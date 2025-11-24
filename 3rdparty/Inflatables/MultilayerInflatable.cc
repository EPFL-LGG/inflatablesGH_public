#include "MultilayerInflatable.hh"
#include "MeshFEM/unused.hh"

#include <MeshFEM/MSHFieldWriter.hh>
#include <MeshFEM/ParallelAssembly.hh>
#include <MeshFEM/filters/remove_dangling_vertices.hh>

// #include <unsupported/Eigen/MPRealSupport>

// Construct from a triangle mesh representing the "top" sheet.
// The "bottom" sheet is an oppositely oriented copy.
MultilayerInflatable::MultilayerInflatable(const std::shared_ptr<Mesh> &inMesh, const size_t num_sheets, const std::vector<Real> &pressures, const Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> &reducedVarIdxForVertexOnSheet)
        : m_sheetMesh(inMesh), m_num_sheets(num_sheets), m_reducedVarIdxForVertexOnSheet(reducedVarIdxForVertexOnSheet), m_pressures(pressures) {
    const auto &m = mesh();
    const size_t nv = m.numVertices(),
                 nt = m.numTris();

    if (num_sheets < 2) throw std::runtime_error("MultilayerInflatable requires at least 2 sheets");
    if (pressures.size() != num_sheets - 1) throw std::runtime_error("Incorrect number of pressures");
    if (reducedVarIdxForVertexOnSheet.rows() != nv) throw std::runtime_error("Incorrect reducedVarIdxForVertexOnSheet size");
    if (reducedVarIdxForVertexOnSheet.cols() != num_sheets) throw std::runtime_error("Incorrect reducedVarIdxForVertexOnSheet size");

    m_useHessianProjectedEnergy = std::vector<bool>(num_sheets * inMesh->numTris(), false);
    // Get the max reduced var index
    m_numReducedVertices = 0;
    std::vector<size_t> vertex_multiplicity(nv, 0);
    for (const auto &v : m.vertices()) {
        for (size_t sheetIdx = 0; sheetIdx < num_sheets; ++sheetIdx) {
            m_numReducedVertices = std::max(m_numReducedVertices, m_reducedVarIdxForVertexOnSheet(v.index(), sheetIdx));
            if (m_reducedVarIdxForVertexOnSheet(v.index(), sheetIdx) < nv)
                vertex_multiplicity[m_reducedVarIdxForVertexOnSheet(v.index(), sheetIdx)] += 1;
        }
    }
    m_numReducedVertices += 1;

    m_updateB();

    setIdentityDeformation();
    m_updateMaterialProperties();
    m_updateGravityCache();

    {
        // Pick "c" and place it at the origin.
        int c_idx;
        // get the center of mass of the mesh
        V3d cm = V3d::Zero();
        for (const auto v : m.vertices()) {
            cm += v.node()->p;
        }
        cm /= m.numVertices();
        // find the index of the vertex closest to the center of mass that's not a fuse vertex.
        Real curr_dist = std::numeric_limits<Real>::max();
        for (const auto v : m.vertices()) {
            if (vertex_multiplicity[v.index()] > 1) continue;
            if ((v.node()->p - cm).norm() < curr_dist) {
                curr_dist = (v.node()->p - cm).norm();
                c_idx = v.index();
            }
        }
        // use varIdx to compute the variable indices corresponding to the vertex c_idx in the top sheet for all x, y, z components.
        m_center_non_fused_vx_idx = varIdx(0, c_idx, 0) / 3;
    }

    m_referenceVolume = std::vector<Real>(num_sheets - 1, 0.0);
    m_clearCache();

}

void MultilayerInflatable::m_updateB() {
    // Generate an orthonormal basis for the tangent plane of each triangle.
    // (All sheets share the same basis)
    const auto &m = mesh();
    const size_t nt = m.numTris();
    m_B.reserve(m_num_sheets * nt);

    // First, check if we actually have a plate in the z = 0 plane; in this
    // case we use the global 2D coordinate system's axis vectors as our
    // orthonormal basis to ease specification of anisotropic materials.
    if (std::abs(m.boundingBox().dimensions()[2]) < 1e-16) {
        M32d globalB(M32d::Identity());
        m_B.assign(nt, globalB);
    }
    else {
        m_B.resize(nt);
        for (auto tri : m.elements()) {
            V3d b0 = (tri.node(1)->p - tri.node(0)->p).cast<Real>().normalized();
            V3d b1 = tri->normal().cast<Real>().cross(b0);
            const size_t ti = tri.index();
            m_B[ti] << b0, b1;
        }
    }
    for (size_t sheetIdx = 1; sheetIdx < m_num_sheets; ++sheetIdx) {
        for (size_t ti = 0; ti < nt; ++ti)
            m_B.push_back(m_B[ti]);
    }
}

void MultilayerInflatable::setVars(Eigen::Ref<const VXd> vars) {
    BENCHMARK_START_TIMER_SECTION("MultilayerInflatable setVars");
    if (size_t(vars.size()) != numVars()) throw std::runtime_error("Invalid variable size");
    m_currVars = vars;

    // Compute Jacobian and per-triangle energy under the new deformation
    const size_t nt = mesh().numTris();
    size_t total_triangles = m_num_sheets * nt;
    m_J .resize(total_triangles);
    m_JB.resize(total_triangles);
    m_triEnergyDensity.resize(total_triangles);
    m_BtGradLambda.resize(total_triangles);
    m_deformed_normals.resize(3, total_triangles);
    m_deformed_normals_scaled_by_areas.resize(3, total_triangles);
    m_deformed_areas  .resize(total_triangles);
    m_projectedTriEnergyDensity.resize(total_triangles);
    auto process_tri = [&](const size_t ti) {
        const auto &tri = mesh().element(ti);
        // const auto &gradLambda = tri->gradBarycentric();
        const auto gradLambda = tri->gradBarycentric().cast<Real>();
        for (size_t sheetIdx = 0; sheetIdx < m_num_sheets; ++sheetIdx) {
            size_t sheet_tri_idx = sheetTriIdx(sheetIdx, tri.index());
            M3d triCornerPos;
            getDeformedTriCornerPositions(ti, sheetIdx, triCornerPos);
            M32d &JB = m_JB[sheet_tri_idx];
            M3d  &J  = m_J [sheet_tri_idx];

            m_BtGradLambda[sheet_tri_idx] = m_B[sheet_tri_idx].transpose() * gradLambda;
            J  = triCornerPos * gradLambda.transpose();
            JB = J * m_B[sheet_tri_idx];
            m_triEnergyDensity[sheet_tri_idx].setDeformationGradient(JB);
            if (m_useHessianProjectedEnergy[sheet_tri_idx])
                m_projectedTriEnergyDensity[sheet_tri_idx].setF(JB);

            const V3d n = (triCornerPos.col(1) - triCornerPos.col(0)).cross(triCornerPos.col(2) - triCornerPos.col(0));
            const Real dblA = n.norm();
            m_deformed_areas[sheet_tri_idx] = 0.5 * dblA;
            m_deformed_normals.col(sheet_tri_idx) = n / dblA;
            m_deformed_normals_scaled_by_areas.col(sheet_tri_idx) = 0.5 * n;
        }
    };

#if MESHFEM_WITH_TBB
    tbb::parallel_for(tbb::blocked_range<size_t>(0, nt), [&](const tbb::blocked_range<size_t> &b) { for (size_t ti = b.begin(); ti < b.end(); ++ti) process_tri(ti); });
#else
    for (size_t ti = 0; ti < nt; ++ti) process_tri(ti);
#endif
    BENCHMARK_STOP_TIMER_SECTION("MultilayerInflatable setVars");
}

// Set the current deformed configuration equal to a rigid transformation
// of the top sheet mesh vertex positions "P". This rigid transformation is
// chosen to place the vertex "c" closest to the center of mass at the origin,
// place the furthest vertex "p" from "c" at (p_x, 0, 0) (defining the x axis)
// and place the furthest vertex "q" from the new x axis at "(q_x, q_y, 0)".
// This allows us to efficiently constrain the sheet's rigid motion with 6
// variable pin constraints (c = 0, p_y = p_z = q_z = 0).
void MultilayerInflatable::setUninflatedDeformation(M3Xd P /* copy modified inside */, bool prepareRigidMotionPinConstraints) {
    VXd vars(numVars());

    if (prepareRigidMotionPinConstraints) {
        // Pick "c" and place it at the origin.
        int c_idx;
        V3d cm = P.rowwise().mean();
        (P.colwise() - cm).colwise().squaredNorm().minCoeff(&c_idx);
        P.colwise() -= P.col(c_idx).eval();

        // Pick "p", defining the unit x axis vector "x_hat"
        int p_idx;
        P.colwise().squaredNorm().maxCoeff(&p_idx);
        V3d x_hat = P.col(p_idx).normalized();

        // Pick "q", defining the unit y axis vector "y_hat"
        int q_idx;
        P.colwise().cross(x_hat).colwise().squaredNorm().maxCoeff(&q_idx);
        // (P - x_hat * (x_hat.transpose() * P)).colwise().squaredNorm().maxCoeff(q_idx);
        V3d y_hat = (P.col(q_idx) - x_hat.dot(P.col(q_idx)) * x_hat).normalized();
        V3d z_hat = x_hat.cross(y_hat);

        M3d R;
        R << x_hat.transpose(),
             y_hat.transpose(),
             z_hat.transpose();

        P = R * P;

        m_rigidMotionPinVars[0] = varIdx(0, c_idx, 0);
        m_rigidMotionPinVars[1] = varIdx(0, c_idx, 1);
        m_rigidMotionPinVars[2] = varIdx(0, c_idx, 2);

        m_rigidMotionPinVars[3] = varIdx(0, p_idx, 1);
        m_rigidMotionPinVars[4] = varIdx(0, p_idx, 2);

        m_rigidMotionPinVars[5] = varIdx(0, q_idx, 2);
    }
    else {
        m_rigidMotionPinVars.fill(0);
    }

    for (const auto v : mesh().vertices()) {
        for (size_t sheetIdx = 0; sheetIdx < m_num_sheets; ++sheetIdx) {
            vars.segment<3>(varIdx(sheetIdx, v.index())) = P.col(v.index()).transpose();
        }
    }

    setVars(vars);
}

// enclosed volume
std::vector<MultilayerInflatable::Real> MultilayerInflatable::volume() const {
    // For each layer, consider the volume between its top sheet and bottom sheet.
    std::vector<Real> result(m_num_sheets - 1, 0);
    for (size_t layer_idx = 0; layer_idx < m_num_sheets - 1; ++layer_idx) {
        Real vol_6 = 0;
        M3d triCornerPos;

        // We expect the reference volume to be close to the current volume and therefore
        // center our volume calculation around it to reduce floating point cancellation.
        vol_6 = -6 * m_referenceVolume[layer_idx];

        MultilayerInflatable::NeumaierSum<Real> sum(vol_6);

        for (const auto tri : mesh().elements()) {

#if 0
            M3d triCornerPosBot;
            // Attempt at a more numerically robust formula--doesn't seem to make much difference,
            // and if anything is slightly less accurate.
            getDeformedTriCornerPositions(tri.index(), layer_idx, triCornerPos);
            getDeformedTriCornerPositions(tri.index(), layer_idx + 1, triCornerPosBot);

            V3d double_nA_top, double_nA_bot;
            double_nA_top = (triCornerPos   .col(1) - triCornerPos   .col(0)).cross(triCornerPos   .col(2) - triCornerPos   .col(0));
            double_nA_bot = (triCornerPosBot.col(2) - triCornerPosBot.col(0)).cross(triCornerPosBot.col(1) - triCornerPosBot.col(0));
            V3d three_c   = (triCornerPos   .rowwise().sum() + triCornerPosBot.rowwise().sum()) / 2;
            Real contrib  = (triCornerPos   .rowwise().sum() - three_c).dot(double_nA_top)
                        + (triCornerPosBot.rowwise().sum() - three_c).dot(double_nA_bot)
                        + three_c.dot(double_nA_bot + double_nA_top);
            vol_6 += contrib / 3;
#else
            getDeformedTriCornerPositions(tri.index(), layer_idx, triCornerPos);
            Real triContrib = triCornerPos.determinant(); // Sum top/bottom sheet contrib first to reduce floating point error
            getDeformedTriCornerPositions(tri.index(), layer_idx + 1, triCornerPos);
            triContrib     -= triCornerPos.determinant();
            sum.accumulate(triContrib);
#endif
        }
        result[layer_idx] = sum.result() / 6.0 + m_referenceVolume[layer_idx];
        // return vol_6 / 6.0 + m_referenceVolume;
    }
    return result;
}

MultilayerInflatable::Real MultilayerInflatable::energyPressurePotential() const {
    Real result = 0;
    std::vector<Real> curr_volume = this->volume();
    for (size_t layer_idx = 0; layer_idx < m_num_sheets - 1; ++layer_idx) {
        result += (curr_volume[layer_idx] - m_referenceVolume[layer_idx]) * m_pressures[layer_idx];
    }
    return -result;
}

MultilayerInflatable::VXd MultilayerInflatable::gradientPressurePotential(bool handleOpenBoundary) const {
    VXd result(VXd::Zero(numVars()));
    for (size_t layer_idx = 0; layer_idx < m_num_sheets - 1; ++layer_idx) {
        Real pressure = m_pressures[layer_idx];
        for (int idx = 0; idx < 2; ++idx) {
            for (const auto tri : mesh().elements()) {
    #if 1    
                const size_t sheet_tri_idx = sheetTriIdx((layer_idx + idx), tri.index());
                Real normalSign = (idx == 0) ? 1.0 : -1.0;
                V3d contrib = (-pressure / 3.0) * normalSign * m_deformed_normals_scaled_by_areas.col(sheet_tri_idx);
                for (const auto v : tri.vertices())
                    result.segment<3>(varIdx((layer_idx + idx), v.index())) += contrib;

                if (handleOpenBoundary) {
                    // Need to handle boundary verties where the one ring is not complete.
                    // These are not computed if the flag is set to false by the caller where the boundary volumes are handled there. 
                    for (size_t i = 0; i < 3; ++i) {
                        const auto edge = tri.halfEdge(i);
                        if (edge.isBoundary()) {
                            const size_t v0_idx = edge.boundaryEdge().node(0).volumeNode().index();
                            const size_t v1_idx = edge.boundaryEdge().node(1).volumeNode().index();
                            V3d bdryContrib = (-pressure) * (1 - (idx) * 2) * getDeformedVtxPosition(v0_idx, (layer_idx + idx)).cross(getDeformedVtxPosition(v1_idx, (layer_idx + idx))) / 6.0; 
                            result.segment<3>(varIdx((layer_idx + idx), v0_idx)) += bdryContrib; 
                            result.segment<3>(varIdx((layer_idx + idx), v1_idx)) += bdryContrib;
                        }
                        
                    }
                }

    #else // equivalent version derived more directly from the signed volume pressure potential
                const double normalSign = (idx == 0) ? 1.0 : -1.0;
                const double signed_pressure_div_6 = normalSign * pressure / 6.0;
                M3d triCornerPos;
                getDeformedTriCornerPositions(tri.index(), layer_idx + idx, triCornerPos);
                for (const auto v : tri.vertices()) {
                    result.segment<3>(varIdx(layer_idx + idx, v.index())) -=
                        signed_pressure_div_6 * triCornerPos.col((v.localIndex() + 1) % 3)
                                        .cross(triCornerPos.col((v.localIndex() + 2) % 3));
                }
    #endif
            }
        }
    }
    return result;

}

MultilayerInflatable::Real MultilayerInflatable::energy(EnergyType etype) const {
    BENCHMARK_START_TIMER_SECTION("MultilayerInflatable energy");
    MultilayerInflatable::NeumaierSum<Real> sum;
    if ((etype == EnergyType::Full) || (etype == EnergyType::Pressure))
        sum.accumulate(energyPressurePotential());

    if ((etype == EnergyType::Full) || (etype == EnergyType::Elastic)) {
        for (const auto tri : mesh().elements()) {
            for (size_t sheet_idx = 0; sheet_idx < m_num_sheets; ++sheet_idx) {
                size_t sheet_tri_idx = sheetTriIdx(sheet_idx, tri.index());
                sum.accumulate(tri->volume() * m_triEnergyDensity[sheet_tri_idx].energy());
            }
        }
    }
    if ((etype == EnergyType::Full) || (etype == EnergyType::Gravity))
        sum.accumulate(m_g_grad.dot(getVars()));

    BENCHMARK_STOP_TIMER_SECTION("MultilayerInflatable energy");

    return sum.result();
}

MultilayerInflatable::VXd MultilayerInflatable::gradient(EnergyType etype, bool handleOpenBoundary) const {
    BENCHMARK_START_TIMER_SECTION("MultilayerInflatable gradient");
    VXd result(VXd::Zero(numVars()));

    if ((etype == EnergyType::Full) || (etype == EnergyType::Pressure))
        result += gradientPressurePotential(handleOpenBoundary);

    if ((etype == EnergyType::Full) || (etype == EnergyType::Elastic)) {
        auto accumulatePerTriContrib = [this](size_t tri_idx, VXd &out) {
            const auto &tri = mesh().tri(tri_idx);
            for (size_t sheetIdx = 0; sheetIdx < m_num_sheets; ++sheetIdx) {
                const size_t sheet_tri_idx = sheetTriIdx(sheetIdx, tri.index());
                const auto &ted            = m_triEnergyDensity[sheet_tri_idx];
                const auto &BtGradLambda   = m_BtGradLambda    [sheet_tri_idx];

                M3d dE_dv = ted.denergy() * BtGradLambda;
                for (const auto v : tri.vertices())
                    out.segment<3>(varIdx(sheetIdx, v.index(), 0)) += tri->volume() * dE_dv.col(v.localIndex());
            }
        };

        // assemble_parallel(accumulatePerTriContrib, result, mesh().numElements());
        // The serial version actually seems to be faster... (not enough work is done for each tri).
        const size_t ntri = mesh().numElements();
        for (size_t i = 0; i < ntri; ++i)
            accumulatePerTriContrib(i, result);
    }

    if ((etype == EnergyType::Full) || (etype == EnergyType::Gravity))
        result += m_g_grad;

    BENCHMARK_STOP_TIMER_SECTION("MultilayerInflatable gradient");

    return result;
}


template <typename MatrixType>
void MultilayerInflatable::hessian(MatrixType &H, EnergyType etype, bool get_sparsity) const {
    BENCHMARK_SCOPED_TIMER_SECTION timer("MultilayerInflatable.hessian");
    auto &varLocks = m_getVarLocks();
    auto assemblePerTriContrib = [&](const size_t ti) {
        const auto &tri = mesh().element(ti);
        for (size_t layerIdx = 0; layerIdx < m_num_sheets - 1; ++layerIdx) {
            Real pressure = m_pressures[layerIdx];
            for (size_t idx = 0; idx < 2; ++idx) {
                size_t sheetIdx = layerIdx + idx;
                const size_t sheet_tri_idx = sheetTriIdx(sheetIdx, tri.index());
                Eigen::Matrix<Real, 9, 9> elemH;
                elemH.setZero();

                if ((etype == EnergyType::Full) || (etype == EnergyType::Pressure)) {
                    M3d triCornerPos;
                    const double normalSign = (idx == 0) ? 1.0 : -1.0;
                    const double signed_pressure_div_6 = normalSign * pressure / 6.0;
                    getDeformedTriCornerPositions(tri.index(), sheetIdx, triCornerPos);
                    for (size_t vlb = 0; vlb < 3; ++vlb) {
                        for (size_t vla = 0; vla < vlb; ++vla) { // strict upper triangle only (no vertex self-interaction)
                            // Gradient wrt v1 of a triangle's signed volume contribution is:
                            //      d vol / d v1 = v_2 x  v_3
                            // so differentiating again with respect to v_2 or v_3
                            // gives a cross product matrix -[v_3]_x or [v_2]_x, respectively.
                            // The sign here is referred to as ordering_sign below.
                            const size_t vlother = 3 - (vla + vlb);
                            const double ordering_sign = (vlb == ((vla + 1) % 3)) ? -1.0 : 1.0;
                            V3d contrib = (-signed_pressure_div_6 * ordering_sign) * triCornerPos.col(vlother);

                            elemH(3 * vla + 1, 3 * vlb + 0) +=  contrib[2];
                            elemH(3 * vla + 2, 3 * vlb + 0) += -contrib[1];
                            elemH(3 * vla + 0, 3 * vlb + 1) += -contrib[2];
                            elemH(3 * vla + 2, 3 * vlb + 1) +=  contrib[0];
                            elemH(3 * vla + 0, 3 * vlb + 2) +=  contrib[1];
                            elemH(3 * vla + 1, 3 * vlb + 2) += -contrib[0];
                        }
                    }
                }

                // For the first layer, consider top and bottom sheet; for the other layers, only count the bottom sheet.
                if (((etype == EnergyType::Full) || (etype == EnergyType::Elastic)) && (layerIdx == 0 || idx == 1)) {
                    // Accumulate contribution from sheet triangle (tri, sheetIdx)
                    // Note: we assume that the variables for the components of a sheet vertex position
                    // are contiguous.
                    const auto &BtGradLambda = m_BtGradLambda[sheet_tri_idx];

                    for (size_t vlb = 0; vlb < 3; ++vlb) {
                        VSFJ vol_dF_b(0, tri->volume() * BtGradLambda.col(vlb));
                        for (size_t comp_b = 0; comp_b < 3; ++comp_b) {
                            vol_dF_b.c = comp_b;
                            const bool useHPE = m_useHessianProjectedEnergy[sheet_tri_idx];
                            M32d delta_de = useHPE ? m_projectedTriEnergyDensity[sheet_tri_idx].delta_denergy(vol_dF_b)
                                                :          m_triEnergyDensity[sheet_tri_idx].delta_denergy(vol_dF_b);

                            Eigen::Map<M3d>(elemH.col(3 * vlb + comp_b).data()) += delta_de * BtGradLambda;
                        }
                    }
                }

                elemH.triangularView<Eigen::StrictlyLower>() = elemH.triangularView<Eigen::StrictlyUpper>().transpose();

                for (const auto v_b : tri.vertices()) {
                    size_t b = varIdx(sheetIdx, v_b.index(), 0);
                    for (size_t comp_b = 0; comp_b < 3; ++comp_b) {
                        size_t col = b + comp_b;
                        while (varLocks[col].exchange(true, std::memory_order_acquire)); // lock column
                        for (const auto v_a : tri.vertices()) {
                            size_t a = varIdx(sheetIdx, v_a.index(), 0);
                            if (a > b) continue; // upper triangle only
                            size_t len = std::min<size_t>(3, col - a + 1);
                            int vlb = v_b.localIndex();
                            int vla = v_a.localIndex();
                            H.addNZStrip(a, col, elemH.col(3 * vlb + comp_b).middleRows(3 * vla, len));
                        }
                        varLocks[col].store(false, std::memory_order_release); // unlock column
                    }
                }
            }
        }
    };

    if (get_sparsity) {
        for (size_t ti = 0; ti < mesh().numElements(); ++ti)
            assemblePerTriContrib(ti);
    } else {
        get_hessian_assembly_arena().execute([&assemblePerTriContrib, this]() {
            parallel_for_range(mesh().numElements(), assemblePerTriContrib);
        });
    }
}

std::shared_ptr<MultilayerInflatable::Mesh> MultilayerInflatable::visualizationMesh(bool duplicateFusedTris) const {
    const size_t nv = mesh().numVertices();
    const size_t nt = mesh().numTris();
    std::vector<V3d> vertices(m_num_sheets * nv);
    std::vector<MeshIO::IOElement> elements;
    elements.reserve(m_num_sheets * nt);

    for (size_t sheetIdx = 0; sheetIdx < m_num_sheets; ++sheetIdx) {
        for (size_t vi = 0; vi < nv; ++vi)
            vertices[vi + sheetIdx * nv] = getDeformedVtxPosition(vi, sheetIdx);
        for (const auto tri : mesh().elements()) {
            elements.emplace_back(tri.vertex(0).index() + sheetIdx * nv,
                                  tri.vertex(1).index() + sheetIdx * nv,
                                  tri.vertex(2).index() + sheetIdx * nv);
            if (sheetIdx == m_num_sheets) std::swap(elements.back()[0], elements.back()[1]); // flip orientation of bottom sheet tris
        }
    }
    if (!duplicateFusedTris)
        remove_dangling_vertices(vertices, elements);

    // Note: deduplicating the fused trianges can result in a non-manifold mesh
    return std::make_shared<Mesh>(elements, vertices, /* suppressNonmanifoldWarning = */ !duplicateFusedTris);
}

Eigen::MatrixXd MultilayerInflatable::visualizationField(Eigen::MatrixXd field, bool duplicateFusedTris) {
    const size_t nt        = mesh().numTris();
    const size_t nv        = mesh().numVertices();
          size_t in_size   = field.rows();
    const size_t field_dim = field.cols();
    // std::cout << "Running visualizationField with" << std::endl;
    // std::cout << "nt        = " << nt        << std::endl;
    // std::cout << "nv        = " << nv        << std::endl;
    // std::cout << "in_size   = " << in_size   << std::endl;
    // std::cout << "field_dim = " << field_dim << std::endl;
    // std::cout << "num reduced vertices = " << m_numReducedVertices << std::endl;

    // Duplicate data defined on just the top sheet to the bottom sheet.
    if ((in_size == nt) || (in_size == nv)) {
        field.conservativeResize(m_num_sheets * in_size, field_dim);
        field.bottomRows(in_size) = field.topRows(in_size);
    }

    // Decode "DoF field" into a per-vertex field on top/bottom sheet
    // (e.g., output of `sheet.gradient()` after reshaping into an `N x 3` matrix)
    if (in_size == m_numReducedVertices) { // "DoF field" 
        Eigen::MatrixXd decodedField(m_num_sheets * nv, field_dim);
        for (size_t sheetIdx = 0; sheetIdx < m_num_sheets; ++sheetIdx) {
            for (size_t vi = 0; vi < nv; ++vi)
                decodedField.row(sheetIdx * nv + vi) = field.row(m_reducedVarIdxForVertexOnSheet(vi, sheetIdx));
        }
        field.swap(decodedField);
    }

    in_size = field.rows(); // possibly updated!

    // Deduplicate fields defined on top/bottom sheet when duplicateFusedTris = false.
    if (in_size == m_num_sheets * nt) {
    }
    else if (in_size == m_num_sheets * nv) {
        if (!duplicateFusedTris) {
            size_t outVtx = nv;
            std::vector<bool> include(nv); // whether to include the bottom sheet copy of a vertex
            for (const auto t: mesh().elements()) {
                include[t.vertex(0).index()] = true;
                include[t.vertex(1).index()] = true;
                include[t.vertex(2).index()] = true;
            }
            for (size_t vi = 0; vi < nv; ++vi) {
                if (!include[vi]) {
                    if (field.row(nv + vi) != field.row(vi))
                        throw std::runtime_error("Inconsistent data on deduplicated top/bottom vertices");
                    continue;
                }
                field.row(outVtx++) = field.row(nv + vi);
            }
            field.conservativeResize(outVtx, field_dim);
        }
    }
    else throw std::runtime_error("Unimplemented/unsupported field type");

    return field;
}

void MultilayerInflatable::writeDebugMesh(const std::string &path) const {
    const size_t nv = mesh().numVertices();
    const size_t nt = mesh().numTris();
    std::vector<MeshIO::IOVertex > vertices(m_num_sheets * nv);
    std::vector<MeshIO::IOElement> elements;
    elements.reserve(m_num_sheets * nt);

    VectorField<double, 3> N(m_num_sheets * nt);
    SymmetricMatrixField<double, 2> strain(m_num_sheets * nt); // *2D* rank-2 tensor field

    for (size_t sheetIdx = 0; sheetIdx < m_num_sheets; ++sheetIdx) {
        for (size_t vi = 0; vi < nv; ++vi)
            vertices[vi + sheetIdx * nv].point = getDeformedVtxPosition(vi, sheetIdx).cast<double>();
        for (const auto tri : mesh().elements()) {
            elements.emplace_back(tri.vertex(0).index() + sheetIdx * nv,
                                  tri.vertex(1).index() + sheetIdx * nv,
                                  tri.vertex(2).index() + sheetIdx * nv);
            if (sheetIdx == 1) std::swap(elements.back()[0], elements.back()[1]); // flip orientation of bottom sheet tris
            size_t sheet_tri_idx = sheetTriIdx(sheetIdx, tri.index());
            N(sheet_tri_idx) = m_deformed_normals.col(sheet_tri_idx).cast<double>();
            strain(sheet_tri_idx) = SymmetricMatrixValue<double, 2>(greenLagrangianStrain(sheetIdx, tri.index()).cast<double>());
        }
    }

    MSHFieldWriter writer(path, vertices, elements);
    writer.addField("Normal",          N, DomainType::PER_ELEMENT);
    writer.addField("G-L Strain", strain, DomainType::PER_ELEMENT);
}

template void MultilayerInflatable::hessian<TripletMatrix<Triplet<double>>>(TripletMatrix<Triplet<double>>& , MultilayerInflatable::EnergyType, bool) const;

template void MultilayerInflatable::hessian<SuiteSparseMatrix>(SuiteSparseMatrix& , MultilayerInflatable::EnergyType, bool) const;
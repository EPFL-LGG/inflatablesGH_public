////////////////////////////////////////////////////////////////////////////////
// MultilayerInflatable.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//  An "inflatable sheet" is a structure formed by fusing together two
//  identical sheets of inextensible material along their boundaries and along
//  some internal curves to form air channels.
//
//  The "top" and "bottom" sheet are two oppositely oriented copies of a single
//  planar triangle mesh. The sheets are fused together by making each copy
//  share variables controlling the "fused vertices." This is done by
//  introducing a "reduced vertex set" whose positions determine the positions
//  of all vertices on the top and bottom sheets.
*/
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Created:  04/05/2019 17:46:33
////////////////////////////////////////////////////////////////////////////////
#ifndef MULTILAYERINFLATABLE_HH
#define MULTILAYERINFLATABLE_HH

#include <MeshFEM/FEMMesh.hh>
#include <MeshFEM/SparseMatrices.hh>
#include <MeshFEM/Utilities/ArrayPadder.hh>
#include <MeshFEM/Utilities/MeshConversion.hh>
#include <memory>
#include <string>
#include <atomic>

#include <MeshFEM/EnergyDensities/NeoHookeanEnergy.hh>
#include <MeshFEM/EnergyDensities/StVenantKirchhoff.hh>
#include <MeshFEM/EnergyDensities/TensionFieldTheory.hh>
#include <MeshFEM/EnergyDensities/TangentElasticityTensor.hh>
#include <MeshFEM/EnergyDensities/IsoCRLETensionFieldMembrane.hh>
#include <MeshFEM/EnergyDensities/TensionFieldNeoHookean.hh>

// #include "TensionFieldEnergy.hh"
#include "IncompressibleBalloonEnergyWithHessProjection.hh"

struct MultilayerInflatable {
#if INFLATABLES_LONG_DOUBLE
    using Real = long double;
#else
    using Real = double;
#endif

#if 0
    using INeo_TFT_CBased = RelaxedEnergyDensity<IncompressibleNeoHookeanEnergyCBased<Real>>;
#else
    using INeo_TFT_CBased = OptionalTensionFieldEnergy<MultilayerInflatable::Real>;
#endif
    using StVk_TFT_CBased = RelaxedEnergyDensity<StVenantKirchhoffEnergyCBased<Real, 2>>;
#if 1
    using EnergyDensityCBased = INeo_TFT_CBased;
#else
    using EnergyDensityCBased = StVk_TFT_CBased;
#endif
#if 1
    using EnergyDensity = EnergyDensityFBasedFromCBased<EnergyDensityCBased, 3>;
#else
    using EnergyDensity = IsoCRLETensionFieldMembrane<Real>;
#endif

    using  V2d = Eigen::Matrix<Real, 2, 1>;
    using  V3d = Eigen::Matrix<Real, 3, 1>;
    using  V4d = Eigen::Matrix<Real, 4, 1>;
    using  VXd = Eigen::Matrix<Real, Eigen::Dynamic, 1>;
    using  M2d = Eigen::Matrix<Real, 2, 2>;
    using  M3d = Eigen::Matrix<Real, 3, 3>;
    using M23d = Eigen::Matrix<Real, 2, 3>;
    using M32d = Eigen::Matrix<Real, 3, 2>;
    using MX2d = Eigen::Matrix<Real, Eigen::Dynamic, 2>;
    using MX3d = Eigen::Matrix<Real, Eigen::Dynamic, 3>;
    using M3Xd = Eigen::Matrix<Real, 3, Eigen::Dynamic>;
    using M34d = Eigen::Matrix<Real, 3, 4>;
    using  M4d = Eigen::Matrix<Real, 4, 4>;
    using VSFJ = VectorizedShapeFunctionJacobian<3, V2d>;
    using ETensor = ElasticityTensor<Real, 2>;

    using TMatrix = TripletMatrix<Triplet<Real>>;


    // Neumaier sum adapted from Wikipedia
    template<typename T>
    struct NeumaierSum {
        NeumaierSum(T val = 0) : sum(val) { }

        void accumulate(T term) {
            T newSum = sum + term;
            if (std::abs(sum) >= std::abs(term))
                c += term  + (sum - newSum); // If sum is bigger, low-order digits of "term" are lost.
            else
                c += sum + (term - newSum);  // Else low-order digits of sum are lost
            sum = newSum;
        }

        T result() {  return sum + c; }

        T c = 0; // roundoff error correction
        T sum = 0;
    };

    using Mesh = FEMMesh<2, 1, V3d>; // Piecewise linear triangle mesh embedded in R^3

    // Be careful when changing this since the periodic classes (InflatablePeriodicUnit and InflateableMidSurfacePeriodicUnit) also uses these.
    enum class EnergyType { Full, Elastic, Pressure, Gravity };

    // Build from a triangle mesh, number of sheets, pressure values for adjacent sheets, and a mapping from (mesh_vx_idx, sheet_idx) to var_idx.
    // The "bottom" sheet is an oppositely oriented copy.
    MultilayerInflatable(const std::shared_ptr<Mesh> &inMesh, const size_t num_sheets, const std::vector<Real> &pressures, const Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> &reducedVarIdxForVertexOnSheet);

    void setMaterial(const EnergyDensity &psi) {
        for (auto &ted : m_triEnergyDensity)
            ted.copyMaterialProperties(psi);
    }

    // Update the sheet design (repositioning rest vertices).
    // Also optionally update the equilibrium variables (i.e. the deformed configuration).
    template<class Derived>
    void setRestVertexPositions(const Eigen::MatrixBase<Derived> &X) {
        mesh().setNodePositions(pad_columns<1>(X));
        m_updateB();
        setVars(getVars()); // Also update the gradient quantities
    }
    template<class Derived>
    void setRestVertexPositions(const Eigen::MatrixBase<Derived> &X, Eigen::Ref<const VXd> vars) {
        mesh().setNodePositions(pad_columns<1>(X));
        m_updateB();
        setVars(vars); // Also update the gradient quantities
    }

    size_t numVars() const { return 3 * m_numReducedVertices; }
    const VXd &getVars() const { return m_currVars; }
    void setVars(Eigen::Ref<const VXd> vars);
    size_t numSheetTris() const { return m_num_sheets * mesh().numTris(); }

    // Use a rigid transformation of the passed (top) sheet vertex positions as
    // an initial deformed configuration for both the top and bottom sheets
    // (resulting in an uninflated structure).
    // This rigid transformation is chosen to enable pinning rigid motion with
    // 6 variable pin constraints (subsequently accessed by rigidMotionPinVars()).
    // If "prepareRigidMotionPinConstraints" is false, then the unmodified P is used,
    // but rigid motion pin constraints are not set up.
    void setUninflatedDeformation(M3Xd P /* copy modified inside */, bool prepareRigidMotionPinConstraints = false);
    void setIdentityDeformation(bool prepareRigidMotionPinConstraints = false) {
        M3Xd P(3, mesh().numVertices());
        for (const auto v : mesh().vertices())
            P.col(v.index()) = v.node()->p.cast<Real>();
        setUninflatedDeformation(P, prepareRigidMotionPinConstraints);
    }

    void setUseTensionFieldEnergy(bool useTFE) {
        for (auto &ted : m_triEnergyDensity)
            ted.setRelaxationEnabled(useTFE);
    }

    // Note: enabling the Hessian projected energy necessarily disables the tension field energy.
    void setUseHessianProjectedEnergy(bool useHPE) {
        for (auto &ted : m_projectedTriEnergyDensity)
            ted.applyHessianProjection = useHPE;
        m_useHessianProjectedEnergy.assign(numSheetTris(), useHPE);
        if (useHPE) setVars(getVars()); // The hessian-projected energy density has not necessarily been updated for the current variables...
    }

    // Note: the behavior here is undefined if the sheets' energy density types are inhomogeneous
    bool usingHessianProjectedEnergy(size_t i) const { return m_useHessianProjectedEnergy.at(i); }
    bool usingTensionFieldEnergy(size_t i)     const { return !m_useHessianProjectedEnergy.at(i)
                                                            && m_triEnergyDensity.at(i).getRelaxationEnabled(); }

    std::vector<bool> usingHessianProjectedEnergy() const { return m_useHessianProjectedEnergy; }

    void setRelaxedStiffnessEpsilon(Real val) {
        for (auto &ted : m_triEnergyDensity)
            ted.setRelaxedStiffnessEpsilon(val);
    }

    std::array<size_t, 3> tensionStateHistogram() const {
        std::array<size_t, 3> counts = {{0, 0, 0}};
        for (const auto &ted : m_triEnergyDensity)
            ++counts[ted.tensionState()];
        return counts;
    }

    const std::array<size_t, 6> &rigidMotionPinVars() const { return m_rigidMotionPinVars; }
    void setRigidMotionPinVars(const std::array<size_t, 6> &pinVars) { m_rigidMotionPinVars = pinVars; }

    void setPressure(std::vector<Real> p) {
        if (p.size() != m_num_sheets - 1)
            throw std::runtime_error("Pressure vector size mismatch"); 
        m_pressures = p; 
    }
    std::vector<Real> getPressure() const { return m_pressures; }

    void setThickness   (Real h) { m_thickness    = h; m_updateMaterialProperties(); }
    void setYoungModulus(Real E) { m_youngModulus = E; m_updateMaterialProperties(); }
    Real getThickness()    const { return m_thickness; }
    Real getYoungModulus() const { return m_youngModulus; }

    void setRho (Real rho) { m_g_rho  = rho; m_updateGravityCache(); }
    void setGravity (const V3d &g) { m_g = g; m_updateGravityCache(); }
    Real getRho () const { return m_g_rho; }
    const V3d &getGravity () const { return m_g; }

    // Volume enclosed by the sheet's tubes
    std::vector<Real> volume() const;
    std::vector<Real> referenceVolume() const { return m_referenceVolume; }
    void setReferenceVolume(std::vector<Real> V0) { m_referenceVolume = V0; }

    Real energy(EnergyType etype = EnergyType::Full) const;
    Real energyPressurePotential() const;

    Real systemEnergy() const { return energy(EnergyType::Elastic); }

    VXd gradientPressurePotential(bool handleOpenBoundary = true) const;
    VXd gradient(EnergyType etype = EnergyType::Full, bool handleOpenBoundary = true) const;

    size_t hessianNNZ() const { return hessianSparsityPattern().nz; } // TODO: predict without constructing
    SuiteSparseMatrix hessianSparsityPattern(Real val = 0.0) const {
        if (m_cachedHessianSparsity.m == 0) {
            TMatrix H(numVars(), numVars());
            H.symmetry_mode = TMatrix::SymmetryMode::UPPER_TRIANGLE;
            H.pruneTol = -1.0;
            hessian<TMatrix>(H, EnergyType::Full, true);
            m_cachedHessianSparsity = SuiteSparseMatrix(H);
        }
        m_cachedHessianSparsity.fill(val);
        return m_cachedHessianSparsity;
    }

    template <typename MatrixType>
    void              hessian(MatrixType &H, EnergyType etype = EnergyType::Full, bool get_sparsity = false) const; // accumulate Hessian to H


    SuiteSparseMatrix hessian(EnergyType etype) const {
        SuiteSparseMatrix H(hessianSparsityPattern(0.0));
        hessian(H, etype);
        return H;
    }


          Mesh &mesh()       { return *m_sheetMesh; }
    const Mesh &mesh() const { return *m_sheetMesh; }

    // Access the mesh shared pointer from this instance
    std::shared_ptr<Mesh>       meshPtr()       { return m_sheetMesh; }
    std::shared_ptr<const Mesh> meshPtr() const { return m_sheetMesh; }

    auto getDeformedVtxPosition(size_t vi, size_t sheetIdx) const {
        return m_currVars.segment<3>(3 * m_reducedVarIdxForVertexOnSheet(vi, sheetIdx));
    }

    void getDeformedTriCornerPositions(size_t ti, size_t sheetIdx, M3d &out) const {
        if (sheetIdx > m_num_sheets - 1) throw std::runtime_error("sheetIdx out of bounds");
        const auto &tri = mesh().element(ti);
        for (const auto v : tri.vertices())
            out.col(v.localIndex()) = getDeformedVtxPosition(v.index(), sheetIdx);
    }

    void getDeformedTriCornerDisplacement(size_t ti, size_t sheetIdx, M3d &out) const {
        getDeformedTriCornerPositions(ti, sheetIdx, out);
        const auto &tri = mesh().element(ti);
        for (const auto v : tri.vertices())
            out.col(v.localIndex()) -= v.node()->p.cast<Real>();
    }

    // System variable corresponding to component "compIdx" of vertex "vtxIdx"
    // on top/bottom sheet "sheetIdx"
    size_t varIdx(size_t sheetIdx, size_t vtxIdx, size_t compIdx = 0) const {
        return 3 * m_reducedVarIdxForVertexOnSheet(vtxIdx, sheetIdx) + compIdx;
    }

    size_t sheetTriIdx(size_t sheetIdx, size_t triIdx) const {
        return mesh().numTris() * sheetIdx + triIdx;
    }

    MX3d restVertexPositions() const {
        const auto &m = mesh();
        MX3d result(m.numVertices(), 3);
        for (const auto v : m.vertices())
            result.row(v.index()) = v.node()->p.transpose().cast<Real>();
        return result;
    }

    M2d greenLagrangianStrain(size_t sheetIdx, size_t triIdx) const {
        const auto &JB = m_JB[triIdx + sheetIdx * mesh().numTris()];
        return 0.5 * (JB.transpose() * JB - M2d::Identity());
    }

#if 0
    std::vector<ETensor> tangentElasticityTensors() const {
        std::vector<ETensor> result;
        for (const auto &ted : m_triEnergyDensity) {
            EnergyDensityCBased psi_C;
            psi_C.copyMaterialProperties(ted);
            const M32d &F = ted.getDeformationGradient();
            result.push_back(tangentElasticityTensor(psi_C, F.transpose() * F));
        }
        return result;
    }
#endif

    const aligned_std_vector<EnergyDensity> &triEnergyDensities() const { return m_triEnergyDensity; }
    const M3d &deformationGradient3D(size_t sheet_tri_idx) const { return m_J.at(sheet_tri_idx); }
    // Gradients of the shape functions expressed in the triangle's 2D tangent plane basis
    // (one gradient per column).
    const std::vector<M23d> &shapeFunctionGradients() const { return m_BtGradLambda; }

    const std::vector<M2d> cauchyGreenDeformationTensors() const {
        std::vector<M2d> result(m_JB.size());
        for (size_t sti = 0; sti < m_JB.size(); ++sti) {
            result[sti] = m_JB[sti].transpose() * m_JB[sti];
        }
        return result;
    }
    const VXd &deformedAreas() const { return m_deformed_areas; }
    VXd      undeformedAreas() const {
        const size_t nt = mesh().numTris();
        VXd result(2 * nt);
        for (const auto tri : mesh().tris())
            result[tri.index() + nt] = result[tri.index()] = tri->volume();
        return result;
    }

    std::shared_ptr<Mesh> visualizationMesh(bool duplicateFusedTris = false) const;

    Eigen::MatrixXd visualizationField(Eigen::MatrixXd field, bool duplicateFusedTris = false);

    void writeDebugMesh(const std::string &path) const;

    // Helper routines for serialization/restore
    using MaterialConfiguration = std::tuple<Real, // m_triEnergyDensity          stiffness
                                             bool, // m_triEnergyDensity          useTensionField
                                             Real, // m_projectedTriEnergyDensity stiffness
                                             bool>;// m_projectedTriEnergyDensity applyHessianProjection
    std::vector<MaterialConfiguration> getMaterialConfiguration() const {
        std::vector<MaterialConfiguration> result;
        const size_t nst = numSheetTris();
        if (m_triEnergyDensity.size() != nst) throw std::runtime_error("Material configuration size mismatch - triEnergyDensity");
        if (m_projectedTriEnergyDensity.size() != nst) throw std::runtime_error("Material configuration size mismatch - projectedTriEnergyDensity");
        for (size_t i = 0; i < nst; ++i) {
            result.emplace_back(m_triEnergyDensity[i].stiffness(), m_triEnergyDensity[i].getRelaxationEnabled(),
                                m_projectedTriEnergyDensity[i].stiffness, m_projectedTriEnergyDensity[i].applyHessianProjection);
        }
        return result;
    }
    void applyMaterialConfiguration(const std::vector<MaterialConfiguration> &c) {
        const size_t nst = numSheetTris();
        if ((m_triEnergyDensity.size() != nst) || (m_projectedTriEnergyDensity.size() != nst))
            throw std::runtime_error("Material configuration size mismatch");
        for (size_t i = 0; i < nst; ++i) {
            m_triEnergyDensity[i].setStiffness(                     std::get<0>(c[i]));
            m_triEnergyDensity[i].setRelaxationEnabled(             std::get<1>(c[i]));
            m_projectedTriEnergyDensity[i].stiffness              = std::get<2>(c[i]);
            m_projectedTriEnergyDensity[i].applyHessianProjection = std::get<3>(c[i]);
        }
    }

    const std::vector<M32d> &getJB() const { return m_JB; }

    size_t center_non_fused_vx_idx() const { return m_center_non_fused_vx_idx; }
    
    void m_clearCache() { m_cachedHessianSparsity.clear(); }

private:
    std::shared_ptr<Mesh> m_sheetMesh;
    size_t m_num_sheets;

    int m_numReducedVertices;
    Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> m_reducedVarIdxForVertexOnSheet;


    VXd m_currVars;
    std::vector<Real> m_pressures = std::vector<Real>(); // current inflation pressure
    std::vector<Real> m_referenceVolume = std::vector<Real>(); // V_0 used for defining the pressure potential.
    Real m_thickness = 0.075, m_youngModulus = 300; // Material properties used to set the energy density stiffness.

    // Gravity
    Real m_g_rho = 0.0; // mass density of the sheet
    V3d m_g = V3d(0.0, 0.0, -9.80635);
    VXd m_g_grad;

    void m_updateGravityCache() {
        VXd result;
        result.setZero(numVars());
        const auto &m = mesh();
        auto integratedPhis = integratedShapeFunctions<1, 2>();
        for (const auto e : m.elements()) {
            for (const auto n : e.nodes()) {
                auto vi_top = varIdx(0, n.index(), 0);
                auto vi_bot = varIdx(1, n.index(), 0);
                auto contrib = m_g * (integratedPhis[n.localIndex()] * e->volume());
                result.template segment<3>(vi_top) += contrib;
                result.template segment<3>(vi_bot) += contrib;
            }
        }
        result *= -m_g_rho;
        m_g_grad = result;
    }

    // Set all energy densities' stiffness parameters based on the thickness,
    // Young's modulus parameters configured for this sheet.
    void m_updateMaterialProperties() {
        const size_t nst = numSheetTris();
        Real stiffness = m_youngModulus * m_thickness / 6.0;
        for (size_t i = 0; i < nst; ++i) {
            m_triEnergyDensity[i].setStiffness(stiffness);
            m_projectedTriEnergyDensity[i].stiffness = stiffness;
        }
    }

    // Orthonormal basis for each triangle's tangent space (both top and bottom sheet)
    std::vector<M32d> m_B;

    // Method to update the tangent space basis for each triangle (call when rest positions change)
    void m_updateB();
    
    mutable SuiteSparseMatrix m_cachedHessianSparsity;

    ////////////////////////////////////////////////////////////////////////////
    // Quantities computed from the current deformation
    ////////////////////////////////////////////////////////////////////////////
    // Jacobian for each triangle (mapping from the triangle's 2D tangent space to 3D)
    // in the top sheet version (first) and bottom sheet version (after) for each sheet.
    std::vector<M32d> m_JB;
    std::vector<M3d > m_J;  // mapping from 3D to 3D (with J n = 0)
    std::vector<M23d> m_BtGradLambda;
    M3Xd m_deformed_normals_scaled_by_areas;
    M3Xd m_deformed_normals;
    VXd m_deformed_areas;

    aligned_std_vector<EnergyDensity> m_triEnergyDensity;
    aligned_std_vector<IncompressibleBalloonEnergyWithHessProjection<Real>> m_projectedTriEnergyDensity;
    std::vector<bool> m_useHessianProjectedEnergy;

    std::array<size_t, 6> m_rigidMotionPinVars;

    size_t m_center_non_fused_vx_idx; 
    

    // Spin locks used for parallel Hessian assembly.
    mutable std::unique_ptr<std::vector<std::atomic<bool>>> m_varLocks;
    auto &m_getVarLocks() const {
        if (!m_varLocks) {
            const size_t nv = numVars();
            m_varLocks = std::make_unique<std::vector<std::atomic<bool>>>(nv);
            for (size_t i = 0; i < nv; ++i)
                atomic_init(&(*m_varLocks)[i], false);
        }
        return *m_varLocks;
    }
};

#endif /* end of include guard: MULTILAYERINFLATABLE_HH */

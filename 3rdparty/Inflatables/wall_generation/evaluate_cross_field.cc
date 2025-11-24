#include "evaluate_stripe_field.hh"
#include "../subdivide_triangle.hh"

// DDG includes
#include <Mesh.h>
#include <DenseMatrix.h>

double get_distance_to_line_segment(const Eigen::Vector2d &p, const Eigen::Vector2d &a, const Eigen::Vector2d &b) {
    Eigen::Vector2d v1 = p - a;
    Eigen::Vector2d v2 = b - a;
    double dot = v1.dot(v2);
    double len_sq = v2.dot(v2);
    double param = -1.;
    if (len_sq != 0) param = dot / len_sq;
    Eigen::Vector2d closest_pt;
    if (param < 0) {
        closest_pt = a;
    } else if (param > 1) {
        closest_pt = b;
    } else {
        closest_pt = a + param * v2;
    }
    Eigen::Vector2d d = p - closest_pt;
    return d.norm();
}

// The pattern needs to have reflectional symmetry. This function returns the gamma and theta in the base lower left corner square.
std::array<double, 2> get_theta_gamma(const std::array<DDG::Complex, 3> &g, double nu, double nv, const Eigen::Vector3d &b /* barycentric coordinates */) {
    double theta = g[0].re * b[0] +
                   g[1].re * b[1] +
                   g[2].re * b[2];

    double gamma = g[0].im * b[0] +
                   g[1].im * b[1] +
                   g[2].im * b[2];
    // compute lArg_n
    double lArg = 0.;
    if      (b[2] <= b[0] && b[2] <= b[1]) lArg = (M_PI / 3.0) * (1.0 + (b[1] - b[0]) / (1.0 - 3.0 * b[2]));
    else if (b[0] <= b[1] && b[0] <= b[2]) lArg = (M_PI / 3.0) * (3.0 + (b[2] - b[1]) / (1.0 - 3.0 * b[0]));
    else                                   lArg = (M_PI / 3.0) * (5.0 + (b[0] - b[2]) / (1.0 - 3.0 * b[1]));

    gamma += lArg * nv;
    // adjust texture coordinates
    theta += lArg * nu;
    
    theta = std::fmod(theta, 2 * M_PI);
    if (theta < 0) theta += 2 * M_PI;
    // theta now in [0, 2 pi].
    // Over [0, pi], distance should interpolate down from pi to 0,
    // then back up from 0 to pi over the range [pi, 2 * pi].

    if (theta < M_PI) theta = M_PI - theta;
    else theta = theta - M_PI;
    
    // return  theta + (M_PI - theta) * (M_PI - theta) * 0.1 * sin(beta / 2);
    // if (theta > M_PI / 2) return theta;

    gamma = std::fmod(gamma, 2 * M_PI);
    if (gamma < 0) gamma += 2 * M_PI;
    if (gamma < M_PI) gamma = M_PI - gamma;
    else gamma = gamma - M_PI;
    return std::array<double, 2>({theta, gamma});
}

using PatternInfoFunction = std::function<double(const double theta, const double gamma, std::vector<double> patternParams, double margin)>;
using TriangleTextureCoords = std::array<std::array<double, 2>, 3>;
using PolylineTextureBarycentricCoords = std::vector<std::array<Eigen::Vector3d, 2>>;
using PatternPolylineIntersectionFunction = std::function<PolylineTextureBarycentricCoords(const TriangleTextureCoords ttc, std::vector<double> patternParams)>;

// Rewrite the evaluate_cross_field function to use a new argument that specifies the pattern function
// the argument should be of type PatternInfoFunction
// also replace p1 p2 p3 with a single argument that is a vector of pattern parameters of the type const std::vector<std::vector<double>> &patternParams
void evaluate_cross_field_custom_pattern(const Eigen::MatrixX3d &vertices,
                                         const Eigen::MatrixX3i &elements,
                                         const std::vector<double> &stretchAngles,
                                         const std::vector<std::vector<double>> &patternParams,
                                         PatternInfoFunction patternFunc,
                                         PatternPolylineIntersectionFunction patternPolylineFunc,
                                         const double frequency,
                                         const double margin,
                                         Eigen::MatrixX3d &outVerticesEigen,
                                         Eigen::MatrixX3i &outTrianglesEigen,
                                         std::vector<double> &stripeField,
                                         PolylineTextureBarycentricCoords &patternPolylineSoup,
                                         const size_t nsubdiv, 
                                         const bool glue) {
    DDG::Mesh m;
    m.import(vertices, elements);
    std::vector<MeshIO::IOVertex > outVertices;
    std::vector<MeshIO::IOElement> outTriangles;
    m.nCoordinateFunctions = 2;
    if (stretchAngles.size() == 0)
        m.computeCurvatureAlignedSection();
    else {
        const size_t nv = m.vertices.size();
        if (stretchAngles.size() != nv) throw std::runtime_error("Field size mismatch");
        // Compute angle in the tangent plane between a vertex's reference halfedge and the input field vector.
        for (size_t i = 0; i < nv; ++i) {
            auto &vtx = m.vertices[i];
            const DDG::Vector dir(std::cos(stretchAngles[i]), std::sin(stretchAngles[i]), 0.0);
            auto n = vtx.normal();
            auto dir_tangent    = (dir - n * dot(n, dir)).unit();
            // std::cout << dir << "\t" << dir_tangent << "\t" << n << std::endl;
            auto refdir_tangent = (vtx.he->vector() - n * dot(n, vtx.he->vector())).unit();
            DDG::Complex f(dot(dir_tangent, refdir_tangent), dot(n, cross(refdir_tangent, dir_tangent)));
            vtx.directionField = f * f; // Convention is to work with twice the angle
        }
        m.lambda = frequency;
    }

    m.parameterize();

    PointGluingMap indexForPoint;

    stripeField.clear();
    
    std::vector<std::array<double, 2>> theta_gamma_field;
    std::vector<std::vector<double>> patternParamsField;
    patternParamsField.clear();
    theta_gamma_field.clear();

    size_t num_pattern_params = patternParams.size();
    // Replicate the interpolation performed by the stripe shader,
    // evaluating on a refined mesh...
    for (const auto &f : m.faces) {
        if (f.isBoundary()) continue;

        double k  = f.fieldIndex(2.0);
        // k = 0; // FOR DEBUGGING

        if (k == 0) {
            double nu = f.paramIndex[0];
            double nv = f.paramIndex[1];
            // nu = 0; // FOR DEBUGGING
            int i = 0;

            auto he = f.he;
            std::array<DDG::Complex,    3> g;
            std::array<Eigen::Vector3d, 3> vx_pos;
            std::vector<std::array<double, 3>> curr_pattern_vals;
            curr_pattern_vals.resize(num_pattern_params);

            do {
                g[i] = he->texcoord;
                const auto &pos = he->vertex->position;
                vx_pos[i] << pos[0], pos[1], pos[2];
                for (size_t j = 0; j < num_pattern_params; ++j) {
                    curr_pattern_vals[j][i] = patternParams[j].at(he->vertex->index);
                }

                i++;
                he = he->next;
            } while (he != f.he);

            if (!glue) indexForPoint.clear();

            auto newPt = [&](const Eigen::Vector3d &p_sub, double b0, double b1, double b2) {
                outVertices.emplace_back(p_sub);
                std::vector<double> curr_pattern_vals_interp;
                for (size_t j = 0; j < num_pattern_params; ++j) {
                    curr_pattern_vals_interp.push_back(b0 * curr_pattern_vals[j][0] + b1 * curr_pattern_vals[j][1] + b2 * curr_pattern_vals[j][2]);
                }
                patternParamsField.push_back(curr_pattern_vals_interp);
                auto [theta, gamma] = get_theta_gamma(g, nu, nv, Eigen::Vector3d(b0, b1, b2));
                theta_gamma_field.push_back(std::array<double, 2>({theta, gamma}));
                stripeField.push_back(patternFunc(theta, gamma, curr_pattern_vals_interp, margin)); //linearly interpolated pattern parameters
                return outVertices.size() - 1;
            };
            auto newTri = [&](size_t i0, size_t i1, size_t i2) { 
                // Average vertex pattern parameters to get triangle pattern parameters
                // The three indices are the index of the vertices of the triangle in the outVertices vector. 
                // Right now we have different pattern parameters at each vertices, if we switch to using a single pattern parameter at a triangle, how do we stitch the triangles together 
                std::vector<double> curr_pattern_vals_interp;
                for (size_t j = 0; j < num_pattern_params; ++j) {
                    curr_pattern_vals_interp.push_back((patternParamsField[i0][j] + patternParamsField[i1][j] + patternParamsField[i2][j]) / 3.);
                }
                TriangleTextureCoords theta_gamma;
                theta_gamma[0] = theta_gamma_field[i0];
                theta_gamma[1] = theta_gamma_field[i1];
                theta_gamma[2] = theta_gamma_field[i2];
                PolylineTextureBarycentricCoords barycentricSoup = patternPolylineFunc(theta_gamma, curr_pattern_vals_interp);

                // Convert polyline soup in barycentric coords to absolute coords. The first and second points are the barycentric coordinates of two end points of the polyline. The barycentric coordinates are defined for the triangle whose end vertices have index i0, i1, i2.
                Eigen::Vector3d tri_p0 = outVertices[i0].point;
                Eigen::Vector3d tri_p1 = outVertices[i1].point;
                Eigen::Vector3d tri_p2 = outVertices[i2].point;
                for (const auto &polyline : barycentricSoup) {
                    Eigen::Vector3d p0 = tri_p0 * polyline[0][0] + tri_p1 * polyline[0][1] + tri_p2 * polyline[0][2];
                    Eigen::Vector3d p1 = tri_p0 * polyline[1][0] + tri_p1 * polyline[1][1] + tri_p2 * polyline[1][2];
                    patternPolylineSoup.emplace_back(std::array<Eigen::Vector3d, 2>({p0, p1}));
                }

                outTriangles.emplace_back(i0, i1, i2); 
            };

            subdivide_triangle(nsubdiv, vx_pos[0], vx_pos[1], vx_pos[2], indexForPoint, newPt, newTri);

         }
         else // singular triangle
         {
            // Get the three half edges.
            auto hij = f.he;
            auto hjk = hij->next;
            auto hkl = hjk->next;

            // Get the three vertices.
            auto vi = hij->vertex;
            auto vj = hjk->vertex;
            auto vk = hkl->vertex;

            // Get the three parameter values---for clarity, let "l"
            // denote the other point in the same fiber as "i".  Purely
            // for clarity, we will explicitly define the value of psi
            // at l, which of course is always just the conjugate of the
            // value of psi at i.
            DDG::Complex psiI = vi->parameterization;
            DDG::Complex psiJ = vj->parameterization;
            DDG::Complex psiK = vk->parameterization;
            DDG::Complex psiL = psiI.bar();

            double cIJ = ( hij->edge->he != hij ? -1. : 1. );
            double cJK = ( hjk->edge->he != hjk ? -1. : 1. );
            double cKL = ( hkl->edge->he != hkl ? -1. : 1. );

            // Get the three omegas, which were used to define our energy.
            double omegaIJ = hij->omega();
            double omegaJK = hjk->omega();
            double omegaKL = hkl->omega();

            // Here's the trickiest part.  If the canonical orientation of
            // this last edge is from l to k (rather than from k to l)...
            omegaKL *= cKL;
            // SIMPLER // if( cKL == -1. )
            // SIMPLER // {
            // SIMPLER //    // ...then the value of omega needs to be negated, since the
            // SIMPLER //    // original value we computed represents transport away from
            // SIMPLER //    // vertex i rather than the corresponding vertex l
            // SIMPLER //    omegaKL = -omegaKL;
            // SIMPLER // }
            // Otherwise we're ok, because the original value was computed
            // starting at k, which is exactly where we want to start anyway.

            // Now we just get consecutive values along the curve from i to j to k to l.
            // (The following logic was already described in our routine for finding
            // zeros of the parameterization.)
            if( hij->crossesSheets() )
            {
                psiJ = psiJ.bar();
                omegaIJ =  cIJ * omegaIJ;
                omegaJK = -cJK * omegaJK;
            }

            // Note that the flag hkl->crossesSheets() is the opposite of what we want here:
            // based on the way it was originally computed, it flags whether the vectors at
            // Xk and Xi have a negative dot product.  But here, we instead want to know if
            // the vectors at Xk and Xl have a negative dot product.  (And since Xi=-Xl, this
            // flag will be reversed.)
            if( !hkl->crossesSheets() )
            {
                psiK = psiK.bar();
                omegaKL = -cKL * omegaKL;
                omegaJK =  cJK * omegaJK;
            }

            // From here, everthing gets computed as usual.
            DDG::Complex rij( cos(omegaIJ), sin(omegaIJ) );
            DDG::Complex rjk( cos(omegaJK), sin(omegaJK) );
            DDG::Complex rkl( cos(omegaKL), sin(omegaKL) );

            double sigmaIJ = omegaIJ - ((rij*psiI)/psiJ).arg();
            double sigmaJK = omegaJK - ((rjk*psiJ)/psiK).arg();
            double sigmaKL = omegaKL - ((rkl*psiK)/psiL).arg();
            //double xi = sigmaIJ + sigmaJK + sigmaKL;

            double betaI = psiI.arg();
            double betaJ = betaI + sigmaIJ;
            double betaK = betaJ + sigmaJK;
            double betaL = betaK + sigmaKL;
            double betaM = betaI + (betaL-betaI)/2.;

            std::array<Eigen::Vector3d, 3> p_ijk;
            p_ijk[0] << vi->position[0], vi->position[1], vi->position[2];
            p_ijk[1] << vj->position[0], vj->position[1], vj->position[2];
            p_ijk[2] << vk->position[0], vk->position[1], vk->position[2];
            Eigen::Vector3d pm = (p_ijk[0] + p_ijk[1] + p_ijk[2]) / 3.;

            std::vector<std::array<double, 3>> patternParams_ijk;
            std::vector<double> patternParams_m;
            for (size_t j = 0; j < num_pattern_params; ++j) {
                std::array<double, 3> curr_pattern_vals;
                curr_pattern_vals[0] = patternParams[j].at(vi->index);
                curr_pattern_vals[1] = patternParams[j].at(vj->index);
                curr_pattern_vals[2] = patternParams[j].at(vk->index);
                patternParams_ijk.push_back(curr_pattern_vals);
                patternParams_m.push_back((curr_pattern_vals[0] + curr_pattern_vals[1] + curr_pattern_vals[2]) / 3.);
            }

            const double nu = 0, nv = 0;

            std::array<DDG::Complex, 3> g;
            std::array<DDG::Complex, 4> beta_ijkl = {{ DDG::Complex(betaI), DDG::Complex(betaJ), DDG::Complex(betaK), DDG::Complex(betaL) }};

            for (int offset = 0; offset < 3; ++offset) {
                int next = (offset + 1) % 3;
                auto newPt = [&](const Eigen::Vector3d &p, double b0, double b1, double b2) {
                    outVertices.emplace_back(p);
                    std::vector<double> curr_pattern_vals_interp;
                    for (size_t j = 0; j < num_pattern_params; ++j) {
                        curr_pattern_vals_interp.push_back(b0 * patternParams_ijk[j][offset] + b1 * patternParams_ijk[j][next] + b2 * patternParams_m[j]);
                    }
                    patternParamsField.push_back(curr_pattern_vals_interp);
                    auto [theta, gamma] = get_theta_gamma(g, nu, nv, Eigen::Vector3d(b0, b1, b2));
                    theta_gamma_field.push_back(std::array<double, 2>({theta, gamma}));
                    stripeField.push_back(patternFunc(theta, gamma, curr_pattern_vals_interp, margin)); //linearly interpolated pattern parameters
                    return outVertices.size() - 1;
                };

                // Repeat of the function in the branch above, but uses different g, nu, nv.
                auto newTri = [&](size_t i0, size_t i1, size_t i2) { 
                    // Average vertex pattern parameters to get triangle pattern parameters
                    // The three indices are the index of the vertices of the triangle in the outVertices vector. 
                    // Right now we have different pattern parameters at each vertices, if we switch to using a single pattern parameter at a triangle, how do we stitch the triangles together 
                    std::vector<double> curr_pattern_vals_interp;
                    for (size_t j = 0; j < num_pattern_params; ++j) {
                        curr_pattern_vals_interp.push_back((patternParamsField[i0][j] + patternParamsField[i1][j] + patternParamsField[i2][j]) / 3.);
                    }
                    TriangleTextureCoords theta_gamma;
                    theta_gamma[0] = theta_gamma_field[i0];
                    theta_gamma[1] = theta_gamma_field[i1];
                    theta_gamma[2] = theta_gamma_field[i2];
                    PolylineTextureBarycentricCoords barycentricSoup = patternPolylineFunc(theta_gamma, curr_pattern_vals_interp);

                    // Convert polyline soup in barycentric coords to absolute coords. The first and second points are the barycentric coordinates of two end points of the polyline. The barycentric coordinates are defined for the triangle whose end vertices have index i0, i1, i2.
                    Eigen::Vector3d tri_p0 = outVertices[i0].point;
                    Eigen::Vector3d tri_p1 = outVertices[i1].point;
                    Eigen::Vector3d tri_p2 = outVertices[i2].point;
                    for (const auto &polyline : barycentricSoup) {
                        Eigen::Vector3d p0 = tri_p0 * polyline[0][0] + tri_p1 * polyline[0][1] + tri_p2 * polyline[0][2];
                        Eigen::Vector3d p1 = tri_p0 * polyline[1][0] + tri_p1 * polyline[1][1] + tri_p2 * polyline[1][2];
                        patternPolylineSoup.emplace_back(std::array<Eigen::Vector3d, 2>({p0, p1}));
                    }

                    outTriangles.emplace_back(i0, i1, i2); 
                };
                

                g[0] = beta_ijkl[offset];
                g[1] = beta_ijkl[offset + 1];
                g[2] = DDG::Complex(betaM);
                if (!glue) indexForPoint.clear();
                subdivide_triangle(nsubdiv, p_ijk[offset], p_ijk[next], pm, indexForPoint, newPt, newTri);
            }
        }
    }

    outVerticesEigen .resize(outVertices .size(), 3);
    outTrianglesEigen.resize(outTriangles.size(), 3);
    for (size_t i = 0; i < outVertices.size(); ++i)
        outVerticesEigen.row(i) = outVertices[i].point;
    for (size_t i = 0; i < outTriangles.size(); ++i)
        outTrianglesEigen.row(i) = Eigen::Vector3i(outTriangles[i][0], outTriangles[i][1], outTriangles[i][2]);
}
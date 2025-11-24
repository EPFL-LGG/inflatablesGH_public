////////////////////////////////////////////////////////////////////////////////
// evaluate_stripe_field.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//  Sampling of the stripe pattern on a subdivided mesh.
*/
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Created:  04/27/2019 15:36:32
////////////////////////////////////////////////////////////////////////////////
#ifndef EVALUATE_STRIPE_FIELD_HH
#define EVALUATE_STRIPE_FIELD_HH

#include <MeshFEM/MeshIO.hh>
#include <Eigen/Dense>

void evaluate_stripe_field(const Eigen::MatrixX3d &vertices,
                           const Eigen::MatrixX3i &elements,
                           const std::vector<double> &stretchAngles,
                           const std::vector<double> &wallWidths,
                           const double frequency,
                           Eigen::MatrixX3d &outVerticesEigen,
                           Eigen::MatrixX3i &outTrianglesEigen,
                           std::vector<double> &stripeField,
                           const size_t nsubdiv = 3,
                           const bool glue = true);

using PatternInfoFunction = std::function<double(const double theta, const double gamma, std::vector<double> patternParams, double margin)>;
using TriangleTextureCoords = std::array<std::array<double, 2>, 3>;
// List of pairs of texture coordinates in barycentric coordinates of the triangle.
using PolylineTextureBarycentricCoords = std::vector<std::array<Eigen::Vector3d, 2>>;
// Take texture coordinates of a triangle and the pattern parameters and return the pattern polylines contained within the triangle in barycentric coordinates of the triangle.
using PatternPolylineIntersectionFunction = std::function<PolylineTextureBarycentricCoords(const TriangleTextureCoords ttc, std::vector<double> patternParams)>;



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
                                         const bool glue);

#endif /* end of include guard: EVALUATE_STRIPE_FIELD_HH */

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include "custom_Stream_lines_2.h"
#include <CGAL/Runge_kutta_integrator_2.h>
#include <CGAL/Triangular_field_2.h>

#include <MeshFEM/MeshIO.hh>

#include <iostream>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel         K;
typedef K::Point_2                                                  Point;
typedef K::Vector_2                                                 Vector;
typedef CGAL::Triangular_field_2<K>                                 Field;
typedef CGAL::Runge_kutta_integrator_2<Field>                       Runge_kutta_integrator;
typedef CGAL::Custom_Stream_lines_2<Field, Runge_kutta_integrator>  Strl;

int main(int argc, const char *argv[]) {
    if (argc != 5) {
        std::cout << "usage: streamlines points vectors separation_distance saturation_ratio" << std::endl;
        std::cout << "example separation_distance: 60.0, saturation_ratio: 1.6" << std::endl;
        exit(-1);
    }

    Runge_kutta_integrator runge_kutta_integrator(1);
    /*datap.tri.cin and datav.tri.cin are ascii files where are stored the vector values*/
    std::ifstream inp(argv[1]);
    std::ifstream inv(argv[2]);
    std::istream_iterator<Point> beginp(inp);
    std::istream_iterator<Vector> beginv(inv);
    std::istream_iterator<Point> endp;
    Field triangular_field(beginp, endp, beginv);
    /* the placement of streamlines */
    std::cout << "processing...\n";
    double dSep = std::stod(argv[3]);
    double dRat = std::stod(argv[4]);
    Strl streamlines(triangular_field, runge_kutta_integrator, dSep, dRat);
    std::cout << "placement generated\n";

    std::vector<MeshIO::IOVertex > vertices;
    std::vector<MeshIO::IOElement> elements;
    for (const auto &polyline : streamlines.stl_container) {
        size_t offset = vertices.size();
        size_t numPts = polyline.size();
        for (const auto &pt : polyline)
            vertices.emplace_back(Point3D(pt.x(), pt.y(), 0.0));
        for (size_t i = 0; i < numPts - 1; ++i)
            elements.emplace_back(offset + i, offset + i + 1);
    }
    MeshIO::save("streamlines.msh", vertices, elements);
}

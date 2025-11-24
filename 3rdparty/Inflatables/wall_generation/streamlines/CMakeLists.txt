find_package(CGAL QUIET)
if (NOT TARGET CGAL::CGAL)
    message("CGAL not found; disabling streamlines binary")
endif()


if (TARGET CGAL::CGAL)
    add_executable(streamlines streamlines.cc)
    target_link_libraries(streamlines CGAL::CGAL MeshFEM)
endif()

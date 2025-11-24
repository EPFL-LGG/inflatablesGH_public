#include <igl/opengl/glfw/Viewer.h>
#include <GLFW/glfw3.h>

#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <igl/file_dialog_open.h>
#include <igl/jet.h>
#include <igl/hsv_to_rgb.h>

#include <igl/barycenter.h>
#include <igl/avg_edge_length.h>
#include <imgui/imgui.h>

#include "examples/paraboloid.hh"

using Viewer     = igl::opengl::glfw::Viewer;
using ViewerCore = igl::opengl::ViewerCore;
using ViewerData = igl::opengl::ViewerData;
using RGBColor = Eigen::Vector3d;

#include "parametrization.hh"
#include "parametrization.hh"
#include "local_global_parametrization.hh"
#include "regularized_parametrization.hh"
#include <MeshFEM/MeshIO.hh>
#include <MeshFEM/MSHFieldWriter.hh>
#include <MeshFEM/GlobalBenchmark.hh>
#include <thread>

using LGParam = parametrization::LocalGlobalParametrizer;

void meshio_to_igl(const std::vector<MeshIO::IOVertex > &vertices, const std::vector<MeshIO::IOElement> &elements,
                   Eigen::MatrixXd &V, Eigen::MatrixXi &F) {
    const size_t nv = vertices.size();
    V.resize(nv, 3);
    for (size_t i = 0; i < nv; ++i)
        V.row(i) = vertices[i].point;

    size_t ne = elements.size();
    if (ne == 0) return;
    size_t es = elements[0].size();

    if (es == 4) {
        ne *= 2; // triangulate quads
    }

    F.resize(ne, 3);

    for (size_t i = 0; i < elements.size(); ++i) {
        const auto &e = elements[i];
        if (es != e.size()) throw std::runtime_error("Mixed element types not supported");
        if (es == 4) {
            F.row(2 * i + 0) << e[0], e[1], e[2];
            F.row(2 * i + 1) << e[0], e[2], e[3];
        }
        else if (es == 3) {
            F.row(i) << e[0], e[1], e[2];
        }
        else throw std::runtime_error("Unsupported mesh type");
    }
}

void plot_vectorfield(ViewerData &data, const Eigen::MatrixXd &points, const Eigen::MatrixXd &vectors, const RGBColor &color, bool append, Real scale = 1.0) {
    int n = vectors.rows();
    if (points.rows() != n) throw std::runtime_error("Tangent vector base point size mismatch");
    if (append) data.lines.conservativeResize(n + data.lines.rows(), 9);
    else        data.lines.resize(n, 9);

    // (Every row contains 9 doubles in the following format S_x, S_y, S_z, T_x, T_y, T_z, C_r, C_g, C_b),
    auto dst = data.lines.bottomRows(n);
    dst.leftCols<6>() << points, (points + scale * vectors);
    dst.rightCols<3>().rowwise() = color.transpose();
    data.dirty |= igl::opengl::MeshGL::DIRTY_OVERLAY_LINES;
}

enum class ViewingField : int { NONE=0, ALPHA=1, PHI=2 };

int main(int argc, char * argv[])
{
    if (argc != 2) {
        std::cerr << "Usage: viewer [surface.obj]" << std::endl;
        exit(-1);
    }

    const size_t width = 1280;
    const size_t height = 800;

    // Create viewer with two empty mesh slots.
    Viewer viewer;
    viewer.append_mesh();
    const int uv_mesh_id = 0, surface_mesh_id = 1;

    // State of the mesh/parametrization.
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    Eigen::MatrixXd barycenters3D, barycenters2D;
    parametrization::UVMap uv;
    std::unique_ptr<LGParam> param;

    double gui_alpha_min = 1.0;
    double gui_alpha_max = 1.0;
    double global_scale = 1.0;

    bool data_updated = false;
    bool show_frame = true;
    ViewingField view_field = ViewingField::NONE;

    auto update_data = [&]() {
        V.leftCols(2) = uv;
        V.col(2).setZero();
        // Scale flattened mesh into [-1, 1]
        V.rowwise() -= V.colwise().mean().eval();
        V *= 1.0 / V.cwiseAbs().maxCoeff();

        viewer.data(uv_mesh_id).set_mesh(V, F);

        igl::barycenter(V, F, barycenters2D);

        auto stretchFrame2D = param->scaledPrincipalDirections(parametrization::Domain::UV);
        auto stretchFrame3D = param->scaledPrincipalDirections(parametrization::Domain::XYZ);

        if (show_frame) {
            plot_vectorfield(viewer.data(uv_mesh_id), barycenters2D, stretchFrame2D.first,  RGBColor(1.0, 0.0, 0.0), false, global_scale);
            plot_vectorfield(viewer.data(uv_mesh_id), barycenters2D, stretchFrame2D.second, RGBColor(0.0, 0.0, 1.0),  true, global_scale);

            plot_vectorfield(viewer.data(surface_mesh_id), barycenters3D, stretchFrame3D.first,  RGBColor(1.0, 0.0, 0.0), false, global_scale);
            plot_vectorfield(viewer.data(surface_mesh_id), barycenters3D, stretchFrame3D.second, RGBColor(0.0, 0.0, 1.0),  true, global_scale);
        }
        else {
            viewer.data(uv_mesh_id).lines.resize(0, 9);
            viewer.data(surface_mesh_id).lines.resize(0, 9);
        }

        if (view_field == ViewingField::ALPHA) {
            const auto &alphas = param->getAlphas();

            Eigen::MatrixXd C;
            // map [1, pi / 2] to [0, 1]
            igl::jet(((Eigen::Map<const Eigen::VectorXd>(alphas.data(), alphas.size()) - Eigen::VectorXd::Ones(alphas.size())) /  (M_PI / 2.0 - 1.0)).eval(), false, C);
            viewer.data(uv_mesh_id).set_colors(C);
        }
        if (view_field == ViewingField::PHI) {
            const size_t nt = param->mesh().numTris();

            Eigen::MatrixXd C(nt, 3);
            for (size_t i = 0; i < nt; ++i) {
                double double_phi = std::fmod(2.0 * Eigen::Rotation2D<Real>(param->getU(i) * param->getR(i)).angle(), 2 * M_PI);
                if (double_phi < 0) double_phi += 2 * M_PI;

                // std::cout << "double_phi: " << double_phi << std::endl;
                double s = 1.0, v = 1.0;
                double hue = double_phi * (180.0 / M_PI);
                igl::hsv_to_rgb(hue, s, v, C(i, 0), C(i, 1), C(i, 2));
            }

            // Manually map orientations to colors using the HSV wheel
            viewer.data(uv_mesh_id).set_colors(C);
        }

        data_updated = false;
    };

    std::thread optimization_thread;
    bool optimization_cancelled = false, optimization_running = false;
    auto arap_thread = [&] {
        optimization_running = true;
        while (!optimization_cancelled) {
            if (gui_alpha_min != param->alphaMin()) param->setAlphaMin(gui_alpha_min);
            if (gui_alpha_max != param->alphaMax()) param->setAlphaMax(gui_alpha_max);
            param->runIteration();
            uv = param->uv();
            data_updated = true;
            glfwPostEmptyEvent(); // Run another event loop iteration so the viewer redraws.
        }
        optimization_cancelled = false;
        optimization_running = false;
    };

    auto launch_optimization = [&] {
        if (!optimization_running) {
            if (optimization_thread.joinable())
                optimization_thread.join(); // shouldn't actually happen...
            optimization_cancelled = false;
            optimization_thread = std::thread(arap_thread);
        }
    };

    auto initialize_with_mesh = [&](std::vector<MeshIO::IOVertex > &vertices,
                                    std::vector<MeshIO::IOElement> &triangles) {
        // Cancel the old optimization and wait it to finish (avoid race conditions)
        {
            optimization_cancelled = true;
            if (optimization_thread.joinable())
                optimization_thread.join();
        }

        {
            auto mesh = std::make_shared<parametrization::Mesh>(triangles, vertices);
            uv = parametrization::lscm(*mesh);
            param = std::make_unique<LGParam>(mesh, uv);
        }

        meshio_to_igl(vertices, triangles, V, F);
        igl::barycenter(V, F, barycenters3D);

        // Libigl requires clearing old meshes before setting new ones (of different sizes).
        viewer.data(surface_mesh_id).clear();
        viewer.data(     uv_mesh_id).clear();

        viewer.data(surface_mesh_id).set_mesh(V, F);

        param->setAlphaMin(gui_alpha_min);
        param->setAlphaMax(gui_alpha_max);
        global_scale = 0.5 * igl::avg_edge_length(V, F);

        update_data();

        for (int id = 0; id < 2; ++id) {
            viewer.data(id).face_based = true;
            viewer.data(id).set_colors(Eigen::RowVector3d(0.65, 0.65, 0.65));
        }
    };

    auto load_mesh = [&](const std::string &path) {
        std::vector<MeshIO::IOVertex > vertices;
        std::vector<MeshIO::IOElement> triangles;
        MeshIO::load(path, vertices, triangles);

        initialize_with_mesh(vertices, triangles);
    };

    // Load a surface mesh if one was specified; otherwise generate a paraboloid surface.
    double k1 = 1.0, k2 = 1.0, triArea = 0.01;

    auto generate_paraboloid = [&]() {
        std::vector<MeshIO::IOVertex > vertices;
        std::vector<MeshIO::IOElement> triangles;
        paraboloid(triArea, k1, k2, vertices, triangles);
        initialize_with_mesh(vertices, triangles);
    };

    if (argc > 1) load_mesh(argv[1]);
    else { generate_paraboloid(); }

    int left_view, right_view;
    viewer.callback_init = [&](Viewer &)
    {
        glfwSetWindowTitle(viewer.window, "Local-Global Parametrizer");
        left_view = viewer.core_list[0].id;
        right_view = viewer.append_core(Eigen::Vector4f(0, 0, width, height)); // will be resized by callback_post_resize
        viewer.core( left_view).background_color << 0.95, 0.95, 0.95, 1.0;
        viewer.core(right_view).background_color << 0.92, 0.92, 0.92, 1.0;

        viewer.core( left_view).rotation_type = ViewerCore::ROTATION_TYPE_TRACKBALL;
        viewer.core(right_view).rotation_type = ViewerCore::ROTATION_TYPE_TRACKBALL;

        viewer.core(left_view).camera_dnear = viewer.core(right_view).camera_dnear = 0.005;
        viewer.core(left_view).camera_dfar  = viewer.core(right_view).camera_dfar  = 50;

        viewer.data(surface_mesh_id).set_visible(false,  left_view);
        viewer.data(     uv_mesh_id).set_visible(false, right_view);

        // Initialize the split views' viewports
        viewer.callback_post_resize(viewer, width, height);

        return false; // also init the plugins
    };

    // Update viewew data if new ARAP iterations have been run.
    viewer.callback_pre_draw = [&](Viewer &/* v */) {
        if (data_updated) update_data();
        return false;
    };

    viewer.callback_key_pressed = [&](Viewer &, unsigned int key, int /* mod */)
    {
        if ((key == 'g') || (key == 'G')) {
            launch_optimization();
            return true;
        }

        if ((key == 'c') || (key == 'C')) {
            optimization_cancelled = true;
            return true;
        }

        return false;
    };

    viewer.callback_post_resize = [&](Viewer &v, int w, int h) {
        // v.core( left_view).viewport = Eigen::Vector4f(0, h / 2, w, h - (h / 2));
        // v.core(right_view).viewport = Eigen::Vector4f(0, 0, w, h / 2);
        v.core( left_view).viewport = Eigen::Vector4f(0, 0, w / 2, h);
        v.core(right_view).viewport = Eigen::Vector4f(w / 2, 0, w / 2, h);
        return true;
    };

    ////////////////////////////////////////////////////////////////////////////
    // IMGui UI
    ////////////////////////////////////////////////////////////////////////////
    igl::opengl::glfw::imgui::ImGuiMenu menu;
    viewer.plugins.push_back(&menu);

    menu.callback_draw_viewer_window = [&]() { };

    // Draw additional windows
    menu.callback_draw_custom_window = [&]() { // Define next window position + size
        ImGui::SetNextWindowPos(ImVec2(10, 10),    ImGuiCond_FirstUseEver);
        ImGui::SetNextWindowSize(ImVec2(200, 200), ImGuiCond_FirstUseEver);
        ImGui::Begin("Parametrization", nullptr,         ImGuiWindowFlags_NoSavedSettings);
        if (optimization_running) ImGui::StyleColorsLight();
        else                      ImGui::StyleColorsDark();

        ImGui::PushItemWidth(-80);

        const double alpha_min_min = 1.0;
        const double alpha_min_max = gui_alpha_max;
        const double alpha_max_min = gui_alpha_min;
        const double alpha_max_max = M_PI / 2.0;
        if (ImGui::DragScalar("alpha_min",  ImGuiDataType_Double, &gui_alpha_min,  0.01, &alpha_min_min, &alpha_min_max, "%.4f")) {
            // To avoid a race condition, only update the parametrizer's alpha if the optimization isn't running.
            // If the optimization is running, the optimization thread will read the updated alpha at the next iteration.
            if (!optimization_running) param->setAlphaMin(gui_alpha_min);
        }
        if (ImGui::DragScalar("alpha_max",  ImGuiDataType_Double, &gui_alpha_max,  0.01, &alpha_max_min, &alpha_max_max, "%.4f")) {
            // To avoid a race condition, only update the parametrizer's alpha if the optimization isn't running.
            // If the optimization is running, the optimization thread will read the updated alpha at the next iteration.
            if (!optimization_running) param->setAlphaMax(gui_alpha_max);
        }

        if (optimization_running) {
            if (ImGui::Button("Pause", ImVec2(-1,0)))
                optimization_cancelled = true;
        }
        else {
            if (ImGui::Button("Run", ImVec2(-1,0)))
                launch_optimization();
        }

        ImGui::Checkbox("Show frame", &show_frame);
        ImGui::Combo("Scalar field visualization", (int *)(&view_field), "None\0Alpha\0Phi\0\0");

        const double kmin = -2, kmax = 2, triAreaMin = 1e-5;
        ImGui::DragScalar("k1",  ImGuiDataType_Double, &k1,  0.1, &kmin, &kmax, "%.4f");
        ImGui::DragScalar("k2",  ImGuiDataType_Double, &k2,  0.1, &kmin, &kmax, "%.4f");
        ImGui::InputScalar("tri area",  ImGuiDataType_Double, &triArea, &triAreaMin, 0, "%.7f");
        if (ImGui::Button("Generate Paraboloid", ImVec2(-1,0))) generate_paraboloid();

        ImGui::Text("Energy: %f", param->energy());
        ImGui::Text("Num flips: %i", int(param->numFlips()));

        ImGui::PopItemWidth();
        ImGui::End();
    };

    launch_optimization(); // Start with ARAP running...

    viewer.launch(/* resizeable = */ true, /* fullscreen = */ false, width, height);

    optimization_cancelled = true;
    if (optimization_thread.joinable())
        optimization_thread.join();

    return EXIT_SUCCESS;
}

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import SmoothBivariateSpline
from scipy.interpolate import griddata
from scipy.interpolate import RectBivariateSpline
from scipy.interpolate import UnivariateSpline

import ndsplines
import os

def get_stiffness_coefficient_splines(x_scale_factors, y_scale_factors, stiffness_coefficients, inequality_constraints):
    points = np.concatenate((x_scale_factors.reshape(-1, 1), y_scale_factors.reshape(-1, 1)), axis = 1)

    # plot the feasible region
    d = np.linspace(0.8,1.8,1000)
    mesh_x,mesh_y = np.meshgrid(d,d)
    inside_point = ((inequality_constraints[0][0] * mesh_x + inequality_constraints[0][1] * mesh_y + inequality_constraints[0][2]<=0) & 
                    (inequality_constraints[1][0] * mesh_x + inequality_constraints[1][1] * mesh_y + inequality_constraints[1][2]<=0) &
                    (inequality_constraints[2][0] * mesh_x + inequality_constraints[2][1] * mesh_y + inequality_constraints[2][2]<=0) &
                    (inequality_constraints[3][0] * mesh_x + inequality_constraints[3][1] * mesh_y + inequality_constraints[3][2]<=0) & 
                    (inequality_constraints[4][0] * mesh_x + inequality_constraints[4][1] * mesh_y + inequality_constraints[4][2]<=0)).astype(int)
    
    grid_x = mesh_x[np.where(inside_point == 1)]
    grid_y = mesh_y[np.where(inside_point == 1)]
    sc_griddata = [griddata(points, stiffness_coefficients[:, i], (grid_x, grid_y), method='cubic') for i in range(5)]
    nan_check = np.sum([np.sum(np.isnan(sc_griddata[i])) for i in range(5)])
    if (nan_check > 0):
        print("Has NAN! Update inequality constraints!")
        return None 
    
    coefficient_splines = [SmoothBivariateSpline(points[:, 0], points[:, 1], stiffness_coefficients[:, i], s = len(points) * 100) for i in range(5)]

    coefficient_gradient_alpha = [smth_func.partial_derivative(1, 0) for smth_func in coefficient_splines]
    coefficient_gradient_beta = [smth_func.partial_derivative(0, 1) for smth_func in coefficient_splines]
    coefficient_hessian_alpha = [smth_func.partial_derivative(2, 0) for smth_func in coefficient_splines]
    coefficient_hessian_beta = [smth_func.partial_derivative(0, 2) for smth_func in coefficient_splines]
    coefficient_hessian_alpha_beta = [smth_func.partial_derivative(1, 1) for smth_func in coefficient_splines]

    splines = coefficient_splines + coefficient_gradient_alpha + coefficient_gradient_beta + coefficient_hessian_alpha + coefficient_hessian_beta + coefficient_hessian_alpha_beta

    # spline_functions = [(lambda capture_spline: lambda alpha, beta: capture_spline([alpha], [beta], grid = False)[0] if np.all(np.array(inequality_constraints) @ np.array([alpha, beta, 1.]) <= 0) else 0)(spline) for spline in splines]

    spline_functions = [(lambda capture_spline: lambda alphas, betas: capture_spline(alphas, betas, grid = False))(spline) for spline in splines]


    # for spline in splines:
    #     def func(alpha, beta):
    #         return spline([alpha], [beta], grid = False)[0]
    #     spline_functions.append(func)

    return spline_functions

def scipy_get_mat_params_over_one_pattern_params_grid_interpolation(grid_pattern, grid_data, smoothing = None):
    objective_splines = [UnivariateSpline(grid_pattern, grid_data[i], k = 5, s = smoothing) for i in range(grid_data.shape[0])]

    splines = []
    for spline in objective_splines:
        # Objective
        splines.append([spline])
        # Gradient
        splines.append([spline.derivative(n=1)])
        # Hessian
        splines.append([spline.derivative(n=2)])

    spline_functions = [(lambda capture_spline: lambda p: np.concatenate([unit_spline(p) for unit_spline in capture_spline]).flatten())(spline) for spline in splines]
    return spline_functions

def scipy_get_mat_params_over_pattern_params_grid_interpolation(grid_pattern_1, grid_pattern_2, grid_data):
    objective_splines = [RectBivariateSpline(grid_pattern_1, grid_pattern_2, grid_data[i], kx = 5, ky = 5) for i in range(grid_data.shape[0])]

    splines = []
    for spline in objective_splines:
        # Objective
        splines.append([spline])
        # Gradient
        splines.append([spline.partial_derivative(1, 0), spline.partial_derivative(0, 1)])
        # Hessian
        splines.append([spline.partial_derivative(2, 0),
                        spline.partial_derivative(1, 1),
                        spline.partial_derivative(1, 1),
                        spline.partial_derivative(0, 2)])

    # spline_functions = [(lambda capture_spline: lambda alpha, beta: capture_spline([alpha], [beta], grid = False)[0] if np.all(np.array(inequality_constraints) @ np.array([alpha, beta, 1.]) <= 0) else 0)(spline) for spline in splines]

    spline_functions = [(lambda capture_spline: lambda p: np.concatenate([unit_spline(p[0], p[1], grid = False) for unit_spline in capture_spline]).flatten())(spline) for spline in splines]
    return spline_functions

def ndsplines_get_mat_params_over_pattern_params_grid_interpolation(grid_data, *grid_patterns, spline_degree = 3):
    if (len(grid_patterns) != len(grid_data.shape) - 1):
        # The grid data should be gather over a grid over the number of pattern parameters for each of the material parameters.
        print("Wrong number of grid pattern parameters!")
        return None

    num_pattern_params = len(grid_patterns)
    
    def build_bspline(matIdx):
        grid_points = np.meshgrid(*grid_patterns, indexing='ij')
        x = np.stack(grid_points, axis=-1)
        y = grid_data[matIdx]

        bspline = ndsplines.make_interp_spline(x, y, degrees = [spline_degree for i in range(num_pattern_params)])
        return bspline
    
    objective_splines = [build_bspline(i) for i in range(len(grid_data))]

    splines = []

    for spline in objective_splines:
        splines.append([spline])
        # Gradient
        grad_splines = []
        for i in range(num_pattern_params):
            grad_splines.append(spline.derivative(dim = i, nu = 1))
        splines.append(grad_splines)
        # Hessian
        hess_splines = []
        for i in range(num_pattern_params):
            for j in range(num_pattern_params):
                hess_splines.append(grad_splines[i].derivative(dim = j, nu = 1))
        splines.append(hess_splines)
    spline_functions = [(lambda capture_spline: lambda p: np.concatenate([unit_spline(np.array(p).transpose()) for unit_spline in capture_spline]).flatten())(spline) for spline in splines]
    return spline_functions

def visualize_scale_factors(lines, alphas, betas, path = None):
    eqns = lines.copy()
    fig, ax = plt.subplots()
    x_coords = np.linspace(min(min(alphas), min(betas)), max(max(alphas), max(betas)), 100)
    for i in range(len(eqns) - 1):
        y = (-eqns[i][2] - x_coords * eqns[i][0]) / eqns[i][1]
        plt.plot(x_coords, y, label="inequality {}".format(i), alpha = 0.2)
    ax.set_aspect('equal')

    # now plot both limits against eachother
    plt.ylim(min(betas) - 0.1,max(betas) + 0.1)
    plt.xlim(min(alphas) - 0.1,max(alphas) + 0.1)

    d = np.linspace(0.8,1.6,600)
    x,y = np.meshgrid(d,d)
    inside_point = (eqns[0][0] * x + eqns[0][1] * y + eqns[0][2]<=0)
    for i in range(len(eqns)):
        inside_point = inside_point & (eqns[i][0] * x + eqns[i][1] * y + eqns[i][2]<=0)
    inside_point = inside_point.astype(int)
    plt.imshow(inside_point, 
                    extent=(x.min(),x.max(),y.min(),y.max()),origin="lower", cmap="Greys", alpha = 0.1);
    # Plot x = y line
    lims = [
    np.min([ax.get_xlim(), ax.get_ylim()]),  # min of both axes
    np.max([ax.get_xlim(), ax.get_ylim()]),  # max of both axes
    ]

    plt.scatter(alphas, betas)


    # now plot both limits against eachother
    ax.plot(lims, lims, 'k-', alpha=0.75, zorder=0)
    ax.set_aspect('equal')
    plt.xlabel(r'$x$')
    plt.ylabel(r'$y$')
    if path is not None:
        plt.savefig(path, dpi = 100)

def get_distance_to_line_segments(p, polyline):
    a = polyline[:-1]
    b = polyline[1:]
    v1 = p - a
    v2 = b - a
    dot = np.sum(v1 * v2, axis=1)
    len_sq = np.sum(v2 * v2, axis=1)
    epsilon = 1e-10  # or any small number you prefer
    param = np.where(len_sq > epsilon, dot / len_sq, -1)
    param = np.clip(param, 0, 1)
    closest_pts = a + np.expand_dims(param, axis=1) * v2
    d = p - closest_pts
    distances = np.linalg.norm(d, axis=1)
    return np.min(distances)

def get_stretch_angle_offset_from_pattern_params(scale_factor_splines, query_pattern_params):
    x_scale_factors = scale_factor_splines[0](query_pattern_params)
    y_scale_factors = scale_factor_splines[3](query_pattern_params)

    result = np.zeros_like(x_scale_factors)
    for i in range(len(x_scale_factors)):
        if (x_scale_factors[i] > y_scale_factors[i]):
            result[i] = 0.0
        else:
            result[i] = np.pi / 2.0
    return result


def compute_barycentric_coordinates(triangle, points):
    # Triangle vertices
    v0, v1, v2 = triangle[0], triangle[1], triangle[2]
    
    # Area of the triangle
    e1 = v1 - v0
    e2 = v2 - v0
    area = np.abs(0.5 * np.cross(e1, e2))
    if (area < 1e-8):
        return None
    
    # Initialize array for barycentric coordinates
    barycentric_coordinates = np.zeros((points.shape[0], 3))
    
    for i, p in enumerate(points):
        # Area of the triangle formed by p, v1 and v2
        e1 = v1 - p
        e2 = v2 - p
        area1 = np.abs(0.5 * np.cross(e1, e2))
        
        # Area of the triangle formed by p, v0 and v2
        e1 = v0 - p
        e2 = v2 - p
        area2 = np.abs(0.5 * np.cross(e1, e2))
        
        # Barycentric coordinates
        barycentric_coordinates[i, 0] = area1 / area
        barycentric_coordinates[i, 1] = area2 / area
        barycentric_coordinates[i, 2] = 1 - barycentric_coordinates[i, 0] - barycentric_coordinates[i, 1]
    
    return barycentric_coordinates

import shapely
def convert_shapely_objs_to_edge_soup(results):
    if isinstance(results, shapely.geometry.LineString):
        polyline = np.array(results.coords)
        polyline = np.concatenate((polyline, np.zeros((len(polyline), 1))), axis = 1)
        # edge soup should be a list of pairs of points
        edgeSoup = [polyline[i:i+2] for i in range(len(polyline)-1)]
    elif isinstance(results, shapely.geometry.MultiLineString):
        polylines = [np.array(line.coords) for line in results.geoms]
        edgeSoup = []
        for polyline in polylines:
            polyline = np.concatenate((polyline, np.zeros((len(polyline), 1))), axis = 1)
            edgeSoup.extend([polyline[i:i+2] for i in range(len(polyline)-1)])
    elif isinstance(results, shapely.geometry.Point):
        coords = np.array(results.coords)
        coords = np.concatenate((coords, np.zeros((len(coords), 1))), axis = 1)[0]
        edgeSoup = [[coords, coords]]
    elif isinstance(results, shapely.geometry.GeometryCollection):
        for obj in results.geoms:
            print(obj)
        edgeSoup = []
        for obj in results.geoms:
            edgeSoup.extend(convert_shapely_objs_to_edge_soup(obj))
    else:
        raise ValueError('Unexpected intersection result type: {}'.format(type(results)))
    return edgeSoup


# Assume the fusing_curve_polyline function return a list of polylines, where each polylines is a list of points that form a connected curve.
def get_pattern_polyline_function(fusing_curve_polyline):
    def pattern_polyline_function(theta_gamma, patternParams):
        default = [[[1, 0, 0], [0, 1, 0]], [[0, 1, 0], [0, 0, 1]], [[1, 0, 0], [0, 0, 1]]]
        default = []

        clip_tri = shapely.geometry.Polygon(theta_gamma)
        polylines = fusing_curve_polyline(patternParams)
        edge_soup = []
        for polyline in polylines:
            line = shapely.geometry.LineString(polyline)
            if line.intersects(clip_tri):
                result = line.intersection(clip_tri)
                edge_soup.extend(convert_shapely_objs_to_edge_soup(result))
               
        edge_soup = np.array(edge_soup)
        flat_edge_soup = edge_soup.reshape((-1, 3))[:, :2]
        barycentric_coords_soup = compute_barycentric_coordinates(np.array(theta_gamma), flat_edge_soup)
        if (barycentric_coords_soup is None):            
            return []
        return default + list(barycentric_coords_soup.reshape((-1, 2, 3)))
    return pattern_polyline_function

def get_pattern_function(fusing_curve_polyline):
    def pattern_function(theta, gamma, patternParams, margin, draw_boundary = True):
        if draw_boundary:
    #         This is for debugging only and shouldn't be used for generating the inflatable mesh.
            if (theta < 0.1):
                return - margin
            if (gamma < 0.1):
                return - margin
        # Gamma is y, theta is x
        # Gamma theta are between 0 and pi
        polylines = fusing_curve_polyline(patternParams)
        polyline_dist = np.inf
        for polyline in polylines:
            polyline_dist = min(polyline_dist, get_distance_to_line_segments(np.array([theta, gamma]), polyline))
        return polyline_dist - margin
    return pattern_function


def form_polylines(edges):
    polylines = []
    
    while edges:
        # Start with the first edge
        polyline = [edges.pop(0)]
        
        while edges:
            # Current edge
            current_edge = polyline[-1]
            
            # Find the next edge
            for i, edge in enumerate(edges):
                if current_edge[1] == edge[0]:
                    next_edge = edge
                    break
                elif current_edge[1] == edge[1]:
                    next_edge = edge[::-1]
                    break
            else:
                # No next edge found, start a new polyline
                break
            
            # Remove the next edge from the list and append it to the polyline
            edges.pop(i)
            polyline.append(next_edge)

        # Extend the polyline in the other direction if there are still edges left 
        while edges:
            # Current edge
            current_edge = polyline[0]
            
            # Find the next edge
            for i, edge in enumerate(edges):
                if current_edge[0] == edge[1]:
                    next_edge = edge
                    break
                elif current_edge[0] == edge[0]:
                    next_edge = edge[::-1]
                    break
            else:
                # No next edge found, start a new polyline
                break
            
            # Remove the next edge from the list and append it to the polyline
            edges.pop(i)
            polyline.insert(0, next_edge)
        
        polylines.append(polyline)
    
    return polylines


def remove_dangling_vertices(v, f, fusing):
    # Find all unique vertices used in the faces
    used_vertices = np.unique(f)

    # Create a mapping from old vertex indices to new ones
    old_to_new = np.full(v.shape[0], -1)
    old_to_new[used_vertices] = np.arange(used_vertices.shape[0])

    # Remove unused vertices from v
    new_v = v[used_vertices]

    # Update the faces to use the new vertex indices
    new_f = old_to_new[f]

    # Update the fusing array to use the new vertex indices
    new_fusing = fusing[used_vertices]

    return new_v, new_f, new_fusing

import numpy as np
from scipy.spatial import cKDTree

def remove_duplicates(points, edges, epsilon, end_points_only = False):
    if len(edges) == 0:
        print("No edges to remove duplicates from!")
        return points, edges
    if (len(points) == 0):
        print("No points to remove duplicates from!")
        return points, edges
    
    if not end_points_only:
        tree = cKDTree(points)
        not_duplicate = np.ones(len(points), dtype=bool)
        old_to_new = np.empty(len(points), dtype=int)
        next_index = 0
        for i in range(len(points)):
            if not_duplicate[i]:
                duplicates = tree.query_ball_point(points[i], epsilon)
                duplicates.remove(i)
                not_duplicate[duplicates] = False
                old_to_new[i] = next_index
                old_to_new[duplicates] = next_index
                next_index += 1
    else:
        # Remove duplicate points only among the valence 1 vertices.
        # This is to ensure that the interior vertices are not removed.
        # First reorder the points and edges to have first non valence 1 vertices followed by valence 1 vertices
        # Calculate the valence of each vertex
        valence = np.zeros(len(points), dtype=int)
        for edge in edges:
            valence[edge[0]] += 1
            valence[edge[1]] += 1

        # Separate the vertices into valence 1 and non-valence 1
        valence_1_vertices = np.where(valence == 1)[0]
        valence_non_1_vertices = np.where(valence != 1)[0]

        # Reorder the points and edges to have non-valence 1 vertices first
        new_order = np.concatenate((valence_non_1_vertices, valence_1_vertices))
        points = points[new_order]

        # Create a mapping from old indices to new indices
        old_to_new = np.empty(len(points), dtype=int)
        old_to_new[new_order] = np.arange(len(points))

        # Update the edges to use the new indices
        edges = old_to_new[edges]
        # Create a KDTree for efficient nearest neighbor search
        tree = cKDTree(points[len(valence_non_1_vertices):])
        old_to_new = np.empty(len(points), dtype=int)
        old_to_new[:len(valence_non_1_vertices)] = np.arange(len(valence_non_1_vertices))
        # Remove duplicate points among the valence 1 vertices
        not_duplicate = np.ones(len(points), dtype=bool)
        next_index = len(valence_non_1_vertices)  # Start the index from the first valence 1 vertex
        for i in range(len(valence_non_1_vertices), len(points)):
            if not_duplicate[i]:
                duplicates = tree.query_ball_point(points[i], epsilon)
                duplicates.remove(i - len(valence_non_1_vertices))
                duplicates = [d + len(valence_non_1_vertices) for d in duplicates]
                not_duplicate[duplicates] = False
                old_to_new[i] = next_index
                old_to_new[duplicates] = next_index
                next_index += 1

    new_points = points[not_duplicate]
    new_edges = old_to_new[edges]

    if len(new_edges) != 0:
        new_edges = new_edges[new_edges[:, 0] != new_edges[:, 1]]
    # Remove duplicate edges
    edge_set = set()
    new_edges_list = []
    for edge in new_edges:
        edge_tuple = tuple(sorted(edge))
        if edge_tuple not in edge_set:
            edge_set.add(edge_tuple)
            new_edges_list.append(edge)
    new_edges = np.array(new_edges_list)

    return new_points, new_edges

import svgwrite

def save_to_svg(boundary, polylines, filename):
    dwg = svgwrite.Drawing(filename, profile='tiny')

    # Add boundary to the SVG
    dwg.add(dwg.polyline(points=boundary, fill='none', stroke='black'))

    # Add polylines to the SVG
    for polyline in polylines:
        dwg.add(dwg.polyline(points=polyline, fill='none', stroke='black'))

    # Save the SVG
    dwg.save()

def save_to_obj(boundary, polylines, filename, use_3d = False):
    with open(filename, 'w') as f:
        # Initialize the index counter
        index = 1
        # Initialize the list of lines
        lines = []
        # Write the boundary to the .obj file and accumulate the indices
        boundary_indices = []
        for point in boundary:
            if use_3d:
                f.write(f'v {point[0]} {point[1]} {point[2]}\n')
            else:
                f.write(f'v {point[0]} {point[1]} 0\n')
            boundary_indices.append(index)
            index += 1
        for i in range(len(boundary_indices)):
            lines.append([boundary_indices[i], boundary_indices[(i + 1) % len(boundary_indices)]])

        # Write the polylines to the .obj file and accumulate the indices
        for polyline in polylines:
            polyline_indices = []
            for point in polyline:
                if use_3d:
                    f.write(f'v {point[0]} {point[1]} {point[2]}\n')
                else:
                    f.write(f'v {point[0]} {point[1]} 0\n')
                polyline_indices.append(index)
                index += 1
            for i in range(len(polyline_indices) - 1):
                lines.append([polyline_indices[i], polyline_indices[i + 1]])

        # Write the lines to the .obj file
        for line in lines:
            f.write(f'l {" ".join(str(i) for i in line)}\n')

def save_line_segments_as_obj(pts, edges, path):
    with open(path, 'w') as f:
        for point in pts:
            f.write(f'v {point[0]} {point[1]} {point[2]}\n')
        for edge in edges:
            f.write(f'l {edge[0] + 1} {edge[1] + 1}\n')

# import json
# # Save to json file
# with open('sheet_pattern.json', 'w') as fp:
#     data = {}
#     data['boundary'] = boundaryVxs[boundaryEdges[:, 0]][:, :2].tolist()
#     data['polylines'] = []
#     for polyline in sheet_edges_polylines:
#         polyline = np.array(polyline)
#         data['polylines'].append(sheet_vxs[np.array(list(polyline[:, 0]) + list([polyline[-1, 1]]))][:, :2].tolist())
#     json.dump(data, fp, indent = 4)
    
# Old code for clipping polylines
# import shapely

# from shapely.geometry import Polygon, LineString
# from shapely.ops import unary_union

# # Create the boundary polygon
# boundary_polygon = Polygon(boundaryVxs[boundaryEdges[:, 0]][:, :2])

# # Create a list to store the clipped edges
# clipped_edge_soup = []

# # Iterate over the polylines in sheet_edges_polylines
# for polyline in sheet_edges_polylines:
#     # Create a LineString from the polyline
#     line = LineString(sheet_vxs[np.array(list(np.array(polyline)[:, 0]) + [polyline[-1][1]])][:, :2])
#     # Check if the line intersects the boundary polygon
#     if line.intersects(boundary_polygon):
#         # Compute the intersection of the line with the boundary polygon
#         intersection = line.intersection(boundary_polygon)
#         clipped_edge_soup.extend(convert_shapely_objs_to_edge_soup(intersection))
# # Convert the list of clipped edges to a numpy array
# clipped_edge_soup = np.array(clipped_edges)

# clipped_pts = clipped_edge_soup.reshape((-1, 3))

# clipped_edges = np.array([[2 * i, 2 * i + 1] for i in range(len(clipped_edge_soup))])

# clipped_SV, clipped_SE = remove_duplicates(clipped_pts, clipped_edges, 1e-7)

# clipped_sheet_vxs = clipped_SV
# clipped_sheet_edges = clipped_SE
# clipped_sheet_edges_polylines = form_polylines(clipped_sheet_edges.tolist())
    
def get_mesh_boundary_info(m):
    boundaryLoop = m.boundaryElements()
    vertices = m.vertices()

    boundaryVxsIdxs = m.boundaryVertices()
    boundaryVxs = vertices[boundaryVxsIdxs]
    bdry_vxs_index_map = np.zeros(len(vertices), dtype = np.int64)
    for i, idx in enumerate(boundaryVxsIdxs):
        bdry_vxs_index_map[idx] = i
    boundaryEdges = bdry_vxs_index_map[boundaryLoop]
    boundaryEdges = form_polylines(boundaryEdges.tolist())
    return boundaryVxs, boundaryEdges

def get_polylines_from_edge_soup(edge_soup, duplicates_removable_threshold):
    new_edge_soup = np.array(edge_soup)
    pts = new_edge_soup.reshape((-1, 3))
    edges = np.array([[2 * i, 2 * i + 1] for i in range(len(new_edge_soup))])

    if duplicates_removable_threshold is not None:
        SV, SE = remove_duplicates(pts, edges, duplicates_removable_threshold[0])
        for threshold in duplicates_removable_threshold[1:-1]:
            SV, SE = remove_duplicates(SV, SE, threshold)
        SV, SE = remove_duplicates(SV, SE, duplicates_removable_threshold[-1], end_points_only = True)
        return SV, SE, form_polylines(SE.tolist())
    else:
        return pts, edges, form_polylines(edges.tolist())

def get_mesh_interior_polyline(m, fusing_data):
    fused_edge_soup = []
    boundaryVxsIdxs = m.boundaryVertices()
    boundary_vx_data = np.zeros(m.numVertices(), dtype = bool)
    boundary_vx_data[boundaryVxsIdxs] = True

    def get_fused_edge_soup(edge, i):
        v1, v2 = edge
        if ((fusing_data[v1] and not boundary_vx_data[v1]) and (fusing_data[v2] and not boundary_vx_data[v2])):
            fused_edge_soup.append([v1, v2])
    m.visitEdges(get_fused_edge_soup)
    return form_polylines(fused_edge_soup)

def get_fusing_lines_from_mesh(m, fusing_data):
    interior_fusing_lines = get_mesh_interior_polyline(m, fusing_data)
    boundary_fusing_lines = form_polylines(m.boundaryElements().tolist())
    # The curve smoothing class that uses this output assume for closed curves the first and last vertices are repeated.
    processed_interior_polylines = []
    for polyline in interior_fusing_lines:
        polyline = np.array(polyline)
        single_curve = polyline[:, 0]
        single_curve = np.concatenate((single_curve, [polyline[-1][1]]))
        processed_interior_polylines.append(single_curve)

    processed_boundary_polylines = []
    for polyline in boundary_fusing_lines:
        polyline = np.array(polyline)
        single_curve = polyline[:, 0]
        single_curve = np.concatenate((single_curve, [polyline[-1][1]]))
        processed_boundary_polylines.append(single_curve)

    # reorder the boundary polylines by the length of the polylines
    processed_boundary_polylines = sorted(processed_boundary_polylines, key = lambda x: len(x), reverse = True)
    return interior_fusing_lines, boundary_fusing_lines, processed_interior_polylines + processed_boundary_polylines[1:], processed_boundary_polylines[:1]

def is_counter_clockwise(boundary_curve):
    signed_area = 0.0
    x1, y1 = boundary_curve[0][:2]
    for i in range(1, len(boundary_curve)):
        x2, y2 = boundary_curve[i][:2]
        signed_area += (x1*y2 - x2*y1)
        x1, y1 = x2, y2
    return signed_area >= 0

import sys
sys.path.append("../")
import wall_generation
import MeshFEM
import pickle

def save_base_data(path, upsampledMesh, boundary_upsampledMesh, upsampledAngles, upsampledPatternParams, rparam, boundaryVxs, boundaryEdges):
    upsampledMesh.save(path + "/upsampledMesh.obj")
    boundary_upsampledMesh.save(path + "/boundary_upsampleMesh.obj")
    np.save(path + "/upsampledAngles.npy", upsampledAngles)
    np.save(path + "/upsampledPatternParams.npy", upsampledPatternParams)
    np.save(path + "/rparam_uv.npy", rparam.uv())
    np.save(path + "/boundaryVxs.npy", boundaryVxs)
    with open(path + "/boundaryEdges.pkl", 'wb') as f:
        pickle.dump(boundaryEdges, f)

def load_base_data(path):
    boundary_upsampledMesh = MeshFEM.mesh.Mesh(path + "/boundary_upsampleMesh.obj")
    upsampledMesh = MeshFEM.mesh.Mesh(path + "/upsampledMesh.obj")
    boundary_upsampleMesh_vertices = boundary_upsampledMesh.vertices()
    boundary_upsampleMesh_vertices = np.concatenate((boundary_upsampleMesh_vertices, np.zeros((len(boundary_upsampleMesh_vertices), 1))), axis = 1)
    boundary_upsampleMesh_triangles = boundary_upsampledMesh.triangles()
    upsampleMesh_vertices = upsampledMesh.vertices()
    upsampleMesh_vertices = np.concatenate((upsampleMesh_vertices, np.zeros((len(upsampleMesh_vertices), 1))), axis = 1)
    upsampleMesh_triangles = upsampledMesh.triangles()
    upsampledAngles = np.load(path + "/upsampledAngles.npy")
    upsampledPatternParams = np.load(path + "/upsampledPatternParams.npy")
    boundaryVxs = np.load(path + "/boundaryVxs.npy")
    # load boundaryEdges using pickle
    with open(path + "/boundaryEdges.pkl", 'rb') as f:
        boundaryEdges = pickle.load(f)
    return boundary_upsampledMesh, upsampledMesh, boundary_upsampleMesh_vertices, boundary_upsampleMesh_triangles, upsampleMesh_vertices, upsampleMesh_triangles, upsampledAngles, upsampledPatternParams, boundaryVxs, boundaryEdges

def save_sdf_data(path, sdfVertices, sdfTris, sdf, sheet_vxs, sheet_edges, sheet_edges_polylines, concatenated_polylines):
    np.save(path + "/sdfVertices.npy", sdfVertices)
    np.save(path + "/sdfTris.npy", sdfTris)
    np.save(path + "/sdf.npy", sdf)
    np.save(path + "/sheet_vxs.npy", sheet_vxs)
    np.save(path + "/sheet_edges.npy", sheet_edges)
    with open(path + "/sheet_edges_polylines.pkl", 'wb') as f:
        pickle.dump(sheet_edges_polylines, f)
    np.save(path + "/concatenated_polylines.npy", concatenated_polylines)

def load_sdf_data(path):
    sdfVertices = np.load(path + "sdfVertices.npy")
    sdfTris = np.load(path + "sdfTris.npy")
    sdf = np.load(path + "sdf.npy")
    sheet_vxs = np.load(path + "/sheet_vxs.npy")
    sheet_edges = np.load(path + "/sheet_edges.npy")
    # load sheet_edges_polylines using pickle
    with open(path + "/sheet_edges_polylines.pkl", 'rb') as f:
        sheet_edges_polylines = pickle.load(f)
    concatenated_polylines = np.load(path + "/concatenated_polylines.npy")
    return sdfVertices, sdfTris, sdf, sheet_vxs, concatenated_polylines, sheet_edges_polylines

# load_data = None, "param", "final_results"
def get_polyline_from_pattern_parameters(rparam, fusing_curve_polyline, nsubdiv = 7, frequency = 0.1, duplicates_removable_threshold = [1e-4, 1e-2, 1e-1, 2e0, 3e0], path = None, load_data = None):
    if load_data not in [None, "param", "final_results"]:
        print("Wrong load data type!")
        return None
    
    if load_data is None:
        boundary_upsampledMesh, upsampledAngles, upsampledPatternParams = rparam.upsampledVertexLeftStretchAnglesAndPatternParameters(2)
        boundary_upsampleMesh_vertices = boundary_upsampledMesh.vertices()
        boundary_upsampleMesh_triangles = boundary_upsampledMesh.triangles()

        nsubdiv=nsubdiv
        upsampledMesh, upsampledAngles, upsampledPatternParams = rparam.upsampledVertexLeftStretchAnglesAndPatternParameters(nsubdiv)
        upsampleMesh_vertices = upsampledMesh.vertices()
        upsampleMesh_triangles = upsampledMesh.triangles()

        #### Get boundary edges for meshing
        flatten_mesh = MeshFEM.mesh.Mesh(boundary_upsampleMesh_vertices, boundary_upsampleMesh_triangles)
        boundaryVxs, boundaryEdges = get_mesh_boundary_info(flatten_mesh)
            
        if path is not None:
            save_base_data(path, upsampledMesh, boundary_upsampledMesh, upsampledAngles, upsampledPatternParams, rparam, boundaryVxs, boundaryEdges)

    else:
        boundary_upsampledMesh, upsampledMesh, boundary_upsampleMesh_vertices, boundary_upsampleMesh_triangles, upsampleMesh_vertices, upsampleMesh_triangles, upsampledAngles, upsampledPatternParams, boundaryVxs, boundaryEdges = load_base_data(path)
        if load_data == "final_results":
            sdfVertices, sdfTris, sdf, sheet_vxs, concatenated_polylines, sheet_edges_polylines = load_sdf_data(path)
            return sdfVertices, sdfTris, sdf, sheet_vxs, concatenated_polylines, sheet_edges_polylines,  boundaryVxs, boundaryEdges, upsampleMesh_vertices, upsampleMesh_triangles, upsampledAngles, upsampledPatternParams

    pattern_function = get_pattern_function(fusing_curve_polyline)
    pattern_polyline_function = get_pattern_polyline_function(fusing_curve_polyline)
    import time
    start_time = time.time()
    (sdfVertices, sdfTris, sdf, edge_soup) = wall_generation.evaluate_cross_field_custom_pattern(upsampleMesh_vertices, upsampleMesh_triangles, upsampledAngles, upsampledPatternParams, pattern_function, pattern_polyline_function, frequency=frequency, margin = 0.0, nsubdiv = 0)
    print("Compute polyline field takes: ", time.time() - start_time)

    start_time = time.time()

    sheet_vxs, sheet_edges, sheet_edges_polylines = get_polylines_from_edge_soup(edge_soup, duplicates_removable_threshold)
    concatenated_polylines = np.concatenate(sheet_edges_polylines)

    print("Merge edge soup takes: ", time.time() - start_time)
    if path is not None:
        save_sdf_data(path, sdfVertices, sdfTris, sdf, sheet_vxs, sheet_edges, sheet_edges_polylines, concatenated_polylines)

    return sdfVertices, sdfTris, sdf, sheet_vxs, concatenated_polylines, sheet_edges_polylines,  boundaryVxs, boundaryEdges, upsampleMesh_vertices, upsampleMesh_triangles, upsampledAngles, upsampledPatternParams

def get_fabrication_file_from_mesh(sheet, m, fusing_data, channelMargin, duplicates_removable_threshold, path, use_holes = False, use_obj = True):
    fused_edge_soup = []
    boundaryVxsIdxs = m.boundaryVertices()
    boundary_vx_data = np.zeros(m.numVertices(), dtype = bool)
    boundary_vx_data[boundaryVxsIdxs] = True

    def get_fused_edge_soup(edge, i):
        v1, v2 = edge
        if ((fusing_data[v1] and not boundary_vx_data[v1]) and (fusing_data[v2] and not boundary_vx_data[v2])):
            fused_edge_soup.append([v1, v2])
    m.visitEdges(get_fused_edge_soup)
    if fused_edge_soup == []:
        sheet_vxs, sheet_edges, sheet_edges_polylines = [], [], []
    else:
        sheet_vxs, sheet_edges, sheet_edges_polylines = get_polylines_from_edge_soup(m.vertices()[fused_edge_soup], duplicates_removable_threshold)

    if channelMargin > 0:
        boundaryVxs, boundaryEdges = get_mesh_boundary_info(m)
        V = m.vertices()
        polylines = sheet.fusedRegionBooleanIntersectSheetBoundary()
        shapely_boundaryEdges = shapely.MultiLineString([V[p] for p in polylines])
        bypasses = shapely_boundaryEdges.buffer(channelMargin)
        if bypasses.geom_type == 'Polygon': bypasses = [bypasses] # we generally expect a multipolygon...
        outerAirChannelPolygons = [shapely.ops.unary_union([shapely.Polygon(boundaryVxs[np.array(boundaryEdges[0])[:, 0]])] + list(bypasses.geoms))]
        smart_polygon = outerAirChannelPolygons[0]
        coords = np.array(smart_polygon.exterior.coords)
    else:
        boundaryVxs, boundaryEdges = get_mesh_boundary_info(m)

    polylines = []
    for polyline in sheet_edges_polylines:
        polyline = np.array(polyline)
        polylines.append(sheet_vxs[np.array(list(polyline[:, 0]) + list([polyline[-1, 1]]))][:, :2].tolist())
    
    if channelMargin > 0:
        if use_obj:
            save_to_obj(coords[:-1], polylines, path)
        else:
            save_to_svg(coords[:-1], polylines, path)
    else:
        for be in boundaryEdges[1:]:
            be = np.array(be)
            polylines.append(boundaryVxs[np.array(list(be[:, 0]) + list([be[-1, 1]]))][:, :2].tolist())

        if use_obj:
            save_to_obj(boundaryVxs[np.array(boundaryEdges[0])[:, 0]][:, :2], polylines, path)
        else:
            save_to_svg(boundaryVxs[np.array(boundaryEdges[0])[:, 0]][:, :2], polylines, path)
    
    if channelMargin > 0:
        final_vertices = np.concatenate((sheet_vxs, coords))
        concatenated_polylines = list(np.concatenate(sheet_edges_polylines))
        for i in range(coords.shape[0]):
            concatenated_polylines.append([len(sheet_vxs) + i, len(sheet_vxs) + (i + 1) % (coords.shape[0])])
        return final_vertices, concatenated_polylines
    else:
        if len(sheet_edges_polylines) > 0:
            return sheet_vxs, np.concatenate(sheet_edges_polylines)
        return [], []


def get_non_manifold_boundary_vertices(m):
    edge_count = np.zeros(m.numVertices())
    for edge in m.boundaryElements():
        p1 = edge[0]
        p2 = edge[1]
        if p1 == p2:
            print("warning!")
            continue
        edge_count[p1] += 1
        edge_count[p2] += 1
    nonmanifold_vxs = []
    for i in range(m.numVertices()):
        if edge_count[i] > 2:
            nonmanifold_vxs.append(i)

    return nonmanifold_vxs


def export_top_bottom_mesh(isheet, export_path, shape_name, pattern_name):
    mesh_3d = isheet.visualizationMesh(True)
    mesh_2d = isheet.mesh()
    vx_3d = mesh_3d.vertices()
    elements_3d = mesh_3d.elements()

    new_mesh_3d = MeshFEM.Mesh(vx_3d[:mesh_2d.numVertices()], elements_3d[:mesh_2d.numElements()])

    new_mesh_3d.save(export_path + '{}_{}_mesh_3d_top.obj'.format(shape_name, pattern_name))

    new_mesh_3d = MeshFEM.Mesh(vx_3d[mesh_2d.numVertices():], elements_3d[mesh_2d.numElements():] - mesh_2d.numVertices())
    new_mesh_3d.save(export_path + '{}_{}_mesh_3d_bottom.obj'.format(shape_name, pattern_name))

from scipy.interpolate import splprep, splev
def smooth_polyline(polyline, s=10.0, segment_length=0.1):
    # Separate the polyline into x and y coordinates
    x, y, _ = polyline.T

    # Calculate the total length of the polyline
    total_length = np.sum(np.sqrt(np.diff(x)**2 + np.diff(y)**2))

    # Calculate the number of points to use
    num_points = int(np.ceil(total_length / segment_length))

    # Fit a spline to the polyline
    tck, u = splprep([x, y], s=s, per=1)

    # Generate new points on the spline
    unew = np.linspace(0, 1.0, num_points)
    out = splev(unew, tck)

    # Create a 2D array of zeros with the same shape as out
    zeros = np.zeros_like(out[0])

    # Stack the smoothed polyline with the zeros
    return np.column_stack((out[0], out[1], zeros))
    
import numpy as np
from scipy.linalg import eig

import numpy as np
from scipy.linalg import eig

def find_bounding_box(points):
    # Calculate the centroid of the points
    center = points.mean(axis=0)

    # Center the points
    C = points - center

    # Calculate the covariance matrix
    cov = np.cov(C.T)

    # Find the eigenvalues and eigenvectors
    eigenvalues, eigenvectors = eig(cov)

    # The principal axis is the eigenvector corresponding to the largest eigenvalue
    principal_axis = eigenvectors[:, np.argmax(eigenvalues)]

    # Project the points onto the principal axis
    projected_points = C @ principal_axis

    # Find the minimum and maximum coordinates of the projected points
    min_coord, max_coord = np.min(projected_points), np.max(projected_points)

    # The semi-major axis of the bounding box is half the length of the box
    semi_major_axis = (max_coord - min_coord) / 2

    # The minor axis is the range of the points projected onto the direction orthogonal to the principal axis
    orthogonal_direction = np.array([-principal_axis[1], principal_axis[0]])
    projected_points_orthogonal = C @ orthogonal_direction
    min_coord_orthogonal, max_coord_orthogonal = np.min(projected_points_orthogonal), np.max(projected_points_orthogonal)
    semi_minor_axis = (max_coord_orthogonal - min_coord_orthogonal) / 2

    return center, semi_major_axis, semi_minor_axis, principal_axis

import os
script_dir = os.path.dirname(os.path.realpath(__file__))

sys.path.append(os.path.join(script_dir, '../'))

import inflatables_parametrization as parametrization

def ellipsify_holes(boundaryVxs, boundaryEdges, sheet_vxs, sheet_edges_polylines, smoothing = 10.0, area_threshold = 1.0, avg_len = 4, duplicates_removable_threshold = [], distance_threshod = 1e-2):
    selected_elements = [np.array(sublist)[:, 0] for sublist in boundaryEdges[1:]]
    holes_vxs = list(boundaryVxs[selected_elements])
    removed_holes_vxs = []

    # Gather all boundaryEdges
    all_boundary_edges = []
    for sublist in boundaryEdges:
        all_boundary_edges.extend(sublist)
    all_boundary_edges = np.array(all_boundary_edges)
        
    mesh_vertices = boundaryVxs[:, :2]
    mesh_edges = all_boundary_edges

    # Construct an AABB tree for the edges
    tree = igl.AABB_f64_2()
    tree.init(mesh_vertices, mesh_edges)
    # Find the closest point on the edges for each vertex

    print("number of input holes: ", len(sheet_edges_polylines))
    hole_marker = np.zeros(len(sheet_edges_polylines))

    from scipy import optimize
    from scipy.special import ellipe
    from scipy.special import ellipeinc

    def get_ellipse_using_box(center, semi_major_axis, semi_minor_axis, principal_axis):

        iw =semi_minor_axis
        ih = semi_major_axis  

        # Center of the ellipse
        h = center[0]
        k = center[1]

        # Number of points in the fusing line
        num_points = 15

        # Eccentricity of the ellipse
        e = np.sqrt(1 - (min(iw, ih) / max(iw, ih))**2)

        # Total arc length of the ellipse
        total_arclength = 4 * max(iw, ih) * ellipe(e**2)

        num_points = min(num_points, int(total_arclength / avg_len))
        print("num_points: ", num_points)

        # Compute evenly spaced arc lengths
        arc_lengths = np.linspace(0, total_arclength, num_points)

        # Function to compute the difference between the target and actual arc length
        def func(t, target_arclength, e):
            return target_arclength - max(iw, ih) * ellipeinc(t, e**2)

        # Compute the angles for the points
        t = np.array([optimize.root(func, [0], args=(l, e)).x[0] for l in arc_lengths])
        # Compute the points on the ellipse
        x_unrotated = iw * np.cos(t)
        y_unrotated = ih * np.sin(t)

        theta = np.arctan2(principal_axis[1], principal_axis[0]) + np.pi / 2

        # Rotate the ellipse by angle
        x = h + (x_unrotated) * np.cos(theta) - (y_unrotated) * np.sin(theta)
        y = k + (x_unrotated) * np.sin(theta) + (y_unrotated) * np.cos(theta)
        return np.column_stack((x, y, np.zeros((len(x), 1))))[:-1]
    
    def get_ellipse(points):
        # Find the best oriented bounding ellipse for these points
        center, semi_major_axis, semi_minor_axis, principal_axis = find_bounding_box(points[:, :2])

        return get_ellipse_using_box(center, semi_major_axis, semi_minor_axis, principal_axis)

    
    def merge_polylines(open_polylines, duplicates_removable_threshold, end_points_only = False):
        # Turn holes_vxs into edge soup
        edge_soup = []
        for hole in open_polylines:
            for i in range(len(hole)):
                edge_soup.append((hole[i], hole[(i+1)%len(hole)]))
        edge_soup = np.array(edge_soup)

        pts = edge_soup.reshape((-1, 3))

        edges = np.array([[2 * i, 2 * i + 1] for i in range(len(edge_soup))])
        print(len(pts))
        print(len(pts), len(edges))

        SV, SE = remove_duplicates(pts, edges, duplicates_removable_threshold[0], end_points_only = end_points_only)
        for threshold in duplicates_removable_threshold[1:-1]:
            SV, SE = remove_duplicates(SV, SE, threshold, end_points_only = end_points_only)
            print(len(SV), len(SE))
        SV, SE = remove_duplicates(SV, SE, duplicates_removable_threshold[-1], end_points_only=True)
        print(len(SV), len(SE))
        edges_polylines = form_polylines(SE.tolist())
        selected_elements = [np.array(sublist)[:, 0] for sublist in edges_polylines]
        open_polylines = [SV[element] for element in selected_elements]
        print(len(open_polylines))
        return open_polylines
    


    threshold = 4
    open_polylines = []
    fusing_lines = []

    for index, polyline in enumerate(sheet_edges_polylines):
        polyline = np.array(polyline)

        curr_hole_point = np.array(sheet_vxs[polyline[:, 0]])[:, :2]

        indicator = tree.squared_distance(mesh_vertices, mesh_edges, curr_hole_point)
        adj_boundary = np.min(indicator) < distance_threshod
        closed_loop = (polyline[0, 0] == polyline[-1, 1]) 

        # Check if the polyline forms a closed loop
        if (not adj_boundary) and closed_loop:
            points = np.array(sheet_vxs[polyline[:, 0]])
            center, semi_major_axis, semi_minor_axis, principal_axis = find_bounding_box(points[:, :2])
            if (semi_minor_axis * semi_major_axis < area_threshold):
                fusing_lines.append(get_ellipse(points))
            else:
                holes_vxs.append(get_ellipse(points))

        elif (not adj_boundary) and (not closed_loop):
            open_polylines.append(np.array(sheet_vxs[polyline[:, 0]]))
        elif adj_boundary and (not closed_loop):
            points = np.array(sheet_vxs[polyline[:, 0]])
            if len(points) < threshold:
                fusing_lines.append(points)
                continue
            center, semi_major_axis, semi_minor_axis, principal_axis = find_bounding_box(points[:, :2])
            finished = False
            while not finished:
                semi_major_axis *= 0.8
                if semi_major_axis < semi_minor_axis:
                    finished = True
                elif (semi_minor_axis * semi_major_axis < area_threshold):
                    fusing_lines.append(get_ellipse(points))
                    finished = True
                else:
                    try:
                        print("Trying smaller ellipse")
                        ellipse = get_ellipse_using_box(center, semi_major_axis, semi_minor_axis, principal_axis)
                        indicator = tree.squared_distance(mesh_vertices, mesh_edges, ellipse)
                        adj_boundary = np.min(indicator) < distance_threshod
                        if not adj_boundary:
                            holes_vxs.append(ellipse)
                            print("Got new ellipse")
                            finished = True
                    except:
                        finished = True


    duplicates_removable_threshold = np.array([1e-4, 1e-2, 1e-1, 1e0])
    # open_polylines = merge_polylines(open_polylines, duplicates_removable_threshold)

    for points in open_polylines:
        if len(points) < threshold:
            fusing_lines.append(points)
            continue
        # Find the best oriented bounding ellipse for these points
        try:
            center, semi_major_axis, semi_minor_axis, principal_axis = find_bounding_box(points[:, :2])
            if (semi_major_axis < distance_threshod):
                # return a line segment along the major axis
                # Calculate the endpoints of the line segment
                start_point = center - semi_major_axis * principal_axis
                end_point = center + semi_major_axis * principal_axis

                # The line segment is then the pair of points (start_point, end_point)
                line_segment = np.array([start_point, end_point])
                fusing_lines.append(np.column_stack((line_segment, np.zeros((len(line_segment), 1)))))

            elif (semi_minor_axis * semi_major_axis < area_threshold):
                fusing_lines.append(get_ellipse(points))
            else:
                holes_vxs.append(get_ellipse(points))
        except:
            fusing_lines.append(points)

    #     # elif adj_boundary and closed_loop:
    # duplicates_removable_threshold = np.array([1e-4, 1e-2, 1e-1])
    # print("fusing_lines", fusing_lines)
    # if len(fusing_lines) > 10:
    #     fusing_lines = merge_polylines(fusing_lines, duplicates_removable_threshold)


    #     # else:
    boundary_holes_vxs = []
    boundary_polygon = shapely.geometry.Polygon(boundaryVxs[np.array(boundaryEdges[0])[:, 0]])


    print("number of output holes: ", len(boundary_holes_vxs) + len(holes_vxs))
    print("number of boundary holes: ", len(boundary_holes_vxs))
    print("number of non-boundary holes: ", len(holes_vxs))
    return fusing_lines, holes_vxs, boundary_polygon



from shapely import LineString
from scipy.spatial import ConvexHull, convex_hull_plot_2d
import igl
import numpy.linalg as la
def post_process_holes(boundaryVxs, boundaryEdges, sheet_vxs, sheet_edges_polylines, smoothing = 10.0, area_threshold = 1.0, avg_len = 4, duplicates_removable_threshold = [], distance_threshod = 1e-2):
    selected_elements = [np.array(sublist)[:, 0] for sublist in boundaryEdges[1:]]
    holes_vxs = list(boundaryVxs[selected_elements])
    removed_holes_vxs = []

    # Gather all boundaryEdges
    all_boundary_edges = []
    for sublist in boundaryEdges:
        all_boundary_edges.extend(sublist)
    all_boundary_edges = np.array(all_boundary_edges)
        
    mesh_vertices = boundaryVxs[:, :2]
    mesh_edges = all_boundary_edges

    # Construct an AABB tree for the edges
    tree = igl.AABB_f64_2()
    tree.init(mesh_vertices, mesh_edges)
    # Find the closest point on the edges for each vertex

    print("number of input holes: ", len(sheet_edges_polylines))
    hole_marker = np.zeros(len(sheet_edges_polylines))

    threshold = 4
    for index, polyline in enumerate(sheet_edges_polylines):
        polyline = np.array(polyline)
        
        if len(polyline) < threshold:
            removed_holes_vxs.append(np.array(sheet_vxs[polyline[:, 0]]))
            hole_marker[index] = 1
            continue

        curr_hole_point = np.array(sheet_vxs[polyline[:, 0]])[:, :2]

        indicator = tree.squared_distance(mesh_vertices, mesh_edges, curr_hole_point)
        adj_boundary = np.min(indicator) < distance_threshod
        closed_loop = (polyline[0, 0] == polyline[-1, 1]) 

        # Check if the polyline forms a closed loop
        if (not adj_boundary) and closed_loop:
            # If it's a closed loop, append it to holes_vxs
            holes_vxs.append(np.array(sheet_vxs[polyline[:, 0]]))
            hole_marker[index] = 2
        elif adj_boundary and (not closed_loop):
            # If it's not a closed loop, connect the end points with a half circular arc that points away from the polyline
            # Convert polyline to 2D points
            points = np.array(sheet_vxs[polyline[:, 0]])

            # Calculate the midpoint, direction, and length of the last segment
            p1 = points[0]
            p2 = points[-1]
            midpoint = (p1 + p2) / 2
            # Calculate the radius and center of the semi-circle
            radius = np.linalg.norm(p2 - p1) / 2
            center = midpoint

            # Calculate the direction from p1 to p2
            direction_vector = p2 - p1
            direction = np.arctan2(direction_vector[1], direction_vector[0])
            # Calculate the middle point of the points
            middle_point = points[len(points) // 2]

            # Calculate the vector from the first point to the middle point
            vector = middle_point - points[0]

            # Calculate the z-component of the cross product of the direction vector and the vector
            cross_product_z = direction_vector[0] * vector[1] - direction_vector[1] * vector[0]


            # Generate points on the semi-circle from p2 to p1
            if cross_product_z > 0:
                t = np.linspace(direction, direction - np.pi, 6)[1:-1]
            else:
                t = np.linspace(direction, direction + np.pi, 6)[1:-1]
            semi_circle_points = np.empty((4, 3))
            semi_circle_points[:, 0] = center[0] + radius * np.cos(t)
            semi_circle_points[:, 1] = center[1] + radius * np.sin(t)
            semi_circle_points[:, 2] = 0
            # Append the polyline points and the semi-circle points to holes_vxs
            final_points = np.concatenate([points, semi_circle_points])
            new_points, _ = remove_duplicates(final_points, [], 1e0)
            if (len(new_points) < threshold):
                removed_holes_vxs.append(np.array(sheet_vxs[polyline[:, 0]]))
                hole_marker[index] = 1
                continue
            hull = ConvexHull(new_points[:, :2])  # We only consider the first two columns (x, y coordinates)
            
            hull_points = np.concatenate((new_points[hull.vertices], np.zeros((len(hull.vertices), 1))), axis = 1)

            indicator = tree.squared_distance(mesh_vertices, mesh_edges, hull_points[:, :2])
            adj_boundary = np.min(indicator) < distance_threshod
            if adj_boundary:
                removed_holes_vxs.append(np.array(sheet_vxs[polyline[:, 0]]))
                hole_marker[index] = -1
                continue
            
            holes_vxs.append(hull_points[:, :3])
            hole_marker[index] = 2

        elif adj_boundary and closed_loop:
            removed_holes_vxs.append(np.array(sheet_vxs[polyline[:, 0]]))
            hole_marker[index] = -1
            continue
        else:
            points = np.array(sheet_vxs[polyline[:, 0]])
            new_points, _ = remove_duplicates(points, [], 1e0)
            if (len(new_points) < threshold):
                removed_holes_vxs.append(np.array(sheet_vxs[polyline[:, 0]]))
                hole_marker[index] = 1
                continue
            hull = ConvexHull(new_points[:, :2])  # We only consider the first two columns (x, y coordinates)
            hull_points = np.concatenate((new_points[hull.vertices], np.zeros((len(hull.vertices), 1))), axis = 1)
            holes_vxs.append(hull_points[:, :3])
            hole_marker[index] = 2

    # if smoothing != 0:
    #     holes_vxs = [smooth_polyline(hole, smoothing, avg_len) for hole in holes_vxs if len(hole) > threshold]
    #     short_holes_vxs = [hole for hole in holes_vxs if len(hole) <= threshold]
    #     removed_holes_vxs.extend(short_holes_vxs)

    # def calculate_area(polyline):
    #     x, y, _ = polyline.T
    #     return 0.5 * np.abs(np.dot(x, np.roll(y, 1)) - np.dot(y, np.roll(x, 1)))

    # # Calculate the area of each hole and filter out the ones below the threshold
    # holes_vxs = [hole[:-1] for hole in holes_vxs if calculate_area(hole) > area_threshold]

    # small_holes_vxs = [hole[:-1] for hole in holes_vxs if calculate_area(hole) <= area_threshold]
    # removed_holes_vxs.extend(small_holes_vxs)

    # # filtered_fusing_vxs =  [hole[:-1] for hole in smoothed_holes_vxs if calculate_area(hole) <= area_threshold ]


    # def check_self_intersection_and_replace_with_hull(polyline):
    #     curr_points = np.array(polyline)[:, :2]
    #     line = LineString(curr_points)

    #     # Check if the line has self-intersections
    #     if not line.is_simple:
    #         hull = ConvexHull(curr_points)  # We only consider the first two columns (x, y coordinates)
    #         hull_points = np.concatenate((curr_points[hull.vertices], np.zeros((len(hull.vertices), 1))), axis = 1)
    #         return hull_points
    #     else:
    #         return polyline

    # # Apply the function to each hole in filtered_holes_vxs
    # holes_vxs = [check_self_intersection_and_replace_with_hull(hole) for hole in holes_vxs]

    # if len(duplicates_removable_threshold) > 1:
    #     # Turn holes_vxs into edge soup
    #     edge_soup = []
    #     for hole in holes_vxs:
    #         for i in range(len(hole)):
    #             edge_soup.append((hole[i], hole[(i+1)%len(hole)]))
    #     edge_soup = np.array(edge_soup)

    #     pts = edge_soup.reshape((-1, 3))

    #     edges = np.array([[2 * i, 2 * i + 1] for i in range(len(edge_soup))])
    #     print(len(pts))

    #     SV, SE = remove_duplicates(pts, edges, duplicates_removable_threshold[0])
    #     for threshold in duplicates_removable_threshold[1:-1]:
    #         SV, SE = remove_duplicates(SV, SE, threshold)
    #     SV, SE = remove_duplicates(SV, SE, duplicates_removable_threshold[-1], end_points_only=True)
    #     print(len(SV))
    #     edges_polylines = form_polylines(SE.tolist())
    #     selected_elements = [np.array(sublist)[:, 0] for sublist in edges_polylines]
    #     holes_vxs = [SV[element] for element in selected_elements]

    # Separate holes into the ones that intersect the boundary vs the ones that don't
    boundary_holes_vxs = []
    non_boundary_holes_vxs = []
    boundary_polygon = shapely.geometry.Polygon(boundaryVxs[np.array(boundaryEdges[0])[:, 0]])
    for hole in holes_vxs:
        hole = np.array(hole)
        # check whether the hole intersects the boundary using shapely
        hole_polygon = shapely.geometry.Polygon(hole[:, :2])

        if not boundary_polygon.contains(hole_polygon):
            boundary_holes_vxs.append(hole)
        else:
            non_boundary_holes_vxs.append(hole)


    print("number of output holes: ", len(boundary_holes_vxs) + len(non_boundary_holes_vxs))
    print("number of boundary holes: ", len(boundary_holes_vxs))
    print("number of non-boundary holes: ", len(non_boundary_holes_vxs))
    return boundary_holes_vxs, non_boundary_holes_vxs, boundary_polygon, removed_holes_vxs, hole_marker

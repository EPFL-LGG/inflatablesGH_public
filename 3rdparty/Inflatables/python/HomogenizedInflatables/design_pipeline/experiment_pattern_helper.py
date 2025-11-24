import numpy as np 

import os 
import sys
# Get the directory of the current script
script_dir = os.path.dirname(os.path.realpath(__file__))

# Append the parent directory to the system path
sys.path.append(os.path.join(script_dir, '..'))
from experiment_configuration import *

sys.path.append(os.path.join(script_dir, '../../..'))

import inflation
import inflatables_parametrization as parametrization

def cosine_fusing_curve_polyline(patternParams):
#     Draw cosine curves.
    amp = patternParams[0]
    def get_y_from_x(x):
        return amp * np.cos(x) * 0.5 * np.pi + np.pi / 2

    x_coords = np.linspace(-np.pi, np.pi, 15)
    y_coords = get_y_from_x(x_coords)
    x_coords += np.pi
    x_coords /= 2
    polyline = np.concatenate(((y_coords).reshape(-1, 1), (x_coords).reshape(-1, 1)), axis = 1)
    return [polyline]

def dash_line_fusing_curve_polyline(patternParams):
#     Draw dash_line
    radius = patternParams[0]
    angle = patternParams[1]
    
    dash_point = np.array([np.cos(angle / 180 * np.pi), np.sin(angle / 180 * np.pi)]) * radius + np.array([0, 0])
    mid_point = np.array([0, 0])
    dash_line = (np.array([dash_point, mid_point * 2 - dash_point]) + [2.5, 2.5]) / 5 * np.pi
    return [dash_line]

mid_point = np.array([0, 5 / 4])
from shapely.geometry import LineString


def line_segment_intersection(segment1, segment2):
    line1 = LineString(segment1)
    line2 = LineString(segment2)
    intersection = line1.intersection(line2)

    if intersection.is_empty:
        return None
    else:
        return intersection.x, intersection.y
    
def double_zigzag_fusing_curve_polyline(patternParams):
#     Draw dash_line
    r = patternParams[0]
    angle = patternParams[1]
    
    dash_point = np.array([np.cos(angle / 180 * np.pi), np.sin(angle / 180 * np.pi)]) * r + mid_point

    dash1 = np.array([dash_point, mid_point * 2 - dash_point])
    
    dash2 = np.array(dash1)
    dash2[:, 1] = - dash2[:, 1]

    # Check if the dashes intersect
    intersection = line_segment_intersection(dash1, dash2)
    if intersection is not None:
        # If they intersect, clip them at the intersection point
        dash1[1] = intersection
        dash2[1] = intersection
        
    return np.array([(dash1 + np.array([2.5, 2.5])) / 5 * np.pi, (dash2 + np.array([2.5, 2.5])) / 5 * np.pi])

def elliptic_hole_fusing_curve_polyline(patternParams):
    # Draw dash_line
    angle = patternParams[0]
    theta = np.radians(angle)  # convert angle to radians

    iw = patternParams[1]
    ih = patternParams[2]
    
    # Center of the ellipse
    h = k = 0

    t, area = parametrization.ellipsePointParameters(0.3, ih, iw)
    t = list(t) + [t[0]]
    # Compute the points on the ellipse
    x_unrotated = iw * np.cos(t)
    y_unrotated = ih * np.sin(t)

    # Rotate the ellipse by angle
    x = h + (x_unrotated) * np.cos(theta) - (y_unrotated) * np.sin(theta)
    y = k + (x_unrotated) * np.sin(theta) + (y_unrotated) * np.cos(theta)

    # Combine x and y into an array of points
    points = np.column_stack((x, y))
    return [(points + np.array([2.5, 2.5])) / 5 * np.pi]


Pattern_data = [
    {
        'name': 'dashed_line',
        'experiment_file': '{}/dashed_line/experiment_result.json'.format(pattern_data_base_folder),
        'stiffness_path': '{}/dashed_line'.format(pattern_data_base_folder),
        'label': 'Dashed line',
        'num_pattern_params': 2,
        'param_index': [0, 1],
        'default_param': [1, 50],
        'param_range': [[0.5, 2.5], [45, 75]],
        'param_normalization_factor': [2.5 - 0.5, 30],
        'fusing_curve_polyline_function': dash_line_fusing_curve_polyline,
        'use_holes': False
    },
    {
        'name': 'cosine_curve_amplitude_full_period',
        'experiment_file': '{}/cosine_curve_amplitude_full_period/experiment_result.json'.format(pattern_data_base_folder),
        'stiffness_path': '{}/cosine_curve_amplitude_full_period/'.format(pattern_data_base_folder),
        'label': 'Cosine curve',
        'num_pattern_params': 1,
        'param_index': [0],
        'default_param': [0.2],
        'param_range': [[0.03, 0.4]],
        'param_normalization_factor': [1],
        'fusing_curve_polyline_function': cosine_fusing_curve_polyline,
        'use_holes': False
    },
    {
        'name': 'square_with_ellipse_hole_angle_width_height',
        'label': 'Elliptic holes',
        'experiment_file': '{}/square_with_ellipse_hole_angle_width_height/experiment_result.json'.format(pattern_data_base_folder),
        'stiffness_path': '{}/square_with_ellipse_hole_angle_width_height'.format(pattern_data_base_folder),
        'num_pattern_params': 3,
        'param_index': [0, 1, 2],
        'default_param': [36, 0.6, 1.65],
        'param_range': [[6, 45], [0.2, 0.9], [1.3, 1.8]],
        'param_normalization_factor': [45 - 6, 0.9 - 0.2, 1.8 - 1.3],
        'fusing_curve_polyline_function': elliptic_hole_fusing_curve_polyline,
        'use_holes': True
    },
    {
        'name': 'voronoi_5',
        'label': 'Random Voronoi diagrams',
        'experiment_file': '{}/voronoi_5/experiment_result.json'.format(pattern_data_base_folder),
        'stiffness_path': '{}/voronoi_5'.format(pattern_data_base_folder),
        'num_pattern_params': 0,
        'param_index': [],
        'default_param': [],
        'param_range': [],
        'param_normalization_factor': [],
        'fusing_curve_polyline_function': None,
        'use_holes': False
    },
]

Shape_data = [
    {
        'name': 'taller_hill',
        'path': '{}/examples/taller_hill.obj'.format(shape_data_base_folder),
    },
    {
        'name': 'igloo',
        'path': '{}/examples/igloo.obj'.format(shape_data_base_folder),
    },
    {
        'name': 'vest_coarse',
        'path': '{}/examples/vest_coarse.obj'.format(shape_data_base_folder),
    },
    {
        'name': 'new_neck_brace_coarse',
        'path': '{}/examples/new_neck_brace_coarse.obj'.format(shape_data_base_folder),
    },
    {
        'name': 'cashew_planar_coarse',
        'path': '{}/examples/cashew_planar_coarse.obj'.format(shape_data_base_folder),
    },
    {
        'name': 'squiward',
        'path': '{}/examples/squidward_remesh.obj'.format(shape_data_base_folder),
    },
    {
        'name': 'two_rings',
        'path': '{}/examples/two_rings.obj'.format(shape_data_base_folder),
    },
    {
        'name': 'lemonade_stand_modular',
        'path': '{}/examples/lemonade_stand_modular.obj'.format(shape_data_base_folder),
    },
    {
        'name': 'helmet_with_vents_coarse',
        'path': '{}/examples/helmet_with_vents_coarse.obj'.format(shape_data_base_folder),
    },
    {
        'name': 'half_sphere_iso',
        'path': '{}/examples/half_sphere_iso.obj'.format(shape_data_base_folder),
    },
    {
        'name': 'lilium',
        'path': '{}/examples/lilium.msh'.format(shape_data_base_folder),
    },
    {
        'name': 'hive',
        'path': '{}/examples/hive.obj'.format(shape_data_base_folder),
    },
    {
        'name': 'mask',
        'path': '{}/python/SiggraphExamples/Meshes/mask.remeshed.obj'.format(shape_data_base_folder),
    },
    {
        'name': 'cast',
        'path': '{}/python/SiggraphExamples/Meshes/20200116_cast_velcro_strip.remeshed.obj'.format(shape_data_base_folder),
    },
    {
        'name': 'donut_pavilion',
        'path': '{}/examples/donut.obj'.format(shape_data_base_folder),
    },
    {
        'name': 'snail_shell',
        'path': '{}/examples/snail_shell_low_res.obj'.format(shape_data_base_folder),
    },
    {
        'name': 'cashew',
        'path': '{}/examples/cashew.obj'.format(shape_data_base_folder),
    },
]

def parse_input(shape_index, pattern_index):
    pattern = Pattern_data[pattern_index]

    experiment_file = pattern['experiment_file']
    stiffness_path = pattern['stiffness_path']
    pattern_name = pattern['name']
    num_pattern_params = pattern['num_pattern_params']
    param_index = pattern['param_index']
    default_param = pattern['default_param']
    param_range = pattern['param_range']
    param_normalization_factor = pattern['param_normalization_factor']
    fusing_curve_polyline = pattern['fusing_curve_polyline_function']
    use_holes = pattern['use_holes']
    shape = Shape_data[shape_index]

    shape_name = shape['name']
    shape_path = shape['path']
    return experiment_file, stiffness_path, pattern_name, num_pattern_params, param_index, default_param, param_range, param_normalization_factor, fusing_curve_polyline, shape_name, shape_path, use_holes
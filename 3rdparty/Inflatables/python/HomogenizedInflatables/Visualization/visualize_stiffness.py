import matplotlib.pyplot as plt
import numpy as np
import numpy.linalg as la
import os

change_ticks = True
def plot_kappa(data, name):
    fig, ax = plt.subplots()
    plt.plot(data)
    fig.set_size_inches(12, 4)
    ax.title.set_text("Kappa")
    fig.tight_layout()
    plt.savefig('kappa_{}.png'.format(name), dpi = 300)
    plt.show()

def plot_stiffness(data, name, merge_data = None, position = None, bending = True, sort_data = False):
    fig, ax = plt.subplots()
    if sort_data:
        data = sorted(data, key=np.min)
    plt.violinplot(data, showmedians=True)
    if (merge_data is not None):
        plt.violinplot(merge_data, positions = [position], showmedians=True)
    fig.set_size_inches(12, 4)
    title = "Bending stiffness values" if bending else "Stretching stiffness values"
    ax.title.set_text(title)

    if change_ticks:
        fig.canvas.draw()
        labels = [(item.get_text()) for item in ax.get_xticklabels()]
        select_theta = [int(x) for x in labels[1:-1]]
        def get_angle(i):
            return 80 / 49 * i + 10
        labels[1:-1] = np.round([get_angle(x) for x in select_theta], 2)
        ax.set_xticklabels(labels)

    fig.tight_layout()
    plt.savefig('{}_stiffness_values_{}.png'.format("bending" if bending else "stretching", name), dpi = 300)
    plt.show()

def plot_strain(stiffness_path, name, iteration_tags, high_pressure):
    data = []
    for tag in iteration_tags:
        data.append(np.load("{}/{}_strain_values_{}_{}.npy".format(stiffness_path, "high_pressure" if high_pressure else "low_pressure", name, tag), allow_pickle = True))
            
    fig, ax = plt.subplots()
    plt.violinplot(data, showmedians=True)
    fig.set_size_inches(12, 4)
    ax.title.set_text("Strain values")
    fig.tight_layout()
    plt.savefig('{}_strain_values_{}.png'.format("high_pressure" if high_pressure else "low_pressure", name), dpi = 300)
    plt.show()

def plot_scale_factors(data, name, merge_data = None, position = None, sort_data = False):
    fig, ax = plt.subplots()
    if sort_data:
        data = sorted(data, key=np.min)
    violin_parts = plt.violinplot(data)
    for pc in violin_parts['bodies']:
        pc.set_facecolor('white')
    if (merge_data is not None):
        merge_violin_parts = plt.violinplot(merge_data, positions = [position])
        for pc in merge_violin_parts['bodies']:
            pc.set_facecolor('white')
    fig.set_size_inches(12, 4)
    ax.title.set_text("Max and min scale factors")

    if change_ticks:
        fig.canvas.draw()
        labels = [(item.get_text()) for item in ax.get_xticklabels()]
        select_theta = [int(x) for x in labels[1:-1]]
        def get_angle(i):
            return 80 / 49 * i + 10
        labels[1:-1] = np.round([get_angle(x) for x in select_theta], 2)
        ax.set_xticklabels(labels)

    fig.tight_layout()
    plt.savefig('scale_factors_{}.png'.format(name), dpi = 300)
    plt.show()


def plot_extreme_stiffness_data(data, name, use_min=True, bending = True, sort_data = False):
    fig, ax = plt.subplots()
    if sort_data:
        data = sorted(data, key=np.min)
    plt.plot(data)
    # ax.set_yscale('log')
    fig.set_size_inches(12, 4)
    ax.title.set_text("{} stiffness values".format("Min" if use_min else "Max"))

    if change_ticks:
        fig.canvas.draw()
        labels = [(item.get_text()) for item in ax.get_xticklabels()]
        select_theta = [int(x) for x in labels[1:-1]]
        def get_angle(i):
            return 80 / 49 * i + 10
        labels[1:-1] = np.round([get_angle(x) for x in select_theta], 2)
        ax.set_xticklabels(labels)
    
    fig.tight_layout()
    plt.savefig('{}_{}_stiffness_values_{}.png'.format("min" if use_min else "max", "bending" if bending else "stretching", name), dpi = 300)
    plt.show()

    return data

def plot_all_data(kappa_path, stiffness_path, name, iteration_tags, merge_data_path = None, plot_data = True, sort_data = False):
    if (kappa_path is not None):
        data = []
        for tag in iteration_tags:
            path = "{}/kappa_{}_{}.npy".format(kappa_path, name, tag)
            if (os.path.isfile(path)):
                data.append(np.load(path, allow_pickle = True))
        plot_kappa(data, name)

    used_tags = []

    if (stiffness_path is not None):
        bending_stiffness_data = []

        for tag in iteration_tags:
            path = "{}/{}/bending_stiffness_values_{}_{}.npy".format(stiffness_path, tag, name, tag)
            if (os.path.isfile(path)):
                bending_stiffness = np.load(path, allow_pickle = True)
                bending_stiffness_data.append(bending_stiffness)
                used_tags.append(tag)
        
        if plot_data: plot_stiffness(bending_stiffness_data, name, bending = True, sort_data = sort_data)

        stretching_stiffness_data = []
        for tag in iteration_tags:
            path = "{}/{}/stretching_stiffness_{}_{}.npy".format(stiffness_path, tag, name, tag)
            if (os.path.isfile(path)):
                stretching_stiffness_data.append(np.load(path, allow_pickle = True))

        if plot_data: plot_stiffness(stretching_stiffness_data, name, bending = False, sort_data = sort_data)
            
        scale_factor_data = []

        for tag in iteration_tags:
            path = "{}/{}/scale_factors_{}_{}.npy".format(stiffness_path, tag, name, tag)
            if (os.path.isfile(path)):
                scale_factor_data.append(np.load(path, allow_pickle = True))
        if plot_data: 
            plot_scale_factors(scale_factor_data, name, sort_data = sort_data)

            plot_extreme_stiffness_data(np.max(bending_stiffness_data, axis = 1), name, use_min=False, bending = True, sort_data = sort_data)
            plot_extreme_stiffness_data(np.min(bending_stiffness_data, axis = 1), name, use_min=True, bending = True, sort_data = sort_data)

            plot_extreme_stiffness_data(np.max(stretching_stiffness_data, axis = 1), name, use_min=False, bending = False, sort_data = sort_data)
            plot_extreme_stiffness_data(np.min(stretching_stiffness_data, axis = 1), name, use_min=True, bending = False, sort_data = sort_data)

    return bending_stiffness_data, stretching_stiffness_data, scale_factor_data, used_tags


def get_axis_scale_factors(scale_factor_path, name, iteration_tags):
    x_scale_factors = []
    y_scale_factors = []

    for tag in iteration_tags:
        path = "{}/{}/average_deformation_gradient_matrix_{}_{}.npy".format(scale_factor_path, tag, name, tag)
        if (os.path.isfile(path)):
            deformation_gradient_matrix = np.load(path, allow_pickle = True)
            x_scale_factors.append(1. / la.norm(deformation_gradient_matrix @ np.array([1, 0])))
            y_scale_factors.append(1. / la.norm(deformation_gradient_matrix @ np.array([0, 1])))
    return np.array(x_scale_factors), np.array(y_scale_factors)
    
def get_max_flattening_factor_offset(scale_factor_path, name, iteration_tags):
    max_flattening_factor_offset = []
    for tag in iteration_tags:
        path = "{}/{}/average_deformation_gradient_matrix_{}_{}.npy".format(scale_factor_path, tag, name, tag)
        if (os.path.isfile(path)):
            deformation_gradient_matrix = np.load(path, allow_pickle = True)
            # this is only true if we are using patches with reflectional symmetry.
            if deformation_gradient_matrix[0, 0] < deformation_gradient_matrix[1, 1]:
                max_flattening_factor_offset.append(0)
            else:
                max_flattening_factor_offset.append(np.pi / 2)
    
    return max_flattening_factor_offset


def get_extreme_stiffness(stiffness_path, name, iteration_tags, use_min=True):
    data = []
    selection = np.min if use_min else np.max
    for tag in iteration_tags:
        path = "{}/{}/stiffness_values_{}_{}.npy".format(stiffness_path, tag, name, tag)
        if (os.path.isfile(path)):
            data.append(selection(np.load(path, allow_pickle = True)))
    return data

def get_stiffness_coefficients(stiffness_path, name, iteration_tags):
    data = []
    for tag in iteration_tags:
        if (os.path.isfile("{}/{}/stiffness_coefficient_{}_{}.npy".format(stiffness_path, tag, name, tag))):
            data.append(np.load("{}/{}/stiffness_coefficient_{}_{}.npy".format(stiffness_path, tag, name, tag), allow_pickle = True))
    return data

def plot_extreme_stiffness(stiffness_path, name, iteration_tags, use_min=True):
    data = get_extreme_stiffness(stiffness_path, name, iteration_tags, use_min)
    fig, ax = plt.subplots()
    plt.plot(data)
    ax.set_yscale('log')
    fig.set_size_inches(12, 4)
    ax.title.set_text("{} stiffness values".format("Min" if use_min else "Max"))

    if change_ticks:
        fig.canvas.draw()
        labels = [(item.get_text()) for item in ax.get_xticklabels()]
        select_theta = [int(x) for x in labels[1:-1]]
        def get_angle(i):
            return 80 / 49 * i + 10
        labels[1:-1] = np.round([get_angle(x) for x in select_theta], 2)
        ax.set_xticklabels(labels)
    
    fig.tight_layout()
    plt.savefig('{}_stiffness_values_{}.png'.format("min" if use_min else "max", name), dpi = 300)
    plt.show()

    return data

def compare_stiffness_across_pressure_values(stiffness_paths, name, iteration_tags, legends = None):
    all_stiffness_data = []
    for stiffness_path in stiffness_paths:        
        stiffness_data = []
        for tag in iteration_tags:
            stiffness_data.append(max(np.load("{}/stiffness_values_{}_{}.npy".format(stiffness_path, name, tag), allow_pickle = True)))
        all_stiffness_data.append(stiffness_data)

    fig, ax = plt.subplots()
    for i in range(len(stiffness_paths)):
        plt.plot(all_stiffness_data[i], label="Simulation {}".format(legends[i]))
    plt.legend(loc="upper right")

    fig.set_size_inches(12, 4)
    ax.title.set_text("Stiffness values")
    fig.tight_layout()
    plt.savefig('stiffness_values_comparison_{}_{}.png'.format(name, len(stiffness_paths)), dpi = 300)
    plt.show()

def plot_min_stretching_factor(scale_factor_paths, name, iteration_tags, legends = None):
    all_scale_factor_data = []
    for scale_factor_path in scale_factor_paths:        
        scale_factor_data = []
        for tag in iteration_tags:
            scale_factor_data.append(np.load("{}/scale_factors_{}_{}.npy".format(scale_factor_path, name, tag), allow_pickle = True)[0])
        all_scale_factor_data.append(scale_factor_data)

    theoretical_data = [(2 * (5 - w) / np.pi + w) / 5 for w in iteration_tags]
    fig, ax = plt.subplots()
    for i in range(len(scale_factor_paths)):
        plt.plot(all_scale_factor_data[i], label="Simulation {}".format(legends[i]), color = "tab:blue")
    plt.plot(theoretical_data, label="Analytical", color = "tab:orange")
    plt.legend(loc="upper left")

    fig.set_size_inches(12, 4)
    ax.title.set_text("Min scale factors")
    plt.xlabel("Fusing width")

    fig.tight_layout()
    plt.savefig('scale_factors_comparison_{}_{}.png'.format(name, len(scale_factor_paths)), dpi = 300)
    plt.show()


def analytical_zigzag_scale_factors(h1, h2, h3, parallel_tube_scale_factors = None, total_area = None):

    # parallel_tube_2_over_pi = 0.6670760889008647
    parallel_tube_2_over_pi = 2 / np.pi

    aux_len = np.sqrt((h3 / 4.0 )**2 + (h2 - h1) ** 2)
    cos_theta = (h2 - h1) / aux_len
    sin_theta = h3 / 4.0 / aux_len
    theta = np.arccos(cos_theta) / np.pi * 180

    if (parallel_tube_scale_factors is None):
        parallel_tube_scale_factors = (2 * parallel_tube_2_over_pi * h2 + h3 - 2 * h2) / h3
        temp = np.ones((len(parallel_tube_scale_factors), 2))
        temp[:, 0] = parallel_tube_scale_factors
        parallel_tube_scale_factors = temp

    scale_factor_1 = np.sqrt(parallel_tube_scale_factors[:, 1] ** 2 * cos_theta ** 2 + parallel_tube_scale_factors[:, 0] ** 2 * sin_theta ** 2)
    if (total_area is not None):
        scale_factor_2 = total_area / scale_factor_1
    else:
        scale_factor_2 = parallel_tube_scale_factors[:, 0] * parallel_tube_scale_factors[:, 1] / scale_factor_1

    return scale_factor_1, scale_factor_2, theta


def analytical_zigzag_scale_factors_angles(angles, parallel_tube_scale_factors, total_area = None):
    cos_theta = np.cos(angles)
    sin_theta = np.sin(angles)

    scale_factor_1 = np.sqrt(parallel_tube_scale_factors[:, 1] ** 2 * cos_theta ** 2 + parallel_tube_scale_factors[:, 0] ** 2 * sin_theta ** 2)
    if (total_area is not None):
        scale_factor_2 = total_area / scale_factor_1
    else:
        scale_factor_2 = parallel_tube_scale_factors[:, 0] * parallel_tube_scale_factors[:, 1] / scale_factor_1

    return scale_factor_1, scale_factor_2


def get_scale_factor(scale_factor_path, name, iteration_tags):
    scale_factor_data = []
    for tag in iteration_tags:
        path = "{}/{}/scale_factors_{}_{}.npy".format(scale_factor_path, tag, name, tag)
        if os.path.isfile(path):
            scale_factor_data.append(np.load(path, allow_pickle = True))
        else:
            print("Missing {}".format(tag))
    return scale_factor_data

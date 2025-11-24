import json
import matplotlib.pyplot as plt
import numpy as np


def plot_stretching_graph(values : float, filename : str, show_graph : bool):
    plot_min_r = 0
    plot_max_r = None
    r = list(values)

    betas = np.linspace(0, 2*np.pi, 1000)
    theta = list(betas)

    fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
    ax.plot(theta, r)
    ax.set_rmax(max(values) if plot_max_r is None else plot_max_r)
    ax.set_rmin(min(values) - 0.2 * (max(values) - min(values)) if plot_min_r is None else plot_min_r)
    # ax.set_rticks([0.5, 1, 1.5, 2])  # Less radial ticks
    # ax.set_rlabel_position(-22.5)  # Move radial labels away from plotted line
    ax.grid(True)

    ax.set_title("Stretching Stiffness", va='bottom')
    plt.tight_layout()
    if(show_graph): plt.show()
    plt.savefig(filename, dpi = 300)

def plot_bending_graph(sampled_stiffness : float, filename : str, show_graph : bool):
    plot_min_r = None
    plot_max_r = None
    sampled_alpha = np.linspace(0, np.pi, 1000)
    r = list(sampled_stiffness) + list(sampled_stiffness)
    theta = list(sampled_alpha) + list(np.pi + np.array(sampled_alpha))

    fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
    ax.plot(theta, r)
    ax.set_rmax(max(sampled_stiffness) if plot_max_r is None else plot_max_r)
    ax.set_rmin(min(sampled_stiffness) - 0.2 * (max(sampled_stiffness) - min(sampled_stiffness)) if plot_min_r is None else plot_min_r)
    # ax.set_rticks([0.5, 1, 1.5, 2])  # Less radial ticks
    # ax.set_rlabel_position(-22.5)  # Move radial labels away from plotted line
    ax.grid(True)

    ax.set_title("Bending Stiffness", va='bottom')
    plt.tight_layout()
    if(show_graph): plt.show()
    plt.savefig(filename, dpi = 300)
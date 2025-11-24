
import pandas as pd
import numpy as np

def get_largest_box_of_valid_data_point(data):
    grid_size = [0, 0, 0]
    for i, entry in enumerate(data['pattern_parameters']):
        grid_size[i] = len(entry['values'])

    df = pd.DataFrame(data['data'])
    param_grid = np.multiply(np.array(df['Planar equilibrium']).reshape(grid_size), np.array(df['Ipu simulation succeed']).reshape(grid_size))
    import time

    box_sizes = []
    box_entries = []
    curr_max = 0
    for i in range(grid_size[0]):
        for j in range(grid_size[1]):
            for k in range(grid_size[2]):
                for i1 in range(grid_size[0] + 1)[i+1:]:
                    for j1 in range(grid_size[1] + 1)[j+1:]:
                        for k1 in range(grid_size[2] + 1)[k+1:]:
                            if 0 in param_grid[i:i1, j:j1, k:k1]:
                                box_size = 0
                                break
                            else:
                                box_size = (i1 - i) * (j1 - j) * (k1 - k)
                                box_entries.append([i, i1, j, j1, k, k1])
                                box_sizes.append(box_size)
                                curr_max = max(curr_max, box_size)

    largest_box = box_entries[np.argmax(box_sizes)]

    box_sizes[np.argmax(box_sizes)]
    return largest_box
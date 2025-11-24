import experiment_helper
import igl
from periodic_simulation_setup import *
import json

allowBending = False
useTFT = True


name = 'zigzag_line'
time_stamp = time.strftime("%Y_%m_%d_%H_%M")
result_folder = 'output/{}/{}'.format(name, time_stamp)
if not os.path.exists(result_folder):
    os.makedirs(result_folder)  

# pressure = 0.8
stiffness_pressure = 0.4
scale_factor_pressure = 0.01

labels = np.arange(51)
for label in labels:
    label = 20
    with open('../data/ZigZag_Line/ZigZag_Line{}.json'.format(label), 'r') as f:
        data = json.load(f)
    fusedVertices = data['FusedVertices']
    fusedVertices = [True if vx == 1 else False for vx in fusedVertices]
    V = data['Vertices']
    F = data['Faces']
    m = MeshFEM.Mesh(V, F)
    for vx in m.boundaryVertices():
        fusedVertices[vx] = False
    ipu = inflation.InflatablePeriodicUnit(m, fusedVertices)

    finalMarkers = np.where(np.array(fusedVertices) == 1)[0]
    m, finalMarkers = periodic_unit_helper.shift_and_merge_2D_periodic_mesh(m, finalMarkers)
    m, finalMarkers = periodic_unit_helper.shift_and_merge_2D_periodic_mesh(m, finalMarkers, axis = 1)

    fusedVtx = get_fusedVtx_using_markers(len(m.vertices()), finalMarkers)
    ipu = inflation.InflatablePeriodicUnit(m, fusedVtx = fusedVtx, epsilon = 1e-9)

    experiment_helper.run_experiment(ipu, m, fusedVtx, stiffness_pressure, scale_factor_pressure, name, label, allowBending, result_folder, useTFT = useTFT)
    break
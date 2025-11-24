import sys; sys.path.append('../..')
import MeshFEM, mesh
import numpy as np
import reinflate

reinflate.inflation.set_max_num_tbb_threads(6)

m = mesh.Mesh('SiggraphExamples/Mask/high_res_fw1e-7/optDesign.msh', embeddingDimension=3)
iwv = np.loadtxt('SiggraphExamples/Mask/high_res_fw1e-7/iwv.txt.gz')
liftedSheetPositions = np.loadtxt('SiggraphExamples/Mask/high_res_fw1e-7/liftedSheetPositions.txt.gz') 

frameOutDir = 'test_frames'
reinflate.reinflate(frameOutDir, m, iwv, liftedSheetPositions)

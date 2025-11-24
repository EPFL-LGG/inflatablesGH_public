# Utilities for comparing the accuracy of ElasticSheet and ElasticSolid
# simulations of a simple semi-sphere(hemisphere) plate under different gravitites, mesh refinement, and at various
# thicknesses.

import os
import sys; sys.path.extend(['..', '../validations/'])
import numpy as np, time
import copy
import mesh, elastic_sheet, elastic_solid, energy, tensors
import meshing, triangulation, py_newton_optimizer, mesh_operations
import loads
from io_redirection import suppress_stdout as so
from sim_utils import getBBoxVars, BBoxFace
from tri_mesh_viewer import TriMeshViewer

import meshpy # We use `meshpy`'s tetgen wrapper
# import meshpy.triangle as meshtriangle
from meshpy.tet import MeshInfo, build, Options

# Test geometry: rectangular strip
def stripBoundary(L = 4, H = 1):
    pts = [[0, 0], [0, H], [L, H], [L, 0]]
    edges = [[0, 1], [1, 2], [2, 3], [3, 0]]
    return pts, edges

# Test geometry: creased rectangular strip
def creasedBoundary(L = 4, H = 1):
    pts = [[0, 0], [0, H], [L/2, H],[L, H], [L, 0], [L/2, 0]]
    edges = [[0, 1], [1, 3], [3, 4], [4, 0], [2, 5]]
    return pts, edges

def RotateMat3D(ux,uy,uz,theta):
    if ux*ux+uy*uy+uz*uz != 1:
        print('ERROR! (ux,uy,uz) do not represent a unit vector!')
        return -1
    else:
        c = np.cos(theta)
        s = np.sin(theta)
        C = 1 - c
        # construct Rotate Matrix using Axis-Angle representation
        R = np.zeros((3,3),dtype=float)
        R[0,0] = ux*ux*C + c
        R[0,1] = ux*uy*C - uz*s
        R[0,2] = ux*uz*C + uy*s
        R[1,0] = uy*ux*C + uz*s
        R[1,1] = uy*uy*C + c
        R[1,2] = uy*uz*C - ux*s
        R[2,0] = uz*ux*C - uy*s
        R[2,1] = uz*uy*C + ux*s
        R[2,2] = uz*uz*C + c
        return R
    
# create simple strip rectangular mesh
def getRecSheetMesh(L = 4, H = 1, maxArea = 1e-2):
    m = mesh.Mesh(*triangulation.triangulate(*stripBoundary(L,H), triArea=maxArea)[0:2], embeddingDimension=3)
    return m

def getRecTetMesh(thickness, L = 4, H = 1, maxVol = 2e-2, tetdegree = 2):
    pts, _ = stripBoundary(L,H)
    numP = len(pts)
    pts3d_array = np.zeros((numP,3),dtype=float)
    pts3d_array[:,0:2] = pts
    # outer shell z + thickness/2
    pts_outer = copy.deepcopy(pts3d_array)
    pts_outer[:,2] += thickness/2
    #inner shell z - thickness/2
    pts_inner = copy.deepcopy(pts3d_array)
    pts_inner[:,2] -= thickness/2
    
    V = np.vstack((pts_outer,pts_inner))
    F = [[0,3,2,1],
         [0,1,5,4],
         [2,3,7,6],
         [3,0,4,7],
         [1,2,6,5],
         [5,6,7,4]]
    
    #build tet mesh
    tm_info = MeshInfo()
    tm_info.set_points(V)
    tm_info.set_facets(F)
    tmesh = build(tm_info,max_volume=maxVol)
    tm = mesh.Mesh(tmesh.points,tmesh.elements,degree=tetdegree)   
    return tm

# create triangular mesh of creased rectangle
def getCreasedRecSheetMesh(alpha, L = 4, H = 1, maxArea = 1e-2):
    rotangle = np.pi/2 - alpha/2
    # create flat rectangle mesh 
    V, F, edgeMarkers = triangulation.triangulate(*creasedBoundary(L,H), triArea=maxArea, outputPointMarkers=False, outputEdgeMarkers=True)
    m = mesh.Mesh(V, F,embeddingDimension=3)
    # get creases list
    isBoundary = np.zeros(m.numVertices(), dtype=bool)
    isBoundary[m.boundaryVertices()] = True
    creases = np.array([em for em in edgeMarkers if not isBoundary[em].all()])
    
    # copy of original vertices
    vertex_ori = copy.deepcopy(V)
    vertex_rotate = np.empty_like(vertex_ori)
    vertex_rotate[:] = vertex_ori
    # left vertices
    vertex_left_ind = np.where(vertex_ori[:,0]<(L/2))
    vertex_left = vertex_ori[vertex_left_ind]
    # left vertex translate
    vertex_left[:,0] -= L/2
    # left vertices rotate
    RotMatL = RotateMat3D(0,1,0,-rotangle)
    vertex_left_rotate_T = np.matmul(RotMatL,np.transpose(vertex_left))
    vertex_left_rotate = np.transpose(vertex_left_rotate_T)
    vertex_left_rotate[:,0] += L/2 # translate back
    vertex_rotate[vertex_left_ind] = vertex_left_rotate
    # right vertices
    vertex_right_ind = np.where(vertex_ori[:,0]>(L/2))
    vertex_right = vertex_ori[vertex_right_ind]
    # right vertex translate
    vertex_right[:,0] -= L/2
    # right vertex rotate
    RotMatR = RotateMat3D(0,1,0,rotangle)
    vertex_right_rotate_T = np.matmul(RotMatR,np.transpose(vertex_right))
    vertex_right_rotate = np.transpose(vertex_right_rotate_T)
    vertex_right_rotate[:,0] += L/2 # translate back
    vertex_rotate[vertex_right_ind] = vertex_right_rotate
    
    m = mesh.Mesh(vertex_rotate,F)   
    return m, creases
    
def getCreasedRecTetMesh(alpha, thickness, L=4, H=1, maxVol = 2e-2, tetdegree = 2):   
    rotangle = np.pi/2 - alpha/2
    h = thickness/np.tan(alpha/2)
    
    pts, edges = creasedBoundary(L,H)
    numP = len(pts)    
    
    # for now only consider alpha > 0 and <= pi
    if alpha <= 0 or alpha > np.pi:
        print("Error: alpha shoule be in the range (0,pi] !")
        return -1
    elif alpha == np.pi:
        pts3d_array = np.zeros((numP,3),dtype=float)
        pts3d_array[:,0:2] = pts
        # outer shell z + thickness/2
        pts_outer = copy.deepcopy(pts3d_array)
        pts_outer[:,2] += thickness/2
        #inner shell z - thickness/2
        pts_inner = copy.deepcopy(pts3d_array)
        pts_inner[:,2] -= thickness/2
        
    else:
        h = thickness/np.tan(alpha/2)     
        pts3d_array = np.zeros((numP,3),dtype=float)
        pts3d_array[:,0:2] = pts
        edges_array = np.array(edges,dtype=int)
        pts3d_rotate = copy.deepcopy(pts3d_array)

        # left two corner points
        pts_left_ind = np.where(pts3d_array[:,0]<(L/2))
        pts_left = pts3d_array[pts_left_ind]
        pts_left[:,0] -= L/2
        RotMatL = RotateMat3D(0,1,0,-rotangle)
        pts_left_rotate = np.transpose(np.matmul(RotMatL,np.transpose(pts_left)))
        pts_left_rotate[:,0] += L/2
        pts3d_rotate[pts_left_ind] = pts_left_rotate
        # right two corner points
        pts_right_ind = np.where(pts3d_array[:,0]>(L/2))
        pts_right = pts3d_array[pts_right_ind]
        pts_right[:,0] -= L/2
        RotMatR = RotateMat3D(0,1,0,rotangle)
        pts_right_rotate = np.transpose(np.matmul(RotMatR,np.transpose(pts_right)))
        pts_right_rotate[:,0] += L/2
        pts3d_rotate[pts_right_ind] = pts_right_rotate

        pts_middle_ind = np.where(pts3d_array[:,0] == (L/2))
        #outershell
        pts_outer = copy.deepcopy(pts3d_rotate)
        pts_temp_ps = pts_outer[pts_middle_ind]
        pts_temp_ps[:,2] += h/2
        pts_outer[pts_middle_ind] = pts_temp_ps
        pts_temp_ps = pts_outer[pts_left_ind]
        pts_temp_ps[:,0] -= thickness/2
        pts_outer[pts_left_ind] = pts_temp_ps
        pts_temp_ps = pts_outer[pts_right_ind]
        pts_temp_ps[:,0] += thickness/2
        pts_outer[pts_right_ind] = pts_temp_ps
        #inner shell
        pts_inner = copy.deepcopy(pts3d_rotate)
        pts_temp_ps = pts_inner[pts_middle_ind]
        pts_temp_ps[:,2] -= h/2
        pts_inner[pts_middle_ind] = pts_temp_ps
        pts_temp_ps = pts_inner[pts_left_ind]
        pts_temp_ps[:,0] += thickness/2
        pts_inner[pts_left_ind] = pts_temp_ps
        pts_temp_ps = pts_inner[pts_right_ind]
        pts_temp_ps[:,0] -= thickness/2
        pts_inner[pts_right_ind] = pts_temp_ps
    
    V = np.vstack((pts_outer,pts_inner))
    F = [[0,5,2,1],
          [0,1,7,6],
          [0,6,11,5],
          [11,6,7,8],
          [1,2,8,7],
          [5,4,3,2],
          [5,11,10,4],
          [10,9,3,4],
          [10,11,8,9],
          [9,8,2,3]]
    
    #build tet mesh
    tm_info = MeshInfo()
    tm_info.set_points(V)
    tm_info.set_facets(F)
    tmesh = build(tm_info,max_volume=maxVol)
    tm = mesh.Mesh(tmesh.points,tmesh.elements,degree=tetdegree)   
    return tm

# get hemisphere sheet mesh and creases
def getCreasedHemisphereSheetMesh(thickness, resolution = 5):
    # open obj mesh
    thisfilepath = os.path.dirname(os.path.realpath(__file__))
    filepath_parent = os.path.abspath(os.path.join(thisfilepath, os.pardir))
    basepath = os.path.join(filepath_parent,'demos','Data','creased_shell_')
    filepath = basepath + str(resolution) + '.obj'
    m = mesh.Mesh(filepath)
    
    # determine creases (exclude boundary edges)
    edges = []
    edgeFaceIdx = []
    m.visitEdges(lambda e, i: edges.append(e))
    m.visitEdgesFaceIdxPair(lambda e,i: edgeFaceIdx.append(e))
    boundaryEdgeIdx = []
    for i in range(len(edgeFaceIdx)):
        if edgeFaceIdx[i][0] == -1 or edgeFaceIdx[i][1] == -1 :
            boundaryEdgeIdx.append(i)
    
    # print(boundaryEdgeIdx)
    # delete boundaryEdge in cedges[] list
    cedges = copy.deepcopy(edges)        
    for ele in sorted(boundaryEdgeIdx,reverse=True):
        del cedges[ele]
    # print(cedges)
    creases = np.array(cedges,dtype=np.int32)
    return m, creases

def getCreasedHemisphereTetMesh(thickness, maxVol = 2e-3, tetdegree = 2):
    # open obj mesh
    thisfilepath = os.path.dirname(os.path.realpath(__file__))
    filepath_parent = os.path.abspath(os.path.join(thisfilepath, os.pardir))
    meshpath = os.path.join(filepath_parent,'demos','Data','creased_shell_1.obj')
    m_ori = mesh.Mesh(meshpath)
    
    # deep copy of original mesh
    m = mesh.Mesh(m_ori.vertices(),m_ori.elements())
    
    # generate outer shell and inner shell according to thinkness/2
    m_outer = mesh.Mesh(m.vertices() * (1+(thickness / 2)), m.elements())
    m.setVertices(m.vertices() * (1-(thickness / 2)))
    # get boundary node sequences
    m_bdryloops = m.boundaryNodes()[m.boundaryLoops()]
    m_bdryloops_c = m_bdryloops[0]
    
    # concatenateMeshes m and m_outer
    m_tris = m.elements()
    # change orientation of the outer shell
    m_tris[:,[1,0]] = m_tris[:,[0,1]]
    V_combined, F_combined = mesh_operations.concatenateMeshes([(m.vertices(),m.elements()), (m_outer.vertices(), m_tris)])
    # generate tri mesh of bottom surface: we have to stitch the outer shell and the inner shell
    m_boundnode_idx = m_bdryloops_c
    V_offset_val = m.vertices().shape[0]

    stitchTris = []
    for i in range(m_boundnode_idx.shape[0]):
        vind_1 = m_boundnode_idx[i]
        next_i = i + 1
        if next_i == m_boundnode_idx.shape[0]:
            next_i = 0
        vind_2 = m_boundnode_idx[next_i]
        stitchTris.append([vind_1, vind_2, vind_1 + V_offset_val])
        stitchTris.append([vind_2, vind_2 + V_offset_val, vind_1 + V_offset_val])

    stTris = np.array(stitchTris)
    F_combined_stitch = np.vstack((F_combined, stTris))
    m_concatenate = mesh.Mesh(V_combined, F_combined_stitch)
    
    # Generate Tet Mesh
    tetmesh_info = MeshInfo()
    tetmesh_info.set_points(m_concatenate.vertices())
    tetmesh_info.set_facets(m_concatenate.elements())
    tmesh = build(tetmesh_info, max_volume=maxVol)
    V_tmesh = np.array(tmesh.points)
    F_tmesh = np.array(tmesh.elements)
    tm = mesh.Mesh(V_tmesh,F_tmesh, degree=tetdegree)
    return tm
    
    
def getElasticSheet(smesh,thickness,creases = [], useNeoHookean = True,useCreases = True):
    if smesh.is_tet_mesh():
        print('ERROR: INPUT Tet Mesh!')
        return -1
    if useNeoHookean: psi = energy.NeoHookeanYoungPoisson (2, 2000, 0.3)
    else:             psi = energy.StVenantKirchhoffCBased(tensors.ElasticityTensor2D(2000, 0.3))
    
    if useCreases: plate = elastic_sheet.ElasticSheet(smesh,psi,creases)
    else: plate = elastic_sheet.ElasticSheet(smesh,psi)
    plate.thickness = thickness
    return plate

def getElasticSolid(mesh,useNeoHookean = True):
    if not mesh.is_tet_mesh():
        print('ERROR: Input Mesh is not Tet Mesh!')
        return -1
    if useNeoHookean: psi = energy.NeoHookeanYoungPoisson(3,2000,0.3)
    else: psi = energy.IsotropicLinearElastic(3, 2000, 0.3)
    esolid = elastic_solid.ElasticSolid(mesh, psi)
    return esolid
    
def gravitySimulation(obj, gravity = -9.80635e-3, opts=None):
    
    # Optimizer
    if opts is None:
        opts = py_newton_optimizer.NewtonOptimizerOptions()
        opts.niter = 100
        opts.gradTol = 1e-12
        
    bottomThetas = getBBoxVars(obj, BBoxFace.MIN_Z, displacementComponents=[])
    bottomEdgeVars = getBBoxVars(obj, BBoxFace.MIN_Z)
    
    start_t = time.time()
    
    # sheet simulation
    if hasattr(obj, 'numCreases'):
        obj.programFlatRestCurvature()
        obj.computeEquilibrium(loads = [], fixedVars=np.arange(obj.creaseAngleOffset()), opts=opts)
        obj.programRestCurvature()      
        g = loads.Gravity(obj, 1.06e-3*obj.thickness, [0, 0, gravity])
        # Compute crease angles between each two faces
        creaseVars = np.arange(obj.numCreases()) + obj.creaseAngleOffset()
        # constrain all crease angles
        obj.computeEquilibrium(loads = [g], fixedVars=bottomEdgeVars+bottomThetas+list(creaseVars), opts=opts)
    
    # solid simulation
    else:
        g = loads.Gravity(obj, 1.06e-3, [0, 0, gravity])
        obj.computeEquilibrium(loads = [g], fixedVars=bottomEdgeVars+bottomThetas, opts=opts)
    
    return obj, time.time()-start_t

def cantilverGraSimulation(obj, gravity = -9.80635e-3, opts = None):
    
    # Optimizer
    if opts is None:
        opts = py_newton_optimizer.NewtonOptimizerOptions()
        opts.niter = 100
        opts.gradTol = 1e-12
    start_t = time.time()
    
    leftThetas = getBBoxVars(obj, BBoxFace.MIN_X, displacementComponents=[])
    leftEdgeVars = getBBoxVars(obj, BBoxFace.MIN_X, displacementComponents=[0,2])
    leftEdgeVars.append(leftEdgeVars[0] + 1)
    
    # sheet simulation
    if hasattr(obj, 'numCreases'):      
        g = loads.Gravity(obj, 1.06e-3*obj.thickness, [0, 0, gravity])
        obj.computeEquilibrium(loads = [g], fixedVars=leftThetas+leftEdgeVars, opts=opts)
    
    # solid simulation
    else:
        g = loads.Gravity(obj, 1.06e-3, [0, 0, gravity])
        obj.computeEquilibrium(loads = [g], fixedVars=leftThetas+leftEdgeVars, opts=opts)
    
    return obj, time.time()-start_t

# Compute the vertex-averaged maximum principal strains
import field_sampler
maxEigenvalue = lambda e: np.linalg.eigh(e)[0].max()
minEigenvalue = lambda e: np.linalg.eigh(e)[0].min()
frobeniusNorm = lambda e: np.linalg.norm(e)
fullTensor    = lambda e: e

def sampleStrains(esolid, esheet, h):
    """
    Sample the volumetric strain fields computed by the sheet and solid simulations
    at the vertices of the sheet mesh, and at offset `h` away from the midsurface.

    Returns a pair holding the sheet's strain field, followed by the strains sampled from the solid.
    """
    if not hasattr(esheet, 'numCreases'):
        print('ERROR! Wrong Simulation Object for esheet!')
        return -1

    sheetStrains = esheet.getVertexVolumetricStrains(h)
    fs = field_sampler.FieldSampler(esolid.mesh())

    samplePts = esheet.mesh().vertices() + h * esheet.mesh().vertexNormals()
    sampledSolidStrains = [esolid.greenStrain(ei, bc) for ei, bc in zip(*fs.closestElementAndBaryCoords(samplePts))]
    return sheetStrains, sampledSolidStrains

def sampleStresses(esolid, esheet, h):
    """
    Sample the volumetric stress fields computed by the sheet and solid simulations
    at the vertices of the sheet mesh, and at offset `h` away from the midsurface.

    Returns a pair holding the sheet's stress field, followed by the stresses sampled from the solid.
    """
    if not hasattr(esheet, 'numCreases'):
        print('ERROR! Wrong Simulation Object for esheet!')
        return -1

    sheetStresses = esheet.getVertexCauchyStresses(h)
    fs = field_sampler.FieldSampler(esolid.mesh())

    samplePts = esheet.mesh().vertices() + h * esheet.mesh().vertexNormals()
    sampledSolidStresses = [esolid.cauchyStress(ei, bc) for ei, bc in zip(*fs.closestElementAndBaryCoords(samplePts))]
    return sheetStresses, sampledSolidStresses

def computeVertexStrain(esolid_sim_obj, esheet_sim_obj, h = None, scalarMeasure=maxEigenvalue):
    thickness = esheet_sim_obj.thickness
    if h is None: h = thickness / 2

    sheetStrains, sampledSolidStrains = sampleStrains(esolid_sim_obj, esheet_sim_obj, h)

    strainMeasure_sheet = np.array([scalarMeasure(e) for e in sheetStrains])
    strainMeasure_solid = np.array([scalarMeasure(e) for e in esolid_sim_obj.vertexGreenStrains()])
    strainMeasure_solid_on_sheet = np.array([scalarMeasure(e) for e in sampledSolidStrains])

    rel_error = np.linalg.norm(strainMeasure_sheet - strainMeasure_solid_on_sheet) / np.linalg.norm(strainMeasure_solid_on_sheet)

    return rel_error, strainMeasure_sheet, strainMeasure_solid, strainMeasure_solid_on_sheet, sheetStrains, sampledSolidStrains

def computeVertexStresses(esolid_sim_obj, esheet_sim_obj, h = None, scalarMeasure=maxEigenvalue):
    thickness = esheet_sim_obj.thickness
    if h is None: h = thickness / 2

    sheetStresses, sampledSolidStresses = sampleStresses(esolid_sim_obj, esheet_sim_obj, h)

    stressMeasure_sheet = np.array([scalarMeasure(e) for e in sheetStresses])
    stressMeasure_solid = np.array([scalarMeasure(e) for e in esolid_sim_obj.vertexCauchyStresses()])
    stressMeasure_solid_on_sheet = np.array([scalarMeasure(e) for e in sampledSolidStresses])

    rel_error = np.linalg.norm(stressMeasure_sheet - stressMeasure_solid_on_sheet) / np.linalg.norm(stressMeasure_solid_on_sheet)

    return rel_error, stressMeasure_sheet, stressMeasure_solid, stressMeasure_solid_on_sheet, sheetStresses, sampledSolidStresses

# visualize strain field
import tri_mesh_viewer
import vis
def visVertexStrainOnSheetAndSolid(esolid_sim_obj, esheet_sim_obj, h = None, scalarMeasure=maxEigenvalue):
    _,emax_sheet, emax_solid, emax_solid_on_sheet, _, _  = computeVertexStrain(esolid_sim_obj, esheet_sim_obj, h, scalarMeasure)
    vmin = min(np.min(emax_sheet), np.min(emax_solid_on_sheet))
    vmax = max(np.max(emax_sheet), np.max(emax_solid_on_sheet))
    vmin = np.min(emax_solid)
    vmax = np.max(emax_solid) 
    vdefo = tri_mesh_viewer.Viewer(esolid_sim_obj, scalarField={'data': emax_solid, 'vmin': vmin, 'vmax': vmax})
    sheetView = TriMeshViewer(esheet_sim_obj, scalarField={'data': emax_sheet, 'vmin': vmin, 'vmax': vmax})
    sampledSolidView = TriMeshViewer(esheet_sim_obj, scalarField={'data': emax_solid_on_sheet, 'vmin': vmin, 'vmax': vmax})
    # sampledSolidView.update(scalarField={'data': np.abs(emax_solid_on_sheet - emax_sheet), 'vmin': vmin, 'vmax': vmax})
    return sheetView, sampledSolidView, vdefo

def visVertexStressOnSheetAndSolid(esolid_sim_obj, esheet_sim_obj, h = None, scalarMeasure=maxEigenvalue):
    _,emax_sheet, emax_solid, emax_solid_on_sheet, _, _  = computeVertexStresses(esolid_sim_obj, esheet_sim_obj, h, scalarMeasure)
    vmin = min(np.min(emax_sheet), np.min(emax_solid_on_sheet))
    vmax = max(np.max(emax_sheet), np.max(emax_solid_on_sheet))
    vmin = np.min(emax_solid)
    vmax = np.max(emax_solid) 
    vdefo = tri_mesh_viewer.Viewer(esolid_sim_obj, scalarField={'data': emax_solid, 'vmin': vmin, 'vmax': vmax})
    sheetView = TriMeshViewer(esheet_sim_obj, scalarField={'data': emax_sheet, 'vmin': vmin, 'vmax': vmax})
    sampledSolidView = TriMeshViewer(esheet_sim_obj, scalarField={'data': emax_solid_on_sheet, 'vmin': vmin, 'vmax': vmax})
    # sampledSolidView.update(scalarField={'data': np.abs(emax_solid_on_sheet - emax_sheet), 'vmin': vmin, 'vmax': vmax})
    
    return sheetView, sampledSolidView, vdefo

    
# Created By Xinzhuo (Johnson) Hu 2023-01-10 04:31
def thickConvergenceSweep(solidMeshFunc, shellMeshFunc, simFunc,thicks = np.linspace(0.01,0.1,10)):
    
    result = { 'thickness': [],
               'solidEnergies': [],
               'shellEnergies': [],
               'solidSimTimes': [],
               'shellSimTimes': [],
               'relativeErrors': []}
    
    thicksr = np.flip(thicks)
    for i, thick in enumerate(thicksr):
        print(f'Sim Iter {i+1}/{len(thicksr)}',end = '\r', flush=True)
        result['thickness'].append(thick)
        
        tm = solidMeshFunc(thick)
        esolid, t1 = simFunc(getElasticSolid(tm)) 
        esolid_energy = esolid.energy()
        result['solidEnergies'].append(esolid_energy)
        result['solidSimTimes'].append(t1)
        
        sm = shellMeshFunc()  
        eshell,t2 = simFunc(getElasticSheet(sm,thick,useCreases=False))
        eshell_energy = eshell.energy()
        result['shellEnergies'].append(eshell_energy)
        result['shellSimTimes'].append(t2)
        
        rel_error = (np.abs(esolid_energy-eshell_energy))/esolid_energy
        result['relativeErrors'].append(rel_error)
        
    return result
    


    
# !!! The following functions need to be reorganized !!!



# function to simulate semi-sphere tet mesh under user-defined gravity
def semisphereTetGraSimulation(thickness, resolution = 1, maxVol = 0.02, tetdegree = 2, gravity = -9.80635e-3):
    # open obj mesh
    basepath = '/Users/xinzhuohu/PersonalFiles/Projects/MeshFEM_dev/python/demos/Data/creased_shell_'
    filepath = basepath + str(resolution) + '.obj'
    m_ori = mesh.Mesh(filepath)
    
    # deep copy of original mesh
    m = mesh.Mesh(m_ori.vertices(),m_ori.elements())
    
    # generate outer shell and inner shell according to thinkness/2
    m_outer = mesh.Mesh(m.vertices() * (1+(thickness / 2)), m.elements())
    m.setVertices(m.vertices() * (1-(thickness / 2)))
    # get boundary node sequences
    m_bdryloops = m.boundaryNodes()[m.boundaryLoops()]
    m_bdryloops_c = m_bdryloops[0]
    
    # concatenateMeshes m and m_outer
    m_tris = m.elements()
    # change orientation of the outer shell
    m_tris[:,[1,0]] = m_tris[:,[0,1]]
    V_combined, F_combined = mesh_operations.concatenateMeshes([(m.vertices(),m.elements()), (m_outer.vertices(), m_tris)])
    # generate tri mesh of bottom surface: we have to stitch the outer shell and the inner shell
    m_boundnode_idx = m_bdryloops_c
    V_offset_val = m.vertices().shape[0]

    stitchTris = []
    for i in range(m_boundnode_idx.shape[0]):
        vind_1 = m_boundnode_idx[i]
        next_i = i + 1
        if next_i == m_boundnode_idx.shape[0]:
            next_i = 0
        vind_2 = m_boundnode_idx[next_i]
        stitchTris.append([vind_1, vind_2, vind_1 + V_offset_val])
        stitchTris.append([vind_2, vind_2 + V_offset_val, vind_1 + V_offset_val])

    stTris = np.array(stitchTris)
    F_combined_stitch = np.vstack((F_combined, stTris))
    m_concatenate = mesh.Mesh(V_combined, F_combined_stitch)
    
    # Generate Tet Mesh
    tetmesh_info = MeshInfo()
    tetmesh_info.set_points(m_concatenate.vertices())
    tetmesh_info.set_facets(m_concatenate.elements())
    tmesh = build(tetmesh_info, max_volume=maxVol)
    V_tmesh = np.array(tmesh.points)
    F_tmesh = np.array(tmesh.elements)
    tm = mesh.Mesh(V_tmesh,F_tmesh, degree=tetdegree)
    
    # Generate Elastic Object
    # psi = energy.NeoHookeanYoungPoisson (3, 0.2587, 1.5828e-3) # material properity
    psi = energy.NeoHookeanYoungPoisson (3, 0.07, 0.4)
    semiss = elastic_solid.ElasticSolid(tm,psi)
    
    # Dirichlet BC: constrain bottom Z nodes and normals
    bottomThetas = getBBoxVars(semiss, BBoxFace.MIN_Z, displacementComponents=[])
    bottomEdgeVars = getBBoxVars(semiss, BBoxFace.MIN_Z)
    # Neuman BC: gravity load
    # Gravity(obj, rho, g-vector)
    g = loads.Gravity(semiss, 1.06e-3, [0, 0, gravity])
    
    # Simulation
    opts = py_newton_optimizer.NewtonOptimizerOptions()
    opts.niter = 100
    opts.gradTol = 1e-10
    start_t = time.time()
    semiss.computeEquilibrium(loads = [g], fixedVars=bottomEdgeVars+bottomThetas, opts=opts)
    return semiss, time.time()-start_t

# function to simulate semi-sphere shell mesh under user-defined gravity
def semisphereCreaseGraSimulation(thickness, resolution = 1, gravity = -9.80635e-3):
    # open obj mesh
    basepath = '/Users/xinzhuohu/PersonalFiles/Projects/MeshFEM_dev/python/demos/Data/creased_shell_'
    filepath = basepath + str(resolution) + '.obj'
    m = mesh.Mesh(filepath)
    
    # determine creases (exclude boundary edges)
    edges = []
    edgeFaceIdx = []
    m.visitEdges(lambda e, i: edges.append(e))
    m.visitEdgesFaceIdxPair(lambda e,i: edgeFaceIdx.append(e))
    boundaryEdgeIdx = []
    for i in range(len(edgeFaceIdx)):
        if edgeFaceIdx[i][0] == -1 or edgeFaceIdx[i][1] == -1 :
            boundaryEdgeIdx.append(i)
    
    # print(boundaryEdgeIdx)
    # delete boundaryEdge in cedges[] list
    cedges = copy.deepcopy(edges)        
    for ele in sorted(boundaryEdgeIdx,reverse=True):
        del cedges[ele]
    # print(cedges)
    creases = np.array(cedges,dtype=np.int32)
    
    # Generate elastic shell object
    # rubber material
    # psi = energy.NeoHookeanYoungPoisson (2, 0.2587, 1.5828e-3)
    psi = energy.NeoHookeanYoungPoisson (2, 0.07, 0.4)
    # constructor with creases
    es = elastic_sheet.ElasticSheet(m,psi,creases)
    es.thickness = thickness
    
    # find optimal crease angle by optimizing from flat configuration
    opts = py_newton_optimizer.NewtonOptimizerOptions()
    opts.niter = 100
    opts.gradTol = 1e-10
    start_t = time.time()
    es.programFlatRestCurvature()
    es.computeEquilibrium(loads = [], fixedVars=np.arange(es.creaseAngleOffset()), opts=opts)
    print("Shell Energy after Finding Optimal Crease Angle:  Total Energy | Membrane Energy | Bending Energy")
    print(es.energy(), es.energy(etype=es.EnergyType.Membrane), es.energy(etype=es.EnergyType.Bending))
    es.programRestCurvature()
    # print("Shell Energy after Finding Optimal Crease Angle:  Total Energy | Membrane Energy | Bending Energy")
    # print(es.energy(), es.energy(etype=es.EnergyType.Membrane), es.energy(etype=es.EnergyType.Bending))
    
    # Dirichlet BC: constrain bottom Z nodes and normals
    bottomThetas = getBBoxVars(es, BBoxFace.MIN_Z, displacementComponents=[])
    bottomEdgeVars = getBBoxVars(es, BBoxFace.MIN_Z)
    # Neuman BC: gravity
    # density should multiply thickness
    g = loads.Gravity(es, 1.06e-3*thickness, [0, 0, gravity])
    # Compute crease angles between each two faces
    creaseVars = np.arange(es.numCreases()) + es.creaseAngleOffset()
    # constrain all crease angles
    es.computeEquilibrium(loads = [g], fixedVars=bottomEdgeVars+bottomThetas+list(creaseVars), opts=opts)
    
    print("Shell Energy after Gravity Simulation:  Total Energy | Membrane Energy | Bending Energy")
    print(es.energy(), es.energy(etype=es.EnergyType.Membrane), es.energy(etype=es.EnergyType.Bending))
    
    return es, time.time()-start_t
    
def tetGraSimConvergenceSweep(thickness, maxVols = np.logspace(-2, -9, 30,base = 5.0), includeDeg1 = False):
    energies = {1 : [], 2 : []}
    times = {1 :[], 2 : []}
    elements = {1 : [], 2 : []}
    edgeLens = {1 : [], 2 : []}
    minedgeLen = {1 : [], 2 : []}
    for i, maxVolScale in enumerate(maxVols):
        for deg in [1, 2] if includeDeg1 else [2]:
            print(f'Tet Sum {i + 1}/{len(maxVols)} deg {deg}', end='\r', flush = True)
            maxVol = maxVolScale * (1 if deg == 2 else 0.1)
            semiss, t = semisphereTetGraSimulation(thickness, resolution=1, maxVol=maxVol, tetdegree=deg)
            energies[deg].append(semiss.energy())
            times[deg].append(t)
            elements[deg].append(semiss.mesh().numElements())
            edgeLens[deg].append(np.median(semiss.mesh().edgeLengths()))
            minedgeLen[deg].append(semiss.mesh().edgeLengths().min())
    return {'energies': energies,
            'times': times,
            'elements': elements,
            'edgeLens': edgeLens,
            'minedgeLen':minedgeLen}

def shellGraSimConvergenceSweep(thickness, ResList):
    result = { 'times':            [],
               'energies':         [],
               'bendingEnergies':  [],
               'membraneEnergies': [],
               'edgeLens':         [],
               'elements':         [],
               'minedgeLen':       []}
    for i, res in enumerate(ResList):
        print(f'Sheet sim {i + 1}/{len(ResList)}', end='\r', flush=True)
        eshell,t = semisphereCreaseGraSimulation(thickness, res)
        result['times'].append(t)
        result['energies'].append(eshell.energy())
        result['bendingEnergies'].append(eshell.energy(etype = eshell.EnergyType.Bending))
        result['membraneEnergies'].append(eshell.energy(etype = eshell.EnergyType.Membrane))
        result['elements'].append(eshell.mesh().numElements())
        result['edgeLens'].append(np.median(eshell.mesh().edgeLengths()))
        result['minedgeLen'].append(eshell.mesh().edgeLengths().min())
    return result

# Created By Xinzhuo (Johnson) Hu 2022-12-15 06:00 pm   
def thicknessConvergenceTest_OLD(solidSimFunc, shellSimFunc,thicks = np.linspace(0.01,0.1,10), maxV = 2e-5, Res = 7):
    
    result = { 'thickness': [],
               'solidEnergies': [],
               'shellEnergies': [],
               'relativeErrors': []}
    
    thicksr = np.flip(thicks)
    for i, thick in enumerate(thicksr):
        print(f'Sim Iter {i+1}/{len(thicksr)}',end = '\r', flush=True)
        result['thickness'].append(thick)
        esolid,_ = solidSimFunc(thick, maxVol = maxV) #default degree 2
        esolid_energy = esolid.energy()
        result['solidEnergies'].append(esolid_energy)
        eshell,_ = shellSimFunc(thick, Res)
        eshell_energy = eshell.energy()
        result['shellEnergies'].append(eshell_energy)
        rel_error = (np.abs(esolid_energy-eshell_energy))/esolid_energy
        result['relativeErrors'].append(rel_error)
        
    return result

from matplotlib import pyplot as plt
def convergencePlot(thickness, shellResult, tetResult, includeDeg1 = False):
    plt.figure(figsize=(14, 8))
    for i, xaxis in enumerate(['times', 'edgeLens']):
        plt.subplot(2, 1, i+1)
        plt.semilogx(shellResult[xaxis], shellResult['energies'], label = 'semisphere-shell')
        plt.semilogx(tetResult[xaxis][2], tetResult['energies'][2], label = 'tet (degree 2)')
        if includeDeg1:
            plt.semilogx(tetResult[xaxis][1], tetResult['energies'][1], label = 'tet (degree 1)')
        plt.xlabel(xaxis)
        plt.legend()
        plt.title(f'Elastic Shell Gravity Simulation Convergence -- Thickness {thickness}')
        plt.ylabel('Elastic Energy')
        plt.grid()
    plt.tight_layout()

def thickConvergencePlot(result):
    plt.figure(figsize=(14,8))
    plt.subplot(2,1,1)
    plt.plot(result['thickness'],result['solidEnergies'],label = 'Solid Energy')
    plt.plot(result['thickness'],result['shellEnergies'],label = 'Shell Energy')
    plt.xlabel('Thickness')
    plt.legend()
    plt.ylabel('Elastic Energy')
    plt.grid()
    plt.title('Solid and Shell Simulation Energy')
    
    plt.subplot(2,1,2)
    plt.plot(result['thickness'],result['relativeErrors'],label = 'Relative Error')
    plt.xlabel('Thickness')
    plt.ylabel('Relative Errors')
    plt.grid()
    plt.title('Energy Comparsion Between Solid and Shell Simulation by Different Thickness')
    plt.tight_layout()


# ### Save pattern
isheet = inflation.InflatableSheet(m, fusing_data)
uv = np.load(output_data_path + '/rparam_uv.npy')

paramSampler = SurfaceSampler(np.pad(uv, [(0, 0), (0, 1)], 'constant'), target_surf.triangles())
liftedSheetPositions = paramSampler.sample(m.vertices(), target_surf.vertices())
isheet.setUninflatedDeformation(liftedSheetPositions.transpose())

opts = py_newton_optimizer.NewtonOptimizerOptions()
opts.useIdentityMetric = True
opts.beta = 1e-4
opts.gradTol = 1e-10

sys.path.append("../../")

# Are the flat region causing a problem? They might not actually control the metric...
# Try replacing them with single wall...
# Analyze the actual stretching factor (much easier to do with skeleton walls)

bdryVars = boundaries.getOuterBoundaryVars(isheet)
isheet.setUseTensionFieldEnergy(True)
isheet.setUseHessianProjectedEnergy(False)
# First generate results with fixed boundary 
viewer = TriMeshViewerWithSurface(isheet, target_surf, width=768, height=640)
viewer.showWireframe(True)

# viewer.setCameraParams(((1.613494603240345, -3.9332708615926393, 1.4922998234349831),
# (-0.05948468564942635, 0.33267929672385665, 0.941162078339598),
# (0.0, 0.0, 0.0)))


framerate = 10
def cb(it):
    return
    if it % framerate == 0:
        viewer.update(scalarField=utils.getStrains(isheet)[:, 0])    


# ### First solve with low pressure to get out of indefinite state
isheet.pressure = 1e-7
opts.niter = 20
cr = inflation.inflation_newton(isheet, bdryVars, opts, hessianShift = 0, callback = cb)

isheet.pressure = 1e-3
opts.niter = 20
cr = inflation.inflation_newton(isheet, bdryVars, opts, hessianShift = 0, callback = cb)

isheet.pressure = 1e-2
opts.niter = 20
cr = inflation.inflation_newton(isheet, bdryVars, opts, hessianShift = 0, callback = cb)

fixedVars_list = [bdryVars, []]
tag_name = ['fixed_boundary', 'free_boundary']
hessian_shifts = [0, 1e-6]
viewer.update(scalarField=utils.getStrains(isheet)[:, 0])    

for i in range(len(fixedVars_list)):
    fixedVars = fixedVars_list[i]
    tag = tag_name[i]
    hessian_shift = hessian_shifts[i]
    if os.path.exists(output_data_path + '/{}_inflated_sheet_vars.npy'.format(tag)):
        isheet.setVars(np.load(output_data_path + '/{}_inflated_sheet_vars.npy'.format(tag)))
        viewer.update(scalarField=utils.getStrains(isheet)[:, 0])    

    else:
        # ### Then inflate
        isheet.pressure = 0.025
        opts.niter = 100
        opts.gradTol = 1e-7

        benchmark.reset()
        cr = inflation.inflation_newton(isheet, fixedVars, opts, hessianShift = hessian_shift, callback = cb)
        benchmark.report()

        np.save(output_data_path + '/{}_inflated_sheet_vars.npy'.format(tag), isheet.getVars())

    orender = viewer.offscreenRenderer(width=1024,height=1024)
    orender.render()
    orender.save(output_data_path + '/{}_parametrized_mesh_inflated.png'.format(tag))

    def export_top_bottom_mesh(isheet, export_path, shape_name, pattern_name):
        mesh_3d = isheet.visualizationMesh(True)
        mesh_2d = isheet.mesh()
        vx_3d = mesh_3d.vertices()
        elements_3d = mesh_3d.elements()

        new_mesh_3d = MeshFEM.Mesh(vx_3d[:mesh_2d.numVertices()], elements_3d[:mesh_2d.numElements()])

        new_mesh_3d.save(export_path + '/{}_{}_{}_mesh_3d_top.obj'.format(tag, shape_name, pattern_name))

        new_mesh_3d = MeshFEM.Mesh(vx_3d[m.numVertices():], elements_3d[mesh_2d.numElements():] - mesh_2d.numVertices())
        new_mesh_3d.save(export_path + '/{}_{}_{}_mesh_3d_bottom.obj'.format(tag, shape_name, pattern_name))

    export_top_bottom_mesh(isheet, output_data_path, shape_name, pattern_name)

inflation_time = time.time()
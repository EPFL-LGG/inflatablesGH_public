import sys; sys.path.append('../python')
import re, glob, utils, MeshFEM, mesh, os, subprocess

mshDir, tgtMesh = None, None

if len(sys.argv) == 3:
    mshDir, tgtMesh = sys.argv[1:]
else:
    mshDir,         = sys.argv[1:]

print(f"mshDir: {mshDir}")
print(f"tgtMesh: {tgtMesh}")

def natural_sort(l):
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ]
    return sorted(l, key = alphanum_key)

meshes = natural_sort(glob.glob(mshDir + '/step_*.msh'))
if len(meshes) == 0:
    meshes = natural_sort(glob.glob(mshDir + '/frame_*.msh'))
print(meshes)

fullyInflated = mesh.Mesh(meshes[-1])
xf = utils.renderingNormalization(fullyInflated.vertices(), placeAtopFloor = True)

blender = '/Applications/blender_2.79.app/Contents/MacOS/blender'

mesh.save('xf.obj', xf(fullyInflated.vertices()), fullyInflated.triangles())

if tgtMesh is not None:
    tgt = mesh.Mesh(tgtMesh)
    mesh.save('xf_tgt.obj', xf(tgt.vertices()), tgt.triangles())
    subprocess.call([blender, '../rendering/inflatable_render_blank_zoomout.blend', '--python', 'inflatable_render_scene_setup_only.py', '--', 'xf_tgt.obj', 'xf.obj'])
else:
    subprocess.call([blender, '../rendering/inflatable_render_blank_zoomout.blend', '--python', 'inflatable_render_scene_setup_only.py', '--', 'xf.obj'])

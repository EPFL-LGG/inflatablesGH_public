import sys; sys.path.append('../python')
import re, glob, utils, MeshFEM, mesh, os, subprocess
import uuid

mshDir, tgtMesh = None, None

if len(sys.argv) == 4:
    mshDir, tgtMesh, outPngPath = sys.argv[1:]
else:
    mshDir, outPngPath = sys.argv[1:]

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

# blender = '/Applications/blender_2.79.app/Contents/MacOS/blender'
blender = 'blender'

xf_mesh = f'tmp.{uuid.uuid4()}.obj'

mesh.save(xf_mesh, xf(fullyInflated.vertices()), fullyInflated.triangles())

if tgtMesh is not None:
    tgt = mesh.Mesh(tgtMesh)
    xf_tgt = f'tgt.{uuid.uuid4()}.obj'
    mesh.save(xf_tgt, xf(tgt.vertices()), tgt.triangles())
    subprocess.call([blender, '../rendering/inflatable_render_blank_zoomout.blend', '-b', '--python', 'inflatable_render_orbit_render.py', '--', xf_tgt, xf_mesh, outPngPath])
    os.unlink(xf_tgt)
else:
    subprocess.call([blender, '../rendering/inflatable_render_blank_zoomout.blend', '-b', '--python', 'inflatable_render_orbit_render.py', '--', xf_mesh, outPngPath])

os.unlink(xf_mesh)

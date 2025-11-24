import sys; sys.path.append('../python')
import re, glob, utils, MeshFEM, mesh, os, subprocess

mshDir, outDir = sys.argv[1:]

os.makedirs(outDir, exist_ok=True)

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

for i, mpath in enumerate(meshes):
    m = mesh.Mesh(mpath)
    mesh.save('tmp.obj', xf(m.vertices()), m.triangles())
    # subprocess.call([blender, '-b', '../rendering/inflatable_render_blank.blend', '--python', '../rendering/inflatable_render.py', '--', tgtMesh, eqiMesh, flipDataDir + '/equilibrium_with_target.png'])
    subprocess.call([blender, '-b', '../rendering/inflatable_render_blank.blend', '--python', '../rendering/inflatable_render.py', '--',          'tmp.obj', f'{outDir}/step_{i}.png'])

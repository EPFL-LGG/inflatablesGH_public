import bpy, sys

argv = []
if '--' in sys.argv:
    argv = sys.argv[sys.argv.index('--') + 1:]

targetSurfMeshPath = None
if (len(argv) == 1):
    sheetMeshPath, = argv
elif (len(argv) == 2):
    targetSurfMeshPath, sheetMeshPath = argv
else:
    raise Exception('usage: blender inflatable_render_blank.blend --python inflatable_render_scene_setup_only.py -- [targetSurfMeshPath] sheetMeshPath')

def setSmoothShading(blenderObject):
    mesh = blenderObject.data
    for f in mesh.polygons:
        f.use_smooth = True

def importMesh(path, materialName):
    bpy.ops.import_scene.obj(filepath=path, axis_forward='Y', axis_up='Z')
    #importedObject = bpy.context.object
    importedObject = bpy.context.selected_objects[-1]
    importedObject.data.materials.clear()
    importedObject.data.materials.append(bpy.data.materials[materialName])
    setSmoothShading(importedObject)

if (targetSurfMeshPath):
    importMesh(targetSurfMeshPath, 'TargetSurfMaterial')
importMesh(sheetMeshPath, 'SheetMaterial')

scene = bpy.context.scene
# scene.render.resolution_x = 1920
# scene.render.resolution_y = 1080
scene.render.resolution_x = 1024
scene.render.resolution_y = 576
scene.render.resolution_percentage = 100

# scene.render.threads_mode = 'FIXED'
# scene.render.threads = 36

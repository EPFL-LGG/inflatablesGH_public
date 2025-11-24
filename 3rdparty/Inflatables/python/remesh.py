import subprocess, os, MeshFEM, mesh, uuid

def remesh(m):
    env = os.environ.copy()
    env['DISPLAY'] = ':1'

    infile = f'{uuid.uuid4()}.obj'
    remeshedfile = f'{uuid.uuid4()}.obj'

    m.save(infile)
    subprocess.call(['meshlabserver', '-s', env['INFLATABLES'] + '/scripts/remesh_filter2.mlx', '-i', infile, '-o', remeshedfile], env=env)
    remeshed = mesh.Mesh(remeshedfile, embeddingDimension=3)

    os.unlink(infile)
    os.unlink(remeshedfile)

    return remeshed

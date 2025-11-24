import gmsh
import sys; sys.path.append('../')
import numpy as np

def generate_mesh_from_embeddings(filename_boundary, filename_embeddings=None, embed_pts=None, mesh_size_at_boundary=1.0, mesh_size_at_embeddings=1.0, gui=False):
    """Generate mesh from embeddings"""

    file_boundary = open(filename_boundary)

    # Initialize gmsh
    gmsh.initialize()
    gmsh.model.add("mesher_from_embeddings")

    # Read file
    dim = 2
    for line in file_boundary.readlines():
        data = [d.strip() for d in line.split(' ')]
        
        # Generate points
        if(data[0]=="v"):
            last_point_tag = gmsh.model.geo.addPoint(float(data[1]),float(data[2]),float(data[3]), mesh_size_at_boundary)

        # Generate lines
        if(data[0]=="l"):
            gmsh.model.geo.addLine(int(data[1]),int(data[2]))

    # Generate planar surface
    loop_tag = gmsh.model.geo.addCurveLoop([1,2,3,4])
    srf_tag = gmsh.model.geo.addPlaneSurface([loop_tag])

    # Synchronize
    gmsh.model.geo.synchronize()

    # Build embeddings
    if filename_embeddings is not None: 
        file_embeddings= open(filename_embeddings)
        embed_tags = []
        temp_pts = []
        embed_edges = []
        for line in file_embeddings.readlines():
            data = [d.strip() for d in line.split(' ')]
            if(data[0]=="v"):
                gmsh.model.geo.addPoint(float(data[1]),float(data[2]),float(data[3]), mesh_size_at_embeddings)
                temp_pts.append( [float(data[1]),float(data[2]),float(data[3])] )
            
            if(data[0]=="l"):
                embed_tags.append( gmsh.model.geo.addLine(last_point_tag + int(data[1]), last_point_tag + int(data[2])) )
                embed_edges.append([temp_pts[int(data[1])-1], temp_pts[int(data[2])-1]])

        # Synchronize
        gmsh.model.geo.synchronize()

        # Add embeddings to surface
        gmsh.model.mesh.embed(dim-1, embed_tags, dim, srf_tag)
    
    if embed_pts is not None:
        embed_tags = []
        for p in embed_pts:
            tag = gmsh.model.geo.addPoint(float(p[0]),float(p[1]),float(p[2]), mesh_size_at_embeddings)
            embed_tags.append(tag)

         # Synchronize
        gmsh.model.geo.synchronize()

        # Add embeddings to surface
        gmsh.model.mesh.embed(dim-2, embed_tags, dim, srf_tag)


    # Generate mesh
    if mesh_size_at_boundary<mesh_size_at_embeddings: mesh_size = mesh_size_at_boundary
    else : mesh_size = mesh_size_at_embeddings
    gmsh.option.setNumber("Mesh.MeshSizeMin", mesh_size)
    gmsh.option.setNumber("Mesh.MeshSizeMax", mesh_size)
    gmsh.model.mesh.generate(dim)

    node_tags, node_coords, node_param = gmsh.model.mesh.getNodes()
    element_tags, elements_node_tags = gmsh.model.mesh.getElementsByType(dim)

    # Get new nodes and faces
    v = node_coords.reshape((len(node_tags),3))
    f = order_faces( v, elements_node_tags.reshape((len(element_tags),3)) )
    fusing_data = [0] * len(v)
    
    if filename_embeddings is not None:
        fusing_data = generate_fusing_data_from_lines(v, embed_edges, 1e-3)
    if embed_pts is not None:
        for p in embed_pts:
            idx = find_point_closest_point(v, p)
            fusing_data[idx] = 1

    # Open gmsh GUI
    if(gui): gmsh.fltk.run()

    # End gmsh
    gmsh.finalize()
    
    return v, f, fusing_data


def generate_mesh_from_regions(filename_regions, filename_edges, mesh_size=1.0, gui=False):
    ''''Generate mesh from regions'''
    # Import mesh
    file_regions = open(filename_regions)

    # Initialize gmsh
    gmsh.initialize()
    gmsh.model.add("mesher_from_regions")

    # Tags
    region_tags = []
    dim = 2
    # Read file
    for line in file_regions.readlines():
        data = [d.strip() for d in line.split(' ')]
        
        # Generate points
        if(data[0]=="v"):
            gmsh.model.geo.addPoint(float(data[1]),float(data[2]),float(data[3]), mesh_size)

        # Generate lines
        if(data[0]=="l"):
            gmsh.model.geo.addLine(int(data[1]),int(data[2]))

        # Generate planar surfaces
        if(data[0]=="pl"):
            loop_tag = gmsh.model.geo.addCurveLoop([int(d) for d in data[1:]])
            srf_tag = gmsh.model.geo.addPlaneSurface([loop_tag])
            region_tags.append(srf_tag)
    
    # Build fusing lines
    file_edges= open(filename_edges)
    temp_pts = []
    edges = []
    for line in file_edges.readlines():
        data = [d.strip() for d in line.split(' ')]
        if(data[0]=="v"):
            temp_pts.append( [float(data[1]),float(data[2]),float(data[3])] )
            
        if(data[0]=="l"):
            edges.append([temp_pts[int(data[1])-1], temp_pts[int(data[2])-1]])

    # Generate surface loop
    gmsh.model.geo.addSurfaceLoop(region_tags)
    gmsh.model.geo.synchronize()

    # Generate mesh
    gmsh.option.setNumber("Mesh.MeshSizeMin", mesh_size)
    gmsh.option.setNumber("Mesh.MeshSizeMax", mesh_size)
    gmsh.model.mesh.generate(dim)

    node_tags, node_coords, node_param = gmsh.model.mesh.getNodes()
    element_tags, elements_node_tags = gmsh.model.mesh.getElementsByType(dim)

    # Get nodes and faces
    v = node_coords.reshape((len(node_tags),3))
    f = order_faces( v, elements_node_tags.reshape((len(element_tags),3)) )

    fusing_data = fusing_data = generate_fusing_data_from_lines(v, edges, 1e-3)

    # Open gmsh GUI
    if(gui): gmsh.fltk.run()

    # End gmsh
    gmsh.finalize()
    
    return v, f, fusing_data

def write_obj(vertices, faces, filename):
    """Write obj file"""
    file = open(filename, 'w')
    for v in vertices:
        file.write('v ' + str(v[0]) + ' ' + str(v[1]) + ' ' + str(v[2]) + '\n')
    for f in faces:
        file.write('f ' + str(f[0]) + ' ' + str(f[1]) + ' ' + str(f[2]) + '\n')
    file.close()

def order_faces(vertices, faces):
    """Order faces in a counter-clockwise fashion"""
    new_faces = []
    for f in faces:
        v1 = vertices[int(f[0]-1)]
        v2 = vertices[int(f[1]-1)]
        v3 = vertices[int(f[2]-1)]
        v12 = v2-v1
        v13 = v3-v1
        n = np.cross(v12, v13)
        if(n[2]>0):
            new_faces.append([int(f[0]), int(f[1]), int(f[2])])
        else:
            new_faces.append([int(f[0]), int(f[2]), int(f[1])])
    return np.array(new_faces)

def generate_fusing_data_from_lines(vertices, edges, min_distance = float("inf")):
    """
    Generate fusing data
    """
    fusing_data = [0] * len(vertices)
    for e in edges:
        idx = find_line_closest_points(vertices, e, min_distance)#generate_line_in_slope_intercept_form(e[0], e[1]), min_distance)
        for i in idx: 
            fusing_data[i] = 1
    return fusing_data                   

def find_line_closest_points(point_set, line, min_distance):
    """
    Finds the closest points in a given point set to a given line, and returns their indexes.
    """

    slope = (line[1][1] - line[0][1]) / (line[1][0] - line[0][0])
    y_intercept = line[0][1] - slope * line[0][0]

    closest_points = []
    for i, point in enumerate(point_set):
        distance = abs(slope * point[0] - point[1] + y_intercept) / (slope**2 + 1)**0.5
        if distance < min_distance:
            closest_points.append(i)
    
    return closest_points

def find_point_closest_point(point_set, point):
    """
    Finds the closest point in a given point set to a given point, and returns its index.
    """
    distances = np.sqrt(np.sum((point_set - point)**2, axis=1)) 
    return np.argmin(distances)
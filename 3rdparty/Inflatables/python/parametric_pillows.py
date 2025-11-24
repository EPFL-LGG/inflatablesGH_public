import numpy as np
def squareWithVerticalChannels(channelHeightFrac, numChannels):
    # border
    pts = [[0, 0], [1, 0], [1, 1], [0, 1]]
    edges = [(0, 1), (1, 2), (2, 3), (3, 0)]
    channelSpacing = 1.0 / numChannels
    ybot = (1.0 - channelHeightFrac) * 0.5
    ytop = 1.0 - ybot
    for wall in range(1, numChannels):
        x = wall * channelSpacing
        edges.append(tuple([len(pts), len(pts) + 1]))
        pts += [[x, ybot], [x, ytop]]
    return pts, edges

def parallelTubes(N, h, d, w, triArea):
    """
    Generate a rectangular sheet mesh with N parallel tubes of height h with
    tube width w and fusing curve width d.
    Unlike the other functions in this module, this creates an actual triangle mesh
    """
    tubeXCoords = [0, w] # tube 0
    nonemptyFuseRegion = d > 1e-10
    for i in range(1, N):
        l_old = tubeXCoords[-1]
        if nonemptyFuseRegion: tubeXCoords += [l_old + d, l_old + d + w] # fuse gap + tube i 
        else:                  tubeXCoords += [l_old + w] # just tube i
    pts = [[x, y] for y in [0, h] for x in tubeXCoords]
    numBottomPts = len(tubeXCoords)
    edges  = [(i, i + 1) for i in range(numBottomPts - 1)]                      # bottom edges
    edges += [(i, i + 1) for i in range(numBottomPts, 2 * numBottomPts - 1)]    # top edges
    edges += [(i, numBottomPts + i) for i in range(numBottomPts)]               # vertical edges

    import wall_generation
    m, fuseMarkers, fuseSegments = wall_generation.triangulate_channel_walls(pts, edges, triArea)

    if (nonemptyFuseRegion):
        import sheet_meshing
        fakeSignedDistance = lambda x: w - np.fmod(x, d + w) # will be negative inside the wall region
        pointSetLiesInWall = lambda X: np.mean(fakeSignedDistance(X[:, 0])) < 0.0
        m, iwv, iwbv = sheet_meshing.remeshWallRegions(m, fuseMarkers, fuseSegments, pointSetLiesInWall)
    else:
        iwv = iwbv = fuseMarkers

    return m, iwv, iwbv

def concentricCircles(numChannels, numSegments):
    pts, edges = [], []
    for wall in range(1, numChannels + 1): # last wall is actually the pillow boundary
        # The innermost circle's diameter should match the channel width, so its radius
        # is half the channel spacing (treat it like a half-channel).
        r = (0.5 + (wall - 1)) / (numChannels - 0.5)
        
        angles = np.linspace(0, 2 * np.pi, numSegments, endpoint = False)
        newPtIdxs = list(range(len(pts), len(pts) + len(angles)))
        pts += [[r * np.cos(t), r * np.sin(t)] for t in angles]
        edges += zip(newPtIdxs, newPtIdxs[1:] + [newPtIdxs[0]])
    return pts, edges

def circle(resolution):
    t = np.linspace(0, 2 * np.pi, resolution, endpoint = False)
    pts   = list(np.column_stack((np.cos(t), np.sin(t))))
    edges = list(zip(range(resolution), list(range(1, resolution)) + [0]))
    return pts, edges

def radialChannels(numSectors, numSectorSubdivisions, circleResolution):
    pts, edges = circle(circleResolution)

    def subdivide(nsubdiv, angleStart, angleEnd, rStart, rEnd):
        if (nsubdiv == 0): return
        angleMid = 0.5 * (angleStart + angleEnd)
        alpha = 0.5
        rMid     = (1 - alpha) * rStart + alpha * rEnd

        edges.append((len(pts), len(pts) + 1))
        pts.append([rMid * np.cos(angleMid), rMid * np.sin(angleMid)])
        pts.append([rEnd * np.cos(angleMid), rEnd * np.sin(angleMid)])

        subdivide(nsubdiv - 1, angleStart, angleMid, rMid, rEnd)
        subdivide(nsubdiv - 1, angleMid  , angleEnd, rMid, rEnd)

    for sector in range(numSectors):
        sectorBegin = (2 * np.pi * sector) / numSectors
        sectorEnd =   (2 * np.pi * (sector + 1)) / numSectors

        rStart = 0.06
        rEnd = 0.95

        edges.append((len(pts), len(pts) + 1))
        pts.append([rStart * np.cos(sectorBegin), rStart * np.sin(sectorBegin)])
        pts.append([rEnd   * np.cos(sectorBegin), rEnd   * np.sin(sectorBegin)])

        subdivide(numSectorSubdivisions, sectorBegin, sectorEnd, rStart, rEnd)

    return pts, edges

from matplotlib import pyplot as plt

# alpha: angle between spiral tangent and axis vector dtheta (not dr)
def logSpiralPlot(alpha = 70.0, radius = 1.0, minDist = 0.05, margin = 0.05, edgeLength = 0.02):
    alpha_rad = np.deg2rad(alpha)
    b = np.tan(alpha_rad)

    # Note: the logarithmic spiral is self-similar, so its scale is irrelevant.
    # We therefore use the unit-scaled logarithmic spiral, given in polar coordinates by:
    #   r = e^(b theta)
    # We evaluate the spiral at evenly spaced points along its arclength.
    sqrtTerm = np.sqrt(1 + 1 / (b * b))
    rForTheta      = lambda th: np.exp(b * th)
    thetaForR      = lambda r: np.log(r) / b
    thetaForArclen = lambda s: (1.0 / b) * (np.log(s + sqrtTerm) - np.log((sqrtTerm)))
    arclenForTheta = lambda th: sqrtTerm * (np.exp(b * th) - 1.0)

    def thetasForRadiusInterval(rmin, rmax):
        smin, smax = map(lambda r: arclenForTheta(thetaForR(r)), [rmin, rmax])
        nsubdiv = int(np.round((smax - smin) / edgeLength))
        return thetaForArclen(np.linspace(smin, smax, nsubdiv))

    pts, edges = circle(int(np.round(2 * radius * np.pi / edgeLength)))

    def generate_points(rs, thetas, rotation = 0):
        return np.column_stack((rs * np.cos(thetas + rotation), rs * np.sin(thetas + rotation)))

    numSectors = 1
    # Draw walls (spiral arms) dividing the circle into numSectors sectors
    # for numSectors in 2, 4, 8, ...
    while True:
        # We approximate the channel thickness by multiplying the channel's
        # sector angle by the normal velocity of spiral arm (wall) points as
        # the arms are rotated at unit angular velocity.
        #       thickness ~= (2 * pi / numSectors) * r * sin(alpha)
        # where sin(alpha) is the (constant) angle between the curve's normal
        # and the radial axis vector dr.
        # Then we can solve for the minimum radius such that the thickness is >= minDist:
        #       (2 * pi / numSectors) * r * sin(alpha) >= minDist   ==>
        #       r >= (minDist * numSectors) / (2 * pi * sin(alpha))
        rmin = (minDist * numSectors) / (2 * np.pi * np.sin(alpha_rad))
        # The following version reproduces the original Matlab behavior (but
        # leads to tightly spaced channels for small alpha)
        # rmin = max(minDist / 2, rmin) if numSectors > 2 else minDist / 2
        rmin = max(minDist / 2, rmin)
        rmax = radius - margin
        if (rmin > rmax - edgeLength): break # admissible channel walls have shrunk below the target edge length
        thetas = thetasForRadiusInterval(rmin, rmax)
        rvalues = rForTheta(thetas)

        for arm in range(numSectors):
            if ((numSectors > 1) and (arm % 2 == 0)): continue # even arms have already been drawn by previous passes
            newPts = list(generate_points(rvalues, thetas, arm * (2 * np.pi) / numSectors))
            if (len(newPts) >= 2): # At least two points must be added to form a segment
                ptOffset = len(pts)
                pts += newPts
                for i in range(len(newPts) - 1):
                    edges.append((ptOffset + i, ptOffset + i + 1))
        numSectors *= 2

    return pts, edges

def bentArc(length, width, curvature, numArcSegments, includeStart = True, includeEnd = True):
    r = 1.0 / curvature
    thetaLen = length / r
    
    thetaStart, thetaEnd = 0, thetaLen
    if (not includeStart): thetaStart += thetaLen / numArcSegments
    if (not includeEnd):   thetaEnd   -= thetaLen / numArcSegments
    numArcSegments -= includeStart + includeEnd
    
    thetas = np.linspace(thetaStart, thetaEnd, numArcSegments + 1)[:, np.newaxis]
    rOuter, rInner = r + width / 2, r - width / 2
    circlePts = np.concatenate((np.cos(thetas),
                                np.sin(thetas)), axis=1)
    pts = np.concatenate((rOuter * circlePts, rInner * circlePts[::-1]), axis=0)
    pts -= [r, 0]
    
    arcPtIdxs = np.arange(numArcSegments + 1)
    edges = list(zip(arcPtIdxs, arcPtIdxs[1:]))

    if includeEnd: edges.append((numArcSegments, numArcSegments + 1))
    innerPtOffset = numArcSegments + 1
    edges += zip(innerPtOffset + arcPtIdxs, innerPtOffset + arcPtIdxs[1:])
    if includeStart: edges.append((len(pts) - 1, 0))

    return pts, edges

# Get a vector field by sampling the contraction direction at "pts"
def bentArcContractionDirection(curvature, pts):
    r = 1.0 / curvature;
    pts = pts + [r, 0] # Undo translation
    return pts / np.linalg.norm(pts, axis=1)[:, np.newaxis]

# numLinks: number of links in the chain of circular arcs
# (if numLinks == 1, we get an s-shaped curve)
def bentArchChain(numLinks, length, width, curvature, numArcSegments):
    r = 1.0 / curvature
    thetaLen = length / r
    pts, edges = np.empty((0, 2)), []
    rotMat = lambda x: np.array([[np.cos(x), -np.sin(x)],[np.sin(x), np.cos(x)]])
    terminals = None
    for l in range(numLinks):
        # Transform the old points so that the new chain will connect at the origin
        pts = (pts + [r, 0]) @ rotMat(-thetaLen).transpose() - [r, 0]
        # Reflect the chain for every added link to curve in the opposite direction.
        # Note: we make no guarantees on the curve orientation, so we needn't reverse
        #       any edges!
        pts *= [-1, 1]
        link_pts, link_edges = bentArc(length, width, curvature, numArcSegments, l == 0, l == numLinks - 1)
        offset = pts.shape[0]
        pts = np.concatenate((pts, link_pts), axis=0)
        # Connect with the terminals of the previous chain link (if one exists)
        if (terminals is not None):
            edges += [(terminals[1], offset), (terminals[0], offset + len(link_pts) - 1)]

        edges += [(i + offset, j + offset) for i, j in link_edges]
        terminals = (offset + numArcSegments - 1 + (l > 0), offset + numArcSegments + (l > 0))
        
    return pts, edges

from scipy.special import comb
def tilted_sine_function(x, n, freq, amplitude):
    k = (np.arange(n) + 1).reshape(n, 1)
    top = np.array([comb(2*n, n - i) for i in k])
    factor = top / comb(2*n, n) / k
    kx = np.einsum('ik,jk->ij', x.reshape(len(x), 1), k)
    return (np.einsum('ij,jk->ik', np.sin(kx * freq) * amplitude, factor)).flatten()

epsilon = 1e-7

def trim_vertices(pts, edges, lower_x, upper_x):
    import networkx as nx
    G = nx.Graph() 
    G.add_nodes_from(np.arange(len(pts)))
    G.add_edges_from(edges)
    for i, pt in enumerate(pts):
        if pt[0] < lower_x - epsilon or pt[0] > upper_x + epsilon:
            G.remove_node(i)
    Gcc = sorted(nx.connected_components(G), key=len, reverse=True)
    if (nx.number_connected_components(G) > 1):
        G0 = G.subgraph(Gcc[0])
    else:
        G0 = G
    updated_edges = list(G0.edges)
    new_pidx = list(G0.nodes)
    index_map = {new_pidx[i] : i for i in range(len(new_pidx))}
    new_edges = [(index_map[e[0]], index_map[e[1]]) for e in (updated_edges)]
    new_pts = np.array(pts)[new_pidx]
    return new_pts, new_edges

from edge_edge_intersection_helper import Point, doIntersect, line, intersection

import numpy.linalg as la
def get_intermediate_point(sp, ep, target_length):
    sp = np.array(sp)
    ep = np.array(ep)
    total_length = la.norm(ep - sp)
    if (target_length is not None):
        num_edges = int(np.round(total_length / target_length))
    else:
        return [sp, ep]
    point_list = [sp]
    for i in range(num_edges)[1:]:
        lam = 1.0 / num_edges * i
        point_list.append(lam * sp + (1 - lam) * ep)
    point_list.append(ep)

    return point_list

def connect_ps_pt(pts, edges, sp, tp, target_length):
    # Connect point sp and point pt with new points.
    pl = get_intermediate_point(pts[sp], pts[tp], target_length)
    prev = sp
    for j in range(len(pl) - 2):
        edges.append([prev, len(pts) + j])
        prev = len(pts) + j
    edges.append([prev, tp])

    for j in range(len(pl) - 2):
        pts.append(pl[j + 1])
    return pts, edges 

def add_lines(p1, q1, ip1, iq1, pts, edges, target_length):
#     First compute all intersections of the input line segments with all existing edges.
    intersections = []
    for i, e in enumerate(edges):
        p2 = Point(*pts[e[0]])
        q2 = Point(*pts[e[1]])
        if (doIntersect(p1, q1, p2, q2)):
            int_pt = intersection(line(pts[e[0]], pts[e[1]]), line(p1.coord, q1.coord))
            intersections.append(list(int_pt) + [i])

    # Sort the intersections by their y coordinates and remove duplicates as well as the first and last ones (the ones that share the end points of the input line segment)
    intersections.sort(key = lambda x : x[1])
    intersections = np.array(intersections)
    unique_keys, indices = np.unique(np.array(intersections)[:,1].round(6), return_index=True)
    edges = list(edges)
    prev = ip1
    remove_edge = []

    for ie in intersections[indices][1:]:
        e = edges[int(ie[-1])]
        curr_tl = target_length

        int_pt = ie[:2]
        if np.allclose(int_pt, pts[e[0]]):
            pts[e[0]] = int_pt
            pts, edges = connect_ps_pt(pts, edges, prev, e[0], curr_tl)

            prev = e[0]


        elif (np.allclose(int_pt, pts[e[1]])):
            pts[e[1]] = int_pt
            pts, edges = connect_ps_pt(pts, edges, prev, e[1], curr_tl)

            prev = e[1]

        else:
            remove_edge.append(e)

            pl = get_intermediate_point(pts[prev], int_pt, curr_tl)
            for i in range(len(pl) - 1):
                edges.append([prev, len(pts) + i])
                prev = len(pts) + i
            # Need to connect the new intersection point with the original end points of the crossing edges. 
            edges.append([prev, e[0]])
            edges.append([prev, e[1]])

            for i in range(len(pl) - 1):
                pts.append(pl[i + 1])

    # edges.append([prev, iq1])

    for e in remove_edge:
        edges.remove(e)
    return pts, edges


def get_perioidic_mesh(points, edges):
    import numpy.linalg as la
    # Build graph from mesh.

    import networkx as nx
    G = nx.Graph()
    G.add_edges_from(edges)
    for e in G.edges():
        G[e[0]][e[1]]['weight'] = la.norm(points[e[0]] - points[e[1]])

    # Find boundary path segments.
    bot_vx_idxs = np.where(np.isclose(points[:, 1], 0))
    bot_vxs = points[bot_vx_idxs]
    top_vx_idxs = np.where(np.isclose(points[:, 1], 1))
    top_vxs = points[top_vx_idxs]

    top_idxs_sorted = top_vx_idxs[0][top_vxs[:, 0].argsort()]
    bot_idxs_sorted = bot_vx_idxs[0][bot_vxs[:, 0].argsort()]

    top_path = nx.shortest_path(G, top_idxs_sorted[0], top_idxs_sorted[-1], weight='weight')
    bot_path = nx.shortest_path(G, bot_idxs_sorted[0], bot_idxs_sorted[-1], weight='weight')

    left_path = nx.shortest_path(G, bot_idxs_sorted[0], top_idxs_sorted[0], weight='weight')
    right_path = nx.shortest_path(G, bot_idxs_sorted[-1], top_idxs_sorted[-1], weight='weight')

    def copy_path(source_path, dest_path, shift, shift_direction, points):
        def remove_intermediate_vx_in_path(path):
            for idx in path[1:-1]:
                keep_vx = False
                for n in G.neighbors(idx):
                    if n not in path:
                        keep_vx = True
                        break
                if not keep_vx:
                    G.remove_node(idx)
        def replace_old_path_with_new_path(old_path, new_path_points, points):        
            prev = old_path[0]
            for vx in new_path_points[1:]:
                index = np.where(np.isclose(points, vx).all(axis = 1))[0]
                if (len(index) > 0):
                    index = index[0]
                    G.add_edge(prev, index)
                    prev = index
                else:
                    G.add_edge(prev, len(points))
                    prev = len(points)
                    points = np.vstack([points, vx])
            return points

        # Copy source_path to dest_path after shifting.
        remove_intermediate_vx_in_path(dest_path)
        dest_path_points = points[source_path]
        dest_path_points[:, shift_direction] += shift
        points = replace_old_path_with_new_path(dest_path, dest_path_points, points)
        return points

    points = copy_path(left_path, right_path, (points[bot_idxs_sorted[-1]] - points[bot_idxs_sorted[0]])[0], 0, points)
    points = copy_path(bot_path, top_path, 1, 1, points)

    # Export new line segments constraints.
    new_pidx = list(G.nodes)
    updated_edges = list(G.edges)
    index_map = {new_pidx[i] : i for i in range(len(new_pidx))}
    new_edges = [(index_map[e[0]], index_map[e[1]]) for e in (updated_edges)]
    new_pts = np.array(points)[new_pidx]
    return new_pts, new_edges


def sinusoid_raw(N, h, d, w, triArea, numSegments, freq, amplitude, tilt_n, clipping_lN = None, clipping_uN = None, target_length = None):
    """
    Generate a rectangular sheet mesh with N parallel tubes of height h with
    tube width w and fusing curve width d.
    Unlike the other functions in this module, this creates an actual triangle mesh
    """
    nonemptyFuseRegion = d > 1e-10
    tubeXCoords = [0, d, w + d] if nonemptyFuseRegion else [0, w] # tube 0
    for i in range(1, N):
        l_old = tubeXCoords[-1]
        if nonemptyFuseRegion: tubeXCoords += [l_old + d, l_old + d + w] # fuse gap + tube i 
        else:                  tubeXCoords += [l_old + w] # just tube i

    # tubeXCoords += [tubeXCoords[-1] + d]
    y = np.linspace(0, h, numSegments)
    offset = tilted_sine_function(y, tilt_n, freq, amplitude)
    pts = [[x + offset[i], y[i]] for i in np.arange(numSegments) for x in tubeXCoords]

    numBottomPts = len(tubeXCoords)

    if (target_length is not None):
        edges = []

        for i in range(numBottomPts - 1):
            if i % 2 == 0 and nonemptyFuseRegion:
                edges.append((i, i+1))
            else:
                pts, edges = connect_ps_pt(pts, edges, i, i+1, target_length)

        for i in range((numSegments - 1) * numBottomPts, numSegments * numBottomPts - 1):
            if (i % 2 != 0 and nonemptyFuseRegion):
                edges.append((i, i+1))
            else:
                pts, edges = connect_ps_pt(pts, edges, i, i+1, target_length)

    else:
        edges  = [(i, i + 1) for i in range(numBottomPts - 1)]  # bottom edges
        edges += [(i, i + 1) for i in range((numSegments - 1) * numBottomPts, numSegments * numBottomPts - 1)]    # top edges

    for i in range(numBottomPts):
        for j in range(numSegments - 1):
            edges.append((i + j * numBottomPts, i + (j + 1) * numBottomPts))             # vertical edges

    if (clipping_lN is not None and clipping_uN is not None):
        l1 = ((clipping_lN * (2 if nonemptyFuseRegion else 1), clipping_lN * (2 if nonemptyFuseRegion else 1) + numBottomPts * (numSegments - 1)))
        l2 = ((clipping_uN * (2 if nonemptyFuseRegion else 1), clipping_uN * (2 if nonemptyFuseRegion else 1) + numBottomPts * (numSegments - 1)))
        l1ps, l1pe = Point(*pts[l1[0]]), Point(*pts[l1[1]])
        l2ps, l2pe = Point(*pts[l2[0]]), Point(*pts[l2[1]])

        def x_in_range(x):
            return x >= l1ps.x and x <= l2ps.x
        pts, edges = add_lines(l1ps, l1pe, l1[0], l1[1], pts, edges, target_length)

        pts, edges = add_lines(l2ps, l2pe, l2[0], l2[1], pts, edges, target_length)

        pts, edges = trim_vertices(pts, edges, l1ps.x, l2ps.x)
    return pts, edges

def sinusoid(N, h, d, w, triArea, numSegments, freq, amplitude, tilt_n, clipping_lN = None, clipping_uN = None, target_length = None, use_periodic = True, epsilon = 1e-3):
    pts, edges = sinusoid_raw(N, h, d, w, triArea, numSegments, freq, amplitude, tilt_n, clipping_lN, clipping_uN, target_length)
    import wall_generation
    m, fuseMarkers, fuseSegments = wall_generation.triangulate_channel_walls(pts, edges, triArea, flags="")
    if (use_periodic):
        pts, edges = get_perioidic_mesh(m.vertices()[:, :2], fuseSegments)
        m, fuseMarkers, fuseSegments = wall_generation.triangulate_channel_walls(pts, edges, triArea, flags="Y")

    if (d > 1e-10):
        import sheet_meshing
        fakeSignedDistance = lambda x: d - np.fmod(x, d + w) # will be positive inside the wall region
        pointSetLiesInWall = lambda X: np.mean(fakeSignedDistance(X[:, 0] - tilted_sine_function(X[:, 1], tilt_n, freq, amplitude).flatten())) >= -1e-16
        m, iwv, iwbv = sheet_meshing.remeshWallRegions(m, fuseMarkers, fuseSegments, pointSetLiesInWall, markBoundaryAsWall = True)

        pointLiesInWall= lambda X: np.logical_or(fakeSignedDistance(X[:, 0] - tilted_sine_function(X[:, 1], tilt_n, freq, amplitude).flatten()) >= -epsilon, w + fakeSignedDistance(X[:, 0] - tilted_sine_function(X[:, 1], tilt_n, freq, amplitude).flatten()) <= epsilon)

    else:        
        iwv = iwbv = fuseMarkers

        fakeSignedDistance = lambda x: np.fmod(x, w) # will be positive inside the wall region
        pointLiesInWall= lambda X: np.logical_or(fakeSignedDistance(X[:, 0] - tilted_sine_function(X[:, 1], tilt_n, freq, amplitude).flatten()) < epsilon, fakeSignedDistance(X[:, 0] - tilted_sine_function(X[:, 1], tilt_n, freq, amplitude).flatten()) > w - epsilon)


    if (use_periodic):
        points = m.vertices()[:, :2]
        iwv = np.array(iwv)
        iwv = pointLiesInWall(points)
        iwv = pointLiesInWall(points)

        iwbv = iwv

    return m, iwv, iwbv

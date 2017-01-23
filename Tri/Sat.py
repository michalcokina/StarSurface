# implementation of SAP algorithm in 2d
# detection of collision of two triangles
# if triangles are touched, return False too

import numpy as np
PRECISION = 10

def projection(v, t):
    t, v = np.array(t), np.array(v)
    return (np.dot(v, t) / np.dot(t, t)) * t


def to_tangential(tn, nn, pt):
    m = [[np.dot(np.array([1.0, 0.0]), tn), np.dot(np.array([0.0, 1.0]), tn)],
         [np.dot(np.array([1.0, 0.0]), nn), np.dot(np.array([0.0, 1.0]), nn)]]
    return np.dot(m, pt)


def separating_axis_theorem(tri1, tri2):
    # each edge of triangle tri1 and tri2
    edges = np.array([[tri1[i], tri1[i + 1]] for i in range(-1, 2)] + [[tri2[i], tri2[i + 1]] for i in range(-1, 2)])

    for edge in edges:
        # tangent vector of testing line
        tangent = (edge[1] - edge[0]) / np.linalg.norm(edge[1] - edge[0])
        # normal vector of testing line
        normal = -tangent[1], tangent[0]

        # we flip tangent and normal vector (we will project all points of each triangle on the normal line of
        # current edge)
        tangent, normal = normal, tangent

        # projection of each vertex on testing line
        projection_t1 = [to_tangential(tangent, normal, projection(vertex, tangent)).tolist() for vertex in tri1]
        projection_t2 = [to_tangential(tangent, normal, projection(vertex, tangent)).tolist() for vertex in tri2]

        # new x coo (y should be zero now)
        projection_t1_newx, projection_t2_newx = list(zip(*projection_t1))[0], list(zip(*projection_t2))[0]

        # maximal length projected of each triangle
        projection_edge1, projection_edge2 = [round(min(projection_t1_newx), PRECISION),
                                              round(max(projection_t1_newx), PRECISION)], \
                                             [round(min(projection_t2_newx), PRECISION),
                                              round(max(projection_t2_newx), PRECISION)]

        projection_edge = [projection_edge1, projection_edge2]
        projection_edge.sort(key=lambda x: x[0])

        # if intervals connected in point or separated return True
        if projection_edge[0][1] <= projection_edge[1][0]:
            return True

    return False


def intersection(t1, t2):
    return not separating_axis_theorem(t1, t2)

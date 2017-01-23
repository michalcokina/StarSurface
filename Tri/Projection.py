import numpy as np

PRECISION = 10


def line(u, a, b):
    # line equation
    # x = a  + u * (b - a)
    a, b = np.array(a), np.array(b)
    return a + (u * (b - a))


def line_t(u, a, t):
    """
    :rtype: ndarray
    """
    a, t = np.array(a), np.array(t)
    return a + (u * t)


def plane(n, p, a, b):
    # plane equation
    # (x - p)  n = 0
    # coo of plane point x = [a, b, c]
    # get a, b and return [a, b, c] if n[2] != 0 etc.
    n, p = np.array(n), np.array(p)
    if n[2] != 0:
        c = (np.dot(p, n) - (a * n[0]) - (b * n[1])) / n[2]
        return a, b, c
    elif n[1] != 0:
        c = (np.dot(p, n) - (a * n[0]) - (b * n[2])) / n[1]
        return a, c, b
    elif n[0] != 0:
        c = (np.dot(p, n) - (a * n[1]) - (b * n[2])) / n[0]
        return c, a, b
    else:
        return False


def planeline_intersection(a, b, p, n):
    """
    :rtype: ndarray
    """
    a, b, p, n = np.array(a), np.array(b), np.array(p), np.array(n)
    if np.dot((b - a), n) != 0:
        u_intersection = np.dot((p - a), n) / np.dot((b - a), n)
        return line(u_intersection, a, b)
    return False


def to_2d_plane(s1, s2, s3, pt):
    # s1, s2, s3 - orthonormal vectors of reference frame
    s1, s2, s3, pt = np.array(s1), np.array(s2), np.array(s3), np.array(pt)
    m = [[np.dot(np.array([1., 0., 0.]), s1), np.dot(np.array([0., 1., 0.]), s1), np.dot(np.array([0., 0., 1.]), s1)],
         [np.dot(np.array([1., 0., 0.]), s2), np.dot(np.array([0., 1., 0.]), s2), np.dot(np.array([0., 0., 1.]), s2)],
         [np.dot(np.array([1., 0., 0.]), s3), np.dot(np.array([0., 1., 0.]), s3), np.dot(np.array([0., 0., 1.]), s3)]]
    return np.dot(m, pt)


def get_projection_plane(n, p):
    projection_plane = [p]
    a, b = np.array([[1.5, 1.5], [0.0, 1.0]])
    # a[i], b[i] - values in equation (p - p0) . n = 0
    # p . n - p0 . n = 0
    # p0 . n - (a[i] * n[0] + b[i] * n[1] + c[i] * n[2]) = 0
    # we have to find c[i], then point in plane defined by (n, p) is [a[i], b[i], c[i]]
    for i in range(0, 2):
        projection_plane.append(np.array(plane(n, p, a[i], b[i])))
    return np.array(projection_plane)


def projection(n, p, pp, x, dim="2d"):
    # n - normal vector define plane (outer vector of surface)
    # p - point define plane
    # pp - projection point parameter (define projection point), # parameter "u" in eq. [x, y] = A + u (B - A)
    # x - array like, points for projections to the plane
    # dim - if 2d return only first two coordinate of points (anyway third should be close to zero, in ideal case: 0)

    # return variable
    inplane = []
    # trasnform to np.array
    n, p, x = np.array(n), np.array(p), np.array(x)

    # create distant prejection point
    projection_point = line_t(pp, p, n)
    # create projection plane
    projection_plane = get_projection_plane(n, p)
    # unit vector of projection plane
    # normal vector of plane will be still "z" - axis of new coo system
    ref1 = (projection_plane[0] - projection_plane[1]) / (np.linalg.norm(projection_plane[0] - projection_plane[1]))
    ref3 = n / np.linalg.norm(n)
    ref2 = np.cross(ref1, ref3) / np.linalg.norm(np.cross(ref1, ref3))

    for x_val in x:
        intersection = planeline_intersection(projection_point, x_val, p, n)
        p_inplane = to_2d_plane(ref1, ref2, ref3, intersection - p)
        if dim == "2d":
            inplane.append([p_inplane[0], p_inplane[1]])
            continue
        inplane.append(p_inplane)

    return inplane

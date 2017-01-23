import inspect
import numpy as np
import math as m
import warnings
import operator

def center_of_mass(faces=None):
    return np.array([(tri[0] + tri[1] + tri[2]) / 3.0 for tri in faces])

def average_spacing(data=None, neighbours=6):
    if type(data) != type(np.array([])):
        data = np.array(data)

    # ########################## !!!!!!!!!!!!!!!!! IMPORTANT
    # neighbours is variable match same variable in v cgal function
    from scipy.spatial import distance
    dis = distance.cdist(data, data, 'euclidean')
    tot = 0
    for line in dis:
        tot += np.sort(line)[1:1 + neighbours].sum() / (neighbours + 1)
    return tot / dis.shape[0]


def distances_matrix(vertices):
    from scipy.spatial import KDTree
    tree = KDTree(vertices)
    dist = tree.sparse_distance_matrix(tree, max_distance=10.0)
    dist_array = {}

    for i in range(0, len(vertices)):
        for j in range(0, len(vertices)):
            if i == j:
                continue
            dist_array[(i, j)] = dist[(i, j)]

    sorted_dist_array = dict(sorted(dist_array.items(), key=operator.itemgetter(1)))
    return sorted_dist_array

def distance_to(to, matrix):
    dist, ml = {}, int(np.ceil(np.sqrt(len(matrix))))
    for i in range(ml):
        if i == to:
            continue
        dist[(to, i)] = matrix[(to, i)]
    sorted_dist_array = sorted(dist.items(), key=operator.itemgetter(1))
    return sorted_dist_array


def lineno():
    """Returns the current line number in our program."""
    return inspect.currentframe().f_back.f_lineno


def cartesian_to_spheric(vector=np.arange(3, dtype=np.float), degrees=False):
    spherical = np.arange(3, dtype=np.float)
    spherical[0] = np.sqrt(vector[0] ** 2 + vector[1] ** 2 + vector[2] ** 2)  # vypocet vzdialenosti "r"

    np.seterr(divide='raise', invalid='raise')
    # vypocet a osetrenie vypoctu pre azimutalny uhol
    try:
        spherical[1] = np.arcsin(
            vector[1] / (np.sqrt(vector[0] ** 2 + vector[1] ** 2)))  # vypocet azimutalneho (rovinneho) uhla
    except:
        spherical[1] = 0.0
    # vypocet a osetrenie vypoctu pre polarny uhol
    try:
        spherical[2] = np.arccos(vector[2] / (spherical[0]))  # vypocet polarneho (elevacneho) uhla
    except:
        spherical[2] = 0.0

    np.seterr(divide='print', invalid='print')

    if vector[0] < 0:
        spherical[1] = np.pi - spherical[1]

    if not degrees:
        return spherical
    else:
        return np.array([spherical[0], np.degrees(spherical[1]), np.degrees(spherical[2])])


def spheric_to_cartesian(vector=np.arange(3, dtype=np.float)):
    cartesian = np.arange(3, dtype=np.float)
    cartesian[0] = vector[0] * m.cos(vector[1]) * m.sin(vector[2])
    cartesian[1] = vector[0] * m.sin(vector[1]) * m.sin(vector[2])
    cartesian[2] = vector[0] * m.cos(vector[2])
    return cartesian


def rotate(
        angle=0,
        vector=np.arange(3, dtype=np.float),
        inverse=False,
        axis="x"):
    matrix = np.arange(9, dtype=np.float).reshape((3, 3))
    vector = np.array(vector)
    if axis == "x":
        matrix[0][0], matrix[0][1], matrix[0][2] = 1, 0, 0
        matrix[1][0], matrix[1][1], matrix[1][2] = 0, np.cos(angle), - np.sin(angle)
        matrix[2][0], matrix[2][1], matrix[2][2] = 0, np.sin(angle), np.cos(angle)
        if inverse:
            matrix[1][2], matrix[2][1] = np.sin(angle), - np.sin(angle)
    if axis == "y":
        matrix[0][0], matrix[0][1], matrix[0][2] = np.cos(angle), 0, np.sin(angle)
        matrix[1][0], matrix[1][1], matrix[1][2] = 0, 1, 0
        matrix[2][0], matrix[2][1], matrix[2][2] = - np.sin(angle), 0, np.cos(angle)
        if inverse:
            matrix[2][0], matrix[0][2] = + np.sin(angle), - np.sin(angle)
    if axis == "z":
        matrix[0][0], matrix[0][1], matrix[0][2] = np.cos(angle), - np.sin(angle), 0
        matrix[1][0], matrix[1][1], matrix[1][2] = np.sin(angle), np.cos(angle), 0
        matrix[2][0], matrix[2][1], matrix[2][2] = 0, 0, 1
        if inverse:
            matrix[0][1], matrix[1][0] = + np.sin(angle), - np.sin(angle)
    return np.dot(matrix, vector)


def arbitrary_rotate(theta, omega=None, vector=None, degrees=False):
    # omega - lubovolny vektor okolo ktoreho sa ma rotovat
    # theta - uhol o kolko rotovat
    # vector - vector ktory sa bude rotovat okolo omega

    omega = np.array(omega) / np.linalg.norm(np.array(omega))
    if degrees:
        theta = np.radians(theta)

    matrix = np.arange(9, dtype=np.float).reshape((3, 3))

    matrix[0][0] = (np.cos(theta)) + (omega[0] ** 2 * (1. - np.cos(theta)))
    matrix[0][1] = (omega[0] * omega[1] * (1. - np.cos(theta))) - (omega[2] * np.sin(theta))
    matrix[0][2] = (omega[1] * np.sin(theta)) + (omega[0] * omega[2] * (1. - np.cos(theta)))

    matrix[1][0] = (omega[2] * np.sin(theta)) + (omega[0] * omega[1] * (1. - np.cos(theta)))
    matrix[1][1] = (np.cos(theta)) + (omega[1] ** 2 * (1. - np.cos(theta)))
    matrix[1][2] = (- omega[0] * np.sin(theta)) + (omega[1] * omega[2] * (1. - np.cos(theta)))

    matrix[2][0] = (- omega[1] * np.sin(theta)) + (omega[0] * omega[2] * (1. - np.cos(theta)))
    matrix[2][1] = (omega[0] * np.sin(theta)) + (omega[1] * omega[2] * (1. - np.cos(theta)))
    matrix[2][2] = (np.cos(theta)) + (omega[2] ** 2 * (1. - np.cos(theta)))

    return np.dot(matrix, vector)


def empty(var, debug=False):
    # if var is numpy arr type
    if type(var) == type(np.array([])):
        if debug:
            print("Variable type is <numpy array>")
        if var.size == 0:
            return True
    # if var is python tuple
    elif type(var) == type(()):
        if debug:
            print("Variable type is <python tuple>")
        if len(var) == 0:
            return True
    # if var is python list type
    elif type(var) == type([]):
        if debug:
            print("Variable type is <python list>")
        if np.array(var).size == 0:
            return True
    # if var is dictionary type
    elif type(var) == type({}):
        if debug:
            print("Variable type is <python dict>")
        if var is {}:
            return True
    elif type(var) == type(True):
        if debug:
            print("Variable type is <bool>")
        if not var:
            return True
    elif var is None:
        if debug:
            print("Variable type is <NoneType>")
        return True
    elif type(var) == type("foo"):
        if debug:
            print("Variable type is <string>")
        if var is "" or var == "0":
            return True
    else:
        try:
            if np.isnan(var):
                if debug:
                    print("Variable type is <numpy.nan>")
                return True
            else:
                if debug:
                    print("Variable type is <number>")
                if var == 0:
                    return True
        except:
            print("Variable type is invalid")
            return True
    return False


def numpy_array_is_empty_verbose(arr, function_name, class_name, var_name, verbose, line):
    warnings.simplefilter(action="ignore", category=FutureWarning)

    try:
        if arr is None:
            if verbose:
                print("ValueError: In reference: " + str(
                      class_name) + ", function: " + str(function_name) + ", line: " + line + ". Variable `" + str(
                      var_name) + "` is None.")
            return True

        if arr.size == 0:
            if verbose:
                print("EmptyVariableError: In class: " + str(
                      class_name) + ", function: ",
                      str(function_name) + ", line: " + line + ". Variable `" + str(var_name) + "` is empty.")
            return True
        else:
            return False
    except:
        if verbose:
            print("ValueError: In class: " + str(
                  class_name) + ", function: " + str(function_name) + ", line: " + line + ". Variable `" + str(
                  var_name) + "` is invalid.")
        return True


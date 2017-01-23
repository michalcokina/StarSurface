import scipy
from scipy.spatial import KDTree
import operator
import numpy as np
import Tri.Projection as Projection
import Tri.Sat as Sat
import copy
import matplotlib.path as mpltpath
import os
import Tri.Iostream as Io
from Star.StarCreator import Geometry

def separate(faces=None, separation=None):
    var = {"primary": [], "secondary": []}
    com = center_of_mass(faces=np.array(faces))

    for f in faces:
        pass

    return True


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


def average_spacing(data=None, neighbours=6):
    if type(data) != type(np.array([])):
        data = np.array(data)

    # ########################## !!!!!!!!!!!!!!!!! IMPORTANT
    # neighbours is variable match same variable in v cgal function
    from scipy.spatial import distance
    dis = scipy.spatial.distance.cdist(data, data, 'euclidean')
    tot = 0
    for line in dis:
        tot += np.sort(line)[1:1 + neighbours].sum() / (neighbours + 1)
    return tot / dis.shape[0]


def center_of_mass(faces=None):
    return np.array([(tri[0] + tri[1] + tri[2]) / 3.0 for tri in faces])


def distances_matrix(vertices):

    avsp = average_spacing(data=vertices, neighbours=5)
    tree = KDTree(vertices)
    dist = tree.sparse_distance_matrix(tree, max_distance=avsp * 4)
    dist_array = {}

    # import Plot.Plot as Plt
    # Plt.plot_3d(vertices=[[vertices[354], vertices[528]]], faces=None, normals_view=False, points_view=True,
    #             faces_view=False, point_color="r", point_size=5.0, face_color="c", edge_color="k")

    for i in range(0, len(vertices)):
        for j in range(0, len(vertices)):
            if i == j:
                continue
            if dist[(i, j)] == 0.0:
                # dist_array[(i, j)] = np.inf
                continue
            dist_array[(i, j)] = dist[(i, j)]

    sorted_dist_array = dict(sorted(dist_array.items(), key=operator.itemgetter(1)))
    del(dist_array, dist)

    return sorted_dist_array


# def distance_to(to, matrix):
#     dist, ml = {}, int(np.ceil(np.sqrt(len(matrix))))
#     for i in range(ml):
#         try:
#             if i == to:
#                 continue
#             dist[(to, i)] = matrix[(to, i)]
#         except KeyError:
#             continue
#
#     sorted_dist_array = sorted(dist.items(), key=operator.itemgetter(1))
#
#     return sorted_dist_array


# def distance_to(to, matrix):
#     dist = {}
#     for v in matrix:
#         if v[0] == to:
#             dist[(to, v[1])] = matrix[(to, v[1])]
#     return sorted(dist.items(), key=operator.itemgetter(1))


def array_mask(array, mask):
    return np.array([array[idx] for idx in mask])


def sim2tri(vertices, simplices):
    return [[vertices[sim[0]], vertices[sim[1]], vertices[sim[2]]] for sim in simplices]


def find_common_simplices(tsim1, tsim2):
    return [s for s in tsim1 if s in tsim2]


def find_not_common_simplices(simplices, common):
    return [s for s in simplices if s not in common]


def delaunay2d(vertices):
    from scipy.spatial import Delaunay as Del
    return Del(np.array(vertices, dtype="float64")).simplices


def delaunay3d(vertices):
    from scipy.spatial import Delaunay as Del
    vertices = np.array(vertices)
    convex_hull = Del(vertices).convex_hull
    return convex_hull, vertices[convex_hull]


# noinspection PyTypeChecker
def trispot(vertices=None, norms=None, spots=None, binary_morph=None):
    simplices_map = {"primary": {}, "secondary": {}}
    vertices_t, norms_t = copy.copy(vertices), copy.copy(norms)

    # import Plot.Plot as Plt
    # Plt.plot_3d(vertices=[spots["primary"][0]["vertices"]], faces=None, normals_view=False, points_view=True,
    #             faces_view=False, point_color="r", point_size=5.0, face_color="c", edge_color="k")

    # Plt.plot_3d(vertices=[np.concatenate((spots["primary"][0]["vertices"], vertices["primary"]), 0)],
    #             faces=None, normals_view=False, points_view=True,
    #             faces_view=False, point_color=
    #             [np.concatenate((["b"] * len(spots["primary"][0]["vertices"]), ["r"] * len(vertices["primary"])), 0)],
    #             point_size=5.0, face_color="r", edge_color="k")

    for t_object in ["primary", "secondary"]:
        for i in range(0, len(vertices[t_object])):
            simplices_map[t_object][i] = {"type": "t_object", "clearance": True, "enum": -1}

    if "t_object" in locals():
        del t_object

    for t_object in ["primary", "secondary"]:
        # average spacing of stars
        avsp = average_spacing(data=vertices_t[t_object], neighbours=6)
        for spot, spot_num in list(zip(spots[t_object], range(0, len(spots[t_object])))):
            avsp_spot = average_spacing(data=spot["vertices"], neighbours=6)
            simplices_test, vertices_test = [], []

            try:
                # find nerest points to spot alt center
                tree = KDTree(vertices_t[t_object])
                distance, ndx = tree.query(spot["alt_center"], k=len(vertices_t[t_object]))

                for dist, i in list(zip(distance, ndx)):
                    if dist > spot["size"] + (0.25 * avsp):
                        # break, because distancies are ordered by size
                        break

                    # distance exeption for spots points
                    if simplices_map[t_object][i]["type"] == "spot":
                        if dist > spot["size"] + (0.05 * avsp_spot):
                            continue

                    # norms 0 belong to alt center
                    if np.dot(spot["norms"][0], norms_t[t_object][i]) > 0:
                        simplices_test.append(i)
                # simplices of target object for testing whether point lying inside or not of spot boundary
                simplices_test = list(set(simplices_test))

                # test if index to remove from all current vertices belong to any spot
                clearance_enum, spot_indices, star_indices = [], [], []
                for st in simplices_test:
                    if simplices_map[t_object][st]["type"] == "spot":
                        clearance_enum.append(simplices_map[t_object][st]["enum"])
                        spot_indices.append(st)
                    else:
                        star_indices.append(st)
                clearance_enum = list(set(clearance_enum))

                # tento kod je momentalne odstrihnuty
                # myslienka tohto kodu je taka, ze sa spravi priemet do normalovej roviny
                # a odstrania sa len body, ktory su vo vnutri n-uholnika vytvoreneho z hranicnych
                # bodov skrvny;
                # je to vypnute, pretoze cez vzdialenosti ked sa to riesi, tak to skvrny vyzeraju krajsie
                if len(clearance_enum) > 0 and False:
                    vertices_test = [vertices_t[t_object][i] for i in spot_indices]
                    # projection to tangential plane
                    vp = Projection.projection(n=spot["norms"][0], p=spot["alt_center"], pp=-20, x=vertices_test,
                                               dim="2d")
                    sp = Projection.projection(n=spot["norms"][0], p=spot["alt_center"], pp=-20, x=spot["boundary"],
                                               dim="2d")

                    # import matplotlib.pyplot as plt
                    # plt.scatter(list(zip(*sp))[0], list(zip(*sp))[1])
                    # plt.scatter(list(zip(*vp))[0], list(zip(*vp))[1], c="r")
                    # plt.axis('equal')
                    # plt.show()

                    # remove star points lying inside spot boundary
                    bb_path = mpltpath.Path(np.array(sp))
                    p_inside = bb_path.contains_points(vp)

                    exclude = []
                    if True in p_inside:
                        exclude = [ndx for statement, ndx in list(zip(p_inside, spot_indices)) if statement == True]
                    simplices_test = np.concatenate((star_indices, exclude), axis=0)

                    # repeated clearance enum (opakuje sa to preto, lebo ak pouzijem toto projektovanie, tak sa
                    # v pohode stane, ze tie vertexy, ktore by sa mali vyhodit sa vyhodit uz nemaju
                    # a moze sa stat nahodou, ze zrazu uz z nejakej skvrny nebude treba odstranovat ziaden bod)
                    clearance_enum, spot_indices = [], []
                    for st in simplices_test:
                        if simplices_map[t_object][st]["type"] == "spot":
                            clearance_enum.append(simplices_map[t_object][st]["enum"])
                            spot_indices.append(st)
                    clearance_enum = list(set(clearance_enum))

                # change clearance status (if any vertex is removed from spot, then spot clearance status is
                # set to False, otherwise, if entire spot is used, clearance status is set to default True)
                for enum in clearance_enum:
                    for i in simplices_map[t_object]:
                        if simplices_map[t_object][i]["enum"] == enum:
                            simplices_map[t_object][i]["clearance"] = False

                # vertices, norms and simplices_map update
                if len(simplices_test) != 0:
                    vertices_temp, norms_temp = [], []
                    map_temp, j = {}, 0
                    for i, vertex, norm in list(zip(range(0, len(vertices_t[t_object])), vertices_t[t_object],
                                                    norms_t[t_object])):
                        if i in simplices_test:
                            continue
                        vertices_temp.append(vertex)
                        norms_temp.append(norm)

                        map_temp[j] = {"type": simplices_map[t_object][i]["type"],
                                       "clearance": simplices_map[t_object][i]["clearance"],
                                       "enum": simplices_map[t_object][i]["enum"]}
                        j += 1

                    shift = len(vertices_temp)
                    for i, vertex, norm in list(zip(range(shift, shift + len(spot["vertices"])), spot["vertices"],
                                                    spot["norms"])):
                        vertices_temp.append(vertex)
                        norms_temp.append(norm)
                        map_temp[i] = {"type": "spot", "clearance": True, "enum": spot_num}

                    vertices_t[t_object] = copy.deepcopy(vertices_temp)
                    simplices_map[t_object] = copy.deepcopy(map_temp)
                    norms_t[t_object] = copy.deepcopy(norms_temp)

                    del (vertices_temp, map_temp, norms_temp)

            except IndexError:
                print("IndexError")
                continue
    if "t_object" in locals():
        del t_object

    # # visualization
    # v, s = [], []
    # v = [x for x, i in list(zip(vertices_t["primary"], simplices_map["primary"])) if
    #      simplices_map["primary"][i]["type"] == "t_object" and x[0] > 0.1]
    # s = [x for x, i in list(zip(vertices_t["primary"], simplices_map["primary"])) if
    #      simplices_map["primary"][i]["type"] == "spot"]
    #
    # import Plot.Plot as Plt
    # Plt.plot_3d(vertices=[v, s], faces=None, normals_view=False, points_view=True, faces_view=False,
    #             point_color=[["r"] * len(v), ["b"] * len(s)], point_size=30.0, face_color=None, edge_color="k",
    #             elev=45, azim=45, save=False)

    # triangulation process
    if binary_morph == "detached" or binary_morph == "semi-contact":
        tri = {"primary": delaunay3d(vertices=vertices_t["primary"]),
               "secondary": delaunay3d(vertices=vertices_t["secondary"])}

    elif binary_morph == "over-contact":
        import pickle

        # !!! vertices_t and similar (norms_t...): key "system" do not contain correct values,
        # use "primary" and "secondary"

        if True:

            tri_class = {"primary": Tri(vertices=vertices_t["primary"], norms=norms_t["primary"]),
                         "secondary": Tri(vertices=vertices_t["secondary"], norms=norms_t["secondary"])}

            tri_class["primary"].triangulate()
            tri_class["secondary"].triangulate()

            tri = {"primary": [tri_class["primary"].simplices(), tri_class["primary"].hull()],
                   "secondary": [tri_class["secondary"].simplices(), tri_class["secondary"].hull()]}

            del(tri_class["primary"], tri_class["secondary"])

            # pickle.dump(tri, open("tri.pickle", "wb"))

        # tri = pickle.load(open("tri.pickle", "rb"))

        face_orientation = {"primary": Geometry.face_orientation_a(face=np.array(tri["primary"][1]), t_object="primary",
                                                                   actual_distance=0.0),
                            "secondary": Geometry.face_orientation_a(face=np.array(tri["secondary"][1]),
                                                                     t_object="secondary",
                                                                     actual_distance=1.0)}

        # import Plot.Plot as Plt
        # Plt.plot_3d(vertices=None, faces=[tri["primary"][1]], normals_view=False, points_view=False, faces_view=True,
        #             point_color=None, point_size=30.0, face_color="c", edge_color="k",
        #             elev=45, azim=45, save=False)

        rm_indices = {"primary": [], "secondary": []}

        for t_object in ["primary", "secondary"]:
            for f, fo, ndx in list(zip(tri[t_object][1], face_orientation[t_object], range(0, len(tri[t_object][1])))):
                # if x coo is more then center coo of primary star
                if t_object == "primary" and (f[0][0] > 0.0 and f[1][0] > 0.0 and f[2][0] > 0.0):
                    # we have to test face orientation, because should happpend, there is wrong face on neck (is inside
                    # of neck, (create cover of it))
                    if abs(Geometry.angle(u=fo, v=np.array([1., 0., 0.]))) <= np.pi / 180.0:
                        rm_indices[t_object].append(ndx)
                elif t_object == "secondary" and (f[0][0] < 1.0 and f[1][0] < 1.0 and f[2][0] < 1.0):
                    if abs(Geometry.angle(u=fo, v=np.array([-1., 0., 0.]))) <= np.pi / 180.0:
                        rm_indices[t_object].append(ndx)

            s = [i[1] for i in enumerate(tri[t_object][0]) if i[0] not in rm_indices[t_object]]
            t = [i[1] for i in enumerate(tri[t_object][1]) if i[0] not in rm_indices[t_object]]

            tri[t_object][0], tri[t_object][1] = s, t

    # separation spot and t_object
    model = {"primary":
                 {"t_object": [], "spots": {}},
             "secondary":
                 {"t_object": [], "spots": {}}}

    # spots enumeration and structure preparation
    # we need this tlanslator, because spot "enum" should be different to list index of this spot in variable
    # "model", e.g. if entire spot lies inside of another, it is removed;
    # model variable is based on existence of spots in simplices_map and if spot is removed, there is missing
    # index in simplices_map, e.g. we have spots with enum (0, 1, 3) because 2 is for examle inside of 3 and
    # then spots in model have anyway indices (0, 1, 2) so we have to trransale 0 -> 0, 1 -> 1 and 2 -> 3 or
    # in reverse way depend on current situation we handle with.
    spots_indices_map = {"primary": list(set([simplices_map["primary"][i]["enum"] for i in simplices_map["primary"]
                                              if simplices_map["primary"][i]["type"] == "spot"])),
                         "secondary": list(
                             set([simplices_map["secondary"][i]["enum"] for i in simplices_map["secondary"]
                                  if simplices_map["secondary"][i]["type"] == "spot"]))}

    rev_spots_indices_map = {"primary": None, "secondary": None}
    spot_candidate = {"primary": {"vertices": {}, "com": {}, "3rd_enum": {}},
                      "secondary": {"vertices": {}, "com": {}, "3rd_enum": {}}}

    # i[0] - num, i[1] - spot enum
    for t_object in ["primary", "secondary"]:
        spots_indices_map[t_object] = {i[1]: i[0] for i in enumerate(spots_indices_map[t_object])}
        rev_spots_indices_map[t_object] = {spots_indices_map[t_object][i]: i for i in spots_indices_map[t_object]}

        if len(spots_indices_map[t_object]) != 0:
            for i in range(len(spots_indices_map[t_object])):
                model[t_object]["spots"][i] = []
                spot_candidate[t_object]["vertices"][i] = []
                spot_candidate[t_object]["com"][i] = []
                spot_candidate[t_object]["3rd_enum"][i] = []
    if "t_object" in locals():
        del t_object

    # next will be t_o used as t_object for shorter notation, so if somewhere in comment is mentioned t_object as
    # variable, I mean t_o variable
    for t_o in ["primary", "secondary"]:
        for s, t in list(zip(tri[t_o][0], tri[t_o][1])):
            # test if each point belongs to spot
            if simplices_map[t_o][s[0]]["type"] == "spot" and \
                            simplices_map[t_o][s[1]]["type"] == "spot" and \
                            simplices_map[t_o][s[2]]["type"] == "spot":

                # if each point belongs to one spot, then it is for sure face of that spot
                if simplices_map[t_o][s[0]]["enum"] == \
                        simplices_map[t_o][s[1]]["enum"] == \
                        simplices_map[t_o][s[2]]["enum"]:
                    model[t_o]["spots"][spots_indices_map[t_o][simplices_map[t_o][s[0]]["enum"]]].append(
                        np.array(t).tolist())

                else:
                    # if at least one of points of face belongs to different spot, we have to test
                    # which one of those spots current face belongs to
                    reference = None

                    # variable trd_enum is enum index of 3rd corner of face;

                    if simplices_map[t_o][s[-1]]["enum"] == simplices_map[t_o][s[0]]["enum"]:
                        reference = simplices_map[t_o][s[-1]]["enum"]
                        trd_enum = simplices_map[t_o][s[1]]["enum"]
                    elif simplices_map[t_o][s[0]]["enum"] == simplices_map[t_o][s[1]]["enum"]:
                        reference = simplices_map[t_o][s[0]]["enum"]
                        trd_enum = simplices_map[t_o][s[-1]]["enum"]
                    elif simplices_map[t_o][s[1]]["enum"] == simplices_map[t_o][s[-1]]["enum"]:
                        reference = simplices_map[t_o][s[1]]["enum"]
                        trd_enum = simplices_map[t_o][s[0]]["enum"]

                    if reference is not None:
                        spot_candidate[t_o]["vertices"][spots_indices_map[t_o][reference]].append(t)
                        spot_candidate[t_o]["com"][spots_indices_map[t_o][reference]].append(center_of_mass([t])[0])
                        spot_candidate[t_o]["3rd_enum"][spots_indices_map[t_o][reference]].append(trd_enum)

            # if at least one of points belongs to t_object, then it is for sure t_object face
            elif simplices_map[t_o][s[0]]["type"] == "t_object" or \
                simplices_map[t_o][s[1]]["type"] == "t_object" or \
                simplices_map[t_o][s[2]]["type"] == "t_object":

                model[t_o]["t_object"].append(np.array(t).tolist())
            else:
                model[t_o]["t_object"].append(np.array(t).tolist())

        # test unclear spots faces
        # ! enum of spot is same as list index in spots[t_object]
        # ! spots_indices_map{spot_enum (same as spots[t_object] index of spot): model_number (model variable)}
        # ! rev_spots_indices_map{model_number: spot_enum (same as spots[t_object] index of spot)}

        if len(spot_candidate[t_o]["com"]) != 0:
            for spot_reference in spot_candidate[t_o]["com"]:
                # translate spot_reference from model list order to spot_number (index from original spots list)
                rev_spot_reference = rev_spots_indices_map[t_o][spot_reference]

                # get center and size of current spot candidate
                center, size = spots[t_o][rev_spot_reference]["alt_center"], spots[t_o][rev_spot_reference]["size"]

                # compute distance of all center of mass of triangles of current spot candidate to center of
                # this candidate
                dists = [np.linalg.norm(np.array(com) - np.array(center))
                         for com in spot_candidate[t_o]["com"][spot_reference]]

                # test if dist is smaller as size;
                # if dist is smaller, then current face belongs to spots
                # otherwise face belongs to t_object itself

                for dist in enumerate(dists):
                    dist = list(dist)
                    app = np.array(spot_candidate[t_o]["vertices"][spot_reference][dist[0]]).tolist()
                    if dist[1] < size:
                        model[t_o]["spots"][spot_reference].append(app)
                        continue

                    # make the same computation for 3rd vertex of face
                    # it might be confusing, but spot candidate is spot where 2 of 3 vertex of face belong to,
                    # the 3rd index belongs to another (neighbour) spot and it has to be tested
                    # too, if face finally do not belongs to spot candidate;
                    # reference to this another spot is saved in 3rd_enum of spot_candidate[t_o][3rd_enum][ref]
                    # dictionary, referenced by "ref", where "ref" is "spot_reference" variable from above
                    else:
                        # in 3rd_enum is directly stored reference to default spots list which comes as
                        # function parameterclear

                        trd_rev_spot_reference = spot_candidate[t_o]["3rd_enum"][spot_reference][dist[0]]
                        trd_spot_reference = spots_indices_map[t_o][trd_rev_spot_reference]

                        trd_center = spots[t_o][trd_rev_spot_reference]["alt_center"]
                        trd_size = spots[t_o][trd_rev_spot_reference]["size"]
                        com = spot_candidate[t_o]["com"][spot_reference][dist[0]]
                        dist[1] = np.linalg.norm(np.array(com) - np.array(trd_center))

                        if dist[1] < trd_size:
                            model[t_o]["spots"][trd_spot_reference].append(app)
                            continue
                        else:
                            model[t_o]["t_object"].append(app)

    if "model" in locals():
        return model

def cgal_triangulation(
        normals=None,
        points=None,
        min_triangle_angle=0.349066,  # a lower bound on the minimum angle in degrees of the surface mesh facets.
        # an upper bound on the radius of surface Delaunay balls. A surface Delaunay ball is a ball circumscribing
        # a facet, centered on the surface and empty of vertices. Such a ball exists for each facet of the current
        # surface mesh. Indeed the current surface mesh is the Delaunay triangulation of the current sampling restricted
        # to the surface which is just the set of facets in the three dimensional Delaunay triangulation of the sampling
        # that have a Delaunay surface ball.
        max_triangle_size=1,
        # an upper bound on the center-center distances of the surface mesh facets. The center-center distance of
        # a surface mesh facet is the distance between the facet circumcenter and the center of its surface
        # Delaunay ball.
        surface_aproximation_error=0.375,
        to_average_spacing=5,
        # to_average_spacing: k nearest neighbours used in computing of point cloud average spacing
        # best value used to be 5

        # zasadne pouizvat hodnotu 5 pre average spacing, pretoze to zodpoveda rozumenj vzdialenosti
        # dvoch blizkych bodoch na povrchu
        input_path="input.xyz"
):
    # if error osccured during pcd file saving, function will be terminated and will return False
    if not Io.save_cgal_pcd_with_normals(filename="input.xyz", filedir=input_path + "/", points=points,
                                         normals=normals):
        return False

    # Min triangle angle
    # Max triangle size w.r.t. point set average spacing;
    # Surface Approximation error w.r.t. point set average spacing
    # to average spacing

    # tato krkolomna konstrukcia ''np.degrees([min_triangle_angle])[0]'' tu je preto, aby pycharm nepapuloval

    terminal = os.path.dirname(os.path.realpath(__file__)) + '/poisson_reconstruction ' + \
        str(np.degrees([min_triangle_angle])[0]) + ' ' + str(max_triangle_size) + ' ' + \
        str(surface_aproximation_error) + ' ' + str(to_average_spacing) + ' ' + str(input_path) + "/input.xyz"

    # # same as above, but without verbosity
    # terminal = os.path.dirname(os.path.realpath(__file__)) + '/poisson_reconstruction ' + \
    #            str(np.degrees([min_triangle_angle])[0]) + ' ' + str(max_triangle_size) + ' ' + \
    #            str(surface_aproximation_error) + ' ' + str(to_average_spacing) + ' ' + str(input_path) +
    #            "/input.xyz" + ' > /dev/null 2>&1'

    os.system(terminal)  # returns the exit status
    os.remove(input_path + '/input.xyz')

    if os.path.isfile(input_path + '/output.off'):
        # load_cgal_3d_mesh will return tuple of (faces, vertices, simplices)
        tri = Io.load_cgal_3d_mesh(filename='output.off', filedir=input_path + "/")
        os.remove(input_path + '/output.off')
    else:
        return False

    if not empty(var=tri):
        return tri
    else:
        return False


class Tri:
    def __init__(self, vertices=None, norms=None):
        self.vertices = np.array(vertices)
        self.norms = np.array(norms)
        self._simplices = None
        self._hull = None

    def triangulate(self, srange=14, sindex=0):

        # srange (simplex range) (int)
        # sindex (start index) (int)

        # distances matrix
        d_mtx = distances_matrix(vertices=self.vertices)
        print("d_mtx done...")
        # first face
        fface = self.get_face(srange=srange, sindex=sindex, dmatrix=d_mtx)[0]

        # variables declar.
        tsimplices, tsimplices_queue, simplices_list = [fface], [fface], []

        while True:
            # break if computing queue is empty
            if len(tsimplices_queue) < 1:
                break

            # current vertices of triangle for finding next one
            current_tsimplices = tsimplices_queue[0]
            # run over all vertices in current triangle
            for simplex in current_tsimplices:
                new_tsimplices, intersc = None, None
                # trying to optim. doesn't work
                # if simplex not in simplices_list:
                #     simplices_list.append(simplex)
                # else:
                #     continue

                # new simplices
                nwsimplices = self.get_face(srange=srange, sindex=simplex, dmatrix=d_mtx)

                # iterate through every triangle created in current projection
                for simm in nwsimplices:
                    if simm in tsimplices:
                        # such simm already exists
                        continue

                    intersc = self.intersection(swhat=simm, swhere=tsimplices, srange=srange, dmatrix=d_mtx)

                    if intersc[0] is True:
                        use = False

                        # in many case happend, if triangles are in intersection, they have a two common points
                        # (see pic.) and then we should test, if we flip them, if there will be solution
                        #           D                      ________D
                        # A\\     /\                     A \\      \
                        #   \  \ /  \           ===>        \  \    \
                        #    \  / \  \                       \    \  \
                        #    B\/______\ C                    B\_______\ C

                        common = find_common_simplices(simm, intersc[1])
                        # if triangles have common simplices different to two, than continue,
                        # something like that cannot be flipped
                        if len(common) != 2:
                            continue

                        # otherwise, we have to create flipped case;
                        # in this case, there should be comnon only one point, so we are able to use [0]
                        ncmn_1 = find_not_common_simplices(simm, common)[0]
                        ncmn_2 = find_not_common_simplices(intersc[1], common)[0]
                        new_tsimplices = [sorted([ncmn_1, ncmn_2, common[0]]), sorted([ncmn_1, ncmn_2, common[1]])]

                        # we have to check, if this new_tsimplices is not in intersection with other,
                        # already existing triangles
                        for new_ts in new_tsimplices:
                            # there exist 2 triangles of course, but one of them should already exists in tsimplices
                            # array,
                            # should already exist, if there is intersection with it
                            if new_ts in tsimplices:
                                continue

                            # the other one, has to be checked "manualy";
                            # we will use the same code as above
                            intersc = self.intersection(swhat=new_ts, swhere=tsimplices, srange=srange, dmatrix=d_mtx)
                            if intersc[0] is False:
                                use, simm = True, new_ts
                                # we should break, because second one already exists
                                break

                    else:
                        use = True

                    if use is True:
                        tsimplices.append(simm)
                        tsimplices_queue.append(simm)
                        break
                del (intersc, nwsimplices, new_tsimplices)

            idx = tsimplices_queue.index(current_tsimplices)
            del tsimplices_queue[idx]

        self._simplices = tsimplices
        del (d_mtx, fface, tsimplices, tsimplices_queue, simplices_list)

    def nearest_to(self, to, srange, matrix):
        # closer_map is link between array closer and array self.vertices (in meaning of indices)
        closer, runner = [to], srange
        # get distance from self.vertices[sindex] to all other
        dist = self.distance_to(to, matrix)

        # append to closer srange closest points if dot product of normals are greather than 0
        for i in range(0, runner):
            # vertex of nearest point
            nsimplex = dist[i][0][1]
            if np.dot(self.norms[to], self.norms[nsimplex]) > 0:
                closer.append(nsimplex)
                continue
            # if dot product is negative, runner is incremented
            runner += 1

        return closer

    def intersection(self, swhat, swhere, srange, dmatrix):
        # vyhlada najblizsie trojuholniky z swhere k trojuholniku swhat a zisti ci sa medzi sebou tieto
        # najblizsie nepretinaju s swhat trojuholnikom

        # vezmeme jeden z bodov na trojuholniku swhat
        s = swhat[0]
        # nearest simplices to s point
        # nsimplices = sorted([dist[i][0][1] for i in range(0, srange)])
        closer = self.nearest_to(to=s, srange=srange, matrix=dmatrix)
        # simplices of nearest triangle from swhere to swhat
        nsimplices = []
        # iterate each close point
        for simplex in closer:
            # iterate each triangle already appended to final list (swhere)
            for face in swhere:
                if simplex in face and face not in nsimplices and face is not swhat:
                    nsimplices.append(face)

        # projection of triangles to 2d plane
        unique = np.unique(nsimplices)

        ps_projection = Projection.projection(n=self.norms[s],
                                              p=self.vertices[s],
                                              pp=-10,
                                              x=array_mask(array=self.vertices, mask=unique),
                                              dim="2d")

        # uprava aby sa to dalo podla indexov z povodneho pola precitat, ked sa pouzite funkcia sim2tri
        ps_projection_dict = {}
        for u, item in list(zip(unique, ps_projection)):
            ps_projection_dict[u] = item

        # project swhat triangle only to the same plane as other points (ps_projection)
        s_projection = Projection.projection(n=self.norms[s],
                                             p=self.vertices[s],
                                             pp=-10,
                                             x=array_mask(array=self.vertices, mask=swhat),
                                             dim="2d")

        # triangle 2d coo to comparsion;
        # nsimplices is in same order as twhere, so if twhere[1] is in intersect with twhat, than simplices for
        # twhere[1] are stored in nsimplices[1]
        twhere = sim2tri(vertices=ps_projection_dict, simplices=nsimplices)
        twhat = s_projection

        for i, t in list(zip(range(len(twhere)), twhere)):
            # tu pokracovat, treba s novym (mozno dalsim v zozname) trojuholnikom spravit intersekcny test s otanymi
            # teda s t-ckami v tomto loope a v momente, ked sa najde intersekcia vratit False, ze sa pretina a teda
            # je zlym kadnidatom na zaradenie, vlastne nemoze byt urcite zaradeny

            if Sat.intersection(twhat, t):
                return True, nsimplices[i]

        return False,

    def get_face(self, srange, sindex, dmatrix):
        closer = self.nearest_to(to=sindex, srange=srange, matrix=dmatrix)
        # projection of closer points to sindex normal plane
        ps_projection = Projection.projection(n=self.norms[sindex],
                                              p=self.vertices[sindex],
                                              pp=-10,
                                              x=array_mask(array=self.vertices, mask=closer),
                                              dim="2d")

        # for plotting chage dim = 2d to dim = 3d
        # import objects.Plot as Plt
        # Plt.plot_3d(vertices=[ps_projection], normals_view=False, points_view=True, faces_view=False,
        #             point_color="r", point_size=5.0, face_color="c", edge_color="k", verbose=True)
        #
        # import matplotlib.pyplot as plt
        # plt.scatter(list(zip(*ps_projection))[0], list(zip(*ps_projection))[1], s=50.0)
        # plt.axis('equal')
        # plt.savefig(filename="o.pdf")

        # first simplices (first triangle)
        fsimplices, return_simplices = delaunay2d(vertices=ps_projection), []

        for item in fsimplices:
            # zero in item, because sindex is added as firts value (zero valu in programmer term) to closer
            if 0 in item:
                return_simplices.append(sorted([closer[item[0]], closer[item[1]], closer[item[2]]]))
        return return_simplices

    def distance_to(self, to, matrix):
        dist, ml = {}, int(np.ceil(np.sqrt(len(self.vertices) ** 2 - len(self.vertices))))
        # povodne bolo toto za ml, funkcia bola staticmethod a v distance_matrix bolo odkomentovane np.inf riadok
        # int(np.ceil(np.sqrt(len(matrix))))
        for i in range(ml):
            try:
                if i == to:
                    continue
                dist[(to, i)] = matrix[(to, i)]
            except KeyError:
                continue

        sorted_dist_array = sorted(dist.items(), key=operator.itemgetter(1))

        return sorted_dist_array

    def load_csv(self, filename, separator, inplace=False):
        vertices, norms = [], []
        with open(filename, "r") as f:
            while True:
                data = f.readline()
                if not data:
                    break
                data = data.split(separator)
                try:
                    vertices.append([data[i] for i in range(0, 3)])
                    norms.append([data[i] for i in range(3, 6)])
                except IndexError:
                    break
        vertices, norms = np.array(vertices, dtype="float64"), np.array(norms, dtype="float64")

        if inplace:
            self.vertices = vertices
            self.norms = norms
        return vertices, norms

    def simplices(self):
        return self._simplices

    def hull(self):
        h = [[self.vertices[simplex[0]], self.vertices[simplex[1]], self.vertices[simplex[2]]] for simplex in
             self._simplices]
        self._hull = np.array(h[:]).tolist()
        return self._hull
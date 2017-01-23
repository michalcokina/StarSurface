"""
###   In this document we show some examples how to create a spots on star surface.                                  ###
###   !!! For runing of any part of code, please change False to True in condition.                                  ###
"""

# imports ##############################################################################################################

import numpy as np
import os
from Star.StarCreator import StarCreator
import Tri.Triangulation as T
import Plot.Plot as Plt

# /imports #############################################################################################################

# colors list ##########################################################################################################
colors = ["#eae310", "#ea8110", "y", "r", "g", "#eabb10"]
# /colors list #########################################################################################################

# detached system ######################################################################################################
primary = {"mass": 2.0, "potential": 4.0, "phis": 5, "thetas": 10}
secondary = {"mass": 1.0, "potential": 4.0, "phis": 5, "thetas": 10}
orbit = {"period": 2.0, "eccentricity": 0.0, "inclination": np.pi, "argument_of_periastron": np.pi / 2.0}
system = StarCreator(primary=primary, secondary=secondary, orbit=orbit)
system.create_model()
binary_vertices = system.get_vertices()
binary_norms = system.get_norms()
from Star.StarCreator import Geometry


# single spot on the primary component
# ------------------------------------
if False:
    system.create_spots(meta=[
        # {"lon": 0.5, "lat": 0.5, "diameter": 0.5, "steps_azimuthal": 5,
        #  "steps_radial": 5, "t_object": "primary"},
        #
        # {"lon": 0.5, "lat": 0.5, "diameter": 0.2, "steps_azimuthal": 5,
        #  "steps_radial": 5, "t_object": "primary"},

        {"lon": 0.5, "lat": 1.1, "diameter": 0.8, "steps_azimuthal": 4,
         "steps_radial": 4, "t_object": "primary"}
    ])

    spots = system.get_spots()

    # v = [x for x in binary_vertices["primary"] if x[0] > 0.1]
    # s = spots["primary"][0]["vertices"]
    #
    # Plt.plot_3d(vertices=[v, s], faces=None, normals_view=False, points_view=True, faces_view=False,
    #             point_color=[["r"] * len(v), ["b"] * len(s)], point_size=30.0, face_color=None, edge_color="k",
    #             elev=45, azim=45, save=False)

    model = T.trispot(vertices=binary_vertices, norms=binary_norms, spots=spots,
                      binary_morph=system.get_binary_morph())

    # visualization
    faces = [model["primary"]["t_object"]]
    color = [[colors[-1]] * len(model["primary"]["t_object"])]

    for i in model["primary"]["spots"]:
        faces.append(model["primary"]["spots"][i])
        # noinspection PyTypeChecker
        color.append([colors[i]] * len(model["primary"]["spots"][i]))

    Plt.plot_3d(vertices=None, faces=faces, normals_view=False, points_view=False, faces_view=True,
                point_color=None, point_size=5.0, face_color=color, edge_color="k",
                elev=45, azim=45, save=False)


# two spots on the secondary component
# ------------------------------------
if False:
    system.create_spots(meta=[
        {"lon": 1.5, "lat": 1.0, "diameter": 0.4, "steps_azimuthal": 4,
         "steps_radial": 4, "t_object": "secondary"},

        {"lon": 0.5, "lat": 1.1, "diameter": 0.5, "steps_azimuthal": 4,
         "steps_radial": 4, "t_object": "secondary"}
    ])

    spots = system.get_spots()
    model = T.trispot(vertices=binary_vertices, norms=binary_norms, spots=spots,
                      binary_morph=system.get_binary_morph())

    # visualization
    faces = [model["secondary"]["t_object"]]
    color = [[colors[-1]] * len(model["primary"]["t_object"])]

    for i in model["secondary"]["spots"]:
        faces.append(model["secondary"]["spots"][i])
        # noinspection PyTypeChecker
        color.append([colors[i]] * len(model["secondary"]["spots"][i]))

    Plt.plot_3d(vertices=None, faces=faces, normals_view=False, points_view=False, faces_view=True,
                point_color=None, point_size=5.0, face_color=color, edge_color="k",
                elev=45, azim=45, save=False, x_range=[0.5, 1.5], y_range=[-0.5, 0.5], z_range=[-0.5, 0.5])


# two spots on the primary component, overlaping spots
# ----------------------------------------------------
if False:
    system.create_spots(meta=[
        {"lon": 0.5, "lat": 1.1, "diameter": 0.7, "steps_azimuthal": 5,
         "steps_radial": 5, "t_object": "primary"},
        {"lon": 0.8, "lat": 1.4, "diameter": 0.5, "steps_azimuthal": 5,
         "steps_radial": 5, "t_object": "primary"}
    ])

    spots = system.get_spots()
    model = T.trispot(vertices=binary_vertices, norms=binary_norms, spots=spots,
                      binary_morph=system.get_binary_morph())

    # visualization
    faces = [model["primary"]["t_object"]]
    color = [[colors[-1]] * len(model["primary"]["t_object"])]

    for i in model["primary"]["spots"]:
        faces.append(model["primary"]["spots"][i])
        # noinspection PyTypeChecker
        color.append([colors[i]] * len(model["primary"]["spots"][i]))

    Plt.plot_3d(vertices=None, faces=faces, normals_view=False, points_view=False, faces_view=True,
                point_color=None, point_size=5.0, face_color=color, edge_color="k",
                elev=45, azim=45, save=False)


# two spots on the primary component, overlaping spots (posibility for umbra and penumbra)
# bigger has to be first
# ----------------------------------------------------------------------------------------
if False:
    system.create_spots(meta=[
        {"lon": 0.5, "lat": 1.1, "diameter": 0.7, "steps_azimuthal": 5,
         "steps_radial": 5, "t_object": "primary"},
        {"lon": 0.5, "lat": 1.1, "diameter": 0.4, "steps_azimuthal": 5,
         "steps_radial": 5, "t_object": "primary"}
    ])

    spots = system.get_spots()
    model = T.trispot(vertices=binary_vertices, norms=binary_norms, spots=spots,
                      binary_morph=system.get_binary_morph())

    # visualization
    faces = [model["primary"]["t_object"]]
    color = [[colors[-1]] * len(model["primary"]["t_object"])]

    for i in model["primary"]["spots"]:
        faces.append(model["primary"]["spots"][i])
        # noinspection PyTypeChecker
        color.append([colors[i]] * len(model["primary"]["spots"][i]))

    Plt.plot_3d(vertices=None, faces=faces, normals_view=False, points_view=False, faces_view=True,
                point_color=None, point_size=5.0, face_color=color, edge_color="k",
                elev=45, azim=45, save=False)

# /detached system #####################################################################################################

# over-contact system ##################################################################################################
primary = {"mass": 3.0, "potential": 2.48, "phis": 6, "thetas": 12}
secondary = {"mass": 1.2, "potential": 2.48, "phis": 6, "thetas": 12}
orbit = {"period": 0.7, "eccentricity": 0.0, "inclination": np.pi, "argument_of_periastron": np.pi / 2.0}
system = StarCreator(primary=primary, secondary=secondary, orbit=orbit)
system.create_model()
binary_vertices = system.get_vertices()
binary_norms = system.get_norms()

# dot surface visualization
# -------------------------
if False:
    Plt.plot_3d(vertices=[binary_vertices["system"]], faces=None, normals_view=False, points_view=True,
                faces_view=False, point_color="r", point_size=5.0, face_color="c", edge_color="k")
# single spot on the primary component
# ------------------------------------
if True:
    system.create_spots(meta=[
        {"lon": np.pi / 2.0, "lat": 0.8, "diameter": 0.7, "steps_azimuthal": 5,
         "steps_radial": 5, "t_object": "primary"}
    ])

    spots = system.get_spots()
    model = T.trispot(vertices=binary_vertices, norms=binary_norms, spots=spots,
                      binary_morph=system.get_binary_morph())


# visualization
    faces = [model["primary"]["t_object"]]
    color = [[colors[-1]] * len(model["primary"]["t_object"])]

    # visualization
    for i in model["primary"]["spots"]:
        faces.append(model["primary"]["spots"][i])
        # noinspection PyTypeChecker
        color.append([colors[i]] * len(model["primary"]["spots"][i]))

    faces.append(model["secondary"]["t_object"])
    color.append([colors[i + 1]] * len(model["secondary"]["t_object"]))

    Plt.plot_3d(vertices=None, faces=faces, normals_view=False, points_view=False, faces_view=True,
                point_color=None, point_size=5.0, face_color=color, edge_color="k",
                elev=45, azim=45, save=False)

# /over-contact system #################################################################################################

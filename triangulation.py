"""
###   In this document we show some examples of using 2d Delaunay triangulation in 3d space.                         ###
###   Mentioned method is created especially for triangulation of points describing W UMa                            ###
###   binary star system.                                                                                            ###
###   For W UMa surface triangulation we show also Poisson surface reconstruction                                    ###
###   implemented in The Computational Geometry Algorithms Library (CGAL)                                            ###
###   We also demonstrate using of classic 3d Delaunay triangulation in problem of Algol                             ###
###   binary stars.                                                                                                  ###
###                                                                                                                  ###
###   !!! For runing of any part of code, please change False to True in condition.                                  ###
"""
# imports ##############################################################################################################

import numpy as np
import Plot.Plot as Plt
from Star.StarCreator import StarCreator
import Tri.Triangulation as T
from Tri.Triangulation import delaunay3d as d3d
from Tri.Triangulation import Tri as Tri
import os

# /imports #############################################################################################################


# detached system ######################################################################################################

primary = {"mass": 3.0, "potential": 4.0, "phis": 5, "thetas": 10}
secondary = {"mass": 1.0, "potential": 4.0, "phis": 5, "thetas": 10}
orbit = {"period": 2.7, "eccentricity": 0.0, "inclination": np.pi, "argument_of_periastron": np.pi / 2.0}
system = StarCreator(primary=primary, secondary=secondary, orbit=orbit)
system.create_model()
binary_vertices = system.get_vertices()
binary_norms = system.get_norms()


# dot surface visualization
# -------------------------
if False:
    Plt.plot_3d(vertices=[binary_vertices["system"]], faces=None, normals_view=False, points_view=True,
                faces_view=False, point_color="r", point_size=10.0, face_color="c", edge_color="k", azim=45, elev=45,
                save=False, s_format="pdf", x_range=[-0.2, 1.2], y_range=[-0.5, 0.5], z_range=[-0.5, 0.5], dpi=1000)


# application of 3d Delaunay convex hull triangulation to primary component of star system
# ----------------------------------------------------------------------------------------
if False:
    primary_faces = d3d(vertices=binary_vertices["primary"])[1]
    # vizualisation
    # for plotting without vertices, set variable "points_view" in function "plot_3d" to False
    Plt.plot_3d(vertices=[binary_vertices["primary"]], faces=[primary_faces], normals_view=False, points_view=True,
                faces_view=True, point_color="r", point_size=5.0, face_color="c", edge_color="k", save=False,
                s_format="pdf", elev=45, azim=45, grid=True)


# application of 2d Delaunay triangulation in 3d space
# ----------------------------------------------------
if False:
    # INFO: it might take some times until operation will finish
    tri = Tri(vertices=binary_vertices["primary"], norms=binary_norms["primary"])
    tri.triangulate()
    hull = tri.hull()
    # for plotting without vertices, set variable "points_view" in function "plot_3d" to False
    Plt.plot_3d(vertices=[binary_vertices["primary"]], faces=[hull], normals_view=False, points_view=True,
                faces_view=True, point_color="r", point_size=5.0, face_color="c", edge_color="k")

# /detached system #####################################################################################################


# over-contact system ##################################################################################################
primary = {"mass": 3.0, "potential": 2.485, "phis": 5, "thetas": 10}
secondary = {"mass": 1.2, "potential": 2.485, "phis": 5, "thetas": 10}
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

# cgal triangulation
# -------------------------
if False:
    triangulation = T.cgal_triangulation(normals=binary_norms["system"], points=binary_vertices['system'],
                                         min_triangle_angle=np.radians([20.0])[0], max_triangle_size=1.0,
                                         surface_aproximation_error=0.13, to_average_spacing=5,
                                         input_path=os.path.dirname(os.path.realpath(__file__)))

    # for plotting without vertices, set variable "points_view" in function "plot_3d" to False
    Plt.plot_3d(vertices=[binary_vertices["system"]], faces=[triangulation[0]], normals_view=False, points_view=False,
                faces_view=True, point_color="r", point_size=5.0, face_color="c", edge_color="k",
                azim=45, elev=45, x_range=[-0.5, 1.4], y_range=[-0.5, 0.5], z_range=[-0.5, 0.5])

# triangulation of concave object using 2d Delaunay triangulation id 3d space
# ---------------------------------------------------------------------------
if False:
    # INFO: it might take some times until operation will finish
    tri = Tri(vertices=binary_vertices["system"], norms=binary_norms["system"])
    tri.triangulate()
    hull = tri.hull()
    # for plotting without vertices, set variable "points_view" in function "plot_3d" to False
    Plt.plot_3d(vertices=[binary_vertices["system"]], faces=[hull], normals_view=False, points_view=True,
                faces_view=True, point_color="r", point_size=5.0, face_color="c", edge_color="k",
                save=False, s_format="pdf", azim=45, elev=45, grid=True)

# /over-contact system #################################################################################################

if False:
    # INFO: it might take some times until operation will finish
    hull = d3d(vertices=binary_vertices["primary"])[1]
    # for plotting without vertices, set variable "points_view" in function "plot_3d" to False
    Plt.plot_3d(vertices=[binary_vertices["primary"]], faces=[hull], normals_view=False, points_view=True,
                faces_view=True, point_color="r", point_size=5.0, face_color="c", edge_color="k",
                save=False, s_format="pdf", azim=45, elev=45, grid=True)

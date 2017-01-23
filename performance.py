"""
###   Performence test                                                                                               ###
"""

import numpy as np
import os
from Star.StarCreator import StarCreator
from Tri.Triangulation import Tri as Tri
from Tri.Triangulation import delaunay3d as d3d
from time import time
import matplotlib.pyplot as plt
import pickle
plt.style.use('seaborn-white')
print(plt.style.available)

# detached system ######################################################################################################
if True:

    phis = list(range(3, 19, 1))
    try:
        xs = pickle.load(open("xs.pickle", "rb"))
        ys = pickle.load(open("ys.pickle", "rb"))
        phis = phis[len(xs):]
    except:
        xs, ys = [], [[], []]

    # for p in phis:
    #     print("phi steps:", p)
    #     primary = {"mass": 2.0, "potential": 4.0, "phis": p, "thetas": p * 2}
    #     secondary = {"mass": 1.0, "potential": 4.0, "phis": p, "thetas": p * 2}
    #     orbit = {"period": 2.0, "eccentricity": 0.0, "inclination": np.pi, "argument_of_periastron": np.pi / 2.0}
    #     system = StarCreator(primary=primary, secondary=secondary, orbit=orbit)
    #     system.create_model()
    #     binary_vertices = system.get_vertices()
    #     binary_norms = system.get_norms()
    #
    #     xs.append(len(binary_vertices["primary"]))
    #
    #     t = time()
    #     d3d(vertices=binary_vertices["primary"])
    #     ys[0].append(time() - t)
    #
    #     t = time()
    #     tri = Tri(vertices=binary_vertices["primary"], norms=binary_norms["primary"])
    #     tri.triangulate()
    #     ys[1].append(time() - t)
    #
    #     pickle.dump(xs, open("xs.pickle", "wb"))
    #     pickle.dump(ys, open("ys.pickle", "wb"))
    #
    #     del(primary, secondary, orbit, system)

    plt.scatter(xs, ys[0], c="b", s=100, label="3D Delaunay")
    plt.scatter(xs, ys[1], c="r", s=100, label="2D to 3D Delaunay")
    plt.legend(loc="upper left", markerscale=1.5, scatterpoints=1, fontsize=17)
    plt.xlabel('points [-]', fontsize=17)
    plt.ylabel('time [s]', fontsize=17)

    plt.tick_params(axis='both', which='major', labelsize=17)
    plt.tick_params(axis='both', which='minor', labelsize=17)

    plt.grid(True)
    plt.savefig("15_benchmark.pdf", dpi=900)
    plt.show()
# /detached system #####################################################################################################

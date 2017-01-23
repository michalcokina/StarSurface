import numpy as np
from os import environ
import matplotlib
import warnings

if 'DISPLAY' not in environ:
    matplotlib.use('pdf')
    import matplotlib.pyplot as plt
else:
    import matplotlib.pyplot as plt

warnings.simplefilter(action="ignore", category=FutureWarning)


def axis_equal_3d(ax):
    extents = np.array([getattr(ax, 'get_{}lim'.format(dim))() for dim in 'xyz'])
    sz = extents[:, 1] - extents[:, 0]
    centers = np.mean(extents, axis=1)
    maxsize = max(abs(sz))
    r = maxsize / 2
    for ctr, dim in list(zip(centers, 'xyz')):
        getattr(ax, 'set_{}lim'.format(dim))(ctr - r, ctr + r)


def figure_axis_range(arr):
    if type(arr) is type([]):
        arr = np.array(arr)
    maxim = - np.inf
    minim = np.inf
    for obj in arr:
        ma = np.amax(obj)
        mi = np.amin(obj)
        if ma > maxim:
            maxim = ma
        if mi < minim:
            minim = mi
        # if abs(mi) > maxim: maxim = abs(mi)
    return np.array([minim, maxim])


def plot_2d(point_color="r", point_marker="o", points=None, x_label="x", y_label="y", point_size=1., save=False,
            filename="plot", aspect="equal", line=False, grid=False, s_format="png"):
    zip_data = list(zip(*points))
    xs, ys = zip_data[0], zip_data[1]
    fig = plt.figure()
    ax = fig.add_subplot(111, aspect=aspect)
    ax.scatter([xs], [ys], color=str(point_color), marker=point_marker, s=point_size)

    if line: ax.plot(xs, ys, "-", c=str(point_color))

    ax.set_xlabel(str(x_label))
    ax.set_ylabel(str(y_label))

    ax.grid(grid)

    if not save:
        plt.show()
    else:
        if s_format == "png":
            plt.savefig(str(filename) + ".png")
        elif s_format == "pdf":
            plt.savefig(str(filename) + ".pdf")


def plot_3d(faces=None, normals=None, vertices=None, edge_color="k",
            face_color="w",  # matplotlib string or list
            normal_color="r", x_label="x", y_label="y", z_label="z", point_color="r", point_size=1.0, axis_off=False,
            faces_view=True, normals_view=True, points_view=True, azim=0, elev=0, face_alpha=1.0, save=False,
            filename="untitled", x_range=None, y_range=None, z_range=None, dpi=300, s_format="png", grid=False):
    import mpl_toolkits.mplot3d as a3
    warnings.simplefilter(action="ignore", category=FutureWarning)

    ax = a3.Axes3D(plt.figure())
    axis_range = 0.0

    # faces plotting
    # --------------
    idx_inner, idx_outer = 0, 0
    if faces_view:
        for obj in faces:
            idx_inner = 0
            for vtx in obj:
                # face color
                if type(face_color) == type(""):
                    fce_color = str(face_color)
                else:
                    fce_color = str(face_color[idx_outer][idx_inner])
                # /face color
                tri = a3.art3d.Poly3DCollection([vtx], alpha=face_alpha, edgecolors=str(edge_color))
                tri.set_facecolor(fce_color)
                ax.add_collection3d(tri)
                idx_inner += 1
            idx_outer += 1
        axis_range = figure_axis_range(arr=faces)

    # points ploting
    # --------------
    if points_view:

        for idx, obj in list(zip(range(len(vertices)), vertices)):
            if type(obj) is type([]):
                obj = np.array(obj)

            zip_points = obj.T  # zip(*obj)
            xs, ys, zs = zip_points[0], zip_points[1], zip_points[2]
            if type(point_color) == type(""):
                ax.scatter([xs], [ys], [zs], color=point_color, marker="o", s=point_size)
            else:
                ax.scatter([xs], [ys], [zs], color=[color for color in point_color[idx]], marker="o", s=point_size)
        axis_range = figure_axis_range(arr=vertices)

    # normals ploting
    # ---------------
    if normals_view:
        for obj_idx in range(0, len(normals)):
            for idx in range(0, len(normals[obj_idx])):
                # normal color
                norm_color = str(normal_color) if type(normal_color) == type("") else str(normal_color[obj_idx][idx])
                # /normal color

                segment = [vertices[obj_idx][idx], normals[obj_idx][idx]]
                line = a3.art3d.Line3DCollection([segment], alpha=1.0, zorder=1, colors=norm_color)
                ax.add_collection3d(line)
        axis_range = figure_axis_range(arr=vertices)

    # additional plot options
    # -----------------------

    if True:
        if empty(x_range) or empty(y_range) or empty(z_range):
            ax.set_xlim3d([-axis_range[1], axis_range[1]])
            ax.set_ylim3d([-axis_range[1], axis_range[1]])
            ax.set_zlim3d([-axis_range[1], axis_range[1]])
        else:
            ax.set_xlim3d(x_range)
            ax.set_ylim3d(y_range)
            ax.set_zlim3d(z_range)

    ax.set_xlabel(str(x_label), fontsize=15)
    ax.set_ylabel(str(y_label), fontsize=15)
    ax.set_zlabel(str(z_label), fontsize=15)

    ax.view_init(azim=azim, elev=elev)

    ax.set_autoscale_on(True)
    plt.gca().set_aspect('equal')
    axis_equal_3d(ax)

    ax.grid(True)
    # ax.set_axis_off()

    plt.tick_params(axis='both', which='major', labelsize=15)
    # plt.tick_params(axis='both', which='minor', labelsize=10)
    # ax.set_xticks()
    # ax.set_yticks()
    # ax.set_zticks()

    if axis_off:
        plt.axis('off')

    if not save:
        plt.show()
    else:
        if s_format == "png":
            plt.savefig(str(filename) + ".png", dpi=dpi)
        elif s_format == "pdf":
            plt.savefig(str(filename) + ".pdf", dpi=dpi, transparent=True)


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

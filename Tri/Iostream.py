import numpy as np
import os


# save point cloud with normals
def save_cgal_pcd_with_normals(filedir='tmp/', filename='3DModel.xyz', points=None, normals=None):
    filename, filedir = str(filename), str(filedir)

    f = open(filedir + filename, 'w')
    if os.path.isfile(filedir + filename):
        for idx in range(0, len(normals)):
            f.write(str(points[idx][0]) + " " + str(points[idx][1]) + " " + str(points[idx][2]) + " " + str(
                normals[idx][0]) + " " + str(normals[idx][1]) + " " + str(normals[idx][2]) + "\n")
    else:
        try:
            f.close()
        except:
            pass
        return False
    f.close()
    os.chmod(filedir + filename, 0o755)
    return True

# load cgal output mesh file
def load_cgal_3d_mesh(filename="out.off", filedir="tmp/"):
    path = str(filedir) + str(filename)

    if os.path.isfile(path):
        vertices, indices, faces = [], [], []  # real points, index of triangle points
        try:
            f = open(path, 'r')
            for _ in range(3):
                next(f)  # skip first 3 header lines
            for line in f:
                line = str(line).strip()
                if line is not "":
                    # if line start with number 3, line represent simplices
                    if line[:2] == "3 ":
                        line = line[3:]
                        num = line.split(" ")
                        indices.append([int(num[0]), int(num[1]), int(num[2])])
                    # otherwise there are vertices
                    else:
                        num = line.split(" ")
                        vertices.append([float(num[0]), float(num[1]), float(num[2])])
            f.close()

            # save faces to same shape as Delaunay triangulation using convex_hull
            for index in indices:
                faces.append([
                    vertices[index[0]],
                    vertices[index[1]],
                    vertices[index[2]]
                ])
            return np.array(faces), np.array(vertices), np.array(indices)
        # in case of any problem, return False
        except:
            return False
    else:
        return False



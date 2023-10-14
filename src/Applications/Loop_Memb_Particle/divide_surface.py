#from mpl_toolkits import mplot3d
#import matplotlib.pyplot as plt
import numpy as np
import math

# This file would output a dictionary of cubes:[triangles] that associates
# all triangles with a particular cube.


points = list()
connectivity = list()
nx = 5
ny = 5
nz = 5
ax = 0; ay=0; az = 0
bx = 0; by = 0; bz = 0;


def setup():
    '''
    Basic set up to read the points and connectivity from vtk file
    '''

    file = open("sphere_212.vtk", "r") #open("sphere_42_no_northpole.vtk", "r")
    #file = open("sphere_642.vtk", "r") #open("sphere_42_no_northpole.vtk", "r")
    #file = open("sphere_2562.vtk", "r") #open("sphere_42_no_northpole.vtk", "r")
    for i in range(0, 4):
        file.readline()

    num_points = int(file.readline().split()[1])  # Get the points data to variable - points
    for i in range(0, num_points):
        str_p = file.readline().split()
        float_p = [float(x) for x in str_p]
        points.append(float_p)
    file.readline()  # skip the empty line

    num_polys = int(file.readline().split()[1])  # Get the connectivity data to variable - connectivity
    for i in range(0, num_polys):
        str_poly = file.readline().split()[1:]
        int_poly = [int(x) for x in str_poly]
        connectivity.append(int_poly)


def gen_points():
    x = [p[0] for p in points]
    y = [p[1] for p in points]
    z = [p[2] for p in points]
    global ax, ay, az
    global bx, by, bz
    ep = 0.5;
    ax = (nx-2*ep) / (max(x) - min(x))
    bx = -ax*min(x) + ep
    ay = (ny-2*ep) / (max(y) - min(y))
    by = -ay * min(y) + ep
    az = (nz-2*ep) / (max(z) - min(z))
    bz = -az * min(z)+ep 
    print(ax)
    print(ay)
    print(az)
    print(bx)
    print(by)
    print(bz)    
    global scaled_x, scaled_y, scaled_z
    scaled_x = [i * ax + bx for i in x]
    scaled_y = [i * ay + by for i in y]
    scaled_z = [i * az + bz for i in z]

    # # Plot the scaled points
    # ax = plt.axes(projection='3d')
    # xticks = np.arange(math.floor(min(scaled_x)), math.ceil(max(scaled_x)) + 1, 1)
    # yticks = np.arange(math.floor(min(scaled_y)), math.ceil(max(scaled_y)) + 1, 1)
    # zticks = np.arange(math.floor(min(scaled_z)), math.ceil(max(scaled_z)) + 1, 1)
    # ax.set_xticks(xticks)
    # ax.set_yticks(yticks)
    # ax.set_zticks(zticks)
    # ax.scatter3D(scaled_x, scaled_y, scaled_z, c="blue")
    # for c in connectivity:
    #     plot_triangle_out_of_con(c, ax)
    # ax.grid()
    # plt.show()

    # Query line segment
    analyze_cubes()


def analyze_cubes():
    '''
    Returned cubes are defined by (floor_x, floor_y) as the lower left corner.
    Each line segment is represented by cubes that might be duplicated.
    Example: (2,7) refers to the cube with 4 vertices (x=2, x=3, y=7, y=8)
    '''
    face_dict = dict() # Output all cubes associated with a line segment.
    cube_dict = dict() # Output all line segments associated with a cube.

    for i in range(len(connectivity)):
        face_dict[i] = find_all_cubes(connectivity[i])

    all_cubes = list()
    for x in range(math.floor(min(scaled_x)), math.ceil(max(scaled_x))):
        for y in range(math.floor(min(scaled_y)), math.ceil(max(scaled_y))):
            for z in range(math.floor(min(scaled_z)), math.ceil(max(scaled_z))):
                all_cubes.append((x,y,z))

    for c in all_cubes:
        temp = list()
        for face in face_dict:
            if c in face_dict[face]:
                temp.append(face)
        if temp: # eliminate empty cubes
            cube_dict[c] = temp
    fid = open("dict_212.txt","w")
    fid.write(str(ax)+"\n")
    fid.write(str(ay)+"\n")
    fid.write(str(az)+"\n")
    fid.write(str(bx)+"\n")
    fid.write(str(by)+"\n")
    fid.write(str(bz)+"\n")
    for d in cube_dict:
        fid.write(str(d[0])+" "+str(d[1])+" "+str(d[2]))
        fid.write("\t")
        for ent in cube_dict[d]:
            fid.write(str(ent)+" ")
        fid.write("\n")
    fid.close()
    print("Keys are each faces, and values are the associated cubes.")
    print(face_dict)
    #test_cube_dict((2,-2,3), cube_dict) # Can change to random cubes.
    print("Keys are each cubes, and values are the faces around it.")
    print(cube_dict)

def find_all_cubes(one_connectivity):
    '''
    Given the information about a face, find all cubes that it lies on
    '''
    result = list()  # a list of tuples in form of (#, #, #)
    xs = [scaled_x[one_connectivity[0]], scaled_x[one_connectivity[1]], scaled_x[one_connectivity[2]]]
    ys = [scaled_y[one_connectivity[0]], scaled_y[one_connectivity[1]], scaled_y[one_connectivity[2]]]
    zs = [scaled_z[one_connectivity[0]], scaled_z[one_connectivity[1]], scaled_z[one_connectivity[2]]]
    lx = math.floor(min(xs))
    hx = math.ceil(max(xs))
    ly = math.floor(min(ys))
    hy = math.ceil(max(ys))
    lz = math.floor(min(zs))
    hz = math.ceil(max(zs))
    for i in range(lx, hx):
        for j in range(ly, hy):
            for k in range(lz, hz):
                temp = (i, j, k)
                result.append(temp)
    return result


def test_cubes(one_connectivity):
    '''
    Project the face onto x-y plane, x-z plane, and y-z plane to check if
    the associated cubes are correct.
    '''
    p1 = [scaled_x[one_connectivity[0]], scaled_y[one_connectivity[0]], scaled_z[one_connectivity[0]]]
    p2 = [scaled_x[one_connectivity[1]], scaled_y[one_connectivity[1]], scaled_z[one_connectivity[1]]]
    p3 = [scaled_x[one_connectivity[2]], scaled_y[one_connectivity[2]], scaled_z[one_connectivity[2]]]

    result = find_all_cubes(one_connectivity)
    print(result)
    print("x ranges from: {:d} to {:d}".format(min([r[0] for r in result]), max([r[0] for r in result])))
    print("y ranges from: {:d} to {:d}".format(min([r[1] for r in result]), max([r[1] for r in result])))
    print("z ranges from: {:d} to {:d}".format(min([r[2] for r in result]), max([r[2] for r in result])))

    #TODO get them into 1 figure
    #plt.subplot(3, 1, 1)
    plot_triangle(p1, p2, p3, 2)
    #plt.subplot(3, 1, 2)
    plot_triangle(p1, p2, p3, 1)
    #plt.subplot(3, 1, 3)
    plot_triangle(p1, p2, p3, 0)

    plt.show()

def test_cube_dict(rand_cube, cube_dict):
    # Test the associated face for a random cube, rand_cube is in
    # terms of (2, -2, 3).
    print(cube_dict[rand_cube])

def plot_triangle_out_of_con(c, ax):
    x_ls = [scaled_x[c[0]], scaled_x[c[1]], scaled_x[c[2]]]
    y_ls = [scaled_y[c[0]], scaled_y[c[1]], scaled_y[c[2]]]
    z_ls = [scaled_z[c[0]], scaled_z[c[1]], scaled_z[c[2]]]
    ax.plot(x_ls, y_ls, z_ls, color='grey')


def plot_triangle(p1, p2, p3, direction):
    # The direction is between {0,1,2}, and if the triangle is in
    # x-y plane, then the z-axis is zero, so "direction" = 2.
    x_ls = [p1[0], p2[0], p3[0]]
    y_ls = [p1[1], p2[1], p3[1]]
    z_ls = [p1[2], p1[2], p3[2]]
    if direction == 0:
        plt.plot(y_ls, z_ls, "blue")
        plt.plot([p1[1], p3[1]], [p1[2], p3[2]], "blue")
        plt.title("Projection onto y-z plane")
    if direction == 1:
        plt.plot(x_ls, z_ls, "blue")
        plt.plot([p1[0], p3[0]], [p1[2], p3[2]], "blue")
        plt.title("Projection onto x-z plane")
    if direction == 2:
        plt.plot(x_ls, y_ls, "blue")
        plt.plot([p1[0], p3[0]], [p1[1], p3[1]], "blue")
        plt.title("Projection onto x-y plane")
    plt.show()


if __name__ == '__main__':
    setup()
    gen_points()
    # find_all_cubes((3, 20, 19))
    # print(scaled_x[3], scaled_y[3], scaled_z[3])
    # print(scaled_x[20], scaled_y[20], scaled_z[20])
    # print(scaled_x[19], scaled_y[19], scaled_z[19])
    # test_cubes((3, 20, 19))

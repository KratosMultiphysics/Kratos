from scipy import spatial
import math
import numpy as np
import os
import sys

from terrain_data_container import data_container

from KratosMultiphysics import *


def mesh_in_concave_hull ( data ):

    ###	  P R E P A R A T I O N   ###
    # evaluation of imported data
    X = []; Y = []; Z = []      # lists are generated to use predefined operations

    for k, v in data.vertex.items():
        if (k[2] == 0.0):	# if (Z == 0.0)
            continue        # fill just coord_3D and not coord_2D
        X.append(k[0])
        Y.append(k[1])
        Z.append(k[2])

    # calculate minimum and maximum value of X, Y and Z coordinate
    x_min = min(X); x_max = max(X)
    y_min = min(Y); y_max = max(Y)
    z_min = min(Z); z_max = max(Z)

    vertex_id_list = data.id_vertex		# filling the id list

    import triangle

    # fill x and y lists with 2D coordinate
    x = []; y = []


    # fill pts list with the couple x, y
    pts = np.vstack((x, y)).T
    A = dict(vertices = pts)
    triangle_dict = triangle.triangulate(A)		# triangle_dict {vertices: [...], triangles: [...], vertex_markers: [...]}












def mesh_in_circle ( data ):

    ###	  P R E P A R A T I O N   ###
    # evaluation of imported data
    X = []; Y = []; Z = []      # lists are generated to use predefined operations

    for k, v in data.vertex.items():
        if (k[2] == 0.0):	# if (Z == 0.0)
            continue        # fill just coord_3D and not coord_2D
        X.append(k[0])
        Y.append(k[1])
        Z.append(k[2])

    # calculate minimum and maximum value of X, Y and Z coordinate
    x_min = min(X); x_max = max(X)
    y_min = min(Y); y_max = max(Y)
    z_min = min(Z); z_max = max(Z)

    # calculate the centre of circle
    x_center = (x_min + x_max)/2.0
    y_center = (y_min + y_max)/2.0

    # calculate the radius considering the smaller side of the rectangle
    r_boundary = min((x_max-x_min), (y_max-y_min))/2.0	    # radius of the domain
    r_ground = r_boundary * 2.0/3.0							# 2/3 of the radius of the domain
    r_buildings = r_boundary * 1.0/3.0						# 1/3 of the radius of the domain

    # buffering zone
    delta = r_boundary / 1000.0			# evaluate this value!!!
    print("PREPARATION DONE!")

    ###	  S M O O T H I N G   ###
    # different smoothing levels in defined circular regions

    vertex_id_list = data.id_vertex		# filling the id list

    # calculate id nodes that are inside the r_ground
    points_no_smooth = np.array(data.coord_2D)		# points that do not require smooth
    tree = spatial.KDTree(points_no_smooth)
    idx_inside = tree.query_ball_point([x_center, y_center], r = r_ground)		# node ids that are inside 2/3 of the circle

    # here there are ids that must to be smooth
    idx_to_smooth = list( set(vertex_id_list) - set(idx_inside) )

    del_id = []						    # list where there are the ids to be deleted because are outside of r_boundary
    for index in idx_to_smooth:
        x_curr = data.coord_3D[index][0]	    # current X coordinate
        y_curr = data.coord_3D[index][1]	    # current Y coordinate
        dist = math.sqrt(((x_center-x_curr)**2)+((y_center-y_curr)**2))		# calculate the distance between current node and center of circle

        if (dist > r_boundary-delta):
            del_id.append(index)
            continue 							# here beta = 0 because the node is outside the circle boundary
        elif (dist > r_ground):
            Z_beta = (-(dist-r_ground) / (r_boundary-r_ground)) + 1.0		# calculate beta
        else:
            # dist <= r_ground
            Z_beta = 1.0

        data.coord_3D[index][2] = (data.coord_3D[index][2] - z_min) * Z_beta + z_min

    # delete the nodes from coord_3D that are outside the r_boundary
    for i in del_id:
        del data.coord_3D[i]

    # create an array with 360 degree (from 0 to 360) -> 2*pi / 360 = pi /180
    theta = np.arange(0, 2*np.pi, np.pi/90.0)			# old value -->  theta = np.arange(0, 2*np.pi, np.pi/180)

    # define the nodes on circle boundary
    x_circle = []
    y_circle = []
    for th in theta:
        x_circle.append(r_boundary * np.cos(th) + x_center)	# move the node along X to consider that the circle have not the centre in (0.0, 0.0)
        y_circle.append(r_boundary * np.sin(th) + y_center)	# move the node along Y to consider that the circle have not the centre in (0.0, 0.0)
    print("SMOOTHING DONE!")

    ###	  T R I A N G L E	###
    # compute the terrain mesh (surface mesh with only terrain)
    import triangle

    # fill x and y lists with 2D coordinate
    x = []; y = []
    for _, coords in data.coord_3D.items():
        x.append(coords[0])
        y.append(coords[1])

    # add circumference's nodes into x and y lists
    for i in range(len(x_circle)):
        x.append(x_circle[i])
        y.append(y_circle[i])

    # fill pts list with the couple x, y
    pts = np.vstack((x, y)).T
    A = dict(vertices = pts)
    triangle_dict = triangle.triangulate(A)		# triangle_dict {vertices: [...], triangles: [...], vertex_markers: [...]}

    print( triangle_dict )

    all_points = []		# vector with all coordinates of nodes
    all_facets = []		# vector with all node ids of faces

    # fill all_points vector
    for _, coords in data.coord_3D.items():
        all_points.append((coords[0], coords[1], coords[2]))

    inner_id_terrain = len(all_points)		# id up to this moment
    list_id_circle_bot = []					# array with only ids of the bottom circle
    for i in range(len(x_circle)):
        all_points.append((x_circle[i], y_circle[i], z_min))
        list_id_circle_bot.append(i+inner_id_terrain)

    for f in triangle_dict["triangles"]:
        all_facets.append([f[0], f[1], f[2]])

    number_node_terrain = len(all_points)

    volume_height = z_max + 200
    # fill all_points vector to "topper"
    for _, coords in data.coord_3D.items():
        all_points.append((coords[0], coords[1], volume_height))

    inner_id_topper	= len(all_points)		# id up to this moment
    list_id_circle_top = []					# array with only ids of the top circle
    for i in range(len(x_circle)):
        all_points.append((x_circle[i], y_circle[i], volume_height))
        list_id_circle_top.append(i+inner_id_topper)

    for f in triangle_dict["triangles"]:
        all_facets.append([f[0]+number_node_terrain, f[1]+number_node_terrain, f[2]+number_node_terrain])		# list of array

    # fill the information about the lateral surface
    for i in range(len(list_id_circle_bot)):
        if (i == len(list_id_circle_bot)-1):
            all_facets.append([list_id_circle_bot[i], list_id_circle_bot[0], list_id_circle_top[0], list_id_circle_top[i]])
            break
        all_facets.append([list_id_circle_bot[i], list_id_circle_bot[i+1], list_id_circle_top[i+1], list_id_circle_top[i]])

    print("MESH WITH TRIANGLE DONE!")

    return ( all_points, all_facets )



def mesh_in_rectangel ( data ):

    ###	  P R E P A R A T I O N   ###
    # evaluation of imported data
    X = []; Y = []; Z = []      # lists are generated to use predefined operations

    for k, v in data.vertex.items():
        if (k[2] == 0.0):	# if (Z == 0.0)
            continue        # fill just coord_3D and not coord_2D
        X.append(k[0])
        Y.append(k[1])
        Z.append(k[2])

    # calculate minimum and maximum value of X, Y and Z coordinate
    x_min = min(X); x_max = max(X)
    y_min = min(Y); y_max = max(Y)
    z_min = min(Z); z_max = max(Z)

    # calculate the centre of circle
    x_center = (x_min + x_max)/2.0
    y_center = (y_min + y_max)/2.0

    print("PREPARATION DONE!")

    ###	  S M O O T H I N G   ###
    # different smoothing levels in defined circular regions

    vertex_id_list = data.id_vertex		# filling the id list

    # calculate id nodes that are inside the r_ground
    points_no_smooth = np.array(data.coord_2D)		# points that do not require smooth
    tree = spatial.KDTree(points_no_smooth)
    idx_inside = tree.query_ball_point([x_center, y_center], r = r_ground)		# node ids that are inside 2/3 of the circle

    # here there are ids that must to be smooth
    idx_to_smooth = list( set(vertex_id_list) - set(idx_inside) )

    del_id = []						    # list where there are the ids to be deleted because are outside of r_boundary
    for index in idx_to_smooth:
        x_curr = data.coord_3D[index][0]	    # current X coordinate
        y_curr = data.coord_3D[index][1]	    # current Y coordinate
        dist = math.sqrt(((x_center-x_curr)**2)+((y_center-y_curr)**2))		# calculate the distance between current node and center of circle

        if (dist > r_boundary-delta):
            del_id.append(index)
            continue 							# here beta = 0 because the node is outside the circle boundary
        elif (dist > r_ground):
            Z_beta = (-(dist-r_ground) / (r_boundary-r_ground)) + 1.0		# calculate beta
        else:
            # dist <= r_ground
            Z_beta = 1.0

        data.coord_3D[index][2] = (data.coord_3D[index][2] - z_min) * Z_beta + z_min

    # delete the nodes from coord_3D that are outside the r_boundary
    for i in del_id:
        del data.coord_3D[i]

    # create an array with 360 degree (from 0 to 360) -> 2*pi / 360 = pi /180
    theta = np.arange(0, 2*np.pi, np.pi/90.0)			# old value -->  theta = np.arange(0, 2*np.pi, np.pi/180)

    # define the nodes on circle boundary
    x_circle = []
    y_circle = []
    for th in theta:
        x_circle.append(r_boundary * np.cos(th) + x_center)	# move the node along X to consider that the circle have not the centre in (0.0, 0.0)
        y_circle.append(r_boundary * np.sin(th) + y_center)	# move the node along Y to consider that the circle have not the centre in (0.0, 0.0)
    print("SMOOTHING DONE!")

    ###	  T R I A N G L E	###
    # compute the terrain mesh (surface mesh with only terrain)
    import triangle

    # fill x and y lists with 2D coordinate
    x = []; y = []
    for _, coords in data.coord_3D.items():
        x.append(coords[0])
        y.append(coords[1])

    # add circumference's nodes into x and y lists
    for i in range(len(x_circle)):
        x.append(x_circle[i])
        y.append(y_circle[i])

    # fill pts list with the couple x, y
    pts = np.vstack((x, y)).T
    A = dict(vertices = pts)
    triangle_dict = triangle.triangulate(A)		# triangle_dict {vertices: [...], triangles: [...], vertex_markers: [...]}

    # print( triangle_dict )

    all_points = []		# vector with all coordinates of nodes
    all_facets = []		# vector with all node ids of faces

    # fill all_points vector
    for _, coords in data.coord_3D.items():
        all_points.append((coords[0], coords[1], coords[2]))

    inner_id_terrain = len(all_points)		# id up to this moment
    list_id_circle_bot = []					# array with only ids of the bottom circle
    for i in range(len(x_circle)):
        all_points.append((x_circle[i], y_circle[i], z_min))
        list_id_circle_bot.append(i+inner_id_terrain)

    for f in triangle_dict["triangles"]:
        all_facets.append([f[0], f[1], f[2]])

    number_node_terrain = len(all_points)

    volume_height = z_max + 200
    # fill all_points vector to "topper"
    for _, coords in data.coord_3D.items():
        all_points.append((coords[0], coords[1], volume_height))

    inner_id_topper	= len(all_points)		# id up to this moment
    list_id_circle_top = []					# array with only ids of the top circle
    for i in range(len(x_circle)):
        all_points.append((x_circle[i], y_circle[i], volume_height))
        list_id_circle_top.append(i+inner_id_topper)

    for f in triangle_dict["triangles"]:
        all_facets.append([f[0]+number_node_terrain, f[1]+number_node_terrain, f[2]+number_node_terrain])		# list of array

    # fill the information about the lateral surface
    for i in range(len(list_id_circle_bot)):
        if (i == len(list_id_circle_bot)-1):
            all_facets.append([list_id_circle_bot[i], list_id_circle_bot[0], list_id_circle_top[0], list_id_circle_top[i]])
            break
        all_facets.append([list_id_circle_bot[i], list_id_circle_bot[i+1], list_id_circle_top[i+1], list_id_circle_top[i]])

    print("MESH WITH TRIANGLE DONE!")

    return ( all_points, all_facets )
'''
This demo demonstrates how to perform closest point queries using FCPW.
The full FCPW API can be viewed using the following commands in the Python console:
>>> import fcpw
>>> help(fcpw)
'''

import numpy as np
import polyscope as ps
import fcpw
import argparse
from pathlib import Path
from typing import Callable, Optional

def load_obj(obj_file_path):
    positions = []
    indices = []

    with open(obj_file_path, 'r') as file:
        for line in file:
            if line.startswith('v '):
                position = list(map(float, line.strip().split()[1:]))
                positions.append(np.array(position, dtype=np.float32, order='C'))

            elif line.startswith('f '):
                index = [int(idx.split('/')[0]) - 1 for idx in line.strip().split()[1:]]
                indices.append(np.array(index, dtype=np.int32, order='C'))

    return positions, indices

def load_fcpw_scene(positions, indices, build_vectorized_cpu_bvh):
    # load positions and indices
    scene = fcpw.scene_3D()
    scene.set_object_count(1)
    scene.set_object_vertices(positions, 0)
    scene.set_object_triangles(indices, 0)

    # build scene on CPU
    aggregate_type = fcpw.aggregate_type.bvh_surface_area
    print_stats = False
    reduce_memory_footprint = False
    scene.build(aggregate_type, build_vectorized_cpu_bvh,
                print_stats, reduce_memory_footprint)

    return scene

def perform_closest_point_queries(scene, query_points):
    # initialize bounding spheres
    bounding_spheres = fcpw.bounding_sphere_3D_list()
    for q in query_points:
        bounding_spheres.append(fcpw.bounding_sphere_3D(q, np.inf))

    # perform cpqs
    interactions = fcpw.interaction_3D_list()
    scene.find_closest_points(bounding_spheres, interactions)

    # extract closest points
    closest_points = np.array([i.p for i in interactions])

    return closest_points

def perform_gpu_closest_point_queries(gpu_scene, query_points):
    # initialize bounding spheres
    bounding_spheres = fcpw.gpu_bounding_sphere_list()
    for q in query_points:
        gpu_query_point = fcpw.float_3D(q[0], q[1], q[2])
        bounding_spheres.append(fcpw.gpu_bounding_sphere(gpu_query_point, np.inf))

    # perform cpqs on GPU
    interactions = fcpw.gpu_interaction_list()
    gpu_scene.find_closest_points(bounding_spheres, interactions)

    # extract closest points
    closest_points = np.array([np.array([i.p.x, i.p.y, i.p.z], dtype=np.float32, order='C') for i in interactions])

    return closest_points

def gui_callback(scene, query_points, use_gpu):
    # animate query points
    for q in query_points:
        q[0] += 0.001 * np.sin(10.0 * q[1])
        q[1] += 0.001 * np.cos(10.0 * q[0])

    # perform closest point queries
    closest_points = None
    if use_gpu:
        closest_points = perform_gpu_closest_point_queries(scene, query_points)
    else:
        closest_points = perform_closest_point_queries(scene, query_points)

    # plot results
    ps.register_point_cloud("query points", query_points)
    ps.register_point_cloud("closest points", closest_points)
    edge_positions = np.concatenate([query_points, closest_points], axis=0)
    edge_indices = np.array([[i, i + len(query_points)] for i in range(len(query_points))])
    network = ps.register_curve_network("edges", edge_positions, edge_indices)
    network.set_radius(0.005, relative=False)

def visualize(scene, positions, indices, query_points, use_gpu):
    # initialize polyscope
    ps.init()
    ps.set_ground_plane_mode("none")

    # register mesh and callback
    ps.register_surface_mesh("mesh", np.array(positions), np.array(indices))
    gui_callback_no_args = lambda: gui_callback(scene, query_points, use_gpu)
    ps.set_user_callback(gui_callback_no_args)

    # give control to polyscope gui
    ps.show()

def main():
    # parse arguments
    parser = argparse.ArgumentParser(description="fcpw demo")
    parser.add_argument("--use_gpu", action="store_true", help="use GPU")
    args = parser.parse_args()

    # load obj file
    positions, indices = load_obj("dragon.obj")

    # generate random query points for closest point queries
    num_query_points = 100
    box_min = np.min(positions, axis=0)
    box_max = np.max(positions, axis=0)
    query_points = np.random.uniform(box_min, box_max, (num_query_points, 3)).astype(np.float32)

    if args.use_gpu:
        # load fcpw scene
        scene = load_fcpw_scene(positions, indices, False) # NOTE: must build non-vectorized CPU BVH

        # transfer scene to GPU
        fcpw_directory_path = str(Path.cwd().parent)
        print_stats = False
        gpu_scene = fcpw.gpu_scene_3D(fcpw_directory_path, print_stats)
        gpu_scene.transfer_to_gpu(scene)

        # visualize scene
        visualize(gpu_scene, positions, indices, query_points, True)

    else:
        # load fcpw scene
        scene = load_fcpw_scene(positions, indices, True)

        # visualize scene
        visualize(scene, positions, indices, query_points, False)

if __name__ == "__main__":
    main()
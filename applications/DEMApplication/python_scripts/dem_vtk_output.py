import os
from pyevtk import hl
from pyevtk import vtk
import numpy as np
import weakref
import shutil
from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *

class VtkOutput():
    def __init__(self, main_path, problem_name, spheres_model_part, rigid_face_model_part):
        self.problem_name = problem_name

        self.particles_X = np.empty(0)
        self.particles_Y = np.empty(0)
        self.particles_Z = np.empty(0)
        self.particles_R = np.empty(0)
        self.particles_material = np.empty(0)
        self.velocities_X = np.empty(0)
        self.velocities_Y = np.empty(0)
        self.velocities_Z = np.empty(0)
        self.spheres_model_part = weakref.proxy(spheres_model_part)

        self.walls_X = np.empty(0)
        self.walls_Y = np.empty(0)
        self.walls_Z = np.empty(0)
        self.walls_connectivity = np.empty(0)
        self.walls_offsets = np.empty(0)
        self.walls_cell_types = np.empty(0)
        self.rigid_face_model_part = weakref.proxy(rigid_face_model_part)

        self.vtk_post_path_directory = os.path.join(main_path, problem_name + "_Post_VTK_Files")
        shutil.rmtree(self.vtk_post_path_directory, ignore_errors=True)
        os.makedirs(str(self.vtk_post_path_directory))

        self.counter = 0

    def ConvertParticlesToNumpyArrays(self):
        number_of_nodes = self.spheres_model_part.NumberOfNodes(0)

        self.particles_X = np.empty(number_of_nodes)
        self.particles_Y = np.empty(number_of_nodes)
        self.particles_Z = np.empty(number_of_nodes)
        self.particles_R = np.empty(number_of_nodes)
        self.particles_material = np.empty(number_of_nodes)
        self.velocities_X = np.empty(number_of_nodes)
        self.velocities_Y = np.empty(number_of_nodes)
        self.velocities_Z = np.empty(number_of_nodes)

        i = 0
        for node in self.spheres_model_part.Nodes:
            self.particles_X[i] = node.X
            self.particles_Y[i] = node.Y
            self.particles_Z[i] = node.Z
            self.particles_R[i] = node.GetSolutionStepValue(RADIUS)
            self.particles_material[i] = node.GetSolutionStepValue(PARTICLE_MATERIAL)
            self.velocities_X[i] = node.GetSolutionStepValue(VELOCITY_X)
            self.velocities_Y[i] = node.GetSolutionStepValue(VELOCITY_Y)
            self.velocities_Z[i] = node.GetSolutionStepValue(VELOCITY_Z)
            i += 1


    def ConvertWallsToNumpyArrays(self):
        number_of_nodes = self.rigid_face_model_part.NumberOfNodes(0)
        self.walls_X = np.empty(number_of_nodes)
        self.walls_Y = np.empty(number_of_nodes)
        self.walls_Z = np.empty(number_of_nodes)

        position_of_each_id = {}
        i = 0
        for node in self.rigid_face_model_part.Nodes:
            self.walls_X[i] = node.X
            self.walls_Y[i] = node.Y
            self.walls_Z[i] = node.Z
            position_of_each_id[node.Id] = i
            i += 1

        number_of_conditions = self.rigid_face_model_part.NumberOfConditions(0)
        self.walls_connectivity = np.empty(number_of_conditions * 3)
        self.walls_offsets = np.empty(number_of_conditions)
        self.walls_cell_types = np.empty(number_of_conditions)
        i = 0
        j = 0
        for cond in self.rigid_face_model_part.Conditions:
            list_of_vertices = cond.GetNodes()
            number_of_vertices = len(list_of_vertices)
            for k in range(number_of_vertices):
                index_in_array = position_of_each_id[list_of_vertices[k].Id]
                self.walls_connectivity[i] = index_in_array
                self.walls_offsets[j] = i+1
                i += 1

            self.walls_cell_types[j] = vtk.VtkTriangle.tid
            j += 1


    def WriteResults(self, time):

        self.ConvertParticlesToNumpyArrays()
        particles_filename = self.problem_name + "_DEM_" + str(self.counter)
        path = os.path.join(self.vtk_post_path_directory, particles_filename)
        hl.pointsToVTK(path, self.particles_X, self.particles_Y, self.particles_Z, {'radius':self.particles_R, 'material':self.particles_material, 'velocity':(self.velocities_X, self.velocities_Y, self.velocities_Z)})


        self.ConvertWallsToNumpyArrays()
        walls_filename = self.problem_name + "_Walls_" + str(self.counter)
        path = os.path.join(self.vtk_post_path_directory, walls_filename)
        hl.unstructuredGridToVTK(path, self.walls_X, self.walls_Y, self.walls_Z, self.walls_connectivity, self.walls_offsets, self.walls_cell_types, cellData=None, pointData=None)

        self.counter += 1


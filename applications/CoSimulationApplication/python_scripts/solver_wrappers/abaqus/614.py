import os
from os.path import join
import subprocess
import time
import numpy as np
import copy


import KratosMultiphysics as KM
from KratosMultiphysics.CoSimulationApplication.co_simulation_component import CoSimulationComponent
from KratosMultiphysics.CoSimulationApplication.co_simulation_interface import CoSimulationInterface
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools
cs_data_structure = cs_tools.cs_data_structure


def Create(parameters):
    return SolverWrapperAbaqus614(parameters)

class SolverWrapperAbaqus614(CoSimulationComponent):
    def __init__(self, parameters):
        super().__init__()

        #settings
        """
               settings of solver_wrappers.abaqus.614:

                   working_directory       absolute path to working directory
                                           or relative path w.r.t current directory
        """

        self.settings = parameters['settings']
        self.dir_csm = join(os.getcwd(), self.settings['working_directory'].GetString())  # *** alternative for getcwd?
        path_src = os.path.realpath(os.path.dirname(__file__))

        self.remove_all_messages()

        self.cores = self.settings['cores'].GetInt() #  number of cpus Abaqus has to use
        self.dimensions = self.settings['dimensions'].GetInt()
        self.array_size = self.settings["arraysize"].GetInt()
        self.surfaces = self.settings["surfaces"].GetInt()

        #prepare Abaqus USR
        path_usr = join(path_src, "USR.f")
        with open(path_usr)as infile:
            for line in infile:
                line = line.replace("|dimension|", self.dimensions)
                line = line.replace("|arraySize|", self.array_size)
                line = line.replace("|surfaces|", self.surfaces)
                line = line.replace("|cpus|", self.cores)
                line = line.replace("|PWD|", os.pardir(self.dir_csm))
                line = line.replace("|csmdir|", self.settings["working directory"])



        #compile Fortran and C++ codes
        path_libusr = join(self.dir_csm, "libusr")
        os.system("rm -r " + path_libusr)
        os.system("mdkir " + path_libusr)
        cmd = "abaqus make library=" + path_usr + " directory=" + path_libusr + " >> AbaqusSolver.log 2>&1"
        print(cmd)
        self.abaqus_process = subprocess.Popen(cmd, executable='/bin/bash', shell=True, cwd=self.dir_csm)

        # TODO:
        #   Read settings
        #   Prepare Abaqus usr (if necessary), probably some keyword substitutions
        #   Abaqus needs to print a file with the location of its load points (usr_init.f ?)
        #   Create Model with ModelParts
        #   Add variables to ModelParts
        #   Add Nodes to input ModelParts and output ModelParts, based on file written by Abaqus


    def Initialize(self):
        # super().Initialize()
        print('\nInitialize')

    def InitializeSolutionStep(self):
        # super().InitializeSolutionStep()
        #
        # self.iteration = 0
        # self.timestep += 1
        # print(f'\tTimestep {self.timestep}')
        #
        # self.send_message('next')
        # self.wait_message('next_ready')
        print('\nInitializeSolutionStep')

    def SolveSolutionStep(self, interface_input):
        print('\nSolveSolutionStep')
        return 0
    def FinalizeSolutionStep(self):
        # super().FinalizeSolutionStep()
        #
        # if not self.timestep % self.settings['save_iterations'].GetInt():
        #     self.send_message('save')
        #     self.wait_message('save_ready')
        print('\nFinalizeSolutionStep')

    def Finalize(self):
        # super().Finalize()
        # self.send_message('stop')
        # self.wait_message('stop_ready')
        # self.fluent_process.wait()
        print('Finalize')

    def GetInterfaceInput(self):
        return self.interface_input.deepcopy()

    def SetInterfaceInput(self):
        Exception("This solver interface provides no mapping.")

    def GetInterfaceOutput(self):
        return self.interface_output.deepcopy()

    def SetInterfaceOutput(self):
        Exception("This solver interface provides no mapping.")

    # def get_unique_face_ids(self, data):
    #     """
    #     Construct unique face IDs based on the face's node IDs.
    #
    #     Parameter data contains a 2D ndarray of node IDs.
    #     Each row corresponds to the unique node IDs corresponding
    #     to a certain face, supplemented with -1-values.
    #     The row is sorted, the -1-values are removed, and then a
    #     string is made by adding the unique node IDs together.
    #         e.g. for a row [5, 9, 7, -1, -1]
    #              the face ID is "5-7-9"
    #     """
    #     data = data.astype(int)
    #     ids = np.zeros(data.shape[0], dtype='U256')  # array is flattened
    #     for j in range(ids.size):
    #         tmp = np.unique(data[j, :])
    #         if tmp[0] == -1:
    #             tmp = tmp[1:]
    #         ids[j] = '-'.join(tuple(tmp.astype(str)))
    #     return ids
    #
    # def set_node_coordinates_test(self, f):
    #     for key in self.settings['interface_input'].keys():
    #         for node in self.model[key].Nodes:
    #             node.Y += (1 - np.cos(2 * np.pi * node.X)) * 0.5 * f
    #
    # def write_node_positions(self):
    #     for key in self.settings['interface_input'].keys():
    #         mp = self.model[key]
    #         tmp = f'nodes_update_timestep{self.timestep}_thread{mp.thread_id}.dat'
    #         file_name = join(self.dir_cfd, tmp)
    #         with open(file_name, 'w') as file:
    #             file.write(f'{mp.NumberOfNodes()}\n')
    #             for node in mp.Nodes:
    #                 if self.dimensions == 2:
    #                     file.write(f'{node.X:27.17e} {node.Y:27.17e} {node.Id:>27}\n')
    #                 else:
    #                     file.write(f'{node.X:27.17e} {node.Y:27.17e} {node.Z:27.17e} {node.Id:>27}\n')

    def send_message(self, message):
        file = join(self.dir_csm, message + ".msg")
        open(file, 'w').close()
        return

    def wait_message(self, message):
        file = join(self.dir_csm, message + ".msg")
        while not os.path.isfile(file):
            time.sleep(0.01)
        os.remove(file)
        return

    def check_message(self, message):
        file = join(self.dir_csm, message + ".msg")
        if os.path.isfile(file):
            os.remove(file)
            return True
        return False

    def remove_all_messages(self):
        for file_name in os.listdir(self.dir_csm):
            if file_name.endswith('.msg'):
                file = join(self.dir_csm, file_name)
                os.remove(file)

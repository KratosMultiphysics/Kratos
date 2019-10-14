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
    return SolverWrapperFluent2019R1(parameters)


class SolverWrapperFluent2019R1(CoSimulationComponent):
    def __init__(self, parameters):
        super().__init__()

        # settings
        self.settings = parameters['settings']
        self.dir_cfd = os.path.join(os.getcwd(), self.settings['working_directory'].GetString())  # *** alternative for getcwd?
        path_src = os.path.realpath(os.path.dirname(__file__))

        self.cores = self.settings['cores'].GetInt()
        self.case_file = self.settings['case_file'].GetString()
        self.mnpf = self.settings['max_nodes_per_face'].GetInt()
        self.dimensions = self.settings['dimensions'].GetInt()

        self.thread_names = [_.GetString() for _ in self.settings['thread_names'].list()]
        self.n_threads = len(self.thread_names)
        self.thread_ids = [None] * self.n_threads

        # prepare Fluent input journal
        journal = '2019R1.jou'
        thread_names_str = ''
        for key in self.thread_names:
            thread_names_str += ' "' + key + '"'
        if self.settings['hybrid_initialization'].GetBool():
            hybrid_initialization = '#t'
        else:
            hybrid_initialization = '#f'
        with open(os.path.join(path_src, journal), 'r') as infile:
            with open(os.path.join(self.dir_cfd, journal), 'w') as outfile:
                for line in infile:
                    line = line.replace('|case|', os.path.join(self.dir_cfd, self.case_file))  #*** change back to capitals, that's clearer...
                    line = line.replace('|thread_names|', thread_names_str)
                    line = line.replace('|hybrid_initialization|', hybrid_initialization)
                    line = line.replace('|iterations|', str(50))
                    outfile.write(line)

        # prepare Fluent UDF
        udf = '2019R1.c'
        with open(os.path.join(path_src, udf), 'r') as infile:
            with open(os.path.join(self.dir_cfd, udf), 'w') as outfile:
                for line in infile:
                    line = line.replace('|max_nodes_per_face|', str(self.mnpf))
                    outfile.write(line)

        # start Fluent with journal
        fluent_gui = self.settings['fluent_gui'].GetBool()
        gui = ''
        if not fluent_gui:
            gui = ' -gu'
        # subprocess.Popen(f'fluent 2ddp{gui} -t{self.cores} -i {journal}',  # *** ON/OFF
        #                      shell=True, executable='/bin/bash', cwd=self.dir_cfd)  # *** ON/OFF

        # get surface thread ID's from report.sum and write them to bcs.txt
        # self.wait_message('surface_info_exported')  # *** ON/OFF
        report = os.path.join(self.dir_cfd, 'report.sum')
        check = 0
        info = []
        with open(report, 'r') as file:
            for line in file:
                if check == 3 and line.islower():
                    name, id, _ = line.strip().split()
                    if name in self.thread_names:
                        info.append(' '.join((name, id)))
                        self.thread_ids[self.thread_names.index(name)] = id
                if check == 3 and not line.islower():
                    break
                if check == 2:  # skip 1 line
                    check = 3
                if 'name' in line and check == 1:
                    check = 2
                if 'Boundary Conditions' in line:
                    check = 1

        with open(os.path.join(self.dir_cfd, 'bcs.txt'), 'w') as file:
            file.write(str(len(info)) + '\n')
            for line in info:
                file.write(line + '\n')
        self.send_message('thread_ids_written_to_file')

        # import node and face information
        # self.wait_message('nodes_and_faces_stored')  # *** ON/OFF

        # create Model
        self.model = cs_data_structure.Model()

        # create ModelParts
        for key, value in (self.settings['interface_input'].items() +
                           self.settings['interface_output'].items()):
            # add ModelPart to Model
            self.model.CreateModelPart(key)
            mp = self.model[key]

            # add historical variables to ModelPart
            for var_name in value.list():
                var = vars(KM)[var_name.GetString()]
                mp.AddNodalSolutionStepVariable(var)

            # add information to ModelPart
            for i in range(self.n_threads):
                if self.thread_names[i] in key:
                    mp.thread_name = self.thread_names[i]
                    mp.thread_id = self.thread_ids[i]
            if 'thread_id' not in dir(mp):
                raise AttributeError('could not find thread name corresponding to key')

        # add Nodes to input ModelParts (nodes)
        for key in self.settings['interface_input'].keys():
            mp = self.model[key]

            # read in datafile
            tmp = 'nodes_thread' + str(mp.thread_id) + '.dat'
            file_name = os.path.join(self.dir_cfd, tmp)
            data = np.loadtxt(file_name, skiprows=1)
            if data.shape[1] != self.dimensions + 1:
                raise ValueError(f'given dimension does not match coordinates')

            # get node coordinates and ids
            coords_tmp = np.zeros((data.shape[0], 3)) * 0.
            coords_tmp[:, :self.dimensions] = data[:, :-1]  # add column z if 2D
            ids_tmp = data[:, -1].astype(int).astype(str)  # array is flattened

            # sort and remove doubles
            args = np.unique(ids_tmp, return_index=True)[1].tolist()
            coords_tmp = coords_tmp[args, :]
            ids_tmp = ids_tmp[args]

            # create Nodes
            for i in range(ids_tmp.size):
                mp.CreateNewNode(ids_tmp[i],
                    coords_tmp[i, 0], coords_tmp[i, 1], coords_tmp[i, 2])

        # add Nodes to output ModelParts (faces)
        for key in self.settings['interface_output'].keys():
            mp = self.model[key]

            # read in datafile
            tmp = 'faces_thread' + str(mp.thread_id) + '.dat'
            file_name = os.path.join(self.dir_cfd, tmp)
            data = np.loadtxt(file_name, skiprows=1)
            if data.shape[1] != self.dimensions + self.mnpf:
                raise ValueError(f'given dimension does not match coordinates')

            # get face coordinates and ids
            coords_tmp = np.zeros((data.shape[0], 3)) * 0.
            coords_tmp[:, :self.dimensions] = data[:, :-self.mnpf]  # add column z if 2D
            ids_tmp = self.get_unique_face_ids(data[:, -self.mnpf:])

            # sort and remove doubles
            args = np.unique(ids_tmp, return_index=True)[1].tolist()
            coords_tmp = coords_tmp[args, :]
            ids_tmp = ids_tmp[args]

            # create Nodes
            for i in range(ids_tmp.size):
                mp.CreateNewNode(ids_tmp[i],
                    coords_tmp[i, 0], coords_tmp[i, 1], coords_tmp[i, 2])

        # create CoSimulationInterfaces
        self.interface_input = CoSimulationInterface(self.model, self.settings['interface_input'])
        self.interface_output = CoSimulationInterface(self.model, self.settings['interface_output'])

        # *** HOW TO STORE absolute coordinates at all times?
        """
        adapt Node.X,Y,Z in every step
            in interface object: add DISPLACEMENT to X, Y, Z
        don't make an other object locally in flow solver! overkill
        """

        # create Variables
        self.pressure = vars(KM)['PRESSURE']
        self.traction = vars(KM)['TRACTION']
        self.displacement = vars(KM)['DISPLACEMENT']

        # test simple FSI loop
        for i in range(0):  # *** ON/OFF
            self.set_node_coordinates_test(0.01)
            self.write_node_positions()
            self.send_message('continue')
            self.wait_message('fluent_ready')
        self.send_message('stop')

        print('FINISHED TEST')

    def Initialize(self):
        super().Initialize()
        # *** should sth happen in here?
        # I think most of the init already happends in __init__?
        # when using restart, __init__ is called again??

        self.n = 0  # time step


    def InitializeSolutionStep(self):
        super().InitializeSolutionStep()
        # *** when is this called? also at start??

        # *** should sth happen in here?
        # update index in Fluent? and let Fluent do some more stuff?

        self.k = 0

        self.send_message('next')

    def SolveSolutionStep(self, interface_input):
        # *** do sth about filenames! what should be included in the name??

        # store incoming displacements
        self.interface_input.SetPythonList(interface_input.GetPythonList())

        # update iteration index
        # *** TODO

        # update X,Y,Z in interface
        for key in [_[0] for _ in self.interface_input.model_parts_variables]:
            for node in self.model[key].Nodes:
                disp = node.GetSolutionStepValue(self.displacement)
                node.X += disp[0]
                node.Y += disp[1]
                node.Z += disp[2]

        # write interface data
        self.write_node_positions()

        # let Fluent run, wait for data
        self.send_message('continue')
        # self.wait_message('fluent_ready')

        # read data from Fluent
        for key in self.settings['interface_output'].keys():
            mp = self.model[key]

            # read in datafile
            tmp = 'pressure_traction_thread' + str(mp.thread_id) + '.dat'
            file_name = os.path.join(self.dir_cfd, tmp)
            data = np.loadtxt(file_name, skiprows=1)
            if data.shape[1] != self.dimensions + 1 + self.mnpf:
                raise ValueError('given dimension does not match coordinates')

            # get face coordinates and ids
            traction_tmp = np.zeros((data.shape[0], 3)) * 0.
            traction_tmp[:, :self.dimensions] = data[:, :-1 - self.mnpf]
            pressure_tmp = data[:, self.dimensions]
            ids_tmp = self.get_unique_face_ids(data[:, -self.mnpf:])

            # sort and remove doubles
            args = np.unique(ids_tmp, return_index=True)[1].tolist()
            traction_tmp = traction_tmp[args, :]
            pressure_tmp = pressure_tmp[args]
            ids_tmp = ids_tmp[args]

            # store pressure and traction in Nodes
            if ids_tmp.size != mp.NumberOfNodes():
                raise ValueError('number of nodes does not match size of data')
            index = 0
            for node in mp.Nodes:
                if ids_tmp[index] != node.Id:
                    raise ValueError(f'node IDs do not match: {ids_tmp[index]}, {node.Id}')
                node.SetSolutionStepValue(self.traction, 0, traction_tmp[index, :].tolist())
                node.SetSolutionStepValue(self.pressure, 0, pressure_tmp[index])
                index += 1

        # return interface_output object
        return self.interface_output

    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()

    def Finalize(self):
        super().Finalize()

    def GetInterfaceInput(self):
        return self.interface_input  # *** protect with deepcopy?

    def SetInterfaceInput(self):
        Exception("This solver interface provides no mapping.")

    def GetInterfaceOutput(self):
        return self.interface_output  # *** protect with deepcopy?

    def SetInterfaceOutput(self):
        Exception("This solver interface provides no mapping.")

    def get_unique_face_ids(self, data):
        """
        Construct unique face IDs based on the face's node IDs.

        Parameter data contains a 2D ndarray of node IDs.
        Each row corresponds to the unique node IDs corresponding
        to a certain face, supplemented with -1-values.
        The row is sorted, the -1-values are removed, and then a
        string is made by adding the unique node IDs together.
            e.g. for a row [5, 9, 7, -1, -1]
                 the face ID is "5-7-9"
        """
        data = data.astype(int)
        ids = np.zeros(data.shape[0], dtype='U256')  # array is flattened
        for j in range(ids.size):
            tmp = np.unique(data[j, :])
            if tmp[0] == -1:
                tmp = tmp[1:]
            ids[j] = '-'.join(tuple(tmp.astype(str)))
        return ids

    def set_node_coordinates_test(self, f):
        for key in self.settings['interface_input'].keys():
            for node in self.model[key].Nodes:
                node.Y += (1 - np.cos(2 * np.pi * node.X)) * 0.5 * f

    def write_node_positions(self):
        # *** TODO: add formatter 27.18e or sth...
        for key in self.settings['interface_input'].keys():
            mp = self.model[key]
            file_name = join(self.dir_cfd, f'new_node_coords_thread{mp.thread_id}.dat')
            with open(file_name, 'w') as file:
                file.write(f'{mp.NumberOfNodes()}\n')
                for node in mp.Nodes:
                    if self.dimensions == 2:
                        file.write(f'{node.X} {node.Y} {node.Id}\n')
                    else:
                        file.write(f'{node.X} {node.Y} {node.Z} {node.Id}\n')

    def send_message(self, message):
        file = os.path.join(self.dir_cfd, message + ".msg")
        open(file, 'w').close()
        return

    def wait_message(self, message):
        file = os.path.join(self.dir_cfd, message + ".msg")
        while not os.path.isfile(file):
            time.sleep(0.01)
        os.remove(file)
        return

    def check_message(self, message):
        file = os.path.join(self.dir_cfd, message + ".msg")
        if os.path.isfile(file):
            os.remove(file)
            return True
        return False




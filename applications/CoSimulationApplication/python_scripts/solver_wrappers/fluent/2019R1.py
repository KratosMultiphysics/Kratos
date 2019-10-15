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
        self.dir_cfd = join(os.getcwd(), self.settings['working_directory'].GetString())  # *** alternative for getcwd?
        path_src = os.path.realpath(os.path.dirname(__file__))

        self.remove_all_messages()

        self.cores = self.settings['cores'].GetInt()
        self.case_file = self.settings['case_file'].GetString()
        self.mnpf = self.settings['max_nodes_per_face'].GetInt()
        self.dimensions = self.settings['dimensions'].GetInt()
        self.unsteady = self.settings['unsteady'].GetBool()
        self.hybrid_initialization = self.settings['hybrid_initialization'].GetBool()
        self.flow_iterations = self.settings['flow_iterations'].GetInt()

        self.thread_names = [_.GetString() for _ in self.settings['thread_names'].list()]
        self.n_threads = len(self.thread_names)
        self.thread_ids = [None] * self.n_threads

        # prepare Fluent input journal
        journal = '2019R1.jou'
        thread_names_str = ''
        for key in self.thread_names:
            thread_names_str += ' "' + key + '"'
        unsteady = '#f'
        if self.unsteady:
            unsteady = '#t'
        hybrid_initialization = '#f'
        if self.hybrid_initialization:
            hybrid_initialization = '#t'
        with open(join(path_src, journal), 'r') as infile:
            with open(join(self.dir_cfd, journal), 'w') as outfile:
                for line in infile:
                    line = line.replace('|CASE|', join(self.dir_cfd, self.case_file))
                    line = line.replace('|THREAD_NAMES|', thread_names_str)
                    line = line.replace('|UNSTEADY|', unsteady)
                    line = line.replace('|HYBRID_INITIALIZATION|', hybrid_initialization)
                    line = line.replace('|FLOW_ITERATIONS|', str(self.flow_iterations))
                    outfile.write(line)

        # prepare Fluent UDF
        udf = '2019R1.c'
        with open(join(path_src, udf), 'r') as infile:
            with open(join(self.dir_cfd, udf), 'w') as outfile:
                for line in infile:
                    line = line.replace('|MAX_NODES_PER_FACE|', str(self.mnpf))
                    outfile.write(line)

        # start Fluent with journal
        fluent_gui = self.settings['fluent_gui'].GetBool()
        gui = ''
        if not fluent_gui:
            gui = ' -gu'
        subprocess.Popen(f'fluent 2ddp{gui} -t{self.cores} -i {journal}',  # *** ON/OFF
                             shell=True, executable='/bin/bash', cwd=self.dir_cfd)  # *** ON/OFF

        # get general simulation info from report.sum
        self.wait_message('surface_info_exported')  # *** ON/OFF
        report = join(self.dir_cfd, 'report.sum')
        check = 0
        with open(report, 'r') as file:
            for line in file:
                if check == 2 and 'Time' in line:
                    if 'Steady' in line and self.unsteady:
                        raise ValueError('steady in JSON does not match unsteady Fluent')
                    elif 'Unsteady' in line and not self.unsteady:
                        raise ValueError('unsteady in JSON does not match steady Fluent')
                    break
                if check == 1 and 'Space' in line:
                    if str(self.dimensions) not in line:
                        raise ValueError(f'dimension in JSON does not match Fluent')
                    check = 2
                if 'Model' in line and 'Settings' in line:
                    check = 1

        # get surface thread ID's from report.sum and write them to bcs.txt
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
        with open(join(self.dir_cfd, 'bcs.txt'), 'w') as file:
            file.write(str(len(info)) + '\n')
            for line in info:
                file.write(line + '\n')
        self.send_message('thread_ids_written_to_file')

        # import node and face information
        self.wait_message('nodes_and_faces_stored')  # *** ON/OFF

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
            file_name = join(self.dir_cfd, tmp)
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
            file_name = join(self.dir_cfd, tmp)
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

        # create Variables
        self.pressure = vars(KM)['PRESSURE']
        self.traction = vars(KM)['TRACTION']
        self.displacement = vars(KM)['DISPLACEMENT']

        # test simple FSI loop
        if 0:  # *** ON/OFF
            for i in range(5):
                self.set_node_coordinates_test(0.01)
                self.write_node_positions()
                self.send_message('continue')
                self.wait_message('fluent_ready')
            self.send_message('stop')

    def Initialize(self):
        super().Initialize()
        # TODO I think most of the init already happends in __init__?
        # when using restart, __init__ is called again??

        self.timestep = 0


    def InitializeSolutionStep(self):
        super().InitializeSolutionStep()
        self.iteration = 0
        self.timestep += 1
        print(f'Timestep {self.timestep}')

        self.send_message('next')
        self.wait_message('next_ready')

    def SolveSolutionStep(self, interface_input):
        self.iteration += 1
        print(f'\tIteration {self.iteration}')

        # store incoming displacements
        self.interface_input.SetPythonList(interface_input.GetPythonList())

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
        self.wait_message('continue_ready')

        # read data from Fluent
        for key in self.settings['interface_output'].keys():
            mp = self.model[key]

            # read in datafile
            tmp = f'pressure_traction_timestep{self.timestep}_thread{mp.thread_id}.dat'
            file_name = join(self.dir_cfd, tmp)
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

        if not self.timestep % self.settings['save_iterations'].GetInt():
            self.send_message('save')
            self.wait_message('save_ready')

    def Finalize(self):
        super().Finalize()
        self.send_message('stop')

    def GetInterfaceInput(self):
        return self.interface_input.deepcopy()

    def SetInterfaceInput(self):
        Exception("This solver interface provides no mapping.")

    def GetInterfaceOutput(self):
        return self.interface_output.deepcopy()

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
        for key in self.settings['interface_input'].keys():
            mp = self.model[key]
            tmp = f'nodes_update_timestep{self.timestep}_thread{mp.thread_id}.dat'
            file_name = join(self.dir_cfd, tmp)
            with open(file_name, 'w') as file:
                file.write(f'{mp.NumberOfNodes()}\n')
                for node in mp.Nodes:
                    if self.dimensions == 2:
                        file.write(f'{node.X:27.17e} {node.Y:27.17e} {node.Id:>27}\n')
                    else:
                        file.write(f'{node.X:27.17e} {node.Y:27.17e} {node.Z:27.17e} {node.Id:>27}\n')

    def send_message(self, message):
        file = join(self.dir_cfd, message + ".msg")
        open(file, 'w').close()
        return

    def wait_message(self, message):
        file = join(self.dir_cfd, message + ".msg")
        while not os.path.isfile(file):
            time.sleep(0.01)
        os.remove(file)
        return

    def check_message(self, message):
        file = join(self.dir_cfd, message + ".msg")
        if os.path.isfile(file):
            os.remove(file)
            return True
        return False

    def remove_all_messages(self):
        for file_name in os.listdir(self.dir_cfd):
            if file_name.endswith('.msg'):
                file = join(self.dir_cfd, file_name)
                os.remove(file)

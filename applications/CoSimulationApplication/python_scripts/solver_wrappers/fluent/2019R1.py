import os
import subprocess
import time
import numpy as np
import copy

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

        self.cores = self.settings['cores'].GetInt()
        self.case_file = self.settings['case_file'].GetString()
        self.thread_names = self.settings['thread_names'].parameters  # ***
        self.n_threads = len(self.thread_names)
        self.thread_ids = [None] * self.n_threads
        self.mnpf = self.settings['max_nodes_per_face'].GetInt()
        self.dimensions = self.settings['dimensions'].GetInt()

        path_src = os.path.realpath(os.path.dirname(__file__))

        # prepare Fluent input journal
        journal = '2019R1.jou'
        thread_names_str = ''
        for ft in self.thread_names:
            thread_names_str += ' "' + ft + '"'
        if self.settings['hybrid_initialization'].parameters['a_py_kratos']:  # ***
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

        # *** Parameters file should be updated! can't extract boolean (nor list)

        # start Fluent with journal
        fluent_gui = self.settings['fluent_gui'].parameters['a_py_kratos']  # ***
        gui = ''
        if not fluent_gui:
            gui = ' -gu'
        subprocess.Popen(f'fluent 2ddp{gui} -t{self.cores} -i {journal}',  # *** ON/OFF
                             shell=True, executable='/bin/bash', cwd=self.dir_cfd)  # *** ON/OFF

        # get surface thread ID's from report.sum and write them to bcs.txt
        self.wait_message('surface_info_exported')  # *** ON/OFF
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
        self.wait_message('nodes_and_faces_stored')  # *** ON/OFF

        # import node data, unique sort on ID-string
        self.node_coords = [None] * self.n_threads
        self.node_ids = [None] * self.n_threads
        for t in range(self.n_threads):
            id = self.thread_ids[t]
            file = os.path.join(self.dir_cfd, f'nodes_thread{id}.dat')
            data = np.loadtxt(file, skiprows=1)
            if data.shape[1] != self.dimensions + 1:
                raise ValueError(f'given dimension does not match coordinates')

            coords_tmp = np.zeros((data.shape[0], 3))
            coords_tmp[:, :self.dimensions] = data[:, :-1]  # add column z if 2D
            ids_tmp = data[:, -1].astype(int).astype(str)  # array is flattened

            args = np.unique(ids_tmp, return_index=True)[1].tolist()
            self.node_coords[t] = coords_tmp[args, :]
            self.node_ids[t] = ids_tmp[args]

        # import face data, unique sort on ID-string
        self.face_coords = [None] * self.n_threads
        self.face_ids = [None] * self.n_threads
        for t in range(self.n_threads):
            id = self.thread_ids[t]
            file = os.path.join(self.dir_cfd, f'faces_thread{id}.dat')
            data = np.loadtxt(file, skiprows=1)
            if data.shape[1] != self.dimensions + self.mnpf:
                raise ValueError(f'given dimension does not match coordinates')

            coords_tmp = np.zeros((data.shape[0], 3))
            coords_tmp[:, :self.dimensions] = data[:, :-self.mnpf]  # add column z if 2D
            ids_tmp_all = data[:, -self.mnpf:].astype(int)
            ids_tmp = np.zeros(data.shape[0], dtype='U256')  # array is flattened
            for j in range(ids_tmp.size):
                tmp = np.unique(ids_tmp_all[j, :])
                if tmp[0] == -1:
                    tmp = tmp[1:]
                ids_tmp[j] = '-'.join(tuple(tmp.astype(str)))

            args = np.unique(ids_tmp, return_index=True)[1].tolist()
            self.face_coords[t] = coords_tmp[args, :]
            self.face_ids[t] = ids_tmp[args]

        # create Model, Modelparts, Interfaces
        self.model = cs_data_structure.Model()
        self.model_part_nodes = [None] * self.n_threads
        self.model_part_faces = [None] * self.n_threads
        for t in range(self.n_threads):
            self.model_part_nodes[t] = self.model.CreateModelPart(f'nodes_{self.thread_ids[t]}')
            self.model_part_faces[t] = self.model.CreateModelPart(f'faces_{self.thread_ids[t]}')

            # *** can I use vectors as NodalSolutionStepVariable?? not that I see...

            self.model_part_nodes[t].AddNodalSolutionStepVariable("DISPLACEMENT_X")
            self.model_part_nodes[t].AddNodalSolutionStepVariable("DISPLACEMENT_Y")
            self.model_part_nodes[t].AddNodalSolutionStepVariable("DISPLACEMENT_Z")

            self.model_part_faces[t].AddNodalSolutionStepVariable("TRACTION_X")
            self.model_part_faces[t].AddNodalSolutionStepVariable("TRACTION_Y")
            self.model_part_faces[t].AddNodalSolutionStepVariable("TRACTION_Z")

            self.model_part_faces[t].AddNodalSolutionStepVariable("pressure")

            for j in range(self.node_ids[t].size):
                self.model_part_nodes[t].CreateNewNode(
                    self.node_ids[t][j],
                    self.node_coords[t][j, 0],
                    self.node_coords[t][j, 1],
                    self.node_coords[t][j, 2])

            for j in range(self.face_ids[t].size):
                self.model_part_faces[t].CreateNewNode(
                    self.face_ids[t][j],
                    self.face_coords[t][j, 0],
                    self.face_coords[t][j, 1],
                    self.face_coords[t][j, 2])

            for node in self.model_part_nodes[t].Nodes:
                node.SetSolutionStepValue("DISPLACEMENT_X", 0, 0.0)  # *** init values??
                node.SetSolutionStepValue("DISPLACEMENT_Y", 0, 0.0)
                node.SetSolutionStepValue("DISPLACEMENT_Z", 0, 0.0)

            for node in self.model_part_faces[t].Nodes:
                node.SetSolutionStepValue("TRACTION_X", 0, 0.0)
                node.SetSolutionStepValue("TRACTION_Y", 0, 0.0)
                node.SetSolutionStepValue("TRACTION_Z", 0, 0.0)
                node.SetSolutionStepValue("PRESSURE", 0, 0.0)


            # self.interface_input = CoSimulationInterface(self.model, self.settings["interface_input"])
            # self.interface_output = CoSimulationInterface(self.model, self.settings["interface_output"])

            # *** can I ask for some data about ModelParts? perhaps print it? to check stuff...

            # *** TO DO: define CoSimulationInterfaces... not sure how


        # *** reading in pressure and traction data is analogous to reading in face data
        # *** perhaps put the overlapping code in some kind of read-function??



        # create and write test data to deform threads in Fluent
        node_coords_new = copy.deepcopy(self.node_coords)
        for t in range(self.n_threads):
            node_coords_new[t][:, 1] += self.node_coords[t][:, 0] * 0.
            with open(os.path.join(self.dir_cfd, f'new_node_coords_thread{self.thread_ids[t]}.dat'), 'w') as file:
                file.write(f'{self.node_ids[t].size}\n')
                for i in range(self.node_ids[t].size):
                    for d in range(self.dimensions):
                        file.write(f'{node_coords_new[t][i, d]}\t')
                    file.write(self.node_ids[t][i] + '\n')

        self.send_message('test_position_updates')
        self.wait_message('test_finished')

        # test simple FSI loop
        for i in range(5):
            self.write_position_update_test(0.005)
            self.send_message('continue')
            self.wait_message('fluent_ready')
        self.send_message('stop')

        print('FINISHED TEST')


    def write_position_update_test(self, f):
        for t in range(self.n_threads):
            self.node_coords[t][:, 1] += (1 - np.cos(2 * np.pi * self.node_coords[t][:, 0])) * 0.5 * f
            with open(os.path.join(self.dir_cfd, f'new_node_coords_thread{self.thread_ids[t]}.dat'), 'w') as file:
                file.write(f'{self.node_ids[t].size}\n')
                for i in range(self.node_ids[t].size):
                    for d in range(self.dimensions):
                        file.write(f'{self.node_coords[t][i, d]}\t')
                    file.write(self.node_ids[t][i] + '\n')

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




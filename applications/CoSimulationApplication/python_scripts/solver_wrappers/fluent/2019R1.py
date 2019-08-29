import os
import subprocess
import time
import numpy as np

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

        cores = self.settings['cores'].GetInt()
        case_file = self.settings['case_file'].GetString()
        thread_names = self.settings['face_threads'].parameters  # ***
        n_threads = len(thread_names)
        thread_ids = [-1] * n_threads
        mnpf = self.settings['max_nodes_per_face'].GetInt()

        journal = '2019R1.jou'
        udf = '2019R1.c'
        path_src = os.path.realpath(os.path.dirname(__file__))

        # prepare Fluent input journal
        thread_names_str = ''
        for ft in thread_names:
            thread_names_str += ' "' + ft + '"'
        with open(os.path.join(path_src, journal), 'r') as infile:
            with open(os.path.join(self.dir_cfd, journal), 'w') as outfile:
                for line in infile:
                    line = line.replace('|case|', os.path.join(self.dir_cfd, case_file))
                    line = line.replace('|thread_names|', thread_names_str)
                    outfile.write(line)

        # prepare Fluent UDF
        with open(os.path.join(path_src, udf), 'r') as infile:
            with open(os.path.join(self.dir_cfd, udf), 'w') as outfile:
                for line in infile:
                    line = line.replace('|max_nodes_per_face|', str(mnpf))
                    outfile.write(line)

        # *** Parameters file should be updated! can't extract boolean (nor list)

        # start Fluent with journal
        fluent_gui = self.settings['fluent_gui'].parameters['a_py_kratos']  # ***
        gui = ''
        if not fluent_gui:
            gui = ' -gu'
        subprocess.Popen(f'fluent 2ddp{gui} -t{cores} -i {journal}',
                             shell=True, executable='/bin/bash', cwd=self.dir_cfd)

        # get surface thread ID's from report.sum and write them to bcs.txt
        self.wait_message('surface_info_exported')
        report = os.path.join(self.dir_cfd, 'report.sum')
        check = 0
        info = []
        with open(report, 'r') as file:
            for line in file:
                if check == 3 and line.islower():
                    name, id, _ = line.strip().split()
                    if name in thread_names:
                        info.append(' '.join((name, id)))
                        thread_ids[thread_names.index(name)] = id
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
        self.wait_message('nodes_and_faces_stored')

        # import node data, sort unique on ID-string
        node_coords = [None] * n_threads
        node_ids = [None] * n_threads
        for i in range(n_threads):
            id = thread_ids[i]
            file = os.path.join(self.dir_cfd, f'nodes_thread{id}.dat')
            data = np.loadtxt(file, skiprows=1)

            coords_tmp = data[:, :-1]
            ids_tmp = data[:, -1:].astype(int).astype(str)

            args = np.unique(ids_tmp.flatten(), return_index=True)[1].tolist()

            node_coords[i] = coords_tmp[args, :]
            node_ids[i] = ids_tmp[args, :]

        # import face data, sort unique on ID-string
        face_coords = [None] * n_threads
        face_ids = [None] * n_threads
        for i in range(n_threads):
            id = thread_ids[i]
            file = os.path.join(self.dir_cfd, f'faces_thread{id}.dat')
            data = np.loadtxt(file, skiprows=1)

            coords_tmp = data[:, :-mnpf]
            ids_tmp_all = data[:, -mnpf:].astype(int)
            ids_tmp = np.zeros((ids_tmp_all.shape[0], 1), dtype='U256')

            for j in range(ids_tmp_all.shape[0]):
                tmp = np.unique(ids_tmp_all[j, :])
                if tmp[0] == -1:
                    tmp = tmp[1:]
                ids_tmp[j] = '-'.join(tuple(tmp.astype(str)))

            args = np.unique(ids_tmp.flatten(), return_index=True)[1].tolist()
            face_coords[i] = coords_tmp[args, :]
            face_ids[i] = ids_tmp[args, :]


        for i in range(n_threads):
            print(f'\n\tTHREAD{i}')
            print(node_coords[i][:10, :])
            print(node_ids[i][:10, :])
            print(face_coords[i][:10, :])
            print(face_ids[i][:10, :])




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




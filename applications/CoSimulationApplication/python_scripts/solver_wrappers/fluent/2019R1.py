import os
import subprocess
import time

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
        face_threads = self.settings['face_threads'].parameters  # ***
        max_nodes_per_face = self.settings['max_nodes_per_face'].GetInt()

        journal = '2019R1.jou'
        udf = '2019R1.c'
        path_src = os.path.realpath(os.path.dirname(__file__))


        # prepare Fluent input journal
        face_threads_str = ''
        for ft in face_threads:
            face_threads_str += ' "' + ft + '"'
        with open(os.path.join(path_src, journal), 'r') as infile:
            with open(os.path.join(self.dir_cfd, journal), 'w') as outfile:
                for line in infile:
                    line = line.replace('|case|', os.path.join(self.dir_cfd, case_file))
                    line = line.replace('|face_threads|', face_threads_str)
                    outfile.write(line)

        # prepare Fluent UDF
        with open(os.path.join(path_src, udf), 'r') as infile:
            with open(os.path.join(self.dir_cfd, udf), 'w') as outfile:
                for line in infile:
                    line = line.replace('|max_nodes_per_face|', str(max_nodes_per_face))
                    outfile.write(line)

        # *** Parameters file should be updated! can't extract boolean (nor list)

        # start Fluent with journal
        fluent_gui = self.settings['fluent_gui'].parameters['a_py_kratos']  # ***
        gui = ''
        if not fluent_gui:
            gui = ' -gu'

        if 1:  # *** rm
            subprocess.Popen(f'fluent 2ddp{gui} -t{cores} -i {journal}',
                                 shell=True, executable='/bin/bash', cwd=self.dir_cfd)

        time.sleep(10)
        print('python start message checks')
        self.send_message('python_1')
        self.wait_message('fluent_2')
        print('received fluent_2')
        time.sleep(2)
        self.send_message('python_3')
        while True:
            time.sleep(0.1)
            b = self.check_message('fluent_4')
            print(b)
            if b:
                break
        print('python finished')



        # get surface thread ID's from report.sum and write them to bcs.txt
        report = os.path.join(self.dir_cfd, 'report.sum')
        while not os.path.isfile(report):
            time.sleep(0.01)

        check = 0
        info = []
        with open(report, 'r') as file:
            for line in file:
                if check == 3 and line.islower():
                    tmp = line.strip().split()
                    if tmp[0] in face_threads:
                        info.append(' '.join(tmp[:2]))
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



    # *** test these message functions, communicate a bit with journal

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




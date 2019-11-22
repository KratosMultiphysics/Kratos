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
                    dimension               dimensionality of the problem 2 or 3 
                    arraySize               declare a sufficiently large array size for load array in FORTRAN
                    surfaces                number of interface surfaces
                    cpus                    number of cpus to be used for Abaqus 
                    CSM_dir                 relative path to directory for the files and execution of the flow solver 
                    ramp                    0 for step load, 1 for ramp load in Abaqus
                    deltaT                  time step size
                    timestep_start          time step from which is to be started (initial = 0) 
        """

        self.settings = parameters['settings']
        self.dir_csm = join(os.getcwd(), self.settings['working_directory'].GetString())  # *** alternative for getcwd?
        path_src = os.path.realpath(os.path.dirname(__file__))

        self.remove_all_messages()

        self.cores = self.settings['cores'].GetInt() #  number of cpus Abaqus has to use
        self.dimensions = self.settings['dimensions'].GetInt()
        self.array_size = self.settings["arraysize"].GetInt()
        self.surfaces = self.settings["surfaces"].GetInt()
        self.ramp = self.settings["ramp"].GetInt()
        self.delta_T = self.settings["delta_T"].GetDouble()  #TODO: move to higher-level parameter file?
        self.timestep_start = self.settings["timestep_start"].GetDouble()  #TODO: move to higher-level parameter file?
        self.surfaceIDs = self.settings["surfaceIDs"].GetString()
        self.mp_mode = self.settings["mp_mode"].GetString()
        # Upon(re)starting Abaqus needs to run USRInit.f
        # A restart requires Abaqus to be booted with a restart file

        # prepare abaqus_v6.env
        self.hostnames = []
        self.hostnames_unique = []
        with open(join(self.dir_csm, "AbaqusHosts.txt"), "r") as hostfile:
            for line in hostfile:
                self.hostnames.append(line.rstrip())
                if not line.rstrip() in self.hostnames_unique:
                    self.hostnames_unique.append(line.rstrip())
        self.hostname_replace = ""
        for hostname in self.hostnames_unique:
            self.hostname_replace += "[\'" + hostname + "\', " + str(self.hostnames.count(hostname)) +"], "
        self.hostname_replace = self.hostname_replace.rstrip(", ")
        with open(join(path_src, "abaqus_v6.env"), "r") as infile:
            with open(join(self.dir_csm, "abaqus_v6.env"), "w") as outfile:
                for line in infile:
                    line = line.replace("|HOSTNAME|", self.hostname_replace)
                    line = line.replace("|MP_MODE|", self.mp_mode)
                    line = line.replace("|PID|", str(os.getpid()))
                    line = line.replace("|PWD|", os.path.abspath(os.path.join(self.dir_csm, os.pardir)))
                    line = line.replace("|CSM_dir|", self.settings["working_directory"].GetString())
                    if "|" in line:
                        raise ValueError(f"The following line in abaqus_v6.env still contains a \"|\" after substitution: \n \t{line} \n Probably a parameter was not subsituted")
                    outfile.write(line)

        # prepare Abaqus USRInit.f
        usr = "USRInit.f"
        with open(join(path_src, usr), "r") as infile:
            with open(join(self.dir_csm, "usrInit.f"), "w") as outfile:
                for line in infile:
                    line = line.replace("|dimension|", str(self.dimensions))
                    line = line.replace("|surfaces|", str(self.surfaces))
                    line = line.replace("|cpus|", str(self.cores))

                    #if PWD is too long then FORTRAN code can not compile so this needs special treatment
                    line = self.FORT_replace(line,"|PWD|", os.path.abspath(os.path.join(self.dir_csm, os.pardir)))
                    line = self.FORT_replace(line,"|CSM_dir|", self.settings["working_directory"].GetString())

                    if "|" in line:
                        raise ValueError(f"The following line in USRInit.f still contains a \"|\" after substitution: \n \t{line} \n Probably a parameter was not subsituted")
                    outfile.write(line)

        # compile Abaqus USRInit.f
        cmd_settings1 = "export INTEL_LICENSE_FILE=28518@157.193.126.6"
        cmd_settings2 = "source /apps/SL6.3/Intel/compiler/2015.3.187/bin/compilervars.sh intel64"
        cmd_settings3 = "module load ABAQUS/6.14"
        path_libusr = join(self.dir_csm, "libusr/")
        os.system("rm -r " + path_libusr)
        os.system("mkdir " + path_libusr)
        cmd = "abaqus make library=usrInit.f directory=" + path_libusr + " >> AbaqusSolver.log 2>&1"
        commands = [cmd_settings1, cmd_settings2, cmd_settings3, cmd]
        self.run_shell(self.dir_csm, commands, name='Compile_USRInit')

        # Get loadpoints from usrInit.f
        if(self.timestep_start == 0):
            cmd1 = f"export PBS_NODEFILE=AbaqusHosts.txt && unset SLURM_GTIDS" #To get this to work on HPC?
            cmd2 = f"abaqus job=CSM_Time{self.timestep_start+1} input=CSM_Time{self.timestep_start} cpus=1 user=usrInit.f" \
                f"output_precision=full interactive >> AbaqusSolver.log 2>&1"
            commands = [cmd1, cmd2]
            self.run_shell(self.dir_csm, commands, name='Abaqus_USRInit_Time0')
        else:
            cmd1 = f"export PBS_NODEFILE=AbaqusHosts.txt && unset SLURM_GTIDS" #To get this to work on HPC?
            cmd2 = f"abaqus job=CSM_Time{self.timestep_start+1} oldjob=CSM_Time{self.timestep_start} input=CSM_Restart cpus=1 user=usrInit.f" \
                f"output_precision=full interactive >> AbaqusSolver.log 2>&1"
            commands = [cmd1, cmd2]
            self.run_shell(self.dir_csm, commands, name=f'Abaqus_USRInit_Restart')

        # prepare GetOutput.cpp
        get_output = "GetOutput.cpp"
        with open(join(path_src, get_output), "r") as infile:
            with open(join(self.dir_csm, get_output), "w") as outfile:
                for line in infile:
                    line = line.replace("|surfaces|", str(self.surfaces))
                    line = line.replace("|surfaceIDs|", str(self.surfaceIDs))
                    line = line.replace("|dimension|", str(self.dimensions))
                    if "|" in line:
                        raise ValueError(f"The following line in GetOutput.cpp still contains a \"|\" after substitution: \n \t{line} \n Probably a parameter was not subsituted")

                    outfile.write(line)

        # compile GetOutput.cpp
        cmd = "abaqus make job=GetOutput user=GetOutput.cpp >> AbaqusSolver.log 2>&1"
        commands = [cmd_settings1, cmd_settings2, cmd_settings3, cmd]
        self.run_shell(self.dir_csm, commands, name='Compile_GetOutput')

        # Get node positions (not load points) at startTimeStep

        # prepare Abaqus USR.f
        usr = "USR.f"
        with open(join(path_src, usr), "r") as infile:
            with open(join(self.dir_csm, "usr.f"), "w") as outfile:
                for line in infile:
                    line = line.replace("|dimension|", str(self.dimensions))
                    line = line.replace("|arraySize|", str(self.array_size))
                    line = line.replace("|surfaces|", str(self.surfaces))
                    line = line.replace("|cpus|", str(self.cores))
                    line = line.replace("|ramp|", str(self.ramp))
                    line = line.replace("|deltaT|", str(self.delta_T))

                    # if PWD is too ling then FORTRAN code can not compile so this needs special treatment
                    line = self.FORT_replace(line, "|PWD|", os.path.abspath(os.path.join(self.dir_csm, os.pardir)))
                    line = self.FORT_replace(line, "|CSM_dir|", self.settings["working_directory"].GetString())
                    if "|" in line:
                        raise ValueError(f"The following line in USR.f still contains a \"|\" after substitution: \n \t{line} \n Probably a parameter was not subsituted")

                    outfile.write(line)

        # compile Abaqus USR.f
        os.system("rm -r " + path_libusr)  # remove libusr containing compiled USRInit.f
        os.system("mkdir " + path_libusr)
        cmd = "abaqus make library=usr.f directory=" + path_libusr + " >> AbaqusSolver.log 2>&1"
        commands = [cmd_settings1, cmd_settings2, cmd_settings3, cmd]
        self.run_shell(self.dir_csm, commands, name='Compile_USR')



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

    def run_shell(self, work_dir, commands, wait=True, name='script'):
        script = f'{work_dir}/{name}.sh'
        with open(script, 'w') as file:
            file.write('#!/bin/bash\n')
            file.write(f'cd {work_dir}\n')
            for line in commands:
                file.write(line + '\n')
        os.chmod(script, 0o700)
        if wait:
            p = subprocess.call(script, shell=True)
        else:
            p = subprocess.Popen(script, shell=True)
        return p

    def FORT_replace(self, line, orig, new):
        '''The length of a line in FORTRAN 77 is limited, replacing working directories can exceed this limiet
        This functions splits these strings over multiple lines'''

        ampersand_location = 6
        char_limit = 72

        if "|" in line:
            temp = line.replace(orig, new)
            N = len(temp)

            if N > char_limit:
                count = 0
                line = ""
                line += temp[0:char_limit] + "\n"
                count += char_limit
                while count < N:
                    temp_string = temp[count:count + char_limit - 12]
                    n = len(temp_string)
                    count += n
                    if count < N:  # need to append an additional new line
                        line += "     &" + "      " + temp_string + "\n"
                    else:
                        line += "     &" + "      " + temp_string
            else:
                line = temp

        return line
import os
from os.path import join
import subprocess
import time
import numpy as np
import copy
import re


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
                    input_file              Name of the Abaqus input file (Located in the directory where Python is launched)
                    dimension               dimensionality of the problem 2 or 3 
                    arraySize               declare a sufficiently large array size for load array in FORTRAN
                    surfaces                number of interface surfaces
                    cpus                    number of cpus to be used for Abaqus 
                    CSM_dir                 relative path to directory for the files and execution of the structural solver 
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
        # self.surfaceIDs = self.settings["surfaceIDs"].GetString()
        self.surfaceIDs = [_.GetString() for _ in self.settings['surfaceIDs'].list()]
        self.n_surfaces = len(self.surfaceIDs)
        self.mp_mode = self.settings["mp_mode"].GetString()
        self.input_file = self.settings["input_file"].GetString()

        # Upon (re)starting Abaqus needs to run USRInit.f
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

        #Create start and restart file
        self.write_start_and_restart_inp(self.input_file, self.dir_csm+"/CSM_Time0.inp", self.dir_csm+"/CSM_Restart.inp")

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
        path_libusr = join(self.dir_csm, "libusr/")
        os.system("rm -r " + path_libusr)
        os.system("mkdir " + path_libusr)
        cmd = "abaqus make library=usrInit.f directory=" + path_libusr + " >> AbaqusSolver.log 2>&1"
        commands = [cmd]
        self.run_shell(self.dir_csm, commands, name='Compile_USRInit')

        # Get loadpoints from usrInit.f
        if(self.timestep_start == 0):
            cmd1 = f"export PBS_NODEFILE=AbaqusHosts.txt && unset SLURM_GTIDS" #To get this to work on HPC?
            cmd2 = f"rm CSM_Time{self.timestep_start}Surface*Faces.dat CSM_Time{self.timestep_start}Surface*FacesBis.dat"
            cmd3 = f"abaqus job=CSM_Time{self.timestep_start+1} input=CSM_Time{self.timestep_start} cpus=1 user=usrInit.f" \
                f" output_precision=full interactive >> AbaqusSolver.log 2>&1"
            commands = [cmd1, cmd2, cmd3]
            self.run_shell(self.dir_csm, commands, name='Abaqus_USRInit_Time0')
        else:
            cmd1 = f"export PBS_NODEFILE=AbaqusHosts.txt && unset SLURM_GTIDS" #To get this to work on HPC?
            cmd2 = f"rm CSM_Time{self.timestep_start}Surface*Faces.dat CSM_Time{self.timestep_start}Surface*FacesBis.dat"
            cmd3 = f"abaqus job=CSM_Time{self.timestep_start+1} oldjob=CSM_Time{self.timestep_start} input=CSM_Restart cpus=1 user=usrInit.f" \
                f" output_precision=full interactive >> AbaqusSolver.log 2>&1"
            commands = [cmd1, cmd2, cmd3]
            self.run_shell(self.dir_csm, commands, name=f'Abaqus_USRInit_Restart')

        # prepare GetOutput.cpp
        get_output = "GetOutput.cpp"
        temp_str = ""
        for j in range(0, self.n_surfaces - 1):
            temp_str += f"\"{self.surfaceIDs[j]}\", "
        temp_str += f"\"{self.surfaceIDs[self.n_surfaces-1]}\""

        with open(join(path_src, get_output), "r") as infile:
            with open(join(self.dir_csm, get_output), "w") as outfile:
                for line in infile:
                    line = line.replace("|surfaces|", str(self.surfaces))
                    line = line.replace("|surfaceIDs|", temp_str)
                    line = line.replace("|dimension|", str(self.dimensions))
                    if "|" in line:
                        raise ValueError(f"The following line in GetOutput.cpp still contains a \"|\" after substitution: \n \t{line} \n Probably a parameter was not subsituted")

                    outfile.write(line)

        # compile GetOutput.cpp
        cmd = "abaqus make job=GetOutput user=GetOutput.cpp >> AbaqusSolver.log 2>&1"
        commands = [cmd]
        self.run_shell(self.dir_csm, commands, name='Compile_GetOutput')

        # Get node positions (not load points) at startTimeStep
        cmd = f"abaqus ./GetOutput.exe CSM_Time{self.timestep_start+1} 0 >> AbaqusSolver.log 2>&1"
        commands=[cmd]
        self.run_shell(self.dir_csm, commands, name='GetOutput_Start')

        for i in range(0,self.surfaces):
            path_output=f"CSM_Time{self.timestep_start+1}Surface{i}Output.dat"
            path_nodes=f"CSM_Time{self.timestep_start}Surface{i}Nodes.dat"
            cmd=f"mv {path_output} {path_nodes}"
            commands = [cmd]
            self.run_shell(self.dir_csm, commands, name='Move_Output_File_To_Node')

            # Create elements file per surface
            face_file = os.path.join(self.dir_csm,f"CSM_Time{self.timestep_start}Surface{i}Cpu0Faces.dat")
            output_file = os.path.join(self.dir_csm,f"CSM_Time{self.timestep_start}Surface{i}Elements.dat")
            self.makeElements(face_file, output_file)

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
        commands = [cmd]
        self.run_shell(self.dir_csm, commands, name='Compile_USR')

        # ### --- Create Model --- ###
        # self.model = cs_data_structure.Model()
        #
        # # create ModelParts
        # for key, value in (self.settings['interface_input'].items() +
        #                    self.settings['interface_output'].items()):
        #     # add ModelPart to Model
        #     self.model.CreateModelPart(key)
        #     mp = self.model[key]
        #
        #     # add historical variables to ModelPart
        #     for var_name in value.list():
        #         var = vars(KM)[var_name.GetString()]
        #         mp.AddNodalSolutionStepVariable(var)
        #
        #     # add information to ModelPart
        #     for i in range(self.n_threads):
        #         if self.thread_names[i] in key:
        #             mp.thread_name = self.thread_names[i]
        #             mp.thread_id = self.thread_ids[i]
        #             if 'thread_id' not in dir(mp):
        #                 raise AttributeError('could not find thread name corresponding to key')

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

    def makeElements(self,face_file, output_file):
        firstLoop = 1
        element = 0
        element_prev = -1
        point = 0
        point_prev =-1
        element_0 = -1
        point_0 = -1
        count=0
        element_str = ""

        with open(face_file, 'r') as file:
            for line in file:
                values = line.strip().split()
                element = int(values[0])
                point = int(values[1])
                if element == element_0 and point == point_0:
                    break
                if element == element_prev:
                    if point == point_prev + 1:
                        point_prev = point
                    else:
                        raise ValueError(f"loadpoint number increases by more than one per line for element {element}")
                else:
                    if point == 1:
                        point_prev = point
                        element_prev = element
                        element_str += str(element) + "\n"
                        count+=1
                        if firstLoop: #Faces contains all values multiple times, but we only want it once
                            element_0=element
                            point_0 = point
                            firstLoop=0
                    else:
                        raise ValueError(
                            f"loadpoint number does not start at 1 for element {element}")

        element_str = f"{count}\n{point_prev}\n"+element_str
        with open(output_file,"w") as file:
            file.write(element_str)

    def run_shell(self, work_dir, commands, wait=True, name='script', delete=True):
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
        if delete:
            os.system("rm " + script)
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
                    temp_string = temp[count:count + char_limit - 6]
                    n = len(temp_string)
                    count += n
                    if count < N:  # need to append an additional new line
                        line += "     &" + temp_string + "\n"
                    else:
                        line += "     &" + temp_string
            else:
                line = temp

        return line

    def write_start_and_restart_inp(self,input_file,output_file,restart_file):
        bool_restart = 0

        rf = open(restart_file,"w")
        of = open(output_file,"w")

        rf.write("*HEADING \n")
        rf.write("*RESTART, READ \n")

        with open(input_file) as f:
            line=f.readline()
            while line:
                if "*step" in line.lower():
                    contents=line.split(",") #Split string on commas
                    for s in contents:
                        if s.strip().startswith("inc="):
                            numbers=re.findall("\d+",s)
                            if int(numbers[0])!=1:
                                raise NotImplementedError(f"inc={numbers[0]}: currently only single increment steps are implemented for the Abaqus wrapper")
                    of.write(line)
                    if bool_restart:
                        rf.write(line)
                    line = f.readline()
                elif "*dynamic" in line.lower():
                    contents = line.split(",")
                    for s in contents:
                        if "application" in s.lower():
                            contents_B = s.split("=")
                            if contents_B[1].lower().strip() != "quasi-static":
                                raise NotImplementedError(
                                    f"{contents_B[1]} not available: Currently only quasi-static is implemented for the Abaqus wrapper")
                    of.write(line)
                    if bool_restart:
                        rf.write(line)
                    line=f.readline() #need to skip the next line
                    of.write(f"{self.delta_T}, {self.delta_T},\n") #Change the time step in the Abaqus step
                    if bool_restart:
                        rf.write(f"{self.delta_T}, {self.delta_T},\n") #Change the time step in the Abaqus step (restart-file)
                    line = f.readline()
                else:
                    of.write(line)
                    if bool_restart:
                        rf.write(line)
                    line=f.readline()
                if "** --"in line:
                    bool_restart=1
        rf.close()
        of.close()
# CoSimulation imports
from KratosMultiphysics.CoSimulationApplication.function_callback_utility import GenericCallFunction
import KratosMultiphysics as KM

# Other imports
import numpy as np
import json
import os
from pathlib import Path

class SDoFSolver(object):
    """ This class implements an SDof solver, independent of Kratos
    Several types of applying loads are available
    """
    def __init__(self, input_name, model_part):

        # mimicing two constructors
        if isinstance(input_name, dict):
            parameters = input_name

        elif isinstance(input_name, str):
            if not input_name.endswith(".json"):
                input_name += ".json"

            with open(input_name,'r') as ProjectParameters:
                parameters = json.load(ProjectParameters)

        else:
            raise Exception("The input has to be provided as a dict or a string")

        default_settings = {
                "system_parameters":{
                    "mass"      : 100.0,
                    "stiffness" : 4000.0,
                    "damping"   : 0.0,
                    "modulus_self_weight": 9.81
                },
                "time_integration_parameters":{
                    "alpha_m"   : -0.3,
                    "start_time": 0.0,
                    "time_step" : 0.05,
                },
                "initial_values":{
                    "displacement"  : 0.0,
                    "velocity"      : 0.0,
                    "acceleration"  : 0.0
                },
                "boundary_conditions":{
                    "load_impulse" : 0.0,
                    "omega_force"        : 0.0,
                    "omega_root_point_displacement"        : 0.0,
                    "excitation_function_force": "A * sin(omega * t)",
                    "excitation_function_root_point_displacement": "A * sin(omega * t)",
                    "amplitude_root_point_displacement": 0.0,
                    "amplitude_force": 0.0
                },
                "solver_parameters": {
                    "buffer_size"   : 2
                },
                "output_parameters":{
                    "write_output_file": True,
                    "file_name" : "sdof_solver/results_sdof.dat"
                },
                "restart_parameters":{
                    "load_restart_file": False,
                    "input_file_name": "",
                    "write_restart_file": False,
                    "restart_save_frequency": 0.0,
                    "restart_control_type": "time",
                    "save_restart_files_in_folder"   : True,
                    "output_path"                  : "sdof_solver",
                    "max_files_to_keep"            : -1
                }
                }

        RecursivelyValidateAndAssignDefaults(default_settings, parameters)

        self.mass = parameters["system_parameters"]["mass"]
        self.stiffness = parameters["system_parameters"]["stiffness"]
        self.damping = parameters["system_parameters"]["damping"]
        self.modulus_self_weight = parameters["system_parameters"]["modulus_self_weight"]


        self.alpha_m = parameters["time_integration_parameters"]["alpha_m"]
        self.delta_t = parameters["time_integration_parameters"]["time_step"]
        self.start_time = parameters["time_integration_parameters"]["start_time"]

        self.initial_displacement = parameters["initial_values"]["displacement"]
        self.initial_velocity = parameters["initial_values"]["velocity"]

        self.excitation_function_force = parameters["boundary_conditions"]["excitation_function_force"]
        self.excitation_function_root_point_displacement = parameters["boundary_conditions"]["excitation_function_root_point_displacement"]
        self.load_impulse = parameters["boundary_conditions"]["load_impulse"]
        self.omega_force = parameters["boundary_conditions"]["omega_force"]
        self.omega_root_point_displacement = parameters["boundary_conditions"]["omega_root_point_displacement"]
        self.amplitude_root_point_displacement = parameters["boundary_conditions"]["amplitude_root_point_displacement"]
        self.amplitude_force = parameters["boundary_conditions"]["amplitude_force"]

        #calculate initial acceleration
        factor = self.load_impulse - self.stiffness * self.initial_displacement
        self.initial_acceleration = (1/self.mass) * factor

        self.beta = 0.25 * (1- self.alpha_m)**2
        self.gamma =  0.50 - self.alpha_m

        self.LHS = np.array([[1.0, 0.0, -self.delta_t**2 * self.beta],
                             [0.0, 1.0, -self.delta_t * self.gamma],
                             [self.stiffness,
                              self.damping,
                               (1-self.alpha_m) * self.mass]])

        self.RHS_matrix = np.array([[1.0, self.delta_t, self.delta_t**2 * (0.5 - self.beta)],
                                    [0.0, 1.0, self.delta_t*(1-self.gamma)],
                                    [0.0,
                                     0.0,
                                     -self.alpha_m * self.mass]])

        self.buffer_size = parameters["solver_parameters"]["buffer_size"]
        self.output_file_name = parameters["output_parameters"]["file_name"]
        self.write_output_file = parameters["output_parameters"]["write_output_file"]

        ## Restart
        self.model_part = model_part
        self.load_restart_file = parameters["restart_parameters"]["load_restart_file"]
        if self.load_restart_file:
            self.input_file_name = parameters["restart_parameters"]["input_file_name"]

        self.write_restart_file = parameters["restart_parameters"]["write_restart_file"]

        if self.write_restart_file:
            self.restart_output_path = parameters["restart_parameters"]["output_path"]
            self.restart_save_frequency = parameters["restart_parameters"]["restart_save_frequency"]

            self.save_restart_files_in_folder   = parameters["restart_parameters"]["save_restart_files_in_folder"]


            restart_control_type = parameters["restart_parameters"]["restart_control_type"]
            if restart_control_type == "time":
                self.restart_control_type_is_time = True
            elif restart_control_type == "step":
                self.restart_control_type_is_time = False
                self.restart_save_frequency = int(self.restart_save_frequency) # STEP is an integer
            else:
                err_msg =  'The requested restart_control_type "' + restart_control_type + '" is not available!\n'
                err_msg += 'Available options are: "time", "step"'
                raise Exception(err_msg)

            self.next_output = self.restart_save_frequency # Schedule the first output to avoid printing in first step

            self.max_files_to_keep = parameters["restart_parameters"]["max_files_to_keep"]

            self.raw_path = os.getcwd()
       
    def Initialize(self):
        #solution buffer
        self.x = np.zeros((3, self.buffer_size))
        #values at the root point buffer
        self.x_f = np.zeros((3, self.buffer_size))

        initial_values = np.array([self.initial_displacement,
                                   self.initial_velocity,
                                   self.initial_acceleration])
        self.dx = initial_values
        self.dx_f = np.zeros(3)
        self.time = self.start_time

        self.root_point_displacement = 0.0

        #restart
        if self.load_restart_file:
            self.LoadRestart()

        #x and dx contain: [displacement, velocity, acceleration]
        if self.write_output_file:            
            self.InitializeOutput()

        #apply external load as an initial impulse
        self.load_vector = np.array([0,
                                     0,
                                     self.load_impulse])
    
    def InitializeOutput(self):
        data_comm = KM.DataCommunicator.GetDefault()
        if data_comm.Rank()==0:
            
            # Create the directory where the results will be stored
            output_file_path = os.path.dirname()
            if not os.path.exists(output_file_path):
                os.makedirs(output_file_path)

            # erase results files from other runs
            if os.path.isfile(self.output_file_name):
                os.remove(self.output_file_name)

            with open(self.output_file_name, "w") as results_sdof:
                results_sdof.write("time"+ " " +
                                   "displacement" + " " +
                                   "velocity" + " " +
                                   "acceleration" + " " +
                                   "root point displacement" + " " +
                                   "root point velocity" + " " +
                                   "root point acceleration" + " " +
                                   "relative displacement" + " " +
                                   "relative velocity" + " " +
                                   "relative accleration" + " " +
                                   "reaction" + "\n")
            
            self.OutputSolutionStep()

            #restart
            if self.write_restart_file:
                self.CreateOutputFolder()

    def OutputSolutionStep(self):
        data_comm = KM.DataCommunicator.GetDefault()
        if data_comm.Rank()==0:
            reaction = self.CalculateReaction()
            if self.write_output_file:
                with open(self.output_file_name, "a") as results_sdof:
                    #outputs results
                    results_sdof.write(str(np.around(self.time, 3)) + " " +
                                    str(self.dx[0]) + " " +
                                    str(self.dx[1]) + " " +
                                    str(self.dx[2]) + " " +
                                    str(self.dx_f[0]) + " " +
                                    str(self.dx_f[1]) + " " +
                                    str(self.dx_f[2]) + " " +
                                    str(self.dx[0] - self.dx_f[0]) + " " +
                                    str(self.dx[1] - self.dx_f[1]) + " " +
                                    str(self.dx[2] - self.dx_f[2]) + " " +
                                    str(reaction) + "\n")

            #restart
            if self.write_restart_file:
                if self.IsRestartOutputStep():
                    self.SaveRestart()

    def AdvanceInTime(self, current_time):
        # similar to the Kratos CloneTimeStep function
        # advances values along the buffer axis (so rolling columns) using numpy's roll
        self.x = np.roll(self.x,1,axis=1)
        self.x_f = np.roll(self.x_f,1,axis=1)
        # overwriting at the buffer_idx=0 the newest values
        buffer_idx = 0
        self.x[:,buffer_idx] = self.dx
        self.x_f[:,buffer_idx] = self.dx_f

        self.time = current_time + self.delta_t
        return self.time

    def CalculateEquivalentForceFromRootPointExcitation(self, d_f):
        #d_f = self.root_point_displacement
        v_f = self.x_f[1,0] + self.delta_t * (self.gamma * d_f + (1-self.gamma) * self.x_f[2,0])
        a_f = 1/(self.delta_t**2 * self.beta) * (d_f - self.x_f[0,1])\
            - 1/(self.delta_t * self.beta) * self.x_f[1,0]\
            + (1-1/(2*self.beta)) * self.x_f[2,0]
        self.dx_f = np.array([d_f, v_f, a_f])
        b_f = np.array([0.0, 0.0, d_f * self.stiffness + v_f * self.damping])
        return b_f

    def ApplyRootPointExcitation(self):
        scope_vars = {'t' : self.time, 'omega': self.omega_root_point_displacement, 'A': self.amplitude_root_point_displacement}
        return GenericCallFunction(self.excitation_function_root_point_displacement, scope_vars, check=False)

    def ApplyForceExcitation(self):
        scope_vars = {'t' : self.time, 'omega': self.omega_force, 'A': self.amplitude_force}
        return GenericCallFunction(self.excitation_function_force, scope_vars, check=False)

    def SolveSolutionStep(self):
        b = self.RHS_matrix @ self.x[:,0]
        #external load
        excitation_load = self.ApplyForceExcitation()
        self.load_vector[-1] += excitation_load
        # print("external load= ", self.load_vector[-1])
        b += self.load_vector
        #root point displacement
        d_f_excitation = self.ApplyRootPointExcitation()
        d_f = d_f_excitation + self.root_point_displacement
        b_f = self.CalculateEquivalentForceFromRootPointExcitation(d_f)
        b += b_f
        self.dx = np.linalg.solve(self.LHS, b)

    def CalculateReaction(self, buffer_idx=0):
        reaction = self.damping * ( self.dx[1] - self.dx_f[1]) \
                 + self.stiffness * (self.dx[0] - self.dx_f[0])
        return reaction

    def CalculateSelfWeight(self):
        self_weight = self.mass * self.modulus_self_weight
        return self_weight

    def GetSolutionStepValue(self, identifier, buffer_idx=0):
        if identifier == "DISPLACEMENT":
            return self.x[:,buffer_idx][0]
        elif identifier == "VELOCITY":
            return self.x[:,buffer_idx][1]
        elif identifier == "ACCELERATION":
            return self.x[:,buffer_idx][2]
        elif identifier == "REACTION":
            return self.CalculateReaction()
        elif identifier == "VOLUME_ACCELERATION":
            return self.CalculateSelfWeight()
        else:
            raise Exception("Identifier is unknown!")

    def SetSolutionStepValue(self, identifier, value, buffer_idx=0):
        if identifier == "DISPLACEMENT":
            self.x[:,buffer_idx][0] = value
        elif identifier == "VELOCITY":
            self.x[:,buffer_idx][1] = value
        elif identifier == "ACCELERATION":
            self.x[:,buffer_idx][2] = value
        elif identifier == "LOAD":
            self.load_vector[-1] = 0.0
            self.load_vector[-1] = value
        elif identifier == "ROOT_POINT_DISPLACEMENT":
            self.root_point_displacement = 0.0
            self.root_point_displacement = value
        else:
            raise Exception("Identifier is unknown!")

    ####### RESTART
    def LoadRestart(self):
        input_file = self.__GetFolderPathLoad()

        data = json.load(open(input_file + '.json'))

        self.time = data["time"]
        
        self.dx[0] = data["displacement"]
        self.dx[1] = data ["velocity"]
        self.dx[2] = data ["acceleration"]
        
        self.dx_f[0] = data["displacement_r"]
        self.dx_f[1] = data ["velocity_r"]
        self.dx_f[2] = data ["acceleration_r"]

        print(self.dx)
        print(self.dx_f)

        self.model_part.ProcessInfo[KM.IS_RESTARTED] = True

        KM.Logger.PrintInfo("Restart Utility", "Finished loading model part from restart file.")       

    def SaveRestart(self):
        """
        This function saves the restart file. It should be called at the end of a time-step.
        Use "IsRestartOutputStep" to check if a restart file should be written in this time-step
        """

        if self.restart_control_type_is_time:
            time = self.time
            control_label = self.__GetPrettyTime(time)
        else:
            pass
            # control_label = self.model_part.ProcessInfo[KM.STEP]

        self.time = time
        self.displacement = self.dx[0]
        self.velocity = self.dx[1]
        self.acceleration = self.dx[2]

        self.displacement_root = self.dx_f[0]
        self.velocity_root = self.dx_f[1]
        self.acceleration_root = self.dx_f[2]
        
        self.dictionary = {
            "time": self.time, "displacement": self.displacement,"displacement_r": self.displacement_root,
            "velocity": self.velocity,"velocity_r": self.velocity_root, "acceleration": self.acceleration,"acceleration_r": self.acceleration_root}
        with open(self.__GetFolderPathSave()+'//Sdof_' + str(self.__GetPrettyTime(time)) + '.json', 'w') as outfile:
            json.dump(self.dictionary, outfile)
        
        # Schedule next output
        if self.restart_save_frequency > 0.0: # Note: if == 0, we'll just always print
            while self.next_output <= control_label:
                self.next_output += self.restart_save_frequency      

        # Cleanup
        self._ClearObsoleteRestartFiles()
        
    def IsRestartOutputStep(self):
        """
        This function checks and returns whether a restart file should be written in this time-step
        """
        if self.restart_control_type_is_time:
            return (self.time > self.next_output)
        else:
            pass
            # return (self.model_part.ProcessInfo[KM.STEP] >= self.next_output)

    def CreateOutputFolder(self):
        if self.save_restart_files_in_folder:
            folder_path = self.__GetFolderPathSave()
            if not os.path.isdir(folder_path) and self.model_part.GetCommunicator().MyPID() == 0:
                os.makedirs(folder_path)
            self.model_part.GetCommunicator().GetDataCommunicator().Barrier()

    def _ClearObsoleteRestartFiles(self):
        """Delete restart files that are no longer needed."""
        if self.max_files_to_keep > -1:
            self.restart_files = self._GetRestartFiles()

            number_of_obsolete_files = len(self.restart_files) - self.max_files_to_keep
            restart_file_keys = sorted(self.restart_files)

            for i in range(number_of_obsolete_files):
                i_key = restart_file_keys[i]
                file_path = os.path.join(self.__GetFolderPathSave(), self.restart_files[i_key])
                KM.kratos_utilities.DeleteFileIfExisting(file_path)

    def _GetFileNamePattern(self):
        """Return the pattern of flags in the file name for FileNameDataCollector."""
        file_name_pattern = "<model_part_name>"
        if self.restart_control_type_is_time:
            file_name_pattern += "_<time>"
        else:
            file_name_pattern += "_<step>"
        file_name_pattern += ".json"
        return file_name_pattern

    def __GetFolderPathLoad(self):
        return os.path.join(self.raw_path, self.input_file_name)
        
    def __GetFolderPathSave(self):
        if self.save_restart_files_in_folder:
            return os.path.join(self.raw_path, self.restart_output_path)
        else:
            return self.raw_path

    def __GetPrettyTime(self, time):
        """This functions reduces the digits of a number to a relevant precision
        """
        pretty_time = "{0:.12g}".format(time)
        pretty_time = float(pretty_time)
        return pretty_time

    def _GetRestartFiles(self):
        """Return a dictionary of stepID - restart_file_list dictionary that stores sets of restart files for each step."""
        restart_path    = Path(self.__GetFolderPathSave())
        restart_files   = {}
        if restart_path.is_dir():
            file_name_data_collector = KM.FileNameDataCollector(self.model_part, os.path.join(self.__GetFolderPathSave(), self._GetFileNamePattern()), {})

            for file_name_data in file_name_data_collector.GetFileNameDataList():
                # Get step id
                if self.restart_control_type_is_time:
                    step_id = file_name_data.GetTime()
                else:
                    step_id = file_name_data.GetStep()

                self._UpdateRestartFilesMap(restart_files, step_id, file_name_data)

            # barrier is necessary to avoid having some ranks deleting files while other ranks still detect them in the same directory
            self.model_part.GetCommunicator().GetDataCommunicator().Barrier()

        return restart_files

    def _UpdateRestartFilesMap(self, restart_files_map, step_id, file_name_data):
        restart_files_map[step_id] = file_name_data.GetFileName()

def ValidateAndAssignDefaults(defaults, settings, recursive=False):
    for key, val in settings.items():
        # check if the current entry also exists in the defaults
        if not key in defaults.keys():
            err_msg  = 'The item with name "' + key + '" is present in this '
            err_msg += 'settings\nbut NOT in the defaults!\n'
            err_msg += 'settings are:\n'
            err_msg += json.dumps(settings, indent=4)
            err_msg += '\ndefaults are:\n'
            err_msg += json.dumps(defaults, indent=4)
            raise Exception(err_msg)

        # check if the type is the same in the defaults
        if type(settings[key]) != type(defaults[key]):
            err_msg  = 'The type of the item with name "' + key + '" (type: "'
            err_msg += str(type(settings[key]).__name__)+'") in this '
            err_msg += 'settings\nis NOT the same as in the defaults (type: "'
            err_msg += str(type(defaults[key]).__name__)+'")!\n'
            err_msg += 'settings are:\n'
            err_msg += json.dumps(settings, indent=4)
            err_msg += '\ndefaults are:\n'
            err_msg += json.dumps(defaults, indent=4)
            raise Exception(err_msg)

    # loop the defaults and add the missing entries
    for key_d, val_d in defaults.items():
        if key_d not in settings: # add the default in case the setting is not present
            settings[key_d] = val_d
        elif recursive and type(val_d) is dict:
            RecursivelyValidateAndAssignDefaults(val_d, settings[key_d])

def RecursivelyValidateAndAssignDefaults(defaults, settings):
    ValidateAndAssignDefaults(defaults, settings, recursive=True)

def is_restarted(self):
        # this function avoids the long call to ProcessInfo and is also safer
        # in case the detection of a restart is changed later
        return self.main_model_part.ProcessInfo[KM.IS_RESTARTED]
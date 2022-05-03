import os
from pyevtk import hl
from pyevtk import vtk
import numpy as np
import shutil
import KratosMultiphysics
import KratosMultiphysics.ParticleMechanicsApplication as KratosParticle

# Import time library
from time import time


def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    model_part = Model[settings["Parameters"]["model_part_name"].GetString()]
    parameters = settings["Parameters"]
    return ParticleVTKOutputProcess(model_part, parameters)

class ParticleVTKOutputProcess(KratosMultiphysics.Process):
    defaults = KratosMultiphysics.Parameters("""{
        "model_part_name"                    : "MPM_Material",
        "output_control_type"                : "step",
        "output_frequency"                   : 1,
        "file_format"                        : "ascii",
        "output_precision"                   : 7,
        "folder_name"                        : "vtk_output",
        "output_sub_model_parts"             : false,
        "save_output_files_in_folder"        : true,
        "gauss_point_results" : []
    }""")

    def __init__(self, model_part, param):
        KratosMultiphysics.Process.__init__(self)

        if param is None:
            param = self.defaults
        else:
            param.ValidateAndAssignDefaults(self.defaults)

        self.param = param
        self.model_part = model_part
        self.problem_name = param["model_part_name"].GetString()

        self.step_count = 0
        self.printed_step_count = 0
        self.next_output = 0.0

        self.coords_X = np.empty(0)
        self.coords_Y = np.empty(0)
        self.coords_Z = np.empty(0)
        self.temp_results = np.empty([0,3])
        self.result_names = []
        self.result_dict = {}

        self.vtk_post_path_directory = os.path.join(param["folder_name"].GetString())
        shutil.rmtree(self.vtk_post_path_directory, ignore_errors=True)
        os.makedirs(str(self.vtk_post_path_directory))


        # Public Functions
    def ExecuteInitialize(self):
        # Set up output frequency and format
        output_control_type = self.param["output_control_type"].GetString()
        if output_control_type == "time":
            self.output_control_is_time = True
        elif output_control_type == "step":
            self.output_control_is_time = False
        else:
            msg = "{0} Error: Unknown value \"{1}\" read for parameter \"{2}\"".format(self.__class__.__name__,output_control_type,"file_label")
            raise Exception(msg)

        self.output_frequency = self.param["output_frequency"].GetDouble()

        # Set Variable list to print
        self.variable_name_list = self.param["gauss_point_results"]
        self.variable_list      = []
        for i in range(self.variable_name_list.size()):
            var_name = self.variable_name_list[i].GetString()
            variable = self._get_variable(var_name)
            self.variable_list.append(variable)

    def ExecuteBeforeSolutionLoop(self): pass

    def ExecuteInitializeSolutionStep(self):
        self.step_count += 1

    def ExecuteFinalizeSolutionStep(self): pass

    def ExecuteFinalize(self): pass

    def PrintOutput(self):
        # Print the output
        self.printed_step_count += 1
        self._get_mp_coords()
        self._get_mp_results()
        self._write_vtk()

        # Schedule next output
        time = self._get_pretty_time(self.model_part.ProcessInfo[KratosMultiphysics.TIME])
        if self.output_frequency > 0.0: # Note: if == 0, we'll just always print
            if self.output_control_is_time:
                while self._get_pretty_time(self.next_output) <= time:
                    self.next_output += self.output_frequency
            else:
                while self.next_output <= self.step_count:
                    self.next_output += self.output_frequency

    def IsOutputStep(self):
        if self.output_control_is_time:
            time = self._get_pretty_time(self.model_part.ProcessInfo[KratosMultiphysics.TIME])
            return (time >= self._get_pretty_time(self.next_output))
        else:
            return ( self.step_count >= self.next_output )

    # Private Functions
    def _get_pretty_time(self,time):
        pretty_time = "{0:.12g}".format(time)
        pretty_time = float(pretty_time)
        return pretty_time

    def _get_attribute(self, my_string, function_pointer, attribute_type):
        """Return the python object named by the string argument.

        To be used with functions from KratosGlobals

        Examples:
        variable = self._get_attribute("DISPLACEMENT",
                                       KratosMultiphysics.ParticleMechanicsApplication.GetVariable,
                                       "Variable")
        """
        splitted = my_string.split(".")

        if len(splitted) == 0:
            raise Exception("Something wrong. Trying to split the string " + my_string)
        if len(splitted) > 3:
            raise Exception("Something wrong. String " + my_string + " has too many arguments")

        attribute_name = splitted[-1]

        if len(splitted) == 2 or len(splitted) == 3:
            warning_msg =  "Ignoring \"" +  my_string.rsplit(".",1)[0]
            warning_msg += "\" for " + attribute_type +" \"" + attribute_name + "\""
            KratosMultiphysics.Logger.PrintInfo("Warning in mpm gid output", warning_msg)

        return function_pointer(attribute_name) # This also checks if the application has been imported

    def _get_variable(self, my_string):
        """Return the python object of a Variable named by the string argument.

        Examples:
        recommended usage:
        variable = self._get_variable("MP_VELOCITY")
        deprecated:
        variable = self._get_variables("KratosMultiphysics.ParticleMechanicsApplication.MP_VELOCITY")
        """
        return self._get_attribute(my_string, KratosMultiphysics.KratosGlobals.GetVariable, "Variable")

    def _get_mp_coords(self):
        number_of_mps = self.model_part.NumberOfElements(0)
        if len(self.coords_X) != number_of_mps:
            self.coords_X = np.empty(number_of_mps)
            self.coords_Y = np.empty(number_of_mps)
            self.coords_Z = np.empty(number_of_mps)

        i = 0
        for mpm in self.model_part.Elements:
            coord = mpm.CalculateOnIntegrationPoints(KratosParticle.MP_COORD,self.model_part.ProcessInfo)[0]
            self.coords_X[i] = coord[0]
            self.coords_Y[i] = coord[1]
            self.coords_Z[i] = coord[2]
            i += 1

    def _get_mp_results(self):
        clock_time = self._start_time_measure()
        number_of_results = self.variable_name_list.size()
        number_of_mps = self.model_part.NumberOfElements(0)

        if len(self.result_names) != number_of_results:
            self.result_names = ["dummy"]*number_of_results
        if len(self.temp_results) != number_of_results:
            self.temp_results = np.empty([number_of_mps,3])

        for result_index in range(number_of_results):
            var_name = self.variable_name_list[result_index].GetString()
            self.result_names[result_index] = var_name

            variable = self.variable_list[result_index]
            var_size = 0
            is_scalar = self._is_scalar(variable)

            # Write in result file
            mpm_index = 0
            for mpm in self.model_part.Elements:
                print_variable = mpm.CalculateOnIntegrationPoints(variable,self.model_part.ProcessInfo)[0]
                if is_scalar:
                    self.temp_results[mpm_index,0] = print_variable
                else:
                    var_size =  print_variable.Size()
                    if var_size == 1 or var_size == 3:
                        for i in range(var_size):
                            self.temp_results[mpm_index,i] = print_variable[i]
                    else:
                        KratosMultiphysics.Logger.PrintInfo("Warning in mpm vtk output", "Printing format is not defined for variable: ", var_name, "with size: ", var_size)

                mpm_index += 1

            # store in dictionary
            if var_size == 1 or is_scalar:
                self.result_dict[var_name] = self._GetCSlice(0)
            elif var_size == 3:
                #self.result_dict[var_name] = (self.temp_results[:,0],self.temp_results[:,1],self.temp_results[:,2])
                self.result_dict[var_name] = (self._GetCSlice(0),self._GetCSlice(1),self._GetCSlice(2))
            else:
                KratosMultiphysics.Logger.PrintInfo("Warning in mpm vtk output", "Printing format is not defined for variable: ", var_name, "with size: ", var_size)

        self._stop_time_measure(clock_time)

    def _GetCSlice(self,index):
        # returns a contiguous slice, needed for vtk writing
        return np.asarray(self.temp_results[:,index], order = 'C')

    def _write_vtk(self):
        particles_filename = self.problem_name + str(self.printed_step_count)
        path = os.path.join(self.vtk_post_path_directory, particles_filename)
        hl.pointsToVTK(path, self.coords_X, self.coords_Y, self.coords_Z, self.result_dict)

    def _start_time_measure(self):
        return time()

    def _stop_time_measure(self, time_ip):
        time_fp = time()
        KratosMultiphysics.Logger.PrintInfo("::[Particle VTK Output Process]:: ", "[Spent time for output = ", time_fp - time_ip, "sec]")

    def _is_scalar(self,variable):
        is_scalar = False
        if (isinstance(variable,KratosMultiphysics.IntegerVariable) or isinstance(variable,KratosMultiphysics.DoubleVariable) or isinstance(variable,KratosMultiphysics.BoolVariable)):
            is_scalar = True
        elif (isinstance(variable,KratosMultiphysics.StringVariable)):
            raise Exception("String variable cant be printed.")
        return is_scalar
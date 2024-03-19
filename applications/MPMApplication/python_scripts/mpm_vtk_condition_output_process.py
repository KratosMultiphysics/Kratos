import numpy as np
import KratosMultiphysics
from KratosMultiphysics.kratos_utilities import IssueDeprecationWarning
import KratosMultiphysics.MPMApplication as KratosMPM
from KratosMultiphysics.MPMApplication.mpm_vtk_output_process import MPMVtkOutputProcess

# Import time library
from time import time

def Factory(settings: KratosMultiphysics.Parameters, model: KratosMultiphysics.Model) -> KratosMultiphysics.OutputProcess:
    if not isinstance(model, KratosMultiphysics.Model):
        raise Exception("expected input shall be a Model object, encapsulating a json string")
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    IssueDeprecationWarning("MPMApplication:","`MPMConditionVTKOutputProcess` is deprecated and replaced with `MPMVtkOutputProcess`")
    return MPMVtkOutputProcess(model, settings["Parameters"])

class LegacyMPMConditionVTKOutputProcess(MPMVtkOutputProcess):
    defaults = KratosMultiphysics.Parameters("""{
        "model_part_name"                    : "MPM_Material",
        "output_control_type"                : "step",
        "output_interval"                   : 1,
        "file_format"                        : "ascii",
        "output_precision"                   : 7,
        "folder_name"                        : "vtk_condition_output",
        "output_sub_model_parts"             : false,
        "save_output_files_in_folder"        : true,
        "gauss_point_results" : []
    }""")

    def __init__(self, model_part, param, Model):

        if param is None:
            param = self.defaults
        else:
            param.ValidateAndAssignDefaults(self.defaults)

        self.model = Model
        self.param = param
        self.model_part = model_part

        # Initiate base class - material point output
        MPMVTKOutputProcess.__init__(self, model_part, param)

    def ExecuteBeforeSolutionLoop(self): 
        if (self.problem_name.startswith('Background_Grid.')):
            self.problem_name = self.problem_name.replace('Background_Grid.','')
        mpm_material_model_part_name = "MPM_Material." + self.problem_name
        self.model_part = self.model[mpm_material_model_part_name]


    def _get_mp_coords(self):
        number_of_mps = self.model_part.NumberOfConditions(0)
        if len(self.coords_X) != number_of_mps:
            self.coords_X = np.empty(number_of_mps)
            self.coords_Y = np.empty(number_of_mps)
            self.coords_Z = np.empty(number_of_mps)

        i = 0
        for mpc in self.model_part.Conditions:
            coord = mpc.CalculateOnIntegrationPoints(KratosMPM.MPC_COORD,self.model_part.ProcessInfo)[0]
            self.coords_X[i] = coord[0]
            self.coords_Y[i] = coord[1]
            self.coords_Z[i] = coord[2]
            i += 1

    def _get_mp_results(self):
        clock_time = self._start_time_measure()
        number_of_results = self.variable_name_list.size()
        number_of_mps = self.model_part.NumberOfConditions(0)

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
            mpc_index = 0
            for mpc in self.model_part.Conditions:
                print_variable = mpc.CalculateOnIntegrationPoints(variable,self.model_part.ProcessInfo)[0]
                if is_scalar:
                    self.temp_results[mpc_index,0] = print_variable
                else:
                    var_size =  print_variable.Size()
                    if var_size == 1 or var_size == 3:
                        for i in range(var_size):
                            self.temp_results[mpc_index,i] = print_variable[i]
                    else:
                        KratosMultiphysics.Logger.PrintInfo("Warning in mpm vtk condition output", "Printing format is not defined for variable: ", var_name, "with size: ", var_size)

                mpc_index += 1

            # store in dictionary
            if var_size == 1 or is_scalar:
                self.result_dict[var_name] = self._GetCSlice(0)
            elif var_size == 3:
                #self.result_dict[var_name] = (self.temp_results[:,0],self.temp_results[:,1],self.temp_results[:,2])
                self.result_dict[var_name] = (self._GetCSlice(0),self._GetCSlice(1),self._GetCSlice(2))
            else:
                KratosMultiphysics.Logger.PrintInfo("Warning in mpm vtk condition output", "Printing format is not defined for variable: ", var_name, "with size: ", var_size)

        self._stop_time_measure(clock_time)

    def _stop_time_measure(self, time_ip):
        time_fp = time()
        KratosMultiphysics.Logger.PrintInfo("::[Material Point Condition VTK Output Process]:: ", "[Spent time for output = ", time_fp - time_ip, "sec]")

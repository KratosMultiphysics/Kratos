# Importing the Kratos Library
import KratosMultiphysics as KM

# other imports
from KratosMultiphysics.time_based_ascii_file_writer_utility import TimeBasedAsciiFileWriterUtility

def Factory(settings, model):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ProcessInfoOutputProcess(model, settings["Parameters"])

class ProcessInfoOutputProcess(KM.OutputProcess):
    """This output process writes variables from ProcessInfo to a file.

    By default it outputs NL_ITERATION_NUMBER, which is used by the ResidualBasedNewtonRaphsonStrategy.

    This process provides the possibility to get the information in each coupling iteration,
    which can be used for coupled simulations by the CoSimulationApplication so far.
    The user may define its execution point.

    This process is not tested for MPI or restarts so far.
    """
    def __init__(self, model, params):
        super().__init__()

        default_settings = KM.Parameters('''{
            "help"                 : "This output process writes variables from ProcessInfo to a file. By default it outputs NL_ITERATION_NUMBER, which is used by the ResidualBasedNewtonRaphsonStrategy. This process provides the possibility to get the information in each coupling iteration, which can be used for coupled simulations by the CoSimulationApplication so far. The user may define its execution point. This process is not tested for MPI or restarts so far.",
            "model_part_name"      : "",
            "interval"             : [0.0, 1e30],
            "execution_point"      : "finalize_solution_step",
            "output_variables"     : ["NL_ITERATION_NUMBER"],
            "print_format"         : "",
            "output_file_settings" : {}
        }''')
        params.ValidateAndAssignDefaults(default_settings)
        self.model = model
        self.interval = KM.IntervalUtility(params)
        self.params = params
        self.format = self.params["print_format"].GetString()

        # validate execution point
        self.execution_point = self.params["execution_point"].GetString()
        list_of_execution_points = ["initialize", "before_solution_loop", "initialize_solution_step", "initialize_coupling_step", "finalize_coupling_step", "finalize_solution_step", "before_output_step", "after_output_step"]
        if not self.execution_point in list_of_execution_points:
            err_msg = 'Execution point "' + self.execution_point + '" is not known!\n'
            err_msg += 'Available execution points are: "' + '", "'.join(list_of_execution_points) + '"'
            raise Exception(err_msg)

        # getting the ModelPart from the Model
        model_part_name = self.params["model_part_name"].GetString()
        if model_part_name == "":
            raise Exception('No "model_part_name" was specified!')
        self.model_part = self.model[model_part_name]

        # retrieving the output variables
        output_var_names = self.params["output_variables"]
        variable_names = [ output_var_names[i].GetString() for i in range( output_var_names.size() ) ]
        self.output_variables = [ KM.KratosGlobals.GetVariable( var ) for var in variable_names ]
        if len(self.output_variables ) == 0:
            raise Exception('No variables specified for output!')
        # initialize bool for checking that variables are in ProcessInfo
        self.variables_of_process_info_unchecked = True

        # create and open file with file header
        file_handler_params = KM.Parameters(self.params["output_file_settings"])
        self.ascii_writer = TimeBasedAsciiFileWriterUtility(self.model_part, file_handler_params, self.__GetFileHeader()).file

        # initialize variable collecting output string
        self.output = ""

    def ExecuteInitialize(self):
        if self.execution_point == "initialize":
            self.__addVariablesToOutput()

    def ExecuteBeforeSolutionLoop(self):
        if self.execution_point == "before_solution_loop":
            self.__addVariablesToOutput()

    def ExecuteInitializeSolutionStep(self):
        if self.execution_point == "initialize_solution_step":
            self.__addVariablesToOutput()

    def ExecuteInitializeCouplingStep(self):
        if self.execution_point == "initialize_coupling_step":
            self.__addVariablesToOutput()

    def ExecuteFinalizeCouplingStep(self):
        if self.execution_point == "finalize_coupling_step":
            self.__addVariablesToOutput()

    def ExecuteFinalizeSolutionStep(self):
        if self.execution_point == "finalize_solution_step":
            self.__addVariablesToOutput()

    def ExecuteBeforeOutputStep(self):
        if self.execution_point == "before_output_step":
            self.__addVariablesToOutput()

    def IsOutputStep(self):
        time = self.model_part.ProcessInfo[KM.TIME]
        if self.interval.IsInInterval(time):
            return True
        else:
            return False

    def PrintOutput(self):
        # write output and reset variable collecting output
        self.ascii_writer.write(self.output)
        self.output = ""

    def ExecuteAfterOutputStep(self):
        if self.execution_point == "after_output_step":
            self.__addVariablesToOutput()

    def ExecuteFinalize(self):
        self.ascii_writer.close()

    def __GetFileHeader(self):
        header = "# ProcessInfo over time"
        if (self.execution_point == "initialize_coupling_step") or (self.execution_point == "finalize_coupling_step"):
            header += " including coupling iterations"
        header += "\n# time"
        for var in self.output_variables:
            header += " " + var.Name()
        header += "\n"
        return header

    def __addVariablesToOutput(self):
        if self.IsOutputStep():
            # check whether variables are in Process Info before they are output for the first time
            if self.variables_of_process_info_unchecked:
                for var in self.output_variables:
                    if not self.model_part.ProcessInfo.Has(var):
                        raise Exception('Given output variable "' + var.Name() + '" does not exist in ProcessInfo!')
                self.variables_of_process_info_unchecked = False
            # add all variables to output with time step
            # if this is called refering to a coupling step, multiple entries in the file will have the same time stemp
            self.output += str(self.model_part.ProcessInfo[KM.TIME])
            for var in self.output_variables:
                if self.variables_of_process_info_unchecked:
                    if not self.model_part.ProcessInfo.Has(var):
                        raise Exception('Given output variable "' + var.Name() + '" does not exist in ProcessInfo!')
                value = str(self.model_part.ProcessInfo[var])
                self.output += " " + format(value,self.format)
            self.output += "\n"

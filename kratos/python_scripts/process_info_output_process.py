# Importing the Kratos Library
import KratosMultiphysics as KM

# other imports
from KratosMultiphysics.time_based_ascii_file_writer_utility import TimeBasedAsciiFileWriterUtility

def Factory(settings, model):
    if(type(settings) != KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ProcessInfoOutputProcess(model, settings["Parameters"])

class ProcessInfoOutputProcess(KM.OutputProcess):
    """This output process writes variables from ProcessInfo to a file.
    By default it outputs NL_ITERATION_NUMBER, which is used by the ResidualBasedNewtonRaphsonStrategy.

    This process provides the possibility to get the information in each coupling iteration,
    which can be used for coupled simulations by the CoSimulationApplication so far.

    This process is not tested for MPI or restarts so far.
    """
    def __init__(self, model, params):
        super().__init__()
        # validate and assign default, create class variables etc. --> see point_output_process

        default_settings = KM.Parameters('''{
            "help"                        : "This output process writes variables from ProcessInfo to a file. By default it outputs NL_ITERATION_NUMBER, which is used by the ResidualBasedNewtonRaphsonStrategy. This process provides the possibility to get the information in each coupling iteration, which can be used for coupled simulations by the CoSimulationApplication so far. This process is not tested for MPI or restarts so far.",
            "model_part_name"             : "",
            "include_coupling_iterations" : false,
            "output_variables"            : ["NL_ITERATION_NUMBER"],
            "print_format"                : "",
            "output_file_settings"        : {}
        }''')
        params.ValidateAndAssignDefaults(default_settings)
        self.model = model
        self.params = params
        self.include_coupling_iterations = self.params["include_coupling_iterations"].GetBool()
        self.format = self.params["print_format"].GetString()

        self.output = ""
        self.latest_coupling_output = ""

    def ExecuteInitialize(self):
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
        # validate types of variables and that variables are in ProcessInfo
        for var in self.output_variables:
            #if not self.model_part.ProcessInfo.Has(var):
                #raise Exception('Given output variable "' + var.Name() + '" does not exist in ProcessInfo!')
            if type(var) == KM.IntegerVariable:
                continue
            elif type(var) == KM.DoubleVariable:
                continue
            else:
                err_msg  = 'Type of variable "' + var.Name() + '" is not valid\n'
                err_msg += 'It can only be integer or double!'
                raise Exception(err_msg)

        # create and open file with file header
        file_handler_params = KM.Parameters(self.params["output_file_settings"])
        self.ascii_writer = TimeBasedAsciiFileWriterUtility(self.model_part, file_handler_params, self.__GetFileHeader()).file

    def ExecuteInitializeSolutionStep(self):
        # empty output strings
        # both variables necessary, so that latest output of coupling iteration is not added to final output
        # in order to omit if condition in ExecuteFinalizeSolutionStep
        self.output = ""
        self.latest_coupling_output = ""

    def ExecuteFinalizeCouplingStep(self):
        if self.include_coupling_iterations:
            self.output += self.latest_coupling_output
            self.latest_coupling_output = str(self.model_part.ProcessInfo[KM.TIME])
            for var in self.output_variables:
                value = self.model_part.ProcessInfo[var]
                self.latest_coupling_output += " " + format(value,self.format)
            self.latest_coupling_output += "\n"

    def ExecuteFinalizeSolutionStep(self):
        # add final output of current time step
        self.output += str(self.model_part.ProcessInfo[KM.TIME])
        for var in self.output_variables:
            value = self.model_part.ProcessInfo[var]
            self.output += " " + format(value,self.format)
        self.output += "\n"

    def IsOutputStep(self):
        return True

    def PrintOutput(self):
        self.ascii_writer.write(self.output)

    def ExecuteFinalize(self):
        self.ascii_writer.close()

    def __GetFileHeader(self):
        header = "# ProcessInfo over time"
        if self.include_coupling_iterations:
            header += " including coupling iterations"
        header += "\n# time"
        for var in self.output_variables:
            header += " " + var.Name()
        header += "\n"

        return header

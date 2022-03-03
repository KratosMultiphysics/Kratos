import KratosMultiphysics as KM

def Factory(settings, model):
    if(type(settings) != KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return CoupledIterationOutputProcess(model, settings["Parameters"])

class CoupledIterationOutputProcess(KM.OutputProcess):
    """describe class

    """
    def __init__(self, model, params):
        super().__init__()
        # validate and assign default, create class variables etc. --> see point_output_process

        default_settings = KM.Parameters('''{
            "help"                   : "This process ...",
            "print_iteration_number" : true,
            "model_part_name"        : "",
            "interval"               : [0.0, 1e30],
            "output_variables"       : [],
            "output_file_settings"   : {}
        }''')
        params.ValidateAndAssignDefaults(default_settings)

        self.model = model
        self.interval = KM.IntervalUtility(params)
        self.params = params

    def ExecuteInitialize(self):
        # get name of solver
        # depending on solver get correct variable in which iteration number is stored
        # GaussSeidelStrongCoupledSolver - KratosCoSim.COUPLING_ITERATION_NUMBER (co-simulation)
        # ResidualBasedNewtonRaphsonStrategy - r_model_part.GetProcessInfo()[NL_ITERATION_NUMBER] (fluid and structure of Mok test case)

        # create/ open file with file header --> see point_output_process
        pass

    def ExecuteBeforeSolutionLoop(self):
        pass

    def ExecuteInitializeSolutionStep(self):
        pass

    def ExecuteFinalizeCouplingStep(self):
        # get iteration number or displacement of a point, etc.
        # get value from corresponding variable

        # ResidualBasedNewtonRaphsonStrategy if ...
        # if fluid structural mechanics
            # if self.settings["analysis_type"].GetString() == "non_linear"
            # and:
            # solver_type = solver_settings["solver_type"].GetString()
            # if solver_type == "monolithic" or solver_type == "Monolithic"
            # elif solver_type == "monolithic_stokes" or solver_type == "MonolithicStokes":
            # elif (solver_type == "Embedded"):
            # elif (solver_type == "Compressible"):
            # elif (solver_type == "CompressibleExplicit"):
            # elif (solver_type == "ConjugateHeatTransfer"):
            # elif solver_type == "two_fluids" or solver_type == "TwoFluids":
            # INFO: fractional_step doesn't save iteration number to any variable, it's only printed
        # if structural mechanics
            # if self.settings["analysis_type"].GetString() == "non_linear"
            # and (self.settings["line_search"].GetBool() == False):
            # if solver_settings.Has("time_integration_method"):
            #     time_integration_method = solver_settings["time_integration_method"].GetString()
            # else:
            #     time_integration_method = "implicit" # defaulting to implicit time-integration
            # and if time_integration_method == "implicit"

        # try exept statement or check if empty etc., check what happens if variable is not used in solving strategy
        # --> use r_model_part.GetProcessInfo()[NL_ITERATION_NUMBER]

        model_part_name = self.params["model_part_name"].GetString()
        if model_part_name == "":
            raise Exception('No "model_part_name" was specified!')
        self.model_part = self.model[model_part_name]
        num_iter = self.model_part.ProcessInfo[KM.NL_ITERATION_NUMBER]

        print("[TEST ITERATION OUTPUT] Number of iterations for " + model_part_name + ": " + str(num_iter))

    def ExecuteFinalizeSolutionStep(self):
        pass

    def ExecuteBeforeOutputStep(self):
        pass

    def ExecuteAfterOutputStep(self):
        pass

    def ExecuteFinalize(self):
        # close file --> see vtk_output_process and point_output_process
        pass

    def IsOutputStep(self):
        # should return "True" when inside the given interval --> see vtk_output_process
        return True

    def PrintOutput(self):
        # will be called once per time step, so data of all coupling iterations during the time step
        # should be collected and then print to the output file here
        pass

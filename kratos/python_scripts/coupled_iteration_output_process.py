import KratosMultiphysics as KM

def Factory(settings, model):
    if(type(settings) != KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return CoupledIterationOutputProcess(model, settings["Parameters"])

class CoupledIterationOutputProcess(KM.OutputProcess):
    """

    """
    def __init__(self, model, params):
        super().__init__()

        default_settings = KM.Parameters('''{
            "help"                 : "This process ...
            "output_variables"     : [],
            "output_file_settings" : {}
        }''')
        params.ValidateAndAssignDefaults(default_settings)

        # get name of solver
        # depending on solver get correct variable in which iteration number is stored
        # GaussSeidelStrongCoupledSolver - KratosCoSim.COUPLING_ITERATION_NUMBER (co-simulation)
        # ResidualBasedNewtonRaphsonStrategy - r_model_part.GetProcessInfo()[NL_ITERATION_NUMBER] (fluid and structure of Mok test case)

    def ExecuteInitialize(self):
        # create/ open file
        pass

    def ExecuteBeforeSolutionLoop(self):
        pass

    def ExecuteInitializeSolutionStep(self):
        pass

    def ExecuteFinalizeCouplingStep(self):
        # get iteration number or displacement of a point, etc.
        # get value from corresponding variable
        print("[TEST ITERATION OUTPUT] another coupling iteration done!")

    def ExecuteFinalizeSolutionStep(self):
        pass

    def ExecuteBeforeOutputStep(self):
        pass

    def ExecuteAfterOutputStep(self):
        pass

    def ExecuteFinalize(self):
        # close file
        pass

    def IsOutputStep(self):
        # True or intervall --> vtk output
        return True

    def PrintOutput(self):
        # print to file --> vtk output
        pass

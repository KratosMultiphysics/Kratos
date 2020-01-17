from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.MeshingApplication as KratosMesh
import KratosMultiphysics.ConvectionDiffusionApplication as KratosConvDiff
import KratosMultiphysics.FluidDynamicsApplication as KratosFluid

def Factory(settings, Model): 
    if(type(settings) != KratosMultiphysics.Parameters): 
        raise Exception("expected input shall be a Parameters object, encapsulating a json string") 
    return ApplyErrorCalculatorProcess(Model, settings["Parameters"])

## All the processes python should be derived from "Process" 
class ApplyErrorCalculatorProcess(KratosMultiphysics.Process):

    def __init__(self, Model, settings):
        KratosMultiphysics.Process.__init__(self)

        default_settings = KratosMultiphysics.Parameters(""" 
            {
                "model_part_name"                       : "MainModelPart",
                "minimal_size"                        : 0.005,
                "maximal_size"                        : 2.0,
                "nodal_averaging"                     : true,
                "reference_variable_name"             : "ERROR_RATIO",
                "historical_results"                  : true,
                "echo_level"                          : 0
            }
            """
        )
        settings.ValidateAndAssignDefaults(default_settings)

        self.model_part = Model[settings["model_part_name"].GetString()]
        self.historical_results = settings["historical_results"].GetBool()
        print(self.model_part.NumberOfNodes())
        print(self.model_part.NumberOfElements())


        #Process Parameters
        error_calc_params = KratosMultiphysics.Parameters("""{}""")
        error_calc_params.AddValue("minimal_size", settings["minimal_size"])
        error_calc_params.AddValue("maximal_size", settings["maximal_size"])
        error_calc_params.AddValue("nodal_averaging", settings["nodal_averaging"])
        error_calc_params.AddValue("reference_variable_name", settings["reference_variable_name"])
        error_calc_params.AddValue("echo_level", settings["echo_level"])
        error_calc_params.AddValue("historical_results", settings["historical_results"])

        if (self.model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 2):
            self.error_calc = KratosConvDiff.MetricsTemperatureGradientProcess2D(self.model_part, error_calc_params)
        else:
            self.error_calc = KratosConvDiff.MetricsTemperatureGradientProcess3D(self.model_part, error_calc_params)
        
        remesh_params = KratosMultiphysics.Parameters("""{}""")
        if (self.model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 2):
            self.remesh_process = KratosMesh.MmgProcess2D(self.model_part, remesh_params)
        else:
            self.remesh_process = KratosMesh.MmgProcess3D(self.model_part, remesh_params)
    
    def ExecuteInitialize(self):
        pass

    def ExecuteBeforeSolutionLoop(self):
        pass

    def ExecuteInitializeSolutionStep(self):
        pass

    def ExecuteFinalizeSolutionStep(self):
        if self.historical_results:
            # We execute our process
            self.error_calc.Execute()
        # if self.model_part.ProcessInfo[KratosMultiphysics.TIME] == 10.0:
        #     self.remesh_process.Execute()

    def ExecuteBeforeOutputStep(self):
        pass

    def ExecuteAfterOutputStep(self):
        pass

    def ExecuteFinalize(self):
        if not self.historical_results:
            # We execute our process
            self.error_calc.Execute()
        self.remesh_process.Execute()
        print(self.model_part)
        KratosMultiphysics.ModelPartIO("Revised_fluid_buoyancy", KratosMultiphysics.IO.WRITE).WriteModelPart(self.model_part)
        

    def Clear(self):
        pass




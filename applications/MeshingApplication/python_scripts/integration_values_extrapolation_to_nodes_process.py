from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Check that applications were imported in the main script
KratosMultiphysics.CheckRegisteredApplications("MeshingApplication")

# Import applications
import KratosMultiphysics.MeshingApplication as MeshingApplication

def Factory(settings, Model): 
    if(type(settings) != KratosMultiphysics.Parameters): 
        raise Exception("expected input shall be a Parameters object, encapsulating a json string") 
    return IntegrationValuesExtrapolationToNodesProcess(Model, settings["Parameters"])
 
## All the processes python should be derived from "Process" 
class IntegrationValuesExtrapolationToNodesProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings ): 
        KratosMultiphysics.Process.__init__(self) 
 
        default_settings = KratosMultiphysics.Parameters(""" 
        {
            "help"                       : "This process extrapolates the values from integration points to the mesh nodes",
            "model_part_name"            : "",
            "echo_level"                 : 0,
            "average_variable"           : "NODAL_AREA",
            "area_average"               : true,
            "list_of_variables"          : [],
            "extrapolate_non_historical" : true
        }
        """
        )

        settings.ValidateAndAssignDefaults(default_settings)

        self.model_part = Model[settings["model_part_name"].GetString()]

        extrapolation_parameters = KratosMultiphysics.Parameters("""{}""")
        extrapolation_parameters.AddValue("echo_level", settings["echo_level"])
        extrapolation_parameters.AddValue("average_variable", settings["average_variable"])
        extrapolation_parameters.AddValue("area_average", settings["area_average"])
        extrapolation_parameters.AddValue("list_of_variables", settings["list_of_variables"])
        extrapolation_parameters.AddValue("extrapolate_non_historical", settings["extrapolate_non_historical"])
        self.integration_values_extrapolation_to_nodes_process = MeshingApplication.IntegrationValuesExtrapolationToNodesProcess(self.model_part, extrapolation_parameters)

    def ExecuteInitialize(self):
        pass

    def ExecuteBeforeSolutionLoop(self):
        self.integration_values_extrapolation_to_nodes_process.ExecuteBeforeSolutionLoop()

    def ExecuteInitializeSolutionStep(self):
        pass

    def ExecuteFinalizeSolutionStep(self):
        self.integration_values_extrapolation_to_nodes_process.ExecuteFinalizeSolutionStep()

    def ExecuteBeforeOutputStep(self):
        pass

    def ExecuteAfterOutputStep(self):
        pass

    def ExecuteFinalize(self):
        self.integration_values_extrapolation_to_nodes_process.ExecuteFinalize()

    def Clear(self):
        pass

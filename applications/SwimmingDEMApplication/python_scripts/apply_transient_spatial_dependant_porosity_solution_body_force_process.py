# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.SwimmingDEMApplication as KratosSDEM

from importlib import import_module

def Factory(settings, Model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyTransientSpatialDependantPorositySolutionBodyForceProcess(Model, settings["Parameters"])

def CreateManufacturedSolution(custom_settings):
    return ApplyTransientSpatialDependantPorositySolutionBodyForceProcess(model_part, custom_settings)

## All the processes python should be derived from "Process"
class ApplyTransientSpatialDependantPorositySolutionBodyForceProcess(KratosMultiphysics.Process):
    def __init__(self, model, settings ):

        KratosMultiphysics.Process.__init__(self)

        default_settings = KratosMultiphysics.Parameters("""
            {
                "model_part_name"          : "please_specify_model_part_name",
                "variable_name"            : "BODY_FORCE",
                "benchmark_name"           : "custom_body_force.vortex",
                "benchmark_parameters"     : {},
                "compute_nodal_error"      : true,
                "print_convergence_output" : false,
                "output_parameters"        : {}
            }
            """
            )

        self.settings = settings
        self.settings.ValidateAndAssignDefaults(default_settings)

        self.model_part = model[self.settings["model_part_name"].GetString()]
        self.variable = KratosMultiphysics.KratosGlobals.GetVariable(self.settings["variable_name"].GetString())

        self.TransientSpatialDependantPorositySolutionBodyForceProcess = KratosSDEM.TransientSpatialDependantPorositySolutionBodyForceProcess(self.model_part, self.settings)

    def ExecuteBeforeSolutionLoop(self):
        self.TransientSpatialDependantPorositySolutionBodyForceProcess.ExecuteBeforeSolutionLoop()


    def ExecuteInitializeSolutionStep(self):
        self.TransientSpatialDependantPorositySolutionBodyForceProcess.ExecuteInitializeSolutionStep()

    def ExecuteFinalizeSolutionStep(self):
        pass
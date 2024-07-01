import KratosMultiphysics
import KratosMultiphysics.ConvectionDiffusionApplication as KratosConvDiff

import numpy as np

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyMicrofluidicTubeProcess(Model, settings["Parameters"])


class ApplyMicrofluidicTubeProcess(KratosMultiphysics.Process):
    def __init__(self, model, settings):
        KratosMultiphysics.Process.__init__(self)

        # Compare with default settings
        default_settings = KratosMultiphysics.Parameters(r'''{
            "model_part_name"              : "please_specify_model_part_name",
            "outlet_model_part_name"       : "outlet_model_part_name",
            "geometry_type"                : "circular",
            "pressure_outlet"              : 0.0,
            "dimensions"                   : [1.0, 1.0],
            "materials_filename"           : "json_file",
            "diffusion_coefficient"        : 1.0
        }''')

        self.settings = settings
        self.settings.ValidateAndAssignDefaults(default_settings)

        self.model = model
        self.model_part = model[self.settings["model_part_name"].GetString()]

        # Check input data
        if (settings["model_part_name"].GetString() == ""):
            raise Exception("Empty model part name string. Set a valid model part name.")

        

        self.ApplyMicrofluidicTubeProcess_ = KratosConvDiff.MicrofluidicTubeProcess(self.model_part, self.settings)

        # Nodal area
        self.area_calculator = KratosMultiphysics.CalculateNodalAreaProcess(self.model_part, 3)
        self.area_calculator.Execute()

    def check_geometry(self, geometry):
        if geometry not in ("circular", "square"):
            raise ValueError(f"MicrofluidicTubeProcess Error: Invalid geometry '{geometry}'. Valid geometries are: 'circular', 'square'.")

    def ExecuteBeforeSolutionLoop(self):
        self.ApplyMicrofluidicTubeProcess_.ExecuteBeforeSolutionLoop()

    def ExecuteInitializeSolutionStep(self):
        self.ApplyMicrofluidicTubeProcess_.ExecuteInitializeSolutionStep()
        # for node in self.model_part.Nodes:
        #     if node.Z >= 0.:
        #         node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, 1.)
        #     else:
        #         node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, 0.)

    def ExecuteFinalizeSolutionStep(self):
        self.ApplyMicrofluidicTubeProcess_.ExecuteFinalizeSolutionStep()

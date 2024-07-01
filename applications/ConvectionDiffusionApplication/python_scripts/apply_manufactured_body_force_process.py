import KratosMultiphysics
import KratosMultiphysics.ConvectionDiffusionApplication as KratosConvDiff

import numpy as np

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyManufacturedBodyForceProcess(Model, settings["Parameters"])


class ApplyManufacturedBodyForceProcess(KratosMultiphysics.Process):
    def __init__(self, model, settings):
        KratosMultiphysics.Process.__init__(self)

        # Compare with default settings
        default_settings = KratosMultiphysics.Parameters(r'''{
            "model_part_name"              : "please_specify_model_part_name",
            "impose_inlet_conditions"      : true,
            "write_exact_solution_output"  : true,
            "cavity"                       : false,
            "inlet_model_part_name"        : "inlet_model_part_name",
            "outlet_model_part_name"        : "outlet_model_part_name",
            "central_velocity"             : 1.0,
            "radius"                       : 1.0,
            "length"                       : 1.0,
            "materials_filename"           : "materials_filename.json",
            "add_heat_flux_source"         : false,
            "add_body_force_source"        : true,
            "temperature_gradient"         : 1.0,
            "temperature_inlet"            : 1.0
        }''')

        self.settings = settings
        self.settings.ValidateAndAssignDefaults(default_settings)

        self.model = model
        self.model_part = model[self.settings["model_part_name"].GetString()]

        # Check input data
        if (settings["model_part_name"].GetString() == ""):
            raise Exception("Empty model part name string. Set a valid model part name.")

        self.ApplyManufacturedBodyForceProcess_ = KratosConvDiff.ManufacturedBodyForceProcess(self.model_part, self.settings)

        # Nodal area
        self.area_calculator = KratosMultiphysics.CalculateNodalAreaProcess(self.model_part, 3)
        self.area_calculator.Execute()

    def ExecuteBeforeSolutionLoop(self):
        self.ApplyManufacturedBodyForceProcess_.ExecuteBeforeSolutionLoop()

    def ExecuteInitializeSolutionStep(self):
        self.ApplyManufacturedBodyForceProcess_.ExecuteInitializeSolutionStep()

    def ExecuteFinalizeSolutionStep(self):
        self.ApplyManufacturedBodyForceProcess_.ExecuteFinalizeSolutionStep()

        # Flux
        # from KratosMultiphysics.ConvectionDiffusionApplication import ComputeFlux
        # print(f"\n{'#' * 10} FLUX {'#' * 10}")
        # inlet_model_part = self.model.GetModelPart("FluidModelPart.ImposedTemperature3D_Automatic_inlet_velocity_Auto1")
        # outlet_model_part = self.model.GetModelPart("FluidModelPart.Outlet3D_Outlet_pressure_Auto1")

        # inlet_flux = abs(ComputeFlux(inlet_model_part, KratosMultiphysics.VELOCITY).ComputeSurfaceIntegral())
        # outlet_flux = abs(ComputeFlux(outlet_model_part, KratosMultiphysics.VELOCITY).ComputeSurfaceIntegral())

        # exact_flux = .5 * np.pi * pow(.5e-3, 2) * 5e-5

        # print(f"Inlet:  {inlet_flux} m^3 / s")
        # print(f"Outlet: {outlet_flux} m^3 / s\n")
        # print(f"Exact: {exact_flux} m^3 / s\n")
        # print(f"Relative error flux: {abs(inlet_flux - outlet_flux) / inlet_flux}")
        # print(f"{'#' * 10}######{'#' * 10}")

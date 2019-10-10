from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Import applications

# Importing the base class
from KratosMultiphysics.ConvectionDiffusionApplication.coupled_fluid_thermal_solver import CoupledFluidThermalSolver
from KratosMultiphysics.python_solver import PythonSolver

from importlib import import_module

def CreateSolver(main_model_part, custom_settings):
    return CoupledFluidTransportSolver(main_model_part, custom_settings)

class CoupledFluidTransportSolver(CoupledFluidThermalSolver):

    def __init__(self, model, custom_settings):

        self._validate_settings_in_baseclass=True # To be removed eventually

        PythonSolver.__init__(self, model, custom_settings)

        ## Get domain size
        self.domain_size = self.settings["domain_size"].GetInt()

        from KratosMultiphysics.FluidDynamicsApplication import python_solvers_wrapper_fluid
        self.fluid_solver = python_solvers_wrapper_fluid.CreateSolverByParameters(self.model, self.settings["fluid_solver_settings"],"OpenMP")

        python_module_name = "KratosMultiphysics.FluidTransportApplication"
        full_module_name = python_module_name + "." + self.settings["thermal_solver_settings"]["solver_type"].GetString()
        solver_module = import_module(full_module_name)
        self.thermal_solver = solver_module.CreateSolver(self.model, self.settings["thermal_solver_settings"])

    @classmethod
    def GetDefaultSettings(cls):

        # TODO adapt thermal solver settings to FluidTransportReplaceSolver settings
        this_defaults = KratosMultiphysics.Parameters("""
        {
            "solver_type" : "coupled_fluid_transport_solver",
            "domain_size" : -1,
            "echo_level": 0,
            "fluid_solver_settings": {
                "solver_type": "navier_stokes_solver_vmsmonolithic",
                "model_import_settings": {
                    "input_type": "mdpa",
                    "input_filename": "unknown_name"
                }
            },
            "thermal_solver_settings": {
                "solver_type": "fluid_transport_replace_solver",
                "solution_type":    "Steady",
                "scheme_type":      "Implicit",
                "newmark_theta":     0.5,
                "strategy_type":   "Linear",
                "time_step": 0.1,
                "model_import_settings": {
                    "input_type": "use_input_model_part"
                },
                "material_import_settings": {
                        "materials_filename": "ThermalMaterials.json"
                }
            }
        }
        """)

        this_defaults.AddMissingParameters(PythonSolver.GetDefaultSettings())
        return this_defaults

    # def ImportModelPart(self):
    #     # Call the fluid solver to import the model part from the mdpa
    #     self.fluid_solver.ImportModelPart()

    #     # Save the convection diffusion settings
    #     convection_diffusion_settings = self.thermal_solver.main_model_part.ProcessInfo.GetValue(KratosMultiphysics.CONVECTION_DIFFUSION_SETTINGS)

    #     # Here the fluid model part is cloned to be thermal model part so that the nodes are shared
    #     modeler = KratosMultiphysics.ConnectivityPreserveModeler()
    #     if self.domain_size == 2:
    #         modeler.GenerateModelPart(self.fluid_solver.main_model_part,
    #                                   self.thermal_solver.main_model_part,
    #                                   "EulerianConvDiff2D",
    #                                   "ThermalFace2D2N")
    #     else:
    #         modeler.GenerateModelPart(self.fluid_solver.main_model_part,
    #                                   self.thermal_solver.main_model_part,
    #                                   "EulerianConvDiff3D",
    #                                   "ThermalFace3D3N")

    #     # Set the saved convection diffusion settings to the new thermal model part
    #     self.thermal_solver.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.CONVECTION_DIFFUSION_SETTINGS, convection_diffusion_settings)

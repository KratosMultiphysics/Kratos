# Importing the Kratos Library
import KratosMultiphysics as KM

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_coupling_operation import CoSimulationCouplingOperation

import os

def Create(*args):
    return CouplingOutput(*args)

class CouplingOutput(CoSimulationCouplingOperation):
    """This operation computes the Normals (NORMAL) on a given ModelPart
    """
    def __init__(self, settings, solver_wrappers, process_info):
        super().__init__(settings, process_info)
        self.model = solver_wrappers[self.settings["solver"].GetString()].model
        self.execution_point = self.settings["execution_point"].GetString()

        available_execution_points = [
            "initialize_solution_step",
            "finalize_solution_step",
            "initialize_coupling_iteration",
            "finalize_coupling_iteration",
        ]

        if self.execution_point not in available_execution_points:
            raise Exception("...")

        self.step = 0 # this should come from self.process_info

    def InitializeSolutionStep(self):
        self.step += 1
        self.coupling_iteration = 0
        if self.execution_point == "initialize_solution_step":
            pass

    def FinalizeSolutionStep(self):
        if self.execution_point == "finalize_solution_step":
            pass

    def InitializeCouplingIteration(self):
        self.coupling_iteration += 1
        if self.execution_point == "initialize_coupling_iteration":
            pass

    def FinalizeCouplingIteration(self):
        if self.execution_point == "finalize_coupling_iteration":
            model_part_name = self.settings["output_parameters"]["model_part_name"].GetString()
            KM.VtkOutput(self.model[model_part_name], self.settings["output_parameters"]).PrintOutput(self.settings["solver"].GetString()+"_"+model_part_name+"_"+str(self.coupling_iteration)) # this also validates the settings


    def PrintInfo(self):
        pass


    @classmethod
    def _GetDefaultParameters(cls):
        this_defaults = KM.Parameters("""{
            "solver"            : "UNSPECIFIED",
            "execution_point"   : "UNSPECIFIED",
            "output_format"     : "vtk",
            "output_parameters" : { }
        }""")
        this_defaults.AddMissingParameters(super()._GetDefaultParameters())
        return this_defaults

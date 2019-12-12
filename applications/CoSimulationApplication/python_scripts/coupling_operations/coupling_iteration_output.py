from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics as KM

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_coupling_operation import CoSimulationCouplingOperation

# CoSimulation imports
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools

def Create(settings, solver_wrappers):
    cs_tools.SettingsTypeCheck(settings)
    return CouplingInterationOutput(settings, solver_wrappers)

class CouplingInterationOutput(CoSimulationCouplingOperation):
    """This operation outputs during the coupling iterations
    """
    def __init__(self, settings, solver_wrappers):
        super(CouplingInterationOutput, self).__init__(settings)
        self.solver_model = solver_wrappers[self.settings["solver"].GetString()].model

    def Initialize(self):
        output_format = self.settings["output_parameters"].GetString()
        if output_format == "vtk":
            from KratosMultiphysics.vtk_output_process import Factory
        elif output_format == "gid":
            from KratosMultiphysics.gid_output_process import Factory
        else:
            raise Exception('Currently only "vtk" and "gid" are supported')

        self.output_process = Factory(self.solver_model, self.settings["output_parameters"])

        self.output_process.ExecuteInitialize()
        self.output_process.ExecuteBeforeSolutionLoop()

    def InitializeSolutionStep(self):
        self.coupling_iteration = 0
        self.output_process.ExecuteInitializeSolutionStep()

    def FinalizeCouplingIteration(self):
        self.coupling_iteration += 1

        # setting variables needed inside the output for setting the correct labels
        process_info = self.output_process.model_part.ProcessInfo
        orig_time = process_info[KM.TIME]
        orig_step = process_info[KM.STEP]
        process_info[KM.TIME] = self.coupling_iteration
        process_info[KM.STEP] = self.coupling_iteration

        self.output_process.PrintOutput()

        # resetting internal variables
        process_info[KM.TIME] = orig_time
        process_info[KM.STEP] = orig_step

    def FinalizeSolutionStep(self):
        self.output_process.ExecuteFinalizeSolutionStep()

    def Finalize(self):
        self.output_process.ExecuteFinalize()

    def PrintInfo(self):
        pass

    def Check(self):
        pass

    @classmethod
    def _GetDefaultSettings(cls):
        this_defaults = KM.Parameters("""{
            "solver"           : "UNSPECIFIED",
            "output_format"    : "vtk",
            "output_parameters : { }
        }""")
        this_defaults.AddMissingParameters(super(CouplingInterationOutput, cls)._GetDefaultSettings())
        return this_defaults




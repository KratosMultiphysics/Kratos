from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics as KM

# CoSimulation imports
from KratosMultiphysics.CoSimulationApplication.factories.convergence_accelerator_factory import CreateConvergenceAccelerator

class ConvergenceAcceleratorWrapper(object):
    """This class wraps the convergence accelerators such that they can be used "automized"
    => this class stores the residual and updates the solutions, such that the
    convergence accelerator can be configured through json
    """
    def __init__(self, settings, solver_wrapper):
        self.interface_data = solver_wrapper.GetInterfaceData(settings["data_name"].GetString())
        settings.RemoveValue("data_name")
        settings.RemoveValue("solver")

        self.conv_acc = CreateConvergenceAccelerator(settings)

    def Initialize(self):
        self.conv_acc.Initialize()

    def Finalize(self):
        self.conv_acc.Finalize()

    def InitializeSolutionStep(self):
        self.conv_acc.InitializeSolutionStep()

    def FinalizeSolutionStep(self):
        self.conv_acc.FinalizeSolutionStep()

    def InitializeNonLinearIteration(self):
        # Saving the previous data for the computation of the residual
        # and the computation of the solution update
        self.input_data = self.interface_data.GetData()

        self.conv_acc.InitializeNonLinearIteration()

    def FinalizeNonLinearIteration(self):
        self.conv_acc.FinalizeNonLinearIteration()

    def ComputeAndApplyUpdate(self):
        current_data = self.interface_data.GetData()
        residual = current_data - self.input_data
        updated_data = self.input_data + self.conv_acc.UpdateSolution(residual, self.input_data)
        self.interface_data.SetData(updated_data)

    def PrintInfo(self):
        self.conv_acc.PrintInfo()

    def Check(self):
        self.conv_acc.Check()

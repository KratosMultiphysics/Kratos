from KratosMultiphysics import Parameters
import os
# Importing the Kratos Library
import KratosMultiphysics
file_path = os.path.abspath(__file__)
dir_path = os.path.dirname(file_path)
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.SwimmingDEMApplication.swimming_DEM_analysis import SwimmingDEMAnalysis

class DragTestAnalysis(SwimmingDEMAnalysis, KratosUnittest.TestCase):
    def __init__(self, model, varying_parameters = Parameters("{}")):
        super().__init__(model, varying_parameters)

    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()

    def CheckValues(self, x_vel):
        tol = 1.0e-18
        x_vel_ref = 0.9156414109471193 #ChienDragLaw

        # Extra results.
        # x_vel_ref = 0.9102456241767058 #SchillerAndNaumannDragLaw
        # x_vel_ref = 0.9929471943824714 #StokesDragLaw
        # x_vel_ref = 0.848557847912186  #BeetstraDragLaw
        # x_vel_ref = 0.9050175862530042 #HaiderAndLevenspielDragLaw
        # x_vel_ref = 0.7623255437033674 #GanserDragLaw

        self.assertAlmostEqual(x_vel, x_vel_ref, delta=tol)


    def Finalize(self):
        for node in self.spheres_model_part.Nodes:
            if node.Id == 111:
                x_vel = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_X)

        self.CheckValues(x_vel)
        self.procedures.RemoveFoldersWithResults(str(dir_path), str('chien_law/chien_drag_test'), '')
        super().Finalize()
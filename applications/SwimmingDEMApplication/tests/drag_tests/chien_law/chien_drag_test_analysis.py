from KratosMultiphysics import Parameters
import os
# Importing the Kratos Library
import KratosMultiphysics
file_path = os.path.abspath(__file__)
dir_path = os.path.dirname(file_path)
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.SwimmingDEMApplication.swimming_DEM_analysis import SwimmingDEMAnalysis

class ChienDragAnalysis(SwimmingDEMAnalysis, KratosUnittest.TestCase):
    def __init__(self, model, varying_parameters = Parameters("{}")):
        super().__init__(model, varying_parameters)

    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()

    def CheckValues(self, x_vel):
        tol = 1.0e-18
        x_vel_ref = 0.9886575480896711 #ChienDragLaw

        # Other results.
        # x_vel_ref = 0.9880047941854828 #SchillerAndNaumannDragLaw
        #StokesDragLaw BeetstraDragLaw HaiderAndLevenspielDragLaw GanserDragLaw

        self.assertAlmostEqual(x_vel, x_vel_ref, delta=tol)


    def Finalize(self):
        for node in self.spheres_model_part.Nodes:
            if node.Id == 111:
                x_vel = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_X)
                print(x_vel)

        self.CheckValues(x_vel)
        self.procedures.RemoveFoldersWithResults(str(dir_path), str('chien_drag_test'), '')
        super().Finalize()
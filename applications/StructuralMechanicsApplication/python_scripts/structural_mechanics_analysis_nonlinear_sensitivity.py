from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing Kratos
import KratosMultiphysics
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication

# Importing the base class
from structural_mechanics_analysis import StructuralMechanicsAnalysis

import numpy as np

class StructuralMechanicsAnalysisNLSensitivity(StructuralMechanicsAnalysis):
    """
    This class is the special-script of the StructuralMechanicsApplication put in a class

    It is used to specifiy the non-linear (over- resp. under-linear ) behaviour of state results (e.g. displacements or stresses).
    """
    def __init__(self, model, project_parameters):
        #solver_settings = project_parameters["solver_settings"]

        super(StructuralMechanicsAnalysisNLSensitivity, self).__init__(model, project_parameters)

    def Initialize(self):
        ##here we initialize user-provided processes
        for process in self._GetListOfProcesses():
            process.ExecuteInitialize()

        for process in self._GetListOfProcesses():
            process.ExecuteBeforeSolutionLoop()

        ## Stepping and time settings
        self.end_time = self.project_parameters["problem_data"]["end_time"].GetDouble()

        if self._GetSolver().GetComputingModelPart().ProcessInfo[KratosMultiphysics.IS_RESTARTED]:
            self.time = self._GetSolver().GetComputingModelPart().ProcessInfo[KratosMultiphysics.TIME]
            self.start_time = self._GetSolver().GetComputingModelPart().ProcessInfo[KratosMultiphysics.TIME]
        else:
            self.time = self.project_parameters["problem_data"]["start_time"].GetDouble()

        if self.is_printing_rank:
            KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(), "Analysis -START- ")

    def FinalizeSolutionStep(self):
        # Check if we are passing already in the last pseudo time step
        if (self.time < self.end_time) == False:
            # Compute displacement E(F)-sensitivity
            model_part = self._GetSolver().GetComputingModelPart()
            time_step = self._GetSolver().ComputeDeltaTime()
            load_factors = np.array([self.start_time, self.start_time + time_step, self.start_time + time_step*2])
            self._ComputeEFDisplacementCurvature(model_part, load_factors)
            self._ComputeFirstAndSecondOrderNLSensitivityFactors(model_part, load_factors)

        super(StructuralMechanicsAnalysisNLSensitivity, self).FinalizeSolutionStep()


    def _ComputeEFDisplacementCurvature(self, model_part, load_factor):
        for node in model_part.Nodes:
            disp_1_x = abs( node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X, 2) )
            disp_2_x = abs( node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X, 1) )
            disp_3_x = abs( node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X, 0) )
            disp_1_y = abs( node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y, 2) )
            disp_2_y = abs( node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y, 1) )
            disp_3_y = abs( node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y, 0) )
            disp_1_z = abs( node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Z, 2) )
            disp_2_z = abs( node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Z, 1) )
            disp_3_z = abs( node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Z, 0) )

            disp_x = np.array([disp_1_x, disp_2_x, disp_3_x])
            px = np.polyfit(load_factor, disp_x, 2)
            disp_y = np.array([disp_1_y, disp_2_y, disp_3_y])
            py = np.polyfit(load_factor, disp_y, 2)
            disp_z = np.array([disp_1_z, disp_2_z, disp_3_z])
            pz = np.polyfit(load_factor, disp_z, 2)

            node.SetValue(StructuralMechanicsApplication.DISPLACEMENT_NL_SENSITIVITY_X, 2 * px[0])
            node.SetValue(StructuralMechanicsApplication.DISPLACEMENT_NL_SENSITIVITY_Y, 2 * py[0])
            node.SetValue(StructuralMechanicsApplication.DISPLACEMENT_NL_SENSITIVITY_Z, 2 * pz[0])

    def _ComputeFirstAndSecondOrderNLSensitivityFactors(self, model_part, load_factor):
        lambda_0 = load_factor[0]
        lambda_1 = load_factor[1]
        lambda_2 = load_factor[2]
        f_1 = lambda_1 / lambda_0
        #f_2 = lambda_2 / lambda_0
        delta_10 = lambda_1 - lambda_0
        delta_20 = lambda_2 - lambda_0
        for node in model_part.Nodes:
            disp_1_x = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X, 2)
            disp_2_x = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X, 1)
            disp_3_x = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X, 0)
            disp_1_y = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y, 2)
            disp_2_y = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y, 1)
            disp_3_y = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y, 0)
            disp_1_z = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Z, 2)
            disp_2_z = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Z, 1)
            disp_3_z = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Z, 0)
            #print("Displacement = ", disp_1_x)
            #print("Displacement = ", disp_2_x)
            #print("Displacement = ", disp_3_x)

            if abs(disp_1_x) > 1e-8:
                sensitivity_first_order_1_x = disp_2_x / ( disp_1_x * f_1 )
                #sensitivity_first_order_2_x = disp_3_x / ( disp_1_x * f_2 )
                #sensitivity_second_order_x = sensitivity_first_order_2_x / sensitivity_first_order_1_x
                slope_10_x = (disp_2_x-disp_1_x)/delta_10
                slope_20_x = (disp_3_x-disp_1_x)/delta_20
                sensitivity_second_order_x = slope_20_x / slope_10_x
                node.SetValue(StructuralMechanicsApplication.NL_SENSITIVITY_FIRST_ORDER_X, sensitivity_first_order_1_x)
                node.SetValue(StructuralMechanicsApplication.NL_SENSITIVITY_SECOND_ORDER_X, sensitivity_second_order_x)
            else:
                node.SetValue(StructuralMechanicsApplication.NL_SENSITIVITY_FIRST_ORDER_X, 0.0)
                node.SetValue(StructuralMechanicsApplication.NL_SENSITIVITY_SECOND_ORDER_X, 0.0)

            if abs(disp_1_y) > 1e-8:
                sensitivity_first_order_1_y = disp_2_y / ( disp_1_y * f_1 )
                #sensitivity_first_order_2_y = disp_3_y / ( disp_1_y * f_2 )
                #sensitivity_second_order_y = sensitivity_first_order_2_y / sensitivity_first_order_1_y
                slope_10_y = (disp_2_y-disp_1_y)/delta_10
                slope_20_y = (disp_3_y-disp_1_y)/delta_20
                sensitivity_second_order_y = slope_20_y / slope_10_y
                node.SetValue(StructuralMechanicsApplication.NL_SENSITIVITY_FIRST_ORDER_Y, sensitivity_first_order_1_y)
                node.SetValue(StructuralMechanicsApplication.NL_SENSITIVITY_SECOND_ORDER_Y, sensitivity_second_order_y)
            else:
                node.SetValue(StructuralMechanicsApplication.NL_SENSITIVITY_FIRST_ORDER_Y, 0.0)
                node.SetValue(StructuralMechanicsApplication.NL_SENSITIVITY_SECOND_ORDER_Y, 0.0)

            if abs(disp_1_z) > 1e-8:
                sensitivity_first_order_1_z = disp_2_z / ( disp_1_z * f_1 )
                #sensitivity_first_order_2_z = disp_3_z / ( disp_1_z * f_2 )
                #sensitivity_second_order_z = sensitivity_first_order_2_z / sensitivity_first_order_1_z
                slope_10_z = (disp_2_z-disp_1_z)/delta_10
                slope_20_z = (disp_3_z-disp_1_z)/delta_20
                sensitivity_second_order_z = slope_20_z / slope_10_z
                node.SetValue(StructuralMechanicsApplication.NL_SENSITIVITY_FIRST_ORDER_Z, sensitivity_first_order_1_z)
                node.SetValue(StructuralMechanicsApplication.NL_SENSITIVITY_SECOND_ORDER_Z, sensitivity_second_order_z)
            else:
                node.SetValue(StructuralMechanicsApplication.NL_SENSITIVITY_FIRST_ORDER_Z, 0.0)
                node.SetValue(StructuralMechanicsApplication.NL_SENSITIVITY_SECOND_ORDER_Z, 0.0)




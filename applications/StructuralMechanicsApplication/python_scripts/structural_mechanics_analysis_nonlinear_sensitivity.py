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

    #def FinalizeSolutionStep(self):
    #    # Check if we are passing already in the last pseudo time step
    #    if (self.time < self.end_time) == False:
    #        # Compute displacement E(F)-sensitivity
    #        model_part = self._GetSolver().GetComputingModelPart()
    #        time_step = self._GetSolver().ComputeDeltaTime()
    #        load_factors = np.array([self.start_time, self.start_time + time_step, self.start_time + time_step*2])
    #        self._ComputeEFDisplacementCurvature(model_part, load_factors)
    #        self._ComputeFirstAndSecondOrderNLDisplacementSensitivityFactors(model_part, load_factors)
#
    #    super(StructuralMechanicsAnalysisNLSensitivity, self).FinalizeSolutionStep()


    #def _ComputeEFDisplacementCurvature(self, model_part, load_factor_array):
    #    for node in model_part.Nodes:
    #        disp_1_x = abs( node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X, 2) )
    #        disp_2_x = abs( node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X, 1) )
    #        disp_3_x = abs( node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X, 0) )
    #        disp_1_y = abs( node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y, 2) )
    #        disp_2_y = abs( node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y, 1) )
    #        disp_3_y = abs( node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y, 0) )
    #        disp_1_z = abs( node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Z, 2) )
    #        disp_2_z = abs( node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Z, 1) )
    #        disp_3_z = abs( node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Z, 0) )
#
    #        # X
    #        disp_x = np.array([disp_1_x, disp_2_x, disp_3_x])
    #        curvature_x = self._ComputeEFCurvature(disp_x, load_factor_array)
    #        node.SetValue(StructuralMechanicsApplication.DISPLACEMENT_NL_SENSITIVITY_X, curvature_x)
    #        # Y
    #        disp_y = np.array([disp_1_y, disp_2_y, disp_3_y])
    #        curvature_y = self._ComputeEFCurvature(disp_y, load_factor_array)
    #        node.SetValue(StructuralMechanicsApplication.DISPLACEMENT_NL_SENSITIVITY_Y, curvature_y)
    #        # Z
    #        disp_z = np.array([disp_1_z, disp_2_z, disp_3_z])
    #        curvature_z = self._ComputeEFCurvature(disp_z, load_factor_array)
    #        node.SetValue(StructuralMechanicsApplication.DISPLACEMENT_NL_SENSITIVITY_Z, curvature_z)
#
    #def _ComputeEFCurvature(self, response_value_array, load_factor_array):
    #    polynom = np.polyfit(load_factor_array, response_value_array, 2)
    #    curvature = 2 * polynom[0]
    #    return curvature
#
#
    #def _ComputeFirstAndSecondOrderNLDisplacementSensitivityFactors(self, model_part, load_factor_array):
    #    for node in model_part.Nodes:
    #        disp_1_x = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X, 2)
    #        disp_2_x = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X, 1)
    #        disp_3_x = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X, 0)
    #        disp_1_y = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y, 2)
    #        disp_2_y = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y, 1)
    #        disp_3_y = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y, 0)
    #        disp_1_z = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Z, 2)
    #        disp_2_z = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Z, 1)
    #        disp_3_z = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Z, 0)
#
    #        # X
    #        disp_x = np.array([disp_1_x, disp_2_x, disp_3_x])
    #        sen_first_x, sen_second_x = self._ComputeFirstAndSecondOrderNLSensitivityFactors(disp_x, load_factor_array)
    #        node.SetValue(StructuralMechanicsApplication.DISPLACEMENT_NL_SENSITIVITY_FIRST_ORDER_X, sen_first_x)
    #        node.SetValue(StructuralMechanicsApplication.DISPLACEMENT_NL_SENSITIVITY_SECOND_ORDER_X, sen_second_x)
    #        # Y
    #        disp_y = np.array([disp_1_y, disp_2_y, disp_3_y])
    #        sen_first_y, sen_second_y = self._ComputeFirstAndSecondOrderNLSensitivityFactors(disp_y, load_factor_array)
    #        node.SetValue(StructuralMechanicsApplication.DISPLACEMENT_NL_SENSITIVITY_FIRST_ORDER_Y, sen_first_y)
    #        node.SetValue(StructuralMechanicsApplication.DISPLACEMENT_NL_SENSITIVITY_SECOND_ORDER_Y, sen_second_y)
    #        # Z
    #        disp_z = np.array([disp_1_z, disp_2_z, disp_3_z])
    #        sen_first_z, sen_second_z = self._ComputeFirstAndSecondOrderNLSensitivityFactors(disp_z, load_factor_array)
    #        node.SetValue(StructuralMechanicsApplication.DISPLACEMENT_NL_SENSITIVITY_FIRST_ORDER_Z, sen_first_z)
    #        node.SetValue(StructuralMechanicsApplication.DISPLACEMENT_NL_SENSITIVITY_SECOND_ORDER_Z, sen_second_z)
#
#
    #def _ComputeFirstAndSecondOrderNLSensitivityFactors(self, response_value_array, load_factor_array):
    #    lambda_0 = load_factor_array[0]
    #    lambda_1 = load_factor_array[1]
    #    lambda_2 = load_factor_array[2]
    #    f_1 = lambda_1 / lambda_0
    #    delta_10 = lambda_1 - lambda_0
    #    delta_20 = lambda_2 - lambda_0
#
    #    response_0 = response_value_array[0]
    #    response_1 = response_value_array[1]
    #    response_2 = response_value_array[2]
#
    #    sensitivity_first_order = 0.0
    #    sensitivity_second_order = 0.0
#
    #    if abs(response_0) > 1e-8:
    #        sensitivity_first_order = response_1 / ( response_0 * f_1 )
    #        slope_10 = (response_1 - response_0) / delta_10
    #        slope_20 = (response_2 - response_1) / delta_20
    #        sensitivity_second_order = slope_20 / slope_10
#
    #    return sensitivity_first_order, sensitivity_second_order
#





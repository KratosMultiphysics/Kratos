from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
import KratosMultiphysics
from KratosMultiphysics import Vector
import KratosMultiphysics.SwimmingDEMApplication as SDEM
#from . import recoverer
import parameters_tools as PT
import numpy as np
import math
import csv
import sys

class L2ErrorProjectionUtility:
    def __init__(self, model):
        self.model_part = model
        self.model = KratosMultiphysics.Model()
        self.element_name = "CalculateErrorL2Projection3D"
        self.error_model_part = self.model.CreateModelPart("ErrorModelPart")
        self.error_model_part.AddNodalSolutionStepVariable(SDEM.VECTORIAL_ERROR)
        self.error_model_part.AddNodalSolutionStepVariable(SDEM.ERROR_X)
        self.error_model_part.AddNodalSolutionStepVariable(SDEM.ERROR_Y)
        self.error_model_part.AddNodalSolutionStepVariable(SDEM.ERROR_Z)
        self.error_model_part.AddNodalSolutionStepVariable(SDEM.SCALAR_ERROR)
        model_part_cloner = KratosMultiphysics.ConnectivityPreserveModeler()
        model_part_cloner.GenerateModelPart(self.model_part, self.error_model_part, self.element_name)
        self.error_model_part.ProcessInfo = self.model_part.ProcessInfo
        self.DOFs = (SDEM.ERROR_X, SDEM.ERROR_Y, SDEM.ERROR_Z, SDEM.SCALAR_ERROR)
        self.AddDofs(self.DOFs)
        self.SetStrategy()

    def ProjectL2(self):
        self.ComputeVelocityError()
        self.ComputePressureError()
        self.Solve()

        self.velocity_error_projected = SDEM.L2ErrorProjection().GetL2VectorProjection(self.error_model_part)
        self.pressure_error_projected = SDEM.L2ErrorProjection().GetL2ScalarProjection(self.error_model_part)
        return self.velocity_error_projected, self.pressure_error_projected, self.error_model_part

    def ComputeVelocityError(self):
        mod_error = 0
        for node in self.error_model_part.Nodes:
            vectorial_error = Vector(node.GetSolutionStepValue(KratosMultiphysics.VELOCITY) - node.GetSolutionStepValue(SDEM.EXACT_VELOCITY))
            mod_error += np.sqrt(vectorial_error[0]**2 + vectorial_error[1]**2 + vectorial_error[2]**2)
            node.SetSolutionStepValue(SDEM.VECTORIAL_ERROR, vectorial_error)

        print("av_mod_error :", mod_error / len(self.error_model_part.Elements))

    def ComputePressureError(self):
        mod_error = 0
        for node in self.error_model_part.Nodes:
            scalar_error = node.GetSolutionStepValue(KratosMultiphysics.PRESSURE) - node.GetSolutionStepValue(SDEM.EXACT_PRESSURE)
            node.SetSolutionStepValue(SDEM.SCALAR_ERROR, scalar_error)

    def SetStrategy(self):
        scheme = KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme()

        linear_solver = KratosMultiphysics.AMGCLSolver()
        self.l2_projector_strategy = KratosMultiphysics.ResidualBasedLinearStrategy(self.error_model_part, scheme, linear_solver, False, False, False, False)

    def AddDofs(self, DOF_variables):
        for node in self.error_model_part.Nodes:
            for var in DOF_variables:
                node.AddDof(var)

    def Solve(self):
        print("\nSolving for the fluid acceleration...")
        sys.stdout.flush()
        self.l2_projector_strategy.Solve()

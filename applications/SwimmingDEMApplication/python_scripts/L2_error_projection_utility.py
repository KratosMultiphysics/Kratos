from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
import KratosMultiphysics
from KratosMultiphysics import Vector
import KratosMultiphysics.SwimmingDEMApplication as SDEM
import sys

class L2ErrorProjectionUtility:
    def __init__(self, model, parameters):
        """The default constructor of the class.

        Keyword arguments:
        self -- It signifies an instance of a class.
        model -- the container of the different model parts.
        """
        self.model_part = model

        self.u_characteristic = parameters["error_projection_parameters"]["u_characteristic"].GetDouble()
        for element in self.model_part.Elements:
            rho = element.Properties.GetValue(KratosMultiphysics.DENSITY)
            break
        self.p_characteristic = (1/2)*rho*self.u_characteristic**2

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
        self.ComputeDofsErrors(self.error_model_part)
        self.Solve()

        self.velocity_error_projected = self.VectorL2Projection(self.error_model_part)
        self.pressure_error_projected = self.ScalarL2Projection(self.error_model_part)

        return self.velocity_error_projected/self.u_characteristic, self.pressure_error_projected/self.p_characteristic, self.error_model_part

    def ComputeDofsErrors(self, error_model_part):
        SDEM.L2ErrorProjection().ComputeDofsErrors(self.error_model_part)

    def VectorL2Projection(self, error_model_part):
        return SDEM.L2ErrorProjection().GetL2VectorProjection(self.error_model_part)

    def ScalarL2Projection(self, error_model_part):
        return SDEM.L2ErrorProjection().GetL2ScalarProjection(self.error_model_part)

    def SetStrategy(self):
        scheme = KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme()
        linear_solver = KratosMultiphysics.AMGCLSolver()
        self.l2_projector_strategy = KratosMultiphysics.ResidualBasedLinearStrategy(self.error_model_part, scheme, linear_solver, False, False, False, False)

    def AddDofs(self, DOF_variables):
        for node in self.error_model_part.Nodes:
            for var in DOF_variables:
                node.AddDof(var)

    def Solve(self):
        print("\nSolving for the L2 error projection...")
        sys.stdout.flush()
        self.l2_projector_strategy.Solve()

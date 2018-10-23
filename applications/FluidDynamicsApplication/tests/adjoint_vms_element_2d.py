from KratosMultiphysics import *
from KratosMultiphysics.FluidDynamicsApplication import *

import KratosMultiphysics.KratosUnittest as KratosUnittest
import random

class AdjointVMSElement2D(KratosUnittest.TestCase):

    def setUp(self):
        self.delta_time = 1.0
        # create test model part
        self.model_part = ModelPart("test")
        self.model_part.AddNodalSolutionStepVariable(VELOCITY)
        self.model_part.AddNodalSolutionStepVariable(ACCELERATION)
        self.model_part.AddNodalSolutionStepVariable(MESH_VELOCITY)
        self.model_part.AddNodalSolutionStepVariable(PRESSURE)
        self.model_part.AddNodalSolutionStepVariable(VISCOSITY)
        self.model_part.AddNodalSolutionStepVariable(DENSITY)
        self.model_part.AddNodalSolutionStepVariable(BODY_FORCE)
        self.model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
        self.model_part.CreateNewNode(2, 1.0, 0.0, 0.0)
        self.model_part.CreateNewNode(3, 1.0, 1.0, 0.0)
        prop = self.model_part.GetProperties()[0]
        self.model_part.CreateNewElement("VMS2D3N", 1, [1, 2, 3], prop)
        self.model_part.CreateNewElement("VMSAdjointElement2D", 2, [1, 2, 3], prop)
        self.model_part.SetBufferSize(2)
        self.model_part.ProcessInfo[OSS_SWITCH] = 0
        self.model_part.ProcessInfo[DELTA_TIME] = self.delta_time
        self.model_part.ProcessInfo[DYNAMIC_TAU] = 1.0

        self.vms_element = self.model_part.GetElement(1)
        self.adjoint_element = self.model_part.GetElement(2)

        self._AssignSolutionStepData1(0)
        self._AssignSolutionStepData2(1)

    def _AssignSolutionStepData1(self, step=0):
        # generate nodal solution step test data
        random.seed(1.0)
        for node in self.model_part.Nodes:
            node.SetSolutionStepValue(DENSITY,step,1.0)
            node.SetSolutionStepValue(VISCOSITY,step,1.0e-5)
            node.SetSolutionStepValue(VELOCITY_X,step,random.random())
            node.SetSolutionStepValue(VELOCITY_Y,step,random.random())
            node.SetSolutionStepValue(ACCELERATION_X,step,random.random())
            node.SetSolutionStepValue(ACCELERATION_Y,step,random.random())
            node.SetSolutionStepValue(PRESSURE,step,random.random())

    def _AssignSolutionStepData2(self, step=0):
        # generate nodal solution step test data
        random.seed(2.0)
        for node in self.model_part.Nodes:
            node.SetSolutionStepValue(DENSITY,step,1.0)
            node.SetSolutionStepValue(VISCOSITY,step,1.0e-5)
            node.SetSolutionStepValue(VELOCITY_X,step,random.random())
            node.SetSolutionStepValue(VELOCITY_Y,step,random.random())
            node.SetSolutionStepValue(ACCELERATION_X,step,random.random())
            node.SetSolutionStepValue(ACCELERATION_Y,step,random.random())
            node.SetSolutionStepValue(PRESSURE,step,random.random())

    def _zeroVector(self,size):
        v = Vector(size)
        for i in range(size):
            v[i] = 0.0
        return v

    def _transpose(self, m):
        tmp = Matrix(m.Size1(), m.Size2())
        for i in range(m.Size1()):
            for j in range(m.Size2()):
                tmp[i,j] = m[j,i]
        return tmp

    def _assertMatrixAlmostEqual(self, matrix1, matrix2, prec=7):
        self.assertEqual(matrix1.Size1(), matrix2.Size1())
        self.assertEqual(matrix1.Size2(), matrix2.Size2())
        for i in range(matrix1.Size1()):
            for j in range(matrix1.Size2()):
                self.assertAlmostEqual(matrix1[i,j], matrix2[i,j], prec)

    def testCalculateSecondDerivativesLHS(self):
        Mass1 = Matrix(9,9)
        self.model_part.ProcessInfo[DELTA_TIME] = self.delta_time
        self.vms_element.CalculateMassMatrix(Mass1,self.model_part.ProcessInfo)
        self.model_part.ProcessInfo[DELTA_TIME] =-self.delta_time
        mass2_trans = Matrix(9,9)
        self.adjoint_element.CalculateSecondDerivativesLHS(mass2_trans,self.model_part.ProcessInfo)
        self._assertMatrixAlmostEqual(Mass1, self._transpose(mass2_trans)*(-1.0))

    def testCalculateFirstDerivativesLHS1(self):
        # test for steady state.
        for node in self.model_part.Nodes:
            for step in range(2):
                node.SetSolutionStepValue(ACCELERATION_X, step, 0.0)
                node.SetSolutionStepValue(ACCELERATION_Y, step, 0.0)
        # unperturbed residual
        LHS = Matrix(9,9)
        RHS = self._zeroVector(9)
        FirstDerivatives = Vector(9)
        self.model_part.ProcessInfo[DELTA_TIME] = self.delta_time
        self.vms_element.CalculateLocalVelocityContribution(LHS,RHS,self.model_part.ProcessInfo)
        self.vms_element.GetFirstDerivativesVector(FirstDerivatives,0)
        res0 = LHS * FirstDerivatives
        # finite difference approximation
        h = 0.0000001
        FDAdjointMatrix = Matrix(9,9)
        row_index = 0
        for node in self.model_part.Nodes:
            # VELOCITY_X
            vx = node.GetSolutionStepValue(VELOCITY_X,0)
            node.SetSolutionStepValue(VELOCITY_X,0,vx+h)
            self.vms_element.CalculateLocalVelocityContribution(LHS,RHS,self.model_part.ProcessInfo)
            self.vms_element.GetFirstDerivativesVector(FirstDerivatives,0)
            node.SetSolutionStepValue(VELOCITY_X,0,vx)
            res = LHS * FirstDerivatives
            for j in range(9):
                FDAdjointMatrix[row_index,j] = -(res[j] - res0[j]) / h
            row_index = row_index + 1
            # VELOCITY_Y
            vy = node.GetSolutionStepValue(VELOCITY_Y,0)
            node.SetSolutionStepValue(VELOCITY_Y,0,vy+h)
            self.vms_element.CalculateLocalVelocityContribution(LHS,RHS,self.model_part.ProcessInfo)
            self.vms_element.GetFirstDerivativesVector(FirstDerivatives,0)
            node.SetSolutionStepValue(VELOCITY_Y,0,vy)
            res = LHS * FirstDerivatives
            for j in range(9):
                FDAdjointMatrix[row_index,j] = -(res[j] - res0[j]) / h
            row_index = row_index + 1
            # PRESSURE
            p = node.GetSolutionStepValue(PRESSURE,0)
            node.SetSolutionStepValue(PRESSURE,0,p+h)
            self.vms_element.CalculateLocalVelocityContribution(LHS,RHS,self.model_part.ProcessInfo)
            self.vms_element.GetFirstDerivativesVector(FirstDerivatives,0)
            node.SetSolutionStepValue(PRESSURE,0,p)
            res = LHS * FirstDerivatives
            for j in range(9):
                FDAdjointMatrix[row_index,j] = -(res[j] - res0[j]) / h
            row_index = row_index + 1

        # analytical implementation
        self.model_part.ProcessInfo[DELTA_TIME] =-self.delta_time
        AdjointMatrix = Matrix(9,9)
        self.adjoint_element.CalculateFirstDerivativesLHS(AdjointMatrix,self.model_part.ProcessInfo)
        self._assertMatrixAlmostEqual(FDAdjointMatrix, AdjointMatrix)
        # reset test data
        self._AssignSolutionStepData1(0)
        self._AssignSolutionStepData2(1)

    def testCalculateFirstDerivativesLHS2(self):
        # unperturbed residual
        Mass = Matrix(9,9)
        LHS = Matrix(9,9)
        RHS = self._zeroVector(9)
        FirstDerivatives = Vector(9)
        SecondDerivatives = Vector(9)
        self.model_part.ProcessInfo[DELTA_TIME] = self.delta_time
        self.vms_element.CalculateMassMatrix(Mass,self.model_part.ProcessInfo)
        self.vms_element.CalculateLocalVelocityContribution(LHS,RHS,self.model_part.ProcessInfo)
        self.vms_element.GetFirstDerivativesVector(FirstDerivatives,0)
        self.vms_element.GetSecondDerivativesVector(SecondDerivatives,0)
        res0 = LHS * FirstDerivatives + Mass * SecondDerivatives
        # finite difference approximation
        h = 0.0000001
        FDAdjointMatrix = Matrix(9,9)
        row_index = 0
        for node in self.model_part.Nodes:
            # VELOCITY_X
            vx = node.GetSolutionStepValue(VELOCITY_X,0)
            node.SetSolutionStepValue(VELOCITY_X,0,vx+h)
            self.vms_element.CalculateMassMatrix(Mass,self.model_part.ProcessInfo)
            self.vms_element.CalculateLocalVelocityContribution(LHS,RHS,self.model_part.ProcessInfo)
            self.vms_element.GetFirstDerivativesVector(FirstDerivatives,0)
            self.vms_element.GetSecondDerivativesVector(SecondDerivatives,0)
            node.SetSolutionStepValue(VELOCITY_X,0,vx)
            res = LHS * FirstDerivatives + Mass * SecondDerivatives
            for j in range(9):
                FDAdjointMatrix[row_index,j] = -(res[j] - res0[j]) / h
            row_index = row_index + 1
            # VELOCITY_Y
            vy = node.GetSolutionStepValue(VELOCITY_Y,0)
            node.SetSolutionStepValue(VELOCITY_Y,0,vy+h)
            self.vms_element.CalculateMassMatrix(Mass,self.model_part.ProcessInfo)
            self.vms_element.CalculateLocalVelocityContribution(LHS,RHS,self.model_part.ProcessInfo)
            self.vms_element.GetFirstDerivativesVector(FirstDerivatives,0)
            self.vms_element.GetSecondDerivativesVector(SecondDerivatives,0)
            node.SetSolutionStepValue(VELOCITY_Y,0,vy)
            res = LHS * FirstDerivatives + Mass * SecondDerivatives
            for j in range(9):
                FDAdjointMatrix[row_index,j] = -(res[j] - res0[j]) / h
            row_index = row_index + 1
            # PRESSURE
            p = node.GetSolutionStepValue(PRESSURE,0)
            node.SetSolutionStepValue(PRESSURE,0,p+h)
            self.vms_element.CalculateMassMatrix(Mass,self.model_part.ProcessInfo)
            self.vms_element.CalculateLocalVelocityContribution(LHS,RHS,self.model_part.ProcessInfo)
            self.vms_element.GetFirstDerivativesVector(FirstDerivatives,0)
            self.vms_element.GetSecondDerivativesVector(SecondDerivatives,0)
            node.SetSolutionStepValue(PRESSURE,0,p)
            res = LHS * FirstDerivatives + Mass * SecondDerivatives
            for j in range(9):
                FDAdjointMatrix[row_index,j] = -(res[j] - res0[j]) / h
            row_index = row_index + 1

        # analytical implementation
        self.model_part.ProcessInfo[DELTA_TIME] =-self.delta_time
        AdjointMatrix = Matrix(9,9)
        self.adjoint_element.CalculateFirstDerivativesLHS(AdjointMatrix,self.model_part.ProcessInfo)
        self._assertMatrixAlmostEqual(FDAdjointMatrix, AdjointMatrix)

    def testCalculateSensitivityMatrix(self):
        # unperturbed residual
        Mass = Matrix(9,9)
        LHS = Matrix(9,9)
        RHS = self._zeroVector(9)
        FirstDerivatives = Vector(9)
        SecondDerivatives = Vector(9)
        self.model_part.ProcessInfo[DELTA_TIME] = self.delta_time
        self.vms_element.CalculateMassMatrix(Mass,self.model_part.ProcessInfo)
        self.vms_element.CalculateLocalVelocityContribution(LHS,RHS,self.model_part.ProcessInfo)
        self.vms_element.GetFirstDerivativesVector(FirstDerivatives,0)
        self.vms_element.GetSecondDerivativesVector(SecondDerivatives,0)
        res0 = LHS * FirstDerivatives + Mass * SecondDerivatives
        # finite difference approximation
        h = 0.00000001
        FDShapeDerivativeMatrix = Matrix(6,9)
        row_index = 0
        for node in self.model_part.Nodes:
            # X
            x = node.X
            node.X = x+h
            self.vms_element.CalculateMassMatrix(Mass,self.model_part.ProcessInfo)
            self.vms_element.CalculateLocalVelocityContribution(LHS,RHS,self.model_part.ProcessInfo)
            node.X = x
            res = LHS * FirstDerivatives + Mass * SecondDerivatives
            for j in range(9):
                FDShapeDerivativeMatrix[row_index,j] = -(res[j] - res0[j]) / h
            row_index = row_index + 1
            # Y
            y = node.Y
            node.Y = y+h
            self.vms_element.CalculateMassMatrix(Mass,self.model_part.ProcessInfo)
            self.vms_element.CalculateLocalVelocityContribution(LHS,RHS,self.model_part.ProcessInfo)
            node.Y = y
            res = LHS * FirstDerivatives + Mass * SecondDerivatives
            for j in range(9):
                FDShapeDerivativeMatrix[row_index,j] = -(res[j] - res0[j]) / h
            row_index = row_index + 1

        # analytical implementation
        self.model_part.ProcessInfo[DELTA_TIME] =-self.delta_time
        ShapeDerivativeMatrix = Matrix(6,9)
        self.adjoint_element.CalculateSensitivityMatrix(SHAPE_SENSITIVITY,ShapeDerivativeMatrix,self.model_part.ProcessInfo)
        self._assertMatrixAlmostEqual(FDShapeDerivativeMatrix, ShapeDerivativeMatrix)

if __name__ == '__main__':
    KratosUnittest.main()

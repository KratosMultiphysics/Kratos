from KratosMultiphysics import *
from KratosMultiphysics.FluidDynamicsApplication import *
from KratosMultiphysics.AdjointFluidApplication import *
import KratosMultiphysics.KratosUnittest as KratosUnittest
import random

class TestCase(KratosUnittest.TestCase):

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

        self.AssignSolutionStepDataSet1(0)
        self.AssignSolutionStepDataSet2(1)

    def AssignSolutionStepDataSet1(self, step=0):
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

    def AssignSolutionStepDataSet2(self, step=0):
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

    def zeroVector(self,size):
        v = Vector(size)
        for i in range(size):
            v[i] = 0.0
        return v

    def assertMatrixAlmostEqual(self, matrix1, matrix2, prec=7):
        self.assertEqual(matrix1.Size1(), matrix2.Size1())
        self.assertEqual(matrix1.Size2(), matrix2.Size2())
        for i in range(matrix1.Size1()):
            for j in range(matrix1.Size2()):
                self.assertAlmostEqual(matrix1[i,j], matrix2[i,j], prec)

    def test_MASS_MATRIX_0(self):
        Mass1 = Matrix(9,9)
        self.model_part.ProcessInfo[DELTA_TIME] = self.delta_time
        self.vms_element.CalculateMassMatrix(Mass1,self.model_part.ProcessInfo)
        self.model_part.ProcessInfo[DELTA_TIME] =-self.delta_time
        Mass2 = self.adjoint_element.Calculate(MASS_MATRIX_0,self.model_part.ProcessInfo)
        self.assertMatrixAlmostEqual(Mass1, Mass2)

    def test_MASS_MATRIX_1(self):
        Mass1 = Matrix(9,9)
        self.AssignSolutionStepDataSet2(0)
        self.model_part.ProcessInfo[DELTA_TIME] = self.delta_time
        self.vms_element.CalculateMassMatrix(Mass1,self.model_part.ProcessInfo)
        self.AssignSolutionStepDataSet1(0)
        self.model_part.ProcessInfo[DELTA_TIME] =-self.delta_time
        Mass2 = self.adjoint_element.Calculate(MASS_MATRIX_1,self.model_part.ProcessInfo)
        self.assertMatrixAlmostEqual(Mass1, Mass2)

    def test_ADJOINT_MATRIX_1(self):
        # unperturbed residual
        LHS = Matrix(9,9)
        RHS = self.zeroVector(9)
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
        AdjointMatrix = self.adjoint_element.Calculate(ADJOINT_MATRIX_1,self.model_part.ProcessInfo)
        self.assertMatrixAlmostEqual(FDAdjointMatrix, AdjointMatrix)

    def test_ADJOINT_MATRIX_2(self):
        # unperturbed residual
        Mass = Matrix(9,9)
        LHS = Matrix(9,9)
        RHS = self.zeroVector(9)
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
        AdjointMatrix = self.adjoint_element.Calculate(ADJOINT_MATRIX_2,self.model_part.ProcessInfo)
        self.assertMatrixAlmostEqual(FDAdjointMatrix, AdjointMatrix)

    def test_SHAPE_DERIVATIVE_MATRIX_1(self):
        # unperturbed residual
        LHS = Matrix(9,9)
        RHS = self.zeroVector(9)
        FirstDerivatives = Vector(9)
        self.model_part.ProcessInfo[DELTA_TIME] = self.delta_time
        self.vms_element.CalculateLocalVelocityContribution(LHS,RHS,self.model_part.ProcessInfo)
        self.vms_element.GetFirstDerivativesVector(FirstDerivatives,0)
        res0 = LHS * FirstDerivatives
        # finite difference approximation
        h = 0.00000001
        FDShapeDerivativeMatrix = Matrix(6,9)
        row_index = 0
        for node in self.model_part.Nodes:
            # X
            x = node.X
            node.X = x+h
            self.vms_element.CalculateLocalVelocityContribution(LHS,RHS,self.model_part.ProcessInfo)
            node.X = x
            res = LHS * FirstDerivatives
            for j in range(9):
                FDShapeDerivativeMatrix[row_index,j] = -(res[j] - res0[j]) / h
            row_index = row_index + 1
            # Y
            y = node.Y
            node.Y = y+h
            self.vms_element.CalculateLocalVelocityContribution(LHS,RHS,self.model_part.ProcessInfo)
            node.Y = y
            res = LHS * FirstDerivatives
            for j in range(9):
                FDShapeDerivativeMatrix[row_index,j] = -(res[j] - res0[j]) / h
            row_index = row_index + 1
            
        # analytical implementation
        self.model_part.ProcessInfo[DELTA_TIME] =-self.delta_time
        ShapeDerivativeMatrix = self.adjoint_element.Calculate(SHAPE_DERIVATIVE_MATRIX_1,self.model_part.ProcessInfo)
        self.assertMatrixAlmostEqual(FDShapeDerivativeMatrix, ShapeDerivativeMatrix)

    def test_SHAPE_DERIVATIVE_MATRIX_2(self):
        # unperturbed residual
        Mass = Matrix(9,9)
        LHS = Matrix(9,9)
        RHS = self.zeroVector(9)
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
        ShapeDerivativeMatrix = self.adjoint_element.Calculate(SHAPE_DERIVATIVE_MATRIX_2,self.model_part.ProcessInfo)
        self.assertMatrixAlmostEqual(FDShapeDerivativeMatrix, ShapeDerivativeMatrix)

if __name__ == '__main__':
    KratosUnittest.main()

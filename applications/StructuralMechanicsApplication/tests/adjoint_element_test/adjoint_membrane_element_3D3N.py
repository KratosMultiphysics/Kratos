from KratosMultiphysics import *
from KratosMultiphysics.StructuralMechanicsApplication import *
import KratosMultiphysics.StructuralMechanicsApplication as KratosSM

import KratosMultiphysics.KratosUnittest as KratosUnittest
import random

class AdjointMembraneElement3D3N(KratosUnittest.TestCase):

    def setUp(self):
        self.delta_time = 1.0
        # create test model part
        self.model = Model()
        self.model_part = self.model.CreateModelPart("test")
        self.model_part.ProcessInfo.SetValue(DOMAIN_SIZE, 3)
        self.model_part.AddNodalSolutionStepVariable(DISPLACEMENT)
        self.model_part.AddNodalSolutionStepVariable(VELOCITY)
        self.model_part.AddNodalSolutionStepVariable(ACCELERATION)
        self.model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
        self.model_part.CreateNewNode(2, 1.0, 0.0, 0.0)
        self.model_part.CreateNewNode(3, 0.0, 1.0, 0.0)

        prop_parameters = Parameters("""{
            "properties": [{
                "model_part_name": "test",
                "properties_id": 1,
                "Material": {
                    "constitutive_law": {
                        "name": "LinearElasticPlaneStress2DLaw"
                    },
                    "Variables": {
                        "THICKNESS"               : 0.001,
                        "PRESTRESS_VECTOR"        : [4.0e3,4.0e3,0],
                        "DENSITY": 7850.0,
                        "YOUNG_MODULUS": 2e4,
                        "POISSON_RATIO": 0.3
                    }
                }
            }]
        }""")
        ReadMaterialsUtility(self.model).ReadMaterials(prop_parameters)
        prop = self.model_part.GetProperties(1)

        self.model_part.CreateNewElement("MembraneElement3D3N", 1, [1, 2, 3], prop)
        self.model_part.CreateNewElement("AdjointFiniteDifferenceMembraneElement3D3N", 2, [1, 2, 3], prop)
        self.model_part.SetBufferSize(2)
        self.model_part.ProcessInfo[OSS_SWITCH] = 0
        self.model_part.ProcessInfo[DELTA_TIME] = self.delta_time
   
        self.primal_element = self.model_part.GetElement(1)
        self.adjoint_element = self.model_part.GetElement(2)

        self._AssignSolutionStepData1(0)
        self._AssignSolutionStepData2(1)

        force = [100, 0, 0]
        self.model_part.GetNode(1).SetValue(FORCE, force)

        force2 = [0, 100, 0]
        self.model_part.GetNode(2).SetValue(FORCE, force2)

    def _AssignSolutionStepData1(self, step=0):
        # generate nodal solution step test data
        random.seed(1.0)
        for node in self.model_part.Nodes:
            node.SetSolutionStepValue(DISPLACEMENT_X,step,random.random())
            node.SetSolutionStepValue(DISPLACEMENT_Y,step,random.random())
            node.SetSolutionStepValue(DISPLACEMENT_Z,step,random.random())
            node.SetSolutionStepValue(VELOCITY_X,step,random.random())
            node.SetSolutionStepValue(VELOCITY_Y,step,random.random())
            node.SetSolutionStepValue(VELOCITY_Z,step,random.random())
            node.SetSolutionStepValue(ACCELERATION_X,step,random.random())
            node.SetSolutionStepValue(ACCELERATION_Y,step,random.random())
            node.SetSolutionStepValue(ACCELERATION_Z,step,random.random())

    def _AssignSolutionStepData2(self, step=0):
        # generate nodal solution step test data
        random.seed(2.0)
        for node in self.model_part.Nodes:
            node.SetSolutionStepValue(DISPLACEMENT_X,step,random.random())
            node.SetSolutionStepValue(DISPLACEMENT_Y,step,random.random())
            node.SetSolutionStepValue(DISPLACEMENT_Z,step,random.random())
            node.SetSolutionStepValue(VELOCITY_X,step,random.random())
            node.SetSolutionStepValue(VELOCITY_Y,step,random.random())
            node.SetSolutionStepValue(VELOCITY_Z,step,random.random())
            node.SetSolutionStepValue(ACCELERATION_X,step,random.random())
            node.SetSolutionStepValue(ACCELERATION_Y,step,random.random())
            node.SetSolutionStepValue(ACCELERATION_Z,step,random.random())

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

    def _assertVectorAlmostEqual(self, vector1, vector2, prec=7):
        self.assertEqual(vector1.Size(), vector2.Size())
        for i in range(vector1.Size()):
            self.assertAlmostEqual(vector1[i], vector2[i], prec)

    def testCalculateSecondDerivativesLHS(self):
        Mass1 = Matrix(9, 9)
        self.model_part.ProcessInfo[DELTA_TIME] = self.delta_time
        self.primal_element.CalculateMassMatrix(Mass1,self.model_part.ProcessInfo)
        print("Mass from primal ", Mass1)
        self.model_part.ProcessInfo[DELTA_TIME] =-self.delta_time
        mass2_trans = Matrix(9, 9)
        self.adjoint_element.CalculateSecondDerivativesLHS(mass2_trans,self.model_part.ProcessInfo)
        print("CalculateSecondDerivativesLHS ", mass2_trans)
        self._assertMatrixAlmostEqual(Mass1, self._transpose(mass2_trans))

    def testCalculateFirstDerivativesLHS(self):
        Damp1 = Matrix(9, 9)
        self.model_part.ProcessInfo[DELTA_TIME] = self.delta_time
        self.primal_element.CalculateDampingMatrix(Damp1,self.model_part.ProcessInfo)
        print("Damping from primal ", Damp1)
        self.model_part.ProcessInfo[DELTA_TIME] =-self.delta_time
        damp2_trans = Matrix(9, 9)
        self.adjoint_element.CalculateFirstDerivativesLHS(damp2_trans,self.model_part.ProcessInfo)
        print("testCalculateFirstDerivativesLHS  ", damp2_trans)
        self._assertMatrixAlmostEqual(Damp1, self._transpose(damp2_trans))

    def testCalculateSecondDerivativesSensitivityRHS(self):
        # Res1 = Vector(9)
        epsilon = 1e-7
        self.model_part.ProcessInfo[DELTA_TIME] = self.delta_time
        self.primal_element.Initialize(self.model_part.ProcessInfo)
        # self.primal_element.CalculateRightHandSide(Res1,self.model_part.ProcessInfo)
        # print(Res1)

        #unperturbed residual
        Mass = Matrix(9, 9)
        LHS = Matrix(9, 9)
        RHS = self._zeroVector(9)
        Values = Vector(9)
        SecondDerivatives = Vector(9)
        #self.model_part.ProcessInfo[DELTA_TIME] = self.delta_time
        self.primal_element.Initialize(self.model_part.ProcessInfo)
        self.primal_element.CalculateMassMatrix(Mass,self.model_part.ProcessInfo)
        self.primal_element.CalculateLeftHandSide(LHS,self.model_part.ProcessInfo)
        self.primal_element.GetValuesVector(Values,0)
        self.primal_element.GetSecondDerivativesVector(SecondDerivatives,0)
        res0 = LHS * Values + Mass * SecondDerivatives * (1 - (-0.3))

        # self.model_part.ProcessInfo[DELTA_TIME] =self.delta_time
        # res1 = Vector(9)
        # self.adjoint_element.Initialize(self.model_part.ProcessInfo)
        # self.adjoint_element.CalculateRightHandSide(res1,self.model_part.ProcessInfo)
        print(res0)

        dof_list = [ACCELERATION_X, ACCELERATION_Y, ACCELERATION_Z]
        for index, node in enumerate(self.model_part.Nodes):
            node: Node
            for j, dof in enumerate(dof_list):

                nodal_val = self.model_part.GetNode(node.Id).GetSolutionStepValue(dof)
                self.model_part.GetNode(node.Id).SetSolutionStepValue(dof, nodal_val*(1+ epsilon))

                Mass2 = Matrix(9, 9)
                LHS2 = Matrix(9, 9)
                RHS = self._zeroVector(9)
                Values2 = Vector(9)
                SecondDerivatives2 = Vector(9)
                self.model_part.ProcessInfo[DELTA_TIME] = self.delta_time
                self.primal_element.CalculateMassMatrix(Mass2,self.model_part.ProcessInfo)
                self.primal_element.CalculateLeftHandSide(LHS2,self.model_part.ProcessInfo)
                self.primal_element.GetValuesVector(Values2,0)
                self.primal_element.GetSecondDerivativesVector(SecondDerivatives2,0)
                res2 = LHS2 * Values2 + Mass2 * SecondDerivatives2 * (1 - (-0.3) )

                # res2 = Vector(9)
                # self.adjoint_element.Initialize(self.model_part.ProcessInfo)
                # self.adjoint_element.CalculateRightHandSide(res2,self.model_part.ProcessInfo)
                #print(res2)
                FD_n1 = (res2 - res0) / (nodal_val * epsilon)
                self.model_part.GetNode(node.Id).SetSolutionStepValue(dof, nodal_val)

                col1 = Vector(9)
                for i in range(9):
                    col1[i] = Mass[i,index*3+j] * (1 - (-0.3))


                self._assertVectorAlmostEqual(FD_n1, col1*(1.0),6)

                print("---------------------------")
                print(f"For Node {node.Id} and dof {dof} :")
                print("From FD partial R/ partial u.. : ", FD_n1)
                # #print(mass2_trans)
                print("From CalculateMassMatrix       : ", col1)
                print("---------------------------")

    # def testCalculateSensitivityRHS(self):
    #     # Res1 = Vector(9)
    #     # self.model_part.ProcessInfo[DELTA_TIME] = self.delta_time
    #     self.primal_element.Initialize(self.model_part.ProcessInfo)
    #     # self.primal_element.CalculateRightHandSide(Res1,self.model_part.ProcessInfo)
    #     # print(Res1)

    #     #unperturbed residual
    #     Mass = Matrix(9, 9)
    #     LHS = Matrix(9, 9)
    #     RHS = self._zeroVector(9)
    #     Values = Vector(9)
    #     SecondDerivatives = Vector(9)
    #     #self.model_part.ProcessInfo[DELTA_TIME] = self.delta_time
    #     #self.primal_element.Initialize(self.model_part.ProcessInfo)
    #     self.primal_element.CalculateMassMatrix(Mass,self.model_part.ProcessInfo)
    #     self.primal_element.CalculateLeftHandSide(LHS,self.model_part.ProcessInfo)
    #     self.primal_element.GetValuesVector(Values,0)
    #     self.primal_element.GetSecondDerivativesVector(SecondDerivatives,0)
    #     res0 = LHS * Values + Mass * SecondDerivatives * (1 - (-0.3))

    #     # self.model_part.ProcessInfo[DELTA_TIME] =self.delta_time
    #     # stiff2_trans = Matrix(9, 9)
    #     # self.adjoint_element.Initialize(self.model_part.ProcessInfo)
    #     # self.adjoint_element.CalculateLeftHandSide(stiff2_trans,self.model_part.ProcessInfo)
    #     # col1 = Vector(9)
    #     # for i in range(9):
    #     #     col1[i] = stiff2_trans[i,0]

    #     # self.model_part.ProcessInfo[DELTA_TIME] =self.delta_time
    #     # res1 = Vector(9)
    #     # self.adjoint_element.Initialize(self.model_part.ProcessInfo)
    #     # self.adjoint_element.CalculateRightHandSide(res1,self.model_part.ProcessInfo)
    #     #print(res1)

    #     dof_list = [DISPLACEMENT_X, DISPLACEMENT_Y, DISPLACEMENT_Z]
    #     for index, node in enumerate(self.model_part.Nodes):
    #         node: Node
    #         for j, dof in enumerate(dof_list):

    #             nodal_val = self.model_part.GetNode(node.Id).GetSolutionStepValue(dof)
    #             self.model_part.GetNode(node.Id).SetSolutionStepValue(dof, nodal_val*(1+ 1e-5))

    #             Mass2 = Matrix(9, 9)
    #             LHS2 = Matrix(9, 9)
    #             RHS = self._zeroVector(9)
    #             Values2 = Vector(9)
    #             SecondDerivatives2 = Vector(9)
    #             self.model_part.ProcessInfo[DELTA_TIME] = self.delta_time
    #             self.primal_element.CalculateMassMatrix(Mass2,self.model_part.ProcessInfo)
    #             self.primal_element.CalculateLeftHandSide(LHS2,self.model_part.ProcessInfo)
    #             self.primal_element.GetValuesVector(Values2,0)
    #             self.primal_element.GetSecondDerivativesVector(SecondDerivatives2,0)
    #             res2 = LHS2 * Values2 + Mass2 * SecondDerivatives2 * (1 - (-0.3) )


    #             # res2 = Vector(9)
    #             # self.adjoint_element.Initialize(self.model_part.ProcessInfo)
    #             # self.adjoint_element.CalculateRightHandSide(res2,self.model_part.ProcessInfo)
    #             #print(res2)
    #             FD_n1 = (res2 - res0) / (nodal_val * 1e-5)
    #             self.model_part.GetNode(node.Id).SetSolutionStepValue(dof, nodal_val)

    #             col1 = Vector(9)
    #             for i in range(9):
    #                 col1[i] = LHS[i,index*3+j]


    #             self._assertVectorAlmostEqual(FD_n1, col1*(1.0),5)

    #             print("---------------------------")
    #             print(f"For Node {node.Id} and dof {dof} :")
    #             print("From FD partial R/ partial u : ", FD_n1)
    #             # #print(LHS)
    #             print("From CalculateLHS            : ", col1*(1.0))
    #             print("---------------------------")

    def testCalculateSensitivityRHS(self):
        # Res1 = Vector(9)
        # self.model_part.ProcessInfo[DELTA_TIME] = self.delta_time
        # self.primal_element.Initialize(self.model_part.ProcessInfo)
        # self.primal_element.CalculateRightHandSide(Res1,self.model_part.ProcessInfo)
        # print(Res1)
        epsilon = 1e-7
        self.model_part.ProcessInfo[DELTA_TIME] =self.delta_time
        stiff2_trans = Matrix(9, 9)
        self.adjoint_element.Initialize(self.model_part.ProcessInfo)
        self.adjoint_element.CalculateLeftHandSide(stiff2_trans,self.model_part.ProcessInfo)

        self.model_part.ProcessInfo[DELTA_TIME] =self.delta_time
        res1 = Vector(9)
        self.adjoint_element.Initialize(self.model_part.ProcessInfo)
        self.adjoint_element.CalculateRightHandSide(res1,self.model_part.ProcessInfo)
        #print(res1)

        dof_list = [DISPLACEMENT_X, DISPLACEMENT_Y, DISPLACEMENT_Z]
        for index, node in enumerate(self.model_part.Nodes):
            node: Node
            for j, dof in enumerate(dof_list):

                nodal_val = self.model_part.GetNode(node.Id).GetSolutionStepValue(dof)
                self.model_part.GetNode(node.Id).SetSolutionStepValue(dof, nodal_val*(1+ epsilon))
                res2 = Vector(9)
                self.adjoint_element.Initialize(self.model_part.ProcessInfo)
                self.adjoint_element.CalculateRightHandSide(res2,self.model_part.ProcessInfo)
                #print(res2)
                FD_n1 = (res2 - res1) / (nodal_val * epsilon)
                self.model_part.GetNode(node.Id).SetSolutionStepValue(dof, nodal_val)

                col1 = Vector(9)
                for i in range(9):
                    col1[i] = stiff2_trans[i,index*3+j]


                self._assertVectorAlmostEqual(FD_n1, col1*(-1.0),5)

                print("---------------------------")
                print(f"For Node {node.Id} and dof {dof} :")
                print("From FD partial R/ partial u : ", FD_n1)
                # #print(stiff2_trans)
                print("From CalculateLHS            : ", col1*(-1.0))
                print("---------------------------")


if __name__ == '__main__':
    KratosUnittest.main()

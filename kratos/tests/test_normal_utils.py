from __future__ import print_function, absolute_import, division

import os
import math

import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.gid_output_process import GiDOutputProcess
from KratosMultiphysics.testing.utilities import ReadModelPart

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

def PostProcess(model_part):
    gid_output = GiDOutputProcess(model_part,
                                "gid_output",
                                KratosMultiphysics.Parameters("""
                                    {
                                        "result_file_configuration" : {
                                            "gidpost_flags": {
                                                "GiDPostMode": "GiD_PostBinary",
                                                "WriteDeformedMeshFlag": "WriteUndeformed",
                                                "WriteConditionsFlag": "WriteConditions",
                                                "MultiFileFlag": "SingleFile"
                                            },
                                            "nodal_results" : ["NORMAL"]
                                        }
                                    }
                                    """)
                                )

    gid_output.ExecuteInitialize()
    gid_output.ExecuteBeforeSolutionLoop()
    gid_output.ExecuteInitializeSolutionStep()
    gid_output.PrintOutput()
    gid_output.ExecuteFinalizeSolutionStep()
    gid_output.ExecuteFinalize()

def CalculateAnalyticalNormal(node):
    normal = []
    norm = math.sqrt(node.X**2+node.Y**2+node.Z**2)
    normal.append(node.X/norm)
    normal.append(node.Y/norm)
    normal.append(node.Z/norm)
    return normal

def CalculateNormalResidual(analytical_normal, utilities_normal):
    return math.sqrt((utilities_normal[0]-analytical_normal[0])**2+(utilities_normal[1]-analytical_normal[1])**2+(utilities_normal[2]-analytical_normal[2])**2)


class TestNormalUtilsCoarseSphere(KratosUnittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.current_model = KratosMultiphysics.Model()
        cls.model_part = cls.current_model.CreateModelPart("Main")
        cls.model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] = 3
        cls.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL)
        model_part_io = KratosMultiphysics.ModelPartIO(GetFilePath("auxiliar_files_for_python_unittest/mdpa_files/coarse_sphere"))
        model_part_io.ReadModelPart(cls.model_part)

    def setUp(self):
        KratosMultiphysics.VariableUtils().SetHistoricalVariableToZero(KratosMultiphysics.NORMAL, self.model_part.Nodes)

    def test_ComputeSimplexNormalModelPart(self):
        KratosMultiphysics.NormalCalculationUtils().CalculateOnSimplex(self.model_part, 3)

        ## DEBUG
        #PostProcess(self.model_part)

        for node in self.model_part.GetSubModelPart("Skin_Part").Nodes:
            normal = CalculateAnalyticalNormal(node)

            solution_normal = node.GetSolutionStepValue(KratosMultiphysics.NORMAL)
            solution_normal_norm = math.sqrt(solution_normal[0]**2+solution_normal[1]**2+solution_normal[2]**2)
            solution_normal /= solution_normal_norm

            self.assertLess(CalculateNormalResidual(normal, solution_normal), 0.15)

    def test_ComputeUnitNormalModelPart(self):
        KratosMultiphysics.NormalCalculationUtils().CalculateUnitNormals(self.model_part)

        ## DEBUG
        #PostProcess(self.model_part)

        for node in self.model_part.GetSubModelPart("Skin_Part").Nodes:
            normal = CalculateAnalyticalNormal(node)
            solution_normal = node.GetSolutionStepValue(KratosMultiphysics.NORMAL)
            self.assertLess(CalculateNormalResidual(normal, solution_normal), 0.15)

    def test_ComputeNodesMeanNormalModelPart(self):
        KratosMultiphysics.MortarUtilities.ComputeNodesMeanNormalModelPart(self.model_part, True)

        ## DEBUG
        #PostProcess(self.model_part)

        for node in self.model_part.GetSubModelPart("Skin_Part").Nodes:
            normal = CalculateAnalyticalNormal(node)
            solution_normal = node.GetSolutionStepValue(KratosMultiphysics.NORMAL)
            self.assertLess(CalculateNormalResidual(normal, solution_normal), 0.1)

    def test_InvertNormal(self):
        KratosMultiphysics.MortarUtilities.InvertNormal(self.model_part.Conditions)
        KratosMultiphysics.MortarUtilities.ComputeNodesMeanNormalModelPart(self.model_part, True)

        ## DEBUG
        #PostProcess(self.model_part)

        for node in self.model_part.GetSubModelPart("Skin_Part").Nodes:
            normal = CalculateAnalyticalNormal(node)
            solution_normal = node.GetSolutionStepValue(KratosMultiphysics.NORMAL) * -1.0
            self.assertLess(CalculateNormalResidual(normal, solution_normal), 0.1)


class TestNormalUtilsQuadSphere(KratosUnittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.current_model = KratosMultiphysics.Model()
        cls.model_part = cls.current_model.CreateModelPart("Main")
        cls.model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] = 3
        cls.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL)
        model_part_io = KratosMultiphysics.ModelPartIO(GetFilePath("auxiliar_files_for_python_unittest/mdpa_files/quad_sphere"))
        model_part_io.ReadModelPart(cls.model_part)

    def setUp(self):
        KratosMultiphysics.VariableUtils().SetHistoricalVariableToZero(KratosMultiphysics.NORMAL, self.model_part.Nodes)

    def test_ComputeUnitNormalQuadModelPart(self):
        KratosMultiphysics.NormalCalculationUtils().CalculateUnitNormals(self.model_part)

        ## DEBUG
        #PostProcess(self.model_part)

        for node in self.model_part.Nodes:
            normal = CalculateAnalyticalNormal(node)
            solution_normal = node.GetSolutionStepValue(KratosMultiphysics.NORMAL)
            self.assertLess(CalculateNormalResidual(normal, solution_normal), 0.15)

if __name__ == '__main__':
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    KratosUnittest.main()


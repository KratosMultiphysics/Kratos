from __future__ import print_function, absolute_import, division
import KratosMultiphysics
import numpy as np
import os

from .MainKratosROM import TestStructuralMechanicsStaticROM
from .MainKratosHROM import TestStructuralMechanicsStaticHROM
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utilities


class ROMStaticStruct(KratosUnittest.TestCase):
#########################################################################################

    def test_Struct_Static_ROM_2D(self):
        with KratosUnittest.WorkFolderScope(".", __file__):
            with open("ProjectParametersROM.json",'r') as parameter_file:
                parameters = KratosMultiphysics.Parameters(parameter_file.read())
            model = KratosMultiphysics.Model()
            Simulation = TestStructuralMechanicsStaticROM(model,parameters)
            Simulation.Run()
            ObtainedOutput = Simulation.EvaluateQuantityOfInterest()
            ExpectedOutput = np.load('ExpectedOutput.npy')
            NodalArea = Simulation.EvaluateQuantityOfInterest2()

            UP=0
            DOWN=0
            for i in range (len(ObtainedOutput)):
                if ExpectedOutput[i] != 0:
                    UP += (NodalArea[i]*(   (1  - (ObtainedOutput[i] / ExpectedOutput[i] )    )**2)  )
                    DOWN +=  NodalArea[i]
            L2 = (np.sqrt(UP/DOWN)) *100
            self.assertLess(L2, 1.2e-10)#percent
            # Cleaning
            kratos_utilities.DeleteDirectoryIfExisting("__pycache__")
            kratos_utilities.DeleteDirectoryIfExisting("vtk_output")
            for file_name in os.listdir(os.getcwd()):
                if file_name.endswith(".bin") or file_name.endswith(".lst") :
                    kratos_utilities.DeleteFileIfExisting(file_name)

    def test_Struct_Static_HROM_2D(self):
        with KratosUnittest.WorkFolderScope(".", __file__):
            with open("ProjectParametersHROM.json",'r') as parameter_file:
                parameters = KratosMultiphysics.Parameters(parameter_file.read())
            model = KratosMultiphysics.Model()
            Simulation = TestStructuralMechanicsStaticHROM(model,parameters)
            Simulation.Run()
            computing_model_part = Simulation._solver.GetComputingModelPart()
            dimension = Simulation._GetSolver().settings["domain_size"].GetInt()
            area_calculator = KratosMultiphysics.CalculateNodalAreaProcess(computing_model_part , dimension)
            area_calculator.Execute()
            ExpectedOutput = np.load('ExpectedOutput.npy')
            UP=0
            DOWN=0
            for node in computing_model_part.Nodes:
                NodalArea = node.GetSolutionStepValue(KratosMultiphysics.NODAL_AREA)
                if (ExpectedOutput[(2*node.Id)-1] != 0) and (ExpectedOutput[(2*node.Id)-2] != 0):
                    UP += NodalArea*(    (ExpectedOutput[(2*node.Id)-2] - node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X, 0)  )**2)
                    UP += NodalArea*(    (ExpectedOutput[(2*node.Id)-1] - node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y, 0)  )**2)
                    DOWN +=  2*NodalArea
            L2 = (np.sqrt(UP/DOWN))*100
            self.assertLess(L2, 7.2e-8) #percent
            # Cleaning
            kratos_utilities.DeleteDirectoryIfExisting("__pycache__")
            kratos_utilities.DeleteDirectoryIfExisting("vtk_output")
            for file_name in os.listdir(os.getcwd()):
                if file_name.endswith(".bin") or file_name.endswith(".lst") :
                    kratos_utilities.DeleteFileIfExisting(file_name)


##########################################################################################

if __name__ == '__main__':
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    KratosUnittest.main()

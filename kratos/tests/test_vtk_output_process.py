from __future__ import print_function, absolute_import, division

import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest

import KratosMultiphysics.kratos_utilities as kratos_utils
import vtk_output_process

import os

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

class TestVtkOutputProcess(KratosUnittest.TestCase):
    def __SetupModelPart(self, model, model_part_name):
        self.mp = model.CreateModelPart(model_part_name)
        self.mp.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        self.mp.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        #create nodes
        self.mp.CreateNewNode(1, 0.00000, 1.00000, 0.00000)
        self.mp.CreateNewNode(2, 0.00000, 0.50000, 0.00000)
        self.mp.CreateNewNode(3, 0.50000, 1.00000, 0.00000)
        self.mp.CreateNewNode(4, 0.50000, 0.50000, 0.00000)
        self.mp.CreateNewNode(5, 0.00000, 0.00000, 0.00000)
        self.mp.CreateNewNode(6, 1.00000, 1.00000, 0.00000)
        self.mp.CreateNewNode(7, 1.00000, 0.50000, 0.00000)
        self.mp.CreateNewNode(8, 0.50000, 0.00000, 0.00000)
        self.mp.CreateNewNode(9, 1.00000, 0.00000, 0.00000)
        self.mp.CreateNewNode(10, 1.50000, 1.00000, 0.00000)
        self.mp.CreateNewNode(11, 1.50000, 0.50000, 0.00000)
        self.mp.CreateNewNode(12, 1.50000, 0.00000, 0.00000)
        self.mp.CreateNewNode(13, 2.00000, 1.00000, 0.00000)
        self.mp.CreateNewNode(14, 2.00000, 0.50000, 0.00000)
        self.mp.CreateNewNode(15, 2.00000, 0.00000, 0.00000)
        self.mp.CreateNewNode(16, 1.00000, 1.00000, 0.00000)
        self.mp.CreateNewNode(17, 1.00000, 0.50000, 0.00000)
        self.mp.CreateNewNode(18, 1.00000, 0.00000, 0.00000)

        #create Element
        self.mp.CreateNewElement("Element2D4N", 1,
                            [14, 11, 12, 15], self.mp.GetProperties()[1])
        self.mp.CreateNewElement("Element2D4N", 2,
                            [13, 10, 11, 14], self.mp.GetProperties()[1])
        self.mp.CreateNewElement("Element2D4N", 3,
                            [11, 17, 18, 12], self.mp.GetProperties()[1])
        self.mp.CreateNewElement("Element2D4N", 4,
                            [10, 16, 17, 11], self.mp.GetProperties()[1])
        self.mp.CreateNewElement("Element2D4N", 5,
                            [2, 4, 3, 1], self.mp.GetProperties()[1])
        self.mp.CreateNewElement("Element2D4N", 6, [5, 8, 4, 2],
                            self.mp.GetProperties()[1])
        self.mp.CreateNewElement("Element2D4N", 7, [4, 7, 6, 3],
                            self.mp.GetProperties()[1])
        self.mp.CreateNewElement("Element2D4N", 8, [8, 9, 7, 4],
                            self.mp.GetProperties()[1])

        self.mp.CreateNewCondition("LineCondition2D2N", 1, [1,2], self.mp.GetProperties()[1])
        self.mp.CreateNewCondition("LineCondition2D2N", 2, [2,5], self.mp.GetProperties()[1])
        self.mp.CreateNewCondition("LineCondition2D2N", 3, [13,14], self.mp.GetProperties()[1])
        self.mp.CreateNewCondition("LineCondition2D2N", 4, [14,15], self.mp.GetProperties()[1])

        #create a submodelpart for boundary conditions
        bcs = self.mp.CreateSubModelPart("FixedEdgeNodes")
        bcs.AddNodes([1, 2, 5])
        bcs.AddConditions([1,2])

        bcmn = self.mp.CreateSubModelPart("MovingNodes")
        bcmn.AddNodes([13, 14, 15])
        bcs.AddConditions([3,4])

    def __SetSolution(self):
        for node in self.mp.Nodes:
            node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT,0,[node.X,node.Y,node.Z])
            node.SetSolutionStepValue(KratosMultiphysics.VELOCITY,0,[2*node.X,2*node.Y,2*node.Z])

    def __SetupVtkOutputProcess(self, current_model, parameters):
        return vtk_output_process.Factory(parameters, current_model)

    def test_vtk_io(self):
        current_model = KratosMultiphysics.Model()
        model_part_name = "Main"
        self.__SetupModelPart(current_model, model_part_name)
        self.__SetSolution()

        vtk_output_parameters = KratosMultiphysics.Parameters("""
        {
            "model_part_name"                    : "Main",
            "file_format"                        : "ascii",
            "output_control_type"                : "step",
            "output_frequency"                   : 1.0,
            "output_sub_model_parts"             : true,
            "folder_name"                        : "test_output",
            "save_output_files_in_folder"        : true,
            "nodal_solution_step_data_variables" : ["DISPLACEMENT", "VELOCITY"],
            "nodal_data_value_variables"         : [],
            "element_data_value_variables"       : []
        }
        """)
        vtk_output_process = self.__SetupVtkOutputProcess(current_model, vtk_output_parameters)

        time = 0.0
        dt = 0.2
        step = 0
        end_time = 1.0
        vtk_output_process.ExecuteInitialize()

        while (time <= end_time):
            time = time + dt
            step = step + 1
            #print("STEP :: ", step, ", TIME :: ", time)
            vtk_output_process.ExecuteInitializeSolutionStep()
            self.mp.CloneTimeStep(time)
            vtk_output_process.ExecuteFinalizeSolutionStep()
            vtk_output_process.PrintOutput()

        vtk_output_process.ExecuteFinalize()

    def tearDown(self):
        kratos_utils.DeleteDirectoryIfExisting("test_output")

if __name__ == '__main__':
    KratosUnittest.main()
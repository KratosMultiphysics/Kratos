from __future__ import print_function, absolute_import, division
import KratosMultiphysics

import KratosMultiphysics.KratosUnittest as KratosUnittest
import os

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

expected_string = '''-Main- model part
    Buffer Size : 1
    Number of tables : 1
    Number of sub model parts : 2
    Current solution step index : 0

    Mesh 0 : 
        Number of Nodes       : 6
        Number of Properties  : 1
        Number of Elements    : 4
        Number of Conditions  : 5
        Number of Constraints : 0

    -Outlet- model part
        Number of tables : 0
        Number of sub model parts : 0

        Mesh 0 : 
            Number of Nodes       : 0
            Number of Properties  : 1
            Number of Elements    : 0
            Number of Conditions  : 1
            Number of Constraints : 0
    -Inlets- model part
        Number of tables : 1
        Number of sub model parts : 2

        Mesh 0 : 
            Number of Nodes       : 3
            Number of Properties  : 0
            Number of Elements    : 1
            Number of Conditions  : 3
            Number of Constraints : 0
        -Inlet2- model part
            Number of tables : 0
            Number of sub model parts : 0

            Mesh 0 : 
                Number of Nodes       : 0
                Number of Properties  : 0
                Number of Elements    : 0
                Number of Conditions  : 2
                Number of Constraints : 0
        -Inlet1- model part
            Number of tables : 0
            Number of sub model parts : 0

            Mesh 0 : 
                Number of Nodes       : 2
                Number of Properties  : 0
                Number of Elements    : 0
                Number of Conditions  : 2
                Number of Constraints : 0
'''
class TestObjectPrinting(KratosUnittest.TestCase):
    def test_Properties_printing(self):
        pass
    def test_ModelPart_printing(self):
        model_part = KratosMultiphysics.ModelPart("Main")
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VISCOSITY)
        model_part_io = KratosMultiphysics.ModelPartIO(GetFilePath("test_model_part_io_read"))
        model_part_io.ReadModelPart(model_part)

        # print(model_part)



        self.assertMultiLineEqual(str(model_part), expected_string)




if __name__ == '__main__':
    KratosUnittest.main()


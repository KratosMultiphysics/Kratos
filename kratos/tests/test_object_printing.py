from __future__ import print_function, absolute_import, division
import KratosMultiphysics

import KratosMultiphysics.KratosUnittest as KratosUnittest
import os

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

# The expected definitions are here to make the handling of the
# multiline-stings easier (no need to deal with indentation)
prop_str = '''Properties
    VISCOSITY : 5.3
    DENSITY : -95.3
This properties contains 0 tables'''

model_part_str = '''-Main- model part
    AMBIENT_TEMPERATURE : 250
    Buffer Size : 1
    Number of tables : 1
    Number of sub model parts : 2
    Current solution step index : 0

    Number of Geometries  : 0
    Mesh 0 :
        Number of Nodes       : 6
        Number of Properties  : 1
        Number of Elements    : 4
        Number of Conditions  : 5
        Number of Constraints : 0

    -Inlets- model part
        Number of tables : 1
        Number of sub model parts : 2

        Number of Geometries  : 0
        Mesh 0 :
            Number of Nodes       : 3
            Number of Properties  : 0
            Number of Elements    : 1
            Number of Conditions  : 3
            Number of Constraints : 0
        -Inlet1- model part
            Number of tables : 0
            Number of sub model parts : 0

            Number of Geometries  : 0
            Mesh 0 :
                Number of Nodes       : 2
                Number of Properties  : 0
                Number of Elements    : 0
                Number of Conditions  : 2
                Number of Constraints : 0
        -Inlet2- model part
            Number of tables : 0
            Number of sub model parts : 0

            Number of Geometries  : 0
            Mesh 0 :
                Number of Nodes       : 0
                Number of Properties  : 0
                Number of Elements    : 0
                Number of Conditions  : 2
                Number of Constraints : 0
    -Outlet- model part
        Number of tables : 0
        Number of sub model parts : 0

        Number of Geometries  : 0
        Mesh 0 :
            Number of Nodes       : 0
            Number of Properties  : 1
            Number of Elements    : 0
            Number of Conditions  : 1
            Number of Constraints : 0
'''

class TestObjectPrinting(KratosUnittest.TestCase):
    maxDiff = None # to display all the diff

    def test_Properties_str(self):
        prop = KratosMultiphysics.Properties(5)

        prop[KratosMultiphysics.VISCOSITY] = 5.3
        prop[KratosMultiphysics.DENSITY] = -95.3

        self.assertMultiLineEqual(str(prop), prop_str)

    def test_ModelPart_str(self):
        current_model = KratosMultiphysics.Model()
        model_part = current_model.CreateModelPart("Main")
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VISCOSITY)
        model_part_io = KratosMultiphysics.ModelPartIO(GetFilePath("auxiliar_files_for_python_unittest/mdpa_files/test_model_part_io_read"))
        model_part_io.ReadModelPart(model_part)

        self.assertMultiLineEqual(str(model_part), model_part_str)


if __name__ == '__main__':
    KratosUnittest.main()


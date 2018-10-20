from __future__ import print_function, absolute_import, division

import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest

import KratosMultiphysics.kratos_utilities as kratos_utils

try:
    import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
    missing_external_dependencies = False
    missing_application = ''
except ImportError as e:
    missing_external_dependencies = True
    # extract name of the missing application from the error message
    import re
    missing_application = re.search(r'''.*'KratosMultiphysics\.(.*)'.*''',
                                    '{0}'.format(e)).group(1)
import os
import sys

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)


class TestModelPartIO(KratosUnittest.TestCase):

    def setUp(self):
        if (sys.version_info < (3, 2)):
            self.assertRaisesRegex = self.assertRaisesRegexp

    def tearDown(self):
        # Clean up temporary files
        kratos_utils.DeleteFileIfExisting(GetFilePath("test_model_part_io_write.out.mdpa"))
        kratos_utils.DeleteFileIfExisting(GetFilePath("test_model_part_io_write.out.time"))
        kratos_utils.DeleteFileIfExisting(GetFilePath("test_model_part_io_write.time"))

    def test_model_part_io_read_model_part(self):
        current_model = KratosMultiphysics.Model()

        model_part = current_model.CreateModelPart("Main")
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VISCOSITY)
        model_part_io = KratosMultiphysics.ModelPartIO(GetFilePath("test_model_part_io_read"))
        model_part_io.ReadModelPart(model_part)

        self.assertEqual(model_part.NumberOfSubModelParts(), 2)

        self.assertEqual(model_part.NumberOfTables(), 1)
        self.assertEqual(model_part.NumberOfProperties(), 1)
        self.assertEqual(model_part.NumberOfNodes(), 6)
        self.assertEqual(model_part.NumberOfProperties(), 1)
        self.assertEqual(model_part.NumberOfElements(), 4)
        self.assertEqual(model_part.NumberOfConditions(), 5)

        self.assertTrue(model_part.GetNode(1).IsFixed(KratosMultiphysics.DISPLACEMENT_X))
        self.assertTrue(model_part.GetNode(2).IsFixed(KratosMultiphysics.DISPLACEMENT_X))
        self.assertFalse(model_part.GetNode(3).IsFixed(KratosMultiphysics.DISPLACEMENT_X))
        self.assertFalse(model_part.GetNode(972).IsFixed(KratosMultiphysics.DISPLACEMENT_X))
        self.assertTrue(model_part.GetNode(973).IsFixed(KratosMultiphysics.DISPLACEMENT_X))
        self.assertTrue(model_part.GetNode(974).IsFixed(KratosMultiphysics.DISPLACEMENT_X))

        self.assertTrue(model_part.GetNode(1).IsFixed(KratosMultiphysics.DISPLACEMENT_Y))
        self.assertTrue(model_part.GetNode(2).IsFixed(KratosMultiphysics.DISPLACEMENT_Y))
        self.assertFalse(model_part.GetNode(3).IsFixed(KratosMultiphysics.DISPLACEMENT_Y))
        self.assertFalse(model_part.GetNode(972).IsFixed(KratosMultiphysics.DISPLACEMENT_Y))
        self.assertTrue(model_part.GetNode(973).IsFixed(KratosMultiphysics.DISPLACEMENT_Y))
        self.assertTrue(model_part.GetNode(974).IsFixed(KratosMultiphysics.DISPLACEMENT_Y))

        self.assertTrue(model_part.GetNode(1).IsFixed(KratosMultiphysics.DISPLACEMENT_Z))
        self.assertTrue(model_part.GetNode(2).IsFixed(KratosMultiphysics.DISPLACEMENT_Z))
        self.assertFalse(model_part.GetNode(3).IsFixed(KratosMultiphysics.DISPLACEMENT_Z))
        self.assertFalse(model_part.GetNode(972).IsFixed(KratosMultiphysics.DISPLACEMENT_Z))
        self.assertTrue(model_part.GetNode(973).IsFixed(KratosMultiphysics.DISPLACEMENT_Z))
        self.assertTrue(model_part.GetNode(974).IsFixed(KratosMultiphysics.DISPLACEMENT_Z))

        self.assertEqual(model_part.GetNode(1).GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X),0.1)
        self.assertEqual(model_part.GetNode(2).GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X),0.2)
        self.assertEqual(model_part.GetNode(3).GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X),0.0)
        self.assertEqual(model_part.GetNode(972).GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X),0.0)
        self.assertEqual(model_part.GetNode(973).GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X),0.0)
        self.assertEqual(model_part.GetNode(974).GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X),0.0)

        self.assertEqual(model_part.GetNode(1).GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y),0.0)
        self.assertEqual(model_part.GetNode(2).GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y),0.0)
        self.assertEqual(model_part.GetNode(3).GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y),0.0)
        self.assertEqual(model_part.GetNode(972).GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y),0.0)
        self.assertEqual(model_part.GetNode(973).GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y),0.000973)
        self.assertEqual(model_part.GetNode(974).GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y),0.000974)

        self.assertEqual(model_part.GetNode(1).GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Z),0.0)
        self.assertEqual(model_part.GetNode(2).GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Z),0.0)
        self.assertEqual(model_part.GetNode(3).GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Z),0.0)
        self.assertEqual(model_part.GetNode(972).GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Z),0.0)
        self.assertEqual(model_part.GetNode(973).GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Z),0.0)
        self.assertEqual(model_part.GetNode(974).GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Z),0.0)

        self.assertEqual(model_part.GetNode(1).GetSolutionStepValue(KratosMultiphysics.VISCOSITY),0.01)
        self.assertEqual(model_part.GetNode(2).GetSolutionStepValue(KratosMultiphysics.VISCOSITY),0.01)
        self.assertEqual(model_part.GetNode(3).GetSolutionStepValue(KratosMultiphysics.VISCOSITY),0.0)
        self.assertEqual(model_part.GetNode(972).GetSolutionStepValue(KratosMultiphysics.VISCOSITY),0.0)
        self.assertEqual(model_part.GetNode(973).GetSolutionStepValue(KratosMultiphysics.VISCOSITY),0.01)
        self.assertEqual(model_part.GetNode(974).GetSolutionStepValue(KratosMultiphysics.VISCOSITY),0.01)

        self.assertTrue(model_part.HasSubModelPart("Inlets"))

        inlets_model_part = model_part.GetSubModelPart("Inlets")

        self.assertEqual(inlets_model_part.NumberOfTables(), 1)
        self.assertEqual(inlets_model_part.NumberOfProperties(), 0)
        self.assertEqual(inlets_model_part.NumberOfNodes(), 3)
        self.assertEqual(inlets_model_part.NumberOfElements(), 1)
        self.assertEqual(inlets_model_part.NumberOfConditions(), 3)
        self.assertEqual(inlets_model_part.NumberOfSubModelParts(), 2)
        self.assertTrue(inlets_model_part.HasSubModelPart("Inlet1"))
        self.assertTrue(inlets_model_part.HasSubModelPart("Inlet2"))

        inlet1_model_part = inlets_model_part.GetSubModelPart("Inlet1")

        self.assertEqual(inlet1_model_part.NumberOfTables(), 0)
        self.assertEqual(inlet1_model_part.NumberOfProperties(), 0)
        self.assertEqual(inlet1_model_part.NumberOfNodes(), 2)
        self.assertEqual(inlet1_model_part.NumberOfElements(), 0)
        self.assertEqual(inlet1_model_part.NumberOfConditions(), 2)
        self.assertEqual(inlet1_model_part.NumberOfSubModelParts(), 0)

        inlet2_model_part = inlets_model_part.GetSubModelPart("Inlet2")

        self.assertEqual(inlet2_model_part.NumberOfTables(), 0)
        self.assertEqual(inlet2_model_part.NumberOfProperties(), 0)
        self.assertEqual(inlet2_model_part.NumberOfNodes(), 0)
        self.assertEqual(inlet2_model_part.NumberOfElements(), 0)
        self.assertEqual(inlet2_model_part.NumberOfConditions(), 2)
        self.assertEqual(inlet2_model_part.NumberOfSubModelParts(), 0)

        self.assertTrue(model_part.HasSubModelPart("Outlet"))

        outlet_model_part = model_part.GetSubModelPart("Outlet")

        self.assertEqual(outlet_model_part.NumberOfTables(), 0)
        self.assertEqual(outlet_model_part.NumberOfProperties(), 1)
        self.assertEqual(outlet_model_part.NumberOfNodes(), 0)
        self.assertEqual(outlet_model_part.NumberOfElements(), 0)
        self.assertEqual(outlet_model_part.NumberOfConditions(), 1)
        self.assertEqual(outlet_model_part.NumberOfSubModelParts(), 0)

        properties_1 = model_part.GetProperties()[1]
        #Bools
        self.assertTrue(properties_1[KratosMultiphysics.IS_RESTARTED])
        self.assertFalse(properties_1[KratosMultiphysics.COMPUTE_DYNAMIC_TANGENT])
        #Double
        self.assertEqual(properties_1[KratosMultiphysics.DENSITY], 3.4E-5)
        #Array3
        self.assertEqual(properties_1[KratosMultiphysics.VOLUME_ACCELERATION][0], 0.00)
        self.assertEqual(properties_1[KratosMultiphysics.VOLUME_ACCELERATION][1], 0.00)
        self.assertEqual(properties_1[KratosMultiphysics.VOLUME_ACCELERATION][2], 9.8)
        #Matrix3x3
        self.assertEqual(properties_1[KratosMultiphysics.LOCAL_INERTIA_TENSOR][0,0], 0)
        self.assertEqual(properties_1[KratosMultiphysics.LOCAL_INERTIA_TENSOR][0,1], 0.27)
        self.assertEqual(properties_1[KratosMultiphysics.LOCAL_INERTIA_TENSOR][0,2], 0.27)
        self.assertEqual(properties_1[KratosMultiphysics.LOCAL_INERTIA_TENSOR][1,0], 0.087)
        self.assertEqual(properties_1[KratosMultiphysics.LOCAL_INERTIA_TENSOR][1,1], 0)
        self.assertEqual(properties_1[KratosMultiphysics.LOCAL_INERTIA_TENSOR][1,2], 0.27)
        self.assertEqual(properties_1[KratosMultiphysics.LOCAL_INERTIA_TENSOR][2,0], 0.075)
        self.assertEqual(properties_1[KratosMultiphysics.LOCAL_INERTIA_TENSOR][2,1], 0.23)
        self.assertEqual(properties_1[KratosMultiphysics.LOCAL_INERTIA_TENSOR][2,2], 0)

        #SubModelPartData
        self.assertTrue(inlets_model_part[KratosMultiphysics.IS_RESTARTED])
        self.assertFalse(inlets_model_part[KratosMultiphysics.COMPUTE_DYNAMIC_TANGENT])

    def test_model_part_io_read_model_part_mesh_only(self):
        current_model = KratosMultiphysics.Model()

        model_part = current_model.CreateModelPart("Main")
        model_part_io = KratosMultiphysics.ModelPartIO(GetFilePath("test_model_part_io_read"), KratosMultiphysics.IO.MESH_ONLY)
        model_part_io.ReadModelPart(model_part)

        self.assertEqual(model_part.NumberOfSubModelParts(), 2)

        self.assertEqual(model_part.NumberOfTables(), 0)
        self.assertEqual(model_part.NumberOfProperties(), 1)
        self.assertEqual(model_part.NumberOfNodes(), 6)
        self.assertEqual(model_part.NumberOfElements(), 4)
        self.assertEqual(model_part.NumberOfConditions(), 5)

        self.assertTrue(model_part.HasSubModelPart("Inlets"))

        inlets_model_part = model_part.GetSubModelPart("Inlets")

        self.assertEqual(inlets_model_part.NumberOfTables(), 0)
        self.assertEqual(inlets_model_part.NumberOfProperties(), 0)
        self.assertEqual(inlets_model_part.NumberOfNodes(), 3)
        self.assertEqual(inlets_model_part.NumberOfElements(), 1)
        self.assertEqual(inlets_model_part.NumberOfConditions(), 3)
        self.assertEqual(inlets_model_part.NumberOfSubModelParts(), 2)
        self.assertTrue(inlets_model_part.HasSubModelPart("Inlet1"))
        self.assertTrue(inlets_model_part.HasSubModelPart("Inlet2"))

        inlet1_model_part = inlets_model_part.GetSubModelPart("Inlet1")

        self.assertEqual(inlet1_model_part.NumberOfTables(), 0)
        self.assertEqual(inlet1_model_part.NumberOfProperties(), 0)
        self.assertEqual(inlet1_model_part.NumberOfNodes(), 2)
        self.assertEqual(inlet1_model_part.NumberOfElements(), 0)
        self.assertEqual(inlet1_model_part.NumberOfConditions(), 2)
        self.assertEqual(inlet1_model_part.NumberOfSubModelParts(), 0)

        inlet2_model_part = inlets_model_part.GetSubModelPart("Inlet2")

        self.assertEqual(inlet2_model_part.NumberOfTables(), 0)
        self.assertEqual(inlet2_model_part.NumberOfProperties(), 0)
        self.assertEqual(inlet2_model_part.NumberOfNodes(), 0)
        self.assertEqual(inlet2_model_part.NumberOfElements(), 0)
        self.assertEqual(inlet2_model_part.NumberOfConditions(), 2)
        self.assertEqual(inlet2_model_part.NumberOfSubModelParts(), 0)

        self.assertTrue(model_part.HasSubModelPart("Outlet"))

        outlet_model_part = model_part.GetSubModelPart("Outlet")

        self.assertEqual(outlet_model_part.NumberOfTables(), 0)
        self.assertEqual(outlet_model_part.NumberOfProperties(), 1)
        self.assertEqual(outlet_model_part.NumberOfNodes(), 0)
        self.assertEqual(outlet_model_part.NumberOfElements(), 0)
        self.assertEqual(outlet_model_part.NumberOfConditions(), 1)
        self.assertEqual(outlet_model_part.NumberOfSubModelParts(), 0)

        properties_1 = model_part.GetProperties()[1]
        #Bools
        self.assertTrue(properties_1[KratosMultiphysics.IS_RESTARTED])
        self.assertFalse(properties_1[KratosMultiphysics.COMPUTE_DYNAMIC_TANGENT])
        #Double
        self.assertEqual(properties_1[KratosMultiphysics.DENSITY], 3.4E-5)
        #Array3
        self.assertEqual(properties_1[KratosMultiphysics.VOLUME_ACCELERATION][0], 0.00)
        self.assertEqual(properties_1[KratosMultiphysics.VOLUME_ACCELERATION][1], 0.00)
        self.assertEqual(properties_1[KratosMultiphysics.VOLUME_ACCELERATION][2], 9.8)
        #Matrix3x3
        self.assertEqual(properties_1[KratosMultiphysics.LOCAL_INERTIA_TENSOR][0,0], 0)
        self.assertEqual(properties_1[KratosMultiphysics.LOCAL_INERTIA_TENSOR][0,1], 0.27)
        self.assertEqual(properties_1[KratosMultiphysics.LOCAL_INERTIA_TENSOR][0,2], 0.27)
        self.assertEqual(properties_1[KratosMultiphysics.LOCAL_INERTIA_TENSOR][1,0], 0.087)
        self.assertEqual(properties_1[KratosMultiphysics.LOCAL_INERTIA_TENSOR][1,1], 0)
        self.assertEqual(properties_1[KratosMultiphysics.LOCAL_INERTIA_TENSOR][1,2], 0.27)
        self.assertEqual(properties_1[KratosMultiphysics.LOCAL_INERTIA_TENSOR][2,0], 0.075)
        self.assertEqual(properties_1[KratosMultiphysics.LOCAL_INERTIA_TENSOR][2,1], 0.23)
        self.assertEqual(properties_1[KratosMultiphysics.LOCAL_INERTIA_TENSOR][2,2], 0)

    def test_model_part_io_write_model_part(self):
        if (missing_external_dependencies is False):
            current_model = KratosMultiphysics.Model()
            model_part = current_model.CreateModelPart("Main")
            model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
            model_part_io = KratosMultiphysics.ModelPartIO(GetFilePath("test_model_part_io_write"))
            model_part_io.ReadModelPart(model_part)

            model_part_io = KratosMultiphysics.ModelPartIO(GetFilePath("test_model_part_io_write.out"), KratosMultiphysics.IO.WRITE)
            model_part_io.WriteModelPart(model_part)

            import filecmp
            value = filecmp.cmp(GetFilePath("test_model_part_io_write.mdpa"), GetFilePath("test_model_part_io_write.out.mdpa"))
            self.assertEqual(value, True)
        else:
            KratosMultiphysics.Logger.PrintInfo("TestModelPartIO", "Please compile StructuralMechanicsApplication in order to test output in IO")

    @KratosUnittest.expectedFailure
    def test_error_on_wrong_input(self):
        current_model =  KratosMultiphysics.Model()

        model_part = current_model.CreateModelPart("Main")
        model_part_io = KratosMultiphysics.ModelPartIO(GetFilePath("wrong_properties_input"))

        #an error shall be thrown while reading the input since the format is not correct
        try:
            with self.assertRaisesRegex(RuntimeError, "wrong input format while reading Properties"): #ideally a more specific error message shall be devised
                pass #the real line shall be the one below but it segfaults badly
                #model_part_io.ReadModelPart(model_part)
        except:
            raise Exception("a segmentation fault is issued!!")
            self.fail("a segmentation fault is issued!!")




    #def test_model_part_io_properties_block(self):
    #    model_part= current_model.CreateModelPart("Main")
    #    model_part_io = ModelPartIO("test_model_part_io")
    #    model_part_io.ReadProperties(model_part.Properties)

if __name__ == '__main__':
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    KratosUnittest.main()

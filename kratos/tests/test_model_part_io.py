import os
from pathlib import Path

# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as KratosUtils
from KratosMultiphysics.testing.utilities import ReadModelPart

structural_mechanics_is_available = KratosUtils.CheckIfApplicationsAvailable("StructuralMechanicsApplication")
if structural_mechanics_is_available:
    import KratosMultiphysics.StructuralMechanicsApplication


def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)


class TestModelPartIO(KratosUnittest.TestCase):
    def tearDown(self):
        # Clean up temporary files
        KratosUtils.DeleteFileIfExisting(GetFilePath("test_model_part_io_write.out.mdpa"))
        KratosUtils.DeleteFileIfExisting(GetFilePath("test_model_part_io_write.out.time"))
        KratosUtils.DeleteFileIfExisting(GetFilePath("test_model_part_io_write.time"))
        KratosUtils.DeleteFileIfExisting(GetFilePath("test_model_part_io_write_mesh_only.out.mdpa"))
        KratosUtils.DeleteFileIfExisting(GetFilePath("test_model_part_io_write_mesh_only.out.time"))
        KratosUtils.DeleteFileIfExisting(GetFilePath("test_model_part_io_write_mesh_only.time"))

    def test_model_part_io_read_model_part(self):
        input_mdpa = GetFilePath("auxiliar_files_for_python_unittest/mdpa_files/test_model_part_io_read")

        with self.subTest("string"):
            self.execute_test_model_part_io_read_model_part(str(input_mdpa))

        with self.subTest("pathlib.Path"):
            self.execute_test_model_part_io_read_model_part(Path(input_mdpa))

    def test_model_part_io_read_using_submodelpart(self):
        current_model = KratosMultiphysics.Model()

        main_model_part = current_model.CreateModelPart("Main")
        model_part = main_model_part.CreateSubModelPart("submodelpart")
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VISCOSITY)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)

        model_part_io = KratosMultiphysics.ModelPartIO(GetFilePath("auxiliar_files_for_python_unittest/mdpa_files/test_model_part_io_read"))
        model_part_io.ReadModelPart(model_part)

        self.assertEqual(model_part.NumberOfSubModelParts(), 2)
        self.assertEqual(model_part.NumberOfTables(), 1)
        self.assertEqual(model_part.NumberOfProperties(), 1)
        self.assertEqual(model_part.NumberOfNodes(), 6)
        self.assertEqual(model_part.NumberOfProperties(), 1)
        self.assertEqual(model_part.NumberOfGeometries(), 9)
        self.assertEqual(model_part.NumberOfElements(), 4)
        self.assertEqual(model_part.NumberOfConditions(), 5)
        self.assertEqual(model_part.NumberOfMasterSlaveConstraints(), 2)

    def execute_test_model_part_io_read_model_part(self, input_mdpa):
        current_model = KratosMultiphysics.Model()

        model_part = current_model.CreateModelPart("Main")
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VISCOSITY)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)

        model_part_io = KratosMultiphysics.ModelPartIO(input_mdpa)
        model_part_io.ReadModelPart(model_part)

        self.assertEqual(model_part.NumberOfSubModelParts(), 2)

        self.assertEqual(model_part.NumberOfTables(), 1)
        self.assertEqual(model_part.NumberOfProperties(), 1)
        self.assertEqual(model_part.NumberOfNodes(), 6)
        self.assertEqual(model_part.NumberOfProperties(), 1)
        self.assertEqual(model_part.NumberOfGeometries(), 9)
        self.assertEqual(model_part.NumberOfElements(), 4)
        self.assertEqual(model_part.NumberOfConditions(), 5)

        self.assertEqual(model_part[KratosMultiphysics.AMBIENT_TEMPERATURE], 250.0)
        self.assertEqual(model_part[KratosMultiphysics.DISPLACEMENT_X], 2.1)
        self.assertEqual(model_part[KratosMultiphysics.DISPLACEMENT_Y], 3.2)
        self.assertEqual(model_part[KratosMultiphysics.DISPLACEMENT_Z], 4.3)
        self.assertEqual(model_part[KratosMultiphysics.VELOCITY_X], 3.8)
        self.assertEqual(model_part[KratosMultiphysics.VELOCITY_Y], 4.9)
        self.assertEqual(model_part[KratosMultiphysics.VELOCITY_Z], 0.0)

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

        self.assertEqual(model_part.GetNode(1).GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y), 0.0)
        self.assertEqual(model_part.GetNode(2).GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y), 0.0)
        self.assertEqual(model_part.GetNode(3).GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y), 0.0)
        self.assertEqual(model_part.GetNode(972).GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y), 0.0)
        self.assertEqual(model_part.GetNode(973).GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y), 0.000973)
        self.assertEqual(model_part.GetNode(974).GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y), 0.000974)

        self.assertEqual(model_part.GetNode(1).GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Z), 0.0)
        self.assertEqual(model_part.GetNode(2).GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Z), 0.0)
        self.assertEqual(model_part.GetNode(3).GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Z), 0.0)
        self.assertEqual(model_part.GetNode(972).GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Z), 0.0)
        self.assertEqual(model_part.GetNode(973).GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Z), 0.0)
        self.assertEqual(model_part.GetNode(974).GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Z), 0.0)

        self.assertEqual(model_part.GetNode(1).GetSolutionStepValue(KratosMultiphysics.VISCOSITY), 0.01)
        self.assertEqual(model_part.GetNode(2).GetSolutionStepValue(KratosMultiphysics.VISCOSITY), 0.01)
        self.assertEqual(model_part.GetNode(3).GetSolutionStepValue(KratosMultiphysics.VISCOSITY), 0.0)
        self.assertEqual(model_part.GetNode(972).GetSolutionStepValue(KratosMultiphysics.VISCOSITY), 0.0)
        self.assertEqual(model_part.GetNode(973).GetSolutionStepValue(KratosMultiphysics.VISCOSITY), 0.01)
        self.assertEqual(model_part.GetNode(974).GetSolutionStepValue(KratosMultiphysics.VISCOSITY), 0.01)

        self.assertEqual(model_part.GetNode(1).GetSolutionStepValue(KratosMultiphysics.VELOCITY_X), 1.1)
        self.assertEqual(model_part.GetNode(1).GetSolutionStepValue(KratosMultiphysics.VELOCITY_Y), 2.2)
        self.assertEqual(model_part.GetNode(1).GetSolutionStepValue(KratosMultiphysics.VELOCITY_Z), 3.3)

        self.assertEqual(model_part.GetNode(1).Is(KratosMultiphysics.BOUNDARY), True)
        self.assertEqual(model_part.GetNode(2).Is(KratosMultiphysics.BOUNDARY), True)
        self.assertEqual(model_part.GetNode(3).Is(KratosMultiphysics.BOUNDARY), False)
        self.assertEqual(model_part.GetNode(972).Is(KratosMultiphysics.BOUNDARY), False)
        self.assertEqual(model_part.GetNode(973).Is(KratosMultiphysics.BOUNDARY), True)
        self.assertEqual(model_part.GetNode(974).Is(KratosMultiphysics.BOUNDARY), True)

        self.assertTrue(model_part.HasSubModelPart("Inlets"))

        inlets_model_part = model_part.GetSubModelPart("Inlets")

        self.assertEqual(inlets_model_part.NumberOfTables(), 1)
        self.assertEqual(inlets_model_part.NumberOfProperties(), 0)
        self.assertEqual(inlets_model_part.NumberOfNodes(), 3)
        self.assertEqual(inlets_model_part.NumberOfGeometries(), 0)
        self.assertEqual(inlets_model_part.NumberOfElements(), 1)
        self.assertEqual(inlets_model_part.NumberOfConditions(), 3)
        self.assertEqual(inlets_model_part.NumberOfSubModelParts(), 2)
        self.assertTrue(inlets_model_part.HasSubModelPart("Inlet1"))
        self.assertTrue(inlets_model_part.HasSubModelPart("Inlet2"))

        inlet1_model_part = inlets_model_part.GetSubModelPart("Inlet1")

        self.assertEqual(inlet1_model_part.NumberOfTables(), 0)
        self.assertEqual(inlet1_model_part.NumberOfProperties(), 0)
        self.assertEqual(inlet1_model_part.NumberOfNodes(), 2)
        self.assertEqual(inlet1_model_part.NumberOfGeometries(), 0)
        self.assertEqual(inlet1_model_part.NumberOfElements(), 0)
        self.assertEqual(inlet1_model_part.NumberOfConditions(), 2)
        self.assertEqual(inlet1_model_part.NumberOfSubModelParts(), 0)

        inlet2_model_part = inlets_model_part.GetSubModelPart("Inlet2")

        self.assertEqual(inlet2_model_part.NumberOfTables(), 0)
        self.assertEqual(inlet2_model_part.NumberOfProperties(), 0)
        self.assertEqual(inlet2_model_part.NumberOfNodes(), 0)
        self.assertEqual(inlet2_model_part.NumberOfGeometries(), 0)
        self.assertEqual(inlet2_model_part.NumberOfElements(), 0)
        self.assertEqual(inlet2_model_part.NumberOfConditions(), 2)
        self.assertEqual(inlet2_model_part.NumberOfSubModelParts(), 0)

        self.assertTrue(model_part.HasSubModelPart("Outlet"))

        outlet_model_part = model_part.GetSubModelPart("Outlet")

        self.assertEqual(outlet_model_part.NumberOfTables(), 0)
        self.assertEqual(outlet_model_part.NumberOfProperties(), 1)
        self.assertEqual(outlet_model_part.NumberOfNodes(), 0)
        self.assertEqual(outlet_model_part.NumberOfGeometries(), 0)
        self.assertEqual(outlet_model_part.NumberOfElements(), 0)
        self.assertEqual(outlet_model_part.NumberOfConditions(), 1)
        self.assertEqual(outlet_model_part.NumberOfSubModelParts(), 0)

        properties_1 = model_part.GetProperties()[1]
        #Bools
        self.assertTrue(properties_1[KratosMultiphysics.IS_RESTARTED])
        self.assertFalse(properties_1[KratosMultiphysics.COMPUTE_DYNAMIC_TANGENT])
        self.assertFalse(properties_1[KratosMultiphysics.COMPUTE_LUMPED_MASS_MATRIX])
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
        self.assertTrue(inlets_model_part[KratosMultiphysics.COMPUTE_LUMPED_MASS_MATRIX])
        self.assertFalse(inlets_model_part[KratosMultiphysics.COMPUTE_DYNAMIC_TANGENT])

    def test_model_part_io_read_model_part_mesh_only(self):
        current_model = KratosMultiphysics.Model()

        model_part = current_model.CreateModelPart("Main")
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        model_part_io = KratosMultiphysics.ModelPartIO(GetFilePath("auxiliar_files_for_python_unittest/mdpa_files/test_model_part_io_read"), KratosMultiphysics.IO.MESH_ONLY)
        model_part_io.ReadModelPart(model_part)

        self.assertEqual(model_part.NumberOfSubModelParts(), 2)

        self.assertEqual(model_part.NumberOfTables(), 0)
        self.assertEqual(model_part.NumberOfProperties(), 1)
        self.assertEqual(model_part.NumberOfNodes(), 6)
        self.assertEqual(model_part.NumberOfGeometries(), 9)
        self.assertEqual(model_part.NumberOfElements(), 4)
        self.assertEqual(model_part.NumberOfConditions(), 5)
        self.assertEqual(model_part.NumberOfMasterSlaveConstraints(), 2)

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

    @KratosUnittest.skipUnless(structural_mechanics_is_available,"StructuralMechanicsApplication is not available")
    def test_model_part_io_write_model_part(self):
        current_model = KratosMultiphysics.Model()
        model_part = current_model.CreateModelPart("Main")
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        model_part_io = KratosMultiphysics.ModelPartIO(GetFilePath("auxiliar_files_for_python_unittest/mdpa_files/test_model_part_io_write"))
        model_part_io.ReadModelPart(model_part)

        model_part_io = KratosMultiphysics.ModelPartIO(GetFilePath("test_model_part_io_write.out"), KratosMultiphysics.IO.WRITE | KratosMultiphysics.IO.SCIENTIFIC_PRECISION)
        model_part_io.WriteModelPart(model_part)

        import filecmp
        value = filecmp.cmp(GetFilePath("auxiliar_files_for_python_unittest/mdpa_files/test_model_part_io_write.mdpa"), GetFilePath("test_model_part_io_write.out.mdpa"))
        self.assertTrue(value)

        # Writing a second time should yield the same results.
        # (unless writing is allowed to mutate the model part because we're not const-correct ...)
        KratosUtils.DeleteFileIfExisting(GetFilePath("test_model_part_io_write.out.mdpa"))
        KratosUtils.DeleteFileIfExisting(GetFilePath("test_model_part_io_write.out.time"))
        model_part_io.WriteModelPart(model_part)
        value = filecmp.cmp(GetFilePath("auxiliar_files_for_python_unittest/mdpa_files/test_model_part_io_write.mdpa"), GetFilePath("test_model_part_io_write.out.mdpa"))
        self.assertTrue(value)

    @KratosUnittest.skipUnless(structural_mechanics_is_available,"StructuralMechanicsApplication is not available")
    def test_model_part_io_write_model_part_mesh_only(self):
        current_model = KratosMultiphysics.Model()
        model_part = current_model.CreateModelPart("Main")
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        model_part_io = KratosMultiphysics.ModelPartIO(GetFilePath("auxiliar_files_for_python_unittest/mdpa_files/test_model_part_io_write"))
        model_part_io.ReadModelPart(model_part)

        model_part_io = KratosMultiphysics.ModelPartIO(GetFilePath("test_model_part_io_write_mesh_only.out"), KratosMultiphysics.IO.WRITE | KratosMultiphysics.IO.MESH_ONLY)
        model_part_io.WriteModelPart(model_part)

        import filecmp
        value = filecmp.cmp(GetFilePath("auxiliar_files_for_python_unittest/mdpa_files/test_model_part_io_write_mesh_only.mdpa"), GetFilePath("test_model_part_io_write_mesh_only.out.mdpa"))
        self.assertTrue(value)

    @KratosUnittest.expectedFailure
    def test_error_on_wrong_input(self):
        current_model =  KratosMultiphysics.Model()

        model_part = current_model.CreateModelPart("Main")
        model_part_io = KratosMultiphysics.ModelPartIO(GetFilePath("auxiliar_files_for_python_unittest/mdpa_files/wrong_properties_input"))

        #an error shall be thrown while reading the input since the format is not correct
        try:
            with self.assertRaisesRegex(RuntimeError, "wrong input format while reading Properties"): #ideally a more specific error message shall be devised
                pass #the real line shall be the one below but it segfaults badly
                #model_part_io.ReadModelPart(model_part)
        except:
            raise Exception("a segmentation fault is issued!!")
            self.fail("a segmentation fault is issued!!")

class TestModelPartIOMPI(KratosUnittest.TestCase):
    def test_model_part_io_read_entity_data(self):
        # testing if the assignment of entity data works correctly in serial and MPI
        current_model = KratosMultiphysics.Model()

        model_part = current_model.CreateModelPart("Main")
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.BULK_MODULUS)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_VAUX)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.EXTERNAL_FORCES_VECTOR)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.LOCAL_AXES_MATRIX)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PARTITION_INDEX)
        model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] = 3

        ReadModelPart(GetFilePath("auxiliar_files_for_python_unittest/mdpa_files/coarse_sphere_with_conditions"), model_part)

        def GetScalar(entity_id):
            return entity_id**2

        def GetArray3(entity_id):
            return [entity_id+1.5, entity_id-23, entity_id*2]

        def GetVector(entity_id):
            vec = [entity_id+1, entity_id-23, entity_id*2]
            if entity_id > 10:
                vec.append(entity_id*1.5)
            if entity_id > 25:
                vec.append(entity_id**2)
            if entity_id > 60:
                vec.append(entity_id/2.5)
            return vec

        def GetMatrix(entity_id):
            mat = KratosMultiphysics.Matrix(2,3)
            for i in range(mat.Size1()):
                for j in range(mat.Size2()):
                    mat[i,j] = i+j+entity_id
            return mat

        for node in model_part.Nodes:
            if node.Id in [1,10,55,81]:
                self.assertAlmostEqual(node.GetSolutionStepValue(KratosMultiphysics.BULK_MODULUS), GetScalar(node.Id))
                self.assertVectorAlmostEqual(node.GetSolutionStepValue(KratosMultiphysics.NODAL_VAUX), GetArray3(node.Id))
                self.assertVectorAlmostEqual(node.GetSolutionStepValue(KratosMultiphysics.EXTERNAL_FORCES_VECTOR), GetVector(node.Id))
                self.assertMatrixAlmostEqual(node.GetSolutionStepValue(KratosMultiphysics.LOCAL_AXES_MATRIX), GetMatrix(node.Id))

        def CheckEntities(entites, ids_to_check):
            for ent in entites:
                if ent.Id in ids_to_check:
                    self.assertAlmostEqual(ent.GetValue(KratosMultiphysics.TEMPERATURE), GetScalar(ent.Id))
                    self.assertVectorAlmostEqual(ent.GetValue(KratosMultiphysics.MESH_VELOCITY), GetArray3(ent.Id))
                    self.assertVectorAlmostEqual(ent.GetValue(KratosMultiphysics.INITIAL_STRAIN), GetVector(ent.Id))
                    self.assertMatrixAlmostEqual(ent.GetValue(KratosMultiphysics.LOCAL_INERTIA_TENSOR), GetMatrix(ent.Id))

        elem_ids_to_check = [5,64,33,214]
        cond_ids_to_check = [2,13,22,121]

        CheckEntities(model_part.Elements, elem_ids_to_check)
        CheckEntities(model_part.Conditions, cond_ids_to_check)


if __name__ == '__main__':
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    KratosUnittest.main()

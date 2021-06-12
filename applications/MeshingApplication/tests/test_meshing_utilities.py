# We import the libraries
import KratosMultiphysics as Kratos
import KratosMultiphysics.MeshingApplication as MeshingApplication

import KratosMultiphysics.KratosUnittest as KratosUnittest

from KratosMultiphysics.kratos_utilities import DeleteFileIfExisting
from KratosMultiphysics.kratos_utilities import DeleteTimeFiles

class TestMPIMetrics(KratosUnittest.TestCase):
    @classmethod
    def tearDownClass(cls):
        DeleteTimeFiles("element_refinement_utilities")

    def test_RefineElement2D(self):
        self.addCleanup(lambda: DeleteFileIfExisting("element_refinement_utilities/2d_refined_model_part.mdpa"))

        model = Kratos.Model()
        model_part = model.CreateModelPart("test")
        model_part.ProcessInfo[Kratos.DOMAIN_SIZE] = 2

        model_part.AddNodalSolutionStepVariable(Kratos.PRESSURE)
        model_part.AddNodalSolutionStepVariable(Kratos.VELOCITY)
        model_part.AddNodalSolutionStepVariable(Kratos.NORMAL)

        node_1 = model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
        node_2 = model_part.CreateNewNode(2, 0.0, 1.0, 0.0)
        node_3 = model_part.CreateNewNode(3, 1.0, 1.0, 0.0)

        properties = model_part.GetProperties()[0]
        element = model_part.CreateNewElement("Element2D3N", 1, [1, 2, 3], properties)

        p_condition_1 = model_part.CreateNewCondition("LineCondition2D2N", 1, [1, 2], properties)
        p_condition_2 = model_part.CreateNewCondition("LineCondition2D2N", 2, [2, 3], properties)
        p_condition_3 = model_part.CreateNewCondition("LineCondition2D2N", 3, [3, 1], properties)
        p_condition_2.SetValue(Kratos.DISTANCE, 10)

        p_condition_1.Set(Kratos.INLET)
        node_1.Set(Kratos.INLET)
        node_2.Set(Kratos.INLET)
        p_condition_2.Set(Kratos.OUTLET)
        node_2.Set(Kratos.OUTLET)
        node_3.Set(Kratos.OUTLET)
        p_condition_3.Set(Kratos.SLIP)
        node_3.Set(Kratos.SLIP)
        node_1.Set(Kratos.SLIP)

        tmoc = Kratos.TetrahedralMeshOrientationCheck
        throw_errors = False
        flags = tmoc.COMPUTE_NODAL_NORMALS | tmoc.COMPUTE_CONDITION_NORMALS | tmoc.ASSIGN_NEIGHBOUR_ELEMENTS_TO_CONDITIONS
        Kratos.TetrahedralMeshOrientationCheck(model_part, throw_errors, flags).Execute()

        element_refinement_utilities = MeshingApplication.ElementRefinementUtilities(
            model_part,
            "refined_model_part",
            element,
            2,
            "Element2D3N",
            "LineCondition2D2N",
            ["VELOCITY_X", "PRESSURE"],
            ["RESIDUAL_NORM"],
            ["INLET", "OUTLET", "SLIP"])

        refined_model_part = element_refinement_utilities.GetRefinedModelPart()

        node_1.X += 1.0
        node_2.X += 1.0
        node_3.X += 1.0

        node_1.Y += 1.0
        node_2.Y += 1.0
        node_3.Y += 1.0

        node_1.SetSolutionStepValue(Kratos.VELOCITY, Kratos.Array3([1.0, 1.0, 1.0]))
        node_2.SetSolutionStepValue(Kratos.VELOCITY, Kratos.Array3([10.0, 10.0, 10.0]))
        node_3.SetSolutionStepValue(Kratos.VELOCITY, Kratos.Array3([20.0, 20.0, 20.0]))

        node_1.SetSolutionStepValue(Kratos.PRESSURE, 100.0)
        node_2.SetSolutionStepValue(Kratos.PRESSURE, 200.0)
        node_3.SetSolutionStepValue(Kratos.PRESSURE, 300.0)

        node_1.SetValue(Kratos.RESIDUAL_NORM, 1000.0)
        node_2.SetValue(Kratos.RESIDUAL_NORM, 2000.0)
        node_3.SetValue(Kratos.RESIDUAL_NORM, 3000.0)

        element_refinement_utilities.InterpolateToRefinedMeshFromCoarseElement(
            element)

        element_refinement_utilities.SetIds(Kratos.SCALE)
        element_refinement_utilities.SetConditionParentIds(Kratos.FIRST_TIME_STEP)

        Kratos.ModelPartIO("element_refinement_utilities/2d_refined_model_part", Kratos.ModelPartIO.WRITE | Kratos.ModelPartIO.IGNORE_VARIABLES_ERROR.AsFalse()).WriteModelPart(refined_model_part)

        # check both files
        with open("element_refinement_utilities/2d_refined_model_part_ref.mdpa", "r") as file_input:
            ref_lines = file_input.readlines()

        with open("element_refinement_utilities/2d_refined_model_part.mdpa", "r") as file_input:
            lines = file_input.readlines()

        self.assertEqual(ref_lines, lines)


    # def test_RefineElement3D(self):
    #     model = Kratos.Model()
    #     model_part = model.CreateModelPart("test")
    #     model_part.ProcessInfo[Kratos.DOMAIN_SIZE] = 3

    #     model_part.AddNodalSolutionStepVariable(Kratos.PRESSURE)
    #     model_part.AddNodalSolutionStepVariable(Kratos.VELOCITY)
    #     model_part.AddNodalSolutionStepVariable(Kratos.NORMAL)

    #     node_1 = model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
    #     node_2 = model_part.CreateNewNode(2, 0.0, 1.0, 0.0)
    #     node_3 = model_part.CreateNewNode(3, 1.0, 1.0, 0.0)
    #     node_4 = model_part.CreateNewNode(4, 1.0 / 3.0, 2.0 / 3.0, 1.0)

    #     properties = model_part.GetProperties()[0]
    #     element = model_part.CreateNewElement("Element3D4N", 1, [1, 2, 3, 4], properties)

    #     p_condition_1 = model_part.CreateNewCondition("SurfaceCondition3D3N", 1, [1, 2, 3], properties)
    #     p_condition_2 = model_part.CreateNewCondition("SurfaceCondition3D3N", 2, [1, 4, 3], properties)
    #     p_condition_3 = model_part.CreateNewCondition("SurfaceCondition3D3N", 3, [1, 2, 4], properties)
    #     p_condition_4 = model_part.CreateNewCondition("SurfaceCondition3D3N", 4, [2, 3, 4], properties)

    #     p_condition_2.SetValue(Kratos.DISTANCE, 10)

    #     p_condition_1.Set(Kratos.INLET)
    #     node_1.Set(Kratos.INLET)
    #     node_2.Set(Kratos.INLET)
    #     p_condition_2.Set(Kratos.OUTLET)
    #     node_2.Set(Kratos.OUTLET)
    #     node_3.Set(Kratos.OUTLET)
    #     p_condition_3.Set(Kratos.SLIP)
    #     node_3.Set(Kratos.SLIP)
    #     node_1.Set(Kratos.SLIP)

    #     tmoc = Kratos.TetrahedralMeshOrientationCheck
    #     throw_errors = False
    #     flags = tmoc.COMPUTE_NODAL_NORMALS | tmoc.COMPUTE_CONDITION_NORMALS | tmoc.ASSIGN_NEIGHBOUR_ELEMENTS_TO_CONDITIONS
    #     Kratos.TetrahedralMeshOrientationCheck(model_part, throw_errors, flags).Execute()

    #     element_refinement_utilities = MeshingApplication.ElementRefinementUtilities(
    #         model_part,
    #         "refined_model_part",
    #         element,
    #         1,
    #         "Element3D4N",
    #         "SurfaceCondition3D3N",
    #         ["VELOCITY_X", "PRESSURE"],
    #         ["RESIDUAL_NORM"],
    #         ["INLET", "OUTLET", "SLIP"])

    #     refined_model_part = element_refinement_utilities.GetRefinedModelPart()

        # node_1.X += 1.0
        # node_2.X += 1.0
        # node_3.X += 1.0
        # node_4.X += 1.0

        # node_1.Y += 1.0
        # node_2.Y += 1.0
        # node_3.Y += 1.0
        # node_4.Y += 1.0

        # node_1.Z += 1.0
        # node_2.Z += 1.0
        # node_3.Z += 1.0
        # node_4.Z += 1.0

        # node_1.SetSolutionStepValue(Kratos.VELOCITY, Kratos.Array3([1.0, 1.0, 1.0]))
        # node_2.SetSolutionStepValue(Kratos.VELOCITY, Kratos.Array3([10.0, 10.0, 10.0]))
        # node_3.SetSolutionStepValue(Kratos.VELOCITY, Kratos.Array3([20.0, 20.0, 20.0]))
        # node_4.SetSolutionStepValue(Kratos.VELOCITY, Kratos.Array3([30.0, 30.0, 30.0]))

        # node_1.SetSolutionStepValue(Kratos.PRESSURE, 100.0)
        # node_2.SetSolutionStepValue(Kratos.PRESSURE, 200.0)
        # node_3.SetSolutionStepValue(Kratos.PRESSURE, 300.0)
        # node_4.SetSolutionStepValue(Kratos.PRESSURE, 400.0)

        # node_1.SetValue(Kratos.RESIDUAL_NORM, 1000.0)
        # node_2.SetValue(Kratos.RESIDUAL_NORM, 2000.0)
        # node_3.SetValue(Kratos.RESIDUAL_NORM, 3000.0)
        # node_4.SetValue(Kratos.RESIDUAL_NORM, 4000.0)

        # element_refinement_utilities.InterpolateToRefinedMeshFromCoarseElement(
        #     element)

        # element_refinement_utilities.SetIds(Kratos.SCALE)
        # element_refinement_utilities.SetConditionParentIds(Kratos.FIRST_TIME_STEP)

        # Kratos.ModelPartIO("3d_refined_model_part", Kratos.ModelPartIO.WRITE).WriteModelPart(refined_model_part)
        # Kratos.ModelPartIO("3d_coarse_model_part", Kratos.ModelPartIO.WRITE).WriteModelPart(model_part)

        # hdf5_file_settings = Kratos.Parameters("""{
        #     "file_name"       : "3d_coarse.h5",
        #     "file_access_mode": "truncate"
        # }""")
        # hdf5_file = KratosHDF5.HDF5FileSerial(hdf5_file_settings)
        # KratosHDF5.HDF5ModelPartIO(hdf5_file, "/ModelData").WriteModelPart(model_part)
        # nodal_io = KratosHDF5.HDF5NodalSolutionStepDataIO(Kratos.Parameters("""{"list_of_variables": ["ALL_VARIABLES_FROM_VARIABLES_LIST"], "prefix":"/ResultsData" }"""), hdf5_file)
        # nodal_io.WriteNodalResults(model_part, 0)
        # nodal_io = KratosHDF5.HDF5NodalDataValueIO(Kratos.Parameters("""{"list_of_variables": ["RESIDUAL_NORM"], "prefix":"/ResultsData" }"""), hdf5_file)
        # nodal_io.WriteNodalResults(model_part.Nodes)

        # hdf5_file_settings = Kratos.Parameters("""{
        #     "file_name"       : "3d_refined.h5",
        #     "file_access_mode": "truncate"
        # }""")

        # hdf5_file = KratosHDF5.HDF5FileSerial(hdf5_file_settings)
        # KratosHDF5.HDF5ModelPartIO(hdf5_file, "/ModelData").WriteModelPart(refined_model_part)
        # nodal_io = KratosHDF5.HDF5NodalSolutionStepDataIO(Kratos.Parameters("""{"list_of_variables": ["ALL_VARIABLES_FROM_VARIABLES_LIST"], "prefix":"/ResultsData" }"""), hdf5_file)
        # nodal_io.WriteNodalResults(refined_model_part, 0)
        # nodal_io = KratosHDF5.HDF5NodalDataValueIO(Kratos.Parameters("""{"list_of_variables": ["RESIDUAL_NORM"], "prefix":"/ResultsData" }"""), hdf5_file)
        # nodal_io.WriteNodalResults(refined_model_part.Nodes)
        # nodal_io = KratosHDF5.HDF5ConditionDataValueIO(Kratos.Parameters("""{"list_of_variables": ["DISTANCE", "SCALE", "FIRST_TIME_STEP"], "prefix":"/ResultsData" }"""), hdf5_file)
        # nodal_io.WriteConditionResults(refined_model_part.Conditions)
        # nodal_io = KratosHDF5.HDF5ElementDataValueIO(Kratos.Parameters("""{"list_of_variables": ["SCALE"], "prefix":"/ResultsData" }"""), hdf5_file)
        # nodal_io.WriteElementResults(refined_model_part.Elements)
        # nodal_io = KratosHDF5.HDF5NodalFlagValueIO(Kratos.Parameters("""{"list_of_variables": ["ACTIVE", "SLIP", "INLET", "OUTLET"], "prefix":"/ResultsData" }"""), hdf5_file)
        # nodal_io.WriteNodalFlags(refined_model_part.Nodes)
        # nodal_io = KratosHDF5.HDF5ConditionFlagValueIO(Kratos.Parameters("""{"list_of_variables": ["ACTIVE", "SLIP", "INLET", "OUTLET"], "prefix":"/ResultsData" }"""), hdf5_file)
        # nodal_io.WriteConditionFlags(refined_model_part.Conditions)


if __name__ == '__main__':
    KratosUnittest.main()
# We import the libraries
from math import sqrt

import KratosMultiphysics as Kratos
import KratosMultiphysics.FluidDynamicsApplication as KratosFluid

import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.kratos_utilities import DeleteFileIfExisting
from KratosMultiphysics.kratos_utilities import DeleteTimeFiles

class TestMPIMetrics(KratosUnittest.TestCase):
    @classmethod
    def tearDownClass(cls):
        DeleteTimeFiles("element_refinement_utilities")

    def _CheckMDPAEqual(self, mdpa_1_name, mdpa_2_name):
        with open(mdpa_1_name, "r") as file_input:
            lines_1 = file_input.readlines()

        with open(mdpa_2_name, "r") as file_input:
            lines_2 = file_input.readlines()

        self.assertEqual(lines_1, lines_2)

    def _CalculateUnitVector(self, v):
        return v * (1.0 / sqrt(v[0]**2 + v[1]**2 + v[2]**2))

    def _CalculateSurfaceArea(self, refined_model_part, surface_index):
        area = 0.0
        for condition in refined_model_part.GetSubModelPart("Surface_{:d}".format(surface_index)).Conditions:
            area += condition.GetGeometry().Area()

        return area

    def _CalculateVolume(self, refined_model_part, volume_computation_method):
        volume = 0.0
        for element in refined_model_part.Elements:
            volume += volume_computation_method(element.GetGeometry())

        return volume

    def _CheckNormals(self, coarse_condition, refined_surface_model_part):
        coarse_unit_normal = self._CalculateUnitVector(coarse_condition.GetValue(Kratos.NORMAL))
        for condition in refined_surface_model_part.Conditions:
            self.assertVectorAlmostEqual(coarse_unit_normal, self._CalculateUnitVector(condition.GetValue(Kratos.NORMAL)), 9)

    def _ExecuteForAllThreadLocalMeshes(self, model, f):
        for model_part_name in model.GetModelPartNames():
            if "." not in model_part_name and model_part_name.startswith("Thread_"):
                refined_model_part = model.GetModelPart(model_part_name)
                f(refined_model_part)

    def test_ElementRefinementProcess2D(self):
        self.addCleanup(lambda: DeleteFileIfExisting("element_refinement_utilities/2d_refined_model_part.mdpa"))

        model = Kratos.Model()
        model_part = model.CreateModelPart("test")
        model_part.ProcessInfo[Kratos.DOMAIN_SIZE] = 2

        model_part.AddNodalSolutionStepVariable(Kratos.PRESSURE)
        model_part.AddNodalSolutionStepVariable(Kratos.VELOCITY)
        model_part.AddNodalSolutionStepVariable(Kratos.NORMAL)

        node_1 = model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
        node_2 = model_part.CreateNewNode(2, 0.0, 1.5, 0.0)
        node_3 = model_part.CreateNewNode(3, 1.0, 1.0, 0.0)

        properties = model_part.GetProperties()[0]
        element = model_part.CreateNewElement("Element2D3N", 1, [1, 2, 3], properties)

        condition_1 = model_part.CreateNewCondition("LineCondition2D2N", 1, [1, 2], properties)
        condition_2 = model_part.CreateNewCondition("LineCondition2D2N", 2, [2, 3], properties)
        condition_3 = model_part.CreateNewCondition("LineCondition2D2N", 3, [3, 1], properties)
        condition_2.SetValue(Kratos.DISTANCE, 10)

        condition_1.Set(Kratos.INLET)
        node_1.Set(Kratos.INLET)
        node_2.Set(Kratos.INLET)
        condition_2.Set(Kratos.OUTLET)
        node_2.Set(Kratos.OUTLET)
        node_3.Set(Kratos.OUTLET)
        condition_3.Set(Kratos.SLIP)
        node_3.Set(Kratos.SLIP)
        node_1.Set(Kratos.SLIP)

        tmoc = Kratos.TetrahedralMeshOrientationCheck
        throw_errors = False
        flags = tmoc.COMPUTE_NODAL_NORMALS | tmoc.COMPUTE_CONDITION_NORMALS | tmoc.ASSIGN_NEIGHBOUR_ELEMENTS_TO_CONDITIONS
        Kratos.TetrahedralMeshOrientationCheck(model_part, throw_errors, flags).Execute()

        parameters = Kratos.Parameters("""{
            "model_part_name"             : "test",
            "refinement_model_part_prefix": "Thread_",
            "refinement_level"            : 4,
            "echo_level"                  : 0,
            "refined_element_name"        : "Element2D3N",
            "refined_condition_name"      : "LineCondition2D2N",
            "nodal_interpolation_settings": {
                "historical_variables_list"    : ["VELOCITY_X", "PRESSURE"],
                "non_hitsorical_variables_list": ["RESIDUAL_NORM"],
                "flags_list"                   : ["INLET", "OUTLET", "SLIP"]
            }
        }""")

        element_refinement_process = KratosFluid.ElementRefinementProcess(model, parameters)
        element_refinement_process.ExecuteInitialize()

        def check_refined_model_part_area_and_volume(model_part):
            # # check surface areas
            self.assertAlmostEqual(self._CalculateSurfaceArea(model_part, 1), condition_1.GetGeometry().Area(), 12)
            self.assertAlmostEqual(self._CalculateSurfaceArea(model_part, 2), condition_3.GetGeometry().Area(), 12)
            self.assertAlmostEqual(self._CalculateSurfaceArea(model_part, 3), condition_2.GetGeometry().Area(), 12)

            # check volume
            self.assertAlmostEqual(self._CalculateVolume(model_part, lambda x: x.Area()), element.GetGeometry().Area(), 12)

        self._ExecuteForAllThreadLocalMeshes(model, check_refined_model_part_area_and_volume)

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

        element_refinement_process.InterpolateAllRefinedMeshesFromCoarseElement(element)

        element_refinement_process.SetEntityIds(Kratos.SCALE)
        element_refinement_process.SetConditionParentIds(Kratos.FIRST_TIME_STEP)

        def check_refined_model_part_interpolations(model_part):
            Kratos.NormalCalculationUtils().CalculateUnitNormals(model_part)
            self._CheckNormals(condition_1, model_part.GetSubModelPart("Surface_1"))
            self._CheckNormals(condition_3, model_part.GetSubModelPart("Surface_2"))
            self._CheckNormals(condition_2, model_part.GetSubModelPart("Surface_3"))
            Kratos.ModelPartIO("element_refinement_utilities/2d_refined_model_part", Kratos.ModelPartIO.WRITE | Kratos.ModelPartIO.IGNORE_VARIABLES_ERROR.AsFalse()).WriteModelPart(model_part)
            self._CheckMDPAEqual("element_refinement_utilities/2d_refined_model_part_ref.mdpa", "element_refinement_utilities/2d_refined_model_part.mdpa")

        self._ExecuteForAllThreadLocalMeshes(model, check_refined_model_part_interpolations)

    def test_ElementRefinementProcess3D(self):
        self.addCleanup(lambda: DeleteFileIfExisting("element_refinement_utilities/3d_refined_model_part.mdpa"))

        model = Kratos.Model()
        model_part = model.CreateModelPart("test")
        model_part.ProcessInfo[Kratos.DOMAIN_SIZE] = 3

        model_part.AddNodalSolutionStepVariable(Kratos.PRESSURE)
        model_part.AddNodalSolutionStepVariable(Kratos.VELOCITY)
        model_part.AddNodalSolutionStepVariable(Kratos.NORMAL)

        node_1 = model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
        node_2 = model_part.CreateNewNode(2, 0.0, 1.2, 0.0)
        node_3 = model_part.CreateNewNode(3, 1.1, 1.3, 0.0)
        node_4 = model_part.CreateNewNode(4, 1.0 / 3.0, 2.0 / 3.0, 1.0)

        properties = model_part.GetProperties()[0]
        element = model_part.CreateNewElement("Element3D4N", 1, [1, 2, 3, 4], properties)

        condition_1 = model_part.CreateNewCondition("SurfaceCondition3D3N", 1, [1, 2, 3], properties)
        condition_2 = model_part.CreateNewCondition("SurfaceCondition3D3N", 2, [1, 4, 3], properties)
        condition_3 = model_part.CreateNewCondition("SurfaceCondition3D3N", 3, [1, 2, 4], properties)
        condition_4 = model_part.CreateNewCondition("SurfaceCondition3D3N", 4, [2, 3, 4], properties)

        condition_2.SetValue(Kratos.DISTANCE, 10)

        condition_1.Set(Kratos.INLET)
        node_1.Set(Kratos.INLET)
        node_2.Set(Kratos.INLET)
        node_3.Set(Kratos.INLET)
        condition_2.Set(Kratos.OUTLET)
        node_1.Set(Kratos.OUTLET)
        node_4.Set(Kratos.OUTLET)
        node_3.Set(Kratos.OUTLET)
        condition_3.Set(Kratos.SLIP)
        node_1.Set(Kratos.SLIP)
        node_2.Set(Kratos.SLIP)
        node_4.Set(Kratos.SLIP)

        tmoc = Kratos.TetrahedralMeshOrientationCheck
        throw_errors = False
        flags = tmoc.COMPUTE_NODAL_NORMALS | tmoc.COMPUTE_CONDITION_NORMALS | tmoc.ASSIGN_NEIGHBOUR_ELEMENTS_TO_CONDITIONS
        Kratos.TetrahedralMeshOrientationCheck(model_part, throw_errors, flags).Execute()

        parameters = Kratos.Parameters("""{
            "model_part_name"             : "test",
            "refinement_model_part_prefix": "Thread_",
            "refinement_level"            : 1,
            "echo_level"                  : 0,
            "refined_element_name"        : "Element3D4N",
            "refined_condition_name"      : "SurfaceCondition3D3N",
            "nodal_interpolation_settings": {
                "historical_variables_list"    : ["VELOCITY_X", "PRESSURE"],
                "non_hitsorical_variables_list": ["RESIDUAL_NORM"],
                "flags_list"                   : ["INLET", "OUTLET", "SLIP"]
            }
        }""")

        element_refinement_process = KratosFluid.ElementRefinementProcess(model, parameters)
        element_refinement_process.ExecuteInitialize()

        def check_refined_model_part_area_and_volume(model_part):
            # # check surface areas
            self.assertAlmostEqual(self._CalculateSurfaceArea(model_part, 1), condition_3.GetGeometry().Area(), 12)
            self.assertAlmostEqual(self._CalculateSurfaceArea(model_part, 2), condition_2.GetGeometry().Area(), 12)
            self.assertAlmostEqual(self._CalculateSurfaceArea(model_part, 3), condition_4.GetGeometry().Area(), 12)
            self.assertAlmostEqual(self._CalculateSurfaceArea(model_part, 4), condition_1.GetGeometry().Area(), 12)

            # check volume
            self.assertAlmostEqual(self._CalculateVolume(model_part, lambda x: x.Volume()), element.GetGeometry().Volume(), 12)

        self._ExecuteForAllThreadLocalMeshes(model, check_refined_model_part_area_and_volume)

        node_1.X += 1.0
        node_2.X += 1.0
        node_3.X += 1.0
        node_4.X += 1.0

        node_1.Y += 1.0
        node_2.Y += 1.0
        node_3.Y += 1.0
        node_4.Y += 1.0

        node_1.Z += 1.0
        node_2.Z += 1.0
        node_3.Z += 1.0
        node_4.Z += 1.0

        node_1.SetSolutionStepValue(Kratos.VELOCITY, Kratos.Array3([1.0, 1.0, 1.0]))
        node_2.SetSolutionStepValue(Kratos.VELOCITY, Kratos.Array3([10.0, 10.0, 10.0]))
        node_3.SetSolutionStepValue(Kratos.VELOCITY, Kratos.Array3([20.0, 20.0, 20.0]))
        node_4.SetSolutionStepValue(Kratos.VELOCITY, Kratos.Array3([30.0, 30.0, 30.0]))

        node_1.SetSolutionStepValue(Kratos.PRESSURE, 100.0)
        node_2.SetSolutionStepValue(Kratos.PRESSURE, 200.0)
        node_3.SetSolutionStepValue(Kratos.PRESSURE, 300.0)
        node_4.SetSolutionStepValue(Kratos.PRESSURE, 400.0)

        node_1.SetValue(Kratos.RESIDUAL_NORM, 1000.0)
        node_2.SetValue(Kratos.RESIDUAL_NORM, 2000.0)
        node_3.SetValue(Kratos.RESIDUAL_NORM, 3000.0)
        node_4.SetValue(Kratos.RESIDUAL_NORM, 4000.0)

        element_refinement_process.InterpolateAllRefinedMeshesFromCoarseElement(element)

        element_refinement_process.SetEntityIds(Kratos.SCALE)
        element_refinement_process.SetConditionParentIds(Kratos.FIRST_TIME_STEP)

        def check_refined_model_part_interpolations(model_part):
            Kratos.NormalCalculationUtils().CalculateUnitNormals(model_part)
            self._CheckNormals(condition_3, model_part.GetSubModelPart("Surface_1"))
            self._CheckNormals(condition_2, model_part.GetSubModelPart("Surface_2"))
            self._CheckNormals(condition_4, model_part.GetSubModelPart("Surface_3"))
            self._CheckNormals(condition_1, model_part.GetSubModelPart("Surface_4"))
            Kratos.ModelPartIO("element_refinement_utilities/3d_refined_model_part", Kratos.ModelPartIO.WRITE | Kratos.ModelPartIO.IGNORE_VARIABLES_ERROR.AsFalse()).WriteModelPart(model_part)
            self._CheckMDPAEqual("element_refinement_utilities/3d_refined_model_part_ref.mdpa", "element_refinement_utilities/3d_refined_model_part.mdpa")

        self._ExecuteForAllThreadLocalMeshes(model, check_refined_model_part_interpolations)

if __name__ == '__main__':
    KratosUnittest.main()
# We import the libraries
import KratosMultiphysics as Kratos
import KratosMultiphysics.MeshingApplication as MeshingApplication

import KratosMultiphysics.KratosUnittest as KratosUnittest

from math import sqrt
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

    def test_RefineElement2D(self):
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

        element_refinement_utilities = MeshingApplication.ElementRefinementUtilities(
            model_part,
            "refined_model_part",
            element,
            4,
            "Element2D3N",
            "LineCondition2D2N",
            ["VELOCITY_X", "PRESSURE"],
            ["RESIDUAL_NORM"],
            ["INLET", "OUTLET", "SLIP"])

        refined_model_part = element_refinement_utilities.GetRefinedModelPart()

        # # check surface areas
        self.assertAlmostEqual(self._CalculateSurfaceArea(refined_model_part, 1), condition_1.GetGeometry().Area(), 12)
        self.assertAlmostEqual(self._CalculateSurfaceArea(refined_model_part, 2), condition_3.GetGeometry().Area(), 12)
        self.assertAlmostEqual(self._CalculateSurfaceArea(refined_model_part, 3), condition_2.GetGeometry().Area(), 12)

        # check volume
        self.assertAlmostEqual(self._CalculateVolume(refined_model_part, lambda x: x.Area()), element.GetGeometry().Area(), 12)

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

        element_refinement_utilities.ComputeSurfaceMap()
        element_refinement_utilities.InterpolateToRefinedMeshFromCoarseElement(element)

        Kratos.NormalCalculationUtils().CalculateUnitNormals(refined_model_part)
        self._CheckNormals(condition_1, refined_model_part.GetSubModelPart("Surface_1"))
        self._CheckNormals(condition_3, refined_model_part.GetSubModelPart("Surface_2"))
        self._CheckNormals(condition_2, refined_model_part.GetSubModelPart("Surface_3"))

        element_refinement_utilities.SetIds(Kratos.SCALE)
        element_refinement_utilities.SetConditionParentIds(Kratos.FIRST_TIME_STEP)

        Kratos.ModelPartIO("element_refinement_utilities/2d_refined_model_part", Kratos.ModelPartIO.WRITE | Kratos.ModelPartIO.IGNORE_VARIABLES_ERROR.AsFalse()).WriteModelPart(refined_model_part)
        self._CheckMDPAEqual("element_refinement_utilities/2d_refined_model_part_ref.mdpa", "element_refinement_utilities/2d_refined_model_part.mdpa")

    def test_RefineElement3D(self):
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

        element_refinement_utilities = MeshingApplication.ElementRefinementUtilities(
            model_part,
            "refined_model_part",
            element,
            1,
            "Element3D4N",
            "SurfaceCondition3D3N",
            ["VELOCITY_X", "PRESSURE"],
            ["RESIDUAL_NORM"],
            ["INLET", "OUTLET", "SLIP"])

        refined_model_part = element_refinement_utilities.GetRefinedModelPart()

        # check surface areas
        self.assertAlmostEqual(self._CalculateSurfaceArea(refined_model_part, 1), condition_3.GetGeometry().Area(), 12)
        self.assertAlmostEqual(self._CalculateSurfaceArea(refined_model_part, 2), condition_2.GetGeometry().Area(), 12)
        self.assertAlmostEqual(self._CalculateSurfaceArea(refined_model_part, 3), condition_4.GetGeometry().Area(), 12)
        self.assertAlmostEqual(self._CalculateSurfaceArea(refined_model_part, 4), condition_1.GetGeometry().Area(), 12)

        # check volume
        self.assertAlmostEqual(self._CalculateVolume(refined_model_part, lambda x: x.Volume()), element.GetGeometry().Volume(), 12)

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

        element_refinement_utilities.ComputeSurfaceMap()
        element_refinement_utilities.InterpolateToRefinedMeshFromCoarseElement(element)
        Kratos.NormalCalculationUtils().CalculateUnitNormals(refined_model_part)
        self._CheckNormals(condition_3, refined_model_part.GetSubModelPart("Surface_1"))
        self._CheckNormals(condition_2, refined_model_part.GetSubModelPart("Surface_2"))
        self._CheckNormals(condition_4, refined_model_part.GetSubModelPart("Surface_3"))
        self._CheckNormals(condition_1, refined_model_part.GetSubModelPart("Surface_4"))

        element_refinement_utilities.SetIds(Kratos.SCALE)
        element_refinement_utilities.SetConditionParentIds(Kratos.FIRST_TIME_STEP)

        Kratos.ModelPartIO("element_refinement_utilities/3d_refined_model_part", Kratos.ModelPartIO.WRITE | Kratos.ModelPartIO.IGNORE_VARIABLES_ERROR.AsFalse()).WriteModelPart(refined_model_part)
        self._CheckMDPAEqual("element_refinement_utilities/3d_refined_model_part_ref.mdpa", "element_refinement_utilities/3d_refined_model_part.mdpa")


if __name__ == '__main__':
    KratosUnittest.main()
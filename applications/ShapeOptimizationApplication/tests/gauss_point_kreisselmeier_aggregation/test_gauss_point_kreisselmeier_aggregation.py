import KratosMultiphysics as Kratos
import KratosMultiphysics.ShapeOptimizationApplication as KratosSOA
import KratosMultiphysics.StructuralMechanicsApplication as KratosStructural
import KratosMultiphysics.KratosUnittest as UnitTest
from KratosMultiphysics.ShapeOptimizationApplication.response_functions.total_volume import TotalVolume


class TestGaussPointKreisselmeierAggregationResponseFunction(UnitTest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.primal_model = Kratos.Model()
        cls.primal_model_part = cls.primal_model.CreateModelPart("structure")
        cls.primal_model_part.AddNodalSolutionStepVariable(Kratos.DISPLACEMENT)
        Kratos.ModelPartIO("structure").ReadModelPart(cls.primal_model_part)
        material_settings = Kratos.Parameters("""{"Parameters": {"materials_filename": ""}} """)
        material_settings["Parameters"]["materials_filename"].SetString("material.json")
        Kratos.ReadMaterialsUtility(material_settings, cls.primal_model)
        Kratos.EntitiesUtilities.InitializeAllEntities(cls.primal_model_part)
        Kratos.EntitiesUtilities.InitializeSolutionStepAllEntities(cls.primal_model_part)
        Kratos.EntitiesUtilities.InitializeNonLinearIterationAllEntities(cls.primal_model_part)

        parameters = Kratos.Parameters("""{
            "model_part_name": "test",
            "gauss_point_values_scalar_variable": "VON_MISES_STRESS",
            "gauss_point_value_scaling_factor"  : 1.14859e+15,
            "rho": 50.0,
            "gradient_mode": "semi_analytic",
            "gradient_mode_settings": {
                "perturbation_size": 1e-7,
                "perturbation_variable_name": "PERTURBATION_SIZE",
                "list_of_scalar_primal_state_variables": [
                    "DISPLACEMENT_X",
                    "DISPLACEMENT_Y",
                    "DISPLACEMENT_Z"
                ]
            }
        }""")
        cls.adjoint_model = Kratos.Model()
        cls.adjoint_model_part = cls.adjoint_model.CreateModelPart("structure")
        cls.adjoint_model_part.AddNodalSolutionStepVariable(Kratos.DISPLACEMENT)
        cls.adjoint_model_part.AddNodalSolutionStepVariable(KratosStructural.ADJOINT_DISPLACEMENT)
        Kratos.ModelPartIO("structure").ReadModelPart(cls.adjoint_model_part)
        replacement_settings = Kratos.Parameters("""
            {
                "element_name_table" :
                {
                    "ShellThinElement3D3N"           : "AdjointFiniteDifferencingShellThinElement3D3N",
                    "CrLinearBeamElement3D2N"        : "AdjointFiniteDifferenceCrBeamElementLinear3D2N",
                    "TrussLinearElement3D2N"         : "AdjointFiniteDifferenceTrussLinearElement3D2N",
                    "TrussElement3D2N"               : "AdjointFiniteDifferenceTrussElement3D2N",
                    "TotalLagrangianElement2D3N"     : "TotalLagrangianAdjointElement2D3N",
                    "TotalLagrangianElement2D4N"     : "TotalLagrangianAdjointElement2D4N",
                    "TotalLagrangianElement2D6N"     : "TotalLagrangianAdjointElement2D6N",
                    "TotalLagrangianElement3D4N"     : "TotalLagrangianAdjointElement3D4N",
                    "TotalLagrangianElement3D8N"     : "TotalLagrangianAdjointElement3D8N",
                    "SmallDisplacementElement3D4N"   : "AdjointFiniteDifferencingSmallDisplacementElement3D4N",
                    "SmallDisplacementElement3D6N"   : "AdjointFiniteDifferencingSmallDisplacementElement3D6N",
                    "SmallDisplacementElement3D8N"   : "AdjointFiniteDifferencingSmallDisplacementElement3D8N",
                    "SpringDamperElement3D2N"        : "AdjointFiniteDifferenceSpringDamperElement3D2N"
                },
                "condition_name_table" :
                {
                    "PointLoadCondition2D1N"         : "AdjointSemiAnalyticPointLoadCondition2D1N",
                    "PointLoadCondition3D1N"         : "AdjointSemiAnalyticPointLoadCondition3D1N",
                    "SurfaceLoadCondition3D3N"       : "AdjointSemiAnalyticSurfaceLoadCondition3D3N",
                    "SurfaceLoadCondition3D4N"       : "AdjointSemiAnalyticSurfaceLoadCondition3D4N",
                    "SmallDisplacementSurfaceLoadCondition3D3N" : "AdjointSemiAnalyticSmallDisplacementSurfaceLoadCondition3D3N",
                    "SmallDisplacementSurfaceLoadCondition3D4N" : "AdjointSemiAnalyticSmallDisplacementSurfaceLoadCondition3D4N",
                    "LineLoadCondition3D2N"                     : "AdjointSemiAnalyticLineLoadCondition3D2N",
                    "SmallDisplacementLineLoadCondition3D2N"    : "AdjointSemiAnalyticSmallDisplacementLineLoadCondition3D2N"
                },
                "ignore_conditions" : [
                    "SurfaceCondition3D3N",
                    "SurfaceCondition3D4N",
                    "PointCondition3D1N"
                ]
            }
        """) # TODO remove "Condition3D" after issue#4439 is resolved

        # KratosStructural.ReplaceMultipleElementsAndConditionsProcess(cls.adjoint_model_part, replacement_settings).Execute()
        material_settings = Kratos.Parameters("""{"Parameters": {"materials_filename": ""}} """)
        material_settings["Parameters"]["materials_filename"].SetString("material.json")
        Kratos.ReadMaterialsUtility(material_settings, cls.adjoint_model)
        Kratos.EntitiesUtilities.InitializeAllEntities(cls.adjoint_model_part)
        Kratos.EntitiesUtilities.InitializeSolutionStepAllEntities(cls.adjoint_model_part)
        Kratos.EntitiesUtilities.InitializeNonLinearIterationAllEntities(cls.adjoint_model_part)

        cls.element_id = 10
        test_model_part = cls.primal_model_part.CreateSubModelPart("test")
        test_model_part.AddElement(cls.primal_model_part.GetElement(cls.element_id), 0)
        test_model_part = cls.adjoint_model_part.CreateSubModelPart("test")
        test_model_part.AddElement(cls.adjoint_model_part.GetElement(cls.element_id), 0)
        cls.response_function = KratosSOA.GaussPointKreisselmeierAggregationResponseFunction(parameters, cls.adjoint_model_part)
        cls.response_function.Initialize()

        index = 0
        for primal_node, adjoint_node in zip(cls.primal_model_part.Nodes, cls.adjoint_model_part.Nodes):
            v = Kratos.Array3([index, index * 3.5, index * index])
            primal_node.SetSolutionStepValue(Kratos.DISPLACEMENT, v)
            adjoint_node.SetSolutionStepValue(Kratos.DISPLACEMENT, v)
            index += 1

        cls.primal_model_part.ProcessInfo[Kratos.DOMAIN_SIZE] = 3
        cls.adjoint_model_part.ProcessInfo[Kratos.DOMAIN_SIZE] = 3
        cls.ref_value = cls.response_function.CalculateValue(cls.primal_model_part)

    def CalculateFiniteDifferenceSensitivity(self, node, variable, delta):
        if variable.Name() == "SHAPE_SENSITIVITY_X":
            def nodal_perturbation_method(node, delta):
                node.X += delta
                node.X0 += delta
        elif variable.Name() == "SHAPE_SENSITIVITY_Y":
            def nodal_perturbation_method(node, delta):
                node.Y += delta
                node.Y0 += delta
        elif variable.Name() == "SHAPE_SENSITIVITY_Z":
            def nodal_perturbation_method(node, delta):
                node.Z += delta
                node.Z0 += delta
        else:
            def nodal_perturbation_method(node, delta):
                node.SetSolutionStepValue(variable, node.GetSolutionStepValue(variable) + delta)

        nodal_perturbation_method(node, delta)
        fd_shape_sensitivity = (self.response_function.CalculateValue(self.primal_model_part) - self.ref_value) / delta
        nodal_perturbation_method(node, -delta)
        return fd_shape_sensitivity

    def testStateVariableSensitivity(self):
        variables_list = [Kratos.DISPLACEMENT_X, Kratos.DISPLACEMENT_Y, Kratos.DISPLACEMENT_Z]
        delta = 2e-5

        element = self.adjoint_model_part.GetElement(self.element_id)
        residual_gradient = Kratos.Matrix()
        response_gradient = Kratos.Vector()
        self.response_function.CalculateGradient(element, residual_gradient, response_gradient, self.adjoint_model_part.ProcessInfo)

        for i_node, node in enumerate(element.GetGeometry()):
            for i_var, variable in enumerate(variables_list):
                fd_sensitivity = self.CalculateFiniteDifferenceSensitivity(self.primal_model_part.GetNode(node.Id), variable, delta)
                self.assertAlmostEqual(fd_sensitivity, response_gradient[i_node * len(variables_list) + i_var], 8)


    def testShapeVariableSensitivity(self):
        variables_list = [Kratos.SHAPE_SENSITIVITY_X, Kratos.SHAPE_SENSITIVITY_Y, Kratos.SHAPE_SENSITIVITY_Z]
        delta = 1e-7

        element = self.adjoint_model_part.GetElement(self.element_id)
        residual_gradient = Kratos.Matrix()
        response_gradient = Kratos.Vector()
        self.response_function.CalculatePartialSensitivity(element, Kratos.SHAPE_SENSITIVITY, residual_gradient, response_gradient, self.adjoint_model_part.ProcessInfo)

        for i_node, node in enumerate(element.GetGeometry()):
            for i_var, variable in enumerate(variables_list):
                fd_sensitivity = self.CalculateFiniteDifferenceSensitivity(self.primal_model_part.GetNode(node.Id), variable, delta)
                self.assertAlmostEqual(fd_sensitivity, response_gradient[i_node * len(variables_list) + i_var], 8)

if __name__ == '__main__':
    UnitTest.main()
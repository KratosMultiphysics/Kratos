
import abc

import KratosMultiphysics as Kratos
from KratosMultiphysics.testing.utilities import ReadModelPart
import KratosMultiphysics.OptimizationApplication as KratosOA
import KratosMultiphysics.StructuralMechanicsApplication

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as kratos_unittest
from KratosMultiphysics.kratos_utilities import DeleteFileIfExisting
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.controls.shape.vertex_morphing_shape_control import VertexMorphingShapeControl

class TestVMShapeControlBase:
    @classmethod
    @abc.abstractmethod
    def GetImplicitControlParameters(self) -> Kratos.Parameters:
        pass

    @classmethod
    @abc.abstractmethod
    def GetExplicitControlParameters(self) -> Kratos.Parameters:
        pass

    @classmethod
    @abc.abstractmethod
    def GetMdpaFileName(self) -> str:
        pass

    @classmethod
    def setUpClass(cls) -> None:
        cls.model = Kratos.Model()
        cls.model_part = cls.model.CreateModelPart("test")
        cls.model_part.ProcessInfo[Kratos.DOMAIN_SIZE] = 3

        with kratos_unittest.WorkFolderScope(".", __file__):
            cls.optimization_problem = OptimizationProblem()

            cls.implicit_shape_control = VertexMorphingShapeControl("implicit_filter", cls.model, cls.GetImplicitControlParameters(), cls.optimization_problem)
            cls.optimization_problem.AddComponent(cls.implicit_shape_control)

            cls.explicit_shape_control = VertexMorphingShapeControl("explicit_filter", cls.model, cls.GetExplicitControlParameters(), cls.optimization_problem)
            cls.optimization_problem.AddComponent(cls.explicit_shape_control)

            ReadModelPart(cls.GetMdpaFileName(), cls.model_part)
            cls.implicit_shape_control.Initialize()
            cls.explicit_shape_control.Initialize()

        cls.initial_nodal_positions_exp = Kratos.Expression.NodalExpression(cls.model_part)
        Kratos.Expression.NodalPositionExpressionIO.Read(cls.initial_nodal_positions_exp, Kratos.Configuration.Initial)

    def setUp(self) -> None:
        Kratos.Expression.NodalPositionExpressionIO.Write(self.initial_nodal_positions_exp, Kratos.Configuration.Initial)
        Kratos.Expression.NodalPositionExpressionIO.Write(self.initial_nodal_positions_exp, Kratos.Configuration.Current)

    @classmethod
    def tearDownClass(cls):
        cls.implicit_shape_control.Finalize()
        cls.explicit_shape_control.Finalize()
        del cls.explicit_shape_control
        del cls.implicit_shape_control
        with kratos_unittest.WorkFolderScope(".", __file__):
            DeleteFileIfExisting(f"{cls.GetMdpaFileName()}.time")

@kratos_unittest.skipIfApplicationsNotAvailable("StructuralMechanicsApplication", "MeshMovingApplication")
class TestVMShapeControlShell(TestVMShapeControlBase, kratos_unittest.TestCase):
    @classmethod
    def GetImplicitControlParameters(self) -> Kratos.Parameters:
        return Kratos.Parameters("""{
                "controlled_model_part_names": ["test"],
                "filter_settings": {
                    "filter_type"                  : "implicit_filter",
                    "filter_radius"                : 0.2,
                    "filtering_boundary_conditions": {
                        "test.top_edge"    : [true,true,true],
                        "test.edge_support": [true,true,true]
                    }
                }
            }""")

    @classmethod
    def GetExplicitControlParameters(self) -> Kratos.Parameters:
        return Kratos.Parameters("""{
                "controlled_model_part_names": ["test"],
                "filter_settings": {
                    "filter_type"               : "explicit_filter",
                    "filter_function_type"      : "linear",
                    "filter_radius_settings":{
                        "filter_radius_type": "constant",
                        "filter_radius"     : 0.2
                    },
                    "filtering_boundary_conditions": {
                        "damping_type"              : "nearest_entity",
                        "damping_function_type"     : "sigmoidal",
                        "damped_model_part_settings": {
                            "test.top_edge"    : [true,true,true],
                            "test.edge_support": [true,true,true]
                        }
                    }
                }
            }""")

    @classmethod
    def GetMdpaFileName(self) -> str:
        return "shell"

    def test_GetControlField(self):
        self.assertAlmostEqual(Kratos.Expression.Utils.NormL2(self.implicit_shape_control.GetControlField()), 22.4268730405143, 10)
        self.assertAlmostEqual(Kratos.Expression.Utils.NormL2(self.explicit_shape_control.GetControlField()), 0.0, 10)

    def test_GetPhysicalField(self):
        shape_field = self.implicit_shape_control.GetPhysicalField()
        self.assertAlmostEqual(Kratos.Expression.Utils.NormL2(shape_field), 20.976176963410882, 10)
        shape_field = self.explicit_shape_control.GetPhysicalField()
        self.assertAlmostEqual(Kratos.Expression.Utils.NormL2(shape_field), 20.976176963410882, 10)

    def test_MapGradient(self):
        physical_gradient =  self.implicit_shape_control.GetEmptyField()
        ta_physical_gradient = Kratos.TensorAdaptors.DoubleTensorAdaptor(physical_gradient.GetContainer(), Kratos.DoubleNDData(physical_gradient.Evaluate()), copy=False)
        ta_physical_gradient.data[:] = 1.0
        self.explicit_shape_control.filter.filter_utils.GetIntegrationWeights(ta_physical_gradient)
        Kratos.Expression.CArrayExpressionIO.Read(physical_gradient, ta_physical_gradient.data)
        implicit_mapped_gradient = self.implicit_shape_control.MapGradient({KratosOA.SHAPE: physical_gradient})
        self.assertAlmostEqual(Kratos.Expression.Utils.NormL2(implicit_mapped_gradient), 27.84296340221239, 10)
        explicit_mapped_gradient = self.explicit_shape_control.MapGradient({KratosOA.SHAPE: physical_gradient})
        self.assertAlmostEqual(Kratos.Expression.Utils.NormL2(explicit_mapped_gradient), 28.075348639115425, 10)

    def test_UpdateImplicit(self):
        update_field = self.implicit_shape_control.GetEmptyField()
        constant_field_value = Kratos.Array3([0.1, 0.1, 0.1])
        Kratos.Expression.LiteralExpressionIO.SetData(update_field, constant_field_value)
        self.implicit_shape_control.Update(update_field)
        implicit_control_field = self.implicit_shape_control.GetControlField()
        implicit_shape_field = self.implicit_shape_control.GetPhysicalField()
        self.assertAlmostEqual(Kratos.Expression.Utils.NormL2(implicit_control_field), 3.633180424916991, 10)
        self.assertAlmostEqual(Kratos.Expression.Utils.NormL2(implicit_shape_field), 10.365298105786017, 10)

    def test_UpdateExplicit(self):
        update_field = self.explicit_shape_control.GetEmptyField()
        constant_field_value = Kratos.Array3([0.1, 0.1, 0.1])
        Kratos.Expression.LiteralExpressionIO.SetData(update_field, constant_field_value)
        self.explicit_shape_control.Update(update_field)
        explicit_control_field = self.explicit_shape_control.GetControlField()
        explicit_shape_field = self.explicit_shape_control.GetPhysicalField()
        self.assertAlmostEqual(Kratos.Expression.Utils.NormL2(explicit_control_field), 3.633180424916991, 10)
        self.assertAlmostEqual(Kratos.Expression.Utils.NormL2(explicit_shape_field), 22.013379276492945, 10)

@kratos_unittest.skipIfApplicationsNotAvailable("StructuralMechanicsApplication", "MeshMovingApplication")
class TestVMShapeControlSolid(TestVMShapeControlBase, kratos_unittest.TestCase):
    @classmethod
    def GetImplicitControlParameters(self) -> Kratos.Parameters:
        return Kratos.Parameters("""{
                "controlled_model_part_names": ["test.design"],
                "mesh_motion_solver_type"    : "filter_based",
                "filter_settings": {
                    "filter_type"                  : "implicit_filter",
                    "filter_radius"                : 0.5,
                    //"linear_solver_settings": {"solver_type": "amgcl", "verbosity": 4},
                    "filtering_boundary_conditions": {
                        "test.fixed" : [true,true,true]
                    }
                }
            }""")

    @classmethod
    def GetExplicitControlParameters(self) -> Kratos.Parameters:
        return Kratos.Parameters("""{
                "controlled_model_part_names": ["test.design"],
                "mesh_motion_solver_type"    : "mesh_moving_analysis",
                "filter_settings": {
                    "filter_type"               : "explicit_filter",
                    "filter_function_type"      : "linear",
                    "filter_radius_settings":{
                        "filter_radius_type": "constant",
                        "filter_radius"     : 0.5
                    },
                    "filtering_boundary_conditions": {
                        "damping_type"              : "nearest_entity",
                        "damping_function_type"     : "sigmoidal",
                        "damped_model_part_settings": {
                            "test.fixed" : [true,true,true]
                        }
                    }
                },
                "mesh_motion_solver_settings": {
                    "problem_data": {
                        "echo_level"   : 0,
                        "parallel_type": "OpenMP",
                        "start_time"   : 0.0,
                        "end_time"     : 1.0
                    },
                    "solver_settings" : {
                        "domain_size"            : 3,
                        "echo_level"             : 0,
                        "solver_type"            : "structural_similarity",
                        "model_part_name"        : "test",
                        "compute_reactions"      : false,
                        "calculate_mesh_velocity": false,
                        "model_import_settings"              : {
                            "input_type"     : "use_input_model_part"
                        },
                        "time_stepping" : {
                            "time_step"       : 1.1
                        }
                    }
                }
            }""")

    @classmethod
    def GetMdpaFileName(self) -> str:
        return "solid"

    def test_GetControlField(self):
        self.assertAlmostEqual(Kratos.Expression.Utils.NormL2(self.implicit_shape_control.GetControlField()), 213.27526425126743, 4)
        self.assertAlmostEqual(Kratos.Expression.Utils.NormL2(self.explicit_shape_control.GetControlField()), 0.0000, 4)

    def test_GetPhysicalField(self):
        shape_field = self.explicit_shape_control.GetPhysicalField()
        self.assertAlmostEqual(Kratos.Expression.Utils.NormL2(shape_field), 4.24264068711927, 10)
        shape_field = self.implicit_shape_control.GetPhysicalField()
        self.assertAlmostEqual(Kratos.Expression.Utils.NormL2(shape_field), 4.24264068711927, 10)

    def test_MapGradient(self):
        constant_field_value = Kratos.Array3([1.0, 1.0, 1.0])

        solid_physical_gradient = self.explicit_shape_control.GetEmptyField()
        Kratos.Expression.LiteralExpressionIO.SetData(solid_physical_gradient, constant_field_value)
        mapped_gradient = self.explicit_shape_control.MapGradient({KratosOA.SHAPE: solid_physical_gradient})
        self.assertAlmostEqual(Kratos.Expression.Utils.NormL2(mapped_gradient), 13.52774925846869, 10)

        solid_physical_gradient = self.implicit_shape_control.GetEmptyField()
        Kratos.Expression.LiteralExpressionIO.SetData(solid_physical_gradient, constant_field_value)
        mapped_gradient = self.implicit_shape_control.MapGradient({KratosOA.SHAPE: solid_physical_gradient})
        self.assertAlmostEqual(Kratos.Expression.Utils.NormL2(mapped_gradient), 20.108624626448545, 10)

    def test_UpdateExplicit(self):
        constant_field_value = Kratos.Array3([0.1, 0.1, 0.1])

        update_field = self.explicit_shape_control.GetEmptyField()
        Kratos.Expression.LiteralExpressionIO.SetData(update_field, constant_field_value)
        self.explicit_shape_control.Update(update_field)
        control_field = self.explicit_shape_control.GetControlField()
        shape_field = self.explicit_shape_control.GetPhysicalField()
        self.assertAlmostEqual(Kratos.Expression.Utils.NormL2(control_field), 0.6480740698407862, 10)
        self.assertAlmostEqual(Kratos.Expression.Utils.NormL2(shape_field), 4.524378410345448, 10)

    def test_UpdateImplicit(self):
        constant_field_value = Kratos.Array3([0.1, 0.1, 0.1])
        update_field = self.implicit_shape_control.GetEmptyField()
        Kratos.Expression.LiteralExpressionIO.SetData(update_field, constant_field_value)
        self.implicit_shape_control.Update(update_field)
        control_field = self.implicit_shape_control.GetControlField()
        shape_field = self.implicit_shape_control.GetPhysicalField()
        self.assertAlmostEqual(Kratos.Expression.Utils.NormL2(control_field), 0.6480740698407862, 10)
        self.assertAlmostEqual(Kratos.Expression.Utils.NormL2(shape_field), 4.370590586296503, 10)

if __name__ == "__main__":
    kratos_unittest.main()

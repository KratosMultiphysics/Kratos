
import abc, numpy

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

        cls.initial_nodal_positions_ta = Kratos.TensorAdaptors.NodePositionTensorAdaptor(cls.model_part.Nodes, Kratos.Configuration.Initial)
        cls.initial_nodal_positions_ta.CollectData()

    def setUp(self) -> None:
        Kratos.TensorAdaptors.NodePositionTensorAdaptor(self.initial_nodal_positions_ta, Kratos.Configuration.Initial, copy=False).StoreData()
        Kratos.TensorAdaptors.NodePositionTensorAdaptor(self.initial_nodal_positions_ta, Kratos.Configuration.Current, copy=False).StoreData()

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
        self.assertAlmostEqual(numpy.linalg.norm(self.implicit_shape_control.GetControlField().data), 22.4268730405143, 10)
        self.assertAlmostEqual(numpy.linalg.norm(self.explicit_shape_control.GetControlField().data), 0.0, 10)

    def test_GetPhysicalField(self):
        shape_field = self.implicit_shape_control.GetPhysicalField()
        self.assertAlmostEqual(numpy.linalg.norm(shape_field.data), 20.976176963410882, 10)
        shape_field = self.explicit_shape_control.GetPhysicalField()
        self.assertAlmostEqual(numpy.linalg.norm(shape_field.data), 20.976176963410882, 10)

    def test_MapGradient(self):
        physical_gradient =  self.implicit_shape_control.GetEmptyField()
        physical_gradient.data[:] = 1.0
        self.explicit_shape_control.filter.filter_utils.GetIntegrationWeights(physical_gradient)
        implicit_mapped_gradient = self.implicit_shape_control.MapGradient({KratosOA.SHAPE: physical_gradient})
        self.assertAlmostEqual(numpy.linalg.norm(implicit_mapped_gradient.data), 27.84296340221239, 10)
        explicit_mapped_gradient = self.explicit_shape_control.MapGradient({KratosOA.SHAPE: physical_gradient})
        self.assertAlmostEqual(numpy.linalg.norm(explicit_mapped_gradient.data), 28.075348639115425, 10)

    def test_UpdateImplicit(self):
        update_field = self.implicit_shape_control.GetEmptyField()
        update_field.data[:] = 0.1
        self.implicit_shape_control.Update(update_field)
        implicit_control_field = self.implicit_shape_control.GetControlField()
        implicit_shape_field = self.implicit_shape_control.GetPhysicalField()
        self.assertAlmostEqual(numpy.linalg.norm(implicit_control_field.data), 3.633180424916991, 10)
        self.assertAlmostEqual(numpy.linalg.norm(implicit_shape_field.data), 10.365298105786017, 10)

    def test_UpdateExplicit(self):
        update_field = self.explicit_shape_control.GetEmptyField()
        update_field.data[:] = 0.1
        self.explicit_shape_control.Update(update_field)
        explicit_control_field = self.explicit_shape_control.GetControlField()
        explicit_shape_field = self.explicit_shape_control.GetPhysicalField()
        self.assertAlmostEqual(numpy.linalg.norm(explicit_control_field.data), 3.633180424916991, 10)
        self.assertAlmostEqual(numpy.linalg.norm(explicit_shape_field.data), 22.013379276492945, 10)

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
        self.assertAlmostEqual(numpy.linalg.norm(self.implicit_shape_control.GetControlField().data), 213.27526425126743, 4)
        self.assertAlmostEqual(numpy.linalg.norm(self.explicit_shape_control.GetControlField().data), 0.0000, 4)

    def test_GetPhysicalField(self):
        shape_field = self.explicit_shape_control.GetPhysicalField()
        self.assertAlmostEqual(numpy.linalg.norm(shape_field.data), 4.24264068711927, 10)
        shape_field = self.implicit_shape_control.GetPhysicalField()
        self.assertAlmostEqual(numpy.linalg.norm(shape_field.data), 4.24264068711927, 10)

    def test_MapGradient(self):
        solid_physical_gradient = self.explicit_shape_control.GetEmptyField()
        solid_physical_gradient.data[:] = 1.0
        mapped_gradient = self.explicit_shape_control.MapGradient({KratosOA.SHAPE: solid_physical_gradient})
        self.assertAlmostEqual(numpy.linalg.norm(mapped_gradient.data), 13.52774925846869, 10)

        solid_physical_gradient = self.implicit_shape_control.GetEmptyField()
        solid_physical_gradient.data[:] = 1.0
        mapped_gradient = self.implicit_shape_control.MapGradient({KratosOA.SHAPE: solid_physical_gradient})
        self.assertAlmostEqual(numpy.linalg.norm(mapped_gradient.data), 20.108624626448545, 10)

    def test_UpdateExplicit(self):
        update_field = self.explicit_shape_control.GetEmptyField()
        update_field.data[:] = 0.1
        self.explicit_shape_control.Update(update_field)
        control_field = self.explicit_shape_control.GetControlField()
        shape_field = self.explicit_shape_control.GetPhysicalField()
        self.assertAlmostEqual(numpy.linalg.norm(control_field.data), 0.6480740698407862, 10)
        self.assertAlmostEqual(numpy.linalg.norm(shape_field.data), 4.524378410345448, 10)

    def test_UpdateImplicit(self):
        update_field = self.implicit_shape_control.GetEmptyField()
        update_field.data[:] = 0.1
        self.implicit_shape_control.Update(update_field)
        control_field = self.implicit_shape_control.GetControlField()
        shape_field = self.implicit_shape_control.GetPhysicalField()
        self.assertAlmostEqual(numpy.linalg.norm(control_field.data), 0.6480740698407862, 10)
        self.assertAlmostEqual(numpy.linalg.norm(shape_field.data), 4.370590586296503, 10)

if __name__ == "__main__":
    kratos_unittest.main()

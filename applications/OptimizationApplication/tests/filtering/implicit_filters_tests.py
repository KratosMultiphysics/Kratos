import numpy

import KratosMultiphysics as KM
from KratosMultiphysics.OptimizationApplication.filtering.helmholtz_analysis import HelmholtzAnalysis
from KratosMultiphysics.testing.utilities import ReadModelPart
import KratosMultiphysics.OptimizationApplication as KOA

# Additional imports
from KratosMultiphysics.KratosUnittest import TestCase
import KratosMultiphysics.KratosUnittest as kratos_unittest

class HelmholtzAnalysisTest(TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = KM.Model()
        cls.solid_scalar_model_part = cls.model.CreateModelPart("solid_scalar")
        cls.solid_scalar_model_part.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 3)
        cls.solid_vector_model_part = cls.model.CreateModelPart("solid_vector")
        cls.solid_vector_model_part.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 3)
        cls.solid_bulk_surface_model_part = cls.model.CreateModelPart("solid_bulk_surf")
        cls.solid_bulk_surface_model_part.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 3)

        with kratos_unittest.WorkFolderScope(".", __file__):
            solid_scalar_filter_parameters = KM.Parameters("""
                                {
                                    "solver_settings" : {
                                        "domain_size"     : 3,
                                        "echo_level"      : 0,
                                        "filter_type"     : "general_scalar",
                                        "filter_radius"     : 2.0,
                                        "model_part_name" : "solid_scalar.structure",
                                        "model_import_settings"              : {
                                            "input_type"     : "use_input_model_part"
                                        },
                                        "linear_solver_settings": {
                                            "solver_type": "skyline_lu_factorization"
                                        }
                                    },
                                    "problem_data": {
                                        "echo_level"    : 0,
                                        "start_time"    : 0.0,
                                        "end_time"      : 1.0,
                                        "parallel_type" : "OpenMP"
                                    }
                                }""")

            solid_vector_filter_parameters = KM.Parameters("""
                                {
                                    "solver_settings" : {
                                        "domain_size"     : 3,
                                        "echo_level"      : 0,
                                        "filter_type"     : "general_vector",
                                        "filter_radius"     : 2.0,
                                        "model_part_name" : "solid_vector.structure",
                                        "model_import_settings"              : {
                                            "input_type"     : "use_input_model_part"
                                        },
                                        "linear_solver_settings": {
                                            "solver_type": "skyline_lu_factorization"
                                        }
                                    },
                                    "problem_data": {
                                        "echo_level"    : 0,
                                        "start_time"    : 0.0,
                                        "end_time"      : 1.0,
                                        "parallel_type" : "OpenMP"
                                    }
                                }""")

            bulk_surface_filter_parameters = KM.Parameters("""
                                {
                                    "solver_settings" : {
                                        "domain_size"     : 3,
                                        "echo_level"      : 0,
                                        "filter_type"     : "bulk_surface_shape",
                                        "filter_radius"     : 2.0,
                                        "model_part_name" : "solid_bulk_surf.design",
                                        "model_import_settings"              : {
                                            "input_type"     : "use_input_model_part"
                                        },
                                        "linear_solver_settings": {
                                            "solver_type": "skyline_lu_factorization"
                                        }
                                    },
                                    "problem_data": {
                                        "echo_level"    : 0,
                                        "start_time"    : 0.0,
                                        "end_time"      : 1.0,
                                        "parallel_type" : "OpenMP"
                                    }
                                }""")
            shape_bulk_filter_parameters = KM.Parameters("""
                                {
                                    "solver_settings" : {
                                        "domain_size"     : 3,
                                        "echo_level"      : 0,
                                        "filter_type"     : "shape",
                                        "filter_radius"     : 2.0,
                                        "model_part_name" : "solid_bulk_surf.design",
                                        "model_import_settings"              : {
                                            "input_type"     : "use_input_model_part"
                                        },
                                        "linear_solver_settings": {
                                            "solver_type": "skyline_lu_factorization"
                                        }
                                    },
                                    "problem_data": {
                                        "echo_level"    : 0,
                                        "start_time"    : 0.0,
                                        "end_time"      : 1.0,
                                        "parallel_type" : "OpenMP"
                                    }
                                }""")
            shape_shell_vector_filter_parameters = KM.Parameters("""
                                {
                                    "solver_settings" : {
                                        "domain_size"     : 3,
                                        "echo_level"      : 0,
                                        "filter_type"     : "shape",
                                        "filter_radius"     : 500.0,
                                        "model_part_name" : "shell",
                                        "model_import_settings"              : {
                                            "input_type"     : "use_input_model_part"
                                        },
                                        "linear_solver_settings": {
                                            "solver_type": "skyline_lu_factorization"
                                        }
                                    },
                                    "problem_data": {
                                        "echo_level"    : 0,
                                        "start_time"    : 0.0,
                                        "end_time"      : 1.0,
                                        "parallel_type" : "OpenMP"
                                    }
                                }""")
            cls.solid_scalar_filter = HelmholtzAnalysis(cls.model, solid_scalar_filter_parameters)
            cls.solid_vector_filter = HelmholtzAnalysis(cls.model, solid_vector_filter_parameters)
            cls.bulk_surface_filter = HelmholtzAnalysis(cls.model, bulk_surface_filter_parameters)
            cls.shape_bulk_surface_filter = HelmholtzAnalysis(cls.model, shape_bulk_filter_parameters)
            ReadModelPart("solid",cls.solid_scalar_model_part)
            ReadModelPart("solid",cls.solid_vector_model_part)
            ReadModelPart("solid",cls.solid_bulk_surface_model_part)

            cls.shell_model_part = cls.model.CreateModelPart("shell")
            cls.shell_model_part.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 3)
            cls.shape_shell_surface_filter = HelmholtzAnalysis(cls.model, shape_shell_vector_filter_parameters)
            shell_vector_filter_parameters = KM.Parameters("""
                                {
                                    "solver_settings" : {
                                        "domain_size"     : 3,
                                        "echo_level"      : 0,
                                        "filter_type"     : "general_vector",
                                        "filter_radius"     : 500.0,
                                        "model_part_name" : "shell",
                                        "model_import_settings"              : {
                                            "input_type"     : "use_input_model_part"
                                        },
                                        "linear_solver_settings": {
                                            "solver_type": "skyline_lu_factorization"
                                        }
                                    },
                                    "problem_data": {
                                        "echo_level"    : 0,
                                        "start_time"    : 0.0,
                                        "end_time"      : 1.0,
                                        "parallel_type" : "OpenMP"
                                    }
                                }""")
            cls.shell_vector_filter = HelmholtzAnalysis(cls.model, shell_vector_filter_parameters)

            cls.shell_model_part.CreateNewNode(1, 0.0,0.0,0.0)
            cls.shell_model_part.CreateNewNode(2, 1.0,0.0,0.0)
            cls.shell_model_part.CreateNewNode(3, 2.0,0.0,0.0)

            cls.shell_model_part.CreateNewNode(4, 0.0,1.0,0.0)
            cls.shell_model_part.CreateNewNode(5, 1.0,1.0,0.0)
            cls.shell_model_part.CreateNewNode(6, 2.0,1.0,0.0)

            cls.shell_model_part.CreateNewNode(7, 0.0,2.0,0.0)
            cls.shell_model_part.CreateNewNode(8, 1.0,2.0,0.0)
            cls.shell_model_part.CreateNewNode(9, 2.0,2.0,0.0)


            properties_1 = cls.shell_model_part.CreateNewProperties(1)
            cls.shell_model_part.CreateNewElement("HelmholtzVectorSurfaceElement3D3N", 1, [1,2,4], properties_1)
            cls.shell_model_part.CreateNewElement("HelmholtzVectorSurfaceElement3D3N", 2, [2,5,4], properties_1)
            cls.shell_model_part.CreateNewElement("HelmholtzVectorSurfaceElement3D3N", 3, [2,3,5], properties_1)
            cls.shell_model_part.CreateNewElement("HelmholtzVectorSurfaceElement3D3N", 4, [3,6,5], properties_1)

            cls.shell_model_part.CreateNewElement("HelmholtzVectorSurfaceElement3D3N", 5, [4,5,7], properties_1)
            cls.shell_model_part.CreateNewElement("HelmholtzVectorSurfaceElement3D3N", 6, [5,8,7], properties_1)
            cls.shell_model_part.CreateNewElement("HelmholtzVectorSurfaceElement3D3N", 7, [5,6,8], properties_1)
            cls.shell_model_part.CreateNewElement("HelmholtzVectorSurfaceElement3D3N", 8, [6,9,8], properties_1)


            cls.closed_shell_model_part =  cls.model.CreateModelPart("closed_surface_shell")
            cls.closed_shell_model_part.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 3)
            closed_surface_shell_scalar_filter_parameters = KM.Parameters("""
                                {
                                    "solver_settings" : {
                                        "domain_size"     : 3,
                                        "echo_level"      : 0,
                                        "filter_type"     : "general_scalar",
                                        "filter_radius"     : 20.0,
                                        "model_part_name" : "closed_surface_shell",
                                        "model_import_settings"              : {
                                            "input_type"     : "use_input_model_part"
                                        },
                                        "linear_solver_settings": {
                                            "solver_type": "skyline_lu_factorization"
                                        }
                                    },
                                    "problem_data": {
                                        "echo_level"    : 0,
                                        "start_time"    : 0.0,
                                        "end_time"      : 1.0,
                                        "parallel_type" : "OpenMP"
                                    }
                                }""")

            cls.closed_shell_scalar_filter = HelmholtzAnalysis(cls.model, closed_surface_shell_scalar_filter_parameters)

            cls.closed_shell_model_part.CreateNewNode(1, 50.0,0.0,0.0)
            cls.closed_shell_model_part.CreateNewNode(2, 51.0,0.0,0.0)
            cls.closed_shell_model_part.CreateNewNode(3, 51.0,51.0,0.0)
            cls.closed_shell_model_part.CreateNewNode(4, 50.0,51.0,0.0)

            cls.closed_shell_model_part.CreateNewNode(5, 50.0,0.0,1.0)
            cls.closed_shell_model_part.CreateNewNode(6, 51.0,0.0,1.0)
            cls.closed_shell_model_part.CreateNewNode(7, 51.0,51.0,1.0)
            cls.closed_shell_model_part.CreateNewNode(8, 50.0,51.0,1.0)

            properties_1 = cls.closed_shell_model_part.CreateNewProperties(1)
            cls.closed_shell_model_part.CreateNewElement("HelmholtzVectorSurfaceElement3D4N", 1, [1,2,3,4], properties_1)
            cls.closed_shell_model_part.CreateNewElement("HelmholtzVectorSurfaceElement3D4N", 2, [5,6,7,8], properties_1)
            cls.closed_shell_model_part.CreateNewElement("HelmholtzVectorSurfaceElement3D4N", 3, [1,2,6,5], properties_1)
            cls.closed_shell_model_part.CreateNewElement("HelmholtzVectorSurfaceElement3D4N", 4, [1,4,8,5], properties_1)
            cls.closed_shell_model_part.CreateNewElement("HelmholtzVectorSurfaceElement3D4N", 5, [3,2,6,7], properties_1)
            cls.closed_shell_model_part.CreateNewElement("HelmholtzVectorSurfaceElement3D4N", 6, [4,3,7,8], properties_1)

    def setUp(self) -> None:
        for model_part_name in self.model.GetModelPartNames():
            model_part = self.model[model_part_name]
            if model_part.HasNodalSolutionStepVariable(KOA.HELMHOLTZ_VECTOR):
                KM.VariableUtils().SetHistoricalVariableToZero(KOA.HELMHOLTZ_VECTOR, model_part.Nodes)
            if model_part.HasNodalSolutionStepVariable(KOA.HELMHOLTZ_SCALAR):
                KM.VariableUtils().SetHistoricalVariableToZero(KOA.HELMHOLTZ_SCALAR, model_part.Nodes)


    def test_scalar_solid_filter(self):
        # initialization of the filter done here so that filtering radius is set.
        # Since, some of the filters use the same model part, the process info is shared
        # hence it is required to change the filter radius based on the filter type
        # to get the expected result.
        self.solid_scalar_filter.Initialize()

        unfiltered_uniform_field_nodal = KM.TensorAdaptors.VariableTensorAdaptor(self.solid_scalar_model_part.Nodes, KM.PRESSURE)
        unfiltered_uniform_field_nodal.data[:] = 1.0

        filtered_field = self.solid_scalar_filter.FilterField(unfiltered_uniform_field_nodal)
        self.assertAlmostEqual(numpy.linalg.norm(filtered_field.data), 3.741657, 4)

        for node in self.solid_scalar_model_part.Nodes:
            node.SetValue(KM.NODAL_VOLUME, 0)

        for element in self.solid_scalar_model_part.Elements:
            for node in element.GetNodes():
                node.SetValue(KM.NODAL_VOLUME, node.GetValue(KM.NODAL_VOLUME)+element.GetGeometry().Volume()/4.0)

        nodal_volume = KM.TensorAdaptors.VariableTensorAdaptor(self.solid_scalar_model_part.Nodes, KM.NODAL_VOLUME)
        nodal_volume.CollectData()

        filtered_field = self.solid_scalar_filter.FilterIntegratedField(nodal_volume)
        self.assertAlmostEqual(numpy.linalg.norm(filtered_field.data), 3.741657, 4)
        filtered_field = self.solid_scalar_filter.UnFilterField(unfiltered_uniform_field_nodal)
        self.assertAlmostEqual(numpy.linalg.norm(filtered_field.data), 3.741657, 4)

        self.solid_scalar_filter.Finalize()

    def test_vector_solid_filter(self):
        # initialization of the filter done here so that filtering radius is set.
        # Since, some of the filters use the same model part, the process info is shared
        # hence it is required to change the filter radius based on the filter type
        # to get the expected result.
        self.solid_vector_filter.Initialize()

        unfiltered_uniform_field_nodal = KM.TensorAdaptors.VariableTensorAdaptor(self.solid_vector_model_part.Nodes, KM.VELOCITY)
        unfiltered_uniform_field_nodal.data[:] = 1.0

        filtered_field = self.solid_vector_filter.FilterField(unfiltered_uniform_field_nodal)
        self.assertAlmostEqual(numpy.linalg.norm(filtered_field.data), 6.48074, 4)

        for node in self.solid_vector_model_part.Nodes:
            node.SetValue(KM.VELOCITY, KM.Array3([0, 0, 0]))

        for element in self.solid_vector_model_part.Elements:
            for node in element.GetNodes():
                node.SetValue(KM.VELOCITY_X, node.GetValue(KM.VELOCITY_X)+element.GetGeometry().Volume()/4.0)
                node.SetValue(KM.VELOCITY_Y, node.GetValue(KM.VELOCITY_Y)+element.GetGeometry().Volume()/4.0)
                node.SetValue(KM.VELOCITY_Z, node.GetValue(KM.VELOCITY_Z)+element.GetGeometry().Volume()/4.0)

        nodal_volume = KM.TensorAdaptors.VariableTensorAdaptor(self.solid_vector_model_part.Nodes, KM.VELOCITY)
        nodal_volume.CollectData()

        filtered_field = self.solid_vector_filter.FilterIntegratedField(nodal_volume)
        self.assertAlmostEqual(numpy.linalg.norm(filtered_field.data), 6.48074, 4)
        filtered_field = self.solid_vector_filter.UnFilterField(unfiltered_uniform_field_nodal)
        self.assertAlmostEqual(numpy.linalg.norm(filtered_field.data), 6.48074, 4)

    def test_scalar_closed_shell_filter(self):
        # initialization of the filter done here so that filtering radius is set.
        # Since, some of the filters use the same model part, the process info is shared
        # hence it is required to change the filter radius based on the filter type
        # to get the expected result.
        self.closed_shell_scalar_filter.Initialize()

        unfiltered_uniform_field_nodal = KM.TensorAdaptors.VariableTensorAdaptor(self.closed_shell_model_part.Nodes, KM.PRESSURE)
        unfiltered_uniform_field_nodal.data[:] = 1.0

        filtered_field = self.closed_shell_scalar_filter.FilterField(unfiltered_uniform_field_nodal)
        self.assertAlmostEqual(numpy.linalg.norm(filtered_field.data), 2.8284271247, 4)

        for node in self.closed_shell_model_part.Nodes:
            node.SetValue(KM.NODAL_AREA, 0)

        for element in self.closed_shell_model_part.Elements:
            for node in element.GetNodes():
                node.SetValue(KM.NODAL_AREA, node.GetValue(KM.NODAL_AREA)+element.GetGeometry().Area()/4.0)

        nodal_area = KM.TensorAdaptors.VariableTensorAdaptor(self.closed_shell_model_part.Nodes, KM.NODAL_AREA)
        nodal_area.CollectData()

        filtered_field = self.closed_shell_scalar_filter.FilterIntegratedField(nodal_area)
        self.assertAlmostEqual(numpy.linalg.norm(filtered_field.data), 2.8284271247, 4)
        filtered_field = self.closed_shell_scalar_filter.UnFilterField(unfiltered_uniform_field_nodal)
        self.assertAlmostEqual(numpy.linalg.norm(filtered_field.data), 2.8284271247, 4)

    def test_vector_shell_filter(self):
        # initialization of the filter done here so that filtering radius is set.
        # Since, some of the filters use the same model part, the process info is shared
        # hence it is required to change the filter radius based on the filter type
        # to get the expected result.
        self.shell_vector_filter.Initialize()

        unfiltered_uniform_field_nodal = KM.TensorAdaptors.VariableTensorAdaptor(self.shell_model_part.Nodes, KM.VELOCITY)
        unfiltered_uniform_field_nodal.data[:] = 1.0

        filtered_field = self.shell_vector_filter.UnFilterField(unfiltered_uniform_field_nodal)
        filtered_field = self.shell_vector_filter.FilterField(filtered_field)
        self.assertAlmostEqual(numpy.linalg.norm(filtered_field.data), 5.196153341144455, 4)

        filtered_field = self.shell_vector_filter.FilterField(unfiltered_uniform_field_nodal)
        self.assertAlmostEqual(numpy.linalg.norm(filtered_field.data), 5.1961524, 4)

        for node in self.shell_model_part.Nodes:
            node.SetValue(KM.VELOCITY, KM.Array3([0, 0, 0]))

        for element in self.shell_model_part.Elements:
            for node in element.GetNodes():
                node.SetValue(KM.VELOCITY_X, node.GetValue(KM.VELOCITY_X)+element.GetGeometry().Area()/3.0)
                node.SetValue(KM.VELOCITY_Y, node.GetValue(KM.VELOCITY_Y)+element.GetGeometry().Area()/3.0)
                node.SetValue(KM.VELOCITY_Z, node.GetValue(KM.VELOCITY_Z)+element.GetGeometry().Area()/3.0)

        nodal_area = KM.TensorAdaptors.VariableTensorAdaptor(self.shell_model_part.Nodes, KM.VELOCITY)
        nodal_area.CollectData()

        filtered_field = self.shell_vector_filter.FilterIntegratedField(nodal_area)
        self.assertAlmostEqual(numpy.linalg.norm(filtered_field.data), 5.196118198546968, 4)

    def test_bulk_surface_shape(self):
        # initialization of the filter done here so that filtering radius is set.
        # Since, some of the filters use the same model part, the process info is shared
        # hence it is required to change the filter radius based on the filter type
        # to get the expected result.
        self.bulk_surface_filter.Initialize()

        unfiltered_uniform_field_nodal = KM.TensorAdaptors.VariableTensorAdaptor(self.solid_bulk_surface_model_part.Nodes, KM.VELOCITY)
        unfiltered_uniform_field_nodal.data[:] = 1.0

        filtered_field = self.bulk_surface_filter.FilterField(unfiltered_uniform_field_nodal)
        self.assertAlmostEqual(numpy.linalg.norm(filtered_field.data), 6.4807406, 4)

        for node in self.solid_bulk_surface_model_part.Nodes:
            node.SetValue(KM.VELOCITY, KM.Array3([0, 0, 0]))

        for element in self.solid_bulk_surface_model_part.Elements:
            for node in element.GetNodes():
                node.SetValue(KM.VELOCITY_X, node.GetValue(KM.VELOCITY_X)+element.GetGeometry().Volume()/4.0)
                node.SetValue(KM.VELOCITY_Y, node.GetValue(KM.VELOCITY_Y)+element.GetGeometry().Volume()/4.0)
                node.SetValue(KM.VELOCITY_Z, node.GetValue(KM.VELOCITY_Z)+element.GetGeometry().Volume()/4.0)

        nodal_area = KM.TensorAdaptors.VariableTensorAdaptor(self.solid_bulk_surface_model_part.Nodes, KM.VELOCITY)
        nodal_area.CollectData()

        filtered_field = self.bulk_surface_filter.FilterIntegratedField(nodal_area)
        self.assertAlmostEqual(numpy.linalg.norm(filtered_field.data), 6.4807406, 4)
        filtered_field = self.bulk_surface_filter.UnFilterField(unfiltered_uniform_field_nodal)
        self.assertAlmostEqual(numpy.linalg.norm(filtered_field.data), 6.4807406, 4)

    def test_shape_surface_shape(self):
        self.shape_bulk_surface_filter.Initialize()

        unfiltered_uniform_field_nodal = KM.TensorAdaptors.VariableTensorAdaptor(self.solid_bulk_surface_model_part.Nodes, KM.VELOCITY)
        unfiltered_uniform_field_nodal.data[:] = 1.0

        filtered_field = self.shape_bulk_surface_filter.FilterField(unfiltered_uniform_field_nodal)
        self.assertAlmostEqual(numpy.linalg.norm(filtered_field.data), 6.4807406, 4)

        self.shape_shell_surface_filter.Initialize()

        unfiltered_uniform_field_nodal = KM.TensorAdaptors.VariableTensorAdaptor(self.shell_model_part.Nodes, KM.VELOCITY)
        unfiltered_uniform_field_nodal.data[:] = 1.0

        filtered_field = self.shape_shell_surface_filter.UnFilterField(unfiltered_uniform_field_nodal)
        filtered_field = self.shape_shell_surface_filter.FilterField(filtered_field)
        self.assertAlmostEqual(numpy.linalg.norm(filtered_field.data), 5.196153341144455, 4)

        filtered_field = self.shape_shell_surface_filter.FilterField(unfiltered_uniform_field_nodal)
        self.assertAlmostEqual(numpy.linalg.norm(filtered_field.data), 5.1961524, 4)

if __name__ == '__main__':
    KM.KratosUnittest.main()

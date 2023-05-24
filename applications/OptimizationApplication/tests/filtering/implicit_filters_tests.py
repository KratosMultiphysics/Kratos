import KratosMultiphysics as KM
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.OptimizationApplication.filtering.helmholtz_analysis import HelmholtzAnalysis
from KratosMultiphysics.testing.utilities import ReadModelPart
import KratosMultiphysics.OptimizationApplication as KOA
from numpy import linalg as LA
import os

# Additional imports
from KratosMultiphysics.KratosUnittest import TestCase

class HelmholtzAnalysisTest(TestCase):

    def test_solid(self):
        with KM.KratosUnittest.WorkFolderScope(".", __file__):
            model = KM.Model()
            model_part = model.CreateModelPart("solid")
            model_part.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 3)
            model_part_io = KM.ModelPartIO("solid")

            general_scalar_parameters = KM.Parameters("""
                            {
                                "solver_settings" : {
                                    "domain_size"     : 3,
                                    "echo_level"      : 0,
                                    "filter_type"     : "general_scalar",
                                    "filter_radius"     : 2.0,
                                    "model_part_name" : "solid",
                                    "time_stepping" : {
                                        "time_step"       : 1.0
                                    },
                                    "model_import_settings"              : {
                                        "input_type"     : "use_input_model_part"
                                    },
                                    "linear_solver_settings" : {
                                        "solver_type" : "amgcl",
                                        "smoother_type":"ilu0",
                                        "krylov_type": "gmres",
                                        "coarsening_type": "aggregation",
                                        "max_iteration": 200,
                                        "verbosity" : 0,
                                        "tolerance": 1e-7
                                    }
                                },
                                "processes" : {
                                    "boundary_conditions_process_list" : [
                                    ]
                                },
                                "problem_data": {
                                    "echo_level"    : 0,
                                    "start_time"    : 0.0,
                                    "end_time"      : 1.0,
                                    "parallel_type" : "OpenMP"
                                }
                            }""")

        my_filter = HelmholtzAnalysis(model, general_scalar_parameters)

        model_part_io.ReadModelPart(model_part)

        unfiltered_uniform_field_nodal = KM.ContainerExpression.NodalNonHistoricalExpression(model_part)
        unfiltered_uniform_field_nodal.SetData(1.0)

        filtered_field = my_filter.FilterField(unfiltered_uniform_field_nodal)
        self.assertAlmostEqual(LA.norm(filtered_field.Evaluate()), 3.741657, 4)

        for node in model_part.Nodes:
            node.SetValue(KM.NODAL_VOLUME, 0)

        for element in model_part.Elements:
            for node in element.GetNodes():
                node.SetValue(KM.NODAL_VOLUME, node.GetValue(KM.NODAL_VOLUME)+element.GetGeometry().Volume()/4.0)

        nodal_volume = KM.ContainerExpression.NodalNonHistoricalExpression(model_part)
        nodal_volume.Read(KM.NODAL_VOLUME)

        filtered_field = my_filter.FilterIntegratedField(nodal_volume)
        self.assertAlmostEqual(LA.norm(filtered_field.Evaluate()), 3.741657, 4)
        filtered_field = my_filter.UnFilterField(unfiltered_uniform_field_nodal)
        self.assertAlmostEqual(LA.norm(filtered_field.Evaluate()), 3.741657, 4)

        my_filter.Finalize()

    def test_vector_solid(self):
        with KM.KratosUnittest.WorkFolderScope(".", __file__):
            model = KM.Model()
            model_part = model.CreateModelPart("solid")
            model_part.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 3)
            model_part_io = KM.ModelPartIO("solid")

            filter_parameters = KM.Parameters("""
                            {
                                "solver_settings" : {
                                    "domain_size"     : 3,
                                    "echo_level"      : 0,
                                    "filter_type"     : "general_vector",
                                    "filter_radius"     : 2.0,
                                    "model_part_name" : "solid",
                                    "time_stepping" : {
                                        "time_step"       : 1.0
                                    },
                                    "model_import_settings"              : {
                                        "input_type"     : "use_input_model_part"
                                    },
                                    "linear_solver_settings" : {
                                        "solver_type" : "amgcl",
                                        "smoother_type":"ilu0",
                                        "krylov_type": "gmres",
                                        "coarsening_type": "aggregation",
                                        "max_iteration": 200,
                                        "verbosity" : 0,
                                        "tolerance": 1e-7
                                    }
                                },
                                "processes" : {
                                    "boundary_conditions_process_list" : [
                                    ]
                                },
                                "problem_data": {
                                    "echo_level"    : 0,
                                    "start_time"    : 0.0,
                                    "end_time"      : 1.0,
                                    "parallel_type" : "OpenMP"
                                }
                            }""")

        my_filter = HelmholtzAnalysis(model, filter_parameters)

        model_part_io.ReadModelPart(model_part)

        unfiltered_uniform_field_nodal = KM.ContainerExpression.NodalNonHistoricalExpression(model_part)
        unfiltered_uniform_field_nodal.SetData(KM.Array3([1, 1, 1]))

        filtered_field = my_filter.FilterField(unfiltered_uniform_field_nodal)
        self.assertAlmostEqual(LA.norm(filtered_field.Evaluate()), 6.48074, 4)

        for node in model_part.Nodes:
            node.SetValue(KM.VELOCITY, KM.Array3([0, 0, 0]))

        for element in model_part.Elements:
            for node in element.GetNodes():
                node.SetValue(KM.VELOCITY_X, node.GetValue(KM.VELOCITY_X)+element.GetGeometry().Volume()/4.0)
                node.SetValue(KM.VELOCITY_Y, node.GetValue(KM.VELOCITY_Y)+element.GetGeometry().Volume()/4.0)
                node.SetValue(KM.VELOCITY_Z, node.GetValue(KM.VELOCITY_Z)+element.GetGeometry().Volume()/4.0)

        nodal_volume = KM.ContainerExpression.NodalNonHistoricalExpression(model_part)
        nodal_volume.Read(KM.VELOCITY)

        filtered_field = my_filter.FilterIntegratedField(nodal_volume)
        self.assertAlmostEqual(LA.norm(filtered_field.Evaluate()), 6.48074, 4)
        filtered_field = my_filter.UnFilterField(unfiltered_uniform_field_nodal)
        self.assertAlmostEqual(LA.norm(filtered_field.Evaluate()), 6.48074, 4)

        my_filter.Finalize()

    def test_shell(self):
        with KM.KratosUnittest.WorkFolderScope(".", __file__):
            model = KM.Model()
            model_part = model.CreateModelPart("shell")
            model_part.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 3)

            general_scalar_parameters = KM.Parameters("""
                            {
                                "solver_settings" : {
                                    "domain_size"     : 3,
                                    "echo_level"      : 0,
                                    "filter_type"     : "general_scalar",
                                    "filter_radius"     : 0.2,
                                    "model_part_name" : "shell",
                                    "model_import_settings"              : {
                                        "input_type"     : "use_input_model_part"
                                    },
                                    "linear_solver_settings" : {
                                        "solver_type" : "amgcl",
                                        "smoother_type":"ilu0",
                                        "krylov_type": "gmres",
                                        "coarsening_type": "aggregation",
                                        "max_iteration": 200,
                                        "verbosity" : 0,
                                        "tolerance": 1e-7
                                    }
                                },
                                "processes" : {
                                    "boundary_conditions_process_list" : [
                                    ]
                                },
                                "problem_data": {
                                    "echo_level"    : 0,
                                    "start_time"    : 0.0,
                                    "end_time"      : 1.0,
                                    "parallel_type" : "OpenMP"
                                }
                            }""")

        my_filter = HelmholtzAnalysis(model, general_scalar_parameters)

        model_part.CreateNewNode(1, 0.0,0.0,0.0)
        model_part.CreateNewNode(2, 1.0,0.0,0.0)
        model_part.CreateNewNode(3, 2.0,0.0,0.0)

        model_part.CreateNewNode(4, 0.0,1.0,0.0)
        model_part.CreateNewNode(5, 1.0,1.0,0.0)
        model_part.CreateNewNode(6, 2.0,1.0,0.0)

        model_part.CreateNewNode(7, 0.0,2.0,0.0)
        model_part.CreateNewNode(8, 1.0,2.0,0.0)
        model_part.CreateNewNode(9, 2.0,2.0,0.0)


        properties_1 = model_part.CreateNewProperties(1)
        model_part.CreateNewElement("HelmholtzVectorSurfaceElement3D3N", 1, [1,2,4], properties_1)
        model_part.CreateNewElement("HelmholtzVectorSurfaceElement3D3N", 2, [2,5,4], properties_1)
        model_part.CreateNewElement("HelmholtzVectorSurfaceElement3D3N", 3, [2,3,5], properties_1)
        model_part.CreateNewElement("HelmholtzVectorSurfaceElement3D3N", 4, [3,6,5], properties_1)

        model_part.CreateNewElement("HelmholtzVectorSurfaceElement3D3N", 5, [4,5,7], properties_1)
        model_part.CreateNewElement("HelmholtzVectorSurfaceElement3D3N", 6, [5,8,7], properties_1)
        model_part.CreateNewElement("HelmholtzVectorSurfaceElement3D3N", 7, [5,6,8], properties_1)
        model_part.CreateNewElement("HelmholtzVectorSurfaceElement3D3N", 8, [6,9,8], properties_1)

        unfiltered_uniform_field_nodal = KM.ContainerExpression.NodalNonHistoricalExpression(model_part)
        unfiltered_uniform_field_nodal.SetData(1.0)

        filtered_field = my_filter.FilterField(unfiltered_uniform_field_nodal)
        self.assertAlmostEqual(LA.norm(filtered_field.Evaluate()), 3.0, 4)

        for node in model_part.Nodes:
            node.SetValue(KM.NODAL_AREA, 0)

        for element in model_part.Elements:
            for node in element.GetNodes():
                node.SetValue(KM.NODAL_AREA, node.GetValue(KM.NODAL_AREA)+element.GetGeometry().Area()/3.0)

        nodal_area = KM.ContainerExpression.NodalNonHistoricalExpression(model_part)
        nodal_area.Read(KM.NODAL_AREA)

        filtered_field = my_filter.FilterIntegratedField(nodal_area)
        self.assertAlmostEqual(LA.norm(filtered_field.Evaluate()), 3.0, 4)
        filtered_field = my_filter.UnFilterField(unfiltered_uniform_field_nodal)
        self.assertAlmostEqual(LA.norm(filtered_field.Evaluate()), 3.0, 4)

        my_filter.Finalize()

    def test_vector_shell(self):
        with KM.KratosUnittest.WorkFolderScope(".", __file__):
            model = KM.Model()
            model_part = model.CreateModelPart("shell")
            model_part.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 3)
            filter_parameters = KM.Parameters("""
                            {
                                "solver_settings" : {
                                    "domain_size"     : 3,
                                    "echo_level"      : 0,
                                    "filter_type"     : "general_vector",
                                    "filter_radius"     : 0.2,
                                    "model_part_name" : "shell",
                                    "model_import_settings"              : {
                                        "input_type"     : "use_input_model_part"
                                    },
                                    "linear_solver_settings" : {
                                        "solver_type" : "amgcl",
                                        "smoother_type":"ilu0",
                                        "krylov_type": "gmres",
                                        "coarsening_type": "aggregation",
                                        "max_iteration": 200,
                                        "verbosity" : 0,
                                        "tolerance": 1e-7
                                    }
                                },
                                "processes" : {
                                    "boundary_conditions_process_list" : [
                                    ]
                                },
                                "problem_data": {
                                    "echo_level"    : 0,
                                    "start_time"    : 0.0,
                                    "end_time"      : 1.0,
                                    "parallel_type" : "OpenMP"
                                }
                            }""")

        my_filter = HelmholtzAnalysis(model, filter_parameters)

        model_part.CreateNewNode(1, 0.0,0.0,0.0)
        model_part.CreateNewNode(2, 1.0,0.0,0.0)
        model_part.CreateNewNode(3, 2.0,0.0,0.0)

        model_part.CreateNewNode(4, 0.0,1.0,0.0)
        model_part.CreateNewNode(5, 1.0,1.0,0.0)
        model_part.CreateNewNode(6, 2.0,1.0,0.0)

        model_part.CreateNewNode(7, 0.0,2.0,0.0)
        model_part.CreateNewNode(8, 1.0,2.0,0.0)
        model_part.CreateNewNode(9, 2.0,2.0,0.0)


        properties_1 = model_part.CreateNewProperties(1)
        model_part.CreateNewElement("HelmholtzVectorSurfaceElement3D4N", 1, [1,2,5,4], properties_1)
        model_part.CreateNewElement("HelmholtzVectorSurfaceElement3D4N", 2, [2,3,6,5], properties_1)
        model_part.CreateNewElement("HelmholtzVectorSurfaceElement3D4N", 3, [4,5,8,7], properties_1)
        model_part.CreateNewElement("HelmholtzVectorSurfaceElement3D4N", 4, [5,6,9,8], properties_1)

        unfiltered_uniform_field_nodal = KM.ContainerExpression.NodalNonHistoricalExpression(model_part)
        unfiltered_uniform_field_nodal.SetData(KM.Array3([1, 1, 1]))

        filtered_field = my_filter.FilterField(unfiltered_uniform_field_nodal)
        self.assertAlmostEqual(LA.norm(filtered_field.Evaluate()), 5.1961524, 4)


        for node in model_part.Nodes:
            node.SetValue(KM.VELOCITY, KM.Array3([0, 0, 0]))

        for element in model_part.Elements:
            for node in element.GetNodes():
                node.SetValue(KM.VELOCITY_X, node.GetValue(KM.VELOCITY_X)+element.GetGeometry().Area()/4.0)
                node.SetValue(KM.VELOCITY_Y, node.GetValue(KM.VELOCITY_Y)+element.GetGeometry().Area()/4.0)
                node.SetValue(KM.VELOCITY_Z, node.GetValue(KM.VELOCITY_Z)+element.GetGeometry().Area()/4.0)

        nodal_area = KM.ContainerExpression.NodalNonHistoricalExpression(model_part)
        nodal_area.Read(KM.VELOCITY)


        filtered_field = my_filter.FilterIntegratedField(nodal_area)
        self.assertAlmostEqual(LA.norm(filtered_field.Evaluate()), 5.1961524, 4)
        filtered_field = my_filter.UnFilterField(unfiltered_uniform_field_nodal)
        self.assertAlmostEqual(LA.norm(filtered_field.Evaluate()), 5.1961524, 4)

        my_filter.Finalize()

    def test_bulk_surface_shape(self):
        with KM.KratosUnittest.WorkFolderScope(".", __file__):
            model = KM.Model()
            model_part = model.CreateModelPart("bulk_surface")
            model_part.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 3)
            model_part_io = KM.ModelPartIO("solid")
            filter_parameters = KM.Parameters("""
                            {
                                "solver_settings" : {
                                    "domain_size"     : 3,
                                    "echo_level"      : 0,
                                    "filter_type"     : "bulk_surface_shape",
                                    "filter_radius"     : 2.0,
                                    "model_part_name" : "bulk_surface",
                                    "model_import_settings"              : {
                                        "input_type"     : "use_input_model_part"
                                    },
                                    "linear_solver_settings" : {
                                        "solver_type" : "amgcl",
                                        "smoother_type":"ilu0",
                                        "krylov_type": "gmres",
                                        "coarsening_type": "aggregation",
                                        "max_iteration": 200,
                                        "verbosity" : 0,
                                        "tolerance": 1e-7
                                    }
                                },
                                "processes" : {
                                    "boundary_conditions_process_list" : [
                                    ]
                                },
                                "problem_data": {
                                    "echo_level"    : 0,
                                    "start_time"    : 0.0,
                                    "end_time"      : 1.0,
                                    "parallel_type" : "OpenMP"
                                }
                            }""")

        my_filter = HelmholtzAnalysis(model, filter_parameters)
        model_part_io.ReadModelPart(model_part)

        unfiltered_uniform_field_nodal = KM.ContainerExpression.NodalNonHistoricalExpression(model_part)
        unfiltered_uniform_field_nodal.SetData(KM.Array3([1, 1, 1]))

        filtered_field = my_filter.FilterField(unfiltered_uniform_field_nodal)
        self.assertAlmostEqual(LA.norm(filtered_field.Evaluate()), 6.4807406, 4)


        for node in model_part.Nodes:
            node.SetValue(KM.VELOCITY, KM.Array3([0, 0, 0]))

        for element in model_part.Elements:
            for node in element.GetNodes():
                node.SetValue(KM.VELOCITY_X, node.GetValue(KM.VELOCITY_X)+element.GetGeometry().Volume()/4.0)
                node.SetValue(KM.VELOCITY_Y, node.GetValue(KM.VELOCITY_Y)+element.GetGeometry().Volume()/4.0)
                node.SetValue(KM.VELOCITY_Z, node.GetValue(KM.VELOCITY_Z)+element.GetGeometry().Volume()/4.0)

        nodal_area = KM.ContainerExpression.NodalNonHistoricalExpression(model_part)
        nodal_area.Read(KM.VELOCITY)

        filtered_field = my_filter.FilterIntegratedField(nodal_area)
        self.assertAlmostEqual(LA.norm(filtered_field.Evaluate()), 6.4807406, 4)
        filtered_field = my_filter.UnFilterField(unfiltered_uniform_field_nodal)
        self.assertAlmostEqual(LA.norm(filtered_field.Evaluate()), 6.4807406, 4)

        my_filter.Finalize()


if __name__ == '__main__':
    KM.KratosUnittest.main()

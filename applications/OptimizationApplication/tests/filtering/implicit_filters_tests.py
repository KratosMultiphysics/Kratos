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

    def test_solid_element(self):
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
        self.assertAlmostEqual(LA.norm(filtered_field.Evaluate()), 16.97056, 4)

        for node in model_part.Nodes:
            node.SetValue(KM.NODAL_VOLUME, 0)

        for element in model_part.Elements:
            for node in element.GetNodes():
                node.SetValue(KM.NODAL_VOLUME, node.GetValue(KM.NODAL_VOLUME)+element.GetGeometry().Volume()/8.0)

        nodal_volume = KM.ContainerExpression.NodalNonHistoricalExpression(model_part)
        nodal_volume.Read(KM.NODAL_VOLUME)

        filtered_field = my_filter.FilterIntegratedField(nodal_volume)
        self.assertAlmostEqual(LA.norm(filtered_field.Evaluate()), 16.97056, 4)
        filtered_field = my_filter.UnFilterField(unfiltered_uniform_field_nodal)
        self.assertAlmostEqual(LA.norm(filtered_field.Evaluate()), 16.97056, 4)

        my_filter.Finalize()

    def test_shell_element(self):
        with KM.KratosUnittest.WorkFolderScope(".", __file__):
            model = KM.Model()
            model_part = model.CreateModelPart("shell")
            model_part.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 3)
            model_part_io = KM.ModelPartIO("shell")

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

        model_part_io.ReadModelPart(model_part)

        unfiltered_uniform_field_nodal = KM.ContainerExpression.NodalNonHistoricalExpression(model_part)
        unfiltered_uniform_field_nodal.SetData(1.0)

        filtered_field = my_filter.FilterField(unfiltered_uniform_field_nodal)
        self.assertAlmostEqual(LA.norm(filtered_field.Evaluate()), 20.97617, 4)
        filtered_field = my_filter.FilterIntegratedField(unfiltered_uniform_field_nodal)
        self.assertAlmostEqual(LA.norm(filtered_field.Evaluate()), 2485.48679, 4)
        filtered_field = my_filter.UnFilterField(unfiltered_uniform_field_nodal)
        self.assertAlmostEqual(LA.norm(filtered_field.Evaluate()), 20.97617, 4)

        my_filter.Finalize()

if __name__ == '__main__':
    KM.KratosUnittest.main()

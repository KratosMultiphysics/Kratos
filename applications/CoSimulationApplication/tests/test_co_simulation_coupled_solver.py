from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

import KratosMultiphysics as KM
import KratosMultiphysics.KratosUnittest as KratosUnittest

from KratosMultiphysics.CoSimulationApplication.co_simulation_analysis import CoSimulationAnalysis
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_solver_wrapper import CoSimulationSolverWrapper
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_coupled_solver import CoSimulationCoupledSolver
from KratosMultiphysics.CoSimulationApplication.helpers.dummy_solver_wrapper import DummySolverWrapper

class TestCoupledSolverGetSolver(KratosUnittest.TestCase):

    def test_GetSolver_zero_levels(self):
        # CoSim contains only one SolverWrapper (not really CoSimulation...)
        params = KM.Parameters("""{
            "problem_data" : {
                "problem_name" : "my_dummy_solver",
                "start_time" : 0.0,
                "end_time" : 3.0,
                "echo_level" : 0,
                "print_colors" : true,
                "parallel_type" : "OpenMP"
            },
            "solver_settings" : {
                "type"        : "helpers.dummy_solver_wrapper",
                "solver_wrapper_settings": {
                    "time_step" : 1.0,
                    "domain_size" : 1,
                    "main_model_part_name" : "dummy"
                }
            }
        }""")

        co_sim_analysis = CoSimulationAnalysis(params)

        top_level_solver_wrapper = co_sim_analysis._GetSolver()

        self.assertIsInstance(top_level_solver_wrapper, CoSimulationSolverWrapper)
        self.assertEqual(top_level_solver_wrapper.name, "my_dummy_solver")

    def test_GetSolver_one_level(self):
        params = KM.Parameters("""{
            "problem_data" : {
                "problem_name" : "coupling_solver_top_level",
                "start_time" : 0.0,
                "end_time" : 3.0,
                "echo_level" : 0,
                "print_colors" : true,
                "parallel_type" : "OpenMP"
            },
            "solver_settings" : {
                "type"        : "coupled_solvers.gauss_seidel_weak",
                "coupling_sequence" : [{
                    "name": "structure",
                    "output_data_list": [],
                    "input_data_list": []
                },{
                    "name": "aux_structure",
                    "output_data_list": [],
                    "input_data_list": []
                },{
                    "name": "fluid",
                    "output_data_list": [],
                    "input_data_list": []
                }],
                "solvers" : {
                    "fluid": {
                        "type"        : "helpers.dummy_solver_wrapper",
                        "solver_wrapper_settings": {
                            "time_step" : 1.0,
                            "domain_size" : 1,
                            "main_model_part_name" : "dummy"
                        }
                    },
                    "structure": {
                        "type"        : "helpers.dummy_solver_wrapper",
                        "solver_wrapper_settings": {
                            "time_step" : 1.0,
                            "domain_size" : 1,
                            "main_model_part_name" : "dummy"
                        }
                    },
                    "aux_structure": {
                        "type"        : "helpers.dummy_solver_wrapper",
                        "solver_wrapper_settings": {
                            "time_step" : 1.0,
                            "domain_size" : 1,
                            "main_model_part_name" : "dummy"
                        }
                    }
                }
            }
        }""")

        co_sim_analysis = CoSimulationAnalysis(params)

        top_level_solver_wrapper  = co_sim_analysis._GetSolver()
        fluid_solver_wrapper      = co_sim_analysis._GetSolver("fluid")
        structural_solver_wrapper = co_sim_analysis._GetSolver("structure")
        aux_structural_solver_wrapper = co_sim_analysis._GetSolver("aux_structure")

        self.assertIsInstance(top_level_solver_wrapper, CoSimulationCoupledSolver)
        self.assertIsInstance(fluid_solver_wrapper, CoSimulationSolverWrapper)
        self.assertIsInstance(structural_solver_wrapper, CoSimulationSolverWrapper)
        self.assertIsInstance(aux_structural_solver_wrapper, CoSimulationSolverWrapper)

        self.assertEqual(top_level_solver_wrapper.name, "coupling_solver_top_level")
        self.assertEqual(fluid_solver_wrapper.name, "fluid")
        self.assertEqual(structural_solver_wrapper.name, "structure")
        self.assertEqual(aux_structural_solver_wrapper.name, "aux_structure")

    def test_GetSolver_two_levels(self):
        params = KM.Parameters("""{
            "problem_data" : {
                "problem_name" : "coupling_solver_top_level",
                "start_time" : 0.0,
                "end_time" : 3.0,
                "echo_level" : 0,
                "print_colors" : true,
                "parallel_type" : "OpenMP"
            },
            "solver_settings" : {
                "type"        : "coupled_solvers.gauss_seidel_weak",
                "solver_wrapper_settings": {
                    "time_step" : 1.0,
                    "domain_size" : 1,
                    "main_model_part_name" : "dummy"
                },
                "coupling_sequence" : [{
                    "name": "fsi",
                    "output_data_list": [],
                    "input_data_list": []
                },{
                    "name": "controller",
                    "output_data_list": [],
                    "input_data_list": []
                }],
                "solvers" : {
                    "fsi": {
                        "type"        : "coupled_solvers.gauss_seidel_weak",
                        "coupling_sequence" : [{
                            "name": "structure",
                            "output_data_list": [],
                            "input_data_list": []
                        },{
                            "name": "fluid",
                            "output_data_list": [],
                            "input_data_list": []
                        }],
                        "solvers" : {
                            "fluid": {
                                "type"        : "helpers.dummy_solver_wrapper",
                                "solver_wrapper_settings": {
                                    "time_step" : 1.0,
                                    "domain_size" : 1,
                                    "main_model_part_name" : "dummy"
                                }
                            },
                            "structure": {
                                "type"        : "helpers.dummy_solver_wrapper",
                                "solver_wrapper_settings": {
                                    "time_step" : 1.0,
                                    "domain_size" : 1,
                                    "main_model_part_name" : "dummy"
                                }
                            }
                        }
                    },
                    "controller": {
                        "type"        : "helpers.dummy_solver_wrapper",
                        "solver_wrapper_settings": {
                            "time_step" : 1.0,
                            "domain_size" : 1,
                            "main_model_part_name" : "dummy"
                        }
                    }
                }
            }
        }""")

        co_sim_analysis = CoSimulationAnalysis(params)

        top_level_solver_wrapper  = co_sim_analysis._GetSolver()
        controller_solver_wrapper = co_sim_analysis._GetSolver("controller")
        fsi_solver_wrapper        = co_sim_analysis._GetSolver("fsi")
        fluid_solver_wrapper      = co_sim_analysis._GetSolver("fsi.fluid")
        structural_solver_wrapper = co_sim_analysis._GetSolver("fsi.structure")

        fsi_fluid_solver_wrapper      = fsi_solver_wrapper._GetSolver("fluid")
        fsi_structural_solver_wrapper = fsi_solver_wrapper._GetSolver("structure")

        self.assertIsInstance(top_level_solver_wrapper, CoSimulationCoupledSolver)
        self.assertIsInstance(controller_solver_wrapper, CoSimulationSolverWrapper)
        self.assertIsInstance(fsi_solver_wrapper, CoSimulationCoupledSolver)
        self.assertIsInstance(fluid_solver_wrapper, CoSimulationSolverWrapper)
        self.assertIsInstance(structural_solver_wrapper, CoSimulationSolverWrapper)
        self.assertIsInstance(fsi_fluid_solver_wrapper, CoSimulationSolverWrapper)
        self.assertIsInstance(fsi_structural_solver_wrapper, CoSimulationSolverWrapper)

        self.assertEqual(top_level_solver_wrapper.name, "coupling_solver_top_level")

        self.assertEqual(controller_solver_wrapper.name, "controller")
        self.assertEqual(fsi_solver_wrapper.name, "fsi")

        self.assertEqual(fluid_solver_wrapper.name, "fluid")
        self.assertEqual(structural_solver_wrapper.name, "structure")

        self.assertEqual(fsi_fluid_solver_wrapper.name, "fluid")
        self.assertEqual(fsi_structural_solver_wrapper.name, "structure")

        self.assertIs(fluid_solver_wrapper, fsi_fluid_solver_wrapper)
        self.assertIs(structural_solver_wrapper, fsi_structural_solver_wrapper)


if __name__ == '__main__':
    KratosUnittest.main()

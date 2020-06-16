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


class TestCoupledSolverModelAccess(KratosUnittest.TestCase):
    def test_model_access_one_level(self):
        params = KM.Parameters("""{
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
                        "main_model_part_name" : "fl_model_part"
                    }
                },
                "structure": {
                    "type"        : "helpers.dummy_solver_wrapper",
                    "solver_wrapper_settings": {
                        "time_step" : 1.0,
                        "domain_size" : 1,
                        "main_model_part_name" : "strct_model_part"
                    }
                }
            }
        }""")

        coupled_solver = CoSimulationCoupledSolver(params, None, "for_testing")

        with self.assertRaisesRegex(Exception, 'No solver_name was specified!'):
            coupled_solver.model[""]

        # checking Models
        fluid_model = coupled_solver.model["fluid"]
        structural_model = coupled_solver.model["structure"]

        self.assertIsInstance(fluid_model, KM.Model)
        self.assertIsInstance(structural_model, KM.Model)

        self.assertIsNot(structural_model, fluid_model)
        self.assertIs(fluid_model, coupled_solver._GetSolver("fluid").model)
        self.assertIs(structural_model, coupled_solver._GetSolver("structure").model)

        # checking ModelParts
        fluid_main_mp = coupled_solver.model["fluid.fl_model_part"]
        self.assertIsInstance(fluid_main_mp, KM.ModelPart)
        self.assertEqual(fluid_main_mp.Name, "fl_model_part")

        # creating some submodelparts to see if getting those also works
        smp_1 = fluid_main_mp.CreateSubModelPart("boundary")
        smp_11 = smp_1.CreateSubModelPart("slip")

        boundary_mp = coupled_solver.model["fluid.fl_model_part.boundary"]
        self.assertIsInstance(boundary_mp, KM.ModelPart)
        self.assertEqual(boundary_mp.Name, "boundary")
        self.assertEqual(boundary_mp.FullName(), "fl_model_part.boundary")

        slip_mp = coupled_solver.model["fluid.fl_model_part.boundary.slip"]
        self.assertIsInstance(slip_mp, KM.ModelPart)
        self.assertEqual(slip_mp.Name, "slip")
        self.assertEqual(slip_mp.FullName(), "fl_model_part.boundary.slip")

    def test_model_access_two_levels(self):
        params = KM.Parameters("""{
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
                                "main_model_part_name" : "fl_model_part"
                            }
                        },
                        "structure": {
                            "type"        : "helpers.dummy_solver_wrapper",
                            "solver_wrapper_settings": {
                                "time_step" : 1.0,
                                "domain_size" : 1,
                                "main_model_part_name" : "strct_model_part"
                            }
                        }
                    }
                },
                "controller": {
                    "type"        : "helpers.dummy_solver_wrapper",
                    "solver_wrapper_settings": {
                        "time_step" : 1.0,
                        "domain_size" : 1,
                        "main_model_part_name" : "interface"
                    }
                }
            }
        }""")

        coupled_solver = CoSimulationCoupledSolver(params, None, "for_testing")

        fluid_main_mp = coupled_solver.model["fsi.fluid.fl_model_part"]
        self.assertIsInstance(fluid_main_mp, KM.ModelPart)
        self.assertEqual(fluid_main_mp.Name, "fl_model_part")

        # creating some submodelparts to see if getting those also works
        smp_1 = fluid_main_mp.CreateSubModelPart("boundary")
        smp_1.CreateSubModelPart("slip")

        boundary_mp = coupled_solver.model["fsi.fluid.fl_model_part.boundary"]
        self.assertIsInstance(boundary_mp, KM.ModelPart)
        self.assertEqual(boundary_mp.Name, "boundary")
        self.assertEqual(boundary_mp.FullName(), "fl_model_part.boundary")

        slip_mp = coupled_solver.model["fsi.fluid.fl_model_part.boundary.slip"]
        self.assertIsInstance(slip_mp, KM.ModelPart)
        self.assertEqual(slip_mp.Name, "slip")
        self.assertEqual(slip_mp.FullName(), "fl_model_part.boundary.slip")


class TestCoupledSolverCouplingInterfaceDataAccess(KratosUnittest.TestCase):
    def test_coupling_interface_data_access_one_level(self):
        pass

    def test_coupling_interface_data_access_two_levels(self):
        pass


class TestCoupledSolverPassingModel(KratosUnittest.TestCase):
    def test_pass_model_to_solver_wrapper(self):
        # CoSim contains only one SolverWrapper (not really CoSimulation...)
        params = KM.Parameters("""{
            "solver_wrapper_settings": {
                "time_step" : 1.0,
                "domain_size" : 1,
                "main_model_part_name" : "dummy"
            }
        }""")

        solver_without_outside_model = DummySolverWrapper(params, None, "dummy_2")
        self.assertIsInstance(solver_without_outside_model.model, KM.Model)

        outside_model = KM.Model()
        solver = DummySolverWrapper(params, outside_model, "dummy")
        self.assertIs(outside_model, solver.model)

        with self.assertRaisesRegex(Exception, 'A solver wrapper can either be passed a Model\nor None, got object of type'):
            DummySolverWrapper(params, "wrongtype", "dummy")



    def test_pass_model_zero_levels(self):
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

        outside_model = KM.Model()

        co_sim_analysis = CoSimulationAnalysis(params, outside_model)

        top_level_solver_wrapper = co_sim_analysis._GetSolver()

        self.assertIs(outside_model, top_level_solver_wrapper.model)

    def test_not_passing_model_zero_levels(self):
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

        self.assertIsInstance(top_level_solver_wrapper.model, KM.Model)


    def test_pass_model_one_level(self):
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
                            "main_model_part_name" : "fl_model_part"
                        }
                    },
                    "structure": {
                        "type"        : "helpers.dummy_solver_wrapper",
                        "solver_wrapper_settings": {
                            "time_step" : 1.0,
                            "domain_size" : 1,
                            "main_model_part_name" : "strct_model_part"
                        }
                    }
                }
            }
        }""")

        structure_model = KM.Model()

        with self.assertRaisesRegex(Exception, 'A coupled solver can either be passed a dict of Models\nor None, got object of type'):
            co_sim_analysis = CoSimulationAnalysis(params, structure_model)

        with self.assertRaisesRegex(Exception, 'A solver wrapper can either be passed a Model\nor None, got object of type'): # this throws in the SolverWrapper
            co_sim_analysis = CoSimulationAnalysis(params, {"structure": "structure_model"})

        with self.assertRaisesRegex(Exception, 'A Model was given for solver "wrong_name" but this solver does not exist!'):
            co_sim_analysis = CoSimulationAnalysis(params, {"wrong_name": structure_model})

        co_sim_analysis = CoSimulationAnalysis(params, {"structure": structure_model}) # name has to match the one in json!

        top_level_solver_wrapper = co_sim_analysis._GetSolver()
        fluid_solver_wrapper = co_sim_analysis._GetSolver("fluid")
        structural_solver_wrapper = co_sim_analysis._GetSolver("structure")

        self.assertIs(structure_model, top_level_solver_wrapper.model["structure"])
        self.assertIsNot(structure_model, top_level_solver_wrapper.model["fluid"])
        self.assertIs(structure_model, structural_solver_wrapper.model)
        self.assertIsNot(structure_model, fluid_solver_wrapper.model)


    def test_pass_model_two_levels(self):
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
                                    "main_model_part_name" : "fl_model_part"
                                }
                            },
                            "structure": {
                                "type"        : "helpers.dummy_solver_wrapper",
                                "solver_wrapper_settings": {
                                    "time_step" : 1.0,
                                    "domain_size" : 1,
                                    "main_model_part_name" : "strct_model_part"
                                }
                            }
                        }
                    },
                    "controller": {
                        "type"        : "helpers.dummy_solver_wrapper",
                        "solver_wrapper_settings": {
                            "time_step" : 1.0,
                            "domain_size" : 1,
                            "main_model_part_name" : "interface"
                        }
                    }
                }
            }
        }""")

        structure_model = KM.Model()
        controller_model = KM.Model()

        models = {
            "controller" : controller_model,
            "fsi" : {
                "structure" : structure_model
            }
        }

        co_sim_analysis = CoSimulationAnalysis(params, models)

        top_level_solver_wrapper = co_sim_analysis._GetSolver()
        fluid_solver_wrapper = co_sim_analysis._GetSolver("fsi.fluid")
        structural_solver_wrapper = co_sim_analysis._GetSolver("fsi.structure")
        controller_solver_wrapper = co_sim_analysis._GetSolver("controller")

        self.assertIs(structure_model, structural_solver_wrapper.model)
        self.assertIs(structure_model, top_level_solver_wrapper.model["fsi.structure"])

        self.assertIs(controller_model, controller_solver_wrapper.model)
        self.assertIs(controller_model, top_level_solver_wrapper.model["controller"])

        self.assertIsNot(structure_model, fluid_solver_wrapper.model)
        self.assertIsNot(controller_model, fluid_solver_wrapper.model)

        self.assertIs(fluid_solver_wrapper.model, top_level_solver_wrapper.model["fsi.fluid"])


if __name__ == '__main__':
    KratosUnittest.main()

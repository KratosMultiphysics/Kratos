from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

import KratosMultiphysics as KM
import KratosMultiphysics.CoSimulationApplication as KMC
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utils

from KratosMultiphysics.CoSimulationApplication.factories import data_transfer_operator_factory
from KratosMultiphysics.CoSimulationApplication.coupling_interface_data import CouplingInterfaceData

from KratosMultiphysics.CoSimulationApplication.co_simulation_tools import UsingPyKratos
using_pykratos = UsingPyKratos()

mapping_app_available = kratos_utils.CheckIfApplicationsAvailable("MappingApplication")


class TestDataTransferOperators(KratosUnittest.TestCase):

    def setUp(self):
        self.model = KM.Model()
        mp_o = self.model.CreateModelPart("mp_origin", 2)
        mp_d_m = self.model.CreateModelPart("mp_destination_matching", 2)
        mp_d_nm = self.model.CreateModelPart("mp_destination_non_matching", 2)
        mp_one_n = self.model.CreateModelPart("mp_single_node", 2)

        mp_o.AddNodalSolutionStepVariable(KM.PRESSURE)
        mp_o.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
        mp_d_m.AddNodalSolutionStepVariable(KM.TEMPERATURE)
        mp_d_m.AddNodalSolutionStepVariable(KM.FORCE)
        mp_d_nm.AddNodalSolutionStepVariable(KM.TEMPERATURE)
        mp_d_nm.AddNodalSolutionStepVariable(KM.FORCE)
        mp_one_n.AddNodalSolutionStepVariable(KMC.SCALAR_DISPLACEMENT)
        mp_one_n.AddNodalSolutionStepVariable(KMC.SCALAR_FORCE)

        mp_o.ProcessInfo[KM.DOMAIN_SIZE] = 2
        mp_d_m.ProcessInfo[KM.DOMAIN_SIZE] = 2
        mp_d_nm.ProcessInfo[KM.DOMAIN_SIZE] = 2
        mp_one_n.ProcessInfo[KM.DOMAIN_SIZE] = 2

        num_nodes_matching = 5
        num_nodes_non_matching = 8

        for i in range(num_nodes_matching):
            node_id = i+1
            node_o = mp_o.CreateNewNode(node_id, 0.0, 0.0, i+1)
            mp_d_m.CreateNewNode(node_id, 0.0, 0.0, i+1)

            node_o.SetSolutionStepValue(KM.PRESSURE, 0, ScalarValueFromId(node_id))
            node_o.SetSolutionStepValue(KM.DISPLACEMENT, 0, VectorValueFromId(node_id))

        for i in range(num_nodes_non_matching-1,-1, -1):
            node_id = i+15
            mp_d_nm.CreateNewNode(node_id, 0.0, 0.0, i+1.1)

        mp_one_n.CreateNewNode(1, 0.0, 0.0, 0.0)

        origin_data_settings_scalar = KM.Parameters("""{
            "model_part_name" : "mp_origin",
            "variable_name"   : "PRESSURE"
        }""")
        origin_data_settings_vector = KM.Parameters("""{
            "model_part_name" : "mp_origin",
            "variable_name"   : "DISPLACEMENT",
            "dimension" : 2
        }""")
        origin_data_settings_single_node = KM.Parameters("""{
            "model_part_name" : "mp_single_node",
            "variable_name"   : "SCALAR_DISPLACEMENT"
        }""")

        self.origin_data_scalar = CouplingInterfaceData(origin_data_settings_scalar, self.model)
        self.origin_data_vector = CouplingInterfaceData(origin_data_settings_vector, self.model)
        self.origin_data_single_node = CouplingInterfaceData(origin_data_settings_single_node, self.model)
        self.origin_data_scalar.Initialize()
        self.origin_data_vector.Initialize()
        self.origin_data_single_node.Initialize()


        destination_matching_data_settings_scalar = KM.Parameters("""{
            "model_part_name" : "mp_destination_matching",
            "variable_name"   : "TEMPERATURE"
        }""")
        destination_matching_data_settings_vector = KM.Parameters("""{
            "model_part_name" : "mp_destination_matching",
            "variable_name"   : "FORCE",
            "dimension" : 2
        }""")
        destination_data_settings_single_node = KM.Parameters("""{
            "model_part_name" : "mp_single_node",
            "variable_name"   : "SCALAR_FORCE"
        }""")

        self.destination_matching_data_scalar = CouplingInterfaceData(destination_matching_data_settings_scalar, self.model)
        self.destination_matching_data_vector = CouplingInterfaceData(destination_matching_data_settings_vector, self.model)
        self.destination_data_single_node = CouplingInterfaceData(destination_data_settings_single_node, self.model)
        self.destination_matching_data_scalar.Initialize()
        self.destination_matching_data_vector.Initialize()
        self.destination_data_single_node.Initialize()


        destination_non_matching_data_settings_scalar = KM.Parameters("""{
            "model_part_name" : "mp_destination_non_matching",
            "variable_name"   : "TEMPERATURE"
        }""")
        destination_non_matching_data_settings_vector = KM.Parameters("""{
            "model_part_name" : "mp_destination_non_matching",
            "variable_name"   : "FORCE",
            "dimension" : 2
        }""")

        self.destination_non_matching_data_scalar = CouplingInterfaceData(destination_non_matching_data_settings_scalar, self.model)
        self.destination_non_matching_data_vector = CouplingInterfaceData(destination_non_matching_data_settings_vector, self.model)
        self.destination_non_matching_data_scalar.Initialize()
        self.destination_non_matching_data_vector.Initialize()


    def test_copy_transfer_operator(self):
        data_transfer_op_settings = KM.Parameters("""{
            "type" : "copy"
        }""")

        data_transfer_op = data_transfer_operator_factory.CreateDataTransferOperator(data_transfer_op_settings)

        self.__TestTransferMatching(data_transfer_op)
        self.__TestTransferMatchingSwapSign(data_transfer_op)
        self.__TestTransferMatchingAddValues(data_transfer_op)
        self.__TestTransferMatchingAddValuesAndSwapSign(data_transfer_op)

        transfer_options_empty = KM.Parameters(""" [] """)
        exp_error = 'The sizes of the data are not matching: {} != {} for interface data "{}" of solver "{}" and interface data "{}" of solver "{}"!'.format(self.origin_data_scalar.Size(), self.destination_non_matching_data_scalar.Size(), self.origin_data_scalar.name, self.origin_data_scalar.solver_name, self.destination_non_matching_data_scalar.name, self.destination_non_matching_data_scalar.solver_name)
        with self.assertRaisesRegex(Exception, exp_error):
            data_transfer_op.TransferData(self.origin_data_scalar, self.destination_non_matching_data_scalar, transfer_options_empty)

    def test_kratos_mapping_transfer_operator(self):
        if using_pykratos:
            self.skipTest("This test cannot be run with pyKratos!")
        if not mapping_app_available:
            self.skipTest("MappingApplication not available!")

        data_transfer_op_settings_missing = KM.Parameters("""{
            "type" : "kratos_mapping"
        }""")

        exp_error = 'No "mapper_settings" provided!'
        with self.assertRaisesRegex(Exception, exp_error):
            data_transfer_operator_factory.CreateDataTransferOperator(data_transfer_op_settings_missing)

        data_transfer_op_settings = KM.Parameters("""{
            "type" : "kratos_mapping",
            "mapper_settings" : {
                "mapper_type" : "nearest_neighbor"
            }
        }""")

        data_transfer_op = data_transfer_operator_factory.CreateDataTransferOperator(data_transfer_op_settings)

        self.__TestTransferMatching(data_transfer_op)
        self.__TestTransferMatchingSwapSign(data_transfer_op)
        self.__TestTransferMatchingAddValues(data_transfer_op)
        self.__TestTransferMatchingAddValuesAndSwapSign(data_transfer_op)

        # with this we make sure that only one mapper is created (and not several ones for each mapping operation!)
        # Hint: requires access to private member
        self.assertEqual(len(data_transfer_op._KratosMappingDataTransferOperator__mappers), 1)

        self.__TestTransferMatchingInverse(data_transfer_op)
        # here we check explicitly the InverseMap fct
        self.assertEqual(len(data_transfer_op._KratosMappingDataTransferOperator__mappers), 1)

        transfer_options_empty = KM.Parameters(""" [] """)
        data_transfer_op.TransferData(self.origin_data_scalar, self.destination_non_matching_data_scalar, transfer_options_empty)
        # here we check explicitly the creation of a second mapper, which is required since the interfaces are not the same this time
        self.assertEqual(len(data_transfer_op._KratosMappingDataTransferOperator__mappers), 2)

        data_settings_model_part = KM.Parameters("""{
            "model_part_name" : "mp_single_node",
            "variable_name"   : "TEMPERATURE",
            "location"        : "model_part"
        }""")

        data_model_part = CouplingInterfaceData(data_settings_model_part, self.model)
        data_model_part.Initialize()

        with self.assertRaisesRegex(Exception, 'Currently only historical nodal values are supported'):
            data_transfer_op.TransferData(self.origin_data_scalar, data_model_part, transfer_options_empty)

    def test_copy_single_to_dist_transfer_operator(self):
        data_transfer_op_settings = KM.Parameters("""{
            "type" : "copy_single_to_distributed"
        }""")

        data_transfer_op = data_transfer_operator_factory.CreateDataTransferOperator(data_transfer_op_settings)
        transfer_options = KM.Parameters(""" [] """)

        for node in self.origin_data_single_node.GetModelPart().Nodes:
            node.SetSolutionStepValue(KMC.SCALAR_DISPLACEMENT, 0, 100.0)

        data_transfer_op.TransferData(self.origin_data_single_node,
                                      self.destination_matching_data_scalar,
                                      transfer_options)

        self.__CompareScalarNodalValues(self.destination_matching_data_scalar.GetModelPart().Nodes,
                                        self.origin_data_single_node.GetModelPart().Nodes,
                                        KM.TEMPERATURE, KMC.SCALAR_DISPLACEMENT,1)

        with self.assertRaisesRegex(Exception, 'Interface data "default" of solver "default_solver" requires to be of size 1, got: 5'):
            data_transfer_op.TransferData(self.origin_data_scalar,
                                          self.destination_matching_data_scalar,
                                          transfer_options)


    def test_copy_single_to_dist_transfer_operator_swap_sign(self):
        data_transfer_op_settings = KM.Parameters("""{
            "type" : "copy_single_to_distributed"
        }""")

        data_transfer_op = data_transfer_operator_factory.CreateDataTransferOperator(data_transfer_op_settings)
        transfer_options = KM.Parameters(""" ["swap_sign"] """)

        for node in self.origin_data_single_node.GetModelPart().Nodes:
            node.SetSolutionStepValue(KMC.SCALAR_DISPLACEMENT, 0, 100.0)

        data_transfer_op.TransferData(self.origin_data_single_node,
                                      self.destination_matching_data_scalar,
                                      transfer_options)

        self.__CompareScalarNodalValues(self.destination_matching_data_scalar.GetModelPart().Nodes,
                                        self.origin_data_single_node.GetModelPart().Nodes,
                                        KM.TEMPERATURE, KMC.SCALAR_DISPLACEMENT,-1)

    def test_copy_single_to_dist_transfer_operator_distribute_values(self):
        data_transfer_op_settings = KM.Parameters("""{
            "type" : "copy_single_to_distributed"
        }""")

        data_transfer_op = data_transfer_operator_factory.CreateDataTransferOperator(data_transfer_op_settings)
        transfer_options = KM.Parameters(""" ["distribute_values"] """)

        for node in self.origin_data_single_node.GetModelPart().Nodes:
            node.SetSolutionStepValue(KMC.SCALAR_DISPLACEMENT, 0, 100.0)

        data_transfer_op.TransferData(self.origin_data_single_node,
                                      self.destination_matching_data_scalar,
                                      transfer_options)

        self.__CompareScalarNodalValues(self.destination_matching_data_scalar.GetModelPart().Nodes,
                                        self.origin_data_single_node.GetModelPart().Nodes,
                                        KM.TEMPERATURE, KMC.SCALAR_DISPLACEMENT, 0.2)

    def test_copy_single_to_dist_transfer_operator_add_values(self):
        data_transfer_op_settings = KM.Parameters("""{
            "type" : "copy_single_to_distributed"
        }""")

        data_transfer_op = data_transfer_operator_factory.CreateDataTransferOperator(data_transfer_op_settings)
        transfer_options_empty = KM.Parameters(""" [] """)
        transfer_options_add = KM.Parameters(""" ["add_values"] """)

        for node in self.origin_data_single_node.GetModelPart().Nodes:
            node.SetSolutionStepValue(KMC.SCALAR_DISPLACEMENT, 0, 100.0)

        data_transfer_op.TransferData(self.origin_data_single_node,
                                      self.destination_matching_data_scalar,
                                      transfer_options_empty)

        self.__CompareScalarNodalValues(self.destination_matching_data_scalar.GetModelPart().Nodes,
                                        self.origin_data_single_node.GetModelPart().Nodes,
                                        KM.TEMPERATURE, KMC.SCALAR_DISPLACEMENT, 1)

        data_transfer_op.TransferData(self.origin_data_single_node,
                                      self.destination_matching_data_scalar,
                                      transfer_options_add)

        self.__CompareScalarNodalValues(self.destination_matching_data_scalar.GetModelPart().Nodes,
                                        self.origin_data_single_node.GetModelPart().Nodes,
                                        KM.TEMPERATURE, KMC.SCALAR_DISPLACEMENT, 2)

    def test_copy_single_to_dist_transfer_operator_distribute_values_swap_sign(self):
        data_transfer_op_settings = KM.Parameters("""{
            "type" : "copy_single_to_distributed"
        }""")

        data_transfer_op = data_transfer_operator_factory.CreateDataTransferOperator(data_transfer_op_settings)
        transfer_options = KM.Parameters(""" ["distribute_values", "swap_sign"] """)

        for node in self.origin_data_single_node.GetModelPart().Nodes:
            node.SetSolutionStepValue(KMC.SCALAR_DISPLACEMENT, 0, 100.0)

        data_transfer_op.TransferData(self.origin_data_single_node,
                                      self.destination_matching_data_scalar,
                                      transfer_options)

        self.__CompareScalarNodalValues(self.destination_matching_data_scalar.GetModelPart().Nodes,
                                        self.origin_data_single_node.GetModelPart().Nodes,
                                        KM.TEMPERATURE, KMC.SCALAR_DISPLACEMENT, -0.2)

    def test_copy_single_to_dist_transfer_operator_distribute_values_swap_sign_add_values(self):
        data_transfer_op_settings = KM.Parameters("""{
            "type" : "copy_single_to_distributed"
        }""")

        data_transfer_op = data_transfer_operator_factory.CreateDataTransferOperator(data_transfer_op_settings)
        transfer_options = KM.Parameters(""" ["distribute_values", "swap_sign"] """)
        transfer_options_with_add_vals = KM.Parameters(""" ["distribute_values", "swap_sign", "add_values"] """)

        for node in self.origin_data_single_node.GetModelPart().Nodes:
            node.SetSolutionStepValue(KMC.SCALAR_DISPLACEMENT, 0, 100.0)

        data_transfer_op.TransferData(self.origin_data_single_node,
                                      self.destination_matching_data_scalar,
                                      transfer_options)

        self.__CompareScalarNodalValues(self.destination_matching_data_scalar.GetModelPart().Nodes,
                                        self.origin_data_single_node.GetModelPart().Nodes,
                                        KM.TEMPERATURE, KMC.SCALAR_DISPLACEMENT, -0.2)

        data_transfer_op.TransferData(self.origin_data_single_node,
                                      self.destination_matching_data_scalar,
                                      transfer_options_with_add_vals)

        self.__CompareScalarNodalValues(self.destination_matching_data_scalar.GetModelPart().Nodes,
                                        self.origin_data_single_node.GetModelPart().Nodes,
                                        KM.TEMPERATURE, KMC.SCALAR_DISPLACEMENT, -0.4)

    def test_sum_dist_to_single(self):
        data_transfer_op_settings = KM.Parameters("""{
            "type" : "sum_distributed_to_single"
        }""")

        data_transfer_op = data_transfer_operator_factory.CreateDataTransferOperator(data_transfer_op_settings)
        transfer_options = KM.Parameters(""" [] """)

        for node in self.origin_data_scalar.GetModelPart().Nodes:
            node.SetSolutionStepValue(KM.PRESSURE, 0, 100.0)

        data_transfer_op.TransferData(self.origin_data_scalar,
                                      self.destination_data_single_node,
                                      transfer_options)

        self.__CompareScalarNodalValues(self.destination_data_single_node.GetModelPart().Nodes,
                                        self.origin_data_scalar.GetModelPart().Nodes,
                                        KMC.SCALAR_FORCE, KM.PRESSURE, 5)

        with self.assertRaisesRegex(Exception, 'Interface data "default" of solver "default_solver" requires to be of size 1, got: 5'):
            data_transfer_op.TransferData(self.origin_data_scalar,
                                        self.destination_matching_data_scalar,
                                        transfer_options)


    def test_sum_dist_to_single_swap_sign(self):
        data_transfer_op_settings = KM.Parameters("""{
            "type" : "sum_distributed_to_single"
        }""")

        data_transfer_op = data_transfer_operator_factory.CreateDataTransferOperator(data_transfer_op_settings)
        transfer_options = KM.Parameters(""" ["swap_sign"] """)

        for node in self.origin_data_scalar.GetModelPart().Nodes:
            node.SetSolutionStepValue(KM.PRESSURE, 0, 100.0)

        data_transfer_op.TransferData(self.origin_data_scalar,
                                      self.destination_data_single_node,
                                      transfer_options)

        self.__CompareScalarNodalValues(self.destination_data_single_node.GetModelPart().Nodes,
                                        self.origin_data_scalar.GetModelPart().Nodes,
                                        KMC.SCALAR_FORCE, KM.PRESSURE, -5)


    def test_sum_dist_to_single_add_values(self):
        data_transfer_op_settings = KM.Parameters("""{
            "type" : "sum_distributed_to_single"
        }""")

        data_transfer_op = data_transfer_operator_factory.CreateDataTransferOperator(data_transfer_op_settings)
        transfer_options = KM.Parameters(""" ["add_values"] """)

        for node in self.origin_data_scalar.GetModelPart().Nodes:
            node.SetSolutionStepValue(KM.PRESSURE, 0, 100.0)

        for node in self.destination_data_single_node.GetModelPart().Nodes:
            node.SetSolutionStepValue(KMC.SCALAR_FORCE, 0, 500.0)

        data_transfer_op.TransferData(self.origin_data_scalar,
                                      self.destination_data_single_node,
                                      transfer_options)

        self.__CompareScalarNodalValues(self.destination_data_single_node.GetModelPart().Nodes,
                                        self.origin_data_scalar.GetModelPart().Nodes,
                                        KMC.SCALAR_FORCE, KM.PRESSURE, 10.0)

    def test_sum_dist_to_single_add_values_swap_sign(self):
        data_transfer_op_settings = KM.Parameters("""{
            "type" : "sum_distributed_to_single"
        }""")

        data_transfer_op = data_transfer_operator_factory.CreateDataTransferOperator(data_transfer_op_settings)
        transfer_options = KM.Parameters(""" ["add_values", "swap_sign"] """)

        for node in self.origin_data_scalar.GetModelPart().Nodes:
            node.SetSolutionStepValue(KM.PRESSURE, 0, 100.0)

        for node in self.destination_data_single_node.GetModelPart().Nodes:
            node.SetSolutionStepValue(KMC.SCALAR_FORCE, 0, 500.0)

        data_transfer_op.TransferData(self.origin_data_scalar,
                                      self.destination_data_single_node,
                                      transfer_options)

        self.__CompareScalarNodalValues(self.destination_data_single_node.GetModelPart().Nodes,
                                        self.origin_data_scalar.GetModelPart().Nodes,
                                        KMC.SCALAR_FORCE, KM.PRESSURE, 0.0)

    def test_sum_dist_to_single_check_var(self):
        data_transfer_op_settings = KM.Parameters("""{
            "type" : "sum_distributed_to_single"
        }""")

        data_transfer_op = data_transfer_operator_factory.CreateDataTransferOperator(data_transfer_op_settings)
        transfer_options = KM.Parameters(""" [] """)

        with self.assertRaisesRegex(Exception, 'Variable of interface data "default" of solver "default_solver" has to be a scalar!'):
            data_transfer_op.TransferData(self.destination_matching_data_vector,
                                        self.destination_data_single_node,
                                        transfer_options)

    def test_copy_single_to_dist_check_var(self):
        data_transfer_op_settings = KM.Parameters("""{
            "type" : "copy_single_to_distributed"
        }""")

        data_transfer_op = data_transfer_operator_factory.CreateDataTransferOperator(data_transfer_op_settings)
        transfer_options = KM.Parameters(""" [] """)

        with self.assertRaisesRegex(Exception, 'Variable of interface data "default" of solver "default_solver" has to be a scalar!'):
            data_transfer_op.TransferData(self.destination_data_single_node,
                                        self.destination_matching_data_vector,
                                        transfer_options)

    def __TestTransferMatching(self, data_transfer_op):
        transfer_options_empty = KM.Parameters(""" [] """)
        data_transfer_op.TransferData(self.origin_data_scalar, self.destination_matching_data_scalar, transfer_options_empty)
        self.__CompareScalarNodalValues(self.destination_matching_data_scalar.GetModelPart().Nodes, self.origin_data_scalar.GetModelPart().Nodes, KM.TEMPERATURE, KM.PRESSURE)
        data_transfer_op.TransferData(self.origin_data_vector, self.destination_matching_data_vector, transfer_options_empty)
        self.__CompareVectorNodalValues(self.destination_matching_data_vector.GetModelPart().Nodes, self.origin_data_vector.GetModelPart().Nodes, KM.FORCE, KM.DISPLACEMENT, 2)

    def __TestTransferMatchingInverse(self, data_transfer_op):
        transfer_options_empty = KM.Parameters(""" [] """)
        data_transfer_op.TransferData(self.destination_matching_data_scalar, self.origin_data_scalar, transfer_options_empty)
        self.__CompareScalarNodalValues(self.destination_matching_data_scalar.GetModelPart().Nodes, self.origin_data_scalar.GetModelPart().Nodes, KM.TEMPERATURE, KM.PRESSURE)
        data_transfer_op.TransferData(self.destination_matching_data_vector, self.origin_data_vector, transfer_options_empty)
        self.__CompareVectorNodalValues(self.destination_matching_data_vector.GetModelPart().Nodes, self.origin_data_vector.GetModelPart().Nodes, KM.FORCE, KM.DISPLACEMENT, 2)

        transfer_options_fail = KM.Parameters(""" ["thisWillHopefullyNeverBeImplementedOtherWiseThisTestWillFail"] """)
        with self.assertRaisesRegex(Exception, ' not recognized for '):
            data_transfer_op.TransferData(self.origin_data_scalar, self.destination_matching_data_scalar, transfer_options_fail)

    def __TestTransferMatchingSwapSign(self, data_transfer_op):
        transfer_options_swap_sign = KM.Parameters(""" ["swap_sign"] """)
        data_transfer_op.TransferData(self.origin_data_scalar, self.destination_matching_data_scalar, transfer_options_swap_sign)
        self.__CompareScalarNodalValues(self.destination_matching_data_scalar.GetModelPart().Nodes, self.origin_data_scalar.GetModelPart().Nodes, KM.TEMPERATURE, KM.PRESSURE, -1.0)
        data_transfer_op.TransferData(self.origin_data_vector, self.destination_matching_data_vector, transfer_options_swap_sign)
        self.__CompareVectorNodalValues(self.destination_matching_data_vector.GetModelPart().Nodes, self.origin_data_vector.GetModelPart().Nodes, KM.FORCE, KM.DISPLACEMENT, 2, -1.0)

    def __TestTransferMatchingAddValues(self, data_transfer_op):
        transfer_options_empty = KM.Parameters(""" [] """)
        transfer_options_add_values = KM.Parameters(""" ["add_values"] """)
        data_transfer_op.TransferData(self.origin_data_scalar, self.destination_matching_data_scalar, transfer_options_empty) # "resetting"
        data_transfer_op.TransferData(self.origin_data_scalar, self.destination_matching_data_scalar, transfer_options_add_values)
        self.__CompareScalarNodalValues(self.destination_matching_data_scalar.GetModelPart().Nodes, self.origin_data_scalar.GetModelPart().Nodes, KM.TEMPERATURE, KM.PRESSURE, 2.0)
        data_transfer_op.TransferData(self.origin_data_vector, self.destination_matching_data_vector, transfer_options_empty) # "resetting"
        data_transfer_op.TransferData(self.origin_data_vector, self.destination_matching_data_vector, transfer_options_add_values)
        self.__CompareVectorNodalValues(self.destination_matching_data_vector.GetModelPart().Nodes, self.origin_data_vector.GetModelPart().Nodes, KM.FORCE, KM.DISPLACEMENT, 2, 2.0)

    def __TestTransferMatchingAddValuesAndSwapSign(self, data_transfer_op):
        transfer_options_empty = KM.Parameters(""" [] """)
        transfer_options_add_values = KM.Parameters(""" ["add_values"] """)
        transfer_options_add_values_swap_sign = KM.Parameters(""" ["add_values", "swap_sign"] """)
        data_transfer_op.TransferData(self.origin_data_scalar, self.destination_matching_data_scalar, transfer_options_empty) # "resetting"
        data_transfer_op.TransferData(self.origin_data_scalar, self.destination_matching_data_scalar, transfer_options_add_values)
        data_transfer_op.TransferData(self.origin_data_scalar, self.destination_matching_data_scalar, transfer_options_add_values)
        data_transfer_op.TransferData(self.origin_data_scalar, self.destination_matching_data_scalar, transfer_options_add_values_swap_sign)
        self.__CompareScalarNodalValues(self.destination_matching_data_scalar.GetModelPart().Nodes, self.origin_data_scalar.GetModelPart().Nodes, KM.TEMPERATURE, KM.PRESSURE, 2.0)
        data_transfer_op.TransferData(self.origin_data_vector, self.destination_matching_data_vector, transfer_options_empty) # "resetting"
        data_transfer_op.TransferData(self.origin_data_vector, self.destination_matching_data_vector, transfer_options_add_values)
        data_transfer_op.TransferData(self.origin_data_vector, self.destination_matching_data_vector, transfer_options_add_values)
        data_transfer_op.TransferData(self.origin_data_vector, self.destination_matching_data_vector, transfer_options_add_values_swap_sign)
        self.__CompareVectorNodalValues(self.destination_matching_data_vector.GetModelPart().Nodes, self.origin_data_vector.GetModelPart().Nodes, KM.FORCE, KM.DISPLACEMENT, 2, 2.0)


    def __CompareScalarNodalValues(self, nodes, nodes_ref, var, var_ref, factor=1.0):
        for node, node_ref in zip(nodes, nodes_ref):
            self.assertAlmostEqual(node.GetSolutionStepValue(var), node_ref.GetSolutionStepValue(var_ref)*factor, 10)

    def __CompareVectorNodalValues(self, nodes, nodes_ref, var, var_ref, dimension, factor=1.0):
        for node, node_ref in zip(nodes, nodes_ref):
            val = node.GetSolutionStepValue(var)
            val_ref = node_ref.GetSolutionStepValue(var_ref)
            for i in range(dimension):
                self.assertAlmostEqual(val[i], val_ref[i]*factor, 10)


def ScalarValueFromId(the_id):
    return the_id*1.5

def VectorValueFromId(the_id):
    return [the_id*1.1, the_id*1.1+2.5, the_id*2.3]

if __name__ == '__main__':
    KratosUnittest.main()

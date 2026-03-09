import KratosMultiphysics as KM
import KratosMultiphysics.CoSimulationApplication as KMC
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utils

from KratosMultiphysics.CoSimulationApplication.factories import data_transfer_operator_factory
from KratosMultiphysics.CoSimulationApplication.coupling_interface_data import CouplingInterfaceData

mapping_app_available = kratos_utils.CheckIfApplicationsAvailable("MappingApplication")

if KM.IsDistributedRun():
    import KratosMultiphysics.mpi as KratosMPI

class TestDataTransferOperators(KratosUnittest.TestCase):

    def setUp(self):
        self.model = KM.Model()
        mp_o = self.model.CreateModelPart("mp_origin", 2)
        mp_d_m = self.model.CreateModelPart("mp_destination_matching", 2)
        mp_d_nm = self.model.CreateModelPart("mp_destination_non_matching", 2)
        mp_one_n = self.model.CreateModelPart("mp_single_node", 2)

        self.data_comm = KM.Testing.GetDefaultDataCommunicator()
        self.my_pid = self.data_comm.Rank()

        mp_o.AddNodalSolutionStepVariable(KM.PRESSURE)
        mp_o.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
        mp_o.AddNodalSolutionStepVariable(KM.PARTITION_INDEX)

        mp_d_m.AddNodalSolutionStepVariable(KM.TEMPERATURE)
        mp_d_m.AddNodalSolutionStepVariable(KM.FORCE)
        mp_d_m.AddNodalSolutionStepVariable(KM.PARTITION_INDEX)

        mp_d_nm.AddNodalSolutionStepVariable(KM.TEMPERATURE)
        mp_d_nm.AddNodalSolutionStepVariable(KM.FORCE)
        mp_d_nm.AddNodalSolutionStepVariable(KM.PARTITION_INDEX)

        mp_one_n.AddNodalSolutionStepVariable(KMC.SCALAR_DISPLACEMENT)
        mp_one_n.AddNodalSolutionStepVariable(KMC.SCALAR_FORCE)
        mp_one_n.AddNodalSolutionStepVariable(KM.PARTITION_INDEX)

        mp_o.ProcessInfo[KM.DOMAIN_SIZE] = 2
        mp_d_m.ProcessInfo[KM.DOMAIN_SIZE] = 2
        mp_d_nm.ProcessInfo[KM.DOMAIN_SIZE] = 2
        mp_one_n.ProcessInfo[KM.DOMAIN_SIZE] = 2

        num_nodes_matching = self.my_pid % 5 + 3 # num_nodes in range (3 ... 7)
        num_nodes_non_matching = self.my_pid % 5 + 5 # num_nodes in range (3 ... 7)

        if self.my_pid == 3:
            num_nodes_matching = 0 # in order to emulate one partition not having local nodes
        if self.my_pid == 2:
            num_nodes_non_matching = 0 # in order to emulate one partition not having local nodes

        for i in range(num_nodes_matching):
            node_id = i+1
            node_o = mp_o.CreateNewNode(node_id, 0.0, 0.0, i+1)
            mp_d_m.CreateNewNode(node_id, 0.0, 0.0, i+1)

            node_o.SetSolutionStepValue(KM.PRESSURE, 0, ScalarValueFromId(node_id))
            node_o.SetSolutionStepValue(KM.DISPLACEMENT, 0, VectorValueFromId(node_id))
            node_o.SetValue(KM.LAMBDA, ScalarValueFromId(node_id)*1.47)

        KM.VariableUtils().SetNonHistoricalVariableToZero(KM.VISCOSITY, mp_d_m.Nodes)

        for i in range(num_nodes_non_matching-1,-1, -1):
            node_id = i+15
            mp_d_nm.CreateNewNode(node_id, 0.0, 0.0, i+1.1)

        if self.my_pid == 0:
            mp_one_n.CreateNewNode(1, 0.0, 0.0, 0.0)

        for node in mp_o.Nodes:
            node.SetSolutionStepValue(KM.PARTITION_INDEX, self.my_pid)
        for node in mp_d_m.Nodes:
            node.SetSolutionStepValue(KM.PARTITION_INDEX, self.my_pid)
        for node in mp_d_nm.Nodes:
            node.SetSolutionStepValue(KM.PARTITION_INDEX, self.my_pid)
        for node in mp_one_n.Nodes:
            node.SetSolutionStepValue(KM.PARTITION_INDEX, self.my_pid)

        origin_data_settings_scalar = KM.Parameters("""{
            "model_part_name" : "mp_origin",
            "variable_name"   : "PRESSURE"
        }""")
        origin_data_settings_scalar_non_hist = KM.Parameters("""{
            "model_part_name" : "mp_origin",
            "variable_name"   : "LAMBDA",
            "location"        : "node_non_historical"
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
        self.origin_data_scalar_non_hist = CouplingInterfaceData(origin_data_settings_scalar_non_hist, self.model)
        self.origin_data_vector = CouplingInterfaceData(origin_data_settings_vector, self.model)
        self.origin_data_single_node = CouplingInterfaceData(origin_data_settings_single_node, self.model)


        destination_matching_data_settings_scalar = KM.Parameters("""{
            "model_part_name" : "mp_destination_matching",
            "variable_name"   : "TEMPERATURE"
        }""")
        destination_matching_data_settings_scalar_non_hist = KM.Parameters("""{
            "model_part_name" : "mp_destination_matching",
            "variable_name"   : "VISCOSITY",
            "location"        : "node_non_historical"
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
        self.destination_matching_data_scalar_non_hist = CouplingInterfaceData(destination_matching_data_settings_scalar_non_hist, self.model)
        self.destination_matching_data_vector = CouplingInterfaceData(destination_matching_data_settings_vector, self.model)
        self.destination_data_single_node = CouplingInterfaceData(destination_data_settings_single_node, self.model)


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

        if KM.IsDistributedRun():
            KratosMPI.ParallelFillCommunicator(mp_o, self.data_comm).Execute()
            KratosMPI.ParallelFillCommunicator(mp_d_m, self.data_comm).Execute()
            KratosMPI.ParallelFillCommunicator(mp_d_nm, self.data_comm).Execute()
            KratosMPI.ParallelFillCommunicator(mp_one_n, self.data_comm).Execute()

    def test_copy_transfer_operator(self):
        data_transfer_op_settings = KM.Parameters("""{
            "type" : "copy"
        }""")

        data_transfer_op = data_transfer_operator_factory.CreateDataTransferOperator(data_transfer_op_settings, KM.Testing.GetDefaultDataCommunicator())

        self.__TestTransferMatching(data_transfer_op)
        self.__TestTransferMatchingSwapSign(data_transfer_op)
        self.__TestTransferMatchingAddValues(data_transfer_op)
        self.__TestTransferMatchingAddValuesAndSwapSign(data_transfer_op)

        transfer_options_empty = KM.Parameters(""" [] """)
        exp_error = 'The sizes of the data are not matching: {} != {} for interface data "{}" of solver "{}" and interface data "{}" of solver "{}"!'.format(self.origin_data_scalar.Size(), self.destination_non_matching_data_scalar.Size(), self.origin_data_scalar.name, self.origin_data_scalar.solver_name, self.destination_non_matching_data_scalar.name, self.destination_non_matching_data_scalar.solver_name)
        with self.assertRaisesRegex(Exception, exp_error):
            data_transfer_op.TransferData(self.origin_data_scalar, self.destination_non_matching_data_scalar, transfer_options_empty)

    def test_kratos_mapping_transfer_operator(self):
        if not mapping_app_available:
            self.skipTest("MappingApplication not available!")

        data_transfer_op_settings_missing = KM.Parameters("""{
            "type" : "kratos_mapping"
        }""")

        exp_error = 'No "mapper_settings" provided!'
        with self.assertRaisesRegex(Exception, exp_error):
            data_transfer_operator_factory.CreateDataTransferOperator(data_transfer_op_settings_missing, KM.Testing.GetDefaultDataCommunicator())

        data_transfer_op_settings = KM.Parameters("""{
            "type" : "kratos_mapping",
            "mapper_settings" : {
                "mapper_type" : "nearest_neighbor"
            }
        }""")

        data_transfer_op = data_transfer_operator_factory.CreateDataTransferOperator(data_transfer_op_settings, KM.Testing.GetDefaultDataCommunicator())

        self.__TestTransferMatching(data_transfer_op)
        self.__TestTransferMatchingSwapSign(data_transfer_op)
        self.__TestTransferMatchingAddValues(data_transfer_op)
        self.__TestTransferMatchingAddValuesAndSwapSign(data_transfer_op)

        # with this we make sure that only one mapper is created (and not several ones for each mapping operation!)
        # Hint: requires access to private member
        self.assertEqual(len(data_transfer_op._mappers), 1)

        self.__TestTransferMatchingInverse(data_transfer_op)
        # here we check explicitly the InverseMap fct
        self.assertEqual(len(data_transfer_op._mappers), 1)

        transfer_options_empty = KM.Parameters(""" [] """)
        data_transfer_op.TransferData(self.origin_data_scalar, self.destination_non_matching_data_scalar, transfer_options_empty)
        # here we check explicitly the creation of a second mapper, which is required since the interfaces are not the same this time
        self.assertEqual(len(data_transfer_op._mappers), 2)

        data_settings_model_part = KM.Parameters("""{
            "model_part_name" : "mp_single_node",
            "variable_name"   : "TEMPERATURE",
            "location"        : "model_part"
        }""")

        data_model_part = CouplingInterfaceData(data_settings_model_part, self.model)

        with self.assertRaisesRegex(Exception, 'Mapping only supports nodal values!'):
            data_transfer_op.TransferData(self.origin_data_scalar, data_model_part, transfer_options_empty)

    def test_kratos_mapping_transfer_operator_non_historical(self):
        if not mapping_app_available:
            self.skipTest("MappingApplication not available!")

        data_transfer_op_settings = KM.Parameters("""{
            "type" : "kratos_mapping",
            "mapper_settings" : {
                "mapper_type" : "nearest_neighbor"
            }
        }""")

        data_transfer_op = data_transfer_operator_factory.CreateDataTransferOperator(data_transfer_op_settings, KM.Testing.GetDefaultDataCommunicator())

        transfer_options_empty = KM.Parameters(""" [] """)
        # origin is non-hist
        data_transfer_op.TransferData(self.origin_data_scalar_non_hist, self.destination_matching_data_scalar, transfer_options_empty)
        self.assertVectorAlmostEqual(self.origin_data_scalar_non_hist.GetData(), self.destination_matching_data_scalar.GetData())

        # destination is non-hist
        data_transfer_op.TransferData(self.origin_data_scalar, self.destination_matching_data_scalar_non_hist, transfer_options_empty)
        self.assertVectorAlmostEqual(self.origin_data_scalar.GetData(), self.destination_matching_data_scalar_non_hist.GetData())

        # both are non-hist
        data_transfer_op.TransferData(self.origin_data_scalar_non_hist, self.destination_matching_data_scalar_non_hist, transfer_options_empty)
        self.assertVectorAlmostEqual(self.origin_data_scalar_non_hist.GetData(), self.destination_matching_data_scalar_non_hist.GetData())

        # with this we make sure that only one mapper is created (and not several ones for each mapping operation!)
        # Hint: requires access to private member
        self.assertEqual(len(data_transfer_op._mappers), 1)

        self.destination_matching_data_scalar.SetData([i+9.47 for i in range(self.destination_matching_data_scalar.Size())]) # reset values
        self.destination_matching_data_scalar_non_hist.SetData([i-89.14 for i in range(self.destination_matching_data_scalar_non_hist.Size())]) # reset values

        # Now also check InverseMap
        # origin is non-hist
        data_transfer_op.TransferData(self.destination_matching_data_scalar, self.origin_data_scalar_non_hist, transfer_options_empty)
        self.assertVectorAlmostEqual(self.destination_matching_data_scalar.GetData(), self.origin_data_scalar_non_hist.GetData())

        # destination is non-hist
        data_transfer_op.TransferData(self.destination_matching_data_scalar_non_hist, self.origin_data_scalar, transfer_options_empty)
        self.assertVectorAlmostEqual(self.destination_matching_data_scalar_non_hist.GetData(), self.origin_data_scalar.GetData())

        # both are non-hist
        data_transfer_op.TransferData(self.destination_matching_data_scalar_non_hist, self.origin_data_scalar_non_hist, transfer_options_empty)
        self.assertVectorAlmostEqual(self.destination_matching_data_scalar_non_hist.GetData(), self.origin_data_scalar_non_hist.GetData())

        # with this we make sure that only one mapper is created (and not several ones for each mapping operation!)
        # Hint: requires access to private member
        self.assertEqual(len(data_transfer_op._mappers), 1)


    def test_copy_single_to_dist_transfer_operator(self):
        data_transfer_op_settings = KM.Parameters("""{
            "type" : "sum_many_to_many"
        }""")

        data_transfer_op = data_transfer_operator_factory.CreateDataTransferOperator(data_transfer_op_settings, KM.Testing.GetDefaultDataCommunicator())
        transfer_options = KM.Parameters(""" [] """)

        for node in self.origin_data_single_node.GetModelPart().Nodes:
            node.SetSolutionStepValue(KMC.SCALAR_DISPLACEMENT, 0, 100.0)

        data_transfer_op.TransferData(self.origin_data_single_node,
                                      self.destination_matching_data_scalar,
                                      transfer_options)

        self.__CompareScalarNodalValues(self.destination_matching_data_scalar.GetModelPart().Nodes,
                                        self.origin_data_single_node.GetModelPart().Nodes,
                                        KM.TEMPERATURE, KMC.SCALAR_DISPLACEMENT,1)

        # Test with false input
        # with self.assertRaisesRegex(Exception, 'Interface data "default" of solver "default_solver" requires to be of size 1, got: 5'):
        #     data_transfer_op.TransferData(self.origin_data_scalar,
        #                                   self.destination_matching_data_scalar,
        #                                   transfer_options)


    def test_copy_single_to_dist_transfer_operator_swap_sign(self):
        data_transfer_op_settings = KM.Parameters("""{
            "type" : "sum_many_to_many"
        }""")

        data_transfer_op = data_transfer_operator_factory.CreateDataTransferOperator(data_transfer_op_settings, KM.Testing.GetDefaultDataCommunicator())
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
            "type" : "sum_many_to_many"
        }""")

        data_transfer_op = data_transfer_operator_factory.CreateDataTransferOperator(data_transfer_op_settings, KM.Testing.GetDefaultDataCommunicator())
        transfer_options = KM.Parameters(""" ["distribute_values"] """)

        for node in self.origin_data_single_node.GetModelPart().Nodes:
            node.SetSolutionStepValue(KMC.SCALAR_DISPLACEMENT, 0, 100.0)

        data_transfer_op.TransferData(self.origin_data_single_node,
                                      self.destination_matching_data_scalar,
                                      transfer_options)

        self.__CompareScalarNodalValues(self.destination_matching_data_scalar.GetModelPart().Nodes,
                                        self.origin_data_single_node.GetModelPart().Nodes,
                                        KM.TEMPERATURE, KMC.SCALAR_DISPLACEMENT, 1/self.data_comm.SumAll(self.destination_matching_data_scalar.Size()))

    def test_copy_single_to_dist_transfer_operator_add_values(self):
        data_transfer_op_settings = KM.Parameters("""{
            "type" : "sum_many_to_many"
        }""")

        data_transfer_op = data_transfer_operator_factory.CreateDataTransferOperator(data_transfer_op_settings, KM.Testing.GetDefaultDataCommunicator())
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
            "type" : "sum_many_to_many"
        }""")

        data_transfer_op = data_transfer_operator_factory.CreateDataTransferOperator(data_transfer_op_settings, KM.Testing.GetDefaultDataCommunicator())
        transfer_options = KM.Parameters(""" ["distribute_values", "swap_sign"] """)

        for node in self.origin_data_single_node.GetModelPart().Nodes:
            node.SetSolutionStepValue(KMC.SCALAR_DISPLACEMENT, 0, 100.0)

        data_transfer_op.TransferData(self.origin_data_single_node,
                                      self.destination_matching_data_scalar,
                                      transfer_options)

        self.__CompareScalarNodalValues(self.destination_matching_data_scalar.GetModelPart().Nodes,
                                        self.origin_data_single_node.GetModelPart().Nodes,
                                        KM.TEMPERATURE, KMC.SCALAR_DISPLACEMENT, -1/self.data_comm.SumAll(self.destination_matching_data_scalar.Size()))

    def test_copy_single_to_dist_transfer_operator_distribute_values_swap_sign_add_values(self):
        data_transfer_op_settings = KM.Parameters("""{
            "type" : "sum_many_to_many"
        }""")

        data_transfer_op = data_transfer_operator_factory.CreateDataTransferOperator(data_transfer_op_settings, KM.Testing.GetDefaultDataCommunicator())
        transfer_options = KM.Parameters(""" ["distribute_values", "swap_sign"] """)
        transfer_options_with_add_vals = KM.Parameters(""" ["distribute_values", "swap_sign", "add_values"] """)

        for node in self.origin_data_single_node.GetModelPart().Nodes:
            node.SetSolutionStepValue(KMC.SCALAR_DISPLACEMENT, 0, 100.0)

        data_transfer_op.TransferData(self.origin_data_single_node,
                                      self.destination_matching_data_scalar,
                                      transfer_options)

        self.__CompareScalarNodalValues(self.destination_matching_data_scalar.GetModelPart().Nodes,
                                        self.origin_data_single_node.GetModelPart().Nodes,
                                        KM.TEMPERATURE, KMC.SCALAR_DISPLACEMENT, -1/self.data_comm.SumAll(self.destination_matching_data_scalar.Size()))

        data_transfer_op.TransferData(self.origin_data_single_node,
                                      self.destination_matching_data_scalar,
                                      transfer_options_with_add_vals)

        self.__CompareScalarNodalValues(self.destination_matching_data_scalar.GetModelPart().Nodes,
                                        self.origin_data_single_node.GetModelPart().Nodes,
                                        KM.TEMPERATURE, KMC.SCALAR_DISPLACEMENT, -2/self.data_comm.SumAll(self.destination_matching_data_scalar.Size()))

    def test_sum_dist_to_single(self):
        data_transfer_op_settings = KM.Parameters("""{
            "type" : "sum_many_to_many"
        }""")

        data_transfer_op = data_transfer_operator_factory.CreateDataTransferOperator(data_transfer_op_settings, KM.Testing.GetDefaultDataCommunicator())
        transfer_options = KM.Parameters(""" [] """)

        for node in self.origin_data_scalar.GetModelPart().Nodes:
            node.SetSolutionStepValue(KM.PRESSURE, 0, 100.0)

        data_transfer_op.TransferData(self.origin_data_scalar,
                                      self.destination_data_single_node,
                                      transfer_options)

        self.__CompareScalarNodalValues(self.destination_data_single_node.GetModelPart().Nodes,
                                        self.origin_data_scalar.GetModelPart().Nodes,
                                        KMC.SCALAR_FORCE, KM.PRESSURE, self.data_comm.SumAll(self.destination_matching_data_scalar.Size()))

        # Test with false input
        # with self.assertRaisesRegex(Exception, 'Interface data "default" of solver "default_solver" requires to be of size 1, got: 5'):
        #     data_transfer_op.TransferData(self.origin_data_scalar,
        #                                 self.destination_matching_data_scalar,
        #                                 transfer_options)


    def test_sum_dist_to_single_swap_sign(self):
        data_transfer_op_settings = KM.Parameters("""{
            "type" : "sum_many_to_many"
        }""")

        data_transfer_op = data_transfer_operator_factory.CreateDataTransferOperator(data_transfer_op_settings, KM.Testing.GetDefaultDataCommunicator())
        transfer_options = KM.Parameters(""" ["swap_sign"] """)

        for node in self.origin_data_scalar.GetModelPart().Nodes:
            node.SetSolutionStepValue(KM.PRESSURE, 0, 100.0)

        data_transfer_op.TransferData(self.origin_data_scalar,
                                      self.destination_data_single_node,
                                      transfer_options)

        self.__CompareScalarNodalValues(self.destination_data_single_node.GetModelPart().Nodes,
                                        self.origin_data_scalar.GetModelPart().Nodes,
                                        KMC.SCALAR_FORCE, KM.PRESSURE, -self.data_comm.SumAll(self.destination_matching_data_scalar.Size()))


    def test_sum_dist_to_single_add_values(self):
        data_transfer_op_settings = KM.Parameters("""{
            "type" : "sum_many_to_many"
        }""")

        data_transfer_op = data_transfer_operator_factory.CreateDataTransferOperator(data_transfer_op_settings, KM.Testing.GetDefaultDataCommunicator())
        transfer_options_empty = KM.Parameters(""" [] """)
        transfer_options_add = KM.Parameters(""" ["add_values"] """)

        for node in self.origin_data_scalar.GetModelPart().Nodes:
            node.SetSolutionStepValue(KM.PRESSURE, 0, 100.0)

        data_transfer_op.TransferData(self.origin_data_scalar,
                                      self.destination_data_single_node,
                                      transfer_options_empty)

        self.__CompareScalarNodalValues(self.destination_data_single_node.GetModelPart().Nodes,
                                        self.origin_data_scalar.GetModelPart().Nodes,
                                        KMC.SCALAR_FORCE, KM.PRESSURE, self.data_comm.SumAll(self.destination_matching_data_scalar.Size()) )

        data_transfer_op.TransferData(self.origin_data_scalar,
                                            self.destination_data_single_node,
                                            transfer_options_add)

        self.__CompareScalarNodalValues(self.destination_data_single_node.GetModelPart().Nodes,
                                        self.origin_data_scalar.GetModelPart().Nodes,
                                        KMC.SCALAR_FORCE, KM.PRESSURE, 2*self.data_comm.SumAll(self.destination_matching_data_scalar.Size()) )


    def test_sum_dist_to_single_add_values_swap_sign(self):
        data_transfer_op_settings = KM.Parameters("""{
            "type" : "sum_many_to_many"
        }""")

        data_transfer_op = data_transfer_operator_factory.CreateDataTransferOperator(data_transfer_op_settings, KM.Testing.GetDefaultDataCommunicator())
        transfer_options_empty = KM.Parameters(""" [] """)
        transfer_options_add = KM.Parameters(""" ["add_values", "swap_sign"] """)

        for node in self.origin_data_scalar.GetModelPart().Nodes:
            node.SetSolutionStepValue(KM.PRESSURE, 0, 100.0)

        data_transfer_op.TransferData(self.origin_data_scalar,
                                      self.destination_data_single_node,
                                      transfer_options_empty)

        self.__CompareScalarNodalValues(self.destination_data_single_node.GetModelPart().Nodes,
                                        self.origin_data_scalar.GetModelPart().Nodes,
                                        KMC.SCALAR_FORCE, KM.PRESSURE, self.data_comm.SumAll(self.destination_matching_data_scalar.Size()))

        data_transfer_op.TransferData(self.origin_data_scalar,
                                      self.destination_data_single_node,
                                      transfer_options_add)

        self.__CompareScalarNodalValues(self.destination_data_single_node.GetModelPart().Nodes,
                                        self.origin_data_scalar.GetModelPart().Nodes,
                                        KMC.SCALAR_FORCE, KM.PRESSURE, 0.0)

    def test_sum_dist_to_dist(self):
        data_transfer_op_settings = KM.Parameters("""{
            "type" : "sum_many_to_many"
        }""")

        data_transfer_op = data_transfer_operator_factory.CreateDataTransferOperator(data_transfer_op_settings, KM.Testing.GetDefaultDataCommunicator())
        transfer_options = KM.Parameters(""" [] """)

        for node in self.origin_data_scalar.GetModelPart().Nodes:
            node.SetSolutionStepValue(KM.PRESSURE, 0, 100.0)

        data_transfer_op.TransferData(self.origin_data_scalar,
                                      self.destination_non_matching_data_scalar,
                                      transfer_options)      

        self.__CompareScalarNodalValues(self.destination_non_matching_data_scalar.GetModelPart().Nodes,
                                        self.origin_data_scalar.GetModelPart().Nodes,
                                        KM.TEMPERATURE, KM.PRESSURE, self.data_comm.SumAll(self.origin_data_scalar.Size()))
     
        # Test with false input
        # with self.assertRaisesRegex(Exception, 'Interface data "default" of solver "default_solver" requires to be of size 1, got: 5'):
        #     data_transfer_op.TransferData(self.origin_data_scalar,
        #                                 self.destination_matching_data_scalar,
        #                                 transfer_options)

    def test_sum_dist_to_dist_distribute_values(self):
        data_transfer_op_settings = KM.Parameters("""{
            "type" : "sum_many_to_many"
        }""")

        data_transfer_op = data_transfer_operator_factory.CreateDataTransferOperator(data_transfer_op_settings, KM.Testing.GetDefaultDataCommunicator())
        transfer_options = KM.Parameters(""" ["distribute_values"] """)

        for node in self.origin_data_scalar.GetModelPart().Nodes:
            node.SetSolutionStepValue(KM.PRESSURE, 0, 100.0)

        data_transfer_op.TransferData(self.origin_data_scalar,
                                      self.destination_non_matching_data_scalar,
                                      transfer_options)      

        self.__CompareScalarNodalValues(self.destination_non_matching_data_scalar.GetModelPart().Nodes,
                                        self.origin_data_scalar.GetModelPart().Nodes,
                                        KM.TEMPERATURE, KM.PRESSURE, self.data_comm.SumAll(self.origin_data_scalar.Size()) / self.data_comm.SumAll(self.destination_non_matching_data_scalar.Size()) )

    def test_sum_dist_to_dist_add_values(self):
        data_transfer_op_settings = KM.Parameters("""{
            "type" : "sum_many_to_many"
        }""")

        data_transfer_op = data_transfer_operator_factory.CreateDataTransferOperator(data_transfer_op_settings, KM.Testing.GetDefaultDataCommunicator())
        transfer_options_empty = KM.Parameters(""" [] """)
        transfer_options_add = KM.Parameters(""" ["add_values"] """)

        for node in self.origin_data_scalar.GetModelPart().Nodes:
            node.SetSolutionStepValue(KM.PRESSURE, 0, 100.0)

        data_transfer_op.TransferData(self.origin_data_scalar,
                                      self.destination_non_matching_data_scalar,
                                      transfer_options_empty)      

        self.__CompareScalarNodalValues(self.destination_non_matching_data_scalar.GetModelPart().Nodes,
                                        self.origin_data_scalar.GetModelPart().Nodes,
                                        KM.TEMPERATURE, KM.PRESSURE, self.data_comm.SumAll(self.origin_data_scalar.Size()) )

        data_transfer_op.TransferData(self.origin_data_scalar,
                                      self.destination_non_matching_data_scalar,
                                      transfer_options_add)      

        self.__CompareScalarNodalValues(self.destination_non_matching_data_scalar.GetModelPart().Nodes,
                                        self.origin_data_scalar.GetModelPart().Nodes,
                                        KM.TEMPERATURE, KM.PRESSURE, 2 * self.data_comm.SumAll(self.origin_data_scalar.Size()) )

    def test_sum_dist_to_dist_swap_sign(self):
        data_transfer_op_settings = KM.Parameters("""{
            "type" : "sum_many_to_many"
        }""")

        data_transfer_op = data_transfer_operator_factory.CreateDataTransferOperator(data_transfer_op_settings, KM.Testing.GetDefaultDataCommunicator())
        transfer_options = KM.Parameters(""" ["swap_sign"] """)

        for node in self.origin_data_scalar.GetModelPart().Nodes:
            node.SetSolutionStepValue(KM.PRESSURE, 0, 100.0)

        data_transfer_op.TransferData(self.origin_data_scalar,
                                      self.destination_non_matching_data_scalar,
                                      transfer_options)      

        self.__CompareScalarNodalValues(self.destination_non_matching_data_scalar.GetModelPart().Nodes,
                                        self.origin_data_scalar.GetModelPart().Nodes,
                                        KM.TEMPERATURE, KM.PRESSURE, -self.data_comm.SumAll(self.origin_data_scalar.Size()))

    def test_sum_dist_to_dist_distribute_values_swap_sign(self):
        data_transfer_op_settings = KM.Parameters("""{
            "type" : "sum_many_to_many"
        }""")

        data_transfer_op = data_transfer_operator_factory.CreateDataTransferOperator(data_transfer_op_settings, KM.Testing.GetDefaultDataCommunicator())
        transfer_options = KM.Parameters(""" ["distribute_values", "swap_sign"] """)

        for node in self.origin_data_scalar.GetModelPart().Nodes:
            node.SetSolutionStepValue(KM.PRESSURE, 0, 100.0)

        data_transfer_op.TransferData(self.origin_data_scalar,
                                      self.destination_non_matching_data_scalar,
                                      transfer_options)      

        self.__CompareScalarNodalValues(self.destination_non_matching_data_scalar.GetModelPart().Nodes,
                                        self.origin_data_scalar.GetModelPart().Nodes,
                                        KM.TEMPERATURE, KM.PRESSURE, -self.data_comm.SumAll(self.origin_data_scalar.Size()) / self.data_comm.SumAll(self.destination_non_matching_data_scalar.Size()) )

    def test_sum_dist_to_dist_distribute_values_swap_sign_add_values(self):
        data_transfer_op_settings = KM.Parameters("""{
            "type" : "sum_many_to_many"
        }""")

        data_transfer_op = data_transfer_operator_factory.CreateDataTransferOperator(data_transfer_op_settings, KM.Testing.GetDefaultDataCommunicator())
        transfer_options_empty = KM.Parameters(""" ["distribute_values", "swap_sign"] """)
        transfer_options_add = KM.Parameters(""" ["distribute_values", "swap_sign", "add_values"] """)

        for node in self.origin_data_scalar.GetModelPart().Nodes:
            node.SetSolutionStepValue(KM.PRESSURE, 0, 100.0)

        data_transfer_op.TransferData(self.origin_data_scalar,
                                      self.destination_non_matching_data_scalar,
                                      transfer_options_empty)      

        self.__CompareScalarNodalValues(self.destination_non_matching_data_scalar.GetModelPart().Nodes,
                                        self.origin_data_scalar.GetModelPart().Nodes,
                                        KM.TEMPERATURE, KM.PRESSURE, -self.data_comm.SumAll(self.origin_data_scalar.Size()) / self.data_comm.SumAll(self.destination_non_matching_data_scalar.Size()) )

        data_transfer_op.TransferData(self.origin_data_scalar,
                                      self.destination_non_matching_data_scalar,
                                      transfer_options_add)      

        self.__CompareScalarNodalValues(self.destination_non_matching_data_scalar.GetModelPart().Nodes,
                                        self.origin_data_scalar.GetModelPart().Nodes,
                                        KM.TEMPERATURE, KM.PRESSURE, -2 * self.data_comm.SumAll(self.origin_data_scalar.Size()) / self.data_comm.SumAll(self.destination_non_matching_data_scalar.Size()) )
        
    def test_sum_dist_to_single_check_var(self):
        data_transfer_op_settings = KM.Parameters("""{
            "type" : "sum_many_to_many"
        }""")

        data_transfer_op = data_transfer_operator_factory.CreateDataTransferOperator(data_transfer_op_settings, KM.Testing.GetDefaultDataCommunicator())
        transfer_options = KM.Parameters(""" [] """)

        with self.assertRaisesRegex(Exception, 'Variable of interface data "default" of solver "default_solver" has to be a scalar!'):
            data_transfer_op.TransferData(self.destination_matching_data_vector,
                                        self.destination_data_single_node,
                                        transfer_options)

    def test_copy_single_to_dist_check_var(self):
        data_transfer_op_settings = KM.Parameters("""{
            "type" : "sum_many_to_many"
        }""")

        data_transfer_op = data_transfer_operator_factory.CreateDataTransferOperator(data_transfer_op_settings, KM.Testing.GetDefaultDataCommunicator())
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

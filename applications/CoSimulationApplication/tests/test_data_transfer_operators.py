from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

import KratosMultiphysics as KM
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

        mp_o.AddNodalSolutionStepVariable(KM.PRESSURE)
        mp_o.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
        mp_d_m.AddNodalSolutionStepVariable(KM.TEMPERATURE)
        mp_d_m.AddNodalSolutionStepVariable(KM.FORCE)
        mp_d_nm.AddNodalSolutionStepVariable(KM.TEMPERATURE)
        mp_d_nm.AddNodalSolutionStepVariable(KM.FORCE)

        mp_o.ProcessInfo[KM.DOMAIN_SIZE] = 2
        mp_d_m.ProcessInfo[KM.DOMAIN_SIZE] = 2
        mp_d_nm.ProcessInfo[KM.DOMAIN_SIZE] = 2

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


        origin_data_settings_scalar = KM.Parameters("""{
            "model_part_name" : "mp_origin",
            "variable_name"   : "PRESSURE"
        }""")
        origin_data_settings_vector = KM.Parameters("""{
            "model_part_name" : "mp_origin",
            "variable_name"   : "DISPLACEMENT",
            "dimension" : 2
        }""")

        self.origin_data_scalar = CouplingInterfaceData(origin_data_settings_scalar, self.model)
        self.origin_data_vector = CouplingInterfaceData(origin_data_settings_vector, self.model)
        self.origin_data_scalar.Initialize()
        self.origin_data_vector.Initialize()


        destination_matching_data_settings_scalar = KM.Parameters("""{
            "model_part_name" : "mp_destination_matching",
            "variable_name"   : "TEMPERATURE"
        }""")
        destination_matching_data_settings_vector = KM.Parameters("""{
            "model_part_name" : "mp_destination_matching",
            "variable_name"   : "FORCE",
            "dimension" : 2
        }""")

        self.destination_matching_data_scalar = CouplingInterfaceData(destination_matching_data_settings_scalar, self.model)
        self.destination_matching_data_vector = CouplingInterfaceData(destination_matching_data_settings_vector, self.model)
        self.destination_matching_data_scalar.Initialize()
        self.destination_matching_data_vector.Initialize()


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
        self.__TestTransferMatchingAddValuesAdSwapSign(data_transfer_op)

        transfer_options_empty = KM.Parameters(""" [] """)
        exp_error = 'The sizes of the data are not matching: {} != {}!'.format(self.origin_data_scalar.Size(), self.destination_non_matching_data_scalar.Size())
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
        self.__TestTransferMatchingAddValuesAdSwapSign(data_transfer_op)


    def __TestTransferMatching(self, data_transfer_op):
        transfer_options_empty = KM.Parameters(""" [] """)
        data_transfer_op.TransferData(self.origin_data_scalar, self.destination_matching_data_scalar, transfer_options_empty)
        self.__CompareScalarNodalValues(self.destination_matching_data_scalar.GetModelPart().Nodes, self.origin_data_scalar.GetModelPart().Nodes, KM.TEMPERATURE, KM.PRESSURE)
        data_transfer_op.TransferData(self.origin_data_vector, self.destination_matching_data_vector, transfer_options_empty)
        self.__CompareVectorNodalValues(self.destination_matching_data_vector.GetModelPart().Nodes, self.origin_data_vector.GetModelPart().Nodes, KM.FORCE, KM.DISPLACEMENT, 2)

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

    def __TestTransferMatchingAddValuesAdSwapSign(self, data_transfer_op):
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

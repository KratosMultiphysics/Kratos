from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

import KratosMultiphysics as KM
import KratosMultiphysics.KratosUnittest as KratosUnittest

from KratosMultiphysics.CoSimulationApplication.coupling_interface_data import CouplingInterfaceData

# The expected definitions are here to make the handling of the
# multiline-stings easier (no need to deal with indentation)
coupling_interface_data_str = '''CouplingInterfaceData:
	ModelPart: "mp_4_test"
	IsDistributed: False
	Variable: "DISPLACEMENT" (Vector with dimension: 2)
	Location: "node_historical"
	Size: 10
'''

class TestCouplingInterfaceData(KratosUnittest.TestCase):

    def setUp(self):
        self.model = KM.Model()
        self.mp = self.model.CreateModelPart("mp_4_test")
        self.mp.AddNodalSolutionStepVariable(KM.PRESSURE)
        self.mp.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
        self.mp.ProcessInfo[KM.DOMAIN_SIZE] = 2
        self.mp.SetBufferSize(2)

        num_nodes = 5
        num_elems = 7
        num_conds = 3

        props = self.mp.CreateNewProperties(2)

        for i in range(num_nodes):
            node_id = i+1
            node = self.mp.CreateNewNode(node_id, 0.0, 0.0, i+1)

            node.SetSolutionStepValue(KM.PRESSURE, node_id*1.5)
            node.SetSolutionStepValue(KM.DISPLACEMENT, [node_id*1.1, node_id*1.1+2.5, node_id*2.3])

            node.SetSolutionStepValue(KM.PRESSURE, 1, node_id*1.8)
            node.SetSolutionStepValue(KM.DISPLACEMENT, 1, [node_id*5.4, node_id+3.8, node_id*24.4])

            node.SetValue(KM.TEMPERATURE, node_id*6.7)
            node.SetValue(KM.VELOCITY, [node_id+11.7, node_id-0.1, node_id*33.1])

        for i in range(num_elems):
            elem_id = i+1
            elem = self.mp.CreateNewElement("Element2D2N", elem_id, [1,2], props)

            elem.SetValue(KM.DENSITY, elem_id*6.7)
            elem.SetValue(KM.FORCE, [elem_id+11.7, elem_id-0.1, elem_id*33.1])

        for i in range(num_conds):
            cond_id = i+1
            cond = self.mp.CreateNewCondition("LineCondition2D2N", cond_id, [1,2], props)

            cond.SetValue(KM.YOUNG_MODULUS, cond_id*6.7)
            cond.SetValue(KM.ROTATION, [cond_id+11.7, cond_id-0.1, cond_id*33.1])

        self.mp[KM.NODAL_MASS] = -55.2
        self.mp.ProcessInfo[KM.TORQUE] = [1.42, 3.7, -15.4]

        self.mp.ProcessInfo[KM.NODAL_ERROR] = 15.4
        self.mp.ProcessInfo[KM.TORQUE] = [101.2, -33.7, 1.4]

    def test_basics(self):
        settings_scal_hist = KM.Parameters("""{
            "model_part_name" : "mp_4_test",
            "variable_name"   : "PRESSURE"
        }""")

        settings_vec_elem = KM.Parameters("""{
            "model_part_name" : "mp_4_test",
            "variable_name"   : "DISPLACEMENT",
            "location"        : "element",
            "dimension"       : 2
        }""")


        coupling_data_scal = CouplingInterfaceData(settings_scal_hist, self.model)
        coupling_data_vec = CouplingInterfaceData(settings_vec_elem, self.model)

        self.assertEqual(coupling_data_scal.GetModelPart().Name, "mp_4_test")
        self.assertEqual(coupling_data_vec.GetModelPart().Name, "mp_4_test")

        self.assertFalse(coupling_data_scal.GetModelPart().IsDistributed())
        self.assertFalse(coupling_data_vec.GetModelPart().IsDistributed())

        self.assertEqual(coupling_data_scal.Size(), 5) # 5 nodes and scalar var
        self.assertEqual(coupling_data_vec.Size(), 14) # 7 elements and vector var with dim==2

        self.assertEqual(coupling_data_scal.GetBufferSize(), 2)
        self.assertEqual(coupling_data_vec.GetBufferSize(), 1)

    def test_printing(self):
        settings = KM.Parameters("""{
            "model_part_name" : "mp_4_test",
            "variable_name"   : "DISPLACEMENT",
            "dimension"       : 2
        }""")

        coupling_data = CouplingInterfaceData(settings, self.model)
        self.assertMultiLineEqual(str(coupling_data), coupling_interface_data_str)

    def test_wrong_input_dim_scalar(self):
        settings = KM.Parameters("""{
            "model_part_name" : "mp_4_test",
            "variable_name"   : "PRESSURE",
            "dimension"       : 2
        }""")

        with self.assertRaisesRegex(Exception, '"dimension" cannot be specifed for scalar variables!'):
            coupling_data = CouplingInterfaceData(settings, self.model)

    def test_wrong_input_no_dim_vector(self):
        settings = KM.Parameters("""{
            "model_part_name" : "mp_4_test",
            "variable_name"   : "DISPLACEMENT"
        }""")

        with self.assertRaisesRegex(Exception, '"dimension" has to be specifed for vector variables!'):
            coupling_data = CouplingInterfaceData(settings, self.model)

    def test_wrong_input_variable_type(self):
        settings = KM.Parameters("""{
            "model_part_name" : "mp_4_test",
            "variable_name"   : "EXTERNAL_FORCES_VECTOR"
        }""")

        exp_error = 'The input for "variable" "EXTERNAL_FORCES_VECTOR" is of variable type "Vector" which is not allowed, only the following variable types are allowed:\nBool, Integer, Unsigned Integer, Double, Component, Array'

        with self.assertRaisesRegex(Exception, exp_error):
            coupling_data = CouplingInterfaceData(settings, self.model)

    def test_wrong_input_location(self):
        settings = KM.Parameters("""{
            "model_part_name" : "mp_4_test",
            "variable_name"   : "PRESSURE",
            "location"        : "dummy"
        }""")

        exp_error = '"dummy" is not allowed as "location", only the following options are possible:\nnode_historical, node_non_historical, element, condition, process_info, model_part'

        with self.assertRaisesRegex(Exception, exp_error):
            coupling_data = CouplingInterfaceData(settings, self.model)







if __name__ == '__main__':
    KratosUnittest.main()

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

model_part_scalar_value = -61.225
model_part_vector_value = [123.5, 54.9, -92.4]

process_info_scalar_value = 9745.34
process_info_vector_value = [-556.3, -334.2, 65.9]

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

            node.SetSolutionStepValue(KM.PRESSURE, NodeScalarHistValueCurrent(node_id))
            node.SetSolutionStepValue(KM.DISPLACEMENT, NodeVectorHistValueCurrent(node_id))

            node.SetSolutionStepValue(KM.PRESSURE, 1, NodeScalarHistValuePrevious(node_id))
            node.SetSolutionStepValue(KM.DISPLACEMENT, 1, NodeVectorHistValuePrevious(node_id))

            node.SetValue(KM.TEMPERATURE, NodeScalarNonHistValue(node_id))
            node.SetValue(KM.VELOCITY, NodeVectorNonHistValue(node_id))

        for i in range(num_elems):
            elem_id = i+1
            elem = self.mp.CreateNewElement("Element2D2N", elem_id, [1,2], props)

            elem.SetValue(KM.DENSITY, ElementScalarValue(elem_id))
            elem.SetValue(KM.FORCE, ElementVectorValue(elem_id))

        for i in range(num_conds):
            cond_id = i+1
            cond = self.mp.CreateNewCondition("LineCondition2D2N", cond_id, [1,2], props)

            cond.SetValue(KM.YOUNG_MODULUS, ConditionScalarValue(cond_id))
            cond.SetValue(KM.ROTATION, ConditionVectorValue(cond_id))

        self.mp[KM.NODAL_MASS] = model_part_scalar_value
        self.mp[KM.TORQUE] = model_part_vector_value

        self.mp.ProcessInfo[KM.NODAL_ERROR] = process_info_scalar_value
        self.mp.ProcessInfo[KM.NORMAL] = process_info_vector_value

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

    def test_wrong_input_set_data(self):
        settings = KM.Parameters("""{
            "model_part_name" : "mp_4_test",
            "variable_name"   : "PRESSURE"
        }""")


        settings_model_part = KM.Parameters("""{
            "model_part_name" : "mp_4_test",
            "variable_name"   : "PRESSURE",
            "location"        : "model_part"
        }""")

        coupling_data = CouplingInterfaceData(settings, self.model)
        coupling_data_mp = CouplingInterfaceData(settings_model_part, self.model)

        wrong_data = [1,2,3]
        correct_data = [1,2,3,4,5]

        with self.assertRaisesRegex(Exception, "The sizes of the data are not matching, got: 3, expected: 5"):
            coupling_data.SetData(wrong_data)

        with self.assertRaisesRegex(Exception, "The buffer-size is not large enough: current buffer size: 2 | requested solution_step_index: 3"):
            coupling_data.SetData(correct_data, 2)

        with self.assertRaisesRegex(Exception, "accessing data from previous steps is only possible with historical nodal data!"):
            coupling_data_mp.SetData(correct_data, 2)

    def test_GetSetNodalHistoricalData(self):
        settings_scal = KM.Parameters("""{
            "model_part_name" : "mp_4_test",
            "variable_name"   : "PRESSURE"
        }""")

        settings_vec = KM.Parameters("""{
            "model_part_name" : "mp_4_test",
            "variable_name"   : "DISPLACEMENT",
            "dimension"       : 2
        }""")

        coupling_data_scal = CouplingInterfaceData(settings_scal, self.model)
        coupling_data_vec = CouplingInterfaceData(settings_vec,  self.model)

        # 1. check the initial values
        exp_data_scal_cur = [NodeScalarHistValueCurrent(node.Id) for node in self.mp.Nodes]
        exp_data_scal_prev = [NodeScalarHistValuePrevious(node.Id) for node in self.mp.Nodes]

        exp_data_vec_cur = [NodeVectorHistValueCurrent(node.Id) for node in self.mp.Nodes]
        exp_data_vec_prev = [NodeVectorHistValuePrevious(node.Id) for node in self.mp.Nodes]

        self.__CheckScalarData(exp_data_scal_cur, coupling_data_scal.GetData())
        self.__CheckVectorData(exp_data_vec_cur, coupling_data_vec.GetData(), 2)

        self.__CheckScalarData(exp_data_scal_prev, coupling_data_scal.GetData(1))
        self.__CheckVectorData(exp_data_vec_prev, coupling_data_vec.GetData(1), 2)

        # 2. check setting and getting works
        set_data_scal_cur = [ElementScalarValue(node.Id) for node in self.mp.Nodes]
        set_data_scal_prev = [ConditionScalarValue(node.Id) for node in self.mp.Nodes]

        set_data_vec_cur = [ElementVectorValue(node.Id) for node in self.mp.Nodes]
        set_data_vec_prev = [ConditionVectorValue(node.Id) for node in self.mp.Nodes]

        self.__CheckSetGetScalarData(set_data_scal_cur, coupling_data_scal)
        self.__CheckSetGetVectorData(set_data_vec_cur, coupling_data_vec, 2)

        self.__CheckSetGetScalarData(set_data_scal_prev, coupling_data_scal, 1)
        self.__CheckSetGetVectorData(set_data_vec_prev, coupling_data_vec, 2, 1)

    def test_GetSetNodalNonHistoricalData(self):
        settings_scal = KM.Parameters("""{
            "model_part_name" : "mp_4_test",
            "location"        : "node_non_historical",
            "variable_name"   : "PRESSURE"
        }""")

        settings_vec = KM.Parameters("""{
            "model_part_name" : "mp_4_test",
            "variable_name"   : "DISPLACEMENT",
            "location"        : "node_non_historical",
            "dimension"       : 3
        }""")

        self.mp.ProcessInfo[KM.DOMAIN_SIZE] = 3

        coupling_data_scal = CouplingInterfaceData(settings_scal, self.model)
        coupling_data_vec = CouplingInterfaceData(settings_vec,  self.model)

        # 1. check the initial values
        exp_data_scal = [NodeScalarHistValueCurrent(node.Id) for node in self.mp.Nodes]
        exp_data_vec = [NodeVectorHistValueCurrent(node.Id) for node in self.mp.Nodes]

        self.__CheckScalarData(exp_data_scal, coupling_data_scal.GetData())
        self.__CheckVectorData(exp_data_vec, coupling_data_vec.GetData(), 3)

        # 2. check setting and getting works
        set_data_scal = [ConditionScalarValue(node.Id) for node in self.mp.Nodes]
        set_data_vec = [ConditionVectorValue(node.Id) for node in self.mp.Nodes]

        self.__CheckSetGetScalarData(set_data_scal, coupling_data_scal)
        self.__CheckSetGetVectorData(set_data_vec, coupling_data_vec, 3)

    def test_GetSetElementalData(self):
        settings_scal = KM.Parameters("""{
            "model_part_name" : "mp_4_test",
            "location"        : "element",
            "variable_name"   : "DENSITY"
        }""")

        settings_vec = KM.Parameters("""{
            "model_part_name" : "mp_4_test",
            "variable_name"   : "FORCE",
            "location"        : "element",
            "dimension"       : 2
        }""")

        coupling_data_scal = CouplingInterfaceData(settings_scal, self.model)
        coupling_data_vec = CouplingInterfaceData(settings_vec,  self.model)

        # 1. check the initial values
        exp_data_scal = [ElementScalarValue(elem.Id) for elem in self.mp.Elements]
        exp_data_vec = [ElementVectorValue(elem.Id) for elem in self.mp.Elements]

        self.__CheckScalarData(exp_data_scal, coupling_data_scal.GetData())
        self.__CheckVectorData(exp_data_vec, coupling_data_vec.GetData(), 2)

        # 2. check setting and getting works
        set_data_scal = [NodeScalarNonHistValue(elem.Id) for elem in self.mp.Elements]
        set_data_vec = [NodeVectorNonHistValue(elem.Id) for elem in self.mp.Elements]

        self.__CheckSetGetScalarData(set_data_scal, coupling_data_scal)
        self.__CheckSetGetVectorData(set_data_vec, coupling_data_vec, 2)

    def test_GetSetConditionalData(self):
        settings_scal = KM.Parameters("""{
            "model_part_name" : "mp_4_test",
            "location"        : "condition",
            "variable_name"   : "YOUNG_MODULUS"
        }""")

        settings_vec = KM.Parameters("""{
            "model_part_name" : "mp_4_test",
            "variable_name"   : "ROTATION",
            "location"        : "condition",
            "dimension"       : 2
        }""")

        coupling_data_scal = CouplingInterfaceData(settings_scal, self.model)
        coupling_data_vec = CouplingInterfaceData(settings_vec,  self.model)

        # 1. check the initial values
        exp_data_scal = [ConditionScalarValue(cond.Id) for cond in self.mp.Conditions]
        exp_data_vec = [ConditionVectorValue(cond.Id) for cond in self.mp.Conditions]

        self.__CheckScalarData(exp_data_scal, coupling_data_scal.GetData())
        self.__CheckVectorData(exp_data_vec, coupling_data_vec.GetData(), 2)

        # 2. check setting and getting works
        set_data_scal = [NodeScalarHistValuePrevious(cond.Id) for cond in self.mp.Conditions]
        set_data_vec = [NodeVectorHistValuePrevious(cond.Id) for cond in self.mp.Conditions]

        self.__CheckSetGetScalarData(set_data_scal, coupling_data_scal)
        self.__CheckSetGetVectorData(set_data_vec, coupling_data_vec, 2)

    def test_GetSetModelPartData(self):
        settings_scal = KM.Parameters("""{
            "model_part_name" : "mp_4_test",
            "location"        : "model_part",
            "variable_name"   : "NODAL_MASS"
        }""")

        settings_vec = KM.Parameters("""{
            "model_part_name" : "mp_4_test",
            "variable_name"   : "TORQUE",
            "location"        : "model_part",
            "dimension"       : 1
        }""")

        self.mp.ProcessInfo[KM.DOMAIN_SIZE] = 1

        coupling_data_scal = CouplingInterfaceData(settings_scal, self.model)
        coupling_data_vec = CouplingInterfaceData(settings_vec,  self.model)

        # 1. check the initial values
        self.__CheckScalarData([model_part_scalar_value], coupling_data_scal.GetData())
        self.__CheckVectorData([model_part_vector_value], coupling_data_vec.GetData(), 1)

        # 2. check setting and getting works
        set_data_scal = [process_info_scalar_value]
        set_data_vec = [process_info_vector_value]

        self.__CheckSetGetScalarData(set_data_scal, coupling_data_scal)
        self.__CheckSetGetVectorData(set_data_vec, coupling_data_vec, 1)

    def test_GetSetProcessInfoData(self):
        settings_scal = KM.Parameters("""{
            "model_part_name" : "mp_4_test",
            "location"        : "process_info",
            "variable_name"   : "NODAL_ERROR"
        }""")

        settings_vec = KM.Parameters("""{
            "model_part_name" : "mp_4_test",
            "variable_name"   : "NORMAL",
            "location"        : "process_info",
            "dimension"       : 2
        }""")

        coupling_data_scal = CouplingInterfaceData(settings_scal, self.model)
        coupling_data_vec = CouplingInterfaceData(settings_vec,  self.model)

        # 1. check the initial values
        self.__CheckScalarData([process_info_scalar_value], coupling_data_scal.GetData())
        self.__CheckVectorData([process_info_vector_value], coupling_data_vec.GetData(), 2)

        # 2. check setting and getting works
        set_data_scal = [model_part_scalar_value]
        set_data_vec = [model_part_vector_value]

        self.__CheckSetGetScalarData(set_data_scal, coupling_data_scal)
        self.__CheckSetGetVectorData(set_data_vec, coupling_data_vec, 2)


    def __CheckScalarData(self, exp_data, data):
        self.assertEqual(len(exp_data), len(data))

        for exp_val, val in zip(exp_data, data):
            self.asserAlmostEqual(exp_val, val)

    def __CheckVectorData(self, exp_data, data, dim):
        self.assertEqual(len(exp_data), len(data))

        for exp_vec, vec in zip(exp_data, data):
            for i in range(dim):
                self.asserAlmostEqual(exp_vec[i], vec[i])

    def __CheckSetGetScalarData(self, the_data, coupling_data, solution_step_index=0):
        # Checking to call fct differently
        if solution_step_index > 0:
            coupling_data.SetData(the_data, solution_step_index)
            extracted_data = coupling_data.GetData(solution_step_index)
        else:
            coupling_data.SetData(the_data)
            extracted_data = coupling_data.GetData()

        self.__CheckScalarData(the_data, extracted_data)

    def __CheckSetGetVectorData(self, the_data, coupling_data, dim, solution_step_index=0):
        # Checking to call fct differently
        if solution_step_index > 0:
            coupling_data.SetData(the_data, solution_step_index)
            extracted_data = coupling_data.GetData(solution_step_index)
        else:
            coupling_data.SetData(the_data)
            extracted_data = coupling_data.GetData()

        self.__CheckVectorData(the_data, extracted_data, dim)


def NodeScalarHistValueCurrent(the_id):
    return the_id*1.5

def NodeScalarHistValuePrevious(the_id):
    return the_id*123.7

def NodeVectorHistValueCurrent(the_id):
    return [the_id*1.1, the_id*1.1+2.5, the_id*2.3]

def NodeVectorHistValuePrevious(the_id):
    return [the_id*1.4-0.4, the_id*10.4+2.5, the_id*5.3]

def NodeScalarNonHistValue(the_id):
    return the_id*2.1+6.2

def NodeVectorNonHistValue(the_id):
    return [the_id*14.7, the_id*19.2-10.5, the_id*303.9]

def ElementScalarValue(the_id):
    return the_id*6.4+1.1

def ElementVectorValue(the_id):
    return [the_id*11.7, the_id*12.2-120.5, the_id*3.9]

def ConditionScalarValue(the_id):
    return the_id*14.5+15.8

def ConditionVectorValue(the_id):
    return [the_id*65.7, the_id*12.1, the_id+3.9]



if __name__ == '__main__':
    KratosUnittest.main()

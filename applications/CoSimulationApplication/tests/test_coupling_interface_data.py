from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

import KratosMultiphysics as KM
import KratosMultiphysics.KratosUnittest as KratosUnittest

from KratosMultiphysics.CoSimulationApplication.coupling_interface_data import CouplingInterfaceData

from KratosMultiphysics.CoSimulationApplication.co_simulation_tools import UsingPyKratos
using_pykratos = UsingPyKratos()

# The expected definitions are here to make the handling of the
# multiline-stings easier (no need to deal with indentation)
coupling_interface_data_str = '''CouplingInterfaceData:
	Name: "default"
	SolverWrapper: "default_solver"
	ModelPart: "mp_4_test"
	IsDistributed: False
	Variable: "DISPLACEMENT" (Vector with dimension: 2)
	Location: "node_historical"
	Size: 10
'''

model_part_scalar_value = -61.225
model_part_vector_value = [123.5, 54.9, -92.4]

model_part_scalar_value_2 = 9745.34
model_part_vector_value_2 = [-556.3, -334.2, 65.9]

class TestCouplingInterfaceData(KratosUnittest.TestCase):

    def setUp(self):
        self.model = KM.Model()
        self.mp = self.model.CreateModelPart("mp_4_test", 2)
        self.mp.AddNodalSolutionStepVariable(KM.PRESSURE)
        self.mp.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
        self.mp.ProcessInfo[KM.DOMAIN_SIZE] = 2

        num_nodes = 5
        num_elems = 7
        num_conds = 3

        props = self.mp.CreateNewProperties(2)

        for i in range(num_nodes):
            node_id = i+1
            node = self.mp.CreateNewNode(node_id, 0.0, 0.0, i+1)

            node.SetSolutionStepValue(KM.PRESSURE, 0, NodeScalarHistValueCurrent(node_id))
            node.SetSolutionStepValue(KM.DISPLACEMENT, 0, NodeVectorHistValueCurrent(node_id))

            node.SetSolutionStepValue(KM.PRESSURE, 1, NodeScalarHistValuePrevious(node_id))
            node.SetSolutionStepValue(KM.DISPLACEMENT, 1, NodeVectorHistValuePrevious(node_id))

            node.SetValue(KM.TEMPERATURE, NodeScalarNonHistValue(node_id))
            node.SetValue(KM.VELOCITY, NodeVectorNonHistValue(node_id))

        for i in range(num_elems):
            elem_id = i+1
            elem = self.mp.CreateNewElement("Element2D2N", elem_id, [1,2], props)

            elem.SetValue(KM.DENSITY, ElementScalarValue(elem_id))
            elem.SetValue(KM.FORCE, ElementVectorValue(elem_id))

        if not using_pykratos: # pyKratos does not have Conditions
            for i in range(num_conds):
                cond_id = i+1
                cond = self.mp.CreateNewCondition("LineCondition2D2N", cond_id, [1,2], props)

                cond.SetValue(KM.YOUNG_MODULUS, ConditionScalarValue(cond_id))
                cond.SetValue(KM.ROTATION, ConditionVectorValue(cond_id))

        self.mp[KM.NODAL_MASS] = model_part_scalar_value
        self.mp[KM.TORQUE] = model_part_vector_value

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
        coupling_data_scal.Initialize()
        coupling_data_vec = CouplingInterfaceData(settings_vec_elem, self.model)
        coupling_data_vec.Initialize()

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
        coupling_data.Initialize()
        self.assertMultiLineEqual(str(coupling_data), coupling_interface_data_str)

    def test_without_initialization(self):
        settings = KM.Parameters("""{
            "model_part_name" : "mp_4_test",
            "variable_name"   : "DISPLACEMENT",
            "dimension"       : 2
        }""")

        coupling_data = CouplingInterfaceData(settings, self.model)
        # coupling_data.Initialize() # intentially commented to raise error
        with self.assertRaisesRegex(Exception, ' can onyl be called after initializing the CouplingInterfaceData!'):
            self.assertMultiLineEqual(str(coupling_data), coupling_interface_data_str)

        with self.assertRaisesRegex(Exception, ' can onyl be called after initializing the CouplingInterfaceData!'):
            coupling_data.PrintInfo()

        with self.assertRaisesRegex(Exception, ' can onyl be called after initializing the CouplingInterfaceData!'):
            coupling_data.GetModelPart()

        with self.assertRaisesRegex(Exception, ' can onyl be called after initializing the CouplingInterfaceData!'):
            coupling_data.IsDistributed()

        with self.assertRaisesRegex(Exception, ' can onyl be called after initializing the CouplingInterfaceData!'):
            coupling_data.Size()

        with self.assertRaisesRegex(Exception, ' can onyl be called after initializing the CouplingInterfaceData!'):
            coupling_data.GetBufferSize()

        with self.assertRaisesRegex(Exception, ' can onyl be called after initializing the CouplingInterfaceData!'):
            coupling_data.GetData()

        with self.assertRaisesRegex(Exception, ' can onyl be called after initializing the CouplingInterfaceData!'):
            coupling_data.SetData([])

    def test_unallowed_names(self):
        settings = KM.Parameters("""{
            "model_part_name" : "mp_4_test",
            "variable_name"   : "PRESSURE"
        }""")

        with self.assertRaisesRegex(Exception, 'The name cannot be empty, contain whitespaces or "."!'):
            CouplingInterfaceData(settings, self.model, "")

        with self.assertRaisesRegex(Exception, 'The name cannot be empty, contain whitespaces or "."!'):
            CouplingInterfaceData(settings, self.model, "aaa.bbbb")

        with self.assertRaisesRegex(Exception, 'The name cannot be empty, contain whitespaces or "."!'):
            CouplingInterfaceData(settings, self.model, "aaa bbb")


    def test_var_does_not_exist(self):
        settings = KM.Parameters("""{
            "model_part_name" : "mp_4_test",
            "variable_name"   : "var_that_hopefully_none_will_ever_create_otherwise_this_test_will_be_wrong"
        }""")

        with self.assertRaisesRegex(Exception, 'does not exist!'):
            CouplingInterfaceData(settings, self.model)

    def test_wrong_input_dim_scalar(self):
        settings = KM.Parameters("""{
            "model_part_name" : "mp_4_test",
            "variable_name"   : "PRESSURE",
            "dimension"       : 2
        }""")

        coupling_data = CouplingInterfaceData(settings, self.model)
        with self.assertRaisesRegex(Exception, '"dimension" cannot be specifed for scalar variables!'):
            coupling_data.Initialize()

    def test_wrong_input_no_dim_vector(self):
        settings = KM.Parameters("""{
            "model_part_name" : "mp_4_test",
            "variable_name"   : "DISPLACEMENT"
        }""")

        coupling_data = CouplingInterfaceData(settings, self.model)
        with self.assertRaisesRegex(Exception, '"dimension" has to be specifed for vector variables!'):
            coupling_data.Initialize()

    def test_wrong_input_variable_type(self):
        settings = KM.Parameters("""{
            "model_part_name" : "mp_4_test",
            "variable_name"   : "EXTERNAL_FORCES_VECTOR"
        }""")

        exp_error = 'The input for "variable" "EXTERNAL_FORCES_VECTOR" is of variable type "Vector" which is not allowed, only the following variable types are allowed:\nBool, Integer, Unsigned Integer, Double, Component, Array'

        with self.assertRaisesRegex(Exception, exp_error):
            CouplingInterfaceData(settings, self.model)

    def test_wrong_input_location(self):
        settings = KM.Parameters("""{
            "model_part_name" : "mp_4_test",
            "variable_name"   : "PRESSURE",
            "location"        : "dummy"
        }""")

        exp_error = '"dummy" is not allowed as "location", only the following options are possible:\nnode_historical, node_non_historical, element, condition, model_part'

        with self.assertRaisesRegex(Exception, exp_error):
            CouplingInterfaceData(settings, self.model)

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
        coupling_data.Initialize()
        coupling_data_mp = CouplingInterfaceData(settings_model_part, self.model)
        coupling_data_mp.Initialize()

        wrong_data = [1,2,3]
        correct_data = [1,2,3,4,5]

        with self.assertRaisesRegex(Exception, "The sizes of the data are not matching, got: 3, expected: 5"):
            coupling_data.SetData(wrong_data)

        with self.assertRaisesRegex(Exception, "The buffer-size is not large enough: current buffer size: 2 | requested solution_step_index: 3"):
            coupling_data.SetData(correct_data, 2)

        with self.assertRaisesRegex(Exception, "accessing data from previous steps is only possible with historical nodal data!"):
            coupling_data_mp.SetData(correct_data, 2)

    def test_wrong_input_dim_array(self):
        settings = KM.Parameters("""{
            "model_part_name" : "mp_4_test",
            "variable_name"   : "DISPLACEMENT",
            "dimension"       : 4
        }""")

        exp_error = '"dimension" can only be 1,2,3 when using variables of type "Array"'

        coupling_data = CouplingInterfaceData(settings, self.model)
        with self.assertRaisesRegex(Exception, exp_error):
            coupling_data.Initialize()

    def test_wrong_input_missing_solutionstepvar(self):
        settings = KM.Parameters("""{
            "model_part_name" : "mp_4_test",
            "variable_name"   : "FORCE",
            "dimension"       : 2
        }""")

        exp_error = '"FORCE" is missing as SolutionStepVariable in ModelPart "mp_4_test"'

        coupling_data = CouplingInterfaceData(settings, self.model)
        with self.assertRaisesRegex(Exception, exp_error):
            coupling_data.Initialize()

    def test_wrong_input_missing_solutionstepvar_component(self):
        settings = KM.Parameters("""{
            "model_part_name" : "mp_4_test",
            "variable_name"   : "FORCE_X"
        }""")

        exp_error = '"FORCE" is missing as SolutionStepVariable in ModelPart "mp_4_test"'

        coupling_data = CouplingInterfaceData(settings, self.model)
        with self.assertRaisesRegex(Exception, exp_error):
            coupling_data.Initialize()

    def test_non_existing_model_part(self):
        settings = KM.Parameters("""{
            "model_part_name" : "something",
            "variable_name"   : "PRESSURE",
            "location"        : "node_non_historical"
        }""")

        coupling_data = CouplingInterfaceData(settings, self.model)
        with self.assertRaisesRegex(Exception, "The specified ModelPart is not in the Model, only the following ModelParts are available:"):
            coupling_data.Initialize()

    def test_GetHistoricalVariableDict(self):
        settings = KM.Parameters("""{
            "model_part_name" : "mp_4_test",
            "variable_name"   : "PRESSURE"
        }""")

        coupling_data = CouplingInterfaceData(settings, self.model)

        exp_dict = {"mp_4_test" : KM.PRESSURE}

        dict_res = coupling_data.GetHistoricalVariableDict()

        self.assertEqual(len(exp_dict), len(dict_res))
        self.assertEqual(exp_dict["mp_4_test"].Name(), dict_res["mp_4_test"].Name())

        # this should not give anything since there are no historical variables in this case
        settings_2 = KM.Parameters("""{
            "model_part_name" : "mp_4_test",
            "variable_name"   : "PRESSURE",
            "location"        : "element"
        }""")

        coupling_data_2 = CouplingInterfaceData(settings_2, self.model)

        exp_dict_2 = {}

        dict_res_2 = coupling_data_2.GetHistoricalVariableDict()

        self.assertEqual(len(exp_dict_2), len(dict_res_2))

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
        coupling_data_scal.Initialize()
        coupling_data_vec.Initialize()

        # 1. check the initial values
        exp_data_scal_cur = [NodeScalarHistValueCurrent(node.Id) for node in self.mp.Nodes]
        exp_data_scal_prev = [NodeScalarHistValuePrevious(node.Id) for node in self.mp.Nodes]

        exp_data_vec_cur = GetVectorValues(self.mp.Nodes, NodeVectorHistValueCurrent, 2)
        exp_data_vec_prev = GetVectorValues(self.mp.Nodes, NodeVectorHistValuePrevious, 2)

        self.__CheckData(exp_data_scal_cur, coupling_data_scal.GetData())
        self.__CheckData(exp_data_vec_cur, coupling_data_vec.GetData())

        self.__CheckData(exp_data_scal_prev, coupling_data_scal.GetData(1))
        self.__CheckData(exp_data_vec_prev, coupling_data_vec.GetData(1))

        # 2. check setting and getting works
        set_data_scal_cur = [ElementScalarValue(node.Id) for node in self.mp.Nodes]
        set_data_scal_prev = [ConditionScalarValue(node.Id) for node in self.mp.Nodes]

        set_data_vec_cur = GetVectorValues(self.mp.Nodes, ElementVectorValue, 2)
        set_data_vec_prev = GetVectorValues(self.mp.Nodes, ConditionVectorValue, 2)

        self.__CheckSetGetData(set_data_scal_cur, coupling_data_scal)
        self.__CheckSetGetData(set_data_vec_cur, coupling_data_vec)

        self.__CheckSetGetData(set_data_scal_prev, coupling_data_scal, 1)
        self.__CheckSetGetData(set_data_vec_prev, coupling_data_vec, 1)

    def test_GetSetNodalNonHistoricalData(self):
        settings_scal = KM.Parameters("""{
            "model_part_name" : "mp_4_test",
            "location"        : "node_non_historical",
            "variable_name"   : "TEMPERATURE"
        }""")

        settings_vec = KM.Parameters("""{
            "model_part_name" : "mp_4_test",
            "variable_name"   : "VELOCITY",
            "location"        : "node_non_historical",
            "dimension"       : 3
        }""")

        self.mp.ProcessInfo[KM.DOMAIN_SIZE] = 3

        coupling_data_scal = CouplingInterfaceData(settings_scal, self.model)
        coupling_data_vec = CouplingInterfaceData(settings_vec,  self.model)
        coupling_data_scal.Initialize()
        coupling_data_vec.Initialize()

        # 1. check the initial values
        exp_data_scal = [NodeScalarNonHistValue(node.Id) for node in self.mp.Nodes]
        exp_data_vec = GetVectorValues(self.mp.Nodes, NodeVectorNonHistValue, 3)

        self.__CheckData(exp_data_scal, coupling_data_scal.GetData())
        self.__CheckData(exp_data_vec, coupling_data_vec.GetData())

        # 2. check setting and getting works
        set_data_scal = [ConditionScalarValue(node.Id) for node in self.mp.Nodes]
        set_data_vec = GetVectorValues(self.mp.Nodes, ConditionVectorValue, 3)

        self.__CheckSetGetData(set_data_scal, coupling_data_scal)
        self.__CheckSetGetData(set_data_vec, coupling_data_vec)

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
        coupling_data_scal.Initialize()
        coupling_data_vec.Initialize()

        # 1. check the initial values
        exp_data_scal = [ElementScalarValue(elem.Id) for elem in self.mp.Elements]
        exp_data_vec = GetVectorValues(self.mp.Elements, ElementVectorValue, 2)

        self.__CheckData(exp_data_scal, coupling_data_scal.GetData())
        self.__CheckData(exp_data_vec, coupling_data_vec.GetData())

        # 2. check setting and getting works
        set_data_scal = [NodeScalarNonHistValue(elem.Id) for elem in self.mp.Elements]
        set_data_vec = GetVectorValues(self.mp.Elements, NodeVectorNonHistValue, 2)

        self.__CheckSetGetData(set_data_scal, coupling_data_scal)
        self.__CheckSetGetData(set_data_vec, coupling_data_vec)

    def test_GetSetConditionalData(self):
        if using_pykratos:
            self.skipTest("This test cannot be run with pyKratos!")
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
        coupling_data_scal.Initialize()
        coupling_data_vec.Initialize()

        # 1. check the initial values
        exp_data_scal = [ConditionScalarValue(cond.Id) for cond in self.mp.Conditions]
        exp_data_vec = GetVectorValues(self.mp.Conditions, ConditionVectorValue, 2)

        self.__CheckData(exp_data_scal, coupling_data_scal.GetData())
        self.__CheckData(exp_data_vec, coupling_data_vec.GetData())

        # 2. check setting and getting works
        set_data_scal = [NodeScalarHistValuePrevious(cond.Id) for cond in self.mp.Conditions]
        set_data_vec = GetVectorValues(self.mp.Conditions, NodeVectorHistValuePrevious, 2)

        self.__CheckSetGetData(set_data_scal, coupling_data_scal)
        self.__CheckSetGetData(set_data_vec, coupling_data_vec)

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
        coupling_data_scal.Initialize()
        coupling_data_vec.Initialize()

        # 1. check the initial values
        self.__CheckData([model_part_scalar_value], coupling_data_scal.GetData())
        self.__CheckData([model_part_vector_value[0]], coupling_data_vec.GetData())

        # 2. check setting and getting works
        set_data_scal = [model_part_scalar_value_2]
        set_data_vec = [model_part_vector_value_2[0]]

        self.__CheckSetGetData(set_data_scal, coupling_data_scal)
        self.__CheckSetGetData(set_data_vec, coupling_data_vec)

    def __CheckData(self, exp_data, data):
        self.assertEqual(len(exp_data), len(data))

        for exp_val, val in zip(exp_data, data):
            self.assertAlmostEqual(exp_val, val)

    def __CheckSetGetData(self, the_data, coupling_data, solution_step_index=0):
        # Checking to call fct differently
        if solution_step_index > 0:
            coupling_data.SetData(the_data, solution_step_index)
            extracted_data = coupling_data.GetData(solution_step_index)
        else:
            coupling_data.SetData(the_data)
            extracted_data = coupling_data.GetData()

        self.__CheckData(the_data, extracted_data)

def GetVectorValues(container, fct_ptr, dim):
    values = []
    for entity in container:
        vector_vals = fct_ptr(entity.Id)
        for i in range(dim):
            values.append(vector_vals[i])
    return values

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

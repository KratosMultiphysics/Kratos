import KratosMultiphysics as KM

import KratosMultiphysics.MappingApplication as KratosMapping
import KratosMultiphysics.KratosUnittest as KratosUnittest

'''
This test is a fast test for testing the mappers, WORKS ONLY IN SERIAL
It stems from the patch tests in the StructuralMechanicsApplication
It covers only 2D!

Setup:
ModelPart Origin:
x-----x-----x---x

Modelpart Destination:
x--x--x----x---x
'''


class TestPatchTestMappers(KratosUnittest.TestCase):
    def setUp(self):
        self.test_model = KM.Model()
        self.mp_origin = self.test_model.CreateModelPart("origin_part")
        self.mp_destination = self.test_model.CreateModelPart("destination_part")

        self._add_variables()
        self._create_nodes()
        self._create_elements()


    def _add_variables(self):
        self.mp_origin.AddNodalSolutionStepVariable(KM.PRESSURE)
        self.mp_origin.AddNodalSolutionStepVariable(KM.FORCE)

        self.mp_destination.AddNodalSolutionStepVariable(KM.TEMPERATURE)
        self.mp_destination.AddNodalSolutionStepVariable(KM.VELOCITY)


    def _create_nodes(self):
        self.mp_origin.CreateNewNode(1, -5.0,  0.0,  0.0)
        self.mp_origin.CreateNewNode(2,  0.0,  0.0,  0.0)
        self.mp_origin.CreateNewNode(3,  5.0,  0.0,  0.0)
        self.mp_origin.CreateNewNode(4,  8.0,  0.0,  0.0)

        self.mp_destination.CreateNewNode(1, -5.0,  0.0,  0.0)
        self.mp_destination.CreateNewNode(2, -3.0,  0.0,  0.0)
        self.mp_destination.CreateNewNode(3, -1.0,  0.0,  0.0)
        self.mp_destination.CreateNewNode(4,  3.0,  0.0,  0.0)
        self.mp_destination.CreateNewNode(5,  6.0,  0.0,  0.0)


    def _create_elements(self):
        element_name = "Element2D2N"
        # This seems to create properties on the fly
        props =self. mp_origin.GetProperties()[1]
        self.mp_origin.CreateNewElement(element_name, 1, [1,2], props)
        self.mp_origin.CreateNewElement(element_name, 2, [2,3], props)
        self.mp_origin.CreateNewElement(element_name, 3, [3,4], props)

        props =self. mp_destination.GetProperties()[1]
        self.mp_destination.CreateNewElement(element_name, 1, [1,2], props)
        self.mp_destination.CreateNewElement(element_name, 2, [2,3], props)
        self.mp_destination.CreateNewElement(element_name, 3, [3,4], props)
        self.mp_destination.CreateNewElement(element_name, 4, [4,5], props)


    def _set_values_origin(self):
        value = 0
        for node in self.mp_origin.Nodes:
            node.SetSolutionStepValue(KM.PRESSURE, value+0.2)
            node.SetSolutionStepValue(KM.FORCE, [value, value+0.1, value-0.3])
            value += 1


    def _set_values_destination(self):
        value = 0
        for node in self.mp_destination.Nodes:
            node.SetSolutionStepValue(KM.TEMPERATURE, value-0.3)
            node.SetSolutionStepValue(KM.VELOCITY, [value, value-0.1, value+0.4])
            value += 1


    def _set_values_mp_const(self, mp, variable, value):
        for node in mp.Nodes:
            node.SetSolutionStepValue(variable, value)


    def _create_mapper(self, mapper_name):
        mapper_settings = KM.Parameters("""{
            "mapper_type" : \"""" + mapper_name + """\"
         }""")

        self.mapper = KM.MapperFactory.CreateMapper(self.mp_origin,
                                                               self.mp_destination,
                                                               mapper_settings)

    def _check_results_scalar(self, mp, results, variable):
        if len(results) != mp.NumberOfNodes():
            raise RuntimeError("Number of results does not match number of Nodes!")
        for index, node in enumerate(mp.Nodes):
            self.assertAlmostEqual(node.GetSolutionStepValue(variable), results[index], 10)


    def _check_results_vector(self, mp, results, variable):
        if len(results) != mp.NumberOfNodes():
            raise RuntimeError("Number of results does not match number of Nodes!")
        for index, node in enumerate(mp.Nodes):
            self.assertAlmostEqual(node.GetSolutionStepValue(variable)[0], results[index][0], 10)
            self.assertAlmostEqual(node.GetSolutionStepValue(variable)[1], results[index][1], 10)
            self.assertAlmostEqual(node.GetSolutionStepValue(variable)[2], results[index][2], 10)


    def _check_results_scalar_const(self, mp, value, variable):
        for node in mp.Nodes:
            self.assertAlmostEqual(node.GetSolutionStepValue(variable), value)

    def _check_results_vector_const(self, mp, value, variable):
        for node in mp.Nodes:
            self.assertAlmostEqual(node.GetSolutionStepValue(variable)[0], value[0])
            self.assertAlmostEqual(node.GetSolutionStepValue(variable)[1], value[1])
            self.assertAlmostEqual(node.GetSolutionStepValue(variable)[2], value[2])

    def _execute_constant_value_test(self):
        # Check mapping of a constant field and the basic functionalities
        ### Map ###
        # Scalar Mapping
        mapping_value = 1.33
        self._set_values_mp_const(self.mp_origin, KM.PRESSURE, mapping_value)

        self.mapper.Map(KM.PRESSURE, KM.TEMPERATURE)
        self._check_results_scalar_const(self.mp_destination, mapping_value, KM.TEMPERATURE)

        self.mapper.UpdateInterface()

        self.mapper.Map(KM.PRESSURE, KM.TEMPERATURE, KM.Mapper.ADD_VALUES)
        self._check_results_scalar_const(self.mp_destination, 2*mapping_value, KM.TEMPERATURE)

        self.mapper.Map(KM.PRESSURE, KM.TEMPERATURE, KM.Mapper.ADD_VALUES | KM.Mapper.SWAP_SIGN)
        self._check_results_scalar_const(self.mp_destination, mapping_value, KM.TEMPERATURE)

        # Vector Mapping
        mapping_value = [1.443, -5.874, 7.99]
        self._set_values_mp_const(self.mp_origin, KM.FORCE, mapping_value)

        self.mapper.Map(KM.FORCE, KM.VELOCITY)
        self._check_results_vector_const(self.mp_destination, mapping_value, KM.VELOCITY)

        self.mapper.Map(KM.FORCE, KM.VELOCITY, KM.Mapper.ADD_VALUES)
        self._check_results_vector_const(self.mp_destination, [2*x for x in mapping_value], KM.VELOCITY)

        self.mapper.Map(KM.FORCE, KM.VELOCITY, KM.Mapper.ADD_VALUES | KM.Mapper.SWAP_SIGN)
        self._check_results_vector_const(self.mp_destination, mapping_value, KM.VELOCITY)

        ### InverseMap ###
        # Scalar Mapping
        mapping_value = -71.33
        self._set_values_mp_const(self.mp_destination, KM.TEMPERATURE, mapping_value)

        self.mapper.InverseMap(KM.PRESSURE, KM.TEMPERATURE)
        self._check_results_scalar_const(self.mp_origin, mapping_value, KM.PRESSURE)

        self.mapper.InverseMap(KM.PRESSURE, KM.TEMPERATURE, KM.Mapper.ADD_VALUES)
        self._check_results_scalar_const(self.mp_origin, 2*mapping_value, KM.PRESSURE)

        self.mapper.InverseMap(KM.PRESSURE, KM.TEMPERATURE, KM.Mapper.ADD_VALUES | KM.Mapper.SWAP_SIGN)
        self._check_results_scalar_const(self.mp_origin, mapping_value, KM.PRESSURE)

        # Vector Mapping
        mapping_value = [-5.443, 44.874, -7.9779]
        self._set_values_mp_const(self.mp_destination, KM.VELOCITY, mapping_value)

        self.mapper.InverseMap(KM.FORCE, KM.VELOCITY)
        self._check_results_vector_const(self.mp_origin, mapping_value, KM.FORCE)

        self.mapper.UpdateInterface(KM.Mapper.REMESHED)

        self.mapper.InverseMap(KM.FORCE, KM.VELOCITY, KM.Mapper.ADD_VALUES)
        self._check_results_vector_const(self.mp_origin, [2*x for x in mapping_value], KM.FORCE)

        self.mapper.InverseMap(KM.FORCE, KM.VELOCITY, KM.Mapper.ADD_VALUES | KM.Mapper.SWAP_SIGN)
        self._check_results_vector_const(self.mp_origin, mapping_value, KM.FORCE)


    def _execute_non_constant_value_test(self, results, use_transpose=False):
        # Check mapping of a non-constant field

        if use_transpose:
            mapper_flag = KM.Mapper.USE_TRANSPOSE
        else:
            mapper_flag=KM.Flags()

        ### Map ###
        # Scalar Mapping
        self._set_values_origin()

        self.mapper.Map(KM.PRESSURE, KM.TEMPERATURE, mapper_flag)
        self._check_results_scalar(self.mp_destination, results[0], KM.TEMPERATURE)
        if use_transpose:
            self.__CheckValuesSum(self.mp_origin, self.mp_destination, KM.PRESSURE, KM.TEMPERATURE, True)

        # Vector Mapping
        self.mapper.Map(KM.FORCE, KM.VELOCITY, mapper_flag)
        self._check_results_vector(self.mp_destination, results[1], KM.VELOCITY)
        if use_transpose:
            self.__CheckValuesSum(self.mp_origin, self.mp_destination, KM.FORCE, KM.VELOCITY, True)

        ### InverseMap ###
        # Scalar Mapping
        self._set_values_destination()
        self.mapper.InverseMap(KM.PRESSURE, KM.TEMPERATURE, mapper_flag)
        self._check_results_scalar(self.mp_origin, results[2], KM.PRESSURE)
        if use_transpose:
            self.__CheckValuesSum(self.mp_origin, self.mp_destination, KM.PRESSURE, KM.TEMPERATURE, True)

        # Vector Mapping
        self.mapper.InverseMap(KM.FORCE, KM.VELOCITY, mapper_flag)
        self._check_results_vector(self.mp_origin, results[3], KM.FORCE)
        if use_transpose:
            self.__CheckValuesSum(self.mp_origin, self.mp_destination, KM.FORCE, KM.VELOCITY, True)

    def test_nearest_neighbor_mapper(self):
        mapper_name = "nearest_neighbor"

        map_results_scalar = [0.2, 0.2, 1.2, 2.2, 2.2]
        map_results_vector     = [[0.0,0.1,-0.3], [0.0,0.1,-0.3], [1.0,1.1,0.7], [2.0,2.1,1.7], [2.0,2.1,1.7]]
        inverse_map_results_scalar     = [-0.3, 1.7, 3.7, 3.7]
        inverse_map_results_vector     = [[0.0,-0.1,0.4], [2.0,1.9,2.4], [4.0,3.9,4.4], [4.0,3.9,4.4]]

        map_results_scalar_conservative = [0.2, 0.0, 1.2, 0.0, 5.4]
        map_results_vector_conservative     = [[0.0,0.1,-0.3], [0.0,0.0,0.0], [1.0,1.1,0.7], [0.0,0.0,0.0], [5.0,5.2,4.4]]
        inverse_map_results_scalar_conservative = [0.4, 1.7, 6.4, 0.0]
        inverse_map_results_vector_conservative = [[1.0,0.8,1.8], [2.0,1.9,2.4], [7.0,6.8,7.8], [0.0,0.0,0.0]]

        results = [map_results_scalar, map_results_vector]
        results.extend([inverse_map_results_scalar, inverse_map_results_vector])

        results_conservative = [map_results_scalar_conservative, map_results_vector_conservative]
        results_conservative.extend([inverse_map_results_scalar_conservative, inverse_map_results_vector_conservative])

        self._create_mapper(mapper_name)

        self._execute_constant_value_test()

        self._execute_non_constant_value_test(results)

        # Test Mapping with transpose
        self._execute_non_constant_value_test(results_conservative, True)


    def test_nearest_element_mapper(self):
        mapper_name = "nearest_element"

        map_results_scalar = [0.2, 0.6, 1.0, 1.8, 7.6/3]
        map_results_vector     = [[0.0,0.1,-0.3], [0.4,0.5,0.1], [0.8,0.9,0.5], [1.6,1.7,1.3], [7/3,7.3/3,6.1/3]]
        inverse_map_results_scalar     = [-0.3, 1.95, 10.1/3, 3.7]
        inverse_map_results_vector     = [[0.0,-0.1,0.4], [2.25,2.15,2.65], [11/3,10.7/3,12.2/3], [4.0,3.9,4.4]]

        map_results_scalar_conservative = [0.2, 0.0, 0.9, 1.033333333333337, 14/3]
        map_results_vector_conservative     = [[0.0,0.1,-0.3], [0.0,0.0,0.0], [0.75,0.825,0.525], [11/12,0.975,0.741666666667], [13/3,4.5,23/6]]
        inverse_map_results_scalar_conservative     = [0.46, 2.72, 4.0866666666666, 37/30]
        inverse_map_results_vector_conservative     = [[1.0,0.82,1.72], [3.2,3.04,3.84], [67/15,4.34,373/75], [4/3,1.3,22/15]]

        results = [map_results_scalar, map_results_vector]
        results.extend([inverse_map_results_scalar, inverse_map_results_vector])

        results_conservative = [map_results_scalar_conservative, map_results_vector_conservative]
        results_conservative.extend([inverse_map_results_scalar_conservative, inverse_map_results_vector_conservative])

        self._create_mapper(mapper_name)

        self._execute_constant_value_test()

        self._execute_non_constant_value_test(results)

        # # Test conservative Mapping
        self._execute_non_constant_value_test(results_conservative, True)

    def __CheckValuesSum(self, mp1, mp2, var1, var2, value_is_historical=True):
        var_type = KM.KratosGlobals.GetVariableType(var1.Name())
        if var_type != KM.KratosGlobals.GetVariableType(var2.Name()):
            raise TypeError("Variable types-mismatch!")

        if value_is_historical:
            if var_type == "Double":
                val_1 = KM.VariableUtils().SumHistoricalNodeScalarVariable(var1, mp1, 0)
                val_2 = KM.VariableUtils().SumHistoricalNodeScalarVariable(var2, mp2, 0)
                self.assertAlmostEqual(val_1, val_2)
            else:
                val_1 = KM.VariableUtils().SumHistoricalNodeVectorVariable(var1, mp1, 0)
                val_2 = KM.VariableUtils().SumHistoricalNodeVectorVariable(var2, mp2, 0)
                self.assertAlmostEqual(val_1[0], val_2[0])
                self.assertAlmostEqual(val_1[1], val_2[1])
                self.assertAlmostEqual(val_1[2], val_2[2])
        else:
            if var_type == "Double":
                val_1 = KM.VariableUtils().SumNonHistoricalNodeScalarVariable(var1, mp1)
                val_2 = KM.VariableUtils().SumNonHistoricalNodeScalarVariable(var2, mp2)
                self.assertAlmostEqual(val_1, val_2)
            else:
                val_1 = KM.VariableUtils().SumNonHistoricalNodeVectorVariable(var1, mp1)
                val_2 = KM.VariableUtils().SumNonHistoricalNodeVectorVariable(var2, mp2)
                self.assertAlmostEqual(val_1[0], val_2[0])
                self.assertAlmostEqual(val_1[1], val_2[1])
                self.assertAlmostEqual(val_1[2], val_2[2])


if __name__ == '__main__':
    KratosUnittest.main()
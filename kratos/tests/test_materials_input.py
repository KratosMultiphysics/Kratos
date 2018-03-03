from __future__ import print_function, absolute_import, division

import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics
import os

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

class TestMaterialsInput(KratosUnittest.TestCase):

    def test_input(self):
        try:
            import KratosMultiphysics.FluidDynamicsApplication as KratosFluid
        except:
            self.skipTest("KratosMultiphysics.FluidDynamicsApplication is not available")

        try:
            import KratosMultiphysics.StructuralMechanicsApplication
        except:
            self.skipTest("KratosMultiphysics.StructuralMechanicsApplication is not available")

        model_part = KratosMultiphysics.ModelPart("Main")
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VISCOSITY)
        model_part_io = KratosMultiphysics.ModelPartIO(GetFilePath("test_model_part_io_read")) #reusing the file that is already in the directory
        model_part_io.ReadModelPart(model_part)

        #define a Model TODO: replace to use the real Model once available
        Model = {
            "Main" : model_part,
            "Inlets" : model_part.GetSubModelPart("Inlets"),
            "Inlet1" : model_part.GetSubModelPart("Inlets").GetSubModelPart("Inlet1"),
            "Inlet2" : model_part.GetSubModelPart("Inlets").GetSubModelPart("Inlet2"),
            "Outlet" : model_part.GetSubModelPart("Outlet")
            }

        test_settings = KratosMultiphysics.Parameters("""
            {
                    "Parameters": {
                            "materials_filename": "materials.json"
                    }
            }
            """)

        #assign the real path
        test_settings["Parameters"]["materials_filename"].SetString(GetFilePath("materials.json"))

        import read_materials_process
        read_materials_process.Factory(test_settings,Model)

        #test if the element properties are assigned correctly to the elements and conditions
        for elem in Model["Inlets"].Elements:
            self.assertTrue(elem.Properties.Id == 1)
        for cond in Model["Inlets"].Conditions:
            self.assertTrue(cond.Properties.Id == 1)
        for elem in Model["Outlet"].Elements:
            self.assertTrue(elem.Properties.Id == 2)
        for cond in Model["Outlet"].Conditions:
            self.assertTrue(cond.Properties.Id == 2)

        #test that the properties are read correctly
        self.assertTrue(model_part.Properties[1].GetValue(KratosMultiphysics.YOUNG_MODULUS) == 200.0)
        self.assertTrue(model_part.Properties[1].GetValue(KratosMultiphysics.POISSON_RATIO) == 0.3)
        self.assertTrue(model_part.Properties[1].GetValue(KratosMultiphysics.YIELD_STRESS) == 400.0)
        self.assertTrue(model_part.Properties[1].GetValue(KratosFluid.SUBSCALE_PRESSURE) == 0.1)
        self.assertTrue(model_part.Properties[1].GetValue(KratosFluid.VORTICITY_MAGNITUDE) == -5.888)

        self.assertTrue(model_part.Properties[2].GetValue(KratosMultiphysics.YOUNG_MODULUS) == 100.0)
        self.assertTrue(model_part.Properties[2].GetValue(KratosMultiphysics.POISSON_RATIO) == 0.1)
        self.assertTrue(model_part.Properties[2].GetValue(KratosMultiphysics.YIELD_STRESS) == 800.0)
        self.assertTrue(model_part.Properties[2].GetValue(KratosMultiphysics.HTC) == 0.3)
        self.assertTrue(model_part.Properties[2].GetValue(KratosMultiphysics.TIME_STEPS) == 159) # int
        self.assertTrue(model_part.Properties[2].GetValue(KratosMultiphysics.UPDATE_SENSITIVITIES) == True) # bool
        self.assertTrue(model_part.Properties[2].GetValue(KratosMultiphysics.IDENTIFIER) == "MyTestString") # std::string

        mat_vector = model_part.Properties[2].GetValue(KratosMultiphysics.CAUCHY_STRESS_VECTOR)
        self.assertAlmostEqual(mat_vector[0],1.5)
        self.assertAlmostEqual(mat_vector[1],0.3)
        self.assertAlmostEqual(mat_vector[2],-2.58)

        mat_matrix = model_part.Properties[2].GetValue(KratosMultiphysics.LOCAL_INERTIA_TENSOR)
        self.assertAlmostEqual(mat_matrix[0,0],1.27)
        self.assertAlmostEqual(mat_matrix[0,1],-22.5)
        self.assertAlmostEqual(mat_matrix[1,0],2.01)
        self.assertAlmostEqual(mat_matrix[1,1],0.257)

        table = model_part.Properties[2].GetTable(KratosMultiphysics.TEMPERATURE, KratosMultiphysics.YOUNG_MODULUS)
        self.assertAlmostEqual(table.GetValue(1.5),11.0)
        self.assertAlmostEqual(table.GetNearestValue(1.1),10.0)
        self.assertAlmostEqual(table.GetDerivative(1.2),2.0)

    def test_input_interpolative(self):
        try:
            import KratosMultiphysics.FluidDynamicsApplication as KratosFluid
        except:
            self.skipTest("KratosMultiphysics.FluidDynamicsApplication is not available")

        try:
            import KratosMultiphysics.StructuralMechanicsApplication
        except:
            self.skipTest("KratosMultiphysics.StructuralMechanicsApplication is not available")

        model_part = KratosMultiphysics.ModelPart("Main")
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.AUX_INDEX)
        model_part_io = KratosMultiphysics.ModelPartIO(GetFilePath("test_processes")) #reusing the file that is already in the directory
        model_part_io.ReadModelPart(model_part)

        #define a Model TODO: replace to use the real Model once available
        Model = {
            "Main" : model_part,
            "Main_domain" : model_part.GetSubModelPart("Main_domain"),
            "Left_side" : model_part.GetSubModelPart("Left_side")
            }

        test_settings = KratosMultiphysics.Parameters("""
            {
                    "Parameters": {
                            "materials_filename": "materials_interpolative.json"
                    }
            }
            """)


        #assign the real path
        test_settings["Parameters"]["materials_filename"].SetString(GetFilePath("materials_interpolative.json"))

        # Populate the Entities with values, usually these are coming from mdpa
        for node in model_part.Nodes:
            node.SetSolutionStepValue(KratosMultiphysics.AUX_INDEX, self._get_value_for_entity(node))
        for element in model_part.Elements:
            element.SetValue(KratosMultiphysics.TEMPERATURE, self._get_value_for_entity(element))

        import read_materials_process
        read_materials_process.Factory(test_settings, Model)

        self.assertEqual(model_part.NumberOfTables(), 4)
        table_1 = model_part.GetTable(1)
        table_2 = model_part.GetTable(2)
        table_3 = model_part.GetTable(3)
        table_4 = model_part.GetTable(4)
        self.assertAlmostEqual(table_1.GetValue(100),300.0)
        self.assertAlmostEqual(table_2.GetValue(150),200.0)
        self.assertAlmostEqual(table_3.GetValue(110),23.75)
        self.assertAlmostEqual(table_4.GetValue(75),121.25)

        self.assertEqual(model_part.NumberOfProperties(), 23) # 2 props are in the mdpa already

        for elem in model_part.GetSubModelPart("Main_domain").Elements:
            self.assertAlmostEqual(elem.Properties.GetValue(KratosMultiphysics.POISSON_RATIO), 0.39)
            self.assertAlmostEqual(elem.Properties.GetValue(KratosMultiphysics.YOUNG_MODULUS),
                                    table_1.GetValue(elem.GetValue(KratosMultiphysics.TEMPERATURE)))


            local_inertia_tensor = elem.Properties.GetValue(KratosMultiphysics.LOCAL_INERTIA_TENSOR)
            local_inertia_tensor_expected = KratosMultiphysics.Matrix(2,2)
            local_inertia_tensor_expected[0,0] = 1.27
            local_inertia_tensor_expected[0,1] = table_3.GetValue(elem.GetValue(KratosMultiphysics.TEMPERATURE))
            local_inertia_tensor_expected[1,0] = table_2.GetValue(self._get_value_for_entity(elem))
            local_inertia_tensor_expected[1,1] = 0.257
            for i in range(2):
                for j in range(2):
                    self.assertAlmostEqual(local_inertia_tensor[i,j], local_inertia_tensor_expected[i,j])

        for cond in model_part.GetSubModelPart("Left_side").Conditions:
            self.assertAlmostEqual(cond.Properties.GetValue(KratosMultiphysics.POISSON_RATIO), 0.55)

            cauchy_stress_vector = cond.Properties.GetValue(KratosMultiphysics.CAUCHY_STRESS_VECTOR)
            cauchy_stress_vector_expected = KratosMultiphysics.Vector(3)
            cauchy_stress_vector_expected[0] = table_4.GetValue(self._get_value_from_nodes(cond))
            cauchy_stress_vector_expected[1] = 0.3
            cauchy_stress_vector_expected[2] = -2.58
            for i in range(3):
                self.assertAlmostEqual(cauchy_stress_vector[i], cauchy_stress_vector_expected[i])


            # TODO test also the tables => should be the same as in the ModelPart

    def _get_value_for_entity(self, entity):
        # assign a value randomly
        return (20.0*entity.Id)%200

    def _get_value_from_nodes(self, entity):
        cond_nodes = entity.GetNodes()
        num_nodes = len(cond_nodes)
        value = 0.0

        for node in cond_nodes:
            value += self._get_value_for_entity(node)

        return value / num_nodes


if __name__ == '__main__':
    KratosUnittest.main()

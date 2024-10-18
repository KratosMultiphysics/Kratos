# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest

import os

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

class TestVariableUtils(KratosUnittest.TestCase):

    def test_copy_model_part_nodal_var(self):
        current_model = KratosMultiphysics.Model()

        ##set the origin model part
        origin_model_part = current_model.CreateModelPart("OriginModelPart")
        origin_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VISCOSITY)
        origin_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        origin_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        origin_model_part.SetBufferSize(2)
        model_part_io = KratosMultiphysics.ModelPartIO(GetFilePath("auxiliar_files_for_python_unittest/mdpa_files/test_model_part_io_read"))
        model_part_io.ReadModelPart(origin_model_part)

        ##set the destination model part
        destination_model_part = current_model.CreateModelPart("DestinationModelPart")
        destination_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VISCOSITY)
        destination_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        destination_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        destination_model_part.SetBufferSize(2)
        model_part_io = KratosMultiphysics.ModelPartIO(GetFilePath("auxiliar_files_for_python_unittest/mdpa_files/test_model_part_io_read"))
        model_part_io.ReadModelPart(destination_model_part)

        ##set the values in the origin model part
        for node in origin_model_part.Nodes:
            node.SetSolutionStepValue(KratosMultiphysics.VISCOSITY, 0, node.X + node.Y)
            node.SetSolutionStepValue(KratosMultiphysics.VISCOSITY, 1, 2.0 * node.X + 3.0 * node.Y)
            node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT, 0, [node.X ** 2, 0.0, 0.0])
            node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT, 1, [node.X, node.Y, node.Z])

        # ##copy the values to the destination model part
        KratosMultiphysics.VariableUtils().CopyModelPartNodalVar(KratosMultiphysics.VISCOSITY, origin_model_part, destination_model_part, 0)
        KratosMultiphysics.VariableUtils().CopyModelPartNodalVar(KratosMultiphysics.VISCOSITY, origin_model_part, destination_model_part, 1)
        KratosMultiphysics.VariableUtils().CopyModelPartNodalVar(KratosMultiphysics.DISPLACEMENT, origin_model_part, destination_model_part, 0)
        KratosMultiphysics.VariableUtils().CopyModelPartNodalVar(KratosMultiphysics.DISPLACEMENT, origin_model_part, destination_model_part, 1)

        ##check the copied values
        for node in destination_model_part.Nodes:
            self.assertEqual(node.GetSolutionStepValue(KratosMultiphysics.VISCOSITY, 0), node.X + node.Y)
            self.assertEqual(node.GetSolutionStepValue(KratosMultiphysics.VISCOSITY, 1), 2.0 * node.X + 3.0 * node.Y)
            self.assertEqual(node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X, 0), node.X ** 2)
            self.assertEqual(node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X, 1), node.X)

        ##copy the values to the destination model part in different buffers
        KratosMultiphysics.VariableUtils().CopyModelPartNodalVar(KratosMultiphysics.VISCOSITY, KratosMultiphysics.VISCOSITY ,origin_model_part, destination_model_part, 0, 1)
        KratosMultiphysics.VariableUtils().CopyModelPartNodalVar(KratosMultiphysics.DISPLACEMENT, KratosMultiphysics.DISPLACEMENT,origin_model_part, destination_model_part, 0, 1)

        ##check the copied values
        for node in destination_model_part.Nodes:
            self.assertEqual(node.GetSolutionStepValue(KratosMultiphysics.VISCOSITY, 1), node.X + node.Y)
            self.assertVectorAlmostEqual(node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT,1),[node.X ** 2, 0.0, 0.0])


    def test_copy_model_part_nodal_var_to_non_historical_var(self):
        ##set the origin model part
        current_model = KratosMultiphysics.Model()
        origin_model_part = current_model.CreateModelPart("OriginModelPart")
        origin_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VISCOSITY)
        origin_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        origin_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        origin_model_part.SetBufferSize(2)
        model_part_io = KratosMultiphysics.ModelPartIO(GetFilePath("auxiliar_files_for_python_unittest/mdpa_files/test_model_part_io_read"))
        model_part_io.ReadModelPart(origin_model_part)

        ##set the destination model part
        destination_model_part = current_model.CreateModelPart("DestinationModelPart")
        destination_model_part.SetBufferSize(2)
        destination_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VISCOSITY)
        destination_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        destination_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        model_part_io = KratosMultiphysics.ModelPartIO(GetFilePath("auxiliar_files_for_python_unittest/mdpa_files/test_model_part_io_read"))
        model_part_io.ReadModelPart(destination_model_part)

        ##set the values in the origin model part
        for node in origin_model_part.Nodes:
            node.SetSolutionStepValue(KratosMultiphysics.VISCOSITY, 0, node.X + node.Y)
            node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT, 1, [node.X, node.Y, node.Z])

        ##  initialize the containers in destination model part (otherwise the operation is not threadsafe!)
        for node in destination_model_part.Nodes:
            node.SetValue(KratosMultiphysics.VISCOSITY, 0)
            node.SetValue(KratosMultiphysics.DISPLACEMENT, KratosMultiphysics.Array3())
            node.SetValue(KratosMultiphysics.VELOCITY, KratosMultiphysics.Array3())

        ##copy the values to the destination model part
        KratosMultiphysics.VariableUtils().CopyModelPartNodalVarToNonHistoricalVar(KratosMultiphysics.VISCOSITY, origin_model_part, destination_model_part, 0)
        KratosMultiphysics.VariableUtils().CopyModelPartNodalVarToNonHistoricalVar(KratosMultiphysics.DISPLACEMENT, origin_model_part, destination_model_part, 1)
        KratosMultiphysics.VariableUtils().CopyModelPartNodalVarToNonHistoricalVar(KratosMultiphysics.DISPLACEMENT, KratosMultiphysics.VELOCITY, origin_model_part, destination_model_part, 1)

        ##check the copied values
        for node in destination_model_part.Nodes:
            self.assertEqual(node.GetValue(KratosMultiphysics.VISCOSITY), node.X + node.Y)
            self.assertEqual(node.GetValue(KratosMultiphysics.VELOCITY_X), node.X)
            self.assertEqual(node.GetValue(KratosMultiphysics.DISPLACEMENT_X), node.X)

    def test_copy_model_part_flagged_nodal_historical_var_to_historical_var(self):
        def test_method(m_orig, m_dest, v_orig, v_dest, read_buffer, write_buffer):
            KratosMultiphysics.VariableUtils().CopyModelPartFlaggedNodalHistoricalVarToHistoricalVar(
                v_orig, v_dest, m_orig, m_dest, KratosMultiphysics.INLET, True, read_buffer, write_buffer)

        def test_case(read_buffer, write_buffer):
            self.__test_copy_model_part_flagged_methods(
                lambda m_orig, m_dest, v_orig, v_dest: test_method(m_orig, m_dest, v_orig, v_dest, read_buffer, write_buffer),
                lambda model_part: model_part.Nodes,
                lambda node, variable, value: node.SetSolutionStepValue(variable, read_buffer, value),
                lambda node, variable: node.GetSolutionStepValue(variable, read_buffer),
                lambda node, variable, value: node.SetSolutionStepValue(variable, write_buffer, value),
                lambda node, variable: node.GetSolutionStepValue(variable, write_buffer)
            )

        test_case(0, 0)
        test_case(0, 1)
        test_case(1, 0)
        test_case(1, 1)

    def test_copy_model_part_flagged_nodal_historical_var_to_non_historical_var(self):
        def test_method(m_orig, m_dest, v_orig, v_dest, read_buffer):
            KratosMultiphysics.VariableUtils().CopyModelPartFlaggedNodalHistoricalVarToNonHistoricalVar(
                v_orig, v_dest, m_orig, m_dest, KratosMultiphysics.INLET, True, read_buffer)
        def test_case(read_buffer):
            self.__test_copy_model_part_flagged_methods(
                lambda m_orig, m_dest, v_orig, v_dest: test_method(m_orig, m_dest, v_orig, v_dest, read_buffer),
                lambda model_part: model_part.Nodes,
                lambda node, variable, value: node.SetSolutionStepValue(variable, read_buffer, value),
                lambda node, variable: node.GetSolutionStepValue(variable, read_buffer),
                lambda node, variable, value: node.SetValue(variable, value),
                lambda node, variable: node.GetValue(variable)
            )

        test_case(0)
        test_case(1)

    def test_copy_model_part_flagged_nodal_non_historical_var_to_historical_var(self):
        def test_method(m_orig, m_dest, v_orig, v_dest, write_buffer):
            KratosMultiphysics.VariableUtils().CopyModelPartFlaggedNodalNonHistoricalVarToHistoricalVar(
                v_orig, v_dest, m_orig, m_dest, KratosMultiphysics.INLET, True, write_buffer)
        def test_case(write_buffer):
            self.__test_copy_model_part_flagged_methods(
                lambda m_orig, m_dest, v_orig, v_dest: test_method(m_orig, m_dest, v_orig, v_dest, write_buffer),
                lambda model_part: model_part.Nodes,
                lambda node, variable, value: node.SetValue(variable, value),
                lambda node, variable: node.GetValue(variable),
                lambda node, variable, value: node.SetSolutionStepValue(variable, write_buffer, value),
                lambda node, variable: node.GetSolutionStepValue(variable, write_buffer)
            )

        test_case(0)
        test_case(1)

    def test_copy_model_part_flagged_nodal_non_historical_var_to_non_historical_var(self):
        self.__test_copy_model_part_flagged_methods(
            lambda m_orig, m_dest, v_orig, v_dest: KratosMultiphysics.VariableUtils().CopyModelPartFlaggedNodalNonHistoricalVarToNonHistoricalVar(v_orig, v_dest, m_orig, m_dest, KratosMultiphysics.INLET, True),
            lambda model_part: model_part.Nodes,
            lambda node, variable, value: node.SetValue(variable, value),
            lambda node, variable: node.GetValue(variable),
            lambda node, variable, value: node.SetValue(variable, value),
            lambda node, variable: node.GetValue(variable)
        )

    def test_copy_model_part_flagged_element_var(self):
        self.__test_copy_model_part_flagged_methods(
            lambda m_orig, m_dest, v_orig, v_dest: KratosMultiphysics.VariableUtils().CopyModelPartFlaggedElementVar(v_orig, v_dest, m_orig, m_dest, KratosMultiphysics.INLET, True),
            lambda model_part: model_part.Elements,
            lambda element, variable, value: element.SetValue(variable, value),
            lambda element, variable: element.GetValue(variable),
            lambda element, variable, value: element.SetValue(variable, value),
            lambda element, variable: element.GetValue(variable)
        )

    def test_copy_model_part_flagged_condition_var(self):
        self.__test_copy_model_part_flagged_methods(
            lambda m_orig, m_dest, v_orig, v_dest: KratosMultiphysics.VariableUtils().CopyModelPartFlaggedConditionVar(v_orig, v_dest, m_orig, m_dest, KratosMultiphysics.INLET, True),
            lambda model_part: model_part.Conditions,
            lambda condition, variable, value: condition.SetValue(variable, value),
            lambda condition, variable: condition.GetValue(variable),
            lambda condition, variable, value: condition.SetValue(variable, value),
            lambda condition, variable: condition.GetValue(variable)
        )

    def test_copy_model_part_elemental_var(self):
        current_model = KratosMultiphysics.Model()

        ##set the origin model part
        origin_model_part = current_model.CreateModelPart("OriginModelPart")
        origin_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VISCOSITY)
        origin_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        origin_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        model_part_io = KratosMultiphysics.ModelPartIO(GetFilePath("auxiliar_files_for_python_unittest/mdpa_files/test_model_part_io_read"))
        model_part_io.ReadModelPart(origin_model_part)

        ##set the destination model part
        destination_model_part = current_model.CreateModelPart("DestinationModelPart")
        destination_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VISCOSITY)
        destination_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        destination_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        model_part_io = KratosMultiphysics.ModelPartIO(GetFilePath("auxiliar_files_for_python_unittest/mdpa_files/test_model_part_io_read"))
        model_part_io.ReadModelPart(destination_model_part)
        ##set the values in the destination model part
        for element in origin_model_part.Elements:
            element.SetValue(KratosMultiphysics.DENSITY, element.Id*100)
            element.SetValue(KratosMultiphysics.VOLUME_ACCELERATION, [element.Id*100, 0.0, 0.0])

        ##copy the values to the destination model part
        KratosMultiphysics.VariableUtils().CopyModelPartElementalVar(KratosMultiphysics.DENSITY, origin_model_part, destination_model_part)
        KratosMultiphysics.VariableUtils().CopyModelPartElementalVar(KratosMultiphysics.VOLUME_ACCELERATION, origin_model_part, destination_model_part)

        ##check the copied values
        for element in destination_model_part.Elements:
            self.assertEqual(element.GetValue(KratosMultiphysics.DENSITY), element.Id*100)
            self.assertEqual(element.GetValue(KratosMultiphysics.VOLUME_ACCELERATION)[0], element.Id*100)


        current_model = KratosMultiphysics.Model()

        ##set the model part
        model_part = current_model.CreateModelPart("Main")
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VISCOSITY)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PARTITION_INDEX)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ORIENTATION)
        model_part_io = KratosMultiphysics.ModelPartIO(GetFilePath("auxiliar_files_for_python_unittest/mdpa_files/test_model_part_io_read"))
        model_part_io.ReadModelPart(model_part)

        ##set the variable values
        viscosity = 0.1
        partition_index = 1
        velocity = KratosMultiphysics.Vector(3)
        velocity[0] = 2.0
        velocity[1] = 4.0
        velocity[2] = 8.0
        displacement = KratosMultiphysics.Vector(3)
        displacement[0] = 1.0
        displacement[1] = 2.0
        displacement[2] = 3.0

        orientation = KratosMultiphysics.Quaternion()
        orientation.X = 1.0
        orientation.Y = 2.0
        orientation.Z = 3.0
        orientation.W = 4.0

        KratosMultiphysics.VariableUtils().SetScalarVar(KratosMultiphysics.VISCOSITY, viscosity, model_part.Nodes)
        KratosMultiphysics.VariableUtils().SetVariable(KratosMultiphysics.VELOCITY_X, velocity[0], model_part.Nodes)
        KratosMultiphysics.VariableUtils().SetVariable(KratosMultiphysics.VELOCITY_Y, velocity[1], model_part.Nodes)
        KratosMultiphysics.VariableUtils().SetVariable(KratosMultiphysics.VELOCITY_Z, velocity[2], model_part.Nodes)
        KratosMultiphysics.VariableUtils().SetVariable(KratosMultiphysics.DISPLACEMENT, displacement, model_part.Nodes)
        KratosMultiphysics.VariableUtils().SetVariable(KratosMultiphysics.PARTITION_INDEX, partition_index, model_part.Nodes)
        KratosMultiphysics.VariableUtils().SetVariable(KratosMultiphysics.ORIENTATION, orientation, model_part.Nodes)

        ##verify the result
        for node in model_part.Nodes:
            self.assertEqual(node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X), 1.0)
            self.assertEqual(node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y), 2.0)
            self.assertEqual(node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Z), 3.0)
            self.assertEqual(node.GetSolutionStepValue(KratosMultiphysics.VISCOSITY), viscosity)
            self.assertEqual(node.GetSolutionStepValue(KratosMultiphysics.ORIENTATION).X, 1.0)
            self.assertEqual(node.GetSolutionStepValue(KratosMultiphysics.ORIENTATION).W, 4.0)

    def test_set_nonhistorical_variable(self):
        current_model = KratosMultiphysics.Model()

        ##set the model part
        model_part = current_model.CreateModelPart("Main")
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VISCOSITY)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ORIENTATION)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        model_part_io = KratosMultiphysics.ModelPartIO(GetFilePath("auxiliar_files_for_python_unittest/mdpa_files/test_model_part_io_read"))
        model_part_io.ReadModelPart(model_part)

        ##set the variable values
        viscosity = 0.1
        displacement = KratosMultiphysics.Vector(3)
        displacement[0] = 1.0
        displacement[1] = 2.0
        displacement[2] = 3.0

        orientation = KratosMultiphysics.Quaternion()
        orientation.X = 1.0
        orientation.Y = 2.0
        orientation.Z = 3.0
        orientation.W = 4.0

        # First for nodes
        KratosMultiphysics.VariableUtils().SetNonHistoricalVariable(KratosMultiphysics.VISCOSITY, viscosity, model_part.Nodes)
        KratosMultiphysics.VariableUtils().SetNonHistoricalVariable(KratosMultiphysics.DISPLACEMENT, displacement, model_part.Nodes)
        KratosMultiphysics.VariableUtils().SetNonHistoricalVariable(KratosMultiphysics.ORIENTATION, orientation, model_part.Nodes)


        ##verify the result
        for node in model_part.Nodes:
            self.assertEqual(node.GetValue(KratosMultiphysics.DISPLACEMENT_X), 1.0)
            self.assertEqual(node.GetValue(KratosMultiphysics.DISPLACEMENT_Y), 2.0)
            self.assertEqual(node.GetValue(KratosMultiphysics.DISPLACEMENT_Z), 3.0)
            self.assertEqual(node.GetValue(KratosMultiphysics.VISCOSITY), viscosity)
            self.assertEqual(node.GetValue(KratosMultiphysics.ORIENTATION).X, 1.0)
            self.assertEqual(node.GetValue(KratosMultiphysics.ORIENTATION).W, 4.0)

        # Now for conditions (it will work for elements too)
        KratosMultiphysics.VariableUtils().SetNonHistoricalVariable(KratosMultiphysics.VISCOSITY, viscosity, model_part.Conditions)
        KratosMultiphysics.VariableUtils().SetNonHistoricalVariable(KratosMultiphysics.DISPLACEMENT, displacement, model_part.Conditions)

        ##verify the result
        for cond in model_part.Conditions:
            disp = cond.GetValue(KratosMultiphysics.DISPLACEMENT)
            self.assertEqual(disp[0], 1.0)
            self.assertEqual(disp[1], 2.0)
            self.assertEqual(disp[2], 3.0)
            self.assertEqual(cond.GetValue(KratosMultiphysics.VISCOSITY), viscosity)

    def test_clear_nonhistorical_values(self):
        current_model = KratosMultiphysics.Model()

        ##set the model part
        model_part = current_model.CreateModelPart("Main")
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VISCOSITY)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        model_part_io = KratosMultiphysics.ModelPartIO(GetFilePath("auxiliar_files_for_python_unittest/mdpa_files/test_model_part_io_read"))
        model_part_io.ReadModelPart(model_part)

        ##set the variable values
        viscosity = 0.1
        displacement = KratosMultiphysics.Vector(3)
        displacement[0] = 1.0
        displacement[1] = 2.0
        displacement[2] = 3.0

        # First for nodes
        KratosMultiphysics.VariableUtils().SetNonHistoricalVariable(KratosMultiphysics.VISCOSITY, viscosity, model_part.Nodes)
        KratosMultiphysics.VariableUtils().SetNonHistoricalVariable(KratosMultiphysics.DISPLACEMENT, displacement, model_part.Nodes)
        KratosMultiphysics.VariableUtils().ClearNonHistoricalData(model_part.Nodes)

        ##verify the result
        for node in model_part.Nodes:
            self.assertFalse(node.Has(KratosMultiphysics.DISPLACEMENT))
            self.assertFalse(node.Has(KratosMultiphysics.VISCOSITY))

        # Now for conditions (it will work for elements too)
        KratosMultiphysics.VariableUtils().SetNonHistoricalVariable(KratosMultiphysics.VISCOSITY, viscosity, model_part.Conditions)
        KratosMultiphysics.VariableUtils().SetNonHistoricalVariable(KratosMultiphysics.DISPLACEMENT, displacement, model_part.Conditions)
        KratosMultiphysics.VariableUtils().ClearNonHistoricalData(model_part.Conditions)

        ##verify the result
        for cond in model_part.Conditions:
            self.assertFalse(cond.Has(KratosMultiphysics.DISPLACEMENT))
            self.assertFalse(cond.Has(KratosMultiphysics.VISCOSITY))

    def test_set_flag(self):
        current_model = KratosMultiphysics.Model()

        ##set the model part
        model_part = current_model.CreateModelPart("Main")
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VISCOSITY)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        model_part_io = KratosMultiphysics.ModelPartIO(GetFilePath("auxiliar_files_for_python_unittest/mdpa_files/test_model_part_io_read"))
        model_part_io.ReadModelPart(model_part)

        KratosMultiphysics.VariableUtils().SetFlag(KratosMultiphysics.VISITED, True, model_part.Nodes)
        KratosMultiphysics.VariableUtils().SetFlag(KratosMultiphysics.VISITED, True, model_part.Conditions)
        KratosMultiphysics.VariableUtils().SetFlag(KratosMultiphysics.VISITED, True, model_part.Elements)

        KratosMultiphysics.VariableUtils().SetFlag(KratosMultiphysics.INLET, True, model_part.GetSubModelPart("Inlets").Nodes)
        KratosMultiphysics.VariableUtils().SetFlag(KratosMultiphysics.INLET, True, model_part.GetSubModelPart("Inlets").Conditions)
        KratosMultiphysics.VariableUtils().SetFlag(KratosMultiphysics.OUTLET, False, model_part.GetSubModelPart("Inlets").Nodes)
        KratosMultiphysics.VariableUtils().SetFlag(KratosMultiphysics.OUTLET, False, model_part.GetSubModelPart("Inlets").Conditions)

        ##verify the main modelpart flags set
        for node in model_part.Nodes:
            self.assertTrue(node.Is(KratosMultiphysics.VISITED))
        for condition in model_part.Conditions:
            self.assertTrue(condition.Is(KratosMultiphysics.VISITED))
        for element in model_part.Elements:
            self.assertTrue(element.Is(KratosMultiphysics.VISITED))
        ##verify the inlet submodelpart flag set
        for node in model_part.GetSubModelPart("Inlets").Nodes:
            self.assertTrue(node.Is(KratosMultiphysics.INLET))
            self.assertTrue(node.IsNot(KratosMultiphysics.OUTLET))
        for condition in model_part.GetSubModelPart("Inlets").Conditions:
            self.assertTrue(condition.Is(KratosMultiphysics.INLET))
            self.assertTrue(condition.IsNot(KratosMultiphysics.OUTLET))

        KratosMultiphysics.VariableUtils().FlipFlag(KratosMultiphysics.VISITED, model_part.Nodes)
        KratosMultiphysics.VariableUtils().FlipFlag(KratosMultiphysics.VISITED, model_part.Conditions)
        KratosMultiphysics.VariableUtils().FlipFlag(KratosMultiphysics.VISITED, model_part.Elements)

        KratosMultiphysics.VariableUtils().FlipFlag(KratosMultiphysics.INLET, model_part.GetSubModelPart("Inlets").Nodes)
        KratosMultiphysics.VariableUtils().FlipFlag(KratosMultiphysics.INLET, model_part.GetSubModelPart("Inlets").Conditions)
        KratosMultiphysics.VariableUtils().FlipFlag(KratosMultiphysics.OUTLET, model_part.GetSubModelPart("Inlets").Nodes)
        KratosMultiphysics.VariableUtils().FlipFlag(KratosMultiphysics.OUTLET, model_part.GetSubModelPart("Inlets").Conditions)

        ##verify the main modelpart flags set
        for node in model_part.Nodes:
            self.assertFalse(node.Is(KratosMultiphysics.VISITED))
        for condition in model_part.Conditions:
            self.assertFalse(condition.Is(KratosMultiphysics.VISITED))
        for element in model_part.Elements:
            self.assertFalse(element.Is(KratosMultiphysics.VISITED))
        ##verify the inlet submodelpart flag set
        for node in model_part.GetSubModelPart("Inlets").Nodes:
            self.assertFalse(node.Is(KratosMultiphysics.INLET))
            self.assertFalse(node.IsNot(KratosMultiphysics.OUTLET))
        for condition in model_part.GetSubModelPart("Inlets").Conditions:
            self.assertFalse(condition.Is(KratosMultiphysics.INLET))
            self.assertFalse(condition.IsNot(KratosMultiphysics.OUTLET))

        KratosMultiphysics.VariableUtils().ResetFlag(KratosMultiphysics.VISITED, model_part.Nodes)
        KratosMultiphysics.VariableUtils().ResetFlag(KratosMultiphysics.VISITED, model_part.Conditions)
        KratosMultiphysics.VariableUtils().ResetFlag(KratosMultiphysics.VISITED, model_part.Elements)

        KratosMultiphysics.VariableUtils().ResetFlag(KratosMultiphysics.INLET, model_part.GetSubModelPart("Inlets").Nodes)
        KratosMultiphysics.VariableUtils().ResetFlag(KratosMultiphysics.INLET, model_part.GetSubModelPart("Inlets").Conditions)
        KratosMultiphysics.VariableUtils().ResetFlag(KratosMultiphysics.OUTLET, model_part.GetSubModelPart("Inlets").Nodes)
        KratosMultiphysics.VariableUtils().ResetFlag(KratosMultiphysics.OUTLET, model_part.GetSubModelPart("Inlets").Conditions)

        ##verify the main modelpart flags unset
        for node in model_part.Nodes:
            self.assertFalse(node.IsDefined(KratosMultiphysics.VISITED))
        for condition in model_part.Conditions:
            self.assertFalse(condition.IsDefined(KratosMultiphysics.VISITED))
        for element in model_part.Elements:
            self.assertFalse(element.IsDefined(KratosMultiphysics.VISITED))
        ##verify the inlet submodelpart flag set
        for node in model_part.GetSubModelPart("Inlets").Nodes:
            self.assertFalse(node.IsDefined(KratosMultiphysics.INLET))
            self.assertFalse(node.IsDefined(KratosMultiphysics.OUTLET))
        for condition in model_part.GetSubModelPart("Inlets").Conditions:
            self.assertFalse(condition.IsDefined(KratosMultiphysics.INLET))
            self.assertFalse(condition.IsDefined(KratosMultiphysics.OUTLET))

    def test_copy_var(self):
        current_model = KratosMultiphysics.Model()

        ##set the model part
        model_part = current_model.CreateModelPart("Main")
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DENSITY)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VISCOSITY)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.FORCE)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)
        model_part_io = KratosMultiphysics.ModelPartIO(GetFilePath("auxiliar_files_for_python_unittest/mdpa_files/test_model_part_io_read"))
        model_part_io.ReadModelPart(model_part)

        ##set the variable values
        viscosity = 0.1
        displacement = KratosMultiphysics.Vector(3)
        displacement[0] = 1.3
        displacement[1] = 2.2
        displacement[2] = 3.1

        KratosMultiphysics.VariableUtils().SetVariable(KratosMultiphysics.VISCOSITY, viscosity, model_part.Nodes)
        KratosMultiphysics.VariableUtils().SetVariable(KratosMultiphysics.DISPLACEMENT, displacement, model_part.Nodes)
        KratosMultiphysics.VariableUtils().SetVariable(KratosMultiphysics.FORCE, displacement, model_part.Nodes)

        ##save the variable values
        KratosMultiphysics.VariableUtils().CopyScalarVar(KratosMultiphysics.VISCOSITY, KratosMultiphysics.DENSITY, model_part.Nodes)
        KratosMultiphysics.VariableUtils().CopyScalarVar(KratosMultiphysics.FORCE_X, KratosMultiphysics.REACTION_Y, model_part.Nodes)
        KratosMultiphysics.VariableUtils().CopyScalarVar(KratosMultiphysics.FORCE_X, KratosMultiphysics.FORCE_Y, model_part.Nodes)
        KratosMultiphysics.VariableUtils().CopyVectorVar(KratosMultiphysics.DISPLACEMENT, KratosMultiphysics.VELOCITY, model_part.Nodes)

        ##verify the result
        for node in model_part.Nodes:
            self.assertEqual(node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X), node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_X))
            self.assertEqual(node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y), node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_Y))
            self.assertEqual(node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Z), node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_Z))
            self.assertEqual(node.GetSolutionStepValue(KratosMultiphysics.FORCE_X), node.GetSolutionStepValue(KratosMultiphysics.REACTION_Y))
            self.assertEqual(node.GetSolutionStepValue(KratosMultiphysics.FORCE_X), node.GetSolutionStepValue(KratosMultiphysics.FORCE_Y))
            self.assertEqual(node.GetSolutionStepValue(KratosMultiphysics.VISCOSITY), node.GetSolutionStepValue(KratosMultiphysics.DENSITY))

    def test_save_var(self):
        current_model = KratosMultiphysics.Model()

        ##set the model part
        model_part = current_model.CreateModelPart("Main")
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DENSITY)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VISCOSITY)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        model_part_io = KratosMultiphysics.ModelPartIO(GetFilePath("auxiliar_files_for_python_unittest/mdpa_files/test_model_part_io_read"))
        model_part_io.ReadModelPart(model_part)

        ##save the variable values
        KratosMultiphysics.VariableUtils().SaveScalarVar(KratosMultiphysics.VISCOSITY, KratosMultiphysics.DENSITY, model_part.Nodes)
        KratosMultiphysics.VariableUtils().SaveVectorVar(KratosMultiphysics.DISPLACEMENT, KratosMultiphysics.VELOCITY, model_part.Nodes)
        KratosMultiphysics.VariableUtils().SaveScalarNonHistoricalVar(KratosMultiphysics.DENSITY, KratosMultiphysics.DISTANCE, model_part.Nodes)
        KratosMultiphysics.VariableUtils().SaveVectorNonHistoricalVar(KratosMultiphysics.VELOCITY, KratosMultiphysics.VOLUME_ACCELERATION, model_part.Nodes)

        ##verify the result
        for node in model_part.Nodes:
            self.assertEqual(node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X), node.GetValue(KratosMultiphysics.VELOCITY_X))
            self.assertEqual(node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y), node.GetValue(KratosMultiphysics.VELOCITY_Y))
            self.assertEqual(node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Z), node.GetValue(KratosMultiphysics.VELOCITY_Z))
            self.assertEqual(node.GetSolutionStepValue(KratosMultiphysics.VISCOSITY), node.GetValue(KratosMultiphysics.DENSITY))
            self.assertEqual(node.GetValue(KratosMultiphysics.VOLUME_ACCELERATION_X), node.GetValue(KratosMultiphysics.VELOCITY_X))
            self.assertEqual(node.GetValue(KratosMultiphysics.VOLUME_ACCELERATION_Y), node.GetValue(KratosMultiphysics.VELOCITY_Y))
            self.assertEqual(node.GetValue(KratosMultiphysics.VOLUME_ACCELERATION_Z), node.GetValue(KratosMultiphysics.VELOCITY_Z))
            self.assertEqual(node.GetValue(KratosMultiphysics.DISTANCE), node.GetValue(KratosMultiphysics.DENSITY))

    def test_set_variable_to_zero(self):
        ## Set the model part
        current_model = KratosMultiphysics.Model()
        model_part = current_model.CreateModelPart("Main")
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VISCOSITY)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        model_part_io = KratosMultiphysics.ModelPartIO(GetFilePath("auxiliar_files_for_python_unittest/mdpa_files/test_model_part_io_read"))
        model_part_io.ReadModelPart(model_part)

        ## Initialize the variable values
        for node in model_part.Elements:
            node.SetValue(KratosMultiphysics.VISCOSITY, node.Id)
            node.SetValue(KratosMultiphysics.DISPLACEMENT, [node.Id, 2 * node.Id, 3.0 * node.Id])
        for elem in model_part.Elements:
            elem.SetValue(KratosMultiphysics.VISCOSITY, elem.Id)
            elem.SetValue(KratosMultiphysics.DISPLACEMENT, [elem.Id, 2 * elem.Id, 3.0 * elem.Id])
        for cond in model_part.Conditions:
            cond.SetValue(KratosMultiphysics.VISCOSITY, cond.Id)
            cond.SetValue(KratosMultiphysics.DISPLACEMENT, [cond.Id, 2 * cond.Id, 3.0 * cond.Id])

        ## Set the variable values to zero
        KratosMultiphysics.VariableUtils().SetNonHistoricalVariableToZero(KratosMultiphysics.VISCOSITY, model_part.Nodes)
        KratosMultiphysics.VariableUtils().SetNonHistoricalVariableToZero(KratosMultiphysics.VISCOSITY, model_part.Elements)
        KratosMultiphysics.VariableUtils().SetNonHistoricalVariableToZero(KratosMultiphysics.VISCOSITY, model_part.Conditions)
        KratosMultiphysics.VariableUtils().SetNonHistoricalVariableToZero(KratosMultiphysics.DISPLACEMENT, model_part.Nodes)
        KratosMultiphysics.VariableUtils().SetNonHistoricalVariableToZero(KratosMultiphysics.DISPLACEMENT, model_part.Elements)
        KratosMultiphysics.VariableUtils().SetNonHistoricalVariableToZero(KratosMultiphysics.DISPLACEMENT, model_part.Conditions)

        ## Verify the results
        for node in model_part.Nodes:
            self.assertEqual(node.GetValue(KratosMultiphysics.VISCOSITY), 0.0)
            aux = node.GetValue(KratosMultiphysics.DISPLACEMENT)
            self.assertEqual(aux[0], 0.0)
            self.assertEqual(aux[1], 0.0)
            self.assertEqual(aux[2], 0.0)
        for elem in model_part.Elements:
            self.assertEqual(elem.GetValue(KratosMultiphysics.VISCOSITY), 0.0)
            aux = elem.GetValue(KratosMultiphysics.DISPLACEMENT)
            self.assertEqual(aux[0], 0.0)
            self.assertEqual(aux[1], 0.0)
            self.assertEqual(aux[2], 0.0)
        for cond in model_part.Conditions:
            self.assertEqual(cond.GetValue(KratosMultiphysics.VISCOSITY), 0.0)
            aux = cond.GetValue(KratosMultiphysics.DISPLACEMENT)
            self.assertEqual(aux[0], 0.0)
            self.assertEqual(aux[1], 0.0)
            self.assertEqual(aux[2], 0.0)

    def test_set_nodal_historical_variable_to_zero(self):
        ## Set the model part
        current_model = KratosMultiphysics.Model()
        model_part = current_model.CreateModelPart("Main")
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VISCOSITY)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        model_part_io = KratosMultiphysics.ModelPartIO(GetFilePath("auxiliar_files_for_python_unittest/mdpa_files/test_model_part_io_read"))
        model_part_io.ReadModelPart(model_part)

        ## Initialize the variable values
        for node in model_part.Nodes:
            node.SetSolutionStepValue(KratosMultiphysics.VISCOSITY, node.Id)
            node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT, [node.Id, 2.0 * node.Id, 3.0 * node.Id])

        ## Set the variable values to zero
        KratosMultiphysics.VariableUtils().SetHistoricalVariableToZero(KratosMultiphysics.VISCOSITY, model_part.Nodes)
        KratosMultiphysics.VariableUtils().SetHistoricalVariableToZero(KratosMultiphysics.DISPLACEMENT, model_part.Nodes)

        ## Verify the result
        for node in model_part.Nodes:
            self.assertEqual(node.GetSolutionStepValue(KratosMultiphysics.VISCOSITY), 0.0)
            self.assertEqual(node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X), 0.0)
            self.assertEqual(node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y), 0.0)
            self.assertEqual(node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Z), 0.0)

    def test_select_node_list(self):
        current_model = KratosMultiphysics.Model()

        ##set the model part
        model_part = current_model.CreateModelPart("Main")
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VISCOSITY)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        model_part_io = KratosMultiphysics.ModelPartIO(GetFilePath("auxiliar_files_for_python_unittest/mdpa_files/test_model_part_io_read"))
        model_part_io.ReadModelPart(model_part)

        ##extract the nodes with KratosMultiphysics.DISPLACEMENT_X equal to 0.0
        ids_list = []
        node_list = KratosMultiphysics.VariableUtils().SelectNodeList(KratosMultiphysics.VISCOSITY, 0.01, model_part.Nodes)
        for node in node_list:
            ids_list.append(node.Id)

        ##verify the result
        self.assertTrue(model_part.Nodes[1].Id in ids_list)
        self.assertTrue(model_part.Nodes[2].Id in ids_list)
        self.assertTrue(model_part.Nodes[973].Id in ids_list)
        self.assertTrue(model_part.Nodes[974].Id in ids_list)

    def test_apply_fixity(self):
        current_model = KratosMultiphysics.Model()

        ##set the model part
        model_part = current_model.CreateModelPart("Main")
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VISCOSITY)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        model_part_io = KratosMultiphysics.ModelPartIO(GetFilePath("auxiliar_files_for_python_unittest/mdpa_files/test_model_part_io_read"))
        model_part_io.ReadModelPart(model_part)

        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.VISCOSITY, model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_X, model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Y, model_part)

        ##apply the fixity
        KratosMultiphysics.VariableUtils().ApplyFixity(KratosMultiphysics.VISCOSITY, True, model_part.Nodes)
        KratosMultiphysics.VariableUtils().ApplyFixity(KratosMultiphysics.DISPLACEMENT_X, True, model_part.Nodes)
        KratosMultiphysics.VariableUtils().ApplyFixity(KratosMultiphysics.DISPLACEMENT_Y, False, model_part.Nodes)

        ##verify the result
        for node in model_part.Nodes:
            self.assertTrue(node.IsFixed(KratosMultiphysics.VISCOSITY))
            self.assertTrue(node.IsFixed(KratosMultiphysics.DISPLACEMENT_X))
            self.assertFalse(node.IsFixed(KratosMultiphysics.DISPLACEMENT_Y))

    def test_apply_flagged_fixity(self):
        current_model = KratosMultiphysics.Model()

        ##set the model part
        model_part = current_model.CreateModelPart("Main")
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VISCOSITY)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        model_part_io = KratosMultiphysics.ModelPartIO(GetFilePath("auxiliar_files_for_python_unittest/mdpa_files/test_model_part_io_read"))
        model_part_io.ReadModelPart(model_part)

        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.VISCOSITY, model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_X, model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Y, model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Z, model_part)

        for node in model_part.Nodes:
            if (node.Id % 3 == 0):
                node.Set(KratosMultiphysics.INLET)
            elif (node.Id % 3 == 1):
                node.Set(KratosMultiphysics.OUTLET)
            elif (node.Id % 3 == 2):
                node.Set(KratosMultiphysics.SLIP)

        ##apply the fixity
        KratosMultiphysics.VariableUtils().ApplyFixity(KratosMultiphysics.VISCOSITY, True, model_part.Nodes, KratosMultiphysics.INLET, True)
        KratosMultiphysics.VariableUtils().ApplyFixity(KratosMultiphysics.DISPLACEMENT_X, True, model_part.Nodes, KratosMultiphysics.OUTLET, False)
        KratosMultiphysics.VariableUtils().ApplyFixity(KratosMultiphysics.DISPLACEMENT_Y, False, model_part.Nodes, KratosMultiphysics.SLIP, True)
        KratosMultiphysics.VariableUtils().ApplyFixity(KratosMultiphysics.DISPLACEMENT_Z, True, model_part.Nodes, KratosMultiphysics.SLIP, False)

        ##verify the result
        for node in model_part.Nodes:
            node_type = node.Id % 3
            if (node_type == 0):
                self.assertTrue(node.IsFixed(KratosMultiphysics.VISCOSITY))
            if (node_type != 1):
                self.assertTrue(node.IsFixed(KratosMultiphysics.DISPLACEMENT_X))
            if (node_type == 2):
                self.assertFalse(node.IsFixed(KratosMultiphysics.DISPLACEMENT_Y))
            if (node_type != 2):
                self.assertTrue(node.IsFixed(KratosMultiphysics.DISPLACEMENT_Z))

    def test_apply_vector(self):
        current_model = KratosMultiphysics.Model()

        ##set the model part
        model_part = current_model.CreateModelPart("Main")
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VISCOSITY)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        model_part_io = KratosMultiphysics.ModelPartIO(GetFilePath("auxiliar_files_for_python_unittest/mdpa_files/test_model_part_io_read"))
        model_part_io.ReadModelPart(model_part)

        ##set the data vector [0,1,2,...]
        data_vector_x1 = KratosMultiphysics.Vector(len(model_part.Nodes))
        data_vector_x2 = KratosMultiphysics.Vector(len(model_part.Nodes))
        for i in range(len(model_part.Nodes)):
            data_vector_x1[i] = i
            data_vector_x2[i] = 2.0*i

        KratosMultiphysics.VariableUtils().ApplyVector(KratosMultiphysics.VISCOSITY, data_vector_x1, model_part.Nodes)
        KratosMultiphysics.VariableUtils().ApplyVector(KratosMultiphysics.DISPLACEMENT_X, data_vector_x2, model_part.Nodes)

        ##verify the result
        i = 0
        for node in model_part.Nodes:
            self.assertEqual(node.GetSolutionStepValue(KratosMultiphysics.VISCOSITY), data_vector_x1[i])
            self.assertEqual(node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X), data_vector_x2[i])
            i+=1

    def test_solutionsteps_vector(self):
        current_model = KratosMultiphysics.Model()

        ##set the model part
        model_part = current_model.CreateModelPart("Main")
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VISCOSITY)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        model_part.CreateNewNode(1,1.0,2.0,3.0)
        model_part.CreateNewNode(2,11.0,25.0,31.0)
        model_part.CreateNewNode(10,12.0,26.0,32.0)

        #reference_results
        ref_initial_coords2d = [1.0,2.0,11.0,25.0,12.0,26.0]
        ref_initial_coords3d = [1.0,2.0,3.0,11.0,25.0,31.0,12.0,26.0,32.0]
        ref_current_coords2d = [2.0,4.0,12.0,27.0,13.0,28.0]
        ref_current_coords3d = [2.0,4.0,6.0,12.0,27.0,34.0,13.0,28.0,35.0]
        ref_values2d = [1.0,2.0,2.0,4.0,10.0,20.0]
        ref_values3d = [1.0,2.0,3.0,2.0,4.0,6.0,10.0,20.0,30.0]

        for node in model_part.Nodes:
            node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X,0,node.Id*1.)
            node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y,0,node.Id*2.)
            node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Z,0,node.Id*3.)
            node.SetSolutionStepValue(KratosMultiphysics.VISCOSITY,0,node.X)
            node.X = node.X0 + 1.0
            node.Y = node.Y0 + 2.0
            node.Z = node.Z0 + 3.0

        init2d = KratosMultiphysics.VariableUtils().GetInitialPositionsVector(model_part.Nodes,2) 
        init3d = KratosMultiphysics.VariableUtils().GetInitialPositionsVector(model_part.Nodes,3)
        current2d = KratosMultiphysics.VariableUtils().GetCurrentPositionsVector(model_part.Nodes,2) 
        current3d = KratosMultiphysics.VariableUtils().GetCurrentPositionsVector(model_part.Nodes,3)

        for v,r in zip(ref_initial_coords2d,init2d):
            self.assertEqual(v,r)

        for v,r in zip(ref_initial_coords3d,init3d):
            self.assertEqual(v,r)

        for v,r in zip(ref_current_coords2d,current2d):
            self.assertEqual(v,r)

        for v,r in zip(ref_current_coords3d,current3d):
            self.assertEqual(v,r)

        values2d = KratosMultiphysics.VariableUtils().GetSolutionStepValuesVector(model_part.Nodes,KratosMultiphysics.DISPLACEMENT,0,2)
        values3d = KratosMultiphysics.VariableUtils().GetSolutionStepValuesVector(model_part.Nodes,KratosMultiphysics.DISPLACEMENT,0,3)

        for v,r in zip(ref_values2d,values2d):
            self.assertEqual(v,r)

        for v,r in zip(ref_values3d,values3d):
            self.assertEqual(v,r)

        input2d = [0.0,1.0,2.0,3.0,4.0,5.0]
        KratosMultiphysics.VariableUtils().SetSolutionStepValuesVector(model_part.Nodes,KratosMultiphysics.DISPLACEMENT,input2d,0)
        verify2d = KratosMultiphysics.VariableUtils().GetSolutionStepValuesVector(model_part.Nodes,KratosMultiphysics.DISPLACEMENT,0,2) 

        for v,r in zip(input2d,verify2d):
            self.assertEqual(v,r)

        input3d = [0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0]
        KratosMultiphysics.VariableUtils().SetSolutionStepValuesVector(model_part.Nodes,KratosMultiphysics.DISPLACEMENT,input2d,0)
        verify3d = KratosMultiphysics.VariableUtils().GetSolutionStepValuesVector(model_part.Nodes,KratosMultiphysics.DISPLACEMENT,0,2) 

        for v,r in zip(input3d,verify3d):
            self.assertEqual(v,r)

    def test_sum_variable(self):
        current_model = KratosMultiphysics.Model()

        ##set the model part
        model_part = current_model.CreateModelPart("Main")
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DENSITY)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VISCOSITY)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        model_part_io = KratosMultiphysics.ModelPartIO(GetFilePath("auxiliar_files_for_python_unittest/mdpa_files/test_model_part_io_read"))
        model_part_io.ReadModelPart(model_part)

        scalar_val = 0.1
        vector_val = KratosMultiphysics.Vector(3)
        vector_val[0] = 1.0
        vector_val[1] = 2.0
        vector_val[2] = 3.0

        ##set non-historical nodal values (historical ones coming from the mdpa)
        for node in model_part.Nodes:
            node.SetValue(KratosMultiphysics.DENSITY, scalar_val)
            node.SetValue(KratosMultiphysics.VELOCITY, vector_val)

        ##set non-historical condition values
        for condition in model_part.Conditions:
            condition.SetValue(KratosMultiphysics.DENSITY, scalar_val)
            condition.SetValue(KratosMultiphysics.VELOCITY, vector_val)

        ##set non-historical element values
        for element in model_part.Elements:
            element.SetValue(KratosMultiphysics.DENSITY, scalar_val)
            element.SetValue(KratosMultiphysics.VELOCITY, vector_val)

        ##sum nodal variables
        sum_hist_scal = KratosMultiphysics.VariableUtils().SumHistoricalNodeScalarVariable(KratosMultiphysics.VISCOSITY, model_part, 0)
        sum_hist_vect = KratosMultiphysics.VariableUtils().SumHistoricalNodeVectorVariable(KratosMultiphysics.DISPLACEMENT, model_part, 0)
        sum_nonhist_scal = KratosMultiphysics.VariableUtils().SumNonHistoricalNodeScalarVariable(KratosMultiphysics.DENSITY, model_part)
        sum_nonhist_vect = KratosMultiphysics.VariableUtils().SumNonHistoricalNodeVectorVariable(KratosMultiphysics.VELOCITY, model_part)

        ##verify the nodal results
        self.assertAlmostEqual(sum_hist_scal, 0.04, delta=1e-6)
        self.assertAlmostEqual(sum_hist_vect[0], 0.3, delta=1e-6)
        self.assertAlmostEqual(sum_hist_vect[1], 0.001947, delta=1e-6)
        self.assertAlmostEqual(sum_hist_vect[2], 0.0, delta=1e-6)
        self.assertAlmostEqual(sum_nonhist_scal, 0.6, delta=1e-6)
        self.assertAlmostEqual(sum_nonhist_vect[0], 6.0, delta=1e-6)
        self.assertAlmostEqual(sum_nonhist_vect[1], 12.0, delta=1e-6)
        self.assertAlmostEqual(sum_nonhist_vect[2], 18.0, delta=1e-6)

        ##sum condition variables
        cond_scal = KratosMultiphysics.VariableUtils().SumConditionScalarVariable(KratosMultiphysics.DENSITY, model_part)
        cond_vect = KratosMultiphysics.VariableUtils().SumConditionVectorVariable(KratosMultiphysics.VELOCITY, model_part)

        ##verify the condition results
        self.assertAlmostEqual(cond_scal, 0.5, delta=1e-6)
        self.assertAlmostEqual(cond_vect[0], 5.0, delta=1e-6)
        self.assertAlmostEqual(cond_vect[1], 10.0, delta=1e-6)
        self.assertAlmostEqual(cond_vect[2], 15.0, delta=1e-6)

        ##sum element variables
        elem_scal = KratosMultiphysics.VariableUtils().SumElementScalarVariable(KratosMultiphysics.DENSITY, model_part)
        elem_vect = KratosMultiphysics.VariableUtils().SumElementVectorVariable(KratosMultiphysics.VELOCITY, model_part)

        ##verify the element results
        self.assertAlmostEqual(elem_scal, 0.4, delta=1e-6)
        self.assertAlmostEqual(elem_vect[0], 4.0, delta=1e-6)
        self.assertAlmostEqual(elem_vect[1], 8.0, delta=1e-6)
        self.assertAlmostEqual(elem_vect[2], 12.0, delta=1e-6)

    def test_UpdateCurrentToInitialConfiguration(self):
        current_model = KratosMultiphysics.Model()

        ##set the model part
        model_part = current_model.CreateModelPart("Main")
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DENSITY)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VISCOSITY)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        model_part_io = KratosMultiphysics.ModelPartIO(GetFilePath("auxiliar_files_for_python_unittest/mdpa_files/test_model_part_io_read"))
        model_part_io.ReadModelPart(model_part)

        ##set the reference model part
        ref_model_part = current_model.CreateModelPart("Reference")
        ref_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DENSITY)
        ref_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        ref_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VISCOSITY)
        ref_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        ref_model_part_io = KratosMultiphysics.ModelPartIO(GetFilePath("auxiliar_files_for_python_unittest/mdpa_files/test_model_part_io_read"))
        ref_model_part_io.ReadModelPart(ref_model_part)

        dx = 0.1
        dy = -0.2
        dz = 0.3

        # update the current configuration, ONLY in master_mp, NOT in the reference
        for node in model_part.Nodes:
            node.X += dx
            node.Y += dy
            node.Z += dz

        KratosMultiphysics.VariableUtils().UpdateCurrentToInitialConfiguration(model_part.Nodes)

        # we set the updated configuration back to the initial configuration
        # therefore testing against the initial configuration of the reference MP
        for node, node_ref in zip(model_part.Nodes, ref_model_part.Nodes):
            self.assertAlmostEqual(node.X0, node_ref.X0)
            self.assertAlmostEqual(node.Y0, node_ref.Y0)
            self.assertAlmostEqual(node.Z0, node_ref.Z0)

            self.assertAlmostEqual(node.X, node_ref.X0)
            self.assertAlmostEqual(node.Y, node_ref.Y0)
            self.assertAlmostEqual(node.Z, node_ref.Z0)


    def test_UpdateInitialToCurrentConfiguration(self):
        current_model = KratosMultiphysics.Model()

        ##set the model part
        model_part = current_model.CreateModelPart("Main")
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DENSITY)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VISCOSITY)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        model_part_io = KratosMultiphysics.ModelPartIO(GetFilePath("auxiliar_files_for_python_unittest/mdpa_files/test_model_part_io_read"))
        model_part_io.ReadModelPart(model_part)

        ##set the reference model part
        ref_model_part = current_model.CreateModelPart("Reference")
        ref_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DENSITY)
        ref_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        ref_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VISCOSITY)
        ref_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        ref_model_part_io = KratosMultiphysics.ModelPartIO(GetFilePath("auxiliar_files_for_python_unittest/mdpa_files/test_model_part_io_read"))
        ref_model_part_io.ReadModelPart(ref_model_part)

        dx = 0.1
        dy = -0.2
        dz = 0.3

        # update the current configuration, in BOTH ModelParts!
        for node, node_ref in zip(model_part.Nodes, ref_model_part.Nodes):
            node.X += dx
            node.Y += dy
            node.Z += dz
            node_ref.X += dx
            node_ref.Y += dy
            node_ref.Z += dz

        KratosMultiphysics.VariableUtils().UpdateInitialToCurrentConfiguration(model_part.Nodes)

        # we set the initial configuration to be the same as the current configuration
        # therefore testing against the current configuration of the reference MP
        for node, node_ref in zip(model_part.Nodes, ref_model_part.Nodes):
            self.assertAlmostEqual(node.X0, node_ref.X)
            self.assertAlmostEqual(node.Y0, node_ref.Y)
            self.assertAlmostEqual(node.Z0, node_ref.Z)

            self.assertAlmostEqual(node.X, node_ref.X)
            self.assertAlmostEqual(node.Y, node_ref.Y)
            self.assertAlmostEqual(node.Z, node_ref.Z)

    def test_UpdateCurrentPosition(self):
        # Set the test model part
        current_model = KratosMultiphysics.Model()
        model_part = current_model.CreateModelPart("Main")
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DENSITY)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VISCOSITY)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.MESH_DISPLACEMENT)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        model_part_io = KratosMultiphysics.ModelPartIO(GetFilePath("auxiliar_files_for_python_unittest/mdpa_files/test_model_part_io_read"))
        model_part_io.ReadModelPart(model_part)

        # Set a fake displacement field
        for node in model_part.Nodes:
            aux_disp = KratosMultiphysics.Vector(3)
            aux_disp[0] = float(node.Id)
            aux_disp[1] = 1.5 * float(node.Id)
            aux_disp[2] = 2.0 * float(node.Id)
            node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT, aux_disp)

        # Update current position
        KratosMultiphysics.VariableUtils().UpdateCurrentPosition(model_part.Nodes)

        # Check current position
        for node in model_part.Nodes:
            self.assertAlmostEqual(node.X, node.X0 + float(node.Id))
            self.assertAlmostEqual(node.Y, node.Y0 + 1.5 * float(node.Id))
            self.assertAlmostEqual(node.Z, node.Z0 + 2.0 * float(node.Id))

        # Set a fake displacement field in another variable
        for node in model_part.Nodes:
            aux_disp = KratosMultiphysics.Vector(3)
            aux_disp[0] = 3.0 * float(node.Id)
            aux_disp[1] = 4.0 * float(node.Id)
            aux_disp[2] = 5.0 * float(node.Id)
            node.SetSolutionStepValue(KratosMultiphysics.MESH_DISPLACEMENT, aux_disp)

        # Update current position using an alternative variable
        KratosMultiphysics.VariableUtils().UpdateCurrentPosition(model_part.Nodes, KratosMultiphysics.MESH_DISPLACEMENT)

        # Check current position
        for node in model_part.Nodes:
            self.assertAlmostEqual(node.X, node.X0 + 3.0 * float(node.Id))
            self.assertAlmostEqual(node.Y, node.Y0 + 4.0 * float(node.Id))
            self.assertAlmostEqual(node.Z, node.Z0 + 5.0 * float(node.Id))

    def testAddDofsList(self):
        # Set the test model part
        current_model = KratosMultiphysics.Model()
        model_part = current_model.CreateModelPart("Main")
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VISCOSITY)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        model_part_io = KratosMultiphysics.ModelPartIO(GetFilePath("auxiliar_files_for_python_unittest/mdpa_files/test_model_part_io_read"))
        model_part_io.ReadModelPart(model_part)

        # Set the DOFs list to be set
        dofs_list = ["VELOCITY_X","VELOCITY_Y","DISPLACEMENT_Z"]

        # Update current position
        KratosMultiphysics.VariableUtils.AddDofsList(dofs_list, model_part)

        # Check that the nodes have the DOFs added
        for node in model_part.Nodes:
            self.assertTrue(node.HasDofFor(KratosMultiphysics.VELOCITY_X))
            self.assertTrue(node.HasDofFor(KratosMultiphysics.VELOCITY_Y))
            self.assertTrue(node.HasDofFor(KratosMultiphysics.DISPLACEMENT_Z))

    def testAddDofsWithReactionList(self):
        # Set the test model part
        current_model = KratosMultiphysics.Model()
        model_part = current_model.CreateModelPart("Main")
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VISCOSITY)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        model_part_io = KratosMultiphysics.ModelPartIO(GetFilePath("auxiliar_files_for_python_unittest/mdpa_files/test_model_part_io_read"))
        model_part_io.ReadModelPart(model_part)

        # Set the DOFs list to be set
        dofs_with_reaction_list = [["VELOCITY_X","REACTION_X"],["VELOCITY_Y","REACTION_Y"],["DISPLACEMENT_Z","REACTION_Z"]]

        # Update current position
        KratosMultiphysics.VariableUtils.AddDofsList(dofs_with_reaction_list, model_part)

        # Check that the nodes have the DOFs added
        for node in model_part.Nodes:
            self.assertTrue(node.HasDofFor(KratosMultiphysics.VELOCITY_X))
            self.assertTrue(node.HasDofFor(KratosMultiphysics.VELOCITY_Y))
            self.assertTrue(node.HasDofFor(KratosMultiphysics.DISPLACEMENT_Z))

    def test_distribute_condition_variable(self):
        current_model = KratosMultiphysics.Model()

        ##set the model part
        model_part = current_model.CreateModelPart("Main")
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VISCOSITY)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        model_part_io = KratosMultiphysics.ModelPartIO(GetFilePath("auxiliar_files_for_python_unittest/mdpa_files/test_model_part_io_read"))
        model_part_io.ReadModelPart(model_part)

        for node in model_part.Nodes:
            node.SetValue(KratosMultiphysics.AUX_MESH_VAR, node.Id)

        for condition in model_part.Conditions:
            condition.SetValue(KratosMultiphysics.DISTANCE, condition.Id)
            vector = KratosMultiphysics.Vector(3)
            vector[0] = condition.Id * 3
            vector[1] = condition.Id * 3 + 1
            vector[2] = condition.Id * 3 + 2
            condition.SetValue(KratosMultiphysics.VELOCITY, vector)

        variable_utils = KratosMultiphysics.VariableUtils()
        variable_utils.WeightedAccumulateConditionVariableOnNodes(model_part, KratosMultiphysics.DISTANCE, KratosMultiphysics.AUX_MESH_VAR, False)
        variable_utils.WeightedAccumulateConditionVariableOnNodes(model_part, KratosMultiphysics.VELOCITY, KratosMultiphysics.AUX_MESH_VAR, False)

        distance_vector = [
            1.0, 3602.0, 10803.0, 3643056.0, 3789835.0, 1897352.0
        ]
        velocity_vector = [
            3.0, 4.0, 5.0, 10806.0, 10810.0, 10814.0, 32409.0, 32415.0, 32421.0, 10929168.0, 10931112.0, 10933056.0, 11369505.0, 11371451.0, 11373397.0, 5692056.0, 5693030.0, 5694004.0
        ]

        local_index = 0
        for node in model_part.Nodes:
            self.assertEqual(node.GetValue(KratosMultiphysics.DISTANCE), distance_vector[local_index])
            self.assertEqual(node.GetValue(KratosMultiphysics.VELOCITY)[0], velocity_vector[local_index * 3])
            self.assertEqual(node.GetValue(KratosMultiphysics.VELOCITY)[1], velocity_vector[local_index * 3 + 1])
            self.assertEqual(node.GetValue(KratosMultiphysics.VELOCITY)[2], velocity_vector[local_index * 3 + 2])

            local_index += 1

    def test_distribute_condition_variable_inverse(self):
        current_model = KratosMultiphysics.Model()

        ##set the model part
        model_part = current_model.CreateModelPart("Main")
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VISCOSITY)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        model_part_io = KratosMultiphysics.ModelPartIO(GetFilePath("auxiliar_files_for_python_unittest/mdpa_files/test_model_part_io_read"))
        model_part_io.ReadModelPart(model_part)

        for node in model_part.Nodes:
            node.SetValue(KratosMultiphysics.AUX_MESH_VAR, 1.0 / node.Id)

        for condition in model_part.Conditions:
            condition.SetValue(KratosMultiphysics.DISTANCE, condition.Id)
            vector = KratosMultiphysics.Vector(3)
            vector[0] = condition.Id * 3
            vector[1] = condition.Id * 3 + 1
            vector[2] = condition.Id * 3 + 2
            condition.SetValue(KratosMultiphysics.VELOCITY, vector)

        variable_utils = KratosMultiphysics.VariableUtils()
        variable_utils.WeightedAccumulateConditionVariableOnNodes(model_part, KratosMultiphysics.DISTANCE, KratosMultiphysics.AUX_MESH_VAR, True)
        variable_utils.WeightedAccumulateConditionVariableOnNodes(model_part, KratosMultiphysics.VELOCITY, KratosMultiphysics.AUX_MESH_VAR, True)

        distance_vector = [
            1.0, 3602.0, 10803.0, 3643056.0, 3789835.0, 1897352.0
        ]
        velocity_vector = [
            3.0, 4.0, 5.0, 10806.0, 10810.0, 10814.0, 32409.0, 32415.0, 32421.0, 10929168.0, 10931112.0, 10933056.0, 11369505.0, 11371451.0, 11373397.0, 5692056.0, 5693030.0, 5694004.0
        ]

        local_index = 0
        for node in model_part.Nodes:
            self.assertAlmostEqual(node.GetValue(KratosMultiphysics.DISTANCE), distance_vector[local_index])
            self.assertAlmostEqual(node.GetValue(KratosMultiphysics.VELOCITY)[0], velocity_vector[local_index * 3])
            self.assertAlmostEqual(node.GetValue(KratosMultiphysics.VELOCITY)[1], velocity_vector[local_index * 3 + 1])
            self.assertAlmostEqual(node.GetValue(KratosMultiphysics.VELOCITY)[2], velocity_vector[local_index * 3 + 2])

            local_index += 1

    def __test_copy_model_part_flagged_methods(
        self,
        method,
        get_entity_container,
        origin_entity_value_setter,
        origin_entity_value_getter,
        destination_entity_value_setter,
        destination_entity_value_getter,
        add_origin_variables_to_solution_step_variables_list = True,
        add_destination_variables_to_solution_step_variables_list = True):
        ##set the origin model part
        current_model = KratosMultiphysics.Model()
        origin_model_part = current_model.CreateModelPart("OriginModelPart")
        if (add_origin_variables_to_solution_step_variables_list):
            origin_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DENSITY)
            origin_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VISCOSITY)
            origin_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
            origin_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        origin_model_part.SetBufferSize(2)
        model_part_io = KratosMultiphysics.ModelPartIO(GetFilePath("auxiliar_files_for_python_unittest/mdpa_files/test_model_part_io_read"), KratosMultiphysics.ModelPartIO.READ | KratosMultiphysics.ModelPartIO.IGNORE_VARIABLES_ERROR)
        model_part_io.ReadModelPart(origin_model_part)

        ##set the destination model part
        destination_model_part = current_model.CreateModelPart("DestinationModelPart")
        destination_model_part.SetBufferSize(2)
        if (add_destination_variables_to_solution_step_variables_list):
            destination_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VISCOSITY)
            destination_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DENSITY)
            destination_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
            destination_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        model_part_io = KratosMultiphysics.ModelPartIO(GetFilePath("auxiliar_files_for_python_unittest/mdpa_files/test_model_part_io_read"), KratosMultiphysics.ModelPartIO.READ | KratosMultiphysics.ModelPartIO.IGNORE_VARIABLES_ERROR)
        model_part_io.ReadModelPart(destination_model_part)

        ##set the values in the origin model part
        for entity in get_entity_container(origin_model_part):
            entity.Set(KratosMultiphysics.INLET, entity.Id % 2)
            origin_entity_value_setter(entity, KratosMultiphysics.VISCOSITY, entity.Id * 3)
            origin_entity_value_setter(entity, KratosMultiphysics.DISPLACEMENT, [entity.Id * 3, entity.Id * 7, entity.Id * 23])

        ##  initialize the containers in destination model part (otherwise the operation is not threadsafe!)
        for entity in get_entity_container(destination_model_part):
            destination_entity_value_setter(entity, KratosMultiphysics.VISCOSITY, 0.0)
            destination_entity_value_setter(entity, KratosMultiphysics.DENSITY, 0.0)
            destination_entity_value_setter(entity, KratosMultiphysics.DISPLACEMENT, KratosMultiphysics.Array3(0.0))
            destination_entity_value_setter(entity, KratosMultiphysics.VELOCITY, KratosMultiphysics.Array3(0.0))

        ##copy the values to the destination model part
        method(origin_model_part, destination_model_part, KratosMultiphysics.VISCOSITY, KratosMultiphysics.DENSITY)
        method(origin_model_part, destination_model_part, KratosMultiphysics.VISCOSITY, KratosMultiphysics.VISCOSITY)
        method(origin_model_part, destination_model_part, KratosMultiphysics.DISPLACEMENT, KratosMultiphysics.DISPLACEMENT)
        method(origin_model_part, destination_model_part, KratosMultiphysics.DISPLACEMENT, KratosMultiphysics.VELOCITY)

        ##check the copied values
        for origin_entity, destination_entity in zip(get_entity_container(origin_model_part), get_entity_container(destination_model_part)):
            if (origin_entity.Is(KratosMultiphysics.INLET)):
                self.assertEqual(origin_entity_value_getter(origin_entity, KratosMultiphysics.VISCOSITY), destination_entity_value_getter(destination_entity, KratosMultiphysics.VISCOSITY))
                self.assertEqual(origin_entity_value_getter(origin_entity, KratosMultiphysics.VISCOSITY), destination_entity_value_getter(destination_entity, KratosMultiphysics.DENSITY))
                self.assertEqual(origin_entity_value_getter(origin_entity, KratosMultiphysics.DISPLACEMENT_X), destination_entity_value_getter(destination_entity, KratosMultiphysics.DISPLACEMENT_X))
                self.assertEqual(origin_entity_value_getter(origin_entity, KratosMultiphysics.DISPLACEMENT_Y), destination_entity_value_getter(destination_entity, KratosMultiphysics.DISPLACEMENT_Y))
                self.assertEqual(origin_entity_value_getter(origin_entity, KratosMultiphysics.DISPLACEMENT_X), destination_entity_value_getter(destination_entity, KratosMultiphysics.VELOCITY_X))
                self.assertEqual(origin_entity_value_getter(origin_entity, KratosMultiphysics.DISPLACEMENT_Y), destination_entity_value_getter(destination_entity, KratosMultiphysics.VELOCITY_Y))
            else:
                self.assertEqual(destination_entity_value_getter(destination_entity, KratosMultiphysics.VISCOSITY), 0.0)
                self.assertEqual(destination_entity_value_getter(destination_entity, KratosMultiphysics.DENSITY), 0.0)
                self.assertEqual(destination_entity_value_getter(destination_entity, KratosMultiphysics.DISPLACEMENT_X), 0.0)
                self.assertEqual(destination_entity_value_getter(destination_entity, KratosMultiphysics.DISPLACEMENT_Y), 0.0)
                self.assertEqual(destination_entity_value_getter(destination_entity, KratosMultiphysics.VELOCITY_X), 0.0)
                self.assertEqual(destination_entity_value_getter(destination_entity, KratosMultiphysics.VELOCITY_Y), 0.0)

if __name__ == '__main__':
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    KratosUnittest.main()

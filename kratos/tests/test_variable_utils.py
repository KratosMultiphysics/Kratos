from __future__ import print_function, absolute_import, division

import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics import *
import math
import os

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)


class TestVariableUtils(KratosUnittest.TestCase):

    def test_copy_model_part_nodal_var(self):
        current_model = Model()

        ##set the origin model part
        origin_model_part = current_model.CreateModelPart("OriginModelPart")
        origin_model_part.AddNodalSolutionStepVariable(VISCOSITY)
        origin_model_part.AddNodalSolutionStepVariable(DISPLACEMENT)
        origin_model_part.SetBufferSize(2)
        model_part_io = ModelPartIO(GetFilePath("auxiliar_files/mdpa_files/test_model_part_io_read"))
        model_part_io.ReadModelPart(origin_model_part)

        ##set the destination model part
        destination_model_part = current_model.CreateModelPart("DestinationModelPart")
        destination_model_part.AddNodalSolutionStepVariable(VISCOSITY)
        destination_model_part.AddNodalSolutionStepVariable(DISPLACEMENT)
        destination_model_part.SetBufferSize(2)
        model_part_io = ModelPartIO(GetFilePath("auxiliar_files/mdpa_files/test_model_part_io_read"))
        model_part_io.ReadModelPart(destination_model_part)

        ##set the values in the destination model part
        for node in origin_model_part.Nodes:
            node.SetSolutionStepValue(VISCOSITY, 0, node.X + node.Y)
            node.SetSolutionStepValue(VISCOSITY, 1, 2.0 * node.X + 3.0 * node.Y)
            node.SetSolutionStepValue(DISPLACEMENT, 0, [node.X ** 2, 0.0, 0.0])
            node.SetSolutionStepValue(DISPLACEMENT, 1, [node.X, node.Y, node.Z])

        ##copy the values to the destination model part
        VariableUtils().CopyModelPartNodalVar(VISCOSITY, origin_model_part, destination_model_part, 0)
        VariableUtils().CopyModelPartNodalVar(VISCOSITY, origin_model_part, destination_model_part, 1)
        VariableUtils().CopyModelPartNodalVar(DISPLACEMENT, origin_model_part, destination_model_part, 0)
        VariableUtils().CopyModelPartNodalVar(DISPLACEMENT, origin_model_part, destination_model_part, 1)

        ##check the copied values
        for node in destination_model_part.Nodes:
            self.assertEqual(node.GetSolutionStepValue(VISCOSITY, 0), node.X + node.Y)
            self.assertEqual(node.GetSolutionStepValue(VISCOSITY, 1), 2.0 * node.X + 3.0 * node.Y)
            self.assertEqual(node.GetSolutionStepValue(DISPLACEMENT_X, 0), node.X ** 2)
            self.assertEqual(node.GetSolutionStepValue(DISPLACEMENT_X, 1), node.X)

    def test_copy_model_part_elemental_var(self):
        current_model = Model()

        ##set the origin model part
        origin_model_part = current_model.CreateModelPart("OriginModelPart")
        origin_model_part.AddNodalSolutionStepVariable(VISCOSITY)
        origin_model_part.AddNodalSolutionStepVariable(DISPLACEMENT)
        model_part_io = ModelPartIO(GetFilePath("auxiliar_files/mdpa_files/test_model_part_io_read"))
        model_part_io.ReadModelPart(origin_model_part)

        ##set the destination model part
        destination_model_part = current_model.CreateModelPart("DestinationModelPart")
        destination_model_part.AddNodalSolutionStepVariable(VISCOSITY)
        destination_model_part.AddNodalSolutionStepVariable(DISPLACEMENT)
        model_part_io = ModelPartIO(GetFilePath("auxiliar_files/mdpa_files/test_model_part_io_read"))
        model_part_io.ReadModelPart(destination_model_part)
        ##set the values in the destination model part
        for element in origin_model_part.Elements:
            element.SetValue(DENSITY, element.Id*100)
            element.SetValue(VOLUME_ACCELERATION, [element.Id*100, 0.0, 0.0])

        ##copy the values to the destination model part
        VariableUtils().CopyModelPartElementalVar(DENSITY, origin_model_part, destination_model_part)
        VariableUtils().CopyModelPartElementalVar(VOLUME_ACCELERATION, origin_model_part, destination_model_part)

        ##check the copied values
        for element in destination_model_part.Elements:
            self.assertEqual(element.GetValue(DENSITY), element.Id*100)
            self.assertEqual(element.GetValue(VOLUME_ACCELERATION)[0], element.Id*100)


    def test_set_variable(self):
        current_model = Model()

        ##set the model part
        model_part = current_model.CreateModelPart("Main")
        model_part.AddNodalSolutionStepVariable(VISCOSITY)
        model_part.AddNodalSolutionStepVariable(DISPLACEMENT)
        model_part_io = ModelPartIO(GetFilePath("auxiliar_files/mdpa_files/test_model_part_io_read"))
        model_part_io.ReadModelPart(model_part)

        ##set the variable values
        viscosity = 0.1
        displacement = Vector(3)
        displacement[0] = 1.0
        displacement[1] = 2.0
        displacement[2] = 3.0

        VariableUtils().SetScalarVar(VISCOSITY, viscosity, model_part.Nodes)
        VariableUtils().SetVectorVar(DISPLACEMENT, displacement, model_part.Nodes)

        ##verify the result
        for node in model_part.Nodes:
            self.assertEqual(node.GetSolutionStepValue(DISPLACEMENT_X), 1.0)
            self.assertEqual(node.GetSolutionStepValue(DISPLACEMENT_Y), 2.0)
            self.assertEqual(node.GetSolutionStepValue(DISPLACEMENT_Z), 3.0)
            self.assertEqual(node.GetSolutionStepValue(VISCOSITY), viscosity)

    def test_set_nonhistorical_variable(self):
        current_model = Model()

        ##set the model part
        model_part = current_model.CreateModelPart("Main")
        model_part.AddNodalSolutionStepVariable(VISCOSITY)
        model_part.AddNodalSolutionStepVariable(DISPLACEMENT)
        model_part_io = ModelPartIO(GetFilePath("auxiliar_files/mdpa_files/test_model_part_io_read"))
        model_part_io.ReadModelPart(model_part)

        ##set the variable values
        viscosity = 0.1
        displacement = Vector(3)
        displacement[0] = 1.0
        displacement[1] = 2.0
        displacement[2] = 3.0

        # First for nodes
        VariableUtils().SetNonHistoricalVariable(VISCOSITY, viscosity, model_part.Nodes)
        VariableUtils().SetNonHistoricalVariable(DISPLACEMENT, displacement, model_part.Nodes)

        ##verify the result
        for node in model_part.Nodes:
            self.assertEqual(node.GetValue(DISPLACEMENT_X), 1.0)
            self.assertEqual(node.GetValue(DISPLACEMENT_Y), 2.0)
            self.assertEqual(node.GetValue(DISPLACEMENT_Z), 3.0)
            self.assertEqual(node.GetValue(VISCOSITY), viscosity)

        # Now for conditions (it will work for elements too)
        VariableUtils().SetNonHistoricalVariable(VISCOSITY, viscosity, model_part.Conditions)
        VariableUtils().SetNonHistoricalVariable(DISPLACEMENT, displacement, model_part.Conditions)

        ##verify the result
        for cond in model_part.Conditions:
            disp = cond.GetValue(DISPLACEMENT)
            self.assertEqual(disp[0], 1.0)
            self.assertEqual(disp[1], 2.0)
            self.assertEqual(disp[2], 3.0)
            self.assertEqual(cond.GetValue(VISCOSITY), viscosity)

    def test_set_flag(self):
        current_model = Model()

        ##set the model part
        model_part = current_model.CreateModelPart("Main")
        model_part.AddNodalSolutionStepVariable(VISCOSITY)
        model_part.AddNodalSolutionStepVariable(DISPLACEMENT)
        model_part_io = ModelPartIO(GetFilePath("auxiliar_files/mdpa_files/test_model_part_io_read"))
        model_part_io.ReadModelPart(model_part)

        VariableUtils().SetFlag(VISITED, True, model_part.Nodes)
        VariableUtils().SetFlag(VISITED, True, model_part.Conditions)
        VariableUtils().SetFlag(VISITED, True, model_part.Elements)

        VariableUtils().SetFlag(INLET, True, model_part.GetSubModelPart("Inlets").Nodes)
        VariableUtils().SetFlag(INLET, True, model_part.GetSubModelPart("Inlets").Conditions)
        VariableUtils().SetFlag(OUTLET, False, model_part.GetSubModelPart("Inlets").Nodes)
        VariableUtils().SetFlag(OUTLET, False, model_part.GetSubModelPart("Inlets").Conditions)

        ##verify the main modelpart flags set
        for node in model_part.Nodes:
            self.assertTrue(node.Is(VISITED))
        for condition in model_part.Conditions:
            self.assertTrue(condition.Is(VISITED))
        for element in model_part.Elements:
            self.assertTrue(element.Is(VISITED))
        ##verify the inlet submodelpart flag set
        for node in model_part.GetSubModelPart("Inlets").Nodes:
            self.assertTrue(node.Is(INLET))
            self.assertTrue(node.IsNot(OUTLET))
        for condition in model_part.GetSubModelPart("Inlets").Conditions:
            self.assertTrue(condition.Is(INLET))
            self.assertTrue(condition.IsNot(OUTLET))

    def test_copy_var(self):
        current_model = Model()

        ##set the model part
        model_part = current_model.CreateModelPart("Main")
        model_part.AddNodalSolutionStepVariable(DENSITY)
        model_part.AddNodalSolutionStepVariable(VELOCITY)
        model_part.AddNodalSolutionStepVariable(VISCOSITY)
        model_part.AddNodalSolutionStepVariable(DISPLACEMENT)
        model_part.AddNodalSolutionStepVariable(FORCE)
        model_part.AddNodalSolutionStepVariable(REACTION)
        model_part_io = ModelPartIO(GetFilePath("auxiliar_files/mdpa_files/test_model_part_io_read"))
        model_part_io.ReadModelPart(model_part)

        ##set the variable values
        viscosity = 0.1
        displacement = Vector(3)
        displacement[0] = 1.3
        displacement[1] = 2.2
        displacement[2] = 3.1

        VariableUtils().SetScalarVar(VISCOSITY, viscosity, model_part.Nodes)
        VariableUtils().SetVectorVar(DISPLACEMENT, displacement, model_part.Nodes)
        VariableUtils().SetVectorVar(FORCE, displacement, model_part.Nodes)

        ##save the variable values
        VariableUtils().CopyScalarVar(VISCOSITY, DENSITY, model_part.Nodes)
        VariableUtils().CopyComponentVar(FORCE_X, REACTION_Y, model_part.Nodes)
        VariableUtils().CopyComponentVar(FORCE_X, FORCE_Y, model_part.Nodes)
        VariableUtils().CopyVectorVar(DISPLACEMENT, VELOCITY, model_part.Nodes)

        ##verify the result
        for node in model_part.Nodes:
            self.assertEqual(node.GetSolutionStepValue(DISPLACEMENT_X), node.GetSolutionStepValue(VELOCITY_X))
            self.assertEqual(node.GetSolutionStepValue(DISPLACEMENT_Y), node.GetSolutionStepValue(VELOCITY_Y))
            self.assertEqual(node.GetSolutionStepValue(DISPLACEMENT_Z), node.GetSolutionStepValue(VELOCITY_Z))
            self.assertEqual(node.GetSolutionStepValue(FORCE_X), node.GetSolutionStepValue(REACTION_Y))
            self.assertEqual(node.GetSolutionStepValue(FORCE_X), node.GetSolutionStepValue(FORCE_Y))
            self.assertEqual(node.GetSolutionStepValue(VISCOSITY), node.GetSolutionStepValue(DENSITY))

    def test_save_var(self):
        current_model = Model()

        ##set the model part
        model_part = current_model.CreateModelPart("Main")
        model_part.AddNodalSolutionStepVariable(DENSITY)
        model_part.AddNodalSolutionStepVariable(VELOCITY)
        model_part.AddNodalSolutionStepVariable(VISCOSITY)
        model_part.AddNodalSolutionStepVariable(DISPLACEMENT)
        model_part_io = ModelPartIO(GetFilePath("auxiliar_files/mdpa_files/test_model_part_io_read"))
        model_part_io.ReadModelPart(model_part)

        ##save the variable values
        VariableUtils().SaveScalarVar(VISCOSITY, DENSITY, model_part.Nodes)
        VariableUtils().SaveVectorVar(DISPLACEMENT, VELOCITY, model_part.Nodes)
        VariableUtils().SaveScalarNonHistoricalVar(DENSITY, DISTANCE, model_part.Nodes)
        VariableUtils().SaveVectorNonHistoricalVar(VELOCITY, VOLUME_ACCELERATION, model_part.Nodes)

        ##verify the result
        for node in model_part.Nodes:
            self.assertEqual(node.GetSolutionStepValue(DISPLACEMENT_X), node.GetValue(VELOCITY_X))
            self.assertEqual(node.GetSolutionStepValue(DISPLACEMENT_Y), node.GetValue(VELOCITY_Y))
            self.assertEqual(node.GetSolutionStepValue(DISPLACEMENT_Z), node.GetValue(VELOCITY_Z))
            self.assertEqual(node.GetSolutionStepValue(VISCOSITY), node.GetValue(DENSITY))
            self.assertEqual(node.GetValue(VOLUME_ACCELERATION_X), node.GetValue(VELOCITY_X))
            self.assertEqual(node.GetValue(VOLUME_ACCELERATION_Y), node.GetValue(VELOCITY_Y))
            self.assertEqual(node.GetValue(VOLUME_ACCELERATION_Z), node.GetValue(VELOCITY_Z))
            self.assertEqual(node.GetValue(DISTANCE), node.GetValue(DENSITY))

    def test_set_to_zero(self):
        current_model = Model()

        ##set the model part
        model_part = current_model.CreateModelPart("Main")
        model_part.AddNodalSolutionStepVariable(VISCOSITY)
        model_part.AddNodalSolutionStepVariable(DISPLACEMENT)
        model_part_io = ModelPartIO(GetFilePath("auxiliar_files/mdpa_files/test_model_part_io_read"))
        model_part_io.ReadModelPart(model_part)

        ##save the variable values
        VariableUtils().SetToZero_ScalarVar(VISCOSITY, model_part.Nodes)
        VariableUtils().SetToZero_VectorVar(DISPLACEMENT, model_part.Nodes)

        ##verify the result
        for node in model_part.Nodes:
            self.assertEqual(node.GetSolutionStepValue(DISPLACEMENT_X), 0.0)
            self.assertEqual(node.GetSolutionStepValue(DISPLACEMENT_Y), 0.0)
            self.assertEqual(node.GetSolutionStepValue(DISPLACEMENT_Z), 0.0)
            self.assertEqual(node.GetSolutionStepValue(VISCOSITY), 0.0)

    def test_select_node_list(self):
        current_model = Model()

        ##set the model part
        model_part = current_model.CreateModelPart("Main")
        model_part.AddNodalSolutionStepVariable(VISCOSITY)
        model_part.AddNodalSolutionStepVariable(DISPLACEMENT)
        model_part_io = ModelPartIO(GetFilePath("auxiliar_files/mdpa_files/test_model_part_io_read"))
        model_part_io.ReadModelPart(model_part)

        ##extract the nodes with DISPLACEMENT_X equal to 0.0
        ids_list = []
        node_list = VariableUtils().SelectNodeList(VISCOSITY, 0.01, model_part.Nodes)
        for node in node_list:
            ids_list.append(node.Id)

        ##verify the result
        self.assertTrue(model_part.Nodes[1].Id in ids_list)
        self.assertTrue(model_part.Nodes[2].Id in ids_list)
        self.assertTrue(model_part.Nodes[973].Id in ids_list)
        self.assertTrue(model_part.Nodes[974].Id in ids_list)

    def test_apply_fixity(self):
        current_model = Model()

        ##set the model part
        model_part = current_model.CreateModelPart("Main")
        model_part.AddNodalSolutionStepVariable(VISCOSITY)
        model_part.AddNodalSolutionStepVariable(DISPLACEMENT)
        model_part_io = ModelPartIO(GetFilePath("auxiliar_files/mdpa_files/test_model_part_io_read"))
        model_part_io.ReadModelPart(model_part)

        ##apply the fixity
        VariableUtils().ApplyFixity(VISCOSITY, True, model_part.Nodes)
        VariableUtils().ApplyFixity(DISPLACEMENT_X, True, model_part.Nodes)
        VariableUtils().ApplyFixity(DISPLACEMENT_Y, False, model_part.Nodes)

        ##verify the result
        for node in model_part.Nodes:
            self.assertTrue(node.IsFixed(VISCOSITY))
            self.assertTrue(node.IsFixed(DISPLACEMENT_X))
            self.assertFalse(node.IsFixed(DISPLACEMENT_Y))

    def test_apply_vector(self):
        current_model = Model()

        ##set the model part
        model_part = current_model.CreateModelPart("Main")
        model_part.AddNodalSolutionStepVariable(VISCOSITY)
        model_part.AddNodalSolutionStepVariable(DISPLACEMENT)
        model_part_io = ModelPartIO(GetFilePath("auxiliar_files/mdpa_files/test_model_part_io_read"))
        model_part_io.ReadModelPart(model_part)

        ##set the data vector [0,1,2,...]
        data_vector_x1 = Vector(len(model_part.Nodes))
        data_vector_x2 = Vector(len(model_part.Nodes))
        for i in range(len(model_part.Nodes)):
            data_vector_x1[i] = i
            data_vector_x2[i] = 2.0*i

        VariableUtils().ApplyVector(VISCOSITY, data_vector_x1, model_part.Nodes)
        VariableUtils().ApplyVector(DISPLACEMENT_X, data_vector_x2, model_part.Nodes)

        ##verify the result
        i = 0
        for node in model_part.Nodes:
            self.assertEqual(node.GetSolutionStepValue(VISCOSITY), data_vector_x1[i])
            self.assertEqual(node.GetSolutionStepValue(DISPLACEMENT_X), data_vector_x2[i])
            i+=1

    def test_sum_variable(self):
        current_model = Model()

        ##set the model part
        model_part = current_model.CreateModelPart("Main")
        model_part.AddNodalSolutionStepVariable(DENSITY)
        model_part.AddNodalSolutionStepVariable(VELOCITY)
        model_part.AddNodalSolutionStepVariable(VISCOSITY)
        model_part.AddNodalSolutionStepVariable(DISPLACEMENT)
        model_part_io = ModelPartIO(GetFilePath("auxiliar_files/mdpa_files/test_model_part_io_read"))
        model_part_io.ReadModelPart(model_part)

        scalar_val = 0.1
        vector_val = Vector(3)
        vector_val[0] = 1.0
        vector_val[1] = 2.0
        vector_val[2] = 3.0

        ##set non-historical nodal values (historical ones coming from the mdpa)
        for node in model_part.Nodes:
            node.SetValue(DENSITY, scalar_val)
            node.SetValue(VELOCITY, vector_val)

        ##set non-historical condition values
        for condition in model_part.Conditions:
            condition.SetValue(DENSITY, scalar_val)
            condition.SetValue(VELOCITY, vector_val)

        ##set non-historical element values
        for element in model_part.Elements:
            element.SetValue(DENSITY, scalar_val)
            element.SetValue(VELOCITY, vector_val)

        ##sum nodal variables
        sum_hist_scal = VariableUtils().SumHistoricalNodeScalarVariable(VISCOSITY, model_part, 0)
        sum_hist_vect = VariableUtils().SumHistoricalNodeVectorVariable(DISPLACEMENT, model_part, 0)
        sum_nonhist_scal = VariableUtils().SumNonHistoricalNodeScalarVariable(DENSITY, model_part)
        sum_nonhist_vect = VariableUtils().SumNonHistoricalNodeVectorVariable(VELOCITY, model_part)

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
        cond_scal = VariableUtils().SumConditionScalarVariable(DENSITY, model_part)
        cond_vect = VariableUtils().SumConditionVectorVariable(VELOCITY, model_part)

        ##verify the condition results
        self.assertAlmostEqual(cond_scal, 0.5, delta=1e-6)
        self.assertAlmostEqual(cond_vect[0], 5.0, delta=1e-6)
        self.assertAlmostEqual(cond_vect[1], 10.0, delta=1e-6)
        self.assertAlmostEqual(cond_vect[2], 15.0, delta=1e-6)

        ##sum element variables
        elem_scal = VariableUtils().SumElementScalarVariable(DENSITY, model_part)
        elem_vect = VariableUtils().SumElementVectorVariable(VELOCITY, model_part)

        ##verify the element results
        self.assertAlmostEqual(elem_scal, 0.4, delta=1e-6)
        self.assertAlmostEqual(elem_vect[0], 4.0, delta=1e-6)
        self.assertAlmostEqual(elem_vect[1], 8.0, delta=1e-6)
        self.assertAlmostEqual(elem_vect[2], 12.0, delta=1e-6)

    def test_UpdateCurrentToInitialConfiguration(self):
        current_model = Model()

        ##set the model part
        model_part = current_model.CreateModelPart("Main")
        model_part.AddNodalSolutionStepVariable(DENSITY)
        model_part.AddNodalSolutionStepVariable(VELOCITY)
        model_part.AddNodalSolutionStepVariable(VISCOSITY)
        model_part.AddNodalSolutionStepVariable(DISPLACEMENT)
        model_part_io = ModelPartIO(GetFilePath("auxiliar_files/mdpa_files/test_model_part_io_read"))
        model_part_io.ReadModelPart(model_part)

        ##set the reference model part
        ref_model_part = current_model.CreateModelPart("Reference")
        ref_model_part.AddNodalSolutionStepVariable(DENSITY)
        ref_model_part.AddNodalSolutionStepVariable(VELOCITY)
        ref_model_part.AddNodalSolutionStepVariable(VISCOSITY)
        ref_model_part.AddNodalSolutionStepVariable(DISPLACEMENT)
        ref_model_part_io = ModelPartIO(GetFilePath("auxiliar_files/mdpa_files/test_model_part_io_read"))
        ref_model_part_io.ReadModelPart(ref_model_part)

        dx = 0.1
        dy = -0.2
        dz = 0.3

        # update the current configuration, ONLY in master_mp, NOT in the reference
        for node in model_part.Nodes:
            node.X += dx
            node.Y += dy
            node.Z += dz

        VariableUtils().UpdateCurrentToInitialConfiguration(model_part.Nodes)

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
        current_model = Model()

        ##set the model part
        model_part = current_model.CreateModelPart("Main")
        model_part.AddNodalSolutionStepVariable(DENSITY)
        model_part.AddNodalSolutionStepVariable(VELOCITY)
        model_part.AddNodalSolutionStepVariable(VISCOSITY)
        model_part.AddNodalSolutionStepVariable(DISPLACEMENT)
        model_part_io = ModelPartIO(GetFilePath("auxiliar_files/mdpa_files/test_model_part_io_read"))
        model_part_io.ReadModelPart(model_part)

        ##set the reference model part
        ref_model_part = current_model.CreateModelPart("Reference")
        ref_model_part.AddNodalSolutionStepVariable(DENSITY)
        ref_model_part.AddNodalSolutionStepVariable(VELOCITY)
        ref_model_part.AddNodalSolutionStepVariable(VISCOSITY)
        ref_model_part.AddNodalSolutionStepVariable(DISPLACEMENT)
        ref_model_part_io = ModelPartIO(GetFilePath("auxiliar_files/mdpa_files/test_model_part_io_read"))
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

        VariableUtils().UpdateInitialToCurrentConfiguration(model_part.Nodes)

        # we set the initial configuration to be the same as the current configuration
        # therefore testing against the current configuration of the reference MP
        for node, node_ref in zip(model_part.Nodes, ref_model_part.Nodes):
            self.assertAlmostEqual(node.X0, node_ref.X)
            self.assertAlmostEqual(node.Y0, node_ref.Y)
            self.assertAlmostEqual(node.Z0, node_ref.Z)

            self.assertAlmostEqual(node.X, node_ref.X)
            self.assertAlmostEqual(node.Y, node_ref.Y)
            self.assertAlmostEqual(node.Z, node_ref.Z)


if __name__ == '__main__':
    Logger.GetDefaultOutput().SetSeverity(Logger.Severity.WARNING)
    KratosUnittest.main()

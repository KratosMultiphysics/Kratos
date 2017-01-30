from __future__ import print_function, absolute_import, division

import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics import *
import math

def GetFilePath(fileName):
    return os.path.dirname(os.path.realpath(__file__)) + "/" + fileName


class TestVariableUtils(KratosUnittest.TestCase):

    def test_set_variable(self):
        ##set the model part
        model_part = ModelPart("Main")
        model_part.AddNodalSolutionStepVariable(VISCOSITY)
        model_part.AddNodalSolutionStepVariable(DISPLACEMENT)
        model_part_io = ModelPartIO(GetFilePath("test_model_part_io"))
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

    def test_copy_var(self):
        ##set the model part
        model_part = ModelPart("Main")
        model_part.AddNodalSolutionStepVariable(DENSITY)
        model_part.AddNodalSolutionStepVariable(VELOCITY)
        model_part.AddNodalSolutionStepVariable(VISCOSITY)
        model_part.AddNodalSolutionStepVariable(DISPLACEMENT)
        model_part_io = ModelPartIO(GetFilePath("test_model_part_io"))
        model_part_io.ReadModelPart(model_part)

        ##save the variable values
        VariableUtils().CopyScalarVar(VISCOSITY, DENSITY, model_part.Nodes)
        VariableUtils().CopyVectorVar(DISPLACEMENT, VELOCITY, model_part.Nodes)

        ##verify the result
        for node in model_part.Nodes:
            self.assertEqual(node.GetSolutionStepValue(DISPLACEMENT_X), node.GetSolutionStepValue(VELOCITY_X))
            self.assertEqual(node.GetSolutionStepValue(DISPLACEMENT_Y), node.GetSolutionStepValue(VELOCITY_Y))
            self.assertEqual(node.GetSolutionStepValue(DISPLACEMENT_Z), node.GetSolutionStepValue(VELOCITY_Z))
            self.assertEqual(node.GetSolutionStepValue(VISCOSITY), node.GetSolutionStepValue(DENSITY))

    def test_save_var(self):
        ##set the model part
        model_part = ModelPart("Main")
        model_part.AddNodalSolutionStepVariable(DENSITY)
        model_part.AddNodalSolutionStepVariable(VELOCITY)
        model_part.AddNodalSolutionStepVariable(VISCOSITY)
        model_part.AddNodalSolutionStepVariable(DISPLACEMENT)
        model_part_io = ModelPartIO(GetFilePath("test_model_part_io"))
        model_part_io.ReadModelPart(model_part)

        ##save the variable values
        VariableUtils().SaveScalarVar(VISCOSITY, DENSITY, model_part.Nodes)
        VariableUtils().SaveVectorVar(DISPLACEMENT, VELOCITY, model_part.Nodes)

        ##verify the result
        for node in model_part.Nodes:
            self.assertEqual(node.GetSolutionStepValue(DISPLACEMENT_X), node.GetValue(VELOCITY_X))
            self.assertEqual(node.GetSolutionStepValue(DISPLACEMENT_Y), node.GetValue(VELOCITY_Y))
            self.assertEqual(node.GetSolutionStepValue(DISPLACEMENT_Z), node.GetValue(VELOCITY_Z))
            self.assertEqual(node.GetSolutionStepValue(VISCOSITY), node.GetValue(DENSITY))

    def test_set_to_zero(self):
        ##set the model part
        model_part = ModelPart("Main")
        model_part.AddNodalSolutionStepVariable(VISCOSITY)
        model_part.AddNodalSolutionStepVariable(DISPLACEMENT)
        model_part_io = ModelPartIO(GetFilePath("test_model_part_io"))
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
        ##set the model part
        model_part = ModelPart("Main")
        model_part.AddNodalSolutionStepVariable(VISCOSITY)
        model_part.AddNodalSolutionStepVariable(DISPLACEMENT)
        model_part_io = ModelPartIO(GetFilePath("test_model_part_io"))
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
        ##set the model part
        model_part = ModelPart("Main")
        model_part.AddNodalSolutionStepVariable(VISCOSITY)
        model_part.AddNodalSolutionStepVariable(DISPLACEMENT)
        model_part_io = ModelPartIO(GetFilePath("test_model_part_io"))
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
        ##set the model part
        model_part = ModelPart("Main")
        model_part.AddNodalSolutionStepVariable(VISCOSITY)
        model_part.AddNodalSolutionStepVariable(DISPLACEMENT)
        model_part_io = ModelPartIO(GetFilePath("test_model_part_io"))
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
        ##set the model part
        model_part = ModelPart("Main")
        model_part.AddNodalSolutionStepVariable(DENSITY)
        model_part.AddNodalSolutionStepVariable(VELOCITY)
        model_part.AddNodalSolutionStepVariable(VISCOSITY)
        model_part.AddNodalSolutionStepVariable(DISPLACEMENT)
        model_part_io = ModelPartIO(GetFilePath("test_model_part_io"))
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


if __name__ == '__main__':
    KratosUnittest.main()

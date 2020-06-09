import KratosMultiphysics as Kratos
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD
from KratosMultiphysics.FluidDynamicsApplication import check_and_prepare_model_process_fluid

import KratosMultiphysics.KratosUnittest as UnitTest


class CFDUtilitiesTest(UnitTest.TestCase):
    def testCalculateLinearLogarithmicWallFunctionBasedYPlusLimit(self):
        kappa = 0.5
        beta = 6.0
        y_plus = KratosCFD.CFDUtilities.CalculateLinearLogarithmicWallFunctionBasedYPlusLimit(
            kappa, beta)
        u_plus = math.log(y_plus) / kappa + beta
        self.assertAlmostEqual(u_plus, y_plus, 6)

    def testCalculateLinearLogarithmicWallFunctionBasedYPlusAndUtau(self):
        u_tau = Kratos.Array3()
        velocity = Kratos.Array3()
        normal = Kratos.Array3([1.0, 0.5, 0.5])
        nu = 10.0
        y = 1.1
        kappa = 0.5
        beta = 6.0

        unit_normal = normal * (1.0 / self.__CalculateArray3DNorm(normal))

        # checking for linear region
        velocity[0] = 0.1
        velocity[1] = 0.2
        velocity[2] = 0.3
        y_plus = KratosCFD.CFDUtilities.CalculateLinearLogarithmicWallFunctionBasedYPlusAndUtau(
            u_tau, velocity, normal, nu, y, kappa, beta)
        tangential_velocity = velocity - unit_normal * (
            velocity[0] * unit_normal[0] + velocity[1] * unit_normal[1] +
            velocity[2] * unit_normal[2])
        u_tau_magnitude = math.sqrt(
            self.__CalculateArray3DNorm(tangential_velocity) * nu / y)
        u_tau_check = tangential_velocity * (
            u_tau_magnitude * -1.0 /
            self.__CalculateArray3DNorm(tangential_velocity))
        y_plus_check = u_tau_magnitude * y / nu

        self.assertVectorAlmostEqual(u_tau, u_tau_check, 7)
        self.assertAlmostEqual(y_plus, y_plus_check, 7)

        # checking for logarithmic region
        velocity[0] = 1000.0
        velocity[1] = 2000.0
        velocity[2] = 3000.0
        y_plus = KratosCFD.CFDUtilities.CalculateLinearLogarithmicWallFunctionBasedYPlusAndUtau(
            u_tau, velocity, normal, nu, y, kappa, beta)
        tangential_velocity = velocity - unit_normal * (
            velocity[0] * unit_normal[0] + velocity[1] * unit_normal[1] +
            velocity[2] * unit_normal[2])
        u_tau_magnitude = y_plus * nu / y
        u_tau_check = tangential_velocity * (
            u_tau_magnitude * -1.0 /
            self.__CalculateArray3DNorm(tangential_velocity))
        y_plus_check = math.exp(
            (self.__CalculateArray3DNorm(tangential_velocity) /
             self.__CalculateArray3DNorm(u_tau) - beta) * kappa)

        self.assertVectorAlmostEqual(u_tau, u_tau_check, 7)
        self.assertAlmostEqual(y_plus, y_plus_check, 7)

    def testCalculateReactionBasedYPlusUTau(self):
        u_tau = Kratos.Array3()
        reaction = Kratos.Array3([1.0, 2.0, 5.0])
        normal = Kratos.Array3([1.0, 0.5, 0.5])
        nu = 10.0
        y = 1.1
        rho = 1.2

        y_plus = KratosCFD.CFDUtilities.CalculateReactionBasedYPlusUTau(
            u_tau, reaction, normal, rho, nu, y)
        surface_area = self.__CalculateArray3DNorm(normal)
        unit_normal = normal * (1.0 / surface_area)

        tangential_reaction_check = reaction - unit_normal * (
            reaction[0] * unit_normal[0] + reaction[1] * unit_normal[1] +
            reaction[2] * unit_normal[2])
        u_tau_magnitude = self.__CalculateArray3DNorm(u_tau)
        tangential_reaction = u_tau * (u_tau_magnitude**2 * rho *
                                       surface_area / u_tau_magnitude)
        y_plus_check = u_tau_magnitude * y / nu
        self.assertVectorAlmostEqual(tangential_reaction,
                                     tangential_reaction_check, 7)
        self.assertAlmostEqual(y_plus, y_plus_check, 7)

    def testCalculateNumberOfNeighbourConditions(self):
        self.__CreateModel()
        model_part = self.model.GetModelPart("test.submodelpart_1")
        KratosCFD.CFDUtilities.CalculateNumberOfNeighbourConditions(model_part)

        check_values = [1, 2, 1]
        for node, value in zip(model_part.Nodes, check_values):
            node_value = node.GetValue(
                KratosCFD.NUMBER_OF_NEIGHBOUR_CONDITIONS)
            self.assertEqual(node_value, value)

    def testCalculateYPlusAndUTauForConditionsBasedOnReaction(self):
        self.__CreateModel()
        model_part = self.model.GetModelPart("test.submodelpart_1")
        KratosCFD.CFDUtilities.CalculateYPlusAndUTauForConditionsBasedOnReaction(
            model_part, Kratos.VISCOSITY, Kratos.REACTION)

        check_y_plus_values = [154.37981820205698, 40.0630232620064]
        check_u_tau_values = [
            Kratos.Array3([1.2995046003231017, 0.0, 1.6293788450205047]),
            Kratos.Array3([1.8308118481697873, 0.0, 0.1722050748278512])
        ]
        for condition, y_plus_check, u_tau_check in zip(
                model_part.Conditions, check_y_plus_values,
                check_u_tau_values):
            self.assertAlmostEqual(condition.GetValue(KratosCFD.Y_PLUS),
                                   y_plus_check, 7)
            self.assertVectorAlmostEqual(
                condition.GetValue(KratosCFD.FRICTION_VELOCITY), u_tau_check,
                7)

    def testCalculateYPlusAndUTauForConditionsBasedOnLinearLogarithmicWallFunction(
        self):
        self.__CreateModel()
        model_part = self.model.GetModelPart("test.submodelpart_1")

        kappa = 0.5
        beta = 10.0
        KratosCFD.CFDUtilities.CalculateYPlusAndUTauForConditionsBasedOnLinearLogarithmicWallFunction(
            model_part, Kratos.VISCOSITY, kappa, beta)

        check_y_plus_values = [42.184008341432026, 14.937830572275276]
        check_u_tau_values = [
            Kratos.Array3([-0.3145718416438325, 0.0, -0.4747175064806927]),
            Kratos.Array3([-0.6091915392914995, 0.0, -0.31463738842527994])
        ]
        for condition, y_plus_check, u_tau_check in zip(
                model_part.Conditions, check_y_plus_values,
                check_u_tau_values):
            self.assertAlmostEqual(condition.GetValue(KratosCFD.Y_PLUS),
                                   y_plus_check, 7)
            self.assertVectorAlmostEqual(
                condition.GetValue(KratosCFD.FRICTION_VELOCITY), u_tau_check,
                7)

    def __CalculateArray3DNorm(self, value):
        return math.sqrt(value[0]**2 + value[1]**2 + value[2]**2)


# test model part is as follows
#       6   5   4
#       .---.---.
#      / |4 /\ 3|
#     / 1| /  \ |
#    /   |/  2 \|
#   .----.------.
#   1    2      3

    def __CreateModel(self,
                      variable_list=[
                          Kratos.VELOCITY, Kratos.REACTION, Kratos.DENSITY,
                          Kratos.VISCOSITY, Kratos.MESH_VELOCITY
                      ]):
        self.model = Kratos.Model()
        self.model_part = self.model.CreateModelPart("test")
        for variable in variable_list:
            self.model_part.AddNodalSolutionStepVariable(variable)

        self.model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
        self.model_part.CreateNewNode(2, 1.0, 0.0, 0.0)
        self.model_part.CreateNewNode(3, 2.0, 0.0, 0.0)
        self.model_part.CreateNewNode(4, 2.0, 1.0, 0.0)
        self.model_part.CreateNewNode(5, 1.5, 1.0, 0.0)
        self.model_part.CreateNewNode(6, 1.0, 1.0, 0.0)
        prop = self.model_part.GetProperties()[0]

        self.model_part.CreateNewElement("Element2D3N", 1, [1, 2, 6], prop)
        self.model_part.CreateNewElement("Element2D3N", 2, [2, 3, 5], prop)
        self.model_part.CreateNewElement("Element2D3N", 3, [3, 4, 5], prop)
        self.model_part.CreateNewElement("Element2D3N", 4, [2, 5, 6], prop)

        self.submodelpart_1 = self.model_part.CreateSubModelPart(
            "submodelpart_1")
        self.submodelpart_1.AddNodes([1, 2, 3])
        self.submodelpart_1.CreateNewCondition("LineCondition2D2N", 1, [1, 2],
                                               prop)
        self.submodelpart_1.CreateNewCondition("LineCondition2D2N", 2, [2, 3],
                                               prop)

        self.submodelpart_2 = self.model_part.CreateSubModelPart(
            "submodelpart_2")
        self.submodelpart_2.AddNodes([3, 4, 5])
        self.submodelpart_2.CreateNewCondition("LineCondition2D2N", 3, [3, 4],
                                               prop)
        self.submodelpart_2.CreateNewCondition("LineCondition2D2N", 4, [4, 5],
                                               prop)

        self.submodelpart_3 = self.model_part.CreateSubModelPart(
            "submodelpart_3")
        self.submodelpart_3.AddNodes([1, 6, 5])
        self.submodelpart_3.CreateNewCondition("LineCondition2D2N", 5, [5, 6],
                                               prop)
        self.submodelpart_3.CreateNewCondition("LineCondition2D2N", 6, [6, 1],
                                               prop)

        prepare_model_part_settings = Kratos.Parameters("{}")

        prepare_model_part_settings.AddEmptyValue("volume_model_part_name")
        prepare_model_part_settings["volume_model_part_name"].SetString("test")

        prepare_model_part_settings.AddEmptyArray("skin_parts")
        prepare_model_part_settings["skin_parts"].Append("submodelpart_1")
        prepare_model_part_settings["skin_parts"].Append("submodelpart_2")
        prepare_model_part_settings["skin_parts"].Append("submodelpart_3")

        prepare_model_part_settings.AddEmptyValue(
            "assign_neighbour_elements_to_conditions")
        prepare_model_part_settings[
            "assign_neighbour_elements_to_conditions"].SetBool(True)

        check_and_prepare_model_process_fluid.CheckAndPrepareModelProcess(
            self.model_part, prepare_model_part_settings).Execute()

        self.model_part.SetBufferSize(2)
        self.model_part.ProcessInfo[Kratos.DOMAIN_SIZE] = 2

        for node in self.model_part.Nodes:
            nid = node.Id
            x = node.X + 0.6
            y = node.Y + 1.2
            z = node.Z + 3.0

            vector = Kratos.Vector(3)
            vector[0] = nid * x + y * z
            vector[1] = nid * x - y * z
            vector[2] = -1.0 * nid * x + y * 1 + 3 * z
            if (node.SolutionStepsDataHas(Kratos.VELOCITY)):
                vector = node.SetSolutionStepValue(Kratos.VELOCITY, 0, vector)

            vector = Kratos.Vector(3)
            vector[0] = nid * x + y * z * 2
            vector[1] = nid * x - y * z * 5
            vector[2] = -2.0 * nid * x + y * 4 + 3 * z
            if (node.SolutionStepsDataHas(Kratos.REACTION)):
                vector = node.SetSolutionStepValue(Kratos.REACTION, 0, vector)

            scalar = nid * y + z
            if (node.SolutionStepsDataHas(Kratos.DENSITY)):
                node.SetSolutionStepValue(Kratos.DENSITY, 0, scalar)

            scalar = (nid * z * x - y) * 1e-3
            if (node.SolutionStepsDataHas(Kratos.VISCOSITY)):
                node.SetSolutionStepValue(Kratos.VISCOSITY, 0, scalar)

if __name__ == '__main__':
    UnitTest.main()

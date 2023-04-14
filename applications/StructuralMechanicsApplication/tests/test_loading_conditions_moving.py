import KratosMultiphysics

import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import KratosMultiphysics.KratosUnittest as KratosUnittest
import math

class TestLoadingConditionsMoving(KratosUnittest.TestCase):


    def __CalculateReaction2n(self, length, angle, load, distance):
        """
        Calculate reaction forces on 2 noded 2d inclined element

        Parameters
        ----------
        length length of the element
        angle angle of the element
        load load vector
        distance local distance of the load within the element

        Returns reaction forces
        -------

        """

        cos = math.cos(angle)
        sin = math.sin(angle)

        p_shear = load[1] * sin + load[0] * -cos
        a = distance
        b = length - distance

        left_shear_reaction = p_shear * b/length
        right_shear_reaction = p_shear * a/length

        x_l = left_shear_reaction * -sin
        y_l = left_shear_reaction * cos

        x_r = right_shear_reaction * -sin
        y_r = right_shear_reaction * cos

        return (x_l, y_l), (x_r, y_r)


    def __CalculateReaction3n(self, length, angle, load, distance):
        """
        Calculate reaction forces on 3 noded 2d inclined element

        Parameters
        ----------
        length length of the element
        angle angle of the element
        load load vector
        distance local distance of the load within the element

        Returns reaction forces
        -------

        """

        cos = math.cos(angle)
        sin = math.sin(angle)

        # calculate local forces
        p_shear = load[1] * sin + load[0] * -cos
        p_normal = load[1] * cos + load[0] * sin

        # calculate local reaction forces
        left_shear_reaction = p_shear * (length**2 - 3*length*distance + 2*distance**2)/length**2
        right_shear_reaction = p_shear * (2*distance**2 - length*distance)/length**2
        mid_shear_reaction = p_shear * 4 * (length * distance - distance ** 2) / length ** 2

        left_normal_reaction = p_normal * (length**2 - 3*length*distance + 2*distance**2)/length**2
        right_normal_reaction = p_normal * (2*distance**2 - length*distance)/length**2
        mid_normal_reaction = p_normal * 4 * (length * distance - distance ** 2) / length ** 2

        # calculate global force
        x_l = left_shear_reaction * -sin + left_normal_reaction * cos
        y_l = left_shear_reaction * cos + left_normal_reaction * sin

        x_mid = mid_shear_reaction * -sin + mid_normal_reaction * cos
        y_mid = mid_shear_reaction * cos + mid_normal_reaction * sin

        x_r = right_shear_reaction * -sin + right_normal_reaction * cos
        y_r = right_shear_reaction * cos + right_normal_reaction * sin

        return (x_l, y_l), (x_mid, y_mid), (x_r, y_r)

    def __CalculateReactionWithRotation(self, length, angle, load, distance):
        """
        Calculate reaction forces including moment on a 2d inclined element

        Parameters
        ----------
        length length of the element
        angle angle of the element
        load load vector
        distance local distance of the load within the element

        Returns reaction forces
        -------

        """

        cos = math.cos(angle)
        sin = math.sin(angle)

        p_shear = load[1] * sin + load[0] * -cos
        a = distance
        b = length - distance

        expected_left_moment = p_shear * a * b ** 2 / length ** 2
        expected_right_moment = -p_shear * a ** 2 * b / length ** 2
        moment_left, moment_right = expected_left_moment, expected_right_moment

        left_shear_reaction = (p_shear * b) / length + (moment_left + moment_right) / length
        right_shear_reaction = (p_shear * a) / length - (moment_left + moment_right) / length

        x_l = left_shear_reaction * -sin
        y_l = left_shear_reaction * cos

        x_r = right_shear_reaction * -sin
        y_r = right_shear_reaction * cos

        return (x_l, y_l, moment_left), (x_r, y_r, moment_right)

    def _MovingLoadCondition3D2NRotDofZ(self):
        """"
        Check moving load on a 3d 2n line element, where the line element is inclined around z axis.
        """

        current_model = KratosMultiphysics.Model()

        mp = current_model.CreateModelPart("solid_part")
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.ROTATION)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION_MOMENT)

        # create nodes
        second_coord = [math.sqrt(2), math.sqrt(2.0), 0.0]
        mp.CreateNewNode(1, 0.0, 0.0, 0.0)
        mp.CreateNewNode(2, second_coord[0], second_coord[1], second_coord[2])

        length = 2.0
        rotation = math.atan(second_coord[1] / second_coord[0])

        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.REACTION_X, mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Y, KratosMultiphysics.REACTION_Y, mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Z, KratosMultiphysics.REACTION_Z, mp)

        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ROTATION_X, KratosMultiphysics.REACTION_MOMENT_X, mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ROTATION_Y, KratosMultiphysics.REACTION_MOMENT_Y, mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ROTATION_Z, KratosMultiphysics.REACTION_MOMENT_Z, mp)

        # create condition
        cond = mp.CreateNewCondition("MovingLoadCondition3D2N", 1, [1, 2], mp.GetProperties()[1])

        # initialise matrices
        lhs = KratosMultiphysics.Matrix(0, 0)
        rhs = KratosMultiphysics.Vector(0)

        # set perpendicular POINT_LOAD to the condition
        load_on_cond = KratosMultiphysics.Vector(3)
        load_on_cond[0] = 2.0
        load_on_cond[1] = -2.0
        load_on_cond[2] = 0.0

        cond.SetValue(StructuralMechanicsApplication.POINT_LOAD, load_on_cond)

        # compare results when load is located at different points
        num_tests = 5
        for i in range(num_tests):
            distance = i/(num_tests-1) * length
            (x_left, y_left, moment_left), (x_right, y_right, moment_right) = \
                self.__CalculateReactionWithRotation(length, rotation, load_on_cond, distance)

            # set load on quarter
            cond.SetValue(StructuralMechanicsApplication.MOVING_LOAD_LOCAL_DISTANCE, distance)
            cond.CalculateLocalSystem(lhs, rhs, mp.ProcessInfo)

            self.assertAlmostEqual(rhs[0], x_left)
            self.assertAlmostEqual(rhs[1], y_left)
            self.assertAlmostEqual(rhs[2], 0)
            self.assertAlmostEqual(rhs[3], 0)
            self.assertAlmostEqual(rhs[4], 0)
            self.assertAlmostEqual(rhs[5], moment_left)
            self.assertAlmostEqual(rhs[6], x_right)
            self.assertAlmostEqual(rhs[7], y_right)
            self.assertAlmostEqual(rhs[8], 0)
            self.assertAlmostEqual(rhs[9], 0)
            self.assertAlmostEqual(rhs[10], 0)
            self.assertAlmostEqual(rhs[11], moment_right)

        # set load outside condition
        cond.SetValue(StructuralMechanicsApplication.MOVING_LOAD_LOCAL_DISTANCE, length * 2)
        cond.CalculateLocalSystem(lhs, rhs, mp.ProcessInfo)
        self.assertAlmostEqual(rhs[0], 0)
        self.assertAlmostEqual(rhs[1], 0)
        self.assertAlmostEqual(rhs[2], 0)
        self.assertAlmostEqual(rhs[3], 0)
        self.assertAlmostEqual(rhs[4], 0)
        self.assertAlmostEqual(rhs[5], 0)
        self.assertAlmostEqual(rhs[6], 0)
        self.assertAlmostEqual(rhs[7], 0)
        self.assertAlmostEqual(rhs[8], 0)
        self.assertAlmostEqual(rhs[9], 0)
        self.assertAlmostEqual(rhs[10], 0)
        self.assertAlmostEqual(rhs[11], 0)

    def _MovingLoadCondition3D2NRotDofY(self):
        """"
        Check moving load on a 3d 2n line element, where the line element is inclined around y axis.
        """

        current_model = KratosMultiphysics.Model()

        mp = current_model.CreateModelPart("solid_part")
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.ROTATION)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION_MOMENT)

        # create nodes
        second_coord = [math.sqrt(2),  0.0, math.sqrt(2.0)]
        mp.CreateNewNode(1, 0.0, 0.0, 0.0)
        mp.CreateNewNode(2, second_coord[0], second_coord[1], second_coord[2])

        length = 2.0
        rotation = math.atan(second_coord[2] / second_coord[0])

        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.REACTION_X, mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Y, KratosMultiphysics.REACTION_Y, mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Z, KratosMultiphysics.REACTION_Z, mp)

        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ROTATION_X, KratosMultiphysics.REACTION_MOMENT_X, mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ROTATION_Y, KratosMultiphysics.REACTION_MOMENT_Y, mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ROTATION_Z, KratosMultiphysics.REACTION_MOMENT_Z, mp)

        # create condition
        cond = mp.CreateNewCondition("MovingLoadCondition3D2N", 1, [1, 2], mp.GetProperties()[1])

        # initialise matrices
        lhs = KratosMultiphysics.Matrix(0, 0)
        rhs = KratosMultiphysics.Vector(0)

        # set perpendicular POINT_LOAD to the condition
        load_on_cond = KratosMultiphysics.Vector(3)
        load_on_cond[0] = 2.0
        load_on_cond[1] = 0.0
        load_on_cond[2] = -2.0

        cond.SetValue(StructuralMechanicsApplication.POINT_LOAD, load_on_cond)

        # compare results when load is located at different points
        num_tests = 5
        for i in range(num_tests):
            distance = i/(num_tests-1) * length
            (x_left, z_left, moment_left), (x_right, z_right, moment_right) = \
                self.__CalculateReactionWithRotation(length, rotation, [load_on_cond[0], load_on_cond[2],
                                                                      load_on_cond[1]], distance)

            # set load on quarter
            cond.SetValue(StructuralMechanicsApplication.MOVING_LOAD_LOCAL_DISTANCE, distance)
            cond.CalculateLocalSystem(lhs, rhs, mp.ProcessInfo)

            self.assertAlmostEqual(rhs[0], x_left)
            self.assertAlmostEqual(rhs[1], 0)
            self.assertAlmostEqual(rhs[2], z_left)
            self.assertAlmostEqual(rhs[3], 0)
            self.assertAlmostEqual(rhs[4], moment_left)
            self.assertAlmostEqual(rhs[5], 0)
            self.assertAlmostEqual(rhs[6], x_right)
            self.assertAlmostEqual(rhs[7], 0)
            self.assertAlmostEqual(rhs[8], z_right)
            self.assertAlmostEqual(rhs[9], 0)
            self.assertAlmostEqual(rhs[10], moment_right)
            self.assertAlmostEqual(rhs[11], 0)

    def _MovingLoadCondition2D2N(self):
        current_model = KratosMultiphysics.Model()
        mp = current_model.CreateModelPart("solid_part")
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)

        # create nodes
        second_coord = [math.sqrt(2), math.sqrt(2.0), 0.0]
        mp.CreateNewNode(1,0.0,0.0,0.0)
        mp.CreateNewNode(2,second_coord[0],second_coord[1],second_coord[2])
        length = 2.0
        rotation = math.atan(second_coord[1] / second_coord[0])

        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.REACTION_X,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Y, KratosMultiphysics.REACTION_Y,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Z, KratosMultiphysics.REACTION_Z,mp)

        # create condition
        cond = mp.CreateNewCondition("MovingLoadCondition2D2N", 1, [1,2], mp.GetProperties()[1])

        # initialise matrices
        lhs = KratosMultiphysics.Matrix(0,0)
        rhs = KratosMultiphysics.Vector(0)

        # set perpendicular POINT_LOAD to the condition
        load_on_cond = KratosMultiphysics.Vector(3)
        load_on_cond[0] = 2.0
        load_on_cond[1] = -2.0
        load_on_cond[2] = 0.0 # note that this is a 2D condition

        num_tests = 5
        for i in range(num_tests):
            distance = i / (num_tests - 1) * length

            (x_left, y_left), (x_right, y_right) = self.__CalculateReaction2n(length, rotation, load_on_cond, distance)

            cond.SetValue(StructuralMechanicsApplication.POINT_LOAD,load_on_cond)

            # set load on node
            cond.SetValue(StructuralMechanicsApplication.MOVING_LOAD_LOCAL_DISTANCE, distance)
            cond.CalculateLocalSystem(lhs, rhs, mp.ProcessInfo)

            self.assertAlmostEqual(rhs[0], x_left)
            self.assertAlmostEqual(rhs[1], y_left)
            self.assertAlmostEqual(rhs[2], x_right)
            self.assertAlmostEqual(rhs[3], y_right)

        # set load outside condition
        cond.SetValue(StructuralMechanicsApplication.MOVING_LOAD_LOCAL_DISTANCE, length * 2)
        cond.CalculateLocalSystem(lhs, rhs, mp.ProcessInfo)
        self.assertAlmostEqual(rhs[0], 0)
        self.assertAlmostEqual(rhs[1], 0)
        self.assertAlmostEqual(rhs[2], 0)
        self.assertAlmostEqual(rhs[3], 0)

    def _MovingLoadCondition2D2NRot(self):
        current_model = KratosMultiphysics.Model()
        mp = current_model.CreateModelPart("solid_part")
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.ROTATION)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION_MOMENT)

        # create nodes
        second_coord = [math.sqrt(2), math.sqrt(2.0), 0.0]
        mp.CreateNewNode(1, 0.0, 0.0, 0.0)
        mp.CreateNewNode(2, math.sqrt(2), math.sqrt(2.0), 0.0)
        length = 2.0
        rotation = math.atan(second_coord[2] / second_coord[0])

        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.REACTION_X, mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Y, KratosMultiphysics.REACTION_Y, mp)

        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ROTATION_X, KratosMultiphysics.REACTION_MOMENT_X, mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ROTATION_Y, KratosMultiphysics.REACTION_MOMENT_Y, mp)

        # create condition
        cond = mp.CreateNewCondition("MovingLoadCondition2D2N", 1, [1, 2], mp.GetProperties()[1])

        # initialise matrices
        lhs = KratosMultiphysics.Matrix(0, 0)
        rhs = KratosMultiphysics.Vector(0)

        # set perpendicular POINT_LOAD to the condition
        load_on_cond = KratosMultiphysics.Vector(3)
        load_on_cond[0] = 2.0
        load_on_cond[1] = -2.0
        load_on_cond[2] = 0.0  # note that this is a 2D condition

        cond.SetValue(StructuralMechanicsApplication.POINT_LOAD, load_on_cond)

        num_tests = 5
        for i in range(num_tests):
            distance = i / (num_tests - 1) * length
            (x_left, y_left, moment_left), (x_right, y_right, moment_right) = \
                self.__CalculateReactionWithRotation(length, rotation, [load_on_cond[0], load_on_cond[2],
                                                                        load_on_cond[1]], distance)

            # set load on node
            cond.SetValue(StructuralMechanicsApplication.MOVING_LOAD_LOCAL_DISTANCE, distance)
            cond.CalculateLocalSystem(lhs, rhs, mp.ProcessInfo)

            self.assertAlmostEqual(rhs[0], x_left)
            self.assertAlmostEqual(rhs[1], y_left)
            self.assertAlmostEqual(rhs[2], moment_left)
            self.assertAlmostEqual(rhs[3], x_right)
            self.assertAlmostEqual(rhs[4], y_right)
            self.assertAlmostEqual(rhs[5], moment_right)

    def _MovingLoadCondition2D3N(self):
        current_model = KratosMultiphysics.Model()
        mp = current_model.CreateModelPart("solid_part")
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)

        #create nodes
        second_coord = [math.sqrt(2), math.sqrt(2.0), 0.0]
        third_coord = [math.sqrt(2)/2, math.sqrt(2.0)/2, 0.0]
        mp.CreateNewNode(1,0.0,0.0,0.0)
        mp.CreateNewNode(2,second_coord[0],second_coord[1],second_coord[2])
        mp.CreateNewNode(3, third_coord[0], third_coord[1], third_coord[2])
        length = 2.0
        rotation = math.atan(second_coord[1] / second_coord[0])

        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.REACTION_X,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Y, KratosMultiphysics.REACTION_Y,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Z, KratosMultiphysics.REACTION_Z,mp)

        # create condition
        cond = mp.CreateNewCondition("MovingLoadCondition2D3N", 1, [1,2,3], mp.GetProperties()[1])

        # initialise matrices
        lhs = KratosMultiphysics.Matrix(0,0)
        rhs = KratosMultiphysics.Vector(0)

        # set perpendicular POINT_LOAD to the condition
        load_on_cond = KratosMultiphysics.Vector(3)
        load_on_cond[0] = 2.0
        load_on_cond[1] = -3.0
        load_on_cond[2] = 0.0 # note that this is a 2D condition

        num_tests = 5
        for i in range(num_tests):
            distance = i / (num_tests - 1) * length

            (x_left, y_left),(x_mid, y_mid), (x_right, y_right) = \
                self.__CalculateReaction3n(length, rotation, load_on_cond, distance)

            cond.SetValue(StructuralMechanicsApplication.POINT_LOAD,load_on_cond)

            # set load on node
            cond.SetValue(StructuralMechanicsApplication.MOVING_LOAD_LOCAL_DISTANCE, distance)
            cond.CalculateLocalSystem(lhs, rhs, mp.ProcessInfo)

            self.assertAlmostEqual(rhs[0], x_left)
            self.assertAlmostEqual(rhs[1], y_left)
            self.assertAlmostEqual(rhs[2], x_right)
            self.assertAlmostEqual(rhs[3], y_right)
            self.assertAlmostEqual(rhs[4], x_mid)
            self.assertAlmostEqual(rhs[5], y_mid)

        # set load outside condition
        cond.SetValue(StructuralMechanicsApplication.MOVING_LOAD_LOCAL_DISTANCE, length * 2)
        cond.CalculateLocalSystem(lhs, rhs, mp.ProcessInfo)
        self.assertAlmostEqual(rhs[0], 0)
        self.assertAlmostEqual(rhs[1], 0)
        self.assertAlmostEqual(rhs[2], 0)
        self.assertAlmostEqual(rhs[3], 0)
        self.assertAlmostEqual(rhs[4], 0)
        self.assertAlmostEqual(rhs[5], 0)

    def _MovingLoadCondition3D3N(self):
        current_model = KratosMultiphysics.Model()
        mp = current_model.CreateModelPart("solid_part")
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)

        # create nodes
        second_coord = [math.sqrt(2), 0.0, math.sqrt(2.0)]
        third_coord = [math.sqrt(2)/2, 0.0, math.sqrt(2.0)/2]
        mp.CreateNewNode(1,0.0,0.0,0.0)
        mp.CreateNewNode(2,second_coord[0],second_coord[1],second_coord[2])
        mp.CreateNewNode(3, third_coord[0], third_coord[1], third_coord[2])
        length = 2.0
        rotation = math.atan(second_coord[2] / second_coord[0])

        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.REACTION_X,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Y, KratosMultiphysics.REACTION_Y,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Z, KratosMultiphysics.REACTION_Z,mp)

        # create condition
        cond = mp.CreateNewCondition("MovingLoadCondition3D3N", 1, [1,2,3], mp.GetProperties()[1])

        # initialise matrices
        lhs = KratosMultiphysics.Matrix(0,0)
        rhs = KratosMultiphysics.Vector(0)

        # set perpendicular POINT_LOAD to the condition
        load_on_cond = KratosMultiphysics.Vector(3)
        load_on_cond[0] = 2.0
        load_on_cond[1] = 0.0
        load_on_cond[2] = -3.0

        num_tests = 5
        for i in range(num_tests):
            distance = i / (num_tests - 1) * length

            (x_left, z_left), (x_mid, z_mid), (x_right, z_right) = \
                self.__CalculateReaction3n(length, rotation, [load_on_cond[0], load_on_cond[2], load_on_cond[1]],
                                           distance)

            cond.SetValue(StructuralMechanicsApplication.POINT_LOAD,load_on_cond)

            # set load on node
            cond.SetValue(StructuralMechanicsApplication.MOVING_LOAD_LOCAL_DISTANCE, distance)
            cond.CalculateLocalSystem(lhs, rhs, mp.ProcessInfo)

            self.assertAlmostEqual(rhs[0], x_left)
            self.assertAlmostEqual(rhs[1], 0)
            self.assertAlmostEqual(rhs[2], z_left)
            self.assertAlmostEqual(rhs[3], x_right)
            self.assertAlmostEqual(rhs[4], 0)
            self.assertAlmostEqual(rhs[5], z_right)
            self.assertAlmostEqual(rhs[6], x_mid)
            self.assertAlmostEqual(rhs[7], 0)
            self.assertAlmostEqual(rhs[8], z_mid)

        # set load outside condition
        cond.SetValue(StructuralMechanicsApplication.MOVING_LOAD_LOCAL_DISTANCE, length * 2)
        cond.CalculateLocalSystem(lhs, rhs, mp.ProcessInfo)
        self.assertAlmostEqual(rhs[0], 0)
        self.assertAlmostEqual(rhs[1], 0)
        self.assertAlmostEqual(rhs[2], 0)
        self.assertAlmostEqual(rhs[3], 0)
        self.assertAlmostEqual(rhs[4], 0)
        self.assertAlmostEqual(rhs[5], 0)
        self.assertAlmostEqual(rhs[6], 0)
        self.assertAlmostEqual(rhs[7], 0)
        self.assertAlmostEqual(rhs[8], 0)

    def test_MovingLoadCondition3D2NRotDofZ(self):
        self._MovingLoadCondition3D2NRotDofZ()

    def test_MovingLoadCondition3D2NRotDofY(self):
        self._MovingLoadCondition3D2NRotDofY()

    def test_MovingLoadCondition2D2N(self):
        self._MovingLoadCondition2D2N()

    def test_MovingLoadCondition2D3N(self):
        self._MovingLoadCondition2D3N()

    def test_MovingLoadCondition3D3N(self):
        self._MovingLoadCondition3D3N()


if __name__ == '__main__':
    KratosUnittest.main()

import KratosMultiphysics

import KratosMultiphysics.StructuralMechanicsApplication as SMA
import KratosMultiphysics.KratosUnittest as KratosUnittest

class TestSetMovingLoadProcess(KratosUnittest.TestCase):

    def _TestSetMovingLoad(self):
        """
        Tests a moving load on 1 condition element, where the nodes of the element are sorted in the direction of the
        moving load
        Returns
        -------

        """

        current_model = KratosMultiphysics.Model()
        mp = current_model.CreateModelPart("solid_part")

        # create nodes
        second_coord = [1, 0, 0.0]
        mp.CreateNewNode(1,0.0,0.0,0.0)
        mp.CreateNewNode(2,second_coord[0],second_coord[1],0.0)

        # create condition
        cond = mp.CreateNewCondition("MovingLoadCondition2D2N", 1, [1,2], mp.GetProperties()[1])

        parameters = KratosMultiphysics.Parameters("""
                {
                    "help"            : "This process applies a moving load condition belonging to a modelpart. The load moves over line elements.",
                    "model_part_name" : "please_specify_model_part_name",
                    "variable_name"   : "POINT_LOAD",
                    "load"            : [0.0, -2.0, 0.0],
                    "direction"       : [1,1,1],
                    "velocity"        : 1,
                    "origin"          : [0,0,0]
                }
                """
                                                         )
        mp.ProcessInfo.SetValue(KratosMultiphysics.TIME, 0)
        mp.ProcessInfo.SetValue(KratosMultiphysics.DELTA_TIME, 0.25)

        process = SMA.SetMovingLoadProcess(mp, parameters)

        # initialise and set load
        process.ExecuteInitialize()
        process.ExecuteInitializeSolutionStep()

        # initialise matrices
        lhs = KratosMultiphysics.Matrix(0, 0)
        rhs = KratosMultiphysics.Vector(0)

        # set load on node
        cond.CalculateLocalSystem(lhs, rhs, mp.ProcessInfo)

        self.assertAlmostEqual(rhs[0], 0)
        self.assertAlmostEqual(rhs[1], -2)
        self.assertAlmostEqual(rhs[2], 0)
        self.assertAlmostEqual(rhs[3], 0)

        # move load
        process.ExecuteFinalizeSolutionStep()
        process.ExecuteInitializeSolutionStep()

        # check if interpolation is done correctly
        cond.CalculateLocalSystem(lhs, rhs, mp.ProcessInfo)

        self.assertAlmostEqual(rhs[0], 0)
        self.assertAlmostEqual(rhs[1], -1.5)
        self.assertAlmostEqual(rhs[2], 0)
        self.assertAlmostEqual(rhs[3], -0.5)

    def _TestSetMovingLoadReverseGeom(self):
        """
        Tests a moving load on 1 condition element, where the nodes of the element are reversed compared to the
        direction of the moving load
        Returns
        -------

        """

        current_model = KratosMultiphysics.Model()
        mp = current_model.CreateModelPart("solid_part")

        # create nodes
        second_coord = [1, 0, 0.0]
        mp.CreateNewNode(1, 0.0, 0.0, 0.0)
        mp.CreateNewNode(2, second_coord[0], second_coord[1], 0.0)

        # create condition
        cond = mp.CreateNewCondition("MovingLoadCondition2D2N", 1, [2, 1], mp.GetProperties()[1])

        parameters = KratosMultiphysics.Parameters("""
                   {
                       "help"            : "This process applies a moving load condition belonging to a modelpart. The load moves over line elements.",
                       "model_part_name" : "please_specify_model_part_name",
                       "variable_name"   : "POINT_LOAD",
                       "load"            : [0.0, -2.0, 0.0],
                       "direction"       : [1,1,1],
                       "velocity"        : 1,
                       "origin"          : [0,0,0]
                   }
                   """
                                                   )
        mp.ProcessInfo.SetValue(KratosMultiphysics.TIME,
                                0)
        mp.ProcessInfo.SetValue(KratosMultiphysics.DELTA_TIME,
                                0.25)

        process = SMA.SetMovingLoadProcess(mp, parameters)

        # initialize and set load
        process.ExecuteInitialize()
        process.ExecuteInitializeSolutionStep()

        # initialise matrices
        lhs = KratosMultiphysics.Matrix(0, 0)
        rhs = KratosMultiphysics.Vector(0)

        # set load on node
        # cond.SetValue(SMA.MOVING_LOAD_LOCAL_DISTANCE, 0)
        cond.CalculateLocalSystem(lhs, rhs, mp.ProcessInfo)

        self.assertAlmostEqual(rhs[0], 0)
        self.assertAlmostEqual(rhs[1], 0)
        self.assertAlmostEqual(rhs[2], 0)
        self.assertAlmostEqual(rhs[3], -2)

        # move load
        process.ExecuteFinalizeSolutionStep()
        process.ExecuteInitializeSolutionStep()

        # check if interpolation is done correctly
        cond.CalculateLocalSystem(lhs, rhs, mp.ProcessInfo)

        self.assertAlmostEqual(rhs[0], 0)
        self.assertAlmostEqual(rhs[1], -0.5)
        self.assertAlmostEqual(rhs[2], 0)
        self.assertAlmostEqual(rhs[3], -1.5)

    def _TestSetMovingLoadMultipleConditions(self):
        """
        Tests a moving load on 2 condition elements, where the order of the elements is sorted in the direction of the
        moving load
        Returns
        -------

        """

        current_model = KratosMultiphysics.Model()
        mp = current_model.CreateModelPart("solid_part")

        #create nodes
        second_coord = [1.0, 0.0, 0.0]
        third_coord = [2.0, 0.0, 0.0]
        mp.CreateNewNode(1, 0.0, 0.0, 0.0)
        mp.CreateNewNode(2, second_coord[0],second_coord[1],second_coord[2])
        mp.CreateNewNode(3, third_coord[0], third_coord[1], third_coord[2])

        # create condition
        conditions=[]
        conditions.append(mp.CreateNewCondition("MovingLoadCondition2D2N", 1, [1, 2], mp.GetProperties()[1]))
        conditions.append(mp.CreateNewCondition("MovingLoadCondition2D2N", 2, [2, 3], mp.GetProperties()[1]))

        parameters = KratosMultiphysics.Parameters("""
                {
                    "help"            : "This process applies a moving load condition belonging to a modelpart. The load moves over line elements.",
                    "model_part_name" : "please_specify_model_part_name",
                    "variable_name"   : "POINT_LOAD",
                    "load"            : [0.0, -2.0, 0.0],
                    "direction"       : [1,1,1],
                    "velocity"        : 1,
                    "origin"          : [0,0,0]
                }
                """
                                                         )
        mp.ProcessInfo.SetValue(KratosMultiphysics.TIME,
                                                  0)
        mp.ProcessInfo.SetValue(KratosMultiphysics.DELTA_TIME,
                                                  0.5)
        process = SMA.SetMovingLoadProcess(mp,parameters)


        # initialize and set load
        process.ExecuteInitialize()
        process.ExecuteInitializeSolutionStep()

        # initialise matrices
        lhs = KratosMultiphysics.Matrix(0,0)
        rhs = KratosMultiphysics.Vector(0)

        # set load on node

        all_rhs = []
        for cond in conditions:
            cond.CalculateLocalSystem(lhs, rhs, mp.ProcessInfo)
            all_rhs.append(list(rhs))

        self.assertAlmostEqual(all_rhs[0][0], 0)
        self.assertAlmostEqual(all_rhs[0][1], -2)
        self.assertAlmostEqual(all_rhs[0][2], 0)
        self.assertAlmostEqual(all_rhs[0][3], 0)

        self.assertAlmostEqual(all_rhs[1][0], 0)
        self.assertAlmostEqual(all_rhs[1][1], 0)
        self.assertAlmostEqual(all_rhs[1][2], 0)
        self.assertAlmostEqual(all_rhs[1][3], 0)

        # move load within first element
        process.ExecuteFinalizeSolutionStep()
        process.ExecuteInitializeSolutionStep()

        # check if interpolation is done correctly
        all_rhs = []
        for cond in conditions:
            cond.CalculateLocalSystem(lhs, rhs, mp.ProcessInfo)
            all_rhs.append(list(rhs))

        self.assertAlmostEqual(all_rhs[0][0], 0)
        self.assertAlmostEqual(all_rhs[0][1], -1)
        self.assertAlmostEqual(all_rhs[0][2], 0)
        self.assertAlmostEqual(all_rhs[0][3], -1)

        self.assertAlmostEqual(all_rhs[1][0], 0)
        self.assertAlmostEqual(all_rhs[1][1], 0)
        self.assertAlmostEqual(all_rhs[1][2], 0)
        self.assertAlmostEqual(all_rhs[1][3], 0)

        # move load to element connection element
        process.ExecuteFinalizeSolutionStep()
        process.ExecuteInitializeSolutionStep()

        all_rhs = []
        for cond in conditions:
            cond.CalculateLocalSystem(lhs, rhs, mp.ProcessInfo)
            all_rhs.append(list(rhs))

        self.assertAlmostEqual(all_rhs[0][0], 0)
        self.assertAlmostEqual(all_rhs[0][1], 0)
        self.assertAlmostEqual(all_rhs[0][2], 0)
        self.assertAlmostEqual(all_rhs[0][3], -2)

        self.assertAlmostEqual(all_rhs[1][0], 0)
        self.assertAlmostEqual(all_rhs[1][1], 0)
        self.assertAlmostEqual(all_rhs[1][2], 0)
        self.assertAlmostEqual(all_rhs[1][3], 0)

        # move load to next element
        process.ExecuteFinalizeSolutionStep()
        process.ExecuteInitializeSolutionStep()

        all_rhs = []
        for cond in conditions:
            cond.CalculateLocalSystem(lhs, rhs, mp.ProcessInfo)
            all_rhs.append(list(rhs))

        self.assertAlmostEqual(all_rhs[0][0], 0)
        self.assertAlmostEqual(all_rhs[0][1], 0)
        self.assertAlmostEqual(all_rhs[0][2], 0)
        self.assertAlmostEqual(all_rhs[0][3], 0)

        self.assertAlmostEqual(all_rhs[1][0], 0)
        self.assertAlmostEqual(all_rhs[1][1], -1)
        self.assertAlmostEqual(all_rhs[1][2], 0)
        self.assertAlmostEqual(all_rhs[1][3], -1)


    def _TestSetMovingLoadMultipleConditionsReversed(self):
        """
        Tests a moving load on 2 condition elements, where the order of the elements is reversed compared to the moving
        direction of the load
        Returns
        -------

        """

        current_model = KratosMultiphysics.Model()
        mp = current_model.CreateModelPart("solid_part")

        # create nodes
        second_coord = [1.0, 0.0, 0.0]
        third_coord = [2.0, 0.0, 0.0]
        mp.CreateNewNode(1, 0.0, 0.0, 0.0)
        mp.CreateNewNode(2, second_coord[0],second_coord[1],second_coord[2])
        mp.CreateNewNode(3, third_coord[0], third_coord[1], third_coord[2])

        # create condition
        conditions = []
        conditions.append(mp.CreateNewCondition("MovingLoadCondition2D2N", 1, [3, 2], mp.GetProperties()[1]))
        conditions.append(mp.CreateNewCondition("MovingLoadCondition2D2N", 2, [2, 1], mp.GetProperties()[1]))

        # set parameters and process info
        parameters = KratosMultiphysics.Parameters("""
                {
                    "help"            : "This process applies a moving load condition belonging to a modelpart. The load moves over line elements.",
                    "model_part_name" : "please_specify_model_part_name",
                    "variable_name"   : "POINT_LOAD",
                    "load"            : [0.0, -2.0, 0.0],
                    "direction"       : [1,1,1],
                    "velocity"        : 1,
                    "origin"          : [0,0,0]
                }
                """
                                                         )
        mp.ProcessInfo.SetValue(KratosMultiphysics.TIME, 0)
        mp.ProcessInfo.SetValue(KratosMultiphysics.DELTA_TIME, 0.5)
        process = SMA.SetMovingLoadProcess(mp,parameters)

        # initialize and set load
        process.ExecuteInitialize()
        process.ExecuteInitializeSolutionStep()

        # initialise matrices
        lhs = KratosMultiphysics.Matrix(0,0)
        rhs = KratosMultiphysics.Vector(0)

        # calculate load
        all_rhs = []
        for cond in conditions:
            cond.CalculateLocalSystem(lhs, rhs, mp.ProcessInfo)
            all_rhs.append(list(rhs))

        self.assertAlmostEqual(all_rhs[0][0], 0)
        self.assertAlmostEqual(all_rhs[0][1], 0)
        self.assertAlmostEqual(all_rhs[0][2], 0)
        self.assertAlmostEqual(all_rhs[0][3], 0)

        self.assertAlmostEqual(all_rhs[1][0], 0)
        self.assertAlmostEqual(all_rhs[1][1], 0)
        self.assertAlmostEqual(all_rhs[1][2], 0)
        self.assertAlmostEqual(all_rhs[1][3], -2)

        # move load within first element
        process.ExecuteFinalizeSolutionStep()
        process.ExecuteInitializeSolutionStep()

        # check if interpolation is done correctly
        all_rhs = []
        for cond in conditions:
            cond.CalculateLocalSystem(lhs, rhs, mp.ProcessInfo)
            all_rhs.append(list(rhs))

        self.assertAlmostEqual(all_rhs[0][0], 0)
        self.assertAlmostEqual(all_rhs[0][1], 0)
        self.assertAlmostEqual(all_rhs[0][2], 0)
        self.assertAlmostEqual(all_rhs[0][3], 0)

        self.assertAlmostEqual(all_rhs[1][0], 0)
        self.assertAlmostEqual(all_rhs[1][1], -1)
        self.assertAlmostEqual(all_rhs[1][2], 0)
        self.assertAlmostEqual(all_rhs[1][3], -1)

        # move load to element connection element
        process.ExecuteFinalizeSolutionStep()
        process.ExecuteInitializeSolutionStep()

        # calculate load
        all_rhs = []
        for cond in conditions:
            cond.CalculateLocalSystem(lhs, rhs, mp.ProcessInfo)
            all_rhs.append(list(rhs))

        self.assertAlmostEqual(all_rhs[0][0], 0)
        self.assertAlmostEqual(all_rhs[0][1], 0)
        self.assertAlmostEqual(all_rhs[0][2], 0)
        self.assertAlmostEqual(all_rhs[0][3], -2)

        self.assertAlmostEqual(all_rhs[1][0], 0)
        self.assertAlmostEqual(all_rhs[1][1], 0)
        self.assertAlmostEqual(all_rhs[1][2], 0)
        self.assertAlmostEqual(all_rhs[1][3], 0)

        # move load to next element, also increase time step
        mp.ProcessInfo.SetValue(KratosMultiphysics.DELTA_TIME, 0.75)
        process.ExecuteFinalizeSolutionStep()
        process.ExecuteInitializeSolutionStep()

        # calculate load
        all_rhs = []
        for cond in conditions:
            cond.CalculateLocalSystem(lhs, rhs, mp.ProcessInfo)
            all_rhs.append(list(rhs))

        self.assertAlmostEqual(all_rhs[0][0], 0)
        self.assertAlmostEqual(all_rhs[0][1], -1.5)
        self.assertAlmostEqual(all_rhs[0][2], 0)
        self.assertAlmostEqual(all_rhs[0][3], -0.5)

        self.assertAlmostEqual(all_rhs[1][0], 0)
        self.assertAlmostEqual(all_rhs[1][1], 0)
        self.assertAlmostEqual(all_rhs[1][2], 0)
        self.assertAlmostEqual(all_rhs[1][3], 0)

    def _TestSetMovingLoadMultipleConditionsDifferentOrigin(self):
        """
        Tests a moving load on 2 condition elements, where the order of the elements is sorted in the direction of the
        moving load
        Returns
        -------

        """

        current_model = KratosMultiphysics.Model()
        mp = current_model.CreateModelPart("solid_part")

        #create nodes
        second_coord = [1.0, 0.0, 0.0]
        third_coord = [2.0, 0.0, 0.0]
        mp.CreateNewNode(1, 0.0, 0.0, 0.0)
        mp.CreateNewNode(2, second_coord[0],second_coord[1],second_coord[2])
        mp.CreateNewNode(3, third_coord[0], third_coord[1], third_coord[2])

        # create condition
        conditions=[]
        conditions.append(mp.CreateNewCondition("MovingLoadCondition2D2N", 1, [1, 2], mp.GetProperties()[1]))
        conditions.append(mp.CreateNewCondition("MovingLoadCondition2D2N", 2, [2, 3], mp.GetProperties()[1]))

        parameters = KratosMultiphysics.Parameters("""
                {
                    "help"            : "This process applies a moving load condition belonging to a modelpart. The load moves over line elements.",
                    "model_part_name" : "please_specify_model_part_name",
                    "variable_name"   : "POINT_LOAD",
                    "load"            : [0.0, -2.0, 0.0],
                    "direction"       : [1,1,1],
                    "velocity"        : 1,
                    "origin"          : [1.25,0,0]
                }
                """
                                                         )
        process = SMA.SetMovingLoadProcess(mp,parameters)

        # initialize and set load
        process.ExecuteInitialize()
        process.ExecuteInitializeSolutionStep()

        # initialise matrices
        lhs = KratosMultiphysics.Matrix(0,0)
        rhs = KratosMultiphysics.Vector(0)

        # set load on node
        all_rhs = []
        for cond in conditions:
            cond.CalculateLocalSystem(lhs, rhs, mp.ProcessInfo)
            all_rhs.append(list(rhs))

        self.assertAlmostEqual(all_rhs[0][0], 0)
        self.assertAlmostEqual(all_rhs[0][1], 0)
        self.assertAlmostEqual(all_rhs[0][2], 0)
        self.assertAlmostEqual(all_rhs[0][3], 0)

        self.assertAlmostEqual(all_rhs[1][0], 0)
        self.assertAlmostEqual(all_rhs[1][1], -1.5)
        self.assertAlmostEqual(all_rhs[1][2], 0)
        self.assertAlmostEqual(all_rhs[1][3], -0.5)


    def _TestSetMovingLoadMultipleConditionsDifferentOriginReversed(self):
        """
        Tests a moving load on 2 condition elements, where the order of the elements is sorted in the direction of the
        moving load
        Returns
        -------

        """

        current_model = KratosMultiphysics.Model()
        mp = current_model.CreateModelPart("solid_part")

        # create nodes
        second_coord = [1.0, 0.0, 0.0]
        third_coord = [2.0, 0.0, 0.0]
        mp.CreateNewNode(1, 0.0, 0.0, 0.0)
        mp.CreateNewNode(2, second_coord[0],second_coord[1],second_coord[2])
        mp.CreateNewNode(3, third_coord[0], third_coord[1], third_coord[2])

        # create condition
        conditions=[]
        conditions.append(mp.CreateNewCondition("MovingLoadCondition2D2N", 1, [3, 2], mp.GetProperties()[1]))
        conditions.append(mp.CreateNewCondition("MovingLoadCondition2D2N", 2, [2, 1], mp.GetProperties()[1]))

        parameters = KratosMultiphysics.Parameters("""
                {
                    "help"            : "This process applies a moving load condition belonging to a modelpart. The load moves over line elements.",
                    "model_part_name" : "please_specify_model_part_name",
                    "variable_name"   : "POINT_LOAD",
                    "load"            : [0.0, -2.0, 0.0],
                    "direction"       : [1,1,1],
                    "velocity"        : 1,
                    "origin"          : [1.25,0,0]
                }
                """
                                                         )
        process = SMA.SetMovingLoadProcess(mp,parameters)

        # initialize and set load
        process.ExecuteInitialize()
        process.ExecuteInitializeSolutionStep()

        # initialise matrices
        lhs = KratosMultiphysics.Matrix(0,0)
        rhs = KratosMultiphysics.Vector(0)

        # set load on node
        all_rhs = []
        for cond in conditions:
            cond.CalculateLocalSystem(lhs, rhs, mp.ProcessInfo)
            all_rhs.append(list(rhs))

        self.assertAlmostEqual(all_rhs[0][0], 0)
        self.assertAlmostEqual(all_rhs[0][1], -0.5)
        self.assertAlmostEqual(all_rhs[0][2], 0)
        self.assertAlmostEqual(all_rhs[0][3], -1.5)

        self.assertAlmostEqual(all_rhs[1][0], 0)
        self.assertAlmostEqual(all_rhs[1][1], 0)
        self.assertAlmostEqual(all_rhs[1][2], 0)
        self.assertAlmostEqual(all_rhs[1][3], 0)

    def _TestSetMovingLoadWithLoadFunction(self):
        """
       Tests a moving load on a condition element, where the load is a function of time
       Returns
       -------

       """

        current_model = KratosMultiphysics.Model()
        mp = current_model.CreateModelPart("solid_part")

        # create nodes
        second_coord = [1, 0, 0.0]
        mp.CreateNewNode(1,0.0,0.0,0.0)
        mp.CreateNewNode(2,second_coord[0],second_coord[1],0.0)

        # create condition
        cond = mp.CreateNewCondition("MovingLoadCondition2D2N", 1, [1,2], mp.GetProperties()[1])

        parameters = KratosMultiphysics.Parameters("""
                {
                    "help"            : "This process applies a moving load condition belonging to a modelpart. The load moves over line elements.",
                    "model_part_name" : "please_specify_model_part_name",
                    "variable_name"   : "POINT_LOAD",
                    "load"            : ["0.0", "-2.0*t", "0.0"],
                    "direction"       : [1,1,1],
                    "velocity"        : 1,
                    "origin"          : [0.5,0,0]
                }
                """
                                                         )
        mp.ProcessInfo.SetValue(KratosMultiphysics.TIME, 0)
        mp.ProcessInfo.SetValue(KratosMultiphysics.DELTA_TIME, 0.25)

        process = SMA.SetMovingLoadProcess(mp, parameters)

        # initialise and set load
        process.ExecuteInitialize()
        process.ExecuteInitializeSolutionStep()

        # initialise matrices
        lhs = KratosMultiphysics.Matrix(0, 0)
        rhs = KratosMultiphysics.Vector(0)

        # set load on node
        cond.CalculateLocalSystem(lhs, rhs, mp.ProcessInfo)

        self.assertAlmostEqual(rhs[0], 0)
        self.assertAlmostEqual(rhs[1], 0)
        self.assertAlmostEqual(rhs[2], 0)
        self.assertAlmostEqual(rhs[3], 0)

        # change time and recalculate load
        mp.ProcessInfo.SetValue(KratosMultiphysics.TIME, 0.5)
        process.ExecuteInitializeSolutionStep()

        # check if interpolation is done correctly
        cond.CalculateLocalSystem(lhs, rhs, mp.ProcessInfo)

        self.assertAlmostEqual(rhs[0], 0)
        self.assertAlmostEqual(rhs[1], -0.5)
        self.assertAlmostEqual(rhs[2], 0)
        self.assertAlmostEqual(rhs[3], -0.5)


    def _TestSetMovingLoadWithVelocityFunction(self):
        """
       Tests a moving load on a condition element, where the load velocity is a function of time.

       Returns
       -------

       """

        current_model = KratosMultiphysics.Model()
        mp = current_model.CreateModelPart("solid_part")

        # create nodes
        second_coord = [1, 0, 0.0]
        mp.CreateNewNode(1,0.0,0.0,0.0)
        mp.CreateNewNode(2,second_coord[0],second_coord[1],0.0)

        # create condition
        cond = mp.CreateNewCondition("MovingLoadCondition2D2N", 1, [1,2], mp.GetProperties()[1])

        parameters = KratosMultiphysics.Parameters("""
                {
                    "help"            : "This process applies a moving load condition belonging to a modelpart. The load moves over line elements.",
                    "model_part_name" : "please_specify_model_part_name",
                    "variable_name"   : "POINT_LOAD",
                    "load"            : [0.0, -2.0, 0.0],
                    "direction"       : [1,1,1],
                    "velocity"        : "2.0*t",
                    "origin"          : [0.0,0,0]
                }
                """
                                                         )
        mp.ProcessInfo.SetValue(KratosMultiphysics.TIME, 0)
        mp.ProcessInfo.SetValue(KratosMultiphysics.DELTA_TIME, 0.25)

        process = SMA.SetMovingLoadProcess(mp, parameters)

        # initialise and set load
        process.ExecuteInitialize()
        process.ExecuteInitializeSolutionStep()

        # initialise matrices
        lhs = KratosMultiphysics.Matrix(0, 0)
        rhs = KratosMultiphysics.Vector(0)

        # set load on node
        cond.CalculateLocalSystem(lhs, rhs, mp.ProcessInfo)

        self.assertAlmostEqual(rhs[0], 0)
        self.assertAlmostEqual(rhs[1], -2)
        self.assertAlmostEqual(rhs[2], 0)
        self.assertAlmostEqual(rhs[3], 0)

        # change time and recalculate load
        process.ExecuteFinalizeSolutionStep()
        process.ExecuteInitializeSolutionStep()

        # check if interpolation is done correctly
        cond.CalculateLocalSystem(lhs, rhs, mp.ProcessInfo)

        self.assertAlmostEqual(rhs[0], 0)
        self.assertAlmostEqual(rhs[1], -2)
        self.assertAlmostEqual(rhs[2], 0)
        self.assertAlmostEqual(rhs[3], 0)

        mp.ProcessInfo.SetValue(KratosMultiphysics.TIME, 0.5)

        process.ExecuteFinalizeSolutionStep()
        process.ExecuteInitializeSolutionStep()

        # check if interpolation is done correctly
        cond.CalculateLocalSystem(lhs, rhs, mp.ProcessInfo)

        self.assertAlmostEqual(rhs[0], 0)
        self.assertAlmostEqual(rhs[1], -1.5)
        self.assertAlmostEqual(rhs[2], 0)
        self.assertAlmostEqual(rhs[3], -0.5)


    def test_SetMovingLoad(self):
        self._TestSetMovingLoad()

    def test_SetMovingLoadReverseGeom(self):
        self._TestSetMovingLoadReverseGeom()

    def test_SetMovingLoadMultipleConditions(self):
        self._TestSetMovingLoadMultipleConditions()

    def test_SetMovingLoadMultipleConditionsReversed(self):
        self._TestSetMovingLoadMultipleConditionsReversed()

    def test_SetMovingLoadMultipleConditionsDifferentOrigin(self):
        self._TestSetMovingLoadMultipleConditionsDifferentOrigin()

    def test_SetMovingLoadMultipleConditionsDifferentOriginReversed(self):
        self._TestSetMovingLoadMultipleConditionsDifferentOriginReversed()

    def test_SetMovingLoadWithLoadFunction(self):
        self._TestSetMovingLoadWithLoadFunction()

    def test_SetMovingLoadWithVelocityFunction(self):
        self._TestSetMovingLoadWithVelocityFunction()


if __name__ == '__main__':
    KratosUnittest.main()

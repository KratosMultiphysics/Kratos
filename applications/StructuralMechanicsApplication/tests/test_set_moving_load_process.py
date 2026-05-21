import KratosMultiphysics

import KratosMultiphysics.StructuralMechanicsApplication as SMA
import KratosMultiphysics.KratosUnittest as KratosUnittest

class TestSetMovingLoadProcess(KratosUnittest.TestCase):

    @staticmethod
    def setup_strategy(mp):
        """
        Setup default strategy for solving the problem

        Parameters
        ----------
        mp: model part

        Returns default strategy
        -------

        """

        linear_solver = KratosMultiphysics.SkylineLUFactorizationSolver()
        builder_and_solver = KratosMultiphysics.ResidualBasedBlockBuilderAndSolver(linear_solver)
        scheme = KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme()
        convergence_criterion = KratosMultiphysics.ResidualCriteria(1e-1, 1e-1)
        convergence_criterion.SetEchoLevel(0)

        max_iters = 10
        compute_reactions = False
        reform_step_dofs = False
        move_mesh_flag = False
        strategy = KratosMultiphysics.ResidualBasedNewtonRaphsonStrategy(mp,
                                                                         scheme,
                                                                         convergence_criterion,
                                                                         builder_and_solver,
                                                                         max_iters,
                                                                         compute_reactions,
                                                                         reform_step_dofs,
                                                                         move_mesh_flag)

        strategy.SetEchoLevel(0)

        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_X, mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Y, mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Z, mp)

        return strategy

    def checkRHS(self, rhs, expected_res):
        """
        routine to check calculation of rhs side within context of testing SetMovingLoad
        Returns
        -------

        """
        self.assertAlmostEqual(rhs[0], expected_res[0])
        self.assertAlmostEqual(rhs[1], expected_res[1])
        self.assertAlmostEqual(rhs[2], expected_res[2])
        self.assertAlmostEqual(rhs[3], expected_res[3])

    def _TestSetMovingLoad(self):
        """
        Tests a moving load on 1 condition element, where the nodes of the element are sorted in the direction of the
        moving load
        Returns
        -------

        """

        current_model = KratosMultiphysics.Model()
        mp = current_model.CreateModelPart("solid_part")
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)

        # create nodes
        second_coord = [1, 0, 0.0]
        mp.CreateNewNode(1, 0.0, 0.0, 0.0)
        mp.CreateNewNode(2, second_coord[0], second_coord[1], 0.0)

        strategy = self.setup_strategy(mp)

        # create condition
        cond = mp.CreateNewCondition("MovingLoadCondition2D2N", 1, [1, 2], mp.GetProperties()[1])

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
        strategy.InitializeSolutionStep()
        cond.CalculateLocalSystem(lhs, rhs, mp.ProcessInfo)

        self.checkRHS(rhs, [0.0, -2.0, 0.0, 0.0])

        # move load
        process.ExecuteFinalizeSolutionStep()
        process.ExecuteInitializeSolutionStep()

        # check if interpolation is done correctly
        cond.CalculateLocalSystem(lhs, rhs, mp.ProcessInfo)

        self.checkRHS(rhs, [0.0, -1.5, 0.0, -0.5])

    def _TestSetMovingLoadOffsetPositive(self):
        """
        Tests a moving load on 1 condition element, where the nodes of the element are sorted in the direction of the
        moving load, including a positive offset along line condition direction in velocity direction.
        Returns
        -------

        """

        current_model = KratosMultiphysics.Model()
        mp = current_model.CreateModelPart("solid_part")
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)

        # create nodes
        second_coord = [1, 0, 0.0]
        mp.CreateNewNode(1, 0.0, 0.0, 0.0)
        mp.CreateNewNode(2, second_coord[0], second_coord[1], 0.0)

        strategy = self.setup_strategy(mp)

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
                    "origin"          : [0,0,0],
                    "offset"          : 0.25
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
        strategy.InitializeSolutionStep()
        cond.CalculateLocalSystem(lhs, rhs, mp.ProcessInfo)

        self.checkRHS(rhs, [0.0, -1.5, 0.0, -0.5])

        # move load
        process.ExecuteFinalizeSolutionStep()
        process.ExecuteInitializeSolutionStep()

        # check if interpolation is done correctly
        cond.CalculateLocalSystem(lhs, rhs, mp.ProcessInfo)

        self.checkRHS(rhs, [0.0, -1.0, 0.0, -1.0])

    def _TestSetMovingLoadOffsetNegative(self):
        """
        Tests a moving load on 1 condition element, where the nodes of the element are sorted in the direction of the
        moving load, , including a negative offset along line condition direction in velocity direction.
        Returns
        -------

        """

        current_model = KratosMultiphysics.Model()
        mp = current_model.CreateModelPart("solid_part")
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)

        # create nodes
        second_coord = [1, 0, 0.0]
        mp.CreateNewNode(1,0.0,0.0,0.0)
        mp.CreateNewNode(2,second_coord[0],second_coord[1],0.0)

        strategy = self.setup_strategy(mp)

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
                    "origin"          : [0,0,0],
                    "offset"          : -0.25
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
        strategy.InitializeSolutionStep()
        cond.CalculateLocalSystem(lhs, rhs, mp.ProcessInfo)

        self.checkRHS(rhs, [0.0, 0.0, 0.0, 0.0])

        # move load
        process.ExecuteFinalizeSolutionStep()
        process.ExecuteInitializeSolutionStep()

        # check if interpolation is done correctly
        strategy.InitializeSolutionStep()
        cond.CalculateLocalSystem(lhs, rhs, mp.ProcessInfo)

        self.checkRHS(rhs, [0.0, -2.0, 0.0, 0.0])

    def _TestSetMovingLoadReverseGeom(self):
        """
        Tests a moving load on 1 condition element, where the nodes of the element are reversed compared to the
        direction of the moving load
        Returns
        -------

        """

        current_model = KratosMultiphysics.Model()
        mp = current_model.CreateModelPart("solid_part")
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)

        # create nodes
        second_coord = [1, 0, 0.0]
        mp.CreateNewNode(1, 0.0, 0.0, 0.0)
        mp.CreateNewNode(2, second_coord[0], second_coord[1], 0.0)

        strategy = self.setup_strategy(mp)

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
        strategy.InitializeSolutionStep()
        cond.CalculateLocalSystem(lhs, rhs, mp.ProcessInfo)

        self.checkRHS(rhs, [0.0, 0.0, 0.0, -2.0])

        # move load
        process.ExecuteFinalizeSolutionStep()
        process.ExecuteInitializeSolutionStep()

        # check if interpolation is done correctly
        cond.CalculateLocalSystem(lhs, rhs, mp.ProcessInfo)

        self.checkRHS(rhs, [0.0, -0.5, 0.0, -1.5])

    def _TestSetMovingLoadReverseGeomOffsetPositive(self):
        """
        Tests a moving load on 1 condition element, where the nodes of the element are reversed compared to the
        direction of the moving load, including a positive offset along line condition direction in velocity direction.
        Returns
        -------

        """

        current_model = KratosMultiphysics.Model()
        mp = current_model.CreateModelPart("solid_part")
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)

        # create nodes
        second_coord = [1, 0, 0.0]
        mp.CreateNewNode(1, 0.0, 0.0, 0.0)
        mp.CreateNewNode(2, second_coord[0], second_coord[1], 0.0)

        strategy = self.setup_strategy(mp)

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
                       "origin"          : [0,0,0],
                       "offset"          : 0.25
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
        strategy.InitializeSolutionStep()
        cond.CalculateLocalSystem(lhs, rhs, mp.ProcessInfo)

        self.checkRHS(rhs, [0.0, -0.5, 0.0, -1.5])

        # move load
        process.ExecuteFinalizeSolutionStep()
        process.ExecuteInitializeSolutionStep()

        # check if interpolation is done correctly
        cond.CalculateLocalSystem(lhs, rhs, mp.ProcessInfo)

        self.checkRHS(rhs, [0.0, -1.0, 0.0, -1.0])

    def _TestSetMovingLoadReverseGeomOffsetNegative(self):
        """
        Tests a moving load on 1 condition element, where the nodes of the element are reversed compared to the
        direction of the moving load, including a negative offset along line condition direction in velocity direction.
        Returns
        -------

        """

        current_model = KratosMultiphysics.Model()
        mp = current_model.CreateModelPart("solid_part")
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)

        # create nodes
        second_coord = [1, 0, 0.0]
        mp.CreateNewNode(1, 0.0, 0.0, 0.0)
        mp.CreateNewNode(2, second_coord[0], second_coord[1], 0.0)

        strategy = self.setup_strategy(mp)
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
                       "origin"          : [0,0,0],
                       "offset"          : -0.25
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
        strategy.InitializeSolutionStep()
        cond.CalculateLocalSystem(lhs, rhs, mp.ProcessInfo)

        self.checkRHS(rhs, [0.0, 0.0, 0.0, 0.0])

        # move load
        process.ExecuteFinalizeSolutionStep()
        process.ExecuteInitializeSolutionStep()

        # check if interpolation is done correctly
        strategy.InitializeSolutionStep()
        cond.CalculateLocalSystem(lhs, rhs, mp.ProcessInfo)

        self.checkRHS(rhs, [0.0, 0.0, 0.0, -2.0])

    def _TestSetMovingLoadMultipleConditions(self):
        """
        Tests a moving load on 2 condition elements, where the order of the elements is sorted in the direction of the
        moving load
        Returns
        -------

        """

        current_model = KratosMultiphysics.Model()
        mp = current_model.CreateModelPart("solid_part")
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)

        #create nodes
        second_coord = [1.0, 0.0, 0.0]
        third_coord = [2.0, 0.0, 0.0]
        mp.CreateNewNode(4, 0.0, 0.0, 0.0)
        mp.CreateNewNode(2, second_coord[0],second_coord[1],second_coord[2])
        mp.CreateNewNode(3, third_coord[0], third_coord[1], third_coord[2])

        strategy = self.setup_strategy(mp)
        # create condition
        conditions = []
        conditions.append(mp.CreateNewCondition("MovingLoadCondition2D2N", 1, [4, 2], mp.GetProperties()[1]))
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
        mp.ProcessInfo.SetValue(KratosMultiphysics.TIME, 0)
        mp.ProcessInfo.SetValue(KratosMultiphysics.DELTA_TIME, 0.5)
        process = SMA.SetMovingLoadProcess(mp, parameters)


        # initialize and set load
        process.ExecuteInitialize()
        process.ExecuteInitializeSolutionStep()

        # initialise matrices
        lhs = KratosMultiphysics.Matrix(0,0)
        rhs = KratosMultiphysics.Vector(0)

        # set load on node
        strategy.InitializeSolutionStep()

        all_rhs = []
        for cond in conditions:
            cond.CalculateLocalSystem(lhs, rhs, mp.ProcessInfo)
            all_rhs.append(list(rhs))

        self.checkRHS(all_rhs[0], [0.0, -2.0, 0.0, 0.0])
        self.checkRHS(all_rhs[1], [0.0, 0.0, 0.0, 0.0])

        # move load within first element
        process.ExecuteFinalizeSolutionStep()
        process.ExecuteInitializeSolutionStep()

        # check if interpolation is done correctly
        strategy.InitializeSolutionStep()

        all_rhs = []
        for cond in conditions:
            cond.CalculateLocalSystem(lhs, rhs, mp.ProcessInfo)
            all_rhs.append(list(rhs))

        self.checkRHS(all_rhs[0], [0.0, -1.0, 0.0, -1.0])
        self.checkRHS(all_rhs[1], [0.0, 0.0, 0.0, 0.0])

        # move load to element connection element
        process.ExecuteFinalizeSolutionStep()
        process.ExecuteInitializeSolutionStep()

        strategy.InitializeSolutionStep()
        all_rhs = []
        for cond in conditions:
            cond.CalculateLocalSystem(lhs, rhs, mp.ProcessInfo)
            all_rhs.append(list(rhs))

        self.checkRHS(all_rhs[0], [0.0, 0.0, 0.0, -2.0])
        self.checkRHS(all_rhs[1], [0.0, 0.0, 0.0, 0.0])

        # move load to next element
        process.ExecuteFinalizeSolutionStep()
        process.ExecuteInitializeSolutionStep()

        strategy.InitializeSolutionStep()
        all_rhs = []
        for cond in conditions:
            cond.CalculateLocalSystem(lhs, rhs, mp.ProcessInfo)
            all_rhs.append(list(rhs))

        self.checkRHS(all_rhs[0], [0.0, 0.0, 0.0, 0.0])
        self.checkRHS(all_rhs[1], [0.0, -1.0, 0.0, -1.0])

    def _TestSetMovingLoadMultipleConditionsOffSetPositive(self):
        """
        Tests a moving load on 2 condition elements, where the order of the elements is sorted in the direction of the
        moving load
        Returns
        -------

        """

        current_model = KratosMultiphysics.Model()
        mp = current_model.CreateModelPart("solid_part")
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)

        #create nodes
        second_coord = [1.0, 0.0, 0.0]
        third_coord = [2.0, 0.0, 0.0]
        mp.CreateNewNode(1, 0.0, 0.0, 0.0)
        mp.CreateNewNode(2, second_coord[0],second_coord[1],second_coord[2])
        mp.CreateNewNode(3, third_coord[0], third_coord[1], third_coord[2])

        strategy = self.setup_strategy(mp)

        # create condition
        conditions = []
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
                    "origin"          : [0,0,0],
                    "offset"          : 0.5
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
        strategy.InitializeSolutionStep()

        all_rhs = []
        for cond in conditions:
            cond.CalculateLocalSystem(lhs, rhs, mp.ProcessInfo)
            all_rhs.append(list(rhs))

        self.checkRHS(all_rhs[0], [0.0, -1.0, 0.0, -1.0])
        self.checkRHS(all_rhs[1], [0.0, 0.0, 0.0, 0.0])

        # move load within first element
        process.ExecuteFinalizeSolutionStep()
        process.ExecuteInitializeSolutionStep()

        # check if interpolation is done correctly
        strategy.InitializeSolutionStep()

        all_rhs = []
        for cond in conditions:
            cond.CalculateLocalSystem(lhs, rhs, mp.ProcessInfo)
            all_rhs.append(list(rhs))

        self.checkRHS(all_rhs[0], [0.0, 0.0, 0.0, -2.0])
        self.checkRHS(all_rhs[1], [0.0, 0.0, 0.0, 0.0])

        # move load to element connection element
        process.ExecuteFinalizeSolutionStep()
        process.ExecuteInitializeSolutionStep()

        strategy.InitializeSolutionStep()

        all_rhs = []
        for cond in conditions:
            cond.CalculateLocalSystem(lhs, rhs, mp.ProcessInfo)
            all_rhs.append(list(rhs))

        self.checkRHS(all_rhs[0], [0.0, 0.0, 0.0, 0.0])
        self.checkRHS(all_rhs[1], [0.0, -1.0, 0.0, -1.0])

        # move load to next element
        process.ExecuteFinalizeSolutionStep()
        process.ExecuteInitializeSolutionStep()

        strategy.InitializeSolutionStep()
        all_rhs = []
        for cond in conditions:
            cond.CalculateLocalSystem(lhs, rhs, mp.ProcessInfo)
            all_rhs.append(list(rhs))

        self.checkRHS(all_rhs[0], [0.0, 0.0, 0.0, 0.0])
        self.checkRHS(all_rhs[1], [0.0, 0.0, 0.0, -2.0])

    def _TestSetMovingLoadMultipleConditionsOffSetNegative(self):
        """
        Tests a moving load on 2 condition elements, where the order of the elements is sorted in the direction of the
        moving load
        Returns
        -------

        """

        current_model = KratosMultiphysics.Model()
        mp = current_model.CreateModelPart("solid_part")
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)

        #create nodes
        second_coord = [1.0, 0.0, 0.0]
        third_coord = [2.0, 0.0, 0.0]
        mp.CreateNewNode(1, 0.0, 0.0, 0.0)
        mp.CreateNewNode(2, second_coord[0],second_coord[1],second_coord[2])
        mp.CreateNewNode(3, third_coord[0], third_coord[1], third_coord[2])

        strategy = self.setup_strategy(mp)
        # create condition
        conditions = []
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
                    "origin"          : [0,0,0],
                    "offset"          : -0.5
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
        strategy.InitializeSolutionStep()

        all_rhs = []
        for cond in conditions:
            cond.CalculateLocalSystem(lhs, rhs, mp.ProcessInfo)
            all_rhs.append(list(rhs))

        self.checkRHS(all_rhs[0], [0.0, 0.0, 0.0, 0.0])
        self.checkRHS(all_rhs[1], [0.0, 0.0, 0.0, 0.0])

        # move load within first element
        process.ExecuteFinalizeSolutionStep()
        process.ExecuteInitializeSolutionStep()

        # check if interpolation is done correctly
        strategy.InitializeSolutionStep()

        all_rhs = []
        for cond in conditions:
            cond.CalculateLocalSystem(lhs, rhs, mp.ProcessInfo)
            all_rhs.append(list(rhs))

        self.checkRHS(all_rhs[0], [0.0, -2.0, 0.0, 0.0])
        self.checkRHS(all_rhs[1], [0.0, 0.0, 0.0, 0.0])

        # move load to element connection element
        process.ExecuteFinalizeSolutionStep()
        process.ExecuteInitializeSolutionStep()

        strategy.InitializeSolutionStep()
        all_rhs = []
        for cond in conditions:
            cond.CalculateLocalSystem(lhs, rhs, mp.ProcessInfo)
            all_rhs.append(list(rhs))

        self.checkRHS(all_rhs[0], [0.0, -1.0, 0.0, -1.0])
        self.checkRHS(all_rhs[1], [0.0, 0.0, 0.0, 0.0])

        # move load to next element
        process.ExecuteFinalizeSolutionStep()
        process.ExecuteInitializeSolutionStep()

        strategy.InitializeSolutionStep()
        all_rhs = []
        for cond in conditions:
            cond.CalculateLocalSystem(lhs, rhs, mp.ProcessInfo)
            all_rhs.append(list(rhs))

        self.checkRHS(all_rhs[0], [0.0, 0.0, 0.0, -2.0])
        self.checkRHS(all_rhs[1], [0.0, 0.0, 0.0, 0.0])

        # move load to next element
        process.ExecuteFinalizeSolutionStep()
        process.ExecuteInitializeSolutionStep()

        strategy.InitializeSolutionStep()
        all_rhs = []
        for cond in conditions:
            cond.CalculateLocalSystem(lhs, rhs, mp.ProcessInfo)
            all_rhs.append(list(rhs))

        self.checkRHS(all_rhs[0], [0.0, 0.0, 0.0, 0.0])
        self.checkRHS(all_rhs[1], [0.0, -1.0, 0.0, -1.0])

    def _TestSetMovingLoadMultipleConditionsReversed(self):
        """
        Tests a moving load on 2 condition elements, where the order of the elements is reversed compared to the moving
        direction of the load
        Returns
        -------

        """

        current_model = KratosMultiphysics.Model()
        mp = current_model.CreateModelPart("solid_part")
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)

        # create nodes
        second_coord = [1.0, 0.0, 0.0]
        third_coord = [2.0, 0.0, 0.0]
        mp.CreateNewNode(1, 0.0, 0.0, 0.0)
        mp.CreateNewNode(2, second_coord[0],second_coord[1],second_coord[2])
        mp.CreateNewNode(3, third_coord[0], third_coord[1], third_coord[2])

        strategy = self.setup_strategy(mp)

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
        strategy.InitializeSolutionStep()

        all_rhs = []
        for cond in conditions:
            cond.CalculateLocalSystem(lhs, rhs, mp.ProcessInfo)
            all_rhs.append(list(rhs))

        self.checkRHS(all_rhs[0], [0.0, 0.0, 0.0, 0.0])
        self.checkRHS(all_rhs[1], [0.0, 0.0, 0.0, -2.0])

        # move load within first element
        process.ExecuteFinalizeSolutionStep()
        process.ExecuteInitializeSolutionStep()

        # check if interpolation is done correctly
        strategy.InitializeSolutionStep()

        all_rhs = []
        for cond in conditions:
            cond.CalculateLocalSystem(lhs, rhs, mp.ProcessInfo)
            all_rhs.append(list(rhs))

        self.checkRHS(all_rhs[0], [0.0, 0.0, 0.0, 0.0])
        self.checkRHS(all_rhs[1], [0.0, -1.0, 0.0, -1.0])

        # move load to element connection element
        process.ExecuteFinalizeSolutionStep()
        process.ExecuteInitializeSolutionStep()


        # calculate load
        strategy.InitializeSolutionStep()
        all_rhs = []
        for cond in conditions:
            cond.CalculateLocalSystem(lhs, rhs, mp.ProcessInfo)
            all_rhs.append(list(rhs))

        self.checkRHS(all_rhs[0], [0.0, 0.0, 0.0, -2.0])
        self.checkRHS(all_rhs[1], [0.0, 0.0, 0.0, 0.0])

        # move load to next element, also increase time step
        mp.ProcessInfo.SetValue(KratosMultiphysics.DELTA_TIME, 0.75)
        process.ExecuteFinalizeSolutionStep()
        process.ExecuteInitializeSolutionStep()

        # calculate load
        strategy.InitializeSolutionStep()
        all_rhs = []
        for cond in conditions:
            cond.CalculateLocalSystem(lhs, rhs, mp.ProcessInfo)
            all_rhs.append(list(rhs))

        self.checkRHS(all_rhs[0], [0.0, -1.5, 0.0, -0.5])
        self.checkRHS(all_rhs[1], [0.0, 0.0, 0.0, 0.0])

    def _TestSetMovingLoadMultipleConditionsReversedOffsetPositive(self):
        """
        Tests a moving load on 2 condition elements, where the order of the elements is reversed compared to the moving
        direction of the load, including a positive offset along line condition direction in velocity direction.
        Returns
        -------

        """

        current_model = KratosMultiphysics.Model()
        mp = current_model.CreateModelPart("solid_part")
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)

        # create nodes
        second_coord = [1.0, 0.0, 0.0]
        third_coord = [2.0, 0.0, 0.0]
        mp.CreateNewNode(1, 0.0, 0.0, 0.0)
        mp.CreateNewNode(2, second_coord[0],second_coord[1],second_coord[2])
        mp.CreateNewNode(3, third_coord[0], third_coord[1], third_coord[2])

        strategy = self.setup_strategy(mp)

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
                    "origin"          : [0,0,0],
                    "offset"          : 0.5
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
        strategy.InitializeSolutionStep()

        all_rhs = []
        for cond in conditions:
            cond.CalculateLocalSystem(lhs, rhs, mp.ProcessInfo)
            all_rhs.append(list(rhs))

        self.checkRHS(all_rhs[0], [0.0, 0.0, 0.0, 0.0])
        self.checkRHS(all_rhs[1], [0.0, -1.0, 0.0, -1.0])

        # move load within first element
        process.ExecuteFinalizeSolutionStep()
        process.ExecuteInitializeSolutionStep()

        # check if interpolation is done correctly
        strategy.InitializeSolutionStep()
        all_rhs = []
        for cond in conditions:
            cond.CalculateLocalSystem(lhs, rhs, mp.ProcessInfo)
            all_rhs.append(list(rhs))

        self.checkRHS(all_rhs[0], [0.0, 0.0, 0.0, -2.0])
        self.checkRHS(all_rhs[1], [0.0, 0.0, 0.0, 0.0])

        # move load to element connection element
        process.ExecuteFinalizeSolutionStep()
        process.ExecuteInitializeSolutionStep()

        # calculate load
        strategy.InitializeSolutionStep()
        all_rhs = []
        for cond in conditions:
            cond.CalculateLocalSystem(lhs, rhs, mp.ProcessInfo)
            all_rhs.append(list(rhs))

        self.checkRHS(all_rhs[0], [0.0, -1.0, 0.0, -1.0])
        self.checkRHS(all_rhs[1], [0.0, 0.0, 0.0, 0.0])

        # move load to next element, also increase time step
        mp.ProcessInfo.SetValue(KratosMultiphysics.DELTA_TIME, 0.75)
        process.ExecuteFinalizeSolutionStep()
        process.ExecuteInitializeSolutionStep()

        # calculate load
        strategy.InitializeSolutionStep()
        all_rhs = []
        for cond in conditions:
            cond.CalculateLocalSystem(lhs, rhs, mp.ProcessInfo)
            all_rhs.append(list(rhs))

        self.checkRHS(all_rhs[0], [0.0, 0.0, 0.0, 0.0])
        self.checkRHS(all_rhs[1], [0.0, 0.0, 0.0, 0.0])

    def _TestSetMovingLoadMultipleConditionsReversedOffsetNegative(self):
        """
        Tests a moving load on 2 condition elements, where the order of the elements is reversed compared to the moving
        direction of the load, including a negative offset along line condition direction in velocity direction.
        Returns
        -------

        """

        current_model = KratosMultiphysics.Model()
        mp = current_model.CreateModelPart("solid_part")
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)

        # create nodes
        second_coord = [1.0, 0.0, 0.0]
        third_coord = [2.0, 0.0, 0.0]
        mp.CreateNewNode(1, 0.0, 0.0, 0.0)
        mp.CreateNewNode(2, second_coord[0],second_coord[1],second_coord[2])
        mp.CreateNewNode(3, third_coord[0], third_coord[1], third_coord[2])

        strategy = self.setup_strategy(mp)

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
                    "origin"          : [0,0,0],
                    "offset"          : -0.5
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
        strategy.InitializeSolutionStep()

        all_rhs = []
        for cond in conditions:
            cond.CalculateLocalSystem(lhs, rhs, mp.ProcessInfo)
            all_rhs.append(list(rhs))

        self.checkRHS(all_rhs[0], [0.0, 0.0, 0.0, 0.0])
        self.checkRHS(all_rhs[1], [0.0, 0.0, 0.0, 0.0])

        # move load within first element
        process.ExecuteFinalizeSolutionStep()
        process.ExecuteInitializeSolutionStep()

        # check if interpolation is done correctly
        strategy.InitializeSolutionStep()
        all_rhs = []
        for cond in conditions:
            cond.CalculateLocalSystem(lhs, rhs, mp.ProcessInfo)
            all_rhs.append(list(rhs))

        self.checkRHS(all_rhs[0], [0.0, 0.0, 0.0, 0.0])
        self.checkRHS(all_rhs[1], [0.0, 0.0, 0.0, -2.0])

        # move load to element connection element
        process.ExecuteFinalizeSolutionStep()
        process.ExecuteInitializeSolutionStep()

        # calculate load
        strategy.InitializeSolutionStep()
        all_rhs = []
        for cond in conditions:
            cond.CalculateLocalSystem(lhs, rhs, mp.ProcessInfo)
            all_rhs.append(list(rhs))

        self.checkRHS(all_rhs[0], [0.0, 0.0, 0.0, 0.0])
        self.checkRHS(all_rhs[1], [0.0, -1.0, 0.0, -1.0])

        # move load to next element, also increase time step
        mp.ProcessInfo.SetValue(KratosMultiphysics.DELTA_TIME, 0.75)
        process.ExecuteFinalizeSolutionStep()
        process.ExecuteInitializeSolutionStep()

        # calculate load
        strategy.InitializeSolutionStep()
        all_rhs = []
        for cond in conditions:
            cond.CalculateLocalSystem(lhs, rhs, mp.ProcessInfo)
            all_rhs.append(list(rhs))

        self.checkRHS(all_rhs[0], [0.0, -0.5, 0.0, -1.5])
        self.checkRHS(all_rhs[1], [0.0, 0.0, 0.0, 0.0])

    def _TestSetMovingLoadMultipleConditionsDifferentOrigin(self):
        """
        Tests a moving load on 2 condition elements, where the order of the elements is sorted in the direction of the
        moving load
        Returns
        -------

        """

        current_model = KratosMultiphysics.Model()
        mp = current_model.CreateModelPart("solid_part")
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)

        #create nodes
        second_coord = [1.0, 0.0, 0.0]
        third_coord = [2.0, 0.0, 0.0]
        mp.CreateNewNode(1, 0.0, 0.0, 0.0)
        mp.CreateNewNode(2, second_coord[0],second_coord[1],second_coord[2])
        mp.CreateNewNode(3, third_coord[0], third_coord[1], third_coord[2])

        strategy = self.setup_strategy(mp)

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
        strategy.InitializeSolutionStep()

        all_rhs = []
        for cond in conditions:
            cond.CalculateLocalSystem(lhs, rhs, mp.ProcessInfo)
            all_rhs.append(list(rhs))

        self.checkRHS(all_rhs[0], [0.0, 0.0, 0.0, 0.0])
        self.checkRHS(all_rhs[1], [0.0, -1.5, 0.0, -0.5])


    def _TestSetMovingLoadMultipleConditionsDifferentOriginOffsetPositive(self):
        """
        Tests a moving load on 2 condition elements, where the order of the elements is sorted in the direction of the
        moving load, including a positive offset along line condition direction in velocity direction.
        Returns
        -------

        """

        current_model = KratosMultiphysics.Model()
        mp = current_model.CreateModelPart("solid_part")
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)

        #create nodes
        second_coord = [1.0, 0.0, 0.0]
        third_coord = [2.0, 0.0, 0.0]
        mp.CreateNewNode(1, 0.0, 0.0, 0.0)
        mp.CreateNewNode(2, second_coord[0],second_coord[1],second_coord[2])
        mp.CreateNewNode(3, third_coord[0], third_coord[1], third_coord[2])

        strategy = self.setup_strategy(mp)

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
                    "origin"          : [1.25,0,0],
                    "offset"          : 0.75
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
        strategy.InitializeSolutionStep()

        all_rhs = []
        for cond in conditions:
            cond.CalculateLocalSystem(lhs, rhs, mp.ProcessInfo)
            all_rhs.append(list(rhs))

        self.checkRHS(all_rhs[0], [0.0, 0.0, 0.0, 0.0])
        self.checkRHS(all_rhs[1], [0.0, 0.0, 0.0, -2.0])


    def _TestSetMovingLoadMultipleConditionsDifferentOriginOffsetNegative(self):
        """
        Tests a moving load on 2 condition elements, where the order of the elements is sorted in the direction of the
        moving load, including a negative offset along line condition direction in velocity direction.
        Returns
        -------

        """

        current_model = KratosMultiphysics.Model()
        mp = current_model.CreateModelPart("solid_part")
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)

        #create nodes
        second_coord = [1.0, 0.0, 0.0]
        third_coord = [2.0, 0.0, 0.0]
        mp.CreateNewNode(1, 0.0, 0.0, 0.0)
        mp.CreateNewNode(2, second_coord[0],second_coord[1],second_coord[2])
        mp.CreateNewNode(3, third_coord[0], third_coord[1], third_coord[2])

        strategy = self.setup_strategy(mp)

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
                    "origin"          : [1.25,0.0,0.0],
                    "offset"          : -0.25
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
        strategy.InitializeSolutionStep()

        all_rhs = []
        for cond in conditions:
            cond.CalculateLocalSystem(lhs, rhs, mp.ProcessInfo)
            all_rhs.append(list(rhs))

        self.checkRHS(all_rhs[0], [0.0, 0.0, 0.0, -2.0])
        self.checkRHS(all_rhs[1], [0.0, 0.0, 0.0, 0.0])


    def _TestSetMovingLoadMultipleConditionsDifferentOriginReversed(self):
        """
        Tests a moving load on 2 condition elements, where the order of the elements is sorted in the direction of the
        moving load
        Returns
        -------

        """

        current_model = KratosMultiphysics.Model()
        mp = current_model.CreateModelPart("solid_part")
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)

        # create nodes
        second_coord = [1.0, 0.0, 0.0]
        third_coord = [2.0, 0.0, 0.0]
        mp.CreateNewNode(1, 0.0, 0.0, 0.0)
        mp.CreateNewNode(2, second_coord[0],second_coord[1],second_coord[2])
        mp.CreateNewNode(3, third_coord[0], third_coord[1], third_coord[2])

        strategy = self.setup_strategy(mp)

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
        strategy.InitializeSolutionStep()

        all_rhs = []
        for cond in conditions:
            cond.CalculateLocalSystem(lhs, rhs, mp.ProcessInfo)
            all_rhs.append(list(rhs))

        self.checkRHS(all_rhs[0], [0.0, -0.5, 0.0, -1.5])
        self.checkRHS(all_rhs[1], [0.0, 0.0, 0.0, 0.0])


    def _TestSetMovingLoadMultipleConditionsDifferentOriginReversedOffsetPositive(self):
        """
        Tests a moving load on 2 condition elements, where the order of the elements is sorted in the direction of the
        moving load, including a positive offset along line condition direction in velocity direction.
        Returns
        -------

        """

        current_model = KratosMultiphysics.Model()
        mp = current_model.CreateModelPart("solid_part")
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)

        # create nodes
        second_coord = [1.0, 0.0, 0.0]
        third_coord = [2.0, 0.0, 0.0]
        mp.CreateNewNode(1, 0.0, 0.0, 0.0)
        mp.CreateNewNode(2, second_coord[0],second_coord[1],second_coord[2])
        mp.CreateNewNode(3, third_coord[0], third_coord[1], third_coord[2])

        strategy = self.setup_strategy(mp)

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
                    "origin"          : [1.25,0,0],
                    "offset"          : 0.25
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
        strategy.InitializeSolutionStep()

        all_rhs = []
        for cond in conditions:
            cond.CalculateLocalSystem(lhs, rhs, mp.ProcessInfo)
            all_rhs.append(list(rhs))

        self.checkRHS(all_rhs[0], [0.0, -1.0, 0.0, -1.0])
        self.checkRHS(all_rhs[1], [0.0, 0.0, 0.0, 0.0])

    def _TestSetMovingLoadMultipleConditionsDifferentOriginReversedOffsetNegative(self):
        """
        Tests a moving load on 2 condition elements, where the order of the elements is sorted in the direction of the
        moving load, including a negative offset along line condition direction in velocity direction.
        Returns
        -------

        """

        current_model = KratosMultiphysics.Model()
        mp = current_model.CreateModelPart("solid_part")
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)

        # create nodes
        second_coord = [1.0, 0.0, 0.0]
        third_coord = [2.0, 0.0, 0.0]
        mp.CreateNewNode(1, 0.0, 0.0, 0.0)
        mp.CreateNewNode(2, second_coord[0],second_coord[1],second_coord[2])
        mp.CreateNewNode(3, third_coord[0], third_coord[1], third_coord[2])

        strategy = self.setup_strategy(mp)

        # create condition
        conditions = []
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
                    "origin"          : [1.25,0,0],
                    "offset"          : -1.25
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
        strategy.InitializeSolutionStep()

        all_rhs = []
        for cond in conditions:
            cond.CalculateLocalSystem(lhs, rhs, mp.ProcessInfo)
            all_rhs.append(list(rhs))

        self.checkRHS(all_rhs[0], [0.0, 0.0, 0.0, 0.0])
        self.checkRHS(all_rhs[1], [0.0, 0.0, 0.0, -2.0])

    def _TestSetMovingLoadWithLoadFunction(self):
        """
       Tests a moving load on a condition element, where the load is a function of time
       Returns
       -------

       """

        current_model = KratosMultiphysics.Model()
        mp = current_model.CreateModelPart("solid_part")
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)

        # create nodes
        second_coord = [1, 0, 0.0]
        mp.CreateNewNode(1,0.0,0.0,0.0)
        mp.CreateNewNode(2,second_coord[0],second_coord[1],0.0)

        strategy = self.setup_strategy(mp)

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
        strategy.InitializeSolutionStep()
        cond.CalculateLocalSystem(lhs, rhs, mp.ProcessInfo)

        self.checkRHS(rhs, [0.0, 0.0, 0.0, 0.0])

        # change time and recalculate load
        mp.ProcessInfo.SetValue(KratosMultiphysics.TIME, 0.5)
        process.ExecuteInitializeSolutionStep()

        # check if interpolation is done correctly
        strategy.InitializeSolutionStep()
        cond.CalculateLocalSystem(lhs, rhs, mp.ProcessInfo)

        self.checkRHS(rhs, [0.0, -0.5, 0.0, -0.5])

    def _TestSetMovingLoadWithLoadFunctionOffsetPositive(self):
        """
       Tests a moving load on a condition element, where the load is a function of time, including a positive offset along line condition direction in velocity direction.
       Returns
       -------

       """

        current_model = KratosMultiphysics.Model()
        mp = current_model.CreateModelPart("solid_part")
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)

        # create nodes
        second_coord = [1, 0, 0.0]
        mp.CreateNewNode(1,0.0,0.0,0.0)
        mp.CreateNewNode(2,second_coord[0],second_coord[1],0.0)

        strategy = self.setup_strategy(mp)

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
                    "origin"          : [0.5,0,0],
                    "offset"          : 0.25
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
        strategy.InitializeSolutionStep()
        cond.CalculateLocalSystem(lhs, rhs, mp.ProcessInfo)

        self.checkRHS(rhs, [0.0, 0.0, 0.0, 0.0])

        # change time and recalculate load
        mp.ProcessInfo.SetValue(KratosMultiphysics.TIME, 0.5)
        process.ExecuteInitializeSolutionStep()

        # check if interpolation is done correctly
        strategy.InitializeSolutionStep()
        cond.CalculateLocalSystem(lhs, rhs, mp.ProcessInfo)

        self.checkRHS(rhs, [0.0, -0.25, 0.0, -0.75])

    def _TestSetMovingLoadWithLoadFunctionOffsetNegative(self):
        """
       Tests a moving load on a condition element, where the load is a function of time, including a negative offset along line condition direction in velocity direction.
       Returns
       -------

       """

        current_model = KratosMultiphysics.Model()
        mp = current_model.CreateModelPart("solid_part")
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)

        # create nodes
        second_coord = [1, 0, 0.0]
        mp.CreateNewNode(1,0.0,0.0,0.0)
        mp.CreateNewNode(2,second_coord[0],second_coord[1],0.0)

        strategy = self.setup_strategy(mp)

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
                    "origin"          : [0.5,0,0],
                    "offset"          : -0.25
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
        strategy.InitializeSolutionStep()
        cond.CalculateLocalSystem(lhs, rhs, mp.ProcessInfo)

        self.checkRHS(rhs, [0.0, 0.0, 0.0, 0.0])

        # change time and recalculate load
        mp.ProcessInfo.SetValue(KratosMultiphysics.TIME, 0.5)
        process.ExecuteInitializeSolutionStep()

        # check if interpolation is done correctly
        strategy.InitializeSolutionStep()
        cond.CalculateLocalSystem(lhs, rhs, mp.ProcessInfo)

        self.checkRHS(rhs, [0.0, -0.75, 0.0, -0.25])

    def _TestSetMovingLoadWithVelocityFunction(self):
        """
       Tests a moving load on a condition element, where the load velocity is a function of time.

       Returns
       -------

       """

        current_model = KratosMultiphysics.Model()
        mp = current_model.CreateModelPart("solid_part")
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)

        # create nodes
        second_coord = [1, 0, 0.0]
        mp.CreateNewNode(1,0.0,0.0,0.0)
        mp.CreateNewNode(2,second_coord[0],second_coord[1],0.0)

        strategy = self.setup_strategy(mp)

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
        strategy.InitializeSolutionStep()
        cond.CalculateLocalSystem(lhs, rhs, mp.ProcessInfo)

        self.checkRHS(rhs, [0.0, -2.0, 0.0, 0.0])

        # change time and recalculate load
        process.ExecuteFinalizeSolutionStep()
        process.ExecuteInitializeSolutionStep()

        # check if interpolation is done correctly
        strategy.InitializeSolutionStep()
        cond.CalculateLocalSystem(lhs, rhs, mp.ProcessInfo)

        self.checkRHS(rhs, [0.0, -2.0, 0.0, 0.0])

        mp.ProcessInfo.SetValue(KratosMultiphysics.TIME, 0.5)

        process.ExecuteFinalizeSolutionStep()
        process.ExecuteInitializeSolutionStep()

        # check if interpolation is done correctly
        strategy.InitializeSolutionStep()
        cond.CalculateLocalSystem(lhs, rhs, mp.ProcessInfo)

        self.checkRHS(rhs, [0.0, -1.5, 0.0, -0.5])

    def _TestSetMovingLoadWithVelocityFunctionOffsetPositive(self):
        """
       Tests a moving load on a condition element, where the load velocity is a function of time, including a positive offset along line condition direction in velocity direction.

       Returns
       -------

       """

        current_model = KratosMultiphysics.Model()
        mp = current_model.CreateModelPart("solid_part")
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)

        # create nodes
        second_coord = [1, 0, 0.0]
        mp.CreateNewNode(1,0.0,0.0,0.0)
        mp.CreateNewNode(2,second_coord[0],second_coord[1],0.0)

        strategy = self.setup_strategy(mp)

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
                    "origin"          : [0.0,0,0],
                    "offset"          : 0.5
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
        strategy.InitializeSolutionStep()
        cond.CalculateLocalSystem(lhs, rhs, mp.ProcessInfo)

        self.checkRHS(rhs, [0.0, -1.0, 0.0, -1.0])

        # change time and recalculate load
        process.ExecuteFinalizeSolutionStep()
        process.ExecuteInitializeSolutionStep()

        # check if interpolation is done correctly
        strategy.InitializeSolutionStep()
        cond.CalculateLocalSystem(lhs, rhs, mp.ProcessInfo)

        self.checkRHS(rhs, [0.0, -1.0, 0.0, -1.0])

        mp.ProcessInfo.SetValue(KratosMultiphysics.TIME, 0.5)

        process.ExecuteFinalizeSolutionStep()
        process.ExecuteInitializeSolutionStep()

        # check if interpolation is done correctly
        strategy.InitializeSolutionStep()
        cond.CalculateLocalSystem(lhs, rhs, mp.ProcessInfo)

        self.checkRHS(rhs, [0.0, -0.5, 0.0, -1.5])

    def _TestSetMovingLoadWithVelocityFunctionOffsetNegative(self):
        """
       Tests a moving load on a condition element, where the load velocity is a function of time, including a negative offset along line condition direction in velocity direction.

       Returns
       -------

       """

        current_model = KratosMultiphysics.Model()
        mp = current_model.CreateModelPart("solid_part")
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)

        # create nodes
        second_coord = [1, 0, 0.0]
        mp.CreateNewNode(1,0.0,0.0,0.0)
        mp.CreateNewNode(2,second_coord[0],second_coord[1],0.0)

        strategy = self.setup_strategy(mp)

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
                    "origin"          : [0.0,0,0],
                    "offset"          : -0.25
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
        strategy.InitializeSolutionStep()
        cond.CalculateLocalSystem(lhs, rhs, mp.ProcessInfo)

        self.checkRHS(rhs, [0.0, 0.0, 0.0, 0.0])

        # change time and recalculate load
        process.ExecuteFinalizeSolutionStep()
        process.ExecuteInitializeSolutionStep()

        # check if interpolation is done correctly
        strategy.InitializeSolutionStep()
        cond.CalculateLocalSystem(lhs, rhs, mp.ProcessInfo)

        self.checkRHS(rhs, [0.0, 0.0, 0.0, 0.0])

        mp.ProcessInfo.SetValue(KratosMultiphysics.TIME, 0.5)

        process.ExecuteFinalizeSolutionStep()
        process.ExecuteInitializeSolutionStep()

        # check if interpolation is done correctly
        strategy.InitializeSolutionStep()
        cond.CalculateLocalSystem(lhs, rhs, mp.ProcessInfo)

        self.checkRHS(rhs, [0.0, -2.0, 0.0, 0.0])


    def test_SetMovingLoad(self):
        self._TestSetMovingLoad()

    def test_SetMovingLoadOffsetPositive(self):
        self._TestSetMovingLoadOffsetPositive()

    def test_SetMovingLoadOffsetNegative(self):
        self._TestSetMovingLoadOffsetNegative()

    def test_SetMovingLoadReverseGeom(self):
        self._TestSetMovingLoadReverseGeom()

    def test_SetMovingLoadReverseGeomOffsetPositive(self):
        self._TestSetMovingLoadReverseGeomOffsetPositive()

    def test_SetMovingLoadReverseGeomOffsetNegative(self):
        self._TestSetMovingLoadReverseGeomOffsetNegative()

    def test_SetMovingLoadMultipleConditions(self):
        self._TestSetMovingLoadMultipleConditions()

    def test_SetMovingLoadMultipleConditionsOffsetPositive(self):
        self._TestSetMovingLoadMultipleConditionsOffSetPositive()

    def test_SetMovingLoadMultipleConditionsOffsetNegative(self):
        self._TestSetMovingLoadMultipleConditionsOffSetNegative()

    def test_SetMovingLoadMultipleConditionsReversed(self):
        self._TestSetMovingLoadMultipleConditionsReversed()

    def test_SetMovingLoadMultipleConditionsReversedOffsetPositive(self):
        self._TestSetMovingLoadMultipleConditionsReversedOffsetPositive()

    def test_SetMovingLoadMultipleConditionsReversedOffsetNegative(self):
        self._TestSetMovingLoadMultipleConditionsReversedOffsetNegative()

    def test_SetMovingLoadMultipleConditionsDifferentOrigin(self):
        self._TestSetMovingLoadMultipleConditionsDifferentOrigin()

    def test_SetMovingLoadMultipleConditionsDifferentOriginOffsetPositive(self):
        self._TestSetMovingLoadMultipleConditionsDifferentOriginOffsetPositive()

    def test_SetMovingLoadMultipleConditionsDifferentOriginOffsetNegative(self):
        self._TestSetMovingLoadMultipleConditionsDifferentOriginOffsetNegative()

    def test_SetMovingLoadMultipleConditionsDifferentOriginReversed(self):
        self._TestSetMovingLoadMultipleConditionsDifferentOriginReversed()

    def test_SetMovingLoadMultipleConditionsDifferentOriginReversedOffsetPositive(self):
        self._TestSetMovingLoadMultipleConditionsDifferentOriginReversedOffsetPositive()

    def test_SetMovingLoadMultipleConditionsDifferentOriginReversedOffsetNegative(self):
        self._TestSetMovingLoadMultipleConditionsDifferentOriginReversedOffsetNegative()

    def test_SetMovingLoadWithLoadFunction(self):
        self._TestSetMovingLoadWithLoadFunction()

    def test_SetMovingLoadWithLoadFunctionOffsetPositive(self):
        self._TestSetMovingLoadWithLoadFunctionOffsetPositive()

    def test_SetMovingLoadWithLoadFunctionOffsetNegative(self):
        self._TestSetMovingLoadWithLoadFunctionOffsetNegative()

    def test_SetMovingLoadWithVelocityFunction(self):
        self._TestSetMovingLoadWithVelocityFunction()

    def test_SetMovingLoadWithVelocityFunctionOffsetPositive(self):
        self._TestSetMovingLoadWithVelocityFunctionOffsetPositive()

    def test_SetMovingLoadWithVelocityFunctionOffsetNegative(self):
        self._TestSetMovingLoadWithVelocityFunctionOffsetNegative()


if __name__ == '__main__':
    KratosUnittest.main()

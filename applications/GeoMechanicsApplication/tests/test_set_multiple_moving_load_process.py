import KratosMultiphysics
import KratosMultiphysics.GeoMechanicsApplication as GMA
import KratosMultiphysics.KratosUnittest as KratosUnittest


class KratosGeoMechanicsSetMultipleMovingLoadProcessTests(KratosUnittest.TestCase):
    time_dependent_load = '["0.0", "-2.0*t", "0.0"]'
    time_dependent_velocity = ' "2.0*t" '

    def set_nodal_displacement_dof(self):
        for node in self.model_part.Nodes:
            node.AddDof(KratosMultiphysics.DISPLACEMENT_X)
            node.AddDof(KratosMultiphysics.DISPLACEMENT_Y)
            node.AddDof(KratosMultiphysics.DISPLACEMENT_Z)

    @staticmethod
    def setup_strategy(model_part):
        """
        Setup default strategy for solving the problem

        Parameters
        ----------
        model_part: model part

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
        strategy = KratosMultiphysics.ResidualBasedNewtonRaphsonStrategy(model_part,
                                                                         scheme,
                                                                         convergence_criterion,
                                                                         builder_and_solver,
                                                                         max_iters,
                                                                         compute_reactions,
                                                                         reform_step_dofs,
                                                                         move_mesh_flag)

        strategy.SetEchoLevel(0)
        return strategy

    def setUp(self):
        self.current_model = KratosMultiphysics.Model()
        self.root_model_part = self.current_model.CreateModelPart("root_part")
        self.model_part = self.root_model_part.CreateSubModelPart("solid_part")
        self.compute_model_part = self.root_model_part.CreateSubModelPart("compute_part")
        self.base_parameters = KratosMultiphysics.Parameters("""
                {
                    "help"                      : "This process applies multiple moving load conditions belonging to a modelpart. The load moves over line elements.",
                    "model_part_name"           : "solid_part",
                    "compute_model_part_name"   : "compute_part",
                    "variable_name"             : "POINT_LOAD",
                    "load"                      : [0.0, -2.0, 0.0],
                    "direction"                 : [1,1,1],
                    "velocity"                  : 1,
                    "origin"                    : [0,0,0]
                }
                """)
        # Add nodal solution space
        self.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        self.strategy = self.setup_strategy(self.compute_model_part)

    def create_two_nodes_condition(self, reverse_nodes=False):
        # create nodes
        self.model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
        self.model_part.CreateNewNode(2, 1.0, 0.0, 0.0)
        self.set_nodal_displacement_dof()

        # create condition
        if reverse_nodes:
            self.model_part.CreateNewCondition("MovingLoadCondition2D2N", 1, [2, 1], self.model_part.GetProperties()[1])
        else:
            self.model_part.CreateNewCondition("MovingLoadCondition2D2N", 1, [1, 2], self.model_part.GetProperties()[1])

    def create_three_nodes_condition(self, reverse_nodes=False):
        # create nodes
        self.model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
        self.model_part.CreateNewNode(2, 1.0, 0.0, 0.0)
        self.model_part.CreateNewNode(3, 2.0, 0.0, 0.0)
        self.set_nodal_displacement_dof()

        # create condition
        if reverse_nodes:
            self.model_part.CreateNewCondition("MovingLoadCondition2D2N", 1, [3, 2], self.model_part.GetProperties()[1])
            self.model_part.CreateNewCondition("MovingLoadCondition2D2N", 2, [2, 1], self.model_part.GetProperties()[1])
        else:
            self.model_part.CreateNewCondition("MovingLoadCondition2D2N", 1, [1, 2], self.model_part.GetProperties()[1])
            self.model_part.CreateNewCondition("MovingLoadCondition2D2N", 2, [2, 3], self.model_part.GetProperties()[1])

    def initialize_test(self, process):
        process.ExecuteInitialize()
        process.ExecuteInitializeSolutionStep()
        self.strategy.InitializeSolutionStep()

    def next_solution_step(self, process):
        process.ExecuteFinalizeSolutionStep()
        process.ExecuteInitializeSolutionStep()
        self.strategy.InitializeSolutionStep()

    def tearDown(self):
        self.root_model_part.Clear()

    def checkRHS(self, rhs, expected_res, tols=None):
        """
        routine to check calculation of rhs side within context of testing SetMovingLoad
        Returns
        -------
        """

        if tols is None:
            tols = [None] * len(rhs)

        for rhs_val, expected_val, tol in zip(rhs, expected_res, tols):
            self.assertAlmostEqual(rhs_val, expected_val, tol)

    def test_SetMultipleMovingLoads(self):
        """
        Tests a moving load on 1 condition element, where the nodes of the element are sorted in the direction of the
        moving load - a single load
        Returns
        -------

        """

        self.create_two_nodes_condition()

        parameters = self.base_parameters
        self.model_part.ProcessInfo.SetValue(KratosMultiphysics.TIME, 0)
        self.model_part.ProcessInfo.SetValue(KratosMultiphysics.DELTA_TIME, 0.25)

        process = GMA.SetMultipleMovingLoadsProcess(self.model_part, parameters)
        cond = self.compute_model_part.GetCondition(2)

        # initialise and set load
        self.initialize_test(process)

        # initialise matrices
        lhs = KratosMultiphysics.Matrix(0, 0)
        rhs = KratosMultiphysics.Vector(0)

        # set load on node
        cond.CalculateLocalSystem(lhs, rhs, self.model_part.ProcessInfo)

        self.checkRHS(rhs, [0.0, -2.0, 0.0, 0.0])

        # move load
        self.next_solution_step(process)

        # check if interpolation is done correctly
        cond.CalculateLocalSystem(lhs, rhs, self.model_part.ProcessInfo)

        self.checkRHS(rhs, [0.0, -1.5, 0.0, -0.5])

    def test_SetMultipleMovingLoadsConfigurationPositive(self):
        """
        Tests a moving load on 1 condition element, where the nodes of the element are sorted in the direction of the
        moving load, including a positive configuration along line condition direction in velocity direction.
        Returns
        -------

        """

        self.create_two_nodes_condition()

        parameters = self.base_parameters
        parameters.AddVector("configuration", [0.25])

        self.model_part.ProcessInfo.SetValue(KratosMultiphysics.TIME, 0)
        self.model_part.ProcessInfo.SetValue(KratosMultiphysics.DELTA_TIME, 0.25)

        process = GMA.SetMultipleMovingLoadsProcess(self.model_part, parameters)
        cond = self.compute_model_part.GetCondition(2)

        # initialise and set load
        self.initialize_test(process)

        # initialise matrices
        lhs = KratosMultiphysics.Matrix(0, 0)
        rhs = KratosMultiphysics.Vector(0)

        # set load on node
        cond.CalculateLocalSystem(lhs, rhs, self.model_part.ProcessInfo)

        self.checkRHS(rhs, [0.0, -1.5, 0.0, -0.5])

        # move load
        self.next_solution_step(process)

        # check if interpolation is done correctly
        cond.CalculateLocalSystem(lhs, rhs, self.model_part.ProcessInfo)

        self.checkRHS(rhs, [0.0, -1.0, 0.0, -1.0])

    def test_SetMultipleMovingLoadsConfigurationNegative(self):
        """
        Tests a moving load on 1 condition element, where the nodes of the element are sorted in the direction of the
        moving load, , including a negative configuration along line condition direction in velocity direction.
        Returns
        -------

        """

        self.create_two_nodes_condition()

        parameters = self.base_parameters
        parameters.AddVector("configuration", [-0.25])

        self.model_part.ProcessInfo.SetValue(KratosMultiphysics.TIME, 0)
        self.model_part.ProcessInfo.SetValue(KratosMultiphysics.DELTA_TIME, 0.25)

        process = GMA.SetMultipleMovingLoadsProcess(self.model_part, parameters)
        cond = self.compute_model_part.GetCondition(2)

        # initialise and set load
        self.initialize_test(process)

        # initialise matrices
        lhs = KratosMultiphysics.Matrix(0, 0)
        rhs = KratosMultiphysics.Vector(0)

        # set load on node
        cond.CalculateLocalSystem(lhs, rhs, self.model_part.ProcessInfo)

        self.checkRHS(rhs, [0.0, 0.0, 0.0, 0.0])

        # move load
        self.next_solution_step(process)

        # check if interpolation is done correctly
        cond.CalculateLocalSystem(lhs, rhs, self.model_part.ProcessInfo)

        self.checkRHS(rhs, [0.0, -2.0, 0.0, 0.0])

    def test_SetMultipleMovingLoadsConfigurationCombined(self):
        """
        Tests a moving load on 1 condition element, where the nodes of the element are sorted in the direction of the
        moving load, , including a combined configuration along line condition direction in velocity direction.
        Returns
        -------

        """

        self.create_two_nodes_condition()

        parameters = self.base_parameters
        parameters.AddVector("configuration", [-0.25, 0, 0.25])

        self.model_part.ProcessInfo.SetValue(KratosMultiphysics.TIME, 0)
        self.model_part.ProcessInfo.SetValue(KratosMultiphysics.DELTA_TIME, 0.25)

        process = GMA.SetMultipleMovingLoadsProcess(self.model_part, parameters)
        conditions = []
        conditions.append(self.compute_model_part.GetCondition(2))
        conditions.append(self.compute_model_part.GetCondition(3))
        conditions.append(self.compute_model_part.GetCondition(4))

        # initialise and set load
        self.initialize_test(process)

        # set load on node
        all_rhs = []
        for cond in conditions:
            # initialise matrices
            lhs = KratosMultiphysics.Matrix(0, 0)
            rhs = KratosMultiphysics.Vector(0)
            cond.CalculateLocalSystem(lhs, rhs, self.model_part.ProcessInfo)
            all_rhs.append(list(rhs))

        self.checkRHS(all_rhs[0], [0.0, 0.0, 0.0, 0.0])
        self.checkRHS(all_rhs[1], [0.0, -2.0, 0.0, 0.0])
        self.checkRHS(all_rhs[2], [0.0, -1.5, 0.0, -0.5])

        # move load within first element
        self.next_solution_step(process)

        # check if interpolation is done correctly
        all_rhs = []
        for cond in conditions:
            lhs = KratosMultiphysics.Matrix(0, 0)
            rhs = KratosMultiphysics.Vector(0)
            cond.CalculateLocalSystem(lhs, rhs, self.model_part.ProcessInfo)
            all_rhs.append(list(rhs))

        self.checkRHS(all_rhs[0], [0.0, -2.0, 0.0, 0.0])
        self.checkRHS(all_rhs[1], [0.0, -1.5, 0.0, -0.5])
        self.checkRHS(all_rhs[2], [0.0, -1.0, 0.0, -1.0])

    def test_SetMultipleMovingLoadsReverseGeom(self):
        """
        Tests a moving load on 1 condition element, where the nodes of the element are reversed compared to the
        direction of the moving load
        Returns
        -------

        """

        self.create_two_nodes_condition(reverse_nodes=True)

        parameters = self.base_parameters

        self.model_part.ProcessInfo.SetValue(KratosMultiphysics.TIME, 0)
        self.model_part.ProcessInfo.SetValue(KratosMultiphysics.DELTA_TIME, 0.25)

        process = GMA.SetMultipleMovingLoadsProcess(self.model_part, parameters)
        cond = self.compute_model_part.GetCondition(2)

        # initialize and set load
        self.initialize_test(process)

        # initialise matrices
        lhs = KratosMultiphysics.Matrix(0, 0)
        rhs = KratosMultiphysics.Vector(0)

        # set load on node
        # cond.SetValue(SMA.MOVING_LOAD_LOCAL_DISTANCE, 0)
        cond.CalculateLocalSystem(lhs, rhs, self.model_part.ProcessInfo)

        self.checkRHS(rhs, [0.0, 0.0, 0.0, -2.0])

        # move load
        self.next_solution_step(process)

        # check if interpolation is done correctly
        cond.CalculateLocalSystem(lhs, rhs, self.model_part.ProcessInfo)

        self.checkRHS(rhs, [0.0, -0.5, 0.0, -1.5])

    def test_SetMultipleMovingLoadsReverseGeomConfigurationPositive(self):
        """
        Tests a moving load on 1 condition element, where the nodes of the element are reversed compared to the
        direction of the moving load, including a positive configuration along line condition direction in velocity direction.
        Returns
        -------

        """

        self.create_two_nodes_condition(reverse_nodes=True)

        parameters = self.base_parameters
        parameters.AddVector("configuration", [0.25])

        self.model_part.ProcessInfo.SetValue(KratosMultiphysics.TIME, 0)
        self.model_part.ProcessInfo.SetValue(KratosMultiphysics.DELTA_TIME, 0.25)

        process = GMA.SetMultipleMovingLoadsProcess(self.model_part, parameters)
        cond = self.compute_model_part.GetCondition(2)

        # initialize and set load
        self.initialize_test(process)

        # initialise matrices
        lhs = KratosMultiphysics.Matrix(0, 0)
        rhs = KratosMultiphysics.Vector(0)

        # set load on node
        cond.CalculateLocalSystem(lhs, rhs, self.model_part.ProcessInfo)

        self.checkRHS(rhs, [0.0, -0.5, 0.0, -1.5])

        # move load
        self.next_solution_step(process)

        # check if interpolation is done correctly
        cond.CalculateLocalSystem(lhs, rhs, self.model_part.ProcessInfo)

        self.checkRHS(rhs, [0.0, -1.0, 0.0, -1.0])

    def test_SetMultipleMovingLoadsReverseGeomConfigurationNegative(self):
        """
        Tests a moving load on 1 condition element, where the nodes of the element are reversed compared to the
        direction of the moving load, including a negative configuration along line condition direction in velocity direction.
        Returns
        -------

        """

        self.create_two_nodes_condition(reverse_nodes=True)

        parameters = self.base_parameters
        parameters.AddVector("configuration", [-0.25])

        self.model_part.ProcessInfo.SetValue(KratosMultiphysics.TIME, 0)
        self.model_part.ProcessInfo.SetValue(KratosMultiphysics.DELTA_TIME, 0.25)

        process = GMA.SetMultipleMovingLoadsProcess(self.model_part, parameters)
        cond = self.compute_model_part.GetCondition(2)

        # initialize and set load
        self.initialize_test(process)

        # initialise matrices
        lhs = KratosMultiphysics.Matrix(0, 0)
        rhs = KratosMultiphysics.Vector(0)

        # set load on node
        # cond.SetValue(SMA.MOVING_LOAD_LOCAL_DISTANCE, 0)
        cond.CalculateLocalSystem(lhs, rhs, self.model_part.ProcessInfo)

        self.checkRHS(rhs, [0.0, 0.0, 0.0, 0.0])

        # move load
        self.next_solution_step(process)

        # check if interpolation is done correctly
        cond.CalculateLocalSystem(lhs, rhs, self.model_part.ProcessInfo)

        self.checkRHS(rhs, [0.0, 0.0, 0.0, -2.0])

    def test_SetMultipleMovingLoadsReverseGeomConfigurationCombined(self):
        """
        Tests a moving load on 1 condition element, where the nodes of the element are reversed compared to the
        direction of the moving load, including a negative configuration along line condition direction in velocity direction.
        Returns
        -------

        """

        self.create_two_nodes_condition(reverse_nodes=True)

        parameters = self.base_parameters
        parameters.AddVector("configuration", [-0.25, 0, 0.25])

        self.model_part.ProcessInfo.SetValue(KratosMultiphysics.TIME, 0)
        self.model_part.ProcessInfo.SetValue(KratosMultiphysics.DELTA_TIME, 0.25)

        process = GMA.SetMultipleMovingLoadsProcess(self.model_part, parameters)

        conditions = []
        conditions.append(self.compute_model_part.GetCondition(2))
        conditions.append(self.compute_model_part.GetCondition(3))
        conditions.append(self.compute_model_part.GetCondition(4))

        # initialize and set load
        self.initialize_test(process)

        # set load on node
        all_rhs = []
        for cond in conditions:
            # initialise matrices
            lhs = KratosMultiphysics.Matrix(0, 0)
            rhs = KratosMultiphysics.Vector(0)
            cond.CalculateLocalSystem(lhs, rhs, self.model_part.ProcessInfo)
            all_rhs.append(list(rhs))

        self.checkRHS(all_rhs[0], [0.0, 0.0, 0.0, 0.0])
        self.checkRHS(all_rhs[1], [0.0, 0.0, 0.0, -2.0])
        self.checkRHS(all_rhs[2], [0.0, -0.5, 0.0, -1.5])

        # move load within first element
        self.next_solution_step(process)

        # check if interpolation is done correctly
        all_rhs = []
        for cond in conditions:
            lhs = KratosMultiphysics.Matrix(0, 0)
            rhs = KratosMultiphysics.Vector(0)
            cond.CalculateLocalSystem(lhs, rhs, self.model_part.ProcessInfo)
            all_rhs.append(list(rhs))

        self.checkRHS(all_rhs[0], [0.0, 0.0, 0.0, -2.0])
        self.checkRHS(all_rhs[1], [0.0, -0.5, 0.0, -1.5])
        self.checkRHS(all_rhs[2], [0.0, -1.0, 0.0, -1.0])

    def test_SetMultipleMovingLoadsMultipleConditions(self):
        """
        Tests a moving load on 2 condition elements, where the order of the elements is sorted in the direction of the
        moving load
        Returns
        -------

        """

        self.create_three_nodes_condition()

        parameters = self.base_parameters

        self.model_part.ProcessInfo.SetValue(KratosMultiphysics.TIME, 0)
        self.model_part.ProcessInfo.SetValue(KratosMultiphysics.DELTA_TIME, 0.5)
        process = GMA.SetMultipleMovingLoadsProcess(self.model_part, parameters)

        # get conditions
        conditions = [self.compute_model_part.GetCondition(3), self.compute_model_part.GetCondition(4)]

        # initialize and set load
        self.initialize_test(process)

        # initialise matrices
        lhs = KratosMultiphysics.Matrix(0, 0)
        rhs = KratosMultiphysics.Vector(0)

        # set load on node
        all_rhs = []
        for cond in conditions:
            cond.CalculateLocalSystem(lhs, rhs, self.model_part.ProcessInfo)
            all_rhs.append(list(rhs))

        self.checkRHS(all_rhs[0], [0.0, -2.0, 0.0, 0.0])
        self.checkRHS(all_rhs[1], [0.0, 0.0, 0.0, 0.0])

        # move load within first element
        self.next_solution_step(process)

        # check if interpolation is done correctly
        all_rhs = []
        for cond in conditions:
            cond.CalculateLocalSystem(lhs, rhs, self.model_part.ProcessInfo)
            all_rhs.append(list(rhs))

        self.checkRHS(all_rhs[0], [0.0, -1.0, 0.0, -1.0])
        self.checkRHS(all_rhs[1], [0.0, 0.0, 0.0, 0.0])

        # move load to element connection element
        self.next_solution_step(process)

        all_rhs = []
        for cond in conditions:
            cond.CalculateLocalSystem(lhs, rhs, self.model_part.ProcessInfo)
            all_rhs.append(list(rhs))

        self.checkRHS(all_rhs[0], [0.0, 0.0, 0.0, -2.0])
        self.checkRHS(all_rhs[1], [0.0, 0.0, 0.0, 0.0])

        # move load to next element
        self.next_solution_step(process)

        all_rhs = []
        for cond in conditions:
            cond.CalculateLocalSystem(lhs, rhs, self.model_part.ProcessInfo)
            all_rhs.append(list(rhs))

        self.checkRHS(all_rhs[0], [0.0, 0.0, 0.0, 0.0])
        self.checkRHS(all_rhs[1], [0.0, -1.0, 0.0, -1.0])

    def test_SetMultipleMovingLoadsMultipleConditionsConfigurationPositive(self):
        """
        Tests a moving load on 2 condition elements, where the order of the elements is sorted in the direction of the
        moving load
        Returns
        -------

        """

        self.create_three_nodes_condition()

        parameters = self.base_parameters
        parameters.AddVector("configuration", [0.5])

        self.model_part.ProcessInfo.SetValue(KratosMultiphysics.TIME, 0)
        self.model_part.ProcessInfo.SetValue(KratosMultiphysics.DELTA_TIME, 0.5)
        process = GMA.SetMultipleMovingLoadsProcess(self.model_part, parameters)

        # get conditions
        conditions = [self.compute_model_part.GetCondition(3), self.compute_model_part.GetCondition(4)]

        # initialize and set load
        self.initialize_test(process)

        # initialise matrices
        lhs = KratosMultiphysics.Matrix(0, 0)
        rhs = KratosMultiphysics.Vector(0)

        # set load on node

        all_rhs = []
        for cond in conditions:
            cond.CalculateLocalSystem(lhs, rhs, self.model_part.ProcessInfo)
            all_rhs.append(list(rhs))

        self.checkRHS(all_rhs[0], [0.0, -1.0, 0.0, -1.0])
        self.checkRHS(all_rhs[1], [0.0, 0.0, 0.0, 0.0])

        # move load within first element
        self.next_solution_step(process)

        # check if interpolation is done correctly
        all_rhs = []
        for cond in conditions:
            cond.CalculateLocalSystem(lhs, rhs, self.model_part.ProcessInfo)
            all_rhs.append(list(rhs))

        self.checkRHS(all_rhs[0], [0.0, 0.0, 0.0, -2.0])
        self.checkRHS(all_rhs[1], [0.0, 0.0, 0.0, 0.0])

        # move load to element connection element
        self.next_solution_step(process)

        all_rhs = []
        for cond in conditions:
            cond.CalculateLocalSystem(lhs, rhs, self.model_part.ProcessInfo)
            all_rhs.append(list(rhs))

        self.checkRHS(all_rhs[0], [0.0, 0.0, 0.0, 0.0])
        self.checkRHS(all_rhs[1], [0.0, -1.0, 0.0, -1.0])

        # move load to next element
        self.next_solution_step(process)

        all_rhs = []
        for cond in conditions:
            cond.CalculateLocalSystem(lhs, rhs, self.model_part.ProcessInfo)
            all_rhs.append(list(rhs))

        self.checkRHS(all_rhs[0], [0.0, 0.0, 0.0, 0.0])
        self.checkRHS(all_rhs[1], [0.0, 0.0, 0.0, -2.0])

    def test_SetMultipleMovingLoadsMultipleConditionsConfigurationNegative(self):
        """
        Tests a moving load on 2 condition elements, where the order of the elements is sorted in the direction of the
        moving load
        Returns
        -------

        """

        self.create_three_nodes_condition()

        parameters = self.base_parameters
        parameters.AddVector("configuration", [-0.5])

        self.model_part.ProcessInfo.SetValue(KratosMultiphysics.TIME, 0)
        self.model_part.ProcessInfo.SetValue(KratosMultiphysics.DELTA_TIME, 0.5)
        process = GMA.SetMultipleMovingLoadsProcess(self.model_part, parameters)

        # get conditions
        conditions = [self.compute_model_part.GetCondition(3), self.compute_model_part.GetCondition(4)]

        # initialize and set load
        self.initialize_test(process)

        # initialise matrices
        lhs = KratosMultiphysics.Matrix(0, 0)
        rhs = KratosMultiphysics.Vector(0)

        # set load on node

        all_rhs = []
        for cond in conditions:
            cond.CalculateLocalSystem(lhs, rhs, self.model_part.ProcessInfo)
            all_rhs.append(list(rhs))

        self.checkRHS(all_rhs[0], [0.0, 0.0, 0.0, 0.0])
        self.checkRHS(all_rhs[1], [0.0, 0.0, 0.0, 0.0])

        # move load within first element
        self.next_solution_step(process)

        # check if interpolation is done correctly
        all_rhs = []
        for cond in conditions:
            cond.CalculateLocalSystem(lhs, rhs, self.model_part.ProcessInfo)
            all_rhs.append(list(rhs))

        self.checkRHS(all_rhs[0], [0.0, -2.0, 0.0, 0.0])
        self.checkRHS(all_rhs[1], [0.0, 0.0, 0.0, 0.0])

        # move load to element connection element
        self.next_solution_step(process)

        all_rhs = []
        for cond in conditions:
            cond.CalculateLocalSystem(lhs, rhs, self.model_part.ProcessInfo)
            all_rhs.append(list(rhs))

        self.checkRHS(all_rhs[0], [0.0, -1.0, 0.0, -1.0])
        self.checkRHS(all_rhs[1], [0.0, 0.0, 0.0, 0.0])

        # move load to next element
        self.next_solution_step(process)

        all_rhs = []
        for cond in conditions:
            cond.CalculateLocalSystem(lhs, rhs, self.model_part.ProcessInfo)
            all_rhs.append(list(rhs))

        self.checkRHS(all_rhs[0], [0.0, 0.0, 0.0, -2.0])
        self.checkRHS(all_rhs[1], [0.0, 0.0, 0.0, 0.0])

        # move load to next element
        self.next_solution_step(process)

        all_rhs = []
        for cond in conditions:
            cond.CalculateLocalSystem(lhs, rhs, self.model_part.ProcessInfo)
            all_rhs.append(list(rhs))

        self.checkRHS(all_rhs[0], [0.0, 0.0, 0.0, 0.0])
        self.checkRHS(all_rhs[1], [0.0, -1.0, 0.0, -1.0])

    def test_SetMultipleMovingLoadsMultipleConditionsReversed(self):
        """
        Tests a moving load on 2 condition elements, where the order of the elements is reversed compared to the moving
        direction of the load
        Returns
        -------

        """

        self.create_three_nodes_condition(reverse_nodes=True)

        parameters = self.base_parameters

        self.model_part.ProcessInfo.SetValue(KratosMultiphysics.TIME, 0)
        self.model_part.ProcessInfo.SetValue(KratosMultiphysics.DELTA_TIME, 0.5)
        process = GMA.SetMultipleMovingLoadsProcess(self.model_part, parameters)

        # get conditions
        conditions = [self.compute_model_part.GetCondition(3), self.compute_model_part.GetCondition(4)]

        # initialize and set load
        self.initialize_test(process)

        # initialise matrices
        lhs = KratosMultiphysics.Matrix(0, 0)
        rhs = KratosMultiphysics.Vector(0)

        # calculate load
        all_rhs = []
        for cond in conditions:
            cond.CalculateLocalSystem(lhs, rhs, self.model_part.ProcessInfo)
            all_rhs.append(list(rhs))

        self.checkRHS(all_rhs[0], [0.0, 0.0, 0.0, 0.0])
        self.checkRHS(all_rhs[1], [0.0, 0.0, 0.0, -2.0])

        # move load within first element
        self.next_solution_step(process)

        # check if interpolation is done correctly
        all_rhs = []
        for cond in conditions:
            cond.CalculateLocalSystem(lhs, rhs, self.model_part.ProcessInfo)
            all_rhs.append(list(rhs))

        self.checkRHS(all_rhs[0], [0.0, 0.0, 0.0, 0.0])
        self.checkRHS(all_rhs[1], [0.0, -1.0, 0.0, -1.0])

        # move load to element connection element
        self.next_solution_step(process)

        # calculate load
        all_rhs = []
        for cond in conditions:
            cond.CalculateLocalSystem(lhs, rhs, self.model_part.ProcessInfo)
            all_rhs.append(list(rhs))

        self.checkRHS(all_rhs[0], [0.0, 0.0, 0.0, -2.0])
        self.checkRHS(all_rhs[1], [0.0, 0.0, 0.0, 0.0])

        # move load to next element, also increase time step
        self.model_part.ProcessInfo.SetValue(KratosMultiphysics.DELTA_TIME, 0.75)

        self.next_solution_step(process)

        # calculate load
        all_rhs = []
        for cond in conditions:
            cond.CalculateLocalSystem(lhs, rhs, self.model_part.ProcessInfo)
            all_rhs.append(list(rhs))

        self.checkRHS(all_rhs[0], [0.0, -1.5, 0.0, -0.5])
        self.checkRHS(all_rhs[1], [0.0, 0.0, 0.0, 0.0])

    def test_SetMultipleMovingLoadsMultipleConditionsReversedConfigurationPositive(self):
        """
        Tests a moving load on 2 condition elements, where the order of the elements is reversed compared to the moving
        direction of the load, including a positive configuration along line condition direction in velocity direction.
        Returns
        -------

        """

        self.create_three_nodes_condition(reverse_nodes=True)

        parameters = self.base_parameters
        parameters.AddVector("configuration", [0.5])

        self.model_part.ProcessInfo.SetValue(KratosMultiphysics.TIME, 0)
        self.model_part.ProcessInfo.SetValue(KratosMultiphysics.DELTA_TIME, 0.5)
        process = GMA.SetMultipleMovingLoadsProcess(self.model_part, parameters)

        # get conditions
        conditions = [self.compute_model_part.GetCondition(3), self.compute_model_part.GetCondition(4)]

        # initialize and set load
        self.initialize_test(process)

        # initialise matrices
        lhs = KratosMultiphysics.Matrix(0, 0)
        rhs = KratosMultiphysics.Vector(0)

        # calculate load
        all_rhs = []
        for cond in conditions:
            cond.CalculateLocalSystem(lhs, rhs, self.model_part.ProcessInfo)
            all_rhs.append(list(rhs))

        self.checkRHS(all_rhs[0], [0.0, 0.0, 0.0, 0.0])
        self.checkRHS(all_rhs[1], [0.0, -1.0, 0.0, -1.0])

        # move load within first element
        self.next_solution_step(process)

        # check if interpolation is done correctly
        all_rhs = []
        for cond in conditions:
            cond.CalculateLocalSystem(lhs, rhs, self.model_part.ProcessInfo)
            all_rhs.append(list(rhs))

        self.checkRHS(all_rhs[0], [0.0, 0.0, 0.0, -2.0])
        self.checkRHS(all_rhs[1], [0.0, 0.0, 0.0, 0.0])

        # move load to element connection element
        self.next_solution_step(process)

        # calculate load
        all_rhs = []
        for cond in conditions:
            cond.CalculateLocalSystem(lhs, rhs, self.model_part.ProcessInfo)
            all_rhs.append(list(rhs))

        self.checkRHS(all_rhs[0], [0.0, -1.0, 0.0, -1.0])
        self.checkRHS(all_rhs[1], [0.0, 0.0, 0.0, 0.0])

        # move load to next element, also increase time step
        self.model_part.ProcessInfo.SetValue(KratosMultiphysics.DELTA_TIME, 0.75)

        self.next_solution_step(process)

        # calculate load
        all_rhs = []
        for cond in conditions:
            cond.CalculateLocalSystem(lhs, rhs, self.model_part.ProcessInfo)
            all_rhs.append(list(rhs))

        self.checkRHS(all_rhs[0], [0.0, 0.0, 0.0, 0.0])
        self.checkRHS(all_rhs[1], [0.0, 0.0, 0.0, 0.0])

    def test_SetMultipleMovingLoadsMultipleConditionsReversedConfigurationNegative(self):
        """
        Tests a moving load on 2 condition elements, where the order of the elements is reversed compared to the moving
        direction of the load, including a negative configuration along line condition direction in velocity direction.
        Returns
        -------

        """

        self.create_three_nodes_condition(reverse_nodes=True)

        # set parameters and process info
        parameters = self.base_parameters
        parameters.AddVector("configuration", [-0.5])

        self.model_part.ProcessInfo.SetValue(KratosMultiphysics.TIME, 0)
        self.model_part.ProcessInfo.SetValue(KratosMultiphysics.DELTA_TIME, 0.5)
        process = GMA.SetMultipleMovingLoadsProcess(self.model_part, parameters)

        # get conditions
        conditions = [self.compute_model_part.GetCondition(3), self.compute_model_part.GetCondition(4)]

        # initialize and set load
        self.initialize_test(process)

        # initialise matrices
        lhs = KratosMultiphysics.Matrix(0, 0)
        rhs = KratosMultiphysics.Vector(0)

        # calculate load
        all_rhs = []
        for cond in conditions:
            cond.CalculateLocalSystem(lhs, rhs, self.model_part.ProcessInfo)
            all_rhs.append(list(rhs))

        self.checkRHS(all_rhs[0], [0.0, 0.0, 0.0, 0.0])
        self.checkRHS(all_rhs[1], [0.0, 0.0, 0.0, 0.0])

        # move load within first element
        self.next_solution_step(process)

        # check if interpolation is done correctly
        all_rhs = []
        for cond in conditions:
            cond.CalculateLocalSystem(lhs, rhs, self.model_part.ProcessInfo)
            all_rhs.append(list(rhs))

        self.checkRHS(all_rhs[0], [0.0, 0.0, 0.0, 0.0])
        self.checkRHS(all_rhs[1], [0.0, 0.0, 0.0, -2.0])

        # move load to element connection element
        self.next_solution_step(process)

        # calculate load
        all_rhs = []
        for cond in conditions:
            cond.CalculateLocalSystem(lhs, rhs, self.model_part.ProcessInfo)
            all_rhs.append(list(rhs))

        self.checkRHS(all_rhs[0], [0.0, 0.0, 0.0, 0.0])
        self.checkRHS(all_rhs[1], [0.0, -1.0, 0.0, -1.0])

        # move load to next element, also increase time step
        self.model_part.ProcessInfo.SetValue(KratosMultiphysics.DELTA_TIME, 0.75)
        self.next_solution_step(process)

        # calculate load
        all_rhs = []
        for cond in conditions:
            cond.CalculateLocalSystem(lhs, rhs, self.model_part.ProcessInfo)
            all_rhs.append(list(rhs))

        self.checkRHS(all_rhs[0], [0.0, -0.5, 0.0, -1.5])
        self.checkRHS(all_rhs[1], [0.0, 0.0, 0.0, 0.0])

    def test_SetMultipleMovingLoadsMultipleConditionsDifferentOrigin(self):
        """
        Tests a moving load on 2 condition elements, where the order of the elements is sorted in the direction of the
        moving load
        Returns
        -------

        """

        self.create_three_nodes_condition()

        parameters = self.base_parameters
        parameters.AddVector("origin", [1.25, 0, 0])

        process = GMA.SetMultipleMovingLoadsProcess(self.model_part, parameters)

        # get conditions
        conditions = [self.compute_model_part.GetCondition(3), self.compute_model_part.GetCondition(4)]

        # initialize and set load
        self.initialize_test(process)

        # initialise matrices
        lhs = KratosMultiphysics.Matrix(0, 0)
        rhs = KratosMultiphysics.Vector(0)

        # set load on node
        all_rhs = []
        for cond in conditions:
            cond.CalculateLocalSystem(lhs, rhs, self.model_part.ProcessInfo)
            all_rhs.append(list(rhs))

        self.checkRHS(all_rhs[0], [0.0, 0.0, 0.0, 0.0])
        self.checkRHS(all_rhs[1], [0.0, -1.5, 0.0, -0.5])

    def test_SetMultipleMovingLoadsMultipleConditionsDifferentOriginConfigurationPositive(self):
        """
        Tests a moving load on 2 condition elements, where the order of the elements is sorted in the direction of the
        moving load, including a positive configuration along line condition direction in velocity direction.
        Returns
        -------

        """

        self.create_three_nodes_condition()

        parameters = self.base_parameters
        parameters.AddVector("origin", [1.25, 0, 0])
        parameters.AddVector("configuration", [0.75])

        process = GMA.SetMultipleMovingLoadsProcess(self.model_part, parameters)

        # get conditions
        conditions = [self.compute_model_part.GetCondition(3), self.compute_model_part.GetCondition(4)]

        # initialize and set load
        self.initialize_test(process)

        # initialise matrices
        lhs = KratosMultiphysics.Matrix(0, 0)
        rhs = KratosMultiphysics.Vector(0)

        # set load on node
        all_rhs = []
        for cond in conditions:
            cond.CalculateLocalSystem(lhs, rhs, self.model_part.ProcessInfo)
            all_rhs.append(list(rhs))

        self.checkRHS(all_rhs[0], [0.0, 0.0, 0.0, 0.0])
        self.checkRHS(all_rhs[1], [0.0, 0.0, 0.0, -2.0])

    def test_SetMultipleMovingLoadsMultipleConditionsDifferentOriginConfigurationNegative(self):
        """
        Tests a moving load on 2 condition elements, where the order of the elements is sorted in the direction of the
        moving load, including a negative configuration along line condition direction in velocity direction.
        Returns
        -------

        """

        self.create_three_nodes_condition()

        parameters = self.base_parameters
        parameters.AddVector("origin", [1.25, 0, 0])
        parameters.AddVector("configuration", [-0.25])

        process = GMA.SetMultipleMovingLoadsProcess(self.model_part, parameters)

        # get conditions
        conditions = [self.compute_model_part.GetCondition(3), self.compute_model_part.GetCondition(4)]

        # initialize and set load
        self.initialize_test(process)

        # initialise matrices
        lhs = KratosMultiphysics.Matrix(0, 0)
        rhs = KratosMultiphysics.Vector(0)

        # set load on node
        all_rhs = []
        for cond in conditions:
            cond.CalculateLocalSystem(lhs, rhs, self.model_part.ProcessInfo)
            all_rhs.append(list(rhs))

        self.checkRHS(all_rhs[0], [0.0, 0.0, 0.0, -2.0])
        self.checkRHS(all_rhs[1], [0.0, 0.0, 0.0, 0.0])

    def test_SetMultipleMovingLoadsMultipleConditionsDifferentOriginReversed(self):
        """
        Tests a moving load on 2 condition elements, where the order of the elements is sorted in the direction of the
        moving load
        Returns
        -------

        """

        self.create_three_nodes_condition(reverse_nodes=True)

        parameters = self.base_parameters
        parameters.AddVector("origin", [1.25, 0, 0])

        process = GMA.SetMultipleMovingLoadsProcess(self.model_part, parameters)

        # get conditions
        conditions = [self.compute_model_part.GetCondition(3), self.compute_model_part.GetCondition(4)]

        # initialize and set load
        self.initialize_test(process)

        # initialise matrices
        lhs = KratosMultiphysics.Matrix(0, 0)
        rhs = KratosMultiphysics.Vector(0)

        # set load on node
        all_rhs = []
        for cond in conditions:
            cond.CalculateLocalSystem(lhs, rhs, self.model_part.ProcessInfo)
            all_rhs.append(list(rhs))

        self.checkRHS(all_rhs[0], [0.0, -0.5, 0.0, -1.5])
        self.checkRHS(all_rhs[1], [0.0, 0.0, 0.0, 0.0])

    def test_SetMultipleMovingLoadsMultipleConditionsDifferentOriginReversedConfigurationPositive(self):
        """
        Tests a moving load on 2 condition elements, where the order of the elements is sorted in the direction of the
        moving load, including a positive configuration along line condition direction in velocity direction.
        Returns
        -------

        """

        self.create_three_nodes_condition(reverse_nodes=True)

        parameters = self.base_parameters
        parameters.AddVector("origin", [1.25, 0, 0])
        parameters.AddVector("configuration", [0.25])

        process = GMA.SetMultipleMovingLoadsProcess(self.model_part, parameters)

        # get conditions
        conditions = [self.compute_model_part.GetCondition(3), self.compute_model_part.GetCondition(4)]

        # initialize and set load
        self.initialize_test(process)

        # initialise matrices
        lhs = KratosMultiphysics.Matrix(0, 0)
        rhs = KratosMultiphysics.Vector(0)

        # set load on node
        all_rhs = []
        for cond in conditions:
            cond.CalculateLocalSystem(lhs, rhs, self.model_part.ProcessInfo)
            all_rhs.append(list(rhs))

        self.checkRHS(all_rhs[0], [0.0, -1.0, 0.0, -1.0])
        self.checkRHS(all_rhs[1], [0.0, 0.0, 0.0, 0.0])

    def test_SetMultipleMovingLoadsMultipleConditionsDifferentOriginReversedConfigurationNegative(self):
        """
        Tests a moving load on 2 condition elements, where the order of the elements is sorted in the direction of the
        moving load, including a negative configuration along line condition direction in velocity direction.
        Returns
        -------

        """

        self.create_three_nodes_condition(reverse_nodes=True)

        parameters = self.base_parameters
        parameters.AddVector("origin", [1.25, 0, 0])
        parameters.AddVector("configuration", [-1.25])

        process = GMA.SetMultipleMovingLoadsProcess(self.model_part, parameters)

        # get conditions
        conditions = [self.compute_model_part.GetCondition(3), self.compute_model_part.GetCondition(4)]

        # initialize and set load
        self.initialize_test(process)

        # initialise matrices
        lhs = KratosMultiphysics.Matrix(0, 0)
        rhs = KratosMultiphysics.Vector(0)

        # set load on node
        all_rhs = []
        for cond in conditions:
            cond.CalculateLocalSystem(lhs, rhs, self.model_part.ProcessInfo)
            all_rhs.append(list(rhs))

        self.checkRHS(all_rhs[0], [0.0, 0.0, 0.0, 0.0])
        self.checkRHS(all_rhs[1], [0.0, 0.0, 0.0, -2.0])

    def test_SetMultipleMovingLoadsWithLoadFunction(self):
        """
       Tests a moving load on a condition element, where the load is a function of time
       Returns
       -------

       """

        self.create_two_nodes_condition()

        parameters = self.base_parameters
        parameters.AddVector("origin", [0.5, 0, 0])
        parameters["load"] = KratosMultiphysics.Parameters(self.time_dependent_load)

        d_time = 0.25
        self.model_part.ProcessInfo.SetValue(KratosMultiphysics.TIME, 0)
        self.model_part.ProcessInfo.SetValue(KratosMultiphysics.DELTA_TIME, d_time)

        process = GMA.SetMultipleMovingLoadsProcess(self.model_part, parameters)

        # get conditions
        cond = self.compute_model_part.GetCondition(2)

        # initialise and set load
        self.initialize_test(process)

        # initialise matrices
        lhs = KratosMultiphysics.Matrix(0, 0)
        rhs = KratosMultiphysics.Vector(0)

        # set load on node
        cond.CalculateLocalSystem(lhs, rhs, self.model_part.ProcessInfo)

        self.checkRHS(rhs, [0.0, 0.0, 0.0, 0.0])

        # change time and recalculate load
        self.model_part.ProcessInfo.SetValue(KratosMultiphysics.TIME, d_time)
        self.next_solution_step(process)

        # check if interpolation is done correctly
        cond.CalculateLocalSystem(lhs, rhs, self.model_part.ProcessInfo)

        self.checkRHS(rhs, [0.0, -0.125, 0.0, -0.375])

    def test_SetMultipleMovingLoadsWithLoadFunctionConfigurationPositive(self):
        """
       Tests a moving load on a condition element, where the load is a function of time, including a positive configuration along line condition direction in velocity direction.
       Returns
       -------

       """

        self.create_two_nodes_condition()

        parameters = self.base_parameters
        parameters.AddVector("origin", [0.5, 0, 0])
        parameters["load"] = KratosMultiphysics.Parameters(self.time_dependent_load)
        parameters.AddVector("configuration", [0.25])

        d_time = 0.25
        self.model_part.ProcessInfo.SetValue(KratosMultiphysics.TIME, 0)
        self.model_part.ProcessInfo.SetValue(KratosMultiphysics.DELTA_TIME, d_time)

        process = GMA.SetMultipleMovingLoadsProcess(self.model_part, parameters)
        cond = self.compute_model_part.GetCondition(2)

        # initialise and set load
        self.initialize_test(process)

        # initialise matrices
        lhs = KratosMultiphysics.Matrix(0, 0)
        rhs = KratosMultiphysics.Vector(0)

        # set load on node
        cond.CalculateLocalSystem(lhs, rhs, self.model_part.ProcessInfo)

        self.checkRHS(rhs, [0.0, 0.0, 0.0, 0.0])

        # change time and recalculate load
        self.model_part.ProcessInfo.SetValue(KratosMultiphysics.TIME, d_time)
        self.next_solution_step(process)

        # check if interpolation is done correctly
        cond.CalculateLocalSystem(lhs, rhs, self.model_part.ProcessInfo)

        self.checkRHS(rhs, [0.0, 0.0, 0.0, -0.5])

    def test_SetMultipleMovingLoadsWithLoadFunctionConfigurationNegative(self):
        """
       Tests a moving load on a condition element, where the load is a function of time, including a negative configuration along line condition direction in velocity direction.
       Returns
       -------

       """

        self.create_two_nodes_condition()

        parameters = self.base_parameters
        parameters.AddVector("origin", [0.5, 0, 0])
        parameters["load"] = KratosMultiphysics.Parameters(self.time_dependent_load)
        parameters.AddVector("configuration", [-0.25])

        d_time = 0.25
        self.model_part.ProcessInfo.SetValue(KratosMultiphysics.TIME, 0)
        self.model_part.ProcessInfo.SetValue(KratosMultiphysics.DELTA_TIME, d_time)

        process = GMA.SetMultipleMovingLoadsProcess(self.model_part, parameters)
        cond = self.compute_model_part.GetCondition(2)

        # initialise and set load
        self.initialize_test(process)

        # initialise matrices
        lhs = KratosMultiphysics.Matrix(0, 0)
        rhs = KratosMultiphysics.Vector(0)

        # set load on node
        cond.CalculateLocalSystem(lhs, rhs, self.model_part.ProcessInfo)

        self.checkRHS(rhs, [0.0, 0.0, 0.0, 0.0])

        # change time and recalculate load
        self.model_part.ProcessInfo.SetValue(KratosMultiphysics.TIME, d_time)
        self.next_solution_step(process)

        # check if interpolation is done correctly
        cond.CalculateLocalSystem(lhs, rhs, self.model_part.ProcessInfo)

        self.checkRHS(rhs, [0.0, -0.25, 0.0, -0.25])

    def test_SetMultipleMovingLoadsWithVelocityFunction(self):
        """
        Tests a moving load on a condition element, where the load velocity is a function of time.

        Returns
        -------

        """

        self.create_two_nodes_condition()

        parameters = self.base_parameters
        parameters["velocity"] = KratosMultiphysics.Parameters(self.time_dependent_velocity)

        self.model_part.ProcessInfo.SetValue(KratosMultiphysics.TIME, 0)
        self.model_part.ProcessInfo.SetValue(KratosMultiphysics.DELTA_TIME, 0.25)

        process = GMA.SetMultipleMovingLoadsProcess(self.model_part, parameters)
        cond = self.compute_model_part.GetCondition(2)

        # initialise and set load
        self.initialize_test(process)

        # initialise matrices
        lhs = KratosMultiphysics.Matrix(0, 0)
        rhs = KratosMultiphysics.Vector(0)

        # set load on node
        cond.CalculateLocalSystem(lhs, rhs, self.model_part.ProcessInfo)

        self.checkRHS(rhs, [0.0, -2.0, 0.0, 0.0])

        # change time and recalculate load
        self.next_solution_step(process)

        # check if interpolation is done correctly
        cond.CalculateLocalSystem(lhs, rhs, self.model_part.ProcessInfo)

        self.checkRHS(rhs, [0.0, -2.0, 0.0, 0.0])

        self.model_part.ProcessInfo.SetValue(KratosMultiphysics.TIME, 0.5)

        self.next_solution_step(process)

        # check if interpolation is done correctly
        cond.CalculateLocalSystem(lhs, rhs, self.model_part.ProcessInfo)

        self.checkRHS(rhs, [0.0, -1.5, 0.0, -0.5])

    def test_SetMultipleMovingLoadsWithVelocityFunctionConfigurationPositive(self):
        """
       Tests a moving load on a condition element, where the load velocity is a function of time, including a positive configuration along line condition direction in velocity direction.

       Returns
       -------

       """

        self.create_two_nodes_condition()

        parameters = self.base_parameters
        parameters["velocity"] = KratosMultiphysics.Parameters(self.time_dependent_velocity)
        parameters.AddVector("configuration", [0.5])

        self.model_part.ProcessInfo.SetValue(KratosMultiphysics.TIME, 0)
        self.model_part.ProcessInfo.SetValue(KratosMultiphysics.DELTA_TIME, 0.25)

        process = GMA.SetMultipleMovingLoadsProcess(self.model_part, parameters)
        cond = self.compute_model_part.GetCondition(2)

        # initialise and set load
        self.initialize_test(process)

        # initialise matrices
        lhs = KratosMultiphysics.Matrix(0, 0)
        rhs = KratosMultiphysics.Vector(0)

        # set load on node
        cond.CalculateLocalSystem(lhs, rhs, self.model_part.ProcessInfo)

        self.checkRHS(rhs, [0.0, -1.0, 0.0, -1.0])

        # change time and recalculate load
        self.next_solution_step(process)

        # check if interpolation is done correctly
        cond.CalculateLocalSystem(lhs, rhs, self.model_part.ProcessInfo)

        self.checkRHS(rhs, [0.0, -1.0, 0.0, -1.0])

        self.model_part.ProcessInfo.SetValue(KratosMultiphysics.TIME, 0.5)

        self.next_solution_step(process)

        # check if interpolation is done correctly
        cond.CalculateLocalSystem(lhs, rhs, self.model_part.ProcessInfo)

        self.checkRHS(rhs, [0.0, -0.5, 0.0, -1.5])

    def test_SetMultipleMovingLoadsWithVelocityFunctionConfigurationNegative(self):
        """
       Tests a moving load on a condition element, where the load velocity is a function of time, including a negative configuration along line condition direction in velocity direction.

       Returns
       -------

       """

        self.create_two_nodes_condition()

        parameters = self.base_parameters
        parameters["velocity"] = KratosMultiphysics.Parameters(self.time_dependent_velocity)
        parameters.AddVector("configuration", [-0.25])

        self.model_part.ProcessInfo.SetValue(KratosMultiphysics.TIME, 0)
        self.model_part.ProcessInfo.SetValue(KratosMultiphysics.DELTA_TIME, 0.25)

        process = GMA.SetMultipleMovingLoadsProcess(self.model_part, parameters)
        cond = self.compute_model_part.GetCondition(2)

        # initialise and set load
        self.initialize_test(process)

        # initialise matrices
        lhs = KratosMultiphysics.Matrix(0, 0)
        rhs = KratosMultiphysics.Vector(0)

        # set load on node
        cond.CalculateLocalSystem(lhs, rhs, self.model_part.ProcessInfo)

        self.checkRHS(rhs, [0.0, 0.0, 0.0, 0.0])

        # change time and recalculate load
        self.next_solution_step(process)

        # check if interpolation is done correctly
        cond.CalculateLocalSystem(lhs, rhs, self.model_part.ProcessInfo)

        self.checkRHS(rhs, [0.0, 0.0, 0.0, 0.0])

        self.model_part.ProcessInfo.SetValue(KratosMultiphysics.TIME, 0.5)

        self.next_solution_step(process)

        # check if interpolation is done correctly
        cond.CalculateLocalSystem(lhs, rhs, self.model_part.ProcessInfo)

        self.checkRHS(rhs, [0.0, -2.0, 0.0, 0.0])


if __name__ == '__main__':
    KratosUnittest.main()

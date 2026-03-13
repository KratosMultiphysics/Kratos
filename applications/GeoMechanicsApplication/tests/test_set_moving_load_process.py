
import KratosMultiphysics

import KratosMultiphysics.GeoMechanicsApplication as KratosGeo
import KratosMultiphysics.StructuralMechanicsApplication as SMA
import KratosMultiphysics.KratosUnittest as KratosUnittest



class KratosGeoMechanicsSetMovingLoadProcessTests(KratosUnittest.TestCase):

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

        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ROTATION_X, mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ROTATION_Y, mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ROTATION_Z, mp)

        return strategy

    def test_SetMovingLoadWithMotionTypeTotal(self):
        """
        tests set moving load with 'motion_type' = total. Which means that the TOTAL_DISPLACEMENT and TOTAL_ROTATION
        at the position of the load is calculated. If 'motion_type' is set to 'base', the DISPLACEMENT and ROTATION
        at the position of the load is calculated.
        """

        current_model = KratosMultiphysics.Model()
        mp = current_model.CreateModelPart("solid_part")
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.ROTATION)
        mp.AddNodalSolutionStepVariable(KratosGeo.TOTAL_DISPLACEMENT)
        mp.AddNodalSolutionStepVariable(KratosGeo.TOTAL_ROTATION)

        # create nodes
        second_coord = [1, 0, 0.0]
        mp.CreateNewNode(1, 0.0, 0.0, 0.0)
        mp.CreateNewNode(2, second_coord[0], second_coord[1], 0.0)

        # create condition
        cond = mp.CreateNewCondition("MovingLoadCondition2D2N", 1, [1, 2], mp.GetProperties()[1])
        nodes = cond.GetNodes()
        nodes[0].SetSolutionStepValue(KratosGeo.TOTAL_DISPLACEMENT_Y, 1.0)
        nodes[1].SetSolutionStepValue(KratosGeo.TOTAL_DISPLACEMENT_Y, 2.0)

        nodes[0].SetSolutionStepValue(KratosGeo.TOTAL_ROTATION_Z, 0.0)
        nodes[1].SetSolutionStepValue(KratosGeo.TOTAL_ROTATION_Z, 1.0)

        strategy = self.setup_strategy(mp)

        parameters = KratosMultiphysics.Parameters("""
                {
                    "help"            : "This process applies a moving load condition belonging to a modelpart. The load moves over line elements.",
                    "model_part_name" : "please_specify_model_part_name",
                    "variable_name"   : "POINT_LOAD",
                    "load"            : [0.0, -2.0, 0.0],
                    "direction"       : [1,1,1],
                    "velocity"        : 1,
                    "origin"          : [0.5,0,0],
                    "motion_type"     : "total"
                }
                """
                                                         )
        mp.ProcessInfo.SetValue(KratosMultiphysics.TIME, 0)
        mp.ProcessInfo.SetValue(KratosMultiphysics.DELTA_TIME, 0.25)

        process = SMA.SetMovingLoadProcess(mp, parameters)

        # initialise and set load
        process.ExecuteInitialize()

        # set load on node
        process.ExecuteInitializeSolutionStep()
        strategy.InitializeSolutionStep()

        # this calculates the total displacement and rotation at the position of the load
        strategy.SolveSolutionStep()

        calculated_tot_displacement = cond.GetValue(KratosGeo.TOTAL_DISPLACEMENT)
        calculated_tot_rotation = cond.GetValue(KratosGeo.TOTAL_ROTATION)

        self.assertVectorAlmostEqual(calculated_tot_displacement, [0.0, 1.375, 0.0], 6)
        self.assertVectorAlmostEqual(calculated_tot_rotation, [0.0, 0.0, 1.25], 6)


if __name__ == '__main__':
    KratosUnittest.main()
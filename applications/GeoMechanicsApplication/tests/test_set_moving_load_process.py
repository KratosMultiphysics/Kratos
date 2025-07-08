import sys
sys.path.append(r"C:\software_development\Kratos2\bin\Debug")
sys.path.append(r"C:\software_development\Kratos2\bin\Debug\libs")

import KratosMultiphysics

import KratosMultiphysics.GeoMechanicsApplication as KratosGeo
import KratosMultiphysics.StructuralMechanicsApplication as SMA
import KratosMultiphysics.KratosUnittest as KratosUnittest



class TestGeoSetMovingLoadProcess(KratosUnittest.TestCase):

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

    def test_SetMovingLoadWithMotionTypeTotal(self):
        """
        tests set moving load with 'motion_type' = total
        """

        current_model = KratosMultiphysics.Model()
        mp = current_model.CreateModelPart("solid_part")
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        mp.AddNodalSolutionStepVariable(KratosGeo.TOTAL_DISPLACEMENT)
        mp.AddNodalSolutionStepVariable(KratosGeo.TOTAL_ROTATION)

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
                    "origin"          : [0,0,0],
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

        # Work in progress
        strategy.SolveSolutionStep()


if __name__ == '__main__':
    KratosUnittest.main()
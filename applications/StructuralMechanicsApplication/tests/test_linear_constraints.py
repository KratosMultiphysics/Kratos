from __future__ import print_function, absolute_import, division
import KratosMultiphysics
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import KratosMultiphysics.KratosUnittest as KratosUnittest
import math

def inner_prod(a,b):
    tmp = 0
    for i in range(len(a)):
        tmp += a[i]*b[i]
    return tmp

class TestLinearConstraints(KratosUnittest.TestCase):
    def setUp(self):
        pass
    
    def _add_variables(self,mp):
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.VOLUME_ACCELERATION)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL)

    def _apply_material_properties(self, mp, dim, small_strain = True):
        #define properties
        mp.GetProperties()[1].SetValue(KratosMultiphysics.YOUNG_MODULUS,210)
        mp.GetProperties()[1].SetValue(KratosMultiphysics.POISSON_RATIO,0.3)
        mp.GetProperties()[1].SetValue(KratosMultiphysics.THICKNESS,1.0)
        mp.GetProperties()[1].SetValue(KratosMultiphysics.DENSITY,1.0)

        g = [0,-10.0,0]
        mp.GetProperties()[1].SetValue(KratosMultiphysics.VOLUME_ACCELERATION,g)

        cl = StructuralMechanicsApplication.LinearElasticPlaneStrain2DLaw()
        mp.GetProperties()[1].SetValue(KratosMultiphysics.CONSTITUTIVE_LAW,cl)

    def _set_buffer(self,mp):
        buffer_size = 3
        mp.SetBufferSize(buffer_size)
        # Cycle the buffer. This sets all historical nodal solution step data to
        # the current value and initializes the time stepping in the process info.
        mp.ProcessInfo[KratosMultiphysics.DELTA_TIME] = 1.0
        delta_time = mp.ProcessInfo[KratosMultiphysics.DELTA_TIME]
        time = mp.ProcessInfo[KratosMultiphysics.TIME]
        step =-buffer_size
        time = time - delta_time * buffer_size
        mp.ProcessInfo.SetValue(KratosMultiphysics.TIME, time)
        for i in range(0, buffer_size):
            step = step + 1
            time = time + delta_time
            mp.ProcessInfo.SetValue(KratosMultiphysics.STEP, step)
            mp.CloneTimeStep(time)

    def _apply_BCs(self,mp,A,b):
        pass

    def _create_strategy(self, mp):
        #define a minimal newton raphson solver
        linear_solver = KratosMultiphysics.SkylineLUFactorizationSolver()
        builder_and_solver = KratosMultiphysics.ResidualBasedBlockBuilderAndSolver(linear_solver)
        scheme = KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme()
#        convergence_criterion = KratosMultiphysics.ResidualCriteria(1e-13,1e-14) #Parallelism problem, not sure where!
        convergence_criterion = KratosMultiphysics.DisplacementCriteria(1e-13,1e-14)
        convergence_criterion.SetEchoLevel(0)

        #max_iters = 1
        max_iters = 20
        compute_reactions = True
        reform_step_dofs = True
        move_mesh_flag = True
        strategy = KratosMultiphysics.ResidualBasedNewtonRaphsonStrategy(mp,
                                                                        scheme,
                                                                        linear_solver,
                                                                        convergence_criterion,
                                                                        builder_and_solver,
                                                                        max_iters,
                                                                        compute_reactions,
                                                                        reform_step_dofs,
                                                                        move_mesh_flag)
        strategy.SetEchoLevel(1)

        return strategy

    def _solve_with_strategy(self, strategy, step):
        strategy.Check()
        strategy.Initialize()
        strategy.InitializeSolutionStep()
        strategy.Predict()
        strategy.SolveSolutionStep()
        strategy.FinalizeSolutionStep()
   

    def test_constraints(self):
        dim = 2

        current_model = KratosMultiphysics.Model()
        mp = current_model.CreateModelPart("tl_solid_part")

        self._add_variables(mp)
        self._apply_material_properties(mp, dim, False)

        #
        #     3
        #  4     2  
        #     1
        # Create nodes
        n1 = mp.CreateNewNode(1, 0.0, 0.0, 0.0)
        n2 = mp.CreateNewNode(2, 2.0, 1.0, 0.0)
        n3 = mp.CreateNewNode(3, 0.0, 2.0, 0.0)
        n4 = mp.CreateNewNode(4, -2.0, 1.0, 0.0)

        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.REACTION_X,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Y, KratosMultiphysics.REACTION_Y,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Z, KratosMultiphysics.REACTION_Z,mp)

        mp.CreateNewElement("TotalLagrangianElement2D4N", 1, [1,2,3,4], mp.GetProperties()[1])

        #node 1 is completely fixed
        n1.Fix(KratosMultiphysics.DISPLACEMENT_X)
        n1.Fix(KratosMultiphysics.DISPLACEMENT_Y)

        ################# apply a constraint on node 2
        #node2 is at allowed to slide normally to normal2 (taken as normal to the side 1-2)
        dx = n2.X - n1.X
        dy = n2.Y - n1.Y
        normal_2 = KratosMultiphysics.Array3([dy,-dx,0])/math.sqrt(dx**2+dy**2)
        inner_prod(normal_2,normal_2)
        master_dofs_2 = [n2.GetDof(KratosMultiphysics.DISPLACEMENT_Y)]
        slave_dofs_2 = [n2.GetDof(KratosMultiphysics.DISPLACEMENT_X)]
        RelationMatrix2 = KratosMultiphysics.Matrix(1,1)
        RelationMatrix2[0,0] = -normal_2[1]/normal_2[0]
        ConstantVector = KratosMultiphysics.Vector([0.0])

        constraint_2 = KratosMultiphysics.LinearMasterSlaveConstraint(2,
 
            master_dofs_2, 
            slave_dofs_2, 
            RelationMatrix2, 
            ConstantVector)
        mp.AddMasterSlaveConstraint(constraint_2)

        ################# apply a slip constraint on node 4 
        # note that SlipConstraint is practically a placeholder for the operations done just above
        #node4 is at allowed to slide normally to normal4
        dx = n1.X - n4.X
        dy = n1.Y - n4.Y
        normal_4 = KratosMultiphysics.Array3([dy,-dx,0])/math.sqrt(dx**2+dy**2)
        n4.SetSolutionStepValue(KratosMultiphysics.NORMAL, normal_4) #artificially set the normal

        constraint_4 = KratosMultiphysics.SlipConstraint(4,
            n4,
            KratosMultiphysics.DISPLACEMENT_X,
            KratosMultiphysics.DISPLACEMENT_Y,
            KratosMultiphysics.DISPLACEMENT_Z,
            KratosMultiphysics.NORMAL)
        mp.AddMasterSlaveConstraint(constraint_4)

        #solve the problem
        strategy = self._create_strategy(mp)
        self._solve_with_strategy(strategy,0)

        ##verify the results
        d1 = n1.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT)
        d2 = n2.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT)
        d3 = n3.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT)
        d4 = n4.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT)

        self.assertAlmostEqual(d2[0]*normal_2[0] + d2[1]*normal_2[1],0.0, 15)
        self.assertAlmostEqual(d4[0]*normal_4[0] + d4[1]*normal_4[1],0.0, 15)
        self.assertAlmostEqual(d3[0], 0.0, 15) #symmetry condition       
        self.assertAlmostEqual(d3[1], 2.0*d2[1], 15)

        R2 = n2.GetSolutionStepValue(KratosMultiphysics.REACTION)
        R4 = n4.GetSolutionStepValue(KratosMultiphysics.REACTION)

        nR2 = inner_prod(normal_2, R2)
        nR4 = inner_prod(normal_4, R4)

        #check that Reactions are in compression
        self.assertTrue(nR2 < 0.0)
        self.assertTrue(nR4 < 0.0)
        
        #check that tangential component is zero
        tang_2 = R2 - nR2*normal_2
        tang_4 = R4 - nR4*normal_4

        self.assertTrue(tang_2.norm_2() < 1e-12)
        self.assertTrue(tang_4.norm_2() < 1e-12)

        self.assertEqual(mp.ProcessInfo[KratosMultiphysics.NL_ITERATION_NUMBER], 5) #4 if using Residual Criteria

    def __post_process(self, main_model_part):
        from gid_output_process import GiDOutputProcess
        self.gid_output = GiDOutputProcess(main_model_part,
                                    "gid_output",
                                    KratosMultiphysics.Parameters("""
                                        {
                                            "result_file_configuration" : {
                                                "gidpost_flags": {
                                                    "GiDPostMode": "GiD_PostBinary",
                                                    "WriteDeformedMeshFlag": "WriteUndeformed",
                                                    "WriteConditionsFlag": "WriteConditions",
                                                    "MultiFileFlag": "SingleFile"
                                                },
                                                "nodal_results"       : ["DISPLACEMENT"],
                                                "gauss_point_results" : ["GREEN_LAGRANGE_STRAIN_TENSOR","CAUCHY_STRESS_TENSOR"]
                                            }
                                        }
                                        """)
                                    )

        self.gid_output.ExecuteInitialize()
        self.gid_output.ExecuteBeforeSolutionLoop()
        self.gid_output.ExecuteInitializeSolutionStep()
        self.gid_output.PrintOutput()
        self.gid_output.ExecuteFinalizeSolutionStep()
        self.gid_output.ExecuteFinalize()

if __name__ == '__main__':
    KratosUnittest.main()

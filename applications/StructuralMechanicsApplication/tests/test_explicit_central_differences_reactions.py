import KratosMultiphysics

import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import KratosMultiphysics.KratosUnittest as KratosUnittest


class TestExplicitCentralDifferencesReactions(KratosUnittest.TestCase):
    def setUp(self):
        pass


    def _add_variables(self,mp):
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.ACCELERATION)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)
        mp.AddNodalSolutionStepVariable(StructuralMechanicsApplication.POINT_LOAD)

        mp.AddNodalSolutionStepVariable(StructuralMechanicsApplication.MIDDLE_VELOCITY)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_MASS)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.FORCE_RESIDUAL)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.RESIDUAL_VECTOR)


    def _add_dofs(self,mp):
        # Adding the dofs AND their corresponding reaction!
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.REACTION_X,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Y, KratosMultiphysics.REACTION_Y,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Z, KratosMultiphysics.REACTION_Z,mp)


    def _create_nodes(self,mp,element_name):
        mp.CreateNewNode(1,  0.0,  1.0,   0.0)
        mp.CreateNewNode(2,  2.0,  1.0,   0.0)
        mp.CreateNewNode(3,  0.0,  0.0,   0.0)
        mp.CreateNewNode(4,  2.0,  0.0,   0.0)
        mp.CreateNewNode(5,  4.0,  1.0,   0.0)
        mp.CreateNewNode(6,  4.0,  0.0,   0.0)


    def _create_elements(self,mp,element_name):
        mp.CreateNewElement(element_name, 1, [4,2,1,3], mp.GetProperties()[1])
        mp.CreateNewElement(element_name, 2, [6,5,2,4], mp.GetProperties()[1])


    def _apply_dirichlet_BCs(self,mp):
        KratosMultiphysics.VariableUtils().ApplyFixity(KratosMultiphysics.DISPLACEMENT_X, True, mp.Nodes)
        KratosMultiphysics.VariableUtils().ApplyFixity(KratosMultiphysics.DISPLACEMENT_Y, True, mp.Nodes)
        KratosMultiphysics.VariableUtils().ApplyFixity(KratosMultiphysics.DISPLACEMENT_Z, True, mp.Nodes)


    def _apply_neumann_BCs(self,mp):
        index = 1
        for node in mp.Nodes:
            node.SetSolutionStepValue(StructuralMechanicsApplication.POINT_LOAD,0,[0.5,0.0,0.0])
            mp.CreateNewCondition("PointLoadCondition3D1N",index,[node.Id],mp.GetProperties()[1])
            index += 1


    def _apply_material_properties(self,mp):
        #define properties
        mp.GetProperties()[1].SetValue(KratosMultiphysics.YOUNG_MODULUS,1000.0)
        mp.GetProperties()[1].SetValue(KratosMultiphysics.POISSON_RATIO,0.0)
        mp.GetProperties()[1].SetValue(KratosMultiphysics.THICKNESS,1.0)
        mp.GetProperties()[1].SetValue(KratosMultiphysics.DENSITY,100.0)

        cl = StructuralMechanicsApplication.LinearElasticPlaneStress2DLaw()

        mp.GetProperties()[1].SetValue(KratosMultiphysics.CONSTITUTIVE_LAW,cl)


    def _solve(self,mp):
        #define a minimal newton raphson solver
        linear_solver = KratosMultiphysics.SkylineLUFactorizationSolver()
        builder_and_solver = KratosMultiphysics.ResidualBasedBlockBuilderAndSolver(linear_solver)
        scheme = StructuralMechanicsApplication.ExplicitCentralDifferencesScheme(0.00,0.00,0.00)
        convergence_criterion = KratosMultiphysics.ResidualCriteria(1e-4,1e-9)

        max_iters = 20
        compute_reactions = True
        reform_step_dofs = True
        move_mesh_flag = False
        
        strategy = StructuralMechanicsApplication.MechanicalExplicitStrategy(mp,
                                                                        scheme,
                                                                        compute_reactions,
                                                                        reform_step_dofs,
                                                                        move_mesh_flag)
        
        strategy.SetEchoLevel(0)
        mp.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] = 2
        strategy.Check()

        #time integration parameters
        time = 0.0
        end_time = 1.0 
        step = 0
        dt = 0.001

        # Solve the problem
        while time <= end_time:
            time = time + dt
            step = step + 1
            mp.CloneTimeStep(time)
            mp.ProcessInfo[KratosMultiphysics.STEP] = step

            strategy.Solve()


    def _check_results(self,node, reaction_results):
        #check that the results are exact on the node
        reaction = node.GetSolutionStepValue(KratosMultiphysics.REACTION, 0)
        self.assertAlmostEqual(reaction[0], reaction_results[0])
        self.assertAlmostEqual(reaction[1], reaction_results[1])
        self.assertAlmostEqual(reaction[2], reaction_results[2])


    def execute_explicit_central_differences_test(self, current_model,element_name, reaction_results):
        mp = current_model.CreateModelPart("Plate")
        mp.SetBufferSize(2)

        self._add_variables(mp)
        self._apply_material_properties(mp)
        self._create_nodes(mp,element_name)
        self._add_dofs(mp)
        self._create_elements(mp,element_name)

        #create a submodelpart for dirichlet boundary conditions
        bcs_dirichlet = mp.CreateSubModelPart("BoundaryCondtionsDirichlet")
        bcs_dirichlet.AddNodes([1,3])

        #create a submodelpart for neumann boundary conditions
        bcs_neumann = mp.CreateSubModelPart("BoundaryCondtionsNeumann")
        bcs_neumann.AddNodes([5,6])

        self._apply_dirichlet_BCs(bcs_dirichlet)
        self._apply_neumann_BCs(bcs_neumann)
        self._solve(mp)

        self._check_results(mp.Nodes[1],reaction_results)


    def test_explicit_central_differences_reactions(self):
        element_name = "SmallDisplacementElement2D4N"
        reaction_results = [-0.18534577 ,0.0 ,0.0]

        current_model = KratosMultiphysics.Model()
        self.execute_explicit_central_differences_test(current_model,
                                element_name,
                                reaction_results) 

if __name__ == '__main__':
    KratosUnittest.main()
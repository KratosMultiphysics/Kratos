from __future__ import print_function, absolute_import, division
import KratosMultiphysics
from KratosMultiphysics import *

import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import KratosMultiphysics.KratosUnittest as KratosUnittest

class TestPatchTestShellsOrthotropic(KratosUnittest.TestCase):
    def setUp(self):
        pass


    def _add_variables(self,mp):
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.ROTATION)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION_MOMENT)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.VOLUME_ACCELERATION)
        mp.AddNodalSolutionStepVariable(StructuralMechanicsApplication.POINT_LOAD)


    def _add_dofs(self,mp):
        # Adding the dofs AND their corresponding reaction!
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.REACTION_X,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Y, KratosMultiphysics.REACTION_Y,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Z, KratosMultiphysics.REACTION_Z,mp)

        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ROTATION_X, KratosMultiphysics.REACTION_MOMENT_X,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ROTATION_Y, KratosMultiphysics.REACTION_MOMENT_Y,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ROTATION_Z, KratosMultiphysics.REACTION_MOMENT_Z,mp)


    def _create_nodes(self,mp,element_name):
        mp.CreateNewNode(1, -0.5, - 0.45,  0.1)
        mp.CreateNewNode(2,  0.7,  -0.5,   0.2)
        mp.CreateNewNode(3,  0.55,  0.6,   0.15)
        mp.CreateNewNode(4, -0.48,  0.65,  0.0)
        mp.CreateNewNode(5,  0.02, -0.01, -0.15)

        if element_name.endswith("4N"): # create aditional nodes needed for quad-setup
            mp.CreateNewNode(6, -0.03, -0.5,   0.0)
            mp.CreateNewNode(7,  0.51,  0.02,  0.03)
            mp.CreateNewNode(8, -0.01,  0.52, -0.05)
            mp.CreateNewNode(9, -0.49, -0.0,   0.0)


    def _create_elements(self,mp,element_name):
        if element_name.endswith("4N"): # Quadrilaterals
            mp.CreateNewElement(element_name, 1, [1,6,5,9], mp.GetProperties()[1])
            mp.CreateNewElement(element_name, 2, [6,2,7,5], mp.GetProperties()[1])
            mp.CreateNewElement(element_name, 3, [5,7,3,8], mp.GetProperties()[1])
            mp.CreateNewElement(element_name, 4, [9,5,8,4], mp.GetProperties()[1])
        else: # Triangles
            mp.CreateNewElement(element_name, 1, [1,2,5], mp.GetProperties()[1])
            mp.CreateNewElement(element_name, 2, [2,3,5], mp.GetProperties()[1])
            mp.CreateNewElement(element_name, 3, [3,4,5], mp.GetProperties()[1])
            mp.CreateNewElement(element_name, 4, [4,1,5], mp.GetProperties()[1])


    def _apply_dirichlet_BCs(self,mp):
        KratosMultiphysics.VariableUtils().ApplyFixity(KratosMultiphysics.DISPLACEMENT_X, True, mp.Nodes)
        KratosMultiphysics.VariableUtils().ApplyFixity(KratosMultiphysics.DISPLACEMENT_Y, True, mp.Nodes)
        KratosMultiphysics.VariableUtils().ApplyFixity(KratosMultiphysics.DISPLACEMENT_Z, True, mp.Nodes)
        KratosMultiphysics.VariableUtils().ApplyFixity(KratosMultiphysics.ROTATION_X, True, mp.Nodes)
        KratosMultiphysics.VariableUtils().ApplyFixity(KratosMultiphysics.ROTATION_Y, True, mp.Nodes)
        KratosMultiphysics.VariableUtils().ApplyFixity(KratosMultiphysics.ROTATION_Z, True, mp.Nodes)


    def _apply_neumann_BCs(self,mp):
        for node in mp.Nodes:
            node.SetSolutionStepValue(StructuralMechanicsApplication.POINT_LOAD,0,[6.1,-5.5,8.9])
            mp.CreateNewCondition("PointLoadCondition3D1N",1,[node.Id],mp.GetProperties()[1])


    def _apply_material_properties(self,mp):
        #define properties
        orthotropic_props = Matrix(4,16)

        # Orthotropic mechanical moduli
        orthotropic_props[0,0] = 0.5 #lamina thickness
        orthotropic_props[0,1] = 0.0 #lamina rotation (deg)
        orthotropic_props[0,2] = 7850 #density
        orthotropic_props[0,3] = 7500 #E1
        orthotropic_props[0,4] = 2000 #E2
        orthotropic_props[0,5] = 0.25 #nu_12
        orthotropic_props[0,6] = 1250 #G_12
        orthotropic_props[0,7] = 625 #G_13
        orthotropic_props[0,8] = 625 #G_23

        # Orthotropic mechanical strengths. (T)ensile, (C)ompression, (S)hear
        # along 1, 2, 3 lamina directions
        orthotropic_props[0,9] = 800 #T1
        orthotropic_props[0,10] = 500 #C1
        orthotropic_props[0,11] = 40 #T2
        orthotropic_props[0,12] = 300 #C2
        orthotropic_props[0,13] = 60 #S12
        orthotropic_props[0,14] = 60 #S13
        orthotropic_props[0,15] = 60 #S23

        for row in range(1,4):
            for col in range(16):
                orthotropic_props[row,col] = orthotropic_props[0,col]
        orthotropic_props[1,1] = 90
        orthotropic_props[2,1] = 90

        mp.GetProperties()[1].SetValue(
            KratosMultiphysics.StructuralMechanicsApplication.SHELL_ORTHOTROPIC_LAYERS,orthotropic_props)

        g = [0,0,0]
        mp.GetProperties()[1].SetValue(KratosMultiphysics.VOLUME_ACCELERATION,g)

        cl = StructuralMechanicsApplication.LinearElasticOrthotropic2DLaw()

        mp.GetProperties()[1].SetValue(KratosMultiphysics.CONSTITUTIVE_LAW,cl)


    def _solve(self,mp):
        #define a minimal newton raphson solver
        linear_solver = KratosMultiphysics.SkylineLUFactorizationSolver()
        builder_and_solver = KratosMultiphysics.ResidualBasedBlockBuilderAndSolver(linear_solver)
        scheme = KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme()
        convergence_criterion = KratosMultiphysics.ResidualCriteria(1e-14,1e-20)
        convergence_criterion.SetEchoLevel(0)

        max_iters = 20
        compute_reactions = True
        reform_step_dofs = True
        calculate_norm_dx = False
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
        strategy.SetEchoLevel(0)

        strategy.Check()
        strategy.Solve()


    def _check_results(self,node,displacement_results, rotation_results):
        ##check that the results are exact on the node
        disp = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT)
        self.assertAlmostEqual(disp[0], displacement_results[0], 10)
        self.assertAlmostEqual(disp[1], displacement_results[1], 10)
        self.assertAlmostEqual(disp[2], displacement_results[2], 10)

        rot = node.GetSolutionStepValue(KratosMultiphysics.ROTATION)
        self.assertAlmostEqual(rot[0], rotation_results[0], 10)
        self.assertAlmostEqual(rot[1], rotation_results[1], 10)
        self.assertAlmostEqual(rot[2], rotation_results[2], 10)


    def _check_results_stress(self,element,stress_variable,reference_stress_results,processInfo):
        ##check that the results are exact on the first gauss point
        ##only upper triangle of stresses are checked due to symmetry
        stress = element.CalculateOnIntegrationPoints(stress_variable, processInfo)[0]
        self.assertAlmostEqual(stress[0,0], reference_stress_results[0], 10)
        self.assertAlmostEqual(stress[0,1], reference_stress_results[1], 10)
        self.assertAlmostEqual(stress[0,2], reference_stress_results[2], 10)
        self.assertAlmostEqual(stress[1,1], reference_stress_results[3], 10)
        self.assertAlmostEqual(stress[1,2], reference_stress_results[4], 10)
        self.assertAlmostEqual(stress[2,2], reference_stress_results[5], 10)


    def execute_shell_test(self, current_model, element_name, displacement_results, rotation_results, shell_stress_top_surface_results, shell_stress_bottom_surface_results, tsai_wu_result,do_post_processing):
        mp = current_model.CreateModelPart("solid_part")
        mp.SetBufferSize(2)

        self._add_variables(mp)
        self._apply_material_properties(mp)
        self._create_nodes(mp,element_name)
        self._add_dofs(mp)
        self._create_elements(mp,element_name)

        #create a submodelpart for dirichlet boundary conditions
        bcs_dirichlet = mp.CreateSubModelPart("BoundaryCondtionsDirichlet")
        bcs_dirichlet.AddNodes([1,2,4])

        #create a submodelpart for neumann boundary conditions
        bcs_neumann = mp.CreateSubModelPart("BoundaryCondtionsNeumann")
        bcs_neumann.AddNodes([3])

        self._apply_dirichlet_BCs(bcs_dirichlet)
        self._apply_neumann_BCs(bcs_neumann)
        self._solve(mp)

        # Check displacements
        self._check_results(mp.Nodes[3],displacement_results, rotation_results)

        # Check stresses at each surface
        self._check_results_stress(mp.Elements[1],
                                   StructuralMechanicsApplication.SHELL_ORTHOTROPIC_STRESS_TOP_SURFACE,
                                   shell_stress_top_surface_results,mp.ProcessInfo)
        self._check_results_stress(mp.Elements[1],
                                   StructuralMechanicsApplication.SHELL_ORTHOTROPIC_STRESS_BOTTOM_SURFACE,
                                   shell_stress_bottom_surface_results,mp.ProcessInfo)

        # Check results of doubles on 2nd element @ Gauss Point [0] only
        self.assertAlmostEqual(mp.Elements[1].CalculateOnIntegrationPoints(StructuralMechanicsApplication.TSAI_WU_RESERVE_FACTOR,
                               mp.ProcessInfo)[0], tsai_wu_result, 9)

        if do_post_processing:
            self.__post_process(mp)


    def test_thin_shell_triangle(self):
        element_name = "ShellThinElementCorotational3D3N"
        displacement_results = [0.0028456068244 , -0.0021804536526 , 0.0014855251225]
        rotation_results     = [0.0028315743508 , -0.000450044246 , -0.0055701845132]
        shell_stress_top_surface_results    = [0.9088110489672 , -0.0570461205561 , 0.0 , 1.7678124328652 , 0.0 , 0.0]
        shell_stress_bottom_surface_results = [-0.4936295259123 , 0.2914348407351 , 0.0 , -0.5256560385672 , 0.0 , 0.0]
        tsai_wu_result = 39.6023549141987

        current_model = KratosMultiphysics.Model()
        self.execute_shell_test(current_model,
                                element_name, 
                                displacement_results, 
                                rotation_results,
                                shell_stress_top_surface_results,
                                shell_stress_bottom_surface_results,
                                tsai_wu_result,
                                False) # Do PostProcessing for GiD?


    def test_thick_shell_triangle(self):
        element_name = "ShellThickElementCorotational3D3N"
        displacement_results = [0.0004043490308 , -0.0016074440019 , 0.0092911008314]
        rotation_results     = [0.0021176894774 , -0.0005954288823 , -0.0015930914838]
        shell_stress_top_surface_results    = [3.4555559859345 , 3.6328430864296 , 0.2347447591457 , 0.1945591765769 , -1.5033148859134 , 0.0]
        shell_stress_bottom_surface_results = [-0.5442976284974 , -0.1011836349433 , 0.2347447591457 , -2.8139010064313 , -1.5033148859134 , 0.0]
        tsai_wu_result = 15.0065495746848

        current_model = KratosMultiphysics.Model()
        self.execute_shell_test(current_model,
                                element_name, 
                                displacement_results, 
                                rotation_results,
                                shell_stress_top_surface_results,
                                shell_stress_bottom_surface_results,
                                tsai_wu_result,
                                False) # Do PostProcessing for GiD?


    def test_thin_shell_quadrilateral(self):
        element_name = "ShellThinElementCorotational3D4N"
        displacement_results = [0.0225804891311 , -0.0233155244988 , 0.0048050841112]
        rotation_results     = [0.0248341724156 , 0.0105468617083 , -0.0691658930497]
        shell_stress_top_surface_results    = [0.284184788186 , -12.2844786822622 , 0.0 , 4.3796427631839 , 0.0 , 0.0]
        shell_stress_bottom_surface_results = [10.340621141106 , 5.6934270260323 , 0.0 , -2.973608875272 , 0.0 , 0.0]
        tsai_wu_result = 3.828332205752

        current_model = KratosMultiphysics.Model()
        self.execute_shell_test(current_model,
                                element_name, 
                                displacement_results, 
                                rotation_results, 
                                shell_stress_top_surface_results,
                                shell_stress_bottom_surface_results,
                                tsai_wu_result,
                                False) # Do PostProcessing for GiD?


    def test_thick_shell_quadrilateral(self):
        element_name = "ShellThickElementCorotational3D4N"
        displacement_results = [0.0035689894826 , -0.0094851917758 , 0.0191734998621]
        rotation_results     = [0.009933211939  , 0.0006068078079  , -0.0174332051568]
        shell_stress_top_surface_results    = [-3.9178477532111 , -4.1074850572552 , -2.4426862077188 , 10.3723187292559 , 1.6354826554283 , 0.0]
        shell_stress_bottom_surface_results = [5.2113212123242 , -0.2324161069908 , -2.4426862077188 , -11.6664322521041 , 1.6354826554283 , 0.0]
        tsai_wu_result = 3.4966651118454

        current_model = KratosMultiphysics.Model()
        self.execute_shell_test(current_model,
                                element_name, 
                                displacement_results, 
                                rotation_results,
                                shell_stress_top_surface_results,
                                shell_stress_bottom_surface_results,
                                tsai_wu_result,
                                False) # Do PostProcessing for GiD?


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
                                                "nodal_results"       : ["DISPLACEMENT", "ROTATION", "POINT_LOAD"],
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

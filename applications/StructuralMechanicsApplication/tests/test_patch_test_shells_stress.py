from __future__ import print_function, absolute_import, division
import KratosMultiphysics 

import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import KratosMultiphysics.KratosUnittest as KratosUnittest


class TestPatchTestShellsStressRec(KratosUnittest.TestCase):
    def setUp(self):
        pass
    

    def _add_variables(self,mp):
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.ROTATION)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.TORQUE)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.VOLUME_ACCELERATION)  
        mp.AddNodalSolutionStepVariable(StructuralMechanicsApplication.POINT_LOAD)      
        
    
    def _add_dofs(self,mp):
        # Adding the dofs AND their corresponding reaction!
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.REACTION_X,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Y, KratosMultiphysics.REACTION_Y,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Z, KratosMultiphysics.REACTION_Z,mp)

        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ROTATION_X, KratosMultiphysics.TORQUE_X,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ROTATION_Y, KratosMultiphysics.TORQUE_Y,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ROTATION_Z, KratosMultiphysics.TORQUE_Z,mp)


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
        mp.GetProperties()[1].SetValue(KratosMultiphysics.YOUNG_MODULUS,100e3)
        mp.GetProperties()[1].SetValue(KratosMultiphysics.POISSON_RATIO,0.3)
        mp.GetProperties()[1].SetValue(KratosMultiphysics.THICKNESS,1.0)
        mp.GetProperties()[1].SetValue(KratosMultiphysics.DENSITY,1.0)
        
        g = [0,0,0]
        mp.GetProperties()[1].SetValue(KratosMultiphysics.VOLUME_ACCELERATION,g)
        
        cl = StructuralMechanicsApplication.LinearElasticPlaneStress2DLaw()

        mp.GetProperties()[1].SetValue(KratosMultiphysics.CONSTITUTIVE_LAW,cl) 
        

    def _solve(self,mp):
        #define a minimal newton raphson solver
        linear_solver = KratosMultiphysics.SkylineLUFactorizationSolver()
        builder_and_solver = KratosMultiphysics.ResidualBasedBlockBuilderAndSolver(linear_solver)
        scheme = KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme()
        convergence_criterion = KratosMultiphysics.ResidualCriteria(1e-14,1e-20)
        
        max_iters = 20
        compute_reactions = True
        reform_step_dofs = True
        calculate_norm_dx = False
        move_mesh_flag = True
        strategy = KratosMultiphysics.ResidualBasedLinearStrategy(mp, 
                                                                  scheme, 
                                                                  linear_solver, 
                                                                  builder_and_solver, 
                                                                  compute_reactions, 
                                                                  reform_step_dofs, 
                                                                  calculate_norm_dx,
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
        # print("stress:", round(stress[0,0], 13), ",", round(stress[0,1], 13), ",", round(stress[0,2], 13), ",", round(stress[1,1], 13), ",", round(stress[1,2], 13), ",", round(stress[2,2], 13))
        self.assertAlmostEqual(stress[0,0], reference_stress_results[0], 10)
        self.assertAlmostEqual(stress[0,1], reference_stress_results[1], 10)
        self.assertAlmostEqual(stress[0,2], reference_stress_results[2], 10)
        self.assertAlmostEqual(stress[1,1], reference_stress_results[3], 10)
        self.assertAlmostEqual(stress[1,2], reference_stress_results[4], 10)
        self.assertAlmostEqual(stress[2,2], reference_stress_results[5], 10)


    def execute_shell_test(self, element_name, displacement_results, rotation_results, shell_stress_middle_surface_results, shell_stress_top_surface_results, shell_stress_bottom_surface_results, shell_von_mises_result,do_post_processing):
        mp = KratosMultiphysics.ModelPart("solid_part")
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
                                   StructuralMechanicsApplication.SHELL_STRESS_MIDDLE_SURFACE,
                                   shell_stress_middle_surface_results,mp.ProcessInfo)
        self._check_results_stress(mp.Elements[1],
                                   StructuralMechanicsApplication.SHELL_STRESS_TOP_SURFACE,
                                   shell_stress_top_surface_results,mp.ProcessInfo)
        self._check_results_stress(mp.Elements[1],
                                   StructuralMechanicsApplication.SHELL_STRESS_BOTTOM_SURFACE,
                                   shell_stress_bottom_surface_results,mp.ProcessInfo)
        
        # Check results of doubles on 2nd element @ Gauss Point [0] only
        self.assertAlmostEqual(mp.Elements[1].CalculateOnIntegrationPoints(StructuralMechanicsApplication.VON_MISES_STRESS, 
                               mp.ProcessInfo)[0], shell_von_mises_result, 10)
                    
        if do_post_processing:
            self.__post_process(mp)


    def test_thin_shell_triangle(self):
        element_name = "ShellThinElementCorotational3D3N"
        displacement_results = [0.0002324779832 , -0.0002233435997 , 0.0002567143455]
        rotation_results     = [0.0003627433341 , -0.0001926662603 , -0.0004682681704]
        shell_stress_middle_surface_results = [0.5441631346531 , 0.9137998870586 , 0.0 , -1.8281753448172 , 0.0 , 0.0]
        shell_stress_top_surface_results    = [-0.9165579881878 , -1.8047024313424 , 0.0 , -8.8090405278075 , 0.0 , 0.0]
        shell_stress_bottom_surface_results = [2.0048842574939 , 3.6323022054597 , 0.0 , 5.1526898381731 , 0.0 , 0.0]
        shell_von_mises_result = 7.93786988370381

        self.execute_shell_test(element_name, 
                                displacement_results, 
                                rotation_results,
                                shell_stress_middle_surface_results,
                                shell_stress_top_surface_results,
                                shell_stress_bottom_surface_results,
                                shell_von_mises_result,
                                False) # Do PostProcessing for GiD?


    def test_thick_shell_triangle(self):
        element_name = "ShellThickElementCorotational3D3N"
        displacement_results = [7.18997182e-05 , -0.0001572802804 , 0.0005263940488]
        rotation_results     = [0.0003316612014 , -0.0002798472414 , 5.141506e-07]
        shell_stress_middle_surface_results = [0.3253823265812 , 2.9161070474844 , 0.2817211197418 , -3.1242915352812 , -1.8041533899431 , 0.0]
        shell_stress_top_surface_results    = [-4.1557154888598 , -3.4345294031775 , 0.0 , -9.9078945174138 , 0.0 , 0.0]
        shell_stress_bottom_surface_results = [4.8064801420222 , 9.2667434981463 , 0.0 , 3.6593114468514 , 0.0 , 0.0]
        shell_von_mises_result = 16.628950929004937

        self.execute_shell_test(element_name, 
                                displacement_results, 
                                rotation_results, 
                                shell_stress_middle_surface_results,
                                shell_stress_top_surface_results,
                                shell_stress_bottom_surface_results,
                                shell_von_mises_result,
                                False) # Do PostProcessing for GiD?


    def test_thin_shell_quadrilateral(self):
        element_name = "ShellThinElementCorotational3D4N"
        displacement_results = [0.0021909310921 , -0.0021683746759 , 0.0007191338749]
        rotation_results     = [0.0028191154606 , 0.0008171818407 , -0.0069146010725]
        shell_stress_middle_surface_results = [3.2873048917874 , -11.2253739604803 , 0.0 , 3.3609024020594 , 0.0 , 0.0]
        shell_stress_top_surface_results    = [20.2302593264534 , 4.5024470621022 , 0.0 , -8.6119092854543 , 0.0 , 0.0]
        shell_stress_bottom_surface_results = [-13.6556495428787 , -26.9531949830627 , 0.0 , 15.333714089573 , 0.0 , 0.0]
        shell_von_mises_result = 53.013352444470044

        self.execute_shell_test(element_name, 
                                displacement_results, 
                                rotation_results, 
                                shell_stress_middle_surface_results,
                                shell_stress_top_surface_results,
                                shell_stress_bottom_surface_results,
                                shell_von_mises_result,
                                False) # Do PostProcessing for GiD?


    def test_thick_shell_quadrilateral(self):
        element_name = "ShellThickElementCorotational3D4N"
        displacement_results = [0.0003572969872 , -0.0006341259132 , 0.00127807995]
        rotation_results     = [0.0012082600485 , -0.0004098356773 , -0.001167379835]
        shell_stress_middle_surface_results = [2.888402606167 , -3.9172901967814 , -10.5504072853098 , 2.3194532358326 , 9.1735173861104 , 0.0]
        shell_stress_top_surface_results    = [2.888402606167 , -3.9172901967814 , 0.0 , 2.3194532358326 , 0.0 , 0.0]
        shell_stress_bottom_surface_results = [2.888402606167 , -3.9172901967814 , 0.0 , 2.3194532358326 , 0.0 , 0.0]
        shell_von_mises_result = 25.2873931232983

        self.execute_shell_test(element_name, 
                                displacement_results, 
                                rotation_results, 
                                shell_stress_middle_surface_results,
                                shell_stress_top_surface_results,
                                shell_stress_bottom_surface_results,
                                shell_von_mises_result,
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
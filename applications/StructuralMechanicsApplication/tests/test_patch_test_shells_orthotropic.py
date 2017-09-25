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
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.VOLUME_ACCELERATION)  
        mp.AddNodalSolutionStepVariable(StructuralMechanicsApplication.POINT_LOAD)      
        
    
    def _add_dofs(self,mp):
        for node in mp.Nodes:
            node.AddDof(KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.REACTION_X)
            node.AddDof(KratosMultiphysics.ROTATION_X)
            node.AddDof(KratosMultiphysics.DISPLACEMENT_Y, KratosMultiphysics.REACTION_Y)
            node.AddDof(KratosMultiphysics.ROTATION_Y)
            node.AddDof(KratosMultiphysics.DISPLACEMENT_Z, KratosMultiphysics.REACTION_Z)
            node.AddDof(KratosMultiphysics.ROTATION_Z)


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
        for node in mp.Nodes:
            node.Fix(KratosMultiphysics.DISPLACEMENT_X)
            node.Fix(KratosMultiphysics.DISPLACEMENT_Y)
            node.Fix(KratosMultiphysics.DISPLACEMENT_Z)
            # Adding rotations does not work for some reason...
            # node.Fix(KratosMultiphysics.ROTATION_X)
            # node.Fix(KratosMultiphysics.ROTATION_Y)
            # node.Fix(KratosMultiphysics.ROTATION_Z)


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
        self.assertAlmostEqual(stress[0,0], reference_stress_results[0], 10)
        self.assertAlmostEqual(stress[0,1], reference_stress_results[1], 10)
        self.assertAlmostEqual(stress[0,2], reference_stress_results[2], 10)
        self.assertAlmostEqual(stress[1,1], reference_stress_results[3], 10)
        self.assertAlmostEqual(stress[1,2], reference_stress_results[4], 10)
        self.assertAlmostEqual(stress[2,2], reference_stress_results[5], 10)


    def execute_shell_test(self, element_name, displacement_results, rotation_results, shell_stress_top_surface_results, shell_stress_bottom_surface_results, tsai_wu_result,do_post_processing):
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
        
        
#        print("\n\n---------------------------------------------------------\n",element_name,"\tDisplacement results:\n")
#        stress_out = mp.Nodes[3].GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT)
#        print("displacement_results = ","[",stress_out[0],", ",stress_out[1],", ",stress_out[2],"]")
#        stress_out = mp.Nodes[3].GetSolutionStepValue(KratosMultiphysics.ROTATION)
#        print("rotation_results = ","[",stress_out[0],", ",stress_out[1],", ",stress_out[2],"]")
#
#        
#        
#        print("\tStress results:\n")
#        
#        #results of the second element taken
#        stress_out = mp.Elements[1].CalculateOnIntegrationPoints(StructuralMechanicsApplication.SHELL_ORTHOTROPIC_STRESS_BOTTOM_SURFACE, mp.ProcessInfo)[0]
#
#        print("Element 2 results\n","SHELL_ORTHOTROPIC_STRESS_BOTTOM_SURFACE = ","[",stress_out[0,0],", ",stress_out[0,1],", ",stress_out[0,2],", ",stress_out[1,1],", ",stress_out[1,2],", ",stress_out[2,2],"]")
#        
#        stress_out = mp.Elements[1].CalculateOnIntegrationPoints(StructuralMechanicsApplication.SHELL_ORTHOTROPIC_STRESS_TOP_SURFACE, mp.ProcessInfo)[0]
#
#        print("Element 2 results\n","SHELL_ORTHOTROPIC_STRESS_TOP_SURFACE = ","[",stress_out[0,0],", ",stress_out[0,1],", ",stress_out[0,2],", ",stress_out[1,1],", ",stress_out[1,2],", ",stress_out[2,2],"]")
#        
#        #results of the second element taken
#        stress_out = mp.Elements[1].CalculateOnIntegrationPoints(StructuralMechanicsApplication.TSAI_WU_RESERVE_FACTOR, mp.ProcessInfo)[0]
#        
#        print("TSAI_WU_RESERVE_FACTOR = ",stress_out)
        
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
                               mp.ProcessInfo)[0], tsai_wu_result, 10)
                    
        if do_post_processing:
            self.__post_process(mp)


    def test_thin_shell_triangle(self):
        element_name = "ShellThinElementCorotational3D3N"
        displacement_results =  [ 0.0030584358098890915 ,  -0.0023435757330126693 ,  0.0031185473927106983 ]
        rotation_results =  [ 0.004121996123937352 ,  -0.0016014578481286985 ,  -0.005659631319086533 ]
        shell_stress_top_surface_results =  [ 0.8677118027303852 ,  5.239970870392824 ,  0.0 ,  5.458479391807029 ,  0.0 ,  0.0 ]
        shell_stress_bottom_surface_results =  [ 2.4424603977815895 ,  -3.941407794942941 ,  0.0 ,  -8.030484439805466 ,  0.0 ,  0.0 ]
        tsai_wu_result =  9.673893059720848

        self.execute_shell_test(element_name, 
                                displacement_results, 
                                rotation_results,
                                shell_stress_top_surface_results,
                                shell_stress_bottom_surface_results,
                                tsai_wu_result,
                                False) # Do PostProcessing for GiD?


    def test_thick_shell_triangle(self):
        element_name = "ShellThickElementCorotational3D3N"
        displacement_results =  [ -0.0020376091655285294 ,  -0.002111797305261568 ,  0.02228716026170336 ]
        rotation_results =  [ 0.012334859280482497 ,  -0.010873020793242233 ,  -0.0001569679156864581 ]
        shell_stress_top_surface_results =  [ 4.10385008882551 ,  7.665054447303106 ,  3.944609072109664 ,  -2.0590647897139545 ,  1.411446166379739 ,  0.0 ]
        shell_stress_bottom_surface_results =  [ 3.2912405481215083 ,  0.7649741019715389 ,  3.944609072109664 ,  -1.3979528569250694 ,  1.411446166379739 ,  0.0 ]
        tsai_wu_result =  7.546175761508479

        self.execute_shell_test(element_name, 
                                displacement_results, 
                                rotation_results,
                                shell_stress_top_surface_results,
                                shell_stress_bottom_surface_results,
                                tsai_wu_result,
                                False) # Do PostProcessing for GiD?


    def test_thin_shell_quadrilateral(self):
        element_name = "ShellThinElementCorotational3D4N"
        displacement_results =  [ 0.028324544197644595 ,  -0.028004212791980197 ,  0.005819183507034744 ]
        rotation_results =  [ 0.02562623966199286 ,  0.010203033530205864 ,  -0.07750730506722403 ]
        shell_stress_top_surface_results =  [ 1.7825139429109182 ,  -5.816144708664789 ,  0.0 ,  6.656606883510429 ,  0.0 ,  0.0 ]
        shell_stress_bottom_surface_results =  [ 17.91814192958843 ,  -9.248506597807062 ,  0.0 ,  -3.0956520381202655 ,  0.0 ,  0.0 ]  
        tsai_wu_result =  4.860951226584079

        self.execute_shell_test(element_name, 
                                displacement_results, 
                                rotation_results, 
                                shell_stress_top_surface_results,
                                shell_stress_bottom_surface_results,
                                tsai_wu_result,
                                False) # Do PostProcessing for GiD?


    def test_thick_shell_quadrilateral(self):
        element_name = "ShellThickElementCorotational3D4N"
        displacement_results =  [ 0.0010921056093753834 ,  -0.010786255165378924 ,  0.0377492219218906 ]
        rotation_results =  [ 0.02486239006907997 ,  -0.012527260276340912 ,  -0.01693998229995565 ]
        shell_stress_top_surface_results =  [ 6.271587609566849 ,  -1.12820161459703 ,  -2.305720987313501 ,  4.863927115578198 ,  2.9325229909436965 ,  0.0 ]
        shell_stress_bottom_surface_results =  [ 6.271587609566865 ,  -1.1282016145970397 ,  -2.305720987313501 ,  4.863927115578207 ,  2.9325229909436965 ,  0.0 ]
        tsai_wu_result =  7.047288919615517

        self.execute_shell_test(element_name, 
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
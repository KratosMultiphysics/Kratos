from __future__ import print_function, absolute_import, division
import KratosMultiphysics
from KratosMultiphysics import * 

import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import KratosMultiphysics.KratosUnittest as KratosUnittest

orthotropic_props = []

class TestPatchTestShells(KratosUnittest.TestCase):
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
#        mp.GetProperties()[1].SetValue(KratosMultiphysics.YOUNG_MODULUS,100e3)
#        mp.GetProperties()[1].SetValue(KratosMultiphysics.POISSON_RATIO,0.3)
#        mp.GetProperties()[1].SetValue(KratosMultiphysics.THICKNESS,1.0)
#        mp.GetProperties()[1].SetValue(KratosMultiphysics.DENSITY,1.0)

        0.005,0,7850,20010000,1301000,0.3,1001000,1001000,1001000,80000,50000,4000,30000,6000,6000,6000
        
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
        
        
        print("\n\n---------------------------------------------------------\n",element_name,"\tDisplacement results:\n")
        stress_out = mp.Nodes[3].GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT)
        print("displacement_results = ","[",stress_out[0],", ",stress_out[1],", ",stress_out[2],"]")
        stress_out = mp.Nodes[3].GetSolutionStepValue(KratosMultiphysics.ROTATION)
        print("rotation_results = ","[",stress_out[0],", ",stress_out[1],", ",stress_out[2],"]")

        
        
        print("\tStress results:\n")
        
        #results of the second element taken
        stress_out = mp.Elements[1].CalculateOnIntegrationPoints(StructuralMechanicsApplication.SHELL_ORTHOTROPIC_STRESS_BOTTOM_SURFACE, mp.ProcessInfo)[0]

        print("Element 2 results\n","SHELL_ORTHOTROPIC_STRESS_BOTTOM_SURFACE = ","[",stress_out[0,0],", ",stress_out[0,1],", ",stress_out[0,2],", ",stress_out[1,1],", ",stress_out[1,2],", ",stress_out[2,2],"]")
        
        #results of the second element taken
        stress_out = mp.Elements[1].CalculateOnIntegrationPoints(StructuralMechanicsApplication.TSAI_WU_RESERVE_FACTOR, mp.ProcessInfo)[0]
        
        print("TSAI_WU_RESERVE_FACTOR = ",stress_out)
        
        # Check displacements
        self._check_results(mp.Nodes[3],displacement_results, rotation_results)
        
        # Check stresses at each surface
#        self._check_results_stress(mp.Elements[1],
#                                   StructuralMechanicsApplication.SHELL_STRESS_MIDDLE_SURFACE,
#                                   shell_stress_middle_surface_results,mp.ProcessInfo)
#        self._check_results_stress(mp.Elements[1],
#                                   StructuralMechanicsApplication.SHELL_STRESS_TOP_SURFACE,
#                                   shell_stress_top_surface_results,mp.ProcessInfo)
#        self._check_results_stress(mp.Elements[1],
#                                   StructuralMechanicsApplication.SHELL_STRESS_BOTTOM_SURFACE,
#                                   shell_stress_bottom_surface_results,mp.ProcessInfo)
#        
#        # Check results of doubles on 2nd element @ Gauss Point [0] only
#        self.assertAlmostEqual(mp.Elements[1].CalculateOnIntegrationPoints(StructuralMechanicsApplication.VON_MISES_STRESS, 
#                               mp.ProcessInfo)[0], shell_von_mises_result, 10)
                    
        if do_post_processing:
            self.__post_process(mp)


    def test_thin_shell_triangle(self):
        element_name = "ShellThinElementCorotational3D3N"
        displacement_results = [0.0001967925754 , -0.0002074508275 , 0.0007102373246]
        rotation_results     = [0.0007200850431 , -0.0005274945235 , -0.0004217630272]
        shell_stress_middle_surface_results = [0.25899539638307356  ,  2.707180029120507  ,  0.0  ,  -11.81966960463972  ,  0.0  ,  0.0]
        shell_stress_top_surface_results =  [ -12.07383381852614 ,  -20.247177695843806 ,  0.0 ,  -70.76152753587681 ,  0.0 ,  0.0 ]
        shell_stress_bottom_surface_results =  [ 12.591824611292282 ,  25.661537754084822 ,  0.0 ,  47.122188326597374 ,  0.0 ,  0.0 ]
        shell_von_mises_result =  64.7311391504582

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
        displacement_results = [-8.37692736e-05 , -0.0001866547032 , 0.0013194833292]
        rotation_results     = [0.0009557422766 , -0.0008172919756 , -0.0001582810935]
        shell_stress_middle_surface_results = [1.283076815157277  ,  7.2488387046371  ,  12.850376963002141  ,  -6.060997105375106  ,  4.164269658572349  ,  0.0]
        shell_stress_top_surface_results =  [ -0.4308457214682144 ,  -7.655312031416704 ,  0.0 ,  -2.9286449474358367 ,  0.0 ,  0.0 ]
        shell_stress_bottom_surface_results =  [ 2.996999351782768 ,  22.152989440690906 ,  0.0 ,  -9.193349263314376 ,  0.0 ,  0.0 ]
        shell_von_mises_result =  39.91637459806218

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
        displacement_results = [0.0025324078566 , -0.0025556964999 , 0.0010347939593]
        rotation_results     = [0.0029227371131 , 0.0005461189484 , -0.0073009476025]
        shell_stress_middle_surface_results = [5.115629891110861  ,  -17.33491344410078  ,  0.0  ,  2.748284695572128  ,  0.0  ,  0.0]
        shell_stress_top_surface_results =  [ 30.069565917702892 ,  -20.363353385977423 ,  0.0 ,  -19.9770634365775 ,  0.0 ,  0.0 ]
        shell_stress_bottom_surface_results =  [ -19.838306135481172 ,  -14.30647350222414 ,  0.0 ,  25.473632827721755 ,  0.0 ,  0.0 ]
        shell_von_mises_result =  56.106701625055585

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
        displacement_results =  [ 0.00027136365152691065 ,  -0.0006682902441643944 ,  0.0021084790619259235 ]
        rotation_results =  [ 0.0017763515742804457 ,  -0.000981713367534383 ,  -0.001165116922497465 ]
        shell_stress_middle_surface_results = [6.921145522180776  ,  -2.785925304755268  ,  -14.526750876637958  ,  15.381301373773823  ,  10.378077258376223  ,  0.0]
        shell_stress_top_surface_results =  [ 6.921145522180776 ,  -2.785925304755268 ,  0.0 ,  15.381301373773823 ,  0.0 ,  0.0 ]
        shell_stress_bottom_surface_results =  [ 6.921145522180776 ,  -2.785925304755268 ,  0.0 ,  15.381301373773823 ,  0.0 ,  0.0 ]
        shell_von_mises_result =  34.02216244466042

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
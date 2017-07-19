from __future__ import print_function, absolute_import, division
import KratosMultiphysics 

import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import KratosMultiphysics.KratosUnittest as KratosUnittest


class TestQuadraticElements(KratosUnittest.TestCase):
    def setUp(self):
        pass
    
    def _add_variables(self,mp):
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.VOLUME_ACCELERATION)        
        
    
    def _apply_BCs(self,mp):
        for node in mp.Nodes:
            node.Fix(KratosMultiphysics.DISPLACEMENT_X)
            node.Fix(KratosMultiphysics.DISPLACEMENT_Y)
            node.Fix(KratosMultiphysics.DISPLACEMENT_Z)
        
        for node in mp.Nodes:
            u = KratosMultiphysics.Vector(3)
            u[0] = node.X0**2
            u[1] = node.Y0**2
            u[2] = node.Z0**2
            
            node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT,0,u)

    def _apply_material_properties(self,mp,dim):
        #define properties
        for prop in mp.GetProperties():
            prop.SetValue(KratosMultiphysics.YOUNG_MODULUS,210e9)
            prop.SetValue(KratosMultiphysics.POISSON_RATIO,0.3)
            prop.SetValue(KratosMultiphysics.THICKNESS,1.0)
            prop.SetValue(KratosMultiphysics.DENSITY,1.0)
            
            g = [0,0,0]
            prop.SetValue(KratosMultiphysics.VOLUME_ACCELERATION,g)
            
            if(dim == 2):
                cl = StructuralMechanicsApplication.LinearElasticPlaneStress2DLaw()
            else:
                cl = StructuralMechanicsApplication.LinearElastic3DLaw()
            prop.SetValue(KratosMultiphysics.CONSTITUTIVE_LAW,cl) 
        
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
        

        #strategy = KratosMultiphysics.ResidualBasedNewtonRaphsonStrategy(mp, 
                                                                        #scheme, 
                                                                        #linear_solver, 
                                                                        #convergence_criterion, 
                                                                        #builder_and_solver, 
                                                                        #max_iters, 
                                                                        #compute_reactions, 
                                                                        #reform_step_dofs, 
                                                                        #move_mesh_flag)
        strategy.SetEchoLevel(0)
        
        strategy.Check()
        strategy.Solve()
            
    def _check_outputs(self,mp,dim):
        for elem in mp.Elements:
            strains = elem.CalculateOnIntegrationPoints(KratosMultiphysics.GREEN_LAGRANGE_STRAIN_VECTOR, mp.ProcessInfo)
            coords = elem.CalculateOnIntegrationPoints(KratosMultiphysics.INTEGRATION_COORDINATES, mp.ProcessInfo)
            for strain,coord in zip(strains, coords):
                for i in range(2):
                    self.assertAlmostEqual((coord[i] - 0.5 * strain[i])/coord[i], 0.0 , 1)

    def test_Quad8(self):
        dim = 2
        mp = KratosMultiphysics.ModelPart("solid_part")
        self._add_variables(mp)
        
        KratosMultiphysics.ModelPartIO("quadratic_test/static_quadratic_quad_test").ReadModelPart(mp)
        
        self._apply_material_properties(mp,dim)
        
        for node in mp.Nodes:
            node.AddDof(KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.REACTION_X)
            node.AddDof(KratosMultiphysics.DISPLACEMENT_Y, KratosMultiphysics.REACTION_Y)
            node.AddDof(KratosMultiphysics.DISPLACEMENT_Z, KratosMultiphysics.REACTION_Z)
            
        self._apply_BCs(mp)
        self._solve(mp)
        self._check_outputs(mp,dim)
        
        #self.__post_process(mp)
        
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

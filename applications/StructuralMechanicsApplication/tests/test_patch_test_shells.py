from __future__ import print_function, absolute_import, division
import KratosMultiphysics 

import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import KratosMultiphysics.KratosUnittest as KratosUnittest

'''
Questions Andreas:
- Tests for Corotational and other version?
- Fix Rotation?
- 2D ???
- Z-Component of node != 0?
'''

class TestPatchTestShells(KratosUnittest.TestCase):
    def setUp(self):
        pass
    

    def _add_variables(self,mp):
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.ROTATION)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.VOLUME_ACCELERATION) 


    def _add_dofs(self,mp):
        for node in mp.Nodes:
            node.AddDof(KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.REACTION_X)
            node.AddDof(KratosMultiphysics.ROTATION_X)
            node.AddDof(KratosMultiphysics.DISPLACEMENT_Y, KratosMultiphysics.REACTION_Y)
            node.AddDof(KratosMultiphysics.ROTATION_Y)
            node.AddDof(KratosMultiphysics.DISPLACEMENT_Z, KratosMultiphysics.REACTION_Z)
            node.AddDof(KratosMultiphysics.ROTATION_Z)       
        
    
    def _apply_BCs(self,mp,A,b):
        for node in mp.Nodes:
            node.Fix(KratosMultiphysics.DISPLACEMENT_X)
            node.Fix(KratosMultiphysics.DISPLACEMENT_Y)
            node.Fix(KratosMultiphysics.DISPLACEMENT_Z)
        
        for node in mp.Nodes:
            xvec = KratosMultiphysics.Vector(3)
            xvec[0] = node.X0
            xvec[1] = node.Y0
            xvec[2] = node.Z0
            
            u = KratosMultiphysics.Vector()
            u = A*xvec
            u += b
            
            node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT,0,u)


    def _apply_material_properties(self,mp,dim):
        #define properties
        mp.GetProperties()[1].SetValue(KratosMultiphysics.YOUNG_MODULUS,210e9)
        mp.GetProperties()[1].SetValue(KratosMultiphysics.POISSON_RATIO,0.3)
        mp.GetProperties()[1].SetValue(KratosMultiphysics.THICKNESS,1.0)
        mp.GetProperties()[1].SetValue(KratosMultiphysics.DENSITY,1.0)
        
        g = [0,0,0]
        mp.GetProperties()[1].SetValue(KratosMultiphysics.VOLUME_ACCELERATION,g)
        
        if(dim == 2):
            cl = StructuralMechanicsApplication.LinearElasticPlaneStress2DLaw()
        else:
            cl = StructuralMechanicsApplication.LinearElastic3DLaw()
        mp.GetProperties()[1].SetValue(KratosMultiphysics.CONSTITUTIVE_LAW,cl) 
            

    def _define_movement(self,dim):
        if(dim == 2):
            #define the applied motion - the idea is that the displacement is defined as u = A*xnode + b
            #so that the displcement is linear and the exact F = I + A
            A = KratosMultiphysics.Matrix(3,3)
            A[0,0] = 1.0e-10;  A[0,1] = 2.0e-10; A[0,2] = 0.0
            A[1,0] = 0.5e-10;  A[1,1] = 0.7e-10; A[1,2] = 0.0
            A[2,1] = 0.0;  A[2,1] = 0.0; A[2,2] = 0.0
                    
            b = KratosMultiphysics.Vector(3)
            b[0] = 0.5e-10
            b[1] = -0.2e-10     
            b[2] = 0.0
            
        else:
            #define the applied motion - the idea is that the displacement is defined as u = A*xnode + b
            #so that the displcement is linear and the exact F = I + A
            A = KratosMultiphysics.Matrix(3,3)
            A[0,0] = 1.0e-10;   A[0,1] = 2.0e-10; A[0,2] = 0.0
            A[1,0] = 0.5e-10;   A[1,1] = 0.7e-10; A[1,2] = 0.1e-10
            A[2,1] = -0.2e-10;  A[2,1] = 0.0;     A[2,2] = -0.3e-10
                    
            b = KratosMultiphysics.Vector(3)
            b[0] = 0.5e-10
            b[1] = -0.2e-10     
            b[2] = 0.7e-10
        
        return A,b
        

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
        
    
    def _check_results(self,mp,A,b):
        ##check that the results are exact on the nodes
        for node in mp.Nodes:
            xvec = KratosMultiphysics.Vector(len(b))
            xvec[0] = node.X0
            xvec[1] = node.Y0
            xvec[2] = node.Z0
            
            u = KratosMultiphysics.Vector(2)
            u = A*xvec
            u += b            
            
            d = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT)
            self.assertAlmostEqual(d[0], u[0])
            self.assertAlmostEqual(d[1], u[1])
            self.assertAlmostEqual(d[2], u[2])
            

    def _check_outputs(self,mp,A,dim):
        E = mp.GetProperties()[1].GetValue(KratosMultiphysics.YOUNG_MODULUS)
        NU =mp.GetProperties()[1].GetValue(KratosMultiphysics.POISSON_RATIO)
        
        #given the matrix A, the analytic deformation graident is F+I
        F = A
        for i in range(3):
            F[i,i] += 1.0
        
        #here compute the Cauchy green strain tensor
        Etensor = KratosMultiphysics.Matrix(3,3)

        for i in range(3):
            for j in range(3):
                Etensor[i,j] = 0.0
                
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    Etensor[i,j] += A[k,i]*A[k,j]
                    
        for i in range(3):
            Etensor[i,i] -= 1.0
                    
        for i in range(3):
            for j in range(3):
                Etensor[i,j] = 0.5*Etensor[i,j]
        
        if(dim == 2):
            #verify strain
            reference_strain = KratosMultiphysics.Vector(3)
            reference_strain[0] = Etensor[0,0]
            reference_strain[1] = Etensor[1,1]
            reference_strain[2] = 2.0*Etensor[0,1]
        else:
            reference_strain = KratosMultiphysics.Vector(6)
            reference_strain[0] = Etensor[0,0]
            reference_strain[1] = Etensor[1,1]
            reference_strain[2] = Etensor[2,2]
            reference_strain[3] = 2.0*Etensor[0,1]
            reference_strain[4] = 2.0*Etensor[1,2]
            reference_strain[5] = 2.0*Etensor[0,2]
            
        for elem in mp.Elements:
            out = elem.CalculateOnIntegrationPoints(KratosMultiphysics.GREEN_LAGRANGE_STRAIN_VECTOR, mp.ProcessInfo)
            for strain in out:
                for i in range(len(reference_strain)):
                    self.assertAlmostEqual(reference_strain[i], strain[i])
                    
        #finally compute stress
        if(dim == 2):
            #here assume plane stress
            c1 = E / (1.00 - NU*NU);
            c2 = c1 * NU;
            c3 = 0.5* E / (1 + NU);
            reference_stress = KratosMultiphysics.Vector(3)
            reference_stress[0] = c1*reference_strain[0] + c2 * (reference_strain[1])	;
            reference_stress[1] = c1*reference_strain[1] + c2 * (reference_strain[0])	;
            reference_stress[2] = c3*reference_strain[2];
        else:
            c1 = E / (( 1.00 + NU ) * ( 1 - 2 * NU ) );
            c2 = c1 * ( 1 - NU );
            c3 = c1 * NU;
            c4 = c1 * 0.5 * ( 1 - 2 * NU );
            reference_stress = KratosMultiphysics.Vector(6)
            reference_stress[0] = c2*reference_strain[0] + c3 * (reference_strain[1] + reference_strain[2])
            reference_stress[1] = c2*reference_strain[1] + c3 * (reference_strain[0] + reference_strain[2])
            reference_stress[2] = c2*reference_strain[2] + c3 * (reference_strain[0] + reference_strain[1])
            reference_stress[3] = c4*reference_strain[3]
            reference_stress[4] = c4*reference_strain[4]
            reference_stress[5] = c4*reference_strain[5]
            
        for elem in mp.Elements:
            out = elem.CalculateOnIntegrationPoints(KratosMultiphysics.PK2_STRESS_VECTOR, mp.ProcessInfo)
            for stress in out:
                for i in range(len(reference_stress)):
                    self.assertAlmostEqual(reference_stress[i], stress[i],2)        


    def test_thin_shell_triangle(self):
        dim = 3
        mp = KratosMultiphysics.ModelPart("solid_part")
        mp.SetBufferSize(2)
        self._add_variables(mp)
        self._apply_material_properties(mp,dim)
        
        #create nodes
        mp.CreateNewNode(1,0.5,0.5,0.0)
        mp.CreateNewNode(2,0.7,0.2,0.0)
        mp.CreateNewNode(3,0.9,0.8,0.0)
        mp.CreateNewNode(4,0.3,0.7,0.0)
        mp.CreateNewNode(5,0.6,0.6,0.0)
        
        self._add_dofs(mp)
            
        #create a submodelpart for boundary conditions
        bcs = mp.CreateSubModelPart("BoundaryCondtions")
        bcs.AddNodes([1,2,3,4])
                
        #create Element
        mp.CreateNewElement("ShellThinElementCorotational3D3N", 1, [1,2,5], mp.GetProperties()[1])
        mp.CreateNewElement("ShellThinElementCorotational3D3N", 2, [2,3,5], mp.GetProperties()[1])
        mp.CreateNewElement("ShellThinElementCorotational3D3N", 3, [3,4,5], mp.GetProperties()[1])
        mp.CreateNewElement("ShellThinElementCorotational3D3N", 4, [4,1,5], mp.GetProperties()[1])
        
        A,b = self._define_movement(dim)
        
        self._apply_BCs(bcs,A,b)
        self._solve(mp)
        self._check_results(mp,A,b)
        self._check_outputs(mp,A,dim)
        
        #checking consistent mass matrix
        M = KratosMultiphysics.Matrix(0,0)
        mp.Elements[1].CalculateMassMatrix(M,mp.ProcessInfo)
        Area = mp.Elements[1].GetArea()
        for i in range(3):
            for j in range(3):
                for k in range(dim):
                    if(i==j):
                        coeff = Area/6.0
                    else:
                        coeff = Area/12.0
                    self.assertAlmostEqual(M[i*dim+k,j*dim+k],coeff)
                    
        self.__post_process(mp)
                    

    def test_thick_shell_triangle(self):
        dim = 3
        mp = KratosMultiphysics.ModelPart("solid_part")
        mp.SetBufferSize(2)
        self._add_variables(mp)
        self._apply_material_properties(mp,dim)
        
        #create nodes
        mp.CreateNewNode(1,0.5,0.5,0.0)
        mp.CreateNewNode(2,0.7,0.2,0.0)
        mp.CreateNewNode(3,0.9,0.8,0.0)
        mp.CreateNewNode(4,0.3,0.7,0.0)
        mp.CreateNewNode(5,0.6,0.6,0.0)
        
        self._add_dofs(mp)
            
        #create a submodelpart for boundary conditions
        bcs = mp.CreateSubModelPart("BoundaryCondtions")
        bcs.AddNodes([1,2,3,4])
                
        #create Element
        mp.CreateNewElement("ShellThickElementCorotational3D3N", 1, [1,2,5], mp.GetProperties()[1])
        mp.CreateNewElement("ShellThickElementCorotational3D3N", 2, [2,3,5], mp.GetProperties()[1])
        mp.CreateNewElement("ShellThickElementCorotational3D3N", 3, [3,4,5], mp.GetProperties()[1])
        mp.CreateNewElement("ShellThickElementCorotational3D3N", 4, [4,1,5], mp.GetProperties()[1])
        
        A,b = self._define_movement(dim)
        
        self._apply_BCs(bcs,A,b)
        self._solve(mp)
        self._check_results(mp,A,b)
        self._check_outputs(mp,A,dim)
        
        #checking consistent mass matrix
        M = KratosMultiphysics.Matrix(0,0)
        mp.Elements[1].CalculateMassMatrix(M,mp.ProcessInfo)
        Area = mp.Elements[1].GetArea()
        for i in range(3):
            for j in range(3):
                for k in range(dim):
                    if(i==j):
                        coeff = Area/6.0
                    else:
                        coeff = Area/12.0
                    self.assertAlmostEqual(M[i*dim+k,j*dim+k],coeff)
                    
        #self.__post_process(mp)
        

    def test_thin_shell_quadrilateral(self): 
        dim = 3
        mp = KratosMultiphysics.ModelPart("solid_part")
        mp.SetBufferSize(2)
        self._add_variables(mp)
        self._apply_material_properties(mp,dim)
        
        #create nodes
        mp.CreateNewNode(1,0.00,3.00,0.00)
        mp.CreateNewNode(2,1.00,2.25,0.00)
        mp.CreateNewNode(3,0.75,1.00,0.00)
        mp.CreateNewNode(4,2.25,2.00,0.00)
        mp.CreateNewNode(5,0.00,0.00,0.00)
        mp.CreateNewNode(6,3.00,3.00,0.00)
        mp.CreateNewNode(7,2.00,0.75,0.00)
        mp.CreateNewNode(8,3.00,0.00,0.00)
        
        self._add_dofs(mp)
            
        #create a submodelpart for boundary conditions
        bcs = mp.CreateSubModelPart("BoundaryCondtions")
        bcs.AddNodes([1,5,6,8])
                
        #create Element
        mp.CreateNewElement("ShellThinElementCorotational3D4N", 1, [8,7,3,5], mp.GetProperties()[1])
        mp.CreateNewElement("ShellThinElementCorotational3D4N", 2, [6,4,7,8], mp.GetProperties()[1])
        mp.CreateNewElement("ShellThinElementCorotational3D4N", 3, [1,2,4,6], mp.GetProperties()[1])
        mp.CreateNewElement("ShellThinElementCorotational3D4N", 4, [4,2,3,7], mp.GetProperties()[1])
        mp.CreateNewElement("ShellThinElementCorotational3D4N", 5, [2,1,5,3], mp.GetProperties()[1])
        
        A,b = self._define_movement(dim)
        
        self._apply_BCs(bcs,A,b)
        self._solve(mp)
        self._check_results(mp,A,b)
        self._check_outputs(mp,A,dim)
        
        #self.__post_process(mp)
        

    def test_thick_shell_quadrilateral(self): 
        dim = 3
        mp = KratosMultiphysics.ModelPart("solid_part")
        mp.SetBufferSize(2)
        self._add_variables(mp)
        self._apply_material_properties(mp,dim)
        
        #create nodes
        mp.CreateNewNode(1,0.00,3.00,0.00)
        mp.CreateNewNode(2,1.00,2.25,0.00)
        mp.CreateNewNode(3,0.75,1.00,0.00)
        mp.CreateNewNode(4,2.25,2.00,0.00)
        mp.CreateNewNode(5,0.00,0.00,0.00)
        mp.CreateNewNode(6,3.00,3.00,0.00)
        mp.CreateNewNode(7,2.00,0.75,0.00)
        mp.CreateNewNode(8,3.00,0.00,0.00)
        
        self._add_dofs(mp)
            
        #create a submodelpart for boundary conditions
        bcs = mp.CreateSubModelPart("BoundaryCondtions")
        bcs.AddNodes([1,5,6,8])
                
        #create Element
        mp.CreateNewElement("ShellThickElementCorotational3D4N", 1, [8,7,3,5], mp.GetProperties()[1])
        mp.CreateNewElement("ShellThickElementCorotational3D4N", 2, [6,4,7,8], mp.GetProperties()[1])
        mp.CreateNewElement("ShellThickElementCorotational3D4N", 3, [1,2,4,6], mp.GetProperties()[1])
        mp.CreateNewElement("ShellThickElementCorotational3D4N", 4, [4,2,3,7], mp.GetProperties()[1])
        mp.CreateNewElement("ShellThickElementCorotational3D4N", 5, [2,1,5,3], mp.GetProperties()[1])
        
        A,b = self._define_movement(dim)
        
        self._apply_BCs(bcs,A,b)
        self._solve(mp)
        self._check_results(mp,A,b)
        self._check_outputs(mp,A,dim)
        
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

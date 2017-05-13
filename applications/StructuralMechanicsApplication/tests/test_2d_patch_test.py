from __future__ import print_function, absolute_import, division
import KratosMultiphysics 

import KratosMultiphysics.StructuralMechanicsApplication
import KratosMultiphysics.KratosUnittest as KratosUnittest


class TestGeometry(KratosUnittest.TestCase):
    def setUp(self):
        pass
    
    def _add_variables(self,mp):
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.VOLUME_ACCELERATION)        
        
    
    def _apply_BCs(self,mp,A,b):
        for node in mp.Nodes:
            node.Fix(KratosMultiphysics.DISPLACEMENT_X)
            node.Fix(KratosMultiphysics.DISPLACEMENT_Y)
        
        for node in mp.Nodes:
            xvec = KratosMultiphysics.Vector(3)
            xvec[0] = node.X0
            xvec[1] = node.Y0
            xvec[2] = node.Z0
            
            u = KratosMultiphysics.Vector()
            u = A*xvec
            u += b
            
            node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT,0,u)
            
    def _solve(self,mp):
        
        #define a minimal newton raphson solver
        linear_solver = KratosMultiphysics.SkylineLUFactorizationSolver()
        builder_and_solver = KratosMultiphysics.ResidualBasedBlockBuilderAndSolver(linear_solver)
        scheme = KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme()
        convergence_criterion = KratosMultiphysics.ResidualCriteria(1e-14,1e-20)
        
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
        strategy.SetEchoLevel(0)
        
        strategy.Check()
        strategy.Solve()
        
    
    def _check_results(self,mp,A,b):
        
        ##check that the results are exact on the nodes
        for node in mp.Nodes:
            xvec = KratosMultiphysics.Vector(len(b))
            xvec[0] = node.X0
            xvec[1] = node.Y0
            
            u = KratosMultiphysics.Vector(2)
            u = A*xvec
            u += b            
            
            d = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT)
            
            self.assertEqual(d[0], u[0])
            self.assertEqual(d[1], u[1])
            self.assertEqual(d[2], u[2])
            
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
            reference_strain = KratosMultiphysics.Vector(4)
            reference_strain[0] = Etensor[0,0]
            reference_strain[1] = Etensor[1,1]
            reference_strain[2] = Etensor[2,2]
            reference_strain[3] = 2.0*Etensor[0,1]
            reference_strain[4] = 2.0*Etensor[1,2]
            reference_strain[5] = 2.0*Etensor[0,2]
            
        for elem in mp.Elements:
            out = elem.CalculateOnIntegrationPoints(KratosMultiphysics.GREEN_LAGRANGE_STRAIN_TENSOR, mp.ProcessInfo)
            for strain in out:
                for i in range(len(reference_strain)):
                    self.assertAlmostEqual(reference_strain[i], strain[0,i])
                    
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
            raise Exception("3D case not yet implemented")
        
        for elem in mp.Elements:
            out = elem.CalculateOnIntegrationPoints(KratosMultiphysics.PK2_STRESS_TENSOR, mp.ProcessInfo)
            for stress in out:
                for i in range(len(reference_stress)):
                    self.assertAlmostEqual(reference_stress[i], stress[0,i],2)        
        
        

    def test_TL_2D_triangle(self):
        mp = KratosMultiphysics.ModelPart("solid_part")
        self._add_variables(mp)
        
        #create nodes
        mp.CreateNewNode(1,0.5,0.5,0.0)
        mp.CreateNewNode(2,0.7,0.2,0.0)
        mp.CreateNewNode(3,0.9,0.8,0.0)
        mp.CreateNewNode(4,0.3,0.7,0.0)
        mp.CreateNewNode(5,0.6,0.6,0.0)
        
        for node in mp.Nodes:
            node.AddDof(KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.REACTION_X)
            node.AddDof(KratosMultiphysics.DISPLACEMENT_Y, KratosMultiphysics.REACTION_Y)
            node.AddDof(KratosMultiphysics.DISPLACEMENT_Z, KratosMultiphysics.REACTION_Z)
            
        #create a submodelpart for boundary conditions
        bcs = mp.CreateSubModelPart("BoundaryCondtions")
        bcs.AddNodes([1,2,3,4])
        
        #define properties
        mp.GetProperties()[1].SetValue(KratosMultiphysics.YOUNG_MODULUS,210e9)
        mp.GetProperties()[1].SetValue(KratosMultiphysics.POISSON_RATIO,0.3)
        mp.GetProperties()[1].SetValue(KratosMultiphysics.THICKNESS,1.0)
        
        cl = KratosMultiphysics.StructuralMechanicsApplication.LinearPlaneStress()
        mp.GetProperties()[1].SetValue(KratosMultiphysics.CONSTITUTIVE_LAW,cl) #TODO: take care, using solid mechanics here
        
        #create Element
        mp.CreateNewElement("TotalLagrangian2D3N", 1, [1,2,5], mp.GetProperties()[1])
        mp.CreateNewElement("TotalLagrangian2D3N", 2, [2,3,5], mp.GetProperties()[1])
        mp.CreateNewElement("TotalLagrangian2D3N", 3, [3,4,5], mp.GetProperties()[1])
        mp.CreateNewElement("TotalLagrangian2D3N", 4, [4,1,5], mp.GetProperties()[1])
        
        #for elem in mp.Elements:
            #elem.Check()
        
        #define the applied motion
        A = KratosMultiphysics.Matrix(3,3)
        A[0,0] = 1.0;  A[0,1] = 2.0; A[0,2] = 0.0
        A[1,0] = 0.5;  A[1,1] = 0.7; A[1,2] = 0.0
        A[2,1] = 0.0;  A[2,1] = 0.0; A[2,2] = 0.0
                
        b = KratosMultiphysics.Vector(3)
        b[0] = 0.5
        b[1] = -0.2        
        b[2] = 0.0
        
        self._apply_BCs(bcs,A,b)
        self._solve(mp)
        self._check_results(mp,A,b)
        self._check_outputs(mp,A,2)


if __name__ == '__main__':
    KratosUnittest.main()

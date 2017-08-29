from __future__ import print_function, absolute_import, division
import KratosMultiphysics 

import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import KratosMultiphysics.KratosUnittest as KratosUnittest
#import KratosMultiphysics.ExternalSolversApplication as ExternalSolversApplication


class TestTruss3D2N(KratosUnittest.TestCase):
    def setUp(self):
        pass

    def _add_dofs(self,mp):
        for node in mp.Nodes:
            node.AddDof(KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.REACTION_X)
            node.AddDof(KratosMultiphysics.DISPLACEMENT_Y, KratosMultiphysics.REACTION_Y)
            node.AddDof(KratosMultiphysics.DISPLACEMENT_Z, KratosMultiphysics.REACTION_Z)
            
    
    def _add_variables(self,mp):
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)
        mp.AddNodalSolutionStepVariable(StructuralMechanicsApplication.POINT_LOAD)       
        
    def _apply_material_properties(self,mp,dim):
        #define properties
        mp.GetProperties()[1].SetValue(KratosMultiphysics.YOUNG_MODULUS,210e9)
        mp.GetProperties()[1].SetValue(KratosMultiphysics.DENSITY,7850)
        mp.GetProperties()[1].SetValue(StructuralMechanicsApplication.CROSS_AREA,0.01)
        mp.GetProperties()[1].SetValue(StructuralMechanicsApplication.TRUSS_PRESTRESS_PK2,0)

        cl = StructuralMechanicsApplication.LinearElastic3DLaw()
        mp.GetProperties()[1].SetValue(KratosMultiphysics.CONSTITUTIVE_LAW,cl)

        
    def _apply_BCs(self,mp,which_dof):
        
        if (which_dof == 'xyz'):
            for node in mp.Nodes:
                node.Fix(KratosMultiphysics.DISPLACEMENT_X)
                node.Fix(KratosMultiphysics.DISPLACEMENT_Y)
                node.Fix(KratosMultiphysics.DISPLACEMENT_Z)
        if (which_dof == 'xz'):
            for node in mp.Nodes:
                node.Fix(KratosMultiphysics.DISPLACEMENT_X)
                node.Fix(KratosMultiphysics.DISPLACEMENT_Z)

    def _apply_Neumann_BCs(self,mp,which_dof,load_size_dir):

        if(which_dof == 'y'):
            for node in mp.Nodes:
                node.SetSolutionStepValue(StructuralMechanicsApplication.
                POINT_LOAD_Y,0,load_size_dir)        


        
    def _solve(self,mp):
        

        print("load on node 2: ", mp.Nodes[2].GetSolutionStepValue(
            StructuralMechanicsApplication.POINT_LOAD,0))
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

        
    
    def _check_results(self,mp):
        
        displacement_node_2 = mp.Nodes[2].GetSolutionStepValue(
            KratosMultiphysics.DISPLACEMENT)
        self.assertAlmostEqual(displacement_node_2[2], -0.0026619856874997507)

        
    def test_truss3D2N_linear(self):
        dim = 3
        mp = KratosMultiphysics.ModelPart("solid_part")
        self._add_variables(mp)
        self._apply_material_properties(mp,dim)
        
        #create nodes
        mp.CreateNewNode(1,0.0,0.0,0.0)
        mp.CreateNewNode(2,2.0,1.0,0.0)
        mp.CreateNewNode(3,4.0,0.0,0.0)

        self._add_dofs(mp)

        #create condition
        mp.CreateNewCondition("PointLoadCondition3D1N",1,[2],mp.GetProperties()[1])

        #create submodelparts for dirichlet boundary conditions

        bcs_xyz = mp.CreateSubModelPart("Dirichlet_XYZ")
        bcs_xyz.AddNodes([1,3])
        bcs_xz = mp.CreateSubModelPart("Dirichlet_XZ")
        bcs_xz.AddNodes([2])

        #create a submodalpart for neumann boundary conditions
        bcs_neumann = mp.CreateSubModelPart("PointLoad3D_neumann")
        bcs_neumann.AddNodes([2])
        bcs_neumann.AddConditions([1])


        #create Element
        mp.CreateNewElement("TrussLinearElement3D2N", 1, [1,2], mp.GetProperties()[1])
        mp.CreateNewElement("TrussLinearElement3D2N", 2, [2,3], mp.GetProperties()[1])

        #apply boundary conditions
        Force_Y = -1000000.00
        self._apply_BCs(bcs_xyz,'xyz')
        self._apply_BCs(bcs_xz,'xz')
        self._apply_Neumann_BCs(bcs_neumann,'y',Force_Y)

        #solve + compare
        self._solve(mp)
        self._check_results(mp)
        
        
if __name__ == '__main__':
    KratosUnittest.main()


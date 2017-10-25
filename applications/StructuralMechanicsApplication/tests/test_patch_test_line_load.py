from __future__ import print_function, absolute_import, division
import KratosMultiphysics 

import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import KratosMultiphysics.KratosUnittest as KratosUnittest


from math import sqrt, sin, cos, pi, exp, atan

class TestLineLoad(KratosUnittest.TestCase):
    def setUp(self):
        pass

    def _add_dofs(self,mp):
        # Adding dofs AND their corresponding reactions
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.REACTION_X,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Y, KratosMultiphysics.REACTION_Y,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Z, KratosMultiphysics.REACTION_Z,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ROTATION_X, KratosMultiphysics.TORQUE_X,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ROTATION_Y, KratosMultiphysics.TORQUE_Y,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ROTATION_Z, KratosMultiphysics.TORQUE_Z,mp)
    
    def _add_variables(self,mp):
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.ROTATION)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.TORQUE)
        mp.AddNodalSolutionStepVariable(StructuralMechanicsApplication.POINT_LOAD)  
        mp.AddNodalSolutionStepVariable(StructuralMechanicsApplication.LINE_LOAD)  
        mp.AddNodalSolutionStepVariable(StructuralMechanicsApplication.POINT_MOMENT)  
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.VOLUME_ACCELERATION) 
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.ANGULAR_VELOCITY)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.ACCELERATION)     
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.ANGULAR_ACCELERATION)  
        
    def _apply_material_properties_beam(self,mp,dim):
        #define properties
        mp.GetProperties()[0].SetValue(KratosMultiphysics.YOUNG_MODULUS,210e9)
        mp.GetProperties()[0].SetValue(KratosMultiphysics.DENSITY,7850)
        mp.GetProperties()[0].SetValue(StructuralMechanicsApplication.CROSS_AREA,0.01)
        mp.GetProperties()[0].SetValue(KratosMultiphysics.POISSON_RATIO,0.30)

        local_inertia_vector = KratosMultiphysics.Vector(3)
        local_inertia_vector[0] = 0.00001
        local_inertia_vector[1] = 0.00001
        local_inertia_vector[2] = 0.00001

        mp.GetProperties()[0].SetValue(StructuralMechanicsApplication.LOCAL_INERTIA_VECTOR,local_inertia_vector)

        g = [0,0,0]  
        mp.GetProperties()[0].SetValue(KratosMultiphysics.VOLUME_ACCELERATION,g)

        cl = StructuralMechanicsApplication.LinearElastic3DLaw()
        mp.GetProperties()[0].SetValue(KratosMultiphysics.CONSTITUTIVE_LAW,cl)

    def _apply_material_properties_truss(self,mp,dim):
        #define properties
        mp.GetProperties()[0].SetValue(KratosMultiphysics.YOUNG_MODULUS,210e9)
        mp.GetProperties()[0].SetValue(KratosMultiphysics.DENSITY,7850)
        mp.GetProperties()[0].SetValue(StructuralMechanicsApplication.CROSS_AREA,0.01)


        g = [0,0,0]  
        mp.GetProperties()[0].SetValue(KratosMultiphysics.VOLUME_ACCELERATION,g)

        cl = StructuralMechanicsApplication.LinearElastic3DLaw()
        mp.GetProperties()[0].SetValue(KratosMultiphysics.CONSTITUTIVE_LAW,cl)
        
    def _apply_BCs(self,mp,which_dof):
        
        if (which_dof == 'xyz'):
            KratosMultiphysics.VariableUtils().ApplyFixity(KratosMultiphysics.DISPLACEMENT_X, True, mp.Nodes)
            KratosMultiphysics.VariableUtils().ApplyFixity(KratosMultiphysics.DISPLACEMENT_Y, True, mp.Nodes)
            KratosMultiphysics.VariableUtils().ApplyFixity(KratosMultiphysics.DISPLACEMENT_Z, True, mp.Nodes)

        if (which_dof == 'xz'):            
            KratosMultiphysics.VariableUtils().ApplyFixity(KratosMultiphysics.DISPLACEMENT_X, True, mp.Nodes)
            KratosMultiphysics.VariableUtils().ApplyFixity(KratosMultiphysics.DISPLACEMENT_Z, True, mp.Nodes)

        if (which_dof == 'yz'):
            KratosMultiphysics.VariableUtils().ApplyFixity(KratosMultiphysics.DISPLACEMENT_Y, True, mp.Nodes)
            KratosMultiphysics.VariableUtils().ApplyFixity(KratosMultiphysics.DISPLACEMENT_Z, True, mp.Nodes)            

        if (which_dof == 'rotXYZ'):
            KratosMultiphysics.VariableUtils().ApplyFixity(KratosMultiphysics.ROTATION_X, True, mp.Nodes)
            KratosMultiphysics.VariableUtils().ApplyFixity(KratosMultiphysics.ROTATION_Y, True, mp.Nodes)
            KratosMultiphysics.VariableUtils().ApplyFixity(KratosMultiphysics.ROTATION_Z, True, mp.Nodes)
              

    def _apply_Neumann_BCs(self,mp,which_dof,load_size_dir):

        if(which_dof == 'y'):
            KratosMultiphysics.VariableUtils().SetScalarVar(StructuralMechanicsApplication.
                LINE_LOAD_Y, load_size_dir, mp.Nodes)
       

        
    def _solve_linear(self,mp):
        
        linear_solver = KratosMultiphysics.SkylineLUFactorizationSolver()
        builder_and_solver = KratosMultiphysics.ResidualBasedBlockBuilderAndSolver(linear_solver)
        scheme = KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme()
        
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
        return strategy
   
           
    
    def _check_results_linear_beam(self,mp,endNode):
        
        #check displacement result
        displacement_cantilever_tip = mp.Nodes[endNode].GetSolutionStepValue(
            KratosMultiphysics.DISPLACEMENT)

        disp_y_analytical = -10000.00*(2.00**4)/(24*210e9*0.00001)*3.00

        self.assertAlmostEqual(0.00, displacement_cantilever_tip[0])
        self.assertAlmostEqual(disp_y_analytical, displacement_cantilever_tip[1])
        self.assertAlmostEqual(0.00, displacement_cantilever_tip[2])

    def _check_results_linear_truss(self,mp):
        
        #check displacement result
        displacement_truss_tip = mp.Nodes[2].GetSolutionStepValue(
            KratosMultiphysics.DISPLACEMENT)

        disp_y_analytical = -10000.00/(210e9*0.01)*2

        self.assertAlmostEqual(0.00, displacement_truss_tip[0])
        self.assertAlmostEqual(disp_y_analytical, displacement_truss_tip[1])
        self.assertAlmostEqual(0.00, displacement_truss_tip[2])
   
    def test_cr_beam_linear(self):
        dim = 3
        nr_nodes = 2
        nr_elements = nr_nodes-1
        mp = KratosMultiphysics.ModelPart("solid_part")
        self._add_variables(mp)
        self._apply_material_properties_beam(mp,dim)

        #create nodes
        dx = 2.00 / nr_elements
        for i in range(nr_nodes):
            mp.CreateNewNode(i+1,i*dx,0.00,0.00)
        #add dofs
        self._add_dofs(mp)
        #create condition
        mp.CreateNewCondition("LineLoadCondition3D2N",1,[1,nr_nodes],mp.GetProperties()[0])  
        #create submodelparts for dirichlet boundary conditions
        bcs_xyz = mp.CreateSubModelPart("Dirichlet_XYZ")
        bcs_xyz.AddNodes([1])
        bcs_rot = mp.CreateSubModelPart("Dirichlet_RotAll")
        bcs_rot.AddNodes([1])    
        #create a submodalpart for neumann boundary conditions
        bcs_neumann = mp.CreateSubModelPart("LineLoad3D_neumann")
        bcs_neumann.AddNodes([1,nr_nodes])
        bcs_neumann.AddConditions([1])             
        #create Element
        for i in range(nr_elements):
            mp.CreateNewElement("CrLinearBeamElement3D2N", i+1, [i+1,i+2],
             mp.GetProperties()[0])

        #apply boundary conditions
        Force_Y = -10000.00
        self._apply_BCs(bcs_xyz,'xyz')
        self._apply_BCs(bcs_rot,'rotXYZ')

        self._apply_Neumann_BCs(bcs_neumann,'y',Force_Y)

        
        #solve + compare
        strategy = self._solve_linear(mp)  
        strategy.Solve()    
        self._check_results_linear_beam(mp,nr_nodes)    


    def test_truss_linear(self):
        dim = 3
        nr_nodes = 3
        nr_elements = 2
        mp = KratosMultiphysics.ModelPart("solid_part")
        self._add_variables(mp)
        self._apply_material_properties_truss(mp,dim)

        #create nodes
        mp.CreateNewNode(1,0.00,0.00,0.00)
        mp.CreateNewNode(2,1.00,1.00,0.00)
        mp.CreateNewNode(3,2.00,0.00,0.00)

        #add dofs
        self._add_dofs(mp)
        #create condition
        mp.CreateNewCondition("LineLoadCondition3D2N",1,[1,2],mp.GetProperties()[0])  
        mp.CreateNewCondition("LineLoadCondition3D2N",2,[2,3],mp.GetProperties()[0])  
        #create submodelparts for dirichlet boundary conditions
        bcs_xyz = mp.CreateSubModelPart("Dirichlet_XYZ")
        bcs_xyz.AddNodes([1,3])  
        bcs_xz = mp.CreateSubModelPart("Dirichlet_XZ")
        bcs_xz.AddNodes([2])  
        #create a submodalpart for neumann boundary conditions
        bcs_neumann = mp.CreateSubModelPart("LineLoad3D_neumann")
        bcs_neumann.AddNodes([1,2,3])
        bcs_neumann.AddConditions([1,2])             
        #create Element
        for i in range(nr_elements):
            mp.CreateNewElement("TrussLinearElement3D2N", i+1, [i+1,i+2],
             mp.GetProperties()[0])

        #apply boundary conditions
        Force_Y = -10000.00
        self._apply_BCs(bcs_xyz,'xyz')
        self._apply_BCs(bcs_xz,'xz')

        self._apply_Neumann_BCs(bcs_neumann,'y',Force_Y)

        
        #solve + compare
        strategy = self._solve_linear(mp)  
        strategy.Solve()   
        self._check_results_linear_truss(mp)    

if __name__ == '__main__':
    KratosUnittest.main()


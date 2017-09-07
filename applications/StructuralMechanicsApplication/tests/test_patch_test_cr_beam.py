from __future__ import print_function, absolute_import, division
import KratosMultiphysics 

import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import KratosMultiphysics.KratosUnittest as KratosUnittest


from math import sqrt, sin, cos, pi, exp, atan

class TestCrBeam3D2N(KratosUnittest.TestCase):
    def setUp(self):
        pass

    def _add_dofs(self,mp):
        for node in mp.Nodes:
            node.AddDof(KratosMultiphysics.DISPLACEMENT_X)
            node.AddDof(KratosMultiphysics.DISPLACEMENT_Y)
            node.AddDof(KratosMultiphysics.DISPLACEMENT_Z)
            node.AddDof(KratosMultiphysics.REACTION_X)
            node.AddDof(KratosMultiphysics.REACTION_Y)  
            node.AddDof(KratosMultiphysics.REACTION_Z)
            node.AddDof(KratosMultiphysics.ROTATION_X)
            node.AddDof(KratosMultiphysics.ROTATION_Y)
            node.AddDof(KratosMultiphysics.ROTATION_Z)
            node.AddDof(KratosMultiphysics.TORQUE_X)
            node.AddDof(KratosMultiphysics.TORQUE_Y)
            node.AddDof(KratosMultiphysics.TORQUE_Z)
            
    
    def _add_variables(self,mp):
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.ROTATION)
        mp.AddNodalSolutionStepVariable(StructuralMechanicsApplication.POINT_LOAD)  
        mp.AddNodalSolutionStepVariable(StructuralMechanicsApplication.POINT_MOMENT)  
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.VOLUME_ACCELERATION) 
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.ACCELERATION)     
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.ANGULAR_VELOCITY)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.ANGULAR_ACCELERATION)    
        mp.AddNodalSolutionStepVariable(StructuralMechanicsApplication.LOCAL_POINT_MOMENT)
        mp.AddNodalSolutionStepVariable(StructuralMechanicsApplication.POINT_TORQUE)
        mp.AddNodalSolutionStepVariable(StructuralMechanicsApplication.LOCAL_POINT_MOMENT)


        
    def _apply_material_properties(self,mp,dim):
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
        if (which_dof == 'yz'):
            for node in mp.Nodes:
                node.Fix(KratosMultiphysics.DISPLACEMENT_Y)
                node.Fix(KratosMultiphysics.DISPLACEMENT_Z)
        if (which_dof == 'rotXYZ'):
            for node in mp.Nodes:
                node.Fix(KratosMultiphysics.ROTATION_X)
                node.Fix(KratosMultiphysics.ROTATION_Y)
                node.Fix(KratosMultiphysics.ROTATION_Z)                

    def _apply_Neumann_BCs(self,mp,which_dof,load_size_dir):

        if(which_dof == 'y'):
            for node in mp.Nodes:
                node.SetSolutionStepValue(StructuralMechanicsApplication.
                POINT_LOAD_Y,0,load_size_dir)        
        if(which_dof == 'x'):
            for node in mp.Nodes:
                node.SetSolutionStepValue(StructuralMechanicsApplication.
                POINT_LOAD_X,0,load_size_dir)  

        
    def _solve_linear(self,mp):
        
        linear_solver = KratosMultiphysics.SkylineLUFactorizationSolver()
        builder_and_solver = KratosMultiphysics.ResidualBasedBlockBuilderAndSolver(linear_solver)
        scheme = KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme()
        
        compute_reactions = False  #False --> working
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

    def _solve_nonlinear(self,mp):
        linear_solver = KratosMultiphysics.SkylineLUFactorizationSolver()
        builder_and_solver = KratosMultiphysics.ResidualBasedBlockBuilderAndSolver(linear_solver)
        scheme = KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme()
        convergence_criterion = StructuralMechanicsApplication.DisplacementAndOtherDoFCriteria(1e-15,1e-15)
        convergence_criterion.SetEchoLevel(0)
        
        max_iters = 1000
        compute_reactions = False
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
    def _solve_dynamic(self,mp):
        
        #define a minimal newton raphson solver
        linear_solver = KratosMultiphysics.SkylineLUFactorizationSolver()
        builder_and_solver = KratosMultiphysics.ResidualBasedBlockBuilderAndSolver(linear_solver)
        scheme = KratosMultiphysics.ResidualBasedBossakDisplacementScheme(0.00)
        convergence_criterion = KratosMultiphysics.ResidualCriteria(1e-8,1e-8)
        convergence_criterion.SetEchoLevel(0)


        max_iters = 1000
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
        
           
    
    def _check_results_linear(self,mp,endNode):
        
        #check displacement result
        displacement_cantilever_tip = mp.Nodes[endNode].GetSolutionStepValue(
            KratosMultiphysics.DISPLACEMENT)

        disp_y_analytical = -400000.00*(1.2**3)/(3*210e9*0.00001)

        self.assertAlmostEqual(0.00, displacement_cantilever_tip[0])
        self.assertAlmostEqual(disp_y_analytical, displacement_cantilever_tip[1])
        self.assertAlmostEqual(0.00, displacement_cantilever_tip[2])

    def _check_results_nonlinear(self,mp,timestep,Moment_i,endNode):
        ##node at cantilever tip
        node_temp = mp.Nodes[endNode]

        displacement_x = node_temp.GetSolutionStepValue(
            KratosMultiphysics.DISPLACEMENT_X)
        displacement_y = node_temp.GetSolutionStepValue(
            KratosMultiphysics.DISPLACEMENT_Y)
        moment_z = node_temp.GetSolutionStepValue(StructuralMechanicsApplication.
                POINT_MOMENT_Z)

        #check moment z
        self.assertAlmostEqual(moment_z, Moment_i)
        #check displacement as soon as a total circle is formed
        #M = EI * 2 * pi / L ---> 13200000.0 reached at t_step = 527
        if (timestep == 527):
            self.assertAlmostEqual(displacement_x, -1.00,2)
            self.assertAlmostEqual(displacement_y, 0.00,5)

    def _check_results_dynamic(self,mp,time_i):
        
        #analaytical free-vibration node 3
        we1 = 7917.25
        we2 = 19113.94
        y1 = 1.4142*2.874e-5
        y2 = -1.4142*4.93107e-6
        test_disp_temp = y1*(1-cos(we1*time_i))-y2*(1-cos(we2*time_i))
        simulated_disp_temp = mp.Nodes[3].GetSolutionStepValue(
            KratosMultiphysics.DISPLACEMENT_X)  

        self.assertAlmostEqual(simulated_disp_temp, test_disp_temp,6)

        #analaytical free-vibration node 2
        we1 = 7917.25
        we2 = 19113.94
        y1 = 1.000*2.874e-5
        y2 = 1.000*4.93107e-6
        test_disp_temp = y1*(1-cos(we1*time_i))-y2*(1-cos(we2*time_i))
        simulated_disp_temp = mp.Nodes[2].GetSolutionStepValue(
            KratosMultiphysics.DISPLACEMENT_X)  

        self.assertAlmostEqual(simulated_disp_temp, test_disp_temp,6)
    def _set_and_fill_buffer(self,mp,buffer_size,delta_time):
        # Set buffer size
        mp.SetBufferSize(buffer_size)

        # Fill buffer
        time = mp.ProcessInfo[KratosMultiphysics.TIME]
        time = time - delta_time * (buffer_size)
        mp.ProcessInfo.SetValue(KratosMultiphysics.TIME, time)
        for size in range(0, buffer_size):
            step = size - (buffer_size -1)
            mp.ProcessInfo.SetValue(KratosMultiphysics.STEP, step)
            time = time + delta_time
            #delta_time is computed from previous time in process_info
            mp.CloneTimeStep(time)

        mp.ProcessInfo[KratosMultiphysics.IS_RESTARTED] = False        

    def test_cr_beam_linear(self):
        dim = 3
        nr_nodes = 11
        nr_elements = nr_nodes-1
        mp = KratosMultiphysics.ModelPart("solid_part")
        self._add_variables(mp)
        self._apply_material_properties(mp,dim)

        #create nodes
        dx = 1.20 / nr_elements
        for i in range(nr_nodes):
            mp.CreateNewNode(i+1,i*dx,0.00,0.00)
        #add dofs
        self._add_dofs(mp)
        #create condition
        mp.CreateNewCondition("PointLoadCondition3D1N",1,[nr_nodes],mp.GetProperties()[0])  
        #create submodelparts for dirichlet boundary conditions
        bcs_xyz = mp.CreateSubModelPart("Dirichlet_XYZ")
        bcs_xyz.AddNodes([1])
        bcs_rot = mp.CreateSubModelPart("Dirichlet_RotAll")
        bcs_rot.AddNodes([1])    
        #create a submodalpart for neumann boundary conditions
        bcs_neumann = mp.CreateSubModelPart("PointLoad3D_neumann")
        bcs_neumann.AddNodes([nr_nodes])
        bcs_neumann.AddConditions([1])             
        #create Element
        for i in range(nr_elements):
            mp.CreateNewElement("CrLinearBeamElement3D2N", i+1, [i+1,i+2],
             mp.GetProperties()[0])

        #apply boundary conditions
        Force_Y = -400000.00
        self._apply_BCs(bcs_xyz,'xyz')
        self._apply_BCs(bcs_rot,'rotXYZ')

        self._apply_Neumann_BCs(bcs_neumann,'y',Force_Y)

        
        #solve + compare
        self._solve_linear(mp)    
        self._check_results_linear(mp,nr_nodes)
    def test_cr_beam_nonlinear(self):
        dim = 3
        nr_nodes = 21
        nr_elements = nr_nodes-1
        mp = KratosMultiphysics.ModelPart("solid_part")
        self._add_variables(mp)
        self._apply_material_properties(mp,dim)

        #create nodes
        dx = 1.00 / nr_elements
        for i in range(nr_nodes):
            mp.CreateNewNode(i+1,i*dx,0.00,0.00)
        #add dofs
        self._add_dofs(mp)
        #create condition
        mp.CreateNewCondition("PointMomentCondition3D1N",1,[nr_nodes],mp.GetProperties()[0])  
        #create submodelparts for dirichlet boundary conditions
        bcs_xyz = mp.CreateSubModelPart("Dirichlet_XYZ")
        bcs_xyz.AddNodes([1])
        bcs_rot = mp.CreateSubModelPart("Dirichlet_RotAll")
        bcs_rot.AddNodes([1])               
        #create Element
        for i in range(nr_elements):
            mp.CreateNewElement("CrBeamElement3D2N", i+1, [i+1,i+2],
             mp.GetProperties()[0])

        #apply constant boundary conditions
        self._apply_BCs(bcs_xyz,'xyz')
        self._apply_BCs(bcs_rot,'rotXYZ')


        #incrementally increase load -> nonlinear case
        Moment_Z = 25000.00
        time_start = 0.00
        time_end = 528
        time_delta = 1
        time_i = time_start
        time_step = 0
        while (time_i <= time_end):
            time_i += time_delta
            #apply non-constant boundary conditions
            Moment_i = Moment_Z*time_i
            mp.Nodes[nr_nodes].SetSolutionStepValue(StructuralMechanicsApplication.
                POINT_MOMENT_Z,0,Moment_i)  
            #solve + compare
            self._solve_nonlinear(mp)
            self._check_results_nonlinear(mp,time_step,Moment_i,nr_nodes)
            time_step += 1

  
        
  

if __name__ == '__main__':
    KratosUnittest.main()


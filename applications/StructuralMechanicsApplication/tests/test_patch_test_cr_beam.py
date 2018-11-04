from __future__ import print_function, absolute_import, division
import KratosMultiphysics

import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import KratosMultiphysics.KratosUnittest as KratosUnittest


from math import sqrt, sin, cos, pi, exp, atan

class TestCrBeam3D2N(KratosUnittest.TestCase):
    def setUp(self):
        pass

    def _add_dofs(self,mp):
        # Adding dofs AND their corresponding reactions
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.REACTION_X,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Y, KratosMultiphysics.REACTION_Y,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Z, KratosMultiphysics.REACTION_Z,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ROTATION_X, KratosMultiphysics.REACTION_MOMENT_X,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ROTATION_Y, KratosMultiphysics.REACTION_MOMENT_Y,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ROTATION_Z, KratosMultiphysics.REACTION_MOMENT_Z,mp)


    def _add_variables(self,mp):
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.ROTATION)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION_MOMENT)
        mp.AddNodalSolutionStepVariable(StructuralMechanicsApplication.POINT_LOAD)
        mp.AddNodalSolutionStepVariable(StructuralMechanicsApplication.POINT_MOMENT)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.VOLUME_ACCELERATION)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.ANGULAR_VELOCITY)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.ACCELERATION)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.ANGULAR_ACCELERATION)

    def _apply_material_properties(self,mp,dim):
        #define properties
        mp.GetProperties()[0].SetValue(KratosMultiphysics.YOUNG_MODULUS,210e9)
        mp.GetProperties()[0].SetValue(KratosMultiphysics.DENSITY,7850)
        mp.GetProperties()[0].SetValue(StructuralMechanicsApplication.CROSS_AREA,0.01)
        mp.GetProperties()[0].SetValue(KratosMultiphysics.POISSON_RATIO,0.30)
        mp.GetProperties()[0].SetValue(StructuralMechanicsApplication.TORSIONAL_INERTIA,0.00001)
        mp.GetProperties()[0].SetValue(StructuralMechanicsApplication.I22,0.00001)
        mp.GetProperties()[0].SetValue(StructuralMechanicsApplication.I33,0.00001)

        g = [0,0,0]
        mp.GetProperties()[0].SetValue(KratosMultiphysics.VOLUME_ACCELERATION,g)

        cl = StructuralMechanicsApplication.LinearElastic3DLaw()
        mp.GetProperties()[0].SetValue(KratosMultiphysics.CONSTITUTIVE_LAW,cl)

    def _apply_elemental_data(self,element):
        # Adding LOCAL_AXIS_2
        element.SetValue(KratosMultiphysics.LOCAL_AXIS_2,[0,1,0])

    def _apply_3D_moment_hinge_z(self,element):
        # Adding LOCAL_AXIS_2
        element.SetValue(StructuralMechanicsApplication.CONDENSED_DOF_LIST,[11])


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
                POINT_LOAD_Y, load_size_dir, mp.Nodes)
            # for node in mp.Nodes:
            #     node.SetSolutionStepValue(StructuralMechanicsApplication.
            #     POINT_LOAD_Y,0,load_size_dir)
        if(which_dof == 'x'):
            KratosMultiphysics.VariableUtils().SetScalarVar(StructuralMechanicsApplication.
                POINT_LOAD_X, load_size_dir, mp.Nodes)
            # for node in mp.Nodes:
            #     node.SetSolutionStepValue(StructuralMechanicsApplication.
            #     POINT_LOAD_X,0,load_size_dir)

    def _solve_linear(self,mp):
        linear_solver = KratosMultiphysics.SkylineLUFactorizationSolver()
        builder_and_solver = KratosMultiphysics.ResidualBasedBlockBuilderAndSolver(linear_solver)
        scheme = KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme()

        compute_reactions = True  #Now the rotation reactions (REACTION_MOMENT) is added, so it works
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

        self.assertAlmostEqual(0.00, displacement_cantilever_tip[0],6)
        self.assertAlmostEqual(disp_y_analytical, displacement_cantilever_tip[1],6)
        self.assertAlmostEqual(0.00, displacement_cantilever_tip[2],6)

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
        if (timestep == 5):
            self.assertAlmostEqual(displacement_x, -0.0008214481371826507,5)
            self.assertAlmostEqual(displacement_y, 0.03560205297258129,5)


    def _check_results_dynamic(self,mp,time_i,nr_nodes,time_step):
        #check free vibration of cantilever tip
        disp_y_simulated = mp.Nodes[nr_nodes].GetSolutionStepValue(
            KratosMultiphysics.DISPLACEMENT_Y)
        disp_y_analytical = [-4.4017262561983686e-05,-0.00018621779467051006,
        -0.00040709297076666834,-0.0006775988708011861,-0.0009923249270175282]

        self.assertAlmostEqual(disp_y_analytical[time_step], disp_y_simulated)

    def _check_results_dynamic_lumped(self,mp,time_i,nr_nodes,time_step):
        #check free vibration of cantilever tip
        disp_y_simulated = mp.Nodes[nr_nodes].GetSolutionStepValue(
            KratosMultiphysics.DISPLACEMENT_Y)
        disp_y_analytical = [-4.162518390580818e-05,-0.00017969144438005632,
        -0.00039846788371390653,-0.0006674048593190372,-0.000980511641724115]

        self.assertAlmostEqual(disp_y_analytical[time_step], disp_y_simulated)

    def _check_results_dynamic_explicit(self,mp,time_i,nr_nodes,time_step):
        #check free vibration of cantilever tip
        disp_y_simulated = mp.Nodes[nr_nodes].GetSolutionStepValue(
            KratosMultiphysics.DISPLACEMENT_Y)
        disp_y_analytical = [-2.8662420382165618e-06,-6.5284595740263744e-06,
            -8.667247023317519e-06,-1.5193910085875024e-05,-2.4013659391608225e-05,
            -3.2183559803816915e-05,-4.43709426806417e-05,-5.9041345875279067e-05,
            -7.330044103901265e-05,-9.107969267628107e-05,-0.00011105794571720542,
            -0.0001304810033530728,-0.00015308391062689988,-0.00017782947133379302,
            -0.00020200606491067104,-0.00022918489563129963,-0.00025850356356149775,
            -0.00028714489980093377,-0.0003185138569101933,-0.0003519411006451715,
            -0.0003845236264241205,-0.000419580648891076,-0.0004566985915125524,
            -0.0004929004674252967,-0.0005314362205181631,-0.0005721203453514332,
            -0.0006117811911574032]

        self.assertAlmostEqual(disp_y_analytical[time_step], disp_y_simulated,6)


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
        current_model = KratosMultiphysics.Model()
        mp = current_model.CreateModelPart("solid_part")
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

    def test_cr_beam_linear_local_axis2(self):
            dim = 3
            nr_nodes = 11
            nr_elements = nr_nodes-1
            current_model = KratosMultiphysics.Model()
            mp = current_model.CreateModelPart("solid_part")
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

            #apply local_axis_2 elemental data
            for i in range(nr_elements): self._apply_elemental_data(mp.GetElement(i+1))
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
        current_model = KratosMultiphysics.Model()
        mp = current_model.CreateModelPart("solid_part")
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
        time_end = 5
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

    def test_cr_beam_dynamic_lumped_mass_matrix(self):
        dim = 3
        nr_nodes = 11
        nr_elements = nr_nodes-1
        current_model = KratosMultiphysics.Model()
        mp = current_model.CreateModelPart("solid_part")
        self._add_variables(mp)
        self._apply_material_properties(mp,dim)
        mp.GetProperties()[0].SetValue(StructuralMechanicsApplication.USE_CONSISTENT_MASS_MATRIX,False)

        #create nodes
        dx = 1.00 / nr_elements
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
            mp.CreateNewElement("CrBeamElement3D2N", i+1, [i+1,i+2],
             mp.GetProperties()[0])

        #apply constant boundary conditions
        self._apply_BCs(bcs_xyz,'xyz')
        self._apply_BCs(bcs_rot,'rotXYZ')
        Force_Y = -100000.000
        self._apply_Neumann_BCs(bcs_neumann,'y',Force_Y)

        #loop over time
        time_start = 0.00
        time_end = 0.0004
        # time_delta = 0.001
        time_delta = 0.0001
        time_i = time_start
        time_step = 0
        self._set_and_fill_buffer(mp,2,time_delta)

        x = []
        y = []
        y_1 = []
        while (time_i <= time_end):

            time_i += time_delta
            mp.CloneTimeStep(time_i)
            #solve + compare
            self._solve_dynamic(mp)
            self._check_results_dynamic_lumped(mp,time_i,nr_nodes,time_step)
            time_step += 1

    def test_cr_beam_dynamic_consistent_mass_matrix(self):
        dim = 3
        nr_nodes = 11
        nr_elements = nr_nodes-1
        current_model = KratosMultiphysics.Model()
        mp = current_model.CreateModelPart("solid_part")
        self._add_variables(mp)
        self._apply_material_properties(mp,dim)
        mp.GetProperties()[0].SetValue(StructuralMechanicsApplication.USE_CONSISTENT_MASS_MATRIX,True)

        #create nodes
        dx = 1.00 / nr_elements
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
            mp.CreateNewElement("CrBeamElement3D2N", i+1, [i+1,i+2],
             mp.GetProperties()[0])

        #apply constant boundary conditions
        self._apply_BCs(bcs_xyz,'xyz')
        self._apply_BCs(bcs_rot,'rotXYZ')
        Force_Y = -100000.000
        self._apply_Neumann_BCs(bcs_neumann,'y',Force_Y)

        #loop over time
        time_start = 0.00
        time_end = 0.0004
        # time_delta = 0.001
        time_delta = 0.0001
        time_i = time_start
        time_step = 0
        self._set_and_fill_buffer(mp,2,time_delta)

        x = []
        y = []
        y_1 = []
        while (time_i <= time_end):

            time_i += time_delta
            mp.CloneTimeStep(time_i)
            #solve + compare
            self._solve_dynamic(mp)
            self._check_results_dynamic(mp,time_i,nr_nodes,time_step)
            time_step += 1

    def test_cr_beam_dynamic_explicit(self):
            dim = 3
            nr_nodes = 11
            nr_elements = nr_nodes-1
            current_model = KratosMultiphysics.Model()
            mp = current_model.CreateModelPart("solid_part")
            self._add_variables(mp)
            _add_explicit_variables(mp)
            self._apply_material_properties(mp,dim)
            mp.GetProperties()[0].SetValue(StructuralMechanicsApplication.USE_CONSISTENT_MASS_MATRIX,True)

            #create nodes
            dx = 1.00 / nr_elements
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

            #apply constant boundary conditions
            self._apply_BCs(bcs_xyz,'xyz')
            self._apply_BCs(bcs_rot,'rotXYZ')
            Force_Y = -100000.000
            self._apply_Neumann_BCs(bcs_neumann,'y',Force_Y)

            #loop over time
            time_start = 0.00
            time_delta = 0.000015
            time_end = time_delta*12
            # time_delta = 0.001
            time_i = time_start
            time_step = 0
            self._set_and_fill_buffer(mp,2,time_delta)
            strategy_expl = _create_dynamic_explicit_strategy(mp)
            while (time_i <= time_end):

                time_i += time_delta
                mp.CloneTimeStep(time_i)
                #solve + compare
                strategy_expl.Solve()
                self._check_results_dynamic_explicit(mp,time_i,nr_nodes,time_step)

                time_step += 1

    def test_cr_beam_linear_moment_hinge(self):
        dim = 2
        nr_nodes = 3
        nr_elements = nr_nodes-1
        current_model = KratosMultiphysics.Model()
        mp = current_model.CreateModelPart("solid_part")
        self._add_variables(mp)
        self._apply_material_properties(mp,dim)

        #create nodes
        dx = 2.20 / nr_elements
        for i in range(nr_nodes):
            mp.CreateNewNode(i+1,i*dx,0.00,0.00)
        #add dofs
        self._add_dofs(mp)
        #create condition
        mp.CreateNewCondition("PointLoadCondition3D1N",1,[nr_nodes-1],mp.GetProperties()[0])
        #create submodelparts for dirichlet boundary conditions
        bcs_xyz = mp.CreateSubModelPart("Dirichlet_XYZ")
        bcs_xyz.AddNodes([1,nr_nodes])
        bcs_rot = mp.CreateSubModelPart("Dirichlet_RotAll")
        bcs_rot.AddNodes([1,nr_nodes])
        #create a submodalpart for neumann boundary conditions
        bcs_neumann = mp.CreateSubModelPart("PointLoad3D_neumann")
        bcs_neumann.AddNodes([nr_nodes-1])
        bcs_neumann.AddConditions([1])
        #create Element
        for i in range(nr_elements):
            mp.CreateNewElement("CrLinearBeamElement3D2N", i+1, [i+1,i+2],
             mp.GetProperties()[0])

        #apply condensation elemental data
        self._apply_3D_moment_hinge_z(mp.GetElement(1))

        #apply boundary conditions
        Force_Y = -400000.00
        self._apply_BCs(bcs_xyz,'xyz')
        self._apply_BCs(bcs_rot,'rotXYZ')

        self._apply_Neumann_BCs(bcs_neumann,'y',Force_Y)

        #solve + compare
        self._solve_linear(mp)

        displacement_mid_field = mp.Nodes[nr_nodes-1].GetSolutionStepValue(
            KratosMultiphysics.DISPLACEMENT_Y)
        self.assertAlmostEqual(displacement_mid_field, -0.04225396825396825)


        out1 = mp.Elements[1].CalculateOnIntegrationPoints(KratosMultiphysics.MOMENT,mp.ProcessInfo)
        out2 = mp.Elements[2].CalculateOnIntegrationPoints(KratosMultiphysics.MOMENT,mp.ProcessInfo)
        self.assertAlmostEqual(out1[0][2], 165000.0)
        self.assertAlmostEqual(out1[1][2], 110000.0)
        self.assertAlmostEqual(out1[2][2], 55000.0)
        self.assertAlmostEqual(out2[2][2], 165000.0)
        self.assertAlmostEqual(out2[1][2], 110000.0)
        self.assertAlmostEqual(out2[0][2], 55000.0)


class TestCrBeam2D2N(KratosUnittest.TestCase):
    def setUp(self):
        pass

    def _add_dofs(self,mp):
        # Adding dofs AND their corresponding reactions
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.REACTION_X,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Y, KratosMultiphysics.REACTION_Y,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Z, KratosMultiphysics.REACTION_Y,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ROTATION_X, KratosMultiphysics.REACTION_MOMENT_X,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ROTATION_Y, KratosMultiphysics.REACTION_MOMENT_Y,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ROTATION_Z, KratosMultiphysics.REACTION_MOMENT_Z,mp)


    def _add_variables(self,mp):
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.ROTATION)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION_MOMENT)
        mp.AddNodalSolutionStepVariable(StructuralMechanicsApplication.POINT_LOAD)
        mp.AddNodalSolutionStepVariable(StructuralMechanicsApplication.POINT_MOMENT)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.VOLUME_ACCELERATION)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.ANGULAR_VELOCITY)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.ACCELERATION)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.ANGULAR_ACCELERATION)

    def _apply_material_properties(self,mp,dim):
        #define properties
        mp.GetProperties()[0].SetValue(KratosMultiphysics.YOUNG_MODULUS,210e9)
        mp.GetProperties()[0].SetValue(KratosMultiphysics.DENSITY,7850)
        mp.GetProperties()[0].SetValue(StructuralMechanicsApplication.CROSS_AREA,0.01)
        mp.GetProperties()[0].SetValue(KratosMultiphysics.POISSON_RATIO,0.30)
        mp.GetProperties()[0].SetValue(StructuralMechanicsApplication.I33,0.00001)

        g = [0,0,0]
        mp.GetProperties()[0].SetValue(KratosMultiphysics.VOLUME_ACCELERATION,g)

        cl = StructuralMechanicsApplication.LinearElastic3DLaw()
        mp.GetProperties()[0].SetValue(KratosMultiphysics.CONSTITUTIVE_LAW,cl)

    def _apply_elemental_data(self,element):
        # Adding LOCAL_AXIS_2
        element.SetValue(KratosMultiphysics.LOCAL_AXIS_2,(0,1,0))

    def _apply_2D_moment_hinge_z(self,element):
        # Adding LOCAL_AXIS_2
        element.SetValue(StructuralMechanicsApplication.CONDENSED_DOF_LIST,[5])

    def _apply_BCs(self,mp,which_dof):
        if (which_dof == 'xy'):
            KratosMultiphysics.VariableUtils().ApplyFixity(KratosMultiphysics.DISPLACEMENT_X, True, mp.Nodes)
            KratosMultiphysics.VariableUtils().ApplyFixity(KratosMultiphysics.DISPLACEMENT_Y, True, mp.Nodes)
        if (which_dof == 'x'):
            KratosMultiphysics.VariableUtils().ApplyFixity(KratosMultiphysics.DISPLACEMENT_X, True, mp.Nodes)
        if (which_dof == 'y'):
            KratosMultiphysics.VariableUtils().ApplyFixity(KratosMultiphysics.DISPLACEMENT_Y, True, mp.Nodes)
        if (which_dof == 'rotZ'):
            KratosMultiphysics.VariableUtils().ApplyFixity(KratosMultiphysics.ROTATION_Z, True, mp.Nodes)

    def _apply_Neumann_BCs(self,mp,which_dof,load_size_dir):
        if(which_dof == 'y'):
            KratosMultiphysics.VariableUtils().SetScalarVar(StructuralMechanicsApplication.
                POINT_LOAD_Y, load_size_dir, mp.Nodes)
            # for node in mp.Nodes:
            #     node.SetSolutionStepValue(StructuralMechanicsApplication.
            #     POINT_LOAD_Y,0,load_size_dir)
        if(which_dof == 'x'):
            KratosMultiphysics.VariableUtils().SetScalarVar(StructuralMechanicsApplication.
                POINT_LOAD_X, load_size_dir, mp.Nodes)
            # for node in mp.Nodes:
            #     node.SetSolutionStepValue(StructuralMechanicsApplication.
            #     POINT_LOAD_X,0,load_size_dir)

    def _solve_linear(self,mp):
        linear_solver = KratosMultiphysics.SkylineLUFactorizationSolver()
        builder_and_solver = KratosMultiphysics.ResidualBasedBlockBuilderAndSolver(linear_solver)
        scheme = KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme()

        compute_reactions = True  #Now the rotation reactions (REACTION_MOMENT) is added, so it works
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
        if (timestep == 5):
            self.assertAlmostEqual(displacement_x, -0.0008258150820168632)
            self.assertAlmostEqual(displacement_y, 0.03569588131158605)

    def _check_results_dynamic(self,mp,time_i,nr_nodes,time_step):
        #check free vibration of cantilever tip
        disp_y_analytical = [-4.769300755191701e-05,-0.00019510832494467397,
        -0.0004182462789190604,-0.0006910029789562587,-0.0010073938697706265]

        disp_y_simulated = mp.Nodes[nr_nodes].GetSolutionStepValue(
            KratosMultiphysics.DISPLACEMENT_Y)

        self.assertAlmostEqual(disp_y_analytical[time_step], disp_y_simulated)


    def _check_results_dynamic_explicit(self,mp,time_i,nr_nodes,time_step):
        #check free vibration of cantilever tip
        disp_y_simulated = mp.Nodes[nr_nodes].GetSolutionStepValue(
            KratosMultiphysics.DISPLACEMENT_Y)
        disp_y_analytical = [-2.8662420382165618e-06,-6.5284598580628135e-06,
            -8.661766714284483e-06,-1.5166844659226568e-05,-2.3974692146734503e-05,
            -3.213060541808351e-05,-4.4284460087458615e-05,-5.8937873140393946e-05,
            -7.318106720395262e-05,-9.09215604293434e-05,-0.00011088212163774735,
            -0.00013029114264297788,-0.00015285402390510065,-0.00017758464802435417,
            -0.0002017526913999827,-0.00022889337656467832,-0.00025820171682640884,
            -0.0002868413100380328,-0.000318174460796471,-0.0003515951807656998,
            -0.00038418182365403656,-0.0004192072728098697,-0.0004563324503910268,
            -0.0004925439930684863,-0.0005309736146781797,-0.0005715335001692381,
            -0.0006112486419960861]


        self.assertAlmostEqual(disp_y_analytical[time_step], disp_y_simulated)

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
        dim = 2
        nr_nodes = 11
        nr_elements = nr_nodes-1
        current_model = KratosMultiphysics.Model()
        mp = current_model.CreateModelPart("solid_part")
        self._add_variables(mp)
        self._apply_material_properties(mp,dim)

        #create nodes
        dx = 1.20 / nr_elements
        for i in range(nr_nodes):
            mp.CreateNewNode(i+1,i*dx,0.00,0.00)
        #add dofs
        self._add_dofs(mp)
        #create condition
        mp.CreateNewCondition("PointLoadCondition2D1N",1,[nr_nodes],mp.GetProperties()[0])
        #create submodelparts for dirichlet boundary conditions
        bcs_xyz = mp.CreateSubModelPart("Dirichlet_XY")
        bcs_xyz.AddNodes([1])
        bcs_rot = mp.CreateSubModelPart("Dirichlet_Rot")
        bcs_rot.AddNodes([1])
        #create a submodalpart for neumann boundary conditions
        bcs_neumann = mp.CreateSubModelPart("PointLoad3D_neumann")
        bcs_neumann.AddNodes([nr_nodes])
        bcs_neumann.AddConditions([1])
        #create Element
        for i in range(nr_elements):
            mp.CreateNewElement("CrLinearBeamElement2D2N", i+1, [i+1,i+2],
             mp.GetProperties()[0])

        #apply boundary conditions
        Force_Y = -400000.00
        self._apply_BCs(bcs_xyz,'xy')
        self._apply_BCs(bcs_rot,'rotZ')

        self._apply_Neumann_BCs(bcs_neumann,'y',Force_Y)

        #solve + compare
        self._solve_linear(mp)
        self._check_results_linear(mp,nr_nodes)



    def test_cr_beam_dynamic_consistent_mass_matrix(self):
        dim = 2
        nr_nodes = 11
        nr_elements = nr_nodes-1
        current_model = KratosMultiphysics.Model()
        mp = current_model.CreateModelPart("solid_part")
        self._add_variables(mp)
        self._apply_material_properties(mp,dim)

        #create nodes
        dx = 1.00 / nr_elements
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
            mp.CreateNewElement("CrBeamElement2D2N", i+1, [i+1,i+2],
             mp.GetProperties()[0])

        #apply constant boundary conditions
        self._apply_BCs(bcs_xyz,'xy')
        self._apply_BCs(bcs_rot,'rotZ')
        Force_Y = -100000.000
        self._apply_Neumann_BCs(bcs_neumann,'y',Force_Y)

        #loop over time
        time_start = 0.00
        time_end = 0.0004
        # time_delta = 0.001
        time_delta = 0.0001
        time_i = time_start
        time_step = 0
        self._set_and_fill_buffer(mp,2,time_delta)

        x = []
        y = []
        y_1 = []
        while (time_i <= time_end):

            time_i += time_delta
            mp.CloneTimeStep(time_i)
            #solve + compare
            self._solve_dynamic(mp)
            self._check_results_dynamic(mp,time_i,nr_nodes,time_step)
            time_step += 1

    def test_cr_beam_nonlinear(self):
        dim = 2
        nr_nodes = 21
        nr_elements = nr_nodes-1
        current_model = KratosMultiphysics.Model()
        mp = current_model.CreateModelPart("solid_part")
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
            mp.CreateNewElement("CrBeamElement2D2N", i+1, [i+1,i+2],
             mp.GetProperties()[0])

        #apply constant boundary conditions
        self._apply_BCs(bcs_xyz,'xy')
        self._apply_BCs(bcs_rot,'rotZ')

        #incrementally increase load -> nonlinear case
        Moment_Z = 25000.00
        time_start = 0.00
        time_end = 5
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


    def test_cr_beam_linear_moment_hinge(self):
        dim = 2
        nr_nodes = 3
        nr_elements = nr_nodes-1
        current_model = KratosMultiphysics.Model()
        mp = current_model.CreateModelPart("solid_part")
        self._add_variables(mp)
        self._apply_material_properties(mp,dim)

        #create nodes
        dx = 2.20 / nr_elements
        for i in range(nr_nodes):
            mp.CreateNewNode(i+1,i*dx,0.00,0.00)
        #add dofs
        self._add_dofs(mp)
        #create condition
        mp.CreateNewCondition("PointLoadCondition2D1N",1,[nr_nodes-1],mp.GetProperties()[0])
        #create submodelparts for dirichlet boundary conditions
        bcs_xyz = mp.CreateSubModelPart("Dirichlet_XY")
        bcs_xyz.AddNodes([1,nr_nodes])
        bcs_rot = mp.CreateSubModelPart("Dirichlet_Rot")
        bcs_rot.AddNodes([1,nr_nodes])
        #create a submodalpart for neumann boundary conditions
        bcs_neumann = mp.CreateSubModelPart("PointLoad3D_neumann")
        bcs_neumann.AddNodes([nr_nodes-1])
        bcs_neumann.AddConditions([1])
        #create Element
        for i in range(nr_elements):
            mp.CreateNewElement("CrLinearBeamElement2D2N", i+1, [i+1,i+2],
             mp.GetProperties()[0])

        #apply condensation elemental data
        self._apply_2D_moment_hinge_z(mp.GetElement(1))

        #apply boundary conditions
        Force_Y = -400000.00
        self._apply_BCs(bcs_xyz,'xy')
        self._apply_BCs(bcs_rot,'rotZ')

        self._apply_Neumann_BCs(bcs_neumann,'y',Force_Y)

        #solve + compare
        self._solve_linear(mp)

        displacement_mid_field = mp.Nodes[nr_nodes-1].GetSolutionStepValue(
            KratosMultiphysics.DISPLACEMENT_Y)
        self.assertAlmostEqual(displacement_mid_field, -0.04225396825396825)


        out1 = mp.Elements[1].CalculateOnIntegrationPoints(KratosMultiphysics.MOMENT,mp.ProcessInfo)
        out2 = mp.Elements[2].CalculateOnIntegrationPoints(KratosMultiphysics.MOMENT,mp.ProcessInfo)
        self.assertAlmostEqual(out1[0][2], 165000.0)
        self.assertAlmostEqual(out1[1][2], 110000.0)
        self.assertAlmostEqual(out1[2][2], 55000.0)
        self.assertAlmostEqual(out2[2][2], 165000.0)
        self.assertAlmostEqual(out2[1][2], 110000.0)
        self.assertAlmostEqual(out2[0][2], 55000.0)

    def test_cr_beam_dynamic_explicit(self):
            dim = 2
            nr_nodes = 11
            nr_elements = nr_nodes-1
            current_model = KratosMultiphysics.Model()
            mp = current_model.CreateModelPart("solid_part")
            self._add_variables(mp)
            _add_explicit_variables(mp)
            self._apply_material_properties(mp,dim)
            mp.GetProperties()[0].SetValue(StructuralMechanicsApplication.USE_CONSISTENT_MASS_MATRIX,True)

            #create nodes
            dx = 1.00 / nr_elements
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
                mp.CreateNewElement("CrBeamElement2D2N", i+1, [i+1,i+2],
                mp.GetProperties()[0])

            #apply constant boundary conditions
            self._apply_BCs(bcs_xyz,'xyz')
            self._apply_BCs(bcs_rot,'rotXYZ')
            Force_Y = -100000.000
            self._apply_Neumann_BCs(bcs_neumann,'y',Force_Y)

            #loop over time
            time_start = 0.00
            time_end = 0.0004
            # time_delta = 0.001
            time_delta = 0.000015
            time_i = time_start
            time_step = 0
            self._set_and_fill_buffer(mp,2,time_delta)
            strategy_expl = _create_dynamic_explicit_strategy(mp)
            while (time_i <= time_end):

                time_i += time_delta
                mp.CloneTimeStep(time_i)
                #solve + compare
                strategy_expl.Solve()
                self._check_results_dynamic_explicit(mp,time_i,nr_nodes,time_step)

                time_step += 1

def _add_explicit_variables(mp):
    mp.AddNodalSolutionStepVariable(StructuralMechanicsApplication.MIDDLE_VELOCITY)
    mp.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_MASS)
    mp.AddNodalSolutionStepVariable(KratosMultiphysics.FORCE_RESIDUAL)
    mp.AddNodalSolutionStepVariable(KratosMultiphysics.RESIDUAL_VECTOR)
    mp.AddNodalSolutionStepVariable(StructuralMechanicsApplication.MIDDLE_ANGULAR_VELOCITY)
    mp.AddNodalSolutionStepVariable(StructuralMechanicsApplication.NODAL_INERTIA)
    mp.AddNodalSolutionStepVariable(KratosMultiphysics.MOMENT_RESIDUAL)

def _create_dynamic_explicit_strategy(mp):
    scheme = StructuralMechanicsApplication.ExplicitCentralDifferencesScheme(0.00,0.00,0.00)

    strategy = StructuralMechanicsApplication.MechanicalExplicitStrategy(mp,scheme,0,0,1)
    strategy.SetEchoLevel(0)
    return strategy


if __name__ == '__main__':
    KratosUnittest.main()


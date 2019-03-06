from __future__ import print_function, absolute_import, division
import KratosMultiphysics

import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import KratosMultiphysics.KratosUnittest as KratosUnittest


class TestPatchTestLargeStrain(KratosUnittest.TestCase):
    def setUp(self):
        pass
    
    def _add_variables(self,mp):
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.VOLUME_ACCELERATION)

    def _apply_material_properties(self, mp, dim, small_strain = True):
        #define properties
        mp.GetProperties()[1].SetValue(KratosMultiphysics.YOUNG_MODULUS,210e9)
        mp.GetProperties()[1].SetValue(KratosMultiphysics.POISSON_RATIO,0.3)
        mp.GetProperties()[1].SetValue(KratosMultiphysics.THICKNESS,1.0)

        g = [0,0,0]
        mp.GetProperties()[1].SetValue(KratosMultiphysics.VOLUME_ACCELERATION,g)

        if(dim == 2):
            if (small_strain == True):
                cl = StructuralMechanicsApplication.LinearElasticPlaneStress2DLaw()
            else:
                cl = StructuralMechanicsApplication.HyperElasticPlaneStrain2DLaw()
        else:
            if (small_strain == True):
                cl = StructuralMechanicsApplication.LinearElastic3DLaw()
            else:
                cl = StructuralMechanicsApplication.HyperElastic3DLaw()
        mp.GetProperties()[1].SetValue(KratosMultiphysics.CONSTITUTIVE_LAW,cl)

    def _set_buffer(self,mp):
        buffer_size = 3
        mp.SetBufferSize(buffer_size)
        # Cycle the buffer. This sets all historical nodal solution step data to
        # the current value and initializes the time stepping in the process info.
        mp.ProcessInfo[KratosMultiphysics.DELTA_TIME] = 1.0
        delta_time = mp.ProcessInfo[KratosMultiphysics.DELTA_TIME]
        time = mp.ProcessInfo[KratosMultiphysics.TIME]
        step =-buffer_size
        time = time - delta_time * buffer_size
        mp.ProcessInfo.SetValue(KratosMultiphysics.TIME, time)
        for i in range(0, buffer_size):
            step = step + 1
            time = time + delta_time
            mp.ProcessInfo.SetValue(KratosMultiphysics.STEP, step)
            mp.CloneTimeStep(time)

    def _apply_BCs(self,mp,A,b):
        for node in mp.Nodes:
            node.Fix(KratosMultiphysics.DISPLACEMENT_X)
            node.Fix(KratosMultiphysics.DISPLACEMENT_Y)
            node.Fix(KratosMultiphysics.DISPLACEMENT_Z)

        for node in mp.Nodes:
            xvec = KratosMultiphysics.Vector(3)
            u = KratosMultiphysics.Vector()

            xvec[0] = node.X0
            xvec[1] = node.Y0
            xvec[2] = node.Z0

            u = A*xvec
            u += b

            node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT,0,u)

    def _define_movement(self,dim):
        if(dim == 2):
            #define the applied motion - the idea is that the displacement is defined as u = A*xnode + b
            #so that the displcement is linear and the exact F = I + A
            A = KratosMultiphysics.Matrix(3,3)
            A[0,0] =   0.10;  A[0,1] = 0.12; A[0,2] = 0.0
            A[1,0] = - 0.05;  A[1,1] = 0.07; A[1,2] = 0.0
            A[2,0] =   0.00;  A[2,1] = 0.0;  A[2,2] = 0.0

            b = KratosMultiphysics.Vector(3)
            b[0] =  0.05
            b[1] = -0.02
            b[2] =  0.00

        else:
            #define the applied motion - the idea is that the displacement is defined as u = A*xnode + b
            #so that the displcement is linear and the exact F = I + A
            A = KratosMultiphysics.Matrix(3,3)
            A[0,0] =   0.10; A[0,1] = 0.12; A[0,2] = 0.0
            A[1,0] = - 0.05; A[1,1] = 0.07; A[1,2] = 0.1
            A[2,0] = - 0.02; A[2,1] = 0.0;  A[2,2] = -0.3

            b = KratosMultiphysics.Vector(3)
            b[0] =  0.05
            b[1] = -0.02
            b[2] =  0.07

        return A,b

    def _solve(self,mp):

        #define a minimal newton raphson solver
        linear_solver = KratosMultiphysics.SkylineLUFactorizationSolver()
        builder_and_solver = KratosMultiphysics.ResidualBasedEliminationBuilderAndSolver(linear_solver)
        scheme = KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme()
        convergence_criterion = KratosMultiphysics.ResidualCriteria(1e-14,1e-20)
        convergence_criterion.SetEchoLevel(0)

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

    def _create_strategy(self, mp):
        #define a minimal newton raphson solver
        linear_solver = KratosMultiphysics.SkylineLUFactorizationSolver()
        builder_and_solver = KratosMultiphysics.ResidualBasedEliminationBuilderAndSolver(linear_solver)
        scheme = KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme()
        convergence_criterion = KratosMultiphysics.ResidualCriteria(1e-4,1e-9)
        convergence_criterion.SetEchoLevel(0)

        #max_iters = 1
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

        return strategy

    def _solve_with_strategy(self, strategy, lhs, step):
        strategy.Check()
        strategy.Initialize()
        strategy.InitializeSolutionStep()
        strategy.Predict()
        strategy.SolveSolutionStep()
        lhs = strategy.GetSystemMatrix()
        strategy.FinalizeSolutionStep()

    def _check_results(self,mp,A,b):

        ##check that the results are exact on the nodes
        for node in mp.Nodes:
            xvec = KratosMultiphysics.Vector(3)
            xvec[0] = node.X0
            xvec[1] = node.Y0
            xvec[2] = node.Z0

            u = A*xvec
            u += b

            d = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT)
            self.assertAlmostEqual(d[0], u[0])
            self.assertAlmostEqual(d[1], u[1])
            self.assertAlmostEqual(d[2], u[2])

    def _check_outputs(self,mp,A,dim):

        E = mp.GetProperties()[1].GetValue(KratosMultiphysics.YOUNG_MODULUS)
        NU =mp.GetProperties()[1].GetValue(KratosMultiphysics.POISSON_RATIO)

        #given the matrix A, the analytic deformation gradient is F+I
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
                    self.assertTrue((abs((reference_strain[i] - strain[i])/strain[i]) < 1.0e-4))

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
                    self.assertTrue((abs((reference_stress[i] - stress[i])/stress[i]) < 1.0e-4))

    def test_compare_TL_UL_2D_triangle(self):
        dim = 2

        bc_nodes = [1, 2]
        load_nodes = [3, 4]

        current_model = KratosMultiphysics.Model()
        tl_mp = current_model.CreateModelPart("tl_solid_part")

        self._add_variables(tl_mp)
        self._apply_material_properties(tl_mp, dim, False)

        # Create nodes
        tl_mp.CreateNewNode(1, 0.0, 0.0, 0.0)
        tl_mp.CreateNewNode(2, 1.0, 0.0, 0.0)
        tl_mp.CreateNewNode(3, 1.0, 1.0, 0.0)
        tl_mp.CreateNewNode(4, 0.0, 1.0, 0.0)

        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.REACTION_X,tl_mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Y, KratosMultiphysics.REACTION_Y,tl_mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Z, KratosMultiphysics.REACTION_Z,tl_mp)

        # Create a submodelpart for boundary conditions
        tl_bcs = tl_mp.CreateSubModelPart("BoundaryCondtions")
        tl_bcs.AddNodes(bc_nodes)
        for node in tl_bcs.Nodes:
            node.Fix(KratosMultiphysics.DISPLACEMENT_X)
            node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X, 0.0)
            node.Fix(KratosMultiphysics.DISPLACEMENT_Y)
            node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y, 0.0)

        # Create Element and condition
        tl_elem = tl_mp.CreateNewElement("TotalLagrangianElement2D4N", 1, [1,2,3,4], tl_mp.GetProperties()[1])
        tl_load = tl_mp.CreateSubModelPart("LoadCondtions")
        tl_load.AddNodes(load_nodes)
        tl_cond = tl_mp.CreateNewCondition("LineLoadCondition2D2N", 1, load_nodes, tl_mp.GetProperties()[1])

        self._set_buffer(tl_mp)
        tl_lhs = KratosMultiphysics.CompressedMatrix()

        current_model = KratosMultiphysics.Model()
        ul_mp = current_model.CreateModelPart("ul_solid_part")
        self._add_variables(ul_mp)
        self._apply_material_properties(ul_mp, dim, False)

        # Create nodes
        ul_mp.CreateNewNode(1, 0.0, 0.0, 0.0)
        ul_mp.CreateNewNode(2, 1.0, 0.0, 0.0)
        ul_mp.CreateNewNode(3, 1.0, 1.0, 0.0)
        ul_mp.CreateNewNode(4, 0.0, 1.0, 0.0)

        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.REACTION_X,ul_mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Y, KratosMultiphysics.REACTION_Y,ul_mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Z, KratosMultiphysics.REACTION_Z,ul_mp)

        # Create a submodelpart for boundary conditions
        ul_bcs = ul_mp.CreateSubModelPart("BoundaryCondtions")
        ul_bcs.AddNodes(bc_nodes)
        for node in ul_bcs.Nodes:
            node.Fix(KratosMultiphysics.DISPLACEMENT_X)
            node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X, 0.0)
            node.Fix(KratosMultiphysics.DISPLACEMENT_Y)
            node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y, 0.0)

        # Create Element
        ul_elem = ul_mp.CreateNewElement("UpdatedLagrangianElement2D4N", 1, [1,2,3,4], ul_mp.GetProperties()[1])
        ul_load = ul_mp.CreateSubModelPart("LoadCondtions")
        ul_load.AddNodes(load_nodes)
        ul_cond = ul_mp.CreateNewCondition("LineLoadCondition2D2N", 1, load_nodes, tl_mp.GetProperties()[1])

        self._set_buffer(ul_mp)
        ul_lhs = KratosMultiphysics.CompressedMatrix()

        # Now we solve
        load = KratosMultiphysics.Vector(3)
        load[0] = 0.0
        load[1] = 0.0
        load[2] = 0.0

        delta_time = ul_mp.ProcessInfo[KratosMultiphysics.DELTA_TIME]
        time = ul_mp.ProcessInfo[KratosMultiphysics.TIME]

        tl_strategy = self._create_strategy(tl_mp)
        ul_strategy = self._create_strategy(ul_mp)

        for iter in range(1, 4):

            time += iter * delta_time
            tl_mp.CloneTimeStep(time)
            ul_mp.CloneTimeStep(time)

            #load[1] = iter * 1.0e10
            #tl_cond.SetValue(StructuralMechanicsApplication.LINE_LOAD, load)
            #ul_cond.SetValue(StructuralMechanicsApplication.LINE_LOAD, load)

            for node in tl_load.Nodes:
                node.Fix(KratosMultiphysics.DISPLACEMENT_X)
                node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X, iter * 5.0e-1)
            for node in ul_load.Nodes:
                node.Fix(KratosMultiphysics.DISPLACEMENT_X)
                node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X, iter * 5.0e-1)
            #for node in tl_load.Nodes:
                #node.Fix(KratosMultiphysics.DISPLACEMENT_Y)
                #node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y, iter * 5.0e-1)
            #for node in ul_load.Nodes:
                #node.Fix(KratosMultiphysics.DISPLACEMENT_Y)
                #node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y, iter * 5.0e-1)

            self._solve_with_strategy(tl_strategy, tl_lhs, iter)
            self._solve_with_strategy(ul_strategy, ul_lhs, iter)

            # Check displacement
            for i in range(2, 4):
                tl_dx = tl_mp.Nodes[i].GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X)
                tl_dy = tl_mp.Nodes[i].GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y)
                ul_dx = ul_mp.Nodes[i].GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X)
                ul_dy = ul_mp.Nodes[i].GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y)

                if (ul_dx > 0.0):
                    self.assertLess((tl_dx - ul_dx) / ul_dx, 1.0e-3)
                if (ul_dy > 0.0):
                    self.assertLess((tl_dy - ul_dy) / ul_dy, 1.0e-3)

            # Compare matrices
            for i in range(ul_lhs.Size1()):
                for j in range(ul_lhs.Size2()):
                    self.assertLess((ul_lhs[i, j] - tl_lhs[i, j]) / tl_lhs[i, j], 1.0e-3)

        #self.__post_process(tl_mp)
        #self.__post_process(ul_mp)

    def test_TL_2D_triangle(self):
        dim = 2
        current_model = KratosMultiphysics.Model()
        mp = current_model.CreateModelPart("solid_part")
        self._add_variables(mp)
        self._apply_material_properties(mp,dim)

        #create nodes
        mp.CreateNewNode(1,0.5,0.5,0.0)
        mp.CreateNewNode(2,0.7,0.2,0.0)
        mp.CreateNewNode(3,0.9,0.8,0.0)
        mp.CreateNewNode(4,0.3,0.7,0.0)
        mp.CreateNewNode(5,0.6,0.6,0.0)

        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.REACTION_X,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Y, KratosMultiphysics.REACTION_Y,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Z, KratosMultiphysics.REACTION_Z,mp)

        #create a submodelpart for boundary conditions
        bcs = mp.CreateSubModelPart("BoundaryCondtions")
        bcs.AddNodes([1,2,3,4])

        #create Element
        mp.CreateNewElement("TotalLagrangianElement2D3N", 1, [1,2,5], mp.GetProperties()[1])
        mp.CreateNewElement("TotalLagrangianElement2D3N", 2, [2,3,5], mp.GetProperties()[1])
        mp.CreateNewElement("TotalLagrangianElement2D3N", 3, [3,4,5], mp.GetProperties()[1])
        mp.CreateNewElement("TotalLagrangianElement2D3N", 4, [4,1,5], mp.GetProperties()[1])

        A,b = self._define_movement(dim)

        self._apply_BCs(bcs,A,b)
        self._solve(mp)
        self._check_results(mp,A,b)
        self._check_outputs(mp,A,dim)

        #self.__post_process(mp)

    def test_TL_2D_quadrilateral(self):
        dim = 2
        current_model = KratosMultiphysics.Model()
        mp = current_model.CreateModelPart("solid_part")
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

        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.REACTION_X,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Y, KratosMultiphysics.REACTION_Y,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Z, KratosMultiphysics.REACTION_Z,mp)

        #create a submodelpart for boundary conditions
        bcs = mp.CreateSubModelPart("BoundaryCondtions")
        bcs.AddNodes([1,5,6,8])

        #create Element
        mp.CreateNewElement("TotalLagrangianElement2D4N", 1, [8,7,3,5], mp.GetProperties()[1])
        mp.CreateNewElement("TotalLagrangianElement2D4N", 2, [6,4,7,8], mp.GetProperties()[1])
        mp.CreateNewElement("TotalLagrangianElement2D4N", 3, [1,2,4,6], mp.GetProperties()[1])
        mp.CreateNewElement("TotalLagrangianElement2D4N", 4, [4,2,3,7], mp.GetProperties()[1])
        mp.CreateNewElement("TotalLagrangianElement2D4N", 5, [2,1,5,3], mp.GetProperties()[1])

        A,b = self._define_movement(dim)

        self._apply_BCs(bcs,A,b)
        self._solve(mp)
        self._check_results(mp,A,b)
        self._check_outputs(mp,A,dim)

        #self.__post_process(mp)

    def test_TL_3D_tetra(self):
        dim = 3
        current_model = KratosMultiphysics.Model()
        mp = current_model.CreateModelPart("solid_part")
        self._add_variables(mp)
        self._apply_material_properties(mp,dim)

        #create nodes
        mp.CreateNewNode(1,0.0, 1.0, 0.0)
        mp.CreateNewNode(2,0.0, 1.0, 0.1)
        mp.CreateNewNode(3, 0.28739360416666665, 0.27808503701741405, 0.05672979583333333)
        mp.CreateNewNode(4, 0.0, 0.1, 0.0)
        mp.CreateNewNode(5, 0.1, 0.1, 0.1)
        mp.CreateNewNode(6, 1.0, 0.0, 0.0)
        mp.CreateNewNode(7, 1.2, 0.0, 0.1)

        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.REACTION_X,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Y, KratosMultiphysics.REACTION_Y,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Z, KratosMultiphysics.REACTION_Z,mp)

        #create a submodelpart for boundary conditions
        bcs = mp.CreateSubModelPart("BoundaryCondtions")
        bcs.AddNodes([1,2,4,5,6,7])

        #create Element
        mp.CreateNewElement("TotalLagrangianElement3D4N", 1,[5,3,1,2], mp.GetProperties()[1])
        mp.CreateNewElement("TotalLagrangianElement3D4N", 2,[3,1,2,6], mp.GetProperties()[1])
        mp.CreateNewElement("TotalLagrangianElement3D4N", 3,[6,4,7,3], mp.GetProperties()[1])
        mp.CreateNewElement("TotalLagrangianElement3D4N", 4,[5,4,1,3], mp.GetProperties()[1])
        mp.CreateNewElement("TotalLagrangianElement3D4N", 5,[4,1,3,6], mp.GetProperties()[1])
        mp.CreateNewElement("TotalLagrangianElement3D4N", 6,[5,4,3,7], mp.GetProperties()[1])
        mp.CreateNewElement("TotalLagrangianElement3D4N", 7,[3,5,7,2], mp.GetProperties()[1])
        mp.CreateNewElement("TotalLagrangianElement3D4N", 8,[6,7,2,3], mp.GetProperties()[1])

        A,b = self._define_movement(dim)

        self._apply_BCs(bcs,A,b)
        self._solve(mp)
        self._check_results(mp,A,b)
        self._check_outputs(mp,A,dim)

        #self.__post_process(mp)

    def test_TL_3D_prism(self):
        dim = 3
        current_model = KratosMultiphysics.Model()
        mp = current_model.CreateModelPart("solid_part")
        self._add_variables(mp)
        self._apply_material_properties(mp,dim)

        #create nodes
        mp.CreateNewNode(1,0.5,0.5,0.0)
        mp.CreateNewNode(2,0.7,0.2,0.0)
        mp.CreateNewNode(3,0.9,0.8,0.0)
        mp.CreateNewNode(4,0.3,0.7,0.0)
        mp.CreateNewNode(5,0.6,0.6,0.0)
        mp.CreateNewNode(6,0.5,0.5,0.1)
        mp.CreateNewNode(7,0.7,0.2,0.1)
        mp.CreateNewNode(8,0.9,0.8,0.1)
        mp.CreateNewNode(9,0.3,0.7,0.1)
        mp.CreateNewNode(10,0.6,0.6,0.1)
        mp.CreateNewNode(11,0.5,0.5,0.2)
        mp.CreateNewNode(12,0.7,0.2,0.2)
        mp.CreateNewNode(13,0.9,0.8,0.2)
        mp.CreateNewNode(14,0.3,0.7,0.2)
        mp.CreateNewNode(15,0.6,0.6,0.2)

        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.REACTION_X,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Y, KratosMultiphysics.REACTION_Y,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Z, KratosMultiphysics.REACTION_Z,mp)

        #create a submodelpart for boundary conditions
        bcs = mp.CreateSubModelPart("BoundaryCondtions")
        bcs.AddNodes([1,2,3,4,5,6,7,8,9,11,12,13,14,15])

        #create Element
        mp.CreateNewElement("TotalLagrangianElement3D6N", 1, [1,2,5,6,7,10], mp.GetProperties()[1])
        mp.CreateNewElement("TotalLagrangianElement3D6N", 2, [2,3,5,7,8,10], mp.GetProperties()[1])
        mp.CreateNewElement("TotalLagrangianElement3D6N", 3, [3,4,5,8,9,10], mp.GetProperties()[1])
        mp.CreateNewElement("TotalLagrangianElement3D6N", 4, [4,1,5,9,6,10], mp.GetProperties()[1])
        mp.CreateNewElement("TotalLagrangianElement3D6N", 5, [6,7,10,11,12,15], mp.GetProperties()[1])
        mp.CreateNewElement("TotalLagrangianElement3D6N", 6, [7,8,10,12,13,15], mp.GetProperties()[1])
        mp.CreateNewElement("TotalLagrangianElement3D6N", 7, [8,9,10,13,14,15], mp.GetProperties()[1])
        mp.CreateNewElement("TotalLagrangianElement3D6N", 8, [9,6,10,14,11,15], mp.GetProperties()[1])

        A,b = self._define_movement(dim)

        self._apply_BCs(bcs,A,b)
        self._solve(mp)
        self._check_results(mp,A,b)
        self._check_outputs(mp,A,dim)

        #self.__post_process(mp)

    def test_TL_3D_hexa(self):
        dim = 3
        current_model = KratosMultiphysics.Model()
        mp = current_model.CreateModelPart("solid_part")
        self._add_variables(mp)
        self._apply_material_properties(mp,dim)

        #create nodes
        mp.CreateNewNode(1, 0.00000,  1.00000,  1.00000)
        mp.CreateNewNode(2, 0.16500,  0.74500,  0.70200)
        mp.CreateNewNode(3, 0.27300,  0.75000,  0.23000)
        mp.CreateNewNode(4, 0.78800,  0.69300,  0.64400)
        mp.CreateNewNode(5, 0.32000,  0.18600,  0.64300)
        mp.CreateNewNode(6, 0.00000,  1.00000,  0.00000)
        mp.CreateNewNode(7, 0.00000,  0.00000,  1.00000)
        mp.CreateNewNode(8, 1.00000,  1.00000,  1.00000)
        mp.CreateNewNode(9, 0.67700,  0.30500,  0.68300)
        mp.CreateNewNode(10, 0.24900,  0.34200,  0.19200)
        mp.CreateNewNode(11, 0.85000,  0.64900,  0.26300)
        mp.CreateNewNode(12, 0.82600,  0.28800,  0.28800)
        mp.CreateNewNode(13, 0.00000,  0.00000,  0.00000)
        mp.CreateNewNode(14, 1.00000,  1.00000,  0.00000)
        mp.CreateNewNode(15, 1.00000,  0.00000,  1.00000)
        mp.CreateNewNode(16, 1.00000,  0.00000,  0.00000)

        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.REACTION_X,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Y, KratosMultiphysics.REACTION_Y,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Z, KratosMultiphysics.REACTION_Z,mp)

        #create a submodelpart for boundary conditions
        bcs = mp.CreateSubModelPart("BoundaryCondtions")
        bcs.AddNodes([1,6,7,8,13,14,15,16])


        #create Element
        mp.CreateNewElement("TotalLagrangianElement3D8N", 1,[10,5,2,3,13,7,1,6], mp.GetProperties()[1])
        mp.CreateNewElement("TotalLagrangianElement3D8N", 2,[12,9,5,10,16,15,7,13], mp.GetProperties()[1])
        mp.CreateNewElement("TotalLagrangianElement3D8N", 3,[12,11,3,10,9,4,2,5], mp.GetProperties()[1])
        mp.CreateNewElement("TotalLagrangianElement3D8N", 4,[9,4,2,5,15,8,1,7], mp.GetProperties()[1])
        mp.CreateNewElement("TotalLagrangianElement3D8N", 5,[4,11,3,2,8,14,6,1], mp.GetProperties()[1])
        mp.CreateNewElement("TotalLagrangianElement3D8N", 6,[11,4,9,12,14,8,15,16], mp.GetProperties()[1])
        mp.CreateNewElement("TotalLagrangianElement3D8N", 7,[11,12,10,3,14,16,13,6], mp.GetProperties()[1])

        A,b = self._define_movement(dim)

        self._apply_BCs(bcs,A,b)
        self._solve(mp)
        self._check_results(mp,A,b)
        self._check_outputs(mp,A,dim)

        #self.__post_process(mp)

    def test_UL_2D_triangle(self):
        dim = 2
        current_model = KratosMultiphysics.Model()
        mp = current_model.CreateModelPart("solid_part")
        self._add_variables(mp)
        self._apply_material_properties(mp,dim)

        #create nodes
        mp.CreateNewNode(1,0.5,0.5,0.0)
        mp.CreateNewNode(2,0.7,0.2,0.0)
        mp.CreateNewNode(3,0.9,0.8,0.0)
        mp.CreateNewNode(4,0.3,0.7,0.0)
        mp.CreateNewNode(5,0.6,0.6,0.0)

        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.REACTION_X,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Y, KratosMultiphysics.REACTION_Y,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Z, KratosMultiphysics.REACTION_Z,mp)

        #create a submodelpart for boundary conditions
        bcs = mp.CreateSubModelPart("BoundaryCondtions")
        bcs.AddNodes([1,2,3,4])

        #create Element
        mp.CreateNewElement("UpdatedLagrangianElement2D3N", 1, [1,2,5], mp.GetProperties()[1])
        mp.CreateNewElement("UpdatedLagrangianElement2D3N", 2, [2,3,5], mp.GetProperties()[1])
        mp.CreateNewElement("UpdatedLagrangianElement2D3N", 3, [3,4,5], mp.GetProperties()[1])
        mp.CreateNewElement("UpdatedLagrangianElement2D3N", 4, [4,1,5], mp.GetProperties()[1])

        A,b = self._define_movement(dim)

        self._set_buffer(mp)
        self._apply_BCs(bcs,A,b)
        self._solve(mp)
        self._check_results(mp,A,b)
        self._check_outputs(mp,A,dim)

        #self.__post_process(mp)

    def test_UL_2D_quadrilateral(self):
        dim = 2
        current_model = KratosMultiphysics.Model()
        mp = current_model.CreateModelPart("solid_part")
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

        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.REACTION_X,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Y, KratosMultiphysics.REACTION_Y,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Z, KratosMultiphysics.REACTION_Z,mp)

        #create a submodelpart for boundary conditions
        bcs = mp.CreateSubModelPart("BoundaryCondtions")
        bcs.AddNodes([1,5,6,8])

        #create Element
        mp.CreateNewElement("UpdatedLagrangianElement2D4N", 1, [8,7,3,5], mp.GetProperties()[1])
        mp.CreateNewElement("UpdatedLagrangianElement2D4N", 2, [6,4,7,8], mp.GetProperties()[1])
        mp.CreateNewElement("UpdatedLagrangianElement2D4N", 3, [1,2,4,6], mp.GetProperties()[1])
        mp.CreateNewElement("UpdatedLagrangianElement2D4N", 4, [4,2,3,7], mp.GetProperties()[1])
        mp.CreateNewElement("UpdatedLagrangianElement2D4N", 5, [2,1,5,3], mp.GetProperties()[1])

        A,b = self._define_movement(dim)

        self._set_buffer(mp)
        self._apply_BCs(bcs,A,b)
        self._solve(mp)
        self._check_results(mp,A,b)
        self._check_outputs(mp,A,dim)

        #self.__post_process(mp)

    def test_UL_3D_hexa(self):
        dim = 3
        current_model = KratosMultiphysics.Model()
        mp = current_model.CreateModelPart("solid_part")
        self._add_variables(mp)
        self._apply_material_properties(mp,dim)

        #create nodes
        mp.CreateNewNode(1, 0.00000,  1.00000,  1.00000)
        mp.CreateNewNode(2, 0.16500,  0.74500,  0.70200)
        mp.CreateNewNode(3, 0.27300,  0.75000,  0.23000)
        mp.CreateNewNode(4, 0.78800,  0.69300,  0.64400)
        mp.CreateNewNode(5, 0.32000,  0.18600,  0.64300)
        mp.CreateNewNode(6, 0.00000,  1.00000,  0.00000)
        mp.CreateNewNode(7, 0.00000,  0.00000,  1.00000)
        mp.CreateNewNode(8, 1.00000,  1.00000,  1.00000)
        mp.CreateNewNode(9, 0.67700,  0.30500,  0.68300)
        mp.CreateNewNode(10, 0.24900,  0.34200,  0.19200)
        mp.CreateNewNode(11, 0.85000,  0.64900,  0.26300)
        mp.CreateNewNode(12, 0.82600,  0.28800,  0.28800)
        mp.CreateNewNode(13, 0.00000,  0.00000,  0.00000)
        mp.CreateNewNode(14, 1.00000,  1.00000,  0.00000)
        mp.CreateNewNode(15, 1.00000,  0.00000,  1.00000)
        mp.CreateNewNode(16, 1.00000,  0.00000,  0.00000)

        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.REACTION_X,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Y, KratosMultiphysics.REACTION_Y,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Z, KratosMultiphysics.REACTION_Z,mp)

        #create a submodelpart for boundary conditions
        bcs = mp.CreateSubModelPart("BoundaryCondtions")
        bcs.AddNodes([1,6,7,8,13,14,15,16])


        #create Element
        mp.CreateNewElement("UpdatedLagrangianElement3D8N", 1,[10,5,2,3,13,7,1,6], mp.GetProperties()[1])
        mp.CreateNewElement("UpdatedLagrangianElement3D8N", 2,[12,9,5,10,16,15,7,13], mp.GetProperties()[1])
        mp.CreateNewElement("UpdatedLagrangianElement3D8N", 3,[12,11,3,10,9,4,2,5], mp.GetProperties()[1])
        mp.CreateNewElement("UpdatedLagrangianElement3D8N", 4,[9,4,2,5,15,8,1,7], mp.GetProperties()[1])
        mp.CreateNewElement("UpdatedLagrangianElement3D8N", 5,[4,11,3,2,8,14,6,1], mp.GetProperties()[1])
        mp.CreateNewElement("UpdatedLagrangianElement3D8N", 6,[11,4,9,12,14,8,15,16], mp.GetProperties()[1])
        mp.CreateNewElement("UpdatedLagrangianElement3D8N", 7,[11,12,10,3,14,16,13,6], mp.GetProperties()[1])

        A,b = self._define_movement(dim)

        self._set_buffer(mp)
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

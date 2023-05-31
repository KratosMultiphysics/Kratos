import KratosMultiphysics


import KratosMultiphysics.CableNetApplication as CableNetApplication
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import KratosMultiphysics.KratosUnittest as KratosUnittest

from KratosMultiphysics.CableNetApplication import empirical_spring_element_process

numpy_polyfit_available = False

import numpy as np

try:
    from numpy import polyfit
    numpy_polyfit_available = True
except: pass


import sys

class EmpiricalSpringTests(KratosUnittest.TestCase):
    def setUp(self):
        pass

    def _add_dofs(self,mp):
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.REACTION_X,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Y, KratosMultiphysics.REACTION_Y,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Z, KratosMultiphysics.REACTION_Z,mp)
    def _add_explicit_variables(self,mp):
        mp.AddNodalSolutionStepVariable(StructuralMechanicsApplication.MIDDLE_VELOCITY)
        mp.AddNodalSolutionStepVariable(StructuralMechanicsApplication.FRACTIONAL_ACCELERATION)
        mp.AddNodalSolutionStepVariable(StructuralMechanicsApplication.FRACTIONAL_ANGULAR_ACCELERATION)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_MASS)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.FORCE_RESIDUAL)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.RESIDUAL_VECTOR)
        mp.AddNodalSolutionStepVariable(StructuralMechanicsApplication.MIDDLE_ANGULAR_VELOCITY)
        mp.AddNodalSolutionStepVariable(StructuralMechanicsApplication.NODAL_INERTIA)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.MOMENT_RESIDUAL)
    def _add_variables(self,mp):
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)
        mp.AddNodalSolutionStepVariable(StructuralMechanicsApplication.POINT_LOAD)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.VOLUME_ACCELERATION)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.ACCELERATION)
        self._add_explicit_variables(mp)

    def _add_constitutive_law(self,mp):
        mp.GetProperties()[0].SetValue(KratosMultiphysics.CONSTITUTIVE_LAW,StructuralMechanicsApplication.TrussConstitutiveLaw())
    def _apply_material_properties(self,mp,dim,is_process_test):
        #define properties
        mp.GetProperties()[0].SetValue(KratosMultiphysics.YOUNG_MODULUS,206900000000.0)
        mp.GetProperties()[0].SetValue(KratosMultiphysics.DENSITY,7850)
        mp.GetProperties()[0].SetValue(StructuralMechanicsApplication.CROSS_AREA,40.0)
        mp.GetProperties()[0].SetValue(StructuralMechanicsApplication.RAYLEIGH_ALPHA,0.1)
        mp.GetProperties()[0].SetValue(StructuralMechanicsApplication.RAYLEIGH_BETA,0.001)
        if not is_process_test:
            mp.GetProperties()[0].SetValue(CableNetApplication.SPRING_DEFORMATION_EMPIRICAL_POLYNOMIAL,[-9.52380952e+01  ,1.16190476e+03 ,-5.34761905e+03  ,1.23952381e+04 ,1.30156287e-11])
    def _apply_BCs(self,mp,which_dof):
        if (which_dof == 'xyz'):
            KratosMultiphysics.VariableUtils().ApplyFixity(KratosMultiphysics.DISPLACEMENT_X, True, mp.Nodes)
            KratosMultiphysics.VariableUtils().ApplyFixity(KratosMultiphysics.DISPLACEMENT_Y, True, mp.Nodes)
            KratosMultiphysics.VariableUtils().ApplyFixity(KratosMultiphysics.DISPLACEMENT_Z, True, mp.Nodes)
        if (which_dof == 'yz'):
            KratosMultiphysics.VariableUtils().ApplyFixity(KratosMultiphysics.DISPLACEMENT_Y, True, mp.Nodes)
            KratosMultiphysics.VariableUtils().ApplyFixity(KratosMultiphysics.DISPLACEMENT_Z, True, mp.Nodes)
    def _apply_Neumann_BCs(self,mp,which_dof,load_size_dir):
        if(which_dof == 'x'):
            KratosMultiphysics.VariableUtils().SetScalarVar(StructuralMechanicsApplication.
                POINT_LOAD_X, load_size_dir, mp.Nodes)
        else:
            sys.exit('cannot apply given neumann boundary condition')

    def _create_nonlinear_static_strategy(self,mp):
        linear_solver = KratosMultiphysics.SkylineLUFactorizationSolver()
        builder_and_solver = KratosMultiphysics.ResidualBasedBlockBuilderAndSolver(linear_solver)
        scheme = KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme()
        convergence_criterion_1 = KratosMultiphysics.ResidualCriteria(1e-8,1e-8)
        convergence_criterion_2 = KratosMultiphysics.DisplacementCriteria(1e-8,1e-8)
        convergence_criterion   = KratosMultiphysics.AndCriteria(convergence_criterion_1,convergence_criterion_2)
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
        strategy.Initialize()
        strategy.Check()
        return strategy
    def _create_dynamic_implicit_strategy(self,mp):
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

        strategy.Initialize()
        strategy.Check()

        return strategy
    def _create_dynamic_explicit_strategy(self,mp):
        scheme = StructuralMechanicsApplication.ExplicitCentralDifferencesScheme(0.00,0.00,0.00)
        strategy = StructuralMechanicsApplication.MechanicalExplicitStrategy(mp,scheme,0,0,1)
        strategy.SetEchoLevel(0)
        return strategy
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

    def _evaluate_function(self,x,func):
        res = 0.0
        for i,num in enumerate(func):
            res += num*np.power(x,len(func)-i-1)
        return res
    def _set_up_standard_test(self,is_process_test):
        dim = 3
        current_model = KratosMultiphysics.Model()
        mp = current_model.CreateModelPart("solid_part")
        mp.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, dim)
        self._add_variables(mp)
        self._apply_material_properties(mp,dim,is_process_test)
        self._add_constitutive_law(mp)
        #create nodes
        mp.CreateNewNode(1,1.2,0.0,0.0)
        mp.CreateNewNode(2,0.0,0.0,0.0)
        #add dofs
        self._add_dofs(mp)
        #create condition
        mp.CreateNewCondition("PointLoadCondition3D1N",1,[1],mp.GetProperties()[0])
        #create submodelparts for dirichlet boundary conditions
        bcs_xyz = mp.CreateSubModelPart("Dirichlet_XYZ")
        bcs_xyz.AddNodes([2])
        bcs_xz = mp.CreateSubModelPart("Dirichlet_YZ")
        bcs_xz.AddNodes([1])
        #create a submodalpart for neumann boundary conditions
        bcs_neumann = mp.CreateSubModelPart("PointLoad3D_neumann")
        bcs_neumann.AddNodes([1])
        bcs_neumann.AddConditions([1])
        if not is_process_test:
            #create Element
            mp.CreateNewElement("EmpiricalSpringElement3D2N", 1, [2,1], mp.GetProperties()[0])
        #apply boundary conditions
        self._apply_BCs(bcs_xyz,'xyz')
        self._apply_BCs(bcs_xz,'xz')

        if not is_process_test:
            return mp,bcs_neumann
        return current_model,mp,bcs_neumann

    def test_quasi_static(self):
        mp,bcs_neumann = self._set_up_standard_test(is_process_test=False)

        #apply boundary conditions
        force_x = 100.0
        #createfunction
        f = [-9.52380952e+01  ,1.16190476e+03 ,-5.34761905e+03  ,1.23952381e+04 ,1.30156287e-11]

        #loop over time
        time_end = 140.0
        time_delta = 10.0
        time_i = 0.0
        self._set_and_fill_buffer(mp,2,time_delta)

        strategy = self._create_nonlinear_static_strategy(mp)
        while (time_i < time_end):
            time_i += time_delta
            mp.CloneTimeStep(time_i)

            self._apply_Neumann_BCs(bcs_neumann,'x',force_x*time_i)
            #solve + compare
            strategy.Solve()

            disp_solution = bcs_neumann.Nodes[1].GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X)
            force_solution = self._evaluate_function(disp_solution,f)

            self.assertAlmostEqual(force_solution, force_x*time_i)
    def test_quasi_static_process(self):
        if not numpy_polyfit_available: self.skipTest("python package \"numpy.polyfit\" not available")
        model,mp,bcs_neumann = self._set_up_standard_test(is_process_test=True)

        #apply boundary conditions
        force_x = 100.0
        #createfunction
        f = [-9.52380952e+01  ,1.16190476e+03 ,-5.34761905e+03  ,1.23952381e+04 ,1.30156287e-11]

        #loop over time
        time_end = 140.0
        time_delta = 10.0
        time_i = 0.0
        self._set_and_fill_buffer(mp,2,time_delta)

        custom_settings = KratosMultiphysics.Parameters("""
        {
            "model_part_name"           : "solid_part.element_mp",
            "computing_model_part_name" : "solid_part",
            "node_ids"                  : [1,2],
            "element_id"                : 1,
            "property_id"               : 0,
            "displacement_data"         : [0.0,0.5,1.5,2.5,5.0],
            "force_data"                : [0.0,5000.0,10000.0,12000.0,14000.0],
            "polynomial_order"          : 4
        }
        """)

        element_mp = mp.CreateSubModelPart("element_mp")
        element_mp.AddNodes([1,2])

        test_process = empirical_spring_element_process.EmpiricalSpringElementProcess(model,custom_settings)
        test_process.ExecuteInitialize()


        strategy = self._create_nonlinear_static_strategy(mp)
        while (time_i < time_end):
            time_i += time_delta
            mp.CloneTimeStep(time_i)

            self._apply_Neumann_BCs(bcs_neumann,'x',force_x*time_i)
            #solve + compare
            strategy.Solve()

            disp_solution = bcs_neumann.Nodes[1].GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X)
            force_solution = self._evaluate_function(disp_solution,f)
            self.assertAlmostEqual(force_solution, force_x*time_i, 2)
    def test_implicit_dynamic(self):
        mp,bcs_neumann = self._set_up_standard_test(is_process_test=False)

        #apply boundary conditions
        force_x = 14000.0
        displacement_expected = 4.0

        #loop over time
        time_end = 140.0
        time_delta = 6.0
        time_i = 0.0
        self._set_and_fill_buffer(mp,2,time_delta)

        strategy = self._create_dynamic_implicit_strategy(mp)
        while (time_i < time_end):
            time_i += time_delta
            mp.CloneTimeStep(time_i)

            self._apply_Neumann_BCs(bcs_neumann,'x',force_x)
            #solve + compare
            strategy.Solve()

        disp_solution = bcs_neumann.Nodes[1].GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X)
        self.assertAlmostEqual(disp_solution, displacement_expected, 1)
    def test_explicit_dynamic(self):
        mp,bcs_neumann = self._set_up_standard_test(is_process_test=False)


        #apply boundary conditions
        force_x = 14000.0
        displacement_expected = 4.0

        #loop over time
        time_end = 140.0
        time_delta = 6.0
        time_i = 0.0
        self._set_and_fill_buffer(mp,2,time_delta)

        strategy = self._create_dynamic_explicit_strategy(mp)
        while (time_i < time_end):
            time_i += time_delta
            mp.CloneTimeStep(time_i)

            self._apply_Neumann_BCs(bcs_neumann,'x',force_x)
            #solve + compare
            strategy.Solve()

        disp_solution = bcs_neumann.Nodes[1].GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X)
        self.assertAlmostEqual(disp_solution, displacement_expected, 1)


if __name__ == '__main__':
    KratosUnittest.main()

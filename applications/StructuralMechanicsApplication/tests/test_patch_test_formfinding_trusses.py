from __future__ import print_function, absolute_import, division
import KratosMultiphysics

import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import KratosMultiphysics.KratosUnittest as KratosUnittest



class TestPatchTestFormfinding(KratosUnittest.TestCase):
    def setUp(self):
        pass

    def _add_dofs(self,mp):
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.REACTION_X,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Y, KratosMultiphysics.REACTION_Y,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Z, KratosMultiphysics.REACTION_Z,mp)

    def _add_variables(self,mp):
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)
        mp.AddNodalSolutionStepVariable(StructuralMechanicsApplication.POINT_LOAD)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.VOLUME_ACCELERATION)


    def _add_constitutive_law(self,mp,elastic_flag):
        cl = StructuralMechanicsApplication.TrussPlasticityConstitutiveLaw()
        if elastic_flag:
            cl = StructuralMechanicsApplication.TrussConstitutiveLaw()
        mp.GetProperties()[0].SetValue(KratosMultiphysics.CONSTITUTIVE_LAW,cl)


    def _apply_material_properties(self,mp,dim,pre_stress,area):
        #define properties
        mp.GetProperties()[0].SetValue(KratosMultiphysics.YOUNG_MODULUS,210e9)
        mp.GetProperties()[0].SetValue(KratosMultiphysics.DENSITY,7850)
        mp.GetProperties()[0].SetValue(StructuralMechanicsApplication.CROSS_AREA,area)
        mp.GetProperties()[0].SetValue(StructuralMechanicsApplication.TRUSS_PRESTRESS_PK2,pre_stress)
        g = [0,0,0]
        mp.GetProperties()[0].SetValue(KratosMultiphysics.VOLUME_ACCELERATION,g)

    def _apply_material_properties_formfinding(self,mp,dim):
        #define properties
        mp.GetProperties()[0].SetValue(KratosMultiphysics.YOUNG_MODULUS,0.0)
        mp.GetProperties()[0].SetValue(KratosMultiphysics.DENSITY,0.0)
        mp.GetProperties()[0].SetValue(StructuralMechanicsApplication.CROSS_AREA,1.0)
        mp.GetProperties()[0].SetValue(StructuralMechanicsApplication.TRUSS_PRESTRESS_PK2,1.0)
        g = [0,0,0]
        mp.GetProperties()[0].SetValue(KratosMultiphysics.VOLUME_ACCELERATION,g)


    def _apply_BCs(self,mp,which_dof):
        KratosMultiphysics.VariableUtils().ApplyFixity(KratosMultiphysics.DISPLACEMENT_X, True, mp.Nodes)
        KratosMultiphysics.VariableUtils().ApplyFixity(KratosMultiphysics.DISPLACEMENT_Y, True, mp.Nodes)
        KratosMultiphysics.VariableUtils().ApplyFixity(KratosMultiphysics.DISPLACEMENT_Z, True, mp.Nodes)


    def _solve_nonlinear(self,mp):
        linear_solver = KratosMultiphysics.SkylineLUFactorizationSolver()
        builder_and_solver = KratosMultiphysics.ResidualBasedBlockBuilderAndSolver(linear_solver)
        scheme = KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme()
        convergence_criterion = KratosMultiphysics.ResidualCriteria(1e-12,1e-12)
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
        strategy.Solve()

    def _setup_formfinding(self,mp):
        linear_solver = KratosMultiphysics.SkylineLUFactorizationSolver()
        builder_and_solver = KratosMultiphysics.ResidualBasedBlockBuilderAndSolver(linear_solver)
        scheme = KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme()
        convergence_criterion = KratosMultiphysics.ResidualCriteria(1e-12,1e-12)
        convergence_criterion.SetEchoLevel(0)

        max_iters = 1000
        compute_reactions = True
        reform_step_dofs = True
        move_mesh_flag = True
        proj_settings = KratosMultiphysics.Parameters("""{
                "model_part_name"  : "Structure",
                "echo_level"       : 0,
                "projection_type"  : "planar",
                "global_direction" : [1,0,0],
                "variable_name"    : "LOCAL_MATERIAL_AXIS_1",
                "visualize_in_vtk" : false,
                "method_specific_settings" : { },
                "check_local_space_dimension" : false
        }""")
        strategy = StructuralMechanicsApplication.FormfindingStrategy(mp,
                                                        scheme,
                                                        linear_solver,
                                                        convergence_criterion,
                                                        builder_and_solver,
                                                        mp,
                                                        False,
                                                        "none",
                                                        proj_settings,
                                                        max_iters,
                                                        compute_reactions,
                                                        reform_step_dofs,
                                                        move_mesh_flag)
        strategy.SetEchoLevel(0)
        strategy.Initialize()
        return strategy

    def _do_formfinding(self,strategy):
        strategy.Solve()



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

    def test_formfinding_trusses(self):
        dim = 3
        current_model = KratosMultiphysics.Model()
        mp = current_model.CreateModelPart("solid_part")
        mp.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, dim)
        self._add_variables(mp)
        self._apply_material_properties_formfinding(mp,dim)
        self._add_constitutive_law(mp,True)

        #create nodes
        mp.CreateNewNode(1,0.0,0.0,0.0)
        mp.CreateNewNode(2,2.0,0.0,0.0)
        mp.CreateNewNode(3,2.0,1.0,0.0)
        mp.CreateNewNode(4,0.0,1.0,0.0)
        mp.CreateNewNode(5,1.0,0.0,0.0)
        #add dofs
        self._add_dofs(mp)

        #create submodelparts for dirichlet boundary conditions
        bcs_xyz = mp.CreateSubModelPart("Dirichlet_XYZ")
        bcs_xyz.AddNodes([1,2,3,4])

        #create Element
        mp.CreateNewElement("TrussElement3D2N", 1, [1,5], mp.GetProperties()[0])
        mp.CreateNewElement("TrussElement3D2N", 2, [2,5], mp.GetProperties()[0])
        mp.CreateNewElement("TrussElement3D2N", 3, [3,5], mp.GetProperties()[0])
        mp.CreateNewElement("TrussElement3D2N", 4, [4,5], mp.GetProperties()[0])


        self._apply_BCs(bcs_xyz,'xyz')

        #find form suitable for given pre-stress state
        formfinding_strategy = self._setup_formfinding(mp)
        self._do_formfinding(formfinding_strategy)

        final_shape = [mp.Nodes[5].X0,mp.Nodes[5].Y0,mp.Nodes[5].Z0]
        self.assertAlmostEqual(final_shape[0],1.0)
        self.assertAlmostEqual(final_shape[1],0.5)
        self.assertAlmostEqual(final_shape[2],0.0)


        ## now do the analysis on the form-found shape
        pre_stress,area = 1.0e+11,1.0e-2

        self._apply_material_properties(mp,dim,pre_stress,area)
        self._solve_nonlinear(mp)

        self.assertAlmostEqual(mp.Nodes[5].GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X),0.0)
        self.assertAlmostEqual(mp.Nodes[5].GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y),0.0)
        self.assertAlmostEqual(mp.Nodes[5].GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Z),0.0)

        for element in mp.Elements:
            out = element.CalculateOnIntegrationPoints(KratosMultiphysics.FORCE,mp.ProcessInfo)
            self.assertAlmostEqual(out[0][0],area*pre_stress,0)
            self.assertAlmostEqual(out[0][1],0.0)
            self.assertAlmostEqual(out[0][2],0.0)






if __name__ == '__main__':
    KratosUnittest.main()


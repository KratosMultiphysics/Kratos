from __future__ import print_function, absolute_import, division
import KratosMultiphysics

import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import KratosMultiphysics.KratosUnittest as KratosUnittest


class TestPatchTestMembrane(KratosUnittest.TestCase):
    def setUp(self):
        pass
    
    def _add_variables(self,mp):
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)

    def _add_dofs(self,mp):
        # Adding the dofs AND their corresponding reaction!
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.REACTION_X,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Y, KratosMultiphysics.REACTION_Y,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Z, KratosMultiphysics.REACTION_Z,mp)

    def _create_nodes(self,mp):
        mp.CreateNewNode(1, 0.0, 0.0, 0.0)
        mp.CreateNewNode(2, 0.25, 0.0, 0.0)
        mp.CreateNewNode(3, 0.50, 0.0, 0.0)
        mp.CreateNewNode(4, 0.0, 0.25, 0.01)
        mp.CreateNewNode(5, 0.25, 0.25, 0.01)
        mp.CreateNewNode(6, 0.50, 0.25, 0.01)
        mp.CreateNewNode(7, 0.0, 0.50, 0.02)
        mp.CreateNewNode(8, 0.25, 0.50, 0.02)
        mp.CreateNewNode(9, 0.50, 0.50, 0.02)
        mp.CreateNewNode(10, 0.0, 0.75, 0.01)
        mp.CreateNewNode(11, 0.25, 0.75, 0.01)
        mp.CreateNewNode(12, 0.50, 0.75, 0.01)
        mp.CreateNewNode(13, 0.0, 1.0, 0.0)
        mp.CreateNewNode(14, 0.25, 1.0, 0.0)
        mp.CreateNewNode(15, 0.50, 1.0, 0.0)

    def _create_elements_3d3n(self,mp):
        element_name = "PreStressMembraneElement3D3N"
        mp.CreateNewElement(element_name, 1, [1,2,5], mp.GetProperties()[1])
        mp.CreateNewElement(element_name, 2, [1,5,4], mp.GetProperties()[1])
        mp.CreateNewElement(element_name, 3, [2,3,6], mp.GetProperties()[1])
        mp.CreateNewElement(element_name, 4, [2,6,5], mp.GetProperties()[1])

        mp.CreateNewElement(element_name, 5, [4,5,8], mp.GetProperties()[1])
        mp.CreateNewElement(element_name, 6, [4,8,7], mp.GetProperties()[1])
        mp.CreateNewElement(element_name, 7, [5,6,9], mp.GetProperties()[1])
        mp.CreateNewElement(element_name, 8, [5,9,8], mp.GetProperties()[1])

        mp.CreateNewElement(element_name, 9, [7,8,11], mp.GetProperties()[1])
        mp.CreateNewElement(element_name, 10, [7,11,10], mp.GetProperties()[1])
        mp.CreateNewElement(element_name, 11, [8,9,12], mp.GetProperties()[1])
        mp.CreateNewElement(element_name, 12, [8,12,11], mp.GetProperties()[1])

        mp.CreateNewElement(element_name, 13, [10,11,14], mp.GetProperties()[1])
        mp.CreateNewElement(element_name, 14, [10,14,13], mp.GetProperties()[1])
        mp.CreateNewElement(element_name, 15, [11,12,15], mp.GetProperties()[1])
        mp.CreateNewElement(element_name, 16, [11,15,14], mp.GetProperties()[1])

    def _create_elements_3d4n(self,mp):
        element_name = "PreStressMembraneElement3D4N"
        mp.CreateNewElement(element_name, 1, [1,2,5,4], mp.GetProperties()[1])
        mp.CreateNewElement(element_name, 2, [2,3,6,5], mp.GetProperties()[1])

        mp.CreateNewElement(element_name, 3, [4,5,8,7], mp.GetProperties()[1])
        mp.CreateNewElement(element_name, 4, [5,6,9,8], mp.GetProperties()[1])

        mp.CreateNewElement(element_name, 5, [7,8,11,10], mp.GetProperties()[1])
        mp.CreateNewElement(element_name, 6, [8,9,12,11], mp.GetProperties()[1])

        mp.CreateNewElement(element_name, 7, [10,11,14,13], mp.GetProperties()[1])
        mp.CreateNewElement(element_name, 8, [11,12,15,14], mp.GetProperties()[1])

    def _apply_dirichlet_BCs(self,mp):
        KratosMultiphysics.VariableUtils().ApplyFixity(KratosMultiphysics.DISPLACEMENT_X, True, mp.Nodes)
        KratosMultiphysics.VariableUtils().ApplyFixity(KratosMultiphysics.DISPLACEMENT_Y, True, mp.Nodes)
        KratosMultiphysics.VariableUtils().ApplyFixity(KratosMultiphysics.DISPLACEMENT_Z, True, mp.Nodes)


    def _apply_material_properties(self,mp):
        #define properties
        mp.GetProperties()[1].SetValue(KratosMultiphysics.YOUNG_MODULUS,1000.0)
        mp.GetProperties()[1].SetValue(KratosMultiphysics.POISSON_RATIO,0.20)
        mp.GetProperties()[1].SetValue(KratosMultiphysics.THICKNESS,0.001)
        mp.GetProperties()[1].SetValue(KratosMultiphysics.DENSITY,700.0)
        prestress = KratosMultiphysics.Vector(3)
        prestress[0]=0.0
        prestress[1]=20.0
        prestress[2]=0.0
        mp.GetProperties()[1].SetValue(StructuralMechanicsApplication.PRESTRESS_VECTOR,prestress)

        gravity = [0,0,-9.81]
        mp.GetProperties()[1].SetValue(KratosMultiphysics.VOLUME_ACCELERATION,gravity)

        constitutive_law = StructuralMechanicsApplication.LinearElasticPlaneStress2DLaw()

        mp.GetProperties()[1].SetValue(KratosMultiphysics.CONSTITUTIVE_LAW,constitutive_law)
        mp.GetProperties()[1].SetValue(StructuralMechanicsApplication.PROJECTION_TYPE_COMBO,"planar")

    def _solve(self,mp):
        linear_solver = KratosMultiphysics.SkylineLUFactorizationSolver()
        builder_and_solver = KratosMultiphysics.ResidualBasedBlockBuilderAndSolver(linear_solver)
        scheme = KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme()
        convergence_criterion = StructuralMechanicsApplication.DisplacementAndOtherDoFCriteria(1e-15,1e-15)
        convergence_criterion.SetEchoLevel(0)

        max_iters = 1000
        compute_reactions = True
        reform_step_dofs = False
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



    def _check_results(self,node,displacement_results):
        #check that the results are exact on the node
        displacement = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT)
        self.assertAlmostEqual(displacement[0], displacement_results[0], 4)
        self.assertAlmostEqual(displacement[1], displacement_results[1], 4)
        self.assertAlmostEqual(displacement[2], displacement_results[2], 4)


    
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

    '''
    def execute_membrane_test_3d3n(self, current_model, displacement_results, do_post_processing):
        mp = current_model.CreateModelPart("Structure")
        mp.SetBufferSize(2)

        self._add_variables(mp)
        self._apply_material_properties(mp)
        self._create_nodes(mp)
        self._add_dofs(mp)
        self._create_elements_3d3n(mp)

        #create a submodelpart for dirichlet boundary conditions
        bcs_dirichlet = mp.CreateSubModelPart("BoundaryCondtionsDirichlet")
        bcs_dirichlet.AddNodes([1,2,3,13,14,15])

        self._apply_dirichlet_BCs(bcs_dirichlet)
        self._solve(mp)

        self._check_results(mp.Nodes[8],displacement_results)

        if do_post_processing == True:
            self.__post_process(mp)
    '''

    def _set_up_base_system(self,current_model):
        mp = current_model.CreateModelPart("Structure")
        mp.SetBufferSize(2)

        self._add_variables(mp)
        self._apply_material_properties(mp)
        self._create_nodes(mp)
        self._add_dofs(mp)

        #create a submodelpart for dirichlet boundary conditions
        bcs_dirichlet = mp.CreateSubModelPart("BoundaryCondtionsDirichlet")
        bcs_dirichlet.AddNodes([1,2,3,13,14,15])
        self._apply_dirichlet_BCs(bcs_dirichlet)

        return mp

    def set_up_membrane_system_3d3n(self,current_model):
        mp = self._set_up_base_system(current_model)
        self._create_elements_3d3n(mp)

    def set_up_membrane_system_3d4n(self,current_model):
        mp = self._set_up_base_system(current_model)
        self._create_elements_3d4n(mp)

    def test_membrane_3d3n_static(self):
        displacement_results = [-0.3853903940829765 , -0.2299393888361787 , -2.213110569935068]

        current_model = KratosMultiphysics.Model()
        '''
        self.execute_membrane_test_3d3n(current_model,
                                displacement_results,
                                False) # Do PostProcessing for GiD?
        '''

        self.set_up_membrane_system_3d3n(current_model)

        self._solve(mp)

        self._check_results(mp.Nodes[8],displacement_results)

        if do_post_processing == True:
            self.__post_process(mp)

    def test_membrane_3d3n_dynamic(self):
        displacement_results = [-0.3853903940829765 , -0.2299393888361787 , -2.213110569935068]

        current_model = KratosMultiphysics.Model()
        '''
        self.execute_membrane_test_3d3n(current_model,
                                displacement_results,
                                False) # Do PostProcessing for GiD?
        '''
        self.set_up_membrane_system_3d3n(current_model)        
        #time integration parameters
        dt = 0.01
        time = 0.0
        end_time = 10.0
        step = 0

        self._set_and_fill_buffer(mp,2,dt)

        while(time <= end_time):
            time = time + dt
            step = step + 1
            mp.CloneTimeStep(time)
        
            self._solve(mp)

            self._check_results(mp.Nodes[8],displacement_results)

        if do_post_processing == True:
            self.__post_process(mp)


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
                                                "nodal_results"       : ["DISPLACEMENT", "REACTIONS"],
                                                "gauss_point_results" : []
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

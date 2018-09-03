from __future__ import print_function, absolute_import, division
import KratosMultiphysics

import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import KratosMultiphysics.KratosUnittest as KratosUnittest


class TestPatchTestFormfinding(KratosUnittest.TestCase):
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

    def _create_nodes(self,mp,element_name):
        mp.CreateNewNode(1, 0.0, 0.0, 0.0)
        mp.CreateNewNode(2, 0.0, 10.0, 5.0)
        mp.CreateNewNode(3, 10.0, 0.0, 5.0)
        mp.CreateNewNode(4, 10.0,  10.0,  0.0)
        mp.CreateNewNode(5, 6.0,  6.0,  5.0)

    def _create_elements(self,mp,element_name):
        mp.CreateNewElement(element_name, 1, [1,2,5], mp.GetProperties()[1])
        mp.CreateNewElement(element_name, 2, [1,5,3], mp.GetProperties()[1])
        mp.CreateNewElement(element_name, 3, [3,4,5], mp.GetProperties()[1])
        mp.CreateNewElement(element_name, 4, [2,4,5], mp.GetProperties()[1])


    def _apply_dirichlet_BCs(self,mp):
        KratosMultiphysics.VariableUtils().ApplyFixity(KratosMultiphysics.DISPLACEMENT_X, True, mp.Nodes)
        KratosMultiphysics.VariableUtils().ApplyFixity(KratosMultiphysics.DISPLACEMENT_Y, True, mp.Nodes)
        KratosMultiphysics.VariableUtils().ApplyFixity(KratosMultiphysics.DISPLACEMENT_Z, True, mp.Nodes)


    def _apply_material_properties(self,mp):
        #define properties
        mp.GetProperties()[1].SetValue(KratosMultiphysics.YOUNG_MODULUS,0.0)
        mp.GetProperties()[1].SetValue(KratosMultiphysics.POISSON_RATIO,0.0)
        mp.GetProperties()[1].SetValue(KratosMultiphysics.THICKNESS,1.0)
        mp.GetProperties()[1].SetValue(KratosMultiphysics.DENSITY,1.0)
        prestress = KratosMultiphysics.Vector(3)
        prestress[0]=2.0
        prestress[1]=1.0
        prestress[2]=0.0
        mp.GetProperties()[1].SetValue(StructuralMechanicsApplication.PRESTRESS_VECTOR,prestress)


        cl = StructuralMechanicsApplication.LinearElasticPlaneStress2DLaw()

        mp.GetProperties()[1].SetValue(KratosMultiphysics.CONSTITUTIVE_LAW,cl)
        mp.GetProperties()[1].SetValue(StructuralMechanicsApplication.PROJECTION_TYPE_COMBO,"planar")

    def _solve(self,mp):
        #define a minimal newton raphson solver
        linear_solver = KratosMultiphysics.SkylineLUFactorizationSolver()
        builder_and_solver = KratosMultiphysics.ResidualBasedBlockBuilderAndSolver(linear_solver)
        scheme = KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme()
        convergence_criterion = KratosMultiphysics.DisplacementCriteria(1e-3,1e-3)
        convergence_criterion.SetEchoLevel(0)

        max_iters = 50
        compute_reactions = False
        reform_step_dofs = False
        calculate_norm_dx = False
        move_mesh_flag = True
        print_after_formfinding = False
        line_search = True
        strategy = KratosMultiphysics.StructuralMechanicsApplication.FormfindingUpdatedReferenceStrategy(mp,
                                                                  scheme,
                                                                  linear_solver,
                                                                  convergence_criterion,
                                                                  max_iters,
                                                                  compute_reactions,
                                                                  reform_step_dofs,
                                                                  move_mesh_flag,
                                                                  print_after_formfinding,
                                                                  line_search)
        strategy.SetEchoLevel(0)
        strategy.Check()
        strategy.Solve()



    def _check_results(self,node,displacement_results):
        #check that the results are exact on the node
        disp = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT)
        self.assertAlmostEqual(disp[0], displacement_results[0], 4)
        self.assertAlmostEqual(disp[1], displacement_results[1], 4)
        self.assertAlmostEqual(disp[2], displacement_results[2], 4)




    def execute_formfinding_test(self, element_name, displacement_results, do_post_processing):
        mp = KratosMultiphysics.ModelPart("solid_part")
        mp.SetBufferSize(2)

        self._add_variables(mp)
        self._apply_material_properties(mp)
        self._create_nodes(mp,element_name)
        self._add_dofs(mp)
        self._create_elements(mp,element_name)

        #create a submodelpart for dirichlet boundary conditions
        bcs_dirichlet = mp.CreateSubModelPart("BoundaryCondtionsDirichlet")
        bcs_dirichlet.AddNodes([1,2,3,4])

        self._apply_dirichlet_BCs(bcs_dirichlet)
        self._solve(mp)


        self._check_results(mp.Nodes[5],displacement_results)

        if do_post_processing:
            self.__post_process(mp)


    def test_formfinding(self):
        element_name = "PreStressMembraneElement3D3N"
        displacement_results = [-0.3853903940829765 , -0.2299393888361787 , -2.213110569935068]

        self.execute_formfinding_test(element_name,
                                displacement_results,
                                False) # Do PostProcessing for GiD?




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

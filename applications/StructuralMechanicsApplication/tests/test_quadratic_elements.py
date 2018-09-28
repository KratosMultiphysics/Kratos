from __future__ import print_function, absolute_import, division
import KratosMultiphysics

import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import KratosMultiphysics.KratosUnittest as KratosUnittest


class TestQuadraticElements(KratosUnittest.TestCase):
    def setUp(self):
        pass

    def _add_variables(self,mp):
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.VOLUME_ACCELERATION)


    def _apply_BCs(self,mp, coeff = 1.0):
        for node in mp.Nodes:
            node.Fix(KratosMultiphysics.DISPLACEMENT_X)
            node.Fix(KratosMultiphysics.DISPLACEMENT_Y)
            node.Fix(KratosMultiphysics.DISPLACEMENT_Z)

        for node in mp.Nodes:
            u = KratosMultiphysics.Vector(3)
            u[0] = coeff * node.X0**2
            u[1] = coeff * node.Y0**2
            u[2] = coeff * node.Z0**2

            node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT,0,u)

    def _apply_material_properties(self,mp,dim):
        #define properties
        mp.GetProperties()[0].SetValue(KratosMultiphysics.YOUNG_MODULUS,210e9)
        mp.GetProperties()[0].SetValue(KratosMultiphysics.POISSON_RATIO,0.3)
        mp.GetProperties()[0].SetValue(KratosMultiphysics.THICKNESS,1.0)
        mp.GetProperties()[0].SetValue(KratosMultiphysics.DENSITY,1.0)

        g = [0,0,0]
        mp.GetProperties()[0].SetValue(KratosMultiphysics.VOLUME_ACCELERATION,g)

        if(dim == 2):
            cl = StructuralMechanicsApplication.LinearElasticPlaneStress2DLaw()
        else:
            cl = StructuralMechanicsApplication.LinearElastic3DLaw()
        mp.GetProperties()[0].SetValue(KratosMultiphysics.CONSTITUTIVE_LAW,cl)

    def _solve(self,mp):

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


        #strategy = KratosMultiphysics.ResidualBasedNewtonRaphsonStrategy(mp,
                                                                        #scheme,
                                                                        #linear_solver,
                                                                        #convergence_criterion,
                                                                        #builder_and_solver,
                                                                        #max_iters,
                                                                        #compute_reactions,
                                                                        #reform_step_dofs,
                                                                        #move_mesh_flag)
        strategy.SetEchoLevel(0)

        strategy.Check()
        strategy.Solve()

    def _check_outputs(self,mp,dim, coeff = 1.0):
        for elem in mp.Elements:
            strains = elem.CalculateOnIntegrationPoints(KratosMultiphysics.GREEN_LAGRANGE_STRAIN_VECTOR, mp.ProcessInfo)
            coords = elem.CalculateOnIntegrationPoints(KratosMultiphysics.INTEGRATION_COORDINATES, mp.ProcessInfo)
            for strain,coord in zip(strains, coords):
                for i in range(2):
                    self.assertLessEqual(abs((coord[i] - 0.5 * strain[i]/coeff)/coord[i]), 5.0e-2)

    def test_Quad8(self):
        dim = 2
        mp = KratosMultiphysics.ModelPart("solid_part")
        self._add_variables(mp)

        #KratosMultiphysics.ModelPartIO("quadratic_test/static_quadratic_quad_test").ReadModelPart(mp)

        self._apply_material_properties(mp,dim)

        # Create nodes
        mp.CreateNewNode(1,0.0000000000,2.0000000000,0.0000000000)
        mp.CreateNewNode(2,0.0000000000,1.7500000000,0.0000000000)
        mp.CreateNewNode(3,0.2500000000,1.9375000000,0.0000000000)
        mp.CreateNewNode(4,0.0000000000,1.5000000000,0.0000000000)
        mp.CreateNewNode(5,0.5000000000,1.8750000000,0.0000000000)
        mp.CreateNewNode(6,0.2500000000,1.4610595703,0.0000000000)
        mp.CreateNewNode(7,0.5000000000,1.6485595703,0.0000000000)
        mp.CreateNewNode(8,0.0000000000,1.2500000000,0.0000000000)
        mp.CreateNewNode(9,0.5000000000,1.4287109375,0.0000000000)
        mp.CreateNewNode(10,0.7500000000,1.8125000000,0.0000000000)
        mp.CreateNewNode(11,0.5000000000,1.2132568359,0.0000000000)
        mp.CreateNewNode(12,0.7500000000,1.4007568359,0.0000000000)
        mp.CreateNewNode(13,0.0000000000,1.0000000000,0.0000000000)
        mp.CreateNewNode(14,1.0000000000,1.7500000000,0.0000000000)
        mp.CreateNewNode(15,0.2500000000,1.0000000000,0.0000000000)
        mp.CreateNewNode(16,1.0000000000,1.5625000000,0.0000000000)
        mp.CreateNewNode(17,0.5000000000,1.0000000000,0.0000000000)
        mp.CreateNewNode(18,1.0000000000,1.3750000000,0.0000000000)
        mp.CreateNewNode(19,0.0000000000,0.7500000000,0.0000000000)
        mp.CreateNewNode(20,0.7500000000,1.0000000000,0.0000000000)
        mp.CreateNewNode(21,1.2500000000,1.6875000000,0.0000000000)
        mp.CreateNewNode(22,1.0000000000,1.1875000000,0.0000000000)
        mp.CreateNewNode(23,0.5000000000,0.7867431641,0.0000000000)
        mp.CreateNewNode(24,1.2500000000,1.3492431641,0.0000000000)
        mp.CreateNewNode(25,1.0000000000,1.0000000000,0.0000000000)
        mp.CreateNewNode(26,0.2500000000,0.5389404297,0.0000000000)
        mp.CreateNewNode(27,0.0000000000,0.5000000000,0.0000000000)
        mp.CreateNewNode(28,0.5000000000,0.5712890625,0.0000000000)
        mp.CreateNewNode(29,1.5000000000,1.6250000000,0.0000000000)
        mp.CreateNewNode(30,1.0000000000,0.8125000000,0.0000000000)
        mp.CreateNewNode(31,1.5000000000,1.4764404297,0.0000000000)
        mp.CreateNewNode(32,0.7500000000,0.5992431641,0.0000000000)
        mp.CreateNewNode(33,1.2500000000,1.0000000000,0.0000000000)
        mp.CreateNewNode(34,1.5000000000,1.3212890625,0.0000000000)
        mp.CreateNewNode(35,1.0000000000,0.6250000000,0.0000000000)
        mp.CreateNewNode(36,1.5000000000,1.1617431641,0.0000000000)
        mp.CreateNewNode(37,0.5000000000,0.3514404297,0.0000000000)
        mp.CreateNewNode(38,0.0000000000,0.2500000000,0.0000000000)
        mp.CreateNewNode(39,1.5000000000,1.0000000000,0.0000000000)
        mp.CreateNewNode(40,1.7500000000,1.5625000000,0.0000000000)
        mp.CreateNewNode(41,1.2500000000,0.6507568359,0.0000000000)
        mp.CreateNewNode(42,1.0000000000,0.4375000000,0.0000000000)
        mp.CreateNewNode(43,1.7500000000,1.2889404297,0.0000000000)
        mp.CreateNewNode(44,1.5000000000,0.8382568359,0.0000000000)
        mp.CreateNewNode(45,0.5000000000,0.1250000000,0.0000000000)
        mp.CreateNewNode(46,0.2500000000,0.0625000000,0.0000000000)
        mp.CreateNewNode(47,0.7500000000,0.1875000000,0.0000000000)
        mp.CreateNewNode(48,1.5000000000,0.6787109375,0.0000000000)
        mp.CreateNewNode(49,0.0000000000,0.0000000000,0.0000000000)
        mp.CreateNewNode(50,1.0000000000,0.2500000000,0.0000000000)
        mp.CreateNewNode(51,1.7500000000,1.0000000000,0.0000000000)
        mp.CreateNewNode(52,2.0000000000,1.5000000000,0.0000000000)
        mp.CreateNewNode(53,2.0000000000,1.3750000000,0.0000000000)
        mp.CreateNewNode(54,1.2500000000,0.3125000000,0.0000000000)
        mp.CreateNewNode(55,1.5000000000,0.5235595703,0.0000000000)
        mp.CreateNewNode(56,2.0000000000,1.2500000000,0.0000000000)
        mp.CreateNewNode(57,1.7500000000,0.7110595703,0.0000000000)
        mp.CreateNewNode(58,2.0000000000,1.1250000000,0.0000000000)
        mp.CreateNewNode(59,1.5000000000,0.3750000000,0.0000000000)
        mp.CreateNewNode(60,2.0000000000,1.0000000000,0.0000000000)
        mp.CreateNewNode(61,2.0000000000,0.8750000000,0.0000000000)
        mp.CreateNewNode(62,1.7500000000,0.4375000000,0.0000000000)
        mp.CreateNewNode(63,2.0000000000,0.7500000000,0.0000000000)
        mp.CreateNewNode(64,2.0000000000,0.6250000000,0.0000000000)
        mp.CreateNewNode(65,2.0000000000,0.5000000000,0.0000000000)

        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.REACTION_X,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Y, KratosMultiphysics.REACTION_Y,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Z, KratosMultiphysics.REACTION_Z,mp)

        # Create Element
        mp.CreateNewElement("SmallDisplacementElement2D8N",1, [45,28,27,49,37,26,38,46], mp.GetProperties()[0])
        mp.CreateNewElement("SmallDisplacementElement2D8N",2, [50,35,28,45,42,32,37,47], mp.GetProperties()[0])
        mp.CreateNewElement("SmallDisplacementElement2D8N",3, [59,48,35,50,55,41,42,54], mp.GetProperties()[0])
        mp.CreateNewElement("SmallDisplacementElement2D8N",4, [65,63,48,59,64,57,55,62], mp.GetProperties()[0])
        mp.CreateNewElement("SmallDisplacementElement2D8N",5, [28,17,13,27,23,15,19,26], mp.GetProperties()[0])
        mp.CreateNewElement("SmallDisplacementElement2D8N",6, [35,25,17,28,30,20,23,32], mp.GetProperties()[0])
        mp.CreateNewElement("SmallDisplacementElement2D8N",7, [48,39,25,35,44,33,30,41], mp.GetProperties()[0])
        mp.CreateNewElement("SmallDisplacementElement2D8N",8, [63,60,39,48,61,51,44,57], mp.GetProperties()[0])
        mp.CreateNewElement("SmallDisplacementElement2D8N",9, [17,9,4,13,11,6,8,15], mp.GetProperties()[0])
        mp.CreateNewElement("SmallDisplacementElement2D8N",10, [25,18,9,17,22,12,11,20], mp.GetProperties()[0])
        mp.CreateNewElement("SmallDisplacementElement2D8N",11, [39,34,18,25,36,24,22,33], mp.GetProperties()[0])
        mp.CreateNewElement("SmallDisplacementElement2D8N",12, [60,56,34,39,58,43,36,51], mp.GetProperties()[0])
        mp.CreateNewElement("SmallDisplacementElement2D8N",13, [9,5,1,4,7,3,2,6], mp.GetProperties()[0])
        mp.CreateNewElement("SmallDisplacementElement2D8N",14, [18,14,5,9,16,10,7,12], mp.GetProperties()[0])
        mp.CreateNewElement("SmallDisplacementElement2D8N",15, [34,29,14,18,31,21,16,24], mp.GetProperties()[0])
        mp.CreateNewElement("SmallDisplacementElement2D8N",16, [56,52,29,34,53,40,31,43], mp.GetProperties()[0])

        coeff = 1.0e-2
        self._apply_BCs(mp, coeff)
        self._solve(mp)
        self._check_outputs(mp,dim, coeff)

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

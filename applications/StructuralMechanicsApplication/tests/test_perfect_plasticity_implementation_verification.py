from __future__ import print_function, absolute_import, division
import KratosMultiphysics

import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import KratosMultiphysics.KratosUnittest as KratosUnittest

import math

class TestPerfectPlasticityImplementationVerification(KratosUnittest.TestCase):
    def setUp(self):
        pass
    
    def _add_variables(self,mp):
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.VOLUME_ACCELERATION)

    def _apply_fix(self,mp):
        for node in mp.Nodes:
            node.Fix(KratosMultiphysics.DISPLACEMENT_X)
            node.Fix(KratosMultiphysics.DISPLACEMENT_Y)
            node.Fix(KratosMultiphysics.DISPLACEMENT_Z)

        u = KratosMultiphysics.Vector(3)
        u[0] = 0.0
        u[1] = 0.0
        u[2] = 0.0

        for node in mp.Nodes:
            node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT,0,u)

    def _apply_BCs(self,mp,A,B,b,time, w):
        for node in mp.Nodes:
            node.Fix(KratosMultiphysics.DISPLACEMENT_X)
            node.Fix(KratosMultiphysics.DISPLACEMENT_Y)
            node.Fix(KratosMultiphysics.DISPLACEMENT_Z)

        for node in mp.Nodes:
            xvec = KratosMultiphysics.Vector(3)
            xvec[0] = node.X0
            xvec[1] = node.Y0
            xvec[2] = node.Z0

            xvec1 = KratosMultiphysics.Vector(3)
            xvec1 = time*math.sin(w*time)*xvec

            xvec2 = KratosMultiphysics.Vector(3)
            xvec2 = time*math.sin(2*w*time)*xvec

            u = A*xvec1 + B*xvec2
            u += b

            node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT,0,u)

    def _apply_material_properties(self,mp, constitutive_law_type = "LinearJ2Plasticity3DLaw"):
        # Define properties
        mp.GetProperties()[1].SetValue(KratosMultiphysics.YOUNG_MODULUS, 2.0e11)
        mp.GetProperties()[1].SetValue(KratosMultiphysics.POISSON_RATIO, 0.3)
        mp.GetProperties()[1].SetValue(KratosMultiphysics.YIELD_STRESS, 3.0)

        if constitutive_law_type == "LinearJ2Plasticity3DLaw" or constitutive_law_type == "PlasticityIsotropicKinematicJ2Law":
            mp.GetProperties()[1].SetValue(KratosMultiphysics.REFERENCE_HARDENING_MODULUS, 1.0)
            mp.GetProperties()[1].SetValue(KratosMultiphysics.ISOTROPIC_HARDENING_MODULUS, 0.0)
            mp.GetProperties()[1].SetValue(KratosMultiphysics.INFINITY_HARDENING_MODULUS, 0.0)
            mp.GetProperties()[1].SetValue(KratosMultiphysics.HARDENING_EXPONENT, 1.0)
        else:
            mp.GetProperties()[1].SetValue(KratosMultiphysics.FRACTURE_ENERGY, 1.0e16) # Perfect plasticity
            mp.GetProperties()[1].SetValue(StructuralMechanicsApplication.HARDENING_CURVE, 3) # Perfect plasticity

        g = [0,0,0]
        mp.GetProperties()[1].SetValue(KratosMultiphysics.VOLUME_ACCELERATION,g)

        cl = KratosMultiphysics.KratosGlobals.GetConstitutiveLaw(constitutive_law_type).Clone()
        mp.GetProperties()[1].SetValue(KratosMultiphysics.CONSTITUTIVE_LAW,cl)

    def _define_movement(self):
        #define the applied motion - the idea is that the displacement is defined as u = A*xnode + b
        #so that the displcement is linear and the exact F = I + A
        A = KratosMultiphysics.Matrix(3,3)
        A[0,0] = 1.0e-10;   A[0,1] = 2.0e-10; A[0,2] = 0.0
        A[1,0] = 0.5e-10;   A[1,1] = 0.7e-10; A[1,2] = 0.1e-10
        A[2,0] = -0.2e-10;  A[2,1] = 0.0;     A[2,2] = -0.3e-10

        B = KratosMultiphysics.Matrix(3,3)
        B[0,0] = 1.0e-10;   B[0,1] = 2.0e-10; B[0,2] = 0.0
        B[1,0] = 0.5e-10;   B[1,1] = 0.7e-10; B[1,2] = 0.1e-10
        B[2,0] = -0.2e-10;  B[2,1] = 0.0;     B[2,2] = -0.3e-10

        b = KratosMultiphysics.Vector(3)
        b[0] = 0.0
        b[1] = 0.0
        b[2] = 0.0
        #b[0] = 0.5e-10
        #b[1] = -0.2e-10
        #b[2] = 0.7e-10

        return A,B,b

    def _solve(self,mp):

        # Define a minimal newton raphson solver
        linear_solver = KratosMultiphysics.SkylineLUFactorizationSolver()
        builder_and_solver = KratosMultiphysics.ResidualBasedBlockBuilderAndSolver(linear_solver)
        scheme = KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme()
        convergence_criterion = KratosMultiphysics.ResidualCriteria(1e-14,1e-20)

        max_iters = 20
        compute_reactions = True
        reform_step_dofs = True
        calculate_norm_dx = False
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
        convergence_criterion.SetEchoLevel(0)
        strategy.SetEchoLevel(0)

        strategy.Check()
        strategy.Solve()

    def _create_check_outputs(self,current_model):
        import from_json_check_result_process

        check_parameters = KratosMultiphysics.Parameters("""
        {
            "gauss_points_check_variables": ["VON_MISES_STRESS"],
            "input_file_name"      : "",
            "time_frequency"       : 0.01,
            "model_part_name"      : "solid_part",
            "sub_model_part_name"  : "Body"
        }
        """)

        check_parameters["input_file_name"].SetString("cl_test/test_perfect_plasticity_implementation_verification_reference.json")

        check = from_json_check_result_process.FromJsonCheckResultProcess(current_model, check_parameters)
        check.ExecuteInitialize()
        check.ExecuteBeforeSolutionLoop()

        return check

    def _create_reference_solution(self, current_model):
        # The following is used to create the solution database
        import json_output_process

        out_parameters = KratosMultiphysics.Parameters("""
        {
            "gauss_points_output_variables": ["VON_MISES_STRESS"],
            "output_file_name"     : "",
            "time_frequency"       : 0.01,
            "model_part_name"      : "solid_part",
            "sub_model_part_name"  : "Body"
        }
        """)

        out_parameters["output_file_name"].SetString("cl_test/test_perfect_plasticity_implementation_verification_reference.json")

        out = json_output_process.JsonOutputProcess(current_model, out_parameters)
        out.ExecuteInitialize()
        out.ExecuteBeforeSolutionLoop()

        return out

    def _base_test(self, constitutive_law_type, debug = False):
        current_model = KratosMultiphysics.Model()
        mp = current_model.CreateModelPart("solid_part")
        self._add_variables(mp)
        self._apply_material_properties(mp, constitutive_law_type)

        # Create nodes
        mp.CreateNewNode(1, 0.0000000000, 0.0000000000, 0.0000000000)
        mp.CreateNewNode(2, 1.0000000000, 0.0000000000, 0.0000000000)
        mp.CreateNewNode(3, 1.0000000000, 1.0000000000, 0.0000000000)
        mp.CreateNewNode(4, 0.0000000000, 1.0000000000, 0.0000000000)
        mp.CreateNewNode(5, 0.0000000000, 0.0000000000, 1.0000000000)
        mp.CreateNewNode(6, 1.0000000000, 0.0000000000, 1.0000000000)
        mp.CreateNewNode(7, 1.0000000000, 1.0000000000, 1.0000000000)
        mp.CreateNewNode(8, 0.0000000000, 1.0000000000, 1.0000000000)

        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.REACTION_X,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Y, KratosMultiphysics.REACTION_Y,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Z, KratosMultiphysics.REACTION_Z,mp)

        # Create a submodelpart for boundary conditions
        bcs = mp.CreateSubModelPart("BoundaryCondtions")
        bcs.AddNodes([5,6,7,8])
        fix = mp.CreateSubModelPart("FixCondtions")
        fix.AddNodes([1,2,3,4])

        # Create Element
        mp.CreateNewElement("SmallDisplacementBbarElement3D8N", 1,[1,2,3,4,5,6,7,8], mp.GetProperties()[1])

        # Create model part
        body = mp.CreateSubModelPart("Body")
        body.AddElements([1])

        A,B,b = self._define_movement()
        self._apply_fix(fix)

        if debug:
            gid_output = self.__create_post_process(mp, constitutive_law_type)

        w = 5.0 # Frequency  movement
        dt = 0.01 # Delta time
        step = 0
        time = 0.0;
        end_time = 1.0

        while time < end_time:
            time = time + dt
            step += 1
            mp.ProcessInfo[KratosMultiphysics.STEP] = step
            mp.CloneTimeStep(time)
            self._apply_BCs(bcs,A,B,b,time,w)
            self._solve(mp)
            if step == 1:
                check = self._create_check_outputs(current_model)
            check.ExecuteFinalizeSolutionStep()

            #if step == 1:
                #out = self._create_reference_solution(current_model)
            #out.ExecuteFinalizeSolutionStep()

            if debug:
                gid_output.PrintOutput()

        if debug:
            gid_output.ExecuteFinalizeSolutionStep()
            gid_output.ExecuteFinalize()

    def test_PerfectPlasticityLinearJ2Plasticity3DLaw(self):
        self._base_test("LinearJ2Plasticity3DLaw")

    def test_PerfectPlasticityPlasticityIsotropicKinematicJ2Law(self):
        self._base_test("PlasticityIsotropicKinematicJ2Law")

    def test_PerfectPlasticitySmallStrainIsotropicPlasticity3DVonMisesVonMises(self):
        self._base_test("SmallStrainIsotropicPlasticity3DVonMisesVonMises")

    def __create_post_process(self, main_model_part, constitutive_law_type):
        from gid_output_process import GiDOutputProcess
        gid_output = GiDOutputProcess(main_model_part,
                                    "gid_output_" + constitutive_law_type,
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
                                                "gauss_point_results" : ["VON_MISES_STRESS","GREEN_LAGRANGE_STRAIN_VECTOR","PK2_STRESS_VECTOR","PLASTIC_DISSIPATION","PLASTIC_STRAIN","EQUIVALENT_PLASTIC_STRAIN","PLASTIC_STRAIN_VECTOR","UNIAXIAL_STRESS"]
                                            }
                                        }
                                        """)
                                    )

        gid_output.ExecuteInitialize()
        gid_output.ExecuteBeforeSolutionLoop()
        gid_output.ExecuteInitializeSolutionStep()

        return gid_output

if __name__ == '__main__':
    KratosUnittest.main()

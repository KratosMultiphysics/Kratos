import KratosMultiphysics

import KratosMultiphysics.KratosUnittest as KratosUnittest

from KratosMultiphysics.kratos_utilities import CheckIfApplicationsAvailable
if CheckIfApplicationsAvailable("StructuralMechanicsApplication"):
    from KratosMultiphysics import StructuralMechanicsApplication

from KratosMultiphysics import ConstitutiveLawsApplication

import os
import math

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

@KratosUnittest.skipIfApplicationsNotAvailable("StructuralMechanicsApplication")
class TestPerfectPlasticityImplementationVerification(KratosUnittest.TestCase):
    def setUp(self):
        pass

    def _base_test(self, constitutive_law_type, debug = False):
        current_model = KratosMultiphysics.Model()
        mp = current_model.CreateModelPart("solid_part")
        _add_variables(mp)
        _apply_material_properties(mp, constitutive_law_type)

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

        A,B,b = _define_movement()
        _apply_fix(fix)

        if debug:
            output = _create_post_process(mp, constitutive_law_type)

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
            _apply_BCs(bcs,A,B,b,time,w)
            _solve(mp)
            if step == 1:
                check = _create_check_outputs(current_model)
            check.ExecuteFinalizeSolutionStep()

            #if step == 1:
                #out = _create_reference_solution(current_model)
            #out.ExecuteFinalizeSolutionStep()

            if debug:
                output.PrintOutput()

        if debug:
            output.ExecuteFinalizeSolutionStep()
            output.ExecuteFinalize()

    def test_PerfectPlasticitySmallStrainJ2Plasticity3DLaw(self):
        self._base_test("SmallStrainJ2Plasticity3DLaw")

    def test_PerfectPlasticityPlasticityIsotropicKinematicJ2Law(self):
        self._base_test("PlasticityIsotropicKinematicJ2Law")

    def test_PerfectPlasticitySmallStrainIsotropicPlasticity3DVonMisesVonMises(self):
        self._base_test("SmallStrainIsotropicPlasticity3DVonMisesVonMises")

    def test_PerfectPlasticityFiniteStrainIsotropicPlasticity3DVonMisesVonMises(self):
        self._base_test("FiniteStrainIsotropicPlasticity3DVonMisesVonMises")

def _apply_BCs(mp,A,B,b,time, w):
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

def _add_variables(mp):
    mp.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
    mp.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)
    mp.AddNodalSolutionStepVariable(KratosMultiphysics.VOLUME_ACCELERATION)

def _apply_fix(mp):
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

def _apply_material_properties(mp, constitutive_law_type = "SmallStrainJ2Plasticity3DLaw"):
    # Define properties
    mp.GetProperties()[1].SetValue(KratosMultiphysics.YOUNG_MODULUS, 2.0e11)
    mp.GetProperties()[1].SetValue(KratosMultiphysics.POISSON_RATIO, 0.3)
    mp.GetProperties()[1].SetValue(KratosMultiphysics.YIELD_STRESS, 9.0)
    mp.GetProperties()[1].SetValue(ConstitutiveLawsApplication.MAX_NUMBER_NL_CL_ITERATIONS, 100)

    if constitutive_law_type == "SmallStrainJ2Plasticity3DLaw" or constitutive_law_type == "PlasticityIsotropicKinematicJ2Law":
        mp.GetProperties()[1].SetValue(KratosMultiphysics.ISOTROPIC_HARDENING_MODULUS, 0.0)
        mp.GetProperties()[1].SetValue(ConstitutiveLawsApplication.EXPONENTIAL_SATURATION_YIELD_STRESS, 9.0)
        mp.GetProperties()[1].SetValue(KratosMultiphysics.HARDENING_EXPONENT, 1.0)
    else:
        mp.GetProperties()[1].SetValue(KratosMultiphysics.FRACTURE_ENERGY, 1.0e16) # Perfect plasticity
        mp.GetProperties()[1].SetValue(ConstitutiveLawsApplication.HARDENING_CURVE, 3) # Perfect plasticity

    g = [0,0,0]
    mp.GetProperties()[1].SetValue(KratosMultiphysics.VOLUME_ACCELERATION,g)

    cl = KratosMultiphysics.KratosGlobals.GetConstitutiveLaw(constitutive_law_type).Clone()
    mp.GetProperties()[1].SetValue(KratosMultiphysics.CONSTITUTIVE_LAW,cl)

def _define_movement():
    #define the applied motion - the idea is that the displacement is defined as u = A*xnode + b
    #so that the displcement is linear and the exact F = I + A
    A = KratosMultiphysics.Matrix(3,3)
    A[0,0] = 1.0e-10;   A[0,1] = 2.0e-10; A[0,2] = 0.0
    A[1,0] = 0.5e-10;   A[1,1] = 0.7e-10; A[1,2] = 0.1e-10
    A[2,0] = -0.2e-10;  A[2,1] = 0.0;     A[2,2] = -0.3e-10

    B = KratosMultiphysics.Matrix(3,3)
    B[0,0] = 0.0;     B[0,1] = 7.0e-10; B[0,2] = 0.0
    B[1,0] = 2.5e-10; B[1,1] = 1.7e-10; B[1,2] = 0.1e-10
    B[2,0] = 0.0;     B[2,1] = 1.0e-10;     B[2,2] = -0.3e-10

    b = KratosMultiphysics.Vector(3)
    b[0] = 0.0
    b[1] = 0.0
    b[2] = 0.0
    #b[0] = 0.5e-10
    #b[1] = -0.2e-10
    #b[2] = 0.7e-10

    return A,B,b

def _solve(mp):

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

def _create_check_outputs(current_model):
    from KratosMultiphysics.from_json_check_result_process import FromJsonCheckResultProcess

    check_parameters = KratosMultiphysics.Parameters("""
    {
        "gauss_points_check_variables": ["VON_MISES_STRESS"],
        "input_file_name"      : "",
        "time_frequency"       : 0.01,
        "model_part_name"      : "solid_part",
        "sub_model_part_name"  : "Body"
    }
    """)

    check_parameters["input_file_name"].SetString(GetFilePath("test_perfect_plasticity_implementation_verification_reference.json"))

    check = FromJsonCheckResultProcess(current_model, check_parameters)
    check.ExecuteInitialize()
    check.ExecuteBeforeSolutionLoop()

    return check

def _create_reference_solution(current_model):
    # The following is used to create the solution database
    from KratosMultiphysics.json_output_process import JsonOutputProcess

    out_parameters = KratosMultiphysics.Parameters("""
    {
        "gauss_points_output_variables": ["VON_MISES_STRESS"],
        "output_file_name"     : "",
        "time_frequency"       : 0.01,
        "model_part_name"      : "solid_part",
        "sub_model_part_name"  : "Body"
    }
    """)

    out_parameters["output_file_name"].SetString(GetFilePath("test_perfect_plasticity_implementation_verification_reference.json"))

    out = JsonOutputProcess(current_model, out_parameters)
    out.ExecuteInitialize()
    out.ExecuteBeforeSolutionLoop()

    return out

def _create_post_process(main_model_part, constitutive_law_type, debug = "gid"):
    if debug == "gid":
        from KratosMultiphysics.gid_output_process import GiDOutputProcess
        output = GiDOutputProcess(main_model_part,
                                    "output_" + constitutive_law_type,
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
    elif debug == "vtk":
        from KratosMultiphysics.vtk_output_process import VtkOutputProcess
        output = VtkOutputProcess(main_model_part.GetModel(),
                                    KratosMultiphysics.Parameters("""{
                                            "model_part_name"                    : "solid_part",
                                            "nodal_solution_step_data_variables" : ["DISPLACEMENT"],
                                            "gauss_point_variables"              : ["VON_MISES_STRESS"]
                                        }
                                        """)
                                    )

    output.ExecuteInitialize()
    output.ExecuteBeforeSolutionLoop()
    output.ExecuteInitializeSolutionStep()

    return output

if __name__ == '__main__':
    KratosUnittest.main()

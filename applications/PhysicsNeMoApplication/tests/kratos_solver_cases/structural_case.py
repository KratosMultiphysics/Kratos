"""Real linear-static cantilever case built entirely in memory.

A thin rectangular strip (length 1, height 0.1) of plane-strain linear
elastic material, clamped at x = 0 and loaded by a vertical tip force
spread over the x = 1 edge nodes: the classic cantilever solved by
StructuralMechanicsApplication's static solver. The structural counterpart
of thermal_case.py - everything parametrized, no .mdpa fixtures.

Requires StructuralMechanicsApplication (and LinearSolversApplication for
the default linear solver); callers gate on availability via
kratos_utilities.CheckIfApplicationsAvailable.
"""

import KratosMultiphysics as Kratos

_LENGTH = 1.0
_HEIGHT = 0.1

# The complete variable set the static solver's AddVariables adds (rotation,
# volumetric-strain and strain dofs off) - with use_input_model_part the
# nodes exist before AddVariables runs, so everything must be pre-added.
_CORE_HISTORICAL_VARIABLES = (
    "DISPLACEMENT", "REACTION", "POSITIVE_FACE_PRESSURE",
    "NEGATIVE_FACE_PRESSURE", "VOLUME_ACCELERATION",
    # the implicit dynamic solver's _add_dynamic_variables(); like the rest
    # they must exist before the nodes do (use_input_model_part)
    "VELOCITY", "ACCELERATION",
)
_APP_HISTORICAL_VARIABLES = ("POINT_LOAD", "LINE_LOAD", "SURFACE_LOAD")


def _CreateProjectParameters(echo_level=0, dynamic_settings=None):
    """Static by default; dynamic_settings (dict with time_step,
    end_time, damp_factor_m, newmark_beta) switches to the implicit
    dynamic solver (Bossak)."""
    parameters = Kratos.Parameters("""{
        "problem_data": {
            "problem_name"  : "physics_nemo_structural_case",
            "parallel_type" : "OpenMP",
            "start_time"    : 0.0,
            "end_time"      : 0.99,
            "echo_level"    : %d
        },
        "solver_settings": {
            "solver_type"                        : "Static",
            "analysis_type"                      : "linear",
            "model_part_name"                    : "StructuralModelPart",
            "domain_size"                        : 2,
            "model_import_settings"              : { "input_type" : "use_input_model_part" },
            "material_import_settings"           : { "materials_filename" : "" },
            "echo_level"                         : %d,
            "rotation_dofs"                      : false,
            "time_stepping"                      : { "time_step" : 1.0 }
        },
        "processes"        : {},
        "output_processes" : {}
    }""" % (echo_level, echo_level))

    if dynamic_settings is None:
        return parameters

    solver_settings = parameters["solver_settings"]
    solver_settings["solver_type"].SetString("dynamic")
    solver_settings["time_stepping"]["time_step"].SetDouble(
        float(dynamic_settings.get("time_step", 0.005)))
    parameters["problem_data"]["end_time"].SetDouble(
        float(dynamic_settings.get("end_time", 0.05)))
    solver_settings.AddEmptyValue("time_integration_method").SetString("implicit")
    solver_settings.AddEmptyValue("scheme_type").SetString(
        dynamic_settings.get("scheme_type", "bossak"))
    solver_settings.AddEmptyValue("damp_factor_m").SetDouble(
        float(dynamic_settings.get("damp_factor_m", -0.3)))
    solver_settings.AddEmptyValue("newmark_beta").SetDouble(
        float(dynamic_settings.get("newmark_beta", 0.25)))
    return parameters


def CreateStructuralModelPart(model: Kratos.Model, divisions: int = 6) -> Kratos.ModelPart:
    """Creates the meshed cantilever model part (variables added before the
    mesh; elements are real SmallDisplacementElement2D3N, no replacement)."""
    import KratosMultiphysics.StructuralMechanicsApplication as SMA

    model_part = model.CreateModelPart("StructuralModelPart")
    model_part.ProcessInfo[Kratos.DOMAIN_SIZE] = 2
    model_part.SetBufferSize(2)
    for name in _CORE_HISTORICAL_VARIABLES:
        model_part.AddNodalSolutionStepVariable(Kratos.KratosGlobals.GetVariable(name))
    for name in _APP_HISTORICAL_VARIABLES:
        model_part.AddNodalSolutionStepVariable(getattr(SMA, name))

    generator_geometry = Kratos.Quadrilateral2D4(
        Kratos.Node(1, 0.0, 0.0, 0.0), Kratos.Node(2, 0.0, _HEIGHT, 0.0),
        Kratos.Node(3, _LENGTH, _HEIGHT, 0.0), Kratos.Node(4, _LENGTH, 0.0, 0.0))
    mesh_parameters = Kratos.Parameters("""{
        "number_of_divisions"        : %d,
        "element_name"               : "SmallDisplacementElement2D3N",
        "condition_name"             : "LineCondition",
        "create_skin_sub_model_part" : false
    }""" % divisions)
    domain = model_part.CreateSubModelPart("Domain")
    Kratos.StructuredMeshGeneratorProcess(generator_geometry, domain, mesh_parameters).Execute()

    # material data goes onto the properties the generated elements reference
    properties = next(iter(domain.Elements)).Properties
    properties.SetValue(Kratos.YOUNG_MODULUS, 210.0e9)
    properties.SetValue(Kratos.POISSON_RATIO, 0.3)
    properties.SetValue(Kratos.CONSTITUTIVE_LAW, SMA.LinearElasticPlaneStrain2DLaw())
    return model_part


def ApplyCaseData(model_part: Kratos.ModelPart, tip_load: float, tolerance: float = 1e-8) -> None:
    """Clamps the x = 0 edge and spreads a downward tip_load (total, in N)
    over PointLoadCondition2D1N conditions on the x = LENGTH edge nodes."""
    import KratosMultiphysics.StructuralMechanicsApplication as SMA

    tip_nodes = [node for node in model_part.Nodes if abs(node.X - _LENGTH) < tolerance]
    properties = next(iter(model_part.Elements)).Properties
    condition_id = model_part.NumberOfConditions()
    for node in model_part.Nodes:
        if abs(node.X) < tolerance:
            node.Fix(Kratos.DISPLACEMENT_X)
            node.Fix(Kratos.DISPLACEMENT_Y)
            node.SetSolutionStepValue(Kratos.DISPLACEMENT, [0.0, 0.0, 0.0])
    domain = model_part.GetSubModelPart("Domain")
    for node in tip_nodes:
        node.SetSolutionStepValue(SMA.POINT_LOAD, [0.0, -tip_load / len(tip_nodes), 0.0])
        condition_id += 1
        domain.CreateNewCondition("PointLoadCondition2D1N", condition_id, [node.Id], properties)


def CreateStructuralAnalysis(model: Kratos.Model,
                             tip_load: float = 1.0e6,
                             divisions: int = 6,
                             echo_level: int = 0):
    """Builds the meshed cantilever and returns a ready
    StructuralMechanicsAnalysis.

    Call .Run() on the result; the solved DISPLACEMENT field lives on
    model["StructuralModelPart"].
    """
    from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import StructuralMechanicsAnalysis

    model_part = CreateStructuralModelPart(model, divisions)
    ApplyCaseData(model_part, tip_load)
    return StructuralMechanicsAnalysis(model, _CreateProjectParameters(echo_level))


def GetTipDeflection(model_part: Kratos.ModelPart, tolerance: float = 1e-8) -> float:
    """Mean vertical displacement of the x = LENGTH edge nodes (selected by
    the initial position X0: the static solver moves the mesh)."""
    tip_values = [node.GetSolutionStepValue(Kratos.DISPLACEMENT_Y)
                  for node in model_part.Nodes if abs(node.X0 - _LENGTH) < tolerance]
    return sum(tip_values) / len(tip_values)


def CreateTransientStructuralAnalysis(model: Kratos.Model,
                                      tip_load: float = 1.0e6,
                                      divisions: int = 6,
                                      time_step: float = 0.005,
                                      end_time: float = 0.05,
                                      damp_factor_m: float = -0.3,
                                      echo_level: int = 0):
    """The same cantilever run by the IMPLICIT DYNAMIC solver (Bossak).

    The suddenly applied tip load makes the strip oscillate, so the
    collected trajectory is a genuine transient - the crash/deformation
    surrogate setting. The scheme is a real displacement time-integration
    scheme (StructuralMechanicsBossakScheme), so residual assembly on this
    case needs differentiable_residual's scheme= path.
    """
    from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import StructuralMechanicsAnalysis

    model_part = CreateStructuralModelPart(model, divisions)
    # DENSITY is set only here: the dynamic mass matrix needs it, while the
    # static case is left exactly as its own tests pinned it.
    next(iter(model_part.Elements)).Properties.SetValue(Kratos.DENSITY, 7850.0)
    ApplyCaseData(model_part, tip_load)
    return StructuralMechanicsAnalysis(model, _CreateProjectParameters(
        echo_level, {"time_step": time_step, "end_time": end_time,
                     "damp_factor_m": damp_factor_m}))


def CollectDisplacements(model_part: Kratos.ModelPart):
    """(N, 2) planar DISPLACEMENT state, for RunTransientAnalysis(collect=...)."""
    return [[node.GetSolutionStepValue(Kratos.DISPLACEMENT_X),
             node.GetSolutionStepValue(Kratos.DISPLACEMENT_Y)]
            for node in model_part.Nodes]

"""Real thermo-mechanical sintering case built entirely in memory.

A unit square held at a uniform reference temperature is cooled through its
boundary; thermal contraction pulls the body inward. That is sintering's
defining signal - a temperature schedule driving progressive, spatially
non-uniform shrinkage - and it is what the VFGN surrogate consumes, in place
of the synthetic shrinking clouds the bridge's own tests use.

Solved by ConvectionDiffusionApplication's CoupledThermoMechanicalSolver: a
one-pass staggered step (thermal, then structural). There is no mapper
between them. ConnectivityPreserveModeler makes the thermal model part share
the structural part's NODES, PROPERTIES and PROCESSINFO, so the coupling is
simply the shared nodal TEMPERATURE, read at Gauss points by a thermal
constitutive law and turned into strain through THERMAL_EXPANSION_COEFFICIENT
and REFERENCE_TEMPERATURE.

Requires ConvectionDiffusionApplication, StructuralMechanicsApplication,
ConstitutiveLawsApplication (which registers the Thermal* laws - Structural
registers none) and LinearSolversApplication; callers gate on availability
via kratos_utilities.CheckIfApplicationsAvailable.
"""

import KratosMultiphysics as Kratos

# The UNION of both sub-solvers' AddVariables lists. This case is stricter
# than the single-physics ones: MergeVariableListsUtility::Merge unions the
# two lists by calling AddNodalSolutionStepVariable on the already-meshed
# structural part, so anything missing here raises "Attempting to add the
# variable ... to the model part ... which is not empty" - and it fires at
# ANALYSIS CONSTRUCTION, before Initialize().
_STRUCTURAL_HISTORICAL_VARIABLES = (
    "DISPLACEMENT", "REACTION", "POSITIVE_FACE_PRESSURE", "NEGATIVE_FACE_PRESSURE",
    "VOLUME_ACCELERATION", "POINT_LOAD", "LINE_LOAD", "SURFACE_LOAD",
)
_THERMAL_HISTORICAL_VARIABLES = (
    "TEMPERATURE", "DENSITY", "SPECIFIC_HEAT", "CONDUCTIVITY", "HEAT_FLUX",
    "FACE_HEAT_FLUX", "PROJECTED_SCALAR1", "CONVECTION_VELOCITY",
    "TEMPERATURE_GRADIENT", "TRANSFER_COEFFICIENT", "VELOCITY", "MESH_VELOCITY",
    "REACTION_FLUX",
)

# The thermal part must NOT be created by the caller: the modeler needs an
# empty destination and the thermal sub-solver creates it itself.
_STRUCTURAL_MODEL_PART_NAME = "Structure"
_COOLED_SKIN_NAME = "CooledSkin"


def _CreateProjectParameters(reference_temperature: float,
                             cooling_rate: float,
                             time_step: float,
                             end_time: float,
                             echo_level: int = 0) -> Kratos.Parameters:
    """The coupled configuration, with a time-dependent cooling ramp.

    domain_size is needed in THREE places: at the top level (it selects
    EulerianConvDiff2D vs 3D for the modeler) and in both sub-blocks - the
    thermal default of -1 raises.
    """
    return Kratos.Parameters("""{
        "problem_data": {
            "problem_name"  : "physics_nemo_sintering_case",
            "parallel_type" : "OpenMP",
            "start_time"    : 0.0,
            "end_time"      : %f,
            "echo_level"    : %d
        },
        "solver_settings": {
            "solver_type" : "ThermoMechanicallyCoupled",
            "domain_size" : 2,
            "echo_level"  : %d,
            "structural_solver_settings": {
                "solver_type"              : "Static",
                "model_part_name"          : "%s",
                "domain_size"              : 2,
                "echo_level"               : %d,
                "analysis_type"            : "linear",
                "model_import_settings"    : { "input_type" : "use_input_model_part" },
                "material_import_settings" : { "materials_filename" : "" },
                "time_stepping"            : { "time_step" : %f },
                "rotation_dofs"            : false
            },
            "thermal_solver_settings": {
                "solver_type"              : "transient",
                "analysis_type"            : "linear",
                "model_part_name"          : "ThermalModelPart",
                "domain_size"              : 2,
                "echo_level"               : %d,
                "model_import_settings"    : { "input_type" : "use_input_model_part" },
                "material_import_settings" : { "materials_filename" : "" },
                "time_stepping"            : { "time_step" : %f },
                "transient_parameters"     : {
                    "theta"                          : 0.5,
                    "dynamic_tau"                    : 0.0,
                    "cross_wind_stabilization_factor": 0.0
                }
            }
        },
        "processes": {
            "constraints_process_list": [ {
                "python_module" : "assign_scalar_variable_process",
                "kratos_module" : "KratosMultiphysics",
                "Parameters"    : {
                    "model_part_name" : "%s.%s",
                    "variable_name"   : "TEMPERATURE",
                    "constrained"     : true,
                    "value"           : "%f - %f*t",
                    "interval"        : [0.0, "End"]
                }
            } ]
        },
        "output_processes" : {}
    }""" % (end_time, echo_level, echo_level,
            _STRUCTURAL_MODEL_PART_NAME, echo_level, time_step,
            echo_level, time_step,
            _STRUCTURAL_MODEL_PART_NAME, _COOLED_SKIN_NAME,
            reference_temperature, cooling_rate))


def CreateSinteringModelPart(model: Kratos.Model, divisions: int = 6,
                             reference_temperature: float = 1000.0,
                             thermal_expansion: float = 2.0e-4) -> Kratos.ModelPart:
    """Builds the meshed structural part (the thermal one is generated later)."""
    import KratosMultiphysics.ConstitutiveLawsApplication as ConstitutiveLaws
    import KratosMultiphysics.StructuralMechanicsApplication  # noqa: F401 (registers elements)

    model_part = model.CreateModelPart(_STRUCTURAL_MODEL_PART_NAME)
    model_part.ProcessInfo[Kratos.DOMAIN_SIZE] = 2
    model_part.SetBufferSize(2)
    for name in _STRUCTURAL_HISTORICAL_VARIABLES + _THERMAL_HISTORICAL_VARIABLES:
        model_part.AddNodalSolutionStepVariable(Kratos.KratosGlobals.GetVariable(name))

    generator_geometry = Kratos.Quadrilateral2D4(
        Kratos.Node(1, 0.0, 0.0, 0.0), Kratos.Node(2, 0.0, 1.0, 0.0),
        Kratos.Node(3, 1.0, 1.0, 0.0), Kratos.Node(4, 1.0, 0.0, 0.0))
    mesh_parameters = Kratos.Parameters("""{
        "number_of_divisions"        : %d,
        "element_name"               : "SmallDisplacementElement2D3N",
        "condition_name"             : "LineCondition",
        "create_skin_sub_model_part" : false
    }""" % divisions)
    domain = model_part.CreateSubModelPart("Domain")
    Kratos.StructuredMeshGeneratorProcess(
        generator_geometry, domain, mesh_parameters).Execute()

    properties = next(iter(domain.Elements)).Properties
    properties.SetValue(Kratos.YOUNG_MODULUS, 2.0e9)
    properties.SetValue(Kratos.POISSON_RATIO, 0.3)
    properties.SetValue(Kratos.DENSITY, 7800.0)
    # must be NON-NEGATIVE: the law's Check rejects negatives, so shrinkage
    # comes from cooling BELOW the reference temperature, not from a sign flip
    properties.SetValue(Kratos.THERMAL_EXPANSION_COEFFICIENT, thermal_expansion)
    properties.SetValue(Kratos.REFERENCE_TEMPERATURE, reference_temperature)
    properties.SetValue(Kratos.CONSTITUTIVE_LAW, ConstitutiveLaws.ThermalLinearPlaneStrain())
    return model_part


def ApplyCaseData(model_part: Kratos.ModelPart, reference_temperature: float = 1000.0,
                  tolerance: float = 1e-8) -> None:
    """Uniform hot start, a cooled skin, and the minimum restraint.

    The body must be free to contract, so only rigid-body motion is removed:
    one node pinned in both directions and one more in Y.
    """
    skin = model_part.CreateSubModelPart(_COOLED_SKIN_NAME)
    skin_ids = []
    for node in model_part.Nodes:
        # both buffer steps: the transient solver differences in time
        node.SetSolutionStepValue(Kratos.TEMPERATURE, 0, reference_temperature)
        node.SetSolutionStepValue(Kratos.TEMPERATURE, 1, reference_temperature)
        node.SetSolutionStepValue(Kratos.DENSITY, 1.0)
        node.SetSolutionStepValue(Kratos.SPECIFIC_HEAT, 1.0)
        node.SetSolutionStepValue(Kratos.CONDUCTIVITY, 1.0)
        on_boundary = (abs(node.X) < tolerance or abs(node.X - 1.0) < tolerance
                       or abs(node.Y) < tolerance or abs(node.Y - 1.0) < tolerance)
        if on_boundary:
            skin_ids.append(node.Id)
    skin.AddNodes(skin_ids)

    centre = min(model_part.Nodes,
                 key=lambda node: (node.X - 0.5) ** 2 + (node.Y - 0.5) ** 2)
    centre.Fix(Kratos.DISPLACEMENT_X)
    centre.Fix(Kratos.DISPLACEMENT_Y)
    for node in model_part.Nodes:
        if node.Id != centre.Id and abs(node.Y - centre.Y) < tolerance:
            node.Fix(Kratos.DISPLACEMENT_Y)   # kills the remaining rotation
            break


def CreateSinteringAnalysis(model: Kratos.Model,
                            divisions: int = 6,
                            reference_temperature: float = 1000.0,
                            cooling_rate: float = 1600.0,
                            thermal_expansion: float = 2.0e-4,
                            time_step: float = 0.05,
                            end_time: float = 0.5,
                            echo_level: int = 0):
    """Builds the meshed body and returns a ready ConvectionDiffusionAnalysis.

    There is no coupled analysis stage upstream - the plain
    ConvectionDiffusionAnalysis drives the coupled solver. Run it with
    transient_harness.RunTransientAnalysis(analysis, collect=CollectPositions)
    to get a (T, N, 3) trajectory.
    """
    from KratosMultiphysics.ConvectionDiffusionApplication.convection_diffusion_analysis import (
        ConvectionDiffusionAnalysis)

    model_part = CreateSinteringModelPart(
        model, divisions, reference_temperature, thermal_expansion)
    ApplyCaseData(model_part, reference_temperature)
    return ConvectionDiffusionAnalysis(model, _CreateProjectParameters(
        reference_temperature, cooling_rate, time_step, end_time, echo_level))


def CollectPositions(model_part: Kratos.ModelPart):
    """(N, 3) deformed node positions, for RunTransientAnalysis(collect=...).

    The structural solver moves the mesh (move_mesh_flag defaults true), so
    node.X/Y/Z are the deformed coordinates VFGN consumes; X0/Y0/Z0 remain
    the reference configuration.
    """
    return [[node.X, node.Y, node.Z] for node in model_part.Nodes]


def CollectTemperatures(model_part: Kratos.ModelPart):
    """(N,) nodal TEMPERATURE - the schedule driving the shrinkage."""
    return [node.GetSolutionStepValue(Kratos.TEMPERATURE) for node in model_part.Nodes]

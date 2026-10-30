"""A SIMP compliance case: density in, stiffness and compliance out.

The ground truth behind generative topology design. Each element of a
plane-stress square carries its own Young's modulus through the SIMP
interpolation E = E_min + rho^p (E_0 - E_min), so a DENSITY FIELD is the
design variable and the structure's compliance f.u is the thing a design is
judged by. That gives a generative model something real to be conditioned
on and scored against, with no optimization application required.

Built entirely in memory, like the other cases here, so every solution-step
variable is added before the mesh exists.

Per-element properties are the mechanism: the structured mesh generator
gives every element ONE shared Properties object, so each element is given
a fresh one (the constitutive law included, since the law lives on the
properties) before its modulus is set.
"""

import numpy

import KratosMultiphysics as Kratos

_LENGTH = 1.0
_HEIGHT = 1.0
_YOUNG_MODULUS = 1.0e9
_MINIMUM_MODULUS = 1.0e6      # keeps the void stiff enough to stay solvable
_POISSON_RATIO = 0.3
_PENALIZATION = 3.0           # the "p" of SIMP
_TIP_LOAD = -1.0e4

# the UNION of what the static solver's AddVariables asks for: with
# use_input_model_part they must all exist before the mesh does, or
# construction fails with "which is not empty" before Initialize() runs
_CORE_HISTORICAL_VARIABLES = (
    "DISPLACEMENT", "REACTION", "POSITIVE_FACE_PRESSURE",
    "NEGATIVE_FACE_PRESSURE", "VOLUME_ACCELERATION", "VELOCITY", "ACCELERATION",
)
_APP_HISTORICAL_VARIABLES = ("POINT_LOAD", "LINE_LOAD", "SURFACE_LOAD")


def _CreateProjectParameters(echo_level: int = 0) -> Kratos.Parameters:
    return Kratos.Parameters("""{
        "problem_data" : {
            "problem_name"  : "compliance",
            "parallel_type" : "OpenMP",
            "echo_level"    : %d,
            "start_time"    : 0.0,
            "end_time"      : 1.0
        },
        "solver_settings" : {
            "solver_type"              : "Static",
            "model_part_name"          : "ComplianceModelPart",
            "domain_size"              : 2,
            "echo_level"               : %d,
            "analysis_type"            : "linear",
            "model_import_settings"    : { "input_type" : "use_input_model_part" },
            "material_import_settings" : { "materials_filename" : "" },
            "time_stepping"            : { "time_step" : 1.0 },
            "rotation_dofs"            : false
        },
        "processes" : {}
    }""" % (echo_level, echo_level))


def SmoothRandomDensities(divisions: int, seed: int = 0, threshold: float = 0.5,
                          minimum_density: float = 0.2):
    """A smooth random density field on the element grid, in [rho_min, 1].

    Smoothed rather than white, because a per-element coin flip is not a
    structure and gives a generative model nothing to learn.

    Returns:
        (2 * divisions * divisions,) array in element order - two triangles
        per structured cell, both taking the cell's density.
    """
    generator = numpy.random.default_rng(seed)
    field = generator.standard_normal((divisions, divisions))
    for _ in range(3):
        field = (field
                 + numpy.roll(field, 1, 0) + numpy.roll(field, -1, 0)
                 + numpy.roll(field, 1, 1) + numpy.roll(field, -1, 1)) / 5.0
    field = (field - field.mean()) / (field.std() + 1e-12)
    densities = numpy.where(field > threshold, 1.0, minimum_density)
    return numpy.repeat(densities.reshape(-1), 2)   # two triangles per cell


def CreateComplianceModelPart(model: Kratos.Model, divisions: int = 8) -> Kratos.ModelPart:
    """The meshed plane-stress square, one Properties object per element."""
    import KratosMultiphysics.StructuralMechanicsApplication as SMA

    model_part = model.CreateModelPart("ComplianceModelPart")
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

    # the generator hands every element the SAME properties; SIMP needs one
    # modulus per element, so each gets its own (with its own law)
    for element in model_part.Elements:
        properties = model_part.CreateNewProperties(1000 + element.Id)
        properties.SetValue(Kratos.YOUNG_MODULUS, _YOUNG_MODULUS)
        properties.SetValue(Kratos.POISSON_RATIO, _POISSON_RATIO)
        properties.SetValue(Kratos.CONSTITUTIVE_LAW, SMA.LinearElasticPlaneStress2DLaw())
        element.Properties = properties
    return model_part


def ApplyDensities(model_part: Kratos.ModelPart, densities) -> None:
    """Writes a density field onto the elements through SIMP.

    E = E_min + rho^p (E_0 - E_min), the standard penalization: intermediate
    densities are made expensive so the design tends to black and white.
    """
    densities = numpy.asarray(densities, dtype=float).reshape(-1)
    if densities.size != model_part.NumberOfElements():
        raise ValueError(
            f"Got {densities.size} densities for {model_part.NumberOfElements()} elements.")
    for element, density in zip(model_part.Elements, densities):
        modulus = _MINIMUM_MODULUS + (density ** _PENALIZATION) * (
            _YOUNG_MODULUS - _MINIMUM_MODULUS)
        element.Properties.SetValue(Kratos.YOUNG_MODULUS, float(modulus))
        element.SetValue(Kratos.DENSITY, float(density))


def ApplyCaseData(model_part: Kratos.ModelPart, tip_load: float = _TIP_LOAD,
                  tolerance: float = 1e-8) -> None:
    """Clamps the left edge and hangs a point load at the right mid-height."""
    import KratosMultiphysics.StructuralMechanicsApplication as SMA

    for node in model_part.Nodes:
        if node.X0 < tolerance:
            node.Fix(Kratos.DISPLACEMENT_X)
            node.Fix(Kratos.DISPLACEMENT_Y)

    loaded = min((node for node in model_part.Nodes if node.X0 > _LENGTH - tolerance),
                 key=lambda node: abs(node.Y0 - 0.5 * _HEIGHT), default=None)
    if loaded is None:
        raise ValueError("No node on the loaded edge; the mesh is not what this case expects.")
    loaded.SetSolutionStepValue(SMA.POINT_LOAD, [0.0, tip_load, 0.0])
    # the nodal POINT_LOAD value alone does nothing: it enters the system
    # through a point-load CONDITION, without which the solve returns a zero
    # displacement field and a compliance of exactly zero
    domain = model_part.GetSubModelPart("Domain")
    domain.CreateNewCondition(
        "PointLoadCondition2D1N", model_part.NumberOfConditions() + 1,
        [loaded.Id], next(iter(model_part.Elements)).Properties)
    model_part.ProcessInfo[Kratos.STEP] = 0


def CreateComplianceAnalysis(model: Kratos.Model, densities=None, divisions: int = 8,
                             tip_load: float = _TIP_LOAD, echo_level: int = 0):
    """A ready-to-Run() analysis of one density design."""
    from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import (
        StructuralMechanicsAnalysis)

    model_part = CreateComplianceModelPart(model, divisions)
    if densities is None:
        densities = numpy.ones(model_part.NumberOfElements())
    ApplyDensities(model_part, densities)
    ApplyCaseData(model_part, tip_load)
    return StructuralMechanicsAnalysis(model, _CreateProjectParameters(echo_level))


def ComputeCompliance(model_part: Kratos.ModelPart, tolerance: float = 1e-8) -> float:
    """The structure's compliance, f.u over the loaded nodes.

    Lower is stiffer. With one point load this is just the load times the
    displacement under it, which is what a design is scored by.
    """
    import KratosMultiphysics.StructuralMechanicsApplication as SMA

    compliance = 0.0
    for node in model_part.Nodes:
        load = node.GetSolutionStepValue(SMA.POINT_LOAD)
        if abs(load[0]) < tolerance and abs(load[1]) < tolerance:
            continue
        displacement = node.GetSolutionStepValue(Kratos.DISPLACEMENT)
        compliance += sum(load[i] * displacement[i] for i in range(3))
    return float(compliance)


def SolveCompliance(model: Kratos.Model, densities, divisions: int = 8,
                    tip_load: float = _TIP_LOAD):
    """Solves one design and returns (compliance, model_part)."""
    analysis = CreateComplianceAnalysis(model, densities, divisions, tip_load)
    analysis.Run()
    model_part = model["ComplianceModelPart"]
    return ComputeCompliance(model_part), model_part


def DensityGrid(densities, divisions: int):
    """The element densities as the (divisions, divisions) image a grid model
    consumes - the two triangles of a cell share its density."""
    densities = numpy.asarray(densities, dtype=float).reshape(-1, 2)
    return densities[:, 0].reshape(divisions, divisions)


def ConstraintChannels(divisions: int, volume_fraction: float):
    """The conditioning a generative design model is given: where the
    structure is held, where it is loaded, and how much material it may use.

    The layout is DensityGrid's: the FIRST index runs along x and the second
    along y, so the clamped edge (x = 0, see ApplyCaseData) is a row of the
    first index and the loaded cell sits at the far x with y at mid-height.
    Writing them the other way round transposes the conditioning against the
    design image without changing either channel's sum, which is exactly the
    kind of mistake a sum-only test cannot see.

    Returns:
        (3, divisions, divisions) float array - support mask, load mask and
        a constant volume-fraction plane.
    """
    supports = numpy.zeros((divisions, divisions))
    supports[0, :] = 1.0                       # the clamped edge, x = 0
    loads = numpy.zeros((divisions, divisions))
    loads[-1, divisions // 2] = 1.0            # the loaded cell, far x at mid-height
    fraction = numpy.full((divisions, divisions), float(volume_fraction))
    return numpy.stack([supports, loads, fraction])

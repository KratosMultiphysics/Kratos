"""Real lid-driven cavity case built entirely in memory.

Unit-square domain, structured triangle mesh, no-slip on three walls and a
prescribed tangential velocity on the lid: the classic incompressible cavity
solved by FluidDynamicsApplication's monolithic VMS solver. The fluid
counterpart of thermal_case.py / structural_case.py - everything
parametrized, no .mdpa fixtures.

Two setup details the solver enforces and that are easy to get wrong:
"volume_model_part_name" and every "skin_parts" entry must be FULL model part
names, and the default "skin_parts" of [""] raises - an empty list must be
passed explicitly when no wall conditions are needed (velocity is fixed
nodally here instead).

Requires FluidDynamicsApplication (and LinearSolversApplication for the
default linear solver); callers gate on availability via
kratos_utilities.CheckIfApplicationsAvailable.
"""

import KratosMultiphysics as Kratos

# The complete variable set the monolithic solver's AddVariables adds - with
# use_input_model_part the nodes exist before AddVariables runs, so everything
# must be pre-added (DENSITY/VISCOSITY come from the vms element's nodal
# properties list).
_CORE_HISTORICAL_VARIABLES = (
    "VELOCITY", "ACCELERATION", "MESH_VELOCITY", "PRESSURE", "IS_STRUCTURE",
    "DISPLACEMENT", "BODY_FORCE", "NODAL_AREA", "NODAL_H", "ADVPROJ", "DIVPROJ",
    "REACTION", "REACTION_WATER_PRESSURE", "EXTERNAL_PRESSURE", "NORMAL",
    "Y_WALL", "DENSITY", "VISCOSITY",
)

# FluidDynamicsApplication's own variables, added by the same AddVariables call
_APP_HISTORICAL_VARIABLES = ("Q_VALUE",)


def _CreateProjectParameters(time_step=0.1, end_time=0.1, echo_level=0):
    return Kratos.Parameters("""{
        "problem_data": {
            "problem_name"  : "physics_nemo_fluid_case",
            "parallel_type" : "OpenMP",
            "start_time"    : 0.0,
            "end_time"      : %f,
            "echo_level"    : %d
        },
        "solver_settings": {
            "solver_type"             : "monolithic",
            "model_part_name"         : "FluidModelPart",
            "domain_size"             : 2,
            "model_import_settings"   : { "input_type" : "use_input_model_part" },
            "material_import_settings": { "materials_filename" : "" },
            "echo_level"              : %d,
            "formulation"             : { "element_type" : "vms" },
            "volume_model_part_name"  : "FluidModelPart.Domain",
            "skin_parts"              : [],
            "no_skin_parts"           : [],
            "time_stepping"           : { "automatic_time_step" : false, "time_step" : %f },
            "linear_solver_settings"  : { "solver_type" : "skyline_lu_factorization" },
            "maximum_iterations"      : 10,
            "relative_velocity_tolerance" : 1e-5,
            "absolute_velocity_tolerance" : 1e-7,
            "relative_pressure_tolerance" : 1e-5,
            "absolute_pressure_tolerance" : 1e-7
        },
        "processes"        : {},
        "output_processes" : {}
    }""" % (end_time, echo_level, echo_level, time_step))


def CreateFluidModelPart(model: Kratos.Model, divisions: int = 8) -> Kratos.ModelPart:
    """Creates the meshed cavity model part (variables added before the mesh)."""
    import KratosMultiphysics.FluidDynamicsApplication as KratosCFD

    model_part = model.CreateModelPart("FluidModelPart")
    model_part.ProcessInfo[Kratos.DOMAIN_SIZE] = 2
    model_part.SetBufferSize(3)  # the fluid solver's minimum
    for name in _CORE_HISTORICAL_VARIABLES:
        model_part.AddNodalSolutionStepVariable(Kratos.KratosGlobals.GetVariable(name))
    for name in _APP_HISTORICAL_VARIABLES:
        model_part.AddNodalSolutionStepVariable(getattr(KratosCFD, name))

    generator_geometry = Kratos.Quadrilateral2D4(
        Kratos.Node(1, 0.0, 0.0, 0.0), Kratos.Node(2, 0.0, 1.0, 0.0),
        Kratos.Node(3, 1.0, 1.0, 0.0), Kratos.Node(4, 1.0, 0.0, 0.0))
    mesh_parameters = Kratos.Parameters("""{
        "number_of_divisions"        : %d,
        "element_name"               : "Element2D3N",
        "condition_name"             : "LineCondition",
        "create_skin_sub_model_part" : false
    }""" % divisions)
    domain = model_part.CreateSubModelPart("Domain")
    Kratos.StructuredMeshGeneratorProcess(generator_geometry, domain, mesh_parameters).Execute()

    # the solver replaces Element2D3N with VMS2D3N and reads these off the
    # element properties to fill the nodal DENSITY/VISCOSITY
    properties = next(iter(domain.Elements)).Properties
    properties.SetValue(Kratos.DENSITY, 1.0)
    properties.SetValue(Kratos.DYNAMIC_VISCOSITY, 0.01)
    return model_part


def ApplyCaseData(model_part: Kratos.ModelPart, lid_velocity: float,
                  tolerance: float = 1e-8) -> None:
    """No-slip on three walls, a tangential lid, and one pinned pressure DOF."""
    for node in model_part.Nodes:
        on_lid = abs(node.Y - 1.0) < tolerance
        on_wall = (abs(node.X) < tolerance or abs(node.X - 1.0) < tolerance
                   or abs(node.Y) < tolerance)
        if on_lid or on_wall:
            node.Fix(Kratos.VELOCITY_X)
            node.Fix(Kratos.VELOCITY_Y)
            # the lid's corners belong to the walls: no-slip wins there
            velocity_x = lid_velocity if (on_lid and not on_wall) else 0.0
            node.SetSolutionStepValue(Kratos.VELOCITY, [velocity_x, 0.0, 0.0])

    # incompressible flow with velocity fixed everywhere on the boundary
    # determines pressure only up to a constant: pin one node
    corner = min(model_part.Nodes, key=lambda node: (node.X**2 + node.Y**2))
    corner.Fix(Kratos.PRESSURE)
    corner.SetSolutionStepValue(Kratos.PRESSURE, 0.0)


def CreateFluidAnalysis(model: Kratos.Model,
                        lid_velocity: float = 1.0,
                        divisions: int = 8,
                        time_step: float = 0.1,
                        end_time: float = 0.1,
                        echo_level: int = 0):
    """Builds the meshed cavity and returns a ready FluidDynamicsAnalysis.

    Call .Run() on the result (the default end_time runs a single step); the
    solved VELOCITY/PRESSURE fields live on model["FluidModelPart"].
    """
    from KratosMultiphysics.FluidDynamicsApplication.fluid_dynamics_analysis import FluidDynamicsAnalysis

    model_part = CreateFluidModelPart(model, divisions)
    ApplyCaseData(model_part, lid_velocity)
    return FluidDynamicsAnalysis(
        model, _CreateProjectParameters(time_step, end_time, echo_level))

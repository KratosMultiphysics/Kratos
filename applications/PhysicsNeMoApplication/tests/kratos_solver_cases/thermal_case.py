"""Real stationary heat-conduction case built entirely in memory.

Unit-square domain, structured triangle mesh, TEMPERATURE = 0 on the whole
boundary, uniform conductivity k and volumetric source f: the classic
-k lap(u) = f Poisson problem solved by ConvectionDiffusionApplication's
stationary solver. Used by the real-solver integration tests and the
active-learning template — everything parametrized, no .mdpa fixtures.

Requires ConvectionDiffusionApplication (and LinearSolversApplication for the
default linear solver); callers gate on availability via
kratos_utilities.CheckIfApplicationsAvailable.
"""

import KratosMultiphysics as Kratos

# The complete variable set the stationary solver's AddVariables adds (its
# default convection_diffusion_variables block) — with use_input_model_part
# the nodes exist before AddVariables runs, so everything must be pre-added.
_HISTORICAL_VARIABLES = (
    "TEMPERATURE", "DENSITY", "SPECIFIC_HEAT", "CONDUCTIVITY", "HEAT_FLUX",
    "FACE_HEAT_FLUX", "PROJECTED_SCALAR1", "CONVECTION_VELOCITY",
    "TEMPERATURE_GRADIENT", "TRANSFER_COEFFICIENT", "REACTION",
    "VELOCITY", "MESH_VELOCITY", "REACTION_FLUX",
)


def _CreateProjectParameters(echo_level=0, transient_settings=None):
    """Stationary by default; transient_settings (dict with time_step,
    end_time, theta, dynamic_tau) switches to the transient solver."""
    parameters = Kratos.Parameters("""{
        "problem_data": {
            "problem_name"  : "physics_nemo_thermal_case",
            "parallel_type" : "OpenMP",
            "start_time"    : 0.0,
            "end_time"      : 0.99,
            "echo_level"    : %d
        },
        "solver_settings": {
            "solver_type"                        : "stationary",
            "analysis_type"                      : "linear",
            "model_part_name"                    : "ThermalModelPart",
            "domain_size"                        : 2,
            "model_import_settings"              : { "input_type" : "use_input_model_part" },
            "material_import_settings"           : { "materials_filename" : "" },
            "echo_level"                         : %d,
            "problem_domain_sub_model_part_list" : ["Domain"],
            "processes_sub_model_part_list"      : [],
            "convection_diffusion_variables"     : {
                "unknown_variable"       : "TEMPERATURE",
                "density_variable"       : "DENSITY",
                "specific_heat_variable" : "SPECIFIC_HEAT",
                "diffusion_variable"     : "CONDUCTIVITY",
                "volume_source_variable" : "HEAT_FLUX",
                "velocity_variable"      : "VELOCITY",
                "mesh_velocity_variable" : "MESH_VELOCITY",
                "reaction_variable"      : "REACTION_FLUX"
            },
            "time_stepping"                      : { "time_step" : 1.0 }
        },
        "processes"        : {},
        "output_processes" : {}
    }""" % (echo_level, echo_level))

    if transient_settings is None:
        return parameters

    solver_settings = parameters["solver_settings"]
    solver_settings["solver_type"].SetString("transient")
    solver_settings["time_stepping"]["time_step"].SetDouble(
        float(transient_settings.get("time_step", 0.05)))
    parameters["problem_data"]["end_time"].SetDouble(
        float(transient_settings.get("end_time", 0.5)))
    transient_parameters = solver_settings.AddEmptyValue("transient_parameters")
    transient_parameters.AddEmptyValue("theta").SetDouble(
        float(transient_settings.get("theta", 0.5)))
    transient_parameters.AddEmptyValue("dynamic_tau").SetDouble(
        float(transient_settings.get("dynamic_tau", 0.0)))
    # the solver reads every key of this block directly (no defaults merge)
    transient_parameters.AddEmptyValue("cross_wind_stabilization_factor").SetDouble(
        float(transient_settings.get("cross_wind_stabilization_factor", 0.0)))
    return parameters


def CreateThermalModelPart(model: Kratos.Model, divisions: int = 10) -> Kratos.ModelPart:
    """Creates the meshed model part (variables added before the mesh)."""
    model_part = model.CreateModelPart("ThermalModelPart")
    model_part.ProcessInfo[Kratos.DOMAIN_SIZE] = 2
    model_part.SetBufferSize(2)
    for name in _HISTORICAL_VARIABLES:
        model_part.AddNodalSolutionStepVariable(Kratos.KratosGlobals.GetVariable(name))

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
    return model_part


def ApplyCaseData(model_part: Kratos.ModelPart, conductivity: float, heat_flux: float, tolerance: float = 1e-8) -> None:
    """Sets material data everywhere and TEMPERATURE = 0 on the boundary."""
    for node in model_part.Nodes:
        node.SetSolutionStepValue(Kratos.DENSITY, 1.0)
        node.SetSolutionStepValue(Kratos.SPECIFIC_HEAT, 1.0)
        node.SetSolutionStepValue(Kratos.CONDUCTIVITY, conductivity)
        node.SetSolutionStepValue(Kratos.HEAT_FLUX, heat_flux)
        on_boundary = (abs(node.X) < tolerance or abs(node.X - 1.0) < tolerance or
                       abs(node.Y) < tolerance or abs(node.Y - 1.0) < tolerance)
        if on_boundary:
            node.Fix(Kratos.TEMPERATURE)
            node.SetSolutionStepValue(Kratos.TEMPERATURE, 0.0)


def CreateThermalAnalysis(model: Kratos.Model,
                          conductivity: float = 1.0,
                          heat_flux: float = 1.0,
                          divisions: int = 10,
                          echo_level: int = 0):
    """Builds the meshed case and returns a ready ConvectionDiffusionAnalysis.

    Call .Run() on the result; the solved TEMPERATURE field lives on
    model["ThermalModelPart"].
    """
    from KratosMultiphysics.ConvectionDiffusionApplication.convection_diffusion_analysis import ConvectionDiffusionAnalysis

    model_part = CreateThermalModelPart(model, divisions)
    ApplyCaseData(model_part, conductivity, heat_flux)
    return ConvectionDiffusionAnalysis(model, _CreateProjectParameters(echo_level))


def CreateTransientThermalAnalysis(model: Kratos.Model,
                                   conductivity: float = 1.0,
                                   heat_flux: float = 1.0,
                                   divisions: int = 10,
                                   time_step: float = 0.05,
                                   end_time: float = 0.5,
                                   theta: float = 0.5,
                                   echo_level: int = 0):
    """The same case run by ConvectionDiffusion's TRANSIENT solver.

    Time integration is done by the ELEMENT (the solver installs the static
    scheme as a "fake" scheme and passes theta/dynamic_tau/DELTA_TIME
    through ProcessInfo), so residual assembly needs no special scheme -
    only the solution-step buffer the solver already maintains. The initial
    condition is TEMPERATURE = 0 everywhere.
    """
    from KratosMultiphysics.ConvectionDiffusionApplication.convection_diffusion_analysis import ConvectionDiffusionAnalysis

    model_part = CreateThermalModelPart(model, divisions)
    ApplyCaseData(model_part, conductivity, heat_flux)
    return ConvectionDiffusionAnalysis(model, _CreateProjectParameters(
        echo_level, {"time_step": time_step, "end_time": end_time, "theta": theta}))

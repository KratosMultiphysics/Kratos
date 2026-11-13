"""Process closing the residual-scoring -> mesh-adaptation loop.

At a configurable step interval, assembles the PDE residual of the model
part's CURRENT state (whatever a solver or a deployed surrogate last wrote)
through solver_residuals.ResidualEvaluator, turns the per-node residual into
a target size field, and remeshes with MMG. Pair it with any of the
deployment processes to refine where the surrogate's physics error
concentrates.

Requires MeshingApplication with MMG (checked lazily at first remesh).
"""

import KratosMultiphysics as Kratos

from KratosMultiphysics.PhysicsNeMoApplication.bridges.mesh_bridge import adaptive_remeshing
from KratosMultiphysics.PhysicsNeMoApplication.physics import solver_residuals
def Factory(settings: Kratos.Parameters, model: Kratos.Model):
    if not isinstance(settings, Kratos.Parameters):
        raise TypeError("expected input shall be a Parameters object, encapsulating a json string")
    return AdaptiveRemeshProcess(model, settings["Parameters"])


class AdaptiveRemeshProcess(Kratos.Process):
    """Residual-driven MMG adaptation at a step interval.

    Settings:
        {
            "model_part_name" : "",
            "remesh_interval" : 1,
            "size_settings"   : { "target_error": 1e-3, "exponent": 0.5,
                                  "minimal_size": 1e-4, "maximal_size": 1.0 },
            "mmg_parameters"  : { },
            "echo_level"      : 0
        }

    "mmg_parameters" is passed through to MmgProcess (merged over the
    module defaults in adaptive_remeshing.RunMmgAdaptation).
    """

    def __init__(self, model: Kratos.Model, settings: Kratos.Parameters):
        super().__init__()

        # free-form pass-through block: split off before validation
        self.mmg_parameters = None
        if settings.Has("mmg_parameters"):
            self.mmg_parameters = settings["mmg_parameters"].Clone()
            settings.RemoveValue("mmg_parameters")

        defaults = Kratos.Parameters("""{
            "model_part_name" : "",
            "remesh_interval" : 1,
            "size_settings"   : {},
            "echo_level"      : 0
        }""")
        settings.ValidateAndAssignDefaults(defaults)

        self.model_part = model[settings["model_part_name"].GetString()]
        self.remesh_interval = settings["remesh_interval"].GetInt()
        if self.remesh_interval < 1:
            raise ValueError(f"remesh_interval must be >= 1, got {self.remesh_interval}.")
        self.size_settings = settings["size_settings"].Clone()
        self.echo_level = settings["echo_level"].GetInt()

    def ExecuteFinalizeSolutionStep(self):
        step = self.model_part.ProcessInfo[Kratos.STEP]
        if step % self.remesh_interval != 0:
            return
        self.Execute()

    def Execute(self):
        # The DOF set and the mesh both change across remeshes: build a
        # fresh evaluator every time.
        evaluator = solver_residuals.BuildResidualEvaluator(self.model_part)
        nodal_error = adaptive_remeshing.NodalErrorArray(
            self.model_part, evaluator.ComputeNodalResiduals())

        size_field = adaptive_remeshing.ComputeTargetSizeField(
            self.model_part, nodal_error, self.size_settings.Clone())
        nodes_before = self.model_part.NumberOfNodes()
        adaptive_remeshing.RunMmgAdaptation(self.model_part, size_field, self.mmg_parameters)

        if self.echo_level > 0:
            Kratos.Logger.PrintInfo(
                "AdaptiveRemeshProcess",
                f"remeshed at step {self.model_part.ProcessInfo[Kratos.STEP]}: "
                f"{nodes_before} -> {self.model_part.NumberOfNodes()} nodes")

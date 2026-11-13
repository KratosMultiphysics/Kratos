"""Process running a pure-PINN solve on a model part's nodes.

A physics-informed forward solve without training data: a coordinate MLP
(physicsnemo.models.mlp.fully_connected.FullyConnected) is trained against
the differentiable PDE residual (physics_informed / physicsnemo.sym's
PhysicsInformer, autodiff at the node coordinates plus optional random
collocation points) and the Dirichlet data of the model part (fixed DOFs
keep their current solution-step values as boundary targets). The
converged field is written into the output fields at
ExecuteBeforeSolutionLoop (or on demand via Solve()).

"inverse" mode recovers PDE coefficients from observations instead: the
coefficients named in "inverse_parameters" become trainable scalars fed to
the PDE as spatial inputs, the observation fields provide the (fixed)
field values, and detach_names blocks gradients through the observed
fields - upstream's documented inverse-problem mechanism. The recovered
values are exposed as process.inverse_values.

SystemIdentificationApplication's adjoint-based identification is the
classical baseline for the inverse mode.

torch/physicsnemo(.sym) are imported lazily at solve time.
"""

import numpy

import KratosMultiphysics as Kratos
from KratosMultiphysics.PhysicsNeMoApplication.deployment import model_registry
from KratosMultiphysics.PhysicsNeMoApplication.physics import physics_informed
from KratosMultiphysics.PhysicsNeMoApplication.bridges import torch_bridge
from KratosMultiphysics.PhysicsNeMoApplication.processes.inference.inference_process import WriteOutputFields
from KratosMultiphysics.PhysicsNeMoApplication.processes.inference.point_cloud_inference_process import GatherPointCloudCoordinates

_MODES = ("forward", "inverse")


def Factory(settings: Kratos.Parameters, model: Kratos.Model) -> "PinnSolveProcess":
    if not isinstance(settings, Kratos.Parameters):
        raise TypeError("Expected input shall be a Parameters object, encapsulating a json string")
    return PinnSolveProcess(model, settings["Parameters"])


class _NormalizingNetwork:
    """Min-max normalizes the network's INPUT, leaving derivatives physical.

    Built lazily so this module keeps importing without torch. See
    _Build for the actual nn.Module.
    """

    def __new__(cls, network, model_part, device):
        return _BuildNormalizingNetwork(network, model_part, device)


def _BuildNormalizingNetwork(network, model_part, device):
    import torch

    coordinates = GatherPointCloudCoordinates(model_part, normalize=False)
    low = coordinates.min(axis=0)
    extent = coordinates.max(axis=0) - low
    extent[extent == 0.0] = 1.0          # degenerate axes stay at 0

    class NormalizedInput(torch.nn.Module):
        """Applies (x - low) / extent before the network.

        Normalizing here rather than in the coordinates handed to
        PhysicsInformer is the whole point: the network still sees
        well-conditioned inputs in [0, 1], but autodiff is taken with
        respect to PHYSICAL coordinates, so du/dx is a physical derivative
        and the residual is the physical PDE.
        """

        def __init__(self, wrapped):
            super().__init__()
            self.wrapped = wrapped
            self.register_buffer("low", torch.tensor(low, dtype=torch.float32))
            self.register_buffer("extent", torch.tensor(extent, dtype=torch.float32))

        def forward(self, points):
            return self.wrapped((points - self.low) / self.extent)

    return NormalizedInput(network).to(device)


class PinnSolveProcess(Kratos.Process):
    """Trains a coordinate network against the PDE residual and writes the field."""

    def __init__(self, model: Kratos.Model, settings: Kratos.Parameters) -> None:
        super().__init__()

        # free-form: {"coefficient_name": initial_value}
        inverse_parameters = Kratos.Parameters("{}")
        if settings.Has("inverse_parameters"):
            inverse_parameters = settings["inverse_parameters"].Clone()
            settings.RemoveValue("inverse_parameters")
        self.inverse_parameter_init = {
            name: inverse_parameters[name].GetDouble() for name in inverse_parameters.keys()}

        default_settings = Kratos.Parameters("""{
            "model_part_name"    : "PLEASE_SPECIFY_MODEL_PART_NAME",
            "mode"               : "forward",
            "physics"            : {},
            "fields"             : [ { "name" : "u", "width" : 1 } ],
            "solution_fields"    : [
                {
                    "variable_name" : "TEMPERATURE",
                    "data_location" : "node_historical"
                }
            ],
            "observation_fields" : [],
            "output_fields"      : [
                {
                    "variable_name" : "TEMPERATURE",
                    "data_location" : "node_historical"
                }
            ],
            "network"            : {
                "layer_size"    : 64,
                "num_layers"    : 4,
                "activation_fn" : "silu"
            },
            "training"           : {
                "epochs"             : 500,
                "learning_rate"      : 1e-3,
                "collocation_points" : 0,
                "physics_weight"     : 1.0,
                "boundary_weight"    : 10.0,
                "data_weight"        : 1.0,
                "echo_interval"      : 0,
                "seed"               : 0
            },
            "device"             : "auto",
            "normalize_coordinates" : true
        }""")
        settings.ValidateAndAssignDefaults(default_settings)
        settings["network"].ValidateAndAssignDefaults(default_settings["network"])
        settings["training"].ValidateAndAssignDefaults(default_settings["training"])
        field_defaults = default_settings["solution_fields"][0]
        for key in ("solution_fields", "observation_fields", "output_fields"):
            for i in range(settings[key].size()):
                settings[key][i].ValidateAndAssignDefaults(field_defaults)

        self.model_part = model[settings["model_part_name"].GetString()]
        self.mode = settings["mode"].GetString()
        if self.mode not in _MODES:
            raise ValueError(f"Unsupported mode \"{self.mode}\". Use one of {_MODES}.")
        if self.mode == "inverse" and not self.inverse_parameter_init:
            raise ValueError("\"inverse\" mode needs a non-empty \"inverse_parameters\" block.")

        self.physics_settings = settings["physics"].Clone()
        if not self.physics_settings.Has("grad_method"):
            self.physics_settings.AddEmptyValue("grad_method").SetString("autodiff")
        elif self.physics_settings["grad_method"].GetString() != "autodiff":
            raise ValueError("PinnSolveProcess supports grad_method \"autodiff\" only.")

        def read_specs(key):
            return [(settings[key][i]["variable_name"].GetString(),
                     settings[key][i]["data_location"].GetString())
                    for i in range(settings[key].size())]

        self.solution_specs = read_specs("solution_fields")
        self.observation_specs = read_specs("observation_fields")
        self.output_specs = read_specs("output_fields")
        if self.mode == "inverse" and not self.observation_specs:
            raise ValueError("\"inverse\" mode needs \"observation_fields\".")

        # (name, width, components) triples - shared format with the
        # physics_informed loss terms (components name the sympy functions)
        self.field_specs = physics_informed._ReadFieldSpecs(settings)

        self.network_settings = settings["network"].Clone()
        self.training_settings = settings["training"].Clone()
        self.device_name = settings["device"].GetString()
        self.normalize_coordinates = settings["normalize_coordinates"].GetBool()

        self.loss_history = []
        self.inverse_values = {}
        self._network = None

    def ExecuteBeforeSolutionLoop(self) -> None:
        self.Solve()

    # --- assembly helpers ----------------------------------------------------

    def _GatherFieldMatrix(self, specs):
        """(N, total width) float64 numpy over the model part's nodes.

        Blocks are cut to the PDE FIELD widths, not the Kratos variable
        widths: a 2-D velocity is a 3-component Kratos array but only two
        PDE functions, and the network predicts two. Gathering all three
        would leave the targets a column wider than the prediction.
        """
        from KratosMultiphysics.PhysicsNeMoApplication.processes.inference.inference_process import GatherInputFields
        blocks, _ = GatherInputFields(self.model_part, specs)
        import torch
        widths = self._FieldWidths(len(blocks))
        blocks = [block[:, :width] for block, width in zip(blocks, widths)]
        return torch.cat(blocks, dim=-1).numpy()

    def _ExpandToKratosWidths(self, solution):
        """The network's PDE-width output padded back to Kratos widths.

        The mirror of the cut _GatherFieldMatrix makes on the way in: a 2-D
        velocity is two PDE functions but a three-component Kratos array, so
        the missing out-of-plane column is written as zero rather than left
        for the writer to trip over.
        """
        torch = torch_bridge._TryImportTorch()

        widths = self._FieldWidths(len(self.output_specs))
        blocks = []
        offset = 0
        for (variable_name, _), width in zip(self.output_specs, widths):
            block = solution[:, offset:offset + width]
            offset += width
            variable = Kratos.KratosGlobals.GetVariable(variable_name)
            kratos_width = 1 if isinstance(variable, Kratos.DoubleVariable) else 3
            if block.shape[1] < kratos_width:
                padding = torch.zeros(
                    (block.shape[0], kratos_width - block.shape[1]), dtype=block.dtype)
                block = torch.cat([block, padding], dim=1)
            blocks.append(block)
        return torch.cat(blocks, dim=1)

    def _FieldWidths(self, count: int):
        """The PDE width of each solution field, in order."""
        if len(self.field_specs) != count:
            raise ValueError(
                f"\"fields\" declares {len(self.field_specs)} field(s) but "
                f"{count} solution field(s) were given; they are paired in order, "
                "so a 2-D velocity needs {\"name\": \"velocity\", \"width\": 2} "
                "against the VELOCITY variable.")
        return [width for _, width, _ in self.field_specs]

    def _DirichletMask(self):
        """(N, total width) bool: which solution DOFs are fixed.

        Cut to the PDE field widths for the same reason _GatherFieldMatrix
        is: a 2-D case fixes VELOCITY_X and VELOCITY_Y, and a third column
        would not line up with anything the network predicts.
        """
        columns = []
        widths = self._FieldWidths(len(self.solution_specs))
        for (variable_name, _), width in zip(self.solution_specs, widths):
            variable = Kratos.KratosGlobals.GetVariable(variable_name)
            if isinstance(variable, Kratos.DoubleVariable):
                components = [variable]
            else:  # array variable: component-wise fixity, to the PDE's width
                components = [Kratos.KratosGlobals.GetVariable(f"{variable_name}_{axis}")
                              for axis in "XYZ"[:width]]
            block = numpy.zeros((self.model_part.NumberOfNodes(), len(components)), dtype=bool)
            for row, node in enumerate(self.model_part.Nodes):
                for column, component in enumerate(components):
                    block[row, column] = node.IsFixed(component)
            columns.append(block)
        return numpy.concatenate(columns, axis=1)

    # --- the solve -----------------------------------------------------------

    def Solve(self) -> None:
        torch = torch_bridge._TryImportTorch()

        training = self.training_settings
        seed = training["seed"].GetInt()
        if seed >= 0:
            torch.manual_seed(seed)
            numpy.random.seed(seed if seed > 0 else None)
        device = model_registry.ResolveDevice(self.device_name)

        total_width = sum(width for _, width, _ in self.field_specs)
        from physicsnemo.models.mlp.fully_connected import FullyConnected
        network = FullyConnected(
            in_features=3, out_features=total_width,
            layer_size=self.network_settings["layer_size"].GetInt(),
            num_layers=self.network_settings["num_layers"].GetInt(),
            activation_fn=self.network_settings["activation_fn"].GetString()).to(device)
        if self.normalize_coordinates:
            network = _NormalizingNetwork(
                network, self.model_part, device).to(device)
        self._network = network

        # ALWAYS gather physical coordinates. normalize_coordinates conditions
        # the NETWORK's inputs (below); it must not touch the coordinates the
        # residual is differentiated against, or the PDE being enforced becomes
        # sum (1/L_i^2) d2u/dx_i^2 - a different equation, and a per-axis
        # different one on an anisotropic domain. A unit cube hides this
        # entirely because every L_i is 1.
        coordinates = GatherPointCloudCoordinates(self.model_part, normalize=False)
        n_collocation = training["collocation_points"].GetInt()
        if n_collocation > 0:
            low = coordinates.min(axis=0)
            high = coordinates.max(axis=0)
            extra = numpy.random.uniform(low, high, size=(n_collocation, 3))
            all_points = numpy.concatenate([coordinates, extra])
        else:
            all_points = coordinates
        points = torch.tensor(all_points, dtype=torch.float32, device=device)
        n_nodes = coordinates.shape[0]

        physics = self.physics_settings.Clone()
        parameters = list(network.parameters())
        inverse_tensors = {}
        if self.mode == "inverse":
            # coefficients become spatial PDE inputs fed from trainable
            # scalars; the field network trains on the observation data only
            # - detach_names blocks the physics gradients through the field
            # values AND their derivatives (upstream's documented mechanism)
            for name, initial in self.inverse_parameter_init.items():
                inverse_tensors[name] = torch.tensor(
                    float(initial), device=device, requires_grad=True)
            parameters += list(inverse_tensors.values())
            detach = set(physics["detach_names"].GetStringArray()) if physics.Has("detach_names") else set()
            for _, _, components in self.field_specs:  # components = informer input names
                for component in components:
                    detach.add(component)
                    for a in "xyz":
                        detach.add(f"{component}__{a}")
                        detach.add(f"{component}__{a}__{a}")
            if physics.Has("detach_names"):
                physics.RemoveValue("detach_names")
            physics.AddEmptyArray("detach_names")
            for name in sorted(detach):
                physics["detach_names"].Append(name)
        informer, residual_names, _ = physics_informed.CreatePhysicsInformer(physics)

        boundary_mask = self._DirichletMask()
        boundary_targets = torch.tensor(
            self._GatherFieldMatrix(self.solution_specs)[boundary_mask],
            dtype=torch.float32, device=device)
        mask = torch.from_numpy(boundary_mask).to(device)
        observations = (torch.tensor(self._GatherFieldMatrix(self.observation_specs),
                                     dtype=torch.float32, device=device)
                        if self.observation_specs else None)

        physics_weight = training["physics_weight"].GetDouble()
        boundary_weight = training["boundary_weight"].GetDouble()
        data_weight = training["data_weight"].GetDouble()
        echo_interval = training["echo_interval"].GetInt()
        optimizer = torch.optim.Adam(parameters, lr=training["learning_rate"].GetDouble())

        self.loss_history = []
        for epoch in range(training["epochs"].GetInt()):
            optimizer.zero_grad()
            # the (1, 3, P) coordinates leaf must be what the network input
            # derives from, or the autodiff graph breaks at the transpose
            solve_points = points.detach().clone().T[None]
            solve_points.requires_grad_(True)
            prediction = network(solve_points[0].T)

            informer_inputs = physics_informed._SplitFields(prediction, self.field_specs)
            informer_inputs["coordinates"] = solve_points
            for name, tensor in inverse_tensors.items():
                informer_inputs[name] = tensor.expand(1, 1, prediction.shape[0])
            residuals = informer.forward(informer_inputs)
            loss = None
            for name in residual_names:
                contribution = residuals[name].square().mean()
                loss = contribution if loss is None else loss + contribution
            loss = physics_weight * loss

            nodal_prediction = prediction[:n_nodes]
            if self.mode == "forward" and bool(mask.any()):
                loss = loss + boundary_weight * (
                    nodal_prediction[mask] - boundary_targets).square().mean()
            if observations is not None and data_weight > 0.0:
                loss = loss + data_weight * (nodal_prediction - observations).square().mean()

            loss.backward()
            optimizer.step()
            loss_value = loss.item()
            self.loss_history.append(loss_value)
            if echo_interval > 0 and (epoch + 1) % echo_interval == 0:
                Kratos.Logger.PrintInfo(
                    type(self).__name__,
                    f"epoch {epoch + 1}/{training['epochs'].GetInt()}: loss = {loss_value:.6e}")

        self.inverse_values = {name: float(tensor) for name, tensor in inverse_tensors.items()}
        if self.mode == "inverse":
            Kratos.Logger.PrintInfo(
                type(self).__name__, f"recovered coefficients: {self.inverse_values}")
            return

        with torch.no_grad():
            solution = network(points[:n_nodes]).cpu().to(torch.float64)
        WriteOutputFields(self.model_part, self.output_specs,
                          self._ExpandToKratosWidths(solution), n_nodes)

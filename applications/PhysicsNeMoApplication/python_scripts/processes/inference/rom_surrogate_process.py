"""Process deploying a neural ROM surrogate: parameters in, POD field out.

Neural-augmented reduced bases on top of RomApplication's POD output: a
trained model maps the case parameters to the reduced coordinates q, and
the full-order field is reconstructed as u = phi q through rom_bridge and
written onto the basis's nodal unknowns - a full-field prediction at the
cost of an n_modes-sized network.

The parameters travel as ordinary input fields (any constant nodal carrier
works): they are gathered like InferenceProcess inputs and MEAN-reduced over
the entities to one (1, C_in) vector; the model returns (1, n_modes).
``output_fields`` is derived from the basis's ``nodal_unknowns`` (never set
by the user), so the advisory model-card check validates the ROM contract
automatically.

torch is imported lazily on first execution; the basis loads with numpy at
construction. Producing a basis needs the compiled RomApplication - see
rom_bridge.
"""

import KratosMultiphysics as Kratos
from KratosMultiphysics.PhysicsNeMoApplication.bridges import rom_bridge
from KratosMultiphysics.PhysicsNeMoApplication.bridges import torch_bridge
from KratosMultiphysics.PhysicsNeMoApplication.processes.inference.inference_process import InferenceProcess


def Factory(settings: Kratos.Parameters, model: Kratos.Model) -> "RomSurrogateProcess":
    if not isinstance(settings, Kratos.Parameters):
        raise TypeError("Expected input shall be a Parameters object, encapsulating a json string")
    return RomSurrogateProcess(model, settings["Parameters"])


class RomSurrogateProcess(InferenceProcess):
    """Runs the parameter -> q -> phi q reconstruction each due step."""

    def __init__(self, model: Kratos.Model, settings: Kratos.Parameters) -> None:
        # Split the subclass keys off (with their defaults) before the parent
        # validates the shared InferenceProcess settings.
        rom_basis_folder = "rom_data"
        rom_basis_name = "RomParameters"
        input_reduction = "mean"
        if settings.Has("rom_basis_folder"):
            rom_basis_folder = settings["rom_basis_folder"].GetString()
            settings.RemoveValue("rom_basis_folder")
        if settings.Has("rom_basis_name"):
            rom_basis_name = settings["rom_basis_name"].GetString()
            settings.RemoveValue("rom_basis_name")
        if settings.Has("input_reduction"):
            input_reduction = settings["input_reduction"].GetString()
            settings.RemoveValue("input_reduction")
        if input_reduction != "mean":
            raise ValueError(
                f"Unsupported input reduction \"{input_reduction}\". Only \"mean\" is supported.")

        if settings.Has("output_fields"):
            raise ValueError(
                "RomSurrogateProcess derives \"output_fields\" from the basis's nodal_unknowns; "
                "remove the setting.")

        self.rom_basis = rom_bridge.LoadRomBasis(rom_basis_folder, rom_basis_name)

        # Synthesize output_fields = the basis's nodal unknowns, so the parent
        # (and the advisory model card check) validate the ROM contract.
        settings.AddEmptyArray("output_fields")
        for name in self.rom_basis.nodal_unknowns:
            entry = Kratos.Parameters("""{ "variable_name": "", "data_location": "node_historical" }""")
            entry["variable_name"].SetString(name)
            settings["output_fields"].Append(entry)

        super().__init__(model, settings)
        self._permutation = None
        self._permutation_nodes = None

    def RunInference(self) -> None:
        model = self._GetModel()
        torch = torch_bridge._TryImportTorch()

        features, _ = self._GatherFeatures()
        features = features.mean(dim=0, keepdim=True)  # (1, C_in)
        parameter = next(model.parameters(), None)
        if parameter is not None:
            features = features.to(parameter.dtype)

        with torch.no_grad():
            q = model(features.to(self._device)).cpu()

        if tuple(q.shape) != (1, self.rom_basis.n_modes):
            raise ValueError(
                f"The model must return a (1, {self.rom_basis.n_modes}) reduced-coordinates "
                f"prediction (the basis has {self.rom_basis.n_modes} mode(s)); got "
                f"{list(q.shape)}.")

        u = rom_bridge.ReconstructFromReducedSpace(
            self.rom_basis, q[0].to(torch.float64).numpy())
        # the basis row ordering owns the write-back (not the parent's
        # _WriteOutputs, which knows nothing about the interleaving); the
        # permutation is topology, built once and reused until the node
        # count changes (AdaptiveRemeshProcess ships in this application)
        n_nodes = self.model_part.NumberOfNodes()
        if self._permutation is None or self._permutation_nodes != n_nodes:
            self._permutation = rom_bridge.NodePermutation(self.model_part, self.rom_basis)
            self._permutation_nodes = n_nodes
        rom_bridge.ScatterUnknownsVector(self.model_part, self.rom_basis, u,
                                         permutation=self._permutation)

import pathlib

import numpy

import KratosMultiphysics
from KratosMultiphysics.RomApplication.future.rom_decoders import GetBasisRows, GetRomParameters, LinearDecoder, RegisterDecoder, RomDecoder


class LocalRomDecoder(RomDecoder):
    """Decoder made of several local models (each a RomDecoder), one of them active at a time.

    Decode and Jacobian are the ones of the active model, so the solver and the projector work as with a single decoder. At the start
    of each step (UpdateReducedCoordinates) the selector picks the model of the step; on a switch, the current state is encoded in
    the new model (q = E_k(u)), whose number of reduced coordinates may differ. With one model it behaves as that model.

    The state of the first step is always encoded in the model chosen for it, so a model whose origin is not the initial state
    starts from it too.
    """

    def __init__(self, models, selector):
        if len(models) == 0:
            raise Exception("A LocalRomDecoder needs at least one model.")
        self.models = list(models)
        self.selector = selector
        self.active_model = 0
        self.is_initialized = False
        # Model of each step, in the order they were solved
        self.model_history = []

    @classmethod
    def FromRomData(cls, rom_data_folder, dof_set):
        """From the "local_roms" of the RomParameters.json: one linear decoder per local model, each read from its own subfolder
        ("models") with the files of a single one (RightBasisMatrix.npy, and optionally ReferenceSnapshot.npy as its origin, else the
        initial state), with the "number_of_rom_dofs" of each model (else the one of the "rom_settings"), and the "selector" of the
        model of each step (see CreateModelSelector). The rows of the bases are given by the NodeIds.npy of rom_data_folder."""
        folder = pathlib.Path(rom_data_folder)
        rom_parameters = GetRomParameters(folder)
        local_roms = rom_parameters["local_roms"]
        rows = GetBasisRows(folder, dof_set)
        u0 = numpy.array(dof_set.GetValues())
        models = local_roms["models"]
        number_of_rom_dofs = local_roms.get("number_of_rom_dofs", [rom_parameters["rom_settings"]["number_of_rom_dofs"]] * len(models))
        decoders = []
        for name, n in zip(models, number_of_rom_dofs):
            right_modes = numpy.load(folder / name / "RightBasisMatrix.npy")
            if right_modes.ndim == 1:
                right_modes = right_modes.reshape(-1, 1)
            origin_file = folder / name / "ReferenceSnapshot.npy"
            origin = numpy.load(origin_file)[rows] if origin_file.exists() else u0
            decoders.append(LinearDecoder(right_modes[rows, :n], origin))
        return cls(decoders, CreateModelSelector(local_roms["selector"], folder, rows))

    def WithPrescribedFixedDofs(self, dof_set):
        """Each model gives the prescribed values on the fixed DOFs, so that a switch encodes the state with them."""
        self.models = [model.WithPrescribedFixedDofs(dof_set) for model in self.models]
        return self

    def GetActiveModel(self):
        return self.active_model

    def NumberOfRomDofs(self):
        return self.models[self.active_model].NumberOfRomDofs()

    def Decode(self, q):
        return self.models[self.active_model].Decode(q)

    def Jacobian(self, q):
        return self.models[self.active_model].Jacobian(q)

    def Encode(self, u, **kwargs):
        return self.models[self.active_model].Encode(u, **kwargs)

    def InitializeSolutionStep(self):
        # All of them, so that any model can take over with the conditions of the step (e.g. a FixedDofsDecoder)
        for model in self.models:
            model.InitializeSolutionStep()

    def UpdateReducedCoordinates(self, q, dof_set, process_info):
        u = numpy.array(dof_set.GetValues())
        model = self.selector.Select(self, u, process_info)
        if not 0 <= model < len(self.models):
            raise Exception("The selector chose model {}, but there are {} models.".format(model, len(self.models)))
        if model != self.active_model or not self.is_initialized:
            self.active_model = model
            q = self.models[model].Encode(u)
            self.is_initialized = True
        self.model_history.append(self.active_model)
        return q


class RomModelSelector:
    """Base class of the selectors of a LocalRomDecoder: Select returns the model of the step."""

    def Select(self, local_decoder, u, process_info):
        """Model of the step, given the decoder (its active model and models), the current state u (effective DOF set order) and the
        process info of the step."""
        raise NotImplementedError("Calling base class 'Select'.")


class ConstantSelector(RomModelSelector):
    """Always the same model."""

    def __init__(self, model=0):
        self.model = model

    def Select(self, local_decoder, u, process_info):
        return self.model


class ProcessInfoSelector(RomModelSelector):
    """Parametric selection: the model is given by a variable of the process info of the step (e.g. TIME, or a parameter),
    compared with increasing thresholds: model k for thresholds[k-1] <= value < thresholds[k]."""

    def __init__(self, variable, thresholds):
        self.variable = KratosMultiphysics.KratosGlobals.GetVariable(variable) if isinstance(variable, str) else variable
        self.thresholds = numpy.array(thresholds, dtype=numpy.float64)
        if numpy.any(numpy.diff(self.thresholds) <= 0.0):
            raise Exception("The thresholds of a ProcessInfoSelector must be increasing.")

    def Select(self, local_decoder, u, process_info):
        return int(numpy.searchsorted(self.thresholds, process_info[self.variable], side="right"))


class NearestCentreSelector(RomModelSelector):
    """State-based selection: the model whose centre (e.g. the mean of its snapshots) is closest to the current state.
    The active model is kept unless another centre is closer by more than a relative 'margin' (hysteresis), so that the
    selection does not alternate at the boundary between two models. centres: one column per model, rows in the effective DOF set order."""

    def __init__(self, centres, margin=0.0):
        self.centres = numpy.asarray(centres, dtype=numpy.float64)
        if self.centres.ndim != 2:
            raise Exception("The centres must be a matrix with one row per DOF and one column per model.")
        self.margin = float(margin)
        self.distances = []

    def Select(self, local_decoder, u, process_info):
        distances = numpy.linalg.norm(self.centres - numpy.asarray(u)[:, None], axis=0)
        self.distances.append(distances)
        nearest = int(numpy.argmin(distances))
        active = local_decoder.GetActiveModel()
        if local_decoder.is_initialized and distances[nearest] >= (1.0 - self.margin) * distances[active]:
            return active
        return nearest


def CreateModelSelector(settings, rom_data_folder, rows):
    """The selector of a LocalRomDecoder from its settings: {"type": "process_info", "variable": "TIME", "thresholds": [...]} or
    {"type": "nearest_centre", "centres": "<file>.npy", "margin": 0.0}, with the centres saved in rom_data_folder as the bases (one row
    per node and nodal unknown, one column per model) and 'rows' their rows in the order of the effective DOF set."""
    if settings["type"] == "process_info":
        return ProcessInfoSelector(settings["variable"], settings["thresholds"])
    if settings["type"] == "nearest_centre":
        return NearestCentreSelector(numpy.load(pathlib.Path(rom_data_folder) / settings["centres"])[rows, :], settings.get("margin", 0.0))
    raise Exception("Unknown local ROM selector '{}'. Available options are 'process_info' and 'nearest_centre'.".format(settings["type"]))


RegisterDecoder("local_roms", LocalRomDecoder)

"""A trained surrogate as a Kratos response function.

The mirror image of ``cosim_surrogate_solver_wrapper``. That one lets a
model take a *solver*'s place in a co-simulation; this one lets a model take
a *response function*'s place, which is the role Kratos's optimization
tooling drives: value in, gradient out, both through
``ResponseFunctionInterface``.

    kratos/python_scripts/response_functions/response_function_interface.py

Because that interface lives in the core and is duck-typed, anything that
consumes a Kratos response function consumes this one - it is created by the
same ``CreateResponseFunction(response_id, response_settings, model)``
signature the applications' own factories use, so a driver that resolves
responses by module path needs no special case.

Two gradient modes, and the difference between them is exactly the
difference between the two dJ/dX paths this application already ships:

``"surrogate"``
    Autograd straight through the model's forward
    (``sensitivity_utils.ComputeSurrogateSensitivities``). Cheap - one
    backward pass, no assembly - and as accurate as the surrogate is. This is
    the mode that makes a surrogate worth having in an optimization loop.
``"exact"``
    The FEM adjoint (``sensitivity_utils.ComputeShapeSensitivityField``) on
    the state the surrogate just wrote. The surrogate replaces the *solve*,
    not the sensitivity analysis, so the gradient is discretely exact for the
    state it is given - and therefore only as trustworthy as that state. Use
    it to check the surrogate mode, or when the prediction is a warm start
    good enough to differentiate around.

torch and the model are loaded lazily, at the first evaluation.
"""

import numpy

import KratosMultiphysics as Kratos
from KratosMultiphysics.response_functions.response_function_interface import (
    ResponseFunctionInterface)

from KratosMultiphysics.PhysicsNeMoApplication.bridges import adjoint_bridge
from KratosMultiphysics.PhysicsNeMoApplication.bridges import torch_bridge
from KratosMultiphysics.PhysicsNeMoApplication.deployment import model_registry
from KratosMultiphysics.PhysicsNeMoApplication.physics import differentiable_residual
from KratosMultiphysics.PhysicsNeMoApplication.physics import sensitivity_utils
from KratosMultiphysics.PhysicsNeMoApplication.processes.inference.inference_process import (
    GatherInputFields, WriteOutputFields)
from KratosMultiphysics.PhysicsNeMoApplication.processes.inference.point_cloud_inference_process import (
    GatherPointCloudCoordinates, RunPointCloudForward, _MODEL_INTERFACES)

_GRADIENT_MODES = ("surrogate", "exact")
_RESPONSE_INTERFACES = ("flat",) + _MODEL_INTERFACES


def CreateResponseFunction(response_id, response_settings, model):
    """Factory with the signature Kratos's response-function factories use.

    Point ``"response_module"`` at this module (in
    ``adjoint_bridge.CreateResponseFunction``, or in any other driver that
    resolves a response by module path) to put a surrogate where an adjoint
    analysis would go.
    """
    return SurrogateResponseFunction(response_id, response_settings, model)


class SurrogateResponseFunction(ResponseFunctionInterface):
    """J and dJ/dX from a trained model.

    Settings:
        {
            "model_part_name"     : "",
            "model_settings"      : { },
            "input_fields"        : [ { "variable_name" : "", "data_location" : "node_historical" } ],
            "output_fields"       : [ { "variable_name" : "", "data_location" : "node_historical" } ],
            "objective"           : { },
            "gradient_mode"       : "surrogate",
            "model_interface"     : "flat",
            "normalize_coordinates" : false,
            "pass_geometry"       : true,
            "sensitivity_variable": "SHAPE_SENSITIVITY",
            "dof_fields"          : [ ],
            "fd_step"             : 1e-6,
            "design_sub_model_part_name" : ""
        }

    ``objective`` is an ``adjoint_bridge.MakeObjectiveWeights`` block - the
    same one the sensitivity process uses, deliberately, so "the objective"
    means one thing on both sides of a comparison. ``dof_fields`` is needed
    only by ``"exact"``.
    """

    def __init__(self, response_id, response_settings: Kratos.Parameters, model: Kratos.Model):
        self.identifier = response_id
        self.response_settings = response_settings

        defaults = Kratos.Parameters("""{
            "response_type"              : "surrogate",
            "model_part_name"            : "PLEASE_SPECIFY_MODEL_PART_NAME",
            "model_settings"             : {},
            "input_fields"               : [],
            "output_fields"              : [],
            "objective"                  : {},
            "gradient_mode"              : "surrogate",
            "model_interface"            : "flat",
            "normalize_coordinates"      : false,
            "pass_geometry"              : true,
            "sensitivity_variable"       : "SHAPE_SENSITIVITY",
            "dof_fields"                 : [],
            "fd_step"                    : 1e-6,
            "design_sub_model_part_name" : ""
        }""")
        settings = response_settings.Clone()
        settings.ValidateAndAssignDefaults(defaults)

        self.model_part = model[settings["model_part_name"].GetString()]
        self.model_settings = settings["model_settings"].Clone()
        self.objective_settings = settings["objective"].Clone()

        self.input_specs = self._ReadFieldSpecs(settings["input_fields"])
        self.output_specs = self._ReadFieldSpecs(settings["output_fields"])
        self.dof_field_specs = self._ReadFieldSpecs(settings["dof_fields"])

        self.gradient_mode = settings["gradient_mode"].GetString()
        if self.gradient_mode not in _GRADIENT_MODES:
            raise ValueError(
                f"Unsupported \"gradient_mode\" \"{self.gradient_mode}\". Use "
                f"{' or '.join(chr(34) + m + chr(34) for m in _GRADIENT_MODES)}.")

        self.model_interface = settings["model_interface"].GetString()
        if self.model_interface not in _RESPONSE_INTERFACES:
            raise ValueError(
                f"Unsupported \"model_interface\" \"{self.model_interface}\". Supported: "
                f"{', '.join(_RESPONSE_INTERFACES)}.")
        if self.gradient_mode == "surrogate" and self.model_interface == "flat":
            raise ValueError(
                "\"gradient_mode\" : \"surrogate\" differentiates the prediction with respect "
                "to the node COORDINATES, which the \"flat\" interface never sees - it feeds "
                "the model field values only. Use a point-cloud interface "
                f"({', '.join(_MODEL_INTERFACES)}) or \"gradient_mode\" : \"exact\".")

        self.normalize_coordinates = settings["normalize_coordinates"].GetBool()
        self.pass_geometry = settings["pass_geometry"].GetBool()
        self.sensitivity_variable = Kratos.KratosGlobals.GetVariable(
            settings["sensitivity_variable"].GetString())
        self.fd_step = settings["fd_step"].GetDouble()

        design_name = settings["design_sub_model_part_name"].GetString()
        self.design_node_ids = None
        if design_name:
            self.design_node_ids = [node.Id
                                    for node in self.model_part.GetSubModelPart(design_name).Nodes]

        if self.gradient_mode == "exact" and not self.dof_field_specs:
            raise ValueError(
                "\"gradient_mode\" : \"exact\" needs \"dof_fields\" - the solved unknowns the "
                "FEM adjoint system is posed on.")

        self._model = None
        self._device = None
        self._assembler = None
        self._dof_map = None
        self._value = None
        self._nodal_gradient = None

    @staticmethod
    def _ReadFieldSpecs(fields: Kratos.Parameters):
        return [(fields[i]["variable_name"].GetString(),
                 fields[i]["data_location"].GetString())
                for i in range(fields.size())]

    def _GetModel(self):
        if self._model is None:
            self._model, self._device = model_registry.LoadModelWithCardCheck(
                self.model_settings, self.input_specs, self.output_specs, type(self).__name__)
        return self._model

    def _GetNormalization(self):
        """The card's output de-normalization, or None.

        J is a function of the PHYSICAL prediction, so this is applied on
        the write path and inside the autograd objective alike - a
        de-normalization applied only where the field is written would
        leave dJ/dX wrong by exactly the training scale.
        """
        if not hasattr(self, "_normalization"):
            self._normalization = model_registry.LoadOutputNormalization(self.model_settings)
        return self._normalization

    def _GetInputNormalization(self):
        """The card's input normalization, or None (the symmetric half).
        Applied to the FIELD features only - the coordinates keep their own
        convention (``normalize_coordinates``), so the chain rule above is
        untouched by it."""
        if not hasattr(self, "_input_normalization"):
            self._input_normalization = model_registry.LoadInputNormalization(self.model_settings)
        return self._input_normalization

    def _ObjectiveWeights(self):
        return adjoint_bridge.MakeObjectiveWeights(self.objective_settings, self.model_part)

    def _ObjectiveVariable(self):
        return Kratos.KratosGlobals.GetVariable(
            self.objective_settings["variable_name"].GetString())

    # ------------------------------------------------------------------ value

    def CalculateValue(self):
        """Runs the surrogate and reduces its prediction to J.

        The prediction is written onto the model part first, so the state the
        gradient is then taken around is the surrogate's - and so any process
        or output attached afterwards sees it.
        """
        model = self._GetModel()
        torch = torch_bridge._TryImportTorch()

        inputs, n_entities = GatherInputFields(self.model_part, self.input_specs)
        features = model_registry.ApplyInputNormalization(
            torch.cat(inputs, dim=-1), self._GetInputNormalization())

        if self.model_interface == "flat":
            with torch.no_grad():
                prediction = model(features.to(self._device)).cpu()
        else:
            coordinates = torch.from_numpy(GatherPointCloudCoordinates(
                self.model_part, self.normalize_coordinates))
            prediction, _ = RunPointCloudForward(
                model, self._device, self.model_interface, features, coordinates,
                self.pass_geometry)

        WriteOutputFields(self.model_part, self.output_specs, prediction, n_entities,
                          normalization=self._GetNormalization())

        weights = self._ObjectiveWeights()
        variable = self._ObjectiveVariable()
        location = self._LocationOf(variable.Name())
        self._value = adjoint_bridge.EvaluateObjective(
            self.model_part, weights, variable, location)
        return self._value

    def _LocationOf(self, variable_name: str) -> str:
        for name, data_location in self.output_specs:
            if name == variable_name:
                return data_location
        raise ValueError(
            f"The objective reads \"{variable_name}\", which is not among the "
            f"\"output_fields\" {[spec[0] for spec in self.output_specs]} - the objective "
            "must be a function of what the surrogate predicts.")

    def GetValue(self):
        if self._value is None:
            raise RuntimeError(
                f"Response \"{self.identifier}\" has no value yet - call CalculateValue().")
        return self._value

    # --------------------------------------------------------------- gradient

    def CalculateGradient(self):
        if self.gradient_mode == "surrogate":
            self._nodal_gradient = self._SurrogateGradient()
        else:
            self._nodal_gradient = self._ExactGradient()
        adjoint_bridge.WriteSensitivityField(
            self.model_part, self._nodal_gradient, self.sensitivity_variable,
            "node_non_historical")

    def _SurrogateGradient(self):
        """dJ/dX by autograd through the surrogate's own forward.

        With ``normalize_coordinates`` the model is fed x_norm =
        (x - min)/extent, so autograd returns dJ/dx_norm and the physical
        gradient needs the 1/extent chain rule applied here. What that
        chain rule NEGLECTS is that the bounding box is itself recomputed
        from the design every call, so min and extent depend on the design
        too: the result is exact for a design that does not move the
        bounding box (an interior design surface, the usual case) and
        approximate for one that does. Turning normalization off removes
        the approximation entirely, which is why it is off by default here
        and on by default in the deployment processes.
        """
        model = self._GetModel()
        torch = torch_bridge._TryImportTorch()

        inputs, _ = GatherInputFields(self.model_part, self.input_specs)
        features = model_registry.ApplyInputNormalization(
            torch.cat(inputs, dim=-1), self._GetInputNormalization())
        raw_coordinates = GatherPointCloudCoordinates(self.model_part, normalize=False)
        extent = numpy.ones(3, dtype=numpy.float64)
        if self.normalize_coordinates:
            low = raw_coordinates.min(axis=0)
            extent = raw_coordinates.max(axis=0) - low
            extent[extent == 0.0] = 1.0
            coordinates = (raw_coordinates - low) / extent
        else:
            coordinates = raw_coordinates

        weights = self._ObjectiveWeights()
        variable_name = self._ObjectiveVariable().Name()
        offset, width = self._PredictionColumns(variable_name)
        weight_tensor = torch.as_tensor(weights, dtype=torch.float64)
        normalization = self._GetNormalization()

        def Objective(prediction):
            # de-normalized INSIDE the graph (torch-native), so the gradient
            # carries the training scale, and the value matches the field
            # CalculateValue wrote
            prediction = model_registry.ApplyOutputNormalization(prediction, normalization)
            block = prediction[..., offset:offset + width].reshape(weight_tensor.shape)
            return (block.to(torch.float64) * weight_tensor).sum()

        result = sensitivity_utils.ComputeSurrogateSensitivities(
            model, features, coordinates, Objective, model_interface=self.model_interface,
            device=self._device, pass_geometry=self.pass_geometry, wrt=("coordinates",))
        self._value = result["objective"]
        return result["coordinates"] / extent[None, :]

    def _PredictionColumns(self, variable_name: str):
        """Column offset and width of an output field in the prediction.

        The widths come from the tensor adaptors rather than from the
        variable type, which is how WriteOutputFields splits the same
        prediction - guessing them separately would let the two disagree
        silently on anything but a scalar or a 3-component array.
        """
        from KratosMultiphysics.PhysicsNeMoApplication.utilities.tensor_adaptor_dataset_utils import (
            GetTensorAdaptor)
        offset = 0
        for name, data_location in self.output_specs:
            variable = Kratos.KratosGlobals.GetVariable(name)
            adaptor = GetTensorAdaptor(self.model_part, data_location, variable)
            width = int(numpy.prod(adaptor.data.shape[1:], dtype=int))
            if name == variable_name:
                return offset, width
            offset += width
        raise ValueError(
            f"The objective reads \"{variable_name}\", which is not among the "
            f"\"output_fields\" {[spec[0] for spec in self.output_specs]}.")

    def _ExactGradient(self):
        """dJ/dX by the FEM adjoint, around the state the surrogate wrote."""
        if self._assembler is None:
            self._assembler = differentiable_residual.TangentAssembler(self.model_part)
            self._dof_map = differentiable_residual.DofFieldMap(
                self._assembler, self.dof_field_specs)

        weights = self._ObjectiveWeights()
        variable_name = self._ObjectiveVariable().Name()
        selector = numpy.zeros((self._dof_map.n_nodes, self._dof_map.total_width),
                               dtype=numpy.float64)
        columns = ([variable_name] if weights.shape[1] == 1
                   else [f"{variable_name}_{axis}" for axis in "XYZ"])
        for column_index, name in enumerate(columns):
            if name not in self._dof_map.column_of:
                raise ValueError(
                    f"The objective reads \"{name}\", which is not among the \"dof_fields\" "
                    f"{[spec[0] for spec in self.dof_field_specs]}.")
            selector[:, self._dof_map.column_of[name]] = weights[:, column_index]

        return sensitivity_utils.ComputeShapeSensitivityField(
            self._assembler, self._dof_map, self._dof_map.FieldsToDofVector(selector),
            fd_step=self.fd_step, design_node_ids=self.design_node_ids)

    def GetNodalGradient(self, variable):
        """{node_id: [dx, dy, dz]} - Kratos's own gradient contract."""
        if variable != self.sensitivity_variable:
            raise RuntimeError(
                f"GetNodalGradient: no gradient for \"{variable.Name()}\"; this response "
                f"reports \"{self.sensitivity_variable.Name()}\" (its "
                "\"sensitivity_variable\" setting).")
        if self._nodal_gradient is None:
            raise RuntimeError(
                f"Response \"{self.identifier}\" has no gradient yet - call CalculateGradient().")
        return {node.Id: [float(v) for v in row]
                for node, row in zip(self.model_part.Nodes, self._nodal_gradient)}

    def GetElementalGradient(self, variable):
        raise RuntimeError(
            f"GetElementalGradient: no gradient for \"{variable.Name()}\". A surrogate "
            "response differentiates with respect to the node coordinates only; elemental "
            "sensitivities (thickness, Young's modulus) need Kratos's own adjoint elements - "
            "reach them through adjoint_bridge.CreateResponseFunction.")

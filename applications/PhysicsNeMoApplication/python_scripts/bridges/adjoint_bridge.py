"""Bridge to Kratos's own adjoint sensitivity stack.

The counterpart of ``rom_bridge`` for gradients. RomApplication is consumed
through its *file format*; the adjoint stack is consumed through its
*Python contract*, which lives in the core and not in any one application:

    kratos/python_scripts/response_functions/response_function_interface.py

``ResponseFunctionInterface`` is what every Kratos response function
implements - ``CalculateValue``/``CalculateGradient`` then ``GetValue`` and
``GetNodalGradient``/``GetElementalGradient``. Three compiled applications
here provide implementations (StructuralMechanics, ConvectionDiffusion,
CompressiblePotentialFlow), and the gradients themselves are produced by the
core ``SensitivityBuilder`` writing ``SHAPE_SENSITIVITY`` and friends. So this
module needs the core plus whichever application owns the response - never a
compiled optimization application.

Row-ordering contract (the same one every gather in this application uses,
shared with ``graph_bridge.NodePositions``, ``differentiable_residual.DofFieldMap``
and ``torch_bridge``):

    row r of a nodal field     <->  ``model_part.Nodes`` iteration order
    row r of an elemental field <-> ``model_part.Elements`` iteration order

Kratos's response functions return ``{entity_id: value}`` dicts instead, and a
dict is not an array: it carries no order, and its iteration order is the
adjoint model part's, which need not be the primal's. ``EvaluateResponse``
converts by *id*, so the arrays line up with everything else in the
application by construction.

**Which model part the gradient lives on differs by application** and is the
first thing that surprises people:

- ``StructuralMechanicsApplication`` builds a *separate* ``Kratos.Model`` and
  re-reads the mdpa for its adjoint part. The response therefore has to run
  with the case files in the working directory, and its dict is keyed by that
  second part's node ids.
- ``ConvectionDiffusionApplication`` puts primal and adjoint parts in the
  same model, replacing the elements with ``AdjointDiffusionElement``.

Both are handled the same way here because the conversion keys off ids.

**Do not call ``analysis.Finalize()`` before assembling anything adjoint.**
The boundary-condition processes release the DOFs they fixed in their
``ExecuteFinalizeSolutionStep``, leaving an unconstrained (singular) tangent.
The failure is quiet: you get finite sensitivities that are wrong by about
six orders of magnitude.

Pure Kratos + numpy; torch is imported lazily, exactly as in ``rom_bridge``.
"""

import dataclasses
import importlib

import numpy

import KratosMultiphysics as Kratos
import KratosMultiphysics.kratos_utilities as kratos_utilities

from KratosMultiphysics.PhysicsNeMoApplication.utilities.tensor_adaptor_dataset_utils import (
    GetTensorAdaptor)

# response_application -> (required applications, factory module path)
_RESPONSE_FACTORIES = {
    "structural_mechanics": (
        ("StructuralMechanicsApplication",),
        "KratosMultiphysics.StructuralMechanicsApplication."
        "structural_response_function_factory"),
    "convection_diffusion": (
        ("ConvectionDiffusionApplication",),
        "KratosMultiphysics.ConvectionDiffusionApplication.response_functions."
        "convection_diffusion_response_function_factory"),
    "compressible_potential_flow": (
        ("CompressiblePotentialFlowApplication",),
        "KratosMultiphysics.CompressiblePotentialFlowApplication."
        "potential_flow_response_function_factory"),
}

_INTERFACE_METHODS = ("CalculateValue", "CalculateGradient", "GetValue")


@dataclasses.dataclass(frozen=True)
class SensitivityFields:
    """One response evaluation, in the application's row order.

    Attributes:
        value: The response value J.
        nodal: {variable_name: (n_nodes, width) float64}, rows in
            ``model_part.Nodes`` order.
        elemental: {variable_name: (n_elements, width) float64}, rows in
            ``model_part.Elements`` order.
        node_ids: (n_nodes,) int64 ids in that same row order.
        element_ids: (n_elements,) int64 ids in that same row order.
    """
    value: float
    nodal: dict
    elemental: dict
    node_ids: numpy.ndarray
    element_ids: numpy.ndarray

    def AsTorchTensors(self, device=None) -> dict:
        """Every field as a torch tensor, keyed "<VARIABLE>__<location>".

        The keys match the ``.npz`` export convention, so a gradient read
        here and a gradient read from an exported sample are named alike.
        """
        from KratosMultiphysics.PhysicsNeMoApplication.utilities import array_backend_utils
        tensors = {}
        for name, array in self.nodal.items():
            tensors[f"{name}__node_historical"] = array_backend_utils.AsTorchTensor(array, device)
        for name, array in self.elemental.items():
            tensors[f"{name}__element"] = array_backend_utils.AsTorchTensor(array, device)
        return tensors


def CreateResponseFunction(settings: Kratos.Parameters, model: Kratos.Model):
    """Builds a Kratos response function from Parameters.

    Args:
        settings: A Parameters block with

            - ``response_id``: the identifier Kratos labels the response with.
            - ``response_application``: one of ``structural_mechanics``,
              ``convection_diffusion``, ``compressible_potential_flow`` -
              dispatched to that application's ``CreateResponseFunction``.
            - ``response_module``: an explicit module path exposing
              ``CreateResponseFunction(response_id, response_settings, model)``,
              overriding ``response_application``. The escape hatch for an
              application this table does not know, and the same module-path
              convention ``cosim_surrogate_solver_wrapper`` uses.
            - ``response_settings``: passed through verbatim to that factory
              (``response_type``, ``primal_settings``, ``sensitivity_settings``,
              ...). Not validated here: its schema is the application's.

        model: The Kratos Model the response builds its parts in.

    Returns:
        The response function object (a ``ResponseFunctionInterface``).
    """
    defaults = Kratos.Parameters("""{
        "response_id"           : "physicsnemo_response",
        "response_application"  : "",
        "response_module"       : "",
        "response_settings"     : {}
    }""")
    settings = settings.Clone()
    # response_settings is the foreign application's schema - validate around it
    settings.ValidateAndAssignDefaults(defaults)

    module_path = settings["response_module"].GetString()
    application_key = settings["response_application"].GetString()
    if module_path and application_key:
        raise ValueError(
            "Set either \"response_application\" or \"response_module\", not both "
            f"(got \"{application_key}\" and \"{module_path}\").")
    if not module_path:
        if application_key not in _RESPONSE_FACTORIES:
            raise ValueError(
                f"Unsupported \"response_application\" \"{application_key}\". Known: "
                f"{', '.join(sorted(_RESPONSE_FACTORIES))}. For any other application, name "
                "its factory module in \"response_module\" instead.")
        required_applications, module_path = _RESPONSE_FACTORIES[application_key]
        if not kratos_utilities.CheckIfApplicationsAvailable(*required_applications):
            raise RuntimeError(
                f"\"response_application\" \"{application_key}\" needs "
                f"{', '.join(required_applications)}, which is not compiled in this Kratos "
                "build. Add it to the configure script's application list and rebuild.")

    try:
        factory_module = importlib.import_module(module_path)
    except ImportError as e:
        raise ImportError(
            f"Could not import the response-function factory module \"{module_path}\".") from e
    if not hasattr(factory_module, "CreateResponseFunction"):
        raise AttributeError(
            f"Module \"{module_path}\" has no \"CreateResponseFunction\" - a response factory "
            "module must expose CreateResponseFunction(response_id, response_settings, model).")

    response = factory_module.CreateResponseFunction(
        settings["response_id"].GetString(), settings["response_settings"], model)

    missing = [name for name in _INTERFACE_METHODS if not callable(getattr(response, name, None))]
    if missing:
        raise TypeError(
            f"\"{module_path}\" returned {type(response).__name__}, which does not implement "
            f"{', '.join(missing)} - it is not a ResponseFunctionInterface.")
    return response


def _EntityIds(container, count: int) -> numpy.ndarray:
    return numpy.fromiter((entity.Id for entity in container), dtype=numpy.int64, count=count)


def _RowsForIds(row_ids: numpy.ndarray, wanted_ids: numpy.ndarray, what: str) -> numpy.ndarray:
    """Row of each wanted id inside row_ids, via searchsorted.

    The obvious ``{id: row}`` dict is what this replaced everywhere else in
    the application; searchsorted measured 6-15x faster on the gather paths
    and the reason is the same here - the dict is built in the interpreter,
    one Python object per entity.
    """
    order = numpy.argsort(row_ids)
    sorted_ids = row_ids[order]
    positions = numpy.searchsorted(sorted_ids, wanted_ids)
    numpy.clip(positions, 0, max(len(sorted_ids) - 1, 0), out=positions)
    if len(sorted_ids) == 0 or not numpy.array_equal(sorted_ids[positions], wanted_ids):
        unknown = wanted_ids[sorted_ids[positions] != wanted_ids] if len(sorted_ids) else wanted_ids
        raise RuntimeError(
            f"The response's gradient names {what} id(s) {sorted(unknown.tolist())[:8]} that are "
            "not in the given model part - the response was built on a different mesh.")
    return order[positions]


def _GradientToArray(gradient: dict, row_ids: numpy.ndarray, what: str) -> numpy.ndarray:
    """A Kratos {entity_id: value} gradient dict as a row-ordered array."""
    n_rows = len(row_ids)
    if not gradient:
        return numpy.zeros((n_rows, 1), dtype=numpy.float64)

    wanted_ids = numpy.fromiter((int(key) for key in gradient), dtype=numpy.int64,
                                count=len(gradient))
    rows = _RowsForIds(row_ids, wanted_ids, what)

    values = [numpy.atleast_1d(numpy.asarray(value, dtype=numpy.float64)).ravel()
              for value in gradient.values()]
    width = values[0].size
    if any(value.size != width for value in values):
        raise RuntimeError(
            f"The response returned {what} gradient entries of differing widths "
            f"({sorted({value.size for value in values})}).")

    # Entities the response says nothing about keep a zero row: a restricted
    # sensitivity model part is the normal case, not an error.
    array = numpy.zeros((n_rows, width), dtype=numpy.float64)
    array[rows] = numpy.stack(values)
    return array


def _CallGradientGetter(response, getter_name: str, variable, what: str):
    getter = getattr(response, getter_name, None)
    if not callable(getter):
        raise AttributeError(
            f"{type(response).__name__} has no \"{getter_name}\", so it cannot report "
            f"{what} sensitivities.")
    try:
        return getter(variable)
    except RuntimeError as e:
        # Kratos's own implementations raise a plain RuntimeError naming only
        # the variable; say what the caller can actually do about it.
        raise RuntimeError(
            f"{type(response).__name__}.{getter_name} refuses \"{variable.Name()}\" "
            f"({e}). A Kratos response reports only the sensitivity variables its "
            "\"sensitivity_settings\" asked SensitivityBuilder to compute - nodal gradients "
            "are normally SHAPE_SENSITIVITY alone.") from e


def EvaluateResponse(response, model_part: Kratos.ModelPart,
                     nodal_variables=(Kratos.SHAPE_SENSITIVITY,),
                     elemental_variables=(),
                     run_lifecycle: bool = True,
                     initialize: bool = True) -> SensitivityFields:
    """Runs a Kratos response function and returns its gradients as arrays.

    Args:
        response: A ``ResponseFunctionInterface`` (see CreateResponseFunction).
        model_part: The part whose entity order defines the rows. This is the
            *primal*/design part; the response may internally own a different
            adjoint part, which is why the conversion goes through ids.
        nodal_variables: Nodal sensitivity variables to read
            (``SHAPE_SENSITIVITY`` is what SensitivityBuilder writes and what
            Kratos's optimization tooling reads).
        elemental_variables: Elemental sensitivity variables to read, e.g.
            ``StructuralMechanicsApplication.YOUNG_MODULUS_SENSITIVITY``.
        run_lifecycle: Drive ``InitializeSolutionStep``/``CalculateValue``/
            ``CalculateGradient``. Pass False to read a response that has
            already been evaluated.
        initialize: Call ``Initialize`` first. Pass False when it has been
            called already - a second call re-initializes both analyses.

    Returns:
        A SensitivityFields.
    """
    if run_lifecycle:
        if initialize:
            response.Initialize()
        response.InitializeSolutionStep()
        response.CalculateValue()
        response.CalculateGradient()

    node_ids = _EntityIds(model_part.Nodes, model_part.NumberOfNodes())
    element_ids = _EntityIds(model_part.Elements, model_part.NumberOfElements())

    nodal = {}
    for variable in nodal_variables:
        gradient = _CallGradientGetter(response, "GetNodalGradient", variable, "nodal")
        nodal[variable.Name()] = _GradientToArray(gradient, node_ids, "node")

    elemental = {}
    for variable in elemental_variables:
        gradient = _CallGradientGetter(response, "GetElementalGradient", variable, "elemental")
        elemental[variable.Name()] = _GradientToArray(gradient, element_ids, "element")

    return SensitivityFields(value=float(response.GetValue()), nodal=nodal, elemental=elemental,
                             node_ids=node_ids, element_ids=element_ids)


def ReadSensitivityField(model_part: Kratos.ModelPart, variable,
                         data_location: str = "node_historical") -> numpy.ndarray:
    """Reads a sensitivity field a SensitivityBuilder already wrote.

    The response-function-free path: when the adjoint solver ran inside the
    analysis (rather than behind a response object), the sensitivities are
    simply Kratos variables. Goes through the shared tensor-adaptor factory,
    so the rows and the layout match every other gather in the application -
    and the key ``"<VARIABLE>__<data_location>"`` matches the ``.npz`` export.

    Returns:
        (n_entities, width) float64, with a trailing axis even for scalars.
    """
    tensor_adaptor = GetTensorAdaptor(model_part, data_location, variable)
    data = numpy.array(tensor_adaptor.data, dtype=numpy.float64)
    return data.reshape(data.shape[0], -1)


def WriteSensitivityField(model_part: Kratos.ModelPart, field, variable,
                          data_location: str = "node_non_historical") -> None:
    """Writes a sensitivity field into a Kratos variable.

    The direction that lets a surrogate's gradient be read by tooling that
    knows nothing about this application: ``SHAPE_SENSITIVITY`` is where
    Kratos's optimization drivers look.

    ``node_non_historical`` is the default because it needs no
    pre-allocation; the historical database only accepts variables added to
    the solution-step list before the mesh was read, which a primal analysis
    has no reason to have done for a sensitivity variable.
    """
    tensor_adaptor = GetTensorAdaptor(model_part, data_location, variable)
    expected = tensor_adaptor.data.shape
    field = numpy.ascontiguousarray(numpy.asarray(field, dtype=numpy.float64).reshape(expected))
    tensor_adaptor.data[:] = field
    tensor_adaptor.StoreData()


def _VariableWidth(variable) -> int:
    """Component count of a nodal variable: 1 for a scalar, 3 for an array."""
    if isinstance(variable, Kratos.DoubleVariable):
        return 1
    if isinstance(variable, Kratos.Array1DVariable3):
        return 3
    raise ValueError(
        f"Objective variable \"{variable.Name()}\" is neither a double nor a 3-component "
        "array variable; objectives are defined over nodal DOF fields.")


def MakeObjectiveWeights(settings: Kratos.Parameters, model_part: Kratos.ModelPart):
    """Builds the (n_nodes, width) weights w of a linear objective J = w . u.

    Both consumers of an objective in this application need the same thing
    and must agree on it exactly: the sensitivity process contracts the
    weights with the FEM state to get dJ/du, and ``SurrogateResponseFunction``
    contracts them with a surrogate's prediction to get J. Sharing one
    builder is what keeps "the objective" from meaning two different things
    on the two sides of a cross-validation.

    Rows follow ``model_part.Nodes`` order; columns are the variable's
    components (X, Y, Z for an array variable).

    Args:
        settings: Parameters with

            - ``type``: ``traced_node`` - one node's field, projected on
              ``direction``; or ``weighted_sum`` - the field summed over a
              node set, optionally scaled by a nodal weight variable.
            - ``variable_name``: the nodal field the objective reads.
            - ``node_id``: (traced_node) the traced node.
            - ``direction``: (traced_node) the projection direction; defaults
              to the Z axis for an array variable, ignored for a scalar.
            - ``model_part_name``: (weighted_sum) restrict the sum to this
              sub-model-part; empty means the whole part.
            - ``weight_variable_name``: (weighted_sum) a nodal double
              variable scaling each node's contribution - a nodal volume or
              area turns the sum into an integral.

        model_part: The part whose node order defines the rows.

    Returns:
        (n_nodes, width) float64 weights.
    """
    defaults = Kratos.Parameters("""{
        "type"                 : "traced_node",
        "variable_name"        : "PLEASE_SPECIFY_VARIABLE_NAME",
        "node_id"              : 0,
        "direction"            : [0.0, 0.0, 1.0],
        "model_part_name"      : "",
        "weight_variable_name" : ""
    }""")
    settings = settings.Clone()
    settings.ValidateAndAssignDefaults(defaults)

    variable = Kratos.KratosGlobals.GetVariable(settings["variable_name"].GetString())
    width = _VariableWidth(variable)
    node_ids = _EntityIds(model_part.Nodes, model_part.NumberOfNodes())
    weights = numpy.zeros((len(node_ids), width), dtype=numpy.float64)

    objective_type = settings["type"].GetString()
    if objective_type == "traced_node":
        node_id = settings["node_id"].GetInt()
        rows = _RowsForIds(node_ids, numpy.array([node_id], dtype=numpy.int64), "node")
        if width == 1:
            weights[rows[0], 0] = 1.0
        else:
            direction = numpy.asarray(settings["direction"].GetVector(), dtype=numpy.float64)
            if direction.size != 3:
                raise ValueError(
                    f"\"direction\" must have 3 entries, got {direction.size}.")
            weights[rows[0], :] = direction
    elif objective_type == "weighted_sum":
        sub_name = settings["model_part_name"].GetString()
        summed_part = model_part.GetSubModelPart(sub_name) if sub_name else model_part
        summed_ids = _EntityIds(summed_part.Nodes, summed_part.NumberOfNodes())
        rows = _RowsForIds(node_ids, summed_ids, "node")

        weight_name = settings["weight_variable_name"].GetString()
        if weight_name:
            weight_variable = Kratos.KratosGlobals.GetVariable(weight_name)
            if _VariableWidth(weight_variable) != 1:
                raise ValueError(
                    f"\"weight_variable_name\" \"{weight_name}\" must be a scalar variable.")
            per_node = numpy.fromiter(
                (node.GetValue(weight_variable) for node in summed_part.Nodes),
                dtype=numpy.float64, count=len(summed_ids))
        else:
            per_node = numpy.ones(len(summed_ids), dtype=numpy.float64)

        if width == 1:
            weights[rows, 0] = per_node
        else:
            direction = numpy.asarray(settings["direction"].GetVector(), dtype=numpy.float64)
            weights[rows, :] = per_node[:, None] * direction[None, :]
    else:
        raise ValueError(
            f"Unsupported objective \"type\" \"{objective_type}\". Use \"traced_node\" or "
            "\"weighted_sum\".")
    return weights


def EvaluateObjective(model_part: Kratos.ModelPart, weights, variable,
                      data_location: str = "node_historical") -> float:
    """J = sum(weights * field) for the weights MakeObjectiveWeights built."""
    field = ReadSensitivityField(model_part, variable, data_location)
    weights = numpy.asarray(weights, dtype=numpy.float64)
    if field.shape != weights.shape:
        raise ValueError(
            f"Objective weights have shape {list(weights.shape)} but the field is "
            f"{list(field.shape)}.")
    return float((weights * field).sum())

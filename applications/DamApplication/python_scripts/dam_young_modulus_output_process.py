import KratosMultiphysics as KM
import KratosMultiphysics.DamApplication as KratosDam


def Factory(settings, model):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a JSON string")
    return DamYoungModulusOutputProcess(model, settings["Parameters"])


class DamYoungModulusOutputProcess(KM.Process):
    """Prepare nonhistorical NODAL_YOUNG_MODULUS for post-processing.

    Young's modulus is usually stored as a constant element material property,
    while the available GiD/VTK nodal result is the nonhistorical
    NODAL_YOUNG_MODULUS field. This process projects the element property onto
    the element nodes so that the field can be visualized.

    Projection: every element whose Properties define a constant YOUNG_MODULUS
    contributes its value to each of its nodes. The nodal value is the
    arithmetic mean over the contributing elements. Consequently, a node shared
    by elements with different materials receives an averaged value. The
    resulting field is a visualization projection, not an exact elemental
    material property.

    The projection is intentionally limited to constant properties:

    - Accessor-backed properties are skipped. A genuine spatially varying
      Young's modulus must be written through the historical YOUNG_MODULUS
      field, which this process does not modify.
    - Elements without YOUNG_MODULUS (for instance joint/interface elements)
      do not contribute. Nodes without any contributing element are left
      untouched, so standard writers may still display their default zero.

    This process does not modify Properties, constitutive laws, element
    formulations, accessors, or the mechanical solution.
    """

    def __init__(self, model, settings):
        super().__init__()

        settings.ValidateAndAssignDefaults(KM.Parameters('''{
            "model_part_name" : "MainModelPart"
        }'''))

        self.model = model
        self.model_part_name = settings["model_part_name"].GetString()
        self.model_part = None

    def ExecuteBeforeSolutionLoop(self):
        self._PrepareNodalValues()

    def ExecuteBeforeOutputStep(self):
        self._PrepareNodalValues()

    def _GetModelPart(self):
        # Resolve lazily: the process may be constructed before the solver has
        # imported the model part.
        if self.model_part is None:
            self.model_part = self.model[self.model_part_name]
        return self.model_part

    def _PrepareNodalValues(self):
        model_part = self._GetModelPart()

        value_sum = {}
        contribution_count = {}

        for element in model_part.Elements:
            properties = element.Properties
            if properties.HasAccessor(KM.YOUNG_MODULUS):
                continue
            if not properties.Has(KM.YOUNG_MODULUS):
                continue

            young_modulus = properties[KM.YOUNG_MODULUS]
            for node in element.GetGeometry():
                node_id = node.Id
                value_sum[node_id] = value_sum.get(node_id, 0.0) + young_modulus
                contribution_count[node_id] = contribution_count.get(node_id, 0) + 1

        for node in model_part.Nodes:
            count = contribution_count.get(node.Id, 0)
            if count:
                node.SetValue(
                    KratosDam.NODAL_YOUNG_MODULUS,
                    value_sum[node.Id] / count)

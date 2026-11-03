"""Focused coordination tests for multi-part initial-state imports."""
from pathlib import Path
import KratosMultiphysics as KM
import KratosMultiphysics.DamApplication as Dam
import KratosMultiphysics.KratosUnittest as KratosUnittest

# Make the source-tree Dam scripts visible when running from the build tree.
Dam.__path__.insert(0, str(Path(__file__).resolve().parents[1] / "python_scripts"))
from KratosMultiphysics.DamApplication import dam_initial_state_file


class _Element:
    def __init__(self, element_id):
        self.Id = element_id


class _Part:
    def __init__(self, name, element_ids):
        self.Name = name
        self.Elements = [_Element(element_id) for element_id in element_ids]


class _Process:
    """Minimal stand-in for LoadDamInitialStateProcess."""

    def __init__(self, name, file_name, element_ids, temperature=None,
                 reference=None, stresses=True, transfer_temperature=False,
                 transfer_reference=False, policy="preserve"):
        self.input_file_name = file_name
        self.model_part = _Part(name, element_ids)
        self.transfer_initial_stresses = stresses
        self.transfer_initial_temperature = transfer_temperature
        self.transfer_reference_temperature = transfer_reference
        self.reference_temperature_policy = policy
        fields = []
        if stresses:
            fields.append("initial_stresses")
        if transfer_temperature or policy == "initial_temperature":
            fields.append("temperature")
        if transfer_reference:
            fields.append("reference_temperature")
        self._data = {
            "fields": fields,
            "temperature": temperature or {},
            "reference_temperature": reference or {}}
        if stresses:
            self._data["initial_stresses"] = {
                element_id: [[0.0, 0.0, 0.0]] for element_id in element_ids}

    def PrevalidateInitialState(self):
        return self._data


class TestDamInitialStateMultipart(KratosUnittest.TestCase):

    def test_disjoint_stress_elements_are_allowed(self):
        concrete = _Process("Parts_Concrete", "concrete.json", [1, 2])
        soil = _Process("Parts_Soil", "soil.json", [3, 4])
        dam_initial_state_file.ValidateCombinedImports([concrete, soil])

    def test_overlapping_stress_elements_are_rejected(self):
        concrete = _Process("Parts_Concrete", "concrete.json", [1, 2])
        overlap = _Process("Parts_Soil", "soil.json", [2, 3])
        with self.assertRaisesRegex(RuntimeError, "element 2"):
            dam_initial_state_file.ValidateCombinedImports([concrete, overlap])

    def test_compatible_shared_temperature_is_allowed(self):
        concrete = _Process(
            "Parts_Concrete", "concrete.json", [1],
            temperature={2: 20.0}, transfer_temperature=True)
        soil = _Process(
            "Parts_Soil", "soil.json", [2],
            temperature={2: 20.0}, transfer_temperature=True)
        dam_initial_state_file.ValidateCombinedImports([concrete, soil])

    def test_conflicting_shared_temperature_is_rejected(self):
        concrete = _Process(
            "Parts_Concrete", "concrete.json", [1],
            temperature={2: 20.0}, transfer_temperature=True)
        soil = _Process(
            "Parts_Soil", "soil.json", [2],
            temperature={2: 25.0}, transfer_temperature=True)
        with self.assertRaisesRegex(RuntimeError, "shared node 2 receives TEMPERATURE"):
            dam_initial_state_file.ValidateCombinedImports([concrete, soil])

    def test_initial_temperature_policy_conflicts_with_reference_field(self):
        concrete = _Process(
            "Parts_Concrete", "concrete.json", [1],
            temperature={2: 20.0}, transfer_temperature=True,
            policy="initial_temperature")
        soil = _Process(
            "Parts_Soil", "soil.json", [2],
            reference={2: 15.0}, stresses=False, transfer_reference=True)
        with self.assertRaisesRegex(
                RuntimeError, "NODAL_REFERENCE_TEMPERATURE"):
            dam_initial_state_file.ValidateCombinedImports([concrete, soil])

    def test_different_fields_on_shared_node_do_not_conflict(self):
        concrete = _Process(
            "Parts_Concrete", "concrete.json", [1],
            temperature={2: 20.0}, transfer_temperature=True)
        soil = _Process(
            "Parts_Soil", "soil.json", [2],
            reference={2: 15.0}, stresses=False, transfer_reference=True)
        dam_initial_state_file.ValidateCombinedImports([concrete, soil])


if __name__ == "__main__":
    KratosUnittest.main()

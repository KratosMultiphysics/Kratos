import json


class MaterialEditor:
    def __init__(self, json_path):
        self.json_path = json_path
        self.data = self.load_json()

    def _convert_type(self, value_string):
        try:
            return float(value_string)
        except ValueError:
            return value_string

    def load_json(self):
        try:
            with open(self.json_path, 'r') as f:
                return json.load(f)
        except (FileNotFoundError, json.JSONDecodeError) as e:
            raise RuntimeError(f"Failed to load JSON file: {e}")

    def update_material_properties(self, entries: dict):
        variables = self.data["properties"][0]["Material"]["Variables"]
        for key, entry in entries.items():
            value = entry
            if isinstance(entry, list):
                value_str = [str(x).strip() for x in entry]
                value = [self._convert_type(x) for x in value_str]
            variables[key] = value
        self.save()

    def set_constitutive_law(self, law_name: str):
        self.data["properties"][0]["Material"]["constitutive_law"]["name"] = law_name
        self.save()

    def save(self):
        with open(self.json_path, 'w') as f:
            json.dump(self.data, f, indent=4)

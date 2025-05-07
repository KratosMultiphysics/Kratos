import json


class MaterialEditor:
    def __init__(self, json_path):
        self.json_path = json_path
        self.data = self._load_json()

    def _load_json(self):
        with open(self.json_path, 'r') as f:
            data = json.load(f)
        return data

    def _update_material_and_save(self, entries: dict):
        variables = self.data["properties"][0]["Material"]["Variables"]
        for key, entry in entries.items():
            value = entry
            if isinstance(entry, list):
                value_str = [str(x).strip() for x in entry]
                value = [self._convert_type(x) for x in value_str]
            variables[key] = value

        with open(self.json_path, 'w') as f:
            json.dump(self.data, f, indent=4)

    def _convert_type(self, value):
        try:
            if '.' in value or 'e' in value.lower():
                return float(value)
            else:
                return int(value)
        except ValueError:
            return value

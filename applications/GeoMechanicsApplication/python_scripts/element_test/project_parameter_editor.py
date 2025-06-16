import re
import json
from tkinter import messagebox


class ProjectParameterEditor:
    def __init__(self, json_path):
        self.json_path = json_path
        with open(self.json_path, 'r') as f:
            self.raw_text = f.read()

    def _write_back(self):
        with open(self.json_path, 'w') as f:
            f.write(self.raw_text)

    def update_nested_value(self, module_name, key, new_list):
        try:
            data = json.loads(self.raw_text)

            loads_list = data.get("processes", {}).get("loads_process_list", [])
            for process in loads_list:
                if (process.get("python_module") == module_name
                    and key in process.get("Parameters", {})
                ):
                    process["Parameters"][key] = new_list
                    break
            else:
                messagebox.showwarning("Warning", f"Could not find '{key}' under '{module_name}'.")
                return

            self.raw_text = json.dumps(data, indent=4)
            self._write_back()
        except Exception as e:
            messagebox.showerror("Error", f"Failed to update stress vector: {e}")

    def update_property(self, property_name, new_value):
        pattern = rf'("{property_name}"\s*:\s*)([0-9eE+\.\-]+)'
        replacement = rf'\g<1>{new_value}'
        self.raw_text, count = re.subn(pattern, replacement, self.raw_text)
        if count == 0:
            messagebox.showwarning("Warning", f"Could not find '{property_name}' to update.")
        elif count > 1:
            messagebox.showwarning("Warning", f"Multiple occurrences of '{property_name}' found. Updated all {count}.")
        self._write_back()

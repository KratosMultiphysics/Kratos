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

    def _update_property(self, property_name, new_value):
        pattern = rf'("{property_name}"\s*:\s*)([0-9eE+\.\-]+)'
        replacement = rf'\g<1>{new_value}'
        self.raw_text, count = re.subn(pattern, replacement, self.raw_text)
        if count == 0:
            messagebox.showwarning("Warning", f"Could not find '{property_name}' to update.")
        else:
            self._write_back()

    def _update_nested_value(self, module_name, key, new_list):
        try:
            data = json.loads(self.raw_text)

            loads_list = data.get("processes", {}).get("loads_process_list", [])
            for process in loads_list:
                if (
                        process.get("python_module") == module_name
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

    def update_time_step_properties(self, number_of_step, end_time):
        new_time_step = end_time / number_of_step
        ProjectParameterEditor._update_property(self, 'time_step', new_time_step)

    def update_end_time_properties(self, end_time):
        ProjectParameterEditor._update_property(self, 'end_time', end_time)

    def update_initial_stress_vector(self,  initial_effective_cell_pressure):

        new_stress_values = [-initial_effective_cell_pressure] * 3 + [0.0]
        self._update_nested_value("apply_initial_uniform_stress_field", "value", new_stress_values)


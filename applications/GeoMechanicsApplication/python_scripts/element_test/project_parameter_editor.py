import re
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

    def update_time_step_properties(self, number_of_step, end_time):
        new_time_step = end_time / number_of_step
        ProjectParameterEditor._update_property(self, 'time_step', new_time_step)

    def update_end_time_properties(self, end_time):
        ProjectParameterEditor._update_property(self, 'end_time', end_time)


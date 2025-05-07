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

    def update_time_step_properties(self, number_of_step, end_time):
        new_time_step = end_time / number_of_step
        pattern = r'("time_step"\s*:\s*)([0-9eE+\.\-]+)'
        replacement = r'\g<1>' + str(new_time_step)
        self.raw_text, count = re.subn(pattern, replacement, self.raw_text)
        if count == 0:
            messagebox.showwarning("Warning", "Could not find 'time_step' to update.")
        else:
            self._write_back()

    def update_end_time_properties(self, end_time):
        pattern = r'("end_time"\s*:\s*)([0-9eE+\.\-]+)'
        replacement = r'\g<1>' + str(end_time)
        self.raw_text, count = re.subn(pattern, replacement, self.raw_text)
        if count == 0:
            messagebox.showwarning("Warning", "Could not find 'end_time' to update.")
        else:
            self._write_back()

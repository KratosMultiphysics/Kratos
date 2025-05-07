import re
from tkinter import messagebox

class MdpaEditor:
    def __init__(self, mdpa_path):
        self.mdpa_path = mdpa_path
        with open(self.mdpa_path, 'r') as f:
            self.raw_text = f.read()

    def update_maximum_strain(self, maximum_strain):
        pattern = r'(\s*)\$maximum_strain(\s*)'
        prescribed_displacement = -maximum_strain / 100

        def replacer(match):
            leading_ws = match.group(1)
            trailing_ws = match.group(2)
            return f"{leading_ws}{prescribed_displacement}{trailing_ws}"

        new_text, count = re.subn(pattern, replacer, self.raw_text)
        if count == 0:
            messagebox.showwarning("Warning", "Could not find '$maximum_strain' to update.")
        else:
            self.raw_text = new_text
        with open(self.mdpa_path, 'w') as f:
            f.write(self.raw_text)
        # messagebox.showinfo("Success", f"'$maximum_strain' replaced with {prescribed_displacement}")

    def update_initial_effective_cell_pressure(self, initial_effective_cell_pressure):
        pattern = r'(\s*)\$initial_effective_cell_pressure(\s*)'

        def replacer(match):
            leading_ws = match.group(1)
            trailing_ws = match.group(2)
            return f"{leading_ws}{initial_effective_cell_pressure}{trailing_ws}"

        new_text, count = re.subn(pattern, replacer, self.raw_text)
        if count == 0:
            messagebox.showwarning("Warning", "Could not find '$initial_effective_cell_pressure' to update.")
        else:
            self.raw_text = new_text
        with open(self.mdpa_path, 'w') as f:
            f.write(self.raw_text)
        #     messagebox.showinfo("Success", f"'$initial_effective_cell_pressure' replaced with {initial_effective_cell_pressure}")

    def update_first_timestep(self, num_steps):
        first_timestep = 1.0 / num_steps
        pattern = r'(\s*)\$first_timestep(\s*)'

        def replacer(match):
            leading_ws = match.group(1)
            trailing_ws = match.group(2)
            return f"{leading_ws}{first_timestep}{trailing_ws}"

        new_text, count = re.subn(pattern, replacer, self.raw_text)
        if count == 0:
            messagebox.showwarning("Warning", "Could not find '$first_timestep' to update.")
        else:
            with open(self.mdpa_path, 'w') as f:
                f.write(new_text)
        #     messagebox.showinfo("Success", f"'$first_timestep' replaced with {first_timestep}")

    def update_end_time(self, end_time):
        pattern = r'\$end_time\b'
        replacement = str(end_time)

        new_text, count = re.subn(pattern, replacement, self.raw_text)

        if count == 0:
            messagebox.showwarning("Warning", "Could not find '$end_time' to update.")
        else:
            self.raw_text = new_text
            with open(self.mdpa_path, 'w') as f:
                f.write(self.raw_text)
            # messagebox.showinfo("Success", f"'$end_time' replaced with {end_time}")

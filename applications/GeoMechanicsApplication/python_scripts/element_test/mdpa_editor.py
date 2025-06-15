import re
from ui_logger import log_message


class MdpaEditor:
    def __init__(self, mdpa_path):
        self.mdpa_path = mdpa_path
        try:
            with open(self.mdpa_path, 'r') as f:
                self.raw_text = f.read()
        except FileNotFoundError:
            raise RuntimeError(f"File not found: {self.mdpa_path}")

    def save(self):
        with open(self.mdpa_path, 'w') as f:
            f.write(self.raw_text)

    def _replacer_factory(self, variable_value):
        def replacer(match):
            return f"{variable_value:.4f}"
        return replacer

    def update_maximum_strain(self, maximum_strain):
        pattern = r'\$maximum_strain\b'
        prescribed_displacement = -maximum_strain / 100

        replacer = self._replacer_factory(prescribed_displacement)
        new_text, count = re.subn(pattern, replacer, self.raw_text)
        if count == 0:
            log_message("Could not update maximum strain.", "Warning")
        else:
            self.raw_text = new_text
            MdpaEditor.save(self)

    def update_initial_effective_cell_pressure(self, initial_effective_cell_pressure):
        pattern = r'\$initial_effective_cell_pressure\b'

        replacer = self._replacer_factory(initial_effective_cell_pressure)
        new_text, count = re.subn(pattern, replacer, self.raw_text)
        if count == 0:
            log_message("Could not update initial effective cell pressure.", "Warning")
        else:
            self.raw_text = new_text
            MdpaEditor.save(self)

    def update_first_timestep(self, num_steps, end_time):
        first_timestep = end_time / num_steps
        pattern = r'\$first_timestep\b'

        replacer = self._replacer_factory(first_timestep)
        new_text, count = re.subn(pattern, replacer, self.raw_text)
        if count == 0:
            log_message("Could not apply the first time step.", "Warning")
        else:
            self.raw_text = new_text
            MdpaEditor.save(self)

    def update_end_time(self, end_time):
        pattern = r'\$end_time\b'
        replacement = str(end_time)

        new_text, count = re.subn(pattern, replacement, self.raw_text)

        if count == 0:
            log_message("Could not update the end time.", "Warning")
        else:
            self.raw_text = new_text
            MdpaEditor.save(self)

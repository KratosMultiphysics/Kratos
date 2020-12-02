class FiniteDifferenceDragShapeSensitivityAnalysis:
    def __init__(self, domain_size, check_node_ids, model_part_file_name, primal_analysis):
        self.domain_size = domain_size
        self.model_part_file_name = model_part_file_name
        self.check_node_ids = check_node_ids
        self.primal_analysis = primal_analysis
        self.check_node_id_coordinates = []

    def ComputeSensitivity(self, drag_direction, drag_file_name, delta):
        self._ReadNodalCoordinates()

        finite_difference_drag_sensitivities = []
        for index, check_node_id in enumerate(self.check_node_ids):
            coordinates = self.check_node_id_coordinates[index]
            fd_sensitivity = [0.0, 0.0, 0.0]
            for k in range(self.domain_size):
                coordinates[k] += delta
                self._WriteNodalCoordinates(check_node_id, coordinates)

                self.primal_analysis().Run()
                perturbed_drag = FiniteDifferenceDragShapeSensitivityAnalysis._GetTimeAveragedDrag(drag_direction, drag_file_name)
                fd_sensitivity[k] = perturbed_drag / delta

                coordinates[k] -= delta
                self._WriteNodalCoordinates(check_node_id, coordinates)

            finite_difference_drag_sensitivities.append(fd_sensitivity)

        # do the unperturbed state last, so the same results can be used for adjoint method
        self.primal_analysis().Run()
        ref_drag = FiniteDifferenceDragShapeSensitivityAnalysis._GetTimeAveragedDrag(drag_direction, drag_file_name)

        for index, _ in enumerate(self.check_node_ids):
            for k in range(self.domain_size):
                finite_difference_drag_sensitivities[index][k] -= ref_drag / delta

        return finite_difference_drag_sensitivities


    def _ReadNodalCoordinates(self):
        with open(self.model_part_file_name + '.mdpa', 'r') as model_part_file:
            lines = model_part_file.readlines()

        lines = lines[lines.index('Begin Nodes\n'):lines.index('End Nodes\n')]
        for check_node_id in self.check_node_ids:
            line = lines[check_node_id] # assumes consecutive node numbering starting with 1
            components = line.split()
            if int(components[0]) != check_node_id:
                raise RuntimeError('Error parsing file ' + self.model_part_file_name)
            self.check_node_id_coordinates.append([float(components[i]) for i in range(1,4)])

    def _WriteNodalCoordinates(self, node_id, new_coordinates):
        with open(self.model_part_file_name + '.mdpa', 'r') as model_part_file:
            lines = model_part_file.readlines()
        node_lines = lines[lines.index('Begin Nodes\n'):lines.index('End Nodes\n')]
        old_line = node_lines[node_id] # assumes consecutive node numbering starting with 1
        components = old_line.split()
        if int(components[0]) != node_id:
            raise RuntimeError('Error parsing file ' + self.model_part_file_name)
        new_line = '{:5d}'.format(node_id) + ' ' \
             + '{:19.10f}'.format(new_coordinates[0]) + ' ' \
             + '{:19.10f}'.format(new_coordinates[1]) + ' ' \
             + '{:19.10f}'.format(new_coordinates[2]) + '\n'
        lines[lines.index(old_line)] = new_line
        with open(self.model_part_file_name + '.mdpa', 'w') as model_part_file:
            model_part_file.writelines(lines)

    @staticmethod
    def _GetTimeAveragedDrag(direction, drag_file_name):
        time_steps, reactions = FiniteDifferenceDragShapeSensitivityAnalysis._ReadDrag(drag_file_name)
        total_drag = 0.0
        for reaction in reactions:
            total_drag += reaction[0]*direction[0]+reaction[1]*direction[1]+reaction[2]*direction[2]
        if len(time_steps) > 1:
            delta_time = time_steps[1] - time_steps[0]
            total_drag *= delta_time
        return total_drag

    @staticmethod
    def _ReadDrag(file_name):
        with open(file_name, "r") as file_input:
            lines = file_input.readlines()
        time_steps = []
        reaction = []
        for line in lines:
            line = line.strip()
            if len(line) == 0 or line[0] == '#':
                continue
            time_step_data = [float(v) for v in line.split()]
            time, fx, fy, fz = time_step_data
            time_steps.append(time)
            reaction.append([fx, fy, fz])
        return time_steps, reaction


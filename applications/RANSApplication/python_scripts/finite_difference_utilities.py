import os, math
import KratosMultiphysics as Kratos

from KratosMultiphysics.kratos_utilities import CheckIfApplicationsAvailable
if CheckIfApplicationsAvailable("FluidDynamicsApplication"):
    from KratosMultiphysics.FluidDynamicsApplication.fluid_dynamics_analysis import FluidDynamicsAnalysis
else:
    msg = "RANSApplication requires FluidDynamicsApplication which is not found."
    msg += " Please install/compile it and try again."
    raise Exception(msg)


def _is_different(line1, line2, tol=1e-10):
    for e1, e2 in zip(line1.split(), line2.split()):
        if not math.AlmostEqual(float(e1), float(e2), tol):
            return True
    return False


def _get_position(node_id, mdpa_lines):
    index_begin_nodes = mdpa_lines.index('Begin Nodes\n')
    index_end_nodes = mdpa_lines.index('End Nodes\n')
    for i in range(index_begin_nodes + 1, index_end_nodes):
        if mdpa_lines[i].split()[0] == str(node_id):
            return i
    raise RuntimeError("Node not found: " + str(node_id))


def _read_drag_file(filename):
    with open(filename, "r") as file_input:
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


class MDPAUtilities():
    def __init__(self, file_name):
        self.file_name = file_name
        with open(self.file_name, "r") as file_input:
            mdpa_lines = file_input.readlines()
        self.mdpa_lines = mdpa_lines

    def GetNodalCoordinates(self, node_id):
        _, x, y, z = self.mdpa_lines[_get_position(node_id,
                                                   self.mdpa_lines)].split()
        return [float(x), float(y), float(z)]

    def WritePerturbedMdpaFile(self, node_id, perturbation, output_file_name):
        i = _get_position(node_id, self.mdpa_lines)
        coords = self.GetNodalCoordinates(node_id)
        new_coords = [(coords[i_dim] + perturbation[i_dim])
                      for i_dim in range(3)]

        new_line = '{} {:19.10f} {:19.10f} {:19.10f}\n'.format(
            str(node_id), new_coords[0], new_coords[1], new_coords[2])

        new_mdpa_lines = [mdpa_line for mdpa_line in self.mdpa_lines]

        if (perturbation[0] != 0.0 or perturbation[1] != 0.0
                or perturbation[2] != 0.0):
            new_mdpa_lines[i] = new_line

        with open(output_file_name, "w") as new_file:
            new_file.writelines(new_mdpa_lines)


class FiniteDifferenceDragSensitivities():
    def __init__(
            self,
            original_mdpa_file_name,
            output_mdpa_file_name,
            project_parameters_file_name,
            drag_file_name,
            drag_direction,
            working_folder=".",
            sensitivities_file_name="finite_difference_drag_sensitivities.dat"
    ):
        self.original_mdpa_file_name = original_mdpa_file_name
        self.output_mdpa_file_name = output_mdpa_file_name
        self.project_parameters_file_name = project_parameters_file_name
        self.working_folder = working_folder
        self.drag_file_name = drag_file_name
        self.drag_direction = drag_direction
        self.sensitivities_file_name = sensitivities_file_name

        current_path = os.getcwd()
        working_dir = os.path.join(current_path, self.working_folder)
        os.chdir(working_dir)
        self.mdpa_utilities = MDPAUtilities(self.original_mdpa_file_name)
        os.chdir(current_path)

    def _GetTimeAveraged(self):
        time_steps, reactions = _read_drag_file(self.drag_file_name)
        total_drag = 0.0
        for reaction in reactions:
            total_drag += reaction[0] * self.drag_direction[0] + reaction[
                1] * self.drag_direction[1] + reaction[
                    2] * self.drag_direction[2]
        if len(time_steps) > 1:
            delta_time = time_steps[1] - time_steps[0]
            total_drag *= delta_time
        return total_drag

    def ComputeFiniteDifferenceSensitivity(self, node_ids, step_size):
        current_path = os.getcwd()
        working_dir = os.path.join(current_path, self.working_folder)
        os.chdir(working_dir)

        print("Calculating finite difference sensitivity in " + working_dir)

        perturbed_sensitivities = []
        for node_id in node_ids:
            perturbed_drag = []

            # X perturbation
            perturbation = [step_size, 0.0, 0.0]
            self.mdpa_utilities.WritePerturbedMdpaFile(
                node_id, perturbation, self.output_mdpa_file_name)
            print("Running " + str(node_id) + " node x perturbation...")
            self._Solve()
            perturbed_drag.append(self._GetTimeAveraged())

            # Y perturbation
            perturbation = [0.0, step_size, 0.0]
            self.mdpa_utilities.WritePerturbedMdpaFile(
                node_id, perturbation, self.output_mdpa_file_name)
            print("Running " + str(node_id) + " node y perturbation...")
            self._Solve()
            perturbed_drag.append(self._GetTimeAveraged())

            perturbed_sensitivities.append(perturbed_drag)

        # running the unperturbed sensitivity last, then same can be used for adjoints
        perturbation = [0.0, 0.0, 0.0]
        print("Running unperturbed...")
        self.mdpa_utilities.WritePerturbedMdpaFile(1, perturbation,
                                                   self.output_mdpa_file_name)
        self._Solve()
        drag_ref = self._GetTimeAveraged()

        for i, _ in enumerate(perturbed_sensitivities):
            for j in range(2):
                perturbed_sensitivities[i][j] = (
                    perturbed_sensitivities[i][j] - drag_ref) / step_size

        os.chdir(current_path)

        fd_sensitivities = []

        with open(self.sensitivities_file_name, "w") as file_output:
            file_output.write("# Finite difference drag sensitivities\n")
            file_output.write(
                "# NodeId, SHAPE_SENSITIVITY_X, SHAPE_SENSITIVITY_Y\n")
            for node_id, sensitivity in zip(node_ids, perturbed_sensitivities):
                file_output.write("{0:4d},{1:1.16e},{2:1.16e}\n".format(
                    node_id, sensitivity[0], sensitivity[1]))
                fd_sensitivities.append([sensitivity[0], sensitivity[1]])

        return fd_sensitivities

    def _CreateFluidTest(self):
        with open(self.project_parameters_file_name, 'r') as parameter_file:
            project_parameters = Kratos.Parameters(parameter_file.read())
            parameter_file.close()
        test = FluidDynamicsAnalysis(Kratos.Model(), project_parameters)
        return test

    def _Solve(self):
        test = self._CreateFluidTest()
        test.Run()

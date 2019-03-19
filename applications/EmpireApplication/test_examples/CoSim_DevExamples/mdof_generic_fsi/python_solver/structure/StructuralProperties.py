# Class for reading structural properties:
import numpy as np


class StructuralProperties:

    'read defined parameters in project parameter file'

    def __init__(self, ProjectParameters):

        # Geometric properties
        self.type = ProjectParameters["structure_data"]["type"].GetString()
        self.height = ProjectParameters["structure_data"]["height"].GetDouble()
        self.length = ProjectParameters["structure_data"]["length"].GetDouble()
        self.width = ProjectParameters["structure_data"]["width"].GetDouble()
        self.levels = int(
            ProjectParameters["structure_data"]["levels"].GetDouble())

        self.rot_inertia = ProjectParameters[
            "structure_data"]["rot_inertia"].GetDouble()

        # Time integration properties
        self.dt = ProjectParameters["problem_data"]["time_step"].GetDouble()

        # Material properties
        self.density = ProjectParameters[
            "structure_data"]["density"].GetDouble()
        self.mass = ProjectParameters["structure_data"]["mass"].GetDouble()
        self.elast_modulus = ProjectParameters[
            "structure_data"]["elastic_modulus"].GetDouble()

        # Desired eigenfrequencies
        self.eigen_freq_X = ProjectParameters[
            "structure_data"]["eigen_frequencies"][0].GetDouble()
        self.eigen_freq_Y = ProjectParameters[
            "structure_data"]["eigen_frequencies"][1].GetDouble()
        self.eigen_freq_R = ProjectParameters[
            "structure_data"]["eigen_frequencies"][2].GetDouble()

        # Generalized alpha parameters
        self.zeta_X = ProjectParameters["structure_data"]["zeta"][0].GetDouble(
        )
        self.zeta_Y = ProjectParameters["structure_data"]["zeta"][1].GetDouble(
        )
        self.zeta_R = ProjectParameters["structure_data"]["zeta"][2].GetDouble(
        )
        self.rho_inf = ProjectParameters[
            "structure_data"]["rho_inf"].GetDouble()

        # Initial conditions of the structure
        # Displacement
        self.disp_X = np.zeros((self.levels)*2)
        self.disp_Y = np.zeros((self.levels)*2)
        self.disp_R = np.zeros((self.levels))

        # Velocity
        self.vel_X = np.zeros((self.levels)*2)
        self.vel_Y = np.zeros((self.levels)*2)
        self.vel_R = np.zeros((self.levels))

        # Acceleration
        self.acc_X = np.zeros((self.levels)*2)
        self.acc_Y = np.zeros((self.levels)*2)
        self.acc_R = np.zeros((self.levels))

        # Output File
        self.output_filename_X = ProjectParameters[
            "structure_data"]["output_filename_X"].GetString()
        self.output_filename_Y = ProjectParameters[
            "structure_data"]["output_filename_Y"].GetString()
        self.output_filename_R = ProjectParameters[
            "structure_data"]["output_filename_R"].GetString()
        self.output_filename_Result = ProjectParameters[
            "structure_data"]["output_filename_Result"].GetString()

        # FSI parameters
        self.fsi_abs_res = ProjectParameters[
            "FSI_parameters"]["abs_residual"].GetDouble()
        self.fsi_rel_res = ProjectParameters[
            "FSI_parameters"]["rel_residual"].GetDouble()
        self.fsi_relax_coef = ProjectParameters[
            "FSI_parameters"]["relax_coef"].GetDouble()
        self.fsi_max_iter = int(
            ProjectParameters["FSI_parameters"]["max_FSI_iteration"].GetDouble())

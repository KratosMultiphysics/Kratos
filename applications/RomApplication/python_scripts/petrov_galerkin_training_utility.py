# Import Python modules
import json
import numpy as np
from pathlib import Path
from glob import glob

# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.RomApplication as KratosROM
from KratosMultiphysics.RomApplication.randomized_singular_value_decomposition import RandomizedSingularValueDecomposition
import KratosMultiphysics.kratos_utilities as kratos_utils

class PetrovGalerkinTrainingUtility(object):
    """Auxiliary utility for the Petrov Galerkin training.
    This class encapsulates all the functions required for the Petrov Galerkin training.
    These are snapshots collection for the basis Psi used for solving a Petrov Galerkin ROM.
    The snapshots depends on the basis strategy (i.e. Jacobian or Residuals).
    """

    def __init__(self, solver, custom_settings):
        # Validate and assign the HROM training settings
        settings = custom_settings["train_petrov_galerkin"]
        settings.ValidateAndAssignDefaults(self.__GetPetrovGalerkinTrainingDefaultSettings())

        # Auxiliary member variables
        self.solver = solver
        self.time_step_snapshots_matrix_container = [] ## J@Phi or R
        self.echo_level = settings["echo_level"].GetInt()
        self.rom_settings = custom_settings["rom_settings"]
        self.basis_strategy = settings["basis_strategy"].GetString()
        self.include_phi = settings["include_phi"].GetBool()
        self.svd_truncation_tolerance = settings["svd_truncation_tolerance"].GetDouble()
        self.rom_basis_output_name = Path(custom_settings["rom_basis_output_name"].GetString())
        self.rom_basis_output_folder = Path(custom_settings["rom_basis_output_folder"].GetString())

        self.rom_format =  custom_settings["rom_format"].GetString()
        available_rom_format = ["json", "numpy"]
        if self.rom_format not in available_rom_format:
            err_msg = "Provided \'rom format\' is {}. Available options are {}.".format(self.rom_format, available_rom_format)
            raise Exception(err_msg)


    def AppendCurrentStepProjectedSystem(self):
        # Get the computing model part from the solver implementing the problem physics
        computing_model_part = self.solver.GetComputingModelPart()

        # If not created yet, create the ROM residuals utility
        # Note that this ensures that the residuals utility is created in the first residuals append call
        # If not, it might happen that the solver scheme is created by the ROM residuals call rather than by the solver one
        if not hasattr(self, '__rom_residuals_utility'):
            self.__rom_residuals_utility = KratosROM.RomResidualsUtility(
                computing_model_part,
                self.rom_settings,
                self.solver._GetScheme())

            if self.echo_level > 0 : KratosMultiphysics.Logger.PrintInfo("PetrovGalerkinTrainingUtility","RomResidualsUtility created.")

        # Generate the matrix of residuals or projected Jacobians.
        if self.basis_strategy=="jacobian":
            snapshots_matrix = self.__rom_residuals_utility.GetProjectedGlobalLHS()
            if self.echo_level > 0 : KratosMultiphysics.Logger.PrintInfo("PetrovGalerkinTrainingUtility","Generated matrix of projected Jacobian.")
        elif self.basis_strategy=="residuals":
            snapshots_matrix = []
            files_to_read_and_delete = glob('*.res.mm')#TODO: Stop writing to disk.
            for to_erase_file in files_to_read_and_delete:
                non_converged_iteration_snapshot = KratosMultiphysics.Vector()
                KratosMultiphysics.ReadMatrixMarketVector(to_erase_file, non_converged_iteration_snapshot)
                snapshots_matrix.append(non_converged_iteration_snapshot)
            snapshots_matrix = np.array(snapshots_matrix).T
            if self.echo_level > 0 : KratosMultiphysics.Logger.PrintInfo("PetrovGalerkinTrainingUtility","Generating matrix of residuals.")
        else:
            err_msg = "\'self.basis_strategy\' is not available. Select either 'jacobian' or 'residuals'."
            raise Exception(err_msg)

        np_snapshots_matrix = np.array(snapshots_matrix, copy=False)
        self.time_step_snapshots_matrix_container.append(np_snapshots_matrix)

    def CalculateAndSaveBasis(self, snapshots_matrix = None):
        # Calculate the new basis and save
        snapshots_basis = self.__CalculateResidualBasis(snapshots_matrix)
        self.__AppendNewBasisToRomParameters(snapshots_basis)


    @classmethod
    def __GetPetrovGalerkinTrainingDefaultSettings(cls):
        default_settings = KratosMultiphysics.Parameters("""{
                "train": false,
                "basis_strategy": "residuals",
                "include_phi": false,
                "svd_truncation_tolerance": 1.0e-6,
                "echo_level": 0
        }""")
        return default_settings

    def __CalculateResidualBasis(self, snapshots_matrix):
        if snapshots_matrix is None:
            snapshots_matrix = self._GetSnapshotsMatrix()
        u_left,s_left,_,_ = RandomizedSingularValueDecomposition(COMPUTE_V=False).Calculate(
            snapshots_matrix,
            self.svd_truncation_tolerance)
        del(snapshots_matrix)
        #Include Phi in the span of Psi
        if self.include_phi:
            #Read Galerkin modes from RomParameters.json
            u_right = self.__GetGalerkinBasis()
            u , _= np.linalg.qr(np.c_[u_right,u_left])
        else:
            u = u_left

        return u

    def __AppendNewBasisToRomParameters(self, u):
        petrov_galerkin_number_of_rom_dofs= np.shape(u)[1]
        n_nodal_unknowns = len(self.rom_settings["nodal_unknowns"].GetStringArray())
        petrov_galerkin_nodal_modes = {}
        computing_model_part = self.solver.GetComputingModelPart()

        if self.rom_format == "json":
            i = 0
            for node in computing_model_part.Nodes:
                petrov_galerkin_nodal_modes[node.Id] = u[i:i+n_nodal_unknowns].tolist()
                i += n_nodal_unknowns

        elif self.rom_format == "numpy":
            # Storing Petrov-Galerkin modes in Numpy format
            np.save(self.rom_basis_output_folder / "LeftBasisMatrix.npy", u)

        with (self.rom_basis_output_folder / self.rom_basis_output_name).with_suffix('.json').open('r') as f:
            updated_rom_parameters = json.load(f)
            updated_rom_parameters["rom_settings"]["petrov_galerkin_number_of_rom_dofs"] = petrov_galerkin_number_of_rom_dofs
            updated_rom_parameters["petrov_galerkin_nodal_modes"] = petrov_galerkin_nodal_modes

        with (self.rom_basis_output_folder / self.rom_basis_output_name).with_suffix('.json').open('w') as f:
            json.dump(updated_rom_parameters, f, indent=4)

        if self.echo_level > 0 : KratosMultiphysics.Logger.PrintInfo("PetrovGalerkinTrainingUtility","\'RomParameters.json\' file updated with HROM weights.")

    def __GetPrettyFloat(self, number):
        float_format = "{:.12f}"
        pretty_number = float(float_format.format(number))
        return pretty_number

    def __GetGalerkinBasis(self):
        with open(self.rom_basis_output_folder / self.rom_basis_output_name.with_suffix(".json"), 'r') as f:
            galerkin_rom_parameters = json.load(f)
            N_Dof_per_node = len(galerkin_rom_parameters["rom_settings"]["nodal_unknowns"])
            N_nodes = len(galerkin_rom_parameters["nodal_modes"])
            N_Dofs = int(N_Dof_per_node*N_nodes)
            N_Dofs_rom = galerkin_rom_parameters["rom_settings"]["number_of_rom_dofs"]
            u = np.zeros((N_Dofs,N_Dofs_rom))
            counter_in = 0
            for key in galerkin_rom_parameters["nodal_modes"].keys():
                counter_fin = counter_in + N_Dof_per_node
                u[counter_in:counter_fin,:] = np.array(galerkin_rom_parameters["nodal_modes"][key])
                counter_in = counter_fin

        return u


    def _GetSnapshotsMatrix(self):
        # Set up the snapshots matrix for new basis
        n_steps = len(self.time_step_snapshots_matrix_container)
        snapshots_matrix = self.time_step_snapshots_matrix_container[0]
        for i in range(1,n_steps):
            del self.time_step_snapshots_matrix_container[0] # Avoid having two matrices, numpy does not concatenate references.
            snapshots_matrix = np.c_[snapshots_matrix,self.time_step_snapshots_matrix_container[0]]
        return snapshots_matrix


    @classmethod
    def __OrthogonalProjector(self, A, B):
        # A - B @(B.T @ A)
        BtA = B.T@A
        A -= B @ BtA
        return A
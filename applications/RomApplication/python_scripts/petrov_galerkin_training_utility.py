# Import Python modules
import json
import numpy as np
import os
import scipy.linalg

# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.RomApplication as KratosROM
from KratosMultiphysics.RomApplication.randomized_singular_value_decomposition import RandomizedSingularValueDecomposition
import KratosMultiphysics.kratos_utilities as kratos_utils

def GetFilePath(fileName):
    return os.path.join(os.getcwd(), fileName)

class PetrovGalerkinTrainingUtility(object):
    """Auxiliary utility for the Petrov Galerkin training.
    This class encapsulates all the functions required for the Petrov Galerkin training.
    These are snapshots collection for the basis Psi used for solving a Petrov Galerkin ROM.
    The snapshots depends on the basis strategy (i.e. Projected_System or Assembled_Residuals).
    """

    def __init__(self, solver, custom_settings):
        # Validate and assign the HROM training settings
        settings = custom_settings["train_petrov_galerkin"]
        settings.ValidateAndAssignDefaults(self.__GetPetrovGalerkinTrainingDefaultSettings())

        # Auxiliary member variables
        self.solver = solver
        self.time_step_snapshots_matrix_container = [] ## K@Phi
        self.echo_level = settings["echo_level"].GetInt()
        self.rom_settings = custom_settings["rom_settings"]
        # self.basis_strategy = "Assembled_Residuals" if self.rom_settings["solving_strategy"].GetString()=="QR" else settings["basis_strategy"].GetString() #FIXME: Either stop writing to memory or make sure that we can erase after SVD. When solving with QR, force this option to erase the matrices that are written to build the assembled residual option.
        self.basis_strategy = settings["basis_strategy"].GetString()
        self.include_phi = settings["include_phi"].GetBool()
        self.svd_truncation_tolerance = settings["svd_truncation_tolerance"].GetDouble()

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

        # Generate the matrix of residuals
        if self.echo_level > 0 : KratosMultiphysics.Logger.PrintInfo("PetrovGalerkinTrainingUtility","Generating matrix of projected system.")
        
        if self.basis_strategy=="Projected_System":
            snapshots_matrix = self.__rom_residuals_utility.GetProjectedGlobalLHS()
        elif self.basis_strategy=="Assembled_Residuals":
            snapshots_matrix = []
            current_time = self.__GetPrettyFloat(computing_model_part.ProcessInfo[KratosMultiphysics.TIME])
            if current_time.is_integer():
                current_time = int(current_time) #FIXME: I dont like this, c++ writes step with no decimals when is integer and python does. (i.e. c++ step 1=1, python step 1=1.0)
            number_of_iterations = computing_model_part.ProcessInfo[KratosMultiphysics.NL_ITERATION_NUMBER]
            # This loop reads NON-CONVERGED assembled residuals that had been written into disk
            for i in range(number_of_iterations+1):
                non_converged_iteration_snapshot = KratosMultiphysics.Vector()
                snapshot_name = "R_"+str(current_time)+"_"+str(i)+".res.mm"
                KratosMultiphysics.ReadMatrixMarketVector(snapshot_name, non_converged_iteration_snapshot)
                snapshots_matrix.append(non_converged_iteration_snapshot)
                # kratos_utils.DeleteFileIfExisting(GetFilePath(snapshot_name))
            write_to_disk_converged_residual = False
            converged_snapshot = self.__rom_residuals_utility.GetAssembledResiduals(write_to_disk_converged_residual)
            snapshots_matrix.append(converged_snapshot)
            snapshots_matrix = np.array(snapshots_matrix).T
        else:
            err_msg = "\'self.basis_strategy\' is not available. Select either Projected_System or Assembled_Residuals."
            raise Exception(err_msg)

        np_snapshots_matrix = np.array(snapshots_matrix, copy=False)
        self.time_step_snapshots_matrix_container.append(np_snapshots_matrix)

    def CalculateAndSaveBasis(self):
        # Calculate the new basis and save
        snapshots_basis = self.__CalculateResidualBasis()
        self.__AppendNewBasisToRomParameters(snapshots_basis)
    

    @classmethod
    def __GetPetrovGalerkinTrainingDefaultSettings(cls):
        default_settings = KratosMultiphysics.Parameters("""{
                "train": false,
                "basis_strategy": "Assembled_Residuals",
                "include_phi": true,
                "svd_truncation_tolerance": 1.0e-4,
                "echo_level": 0
        }""")
        return default_settings

    def __CalculateResidualBasis(self):
        # Set up the snapshots matrix for new basis
        n_steps = len(self.time_step_snapshots_matrix_container)
        snapshots_matrix = self.time_step_snapshots_matrix_container[0]
        for i in range(1,n_steps):
            del self.time_step_snapshots_matrix_container[0] # Avoid having two matrices, numpy does not concatenate references.
            snapshots_matrix = np.c_[snapshots_matrix,self.time_step_snapshots_matrix_container[0]]
        u_left,s_left,_,_ = RandomizedSingularValueDecomposition(COMPUTE_V=False).Calculate(
            snapshots_matrix,
            self.svd_truncation_tolerance)
        del(snapshots_matrix)
        if self.include_phi:
            #Read Galerkin modes from RomParameters.json
            u_right = self.__GetGalerkinBasis()
            u , _= np.linalg.qr(np.c_[u_right,u_left])
        else:
            u = u_left
        # u,s,_,_ = RandomizedSingularValueDecomposition(COMPUTE_V=False).Calculate(
        #     np.c_[u_right,u_left],
        #     self.svd_truncation_tolerance)
        # del(u_right)
        # del(u_left)
        # u,R, _ = scipy.linalg.qr(np.c_[u_right,u_left],pivoting=True)
        # del(u_right)
        # del(u_left)
        # tol = np.sqrt(u.shape[1])*self.svd_truncation_tolerance # Norm of u is equal to the sqrt of the number of nodes because of its orthogonality.
        # rank = R.shape[1]
        # for i in range(rank):
        #     if abs(R[i,i])<tol:
        #         rank = i
        #         break
        # u = u[:,:rank]

        return u
    
    def __AppendNewBasisToRomParameters(self, u):
        petrov_galerkin_number_of_rom_dofs= np.shape(u)[1]
        n_nodal_unknowns = len(self.rom_settings["nodal_unknowns"].GetStringArray())
        petrov_galerkin_nodal_modes = {}
        computing_model_part = self.solver.GetComputingModelPart()
        i = 0
        for node in computing_model_part.Nodes:
            petrov_galerkin_nodal_modes[node.Id] = u[i:i+n_nodal_unknowns].tolist()
            i += n_nodal_unknowns

        with open('RomParameters.json','r') as f:
            updated_rom_parameters = json.load(f)
            updated_rom_parameters["rom_settings"]["petrov_galerkin_number_of_rom_dofs"] = petrov_galerkin_number_of_rom_dofs
            updated_rom_parameters["petrov_galerkin_nodal_modes"] = petrov_galerkin_nodal_modes

        with open('RomParameters.json','w') as f:
            json.dump(updated_rom_parameters, f, indent = 4)

        if self.echo_level > 0 : KratosMultiphysics.Logger.PrintInfo("PetrovGalerkinTrainingUtility","\'RomParameters.json\' file updated with HROM weights.")

    def __GetPrettyFloat(self, number):
        float_format = "{:.12f}"
        pretty_number = float(float_format.format(number))
        return pretty_number
    
    def __GetGalerkinBasis(self):
        with open('RomParameters.json','r') as f:
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
    
    def __OrthogonalProjector(self, A, B):
        # A - B @(B.T @ A)
        BtA = B.T@A
        A -= B @ BtA
        return A
import json

import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsApplication
import KratosMultiphysics.RomApplication as romapp

from KratosMultiphysics.RomApplication import networks.check_gradient as check_gradient
from KratosMultiphysics.RomApplication import networks.shallow_autoencoder as shallow_autoencoder
from KratosMultiphysics.RomApplication import networks.cluster_autoencoder as cluster_autoencoder
from KratosMultiphysics.RomApplication import networks.gradient_shallow_jacobi as gradient_shallow
from KratosMultiphysics.RomApplication import python_solvers_wrapper_rom as solver_wrapper

from KratosMultiphysics.FluidDynamicsApplication.fluid_dynamics_analysis import FluidDynamicsAnalysis

class FluidDynamicsAnalysisNNROM(FluidDynamicsAnalysis):

    def __init__(self,model,project_parameters):
        KratosMultiphysics.Logger.PrintWarning('\x1b[1;31m[DEPRECATED CLASS] \x1b[0m',"\'FluidDynamicsAnalysisNNROM\'", "class is deprecated. Use the \'RomAnalysis\' one instead.")
        super().__init__(model,project_parameters)

    ### NN ##
    def _CreateNNModel(self):
        # List of files to read from
        data_inputs = [
            "hdf5_output/result_10.h5",
            "hdf5_output/result_15.h5",
            "hdf5_output/result_20.h5",
            "hdf5_output/result_25.h5",
            "hdf5_output/result_30.h5",
        ]

        # List of variables to run be used (VELOCITY, PRESSURE, etc...)
        S = kratos_io.build_snapshot_grid(
            data_inputs, 
            [
                "VELOCITY"
            ]
        )

        # Some configuration
        config = {
            "load_model":       False,
            "train_model":      True,
            "save_model":       False,
            "print_results":    True,
            "use_reduced":      True,
        }

        # Select the netwrok to use
        self.kratos_network = gradient_shallow_ae.GradientShallow()

        print("=== Calculating Randomized Singular Value Decomposition ===")
        with contextlib.redirect_stdout(None):
        # Note: we redirect the output here because there is no option to reduce the log level of this method.
            U,sigma,_,error = RandomizedSingularValueDecomposition().Calculate(S,0)

            SPri = U @ U.T @ S

        # Select the reduced snapshot or the full input
        if config["use_reduced"]:
            SReduced = U.T @ S
        else:
            SReduced = S

        kratos_network.calculate_data_limits(SReduced)

        # Set the properties for the clusters
        num_clusters=1                      # Number of different bases chosen
        num_cluster_col=5                   # If I use num_cluster_col = num_variables result should be exact.
        num_encoding_var=num_cluster_col

        # Load the model or train a new one. TODO: Fix custom_loss not being saved correctly
        if config["load_model"]:
            self.autoencoder = tf.keras.models.load_model(kratos_network.model_name, custom_objects={"custom_loss":custom_loss})
        else:
            self.autoencoder, self.autoencoder_err = kratos_network.define_network(SReduced, custom_loss, num_encoding_var)

        # Prepare the data
        SReducedNormalized = (SReduced - kratos_network.data_min) / (kratos_network.data_max - kratos_network.data_min)

        # Obtain Clusters
        # -- Currenctly this function is quite useless because we are always using the classical SVD with 1 cluster.
        print("=== Calculating Cluster Bases ===")
        with contextlib.redirect_stdout(None):
            CB, kmeans_object = clustering.calcualte_snapshots_with_columns(
                snapshot_matrix=SReducedNormalized,
                number_of_clusters=num_clusters,
                number_of_columns_in_the_basis=num_cluster_col,
                truncation_tolerance=1e-5
            )
        
        # Obtain the reduced representations (that would be our G's)
        q, s, c, b = {}, {}, {}, {}
        total_respresentations = 0

        for i in range(num_clusters):
            q[i] = CB[i] @ CB[i].T @ SReduced[:,kmeans_object.labels_==i]
            b[i] = CB[i] @ CB[i].T
            c[i] = np.array([i for _ in range(q[i].shape[1])])
            s[i] = SReduced[:,kmeans_object.labels_==i]

        print("Total number of representations:", total_respresentations)

        temp_size = num_cluster_col
        grad_size = num_cluster_col

        qc = np.empty(shape=(data_rows,0))
        sc = np.empty(shape=(data_rows,0))
        cc = np.empty(shape=(0))

        for i in range(num_clusters):
            qc = np.concatenate((qc, q[i]), axis=1)
            sc = np.concatenate((sc, s[i]), axis=1)
            cc = np.concatenate((cc, c[i]), axis=0)

        nqc = kratos_network.normalize_data(qc)
        nsc = kratos_network.normalize_data(sc)

        print(f"{cc.shape=}")

        if config["train_model"]:
            self.autoencoder.set_m_grad(b)
            self.autoencoder.set_g_weight(1)
            kratos_network.train_network(self.autoencoder, nsc, cc, len(data_inputs), 400)
            print("Model Initialized")

        if config["save_model"]:
            self.autoencoder.save(kratos_network.model_name)

    #### Internal functions ####
    def _CreateSolver(self):
        """ Create the Solver (and create and import the ModelPart if it is not alread in the model) """
        ## Solver construction
        with open('RomParameters.json') as rom_parameters:
            rom_settings = KratosMultiphysics.Parameters(rom_parameters.read())
            self.project_parameters["solver_settings"].AddValue("rom_settings", rom_settings["rom_settings"])
        return solver_wrapper.CreateSolverByParameters(self.model, self.project_parameters["solver_settings"],self.project_parameters["problem_data"]["parallel_type"].GetString())

    def _GetSimulationName(self):
        return "::[ROM Simulation]:: "

    def ModifyAfterSolverInitialize(self):
        """Here is where the ROM_BASIS is imposed to each node"""
        super().ModifyAfterSolverInitialize()
        computing_model_part = self._solver.GetComputingModelPart()
        with open('RomParameters.json') as f:
            data = json.load(f)
            nodal_dofs = len(data["rom_settings"]["nodal_unknowns"])
            nodal_modes = data["nodal_modes"]
            counter = 0
            rom_dofs= self.project_parameters["solver_settings"]["rom_settings"]["number_of_rom_dofs"].GetInt()
            for node in computing_model_part.Nodes:
                aux = KratosMultiphysics.Matrix(nodal_dofs, rom_dofs)
                for j in range(nodal_dofs):
                    Counter=str(node.Id)
                    for i in range(rom_dofs):
                        aux[j,i] = nodal_modes[Counter][j][i]
                node.SetValue(romapp.ROM_BASIS, aux ) # ROM basis
                counter+=1

    def InitializeSolutionStep(self):
        """This function performs all the required operations that should be executed
        (for each step) BEFORE solving the solution step.
        """

        # Resotre with the coarse model
        fine_model = self.kratos_network.predict_snapshot(autoencoder, nsc)

        # Predict with the fine model
        coarse_model = self.kratos_network.predict_snapshot(autoencoder, nsc)

        # Check the gradients:
        grad_pred = np.array([coarse_model.T])

        grad_chk = check_gradient.GradientCheck(
            kratos_network.encoder_model, 
            kratos_network.get_gradients(
                kratos_network.encoder_model, 
                [kratos_network.encoder_model],
                grad_pred,
                grad_expt
            )
        )

        super().InitializeSolutionStep()

    
    def FinalizeSolutionStep(self):
        """This function performs all the required operations that should be executed
        (for each step) AFTER solving the solution step.
        """
        super().FinalizeSolutionStep()
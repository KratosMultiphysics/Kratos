import os
import sys
import logging

# logging.disable(logging.WARNING)
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3' 

import h5py
import numpy as np

import keras
import tensorflow as tf

from keras import layers
from itertools import repeat

import importlib
import contextlib

import KratosMultiphysics
import KratosMultiphysics.RomApplication as KratosROM

from KratosMultiphysics.FluidDynamicsApplication.fluid_dynamics_analysis import FluidDynamicsAnalysis

from KratosMultiphysics.RomApplication.networks import gradient_shallow as gradient_shallow_ae
from KratosMultiphysics.RomApplication import kratos_io as kratos_io
from KratosMultiphysics.RomApplication import python_solvers_wrapper_rom as solver_wrapper

from KratosMultiphysics.RomApplication import new_python_solvers_wrapper_rom
from KratosMultiphysics.RomApplication.hrom_training_utility import HRomTrainingUtility
from KratosMultiphysics.RomApplication.calculate_rom_basis_output_process import CalculateRomBasisOutputProcess

from KratosMultiphysics.RomApplication.randomized_singular_value_decomposition import RandomizedSingularValueDecomposition

def custom_loss(y_true, y_pred):
    y_diff = y_true-y_pred
    y_diff = y_diff ** 2

    return  y_diff

def CreateRomAnalysisInstance(cls, global_model, parameters):
    class RomAnalysis(cls):

        def __init__(self,global_model, parameters):
            super().__init__(global_model, parameters)

        ### NN ##
        def _CreateModel(self):
            # List of files to read from
            data_inputs = [
                "hdf5_bases/result_10.h5",
                "hdf5_bases/result_15.h5",
                "hdf5_bases/result_20.h5",
                "hdf5_bases/result_25.h5",
                "hdf5_bases/result_30.h5",
            ]

            # List of variables to run be used (VELOCITY, PRESSURE, etc...)
            S = kratos_io.build_snapshot_grid(
                data_inputs, 
                [
                    "VELOCITY",
                    "PRESSURE"
                ]
            )

            # Some configuration
            config = {
                "load_model":       True,
                "train_model":      False,
                "save_model":       False,
                "print_results":    False,
                "use_reduced":      True,
            }

            # Select the netwrok to use
            self.kratos_network = gradient_shallow_ae.GradientShallow()

            print("=== Calculating Randomized Singular Value Decomposition ===")
            with contextlib.redirect_stdout(None):
            # Note: we redirect the output here because there is no option to reduce the log level of this method.
                self.U,sigma,_,error = RandomizedSingularValueDecomposition().Calculate(S,1e-5)

                SPri = self.U @ self.U.T @ S

            # Select the reduced snapshot or the full input
            if config["use_reduced"]:
                SReduced = self.U.T @ S
            else:
                SReduced = S

            self.kratos_network.calculate_data_limits(SReduced)

            # Set the properties for the clusters
            num_clusters=1                      # Number of different bases chosen
            num_cluster_col=5                   # If I use num_cluster_col = num_variables result should be exact.
            num_encoding_var=num_cluster_col

            # Load the model or train a new one. TODO: Fix custom_loss not being saved correctly
            if config["load_model"]:
                self.autoencoder = tf.keras.models.load_model(self.kratos_network.model_name, custom_objects={
                    "custom_loss":custom_loss,
                    "set_m_grad":gradient_shallow_ae.GradModel2.set_m_grad
                })
                self.kratos_network.encoder_model = tf.keras.models.load_model(self.kratos_network.model_name+"_encoder", custom_objects={
                    "custom_loss":custom_loss,
                    "set_m_grad":gradient_shallow_ae.GradModel2.set_m_grad
                })
                self.kratos_network.decoder_model = tf.keras.models.load_model(self.kratos_network.model_name+"_decoder", custom_objects={
                    "custom_loss":custom_loss,
                    "set_m_grad":gradient_shallow_ae.GradModel2.set_m_grad
                })
            else:
                self.autoencoder, self.autoencoder_err = self.kratos_network.define_network(SReduced, custom_loss, num_encoding_var)

            # Prepare the data
            SReducedNormalized = (SReduced - self.kratos_network.data_min) / (self.kratos_network.data_max - self.kratos_network.data_min)

            # Obtain Clusters
            nsc = self.kratos_network.normalize_data(S)

            if config["train_model"]:
                self.autoencoder.set_m_grad(b)
                self.autoencoder.set_g_weight(1)
                self.kratos_network.train_network(self.autoencoder, nsc, cc, len(data_inputs), 400)
                print("Model Initialized")

            if config["save_model"]:
                self.autoencoder.save(self.kratos_network.model_name)

        def _CreateSolver(self):
            """ Create the Solver (and create and import the ModelPart if it is not alread in the model) """

            # Get the ROM settings from the RomParameters.json input file
            # with open('RomParameters.json') as rom_parameters:
            #     self.rom_parameters = KratosMultiphysics.Parameters(rom_parameters.read())

            # Set the ROM settings in the "solver_settings" of the solver introducing the physics
            # self.project_parameters["solver_settings"].AddValue("rom_settings", self.rom_parameters["rom_settings"])

            # HROM operations flags
            self.rom_basis_process_list_check = False
            self.rom_basis_output_process_check = False
            self.run_hrom = False # self.rom_parameters["run_hrom"].GetBool() if self.rom_parameters.Has("run_hrom") else False
            self.train_hrom = False # self.rom_parameters["train_hrom"].GetBool() if self.rom_parameters.Has("train_hrom") else False
            if self.run_hrom and self.train_hrom:
                # Check that train an run HROM are not set at the same time
                err_msg = "\'run_hrom\' and \'train_hrom\' are both \'true\'. Select either training or running (if training has been already done)."
                raise Exception(err_msg)

            # Create the ROM solver
            return new_python_solvers_wrapper_rom.CreateSolver(
                self.model,
                self.project_parameters)

        def _GetListOfProcesses(self):
            # Get the already existent processes list
            list_of_processes = super()._GetListOfProcesses()

            # Check if there is any instance of ROM basis output
            if self.rom_basis_process_list_check:
                for process in list_of_processes:
                    if isinstance(process, KratosROM.calculate_rom_basis_output_process.CalculateRomBasisOutputProcess):
                        warn_msg = "\'CalculateRomBasisOutputProcess\' instance found in ROM stage. Basis must be already stored in \'RomParameters.json\'. Removing instance from processes list."
                        KratosMultiphysics.Logger.PrintWarning("RomAnalysis", warn_msg)
                        list_of_processes.remove(process)
                self.rom_basis_process_list_check = False

            return list_of_processes

        def _GetListOfOutputProcesses(self):
            # Get the already existent output processes list
            list_of_output_processes = super()._GetListOfOutputProcesses()

            # Check if there is any instance of ROM basis output
            if self.rom_basis_output_process_check:
                for process in list_of_output_processes:
                    if isinstance(process, KratosROM.calculate_rom_basis_output_process.CalculateRomBasisOutputProcess):
                        warn_msg = "\'CalculateRomBasisOutputProcess\' instance found in ROM stage. Basis must be already stored in \'RomParameters.json\'. Removing instance from output processes list."
                        KratosMultiphysics.Logger.PrintWarning("RomAnalysis", warn_msg)
                        list_of_output_processes.remove(process)
                self.rom_basis_output_process_check = False

            return list_of_output_processes

        def _GetSimulationName(self):
            return "::[ROM Simulation]:: "

        def InitializeSolutionStep(self):
            # Generate the input snapshot

            sc = []
            computing_model_part = self._GetSolver().GetComputingModelPart().GetRootModelPart()
            for node in computing_model_part.Nodes:
                sc.append(node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_X))
                sc.append(node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_Y))
                sc.append(node.GetSolutionStepValue(KratosMultiphysics.PRESSURE))
            sc = self.U.T @ np.array([sc]).T

            # nsc = self.kratos_network.normalize_data(sc)
            nsc = sc
            NSPredict_e = self.kratos_network.predict_snapshot(self.autoencoder, nsc)
            NSPredict_d = self.kratos_network.predict_snapshot(self.kratos_network.encoder_model, nsc)

            grad_enc_input = np.array([nsc.T[0]], dtype="float32")
            grad_dec_input = np.array([NSPredict_d.T[0]])

            nn_nodal_modes_enc = self.kratos_network.get_gradients(
                self.kratos_network.encoder_model, 
                [self.kratos_network.encoder_model],
                grad_enc_input,
                None
            )

            nn_nodal_modes_dec = self.kratos_network.get_gradients(
                self.kratos_network.decoder_model, 
                [self.kratos_network.decoder_model],
                grad_dec_input,
                None
            )

            print("U        shapes:", self.U.shape)
            print("gradient shapes:", nn_nodal_modes_enc.T.shape)
            print("gradient shapes:", nn_nodal_modes_dec.T.shape)

            # Set the nodal ROM basis
            nodal_modes_e = self.U @ nn_nodal_modes_enc.T 
            nodal_modes_d = self.U @ nn_nodal_modes_dec    # This is the same as (nn_nodal_modes_dec.T @ self.U.T).T
            # nodal_modes_d = nodal_modes_d.T

            # print("nodal_modes shapes:", nodal_modes.shape)

            nodal_dofs = 3
            rom_dofs = 5

            aux_e = KratosMultiphysics.Matrix(nodal_dofs, rom_dofs)
            aux_d = KratosMultiphysics.Matrix(nodal_dofs, rom_dofs)
            for node in computing_model_part.Nodes:
                node_id = node.Id-1
                for j in range(nodal_dofs):
                    for i in range(rom_dofs):
                        aux_e[j,i] = nodal_modes_e[node_id*nodal_dofs+j][i]
                        aux_d[j,i] = nodal_modes_d[node_id*nodal_dofs+j][i]
                node.SetValue(KratosROM.ROM_BASIS,     aux_e)
                node.SetValue(KratosROM.ROM_BASIS_DEC, aux_e)

            super().InitializeSolutionStep()

        def ModifyAfterSolverInitialize(self):
            """Here is where the ROM_BASIS is imposed to each node"""
            super().ModifyAfterSolverInitialize()

            # Get the model part where the ROM is to be applied
            computing_model_part = self._GetSolver().GetComputingModelPart().GetRootModelPart()
            # computing_model_part = self._GetSolver().GetComputingModelPart()

            # Initialize NN
            self._CreateModel()

            # # Set ROM basis
            # nodal_modes = self.rom_parameters["nodal_modes"]
            # nodal_dofs = len(self.project_parameters["solver_settings"]["rom_settings"]["nodal_unknowns"].GetStringArray())
            # rom_dofs = self.project_parameters["solver_settings"]["rom_settings"]["number_of_rom_dofs"].GetInt()

            # # Set the nodal ROM basis
            # aux = KratosMultiphysics.Matrix(nodal_dofs, rom_dofs)
            # for node in computing_model_part.Nodes:
            #     node_id = str(node.Id)
            #     for j in range(nodal_dofs):
            #         for i in range(rom_dofs):
            #             aux[j,i] = nodal_modes[node_id][j][i].GetDouble()
            #     node.SetValue(KratosROM.ROM_BASIS, aux)

            # # Check for HROM stages
            # if self.train_hrom:
            #     # Create the training utility to calculate the HROM weights
            #     self.__hrom_training_utility = HRomTrainingUtility(
            #         self._GetSolver(),
            #         self.rom_parameters)
            # elif self.run_hrom:
            #     # Set the HROM weights in elements and conditions
            #     hrom_weights_elements = self.rom_parameters["elements_and_weights"]["Elements"]
            #     for key,value in zip(hrom_weights_elements.keys(), hrom_weights_elements.values()):
            #         computing_model_part.GetElement(int(key)+1).SetValue(KratosROM.HROM_WEIGHT, value.GetDouble()) #FIXME: FIX THE +1

            #     hrom_weights_condtions = self.rom_parameters["elements_and_weights"]["Conditions"]
            #     for key,value in zip(hrom_weights_condtions.keys(), hrom_weights_condtions.values()):
            #         computing_model_part.GetCondition(int(key)+1).SetValue(KratosROM.HROM_WEIGHT, value.GetDouble()) #FIXME: FIX THE +1

        def FinalizeSolutionStep(self):
            # Call the HROM training utility to append the current step residuals
            # Note that this needs to be done prior to the other processes to avoid unfixing the BCs
            if self.train_hrom:
                self.__hrom_training_utility.AppendCurrentStepResiduals()

            # #FIXME: Make this optional. This must be a process
            # # Project the ROM solution onto the visualization modelparts
            # if self.run_hrom:
            #     model_part_name = self._GetSolver().settings["model_part_name"].GetString()
            #     visualization_model_part = self.model.GetModelPart("{}.{}Visualization".format(model_part_name, model_part_name))
            #     KratosROM.RomAuxiliaryUtilities.ProjectRomSolutionIncrementToNodes(
            #         self.rom_parameters["rom_settings"]["nodal_unknowns"].GetStringArray(),
            #         visualization_model_part)

            # This calls the physics FinalizeSolutionStep (e.g. BCs)
            super().FinalizeSolutionStep()

        def Finalize(self):
            # This calls the physics Finalize
            super().Finalize()

            # Once simulation is completed, calculate and save the HROM weights
            if self.train_hrom:
                self.__hrom_training_utility.CalculateAndSaveHRomWeights()
                self.__hrom_training_utility.CreateHRomModelParts()

    return RomAnalysis(global_model, parameters)

class FluidDynamicsAnalysisWithFlush(FluidDynamicsAnalysis):

    def __init__(self, model, project_parameters):
        super().__init__(model, project_parameters)
        self.w = 1 # original narrowing size
        time_step_size = self.project_parameters["solver_settings"]["fluid_solver_settings"]["time_stepping"]["time_step"].GetDouble()
        self.delta_w = 0.025 * time_step_size # this ensures to obtain a maximum narrowing size of 2.9 and a minimum of 0.1
        self.control_point = 538 #a node around the middle of the geometry to capture the bufurcation
        self.velocity_y_at_control_point = []
        self.narrowing_width = []

    def MovePart(self, part_name, jacobian, centering_vector, extra_centering):
        x_original = []
        y_original = []
        # first loop
        for node in self.model.GetModelPart(f"FluidModelPart.{part_name}").Nodes:
            if not node.IsFixed(KratosMultiphysics.MESH_DISPLACEMENT_X):
                x_original.append(node.X0)
            if not node.IsFixed(KratosMultiphysics.MESH_DISPLACEMENT_Y):
                y_original.append(node.Y0)
        x_original = np.array(x_original).reshape(1,-1)
        y_original = np.array(y_original).reshape(1,-1)
        matrix_of_coordinates = np.r_[x_original, y_original]
        modified_matrix_of_coordinates = np.linalg.inv(jacobian) @ (matrix_of_coordinates - centering_vector)
        modified_matrix_of_coordinates += centering_vector + extra_centering #re-locating
        # second loop
        i = 0
        for node in self.model.GetModelPart(f"FluidModelPart.{part_name}").Nodes:
            if not node.IsFixed(KratosMultiphysics.MESH_DISPLACEMENT_X):
                x_disp = modified_matrix_of_coordinates[0,i] - node.X0
                node.SetSolutionStepValue(KratosMultiphysics.MESH_DISPLACEMENT_X,0, x_disp )
            if not node.IsFixed(KratosMultiphysics.MESH_DISPLACEMENT_Y):
                y_disp = modified_matrix_of_coordinates[1,i] - node.Y0
                node.SetSolutionStepValue(KratosMultiphysics.MESH_DISPLACEMENT_Y,0, y_disp )
                i +=1
            node.Fix(KratosMultiphysics.MESH_DISPLACEMENT_X)
            node.Fix(KratosMultiphysics.MESH_DISPLACEMENT_Y)


    def StoreBifurcationData(self):
        node =  self.model.GetModelPart("FluidModelPart").GetNode(self.control_point)
        self.velocity_y_at_control_point.append(node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_Y))
        self.narrowing_width.append(self.w)


    def InitializeSolutionStep(self):
        print("Print Executing this")
        super().InitializeSolutionStep()

        if self.time>10.0: # start modifying narrowing from 10 seconds onwards

            if self.time<86: #expansion
                self.w += self.delta_w
            elif self.time>86 and self.time<197: #contraction
                self.w -= self.delta_w
            elif  self.time>197: #expansion again
                self.w += self.delta_w

            print('narrowing size: ', self.w)

            #############################
            ####    FREE ALL NODES   ####
            #############################
            for node in self.model.GetModelPart("FluidModelPart").Nodes:
                node.Free(KratosMultiphysics.MESH_DISPLACEMENT_X)
                node.Free(KratosMultiphysics.MESH_DISPLACEMENT_Y)

            #############################
            #### FIXING OUTSIDE PART ####
            #############################
            for node in self.model.GetModelPart("FluidModelPart.GENERIC_not_moving").Nodes:
                node.SetSolutionStepValue(KratosMultiphysics.MESH_DISPLACEMENT_X,0, 0)
                node.Fix(KratosMultiphysics.MESH_DISPLACEMENT_X)
                node.SetSolutionStepValue(KratosMultiphysics.MESH_DISPLACEMENT_Y,0, 0)
                node.Fix(KratosMultiphysics.MESH_DISPLACEMENT_Y)

            self.MovePart('GENERIC_green', np.array([[1,0],[0,1/self.w]]), np.array([[0],[1.5]]), np.array([[0],[0]]))
            self.MovePart('GENERIC_yellow_up', np.array([[1,0],[0, (2/(3-self.w))]]), np.array([[0],[3]]), np.array([[0],[0]]))
            self.MovePart('GENERIC_yellow_down', np.array([[1,0],[0, 2/(3-self.w)]]), np.array([[0],[0]]), np.array([[0],[0]]))
            self.MovePart('GENERIC_blue', np.array([[1,0],[(self.w-1)/2, 1]]), np.array([[0],[0]]), np.array([[0],[(self.w-1)/4]]))
            self.MovePart('GENERIC_grey', np.array([[1,0],[(1-self.w)/2, 1]]), np.array([[0],[0]]), np.array([[0],[- (self.w-1)/4]]))

            #############################
            ###  MOVE EACH SUB-PART   ###
            #############################
            self.MovePart('GENERIC_green', np.array([[1,0],[0,1/self.w]]), np.array([[0],[1.5]]), np.array([[0],[0]]))
            self.MovePart('GENERIC_yellow_up', np.array([[1,0],[0, (2/(3-self.w))]]), np.array([[0],[3]]), np.array([[0],[0]]))
            self.MovePart('GENERIC_yellow_down', np.array([[1,0],[0, 2/(3-self.w)]]), np.array([[0],[0]]), np.array([[0],[0]]))
            self.MovePart('GENERIC_blue', np.array([[1,0],[(self.w-1)/2, 1]]), np.array([[0],[0]]), np.array([[0],[(self.w-1)/4]]))
            self.MovePart('GENERIC_grey', np.array([[1,0],[(1-self.w)/2, 1]]), np.array([[0],[0]]), np.array([[0],[- (self.w-1)/4]]))

    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()
        self.StoreBifurcationData()

    def GetBifuracationData(self):
        return self.velocity_y_at_control_point ,  self.narrowing_width


if __name__ == "__main__":

    with open("ProjectParameters.json", 'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())

    analysis_stage_module_name = parameters["analysis_stage"].GetString()
    analysis_stage_class_name = analysis_stage_module_name.split('.')[-1]
    analysis_stage_class_name = ''.join(x.title() for x in analysis_stage_class_name.split('_'))

    # analysis_stage_module = importlib.import_module(analysis_stage_module_name)
    # analysis_stage_class = getattr(analysis_stage_module, analysis_stage_class_name)

    analysis_stage_class = FluidDynamicsAnalysisWithFlush

    global_model = KratosMultiphysics.Model()
    # simulation = CreateRomAnalysisInstance(analysis_stage_class, global_model, parameters)
    
    simulation = analysis_stage_class(global_model, parameters)
    simulation.Run()

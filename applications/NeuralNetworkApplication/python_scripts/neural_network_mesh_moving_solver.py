from importlib import import_module
from xml.etree.ElementTree import ProcessingInstruction
from KratosMultiphysics import python_solver
import KratosMultiphysics.NeuralNetworkApplication.input_dataclasses as InputDataclasses
import KratosMultiphysics.NeuralNetworkApplication.data_loading_utilities
import KratosMultiphysics
from KratosMultiphysics.NeuralNetworkApplication.input_dataclasses import ListDataWithLookback, ListNeuralNetworkData, NeuralNetworkData
from KratosMultiphysics.NeuralNetworkApplication.neural_network_solver import NeuralNetworkSolver
from KratosMultiphysics.NeuralNetworkApplication.neural_network_process_factory import NeuralNetworkProcessFactory
from importlib import import_module
import KratosMultiphysics.MeshMovingApplication.python_solvers_wrapper_mesh_motion as mesh_mothion_solvers_wrapper

import tensorflow.keras as keras
import numpy as np

def CreateSolver(project_parameters, model = None):
    return NeuralNetworkMeshMovingSolver(project_parameters,model)

class NeuralNetworkMeshMovingSolver(NeuralNetworkSolver):
    """
    This class implements a simple moving mesh variation of the neural network solver.
    """

    # TODO: Addapt solver to PythonSolver structure (project parameters)
    def __init__(self, project_parameters, model = None):      
        
        self.settings = project_parameters
        super().__init__(project_parameters,model)
        self.mesh_motion_solver_settings = project_parameters["mesh_motion_solver_settings"]
        parallelism = project_parameters["parallel_type"].GetString()
        self.parallelism = parallelism
         # Creating the mesh-motion solver
        if not self.mesh_motion_solver_settings.Has("echo_level"):
            self.mesh_motion_solver_settings.AddValue("echo_level", self.settings["echo_level"])

        # Making sure the settings are consistent btw fluid and mesh-motion
        if self.mesh_motion_solver_settings.Has("model_part_name"):
            if not self.model_geometry_name == self.mesh_motion_solver_settings["model_part_name"].GetString():
                err_msg =  'NeuralNetwork- and Mesh-Solver have to use the same MainModelPart ("model_part_name")!\n'
                err_msg += 'Use "mesh_motion_parts" for specifying mesh-motion on sub-model-parts'
                raise Exception(err_msg)
        else:
            self.mesh_motion_solver_settings.AddValue("model_part_name", project_parameters["model_part_name"])

        domain_size = project_parameters["domain_size"].GetInt()
        if self.mesh_motion_solver_settings.Has("domain_size"):
            mesh_motion_domain_size = self.mesh_motion_solver_settings["domain_size"].GetInt()
            if not domain_size == mesh_motion_domain_size:
                raise Exception('NeuralNetwork- and Mesh-Solver have to use the same "domain_size"!')
        else:
            self.mesh_motion_solver_settings.AddValue("domain_size", project_parameters["domain_size"])
        
        self.mesh_motion_solver_full_mesh = mesh_mothion_solvers_wrapper.CreateSolverByParameters(
            model, self.mesh_motion_solver_settings, parallelism)

    def Initialize(self):
        # Saving the ALE-interface-parts for later
        # this can only be done AFTER reading the ModelPart
        main_model_part_name = self.settings["model_part_name"].GetString()

        mesh_boundary_parts_params = self.settings["mesh_boundary_parts"]
        self.mesh_boundary_parts = []
        for i_name in range(mesh_boundary_parts_params.size()):
            sub_model_part_name = mesh_boundary_parts_params[i_name].GetString()
            full_model_part_name = main_model_part_name + "." + sub_model_part_name
            self.mesh_boundary_parts.append(self.model[full_model_part_name])

        mesh_motion_parts_params = self.settings["mesh_motion_parts"]
        self.mesh_motion_solvers = []
        if mesh_motion_parts_params.size() == 0:
            # the entire Main-ModelPart is used in the Mesh-Solver
            self.mesh_motion_solvers.append(self.mesh_motion_solver_full_mesh)
        else:
            # SubModelParts of the Main-ModelPart are used in the Mesh-Solver
            # each SubModelPart has its own mesh-solver
            # Note that these solvers do NOT need to call AddVariables and AddDofs
            # since this is done already for the MainModelPart
            for i_name in range(mesh_motion_parts_params.size()):
                sub_model_part_name = mesh_motion_parts_params[i_name].GetString()
                if sub_model_part_name == main_model_part_name:
                    err_msg =  'The MainModelPart cannot be used as one of the Sub-Mesh-Solvers!\n'
                    err_msg += 'Remove "mesh_motion_parts" for specifying mesh-motion on the MainModelPart'
                    raise Exception(err_msg)
                full_model_part_name = main_model_part_name + "." + sub_model_part_name
                sub_mesh_solver_settings = self.settings["mesh_motion_solver_settings"].Clone()
                sub_mesh_solver_settings["model_part_name"].SetString(full_model_part_name)

                self.mesh_motion_solvers.append(mesh_mothion_solvers_wrapper.CreateSolverByParameters(
                    self.model, sub_mesh_solver_settings, self.parallelism))

        for mesh_solver in self.mesh_motion_solvers:
            mesh_solver.Initialize()

        super().Initialize()
        
    def InitializeSolutionStep(self):
        for mesh_solver in self.mesh_motion_solvers:
            mesh_solver.InitializeSolutionStep()
        super().InitializeSolutionStep()

    def Predict(self):
        for mesh_solver in self.mesh_motion_solvers:
            mesh_solver.Predict()
        super().Predict()

    def FinalizeSolutionStep(self):
        for mesh_solver in self.mesh_motion_solvers:
            mesh_solver.FinalizeSolutionStep()
        super().FinalizeSolutionStep()

    def SolveSolutionStep(self):
        is_converged = True
        for mesh_solver in self.mesh_motion_solvers:
            is_converged &= mesh_solver.SolveSolutionStep()
            
        super().SolveSolutionStep()

    def Finalize(self):
        for mesh_solver in self.mesh_motion_solvers:
            mesh_solver.Finalize()
        
        super().Finalize()

    def GetMeshMotionSolver(self):
        if len(self.mesh_motion_solvers) > 1:
            raise Exception('More than one mesh-motion-solver exists, please use "GetMeshMotionSolvers"')
        return self.mesh_motion_solvers[0]

    def GetMeshMotionSolvers(self):
        return self.mesh_motion_solvers

    def MoveMesh(self):
        for mesh_solver in self.mesh_motion_solvers:
            mesh_solver.MoveMesh()

    def Check(self):
        for mesh_solver in self.mesh_motion_solvers:
            mesh_solver.Check()
        super().Check()

    def Clear(self):
        for mesh_solver in self.mesh_motion_solvers:
            mesh_solver.Clear()
        super().Clear()

    def AddVariables(self):
        self.mesh_motion_solver_full_mesh.AddVariables()

    def AddDofs(self):
        self.mesh_motion_solver_full_mesh.AddDofs()


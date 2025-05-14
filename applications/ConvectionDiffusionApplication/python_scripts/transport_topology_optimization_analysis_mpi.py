from sys import argv

import numpy as np #import the numpy library
import scipy as sp #import the scipy library
from scipy.spatial import KDTree
from scipy.sparse import dok_matrix, lil_matrix
import mmapy as MMA # import the MMA subroutines python library: https://github.com/arjendeetman/GCMMA-MMA-Python
import time

import KratosMultiphysics as KratosMultiphysics
# Import applications
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD
import KratosMultiphysics.ConvectionDiffusionApplication as KratosCD
import KratosMultiphysics.MeshingApplication as KratosMMG

from KratosMultiphysics.FluidDynamicsApplication.fluid_dynamics_analysis import FluidDynamicsAnalysis
from KratosMultiphysics.ConvectionDiffusionApplication.convection_diffusion_analysis import ConvectionDiffusionAnalysis
from KratosMultiphysics.FluidDynamicsApplication.fluid_topology_optimization_analysis import FluidTopologyOptimizationAnalysis
from KratosMultiphysics.ConvectionDiffusionApplication.transport_topology_optimization_analysis import TransportTopologyOptimizationAnalysis
from KratosMultiphysics.ConvectionDiffusionApplication import transport_topology_optimization_solver_mpi

class TransportTopologyOptimizationAnalysisMpi(TransportTopologyOptimizationAnalysis):

    def _CreateTopologyOptimizationSolvers(self):
        """
        This method creates the NS and ADJ_NS solvers
        """
        self.physics_solver = transport_topology_optimization_solver_mpi.CreateSolver(self.model, self.project_parameters)
        self.adjoint_solver = transport_topology_optimization_solver_mpi.CreateSolver(self.model, self.project_parameters, isAdjointSolver=True)

    def _CreateSolver(self, isAdjointSolver = False):
        """
        This method creates a solver
        isAdjointSolver == False --> physics_solver
        isAdjointSolver == True  --> adjoint_solver
        """
        return transport_topology_optimization_solver_mpi.CreateSolver(self.model, self.project_parameters, isAdjointSolver)
    
    def PrepareAdjointSolver(self):
        # pass
        """This method prepares the Adjoint Navier-Stokes Solver in the AnalysisStage 
        Usage: It is designed to be called ONCE, BEFORE the execution of the solution-loop
        Prepare Solver : ImportModelPart -> PrepareModelPart -> AddDofs
        """
        # Modelers:
        self._GetAdjointSolver().ImportModelPart(model_parts=self._GetPhysicsMainModelPartsList(), physics_solver_distributed_model_part_importer=self.physics_solver.distributed_model_part_importer)
        self._GetAdjointSolver().PrepareModelPart()
        self._GetAdjointSolver().AddDofs()

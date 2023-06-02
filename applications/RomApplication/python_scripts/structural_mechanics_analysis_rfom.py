import numpy
import KratosMultiphysics
import KratosMultiphysics.RomApplication as romapp
import KratosMultiphysics.StructuralMechanicsApplication
from KratosMultiphysics.RomApplication.empirical_cubature_method import EmpiricalCubatureMethod
from KratosMultiphysics.RomApplication import python_solvers_wrapper_rom as solver_wrapper
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import StructuralMechanicsAnalysis
from KratosMultiphysics.RomApplication.randomized_singular_value_decomposition import RandomizedSingularValueDecomposition

import json
import numpy as np

class StructuralMechanicsAnalysisRFom(StructuralMechanicsAnalysis):

    def __init__(self,model,rom_project_parameters,fom_project_parameters, hyper_reduction_element_selector = None):
        KratosMultiphysics.Logger.PrintWarning('\x1b[1;31m[DEPRECATED CLASS] \x1b[0m',"\'StructuralMechanicsAnalysisRFom\'", "class is deprecated. Use the \'RomAnalysis\' one instead.")
        self.svd_truncation_tolerance = rom_project_parameters["svd_truncation_tolerance"].GetDouble()
        rom_project_parameters.RemoveValue("svd_truncation_tolerance")
        self.rom_project_params = KratosMultiphysics.Parameters()
        self.rom_project_params.AddValue("solver_settings", fom_project_parameters["solver_settings"])
        self.rom_project_params["solver_settings"].AddValue("rom_settings", rom_project_parameters)
        super().__init__(model,fom_project_parameters)
        self._is_rom = False

    #### Internal functions ####
    def _CreateSolver(self):
        self._CreateRomSolver()
        self._CreateFomSolver()
        return self.fom_solver

    #### Internal functions ####
    def _CreateRomSolver(self):
        """ Create the Solver (and create and import the ModelPart if it is not alread in the model) """
        ## Solver construction
        self.rom_solver = solver_wrapper.CreateSolverByParameters(self.model, self.rom_project_params["solver_settings"],self.project_parameters["problem_data"]["parallel_type"].GetString())
        self.n_nodal_unknowns = len(self.rom_project_params["solver_settings"]["rom_settings"]["nodal_unknowns"].GetStringArray())

    def _CreateFomSolver(self):
        self.fom_solver = super()._CreateSolver()

    def _GetSimulationName(self):
        if self._is_rom:
            return "::[ROM Simulation]:: "
        else:
            return "::[FOM Simulation]" + super()._GetSimulationName()

    def _GetComputingModelPart(self):
        if self._is_rom:
            return self.rom_solver.GetComputingModelPart()
        else:
            return self.fom_solver.GetComputingModelPart()

    def SetSolverToRom(self):
        self._is_rom = True
        self._solver = self.rom_solver

    def SetSolverToFom(self):
        self._is_rom = False
        self._solver = self.fom_solver

    def UpdateRomBases(self,snapshots_matrix):
        # Calculate the randomized SVD of the snapshots matrix
        u,_,_,eSVD = RandomizedSingularValueDecomposition().Calculate(snapshots_matrix, self.svd_truncation_tolerance)
        num_rom_dofs = numpy.shape(u)[1]

        node_counter = 0
        for node in self._GetComputingModelPart().Nodes:
            aux = KratosMultiphysics.Matrix(self.n_nodal_unknowns, num_rom_dofs)
            for j in range(self.n_nodal_unknowns):
                for i in range(num_rom_dofs):
                    aux[j,i] = u.item((node_counter * self.n_nodal_unknowns +j,i))
            node.SetValue(romapp.ROM_BASIS, aux ) # ROM basis
            node_counter += 1

        self._GetComputingModelPart().ProcessInfo.SetValue(romapp.NUM_ROM_BASIS,num_rom_dofs)

    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()

    def Finalize(self):
        super().Finalize()
        if self._is_rom:
            self.fom_solver.Finalize()
        else:
            self.rom_solver.Finalize()

from importlib import import_module
import numpy
from numpy import linalg as LA
import json

import KratosMultiphysics as Kratos
from KratosMultiphysics.analysis_stage import AnalysisStage
from KratosMultiphysics.OptimizationApplication.execution_policies.execution_policy import ExecutionPolicy
from KratosMultiphysics.OptimizationApplication.utilities.helper_utilities import GetClassModuleFromKratos
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
import KratosMultiphysics.RomApplication as ROM
from KratosMultiphysics.RomApplication.randomized_singular_value_decomposition import RandomizedSingularValueDecomposition

def Factory(model: Kratos.Model, parameters: Kratos.Parameters, _) -> ExecutionPolicy:
    if not parameters.Has("name"):
        raise RuntimeError(f"SteppingRomAnalysisExecutionPolicy instantiation requires a \"name\" in parameters [ parameters = {parameters}].")
    if not parameters.Has("settings"):
        raise RuntimeError(f"SteppingRomAnalysisExecutionPolicy instantiation requires a \"settings\" in parameters [ parameters = {parameters}].")

    return SteppingRomAnalysisExecutionPolicy(parameters["name"].GetString(), model, parameters["settings"])

class SteppingRomAnalysisExecutionPolicy(ExecutionPolicy):
    def __init__(self, name: str, model: Kratos.Model, parameters: Kratos.Parameters):
        super().__init__(name)

        default_settings = Kratos.Parameters("""{
            "model_part_names" : [],
            "analysis_module"  : "KratosMultiphysics",
            "analysis_type"    : "",
            "analysis_settings": {}
        }""")
        self.model = model
        self.parameters = parameters
        self.parameters.ValidateAndAssignDefaults(default_settings)

        analysis_module = parameters["analysis_module"].GetString()
        rom_analysis_type = parameters["analysis_type"].GetString() + "RFom"
        analysis_settings = parameters["analysis_settings"]

        if analysis_module == "KratosMultiphysics":
            rom_analysis_module = GetClassModuleFromKratos(rom_analysis_type)

        self.model_parts = []
        rom_analysis_full_module = f"{rom_analysis_module}.{Kratos.StringUtilities.ConvertCamelCaseToSnakeCase(rom_analysis_type)}"

        self.rom_analysis: AnalysisStage = getattr(import_module(rom_analysis_full_module), rom_analysis_type)(self.model, analysis_settings.Clone())
        self.analysis = self.rom_analysis

        #TODO: get optimization iteration number from OptimizationProblem
        self.optimization_iteration = 0

        # Initialize the snapshots data list
        self.snapshots_data_list = []
        self.snapshot_variables_list = ["DISPLACEMENT_X","DISPLACEMENT_Y","DISPLACEMENT_Z","ROTATION_X","ROTATION_Y","ROTATION_Z"]
        self.n_nodal_unknowns = len(self.snapshot_variables_list)
        self.svd_truncation_tolerance = 1e-6
        self.update_rom_bases_frequency = 10
        self.run_rom = False

    def GetAnalysisModelPart(self):
        return self.rom_analysis._GetSolver().GetComputingModelPart()

    def Initialize(self):
        self.rom_analysis.Initialize()

        # initialize model parts
        self.model_parts = [self.model[model_part_name] for model_part_name in self.parameters["model_part_names"].GetStringArray()]

    def Check(self) -> None:
        pass

    def Finalize(self) -> None:
        pass

    def Execute(self):

        self.optimization_iteration += 1

        self.RunAnalysis()

        self.Update()

    def Update(self):

        if self.optimization_iteration < self.update_rom_bases_frequency:
            self.AddToSnapshotMatrix()
            self.run_rom = False
        elif self.optimization_iteration % self.update_rom_bases_frequency == 0.0:
            self.AddToSnapshotMatrix()
            self.rom_analysis.UpdateRomBases(self._GetSnapshotsMatrix())
            self.run_rom = True
        elif self.update_rom_bases_frequency - (self.optimization_iteration % self.update_rom_bases_frequency) == 1:
            self.run_rom = False
        else:
            self.run_rom = True

        if self.run_rom:
            self.rom_analysis.SetSolverToRom()
        else:
            self.rom_analysis.SetSolverToFom()

    def RunAnalysis(self):

        time_before_analysis = []
        step_before_analysis = []
        delta_time_before_analysis = []

        for model_part in self.model_parts:
            time_before_analysis.append(model_part.ProcessInfo[Kratos.TIME])
            step_before_analysis.append(model_part.ProcessInfo[Kratos.STEP])
            delta_time_before_analysis.append(model_part.ProcessInfo[Kratos.DELTA_TIME])

        # Reset step/time iterators such that they match the optimization iteration after calling CalculateValue (which internally calls CloneTimeStep)
        for index, model_part in enumerate(self.model_parts):
            model_part.ProcessInfo.SetValue(Kratos.STEP, step_before_analysis[index] - 1)
            model_part.ProcessInfo.SetValue(Kratos.TIME, time_before_analysis[index] - 1)
            model_part.ProcessInfo.SetValue(Kratos.DELTA_TIME, 0)

        self.rom_analysis.time = self.rom_analysis._GetSolver().AdvanceInTime(self.rom_analysis.time)
        self.rom_analysis.InitializeSolutionStep()
        self.rom_analysis._GetSolver().Predict()
        self.rom_analysis._GetSolver().SolveSolutionStep()

        self.rom_analysis.FinalizeSolutionStep()
        self.rom_analysis.OutputSolutionStep()

        # Clear results or modifications on model parts
        for index, model_part in enumerate(self.model_parts):
            model_part.ProcessInfo.SetValue(Kratos.STEP, step_before_analysis[index])
            model_part.ProcessInfo.SetValue(Kratos.TIME, time_before_analysis[index])
            model_part.ProcessInfo.SetValue(Kratos.DELTA_TIME, delta_time_before_analysis[index])


    def AddToSnapshotMatrix(self):
        # Save the data in the snapshots data list
        aux_data_array = []
        for node in self.GetAnalysisModelPart().Nodes:
            for snapshot_var in self.snapshot_variables_list:
                aux_data_array.append(node.GetSolutionStepValue(Kratos.KratosGlobals.GetVariable(snapshot_var)))
        if len(self.snapshots_data_list) < self.update_rom_bases_frequency:
            self.snapshots_data_list.append(aux_data_array)
        else:
            self.snapshots_data_list.pop(0)
            self.snapshots_data_list.append(aux_data_array)

    def _GetSnapshotsMatrix(self):
        self.n_nodes = self.GetAnalysisModelPart().NumberOfNodes()
        self.n_data_cols = len(self.snapshots_data_list)
        snapshots_matrix = numpy.empty((self.n_nodal_unknowns*self.n_nodes,self.n_data_cols))
        for i_col in range(self.n_data_cols):
            aux_col = numpy.array(self.snapshots_data_list[i_col])
            snapshots_matrix[:,i_col] = aux_col.transpose()
        return snapshots_matrix



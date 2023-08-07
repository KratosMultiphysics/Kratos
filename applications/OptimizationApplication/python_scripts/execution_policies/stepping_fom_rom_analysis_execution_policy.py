from importlib import import_module
import numpy
from numpy import linalg as LA
import json

import KratosMultiphysics as Kratos
from KratosMultiphysics.OptimizationApplication.execution_policies.execution_policy import ExecutionPolicy
from KratosMultiphysics.OptimizationApplication.utilities.helper_utilities import GetClassModuleFromKratos
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
import KratosMultiphysics.RomApplication as ROM
import KratosMultiphysics.RomApplication.fom_rom_analysis

def Factory(model: Kratos.Model, parameters: Kratos.Parameters, _) -> ExecutionPolicy:
    if not parameters.Has("name"):
        raise RuntimeError(f"SteppingFomRomAnalysisExecutionPolicy instantiation requires a \"name\" in parameters [ parameters = {parameters}].")
    if not parameters.Has("settings"):
        raise RuntimeError(f"SteppingFomRomAnalysisExecutionPolicy instantiation requires a \"settings\" in parameters [ parameters = {parameters}].")

    return SteppingFomRomAnalysisExecutionPolicy(parameters["name"].GetString(), model, parameters["settings"])

class SteppingFomRomAnalysisExecutionPolicy(ExecutionPolicy):
    def __init__(self, name: str, model: Kratos.Model, parameters: Kratos.Parameters):
        super().__init__(name)

        default_settings = Kratos.Parameters("""{
            "model_part_names" : [],
            "analysis_module"  : "KratosMultiphysics",
            "fom_analysis_type"    : "",
            "fom_analysis_settings": {},
            "rom_analysis_settings": {},
            "update_rom_bases_frequency": 10,
            "svd_truncation_tolerance": 1.0e-6
        }""")
        self.model = model
        self.parameters = parameters
        self.parameters.ValidateAndAssignDefaults(default_settings)

        self.rom_parameters = self.parameters["rom_analysis_settings"].Clone()

        self.svd_truncation_tolerance = self.parameters["svd_truncation_tolerance"].GetDouble()
        self.update_rom_bases_frequency = self.parameters["update_rom_bases_frequency"].GetInt()

        # Initialize the snapshots data list
        self.snapshots_data_list = []
        self.snapshot_variables_list = self.rom_parameters["rom_settings"]["nodal_unknowns"].GetStringArray()
        self.n_nodal_unknowns = len(self.snapshot_variables_list)

        analysis_module = self.parameters["analysis_module"].GetString()
        rom_analysis_type = self.parameters["fom_analysis_type"].GetString()
        fom_analysis_settings = self.parameters["fom_analysis_settings"]

        if analysis_module == "KratosMultiphysics":
            rom_analysis_module = GetClassModuleFromKratos(rom_analysis_type)

        rom_analysis_full_module = f"{rom_analysis_module}.{Kratos.StringUtilities.ConvertCamelCaseToSnakeCase(rom_analysis_type)}"

        analysis_stage_class = getattr(import_module(rom_analysis_full_module), rom_analysis_type)

        instance_factory = Kratos.RomApplication.fom_rom_analysis.CreateFomRomAnalysisInstance
        self.analysis = instance_factory(analysis_stage_class, self.model, fom_analysis_settings.Clone(), self.rom_parameters)

        #some initializations
        #TODO: get optimization iteration number from OptimizationProblem
        self.optimization_iteration = 0
        self.run_rom = False
        self.model_parts = []

    def GetAnalysisModelPart(self):
        return self.analysis._GetSolver().GetComputingModelPart()

    def Initialize(self):
        self.analysis.Initialize()

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
            self.analysis.UpdateRomBases(self._GetSnapshotsMatrix(),self.svd_truncation_tolerance)
            self.run_rom = True
        elif self.update_rom_bases_frequency - (self.optimization_iteration % self.update_rom_bases_frequency) == 1:
            self.run_rom = False
        else:
            self.run_rom = True

        if self.run_rom:
            self.analysis.SetSolverToRom()
        else:
            self.analysis.SetSolverToFom()

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

        self.analysis.time = self.analysis._GetSolver().AdvanceInTime(self.analysis.time)
        self.analysis.InitializeSolutionStep()
        self.analysis._GetSolver().Predict()
        self.analysis._GetSolver().SolveSolutionStep()

        self.analysis.FinalizeSolutionStep()
        self.analysis.OutputSolutionStep()

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



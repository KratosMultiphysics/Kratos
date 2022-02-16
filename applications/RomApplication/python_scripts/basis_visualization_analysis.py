import KratosMultiphysics
from KratosMultiphysics import gid_output_process
import KratosMultiphysics.RomApplication as RomApp
from KratosMultiphysics.analysis_stage import AnalysisStage
from KratosMultiphysics.python_solver import PythonSolver

class FauxSolver(PythonSolver):
    """
    Pretends to be a solver, but actually just reads the results from the RomParameters.json
    """

    def __init__(self, model, settings):
        super().__init__(model, settings)

        model_part_name = settings["model_part_name"].GetString()

        if self.model.HasModelPart(model_part_name):
            self.main_model_part = self.model.GetModelPart(model_part_name)
        else:
            self.main_model_part = self.model.CreateModelPart(model_part_name)

        settings_filename = self.settings["rom_parameters_file"].GetString()
        with open(settings_filename, "r") as f:
            rom_params = KratosMultiphysics.Parameters(f.read())

        self.rom_settings = rom_params

        self.number_of_modes = rom_params["rom_settings"]["number_of_rom_dofs"].GetInt()

        self.variable_names = [p.GetString() for p in rom_params["rom_settings"]["nodal_unknowns"]]
        self.variables = None

    def GetVariableNames(self):
        return self.variable_names

    def GetVariables(self):
        return self.variables

    def GetDefaultParameters(self):
        return KratosMultiphysics.Parameters("""{
            "echo_level" : 0,
            "model_part_name" : "main",
            "rom_parameters_file" : "RomParameters.json",
            "model_import_settings" : {}
        }""")

    def Check(self):
        pass

    def ImportModelPart(self):
        self._ImportModelPart(self.main_model_part,self.settings["model_import_settings"])

    def GetComputingModelPart(self):
        return self.main_model_part

    def AdvanceInTime(self, current_time):
        t = current_time + 1
        self.GetComputingModelPart().ProcessInfo[KratosMultiphysics.STEP] += 1
        self.GetComputingModelPart().ProcessInfo[KratosMultiphysics.TIME] = t
        return t

    def Initialize(self):
        self.variables = []
        for name in self.GetVariableNames():
            variable = KratosMultiphysics.KratosGlobals.GetVariable(name)
            self.variables.append(variable)

        for variable in self.variables:
            KratosMultiphysics.VariableUtils().SetNonHistoricalVariableToZero(variable, self.GetComputingModelPart().Nodes)

    def GetCurrentMode(self):
        return self.GetComputingModelPart().ProcessInfo[KratosMultiphysics.STEP]

    def GetNumberOfModes(self):
        return self.number_of_modes

    def SolveSolutionStep(self):
        mode = self.GetCurrentMode()
        for node in self.GetComputingModelPart().Nodes:
            variable_idx = 0
            for variable in self.GetVariables():
                value = self.rom_settings["nodal_modes"][str(node.Id)][variable_idx][mode].GetDouble()
                node.SetValue(variable, value)
                variable_idx += 1
        return 0


class BasisVisualizationAnalysis(AnalysisStage):

    def _CreateSolver(self):
        return FauxSolver(self.model, self.project_parameters["solver_settings"])

    def KeepAdvancingSolutionLoop(self):
        return self._GetSolver().GetCurrentMode() + 1 != self._GetSolver().GetNumberOfModes()

    def __init__(self, model, project_parameters):

        if "problem_data" not in project_parameters:
            # Adding dummy values used by base analysis stage
            project_parameters.AddEmptyValue("problem_data")

            project_parameters["problem_data"] = KratosMultiphysics.Parameters("""
                {
                    "parallel_type": "OpenMP",
                    "start_time": 0.0,
                    "time_step": 1.0,
                    "end_time": 1.0,
                    "echo_level": 0
                }
            """)

        super().__init__(model, project_parameters)

        gid_settings = KratosMultiphysics.Parameters("""
        {
            "python_module" : "gid_output_process",
            "kratos_module" : "KratosMultiphysics",
            "process_name"  : "GiDOutputProcess",
            "Parameters" : {
                "model_part_name"        : "",
                "output_name"            : "",
                "postprocess_parameters" : {
                    "result_file_configuration" : {
                        "gidpost_flags"               : {
                            "GiDPostMode"           : "GiD_PostBinary",
                            "WriteDeformedMeshFlag" : "WriteDeformed",
                            "WriteConditionsFlag"   : "WriteConditions",
                            "MultiFileFlag"         : "SingleFile"
                        },
                        "file_label"                  : "step",
                        "output_control_type"         : "step",
                        "output_interval"             : 0,
                        "body_output"                 : true,
                        "node_output"                 : false,
                        "skin_output"                 : false,
                        "plane_output"                : [],
                        "nodal_results"               : [],
                        "nodal_nonhistorical_results" : [],
                        "gauss_point_results"         : []
                    },
                    "point_data_configuration"  : []
                }
            }
        }""")

        gid_settings["Parameters"]["model_part_name"].SetString(self.project_parameters["solver_settings"]["model_part_name"].GetString())
        gid_settings["Parameters"]["output_name"].SetString(self.project_parameters["solver_settings"]["model_import_settings"]["input_filename"].GetString())

        for varname in self._GetSolver().GetVariableNames():
            gid_settings["Parameters"]["postprocess_parameters"]["result_file_configuration"]["nodal_nonhistorical_results"].Append(varname)


        if not "output_processes" in self.project_parameters:
            self.project_parameters.AddEmptyValue("output_processes")

        if not "gid_output" in self.project_parameters["output_processes"]:
            self.project_parameters["output_processes"].AddEmptyArray("gid_output")

        self.project_parameters["output_processes"]["gid_output"].Append(gid_settings)


if __name__ == "__main__":

    with open("ProjectParameters.json", 'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())

    Model = KratosMultiphysics.Model()
    simulation = BasisVisualizationAnalysis(Model, parameters)
    simulation.Run()

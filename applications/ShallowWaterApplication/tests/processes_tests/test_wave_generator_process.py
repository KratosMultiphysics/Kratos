import KratosMultiphysics as KM
import KratosMultiphysics.ShallowWaterApplication as SW
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.process_factory import KratosProcessFactory
import os

try:
    import scipy
    missing_scipy = False
except ImportError:
    missing_scipy = True

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)


class TestWaveGeneratorProcess(KratosUnittest.TestCase):

    @KratosUnittest.skipIf(missing_scipy, "Missing python libraries (scipy)")
    def test_wave_generator_process_by_direction(self):
        model = KM.Model()
        settings = KM.Parameters("""{
            "process_list" : [{
                "python_module"  : "wave_generator_process",
                "kratos_module"  : "KratosMultiphysics.ShallowWaterApplication",
                "Parameters"            : {
                    "model_part_name"          : "model_part.Condition1",
                    "interval"                 : [0.0, "End"],
                    "direction"                : [1.0, 0.0, 0.0],
                    "normal_positive_outwards" : true,
                    "smooth_time"              : 0.0,
                    "wave_specifications"      : {
                        "wave_theory"               : "boussinesq",
                        "period"                    : 0.0,
                        "wavelength"                : 20.0,
                        "amplitude"                 : 0.5,
                        "depth"                     : 2.0,
                        "get_depth_from_model_part" : false
                    }
                }
            },{
                "kratos_module"   : "KratosMultiphysics",
                "python_module"   : "from_json_check_result_process",
                "Parameters"      : {
                    "model_part_name"  : "model_part",
                    "check_variables"  : ["HEIGHT","VELOCITY"],
                    "input_file_name"  : "wave_generator_process_reference.json",
                    "time_frequency"   : 10.0,
                    "tolerance"        : 1e-12
                }
                // "kratos_module"   : "KratosMultiphysics",
                // "python_module"   : "json_output_process",
                // "Parameters"      : {
                //     "model_part_name"  : "model_part",
                //     "output_variables" : ["HEIGHT","VELOCITY"],
                //     "output_file_name" : "wave_generator_process_reference.json",
                //     "time_frequency"   : 10.0
                // }
            }]
        }""")
        reference_file_name = GetFilePath("wave_generator_process_reference.json")
        settings["process_list"][1]["Parameters"]["input_file_name"].SetString(reference_file_name)

        end_time = 10
        CreateAndSetModelPartForWaveGeneratorProcess(model)
        SolutionLoopForProcesses(model, settings["process_list"], end_time)


    @KratosUnittest.skipIf(missing_scipy, "Missing python libraries (scipy)")
    def test_wave_generator_process_by_normal(self):
        model = KM.Model()
        settings = KM.Parameters("""{
            "process_list" : [{
                "python_module"  : "wave_generator_process",
                "kratos_module"  : "KratosMultiphysics.ShallowWaterApplication",
                "Parameters"     : {
                    "model_part_name"          : "model_part.Condition1",
                    "interval"                 : [0.0, "End"],
                    "direction"                : "normal",
                    "normal_positive_outwards" : true,
                    "smooth_time"              : 0.0,
                    "wave_specifications"      : {
                        "wave_theory"               : "boussinesq",
                        "period"                    : 0.0,
                        "wavelength"                : 20.0,
                        "amplitude"                 : 0.5,
                        "depth"                     : 2.0,
                        "get_depth_from_model_part" : false
                    }
                }
            },{
                "kratos_module"   : "KratosMultiphysics",
                "python_module"   : "from_json_check_result_process",
                "Parameters"      : {
                    "model_part_name"  : "model_part",
                    "check_variables"  : ["HEIGHT","VELOCITY"],
                    "input_file_name"  : "wave_generator_process_reference.json",
                    "time_frequency"   : 10.0,
                    "tolerance"        : 1e-12
                }
            }]
        }""")
        reference_file_name = GetFilePath("wave_generator_process_reference.json")
        settings["process_list"][1]["Parameters"]["input_file_name"].SetString(reference_file_name)

        end_time = 10
        CreateAndSetModelPartForWaveGeneratorProcess(model)
        SolutionLoopForProcesses(model, settings["process_list"], end_time)


    @KratosUnittest.skipIf(missing_scipy, "Missing python libraries (scipy)")
    def test_wave_generator_process_with_topography(self):
        model = KM.Model()
        settings = KM.Parameters("""{
            "process_list" : [{
                "python_module"  : "process_factory",
                "kratos_module"  : "KratosMultiphysics",
                "process_name"   : "ApplyConstantScalarValueProcess",
                "Parameters"     : {
                    "model_part_name"      : "model_part",
                    "variable_name"        : "TOPOGRAPHY",
                    "value"                : -2.0
                }
            },{
                "python_module"  : "wave_generator_process",
                "kratos_module"  : "KratosMultiphysics.ShallowWaterApplication",
                "Parameters"     : {
                    "model_part_name"          : "model_part.Condition1",
                    "interval"                 : [0.0, "End"],
                    "direction"                : "normal",
                    "normal_positive_outwards" : true,
                    "smooth_time"              : 0.0,
                    "wave_specifications"      : {
                        "wave_theory"               : "boussinesq",
                        "period"                    : 0.0,
                        "wavelength"                : 20.0,
                        "amplitude"                 : 0.5,
                        "depth"                     : 2.0,
                        "get_depth_from_model_part" : true
                    }
                }
            },{
                "kratos_module"   : "KratosMultiphysics",
                "python_module"   : "from_json_check_result_process",
                "Parameters"      : {
                    "model_part_name"  : "model_part",
                    "check_variables"  : ["HEIGHT","VELOCITY"],
                    "input_file_name"  : "wave_generator_process_reference.json",
                    "time_frequency"   : 10.0,
                    "tolerance"        : 1e-12
                }
            }]
        }""")
        reference_file_name = GetFilePath("wave_generator_process_reference.json")
        settings["process_list"][2]["Parameters"]["input_file_name"].SetString(reference_file_name)

        end_time = 10
        CreateAndSetModelPartForWaveGeneratorProcess(model)
        SolutionLoopForProcesses(model, settings["process_list"], end_time)


def CreateAndSetModelPartForWaveGeneratorProcess(model):
    model_part = model.CreateModelPart("model_part")

    model_part.AddNodalSolutionStepVariable(KM.VELOCITY)
    model_part.AddNodalSolutionStepVariable(SW.HEIGHT)
    model_part.AddNodalSolutionStepVariable(SW.TOPOGRAPHY)

    model_part_io = KM.ModelPartIO(GetFilePath("model_part"))
    model_part_io.ReadModelPart(model_part)

    KM.VariableUtils().AddDof(KM.VELOCITY_X, model_part)
    KM.VariableUtils().AddDof(KM.VELOCITY_Y, model_part)
    KM.VariableUtils().AddDof(SW.HEIGHT, model_part)

    model_part.ProcessInfo[KM.GRAVITY_Z] = 9.81
    model_part.ProcessInfo[KM.DELTA_TIME] = 15
    model_part.ProcessInfo[KM.TIME] = 0.0


def SolutionLoopForProcesses(model, process_list, end_time):
    model_part = model.GetModelPart("model_part")
    element_replace_settings = KM.Parameters("""{
        "element_name"   : "WaveElement2D3N",
        "condition_name" : "WaveCondition2D2N"
    }""")
    KM.ReplaceElementsAndConditionsProcess(model_part, element_replace_settings).Execute()
    list_of_processes = KratosProcessFactory(model).ConstructListOfProcesses(process_list)

    for process in list_of_processes:
        process.ExecuteInitialize()

    for process in list_of_processes:
        process.Check()

    for process in list_of_processes:
        process.ExecuteBeforeSolutionLoop()

    while model_part.ProcessInfo[KM.TIME] < end_time:
        model_part.ProcessInfo[KM.TIME] += model_part.ProcessInfo[KM.DELTA_TIME]

        for process in list_of_processes:
            process.ExecuteInitializeSolutionStep()

        for process in list_of_processes:
            process.ExecuteBeforeOutputStep()

        for process in list_of_processes:
            process.ExecuteAfterOutputStep()

        for process in list_of_processes:
            process.ExecuteFinalizeSolutionStep()

    for process in list_of_processes:
        process.ExecuteFinalize()


if __name__ == '__main__':
    KM.Logger.GetDefaultOutput().SetSeverity(KM.Logger.Severity.WARNING)
    KratosUnittest.main()

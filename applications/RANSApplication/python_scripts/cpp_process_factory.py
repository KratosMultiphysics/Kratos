import KratosMultiphysics as Kratos
from KratosMultiphysics import IsDistributedRun
from KratosMultiphysics.kratos_utilities import CheckIfApplicationsAvailable
import KratosMultiphysics.RANSApplication as KratosRANS

if (IsDistributedRun() and CheckIfApplicationsAvailable("TrilinosApplication")):
    from KratosMultiphysics.RANSApplication.TrilinosExtension import TrilinosRansWallDistanceCalculationProcess as wall_distance_calculation_process
elif (not IsDistributedRun()):
    from KratosMultiphysics.RANSApplication import RansWallDistanceCalculationProcess as wall_distance_calculation_process
else:
    raise Exception("Distributed run requires TrilinosApplication")


def Factory(settings, Model):
    if (not isinstance(settings, Kratos.Parameters)):
        raise Exception(
            "expected input shall be a Parameters object, encapsulating a json string"
        )
    if (not isinstance(Model, Kratos.Model)):
        raise Exception("expected input shall be a Model object")

    allowed_processes_list = [
        [
            "KTurbulentIntensityInletProcess",
            KratosRANS.RansKTurbulentIntensityInletProcess
        ],
        [
            "EpsilonTurbulentMixingLengthInletProcess",
            KratosRANS.RansEpsilonTurbulentMixingLengthInletProcess
        ],
        [
            "OmegaTurbulentMixingLengthInletProcess",
            KratosRANS.RansOmegaTurbulentMixingLengthInletProcess
        ],
        [
            "ApplyExactNodalPeriodicConditionProcess",
            KratosRANS.RansApplyExactNodalPeriodicConditionProcess
        ],
        [
            "ApplyFlagProcess",
            KratosRANS.RansApplyFlagToSkinProcess
        ],
        [
            "ClipScalarVariableProcess",
            KratosRANS.RansClipScalarVariableProcess
        ],
        [
            "LineOutputProcess",
            KratosRANS.RansLineOutputProcess
        ],
        [
            "WallDistanceCalculationProcess",
            wall_distance_calculation_process
        ],
        [
            "NutKEpsilonUpdateProcess",
            KratosRANS.RansNutKEpsilonUpdateProcess
        ],
        [
            "NutKOmegaUpdateProcess",
            KratosRANS.RansNutKOmegaUpdateProcess
        ],
        [
            "NutKOmegaSSTUpdateProcess",
            KratosRANS.RansNutKOmegaSSTUpdateProcess
        ],
        [
            "NutYPlusWallFunctionUpdateProcess",
            KratosRANS.RansNutYPlusWallFunctionUpdateProcess
        ],
        [
            "WallFunctionUpdateProcess",
            KratosRANS.RansWallFunctionUpdateProcess
        ],
        [
            "ComputeReactionsProcess",
            KratosRANS.RansComputeReactionsProcess
        ],
        [
            "CheckScalarBoundsProcess",
            RansCheckScalarBoundsProcess
        ]
    ]

    process_name = settings["process_name"].GetString()

    process_names_list = [
        allowed_processes_list[i][0]
        for i in range(len(allowed_processes_list))
    ]
    process_list = [
        allowed_processes_list[i][1]
        for i in range(len(allowed_processes_list))
    ]

    if (process_name not in process_names_list):
        msg = "Unknown process_name=\"" + process_name + \
            "\". \nFollowing process names are allowed:\n    "
        msg += "\n    ".join(sorted(process_names_list))
        raise Exception(msg + "\n")

    current_process = process_list[process_names_list.index(process_name)](
        Model, settings["Parameters"])

    Kratos.Logger.PrintInfo("RANSApplyCustomProcess",
                            "Created " + process_name + " with following properties...\n" + str(settings["Parameters"]))

    return current_process


class RansCheckScalarBoundsProcess(KratosRANS.RansFormulationProcess):
    """
    Checks bounds of a scalar variable for given model part

    Args:
        model (Kratos.Model): Kratos model
        settings (Kratos.Parameters): Settings for process
    """
    def __init__(self, model, settings):
        super().__init__()

        default_parameters = Kratos.Parameters("""
        {
            "model_part_name" : "PLEASE_SPECIFY_MODEL_PART_NAME",
            "variable_name"   : "PLEASE_SPECIFY_SCALAR_VARIABLE"
        }""")

        settings.ValidateAndAssignDefaults(default_parameters)

        self.variable = Kratos.KratosGlobals.GetVariable(settings["variable_name"].GetString())
        self.model_part = model[settings["model_part_name"].GetString()]

    def Execute(self):
        min_value = KratosRANS.RansVariableUtilities.GetMinimumScalarValue(self.model_part, self.variable)
        max_value = KratosRANS.RansVariableUtilities.GetMaximumScalarValue(self.model_part, self.variable)

        Kratos.Logger.PrintInfo(self.__class__.__name__, "{:s} is bounded between [ {:f}, {:f} ] in {:s}.".format(
            self.variable.Name(), min_value, max_value, self.model_part.Name))

    def ExecuteAfterCouplingSolveStep(self):
        self.Execute()

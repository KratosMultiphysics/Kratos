import KratosMultiphysics as Kratos
import KratosMultiphysics.RANSApplication as KratosRANS

from KratosMultiphysics import IsDistributedRun
from KratosMultiphysics.kratos_utilities import CheckIfApplicationsAvailable

if (IsDistributedRun()
        and CheckIfApplicationsAvailable("TrilinosApplication")):
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
        ["ApplyFlagProcess", KratosRANS.RansApplyFlagProcess],
        ["FindNodalNeighboursProcess", FindNodalNeighboursProcess],
        [
            "FindConditionParentProcess",
            KratosRANS.RansFindConditionParentProcess
        ],
        [
            "ScalarCellCenteredAveragingProcess",
            KratosRANS.RansScalarCellCenterAveragingProcess
        ],
        [
            "VectorCellCenteredAveragingProcess",
            KratosRANS.RansVectorCellCenterAveragingProcess
        ],
        [   "VectorAlignProcess",
            KratosRANS.RansVectorAlignProcess
        ],
        [
            "CalculateNormalsProcess",
            CalculateNormalsProcess
        ],
        [
            "WallDistanceCalculationProcess",
            wall_distance_calculation_process
        ],
        [
            "LogarithmicYPlusCalculationProcess",
            KratosRANS.RansLogarithmicYPlusCalculationProcess
        ],
        [
            "CheckScalarBoundsProcess",
            KratosRANS.RansCheckScalarBoundsProcess
        ],
        [
            "NutKEpsilonHighReCalculationProcess",
            KratosRANS.RansNutKEpsilonHighReCalculationProcess
        ],
        [
            "NutKEpsilonHighReUpdateProcess",
            KratosRANS.RansNutKEpsilonHighReUpdateProcess
        ],
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
            "ClipScalarVariableProcess",
            KratosRANS.RansClipScalarVariableProcess
        ],
        [
            "ApplyExactNodalPeriodicConditionProcess",
            KratosRANS.RansApplyExactNodalPeriodicConditionProcess
        ],
        [
            "NutYPlusWallFunctionProcess",
            KratosRANS.RansNutYPlusWallFunctionProcess
        ],
        [
            "NutYPlusWallFunctionUpdateProcess",
            KratosRANS.RansNutYPlusWallFunctionUpdateProcess
        ],
        [
            "NuTLowReCalculationProcess",
            KratosRANS.RansNutLowReCalculationProcess
        ],
        [
            "LogarithmicYPlusVelocitySensitivitiesProcess",
            KratosRANS.RansLogarithmicYPlusVelocitySensitivitiesProcess
        ],
        [
            "NutKEpsilonHighReSensitivitiesProcess",
            KratosRANS.RansNutKEpsilonHighReSensitivitiesProcess
        ],
        [
            "LineOutputProcess",
            KratosRANS.RansLineOutputProcess
        ],
        [
            "WallFunctionUpdateProcess",
            KratosRANS.RansWallFunctionUpdateProcess
        ],
        [
            "NutYPlusWallFunctionSensitivitiesProcess",
            KratosRANS.RansNutYPlusWallFunctionSensitivitiesProcess
        ],
        [
            "CheckVectorBoundsProcess",
            KratosRANS.RansCheckVectorBoundsProcess
        ],
        [
            "ComputeReactionsProcess",
            KratosRANS.RansComputeReactionsProcess
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
        msg = "Unknown process_name=\"" + process_name + "\". \nFollowing process names are allowed:\n    "
        msg += "\n    ".join(sorted(process_names_list))
        raise Exception(msg + "\n")

    current_process = process_list[process_names_list.index(process_name)](
        Model, settings["Parameters"])

    Kratos.Logger.PrintInfo("RANSApplyCustomProcess",
                            "Created " + process_name + " with following properties...\n" + str(settings["Parameters"]))

    return current_process


class FindNodalNeighboursProcess(Kratos.Process):
    def __init__(self, Model, settings):
        Kratos.Process.__init__(self)

        default_parameters = Kratos.Parameters("""
            {
                "model_part_name"          : "PLEASE_CHOOSE_MODEL_PART_NAME",
                "average_neighbour_nodes"  : 10,
                "average_neighbour_element": 10
            }  """)

        settings.ValidateAndAssignDefaults(default_parameters)
        self.settings = settings

        self.model_part = Model[settings["model_part_name"].GetString()]
        self.process = Kratos.FindNodalNeighboursProcess(
            self.model_part,
            self.settings["average_neighbour_element"].GetInt(),
            self.settings["average_neighbour_nodes"].GetInt())

    def Check(self):
        self.process.Check()

    def ExecuteInitialize(self):
        self.process.Execute()
        Kratos.Logger.PrintInfo(
            "FindNodalNeighboursProcess",
            "Nodal neighbours found for nodes in " + self.model_part.Name +
            ".")

    def Execute(self):
        self.process.Execute()
        Kratos.Logger.PrintInfo(
            "FindNodalNeighboursProcess",
            "Nodal neighbours found for nodes in " + self.model_part.Name +
            ".")


class CalculateNormalsProcess(Kratos.Process):
    def __init__(self, Model, settings):
        Kratos.Process.__init__(self)

        default_parameters = Kratos.Parameters("""
            {
                "model_part_name"          : "PLEASE_CHOOSE_MODEL_PART_NAME"
            }  """)

        settings.ValidateAndAssignDefaults(default_parameters)
        self.settings = settings

        self.model_part = Model[settings["model_part_name"].GetString()]
        self.domain_size = self.model_part.ProcessInfo[Kratos.DOMAIN_SIZE]
        self.process = Kratos.NormalCalculationUtils()

    def ExecuteInitialize(self):
        self.process.CalculateOnSimplex(self.model_part, self.domain_size)
        Kratos.Logger.PrintInfo(
            "NormalCalculationUtils", "Nodal normals calculated for nodes in "
            + self.model_part.Name + ".")

    def Execute(self):
        self.process.CalculateOnSimplex(self.model_part, self.domain_size)
        Kratos.Logger.PrintInfo(
            "NormalCalculationUtils", "Nodal normals calculated for nodes in "
            + self.model_part.Name + ".")

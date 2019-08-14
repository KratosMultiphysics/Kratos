import KratosMultiphysics as Kratos
import KratosMultiphysics.RANSModellingApplication as KratosRANS


def Factory(settings, Model):
    if (type(settings) != Kratos.Parameters):
        raise Exception(
            "expected input shall be a Parameters object, encapsulating a json string"
        )
    if (type(Model) != Kratos.Model):
        raise Exception("expected input shall be a Model object")

    allowed_process_names_list = [
        "ApplyFlagProcess", "FindNodalNeighboursProcess",
        "FindConditionParentProcess",
        "ApplyScalarCellCenteredAveragingProcess",
        "ApplyVectorCellCenteredAveragingProcess", "ApplyVectorAlignProcess",
        "CalculateNormalsProcess", "WallDistanceCalculationProcess",
        "LogarithmicYPlusCalculationProcess", "CheckScalarBoundsProcess",
        "EpsilonWallFunctionProcess", "NuTKWallFunctionProcess",
        "NuTHighReCalculationProcess", "ApplyKTurbulentIntensityInletProcess",
        "ApplyEpsilonTurbulentMixingLengthInletProcess",
        "ApplyKWallFrictionVelocityProcess",
        "ApplyEpsilonWallFrictionVelocityProcess", "ClipScalarVariableProcess",
        "ClipScalarVariableByNeighbourAveragingProcess",
        "ApplyExactNodalPeriodicConditionProcess",
        "ApplyYPlusKCalculationProcess", "ApplyNutYPlusWallFunctionProcess",
        "NuTLowReCalculationProcess", "CheckScalarConditionBoundsProcess"
    ]

    process_name = settings["process_name"].GetString()
    if (process_name not in allowed_process_names_list):
        msg = "Unknown process_name=\"" + process_name + "\". Following process names are allowed:\n    "
        msg += "\n    ".join(allowed_process_names_list)
        raise Exception(msg + "\n")

    Kratos.Logger.PrintInfo("RANSApplyCustomProcess",
                            "Creating " + process_name)

    if (process_name == "ApplyFlagProcess"):
        return KratosRANS.RansApplyFlagProcess(Model, settings["Parameters"])
    elif (process_name == "FindNodalNeighboursProcess"):
        return FindNodalNeighboursProcess(Model, settings["Parameters"])
    elif (process_name == "FindConditionParentProcess"):
        return KratosRANS.RansFindConditionParentProcess(
            Model, settings["Parameters"])
    elif (process_name == "ApplyScalarCellCenteredAveragingProcess"):
        return KratosRANS.RansScalarCellCenterAveragingProcess(
            Model, settings["Parameters"])
    elif (process_name == "ApplyVectorCellCenteredAveragingProcess"):
        return KratosRANS.RansVectorCellCenterAveragingProcess(
            Model, settings["Parameters"])
    elif (process_name == "ApplyVectorAlignProcess"):
        return KratosRANS.RansVectorAlignProcess(Model, settings["Parameters"])
    elif (process_name == "CalculateNormalsProcess"):
        return CalculateNormalsProcess(Model, settings["Parameters"])
    elif (process_name == "WallDistanceCalculationProcess"):
        return KratosRANS.RansWallDistanceCalculationProcess(
            Model, settings["Parameters"])
    elif (process_name == "LogarithmicYPlusCalculationProcess"):
        return KratosRANS.RansLogarithmicYPlusCalculationProcess(
            Model, settings["Parameters"])
    elif (process_name == "CheckScalarBoundsProcess"):
        return KratosRANS.RansCheckScalarBoundsProcess(Model,
                                                       settings["Parameters"])
    elif (process_name == "EpsilonWallFunctionProcess"):
        return KratosRANS.RansEpsilonWallFunctionProcess(
            Model, settings["Parameters"])
    elif (process_name == "NuTKWallFunctionProcess"):
        return KratosRANS.RansNutKWallFunctionProcess(Model,
                                                      settings["Parameters"])
    elif (process_name == "NuTHighReCalculationProcess"):
        return KratosRANS.RansNutHighReCalculationProcess(
            Model, settings["Parameters"])
    elif (process_name == "ApplyKTurbulentIntensityInletProcess"):
        return KratosRANS.RansKTurbulentIntensityInletProcess(
            Model, settings["Parameters"])
    elif (process_name == "ApplyEpsilonTurbulentMixingLengthInletProcess"):
        return KratosRANS.RansEpsilonTurbulentMixingLengthInletProcess(
            Model, settings["Parameters"])
    elif (process_name == "ApplyKWallFrictionVelocityProcess"):
        return KratosRANS.RansKWallFrictionVelocityProcess(
            Model, settings["Parameters"])
    elif (process_name == "ApplyEpsilonWallFrictionVelocityProcess"):
        return KratosRANS.RansEpsilonWallFrictionVelocityProcess(
            Model, settings["Parameters"])
    elif (process_name == "ClipScalarVariableProcess"):
        return KratosRANS.RansClipScalarVariableProcess(
            Model, settings["Parameters"])
    elif (process_name == "ClipScalarVariableByNeighbourAveragingProcess"):
        return KratosRANS.RansClipScalarVariableByNeighbourAveragingProcess(
            Model, settings["Parameters"])
    elif (process_name == "ApplyExactNodalPeriodicConditionProcess"):
        return KratosRANS.RansApplyExactNodalPeriodicConditionProcess(
            Model, settings["Parameters"])
    elif (process_name == "ApplyYPlusKCalculationProcess"):
        return KratosRANS.RansYPlusKCalculationProcess(Model,
                                                       settings["Parameters"])
    elif (process_name == "ApplyNutYPlusWallFunctionProcess"):
        return KratosRANS.RansNutYPlusWallFunctionProcess(
            Model, settings["Parameters"])
    elif (process_name == "NuTLowReCalculationProcess"):
        return KratosRANS.RansNutLowReCalculationProcess(
            Model, settings["Parameters"])
    elif (process_name == "CheckScalarConditionBoundsProcess"):
        return KratosRANS.RansCheckScalarConditionBoundsProcess(
            Model, settings["Parameters"])


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

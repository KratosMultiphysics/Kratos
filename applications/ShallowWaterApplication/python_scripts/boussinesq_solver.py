# importing the Kratos Library
import KratosMultiphysics as KM
import KratosMultiphysics.ShallowWaterApplication as SW

## Import base class file
from KratosMultiphysics.ShallowWaterApplication.shallow_water_base_solver import ShallowWaterBaseSolver

def CreateSolver(model, custom_settings):
    return BoussinesqSolver(model, custom_settings)

class BoussinesqSolver(ShallowWaterBaseSolver):
    def __init__(self, model, settings):
        super().__init__(model, settings)
        self.element_name = "BoussinesqElement"
        self.condition_name = "BoussinesqCondition"
        if self.settings["time_integration_scheme"].GetString().lower() == "bdf":
            self.min_buffer_size = self.settings["time_integration_order"].GetInt() + 1
        else:
            self.min_buffer_size = 4
            if self.settings["time_integration_order"].GetInt() != 4:
                msg = "The order considered for the Adams-Moulton scheme is 4. The user provided order is {}"
                KM.Logger.PrintWarning(self.__class__.__name__, msg.format(self.settings["time_integration_order"].GetInt()))

    def AddDofs(self):
        KM.VariableUtils().AddDof(KM.VELOCITY_X, self.main_model_part)
        KM.VariableUtils().AddDof(KM.VELOCITY_Y, self.main_model_part)
        KM.VariableUtils().AddDof(SW.HEIGHT, self.main_model_part)
        KM.Logger.PrintInfo(self.__class__.__name__, "Boussinesq equations DOFs added correctly.")

    def AddVariables(self):
        super().AddVariables()
        self.main_model_part.AddNodalSolutionStepVariable(SW.DISPERSION_H)      # Intermediate field
        self.main_model_part.AddNodalSolutionStepVariable(SW.DISPERSION_V)      # Intermediate field
        self.main_model_part.AddNodalSolutionStepVariable(KM.RHS)               # This is used by the predictor in ABM scheme
        self.main_model_part.AddNodalSolutionStepVariable(KM.NODAL_AREA)        # This is used to assemble the RHS by the predictor and the dispersive fields

    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()
        SW.ShallowWaterUtilities().ComputeFreeSurfaceElevation(self.main_model_part)

    def _SetProcessInfo(self):
        super()._SetProcessInfo()
        self.main_model_part.ProcessInfo.SetValue(SW.RELATIVE_DRY_HEIGHT, self.settings["relative_dry_height"].GetDouble())
        self.main_model_part.ProcessInfo.SetValue(KM.STABILIZATION_FACTOR, self.settings["stabilization_factor"].GetDouble())
        self.main_model_part.ProcessInfo.SetValue(SW.SHOCK_STABILIZATION_FACTOR, self.settings["shock_capturing_factor"].GetDouble())

    def _CreateScheme(self):
        if self.settings["time_integration_scheme"].GetString().lower() == "bdf":
            scheme_settings = KM.Parameters()
            scheme_settings.AddStringArray("solution_variables", ["VELOCITY","HEIGHT"])
            scheme_settings.AddValue("integration_order", self.settings["time_integration_order"])
            scheme_settings.AddBool("project_dispersive_field", True)
            return SW.ShallowWaterResidualBasedBDFScheme(scheme_settings)
        elif self.settings["time_integration_scheme"].GetString().lower() == "adams-moulton":
            return SW.ResidualBasedAdamsMoultonScheme()
        else:
            raise Exception("Unknown time scheme, possible options are 'Adams-Moulton' or 'BDF'")

    @classmethod
    def GetDefaultParameters(cls):
        default_settings = KM.Parameters("""{
            "relative_dry_height"        : 0.1,
            "stabilization_factor"       : 0.01,
            "shock_capturing_factor"     : 0.0,
            "time_integration_scheme"    : "Adams-Moulton",
            "time_integration_order"     : 4
        }""")
        default_settings.AddMissingParameters(super().GetDefaultParameters())
        return default_settings

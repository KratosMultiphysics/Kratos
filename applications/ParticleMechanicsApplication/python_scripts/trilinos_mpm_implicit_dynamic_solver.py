# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.TrilinosApplication as TrilinosApplication

# Import base class file
from KratosMultiphysics.ParticleMechanicsApplication.trilinos_mpm_solver import TrilinosMPMSolver
from KratosMultiphysics.ParticleMechanicsApplication import TrilinosExtension as TrilinosExtension

def CreateSolver(model, custom_settings):
    return TrilinosMPMImplicitDynamicSolver(model, custom_settings)

class TrilinosMPMImplicitDynamicSolver(TrilinosMPMSolver):
    """The trilinos mpm static solver.

    For more information see:
    mpm_solver.py
    trilinos_mpm_solver.py
    """
    def __init__(self, model, custom_settings):
        # Construct the base solver.
        super(TrilinosMPMImplicitDynamicSolver, self).__init__(model, custom_settings)
        KratosMultiphysics.Logger.PrintInfo("::[TrilinosMPMStaticSolver]:: ", "Construction finished")

    @classmethod
    def GetDefaultParameters(cls):
        this_defaults = KratosMultiphysics.Parameters("""{
            "scheme_type"   : "bossak",
            "damp_factor_m" : -0.3,
            "newmark_beta"  : 0.25
        }""")
        this_defaults.AddMissingParameters(super(TrilinosMPMImplicitDynamicSolver, cls).GetDefaultParameters())
        return this_defaults

    def AddVariables(self):
        super(TrilinosMPMImplicitDynamicSolver, self).AddVariables()
        self._AddDynamicVariables(self.grid_model_part)
        KratosMultiphysics.Logger.PrintInfo("::[MPMImplicitDynamicSolver]:: ", "Variables are all added.")

    ### Protected functions ###
    def _CreateSolutionScheme(self):
        grid_model_part = self.GetGridModelPart()
        domain_size = self._GetDomainSize()
        block_size  = domain_size
        if (self.settings["pressure_dofs"].GetBool()):
            block_size += 1

        # Setting the time integration schemes
        scheme_type = self.settings["scheme_type"].GetString()
        if(scheme_type == "newmark"):
            damp_factor_m = 0.0
            newmark_beta = self.settings["newmark_beta"].GetDouble()
        elif(scheme_type == "bossak"):
            damp_factor_m = self.settings["damp_factor_m"].GetDouble()
            newmark_beta = self.settings["newmark_beta"].GetDouble()
        else:
            err_msg = "The requested scheme type \"" + scheme_type + "\" is not available!\n"
            err_msg += "Available options are: \"newmark\", \"bossak\""
            raise Exception(err_msg)

        is_dynamic = self._IsDynamic()

        return TrilinosExtension.TrilinosMPMResidualBasedBossakScheme( grid_model_part,
                                                                       domain_size,
                                                                       block_size,
                                                                       damp_factor_m,
                                                                       newmark_beta,
                                                                       is_dynamic)
    def _IsDynamic(self):
        return True
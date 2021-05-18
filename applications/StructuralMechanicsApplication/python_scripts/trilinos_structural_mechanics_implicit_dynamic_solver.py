# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import KratosMultiphysics.TrilinosApplication as TrilinosApplication

# Import base class file
from KratosMultiphysics.StructuralMechanicsApplication.trilinos_structural_mechanics_solver import TrilinosMechanicalSolver

from KratosMultiphysics.StructuralMechanicsApplication import auxiliar_methods_solvers

def CreateSolver(model, custom_settings):
    return TrilinosImplicitMechanicalSolver(model, custom_settings)

class TrilinosImplicitMechanicalSolver(TrilinosMechanicalSolver):
    """The trilinos structural mechanics implicit dynamic solver.

    For more information see:
    structural_mechanics_solver.py
    trilinos_structural_mechanics_solver.py
    """
    def __init__(self, model, custom_settings):
        # Construct the base solver.
        super().__init__(model, custom_settings)
        KratosMultiphysics.Logger.PrintInfo("::[TrilinosImplicitMechanicalSolver]:: ", "Construction finished")

    @classmethod
    def GetDefaultParameters(cls):
        this_defaults = KratosMultiphysics.Parameters("""{
            "time_integration_method" : "implicit",
            "scheme_type"             : "bossak",
            "damp_factor_m"           :-0.3,
            "rayleigh_alpha"          : 0.0,
            "rayleigh_beta"           : 0.0
        }""")
        this_defaults.AddMissingParameters(super().GetDefaultParameters())
        return this_defaults

    def AddVariables(self):
        super().AddVariables()
        self._add_dynamic_variables()
        KratosMultiphysics.Logger.PrintInfo("::[TrilinosImplicitMechanicalSolver]:: Variables ADDED")

    def AddDofs(self):
        super().AddDofs()
        self._add_dynamic_dofs()
        KratosMultiphysics.Logger.PrintInfo("::[TrilinosImplicitMechanicalSolver]:: DOF's ADDED")

    def GetMinimumBufferSize(self):
        base_min_buffer_size = super().GetMinimumBufferSize()

        scheme_type = self.settings["scheme_type"].GetString()
        if "bdf" in scheme_type or scheme_type == "backward_euler":
            return max(base_min_buffer_size, auxiliar_methods_solvers.GetBDFIntegrationOrder(scheme_type)+1)
        else:
            return base_min_buffer_size

    #### Private functions ####

    def _create_solution_scheme(self):
        scheme_type = self.settings["scheme_type"].GetString()
        process_info = self.main_model_part.ProcessInfo
        process_info[StructuralMechanicsApplication.RAYLEIGH_ALPHA] = self.settings["rayleigh_alpha"].GetDouble()
        process_info[StructuralMechanicsApplication.RAYLEIGH_BETA] = self.settings["rayleigh_beta"].GetDouble()
        if scheme_type == "newmark":
            damp_factor_m = 0.0
            mechanical_scheme = TrilinosApplication.TrilinosResidualBasedBossakDisplacementScheme(damp_factor_m)
        elif scheme_type == "bossak":
            damp_factor_m = self.settings["damp_factor_m"].GetDouble()
            mechanical_scheme = TrilinosApplication.TrilinosResidualBasedBossakDisplacementScheme(damp_factor_m)
        elif scheme_type.startswith("bdf") or scheme_type == "backward_euler" :
            order = auxiliar_methods_solvers.GetBDFIntegrationOrder(scheme_type)
            # In case of rotation dof we declare the dynamic variables
            if self.settings["rotation_dofs"].GetBool():
                bdf_parameters = KratosMultiphysics.Parameters(""" {
                    "domain_size"           : 3,
                    "integration_order"     : 2,
                    "solution_variables"    : ["DISPLACEMENT","ROTATION"]
                } """)
                bdf_parameters["domain_size"].SetInt(process_info[KratosMultiphysics.DOMAIN_SIZE])
                mechanical_scheme = TrilinosApplication.TrilinosResidualBasedBDFCustomScheme(order, bdf_parameters)
            else:
                mechanical_scheme = TrilinosApplication.TrilinosResidualBasedBDFDisplacementScheme(order)
        else:
            err_msg =  "The requested scheme type \"" + scheme_type + "\" is not available!\n"
            err_msg += "Available options are: \"newmark\", \"bossak\""
            raise Exception(err_msg)
        return mechanical_scheme

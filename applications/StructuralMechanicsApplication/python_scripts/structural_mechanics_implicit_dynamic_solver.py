from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication

# Import base class file
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_solver import MechanicalSolver

def CreateSolver(model, custom_settings):
    return ImplicitMechanicalSolver(model, custom_settings)

class ImplicitMechanicalSolver(MechanicalSolver):
    """The structural mechanics implicit dynamic solver.

    This class creates the mechanical solvers for implicit dynamic analysis.
    It currently supports Newmark, Bossak and dynamic relaxation schemes.

    See structural_mechanics_solver.py for more information.
    """
    def __init__(self, model, custom_settings):
        # Construct the base solver.
        super(ImplicitMechanicalSolver, self).__init__(model, custom_settings)
        KratosMultiphysics.Logger.PrintInfo("::[ImplicitMechanicalSolver]:: ", "Construction finished")

    @classmethod
    def GetDefaultSettings(cls):
        this_defaults = KratosMultiphysics.Parameters("""{
            "time_integration_method" : "implicit",
            "scheme_type"             : "bossak",
            "damp_factor_m"           :-0.3,
            "newmark_beta"            : 0.25,
            "rayleigh_alpha"          : 0.0,
            "rayleigh_beta"           : 0.0
        }""")
        this_defaults.AddMissingParameters(super(ImplicitMechanicalSolver, cls).GetDefaultSettings())
        return this_defaults

    def GetMinimumBufferSize(self):
        base_min_buffer_size = super(ImplicitMechanicalSolver, self).GetMinimumBufferSize()

        scheme_type = self.settings["scheme_type"].GetString()
        if "bdf" in scheme_type or scheme_type == "backward_euler":
            return max(base_min_buffer_size, self._bdf_integration_order()+1)
        else:
            return base_min_buffer_size

    def AddVariables(self):
        super(ImplicitMechanicalSolver, self).AddVariables()
        self._add_dynamic_variables()
        KratosMultiphysics.Logger.PrintInfo("::[ImplicitMechanicalSolver]:: ", "Variables ADDED")

    def AddDofs(self):
        super(ImplicitMechanicalSolver, self).AddDofs()
        self._add_dynamic_dofs()
        KratosMultiphysics.Logger.PrintInfo("::[ImplicitMechanicalSolver]:: ", "DOF's ADDED")

    def InitializeSolutionStep(self):
        # Using the base InitializeSolutionStep
        super(ImplicitMechanicalSolver, self).InitializeSolutionStep()

        # Some pre-processes may affect the system of equations, we rebuild the equation ids
        process_info = self.main_model_part.ProcessInfo
        if process_info[KratosMultiphysics.STEP] == 1 and process_info[StructuralMechanicsApplication.RESET_EQUATION_IDS]:
            # Resetting the global equations ids
            self.get_builder_and_solver().SetUpSystem(self.GetComputingModelPart())

    #### Private functions ####

    def _create_solution_scheme(self):
        scheme_type = self.settings["scheme_type"].GetString()

        # Setting the Rayleigh damping parameters
        process_info = self.main_model_part.ProcessInfo
        process_info[StructuralMechanicsApplication.RAYLEIGH_ALPHA] = self.settings["rayleigh_alpha"].GetDouble()
        process_info[StructuralMechanicsApplication.RAYLEIGH_BETA] = self.settings["rayleigh_beta"].GetDouble()

        # Setting the time integration schemes
        if(scheme_type == "newmark"):
            damp_factor_m = 0.0
            newmark_beta = self.settings["newmark_beta"].GetDouble()
            mechanical_scheme = KratosMultiphysics.ResidualBasedBossakDisplacementScheme(damp_factor_m, newmark_beta)
        elif(scheme_type == "bossak"):
            damp_factor_m = self.settings["damp_factor_m"].GetDouble()
            newmark_beta = self.settings["newmark_beta"].GetDouble()
            mechanical_scheme = KratosMultiphysics.ResidualBasedBossakDisplacementScheme(damp_factor_m, newmark_beta)
        elif(scheme_type == "pseudo_static"):
            mechanical_scheme = KratosMultiphysics.ResidualBasedPseudoStaticDisplacementScheme(StructuralMechanicsApplication.RAYLEIGH_BETA)
        elif(scheme_type.startswith("bdf") or scheme_type == "backward_euler"):
            order = self._bdf_integration_order()
            # In case of rotation dof we declare the dynamic variables
            if self.settings["rotation_dofs"].GetBool():
                bdf_parameters = KratosMultiphysics.Parameters(""" {
                    "domain_size"           : 3,
                    "integration_order"     : 2,
                    "variable"              : ["DISPLACEMENT","ROTATION"],
                    "first_derivative"      : ["VELOCITY","ANGULAR_VELOCITY"],
                    "second_derivative"     : ["ACCELERATION","ANGULAR_ACCELERATION"]
                } """)
                bdf_parameters["domain_size"].SetInt(process_info[KratosMultiphysics.DOMAIN_SIZE])
                mechanical_scheme = KratosMultiphysics.ResidualBasedBDFCustomScheme(order, bdf_parameters)
            else:
                mechanical_scheme = KratosMultiphysics.ResidualBasedBDFDisplacementScheme(order)
        elif(scheme_type == "relaxation"):
            damp_factor_f =-0.3
            dynamic_factor_m = 10.0
            mechanical_scheme = StructuralMechanicsApplication.ResidualBasedRelaxationScheme(damp_factor_f, dynamic_factor_m)
        else:
            err_msg = "The requested scheme type \"" + scheme_type + "\" is not available!\n"
            err_msg += "Available options are: \"newmark\", \"bossak\", \"pseudo_static\", \"backward_euler\", \"bdf1\", \"bdf2\", \"bdf3\", \"bdf4\", \"bdf5\", \"relaxation\""
            raise Exception(err_msg)
        return mechanical_scheme

    def _bdf_integration_order(self):
        scheme_type = self.settings["scheme_type"].GetString()
        if scheme_type == "backward_euler":
            order = 1
        else:
            if scheme_type == "bdf":
                raise Exception('Wrong input for scheme type: "bdf"! Please append the order to the bdf-scheme, e.g. "bdf2"')
            # BDF schemes can be from 1 to 5 order, so in order to detect the integration order from the scheme_type we remove the "bdf" string, that is, if the user tells bdf3 only 3 will remain when we remove bdf which corresponds to the method of choice
            order = int(scheme_type.replace("bdf", ""))

        # Warning
        if (order > 2):
            KratosMultiphysics.Logger.PrintWarning("WARNING:: BDF Order: ", str(order) + " constant time step must be considered")

        return order

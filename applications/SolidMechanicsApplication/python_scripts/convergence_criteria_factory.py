## This script collects the available convergence criteria to be used in the SolidMechanicsApplication

from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.SolidMechanicsApplication as KratosSolid

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

# Convergence criteria class
class ConvergenceCriterion:

    def __init__(self, custom_settings):
        """
        Create a convergence criterion from json parameters.
        """
        default_settings = KratosMultiphysics.Parameters("""
        {
               "convergence_criterion": "Residual_criterion",
               "variable_relative_tolerance": 1.0e-4,
               "variable_absolute_tolerance": 1.0e-9,
               "residual_relative_tolerance": 1.0e-4,
               "residual_absolute_tolerance": 1.0e-9,
               "echo_level": 0
        }
        """)

        # Overwrite the default settings with user-provided parameters
        self.settings = custom_settings
        self.settings.ValidateAndAssignDefaults(default_settings)


        # Note that all the convergence settings are introduced via a Kratos parameters object.
        V_RT = self.settings["variable_relative_tolerance"].GetDouble()
        V_AT = self.settings["variable_absolute_tolerance"].GetDouble()
        R_RT = self.settings["residual_relative_tolerance"].GetDouble()
        R_AT = self.settings["residual_absolute_tolerance"].GetDouble()

        echo_level = self.settings["echo_level"].GetInt()

        if(echo_level >= 1):
            print("::[Mechanical_Solver]:: Convergence criterion [", self.settings["convergence_criterion"].GetString(),"]")

        self.convergence_criterion = None
        if(self.settings["convergence_criterion"].GetString() == "Variable_criterion"):
            self.convergence_criterion = KratosSolid.VariableCriterion(V_RT, V_AT)
            self.convergence_criterion.SetEchoLevel(echo_level)
            self.convergence_criterion.Set(KratosSolid.CriterionLocalFlags.INCREMENTAL)
        elif(self.settings["convergence_criterion"].GetString() == "Residual_criterion"):
            self.convergence_criterion = KratosSolid.ResidualCriterion(R_RT, R_AT)
            self.convergence_criterion.SetEchoLevel(echo_level)
        elif(self.settings["convergence_criterion"].GetString() == "And_criterion"):
            Variable = KratosSolid.VariableCriterion(V_RT, V_AT)
            Variable.SetEchoLevel(echo_level)
            Variable.Set(KratosSolid.CriterionLocalFlags.INCREMENTAL)
            Residual = KratosSolid.ResidualCriterion(R_RT, R_AT)
            Residual.SetEchoLevel(echo_level)
            self.convergence_criterion = KratosSolid.CompositeCriterion(Residual, Variable)
            self.convergence_criterion.Set(KratosSolid.CriterionLocalFlags.AND)
        elif(self.settings["convergence_criterion"].GetString() == "Or_criterion"):
            Variable = KratosSolid.VariableCriterion(V_RT, V_AT)
            Variable.SetEchoLevel(echo_level)
            Variable.Set(KratosSolid.CriterionLocalFlags.INCREMENTAL)
            Variable = KratosSolid.ResidualCriterion(R_RT, R_AT)
            Residual.SetEchoLevel(echo_level)
            self.convergence_criterion = KratosSolid.CompositeCriterion(Residual, Variable)
            self.convergence_criterion.Set(KratosSolid.CriterionLocalFlags.OR)
        else:
            raise Exception("Unsupported \"convergence_criterion\" : " + self.settings["convergence_criterion"].GetString())


    #
    def GetConvergenceCriterion(self):
        return self.convergence_criterion

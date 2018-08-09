## This script collects the available convergence criteria to be used in the SolidMechanicsApplication

from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.SolidMechanicsApplication as KratosSolid

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

# Convergence criteria class
class ConvergenceCriterion:

    def __init__(self, custom_settings, dofs_list):
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
               "separate_dofs" : true,
               "echo_level": 0
        }
        """)

        # Overwrite the default settings with user-provided parameters
        self.settings = custom_settings
        self.settings.ValidateAndAssignDefaults(default_settings)

        self.dofs = []
        for i in range(0, dofs_list.size() ):
            self.dofs.append(dofs_list[i].GetString())

        # add default DISPLACEMENT dof
        if( len(self.dofs) == 0 or (len(self.dofs) == 1 and self.dofs[0] =="ROTATION") ):
            self.dofs.append('DISPLACEMENT')
    #
    def GetConvergenceCriterion(self):

        if( len(self.dofs) > 0 and self.settings["separate_dofs"].GetBool() ):
            convergence_criteria = []
            for dof in self.dofs:
                kratos_variable = KratosMultiphysics.KratosGlobals.GetVariable(dof)
                convergence_criteria.append(self._get_dof_criterion(kratos_variable))

            convergence_criterion = KratosSolid.CompositeCriterion(convergence_criteria)
            convergence_criterion.Set(KratosSolid.CriterionLocalFlags.AND)
        else:
            convergence_criterion = self._get_criterion()

        return convergence_criterion
    #
    def _get_dof_criterion(self, kratos_variable):

        # Note that all the convergence settings are introduced via a Kratos parameters object.
        V_RT = self.settings["variable_relative_tolerance"].GetDouble()
        V_AT = self.settings["variable_absolute_tolerance"].GetDouble()
        R_RT = self.settings["residual_relative_tolerance"].GetDouble()
        R_AT = self.settings["residual_absolute_tolerance"].GetDouble()

        echo_level = self.settings["echo_level"].GetInt()

        convergence_criterion = None
        if(self.settings["convergence_criterion"].GetString() == "Variable_criterion"):
            convergence_criterion = KratosSolid.DofsCriterion(kratos_variable,V_RT,V_AT)
            convergence_criterion.SetEchoLevel(echo_level)
            #convergence_criterion.Set(KratosSolid.CriterionLocalFlags.INCREMENTAL) //smaller reference value
        elif(self.settings["convergence_criterion"].GetString() == "Residual_criterion"):
            convergence_criterion = KratosSolid.ResidualCriterion(kratos_variable,R_RT,R_AT)
            convergence_criterion.SetEchoLevel(echo_level)
        elif(self.settings["convergence_criterion"].GetString() == "And_criterion"):
            Variable = KratosSolid.DofsCriterion(kratos_variable,V_RT,V_AT)
            Variable.SetEchoLevel(echo_level)
            #Variable.Set(KratosSolid.CriterionLocalFlags.INCREMENTAL) //smaller reference value
            Residual = KratosSolid.ResidualCriterion(kratos_variable,R_RT,R_AT)
            Residual.SetEchoLevel(echo_level)
            convergence_criterion = KratosSolid.CompositeCriterion(Residual,Variable)
            convergence_criterion.Set(KratosSolid.CriterionLocalFlags.AND)
        elif(self.settings["convergence_criterion"].GetString() == "Or_criterion"):
            Variable = KratosSolid.DofsCriterion(kratos_variable,V_RT,V_AT)
            Variable.SetEchoLevel(echo_level)
            #Variable.Set(KratosSolid.CriterionLocalFlags.INCREMENTAL) //smaller reference value
            Residual = KratosSolid.ResidualCriterion(kratos_variable,R_RT,R_AT)
            Residual.SetEchoLevel(echo_level)
            convergence_criterion = KratosSolid.CompositeCriterion(Residual,Variable)
            convergence_criterion.Set(KratosSolid.CriterionLocalFlags.OR)
        else:
            raise Exception("Unsupported \"convergence_criterion\" : " + self.settings["convergence_criterion"].GetString())

        return convergence_criterion

    #
    def _get_criterion(self):

        # Note that all the convergence settings are introduced via a Kratos parameters object.
        V_RT = self.settings["variable_relative_tolerance"].GetDouble()
        V_AT = self.settings["variable_absolute_tolerance"].GetDouble()
        R_RT = self.settings["residual_relative_tolerance"].GetDouble()
        R_AT = self.settings["residual_absolute_tolerance"].GetDouble()

        echo_level = self.settings["echo_level"].GetInt()

        if(echo_level >= 1):
            print("::[-----Criterion-----]::", self.settings["convergence_criterion"].GetString(),)

        convergence_criterion = None
        if(self.settings["convergence_criterion"].GetString() == "Variable_criterion"):
            convergence_criterion = KratosSolid.DofsCriterion(V_RT,V_AT)
            convergence_criterion.SetEchoLevel(echo_level)
            convergence_criterion.Set(KratosSolid.CriterionLocalFlags.INCREMENTAL)
        elif(self.settings["convergence_criterion"].GetString() == "Residual_criterion"):
            convergence_criterion = KratosSolid.ResidualCriterion(R_RT,R_AT)
            convergence_criterion.SetEchoLevel(echo_level)
        elif(self.settings["convergence_criterion"].GetString() == "And_criterion"):
            Variable = KratosSolid.DofsCriterion(V_RT,V_AT)
            Variable.SetEchoLevel(echo_level)
            Variable.Set(KratosSolid.CriterionLocalFlags.INCREMENTAL)
            Residual = KratosSolid.ResidualCriterion(R_RT,R_AT)
            Residual.SetEchoLevel(echo_level)
            convergence_criterion = KratosSolid.CompositeCriterion(Residual,Variable)
            convergence_criterion.Set(KratosSolid.CriterionLocalFlags.AND)
        elif(self.settings["convergence_criterion"].GetString() == "Or_criterion"):
            Variable = KratosSolid.DofsCriterion(V_RT,V_AT)
            Variable.SetEchoLevel(echo_level)
            Variable.Set(KratosSolid.CriterionLocalFlags.INCREMENTAL)
            Residual = KratosSolid.ResidualCriterion(R_RT,R_AT)
            Residual.SetEchoLevel(echo_level)
            convergence_criterion = KratosSolid.CompositeCriterion(Residual,Variable)
            convergence_criterion.Set(KratosSolid.CriterionLocalFlags.OR)
        else:
            raise Exception("Unsupported \"convergence_criterion\" : " + self.settings["convergence_criterion"].GetString())

        return convergence_criterion

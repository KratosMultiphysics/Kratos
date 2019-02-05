## This script collects the available convergence criteria to be used in the SolidMechanicsApplication

from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.SolidMechanicsApplication as KratosSolid

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

# Convergence criteria class
class convergence_criterion:

    def __init__(self, convergence_criterion_parameters):
        """Create a convergence criterion from json parameters.

        Keyword arguments:
        convergence_criterion_parameters

        If no error is raised, a valid convergence criterion is accessible
        through the member variable mechanical_convergence_criterion after
        return.
        """
        # Note that all the convergence settings are introduced via a Kratos parameters object.
        D_RT = convergence_criterion_parameters["displacement_relative_tolerance"].GetDouble()
        D_AT = convergence_criterion_parameters["displacement_absolute_tolerance"].GetDouble()
        R_RT = convergence_criterion_parameters["residual_relative_tolerance"].GetDouble()
        R_AT = convergence_criterion_parameters["residual_absolute_tolerance"].GetDouble()

        echo_level = convergence_criterion_parameters["echo_level"].GetInt()

        if(echo_level >= 1):
            print("::[Mechanical_Solver]:: Convergence criterion [", convergence_criterion_parameters["convergence_criterion"].GetString(),"]")

        rotation_dofs = False
        if(convergence_criterion_parameters.Has("rotation_dofs")):
            if(convergence_criterion_parameters["rotation_dofs"].GetBool()):
                rotation_dofs = True

        component_wise = False
        if(convergence_criterion_parameters.Has("component_wise")):
            if(convergence_criterion_parameters["component_wise"].GetBool()):
                component_wise = True

        # Convergence criteria if there are rotation DOFs in the problem
        if(rotation_dofs == True):
            if(convergence_criterion_parameters["convergence_criterion"].GetString() == "Displacement_criterion"):
                self.mechanical_convergence_criterion = KratosSolid.DisplacementCriteria(D_RT, D_AT)
                self.mechanical_convergence_criterion.SetEchoLevel(echo_level)
            elif(convergence_criterion_parameters["convergence_criterion"].GetString() == "Residual_criterion"):
                self.mechanical_convergence_criterion = KratosMultiphysics.ResidualCriteria(R_RT, R_AT)
                self.mechanical_convergence_criterion.SetEchoLevel(echo_level)
            elif(convergence_criterion_parameters["convergence_criterion"].GetString() == "And_criterion"):
                Displacement = KratosSolid.DisplacementCriteria(D_RT, D_AT)
                Displacement.SetEchoLevel(echo_level)
                Residual = KratosMultiphysics.ResidualCriteria(R_RT, R_AT)
                Residual.SetEchoLevel(echo_level)
                self.mechanical_convergence_criterion = KratosMultiphysics.AndCriteria(Residual, Displacement)
            elif(convergence_criterion_parameters["convergence_criterion"].GetString() == "Or_criterion"):
                Displacement = KratosSolid.DisplacementCriteria(D_RT, D_AT)
                Displacement.SetEchoLevel(echo_level)
                Residual = KratosMultiphysics.ResidualCriteria(R_RT, R_AT)
                Residual.SetEchoLevel(echo_level)
                self.mechanical_convergence_criterion = KratosMultiphysics.OrCriteria(Residual, Displacement)
            else:
                raise Exception("Unsupported \"convergence_criterion\" : " + convergence_criterion_parameters["convergence_criterion"].GetString())

        # Convergence criteria without rotation DOFs
        else:
            if(convergence_criterion_parameters["convergence_criterion"].GetString() == "Displacement_criterion"):
                self.mechanical_convergence_criterion = KratosSolid.DisplacementConvergenceCriterion(D_RT, D_AT)
                self.mechanical_convergence_criterion.SetEchoLevel(echo_level)
            elif(convergence_criterion_parameters["convergence_criterion"].GetString() == "Residual_criterion"):
                if(component_wise == True):
                    self.mechanical_convergence_criterion = KratosSolid.ComponentWiseResidualConvergenceCriterion(R_RT, R_AT)
                else:
                    self.mechanical_convergence_criterion = KratosMultiphysics.ResidualCriteria(R_RT, R_AT)
                self.mechanical_convergence_criterion.SetEchoLevel(echo_level)
            elif(convergence_criterion_parameters["convergence_criterion"].GetString() == "And_criterion"):
                Displacement = KratosSolid.DisplacementConvergenceCriterion(D_RT, D_AT)
                Displacement.SetEchoLevel(echo_level)
                if(component_wise == True):
                    Residual = KratosSolid.ComponentWiseResidualConvergenceCriterion(R_RT, R_AT)
                else:
                    Residual = KratosMultiphysics.ResidualCriteria(R_RT, R_AT)
                Residual.SetEchoLevel(echo_level)
                self.mechanical_convergence_criterion = KratosMultiphysics.AndCriteria(Residual, Displacement)
            elif(convergence_criterion_parameters["convergence_criterion"].GetString() == "Or_criterion"):
                Displacement = KratosSolid.DisplacementConvergenceCriterion(D_RT, D_AT)
                Displacement.SetEchoLevel(echo_level)
                if(component_wise == True):
                    Residual = KratosSolid.ComponentWiseResidualConvergenceCriterion(R_RT, R_AT)
                else:
                    Residual = KratosMultiphysics.ResidualCriteria(R_RT, R_AT)
                Residual.SetEchoLevel(echo_level)
                self.mechanical_convergence_criterion = KratosMultiphysics.OrCriteria(Residual, Displacement)
            else:
                raise Exception("Unsupported \"convergence_criterion\" : " + convergence_criterion_parameters["convergence_criterion"].GetString())

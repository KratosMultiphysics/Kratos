## This script collects the available convergence criteria to be used in the SolidMechanicsApplication

#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.FemToDemApplication as KratosFemDem

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


        if(convergence_criterion_parameters["convergence_criterion"].GetString() == "FemDem_Residual_criterion"):
            self.mechanical_convergence_criterion = KratosFemDem.FemDemResidualCriteria(R_RT, R_AT)
            self.mechanical_convergence_criterion.SetEchoLevel(echo_level)
        else:
            raise Exception("Unsupported \"convergence_criterion\" : " + convergence_criterion_parameters["convergence_criterion"].GetString())

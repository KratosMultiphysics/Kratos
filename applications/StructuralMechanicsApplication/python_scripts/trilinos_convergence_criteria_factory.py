# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.TrilinosApplication as TrilinosApplication


# Convergence criteria class
class convergence_criterion:
    def __init__(self, convergence_criterion_parameters):
        # Note that all the convergence settings are introduced via a Kratos parameters object.

        D_RT = convergence_criterion_parameters["displacement_relative_tolerance"].GetDouble()
        D_AT = convergence_criterion_parameters["displacement_absolute_tolerance"].GetDouble()
        R_RT = convergence_criterion_parameters["residual_relative_tolerance"].GetDouble()
        R_AT = convergence_criterion_parameters["residual_absolute_tolerance"].GetDouble()

        convergence_crit = convergence_criterion_parameters["convergence_criterion"].GetString()

        echo_level = convergence_criterion_parameters["echo_level"].GetInt()

        rotation_dofs = False
        if convergence_criterion_parameters.Has("rotation_dofs"):
            rotation_dofs = convergence_criterion_parameters["rotation_dofs"].GetBool()

        volumetric_strain_dofs = False
        if convergence_criterion_parameters.Has("volumetric_strain_dofs"):
            volumetric_strain_dofs = convergence_criterion_parameters["volumetric_strain_dofs"].GetBool()

        if(echo_level >= 1):
            KratosMultiphysics.Logger.PrintInfo("::[Mechanical Solver]::", "MPI CONVERGENCE CRITERION : " + convergence_criterion_parameters["convergence_criterion"].GetString())

        if(convergence_crit == "displacement_criterion"):
            if rotation_dofs:
                self.mechanical_convergence_criterion = TrilinosApplication.TrilinosMixedGenericCriteria(
                    [(KratosMultiphysics.DISPLACEMENT, D_RT, D_AT),
                    (KratosMultiphysics.ROTATION, D_RT, D_AT)])
            elif volumetric_strain_dofs:
                self.mechanical_convergence_criterion = TrilinosApplication.TrilinosMixedGenericCriteria(
                    [(KratosMultiphysics.DISPLACEMENT, D_RT, D_AT),
                    (KratosMultiphysics.VOLUMETRIC_STRAIN, D_RT, D_AT)])
            else:
                self.mechanical_convergence_criterion = TrilinosApplication.TrilinosDisplacementCriteria(D_RT, D_AT)
            self.mechanical_convergence_criterion.SetEchoLevel(echo_level)

        elif(convergence_crit == "residual_criterion"):
            self.mechanical_convergence_criterion = TrilinosApplication.TrilinosResidualCriteria(R_RT, R_AT)
            self.mechanical_convergence_criterion.SetEchoLevel(echo_level)

        elif(convergence_crit == "and_criterion"):
            Displacement = TrilinosApplication.TrilinosDisplacementCriteria(D_RT, D_AT)
            Displacement.SetEchoLevel(echo_level)

            Residual = TrilinosApplication.TrilinosResidualCriteria(R_RT, R_AT)
            Residual.SetEchoLevel(echo_level)
            self.mechanical_convergence_criterion = TrilinosApplication.TrilinosAndCriteria(Residual, Displacement)

        elif(convergence_crit == "or_criterion"):
            Displacement = TrilinosApplication.TrilinosDisplacementCriteria(D_RT, D_AT)
            Displacement.SetEchoLevel(echo_level)

            Residual = TrilinosApplication.TrilinosResidualCriteria(R_RT, R_AT)
            Residual.SetEchoLevel(echo_level)
            self.mechanical_convergence_criterion = TrilinosApplication.TrilinosOrCriteria(Residual, Displacement)

        else:
            err_msg =  "The requested convergence criterion \"" + convergence_crit + "\" is not available!\n"
            err_msg += "Available options are: \"displacement_criterion\", \"residual_criterion\", \"and_criterion\", \"or_criterion\""
            raise Exception(err_msg)

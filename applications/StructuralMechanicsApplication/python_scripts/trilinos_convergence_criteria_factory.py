from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication

# MPI
import KratosMultiphysics.mpi as mpi
import KratosMultiphysics.TrilinosApplication as TrilinosApplication
import KratosMultiphysics.MetisApplication as MetisApplication

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

# Convergence criteria class
class convergence_criterion:
    def __init__(self, convergence_criterion_parameters):
        # Note that all the convergence settings are introduced via a Kratos parameters object.
        
        D_RT = convergence_criterion_parameters["displacement_relative_tolerance"].GetDouble()
        D_AT = convergence_criterion_parameters["displacement_absolute_tolerance"].GetDouble()
        R_RT = convergence_criterion_parameters["residual_relative_tolerance"].GetDouble()
        R_AT = convergence_criterion_parameters["residual_absolute_tolerance"].GetDouble()
        
        echo_level = convergence_criterion_parameters["echo_level"].GetInt()
        
        if(echo_level >= 1):
            print("::[Mechanical Solver]:: MPI CONVERGENCE CRITERION : ", convergence_criterion_parameters["convergence_criterion"].GetString())

        if(convergence_criterion_parameters["convergence_criterion"].GetString() == "displacement_criterion"):
            self.mechanical_convergence_criterion = TrilinosApplication.TrilinosDisplacementCriteria(D_RT, D_AT)
            self.mechanical_convergence_criterion.SetEchoLevel(echo_level)
            
        elif(convergence_criterion_parameters["convergence_criterion"].GetString() == "residual_criterion"):
            self.mechanical_convergence_criterion = TrilinosApplication.TrilinosResidualCriteria(R_RT, R_AT)
            self.mechanical_convergence_criterion.SetEchoLevel(echo_level)
                
        elif(convergence_criterion_parameters["convergence_criterion"].GetString() == "and_criterion"):
            Displacement = TrilinosApplication.TrilinosDisplacementCriteria(D_RT, D_AT)
            Displacement.SetEchoLevel(echo_level)

            Residual = TrilinosApplication.TrilinosResidualCriteria(R_RT, R_AT)
            Residual.SetEchoLevel(echo_level)
            self.mechanical_convergence_criterion = TrilinosApplication.TrilinosAndCriteria(Residual, Displacement)
            
        elif(convergence_criterion_parameters["convergence_criterion"].GetString() == "or_criterion"):
            Displacement = TrilinosApplication.TrilinosDisplacementCriteria(D_RT, D_AT)
            Displacement.SetEchoLevel(echo_level)

            Residual = TrilinosApplication.TrilinosResidualCriteria(R_RT, R_AT)
            Residual.SetEchoLevel(echo_level)
            self.mechanical_convergence_criterion = TrilinosApplication.TrilinosOrCriteria(Residual, Displacement)
        


## This script collects the available convergence criteria to be used in the SolidMechanicsApplication

from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.SolidMechanicsApplication as SolidMechanicsApplication

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
        
        if(convergence_criterion_parameters["echo_level"].GetInt() > 1):
                print("::[Mechanical Solver]:: CONVERGENCE CRITERION : ", convergence_criterion_parameters["convergence_criterion"].GetString())

        # Convergence criteria if there are rotation DOFs in the problem
        if(convergence_criterion_parameters["rotation_dofs"].GetBool()):
            if(convergence_criterion_parameters["convergence_criterion"].GetString() == "Displacement_criterion"):
                self.mechanical_convergence_criterion = SolidMechanicsApplication.DisplacementCriteria(D_RT, D_AT)
                
            elif(convergence_criterion_parameters["convergence_criterion"].GetString() == "Residual_criterion"):
                self.mechanical_convergence_criterion = KratosMultiphysics.ResidualCriteria(R_RT, R_AT)
                
            elif(convergence_criterion_parameters["convergence_criterion"].GetString() == "And_criterion"):
                Displacement = SolidMechanicsApplication.DisplacementCriteria(D_RT, D_AT)
                Residual = KratosMultiphysics.ResidualCriteria(R_RT, R_AT)
                self.mechanical_convergence_criterion = KratosMultiphysics.AndCriteria(Residual, Displacement)
                
            elif(convergence_criterion_parameters["convergence_criterion"].GetString() == "Or_criterion"):
                Displacement = SolidMechanicsApplication.DisplacementCriteria(D_RT, D_AT)
                Residual = KratosMultiphysics.ResidualCriteria(R_RT, R_AT)
                self.mechanical_convergence_criterion = KratosMultiphysics.OrCriteria(Residual, Displacement)
            
        # Convergence criteria without rotation DOFs        
        else:
            if(convergence_criterion_parameters["convergence_criterion"].GetString() == "Displacement_criterion"):
                self.mechanical_convergence_criterion = SolidMechanicsApplication.DisplacementConvergenceCriterion(D_RT, D_AT)
                
            elif(convergence_criterion_parameters["convergence_criterion"].GetString() == "Residual_criterion"):
                if(convergence_criterion_parameters["component_wise"].GetBool()):
                    self.mechanical_convergence_criterion = SolidMechanicsApplication.ComponentWiseResidualConvergenceCriterion(R_RT, R_AT)
                else:
                    self.mechanical_convergence_criterion = KratosMultiphysics.ResidualCriteria(R_RT, R_AT)
                    
            elif(convergence_criterion_parameters["convergence_criterion"].GetString() == "And_criterion"):
                Displacement = SolidMechanicsApplication.DisplacementConvergenceCriterion(D_RT, D_AT)
                if(convergence_criterion_parameters["component_wise"].GetBool()):
                    Residual = SolidMechanicsApplication.ComponentWiseResidualConvergenceCriterion(R_RT, R_AT)
                else:
                    Residual = KratosMultiphysics.ResidualCriteria(R_RT, R_AT)
                self.mechanical_convergence_criterion = KratosMultiphysics.AndCriteria(Residual, Displacement)
                
            elif(convergence_criterion_parameters["convergence_criterion"].GetString() == "Or_criterion"):
                Displacement = SolidMechanicsApplication.DisplacementConvergenceCriterion(D_RT, D_AT)
                if(convergence_criterion_parameters["component_wise"].GetBool()):
                    Residual = SolidMechanicsApplication.ComponentWiseResidualConvergenceCriterion(R_RT, R_AT)
                else:
                    Residual = KratosMultiphysics.ResidualCriteria(R_RT, R_AT)
                self.mechanical_convergence_criterion = KratosMultiphysics.OrCriteria(Residual, Displacement)
        


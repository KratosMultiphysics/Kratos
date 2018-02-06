## This script collects the available convergence criteria to be used in the SolidMechanicsApplication

from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import KratosMultiphysics.ContactStructuralMechanicsApplication as ContactStructuralMechanicsApplication
import KratosMultiphysics.ExternalSolversApplication as ExternalSolversApplication

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

# Convergence criteria class
class convergence_criterion:
    def __init__(self, convergence_criterion_parameters):
        # Note that all the convergence settings are introduced via a Kratos parameters object.
        echo_level = convergence_criterion_parameters["echo_level"].GetInt()
        if "contact" in convergence_criterion_parameters["convergence_criterion"].GetString():
            D_RT = convergence_criterion_parameters["displacement_relative_tolerance"].GetDouble()
            D_AT = convergence_criterion_parameters["displacement_absolute_tolerance"].GetDouble()
            R_RT = convergence_criterion_parameters["residual_relative_tolerance"].GetDouble()
            R_AT = convergence_criterion_parameters["residual_absolute_tolerance"].GetDouble()
            CD_RT = convergence_criterion_parameters["contact_displacement_relative_tolerance"].GetDouble()
            CD_AT = convergence_criterion_parameters["contact_displacement_absolute_tolerance"].GetDouble()
            CR_RT = convergence_criterion_parameters["contact_residual_relative_tolerance"].GetDouble()
            CR_AT = convergence_criterion_parameters["contact_residual_absolute_tolerance"].GetDouble()
            condn_convergence_criterion = convergence_criterion_parameters["condn_convergence_criterion"].GetBool()
            fancy_convergence_criterion = convergence_criterion_parameters["fancy_convergence_criterion"].GetBool()
            print_convergence_criterion = convergence_criterion_parameters["print_convergence_criterion"].GetBool()
            ensure_contact = convergence_criterion_parameters["ensure_contact"].GetBool()
            gidio_debug = convergence_criterion_parameters["gidio_debug"].GetBool()
            
            if(echo_level >= 1):
                print("::[Mechanical Solver]:: CONVERGENCE CRITERION : ", convergence_criterion_parameters["convergence_criterion"].GetString())
            
            if (fancy_convergence_criterion == True):
                table = KratosMultiphysics.TableStreamUtility()
            
            if(convergence_criterion_parameters["convergence_criterion"].GetString() == "contact_displacement_criterion"):
                if (fancy_convergence_criterion == True):
                    self.mechanical_convergence_criterion = ContactStructuralMechanicsApplication.DisplacementLagrangeMultiplierContactCriteria(D_RT, D_AT, D_RT, D_AT, ensure_contact, table, print_convergence_criterion)
                else:
                    self.mechanical_convergence_criterion = ContactStructuralMechanicsApplication.DisplacementLagrangeMultiplierContactCriteria(D_RT, D_AT, D_RT, D_AT, ensure_contact)
                self.mechanical_convergence_criterion.SetEchoLevel(echo_level)
                
            elif(convergence_criterion_parameters["convergence_criterion"].GetString() == "contact_residual_criterion"):
                if (fancy_convergence_criterion == True):
                    self.mechanical_convergence_criterion = ContactStructuralMechanicsApplication.DisplacementLagrangeMultiplierResidualContactCriteria(R_RT, R_AT, CR_RT, CR_AT, ensure_contact, table, print_convergence_criterion)
                else:
                    self.mechanical_convergence_criterion = ContactStructuralMechanicsApplication.DisplacementLagrangeMultiplierResidualContactCriteria(R_RT, R_AT, CR_RT, CR_AT, ensure_contact)
                self.mechanical_convergence_criterion.SetEchoLevel(echo_level)
                
            elif(convergence_criterion_parameters["convergence_criterion"].GetString() == "contact_mixed_criterion"):
                if (fancy_convergence_criterion == True):
                    self.mechanical_convergence_criterion = ContactStructuralMechanicsApplication.DisplacementLagrangeMultiplierMixedContactCriteria(R_RT, R_AT, CR_RT, CR_AT, ensure_contact, table, print_convergence_criterion)
                else:
                    self.mechanical_convergence_criterion = ContactStructuralMechanicsApplication.DisplacementLagrangeMultiplierMixedContactCriteria(R_RT, R_AT, CR_RT, CR_AT, ensure_contact)
                self.mechanical_convergence_criterion.SetEchoLevel(echo_level)
                    
            elif(convergence_criterion_parameters["convergence_criterion"].GetString() == "contact_and_criterion"):
                if (fancy_convergence_criterion == True):
                    Displacement = ContactStructuralMechanicsApplication.DisplacementLagrangeMultiplierContactCriteria(D_RT, D_AT, CD_RT, CD_AT, ensure_contact, table, print_convergence_criterion)
                    Residual = ContactStructuralMechanicsApplication.DisplacementLagrangeMultiplierResidualContactCriteria(R_RT, R_AT, CR_RT, CR_AT, ensure_contact, table, print_convergence_criterion)
                else:
                    Displacement = ContactStructuralMechanicsApplication.DisplacementLagrangeMultiplierContactCriteria(D_RT, D_AT, CD_RT, CD_AT, ensure_contact)
                    Residual = ContactStructuralMechanicsApplication.DisplacementLagrangeMultiplierResidualContactCriteria(R_RT, R_AT, CR_RT, CR_AT, ensure_contact)
                
                Displacement.SetEchoLevel(echo_level)
                Residual.SetEchoLevel(echo_level)
                self.mechanical_convergence_criterion = KratosMultiphysics.AndCriteria(Residual, Displacement)
                
            elif(convergence_criterion_parameters["convergence_criterion"].GetString() == "contact_or_criterion"):
                if (fancy_convergence_criterion == True):
                    Displacement = ContactStructuralMechanicsApplication.DisplacementLagrangeMultiplierContactCriteria(D_RT, D_AT, CD_RT, CD_AT, ensure_contact, table, print_convergence_criterion)
                    Residual = ContactStructuralMechanicsApplication.DisplacementLagrangeMultiplierResidualContactCriteria(R_RT, R_AT, CR_RT, CR_AT, ensure_contact, table, print_convergence_criterion)
                else:
                    Displacement = ContactStructuralMechanicsApplication.DisplacementLagrangeMultiplierContactCriteria(D_RT, D_AT, CD_RT, CD_AT,ensure_contact)
                    Residual = ContactStructuralMechanicsApplication.DisplacementLagrangeMultiplierResidualContactCriteria(R_RT, R_AT, CR_RT, CR_AT, ensure_contact)
                
                Displacement.SetEchoLevel(echo_level)
                Residual.SetEchoLevel(echo_level)
                self.mechanical_convergence_criterion = KratosMultiphysics.OrCriteria(Residual, Displacement)
            
            # Adding the mortar criteria
            if  (convergence_criterion_parameters["mortar_type"].GetString() == "ALMContactFrictionless"):
                if (fancy_convergence_criterion == True):
                    Mortar = ContactStructuralMechanicsApplication.ALMFrictionlessMortarConvergenceCriteria(table, print_convergence_criterion, gidio_debug)
                else:
                    Mortar = ContactStructuralMechanicsApplication.ALMFrictionlessMortarConvergenceCriteria()
            elif  (convergence_criterion_parameters["mortar_type"].GetString() == "ALMContactFrictional"):
                if (fancy_convergence_criterion == True):
                    Mortar = ContactStructuralMechanicsApplication.ALMFrictionalMortarConvergenceCriteria(table, print_convergence_criterion, gidio_debug)
                else:
                    Mortar = ContactStructuralMechanicsApplication.ALMFrictionalMortarConvergenceCriteria()
            elif ("MeshTying" in convergence_criterion_parameters["mortar_type"].GetString()):
                if (fancy_convergence_criterion == True):
                    Mortar = ContactStructuralMechanicsApplication.MeshTyingMortarConvergenceCriteria(table)
                else:
                    Mortar = ContactStructuralMechanicsApplication.MeshTyingMortarConvergenceCriteria()
            
            Mortar.SetEchoLevel(echo_level)

            if (fancy_convergence_criterion == True):
                
                if (condn_convergence_criterion == True):
                    # Construct the solver
                    import eigen_solver_factory
                    settings_max = KratosMultiphysics.Parameters("""
                    {
                        "solver_type"             : "power_iteration_highest_eigenvalue_solver",
                        "max_iteration"           : 10000,
                        "tolerance"               : 1e-9,
                        "required_eigen_number"   : 1,
                        "verbosity"               : 0,
                        "linear_solver_settings"  : {
                            "solver_type"             : "SuperLUSolver",
                            "max_iteration"           : 500,
                            "tolerance"               : 1e-9,
                            "scaling"                 : false,
                            "verbosity"               : 0
                        }
                    }
                    """)
                    eigen_solver_max = eigen_solver_factory.ConstructSolver(settings_max)
                    settings_min = KratosMultiphysics.Parameters("""
                    {
                        "solver_type"             : "power_iteration_eigenvalue_solver",
                        "max_iteration"           : 10000,
                        "tolerance"               : 1e-9,
                        "required_eigen_number"   : 1,
                        "verbosity"               : 0,
                        "linear_solver_settings"  : {
                            "solver_type"             : "SuperLUSolver",
                            "max_iteration"           : 500,
                            "tolerance"               : 1e-9,
                            "scaling"                 : false,
                            "verbosity"               : 0
                        }
                    }
                    """)
                    eigen_solver_min = eigen_solver_factory.ConstructSolver(settings_min)
                    
                    condition_number_utility = KratosMultiphysics.ConditionNumberUtility(eigen_solver_max, eigen_solver_min)
                else:
                    condition_number_utility = None
                
                self.mechanical_convergence_criterion = ContactStructuralMechanicsApplication.MortarAndConvergenceCriteria(self.mechanical_convergence_criterion, Mortar, table, print_convergence_criterion, condition_number_utility)
            else:
                self.mechanical_convergence_criterion = ContactStructuralMechanicsApplication.MortarAndConvergenceCriteria(self.mechanical_convergence_criterion, Mortar)
            self.mechanical_convergence_criterion.SetEchoLevel(echo_level)
            self.mechanical_convergence_criterion.SetActualizeRHSFlag(True)
        
        else: # Standard criteria (same as structural mechanics application)
            # Construction of the class convergence_criterion
            import convergence_criteria_factory
            base_mechanical_convergence_criterion = convergence_criteria_factory.convergence_criterion(convergence_criterion_parameters)
        
            # Adding the mortar criteria
            if  (convergence_criterion_parameters["mortar_type"].GetString() == "ALMContactFrictionless"):
                Mortar = ContactStructuralMechanicsApplication.ALMFrictionlessMortarConvergenceCriteria()
                Mortar.SetEchoLevel(echo_level)
                self.mechanical_convergence_criterion = KratosMultiphysics.AndCriteria( base_mechanical_convergence_criterion.mechanical_convergence_criterion, Mortar)
                (self.mechanical_convergence_criterion).SetActualizeRHSFlag(True)
            elif  (convergence_criterion_parameters["mortar_type"].GetString() == "ALMContactFrictional"):
                Mortar = ContactStructuralMechanicsApplication.ALMFrictionalMortarConvergenceCriteria()
                Mortar.SetEchoLevel(echo_level)
                self.mechanical_convergence_criterion = KratosMultiphysics.AndCriteria( base_mechanical_convergence_criterion.mechanical_convergence_criterion, Mortar)
                (self.mechanical_convergence_criterion).SetActualizeRHSFlag(True)
            elif ("MeshTying" in convergence_criterion_parameters["mortar_type"].GetString()):
                Mortar = ContactStructuralMechanicsApplication.MeshTyingMortarConvergenceCriteria()
                Mortar.SetEchoLevel(echo_level)
                self.mechanical_convergence_criterion = KratosMultiphysics.AndCriteria( base_mechanical_convergence_criterion.mechanical_convergence_criterion, Mortar)
                (self.mechanical_convergence_criterion).SetActualizeRHSFlag(True)
        


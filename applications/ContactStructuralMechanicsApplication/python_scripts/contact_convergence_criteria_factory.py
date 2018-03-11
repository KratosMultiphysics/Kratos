## This script collects the available convergence criteria to be used in the ContactStructuralMechanicsApplication

from __future__ import print_function, absolute_import, division  # makes KM backward compatible with python 2.6 and 2.7
#import kratos core and applications
import KratosMultiphysics as KM
import KratosMultiphysics.StructuralMechanicsApplication as SMA
import KratosMultiphysics.ContactStructuralMechanicsApplication as CSMA

# Check that KM was imported in the main script
KM.CheckForPreviousImport()

# Convergence criteria class
class convergence_criterion:
    def __init__(self, convergence_criterion_parameters):
        # Note that all the convergence settings are introduced via a Kratos parameters object.
        echo_level = convergence_criterion_parameters["echo_level"].GetInt()
        convergence_criterion_name = convergence_criterion_parameters["convergence_criterion"].GetString()
        if "contact" in convergence_criterion_name:
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
                KM.Logger.PrintInfo("::[Mechanical Solver]:: ", "CONVERGENCE CRITERION : " + convergence_criterion_name)

            if (fancy_convergence_criterion == True):
                table = KM.TableStreamUtility()

            if(convergence_criterion_name == "contact_displacement_criterion"):
                if (fancy_convergence_criterion == True):
                    self.mechanical_convergence_criterion = CSMA.DisplacementLagrangeMultiplierContactCriteria(D_RT, D_AT, D_RT, D_AT, ensure_contact, table, print_convergence_criterion)
                else:
                    self.mechanical_convergence_criterion = CSMA.DisplacementLagrangeMultiplierContactCriteria(D_RT, D_AT, D_RT, D_AT, ensure_contact)
                self.mechanical_convergence_criterion.SetEchoLevel(echo_level)

            elif(convergence_criterion_name == "contact_residual_criterion"):
                if (fancy_convergence_criterion == True):
                    self.mechanical_convergence_criterion = CSMA.DisplacementLagrangeMultiplierResidualContactCriteria(R_RT, R_AT, CR_RT, CR_AT, ensure_contact, table, print_convergence_criterion)
                else:
                    self.mechanical_convergence_criterion = CSMA.DisplacementLagrangeMultiplierResidualContactCriteria(R_RT, R_AT, CR_RT, CR_AT, ensure_contact)
                self.mechanical_convergence_criterion.SetEchoLevel(echo_level)

            elif(convergence_criterion_name == "contact_mixed_criterion"):
                if (fancy_convergence_criterion == True):
                    self.mechanical_convergence_criterion = CSMA.DisplacementLagrangeMultiplierMixedContactCriteria(R_RT, R_AT, CR_RT, CR_AT, ensure_contact, table, print_convergence_criterion)
                else:
                    self.mechanical_convergence_criterion = CSMA.DisplacementLagrangeMultiplierMixedContactCriteria(R_RT, R_AT, CR_RT, CR_AT, ensure_contact)
                self.mechanical_convergence_criterion.SetEchoLevel(echo_level)

            elif(convergence_criterion_name == "contact_and_criterion"):
                if (fancy_convergence_criterion == True):
                    Displacement = CSMA.DisplacementLagrangeMultiplierContactCriteria(D_RT, D_AT, CD_RT, CD_AT, ensure_contact, table, print_convergence_criterion)
                    Residual = CSMA.DisplacementLagrangeMultiplierResidualContactCriteria(R_RT, R_AT, CR_RT, CR_AT, ensure_contact, table, print_convergence_criterion)
                else:
                    Displacement = CSMA.DisplacementLagrangeMultiplierContactCriteria(D_RT, D_AT, CD_RT, CD_AT, ensure_contact)
                    Residual = CSMA.DisplacementLagrangeMultiplierResidualContactCriteria(R_RT, R_AT, CR_RT, CR_AT, ensure_contact)

                Displacement.SetEchoLevel(echo_level)
                Residual.SetEchoLevel(echo_level)
                self.mechanical_convergence_criterion = KM.AndCriteria(Residual, Displacement)

            elif(convergence_criterion_name == "contact_or_criterion"):
                if (fancy_convergence_criterion == True):
                    Displacement = CSMA.DisplacementLagrangeMultiplierContactCriteria(D_RT, D_AT, CD_RT, CD_AT, ensure_contact, table, print_convergence_criterion)
                    Residual = CSMA.DisplacementLagrangeMultiplierResidualContactCriteria(R_RT, R_AT, CR_RT, CR_AT, ensure_contact, table, print_convergence_criterion)
                else:
                    Displacement = CSMA.DisplacementLagrangeMultiplierContactCriteria(D_RT, D_AT, CD_RT, CD_AT,ensure_contact)
                    Residual = CSMA.DisplacementLagrangeMultiplierResidualContactCriteria(R_RT, R_AT, CR_RT, CR_AT, ensure_contact)

                Displacement.SetEchoLevel(echo_level)
                Residual.SetEchoLevel(echo_level)
                self.mechanical_convergence_criterion = KM.OrCriteria(Residual, Displacement)

            # Adding the mortar criteria
            mortar_type = convergence_criterion_parameters["mortar_type"].GetString()
            if  (mortar_type == "ALMContactFrictionless"):
                if (fancy_convergence_criterion == True):
                    Mortar = CSMA.ALMFrictionlessMortarConvergenceCriteria(table, print_convergence_criterion, gidio_debug)
                else:
                    Mortar = CSMA.ALMFrictionlessMortarConvergenceCriteria()
            elif  (mortar_type == "ALMContactFrictionlessComponents"):
                if (fancy_convergence_criterion == True):
                    Mortar = CSMA.ALMFrictionlessComponentsMortarConvergenceCriteria(table, print_convergence_criterion, gidio_debug)
                else:
                    Mortar = CSMA.ALMFrictionlessComponentsMortarConvergenceCriteria()
            elif  (mortar_type == "ALMContactFrictional"):
                if (fancy_convergence_criterion == True):
                    Mortar = CSMA.ALMFrictionalMortarConvergenceCriteria(table, print_convergence_criterion, gidio_debug)
                else:
                    Mortar = CSMA.ALMFrictionalMortarConvergenceCriteria()
            elif ("MeshTying" in mortar_type):
                if (fancy_convergence_criterion == True):
                    Mortar = CSMA.MeshTyingMortarConvergenceCriteria(table)
                else:
                    Mortar = CSMA.MeshTyingMortarConvergenceCriteria()

            Mortar.SetEchoLevel(echo_level)

            if (fancy_convergence_criterion == True):

                if (condn_convergence_criterion == True):
                    # Construct the solver
                    import eigen_solver_factory
                    settings_max = KM.Parameters("""
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
                    settings_min = KM.Parameters("""
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

                    condition_number_utility = KM.ConditionNumberUtility(eigen_solver_max, eigen_solver_min)
                else:
                    condition_number_utility = None

                self.mechanical_convergence_criterion = CSMA.MortarAndConvergenceCriteria(self.mechanical_convergence_criterion, Mortar, table, print_convergence_criterion, condition_number_utility)
            else:
                self.mechanical_convergence_criterion = CSMA.MortarAndConvergenceCriteria(self.mechanical_convergence_criterion, Mortar)
            self.mechanical_convergence_criterion.SetEchoLevel(echo_level)
            self.mechanical_convergence_criterion.SetActualizeRHSFlag(True)

        else: # Standard criteria (same as structural mechanics application)
            # Construction of the class convergence_criterion
            import convergence_criteria_factory
            base_mechanical_convergence_criterion = convergence_criteria_factory.convergence_criterion(convergence_criterion_parameters)

            # Adding the mortar criteria
            mortar_type = convergence_criterion_parameters["mortar_type"].GetString()
            if  (mortar_type == "ALMContactFrictionless"):
                Mortar = CSMA.ALMFrictionlessMortarConvergenceCriteria()
                Mortar.SetEchoLevel(echo_level)
                self.mechanical_convergence_criterion = KM.AndCriteria( base_mechanical_convergence_criterion.mechanical_convergence_criterion, Mortar)
                (self.mechanical_convergence_criterion).SetActualizeRHSFlag(True)
            elif  (mortar_type == "ALMContactFrictionless"):
                Mortar = CSMA.ALMFrictionlessComponentsMortarConvergenceCriteria()
                Mortar.SetEchoLevel(echo_level)
                self.mechanical_convergence_criterion = KM.AndCriteria( base_mechanical_convergence_criterion.mechanical_convergence_criterion, Mortar)
                (self.mechanical_convergence_criterion).SetActualizeRHSFlag(True)
            elif  (mortar_type == "ALMContactFrictional"):
                Mortar = CSMA.ALMFrictionalMortarConvergenceCriteria()
                Mortar.SetEchoLevel(echo_level)
                self.mechanical_convergence_criterion = KM.AndCriteria( base_mechanical_convergence_criterion.mechanical_convergence_criterion, Mortar)
                (self.mechanical_convergence_criterion).SetActualizeRHSFlag(True)
            elif ("MeshTying" in mortar_type):
                Mortar = CSMA.MeshTyingMortarConvergenceCriteria()
                Mortar.SetEchoLevel(echo_level)
                self.mechanical_convergence_criterion = KM.AndCriteria( base_mechanical_convergence_criterion.mechanical_convergence_criterion, Mortar)
                (self.mechanical_convergence_criterion).SetActualizeRHSFlag(True)

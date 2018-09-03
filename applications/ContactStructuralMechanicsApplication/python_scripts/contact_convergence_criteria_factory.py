from __future__ import print_function, absolute_import, division  # makes KM backward compatible with python 2.6 and 2.7
#import kratos core and applications
import KratosMultiphysics as KM

# Check that applications were imported in the main script
KM.CheckRegisteredApplications("StructuralMechanicsApplication")
KM.CheckRegisteredApplications("ContactStructuralMechanicsApplication")

# Import applications
import KratosMultiphysics.StructuralMechanicsApplication as SMA
import KratosMultiphysics.ContactStructuralMechanicsApplication as CSMA

# Convergence criteria class
class convergence_criterion:
    def __init__(self, convergence_criterion_parameters):
        # Note that all the convergence settings are introduced via a Kratos parameters object.
        self.echo_level = convergence_criterion_parameters["echo_level"].GetInt()
        self.convergence_criterion_name = convergence_criterion_parameters["convergence_criterion"].GetString()
        self.mortar_type = convergence_criterion_parameters["mortar_type"].GetString()
        self.frictional_decomposed = convergence_criterion_parameters["frictional_decomposed"].GetBool()
        self.print_convergence_criterion = convergence_criterion_parameters["print_convergence_criterion"].GetBool()
        self.gidio_debug = convergence_criterion_parameters["gidio_debug"].GetBool()
        if "contact" in self.convergence_criterion_name:
            D_RT = convergence_criterion_parameters["displacement_relative_tolerance"].GetDouble()
            D_AT = convergence_criterion_parameters["displacement_absolute_tolerance"].GetDouble()
            R_RT = convergence_criterion_parameters["residual_relative_tolerance"].GetDouble()
            R_AT = convergence_criterion_parameters["residual_absolute_tolerance"].GetDouble()
            CD_RT = convergence_criterion_parameters["contact_displacement_relative_tolerance"].GetDouble()
            CD_AT = convergence_criterion_parameters["contact_displacement_absolute_tolerance"].GetDouble()
            CR_RT = convergence_criterion_parameters["contact_residual_relative_tolerance"].GetDouble()
            CR_AT = convergence_criterion_parameters["contact_residual_absolute_tolerance"].GetDouble()
            FCD_RT = convergence_criterion_parameters["frictional_contact_displacement_relative_tolerance"].GetDouble()
            FCD_AT = convergence_criterion_parameters["frictional_contact_displacement_absolute_tolerance"].GetDouble()
            FCR_RT = convergence_criterion_parameters["frictional_contact_residual_relative_tolerance"].GetDouble()
            FCR_AT = convergence_criterion_parameters["frictional_contact_residual_absolute_tolerance"].GetDouble()
            condn_convergence_criterion = convergence_criterion_parameters["condn_convergence_criterion"].GetBool()
            ensure_contact = convergence_criterion_parameters["ensure_contact"].GetBool()

            if(self.echo_level >= 1):
                KM.Logger.PrintInfo("::[Mechanical Solver]:: ", "CONVERGENCE CRITERION : " + self.convergence_criterion_name)

            if(self.convergence_criterion_name == "contact_displacement_criterion"):
                if (self.mortar_type == "ALMContactFrictional" and self.frictional_decomposed):
                    self.mechanical_convergence_criterion = CSMA.DisplacementLagrangeMultiplierFrictionalContactCriteria(D_RT, D_AT, CD_RT, CD_AT, FCD_RT, FCD_AT, ensure_contact, self.print_convergence_criterion)
                else:
                    self.mechanical_convergence_criterion = CSMA.DisplacementLagrangeMultiplierContactCriteria(D_RT, D_AT, CD_RT, CD_AT, ensure_contact, self.print_convergence_criterion)
                self.mechanical_convergence_criterion.SetEchoLevel(self.echo_level)

            elif(self.convergence_criterion_name == "contact_residual_criterion"):
                if (self.mortar_type == "ALMContactFrictional" and self.frictional_decomposed):
                    self.mechanical_convergence_criterion = CSMA.DisplacementLagrangeMultiplierResidualFrictionalContactCriteria(R_RT, R_AT, CR_RT, CR_AT, FCR_RT, FCR_AT, ensure_contact, self.print_convergence_criterion)
                else:
                    self.mechanical_convergence_criterion = CSMA.DisplacementLagrangeMultiplierResidualContactCriteria(R_RT, R_AT, CR_RT, CR_AT, ensure_contact, self.print_convergence_criterion)
                self.mechanical_convergence_criterion.SetEchoLevel(self.echo_level)

            elif(self.convergence_criterion_name == "contact_mixed_criterion"):
                if (self.mortar_type == "ALMContactFrictional" and self.frictional_decomposed):
                    self.mechanical_convergence_criterion = CSMA.DisplacementLagrangeMultiplierMixedFrictionalontactCriteria(R_RT, R_AT, CR_RT, CR_AT, FCR_RT, FCR_AT, ensure_contact, self.print_convergence_criterion)
                else:
                    self.mechanical_convergence_criterion = CSMA.DisplacementLagrangeMultiplierMixedContactCriteria(R_RT, R_AT, CR_RT, CR_AT, ensure_contact, self.print_convergence_criterion)
                self.mechanical_convergence_criterion.SetEchoLevel(self.echo_level)

            elif(self.convergence_criterion_name == "contact_and_criterion"):
                Displacement = CSMA.DisplacementLagrangeMultiplierContactCriteria(D_RT, D_AT, CD_RT, CD_AT, ensure_contact, self.print_convergence_criterion)
                Residual = CSMA.DisplacementLagrangeMultiplierResidualContactCriteria(R_RT, R_AT, CR_RT, CR_AT, ensure_contact, self.print_convergence_criterion)

                Displacement.SetEchoLevel(self.echo_level)
                Residual.SetEchoLevel(self.echo_level)
                self.mechanical_convergence_criterion = KM.AndCriteria(Residual, Displacement)

            elif(self.convergence_criterion_name == "contact_or_criterion"):
                Displacement = CSMA.DisplacementLagrangeMultiplierContactCriteria(D_RT, D_AT, CD_RT, CD_AT, ensure_contact, self.print_convergence_criterion)
                Residual = CSMA.DisplacementLagrangeMultiplierResidualContactCriteria(R_RT, R_AT, CR_RT, CR_AT, ensure_contact, self.print_convergence_criterion)

                Displacement.SetEchoLevel(self.echo_level)
                Residual.SetEchoLevel(self.echo_level)
                self.mechanical_convergence_criterion = KM.OrCriteria(Residual, Displacement)

            # Adding the mortar criteria

            Mortar = self.GetMortarCriteria()

            if (condn_convergence_criterion is True):
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

            self.mechanical_convergence_criterion = CSMA.MortarAndConvergenceCriteria(self.mechanical_convergence_criterion, Mortar,  self.print_convergence_criterion, condition_number_utility)

            self.mechanical_convergence_criterion.SetEchoLevel(self.echo_level)
            self.mechanical_convergence_criterion.SetActualizeRHSFlag(True)

        elif self.convergence_criterion_name == "adaptative_remesh_criteria":
            self.mechanical_convergence_criterion = None
        else: # Standard criteria (same as structural mechanics application)
            # Construction of the class convergence_criterion
            import convergence_criteria_factory
            base_mechanical_convergence_criterion = convergence_criteria_factory.convergence_criterion(convergence_criterion_parameters)

            # Adding the mortar criteria
            Mortar = self.GetMortarCriteria(False)
            if ("ALMContact" in self.mortar_type or "MeshTying" in self.mortar_type):
                self.mechanical_convergence_criterion = KM.AndCriteria( base_mechanical_convergence_criterion.mechanical_convergence_criterion, Mortar)
                (self.mechanical_convergence_criterion).SetActualizeRHSFlag(True)
            else:
                self.mechanical_convergence_criterion = base_mechanical_convergence_criterion.mechanical_convergence_criterion

    def GetMortarCriteria(self, include_table = True):
        # Adding the mortar criteria
        if (self.mortar_type == "ALMContactFrictionless"):
            if (include_table is True):
                Mortar = CSMA.ALMFrictionlessMortarConvergenceCriteria(self.print_convergence_criterion, self.gidio_debug)
            else:
                Mortar = CSMA.ALMFrictionlessMortarConvergenceCriteria()
        elif (self.mortar_type == "ALMContactFrictionlessComponents"):
            if (include_table is True):
                Mortar = CSMA.ALMFrictionlessComponentsMortarConvergenceCriteria(self.print_convergence_criterion, self.gidio_debug)
            else:
                Mortar = CSMA.ALMFrictionlessComponentsMortarConvergenceCriteria()
        elif (self.mortar_type == "ALMContactFrictional"):
            if (include_table is True):
                Mortar = CSMA.ALMFrictionalMortarConvergenceCriteria(self.print_convergence_criterion, self.gidio_debug)
            else:
                Mortar = CSMA.ALMFrictionalMortarConvergenceCriteria()
        elif ("MeshTying" in self.mortar_type):
            Mortar = CSMA.MeshTyingMortarConvergenceCriteria()

        Mortar.SetEchoLevel(self.echo_level)

        return Mortar

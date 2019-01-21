from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication

# Import base class file
import structural_mechanics_solver


def CreateSolver(model, custom_settings):
    return ExplicitMechanicalSolver(model, custom_settings)

class ExplicitMechanicalSolver(structural_mechanics_solver.MechanicalSolver):
    """The structural mechanics explicit dynamic solver.

    This class creates the mechanical solvers for explicit dynamic analysis.

    Public member variables:
    dynamic_settings -- settings for the explicit dynamic solvers.

    See structural_mechanics_solver.py for more information.
    """
    def __init__(self, model, custom_settings):
         # Set defaults and validate custom settings.
        self.dynamic_settings = KratosMultiphysics.Parameters("""
        {
            "scheme_type"                : "central_differences",
            "time_step_prediction_level" : 0,
            "max_delta_time"             : 1.0e-5,
            "fraction_delta_time"        : 0.9,
            "rayleigh_alpha"             : 0.0,
            "rayleigh_beta"              : 0.0,
            "determine_rayleigh_damping" : false,
            "determine_rayleigh_damping_settings" : {
                "echo_level"          : 0,
                "write_on_properties" : true,
                "damping_ratio_0"     : 0.0,
                "damping_ratio_1"     : -1.0,
                "eigen_system_settings" : {
                    "solver_type"                : "FEASTSolver",
                    "print_feast_output"         : false,
                    "perform_stochastic_estimate": true,
                    "solve_eigenvalue_problem"   : true,
                    "lambda_min"                 : 0.0,
                    "lambda_max"                 : 4.0e5,
                    "number_of_eigenvalues"      : 2,
                    "search_dimension"           : 15,
                    "linear_solver_settings": {
                        "solver_type": "skyline_lu"
                    }
                }
            }
        }
        """)
        self.validate_and_transfer_matching_settings(custom_settings, self.dynamic_settings)
        # Validate the remaining settings in the base class.

        # Construct the base solver.
        super(ExplicitMechanicalSolver, self).__init__(model, custom_settings)
        # Lumped mass-matrix is necessary for explicit analysis
        self.main_model_part.ProcessInfo[KratosMultiphysics.COMPUTE_LUMPED_MASS_MATRIX] = True
        self.print_on_rank_zero("::[ExplicitMechanicalSolver]:: Construction finished")

    def AddVariables(self):
        super(ExplicitMechanicalSolver, self).AddVariables()
        self._add_dynamic_variables()
        self.main_model_part.AddNodalSolutionStepVariable(StructuralMechanicsApplication.MIDDLE_VELOCITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_MASS)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.FORCE_RESIDUAL)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.RESIDUAL_VECTOR)

        if (self.settings["rotation_dofs"].GetBool()):
            self.main_model_part.AddNodalSolutionStepVariable(StructuralMechanicsApplication.MIDDLE_ANGULAR_VELOCITY)
            self.main_model_part.AddNodalSolutionStepVariable(StructuralMechanicsApplication.NODAL_INERTIA)
            self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.MOMENT_RESIDUAL)

        self.print_on_rank_zero("::[ExplicitMechanicalSolver]:: Variables ADDED")

    def AddDofs(self):
        super(ExplicitMechanicalSolver, self).AddDofs()
        self._add_dynamic_dofs()
        self.print_on_rank_zero("::[ExplicitMechanicalSolver]:: DOF's ADDED")

    def InitializeSolutionStep(self):
        # Using the base InitializeSolutionStep
        super(ExplicitMechanicalSolver, self).InitializeSolutionStep()

        # We compute for the different parts (each body can have it own damping)
        if self.dynamic_settings["determine_rayleigh_damping"].GetBool() and self.main_model_part.ProcessInfo[KratosMultiphysics.STEP] == 1:
            import KratosMultiphysics.kratos_utilities as kratos_utils
            if kratos_utils.IsApplicationAvailable("ExternalSolversApplication"):
                from KratosMultiphysics import ExternalSolversApplication
            else:
                raise Exception("ExternalSolversApplication not available")

            # The general damping ratios
            damping_ratio_0 = self.dynamic_settings["determine_rayleigh_damping_settings"]["damping_ratio_0"].GetDouble()
            damping_ratio_1 = self.dynamic_settings["determine_rayleigh_damping_settings"]["damping_ratio_1"].GetDouble()

            # We get the model parts which divide the problem
            structural_parts = self.__extract_model_parts(self.settings["problem_domain_sub_model_part_list"])
            for part in structural_parts:

                # Reseting values
                self.dynamic_settings["determine_rayleigh_damping_settings"]["damping_ratio_0"].SetDouble(damping_ratio_0)
                self.dynamic_settings["determine_rayleigh_damping_settings"]["damping_ratio_1"].SetDouble(damping_ratio_1)

                # We check if the values are previously defined
                properties = part.GetProperties()
                for prop in properties:
                    if prop.Has(StructuralMechanicsApplication.SYSTEM_DAMPING_RATIO):
                        self.dynamic_settings["determine_rayleigh_damping_settings"]["damping_ratio_0"].SetDouble(prop.GetValue(StructuralMechanicsApplication.SYSTEM_DAMPING_RATIO))
                        break
                for prop in properties:
                    if prop.Has(StructuralMechanicsApplication.SECOND_SYSTEM_DAMPING_RATIO):
                        self.dynamic_settings["determine_rayleigh_damping_settings"]["damping_ratio_1"].SetDouble(prop.GetValue(StructuralMechanicsApplication.SECOND_SYSTEM_DAMPING_RATIO))
                        break

                compute_rayleigh_damping_process = StructuralMechanicsApplication.ComputeRayleighDampingCoefficientsProcess(part, self.dynamic_settings["determine_rayleigh_damping_settings"])
                compute_rayleigh_damping_process.Execute()

            # Reseting values
            self.dynamic_settings["determine_rayleigh_damping_settings"]["damping_ratio_0"].SetDouble(damping_ratio_0)
            self.dynamic_settings["determine_rayleigh_damping_settings"]["damping_ratio_1"].SetDouble(damping_ratio_1)

    #### Specific internal functions ####
    def _create_solution_scheme(self):
        scheme_type = self.dynamic_settings["scheme_type"].GetString()

        # Setting the Rayleigh damping parameters
        self.main_model_part.ProcessInfo[StructuralMechanicsApplication.RAYLEIGH_ALPHA] = self.dynamic_settings["rayleigh_alpha"].GetDouble()
        self.main_model_part.ProcessInfo[StructuralMechanicsApplication.RAYLEIGH_BETA] = self.dynamic_settings["rayleigh_beta"].GetDouble()

        # Setting the time integration schemes
        if(scheme_type == "central_differences"):
            mechanical_scheme = StructuralMechanicsApplication.ExplicitCentralDifferencesScheme(self.dynamic_settings["max_delta_time"].GetDouble(),
                                                                             self.dynamic_settings["fraction_delta_time"].GetDouble(),
                                                                             self.dynamic_settings["time_step_prediction_level"].GetDouble())
        else:
            err_msg =  "The requested scheme type \"" + scheme_type + "\" is not available!\n"
            err_msg += "Available options are: \"central_differences\""
            raise Exception(err_msg)
        return mechanical_scheme

    def _create_mechanical_solution_strategy(self):
        computing_model_part = self.GetComputingModelPart()
        mechanical_scheme = self.get_solution_scheme()

        mechanical_solution_strategy = StructuralMechanicsApplication.MechanicalExplicitStrategy(computing_model_part,
                                            mechanical_scheme,
                                            self.settings["compute_reactions"].GetBool(),
                                            self.settings["reform_dofs_at_each_step"].GetBool(),
                                            self.settings["move_mesh_flag"].GetBool())

        mechanical_solution_strategy.SetRebuildLevel(0)
        return mechanical_solution_strategy

    #### Private functions ####

    def __extract_model_parts(self, params_list):
        """Extracting a list of SubModelParts to be added to the ComputingModelPart
        If the MainModelPart is to be added to the ComputingModelPart, then this is
        done directly, no need to also add the SubModelParts, since they are contained
        in the MainModelpart
        """
        list_model_part_names = [params_list[i].GetString() for i in range(params_list.size())]
        main_model_part_name = self.settings["model_part_name"].GetString()
        model_parts_list = []
        for model_part_name  in list_model_part_names:
            # only SubModelParts of the MainModelPart can be used!
            if model_part_name != main_model_part_name:
                full_model_part_name = main_model_part_name + "." + model_part_name
            else:
                full_model_part_name = main_model_part_name
            model_parts_list.append(self.model[full_model_part_name])
        return model_parts_list


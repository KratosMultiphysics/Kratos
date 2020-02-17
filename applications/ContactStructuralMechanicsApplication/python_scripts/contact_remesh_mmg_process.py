from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics as KratosMultiphysics
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import KratosMultiphysics.ContactStructuralMechanicsApplication as ContactStructuralMechanicsApplication
import KratosMultiphysics.MeshingApplication as MeshingApplication

# Import base process
from KratosMultiphysics.MeshingApplication.mmg_process import MmgProcess

def Factory(settings, Model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ContactRemeshMmgProcess(Model, settings["Parameters"])


class ContactRemeshMmgProcess(MmgProcess):
    """This process remeshes using MMG library. This process uses different utilities and processes. It is adapted to be used for contact problems

    Only the member variables listed below should be accessed directly.

    Public member variables:
    Model -- the container of the different model parts.
    settings -- Kratos parameters containing solver settings.
    """

    def __init__(self, Model, settings ):
        """ The default constructor of the class

        Keyword arguments:
        self -- It signifies an instance of a class.
        Model -- the container of the different model parts.
        settings -- Kratos parameters containing solver settings.
        """

        ## Settings string in json format
        default_parameters = KratosMultiphysics.Parameters("""
        {
            "help"                             : "This process remeshes using MMG library. This process uses different utilities and processes. It is adapted to be used for contact problems",
            "mesh_id"                          : 0,
            "filename"                         : "out",
            "automatic_normalization_factor"   : true,
            "consider_strain_energy"           : false,
            "model_part_name"                  : "PLEASE_SPECIFY_MODEL_PART_NAME",
            "blocking_threshold_size"          : false,
            "threshold_sizes" : {
                "minimal_size"                     : 0.1,
                "maximal_size"                     : 10.0
            },
            "strategy"                             : "Hessian",
            "discretization_type"                  : "Standard",
            "framework"                            : "Lagrangian",
            "internal_variables_parameters"        :
            {
                "allocation_size"                      : 1000,
                "bucket_size"                          : 4,
                "search_factor"                        : 2,
                "interpolation_type"                   : "LST",
                "internal_variable_interpolation_list" :[]
            },
            "hessian_strategy_parameters"              :{
                "metric_variable"                  : ["VON_MISES_STRESS","AUGMENTED_NORMAL_CONTACT_PRESSURE","STRAIN_ENERGY"],
                "non_historical_metric_variable"   : [true, true, true],
                "normalization_factor"             : [1.0, 1.0, 1.0],
                "estimate_interpolation_error"     : false,
                "interpolation_error"              : 0.04,
                "mesh_dependent_constant"          : 0.28125
            },
            "error_strategy_parameters"              :{
                "compute_error_extra_parameters":
                {
                    "stress_vector_variable"              : "CAUCHY_STRESS_VECTOR",
                    "penalty_normal"                      : 1.0e4,
                    "penalty_tangential"                  : 1.0e4
                },
                "error_metric_parameters"                 :
                {
                    "error_threshold"                       : 1.0e-4,
                    "interpolation_error"                   : 0.04
                },
                "set_target_number_of_elements"       : false,
                "target_number_of_elements"           : 1000,
                "perform_nodal_h_averaging"           : false
            },
            "enforce_current"                  : true,
            "remesh_control_type"              : "step",
            "initial_step"                     : 1,
            "step_frequency"                   : 0,
            "interval"                         : [0.0, 1e30],
            "time_stepping"                    : {},
            "automatic_remesh"                 : true,
            "automatic_remesh_parameters"      :{
                "automatic_remesh_type"            : "Ratio",
                "min_size_ratio"                   : 1.0,
                "max_size_ratio"                   : 3.0,
                "refer_type"                       : "Mean",
                "min_size_current_percentage"      : 50.0,
                "max_size_current_percentage"      : 98.0
            },
            "fix_contour_model_parts"          : [],
            "fix_conditions_model_parts"       : [],
            "fix_elements_model_parts"         : [],
            "force_min"                        : false,
            "minimal_size"                     : 0.1,
            "force_max"                        : false,
            "maximal_size"                     : 10.0,
            "advanced_parameters"                  :
            {
                "force_hausdorff_value"               : false,
                "hausdorff_value"                     : 0.0001,
                "no_move_mesh"                        : false,
                "no_surf_mesh"                        : false,
                "no_insert_mesh"                      : false,
                "no_swap_mesh"                        : false,
                "deactivate_detect_angle"             : false,
                "force_gradation_value"               : false,
                "gradation_value"                     : 1.3
            },
            "anisotropy_remeshing"                 : true,
            "enforce_anisotropy_relative_variable" : false,
            "anisotropy_parameters":{
                "reference_variable_name"          : "DISTANCE",
                "hmin_over_hmax_anisotropic_ratio" : 0.01,
                "boundary_layer_max_distance"      : 1.0,
                "boundary_layer_min_size_ratio"    : 2.0,
                "interpolation"                    : "Linear"
            },
            "collapse_prisms_elements"         : false,
            "save_external_files"              : false,
            "save_colors_files"                : false,
            "save_mdpa_file"                   : false,
            "max_number_of_searchs"            : 1000,
            "preserve_flags"                   : false,
            "interpolate_non_historical"       : true,
            "extrapolate_contour_values"       : true,
            "surface_elements"                 : false,
            "search_parameters"                : {
                "allocation_size"                  : 1000,
                "bucket_size"                      : 4,
                "search_factor"                    : 2.0
            },
            "debug_mode"                       : "",
            "debug_result_mesh"                : false,
            "echo_level"                       : 3
        }
        """)

        # Time stepping settings
        self.time_stepping = KratosMultiphysics.Parameters("""{}""")
        if settings.Has("time_stepping"):
            self.time_stepping = settings["time_stepping"].Clone()
            settings.RemoveValue("time_stepping")

        # Validate the settings
        settings.RecursivelyValidateAndAssignDefaults(default_parameters)

        # Parameters of the process
        self.automatic_normalization_factor = settings["automatic_normalization_factor"].GetBool()
        self.consider_strain_energy = settings["consider_strain_energy"].GetBool()

        # Refill missing
        number_of_metric_variable = settings["hessian_strategy_parameters"]["metric_variable"].size()
        number_of_non_historical_metric_variable = settings["hessian_strategy_parameters"]["non_historical_metric_variable"].size()
        number_of_normalization_factor = settings["hessian_strategy_parameters"]["normalization_factor"].size()
        if number_of_non_historical_metric_variable < number_of_metric_variable:
            for i in range(number_of_non_historical_metric_variable, number_of_metric_variable):
                settings["hessian_strategy_parameters"]["non_historical_metric_variable"].Append(True)
        if number_of_normalization_factor < number_of_metric_variable:
            for i in range(number_of_normalization_factor, number_of_metric_variable):
                settings["hessian_strategy_parameters"]["normalization_factor"].Append(1.0)

        # Remove unused
        if not self.consider_strain_energy:
            # Locating index to remove
            index_to_remove = None
            for i in range(number_of_metric_variable):
                if settings["hessian_strategy_parameters"]["metric_variable"][i].GetString() == "STRAIN_ENERGY":
                    index_to_remove = i

            # Copying
            auxiliar_parameters = KratosMultiphysics.Parameters("""{"metric_variable" : [], "non_historical_metric_variable" : [], "normalization_factor" : []}""")
            for i in range(number_of_metric_variable):
                if not i == index_to_remove:
                    auxiliar_parameters["metric_variable"].Append(settings["hessian_strategy_parameters"]["metric_variable"][i])
                    auxiliar_parameters["non_historical_metric_variable"].Append(settings["hessian_strategy_parameters"]["non_historical_metric_variable"][i])
                    auxiliar_parameters["normalization_factor"].Append(settings["hessian_strategy_parameters"]["normalization_factor"][i])

            # Removing old
            settings["hessian_strategy_parameters"].RemoveValue("metric_variable")
            settings["hessian_strategy_parameters"].RemoveValue("non_historical_metric_variable")
            settings["hessian_strategy_parameters"].RemoveValue("normalization_factor")

            # Adding new
            settings["hessian_strategy_parameters"].AddValue("metric_variable", auxiliar_parameters["metric_variable"])
            settings["hessian_strategy_parameters"].AddValue("non_historical_metric_variable", auxiliar_parameters["non_historical_metric_variable"])
            settings["hessian_strategy_parameters"].AddValue("normalization_factor", auxiliar_parameters["normalization_factor"])

        # Auxiliar dictionary with the variables and index
        self.variables_dict = {}
        number_of_metric_variable = settings["hessian_strategy_parameters"]["metric_variable"].size()
        for i in range(number_of_metric_variable):
            self.variables_dict[settings["hessian_strategy_parameters"]["metric_variable"][i].GetString()] = i

        # Avoid conflict with mother class
        settings.RemoveValue("automatic_normalization_factor")
        settings.RemoveValue("consider_strain_energy")

        # Construct the base process.
        super(ContactRemeshMmgProcess, self).__init__(Model, settings)

        # Create model parts
        model_part_name = settings["model_part_name"].GetString()
        self.main_model_part = Model[model_part_name]

        # Create extrapolation process
        extrapolation_parameters = KratosMultiphysics.Parameters("""
        {
            "model_part_name"            : "",
            "echo_level"                 : 0,
            "average_variable"           : "NODAL_AREA",
            "area_average"               : true,
            "list_of_variables"          : [],
            "extrapolate_non_historical" : true
        }
        """)
        if "VON_MISES_STRESS" in self.variables_dict.keys():
            extrapolation_parameters["list_of_variables"].Append("VON_MISES_STRESS")
        if self.consider_strain_energy and "STRAIN_ENERGY" in self.variables_dict.keys():
            extrapolation_parameters["list_of_variables"].Append("STRAIN_ENERGY")
        extrapolation_parameters["model_part_name"].SetString(model_part_name)
        self.integration_values_extrapolation_to_nodes_process = KratosMultiphysics.IntegrationValuesExtrapolationToNodesProcess(self.main_model_part, extrapolation_parameters)

    def ExecuteInitialize(self):
        """ This method is executed at the begining to initialize the process

        Keyword arguments:
        self -- It signifies an instance of a class.
        """

        # Calculation automatically the normalization factors
        if self.automatic_normalization_factor:
            E = 0.0
            mu = 0.0
            for prop in self.main_model_part.GetProperties():
                if prop.Has(KratosMultiphysics.YOUNG_MODULUS):
                    E = prop.GetValue(KratosMultiphysics.YOUNG_MODULUS)
                    break
            for prop in self.main_model_part.GetProperties():
                if prop.Has(KratosMultiphysics.POISSON_RATIO):
                    mu = prop.GetValue(KratosMultiphysics.POISSON_RATIO)
                    break

            normalization_factor = 2.0e1/(mu**2 * E)
            if self.consider_strain_energy and "STRAIN_ENERGY" in self.variables_dict.keys():
                self.settings["hessian_strategy_parameters"]["normalization_factor"][self.variables_dict["STRAIN_ENERGY"]].SetDouble(normalization_factor)
            if "VON_MISES_STRESS" in self.variables_dict.keys():
                self.settings["hessian_strategy_parameters"]["normalization_factor"][self.variables_dict["VON_MISES_STRESS"]].SetDouble(normalization_factor)
            if "AUGMENTED_NORMAL_CONTACT_PRESSURE" in self.variables_dict.keys():
                self.settings["hessian_strategy_parameters"]["normalization_factor"][self.variables_dict["AUGMENTED_NORMAL_CONTACT_PRESSURE"]].SetDouble(normalization_factor)

        # We print the parameters considered
        if self.consider_strain_energy and "STRAIN_ENERGY" in self.variables_dict.keys():
            KratosMultiphysics.Logger.PrintInfo("STRAIN_ENERGY NORMALIZARION FACTOR: ", "{:.2e}".format(self.settings["hessian_strategy_parameters"]["normalization_factor"][self.variables_dict["STRAIN_ENERGY"]].GetDouble()))
        if "VON_MISES_STRESS" in self.variables_dict.keys():
            KratosMultiphysics.Logger.PrintInfo("VON_MISES_STRESS NORMALIZARION FACTOR: ", "{:.2e}".format(self.settings["hessian_strategy_parameters"]["normalization_factor"][self.variables_dict["VON_MISES_STRESS"]].GetDouble()))
        if "AUGMENTED_NORMAL_CONTACT_PRESSURE" in self.variables_dict.keys():
            KratosMultiphysics.Logger.PrintInfo("AUGMENTED_NORMAL_CONTACT_PRESSURE NORMALIZARION FACTOR: ", "{:.2e}".format(self.settings["hessian_strategy_parameters"]["normalization_factor"][self.variables_dict["AUGMENTED_NORMAL_CONTACT_PRESSURE"]].GetDouble()))

        # We initialize the STRAIN_ENERGY, VON_MISES_STRESS and AUGMENTED_NORMAL_CONTACT_PRESSURE
        if self.consider_strain_energy and "STRAIN_ENERGY" in self.variables_dict.keys():
            KratosMultiphysics.VariableUtils().SetNonHistoricalVariableToZero(KratosMultiphysics.STRAIN_ENERGY, self.main_model_part.Nodes)
        if "VON_MISES_STRESS" in self.variables_dict.keys():
            KratosMultiphysics.VariableUtils().SetNonHistoricalVariableToZero(StructuralMechanicsApplication.VON_MISES_STRESS, self.main_model_part.Nodes)
        if "AUGMENTED_NORMAL_CONTACT_PRESSURE" in self.variables_dict.keys():
            KratosMultiphysics.VariableUtils().SetNonHistoricalVariableToZero(ContactStructuralMechanicsApplication.AUGMENTED_NORMAL_CONTACT_PRESSURE, self.main_model_part.Nodes)

        # We call to the base process
        super(ContactRemeshMmgProcess, self).ExecuteInitialize()

    def ExecuteBeforeSolutionLoop(self):
        """ This method is executed before starting the time loop

        Keyword arguments:
        self -- It signifies an instance of a class.
        """

        # We call the extrapolation
        self.integration_values_extrapolation_to_nodes_process.ExecuteBeforeSolutionLoop()
        self.integration_values_extrapolation_to_nodes_process.ExecuteFinalizeSolutionStep()

        # We call to the base process
        super(ContactRemeshMmgProcess, self).ExecuteBeforeSolutionLoop()

    def ExecuteInitializeSolutionStep(self):
        """ This method is executed in order to initialize the current step

        Keyword arguments:
        self -- It signifies an instance of a class.
        """

        # We call the extrapolation
        self.integration_values_extrapolation_to_nodes_process.ExecuteBeforeSolutionLoop()
        self.integration_values_extrapolation_to_nodes_process.ExecuteFinalizeSolutionStep()

        # We call to the base process
        super(ContactRemeshMmgProcess, self).ExecuteInitializeSolutionStep()

    def ExecuteFinalizeSolutionStep(self):
        """ This method is executed in order to finalize the current step

        Keyword arguments:
        self -- It signifies an instance of a class.
        """

        # We call to the base process
        super(ContactRemeshMmgProcess, self).ExecuteFinalizeSolutionStep()

    def ExecuteBeforeOutputStep(self):
        """ This method is executed right before the ouput process computation

        Keyword arguments:
        self -- It signifies an instance of a class.
        """

        # We call to the base process
        super(ContactRemeshMmgProcess, self).ExecuteBeforeOutputStep()

    def ExecuteAfterOutputStep(self):
        """ This method is executed right after the ouput process computation

        Keyword arguments:
        self -- It signifies an instance of a class.
        """

        # We call to the base process
        super(ContactRemeshMmgProcess, self).ExecuteAfterOutputStep()

    def ExecuteFinalize(self):
        """ This method is executed in order to finalize the current computation

        Keyword arguments:
        self -- It signifies an instance of a class.
        """

        ## We call the extrapolation
        #self.integration_values_extrapolation_to_nodes_process.ExecuteFinalize()

        # We call to the base process
        super(ContactRemeshMmgProcess, self).ExecuteFinalize()

    def _AuxiliarCallsBeforeRemesh(self):
        """ This method is executed right before execute the remesh

        Keyword arguments:
        self -- It signifies an instance of a class.
        """

        # We remove the submodelpart
        KratosMultiphysics.VariableUtils().SetFlag(KratosMultiphysics.TO_ERASE, True, self.main_model_part.GetSubModelPart("ComputingContact").Conditions)
        self.main_model_part.GetRootModelPart().RemoveConditionsFromAllLevels(KratosMultiphysics.TO_ERASE)

        # We clean the computing before remesh
        KratosMultiphysics.VariableUtils().SetFlag(KratosMultiphysics.TO_ERASE, True, self.main_model_part.Nodes)
        KratosMultiphysics.VariableUtils().SetFlag(KratosMultiphysics.TO_ERASE, True, self.main_model_part.Conditions)
        KratosMultiphysics.VariableUtils().SetFlag(KratosMultiphysics.TO_ERASE, True, self.main_model_part.Elements)
        KratosMultiphysics.VariableUtils().SetFlag(KratosMultiphysics.TO_ERASE, True, self.main_model_part.MasterSlaveConstraints)

        self.main_model_part.RemoveNodes(KratosMultiphysics.TO_ERASE)
        self.main_model_part.RemoveConditions(KratosMultiphysics.TO_ERASE)
        self.main_model_part.RemoveElements(KratosMultiphysics.TO_ERASE)
        self.main_model_part.RemoveMasterSlaveConstraints(KratosMultiphysics.TO_ERASE)

        # We remove the contact submodelparts
        self.main_model_part.RemoveSubModelPart("Contact")
        self.main_model_part.RemoveSubModelPart("ComputingContact")

        # Ensure properties defined
        MeshingApplication.MeshingUtilities.RecursiveEnsureModelPartOwnsProperties(self.main_model_part.GetRootModelPart())

        # We create the contact submodelparts
        self.main_model_part.CreateSubModelPart("Contact")
        self.main_model_part.CreateSubModelPart("ComputingContact")

    def _AuxiliarCallsAfterRemesh(self):
        """ This method is executed right after execute the remesh

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        KratosMultiphysics.FastTransferBetweenModelPartsProcess(self.main_model_part, self.main_model_part.GetParentModelPart()).Execute()

    def _GenerateErrorProcess(self):
        """ This method creates an erro process to compute the metric

        Keyword arguments:
        self -- It signifies an instance of a class.
        """

        # We compute the error
        error_compute_parameters = KratosMultiphysics.Parameters("""{}""")
        error_compute_parameters.AddValue("stress_vector_variable", self.settings["error_strategy_parameters"]["compute_error_extra_parameters"]["stress_vector_variable"])
        error_compute_parameters.AddValue("penalty_normal", self.settings["error_strategy_parameters"]["compute_error_extra_parameters"]["penalty_normal"])
        error_compute_parameters.AddValue("penalty_tangential", self.settings["error_strategy_parameters"]["compute_error_extra_parameters"]["penalty_tangential"])
        error_compute_parameters.AddValue("echo_level", self.settings["echo_level"])
        if self.domain_size == 2:
            return ContactStructuralMechanicsApplication.ContactSPRErrorProcess2D(self.main_model_part, error_compute_parameters)
        else:
            return ContactStructuralMechanicsApplication.ContactSPRErrorProcess3D(self.main_model_part, error_compute_parameters)

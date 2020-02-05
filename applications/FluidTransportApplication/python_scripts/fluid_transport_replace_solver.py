from __future__ import print_function, absolute_import, division # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics
from KratosMultiphysics.FluidTransportApplication.fluid_transport_solver import FluidTransportSolver

# Import applications

def CreateSolver(main_model_part, custom_settings):

    return FluidTransportReplaceSolver(main_model_part, custom_settings)


class FluidTransportReplaceSolver(FluidTransportSolver):

    def __init__(self, model, custom_settings):

        # TODO check that we can derive from FluidTransportSolver instead of ConvectionDiffusionBaseSolver
        super(FluidTransportReplaceSolver,self).__init__(model, custom_settings)

    @classmethod
    def GetDefaultSettings(cls):
        this_defaults = KratosMultiphysics.Parameters("""{
            "solver_type": "fluid_transport_solver",
            "model_part_name": "FluidTransportDomain",
            "domain_size": 2,
            "start_time": 0.0,
            "time_step": 0.1,
            "model_import_settings":{
                "input_type": "mdpa",
                "input_filename": "unknown_name",
                "input_file_label": 0
            },
            "computing_model_part_name" : "fluid_transport_computing_domain",
            "buffer_size":                        3,
            "echo_level":                         0,
            "clear_storage":                      false,
            "compute_reactions":                  false,
            "move_mesh_flag":                     false,
            "reform_dofs_at_each_step":           false,
            "block_builder":                      true,
            "solution_type":                      "Steady",
            "scheme_type":                        "Implicit",
            "newmark_theta":                      0.5,
            "strategy_type":                      "Linear",
            "convergence_criterion":              "And_criterion",
            "displacement_relative_tolerance":    1.0E-4,
            "displacement_absolute_tolerance":    1.0E-9,
            "residual_relative_tolerance":        1.0E-4,
            "residual_absolute_tolerance":        1.0E-9,
            "max_iteration":                      15,
            "linear_solver_settings":             {
                "solver_type":   "ExternalSolversApplication.super_lu",
                "tolerance": 1.0e-6,
                "max_iteration": 100,
                "scaling": false,
                "verbosity": 0,
                "preconditioner_type": "ilu0",
                "smoother_type": "ilu0",
                "krylov_type": "gmres",
                "coarsening_type": "aggregation"
            },
            "element_replace_settings" : {
                "element_name" : "EulerianConvDiff",
                "condition_name" : "FluxCondition"
            },
            "material_import_settings": {
                "materials_filename": "BuoyancyMaterials.json"
            },
            "problem_domain_sub_model_part_list": [""],
            "processes_sub_model_part_list": [""]
        }""")

        this_defaults.AddMissingParameters(super(FluidTransportReplaceSolver, cls).GetDefaultSettings())
        return this_defaults

    def PrepareModelPart(self):

        # Set ProcessInfo variables
        # self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.TIME,
        #                                           self.settings["start_time"].GetDouble())
        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DELTA_TIME,
                                                  self.settings["time_step"].GetDouble())
        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.STEP, 0)

        if not self.main_model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED]:
            # Check and prepare computing model part and import constitutive laws.
            self._execute_after_reading()
            # self._ExecuteCheckAndPrepare()

            throw_errors = False
            KratosMultiphysics.TetrahedralMeshOrientationCheck(self.main_model_part, throw_errors).Execute()

            KratosMultiphysics.ReplaceElementsAndConditionsProcess(self.main_model_part,self._get_element_condition_replace_settings()).Execute()

            # self._set_and_fill_buffer()
            self._SetBufferSize()

        if (self.settings["echo_level"].GetInt() > 0):
            KratosMultiphysics.Logger.PrintInfo(self.model)

        KratosMultiphysics.Logger.PrintInfo("FluidTransportReplaceSolver", "Model reading finished.")

    def import_materials(self):
        materials_filename = self.settings["material_import_settings"]["materials_filename"].GetString()
        if (materials_filename != ""):
            # Add constitutive laws and material properties from json file to model parts.
            material_settings = KratosMultiphysics.Parameters("""{"Parameters": {"materials_filename": ""}} """)
            material_settings["Parameters"]["materials_filename"].SetString(materials_filename)
            KratosMultiphysics.ReadMaterialsUtility(material_settings, self.model)

            # We set the properties that are nodal
            self._assign_nodally_properties()
            materials_imported = True
        else:
            materials_imported = False
        return materials_imported

    def GetComputingModelPart(self):
        return self.main_model_part.GetSubModelPart(self.settings["computing_model_part_name"].GetString())

    #### Specific internal functions ####

    def _execute_after_reading(self):
        """Prepare computing model part and import constitutive laws. """
        self.computing_model_part_name = "fluid_transport_computing_domain"

        # Auxiliary parameters object for the CheckAndPepareModelProcess
        params = KratosMultiphysics.Parameters("{}")
        params.AddValue("computing_model_part_name",self.settings["computing_model_part_name"])
        params.AddValue("problem_domain_sub_model_part_list",self.settings["problem_domain_sub_model_part_list"])
        params.AddValue("processes_sub_model_part_list",self.settings["processes_sub_model_part_list"])
        # Assign mesh entities from domain and process sub model parts to the computing model part.
        from KratosMultiphysics.ConvectionDiffusionApplication import check_and_prepare_model_process_convection_diffusion as check_and_prepare_model_process
        check_and_prepare_model_process.CheckAndPrepareModelProcess(self.main_model_part, params).Execute()

        # Import constitutive laws.
        materials_imported = self.import_materials()
        if materials_imported:
            KratosMultiphysics.Logger.PrintInfo("::[ConvectionDiffusionBaseSolver]:: ", "Materials were successfully imported.")

        else:
            KratosMultiphysics.Logger.PrintInfo("::[ConvectionDiffusionBaseSolver]:: ", "Materials were not imported.")

    def _assign_nodally_properties(self):

        # We transfer the values of the con.diff variables to the nodes
        with open(self.settings["material_import_settings"]["materials_filename"].GetString(), 'r') as parameter_file:
            materials = KratosMultiphysics.Parameters(parameter_file.read())

        for i in range(materials["properties"].size()):
            model_part = self.model.GetModelPart(materials["properties"][i]["model_part_name"].GetString())
            mat = materials["properties"][i]["Material"]

            for key, value in mat["Variables"].items():
                var = KratosMultiphysics.KratosGlobals.GetVariable(key)
                if (self._check_variable_to_set(var)):
                    if value.IsDouble():
                        KratosMultiphysics.VariableUtils().SetScalarVar(var, value.GetDouble(), model_part.Nodes)
                    elif value.IsVector():
                        KratosMultiphysics.VariableUtils().SetVectorVar(var, value.GetVector(), model_part.Nodes)
                    else:
                        raise ValueError("Type of value is not available")

    def _check_variable_to_set(self, var):
        thermal_settings = self.main_model_part.ProcessInfo[KratosMultiphysics.CONVECTION_DIFFUSION_SETTINGS]
        if (thermal_settings.IsDefinedDensityVariable()):
            if (thermal_settings.GetDensityVariable() == var):
                return True
        if (thermal_settings.IsDefinedDiffusionVariable()):
            if (thermal_settings.GetDiffusionVariable() == var):
                return True
        if (thermal_settings.IsDefinedVolumeSourceVariable()):
            if (thermal_settings.GetVolumeSourceVariable() == var):
                return True
        if (thermal_settings.IsDefinedSurfaceSourceVariable()):
            if (thermal_settings.GetSurfaceSourceVariable() == var):
                return True
        if (thermal_settings.IsDefinedProjectionVariable()):
            if (thermal_settings.GetProjectionVariable() == var):
                return True
        if (thermal_settings.IsDefinedConvectionVariable()):
            if (thermal_settings.GetConvectionVariable() == var):
                return True
        if (thermal_settings.IsDefinedTransferCoefficientVariable()):
            if (thermal_settings.GetTransferCoefficientVariable() == var):
                return True
        if (thermal_settings.IsDefinedSpecificHeatVariable()):
            if (thermal_settings.GetSpecificHeatVariable() == var):
                return True
        else:
            return False

    def _get_element_condition_replace_settings(self):
        domain_size = self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]

        ## Elements
        num_nodes_elements = 0
        if (len(self.main_model_part.Elements) > 0):
            for elem in self.main_model_part.Elements:
                num_nodes_elements = len(elem.GetNodes())
                break

        if domain_size == 2:
            if (self.settings["element_replace_settings"]["element_name"].GetString() == "EulerianConvDiff"):
                if (num_nodes_elements == 3):
                    if(self.settings["solution_type"].GetString() == "Steady"):
                        self.settings["element_replace_settings"]["element_name"].SetString("SteadyConvectionDiffusionFICElement2D3N")
                    elif(self.settings["scheme_type"].GetString() == "Implicit"):
                        self.settings["element_replace_settings"]["element_name"].SetString("TransientConvectionDiffusionFICElement2D3N")
                    else:
                        self.settings["element_replace_settings"]["element_name"].SetString("TransientConvectionDiffusionFICExplicitElement2D3N")
                else:
                    if(self.settings["solution_type"].GetString() == "Steady"):
                        self.settings["element_replace_settings"]["element_name"].SetString("SteadyConvectionDiffusionFICElement2D4N")
                    elif(self.settings["scheme_type"].GetString() == "Implicit"):
                        self.settings["element_replace_settings"]["element_name"].SetString("TransientConvectionDiffusionFICElement2D4N")
                    else:
                        self.settings["element_replace_settings"]["element_name"].SetString("TransientConvectionDiffusionFICExplicitElement2D4N")
        elif domain_size == 3:
            if (self.settings["element_replace_settings"]["element_name"].GetString() == "EulerianConvDiff"):
                if (num_nodes_elements == 4):
                    if(self.settings["solution_type"].GetString() == "Steady"):
                        self.settings["element_replace_settings"]["element_name"].SetString("SteadyConvectionDiffusionFICElement3D4N")
                    elif(self.settings["scheme_type"].GetString() == "Implicit"):
                        self.settings["element_replace_settings"]["element_name"].SetString("TransientConvectionDiffusionFICElement3D4N")
                    else:
                        self.settings["element_replace_settings"]["element_name"].SetString("TransientConvectionDiffusionFICExplicitElement3D4N")
                else:
                    if(self.settings["solution_type"].GetString() == "Steady"):
                        self.settings["element_replace_settings"]["element_name"].SetString("SteadyConvectionDiffusionFICElement3D8N")
                    elif(self.settings["scheme_type"].GetString() == "Implicit"):
                        self.settings["element_replace_settings"]["element_name"].SetString("TransientConvectionDiffusionFICElement3D8N")
                    else:
                        self.settings["element_replace_settings"]["element_name"].SetString("TransientConvectionDiffusionFICExplicitElement3D8N")
        else:
            raise Exception("DOMAIN_SIZE not set")

        ## Conditions
        num_nodes_conditions = 0
        if (len(self.main_model_part.Conditions) > 0):
            for cond in self.main_model_part.Conditions:
                num_nodes_conditions = len(cond.GetNodes())
                break

        if domain_size == 2:
            if (self.settings["element_replace_settings"]["condition_name"].GetString() == "FluxCondition"):
                self.settings["element_replace_settings"]["condition_name"].SetString("FluxCondition2D2N")
        elif domain_size == 3:
            aux_str = "3D" + str(num_nodes_conditions) + "N"
            if (self.settings["element_replace_settings"]["condition_name"].GetString() == "FluxCondition"):
                self.settings["element_replace_settings"]["condition_name"].SetString("FluxCondition" + aux_str)
        else:
            raise Exception("DOMAIN_SIZE not set")

        return self.settings["element_replace_settings"]

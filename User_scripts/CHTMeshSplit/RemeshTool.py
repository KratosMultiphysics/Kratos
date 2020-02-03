from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import os
import KratosMultiphysics
import KratosMultiphysics.MeshingApplication as KratosMesh
import KratosMultiphysics.FluidDynamicsApplication as KratosFluid
import KratosMultiphysics.ConvectionDiffusionApplication as KratosConvDiff
from KratosMultiphysics.gid_output_process import GiDOutputProcess

class CHTWorkflow():

    def __init__(self, model1, model2, strategy, import_settings=KratosMultiphysics.Parameters("""{}""")):
        
        file_path = os.path.dirname(os.path.realpath(__file__))
        modelA = KratosMultiphysics.Model()
        self.model_fluid = modelA.CreateModelPart("FluidPart")
        self.model_fluid.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, 3)
        self.model_fluid.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)
        self.model_fluid.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE_GRADIENT)
        KratosMultiphysics.ModelPartIO(file_path + model1).ReadModelPart(self.model_fluid)
        print(self.model_fluid)

        modelB = KratosMultiphysics.Model()
        self.model_solid = modelB.CreateModelPart("SolidPart")
        self.model_solid.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, 3)
        self.model_solid.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)
        self.model_solid.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE_GRADIENT)
        KratosMultiphysics.ModelPartIO(file_path + model2).ReadModelPart(self.model_solid)
        print(self.model_solid)

        default_params = KratosMultiphysics.Parameters( """
        {
            "shape_functions"                       : "standard",
            "reform_model_part_at_each_time_step"   : false,
            "visualization_variables"               : ["DISTANCE","DISTANCE_GRADIENT"],
            "help"                                  : "This process detects the skin from a given submodelpart and it generates the correspong conditions",
            "model_part_name"                       : "MainModelPart",
            "computing_model_part_name"             : "ComputingDomain",
            "recursive_detection"                   : true,
            "name_auxiliar_model_part"              : "SolidSkinPart",
            "name_auxiliar_condition"               : "Condition",
            "list_model_parts_to_assign_conditions" : [],
            "echo_level"                            : 0,
            "output_configuration"                  : {
                "result_file_configuration" : {
                    "gidpost_flags" : {
                        "GiDPostMode"           : "GiD_PostBinary",
                        "WriteDeformedMeshFlag" : "WriteDeformed",
                        "WriteConditionsFlag"   : "WriteConditions",
                        "MultiFileFlag"         : "SingleFile"
                    },
                    "file_label"          : "time",
                    "output_control_type" : "time",
                    "output_frequency"    : 0.1,
                    "body_output"         : true,
                    "node_output"         : false,
                    "skin_output"         : false,
                    "nodal_results"       : []
                },
                "point_data_configuration"  : []
            }
        } """ )

        import_settings.ValidateAndAssignDefaults(default_params)
        self.settings = import_settings
        self._refine_before_cut = strategy
        self.Solid_Refined = False

    """
    Internal Function : No call from outside the scope of the class
    """
    def _detect_solid_skin(self):
        
        self.recursive_detection = self.settings["recursive_detection"].GetBool()
        self.name_auxiliar_model_part = self.settings["name_auxiliar_model_part"].GetString()
        
        # Process parameters
        detect_skin_parameters = KratosMultiphysics.Parameters("""{}""")
        detect_skin_parameters.AddValue("name_auxiliar_model_part", self.settings["name_auxiliar_model_part"])
        detect_skin_parameters.AddValue("name_auxiliar_condition", self.settings["name_auxiliar_condition"])
        detect_skin_parameters.AddValue("list_model_parts_to_assign_conditions", self.settings["list_model_parts_to_assign_conditions"])
        detect_skin_parameters.AddValue("echo_level", self.settings["echo_level"])

        if (self.model_solid.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 2):
            self.detect_skin = KratosMultiphysics.SkinDetectionProcess2D(self.model_solid, detect_skin_parameters)
        else:
            self.detect_skin = KratosMultiphysics.SkinDetectionProcess3D(self.model_solid, detect_skin_parameters)
        self.detect_skin.Execute()

    def _fluid_distance_to_skin(self):

        #self._detect_solid_skin()
        skin_model_part = self.model_solid.GetSubModelPart(self.name_auxiliar_model_part)
        KratosMultiphysics.CalculateDistanceToSkinProcess3D(self.model_fluid, skin_model_part).Execute()
    
    def _solid_distance_to_skin(self):

        #self._detect_solid_skin()
        skin_model_part = self.model_solid.GetSubModelPart(self.name_auxiliar_model_part)
        KratosMultiphysics.CalculateDistanceToSkinProcess3D(self.model_solid, skin_model_part).Execute()

    """
    External Functions : Called from the Workflow
    """

    def MoveSolidinFluidDomain(self, dist):

        for node in self.model_solid.Nodes:
            node.X += dist[0]
            node.Y += dist[1]
            node.Z += dist[2]

    def RefineSolid(self, minimum_size):

        self.Solid_Refined = True # Set Solid Refinement Check to true

        self._detect_solid_skin() # Detect Solid Skin
        self._solid_distance_to_skin() # Calculate Distance to Skin for LevelSet

        find_nodal_h = KratosMultiphysics.FindNodalHNonHistoricalProcess(self.model_solid)
        find_nodal_h.Execute()
        
        KratosMultiphysics.VariableUtils().SetNonHistoricalVariable(KratosMultiphysics.NODAL_AREA, 0.0, self.model_solid.Nodes)
        local_gradient = KratosMultiphysics.ComputeNodalGradientProcess3D(self.model_solid, KratosMultiphysics.DISTANCE, KratosMultiphysics.DISTANCE_GRADIENT, KratosMultiphysics.NODAL_AREA)
        local_gradient.Execute()

        # We set to zero (or unit) the metric
        ZeroVector = KratosMultiphysics.Vector(6)
        ZeroVector[0] = 0.0; ZeroVector[1] = 0.0; ZeroVector[2] = 0.0
        ZeroVector[3] = 0.0; ZeroVector[4] = 0.0; ZeroVector[5] = 0.0

        for node in self.model_solid.Nodes:
        	node.SetValue(KratosMesh.METRIC_TENSOR_3D, ZeroVector)

        min_size = minimum_size
        max_dist = 1.25 * minimum_size
        # We define a metric using the ComputeLevelSetSolMetricProcess
        level_set_param = KratosMultiphysics.Parameters("""
        	{
        		"minimal_size"                         : """ + str(min_size) + """,
        		"enforce_current"                      : true,
        		"anisotropy_remeshing"                 : true,
        		"anisotropy_parameters": {
        			"hmin_over_hmax_anisotropic_ratio"      : 0.9,
        			"boundary_layer_max_distance"           : """ + str(max_dist) + """,
        			"interpolation"                         : "Linear" }
        	}
        	""")
        metric_process = KratosMesh.ComputeLevelSetSolMetricProcess3D(self.model_solid, KratosMultiphysics.DISTANCE_GRADIENT, level_set_param)
        metric_process.Execute()

        remesh_param = KratosMultiphysics.Parameters("""{ }""")
        MmgProcess = KratosMesh.MmgProcess3D(self.model_solid, remesh_param)
        MmgProcess.Execute()

        #Recompute Distance and Distance Gradient

        self._solid_distance_to_skin() # Calculate Distance to Skin
        find_nodal_h = KratosMultiphysics.FindNodalHNonHistoricalProcess(self.model_solid)
        find_nodal_h.Execute()
        
        KratosMultiphysics.VariableUtils().SetNonHistoricalVariable(KratosMultiphysics.NODAL_AREA, 0.0, self.model_solid.Nodes)
        local_gradient = KratosMultiphysics.ComputeNodalGradientProcess3D(self.model_solid, KratosMultiphysics.DISTANCE, KratosMultiphysics.DISTANCE_GRADIENT, KratosMultiphysics.NODAL_AREA)
        local_gradient.Execute()

    def RefineBeforeCut(self, level_set_size):

        if self._refine_before_cut:
            self._detect_solid_skin() # Detect Solid Skin
            self._fluid_distance_to_skin() # Calculate Distance to Skin
        
        find_nodal_h = KratosMultiphysics.FindNodalHNonHistoricalProcess(self.model_fluid)
        find_nodal_h.Execute()
        
        KratosMultiphysics.VariableUtils().SetNonHistoricalVariable(KratosMultiphysics.NODAL_AREA, 0.0, self.model_fluid.Nodes)
        local_gradient = KratosMultiphysics.ComputeNodalGradientProcess3D(self.model_fluid, KratosMultiphysics.DISTANCE, KratosMultiphysics.DISTANCE_GRADIENT, KratosMultiphysics.NODAL_AREA)
        local_gradient.Execute()

        # We set to zero (or unit) the metric
        ZeroVector = KratosMultiphysics.Vector(6)
        ZeroVector[0] = 0.0; ZeroVector[1] = 0.0; ZeroVector[2] = 0.0
        ZeroVector[3] = 0.0; ZeroVector[4] = 0.0; ZeroVector[5] = 0.0

        for node in self.model_fluid.Nodes:
        	node.SetValue(KratosMesh.METRIC_TENSOR_3D, ZeroVector)

        min_size = level_set_size
        max_dist = 1.25 * level_set_size
        # We define a metric using the ComputeLevelSetSolMetricProcess
        level_set_param = KratosMultiphysics.Parameters("""
        	{
        		"minimal_size"                         : """ + str(min_size) + """,
        		"enforce_current"                      : true,
        		"anisotropy_remeshing"                 : true,
        		"anisotropy_parameters": {
        			"hmin_over_hmax_anisotropic_ratio"      : 0.9,
        			"boundary_layer_max_distance"           : """ + str(max_dist) + """,
        			"interpolation"                         : "Linear" }
        	}
        	""")
        metric_process = KratosMesh.ComputeLevelSetSolMetricProcess3D(self.model_fluid, KratosMultiphysics.DISTANCE_GRADIENT, level_set_param)
        metric_process.Execute()

        remesh_param = KratosMultiphysics.Parameters("""{ }""")
        MmgProcess = KratosMesh.MmgProcess3D(self.model_fluid, remesh_param)
        MmgProcess.Execute()

        if self._refine_before_cut:
            self._fluid_distance_to_skin() # Recalculate distance to skin


    def DistancetoSolidSkin(self):

        self._detect_solid_skin() # Detect Solid Skin
        self._fluid_distance_to_skin() # Calculate Distance to Skin


    def CutFluidMesh(self, single_parameter, min_size, max_size, hausdorff_value):
        
        find_nodal_h = KratosMultiphysics.FindNodalHNonHistoricalProcess(self.model_fluid)
        find_nodal_h.Execute()
        
        KratosMultiphysics.VariableUtils().SetNonHistoricalVariable(KratosMultiphysics.NODAL_AREA, 0.0, self.model_fluid.Nodes)
        local_gradient = KratosMultiphysics.ComputeNodalGradientProcess3D(self.model_fluid, KratosMultiphysics.DISTANCE, KratosMultiphysics.DISTANCE_GRADIENT, KratosMultiphysics.NODAL_AREA)
        local_gradient.Execute()

        # We set to zero (or unit) the metric
        ZeroVector = KratosMultiphysics.Vector(6)
        ZeroVector[0] = 0.0; ZeroVector[1] = 0.0; ZeroVector[2] = 0.0
        ZeroVector[3] = 0.0; ZeroVector[4] = 0.0; ZeroVector[5] = 0.0

        for node in self.model_fluid.Nodes:
        	node.SetValue(KratosMesh.METRIC_TENSOR_3D, ZeroVector)

        min_size = single_parameter
        max_dist = 1.25 * single_parameter
        # We define a metric using the ComputeLevelSetSolMetricProcess
        level_set_param = KratosMultiphysics.Parameters("""
        	{
        		"minimal_size"                         : """ + str(min_size) + """,
        		"enforce_current"                      : true,
        		"anisotropy_remeshing"                 : true,
        		"anisotropy_parameters": {
        			"hmin_over_hmax_anisotropic_ratio"      : 0.9,
        			"boundary_layer_max_distance"           : """ + str(max_dist) + """,
        			"interpolation"                         : "Linear" }
        	}
        	""")
        metric_process = KratosMesh.ComputeLevelSetSolMetricProcess3D(self.model_fluid, KratosMultiphysics.DISTANCE_GRADIENT, level_set_param)
        metric_process.Execute()

        # We create the remeshing process
        remesh_param = KratosMultiphysics.Parameters("""{
            "advanced_parameters": {
                "deactivate_detect_angle": false,
                "force_gradation_value": true,
                "force_hausdorff_value": true,
                "gradation_value": 0.1,
                "hausdorff_value": """ + str( hausdorff_value ) + """,
                "no_insert_mesh": false,
                "no_move_mesh": false,
                "no_surf_mesh": false,
                "no_swap_mesh": false,
                "normal_regularization_mesh": false
            },
            "buffer_size": 0,
            "debug_result_mesh": false,
            "discretization_type": "IsoSurface",
            "echo_level": 3,
            "extrapolate_contour_values": true,
            "filename": "out",
            "force_sizes": {
                "force_max": true,
                "force_min": true,
                "maximal_size": """ + str(max_size) + """,
                "minimal_size": """ + str(min_size) + """
            },
            "framework": "Eulerian",
            "initialize_entities": true,
            "internal_variables_parameters": {
                "allocation_size": 1000,
                "bucket_size": 4,
                "internal_variable_interpolation_list": [],
                "interpolation_type": "LST",
                "search_factor": 2
            },
            "interpolate_non_historical": false,
            "isosurface_parameters": {
                "isosurface_variable": "DISTANCE",
                "nonhistorical_variable": false,
                "remove_internal_regions": true
            },
            "max_number_of_searchs": 1000,
            "remesh_at_non_linear_iteration": false,
            "save_external_files": false,
            "save_mdpa_file": false,
            "search_parameters": {
                "allocation_size": 1000,
                "bucket_size": 4,
                "search_factor": 2.0
            },
            "step_data_size": 0,
            "surface_elements": false
        }""")
        MmgProcess = KratosMesh.MmgProcess3D(self.model_fluid, remesh_param)
        MmgProcess.Execute()

    def RefineAfterCut(self, level_set_size):

        if not self._refine_before_cut:
            self._detect_solid_skin() # Detect Solid Skin
            self._fluid_distance_to_skin() # Calculate Distance to Skin
        
        find_nodal_h = KratosMultiphysics.FindNodalHNonHistoricalProcess(self.model_fluid)
        find_nodal_h.Execute()
        
        KratosMultiphysics.VariableUtils().SetNonHistoricalVariable(KratosMultiphysics.NODAL_AREA, 0.0, self.model_fluid.Nodes)
        local_gradient = KratosMultiphysics.ComputeNodalGradientProcess3D(self.model_fluid, KratosMultiphysics.DISTANCE, KratosMultiphysics.DISTANCE_GRADIENT, KratosMultiphysics.NODAL_AREA)
        local_gradient.Execute()

        # We set to zero (or unit) the metric
        ZeroVector = KratosMultiphysics.Vector(6)
        ZeroVector[0] = 0.0; ZeroVector[1] = 0.0; ZeroVector[2] = 0.0
        ZeroVector[3] = 0.0; ZeroVector[4] = 0.0; ZeroVector[5] = 0.0

        for node in self.model_fluid.Nodes:
        	node.SetValue(KratosMesh.METRIC_TENSOR_3D, ZeroVector)

        min_size = level_set_size
        max_dist = 1.25 * level_set_size
        # We define a metric using the ComputeLevelSetSolMetricProcess
        level_set_param = KratosMultiphysics.Parameters("""
        	{
        		"minimal_size"                         : """ + str(min_size) + """,
        		"enforce_current"                      : true,
        		"anisotropy_remeshing"                 : true,
        		"anisotropy_parameters": {
        			"hmin_over_hmax_anisotropic_ratio"      : 0.9,
        			"boundary_layer_max_distance"           : """ + str(max_dist) + """,
        			"interpolation"                         : "Linear" }
        	}
        	""")
        metric_process = KratosMesh.ComputeLevelSetSolMetricProcess3D(self.model_fluid, KratosMultiphysics.DISTANCE_GRADIENT, level_set_param)
        metric_process.Execute()

        remesh_param = KratosMultiphysics.Parameters("""{ }""")
        MmgProcess = KratosMesh.MmgProcess3D(self.model_fluid, remesh_param)
        MmgProcess.Execute()

        if not self._refine_before_cut:
            self._fluid_distance_to_skin() # Recalculate distance to skin
        

    def RectifyModelPart(self):

        if self.model_fluid.HasSubModelPart("Boussinesq__Boussinesq_hidden_"):
            self.model_fluid.RemoveSubModelPart("Boussinesq__Boussinesq_hidden_")
            self.model_fluid.CreateSubModelPart("Boussinesq__Boussinesq_hidden_")
        
        if self.model_fluid.HasSubModelPart("TEMPERATURE_Fluid"):
            self.model_fluid.RemoveSubModelPart("TEMPERATURE_Fluid")
            self.model_fluid.CreateSubModelPart("TEMPERATURE_Fluid")
            for node in self.model_fluid.Nodes:
                self.model_fluid.GetSubModelPart("Boussinesq__Boussinesq_hidden_").AddNode(node, 0)
                self.model_fluid.GetSubModelPart("TEMPERATURE_Fluid").AddNode(node, 0)
        else:
            for node in self.model_fluid.Nodes:
                self.model_fluid.GetSubModelPart("Boussinesq__Boussinesq_hidden_").AddNode(node, 0)
        
        if self.model_solid.HasSubModelPart("TEMPERATURE_Solid"):
            self.model_solid.RemoveSubModelPart("TEMPERATURE_Solid")
            self.model_solid.CreateSubModelPart("TEMPERATURE_Solid")
            for node in self.model_solid.Nodes:
                self.model_solid.GetSubModelPart("TEMPERATURE_Solid").AddNode(node, 0)

        print(self.model_fluid)
            
    def CreateGIDOutput(self, model_name):

        KratosMultiphysics.ModelPartIO(model_name, KratosMultiphysics.IO.WRITE).WriteModelPart(self.model_fluid)
        
        self.gid_output = GiDOutputProcess(self.model_fluid, model_name, self.settings["output_configuration"])
        self.gid_output.ExecuteInitialize()
        self.gid_output.ExecuteBeforeSolutionLoop()
        self.gid_output.ExecuteInitializeSolutionStep()
        self.gid_output.PrintOutput()
        self.gid_output.ExecuteFinalizeSolutionStep()
        self.gid_output.ExecuteFinalize()

        if self.Solid_Refined:
            KratosMultiphysics.ModelPartIO("Solid_Refined", KratosMultiphysics.IO.WRITE).WriteModelPart(self.model_solid)
    
    def VTKFileOutput(self):
        
        vtk_settings = KratosMultiphysics.Parameters("""{
            "condition_data_value_variables": [],
            "condition_flags": [],
            "custom_name_postfix": "",
            "custom_name_prefix": "",
            "element_data_value_variables": [],
            "element_flags": [],
            "file_format": "ascii",
            "folder_name": "VTK_Output",
            "gauss_point_variables_extrapolated_to_nodes": [],
            "gauss_point_variables_in_elements": [],
            "model_part_name": "FluidPart",
            "nodal_data_value_variables": ["ANISOTROPIC_RATIO"],
            "nodal_flags": [],
            "nodal_solution_step_data_variables": ["DISTANCE","DISTANCE_GRADIENT"],
            "output_control_type": "step",
            "output_frequency": 1.0,
            "output_precision": 7,
            "output_sub_model_parts": true,
            "save_output_files_in_folder": false,
            "write_deformed_configuration": false,
            "write_ids": false
        }""")
        #vtk_settings["model_part_name"] = str(model_input)
        vtk_io = KratosMultiphysics.VtkOutput(self.model_fluid, vtk_settings)
        vtk_io.PrintOutput("fine_fluid")
	
if __name__ == "__main__":
    
    CHT_tool = CHTWorkflow("/mdpa_files/cavity_fluid3", "/mdpa_files/solid3D", True)
    #CHT_tool.RefineSolid(0.02)
    CHT_tool.RefineBeforeCut(1.5)
    CHT_tool.DistancetoSolidSkin()
    CHT_tool.CutFluidMesh(1.5, 2, 5, 0.1)
    #CHT_tool.RefineAfterCut(2.5)
    CHT_tool.RectifyModelPart()
    CHT_tool.VTKFileOutput()
    CHT_tool.CreateGIDOutput("fine_fluid")
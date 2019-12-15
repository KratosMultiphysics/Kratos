from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import os
import KratosMultiphysics
import KratosMultiphysics.MeshingApplication as KratosMesh
import KratosMultiphysics.FluidDynamicsApplication as KratosFluid
import KratosMultiphysics.ConvectionDiffusionApplication as KratosConvDiff
from KratosMultiphysics.gid_output_process import GiDOutputProcess

class CHTWorkflow():

    def __init__(self, model1, model2, import_settings=KratosMultiphysics.Parameters("""{}""")):
        
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
        #print(self.model_solid)

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
    
    def DetectSolidSkin(self):

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

        self.skin_model_part = self.model_solid.GetSubModelPart(self.name_auxiliar_model_part)
        print(self.skin_model_part)
        KratosMultiphysics.CalculateDistanceToSkinProcess3D(self.model_fluid, self.skin_model_part).Execute()

    def RefineSolid(self):

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

        min_size = 2.0
        max_dist = 1.25 * 2.0
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

    def RefineBeforeVisualization(self, level_set_size):

        self.DetectSolidSkin()

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

        self.DetectSolidSkin()


    def CutMeshnearSolid(self, single_parameter, min_size, max_size, hausdorff_value):
        
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
                "gradation_value": 3.0,
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

    def CleanFluidModel(self):

        self.model_fluid.RemoveSubModelPart("Boussinesq__Boussinesq_hidden_")
        self.model_fluid.RemoveSubModelPart("TEMPERATURE_fluid")

        self.model_fluid.CreateSubModelPart("Boussinesq__Boussinesq_hidden_")
        self.model_fluid.CreateSubModelPart("TEMPERATURE_fluid")

        # slip_model_part = self.model_fluid.GetSubModelPart("Slip_3D")
        # skin_isosurface_part = self.model_fluid.GetSubModelPart("SKIN_ISOSURFACE")

        # for node in skin_isosurface_part.Nodes:
        #     slip_model_part.AddNode(node, 0)
        
        # for cond in skin_isosurface_part.Conditions:
        #     slip_model_part.AddCondition(cond, 0)

        for node in self.model_fluid.Nodes:
            self.model_fluid.GetSubModelPart("Boussinesq__Boussinesq_hidden_").AddNode(node, 0)
            self.model_fluid.GetSubModelPart("TEMPERATURE_fluid").AddNode(node, 0)

        print(self.model_fluid)
            
    def CreateGIDOutput(self, model_name):

        KratosMultiphysics.ModelPartIO("Modified_Fluid3", KratosMultiphysics.IO.WRITE).WriteModelPart(self.model_fluid)
        self.gid_output = GiDOutputProcess(self.model_fluid, model_name, self.settings["output_configuration"])
        self.gid_output.ExecuteInitialize()
        self.gid_output.ExecuteBeforeSolutionLoop()
        self.gid_output.ExecuteInitializeSolutionStep()
        self.gid_output.PrintOutput()
        self.gid_output.ExecuteFinalizeSolutionStep()
        self.gid_output.ExecuteFinalize()
    
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
            "output_sub_model_parts": false,
            "save_output_files_in_folder": false,
            "write_deformed_configuration": false,
            "write_ids": false
        }""")
        #vtk_settings["model_part_name"] = str(model_input)
        vtk_io = KratosMultiphysics.VtkOutput(self.model_fluid, vtk_settings)
        vtk_io.PrintOutput()
	
if __name__ == "__main__":
    
    CHT_tool = CHTWorkflow("/mdpa_files/cavity_fluid3", "/mdpa_files/solid3D")
    CHT_tool.RefineSolid()
    #CHT_tool.EmbeddedSkinVisualization()
    CHT_tool.RefineBeforeVisualization(2.5)
    CHT_tool.CutMeshnearSolid(1.5, 2, 5, 0.1)
    CHT_tool.CleanFluidModel()
    CHT_tool.VTKFileOutput()
    CHT_tool.CreateGIDOutput("modified_fluid3")
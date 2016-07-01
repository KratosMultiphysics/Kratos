{
    "problem_data"             : {
        "problem_name"    : "test_interface5",
        "model_part_name" : "Structure",
        "domain_size"     : *GenData(DOMAIN_SIZE,INT),
        "echo_level"      : *GenData(Echo_Level)
    },
    "solver_settings"          : {
        "solver_type"                        : "solid_mechanics_static_solver",
        "echo_level"                         : 0,
        "solution_type"                      : "Static",
        "analysis_type"                      : "Linear",
        "model_import_settings"              : {
            "input_type"     : "mdpa",
            "input_filename" : "test_interface5"
        },
        "line_search"                        : *GenData(LineSearch),
        "convergence_criterion"              : "Residual_criterion",
        "displacement_relative_tolerance"    : *GenData(Convergence_Tolerance),
        "displacement_absolute_tolerance"    : *GenData(Absolute_Tolerance),
        "residual_relative_tolerance"        : *GenData(Convergence_Tolerance),
        "residual_absolute_tolerance"        : *GenData(Absolute_Tolerance),
        "max_iteration"                      : *GenData(Max_Iter,INT),
        "problem_domain_sub_model_part_list" : ["Parts_Parts_Auto2"],
        "processes_sub_model_part_list"      : ["DISPLACEMENT_Displacement_Auto1","SelfWeight2D_Self_weight_Auto1"]
    },
*set var ndomains(int) = 0
*Set cond group_DeformableBodies *groups
*loop groups *OnlyInCond           
    "meshin_domains" : {
    	"mesh_id": *cond(Group_ID),	   
	"domain_size": *GenData(DOMAIN_SIZE,INT),
	"echo_level": *GenData(Echo_Level),
	"alpha_shape": 2.4,
	"offset_factor": *GenData(Offset_Factor),
	"meshing_strategy":{
		"strategy_type": "meshing_strategy",
		"remesh": *cond(Remesh),
		"refine": *cond(Refine),
		"reconnect": *cond(Remesh),
		"transfer": *cond(Transfer),
		"constrained": *cond(Constrained),
		"mesh_smoothing": *cond(MeshSmoothing),
		"variables_smoothing": *cond(JacobiSmoothing),
		"elemental_variables_to_smooth":[ "DETERMINANT_F" ],
		"reference_element": "Element2D3N"
		"reference_condition": "CompositeCondition2D3N"
	}	
	"spatial_bounding_box":{
		"refine_in_box_only": *cond(Refine_on_box_only),
		"radius":  *cond(Radius_box),
		"center": [0.0, 0.0, 0.0],
		"velocity": [0.0, 0.0, 0.0]
	},	
	"refining_parameters":{
		"critical_size": *cond(Critical_Mesh_Size),
		"threshold_variable": *cond(Dissipation_Variable),
		"reference_threshold" : *cond(Critical_Dissipation),
		"error_variable": *cond(Error_Variable),
		"reference_error" : *cond(Critical_Error),
		"add_nodes": True,
		"insert_nodes": False,
		"remove_nodes": {
			"apply_removal": False,
			"on_distance": False,
			"on_threshold": False,
			"on_error": False
		},
		"remove_boundary": {
			"apply_removal": False,
			"on_distance": False,
			"on_threshold": False,
			"on_error": False
		},
		"refine_elements": {
			"apply_refinement": False,
			"on_distance": False,
			"on_threshold": False,
			"on_error": False
		},
		"refine_boundary": {
			"apply_refinement": False,
			"on_distance": False,
			"on_threshold": False,
			"on_error": False
		},              
		"refining_box":{
			"radius": 0.0,
			"center": [0.0, 0.0, 0.0],
			"velocity": [0.0, 0.0, 0.0]
		}
	},            
	"elemental_variables_to_transfer":[ "CAUCHY_STRESS_VECTOR", "DEFORMATION_GRADIENT" ]
    },
*end groups	
    "constraints_process_list" : [{
        "implemented_in_file"   : "apply_displacement_process",
        "implemented_in_module" : "KratosMultiphysics.SolidMechanicsApplication",
        "help"                  : "This process applies a displacement via the apply_constant_vector_value_process in kratos core",
        "process_name"          : "ApplyDisplacementProcess",
        "Parameters"            : {
            "mesh_id"         : 0,
            "model_part_name" : "DISPLACEMENT_Displacement_Auto1",
            "direction"       : [1.0,0.0,0.0]
        }
    }],
    "loads_process_list"       : [{
        "implemented_in_file"   : "apply_volume_acceleration_process",
        "implemented_in_module" : "KratosMultiphysics.SolidMechanicsApplication",
        "help"                  : "This process applies a volume acceleration node by node via the apply_constant_vector_value_process in kratos core",
        "process_name"          : "ApplyVolumeAccelerationProcess",
        "Parameters"            : {
            "mesh_id"         : 0,
            "model_part_name" : "SelfWeight2D_Self_weight_Auto1",
            "factor"          : 9.8,
            "direction"       : [0.0,-1,0.0]
        }
    }],
    "output_configuration"     : {
        "result_file_configuration" : {
            "gidpost_flags"       : {
                "GiDPostMode"           : "GiD_PostBinary",
                "WriteDeformedMeshFlag" : "WriteDeformed",
                "WriteConditionsFlag"   : "WriteConditions",
                "MultiFileFlag"         : "SingleFile"
            },
            "file_label"          : "step",
            "output_control_type" : "step",
            "output_frequency"    : *GenData(Write_Frequency),
            "body_output"         : true,
            "node_output"         : false,
            "skin_output"         : false,
            "plane_output"        : [],
            "nodal_results"       : ["DISPLACEMENT","VELOCITY","ACCELERATION","PRESSURE"],
            "gauss_point_results" : ["VON_MISES_STRESS"]
        },
        "point_data_configuration"  : []
    },
    "restart_options"          : {
        "SaveRestart"      : *GenData(Print_Restart),
        "RestartFrequency" : *GenData(Write_Frequency),
        "LoadRestart"      : *GenData(Load_Restart),
        "Restart_Step"     : *GenData(Load_Step)
    },
    "constraints_data"         : {
        "incremental_load"         : *GenData(Incremental_Load),
        "incremental_displacement" : *GenData(Incremental_Displacement)
    }
}
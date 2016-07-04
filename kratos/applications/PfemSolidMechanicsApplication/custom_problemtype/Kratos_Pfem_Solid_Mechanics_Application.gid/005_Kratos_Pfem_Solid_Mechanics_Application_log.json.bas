{
    "problem_data"             : {
        "problem_name"    : "*tcl(file tail [GiD_Info Project ModelName])",
        "model_part_name" : "Solid Domain",
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
        "line_search"                        : *tcl(string tolower *GenData(LineSearch)),
        "convergence_criterion"              : "Residual_criterion",
        "displacement_relative_tolerance"    : *GenData(Convergence_Tolerance),
        "displacement_absolute_tolerance"    : *GenData(Absolute_Tolerance),
        "residual_relative_tolerance"        : *GenData(Convergence_Tolerance),
        "residual_absolute_tolerance"        : *GenData(Absolute_Tolerance),
        "max_iteration"                      : *GenData(Max_Iter,INT),
        "problem_domain_sub_model_part_list" : ["Parts_Parts_Auto2"],
        "processes_sub_model_part_list"      : ["DISPLACEMENT_Displacement_Auto1","SelfWeight2D_Self_weight_Auto1"]
    },
    "meshing_domains" : [
*set var ndomains(int) = 0
*Set cond group_DeformableBodies *groups
*loop groups *OnlyInCond           
      	{
	"domain_type": "meshing_domain",
        "mesh_id": *cond(Group_ID),
	"domain_size": *GenData(DOMAIN_SIZE,INT),
	"echo_level": *GenData(Echo_Level),
	"alpha_shape": 2.4,
	"offset_factor": *GenData(Offset_Factor),
	"meshing_strategy":{
		"strategy_type": "meshing_strategy",
		"remesh": *tcl(string tolower *cond(Remesh)),
		"refine": *tcl(string tolower *cond(Refine)),
		"reconnect": *tcl(string tolower *cond(Remesh)),
		"transfer": *tcl(string tolower *cond(Transfer)),
		"constrained": *tcl(string tolower *cond(Constrained)),
		"mesh_smoothing": *tcl(string tolower *cond(MeshSmoothing)),
		"variables_smoothing": *tcl(string tolower *cond(JacobiSmoothing)),
		"elemental_variables_to_smooth":[ "DETERMINANT_F" ],
*if(GenData(DOMAIN_SIZE,INT)==2)
		"reference_element": "Element2D3N",
		"reference_condition": "CompositeCondition2D2N"
*elseif(GenData(DOMAIN_SIZE,INT)==3)
		"reference_element": "Element3D4N" ,
		"reference_condition": "CompositeCondition3D3N"
*endif
	},
	"spatial_bounding_box":{
		"radius": 0.0,
		"center": [0.0, 0.0, 0.0],
		"velocity": [0.0, 0.0, 0.0]
	},	
	"refining_parameters":{
		"critical_size": *cond(Critical_Mesh_Size),
		"threshold_variable": "*cond(Dissipation_Variable)",
		"reference_threshold" : *cond(Critical_Dissipation),
		"error_variable": "*cond(Error_Variable)",
		"reference_error" : *cond(Critical_Error),
		"add_nodes": true,
		"insert_nodes": false,
		"remove_nodes": {
			"apply_removal": false,
			"on_distance": true,
			"on_threshold": false,
			"on_error": true
		},
		"remove_boundary": {
			"apply_removal": false,
			"on_distance": true,
			"on_threshold": false,
			"on_error": false
		},
		"refine_elements": {
			"apply_refinement": false,
			"on_distance": true,
			"on_threshold": true,
			"on_error": false
		},
		"refine_boundary": {
			"apply_refinement": false,
			"on_distance": false,
			"on_threshold": false,
			"on_error": false
		},              
		"refining_box":{
			"refine_in_box_only": *tcl(string tolower *cond(Refine_on_box_only)),
			"radius":  *cond(Radius_box),
			"center": [*tcl(JoinByComma *cond(Center_box))],
			"velocity": [*tcl(JoinByComma *cond(Velocity_box))]
		}
	},            
	"elemental_variables_to_transfer":[ "CAUCHY_STRESS_VECTOR", "DEFORMATION_GRADIENT" ]
	}
*end groups
    ],	
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
        "SaveRestart"      : *tcl(string tolower *GenData(Print_Restart)),
        "RestartFrequency" : *GenData(Write_Frequency),
        "LoadRestart"      : *tcl(string tolower *GenData(Load_Restart)),
        "Restart_Step"     : *GenData(Load_Step)
    },
    "constraints_data"         : {
        "incremental_load"         : *tcl(string tolower *GenData(Incremental_Load)),
        "incremental_displacement" : *tcl(string tolower *GenData(Incremental_Displacement))
    }
}
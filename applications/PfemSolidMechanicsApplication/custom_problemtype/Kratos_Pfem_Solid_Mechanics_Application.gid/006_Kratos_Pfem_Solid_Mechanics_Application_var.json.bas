{
    "problem_data"             : {
        "problem_name"    : "*tcl(file tail [GiD_Info Project ModelName])",
        "model_part_name" : "Solid Domain",
        "domain_size"     : *GenData(DOMAIN_SIZE,INT),
	"time_step"       : *GenData(Time_Step),
        "start_time"      : 0.0,
        "end_time"        : 1.0,
        "echo_level"      : *GenData(Echo_Level)
    },
    "solver_settings"          : {
        "solver_type"                        : "solid_mechanics_implicit_dynamic_solver",
        "echo_level"                         : 0,
        "solution_type"                      : "Dynamic",
	"time_integration_method"            : "Implicit",
        "scheme_type"                        : "Bossak",
        "analysis_type"                      : "Non-Linear",
        "model_import_settings"              : {
*if(strcmp(GenData(Load_Restart),"True")==0)
            "input_type"       : "rest",
            "input_filename"   : "*tcl(file tail [GiD_Info Project ModelName])",
	    "input_file_label" : 0
*else
            "input_type"       : "mdpa",
            "input_filename"   : "*tcl(file tail [GiD_Info Project ModelName])",
	    "input_file_label" : *GenData(Load_Step)
*endif
        },
        "line_search"                        : *tcl(string tolower *GenData(LineSearch)),
        "convergence_criterion"              : "Residual_criterion",
        "displacement_relative_tolerance"    : *GenData(Convergence_Tolerance),
        "displacement_absolute_tolerance"    : *GenData(Absolute_Tolerance),
        "residual_relative_tolerance"        : *GenData(Convergence_Tolerance),
        "residual_absolute_tolerance"        : *GenData(Absolute_Tolerance),
        "max_iteration"                      : *GenData(Max_Iter,INT),
        "linear_solver_settings"             : {
             "solver_type" : "Super_LU",
             "scaling"     : true,
             "verbosity"   : 0
         },
        "problem_domain_sub_model_part_list" : [
*Set cond group_DeformableBodies *groups
*if(CondNumEntities > 0)
*set var GroupNum = 0
*loop groups *OnlyInCond
*set var GroupNum=operation(GroupNum+1)
*end groups
*set var Counter = 0
*loop groups *OnlyInCond
*set var Counter=operation(Counter+1)
*if( Counter == GroupNum )
         "*cond(Group_ID)"
*else
	 "*cond(Group_ID)",
*end
*end groups
*endif
      	],
        "processes_sub_model_part_list" : [
*Set cond group_RigidWalls *groups
*if(CondNumEntities > 0)
*set var GroupNum = 0
*loop groups *OnlyInCond
*set var GroupNum=operation(GroupNum+1)
*end groups

*set var Counter = 0
*loop groups *OnlyInCond
*set var Counter=operation(Counter+1)
*if( Counter == GroupNum )
         "*cond(Group_ID)"
*else
	 "*cond(Group_ID)",
*end
*end groups
*endif
      	]
    },
    "problem_process_list" : [{
        "implemented_in_file"   : "remesh_domains_process",
        "implemented_in_module" : "KratosMultiphysics.PfemBaseApplication",
        "help"                  : "This process applies meshing to the problem domains",
        "process_name"          : "RemeshDomainsProcess",
        "Parameters"            : {
	    "model_part_name"       : "Solid Domain",
            "meshing_control_type"  : "step",
            "meshing_frequency"     : 0.0,
            "meshing_before_output" : true,
	    "meshing_domains" : [
*set var ndomains(int) = 0
*Set cond group_DeformableBodies *groups
*loop groups *OnlyInCond           
            {
		"python_file_name": "meshing_domain",
		"mesh_id": *cond(Group_ID),
		"sub_model_part_name": "*cond(Group_ID)",
		"alpha_shape": 2.4,
		"offset_factor": *GenData(Offset_Factor),
		"meshing_strategy":{
		    "python_file_name": "meshing_strategy",
		    "meshing_frequency": *cond(Meshing_Frequency),
 		    "remesh": *tcl(string tolower *cond(Remesh)),
		    "refine": *tcl(string tolower *cond(Refine)),
		    "reconnect": *tcl(string tolower *cond(Remesh)),
		    "transfer": *tcl(string tolower *cond(Transfer)),
		    "constrained": *tcl(string tolower *cond(Constrained)),
		    "mesh_smoothing": *tcl(string tolower *cond(MeshSmoothing)),
		    "variables_smoothing": *tcl(string tolower *cond(JacobiSmoothing)),
		    "elemental_variables_to_smooth":[ "DETERMINANT_F" ],
*if(GenData(DOMAIN_SIZE,INT)==2)
		    "reference_element_type": "Element2D3N",
		    "reference_condition_type": "CompositeCondition2D2N"
*elseif(GenData(DOMAIN_SIZE,INT)==3)
		    "reference_element_type": "Element3D4N" ,
		    "reference_condition_type": "CompositeCondition3D3N"
*endif
	    },
	    "spatial_bounding_box":{
	    	"upper_point": [0.0, 0.0, 0.0],
		"lower_point": [0.0, 0.0, 0.0],
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
            ]
        }
    },{	
        "implemented_in_file"   : "parametric_walls_process",
        "implemented_in_module" : "KratosMultiphysics.ContactMechanicsApplication",
        "help"                  : "This process applies parametric walls and search contact",
        "process_name"          : "ParametricWallsProcess",
        "Parameters"            : {
	    "model_part_name"      : "Solid Domain",
            "search_control_type"  : "step",
            "search_frequency"     : 0.0,
	    "parametric_walls" : [
*Set cond group_RigidWalls *groups
*loop groups *OnlyInCond    
		{
		    "python_file_name": "parametric_wall",
		    "mesh_id": *cond(Group_ID),
		    "sub_model_part_name" : "*cond(Group_ID)",
		    "rigid_body_settings":{
			"rigid_body_element_type": "TranslatoryRigidElement3D1N",
			"fixed_body": true,
			"compute_parameters": true,
			"rigid_body_parameters":{
			    "center_of_gravity": [0.0 ,0.0, 0.0],
			    "mass":0.0,
			    "main_inertias": [0.0, 0.0, 0.0],
			    "main_axes": [ [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0] ]
			}
		    },
		    "bounding_box_settings":{
			"implemented_in_module": "KratosMultiphysics.ContactMechanicsApplication",
*if(strcmp(cond(Wall_Type),"PLANE")==0)
			"bounding_box_type": "PlaneBoundingBox",
			"bounding_box_parameters":{
			    "parameters_list":[
				"point": [*cond(Wall_Plane,1), *cond(Wall_Plane,2), *cond(Wall_Plane,3)],
			    	"normal": [*cond(Wall_Plane,4), *cond(Wall_Plane,5), *cond(Wall_Plane,6)],
				"convexity": *cond(Wall_Circle,5)
			    ],
			    "velocity" : [*tcl(JoinByComma *cond(Linear_Velocity))]
			}
*elseif(strcmp(cond(Wall_Type),"CIRCLE")==0)
			"bounding_box_type": "CircleBoundingBox",
			"bounding_box_parameters":{
			    "parameters_list":[
				"center": [*cond(Wall_Circle,1), *cond(Wall_Circle,2), *cond(Wall_Circle,3)],
				"radius": *cond(Wall_Circle,4),
				"convexity": *cond(Wall_Circle,5)
			    ],
			    "velocity" : [0.0, 0.0, 0.0]
			}
*elseif(strcmp(cond(Wall_Type),"NOSE-WALL")==0)
			"bounding_box_type": "CompoundNosesBoundingBox",
			"bounding_box_parameters":{
			    "parameters_list":[
*for(i=1;i<=cond(Wall_Noses,INT);i=i+7)
				{  
				   "radius": *cond(Wall_Noses,*Operation(i+3)),		       
				   "center": [*cond(Wall_Noses,*i), *cond(Wall_Noses,*Operation(i+1)), *cond(Wall_Noses,*Operation(i+2)))],
			    	   "rake_angle": *cond(Wall_Noses,*Operation(i+4)),
			    	   "clearance_angle": *cond(Wall_Noses,*Operation(i+5)),
			    	   "convexity": *cond(Wall_Noses,*Operation(i+6))
*if( i<cond(Wall_Noses,INT) )
			     	},
*else
				}
*endif
*end
			    ],
			    "velocity" : [0.0, 0.0, 0.0]
			}

*endif
		    },
		    "contact_search_settings":{
			"python_file_name": "parametric_wall_contact_search",
			"search_frequency": 0,            
			"contact_parameters":{
			    "contact_condition_type": "*cond(Contact_Condition)",
			    "friction_law_type": "HardeningCoulombFrictionLaw",
			    "implemented_in_module": "KratosMultiphysics.ContactMechanicsApplication",
			    "variables_of_properties":{
				"FRICTION_ACTIVE": false,
				"MU_STATIC": 0.3,
				"MU_DYNAMIC": 0.2,
				"PENALTY_PARAMETER": *cond(Penalty_Parameter),
				"TANGENTIAL_PENALTY_RATIO": 0.1,
				"TAU_STAB": 1
			    }
			}
		    }
		}
*end groups
	    ]
	}
    },{
        "implemented_in_file"   : "contact_domain_process",
        "implemented_in_module" : "KratosMultiphysics.ContactMechanicsApplication",
        "help"                  : "This process applies contact domain search by remeshing outer boundaries",
        "process_name"          : "ContactDomainProcess",
        "Parameters"            : {
            "mesh_id"               : 0,
	    "model_part_name"       : "Solid Domain",
            "meshing_control_type"  : "step",
            "meshing_frequency"     : 0.0,
            "meshing_before_output" : true,
	    "meshing_domains" : [
		{
		    "python_file_name": "contact_domain",
		    "sub_model_part_name": "contact_domain_domain",
		    "alpha_shape": 1.4,
		    "offset_factor": *GenData(Offset_Factor),
		    "meshing_strategy":{
			"python_file_name": "contact_meshing_strategy",
			"meshing_frequency": *GenData(Contact_Search_Frequency),
			"remesh": true,
			"constrained": *tcl(string tolower *GenData(Constrained_Contact)),
			"contact_parameters":{
			    "contact_condition_type": "*GenData(ContactCondition)",
			    "friction_law_type": "HardeningCoulombFrictionLaw",
			    "implemented_in_module": "KratosMultiphysics.ContactMechanicsApplication",
			    "variables_of_properties":{
				"FRICTION_ACTIVE": *tcl(string tolower *GenData(Friction_Active)),
				"MU_STATIC": 0.3,
				"MU_DYNAMIC": 0.2,
				"PENALTY_PARAMETER": *GenData(Penalty_Parameter),
				"TANGENTIAL_PENALTY_RATIO": 0.1,
				"TAU_STAB": *GenData(Stability_Parameter)
			    }
			}
		    },
		    "elemental_variables_to_transfer":[ "CAUCHY_STRESS_VECTOR", "DEFORMATION_GRADIENT" ]
		}
            ]
        }
    }],
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
    "output_process_list" : [{
        "implemented_in_file"   : "restart_process",
        "implemented_in_module" : "KratosMultiphysics.SolidMechanicsApplication",
        "help"                  : "This process writes restart files",
        "process_name"          : "RestartProcess",
        "Parameters"            : {
            "save_restart"        : *tcl(string tolower *GenData(Print_Restart)),
            "restart_file_name"   : "*tcl(file tail [GiD_Info Project ModelName])",
            "restart_file_label"  : "step",
            "output_control_type" : "step",
            "output_frequency"    : *GenData(Write_Frequency),
            "json_output"         : false
        }
    }],
    "output_configuration"     : {
        "result_file_configuration" : {
            "gidpost_flags"       : {
*if(strcmp(GenData(File_Format),"binary")==0)
                "GiDPostMode"           : "GiD_PostBinary",
*elseif(strcmp(GenData(File_Format),"ascii")==0)
                "GiDPostMode"           : "GiD_PostAscii",
*endif
*if(strcmp(GenData(Write_Mesh),"Deformed")==0)
                "WriteDeformedMeshFlag" : "WriteDeformed",
*else
                "WriteDeformedMeshFlag" : "WriteUndeformed",
*endif
*if(strcmp(GenData(Write_Conditions),"True")==0)
		"WriteConditionsFlag"   : "WriteConditions",
*else
                "WriteConditionsFlag" : "WriteElementsOnly",
*endif
                "MultiFileFlag"         : "MultipleFiles"
            },
            "file_label"          : "step",
            "output_control_type" : "time",
            "output_frequency"    : *GenData(Write_Frequency),
            "body_output"         : true,
            "node_output"         : false,
            "skin_output"         : false,
            "plane_output"        : [],
            "nodal_results"       : ["DISPLACEMENT","VELOCITY","ACCELERATION","PRESSURE"],
            "gauss_point_results" : ["CAUCHY_STRESS_TENSOR","GREEN_LAGRANGE_STRAIN_TENSOR","VON_MISES_STRESS"]
        },
        "point_data_configuration"  : []
    }

}
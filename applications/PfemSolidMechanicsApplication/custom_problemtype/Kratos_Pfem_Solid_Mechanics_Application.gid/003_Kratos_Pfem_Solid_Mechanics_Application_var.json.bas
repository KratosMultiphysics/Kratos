{
    "problem_data"             : {
        "problem_name"    : "*tcl(file tail [GiD_Info Project ModelName])",
        "threads"         : *GenData(Number_of_threads,INT),
        "echo_level"      : *GenData(Echo_Level)
    },
    "time_settings"             : {
         "time_step"  : *GenData(Time_Step),
         "start_time" : *GenData(Start_Time),
         "end_time"   : *GenData(End_Time)
    },
    "model_settings"           : {
        "model_name": "Main_Domain",
        "dimension"       : *GenData(DIMENSION,INT),
	"bodies_list":[
*set cond group_DeformableBodies *groups
*add cond group_RigidBodies *groups
*if(CondNumEntities > 0)
*set var GroupNumber = 0
*loop groups *OnlyInCond
*set var GroupNumber=operation(GroupNumber+1)
*end groups
*set var Counter = 0
*loop groups *OnlyInCond
*set var Counter=operation(Counter+1)
     	 {
     	 "body_type": "*cond(Structural_Type)",
     	 "body_name": "*cond(Structural_Type)_*GroupName",
     	 "parts_list": ["*GroupName"]
*if( Counter == GroupNumber )
     	 }
*else
     	 },
*endif
*end groups
*endif
	],
        "domain_parts_list" : [
*set cond group_DeformableBodies *groups
*add cond group_RigidBodies *groups
*if(CondNumEntities > 0)
*set var GroupNumber = 0
*loop groups *OnlyInCond
*set var GroupNumber=operation(GroupNumber+1)
*end groups
*set var Counter = 0
*loop groups *OnlyInCond
*set var Counter=operation(Counter+1)
*if( Counter == GroupNumber )
     	 "*GroupName"
*else
     	 "*GroupName",
*endif
*end groups
*endif
	],
        "processes_parts_list" : [
*set cond group_LINEAR_MOVEMENT *groups
*add cond group_ANGULAR_MOVEMENT *groups
*add cond group_POINT_LOAD *groups
*add cond group_LINE_LOAD *groups
*add cond group_SURFACE_LOAD *groups
*add cond group_LINE_PRESSURE *groups
*add cond group_SURFACE_PRESSURE *groups
*add cond group_POINT_MOMENT *groups
*add cond group_VOLUME_ACCELERATION *groups
*add cond group_WATER_PRESSURE *groups
*add cond group_WATER_MOVEMENT *groups
*if(CondNumEntities > 0)
*set var GroupNumber = 0
*loop groups *OnlyInCond
*set var GroupNumber=operation(GroupNumber+1)
*end groups
*set var Counter = 0
*loop groups *OnlyInCond
*set var Counter=operation(Counter+1)
*if( Counter == GroupNumber )
     	 "*GroupName"
*else
     	 "*GroupName",
*endif
*end groups
*endif
	],
        "input_file_settings"     : {
*if(strcmp(GenData(Load_Restart),"True")==0)
            "type"       : "rest",
            "name"   : "*tcl(file tail [GiD_Info Project ModelName])",
	    "label" : *GenData(Load_Step)
*else
            "type"       : "mdpa",
            "name"   : "*tcl(file tail [GiD_Info Project ModelName])",
	    "label" : 0
*endif
        }
    },
    "solver_settings"          : {
*if(strcmp(GenData(Solver_Type),"DynamicSolver")==0)
*if(strcmp(GenData(Time_Integration_Method),"Explicit")==0)
        "solver_type" : "solid_mechanics_explicit_dynamic_solver",
*elseif(strcmp(GenData(Time_Integration_Method),"Implicit")==0)
        "solver_type" : "solid_mechanics_implicit_dynamic_solver",
*endif
*else
        "solver_type" : "solid_mechanics_static_solver",
*endif
        "Parameters"  : {       
              "time_integration_settings" : {
*if(strcmp(GenData(Solver_Type),"DynamicSolver")==0)
                   "solution_type"         : "Dynamic",
*if(strcmp(GenData(Time_Integration_Method),"Explicit")==0)
                   "time_integration"      : "Explicit",
                   "integration_method"    : "CentralDifferences"
*elseif(strcmp(GenData(Time_Integration_Method),"Implicit")==0)
*if(strcmp(GenData(DOFS),"U-W")==0)
                   "time_integration"      : "Implicit",
                   "integration_method"    : "Bossak",
		   "lumped_matrix": false,
		   "consistent_mass_matrix": true
*elseif(strcmp(GenData(DOFS),"U-W-wP")==0)
                   "time_integration"      : "Implicit",
                   "integration_method"    : "Bossak",
		   "lumped_matrix": false,
		   "consistent_mass_matrix": true
*else
                   "time_integration"      : "Implicit",
                   "integration_method"    : "Bossak",
		   "lumped_matrix": false,
		   "consistent_mass_matrix": true
*endif
*endif
*else
*if(strcmp(GenData(Solver_Type),"StaticSolver")==0)
                   "solution_type"         : "Static",
                   "integration_method"    : "Static"
*elseif(strcmp(GenData(Solver_Type),"QuasiStaticSolver")==0)
                   "solution_type"         : "Quasi-static",
                   "integration_method"    : "Static"
*endif
*endif
              },
              "solving_strategy_settings" : {
                   "line_search"                 : *tcl(string tolower *GenData(LineSearch)),
                   "implex"                      : *tcl(string tolower *GenData(Implex)),
                   "compute_reactions"           : *tcl(string tolower *GenData(Write_Reactions)),
	           "compute_contact_forces"      : *tcl(string tolower *GenData(Write_Contact_Forces)),
                   "max_iteration"               : *GenData(Max_Iter,INT)
              },
              "convergence_criterion_settings":{
                   "convergence_criterion"       : "*GenData(Convergence_Criteria)",
                   "reform_dofs_at_each_step"    : true,
                   "variable_relative_tolerance" : *GenData(Convergence_Tolerance),
                   "variable_absolute_tolerance" : *GenData(Absolute_Tolerance),
                   "residual_relative_tolerance" : *GenData(Convergence_Tolerance),
                   "residual_absolute_tolerance" : *GenData(Absolute_Tolerance)
              },
              "linear_solver_settings"   : {
                   "solver_type"    : "*GenData(Linear_Solver)",
                   "tolerance"      : 1e-7,
                   "max_iteration"  : *GenData(Linear_Solver_Max_Iteration,INT),
                   "scaling"        : false
              },
              "dofs"                            : [
*if(strcmp(GenData(DOFS),"DISPLACEMENTS")==0)
                                                "DISPLACEMENT"
*endif
*if(strcmp(GenData(DOFS),"ROTATIONS")==0)
                                                "DISPLACEMENT",
                                                "ROTATION"
*endif
*if(strcmp(GenData(DOFS),"U-P")==0)
                                                "DISPLACEMENT",
                                                "PRESSURE"
*endif
*if(strcmp(GenData(DOFS),"U-wP")==0 )
                                                "DISPLACEMENT",
                                                "WATER_PRESSURE"
*endif
*if( strcmp(GenData(DOFS),"U-J-wP")==0 )
                                                "DISPLACEMENT",
                                                "WATER_PRESSURE",
						"JACOBIAN"
*endif
*if( strcmp(GenData(DOFS),"U-J")==0 )
                                                "DISPLACEMENT",
						"JACOBIAN"
*endif
*if(strcmp(GenData(DOFS),"U-W")==0)
                                                "DISPLACEMENT",
						"WATER_DISPLACEMENT"
*endif
*if(strcmp(GenData(DOFS),"U-W-wP")==0)
                                                "DISPLACEMENT",
						"WATER_DISPLACEMENT",
                                                "WATER_PRESSURE"
*endif
*if(strcmp(GenData(DOFS),"U-J-W-wP")==0)
                                                "DISPLACEMENT",
						"JACOBIAN",
						"WATER_DISPLACEMENT",
                                                "WATER_PRESSURE"
*endif
					       ]
        }
    },
    "problem_process_list" : [{
        "help"            : "This process applies meshing to the problem domains",
        "kratos_module"   : "KratosMultiphysics.DelaunayMeshingApplication",
        "python_module"   : "remesh_domains_process",
        "process_name"    : "RemeshDomainsProcess",
        "Parameters"      : {
	    "model_part_name"       : "Main_Domain",
            "meshing_control_type"  : "step",
            "meshing_frequency"     : 0.0,
            "meshing_before_output" : true,
	    "meshing_domains" : [
*set var numberofdomains= 0
*set cond group_DeformableBodies *groups
*loop groups *OnlyInCond
*set var numberofdomains(int)=Operation(numberofdomains(int)+1)
*end groups
*set var Counter = 0
*set cond group_DeformableBodies *groups
*loop groups *OnlyInCond           
*set var Counter=operation(Counter+1)
            {
		"python_module": "meshing_domain",
		"model_part_name": "*cond(Structural_Type)_*GroupName",
		"alpha_shape": 2.4,
		"offset_factor": *GenData(Offset_Factor),
		"meshing_strategy":{
		    "python_module": "meshing_strategy",
		    "meshing_frequency": *cond(Meshing_Frequency),
		    "remesh": *tcl(string tolower *cond(Remesh)),
		    "refine": *tcl(string tolower *cond(Refine)),
		    "reconnect": *tcl(string tolower *cond(Remesh)),
		    "transfer": *tcl(string tolower *cond(Transfer)),
		    "constrained": *tcl(string tolower *cond(Constrained)),
		    "mesh_smoothing": *tcl(string tolower *cond(MeshSmoothing)),
		    "variables_smoothing": *tcl(string tolower *cond(JacobiSmoothing)),
		    "elemental_variables_to_smooth":[ "DETERMINANT_F" ],
*if(GenData(DIMENSION,INT)==2)
		    "reference_element_type": "Element2D3N",
		    "reference_condition_type": "CompositeCondition2D2N"
*elseif(GenData(DIMENSION,INT)==3)
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
			"upper_point": [*tcl(JoinByComma *cond(Upper_Point_rbox))],
			"lower_point": [*tcl(JoinByComma *cond(Lower_Point_rbox))],
			"velocity": [*tcl(JoinByComma *cond(Velocity_rbox))]
		}
	    },            
	    "elemental_variables_to_transfer":[ "CAUCHY_STRESS_VECTOR", "DEFORMATION_GRADIENT" ]
*if( Counter == numberofdomains )
        }
*else
	},
*endif
*end groups
            ]
        }
*if(strcmp(GenData(FindContacts),"True")==0)
    },{
        "help"             : "This process applies contact domain search by remeshing outer boundaries",    
        "kratos_module"    : "KratosMultiphysics.ContactMechanicsApplication",   
        "python_module"    : "contact_domain_process",
        "process_name"     : "ContactDomainProcess",
        "Parameters"       : {
	    "model_part_name"       : "Main_Domain",
            "meshing_control_type"  : "step",
            "meshing_frequency"     : 0.0,
            "meshing_before_output" : true,
	    "meshing_domains" : [
		{
		    "python_module": "contact_domain",
		    "model_part_name": "contact_domain",
		    "alpha_shape": 1.4,
		    "offset_factor": *GenData(Offset_Factor),
		    "meshing_strategy":{
			"python_module": "contact_meshing_strategy",
			"meshing_frequency": *GenData(Contact_Search_Frequency),
			"remesh": true,
			"constrained": *tcl(string tolower *GenData(Constrained_Contact)),
			"contact_parameters":{
			    "contact_condition_type": "*GenData(ContactCondition)",
			    "kratos_module": "KratosMultiphysics.ContactMechanicsApplication",			    
			    "friction_law_type": "HardeningCoulombFrictionLaw",
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
		    "elemental_variables_to_transfer":[ "CAUCHY_STRESS_VECTOR", "DEFORMATION_GRADIENT" ],
                    "contact_bodies_list" : [
*set cond group_DeformableBodies *groups
*if(CondNumEntities > 0)
*set var GroupNumber = 0
*loop groups *OnlyInCond
*set var GroupNumber=operation(GroupNumber+1)
*end groups
*set var Counter = 0
*loop groups *OnlyInCond
*set var Counter=operation(Counter+1)
*if( Counter == GroupNumber )
         "*cond(Structural_Type)_*GroupName"
*else
	 "*cond(Structural_Type)_*GroupName",
*endif
*end groups
*endif
      	            ]
		}
            ]
        }
*endif
*set var numberofwalls= 0
*set cond group_RigidBodies *groups
*loop groups *OnlyInCond
*if(strcmp(cond(Parametric_Wall),"True")==0)
*set var numberofwalls(int)=Operation(numberofwalls(int)+1)
*endif
*end groups
*if( numberofwalls > 0 )
    },{
        "help"            : "This process applies parametric walls and search contact",
	"kratos_module"   : "KratosMultiphysics.ContactMechanicsApplication", 	       
        "python_module"   : "parametric_walls_process",
        "process_name"    : "ParametricWallsProcess",
        "Parameters"      : {
	    "model_part_name"      : "Main_Domain",
            "search_control_type"  : "step",
            "search_frequency"     : 0.0,
	    "parametric_walls" : [
*set var Counter = 0
*set cond group_RigidBodies *groups
*loop groups *OnlyInCond
*if(strcmp(cond(Parametric_Wall),"True")==0)
*set var Counter=operation(Counter+1)
		{
		    "python_module": "parametric_wall",
		    "model_part_name" : "*cond(Structural_Type)_*GroupName",
		    "rigid_body_settings":{
			"rigid_body_element_type": "*cond(Rigid_Element)",
			"fixed_body": true,
			"compute_body_parameters": *tcl(string tolower *cond(Compute_Weight_Centroid_and_Inertia)),
			"rigid_body_model_part_name": "*cond(Structural_Type)_*GroupName",
			"rigid_body_parameters":{
			    "center_of_gravity": [*tcl(JoinByComma *cond(Centroid))],
			    "mass": *cond(Mass),
			    "main_inertias": [*cond(LocalInertiaTensor,1), *cond(LocalInertiaTensor,2), *cond(LocalInertiaTensor,3)],
			    "main_axes": [ [*tcl(JoinByComma *cond(LocalAxisX))], [*tcl(JoinByComma *cond(LocalAxisY))], [*tcl(JoinByComma *cond(LocalAxisZ))] ]
			}
		    },
		    "bounding_box_settings":{
			"kratos_module": "KratosMultiphysics.ContactMechanicsApplication",
*if(strcmp(cond(Wall_Type),"PLANE")==0)
			"bounding_box_type": "PlaneBoundingBox",
			"bounding_box_parameters":{
			    "parameters_list":[
			      {
				"point": [*cond(Wall_Plane,1), *cond(Wall_Plane,2), *cond(Wall_Plane,3)],
				"normal": [*cond(Wall_Plane,4), *cond(Wall_Plane,5), *cond(Wall_Plane,6)],
				"convexity": *cond(Wall_Circle,5)
			      }
			    ],
			    "velocity" : [*tcl(JoinByComma *cond(Linear_Velocity))],
			    "plane_size": *cond(Plane_Size)
			}
*elseif(strcmp(cond(Wall_Type),"CIRCLE")==0)
			"bounding_box_type": "CircleBoundingBox",
			"bounding_box_parameters":{
			    "parameters_list":[
			      {
				"center": [*cond(Wall_Circle,1), *cond(Wall_Circle,2), *cond(Wall_Circle,3)],
				"radius": *cond(Wall_Circle,4),
				"convexity": *cond(Wall_Circle,5)
			      }
			    ],
			    "velocity" : [0.0, 0.0, 0.0]
			}
*elseif(strcmp(cond(Wall_Type),"CYLINDER")==0)
			"bounding_box_type": "CylinderBoundingBox",
			"bounding_box_parameters":{
			    "parameters_list":[
			      {
				"first_center": [*cond(Wall_Cylinder,1), *cond(Wall_Cylinder,2), *cond(Wall_Cylinder,3)],
				"second_center": [*cond(Wall_Cylinder,5), *cond(Wall_Cylinder,6), *cond(Wall_Cylinder,7)],	
				"radius": *cond(Wall_Cylinder,4),
				"convexity": *cond(Wall_Cylinder,8)
			      }
			    ],
			    "velocity" : [0.0, 0.0, 0.0]
			}
*elseif(strcmp(cond(Wall_Type),"SPHERE")==0)
			"bounding_box_type": "SphereBoundingBox",
			"bounding_box_parameters":{
			    "parameters_list":[
			      {
				"center": [*cond(Wall_Sphere,1), *cond(Wall_Sphere,2), *cond(Wall_Sphere,3)],
				"radius": *cond(Wall_Sphere,4),
				"convexity": *cond(Wall_Sphere,5)
			      }
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
				   "center": [*cond(Wall_Noses,*i), *cond(Wall_Noses,*Operation(i+1)), *cond(Wall_Noses,*Operation(i+2))],
				   "rake_angle": *cond(Wall_Noses,*Operation(i+4)),
				   "clearance_angle": *cond(Wall_Noses,*Operation(i+5)),
				   "convexity": *cond(Wall_Noses,*Operation(i+6))
*if( i+7>=cond(Wall_Noses,INT) )
				}
*else
				},
*endif
*end
			    ],
			    "velocity" : [0.0, 0.0, 0.0],
			    "plane_size": *cond(Plane_Size)
			}

*endif
		    },		    
		    "contact_search_settings":{
			"kratos_module": "KratosMultiphysics.ContactMechanicsApplication",
*if(strcmp(cond(Hydraulic_Condition),"True")==0)
			"contact_search_type": "HMParametricWallContactSearch",
*else
			"contact_search_type": "ParametricWallContactSearch",
*endif
			"contact_parameters":{
			    "contact_condition_type": "*cond(Contact_Condition)",
*if(strcmp(cond(Hydraulic_Condition),"True")==0)
			    "hydraulic_condition_type": "*cond(Hydraulic_Contact_Condition)",
*endif
			    "kratos_module": "KratosMultiphysics.ContactMechanicsApplication",			    
			    "friction_law_type": "HardeningCoulombFrictionLaw",
			    "variables_of_properties":{
				"FRICTION_ACTIVE": *tcl(string tolower *cond(Friction_active)),
				"MU_STATIC": *cond(Static_friction_coeffitient),
				"MU_DYNAMIC": *cond(Dynamic_friction_coeffitient),
				"PENALTY_PARAMETER": *cond(Penalty_Parameter),
				"TANGENTIAL_PENALTY_RATIO": *cond(Tangent_penalty_ratio),
				"TAU_STAB": 1
			    }
			}
		    }
*if( Counter == numberofwalls )
		}
*else
		},
*endif
*endif
*end groups
	    ]
	}
*endif	
*set var numberofrigidbodies= 0
*set cond group_RigidBodies *groups
*loop groups *OnlyInCond
*if(strcmp(cond(Parametric_Wall),"False")==0)
*set var numberofrigidbodies(int)=Operation(numberofrigidbodies(int)+1)
*endif
*end groups
*if( numberofrigidbodies > 0 )
    },{
        "help"         : "This process creates the rigid bodies of the model",
        "kratos_module": "KratosMultiphysics.ContactMechanicsApplication",
        "python_module": "rigid_bodies_process",
        "process_name" : "RigidBodyProcess",
	"Parameters"   : {
	    "model_part_name" : "Main_Domain",
	    "rigid_bodies"    : [
*set var Counter = 0
*set cond group_RigidBodies *groups
*loop groups *OnlyInCond
*if(strcmp(cond(Parametric_Wall),"False")==0)
*set var Counter=operation(Counter+1)
		{
		    "python_module"  : "rigid_body",
		    "model_part_name": "*cond(Structural_Type)_*GroupName",
		    "rigid_body_settings":{
			"rigid_body_element_type": "RigidBodyElement2D1N",
			"fixed_body": false,
			"compute_body_parameters": *tcl(string tolower *cond(Compute_Weight_Centroid_and_Inertia)),
			"rigid_body_model_part_name": "*cond(Structural_Type)_*GroupName",
			"rigid_body_parameters":{
			    "center_of_gravity": [*tcl(JoinByComma *cond(Centroid))],
			    "mass": *cond(Mass),
			    "main_inertias": [*cond(LocalInertiaTensor,1), *cond(LocalInertiaTensor,2), *cond(LocalInertiaTensor,3)],
			    "main_axes": [ [*tcl(JoinByComma *cond(LocalAxisX))], [*tcl(JoinByComma *cond(LocalAxisY))], [*tcl(JoinByComma *cond(LocalAxisZ))] ]
			}
		    }
*if( Counter == numberofrigidbodies )
		}
*else
		},
*endif
*endif
*end groups
	    ]
    	}
*endif
    }],
    "constraints_process_list" : [
*set var numberconstraints= 0
*set cond group_WATER_MOVEMENT *groups 
*add cond group_LINEAR_MOVEMENT *groups
*add cond group_ANGULAR_MOVEMENT *groups
*add cond group_WATER_PRESSURE *groups
*loop groups *OnlyInCond
*set var numberconstraints(int)=Operation(numberconstraints(int)+1)
*end groups
*set var Counter = 0
*set cond group_LINEAR_MOVEMENT *groups
*add cond group_ANGULAR_MOVEMENT *groups
*add cond group_WATER_MOVEMENT *groups
*if(strcmp(GenData(Set_initial_state),"True")==0)
*set var numberconstraints(int)=Operation(numberconstraints(int)+1)
*endif
*loop groups *OnlyInCond
*set var Counter=operation(Counter+1)
     	{
        "help"            : "This process imposes a vector constraint",	
        "kratos_module"   : "KratosMultiphysics.SolidMechanicsApplication",
        "python_module"   : "assign_vector_components_to_nodes_process",
        "process_name"    : "AssignVectorComponentsToNodesProcess",
        "Parameters"      : {
            "model_part_name" : "*GroupName",
            "variable_name"   : "*cond(Variable)",
*if(strcmp(cond(Time_Evolution),"INITIAL")==0)
	    "interval"        : [*GenData(Start_Time), *GenData(Start_Time)],
*elseif(strcmp(cond(Time_Evolution),"FULL")==0)
	    "interval"        : [*GenData(Start_Time), "End"],
*elseif(strcmp(cond(Time_Evolution),"INTERVAL")==0)
	    "interval"        : [*cond(Time_Interval,1), *cond(Time_Interval,2)],
*endif
	    "value"           : [
*if(strcmp(cond(Set_X),"True")==0 )
*if(strcmp(cond(by_function_X),"True")==0 )
	    			"*cond(X_Function)",
*else
	    			*cond(X_Value),
*endif				
*else
				null,
*endif				
*if(strcmp(cond(Set_Y),"True")==0 )
*if(strcmp(cond(by_function_Y),"True")==0 )
	    			"*cond(Y_Function)",
*else
	    			*cond(Y_Value),
*endif				
*else
				null,
*endif
*if(strcmp(cond(Set_Z),"True")==0 )
*if(strcmp(cond(by_function_Z),"True")==0 )
	    			"*cond(Z_Function)"
*else
	    			*cond(Z_Value)
*endif				
*else
				null
*endif
				]
	    }
*if( Counter == numberconstraints )
        }
*else
	},
*endif
*end groups
*set cond group_WATER_PRESSURE *groups
*loop groups *OnlyInCond
*set var Counter=operation(Counter+1)
        {
        "help"            : "This process tries to fix water pressure",
        "kratos_module"   : "KratosMultiphysics.SolidMechanicsApplication",
        "python_module"   : "assign_scalar_to_nodes_process",
        "process_name"    : "AssignScalarToNodesProcess",
	"Parameters"      : {
	   "variable_name":    "WATER_PRESSURE",
           "model_part_name":  "*GroupName",
*if(strcmp(cond(Time_Evolution),"INITIAL")==0)
	    "interval"        : [*GenData(Start_Time), *GenData(Start_Time)],
*elseif(strcmp(cond(Time_Evolution),"FULL")==0)
	    "interval"        : [*GenData(Start_Time), "End"],
*elseif(strcmp(cond(Time_Evolution),"INTERVAL")==0)
	    "interval"        : [*cond(Time_Interval,1), *cond(Time_Interval,2)],
*endif
*if(strcmp(cond(by_function),"True")==0 )
	    "value"           :	"*cond(Function)"
*else
	    "value"           :	*cond(Value)
*endif			
	    }
*if( Counter == numberconstraints )
        }
*else
	},
*endif
*end groups
*if(strcmp(GenData(Set_initial_state),"True")==0)
        {
        "help"            : "This process tries to set the initial state of the soil",
        "kratos_module"   : "KratosMultiphysics.PfemSolidMechanicsApplication",
        "python_module"   : "assign_initial_HM_state_process",
        "process_name"    : "SetInitialStateProcess",
        "Parameters"      : {
             "model_part_name": "Main_Domain",
*if(strcmp(GenData(Constant_weight),"True")==0)
             "gravity_active":              false,
             "constant_horizontal_stress":  *GenData(SX), 
             "constant_vertical_stress":    *GenData(SY),
             "constant_water_pressure":     *GenData(WP)
*else 
             "gravity_active":   true
*endif
             }
         }
*endif
    ],
    "loads_process_list"       : [
*set var numberloads= 0
*set cond group_POINT_LOAD *groups
*add cond group_LINE_LOAD *groups
*add cond group_SURFACE_LOAD *groups
*add cond group_LINE_PRESSURE *groups
*add cond group_SURFACE_PRESSURE *groups
*add cond group_POINT_MOMENT *groups
*add cond group_VOLUME_ACCELERATION *groups
*loop groups *OnlyInCond
*set var numberloads(int)=Operation(numberloads(int)+1)
*end groups
*set var Counter = 0
*set cond group_POINT_LOAD *groups
*add cond group_POINT_MOMENT *groups
*loop groups *OnlyInCond       
*set var Counter=operation(Counter+1)
        {
        "help"            : "This process assigns a load value on conditions",	
        "kratos_module"   : "KratosMultiphysics.SolidMechanicsApplication",	
        "python_module"   : "assign_modulus_and_direction_to_conditions_process",
        "process_name"    : "AssignModulusAndDirectionToConditionsProcess",
        "Parameters"      : {
            "model_part_name" : "*GroupName",
            "variable_name"   : "*cond(Variable)",
*if(strcmp(cond(Time_Evolution),"INITIAL")==0)
	    "interval"        : [*GenData(Start_Time), *GenData(Start_Time)],
*elseif(strcmp(cond(Time_Evolution),"FULL")==0)
	    "interval"        : [*GenData(Start_Time), "End"],
*elseif(strcmp(cond(Time_Evolution),"INTERVAL")==0)
	    "interval"        : [*cond(Time_Interval,1), *cond(Time_Interval,2)],
*endif
*if(strcmp(cond(by_function),"True")==0 )
            "modulus"         : "*cond(Function)",
*else
            "modulus"         : *cond(Value),
*endif	    
            "direction"       : [*tcl(JoinByComma *cond(Direction))]
	    }
*if( Counter == numberloads )
        }
*else
	},
*endif
*end groups
*set cond group_LINE_LOAD *groups
*add cond group_SURFACE_LOAD *groups
*loop groups *OnlyInCond       
*set var Counter=operation(Counter+1)
        {
        "help"            : "This process assigns a load value on conditions",	
        "kratos_module"   : "KratosMultiphysics.SolidMechanicsApplication",	
        "python_module"   : "assign_modulus_and_direction_to_conditions_process",
        "process_name"    : "AssignModulusAndDirectionToConditionsProcess",
        "Parameters"      : {
            "model_part_name" : "*GroupName",
*if(strcmp(cond(by_function),"True")==0 )	    
            "variable_name"   : "*cond(Variable)_VECTOR",
*else
            "variable_name"   : "*cond(Variable)",
*endif	
*if(strcmp(cond(Time_Evolution),"INITIAL")==0)
	    "interval"        : [*GenData(Start_Time), *GenData(Start_Time)],
*elseif(strcmp(cond(Time_Evolution),"FULL")==0)
	    "interval"        : [*GenData(Start_Time), "End"],
*elseif(strcmp(cond(Time_Evolution),"INTERVAL")==0)
	    "interval"        : [*cond(Time_Interval,1), *cond(Time_Interval,2)],
*endif
*if(strcmp(cond(by_function),"True")==0 )
            "modulus"         : "*cond(Function)",
*else
            "modulus"         : *cond(Value),
*endif	    
            "direction"       : [*tcl(JoinByComma *cond(Direction))]
	    }
*if( Counter == numberloads )
        }
*else
	},
*endif
*end groups
*set cond group_VOLUME_ACCELERATION *groups
*loop groups *OnlyInCond       
*set var Counter=operation(Counter+1)
        {
        "help"            : "This process applies a volume acceleration",	
        "kratos_module"   : "KratosMultiphysics.SolidMechanicsApplication",
        "python_module"   : "assign_modulus_and_direction_to_nodes_process",
        "process_name"    : "AssignModulusAndDirectionToNodesProcess",
        "Parameters"      : {
            "model_part_name" : "*GroupName",
            "variable_name"   : "*cond(Variable)",
*if(strcmp(cond(Time_Evolution),"INITIAL")==0)
	    "interval"        : [*GenData(Start_Time), *GenData(Start_Time)],
*elseif(strcmp(cond(Time_Evolution),"FULL")==0)
	    "interval"        : [*GenData(Start_Time), "End"],
*elseif(strcmp(cond(Time_Evolution),"INTERVAL")==0)
	    "interval"        : [*cond(Time_Interval,1), *cond(Time_Interval,2)],
*endif
*if(strcmp(cond(by_function),"True")==0 )
            "modulus"         : "*cond(Function)",
*else
            "modulus"         : *cond(Value),
*endif
            "direction"       : [*tcl(JoinByComma *cond(Direction))]
	    }
*if( Counter == numberloads )
        }
*else
	},
*endif
*end groups
*set cond group_LINE_PRESSURE *groups
*add cond group_SURFACE_PRESSURE *groups
*loop groups *OnlyInCond    
*set var Counter=operation(Counter+1)
        }
        "help"            : "This process applies a pressure load",	
        "kratos_module"   : "KratosMultiphysics.SolidMechanicsApplication",
        "python_module"   : "assign_scalar_to_conditions_process",
        "process_name"    : "AssignScalarToConditionsProcess",
        "Parameters"      : {
            "model_part_name" : "*GroupName",
*if(strcmp(cond(by_function),"True")==0 )	    
            "variable_name"   : "*cond(Variable)_VECTOR",
*else
            "variable_name"   : "*cond(Variable)",
*endif	    
*if(strcmp(cond(Time_Evolution),"INITIAL")==0)
	    "interval"        : [*GenData(Start_Time), *GenData(Start_Time)],
*elseif(strcmp(cond(Time_Evolution),"FULL")==0)
	    "interval"        : [*GenData(Start_Time), "End"],
*elseif(strcmp(cond(Time_Evolution),"INTERVAL")==0)
	    "interval"        : [*cond(Time_Interval,1), *cond(Time_Interval,2)],
*endif
*if(strcmp(cond(by_function),"True")==0 )
            "value"         : "*cond(Function)"
*else
            "value"         : *cond(Value)
*endif
	    }
*if( Counter == numberloads )
        }
*else
	},
*endif
*end groups
    ],
    "output_process_list" : [{
        "help"            : "This process writes restart files",    
        "kratos_module"   : "KratosMultiphysics.SolidMechanicsApplication",    
        "python_module"   : "restart_process",
        "process_name"    : "RestartProcess",
        "Parameters"      : {
            "model_part_name"     : "Main_Domain",
            "save_restart"        : *tcl(string tolower *GenData(Print_Restart)),
            "restart_file_name"   : "*tcl(file tail [GiD_Info Project ModelName])",
            "restart_file_label"  : "step",
            "output_control_type" : "time",
            "output_frequency"    : *GenData(Restart_Frequency),
            "json_output"         : false
        }
    }],
    "output_configuration"     : {
        "result_file_configuration" : {
            "gidpost_flags"       : {
*if(strcmp(GenData(File_Format),"Binary")==0)
                "GiDPostMode"           : "GiD_PostBinary",
*elseif(strcmp(GenData(File_Format),"Ascii")==0)
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
            "node_output"         : *tcl(string tolower *GenData(Write_Particles)),
            "skin_output"         : false,
            "plane_output"        : [],
	    "nodal_results"       : [
*if(strcmp(GenData(Solver_Type),"DynamicSolver")==0)
                                      "VELOCITY",
				      "ACCELERATION",
*if(strcmp(GenData(DOFS),"U-W")==0)
                                      "WATER_DISPLACEMENT",
                                      "WATER_VELOCITY",
				      "WATER_ACCELERATION",
*elseif(strcmp(GenData(DOFS),"U-W-wP")==0)
                                      "WATER_DISPLACEMENT",
                                      "WATER_VELOCITY",
				      "WATER_ACCELERATION",
				      "WATER_PRESSURE",
				      "WATER_PRESSURE_VELOCITY",
				      "WATER_PRESSURE_ACCELERATION",
*elseif(strcmp(GenData(DOFS),"U-J-W-wP")==0)
                                      "JACOBIAN",
                                      "WATER_DISPLACEMENT",
                                      "WATER_VELOCITY",
				      "WATER_ACCELERATION",
				      "WATER_PRESSURE",
				      "WATER_PRESSURE_VELOCITY",
				      "WATER_PRESSURE_ACCELERATION",
*endif
*endif
*if(strcmp(GenData(Write_Reactions),"True")==0)
				      "DISPLACEMENT_REACTION",
*endif
*if(strcmp(GenData(Write_Contact_Forces),"True")==0)
				      "NORMAL",
				      "CONTACT_FORCE",
*endif
*if( strcmp(GenData(DOFS),"U-J")==0 )
						"JACOBIAN",
*endif
*if(strcmp(GenData(DOFS),"U-P")==0)
				      "PRESSURE",
*endif
*if(strcmp(GenData(DOFS),"U-J-wP")==0)
				      "WATER_PRESSURE",
				      "JACOBIAN",
*endif
*if(strcmp(GenData(DOFS),"U-wP")==0)
				      "WATER_PRESSURE",
*if(strcmp(GenData(Write_Contact_Forces),"True")==0)
				      "EFFECTIVE_CONTACT_FORCE",
				      "EFFECTIVE_CONTACT_STRESS",
*endif
*endif
				      "DISPLACEMENT"
				    ],
            "gauss_point_results" : [
*if(strcmp(GenData(Problem_Type),"mechanical")==0)
				      "CAUCHY_STRESS_TENSOR",
				      "GREEN_LAGRANGE_STRAIN_TENSOR",
*endif
*if(strcmp(GenData(DOFS),"U-W")==0)
				      "WATER_PRESSURE",
*endif
				      "STRESS_INV_P",
				      "STRESS_INV_J2",
				      "STRESS_INV_THETA",
				      "PLASTIC_STRAIN",
				      "DELTA_PLASTIC_STRAIN",
				      "PS",
				      "PM",
				      "PT",
                                      "INCR_SHEAR_PLASTIC"	
				    ],
	    "additional_list_files": [
*for(i=1;i<=GenData(List_Files,INT);i=i+1)
*if( i<GenData(List_Files,INT) )
			     	  *GenData(List_Files,*i),
*else
			     	  *GenData(List_Files,*i)
*endif
*end
	    			 ]
        },
        "point_data_configuration"  : []
    }

}

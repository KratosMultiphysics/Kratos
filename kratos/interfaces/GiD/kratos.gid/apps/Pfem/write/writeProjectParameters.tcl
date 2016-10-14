
# Project Parameters
proc Pfem::write::getParametersDict { } {
    Pfem::write::CalculateMyVariables
    set projectParametersDict [dict create]
    
    ##### Problem data #####
    # Create section
    set problemDataDict [GetPFEM_ProblemDataDict]
    # Add section to document
    dict set projectParametersDict problem_data $problemDataDict
    
    ##### solver_settings #####
    set solverSettingsDict [GetPFEM_SolverSettingsDict]
    dict set projectParametersDict solver_settings $solverSettingsDict
    
    ##### problem_process_list
    set problemProcessList [GetPFEM_ProblemProcessList]
    dict set projectParametersDict problem_process_list $problemProcessList
    
    ##### constraints_process_list
    dict set projectParametersDict constraints_process_list [write::getConditionsParametersDict SLNodalConditions "Nodal"]
    
    ##### loads_process_list
    dict set projectParametersDict loads_process_list [write::getConditionsParametersDict SLLoads]
    
    ##### output_configuration
    dict set projectParametersDict output_configuration [write::GetDefaultOutputDict]
    # restart options
    set restartDict [dict create ]
    dict set restartDict SaveRestart false
    dict set restartDict RestartFrequency 0
    dict set restartDict LoadRestart false
    dict set restartDict Restart_Step 0
    dict set projectParametersDict restart_options $restartDict
    
    # Constraints data
    set contraintsDict [dict create ]
    dict set contraintsDict incremental_load false
    dict set contraintsDict incremental_displacement false
    dict set projectParametersDict constraints_data $contraintsDict
    
    return $projectParametersDict
}
proc Pfem::write::GetPFEM_ProblemProcessList { } {
    set resultList [list ]
    lappend resultList [GetPFEM_RemeshDict]
    lappend resultList [GetPFEM_ContactDict]
    return $resultList
}

proc Pfem::write::GetPFEM_ContactDict { } {
    variable bodies_list
    set resultDict [dict create ]
    dict set resultDict "python_module" "contact_domain_process"
    dict set resultDict "kratos_module" "KratosMultiphysics.ContactMechanicsApplication"
    dict set resultDict "help" "This process applies contact domain search by remeshing outer boundaries"
    dict set resultDict "process_name" "ContactDomainProcess"
    
    set paramsDict [dict create ]
    dict set paramsDict "mesh_id" 0
    dict set paramsDict "model_part_name" "Main Domain"
    dict set paramsDict "meshing_control_type" "step"
    dict set paramsDict "meshing_frequency" 1.0
    dict set paramsDict "meshing_before_output" true
    
    foreach body $bodies_list {
        set bodyDict [dict create ]
        dict set bodyDict "python_module" "contact_domain"
        dict set bodyDict "model_part_name" [dict get $body body_name]
        dict set bodyDict "alpha_shape" 1.4
        dict set bodyDict "offset_factor" 0.0
        
            set meshing_strategyDict [dict create ]
            dict set meshing_strategyDict "python_module" "contact_meshing_strategy"
            dict set meshing_strategyDict "meshing_frequency" 0
            dict set meshing_strategyDict "remesh" true
            dict set meshing_strategyDict "constrained" false
            
                set contact_parametersDict [dict create ]
                dict set contact_parametersDict "contact_condition_type" "ContactDomainLM2DCondition"
                dict set contact_parametersDict "friction_law_type" "FrictionLaw"
                dict set contact_parametersDict "kratos_module" "KratosMultiphysics.ContactMechanicsApplication"
                
                    set variables_of_propertiesDict [dict create ]
                    dict set variables_of_propertiesDict "FRICTION_ACTIVE" false
                    dict set variables_of_propertiesDict "MU_STATIC" 0.3
                    dict set variables_of_propertiesDict "MU_DYNAMIC" 0.2
                    dict set variables_of_propertiesDict "PENALTY_PARAMETER" 1000
                    dict set variables_of_propertiesDict "TANGENTIAL_PENALTY_RATIO" 0.1
                    dict set variables_of_propertiesDict "TAU_STAB" 1
                dict set contact_parametersDict variables_of_properties $variables_of_propertiesDict
            dict set meshing_strategyDict contact_parameters $contact_parametersDict
        dict set bodyDict meshing_strategy $meshing_strategyDict
        lappend bodies $bodyDict
    }
    dict set paramsDict meshing_domains $bodies
    dict set resultDict Parameters $paramsDict
    return $resultDict
}

proc Pfem::write::GetPFEM_RemeshDict { } {
    variable bodies_list
    set resultDict [dict create ]
    dict set resultDict "python_module" "remesh_domains_process"
    dict set resultDict "kratos_module" "KratosMultiphysics.PfemBaseApplication"
    dict set resultDict "help" "This process applies meshing to the problem domains"
    dict set resultDict "process_name" "RemeshDomainsProcess"
    
    set paramsDict [dict create]
    dict set paramsDict "model_part_name" "Main Domain"
    dict set paramsDict "meshing_control_type" "step"
    dict set paramsDict "meshing_frequency" 1.0
    dict set paramsDict "meshing_before_output" true
    set meshing_domains_list [list ]
    foreach body $bodies_list {
        set bodyDict [dict create ]
        set body_name [dict get $body body_name]
        dict set bodyDict "python_module" "meshing_domain"
        dict set bodyDict "mesh_id" 1
        dict set bodyDict "model_part_name" $body_name
        dict set bodyDict "alpha_shape" 2.4
        dict set bodyDict "offset_factor" 0.0
        set remesh [write::getStringBinaryFromValue [Pfem::write::GetRemeshProperty $body_name "Remesh"]]
        set refine [write::getStringBinaryFromValue [Pfem::write::GetRemeshProperty $body_name "Refine"]]
        set meshing_strategyDict [dict create ]
        dict set meshing_strategyDict "python_module" "meshing_strategy"
        dict set meshing_strategyDict "meshing_frequency" 0
        dict set meshing_strategyDict "remesh" $remesh
        dict set meshing_strategyDict "refine" $refine
        dict set meshing_strategyDict "reconnect" false
        dict set meshing_strategyDict "transfer" false
        dict set meshing_strategyDict "constrained" false
        dict set meshing_strategyDict "mesh_smoothing" false
        dict set meshing_strategyDict "variables_smoothing" false
        dict set meshing_strategyDict "elemental_variables_to_smooth" [list "DETERMINANT_F" ]
        dict set meshing_strategyDict "reference_element_type" "Element2D3N"
        dict set meshing_strategyDict "reference_condition_type" "CompositeCondition2D2N"
        dict set bodyDict meshing_strategy $meshing_strategyDict
        
        set spatial_bounding_boxDict [dict create ]
        set upX [expr 0.0]; set upY [expr 0.0]; set upZ [expr 0.0]
        dict set spatial_bounding_boxDict "upper_point" [list $upX $upY $upZ]
        set lpX [expr 0.0]; set lpY [expr 0.0]; set lpZ [expr 0.0]
        dict set spatial_bounding_boxDict "lower_point" [list $lpX $lpY $lpZ]
        set vlX [expr 0.0]; set vlY [expr 0.0]; set vlZ [expr 0.0]
        dict set spatial_bounding_boxDict "velocity" [list $vlX $vlY $vlZ]
        dict set bodyDict spatial_bounding_box $spatial_bounding_boxDict
        
        set refining_parametersDict [dict create ]
        dict set refining_parametersDict "critical_size" 0.0
        dict set refining_parametersDict "threshold_variable" "PLASTIC_STRAIN"
        dict set refining_parametersDict "reference_threshold" 0.0
        dict set refining_parametersDict "error_variable" "NORM_ISOCHORIC_STRESS"
        dict set refining_parametersDict "reference_error" 0.0
        dict set refining_parametersDict "add_nodes" true
        dict set refining_parametersDict "insert_nodes" false
        
        set remove_nodesDict [dict create]
        dict set remove_nodesDict "apply_removal" false
        dict set remove_nodesDict "on_distance" false
        dict set remove_nodesDict "on_threshold" false
        dict set remove_nodesDict "on_error" false
        dict set refining_parametersDict remove_nodes $remove_nodesDict
        
        set remove_boundaryDict [dict create]
        dict set remove_boundaryDict "apply_removal" false
        dict set remove_boundaryDict "on_distance" false
        dict set remove_boundaryDict "on_threshold" false
        dict set remove_boundaryDict "on_error" false
        dict set refining_parametersDict remove_boundary $remove_boundaryDict
        
        set refine_elementsDict [dict create]
        dict set refine_elementsDict "apply_removal" false
        dict set refine_elementsDict "on_distance" false
        dict set refine_elementsDict "on_threshold" false
        dict set refine_elementsDict "on_error" false
        dict set refining_parametersDict refine_elements $refine_elementsDict
        
        set refine_boundaryDict [dict create]
        dict set refine_boundaryDict "apply_removal" false
        dict set refine_boundaryDict "on_distance" false
        dict set refine_boundaryDict "on_threshold" false
        dict set refine_boundaryDict "on_error" false
        dict set refining_parametersDict refine_boundary $refine_boundaryDict
        
        set refining_boxDict [dict create]
        dict set refining_boxDict "refine_in_box_only" false
        set upX [expr 0.0]; set upY [expr 0.0]; set upZ [expr 0.0]
        dict set refining_boxDict "upper_point" [list $upX $upY $upZ]
        set lpX [expr 0.0]; set lpY [expr 0.0]; set lpZ [expr 0.0]
        dict set refining_boxDict "lower_point" [list $lpX $lpY $lpZ]
        set vlX [expr 0.0]; set vlY [expr 0.0]; set vlZ [expr 0.0]
        dict set refining_boxDict "velocity" [list $vlX $vlY $vlZ]
        dict set refining_parametersDict refining_box $refining_boxDict
        
        dict set bodyDict refining_parameters $refining_parametersDict
        
        dict set bodyDict "elemental_variables_to_transfer" [list "CAUCHY_STRESS_VECTOR" "DEFORMATION_GRADIENT"]
        lappend meshing_domains_list $bodyDict
    }
    dict set paramsDict meshing_domains $meshing_domains_list
    dict set resultDict Parameters $paramsDict
    return $resultDict
}

proc Pfem::write::GetRemeshProperty { body_name property } {
    set ret ""
    set root [customlib::GetBaseRoot]
    set xp1 "[spdAux::getRoute "PFEM_Bodies"]/blockdata"
    set remesh_name ""
    foreach body_node [$root selectNodes $xp1] {
        if {[$body_node @name] eq $body_name} {
            set remesh_name [get_domnode_attribute [$body_node selectNodes ".//value\[@n='MeshingStrategy'\]"] v]
            break
        }
    }
    if {$remesh_name ne ""} {
        variable remesh_domains_dict
        if {[dict exists $remesh_domains_dict ${remesh_name} $property]} {
            set ret [dict get $remesh_domains_dict ${remesh_name} $property]
        }
    }
    if {$ret eq ""} {set ret false}
    return $ret
}

proc Pfem::write::GetPFEM_ProblemDataDict { } {
    set problemDataDict [dict create]
    dict set problemDataDict problem_name [file tail [GiD_Info Project ModelName]]

    dict set problemDataDict model_part_name "Main Domain"
    set nDim $::Model::SpatialDimension
    set nDim [expr [string range [write::getValue nDim] 0 0] ]
    dict set problemDataDict domain_size $nDim
   
    dict set problemDataDict time_step [write::getValue PFEM_TimeParameters DeltaTime]
    dict set problemDataDict start_time [write::getValue PFEM_TimeParameters StartTime]
    dict set problemDataDict end_time [write::getValue PFEM_TimeParameters EndTime]
    dict set problemDataDict echo_level [write::getValue Results EchoLevel]
    dict set problemDataDict threads 1
    
    return $problemDataDict
}

proc Pfem::write::GetPFEM_SolverSettingsDict { } {
    variable bodies_list
    set solverSettingsDict [dict create]
    set currentStrategyId [write::getValue PFEM_SolStrat]
    set strategy_write_name [[::Model::GetSolutionStrategy $currentStrategyId] getAttribute "ImplementedInPythonFile"]
    dict set solverSettingsDict solver_type $strategy_write_name
    #~ dict set solverSettingsDict domain_size [expr $nDim]
    dict set solverSettingsDict echo_level [write::getValue Results EchoLevel]
    dict set solverSettingsDict solution_type [write::getValue PFEM_SolutionType]
    
    dict set solverSettingsDict time_integration_method [write::getValue PFEM_SolStrat]
    dict set solverSettingsDict scheme_type [write::getValue PFEM_Scheme]
    
    # model import settings
    set modelDict [dict create]
    dict set modelDict input_type "mdpa"
    dict set modelDict input_filename [file tail [GiD_Info Project ModelName]]
    dict set modelDict input_file_label 0
    dict set solverSettingsDict model_import_settings $modelDict
    
    # Solution strategy parameters and Solvers
    set solverSettingsDict [dict merge $solverSettingsDict [write::getSolutionStrategyParametersDict] ]
    set solverSettingsDict [dict merge $solverSettingsDict [write::getSolversParametersDict Pfem] ]
    
    set listsubmodelparts [list ]
    foreach part_un [Pfem::write::GetPartsUN] {
        lappend listsubmodelparts {*}[write::getSubModelPartNames $part_un]
    }
    dict set solverSettingsDict bodies_list $bodies_list
    dict set solverSettingsDict problem_domain_sub_model_part_list $listsubmodelparts
    dict set solverSettingsDict processes_sub_model_part_list [write::getSubModelPartNames "PFEM_NodalConditions" "PFEM_Loads"]
    
    return $solverSettingsDict
}

proc Pfem::write::ProcessBodiesList { } {
    customlib::UpdateDocument
    set bodiesList [list ]
    set root [customlib::GetBaseRoot]
    set xp1 "[spdAux::getRoute "PFEM_Bodies"]/blockdata"
    foreach body_node [$root selectNodes $xp1] {
        set body [dict create]
        set name [$body_node @name]
        set body_type_path ".//value\[@n='BodyType'\]"
        set body_type [get_domnode_attribute [$body_node selectNodes $body_type_path] v]
        set parts [list ]
        foreach part_node [$body_node selectNodes "./condition/group"] {
            lappend parts [write::getMeshId "Parts" [$part_node @n]]
        }
        dict set body "body_type" $body_type
        dict set body "body_name" $name
        dict set body "parts_list" $parts
        lappend bodiesList $body
    }
    return $bodiesList
}

proc Pfem::write::GetNodalDataDict { } {
    set doc $gid_groups_conds::doc
    set root [$doc documentElement]
    set NodalData [list ]
    set parts [list "PFEM_Rigid2DParts" "PFEM_Rigid3DParts" "PFEM_Deformable2DParts" "PFEM_Deformable3DParts" "PFEM_Fluid2DParts" "PFEM_Fluid3DParts"]
    
    foreach part $parts {
        set xp1 "[spdAux::getRoute $part]/group"
        set groups [$root selectNodes $xp1]
        foreach group $groups {
            set partid [[$group parent] @n]
            set groupid [$group @n]
            set processDict [dict create]
            dict set processDict process_name "ApplyValuesToNodes"
            dict set processDict kratos_module "KratosMultiphysics.PFEMBaseApplication"
            
            set params [dict create]
            set xp2 "./value"
            set atts [$group selectNodes $xp2]
            #W "$group $groupid $atts"
            foreach att $atts {
                set state [get_domnode_attribute $att state]
                if {$state ne "hidden"} {
                    set paramName [$att @n]
                    set paramValue [get_domnode_attribute $att v]
                    if {$paramName eq "Material"} {
                        set matdict [::write::getAllMaterialParametersDict $paramValue]
                        dict set matdict Name $paramValue
                        dict set params $paramName $matdict
                    } {
                        if {[write::isBoolean $paramValue]} {set paramValue [expr $paramValue]}
                        dict set params $paramName $paramValue
                    }
                }
            }
            dict set params "model_part_name" [::write::getMeshId $partid $groupid]
            dict set processDict "Parameters" $params
            lappend NodalData $processDict
        }
    }
    
    return $NodalData
}

proc Pfem::write::ProcessRemeshDomainsDict { } {
    customlib::UpdateDocument
    set domains_dict [dict create ]
    set root [customlib::GetBaseRoot]
    set xp1 "[spdAux::getRoute "PFEM_meshing_domains"]/blockdata"
    foreach domain_node [$root selectNodes $xp1] {
        set name [$domain_node @name]
        foreach part_node [$domain_node selectNodes "./value"] {
            dict set domains_dict $name [get_domnode_attribute $part_node n] [get_domnode_attribute $part_node v]
        }
    }
    return $domains_dict
}

proc Pfem::write::CalculateMyVariables { } {
    variable bodies_list
    set bodies_list [Pfem::write::ProcessBodiesList]
    variable remesh_domains_dict
    set remesh_domains_dict [Pfem::write::ProcessRemeshDomainsDict]
 }
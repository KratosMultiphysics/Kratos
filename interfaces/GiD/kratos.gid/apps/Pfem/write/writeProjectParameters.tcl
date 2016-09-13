
# Project Parameters
proc Pfem::write::getParametersDict { } {
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

    return $resultList
}

proc Pfem::write::GetPFEM_ProblemDataDict { } {
    set problemDataDict [dict create]
    dict set problemDataDict problem_name [file tail [GiD_Info Project ModelName]]
    
    set nDim $::Model::SpatialDimension
    dict set problemDataDict domain_size $nDim
    
    dict set problemDataDict time_step [write::getValue PFEM_TimeParameters DeltaTime]
    dict set problemDataDict start_time [write::getValue PFEM_TimeParameters StartTime]
    dict set problemDataDict end_time [write::getValue PFEM_TimeParameters EndTime]
    dict set problemDataDict echo_level [write::getValue Results EchoLevel]
    
    return $problemDataDict
}

proc Pfem::write::GetPFEM_SolverSettingsDict { } {
    set solverSettingsDict [dict create]
    set currentStrategyId [write::getValue PFEM_SolStrat]
    set strategy_write_name [[::Model::GetSolutionStrategy $currentStrategyId] getAttribute "ImplementedInPythonFile"]
    dict set solverSettingsDict solver_type $strategy_write_name
    #~ dict set solverSettingsDict domain_size [expr $nDim]
    dict set solverSettingsDict echo_level [write::getValue Results EchoLevel]
    dict set solverSettingsDict solution_type [write::getValue PFEM_SoluType]
    
    dict set solverSettingsDict time_integration_method [write::getValue PFEM_SolStrat]
    dict set solverSettingsDict scheme_type [write::getValue PFEM_Scheme]
    
    # model import settings
    set modelDict [dict create]
    dict set modelDict input_type "mdpa"
    dict set modelDict input_filename [file tail [GiD_Info Project ModelName]]
    dict set solverSettingsDict model_import_settings $modelDict
    
    # Solution strategy parameters and Solvers
    set solverSettingsDict [dict merge $solverSettingsDict [write::getSolutionStrategyParametersDict] ]
    set solverSettingsDict [dict merge $solverSettingsDict [write::getSolversParametersDict Pfem] ]
    
    set listsubmodelparts [list ]
    foreach part_un [Pfem::write::GetPartsUN] {
        lappend listsubmodelparts {*}[write::getSubModelPartNames $part_un]
    }
    dict set solverSettingsDict problem_domain_sub_model_part_list $listsubmodelparts
    dict set solverSettingsDict bodies_list [GetBodiesList]
    dict set solverSettingsDict processes_sub_model_part_list [write::getSubModelPartNames "PFEM_NodalConditions" "PFEM_Loads"]
    
    return $solverSettingsDict
}

proc Pfem::write::GetBodiesList { } {
    customlib::UpdateDocument
    set bodiesList [list ]
    set root [customlib::GetBaseRoot]
    set xp1 "[spdAux::getRoute "PFEM_Bodies"]/blockdata"
    foreach body_node [$root selectNodes $xp1] {
        set body [dict create]
        set name [$body_node @name]
        set parts [list ]
        foreach part_node [$body_node selectNodes "./condition/group"] {
            lappend parts [write::getMeshId "Parts" [$part_node @n]]
        }
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
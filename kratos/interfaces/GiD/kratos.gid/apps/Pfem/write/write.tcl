namespace eval Pfem::write {
}

proc Pfem::write::Init { } {
    Solid::write::AddValidApps "Pfem"
}

# Project Parameters
proc Pfem::write::getParametersDict { } {
    set projectParametersDict [dict create]
    
    ##### Problem data #####
    # Create section
    set problemDataDict [dict create]
    
    # Add items to section
    set model_name [file tail [GiD_Info Project ModelName]]
    dict set problemDataDict problem_name $model_name
    
    set nDim $::Model::SpatialDimension
    dict set problemDataDict domain_size $nDim
    
    dict set problemDataDict time_step [write::getValue PFEM_TimeParameters DeltaTime]
    dict set problemDataDict start_time [write::getValue PFEM_TimeParameters StartTime]
    dict set problemDataDict end_time [write::getValue PFEM_TimeParameters EndTime]
    set echo_level [write::getValue Results EchoLevel]
    dict set problemDataDict echo_level $echo_level

    # Add section to document
    dict set projectParametersDict problem_data $problemDataDict
    
    ##### solver_settings #####
    set solverSettingsDict [dict create]
    set currentStrategyId [write::getValue PFEM_SolStrat]
    set strategy_write_name [[::Model::GetSolutionStrategy $currentStrategyId] getAttribute "ImplementedInPythonFile"]
    dict set solverSettingsDict solver_type $strategy_write_name
    #~ dict set solverSettingsDict domain_size [expr $nDim]
    dict set solverSettingsDict echo_level $echo_level
    dict set solverSettingsDict solution_type [write::getValue PFEM_SoluType]
    
    dict set solverSettingsDict time_integration_method [write::getValue PFEM_SolStrat]
    dict set solverSettingsDict scheme_type [write::getValue PFEM_Scheme]
    
    
    # model import settings
    set modelDict [dict create]
    dict set modelDict input_type "mdpa"
    dict set modelDict input_filename $model_name
    dict set solverSettingsDict model_import_settings $modelDict
    
    # Solution strategy parameters and Solvers
    set solverSettingsDict [dict merge $solverSettingsDict [write::getSolutionStrategyParametersDict] ]
    set solverSettingsDict [dict merge $solverSettingsDict [write::getSolversParametersDict Pfem] ]
    
    dict set solverSettingsDict problem_domain_sub_model_part_list [write::getSubModelPartNames "PFEM_Parts"]
    dict set solverSettingsDict processes_sub_model_part_list [write::getSubModelPartNames "PFEM_NodalConditions" "PFEM_Loads"]
    
    dict set projectParametersDict solver_settings $solverSettingsDict
    
    ##### problem_process_list
    
    ##### constraints_process_list
    ##### loads_process_list
    ##### solver_settings
    ##### output_configuration
    
    return $projectParametersDict
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
proc Pfem::write::writeParametersEvent { } {
    write::WriteJSON [getParametersDict]
}

# Model Part Blocks
proc Pfem::write::writeModelPartEvent { } {
    Solid::write::writeModelPartEvent
}


# Custom files (Copy python scripts, write materials file...)
proc Pfem::write::writeCustomFilesEvent { } {

}


Pfem::write::Init

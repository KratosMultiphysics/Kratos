# Project Parameters
proc Fluid::write::writeParametersEvent { } {
    variable BCUN
    set projectParametersDict [dict create]
    
    # First section -> Problem data
    set problemDataDict [dict create]
    set model_name [file tail [GiD_Info Project ModelName]]
    dict set problemDataDict problem_name $model_name
    dict set problemDataDict model_part_name "MainModelPart"
    set nDim [expr [string range [write::getValue nDim] 0 0]]
    dict set problemDataDict domain_size $nDim
    
   
    # Time Parameters
    dict set problemDataDict start_step [write::getValue FLTimeParameters StartTime] 
    dict set problemDataDict end_time [write::getValue FLTimeParameters EndTime]
    dict set problemDataDict time_step [write::getValue FLTimeParameters DeltaTime]
    #dict set problemDataDict divergence_step [expr [write::getValue FLTimeParameters DivergenceCleareanceStep]]
    
    dict set projectParametersDict problem_data $problemDataDict
    
    # output configuration
    dict set projectParametersDict output_configuration [write::GetDefaultOutputDict]
    
    # restart options
    set restartDict [dict create]
    dict set restartDict SaveRestart False
    dict set restartDict RestartFrequency 0
    dict set restartDict LoadRestart False
    dict set restartDict Restart_Step 0
    dict set projectParametersDict restart_options $restartDict
    
    # Solver settings
    set solverSettingsDict [dict create]
    set currentStrategyId [write::getValue FLSolStrat]
    set strategy_write_name [[::Model::GetSolutionStrategy $currentStrategyId] getAttribute "ImplementedInPythonFile"]
    dict set solverSettingsDict solver_type $strategy_write_name
    
    # model import settings
    set modelDict [dict create]
    dict set modelDict input_type "mdpa"
    dict set modelDict input_filename $model_name
    dict set solverSettingsDict model_import_settings $modelDict
    
    set solverSettingsDict [dict merge $solverSettingsDict [write::getSolutionStrategyParametersDict] ]
    set solverSettingsDict [dict merge $solverSettingsDict [write::getSolversParametersDict] ]
    # Parts
    dict set solverSettingsDict volume_model_part_name {*}[write::getPartsMeshId]
    # Skin parts
    dict set solverSettingsDict skin_parts [getBoundaryConditionMeshId]
        
    dict set projectParametersDict solver_settings $solverSettingsDict
        
    # Boundary conditions processes
    dict set projectParametersDict boundary_conditions_process_list [write::getConditionsParametersDict $BCUN]
    dict set projectParametersDict gravity [getGravityProcessDict]
    
    write::WriteJSON $projectParametersDict
}

# Skin SubModelParts ids
proc Fluid::write::getGravityProcessDict {} {
    set doc $gid_groups_conds::doc
    set root [$doc documentElement]
    
    set xp1 [spdAux::getRoute "FLGravity"]
    set value [get_domnode_attribute [$root selectNodes "$xp1/value\[@n='GravityValue'\]"] v]
    set cx [get_domnode_attribute [$root selectNodes "$xp1/value\[@n='Cx'\]"] v]
    set cy [get_domnode_attribute [$root selectNodes "$xp1/value\[@n='Cy'\]"] v]
    set cz [get_domnode_attribute [$root selectNodes "$xp1/value\[@n='Cz'\]"] v]
    #W "Gravity $value on \[$cx , $cy , $cz\]"
    set pdict [dict create]
    dict set pdict "implemented_in_file" "apply_gravity_process"
    dict set pdict "implemented_in_module" "KratosMultiphysics.FluidDynamicsApplication"
    dict set pdict "process_name" "ApplyGravity"
    set params [dict create]
    dict set params "mesh_id" 0
    set partgroup [write::getPartsMeshId]
    dict set params "model_part_name" $partgroup
    dict set params "value" $value
    dict set params "vector" [list $cx $cy $cz]
    dict set pdict "Parameters" $params
    
    return $pdict
}

# Skin SubModelParts ids
proc Fluid::write::getBoundaryConditionMeshId {} {
    variable BCUN
    set doc $gid_groups_conds::doc
    set root [$doc documentElement]
    
    set listOfBCGroups [list ]
    set xp1 "[spdAux::getRoute $BCUN]/condition/group"
    set groups [$root selectNodes $xp1]    
    foreach group $groups {
        set groupName [$group @n]
        set cid [[$group parent] @n]
        set gname [::write::getMeshId $cid $groupName]
        if {$gname ni $listOfBCGroups} {lappend listOfBCGroups $gname}
    }
    
    return $listOfBCGroups
}


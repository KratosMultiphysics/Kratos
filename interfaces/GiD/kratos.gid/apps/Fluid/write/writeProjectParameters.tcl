# Project Parameters
proc Fluid::write::writeParametersEvent { } {
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
    set outputConfigDict [dict create]
    dict set outputConfigDict output_filename "[file tail [GiD_Info Project ModelName]].out"
    dict set outputConfigDict gid_post_mode [write::getValue FLResults GiDPostMode]
    dict set outputConfigDict gid_multi_file_flag [write::getValue FLResults GiDMultiFileFlag]
    dict set outputConfigDict gid_write_mesh_flag True
    dict set outputConfigDict gid_write_conditions_flag True
    dict set outputConfigDict gid_write_particles_flag False
    dict set outputConfigDict gid_write_frequency [write::getValue FLResults OutputDeltaTime]
    dict set outputConfigDict plot_graphs False
    dict set outputConfigDict plot_frequency 0
    dict set outputConfigDict print_lists True
    dict set outputConfigDict output_time [write::getValue FLResults OutputDeltaTime]
    dict set outputConfigDict volume_output [write::getValue FLResults VolumeOutput]
    dict set outputConfigDict add_skin true
    
    set cut_list [list ]
    set normal [list 0.0 1.0 0.0]
    set point [list 0.0 0.0 0.0]
    lappend cut_list [dict create normal $normal point $point]
    
    dict set outputConfigDict cut_planes [dict create cut_list $cut_list]

    # on nodes
    dict set outputConfigDict nodal_results [write::GetResultsList "FLNodalResults"]
    # on elements
    dict set projectParametersDict output_configuration [write::GetResultsList "FLElementResults"]
    
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
    variable BCUN
    dict set projectParametersDict boundary_conditions_process_list [write::getConditionsParametersDict $BCUN]
    
    write::WriteJSON $projectParametersDict
}

# Skin SubModelParts ids
proc Fluid::write::getBoundaryConditionMeshId {} {
    variable BCUN
    set doc $gid_groups_conds::doc
    set root [$doc documentElement]
    
    set listOfBCGroups [list ]
    set xp1 "[apps::getRoute $BCUN]/condition/group"
    set groups [$root selectNodes $xp1]    
    foreach group $groups {
        set groupName [$group @n]
        set cid [[$group parent] @n]
        set gname [::write::getMeshId $cid $groupName]
        if {$gname ni $listOfBCGroups} {lappend listOfBCGroups $gname}
    }
    
    return $listOfBCGroups
}


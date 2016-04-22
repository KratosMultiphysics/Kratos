# Project Parameters
proc Solid::write::writeParametersEvent { } {
    set projectParametersDict [dict create]
    
    # Problem data
    # Create section
    set problemDataDict [dict create]
    
    # Add items to section
    set model_name [file tail [GiD_Info Project ModelName]]
    dict set problemDataDict problem_name $model_name
    dict set problemDataDict model_part_name "Structure"
    set nDim [expr [string range [write::getValue nDim] 0 0] ]
    dict set problemDataDict domain_size $nDim
        
    # Parallelization
    set paralleltype [write::getValue ParallelType]
    if {$paralleltype eq "OpenMP"} {
        set nthreads [write::getValue Parallelization OpenMPNumberOfThreads]
        dict set problemDataDict NumberofThreads $nthreads
    } else {
        set nthreads [write::getValue Parallelization MPINumberOfProcessors]
        dict set problemDataDict NumberofProcessors $nthreads
    }
    
    # Time Parameters
    dict set problemDataDict time_step [write::getValue SLTimeParameters DeltaTime]
    dict set problemDataDict start_time [write::getValue SLTimeParameters StartTime]
    dict set problemDataDict end_time [write::getValue SLTimeParameters EndTime]
    set echo_level [write::getValue Results EchoLevel]
    dict set problemDataDict echo_level $echo_level
    
    # Add section to document
    dict set projectParametersDict problem_data $problemDataDict
    
    
    # Solution strategy
    set solverSettingsDict [dict create]
    set currentStrategyId [write::getValue SLSolStrat]
    set strategy_write_name [[::Model::GetSolutionStrategy $currentStrategyId] getAttribute "ImplementedInPythonFile"]
    dict set solverSettingsDict solver_type $strategy_write_name
    #~ dict set solverSettingsDict domain_size [expr $nDim]
    dict set solverSettingsDict echo_level $echo_level
    dict set solverSettingsDict solution_type [write::getValue SLSoluType]
    set solutiontype [write::getValue SLSoluType]
    if {$solutiontype eq "Static"} {
        dict set solverSettingsDict analysis_type [write::getValue SLAnalysisType]
    } elseif {$solutiontype eq "Dynamic"} {
        dict set solverSettingsDict time_integration_method [write::getValue SLSolStrat]
        dict set solverSettingsDict scheme_type [write::getValue SLScheme]
    }
    
    # model import settings
    set modelDict [dict create]
    dict set modelDict input_type "mdpa"
    dict set modelDict input_filename $model_name
    dict set solverSettingsDict model_import_settings $modelDict
    
    # Solution strategy parameters and Solvers
    set solverSettingsDict [dict merge $solverSettingsDict [write::getSolutionStrategyParametersDict] ]
    set solverSettingsDict [dict merge $solverSettingsDict [write::getSolversParametersDict] ]
    
    dict set solverSettingsDict problem_domain_sub_model_part_list [getSubModelPartNames "SLParts"]
    dict set solverSettingsDict processes_sub_model_part_list [getSubModelPartNames "SLNodalConditions" "SLLoads"]
    
    dict set projectParametersDict solver_settings $solverSettingsDict
    
    # Lists of processes
    dict set projectParametersDict constraints_process_list [write::getConditionsParametersDict SLNodalConditions "Nodal"]
    
    dict set projectParametersDict loads_process_list [write::getConditionsParametersDict SLLoads]

    # GiD output configuration
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
    
    write::WriteJSON $projectParametersDict
}

proc Solid::write::getSubModelPartNames { args } {
    set doc $gid_groups_conds::doc
    set root [$doc documentElement]
    
    set listOfProcessedGroups [list ]
    set groups [list ]
    foreach un $args {
        set xp1 "[spdAux::getRoute $un]/condition/group"
        set xp2 "[spdAux::getRoute $un]/group"
        set grs [$root selectNodes $xp1]
        if {$grs ne ""} {lappend groups {*}$grs}
        set grs [$root selectNodes $xp2]
        if {$grs ne ""} {lappend groups {*}$grs}
    }
    foreach group $groups {
        set groupName [$group @n]
        set cid [[$group parent] @n]
        set gname [::write::getMeshId $cid $groupName]
        if {$gname ni $listOfProcessedGroups} {lappend listOfProcessedGroups $gname}
    }
    
    return $listOfProcessedGroups
}

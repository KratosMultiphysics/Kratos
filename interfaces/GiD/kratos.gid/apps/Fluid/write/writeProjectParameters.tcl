# Project Parameters
proc ::Fluid::write::getParametersDict { } {
    variable BCUN

    set projectParametersDict [dict create]

    # First section -> Problem data
    set problemDataDict [dict create]
    set model_name [file tail [GiD_Info Project ModelName]]
    dict set problemDataDict problem_name $model_name
    dict set problemDataDict model_part_name "MainModelPart"
    set nDim [expr [string range [write::getValue nDim] 0 0]]
    dict set problemDataDict domain_size $nDim

    # Parallelization
    set paralleltype [write::getValue ParallelType]
    dict set problemDataDict "parallel_type" $paralleltype
    if {$paralleltype eq "OpenMP"} {
        #set nthreads [write::getValue Parallelization OpenMPNumberOfThreads]
        #dict set problemDataDict NumberofThreads $nthreads
    } else {
        #set nthreads [write::getValue Parallelization MPINumberOfProcessors]
        #dict set problemDataDict NumberofProcessors $nthreads
    }

    # Write the echo level in the problem data section
    set echo_level [write::getValue Results EchoLevel]
    dict set problemDataDict echo_level $echo_level

    # Time Parameters
    dict set problemDataDict start_time [write::getValue FLTimeParameters StartTime]
    dict set problemDataDict end_time [write::getValue FLTimeParameters EndTime]
    set automaticDeltaTime [write::getValue FLTimeParameters AutomaticDeltaTime]
    #if {$automaticDeltaTime eq "Yes"} {
    #    dict set problemDataDict "CFL_number" [write::getValue FLTimeParameters CFLNumber]
    #} else {
    #    dict set problemDataDict "time_step" [write::getValue FLTimeParameters DeltaTime]
    #}

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
    set solverSettingsDict [dict merge $solverSettingsDict [write::getSolversParametersDict Fluid] ]
    # Parts
    dict set solverSettingsDict volume_model_part_name {*}[write::getPartsMeshId]
    # Skin parts
    dict set solverSettingsDict skin_parts [getBoundaryConditionMeshId]
    # No skin parts
    dict set solverSettingsDict no_skin_parts [getNoSkinConditionMeshId]
    # Time stepping settings
    set timeSteppingDict [dict create]
    dict set timeSteppingDict automatic_time_step $automaticDeltaTime
    if {$automaticDeltaTime eq "Yes"} {
        dict set timeSteppingDict "CFL_number" [write::getValue FLTimeParameters CFLNumber]
        dict set timeSteppingDict "maximum_delta_time" [write::getValue FLTimeParameters MaximumDeltaTime]
    } else {
        dict set timeSteppingDict "time_step" [write::getValue FLTimeParameters DeltaTime]
    }
    dict set solverSettingsDict time_stepping $timeSteppingDict

    dict set projectParametersDict solver_settings $solverSettingsDict

    # Boundary conditions processes
    dict set projectParametersDict initial_conditions_process_list [write::getConditionsParametersDict "FLNodalConditions" "Nodal"]
    dict set projectParametersDict boundary_conditions_process_list [write::getConditionsParametersDict $BCUN]
    dict set projectParametersDict gravity [list [getGravityProcessDict] ]
    dict set projectParametersDict auxiliar_process_list [getAuxiliarProcessList]

    return $projectParametersDict
}

proc Fluid::write::writeParametersEvent { } {
    set projectParametersDict [getParametersDict]
    write::SetParallelismConfiguration
    write::WriteJSON $projectParametersDict
}

proc Fluid::write::getAuxiliarProcessList {} {
    set process_list [list ]

    foreach process [getDragProcessList] {lappend process_list $process}

    return $process_list
}

proc Fluid::write::getDragProcessList {} {
    set root [customlib::GetBaseRoot]

    set process_list [list ]
    set xp1 "[spdAux::getRoute FLDrags]/group"
    set groups [$root selectNodes $xp1]
    foreach group $groups {
        set groupName [$group @n]
        set cid [[$group parent] @n]
        set submodelpart [::write::getMeshId $cid $groupName]

        set write_output [write::getStringBinaryFromValue [write::getValueByNode [$group selectNodes "./value\[@n='write_drag_output_file'\]"]]]
        set print_screen [write::getStringBinaryFromValue [write::getValueByNode [$group selectNodes "./value\[@n='print_drag_to_screen'\]"]]]
        set interval_name [write::getValueByNode [$group selectNodes "./value\[@n='Interval'\]"]]

        set pdict [dict create]
        dict set pdict "python_module" "compute_drag_process"
        dict set pdict "kratos_module" "KratosMultiphysics.FluidDynamicsApplication"
        dict set pdict "process_name" "ComputeDragProcess"
        set params [dict create]
        dict set params "mesh_id" 0
        dict set params "model_part_name" $submodelpart
        dict set params "write_drag_output_file" $write_output
        dict set params "print_drag_to_screen" $print_screen
        dict set params "interval" [write::getInterval $interval_name]
        dict set pdict "Parameters" $params

        lappend process_list $pdict
    }

    return $process_list
}

# Gravity SubModelParts and Process collection
proc Fluid::write::getGravityProcessDict {} {
    set root [customlib::GetBaseRoot]

    set value [write::getValue FLGravity GravityValue]
    set cx [write::getValue FLGravity Cx]
    set cy [write::getValue FLGravity Cy]
    set cz [write::getValue FLGravity Cz]
    #W "Gravity $value on \[$cx , $cy , $cz\]"
    set pdict [dict create]
    dict set pdict "python_module" "process_factory"
    dict set pdict "kratos_module" "KratosMultiphysics"
    dict set pdict "process_name" "ApplyConstantVectorValueProcess"
    set params [dict create]
    dict set params "mesh_id" 0
    set partgroup [write::getPartsMeshId]
    dict set params "model_part_name" [concat [lindex $partgroup 0]]
    dict set params "variable_name" "BODY_FORCE"
    dict set params "modulus" $value
    dict set params "direction" [list $cx $cy $cz]
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
        set cond [Model::getCondition $cid]
        if {[$cond getAttribute "SkinConditions"] eq "True"} {
            set gname [::write::getMeshId $cid $groupName]
            if {$gname ni $listOfBCGroups} {lappend listOfBCGroups $gname}
        }
    }

    return $listOfBCGroups
}

# No-skin SubModelParts ids
proc Fluid::write::getNoSkinConditionMeshId {} {
    variable BCUN
    set doc $gid_groups_conds::doc
    set root [$doc documentElement]

    set listOfNoSkinGroups [list ]

    # Append drag processes model parts names
    set xp1 "[spdAux::getRoute FLDrags]/group"
    set dragGroups [$root selectNodes $xp1]
    foreach dragGroup $dragGroups {
        set groupName [$dragGroup @n]
        set cid [[$dragGroup parent] @n]
        set submodelpart [::write::getMeshId $cid $groupName]
        if {$submodelpart ni $listOfNoSkinGroups} {lappend listOfNoSkinGroups $submodelpart}
    }

    # Append no skin conditions model parts names
    set xp1 "[spdAux::getRoute $BCUN]/condition/group"
    set groups [$root selectNodes $xp1]
    foreach group $groups {
        set groupName [$group @n]
        set cid [[$group parent] @n]
        set cond [Model::getCondition $cid]
        if {[$cond getAttribute "SkinConditions"] eq "False"} {
            set gname [::write::getMeshId $cid $groupName]
            if {$gname ni $listOfNoSkinGroups} {lappend listOfNoSkinGroups $gname}
        }
    }

    return $listOfNoSkinGroups
}

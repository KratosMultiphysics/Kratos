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
    set solverSettingsDict [dict merge $solverSettingsDict [write::getSolversParametersDict Fluid] ]
    # Parts
    dict set solverSettingsDict volume_model_part_name {*}[write::getPartsMeshId]
    # Skin parts
    dict set solverSettingsDict skin_parts [getBoundaryConditionMeshId]
    # No skin parts
    dict set solverSettingsDict no_skin_parts [getNoSkinConditionMeshId]

    dict set projectParametersDict solver_settings $solverSettingsDict

    # Boundary conditions processes
    dict set projectParametersDict initial_conditions_process_list [write::getConditionsParametersDict "FLNodalConditions" "Nodal"]
    dict set projectParametersDict boundary_conditions_process_list [write::getConditionsParametersDict $BCUN]
    dict set projectParametersDict gravity [list [getGravityProcessDict] ]

    return $projectParametersDict
}

proc Fluid::write::writeParametersEvent { } {
    set projectParametersDict [getParametersDict]
    write::SetParallelismConfiguration
    write::WriteJSON $projectParametersDict
}

# Gravity SubModelParts and Process collection
proc Fluid::write::getGravityProcessDict {} {
    set doc $gid_groups_conds::doc
    set root [$doc documentElement]

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

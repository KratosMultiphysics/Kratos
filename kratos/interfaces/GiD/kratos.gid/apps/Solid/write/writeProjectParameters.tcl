# Project Parameters

proc Solid::write::getParametersDict { } {
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
    dict set problemDataDict "parallel_type" $paralleltype
    if {$paralleltype eq "OpenMP"} {
        #set nthreads [write::getValue Parallelization OpenMPNumberOfThreads]
        #dict set problemDataDict NumberofThreads $nthreads
    } else {
        #set nthreads [write::getValue Parallelization MPINumberOfProcessors]
        #dict set problemDataDict NumberofProcessors $nthreads
    }
    set solutiontype [write::getValue SLSoluType]
    # Time Parameters
    if {$solutiontype eq "Static"} {
        dict set problemDataDict time_step "1.1"
        dict set problemDataDict start_time "0.0"
        dict set problemDataDict end_time "1.0"

    } elseif {$solutiontype eq "Dynamic"} {
        dict set problemDataDict time_step [write::getValue SLTimeParameters DeltaTime]
        dict set problemDataDict start_time [write::getValue SLTimeParameters StartTime]
        dict set problemDataDict end_time [write::getValue SLTimeParameters EndTime]
    }
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
    dict set modelDict input_file_label 0
    dict set solverSettingsDict model_import_settings $modelDict

    # Solution strategy parameters and Solvers
    set solverSettingsDict [dict merge $solverSettingsDict [write::getSolutionStrategyParametersDict] ]
    set solverSettingsDict [dict merge $solverSettingsDict [write::getSolversParametersDict Solid] ]

    dict set solverSettingsDict problem_domain_sub_model_part_list [write::getSubModelPartNames "SLParts"]
    dict set solverSettingsDict processes_sub_model_part_list [write::getSubModelPartNames "SLNodalConditions" "SLLoads"]

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

    return $projectParametersDict
}

proc Solid::write::writeParametersEvent { } {
    write::WriteJSON [getParametersDict]
    write::SetParallelismConfiguration
}

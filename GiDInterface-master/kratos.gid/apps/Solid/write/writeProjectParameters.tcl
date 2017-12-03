# Project Parameters

proc Solid::write::getParametersDict { } {
    set model_part_name "Solid_Domain"
    set projectParametersDict [dict create]

    # Problem data
    # Create section
    set problemDataDict [dict create]

    # Add items to section
    set model_name [file tail [GiD_Info Project ModelName]]
    dict set problemDataDict problem_name $model_name
    dict set problemDataDict model_part_name $model_part_name
    set nDim [expr [string range [write::getValue nDim] 0 0] ]
    dict set problemDataDict dimension $nDim

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
    set strategy_write_name [[::Model::GetSolutionStrategy $currentStrategyId] getAttribute "python_module"]
    dict set solverSettingsDict solver_type $strategy_write_name
    #~ dict set solverSettingsDict dimension [expr $nDim]
    dict set solverSettingsDict echo_level $echo_level
    dict set solverSettingsDict solution_type [write::getValue SLSoluType]

    if {$solutiontype eq "Static"} {
	dict set solverSettingsDict scheme_type [write::getValue SLScheme]
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

    # Dofs
    dict set solverSettingsDict dofs [list {*}[DofsInElements] ]

    # Solution strategy parameters and Solvers
    set solverSettingsDict [dict merge $solverSettingsDict [write::getSolutionStrategyParametersDict] ]
    set solverSettingsDict [dict merge $solverSettingsDict [write::getSolversParametersDict Solid] ]

    dict set solverSettingsDict problem_domain_sub_model_part_list [write::getSubModelPartNames "SLParts"]
    dict set solverSettingsDict processes_sub_model_part_list [write::getSubModelPartNames "SLNodalConditions" "SLLoads"]

    dict set projectParametersDict solver_settings $solverSettingsDict

    # Lists of processes
    set nodal_conditions_dict [Solid::write::getConditionsParametersDict SLNodalConditions "Nodal"] 
    dict set projectParametersDict constraints_process_list $nodal_conditions_dict

    dict set projectParametersDict loads_process_list [Solid::write::getConditionsParametersDict SLLoads]

    # GiD output configuration
    dict set projectParametersDict output_configuration [Solid::write::GetDefaultOutputDict]

    # restart options
    set restartDict [dict create ]
    dict set restartDict SaveRestart false
    dict set restartDict RestartFrequency 0
    dict set restartDict LoadRestart false
    dict set restartDict Restart_Step 0
    dict set projectParametersDict restart_options $restartDict
        
    return $projectParametersDict
}

proc Solid::write::DofsInElements { } {
    set dofs [list ]
    set root [customlib::GetBaseRoot]
    set xp1 "[spdAux::getRoute SLParts]/group/value\[@n='Element'\]"
    set elements [$root selectNodes $xp1]
    foreach element_node $elements {
        set elemid [$element_node @v]
        set elem [Model::getElement $elemid]
        foreach dof [split [$elem getAttribute "Dofs"] ","] {
            if {$dof ni $dofs} {lappend dofs $dof}
        }
    }
    return {*}$dofs
}

proc Solid::write::writeParametersEvent { } {
    write::WriteJSON [getParametersDict]
    write::SetParallelismConfiguration
}

# Project Parameters

proc Structural::write::getOldParametersDict { } {
    set model_part_name "Structure"
    set projectParametersDict [dict create]

    # Problem data
    # Create section
    set problemDataDict [dict create]

    # Add items to section
    set model_name [file tail [GiD_Info Project ModelName]]
    dict set problemDataDict problem_name $model_name
    dict set problemDataDict model_part_name $model_part_name
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
    set solutiontype [write::getValue STSoluType]
    # Time Parameters
    if {$solutiontype eq "Static"} {
        dict set problemDataDict time_step "1.1"
        dict set problemDataDict start_time "0.0"
        dict set problemDataDict end_time "1.0"

    } elseif {$solutiontype eq "Dynamic"} {
        dict set problemDataDict time_step [write::getValue STTimeParameters DeltaTime]
        dict set problemDataDict start_time [write::getValue STTimeParameters StartTime]
        dict set problemDataDict end_time [write::getValue STTimeParameters EndTime]
    }
    set echo_level [write::getValue Results EchoLevel]
    dict set problemDataDict echo_level $echo_level
    # Add section to document
    dict set projectParametersDict problem_data $problemDataDict

    # Solution strategy
    set solverSettingsDict [dict create]
    set currentStrategyId [write::getValue STSolStrat]
    set currentStrategyId [write::getValue STSolStrat]
    # set strategy_write_name [[::Model::GetSolutionStrategy $currentStrategyId] getAttribute "n"]
    dict set solverSettingsDict solver_type $solutiontype
    #~ dict set solverSettingsDict domain_size [expr $nDim]
    dict set solverSettingsDict echo_level $echo_level
    dict set solverSettingsDict analysis_type [write::getValue STAnalysisType]

    if {$solutiontype eq "Dynamic"} {
        dict set solverSettingsDict time_integration_method [write::getValue STSolStrat]
        dict set solverSettingsDict scheme_type [write::getValue STScheme]
    }

    # Model import settings
    set modelDict [dict create]
    dict set modelDict input_type "mdpa"
    dict set modelDict input_filename $model_name
    dict set solverSettingsDict model_import_settings $modelDict

    set materialsDict [dict create]
    dict set materialsDict materials_filename [GetAttribute materials_file]
    dict set solverSettingsDict material_import_settings $materialsDict

    # Solution strategy parameters and Solvers
    set solverSettingsDict [dict merge $solverSettingsDict [write::getSolutionStrategyParametersDict] ]
    set solverSettingsDict [dict merge $solverSettingsDict [write::getSolversParametersDict Structural] ]

    dict set solverSettingsDict problem_domain_sub_model_part_list [write::getSubModelPartNames [GetAttribute parts_un]]
    dict set solverSettingsDict processes_sub_model_part_list [write::getSubModelPartNames [GetAttribute nodal_conditions_un] [GetAttribute conditions_un] ]

    
    if {[usesContact]} {
        
        dict set solverSettingsDict contact_settings mortar_type "ALMContactFrictionless"

        set convergence_criterion [dict get $solverSettingsDict convergence_criterion]
        dict set solverSettingsDict convergence_criterion "contact_$convergence_criterion"
    }

    dict set projectParametersDict solver_settings $solverSettingsDict

    # Lists of processes
    set nodal_conditions_dict [write::getConditionsParametersDict [GetAttribute nodal_conditions_un] "Nodal"]
    lassign [ProcessContacts $nodal_conditions_dict] nodal_conditions_dict contact_conditions_dict
    dict set projectParametersDict constraints_process_list $nodal_conditions_dict
    dict set projectParametersDict contact_process_list $contact_conditions_dict

    dict set projectParametersDict loads_process_list [write::getConditionsParametersDict [GetAttribute conditions_un]]

    # GiD output configuration
    dict set projectParametersDict output_configuration [write::GetDefaultOutputDict]

    # Restart options
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

    set check_list [list "UpdatedLagrangianElementUP2D" "UpdatedLagrangianElementUPAxisym"]
    foreach elem $check_list {
        if {$elem in [Structural::write::GetUsedElements Name]} {
            dict set projectParametersDict pressure_dofs true
            break
        }
    }

    # set materialsDict [dict create]
    # dict set materialsDict materials_filename [GetAttribute materials_file]
    # dict set projectParametersDict material_import_settings $materialsDict

    return $projectParametersDict
}

proc Structural::write::ProcessContacts { nodal_conditions_dict } {
    set process_list [list ]
    set contact_process_list [list ]
    foreach elem $nodal_conditions_dict {
        if {[dict get $elem python_module] in {"alm_contact_process"}} {
            set model_part_name "Structure"
            dict set elem Parameters contact_model_part [dict get $elem Parameters model_part_name]
            dict set elem Parameters model_part_name $model_part_name
            dict set elem Parameters computing_model_part_name "computing_domain"
            lappend contact_process_list $elem
        } else {
            lappend process_list $elem
        }
    }
    return [list $process_list $contact_process_list]
}

proc Structural::write::writeParametersEvent { } {
    write::WriteJSON [getParametersDict]
    
}


# Project Parameters
proc Structural::write::getParametersEvent { } {
    set project_parameters_dict [getOldParametersDict]
    dict set project_parameters_dict solver_settings rotation_dofs [UsingRotationDofElements]
    set solverSettingsDict [dict get $project_parameters_dict solver_settings]
    set solverSettingsDict [dict merge $solverSettingsDict [write::getSolversParametersDict Structural] ]
    dict set project_parameters_dict solver_settings $solverSettingsDict
    return $project_parameters_dict
}
proc Structural::write::writeParametersEvent { } {
    write::WriteJSON [getParametersEvent]
}

proc Structural::write::UsingRotationDofElements { } {
    set root [customlib::GetBaseRoot]
    set xp1 "[spdAux::getRoute [GetAttribute parts_un]]/group/value\[@n='Element'\]"
    set elements [$root selectNodes $xp1]
    set bool false
    foreach element_node $elements {
        set elemid [$element_node @v]
        set elem [Model::getElement $elemid]
        if {[write::isBooleanTrue [$elem getAttribute "RotationDofs"]]} {set bool true; break}
    }

    return $bool
}

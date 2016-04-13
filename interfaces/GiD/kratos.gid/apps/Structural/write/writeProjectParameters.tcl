# Project Parameters
proc Structural::write::writeParametersEvent { } {
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
    set paralleltype [write::getValue STParallelType]
    if {$paralleltype eq "OpenMP"} {
        set nthreads [write::getValue STParallelization OpenMPNumberOfThreads]
        dict set problemDataDict NumberofThreads $nthreads
    } else {
        set nthreads [write::getValue STParallelization MPINumberOfProcessors]
        dict set problemDataDict NumberofProcessors $nthreads
    }
    
    # Time Parameters
    dict set problemDataDict time_step [write::getValue STTimeParameters DeltaTime]
    dict set problemDataDict end_time [write::getValue STTimeParameters EndTime]
    set echo_level [write::getValue STResults EchoLevel]
    dict set problemDataDict EchoLevel $echo_level
    
    # Add section to document
    dict set projectParametersDict problem_data $problemDataDict
    
    
    # Solution strategy
    set solverSettingsDict [dict create]
    set currentStrategyId [write::getValue STSolStrat]
    set strategy_write_name [[::Model::GetSolutionStrategy $currentStrategyId] getAttribute "ImplementedInPythonFile"]
    dict set solverSettingsDict solver_type $strategy_write_name
    dict set solverSettingsDict domain_size [expr $nDim]
    dict set solverSettingsDict echo_level $echo_level
    dict set solverSettingsDict solution_type [write::getValue STSoluType]
    
    # model import settings
    set modelDict [dict create]
    dict set modelDict input_type "mdpa"
    dict set modelDict input_filename $model_name
    dict set solverSettingsDict model_import_settings $modelDict
    
    # Solution strategy parameters and Solvers
    set solverSettingsDict [dict merge $solverSettingsDict [write::getSolutionStrategyParametersDict] ]
    set solverSettingsDict [dict merge $solverSettingsDict [write::getSolversParametersDict] ]
    dict set projectParametersDict solver_settings $solverSettingsDict
    
    # Lists of processes
    dict set projectParametersDict constraints_process_list [write::getConditionsParametersDict STNodalConditions "Nodal"]
    dict set projectParametersDict loads_process_list [write::getConditionsParametersDict STLoads]

    # GiD output configuration
    set outputDict [dict create ]
    
    set GiDPostDict [dict create]
    dict set GiDPostDict GiDPostMode                 [write::getValue STResults GiDPostMode]
    dict set GiDPostDict WriteMeshFlag               [write::getValue STResults GiDWriteMeshFlag]
    dict set GiDPostDict WriteConditionsFlag         [write::getValue STResults GiDWriteConditionsFlag]
    dict set GiDPostDict WriteParticlesFlag          [write::getValue STResults GiDWriteParticlesFlag]
    dict set GiDPostDict gid_write_frequency         [write::getValue STResults OutputDeltaTime]
    dict set outputDict gidpost_flags $GiDPostDict
    
    dict set outputDict write_results "PreMeshing"
    dict set outputDict plot_graphs false
    dict set outputDict plot_frequency 0
    dict set outputDict print_lists true
    dict set outputDict file_list [list ]
    dict set outputDict output_time 0.01
    dict set outputDict volume_output true
    dict set outputDict add_skin true
    
    dict set outputDict nodal_results [write::GetResultsList "STNodalResults"]
    dict set outputDict gauss_points_results [write::GetResultsList "STElementResults"]
    
    dict set projectParametersDict output_configuration $outputDict
    
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


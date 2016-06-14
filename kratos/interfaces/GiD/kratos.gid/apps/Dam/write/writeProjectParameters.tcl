### Project Parameters
proc Dam::write::getParametersDict { } {
    
    set projectParametersDict [dict create]
    
    ### Problem data
    ### Create section
    set generalDataDict [dict create]
    
    ### Add items to section
    set model_name [file tail [GiD_Info Project ModelName]]
    dict set generalDataDict problem_name $model_name
    dict set generalDataDict model_part_name "MainModelPart"
    set nDim [expr [string range [write::getValue nDim] 0 0] ]
    dict set generalDataDict domain_size $nDim
    dict set generalDataDict NumberofThreads [write::getValue DamNumThreads ]
    dict set generalDataDict time_scale [write::getValue DamTimeParameters TimeScale]
    dict set generalDataDict evolution_type [write::getValue DamEvolutionType] 
    dict set generalDataDict delta_time [write::getValue DamTimeParameters DeltaTime]
    dict set generalDataDict ending_time [write::getValue DamTimeParameters EndingTime]
    
    ### Add section to document
    dict set projectParametersDict general_data $generalDataDict
    
    ### Solver Data
    ### Diffusion settings
    set diffusionSolverSettingsDict [dict create]
    dict set diffusionSolverSettingsDict unknown_variable "TEMPERATURE"
    dict set diffusionSolverSettingsDict difussion_variable "CONDUCTIVITY"
    dict set diffusionSolverSettingsDict specific_heat_variable "SPECIFIC_HEAT"
    dict set diffusionSolverSettingsDict density_variable "DENSITY"
    set damTypeofProblem [write::getValue DamTypeofProblem]
    if {$damTypeofProblem eq "Thermo-Mechanical"} {
        dict set diffusionSolverSettingsDict temporal_scheme [write::getValue DamMechanicalSchemeTherm]
        dict set diffusionSolverSettingsDict reference_temperature [write::getValue DamReferenceTemperature]
    } 
    
    ### Add section to document
    dict set projectParametersDict diffusion_settings $diffusionSolverSettingsDict
    
    ### Mechanical Settings
    set mechanicalSolverSettingsDict [dict create]
    dict set mechanicalSolverSettingsDict solution_type [write::getValue DamSoluType]
    dict set mechanicalSolverSettingsDict analysis_type [write::getValue DamAnalysisType]
    dict set mechanicalSolverSettingsDict strategy_type "Newton-Raphson"
    #set mechanicalSolverSettingsDict [dict merge $mechanicalSolverSettingsDict [write::getSolutionStrategyParametersDict] ]
    #dict set mechanicalSolverSettingsDict max_iteration [write::getValue DamSolStrat DamMaxIter]
    #dict set mechanicalSolverSettingsDict dofs_relative_tolerance [write::getValue DofsTol]
    #dict set mechanicalSolverSettingsDict residual_relative_tolerance  [write::getValue RelTol]
    set damTypeofSolver [write::getValue DamTypeofsolver]
    if {$damTypeofSolver eq "Direct"} {
        dict set mechanicalSolverSettingsDict direct_solver [write::getValue DamDirectsolver]
    } elseif {$damTypeofSolver eq "Iterative"} {
        dict set mechanicalSolverSettingsDict direct_solver [write::getValue DamIterativesolver]
    }
    ### Add section to document
    dict set projectParametersDict mechanical_settings $mechanicalSolverSettingsDict
    
    ### Submodels parts
    #set submodelpartDict [dict create]
    #dict set submodelpartDict problem_domain_sub_model_part_list [getSubModelPartNames "DamParts"]
    #dict set submodelpartDict nodal_processes_sub_model_part_list [getSubModelPartNames "DamNodalConditions" ]
    #dict set submodelpartDict load_processes_sub_model_part_list [getSubModelPartNames "DamLoads" ]
    
    ### GiD output configuration
    dict set projectParametersDict output_configuration [write::GetDefaultOutputDict]
        
    return $projectParametersDict
}

proc Dam::write::writeParametersEvent { } {
    set projectParametersDict [getParametersDict]
    write::WriteJSON $projectParametersDict
}

proc write::GetDefaultOutputDict {} {
    set outputDict [dict create]
    set resultDict [dict create]
    
    set GiDPostDict [dict create]
    dict set GiDPostDict GiDPostMode                [getValue Results GiDPostMode]
    dict set GiDPostDict WriteDeformedMeshFlag      [getValue Results GiDWriteMeshFlag]
    dict set GiDPostDict WriteConditionsFlag        [getValue Results GiDWriteConditionsFlag]
    dict set GiDPostDict MultiFileFlag              [getValue Results GiDMultiFileFlag]
    dict set resultDict gidpost_flags $GiDPostDict
    
    dict set resultDict output_frequency [write::getValue Results OutputDeltaTime]   
    dict set resultDict nodal_results [GetResultsList NodalResults]
    dict set resultDict gauss_point_results [GetResultsList ElementResults]
    
    dict set outputDict "result_file_configuration" $resultDict
    return $outputDict
}

#proc write::getSolutionStrategyParametersDict {} {
    #set solStratUN [apps::getCurrentUniqueName SolStrat]
    #set schemeUN [apps::getCurrentUniqueName Scheme]
    
    #set solstratName [write::getValue $solStratUN]
    #set schemeName [write::getValue $schemeUN]
    #set sol [::Model::GetSolutionStrategy $solstratName]
    #set sch [$sol getScheme $schemeName]
    
    #set paramsPath [apps::getCurrentUniqueName StratParams]
    
    #foreach {n in} [$sol getInputs] {
	#dict set mechanicalSolverSettingsDict $n [write::getValue $paramsPath $n ]
    #}
    #foreach {n in} [$sch getInputs] {
	#dict set mechanicalSolverSettingsDict $n [write::getValue $paramsPath $n ]
    #}
    #return $mechanicalSolverSettingsDict
#}

#proc Dam::write::getSubModelPartNames { args } {
    #set doc $gid_groups_conds::doc
    #set root [$doc documentElement]
    
    #set listOfProcessedGroups [list ]
    #set groups [list ]
    #foreach un $args {
        #set xp1 "[spdAux::getRoute $un]/condition/group"
        #set xp2 "[spdAux::getRoute $un]/group"
        #set grs [$root selectNodes $xp1]
        #if {$grs ne ""} {lappend groups {*}$grs}
        #set grs [$root selectNodes $xp2]
        #if {$grs ne ""} {lappend groups {*}$grs}
    #}
    #foreach group $groups {
        #set groupName [$group @n]
        #set cid [[$group parent] @n]
        #set gname [::write::getMeshId $cid $groupName]
        #if {$gname ni $listOfProcessedGroups} {lappend listOfProcessedGroups $gname}
    #}
    
    #return $listOfProcessedGroups
#}


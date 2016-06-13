# Project Parameters
proc ::Dam::write::getParametersDict { } {
    
    set projectParametersDict [dict create]
    
    # Problem data
    # Create section
    set generalDataDict [dict create]
    
    # Add items to section
    set model_name [file tail [GiD_Info Project ModelName]]
    dict set generalDataDict problem_name $model_name
    dict set generalDataDict model_part_name "MainModelPart"
    set nDim [expr [string range [write::getValue nDim] 0 0] ]
    dict set generalDataDict domain_size $nDim
    dict set generalDataDict time_scale [write::getValue DamTimeParameters TimeScale]
    dict set generalDataDict evolution_type [write::getValue DamEvolutionType] 
    dict set generalDataDict delta_time [write::getValue DamTimeParameters DeltaTime]
    dict set generalDataDict ending_time [write::getValue DamTimeParameters EndingTime]
    
    ## Add section to document
    dict set projectParametersDict general_data $generalDataDict
    
    # Solver Data
    # Diffusion settings
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
    
    ## Add section to document
    dict set projectParametersDict diffusion_settings $diffusionSolverSettingsDict
    
    # Mechanical Settings
    set mechanicalSolverSettingsDict [dict create]
    dict set mechanicalSolverSettingsDict solution_type [write::getValue DamSoluType]
    dict set mechanicalSolverSettingsDict analysis_type [write::getValue DamAnalysisType]
    dict set mechanicalSolverSettingsDict strategy_type "Newton-Raphson"
    dict set mechanicalSolverSettingsDict max_iteration [write::getValue MaxIter]
    dict set mechanicalSolverSettingsDict dofs_relative_tolerance [write::getValue DofsTol]
    dict set mechanicalSolverSettingsDict residual_relative_tolerance  [write::getValue RelTol]
    set damTypeofSolver [write::getValue DamTypeofsolver]
    if {$damTypeofSolver eq "Direct"} {
        dict set mechanicalSolverSettingsDict direct_solver [write::getValue DamDirectsolver]
    } elseif {$damTypeofSolver eq "Iterative"} {
        dict set mechanicalSolverSettingsDict direct_solver [write::getValue DamIterativesolver]
    }
    ## Add section to document
    dict set projectParametersDict mechanical_settings $mechanicalSolverSettingsDict
    
    ## model import settings
    set modelDict [dict create]
    dict set modelDict input_type "mdpa"
    dict set modelDict input_filename $model_name
    dict set projectParametersDict model_import_settings $modelDict

    
    return $projectParametersDict
}

proc Dam::write::writeParametersEvent { } {
    set projectParametersDict [getParametersDict]
    write::WriteJSON $projectParametersDict
}

proc Dam::write::getSubModelPartNames { args } {
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

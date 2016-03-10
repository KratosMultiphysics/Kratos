# Project Parameters
proc Fluid::write::writeParametersEvent { } {
    set projectParametersDict [dict create]
    
    # First section -> Problem data
    set problemDataDict [dict create]
    set model_name [file tail [GiD_Info Project ModelName]]
    dict set problemDataDict ProblemName $model_name
    dict set problemDataDict ModelPartName "MainModelPart"
    set nDim [expr [string range [write::getValue nDim] 0 0] ]
    dict set problemDataDict DomainSize $nDim
    
   
    # Time Parameters
    dict set problemDataDict start_step [expr [write::getValue FLTimeParameters StartTime] ]
    dict set problemDataDict end_time [expr [write::getValue FLTimeParameters EndTime]]
    dict set problemDataDict time_step [expr [write::getValue FLTimeParameters DeltaTime]]
    dict set problemDataDict divergence_step [expr [write::getValue FLTimeParameters DivergenceCleareanceStep]]
    
    dict set projectParametersDict problem_data $problemDataDict
    
    # output configuration
    set outputConfigDict [dict create]
    dict set outputConfigDict output_filename "[file tail [GiD_Info Project ModelName]].out"
    dict set outputConfigDict GiDPostMode [write::getValue FLResults GiDPostMode]
    dict set outputConfigDict GiDMultiFileFlag [write::getValue FLResults GiDMultiFileFlag]
    dict set outputConfigDict GiDWriteMeshFlag True
    dict set outputConfigDict GiDWriteConditionsFlag True
    dict set outputConfigDict GiDWriteParticlesFlag False
    dict set outputConfigDict GiDWriteFrequency [expr [write::getValue FLResults OutputDeltaTime]]
    dict set outputConfigDict PlotGraphs False
    dict set outputConfigDict PlotFrequency 0
    dict set outputConfigDict PrintLists True
    dict set outputConfigDict output_time [expr [write::getValue FLResults OutputDeltaTime]]
    dict set outputConfigDict VolumeOutput [expr [write::getValue FLResults VolumeOutput]]
    
    # on nodes
    dict set outputConfigDict nodal_results [list "VELOCITY" "PRESSURE"]
    
    
    #set xp1 "[apps::getRoute FLResults]/container\[@n='OnNodes'\]/value"
    # on elements
    dict set outputConfigDict gauss_points_results [list "VORTICITY"]
    
    dict set projectParametersDict output_configuration $outputConfigDict
    
    # restart options
    set restartDict [dict create]
    dict set restartDict SaveRestart False
    dict set restartDict RestartFrequency 0
    dict set restartDict LoadRestart False
    dict set restartDict Restart_Step 0
    dict set projectParametersDict restart_options $restartDict
    
    
    # Solver settings
    set solverSettingsDict [dict create]
    dict set solverSettingsDict solver_type navier_stokes_solver_fractionalstep
    dict set solverSettingsDict DomainSize $nDim
    dict set solverSettingsDict echo_level 1
    
    # model import settings
    set modelDict [dict create]
    dict set modelDict input_type "mdpa"
    dict set modelDict input_filename $model_name
    dict set solverSettingsDict model_import_settings $modelDict
    
    
    set solverSettingsDict [getSolutionStrategyParameters $solverSettingsDict]
    set solverSettingsDict [getSolversParameters $solverSettingsDict]
    # Parts
    dict set solverSettingsDict volume_model_part_name {*}[write::getPartsMeshId]
    # Skin parts
    dict set solverSettingsDict skin_parts [getBoundaryConditionMeshId]
    
    dict set solverSettingsDict velocity_tolerance 1e-3
    dict set solverSettingsDict pressure_tolerance 1e-2
    dict set solverSettingsDict maximum_velocity_iterations [expr 3]
    dict set solverSettingsDict maximum_pressure_iterations [expr 3]
    dict set solverSettingsDict predictor_corrector false
    dict set solverSettingsDict echo_level 1
    dict set solverSettingsDict DomainSize $nDim 
    dict set solverSettingsDict consider_periodic_conditions false
    dict set solverSettingsDict time_order 2
    dict set solverSettingsDict dynamic_tau 0.001
    dict set solverSettingsDict compute_reactions false
    dict set solverSettingsDict divergence_clearance_steps 0
    dict set solverSettingsDict reform_dofs_at_each_iteration false
    
    dict set projectParametersDict solver_settings $solverSettingsDict
        
    # Boundary conditions processes
    dict set projectParametersDict boundary_conditions_process_list [getBoundaryConditionsParameters]
    
    write::WriteProcess $projectParametersDict
}

proc Fluid::write::getSolutionStrategyParameters {solverSettingsDict} {
    set solstratName [write::getValue FLSolStrat]
    set schemeName [write::getValue FLScheme]
    set sol [::Model::GetSolutionStrategy $solstratName]
    set sch [$sol getScheme $schemeName]
    
    foreach {n in} [$sol getInputs] {
        dict set solverSettingsDict $n [expr [write::getValue FLStratParams $n ]]
    }
    foreach {n in} [$sch getInputs] {
        dict set solverSettingsDict $n [expr [write::getValue FLStratParams $n ]]
    }
    return $solverSettingsDict
}

proc Fluid::write::getSolversParameters {solverSettingsDict} {
    set solstratName [write::getValue FLSolStrat]
    set sol [::Model::GetSolutionStrategy $solstratName]
    foreach se [$sol getSolversEntries] {
        set solverEntryDict [dict create]
        set un "FL$solstratName[$se getName]"
        set solverName [write::getValue $un Solver]
        dict set solverEntryDict solver_type $solverName
          
        foreach {n in} [[::Model::GetSolver $solverName] getInputs] {
            dict set solverEntryDict $n [expr [write::getValue $un $n]]
        }
        dict set solverSettingsDict [$se getName] $solverEntryDict
        unset solverEntryDict
    }
    return $solverSettingsDict
}

proc Fluid::write::printInitialConditions {spacing} {
    set doc $gid_groups_conds::doc
    set root [$doc documentElement]
    set s [write::getSpacing $spacing]

    write::WriteString "#Constraints Data"
    write::WriteString "#################################################"
    write::WriteString "${s}constraints_process_list = \["
    set xp1 "[apps::getRoute "FLNodalConditions"]/condition/group"
    set groups [$root selectNodes $xp1]
    set str ""
    foreach group $groups {
        set groupName [$group @n]
        set cid [[$group parent] @n]
        set groupId [::write::getMeshId $cid $groupName]
        #W [[$group parent] @type]
        if {[[$group parent] @type] eq "vector"} {
            set FixX False
            set FixY False
            set FixZ False
            if {[$group find n FixX] ne ""} {
                set FixX [expr [get_domnode_attribute [$group find n FixX] v] ? True : False]
            }
            if {[$group find n FixY] ne ""} {
                set FixY [expr [get_domnode_attribute [$group find n FixY] v] ? True : False]
            }
            if {[$group find n FixZ] ne ""} {
                set FixZ [expr [get_domnode_attribute [$group find n FixZ] v] ? True : False]
            }
            
            set ValX [get_domnode_attribute [$group find n ValX] v] 
            set ValY [get_domnode_attribute [$group find n ValY] v] 
            set ValZ [get_domnode_attribute [$group find n ValZ] v]
            set factornode [$group find n FACTOR]
            set Factor 1
            if {$factornode ne ""} {
                set Factor [get_domnode_attribute $factornode v]
            }
            set arguments [list "a" "b"]
            #write::WriteProcess "ApplyConstantVectorValueProcess" $arguments
            append str "{ \"process_name\" : \"ApplyConstantVectorValueProcess\",
              \"implemented_in_module\" : \"KratosMultiphysics\",
              \"implemented_in_file\": \"process_factory\",
              
              \"parameters\" : { 
                  \"mesh_id\": $groupId,
                  \"model_part_name\" : ModelPartName,
                  \"variable_name\": \"[[$group parent] @n]\",
                  \"factor\": $Factor,
                  \"value\": \[$ValX,$ValY,$ValZ\],
                  \"is_fixed_x\" : $FixX,
                  \"is_fixed_y\" : $FixY,
                  \"is_fixed_z\" : $FixZ,
                  }
            }, "
        } {
            set Value [get_domnode_attribute [$group firstChild] v] 
            set Fixed True
            
            append str "{ \"process_name\" : \"ApplyConstantScalarValueProcess\",
            \"implemented_in_python\" : True,
            \"implemented_in_module\" : \"KratosMultiphysics\",
            \"implemented_in_file\": \"process_factory\",
                        
            \"parameters\" : { 
                \"mesh_id\": $groupId,
                \"model_part_name\" : ModelPartName,
                \"variable_name\": \"[[$group parent] @n]\",
                \"value\": $Value,
                \"is_fixed\": $Fixed
                }
          }, "
        }
    }
    # Quitar la coma
    if {[llength $groups]} {set str [string range $str 0 end-1]}
    write::WriteString $str
    write::WriteString "]"
}


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


proc Fluid::write::getBoundaryConditionsParameters {} {
    variable BCUN
    set doc $gid_groups_conds::doc
    set root [$doc documentElement]

    set bcCondsDict [list ]
    
    set xp1 "[apps::getRoute $BCUN]/condition/group"
    set groups [$root selectNodes $xp1]    
    foreach group $groups {
        set groupName [$group @n]
        set cid [[$group parent] @n]
        set groupId [::write::getMeshId $cid $groupName]
        set condId [[$group parent] @n]
        set condition [::Model::getCondition $condId]
        set processName [$condition getProcessName]
        set process [::Model::GetProcess $processName]
        set processDict [dict create]
        set paramDict [dict create]
        dict set paramDict mesh_id 0
        dict set paramDict model_part_name $groupId
        
        set process_attributes [$process getAttributes]
        set process_parameters [$process getInputs]
        
        dict set process_attributes process_name [dict get $process_attributes n]
        dict unset process_attributes n
        dict unset process_attributes pn
        
        set processDict [dict merge $processDict $process_attributes]
        foreach {inputName in_obj} $process_parameters {
            set in_type [$in_obj getType]
            if {$in_type eq "vector"} {
                set is_fixed_x [expr False]
                set is_fixed_y [expr False]
                set is_fixed_z [expr False]
                if {[$group find n FixX] ne ""} {
                    set is_fixed_x [expr [get_domnode_attribute [$group find n FixX] v] ? True : False]
                }
                if {[$group find n FixY] ne ""} {
                    set is_fixed_y [expr [get_domnode_attribute [$group find n FixY] v] ? True : False]
                }
                if {[$group find n FixZ] ne ""} {
                    set is_fixed_z [expr [get_domnode_attribute [$group find n FixZ] v] ? True : False]
                }
                set ValX [expr [get_domnode_attribute [$group find n ${inputName}X] v] ]
                set ValY [expr [get_domnode_attribute [$group find n ${inputName}Y] v] ] 
                set ValZ [expr 0.0]
                catch {set ValZ [expr [get_domnode_attribute [$group find n ${inputName}Z] v]]}
                
                dict set paramDict is_fixed_x $is_fixed_x
                dict set paramDict is_fixed_y $is_fixed_y
                dict set paramDict is_fixed_z $is_fixed_z
                dict set paramDict $inputName [list $ValX $ValY $ValZ]
            } elseif {$in_type eq "double"} {
                set value [get_domnode_attribute [$group find n $inputName] v] 
                if {[$group find n Fix] ne ""} {
                    set is_fixed [expr [get_domnode_attribute [$group find n Fix] v] ? True : False]
                    dict set paramDict is_fixed $is_fixed
                }
                dict set paramDict $inputName [expr $value]
            }
        }
        
        dict set processDict Parameters $paramDict
        lappend bcCondsDict $processDict
    }
    return $bcCondsDict
}

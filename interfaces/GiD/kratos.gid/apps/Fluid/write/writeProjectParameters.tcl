# Project Parameters
proc Fluid::write::writeParametersEvent { } {
    write::WriteString "ProblemName =\"[file tail [GiD_Info Project ModelName]]\""
    write::WriteString ""
    write::WriteString "#Problem Data"
    write::WriteString "#################################################"
    write::WriteString "ModelPartName = \"Structure\""
    write::WriteString "DomainSize = [string range [write::getValue nDim] 0 0]"
    
     # Parallelization
    set paralleltype [write::getValue FLParallelization ParallelSolutionType]
    if {$paralleltype eq "OpenMP"} {
        set nthreads [write::getValue FLParallelization OpenMPNumberOfThreads]
        write::WriteString "NumberofThreads = $nthreads"
    } else {
        set nthreads [write::getValue FLParallelization MPINumberOfProcessors]
        write::WriteString "NumberofProcessors = $nthreads"
    }
    
    # Time Parameters
    write::WriteString "start_step = [write::getValue FLTimeParameters StartTime]"
    write::WriteString "end_time = [write::getValue FLTimeParameters EndTime]"
    write::WriteString "time_step = [write::getValue FLTimeParameters DeltaTime]"
    write::WriteString "divergence_step = [write::getValue FLTimeParameters DivergenceCleareanceStep]"
    
    #write::WriteString "EchoLevel = [write::getValue FLResults EchoLevel]"
    
    # Solution strategy
    write::WriteString "#Solver Data"
    write::WriteString "#################################################"
    write::WriteString "class SolverSettings:"
    write::WriteString "    solver_type = \"fluid_python_solver\""
    write::WriteString "    domain_size = DomainSize"
    write::WriteString "    echo_level = EchoLevel"
    #write::WriteString "    solution_type  = \"[write::getValue FLSoluType]\""
    write::WriteString "    "
    
    printSolutionStrategy 4
    printSolvers 4
    
    printInitialConditions 0
    printConds 0
    
    # GiD output configuration
    write::WriteString "class GidOutputConfiguration:"
    write::WriteString "    GiDPostMode = \"[write::getValue FLResults GiDPostMode]\""
    write::WriteString "    GiDMultiFileFlag = \"[write::getValue FLResults GiDMultiFileFlag]\""

    write::WriteString ""   
    write::WriteString "GiDWriteFrequency = [write::getValue FLResults OutputDeltaTime]"
    write::WriteString ""
    write::WriteString "# graph_options"
    write::WriteString "PlotGraphs = \"False\""
    write::WriteString "PlotFrequency = 0"
    write::WriteString ""
    write::WriteString "# list options"
    write::WriteString "PrintLists = \"True\""
    write::WriteString "file_list = \[\] "
    write::WriteString ""
    write::WriteString "# restart options"
    write::WriteString "SaveRestart = False"
    write::WriteString "RestartFrequency = 0"
    write::WriteString "LoadRestart = False"
    write::WriteString "Restart_Step = 0"
    write::WriteString ""
}

proc Fluid::write::printSolutionStrategy {spacing} {
    set solstratName [write::getValue FLSolStrat]
    set schemeName [write::getValue FLScheme]
    set sol [::Model::GetSolutionStrategy $solstratName]
    set sch [$sol getScheme $schemeName]
    
    set spaces [write::getSpacing $spacing]
    write::WriteString "${spaces}scheme_type = \"[$sch getName]\""
    write::WriteString "${spaces}RotationDofs = False"
    write::WriteString "${spaces}PressureDofs = False"
    
    write::WriteString "${spaces}time_integration_method = \"[$sol getName]\""
    foreach {n in} [$sol getInputs] {
        write::WriteString "$spaces$n = [write::getValue FLStratParams $n ]"
    }
    foreach {n in} [$sch getInputs] {
        write::WriteString "$spaces$n = [write::getValue FLStratParams $n ]"
    }
}

proc Fluid::write::printSolvers {spacing} {
    set solstratName [write::getValue FLSolStrat]
    set sol [::Model::GetSolutionStrategy $solstratName]
    set spaces [write::getSpacing $spacing]
    set spaces2 [write::getSpacing [expr $spacing +4]]
    foreach se [$sol getSolversEntries] {
        set un "FL$solstratName[$se getName]"
        set solverName [write::getValue $un Solver]
        write::WriteString "${spaces}class [$se getName]:"
        write::WriteString "${spaces2}solver_type = \"$solverName\""
          
        foreach {n in} [[::Model::GetSolver $solverName] getInputs] {
            write::WriteString "$spaces2$n = [write::getValue $un $n ]"
        }
    }
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


proc Fluid::write::printConds {spacing} {
    set doc $gid_groups_conds::doc
    set root [$doc documentElement]
    set s [write::getSpacing $spacing]

    write::WriteString "#Loads Data"
    write::WriteString "#################################################"
    write::WriteString "${s}loads_process_list = \["
    set xp1 "[apps::getRoute "FLBC"]/condition/group"
    set groups [$root selectNodes $xp1]
    set str ""
    foreach group $groups {
        set groupName [$group @n]
        set cid [[$group parent] @n]
        #W "$cid $groupName"
        set groupId [::write::getMeshId $cid $groupName]
        set condId [[$group parent] @n]
        set condition [::Model::getCondition $condId]
        set type [$condition getProcessFormat]
        set processName [$condition getProcessName]
        #W "type $type pron $processName"
        
        if {$type eq "vector"} {
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
            
            append str "{ \"process_name\" : \"ApplyConstantVectorValueProcess\",
              \"implemented_in_module\" : \"KratosMultiphysics\",
              \"implemented_in_file\": \"process_factory\",
              
              \"parameters\" : { 
                \"mesh_id\": $groupId,
                  \"model_part_name\" : ModelPartName,
                  \"variable_name\": \"$processName\",
                  \"factor\": $Factor,
                  \"value\": \[$ValX,$ValY,$ValZ\],
                  \"is_fixed_x\" : $FixX,
                  \"is_fixed_y\" : $FixY,
                  \"is_fixed_z\" : $FixZ,
                  }
            }, "
        } elseif {$type eq "double"} {
            set Value [get_domnode_attribute [$group firstChild] v] 
            set Fixed True
            
            append str "{ \"process_name\" : \"ApplyConstantScalarValueProcess\",
            \"implemented_in_python\" : True,
            \"implemented_in_module\" : \"KratosMultiphysics\",
            \"implemented_in_file\": \"process_factory\",
                        
            \"parameters\" : { 
				\"mesh_id\": $groupId,
                \"model_part_name\" : ModelPartName,
                \"variable_name\": \"$processName\",
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




proc Fluid::write::printResults {spacing} {
    set doc $gid_groups_conds::doc
    set root [$doc documentElement]
    set s [write::getSpacing $spacing]

    write::WriteString "#PostProcess Data"
    write::WriteString "#################################################"
    
    set xp1 "[apps::getRoute "SMNodalResults"]/value"
    set results [$root selectNodes $xp1]
    set min 0
    set str "nodal_results=\["
    foreach res $results {
        if {[get_domnode_attribute $res v] eq "Yes"} {
            set min 1
            set name [get_domnode_attribute $res n]
            append str "\"$name\","
        }
    }
    
    if {$min} {set str [string range $str 0 end-1]}
    append str "]"
    write::WriteString $str
    
    set xp1 "[apps::getRoute "SMElementResults"]/value"
    set results [$root selectNodes $xp1]
    set min 0
    set str "gauss_points_results=\["
    foreach res $results {
        if {[get_domnode_attribute $res v] in [list "Yes" "True"] && [get_domnode_attribute $res state] eq "normal"} {
            set min 1
            set name [get_domnode_attribute $res n]
            append str "\"$name\","
        }
    }
    
    if {$min} {set str [string range $str 0 end-1]}
    append str "]"
    write::WriteString $str
    write::WriteString ""
}


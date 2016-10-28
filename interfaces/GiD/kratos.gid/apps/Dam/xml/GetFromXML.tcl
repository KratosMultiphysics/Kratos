namespace eval Dam::xml {
     variable dir
}

proc Dam::xml::Init { } {
     variable dir
     Model::InitVariables dir $Dam::dir

     Model::getSolutionStrategies Strategies.xml
     Model::getElements Elements.xml
     Model::getNodalConditions NodalConditions.xml
     Model::getConstitutiveLaws ConstitutiveLaws.xml
     Model::getProcesses Processes.xml
     Model::getConditions Conditions.xml
     Model::getSolvers "../../Common/xml/Solvers.xml"
}

proc Dam::xml::getUniqueName {name} {
    return Dam$name
}


proc ::Dam::xml::MultiAppEvent {args} {
   if {$args eq "init"} {
     spdAux::parseRoutes
     spdAux::ConvertAllUniqueNames SL Dam
   }
}

proc Dam::xml::ProcGetSchemes {domNode args} {
    set type_of_problem [lindex $args 0]
    
    set sol_stratUN "DamSolStratTherm"
    if {$type_of_problem eq "Mechanic"} {
          set sol_stratUN "DamSolStratMech"
    }
    
    set solStratName [write::getValue $sol_stratUN]
    set schemes [::Model::GetAvailableSchemes $solStratName]
    set ids [list ]
    if {[llength $schemes] == 0} {
        if {[get_domnode_attribute $domNode v] eq ""} {$domNode setAttribute v "None"}
        return "None"
    }
    set names [list ]
    set pnames [list ]
    foreach cl $schemes {
        lappend names [$cl getName]
        lappend pnames [$cl getName] 
        lappend pnames [$cl getPublicName]
    }
    
    $domNode setAttribute values [join $names ","]
    if {[get_domnode_attribute $domNode v] eq ""} {$domNode setAttribute v [lindex $names 0]}
    #spdAux::RequestRefresh
    
    return [join $pnames ","]
}


proc Dam::xml::SolStratParamState {outnode} {
    set doc $gid_groups_conds::doc
    set root [$doc documentElement]
    
    set solstrat_un "DamSolStratMech"
    
    #W $solstrat_un
    if {[get_domnode_attribute [$root selectNodes [spdAux::getRoute $solstrat_un]] v] eq ""} {
        get_domnode_attribute [$root selectNodes [spdAux::getRoute $solstrat_un]] dict
    }
    set SolStrat [get_domnode_attribute [$root selectNodes [spdAux::getRoute $solstrat_un]] v]
    set paramName [$outnode @n]
    set ret [::Model::CheckSolStratInputState $SolStrat $paramName]
    if {$ret} {
        lassign [Model::GetSolStratParamDep $SolStrat $paramName] depN depV
        foreach node [[$outnode parent] childNodes] {
            if {[$node @n] eq $depN} {
                if {[get_domnode_attribute $node v] ni [split $depV ,]} {
                    set ret 0
                    break
                }
            }
        }
    }
    return $ret
}

proc Dam::xml::ProcGetConstitutiveLaws {domNode args} {
     set Elementname [$domNode selectNodes {string(../value[@n='Element']/@v)}]
     set Claws [::Model::GetAvailableConstitutiveLaws $Elementname]
     
     set type_of_problem [write::getValue DamTypeofProblem]
     set goodList [list ]
     foreach cl $Claws {
          set type [$cl getAttribute Type]
          if {[string first "Therm" $type] eq -1 && $type_of_problem ne "Thermo-Mechanical"} {
               lappend goodList $cl
          } elseif {[string first "Therm" $type] ne -1 && $type_of_problem eq "Thermo-Mechanical"} {
               lappend goodList $cl
          } elseif {[string first "Interface" $type] ne -1} {lappend goodList $cl}
     }
     set Claws $goodList
     set analysis_type [write::getValue DamAnalysisType]
     set goodList [list ]
     foreach cl $Claws {
          set type [$cl getAttribute AnalysisType]
          if {$analysis_type eq "Non-Linear"} {
               lappend goodList $cl
          }
          if {$type ne "Non-Linear" && $analysis_type eq "Linear"} {
               lappend goodList $cl
          }
     }
     set Claws $goodList
   
     #W "Const Laws que han pasado la criba: $Claws"
     if {[llength $Claws] == 0} {
         if {[get_domnode_attribute $domNode v] eq ""} {$domNode setAttribute v "None"}
         return "None"
     }
     set names [list ]
     foreach cl $Claws {
         lappend names [$cl getName]
     }
     set values [join $names ","]
     if {[get_domnode_attribute $domNode v] eq "" || [get_domnode_attribute $domNode v] ni $names} {$domNode setAttribute v [lindex $names 0]; spdAux::RequestRefresh}
     
     
     return $values
}

Dam::xml::Init

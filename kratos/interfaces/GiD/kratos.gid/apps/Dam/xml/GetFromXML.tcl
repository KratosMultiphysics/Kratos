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

proc Dam::xml::ProcGetSchemes {domNode args} {
    set type_of_problem [lindex $args 0]
    
    set sol_stratUN "DamSolStratTherm"
    if {$type_of_problem eq "Mechanic"} {
          set sol_stratUN "DamSolStratMech"
    }
    
    set solStratName [get_domnode_attribute [$domNode selectNodes [spdAux::getRoute $sol_stratUN]] v]
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
        get_domnode_attribute [$root selectNodes [spdAux::getRoute $solstrat_un]] values
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

Dam::xml::Init

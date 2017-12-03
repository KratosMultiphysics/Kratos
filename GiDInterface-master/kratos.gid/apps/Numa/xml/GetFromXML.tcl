namespace eval Numa::xml {
    variable dir
}

proc Numa::xml::Init { } {
    variable dir
    Model::InitVariables dir $Numa::dir
    
    Model::getSolutionStrategies Strategies.xml
    Model::getElements Elements.xml
    Model::getMaterials Materials.xml
    Model::getNodalConditions NodalConditions.xml
    Model::getConstitutiveLaws ConstitutiveLaws.xml
    Model::getProcesses Processes.xml
    Model::getConditions Conditions.xml
    
}

proc Numa::xml::getUniqueName {name} {
    return Numa$name
}

proc ::Numa::xml::MultiAppEvent {args} {
    if {$args eq "init"} {
        spdAux::parseRoutes
        spdAux::ConvertAllUniqueNames SL Numa
    }
}

proc Numa::xml::ProcCheckNodalConditionState {domNode args} {
    set parts_un NumaParts
    if {[spdAux::getRoute $parts_un] ne ""} {
        set conditionId [$domNode @n]
        set elems [$domNode selectNodes "[spdAux::getRoute $parts_un]/group/value\[@n='Element'\]"]
        set elemnames [list ]
        foreach elem $elems { set a [get_domnode_attribute $elem values]; set a [get_domnode_attribute $elem dict]; lappend elemnames [$elem @v]}
        set elemnames [lsort -unique $elemnames]
        
        # Mirar Type of problem y acceder al contenedor correcto
        set TypeofProblem [get_domnode_attribute [$domNode selectNodes [spdAux::getRoute NumaTypeofProblem]] v]
        switch $TypeofProblem {
            "Mechanical" {
                set solutionType [get_domnode_attribute [$domNode selectNodes "[spdAux::getRoute NumaMechanicalData]/value\[@n='SolutionType'\]"] v]
                set params [list analysis_type $solutionType TypeofProblem $TypeofProblem]
            }
            "Thermo-Mechanical" {
                set solutionType [get_domnode_attribute [$domNode selectNodes "[spdAux::getRoute "NumaThermo-MechanicalData"]/container\[@n='MechanicalPartProblem'\]/value\[@n='SolutionType'\]"] v]
                set params [list analysis_type $solutionType TypeofProblem $TypeofProblem]                    
            }
            default {
                error [= "Check type of problem"]
            }
        }
        if {[::Model::CheckElementsNodalCondition $conditionId $elemnames $params]} {return "normal"} else {return "hidden"}
    } {return "hidden"}
}

proc Numa::xml::ProcCheckConditionState {domNode args} {
    set cond_id [$domNode @n]
    
    # By active parts
    if {[spdAux::ProcActiveIfAnyPartState $domNode $args] eq "hidden"} {return "hidden"}
    
    # By dimension
    set resp [::Model::CheckConditionState $domNode]
    
    # By type of problem
    if {$resp} {
        set TypeofProblem [get_domnode_attribute [$domNode selectNodes [spdAux::getRoute NumaTypeofProblem]] v]
        if {$TypeofProblem ni [dict get [[Model::getCondition $cond_id] getAttributes ] TypeofProblem]} {set resp 0}
    }   
    
    if {$resp} {return "normal"} else {return "hidden"}
}


proc Numa::xml::ProcGetConstitutiveLaws {domNode args} {
    set Elementname [$domNode selectNodes {string(../value[@n='Element']/@v)}]
    set Claws [::Model::GetAvailableConstitutiveLaws $Elementname]
    #W "Round 1 $Claws"
    #foreach cl $Claws {W [$cl getName]}
    set type_of_problem [write::getValue NumaTypeofProblem]
    #W $type_of_problem
    set goodList [list ]
    #W "Pre type problem -> $type_of_problem"
    foreach cl $Claws {
        set type [$cl getAttribute Type]
        #W $type
        #W "cl -> [$cl getName]"
        #W "type -> $type"
        if {[string first "Therm" $type] eq -1 && $type_of_problem ni [list "Thermo-Mechanical"]} {
            lappend goodList $cl
        } elseif {[string first "Therm" $type] ne -1 && $type_of_problem in [list "Thermo-Mechanical"]} {
            lappend goodList $cl
        } 
    }
    #W "good $goodList"
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
    #W $values
    
    return $values
}

Numa::xml::Init

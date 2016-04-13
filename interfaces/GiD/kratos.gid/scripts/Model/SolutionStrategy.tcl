##################################################################################
#   This file is common for all Kratos Applications.
#   Do not change anything here unless it's strictly necessary.
##################################################################################

namespace eval Model {
# Clase Solution Strategey
catch {Scheme destroy}
oo::class create Scheme {
    superclass Entity
    
    variable elementfilters
    variable elementForceIn
    variable elementForceOut
    
    constructor {n} {
        next $n
        variable elementfilters
        variable elementForceIn
        variable elementForceOut
        
        set elementfilters [dict create]
        set elementForceIn [dict create]
        set elementForceOut [dict create]
    }
            
    method addElementFilter {efn efvl} {
        variable elementfilters
        dict set elementfilters $efn $efvl
    }
    method getElementFilters { } {
        variable elementfilters
        return $elementfilters
    }        
    method addElementForceIn {efn efvl} {
        variable elementForceIn
        dict set elementForceIn $efn $efvl
    }
    method getElementForceIn { } {
        variable elementForceIn
        return $elementForceIn
    }        
    method addElementForceOut {efn efvl} {
        variable elementForceOut
        dict set elementForceOut $efn $efvl
    }
    method getElementForceOut { } {
        variable elementForceOut
        return $elementForceOut
    }
    
}
# Clase Solution Strategey
catch {SolStrat destroy}
oo::class create SolStrat {
    superclass Entity
    
    variable solverEntries
    variable elementfilters
    variable schemes
    
    constructor {n} {
        next $n
        variable solverEntries
        variable schemes
        variable elementfilters
        
        set schemes [list ]
        set solverEntries [list ]
        set elementfilters [dict create]
    }
    
    
    method addSolverEntry {se} {
        variable solverEntries
        lappend solverEntries $se
    }
    method getSolversEntries { } {
        variable solverEntries
        return $solverEntries
    }
    method addScheme {se} {
        variable schemes
        
        lappend schemes $se
    }
    method getSchemes { } {
        variable schemes
        return $schemes
    }
    method getScheme { sid } {
        variable schemes
        
        set goodscheme ""
        foreach sc $schemes {
            if {[$sc getName] eq $sid} {set goodscheme $sc; break}
        }
        return $goodscheme
    }
    
    method addElementFilter {efn efvl} {
        variable elementfilters
        dict set elementfilters $efn $efvl
    }
    method getElementFilters { } {
        variable elementfilters
        return $elementfilters
    }
    
}
}

# Parsing
proc Model::ParseSolutionStrategies { doc } {
    variable SolutionStrategies
    
    set SolNodeList [$doc getElementsByTagName StrategyItem]
    foreach SolNode $SolNodeList {
        lappend SolutionStrategies [ParseSolNode $SolNode]
    }
}


proc Model::ParseSolNode { node } {
    set name [$node getAttribute n]
    
    set st [SolStrat new $name]
    $st setPublicName [$node getAttribute pn]
    
    foreach att [$node attributes] {
        $st setAttribute $att [$node getAttribute $att]
    }
    
    set paramListNode [list ]
    foreach nod [$node childNodes] {
        if {[$nod nodeName] eq "parameter_list"} {
            set paramListNode $nod
            break;
        }
    }
    if {[llength $paramListNode]} {
        foreach in [$paramListNode getElementsByTagName parameter] {
            set st [ParseInputParamNode $st $in]
        }
    }
    if {[llength [$node getElementsByTagName linearSolvers]]} {
        foreach se [[$node getElementsByTagName linearSolvers] getElementsByTagName linearSolverItem] {
            set st [ParseSolverEntry $st $se]
        }
    }
    #if {[llength [$node getElementsByTagName element_filters]]} {
    #    foreach ef [[$node getElementsByTagName element_filters] getElementsByTagName filter] {
    #        set st [ParseElementFilter $st $ef]
    #    }
    #}
    if {[llength [$node getElementsByTagName schemes]]} {
        foreach sc [[$node getElementsByTagName schemes] getElementsByTagName scheme] {
            set st [ParseScheme $st $sc]
        }
    }
    return $st
}

proc Model::ParseScheme {st scn} {
    
    set sc [Model::Scheme new [$scn @n]]
    $sc setPublicName [$scn @pn]
    
    if {[llength [$scn getElementsByTagName parameter_list]]} {
        foreach inn [[$scn getElementsByTagName parameter_list] getElementsByTagName parameter] {
            set sc [ParseInputParamNode $sc $inn]
            #W "in sch [$inn asXML]"
        }
    }
    if {[llength [$scn getElementsByTagName element_filters]]} {
        foreach ef [[$scn getElementsByTagName element_filters] getElementsByTagName filter] {
            set sc [ParseElementFilter $sc $ef]
        }
        foreach ef [[$scn getElementsByTagName element_filters] getElementsByTagName forceIn] {
            set sc [ParseElementFilter $sc $ef "in"]
        }
        foreach ef [[$scn getElementsByTagName element_filters] getElementsByTagName forceOut] {
            set sc [ParseElementFilter $sc $ef "out"]
        }
    }
    
    $st addScheme $sc
    return $st
}
proc Model::ParseElementFilter {st ef {forced ""}} {
    set n [$ef @field]
    set v [$ef @value]
    set values [split $v ","]
    
    if {$forced eq ""} {
        $st addElementFilter $n $values
    } elseif {$forced eq "in"} {
        $st addElementForceIn $n $values
    } else {
        $st addElementForceOut $n $values
    }
    return $st
}



# Getters
proc Model::GetSolutionStrategies { args } {
    variable SolutionStrategies
    
    if {$args eq ""} {return $SolutionStrategies}
    set cumplen [list ]
    foreach ss $SolutionStrategies {
        if {[$ss cumple $args]} { lappend cumplen $ss}
    }
    return $cumplen
}

proc Model::GetSolutionStrategy { id } {
    variable SolutionStrategies
    
    foreach ss $SolutionStrategies {
        if {[$ss getName] eq $id} { return $ss}
    }
    W "No solution strategy named $id"
}


proc Model::GetAvailableSchemes {solstrat} {
    set solst [Model::GetSolutionStrategy $solstrat]
    
    return [$solst getSchemes]
}

proc Model::GetAvailableSolvers {solstrat solverentryid} {
    variable Solvers
    
    set goodSolvers [list ]
    
    set solst [Model::GetSolutionStrategy $solstrat]
    foreach sentry [$solst getSolversEntries] {
        if {[$sentry getName] eq $solverentryid} {
            set filters [$sentry getAttributes]
            foreach solver $Solvers {
                if {[$solver cumple $filters]} {lappend goodSolvers $solver}
            }
        }
    }
    return $goodSolvers
}

proc Model::GetAllSolStratParams {} {
    variable SolutionStrategies
    
    set inputs [dict create ]
    foreach st $SolutionStrategies {
        foreach {k v} [$st getInputs] {
            dict set inputs $k $v
        }
    }
    return $inputs
}

proc Model::GetAllSchemeParams {} {
    variable SolutionStrategies
    
    set inputs [dict create ]
    foreach st $SolutionStrategies {
        foreach sc [$st getSchemes] {
            foreach {k v} [$sc getInputs] {
                dict set inputs $k $v
            }
        }
    }
    return $inputs
}

proc Model::GetAvailableElements {solutionStrategyId schemeId} { 
    variable Elements
    variable SolutionStrategies
    
    set cumplen [list ]
    set solst [Model::GetSolutionStrategy $solutionStrategyId]
    set scheme [$solst getScheme $schemeId]
    if {$scheme ne ""} {
        set filters [$scheme getElementFilters]
        set include [$scheme getElementForceIn]
        set exclude [$scheme getElementForceOut]
        
        foreach elem $Elements {
            set f [$elem cumple $filters]
            set i [$elem cumple $include]
            set i 0
            if {[llength $include]} {set i [$elem cumple $include]}
            set o 0
            if {[llength $exclude]} {set o [$elem cumple $exclude]}        
            if {[expr ($f && !$o) || $i]} { lappend cumplen $elem}
        }
    }
    return $cumplen
}

proc Model::GetAvailableConditions {solutionStrategyId schemeId} { 
    variable Conditions
    variable SolutionStrategies

    set cumplen [list ]
    #W $solutionStrategyId
    #W $schemeId
    set solst [Model::GetSolutionStrategy $solutionStrategyId]
    set scheme [$solst getScheme $schemeId]
    set filters [$scheme getElementFilters]
    
    foreach elem $Conditions {
        if {[$elem cumple $filters]} { lappend cumplen $elem} { lappend cumplen $elem}
    }
    
    return $cumplen
}


# State Check
proc Model::CheckSolStratInputState {SolStratName paramName} {
    set SolStrat [Model::GetSolutionStrategy $SolStratName]
    if {[$SolStrat getInputPn $paramName] ne ""} {return 1} {return 0}
}

proc Model::GetSolStratParamDep {SolStratName paramName} {
    set SolStrat [Model::GetSolutionStrategy $SolStratName]
    set in [$SolStrat getInputPn $paramName]
    set depN [$in getDepN]
    set depV [$in getDepV]
    return [list $depN $depV]
}

proc Model::CheckSchemeInputState {SolStratName SchemeName paramName} {
    set SolStrat [Model::GetSolutionStrategy $SolStratName]
    set Scheme [$SolStrat getScheme $SchemeName]
    if {$Scheme ne ""} {
        if {[$Scheme getInputPn $paramName] ne ""} {return 1} {return 0}
    } {return 0}
}
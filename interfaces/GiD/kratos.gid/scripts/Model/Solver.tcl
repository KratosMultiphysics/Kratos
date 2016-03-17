##################################################################################
#   This file is common for all Kratos Applications.
#   Do not change anything here unless it's strictly necessary.
##################################################################################

namespace eval Model {
# Clase Solution Strategey
catch {Solver destroy}
oo::class create Solver {
    superclass Entity
    
    constructor {n} {
        next $n
    }
}

# Clase Solution Strategey
catch {SolverEntry destroy}
oo::class create SolverEntry {
    superclass Entity
    
    constructor {n} {
        next $n
    }
}
}

proc Model::GetSolver { id } {
    variable Solvers
    
    foreach s $Solvers {
        if {[$s getName] eq $id} {return $s}
    }
    
}
proc Model::ParseSolvers { doc } {
    variable Solvers
    
    set SolNodeList [$doc getElementsByTagName solver]
    foreach SolNode $SolNodeList {
        lappend Solvers [ParseSolverNode $SolNode]
    }
}

proc Model::ParseSolverNode { node } {
    set name [$node getAttribute n]
    
    set sl [::Model::Solver new $name]
    $sl setPublicName [$node getAttribute pn]
    
    foreach attr [$node attributes] {
        $sl setAttribute $attr [$node getAttribute $attr]
    }
    
    foreach in [[$node getElementsByTagName inputs] getElementsByTagName parameter] {
        set sl [ParseInputParamNode $sl $in]
    }
    return $sl
}


proc Model::ParseSolverEntry {st sen} {
    set n [$sen @n]
    set pn [$sen @pn]
    
    set se [::Model::SolverEntry new $n]
    $se setPublicName $pn
    foreach f [$sen getElementsByTagName filter] {
        $se addAttribute [$f @field] [$f @value]    
    }
    $st addSolverEntry $se
    
    return $st
}


proc Model::GetAllSolversParams {} {
    variable Solvers
    
    set inputs [dict create ]
    foreach s $Solvers {
        foreach {k v} [$s getInputs] {
            dict set inputs $k $v
        }
    }
    return $inputs
}

proc Model::getSolverParamState {args} {
    variable Solvers
    lassign $args solvid inputid
    foreach solver $Solvers {
        if {[$solver getName] eq $solvid} {
            return [string compare [$solver getInputPn $inputid] ""]
        }
    }
    return 0
}
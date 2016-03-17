##################################################################################
#   This file is common for all Kratos Applications.
#   Do not change anything here unless it's strictly necessary.
##################################################################################

namespace eval Model {

# Clase CLaw
catch {CLaw destroy}
oo::class create CLaw {
    superclass Entity

    constructor {n} {
        next $n
    }
    
    method cumple {filtros} {
        set c [next $filtros]
         
        if {$c} {
            set ptdim $::Model::SpatialDimension
            set eldim [my getAttribute "Dimension"]
            if {$eldim ne ""} {
                if {$ptdim ne $eldim} {set c 0}
            }
        }
        
        return $c
    }
}

}

proc Model::ParseConstitutiveLaws { doc } {
    variable ConstitutiveLaws
    
    set SolNodeList [$doc getElementsByTagName CLaw]
    foreach SolNode $SolNodeList {
        lappend ConstitutiveLaws [ParseClawNode $SolNode]
    }
}

proc Model::ParseClawNode { node } {
    set name [$node getAttribute n]
    
    set cl [::Model::CLaw new $name]
    $cl setPublicName [$node getAttribute pn]
    $cl setHelp [$node getAttribute help]
    
    foreach attr [$node attributes] {
        $cl addAttribute $attr [$node getAttribute $attr]
    }
    
    foreach in [[$node getElementsByTagName inputs] getElementsByTagName parameter] {
        set $cl [ParseInputParamNode $cl $in]
    }
    foreach out [[$node getElementsByTagName outputs] getElementsByTagName parameter] {
        $cl addOutput [$out getAttribute n] [$out getAttribute pn]
    }
    return $cl
}


proc Model::GetConstitutiveLaws { args } {
    variable ConstitutiveLaws
    set cumplen [list ]
    foreach cl $ConstitutiveLaws {
        if {[$cl cumple $args]} { lappend cumplen $cl}
    }
    return $cumplen
}

proc Model::GetAllCLInputs {} {
    variable ConstitutiveLaws
    
    set inputs [dict create]
    foreach cl $ConstitutiveLaws {
        foreach in [dict keys [$cl getInputs]] {
            dict set inputs $in [$cl getInputPn $in]
        }
    }
    return $inputs
}

proc Model::GetAllCLOutputs {} {
    variable ConstitutiveLaws
    
    set outputs [dict create]
    foreach cl $ConstitutiveLaws {
        foreach in [dict keys [$cl getOutputs]] {
            dict set outputs $in [$cl getOutputPn $in]
        }
    }
    return $outputs
}

proc Model::CheckConstLawParamState {args} {
    variable ConstitutiveLaws
    
    lassign $args clid inputid
    foreach claw $ConstitutiveLaws {
        if {[$claw getName] eq $clid} {
            return [string compare [$claw getInputPn $inputid] ""]
        }
    }
    return 0
}


proc Model::CheckConstLawOutputState { constlawactive paramName} {
    variable ConstitutiveLaws
    
    set state 0
    foreach cl $ConstitutiveLaws {
       if {[$cl getName] in $constlawactive} {
           if {$paramName in [$cl getOutputs]} {
                set state 1
                break
           }
       }
    }
    return $state
}

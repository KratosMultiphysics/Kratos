##################################################################################
#   This file is common for all Kratos Applications.
#   Do not change anything here unless it's strictly necessary.
##################################################################################

namespace eval Model {
catch {Condition destroy}
oo::class create Condition {
    superclass Entity
    variable processName
    variable TopologyFeatures
    variable defaults
    
    constructor {n} {
        next $n
        variable TopologyFeatures
        set TopologyFeatures [list ]
        
        variable defaults
        set defaults [dict create ]
        
        variable processName
        set processName ""
    }
    
    method addTopologyFeature {top} {
        variable TopologyFeatures
        lappend TopologyFeatures $top
    }
    
    method getTopologyFeature {geo nod} {
        variable TopologyFeatures
        set ret ""
        foreach top $TopologyFeatures {
            if {[string match [$top getGeometry]* $geo] && [$top getNodes] eq $nod} {set ret $top; break}
        }
        return $ret
    }
    method getTopologyKratosName {geo nod} {
        set top [my getTopologyFeature $geo $nod]
        if {$top ne ""} {return [$top getKratosName]} {return ""}
    }

    method hasTopologyFeatures { } {
        variable TopologyFeatures
        set ret 0
        if {[llength $TopologyFeatures]} {set ret 1}
        return $ret
    }
    
    method setProcessName {pn} {
        variable processName
        set processName $pn
    }
    method getProcessName { } {
        variable processName
        return $processName
    }
    method setDefault {itemName itemField itemValue} {
        variable defaults
        dict set defaults $itemName $itemField $itemValue
    }
    method getDefault {itemName itemField } {
        variable defaults
        set ret ""
        if {[dict exists $defaults $itemName $itemField]} {
            set ret [dict get $defaults $itemName $itemField]
        }
        return $ret
    }
    method getDefaults {itemName} {
        variable defaults
        if {[dict exists $defaults $itemName]} {
            return [dict keys [dict get $defaults $itemName] ]
        }
    }
    
}
}
proc Model::ForgetConditions { } {
    variable Conditions
    set Conditions [list ]
}
proc Model::ForgetCondition { cnd_id } {
    variable Conditions
    set Conditions2 [list ]
    foreach cnd $Conditions {
        if {[$cnd getName] ne $cnd_id} {
            lappend Conditions2 $cnd
        }
    }
    set Conditions $Conditions2
}
proc Model::GetConditions {args} { 
    variable Conditions
    if {$args eq ""} {
        return $Conditions
    }
    
    set cumplen [list ]
    foreach elem $Conditions {
        if {[$elem cumple {*}$args]} { lappend cumplen $elem}
    }
    return $cumplen
}


proc Model::ParseConditions { doc } {
    variable Conditions
    
    set CondNodeList [$doc getElementsByTagName ConditionItem]
    foreach CondNode $CondNodeList {
        lappend Conditions [ParseCondNode $CondNode]
    }
}


proc Model::ParseCondNode { node } {
    set name [$node getAttribute n]
    set cnd [::Model::Condition new $name]
    $cnd setPublicName [$node getAttribute pn]
    if {[$node hasAttribute help]} {$cnd setHelp [$node getAttribute help]}
    if {[$node hasAttribute ProcessName]} {$cnd setProcessName [$node getAttribute ProcessName]}
    
    foreach att [$node attributes] {
        $cnd setAttribute $att [split [$node getAttribute $att] ","]
        #W "$att : [$el getAttribute $att]"
    }
    set topology_base [$node getElementsByTagName TopologyFeatures]
    if {[llength $topology_base] eq 1} {
        foreach top [$topology_base getElementsByTagName item]  {
            set cnd [ParseTopologyNode $cnd $top]
        }
    }
    set defaults_base [$node getElementsByTagName DefaultValues]
    if {[llength $defaults_base] eq 1} {
        foreach def [$defaults_base getElementsByTagName value]  {
            set itemName [$def @n]
            foreach att [$def attributes] {
                if {$att ne "n"} {
                    set itemField $att
                    set itemValue [$def getAttribute $att]
                    $cnd setDefault $itemName $itemField $itemValue
                }
            }
        }
    }
    set inputs_base [$node getElementsByTagName inputs]
    if {[llength $inputs_base] eq 1} {
        foreach in [$inputs_base getElementsByTagName parameter]  {
            set cnd [ParseInputParamNode $cnd $in]
        }
    }
    set outputs_base [$node getElementsByTagName outputs]
    if {[llength $outputs_base] eq 1} {
        foreach out [$outputs_base getElementsByTagName parameter] {
            set n [$out @n]
            set pn [$out @pn]
            $cnd addOutput $n $pn
        }
    }
    return $cnd
}

proc Model::getCondition {cid} { 
    variable Conditions

    foreach elem $Conditions {
        #W "Checking  [$elem getName] -> $cid"
        if {[$elem getName] eq $cid} { return $elem}
    }
    return ""
}


proc Model::GetAllCondOutputs {} {
    variable Conditions
    
    set outputs [dict create]
    foreach el $Conditions {
        foreach in [dict keys [$el getOutputs]] {
            dict set outputs $in [$el getOutputPn $in]
        }
    }
    return $outputs
}

proc Model::GetAllCondInputs {} {
    variable Conditions
    
    set inputs [dict create]
    foreach el $Conditions {
        foreach in [dict keys [$el getInputs]] {
            dict set inputs $in [$el getInputPn $in]
        }
    }
    return $inputs
}

proc Model::CheckConditionState {node args} {
    variable Conditions
    set condid [get_domnode_attribute $node n]
    set cumple 1
    foreach cond $Conditions {
        if {[$cond getName] eq $condid} {
            if {![$cond cumple {*}$args]} {
                set cumple 0
                break
            } {
                set ptdim $::Model::SpatialDimension
                if {[$cond getAttribute "WorkingSpaceDimension"] ne ""} {
                    set eldim [split [$cond getAttribute "WorkingSpaceDimension"] ","]
                    if {$ptdim ni $eldim} {set cumple 0; break}
                }
             
            }
        }
    }
    return $cumple
}


proc Model::CheckCondOutputState {elemsactive paramName} {
    variable Conditions

    set state 0
    foreach el $Conditions {
       if {[$el getName] in $elemsactive} {
           if {$paramName in [$el getOutputs]} {
                set state 1
                break
           }
       }
    }
    return $state
}


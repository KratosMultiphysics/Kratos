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
        catch {
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
    $cnd setHelp [$node getAttribute help]
    $cnd setProcessName [$node getAttribute ProcessName]
    
    foreach att [$node attributes] {
        $cnd setAttribute $att [split [$node getAttribute $att] ","]
        #W "$att : [$el getAttribute $att]"
    }
    catch {
        foreach top [[$node getElementsByTagName TopologyFeatures] getElementsByTagName item]  {
            set cnd [ParseTopologyNode $cnd $top]
        }
    }
    catch {
        foreach def [[$node getElementsByTagName DefaultValues] getElementsByTagName value]  {
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
    catch {
        set inputsNode [$node getElementsByTagName inputs]
        foreach in [$inputsNode getElementsByTagName parameter]  {
            set cnd [ParseInputParamNode $cnd $in]
        }
    }
    catch {
        foreach out [[$node getElementsByTagName outputs] getElementsByTagName parameter] {
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

proc Model::CheckConditionState {node} {
    variable Conditions
    set condid [get_domnode_attribute $node n]
    set cumple 1
    foreach cond $Conditions {
        if {[$cond getName] eq $condid} {
            if {![$cond cumple ""]} {
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


namespace eval Model {
catch {Element destroy}
oo::class create Element {
    superclass Entity
    
    variable TopologyFeatures
    variable ElementNodalCondition
    variable ConstLawFilters
    
    constructor {n} {
        next $n
        variable TopologyFeatures
        variable ElementNodalCondition
        variable ConstLawFilters
        
        set TopologyFeatures [list ]
        set ElementNodalCondition [dict create]
        set constLawFilters [dict create]
    }
    
    method addConstLawFilter {key value} {variable constLawFilters; dict set constLawFilters $key $value}
    method getConstLawFilterValue {key} {
        variable constLawFilters;
        
        set v ""
        catch {
            set v [dict get $constLawFilters $key]
        }
        return $v
    }
    method getConstLawFilters {} {
        variable constLawFilters
        
        set v [list ]
        if {[info exists constLawFilters]} {
            foreach k [dict keys $constLawFilters] {
                lappend v $k [dict get $constLawFilters $k] 
            }
        }
        return $v
    }
    
    method addNodalCondition {key nc} {
        variable ElementNodalCondition
        dict set ElementNodalCondition $key $nc
    }
    method getNodalConditions {} {
        variable ElementNodalCondition
        return $ElementNodalCondition
    }
    method getNodalCondition {key} {
        variable ElementNodalCondition
        set v ""
        catch {
            set v [dict get $ElementNodalCondition $key]
        }
        return $v
    }
    
    
    
    method addTopologyFeature {top} {
        variable TopologyFeatures
        lappend TopologyFeatures $top
    }
    method getTopologyFeature {geo nod} {
        variable TopologyFeatures
        set ret ""
        foreach top $TopologyFeatures {            
            if {[$top getNodes] eq $nod} {set ret $top; break}
            #if {[$top getGeometry] eq $geo && [$top getNodes] eq $nod} {set ret $top; break}
        }
        return $ret
    }
    method cumple {filtros} {
        set c [next $filtros]
         
        if {$c} {
            set ptdim $::Model::SpatialDimension
            set eldim [my getAttribute "WorkingSpaceDimension"]
            append eldim "D"
            if {$ptdim ne $eldim} {set c 0}
        }
        
        return $c
    }
}
catch {NodalCondition destroy}
oo::class create NodalCondition {
    superclass Parameter
    variable reaction
    variable ov
    
    constructor {n pn type} {
        next $n $pn $type "1" "" "" $pn
        set reaction ""
        set ov "point,line,surface,volume"
    }
    
    method setReaction {r} {variable reaction; set reaction $r}
    method getReaction {} {variable reaction; return $reaction}
    method setOv {r} {variable ov; set ov $r}
    method getOv {} { variable ov; return $ov}
    
}
}

proc Model::ParseElements { doc } {
    variable Elements
    
    set ElemNodeList [$doc getElementsByTagName ElementItem]
    foreach ElemNode $ElemNodeList {
        lappend Elements [ParseElemNode $ElemNode]
    }
}


proc Model::ParseElemNode { node } {
    set name [$node getAttribute n]
    
    set el [::Model::Element new $name]
    $el setPublicName [$node getAttribute pn]
    
    foreach att [$node attributes] {
        $el setAttribute $att [split [$node getAttribute $att] ","]
        #W "$att : [$el getAttribute $att]"
    }
    
    foreach top [[$node getElementsByTagName TopologyFeatures] childNodes]  {
        set el [ParseTopologyNode $el $top]
    }
    foreach in [[$node getElementsByTagName inputs] childNodes]  {
        set el [ParseInputParamNode $el $in]
    }
    foreach out [[$node getElementsByTagName outputs] childNodes] {
        set n [$out @n]
        set pn [$out @pn]
        $el addOutput $n $pn
    }
    foreach clf [[$node getElementsByTagName ConstitutiveLaw_FilterFeatures] childNodes]  {
        $el addConstLawFilter [$clf getAttribute field] [$clf getAttribute value]
        #W "CL Filter [$clf nodeName] -> [$el getConstLawFilterValue [$clf nodeName]]"
    }
    foreach ncnode [[$node getElementsByTagName NodalConditions] childNodes]  {
        set n [$ncnode @n]
        set nc [::Model::NodalCondition new $n [$ncnode @pn] [$ncnode @type] ]
        $nc setReaction [$ncnode @reaction]
        if {[$ncnode hasAttribute ov]} {$nc setOv [$ncnode @ov]}
            
        set fi "0"
        catch {set fi [$ncnode @fixity]}
        $nc setFixity $fi
        $nc setUnits [$ncnode @units]
        $nc setUnitMagnitude [$ncnode @unit_magnitude]
        foreach att [$ncnode attributes] {
            $nc addAttribute $att [split [$ncnode getAttribute $att] ","]
        }
        
        $el addNodalCondition $n $nc
    }
    
    return $el
}

# Se usa?
proc Model::GetElements {args} { 
    variable Elements
    #W "Get elements $args"
    set cumplen [list ]
    foreach elem $Elements {
        if {[$elem cumple $args]} { lappend cumplen $elem}
    }
    return $cumplen
}

proc Model::getElement {eid} { 
    variable Elements

    foreach elem $Elements {
        if {[$elem getName] eq $eid} { return $elem}
    }
    return ""
}

proc Model::GetConstitutiveLawFilters {ename} {
    variable Elements
    
    foreach elem $Elements {
        if {[$elem getName] eq $ename} { return [$elem getConstLawFilters]}
    }
}

proc Model::GetAvailableConstitutiveLaws {elementId} { 
    variable ConstitutiveLaws

    set cumplen [list ]
    
    set filters [GetConstitutiveLawFilters $elementId]
    #W "filtros $filters"
    foreach cl $ConstitutiveLaws {
        #W "Cumple [$cl getName]? [$cl cumple $filters]"
        if {[$cl cumple $filters]} { lappend cumplen $cl}
    }
    
    return $cumplen
}

proc Model::GetAllElemOutputs {} {
    variable Elements
    
    set outputs [dict create]
    foreach el $Elements {
        foreach in [dict keys [$el getOutputs]] {
            dict set outputs $in [$el getOutputPn $in]
        }
    }
    return $outputs
}

proc Model::GetAllElemInputs {} {
    variable Elements
    
    set inputs [dict create]
    foreach el $Elements {
        foreach in [dict keys [$el getInputs]] {
            dict set inputs $in [$el getInputPn $in]
        }
    }
    return $inputs
}

proc Model::getAllDOFs {} {
    variable Elements
    
    set dofs [dict create]
    foreach el $Elements {
        foreach in [dict keys [$el getNodalConditions]] {
            dict set dofs $in [$el getNodalCondition $in]
        }
    }
    return $dofs
}

proc Model::getDOFbyId {dofid} {
    return [dict get [getAllNodalConditions] $dofid]
}


proc Model::CheckElemState {elid inputid} {
    variable Elements
    foreach elem $Elements {
        if {[$elem getName] eq $elid} {
            return [string compare [$elem getInputPn $inputid] ""]
        }
    }
    return 0
}

proc Model::CheckElemParamState {node} {
    set id [$node getAttribute n]
    set elem [get_domnode_attribute [[$node parent] selectNodes "./value\[@n='Element'\]"] v]
    if {$elem eq ""} {return 0}

    return [CheckElemState $elem $id]
}

proc Model::CheckElementOutputState {elemsactive paramName} {
    variable Elements

    set state 0
    foreach el $Elements {
       if {[$el getName] in $elemsactive} {
           if {$paramName in [$el getOutputs]} {
                set state 1
                break
           }
       }
    }
    return $state
}

proc Model::CheckElementsCondition {conditionId elemnames {restrictions "" }} {
    set ret 0
    if {[llength $elemnames] < 1} {
        #
    } else {
        foreach eid $elemnames {
            set elem [getElement $eid]
            set dof [$elem getNodalCondition $conditionId]
            
            if {$dof ne ""} {
                set ret 1
                foreach {key value} $restrictions {
                    # JG: Revisar bidireccionalidad
                    if {$value ni [$dof getAttribute $key]} {set ret 0;break}
                }
            }
            
        }
    }
    return $ret
}

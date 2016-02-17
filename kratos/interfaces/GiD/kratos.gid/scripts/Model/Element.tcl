namespace eval Model {
catch {Element destroy}
oo::class create Element {
    superclass Entity
    
    variable TopologyFeatures
    variable ElementDegreeOfFreedom
    variable ConstLawFilters
    
    constructor {n} {
        next $n
        variable TopologyFeatures
        variable ElementDegreeOfFreedom
        variable ConstLawFilters
        
        set TopologyFeatures [list ]
        set ElementDegreeOfFreedom [dict create]
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
    
    method addDOF {key dof} {
        variable ElementDegreeOfFreedom
        dict set ElementDegreeOfFreedom $key $dof
    }
    method getDOFs {} {
        variable ElementDegreeOfFreedom
        return $ElementDegreeOfFreedom
    }
    method getDOF {key} {
        variable ElementDegreeOfFreedom
        set v ""
        catch {
            set v [dict get $ElementDegreeOfFreedom $key]
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
catch {DOF destroy}
oo::class create DOF {
    superclass Parameter
    variable reaction
    
    constructor {n pn type} {
        next $n $pn $type "1" "" "" $pn
        set reaction ""
    }
    
    method setReaction {r} {variable reaction; set reaction $r}
    method getReaction {} {variable reaction; return $reaction}
    
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
    foreach dofnode [[$node getElementsByTagName BoundaryConditions] childNodes]  {
        set n [$dofnode @n]
        set dof [::Model::DOF new $n [$dofnode @pn] [$dofnode @type] ]
        $dof setReaction [$dofnode @reaction]
        $dof setUnits [$dofnode @units]
        $dof setUnitMagnitude [$dofnode @unit_magnitude]
        foreach att [$dofnode attributes] {
            $dof addAttribute $att [split [$dofnode getAttribute $att] ","]
        }
        
        $el addDOF $n $dof
    }
    
    return $el
}

# Se usa?
proc Model::GetElements {args} { 
    variable Elements
    W "Get elements $args"
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
        foreach in [dict keys [$el getDOFs]] {
            dict set dofs $in [$el getDOF $in]
        }
    }
    return $dofs
}

proc Model::getDOFbyId {dofid} {
    return [dict get [getAllDOFs] $dofid]
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

proc Model::CheckElementsCondition {conditionId elemnames solType} {
    set ret 0
    set aparece 0
    set coincideAnalysis 0
    if {[llength $elemnames] < 1} {
        set ret 0
    } else {
        foreach eid $elemnames {
            set elem [getElement $eid]
            set dof [$elem getDOF $conditionId]
            
            if {$dof ne ""} {
                set aparece 1
                if {[$dof getAttribute analysis_type] eq "Dynamic" && $solType eq "Dynamic"} {set coincideAnalysis 1;break}
                if {[$dof getAttribute analysis_type] eq "Static"} {set coincideAnalysis 1;break}
            }
            
        }
    }
    if {[expr $aparece && $coincideAnalysis]} {set ret 1}
    return $ret
}

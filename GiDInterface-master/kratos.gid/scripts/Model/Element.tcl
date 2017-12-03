##################################################################################
#   This file is common for all Kratos Applications.
#   Do not change anything here unless it's strictly necessary.
##################################################################################

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
        set ElementNodalCondition [list ]
        set constLawFilters [dict create]
    }
    
    method addConstLawFilter {key value} {variable constLawFilters; dict set constLawFilters $key $value}
    method getConstLawFilterValue {key} {
        variable constLawFilters;
        
        set v ""
        if {[dict exists $constLawFilters $key]} {
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
    
    method addNodalCondition {nc_name} {
        variable ElementNodalCondition
        lappend ElementNodalCondition $nc_name
    }
    method getNodalConditions {} {
        variable ElementNodalCondition
        set nclist [list ]
        foreach nc [::Model::getAllNodalConditions] {
            if {[$nc getName] in $ElementNodalCondition} {
                lappend nclist $nc
            }
        }
        return $nclist
    }
    method getNodalCondition {key} {
        variable ElementNodalCondition
        set v ""
        foreach nc [::Model::getAllNodalConditions] {
            if {[$nc getName] in $ElementNodalCondition} {
                return $nc
            }
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
    method cumple {args} {
        set c [next {*}$args]
         
        if {$c} {
            set ptdim $::Model::SpatialDimension
            set eldim [split [my getAttribute "WorkingSpaceDimension"] ","]
            if {$ptdim ni $eldim} {set c 0}
        }
        
        return $c
    }
}
catch {NodalCondition destroy}
oo::class create NodalCondition {
    superclass Condition
    variable reaction
    variable ov
    
    constructor {n} {
        next $n
        set reaction ""
        set ov ""
        #set ov "point,line,surface,volume"
    }
    
    method setReaction {r} {variable reaction; set reaction $r}
    method getReaction {} {variable reaction; return $reaction}
    method setOv {r} {variable ov; set ov $r}
    method getOv {} { variable ov; return $ov}
    
}
}

proc Model::ForgetElements { } {
    variable Elements
    set Elements [list ]
}
proc Model::ForgetElement { elem_id } {
    variable Elements
    set Elements2 [list ]
    foreach elem $Elements {
        if {[$elem getName] ne $elem_id} {
            lappend Elements2 $elem
        }
    }
    set Elements $Elements2
}

proc Model::ParseElements { doc } {
    variable Elements
    
    set ElemNodeList [$doc getElementsByTagName ElementItem]
    foreach ElemNode $ElemNodeList {
        lappend Elements [ParseElemNode $ElemNode]
    }
}
proc Model::ParseNodalConditions { doc } {
    variable NodalConditions

    set NCList [$doc getElementsByTagName NodalConditionItem]
    foreach Node $NCList {
        lappend NodalConditions [ParseNodalConditionsNode $Node]
    }
}

proc Model::ForgetNodalConditions { } {
    variable NodalConditions
    set NodalConditions [list ]
}
proc Model::ForgetNodalCondition { cnd_id } {
    variable NodalConditions
    set NodalConditions2 [list ]
    foreach cnd $NodalConditions {
        if {[$cnd getName] ne $cnd_id} {
            lappend NodalConditions2 $cnd
        }
    }
    set NodalConditions $NodalConditions2
}

proc Model::ParseElemNode { node } {
    set name [$node getAttribute n]
    
    set el [::Model::Element new $name]
    $el setPublicName [$node getAttribute pn]
    if {[$node hasAttribute v]} {$el setDv [$node getAttribute v]}
    
    foreach att [$node attributes] {
        $el setAttribute $att [split [$node getAttribute $att] ","]
        #W "$att : [$el getAttribute $att]"
    }
    set top_node [$node getElementsByTagName TopologyFeatures]
    if {$top_node ne "" && [$top_node hasChildNodes]} {
        foreach top [$top_node childNodes]  {
            set el [ParseTopologyNode $el $top]
        }
    }
    set inputs_node [$node getElementsByTagName inputs]
    if {$inputs_node ne "" && [$inputs_node hasChildNodes]} {
        foreach in [$inputs_node childNodes]  {
            set el [ParseInputParamNode $el $in]
        }
    }
    set outputs_node [$node getElementsByTagName outputs]
    if {$outputs_node ne "" && [$outputs_node hasChildNodes]} {
        foreach out [$outputs_node childNodes] {
            set n [$out @n]
            set pn [$out @pn]
            set v false
            if {[$out hasAttribute v]} {set v [$out @v]}
            set outobj [::Model::Parameter new $n $pn bool $v "" "" "" ]
            $el addOutputDone $outobj
        }
    }
    set const_laws_node [$node getElementsByTagName ConstitutiveLaw_FilterFeatures]
    if {$const_laws_node ne "" && [$const_laws_node hasChildNodes]} {
        foreach clf [$const_laws_node childNodes]  {
            $el addConstLawFilter [$clf getAttribute field] [$clf getAttribute value]
        }
    }
    set nodal_conds_node [$node getElementsByTagName NodalConditions]
    if {$nodal_conds_node ne "" && [$nodal_conds_node hasChildNodes]} {
        foreach ncnode [$nodal_conds_node childNodes]  {
            set n [$ncnode @n]
            $el addNodalCondition $n
        }
    }
    
    return $el
}
proc Model::ParseNodalConditionsNode { node } {
    set name [$node getAttribute n]
    #W "Parsing $name"
    set el [::Model::NodalCondition new $name]
    $el setPublicName [$node getAttribute pn]
    
    if {[$node hasAttribute ov]} {
        set ov [$node getAttribute ov]
        $el setOv $ov
    }
    
    
    foreach att [$node attributes] {
        $el setAttribute $att [split [$node getAttribute $att] ","]
        #W "$att : [$el getAttribute $att]"
    }
    set inputNodes [$node getElementsByTagName inputs]
    if {$inputNodes ne ""} {
        foreach in [$inputNodes getElementsByTagName parameter] {
            set el [ParseInputParamNode $el $in]
        }
    }
    
    set outputNodes [$node getElementsByTagName outputs]
    if {$outputNodes ne ""} {
        foreach out [$outputNodes getElementsByTagName parameter] {
            set n [$out @n]
            set pn [$out @pn]
            set v false
            if {[$out hasAttribute v]} {set v [$out @v]}
            set outobj [::Model::Parameter new $n $pn bool $v "" "" "" ]
            $el addOutputDone $outobj
        }
    }
    $el setProcessName [$node getAttribute ProcessName]
    set defaultValuesList [$node getElementsByTagName DefaultValues]
    if {[llength $defaultValuesList]} {
        foreach def [$defaultValuesList getElementsByTagName value]  {
            set itemName [$def @n]
            foreach att [$def attributes] {
                if {$att ne "n"} {
                    set itemField $att
                    set itemValue [$def getAttribute $att]
                    $el setDefault $itemName $itemField $itemValue
                }
            }
        }
    }
    #W "[$el getName]"
    return $el
}

proc Model::GetElements {args} { 
    variable Elements
    #W "Get elements $args"
    set cumplen [list ]
    foreach elem $Elements {
        if {[$elem cumple {*}$args]} { lappend cumplen $elem}
    }
    #W "Elementos buenos $cumplen"
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

proc Model::GetAllElemOutputs {args} {
    set outputs [dict create]
    foreach el [GetElements {*}$args] {
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

proc Model::getAllNodalConditions {} {
    variable NodalConditions
    #W "all"
    return $NodalConditions
}

proc Model::GetNodalConditions {args} {
    #return [Model::getAllNodalConditions]
    variable NodalConditions
    set validlist [list ]
    
    foreach nc $NodalConditions {
        #W "ASK [$nc cumple {*}$args]"
        if {[$nc cumple {*}$args]} {lappend validlist $nc}
    }
    #W "validas $validlist"
    return $validlist
}

proc Model::getNodalConditionbyId {ncid} {
    set ret ""
    foreach nc [getAllNodalConditions] {
        if {[$nc getName] eq $ncid} {set ret $nc; break}
    }
    return $ret
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

proc Model::CheckElementsNodalCondition {conditionId elemnames {restrictions "" }} {
    set ret 0
    if {[llength $elemnames] < 1} {
        #
    } else {
        foreach eid $elemnames {
            if {$eid ne ""} {
                set elem [getElement $eid]
                foreach elemNCNode [$elem getNodalConditions] {
                    set elemNC [$elemNCNode getName]
                    if {$elemNC eq $conditionId} {
                        set ret 1
                        foreach {key value} $restrictions {
                            # JG: Revisar bidireccionalidad
                            set nc [getNodalConditionbyId $conditionId]
                            if {$value ni [$nc getAttribute $key]} {set ret 0;break}
                        }
                    }
                }
            }
        }
    }
    
    return $ret
}
proc Model::CheckNodalConditionOutputState {conditionId outputId {restrictions "" }} {
    set ret 0
    #W "Con $conditionId out $outputId"
    set nc [getNodalConditionbyId $conditionId]
    foreach {name output} [$nc getOutputs] {
        if {$name eq $outputId} {
            set ret 1
            foreach {key value} $restrictions {
                # JG: Revisar bidireccionalidad
                if {$value ni [$nc getAttribute $key]} {set ret 0;break}
            }
        }
    }
    
    return $ret
}

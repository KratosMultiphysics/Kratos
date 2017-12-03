##################################################################################
#   This file is common for all Kratos Applications.
#   Do not change anything here unless it's strictly necessary.
##################################################################################

# Clase Parameter
namespace eval Model {
catch {Parameter destroy}
oo::class create Parameter {
    superclass Entity
    
    variable type
    variable dv
    variable depname
    variable depv
    
    variable actualize
    
    variable values
    variable pvalues
    variable units
    variable unitMagnitude
    
    constructor {n pn t v u um h {vs ""} {pvs "" }} {
        next $n
        variable type
        variable dv
        variable depname
        variable depv
        variable values
        variable pvalues
        
        variable units
        variable unitMagnitude
        
        variable actualize
        
        my setPublicName $pn
        set type $t
        set dv $v
        my setHelp $h
        
        set depname ""
        set depv ""
        
        set units $u
        set unitMagnitude $um
        
        set values [list ]
        set pvalues [list ]
        
        if {$vs ne ""} {
            foreach val [split $vs ,] {
                lappend values $val
            }
        }
        if {$pvs ne ""} {    
            set pvalues [list ]
            foreach val [split $pvs  ,] {
                lappend pvalues $val
            }
        } else {set pvalues $values}
        
        set actualize 0
    }
    
    method setDep {dn dv} {
        variable depname
        variable depv
        
        set depname $dn
        set depv $dv

    }
    
    method getType { } {variable type; return $type}
    method getDv { } {variable dv; return $dv}
    method getValues { } {variable values; return $values}
    method getPValues { } {variable pvalues; return $pvalues}
    
    method setUnits {u} {variable units; set units $u}
    method getUnits { } {variable units; return $units}
    method setUnitMagnitude {u} {variable unitMagnitude; set unitMagnitude $u}
    method getUnitMagnitude { } {variable unitMagnitude; return $unitMagnitude}
    
    method getDepN {} {variable depname; return $depname}
    method getDepV {} {variable depv; return $depv}
    
    method setActualize {u} {variable actualize; set actualize $u}
    method getActualize { } {variable actualize; return $actualize}
}
}

proc Model::ParseInputParamNode {parent_object in} {
    set n [$in @n]
    set pn [$in @pn]
    set t "double"
    if {[$in hasAttribute type]} {set t [$in @type]}
    set v 0
    if {[$in hasAttribute v]} {set v [$in @v]}
    set h ""
    if {[$in hasAttribute help]} {set h [$in @help]}
    set vs ""
    set pvs ""
    if {[$in hasAttribute values]} {set vs [$in @values]}
    if {[$in hasAttribute pvalues]} {set pvs [$in @pvalues]}
    set units ""
    set unitMagnitude ""
    if {[$in hasAttribute units]} {set units [$in @units]}
    if {[$in hasAttribute unit_magnitude]} {set unitMagnitude [$in @unit_magnitude]}
    set input [::Model::Parameter new $n $pn $t $v $units $unitMagnitude $h $vs $pvs]
    if {[$in hasChildNodes]} {$input setActualize 1}
    foreach att [$in attributes] {
        #W "$n $att"
        $input addAttribute $att [$in getAttribute $att]
    }
    $parent_object addInputDone $input
    if {[$in hasAttribute "parent"]} {
        set dn [[$in parent] @n]
        set dv [$in @parent]
        $parent_object addInputDependency $n $dn $dv
    }
    
    return $parent_object
}

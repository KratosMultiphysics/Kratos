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
    
    variable fix
    
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
        variable fix
        
        variable actualize
        
        my setPublicName $pn
        set type $t
        set dv $v
        my setHelp $h
        #if {$h ne ""} {W "$n $h [my getHelp]"}
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
        
        set fix 0
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
    method setFixity {u} {variable fix; set fix $u}
    method getFixity { } {variable fix; return $fix}
    
    method getDepN {} {variable depname; return $depname}
    method getDepV {} {variable depv; return $depv}
    
    
    method setActualize {u} {variable actualize; set actualize $u}
    method getActualize { } {variable actualize; return $actualize}
}
}

proc Model::ParseInputParamNode {st in} {
    set n [$in @n]
    set pn [$in @pn]
    set t "double"
    catch {set t [$in @type]}
    set v 0
    catch {set v [$in @v]}
    set h ""
    catch {set h [$in @help]}
    set vs ""
    set pvs ""
    catch {set vs [$in @values]}
    catch {set pvs [$in @pvalues]}
    set units ""
    set unitMagnitude ""
    catch {set units [$in @units]}
    catch {set unitMagnitude [$in @unit_magnitude]}
    set fi ""
    catch {set fi [$in @fixity]}
    set input [::Model::Parameter new $n $pn $t $v $units $unitMagnitude $h $vs $pvs]
    if {$fi eq ""} {set fi 0}
    $input setFixity $fi
    if {[$in hasChildNodes]} {$input setActualize 1}
    $st addInputDone $input
    if {[$in hasAttribute "parent"]} {
        set dn [[$in parent] @n]
        set dv [$in @parent]
        $st addInputDependency $n $dn $dv
    }
    return $st
}

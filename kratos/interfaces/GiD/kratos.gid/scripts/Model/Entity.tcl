##################################################################################
#   This file is common for all Kratos Applications.
#   Do not change anything here unless it's strictly necessary.
##################################################################################

namespace eval Model {
catch {Entity destroy}
oo::class create Entity {
    variable name
    variable publicname
    variable inputs
    variable outputs
    variable attributes
    variable help
    
    constructor {n} {
        variable name
        variable publicname
        variable inputs
        variable outputs
        variable attributes
        variable help
        
        set name $n
        set publicname $n
        set inputs [dict create]
        set outputs [dict create]
        set attributes [dict create]
        set help ""
    }
    
    method setAttribute {att val} {variable attributes; dict set attributes $att $val}
    method getAttribute {att} {
        variable attributes
        set v ""
        catch {
            set v [dict get $attributes $att]
        }
        return $v
    }
    method addAttribute {key values} {variable attributes; dict set attributes $key $values}
    method getAttributes {} {variable attributes; return $attributes}
    
    method getName { } {variable name; return $name}
    method getPublicName { } {variable publicname; return $publicname}
    method setPublicName { pn } {variable publicname; set publicname $pn}
    
    method getHelp { } {variable help; return $help}
    method setHelp { h } {variable help; set help $h}

    method addInputDone {inp} {
        variable inputs
        dict set inputs [$inp getName] $inp
    }
    method addInput {args} {
        variable inputs
        lassign $args n pn t v h vs pvs
        dict set inputs $n [::Model::Parameter new $n $pn $t $v $h $vs $pvs]
    }
    
    method addOutputDone {out} {
        variable outputs
        dict set outputs [$out getName] $out
    }
    method addOutput {args} {
        variable outputs
        lassign $args n pn t v h vs pvs
        dict set outputs $n [::Model::Parameter new $n $pn $t $v $h $vs $pvs]
    }
        
    method getInputs { } {variable inputs; return $inputs}
    method getInputPn {in} {
        variable inputs
        set pn ""
        catch {
            set pn [dict get $inputs $in]
        }
        return $pn
    }
    
    method getOutputs { } {variable outputs; return $outputs}
    method getOutputPn {in} {
        variable outputs
        set pn ""
        catch {
            set pn [dict get $outputs $in]
        }
        return $pn
    }
    

    method addInputDependency {n dn dv} {
        variable inputs
        [dict get $inputs $n] setDep $dn $dv
    }
    
    
    method cumple {filtros} {
        set c 1
        
        foreach {k listfiltervalues} $filtros {
            set listattributesvalues [my getAttribute $k]
            set b1 [string first $listfiltervalues $listattributesvalues]
            set b2 [string first $listattributesvalues $listfiltervalues] 
            set b3 1
            if {$b1 eq -1 && $b2 eq -1} {
                set c 0
                break
            }
        }
        #if {$c eq 1 } {W "Cumples"}
        return $c
    }
}
}
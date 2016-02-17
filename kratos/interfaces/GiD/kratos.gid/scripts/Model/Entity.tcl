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

    
    method addInput {args} {
        variable inputs
        lassign $args n pn t v h vs pvs
        dict set inputs $n [::Model::Parameter new $n $pn $t $v $h $vs $pvs]
    }
    
    method addOutput {args} {
        variable outputs
        lassign $args n pn t v h vs pvs
        dict set outputs $n [::Model::Parameter new $n $pn $t $v $h $vs $pvs]
    }
    #method addOutput {n pn} {variable outputs; dict set outputs $n $pn}
        
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
            #W "$k : $listfiltervalues -> $listattributesvalues"
            set b1 [string first $listfiltervalues $listattributesvalues]
            set b2 [string first $listattributesvalues $listfiltervalues] 
            set b3 1
            #W "$b1 $b2"
            if {$b1 eq -1 && $b2 eq -1} {
                #W "fail $k"
                set c 0
                break
            }
        }
        #if {$c eq 1 } {W "Cumples"}
        return $c
    }
}
}
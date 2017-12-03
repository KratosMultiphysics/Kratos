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
    
    #method setAttribute {att val} {variable name; W "$name -> adding $att : $val"; variable attributes; dict set attributes $att $val}
    method setAttribute {att val} {variable attributes; dict set attributes $att $val}
    method hasAttribute {att} {
        variable attributes
        return [dict exists $attributes $att]
    }
    method getAttribute {att} {
        variable attributes
        set v ""
        if {[dict exists $attributes $att]} {
            set v [dict get $attributes $att]
        }
        return $v
    }
    method addAttribute {key values} {variable attributes; dict set attributes $key $values}
    method getAttributes {} {variable attributes; return $attributes}
    
    method getName { } {variable name; return $name}
    method getPublicName { } {variable publicname; return $publicname}
    method setPublicName { pn } {variable publicname; set publicname $pn}
    
    method getHelp { } {variable help; if {$help eq ""} {variable publicname; return $publicname} {return $help} }
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
        if {[dict exists $inputs $in]} {
            set pn [dict get $inputs $in]
        }
        return $pn
    }
    method getInputDv {in} {
        variable inputs
        set dv ""
        if {[dict exists $inputs $in]} {
            set i [dict get $inputs $in]
            set dv [$i getDv]
        }
        return $dv
    }
    
    method getOutputs { } {variable outputs; return $outputs}
    method getOutputPn {in} {
        variable outputs
        set pn ""
        if {[dict exists $outputs $in]} {
            set pn [dict get $outputs $in]
        }
        return $pn
    }
    

    method addInputDependency {n dn dv} {
        variable inputs
        [dict get $inputs $n] setDep $dn $dv

    }
    
    
    method cumple {args} {
        #W "Cumplimos con los filtros: $args"
        set c 1
        if {$::Kratos::kratos_private(DevMode) ne "release"} {
            
        } elseif {[my getAttribute "ProductionReady"] ne "" && [my getAttribute "ProductionReady"] ne "ProductionReady"} {
            #W "[my getName] no paso - [my getAttribute "ProductionReady"] "
            return 0
        }
        if {$args ne ""} {
            foreach {k listfiltervalues} {*}$args {
                set listfiltervalues [string map {, " "} $listfiltervalues]
                #W "k: $k vs : $listfiltervalues"
                set listattributesvalues [my getAttribute $k]
                #W "My value $listattributesvalues"
                set b1 0
                foreach possiblevalue $listfiltervalues {
                    if {$possiblevalue in $listattributesvalues} {
                        set b1 1
                        break
                    }
                }
               
                if {$b1 eq 0} {
                    set c 0
                    break
                }
            }
        }
        #if {$c eq 1 } {W "[my getName] Cumples"} {W "[my getName] No Cumples"}
        return $c
    }
}
}
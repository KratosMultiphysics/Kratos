namespace eval Pfem::xml {
    variable dir
}

proc Pfem::xml::Init { } {
    variable dir
    Model::InitVariables dir $Pfem::dir
    
    Model::getSolutionStrategies Strategies.xml
    Model::getElements Elements.xml
    Model::getConstitutiveLaws "../../Solid/xml/ConstitutiveLaws.xml"
    Model::getConstitutiveLaws "../../Pfem/xml/ConstitutiveLaws.xml"
    Model::getProcesses "../../Solid/xml/Processes.xml"
    Model::getNodalConditions "../../Solid/xml/NodalConditions.xml"
    Model::getConditions "../../Solid/xml/Conditions.xml"
    Model::getSolvers "../../Common/xml/Solvers.xml"
    
}

proc Pfem::xml::getUniqueName {name} {
    return PFEM_$name
}


proc ::Pfem::xml::MultiAppEvent {args} {
    if {$args eq "init"} {
        spdAux::parseRoutes
        spdAux::ConvertAllUniqueNames SL PFEM_
    }
}

proc ::Pfem::xml::::CheckElementOutputState { domNode args } {
    set elemsactive [list ]
    foreach parts_un [Pfem::write::GetPartsUN] {
        set parts_path [spdAux::getRoute $parts_un]
        set xp1 "$parts_path/group/value\[@n='Element'\]"
        foreach gNode [[customlib::GetBaseRoot] selectNodes $xp1] {
            lappend elemsactive [get_domnode_attribute $gNode v]
        }
    }
    set paramName [$domNode @n]
    return [::Model::CheckElementOutputState $elemsactive $paramName]
}

proc Pfem::xml::ProcGetElementsDict {domNode args} {
    set names [list ]
    set blockNode [Pfem::xml::FindMyBlocknode $domNode]
    set BodyType [get_domnode_attribute [$blockNode selectNodes "value\[@n='BodyType'\]"] v]
    #W $BodyType
    set argums [list ElementType $BodyType]
    set elems [Pfem::xml::GetElements $domNode $args]
    foreach elem $elems {
        #W [$elem getName]
        if {[$elem cumple $argums]} {
            lappend pnames [$elem getName] 
            lappend pnames [$elem getPublicName]
        }
    }
    set diction [join $pnames ","]
    return $diction
}
proc Pfem::xml::ProcGetElementsValues {domNode args} {
    set names [list ]
    set blockNode [Pfem::xml::FindMyBlocknode $domNode]
    set BodyType [get_domnode_attribute [$blockNode selectNodes "value\[@n='BodyType'\]"] v]

    set argums [list ElementType $BodyType]
    set elems [Pfem::xml::GetElements $domNode $args]
    foreach elem $elems {
        if {[$elem cumple $argums]} {
            lappend names [$elem getName]
        }
    }
    set values [join $names ","]
    
    if {[get_domnode_attribute $domNode v] eq ""} {$domNode setAttribute v [lindex $names 0]}
    if {[get_domnode_attribute $domNode v] ni $names} {$domNode setAttribute v [lindex $names 0]}

    return $values
}

proc Pfem::xml::GetElements {domNode args} {
    
    set nodeApp [spdAux::GetAppIdFromNode $domNode]
    set sol_stratUN [apps::getAppUniqueName $nodeApp SolStrat]
    set schemeUN [apps::getAppUniqueName $nodeApp Scheme]
    
    get_domnode_attribute [$domNode selectNodes [spdAux::getRoute $sol_stratUN]] dict
    get_domnode_attribute [$domNode selectNodes [spdAux::getRoute $schemeUN]] dict
    
    set solStratName [::write::getValue $sol_stratUN]
    set schemeName [write::getValue $schemeUN]
    set elems [::Model::GetAvailableElements $solStratName $schemeName]
    
    return $elems
}

proc Pfem::xml::FindMyBlocknode {domNode} {
    set top 10
    set ret ""
    for {set i 0} {$i < $top} {incr i} {
        if {[$domNode nodeName] eq "blockdata"} {
            set ret $domNode
            break
        } else {
            set domNode [$domNode parent]     
        }
    }
    return $ret
}

proc Pfem::xml::ProcGetMeshingDomains {domNode args} {
    set basepath [spdAux::getRoute "PFEM_meshing_domains"]
    set values [list ]
    foreach meshing_domain [[$domNode selectNodes $basepath] childNodes] {
        #W [$meshing_domain asXML]
        lappend values [get_domnode_attribute $meshing_domain name]
    }
    if {[get_domnode_attribute $domNode v] eq ""} {
        $domNode setAttribute v [lindex $values 0]
    }
    return [join $values ,]
}

proc Pfem::xml::ProcGetContactDomains {domNode args} {
    set basepath [spdAux::getRoute "PFEM_contact_domains"]
    set values [list ]
    set pvalues [list ]
    catch {
        foreach contact_domain [[$domNode selectNodes $basepath] childNodes] {
            #W [$contact_domain asXML]
            lappend values [get_domnode_attribute $contact_domain n]
            lappend pvalues [get_domnode_attribute $contact_domain n]
            lappend pvalues [get_domnode_attribute $contact_domain pn]
        }
    }
    $domNode setAttribute values [join $values ,]
    
    if {[get_domnode_attribute $domNode v] eq "" || [get_domnode_attribute $domNode v] ni $values} {
        $domNode setAttribute v [lindex $values 0]
    }
    return [join $pvalues ,]
}

proc Pfem::xml::ProcCheckNodalConditionStateSolid {domNode args} {
     # Overwritten the base function to add Solution Type restrictions
     set elemsactive [list ]
     foreach parts_un [Pfem::write::GetPartsUN] {
         set parts_path [spdAux::getRoute $parts_un]
         set xp1 "$parts_path/group/value\[@n='Element'\]"
         foreach gNode [[customlib::GetBaseRoot] selectNodes $xp1] {
             lappend elemsactive [get_domnode_attribute $gNode v]
         }
     }
     if {$elemsactive eq ""} {return "hidden"}
     set elemsactive [lsort -unique $elemsactive]
     set conditionId [$domNode @n]
     set solutionType [get_domnode_attribute [$domNode selectNodes [spdAux::getRoute PFEM_SolutionType]] v]
     set params [list analysis_type $solutionType]
     if {[::Model::CheckElementsNodalCondition $conditionId $elemsactive $params]} {return "normal"} else {return "hidden"}
}

proc Pfem::xml::ProcSolutionTypeState {domNode args} {
    set domain_type_un PFEM_DomainType
    set domain_type_route [spdAux::getRoute $domain_type_un]
    set state normal
    if {$domain_type_route ne ""} {
        set domain_type_node [$domNode selectNodes $domain_type_route]
        set domain_type_value [get_domnode_attribute $domain_type_node v]
        
        if {$domain_type_value eq "Fluids"} {
            $domNode setAttribute values Dynamic 
            $domNode setAttribute v Dynamic
            set state disabled
        } {
            $domNode setAttribute values "Dynamic,Static"
            set state normal
        }
    }
    return $state
}

proc Pfem::xml::ProcGetBodyTypeValues {domNode args} {
    set domain_type_un PFEM_DomainType
    set domain_type_route [spdAux::getRoute $domain_type_un]
    set values "Fluid,Solid,Rigid"
    if {$domain_type_route ne ""} {
        set domain_type_node [$domNode selectNodes $domain_type_route]
        set domain_type_value [get_domnode_attribute $domain_type_node v]
        
        if {$domain_type_value eq "Fluids"} {
            set values "Fluid,Rigid"
        }
        if {$domain_type_value eq "Solids"} {
            set values "Solid,Rigid"
        }
    }
    return $values
}

proc Pfem::xml::ProcGetSolutionStrategiesPFEM {domNode args} {
    set names ""
    set pnames ""
    set solutionType [get_domnode_attribute [$domNode selectNodes [spdAux::getRoute PFEM_SolutionType]] v]
    set Sols [::Model::GetSolutionStrategies [list "SolutionType" $solutionType] ]
    set ids [list ]
    set domainType [get_domnode_attribute [$domNode selectNodes [spdAux::getRoute PFEM_DomainType]] v]
    set filter [list Solid Pfem]
    if {$domainType eq "Solids"} {set filter "Solid"}
    if {$domainType eq "Fluids"} {set filter "Pfem"}
    
    foreach ss $Sols {
        if {[$ss getAttribute "App"] in $filter} {
            lappend names [$ss getName]
            lappend pnames [$ss getName]
            lappend pnames [$ss getPublicName]
        }
    }
    
    $domNode setAttribute values [join $names ","]
    set dv [lindex $names 0]
    #W "dv $dv"
    if {[$domNode getAttribute v] eq ""} {$domNode setAttribute v $dv; spdAux::RequestRefresh}
    if {[$domNode getAttribute v] ni $names} {$domNode setAttribute v $dv; spdAux::RequestRefresh}

    return [join $pnames ","]
}

proc Pfem::xml::ProcGetPartUN {domNode args} {
    customlib::UpdateDocument
    set root [customlib::GetBaseRoot]
    set xp1 "[spdAux::getRoute "PFEM_Bodies"]/blockdata/condition"
    set i 0
    foreach part_node [$root selectNodes $xp1] {
        if {$part_node eq $domNode} {
              break
        } {incr i}
    }
    set un "PFEM_Part$i"
    #$domNode setAttribute curr_un $un
    return $un
}

Pfem::xml::Init

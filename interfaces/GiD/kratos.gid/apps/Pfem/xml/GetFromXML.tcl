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
    Model::getSolvers "../../Pfem/xml/Solvers.xml"
    
}

proc Pfem::xml::getUniqueName {name} {
    return PFEM_$name
}

proc Pfem::xml::MultiAppEvent {args} {
    if {$args eq "init"} {
        spdAux::parseRoutes
        spdAux::ConvertAllUniqueNames SL PFEM_
    }
}

proc Pfem::xml::CustomTree { args } {
    # Hide Results Cut planes
    spdAux::SetValueOnTreeItem v time Results FileLabel
    spdAux::SetValueOnTreeItem v time Results OutputControlType
    spdAux::SetValueOnTreeItem v 0.04 Results OutputDeltaTime
    
    for {set i 0} {$i < 3} {incr i} {
        GiD_Groups create "Group $i"
        spdAux::AddConditionGroupOnXPath "[spdAux::getRoute "PFEM_NodalConditions"]/condition\[@n='PRESSURE'\]" "Group $i"
    }
}

proc Pfem::xml::CheckElementOutputState { domNode args } {
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
    set argums [list ElementType $BodyType]
    set elems [Pfem::xml::GetElements $domNode $args]
    set pnames ""
    foreach elem $elems {
        if {[$elem cumple $argums]} {
            lappend pnames [$elem getName] 
            lappend pnames [$elem getPublicName]
        }
    }
    set diction [join $pnames ","]
    if {$diction eq ""} {W "No available elements - Check Solution strategy & scheme - Check Kratos mode (developer)"}
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
    spdAux::setRoute $un [$part_node toXPath]
    #$domNode setAttribute curr_un $un
    return $un
}

proc Pfem::xml::ProcPartsOverWhat {domNode args} {
    set names [list ]
    set blockNode [Pfem::xml::FindMyBlocknode $domNode]
    set BodyType [get_domnode_attribute [$blockNode selectNodes "value\[@n='BodyType'\]"] v]
    if {$BodyType eq "Fluid" || $BodyType eq "Solid"} {
        if {$::Model::SpatialDimension eq "3D"} {
            return "volume"
        } else {
            return "surface"
        }
    } elseif { $BodyType eq "Rigid"} {
        if {$::Model::SpatialDimension eq "3D"} {
            return "surface,volume"
        } else {
            return "line,surface"
        }
    } else {
        return "point,line,surface,volume"
    }
}

proc Pfem::xml::ProcActiveIfAnyPartState {domNode args} {
    set parts ""
    set parts_un [Pfem::xml::ProcGetPartUN $domNode $args]
    catch {
        set parts [$domNode selectNodes "[spdAux::getRoute $parts_un]/group"]
    }
    if {$parts ne ""} {return "normal"} else {return "hidden"}
}

proc Pfem::xml::ProcGetBodiesValues {$domNode $args} {
    customlib::UpdateDocument
    set root [customlib::GetBaseRoot]
    set xp1 "[spdAux::getRoute "PFEM_Bodies"]/blockdata"
    set bodies [list ]
    foreach body_node [$root selectNodes $xp1] {
        lappend bodies [$body_node @name]
    }
    return [join $bodies ","]
}

proc Pfem::xml::StartSortingWindow { } {
    set data_dict [dict create]
    set conds [Pfem::xml::GetConditionsAndGroups PFEM_Loads]
    set nodalconds [Pfem::xml::GetConditionsAndGroups PFEM_NodalConditions]
    if {[dict size $conds]} {dict set data_dict Loads $conds}
    if {[dict size $nodalconds]} {dict set data_dict Constraints $nodalconds}
    SorterWindow::SorterWindow $data_dict "Pfem::xml::GetDataFromSortingWindow"
}
proc Pfem::xml::GetDataFromSortingWindow { data_dict } {
    W $data_dict
}
proc Pfem::xml::GetConditionsAndGroups { cnd_UN } {
    customlib::UpdateDocument
    set data_dict [dict create]
    set root [customlib::GetBaseRoot]
    foreach {cond_type cond_item cond_item_name} {container blockdata name condition group n} {
        set xp1 "[spdAux::getRoute $cnd_UN]/$cond_type"
        foreach cnd_cont_node [$root selectNodes $xp1] {
            set cnd_cont_name [$cnd_cont_node @n]
            set xp2 "./$cond_item"
            foreach cnd_node [$cnd_cont_node selectNodes $xp2] {
                set cnd_name [$cnd_node getAttribute $cond_item_name]
                set num 0
                if {[$cnd_node hasAttribute order]} {set num [$cnd_node @order]}
                dict set data_dict $cnd_cont_name $cnd_name $num
            }
        }
    }
    return $data_dict
}

Pfem::xml::Init

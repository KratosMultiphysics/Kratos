namespace eval Pfem::xml {
    variable dir
    variable bodyNodalCondition
}

proc Pfem::xml::Init { } {
    variable dir
    variable bodyNodalCondition
    
    set bodyNodalCondition [list ]

    Model::InitVariables dir $Pfem::dir
    
    Model::getSolutionStrategies Strategies.xml
    Model::getElements Elements.xml
    Model::getConstitutiveLaws "../../Solid/xml/ConstitutiveLaws.xml"
    Model::getConstitutiveLaws "../../Pfem/xml/ConstitutiveLaws.xml"
    Model::getProcesses "../../Solid/xml/Processes.xml"
    Model::getNodalConditions "../../Solid/xml/NodalConditions.xml"
    Model::getConditions "../../Solid/xml/Conditions.xml"
    Model::getSolvers "../../Pfem/xml/Solvers.xml"
    Pfem::xml::getBodyNodalConditions BodyNodalConditions.xml
    
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
    
    #HOW TO USE THIS FUNCTION:
    #spdAux::SetValueOnTreeItem arg1 arg2 arg3 (arg4)
    #arg1: attribute_to_modify 
    #arg2: value_of_the_attribute 
    #arg3: unique_name_of_the_node  ('unique name is defined by the attribute un=)
    #arg4 (optional): name_of_the_child_we_want_to_modify  ('name'is defined by the attribute n=)
    
    # Hide Results Cut planes  
    
    #intervals
    spdAux::SetValueOnTreeItem icon timeIntervals Intervals
    foreach node [[customlib::GetBaseRoot] selectNodes "[spdAux::getRoute Intervals]/blockdata"] {$node setAttribute icon select}        
    foreach node [[customlib::GetBaseRoot] selectNodes "[spdAux::getRoute Intervals]/blockdata/value"] {$node setAttribute icon data}          
    #results
    spdAux::SetValueOnTreeItem v time Results FileLabel
    spdAux::SetValueOnTreeItem icon results Results
    spdAux::SetValueOnTreeItem icon seeResults Results 
    spdAux::SetValueOnTreeItem icon data Results FileLabel 
    spdAux::SetValueOnTreeItem icon data Results OutputControlType 
    spdAux::SetValueOnTreeItem icon data Results OutputDeltaTime 
    spdAux::SetValueOnTreeItem icon data Results BodyOutput 
    spdAux::SetValueOnTreeItem icon data Results NodeOutput 
    spdAux::SetValueOnTreeItem icon data Results SkinOutput 
    spdAux::SetValueOnTreeItem icon data Results OnElement 
    spdAux::SetValueOnTreeItem icon select Results OnNodes 
    spdAux::SetValueOnTreeItem icon select Results GiDOptions 
    spdAux::SetValueOnTreeItem v time Results OutputControlType
    spdAux::SetValueOnTreeItem v 0.04 Results OutputDeltaTime
    
    #problem settings
    spdAux::SetValueOnTreeItem icon folder PFEM_Implicitlinear_solver_settings 
    spdAux::SetValueOnTreeItem icon folder PFEM_Explicitlinear_solver_settings
    spdAux::SetValueOnTreeItem icon folder PFEM_GenericSolStratlinear_solver_settings
    spdAux::SetValueOnTreeItem icon folder PFEM_TwoStepVPStrategyvelocity_linear_solver_settings
    spdAux::SetValueOnTreeItem icon folder PFEM_TwoStepVPStrategypressure_linear_solver_settings
    
    #restart
    spdAux::SetValueOnTreeItem icon doRestart Restart 
    spdAux::SetValueOnTreeItem icon data Restart SaveRestart
    spdAux::SetValueOnTreeItem icon data Restart LoadRestart
    
    
    #parallelism
    spdAux::SetValueOnTreeItem icon select Parallelization
    spdAux::SetValueOnTreeItem values OpenMP ParallelType 
    spdAux::SetValueOnTreeItem icon data ParallelType
    spdAux::SetValueOnTreeItem icon data Parallelization OpenMPNumberOfThreads
    
    #boundary conditions
    spdAux::SetValueOnTreeItem icon folder PFEM_NodalConditions DISPLACEMENT
    spdAux::SetValueOnTreeItem icon folder PFEM_NodalConditions VELOCITY
    spdAux::SetValueOnTreeItem icon folder PFEM_NodalConditions ACCELERATION
    spdAux::SetValueOnTreeItem icon folder PFEM_NodalConditions PRESSURE
    [[customlib::GetBaseRoot] selectNodes "[spdAux::getRoute PFEM_NodalConditions]/container\[@n='BODYDISPLACEMENT'\]"] setAttribute icon folder
    foreach node [[customlib::GetBaseRoot] selectNodes "[spdAux::getRoute PFEM_NodalConditions]/container\[@n='BODYDISPLACEMENT'\]/blockdata"] {$node setAttribute icon select}
    
    #loads
    spdAux::SetValueOnTreeItem icon setLoad PFEM_Loads 
    spdAux::SetValueOnTreeItem icon folder PFEM_Loads SelfWeight3D
    spdAux::SetValueOnTreeItem icon folder PFEM_Loads SelfWeight2D
    spdAux::SetValueOnTreeItem icon folder PFEM_Loads SelfWeight2Da
    spdAux::SetValueOnTreeItem icon folder PFEM_Loads PointLoad2D
    spdAux::SetValueOnTreeItem icon folder PFEM_Loads PointLoad2DAxisym
    spdAux::SetValueOnTreeItem icon folder PFEM_Loads PointLoad3D
    spdAux::SetValueOnTreeItem icon folder PFEM_Loads LineLoad2D
    spdAux::SetValueOnTreeItem icon folder PFEM_Loads LineLoad2DAxisym
    spdAux::SetValueOnTreeItem icon folder PFEM_Loads SurfaceLoad3D
    spdAux::SetValueOnTreeItem icon folder PFEM_Loads LinePressure2D
    spdAux::SetValueOnTreeItem icon folder PFEM_Loads LinePressure2DAxisym
    spdAux::SetValueOnTreeItem icon folder PFEM_Loads SurfacePressure3D
   
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
        lappend values [get_domnode_attribute $meshing_domain name]
    }
    if {[get_domnode_attribute $domNode v] eq ""} {
	$domNode setAttribute v [lindex $values 0]
    }
    return [join $values ,]
}

proc Pfem::xml::ProcGetContactDomains {domNode args} {
    set basepath [spdAux::getRoute "PFEM_contacts"]
    set values [list "No contact strategy"]
    foreach contact_domain [[$domNode selectNodes $basepath] childNodes] {
        lappend values [get_domnode_attribute $contact_domain name]
    }
        
    if {[get_domnode_attribute $domNode v] eq "" || [get_domnode_attribute $domNode v] ni $values} {
        $domNode setAttribute v [lindex $values 0]
    }
    return [join $values ,]
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
	
	if {$domain_type_value ne "Solids"} {
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
	if {$domain_type_value eq "Coupled"} {
	    set values "Solid,Fluid,Rigid"
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
    if {$domainType eq "Coupled"} {set filter "Pfem"}
    
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

proc Pfem::xml::ProcGetBodiesValues {domNode args} {
    customlib::UpdateDocument
    set root [customlib::GetBaseRoot]
    set xp1 "[spdAux::getRoute "PFEM_Bodies"]/blockdata"
    set bodies [list ]
    foreach body_node [$root selectNodes $xp1] {
        lappend bodies [$body_node @name]
    }
    if {[get_domnode_attribute $domNode v] ni $bodies} {$domNode setAttribute v [lindex $bodies 0]}
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

proc Pfem::xml::getBodyNodalConditionById { id } {
    variable bodyNodalCondition
    
    foreach cnd $bodyNodalCondition {
        if {[$cnd getName] eq $id} {
            return $cnd
        }
    }
    return ""
}
proc Pfem::xml::getBodyNodalConditions { filename } {
    variable bodyNodalCondition
    dom parse [tDOM::xmlReadFile [file join $Pfem::dir xml $filename]] doc
    
    set NCList [$doc getElementsByTagName NodalConditionItem]
    foreach Node $NCList {
        lappend bodyNodalCondition [::Model::ParseNodalConditionsNode $Node]
    }
}
proc Pfem::xml::injectBodyNodalConditions { basenode args} {
    variable bodyNodalCondition
    Pfem::xml::_injectCondsToTree $basenode $bodyNodalCondition nodal
    $basenode delete
}


proc Pfem::xml::_injectCondsToTree {basenode cond_list {cond_type "normal"} } {
    set conds [$basenode parent]
    set AppUsesIntervals [::Pfem::GetAttribute UseIntervals]
    if {$AppUsesIntervals eq ""} {set AppUsesIntervals 0}
    
    foreach cnd $cond_list {
        set n [$cnd getName]
        set pn [$cnd getPublicName]
        set help [$cnd getHelp]
        set units [$cnd getAttribute "units"]
        set um [$cnd getAttribute "unit_magnitude"]
        set process [::Model::GetProcess [$cnd getProcessName]]
        set check [$process getAttribute "check"]
        if {$check eq ""} {set check "UpdateTree"}
        set state "ConditionState"
        if {$cond_type eq "nodal"} {
            set state [$cnd getAttribute state]
            if {$state eq ""} {set state "CheckNodalConditionState"}
        }
        set contNode [gid_groups_conds::addF [$conds toXPath] container [list n $n pn ${pn}s help $help]]
        set blockNode [gid_groups_conds::addF [$contNode toXPath] blockdata [list n $n pn $pn help $help icon shells16 update_proc $check name "$pn 1" sequence 1 editable_name unique sequence_type non_void_disabled]]
        set block_path [$blockNode toXPath]
        set inputs [$process getInputs] 
        foreach {inName in} $inputs {
            set pn [$in getPublicName]
            set type [$in getType]
            set v [$in getDv]
            set help [$in getHelp]
            set state [$in getAttribute "state"]
            if {$state eq ""} {set state "normal"}
            foreach key [$cnd getDefaults $inName] {
                set $key [$cnd getDefault $inName $key]
            }
            
            set has_units [$in getAttribute "has_units"]
            if {$has_units ne ""} { set has_units "units='$units'  unit_magnitude='$um'"}
            if {$type eq "vector"} {
                set vector_type [$in getAttribute "vectorType"]
                lassign [split $v ","] v1 v2 v3
                if {$vector_type eq "bool"} {
                    gid_groups_conds::addF $block_path value [list n ${inName}X wn [concat $n "_X"] pn "X ${pn}" values "1,0"]
                    gid_groups_conds::addF $block_path value [list n ${inName}Y wn [concat $n "_Y"] pn "Y ${pn}" values "1,0"]
                    gid_groups_conds::addF $block_path value [list n ${inName}Z wn [concat $n "_Z"] pn "Z ${pn}" values "1,0" state {[CheckDimension 3D]}]
                } {
                    foreach i [list "X" "Y" "Z"] {
                        set nodev "../value\[@n='${inName}$i'\]"
                        set zstate ""
                        if {$i eq "Z"} { set zstate "state {\[CheckDimension 3D\]}"}
                        if {[$in getAttribute "enabled"] in [list "1" "0"]} {
                            set val [expr [$in getAttribute "enabled"] ? "Yes" : "No"]
                            if {$i eq "Z"} { set val "No" }
                            set valNode [gid_groups_conds::addF $block_path value [list n Enabled_$i pn "$i component" v No values "Yes,No" help "Enables the $i ${inName}" actualize_tree 1 {*}$zstate]]

                            gid_groups_conds::addF [$valNode toXPath] dependencies [list value No node $nodev att1 state v1 hidden]
                            gid_groups_conds::addF [$valNode toXPath] dependencies [list value Yes node $nodev att1 state v1 normal]
                            if {[$in getAttribute "function"] eq "1"} {
                                set fname "${i}function_$inName"
                                set nodef "../value\[@n='$fname'\]"
                                set nodeb "../value\[@n='ByFunction$i'\]"
                                gid_groups_conds::addF [$valNode toXPath] dependencies [list value No node $nodef att1 state v1 hidden]
                                gid_groups_conds::addF [$valNode toXPath] dependencies [list value No node $nodeb att1 state v1 hidden]
                                gid_groups_conds::addF [$valNode toXPath] dependencies [list value Yes node $nodeb att1 state v1 normal att2 v v2 No]
                            }
                        }
                        if {[$in getAttribute "function"] eq "1"} {
                            set fname "${i}function_$inName"
                            set valNode [gid_groups_conds::addF $block_path value [list n ByFunction$i pn "by function -> f(x,y,z,t)" v No values "Yes,No" actualize_tree 1 state hidden]]
                            gid_groups_conds::addF [$valNode toXPath] dependencies [list value No node $nodev att1 state v1 normal]
                            gid_groups_conds::addF [$valNode toXPath] dependencies [list value Yes node $nodev att1 state v1 hidden]
                            gid_groups_conds::addF [$valNode toXPath] dependencies [list value No node $nodef att1 state v1 hidden]
                            gid_groups_conds::addF [$valNode toXPath] dependencies [list value Yes node $nodef att1 state v1 normal]
                            gid_groups_conds::addF $block_path value [list n $fname pn "$i function" state hidden]
                        }
                        gid_groups_conds::addF $block_path value [list n ${inName}$i wn [concat $n "_$i"] pn "$i ${pn}" v $v1 state hidden]
                    }
                }
                
            } elseif { $type eq "combo" } {
                set values [join [$in getValues] ","]
                gid_groups_conds::addF $block_path value [list n $inName pn $pn v $v1 values $values state $state help $help]
            } elseif { $type eq "bool" } {
                set values "1,0"
                gid_groups_conds::addF $block_path value [list n $inName pn $pn v $v1 values $values state $state help $help]
            } elseif { $type eq "file" || $type eq "tablefile" } {
                gid_groups_conds::addF $block_path value [list n $inName pn $pn v $v1 values {[GetFilesValues]} update_proc AddFile type $type state $state help $help]
            } else {
                if {[$in getAttribute "function"] eq "1"} {
                    set fname "function_$inName"
                    set nodev "../value\[@n='$inName'\]"
                    set nodef "../value\[@n='$fname'\]"
                    
                    set valNode [gid_groups_conds::addF $block_path value [list n ByFunction pn "by function -> f(x,y,z,t)" v No values "Yes,No" actualize_tree 1]]
                    gid_groups_conds::addF [$valNode toXPath] dependencies [list value No node $nodev att1 state v1 normal]
                    gid_groups_conds::addF [$valNode toXPath] dependencies [list value Yes node $nodev att1 state v1 hidden]
                    gid_groups_conds::addF [$valNode toXPath] dependencies [list value No node $nodef att1 state v1 hidden]
                    gid_groups_conds::addF [$valNode toXPath] dependencies [list value Yes node $nodef att1 state v1 normal]
                    gid_groups_conds::addF $block_path value [list n $fname pn "Function"]
                }
                append node "<value n='$inName' pn='$pn' v='$v'  units='$units'  unit_magnitude='$um'  help='$help'/>"
                gid_groups_conds::addF $block_path value [list n $inName pn $pn v $v units $units unit_magnitude $um help $help]
            }
        }
        
        set CondUsesIntervals [$cnd getAttribute "Interval"]
        if {$AppUsesIntervals && $CondUsesIntervals ne "False"} {
            gid_groups_conds::addF $block_path value [list n Interval pn "Time interval" v $CondUsesIntervals values {[getIntervals]} help $help]
        }
        gid_groups_conds::addF $block_path value [list n Body pn Body v - values {[GetBodiesValues]} help $help]
    }
    
}

Pfem::xml::Init

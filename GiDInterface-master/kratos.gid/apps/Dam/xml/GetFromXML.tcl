namespace eval Dam::xml {
    variable dir
}

proc Dam::xml::Init { } {
    variable dir
    Model::InitVariables dir $Dam::dir
    
    Model::getSolutionStrategies Strategies.xml
    Model::getElements Elements.xml
    Model::getMaterials Materials.xml
    Model::getNodalConditions NodalConditions.xml
    Model::getConstitutiveLaws ConstitutiveLaws.xml
    Model::getProcesses Processes.xml
    Model::getConditions Conditions.xml
    Model::getSolvers "../../Common/xml/Solvers.xml"
}

proc Dam::xml::getUniqueName {name} {
    return Dam$name
}


proc ::Dam::xml::MultiAppEvent {args} {
    if {$args eq "init"} {
        spdAux::parseRoutes
        spdAux::ConvertAllUniqueNames SL Dam
    }
}


proc Dam::xml::CustomTree { args } {
	
    # Add some nodal results
    set nodal_results_base [[customlib::GetBaseRoot] selectNodes [spdAux::getRoute NodalResults]]
    $nodal_results_base setAttribute state "\[ActiveIfAnyPartState\]"
    if {$nodal_results_base ne ""} {
        set delete_list [list "INITIALTEMPERATURE" "BOFANGTEMPERATURE" "CONSTANTRESERVOIRTEMPERATURE"]
        foreach item $delete_list {
            set i [$nodal_results_base selectNodes "./value\[@n='$item'\]"]
            if {$i ne ""} {$i delete}
        }
        
        ## It has special filter
        set add_special_list [list TEMPERATURE HEAT_FLUX FACE_HEAT_FLUX]
        set add_special_list_pn [list "Temperature" "Heat Source" "Heat Fluxes"]
        foreach it $add_special_list pn $add_special_list_pn {
               gid_groups_conds::addF [$nodal_results_base toXPath] value [list n $it pn $pn v "No" values "Yes,No" state "\[checkStateByUniqueName DamTypeofProblem UP_Thermo-Mechanical DamTypeofProblem Thermo-Mechanical\]"]
        }
        
        set add_list [list ADDED_MASS VOLUME_ACCELERATION POINT_LOAD LINE_LOAD SURFACE_LOAD POSITIVE_FACE_PRESSURE NODAL_CAUCHY_STRESS_TENSOR NODAL_JOINT_WIDTH Vi_POSITIVE Viii_POSITIVE NODAL_YOUNG_MODULUS]
        set add_list_pn [list "Added Mass" "Body Accelerations" "Point Loads" "Line Loads" "Surface Loads" "Normal Loads" "Nodal Total Stress" "Nodal Joint Width" "Traction Principal Stress Vector" "Compression Principal Stress Vector" "Nodal Young Modulus"]
        foreach item $add_list pn $add_list_pn {
               gid_groups_conds::addF [$nodal_results_base toXPath] value [list n $item pn $pn v "No" values "Yes,No" state "\[checkStateByUniqueName DamTypeofProblem UP_Mechanical DamTypeofProblem UP_Thermo-Mechanical DamTypeofProblem Mechanical DamTypeofProblem Thermo-Mechanical\]"]
        }
    }
       
}

proc Dam::xml::ProcGetSchemes {domNode args} {
    set type_of_problem [write::getValue DamTypeofProblem]
    
    set sol_stratUN "DamSolStratTherm"
    if {$type_of_problem ne "Thermal"} {
        set sol_stratUN "DamSolStrat"
    }
    
    set solStratName [write::getValue $sol_stratUN]
    set schemes [::Model::GetAvailableSchemes $solStratName]
    
    set ids [list ]
    if {[llength $schemes] == 0} {
        if {[get_domnode_attribute $domNode v] eq ""} {$domNode setAttribute v "None"}
        return "None"
    }
    set names [list ]
    set pnames [list ]
    foreach cl $schemes {
        lappend names [$cl getName]
        lappend pnames [$cl getName] 
        lappend pnames [$cl getPublicName]
    }
    if {$type_of_problem in [list "Acoustic" "UP_Mechanical" "UP_Thermo-Mechanical"]} {
        set names [list "Newmark"]
        set pnames [list "Newmark" "Newmark"]
    }
    
    if {$type_of_problem in [list "Modal-Analysis"]} {
        set names [list "Eigen-Dynamic-Scheme"]
        set pnames [list "Eigen-Dynamic-Scheme" "Eigen Dynamic Scheme"]
    }
    
    $domNode setAttribute values [join $names ","]
    
    if {[get_domnode_attribute $domNode v] eq ""} {$domNode setAttribute v [lindex $names 0]}
    if {[get_domnode_attribute $domNode v] ni $names} {$domNode setAttribute v [lindex $names 0]}
    spdAux::RequestRefresh
    return [join $pnames ","]
}


proc Dam::xml::SolStratParamState {outnode} {    
    set root [customlib::GetBaseRoot]
    
    set solstrat_un "DamSolStrat"
    
    #W $solstrat_un
    if {[get_domnode_attribute [$root selectNodes [spdAux::getRoute $solstrat_un]] v] eq ""} {
        get_domnode_attribute [$root selectNodes [spdAux::getRoute $solstrat_un]] dict
    }
    set SolStrat [get_domnode_attribute [$root selectNodes [spdAux::getRoute $solstrat_un]] v]
    set paramName [$outnode @n]
    set ret [::Model::CheckSolStratInputState $SolStrat $paramName]
    if {$ret} {
        lassign [Model::GetSolStratParamDep $SolStrat $paramName] depN depV
        foreach node [[$outnode parent] childNodes] {
            if {[$node @n] eq $depN} {
                if {[get_domnode_attribute $node v] ni [split $depV ,]} {
                    set ret 0
                    break
                }
            }
        }
    }
    return $ret
}

proc Dam::xml::ProcGetConstitutiveLaws {domNode args} {
    set Elementname [$domNode selectNodes {string(../value[@n='Element']/@v)}]
    set Claws [::Model::GetAvailableConstitutiveLaws $Elementname]
    #W "Round 1 $Claws"
    #foreach cl $Claws {W [$cl getName]}
    set type_of_problem [write::getValue DamTypeofProblem]
    #W $type_of_problem
    set goodList [list ]
    #W "Pre type problem -> $type_of_problem"
    foreach cl $Claws {
        set type [$cl getAttribute Type]
        #W $type
        #W "cl -> [$cl getName]"
        #W "type -> $type"
        if {[string first "Therm" $type] eq -1 && $type_of_problem ni [list "Thermo-Mechanical" "UP_Thermo-Mechanical"]} {
            lappend goodList $cl
        } elseif {[string first "Therm" $type] ne -1 && $type_of_problem in [list "Thermo-Mechanical" "UP_Thermo-Mechanical"]} {
            lappend goodList $cl
        } elseif {[string first "Interface" $type] ne -1} {
            lappend goodList $cl
        } elseif {[string first "WaveEquationElement" $Elementname] ne -1 && $type eq "Wave"} {
            lappend goodList $cl
        } 
    }
    #W "good $goodList"
    set Claws $goodList
    #W "Round 2 $Claws"
    #foreach cl $Claws {W [$cl getName]}
    #
    set analysis_type ""
    set damage_type ""
    set TypeofProblem [get_domnode_attribute [$domNode selectNodes [spdAux::getRoute DamTypeofProblem]] v]
    switch $TypeofProblem {
        "Mechanical" {
            set analysis_type [get_domnode_attribute [$domNode selectNodes "[spdAux::getRoute DamMechanicalData]/value\[@n='AnalysisType'\]"] v]
            set damage_type [get_domnode_attribute [$domNode selectNodes "[spdAux::getRoute DamMechanicalData]/value\[@n='DamageType'\]"] v]
        }
        "Thermo-Mechanical" {
            set analysis_type [get_domnode_attribute [$domNode selectNodes "[spdAux::getRoute "DamThermo-MechanicalData"]/container\[@n='MechanicalPartProblem'\]/value\[@n='AnalysisType'\]"] v]
            set damage_type [get_domnode_attribute [$domNode selectNodes "[spdAux::getRoute "DamThermo-MechanicalData"]/container\[@n='MechanicalPartProblem'\]/value\[@n='DamageType'\]"] v]
        }
        "UP_Mechanical" {
            set analysis_type [get_domnode_attribute [$domNode selectNodes "[spdAux::getRoute "DamUP_MechanicalData"]/value\[@n='AnalysisType'\]"] v]
            set damage_type [get_domnode_attribute [$domNode selectNodes "[spdAux::getRoute "DamUP_MechanicalData"]/value\[@n='DamageType'\]"] v]
        }
        "UP_Thermo-Mechanical" {
            set analysis_type [get_domnode_attribute [$domNode selectNodes "[spdAux::getRoute "DamUP_Thermo-MechanicalData"]/container\[@n='MechanicalPartProblem'\]/value\[@n='AnalysisType'\]"] v]
            set damage_type [get_domnode_attribute [$domNode selectNodes "[spdAux::getRoute "DamUP_Thermo-MechanicalData"]/container\[@n='MechanicalPartProblem'\]/value\[@n='DamageType'\]"] v]
        }
        "Acoustic" {
            set analysis_type ""
            set damage_type ""
        }
        "Modal-Analysis" {
            set analysis_type [get_domnode_attribute [$domNode selectNodes "[spdAux::getRoute DamUModalData]/container\[@n='StratParams'\]/value\[@n='AnalysisType'\]"] v]
            set damage_type ""
        }
        default {
            error [= "Check type of problem"]
        }
    }
    set goodList [list ]
    #W "Pre analysis -> $analysis_type"
    #W $Claws
    foreach cl $Claws {
        if {$analysis_type eq ""} {
            lappend goodList $cl
        } else {
            set type [split [$cl getAttribute AnalysisType] ","]
            #W "cl -> [$cl getName]"
            #W "type -> $type"
            if {"Non-Linear" in $type && $analysis_type eq "Non-Linear"} {
                set cl_damage_type [split [$cl getAttribute DamageType] ","]
                if {$damage_type in $cl_damage_type} {
                    lappend goodList $cl
                }
            }
            if {"Linear" in $type && $analysis_type eq "Linear"} {
                lappend goodList $cl
            }
        }
    }
    set Claws $goodList
    #W "Round 3 $Claws"
    #foreach cl $Claws {W [$cl getName]}
    
    #W "Const Laws que han pasado la criba: $Claws"
    if {[llength $Claws] == 0} {
        if {[get_domnode_attribute $domNode v] eq ""} {$domNode setAttribute v "None"}
        return "None"
    }
    set names [list ]
    foreach cl $Claws {
        lappend names [$cl getName]
    }
    set values [join $names ","]
    if {[get_domnode_attribute $domNode v] eq "" || [get_domnode_attribute $domNode v] ni $names} {$domNode setAttribute v [lindex $names 0]; spdAux::RequestRefresh}
    #W $values
    
    return $values
}

proc Dam::xml::ProcCheckNodalConditionState {domNode args} {
    set parts_un DamParts
    if {[spdAux::getRoute $parts_un] ne ""} {
        set conditionId [$domNode @n]
        set elems [$domNode selectNodes "[spdAux::getRoute $parts_un]/group/value\[@n='Element'\]"]
        set elemnames [list ]
        foreach elem $elems { set a [get_domnode_attribute $elem values]; set a [get_domnode_attribute $elem dict]; lappend elemnames [$elem @v]}
        set elemnames [lsort -unique $elemnames]
        
        # Mirar Type of problem y acceder al contenedor correcto
        set TypeofProblem [get_domnode_attribute [$domNode selectNodes [spdAux::getRoute DamTypeofProblem]] v]
        switch $TypeofProblem {
            "Mechanical" {
                set solutionType [get_domnode_attribute [$domNode selectNodes "[spdAux::getRoute DamMechanicalData]/value\[@n='SolutionType'\]"] v]
                set params [list analysis_type $solutionType TypeofProblem $TypeofProblem]
            }
            "Thermo-Mechanical" {
                set solutionType [get_domnode_attribute [$domNode selectNodes "[spdAux::getRoute "DamThermo-MechanicalData"]/container\[@n='MechanicalPartProblem'\]/value\[@n='SolutionType'\]"] v]
                set params [list analysis_type $solutionType TypeofProblem $TypeofProblem]                    
            }
            "UP_Mechanical" {
                set solutionType [get_domnode_attribute [$domNode selectNodes "[spdAux::getRoute "DamUP_MechanicalData"]/value\[@n='SolutionType'\]"] v]
                set params [list analysis_type $solutionType TypeofProblem $TypeofProblem]
            }
            "UP_Thermo-Mechanical" {
                set solutionType [get_domnode_attribute [$domNode selectNodes "[spdAux::getRoute "DamUP_Thermo-MechanicalData"]/container\[@n='MechanicalPartProblem'\]/value\[@n='SolutionType'\]"] v]
                set params [list analysis_type $solutionType TypeofProblem $TypeofProblem]                    
            }
            "Acoustic" {
                set params [list TypeofProblem $TypeofProblem]
            }
            "Modal-Analysis" {
                set params [list TypeofProblem $TypeofProblem]
            }
            default {
                error [= "Check type of problem"]
            }
        }
        if {[::Model::CheckElementsNodalCondition $conditionId $elemnames $params]} {return "normal"} else {return "hidden"}
    } {return "hidden"}
}

proc Dam::xml::ProcCheckConditionState {domNode args} {
    set cond_id [$domNode @n]
    
    # By active parts
    if {[spdAux::ProcActiveIfAnyPartState $domNode $args] eq "hidden"} {return "hidden"}
    
    # By dimension
    set resp [::Model::CheckConditionState $domNode]
    
    # By type of problem
    if {$resp} {
        set TypeofProblem [get_domnode_attribute [$domNode selectNodes [spdAux::getRoute DamTypeofProblem]] v]
        if {$TypeofProblem ni [dict get [[Model::getCondition $cond_id] getAttributes ] TypeofProblem]} {set resp 0}
    }
    
    # By active elements
    if {$resp} {
        set elems [$domNode selectNodes "[spdAux::getRoute DamParts]/group/value\[@n='Element'\]"]
        set elemnames [list ]
        foreach elem $elems { set a [get_domnode_attribute $elem values]; set a [get_domnode_attribute $elem dict]; lappend elemnames [$elem @v]}
        set elemnames [lsort -unique $elemnames]
        
        if {$cond_id in [list "UPCondition2D" "UPCondition3D" "FreeSurface2D" "FreeSurface3D" "InfiniteDomain2D" "InfiniteDomain3D"]} {
            if {"WaveEquationElement$::Model::SpatialDimension" ni $elemnames} {set resp 0}
        } else {
            if {"WaveEquationElement$::Model::SpatialDimension" eq $elemnames} {set resp 0}
        }
        set TypeofProblem [get_domnode_attribute [$domNode selectNodes [spdAux::getRoute DamTypeofProblem]] v]
        if {$TypeofProblem ni [dict get [[Model::getCondition [$domNode @n]] getAttributes ] TypeofProblem]} {set resp 0}
    }
    
    
    if {$resp} {return "normal"} else {return "hidden"}
}


proc Dam::xml::ProcGetSolutionStrategies {domNode args} {
    set names ""
    set pnames ""
    set ids [list ]
    set type_of_strategy [lindex $args 0]
    set type_of_problem [write::getValue DamTypeofProblem]
    if {$type_of_strategy eq "Mechanic"} {
        if {$type_of_problem eq "Acoustic"} {
            set n "Newton-Raphson"
        } elseif {$type_of_problem eq "Modal-Analysis"} {
            set n "Eigen-Strategy"
        } else {       
            set n [list "Newton-Raphson" "Arc-length"]
            set type_of_problem [get_domnode_attribute [$domNode selectNodes [spdAux::getRoute DamTypeofProblem]] v]
        }
    } else {
        # Thermal
        set n "Newton-Raphson"
    }
    
    
    set Sols [::Model::GetSolutionStrategies [list n $n] ]
    foreach ss $Sols {
        lappend names [$ss getName]
        lappend pnames [$ss getName]
        lappend pnames [$ss getPublicName] 
    }
    
    $domNode setAttribute values [join $names ","]
    set dv [lindex $names 0]
    #W "dv $dv"
    if {[$domNode getAttribute v] eq ""} {$domNode setAttribute v $dv; spdAux::RequestRefresh}
    if {[$domNode getAttribute v] ni $names} {$domNode setAttribute v $dv; spdAux::RequestRefresh}
    
    return [join $pnames ","]
}

proc Dam::xml::ProcGetElementsValues {domNode args} {
    set TypeofProblem [get_domnode_attribute [$domNode selectNodes [spdAux::getRoute DamTypeofProblem]] v]
    set nodeApp [spdAux::GetAppIdFromNode $domNode]
    set sol_stratUN [apps::getAppUniqueName $nodeApp SolStrat]
    set schemeUN [apps::getAppUniqueName $nodeApp Scheme]
    if {[get_domnode_attribute [$domNode selectNodes [spdAux::getRoute $sol_stratUN]] v] eq ""} {
        get_domnode_attribute [$domNode selectNodes [spdAux::getRoute $sol_stratUN]] dict
    }
    if {[get_domnode_attribute [$domNode selectNodes [spdAux::getRoute $schemeUN]] v] eq ""} {
        get_domnode_attribute [$domNode selectNodes [spdAux::getRoute $schemeUN]] dict
    }
    
    set solStratName [::write::getValue $sol_stratUN]
    set schemeName [write::getValue $schemeUN]
    set elems [::Model::GetAvailableElements $solStratName $schemeName]
    
    set names [list ]
    foreach elem $elems {
        if {[$elem cumple {*}$args]} {
            lappend names [$elem getName]
        }
    }
    if {$TypeofProblem ni [list UP_Mechanical UP_Thermo-Mechanical Acoustic]} {
        set names [lsearch -all -inline -not -exact $names WaveEquationElement2D]
        set names [lsearch -all -inline -not -exact $names WaveEquationElement3D]
    }
    if {$TypeofProblem in [list Acoustic]} {
        set names [list WaveEquationElement2D]
        if {$::Model::SpatialDimension eq "3D"} {
            set names [list WaveEquationElement3D]
        }
    }
    if {$TypeofProblem in [list Modal-Analysis]} {
        set names [list SmallDisplacementElement2D]
        if {$::Model::SpatialDimension eq "3D"} {
            set names [list SmallDisplacementElement3D]
        }
    }
    if {[get_domnode_attribute $domNode v] eq ""} {$domNode setAttribute v [lindex $names 0]}
    if {[get_domnode_attribute $domNode v] ni $names} {$domNode setAttribute v [lindex $names 0]; spdAux::RequestRefresh}
    set values [join $names ","]
    return $values
}

proc Dam::xml::ProcNoorzaiState {domNode args} {
    set SourceType [get_domnode_attribute [$domNode selectNodes [spdAux::getRoute DamSourceType]] v]
    set ActivateConstruction [get_domnode_attribute [$domNode selectNodes [spdAux::getRoute DamActivateConstruction]] v]
    set state hidden
    if {$ActivateConstruction} {
        if {$SourceType == "Adiabatic"} {set state normal}
    }
    return $state
}

proc Dam::xml::ProcAzenhaState {domNode args} {
    set SourceType [get_domnode_attribute [$domNode selectNodes [spdAux::getRoute DamSourceType]] v]
    set ActivateConstruction [get_domnode_attribute [$domNode selectNodes [spdAux::getRoute DamActivateConstruction]] v]
    set state hidden
    if {$ActivateConstruction} {
        if {$SourceType != "Adiabatic"} {set state normal}
    }
    return $state
}


Dam::xml::Init



namespace eval SolidMechanics::xml {
     variable dir
}

proc SolidMechanics::xml::Init { } {
    variable dir
    Model::InitVariables dir $SolidMechanics::dir
    
    getSolutionStrategies
    getElements
    getConditions
    getConstitutiveLaws
    getSolvers
    
    getProcesses
}

proc SolidMechanics::xml::getSolutionStrategies { } {
    Model::InitVariables SolutionStrategyFileName strategydefinition.xml
    Model::getSolutionStrategies
}

proc SolidMechanics::xml::getElements { } {
    Model::InitVariables ElementsFileName Elements.xml
    Model::getElements
}

proc SolidMechanics::xml::getConditions { } {
    Model::InitVariables ConditionsFileName Conditions.xml
    Model::getConditions
}

proc SolidMechanics::xml::getConstitutiveLaws { } {
    Model::InitVariables ConstitutiveLawsFileName ConstitutiveLaws.xml
    Model::getConstitutiveLaws
}

proc SolidMechanics::xml::getSolvers { } {
    Model::InitVariables SolversFileName Solvers.xml
    Model::getSolvers
}

proc SolidMechanics::xml::getProcesses { } {
    Model::InitVariables ProcessesFileName Processes.xml
    Model::getProcesses
}

proc SolidMechanics::xml::injectSolvers {basenode} {
    
    # Get all solvers params
    set paramspuestos [list ]
    set paramsnodes ""
    set params [::Model::GetAllSolversParams]
    
    foreach {parname par} $params {
        if {$parname ni $paramspuestos} {
            lappend paramspuestos $parname
            set pn [$par getPublicName]
            set type [$par getType]
            set dv [$par getDv]
            append paramsnodes "<value n=\"$parname\" pn=\"$pn\" state=\"\[SolverParamState\]\" v=\"$dv\" "
            if {$type eq "bool"} {
                append paramsnodes " values=\"True,False\" "
            }
            if {$type eq "combo"} {
                set vals [join [$par getValues] ,]
                append paramsnodes " values=\"$vals\" "
                set d [list ]
                foreach v [$par getValues] pv [$par getPValues] {
                    lappend d "$v,$pv"
                }
                set dic [join $d ,]
                append paramsnodes " dict=\"$dic\" "
            }
            
            append paramsnodes "/>"
        }
    }
    set contnode [$basenode parent]
    
    # Get All SolversEntry
    set ses [list ]
    foreach st [::Model::GetSolutionStrategies] {
        lappend ses $st [$st getSolversEntries]
    }
    
    # One container per solverEntry 
    foreach {st se} $ses {
        set stn [$st getName]
        set n [$se getName]
        set pn [$se getPublicName]
        set help [$se getHelp]
        set container "<container help=\"$help\" n=\"$n\" pn=\"$pn\" un=\"SM$stn$n\" state=\"\[SolverEntryState\]\" solstratname=\"$stn\">"
        # Inject solvers combo
        append container "<value n=\"Solver\" pn=\"Solver\" v=\"\" values=\"\[GetSolvers\]\" actualize=\"1\">"
        append container "<dependencies node=\"../value\" actualize=\"1\"/> </value>"
        append container $paramsnodes
        append container "</container>"
        $contnode appendXML $container
    }
    $basenode delete
}

proc SolidMechanics::xml::injectSolStratParams {basenode} {
    set contnode [$basenode parent]
    set params [::Model::GetAllSolStratParams]
    foreach {parname par} $params {
        if {[$contnode find n $parname] eq ""} {
            set pn [$par getPublicName]
            set type [$par getType]
            set dv [$par getDv]
            set helptext [$par getHelp]
            set node "<value n=\"$parname\" pn=\"$pn\" state=\"\[SolStratParamState\]\" v=\"$dv\" help=\"$helptext\" "
            if {$type eq "bool"} {
                append node " values=\"True,False\" "
            } 
            append node "/>"
            catch {$contnode appendXML $node}
        }
    }
    
    set params [::Model::GetAllSchemeParams]
    
    foreach {parname par} $params {
        if {[$contnode find n $parname] eq ""} {
            set pn [$par getPublicName]
            set type [$par getType]
            set dv [$par getDv]
            set helptext [$par getHelp]
            set node "<value n=\"$parname\" pn=\"$pn\" state=\"\[SchemeParamState\]\" v=\"$dv\" help=\"$helptext\" "
            if {$type eq "bool"} {
                append node " values=\"True,False\" "
            } 
            append node "/>"
            catch {$contnode appendXML $node}
        }
    }
    $basenode delete
}


proc SolidMechanics::xml::injectProcesses {basenode} {
    
    set procsnode [$basenode parent]
    set procs [::Model::getAllProcs]
    foreach proc $procs {
        set n [$proc getName]
        set pn [$proc getPublicName]
        set help [$proc getHelp]
        
        set node "<condition n=\"$n\" pn=\"$pn\" ov=\"point,line,surface,volume\"  ovm=\"\" icon=\"shells16\" help=\"$help\"  >"
	foreach {n input} [$proc getInputs] {
		set n [$input getName]
		set pn [$input getPublicName]
		set type [$input getType]
		set v [$input getDv]
		set state [$input getAttribute state]
		#W "$n $type $v $state"
		 if {$type eq "vector"} {
			set v1 [lindex [split $v ","] 0]
			set v2 [lindex [split $v ","] 1]
			set v3 [lindex [split $v ","] 2]
		    append node "
		    <value n=\"[concat $n "_X"]\"  pn=\"X Value\" v=\"$v1\" help=\"\" />
		    <value n=\"[concat $n "_Y"]\"  pn=\"Y Value\" v=\"$v2\" help=\"\" />
		    <value n=\"[concat $n "_Z"]\"  pn=\"Z Value\" v=\"$v3\" help=\"\" />
		    "
		} elseif {$type eq "bool"} {
		    append node "<value n=\"$n\" pn=\"$n\" v=\"$v\" values=\"True,False\"  help=\"\"/>"
		} elseif {$type eq "String"} {
		    append node "<value n=\"$n\" pn=\"$n\" v=\"$v\" state=\"$state\"  help=\"\"/>"
		} else {
		    append node "<value n=\"$n\" pn=\"$n\" v=\"$v\"  help=\"\"/>"
		}
	}
        append node "</condition>"
        #W $node
        catch {$procsnode appendXML $node}
    }
    $basenode delete
}

proc SolidMechanics::xml::injectDoFs { basenode } {
    set conds [$basenode parent]
    set dofs [::Model::getAllDOFs]
    foreach n [dict keys $dofs] {
        set dof [dict get $dofs $n]
        set pn [$dof getPublicName]
        set type [$dof getType]
        set reaction [$dof getReaction]
        set units [$dof getUnits]
        set um [$dof getUnitMagnitude]
        set help [$dof getHelp]
        
        set node "<condition n=\"$n\" pn=\"$pn\" ov=\"point,line,surface,volume\" type=\"$type\" ovm=\"\" icon=\"shells16\" help=\"$help\"  state=\"\[CheckConditionContainerState\]\">"
        if {$type eq "vector"} {
            append node "
            <value n=\"FixX\" pn=\"X Fix\" v=\"1\" values=\"1,0\" help=\"\"/>
            <value n=\"FixY\" pn=\"Y Fix\" v=\"1\" values=\"1,0\" help=\"\"/>
            <value n=\"FixZ\" pn=\"Z Fix\" v=\"1\" values=\"1,0\" help=\"\"/>
            <value n=\"ValX\" wn=\"[concat $n "_X"]\" pn=\"X Value\" v=\"0.0\" help=\"\" units=\"$units\" unit_magnitude=\"$um\"/>
            <value n=\"ValY\" wn=\"[concat $n "_Y"]\" pn=\"Y Value\" v=\"0.0\" help=\"\" units=\"$units\" unit_magnitude=\"$um\"/>
            <value n=\"ValZ\" wn=\"[concat $n "_Z"]\" pn=\"Z Value\" v=\"0.0\" help=\"\" units=\"$units\" unit_magnitude=\"$um\"/>
            "
        } {
            append node "<value n=\"Value\" pn=\"Value\" v=\"1\" units=\"$units\"  unit_magnitude=\"$um\" help=\"\"/>"
        }
        append node "</condition>"
        #W $node
        catch {$conds appendXML $node}
    }
    $basenode delete
}

proc SolidMechanics::xml::injectLoads { basenode } {
    set conds [$basenode parent]
    set loads [::Model::getAllConditions]
    foreach n [dict keys $loads] {
        set ld [dict get $loads $n]
        set pn [$ld getPublicName]
        set help [$ld getHelp]
        set etype [string tolower [$ld getAttribute ElementType]]
        
        set node "<condition n=\"$n\" pn=\"$pn\" ov=\"$etype\" ovm=\"\" icon=\"shells16\" help=\"$help\" state=\"\[ConditionState\]\">"
        foreach {inName in} [$ld getInputs] {
            set inPn [$in getPublicName]
            set units [$in getUnits]
            set um [$in getUnitMagnitude]
            set type [$in getType]
            set dv [$in getDv]
            
            if {$type eq "vector"} {
                append node "
                <value n=\"ValX\" wn=\"[concat $n "_X"]\" pn=\"${inPn} X\" v=\"0.0\" help=\"\" units=\"$units\" unit_magnitude=\"$um\"/>
                <value n=\"ValY\" wn=\"[concat $n "_Y"]\" pn=\"${inPn} Y\" v=\"0.0\" help=\"\" units=\"$units\" unit_magnitude=\"$um\"/>
                <value n=\"ValZ\" wn=\"[concat $n "_Z"]\" pn=\"${inPn} Z\" v=\"0.0\" help=\"\" units=\"$units\" unit_magnitude=\"$um\" state=\"\[CheckDimension 3D\]\"/>
                "
            } {
                append node "<value n=\"$inName\" pn=\"$inPn\" v=\"$dv\"  units=\"$units\"  unit_magnitude=\"$um\"  help=\"\"/>"
            }
        }
        append node "</condition>"
        $conds appendXML $node
    }
    $basenode delete
}

proc SolidMechanics::xml::injectElementInputs { basenode } {
    set parts [$basenode parent]
    set inputs [::Model::GetAllElemInputs]
    foreach inName [dict keys $inputs] {
        set in [dict get $inputs $inName] 
        set pn [$in getPublicName]
        set units [$in getUnits]
        set um [$in getUnitMagnitude]
        set node "<value n=\"$inName\" pn=\"$pn\" state=\"\[PartParamState\]\" v=\"\[ElemParamValue\]\" units=\"$units\" unit_magnitude=\"$um\" />"
        catch {$parts appendXML $node}
    }
    $basenode delete
}
    
proc SolidMechanics::xml::injectConstitutiveLawInputs { basenode } {
    set parts [$basenode parent]
    set inputs [::Model::GetAllCLInputs]
    foreach inName [dict keys $inputs] {
        if {[$parts find n $inName] eq ""} {
            set in [dict get $inputs $inName] 
            set pn [$in getPublicName]
            set units [$in getUnits]
            set um [$in getUnitMagnitude]
            set node "<value n=\"$inName\" pn=\"$pn\" state=\"\[PartParamState\]\" v=\"\[ConstLawParamValue\]\" units=\"$units\" unit_magnitude=\"$um\" />"
            catch {$parts appendXML $node}
        }
    }
    $basenode delete
}

proc SolidMechanics::xml::injectElementOutputs { basenode } {
    set parts [$basenode parent]
    set outputs [::Model::GetAllElemOutputs]
    foreach in [dict keys $outputs] {
        set pn [[dict get $outputs $in] getPublicName]
        set node "<value n=\"$in\" pn=\"$pn\" state=\"\[ElementOutputState\]\" v=\"True\" values=\"True,False\" />"
        catch {$parts appendXML $node}
    }
    $basenode delete
}
    
proc SolidMechanics::xml::injectConstitutiveLawOutputs { basenode } {
    set parts [$basenode parent]
    set outputs [::Model::GetAllCLOutputs]
    foreach in [dict keys $outputs] {
        if {[$parts find n $in] eq ""} {
        set pn [[dict get $outputs $in] getPublicName]
            set node "<value n=\"$in\" pn=\"$pn\" state=\"\[ConstLawOutputState\]\" v=\"True\" values=\"True,False\" />"
            catch {$parts appendXML $node}
        }
    }
    $basenode delete
}

proc SolidMechanics::xml::injectProcs { basenode } {
    set nf [file join $::SolidMechanics::dir xml Procs.spd]
    set xml [tDOM::xmlReadFile $nf]
    set newnode [dom parse [string trim $xml]]
    set xmlNode [$newnode documentElement]

    foreach in [$xmlNode getElementsByTagName "proc"] {
        [$basenode parent] appendChild $in
    }
    $basenode delete
}

proc SolidMechanics::xml::CheckConstLawOutputState {outnode} {
    set doc $gid_groups_conds::doc
    set root [$doc documentElement]
    
    set xp1 "[spdAux::getRoute "SMParts"]/group/value\[@n='ConstitutiveLaw'\]"
    set constlawactive [list ]
    foreach gNode [$root selectNodes $xp1] {
        lappend constlawactive [$gNode @v]
    }
    
    set paramName [$outnode @n]
    return [::Model::CheckConstLawOutputState $constlawactive $paramName]
}

proc SolidMechanics::xml::CheckElementOutputState {outnode} {
    set doc $gid_groups_conds::doc
    set root [$doc documentElement]
    
    set xp1 "[spdAux::getRoute "SMParts"]/group/value\[@n='Element'\]"
    set elemsactive [list ]
    foreach gNode [$root selectNodes $xp1] {
        lappend elemsactive [$gNode @v]
    }
    
    set paramName [$outnode @n]
    return [::Model::CheckElementOutputState $elemsactive $paramName]
}

proc SolidMechanics::xml::SolStratParamState {outnode} {
    set doc $gid_groups_conds::doc
    set root [$doc documentElement]
    if {[get_domnode_attribute [$root selectNodes [spdAux::getRoute SMSolStrat]] v] eq ""} {
        get_domnode_attribute [$root selectNodes [spdAux::getRoute SMSolStrat]] values
    }
    set SolStrat [get_domnode_attribute [$root selectNodes [spdAux::getRoute SMSolStrat]] v]
    set paramName [$outnode @n]
    return [::Model::CheckSolStratInputState $SolStrat $paramName]
}

proc SolidMechanics::xml::SchemeParamState {outnode} {
    set doc $gid_groups_conds::doc
    set root [$doc documentElement]
    
    set xp1 "[spdAux::getRoute "SMSolStrat"]"
    set SolStrat [[$root selectNodes $xp1 ] @v]
    set xp1 "[spdAux::getRoute "SMScheme"]"
    set Scheme [[$root selectNodes $xp1 ] @v]
    
    set paramName [$outnode @n]
    return [::Model::CheckSchemeInputState $SolStrat $Scheme $paramName]
}

proc SolidMechanics::xml::getIntervals {} {
    set doc $gid_groups_conds::doc
    set root [$doc documentElement]
    set xp1 "[spdAux::getRoute "SMIntervals"]/blockdata\[@n='Interval'\]"
    set lista [list ]
    foreach node [$root selectNodes $xp1] {
        lappend lista [$node @name]
    }
    
    return $lista
}

proc SolidMechanics::xml::getTimeFunctions {} {
    set doc $gid_groups_conds::doc
    set root [$doc documentElement]
    set xp1 "[spdAux::getRoute "SMFunctions"]/blockdata\[@n='Function'\]"
    set lista [list ]
    foreach node [$root selectNodes $xp1] {
        lappend lista [$node @name]
    }
    
    return $lista
}

proc SolidMechanics::xml::getFields {} {
    set doc $gid_groups_conds::doc
    set root [$doc documentElement]
    set xp1 "[spdAux::getRoute "SMFields"]/blockdata\[@n='Field'\]"
    set lista [list ]
    foreach node [$root selectNodes $xp1] {
        lappend lista [$node @name]
    }
    
    return $lista
}

SolidMechanics::xml::Init

##################################################################################
#   This file is common for all Kratos Applications.
#   Do not change anything here unless it's strictly necessary.
##################################################################################

namespace eval spdAux {
    # Namespace variables declaration
    
    variable uniqueNames
    variable initwind
    
    variable currentexternalfile
    variable refreshTreeTurn
}

proc spdAux::Init { } {
    # Namespace variables inicialization
    variable uniqueNames
    variable initwind
    variable currentexternalfile
    variable refreshTreeTurn
    
    set uniqueNames ""
    dict set uniqueNames "dummy" 0
    set initwind ""
    set  currentexternalfile ""
    set refreshTreeTurn 0
    #spdAux::TryRefreshTree
}

proc spdAux::RequestRefresh {} {
        variable refreshTreeTurn
        set refreshTreeTurn 1
}

proc spdAux::TryRefreshTree { } {
        variable refreshTreeTurn
        #W "HI"
        if {$refreshTreeTurn} {
            #W "there"
            set ::spdAux::refreshTreeTurn 0
            gid_groups_conds::actualize_conditions_window
            set ::spdAux::refreshTreeTurn 0
        }
        after 750 {spdAux::TryRefreshTree}
}

proc spdAux::EndRefreshTree { } {
        variable refreshTreeTurn
        set refreshTreeTurn 0
        after cancel {spdAux::TryRefreshTree}
}

# Includes
proc spdAux::processIncludes { } {
    set doc $gid_groups_conds::doc
    set root [$doc documentElement]
    spdAux::processAppIncludes $root
    spdAux::processDynamicNodes $root
}

proc spdAux::processDynamicNodes { root } {
    foreach elem [$root getElementsByTagName "dynamicnode"] {
        set func [$elem getAttribute command]
        spdAux::${func} $elem
    }
}

proc spdAux::processAppIncludes { root } {
    foreach elem [$root getElementsByTagName "appLink"] {
        set active [$elem getAttribute "active"]
        set appid [$elem getAttribute "appid"]
        set pn [$elem getAttribute "pn"]
        set prefix [$elem getAttribute "prefix"]
        apps::NewApp $appid $pn $prefix
        if {$active} {
            set dir $::Kratos::kratos_private(Path)
            set f [file join $dir apps $appid xml Main.spd]
            # Only keep pn2 when minimum GiD version > 12.1.11d
            set pn1 ""; set pn2 "";
            catch {set pn1 [customlib::processIncludesRecurse $f $dir]; set processedAppnode $pn1}
            catch {set pn2 [customlib::ProcessIncludesRecurse $f $dir]; set processedAppnode $pn2}
            #set processedAppnode [customlib::processIncludesRecurse $f $dir]
            $root insertBefore $processedAppnode $elem
            $elem delete
        }
    }
}

proc spdAux::reactiveApp { } {
    variable initwind
    destroy $initwind
    set doc $gid_groups_conds::doc
    set root [$doc documentElement]
    set ::Model::SpatialDimension [[$root selectNodes "hiddenfield\[@n='nDim'\]"] getAttribute v ]
    set appname [[$root selectNodes "hiddenfield\[@n='activeapp'\]"] @v ]
    apps::setActiveApp $appname
}

proc spdAux::activeApp { appid } {
    variable initwind
    destroy $initwind
    set doc $gid_groups_conds::doc
    set root [$doc documentElement]
    [$root selectNodes "hiddenfield\[@n='activeapp'\]"] setAttribute v $appid
    foreach elem [$root getElementsByTagName "appLink"] {
        if {$appid eq [$elem getAttribute "appid"] && [$elem getAttribute "active"] eq "0"} {
            $elem setAttribute "active" 1
        } else {
            $elem setAttribute "active" 0
        }
    }
    [$root selectNodes "hiddenfield\[@n='nDim'\]"] setAttribute v $::Model::SpatialDimension
    spdAux::processIncludes
    parseRoutes
    catch {apps::ExecuteOnCurrent init MultiAppEvent}
    gid_groups_conds::actualize_conditions_window
}

proc spdAux::CreateWindow {dir} {
    variable initwind
    set doc $gid_groups_conds::doc
    set root [$doc documentElement]
    
    set activeapp [ [$root selectNodes "hiddenfield\[@n='activeapp'\]"] getAttribute v]
        
    if { $activeapp ne "" } {
        reactiveApp
        catch {destroy $initwind}
        return ""
    }
    if { [ winfo exist .gid.win_example]} {
        destroy .gid.win_example        
    }   
    
    toplevel .gid.win_example
    wm withdraw .gid.win_example

    set w .gid.win_example
    
    set x [expr [winfo rootx .gid]+[winfo width .gid]/2-[winfo width $w]/2]
    set y [expr [winfo rooty .gid]+[winfo height .gid]/2-[winfo height $w]/2]
    
    wm geom .gid.win_example +$x+$y
    wm transient .gid.win_example .gid    

    InitWindow $w [_ "Kratos Multiphysics"] Kratos "" "" 1
    set initwind $w
    ttk::frame $w.top
    ttk::label $w.top.title_text -text [_ " Application market"]
   
    ttk::frame $w.information  -relief ridge 
    
    set appsid [::apps::getAllApplicationsID]
    set appspn [::apps::getAllApplicationsName]
    
    
    set col 0
    set row 0
    foreach appname $appspn appid $appsid {
        set img [::apps::getImgFrom $appid]
        ttk::button $w.information.img$appid -image $img -command [list apps::setActiveApp $appid]
        ttk::label $w.information.text$appid -text $appname
        
        grid $w.information.img$appid -column $col -row $row
        grid $w.information.text$appid -column $col -row [expr $row +1]
        
        incr col
        if {$col >= 5} {set col 0; incr row; incr row}
    }
    
    set frsd [ttk::frame $w.sd]
    set sd [ttk::label $frsd.lgl -text "Spatial Dimension:"]
    set vs [list "2D" "3D"]
    set sdcmb [ttk::combobox $frsd.sdcmb -textvariable ::Model::SpatialDimension -values $vs -width 4 -state "readonly"]
    
    
    grid $w.top
    grid $w.top.title_text
    
    grid $w.information
    grid $sd -row 0 -column 0 -padx 20
    grid $sdcmb  -row 0 -column 1 -padx 20
    grid $frsd -sticky we
}
proc spdAux::DestroyWindow {} {
    variable initwind
    
    catch {destroy $initwind}
    return ""
}

# Routes
proc spdAux::getRoute {name} {
    variable uniqueNames
    set v ""
    catch {
        set v [dict get $uniqueNames $name]
    }
    return $v
}
proc spdAux::setRoute {name route} {
    variable uniqueNames
    #if {[dict exists $uniqueNames $name]} {W "Warning: Unique name $name already exists.\n    Previous value: [dict get $uniqueNames $name],\n    Updated value: $route"}
    set uniqueNames [dict set uniqueNames $name $route]
        
        set uniqueNames [dict remove $uniqueNames dummy]
    # W "Add $name $route"
    # set doc $gid_groups_conds::doc
    # set root [$doc documentElement]
    # W "checking [[$root selectNodes $route] asXML]"
}

proc spdAux::parseRoutes { } {
    set doc $gid_groups_conds::doc
    set root [$doc documentElement]
    parseRecurse $root
}

proc spdAux::parseRecurse { root } {
    foreach node [$root childNodes] {
        #catch {W "Nombre: [$node getAttribute n] Tiene uniquename: [$node hasAttribute un]"}
        catch {
            #W [$node asXML]
            if {[$node hasAttribute un] == "1"} {
                setRoute [$node getAttribute un] [$node toXPath]
            }
        }
        if {[$node hasChildNodes]} {
            parseRecurse $node
        }
    }
}


proc spdAux::ExploreAllRoutes { } {
    variable uniqueNames
    
    set root [$gid_groups_conds::doc documentElement]
    W [dict keys $uniqueNames]
    foreach routeName [dict keys $uniqueNames] {
        set route [getRoute $routeName]
        W "Route $routeName $route"
        set node [$root selectNodes $route]
        W "Node $node"
        W "Value [get_domnode_attribute $node values]"
        W "Value [get_domnode_attribute $node v]"
        gid_groups_conds::uncompress_subtree $node
         
        W "Value [get_domnode_attribute $node values]"
        W "Value [get_domnode_attribute $node v]"
    }
        
}
# Dependencies
proc spdAux::insertDependencies { baseNode originUN } {
    set doc $gid_groups_conds::doc
    set root [$doc documentElement]
    set originxpath [$baseNode toXPath]
    set insertxpath [getRoute $originUN]
    set insertonnode [$root selectNodes $insertxpath]
    # a lo bestia, cambiar cuando sepamos inyectar la dependencia, abajo esta a medias
    $insertonnode setAttribute "actualize_tree" 1
    
    # Aun no soy capaz de insertar y que funcione
    #set ready 1
    #foreach c [$insertonnode getElementsByTagName "dependencies"] {
    #    if {[$c getAttribute "node"] eq $originxpath} {set ready 0; break}
    #}
    #
    #if {$ready} {
    #    set str "<dependencies node=\"$originxpath\" actualize=\"1\"/>"
    #    W $str
    #    W $insertxpath
    #    $insertonnode appendChild [[dom parse $str] documentElement]
    #    W [$insertonnode asXML]
    #}
}

proc spdAux::CheckSolverEntryState {domNode} {
        set kw [apps::getCurrentUniqueName SolStrat]
        set nodo [$domNode selectNodes [getRoute $kw]]
        get_domnode_attribute $nodo values
        set currentSolStrat [get_domnode_attribute $nodo v]
        set mySolStrat [get_domnode_attribute $domNode solstratname]
        return [expr [string compare $currentSolStrat $mySolStrat] == 0]
}



proc spdAux::CheckConstLawParamValue {node} {
    set doc $gid_groups_conds::doc
    set root [$doc documentElement]
    
    set id [$node getAttribute n]
    set val 0.0
    
    # JG CHAPUZA A RESOLVER
    set material_name [get_domnode_attribute [[$node parent] selectNodes "./value\[@n='Material'\]"] v]
    set material_name [.gid.central.boundaryconds.gg.data.f0.e2 get]
    
    set mats_un [apps::getCurrentUniqueName Materials]
    set xp3 [spdAux::getRoute $mats_un]
    append xp3 [format_xpath {/blockdata[@n="material" and @name=%s]/value} $material_name]

    foreach valueNode [$root selectNodes $xp3] {
        if {$id eq [$valueNode getAttribute n] } {set val [$valueNode getAttribute v]}
    }   
    #W "mat: $material_name prop $id val $val"
    
    return $val
}


proc spdAux::chk_loads_function_time { domNode load_name } {
    set loads [list [list scalar]]
    lappend loads [list interpolator_func x x T]
    return [join $loads ,]
}

proc spdAux::ViewDoc {} {
    set doc $gid_groups_conds::doc
    set root [$doc documentElement]
    
    W [$root asXML]
}


proc spdAux::CheckElemParamValue {node} {
    set doc $gid_groups_conds::doc
    set root [$doc documentElement]
    
    set id [$node getAttribute n]
    set val 0.0
    set material_name [get_domnode_attribute [[$node parent] selectNodes "./value\[@n='Material'\]"] v]
    set material_name [.gid.central.boundaryconds.gg.data.f0.e2 get]
    
    set mats_un [apps::getCurrentUniqueName Materials]
    set xp3 [spdAux::getRoute $mats_un]
    append xp3 [format_xpath {/blockdata[@n="material" and @name=%s]/value} $material_name]

    foreach valueNode [$root selectNodes $xp3] {
        if {$id eq [$valueNode getAttribute n] } {set val [$valueNode getAttribute v]}
    }
    
    return $val
}

proc spdAux::ConvertAllUniqueNames {oldPrefix newPrefix} {
    variable uniqueNames
    set root [$gid_groups_conds::doc documentElement]

    foreach routeName [dict keys $uniqueNames] {
        if {[string first $oldPrefix $routeName] eq 0} {
            set route [getRoute $routeName]
            set newrouteName [string map [list $oldPrefix $newPrefix] $routeName]
            [$root selectNodes $route] setAttribute un $newrouteName
        }
    }
    spdAux::parseRoutes
}

proc spdAux::ListToValues {lista} {
    set res ""
    foreach elem $lista {
        append res $elem
        append res ","
    }
    return [string range $res 0 end-1]
}

proc spdAux::injectSolvers {basenode} {
    
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
                append paramsnodes " values=\"Yes,No\" "
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
    foreach {st ss} $ses {
        foreach se $ss {
            set stn [$st getName]
            set n [$se getName]
            set pn [$se getPublicName]
            set help [$se getHelp]
            set un [apps::getCurrentUniqueName "$stn$n"]
            set container "<container help=\"$help\" n=\"$n\" pn=\"$pn\" un=\"$un\" state=\"\[SolverEntryState\]\" solstratname=\"$stn\" >"
            # Inject solvers combo
            append container "<value n=\"Solver\" pn=\"Solver\" v=\"\" values=\"\[GetSolvers\]\" actualize=\"1\" update_proc=\"Updateme\"/>"
            #append container "<dependencies node=\"../value\" actualize=\"1\"/>"
            #append container "</value>"
            append container $paramsnodes
            append container "</container>"
            $contnode appendXML $container
        }
    }
    $basenode delete
}

proc spdAux::injectSolStratParams {basenode} {
    set contnode [$basenode parent]
    set params [::Model::GetAllSolStratParams]
    foreach {parname par} $params {
        if {[$contnode find n $parname] eq ""} {
            set pn [$par getPublicName]
            set type [$par getType]
            set dv [$par getDv]
            set helptext [$par getHelp]
            set actualize [$par getActualize]
            set node "<value n=\"$parname\" pn=\"$pn\" state=\"\[SolStratParamState\]\" v=\"$dv\" help=\"$helptext\" "
            
            if {$actualize} {
                append node "actualize_tree=\"1\""
            }
            
            if {$type eq "bool"} {
                append node " values=\"Yes,No\" "
            }
            if {$type eq "combo"} {
                set values [$par getValues]
                set vs [join [$par getValues] ,]
                set pvalues [$par getPValues]
                
                set pv ""
                for {set i 0} {$i < [llength $values]} {incr i} {
                    append pv [lindex $values $i] "," [lindex $pvalues $i] ","
                }
                if {[llength $pv]} {set pv [string range $pv 0 end-1]}
                #set pv [join $values $pvalues ,]
                append node " values=\"$vs\" dict=\"$pv\" "
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
                append node " values=\"Yes,No\" "
            } 
            append node "/>"
            catch {$contnode appendXML $node}
        }
    }
    $basenode delete
}


proc spdAux::injectProcesses {basenode} {
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
                append node "<value n=\"$n\" pn=\"$n\" v=\"$v\" values=\"Yes,No\"  help=\"\"/>"
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

proc spdAux::injectNodalConditions { basenode } {
    set nodalconds [$basenode parent]
    set nodal_conditions [::Model::getAllNodalConditions]
    foreach n [dict keys $nodal_conditions] {
        set nc [dict get $nodal_conditions $n]
        set pn [$nc getPublicName]
        set help [$nc getHelp]
        set ov  [$nc getOv]
        set node "<condition n=\"$n\" pn=\"$pn\" ov=\"$ov\"  ovm=\"\" icon=\"shells16\" help=\"$help\" state=\"\[CheckNodalConditionState\]\">"
        set inputs [$nc getInputs]
        set process [::Model::GetProcess [$nc getProcessName]]
        set unitsnc [$nc getAttribute "units"]
        set umnc [$nc getAttribute "unit_magnitude"]
        foreach processinput [$process getInputs] {lappend inputs $processinput}
        foreach {inName in} $inputs {
            set inPn [$in getPublicName]
            set units [$in getUnits]
            if {$units eq "0"} {set units $unitsnc}
            set um [$in getUnitMagnitude]
            if {$um eq "0"} {set um $umnc}
            set type [$in getType]
            set dv [$in getDv]
            set fix [$in getFixity]
            set help [$in getHelp]
            #W "$pn fix -> $fix"
            if {$type eq "vector"} {
                if {$fix ne 0} {
                    append node "
                        <value n=\"FixX\" pn=\"$fix X\" v=\"1\" values=\"1,0\" help=\"\"/>
                        <value n=\"FixY\" pn=\"$fix Y\" v=\"1\" values=\"1,0\" help=\"\"/>
                        <value n=\"FixZ\" pn=\"$fix Z\" v=\"1\" values=\"1,0\" help=\"\"/>"
                    }
                append node "
                <value n=\"${inName}X\" wn=\"[concat $n "_X"]\" pn=\"${inPn} X\" v=\"0.0\" help=\"$help\" units=\"$units\" unit_magnitude=\"$um\"/>
                <value n=\"${inName}Y\" wn=\"[concat $n "_Y"]\" pn=\"${inPn} Y\" v=\"0.0\" help=\"$help\" units=\"$units\" unit_magnitude=\"$um\"/>
                <value n=\"${inName}Z\" wn=\"[concat $n "_Z"]\" pn=\"${inPn} Z\" v=\"0.0\" help=\"$help\" units=\"$units\" unit_magnitude=\"$um\" state=\"\[CheckDimension 3D\]\"/>
                "
            } {
                if {$fix ne 0} {
                    append node "<value n=\"Fix\" pn=\"$fix\" v=\"1\" values=\"1,0\" help=\"\"/>"
                }
                append node "<value n=\"$inName\" pn=\"$inPn\" v=\"$dv\"  units=\"$units\"  unit_magnitude=\"$um\"  help=\"$help\"/>"
            }
        }
        append node "</condition>"
        $nodalconds appendXML $node
    }
    $basenode delete
}

proc spdAux::injectConditions { basenode } {
    set conds [$basenode parent]
    set loads [::Model::getAllConditions]
    foreach n [dict keys $loads] {
        set ld [dict get $loads $n]
        set pn [$ld getPublicName]
        set help [$ld getHelp]
        set etype [string tolower [$ld getAttribute ElementType]]
        set node "<condition n=\"$n\" pn=\"$pn\" ov=\"$etype\" ovm=\"\" icon=\"shells16\" help=\"$help\" state=\"\[ConditionState\]\">"
        set inputs [$ld getInputs]
        set unitsc [$ld getAttribute "units"]
        set umc [$ld getAttribute "unit_magnitude"]
        set process [::Model::GetProcess [$ld getProcessName]]
        foreach processinput [$process getInputs] {lappend inputs $processinput}
        foreach {inName in} $inputs {
            set inPn [$in getPublicName]
            set units [$in getUnits]
            if {$units eq "0"} {set units $unitsc}
            set um [$in getUnitMagnitude]
            if {$um eq "0"} {set um $umc}
            set type [$in getType]
            set dv [$in getDv]
            set fix [$in getFixity]
            set help [$in getHelp]
            if {$type eq "vector"} {
                if {$fix ne 0} {
                    append node "
                        <value n=\"FixX\" pn=\"X $fix\" v=\"1\" values=\"1,0\" help=\"\"/>
                        <value n=\"FixY\" pn=\"Y $fix\" v=\"1\" values=\"1,0\" help=\"\"/>
                        <value n=\"FixZ\" pn=\"Z $fix\" v=\"1\" values=\"1,0\" help=\"\"/>"
                    }
                append node "
                <value n=\"${inName}X\" wn=\"[concat $n "_X"]\" pn=\"${inPn} X\" v=\"0.0\" help=\"$help\" units=\"$units\" unit_magnitude=\"$um\"/>
                <value n=\"${inName}Y\" wn=\"[concat $n "_Y"]\" pn=\"${inPn} Y\" v=\"0.0\" help=\"$help\" units=\"$units\" unit_magnitude=\"$um\"/>
                <value n=\"${inName}Z\" wn=\"[concat $n "_Z"]\" pn=\"${inPn} Z\" v=\"0.0\" help=\"$help\" units=\"$units\" unit_magnitude=\"$um\" state=\"\[CheckDimension 3D\]\"/>
                "
            } {
                if {$fix ne 0} {
                    append node "<value n=\"Fix\" pn=\"$fix\" v=\"1\" values=\"1,0\" help=\"\"/>"
                }
                append node "<value n=\"$inName\" pn=\"$inPn\" v=\"$dv\"  units=\"$units\"  unit_magnitude=\"$um\"  help=\"$help\"/>"
            }
        }
        append node "</condition>"
        $conds appendXML $node
    }
    $basenode delete
}

proc spdAux::injectElementInputs { basenode } {
    set parts [$basenode parent]
    set inputs [::Model::GetAllElemInputs]
    foreach inName [dict keys $inputs] {
        set in [dict get $inputs $inName] 
        set pn [$in getPublicName]
        set units [$in getUnits]
        set um [$in getUnitMagnitude]
        set help [$in getHelp] 
        set node "<value n=\"$inName\" pn=\"$pn\" state=\"\[PartParamState\]\" v=\"\[ElemParamValue\]\" units=\"$units\" unit_magnitude=\"$um\" help=\"$help\" />"
        catch {$parts appendXML $node}
    }
    $basenode delete
}
    
proc spdAux::injectConstitutiveLawInputs { basenode } {
    set parts [$basenode parent]
    set inputs [::Model::GetAllCLInputs]
    foreach inName [dict keys $inputs] {
        if {[$parts find n $inName] eq ""} {
            set in [dict get $inputs $inName] 
            set pn [$in getPublicName]
            set units [$in getUnits]
            set um [$in getUnitMagnitude]
            set help [$in getHelp]
            set node "<value n=\"$inName\" pn=\"$pn\" state=\"\[PartParamState\]\" v=\"\[ConstLawParamValue\]\" units=\"$units\" unit_magnitude=\"$um\" help=\"$help\" />"
            catch {$parts appendXML $node}
        }
    }
    $basenode delete
}

proc spdAux::injectElementOutputs { basenode } {
    set parts [$basenode parent]
    set outputs [::Model::GetAllElemOutputs]
    foreach in [dict keys $outputs] {
        set pn [[dict get $outputs $in] getPublicName]
        set v [GetBooleanForTree [[dict get $outputs $in] getDv]]
        
        set node "<value n=\"$in\" pn=\"$pn\" state=\"\[ElementOutputState\]\" v=\"$v\" values=\"Yes,No\" />"
        catch {$parts appendXML $node}
    }
    $basenode delete
}

proc spdAux::GetBooleanForTree {raw} {
    set goodList [list "Yes" "1" "yes" "ok" "YES" "Ok" "True" "TRUE" "true"]
    if {$raw in $goodList} {return "Yes" } {return "No"}
}
    
proc spdAux::injectConstitutiveLawOutputs { basenode } {
    set parts [$basenode parent]
    set outputs [::Model::GetAllCLOutputs]
    foreach in [dict keys $outputs] {
        if {[$parts find n $in] eq ""} {
            set pn [[dict get $outputs $in] getPublicName]
            set v [GetBooleanForTree [[dict get $outputs $in] getDv]]
            set node "<value n=\"$in\" pn=\"$pn\" state=\"\[ConstLawOutputState\]\" v=\"$v\" values=\"Yes,No\" />"
            catch {$parts appendXML $node}
        }
    }
    $basenode delete
}

proc spdAux::injectProcs { basenode } {
    set appId [apps::getActiveAppId]
    if {$appId ne ""} {
        set f "::$appId"
        append f "::dir"
        set nf [file join [subst $$f] xml Procs.spd]
        set xml [tDOM::xmlReadFile $nf]
        set newnode [dom parse [string trim $xml]]
        set xmlNode [$newnode documentElement]
    
        foreach in [$xmlNode getElementsByTagName "proc"] {
            # This allows an app to overwrite mandatory procs
            set procn [$in @n]
            catch {
                set pastnode [[$basenode parent] selectNodes "./proc\[@n='$procn'\]"]
                $pastnode delete
            }
            [$basenode parent] appendChild $in
        }
        $basenode delete
    }
}

proc spdAux::CheckConstLawOutputState {outnode} {
    set doc $gid_groups_conds::doc
    set root [$doc documentElement]
    
    set parts_un [apps::getCurrentUniqueName Parts]
    set parts_path [getRoute $parts_un]
    set xp1 "$parts_path/group/value\[@n='ConstitutiveLaw'\]"
    set constlawactive [list ]
    foreach gNode [$root selectNodes $xp1] {
        lappend constlawactive [get_domnode_attribute $gNode v]
    }
    
    set paramName [$outnode @n]
    return [::Model::CheckConstLawOutputState $constlawactive $paramName]
}

proc spdAux::CheckElementOutputState {outnode} {
    set doc $gid_groups_conds::doc
    set root [$doc documentElement]
    
    set parts_un [apps::getCurrentUniqueName Parts]
    set parts_path [getRoute $parts_un]
    set xp1 "$parts_path/group/value\[@n='Element'\]"
    set elemsactive [list ]
    foreach gNode [$root selectNodes $xp1] {
        lappend elemsactive [get_domnode_attribute $gNode v]
    }
    set paramName [$outnode @n]
    return [::Model::CheckElementOutputState $elemsactive $paramName]
}

proc spdAux::SolStratParamState {outnode} {
    set doc $gid_groups_conds::doc
    set root [$doc documentElement]
    set solstrat_un [apps::getCurrentUniqueName SolStrat]
    #W $solstrat_un
    if {[get_domnode_attribute [$root selectNodes [spdAux::getRoute $solstrat_un]] v] eq ""} {
        get_domnode_attribute [$root selectNodes [spdAux::getRoute $solstrat_un]] values
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

proc spdAux::SchemeParamState {outnode} {
    set doc $gid_groups_conds::doc
    set root [$doc documentElement]
    
    set solstrat_un [apps::getCurrentUniqueName SolStrat]
    set schemet_un [apps::getCurrentUniqueName SolStrat]
    set xp1 "[spdAux::getRoute $solstrat_un]"
    set SolStrat [[$root selectNodes $xp1 ] @v]
    set xp1 "[spdAux::getRoute $schemet_un]"
    set Scheme [[$root selectNodes $xp1 ] @v]
    
    set paramName [$outnode @n]
    return [::Model::CheckSchemeInputState $SolStrat $Scheme $paramName]
}

proc spdAux::getIntervals {} {
    set doc $gid_groups_conds::doc
    set root [$doc documentElement]
    set intervals_un [apps::getCurrentUniqueName Intervals]
    set xp1 "[spdAux::getRoute $intervals_un]/blockdata\[@n='Interval'\]"
    set lista [list ]
    foreach node [$root selectNodes $xp1] {
        lappend lista [$node @name]
    }
    
    return $lista
}

proc spdAux::getTimeFunctions {} {
    set doc $gid_groups_conds::doc
    set root [$doc documentElement]
    set functions_un [apps::getCurrentUniqueName Functions]
    set xp1 "[spdAux::getRoute $functions_un]/blockdata\[@n='Function'\]"
    set lista [list ]
    foreach node [$root selectNodes $xp1] {
        lappend lista [$node @name]
    }
    
    return $lista
}

proc spdAux::getFields {} {
    set doc $gid_groups_conds::doc
    set root [$doc documentElement]
    set fields_un [apps::getCurrentUniqueName Fields]
    set xp1 "[spdAux::getRoute $fields_un]/blockdata\[@n='Field'\]"
    set lista [list ]
    foreach node [$root selectNodes $xp1] {
        lappend lista [$node @name]
    }
    
    return $lista
}


spdAux::Init

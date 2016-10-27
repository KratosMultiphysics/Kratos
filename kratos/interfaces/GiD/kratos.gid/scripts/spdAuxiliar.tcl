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
    
    variable TreeVisibility
}

proc spdAux::Init { } {
    # Namespace variables inicialization
    variable uniqueNames
    variable initwind
    variable currentexternalfile
    variable refreshTreeTurn
    variable TreeVisibility
    
    set uniqueNames ""
    dict set uniqueNames "dummy" 0
    set initwind ""
    set  currentexternalfile ""
    set refreshTreeTurn 0
    set TreeVisibility 1
    #spdAux::TryRefreshTree
}

proc spdAux::RequestRefresh {} {
    variable refreshTreeTurn
    set refreshTreeTurn 1
}

proc spdAux::TryRefreshTree { } {
    variable refreshTreeTurn
    #W "HI"
    update
    update idletasks
    if {$refreshTreeTurn} {
        #W "there"
        catch {
            set foc [focus]
            set ::spdAux::refreshTreeTurn 0
            gid_groups_conds::actualize_conditions_window
            gid_groups_conds::actualize_conditions_window
            focus -force $foc
        }
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
        set ar [$elem getAttribute args]
        spdAux::${func} $elem $ar
    }
}

proc spdAux::processAppIncludes { root } {
    foreach elem [$root getElementsByTagName "appLink"] {
        set active [$elem getAttribute "active"]
        set appid [$elem getAttribute "appid"]
        set pn [$elem getAttribute "pn"]
        set prefix [$elem getAttribute "prefix"]
        set public 0
        catch {set public [$elem getAttribute "public"]}
        set app [apps::NewApp $appid $pn $prefix]
        $app setPublic $public
        if {$active} {
            set dir $::Kratos::kratos_private(Path)
            set f [file join $dir apps $appid xml Main.spd]
            set processedAppnode [customlib::ProcessIncludesRecurse $f $dir]
            $root insertBefore $processedAppnode $elem
            $elem delete
        }
    }
}

proc spdAux::reactiveApp { } {
    #W "Reactive"
    variable initwind
    destroy $initwind
    set doc $gid_groups_conds::doc
    set root [$doc documentElement]
    set ::Model::SpatialDimension [[$root selectNodes "value\[@n='nDim'\]"] getAttribute v ]
    set appname [[$root selectNodes "hiddenfield\[@n='activeapp'\]"] @v ]
    #spdAux::activeApp $appname
    apps::setActiveApp $appname
}

proc spdAux::activeApp { appid } {
    #W "Active $appid"
    variable initwind
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
    set nd ""
    catch {set nd [$root selectNodes "value\[@n='nDim'\]"]} 
    if {$nd eq ""} {catch {set nd [$root selectNodes "hiddenfield\[@n='nDim'\]"]}}
    
    if {[$nd getAttribute v] ne "undefined"} {
        spdAux::SwitchDimAndCreateWindow $::Model::SpatialDimension
        #destroy $initwind
        #gid_groups_conds::open_conditions menu
        #gid_groups_conds::actualize_conditions_window
        #spdAux::RequestRefresh
        spdAux::TryRefreshTree
    }
}

proc spdAux::CreateWindow {} {
    variable initwind
    set doc $gid_groups_conds::doc
    set root [$doc documentElement]
    
    set activeapp [ [$root selectNodes "hiddenfield\[@n='activeapp'\]"] getAttribute v]
    
    if { $activeapp ne "" } {
        #W "Reactivando $activeapp"
        apps::setActiveApp $activeapp
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
        if {[apps::isPublic $appid]} {
            set img [::apps::getImgFrom $appid]
            ttk::button $w.information.img$appid -image $img -command [list apps::setActiveApp $appid]
            ttk::label $w.information.text$appid -text $appname
            
            grid $w.information.img$appid -column $col -row $row
            grid $w.information.text$appid -column $col -row [expr $row +1]
            
            incr col
            if {$col >= 5} {set col 0; incr row; incr row}
        }
    }
    
    grid $w.top
    grid $w.top.title_text
    
    grid $w.information
}

proc spdAux::CreateDimensionWindow { } {
    #package require anigif 1.3
    variable initwind
    set doc $gid_groups_conds::doc
    set root [$doc documentElement]
    
    set nd ""
    catch {set nd [ [$root selectNodes "value\[@n='nDim'\]"] getAttribute v]} 
    if {$nd eq ""} {catch {set nd [ [$root selectNodes "hiddenfield\[@n='nDim'\]"] getAttribute v]}}
    if { $nd ne "undefined" } {
        spdAux::SwitchDimAndCreateWindow $nd
    } {
        set dir $::Kratos::kratos_private(Path)
        
        set initwind .gid.win_example
        if { [ winfo exist $initwind]} {
            destroy $initwind
        }
        toplevel $initwind
        wm withdraw $initwind
        
        set w $initwind
        
        set x [expr [winfo rootx .gid]+[winfo width .gid]/2-[winfo width $w]/2]
        set y [expr [winfo rooty .gid]+[winfo height .gid]/2-[winfo height $w]/2]
        
        wm geom $initwind +$x+$y
        wm transient $initwind .gid    
        
        InitWindow $w [_ "Kratos Multiphysics"] Kratos "" "" 1
        set initwind $w
        ttk::frame $w.top
        ttk::label $w.top.title_text -text [_ " Dimension selection"]
        
        ttk::frame $w.information  -relief ridge
        set i 0
        foreach dim $::Model::ValidSpatialDimensions {
            set imagepath [getImagePathDim $dim]
            if {![file exists $imagepath]} {set imagepath [file nativename [file join $dir images "$dim.png"]]}
            set img [gid_themes::GetImageModule $imagepath ""]
            #W [file extension $imagepath]
            set but [ttk::button $w.information.img$dim -image $img -command [list spdAux::SwitchDimAndCreateWindow $dim] ]
            
            grid $w.information.img$dim -column $i -row 0
            #if {[file extension $imagepath] eq ".gif"} {
            #    ::anigif::anigif $imagepath $but
            #    ::anigif::restart $but
            #    W $but
            #}
            incr i
        }
        grid $w.top
        grid $w.top.title_text
        
        grid $w.information
    }
}

proc spdAux::SetSpatialDimmension {ndim} {
    set doc $gid_groups_conds::doc
    set root [$doc documentElement]
    set ::Model::SpatialDimension $ndim
    
    set nd ""
    catch {set nd [$root selectNodes "value\[@n='nDim'\]"]} 
    if {$nd eq ""} {catch {set nd [$root selectNodes "hiddenfield\[@n='nDim'\]"] }}
    
    $nd setAttribute v $::Model::SpatialDimension
}

proc spdAux::SwitchDimAndCreateWindow { ndim } {
    variable TreeVisibility
    SetSpatialDimmension $ndim
    spdAux::DestroyWindow
    
    processIncludes
    parseRoutes
    
    catch {apps::ExecuteOnCurrentXML MultiAppEvent init }
    catch {apps::ExecuteOnCurrentXML CustomTree "" }
    
    if {$TreeVisibility} {
        after 100 [list gid_groups_conds::open_conditions menu ]
        spdAux::PreChargeTree
        spdAux::TryRefreshTree
    }
    
}

proc spdAux::ForceExtremeLoad { } {
    set doc $gid_groups_conds::doc
    set root [$doc documentElement]
    foreach contNode [$root getElementsByTagName "container"] {
        W "Opening [$contNode  @n]"
        $contNode setAttribute tree_state "open"
    }
    gid_groups_conds::actualize_conditions_window
}

proc spdAux::getImagePathDim { dim } {
    set imagepath ""
    set imagepath [apps::getImgPathFrom [apps::getActiveAppId] "$dim.gif"]
    if {[file exists $imagepath]} {return $imagepath}
    set imagepath [apps::getImgPathFrom [apps::getActiveAppId] "$dim.png"]
    if {[file exists $imagepath]} {return $imagepath}
    set imagepath [file nativename [file join $::Model::dir images "$dim.png"] ]
    return $imagepath
}
proc spdAux::DestroyWindow {} {
    variable initwind
    catch {destroy $initwind}
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
    #W "Add $name $route"
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
                foreach u [split [$node getAttribute un] ","] {
                    setRoute $u [$node toXPath]
                }
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
    }
    
}

proc spdAux::GetAppIdFromNode {domNode} {
    set prefix ""
    set prevDomNodeName ""
    while {$prefix eq "" && [$domNode @n] != $prevDomNodeName} {
        set prevDomNode [$domNode @n]
        set domNode [$domNode parent]
        if {[$domNode hasAttribute prefix]} {set prefix [$domNode @prefix]}
    }
    return [$domNode @n]
}

# Dependencies
proc spdAux::insertDependencies { baseNode originUN } {
    set doc $gid_groups_conds::doc
    set root [$doc documentElement]
    catch {
        set originxpath [$baseNode toXPath]
        set insertxpath [getRoute $originUN]
        set insertonnode [$root selectNodes $insertxpath]
        # a lo bestia, cambiar cuando sepamos inyectar la dependencia, abajo esta a medias
        $insertonnode setAttribute "actualize_tree" 1
    }
    ## Aun no soy capaz de insertar y que funcione
    #set ready 1
    #foreach c [$insertonnode getElementsByTagName "dependencies"] {
        #    if {[$c getAttribute "node"] eq $originxpath} {set ready 0; break}
        #}
    #
    #if {$ready} {
        #    set str "<dependencies node='$originxpath' actualize='1'/>"
        #    W $str
        #    W $insertxpath
        #    $insertonnode appendChild [[dom parse $str] documentElement]
        #    W [$insertonnode asXML]
        #}
}
# Dependencies
proc spdAux::insertDependenciesSoft { originxpath relativepath n attn attv} {
    set doc $gid_groups_conds::doc
    set root [$doc documentElement]
    set insertonnode [$root selectNodes $originxpath]
    
    # Aun no soy capaz de insertar y que funcione
    set ready 1
    foreach c [$insertonnode getElementsByTagName "dependencies"] {
        if {[$c getAttribute "node"] eq $originxpath} {set ready 0; break}
    }
    if {$ready} {
        set str "<dependencies n='$n' node='$relativepath' att1='$attn' v1='$attv' actualize='1'/>"
        $insertonnode appendChild [[dom parse $str] documentElement]
    }
}

proc spdAux::CheckSolverEntryState {domNode} {
    set appid [GetAppIdFromNode $domNode]
    set kw [apps::getAppUniqueName $appid SolStrat]
    set nodo [$domNode selectNodes [getRoute $kw]]
    get_domnode_attribute $nodo dict
    set currentSolStrat [get_domnode_attribute $nodo v]
    set mySolStrat [get_domnode_attribute $domNode solstratname]
    return [expr [string compare $currentSolStrat $mySolStrat] == 0]
}

proc spdAux::chk_loads_function_time { domNode load_name } {
    set loads [list [list scalar]]
    lappend loads [list interpolator_func x x T]
    return [join $loads ,]
}

proc spdAux::ViewDoc {} {
    W [[$gid_groups_conds::doc documentElement] asXML]
}

proc spdAux::CheckPartParamValue {node material_name} {
    set doc $gid_groups_conds::doc
    set root [$doc documentElement]
    #W [get_domnode_attribute $node n]
    if {[$node hasAttribute n] || $material_name ne ""} {
        set id [$node getAttribute n]
        set found 0
        set val 0.0
        
        # primero miramos si el material tiene ese campo
        if {$material_name ne ""} {
            set mats_un [apps::getCurrentUniqueName Materials]
            set xp3 [spdAux::getRoute $mats_un]
            append xp3 [format_xpath {/blockdata[@n="material" and @name=%s]/value} $material_name]
            
            foreach valueNode [$root selectNodes $xp3] {
                if {$id eq [$valueNode getAttribute n] } {set val [$valueNode getAttribute v]; set found 1; break}
            }
            #if {$found} {W "mat $material_name value $val"}
        }
        # si no está en el material, miramos en el elemento
        if {!$found} {
            set element_name [get_domnode_attribute [[$node parent] selectNodes "./value\[@n='Element'\]"] v]
            #set claw_name [.gid.central.boundaryconds.gg.data.f0.e1 get]
            set element [Model::getElement $element_name]
            if {$element ne ""} {
                set val [$element getInputDv $id]
                if {$val ne ""} {set found 1}
            }
            #if {$found} {W "element $element_name value $val"}
        }
        # Si no está en el elemento, miramos en la ley constitutiva
        if {!$found} {
            set claw_name [get_domnode_attribute [[$node parent] selectNodes "./value\[@n='ConstitutiveLaw'\]"] v]
            set claw [Model::getConstitutiveLaw $claw_name]
            if {$claw ne ""} {
                set val [$claw getInputDv $id]
                if {$val ne ""} {set found 1}
            }
            #if {$found} {W "claw $claw_name value $val"}
        }
        #if {!$found} {W "Not found $val"}
        if {$val eq ""} {set val 0.0} {return $val}
    }
}

proc spdAux::ConvertAllUniqueNames {oldPrefix newPrefix} {
    variable uniqueNames
    set root [$gid_groups_conds::doc documentElement]
    
    foreach routeName [dict keys $uniqueNames] {
        if {[string first $oldPrefix $routeName] eq 0} {
            set route [getRoute $routeName]
            set newrouteName [string map [list $oldPrefix $newPrefix] $routeName]
            set node [$root selectNodes $route]
            set uns [split [get_domnode_attribute $node un] ","]
            lappend uns $newrouteName
            $node setAttribute un [ListToValues $uns]
        }
    }
    
    spdAux::parseRoutes
}

# Modify the tree: field newValue UniqueName OptionalChild
proc spdAux::SetValueOnTreeItem { field value name {it "" } } {
    set doc $gid_groups_conds::doc
    set root [$doc documentElement]
    #W "$field $value $name $it"
    set node ""
    catch {
        set xp [getRoute $name]
        set node [$root selectNodes $xp]
        if {$it ne ""} {set node [$node find n $it]}
    }
    if {$node ne ""} {
        gid_groups_conds::setAttributes [$node toXPath] [list $field $value]
    } {
        W "$name $it not found - Check GetFromXML.tcl file"
    }
}

proc spdAux::ListToValues {lista} {
    set res ""
    foreach elem $lista {
        append res $elem
        append res ","
    }
    return [string range $res 0 end-1]
}

proc spdAux::injectSolvers {basenode args} {
    
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
            append paramsnodes "<value n='$parname' pn='$pn' state='\[SolverParamState\]' v='$dv' "
            if {$type eq "bool"} {
                append paramsnodes " values='Yes,No' "
            }
            if {$type eq "combo"} {
                set vals [join [$par getValues] ,]
                append paramsnodes " values='$vals' "
                set d [list ]
                foreach v [$par getValues] pv [$par getPValues] {
                    lappend d "$v,$pv"
                }
                set dic [join $d ,]
                append paramsnodes " dict='$dic' "
            }
            
            append paramsnodes "/>"
        }
    }
    set contnode [$basenode parent]
    
    # Get All SolversEntry
    set ses [list ]
    foreach st [::Model::GetSolutionStrategies {*}$args] {
        lappend ses $st [$st getSolversEntries]
    }
    
    # One container per solverEntry 
    foreach {st ss} $ses {
        foreach se $ss {
            set stn [$st getName]
            set n [$se getName]
            set pn [$se getPublicName]
            set help [$se getHelp]
            set appid [GetAppIdFromNode [$basenode parent]]
            set un [apps::getAppUniqueName $appid "$stn$n"]
            set container "<container help='$help' n='$n' pn='$pn' un='$un' state='\[SolverEntryState\]' solstratname='$stn' open_window='0'>"
            set defsolver [lindex [$se getDefaultSolvers] 0]
            append container "<value n='Solver' pn='Solver' v='$defsolver' values='' dict='\[GetSolvers\]' actualize='1' update_proc='UpdateTree'/>"
            #append container "<dependencies node='../value' actualize='1'/>"
            #append container "</value>"
            append container $paramsnodes
            append container "</container>"
            $contnode appendXML $container
        }
    }
    $basenode delete
}

proc spdAux::injectSolStratParams {basenode args} {
    set contnode [$basenode parent]
    set params [::Model::GetSolStratParams {*}$args]
    foreach {parname par} $params {
        #W "$parname [$contnode find n $parname]"
        if {[$contnode find n $parname] eq ""} {
            set pn [$par getPublicName]
            set type [$par getType]
            set dv [$par getDv]
            if {$type eq "bool"} {set dv [GetBooleanForTree $dv]}
            set helptext [$par getHelp]
            set actualize [$par getActualize]
            set node "<value n='$parname' pn='$pn' state='\[SolStratParamState\]' v='$dv' help='$helptext' "
            
            if {$actualize} {
                append node "actualize_tree='1'"
            }
            
            if {$type eq "bool"} {
                
                append node " values='Yes,No' "
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
                append node " values='$vs' dict='$pv' "
            } 
            append node "/>"
            catch {
                $contnode appendXML $node
                set orig [$contnode lastChild]
                set new [$orig cloneNode]
                $orig delete
                $contnode insertBefore $new $basenode
            }
        }
    }
    
    set params [::Model::GetSchemesParams {*}$args]
    
    foreach {parname par} $params {
        #W "$parname [$contnode find n $parname]"
        if {[$contnode find n $parname] eq ""} {
            set pn [$par getPublicName]
            set type [$par getType]
            set dv [$par getDv]
            if {$type eq "bool"} {set dv [GetBooleanForTree $dv]}
            set helptext [$par getHelp]
            set node "<value n='$parname' pn='$pn' state='\[SchemeParamState\]' v='$dv' help='$helptext' "
            if {$type eq "bool"} {
                append node " values='Yes,No' "
            } 
            append node "/>"
            catch {$contnode appendXML $node}
        }
    }
    $basenode delete
}



proc spdAux::injectNodalConditions { basenode args} {
    if {$args eq "{}"} {
        set nodal_conditions [::Model::getAllNodalConditions]
    } {
        set nodal_conditions [::Model::GetNodalConditions {*}$args]
    }
    spdAux::_injectCondsToTree $basenode $nodal_conditions "nodal"
    $basenode delete
}

proc spdAux::injectConditions { basenode args} {
    set conditions [::Model::GetConditions {*}$args]
    spdAux::_injectCondsToTree $basenode $conditions
    $basenode delete
}

proc spdAux::_injectCondsToTree {basenode cond_list {cond_type "normal"} } {
    set conds [$basenode parent]
    set AppUsesIntervals [apps::ExecuteOnApp [GetAppIdFromNode $conds] GetAttribute UseIntervals]
    if {$AppUsesIntervals eq ""} {set AppUsesIntervals 0}
    
    foreach cnd $cond_list {
        set n [$cnd getName]
        set pn [$cnd getPublicName]
        set help [$cnd getHelp]
        set etype ""
        if {$cond_type eq "nodal"} {
            set etype [$cnd getOv]
        } else {
            set etype [join [string tolower [$cnd getAttribute ElementType]] ,]
        }
        if {$etype eq ""} {set etype "point,line,surface,volume"}
        set units [$cnd getAttribute "units"]
        set um [$cnd getAttribute "unit_magnitude"]
        set process [::Model::GetProcess [$cnd getProcessName]]
        set check [$process getAttribute "check"]
        set state "ConditionState"
        if {$cond_type eq "nodal"} {
            set state [$cnd getAttribute state]
            if {$state eq ""} {set state "CheckNodalConditionState"}
        }
        set node "<condition n='$n' pn='$pn' ov='$etype' ovm='' icon='shells16' help='$help' state='\[$state\]' update_proc='$check'>"
        #foreach processinput [$process getInputs] {lappend inputs $processinput}
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
            if {$type eq "vector"} {
                set vector_type [$in getAttribute "vectorType"]
                lassign [split $v ","] v1 v2 v3
                if {$vector_type eq "bool"} {
                    append node "
                        <value n='${inName}X' wn='[concat $n "_X"]' pn='X ${pn}' v='$v1' values='1,0' help=''/>
                        <value n='${inName}Y' wn='[concat $n "_Y"]' pn='Y ${pn}' v='$v2' values='1,0' help=''/>
                        <value n='${inName}Z' wn='[concat $n "_Z"]' pn='Z ${pn}' v='$v3' values='1,0' help='' state='\[CheckDimension 3D\]'/>"
                } {
                    append node "
                    <value n='${inName}X' wn='[concat $n "_X"]' pn='${pn} X' v='$v1' help='$help' />
                    <value n='${inName}Y' wn='[concat $n "_Y"]' pn='${pn} Y' v='$v2' help='$help' />
                    <value n='${inName}Z' wn='[concat $n "_Z"]' pn='${pn} Z' v='$v3' help='$help'  state='\[CheckDimension 3D\]'/>
                    "
                }
                
            } elseif { $type eq "combo" } {
                set values [join [$in getValues] ","]
                append node "<value n='$inName' pn='$pn' v='$v' values='$values' state='$state' help='$help'/>"
            } elseif { $type eq "bool" } {
                set values "1,0"
                append node "<value n='$inName' pn='$pn' v='$v' values='$values'  help='$help' state='$state'/>"
            } elseif { $type eq "file" || $type eq "tablefile" } {
                append node "<value n='$inName' pn='$pn' v='$v' values='\[GetFilesValues\]' update_proc='AddFile' help='$help' state='$state' type='$type'/>"
            } else {
                append node "<value n='$inName' pn='$pn' v='$v'  units='$units'  unit_magnitude='$um'  help='$help'/>"
                if {[$in getAttribute "function"] eq "1"} {
                    set fname "function_$inName"
                    set nodev "../value\[@n='$inName'\]"
                    set nodef "../value\[@n='$fname'\]"
                    append node "<value n='$fname' pn='Function' v='' help='$help'/>"
                    append node "<value n='ByFunction' pn='By function' v='No' values='Yes,No'  actualize_tree='1'>
                                    <dependencies value='No' node=\""
                    append node $nodev
                    append node "\" att1='state' v1='normal'/>
                                    <dependencies value='Yes'  node=\""
                    append node $nodev
                    append node "\" att1='state' v1='hidden'/>
                                    <dependencies value='No' node=\""
                    append node $nodef
                    append node "\" att1='state' v1='hidden'/>
                                    <dependencies value='Yes'  node=\""
                    append node $nodef
                    append node "\" att1='state' v1='normal'/>
                                  </value>"
                }
                #append node "<value n='$inName' pn='$pn' v='$v'   help='$help'/>"
            }
        }
        
        set CondUsesIntervals [$cnd getAttribute "Interval"]
        if {$AppUsesIntervals && $CondUsesIntervals ne "False"} {
            append node "<value n='Interval' pn='Time interval' v='$CondUsesIntervals' values='\[getIntervals\]'  help='$help'/>"
        }
        append node "</condition>"
        $conds appendXML $node
    }
    
}
proc spdAux::injectElementInputs { basenode args} {
    set parts [$basenode parent]
    set inputs [::Model::GetAllElemInputs]
    foreach inName [dict keys $inputs] {
        set in [dict get $inputs $inName] 
        set pn [$in getPublicName]
        set units [$in getUnits]
        set um [$in getUnitMagnitude]
        set help [$in getHelp] 
        set v [$in getDv]
        #set node "<value n='$inName' pn='$pn' state='\[PartParamState\]' v='-' units='$units' unit_magnitude='$um' help='$help' />"
        set node "<value n='$inName' pn='$pn' state='\[PartParamState\]' v='$v' help='$help' />"
        catch {
            $parts appendXML $node
            set orig [$parts lastChild]
            set new [$orig cloneNode]
            $orig delete
            $parts insertBefore $new $basenode
        }
        
        #set originxpath "[$parts toXPath]/value\[@n='Material'\]"
        #set relativexpath "../value\[@n='$inName'\]"
        #spdAux::insertDependenciesSoft $originxpath $relativexpath $inName v "\[PartParamValue\]"
    }
    $basenode delete
}

proc spdAux::injectConstitutiveLawInputs { basenode  args} {
    set parts [$basenode parent]
    set inputs [::Model::GetAllCLInputs]
    foreach inName [dict keys $inputs] {
        if {[$parts find n $inName] eq ""} {
            set in [dict get $inputs $inName] 
            set pn [$in getPublicName]
            set units [$in getUnits]
            set um [$in getUnitMagnitude]
            set help [$in getHelp]
            set v [$in getDv]
            #set node "<value n='$inName' pn='$pn' state='\[PartParamState\]' v='$v' units='$units' unit_magnitude='$um' help='$help' />"
            set node "<value n='$inName' pn='$pn' state='\[PartParamState\]' v='$v' help='$help' />"
            catch {
                $parts appendXML $node
                set orig [$parts lastChild]
                set new [$orig cloneNode]
                $orig delete
                $parts insertBefore $new $basenode
            }
            
            #set originxpath "[$parts toXPath]/value\[@n='Material'\]"
            #set relativexpath "../value\[@n='$inName'\]"
            #spdAux::insertDependenciesSoft $originxpath $relativexpath $inName v "\[PartParamValue\]"
        }
    }
    $basenode delete
}

proc spdAux::injectElementOutputs { basenode args} {
    set parts [$basenode parent]
    set outputs [::Model::GetAllElemOutputs {*}$args]
    foreach in [dict keys $outputs] {
        set pn [[dict get $outputs $in] getPublicName]
        set v [GetBooleanForTree [[dict get $outputs $in] getDv]]
        
        set node "<value n='$in' pn='$pn' state='\[ElementOutputState\]' v='$v' values='Yes,No' />"
        catch {
                $parts appendXML $node
                set orig [$parts lastChild]
                set new [$orig cloneNode]
                $orig delete
                $parts insertBefore $new $basenode
            }
    }
    $basenode delete
}

proc spdAux::injectNodalConditionsOutputs { basenode args} {
    set parts [$basenode parent]
    
    if {$args eq "{}"} {
        set nodal_conditions [::Model::getAllNodalConditions]
    } {
        set nodal_conditions [::Model::GetNodalConditions {*}$args]
    }
    foreach nc $nodal_conditions {
        set n [$nc getName]
        set pn [$nc getPublicName]
        set v [$nc getAttribute v]
        if {$v eq ""} {set v "Yes"}
        set node "<value n='$n' pn='$pn' v='$v' values='Yes,No' state='\[CheckNodalConditionState\]'/>"
        catch {$parts appendXML $node}
        foreach {n1 output} [$nc getOutputs] {
            set nout [$output getName]
            set pn [$output getPublicName]
            set v [$output getAttribute v]
            if {$v eq ""} {set v "Yes"}
            set node "<value n='$nout' pn='$pn' v='$v' values='Yes,No' state='\[CheckNodalConditionOutputState $n\]'/>"
            catch {$parts appendXML $node}
        }
    }
    $basenode delete
}

proc spdAux::GetBooleanForTree {raw} {
    set goodList [list "Yes" "1" "yes" "ok" "YES" "Ok" "True" "TRUE" "true"]
    if {$raw in $goodList} {return "Yes" } {return "No"}
}

proc spdAux::injectConstitutiveLawOutputs { basenode  args} {
    set parts [$basenode parent]
    set outputs [::Model::GetAllCLOutputs {*}$args]
    foreach in [dict keys $outputs] {
        if {[$parts find n $in] eq ""} {
            set pn [[dict get $outputs $in] getPublicName]
            set v [GetBooleanForTree [[dict get $outputs $in] getDv]]
            set node "<value n='$in' pn='$pn' state='\[ConstLawOutputState\]' v='$v' values='Yes,No' />"
            catch {
                $parts appendXML $node
                set orig [$parts lastChild]
                set new [$orig cloneNode]
                $orig delete
                $parts insertBefore $new $basenode
            }
        }
    }
    $basenode delete
}

proc spdAux::injectProcs { basenode  args} {
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
            set pastnode [[$basenode parent] selectNodes "./proc\[@n='$procn'\]"]
            if {$pastnode ne ""} {gid_groups_conds::delete [$pastnode toXPath]}
            
            [$basenode parent] appendChild $in
        }
        $basenode delete
    }
}

proc spdAux::CheckConstLawOutputState {outnode} {
    set doc $gid_groups_conds::doc
    set root [$doc documentElement]
    
    set nodeApp [GetAppIdFromNode $outnode]
    set parts_un [apps::getAppUniqueName $nodeApp Parts]
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
    
    set nodeApp [GetAppIdFromNode $outnode]
    set parts_un [apps::getAppUniqueName $nodeApp Parts]
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
    set nodeApp [GetAppIdFromNode $outnode]
    
    set sol_stratUN [apps::getAppUniqueName $nodeApp SolStrat]
    
    if {[get_domnode_attribute [$root selectNodes [spdAux::getRoute $sol_stratUN]] v] eq ""} {
        get_domnode_attribute [$root selectNodes [spdAux::getRoute $sol_stratUN]] dict
    }
    set SolStrat [::write::getValue $sol_stratUN]
    
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
    set nodeApp [GetAppIdFromNode $outnode]
    
    set sol_stratUN [apps::getAppUniqueName $nodeApp SolStrat]
    set schemeUN [apps::getAppUniqueName $nodeApp Scheme]
    
    if {[get_domnode_attribute [$root selectNodes [spdAux::getRoute $sol_stratUN]] v] eq ""} {
        get_domnode_attribute [$root selectNodes [spdAux::getRoute $sol_stratUN]] dict
    }
    if {[get_domnode_attribute [$root selectNodes [spdAux::getRoute $schemeUN]] v] eq ""} {
         get_domnode_attribute [$root selectNodes [spdAux::getRoute $schemeUN]] dict
    }
    set SolStrat [::write::getValue $sol_stratUN]
    set Scheme [write::getValue $schemeUN]
    
    set paramName [$outnode @n]
    return [::Model::CheckSchemeInputState $SolStrat $Scheme $paramName]
}

proc spdAux::getIntervals { {un ""} } {
    set doc $gid_groups_conds::doc
    set root [$doc documentElement]
    if {$un eq ""} {
        set un "Intervals"
    } 
    set xp1 "[spdAux::getRoute $un]/blockdata\[@n='Interval'\]"
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



############# procs #################

proc spdAux::ProcGetElements { domNode args } {
    set nodeApp [GetAppIdFromNode $domNode]
    set sol_stratUN [apps::getAppUniqueName $nodeApp SolStrat]
    set schemeUN [apps::getAppUniqueName $nodeApp Scheme]
    if {[get_domnode_attribute [$domNode selectNodes [spdAux::getRoute $sol_stratUN]] v] eq ""} {
        get_domnode_attribute [$domNode selectNodes [spdAux::getRoute $sol_stratUN]] dict
    }
    if {[get_domnode_attribute [$domNode selectNodes [spdAux::getRoute $schemeUN]] v] eq ""} {
         get_domnode_attribute [$domNode selectNodes [spdAux::getRoute $schemeUN]] dict
    }
    
    #W "solStrat $sol_stratUN sch $schemeUN"
    set solStratName [::write::getValue $sol_stratUN]
    set schemeName [write::getValue $schemeUN]
    #W "$solStratName $schemeName"
    #W "************************************************************************"
    #W "$nodeApp $solStratName $schemeName"
    set elems [::Model::GetAvailableElements $solStratName $schemeName]
    #W "************************************************************************"
    set names [list ]
    set pnames [list ]
    foreach elem $elems {
        if {[$elem cumple {*}$args]} {
            lappend names [$elem getName]
            lappend pnames [$elem getName] 
            lappend pnames [$elem getPublicName]
        }
    }
    set diction [join $pnames ","]
    set values [join $names ","]
    #W "[get_domnode_attribute $domNode v] $names"
    $domNode setAttribute values $values
    if {[get_domnode_attribute $domNode v] eq ""} {$domNode setAttribute v [lindex $names 0]}
    if {[get_domnode_attribute $domNode v] ni $names} {$domNode setAttribute v [lindex $names 0]; spdAux::RequestRefresh}
    #spdAux::RequestRefresh
    return $diction
}

proc spdAux::ProcGetSolutionStrategies {domNode args} {
    set names [list ]
    set pnames [list ]
    #W $args
    set Sols [::Model::GetSolutionStrategies {*}$args]
    #W $Sols
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

proc spdAux::ProcGetSchemes {domNode args} {
    set nodeApp [GetAppIdFromNode $domNode]
    #W $nodeApp
    set sol_stratUN [apps::getAppUniqueName $nodeApp SolStrat]
    #W "GS sol $sol_stratUN"
    if {[get_domnode_attribute [$domNode selectNodes [spdAux::getRoute $sol_stratUN]] v] eq ""} {
        #W "entra"
        get_domnode_attribute [$domNode selectNodes [spdAux::getRoute $sol_stratUN]] dict
    }
    set solStratName [::write::getValue $sol_stratUN]
    #W "Unique name: $sol_stratUN - Nombre $solStratName"
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
    
    $domNode setAttribute values [join $names ","]
    if {[get_domnode_attribute $domNode v] eq ""} {$domNode setAttribute v [lindex $names 0]}
    if {[get_domnode_attribute $domNode v] ni $names} {$domNode setAttribute v [lindex $names 0]}
    #spdAux::RequestRefresh
    return [join $pnames ","]
}

proc spdAux::SetNoneValue {domNode} {
    $domNode setAttribute v "None"
    $domNode setAttribute values "None"
    #spdAux::RequestRefresh
    return "None,None"
}

#This should go to values
proc spdAux::ProcGetConstitutiveLaws { domNode args } {
    set Elementname [$domNode selectNodes {string(../value[@n='Element']/@v)}]
    set Claws [::Model::GetAvailableConstitutiveLaws $Elementname]
    #W "Const Laws que han pasado la criba: $Claws"
    if {[llength $Claws] == 0} { return [SetNoneValue $domNode] }
    set names [list ]
    foreach cl $Claws {
        lappend names [$cl getName]
    }
    set values [join $names ","]
    if {[get_domnode_attribute $domNode v] eq "" || [get_domnode_attribute $domNode v] ni $names} {$domNode setAttribute v [lindex $names 0]}
    #spdAux::RequestRefresh
    
    return $values
}
#This should go to dict
proc spdAux::ProcGetAllConstitutiveLaws { domNode args } {
    set Claws [Model::GetConstitutiveLaws]
    if {[llength $Claws] == 0} { return [SetNoneValue $domNode] }
    set pnames [list ]
    foreach cl $Claws {
        lappend pnames [$cl getName] 
        lappend pnames [$cl getPublicName]
    }
    set diction [join $pnames ","]
    #spdAux::RequestRefresh

    return $diction
}
proc spdAux::ProcGetSolvers { domNode args } {
    
    set solStrat [get_domnode_attribute [$domNode parent] solstratname]
    set solverEntryId [get_domnode_attribute [$domNode parent] n]
    
    set solvers [Model::GetAvailableSolvers $solStrat $solverEntryId]
    
    set names [list ]
    set pnames [list ]
    foreach slvr $solvers {
        lappend names [$slvr getName]
        lappend pnames [$slvr getName] 
        lappend pnames [$slvr getPublicName]
    }
    $domNode setAttribute values [join $names ","]
    if {[get_domnode_attribute $domNode v] eq ""} {$domNode setAttribute v [lindex $names 0]}
    return [join $pnames ","]
    
}
proc spdAux::ProcConditionState { domNode args } {
    
    set resp [::Model::CheckConditionState $domNode]
    if {$resp} {return "normal"} else {return "hidden"}
}

proc spdAux::ProcCheckNodalConditionState { domNode args } {
    
    set nodeApp [GetAppIdFromNode $domNode]
    set parts_un [apps::getAppUniqueName $nodeApp Parts]
    #W $parts_un
    if {[spdAux::getRoute $parts_un] ne ""} {
        set conditionId [$domNode @n]
        set elems [$domNode selectNodes "[spdAux::getRoute $parts_un]/group/value\[@n='Element'\]"]
        set elemnames [list ]
        foreach elem $elems { lappend elemnames [$elem @v]}
        set elemnames [lsort -unique $elemnames]
        if {$elemnames eq ""} {return "hidden"}
        if {[::Model::CheckElementsNodalCondition $conditionId $elemnames]} {return "normal"} else {return "hidden"}
    } {return "normal"}
}
proc spdAux::ProcCheckNodalConditionOutputState { domNode args } {
    
    set nodeApp [GetAppIdFromNode $domNode]
    set NC_un [apps::getAppUniqueName $nodeApp NodalConditions]
    if {[spdAux::getRoute $NC_un] ne ""} {
        set ncs [$domNode selectNodes "[spdAux::getRoute $NC_un]/condition/group"]
        set ncslist [list ]
        foreach nc $ncs { lappend ncslist [[$nc parent] @n]}
        set ncslist [lsort -unique $ncslist]
        set conditionId [lindex $args 0]
        if {$conditionId ni $ncslist} {return "hidden"} {return "normal"}
        set outputId [$domNode @n]
        if {[::Model::CheckNodalConditionOutputState $conditionId $outputId]} {return "normal"} else {return "hidden"}
    } {return "normal"}
}
proc spdAux::ProcRefreshTree { domNode args } {
    spdAux::RequestRefresh
}

proc spdAux::ProccheckStateByUniqueName { domNode args } {
    
    set total 0
    foreach {un val} {*}$args {
        catch {
            set xpath [spdAux::getRoute $un]
            spdAux::insertDependencies $domNode $un
            set node [$domNode selectNodes $xpath]
            set realval [get_domnode_attribute $node v]
            if {$realval eq ""} {W "Warning: Check unique name $un"}
            if {[lsearch $val $realval] != -1} {
                set total 1
                break
            }
        }
    }
    if {$total} {return "normal"} else {return "hidden"}
}
proc spdAux::ProcSolverParamState { domNode args } {
    
    
    set id [$domNode getAttribute n]
    set nodesolver [[$domNode parent] selectNodes "./value\[@n='Solver'\]"]
    get_domnode_attribute $nodesolver values
    set solverid [get_domnode_attribute $nodesolver v]
    
    if {$solverid eq ""} {set resp 0} {
        set resp [::Model::getSolverParamState $solverid $id]
    }
    
    #spdAux::RequestRefresh
    if {$resp} {return "normal"} else {return "hidden"}
}
proc spdAux::ProcPartParamValue { domNode args } {
    
    set nodename [get_domnode_attribute $domNode n]
    set matname [get_domnode_attribute $domNode v]
    set node [[$domNode parent] selectNode "../value\[@n='$nodename'\]" ]
    set nodevalue [$node @v]
    return [spdAux::CheckPartParamValue $node $matname]
}
proc spdAux::ProcPartParamState { domNode args } {
    set resp [::Model::CheckElemParamState $domNode]
    if {$resp eq "0"} {
        set id [$domNode getAttribute n]
        set constLaw [get_domnode_attribute [[$domNode parent] selectNodes "./value\[@n='ConstitutiveLaw'\]"] v]
        if {$constLaw eq ""} {return hidden}
        set resp [Model::CheckConstLawParamState $constLaw $id]
    }
    
    #W "Calculando estado de [$domNode @pn] : $resp"
    if {$resp} {return "normal"} else {return "hidden"}
}
proc spdAux::ProcSolverEntryState { domNode args } {
    
    set resp [spdAux::CheckSolverEntryState $domNode]
    if {$resp} {return "normal"} else {return "hidden"}
}
proc spdAux::ProcCheckDimension { domNode args } {
    
    set checkdim [lindex $args 0]
    
    if {$checkdim eq $::Model::SpatialDimension} {return "normal"} else {return "hidden"}
}
proc spdAux::ProcgetStateFromXPathValue2 { domNode args } {
    set args {*}$args
    set arglist [split $args " "]
    set xpath {*}[lindex $arglist 0]
    set checkvalue [lindex $arglist 1]
    set pst [$domNode selectNodes $xpath]
    #W "xpath $xpath checkvalue $checkvalue pst $pst"
    if {$pst == $checkvalue} { return "normal"} else {return "hidden"}
}

proc spdAux::ProcgetStateFromXPathValue { domNode args } {
    set args {*}$args
    set arglist [split $args " "]
    set xpath {*}[lindex $arglist 0]
    set checkvalue [split [lindex $arglist 1] ","]
    set pst [$domNode selectNodes $xpath]
    #W "xpath $xpath checkvalue $checkvalue pst $pst"
    if {$pst in $checkvalue} { return "normal"} else {return "hidden"}
}
proc spdAux::ProcSolStratParamState { domNode args } {
    
    set resp [::spdAux::SolStratParamState $domNode]
    if {$resp} {return "normal"} else {return "hidden"}
}
proc spdAux::ProcSchemeParamState { domNode args } {
    
    set resp [::spdAux::SchemeParamState $domNode]
    if {$resp} {return "normal"} else {return "hidden"}
}  
proc spdAux::ProcConstLawOutputState { domNode args } {
    
    set resp [::spdAux::CheckConstLawOutputState $domNode]
    if {$resp} {return "normal"} else {return "hidden"}
}
proc spdAux::ProcElementOutputState { domNode args } {
    
    set resp [::spdAux::CheckElementOutputState $domNode]
    if {$resp} {return "normal"} else {return "hidden"}
}

proc spdAux::ProcActiveIfAnyPartState { domNode args } {
    
    set parts ""
    set nodeApp [GetAppIdFromNode $domNode]
    set parts_un [apps::getAppUniqueName $nodeApp Parts]
    catch {
        set parts [$domNode selectNodes "[spdAux::getRoute $parts_un]/group"]
    }
    if {$parts ne ""} {return "normal"} else {return "hidden"}
}
proc spdAux::ProcActiveIfRestartAvailable { domNode args } {
    
    set active [apps::ExecuteOnApp [GetAppIdFromNode $domNode] GetAttribute UseRestart]
    if {$active ne "" && $active} {return "normal"} else {return "hidden"}
}

proc spdAux::ProcDisableIfUniqueName { domNode args } {
    set total 1
    foreach {un val} {*}$args {
        set xpath [spdAux::getRoute $un]
        spdAux::insertDependencies $domNode $un
        set node [$domNode selectNodes $xpath]
        if {$node eq ""} {
            set total 0
            W "Warning: state of [$domNode @n]"
        } else {
            set realval [get_domnode_attribute $node v]
            if {$realval eq ""} {W "Warning: Check unique name $un"}
            if {[lsearch $val $realval] == -1} {
                set total 0
                break
            } 
        }
    }
    if {!$total} {return "normal"} else {return "disabled"}
}
proc spdAux::ProcCheckGeometry { domNode args } {
    
    set level [lindex $args 0]
    #W $level
    if {$level eq 1} {
        if {$::Model::SpatialDimension eq "3D"} {return volume} {return surface}
    }
    if {$level eq 2} {
        if {$::Model::SpatialDimension eq "3D"} {return surface} {return line}
    }
}
proc spdAux::ProcDirectorVectorNonZero { domNode args } {
    
    set kw [lindex $args 0]
    set update 0
    foreach condgroupnode [$domNode getElementsByTagName group] {
        set valid 0
        foreach dirnode [$condgroupnode getElementsByTagName value] {
            if {[string first $kw [get_domnode_attribute $dirnode n]] eq 0 } {
                if { [get_domnode_attribute $dirnode v] != 0 } {set valid 1; break}
            }
        }
        if {!$valid} {
            $domNode removeChild $condgroupnode
            set update 1
        }
    }
    if {$update} {
        W "Director vector can't be null"
        gid_groups_conds::actualize_conditions_window
    }
}
proc spdAux::ProcShowInMode { domNode args } {
    set kw [lindex $args 0]
    if {$kw ni [list "Release" "Developer"]} {return "hidden"}
    if {$::Kratos::kratos_private(DevMode) eq "dev"} {
        if {$kw eq "Developer"} {return "normal"} {return "hidden"}
    }
    if {$::Kratos::kratos_private(DevMode) eq "release"} {
        if {$kw eq "Developer"} {return "hidden"} {return "normal"}
    }
}


proc spdAux::LoadModelFiles { } {
    customlib::UpdateDocument
    foreach elem [[customlib::GetBaseRoot] getElementsByTagName "file"] {
        FileSelector::AddFile [$elem @n]
    }
}
proc spdAux::SaveModelFile { fileid } {
    customlib::UpdateDocument
    FileSelector::AddFile $fileid
    gid_groups_conds::addF {container[@n='files']} file [list n ${fileid}]
}

proc spdAux::AddFile { domNode } {
    FileSelector::InitWindow "spdAux::UpdateFileField" $domNode
}
proc spdAux::UpdateFileField { fileid domNode} {
    if {$fileid ne ""} {
        $domNode setAttribute v $fileid
        spdAux::SaveModelFile $fileid
        RequestRefresh 
    }
}
proc spdAux::ProcGetFilesValues { } {
    lappend listilla "- No file"
    lappend listilla {*}[FileSelector::GetAllFiles]
    lappend listilla "- Add new file"
    return [join $listilla ","]
}

proc spdAux::ProcGetIntervals {domNode args} {
    set lista [::spdAux::getIntervals]
    if {$lista eq ""} {$domNode setAttribute state "hidden"; spdAux::RequestRefresh}
    if {[$domNode @v] eq "" || [$domNode @v] ni $lista} {
        $domNode setAttribute v [lindex $lista 0]
    }
    set res [spdAux::ListToValues $lista]
    return $res
}

proc spdAux::PreChargeTree { } {
    return ""
    set doc $gid_groups_conds::doc
    set root [$doc documentElement]
    
    foreach field [list value condition container] {
        foreach cndNode [$root getElementsByTagName $field] {
            set a [get_domnode_attribute $cndNode dict]
            set a [get_domnode_attribute $cndNode values]
            set a [get_domnode_attribute $cndNode v]
            #W [get_domnode_attribute $cndNode n]
        }
    }
}
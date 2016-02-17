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
        
        apps::ExecuteOnCurrent $elem $func
    }
}

proc spdAux::processIncludesRecurse { nf basedir} {
    set xml [tDOM::xmlReadFile $nf]
    set newnode [dom parse [string trim $xml]]
    set xmlNode [$newnode documentElement]
    foreach child [$xmlNode childNodes] {
        if {[$child nodeName] eq "include"} {
            set active [$child getAttribute "active"]
            if {$active} {
                set path [$child getAttribute "path"]
                set f [file join $basedir $path]
                set processednode [spdAux::processIncludesRecurse $f $basedir]
                $xmlNode insertBefore $processednode $child
                $child delete
            }
        } else {
            foreach elem [$child getElementsByTagName "include"] {
                set active [$elem getAttribute "active"]
                if {$active} {
                    set path [$elem getAttribute "path"]
                    set f [file join $basedir $path]
                    set processednode [processIncludesRecurse $f $basedir]
                    $child insertBefore $processednode $elem
                    $elem delete
                }
            }
        }
    }
    return $xmlNode
}

proc spdAux::processAppIncludes { root } {
    foreach elem [$root getElementsByTagName "appLink"] {
        set active [$elem getAttribute "active"]
        set appid [$elem getAttribute "appid"]
        set pn [$elem getAttribute "pn"]
        apps::NewApp $appid $pn
        if {$active} {
            set dir $::Kratos::kratos_private(Path)
            set f [file join $dir apps $appid xml Main.spd]
            set processedAppnode [spdAux::processIncludesRecurse $f $dir]
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
    gid_groups_conds::actualize_conditions_window
    parseRoutes
}

proc spdAux::CreateWindow {dir} {
    variable initwind
    set doc $gid_groups_conds::doc
    set root [$doc documentElement]
    
    set activeapp [ [$root selectNodes "hiddenfield\[@n='activeapp'\]"] getAttribute v]
        
    if { $activeapp ne "" } {activeApp $activeapp; catch {destroy $initwind}; return ""}
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

    InitWindow $w [_ "Kratos Multiphysics"] Kratos "" \
        "" 1
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
    set ::Model::SpatialDimension "3D"
    set sdcmb [ttk::combobox $frsd.sdcmb -textvariable ::Model::SpatialDimension -values $vs -width 4 -state "readonly"]
    
    
    grid $w.top
    grid $w.top.title_text
    
    grid $w.information
    grid $sd -row 0 -column 0 -padx 20
    #-sticky w -padx 20
    grid $sdcmb  -row 0 -column 1 -padx 20
    #-sticky e 
    
    #grid rowconfigure $frsd 0 -weight 1
    #grid columnconfigure $frsd 0 -weight 1
    #grid columnconfigure $frsd 1 -weight 1
    grid $frsd -sticky we
    #set pp [ winfo parent $frsd]
    #grid rowconfigure $pp 2 -weight 1
    #grid columnconfigure $pp 0 -weight 1
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
                #set stateparent ""
                #catch {set stateparent [[$node parent] getAttribute "tree_state"]}
                #if {$stateparent eq ""} {set stateparent "open"}
                #W "Padre -> [[$node parent] getAttribute n]"
                #gid_groups_conds::setAttributesN [$node parent] [list "tree_state" "open"]
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


# External File
proc spdAux::getElements {app fil tag val pn} {
    variable currentexternalfile
    
    if {$currentexternalfile ne "$app//$fil"} {
        set dir $::Kratos::kratos_private(Path)
        [xmlGet new ] init "$dir\\apps\\$app\\xml\\$fil.xml"
        set currentexternalfile "$app//$fil"
    }
    return [[xmlGet new] getElementsPretty $tag $val $pn]
}

proc spdAux::CheckSolverEntryState {domNode} {
        set kw "SMSolStrat"
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
    
    set xp3 [spdAux::getRoute "SMMaterials"]
    append xp3 [format_xpath {/blockdata[@n="material" and @name=%s]/value} $material_name]

    foreach valueNode [$root selectNodes $xp3] {
        if {$id eq [$valueNode getAttribute n] } {set val [$valueNode getAttribute v]}
    }     
    
    # W "mat: $material_name prop $id val $val"
    
    
    return $val
}


proc spdAux::chk_loads_function_time { domNode load_name } {
    set loads [list [list scalar]]
    lappend loads [list interpolator_func x x T]
    return [join $loads ,]
}



proc spdAux::CheckElemParamValue {node} {
    set doc $gid_groups_conds::doc
    set root [$doc documentElement]
    
    set id [$node getAttribute n]
    set val 0.0
    set material_name [get_domnode_attribute [[$node parent] selectNodes "./value\[@n='Material'\]"] v]
    set material_name [.gid.central.boundaryconds.gg.data.f0.e2 get]
    
    set xp3 [spdAux::getRoute "SMMaterials"]
    append xp3 [format_xpath {/blockdata[@n="material" and @name=%s]/value} $material_name]

    foreach valueNode [$root selectNodes $xp3] {
        if {$id eq [$valueNode getAttribute n] } {set val [$valueNode getAttribute v]}
    }
    
    return $val
}

proc spdAux::ExecuteOnCurrentApp {cmd args} {
    return [apps::ExecuteOnCurrent $args $cmd]
}

proc spdAux::ListToValues {lista} {
    set res ""
    foreach elem $lista {
        append res $elem
        append res ","
    }
    return [string range $res 0 end-1]
}

spdAux::Init

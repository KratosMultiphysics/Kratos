#####################################################################################
#
#  NAME: kmpropstree.tcl
#
#  PURPOSE: Tree widget used in the kratos main window to manage model/material/curve
#           properties
#
#  QUANTECH ATZ-DEVELOPMENT DEPARTMENT
#
#  AUTHORS : G. Socorro and L. Calvo
#
#  CREATED AT: 29/03/2012
#
#  HISTORY:
# 
#   0.5- 04/07/13- A.Melendo, popup to don't forget to click ok, cancel on bottom frame
#   0.4- 25/06/13- A. Melendo, new List and Draw procs
#   0.3- 17/05/13-G. Socorro, reformat the source code using Emacs tabulations 
#   0.2- 20/09/12-J.Garate, Adaptation for New Kratos Interface Version, including Curves support
#   0.1- 29/03/2012 G. Socorro, create a base source code from the kmprops.tcl script
#
######################################################################################
#                      Procedures that belong to this file
###############################################################################
#         Name                      |        Functionality
#------------------------------------------------------------------------------
# 1. CreateTreeProperties           | Create the treectrl properties 
# 2. RefreshTree                    | Refresh the treectrl tree properties
# 3. CreateTreeMaterials            | Create the treectrl material tree path
# 4. DoubleClickTree                | Manage the double click in the tree
# 5. ClickTree                      | Manage the a click in the tree
# 6.

# if enabled, fulltktree is used and popup helps will show up on the tree leaves
# if not, treectrl is used, with no help baloons
set ::KMProps::UseFullTkTree 1

# ERROR: Using fulltktree gives an error when hidding the tree again!

if { $::KMProps::UseFullTkTree} {
    package require fulltktree
}

proc ::KMProps::CreateTreeProperties {w} {
    # ABSTRACT: Create the treectrl properties 
    # Arguments
    # w => Frame path
    # Return
    # T => The tree path

    variable SystemHighlight
    variable SystemHighlightText
    
    # Scrollbars
    set vsb $w.vsb1
    set hsb $w.hsb1
    # Create the treectrl and set the scrollbars
    if { $::KMProps::UseFullTkTree} {
        if { [ catch {set T [fulltktree $w.t]}]} {
            set T [fulltktree $w.t]
        }
    } else {
        set T [treectrl $w.t -xscrollcommand [list $hsb set] -yscrollcommand [list $vsb set]]
        ttk::scrollbar $vsb -orient vertical -command [list $T yview]
        ttk::scrollbar $hsb -orient horizontal -command [list $T xview]
    }
    
    # Set the height
    set height [font metrics [$T cget -font] -linespace]
    if {$height < 18} {
        set height 18
    }
    
    # Configure the treectrl
    $T configure -indent 15 -itemheight $height -selectmode browse \
        -showroot 0 -showrootbutton 0 -showbuttons 1 -showlines 1 \
        -highlightthickness 0 -borderwidth 0 -height 400 \
        -xscrollincrement 20 -yscrollincrement 20
    
    
    $T column create -text " " -tags C0 -weight 0
    
    # Configure the column weight and arrow
    $T column configure C0 -weight 1 -arrow up
    
    # Configure the column that have the tree
    $T configure -treecolumn C0
    
    # Intento de desactivar el drag
    $T column dragconfigure -enable no
    
    # Create elements
    $T element create elemImgAny image
    $T element create elemTxtRead text -fill [list $SystemHighlightText {selected focus}] -lines 1
    $T element create elemRectSel rect -fill [list $SystemHighlight {selected focus} gray {selected !focus}] -showfocus yes
    $T element create eWindow window
    
    
    # Create styles using the elements
    set S [$T style create styAnyRead]
    $T style elements $S {elemImgAny elemRectSel elemTxtRead }
    $T style layout $S elemImgAny -expand ns
    $T style layout $S elemTxtRead -padx 4 -expand ns -squeeze x
    $T style layout $S elemRectSel -union [list elemTxtRead] -iexpand ns -ipadx 2
    
    set S [$T style create styFrame -orient horizontal]
    $T style elements $S {eWindow}
    $T style layout $S eWindow -squeeze x -padx {0 1} -pady {2 2}
    
    # Items
    set item root
    $T item configure $item -button yes
    $T item style set $item C0 styAnyRead
    $T item element configure $item C0 elemTxtRead -text [= "Root"]
    
    # List of lists: {column style element ...} specifying text elements
    # the user can edit
    TreeCtrl::SetEditable $T {
        
    }
    
    # List of lists: {column style element ...} specifying elements
    # the user can click on or select with the selection rectangle
    TreeCtrl::SetSensitive $T {
        {C0 styAnyRead elemTxtRead}
    }
    
    # Some notify install
    $T notify install <Drag-receive>
    
    # Notify bind
    # TODO

    $T notify bind DragTag <Drag-receive> { ::KMProps::ReceiveDragGroups  %T %l %I }
    #$T notify bind EditTag <Edit-accept> { ::KMProps::SetPropsToRename %T %I %t }

    bind $T <Button-1> [list ::KMProps::ClickTree %x %y $T]
    bind $T <Double-Button-1> [list ::KMProps::DoubleClickTree %x %y $T]
    bind $T <Return> [list ::KMProps::IntroEvent $T]
    bind $T <Button-3> "[list ::KMProps::RightClick %x %y $T]"
    
    # MouseWheel
    if {[string equal "x11" [tk windowingsystem]]} {
        # Support for mousewheels on Linux/Unix commonly comes through mapping
        # the wheel to the extended buttons.  If you have a mousewheel, find
        # Linux configuration info at:
        #        http://www.inria.fr/koala/colas/mouse-wheel-scroll/
        bind $T <4> [list $T yview scroll -3 units ]
        bind $T <5> [list $T yview scroll 3 units ]
    } elseif {[string equal [tk windowingsystem] "aqua"]} {
        bind $T <MouseWheel> [subst -nocommands { $T yview scroll {eval [expr - (%D)]} units } ]
    } else {
        bind $T <MouseWheel> [subst -nocommands { $T yview scroll [expr - (%D / 120) * 3] units } ]    
    }
    
    # Grid the tree
    
    if { $::KMProps::UseFullTkTree} {
        grid $T -sticky nsew
    } else {
        grid $T $vsb -sticky nsew
        grid configure $vsb -sticky ns
        grid $hsb -sticky ew
        grid remove $vsb $hsb
        bind $T <Configure> [list ConfigureListScrollbars $T $hsb $vsb]
    }
    
    grid columnconfigure $w 0 -weight 1
    grid rowconfigure $w 0 -weight 1
    
    bindtags $T [list $T TreeCtrlFileList TreeCtrl [winfo toplevel $T] all]
    
    return $T

}

proc ::KMProps::RefreshTree { {T ""} {onlySave 0} } {
    # ABSTRACT: Refresh the treectrl tree properties 
    # Arguments
    # T         => The tree path
    # onlySave  => Delete the tree and fills it with new items
    variable TreePropsPath; variable winpath
    variable lastSelected
    global KPriv
    if {$T == ""} {
        set T $TreePropsPath
    }
    
    # Primero hay que asegurarse de que exista el árbol
    if {([winfo exists $winpath]) && ([winfo exists $T])} {
        
        set aux [$T item range 0 end]

        foreach item $aux {
            
            set fullname [DecodeName [$T item tag names $item]]
            
            if {($fullname != "") && ([::xmlutils::getXmlNodeName $fullname] != "Item")} {
                
                catch {
                    ::xmlutils::setXml $fullname open "write" [$T item isopen $item]
                }
            }
        }
    }
    
    if { $onlySave == 0 } {
        
        set selectedItem [$T selection get]
        
        # Vacía el árbol y lo vuelve a llenar
        # Delete the material tree
        ::KMat::DeleteTree $T
        # Fill the model tree
        #::KMProps::FillTreePropsRecursive
        ::KMProps::FillTreeProps
        
        
        # Si está disponible añadimos la selección anterior
        catch {
            $T selection add $selectedItem
            $T see [$T selection get]
        }
    }

    # Release the selected item
    set lastSelected {}
}

proc ::KMProps::GetAllChildren { parent} {
    set lst {}
    # lappend lst $parent
    foreach wc [ winfo children $parent] {
        set lwc [ ::KMProps::GetAllChildren $wc]
        lappend lst {*}$lwc
    }
    if { [ llength $lst] == 0} {
        lappend lst $parent
    }
    return $lst
}

proc ::KMProps::LookForWidget { parent widget} {
    set ret ""
    if { [ winfo exists $parent.$widget]} {
        return $parent.$widget
    } else {
        foreach wc [ winfo children $parent] {
            set ret [ ::KMProps::LookForWidget $wc $widget]
            if { $ret != ""} {
                return $ret
            }
        }
    }
    return $ret
}


proc ::KMProps::CreateTreeMaterials { w } {
    # ABSTRACT: Create the treectrl material tree path
    # Arguments
    # w => Frame path
    
    set ::KMat::TreeMatsPath "$w.treeMat"
}

proc ::KMProps::ThereAreChangesPending { f } {
    #Deteccion de si hay frame abierto y popup de aviso que debe dar al ok-cancel antes de cambiar
    #Es una solucion rapida, lo que se deberia hacer es controlar si ha habido cambios, y que saltara un pop up de
    #si quiere guardar o descartar los cambios
    variable selGroup
    if {[winfo exists $f]} {
        # hay ya una bottom frame abierto
        # check presence of bBottomOk and bBottomCancel
        # two layouts: $f.bBottomOk and $f.fbuttons.bBottomOk
        # so let's go through all the children lookking for the buttons
        set but_ok [ ::KMProps::LookForWidget $f bBottomOk]
        set but_cancel [ ::KMProps::LookForWidget $f bBottomCancel]
        if { ( $but_ok != "") && ( $but_cancel != "")} {
            #el bottom frame tiene un ok y un cancel      
            if { [info exists selGroup] && $selGroup!="" } {
                #Not Empty selected group
                WarnWin [= "Click Ok or Cancel in the bottom frame before selecting a new item"] $f
                return 1
            }
        }
    }
    #end Deteccion...
    return 0
}

proc ::KMProps::DoubleClickTree { x y T {item ""}} {
    # ABSTRACT: Manage the double click in the tree
    # Arguments
    # x    => x position
    # y    => y position
    # T    => Tree path 
    # item => Selected item
    variable lastSelected
    
    # Si no llega directamente el item, miramos cual ha sido pulsado
    if { $item == "" } {
        
        set info [$T identify $x $y]
        # wa "info:$info"        
        if {([lindex $info 0] == "item") && ([llength $info] >= 4)} {
            set item [lindex $info 1]
        } else {
            return ""
        }
    }
    # msg $item
    set fullname [DecodeName [$T item tag names $item]]
    set idFull [string map { "." "" "//" ""} $fullname]
    # wa "fullname:$fullname idFull:$idFull"
    if { [::xmlutils::setXml $fullname state] == "disabled" } {
        return ""
    }

    #Deteccion de si hay frame abierto y popup de aviso que debe dar al ok-cancel antes de cambiar
    #Es una solucion rapida, lo que se deberia hacer es controlar si ha habido cambios, y que saltara un pop up de
    #si quiere guardar o descartar los cambios
    variable NbPropsPath
    if { [ ::KMProps::ThereAreChangesPending $NbPropsPath.fBottom ]} {
        return ""
    }
    #end Deteccion...
    
    set id [::xmlutils::setXml $fullname id]
    set dv [::xmlutils::setXml $fullname dv]
    # wa "id:$id dv:$dv"
    # Eliminamos el anterior combo, si aun está visible
    if {([llength $lastSelected]) && ([lindex $lastSelected 0] == $item)} { 
        ::KMProps::cmbSelectChange [lindex $lastSelected 0] $T 1 "anterior"
    }
    
    # Destroy the bottom frame
    ::KMProps::DestroyBottomFrame
    
    if { [::xmlutils::getXmlNodeName $fullname]  == "Item" } {
        
        set f [::KMProps::buildFrame $T $item]  
        $T item style set $item C0 styFrame
        $T item element configure $item C0 eWindow -window $f
        
        # Nos guardamos este item como y el valor seleccionado (es el último visible)
        set f "$T.f$idFull.cmb"
        set selCombo [::xmlutils::getComboDv $f $fullname]
        set selComboText [::xmlutils::getComboDv $f $fullname "text"]
        set lastSelected [list $item $selCombo $selComboText]
        
    } else {
        set ClassType [::xmlutils::setXml $fullname class]
        set idTemplate [::xmlutils::setXml $fullname idTemplate]
        # wa "ClassType:$ClassType idTemplate:$idTemplate fullname:$fullname"
        
        if {$ClassType == "Groups"} {
            # Comprobamos si en esta asignación de Grupos necesitarán propiedades
            #set template [::KMProps::copyTemplate ${idTemplate} $fullname "Nothing" "OnlyGetText"]
            
            # Validamos que exista algún combo de "properties" en el template
            #if {[string match "*GCV=\"Properties*" $template]} {
                
                # Miramos si hay alguna propiedad dada de alta (si la hay,el combo no puede estar vacío)
            #    set props [::KMProps::getProps $fullname]
            #    if { [llength $props] == 0 } {
            #        WarnWin [= "You must define a valid Property for this element type."]
                    
                    # Abrir la edición de propiedades
            #        $T selection clear
            #        set parentItem [$T item parent [$T item parent $item]]
            #        foreach i [$T item children $parentItem] {
            #            if { [$T item text $i 0] == "Properties" } {
            #                ::KMProps::DoubleClickTree 0 0 $T $i
            #            }
            #        }
            #        
            #        return  ""
            #    }
            #}
            ::KMProps::buildGroupsFrame $T $idTemplate $item $fullname
            
        } elseif {$ClassType == "Properties" } { 
            ::KMProps::buildPropertyFrame $T $idTemplate $item $fullname
            
        } elseif {($ClassType == "Tab") || ($ClassType == "Property") || ($ClassType == "Group")} {
            return 0
            ::KMProps::buildTabFrame $T $item $ClassType
            
        } else {
            $T expand -recurse "$item"
        }
    }
}

proc ::KMProps::ClickTree { x y T } {
    # ABSTRACT: Manage the a click in the tree
    # Arguments
    # x    => x position
    # y    => y position
    # T    => Tree path 
    variable lastSelected
    
    #Deteccion de si hay frame abierto y popup de aviso que debe dar al ok-cancel antes de cambiar
    #Es una solucion rapida, lo que se deberia hacer es controlar si ha habido cambios, y que saltara un pop up de
    #si quiere guardar o descartar los cambios
    variable NbPropsPath
    if { [ ::KMProps::ThereAreChangesPending $NbPropsPath.fBottom ]} {
        return ""
    }
    #end Deteccion...

    set info [$T identify $x $y]
    
    if { [lindex $info 0] == "item" && [llength $info] >= 4 } {
        
        set item [lindex $info 1]
        set col [lindex $info 3]
        
        set fullname [DecodeName [$T item tag names $item]]
        set id [::xmlutils::setXml $fullname id]
        
        # Delete the last combo, if it is still visible
        if {([llength $lastSelected]) && ([lindex $lastSelected 0] != $item)} { 
            
            ::KMProps::cmbSelectChange [lindex $lastSelected 0] $T 1 "anterior"
        }
        
        # We get the kind of node to act accordingly
        set ClassType [::xmlutils::setXml $fullname class]
        
        # Let's see if you pressed container or item
        set nodeName [::xmlutils::getXmlNodeName $fullname]
        if { $nodeName == "Container" } {
            # Abrimos el item si tiene hijos y no es de alta de propiedades o grupos
            if { $ClassType == "Groups" || $ClassType == "Properties" } {
                
                # En este caso desplegamos el frame de alta de propiedades o grupos
                # comentado al cambiar muy a menudo el frame de abajo produce que se pierda la informacion introducida por el usuario
                # si se gestiona mejor el control de cambios se puede volver a activar, (relacionado con 
                #::KMProps::DoubleClickTree $x $y $T

            } else {
                # $T item toggle $item
            }
        } elseif {$nodeName == "Item"} {
            
        }
        
    } elseif { [lindex $info 0] == "header" && [lindex $info 1] == "0" } {
        
        if { [$T column cget C0 -arrow] == "up" } {
            $T column configure C0 -arrow down
        } else {
            $T column configure C0 -arrow up
        }
        
        return ""
    } else {
        return ""
    }
    
    if { $col != 0 } {
        if { ![$T selection includes $item] } {
            $T selection clear
            $T selection add $item
        }
    }
    if { $col == 0 } {
    } elseif { $col == 1 } {
        WarnWin "column 1 ???"
    }
    
    if { $col != 0 } {
        return -code break
    }
    return ""
}
proc ::KMProps::RightClick { x y T } {
    $T identify -array clicked $x $y
    if { $clicked(where)== "item" } {
        $T selection clear
        $T selection add $clicked(item)    
        update idletasks
    }
    ::KMProps::MenuContextual $T $x $y  
}

proc ::KMProps::IntroEvent { T } {
    
    ::KMProps::DoubleClickTree 0 0 $T
}

proc ::KMProps::editTag { T item fullname newtext { newPath "" }  } {
    
    if { $newPath == "" } {
        set parts [::KMProps::split2 $fullname //]
        lset parts end $newtext
        set newPath [join $parts //]
    }
    
    if { $fullname != $newPath } {
        
        # Rename from the tree
        $T item tag remove $item [list names [$T item tag names $item]]
        
        $T item tag add $item [EncodeName $newPath]
        $T item element configure $item C0 elemTxtRead -text $newtext
        
        return $newPath
        
    } else {
        
        return $fullname
    }
}

proc ::KMProps::BeginEditGroups { T } {
    # 
    set I [$T item id active]
    set C 0
    set E elemTxtRead
    
    ::TreeCtrl::FileListEdit $T $I $C $E
}

proc ::KMProps::ReceiveDragGroups { T dragged_list dst } {

}

proc ::KMProps::InsertNewProp { node id T {parent ""} {parentitem root} {childs true} {state "normal"} {open 0}} {
    
    set propName [$node getAttribute pid ""]
    # wa "propName:$propName node:$node id:$id T:$T parent:$parent parentitem:$parentitem childs:$childs state:$state open:$open"
   
    if { $state == "hidden" } {
        #No inserta este item pero se sigue el proceso
        return $parentitem
        
    } elseif { $state == -1 } {
        #No inserta este item ni su descendencia
        return -1
    }
    
    if {[string first "NoTree" $id] != -1 } {
        #Solo insertamos determinados nodos (container, item..)
        return "-1"
    }
    
    if {$id == ""} {
        WarnWin [$node getAttribute id ""]
        #set id [$node getAttribute gpid ""]
    }
    
    if {$parent != ""} {
        set fullname $parent$id
    } else {
        set fullname $id
    }

    set tooltip [::xmlutils::setXml $fullname help]
    # wa "fullname:$fullname"
    if {$childs} {
        
        set item [$T item create -button yes -tags [EncodeName $fullname] -open $open]
        if {$::KMProps::UseFullTkTree} {
            $T popup_enter_help $item [= $tooltip]
        }
        $T item lastchild $parentitem $item
        $T item style set $item C0 styAnyRead
        $T item element configure $item C0 elemTxtRead -text [= $propName]
        
    } else {
        
        set item [$T item create -button no -tags [EncodeName $fullname]]
        if {$::KMProps::UseFullTkTree} {
            # Looks like this is the one that creates the annoying tooltip
            $T popup_enter_help $item [= $tooltip]
        }
        $T item lastchild $parentitem $item
        
        if {[string range $id 0 1] == "i."} {
            
            #set f [::KMProps::buildFrame $T $item]
            
            $T item style set $item C0 styFrame                                                
            #$T item element configure $item C0 eWindow -window $f
            
            set dv [::xmlutils::setXml $fullname dvText]
            
            $T item style set $item C0 styAnyRead
            if {$state != "disabled"} {
                
                #Miramos si tiene algun estilo especial
                if { [::xmlutils::setXml $fullname style] == "*" } {
                    $T item element configure $item C0 elemTxtRead -text "[= $propName]* : $dv"
                } else {
                    $T item element configure $item C0 elemTxtRead -text "[= $propName]: $dv"
                }
                
            } else {
                
                if {[gid_themes::GetCurrentTheme] == "GiD_black"} {
                    $T item element configure $item C0 elemTxtRead -text "[= $propName]: $dv" -fill { darkgreen }
                } else {
                    $T item element configure $item C0 elemTxtRead -text "[= $propName]: $dv" -fill { gray }
                }
            }
            
            
            
        } elseif { [string range $id 0 1] == "c." } {                                                                                                                   
            
            $T item style set $item C0 styAnyRead
            $T item element configure $item C0 elemTxtRead -text [= $propName]
        }
    }
    
    # Consultamos el icono en el xml, y si existe en nuestro directorio se lo añadimos                
    set icon [::xmlutils::setXml $fullname icon]
    set imagen [::WinUtils::GetImage $icon]
    if { $imagen != -1 } {
        $T item image $item C0 $imagen
    } else {
        # Si no encuentra la imagen no mostramos ninguna (descomentando mostramos defaultIcon)
        #$T item image $item C0 [::WinUtils::GetImage "defaultIcon.gif"]
    }
    
    # Miramos si el item tiene que estar a disabled (viene del proc ::KMProps::stateNode)
    if {$state == "disabled"} {
        $T item enabled $item 0
    }
    
    return $item
}

proc ::KMProps::changeTypeEvent {T fullname item {changeTo ""} } {
    
    if { $changeTo == "" } {
        set tipo [ set ::KMProps::type${item}]
    } else {
        set tipo "$changeTo"
    }
    
    ::xmlutils::setXml $fullname type "write" $tipo
    
    #Cambiamos también el tipo de todos los item seleccionados
    set items [$T selection get]
    
    for {set i 0} { $i < [llength $items] } {incr i} {
        
        set item [lindex $items $i]
        set ::KMProps::type$item $tipo
        
        set fullname [DecodeName [$T item tag names $item]]
        
        ::xmlutils::setXml $fullname type "write" $tipo
    }
}

proc ::KMProps::MenuContextual { T x y } {
    
    set w $T.menucontextual
    if { [winfo exists $w] } {
        destroy $w
    }
    
    menu $w
    
    #$w add command -label [= "New Group#C#layer"] -command [list ::KMProps::CreateNewGroupId $T]
    
    set item [$T selection get]
    #msg "item$item"
    set numchilds 0
    if { $item != "" } {
        set fullname [DecodeName [$T item tag names $item]]
        set class [::xmlutils::setXml $fullname class]
        set numchilds [llength [$T item children $item]]
        
        if {$class == "Groups" || $class == "Properties"} {
            $w add command -label [= "New"] -command [list ::KMProps::DoubleClickTree 0 0 $T $item] -state normal
            $w add separator
        }
        
        if {$class == "Group"} {
            $w add command -label [= "Delete"] -command [list ::KMProps::deleteGroupCondition $T $item] -state normal
        } elseif {$class == "Property"} {
            $w add command -label [= "Delete"] -command [list ::KMProps::deleteProps $T $item] -state normal
        } else {                
            $w add command -label [= "Delete"] -state disabled
        }
        
        
        if {$class == "Group" || $class == "Property"} {
            $w add command -label [= "Edit"] -command [list ::KMProps::DoubleClickTree 0 0 $T $item] -state normal
        } else {
            $w add command -label [= "Edit"] -state disabled
        }
        
        $w add separator
        
        if {$class == "Groups"} {
            if { $numchilds > 0 } {
                $w add command -label [= "Draw Entities"] -command [list ::KMProps::drawGroupsCondition $T $item] -state normal
                $w add command -label [= "List Entities"] -command [list ::KMProps::listGroupsCondition $T $item] -state normal
            } else {
                $w add command -label [= "Draw Entities"] -state disabled
                $w add command -label [= "List Entities"] -state disabled    
            }
            $w add separator
        } elseif {$class == "Group"} {
            $w add command -label [= "Draw Entities"] -command [list ::KMProps::drawGroupCondition $T $item] -state normal
            $w add command -label [= "List Entities"] -command [list ::KMProps::listGroupCondition $T $item] -state normal
            $w add separator
        }
        
    } else {
        set item "root"                
    }                
    
    
    if { $numchilds > 0 } {
        $w add command -label [= "Collapse All"] -command [list $T collapse -recurse "$item"] -state normal
        $w add command -label [= "Expand All"] -command [list $T expand -recurse "$item"] -state normal
        $w add separator
        $w add command -label [= "List Subtree"] -command [list ::KMProps::listSubtree $T $item] -state normal
    } else {
        $w add command -label [= "Collapse All"] -state disabled
        $w add command -label [= "Expand All"] -state disabled
        $w add separator
        $w add command -label [= "List Subtree"] -state disabled  
    }
    
    set x [expr [winfo rootx $T]+$x+2]
    set y [expr [winfo rooty $T]+$y]
    GiD_PopupMenu $w $x $y
}

proc ::KMProps::deleteItem { T item } {
    
    #Comprobamos q el item exista
    if { [$T item id $item] != "" } {
        
        #msg "T:$T item$item"
        
        #Si el padre se queda sin hijos le quita la imagen
        set padre [$T item parent $item]
        
        #Elimina el item del árbol
        $T item delete $item 
        
        #Si se ha borrado bien, miramos si el padre sigue teniendo hijos
        if { [$T item numchildren $padre] == 0 } {
            #Si no quedan hijos eliminamos el botón de desplegar
            $T item configure $padre -button 0
        }
    }
}

#
# Elimina todas las asignaciones de grupos que existan en el árbol de propiedades
#
proc ::KMProps::deleteProps { T itemProp {editName ""}} {
    
    set propsToDelete {}
    
    # Recorremos tooodo el árbol
    set items [$T item descendants "root"]
    foreach item $items {
        
        #Comprobamos q el item aun exista
        if { [$T item id $item] != "" } {
            
            set fullname [DecodeName [$T item tag names $item]]
            set itemDv [::xmlutils::setXml $fullname dv]
            set propId [$T item text $itemProp 0]
            
            if { $itemDv == $propId } {
                
                
                if {[::xmlutils::setXml $fullname GCV] == "Properties"} {
                    if {$editName == ""} {
                        
                        lappend propsToDelete $item
                    } else {
                        #Renombra el nodo en el xml
                        ::xmlutils::setXml $fullname dv "write" $editName
                        #::xmlutils::setXml $fullname id "write" $editName
                    }
                }
            }
        }
    }
    if { $editName == "" } {
        if { [llength $propsToDelete] > 0 } {
            set aviso [= "This property it is assigned to %s groups.\nDo you want to delete the property and all this groups?" [llength $propsToDelete]]
            set answer [::WinUtils::confirmBox "." "$aviso"]
            if { $answer == "ok" } {
                
                #Borramos todos los grupos
                foreach item $propsToDelete {        
                    
                    #No hay que borrar la propiedad, sino el grupo que la incluye
                    set itemGroup [$T item parent $item]
                    set fullname [DecodeName [$T item tag names $itemGroup]]
                    if {[::xmlutils::setXml $fullname class] == "Group"} {
                        
                        #Elimina el grupo del xml
                        ::xmlutils::unsetXml $fullname
                        
                        #Elimina el grupo del árbol (ya no hace falta, se hace refresh)
                    }
                }
                
                #Borramos la propiedad
                set fullProp [DecodeName [$T item tag names $itemProp]]
                ::xmlutils::unsetXml $fullProp
            }
        } else {
            set aviso "Do you want to remove the property '[$T item text $itemProp 0]'?" 
            set answer [::WinUtils::confirmBox "." "$aviso"]
            if { $answer == "ok" } {
                #Borramos la propiedad
                set fullProp [DecodeName [$T item tag names $itemProp]]
                ::xmlutils::unsetXml $fullProp
            }
        }
    }
    
    #Cierra la ventana de propiedades y la vuelve a cargar para q se renombren
    # los items automáticamente
    ::KMProps::RefreshTree $T
}

#
# Elimina todas las asignaciones de grupos que existan en el árbol de propiedades
#
proc ::KMProps::checkMaterials { materialId {editName ""} } {
    
    global KPriv
    
    # Recorremos el xml
    set nodes [$KPriv(xml) getElementsByTagName "Item"]
    
    foreach node $nodes {
        # wa "GCV:[$node getAttribute GCV ""]"
        # Si encuentra el grupo intenta borrar el nodo y su descendencia
        if { [$node getAttribute GCV ""] == "Materials" } {
            # wa "id:[$node getAttribute id ""]\ndv:[$node getAttribute dv ""]"
            if { [$node getAttribute dv ""] == $materialId } {
                $node setAttribute dv $editName
            }
        }
    }
    
    # Refresh the tree
    ::KMProps::RefreshTree
}

proc ::KMProps::checkFunctions { functionId {newName ""} {action "execute"}} {
    
    global KPriv
    set numChanges 0
    
    #Si no está abierta la ventana, deberemos recorrer el xml a mano
    set nodes [$KPriv(xml) getElementsByTagName "Item"]
    
    foreach node $nodes {
        
        #Si encuentra el grupo intenta borrar el nodo y su descendencia
        if { [$node getAttribute function ""] == 1 } {
            if { [$node getAttribute dv ""] == "$functionId" } {
                
                if { $action == "execute" } {
                    if {$newName == ""} {
                        #Resetea la función
                        $node setAttribute state ""
                        $node setAttribute function 0
                        $node setAttribute dv "0.0"
                    } else {
                        # Rename the function
                        $node setAttribute dv $newName
                    }
                }
                set numChanges [expr $numChanges + 1]
            }
        }
    }
    if { $numChanges > 0 && $action == "execute"} {
        if { [winfo exists $::KMProps::WinPath] } {
            ::KMProps::DestroyBottomFrame
            ::KMProps::RefreshTree
        }
    }
    return $numChanges
}

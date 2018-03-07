###############################################################################
#
#    NAME: curves.tcl
#
#    PURPOSE: Work with curves
#
#    QUANTECH ATZ-DEVELOPMENT DEPARTMENT
#
#    AUTHOR : G. Socorro
#
#    CREATED AT: 21/08/12
#
#    HISTORY:
# 
#     0.4-26/09/12- J. Garate, Curve Get/Set Properties
#     0.3-26/09/12- J. Garate, Curve Renaming functionality. SPD modification functions.// Keyboard binding
#     0.2-20/09/12- J. Garate, Creation of Curves Interface and managment functions
#     0.1-21/08/12-G. Socorro, create the base source code
#
###############################################################################

 namespace eval ::KCurves {
	variable TreeCurvePath
    variable WinLayout;	variable SystemHighlight
	variable SystemHighlightText;
}

proc ::KCurves::Init { } {
	variable WinLayout;	variable SystemHighlight
	variable SystemHighlightText;
	global KPriv
    
	# Get default colors
	set w [listbox .listbox]
	set SystemHighlight [$w cget -selectbackground]
	set SystemHighlightText [$w cget -selectforeground]
	destroy $w
	if { $::tcl_platform(platform) == "unix" } {
		# I hate that gray selection color
		set SystemHighlight #316ac5
		set SystemHighlightText White
    }

    set WinLayout "OUTSIDE"
    
    if {[info exists KPriv(xml)]} {
		set ::KCurves::xml $KPriv(xml)
    }
}

proc ::KCurves::CreateTreeAndToolbar { w } {
    
    # Create the treectrl properties
    set mdf [ttk::frame $w.middle ]
    set T [::KCurves::CreateTreeProperties $w]
    
    grid $w.middle -sticky wens
    
    # Create the frame where set the properties
    set f [ttk::frame $w.fBottom -borderwidth 0 ]
    set ::KCurves::BottomframePath $f

    # Grid for toolbar
    grid $f -row 2 -column 0 -sticky wes
    
    # Create the toolbar frame
    set tbf [ttk::frame $w.tbar -borderwidth 0]
    
    set ::KCurves::Buttonframe $tbf
    
    ttk::button $tbf.newCurve -image [::WinUtils::GetImage add.gif]  -command [list ::KCurves::CreateNewCurve] -style Toolbutton
    tooltip::tooltip $tbf.newCurve [= "Create new curve"]

    ttk::button $tbf.deleteCurveId -image [::WinUtils::GetImage delete_icon.gif]  -command [list ::KCurves::DeleteCurve] -style Toolbutton
    tooltip::tooltip $tbf.deleteCurveId [= "Delete curve"]

    ttk::button $tbf.refreshTree -image [::WinUtils::GetImage refresh.gif]  -command [list ::KCurves::RefreshTree] -style Toolbutton
    tooltip::tooltip $tbf.refreshTree [= "Refresh Tree"]
    
 
    grid $tbf -sticky ews
    grid anchor $tbf w

    # Grid for toolbar 
    grid $tbf.newCurve -sticky we -row 0 -column 0
    grid $tbf.deleteCurveId -sticky we -row 0 -column 1
    grid $tbf.refreshTree -sticky we -row 0 -column 2
    
    focus $T
    return $T
}

proc ::KCurves::CreateTreeProperties {w} {

    variable SystemHighlight
    variable SystemHighlightText
    set SystemHighlight #316a00
    set SystemHighlightText White


    # Scrollbars
    set vsb $w.vsb1
    set hsb $w.hsb1

    # Create the treectrl and set the scrollbars
    set T [treectrl $w.t -xscrollcommand [list $hsb set] -yscrollcommand [list $vsb set]]
    ttk::scrollbar $vsb -orient vertical   -command [list $T yview]
    ttk::scrollbar $hsb -orient horizontal -command [list $T xview]

    # Set the height
    set height [font metrics [$T cget -font] -linespace]
    if {$height < 19} {
	set height 19
    }

    # Configure the treectrl
    $T configure -indent 15 -itemheight $height -selectmode browse \
	-showroot 0 -showrootbutton 0 -showbuttons 1 -showlines 1 \
	-highlightthickness 0 -borderwidth 0 -height 300 \
	-xscrollincrement 20 -yscrollincrement 20 
    
    $T column create -text "Curves" -tags C0 -weight 0
    
    # Configure the column weight and arrow
    $T column configure C0 -weight 1 -arrow up

    # Configure the column that have the tree
    $T configure -treecolumn C0
    
    # Intento de desactivar el drag
    #$T column dragconfigure -enable no

    # Create elements
    $T element create elemImgAny image
    $T element create elemTxtRead text -fill [list $SystemHighlightText {selected focus}] -lines 1
    $T element create elemRectSel rect -fill [list $SystemHighlight {selected focus} gray {selected !focus}] -showfocus yes
    $T element create eWindow window
    
    # Create styles using the elements
    set S [$T style create styAnyRead]
    $T style elements $S {elemImgAny elemRectSel elemTxtRead}
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

    # List of lists: {column style element ...} specifying text elements the user can edit
    TreeCtrl::SetEditable $T {

    }

    # List of lists: {column style element ...} specifying elements
    # the user can click on or select with the selection rectangle
    TreeCtrl::SetSensitive $T {
	{C0 styAnyRead elemTxtRead}
    }

    # List of lists: {column style element ...} specifying elements
    # added to the drag image when dragging selected items
    TreeCtrl::SetDragImage $T {
	
    }
    
    # Some notify install
    $T notify install <Edit-accept>

    # Notify bind
    # TODO

    #		$T notify bind DragTag <Drag-receive> [list ::KMProps::ReceiveDragGroups  %T %l %I]
    $T notify bind EditTag <Edit-accept> [list ::KCurves::SetCurveToRename %T %I %t]

    #bind $T <Button-1> [list ::KCurves::ClickTree %x %y $T]
    bind $T <Double-Button-1> [list ::KCurves::DoubleClickTree %x %y $T]
    bind $T <F2> [list ::KCurves::BeginEditCurve $T]
    bind $T <Key-Delete> [list ::KCurves::DeleteCurve] 

    bind $T <Button-3> "[list ::KCurves::MenuContextualGroup %W %x %y] ; break"

    # Grid the tree
    grid $T $vsb -sticky nsew
    grid configure $vsb -sticky ns
    grid $hsb -sticky ew
    grid remove $vsb $hsb
    bind $T <Configure> [list ConfigureListScrollbars $T $hsb $vsb]

    grid columnconfigure $w 0 -weight 1
    grid rowconfigure $w 0 -weight 1

    bindtags $T [list $T TreeCtrlFileList TreeCtrl [winfo toplevel $T] all]

    return $T
}

proc ::KCurves::SetCurveToRename { T item newtext } {
    
    global KPriv
    if { $newtext == ""} {  return  }
    if { $item == 0 } {
        msg "Root folder can't be edited"
        return
    }
    
    #Controlamos que el nombre no esté ya en el árbol (a no ser q no se haya cambiado)
    set oldId [$T item text $item 0]
    
    if { $oldId == $newtext } {
        return ""
    }
    
    set parent [$T item parent $item]
    set bros [$T item children $parent]
    foreach bro $bros {
        set broId [$T item text $bro 0]
        if { $broId == $newtext } {
            msg "The curve name '$newtext' already exist.\nChoose another, please." 
            return ""
        }
    }

    #Validamos q el nombre no tenga carácteres que vulneran la seguridad y quitamos espacios
    set newtext [::KUtils::parseTreeStr $newtext 1]
    if { $newtext == -1 } {
        msg "You can't use some reservate chars like:\n  :   /   $   .   \\  %  "
        return
    }

    ::KCurves::SetCurveProp "item" $item "pid" $newtext
}

proc ::KCurves::FillTreeCurves { } {
	
	global KPriv
	set T $KPriv(TreeCurvePath)
	
	#Seleccionamos todos los nodos del primer nivel
	set nodes [$::KPriv(xml) selectNodes "/Kratos_Data/RootData\[@id=\"Curves\"\]/Container\[@id]"]
	
	foreach node $nodes {
		# Insertamos cada RootData de 1er nivel en el árbol
		set item [::KCurves::InsertNewProp $node [$node getAttribute id ""] $T "" "root" [$node hasChildNodes] [::KMProps::stateNode $node] [$node getAttribute open "0"]]
		if {$item != -1} {
			set path "[$node getAttribute id 0]//"
		        ::KCurves::FillRecursiveChilds $T $path $node $item
		}
	}
}

proc ::KCurves::FillRecursiveChilds { T path node item} {
	
	set nodes2 [$node childNodes]
	set pathcp $path
	foreach node2 $nodes2 {
		set item2 [::KCurves::InsertNewProp $node2 [::KMProps::splitNode $node2] $T $path "$item" [$node2 hasChildNodes] [::KMProps::stateNode $node2] [$node2 getAttribute open "0"]]
	        if {$item2 != -1} {
			
			#set pathcp [join [concat $path "[::KMProps::splitNode $node2]//"] ""]
      set pathcp [regsub -all "// " [concat $path "[::KMProps::splitNode $node2]//"] "//"]
			
			::KMProps::FillRecursiveChilds $T $pathcp $node2 $item2
			#set pathcp [join [concat $path "[::KMProps::splitNode $node2]//"] ""]
      set pathcp [regsub -all "// " [concat $path "[::KMProps::splitNode $node2]//"] "//"]
      
		}
	}
}

proc ::KCurves::CreateNewCurve { } {
    # Generates a new curve from Template
    
    global KPriv
    ::KTables::CancelButton
    set xml $KPriv(xml)
#   msg "STEP 1: Locate template and transform to Container"
    set nodeTempl [[$xml find id "NewCurveId"] asList]
#   msg "Pre: nodeTempl : $nodeTempl"
    set aux1 $nodeTempl
    set aux1 [lreplace $aux1 [lsearch $aux1 Template] [lsearch $aux1 Template] Container]
#   msg "Post aux1 : $aux1"

#   msg "STEP 2: Change curve names"
    set aux [lindex $aux1 1]
#   msg "Pre: aux $aux"
    set id [::KCurves::NewCurveName]
    set newaux [lreplace $aux 3 3 [lindex $id 1]]
    set newaux [lreplace $newaux 1 1 [lindex $id 0]]
#   msg "Post: newaux $newaux"
	set aux1 [lreplace $aux1 [lsearch $nodeTempl $aux] [lsearch $nodeTempl $aux] $newaux]
    
    set baseCurvePath "/Kratos_Data/RootData\[@id='Curves'\]"
    set baseNode [$xml selectNodes $baseCurvePath]
#   msg "\nGlobalPost appending aux1 : $aux1"
	$baseNode appendFromList $aux1
    
    # Refresh the Tree
    ::KCurves::RefreshTree
    set item [::KCurves::GetItemFromCurveNodeList $aux1]
    set T $KPriv(TreeCurvePath)
    set item [::KCurves::FindMyCurve $T $item]
    $T selection modify $item [$T selection get]
    ::KCurves::DoubleClickTree 0 0 $T $item
}

proc ::KCurves::NewCurveName { } {
    # Returns a valid Curve name
    global KPriv
    
    set i 0
    while { 1 } {
        if { [$KPriv(xml) find id "NewCurve$i"] == "" } { 
            set list "NewCurve$i" 
            lappend list "New Curve $i"
            return $list
        }
        incr i 1
    }
}

proc ::KCurves::GetItemFromCurveNodeList { node } {
    global KPriv
    set T $KPriv(TreeCurvePath)
    set i 0
    set itemfullname [::KCurves::GetFullname $i]

    set pid [lindex $node 1]
    set pid [lindex $pid 3]
    set aux [$T item range 0 end]
    foreach i $aux {
        if {[$T item text $i 0] == $pid} {
            return $i
        }
    }
}

proc ::KCurves::DeleteCurve { } {
    # Delete the selected item's curve
    global KPriv
    set T $KPriv(TreeCurvePath)
    
    set nitems [$T selection count]
    set items [$T selection get]
    ::KTables::CancelButton
    if { $nitems > 0 } {
        set items [lsort -integer $items]
        set pids ""
        set fullnames ""
        foreach item $items {
            set curve [::KCurves::FindMyCurve $T $item]
            set fullname [::KCurves::GetFullname $curve]
            set fullnames [lappend fullnames $fullname]
            set pid [::KCurves::setXml $fullname pid "" ]
            set pids [ append pids $pid ", "]
        }
        set pids [string range $pids 0 end-2]
        set aviso "Are you sure you want to delete the curve/s \"$pids\"?"
        set confirmado [::WinUtils::confirmBox "." "$aviso"]
        
        if {$confirmado == "ok" } {
            foreach fullname $fullnames {
                ::xmlutils::unsetXml $fullname
            }
        }
        ::KCurves::RefreshTree
        
    } else {
        WarnWin "Please select a Curve"
    }
}

proc ::KCurves::FindMyCurve { T item } {
    while { 1 } {
        set fullname [::KCurves::GetFullname $item]
        if {$fullname == ""} {
            return ""
        }
        if { [::KCurves::IsCurve $fullname] } {
            return $item
        } else {
            incr item -1
        }
    }
}

proc ::KCurves::InsertNewProp { node id T {parent ""} {parentitem root} {childs true} {state "normal"} {open 0}} {
	
	set propName [$node getAttribute pid ""]
	
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
	
	if { $id == "" } {
	WarnWin [$node getAttribute id ""]
	}
	
	if { $parent != "" } {
	set fullname $parent$id
	} else {
	set fullname $id
	}
	
	
	if { $childs } {
	
	set item [$T item create -button yes -tags [EncodeName $fullname] -open $open]
	$T item lastchild $parentitem $item
	$T item style set $item C0 styAnyRead
	$T item element configure $item C0 elemTxtRead -text "$propName"
	
	} else {
	
        set item [$T item create -button no -tags [EncodeName $fullname]]
        $T item lastchild $parentitem $item
        
        if { [string range $id 0 1] == "i." } {
            
            #set f [::KMProps::buildFrame $T $item]
            
            $T item style set $item C0 styFrame                                                
            #$T item element configure $item C0 eWindow -window $f

            set prev "Curves//c."
            append prev $fullname
            set fullname $prev
            set dv [::xmlutils::setXml $fullname dv "read" "" "curve"]

            $T item style set $item C0 styAnyRead
            if {$state != "disabled"} {
                
                #Miramos si tiene algun estilo especial
                if { [::xmlutils::setXml $fullname style] == "*" } {
                        $T item element configure $item C0 elemTxtRead -text "$propName* : $dv"
                } else {
                        $T item element configure $item C0 elemTxtRead -text "$propName: $dv"
                }
            
            } else {
                if {[gid_themes::GetCurrentTheme] == "GiD_black"} {
                $T item element configure $item C0 elemTxtRead -text "$propName: $dv" -fill { darkgreen }
                } else {
                    $T item element configure $item C0 elemTxtRead -text "$propName: $dv" -fill { gray }
                }
            }
        } elseif { [string range $id 0 1] == "c." } {                                                                                                                   
            
            $T item style set $item C0 styAnyRead
            $T item element configure $item C0 elemTxtRead -text "$propName"
        } elseif { [string range $id 0 1] == "t." } {

            $T item style set $item C0 styFrame
            set prev "Curves//c."
            append prev $fullname
            set fullname $prev
            set dv [::xmlutils::setXml $fullname dv "read" "" "curve"]

            $T item style set $item C0 styAnyRead
            if {$state != "disabled"} {
            
                #Miramos si tiene algun estilo especial
                if { [::xmlutils::setXml $fullname style] == "*" } {
                        $T item element configure $item C0 elemTxtRead -text "$propName* : $dv"
                } else {
                        $T item element configure $item C0 elemTxtRead -text "$propName: $dv"
                }
            
            } else {
                    
                if {[gid_themes::GetCurrentTheme] == "GiD_black"} {
                    $T item element configure $item C0 elemTxtRead -text "$propName: $dv" -fill { darkgreen }
                } else {
                    $T item element configure $item C0 elemTxtRead -text "$propName: $dv" -fill { gray }
                }
            }
        }
    }
	
	# Consultamos el icono en el xml, y si existe en nuestro directorio se lo añadimos                
	set icon [::xmlutils::setXml $fullname icon]
	set imagen [::WinUtils::GetImage $icon]
	if { $imagen != -1 } {
	$T item image $item C0 $imagen
	} else {
	# Si no encuentra la imagen no mostramos ninguna (descomentando mostramos defaultIcon)

	}
	
	# Miramos si el item tiene que estar a disabled (viene del proc ::KMProps::stateNode)
	if {$state == "disabled"} {
	$T item enabled $item 0
	}
	
	return $item
}

proc ::KCurves::MenuContextualGroup { T x y } {
    
    set w $T.menucontextualgroup
    if { [winfo exists $w] } {
        destroy $w
    }
    menu $w
    
    $w add command -label [= "New Curve"] -command [list ::KCurves::CreateNewCurve]
    $w add command -label [= "Delete Curve"] -command [list ::KCurves::DeleteCurve]
    
    set item [$T selection get 0]	
    if {$item != "" } {	
        set curve [::KCurves::FindMyCurve $T $item]
	    $w add command -label [= "Rename Curve"] -command [list ::KCurves::BeginEditCurve $T ]			
	    $w add command -label [= "Copy Curve"] -state disabled				
	
    }

    set x [expr [winfo rootx $T]+$x+2]
    set y [expr [winfo rooty $T]+$y]
    GiD_PopupMenu $w $x $y
}

proc ::KCurves::BeginEditCurve { T } {
    set I [$T selection get 0]
    set I [::KCurves::FindMyCurve $T $I]
    set C 0
    set E elemTxtRead
    ::TreeCtrl::FileListEdit $T $I $C $E
}

proc ::KCurves::RefreshTree { {T ""} {onlySave 0} } {
	# ABSTRACT: Refresh the treectrl tree properties 
	# Arguments
	# T         => The tree path
	# onlySave  => Delete the tree and fills it with new items

	global KPriv
	if {$T == ""} {
		set T $KPriv(TreeCurvePath)
	}
	
	# Primero hay que asegurarse de que exista el árbol
	if {[winfo exists $T]} {
	
        set aux [$T item range 0 end]

        foreach item $aux {

            set fullname [DecodeName [$T item tag names $item]]
            set prev "Curves//c."
            append prev $fullname
            set fullname $prev

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
        ::KCurves::DeleteTree $T
        # Fill the model tree
        ::KCurves::FillTreeCurves
        
        
        # Si está disponible añadimos la selección anterior
        catch {
            $T selection add $selectedItem
            $T see [$T selection get]
        }
	}

	# Release the selected item
	set lastSelected {}
}

proc ::KCurves::DoubleClickTree {x y T {item ""}} {
    # Llegamos aqui cuando el usuario hace doble click en el arbol

    global KPriv
    
    #Si no llega directamente el item, miramos cual ha sido pulsado
    if { $item == "" } {  
		set info [$T identify $x $y]
		if { [lindex $info 0] == "item" && [llength $info] >= 4 } {
		    set item [lindex $info 1]
		} else {
		    return ""
		}
	}
    
	if {$item == ""} {
        return ""
    }
	# Conseguimos el path, el id y el valor
    set fullname [::KCurves::GetFullname $item]
   
    if {$fullname == ""} {
        return ""
    }
    
    if {[::KCurves::IsCurve $fullname]} {
        ::KTables::CancelButton
    }
    
	set idFull [string map { "." "" "//" ""} $fullname]	
	set id [::KCurves::setXml $fullname id ]
	set value [::KCurves::setXml $fullname value "" ]

	set dv [::xmlutils::setXml $fullname dv "read" "" "curve"]
	


	set f "$T.f$idFull"
	if {[winfo exists $f]} {
		destroy $f
	}

	set f [::KCurves::buildFrame $T $item] 
	
	#Solo es necesario seguir si se trata de un item editable
	if { $f == "" } {
		return ""
	}
	$T item style set $item C0 styFrame
	$T item element configure $item C0 eWindow -window $f		
	::KCurves::insertIcon $item $T $fullname
	
	#Nos guardamos este item y el valor seleccionado
	set f "$T.f$idFull.cmb"
	set xpath "[::KCurves::setXPath $fullname]"

    if {$xpath != ""} {
        set selCombo [::xmlutils::getComboValue $KPriv(xml) $xpath $f]
        set selComboText [::xmlutils::getComboValue $KPriv(xml) $xpath $f "text"]
        set lastSelected [list $item $selCombo $selComboText]

        set path [DecodeName [$T item tag names $item]]   
        set splitted [::KEGroups::split2 $path //]	
           
        if { [llength $splitted] == 3 } {
            # Si tenemos que construir el arbol de Tabs, entraremos por aqui
            set clase "Tab"
            #::KMat::buildTabFrame $T $item $f $clase
            return ""
        }
    }
}

proc ::KCurves::setXPath { path } {
    
    set splitted [::KMProps::split2 $path //]
    set xpath ""
    set i 0
    foreach itemId $splitted {
        if { $i == 0 } {	  
            set xpath "/Kratos_Data/RootData\[@id='$itemId'\]"		 
            
        } else {
            if { [string index $itemId 0] == "t" } {
            set xpath "$xpath/ContainerTable\[@id='[string range $itemId 2 end]'\]"
            }
            if { [string index $itemId 0] == "i" } {
            set xpath "$xpath/Item\[@id='[string range $itemId 2 end]'\]"
            }
            if { [string index $itemId 0] == "c" } {
            set xpath "$xpath/Container\[@id='[string range $itemId 2 end]'\]"
            }
            if { [string index $itemId 0] == "T" } {
            set xpath "$xpath/TItem\[@id='[string range $itemId 2 end]'\]"
            }
        }		
        incr i
        
    }
    return $xpath
}

proc ::KCurves::buildFrame { T item } {
    # Separamos en funcion de que ventana queramos lanzar
    # Si es Container normal, lanzaremos los clasicos Tab de edicion //de momento no se usa ninguno
    # Si es Points, lanzaremos un tipo especial
    # Si es una curva lanzaremos las ventanas nuevas de edicion
    # Si es un item, haremos el combobox editable
    
	global KPriv
    
    set T $KPriv(TreeCurvePath)
    set fullname [::KCurves::GetFullname $item]
    if {$fullname == ""} {
        return ""
    }
	set idFull [string map { "." "" "//" ""} $fullname]		
	#Comprobamos que sea un item

    if { [::KMProps::itemType $fullname] == "i" } {
        # Si es un item, hay que hacer el combo editable
		set id [::KCurves::setXml $fullname id "" ]
		set pid [::KCurves::setXml $fullname pid "" ]
		set unit [::KCurves::setXml $fullname unit "" ]
		set icon [::KCurves::setXml $fullname icon "" ]
		set state [::KCurves::setXml $fullname state "" ]
		set help [::KCurves::setXml $fullname help "" ]
		set xpath "[::KCurves::setXPath $fullname]"
		set comboList [::xmlutils::getValues $KPriv(xml) $xpath "ivalues"]
		set value [::xmlutils::getValueText $KPriv(xml) $xpath "dv"]

		set bg "#F8F8F8"					
		set f [frame "$T.f$idFull" -borderwidth 0 -background $bg]
		if { [llength $comboList] > 0 } {
		    grid [ttk::combobox $f.cmb -values $comboList -state readonly -textvariable "::KCurves::cmb$idFull"] \
			-row 0 -column 0 -padx 3 -sticky ne -in $f
		    
		    ::xmlutils::setComboValue $::KCurves::xml $xpath $f.cmb $value

		    bind $f.cmb <<ComboboxSelected>> [list ::KCurves::cmbSelectChange $item $T 0]
		} else {
		    grid [ttk::combobox $f.cmb -state normal -textvariable "::KCurves::cmb$idFull"] \
			-row 0 -column 0 -padx 3 -sticky nw -in $f
		    set ::KCurves::cmb$idFull $value
		    bind $f.cmb <FocusOut> [list ::KCurves::cmbSelectChange $item $T 0]
		    bind $f.cmb <Escape> [list ::KCurves::cmbCancel $item $T]
		    bind $f.cmb  <<ComboboxSelected>> [list ::KCurves::cmbSelectChange $item $T 0]
		}
        
		# Si pulsan intro también forzamos la salida del combo
		bind $f.cmb <Return> [list KCurves::cmbSelectChange $item $T 1]
		return $f		  
	} elseif { [::KMProps::itemType $fullname] == "t"} {
		# Si es el conjunto de puntos, lanzamos la ventana de edicion de puntos
        # Por desarrollar
        
	} elseif { [::KCurves::IsCurve $fullname] } {
        # Si es una curva lanzaremos las ventanas de edicion de curva
        
        ::KCurves::CurveBottomFrame $T $item
    }
}

proc ::KCurves::CurveBottomFrame { T item } {
    # ABSTRACT: Create the botton frame
    global KPriv
    
    # If the bottom frame exists, destroy it
    set f [::KCurves::DestroyBottomFrame]
    
    # Create the frame where set the properties
    ttk::frame $f -borderwidth 0
    # Grid for toolbar
    grid $f -row 2 -column 0 -sticky wes
    
    set ::KCurves::mycurve [::KCurves::FindMyCurve $T $item]
	::KCurves::AddButtons
    ::KTables::TableFrameCreateDSets $f "crea"
    return ""
}

proc ::KCurves::AddButtons { } {
    
    set tbf $::KCurves::Buttonframe
    
    ttk::label $tbf.split -text "|" 
	
    ttk::button $tbf.acceptCurve -image [::WinUtils::GetImage ok.gif]  -command [list ::KTables::AcceptButton] -style Toolbutton
    tooltip::tooltip $tbf.acceptCurve [= "Accept Curve"]
    
    ttk::button $tbf.cancelCurve -image [::WinUtils::GetImage tableDeleteAll.gif]  -command [list ::KTables::CancelButton] -style Toolbutton
    tooltip::tooltip $tbf.cancelCurve [= "Cancel Curve"]
    
    grid $tbf.split -sticky e -row 0 -column 3
    grid $tbf.acceptCurve -sticky e -row 0 -column 4
    grid $tbf.cancelCurve -sticky e -row 0 -column 5

}

proc ::KCurves::DestroyButtons {} {
    set tbf $::KCurves::Buttonframe
    
        if {[winfo exists $tbf.split] } {
        destroy $tbf.split
    } 
    
    if {[winfo exists $tbf.acceptCurve] } {
        destroy $tbf.acceptCurve
    }    
    if {[winfo exists $tbf.cancelCurve] } {
        destroy $tbf.cancelCurve
    }
}

proc ::KCurves::DestroyBottomFrame { } {
    # ABSTRACT: Destroy the botton frame
    global KPriv
    variable BottomframePath
    set f ${BottomframePath}
    #set f ${BottomframePath}.bottom

    if {[winfo exists $f]} {
        foreach w [winfo children $f] {
            destroy $w
        }
        destroy $f
    }
    return $f
}

proc ::KCurves::IsCurve { correctedFullname } {
    # True if $correctedFullname is a valid path for a curve at KPriv(xml)
    global KPriv
    set splitted [::KMProps::split2 $correctedFullname / ]
    if { [llength $splitted] eq 3 } {
        return 1
    }
    return 0
}

proc ::KCurves::setXml { path property {value ""} } {
    
    global KPriv
    
    set xpath "[::KCurves::setXPath $path]"
    if {$xpath == ""} {
        return ""
    }
    
    if { $value == "" } {
        set value [$KPriv(xml) set "$xpath/@$property" ]				
        #Cuando hay espacios el xml devuelve una lista y si la imprimes tal cual aparecen corchetes
        if { [llength $value] == 1 } {
            set value [lindex $value 0]
        }
        return $value
    } else {
        $KPriv(xml) set "$xpath/@$property" "$value"
        return "1"
    }
}

proc ::KCurves::cmbSelectChange { item T {remove 1} {selectVal current}  } {   
    global KPriv
    
    set fullname [::KCurves::GetFullname $item]
    
    if {$fullname == ""} {
        return ""
    }
    set idFull [string map { "." "" "//" ""} $fullname]
    
    set id [::KCurves::setXml $fullname id "" ]
    set pid [::KCurves::setXml $fullname pid "" ]
    
    set xpath "[::KCurves::setXPath $fullname]"
    set value [::xmlutils::getValueText $KPriv(xml) $xpath "value"]
		
    
    #Gestión de ivalues
    set xpath "[::KCurves::setXPath $fullname]"
    set comboState [::xmlutils::getComboState $KPriv(xml) $xpath]
    if { $comboState == "normal" } {
	
		set selCombo [set ::KCurves::cmb$idFull]
		set selComboText $selCombo
		if {$selectVal != "current"} {
		    set ::KCurves::lastSelected {}
		}
    } else {
	
		if {$selectVal == "current"} {
		    
		    set f "$T.f${idFull}.cmb"
		    set selCombo [::xmlutils::getComboValue $KPriv(xml) $xpath $f]
		    set selComboText [::xmlutils::getComboValue $KPriv(xml) $xpath $f "text"]
		    set ::KCurves::lastSelected [list $item $selCombo $selComboText]
		    
		} else {
		    
		    set selCombo [lindex $::KFun::lastSelected 1]
		    set selComboText [lindex $::KFun::lastSelected 2]
		    set ::KCurves::lastSelected {}
		}
    }
    
    if { $remove } {
		
		set f "$T.f$idFull"
		if {[winfo exists $f]} {
		    destroy $f
		}
		$T item style set $item C0 styAnyRead
		$T item element configure $item C0 elemTxtRead -text "$pid: $selComboText"
    }
    
    #Guarda el nuevo valor en el xml

    #msg "fullname = $fullname"
    ::KCurves::setXml $fullname dv $selCombo
    
    ::KCurves::insertIcon $item $T $fullname
    
    #Volvemos a cargar el árbol de propiedades

    ::KCurves::RefreshTree $T
    
    
}

proc ::KCurves::DeleteTree { {T ""} } {
    global KPriv

    if { $T == "" } {
	set T $KPriv(TreeCurvePath)
    }
    if { [winfo exists $T] } {
	foreach item [$T item range 0 end] {	  
	    # Elimina el item del árbol
	    catch {$T item delete $item}						
	}
    }		
}

proc ::KCurves::cmbCancel { item T } {
    
    set fullname [::KCurves::GetFullname $item]
    if {$fullname == ""} {
        return ""
    }
    set idFull [string map { "." "" "//" ""} $fullname]
    
    set pid [::KCurves::setXml $fullname pid ""]
    set id [::KCurves::setXml $fullname id ""]
    set value [::KCurves::setXml $fullname value ""]
    
    set xpath "[::KCurves::setXPath $fullname]"
    set valueText [::xmlutils::getValueText $KPriv(xml) $xpath "value"]
    
    set ::KCurves::cmb$idFull "$value"
    set f "$T.f$idFull"
    if {[winfo exists $f]} {
        destroy $f
    }
    $T item style set $item C0 styAnyRead
    $T item element configure $item C0 elemTxtRead -text "$pid: $valueText"
    
    ::KCurves::insertIcon $item $T $fullname
    set ::KCurves::lastSelected ""
    ::KCurves::refreshTree
}

proc ::KCurves::insertIcon { item T fullname {type "mat"}} {
    
	set icon [::KCurves::setXml $fullname icon]
    set imagen [::WinUtils::GetImage $icon]
    if { $imagen != -1 } {
        $T item image $item C0 $imagen
    }			
}


proc ::KCurves::GetFullname { item } {
    if {$item != ""} {
        global KPriv
        set T $KPriv(TreeCurvePath)
        set fullname [DecodeName [$T item tag names $item]]
        set prev "Curves//c."
        append prev $fullname
        return $prev
    }
}


proc ::KCurves::GetCurveInfo { curveid {what "all"} } {
    # Returns a list with the info
    # $curveid is the id of the curve to query on
    # $what can be "all" or a combination of "pid" | "xVar" | "yVar" | "CurveType" | "Points" in a list format.
    # Example:
    # set cinfo [::KCurves::GetCurveInfo "Default"]
    # $cinfo will contain [list "Default Curve" "ByPoints" "Time" "YoungModulus" "Points"]
    # Example: 
    # set props [list "xVar" "Yvar" "pid"]
    # set cinfo [::KCurves::GetCurveInfo "Default" $props]
    # $cinfo will contain [list "ByPoints" "Time" "YoungModulus" "Default Curve"]
    
    global KPriv
    
    if {$what == ""} {return ""}
    
    set retlist ""
    
    set curve [::KCurves::GetCurves id $curveid]
    set cpid [$curve getAttribute pid ""]

    set curveprops [$curve childNodes]
    #msg "what $what"
    foreach it $what {
        
        set it [::KCurves::PropertyIdentify $it]
        #msg $it
        if { [string equal $it "pid"] || $what eq "all" } {
            set retlist [lappend retlist $cpid]
        }
        if { $it eq "xVar" || $what eq "all" } {
            set retlist [lappend retlist [[lindex $curveprops 1] getAttribute dv ""]]
        }
        if { $it eq "yVar" || $what eq "all" } {
            set retlist [lappend retlist [[lindex $curveprops 2] getAttribute dv ""]]
        }
        if { $it eq "CurveType" || $what eq "all" } {
            set retlist [lappend retlist [[lindex $curveprops 0] getAttribute dv ""]]
        }
        if { $it eq "Points" || $what eq "all" } {
            set retlist [lappend retlist [::KCurves::GetPointsfromCurveNode $curve]]
        }
    }
    return $retlist
}

proc ::KCurves::PruebaCurvas { } {

    # Ejemplo de uso de Get Curves y GetCurveInfo
    set curve [::KCurves::GetCurves xVar "Time"]
    # msg $curve
    set Timemore [::KCurves::GetCurves yVar "More" $curve "namelist"]
    # msg "Time and more $Timemore"
    set tminfo [::KCurves::GetCurveInfo $Timemore "all"]
    # msg "allinfo $tminfo"
    set xinfo [::KCurves::GetCurveInfo $Timemore "xvar yvar"]
    # msg "xinfo $xinfo"
}

proc ::KCurves::GetCurves { property value {searchon "xml" } {retType "nodelist" } } {
    # Returns a list of nodes according to user specifications at property & value
    # property can be "name" "xvar" "yvar" 
    # value is the queried value for the property
    # searchon can be "xml" to look for all curves in the xml, or a node list
    # retType can be "nodelist" or "namelist", specifies de return data format
    # "nodelist" will return a list of nodes. "namelist" will return a list of names
    # Example: To get all the curves which xValue is Time: set curves [::KCurves::GetCurves xval time]
    # Example: To get the curve called "Young Module on Steel" from the curves list of the previous example: 
    # set c [::KCurves::GetCurves name "Young Module on Steel" $curves
    
    # First identify the property

    set prop [::KCurves::PropertyIdentify $property]

    if { $prop == "" } { return "" }
    
    set finalnodes ""
    
    if { $searchon eq "xml" } {
        # Get the nodes of all the curves at xml
        set curvesBaseNode [::KCurves::CurveBase]
        set nodes [$curvesBaseNode childNodes]
    } else {
        set nodes $searchon
    }
    
    # Get just the curves with the user specifications
    if { $prop eq "pid" } {
        foreach node $nodes {
            set cpid [$node getAttribute pid ""]
            if { $cpid eq $value } {
                set finalnodes [lappend finalnodes $node]
            }
        }
    } elseif { $prop eq "id" } {
        foreach node $nodes {
            set cid [$node getAttribute id ""]
            if { $cid eq $value } {
                set finalnodes [lappend finalnodes $node]
            }
        }
    } else {
        foreach node $nodes {
            set nodeprops [$node childNodes]
            if { $prop eq "xVar" } {
                set subnode [lindex $nodeprops 1]
            } elseif { $prop eq "yVar" } {
                set subnode [lindex $nodeprops 2]
            } elseif { $prop eq "CurveType" } {
                set subnode [lindex $nodeprops 0]
            }
            set val [$subnode getAttribute dv ""]
            if { $val eq $value } {
                set finalnodes [lappend finalnodes $node]
            }
        }
    }
    if { $retType eq "nodelist" } {
        return $finalnodes 
    } elseif {$retType eq "namelist"} {
        set names ""
        foreach node $finalnodes {
            set names [lappend names [$node getAttribute id ""]]
        }
        return $names
    }
}

proc ::KCurves::GetPointsfromCurveNode { node } {
    set curveprops [$node childNodes]
    set pointnodes [[lindex $curveprops 3] childNodes]
    set retlist ""
    foreach point $pointnodes {
        set retlist [lappend retlist [$point getAttribute Xval ""] [$point getAttribute Yval ""]]
    }
    return $retlist
}

proc ::KCurves::PropertyIdentify { prop } {

    set idlist [list "name" "Name" "NAME" "ID" "id" "Id"]
    set namelist [list "pname" "pName" "PNAME" "pid" "Pid" "PID"]
    set xlist [list "x" "xval" "X" "xVal" "Xval" "XVAL" "xvar" "xVar"]
    set ylist [list "y" "yval" "Y" "yVal" "Yval" "YVAL" "yvar" "yVar"]
    set ctlist [list "curve" "Curve" "type" "Type" "CurveType" "CURVETYPE"]
    set pointlist [list "points" "Points" "POINTS" ]
    
    if { $prop in $idlist } {
        return "id"
    }
    if { $prop in $namelist } {
        return "pid"
    }
    if { $prop in $xlist } {
        return "xVar"
    }
    if { $prop in $ylist } {
        return "yVar"
    }
    if { $prop in $ctlist } {
        return "CurveType"
    }
    if { $prop eq "all" } {
        return "all"
    }
    if { $prop in $pointlist } {
        return "Points"
    }
    return ""
    
}


proc ::KCurves::CurveBase { } {
    global KPriv
    
    set baseNodePath "/Kratos_Data/RootData\[@id='Curves'\]"
    set node [$KPriv(xml) selectNodes $baseNodePath]
    return $node
}


proc ::KCurves::SetCurveProp { what curve property value } {
   
    # Modifies in xml the curve which name is $curve, new $property 's value will be $value.
    # $what can be "name" or "item", referrin to "pid" or the item of the tree.
    # $curve will be the pid of the curve [if "name"], or the item of the tree [if "item".
    # $property can be "pid", "CurveType", "xVar", "yVar".
    # $value is the new value.
    
    global KPriv
    set T $KPriv(TreeCurvePath)
    # First identify the property
    set prop [::KCurves::PropertyIdentify $property]
    if { $prop == "" } { return "" }

    # To simplify procedures, we'll transform item to pid
    if {$what eq "item"} {
        set curveitem [::KCurves::FindMyCurve $T $curve]
        set curve [$T item text $curveitem 0]
    }
    # Now, all curves are identified by pid
    
    # Get the nodes of all the curves at xml
    set curvesBaseNode [::KCurves::CurveBase]
    set nodes [$curvesBaseNode childNodes]
    
    # Set just the new values just to the curve with the user specifications
    foreach node $nodes {
        set cpid [$node getAttribute pid ""]
        if { $cpid eq $curve } {
            set path [::xmlutils::getFullnameFromNode $node]
            ::KCurves::setXml $path $prop $value
            break
        }
    }
    ::KCurves::RefreshTree
}

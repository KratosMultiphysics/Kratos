###############################################################################
#
#	NAME: kmprop.tcl
#
#	PURPOSE: Main window to manage model properties
#
#	QUANTECH ATZ-DEVELOPMENT DEPARTMENT
#
#	AUTHOR : G. Socorro
#
#	CREATED AT: 25/02/2010
#
#	HISTORY:
#
#       0.8- 22/06/11-G. Socorro, delete snit, tdom and xmlstruct from the package require
#       0.7- 07/06/11 GS, add composite and plastic material to the structural analysis application
#	0.6- 27/09/10 LC, Correct bugs inserting new materials and double click in containers
#	0.5- 08/06/10 KS, New materials database structure + template for add material.
#	0.4- 11/05/10 KS, varius bugs fixed.
#	0.3- 26/04/10 KS, Material Groups.
#	0.2- 20/04/10 KS, Read, Add, Remove, Rename OK.
#	0.1- 25/02/10 LCA, create a base source code from the kegroups.tcl script
#
###############################################################################

package require treectrl
package require tooltip
package provide KMat 1.0 

# Create a base namespace KMat
namespace eval ::KMat:: {
    
    # Path of the base window 
    variable WinPath ".gid.kmprops"
    # Path of the Materials Tree
    variable TreeMatPath 
    
    variable WinLayout 
    variable SystemHighlight
    variable SystemHighlightText
    
    variable lastSelected {}
    variable abdlist
    variable xml ""
    
    # Se inicializan las clases dinámicamente leyendo del xml
    variable visibilityVars {}

}

proc ::KMat::Init {} {
    
    variable WinLayout;	variable SystemHighlight
    variable SystemHighlightText; variable abdlist
    
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
    
    set abdlist [list A11 A12 A13 A21 A22 A23 A31 A32 A33 \
		     B11 B12 B13 B21 B22 B23 B31 B32 B33 \
		     D11 D12 D13 D21 D22 D23 D31 D32 D33 ]
    
    global KPriv
    if {[info exists KPriv(xmlMat)]} {
	set ::KMat::xml $KPriv(xmlMat)
    }
}

# Init KMat namespace
::KMat::Init

proc ::KMat::CreateTreeAndToolbar { w } {
    
    # Create the treectrl properties
    set mdf [ttk::frame $w.middle ]
    set T [::KMat::CreateTreeProperties $w]
    
    grid $w.middle -sticky wens
    
    # Create the frame where set the properties
    set f [ttk::frame $w.fBottom -borderwidth 0 ]
    
    # Grid for toolbar
    grid $f -row 2 -column 0 -sticky wes
    
    # Create the toolbar frame
    set tbf [ttk::frame $w.tbar -borderwidth 0]
    
    ttk::button $tbf.newmaterial -image [::WinUtils::GetImage add.gif]  -command [list ::KMat::CreateNewMaterial $T] -style Toolbutton
    tooltip::tooltip $tbf.newmaterial [= "Create new material"]

    ttk::button $tbf.deleteMaterialsId -image [::WinUtils::GetImage delete_icon.gif]  -command [list ::KMat::DeleteMaterial $T] -style Toolbutton
    tooltip::tooltip $tbf.deleteMaterialsId [= "Delete material"]

    ttk::button $tbf.refreshTree -image [::WinUtils::GetImage refresh.gif]  -command [list ::KMat::refreshTree] -style Toolbutton
    tooltip::tooltip $tbf.refreshTree [= "Refresh Tree"]

    #	ttk::button $tbf.testbutton -image [::WinUtils::GetImage uc.gif]  -command [list ::KMat::testbutton] -style Toolbutton
    #	tooltip::tooltip $tbf.testbutton [= "tests"]
    
    
    #set tf [ttk::frame $w.close] 
    grid $tbf -sticky ews
    grid anchor $tbf w
    #grid [ttk::button $tf.bClose -text [= "Close"] -command [list destroy $w]]  -sticky ew -padx 5 -pady 3
    
    # Grid for toolbar 
    #grid $w.tbar -row 2 -column 0 -sticky wes
    grid $tbf.newmaterial -sticky we -row 0 -column 0
    grid $tbf.deleteMaterialsId -sticky we -row 0 -column 1
    grid $tbf.refreshTree -sticky we -row 0 -column 2
    #	grid $tbf.testbutton -sticky we -row 0 -column 3

    #grid rowconfigure $w 2 -weight 1
    #grid columnconfigure $w 0 -weight 1
    
    focus $T
    return $T
}

proc ::KMat::CreateTreeProperties {w} {

    variable SystemHighlight
    variable SystemHighlightText

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
    
    $T column create -text "Materials" -tags C0 -weight 0
    
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
    #$T notify install <Drag-receive>
    $T notify install <Edit-accept>

    # Notify bind
    # TODO

    #		$T notify bind DragTag <Drag-receive> { ::KMProps::ReceiveDragGroups  %T %l %I }
    $T notify bind EditTag <Edit-accept> { ::KMat::SetMatToRename %T %I %t }

    bind $T <Button-1> [list ::KMat::ClickTree %x %y $T]
    bind $T <Double-Button-1> [list ::KMat::DoubleClickTree %x %y $T]
    #	 bind $T <Return> [list SetLayersTo TOUSE $T]
    #	 bind $T <Key-Delete> [list SetLayersToDelete $T]
    #	 bind $T <Alt_L> [list InvertSelectionTableList $T]
    #	 bind $T <Alt_R> [list InvertSelectionTableList $T]
    #	 bind $T <Meta_L> [list InvertSelectionTableList $T]
    #	 bind $T <Meta_R> [list InvertSelectionTableList $T]
    bind $T <F2> [list ::KMat::BeginEditMaterial $T]

    bind $T <Button-3> "[list ::KMat::MenuContextualGroup %W %x %y] ; break"

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

#---------------------------------------------------------------------------------------------- 
# Lee el xml y carga el árbol de propiedades de forma iterativa como máximo hasta 7 niveles
#----------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------- 
# Lee el xml y carga el árbol de propiedades de forma iterativa como máximo hasta 7 niveles
#----------------------------------------------------------------------------------------------
proc ::KMat::FillTreeMat { } {
    
    variable abdlist
    global KPriv
    set KPriv(materialsId) {}
    
    set T $::KMat::TreeMatPath
    set nodes [$::KMat::xml selectNodes "/Kratos_KMat_DB/Materials/MaterialGroup\[@id\]"]
    
    foreach node $nodes {			  
	#Nos guardamos todos los Id	
	#		set item [::KMat::InsertNewItem [$node getAttribute pid ""] [$node getAttribute id ""] $T "" "root" [$node hasChildNodes] "normal" [$node getAttribute open "0"]]
	set item [::KMat::InsertNewItem [$node getAttribute pid ""] [$node getAttribute id ""] $T "" "root" [$node hasChildNodes] [::KMat::stateNode $node] [$node getAttribute open "0"]]
	
	set nodes2 [$node childNodes]
	foreach node2 $nodes2 {
	    #			set item2 [::KMat::InsertNewItem [$node2 getAttribute pid ""] [::KMat::splitNode $node2] $T "[$node getAttribute id ""]//" "$item" [$node2 hasChildNodes] "normal" [$node2 getAttribute open "0"]]
	    set item2 [::KMat::InsertNewItem [$node2 getAttribute pid ""] [::KMat::splitNode $node2] $T "[$node getAttribute id ""]//" "$item" [$node2 hasChildNodes] [::KMat::stateNode $node2] [$node2 getAttribute open "0"]]
	    if {$item2 != -1} {						
		lappend KPriv(materialsId) [$node2 getAttribute pid ""]
		#Seleccionamos los hijos (3º nivel)
		set nodes3 [$node2 childNodes]		   
		foreach node3 $nodes3 {
		    if { [$node3 getAttribute pid ""] ni $abdlist } {
			# set item3 [::KMat::InsertNewItem [$node3 getAttribute pid ""] [::KMat::splitNode $node3] $T "[$node getAttribute id ""]//[::KMat::splitNode $node2]//" "$item2" [$node3 hasChildNodes] "normal" [$node3 getAttribute open "0"]]
			set item3 [::KMat::InsertNewItem [$node3 getAttribute pid ""] [::KMat::splitNode $node3] $T "[$node getAttribute id ""]//[::KMat::splitNode $node2]//" "$item2" [$node3 hasChildNodes] [::KMat::stateNode $node3] [$node3 getAttribute open "0"]]
		    }
		    if {$item3 != -1} {
			set nodes4 [$node3 childNodes]
			foreach node4 $nodes4 {
			    # set item4 [::KMat::InsertNewItem [$node4 getAttribute pid ""] [::KMat::splitNode $node4] $T "[$node getAttribute id ""]//[::KMat::splitNode $node2]//[::KMat::splitNode $node3]//" "$item3" [$node4 hasChildNodes] "normal" [$node4 getAttribute open "0"]]
			    set item4 [::KMat::InsertNewItem [$node4 getAttribute pid ""] [::KMat::splitNode $node4] $T "[$node getAttribute id ""]//[::KMat::splitNode $node2]//[::KMat::splitNode $node3]//" "$item3" [$node4 hasChildNodes] [::KMat::stateNode $node4] [$node4 getAttribute open "0"]]
			    if {$item4 != -1} {
				set nodes5 [$node4 childNodes]
				foreach node5 $nodes5 {
				    # set item5 [::KMat::InsertNewItem [$node5 getAttribute pid ""] [::KMat::splitNode $node5] $T "[$node getAttribute id ""]//[::KMat::splitNode $node2]//[::KMat::splitNode $node3]//[::KMat::splitNode $node4]//" "$item4" [$node5 hasChildNodes]  "normal" [$node5 getAttribute open "0"]]
				    set item5 [::KMat::InsertNewItem [$node5 getAttribute pid ""] [::KMat::splitNode $node5] $T "[$node getAttribute id ""]//[::KMat::splitNode $node2]//[::KMat::splitNode $node3]//[::KMat::splitNode $node4]//" "$item4" [$node5 hasChildNodes]  [::KMat::stateNode $node5] [$node5 getAttribute open "0"]]
				    if {$item5 != -1} {
					set nodes6 [$node5 childNodes]				   
					foreach node6 $nodes6 {		
					    set item6 [::KMat::InsertNewItem [$node6 getAttribute pid ""] [::KMat::splitNode $node6] $T "[$node getAttribute id ""]//[::KMat::splitNode $node2]//[::KMat::splitNode $node3]//[::KMat::splitNode $node4]//[::KMat::splitNode $node5]//" "$item5" [$node6 hasChildNodes]]
					    #		set item6 [::KMat::InsertNewItem $node6 [::KMat::splitNode $node6] $T "[$node getAttribute id 0]//[::KMat::splitNode $node2]//[::KMat::splitNode $node3]//[::KMat::splitNode $node4]//[::KMat::splitNode $node5]//" "$item5" [$node6 hasChildNodes] [::KMat::stateNode $node6] [$node6 getAttribute open "0"]]
					    if {$item6 != -1} {
						set nodes7 [$node6 childNodes]			 
						foreach node7 $nodes7 {			
						    set item7 [::KMat::InsertNewItem [$node7 getAttribute pid ""] [::KMat::splitNode $node7] $T "[$node getAttribute id ""]//[::KMat::splitNode $node2]//[::KMat::splitNode $node3]//[::KMat::splitNode $node4]//[::KMat::splitNode $node5]//[::KMat::splitNode $node6]//" "$item6" [$node7 hasChildNodes]]
						    #				set item7 [::KMat::InsertNewItem $node7 [::KMat::splitNode $node7] $T "[$node getAttribute id 0]//[::KMat::splitNode $node2]//[::KMat::splitNode $node3]//[::KMat::splitNode $node4]//[::KMat::splitNode $node5]//[::KMat::splitNode $node6]//" "$item6" [$node7 hasChildNodes] [::KMat::stateNode $node7] [$node7 getAttribute open "0"]]

						}
					    }
					}
				    }
				}
			    }
			}
		    }
		}
	    }
	}		
    }

}


proc ::KMat::DoubleClickTree {x y T {item ""}} {
    
    variable lastSelected
    
    #Si no llega directamente el item, miramos cual ha sido pulsado
    if { $item == "" } {  
	set info [$T identify $x $y]
	if { [lindex $info 0] == "item" && [llength $info] >= 4 } {
	    set item [lindex $info 1]
	} else {
	    return ""
	}
    }
    
    set fullname [DecodeName [$T item tag names $item]]
    
    set idFull [string map { "." "" "//" ""} $fullname]	
    set id [::KMat::setXml $fullname id ]
    set value [::KMat::setXml $fullname value "" ]

    if { $value != "" && [llength $lastSelected] > 0 && [lindex $lastSelected 0] > $item && $id != "ABD"} {
	
	#catch
	::KMat::cmbSelectChange [lindex $lastSelected 0] $T 1 "anterior"
    }
    set f "$T.f$idFull"
    if {[winfo exists $f]} {
	destroy $f
    }
    
    set f [::KMat::buildFrame $T $item]
    
    #Solo es necesario seguir si se trata de un item editable
    if { $f == "" } {
	#msg [$::KMat::xml asXML]
	return ""
    }
    
    $T item style set $item C0 styFrame
    $T item element configure $item C0 eWindow -window $f		
    ::KMat::insertIcon $item $T $fullname
    
    #Nos guardamos este item y el valor seleccionado
    set f "$T.f$idFull.cmb"
    set xpath "[::KMat::setXPath $fullname]"
    set selCombo [::xmlutils::getComboValue $::KMat::xml $xpath $f]
    set selComboText [::xmlutils::getComboValue $::KMat::xml $xpath $f "text"]
    set lastSelected [list $item $selCombo $selComboText]

    #	if {$id == "ABD"} {
    #		WarnWin [_ "Lets make a Tab for ABD...!!!"]
    #		set clase "Tab"
    #		::KMat::buildTabFrameABD $T $item $clase
    #		return ""
    #	}

    set path [DecodeName [$T item tag names $item]]   
    set splitted [::KEGroups::split2 $path //]			   
    if { [llength $splitted] == 2 } {
	set clase "Tab"
	#		::KMat::buildTabFrame $T $item $clase
	return ""
    }


}

proc ::KMat::ClickTree { x y T } {
    
    variable lastSelected
    set info [$T identify $x $y]
    
    if { [lindex $info 0] == "item" && [llength $info] >= 4 } {
	set item [lindex $info 1]
	set col [lindex $info 3]
	set fullname [DecodeName [$T item tag names $item]]
	
	#test getMaterialProperties
	#set matlist [::KMat::getMaterialProperties "p" $fullname]
	
	set id [::KMat::setXml $fullname id "" ]
	
	#Eliminamos el anterior combo, si aun está visible
	if {[llength $lastSelected] > 0 && [lindex $lastSelected 0] != $item} { 
	    
	    ::KMat::cmbSelectChange [lindex $lastSelected 0] $T 1 "anterior"
	}
	
	#Si pinchan en un item con hijos lo abrimos
	#				$T item toggle $item
	
    } elseif { [lindex $info 0] == "header" && [lindex $info 1] == "0" } {
	if { [$T column cget C0 -arrow] == "up" } {
	    $T column configure C0 -arrow down
	    $T item sort 0 -dictionary -increasing
	} else {
	    $T column configure C0 -arrow up
	    $T item sort 0 -dictionary -decreasing
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
	#SetLayersTo TOUSE $T
    } elseif { $col == 1 } {		
	set parent [winfo parent $T]
    }
    
    if { $col != 0 } {
	return -code break
    }
    return ""
}

proc ::KMat::CreateNewMaterial { {T ""} {name ""} } {

    #	if { $name == "" } {
    #				set name [::KMat::GetAutomaticMatName]
    #	} else {
    #				if { ![::KEGroups::isValidGroupName $name] } {
    #					WarnWin [_ "Bad material name, start or end by '//' is not allowed"]
    #					return ""
    #				}
    #	}
    
    # insertamos en el grup adecuado. Si no, no insertamos
    set item [$T selection get 0]	
    if {$item == "" } {
	WarnWin [_ "No materials group selected."]
	return ""
    }

    #Insertamos siempre el nuevo material en el padre del item seleccionado
    set aux 99
    set aux2 [$T item parent $item]
    while {$aux != 1 && $aux2 != 0} {
	set padre [$T item parent $item]
	set item $padre		
	set path [DecodeName [$T item tag names $item]]   
	set splitted [::KEGroups::split2 $path //]
	set aux [llength $splitted]
    }

    #		if { [llength $splitted] != 1 } {
    #		WarnWin [_ "Can not create new material. No material Group selected."]
    #				return ""
    #		}

    set path [DecodeName [$T item tag names $item]]   

    if { $name == "" } {
	set name [::KMat::GetAutomaticMatName "" $path]
    } else {
	if { ![::KEGroups::isValidGroupName $name] } {
	    WarnWin [_ "Bad material name, start or end by '//' is not allowed"]
	    return ""
	}
    }



    ::KMat::insertXml "$path" $name 1 Generic
    
    ::KMat::refreshTree $T
    
    # añadimos el nuevo material en el arbol
    #::KMat::InsertNewMaterial $name $T 1 Generic "" "$item" 0  

    return $name   
}

proc ::KMat::DeleteMaterial { {T ""} {name ""} } {
    
    global KPriv
    
    set items [$T selection get]

    if { $items == "" } {
	WarnWin [_ "No material selected."]
	return ""
    }

    #Obtenemos el path del nodo
    set path [DecodeName [$T item tag names $items]]   

    # si no estamos en el segundo nivel (materiales) no borramos nada
    set splitted [::KEGroups::split2 $path //]		
    
    #		set splitted [::KEGroups::split2 $path //]
    #		set aux [llength $splitted]

    if { [llength $splitted] == 1 } {
	WarnWin [_ "No material selected. Can not delete a material group."]
	return ""
    }
    if { [llength $splitted] == 3 } {
	set padre [$T item parent $items]
	set items $padre		
	set path [DecodeName [$T item tag names $items]]   
	set splitted [::KEGroups::split2 $path //]		
    }

    if { [llength $splitted] != 2 } {
	WarnWin [_ "Error deleting Material."]
	return ""		
    }

    if {[llength $items] > 0 } {

	# Para avisar de qué items se van a borrar
	set tuttoItem {}
	foreach it $items {
	    lappend tuttoItem "[lindex [$T item text $it] 0]"
	}
	set aviso "Are you sure you want to delele $tuttoItem ?"
	set confirmado [::WinUtils::confirmBox "." "$aviso"]
	if { $confirmado == "ok" } {

	    #Buscamos todos los descendientes del item a eliminar
	    set completList {}
	    foreach item $items {						
		lappend completList $item
		foreach i [$T item descendants $item] {
		    set idx [lsearch -exact $items $i]
		    if { $idx == -1 } {
			lappend completList $i		
		    }
		}
	    }
	    for {set i 0} { $i < [llength $completList] } {incr i} {
		set MaterialId [$T item text [lindex $completList $i] 0]
		set KPriv(materialsId) [::KEGroups::listReplace $KPriv(materialsId) $MaterialId]
	    }
	    
	    for {set i 0} { $i < [llength $items] } {incr i} {
		#Elimina el grupo del xml
		catch { ::xmlutils::unsetXml [DecodeName [$T item tag names [lindex $items $i]]] "mat" }			   
		
		#Resetea todas las ocurrencias del material en las propiedades del .spd
		::KMProps::chekMaterials [$T item text $item 0]
		
		#Elimina el grupo del árbol
		::KMProps::deleteItem $T $item
		
	    }
	} else {
	    #WarnWin [= "No material selected"]
	}
    }
}


proc ::KMat::DeleteTree { {T ""} } {

    if { $T == "" } {
	global KPriv
	set T $::KMat::TreeMatPath
    }
    if { [winfo exists $::KMProps::WinPath] } {
	foreach item [$T item range 0 end] {	  
	    #Elimina el item del árbol
	    catch {$T item delete $item}						
	}
    }		
}


proc ::KMat::MenuContextualGroup { T x y } {
    
    set w $T.menucontextualgroup
    if { [winfo exists $w] } {
	destroy $w
    }
    menu $w
    
    $w add command -label [= "New Material"] -command [list ::KMat::CreateNewMaterial $T]
    $w add command -label [= "Delete Material"] -command [list ::KMat::DeleteMaterial $T]
    
    set item [$T selection get 0]	
    if {$item != "" } {		
	set path [DecodeName [$T item tag names $item]]   
	set splitted [::KEGroups::split2 $path //]		
	if { [llength $splitted] == 2 } {
	    $w add command -label [= "Rename material"] -command [list ::KMat::BeginEditMaterial $T]				
	    $w add command -label [= "Copy material"] -command [list ::KMat::CopyMaterial $T]				
	} else {
	    $w add command -label [= "Rename material"] -state disabled				
	    $w add command -label [= "Copy material"] -state disabled				
	}
    }
    
    
    #		$w add command -label [= "Collapse All"] -command [list $T collapse -recurse "$item"] -state normal
    #		$w add command -label [= "Expand All"] -command [list $T expand -recurse "$item"] -state normal
    
    set x [expr [winfo rootx $T]+$x+2]
    set y [expr [winfo rooty $T]+$y]
    GiD_PopupMenu $w $x $y
}

proc ::KMat::insertXml { path id state type } {
    global KPriv

    if { $path == "root" } { 
	set xpath "/Kratos_KMat_DB/Materials"
    } else {
	set xpath "[::KMat::setMatXPath $path]"
    }

    set templatePath "/Kratos_KMat_DB/Templates/Template"
    if { $path == "Metal" } {
	set maticon "grey.gif"
    } elseif { $path == "Fluid" } {
	set maticon "blue.gif"		
    } elseif { $path == "Plastic" } {
	set maticon "red.gif"		
    } elseif { $path == "Composite" } {
	set maticon "green.gif"		
    }

    set attributesArray [ list id=\"$id\" pid=\"$id\" icon=\"$maticon\" help=\"$id\" open=\"1\"]

    ::xmlutils::copyTemplate $::KMat::xml $xpath $templatePath "NewMaterial" "Material" $attributesArray 

    set xmlArray [::xmlutils::replaceTemplate $::KMat::xml $xpath]
    
    set KPriv(xmlDocMat) [lindex $xmlArray 0]
    set KPriv(xmlMat) [lindex $xmlArray 1]
    
    set ::KMat::xml $KPriv(xmlMat)
}

proc ::KMat::insertXmlCopy { path id state type sourcematname} {
    global KPriv

    if { $path == "root" } { 
	set xpath "/Kratos_KMat_DB/Materials"
    } else {
	set xpath "[::KMat::setMatXPath $path]"
    }

    #	set templatePath "/Kratos_KMat_DB/Templates/Template"
    #	set templatePath "/Kratos_KMat_DB/Materials/Material"
    set templatePath "$xpath/Material"
    if { $path == "Metal" } {
	set maticon "grey.gif"
    } elseif { $path == "Fluid" } {
	set maticon "blue.gif"		
    } elseif { $path == "Plastic" } {
	set maticon "red.gif"		
    } elseif { $path == "Composite" } {
	set maticon "green.gif"		
    }

    set attributesArray [ list id=\"$id\" pid=\"$id\" icon=\"$maticon\" help=\"$id\" open=\"1\"]

    #	::xmlutils::copyTemplate  $::KMat::xml $xpath $templatePath "NewMaterial" "Material" $attributesArray 
    ::xmlutils::copyTemplate  $::KMat::xml $xpath $templatePath "$sourcematname" "Material" $attributesArray 

    set xmlArray [::xmlutils::replaceTemplate $::KMat::xml $xpath]
    
    set KPriv(xmlDocMat) [lindex $xmlArray 0]
    set KPriv(xmlMat) [lindex $xmlArray 1]
    
    set ::KMat::xml $KPriv(xmlMat)
}


# Prepara la query para utilizar las funciones de domNOde
proc ::KMat::setMatXPath { path } {
    set splitted [::KEGroups::split2 $path //]		
    if { [llength $splitted] >= 1 } {		
	#				set xpath "/Kratos_KMat_DB/Materials/Material\[@id='[lindex $splitted 0]'\]"
	set xpath "/Kratos_KMat_DB/Materials/MaterialGroup\[@id='[lindex $splitted 0]'\]"
    } 
    if { [llength $splitted] >= 2 } {		
	set xpath "$xpath/Group\[@id='[lindex $splitted 1]'\]"
    } 
    if { [llength $splitted] >= 3 } {		
	set xpath "$xpath/Group\[@id='[lindex $splitted 2]'\]"
    }
    if { [llength $splitted] >= 4 } {		
	set xpath "$xpath/Group\[@id='[lindex $splitted 3]'\]"
    }
    if { [llength $splitted] >= 5 } {		
	set xpath "$xpath/Group\[@id='[lindex $splitted 4]'\]"
    }
    return $xpath
}



proc ::KMat::InsertNewMaterial { MatName T {state 1} {type "Generic"} {parent ""} {parentitem root} {childs true} } {
    
    if { $parent != "" } {
	set fullname $parent$MatName
    } else {
	set fullname $MatName
    }
    
    set item [$T item create -button yes -tags [EncodeName $fullname]]
    $T item lastchild $parentitem $item
    $T item style set $item C0 styAnyRead
    $T item element configure $item C0 elemTxtRead -text "$MatName"

    set item2 [::KMat::InsertNewItem "Density" "p.Density" $T "$fullname" "$item" 0]	
    set item3 [::KMat::InsertNewItem "Young_Modulus" "p.Young_Modulus" $T "$fullname" "$item" 0]	
    set item4 [::KMat::InsertNewItem "Poisson_Ratio" "p.Poisson_Ratio" $T "$fullname" "$item" 0]	
    set item5 [::KMat::InsertNewItem "Viscosity" "p.Viscosity" $T "$fullname" "$item" 0]	
    set item6 [::KMat::InsertNewItem "Description" "p.Description" $T "$fullname" "$item" 0]	

    set item7 [::KMat::InsertNewItem "Description2" "c.General//p.Description2" $T "$fullname" "$item" 0]
    set item8 [::KMat::InsertNewItem "Description" "p.Description" $T "$fullname" "$item" 0]
    set item8 [::KMat::InsertNewItem "Density" "p.Density" $T "$fullname" "$item" 0]
    
    ::KMat::refreshTree $T
    return $item
}


proc ::KMat::BeginEditMaterial { T } {
    set I [$T item id active]
    set C 0
    set E elemTxtRead

    set path [DecodeName [$T item tag names $I]]   
    set splitted [::KEGroups::split2 $path //]		
    if { [llength $splitted] == 2 } {
	::TreeCtrl::FileListEdit $T $I $C $E
    }
    
}

proc ::KMat::CopyMaterial { {T ""} {name ""} } {
    
    #	if { $name == "" } {
    #				set name [::KMat::GetAutomaticMatName]
    #	} else {
    #				if { ![::KEGroups::isValidGroupName $name] } {
    #					WarnWin [_ "Bad material name, start or end by '//' is not allowed"]
    #					return ""
    #				}
    #	}
    
    # insertamos en el grup adecuado. Si no, no insertamos
    set item [$T selection get 0]	
    if {$item == "" } {
	WarnWin [_ "No materials group selected."]
	return ""
    }

    set fullname [DecodeName [$T item tag names $item]]
    set sourcematname [::KMat::setXml $fullname id "" ]

    #Insertamos siempre el nuevo material en el padre del item seleccionado
    set aux 99
    set aux2 [$T item parent $item]
    while {$aux != 1 && $aux2 != 0} {
	set padre [$T item parent $item]
	set item $padre		
	set path [DecodeName [$T item tag names $item]]   
	set splitted [::KEGroups::split2 $path //]
	set aux [llength $splitted]
    }

    set path [DecodeName [$T item tag names $item]]   

    if { $name == "" } {
	set name [::KMat::GetAutomaticMatName "" $sourcematname]
    } else {
	if { ![::KEGroups::isValidGroupName $name] } {
	    WarnWin [_ "Bad material name, start or end by '//' is not allowed"]
	    return ""
	}
    }

    ::KMat::insertXmlCopy "$path" $name 1 Generic $sourcematname
    
    
    ::KMat::refreshTree $T
    
    # añadimos el nuevo material en el arbol
    #::KMat::InsertNewMaterial $name $T 1 Generic "" "$item" 0  
    
    
    return $name   
    
    
}


proc ::KMat::SetMatToRename { T item newtext } {
    
    global KPriv
    if { $newtext == ""} {  return  }
    if { $item == 0 } {
	WarnWin [_ "Root folder can't be edited"]
	return
    }
    set oldId [$T item text $item 0]
    #Controlamos q el nombre no esté ya en el árbol (a no ser q no se haya cambiado)
    if { $oldId != $newtext && $newtext in $KPriv(materialsId) } {
	WarnWin [_ "The group name '$s' already exist.\nChoose another, please." $newtext]
	return
    }
    #Validamos q el nombre no tenga carácteres que vulneran la seguridad y quitamos espacios
    set newtext [::KUtils::parseTreeStr $newtext]
    if { $newtext == -1 } {
	WarnWin [_ "You can't use some reservate chars like:\n  :   /   $   .   \\  %  "]
	return
    }
    
    
    #Guardamos el nivel del item a renombrar
    set fullname [DecodeName [$T item tag names $item]]
    set splitted [::KEGroups::split2 $fullname //]
    set whereRename [llength $splitted]
    
    #Recorremos toda su familia		cambiando cada path
    set items [$T item descendants $item]
    
    foreach i $items {		
	set fullname [DecodeName [$T item tag names $i]]
	set splitted [::KEGroups::split2 $fullname //]
	
	#Sustituimos la posición apropiada del path por el nuevo nombre
	set splitted [::KEGroups::listReplace $splitted [lindex $splitted [expr $whereRename - 1]] m.$newtext]
	
	#Reconstruimos el nuevo path
	set fullNewName ""
	foreach iSplit $splitted {		
	    set fullNewName "$fullNewName$iSplit//"
	}
	set fullNewName [string range $fullNewName 0 end-2]
	
	#Cambiamos la etiqueta (el path) en el arbol
	$T item tag remove $i [list names [$T item tag names $i]]
	$T item tag add $i [EncodeName $fullNewName]
    }
    
    #Cambiamos el item a renombrar 
    set fullname [DecodeName [$T item tag names $item]]
    
    ::KMat::editTag $T $item $fullname $newtext
    
    #Renombramos también todas las ocurrencias del material en el .spd de propiedades
    ::KMProps::chekMaterials $oldId $newtext
    
    
    return ""
}


#
# Renombra un item, modificando su path y reconstruyendo el combo
#
proc ::KMat::editTag { T item fullname newtext } {		
    global KPriv		
    set parts [::KEGroups::split2 $fullname //]
    lset parts end m.$newtext
    set newPath [join $parts //]

    #Renombra en la lista de ID's
    set idItem [lindex [$T item text $item] 0]
    set KPriv(materialsId) [::KEGroups::listReplace $KPriv(materialsId) $idItem m.$newtext]
    
    #Cambiar nombre en el árbol
    $T item tag remove $item [list names [$T item tag names $item]]
    $T item tag add $item [EncodeName $newPath]
    $T item element configure $item C0 elemTxtRead -text $newtext
    
    #Cambiar nombre en el XML
    ::KMat::setXml $fullname pid $newtext 
    ::KMat::setXml $fullname id  $newtext
    
    set childs [$T item children $item]
    foreach child $childs {		
	set fullname [DecodeName [$T item tag names $child]]
	set idFull [string map { "." "" "//" ""} $fullname]
	destroy "$T.f$idFull"		
	set f [::KMat::buildFrame $T $child]					
    }
    return $newPath
}

proc ::KMat::GetAutomaticMatName { {auto ""} { startname "" } } {		
    set name ""
    global KPriv
    set i 0
    foreach grup $KPriv(materialsId) {

	incr $i
    }
    if { [llength $KPriv(materialsId)] > 0 } {		
	for {set i 1} {$i<10000} {incr i} {
	    #			set name ${auto}Material${i}
	    set name ${auto}$startname${i}
	    if { [lsearch -exact $KPriv(materialsId) $name] == -1 } { break }
	}
    } else {
	if { $auto == "" } {
	    set name "Material1"
	} else {
	    set name "${auto}Material1"
	}
    }
    return $name
}

#
# Separa cada node en "inicialNombreTag.idNodo"
#
proc ::KMat::splitNode { node } {
    set id [$node getAttribute id ""]		
    if { [$node tagName] == "Container"} {
	return "c.$id"
    } elseif { [$node tagName] == "Item"} {
	return "i.$id"
    } elseif { [$node tagName] == "Property"} {
	return "p.[$node getAttribute id ""]"
    } elseif { [$node tagName] == "Material"} {
	return "m.[$node getAttribute id ""]"
    } else {
	return "NoTree"
    }
}

proc ::KMat::InsertNewItem { propName id T {parent ""} {parentitem root} {childs true} {state "normal"} {open 0}} {
    
    
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
    if { $parent != "" } {
	# set fullname "$parent//$id"
	set fullname "$parent$id"
    } else {
	set fullname $id
    }
    
    if { $childs } {		
	set item [$T item create -button yes -tags [EncodeName $fullname] -open $open]
	$T item lastchild $parentitem $item
	$T item style set $item C0 styAnyRead
	$T item element configure $item C0 elemTxtRead -text "$propName"						
    } else {
	set item [$T item create -button no -tags [EncodeName $fullname] -open $open]
	$T item lastchild $parentitem $item
	$T item style set $item C0 styFrame
	
	set xpath "[::KMat::setXPath $fullname]"
	set value [::xmlutils::getValueText $::KMat::xml $xpath "value"]
	
	$T item style set $item C0 styAnyRead
	if { $parentitem == "root" } {
	    $T item element configure $item C0 elemTxtRead -text "$propName $value"
	} else {
	    $T item element configure $item C0 elemTxtRead -text "$propName: $value"
	}
    }
    
    # Consultamos el icono en el xml, y si existe en nuestro directorio se lo añadimos		
    set icon [::KMat::setXml $fullname icon "" ]
    set imagen [::WinUtils::GetImage $icon]
    if { $imagen != -1 } {
	$T item image $item C0 $imagen
    } else {
    }	
    #Miramos si el item tiene que estar a disabled (viene del proc ::KMProps::stateNode)
    if {$state == "disabled"} {
	$T item enabled $item 0
    }	
    return $item
}


proc ::KMat::buildFrame { T item } {
    
    global KPriv
    set fullname [DecodeName [$T item tag names $item]]
    set idFull [string map { "." "" "//" ""} $fullname]		
    #Comprobamos que sea un item
    if { [::KMProps::itemType $fullname] == "i" } {
	
    } elseif { [::KMProps::itemType $fullname] == "p" } {
	set id [::KMat::setXml $fullname id "" ]
	set pid [::KMat::setXml $fullname pid "" ]
	set unit [::KMat::setXml $fullname unit "" ]
	set icon [::KMat::setXml $fullname icon "" ]
	set state [::KMat::setXml $fullname state "" ]
	set help [::KMat::setXml $fullname help "" ]
	
	set xpath "[::KMat::setXPath $fullname]"
	set comboList [::xmlutils::getValues $::KMat::xml $xpath]
	set value [::xmlutils::getValueText $::KMat::xml $xpath "value"]
	
	#
	#---------------------------#---------------------------#
	# Configurar frame en función del XML
	#---------------------------#---------------------------#
	#
	set bg "#F8F8F8"					
	set f [frame "$T.f$idFull" -borderwidth 0 -background $bg]
	if { [llength $comboList] > 0 } {
	    grid [ttk::combobox $f.cmb -values $comboList -state readonly -textvariable "::KMat::cmb$idFull"] \
		-row 0 -column 0 -padx 3 -sticky ne -in $f
	    
	    ::xmlutils::setComboValue $::KMat::xml $xpath $f.cmb $value
	    #set selected [::xmlutils::getSelected $value $comboList]
	    #$f.cmb current $selected
	    bind $f.cmb <<ComboboxSelected>> "::KMat::cmbSelectChange $item $T 0 "
	} else {
	    grid [ttk::combobox $f.cmb -state normal -textvariable "::KMat::cmb$idFull"] \
		-row 0 -column 0 -padx 3 -sticky nw -in $f
	    set ::KMat::cmb$idFull $value
	    bind $f.cmb <FocusOut> "::KMat::cmbSelectChange $item $T 0 "
	    bind $f.cmb <Escape> "::KMat::cmbCancel $item $T"						
	}
	# Si pulsan intro o Esc también forzamos la salida del combo (por probar)
	bind $f.cmb <KeyPress> "if { %k == 13  } { ::KMat::cmbSelectChange $item $T 1 }"
	#bind $T <KeyPress> "if { %k == 27   } { ::KMProps::cmbCancel $item $T }"
	return $f		  
    } else {				
	return ""
    }
}

proc ::KMat::cmbSelectChange { item T {remove 1} {selectVal current}  } {   
    
    set fullname [DecodeName [$T item tag names $item]]
    set idFull [string map { "." "" "//" ""} $fullname]
    
    set id [::KMat::setXml $fullname id "" ]
    set pid [::KMat::setXml $fullname pid "" ]
    
    set xpath "[::KMat::setXPath $fullname]"
    set value [::xmlutils::getValueText $::KMat::xml $xpath "value"]
    
    #Antes solo se hacía esto:
    #set selCombo [set ::KMat::cmb$idFull]			
    
    #Gestión de ivalues
    set xpath "[::KMat::setXPath $fullname]"
    set comboState [::xmlutils::getComboState $::KMat::xml $xpath]
    if { $comboState == "normal" } {
	
	set selCombo [set ::KMat::cmb$idFull]
	set selComboText $selCombo
	if {$selectVal != "current"} {
	    set ::KMat::lastSelected {}
	}
    } else {
	
	if {$selectVal == "current"} {
	    
	    set f "$T.f${idFull}.cmb"
	    set selCombo [::xmlutils::getComboValue $::KMat::xml $xpath $f]
	    set selComboText [::xmlutils::getComboValue $::KMat::xml $xpath $f "text"]
	    set ::KMat::lastSelected [list $item $selCombo $selComboText]
	    
	} else {
	    
	    set selCombo [lindex $::KFun::lastSelected 1]
	    set selComboText [lindex $::KFun::lastSelected 2]
	    set ::KMat::lastSelected {}
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
    ::KMat::setXml $fullname value $selCombo
    
    ::KMat::insertIcon $item $T $fullname
    
    #Volvemos a cargar el árbol de propiedades
    if { !$remove } {
	#::KMat::refreshTree $T
    }
    ::KMat::refreshTree $T
}

proc ::KMat::cmbCancel { item T  } {
    
    set fullname [DecodeName [$T item tag names $item]]
    set idFull [string map { "." "" "//" ""} $fullname]
    
    set pid [::KMat::setXml $fullname pid ""]
    set id [::KMat::setXml $fullname id ""]
    set value [::KMat::setXml $fullname value ""]
    
    set xpath "[::KMat::setXPath $fullname]"
    set valueText [::xmlutils::getValueText $::KMat::xml $xpath "value"]
    
    set ::KMat::cmb$idFull "$value"
    set f "$T.f$idFull"
    if {[winfo exists $f]} {
	destroy $f
    }
    $T item style set $item C0 styAnyRead
    $T item element configure $item C0 elemTxtRead -text "$pid: $valueText"
    
    ::KMat::insertIcon $item $T $fullname
    set ::KMat::lastSelected ""
    ::KMat::refreshTree
}

proc ::KMat::insertIcon { item T fullname {type "mat"}} {
    
    # Consultamos el icono en el xml, y si existe en nuestro directorio se lo añadimos		
    if {$type == "mat"} {
	set icon [::KMat::setXml $fullname icon ""]
    } else {
	set icon [::KProps::setXml $fullname icon]
    }
    
    set imagen [::WinUtils::GetImage $icon]
    if { $imagen != -1 } {
	$T item image $item C0 $imagen
    }			
}


proc ::KMat::refreshTree { {T ""} } {
    
    #::KMat::findMaterialParent "S1"
    
    if {$T == ""} {
	set T $::KMat::TreeMatPath
    }
    
    #Primero hay que asegurarse de que exista la ventana
    if { [winfo exists $::KMProps::WinPath] } {
	
	foreach item [$T item range 0 end] {	  
	    
	    set fullname [DecodeName [$T item tag names $item]]						
	    if {$fullname != "" } {		
		catch {
		    ::KMat::setXml $fullname open [$T item isopen $item]
		}
	    }
	}
    }
    ::KMat::DeleteTree		
    ::KMat::FillTreeMat
    
    set ::KMat::lastSelected {}
}

proc ::KMat::findMaterialParent { matid } {

    global KPriv

    set nodes [$::KMat::xml selectNodes "/Kratos_KMat_DB/Materials/MaterialGroup\[@id\]"]
    
    foreach node $nodes {			  
	set nodes2 [$node childNodes]
	foreach node2 $nodes2 {  
	    set aux [$node2 getAttribute id ""]
	    if { $aux == $matid} {
		set parent [$node getAttribute id ""]
		return $parent
	    }
	}		
    }
}

proc ::KMat::getMaterials {{application ""}} {

    global KPriv
    set KPriv(materialsList) {}
    set T $::KMat::TreeMatPath
    
    set nodes [$::KMat::xml selectNodes "/Kratos_KMat_DB/Materials/MaterialGroup\[@id\]"]
    
    if {$application == ""} {
	foreach node $nodes {		
	    set nodes2 [$node childNodes]
	    foreach node2 $nodes2 {  
		set aux [$node2 getAttribute id ""]
		lappend KPriv(materialsList) $aux				
	    }		
	}
    } else {
	if {$application == "StructuralAnalysis"} {
	    set smatlist [list "Metal" "Composite" "Plastic"]
	    foreach node $nodes {		
		if {[$node getAttribute id ""] in $smatlist} {
		    set nodes2 [$node childNodes]
		    foreach node2 $nodes2 {  
			set aux [$node2 getAttribute id ""]
			lappend KPriv(materialsList) $aux				
		    }		
		}
	    }
	} elseif {$application == "Fluid"} {
	    foreach node $nodes {		
		if { [$node getAttribute id ""] == "Fluid"} {
		    set nodes2 [$node childNodes]
		    foreach node2 $nodes2 {  
			set aux [$node2 getAttribute id ""]
			lappend KPriv(materialsList) $aux				
		    }		
		}
	    }
	} elseif {$application == "FluidStructureInteraction"} {
	    foreach node $nodes {		
		if { [$node getAttribute id ""] == "Metal" || [$node getAttribute id ""] == "Fluid"} {
		    set nodes2 [$node childNodes]
		    foreach node2 $nodes2 {  
			set aux [$node2 getAttribute id ""]
			lappend KPriv(materialsList) $aux				
		    }		
		}
	    }
	} elseif {$application == "ConvectionDiffusion"} {
	    foreach node $nodes {		
		if { [$node getAttribute id ""] == "Plastic"} {
		    set nodes2 [$node childNodes]
		    foreach node2 $nodes2 {  
			set aux [$node2 getAttribute id ""]
			lappend KPriv(materialsList) $aux				
		    }		
		}
	    }
	} else {
	    foreach node $nodes {		
		set nodes2 [$node childNodes]
		foreach node2 $nodes2 {  
		    set aux [$node2 getAttribute id ""]
		    lappend KPriv(materialsList) $aux				
		}		
	    }					
	}
    }
    
    return $KPriv(materialsList)
}

proc ::KMat::getMaterialProperty { matid prop } {

    global KPriv
    set KPriv(materialsList) {}
    
    set xpath "/Kratos_KMat_DB/Materials/MaterialGroup/Material\[@id=\"$matid\"]/Property\[@id=\"$prop\"]"
    set node [$::KMat::xml selectNodes $xpath]

    set value [$node @value]
    return $value
}

proc ::KMat::getMaterialProperties { pathtype path } {
    
    # pathtype = c if you pass a container path
    # pathtype = p if you pass a property path
    global KPriv
    if { $pathtype == "c" } {
	set xpath "[::KMat::setXPath $path]/Property\[@id\]"
    } elseif { $pathtype == "p" } {
	set xpath "[::KMat::setXPath $path]"
    }
    
    set matPropsList {}
    
    set nodes [$::KMat::xml selectNodes "$xpath"]
    
    foreach node $nodes {			  
	set aux "[$node getAttribute id ""] [$node @value]"
	lappend matPropsList $aux				
    }

    return $matPropsList

}


proc ::KMat::testbutton { } {
    
    set matlist [::KMat::getMaterialProperties]
}


#
# Editar o extraer propiedades del xml en memoria
#
proc ::KMat::setXml { path property {value ""} } {
    
    global KPriv
    
    set xpath "[::KMat::setXPath $path]"
    
    if { $value == "" } {
	set value [$::KMat::xml set "$xpath/@$property" ]				
	
	
	#Cuando hay espacios el xml devuelve una lista y si la imprimes tal cual aparecen corchetes
	if { [llength $value] == 1 } {
	    set value [lindex $value 0]
	}
	return $value
    } else {
	$::KMat::xml set "$xpath/@$property" "$value"
	return "1"
    }
}


#
# Construye un frame con un tab por cada container que cuelgue de "item"
# y dentro de cada tab, etiquetas y combos para cada item (si los hay)
#
proc ::KMat::buildTabFrame { T item {class "Tab"} } {

    variable abdlist

    set fullname [DecodeName [$T item tag names $item]]
    set id [::KMat::setXml $fullname id "" ]
    set pid [::KMat::setXml $fullname pid "" ]

    set f [::KMat::iniFrameBottom]
    
    # Miramos los descendientes directos y si son container ponemos un tab por cada uno q tenga items
    set children [$T item children $item]

    set listTabs {}
    set listItems {}
    set acceptItems {}

    set nb ${f}.nb
    #	grid [ttk::notebook $nb ] -row 0 -sticky ewn
    grid [ttk::notebook $nb ] -row 0 -column 0 -columnspan 2 -padx 0 -sticky nw -in $f

    # declaramos un tab para el material
    set fTab ${nb}.f$id
    $nb add [ttk::labelframe $fTab -text "[= Properties]" -padding {10 0 10 10}] \
	-text "[string range $pid 0 20]"

    # grid [ttk::label $fTab.lblName -text "[= Property Name:]" ] \
	#		 -row 0 -column 0 -pady 5 -sticky nw -in $fTab
    # grid [ttk::combobox $fTab.cmbPropertyName -state normal -textvariable "::KMProps::propertyName" ] \
	#		 -row 0 -column 1 -padx 3 -pady 5 -sticky nw -in $fTab


    set count 0
    foreach itemChild $children {
	set fullname [DecodeName [$T item tag names $itemChild]]
	set nodeName [::xmlutils::getXmlNodeName $fullname "mat"]
	
	set xpath "[::KMat::setXPath $fullname]"
	set comboList [::xmlutils::getValues $::KMat::xml $xpath]
	
	lappend listItems $itemChild

	set id [::KMat::setXml $fullname id "" ]
	set pid [::KMat::setXml $fullname pid "" ]
	set value [::KMat::setXml $fullname value "" ]
	set help [::KMat::setXml $fullname help "" ]

	if { $id ni $abdlist && $id != "ABD" } {
	    grid [ttk::label $fTab.label$id -text "$pid" ] \
		-row $count -column 0 -pady 5 -sticky nw -in $fTab
	    grid [ttk::combobox $fTab.cmb$id -state normal -textvariable "::KMat::cmb$id" ] \
		-row $count -column 1 -padx 3 -pady 5 -sticky nw -in $fTab
	    set ::KMat::cmb$id $value
	    incr count
	}
	if { $id == "ABD" } {
	    set clase "Tab"
	    grid [ttk::button $fTab.bProp$id -text "ABD Matrix"  -command "::KMat::buildTabFrameABD $T $item $clase" ] \
		-row $count -column 0 -sticky sw  -pady 3 -padx 20  -in $f
	    tooltip::tooltip $fTab.bPropABD [= "ABD Matrix"]

	}
    }

    set clase "Tab"
    grid [ttk::button $f.bPropOk -text "Ok"  -command "::KMat::acceptTabFrame $T $listItems $class $item" ] \
	-row $count -column 0 -sticky sw  -pady 3 -padx 20  -in $f
    tooltip::tooltip $f.bPropOk [= "Confirm values"]
    
    grid [ttk::button $f.bPropCancel -text "Cancel"  -command "::KMat::cancelBottom" ] \
	-row $count -column 0 -sticky sw  -pady 3 -padx 100  -in $f
    tooltip::tooltip $f.bPropCancel [= "Cancel assignation"]

}


proc ::KMat::buildTabFrameABD { T item {class "Tab"} } {

    variable abdlist

    set fullname [DecodeName [$T item tag names $item]]
    set id [::KMat::setXml $fullname id "" ]
    set pid [::KMat::setXml $fullname pid "" ]

    set f [::KMat::iniFrameBottom]
    
    # Miramos los descendientes directos y si son container ponemos un tab por cada uno q tenga items
    set children [$T item children $item]

    set listTabs {}
    set listItems {}
    set acceptItems {}

    set nb ${f}.nb
    #	grid [ttk::notebook $nb ] -row 0 -sticky ewn
    grid [ttk::notebook $nb ] -row 0 -column 0 -columnspan 2 -padx 0 -sticky nw -in $f

    # declaramos un tab para el material
    set fTab ${nb}.f$id
    $nb add [ttk::labelframe $fTab -text "[= Properties]" -padding {10 0 10 10}] \
	-text "[string range $pid 0 20]"

    set count 1
    # foreach itemChild $children {
    #		 set fullname [DecodeName [$T item tag names $itemChild]]
    #		 set nodeName [::xmlutils::getXmlNodeName $fullname "mat"]
    #		 set comboList [::xmlutils::getXMLValues $fullname]

    #		 lappend listItems $itemChild

    #		 set id [::KMat::setXml $fullname id ""]
    #		 set pid [::KMat::setXml $fullname pid "" ]
    #		 set value [::KMat::setXml $fullname value "" ]
    #		 set help [::KMat::setXml $fullname help "" ]

    #		 if { $id ni $abdlist } {
    #			 grid [ttk::label $fTab.label$id -text "$pid" ] \
	#				 -row $count -column 0 -pady 5 -sticky nw -in $fTab
    #			 grid [ttk::combobox $fTab.cmb$id -state normal -textvariable "::KMat::cmb$id" ] \
	#				 -row $count -column 1 -padx 3 -pady 5 -sticky nw -in $fTab
    #			 set ::KMat::cmb$id $value
    #			 incr count
    #		 }
    # }

    grid [ttk::button $f.bPropOk -text "Ok"  -command "::KMat::acceptTabFrameABD $T $listItems $class $item" ] \
	-row $count -column 0 -sticky sw  -pady 3 -padx 20  -in $f
    tooltip::tooltip $f.bPropOk [= "Confirm values"]
    
    grid [ttk::button $f.bPropCancel -text "Cancel"  -command "::KMat::cancelBottom" ] \
	-row $count -column 0 -sticky sw  -pady 3 -padx 100  -in $f
    tooltip::tooltip $f.bPropCancel [= "Cancel assignation"]







    # # Create the matrix frame
    # set tfMatrix [TitleFrame $fpath.tfMatrix${matrixid} \
	#					   -relief groove \
	#					   -bd 2 \
	#					   -side left \
	#					   -text $matrixtxt]
    
    # set tffMatrix [$tfMatrix getframe]

    # Label $tffMatrix.lMatrixId${matrixid} \
	#		 -text $matrixid 

    # # Put the frame
    # grid $fpath.tfMatrix${matrixid} \
	#		 -in $fpath \
	#		 -row $row -column $col \
	#		 -sticky nswe

    # # Put all entry/grid for this matrix
    # set ewidth 12
    # # Row 1
    # set crow 1
    # for {set ii 0} {$ii<=2} {incr ii 1} {
    #		 set ccol [expr $ii+1]
    #		 set comp ${matrixid}${crow}${ccol}
    #		 Entry $tffMatrix.e${comp} \
	#			 -width $ewidth \
	#			 -editable 0 \
	#			 -helptext [= "Current matrix value (%s)" $comp].

    #		 grid $tffMatrix.e${comp} \
	#			 -in $tffMatrix \
	#			 -row 0 -column $ccol \
	#			 -sticky ne
    # }
    
    # grid $tffMatrix.lMatrixId${matrixid} \
	#		 -in $tffMatrix \
	#		 -row 1 -column 0 \
	#		 -sticky ne

    # # Row 2
    # set crow 2
    # for {set ii 0} {$ii<=2} {incr ii 1} {
    #		 set ccol [expr $ii+1]
    #		 set comp ${matrixid}${crow}${ccol}
    #		 Entry $tffMatrix.e${comp} \
	#			 -width $ewidth \
	#			 -editable 0 \
	#			 -helptext [= "Current matrix value (%s)" $comp].

    #		 grid $tffMatrix.e${comp} \
	#			 -in $tffMatrix \
	#			 -row 1 -column $ccol \
	#			 -sticky ne
    # }

    # # Row 3
    # set crow 3
    # for {set ii 0} {$ii<=2} {incr ii 1} {
    #		 set ccol [expr $ii+1]
    #		 set comp ${matrixid}${crow}${ccol}
    #		 Entry $tffMatrix.e${comp} \
	#			 -width $ewidth \
	#			 -editable 0 \
	#			 -helptext [= "Current matrix value (%s)" $comp].

    #		 grid $tffMatrix.e${comp} \
	#			 -in $tffMatrix \
	#			 -row 2 -column $ccol \
	#			 -sticky ne
    # }

    # # Expand
    # grid rowconfigure $tffMatrix {0 1 2} -weight 1
    # grid columnconfigure $tffMatrix {0 1 2 3} -weight 1













}


#
# Crea el frame inferior 
#
proc ::KMat::iniFrameBottom { } {
    
    #Destruye el frame inferior si existía
    set f [::KMat::cancelBottom]
    
    # Create the frame where set the properties
    ttk::frame $f -borderwidth 0
    
    # Grid for toolbar
    grid $f -row 2 -column 0 -sticky wes
    
    return $f
}

#
# Destruye el frame inferior 
#
proc ::KMat::cancelBottom { } {
    
    set f ${::KMat::WinPath}.nb.fMat.fBottom
    if { [winfo exists $f]} {
	foreach w [winfo children $f] {
	    destroy $w
	}
    }
    destroy $f
    
    return $f
}




proc ::KMat::acceptTabFrameABD { T listItems class {itemSel ""}} {
    WarnWin [_ "Accept Tab Frame ABD, save values in XML"]
}


proc ::KMat::acceptTabFrame { T listItems class {itemSel ""}} { 
    WarnWin [_ "Accept Tab Frame, save values in XML"]
}

proc ::KMat::initVisibilityClass { } {
    
    variable visibilityVars
    set visibilityVars {}
    global KPriv

    set classes [$::KMat::xml set "/ Kratos_KMat_DB/Materials/ClassConfiguration/Class" ]
    set classes [split $classes ","]
    
    foreach class $classes {
	lappend visibilityVars $class
	set ::KMat::$class ""
    }
}


#
# Valida varias cosas para cada nodo
#
proc ::KMat::stateNode { node } {
    
    #Validamos para cada nodo si tiene que estar visible 
    #(en función de los valores elegidos en algunos combos)
    if { [$node nodeName] == "Property" } {

	#Leemos la class del nodo para ver si requiere de acciones especiales (ocultar nodos)
	set value [$node getAttribute value ""]
	set class [$node getAttribute class ""]

	#Equivalente a Switch $class
	foreach var $::KMat::visibilityVars {	
	    if {$var == $class} {				
		#Caso especial para el solver de fluidos
		if { $var == "fluidSolvTyp" } {   
		    #El solver de fluidos está duplicado dependiendo de una variable prebia
		    set freeYesOrNo [$node getAttribute freeSurf ""]
		    #msg "FREEEEEEEEEEEEEE:  $freeYesOrNo == \"\" || $freeYesOrNo == $::KMProps::freeSurf"
		    if { $freeYesOrNo == "" || $freeYesOrNo == $::KMProps::freeSurf } {
			set ::KMProps::fluidSolvTyp "$value"
		    }
		} else {
		    #Caso general
		    set ::KMat::$var $value
		}
	    }
	}
    } else {
	
	#Caso especial para application
	set class [$node getAttribute class ""]
	if { $class == "application" } {			
	    set apliState [$node getAttribute state ""]
	    #msg "apliState$apliState"
	    if {$apliState != "hiddenAll" } {
		set ::KMProps::application [$node getAttribute id ""]
	    }
	}
    }
    
    set state [$node getAttribute state "normal"]
    #msg "\n[$node getAttribute id ""]				state:$state"
    
    #Si el estado es hiddenAll se oculta el nodo y toda su descendencia
    if {$state == "hiddenAll"} {
	return "-1"
    }
    
    foreach var $::KMat::visibilityVars {
	set globalVar [set ::KMat::$var]		
	set nodeValuesVar [split [$node getAttribute $var ""] ","]		
	#Si el nodo tiene alguna restriccion de clase (p.ej. del tipo strucType=Shell)
	# y no coincide con el valor seleccionado, ocultamos el nodo
	if { $nodeValuesVar != "" && !($globalVar in $nodeValuesVar) } {			
	    if {$var == "strucType" } {
		if { $globalVar != "Generic"} { 
		    #msg "$nodeValuesVar \"\" !=  && $globalVar in $nodeValuesVar"
		    return -1
		}
	    } else {
		#msg "$nodeValuesVar \"\" !=  && $globalVar in $nodeValuesVar"
		return -1
	    }
	}
    }
    #Devolvemos el estado del nodo ("normal" por defecto)
    return $state	
}

#
#
#
# Prepara la query para utilizar las funciones de domNOde
proc ::KMat::setXPath { path } {
    
    set splitted [::KMProps::split2 $path //]
    
    set i 0
    foreach itemId $splitted {
	if { $i == 0 } {				
	    #set xpath "/Kratos_KMat_DB/Materials/Material\[@id='$itemId'\]"		  
	    set xpath "/Kratos_KMat_DB/Materials/MaterialGroup\[@id='$itemId'\]"		 
	    
	} else {
	    if { [string index $itemId 0] == "p" } {
		set xpath "$xpath/Property\[@id='[string range $itemId 2 end]'\]"
	    }
	    if { [string index $itemId 0] == "m" } {
		set xpath "$xpath/Material\[@id='[string range $itemId 2 end]'\]"
	    }
	    if { [string index $itemId 0] == "c" } {
		set xpath "$xpath/Container\[@id='[string range $itemId 2 end]'\]"
	    }
	}		
	incr i
    }
    return $xpath
}
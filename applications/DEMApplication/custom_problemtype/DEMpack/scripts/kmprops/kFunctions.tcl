##############################################################################################
#
#        NAME: listEntities.tcl
#
#        PURPOSE: List Entities for each group 
#
#        QUANTECH ATZ-DEVELOPMENT DEPARTMENT
#
#        AUTHOR : Luis
#
#        CREATED AT: 20/05/10
#
#        LAST MODIFICATION : 
#
#        VERSION : 0.1
#
#        HISTORY:
#
#         
#
##############################################################################################

package require treectrl
package require tooltip
package provide KEGroups 1.0 
package require snit
package require tdom
package require xmlstruct 1.0

# Create a snit type to store the entity group properties
#snit::type EGroup {}

# Create a base namespace KEGroups
namespace eval ::KFun:: {
    
    # Path of the base window 
    variable WinPath ".gid.kmprops.kfunctions"
    variable TreePath   ;# The tree path
    variable lSort 0
    variable lastSelected {}
    
    variable applyFun 0
    variable propsfTab ""
    variable propsId ""
}

proc ::KFun::InitBaseWindow { {fTab ""} {id ""} } {
    
    if { $fTab != "" } {
	    
	    set ::KFun::applyFun 1
	    set ::KFun::propsfTab $fTab
	    set ::KFun::propsId $id
    }
    
    global KPriv
    
    set w "$::KFun::WinPath"

    # Open the window outside
    ::KFun::OpenWindowOutside $w

    if {![winfo exists $w]} return 
	
	KFun::initXML
	
    set T [::KFun::CreateTreeAndToolbar $w]
    
    ::KFun::fillFunctions $T
    
    #::KFun::buildEditor $T
    
    
    # Binding
    bind $w <Alt-c> "destroy $w"
    #bind $w <Escape> "destroy $w"
    
    #Si está abierta la ventana de propiedades la bloqueamos
    if {[info exists ::KMProps::WinPath]} {
		if {[winfo exists $::KMProps::WinPath]} {
		        #grab $::KFun::WinPath        
		}
	}
    
}

proc ::KFun::initXML { } {
	
	global KPriv

	if { $KPriv(xmlFun) == "" } {
		
		set filePath "$KPriv(dir)/python_functions.xml"
		
		::kfiles::MakeBackupCopyOfSPDFile $filePath ".xml"
		
		set xmlArray [::xmlutils::openFile "." "$filePath" 0]
		   
		   set KPriv(xmlFun) [lindex $xmlArray 0]
		   set KPriv(xmlDocFun) [lindex $xmlArray 2]
		   set KPriv(encrXmlFun) [lindex $xmlArray 1]
	}
}

proc ::KFun::CreateTreeAndToolbar { w } {
    
    variable TreePath   
	
	set T [::KFun::CreateTreeProperties $w]
	
	
	#TREE-FRAME
    # Create the treectrl properties
    set mdf [ttk::frame $w.middle ]
    
    # Set the tree path to a namespace variable
    set TreePath $T
    
    grid $mdf -sticky new
    
	#TOOL-BAR
	set tbf [ttk::frame $w.tbf -borderwidth 3]

	grid [ ttk::button $tbf.newFunction -image [::WinUtils::GetImage new_tree.gif] -command [list ::KFun::CreateNewFunction] -style Toolbutton ] \
	-row 0 -column 0 -sticky wes
	tooltip::tooltip $tbf.newFunction [= "Add a new function"]
	
	grid [ ttk::button $tbf.deleteFunction -image [::WinUtils::GetImage delete_tree.gif] -command [list ::KFun::DeleteFunction] -style Toolbutton ] \
	-row 0 -column 1 -sticky wes
	tooltip::tooltip $tbf.deleteFunction [= "Delete the selected function"]
	
	grid $tbf -sticky sew
	
	#::KEGroups::FillTree
	
	if { $::KFun::applyFun } {
	    grid [ ttk::button $tbf.applyFun -text [= "Apply"] -command [list ::KFun::applyFunction $T]] \
		-row 2 -column 0 -columnspan 2 -sticky wes
		tooltip::tooltip $tbf.applyFun [= "Assign the selection function and close editor"]
	}
		
	#CLOSE BUTTON
	## For lower buttons
	#set tf [ttk::frame $w.buttons] 
	#grid $tf -sticky ews
	#grid anchor $tf center
	#
	#grid [ttk::button $tf.bClose -text [= "Close"] -command [list destroy $w]]  -sticky ew -padx 5 -pady 6
		
		
	#Abre una ventana con el editor de GiD
	#grid [ ttk::button $tbf.openEditor -image [::WinUtils::GetImage list_entities.gif] -command ::KFun::OpenNotes -style Toolbutton ] \
	#-row 0 -column 2 -sticky wes
	#tooltip::tooltip $tbf.openEditor [= "Open the GiD editor"]
	
	
    
    focus $T
    
    return $T
}

proc ::KFun::applyFunction { T } {
	
	set item [$T selection get]
	
	if { $item == ""} {
		WarnWin [= "Please, select the function you want to assign to '%s'." $::KFun::propsId]
		return ""
	} else {
		set idFunction [$T item tag names $item]
		::KMProps::assignFunction $::KFun::propsfTab $::KFun::propsId $idFunction
	}
	destroy $::KFun::WinPath
}

proc ::KFun::ClickTree { x y T } {
    
    set info [$T identify $x $y]
    
    if { [lindex $info 0] == "item" && [llength $info] >= 4 } {
	
	set item [lindex $info 1]
	set col [lindex $info 3]
	
	if { $col == "1" } {
		$T selection clear
		$T selection add $item
	}
	
	#::KFun::buildEditor $T $item
	
    } elseif { [lindex $info 0] == "header" && [lindex $info 1] == "0" } {
	
	if { [$T column cget C0 -arrow] == "up" } {
	    $T column configure C0 -arrow down
	    $T item sort 0 -dictionary -increasing
	} else {
	    $T column configure C0 -arrow up
	    $T item sort 0 -dictionary -decreasing
	}
	#msg "Info:$info"
	return ""
    } else {
	return ""
    }
    
    if { $col != 0 } {
	return -code break
    }
    return ""
    
}

proc ::KFun::DoubleClickTree { x y T } {
    
    
    set info [$T identify $x $y]
    
    if { [lindex $info 0] == "item" && [llength $info] >= 4 } {
	
		set item [lindex $info 1]
		set col [lindex $info 3]
		
		if { $col == "1" || $col == "2" } {
		        $T selection clear
		        $T selection add $item
		}
		
		
		#Eliminamos el anterior combo, si aun está visible
		if {[llength $::KFun::lastSelected] > 0 && [lindex $::KFun::lastSelected 0] != $item} { 
		        global KPriv                
		        #set state [::xmlutils::getComboBoxState $fullname]
		        ::KFun::cmbSelectChange $T $KPriv(xmlDocFun) [lindex $::KFun::lastSelected 0] 1 "anterior"
		}
		
		if { [$T item children $item] != ""} {
		        ::KFun::buildEditor $T $item
		} else {
		        set f [::KFun::buildFrame $T $item]  
		        $T item style set $item C0 styFrame
		        $T item element configure $item C0 eWindow -window $f
		}
	}
}

proc ::KFun::CreateTreeProperties {w} {

    # Scrollbars
    set vsb $w.vsb1
    set hsb $w.hsb1

    # Create the treectrl and set the scrollbars
    set T [treectrl $w.t -xscrollcommand [list $hsb set] -yscrollcommand [list $vsb set]]
    ttk::scrollbar $vsb -orient vertical   -command [list $T yview]
    ttk::scrollbar $hsb -orient horizontal -command [list $T xview]
    
    # Set the height
    set height [font metrics [$T cget -font] -linespace]
    if {$height < 22} {
	set height 22
    }

    # Configure the treectrl
    $T configure -indent 30 -itemheight $height -selectmode browse \
	-showroot 0 -showrootbutton 0 -showbuttons 1 -showlines 1 \
	-highlightthickness 0 -borderwidth 0 \
	-xscrollincrement 20 -yscrollincrement 50
    
    # Create the column identifier list
    set collistid [list [= "Function Id"] [= "Function definition"]]
    
    set i 0
    foreach cid $collistid {
	$T column create -text $cid -tags C$i -weight 0
	incr i
    }
    # Configure the column weight and arrow
    $T column configure C0 -weight 1 -arrow up

    # Configure the column that have the tree
    $T configure -treecolumn C0

    # Create elements
    $T element create elemImgAny image
    $T element create elemTxtRead text -fill [list $::KMProps::SystemHighlightText {selected focus}] -lines 1
    $T element create elemRectSel rect -fill [list $::KMProps::SystemHighlight {selected focus} gray {selected !focus}] -showfocus yes
    $T element create elemRectColor rect -width 30
    $T element create eWindow window
    
    
    # Create styles using the elements
    set S [$T style create styAnyRead]
    $T style elements $S {elemImgAny elemRectSel elemTxtRead }
    $T style layout $S elemImgAny -expand ns
    $T style layout $S elemTxtRead -padx 4 -expand ns -squeeze x
    $T style layout $S elemRectSel -union [list elemTxtRead] -iexpand ns -ipadx 2
    
    set S [$T style create styRectColor]
    $T style elements $S {elemRectColor elemTxtRead}
    $T style layout $S elemRectColor -union [list elemTxtRead] -iexpand ns -ipadx 2 -pady 1

    set S [$T style create styAnyImage]
    $T style elements $S {elemImgAny}
    $T style layout $S elemImgAny -expand ns -pady 1
    
    set S [$T style create styFrame -orient horizontal]
    $T style elements $S {eWindow}
    $T style layout $S eWindow -iexpand x -squeeze x -padx {0 1} -pady {2 2}
    
    # Items
    set item root
    $T item configure $item -button yes
    $T item style set $item C0 styAnyRead
    $T item element configure $item C0 elemTxtRead -text [= "Root"]

    # List of lists: {column style element ...} specifying text elements
    # the user can edit
    TreeCtrl::SetEditable $T {
	}
	#{C0 styAnyRead elemTxtRead}
    
    
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
    #$T notify install <Edit-accept>

    # Notify bind
    # TODO
    #$T notify bind DragTag <Drag-receive> { ::KFun::ReceiveDragGroups  %T %l %I }
    #$T notify bind EditTag <Edit-accept> { ::KFun::SetGroupsToRename %T %I %t }
	
    bind $T <Button-1> [list ::KFun::ClickTree %x %y $T]
    bind $T <Double-Button-1> [list ::KFun::DoubleClickTree %x %y $T]
    #bind $T <Return> [list SetLayersTo TOUSE $T]
    #bind $T <Key-Delete> [list SetLayersToDelete $T]
    #bind $T <Alt_L> [list InvertSelectionTableList $T]
    #bind $T <Alt_R> [list InvertSelectionTableList $T]
    #bind $T <Meta_L> [list InvertSelectionTableList $T]
    #bind $T <Meta_R> [list InvertSelectionTableList $T]
    #bind $T <F2> [list ::KFun::BeginEditGroups $T]

    #bind $T <Button-3> "[list ::KFun::MenuContextualGroup %W %x %y] ; break"
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

proc ::KFun::ReceiveDragGroups { T dragged_list dst } {
    
}

proc ::KFun::OpenWindowOutside { w } {
    
    if { [winfo exists $w] } {
	destroy $w
    }

    # Init the window
    set title [= "Function manager"]
    InitWindow $w $title wKFunctionsGeom ::KFun::InitBaseWindow
    
    return $w
}


proc ::KFun::InsertNewItem { T id pid {definition ""} {parentitem root} {childs 0}} {
    
    if { $childs } {
	
		set item [$T item create -button yes -tags $id]
		$T item lastchild $parentitem $item
		$T item style set $item C0 styAnyRead
		$T item element configure $item C0 elemTxtRead -text $pid
	
    } else {
	
		set item [$T item create -button no -tags $id]
		$T item lastchild $parentitem $item
		$T item style set $item C0 styAnyRead
		$T item element configure $item C0 elemTxtRead -text $pid
    }
    
    #De momento no hay iconos
    set icon ""
    if {$icon != ""} {
		set imagen [::WinUtils::GetImage $icon]
		if { $imagen != -1 } {
		    $T item image $item C0 $imagen
		}
    }
    
    $T item style set $item C1 styAnyRead
    $T item element configure $item C1 elemTxtRead -text $definition
    
    $T collapse $item
    return $item
}

#
# 
#
proc ::KFun::fillFunctions { T } {
    
    global KPriv
    set xml $KPriv(xmlDocFun)
    set xpath "Kratos_Functions/Functions/Function"
    set nodes [$xml selectNodes "$xpath"]
    
    foreach node $nodes {
		
		set id [$node getAttribute id ""]
		
		set itemXPath "${xpath}\[@id='$id'\]/Item\[@id='Definition'\]"
		set definition [::xmlutils::getAttribute $xml $itemXPath sourceCode]
		
		#set definition [::KFun::nodeDefinition $id]
		
		set parentItem [::KFun::InsertNewItem $T $id [$node getAttribute pid ""] $definition root 1]
		
		foreach childNode [$node childNodes] {
		        
		        set id [$childNode getAttribute id ""]
		        set pid [$childNode getAttribute pid ""]
		        
		        set dvText [::xmlutils::getValueText $xml [$childNode toXPath]]
		        
		        if { $id != "Definition" } {
		                
		                ::KFun::InsertNewItem $T $id "$pid: $dvText" "" $parentItem
		        }
		}
	}
}

proc ::KFun::getFunctionAtt { path attr } {
    
    global KPriv
    
    set nodes [$KPriv(xmlDocFun) selectNodes "Kratos_Functions/Functions/Function/Item"]
    
    foreach node $nodes {
		
		::KFun::InsertNewItem $T [$node getAttribute id ""] [$node getAttribute pid ""] [$node text]
	}
}

#
# Codifica el pid elegido por al usuario para que sea un id
#
proc ::KFun::encodeId {text} {
	set texto [string map { " " "_" } $text]
	return $texto
}

proc ::KFun::decodeId {text} {
	set texto [string map { "_" " " } $text]
	return $texto
}

#
# Crea el frame inferior 
#
proc ::KFun::iniFrameBottom { } {
    
    #Destruye el frame inferior si existía
    set f [::KFun::cancelBottom]
    
    # Create the frame where set the properties
    ttk::frame $f -borderwidth 0
    
    # Grid for toolbar
    grid $f -row 3 -column 0 -sticky wes
    
    return $f
}

#
# Destruye el frame inferior 
#
proc ::KFun::cancelBottom { } {
    
    #Si estamos en modo "aplicar función" mostramos el boton Apply
    if { $::KFun::applyFun } {
	    
		grid $::KFun::WinPath.tbf.applyFun
	}
	
    set f ${::KFun::WinPath}.fBottom
    if { [winfo exists $f]} {
	foreach w [winfo children $f] {
	    destroy $w
	}
    }
    destroy $f
    
    return $f
}

proc ::KFun::getComboVars { } {
	
	global KPriv
	
	#COMBO VARIABLES
	set comboVars {}
	set node [$KPriv(xmlDocFun) selectNodes "Kratos_Functions/KratosVariables"]
	if {$node != ""} {
		set comboVars [split "[$node text]" ","]
	}
	
	return $comboVars
}

proc ::KFun::insertText { f varText } {
	
	set text [set $varText]
	$f insert [$f index insert] $text 
}

proc ::KFun::getXMLFunctions { path } {
	
	global KPriv
	
	#COMBO Functions
	set lista {""}
	set nodes [$KPriv(xmlDocFun) selectNodes "Kratos_Functions/${path}/Item"]
	foreach node $nodes {
		
		lappend lista [$node getAttribute pid ""]
	}
	
	return $lista
}

proc ::KFun::nodeDefinition { id {action "read"} {value ""} } {

	global KPriv
	set xml $KPriv(xmlDocFun)
	
	#Function Text
	set xpath "Kratos_Functions/Functions/Function\[@id='$id'\]/Item\[@id='Definition'\]"
	set node [$xml selectNodes $xpath]
	if {$node != ""} {
		if { $action == "read" } {
				return [::xmlutils::getAttribute $xml $xpath sourceCode]
		        #return [::xmlutils::AsXml [$node text] "toTxt"]
		} else {
				#msg "xpath\n[$node toXPath]"
				::xmlutils::getAttribute $xml xpath sourceCode "$value"
			
				return ""
		        if { $value == "" } {
		        	set value " "
		        }
		        	
		        set childNodes [$node childNodes]
		        if { $childNodes != "" } {
		                
		                #$childNodes appendXML "<SourceCode>$value</SourceCode>"
		                $childNodes nodeValue $value
		        } else {
		                set nodetext [$node asXML]
		                #msg $nodetext
		                #$node nodeValue $value
		        }
		        return ""
		}
	}
}

proc ::KFun::buildEditor { T {item ""} } {
	
	global KPriv
	set xml $KPriv(xmlDocFun)
	
	set f [::KFun::iniFrameBottom]
	
	set width 25
	
	if {$item != ""} {
		
		set pid [$T item text $item 0]
		set id [$T item tag names $item]
		
		set xpath "Kratos_Functions/Functions/Function\[@id='$id'\]"
		
		set dv [::xmlutils::getAttribute $xml $xpath dv]        
		set tooltip [::xmlutils::getAttribute $xml $xpath help]
		
	} else {
		set pid ""
		set id ""
		set dv ""
		set xpath "Kratos_Functions/Templates/Function\[@id='Function'\]"
	}
		                            
	#NOM ID
	grid [ttk::label $f.lblNomId -text [= "Name id:"] ] \
	-row 1 -column 0 -pady 2 -padx 3 -sticky nw -in $f
	
	grid [ttk::entry $f.txtNomId -cursor ibeam -textvariable ::KFun::nomId -width [expr $width + 3]  ] \
    -row 1 -column 1 -pady 2 -padx 3 -sticky nw -in $f
	
	set ::KFun::nomId "$pid"
	
	
	## PYTHON-FUN
	#grid [ttk::label $f.lblPythonFun -text [= "Kratos Filter"] ] \
	#-row 2 -column 0 -padx 3 -pady 2 -sticky nw -in $f
	#
	#grid [ttk::combobox $f.cmbPythonFun -values [::KFun::getXMLFunctions "KratosPythonFun"] -state readonly -width $width -textvariable "::KFun::cmbPythonFun"] \
	#-row 2 -column 1 -padx 3 -pady 2 -sticky nw -in $f
	#tooltip::tooltip $f.cmbPythonFun [= "Kratos Filter"]
	
	
	#
	# AppliedTo COMBO
	#

	set cmbXPath "${xpath}/Item\[@id='AppliedTo'\]"
	set comboList [::xmlutils::getValues $xml $cmbXPath]
	set cmbPid [::xmlutils::getAttribute $xml $cmbXPath pid]
	set cmbDv [::xmlutils::getAttribute $xml $cmbXPath dv]
	set cmbTooltip [::xmlutils::getAttribute $xml $cmbXPath tooltip]
	
	grid [ttk::label $f.lblAppliedTo -text [= $cmbPid]:] \
	-row 2 -column 0 -padx 3 -pady 2 -sticky nw -in $f
	
	grid [ttk::combobox $f.cmbAppliedTo -values $comboList -state readonly -width $width -textvariable "::KFun::cmbAppliedTo"] \
	-row 2 -column 1 -padx 3 -pady 2 -sticky nw -in $f
	tooltip::tooltip $f.cmbAppliedTo [= $cmbTooltip]
	
	::xmlutils::setComboValue $xml $cmbXPath $f.cmbAppliedTo $cmbDv
	
	
	#
	# FunctionType COMBO
	#
	set cmbXPath "${xpath}/Item\[@id='FunctionType'\]"
	set comboList [::xmlutils::getValues $xml $cmbXPath]
	set cmbPid [::xmlutils::getAttribute $xml $cmbXPath pid]
	set cmbDv [::xmlutils::getAttribute $xml $cmbXPath dv]
	set cmbTooltip [::xmlutils::getAttribute $xml $cmbXPath tooltip]
	
	grid [ttk::label $f.lblFunctionType -text [= $cmbPid]:] \
	-row 3 -column 0 -padx 3 -pady 2 -sticky nw -in $f
	
	grid [ttk::combobox $f.cmbFunctionType -values $comboList -state readonly -width $width -textvariable "::KFun::cmbFunctionType"] \
	-row 3 -column 1 -padx 3 -pady 2 -sticky nw -in $f
	tooltip::tooltip $f.cmbFunctionType [= $cmbTooltip]
	
	::xmlutils::setComboValue $xml $cmbXPath $f.cmbFunctionType $cmbDv
	
	#
	#  VARIABLES
	#
	grid [ttk::label $f.lblVars -text [= "Kratos variables"]: ] \
	-row 4 -column 0 -pady 3 -pady 2 -sticky nw -in $f
	
	grid [ttk::combobox $f.cmbVars -values [::KFun::getComboVars] -state readonly -width $width -textvariable "::KFun::cmbVars"] \
	-row 4 -column 1 -padx 3 -pady 2 -sticky nw -in $f
	tooltip::tooltip $f.cmbVars [= "Kratos variables"]:
	#Inicializamos la variable del valor seleccionado
	set ::KFun::cmbVars ""
	
	grid [ttk::button $f.addVars -text [= "Add variable"] -command "::KFun::insertText ${f}.ftext.txtEditArea ::KFun::cmbVars" ] \
	-row 4 -column 1 -sticky nw  -pady 2 -padx [expr $width + 180] -in $f
	tooltip::tooltip $f.addVars [= "Add variable"]
	set img [::WinUtils::GetImage "add.gif" ]
	if { $img != -1 } { $f.addVars configure -image $img }
	
	
	#EDITION AREA
	set ft "$f.ftext"
	grid [ttk::frame $ft -width [expr $width + 60] -height 10 ] \
	-row 5 -column 0 -padx 5 -columnspan 4 -sticky nsw -in $f
	
	grid [ text $ft.txtEditArea -yscrollcommand "$ft.srl_y set" -endline 15 \
		-cursor ibeam -width [expr $width + 40] -height 12 -wrap word] \
	-row 0 -column 0 -in $ft
	
	grid [ scrollbar $ft.srl_y -command  "$ft.txtEditArea yview"  ] \
	-row 0 -column 1 -sticky ns -in $ft 
	
	$ft.txtEditArea delete 0.0 end
	#$ft.txtEditArea insert 0.0 [::KFun::nodeDefinition $id]
	
	set cmbXPath "${xpath}/Item\[@id='Definition'\]"
	$ft.txtEditArea insert 0.0 [::xmlutils::getAttribute $xml $cmbXPath sourceCode]
	

	
	#OK CANCEL
  #$f.bBottomOk and $f.bBottomCancel names need to popup click ok cancel before change item
	grid [ttk::button $f.bBottomOk -text [= "Ok"]  -command "::KFun::acceptFunction $ft $T" ] \
	    -row 6 -column 0 -sticky nw -columnspan 2 -pady 3 -padx 30  -in $f
	tooltip::tooltip $f.bBottomOk [= "Add this function"]
	
	grid [ttk::button $f.bBottomCancel -text [= "Cancel"]  -command "::KFun::cancelBottom" ] \
	    -row 6 -column 0 -sticky nw -columnspan 2 -pady 3 -padx 120  -in $f
	tooltip::tooltip $f.bBottomCancel [= "Cancel process"]
	
	if { $::KFun::applyFun } {
		
		grid remove ${::KFun::WinPath}.tbf.applyFun
	}
}

proc ::KFun::acceptFunction { f T } {
	
	global KPriv
	set xml $KPriv(xmlDocFun)
	
	#Caso común editar/insertar
	if {$::KFun::nomId == ""} {
		WarnWin [= "Please, choose a non-empty ID.'"]
		return ""
	}
	
	
	set item [$T selection get]
	
	if { $item == ""} {
		
		#Caso INSERTAR nueva función
		
		#Comprobamos que el id no esté repetido (aun no exista en el xml)
		set id [::KFun::encodeId $::KFun::nomId]
		set xpath "Kratos_Functions/Functions/Function\[@id='$id'\]"
		set node [$xml selectNodes $xpath]
		if {$node != ""} {
		        WarnWin "Ya existe una función con este id:\n'$::KFun::nomId'"
		        return ""
		}
		
		set functionTxt [::xmlutils::AsXml [$f.txtEditArea get 0.0 end]]
		
		set nodeWhere [$xml selectNodes "Kratos_Functions/Functions"]
		
		set nodeWhat [$xml selectNodes "Kratos_Functions/Templates/Function"]
		
		set newNode [$nodeWhat cloneNode -deep]
		
		$newNode setAttribute id $id
		$newNode setAttribute pid $::KFun::nomId
		
		$nodeWhere appendChild $newNode
		
	} else {
		
		#Caso EDITAR función
		
		#Comprobamos que el id no esté repetido (aun no exista en el xml)
		set id [::KFun::encodeId $::KFun::nomId]
		set idAnterior [$T item tag names $item]
		set xpath "Kratos_Functions/Functions/Function\[@id='$id'\]"
		
		if {$id != $idAnterior } {
				
		        set node [$xml selectNodes $xpath]
		        if {$node != ""} {
		                WarnWin "Ya existe una función con este id:\n'$::KFun::nomId'"
		                return ""
		        } else {
		                #Actualizamos el nuevo Id en el XML
		                set node [$xml selectNodes "Kratos_Functions/Functions/Function\[@id='$idAnterior'\]"]
		                $node setAttribute id $id
		                $node setAttribute pid $::KFun::nomId
		        }
		        
		        ::KMProps::checkFunctions $idAnterior $id
		}
		
		set functionTxt [::xmlutils::AsXml [$f.txtEditArea get 0.0 end]]
	}
	
	set itemXPath "${xpath}/Item\[@id='Definition'\]"
	::xmlutils::getAttribute $xml $itemXPath sourceCode $functionTxt
	#::KFun::nodeDefinition $id "write" $functionTxt
	
	#Actualizar los valores de los combos
	set f "${::KFun::WinPath}.fBottom"
	
	set itemXPath "${xpath}/Item\[@id='AppliedTo'\]"
	set fcmb "$f.cmbAppliedTo"
	set selCombo [::xmlutils::getComboValue $xml $itemXPath $fcmb]
	::xmlutils::getAttribute $xml "$itemXPath" dv $selCombo
	
	set itemXPath "${xpath}/Item\[@id='FunctionType'\]"
	set fcmb "$f.cmbFunctionType"
	set selCombo [::xmlutils::getComboValue $xml $itemXPath $fcmb]
	::xmlutils::getAttribute $xml "$itemXPath" dv $selCombo
	
	#msg "id:$id ipath:$itemXPath f:$fcmb"
	
	::KFun::clearFields
	
	::KFun::refreshTree 
	
	if { $item != "" } {
		$T item expand $item
	}
	
	#$nodeWhere appendXML "<Function id=\"$id\" pid=\"$::KFun::nomId\" icon=\"\" help=\"\"> \
	#                <Item id=\"FunctionType\" pid=\"\" dv=\"Initial\" values=\"Initial,Update,Postprocess\" ivalues=\"Initial,Update,Postprocess\" ></Item> \
	#                <Item id=\"Definition\" pid=\"Definition\">$functionTxt</Item> \
	#                <Item id=\"AppliedTo\" pid=\"AppliedTo\" dv=\"onNodes\" values=\"onNodes,onElements,onFaces,onConditions\" ivalues=\"onNodes,onElements,onFaces,onConditions\"></Item> \
	#                </Function>"
}

proc ::KFun::OpenNotes { } {
	
    package require texteditor
    set modelname [GiD_Info Project ModelName]
    if  { $modelname == "UNNAMED" } {
	#WarnWin [_ "A project title is needed. Save project to get it"]
	set filename ""
    } else {
	set name [file tail $modelname]
	set filename [file join $modelname.gid $name.txt]
    }
    set w .gid.texteditor
    if { $filename != "" } {
	if { [file exists $filename] } {
	    TextEditor::Create $filename "" "" 1 $w utf-8
	} else {
	    TextEditor::Create "" "" "" 1 $w utf-8
	    set TextEditor::currentfile $filename ;#to save with this name
	    wm title $w $filename
	}
    } else {
	TextEditor::Create "" "" "" 1 $w utf-8
	TextEditor::ReadFromVariable
	set TextEditor::cansavetovariable 1 ;#to temporary store in a tcl variable
    }
}

#proc ::KFun:CloseNotes { } {
#    GidUtils::CloseWindow NOTES
#}

proc ::KFun::CreateNewFunction { } {
	
	set T $::KFun::TreePath
	
	::KFun::buildEditor $T 
	
	$T selection clear
	
	::KFun::clearFields
	
}

proc ::KFun::DeleteFunction { } {

	global KPriv
	
	set T $::KFun::TreePath
	
	set item [$T selection get]
	
	if { $item == ""} {
		WarnWin [= "No function selected to remove."]
		return ""
	}
	set id [$T item tag names $item]
	set node [$KPriv(xmlDocFun) selectNodes "Kratos_Functions/Functions/Function\[@id='$id'\]"]
	if { $node != "" } {
		
		set changes [::KMProps::checkFunctions $id "" "justCheck"]
		if { $changes > 0 } {
			if {$changes == 1} { 
				set plural ""
				set person "it"
			} else {
				set plural "s"
				set person "them"
			}
			set aviso [= "This function is assigned %s time%s in this project,\ndo you want to unassign %s and continue?" $changes $plural $person]
			set answer [::WinUtils::confirmBox "." "$aviso"]
			if { $answer == "ok" } {
				set changes [::KMProps::checkFunctions $id]
				$node delete
			}
		} else {
			$node delete
		}
		::KFun::refreshTree
	}
	
}

proc ::KFun::clearFields { } {
	
	set w "$::KFun::WinPath"
	${w}.fBottom.ftext.txtEditArea delete 0.0 end
	
	set ::KFun::nomId ""
	#set ::KFun::cmbAppliedTo ""
	#set ::KFun::cmbFunctionType ""
}

proc ::KFun::refreshTree { } {
	
	set T $::KFun::TreePath
	 foreach item [$T item range 0 end] {
	    if { $item != 0 } {
		    $T item remove $item
	    }
    }
    ::KFun::fillFunctions $T
}



#
# Contruye un combo dinámico en el item pulsado del árbol y posteriormente se elimina
#
proc ::KFun::buildFrame { T item } {
	
	global KPriv
	
	set xml $KPriv(xmlDocFun)
	set parent [ $T item parent $item ]
	set parentId [$T item tag names $parent]
	set itemId [$T item tag names $item]
	
	set xpath "Kratos_Functions/Functions/Function\[@id='$parentId'\]/Item\[@id='$itemId'\]"
	
	set pid [::xmlutils::getAttribute $xml $xpath pid]
	set dv [::xmlutils::getAttribute $xml $xpath dv]
		
	set tooltip [::xmlutils::getAttribute $xml $xpath help]
	
	set comboList [::xmlutils::getValues $xml $xpath]
	
	#---------------------------#---------------------------#
	# Configurar frame en función del XML
	#---------------------------#---------------------------#
	set bg "#F8F8F8"
	set idFrame "f${parentId}/${itemId}"
	set f "$T.[list $idFrame]"
	
	if { [winfo exists $f] } {
		destroy $f
	}
	set f [frame "$f" -borderwidth 0 -background $bg]
	
	if { [llength $comboList] > 0 } {
		grid [ttk::combobox $f.cmb -values $comboList -state readonly -width 30 -textvariable "::KFun::cmb$idFrame"] \
		-row 0 -column 0 -padx 3 -sticky ne -in $f
		
		::xmlutils::setComboValue $xml $xpath $f.cmb $dv
		
		bind $f.cmb <<ComboboxSelected>> "::KFun::cmbSelectChange $T $xml $item 1 current"
		
	} else {
		grid [ttk::combobox $f.cmb -state normal -textvariable "::KFun::cmb$idFrame"] \
		-row 0 -column 0 -padx 3 -sticky nw -in $f
		
		set ::KFun::cmb$idFrame $dv
		
	    bind $f.cmb <Leave> [list ::KFun::cmbSelectChange $T $xml $item 0 current]
	    bind $f.cmb <FocusOut> [list ::KFun::cmbSelectChange $T $xml $item 1 current]
	    bind $f.cmb <Escape> [list ::KFun::cmbCancel $T $xml $xpath $item $idFrame]
	}
	# Si pulsan intro o Esc también forzamos la salida del combo
	bind $f.cmb <Return> [list ::KFun::cmbSelectChange $T $xml $item 1 current]
	bind $T <Escape> [list ::KFun::cmbCancel $T $xml $xpath $item $idFrame]
	
	tooltip::tooltip $f.cmb "$tooltip"
	
	return $f
}

proc ::KFun::cmbCancel { T xml xpath item idFrame  } {
	
	set id [::xmlutils::getAttribute $xml $xpath id]
	set pid [::xmlutils::getAttribute $xml $xpath pid]
	set dv [::xmlutils::getAttribute $xml $xpath dv]
	set dvText [::xmlutils::getValueText $xml $xpath]
	
	set ::KFun::cmb$idFrame "$dv"
	
	set f "$T.f$idFrame"
	if {[winfo exists $f]} {
		destroy $f
	}
	#Miramos si tiene algun estilo especial
	$T item style set $item C0 styAnyRead
	$T item element configure $item C0 elemTxtRead -text "$pid: $dvText"
	
	set ::KFun::lastSelected {}
}

proc ::KFun::cmbSelectChange { T xml item {remove 1} {selectVal "current"}} {
	
	set refresh 0
	
	global KPriv
	set xml $KPriv(xmlDocFun)
	set parent [ $T item parent $item ]
	set parentId [$T item tag names $parent]
	set itemId [$T item tag names $item]
	
	set idFrame "f${parentId}/${itemId}"
	
	set xpath "Kratos_Functions/Functions/Function\[@id='$parentId'\]/Item\[@id='$itemId'\]"
	
	set id [::xmlutils::getAttribute $xml $xpath id]
	
	set comboState [::xmlutils::getComboState $xml $xpath]
	if { $comboState == "normal" } {
		
		set selCombo [set ::KFun::cmb$idFrame]
		set selComboText $selCombo
		if {$selectVal != "current"} {
		        set ::KFun::lastSelected {}
		}
	} else {
		
		if {$selectVal == "current"} {
		        
		        set f "$T.[list ${idFrame}].cmb"
		        set selCombo [::xmlutils::getComboValue $xml $xpath $f]
		        set selComboText [::xmlutils::getComboValue $xml $xpath $f "text"]
		        set ::KFun::lastSelected [list $item $selCombo $selComboText]
		        
		} else {  
		        
		        set selCombo [lindex $::KFun::lastSelected 1]
		        set selComboText [lindex $::KFun::lastSelected 2]
		        set ::KFun::lastSelected {}
		}
	}
	
	set pid [::xmlutils::getAttribute $xml $xpath pid]
	
	#msg "id:$id rem:$remove sel:$selCombo"
	if { $remove } {
		
		#Eliminamos el frame
		set f "$T.f$idFrame"
		if {[winfo exists $f]} {
		        destroy $f
		}
		
		#Miramos si tiene algun estilo especial
		$T item style set $item C0 styAnyRead
		$T item element configure $item C0 elemTxtRead -text "$pid: $selComboText"
	}
	
	#Vuelve a poner la imagen por si no recargamos el árbol
	set icon [::xmlutils::getAttribute $xml $xpath icon]
	set imagen [::WinUtils::GetImage $icon]
	if { $imagen != -1 } {
		$T item image $item C0 $imagen
	}
	
	#Miramos si ha cambiado para actualizar en consecuencia
	if { [::xmlutils::getAttribute $xml $xpath dv] != $selCombo } {
		
		#Guarda el nuevo valor en el xml
		::xmlutils::getAttribute $xml $xpath dv $selCombo
		
		#Si el combo tiene un class especial habrá q reconstruir el árbol
		set clase [::xmlutils::getAttribute $xml $xpath class]
		if { $clase != "" } {
		        
		        #Volvemos a cargar el árbol de propiedades
		        set refresh 1
		}
	}
	
	if {$refresh} {
		::KFun::refreshTree $T
	}
}
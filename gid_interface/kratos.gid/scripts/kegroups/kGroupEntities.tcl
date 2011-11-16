##############################################################################################
#
#	NAME: listEntities.tcl
#
#	PURPOSE: List Entities for each group 
#
#	QUANTECH ATZ-DEVELOPMENT DEPARTMENT
#
#	AUTHORS : Luis Calvo, Gerardo Socorro
#
#	CREATED AT: 20/05/10
#
#	HISTORY:
#
#	 0.2- 22/06/11-G. Socorro, delete snit, tdom and xmlstruct from the package require
#        0.1- 22/06/11-L. Calvo, create the base source code
##############################################################################################

package require treectrl
package require tooltip
package provide KEGroups 1.0 

# Create a base namespace KEGroups
namespace eval ::LEntities:: {
    
    # Path of the base window 
    variable WinPath ".gid.kegroups.listEntities"
    # The tree path
    variable TreePath   
    variable lSort 0
}

proc ::LEntities::InitBaseWindow { } {
    
    
    global KPriv
    
    set w "$::LEntities::WinPath"

    # Open the window outside
    ::LEntities::OpenWindowOutside $w

    if {![winfo exists $w]} return 

    set T [::LEntities::CreateTreeAndToolbar $w]
    
    ::LEntities::listEntities $T
    
    # Binding
    bind $w <Alt-c> "destroy $w"
    bind $w <Escape> "destroy $w"
    
}

proc ::LEntities::CreateTreeAndToolbar { w } {
    
    variable TreePath   


    # Create the treectrl properties
    set mdf [ttk::frame $w.middle]
    set T [::LEntities::CreateTreeProperties $w]
    # Set the tree path to a namespace variable
    set TreePath $T
    
    grid $w.middle -sticky wes
    
    
    #grid [ttk::button $tbf.bClose -text [= "Close"] -command [list destroy $w]]  -row 1 -column 3 -sticky e -padx 5 -pady 3
    #grid $tbf -sticky ews
    #grid anchor $tbf.bClose center
    
    focus $T
    
    return $T
}

proc ::LEntities::ClickTree { x y T } {
    
    set info [$T identify $x $y]
    
    if { [lindex $info 0] == "item" && [llength $info] >= 4 } {
	
	set item [lindex $info 1]
	set col [lindex $info 3]
	
	if {[$T item numchildren $item] == 0} {
	    
	    set parentEntity [$T item text [$T item parent $item] 0]
	    set entityId [$T item text $item 0]
	    
	    # Le hacemos un signal a la entidad seleccionada en el árbol
	    GiD_Process render Normal
	    GiD_Process MEscape Utilities SignalEntities $parentEntity $entityId
	    
	    #Esto sería para que solo se muestre la entidad un segundo
	    #GiD_Process Escape 
	} else {
	    $T item toggle $item
	}
	
	
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


proc ::LEntities::CreateTreeProperties {w} {

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
    set collistid [list [= "Group Id"]]
    
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
    $T element create elemTxtRead text -fill [list $::KEGroups::SystemHighlightText {selected focus}] -lines 1
    $T element create elemRectSel rect -fill [list $::KEGroups::SystemHighlight {selected focus} gray {selected !focus}] -showfocus yes
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
    #$T notify bind DragTag <Drag-receive> { ::LEntities::ReceiveDragGroups  %T %l %I }
    #$T notify bind EditTag <Edit-accept> { ::LEntities::SetGroupsToRename %T %I %t }

    bind $T <Button-1> [list ::LEntities::ClickTree %x %y $T]
    #	 bind $T <Double-Button-1> [list DoubleClickTableList_New %x %y $T]
    #	 bind $T <Return> [list SetLayersTo TOUSE $T]
    #	 bind $T <Key-Delete> [list SetLayersToDelete $T]
    #	 bind $T <Alt_L> [list InvertSelectionTableList $T]
    #	 bind $T <Alt_R> [list InvertSelectionTableList $T]
    #	 bind $T <Meta_L> [list InvertSelectionTableList $T]
    #	 bind $T <Meta_R> [list InvertSelectionTableList $T]
    #	  bind $T <F2> [list ::LEntities::BeginEditGroups $T]

    #bind $T <Button-3> "[list ::LEntities::MenuContextualGroup %W %x %y] ; break"

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

proc ::LEntities::ReceiveDragGroups { T dragged_list dst } {
    
}

proc ::LEntities::OpenWindowOutside { w } {
    
    if { [winfo exists $w] } {
	destroy $w
    }

    # Init the window
    set title [= "View group entities "]
    InitWindow $w $title wListEntitiesGeom ::LEntities::InitBaseWindow
    
    return $w
}


proc ::LEntities::InsertNewItem { T text {icon ""} {parentitem root} {childs 1}} {
    
    if { $childs } {
	
	set item [$T item create -button yes]
	$T item lastchild $parentitem $item
	$T item style set $item C0 styAnyRead
	$T item element configure $item C0 elemTxtRead -text $text
	
    } else {
	
	set item [$T item create -button no]
	$T item lastchild $parentitem $item
	$T item style set $item C0 styAnyRead
	$T item element configure $item C0 elemTxtRead -text $text
    }
    
    if {$icon != ""} {
	set imagen [::WinUtils::GetImage $icon]
	if { $imagen != -1 } {
	    $T item image $item C0 $imagen
	}
    }
    
    #$T item style set $item C1 styAnyRead
    #$T item element configure $item C1 elemTxtRead -text $id
    
    return $item
}

#
# 
#
proc ::LEntities::listEntities { T } {
    
    global KPriv
    
    set Tgroups $::KEGroups::TreePath
    
    set itemGroups [$Tgroups selection get]
    
    if {[llength $itemGroups] == 0 } {
	
	#Si no hay ninguno 
	set groups $KPriv(groupsId)
    } else {
	
	set groups {}
	foreach g $itemGroups {
	    lappend groups [$Tgroups item text $g 0]
	}
    }
    
    set w "${::LEntities::WinPath}.progress"
    
    ::LEntities::startProgress $w
    
    set whatuse [GiD_Info Project ViewMode]	
    
    set ::LEntities::totalEntities 0

    foreach groupId $groups {
	
	set gItem [::LEntities::InsertNewItem $T $groupId]
	#set entities [::KEGroups::getAssignedGiDEntities $groupId]
	
	update idletasks
	
	switch $whatuse {
	    
	    GEOMETRYUSE {
		set max [expr [llength $groups] ]
		$w.pro configure -maximum $max
		#set ::LEntities::totalEntities [llength $geomEntities]
		foreach entity [list point line surface volume] title [list Points Lines Surfaces Volumes] {
		    
		    
		    set entityList [::KEGroups::getGroupGiDEntities $groupId $entity]
		    #msg $entityList
		    if { [llength $entityList]} {
			#msg "i:$gItem g:$groupId: e:$entity t:$title"
			::LEntities::insEntity $T $entityList "$title" $gItem "${entity}16x16.gif"
		    }
		    
		    #::LEntities::insEntity $T $points "Points" $gItem "point16x16.gif"
		    #::LEntities::insEntity $T $lines "Lines" $gItem "line16x16.gif"
		    #::LEntities::insEntity $T $surfaces "Surfaces" $gItem "surface16x16.gif"
		    #::LEntities::insEntity $T $volumes "Volumes" $gItem "volume16x16.gif"
		}	
	    }
	    MESHUSE {
		set max [expr [llength $groups] * 2]
		$w.pro configure -maximum $max
		
		foreach entity [list node element] title [list Nodes Elements] {	
		    
		    set entityList [::KEGroups::getGroupGiDEntities $groupId "${entity}s"]	
		    if { [llength $entityList]} {
			
			::LEntities::insEntity $T $entityList "$title" $gItem "${entity}16x16.gif"
		    }
		}
		#set entityList [::KEGroups::getGroupGiDEntities $groupId point]	
		#if { [llength $entityList]} {
		#	
		#	
		#	::LEntities::insEntity $T $entityList "Nodes" $gItem "node16x16.gif"
		#}
		#
		#set entityList {}
		#set i 0
		#foreach entity [list line surface volume] {	
		#	
		#	set list [::KEGroups::getGroupGiDEntities $groupId $entity]
		#	if { [llength $list] } {
		#	 	lappend entityList 
		#	}
		#}
		#
		#	set entityList [::KUtils::FlatEmbeddedLists $entityList]
		#	#msg "entityList $entityList"
		#	
		#	if { [llength $entityList]} {
		#		
		#		::LEntities::insEntity $T $entityList "Elements" $gItem "element16x16.gif"
		#	}
		#
	    }
	}
	
	if {[$T item numchildren $gItem] == 0} {
	    $T item remove $gItem
	}
    }
    
    ::LEntities::stopProgress $w
    
    #$qdb::flic.p1 configure -value $progress
    #set progress [expr $progress + $max]
}

proc ::LEntities::insEntity { T array title parent icon} {
    
    
    if {[llength $array] > 0 } {
	
	dict create word val
	foreach entity $array {
	    dict set word $entity val ""
	}
	
	if { [dict size $word] > 0 } {
	    
	    set eItem [::LEntities::InsertNewItem $T "$title" "$icon" $parent]
	    
	    if { $::LEntities::lSort } {
		
		#Ordenamos y luego recorremos el diccionario insertando en el árbol
		set lista {}
		dict for {id nada} $word {
		    lappend lista $id
		}
		set lista [lsort -integer $lista]
		foreach entity $lista {
		    ::LEntities::InsertNewItem $T "$entity" "" $eItem 0
		}
		$T item collapse $eItem
	    } else {
		
		#Recorremos el diccionario sin ordenar
		dict for {id nada} $word {
		    ::LEntities::InsertNewItem $T "$id" "" $eItem 0
		    
		}
		$T item collapse $eItem
	    }
	}
	
    }
}

proc ::LEntities::insEntity2 { T array title parent icon} {
    
    if {[llength $array] > 0 } {
	set i 0
	set eItem [::LEntities::InsertNewItem $T "$title" "$icon" $parent]
	foreach entity $array {
	    if {$entity != ""} {
		::LEntities::InsertNewItem $T "$entity" "" $eItem 0
	    }
	    #if {$i % 20 == 0 } {
	    #	::LEntities::totalEntities
	    #}
	    
	    #set i [expr $i + 1]
	    #set ::LEntities::totalEntities [expr $::LEntities::totalEntities + 1]
	}
	$T item collapse $eItem
    }
}


proc ::LEntities::startProgress { w } {
    
    set pf [ttk::frame $w ]
    
    # Barra de progreso cuando 
    grid [ttk::progressbar $pf.pro -mode determinate -length 50] \
	-row 0 -column 1 -columnspan 2 -pady 2 -sticky swe -in $pf
    
    grid $pf -row 2 -column 0 -sticky wes
    
    # Configura la barra de progreso en función del número de entidades a mostrar
    #set maximum [expr 100 * ([llength $qdb::inputList])+1]
    #$qdb::flic.p1 configure -maximum $maximum
    #set max [expr $maximum / [llength $qdb::inputList]]
    #set progress $max
    
    grid $pf.pro
    $pf.pro start
}

proc ::LEntities::stopProgress { w } {
    
    $w.pro stop
    destroy $w.pro
}
##############################################################################
#
#    NAME: menus.tcl
#
#    PURPOSE: TCL script to work with menus in preprocess and postprocess in the 
#             Kratos problem type
#
#    QUANTECH ATZ-DEVELOPMENT DEPARTMENT
#
#    AUTHOR : G. Socorro
#
#    CREATED AT: 29/03/06
#
#    LAST MODIFICATION : add more options to the Kratos preprocess toolbar
#
#    VERSION : 0.2
#
#    HISTORY:
#
#     0.2-08/06/10-G. Socorro, add more options to the Kratos preprocess toolbar
#     0.1-01/02/10-G. Socorro, create a base source code
#
###############################################################################

# kmtb => Kratos Menus and ToolBars

namespace eval ::kmtb:: {

}

proc ::kmtb::ChangePreprocessMenu {dir} {
   
    # Add some menu in preprocess menu
    ::kmtb::AddMenuToPreprocessMenu $dir
    
    # Create the process toolbar
    ::kmtb::CreatePreprocessModelTBar $dir
    
    # Update the menu properties
    ::GiDMenu::UpdateMenus
}

proc ::kmtb::AddMenuToPreprocessMenu {dir} {
    global NumMenus MenuNames MenuEntries MenuCommands MenuAcceler MenuIcons
    
    set mname [= "Kratos#C#menu"]
    
    # Insert a interface name entry after meshing menu
    set pos [lsearch -exact $MenuNames [_ "Mesh#C#menu"]]
    set pos2 [lsearch -exact $MenuNames $mname]
    ;# WarnWinText "pos:$pos pos2:$pos2"
    if {($pos!="-1") && ($pos2==-1)} {
     	set MenuNames [linsert $MenuNames [incr pos] $mname]
     	set NumMenus [llength $MenuNames ]
	;# WarnWinText "MenuNames:$MenuNames NumMenus:$NumMenus"
     	# Move all the predefined commands one position
     	for {set ii $NumMenus} {$ii > $pos} {incr ii -1} {
	    set jj [expr $ii - 1]
	    ;# WarnWinText "ii:$ii jj:$jj"
     	    catch { set MenuEntries($ii)  $MenuEntries($jj)  }
     	    catch { set MenuCommands($ii) $MenuCommands($jj) }
     	    catch { set MenuAcceler($ii)  $MenuAcceler($jj)  }
	    catch { set MenuIcons($ii) $MenuIcons($jj)  }
     	    set MenuEntries($jj)  ""
     	    set MenuCommands($jj) ""
     	    set MenuAcceler($jj)  ""
	    set MenuIcons($jj) ""
     	}
    }

    # Disabled all base options
    GidChangeDataLabel "Data units" ""
    GidChangeDataLabel "Interval" ""
    GidChangeDataLabel "Conditions" ""
    GidChangeDataLabel "Materials" ""
    GidChangeDataLabel "Intervals" ""
    GidChangeDataLabel "Problem Data" ""
    GidChangeDataLabel "Local axes" ""
    
    # Create the interface name for workshop menu
    set MenuEntries($pos) [list [= "Entities group editor#C#menu"]... \
			       --- \
			       [= "Model properties#C#menu"] \
			       [= "Material database#C#menu"] \
			       --- \
			       [= "Function editor#C#menu"] \
			       [= "Project setting#C#menu"]]
    
    set MenuCommands($pos) [list [list -np- ::KEGroups::InitBaseWindow] \
				"" \
				[list -np- ::KMProps::InitBaseWindow] \
				[list -np- ::KMProps::InitBaseWindow Materials] \
				"" \
				[list] [list]]
    
    
    set MenuAcceler($pos) {"" "" "" "" ""}
    set MenuIcons($pos) {"" "" "" "" ""}

}

proc ::kmtb::DeleteMenu {pos position {prepost "PRE"}} {
    # For base global  
    global MenuEntries MenuCommands MenuAcceler MenuIcons

    # Delete the menu entry and command
    # WarnWinText "MenuEntries:$MenuEntries($pos)\npos:$pos\nposition:$position"
    # WarnWinText "MenuCommands:$MenuCommands($pos)\npos:$pos\nposition:$position"
    # WarnWinText "MenuAcceler:$MenuAcceler($pos)\npos:$pos\nposition:$position"
    # WarnWinText "MenuIcons:$MenuIcons($pos)\npos:$pos\nposition:$position"
    if {[info exists MenuEntries($pos)]} {
	set llen [expr [llength $MenuEntries($pos)]-1]
	# WarnWinText "MenuEntries length $llen"
	if {$position<=$llen} {
	    set MenuEntries($pos) [lreplace $MenuEntries($pos) $position $position]
	}
    }
    if {[info exists MenuCommands($pos)]} {
	set llen [expr [llength $MenuCommands($pos)]-1]
	# WarnWinText "MenuCommands length $llen"
	if {$position<=$llen} {
	    set MenuCommands($pos) [lreplace $MenuCommands($pos) $position $position]
	}
    }
    if {[info exists MenuAcceler($pos)]} {
	set llen [expr [llength $MenuAcceler($pos)]-1]
	# WarnWinText "MenuAcceler length $llen"
	if {$position<=$llen} {
	    set MenuAcceler($pos) [lreplace $MenuAcceler($pos) $position $position]
	}
    }
    if {[info exists MenuIcons($pos)]} {
	set llen [expr [llength $MenuIcons($pos)]-1]
	# WarnWinText "MenuIcons length $llen"
	if {$position<=$llen} {
	    set MenuIcons($pos) [lreplace $MenuIcons($pos) $position $position]
	}
    }
    # WarnWinText "Despues MenuEntries:$MenuEntries($pos)\npos:$pos\nposition:$position"
}

proc ::kmtb::DeleteSubMenu {pos position {prepost "PRE"}} {
     # For base global  
     global MenuEntries MenuCommands MenuAcceler MenuIcons
    
     if {[info exists MenuEntries($pos,$position)]} {
	 # Reset submenu
	 array unset MenuEntries $pos,$position
     }
     if {[info exists MenuCommands($pos,$position)]} {
	 array unset MenuCommands $pos,$position   
     }
     if {[info exists MenuAcceler($pos,$position)]} {
	 array unset MenuAcceler $pos,$position 
     }
     if {[info exists MenuIcons($pos,$position)]} {
	 array unset MenuIcons $pos,$position
     }
 }

proc ::kmtb::AddToMenu {pos clist {prepost "PRE"}} {
     # For base global  
     global MenuEntries MenuCommands MenuAcceler MenuIcons

     set translationfunc =

     # try to add ---
     if {[lindex $MenuEntries($pos) end] !="---"} {
	 lappend MenuEntries($pos) "---"
	 lappend MenuCommands($pos) ""
	 lappend MenuAcceler($pos) ""
	 lappend MenuIcons($pos) ""
     }
 }

proc ::kmtb::CreatePreprocessModelTBar {dir {type "DEFAULT INSIDELEFT"}} {
    global KBitmapsNames KBitmapsCommands KBitmapsHelp 
    global KPriv
    
    catch { unset KBitmapsNames KBitmapsCommands KBitmapsHelp }
    
    set KBitmapsNames(0) "images/groups.gif images/new_props.gif images/maticon.gif --- images/openrunsim.gif images/runsimulation.gif images/runsiminfo.gif images/stop.gif"
    
    set KBitmapsCommands(0) [list \
				 [list -np- ::KEGroups::InitBaseWindow] \
				 [list -np- ::KMProps::InitBaseWindow] \
				 [list -np- ::KMProps::InitBaseWindow Materials] \
				 "" \
				 [list -np- RunWin] \
				 {Utilities Calculate} \
				 [list -np- PWViewOutput] \
				 {Utilities CancelProcess} ]

    set KBitmapsHelp(0) [list [= "Define the group properties using the group editor"] \
			     [= "Define the model properties"] \
			     [= "Define the material properties"] \
			     "" \
			     [= "Open the process control window"] \
			     [= "Run the simulation"] \
			     [= "View process info"] \
			     [= "Cancel process"]]	     

    # prefix values:
    # Pre        Only active in the preprocessor
    # Post       Only active in the postprocessor
    # PrePost    Active Always
    
     set prefix "Pre"
     set name "KPreprocessModel"
     set winname "::kmtb::CreatePreprocessModelTBar"
     set KPriv(ToolBars,PreprocessModelTBar) [CreateOtherBitmaps ${name}bar \
						   [= "Model definition toolbar"] \
						   KBitmapsNames KBitmapsCommands \
						   KBitmapsHelp $dir "$winname [list $dir]" \
						   $type $prefix]
    
    AddNewToolbar "${name}bar" ${prefix}${winname}WindowGeom "$winname [list $dir]" \
	[= "Model definition toolbar"]
	
}

proc ::kmtb::EndCreatePreprocessTBar {} {
    global KPriv
   
    set name "KPreprocessModel"
    set winname "::kmtb::CreatePreprocessModelTBar"
    catch { 
	ReleaseToolbar "${name}bar"
	rename $winname ""
	destroy $KPriv(ToolBars,PreprocessModelTBar)
	unset KPriv(ToolBars,PreprocessModelTBar)
    }
}

proc ::kmtb::UpdateVersion {} {
    return [GiD_Process MEscape data defaults TransfProblem]
}


 
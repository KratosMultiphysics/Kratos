###############################################################################
#
#    NAME: kmpropswin.tcl
#
#    PURPOSE: Utilities procedures to manage the main kratos model window
#
#    COPYRIGHT
#  
#    QUANTECH ATZ-DEVELOPMENT DEPARTMENT
#
#    AUTHOR : G. Socorro
#
#    CREATED AT: 29/03/12
#
#    HISTORY:
# 
#     0.6-  17/06/13- A. Melendo, layer / tree inside windows, redimension problems solved
#     0.5-  01/10/12- J. Garate, Enable/disable Curves Module
#     0.4-  20/09/12- J. Garate, Adaptation for New Kratos Interface Version, including Curves support
#     0.3-  02/04/12- G. Socorro, correct a bug when the layer and the kratos windows are inside
#     0.2-  30/03/12- G. Socorro, add some new procedures CreatePropertiesTabs, StartBaseWindow,etc.
#     0.1-  29/03/12- G. Socorro, create the base source code from KMProps namespace
#
###############################################################################
#                      Procedures that belong to this file
###############################################################################
#         Name                      |        Functionality
#------------------------------------------------------------------------------
# 1. SwitchWindowInsideOutside      | Switch the model window between inside or outside
# 2. OpenWindowOutside              | Open the model window outside
# 3. CloseWindowOutside             | Close the model window outside
# 4. ChangeKratosModelWindow        | Change the kratos model window
# 5. OpenWindowInside               | Procedure to open window inside, include creation of panedwindow to manage several windows
# 6. WriteKratosModelGeomToVar      | Procedure to write to a GidPriv global variable the geometry of the model window
# 7. CreateWindow                   | Procedure to create the model and material window
# 8. CloseWindowInside              | This procedure take account new panedwindow GidPriv(pwCentral) to manage several windows
# 9. CreatePropertiesFrames         | Procedure to create the properties frames
# 10. CreateCaption                 | Procedure to create the caption frame
# 11. CreatePropertiesTabs          | Procedure to create the properties tabs (model and materials)
# 12. StartBaseWindow               | Start the kratos base window
# 13.

proc ::KMProps::StartBaseWindow {{whattab "Model"} {what "1"} } {
    # ABSTRACT: Start the kratos base window
    # Arguments:
    # whattab => The tab to be raised
    # what    => Option ["1:Inside"|"0:Outside"]
    variable winpath; variable Layout

    # Init some namespace variables
    ::KMProps::Init
    ::KMat::Init
    if { [kipt::CurvesModule ] } {
	::KCurves::Init
    }
    set w "$winpath"
    # wa "w:$w whattab:$whattab what:$what"
    
    # Open the window
    if {[winfo exists $w]} {
	if {[info exists Layout]} {
	    # wa "Layout:$Layout"
	    if {$Layout eq "OUTSIDE"} {
		set what 0
		::KMProps::CreateWindow $w $what $whattab
	    } else {
		set what 1
		::KMProps::CreateWindow $w $what $whattab
	    }
	}
    } else {
	::KMProps::CreateWindow $w $what $whattab
    }
}

proc ::KMProps::CreatePropertiesFrames { w {whattab "Model"}} {
    # ABSTRACT: Procedure to create the properties frames
    # Arguments:
    # w         => The base window path
    # whattab   => What tab to be activate
    global GidPriv
    
    # Create the caption frame
    ::KMProps::CreateCaption $w

    # Create the properties tabs
    ::KMProps::CreatePropertiesTabs $w $whattab

}

proc ::KMProps::CreatePropertiesTabs { w {whattab "Model"}} {
    # ABSTRACT: Procedure to create the properties tabs (model and materials)
    # Arguments:
    # w         => The base window path
    variable NbPropsPath; variable TreePropsPath
    variable TreeCurvePath
   global KPriv
    # Create the tabs
    set nb "$w.nb"
    grid [ttk::notebook $nb] -sticky ewns

    grid columnconfigure $w 0 -weight 1
    grid rowconfigure $w 1 -weight 1
    
    # Model properties
    set txt [= Model]
    set fProp ${nb}.fProp
    $nb add [ttk::frame $fProp -padding {1 1 1 1}] -text "$txt"
    set NbPropsPath $fProp

    # Set the tree properties path
    set TreePropsPath [::KMProps::CreateTreeAndToolbar $fProp]
    # wa "fProp:$fProp"
    
    # Material properties
    set txt [= Materials]
    set fMat ${nb}.fMat
    $nb add [ttk::frame $fMat -padding {1 1 1 1}] -text "$txt"
    set ::KMat::NbMatsPath $fMat
    # Set the tree material properties path
    set ::KMat::TreeMatsPath [::KMat::CreateTreeAndToolbar $fMat]
    
    if { [kipt::CurvesModule ] } {
	# Curve properties
	set txt [= Curves]
	set fCurv ${nb}.fCurv
	$nb add [ttk::frame $fCurv -padding {1 1 1 1}] -text "$txt"
	set ::KMat::NbCurvPath $fCurv
	set KPriv(TreeCurvePath) [::KCurves::CreateTreeAndToolbar $fCurv]
    }
    # Model treeß
    ::KMProps::initVisibilityClass
    ::KMProps::FillTreeProps

    # Material tree
    ::KMat::initVisibilityClass
    ::KMat::FillTreeMat
    
    if { [kipt::CurvesModule ] } {
	# Curves tree
	::KMProps::initVisibilityClass
	::KCurves::FillTreeCurves
    }
    
    # Select the active tab
    switch -exact -- $whattab {
	"Model" {
	    $nb select "$nb.fProp"
	}
	"Materials" {
	    $nb select "$nb.fMat"
	}
	"Curve" {
	if { [kipt::CurvesModule ] } {
	    $nb select "$nb.fCurv"
	}
	}
    }
}

proc ::KMProps::SwitchWindowInsideOutside {{w ".gid.kmprops"}} {
    # ABSTRACT: Switch the model window between inside or outside (1:Inside,0:Outside)
    # Arguments:
    # w => The base window path

    if { ![winfo exists $w] } {
	::KMProps::CreateWindow $w 1
    } elseif { [winfo class $w] eq "Toplevel" } {
	::KMProps::CloseWindowOutside $w
	::KMProps::CreateWindow $w 1
    } else {
	::KMProps::CloseWindowInside $w
	::KMProps::CreateWindow $w 0
    }
}

proc ::KMProps::OpenWindowOutside { w } {
    # ABSTRACT: Open the model window outside
    # Arguments:
    # w => The base window path
    # Return
    # w  => The base window path

    if { [winfo exists $w] } {
	::KMProps::CloseWindowOutside $w
    }
    
    # Init the base model window
    set title [= "Project properties"]
    InitWindow $w $title KratosModelWindowGeom ChangeKratosModelWindow "" 0

    # To close the window
    wm protocol $w WM_DELETE_WINDOW "[list ::KMProps::RefreshTree "" 1]; destroy $w"

    return $w
}

proc ::KMProps::CloseWindowOutside { w } {
    # ABSTRACT: Close the model window outside
    # Arguments:
    # w => The base window path
    
    if {[winfo exists $w]} {
	
	# Refresh the tree
	::KMProps::RefreshTree "" 1

	# Close the window
	destroy $w
    }
}

proc ::KMProps::ChangeKratosModelWindow {{w .gid.kmprops}} {
    # ABSTRACT: Change the kratos model window 
    # Arguments:
    # w => The base window path
    variable Layout

    if { ![info exists Layout] } {
	set Layout OUTSIDE
	#not use INSIDERIGTH because sharing gid.ini cause problems for previous GiD's in ReadDefaultValues
    }
    if { ($Layout eq "INSIDE_RIGHT") || ($Layout eq "INSIDE_LEFT") } {
	set inside 1
    } else {
	set inside 0
    }
    
    # problems with panedwindow: without after idle have bad size when restarting
    # gid with kratos model window inside
    ::KMProps::CreateWindow $w $inside
}

proc ::KMProps::OpenWindowInside { w } {
    # ABSTRACT: Procedure to open window inside, include creation of panedwindow to manage several windows
    # Arguments:
    # w => The base window path
    # Return
    # w => The base window path
    global GidPriv
    variable Layout; variable winpath

    if { [winfo exists $w] } {
	::KMProps::CloseWindowInside $w
    }

    if { [lsearch -exact {GEOMETRYUSE MESHUSE} [GiD_Info Project ViewMode]] == -1 } {
	return
    }
    
    ttk::frame $w

    set focus [focus]

    if { [winfo exists .gid.central.wins] } {
	set grWindow .gid.central.wins
    } else {
	set grWindow .gid.central.s
    }

    # pwCentral is a panedwindow to manage $grWindow with other windows in GiD
    if {![info exists GidPriv(pwCentral)]} {
	set GidPriv(pwCentral) .gid.central.pwCentral
    }
    
    if {![winfo exists $GidPriv(pwCentral)]} {
	set mg [winfo manager $grWindow]
	set mginfo [$mg info $grWindow]        
	$mg forget $grWindow

	panedwindow $GidPriv(pwCentral) -borderwidth 0 -showhandle 0 -sashpad 1 -sashrelief sunken -opaqueresize 0

	$mg $GidPriv(pwCentral) {*}$mginfo
	if {$Layout eq "INSIDE_RIGHT"} {
	    $GidPriv(pwCentral) add $grWindow $w -minsize 60
	} else {
	    $GidPriv(pwCentral) add $w $grWindow -minsize 60
	} 
	raise $grWindow
	if {$grWindow == ".gid.central.wins"} {
	    raise .gid.central.s
	}

	$GidPriv(pwCentral) paneconfigure $grWindow -stretch always 

	raise $w
    } else {
	set panes  [$GidPriv(pwCentral) panes]
	set maxpane [expr [llength $panes]-1]

	# Are there any window after grWindow?
	#if {[lsearch $panes $grWindow] < $maxpane} 
	if { [lsearch $panes $::KMProps::winpath] != -1} {
	    ::KMProps::CloseWindowInside $w
	    ::KMProps::CreateWindow $w 0
	} else {
	    # Save sash position
	    for {set i 0} {$i < $maxpane} {incr i} {
		set _wsash($i) [lindex [$GidPriv(pwCentral) sash coord $i] 0]
	    }
	    if {$Layout eq "INSIDE_RIGHT"} {
		set incre 0
		$GidPriv(pwCentral) add $w -after $grWindow -minsize 60
	    } else {
		$GidPriv(pwCentral) add $w -before $grWindow -minsize 60
		set incre 1
	    }
	    # Restore sash position
	    foreach i [array names _wsash]  {
		update idletask
		$GidPriv(pwCentral) sash place [expr $i+$incre] $_wsash($i) 0
	    }
	}
    }

    if { $focus ne "" } {
	update
	catch { focus -force $focus }
    }
    return $w
}

proc ::KMProps::WriteKratosModelGeomToVar { w what geomname {InitComm ""}} {
    # ABSTRACT: Procedure to write to a GidPriv global variable the geometry of the model window
    # Arguments:
    # w         => The base window path
    # geomname  => Geometry identifier
    # what      => Option ["NONE"|"OPEN"]
    # InitComm  => Window init procedure
    variable Layout
    global GidPriv
    
    set trans 1
    update idletasks
    
    if {$Layout eq "INSIDE_RIGTH"} {
	if {[info exists GidPriv($geomname)]} {
	    # maintain all dimensions except the new width
	    set prevgeom [lindex $GidPriv($geomname) 1]
	    foreach {width height x y} [split $prevgeom x+] break
	    set width [winfo width $w]
	    # wa " INSIDE_RIGTH width:${width} height:${height} x:${x} y:${y}"
	    set geom ${width}x${height}+${x}+${y}
	    set GidPriv($geomname) "$what $geom $trans $InitComm"
	}
    } elseif {$Layout eq "INSIDE_LEFT"} {
	if {[info exists GidPriv($geomname)]} {
	    # maintain all dimensions except the new width
	    set prevgeom [lindex $GidPriv($geomname) 1]
	    foreach {width height x y} [split $prevgeom x+] break
	    set width [winfo width $w]
	    # wa " INSIDE_RIGTH width:${width} height:${height} x:${x} y:${y}"
	    set geom ${width}x${height}+${x}+${y}
	    set GidPriv($geomname) "$what $geom $trans $InitComm"
	}
    } else {
	# "OUTSIDE"
	set width [winfo width $w]
	set height [winfo height $w]
	set x [winfo x $w]
	set y [winfo y $w]
	# wa " OUTSIDE width:${width} height:${height} x:${x} y:${y}"
	set geom ${width}x${height}+${x}+${y}
	set GidPriv($geomname) "$what $geom $trans $InitComm"
    }
    
}

proc ::KMProps::CreateWindow { w {inside 1} {whattab "Model"}} {
    # ABSTRACT: Procedure to create the model and material window
    #           Some changes to take account GidPriv(pwCentral), a panedwindow to manage several windows
    # Arguments:
    # w         => The base window path
    # inside    => Option for inside or outside
    # whattab   => What tab to be activate
    global GidPriv
    variable Layout
    
    if { $inside } {
	::KMProps::OpenWindowInside $w
	if { ![winfo exists $w] } return ;# windows disabled || usemorewindows == 0
	if { [winfo class $w] == "Toplevel" } return 

	# pwCentral is a panedwindow to manage $grWindow with other windows in GiD
	if {![info exists GidPriv(pwCentral)]} return
	if {![winfo exists $GidPriv(pwCentral)]} return
	
	set p [ttk::frame $w.gg]
	
	::KMProps::CreatePropertiesFrames $p $whattab

	# set size
	pack $p -expand 1 -fill both
	update idletasks
	if { [info exists ::GidPriv(KratosModelWindowGeom)] } {
	    set w1 [lindex [split [lindex $::GidPriv(KratosModelWindowGeom) 1] x] 0]
	} else {
	    set w1 [winfo reqwidth $p]
	}
	# wa "Layout:$Layout before w1:$w1 width pwCentral:[winfo width $GidPriv(pwCentral)]"
	if {$Layout eq "INSIDE_RIGHT"} {
	    set sepw 4
	    set w0 [expr [winfo width $GidPriv(pwCentral)]-$w1-$sepw]
	} elseif {$Layout eq "INSIDE_LEFT"} {
	    set w0 $w1
	} else {
	    set w0 $w1
	}
	if { $w0<=0 || $w1<=0 } {
	    set w0 750
	    set w1 250
	}
	# wa "w0:$w0 w1:$w1"
	set wsash [expr [lsearch [$GidPriv(pwCentral) panes] $w]-1]
	if {$wsash < 0} {
	    set wsash 0
	}
	
	# $w add $p ;#raise reshape_cb of togl with unexpected width!!
	if { [winfo class $GidPriv(pwCentral)] == "Panedwindow" } {
	    # problems with panedwindow: without update idletasks have bad size when
	    # reopening again the kratos model window inside
	    update idletasks
	    $GidPriv(pwCentral) sash place $wsash $w0 0 ;#raise reshape_cb of togl
	} elseif { [winfo class $GidPriv(pwCentral)] == "TPanedwindow" } {
	    $GidPriv(pwCentral) sashpos $wsash $w0
	}
	variable StartLayout
	set Layout "$StartLayout"
	bind $p <Configure> [list ::KMProps::WriteKratosModelGeomToVar $p OPEN KratosModelWindowGeom ChangeKratosModelWindow]

    } else {
	
	# Open window outside
	::KMProps::OpenWindowOutside $w

	if { ![winfo exists $w] } return ;# windows disabled || usemorewindows == 0
	set Layout "OUTSIDE"
	
	::KMProps::CreatePropertiesFrames $w $whattab
	
	# Add lower buttons
	ttk::frame $w.buts -style BottomFrame.TFrame
	ttk::button $w.buts.can -text [= "Close"] -command [list destroy $w]
	grid $w.buts.can -sticky ew -padx 5 -pady 6
	grid $w.buts -sticky ews
	grid anchor $w.buts center
	
	bind $w <Alt-c> "destroy $w"
    }
}

proc ::KMProps::CreateCaption { w } {
    # ABSTRACT: Procedure to create the caption frame
    # Arguments:
    # w         => The base window path

    ttk::frame $w.caption -style solid.TFrame -borderwidth 1
    ttk::label $w.caption.title
    
    tooltip::tooltip $w.caption [= "Double-click toggle window or inline visualization"]   
    if { [winfo class $w] != "Toplevel" } {
	$w.caption.title configure -text [= "Project properties"]
	tk::button $w.caption.close -image [gid_themes::GetImage close17.png] -command [list ::KMProps::CloseWindowInside [winfo parent $w]] \
	    -borderwidth 0 -relief flat
	$w.caption.title configure -text [= "Double click here to tear off the window"]
    } else {
	$w.caption.title configure -text [= "Double click here to integrate the window"]
    }

    if { [winfo class $w] != "Toplevel" } {
	grid $w.caption.title $w.caption.close -sticky ew
	grid configure $w.caption.close  -sticky w
    } else {
	grid $w.caption.title -sticky ew
    }
    grid $w.caption -sticky ew
    grid columnconfigure $w.caption 0 -weight 1

    bind $w.caption.title <Double-Button-1> [list ::KMProps::SwitchWindowInsideOutside]
}

proc ::KMProps::CloseWindowInside { p } {
    # ABSTRACT: This procedure take account new panedwindow GidPriv(pwCentral) to manage several windows
    # Arguments
    # p  => Panel window path
    global GidPriv
    
    if {![winfo exists $p]} return

    # pwCentral is a panedwindow to manage $grWindow with other windows in GiD
    if {![info exists GidPriv(pwCentral)]} {
	destroy $p
	return
    }
    
    if {[winfo exists $p.gg]} {
	bind $p.gg <Configure> ""
	::KMProps::WriteKratosModelGeomToVar $p.gg NONE KratosModelWindowGeom ChangeKratosModelWindow
    }

    if { [winfo exists $GidPriv(pwCentral)] } {
	if { [winfo class $GidPriv(pwCentral)] in "Panedwindow TPanedwindow" } {
	    set panes [$GidPriv(pwCentral) panes]
	    set _pos [lsearch $panes $p]
	    if {$_pos >= 0} {
		if { [llength $panes] == 2 } {
		    set oldwin [lreplace $panes $_pos $_pos]
		    set mg [winfo manager $GidPriv(pwCentral)]
		    set mginfo [$mg info $GidPriv(pwCentral)]        
		    $GidPriv(pwCentral) forget $p
		    $mg $oldwin {*}$mginfo
		    update idletask
		    destroy $GidPriv(pwCentral)
		    unset GidPriv(pwCentral)
		}
	    }
	}
    } else {
	unset GidPriv(pwCentral)
    }

    set truco 0
    if { [llength $panes] == 3 } {
	set truco 1
	set x [lindex [$GidPriv(pwCentral) sash coord 1] 0]
    }

    destroy $p 
    update idletasks
    if { $truco == 1 } {        
	$GidPriv(pwCentral) sash place 0 $x 0        
    }

}

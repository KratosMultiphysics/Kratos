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
#     0.2-20/09/12- J. Garate, Creation of Curves Interface and managment functions
#     0.1-21/08/12-G. Socorro, create the base source code
#
###############################################################################

 namespace eval ::KTables {
	
 }

proc ::KTables::TableFrameDSetsInitData { action table } {
	set prpid "Curve"
	set cindex "$prpid,Add"
	if { $action == "crea" } {
		#Borro todos los puntos según la cantidad anterior en la tabla
		$table delete 0 end
	} else {
		set NumberPoints [expr ([llength $::KTables(pairs)] / 2) ]
		set PointLabel 1
		set j 0
		set ::KTables::props($cindex,AllProp) ""
		for { set i 0 } { $i < $NumberPoints } { incr i } {
			lappend ::KTables::props($cindex,AllProp) [list $PointLabel [lindex $::KTables(pairs) $j] [lindex $::KTables(pairs) [expr ($j+1)]]]
			set j [expr ($j+2)]
			incr PointLabel
		}
	}
}

proc ::KTables::TableFrameCreateDSets { w action } {
    package require math
    
    variable props
    global KPriv
    global SCPriv; global GidPriv
    
    
    set ::KTables::w $w
    
    set prpid "Curve"
    set cindex "$prpid,Add"
    
    # Find the selected curve at xml
    set selectedItem $::KCurves::mycurve
    if {$selectedItem == ""} {
        return ""
    }
    set fullname [::KCurves::GetFullname $selectedItem]
    set pid [::KCurves::setXml $fullname pid ""]
    
    # Create the Big Frame
    set BigFrame [ttk::labelframe $w.bigframe -text [= "Curve Properties - $pid"]]
    set ::KTables::BigFrame $BigFrame
    
    # Find the proper XVal and YVal
    set xpath [::xmlutils::setXPath $fullname ]
    set nodes [$::KPriv(xml) selectNodes $xpath]
    set id [$nodes getAttribute id ""]
    set pid [$nodes getAttribute pid ""]
    set xnod [$nodes find id "XVar"]
    set xCoord [$xnod getAttribute dv ""]
    set xival [$xnod getAttribute ivalues ""]
    set xival [::KMProps::split2 $xival "," ]
    set ynod [$nodes find id "YVar"]
    set yCoord [$ynod getAttribute dv ""]
    set yival [$ynod getAttribute ivalues ""]
    set yival [::KMProps::split2 $yival "," ]
    
    set combowidth [::math::max [::KMProps::getCmbWidth $xival] [::KMProps::getCmbWidth $yival]]

    set Comboframe [ttk::frame $BigFrame.comboframe] 
    
    
    # Combo labels
    set Comboframe.lblx [ttk::label $Comboframe.lblx -text [= "XCoord"]]
    set Comboframe.lbly [ttk::label $Comboframe.lbly -text [= "YCoord"]]
    
    # Comboboxes
    set Comboframe.cmbx [ttk::combobox $Comboframe.cmbx -values $xival -textvariable $xCoord -width $combowidth]
    set Comboframe.cmby [ttk::combobox $Comboframe.cmby -values $yival -textvariable $yCoord -width $combowidth]
    
    grid $Comboframe.cmbx -row 0 -column 1 -sticky nw -in $Comboframe
    grid $Comboframe.cmby -row 1 -column 1 -sticky nw -in $Comboframe
    grid $Comboframe.lblx -row 0 -column 0 -sticky nw -in $Comboframe -pady 2 -sticky nw
    grid $Comboframe.lbly -row 1 -column 0 -sticky nw -in $Comboframe -pady 2 -sticky nw
    
    grid $Comboframe -row 0 -columnspan 1 -sticky nsew -in $BigFrame
    
    $Comboframe.cmbx set $xCoord 
    $Comboframe.cmby set $yCoord 
    
    # Frame's Grid
    set Tableframe [ttk::labelframe $BigFrame.table -text [= "Table"]] 
    set Graphframe [ttk::labelframe $BigFrame.graph -text [= "Graph"]]
    
    grid $Tableframe -row 2 -in $BigFrame
    grid $Graphframe -row 3 -in $BigFrame -sticky we
    
    set tam 250
    set windsize [winfo screenmmheight $BigFrame]
	if { $windsize < 250 } {
		set tam 150
	}
    
    # Set the Canvas for the Graph
    set cpath [canvas $Graphframe.c \
	-relief sunken -borderwidth 2 -width 250 -height $tam]

    
    set ::KTables::cpath $cpath

    grid $cpath -in $Graphframe \
        -row 0 -column 0 -sticky nws

    set tffLTable $Tableframe
    
    # Register some widgets from the BWidget package for interactive cell editing
    package require tablelist_tile
    ::tablelist::addBWidgetEntry
    ::tablelist::addBWidgetComboBox
    
    # Create a table list
    set vsb $Tableframe.vsb1
    set hsb $Tableframe.hsb1
    
    set idname [= "Point Id"]

    set clist [list 0 "$idname" center \
		       0 "$xCoord" center \
		       0 "$yCoord" center]
 
    # Create the table 
    set tbl [::tablelist::tablelist $Tableframe.tbl \
  		 -exportselection 0 \
  		 -background $GidPriv(Color,BackgroundListbox) \
  		 -columns $clist \
  		 -stretch all \
		 -listvariable ::KTables::props($cindex,AllProp) \
  		 -selectmode extended \
  		 -selecttype row \
		 -selectbackground #447bcd \
  		 -selectforeground white \
  		 -labelborderwidth 1 \
  		 -labelpady 0 \
  		 -showarrow 0 \
  		 -showseparators yes \
		 -xscrollcommand [list $hsb set] \
		 -yscrollcommand [list $vsb set] \
		 -resizablecolumns yes \
		 -editstartcommand "::KTables::TableFrameDSetsTableSetValues" -height 8 -width 40\
		 -editendcommand "::KTables::TableFrameDSetsTableApplyValues" ]
  
    ttk::scrollbar $vsb -orient vertical   -command [list $tbl yview]
    ttk::scrollbar $hsb -orient horizontal -command [list $tbl xview]
    
    # Binding the upper comboboxes to TableList
    bind $Comboframe.cmbx  <<ComboboxSelected>> "::KTables::cmbSelectChange $Comboframe.cmbx $tbl $fullname"
    bind $Comboframe.cmby  <<ComboboxSelected>> "::KTables::cmbSelectChange $Comboframe.cmby $tbl $fullname"
    
    # Init Table set add Table frame data
    ::KTables::TableFrameDSetsInitData $action $tbl
    
    # Configure column edition
    $tbl columnconfigure 0 -resizable 1 -maxwidth 20 -editable no -editwindow Entry
    $tbl columnconfigure 1 -resizable 1 -maxwidth 10 -editable yes -editwindow Entry
    $tbl columnconfigure 2 -resizable 0 -maxwidth 10 -editable yes -editwindow Entry
  
    # Bottom Buttons
    ButtonBox $Tableframe.buttons \
	    -spacing 0 \
	    -padx 1 

    $Tableframe.buttons add \
	-image [image create photo -file [file join $KPriv(dir) $KPriv(imagesdir) tableArrowDown.png]] \
	-highlightthickness 0 \
	-takefocus 0 -relief link \
	-borderwidth 1 -padx 1 -pady 1 \
	-helptext [= "Add new point to table"]. \
	-command "::KTables::TableFrameDSetsbNewLineEvent $tbl"
    
    $Tableframe.buttons add \
	-image [image create photo -file [file join $::KPriv(dir) $KPriv(imagesdir) tableArrowUp.png]] \
	-highlightthickness 0 \
	-takefocus 0 -relief link \
	-borderwidth 1 -padx 1 -pady 1 \
	-helptext [= "Delete the last point set from the table"]. \
	-command "::KTables::TableFrameDSetsbDeleteLastLineEvent $tbl"

    $Tableframe.buttons add \
	-image [image create photo -file [file join $::KPriv(dir) $KPriv(imagesdir) tableDeleteSel.png]] \
	-highlightthickness 0 \
	-takefocus 0 -relief link \
	-borderwidth 1 -padx 1 \
	-helptext [= "Delete the selected points"]. \
	-command "::KTables::TableFrameDSetsbDeleteEvent $tbl"
    
    $Tableframe.buttons add \
	-image [image create photo -file [file join $::KPriv(dir) $KPriv(imagesdir) tableDeleteAll.png]] \
	-highlightthickness 0 \
	-takefocus 0 -relief link \
	-borderwidth 1 -padx 1 \
	-helptext [= "Delete all points in the table"]. \
	-command "::KTables::TableFrameDSetsbDeleteAllEvent $tbl"
	
 	$Tableframe.buttons add \
  	-image [image create photo -file [file join $::KPriv(dir) $KPriv(imagesdir) PlotGraph.png]] \
  	-highlightthickness 0 \
  	-takefocus 0 -relief link \
  	-borderwidth 1 -padx 1  \
  	-helptext [= "Show/Hide Graph"]. \
  	-command "::KTables::ShowHide_Graph $Graphframe $Tableframe"

    
    # Geometry management
    
    # For Table sets
    grid $BigFrame \
	-in $w \
	-row 1 -column 0 -columnspan 1 \
	-sticky nsw 
    
    grid $tbl $vsb -sticky nsew -row 0
    
    grid $Tableframe \
	-in $BigFrame \
	-row 2 -column 0 \
	-sticky nsw -columnspan 1
    
    grid $Tableframe.buttons -row 1 -sticky w
    
    # Expand
    # For add poimts sets
    grid rowconfigure $w 1 -weight 1
    grid rowconfigure $w 2 -weight 1
    grid columnconfigure $w 0 -weight 1
    grid columnconfigure $w 1 -weight 1
    grid columnconfigure $w 2 -weight 1
  
    # For points set table
    grid rowconfigure $Tableframe 0 -weight 1
    grid columnconfigure $Tableframe 0 -weight 1
  
    grid configure $vsb -sticky ns
   
    set bodyTag [$tbl bodypath]
    catch {
        bind $tbl <<ListboxSelect>> 
        bind $bodyTag <<Button3>>
    }
    #bind $tbl <<ListboxSelect>> [list ::SCUtils::ConfigureListScrollbars $tbl $hsb $vsb]
    # Add some bind for botton-3
    # bind $bodyTag <<Button3>>  [bind TablelistBody <Button-1>]
    # bind $bodyTag <<Button3>> +[bind TablelistBody <ButtonRelease-1>]
    # bind $bodyTag <<Button3>> +[list ::KTables::TableFrameAddShowContextMenu $tbl %X %Y $fpath]
    bind $tbl <<TablelistCellUpdated>> "::KTables::UpdateTable"
   
    # Update frame
    ::KTables::TableFrameUpdateDSets $w
    
    #Fill the table with points saved at xml
    ::KTables::FillExistingPoints $tbl
}

proc ::KTables::ShowHide_Graph { Graphframe Tableframe } {
    # When clicked the graph button, show the graph if hidden, else hide it.
    
    if {[winfo exists $Graphframe]} {
        destroy $Graphframe
        grid configure $Tableframe -sticky nsew
    } else {
        set Graphframe [ttk::labelframe $::KTables::BigFrame.graph -text [= "Graph"]]
        grid $Graphframe -row 3 -in $::KTables::BigFrame -sticky we
        
    set tam 250
    set windsize [winfo screenmmheight $::KTables::BigFrame]
	if { $windsize < 250 } {
		set tam 150
	}
        # Set the Canvas for the Graph
        set cpath [canvas $Graphframe.c -relief sunken -borderwidth 2 -width 250 -height $tam]

        grid $cpath -in $Graphframe -row 0 -column 0 -sticky nws
        ::KTables::CreaByPoints_Graph
    }
}

proc ::KTables::CreaByPoints_Graph { } {
    
    set aux $::KTables::BigFrame
    append aux ".comboframe.cmb"
    set title "Curve definition" 
    set xcombo [append aux "x"]
    set aux [string range $aux 0 end-1]
    set ycombo [append aux "y"]

    set xlabel [$xcombo get]
    set ylabel [$ycombo get] 
    set pvalues ""
    set cindex "Curve,Add"
    
    set NumberPoints [llength $::KTables::props($cindex,AllProp)]
    for {set i 0} { $i < $NumberPoints } { incr i } {
        set x [lindex $::KTables::props($cindex,AllProp) $i 1]
        set y [lindex $::KTables::props($cindex,AllProp) $i 2]
        lappend pvalues [list $x $y]
    }
    
    ::KPlot::Plot $::KTables::cpath $pvalues $title $xlabel $ylabel
}

proc ::KTables::FillExistingPoints { table } {

    global KPriv
    
    # Find at xml the Curve node, to add the points to xml
    set T $KPriv(TreeCurvePath)
	set selected [$T selection get]
    set fullname [::KCurves::GetFullname $selected]
    set xpath [::KCurves::setXPath $fullname]
    set node [$KPriv(xml) selectNodes $xpath]
    # $node contains the proper Curve

    set nodeCont [$node find id "Points"] 
    # $nodeCont is the Point Container
    
    # Get the points
    set points [::KTables::GetPointsFromContNodeXML $nodeCont]
    foreach {id xp yp} $points {
        # Add each point to the table
        set list "$id $xp $yp"
        $table insert end $list
        
    }
    if {[llength $points] >= 4} {
        ::KTables::CreaByPoints_Graph
    }
}

proc ::KTables::GetPointsFromContNodeXML { node } {

    set list ""
    set childs [$node childNodes]
    foreach child $childs {
        set id [$child getAttribute "id" "-1"]
        set x [$child getAttribute "Xval" "-1"]
        set y [$child getAttribute "Yval" "-1"]
        set list [lappend list $id $x $y]
    }
    return $list
}

proc ::KTables::UpdateTable { } {
    # Event called when a Table Cell is updated
    ::KTables::CreaByPoints_Graph
}

proc ::KTables::cmbSelectChange { combo table fullname} {
    global KPriv
    
    set combobase [string range $combo 0 end-1]
    set xcombo [append combobase "x"]
    set combobase [string range $combo 0 end-1]
    set ycombo [append combobase "y"]

    set newvaluex [$xcombo get]
    set newvaluey [$ycombo get]

    set xpath [::xmlutils::setXPath $fullname ]
    set nodes [$::KPriv(xml) selectNodes $xpath]
    set xnod [$nodes find id "XVar"]
    set ynod [$nodes find id "YVar"]

    # xml update
    $xnod setAttribute dv $newvaluex
    $ynod setAttribute dv $newvaluey
    
    # TableList Update
    
    set idname [= "Point Id"]
    set clist [list "$idname" "$newvaluex" "$newvaluey"]
    
    $table configure -columntitles $clist
    ::KCurves::RefreshTree $KPriv(TreeCurvePath)
    ::KTables::CreaByPoints_Graph
}

proc ::KTables::TableFrameAddShowContextMenu {w rootx rooty win} {
    # rootx,rooty are screen coordinates (for knowing where to place the menu)
    
    catch { destroy $w.contextMenu}
    menu $w.contextMenu -tearoff false

    $w.contextMenu add command \
	-label [= "Delete the selected pair set"] \
	-command "::KTables::TableFrameDSetsbDeleteEvent $win"

    $w.contextMenu add command \
	-label [= "Delete all pairs set"] \
	-command "::KTables::TableFrameDSetsbDeleteAllEvent $win"

    tk_popup $w.contextMenu $rootx $rooty
}

proc ::KTables::TableFrameUpdateDSets {fpath} {
     
    if {[winfo exists $fpath]} {
	set tffLTable $fpath.bigframe
	set tbl $tffLTable.tbl
	#set cindex "Curve,Add"
	#WarnWinText que:$::KTables::props($cindex,AllProp) 
	# Init Table set add Table frame data
	#::KTables::TableFrameDSetsInitData $stgid
    }
}

proc ::KTables::TableFrameDSetsChangeWidgetState {fpath what} {
    # Enabled/disabled add Table set widgets
    # Arguments
    # fpath => Frame path
    # what  => Option to enable or disabled widgets

    if {[winfo exists $fpath]} {
        set tffLTable $fpath.bigframe
        set wdlist [list $fpath.bigframe $tffLTable.tbl $fpath.bbox1]
        foreach wd $wdlist {
            $wd configure -state $what
        }
    }
}

proc ::KTables::TableFrameDSetsTableSetValues {tbl row col text} {
    # Applies some configuration options to the edit window; if the latter is a
    # ComboBox, the procedure populates it
    
    # WarnWinText "start=> tbl:$tbl row:$row col:$col text:$text"
    set w [$tbl editwinpath]
    # Fill some combo boxes
    switch $col {
	"0" {
	}
	"1" {
	}
	"2" {
	}
    }
    return $text
}

proc ::KTables::TableFrameDSetsTableApplyValues {tbl row col text} {

    # Performs a final validation of the text contained in the edit window and gets
    # the cell's internal contents
    variable winpath

    # WarnWinText "tbl:$tbl row:$row col:$col text:$text"
    set w [$tbl attrib widget]
    # Get old value
    set rvalue [$tbl cellcget $row,$col -text]
    # Check values
    switch $col {
        0 {
        }
        1 {
            # Frequency value
            set cvalue [$tbl cellcget $row,$col -text]
            # WarnWinText "cvalue:$cvalue"
            if {[regexp {^[-+]?[0-9]*\.?[0-9]*([0-9]\.?[eE][-+]?[0-9]*)?$} $text] ==1} {
            set rvalue $text
            } else {
            set rvalue $cvalue
            WarnWin [= "Define a real value for the frequency"].
            return $rvalue
            }
        } 
        2 {
            # Time value
            set cvalue [$tbl cellcget $row,$col -text]
            # WarnWinText "cvalue:$cvalue"
            if {[regexp {^[-+]?[0-9]*\.?[0-9]*([0-9]\.?[eE][-+]?[0-9]*)?$} $text] ==1} {
            set rvalue $text
            } else {
            set rvalue $cvalue
            WarnWin [= "Define a real value for the time"].
            return $rvalue
            }
        }
    }
    
    return $rvalue
}

proc ::KTables::TableFrameDSetsbDeleteLastLineEvent {tbl} {
       
    if {[winfo exists $tbl]} {

        if {[llength [$tbl get 0 end]]>"0"} { 
            # Clear the last row from the list
            $tbl delete end end
        } else {
            WarnWin [= "All point sets for this stage have been removed"].
            return ""
        }
    }
}

proc ::KTables::AcceptButton { } {
    
    global KPriv
    set prpid "Curve"
    set cindex "$prpid,Add"
    set NumberPoints [llength $::KTables::props($cindex,AllProp)]
    if { $NumberPoints < 1}  {
        return ""
    }
    
    set points $::KTables::props($cindex,AllProp)
    
    # Find at xml the Curve node, to add the points to xml
	set selected $::KCurves::mycurve
    set fullname [::KCurves::GetFullname $selected]
    set xpath [::KCurves::setXPath $fullname]
    set node [$KPriv(xml) selectNodes $xpath]
    # $node contains the proper Curve
    
    # First, update the xcoord and ycoord value, and the curve type
    set aux $::KTables::BigFrame
    append aux ".comboframe.cmb"
    set title "Curve definition" 
    set xcombo [append aux "x"]
    set aux [string range $aux 0 end-1]
    set ycombo [append aux "y"]
    set xval [$xcombo get]
    set yval [$ycombo get] 
    
    set xitem [expr $selected + 2]
    set fullname [::KCurves::GetFullname $xitem]
    ::xmlutils::setXml $fullname "dv" "edit" $xval
    
    set yitem [expr $selected + 3]
    set fullname [::KCurves::GetFullname $yitem]
    ::xmlutils::setXml $fullname "dv" "edit" $yval
    
    # Find the xml container for the points and erase old points information
    set childs [[$node find id "Points"] childNodes]
    foreach child $childs {
        $child delete
    }
    
    set nodeCont [[$node find id "Points"] asList]
    
    # Get the number of current points at the xml
    set np [lindex [lindex $nodeCont 1] 5]
    
    foreach point $points {
        incr np 1
        set id [lindex $point 0]
        set xp [lindex $point 1]
        set yp [lindex $point 2]
        # Adding each new point to a list format node
        set nodeCont [::KTables::AddPointToXml $nodeCont $id $xp $yp]
    }

    # Delete the original node
    [$node find id "Points"] delete
    
    # Adding the new node and updating the number of points
    set nodeaux [lreplace [lindex $nodeCont 1] 5 5 $np]
    set nodeCont [lreplace $nodeCont 1 1 $nodeaux]
    
    # Then, add the points to xml
    $node appendFromList $nodeCont
    
    # Finally, close the bottom frame
    ::KTables::CancelButton
}

proc ::KTables::CancelButton { } {

    ::KCurves::DestroyBottomFrame
    ::KCurves::FillTreeCurves
    ::KCurves::DestroyButtons
    ::KCurves::RefreshTree 
}

proc ::KTables::AddPointToXml { node id xp yp } {

    set node [string range $node 0 end-1]
    set titem " \{TItem \{id "
    set titem [append titem $id]
    set titem [append titem " Xval "]
    set titem [append titem $xp]
    set titem [append titem " Yval "]
    set titem [append titem $yp]
    set titem [append titem " \} \{\}\}\}"]
    #set titem [append titem " \{\}\}\}\}"]
    set node [append node $titem ]

    return $node
}

proc ::KTables::TableFrameDSetsbNewLineEvent {tablepath} {
  
    if {[winfo exists $tablepath]} {
        
        set tbl $tablepath

        # Get a new Table set identifier
        set TotalPoints [llength [$tbl get 0 end]]
        set NewDSetId [expr ($TotalPoints + 1)]
        set newlist [list $NewDSetId 0 0]
        if {[llength $newlist]>0} {
            $tbl insert end $newlist
        }
    }
}

proc ::KTables::TableFrameDSetsbDeleteAllEvent {tbl} {
  
    if {[winfo exists $tbl]} {


	# All element
	set AssigList [$tbl get 0 end]
	# WarnWinText "AssigList:$AssigList"
	if {[llength $AssigList]>"0"} {  
	    set txt1 [= "Do you want to remove all defined points ?"]
	    set w1 .gid.tmpwin
	    set ret [tk_dialogRAM $w1 [= "Warning"] "$txt1" warning 0 [= "Yes"] [= "No"]]
	    if { $ret == 0 } {
            # Clear all the list
            $tbl delete 0 end
	    }
	} else {
	    WarnWin [= "All points have been removed"].
	    return ""
	}
    }
}

proc ::KTables::TableFrameDSetsbDeleteEvent {tbl} {
 
    if {[winfo exists $tbl]} {


        # All select element
        set items [$tbl curselection]
        # WarnWinText "items:$items"
        if {[string length $items] == 0} {
            WarnWin [= "No point have been selected for deletion"].
            return ""
        }

        set firstlist [$tbl get 0 end]
        set count -1
        foreach item $items {
            incr count 1
            if {$count==0} {
            set newlist [lreplace $firstlist $item $item]
            } else {
            set newlist [lreplace $newlist [expr $item-$count] [expr $item-$count]]
            }
        }

        # Delete all element
        $tbl delete 0 end
        # Add all elements
        foreach row $newlist {
            $tbl insert end $row
        }
    }
}

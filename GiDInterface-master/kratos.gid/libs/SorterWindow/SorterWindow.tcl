namespace eval SorterWindow {
    variable winpath
    variable plot
    variable PosX
    variable PosY
    variable Delta
    variable Index
    variable Listado
    variable Scrollposition
    variable Scrollarea
    variable imgs
    variable window_state
    variable data_source
    variable data_source_list
    variable update_proc
    variable data
}

proc SorterWindow::Init {} {
    variable imgs
    variable window_state
    variable winpath
    set winpath ".gid.sortwindow"
    #variable dir
    set dir [file dirname [info script]]
    set imgs(lock) [image create photo -file [file join $dir images lock.png]]
    set imgs(unlock) [image create photo -file [file join $dir images unlock.png]]
    set window_state "unlocked"
    
    variable data_source
    set data_source ""
    variable data_source_list
    set data_source_list [list ]
}

#
# Open a widow for the test dialog. Left half = Scrollbox.
# Right: Some entries and a message for testing the code.
#
proc SorterWindow::Window { } {
    variable PosX
    variable PosY
    variable Delta
    variable winpath
    set w $winpath
    
    if {winfo exists $w} {destroy $w}
    
    ###################
    # CREATING WIDGETS
    ###################
    toplevel $w -class Toplevel -relief groove 
    #wm maxsize $w 500 300
    wm minsize $w 500 300
    wm overrideredirect $w 0
    wm resizable $w 1 1
    wm deiconify $w
    wm title $w [= "Sort conditions window"]
    set Delta 30
    wm attribute $w -topmost 1

    
    SorterWindow::RefreshWindow
    
}

proc SorterWindow::RefreshWindow { } {
    variable winpath
    set w $winpath
    if {[winfo exists $w.fr1]} {destroy $w.fr1}
    if {[winfo exists $w.fr2]} {destroy $w.fr2}
    if {[winfo exists $w.buts]} {destroy $w.buts}
    set fr1 [ttk::frame $w.fr1]
    set fr2 [ttk::labelframe $w.fr2 -text [= "Sort the conditions:"]]
    set buts [ttk::frame $w.buts -style BottomFrame.TFrame]
    
    # One call creates the listbox
    #
    set Canv [SortByDragListbox $fr1 10 10 400 240]
    #
    # Widgets on the right side
    #
    ttk::button $buts.q -text Cancel -command [list destroy $w] -style BottomFrame.TButton
    ttk::button $buts.ok -text Ok -command [list SorterWindow::ReturnData] -style BottomFrame.TButton
    SorterWindow::ConfigurationFrame $fr2
    
    
    grid $buts.ok $buts.q -sticky sew
    
    grid $fr1 -sticky nsew -row 0 -column 0
    grid $fr2 -sticky nw -row 0 -column 1 -padx 20
    grid $buts -sticky sew -columnspan 2
    grid columnconfigure $w 1 -weight 1 
    grid rowconfigure $w 0 -weight 1
    if { $::tcl_version >= 8.5 } { grid anchor $buts center }
}

proc SorterWindow::ConfigurationFrame {w} {
    variable imgs
    variable window_state
    variable data_source_list
    
    set locktext [= "Drag&Drop unlocked"]:
    set unlocktext [= "Drag&Drop locked"]:
    set lt $locktext
    if {$window_state eq "locked"} {set lt $unlocktext}
    set lab1 [ttk::label $w.l1 -text $lt]
    
    set li $imgs(unlock)
    if {$window_state eq "locked"} {set li $imgs(lock)}
    set b1 [ttk::button $w.b1 -image $li -command [list SorterWindow::LockButtonClicked]]
    grid $lab1 $b1 -sticky w
    
    set lab2 [ttk::label $w.l2 -text [= "Data source"]:]
    set cb2 [ttk::combobox $w.cb2 -textvariable SorterWindow::data_source -values $data_source_list -width 10 -state readonly]
    bind $cb2 <<ComboboxSelected>> SorterWindow::DataSourceChanged
    grid $lab2 $cb2 -sticky w
}

proc SorterWindow::LockButtonClicked { } {
    variable window_state
    if {$window_state eq "locked"} {set window_state "unlocked"} {set window_state "locked"}
    SorterWindow::RefreshWindow
}
proc SorterWindow::DataSourceChanged { } {
    SorterWindow::RefreshWindow
}

#
# Create a pseudo-listbox with canvas elements. Looks like a listbox,
# but is really a canvas, and all widgets only pretend to be what they seem.
#
proc SorterWindow::SortByDragListbox { w XNull YNull width height } {
    variable Index
    variable Listado
    variable Scrollposition
    variable Scrollarea
    variable window_state
    if {winfo exists $w.cv} {destroy $w.cv}
    set Canv [canvas $w.cv -borderwidth 0 -highlightthickness 0  -height [expr $height + 2*$YNull] -width [expr $width + 2*$XNull] ]

    $Canv create rectangle $XNull $YNull [expr $XNull + $width] [expr $YNull + $height] -outline black -width 1 -fill white -tags Box
    
    if {$window_state eq "locked"} {
        $Canv create rectangle [expr $XNull - 1] [expr $YNull - 1] [expr $XNull + $width + 1] [expr $YNull + $height + 1] -outline red -width 1 -tags Box -fill seashell
        $Canv configure -state disabled
    } else {
        $Canv create rectangle [expr $XNull - 1] [expr $YNull - 1] [expr $XNull + $width + 1] [expr $YNull + $height + 1] -outline grey50 -width 1 -tags Box
    }
    
    if {winfo exists $w.lbscroll} {destroy $w.lbscroll}
    ttk::scrollbar $w.lbscroll -command "SorterWindow::SortByDragListboxScroll $w" -orient vert
    grid $Canv -row 0 -column 0 -sticky nsew
    grid $w.lbscroll -row 0 -column 1 -sticky wns
    
    SorterWindow::FillListado
    set Scrollposition 0
    set Schrifthoehe 16
    # Cambiar el 20 por el numero de tal
    set Scrollarea [expr $Schrifthoehe * 20]
    SorterWindow::SortByDragListboxScroll $w scroll 0.0 units
    
    return $Canv
}

#
# Scrollbar code.
#
proc SorterWindow::SortByDragListboxScroll { w {was moveto} {Zahl 0.0} {Einheit units} } {
    variable Index
    variable Listado
    variable Scrollposition
    variable Scrollarea
    
    set Canv  $w.cv
    set height [lindex [$Canv configure -height] 4]
    set Schrifthoehe 16
    set Scrollposition 0
    
    if {$was == "scroll"} {
        if {$Einheit == "pages"} {
            incr Scrollposition [expr int($Zahl * $height - 20)]
        } else {
            incr Scrollposition [expr 20 * int($Zahl)]
        }
    } else {
        set Scrollposition [expr int($Zahl * $Scrollarea)]
    }
    
    # Limit the scrollposition to sensible values
    #
    if {$Scrollposition > [expr $Scrollarea - $height]} {
        set Scrollposition [expr $Scrollarea - $height]
    }
    if {$Scrollposition < 0} {set Scrollposition 0}
    #
    # Delete Index and built anew from scratch. In priciple all entries could
    # be moved, but this is messy at the edges.
    #
    set yPos [expr 32 - $Scrollposition]
    for {set i 0} {$i < [array size Listado]} {incr i} {
        $Canv delete ent$i
        if {$yPos < 1} {
            incr yPos $Schrifthoehe
            continue
        }
        #
        if {$yPos < [expr $height - 1]} {
            $Canv create text 24 $yPos -text $Listado($Index($i)) -anchor w -fill black -tags ent$i
            incr yPos $Schrifthoehe
            
            $Canv bind ent$i <1>     "SorterWindow::plotDown $Canv %x %y"
            $Canv bind ent$i <B1-Motion>    "SorterWindow::plotMove $Canv %x %y"
            $Canv bind ent$i <ButtonRelease-1> "SorterWindow::plotCopy $w $Canv %x %y $i"
        }
    }
    #
    $w.lbscroll set [expr double($Scrollposition) / $Scrollarea] [expr double($height + $Scrollposition) / $Scrollarea]
}

#
# plotDown --
# This procedure is invoked when the mouse is pressed over one of the
# data points. It sets up state to allow the point to be dragged.
#
# Arguments:
# w -       The canvas window.
# x, y -    The coordinates of the mouse press.
#
proc SorterWindow::plotDown {w x y} {
    variable plot
    #
    $w dtag selected
    $w addtag selected withtag current
    $w raise current
    set plot(lastX) $x
    set plot(lastY) $y
}

# plotMove --
# This procedure is invoked during mouse motion events. It drags the
# current item.
#
# Arguments:
# w -       The canvas window.
# x, y -    The coordinates of the mouse.
#
proc SorterWindow::plotMove { w x y } {
    variable plot
    variable PosX
    variable PosY
    
    $w move selected [expr $x-$plot(lastX)] [expr $y-$plot(lastY)]
    set plot(lastX) $x
    set plot(lastY) $y
    set PosX        $x
    set PosY        $y
}

#
# When the mouse button is released, this routine determines the new
# position and re-orders the list.
#
proc SorterWindow::plotCopy { w Cv x y i } {
    variable Delta
    variable Index
    variable Scrollposition
    variable Scrollarea
    
    set Schrifthoehe 16
    
    W $i
    set max [expr [array size Index] -1]
    W "max $max"
    set min 0
    set Rang [expr int(($y - $Delta + $Scrollposition) / $Schrifthoehe)]
    set Temp $Index($i)
    if {$Rang > $i} {
        if {$Rang > $max} {set Rang $max}
        for {set j $i} {$j < $Rang && $j} {incr j} {
            set k [expr $j + 1]
            W "AP j $j k $k"
            set Index($j) $Index($k)
        }
    } elseif {$Rang == $i} {
        set number [expr double($Scrollposition) / $Scrollarea]
        SorterWindow::SortByDragListboxScroll $w scroll $number units
        return
    } else {
        set Rang [expr $Rang + 1]
        if {$Rang < $min} {set Rang $min}
        W "Rang $Rang"
        for {set j $i} {$j > $Rang} {incr j -1} {
            set k [expr $j - 1]
            W "DA j $j k $k"

            set Index($j) $Index($k)
        }
    }
    set Index($Rang) $Temp
    #
    # Now scroll the list to the right position.
    #
    set number [expr double($Scrollposition) / $Scrollarea]
    SorterWindow::SortByDragListboxScroll $w scroll $number units
}

proc SorterWindow::FillListado { } {
    variable Index
    variable Listado
    variable data_source
    variable data
    unset -nocomplain Index 
    set Items [list ]
    
    if {$data_source ni [dict keys $data]} {set data_source [lindex [dict keys $data] 0]}
    
    foreach item [dict get $data $data_source] {
        lappend Items "$item"
    }
    set size [llength $Items]
    set i 0
    foreach {item groups} $Items {
        W "$item $groups"
        foreach {group num} $groups {
            W "- $item $group"
            set Listado($i) "$item $group"
            lappend Index($i) $i
            incr i
        }
    }
}

proc SorterWindow::ReturnData { } {
    variable update_proc
    variable data
    variable winpath
    
    $update_proc $data
    destroy $winpath
}

proc SorterWindow::SorterWindow {{datadict ""} {up_pr default}} {
    variable update_proc
    set update_proc $up_pr
    variable data
    if {$datadict eq ""} {set datadict [dict create]}
    set data $datadict
    
    SorterWindow::Window 
    
}

SorterWindow::Init
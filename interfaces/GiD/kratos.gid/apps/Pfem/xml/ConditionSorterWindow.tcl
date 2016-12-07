namespace eval Pfem::xml::CndSortWindow {
        variable plot
        variable PosX
        variable PosY
        variable Delta
        variable Pref
        variable Index
        variable Eintrag
        variable Scrollposition
        variable Scrollbereich
}
    #
    # A good GUI needs only one mouse button
    #
    event add <<Loslassen>> <ButtonRelease-1>
    event add <<Loslassen>> <ButtonRelease-2>
    event add <<Loslassen>> <ButtonRelease-3>
    event add <<Ziehen>>    <B1-Motion>
    event add <<Ziehen>>    <B2-Motion>
    event add <<Ziehen>>    <B3-Motion>
    event add <<Klick>>     <1>
    event add <<Klick>>     <2>
    event add <<Klick>>     <3>

    #
    # Stuff from Visual Tcl. Not pretty, and I don't know if I really need
    # all this, but it works.
    #
    proc Pfem::xml::CndSortWindow::Window {args} {
    #
        set cmd [lindex $args 0]
        set name [lindex $args 1]
        set newname [lindex $args 2]
        set rest [lrange $args 3 end]
        if {$name == "" || $cmd == ""} {return}
        if {$newname == ""} {
            set newname $name
        }
        set exists [winfo exists $newname]
        switch $cmd {
            show {
                if {$exists == "1" && $name != "."} {wm deiconify $name; return}
                if {[info procs vTclWindow(pre)$name] != ""} {
                    eval "vTclWindow(pre)$name $newname $rest"
                }
                if {[info procs vTclWindow$name] != ""} {
                    eval "vTclWindow$name $newname $rest"
                }
                if {[info procs vTclWindow(post)$name] != ""} {
                    eval "vTclWindow(post)$name $newname $rest"
                }
            }
            hide    { if $exists {wm withdraw $newname; return} }
            iconify { if $exists {wm iconify $newname; return} }
            destroy { if $exists {destroy $newname; return} }
        }
    }

    #################################
    # VTCL GENERATED GUI PROCEDURES
    #

    proc Pfem::xml::CndSortWindow::vTclWindow. {base} {
        if {$base == ""} {
            set base .
        }
    ###################
    # CREATING WIDGETS
    ###################
        wm focusmodel $base passive
        wm geometry $base 1x1+25+65
        wm maxsize $base 817 594
        wm minsize $base 1 1
        wm overrideredirect $base 0
        wm resizable $base 1 1
        wm withdraw $base
        wm title $base "Wish"
    }

    #
    # Open a widow for the test dialog. Left half = Scrollbox.
    # Right: Some entries and a message for testing the code.
    #
    proc Pfem::xml::CndSortWindow::vTclWindow.dialog {base} {
        variable Pref
        variable PosX
        variable PosY
        variable Delta
    #
        if {$base == ""} {
            set base .dialog
        }
        if {[winfo exists $base]} {
            wm deiconify $base; return
        }
    ###################
    # CREATING WIDGETS
    ###################
        toplevel $base -class Toplevel -relief groove 
        wm focusmodel $base passive
        wm geometry $base 360x280+100+120
        wm maxsize $base 817 594
        wm minsize $base 1 1
        wm overrideredirect $base 0
        wm resizable $base 1 1
        wm deiconify $base
        wm title $base "Sort conditions window"
        set Delta 30
    #
        set Pref(Fill) yellow
    #
    # fill a list with entries
    #
        set Eintraege [list ]
        for {set i 1} {$i < 21} {incr i} {
            lappend Eintraege "Item $i"
        }
    #
    # One call creates the listbox
    #
        set Canv [Sort-by-Drag_Listbox .dialog $Eintraege 20 20 150 240]
    #
    # Widgets on the right side
    #
        label $base.xl -text X -anchor e
        label $base.yl -text Y -anchor e
        entry $base.x -textvariable PosX -width 12 -justify center
        entry $base.y -textvariable PosY -width 12 -justify center
        button $base.ok -text Quit -command exit -width 12 -default active
        label $base.dl -text Delta -anchor e
        entry $base.d -textvariable Delta -width 12 -justify center
        message $base.m -width 140 \
            -text "Sort the of conditions."
    #
    # Position all widgets
    #
        place $base.xl -x 236 -y  30 -anchor e
        place $base.yl -x 236 -y  60 -anchor e
        place $base.x  -x 240 -y  30 -anchor w
        place $base.y  -x 240 -y  60 -anchor w
        place $base.dl -x 236 -y 110 -anchor e
        place $base.d  -x 240 -y 110 -anchor w
        place $base.m  -x 200 -y 135 -anchor nw
        place $base.ok -x 220 -y 250 -anchor w
    }

    #
    # Create a pseudo-listbox with canvas elements. Looks like a listbox,
    # but is really a canvas, and all widgets only pretend to be what they seem.
    #
    proc Pfem::xml::CndSortWindow::Sort-by-Drag_Listbox { base Eintraege XNull YNull Breite Hoehe } {
        variable Pref
        variable Index
        variable Eintrag
        variable Scrollposition
        variable Scrollbereich
    #
        set Canv [canvas $base.cv -borderwidth 0 -highlightthickness 0 \
            -height [expr $Hoehe + 2*$YNull] -width [expr $Breite + 2*$XNull] \
            -bg $Pref(Fill)]
    #
    # Create the box with a scrollbar on the right
    #
        $Canv create rectangle $XNull $YNull [expr $XNull + $Breite] \
            [expr $YNull + $Hoehe] -outline black -width 1 -fill white -tags Box
        $Canv create rectangle [expr $XNull - 1] [expr $YNull - 1] \
            [expr $XNull + $Breite + 1] [expr $YNull + $Hoehe + 1] \
            -outline grey50 -width 1 -tags Box
        scrollbar $base.lbscroll -command "Pfem::xml::CndSortWindow::Sort-by-Drag_ListboxScroll $base" \
            -borderwidth 0 -orient vert -width 16 -cursor left_ptr
        place $Canv -x 0 -y 0 -anchor nw
        place $base.lbscroll -x [expr $XNull + $Breite - 16] -y $YNull \
            -anchor nw -width 16 -height $Hoehe
    #
    # Fill the box with the list
    #
        for {set i 0} {$i < [llength $Eintraege]} {incr i} {
            set Eintrag($i) "[lindex $Eintraege $i]"
            lappend Index($i) $i
        }
        set Scrollposition 0
        set Schrifthoehe 16
        set Scrollbereich [expr $Schrifthoehe * [llength $Eintraege]]
        Pfem::xml::CndSortWindow::Sort-by-Drag_ListboxScroll $base scroll 0.0 units
    #
        return $Canv
    }

    #
    # Scrollbar code.
    #
    proc Pfem::xml::CndSortWindow::Sort-by-Drag_ListboxScroll { base {was moveto} {Zahl 0.0} {Einheit units} } {
        variable Pref
        variable Index
        variable Eintrag
        variable Scrollposition
        variable Scrollbereich
    #
        set Canv  $base.cv
        set Hoehe [lindex [$Canv configure -height] 4]
        set Schrifthoehe 16
    #
        if {$was == "scroll"} {
            if {$Einheit == "pages"} {
                incr Scrollposition [expr int($Zahl * $Hoehe - 20)]
            } else {
                incr Scrollposition [expr 20 * int($Zahl)]
            }
        } else {
            set Scrollposition [expr int($Zahl * $Scrollbereich)]
        }
    #
    # Limit the scrollposition to sensible values
    #
        if {$Scrollposition > [expr $Scrollbereich - $Hoehe]} {
            set Scrollposition [expr $Scrollbereich - $Hoehe]
        }
        if {$Scrollposition < 0} {set Scrollposition 0}
    #
    # Delete Index and built anew from scratch. In priciple all entries could
    # be moved, but this is messy at the edges.
    #
        set yPos [expr 32 - $Scrollposition]
        for {set i 0} {$i < [array size Eintrag]} {incr i} {
            $Canv delete ent$i
            if {$yPos < 20} {
                incr yPos $Schrifthoehe
                continue
            }
    #
            if {$yPos < [expr $Hoehe - 20]} {
                $Canv create text 24 $yPos -text $Eintrag($Index($i)) \
                    -anchor w -fill black -tags ent$i
                incr yPos $Schrifthoehe
    #
    # Bindings for dragging and dropping of items.
    #
                $Canv bind ent$i <<Klick>>     "Pfem::xml::CndSortWindow::plotDown $Canv %x %y"
                $Canv bind ent$i <<Ziehen>>    "Pfem::xml::CndSortWindow::plotMove $Canv %x %y"
                $Canv bind ent$i <<Loslassen>> "Pfem::xml::CndSortWindow::plotCopy $base $Canv %x %y $i"
            }
        }
    #
        $base.lbscroll set [expr double($Scrollposition) / $Scrollbereich] [expr double($Hoehe + $Scrollposition) / $Scrollbereich]
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
    proc Pfem::xml::CndSortWindow::plotDown {w x y} {
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
    proc Pfem::xml::CndSortWindow::plotMove { w x y } {
        variable plot
        variable PosX
        variable PosY
    #
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
    proc Pfem::xml::CndSortWindow::plotCopy { base Cv x y i } {
        variable Pref
        variable Delta
        variable Index
        variable Scrollposition
        variable Scrollbereich
    #
        set Hoehe [lindex [$Cv configure -height] 4]
        set Schrifthoehe 16
    #
    # Get the new position. Delta is a fudge factor which is different
    # between different operating systems.
    #
        set Rang [expr int(($y - $Delta + $Scrollposition) / $Schrifthoehe)]
        set Speicher $Index($i)
        if {$Rang > $i} {
            for {set j $i} {$j < $Rang} {incr j} {
                set Index($j) $Index([expr $j + 1])
            }
        } elseif {$Rang == $i} {
            set Zahl [expr double($Scrollposition) / $Scrollbereich]
            Pfem::xml::CndSortWindow::Sort-by-Drag_ListboxScroll $base scroll $Zahl units
            return
        } else {
            set Rang [expr $Rang + 1]
            for {set j $i} {$j > $Rang} {incr j -1} {
                set Index($j) $Index([expr $j - 1])
            }
        }
        set Index($Rang) $Speicher
    #
    # Now scroll the list to the right position.
    #
        set Zahl [expr double($Scrollposition) / $Scrollbereich]
        Pfem::xml::CndSortWindow::Sort-by-Drag_ListboxScroll $base scroll $Zahl units
    }

    #
    Pfem::xml::CndSortWindow::Window show .gid.sortcondwind
    Pfem::xml::CndSortWindow::Window show .gid.dialog
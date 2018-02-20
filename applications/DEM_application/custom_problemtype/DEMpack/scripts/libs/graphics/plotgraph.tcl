###############################################################################
#
#    NAME: plotgraph.tcl
#
#    PURPOSE: Namespace to manage the graph plot using Plotchart package
#
#    COPYRIGHT
#  
#    QUANTECH ATZ-DEVELOPMENT DEPARTMENT
#    CIMNE-TTS   
#
#    AUTHOR : G. Socorro
#
#    CREATED AT: 21/08/12
#
#    HISTORY:
# 
#     0.2-20/09/12- J. Garate, Adaptation for Curves Interface
#     0.1-27/02/12-G. Socorro, create the base namespace KPlot
#
###############################################################################
#                      Procedures that belong to this namespace 
###############################################################################
#         Name                  |        Functionality
#------------------------------------------------------------------------------
# 1). 
#******************************************************************************
#                  BEGIN KPlot NAMESPACE
#******************************************************************************

package require Plotchart

namespace eval ::KPlot {
    
}

proc ::KPlot::Plot { c points title xtext ytext {symbol "plus"} } {    
     
    if { $points == "" || [llength $points]<=1 || ![winfo exists $c] } {
        return ""
    }    
    # foreach pt $points {
    # 	foreach coord $pt {
    # 	    if { ![string is double $coord] } {
    # 		WarnWin [= "value '%s' is not a number" $coord]
    # 		return
    # 	    }
    # 	}
    # } 
    
          
    foreach margin {top bottom left right} size { 15 15 0 15 } {
        ::Plotchart::plotconfig xyplot margin $margin $size
    }
    
    $c configure -background [Plotchart::plotconfig xyplot background innercolor]
    
    set xformat "%.3f" ;# -ticklines true -ticklength 3 -scale $xaxis
    set yformat "%.3f" ;# -ticklines true -ticklength 3
    
    set xmin 1e20
    set xmax -1e20
    set ymin 1e20
    set ymax -1e20
    
    # Calculamos los maximos y minimos en X e Y
    foreach pt $points {
        set x [lindex $pt 0]
        if { $xmax < $x } {
            set xmax $x
        }
        if { $xmin > $x } {
            set xmin $x
        }
        foreach y [lrange $pt 1 end] {
            if { $ymax < $y } {
                set ymax $y
            }
            if { $ymin > $y } {
                set ymin $y
            }
        }
    }
    
    if { $xmax <= $xmin } {
        #e.g all points with same x, illegal but...
        set xmin [expr $xmin-1.0]
        set xmax [expr $xmax+1.0]
    }
    if { $ymax <= $ymin } {
        #e.g all points with same y
        set ymin [expr $ymin-1.0]
        set ymax [expr $ymax+1.0]
    }
    
    set nincre [expr [llength $points]/2]
    incr nincre -1
    if { $nincre > 10 } {
        set nincre 10
    } elseif { $nincre < 1 } {
        set nincre 1
    }
    
    set xstep [expr {($xmax-$xmin)/double($nincre)}]
    set ystep [expr {($ymax-$ymin)/double($nincre)}]
    
    set xaxis [list $xmin $xmax $xstep]
    set yaxis [list $ymin $ymax $ystep]
    
    $c delete all
    
    set xyplot [::Plotchart::createXYPlot $c [list $xmin $xmax $xstep] [list $ymin $ymax $ystep]]
    
    $xyplot title $title
    $xyplot xtext $xtext
    $xyplot ytext [join $ytext ,]
    # When TCL 8.6
    #$xyplot vtext [join $ytext ,]
    $xyplot xconfig -format $xformat ;# -ticklines true -ticklength 3 -scale $xaxis
    $xyplot yconfig -format $yformat ;# -ticklines true -ticklength 3
    #$xyplot grid {{1 2 3} {1 2 3}}  {{20 20 20} {21 21 21}}
    #type: line, symbol, or both.
    #symbol: plus, cross, circle, up, down, dot, upfilled or downfilled
    set colors {red green blue yellow magenta cyan orange gray black}
    set ncolors [llength $colors]
    set nseries [expr {[llength [lindex $points 0]]-1}]    
    for {set serie 0} {$serie<$nseries} {incr serie} {
        $xyplot dataconfig $serie -colour [lindex $colors [expr {$serie%$ncolors}]] -type line       
        $xyplot legend $serie [lindex $ytext $serie]
    }
    #$xyplot xconfig -scale $xaxis
    #$xyplot yconfig -scale $yaxis
    
    foreach pt $points {
        set x [lindex $pt 0]
        set serie 0
        foreach y [lrange $pt 1 end] {
            $xyplot plot $serie $x $y
            incr serie
        }
    }
}

#******************************************************************************
#                  END KPlot NAMESPACE
#******************************************************************************

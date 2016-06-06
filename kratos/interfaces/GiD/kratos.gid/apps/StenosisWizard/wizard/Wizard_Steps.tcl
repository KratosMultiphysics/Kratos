
namespace eval StenosisWizard::Wizard {
    # Namespace variables declaration

}

proc StenosisWizard::Wizard::Init { } {
#W "Carga los pasos"
}

proc StenosisWizard::Wizard::Geometry { win } {
# Units
     set mtxt [= "m"]
     set mmtxt [= "mm"]
     set entrywidth 8
     
     set properties [Wizard::GetStepProperties Geometry]
     #W $properties
     
     # Set the widgets
     set imgname [Wizard::GetProperty Geometry ImageGeom,value]
     set img1 [apps::getImgFrom StenosisWizard $imgname]
     
     # Left frame
     set fr1 [ttk::frame $win.fr1 -borderwidth 10]
     # Rigth frame
     set fr2 [ttk::frame $win.fr2 -borderwidth 10]

     set labIm [ttk::label $fr2.lIm -image $img1]
     
     # Geom frame
     set labfr1 [ttk::labelframe $fr1.lfr1 -text [= "Define geometrical data"] -padding 10 -width 200 -height 200 ]
     
     set i 1
     set listids [list ]
     foreach prop $properties {
        if {$prop ni [list "Active" "Visited" "State" "ImageGeom"]} {
            set j [Wizard::GetProperty Geometry $prop,order]
            lappend listids $j
            set txt [= $prop]
            set lab$j [ttk::label $labfr1.l$j -text "${txt}:"]
            set ent$j [ttk::entry $labfr1.e$j -textvariable ::Wizard::wprops(Geometry,$prop,value) -width $entrywidth]
            wcb::callback $labfr1.e$j before insert wcb::checkEntryForReal
            #wcb::callback $labfr1.e$j after insert "StenosisWizard::Wizard::callbackCheckGeom"
            #wcb::callback $labfr1.e$j after delete "StenosisWizard::Wizard::callbackCheckGeom"
            #set labun$i [ttk::label $labfr1.lu$i -text "$mmtxt"]
            set txt [= "Enter a value for $txt"]
            tooltip::tooltip $labfr1.e$j "${txt}."
            incr i
        }
     }
     set listids [lsort $listids]
     # Widget geometry     
     # Grid the frames
     grid $fr1 -column 1 -row 0 -sticky nw
     grid $fr2 -column 2 -row 0 -sticky ne
     
     grid $labIm -column 0 -row 0 -sticky ne -rowspan 3
         
     # Label frames
     grid $labfr1 -column 1 -row 0 -sticky wen -ipadx 2  -columnspan 2
     
     # Label frame 1 => Geometrical data frame

     for {set i 0} {$i < [llength $listids]} {incr i} {
        # Cada uno
        set j [lindex $listids $i]
        set lab "lab$j"
        set ent "ent$j"
        #set labun "labun$i"
        #set units "units$i"
        
        grid [expr $$lab] -column 1 -row $i -sticky w -pady 2
        grid [expr $$ent] -column 2 -row $i -sticky w -pady 2
        #grid [expr $$labun] -column 3 -row $i -sticky w
        #grid [expr $$units] -column 4 -row $i -sticky w
          
     }
    set drawButton [button $labfr1.b2 -text [= "Draw Geometry"] -command [list ::StenosisWizard::Wizard::DrawGeometry]]
    $drawButton configure -bg #74bb92
    grid $drawButton -column 1 -columnspan 4 -row 9 -sticky ew
    
    grid columnconfigure $labfr1 1 -minsize 120
    
}
proc StenosisWizard::Wizard::NextGeometry { } {

}

proc StenosisWizard::Wizard::DrawGeometry {} {
    Kratos::ResetModel
    
    set err [ValidateDraw]
    if {$err ne 0} {
        return ""
    }
    set length [ Wizard::GetProperty Geometry Length,value]
    set radius [ Wizard::GetProperty Geometry Radius,value]
    set start [expr [ Wizard::GetProperty Geometry Z,value] *-1.0]
    set end [ Wizard::GetProperty Geometry Z,value]
    set delta [ Wizard::GetProperty Geometry Delta,value]
    set precision [ Wizard::GetProperty Geometry Precision,value]
    #W "Drawing tube: \nLength $length \nStart $start \nEnd $end \nDelta $delta \nPrecision $precision"
    set points [list]
    
    set layer [GiD_Info Project LayerToUse]
    GiD_Process 'Layers Color $layer 153036015 Transparent $layer 255 escape 

    set zona [expr $end - $start]
    set delta_z [expr double($zona) / double($precision)]
    
    # Initial point
    lappend points [list -$length $radius 0]
    GiD_Geometry create point 1 $layer -$length $radius 0
    
    # first cut
    lappend points [list $start $radius 0]
    #W $points
    
    
    for {set i [expr $start + $delta_z]} {$i < [expr $end - $delta_z]} {set i [expr $i + $delta_z]} {
        set y $radius
        set y [expr double($radius)-((double($delta)/2.0)*(1.0+cos(3.14*$i/double($end))))]
        #W "$i $y"
        lappend points [list $i $y 0]
    }
    
    # last cut
    lappend points [list $end $radius 0]
    # Final point
    GiD_Geometry create point 2 $layer $length $radius 0
    lappend points [list $length $radius 0]
    
    set line [GiD_Geometry create line append nurbsline $layer 1 2 -interpolate [llength $points] {*}$points -tangents {1 0 0} {1 0 0}]

    # Time to Revolute!
    GiD_Process Mescape Utilities Id $line escape escape Mescape Utilities Copy Lines DoExtrude Surfaces MaintainLayers MCopy 2 Rotation FNoJoin -$length,0.0,2.0 FNoJoin $length,0.0,2.0 180 1 escape
    
    # Closing tapas!
    GiD_Process Mescape Geometry Create NurbsSurface 3 5 escape 6 4 escape escape 

    # Volumenizando!
    GiD_Process Mescape Geometry Create volume 1 2 3 4 escape 
    
    # Agrupando
    GiD_Groups create Inlet
    GiD_EntitiesGroups assign Inlet surfaces {3}
    
    GiD_Groups create Outlet
    GiD_EntitiesGroups assign Outlet surfaces {4}
    
    GiD_Groups create NoSlip
    GiD_EntitiesGroups assign NoSlip surfaces {1 2}
    
    GiD_Groups create Fluid
    GiD_EntitiesGroups assign Fluid volumes {1}
    
    GidUtils::UpdateWindow GROUPS
    
    GiD_Process 'Zoom Frame escape
}

proc ValidateDraw { } {
    return 0
}

proc StenosisWizard::Wizard::Material { win } {

}
proc StenosisWizard::Wizard::NextMaterial { } {

}


proc StenosisWizard::Wizard::Fluid { win } {

}
proc StenosisWizard::Wizard::NextFluid { } {

}


proc StenosisWizard::Wizard::Simulation { win } {

}
proc StenosisWizard::Wizard::NextSimulation { } {

}

StenosisWizard::Wizard::Init


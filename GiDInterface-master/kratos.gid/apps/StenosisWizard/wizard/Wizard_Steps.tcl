
namespace eval StenosisWizard::Wizard {
    # Namespace variables declaration

}

proc StenosisWizard::Wizard::Init { } {
#W "Carga los pasos"
}

proc StenosisWizard::Wizard::Geometry { win } {
     Wizard::SetWindowSize 650 500
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
            set pn [Wizard::GetProperty Geometry $prop,name]
            lappend listids $j
            set txt [= $pn]
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
    
    # Partimos las superficies para refinar el mallado en el centro
    GiD_Process Mescape Geometry Edit DivideSurf NumDivisions 2 USense 3 escape escape
    GiD_Process Mescape Geometry Edit DivideSurf NumDivisions 1 USense 3 escape escape

}

proc ValidateDraw { } {
    return 0
}

proc StenosisWizard::Wizard::Material { win } {
     Wizard::SetWindowSize 300 400
     set entrywidth 16
     
     set properties [Wizard::GetStepProperties Material]
     #W $properties
     
     # Set the widgets
     #set imgname [Wizard::GetProperty Geometry ImageGeom,value]
     #set img1 [apps::getImgFrom StenosisWizard $imgname]
     
     # Left frame
     set fr1 [ttk::frame $win.fr1 -borderwidth 10]
     # Rigth frame
     #set fr2 [ttk::frame $win.fr2 -borderwidth 10]

     #set labIm [ttk::label $fr2.lIm -image $img1]
     
     # Geom frame
     set labfr1 [ttk::labelframe $fr1.lfr1 -text [= "Define material data"] -padding 10 -width 200 -height 200 ]
     
     set i 1
     set listids [list ]
     foreach prop $properties {
        if {$prop ni [list "Active" "Visited" "State"]} {
            set j [Wizard::GetProperty Material $prop,order]
            set pn [Wizard::GetProperty Material $prop,name]
            lappend listids $j
            set txt [= $pn]
            set lab$j [ttk::label $labfr1.l$j -text "${txt}:"]
            set type [Wizard::GetProperty Material $prop,type]
            if {$type eq "combo"} {
               set values [Wizard::GetProperty Material $prop,values]
               set ent$j [ttk::combobox $labfr1.e$j -values $values -textvariable ::Wizard::wprops(Material,$prop,value) -width $entrywidth -state readonly]
               bind $labfr1.e$j <<ComboboxSelected>> [list StenosisWizard::Wizard::ChangeFluidType %W] 
            } {
               set ent$j [ttk::entry $labfr1.e$j -textvariable ::Wizard::wprops(Material,$prop,value) -width $entrywidth]
            }
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
     #grid $fr2 -column 2 -row 0 -sticky ne
     
     #grid $labIm -column 0 -row 0 -sticky ne -rowspan 3
         
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
    
    grid columnconfigure $labfr1 1 -minsize 120
}
proc StenosisWizard::Wizard::NextMaterial { } {
     # Quitar parts existentes
     gid_groups_conds::delete {container[@n='StenosisWizard']/condition[@n='Parts']/group}
     
     # Crear una part con los datos que toquen
     set where {container[@n='StenosisWizard']/condition[@n='Parts']} 
     set gnode [spdAux::AddConditionGroupOnXPath $where "Fluid"]
     
     set props [list ConstitutiveLaw DENSITY VISCOSITY YIELD_STRESS POWER_LAW_K POWER_LAW_N]
     foreach prop $props {
          set propnode [$gnode selectNodes "./value\[@n = '$prop'\]"]
          if {$propnode ne "" } {
               $propnode setAttribute v [Wizard::GetProperty Material ${prop},value]
          }
     }
     spdAux::RequestRefresh
}

proc StenosisWizard::Wizard::ChangeFluidType {win} {
     W "Hey $win"
}


proc StenosisWizard::Wizard::Fluid { win } {
     Wizard::SetWindowSize 450 450
     set entrywidth 8
     
     set properties [Wizard::GetStepProperties Fluid]
     #W $properties
     
     # Set the widgets
     set imgname [Wizard::GetProperty Fluid ImageFluid,value]
     set img1 [apps::getImgFrom StenosisWizard $imgname]
     
     # Left frame
     set fr1 [ttk::frame $win.fr1 -borderwidth 10]
     # Rigth frame
     set fr2 [ttk::frame $win.fr2 -borderwidth 10]

     set labIm [ttk::label $fr2.lIm -image $img1]
     
     # Geom frame
     set labfr1 [ttk::labelframe $fr1.lfr1 -text [= "Define condition data"] -padding 10 -width 200 -height 200 ]
     
     set i 1
     set listids [list ]
     foreach prop $properties {
        if {$prop ni [list "Active" "Visited" "State" "ImageFluid"]} {
            set j [Wizard::GetProperty Fluid $prop,order]
            set pn [Wizard::GetProperty Fluid $prop,name]
            lappend listids $j
            set txt [= $pn]
            set lab$j [ttk::label $labfr1.l$j -text "${txt}:"]
            set ent$j [ttk::entry $labfr1.e$j -textvariable ::Wizard::wprops(Fluid,$prop,value) -width $entrywidth]
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
    
    grid columnconfigure $labfr1 1 -minsize 120
}
proc StenosisWizard::Wizard::NextFluid { } {
     # Inlet
     gid_groups_conds::delete {container[@n='StenosisWizard']/container[@n='BoundaryConditions']/condition[@n='Inlet3D']/group}
     set where {container[@n='StenosisWizard']/container[@n='BoundaryConditions']/condition[@n='AutomaticInlet3D']}
     set gnode [spdAux::AddConditionGroupOnXPath $where "Inlet"]
     [$gnode selectNodes "./value\[@n = 'direction'\]"] setAttribute v "automatic_inwards_normal"

     set propnode [$gnode selectNodes "./value\[@n = 'modulus'\]"]
     $propnode setAttribute v [Wizard::GetProperty Fluid Inlet,value]
     
     
          
     # Outlet
     gid_groups_conds::delete {container[@n='StenosisWizard']/container[@n='BoundaryConditions']/condition[@n='Outlet3D']/group}
     set where {container[@n='StenosisWizard']/container[@n='BoundaryConditions']/condition[@n='Outlet3D']}
     set gnode [spdAux::AddConditionGroupOnXPath $where "Outlet"]
     set propnode [$gnode selectNodes "./value\[@n = 'value'\]"]
     $propnode setAttribute v [Wizard::GetProperty Fluid Outlet,value]
     
     
     gid_groups_conds::delete {container[@n='StenosisWizard']/container[@n='BoundaryConditions']/condition[@n='NoSlip3D']/group}
     set where {container[@n='StenosisWizard']/container[@n='BoundaryConditions']/condition[@n='NoSlip3D']}
     set gnode [spdAux::AddConditionGroupOnXPath $where "NoSlip"]
     spdAux::RequestRefresh
}


proc StenosisWizard::Wizard::Simulation { win } {
     Wizard::SetWindowSize 450 450
     set entrywidth 8
     
     set properties [Wizard::GetStepProperties Simulation]
     #W $properties
     
     # Set the widgets
     set imgname [Wizard::GetProperty Simulation ImageSimulation,value]
     set img1 [apps::getImgFrom StenosisWizard $imgname]
     
     # Left frame
     set fr1 [ttk::frame $win.fr1 -borderwidth 10]
     # Rigth frame
     set fr2 [ttk::frame $win.fr2 -borderwidth 10]

     set labIm [ttk::label $fr2.lIm -image $img1]
     
     # Geom frame
     set labfr1 [ttk::labelframe $fr1.lfr1 -text [= "Define simulation data"] -padding 10 -width 200 -height 200 ]
     
     set i 1
     set listids [list ]
     foreach prop $properties {
        if {$prop ni [list "Active" "Visited" "State" "ImageSimulation"]} {
            set j [Wizard::GetProperty Simulation $prop,order]
            set pn [Wizard::GetProperty Simulation $prop,name]
            lappend listids $j
            set txt [= $pn]
            set lab$j [ttk::label $labfr1.l$j -text "${txt}:"]
            set ent$j [ttk::entry $labfr1.e$j -textvariable ::Wizard::wprops(Simulation,$prop,value) -width $entrywidth]
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
     set texto "Mesh"
     set command "::StenosisWizard::Wizard::Mesh"
     if {[GiD_Info Project ModelName] eq "UNNAMED"} {
          set texto "Save Model"
          set command "::StenosisWizard::Wizard::Save"
     }
     
     set drawButton [button $labfr1.b2 -text [= $texto] -command [list $command]]
     $drawButton configure -bg #74bb92
     grid $drawButton -column 1 -columnspan 4 -row [llength $listids] -sticky ew
     
     set state disabled
     if {![GiD_Info Project AreChanges]} {
          set state normal
     }
     
     set RunButton [button $labfr1.b3 -text [= "Run"] -command [list StenosisWizard::Wizard::Run] -state $state]
     $RunButton configure -bg #74bb92
     grid $RunButton -column 1 -columnspan 4 -row [expr [llength $listids] +1] -sticky ew
     
     grid columnconfigure $labfr1 1 -minsize 120
}

proc StenosisWizard::Wizard::Mesh { } {
     if {[lindex [GiD_Info Mesh] 0]>0} {
          #GiD_Process Mescape Meshing reset Yes
          GiD_Process Mescape Meshing CancelMesh PreserveFrozen Yes
     }
     
     set mesh [Wizard::GetProperty Simulation MeshSize,value]
     set meshinner [Wizard::GetProperty Simulation CentralMesh,value]
     GiD_Process Mescape Meshing AssignSizes Surfaces $meshinner 6 9 escape escape
     GiD_Process Mescape Meshing Generate $mesh MeshingParametersFrom=Preferences Mescape Meshing MeshView
}
proc StenosisWizard::Wizard::Run { } {
     set root [customlib::GetBaseRoot]
     set solstrat_un [apps::getCurrentUniqueName SolStrat]
     #W $solstrat_un
     if {[get_domnode_attribute [$root selectNodes [spdAux::getRoute $solstrat_un]] v] eq ""} {
         get_domnode_attribute [$root selectNodes [spdAux::getRoute $solstrat_un]] dict
     }
     set solstrat_un [apps::getCurrentUniqueName Scheme]
     #W $solstrat_un
     if {[get_domnode_attribute [$root selectNodes [spdAux::getRoute $solstrat_un]] v] eq ""} {
         get_domnode_attribute [$root selectNodes [spdAux::getRoute $solstrat_un]] dict
     }
     LastStep
     GiD_Process Mescape Utilities Calculate
}

proc StenosisWizard::Wizard::LastStep { } {
     set i [Wizard::GetProperty Simulation InitialTime,value]
     set e [Wizard::GetProperty Simulation EndTime,value]
     set d [Wizard::GetProperty Simulation DeltaTime,value]
     
     gid_groups_conds::setAttributesF {container[@n='StenosisWizard']/container[@n='SolutionStrat']/container[@n='TimeParameters']/value[@n='InitialTime']} "v $i"
     gid_groups_conds::setAttributesF {container[@n='StenosisWizard']/container[@n='SolutionStrat']/container[@n='TimeParameters']/value[@n='EndTime']} "v $e"
     gid_groups_conds::setAttributesF {container[@n='StenosisWizard']/container[@n='SolutionStrat']/container[@n='TimeParameters']/value[@n='DeltaTime']} "v $d"

     gid_groups_conds::copyNode {container[@n='StenosisWizard']/container[@n='Results']/container[@n='CutPlanes']/blockdata[@n='CutPlane'][1]} {container[@n='StenosisWizard']/container[@n='Results']/container[@n='CutPlanes']}
     
     gid_groups_conds::setAttributesF {container[@n='StenosisWizard']/container[@n='Results']/container[@n='CutPlanes']/blockdata[@n='CutPlane'][2]} {name Main}
     gid_groups_conds::setAttributesF {container[@n='StenosisWizard']/container[@n='Results']/container[@n='CutPlanes']/blockdata[@n='CutPlane'][2]/value[@n='normal']} {v 0.0,1.0,0.0}
     
     set ncuts [Wizard::GetProperty Simulation Cuts,value]
     set length [Wizard::GetProperty Geometry Length,value]
     set delta [expr 2.0*double($length)/(double($ncuts)+1.0)]
     #W "$length $delta"
     for {set i 1} {$i <= $ncuts} {incr i} {
          set x [expr -$length + ($i * $delta)]
          set x [expr double(round(100*$x))/100]
          gid_groups_conds::copyNode {container[@n='StenosisWizard']/container[@n='Results']/container[@n='CutPlanes']/blockdata[@n='CutPlane'][1]} {container[@n='StenosisWizard']/container[@n='Results']/container[@n='CutPlanes']}
          set cutplane "container\[@n='StenosisWizard'\]/container\[@n='Results'\]/container\[@n='CutPlanes'\]/blockdata\[@n='CutPlane'\]\[[expr $i +2]\]"
          gid_groups_conds::setAttributesF $cutplane "name CutPlane$i"
          gid_groups_conds::setAttributesF "$cutplane/value\[@n='normal'\]" "v 1.0,0.0,0.0"
          gid_groups_conds::setAttributesF "$cutplane/value\[@n='point'\]" "v $x,0.0,0.0"
     }

     
     gid_groups_conds::setAttributesF {container[@n='StenosisWizard']/container[@n='SolutionStrat']/container[@n='velocity_linear_solver_settings']/value[@n='Solver']} {v Conjugate_gradient}
     gid_groups_conds::setAttributesF {container[@n='StenosisWizard']/container[@n='SolutionStrat']/container[@n='pressure_linear_solver_settings']/value[@n='Solver']} {v Conjugate_gradient}
     spdAux::RequestRefresh

}
proc StenosisWizard::Wizard::Save { } {
     GiD_Process Mescape Files Save 
     if {[GiD_Info Project ModelName] ne "UNNAMED"} {
          set drawButton .gid.activewizard.w.layoutFrame.wiz.layout.fr1.lfr1.b2
          $drawButton configure -text "Mesh"
          $drawButton configure -command [list "::StenosisWizard::Wizard::Mesh"]
          $drawButton configure -bg #1e90ff
     }
}

proc StenosisWizard::AfterMeshGeneration { fail } {
     GiD_Process Mescape Mescape Mescape
     GiD_Process Mescape Files Save 
     if {$fail eq 0} {
          set Button .gid.activewizard.w.layoutFrame.wiz.layout.fr1.lfr1.b3
          if {[winfo exists $Button]} {$Button configure -state normal}
     }
}

StenosisWizard::Wizard::Init


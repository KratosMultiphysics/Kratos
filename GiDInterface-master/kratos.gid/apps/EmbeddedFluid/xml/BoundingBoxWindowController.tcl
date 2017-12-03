namespace eval EmbeddedFluid::xml::BoundingBox {
    variable winpath
    variable box
    variable boxname
}

proc EmbeddedFluid::xml::BoundingBox::Init {} {
    variable winpath
    set winpath ".gid.bboxwindow"
    variable box
    set box(x1) 1.0
    set box(y1) 1.0
    set box(z1) 1.0
    set box(x2) 1.0
    set box(y2) 1.0
    set box(z2) 1.0
    variable boxname
    set boxname "FluidBox"
}

proc EmbeddedFluid::xml::BoundingBox::CreateWindow { } {
    variable winpath
    #Init
    set w $winpath
    
    if {[winfo exists $w]} {destroy $w}
    
    toplevel $w -class Toplevel -relief groove 
    #wm maxsize $w 500 300
    wm minsize $w 250 250
    wm overrideredirect $w 0
    wm resizable $w 1 1
    wm deiconify $w
    wm title $w [= "Bounding box window"]
    wm attribute $w -topmost 1
    
    EmbeddedFluid::xml::BoundingBox::RefreshWindow 
}

proc EmbeddedFluid::xml::BoundingBox::RefreshWindow { } {
    variable winpath
    set w $winpath
    if {[winfo exists $w.fr1]} {destroy $w.fr1}
    if {[winfo exists $w.fr2]} {destroy $w.fr2}
    if {[winfo exists $w.buts]} {destroy $w.buts}
    set fr1 [ttk::labelframe $w.fr1 -text [= "Model space"]]
    set fr2 [ttk::labelframe $w.fr2 -text [= "Create a bounding box"]]
    set buts [ttk::frame $w.buts -style BottomFrame.TFrame]

    ttk::button $buts.q -text Close -command [list destroy $w] -style BottomFrame.TButton
    ttk::button $buts.ok -text "Create Box" -command [list EmbeddedFluid::xml::BoundingBox::BuildBox] -style BottomFrame.TButton
    EmbeddedFluid::xml::BoundingBox::ModelInfoFrame $fr1
    EmbeddedFluid::xml::BoundingBox::BoxDataFrame $fr2
    
    grid $buts.ok $buts.q -sticky sew
    
    grid $fr1 -sticky nsew -row 0 -column 0
    grid $fr2 -sticky nsew -row 1 -column 0
    grid $buts -sticky sew -columnspan 2
    grid columnconfigure $w 0 -weight 1 
    grid rowconfigure $w 2 -weight 1
    if { $::tcl_version >= 8.5 } { grid anchor $buts center }
}

proc EmbeddedFluid::xml::BoundingBox::ModelInfoFrame {w} {
    lassign [GetCurrentBox] x2 y2 z2 x1 y1 z1
    set x1 [format %.3g $x1]
    set x2 [format %.3g $x2]
    set y1 [format %.3g $y1]
    set y2 [format %.3g $y2]
    set z1 [format %.3g $z1]
    set z2 [format %.3g $z2]
    
    set mintxt [= "Min"]
    set maxtxt [= "Max"]
    
    set tX [= "Model X"]
    set lX [ttk::label $w.lX -text ${tX}:]
    set lminX [ttk::label $w.lminX -text ${mintxt}:]
    set eminX [ttk::label $w.eminX -text $x1]
    set lmaxX [ttk::label $w.lmaxX -text ${maxtxt}:]
    set emaxX [ttk::label $w.emaxX -text $x2]
    
    set tY [= "Model Y"]
    set lY [ttk::label $w.lY -text ${tY}:]
    set lminY [ttk::label $w.lminY -text ${mintxt}:]
    set eminY [ttk::label $w.eminY -text $y1]
    set lmaxY [ttk::label $w.lmaxY -text ${maxtxt}:]
    set emaxY [ttk::label $w.emaxY -text $y2]
    
    set tZ [= "Model Z"]
    set lZ [ttk::label $w.lZ -text ${tZ}:]
    set lminZ [ttk::label $w.lminZ -text ${mintxt}:]
    set eminZ [ttk::label $w.eminZ -text $z1]
    set lmaxZ [ttk::label $w.lmaxZ -text ${maxtxt}:]
    set emaxZ [ttk::label $w.emaxZ -text $z2]
    
    grid $lX $lminX $eminX $lmaxX $emaxX -sticky w
    grid $lY $lminY $eminY $lmaxY $emaxY -sticky w
    grid $lZ $lminZ $eminZ $lmaxZ $emaxZ -sticky w
}

proc EmbeddedFluid::xml::BoundingBox::BoxDataFrame {w} {
    set tDist [= "Distance"]
    set size 8
    set lX [ttk::label $w.lX -text $tDist]
    set lminX [ttk::label $w.lminX -text "-X:"]
    set eminX [ttk::entry $w.eminX -textvariable EmbeddedFluid::xml::BoundingBox::box(x1)]
    wcb::callback $eminX before insert wcb::checkEntryForReal
    set lmaxX [ttk::label $w.lmaxX -text "+X:"]
    set emaxX [ttk::entry $w.emaxX -textvariable EmbeddedFluid::xml::BoundingBox::box(x2)]
    wcb::callback $emaxX before insert wcb::checkEntryForReal
    
    set lY [ttk::label $w.lY -text $tDist]
    set lminY [ttk::label $w.lminY -text "-Y:"]
    set eminY [ttk::entry $w.eminY -textvariable EmbeddedFluid::xml::BoundingBox::box(y1)]
    wcb::callback $eminY before insert wcb::checkEntryForReal
    set lmaxY [ttk::label $w.lmaxY -text "+Y:"]
    set emaxY [ttk::entry $w.emaxY -textvariable EmbeddedFluid::xml::BoundingBox::box(y2)]
    wcb::callback $emaxY before insert wcb::checkEntryForReal
    
    set lZ [ttk::label $w.lZ -text $tDist]
    set lminZ [ttk::label $w.lminZ -text "-Z:"]
    set eminZ [ttk::entry $w.eminZ -textvariable EmbeddedFluid::xml::BoundingBox::box(z1)]
    wcb::callback $eminZ before insert wcb::checkEntryForReal
    set lmaxZ [ttk::label $w.lmaxZ -text "+Z:"]
    set emaxZ [ttk::entry $w.emaxZ -textvariable EmbeddedFluid::xml::BoundingBox::box(z2)]
    wcb::callback $emaxZ before insert wcb::checkEntryForReal
    
    grid $lX $lminX $eminX $lmaxX $emaxX -sticky e
    grid $lY $lminY $eminY $lmaxY $emaxY -sticky e
    grid $lZ $lminZ $eminZ $lmaxZ $emaxZ -sticky e
    
    grid configure $eminX $eminY $eminZ $emaxX $emaxY $emaxZ -sticky ew
    grid columnconfigure $w {2 4} -weight 1
}

proc EmbeddedFluid::xml::BoundingBox::DestroyBox { } {
    variable boxname
    
    if {[GiD_Layers exists $boxname]} {
        GiD_Layers delete $boxname
    }
}

proc EmbeddedFluid::xml::BoundingBox::BuildBox { } {
    variable boxname
    
    GidUtils::DisableGraphics
    
    EmbeddedFluid::xml::BoundingBox::DestroyBox
    
    GiD_Layers create $boxname
    GiD_Layers edit to_use $boxname
    
    EmbeddedFluid::xml::BoundingBox::CreateBoxGeom
    if {![GiD_Groups exists $boxname]} {GiD_Groups create $boxname}
    
    GiD_EntitiesGroups assign $boxname volumes  [GiD_EntitiesLayers get $boxname volumes]
    GiD_EntitiesGroups assign $boxname surfaces [GiD_EntitiesLayers get $boxname surfaces]
    GiD_EntitiesGroups assign $boxname lines    [GiD_EntitiesLayers get $boxname lines]
    GiD_EntitiesGroups assign $boxname points   [GiD_EntitiesLayers get $boxname points]
    
    GidUtils::EnableGraphics
    
    
    GiD_Process 'Zoom Frame
    GidUtils::UpdateWindow GROUPS
    GidUtils::UpdateWindow LAYER
}

proc EmbeddedFluid::xml::BoundingBox::CreateBoxGeom { } {
    variable boxname
    variable box
    lassign [GetCurrentBox] x2 y2 z2 x1 y1 z1
    
    set Xo [expr $x1 - $box(x1)]
    set Yo [expr $y1 - $box(y1)]
    set Zo [expr $z1 - $box(z1)]
    set Xf [expr $x2 + $box(x2)]
    set Yf [expr $y2 + $box(y2)]
    set Zf [expr $z2 + $box(z2)]
    
    set h [expr $Zf - $Zo]
    
    GiD_Process Mescape Geometry Create Object Rectangle $Xo $Yo $Zo $Xf $Yf $Zo Mescape
    set surfid [GiD_EntitiesLayers get $boxname surfaces]
    GiD_Process Mescape Utilities Copy Surfaces Duplicate DoExtrude Volumes MaintainLayers Translation FNoJoin 0.0 0.0 0.0 FNoJoin 0.0 0.0 $h $surfid Mescape 
    GidUtils::SetWarnLine [= "Bounding box created: Xo = %s  Yo = %s  Zo = %s  Xf = %s  Yf = %s  Zf = %s" $Xo $Yo $Zo $Xf $Yf $Zf]
}
    
proc EmbeddedFluid::xml::BoundingBox::GetCurrentBox { } {
    variable boxname
    set layers [lsearch -all -inline -not -exact [GiD_Layers list] $boxname]
    set modelbox [GiD_Info Layers -bbox {*}$layers]
    return {*}$modelbox
}

EmbeddedFluid::xml::BoundingBox::Init
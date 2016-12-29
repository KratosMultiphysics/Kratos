namespace eval ::Fluid {
    # Variable declaration
    variable dir
    variable prefix
    variable attributes
}

proc ::Fluid::Init { } {
    # Variable initialization
    variable dir
    variable prefix
    variable attributes
    
    set dir [apps::getMyDir "Fluid"]
    set attributes [dict create]
    
    set prefix FL
    set ::Model::ValidSpatialDimensions [list 2D 3D]
    
    # Allow to open the tree
    set ::spdAux::TreeVisibility 1
    
    dict set attributes UseIntervals 0
    if {$::Kratos::kratos_private(DevMode) eq "dev"} {dict set attributes UseIntervals 1}
    
    LoadMyFiles
    ::Fluid::FluidAppSelectorWindow
}

proc ::Fluid::LoadMyFiles { } {
    variable dir
    
    uplevel #0 [list source [file join $dir xml GetFromXML.tcl]]
    uplevel #0 [list source [file join $dir write write.tcl]]
    uplevel #0 [list source [file join $dir write writeProjectParameters.tcl]]
}

proc ::Fluid::GetAttribute {name} {
    variable attributes
    set value ""
    catch {set value [dict get $attributes $name]}
    return $value
}

proc ::Fluid::FluidAppSelectorWindow { } {
    set initwind $::spdAux::initwind
    set doc $gid_groups_conds::doc
    set root [$doc documentElement]
    set nd ""
    catch {set nd [ [$root selectNodes "value\[@n='nDim'\]"] getAttribute v]} 
    if {$nd eq ""} {catch {set nd [ [$root selectNodes "hiddenfield\[@n='nDim'\]"] getAttribute v]}}
    if { $nd ne "undefined" } {
        if {[apps::getActiveAppId] eq "Fluid"} {
            spdAux::SwitchDimAndCreateWindow $nd
        }
    } {
        [$root selectNodes "value\[@n='nDim'\]"] setAttribute v wait
        set dir $::Kratos::kratos_private(Path)
        
        set initwind .gid.win_example
        if { [ winfo exist $initwind]} {
            destroy $initwind
        }
        toplevel $initwind
        wm withdraw $initwind
        
        set w $initwind
        
        set x [expr [winfo rootx .gid]+[winfo width .gid]/2-[winfo width $w]/2]
        set y [expr [winfo rooty .gid]+[winfo height .gid]/2-[winfo height $w]/2]
        
        wm geom $initwind +$x+$y
        wm transient $initwind .gid    
        
        InitWindow $w [_ "Fluid applications"] Kratos "" "" 1
        set initwind $w
        ttk::frame $w.top
        ttk::label $w.top.title_text -text [_ "Select a fluid application"]
        
        ttk::frame $w.information  -relief ridge
        set i 0
        set apps [list Fluid EmbeddedFluid]
        foreach app $apps {
            set img [::apps::getImgFrom $app]
            set app_publicname [[::apps::getAppById $app] getPublicName]
            set but [ttk::button $w.information.img$app -image $img -command [list ::Fluid::ChangeAppTo $app] ]
            ttk::label $w.information.text$app -text $app_publicname
            grid $w.information.img$app -column $i -row 0
            grid $w.information.text$app -column $i -row 1
            incr i
        }
        grid $w.top
        grid $w.top.title_text
        
        grid $w.information
    }
}

proc ::Fluid::ChangeAppTo {appid} {
    switch $appid {
        "Fluid" {
            set doc $gid_groups_conds::doc
            set root [$doc documentElement]
            [$root selectNodes "value\[@n='nDim'\]"] setAttribute v undefined
            ::spdAux::CreateDimensionWindow
        }
        "EmbeddedFluid" {
            spdAux::deactiveApp Fluid
            apps::setActiveApp $appid       
        }
        default {
            ::spdAux::CreateDimensionWindow
        }
    }
    
}

::Fluid::Init

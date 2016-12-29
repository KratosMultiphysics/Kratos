namespace eval ::FSI {
    # Variable declaration
    variable dir
    variable prefix
    variable attributes
}

proc ::FSI::Init { } {
    # Variable initialization
    variable dir
    variable prefix
    #W "Sourced FSI"
    set dir [apps::getMyDir "FSI"]
    set prefix FSI
    
    apps::LoadAppById "Fluid"
    apps::LoadAppById "Structural"
    
    # Intervals only in developer mode
    dict set attributes UseIntervals 0
    #if {$::Kratos::kratos_private(DevMode) eq "dev"} {dict set attributes UseIntervals 1}
    # Allow to open the tree
    set ::spdAux::TreeVisibility 1
    
    set ::Model::ValidSpatialDimensions [list 2D 3D]
    LoadMyFiles
    #::spdAux::CreateDimensionWindow
}

proc ::FSI::LoadMyFiles { } {
    variable dir
    
    uplevel #0 [list source [file join $dir examples examples.tcl]]
    uplevel #0 [list source [file join $dir xml GetFromXML.tcl]]
    uplevel #0 [list source [file join $dir write write.tcl]]
    uplevel #0 [list source [file join $dir write writeProjectParameters.tcl]]
}

proc ::FSI::CustomToolbarItems { } {
    Kratos::ToolbarAddItem "Example" "example.png" [list -np- ::FSI::examples::MokChannelFlexibleWall] [= "Example\nMOK - Channel with flexible wall"]   
}

::FSI::Init

namespace eval ::FSI {
    # Variable declaration
    variable dir
    variable prefix
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
    
    set ::Model::ValidSpatialDimensions [list 2D 3D]
    ::spdAux::CreateDimensionWindow
    
    
    GiDMenu::InsertOption "Kratos" [list "---"] 6 PRE "" "" "" replace =
    GiDMenu::InsertOption "Kratos" [list "Mok - Channel with flexible wall" ] 7 PRE [list ::FSI::xml::MokChannelFlexibleWall] "" "" replace =
    GiDMenu::UpdateMenus
}

proc ::FSI::LoadMyFiles { } {
    variable dir
    
    uplevel #0 [list source [file join $dir xml GetFromXML.tcl]]
    uplevel #0 [list source [file join $dir write write.tcl]]
    uplevel #0 [list source [file join $dir write writeProjectParameters.tcl]]
}

::FSI::Init

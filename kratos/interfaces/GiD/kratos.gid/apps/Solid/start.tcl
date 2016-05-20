namespace eval ::Solid {
    # Variable declaration
    variable dir
}

proc ::Solid::Init { } {
    # Variable initialization
    variable dir
    
    set dir [apps::getMyDir "Solid"]
    set ::Model::ValidSpatialDimensions [list 2D 2Da 3D]
    ::spdAux::CreateDimensionWindow
}

proc ::Solid::LoadMyFiles { } {
    variable dir
    
    uplevel #0 [list source [file join $dir xml GetFromXML.tcl]]
    uplevel #0 [list source [file join $dir write write.tcl]]
    uplevel #0 [list source [file join $dir write writeProjectParameters.tcl]]
}

::Solid::Init

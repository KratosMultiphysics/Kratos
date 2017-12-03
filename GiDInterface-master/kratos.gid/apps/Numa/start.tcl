namespace eval ::Numa {
    # Variable declaration
    variable dir
    variable kratos_name
}

proc ::Numa::Init { } {
    # Variable initialization
    variable dir
    variable kratos_name
    set kratos_name "DamApplication"
    
    set dir [apps::getMyDir "Numa"]
    set ::Model::ValidSpatialDimensions [list 2D 3D]
    
    # Allow to open the tree
    set ::spdAux::TreeVisibility 1
    LoadMyFiles
    ::spdAux::CreateDimensionWindow
    
}

proc ::Numa::LoadMyFiles { } {
    variable dir
    
    uplevel #0 [list source [file join $dir xml GetFromXML.tcl]]
    uplevel #0 [list source [file join $dir write write.tcl]]
    uplevel #0 [list source [file join $dir write writeProjectParameters.tcl]]
}

::Numa::Init
namespace eval ::Dam {
    # Variable declaration
    variable dir
}

proc ::Dam::Init { } {
    # Variable initialization
    variable dir
    
    set dir [apps::getMyDir "Dam"]
    set ::Model::ValidSpatialDimensions [list 2D 3D]
    
    # Allow to open the tree
    set ::spdAux::TreeVisibility 1
    LoadMyFiles
    ::spdAux::CreateDimensionWindow
    
}

proc ::Dam::LoadMyFiles { } {
    variable dir
    
    uplevel #0 [list source [file join $dir xml GetFromXML.tcl]]
    uplevel #0 [list source [file join $dir write write.tcl]]
    uplevel #0 [list source [file join $dir write writeProjectParameters.tcl]]
}

::Dam::Init

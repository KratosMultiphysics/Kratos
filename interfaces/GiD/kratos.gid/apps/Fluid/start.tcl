namespace eval ::Fluid {
    # Variable declaration
    variable dir
    variable prefix
}

proc ::Fluid::Init { } {
    # Variable initialization
    variable dir
    variable prefix
    
    set dir [apps::getMyDir "Fluid"]
    set prefix FL
    set ::Model::ValidSpatialDimensions [list 2D 3D]
    ::spdAux::CreateDimensionWindow
}

proc ::Fluid::LoadMyFiles { } {
    variable dir
    
    uplevel #0 [list source [file join $dir xml GetFromXML.tcl]]
    uplevel #0 [list source [file join $dir write write.tcl]]
    uplevel #0 [list source [file join $dir write writeProjectParameters.tcl]]
}

::Fluid::Init

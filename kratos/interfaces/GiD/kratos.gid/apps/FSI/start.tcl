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
    set ::Model::ValidSpatialDimensions [list 2D 3D]
    ::spdAux::CreateDimensionWindow
}

proc ::FSI::LoadMyFiles { } {
    variable dir
    
    uplevel #0 [list source [file join $dir xml GetFromXML.tcl]]
    uplevel #0 [list source [file join $dir write write.tcl]]
    uplevel #0 [list source [file join $dir write writeProjectParameters.tcl]]
}

::FSI::Init

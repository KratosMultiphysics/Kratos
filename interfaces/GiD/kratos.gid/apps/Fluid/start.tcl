namespace eval ::Fluid {
    # Variable declaration
    variable dir
}

proc ::Fluid::Init { } {
    # Variable initialization
    variable dir
    
    set dir [apps::getMyDir "Fluid"]
}

proc ::Fluid::LoadMyFiles { } {
    variable dir
    
    uplevel #0 [list source [file join $dir xml GetFromXML.tcl]]
    uplevel #0 [list source [file join $dir write write.tcl]]
}

::Fluid::Init

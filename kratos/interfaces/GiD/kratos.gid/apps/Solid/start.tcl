namespace eval ::Solid {
    # Variable declaration
    variable dir
}

proc ::Solid::Init { } {
    # Variable initialization
    variable dir
    
    set dir [apps::getMyDir "Solid"]
}

proc ::Solid::LoadMyFiles { } {
    variable dir
    
    uplevel #0 [list source [file join $dir xml GetFromXML.tcl]]
    uplevel #0 [list source [file join $dir write write.tcl]]
}

::Solid::Init

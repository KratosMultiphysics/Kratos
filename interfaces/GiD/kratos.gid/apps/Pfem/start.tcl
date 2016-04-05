namespace eval ::Pfem {
    # Variable declaration
    variable dir
}

proc ::Pfem::Init { } {
    # Variable initialization
    variable dir
    
    set dir [apps::getMyDir "Pfem"]
}

proc ::Pfem::LoadMyFiles { } {
    variable dir
    
    uplevel #0 [list source [file join $dir xml GetFromXML.tcl]]
    uplevel #0 [list source [file join $dir write write.tcl]]
}

::Pfem::Init

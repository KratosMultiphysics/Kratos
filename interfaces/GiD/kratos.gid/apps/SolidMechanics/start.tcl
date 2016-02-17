namespace eval ::SolidMechanics {
    # Variable declaration
    variable dir
}

proc ::SolidMechanics::Init { } {
    # Variable initialization
    variable dir
    
    set dir [apps::getMyDir "SolidMechanics"]
}

proc ::SolidMechanics::LoadMyFiles { } {
    variable dir
    
    uplevel #0 [list source [file join $dir xml GetFromXML.tcl]]
    uplevel #0 [list source [file join $dir write write.tcl]]
}

::SolidMechanics::Init

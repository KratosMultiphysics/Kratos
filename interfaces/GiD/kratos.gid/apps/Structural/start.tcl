namespace eval ::Structural {
    # Variable declaration
    variable dir
}

proc ::Structural::Init { } {
    # Variable initialization
    variable dir
    
    set dir [apps::getMyDir "Structural"]
}

proc ::Structural::LoadMyFiles { } {
    variable dir
    
    uplevel #0 [list source [file join $dir xml GetFromXML.tcl]]
    uplevel #0 [list source [file join $dir write write.tcl]]
    uplevel #0 [list source [file join $dir write writeProjectParameters.tcl]]
}

proc ::Structural::MultiAppEvent {args} {
    W "$args"
}

::Structural::Init

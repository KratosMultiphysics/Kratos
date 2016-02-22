namespace eval ::BaseApp {
    # Variable declaration
    variable dir
}

proc ::BaseApp::Init { } {
    # Variable initialization
    variable dir
    
    set dir [apps::getMyDir "BaseApp"]
}

proc ::BaseApp::LoadMyFiles { } {
    variable dir
    
    uplevel #0 [list source [file join $dir xml GetFromXML.tcl]]
    uplevel #0 [list source [file join $dir write write.tcl]]
}

::BaseApp::Init

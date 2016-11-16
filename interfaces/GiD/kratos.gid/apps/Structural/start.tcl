namespace eval ::Structural {
    # Variable declaration
    variable dir
    variable attributes
}

proc ::Structural::Init { } {
    # Variable initialization
    variable dir
    variable attributes
    
    set dir [apps::getMyDir "Structural"]
    set attributes [dict create]
    # Allow to open the tree
    set ::spdAux::TreeVisibility 1
    apps::LoadAppById "Solid"
    
    LoadMyFiles
}

proc ::Structural::LoadMyFiles { } {
    variable dir
    
    uplevel #0 [list source [file join $dir xml GetFromXML.tcl]]
    uplevel #0 [list source [file join $dir write write.tcl]]
}

proc ::Structural::GetAttribute {name} {
    variable attributes
    set value ""
    catch {set value [dict get $attributes $name]}
    if {$value eq ""} {set value [::Solid::GetAttribute $name]}
    return $value
}

::Structural::Init

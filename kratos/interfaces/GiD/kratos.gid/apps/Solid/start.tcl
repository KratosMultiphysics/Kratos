namespace eval ::Solid {
    # Variable declaration
    variable dir
    variable attributes
}

proc ::Solid::Init { } {
    # Variable initialization
    variable dir
    variable attributes
    
    set dir [apps::getMyDir "Solid"]
    set ::Model::ValidSpatialDimensions [list 2D 2Da 3D]
    set attributes [dict create]
    
    # Intervals only in developer mode
    dict set attributes UseIntervals 0
    if {$::Kratos::kratos_private(DevMode) eq "dev"} {dict set attributes UseIntervals 1}
    
    # Restart available
    dict set attributes UseRestart 1
    # Allow to open the tree
    set ::spdAux::TreeVisibility 1
    LoadMyFiles
    ::spdAux::CreateDimensionWindow
}

proc ::Solid::LoadMyFiles { } {
    variable dir
    
    uplevel #0 [list source [file join $dir xml GetFromXML.tcl]]
    uplevel #0 [list source [file join $dir write write.tcl]]
    uplevel #0 [list source [file join $dir write writeProjectParameters.tcl]]
}

proc ::Solid::GetAttribute {name} {
    variable attributes
    set value ""
    catch {set value [dict get $attributes $name]}
    return $value
}

::Solid::Init

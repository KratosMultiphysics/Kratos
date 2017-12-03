namespace eval ::Structural {
    # Variable declaration
    variable dir
    variable attributes
    variable kratos_name
}

proc ::Structural::Init { } {
    # Variable initialization
    variable dir
    variable attributes
    variable kratos_name
    
    set dir [apps::getMyDir "Structural"]
    set attributes [dict create]
    
    # Allow to open the tree
    set ::spdAux::TreeVisibility 1
    
    # Intervals 
    dict set attributes UseIntervals 1
    if {$::Kratos::kratos_private(DevMode) eq "dev"} {dict set attributes UseIntervals 1}
    
    set kratos_name StructuralMechanicsApplication
    
    LoadMyFiles
}

proc ::Structural::LoadMyFiles { } {
    variable dir
    
    uplevel #0 [list source [file join $dir examples examples.tcl]]
    uplevel #0 [list source [file join $dir xml GetFromXML.tcl]]
    uplevel #0 [list source [file join $dir write write.tcl]]
    uplevel #0 [list source [file join $dir write writeProjectParameters.tcl]]
}


proc ::Structural::CustomToolbarItems { } {
    Kratos::ToolbarAddItem "Example" "example.png" [list -np- ::Structural::examples::TrussCantilever] [= "Example\nTruss cantilever"]   
}

proc ::Structural::GetAttribute {name} {
    variable attributes
    set value ""
    if {[dict exists $attributes $name]} {set value [dict get $attributes $name]}
    return $value
}

::Structural::Init

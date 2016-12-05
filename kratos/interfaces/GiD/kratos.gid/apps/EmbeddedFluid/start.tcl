namespace eval ::EmbeddedFluid {
    # Variable declaration
    variable dir
    variable prefix
    variable attributes
}

proc ::EmbeddedFluid::Init { } {
    # Variable initialization
    variable dir
    variable prefix
    variable attributes
    
    set dir [apps::getMyDir "EmbeddedFluid"]
    set attributes [dict create]
    
    set prefix EMBFL
    set ::Model::ValidSpatialDimensions [list 2D 3D]
    
    # Allow to open the tree
    set ::spdAux::TreeVisibility 1
    
    dict set attributes UseIntervals 0
    if {$::Kratos::kratos_private(DevMode) eq "dev"} {dict set attributes UseIntervals 1}
    
    LoadMyFiles
    ::spdAux::CreateDimensionWindow
}

proc ::EmbeddedFluid::LoadMyFiles { } {
    variable dir
    
    uplevel #0 [list source [file join $dir xml GetFromXML.tcl]]
    uplevel #0 [list source [file join $dir write write.tcl]]
    uplevel #0 [list source [file join $dir write writeProjectParameters.tcl]]
}

proc ::EmbeddedFluid::GetAttribute {name} {
    variable attributes
    set value ""
    catch {set value [dict get $attributes $name]}
    return $value
}

::EmbeddedFluid::Init

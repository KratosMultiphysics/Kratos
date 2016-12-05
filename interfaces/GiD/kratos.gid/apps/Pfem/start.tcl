namespace eval ::Pfem {
    # Variable declaration
    variable dir
    variable attributes
}

proc ::Pfem::Init { } {
    # Variable initialization
    variable dir
    variable attributes
    
    set dir [apps::getMyDir "Pfem"]
    set ::Model::ValidSpatialDimensions [list 2D 2Da 3D]
    # Allow to open the tree
    set ::spdAux::TreeVisibility 1
    set attributes [dict create]
    dict set attributes UseIntervals 1
    #if {$::Kratos::kratos_private(DevMode) eq "dev"} {dict set attributes UseIntervals 1}
    dict set attributes UseRestart 1
    LoadMyFiles
    ::spdAux::CreateDimensionWindow
}

proc ::Pfem::LoadMyFiles { } {
    variable dir
    uplevel #0 [list source [file join $dir xml GetFromXML.tcl]]
    uplevel #0 [list source [file join $dir xml ConditionSorterWindow.tcl]]
    uplevel #0 [list source [file join $dir .. Solid write write.tcl]]
    uplevel #0 [list source [file join $dir write write.tcl]]
    uplevel #0 [list source [file join $dir write writeProjectParameters.tcl]]
}


proc ::Pfem::GetAttribute {name} {
    variable attributes
    set value ""
    catch {set value [dict get $attributes $name]}
    return $value
}

proc ::Pfem::CustomToolbarItems { } {
    Kratos::ToolbarAddItem "Conditions" "list.png" [list -np- Pfem::xml::ConditionSorterWindow] [= "Sort the conditions"]   
}

::Pfem::Init

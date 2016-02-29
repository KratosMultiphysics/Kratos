namespace eval ::Fluid {
    # Variable declaration
    variable dir
}

proc ::Fluid::Init { } {
    # Variable initialization
    variable dir
    
    set dir [apps::getMyDir "Fluid"]
}

proc ::Fluid::LoadMyFiles { } {
    variable dir
    
    uplevel #0 [list source [file join $dir xml GetFromXML.tcl]]
    uplevel #0 [list source [file join $dir write write.tcl]]
}

proc ::Fluid::DrawAutoGeom { } {
    # Draw
    GiD_Process Mescape Geometry Create Object Prism 4 1 1 0 0.0,0.0,1.0 1 3 escape
    # Groups
    GiD_Process 'Groups Create Group0 escape 'Groups Edit Rename Group0 Body escape escape escape escape Utilities EntitiesGroups Assign Body Volumes 1 escape Mescape 'Groups Create Group0 escape 'Groups Edit Rename Group0 Inlet escape escape escape escape Utilities EntitiesGroups Assign Inlet Surfaces 6 escape Mescape 'Groups Create Group0 escape 'Groups Edit Rename Group0 Outlet escape escape escape escape Utilities EntitiesGroups Assign Outlet Surfaces 1 escape Mescape 'Groups Create Group0 escape 'Groups Edit Rename Group0 Slip escape escape escape escape Utilities EntitiesGroups Assign Slip Surfaces 2 2 3 5 2 escape Mescape 'Groups Create Group0 escape 'Groups Edit Rename Group0 NoSlip escape escape escape escape Utilities EntitiesGroups Assign NoSlip Surfaces 4 escape Mescape
}

::Fluid::Init

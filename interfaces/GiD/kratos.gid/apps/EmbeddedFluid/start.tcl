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
    
    apps::LoadAppById "Fluid"
    
    set dir [apps::getMyDir "EmbeddedFluid"]
    set attributes [dict create]
    
    set prefix EMBFL
    
    set ::Model::ValidSpatialDimensions [list 3D]
    spdAux::SetSpatialDimmension "3D"
    
    # Allow to open the tree
    set ::spdAux::TreeVisibility 1
    
    dict set attributes UseIntervals 0
    if {$::Kratos::kratos_private(DevMode) eq "dev"} {dict set attributes UseIntervals 1}
    
    LoadMyFiles
    #::spdAux::CreateDimensionWindow
}

proc ::EmbeddedFluid::LoadMyFiles { } {
    variable dir
    
    uplevel #0 [list source [file join $dir xml GetFromXML.tcl]]
    uplevel #0 [list source [file join $dir xml ImportWindowController.tcl]]
    uplevel #0 [list source [file join $dir write write.tcl]]
    uplevel #0 [list source [file join $dir write writeProjectParameters.tcl]]
}

proc ::EmbeddedFluid::GetAttribute {name} {
    variable attributes
    set value ""
    catch {set value [dict get $attributes $name]}
    return $value
}

proc ::EmbeddedFluid::CustomToolbarItems { } {
    Kratos::ToolbarAddItem "ImportMesh" "Import.png" [list -np- EmbeddedFluid::xml::ImportMeshWindow] [= "Import embedded mesh"]   
}

proc ::EmbeddedFluid::BeforeMeshGeneration {elementsize} {
    catch {file delete -force [file join $::write::dir "[file tail [GiD_Info project modelname] ].post.res"]}
    GiD_Process Escape Escape Utilities Variables EmbeddedMesh Activated 1 escape escape 

}

proc ::EmbeddedFluid::AfterMeshGeneration {fail} {
    GiD_Process Escape Escape Utilities Variables EmbeddedMesh Activated 0 escape escape 
}

::EmbeddedFluid::Init

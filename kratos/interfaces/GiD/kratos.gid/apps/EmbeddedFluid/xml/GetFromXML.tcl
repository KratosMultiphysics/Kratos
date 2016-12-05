namespace eval EmbeddedFluid::xml {
    # Namespace variables declaration
    variable dir
    variable lastImportMeshSize
    variable export_dir

}

proc EmbeddedFluid::xml::Init { } {
    # Namespace variables inicialization
    variable dir
    variable lastImportMeshSize
    set lastImportMeshSize 0
    #Model::DestroyEverything
    Model::ForgetMaterials
    Model::ForgetConstitutiveLaws
    Model::InitVariables dir $EmbeddedFluid::dir
    
    #Model::getSolutionStrategies Strategies.xml
    #Model::getElements Elements.xml
    Model::getMaterials Materials.xml
    #Model::getNodalConditions NodalConditions.xml
    Model::getConstitutiveLaws ConstitutiveLaws.xml
    #Model::getProcesses Processes.xml
    #Model::getConditions Conditions.xml
    #Model::getSolvers "../../Common/xml/Solvers.xml"
}


proc EmbeddedFluid::xml::MultiAppEvent {args} {
   if {$args eq "init"} {
     spdAux::parseRoutes
     spdAux::ConvertAllUniqueNames FL ${::EmbeddedFluid::prefix}
   }
}

proc EmbeddedFluid::xml::getUniqueName {name} {
    return ${::EmbeddedFluid::prefix}${name}
}

proc EmbeddedFluid::xml::CustomTree { args } {
    # Hide Results Cut planes
    spdAux::SetValueOnTreeItem v time Results FileLabel
    spdAux::SetValueOnTreeItem v time Results OutputControlType
    
    # Erase when Fractional step is available
    spdAux::SetValueOnTreeItem v Monolithic EMBFLSolStrat
    spdAux::SetValueOnTreeItem values Monolithic EMBFLSolStrat
    spdAux::SetValueOnTreeItem dict "Monolithic,Navier Stokes - Monolithic" EMBFLSolStrat
    spdAux::SetValueOnTreeItem v MN EMBFLScheme
    spdAux::SetValueOnTreeItem values MN EMBFLScheme
    spdAux::SetValueOnTreeItem dict "MN,Monolitic generic scheme" EMBFLScheme
    
}

proc EmbeddedFluid::xml::ImportMeshWindow { } {
    set filenames [Browser-ramR file read .gid [_ "Read mesh file"] \
		{} [list [list {STL Mesh} {.stl }] [list [_ "All files"] {.*}]] {} 1 \
        [list [_ "Import options"] EmbeddedFluid::xml::MoreImportOptions]]
    if {[llength $filenames] == 0} { return "" }
    set model_name [file rootname [file tail $filenames]]
    if {[GiD_Layers exists $model_name]} {
        set i 0
        set orig $model_name
        while {[GiD_Layers exists $model_name]} {
            set model_name ${orig}$i
            incr i
        }
    }
    GidUtils::DisableGraphics
    GiD_Layers create $model_name
    GiD_Layers edit to_use $model_name
    GiD_Process Mescape Files STLRead $filenames
    if {![GiD_Groups exists $model_name]} {GiD_Groups create $model_name}
    
    GiD_Process Mescape Geometry Edit ReConstruction OneSurfForEachElement Layer:$model_name escape 
    GiD_EntitiesGroups assign $model_name surfaces [GiD_EntitiesLayers get $model_name surfaces]
    
    GidUtils::EnableGraphics
    GidUtils::UpdateWindow GROUPS
    GiD_Process 'Zoom Frame
    
    GidUtils::SetWarnLine [= "STL %s imported!" $model_name]
}

proc EmbeddedFluid::xml::MoreImportOptions { f } {
    #variable _merge 
    ttk::label $f.lblGeometry -text [= "Mesh size"]:
    ttk::entry $f.entGeometry -textvariable EmbeddedFluid::xml::lastImportMeshSize  
    grid columnconfigure $f 1 -weight 1
    grid $f.lblGeometry $f.entGeometry -sticky w 

    variable export_dir
    if { ![info exists export_dir] } {
        #trick to show the more options
        set export_dir open
    }
    return EmbeddedFluid::xml::export_dir
}

EmbeddedFluid::xml::Init

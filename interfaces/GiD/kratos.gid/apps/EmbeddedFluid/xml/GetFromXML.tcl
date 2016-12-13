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
    Model::DestroyEverything
    Model::InitVariables dir $EmbeddedFluid::dir

    Model::getSolutionStrategies Strategies.xml
    Model::getElements "../../Fluid/xml/Elements.xml"
    Model::getMaterials Materials.xml
    Model::getNodalConditions "../../Fluid/xml/NodalConditions.xml"
    Model::getConstitutiveLaws ConstitutiveLaws.xml
    Model::getProcesses "../../Fluid/xml/Processes.xml"
    Model::getConditions "../../Fluid/xml/Conditions.xml"
    Model::getSolvers "../../Common/xml/Solvers.xml"
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
    
    set root [customlib::GetBaseRoot]
    if {[$root selectNodes $xp] eq "[spdAux::getRoute NodalResults]/value\[@n='DISTANCE'\]"} {gid_groups_conds::addF [spdAux::getRoute NodalResults] value [list n DISTANCE pn Distance v Yes values {Yes,No} state normal]}
    if {[$root selectNodes $xp] eq [spdAux::getRoute EMBFLSolutionParameters]} {gid_groups_conds::addF [spdAux::getRoute EMBFLSolutionParameters] include [list n DistanceReading active 1 path {apps/EmbeddedFluid/xml/DistanceReading.spd}]}
    customlib::ProcessIncludes $::Kratos::kratos_private(Path)
    spdAux::parseRoutes
    # Erase when Fractional step is available
    #spdAux::SetValueOnTreeItem v Monolithic EMBFLSolStrat
    #spdAux::SetValueOnTreeItem values Monolithic EMBFLSolStrat
    #spdAux::SetValueOnTreeItem dict "Monolithic,Navier Stokes - Monolithic" EMBFLSolStrat
    #spdAux::SetValueOnTreeItem v MN EMBFLScheme
    #spdAux::SetValueOnTreeItem values MN EMBFLScheme
    #spdAux::SetValueOnTreeItem dict "MN,Monolitic generic scheme" EMBFLScheme

}

EmbeddedFluid::xml::Init

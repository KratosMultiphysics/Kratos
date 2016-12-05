namespace eval EmbeddedFluid::xml {
    # Namespace variables declaration
    variable dir
}

proc EmbeddedFluid::xml::Init { } {
    # Namespace variables inicialization
    variable dir
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

EmbeddedFluid::xml::Init

namespace eval Structural::xml {
     variable dir
}

proc Structural::xml::Init { } {
     variable dir
     Model::InitVariables dir $Structural::dir
    
     Model::getSolutionStrategies strategydefinition.xml
     Model::getProcesses Processes.xml
     Model::getElements Elements.xml
     Model::getConstitutiveLaws ConstitutiveLaws.xml
     Model::getConditions Conditions.xml
     Model::getSolvers Solvers.xml
     
     # As we do load all the Solid app in the Start.tcl file apps::LoadAppById "Solid" we do not need to load it's files
     #Model::getSolutionStrategies  "../../Solid/xml/strategydefinition.xml"
     #Model::getProcesses           "../../Solid/xml/Processes.xml"
     #Model::getElements            "../../Solid/xml/Elements.xml"
     #Model::getConstitutiveLaws    "../../Solid/xml/ConstitutiveLaws.xml"
     #Model::getConditions          "../../Solid/xml/Conditions.xml"
     #Model::getSolvers             "../../Solid/xml/Solvers.xml"
     
}

proc Structural::xml::getUniqueName {name} {
    return ST$name
}

proc ::Structural::xml::MultiAppEvent {args} {
   if {$args eq "init"} {
     spdAux::parseRoutes
     spdAux::ConvertAllUniqueNames SL ST
   }
}

Structural::xml::Init

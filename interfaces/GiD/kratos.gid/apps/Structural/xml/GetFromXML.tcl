namespace eval Structural::xml {
     variable dir
}

proc Structural::xml::Init { } {
     variable dir
     Model::InitVariables dir $Structural::dir
    
     Model::ForgetSolutionStrategies
     Model::getSolutionStrategies Strategies.xml
     Model::getProcesses Processes.xml
     Model::getElements Elements.xml
     Model::getNodalConditions NodalConditions.xml
     Model::getConstitutiveLaws ConstitutiveLaws.xml
     Model::getConditions Conditions.xml
     Model::getSolvers Solvers.xml
     
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

proc Structural::xml::CustomTree { args } {
    Solid::xml::CustomTree $args
}

Structural::xml::Init

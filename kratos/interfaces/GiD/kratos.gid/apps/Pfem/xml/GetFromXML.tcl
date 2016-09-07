namespace eval Pfem::xml {
     variable dir
}

proc Pfem::xml::Init { } {
    variable dir
    Model::InitVariables dir $Pfem::dir
    
    Model::getSolutionStrategies Strategies.xml
    Model::getElements Elements.xml
    Model::getConstitutiveLaws "../../Solid/xml/ConstitutiveLaws.xml"
    Model::getConstitutiveLaws "../../Fluid/xml/ConstitutiveLaws.xml"
    Model::getProcesses "../../Solid/xml/Processes.xml"
    Model::getNodalConditions "../../Solid/xml/NodalConditions.xml"
    Model::getConditions "../../Solid/xml/Conditions.xml"
    Model::getSolvers "../../Common/xml/Solvers.xml"

}

proc Pfem::xml::getUniqueName {name} {
    return PFEM_$name
}


proc ::Pfem::xml::MultiAppEvent {args} {
   if {$args eq "init"} {
     spdAux::parseRoutes
     spdAux::ConvertAllUniqueNames SL PFEM_
   }
}

Pfem::xml::Init

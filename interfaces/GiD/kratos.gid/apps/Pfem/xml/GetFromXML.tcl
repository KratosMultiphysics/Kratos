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

proc ::Pfem::xml::::CheckElementOutputState { domNode args } {
     set elemsactive [list ]
     foreach parts_un [Pfem::write::GetPartsUN] {
          set parts_path [spdAux::getRoute $parts_un]
          set xp1 "$parts_path/group/value\[@n='Element'\]"
          foreach gNode [[customlib::GetBaseRoot] selectNodes $xp1] {
              lappend elemsactive [get_domnode_attribute $gNode v]
          }
     }
     W $elemsactive
     set paramName [$domNode @n]
     return [::Model::CheckElementOutputState $elemsactive $paramName]
}

Pfem::xml::Init

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
     set paramName [$domNode @n]
     return [::Model::CheckElementOutputState $elemsactive $paramName]
}

proc Pfem::xml::ProcGetElements {domNode args} {
     set nodeApp [spdAux::GetAppIdFromNode $domNode]
     set sol_stratUN [apps::getAppUniqueName $nodeApp SolStrat]
     set schemeUN [apps::getAppUniqueName $nodeApp Scheme]
     if {[get_domnode_attribute [$domNode selectNodes [spdAux::getRoute $sol_stratUN]] v] eq ""} {
         get_domnode_attribute [$domNode selectNodes [spdAux::getRoute $sol_stratUN]] dict
     }
     if {[get_domnode_attribute [$domNode selectNodes [spdAux::getRoute $schemeUN]] v] eq ""} {
          get_domnode_attribute [$domNode selectNodes [spdAux::getRoute $schemeUN]] dict
     }
     
     #W "solStrat $sol_stratUN sch $schemeUN"
     set solStratName [::write::getValue $sol_stratUN]
     set schemeName [write::getValue $schemeUN]
     #W "$solStratName $schemeName"
     #W "************************************************************************"
     #W "$nodeApp $solStratName $schemeName"
     set elems [::Model::GetAvailableElements $solStratName $schemeName]
     #W "************************************************************************"
     set names [list ]
     set pnames [list ]
     set BodyType [get_domnode_attribute [[[$domNode parent] parent] selectNodes "value\[@n='BodyType'\]"] v]
     set argums [list ElementType $BodyType]
     
     foreach elem $elems {
         if {[$elem cumple $argums]} {
             lappend names [$elem getName]
             lappend pnames [$elem getName] 
             lappend pnames [$elem getPublicName]
         }
     }
     set diction [join $pnames ","]
     set values [join $names ","]
     
     $domNode setAttribute values $values
     if {[get_domnode_attribute $domNode v] eq ""} {$domNode setAttribute v [lindex $names 0]}
     if {[get_domnode_attribute $domNode v] ni $names} {$domNode setAttribute v [lindex $names 0]}
     #spdAux::RequestRefresh
     return $diction
}

Pfem::xml::Init

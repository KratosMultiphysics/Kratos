namespace eval StenosisWizard::xml {
     variable dir
}

proc StenosisWizard::xml::Init { } {
     variable dir
     Model::InitVariables dir $StenosisWizard::dir
}

proc StenosisWizard::xml::getUniqueName {name} {
    return StenWiz$name
}

proc ::StenosisWizard::xml::MultiAppEvent {args} {
   if {$args eq "init"} {
     spdAux::parseRoutes
     spdAux::ConvertAllUniqueNames FL StenWiz
   }
}

StenosisWizard::xml::Init

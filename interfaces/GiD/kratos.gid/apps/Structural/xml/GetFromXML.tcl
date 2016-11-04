namespace eval Structural::xml {
     variable dir
}

proc Structural::xml::Init { } {
     variable dir
     Model::InitVariables dir $Structural::dir
    
     Model::ForgetSolutionStrategies
     Model::getSolutionStrategies Strategies.xml
     Model::getElements Elements.xml
     Model::getNodalConditions NodalConditions.xml
     Model::getConstitutiveLaws ConstitutiveLaws.xml
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

proc Structural::xml::ProcCheckGeometryStructural {domNode args} {
     set ret "line,surface"
     if {$::Model::SpatialDimension eq "3D"} {
          set ret "line,surface,volume"
     }
     return $ret
}

Structural::xml::Init

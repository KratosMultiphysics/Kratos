# Project Parameters
proc FSI::write::getParametersDict { } {
   set projectParametersDict [dict create]
   UpdateUniqueNames Structural
   apps::setActiveAppSoft Structural
   set StructuralParametersDict [Structural::write::getParametersEvent]
   UpdateUniqueNames Fluid
   apps::setActiveAppSoft Fluid
   set FluidParametersDict [Fluid::write::getParametersDict]
   UpdateUniqueNames FSI
   apps::setActiveAppSoft FSI
   set FSIParametersDict [dict create]
   
   
   
   dict set projectParametersDict structure_solver_settings $StructuralParametersDict
   dict set projectParametersDict fluid_solver_settings $FluidParametersDict
   dict set projectParametersDict coupling_solver_settings $FSIParametersDict
   return $projectParametersDict
}

proc FSI::write::writeParametersEvent { } {
    set projectParametersDict [getParametersDict]
    write::WriteJSON $projectParametersDict
}

proc FSI::write::UpdateUniqueNames { appid } {
    set unList [list "Results"]
    foreach un $unList {
         set current_un [apps::getAppUniqueName $appid $un]
         spdAux::setRoute $un [spdAux::getRoute $current_un]
    }
}

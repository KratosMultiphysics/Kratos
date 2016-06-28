# Project Parameters
proc FSI::write::getParametersDict { } {
   set projectParametersDict [dict create]
   
   # FSI section
   set FSIParametersDict [dict create]
   # Solver settings
   set solverSettingsDict [dict create]
   set currentStrategyId [write::getValue FSISolStrat]
   set strategy_write_name [[::Model::GetSolutionStrategy $currentStrategyId] getAttribute "ImplementedInPythonFile"]
   dict set solverSettingsDict strategy_id $currentStrategyId
   dict set solverSettingsDict solver_type $strategy_write_name
   set solverSettingsDict [dict merge $solverSettingsDict [write::getSolutionStrategyParametersDict] ]
   set solverSettingsDict [dict merge $solverSettingsDict [write::getSolversParametersDict FSI] ]
   dict set FSIParametersDict solver_settings $solverSettingsDict
   
   # Structural section
   UpdateUniqueNames Structural
   apps::setActiveAppSoft Structural
   
   set StructuralParametersDict [Structural::write::getParametersEvent]
   set current [dict get $StructuralParametersDict solver_settings model_import_settings input_filename]
   dict set StructuralParametersDict solver_settings model_import_settings input_filename "${current}_Structural"
   
   # Fluid section
   UpdateUniqueNames Fluid
   apps::setActiveAppSoft Fluid
   
   set FluidParametersDict [Fluid::write::getParametersDict]
   set current [dict get $FluidParametersDict solver_settings model_import_settings input_filename]
   dict set FluidParametersDict solver_settings model_import_settings input_filename "${current}_Fluid"
   UpdateUniqueNames FSI
   apps::setActiveAppSoft FSI
   
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


proc PotentialFluid::write::writeParametersEvent { } {
    set projectParametersDict [Fluid::write::getParametersDict]
    
    set projectParametersDict [dict remove $projectParametersDict gravity]
    set projectParametersDict [dict remove $projectParametersDict auxiliar_process_list]
    dict set projectParametersDict solver_settings [dict remove [dict get $projectParametersDict solver_settings] time_stepping]
    
    write::SetParallelismConfiguration
    write::WriteJSON $projectParametersDict
}

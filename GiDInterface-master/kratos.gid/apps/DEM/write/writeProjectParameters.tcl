
# Project Parameters
proc DEM::write::getParametersEvent { } {
    set project_parameters_dict [dict create]
    return $project_parameters_dict
}
proc DEM::write::writeParametersEvent { } {
    write::WriteJSON [getParametersEvent]
}

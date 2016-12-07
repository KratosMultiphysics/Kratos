# Project Parameters
proc ::EmbeddedFluid::write::getParametersDict { } {
      return [Fluid::write::getParametersDict]
}

proc EmbeddedFluid::write::writeParametersEvent { } {
    set projectParametersDict [getParametersDict]
    write::SetParallelismConfiguration
    write::WriteJSON $projectParametersDict
}

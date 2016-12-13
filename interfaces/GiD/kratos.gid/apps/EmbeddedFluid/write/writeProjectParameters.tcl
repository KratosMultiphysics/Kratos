# Project Parameters
proc ::EmbeddedFluid::write::getParametersDict { } {
      set param_dict [Fluid::write::getParametersDict]
      set bclist [dict get $param_dict boundary_conditions_process_list]
      lappend bclist [GetDistanceModificationDict]
      dict set param_dict boundary_conditions_process_list $bclist
      return $param_dict
}

proc EmbeddedFluid::write::writeParametersEvent { } {
    set projectParametersDict [getParametersDict]
    write::SetParallelismConfiguration
    write::WriteJSON $projectParametersDict
}

proc EmbeddedFluid::write::GetDistanceModificationDict { } {
      set distance_modif_dict [dict create ]
      dict set distance_modif_dict python_module apply_distance_modification_process
      dict set distance_modif_dict kratos_module KratosMultiphysics.FluidDynamicsApplication
      dict set distance_modif_dict process_name ApplyDistanceModificationProcess
            set parameters_dict [dict create ]
            dict set parameters_dict mesh_id 0
            dict set parameters_dict model_part_name [lindex [write::getPartsMeshId] 0]
            dict set parameters_dict check_at_each_time_step false
      dict set distance_modif_dict Parameters $parameters_dict
      return $distance_modif_dict
}

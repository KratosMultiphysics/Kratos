# Project Parameters
proc ::EmbeddedFluid::write::getParametersDict { } {
      set param_dict [Fluid::write::getParametersDict]

      # Set the auxiliar embedded fluid application processes dictionary list
      dict set param_dict auxiliar_process_list [getAuxiliarProcessList]

      # Set the solver settings dictionary
      set solverSettingsDict [dict get $param_dict solver_settings]
      set solverSettingsDict [dict merge $solverSettingsDict [write::getSolversParametersDict EmbeddedFluid] ]

      # Set the distance reading settings dictionary
      set dist_settings_dict [dict create]
      set dist_mode [write::getValue EMBFLDistanceSettings ReadingMode]
      dict set dist_settings_dict import_mode $dist_mode
      if {$dist_mode ne "from_mdpa"} {
            set dist_file [write::getValue EMBFLDistanceSettings distance_file_name]
            dict set dist_settings_dict distance_file_name $dist_file
      }
      dict set solverSettingsDict distance_reading_settings $dist_settings_dict

      dict set param_dict solver_settings $solverSettingsDict
      return $param_dict
}

proc EmbeddedFluid::write::writeParametersEvent { } {
    set projectParametersDict [getParametersDict]
    write::SetParallelismConfiguration
    write::WriteJSON $projectParametersDict
}

proc EmbeddedFluid::write::getAuxiliarProcessList {} {
    set auxiliar_process_list [list]

    # Append the distance modification process
    lappend auxiliar_process_list [getDistanceModificationDict]
    # Append the embedded drag process dictionary
    lappend auxiliar_process_list [getEmbeddedDragProcessDict]

    return $auxiliar_process_list
}

proc EmbeddedFluid::write::getDistanceModificationDict { } {
      set distance_modif_dict [dict create ]
      dict set distance_modif_dict "python_module" apply_distance_modification_process
      dict set distance_modif_dict "kratos_module" KratosMultiphysics.FluidDynamicsApplication
      dict set distance_modif_dict "process_name" ApplyDistanceModificationProcess
            set parameters_dict [dict create ]
            dict set parameters_dict "mesh_id" 0
            dict set parameters_dict "model_part_name" [lindex [write::getPartsMeshId] 0]
            dict set parameters_dict "check_at_each_time_step" [write::getValue EMBFLDistanceSettings correct_distance_at_each_step]
      dict set distance_modif_dict "Parameters" $parameters_dict
      return $distance_modif_dict
}

proc EmbeddedFluid::write::getEmbeddedDragProcessDict {} {
    set pdict [dict create]
    dict set pdict "python_module" "compute_embedded_drag_process"
    dict set pdict "kratos_module" "KratosMultiphysics.FluidDynamicsApplication"
    dict set pdict "process_name" "ComputeEmbeddedDragProcess"
        set params [dict create]
        dict set params "mesh_id" 0
        dict set params "fluid_model_part" [lindex [write::getPartsMeshId] 0]
        dict set params "write_drag_output_file" [write::getValue EMBFLEmbeddedDrag write_drag_output_file]
        dict set params "print_drag_to_screen" [write::getValue EMBFLEmbeddedDrag print_drag_to_screen]
        dict set params "interval" [write::getInterval [write::getValue EMBFLEmbeddedDrag Interval] ]
    dict set pdict "Parameters" $params

    return $pdict
}

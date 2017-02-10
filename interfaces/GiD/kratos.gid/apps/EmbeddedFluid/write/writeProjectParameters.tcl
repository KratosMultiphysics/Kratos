# Project Parameters
proc ::EmbeddedFluid::write::getParametersDict { } {
      set param_dict [Fluid::write::getParametersDict]

      ## Set the meshing adaptivity settings
      set mesh_adaptivity [write::getValue EMBFLAdaptivitySettings mesh_adaptivity]
      # Meshing adaptivity switch on/off in problem data dict
      set new_problem_data [dict get $param_dict problem_data]
      dict set new_problem_data "mesh_adaptivity" $mesh_adaptivity
      dict set param_dict problem_data $new_problem_data
      # Set the meshing adaptivity process list
      if {$mesh_adaptivity eq "Yes"} {
          dict set param_dict mesh_adaptivity_process_list [list [getMeshAdaptivityProcessDict]]
      } else {
          dict set param_dict mesh_adaptivity_process_list [list]
      }

      ## Set the auxiliar embedded fluid application processes dictionary list
      dict set param_dict auxiliar_process_list [getAuxiliarProcessList]

      ## Set the solver settings dictionary
      set solverSettingsDict [dict get $param_dict solver_settings]
      # If drag has to be computed, ensure that "compute_reactions" is set to true
      set compute_embedded_drag [write::getValue EMBFLEmbeddedDrag compute_embedded_drag]
      if {$compute_embedded_drag eq "Yes"} {
          dict set solverSettingsDict "compute_reactions" true
      }
      set solverSettingsDict [dict merge $solverSettingsDict [write::getSolversParametersDict EmbeddedFluid] ]

      ## Set the distance reading settings dictionary
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

    # If required, append the embedded drag process dictionary
    set compute_embedded_drag [write::getValue EMBFLEmbeddedDrag compute_embedded_drag]
    if {$compute_embedded_drag eq "Yes"} {
        lappend auxiliar_process_list [getEmbeddedDragProcessDict]
    }

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
        dict set params "model_part_name" [lindex [write::getPartsMeshId] 0]
        dict set params "write_drag_output_file" [write::getValue EMBFLEmbeddedDrag write_drag_output_file]
        dict set params "print_drag_to_screen" [write::getValue EMBFLEmbeddedDrag print_drag_to_screen]
        dict set params "interval" [write::getInterval [write::getValue EMBFLEmbeddedDrag Interval] ]
    dict set pdict "Parameters" $params

    return $pdict
}

proc EmbeddedFluid::write::getMeshAdaptivityProcessDict {} {
    set pdict [dict create]
    dict set pdict "python_module" "mmg_process"
    dict set pdict "kratos_module" "KratosMultiphysics.MeshingApplication"
    dict set pdict "process_name" "MmgProcess"
        # MmgProcess settings dictionary
        set params_dict [dict create]
        dict set params_dict "mesh_id" 0
        dict set params_dict "model_part_name" "MainModelPart"
        dict set params_dict "initial_step" [write::getValue EMBFLAdaptivitySettings initial_step]
        dict set params_dict "step_frequency" [write::getValue EMBFLAdaptivitySettings step_frequency]
        dict set params_dict "initial_remeshing" [write::getValue EMBFLAdaptivitySettings initial_remeshing]
            # Set anisotropic parameters section
            set anisotropy_parameters_dict [dict create]
            dict set anisotropy_parameters_dict "hmin_over_hmax_anisotropic_ratio" [write::getValue EMBFLAdaptivitySettings hmin_over_hmax_anisotropic_ratio]
            dict set anisotropy_parameters_dict "boundary_layer_min_size_ratio" [write::getValue EMBFLAdaptivitySettings boundary_layer_min_size_ratio]
        dict set params_dict "anisotropy_parameters" $anisotropy_parameters_dict
        dict set params_dict "echo_level" [write::getValue EMBFLAdaptivitySettings echo_level]
    dict set pdict "Parameters" $params_dict

    return $pdict
}

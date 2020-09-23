proc WriteProjectParameters { basename dir problemtypedir } {

    ## Source auxiliar procedures
    source [file join $problemtypedir ProjectParametersAuxProcs.tcl]


    ## Start ProjectParameters.json file
    set filename [file join $dir ProjectParameters.json]
    set FileVar [open $filename w]

    puts $FileVar "\{"


    ## problem_data
    puts $FileVar "    \"problem_data\"             : \{"
    puts $FileVar "        \"problem_name\"             : \"$basename\","
    puts $FileVar "        \"echo_level\"               : [GiD_AccessValue get gendata Echo_Level],"
    puts $FileVar "        \"start_time\"               : [GiD_AccessValue get gendata Start_Time],"
    puts $FileVar "        \"end_time\"                 : [GiD_AccessValue get gendata End_Time],"
    puts $FileVar "        \"parallel_type\"            : \"[GiD_AccessValue get gendata Parallel_Configuration]\""
    puts $FileVar "    \},"


    ## solver_settings
    puts $FileVar "    \"solver_settings\"            : \{"
    puts $FileVar "        \"solver_type\"                : \"stabilized_shallow_water_solver\","
    puts $FileVar "        \"model_part_name\"            : \"model_part\","
    puts $FileVar "        \"domain_size\"                : 2,"
    puts $FileVar "        \"gravity\"                    : [GiD_AccessValue get gendata Gravity],"
    puts $FileVar "        \"model_import_settings\"      : \{"
    puts $FileVar "            \"input_type\"                 : \"mdpa\","
    puts $FileVar "            \"input_filename\"             : \"$basename\""
    puts $FileVar "        \},"
    puts $FileVar "        \"echo_level\"                 : [GiD_AccessValue get gendata Echo_Level],"
    puts $FileVar "        \"buffer_size\"                : 2,"
    puts $FileVar "        \"stabilization_factor\"       : [GiD_AccessValue get gendata Stabilization_parameter],"
    puts $FileVar "        \"shock_stabilization_factor\" : [GiD_AccessValue get gendata Shock_stabilization_parameter],"
    puts $FileVar "        \"relative_tolerance\"         : [GiD_AccessValue get gendata Relative_tolerance],"
    puts $FileVar "        \"absolute_tolerance\"         : [GiD_AccessValue get gendata Absolute_tolerance],"
    puts $FileVar "        \"maximum_iterations\"         : [GiD_AccessValue get gendata Maximum_iterations],"
    puts $FileVar "        \"compute_reactions\"          : [GiD_AccessValue get gendata Compute_Reactions],"
    puts $FileVar "        \"reform_dofs_at_each_step\"   : [GiD_AccessValue get gendata Reform_Dofs_At_Each_Step],"
    puts $FileVar "        \"move_mesh_flag\"             : [GiD_AccessValue get gendata Move_Mesh],"
    # linear_solver_settings
    puts $FileVar "        \"linear_solver_settings\"   : \{"
        puts $FileVar "            \"solver_type\"      : \"[GiD_AccessValue get gendata Solver_Type]\","
    puts $FileVar "        \},"
    puts $FileVar "        \"time_stepping\"            : \{"
    if {[GiD_AccessValue get gendata Time_stepping] eq "Automatic"} {
        puts $FileVar "            \"automatic_time_step\"      : true,"
        puts $FileVar "            \"courant_number\"           : [GiD_AccessValue get gendata Courant_number]"
    } else {
        puts $FileVar "            \"automatic_time_step\"      : false,"
        puts $FileVar "            \"time_step\"                : [GiD_AccessValue get gendata Time_step]"
    }
    puts $FileVar "        \}"
    puts $FileVar "    \},"


    ## output processes
    puts $FileVar "    \"output_processes\" : \{"
    puts $FileVar "        \"output_process_list\" : \[\{"
    puts $FileVar "            \"kratos_module\"        : \"KratosMultiphysics\","
    puts $FileVar "            \"python_module\"        : \"gid_output_process\","
    puts $FileVar "            \"process_name\"         : \"GidOutputProcess\","
    puts $FileVar "            \"Parameters\"           : \{"
    puts $FileVar "                \"model_part_name\"        : \"model_part\","
    puts $FileVar "                \"output_name\"            : \"$basename\","
    puts $FileVar "                \"postprocess_parameters\" : \{"
    puts $FileVar "                    \"result_file_configuration\" : \{"
    puts $FileVar "                        \"gidpost_flags\"         : \{"
    puts $FileVar "                            \"GiDPostMode\"           : \"[GiD_AccessValue get gendata GiD_post_mode]\","
    puts $FileVar "                            \"WriteDeformedMeshFlag\" : \"[GiD_AccessValue get gendata Write_deformed_mesh]\","
    puts $FileVar "                            \"WriteConditionsFlag\"   : \"[GiD_AccessValue get gendata Write_conditions]\","
    puts $FileVar "                            \"MultiFileFlag\"         : \"[GiD_AccessValue get gendata Multi_file_flag]\""
    puts $FileVar "                        \},"
    puts $FileVar "                        \"output_control_type\"   : \"[GiD_AccessValue get gendata Output_control_type]\","
    puts $FileVar "                        \"output_interval\"       : [GiD_AccessValue get gendata Output_frequency],"
    puts $FileVar "                        \"body_output\"           : [GiD_AccessValue get gendata Body_output],"
    puts $FileVar "                        \"node_output\"           : [GiD_AccessValue get gendata Node_output],"
    puts $FileVar "                        \"skin_output\"           : [GiD_AccessValue get gendata Skin_output],"
    puts $FileVar "                        \"plane_output\"          : \[\],"
    # nodal_results
    set PutStrings \[
    append PutStrings \"MOMENTUM\",\"VELOCITY\",\"HEIGHT\",\"FREE_SURFACE_ELEVATION\",\"TOPOGRAPHY\"
    set PutStrings [string trimright $PutStrings ,]
    append PutStrings \]
    puts $FileVar "                        \"nodal_results\"         : $PutStrings,"
    puts $FileVar "                        \"gauss_point_results\"   : \[\]"
    puts $FileVar "                    \},"
    puts $FileVar "                \}"
    puts $FileVar "            \}"
    puts $FileVar "        \}\]"
    puts $FileVar "    \},"


    ## processes
    puts $FileVar "    \"processes\"    : \{"

    # initial conditions
    set Groups [GiD_Info conditions Initial_water_level groups]
    set NumGroups [llength $Groups]
    set iGroup 0
    puts $FileVar "        \"initial_conditions_process_list\"   : \[\{"
    WriteInitialWaterLevelProcess FileVar iGroup $Groups surfaces $NumGroups

    # boundary conditions
    set Groups [GiD_Info conditions Slip_condition groups]
    set NumGroups [llength $Groups]
    set Groups [GiD_Info conditions Water_height groups]
    incr NumGroups [llength $Groups]
    set Groups [GiD_Info conditions Imposed_flux groups]
    incr NumGroups [llength $Groups]
    set iGroup 0
    puts $FileVar "        \"boundary_conditions_process_list\"  : \[\{"
    ## Slip conditions
    set Groups [GiD_Info conditions Slip_condition groups]
    WriteSlipConditionProcess FileVar iGroup $Groups lines $NumGroups
    ## Imposed water height
    set Groups [GiD_Info conditions Water_height groups]
    WriteConstantScalarConditionProcess FileVar iGroup $Groups lines $NumGroups
    ## Imposed water flux
    set Groups [GiD_Info conditions Imposed_flux groups]
    WriteConstantVectorConditionProcess FileVar iGroup $Groups lines $NumGroups

    # bathymetry
    set Groups [GiD_Info conditions Body_Part groups]
    set NumGroups [llength $Groups]
    set iGroup 0
    puts $FileVar "        \"bathymetry_process_list\"     : \[\{"
    WriteBathymetryProcess FileVar iGroup $Groups surfaces $NumGroups

    # end of processes
    puts $FileVar "    \}"


    ## Finish ProjectParameters.json file
    puts $FileVar "\}"

    close $FileVar
}

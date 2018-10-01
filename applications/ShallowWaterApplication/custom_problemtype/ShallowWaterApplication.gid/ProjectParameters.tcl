proc WriteProjectParameters { basename dir problemtypedir } {

    ## Source auxiliar procedures
    source [file join $problemtypedir ProjectParametersAuxProcs.tcl]

    # Start ProjectParameters.json file
    set filename [file join $dir ProjectParameters.json]
    set FileVar [open $filename w]

    puts $FileVar "\{"


    # problem_data
    puts $FileVar "    \"problem_data\"             : \{"
    puts $FileVar "        \"problem_name\"             : \"$basename\","
    puts $FileVar "        \"domain_size\"              : 2,"
    puts $FileVar "        \"start_time\"               : [GiD_AccessValue get gendata Start_Time],"
    puts $FileVar "        \"end_time\"                 : [GiD_AccessValue get gendata End_Time],"
    puts $FileVar "        \"parallel_type\"            : \"[GiD_AccessValue get gendata Parallel_Configuration]\","
    puts $FileVar "        \"number_of_threads\"        : [GiD_AccessValue get gendata Number_of_threads]"
    puts $FileVar "    \},"


    # solver_settings
    puts $FileVar "    \"solver_settings\"          : \{"
    if {[GiD_AccessValue get gendata Framework] eq "Pfem2"} {
    if {[GiD_AccessValue get gendata Variables] eq "Primitive"} {
    puts $FileVar "        \"solver_type\"              : \"pfem2_primitive_var_solver\","
    } else {
    puts $FileVar "        \"solver_type\"              : \"pfem2_conserved_var_solver\","
    }
    } else {
    if {[GiD_AccessValue get gendata Variables] eq "Primitive"} {
    puts $FileVar "        \"solver_type\"              : \"eulerian_primitive_var_solver\","
    } else {
    puts $FileVar "        \"solver_type\"              : \"eulerian_conserved_var_solver\","
    }
    }
    puts $FileVar "        \"model_part_name\"          : \"main_model_part\","
    puts $FileVar "        \"domain_size\"              : 2,"
    puts $FileVar "        \"gravity\"                  : [GiD_AccessValue get gendata Gravity],"
    puts $FileVar "        \"time_scale\"               : \"seconds\","
    puts $FileVar "        \"water_height_scale\"       : \"meters\","
    puts $FileVar "        \"model_import_settings\"    : \{"
    puts $FileVar "            \"input_type\"               : \"mdpa\","
    puts $FileVar "            \"input_filename\"           : \"$basename\""
    puts $FileVar "        \},"
    puts $FileVar "        \"echo_level\"               : [GiD_AccessValue get gendata Echo_Level],"
    puts $FileVar "        \"buffer_size\"              : 2,"
    puts $FileVar "        \"dynamic_tau\"              : [GiD_AccessValue get gendata Stabilization_parameter],"
    puts $FileVar "        \"relative_tolerance\"       : [GiD_AccessValue get gendata Relative_tolerance],"
    puts $FileVar "        \"absolute_tolerance\"       : [GiD_AccessValue get gendata Absolute_tolerance],"
    puts $FileVar "        \"maximum_iterations\"       : [GiD_AccessValue get gendata Maximum_iterations],"
    puts $FileVar "        \"compute_reactions\"        : [GiD_AccessValue get gendata Compute_Reactions],"
    puts $FileVar "        \"reform_dofs_at_each_step\" : [GiD_AccessValue get gendata Reform_Dofs_At_Each_Step],"
    puts $FileVar "        \"move_mesh_flag\"           : [GiD_AccessValue get gendata Move_Mesh],"
    puts $FileVar "        \"volume_model_part_name\"   : \"shallow_water_model\","
    # linear_solver_settings
    puts $FileVar "        \"linear_solver_settings\"   : \{"
    if {[GiD_AccessValue get gendata Parallel_Configuration] eq "MPI"} {
    if {[GiD_AccessValue get gendata Solver_Type] eq "AmgclMPISolver"} {
        puts $FileVar "            \"solver_type\"      : \"AmgclMPISolver\","
        puts $FileVar "            \"krylov_type\"      : \"fgmres\","
        puts $FileVar "            \"max_iteration\"    : 100,"
        puts $FileVar "            \"tolerance\"        : 1.0e-6,"
    } elseif {[GiD_AccessValue get gendata Solver_Type] eq "AztecSolver"} {
        puts $FileVar "            \"solver_type\"      : \"AztecSolver\","
        puts $FileVar "            \"tolerance\"        : 1.0e-6,"
        puts $FileVar "            \"max_iteration\"    : 200,"
        puts $FileVar "            \"preconditioner_type\": \"None\""
    } elseif {([GiD_AccessValue get gendata Solver_Type] eq "Klu") || ([GiD_AccessValue get gendata Solver_Type] eq "MultiLevelSolver")} {
        puts $FileVar "            \"solver_type\"      : \"[GiD_AccessValue get gendata Solver_Type]\","
        puts $FileVar "            \"scaling\"          : false"
    } else {
        puts $FileVar "            \"solver_type\"      : \"Klu\","
        puts $FileVar "            \"scaling\"          : false"
    }
    } else {
    if {[GiD_AccessValue get gendata Solver_Type] eq "AMGCL"} {
        puts $FileVar "            \"solver_type\"      : \"AMGCL\","
        puts $FileVar "            \"smoother_type\"    : \"ilu0\","
        puts $FileVar "            \"krylov_type\"      : \"gmres\","
        puts $FileVar "            \"coarsening_type\"  : \"aggregation\","
        puts $FileVar "            \"max_iteration\"    : 100,"
        puts $FileVar "            \"tolerance\"        : 1.0e-6,"
        puts $FileVar "            \"scaling\"          : false"
    } elseif {[GiD_AccessValue get gendata Solver_Type] eq "BICGSTABSolver"} {
        puts $FileVar "            \"solver_type\"      : \"BICGSTABSolver\","
        puts $FileVar "            \"tolerance\"        : 1.0e-6,"
        puts $FileVar "            \"max_iteration\"    : 100,"
        puts $FileVar "            \"scaling\"          : false,"
        puts $FileVar "            \"preconditioner_type\": \"ILU0Preconditioner\""
    } elseif {([GiD_AccessValue get gendata Solver_Type] eq "SkylineLUFactorizationSolver") || ([GiD_AccessValue get gendata Solver_Type] eq "SuperLUSolver")} {
        puts $FileVar "            \"solver_type\"      : \"[GiD_AccessValue get gendata Solver_Type]\""
    } else {
        puts $FileVar "            \"solver_type\"      : \"SuperLUSolver\""
    }
    }
    puts $FileVar "        \},"
    puts $FileVar "        \"time_stepping\"            : \{"
    puts $FileVar "            \"automatic_time_step\"      : false,"
    puts $FileVar "            \"time_step\"                : [GiD_AccessValue get gendata Delta_Time]"
    puts $FileVar "        \}"
    puts $FileVar "    \},"


    # output_configuration
    puts $FileVar "    \"output_configuration\" : \{"
    puts $FileVar "        \"result_file_configuration\" : \{"
    puts $FileVar "            \"gidpost_flags\"         : \{"
    puts $FileVar "                \"GiDPostMode\"           : \"[GiD_AccessValue get gendata GiD_post_mode]\","
    puts $FileVar "                \"WriteDeformedMeshFlag\" : \"[GiD_AccessValue get gendata Write_deformed_mesh]\","
    puts $FileVar "                \"WriteConditionsFlag\"   : \"[GiD_AccessValue get gendata Write_conditions]\","
    puts $FileVar "                \"MultiFileFlag\"         : \"[GiD_AccessValue get gendata Multi_file_flag]\""
    puts $FileVar "            \},"
    puts $FileVar "            \"output_control_type\"   : \"[GiD_AccessValue get gendata Output_control_type]\","
    puts $FileVar "            \"output_frequency\"      : [GiD_AccessValue get gendata Output_frequency],"
    puts $FileVar "            \"body_output\"           : [GiD_AccessValue get gendata Body_output],"
    puts $FileVar "            \"node_output\"           : [GiD_AccessValue get gendata Node_output],"
    puts $FileVar "            \"skin_output\"           : [GiD_AccessValue get gendata Skin_output],"
    puts $FileVar "            \"plane_output\"          : \[\],"
    # nodal_results
    set PutStrings \[
    if {[GiD_AccessValue get gendata Variables] eq "Conserved"} {
    append PutStrings \"MOMENTUM\" ,
    }
    append PutStrings \"VELOCITY\",\"HEIGHT\",\"FREE_SURFACE_ELEVATION\",\"BATHYMETRY\"
    set PutStrings [string trimright $PutStrings ,]
    append PutStrings \]
    puts $FileVar "            \"nodal_results\"         : $PutStrings,"
    puts $FileVar "            \"gauss_point_results\"   : \[\]"
    puts $FileVar "        \},"
    puts $FileVar "        \"point_data_configuration\"  :  \[\]"
    puts $FileVar "    \},"


    # processes
    puts $FileVar "    \"processes\"    : \{"

    # initial conditions
    set Groups [GiD_Info conditions Initial_water_level groups]
    set NumGroups [llength $Groups]
    set iGroup 0
    puts $FileVar "    \"initial_conditions_process_list\"   : \[\{"
    WriteInitialWaterLevelProcess FileVar iGroup $Groups surfaces $NumGroups

    # boundary conditions
    set Groups [GiD_Info conditions Slip_condition groups]
    set NumGroups [llength $Groups]
    set Groups [GiD_Info conditions Water_height groups]
    incr NumGroups [llength $Groups]
    set Groups [GiD_Info conditions Imposed_flux groups]
    incr NumGroups [llength $Groups]
    set iGroup 0
    puts $FileVar "    \"boundary_conditions_process_list\"  : \[\{"
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
    puts $FileVar "    \"bathymetry_process_list\"     : \[\{"
    WriteBathymetryProcess FileVar iGroup $Groups surfaces $NumGroups

    # end of processes
    puts $FileVar "    \}"


    # Finish ProjectParameters.json file
    puts $FileVar "\}"

    close $FileVar
}

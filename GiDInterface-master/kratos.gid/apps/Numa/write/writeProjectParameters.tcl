### Project Parameters
proc Numa::write::getParametersDict { } {
    
    set projectParametersDict [dict create]
    
    ### Problem data
    ### Create section
    set problemDataDict [dict create]
    set numaTypeofProblem [write::getValue NumaTypeofProblem]
    
    ### Add items to section
    set model_name [file tail [GiD_Info Project ModelName]]
    dict set problemDataDict problem_name $model_name
    dict set problemDataDict model_part_name "MainModelPart"
    set nDim [expr [string range [write::getValue nDim] 0 0] ]
    dict set problemDataDict domain_size $nDim
    dict set problemDataDict parallel_type "OpenMP"
    dict set problemDataDict number_of_threads 1
    dict set problemDataDict start_time [write::getValue NumaTimeParameters StartTime]
    dict set problemDataDict end_time [write::getValue NumaTimeParameters EndTime]
    dict set problemDataDict time_step [write::getValue NumaTimeParameters DeltaTime]
    dict set problemDataDict time_scale [write::getValue NumaTimeParameters TimeScale]
    dict set problemDataDict streamlines_utility [Numa::write::StremalinesUtility]

    ### Add section to document
    dict set projectParametersDict problem_data $problemDataDict

    ### Solver Data   
    set solversettingsDict [dict create]
    
    ### Preguntar el solver haciendo los ifs correspondientes
	if {$numaTypeofProblem eq "Thermo-Mechanical"} {
		dict set solversettingsDict solver_type "dam_thermo_mechanic_solver"
	} else {
		dict set solversettingsDict solver_type "dam_mechanical_solver"
	} 

    set modelDict [dict create]
    dict set modelDict input_type "mdpa"
    dict set modelDict input_filename $model_name
    dict set modelDict input_file_label 0
    dict set solversettingsDict model_import_settings $modelDict
    dict set solversettingsDict echo_level 1
    dict set solversettingsDict buffer_size 2
    
    dict set solversettingsDict processes_sub_model_part_list [write::getSubModelPartNames "NumaNodalConditions" "NumaLoads" "NumaCalibration"]
        
    if {$numaTypeofProblem eq "Thermo-Mechanical" } {
            
        ## Thermal part

        dict set solversettingsDict reference_temperature [write::getValue NumaThermalReferenceTemperature]
        dict set solversettingsDict processes_sub_model_part_list [write::getSubModelPartNames "NumaNodalConditions" "NumaLoads" "NumaCalibration"]
        
        set thermalsettingDict [dict create]
        dict set thermalsettingDict echo_level 1
        dict set thermalsettingDict reform_dofs_at_each_step false
        dict set thermalsettingDict clear_storage false
        dict set thermalsettingDict compute_reactions false
        dict set thermalsettingDict move_mesh_flag false
        dict set thermalsettingDict compute_norm_dx_flag false
        dict set thermalsettingDict theta_scheme [write::getValue NumaThermalScheme]
        dict set thermalsettingDict block_builder false
            
        ## Adding linear solver for thermal part
        set thermal_linear_solver [dict create]
        set thermal_solver_type [write::getValue NumaThermoMechaSolverType]
        if {$thermal_solver_type eq "Direct"} {
            dict set thermal_linear_solver solver_type SuperLUSolver
            dict set thermal_linear_solver scaling false
        } else {
            dict set thermal_linear_solver solver_type AMGCL
            dict set thermal_linear_solver max_iteration 200
            dict set thermal_linear_solver tolerance 1e-7
            dict set thermal_linear_solver provide_coordinates false
            dict set thermal_linear_solver smoother_type ilu0
            dict set thermal_linear_solver krylov_type lgmres
            dict set thermal_linear_solver coarsening_type aggregation
            dict set thermal_linear_solver scaling false
        }

        ## Adding thermal solver settings to solver settings
        dict set thermalsettingDict linear_solver_settings $thermal_linear_solver
        dict set thermalsettingDict problem_domain_sub_model_part_list [Numa::write::getSubModelPartThermalNames]
            
        ## Adding thermal solver settings to solver settings
        dict set solversettingsDict thermal_solver_settings $thermalsettingDict
            

        # Mechanical Part
        set mechanicalSolverSettingsDict [dict create]
        dict set mechanicalSolverSettingsDict solution_type [write::getValue NumaThermoMechaSoluType]
        dict set mechanicalSolverSettingsDict strategy_type Newton-Raphson
        dict set mechanicalSolverSettingsDict scheme_type Newmark
        dict set mechanicalSolverSettingsDict convergence_criterion And_criterion
        dict set mechanicalSolverSettingsDict displacement_relative_tolerance 0.0001
        dict set mechanicalSolverSettingsDict displacement_absolute_tolerance 1e-9
        dict set mechanicalSolverSettingsDict residual_relative_tolerance 0.0001
        dict set mechanicalSolverSettingsDict residual_absolute_tolerance 1e-9
        dict set mechanicalSolverSettingsDict max_iteration 10
        dict set mechanicalSolverSettingsDict echo_level 1
        dict set mechanicalSolverSettingsDict buffer_size 2
        dict set mechanicalSolverSettingsDict compute_reactions false
        dict set mechanicalSolverSettingsDict reform_dofs_at_each_step false
        dict set mechanicalSolverSettingsDict move_mesh_flag false
        dict set mechanicalSolverSettingsDict block_builder false
        dict set mechanicalSolverSettingsDict clear_storage false
        dict set mechanicalSolverSettingsDict rayleigh_m [write::getValue NumaThermoMechaDampMass]
        dict set mechanicalSolverSettingsDict rayleigh_k [write::getValue NumaThermoMechaDampStiff]
        dict set mechanicalSolverSettingsDict nonlocal_damage false

         ## Adding linear solver for mechanical part
        set linear_solver [dict create]
        set mechanical_solver_type [write::getValue NumaThermoMechaSolverType]
        if {$mechanical_solver_type eq "Direct"} {
            dict set linear_solver solver_type SuperLUSolver
            dict set linear_solver scaling false
        } else {
            dict set linear_solver solver_type AMGCL
            dict set linear_solver max_iteration 200
            dict set linear_solver tolerance 1e-7
            dict set linear_solver provide_coordinates false
            dict set linear_solver smoother_type ilu0
            dict set linear_solver krylov_type lgmres
            dict set linear_solver coarsening_type aggregation
            dict set linear_solver scaling false
        }

        ## Adding mechanical solver settings to solver settings
        dict set mechanicalSolverSettingsDict linear_solver_settings $linear_solver

        ### Add section to document
        set mechanicalSolverSettingsDict [dict merge $mechanicalSolverSettingsDict [Numa::write::DefinitionDomains] ]
        ### Add section to document
        dict set solversettingsDict mechanical_solver_settings $mechanicalSolverSettingsDict

    } else {
        
        set mechanicalSolverSettingsDict [dict create]
        dict set mechanicalSolverSettingsDict solution_type [write::getValue NumaMechaSoluType]
        dict set mechanicalSolverSettingsDict strategy_type Newton-Raphson
        dict set mechanicalSolverSettingsDict scheme_type Newmark
        dict set mechanicalSolverSettingsDict convergence_criterion And_criterion
        dict set mechanicalSolverSettingsDict displacement_relative_tolerance 0.0001
        dict set mechanicalSolverSettingsDict displacement_absolute_tolerance 1e-9
        dict set mechanicalSolverSettingsDict residual_relative_tolerance 0.0001
        dict set mechanicalSolverSettingsDict residual_absolute_tolerance 1e-9
        dict set mechanicalSolverSettingsDict max_iteration 10
        dict set mechanicalSolverSettingsDict echo_level 1
        dict set mechanicalSolverSettingsDict buffer_size 2
        dict set mechanicalSolverSettingsDict compute_reactions false
        dict set mechanicalSolverSettingsDict reform_dofs_at_each_step false
        dict set mechanicalSolverSettingsDict move_mesh_flag false
        dict set mechanicalSolverSettingsDict block_builder false
        dict set mechanicalSolverSettingsDict clear_storage false
        dict set mechanicalSolverSettingsDict rayleigh_m [write::getValue NumaMechaDampMass]
        dict set mechanicalSolverSettingsDict rayleigh_k [write::getValue NumaMechaDampStiff]
        dict set mechanicalSolverSettingsDict nonlocal_damage false

         ## Adding linear solver for mechanical part
        set linear_solver [dict create]
        set mechanical_solver_type [write::getValue NumaMechaSolverType]
        if {$mechanical_solver_type eq "Direct"} {
            dict set linear_solver solver_type SuperLUSolver
            dict set linear_solver scaling false
        } else {
            dict set linear_solver solver_type AMGCL
            dict set linear_solver max_iteration 200
            dict set linear_solver tolerance 1e-7
            dict set linear_solver provide_coordinates false
            dict set linear_solver smoother_type ilu0
            dict set linear_solver krylov_type lgmres
            dict set linear_solver coarsening_type aggregation
            dict set linear_solver scaling false
        }

        ## Adding mechanical solver settings to solver settings
        dict set mechanicalSolverSettingsDict linear_solver_settings $linear_solver

        ### Add section to document
        set mechanicalSolverSettingsDict [dict merge $mechanicalSolverSettingsDict [Numa::write::DefinitionDomains] ]
        ### Add section to document
        dict set solversettingsDict mechanical_solver_settings $mechanicalSolverSettingsDict
            
    }
    
    dict set projectParametersDict solver_settings $solversettingsDict
    
    dict set projectParametersDict output_configuration [write::GetDefaultOutputDict]
    
    set nodal_process_list [write::getConditionsParametersDict NumaNodalConditions "Nodal"]
    set load_process_list [write::getConditionsParametersDict NumaLoads ]
    set calibration_process_list [write::getConditionsParametersDict NumaCalibration ]
    
    dict set projectParametersDict constraints_process_list [Numa::write::ChangeFileNameforTableid $nodal_process_list]
    set loads_process_list_prev [Numa::write::ChangeFileNameforTableid $load_process_list]
    dict set projectParametersDict loads_process_list [concat $loads_process_list_prev $calibration_process_list]
    
    dict set projectParametersDict temperature_by_device_list [Numa::write::TemperaturebyDevices]
    dict set projectParametersDict output_device_list [Numa::write::DevicesOutput]   
       
    return $projectParametersDict

}

proc Numa::write::DefinitionDomains { } {
    
 ### Boundary conditions processes
    set domainsDict [dict create]
 
    set body_part_list [list ]
    set joint_part_list [list ]
    set mat_dict [write::getMatDict]
    foreach part_name [dict keys $mat_dict] {
        if {[[Model::getElement [dict get $mat_dict $part_name Element]] getAttribute "ElementType"] eq "Solid"} {
            lappend body_part_list [write::getMeshId Parts $part_name]
            #~ W $body_part_list
        }
    }
    dict set domainsDict problem_domain_sub_model_part_list [write::getSubModelPartNames "NumaParts"]
    dict set domainsDict body_domain_sub_model_part_list $body_part_list

    set loads_sub_model_part_list [list]
    set loads_variable_list [list]
    dict set domainsDict loads_sub_model_part_list $loads_sub_model_part_list
    dict set domainsDict loads_variable_list $loads_variable_list

    return $domainsDict
   
}

proc Numa::write::ChangeFileNameforTableid { processList } {
    set returnList [list ]
    foreach nodalProcess $processList {
        set processName [dict get $nodalProcess process_name]
        set process [::Model::GetProcess $processName]
        set params [$process getInputs]
        foreach {paramName param} $params {
            if {[$param getType] eq "tablefile" && [dict exists $nodalProcess Parameters $paramName] } {
                set filename [dict get $nodalProcess Parameters $paramName]
                set value [Numa::write::GetTableidFromFileid $filename]
                dict set nodalProcess Parameters $paramName $value
            }
            if {[$param getType] eq "vector" && [$param getAttribute vectorType] eq "tablefile" && [dict exists $nodalProcess Parameters $paramName] } {
                for {set i 0} {$i < [llength [dict get $nodalProcess Parameters $paramName]]} {incr i} {
                    set filename [lindex [dict get $nodalProcess Parameters $paramName] $i]
                    set value [Numa::write::GetTableidFromFileid $filename]
                    set values_list [dict get $nodalProcess Parameters $paramName]
                    set values_list [lreplace $values_list $i $i $value]
                    dict set nodalProcess Parameters $paramName $values_list
                }
            }
        }
        lappend returnList $nodalProcess
    }
    return $returnList
}

proc Numa::write::writeParametersEvent { } {
    set projectParametersDict [getParametersDict]
    write::WriteJSON $projectParametersDict
}

proc Numa::write::StremalinesUtility { } {

    set nodalList [write::GetResultsList NodalResults]
    if {[lsearch $nodalList Vi_POSITIVE] >= 0 || [lsearch $nodalList Viii_POSITIVE] >= 0} {
        set streamlines true
    } {
        set streamlines false
    }
    return $streamlines
}

proc Numa::write::DevicesOutput { } {
        
    set output_state [write::getValue NumaOutputState]
    set lista [list ]

    if { $output_state == true } {

        set root [customlib::GetBaseRoot]
        set xp1 "[spdAux::getRoute NumaDevices]/blockdata\[@n='device'\]"
        set nodes [$root selectNodes $xp1]

        foreach node $nodes {

            set deviceDict [dict create]
            dict set deviceDict python_module "point_output_process"
            dict set deviceDict kratos_module "KratosMultiphysics"
            dict set deviceDict help "This process print the selected value according its position"
            dict set deviceDict process_name "PointOutputProcess"

	    	set name [$node @name]
	    	set extension ".txt"
    
            set xp2 "[spdAux::getRoute NumaDevices]/blockdata\[@name='$name'\]/value\[@n='Variable'\]"
	    	set node_xp2 [$root selectNodes $xp2]
	    	set variable [get_domnode_attribute $node_xp2 v]
    
	    	set xp3 "[spdAux::getRoute NumaDevices]/blockdata\[@name='$name'\]/value\[@n='XPosition'\]"
	    	set node_xp3 [$root selectNodes $xp3]
	    	set xposition [write::getValueByNode $node_xp3]
    
	    	set xp4 "[spdAux::getRoute NumaDevices]/blockdata\[@name='$name'\]/value\[@n='YPosition'\]"
	    	set node_xp4 [$root selectNodes $xp4]
	    	set yposition [write::getValueByNode $node_xp4]

	    	set xp5 "[spdAux::getRoute NumaDevices]/blockdata\[@name='$name'\]/value\[@n='ZPosition'\]"
	    	set node_xp5 [$root selectNodes $xp5]
	    	set zposition [write::getValueByNode $node_xp5]

            set parameterDict [dict create]
            set positionList [list ]
            lappend positionList $xposition $yposition $zposition
            dict set parameterDict position $positionList
            dict set parameterDict model_part_name "MainModelPart"
            dict set parameterDict output_file_name $name$extension
            set outputlist [list ]
            lappend outputlist $variable
            dict set parameterDict output_variables $outputlist
            dict set deviceDict Parameters $parameterDict     

            lappend lista $deviceDict
        }
    }

    return $lista

}

proc Numa::write::TemperaturebyDevices { } {

    set device_temp_state [write::getValue NumaTemperatureState]
    set lista [list ]

    if { $device_temp_state == true} {

        set root [customlib::GetBaseRoot]
        set xp1 "[spdAux::getRoute NumaTempDevice]/blockdata\[@n='device'\]"
        set nodes [$root selectNodes $xp1]

        foreach node $nodes {

            set TempdeviceDict [dict create]
            dict set TempdeviceDict  python_module "impose_temperature_by_device_process"
            dict set TempdeviceDict  kratos_module "KratosMultiphysics.DamApplication"
            dict set TempdeviceDict  help "This process assigns a nodal temperature value according the device spatial position"
            dict set TempdeviceDict  process_name "ImposeTemperaturebyDeviceProcess"

	    	set name [$node @name]
    
            set xp2 "[spdAux::getRoute NumaTempDevice]/blockdata\[@name='$name'\]/value\[@n='is_fixed'\]"
	    	set node_xp2 [$root selectNodes $xp2]
	    	set isfixed [get_domnode_attribute $node_xp2 v]

            set xp3 "[spdAux::getRoute NumaTempDevice]/blockdata\[@name='$name'\]/value\[@n='value'\]"
	    	set node_xp3 [$root selectNodes $xp3]
	    	set value [write::getValueByNode $node_xp3]
    
	    	set xp4 "[spdAux::getRoute NumaTempDevice]/blockdata\[@name='$name'\]/value\[@n='table'\]"
	    	set node_xp4 [$root selectNodes $xp4]
	    	set table [write::getValueByNode $node_xp4]

	    	set xp5 "[spdAux::getRoute NumaTempDevice]/blockdata\[@name='$name'\]/value\[@n='XPosition'\]"
	    	set node_xp5 [$root selectNodes $xp5]
	    	set xposition [write::getValueByNode $node_xp5]
    
	    	set xp6 "[spdAux::getRoute NumaTempDevice]/blockdata\[@name='$name'\]/value\[@n='YPosition'\]"
	    	set node_xp6 [$root selectNodes $xp6]
	    	set yposition [write::getValueByNode $node_xp6]

	    	set xp7 "[spdAux::getRoute NumaTempDevice]/blockdata\[@name='$name'\]/value\[@n='ZPosition'\]"
	    	set node_xp7 [$root selectNodes $xp7]
	    	set zposition [write::getValueByNode $node_xp7]

            set parameterDict [dict create]
            set positionList [list ]
            dict set parameterDict mesh_id 0
            dict set parameterDict model_part_name "MainModelPart"
            dict set parameterDict variable_name "TEMPERATURE"
            dict set parameterDict is_fixed $isfixed
            dict set parameterDict value $value
            dict set parameterDict table 0
            lappend positionList $xposition $yposition $zposition
            dict set parameterDict position $positionList

            dict set TempdeviceDict  Parameters $parameterDict     

            lappend lista $TempdeviceDict 
        }
    }

    return $lista
}

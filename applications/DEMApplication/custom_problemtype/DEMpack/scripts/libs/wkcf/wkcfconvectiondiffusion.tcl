###############################################################################
#
#    NAME: wkcfconvectiondiffusion.tcl
#
#    PURPOSE: Useful procedures to work with the convection-diffusion application
#
#    QUANTECH ATZ-DEVELOPMENT DEPARTMENT
#
#    AUTHOR : G. Socorro
#
#    CREATED AT: 12/09/13
#
#    HISTORY:
#
#     0.2- 27/09/13-G. Socorro, create some new proc GetConvectionDiffusionMaterialProperties,WriteConvectionDiffusionPrescribedHeatFluxBC
#                               and WriteConvectionDiffusionProjectParameters
#     0.1- 12/09/13-G. Socorro, create a base source code from wkcf.tcl
#
###############################################################################


proc ::wkcf::GetPropertyIdForFaceHeatFluxBC {AppId GroupId} {
    # ABSTRACT: Get the property identifier from the group identifier
    variable dprops

    set RPropertyId 1
    foreach PropertyId $dprops($AppId,GBCKProps,AllPropertyId) {
	# Get the propertyid-group identifier
	set PropertyGroupId $dprops($AppId,GBCKProps,$PropertyId,PropertyGroupId)
	if {$PropertyGroupId eq $GroupId} {
	    set RPropertyId $PropertyId
	}
    }
    return $RPropertyId
}

proc ::wkcf::GetPropertyDataFromFaceHeatFluxBC {AppId} {
    # ABSTRACT: Get the mapping between property and face heat flux boundary condition
    variable dprops
   
    # For debug
    if {!$::wkcf::pflag} {
        set inittime [clock seconds]
    }

    # Init the property to BC counter
    set bcgroupid -1

    # Solver type
    set cproperty "dv"
    set rootid $AppId
    set cxpath "$rootid//c.AnalysisData//i.AnalysisType"
    set AnalysisType [::xmlutils::setXml $cxpath $cproperty]
    if {$AnalysisType eq "Linear"} {
	set bcgroupid 0
    } 

    # Init the global boundary condition identifier list
    set dprops($AppId,GBCKProps,AllPropertyId) [list]
   
    # Set the rootid
    set rootid "$AppId"
    set cproperty "dv"
    
    # Get the values
    set basexpath "$rootid//c.Conditions//c.FaceHeatFlux"
    set gproplist [::xmlutils::setXmlContainerIds $basexpath]
    # wa "gproplist:$gproplist"
    foreach cgroupid $gproplist {
	# Get the group properties
	set cxpath "${basexpath}//c.[list ${cgroupid}]//c.MainProperties"
	set allgprop [::xmlutils::setXmlContainerPairs $cxpath "" "dv"]
	# wa "allgprop:$allgprop"
	if {[llength $allgprop]} {
	    if {[GiD_EntitiesGroups get $cgroupid nodes -count]} {
		incr bcgroupid 1
		# Create the propertyid-group identifier mapping
		lappend dprops($AppId,GBCKProps,AllPropertyId) $bcgroupid
	
		# Set the properties
		# Modify some properties
		set newprops [list]
	       	foreach propid $allgprop {
		    lassign $propid key value
		    if {($key eq "VFaceHeatFlux")} {
			continue
		    } else {
			lappend newprops [list $key $value]	
		    }
		}
		# wa "newprops:$newprops"
		if {[llength $newprops]} {
		    set dprops($AppId,GBCKProps,$bcgroupid,PropertyList) $newprops
		    set dprops($AppId,GBCKProps,$bcgroupid,PropertyGroupId) $cgroupid
		}
	    }
	}
    }
    
    # For debug
    if {!$::wkcf::pflag} {
	set endtime [clock seconds]
	set ttime [expr $endtime-$inittime]
	# WarnWinText "endtime:$endtime ttime:$ttime"
	WarnWinText "Get the mapping between property and face heat flux boundary condition: [::KUtils::Duration $ttime]"
    }
}

proc ::wkcf::GetConvectionDiffusionMaterialProperties {AppId} {
    # ABSTRACT : Return the material properties for the convection diffusion application
    variable dprops
   
    set Density 0.0; set Conductivity 0.0; set ExpansionCoeff 0.0; set SpecificHeat 0.0
    set flag [expr {([info exists dprops($AppId,AllKElemId)]) && ([llength $dprops($AppId,AllKElemId)])}]
    # wa "flag:$flag"
    
    if {$flag} {
	set cproplist [list "Density" "ExpansionCoeff" "Conductivity" "SpecificHeat"]
	foreach PropertyId $dprops($AppId,GKProps,AllPropertyId) {
	    # Get the material identifier for this property 
	    set MatId $dprops($AppId,Property,$PropertyId,MatId) 
	    # Get the group identifier
	    set GroupId $dprops($AppId,Property,$PropertyId,GroupId)
	    # Get all material properties
	    set mpxpath "[::KMat::findMaterialParent $MatId]//m.[list ${MatId}]"
	    # WarnWinText "mpxpath:$mpxpath"
	    # Get the material properties
	    foreach pid $cproplist {
		if {$pid =="Density"} {
		    set xpath "c.General"
		    # Get the current value for this properties
		    set cvalue [lindex [::KMat::getMaterialProperties "p" "$mpxpath//$xpath//p.[list $pid]"] 0 1]
		    set Density [GiD_FormatReal "%10.5e" $cvalue]
		} elseif {$pid =="ExpansionCoeff"} {
		    set xpath "c.Thermal//c.Isotropic"
		    # Get the current value for this properties
		    set cvalue [lindex [::KMat::getMaterialProperties "p" "$mpxpath//$xpath//p.[list $pid]"] 0 1]
		    set ExpansionCoeff [GiD_FormatReal "%10.5e" $cvalue]
		} elseif {$pid =="Conductivity"} {
		    set xpath "c.Thermal//c.Isotropic"
		    # Get the current value for this properties
		    set cvalue [lindex [::KMat::getMaterialProperties "p" "$mpxpath//$xpath//p.[list $pid]"] 0 1]
		    set Conductivity [GiD_FormatReal "%10.5e" $cvalue]
		    # wa "Conductivity:$Conductivity"
		} elseif {$pid =="SpecificHeat"} {
		    set xpath "c.Thermal//c.Isotropic"
		    # Get the current value for this properties
		    set cvalue [lindex [::KMat::getMaterialProperties "p" "$mpxpath//$xpath//p.[list $pid]"] 0 1]
		    set SpecificHeat [GiD_FormatReal "%10.5e" $cvalue]
		}
	    }
	    # Only the first property
	    break 
	}
    }
  
    return [list $Density $Conductivity $ExpansionCoeff $SpecificHeat]
}

proc ::wkcf::WriteConvectionDiffusionPropertyAtNodes {AppId} {
    # ABSTRACT: Write some properties at the nodal level for Convecion-Diffusion application
    variable dprops
    variable filechannel

    # For debug
    if {!$::wkcf::pflag} {
	set inittime [clock seconds]
    }

    set cproperty "dv"
  
    set flag [expr {([info exists dprops($AppId,AllKElemId)]) && ([llength $dprops($AppId,AllKElemId)])}]
    # wa "flag:$flag"
    # Check for all defined kratos elements
    if {$flag} {
	
	set matkwordlist [list Density Conductivity SpecificHeat]
	set matvarlist [list Density Conductivity SpecificHeat]

	# Write mesh velocity and velocity for each node identifier
	set MeshVelocity 0.0; set Velocity 0.0

	# Get the convection diffusion material properties
	lassign [::wkcf::GetConvectionDiffusionMaterialProperties $AppId] Density Conductivity ExpansionCoeff SpecificHeat

	# wa "Density:$Density Conductivity:$Conductivity ExpansionCoeff:$ExpansionCoeff SpecificHeat:$SpecificHeat"

	set kxpath "Applications/$AppId"
	set cpropid "0"        

	# Write the group nodal properties
	foreach celemid $dprops($AppId,AllKElemId) {
	    # Check for all defined group identifier for this element
	    if {([info exists dprops($AppId,KElem,$celemid,AllGroupId)]) && ([llength $dprops($AppId,KElem,$celemid,AllGroupId)])} {
		# For all defined group identifier for this element
		foreach cgroupid $dprops($AppId,KElem,$celemid,AllGroupId) {                    
		    if {[GiD_EntitiesGroups get $cgroupid nodes -count]} {
			set strnodes ""
			foreach node_id [GiD_EntitiesGroups get $cgroupid nodes] {
			    append  strnodes "[format "%10i%4i" $node_id $cpropid]  [GiD_FormatReal "%10.5e" $MeshVelocity]\n"
			}

			# Write mesh velocity and velocity value for this group
			foreach item [list "MeshVelocity" "Velocity"] {
			    foreach comp [list "X" "Y" "Z"] {
				set vkword [::xmlutils::getKKWord $kxpath "${item}_${comp}" "kkword"]
				GiD_File fprintf $filechannel "%s" "Begin NodalData $vkword \/\/ GUI group identifier: $cgroupid"
				GiD_File fprintf $filechannel "%s"  "$strnodes"
				GiD_File fprintf $filechannel "%s " "End NodalData"
				GiD_File fprintf $filechannel ""
			    }
			}
			
		
			# Write others material properties
			foreach kwordid $matkwordlist varid $matvarlist {
			    set vkword [::xmlutils::getKKWord $kxpath "$kwordid" "kkword"]
			    GiD_File fprintf $filechannel "%s" "Begin NodalData $vkword \/\/ GUI group identifier: $cgroupid"
			    foreach node_id [GiD_EntitiesGroups get $cgroupid nodes] {
				GiD_File fprintf $filechannel "%10i %4i %10.5e" $node_id $cpropid [set $varid]
			    }
			    GiD_File fprintf $filechannel "%s" "End NodalData"
			    GiD_File fprintf $filechannel "%s" ""
			}
		    }
		}
	    }
	}
    }

    # For debug
    if {!$::wkcf::pflag} {
	set endtime [clock seconds]
	set ttime [expr $endtime-$inittime]
	# WarnWinText "endtime:$endtime ttime:$ttime"
	WarnWinText "Write property at nodes: [::KUtils::Duration $ttime]"
    }
}


proc ::wkcf::WriteConvectionDiffusionPrescribedTemperatureBC {AppId ccondid kwordlist} {
    variable ndime; variable dprops; variable filechannel

    set cpropid "1"
    
    # Convection-diffusion prescribed temperature
    set xitem [lindex $kwordlist 0]
  
    # For each group in the this condition
    foreach cgroupid $dprops($AppId,BC,$ccondid,AllGroupId) {
	# wa "cgroupid:$cgroupid"
	# Get the condition properties
	lassign $dprops($AppId,BC,$ccondid,$cgroupid,GProps) ixvar ixval
	# wa "x: $ixvar $ixval"
	# Fixed temperature
	if {$ixvar} {
	    if { [GiD_EntitiesGroups get $cgroupid nodes -count] } {
		GiD_File fprintf $filechannel "%s" "Begin NodalData $xitem \/\/ Convection-Diffusion temperature condition GUI group identifier: $cgroupid"
		foreach node_id [GiD_EntitiesGroups get $cgroupid nodes] {
		    GiD_File fprintf $filechannel "%10i %8i %10.5e" $node_id $cpropid $ixval
		}
		GiD_File fprintf $filechannel "End NodalData"
		GiD_File fprintf $filechannel ""
	    }
	}
    }
}

proc ::wkcf::WriteConvectionDiffusionPrescribedHeatFluxBC {AppId ccondid kwordlist} {
    variable ndime; variable dprops; variable filechannel

    set cpropid "0"
 
    # Convection-diffusion heat flux boundary condition
    set xitem [lindex $kwordlist 0]
  
    # For each group in the this condition
    foreach cgroupid $dprops($AppId,BC,$ccondid,AllGroupId) {
	# wa "cgroupid:$cgroupid"
	# Get the condition properties
	lassign $dprops($AppId,BC,$ccondid,$cgroupid,GProps) heatfluxvalue
	# wa "heatfluxvalue:$heatfluxvalue"
	# Heat flux value
	if { [GiD_EntitiesGroups get $cgroupid nodes -count] } {
	    GiD_File fprintf $filechannel "%s" "Begin NodalData $xitem \/\/ Convection-Diffusion heat flux condition GUI group identifier: $cgroupid"
	    foreach node_id [GiD_EntitiesGroups get $cgroupid nodes] {
		GiD_File fprintf $filechannel "%10i %8i   %10.5e" $node_id $cpropid $heatfluxvalue
	    }
	    GiD_File fprintf $filechannel "End NodalData"
	    GiD_File fprintf $filechannel ""
	}
    }
}

proc ::wkcf::WriteConvectionDiffusionPrescribedFaceHeatFluxBC {AppId ccondid kwordlist} {
    variable ndime; variable dprops; variable filechannel

    set cpropid "1"
 
    # Convection-diffusion face heat flux boundary condition
    set xitem [lindex $kwordlist 0]
  
    # For each group in the this condition
    foreach cgroupid $dprops($AppId,BC,$ccondid,AllGroupId) {
	# wa "cgroupid:$cgroupid"
	# Get the condition properties
	lassign $dprops($AppId,BC,$ccondid,$cgroupid,GProps) faceheatfluxvalue
	# wa "heatfluxvalue:$heatfluxvalue"
	# Face heat flux value
	if { [GiD_EntitiesGroups get $cgroupid nodes -count] } {
	    # Get the property identifier from the group property mapping
	    set cpropid [::wkcf::GetPropertyIdForFaceHeatFluxBC $AppId $cgroupid]
	    GiD_File fprintf $filechannel "%s" "Begin NodalData $xitem \/\/ Convection-Diffusion face heat flux condition GUI group identifier: $cgroupid"
	    foreach node_id [GiD_EntitiesGroups get $cgroupid nodes] {
		GiD_File fprintf $filechannel "%10i %8i   %10.5e" $node_id $cpropid $faceheatfluxvalue
	    }
	    GiD_File fprintf $filechannel "End NodalData"
	    GiD_File fprintf $filechannel ""
	}
    }
}

proc ::wkcf::WriteInitialTemperatureIC {AppId ccondid kwordlist} {
    variable ndime; variable dprops; variable filechannel

    set cpropid "0"
 
    # Convection-diffusion initial temperature initial condition
    set xitem [lindex $kwordlist 0]
  
    # For each group in the this condition
    foreach cgroupid $dprops($AppId,IC,$ccondid,AllGroupId) {
	# wa "cgroupid:$cgroupid"
	# Get the condition properties
	lassign $dprops($AppId,IC,$ccondid,$cgroupid,GProps) initialtempvalue
	# wa "initialtempvalue:$initialtempvalue"
	# Initial temperature value
	if { [GiD_EntitiesGroups get $cgroupid nodes -count] } {
	    GiD_File fprintf $filechannel "%s" "Begin NodalData $xitem \/\/ Convection-Diffusion initial temperature condition GUI group identifier: $cgroupid"
	    foreach node_id [GiD_EntitiesGroups get $cgroupid nodes] {
		GiD_File fprintf $filechannel "%10i %8i   %10.5e" $node_id $cpropid $initialtempvalue
	    }
	    GiD_File fprintf $filechannel "End NodalData"
	    GiD_File fprintf $filechannel ""
	}
    }
}

proc ::wkcf::WriteConvectionDiffusionProjectParameters {AppId fileid PDir} {
    variable ndime; variable dprops
   
    set rootid "$AppId"
    
    set trailing_spaces  "    "

    # Kratos key word xpath
    set kxpath "Applications/$rootid"
    
    # Write linear solver method
    # 0 => Old format
    # 1 => New format
    set wfsmethod 1

    # Get others properties
    set cproperty "dv"

    set domain_size 2
    if {$ndime =="2D"} {
	set domain_size 2
    } elseif {$ndime =="3D"} {
	set domain_size 3
    }
    
    # Domain size
    puts $fileid "domain_size = $domain_size"
    puts $fileid ""
    
    # Solver type
    set cxpath "$rootid//c.AnalysisData//i.AnalysisType"
    set AnalysisType [::xmlutils::setXml $cxpath $cproperty]
    if {$AnalysisType eq "Linear"} {
	puts $fileid "Linear = True"
    } else {
	puts $fileid "Linear = False"
    }
    puts $fileid ""

    puts $fileid "class SolverSettings1:"
    puts $fileid "    solver_type = \"convection_diffusion_solver\""
    puts $fileid "    domain_size= $domain_size"
    # Time order
    set cxpath "$rootid//c.SolutionStrategy//c.Advanced//i.TimeOrder"
    set TimeOrder [::xmlutils::setXml $cxpath $cproperty]
    puts $fileid "${trailing_spaces}time_order = $TimeOrder"
    
    # Predictor corrector
    set cxpath "$rootid//c.SolutionStrategy//c.Advanced//i.PredictorCorrector"
    set PredictorCorrector [::xmlutils::setXml $cxpath $cproperty]
    puts $fileid "${trailing_spaces}predictor_corrector = $PredictorCorrector"	
    puts $fileid "    ReformDofAtEachIteration = False"
    puts $fileid "    echo_level=0"
    
    puts $fileid ""
    puts $fileid ""

    puts $fileid "    ###set variables to be used"
    puts $fileid "    unknown_variable = \"TEMPERATURE\""
    puts $fileid "    density_variable= \"DENSITY\""
    puts $fileid "    projection_variable= \"TEMP_CONV_PROJ\""
    puts $fileid "    volume_source_variable= \"HEAT_FLUX\""
    puts $fileid "    diffusion_variable= \"CONDUCTIVITY\""
    puts $fileid "    surface_source_variable= \"FACE_HEAT_FLUX\""
    puts $fileid "    mesh_velocity_variable= \"MESH_VELOCITY\""
    puts $fileid "    velocity_variable= \"VELOCITY\""
    puts $fileid "    specific_heat_variable= \"SPECIFIC_HEAT\""

    puts $fileid ""
    puts $fileid ""

    # Linear solvers
    # Temperature
    ::wkcf::WriteLinearSolvers $rootid $fileid "Temperature" $wfsmethod $trailing_spaces "linear_solver_config"
    
    puts $fileid ""
    puts $fileid ""
    
    puts $fileid "class SolverSettings2:"
    puts $fileid "    solver_type = \"nonlinear_convection_diffusion_solver\""
    puts $fileid "    domain_size= $domain_size"
    
    # Time order
    set cxpath "$rootid//c.SolutionStrategy//c.Advanced//i.TimeOrder"
    set TimeOrder [::xmlutils::setXml $cxpath $cproperty]
    puts $fileid "${trailing_spaces}time_order = $TimeOrder"
    
    # Predictor corrector
    set cxpath "$rootid//c.SolutionStrategy//c.Advanced//i.PredictorCorrector"
    set PredictorCorrector [::xmlutils::setXml $cxpath $cproperty]
    puts $fileid "${trailing_spaces}predictor_corrector = $PredictorCorrector"	
    puts $fileid "    ReformDofAtEachIteration = False"
    puts $fileid "    echo_level=0"
    
    # Maximum temperature iterations
    puts $fileid "$trailing_spaces"
    puts $fileid "${trailing_spaces}#convergence criteria settings"
    set cxpath "$rootid//c.SolutionStrategy//c.Advanced//i.MaximumTemperatureIterations"
    set cvalue [::xmlutils::setXml $cxpath $cproperty]
    puts $fileid "${trailing_spaces}max_iter = $cvalue"
    
    # Relative temperature tolerance
    set cxpath "$rootid//c.SolutionStrategy//c.Advanced//i.RelativeTemperatureTolerance"
    set cvalue [::xmlutils::setXml $cxpath $cproperty]
    puts $fileid "${trailing_spaces}toll = $cvalue"
    
    puts $fileid ""
    puts $fileid ""

    puts $fileid "    ###set variables to be used"
    puts $fileid "    unknown_variable = \"TEMPERATURE\""
    puts $fileid "    density_variable= \"DENSITY\""
    puts $fileid "    projection_variable= \"TEMP_CONV_PROJ\""
    puts $fileid "    volume_source_variable= \"HEAT_FLUX\""
    puts $fileid "    diffusion_variable= \"CONDUCTIVITY\""
    puts $fileid "    surface_source_variable= \"FACE_HEAT_FLUX\""
    puts $fileid "    mesh_velocity_variable= \"MESH_VELOCITY\""
    puts $fileid "    velocity_variable= \"VELOCITY\""
    puts $fileid "    specific_heat_variable= \"SPECIFIC_HEAT\""

    puts $fileid ""
    puts $fileid ""

    # Linear solvers
    # Temperature
    ::wkcf::WriteLinearSolvers $rootid $fileid "Temperature" $wfsmethod $trailing_spaces "linear_solver_config"

    puts $fileid ""

    puts $fileid "nodal_results=\[\"TEMPERATURE\"\]"
    puts $fileid "gauss_points_results=\[\]"

    # GiD post mode variables
    ::wkcf::WriteGiDPostMode $AppId $fileid 

    puts $fileid ""

    puts $fileid "AutomaticDeltaTime = \"Fixed\""
    # Start time
    set cxpath "$rootid//c.SolutionStrategy//c.ProblemParameters//i.StartTime"
    set StartTime [::xmlutils::setXml $cxpath $cproperty]
    # End time
    set cxpath "$rootid//c.SolutionStrategy//c.ProblemParameters//i.EndTime"
    set EndTime [::xmlutils::setXml $cxpath $cproperty]
    # Delta time
    set cxpath "$rootid//c.SolutionStrategy//c.ProblemParameters//i.DeltaTime"
    set DeltaTime [::xmlutils::setXml $cxpath $cproperty]
    
     # WarnWinText "StartTime:$StartTime EndTime:$EndTime DeltaTime:$DeltaTime"
    puts $fileid "Dt = $DeltaTime"
    puts $fileid "Start_time = $StartTime"
    puts $fileid "max_time = $EndTime"
    set nsteps [expr int(double($EndTime-$StartTime)/double($DeltaTime))]
    puts $fileid "nsteps = $nsteps"
    
    puts $fileid ""

    puts $fileid "#output settings"
    # Output step 
    set cxpath "$rootid//c.Results//i.OutputDeltaTime"
    set OutputDeltaTime [::xmlutils::setXml $cxpath $cproperty]
    puts $fileid "output_time = $OutputDeltaTime"
    # WarnWinText "OutputDeltaTime:$OutputDeltaTime"
    set output_step [expr int($OutputDeltaTime/double($DeltaTime))]
    # WarnWinText "output_step:$output_step"
    puts $fileid "output_step = $output_step"
    # Volumen Output
    set cxpath "$rootid//c.Results//i.VolumeOutput"
    set VolumeOutput [::xmlutils::setXml $cxpath $cproperty]
    if {$VolumeOutput eq "Yes"} {
	puts $fileid "VolumeOutput = True"
    } else {
	puts $fileid "VolumeOutput = False"
    }

    puts $fileid ""
    set PName [::KUtils::GetPaths "PName"]
    puts $fileid "problem_name=\"${PName}${AppId}\"" 
    puts $fileid "problem_path=\"$PDir\"" 
    
    # Get the kratos path 
    set cxpath "GeneralApplicationData//c.ProjectConfiguration//i.KratosPath"
    set cproperty "dv"
    set KratosPath [::xmlutils::setXml $cxpath $cproperty]
    set KratosPath [file native $KratosPath]
    
    # Write the kratos path
    puts $fileid "kratos_path=\"${KratosPath}\"" 
}

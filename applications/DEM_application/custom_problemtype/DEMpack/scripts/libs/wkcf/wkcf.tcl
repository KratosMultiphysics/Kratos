package require tdom
package require xpathq

namespace eval ::wkcf:: {
    # Spatial dimension
    variable ndime
    # Structural analysis
    variable StructuralAnalysis
    # Fluid application
    variable FluidApplication
    # FSI application
    variable FSIApplication
    # ConvectionDiffusion application
    variable ConvectionDiffusionApplication
    # DEM application
    variable DEMApplication
    # Active application list
    variable ActiveAppList

    variable gidetype [list Linear Triangle Quadrilateral Tetrahedra Hexahedra]
    # Check for use quadratic elements
    variable useqelem
    # GiD entity list
    variable gidentitylist [list "point" "line" "surface" "volume"]

    # Properties array
    variable dprops

    # Conditions to BC link
    variable ctbclink

    # Debug/Release variable
    variable pflag

    # To write the bat file
    variable wbatfile

    # filechannel -> output file .mdpa
    variable filechannel
    variable deminletchannel
    variable demfemchannel
    variable clusterschannel

    # Structural analysis condition counter
    variable sa_icondid
}

proc ::wkcf::GetAnalysisDataContainer {} {
    set cxpath "GeneralApplicationData//c.CouplingParameters//i.CouplingLevel"
    set coupling_level_type [::xmlutils::setXml $cxpath "dv"]
    if {$coupling_level_type == "0"} {
	return "AnalysisDataCouplingL0"
    } elseif {$coupling_level_type == "1"} {
	return "AnalysisDataCouplingL1"
    } elseif {$coupling_level_type == "2"} {
	return "AnalysisDataCouplingL2"
    } else {
	return "AnalysisData"
    }
}

proc ::wkcf::GetOptionsContainer {} {
    set analysis_data_path [GetAnalysisDataContainer]
    set cxpath "Fluid//c.$analysis_data_path//i.SolverType"
    set solver_type [::xmlutils::setXml $cxpath "dv"]
    if {$solver_type == "ElementBased"} {
	return "OptionsFractionalStep"
    } else {
	return "Options"
    }
}

proc ::wkcf::WriteCalculationFiles {filename} {
    variable FluidApplication
    variable StructuralAnalysis
    variable ActiveAppList
    variable DSOLID
    variable filechannel
    variable ConvectionDiffusionApplication
    variable DEMApplication
    variable deminletchannel
    variable demfemchannel
    variable property_number
    set property_number 1
    variable dem_props_number
    set dem_props_number  1
    variable dem_ref_to_props_number
    set dem_ref_to_props_number  1
    variable dem_group_mesh_property_number
    set dem_group_mesh_property_number 0
    global KPriv

    set KPriv(list_of_cluster_files) [list]

    set fail 0
    set msg ""

    # Unset some local variables
    ::wkcf::UnsetLocalVariables

    # Init some namespace global variables
    ::wkcf::Preprocess

    # Write each block of the files *.mdpa
    # Rename the file name => Change .dat by .mdpa
    set basefilename [file rootname $filename]

    # Write a mdpa file for each application
    foreach AppId $ActiveAppList {
	  if {$AppId eq "DEM"} {
	    # Inlet Channel
	    set filenameInlet ${basefilename}${AppId}_Inlet.mdpa
	    set deminletchannel [GiD_File fopen $filenameInlet]

	    # MDPA standard channel
	    set filename ${basefilename}${AppId}.mdpa
	    set filechannel [GiD_File fopen $filename]

	    # FEM Boundary channel
	    set demfilenameBoun ${basefilename}${AppId}_FEM_boundary.mdpa
	    set demfemchannel [GiD_File fopen $demfilenameBoun]

	    # Clusters Channel
	    set filenameClusters ${basefilename}${AppId}_Clusters.mdpa
	    set clusterschannel [GiD_File fopen $filenameClusters]

	  } else {
	    set filename ${basefilename}${AppId}.mdpa
	    # Open the file
	    set filechannel [GiD_File fopen $filename]
	}
	# wa "filename:$filename"

	# Write model part data
	if {$AppId ne "DSOLID"} {
	    ::wkcf::WriteModelPartData $AppId
	}

	if {$AppId eq "DEM"} {
	    # Write properties block
	    lassign [::wkcf::WriteDSOLIDProperties $AppId] fail msg
	}

	if {$AppId eq "DEM"} {
	    # Write properties block
	    lassign [::wkcf::WriteDSOLIDNodalCoordinates $AppId] fail msg
	    lassign [::wkcf::WriteDSOLIDVolumeAccelerationOnNodes $AppId] fail msg
	    lassign [::wkcf::WriteDSOLIDElementConnectivities $AppId] fail msg
	    lassign [::wkcf::WriteDSOLIDContactProperties $AppId] fail msg
	    lassign [::wkcf::WriteDSOLIDContactElements $AppId] fail msg
	    lassign [::wkcf::WriteDSOLIDContactKinematics $AppId] fail msg
	}

	# Write properties block
	if {$AppId ne "DSOLID"} {
	    lassign [::wkcf::WriteProperties $AppId] fail msg

	    # Write nodes block
	    if { !$fail } {
		lassign [::wkcf::WriteNodalCoordinates $AppId] fail msg
	    }
	    # Write elements block
	    if { !$fail } {
		::wkcf::WriteElementConnectivities $AppId

		# For fluid or convection-diffusion applications
		if {($AppId == "Fluid")} {
		    ::wkcf::WriteConditions $AppId
		} elseif {($AppId == "ConvectionDiffusion")} {
		    set cproperty "dv"
		    set cxpath "$AppId//c.AnalysisData//i.AnalysisType"
		    set AnalysisType [::xmlutils::setXml $cxpath $cproperty]
		    if {$AnalysisType eq "Non-Linear"} {
		        ::wkcf::WriteConditions $AppId
		    }
		    ::wkcf::WriteInitialConditions $AppId
		}
		if {$KPriv(what_dempack_package) ne "C-DEMPack"} {
		    ::wkcf::WriteBoundaryConditions $AppId
		}

		if {$AppId == "Fluid"} {
		    ::wkcf::WritePropertyAtNodes $AppId
		    ::wkcf::WriteGroupMeshProperties $AppId
		    ::wkcf::WriteCutAndGraph $AppId
		} elseif {$AppId == "DEM"} {
		    ::wkcf::WriteDemNodalVariables2 $AppId $filechannel
		    ::wkcf::WriteGroupMeshProperties $AppId
		    ::wkcf::WriteInletGroupMeshProperties $AppId
			::wkcf::WriteCustomSubModelParts $AppId
		}
	    }
	    ::wkcf::CloseChannels $AppId
	}

    }
    ::wkcf::WriteProjectParameters
    if {$KPriv(what_dempack_package) eq "S-DEMPack"} {
    ::wkcf::WriteDSOLIDMaterialsFile
    }
    ::wkcf::SelectPythonScript
    ::wkcf::CopyClusterDefinitionFiles

    if {$AppId == "DEM"} {
	::wkcf::WriteExplicitSolverVariables
	::wkcf::WriteExplicitSolverVariablesInJsonFile
    }

    variable wbatfile
    if {$wbatfile} {
	# Write bat file only for the Linux OS
	if {($::tcl_platform(os) eq "Linux")} {
	    ::wkcf::WriteBatFile $AppId
	}
    }

    # Unset some local variables
    ::wkcf::UnsetLocalVariables

    return [list $fail $msg]
}

proc ::wkcf::SelectPythonScript {} {
    # Select the correct Python script
    variable FluidApplication
    variable StructuralAnalysis
    variable DEMApplication
    variable DSOLID
    set endfilename "KratosOpenMP.py"
    set mpiendfilename "KratosMPI.py"
    # Get the application root identifier
    set rootdataid [::wkcf::GetApplicationRootId]
    set cproperty "dv"
    set PDir [file native [::KUtils::GetPaths "PDir"]]
    set PTDir [::KUtils::GetPaths "PTDir"]

    # For DEM application
    if {$DEMApplication =="Yes" && $FluidApplication =="No" && $DSOLID =="No"} {
	set endfilename "KratosDEM.py"
	set ppfilename "KratosDEM.py"
	set tofname [file native [file join $PDir $endfilename]]
	set fromfname [file native [file join "$PTDir/python" $ppfilename]]
	if {[catch {file copy -force "$fromfname" "$tofname"} error]} {
	    WarnWin [= "Could not copy the Kratos Python script (%s) to (%s): Error (%s)" $fromfname $tofname $error ]
	    return ""
	}
    }
    if {$DEMApplication =="Yes" && $DSOLID =="Yes"} {
	set ppfilename "KratosDEMFEM.py"
	set endfilename "KratosDEMFEM.py"
	set tofname [file native [file join $PDir $endfilename]]
	set fromfname [file native [file join "$PTDir/python" $ppfilename]]
	if {[catch {file copy -force "$fromfname" "$tofname"} error]} {
	    WarnWin [= "Could not copy the Kratos Python script (%s) to (%s): Error (%s)" $fromfname $tofname $error ]
	    return ""
	}
    }

    # For Swimming DEM application
    if {$DEMApplication =="Yes" && $FluidApplication =="Yes"} {

	set ppfilename "KratosSwimmingDEM.py"
	set cxpath "DEM//c.DEM-Cohesivegroup"
	set glist [::xmlutils::setXmlContainerIds $cxpath]

	if { [llength $glist] } {
	    WarnWin [= "DEM continuum simulations are not available when interacting with fluid."]
	}
	set endfilename "KratosSwimmingDEM.py"
	set tofname [file native [file join $PDir $endfilename]]
	set fromfname [file native [file join "$PTDir/python" $ppfilename]]

	if {[catch {file copy -force "$fromfname" "$tofname"} error]} {
	    WarnWin [= "Could not copy the Kratos Python script (%s) to (%s): Error (%s)" $fromfname $tofname $error ]
	    return ""
	}
    }
}

proc ::wkcf::CopyClusterDefinitionFiles {} {
    global KPriv

    set PDir [file native [::KUtils::GetPaths "PDir"]]
    set PTDir [::KUtils::GetPaths "PTDir"]

    foreach f $KPriv(list_of_cluster_files) {
	set fromfname [file native [file join "$PTDir/clusters" $f]]
	set tofname [file native [file join $PDir $f]]

	if {[catch {file copy -force "$fromfname" "$tofname"} error]} {
	    WarnWin [= "Could not copy the cluster information file (%s) to (%s): Error (%s)" $fromfname $tofname $error ]
	    return ""
	}
    }
}

proc ::wkcf::WriteModelPartData {AppId} {
    # Write the model part data
    # Arguments
    # AppId => Application identifier
    variable filechannel
    variable deminletchannel
    variable demfemchannel
    variable clusterschannel
    global KPriv

    if {$AppId eq "DEM"} {
	foreach channel [list $filechannel $deminletchannel $demfemchannel] {
	    GiD_File fprintf $channel "%s" "Begin ModelPartData"
	    GiD_File fprintf $channel "%s" "// VARIABLE_NAME value"
	    GiD_File fprintf $channel "%s" "End ModelPartData"
	    GiD_File fprintf $channel "%s" ""
	}
    } else {
	GiD_File fprintf $filechannel "%s" "Begin ModelPartData"
	GiD_File fprintf $filechannel "%s" "// VARIABLE_NAME value"
	GiD_File fprintf $filechannel "%s" "End ModelPartData"
	GiD_File fprintf $filechannel "%s" ""
    }
}

proc ::wkcf::WriteProperties {AppId} {
    # Write the properties block
    # Arguments
    # AppId => Application identifier
    variable dprops
    variable ndime
    variable dem_props_number
    global ref_props_counter
    variable filechannel
    variable ActiveAppList
    global KPriv

    if {$KPriv(what_dempack_package) ne "S-DEMPack"} {
    set ref_props_counter 0
    }

    # For fluid application
    if {$AppId =="Fluid"} {
	# For fluid application

	GiD_File fprintf $filechannel "%s" "Begin Properties 0"
	GiD_File fprintf $filechannel "%s" "End Properties"
	GiD_File fprintf $filechannel ""

    } elseif {$AppId =="ConvectionDiffusion"} {
	# For convection-diffusion application
	# Kratos key word xpath
	set kxpath "Applications/$AppId"
	set cproperty "dv"

	# Only in the case of linear analysis
	# Analysis type
	set cxpath "$AppId//c.AnalysisData//i.AnalysisType"
	set AnalysisType [::xmlutils::setXml $cxpath $cproperty]
	if {$AnalysisType eq "Linear"} {

	    # Write the base properties
	    GiD_File fprintf $filechannel "%s" "Begin Properties 0"
	    # Write default values
	    set plist [list "Conductivity" "SpecificHeat" "VEmissivity"]
	    set cvalue "0.0"
	    foreach pid $plist {
		# Get the kratos keyword
		set ckword [::xmlutils::getKKWord $kxpath $pid]
		if {$ckword !=""} {
		    GiD_File fprintf $filechannel "%s" " $ckword $cvalue"
		}
	    }
	    GiD_File fprintf $filechannel "%s" "End Properties"
	    GiD_File fprintf $filechannel ""
	}
    } elseif {$AppId =="DEM"} {
	variable deminletchannel
	variable demfemchannel
	variable dem_props_number
	set basexpath "$AppId//c.DEM-Elements"
	set gelemlist [::wkcf::GetDEMActiveElements]
	set props_number 0
	foreach celemid $gelemlist {
	    set elembasexpath "$basexpath//c.DEM-Element"
	    # Get the group node list
	    set cgroupidlist [::xmlutils::setXmlContainerIds "$AppId//c.DEM-Elements//c.DEM-Element"]
	    foreach cgroupid $cgroupidlist {
		set using_dem_kdem 0
		set active_or_not [::xmlutils::setXml "$AppId//c.DEM-Elements//c.DEM-Element//c.$cgroupid//c.Properties//i.SetActive" "dv"]
		if {$active_or_not=="No"} {
		    continue
		}
		if {$cgroupid != ""} {
		    set elist [GiD_EntitiesGroups get $cgroupid elements]
		    if {[llength $elist]} {
		        # Write all nodes for this group in increasing orden
		        GiD_File fprintf $filechannel "%s" "Begin Properties $dem_props_number"
		        set matxpath "$elembasexpath//c.$cgroupid//c.Properties//i.Material"
		        set material [::xmlutils::setXml $matxpath "dv" ]

		        set cxpath "DEMMaterial//m.$material//p.Density"
		        set propvalue [::xmlutils::setXml $cxpath "dv" "read" "" "mat"]
		        GiD_File fprintf $filechannel "PARTICLE_DENSITY $propvalue"

		        set cxpath "DEMMaterial//m.$material//p.YoungModulus"
		        set propvalue [::xmlutils::setXml $cxpath "dv" "read" "" "mat"]
		        GiD_File fprintf $filechannel "YOUNG_MODULUS $propvalue"

		        set cxpath "DEMMaterial//m.$material//p.PoissonRatio"
		        set propvalue [::xmlutils::setXml $cxpath "dv" "read" "" "mat"]
		        GiD_File fprintf $filechannel "POISSON_RATIO $propvalue"

		        set cxpath "DEMMaterial//m.$material//p.ParticleFrictionAngle"
		        set propvalue [::xmlutils::setXml $cxpath "dv" "read" "" "mat"]
		        set pi 3.1415926535897931
		        set propvalue [expr {tan($propvalue*$pi/180.0)}]
		        GiD_File fprintf $filechannel "FRICTION $propvalue"

		        set cxpath "DEMMaterial//m.$material//p.ParticleCohesion"
		        set propvalue [::xmlutils::setXml $cxpath "dv" "read" "" "mat"]
		        GiD_File fprintf $filechannel "PARTICLE_COHESION $propvalue"

		        set cxpath "DEMMaterial//m.$material//p.CoefficientOfRestitution"
		        set propvalue [::xmlutils::setXml $cxpath "dv" "read" "" "mat"]
		        GiD_File fprintf $filechannel "COEFFICIENT_OF_RESTITUTION $propvalue"

		        set cxpath "DEMMaterial//m.$material//p.Color"
		        set propvalue [::xmlutils::setXml $cxpath "dv" "read" "" "mat"]
		        GiD_File fprintf $filechannel "PARTICLE_MATERIAL $propvalue"

		        set cxpath "DEMMaterial//m.$material//p.RollingFriction"
		        set propvalue [::xmlutils::setXml $cxpath "dv" "read" "" "mat"]
		        GiD_File fprintf $filechannel "ROLLING_FRICTION $propvalue"

		        set cxpath "DEMMaterial//m.$material//p.RollingFrictionWithWalls"
		        set propvalue [::xmlutils::setXml $cxpath "dv" "read" "" "mat"]
		        GiD_File fprintf $filechannel "ROLLING_FRICTION_WITH_WALLS $propvalue"

		        if {$KPriv(what_dempack_package) eq "F-DEMPack"} {
		            set cxpath "DEMMaterial//m.$material//p.ParticleSphericity"
		            set propvalue [::xmlutils::setXml $cxpath "dv" "read" "" "mat"]
		            GiD_File fprintf $filechannel "PARTICLE_SPHERICITY $propvalue"
		        }

		        #Constitutive Law names
		        if {$KPriv(what_dempack_package) eq "C-DEMPack"} {
		            set cxpath "DEMMaterial//m.$material//p.DEM_ConstitutiveLaw"
		            set propvalue [::xmlutils::setXml $cxpath "dv" "read" "" "mat"]
		            if {$ndime =="2D"} {
		                if {$propvalue == "DEMPack"} {
		                    GiD_File fprintf $filechannel "DEM_CONTINUUM_CONSTITUTIVE_LAW_NAME DEM_Dempack2D"
		                } elseif {$propvalue == "KDEM"} {
		                    set using_dem_kdem 1
		                    GiD_File fprintf $filechannel "DEM_CONTINUUM_CONSTITUTIVE_LAW_NAME DEM_KDEM2D"
		                } elseif {$propvalue == "KDEMFabric"} {
		                    set using_dem_kdem 1
		                    GiD_File fprintf $filechannel "FABRIC_COEFFICIENT [::xmlutils::setXml "DEMMaterial//m.$material//p.DEM_Fabric_Coefficient" dv read {} mat]"
		                    GiD_File fprintf $filechannel "DEM_CONTINUUM_CONSTITUTIVE_LAW_NAME DEM_KDEMFabric2D"
		                } else {
		                    WarnWin [= "Unknown Continuum Constitutive Law for material $material"]
		                }
		            } else {
		                if {$propvalue == "DEMPack"} {
		                    GiD_File fprintf $filechannel "DEM_CONTINUUM_CONSTITUTIVE_LAW_NAME DEM_Dempack"
		                } elseif {$propvalue == "KDEM"} {
		                    set using_dem_kdem 1
		                    GiD_File fprintf $filechannel "DEM_CONTINUUM_CONSTITUTIVE_LAW_NAME DEM_KDEM"
		        		} elseif {$propvalue == "KDEM_Rankine"} {
		                    set using_dem_kdem 1
		                    GiD_File fprintf $filechannel "DEM_CONTINUUM_CONSTITUTIVE_LAW_NAME DEM_KDEM_Rankine"
		                } elseif {$propvalue == "KDEM_Mohr_Coulomb"} {
		                    set using_dem_kdem 1
		                    GiD_File fprintf $filechannel "DEM_CONTINUUM_CONSTITUTIVE_LAW_NAME DEM_KDEM_Mohr_Coulomb"
		                    GiD_File fprintf $filechannel "INTERNAL_COHESION [::xmlutils::setXml "DEMMaterial//m.$material//p.DEM_internal_cohesion" dv read {} mat]"
		                    GiD_File fprintf $filechannel "INTERNAL_FRICTION_ANGLE [::xmlutils::setXml "DEMMaterial//m.$material//p.DEM_internal_friction_angle" dv read {} mat]"
						} elseif {$propvalue == "KDEM_Fissured_Rock"} {
		                    set using_dem_kdem 1
		                    GiD_File fprintf $filechannel "DEM_CONTINUUM_CONSTITUTIVE_LAW_NAME DEM_KDEM_Fissured_Rock"
		                    GiD_File fprintf $filechannel "TENSION_LIMIT_INCREASE_SLOPE [::xmlutils::setXml "DEMMaterial//m.$material//p.DEM_stress_limit_growth_slope" dv read {} mat]"
		                } elseif {$propvalue == "KDEMFabric"} {
		                    set using_dem_kdem 1
		                    GiD_File fprintf $filechannel "FABRIC_COEFFICIENT [::xmlutils::setXml "DEMMaterial//m.$material//p.DEM_Fabric_Coefficient" dv read {} mat]"
		                    GiD_File fprintf $filechannel "DEM_CONTINUUM_CONSTITUTIVE_LAW_NAME DEM_KDEMFabric"
		                } else {
		                    WarnWin [= "Unknown Continuum Constitutive Law for material $material"]
		                }
		            }
		        }
		        set cxpath "DEMMaterial//m.$material//p.DEM_ContactLaw"
		        set contact_law [::xmlutils::setXml $cxpath "dv" "read" "" "mat"]
		        if {$KPriv(what_dempack_package) ne "C-DEMPack"} {
		            set cxpath "DEMMaterial//m.$material//p.DEM_CohesionLaw"
		            set cohesion_law [::xmlutils::setXml $cxpath "dv" "read" "" "mat"]

		            if {($contact_law eq "Hertz") && ($cohesion_law eq "None")} {
		                GiD_File fprintf $filechannel "DEM_DISCONTINUUM_CONSTITUTIVE_LAW_NAME DEM_D_Hertz_viscous_Coulomb"
		            } elseif {($contact_law eq "Linear") && ($cohesion_law eq "None")} {
		                GiD_File fprintf $filechannel "DEM_DISCONTINUUM_CONSTITUTIVE_LAW_NAME DEM_D_Linear_viscous_Coulomb"
		            } elseif {($contact_law eq "Bentonite") && ($cohesion_law eq "None")} {
		                GiD_File fprintf $filechannel "DEM_DISCONTINUUM_CONSTITUTIVE_LAW_NAME DEM_D_Bentonite_Colloid"
		            } elseif {($contact_law eq "Hertz") && ($cohesion_law eq "JKR")} {
		                GiD_File fprintf $filechannel "DEM_DISCONTINUUM_CONSTITUTIVE_LAW_NAME DEM_D_Hertz_viscous_Coulomb_JKR"
		            } elseif {($contact_law eq "Linear") && ($cohesion_law eq "JKR")} {
		                GiD_File fprintf $filechannel "DEM_DISCONTINUUM_CONSTITUTIVE_LAW_NAME DEM_D_Linear_viscous_Coulomb_JKR"
		            } elseif {($contact_law eq "Hertz") && ($cohesion_law eq "DMT")} {
		                GiD_File fprintf $filechannel "DEM_DISCONTINUUM_CONSTITUTIVE_LAW_NAME DEM_D_Hertz_viscous_Coulomb_DMT"
		            } elseif {($contact_law eq "Linear") && ($cohesion_law eq "DMT")} {
		                GiD_File fprintf $filechannel "DEM_DISCONTINUUM_CONSTITUTIVE_LAW_NAME DEM_D_Linear_viscous_Coulomb_DMT"
		            } elseif {($contact_law eq "ConicalDamage") && ($cohesion_law eq "None")} {
		                GiD_File fprintf $filechannel "DEM_DISCONTINUUM_CONSTITUTIVE_LAW_NAME DEM_D_Conical_damage"
		            } elseif {$contact_law eq "LinearCustomized"} {
		                GiD_File fprintf $filechannel "K_NORMAL [::xmlutils::setXml DEMMaterial//m.$material//p.KNormal dv read {} mat]"
		                GiD_File fprintf $filechannel "K_TANGENTIAL [::xmlutils::setXml DEMMaterial//m.$material//p.KTangential dv read {} mat]"
		                GiD_File fprintf $filechannel "DEM_DISCONTINUUM_CONSTITUTIVE_LAW_NAME DEM_D_Linear_Custom_Constants"
		            } else {
		                GiD_File fprintf $filechannel "DEM_DISCONTINUUM_CONSTITUTIVE_LAW_NAME DEMDiscontinuumConstitutiveLaw"
		            }
		            GiD_File fprintf $filechannel "DEM_CONTINUUM_CONSTITUTIVE_LAW_NAME DEMContinuumConstitutiveLaw"
		        }
		        if {$KPriv(what_dempack_package) eq "C-DEMPack"} {
		            if {$contact_law eq "Hertz"} {
		                if {$ndime =="2D"} {
		                    GiD_File fprintf $filechannel "DEM_DISCONTINUUM_CONSTITUTIVE_LAW_NAME DEM_D_Hertz_viscous_Coulomb2D"
		                } else {
		                    GiD_File fprintf $filechannel "DEM_DISCONTINUUM_CONSTITUTIVE_LAW_NAME DEM_D_Hertz_viscous_Coulomb"
		                }
		            } elseif {$contact_law eq "Linear"} {
		                if {$ndime =="2D"} {
		                    GiD_File fprintf $filechannel "DEM_DISCONTINUUM_CONSTITUTIVE_LAW_NAME DEM_D_Linear_viscous_Coulomb2D"
		                } else {
		                    GiD_File fprintf $filechannel "DEM_DISCONTINUUM_CONSTITUTIVE_LAW_NAME DEM_D_Linear_viscous_Coulomb"
		                }

					} elseif {$contact_law eq "Linear_HighStiffness"} {
		                    GiD_File fprintf $filechannel "DEM_DISCONTINUUM_CONSTITUTIVE_LAW_NAME DEM_D_Linear_HighStiffness"

		            } elseif {$contact_law eq "LinearCustomized"} {
		                GiD_File fprintf $filechannel "K_NORMAL [::xmlutils::setXml DEMMaterial//m.$material//p.KNormal dv read {} mat]"
		                GiD_File fprintf $filechannel "K_TANGENTIAL [::xmlutils::setXml DEMMaterial//m.$material//p.KTangential dv read {} mat]"
		                GiD_File fprintf $filechannel "DEM_DISCONTINUUM_CONSTITUTIVE_LAW_NAME DEM_D_Linear_Custom_Constants"
		            } else {
		                GiD_File fprintf $filechannel "DEM_DISCONTINUUM_CONSTITUTIVE_LAW_NAME DEMDiscontinuumConstitutiveLaw"
		            }

		            set cxpath "DEM//c.DEM-Options//c.DEM-AdvancedOptions//i.DEM-Thermal"
		            set thermal_option [::xmlutils::setXml $cxpath "dv" "read" "" "props"]
		            if {$thermal_option eq "Yes"} {
		                set cxpath "DEMMaterial//m.$material//c.ThermalParameters//p.ThermalConductivity"
		                set propvalue [::xmlutils::setXml $cxpath "dv" "read" "" "mat"]
		                GiD_File fprintf $filechannel "THERMAL_CONDUCTIVITY $propvalue"
		                set cxpath "DEMMaterial//m.$material//c.ThermalParameters//p.SpecificHeat"
		                set propvalue [::xmlutils::setXml $cxpath "dv" "read" "" "mat"]
		                GiD_File fprintf $filechannel "SPECIFIC_HEAT $propvalue"
		            }
		            if {$using_dem_kdem == 0} {
		                set cxpath "DEMMaterial//m.$material//p.LCS1"
		                set propvalue [::xmlutils::setXml $cxpath "dv" "read" "" "mat"]
		                GiD_File fprintf $filechannel "SLOPE_LIMIT_COEFF_C1 $propvalue"
		                set cxpath "DEMMaterial//m.$material//p.LCS2"
		                set propvalue [::xmlutils::setXml $cxpath "dv" "read" "" "mat"]
		                GiD_File fprintf $filechannel "SLOPE_LIMIT_COEFF_C2 $propvalue"
		                set cxpath "DEMMaterial//m.$material//p.LCS3"
		                set propvalue [::xmlutils::setXml $cxpath "dv" "read" "" "mat"]
		                GiD_File fprintf $filechannel "SLOPE_LIMIT_COEFF_C3 $propvalue"
		                set cxpath "DEMMaterial//m.$material//p.YRC1"
		                set propvalue [::xmlutils::setXml $cxpath "dv" "read" "" "mat"]
		                GiD_File fprintf $filechannel "SLOPE_FRACTION_N1 $propvalue"
		                set cxpath "DEMMaterial//m.$material//p.YRC2"
		                set propvalue [::xmlutils::setXml $cxpath "dv" "read" "" "mat"]
		                GiD_File fprintf $filechannel "SLOPE_FRACTION_N2 $propvalue"
		                set cxpath "DEMMaterial//m.$material//p.YRC3"
		                set propvalue [::xmlutils::setXml $cxpath "dv" "read" "" "mat"]
		                GiD_File fprintf $filechannel "SLOPE_FRACTION_N3 $propvalue"
		                set cxpath "DEMMaterial//m.$material//p.PlasticYoungModulus"
		                set propvalue [::xmlutils::setXml $cxpath "dv" "read" "" "mat"]
		                GiD_File fprintf $filechannel "YOUNG_MODULUS_PLASTIC $propvalue"
		                set cxpath "DEMMaterial//m.$material//p.PlasticYieldStress"
		                set propvalue [::xmlutils::setXml $cxpath "dv" "read" "" "mat"]
		                GiD_File fprintf $filechannel "PLASTIC_YIELD_STRESS $propvalue"
		                set cxpath "DEMMaterial//m.$material//p.DamageDeformationFactor"
		                set propvalue [::xmlutils::setXml $cxpath "dv" "read" "" "mat"]
		                GiD_File fprintf $filechannel "DAMAGE_FACTOR $propvalue"
		                set cxpath "DEMMaterial//m.$material//p.ShearEnergyCoeff"
		                set propvalue [::xmlutils::setXml $cxpath "dv" "read" "" "mat"]
		                GiD_File fprintf $filechannel "SHEAR_ENERGY_COEF $propvalue"
		            }
		            set cxpath "DEMMaterial//m.$material//p.TangentialStrength"
		            set propvalue [::xmlutils::setXml $cxpath "dv" "read" "" "mat"]
		            GiD_File fprintf $filechannel "CONTACT_TAU_ZERO $propvalue"

		            set cxpath "DEMMaterial//m.$material//p.NormalTensileStrength"
		            set propvalue [::xmlutils::setXml $cxpath "dv" "read" "" "mat"]
		            GiD_File fprintf $filechannel "CONTACT_SIGMA_MIN $propvalue"

		            set cxpath "DEMMaterial//m.$material//p.InternalFrictionAngleCoeff"
		            set propvalue [::xmlutils::setXml $cxpath "dv" "read" "" "mat"]
		            GiD_File fprintf $filechannel "CONTACT_INTERNAL_FRICC $propvalue"
		        }
		        GiD_File fprintf $filechannel "%s" "End Properties"
		        GiD_File fprintf $filechannel "%s" ""
		    }
		}
		incr dem_props_number
	    }
	}

	set rootid $AppId
	set basexpath "$rootid//c.DEM-Conditions//c.DEM-FEM-Wall"
	set gproplist [::xmlutils::setXmlContainerIds $basexpath]
	foreach cgroupid $gproplist {
	    variable demwall_element_props
	    set demwall_element_props($cgroupid) $ref_props_counter

	    set cproperty "dv"
	    set cxpath "$rootid//c.DEM-Conditions//c.DEM-FEM-Wall//c.[list ${cgroupid}]//i.SetActive"
	    set active_or_not [::xmlutils::setXml $cxpath $cproperty]
	    if {$active_or_not=="No"} {
		  continue
	    }
	    GiD_File fprintf $demfemchannel "%s" "Begin Properties $ref_props_counter"
	    set cproperty "dv"
	    set options_container Options
	    if {"Fluid" in $ActiveAppList} {
		set options_container [GetOptionsContainer]
	    }
	    set cxpath "$rootid//c.DEM-Conditions//c.DEM-FEM-Wall//c.[list ${cgroupid}]//c.$options_container//i.friction_coeff"
	    set friction_value [::xmlutils::setXml $cxpath $cproperty]
	    set pi 3.1415926535897931
	    set friction_value [expr {tan($friction_value*$pi/180.0)}]
	    GiD_File fprintf $demfemchannel "FRICTION $friction_value"
	    set cxpath "$rootid//c.DEM-Conditions//c.DEM-FEM-Wall//c.[list ${cgroupid}]//c.$options_container//i.WallCohesion"
	    set cohesive_wall [::xmlutils::setXml $cxpath $cproperty]
	    GiD_File fprintf $demfemchannel "WALL_COHESION $cohesive_wall"
	    set cxpath "$rootid//c.DEM-Conditions//c.DEM-FEM-Wall//c.[list ${cgroupid}]//c.$options_container//i.DEM_Wear"
	    set compute_wear [::xmlutils::setXml $cxpath $cproperty]
	    if {$compute_wear == "Yes"} {
		set compute_wear 1
	    } else {
		set compute_wear 0
	    }
	    GiD_File fprintf $demfemchannel "COMPUTE_WEAR $compute_wear"
	    set cxpath "$rootid//c.DEM-Conditions//c.DEM-FEM-Wall//c.[list ${cgroupid}]//c.$options_container//i.K_Abrasion"
	    set archard_value [::xmlutils::setXml $cxpath $cproperty]
	    GiD_File fprintf $demfemchannel "SEVERITY_OF_WEAR $archard_value"
	    set cxpath "$rootid//c.DEM-Conditions//c.DEM-FEM-Wall//c.[list ${cgroupid}]//c.$options_container//i.K_Impact"
	    set k_impact [::xmlutils::setXml $cxpath $cproperty]
	    GiD_File fprintf $demfemchannel "IMPACT_WEAR_SEVERITY $k_impact"
	    set cxpath "$rootid//c.DEM-Conditions//c.DEM-FEM-Wall//c.[list ${cgroupid}]//c.$options_container//i.H_Brinell"
	    set brinell_value [::xmlutils::setXml $cxpath $cproperty]
	    GiD_File fprintf $demfemchannel "BRINELL_HARDNESS $brinell_value"
	    set cxpath "$rootid//c.DEM-Conditions//c.DEM-FEM-Wall//c.[list ${cgroupid}]//c.$options_container//i.YoungModulus"
	    set young_modulus [::xmlutils::setXml $cxpath $cproperty]
	    set path_to_rigid_or_not "$rootid//c.DEM-Conditions//c.DEM-FEM-Wall//c.[list ${cgroupid}]//c.$options_container//i.RigidPlane"
	    set rigid_or_not [::xmlutils::setXml $path_to_rigid_or_not $cproperty]
	    if {$rigid_or_not == "Yes"} {
		  set young_modulus 1e20
	    }
	    GiD_File fprintf $demfemchannel "YOUNG_MODULUS $young_modulus"
	    set cxpath "$rootid//c.DEM-Conditions//c.DEM-FEM-Wall//c.[list ${cgroupid}]//c.$options_container//i.PoissonRatio"
	    set poisson_ratio [::xmlutils::setXml $cxpath $cproperty]
	    GiD_File fprintf $demfemchannel "POISSON_RATIO $poisson_ratio"
	    GiD_File fprintf $demfemchannel "%s" "End Properties"
	    GiD_File fprintf $demfemchannel "%s" ""
	    incr ref_props_counter
	}
	set kxpath "Applications/$AppId"
	set rootid $AppId
	set basexpath "$rootid//c.DEM-Conditions//c.DEM-Inlet"
	set gproplist [::xmlutils::setXmlContainerIds $basexpath]
	set inlet_element_type SphericParticle3D
	variable ActiveAppList
	if {$KPriv(what_dempack_package) eq "C-DEMPack"} {
	    set inlet_element_type SphericContinuumParticle3D
	}
	foreach cgroupid $gproplist {
	    set properties_path "$basexpath//c.$cgroupid//c.MainProperties"
	    set material [::xmlutils::setXml "$properties_path//i.Material" "dv"]
	    set constitutive_law ""
	    set contains_clusters 0
	    set random_euler_angles 0
	    set using_dem_kdem 0
	    set active_or_not [::xmlutils::setXml "$properties_path//i.SetActive" "dv"]
	    if {$active_or_not=="No"} {
	    continue
	    }
	    GiD_File fprintf $deminletchannel "%s" "Begin Properties $dem_props_number // Inlet group: $cgroupid"
	    ##MATERIAL
	set density [::xmlutils::setXml "DEMMaterial//m.$material//p.Density" "dv" "read" "" "mat"]
	GiD_File fprintf $deminletchannel "  PARTICLE_DENSITY $density"
	set young_modulus [::xmlutils::setXml "DEMMaterial//m.$material//p.YoungModulus" "dv" "read" "" "mat"]
	GiD_File fprintf $deminletchannel "  YOUNG_MODULUS $young_modulus"
	set poisson_ratio [::xmlutils::setXml "DEMMaterial//m.$material//p.PoissonRatio" "dv" "read" "" "mat"]
	GiD_File fprintf $deminletchannel "  POISSON_RATIO $poisson_ratio"
	set friction_angle [::xmlutils::setXml "DEMMaterial//m.$material//p.ParticleFrictionAngle" "dv" "read" "" "mat"]
	set friction_angle [expr {tan(3.1415926535897931*$friction_angle/180.0)}]
	GiD_File fprintf $deminletchannel "  FRICTION $friction_angle"
	set cohesion [::xmlutils::setXml "DEMMaterial//m.$material//p.ParticleCohesion" "dv" "read" "" "mat"]
	GiD_File fprintf $deminletchannel "  PARTICLE_COHESION $cohesion"
	set coeff_of_rest [::xmlutils::setXml "DEMMaterial//m.$material//p.CoefficientOfRestitution" "dv" "read" "" "mat"]
	GiD_File fprintf $deminletchannel "  COEFFICIENT_OF_RESTITUTION $coeff_of_rest"
	set particle_material [::xmlutils::setXml "DEMMaterial//m.$material//p.Color" "dv" "read" "" "mat"]
	GiD_File fprintf $deminletchannel "  PARTICLE_MATERIAL $particle_material"
	set rolling_friction [::xmlutils::setXml "DEMMaterial//m.$material//p.RollingFriction" "dv" "read" "" "mat"]
	GiD_File fprintf $deminletchannel "  ROLLING_FRICTION $rolling_friction"
	set rolling_friction_with_walls [::xmlutils::setXml "DEMMaterial//m.$material//p.RollingFrictionWithWalls" "dv" "read" "" "mat"]
	GiD_File fprintf $deminletchannel "  ROLLING_FRICTION_WITH_WALLS $rolling_friction_with_walls"
	set sphericity [::xmlutils::setXml "DEMMaterial//m.$material//p.ParticleSphericity" "dv" "read" "" "mat"]
	GiD_File fprintf $deminletchannel "  PARTICLE_SPHERICITY $sphericity"
	set contact_law [::xmlutils::setXml "DEMMaterial//m.$material//p.DEM_ContactLaw" "dv" "read" "" "mat"]
	if {$contact_law eq "Hertz"} {
	    GiD_File fprintf $deminletchannel "  DEM_DISCONTINUUM_CONSTITUTIVE_LAW_NAME DEM_D_Hertz_viscous_Coulomb"
	} elseif {$contact_law eq "Linear"} {
	    GiD_File fprintf $deminletchannel "  DEM_DISCONTINUUM_CONSTITUTIVE_LAW_NAME DEM_D_Linear_viscous_Coulomb"
	} elseif {$contact_law eq "LinearCustomized"} {
	    GiD_File fprintf $deminletchannel "  K_NORMAL [::xmlutils::setXml DEMMaterial//m.$material//p.KNormal dv read {} mat]"
	    GiD_File fprintf $deminletchannel "  K_TANGENTIAL [::xmlutils::setXml DEMMaterial//m.$material//p.KTangential dv read {} mat]"
	    GiD_File fprintf $deminletchannel "  DEM_DISCONTINUUM_CONSTITUTIVE_LAW_NAME DEM_D_Linear_Custom_Constants"
	} else {
	    GiD_File fprintf $deminletchannel "  DEM_DISCONTINUUM_CONSTITUTIVE_LAW_NAME DEMDiscontinuumConstitutiveLaw"
	    W "No Discontinuum Law was chosen!"
	}
	if {$KPriv(what_dempack_package) eq "C-DEMPack"} {
	    set constitutive_law [::xmlutils::setXml "DEMMaterial//m.$material//p.DEM_ConstitutiveLaw" "dv" "read" "" "mat"]
	}
	if {$ndime =="2D"} {
	    GiD_File fprintf $deminletchannel "  DEM_CONTINUUM_CONSTITUTIVE_LAW_NAME DEM_Dempack2D"
	} elseif {$constitutive_law == "DEMPack"} {
	    GiD_File fprintf $deminletchannel "  DEM_CONTINUUM_CONSTITUTIVE_LAW_NAME DEM_Dempack"
	} elseif {$constitutive_law == "KDEM"} {
	    GiD_File fprintf $deminletchannel "  DEM_CONTINUUM_CONSTITUTIVE_LAW_NAME DEM_KDEM"
	set using_dem_kdem 1
	} else {
	    GiD_File fprintf $deminletchannel "  DEM_CONTINUUM_CONSTITUTIVE_LAW_NAME DEMContinuumConstitutiveLaw"
	    if {$KPriv(what_dempack_package) eq "C-DEMPack"} {
	    W "No Continuum Law was chosen!"
	    }
	}
	if {$KPriv(what_dempack_package) eq "C-DEMPack"} {
		set cohesive [::xmlutils::setXml "$properties_path//i.SetBreakable" "dv"]
		if {$cohesive == "Yes"} {
	    GiD_File fprintf $deminletchannel "  BREAKABLE_CLUSTER 1"
		} else {
	    GiD_File fprintf $deminletchannel "  BREAKABLE_CLUSTER 0"
		}
	    if {$using_dem_kdem == 0} {
	    set slope_limit_coeff_1 [::xmlutils::setXml "DEMMaterial//m.$material//p.LCS1" "dv" "read" "" "mat"]
	    GiD_File fprintf $deminletchannel "  SLOPE_LIMIT_COEFF_C1 $slope_limit_coeff_1"
	    set slope_limit_coeff_2 [::xmlutils::setXml "DEMMaterial//m.$material//p.LCS2" "dv" "read" "" "mat"]
	    GiD_File fprintf $deminletchannel "  SLOPE_LIMIT_COEFF_C2 $slope_limit_coeff_2"
	    set slope_limit_coeff_3 [::xmlutils::setXml "DEMMaterial//m.$material//p.LCS3" "dv" "read" "" "mat"]
	    GiD_File fprintf $deminletchannel "  SLOPE_LIMIT_COEFF_C3 $slope_limit_coeff_3"
	    set slope_fraction_1 [::xmlutils::setXml "DEMMaterial//m.$material//p.YRC1" "dv" "read" "" "mat"]
	    GiD_File fprintf $deminletchannel "  SLOPE_FRACTION_N1 $slope_fraction_1"
	    set slope_fraction_2 [::xmlutils::setXml "DEMMaterial//m.$material//p.YRC2" "dv" "read" "" "mat"]
	    GiD_File fprintf $deminletchannel "  SLOPE_FRACTION_N2 $slope_fraction_2"
	    set slope_fraction_3 [::xmlutils::setXml "DEMMaterial//m.$material//p.YRC3" "dv" "read" "" "mat"]
	    GiD_File fprintf $deminletchannel "  SLOPE_FRACTION_N3 $slope_fraction_3"
	    set plastic_young_modulus [::xmlutils::setXml "DEMMaterial//m.$material//p.PlasticYoungModulus" "dv" "read" "" "mat"]
	    GiD_File fprintf $deminletchannel "  YOUNG_MODULUS_PLASTIC $plastic_young_modulus"
	    set plastic_yield_stress [::xmlutils::setXml "DEMMaterial//m.$material//p.PlasticYieldStress" "dv" "read" "" "mat"]
	    GiD_File fprintf $deminletchannel "  PLASTIC_YIELD_STRESS $plastic_yield_stress"
	    set damage_factor [::xmlutils::setXml "DEMMaterial//m.$material//p.DamageDeformationFactor" "dv" "read" "" "mat"]
	    GiD_File fprintf $deminletchannel "  DAMAGE_FACTOR $damage_factor"
	    set shear_energy_coeff [::xmlutils::setXml "DEMMaterial//m.$material//p.ShearEnergyCoeff" "dv" "read" "" "mat"]
	    GiD_File fprintf $deminletchannel "  SHEAR_ENERGY_COEF $shear_energy_coeff"
	    }
	    set contact_tau_zero [::xmlutils::setXml "DEMMaterial//m.$material//p.TangentialStrength" "dv" "read" "" "mat"]
	    GiD_File fprintf $deminletchannel "  CONTACT_TAU_ZERO $contact_tau_zero"
	    set contact_sigma_min [::xmlutils::setXml "DEMMaterial//m.$material//p.NormalTensileStrength" "dv" "read" "" "mat"]
	    GiD_File fprintf $deminletchannel "  CONTACT_SIGMA_MIN $contact_sigma_min"
	    set internal_friction [::xmlutils::setXml "DEMMaterial//m.$material//p.InternalFrictionAngleCoeff" "dv" "read" "" "mat"]
	    GiD_File fprintf $deminletchannel "  CONTACT_INTERNAL_FRICC $internal_friction"
	}

	if {"Fluid" in $ActiveAppList} {
	    set list_of_cluster_files [list]
	} else {
	    if {[::xmlutils::setXml "${properties_path}//i.InletElementType" "dv"] eq "Cluster3D"} {
	    set inlet_element_type [::xmlutils::setXml "${properties_path}//i.Cluster3D" dv]
	    set contains_clusters 1
	    lassign [::wkcf::GetClusterFileNameAndReplaceInletElementType $inlet_element_type] inlet_element_type cluster_file_name
	    }
	    if {$inlet_element_type eq "Cluster3D"} {
	    lappend KPriv(list_of_cluster_files) $cluster_file_name
	    GiD_File fprintf $deminletchannel "  CLUSTER_FILE_NAME $cluster_file_name"
	    }
	}

	GiD_File fprintf $deminletchannel "%s" "End Properties"
	GiD_File fprintf $deminletchannel ""
	incr dem_props_number
    }

	set KPriv(list_of_cluster_files) [lsort -unique $KPriv(list_of_cluster_files)]

    }
    return [list 0 ""]
}

proc ::wkcf::WriteDSOLIDProperties {AppId} {
  variable demfemchannel
  global ref_props_counter
  set ref_props_counter 0

  set rootid DSOLID
  set basexpath "$rootid//c.Solid-Elements//c.Solid-Element"
  set gproplist [::xmlutils::setXmlContainerIds $basexpath]

  foreach cgroupid $gproplist {
    global solid_element_props
    set solid_element_props($cgroupid) $ref_props_counter

	set cproperty "dv"
	set cxpath "$rootid//c.Solid-Elements//c.Solid-Element//c.[list ${cgroupid}]//c.Properties//i.SetActive"
	set active_or_not [::xmlutils::setXml $cxpath $cproperty]
	if {$active_or_not=="No"} {
	    continue
	}

	GiD_File fprintf $demfemchannel "%s" "Begin Properties $ref_props_counter  // Write DSOLID Properties"

	set cproperty "dv"
	set cxpath "SolidMaterial//m.Solid-DefaultMaterial//c.General//p.Density"
	set propvalue [::xmlutils::setXml $cxpath "dv" "read" "" "mat"]
	GiD_File fprintf $demfemchannel "DENSITY $propvalue"

	set cxpath "SolidMaterial//m.Solid-DefaultMaterial//c.Structural//c.Elastic//c.Isotropic//p.YoungModulus"
	set propvalue [::xmlutils::setXml $cxpath "dv" "read" "" "mat"]
	GiD_File fprintf $demfemchannel "YOUNG_MODULUS $propvalue"

	set cxpath "SolidMaterial//m.Solid-DefaultMaterial//c.Structural//c.Elastic//c.Isotropic//p.PoissonRatio"
	set propvalue [::xmlutils::setXml $cxpath "dv" "read" "" "mat"]
	GiD_File fprintf $demfemchannel "POISSON_RATIO $propvalue"

	set cxpath "SolidMaterial//m.Solid-DefaultMaterial//c.Structural//c.ElastoPlastic//p.YieldStress"
	set propvalue [::xmlutils::setXml $cxpath "dv" "read" "" "mat"]
	GiD_File fprintf $demfemchannel "YIELD_STRESS $propvalue"

	set cxpath "SolidMaterial//m.Solid-DefaultMaterial//c.Structural//c.ElastoPlastic//p.IsotropicHardeningModulus"
	set propvalue [::xmlutils::setXml $cxpath "dv" "read" "" "mat"]
	GiD_File fprintf $demfemchannel "ISOTROPIC_HARDENING_MODULUS $propvalue"
	GiD_File fprintf $demfemchannel "%s" "BODY_FORCE \[3\] (0.0, 0.0, 0.0\)"
	GiD_File fprintf $demfemchannel "%s" "RAYLEIGH_ALPHA 0.0"
	GiD_File fprintf $demfemchannel "%s" "RAYLEIGH_BETA  0.0"
	GiD_File fprintf $demfemchannel "%s" "End Properties"
	GiD_File fprintf $demfemchannel "%s" ""
	incr ref_props_counter
    }
    return [list 0 ""]
}

proc ::wkcf::WriteDSOLIDContactProperties {AppId} {
  variable demfemchannel
  global ref_props_counter

  set rootid DSOLID
  set basexpath "$rootid//c.DEM-Conditions//c.DEM-FEM-Solid-Wall"
  set gproplist [::xmlutils::setXmlContainerIds $basexpath]

  foreach cgroupid $gproplist {
    variable solidcontact_element_props
    set solidcontact_element_props($cgroupid) $ref_props_counter

	set cproperty "dv"
	set cxpath "$rootid//c.DEM-Conditions//c.DEM-FEM-Solid-Wall//c.[list ${cgroupid}]//i.SetActive"
	set active_or_not [::xmlutils::setXml $cxpath $cproperty]
	if {$active_or_not=="No"} {
	    continue
	}

	GiD_File fprintf $demfemchannel "%s" "Begin Properties $ref_props_counter // Write DSOLID-Contact Properties"
	set cxpath "$rootid//c.DEM-Conditions//c.DEM-FEM-Solid-Wall//c.[list ${cgroupid}]//i.friction_coeff"
	set friction_value [::xmlutils::setXml $cxpath $cproperty]
	set cxpath "$rootid//c.DEM-Conditions//c.DEM-FEM-Solid-Wall//c.[list ${cgroupid}]//i.young_modulus"
	set young_modulus [::xmlutils::setXml $cxpath $cproperty]
	set cxpath "$rootid//c.DEM-Conditions//c.DEM-FEM-Solid-Wall//c.[list ${cgroupid}]//i.poisson_ratio"
	set poisson_ratio [::xmlutils::setXml $cxpath $cproperty]
	set pi 3.1415926535897931
	set friction_value [expr {tan($friction_value*$pi/180.0)}]
	GiD_File fprintf $demfemchannel "FRICTION $friction_value"
	GiD_File fprintf $demfemchannel "WALL_COHESION 0.0"
	GiD_File fprintf $demfemchannel "COMPUTE_WEAR 0.0"
	GiD_File fprintf $demfemchannel "SEVERITY_OF_WEAR 0.0"
	GiD_File fprintf $demfemchannel "IMPACT_WEAR_SEVERITY 0.0"
	GiD_File fprintf $demfemchannel "BRINELL_HARDNESS 0.0"
  GiD_File fprintf $demfemchannel "YOUNG_MODULUS $young_modulus"
	GiD_File fprintf $demfemchannel "POISSON_RATIO $poisson_ratio"
	GiD_File fprintf $demfemchannel "%s" "End Properties"
	GiD_File fprintf $demfemchannel "%s" ""
	incr ref_props_counter
    }
    return [list 0 ""]
}


proc ::wkcf::WriteDSOLIDNodalCoordinates {AppId} {
  variable demfemchannel

	set rootid DSOLID
	set basexpath "$rootid//c.Solid-Elements//c.Solid-Element"
	set gproplist [::xmlutils::setXmlContainerIds $basexpath]
	foreach cgroupid $gproplist {
    set cproperty "dv"
    set cxpath "$rootid//c.Solid-Elements//c.Solid-Element//c.[list ${cgroupid}]//c.Properties//i.SetActive"
    set active_or_not [::xmlutils::setXml $cxpath $cproperty]
    if {$active_or_not=="No"} {
	  continue
    }
    # Get the group node list
    set nlist [list]
    lassign [::wkcf::GetDSOLIDGroupNodes $AppId $cgroupid] fail msg nlist
    if { $fail == 1 } {
	  return [list 1 $msg]
    }
    if {[llength $nlist]} {
		# Write all nodes for this group in increasing orden
		GiD_File fprintf $demfemchannel "Begin Nodes // Write DSOLID Nodes"
		foreach nodeid $nlist {
		    lassign [GiD_Mesh get node $nodeid] layer x y z
		    GiD_File fprintf $demfemchannel "$nodeid [format {%#.6g %#.6g %#.6g} $x $y $z]"
		}
		GiD_File fprintf $demfemchannel "End Nodes"
		GiD_File fprintf $demfemchannel ""
	    }
	}
  return [list 0 ""]
}


proc ::wkcf::WriteNodalCoordinates {AppId} {
    # Write the nodal coordinates block
    # Nodes block format
    # Begin Nodes
    # // id          X        Y        Z
    # End Nodes
    # Arguments
    # AppId  => Application identifier
    variable ndime; variable dprops
    variable filechannel

    # For debug
    if {!$::wkcf::pflag} {
	set inittime [clock seconds]
    }

    if {$AppId ne "DEM"} {
	# Check for all defined kratos elements
	if {([info exists dprops($AppId,AllKElemId)]) && ([llength $dprops($AppId,AllKElemId)])} {

	    # For all defined kratos nodes
	    foreach celemid $dprops($AppId,AllKElemId) {

		# Check for all defined group identifier for this element
		if {[info exists dprops($AppId,KElem,$celemid,AllGroupId)] && [llength $dprops($AppId,KElem,$celemid,AllGroupId)]} {
		    # For all defined group identifier for this element
		    foreach cgroupid $dprops($AppId,KElem,$celemid,AllGroupId) {
		        # Get the GiD entity type, element type and property identifier
		        lassign $dprops($AppId,KElem,$celemid,$cgroupid,GProps) GiDEntity GiDElemType PropertyId KEKWord nDim
		        if { [GiD_EntitiesGroups get $cgroupid nodes -count] } {
		            # Write all nodes for this group in increasing orden
		            GiD_File fprintf $filechannel "Begin Nodes // GUI group identifier: $cgroupid celemid $celemid"
		            if {$ndime == "2D"} {
		                foreach nodeid [GiD_EntitiesGroups get $cgroupid nodes] {
		                    lassign [GiD_Mesh get node $nodeid] layer x y
		                    GiD_File fprintf $filechannel "$nodeid [format {%12g %12g} $x $y] 0"
		                }
		            } else {
		                foreach nodeid [GiD_EntitiesGroups get $cgroupid nodes] {
		                    lassign [GiD_Mesh get node $nodeid] layer x y z
		                    GiD_File fprintf $filechannel "$nodeid [format {%12g %12g %12g} $x $y $z]"
		                }
		            }
		            GiD_File fprintf $filechannel "End Nodes"
		            GiD_File fprintf $filechannel ""
		        }
		    }
		}
	    }
	}
    }
    # Special case of DEM application
    if {$AppId eq "DEM"} {

	variable deminletchannel
	variable demfemchannel

	# General MDPA
	set rootid $AppId
	set basexpath "$rootid//c.DEM-Elements//c.DEM-Element"
	set gproplist [::xmlutils::setXmlContainerIds $basexpath]
	set celemid [::wkcf::GetDEMActiveElements]
	foreach cgroupid $gproplist {

	    set cproperty "dv"
	    set cxpath "${basexpath}//c.[list ${cgroupid}]//c.Properties//i.SetActive"
	    set active_or_not [::xmlutils::setXml $cxpath $cproperty]

	    if {$active_or_not=="No"} {
		continue
	    }

	    # Get the group node list
	    set elist [GiD_EntitiesGroups get $cgroupid elements]
	    if {[llength $elist]} {

		# Write all nodes for this group in increasing order
		GiD_File fprintf $filechannel "Begin Nodes // GUI group identifier: $cgroupid celemid $celemid"
		if {$ndime == "2D"} {
		    foreach eid $elist {
		        # TODO: JG cambiar a GiD_Mesh
		        set nodeid [lindex [GiD_Info Mesh Elements circle $eid $eid] 1]
		        if {$nodeid==""} {
		            return [list 0 "An element which is not a 2D circle was found in group $cgroupid! Please remove all entities that are not circles from this group"]
		        }
		        lassign [GiD_Mesh get node $nodeid] layer x y

		        GiD_File fprintf $filechannel "$nodeid [format {%#.19g %#.19g} $x $y] 0"
		    }
		} else {
		    set nodeslist [list]
		    foreach eid $elist {
		        # TODO: JG cambiar a GiD_Mesh
		        set nodeid [lindex [GiD_Info Mesh Elements sphere $eid $eid] 1]
		        if {$nodeid==""} {
		            return [list 1 "An element which is not a 3D sphere was found in group $cgroupid! Please remove all entities that are not spheres from this group"]
		        }
		        lappend nodeslist $nodeid
		    }
		    set nodeslist [lsort -integer -unique $nodeslist]
		    foreach nid $nodeslist {
		        lassign [GiD_Mesh get node $nid] layer x y z
		        GiD_File fprintf $filechannel "$nid [format {%#.19g %#.19g %#.19g} $x $y $z]"
		    }
		}
		GiD_File fprintf $filechannel "End Nodes"
		GiD_File fprintf $filechannel ""
	    }
	}

	# Inlet MDPA
	# Get the values
	set rootid $AppId
	set basexpath "$rootid//c.DEM-Conditions//c.DEM-Inlet"
	set gproplist [::xmlutils::setXmlContainerIds $basexpath]
	foreach cgroupid $gproplist {
	    #
	    set cproperty "dv"
	    set cxpath "${basexpath}//c.[list ${cgroupid}]//c.MainProperties//i.SetActive"
	    set active_or_not [::xmlutils::setXml $cxpath $cproperty]
	    if {$active_or_not=="No"} {
		continue
	    }
	    # Get the group node list
	    set nlist [::wkcf::GetInletGroupNodes $AppId $cgroupid]
	    if {[llength $nlist]} {
		# Write all nodes for this group in increasing orden
		GiD_File fprintf $deminletchannel "Begin Nodes // GUI Inlet group identifier: $cgroupid"
		if {$ndime == "2D"} {
		    foreach nodeid $nlist {
		        lassign [GiD_Mesh get node $nodeid] layer x y
		        GiD_File fprintf $deminletchannel "$nodeid [format {%#.6g %#.6g} $x $y] 0"
		    }
		} else {
		    foreach nodeid $nlist {
		        lassign [GiD_Mesh get node $nodeid] layer x y z
		        GiD_File fprintf $deminletchannel "$nodeid [format {%#.6g %#.6g %#.6g} $x $y $z]"
		    }
		}
		GiD_File fprintf $deminletchannel "End Nodes"
		GiD_File fprintf $deminletchannel ""
	    }
	}

	# DEM FEM MDPA
	set rootid $AppId
	set basexpath "$rootid//c.DEM-Conditions//c.DEM-FEM-Wall"
	set gproplist [::xmlutils::setXmlContainerIds $basexpath]

	foreach cgroupid $gproplist {
	    #
	    set cproperty "dv"
	    set cxpath "$rootid//c.DEM-Conditions//c.DEM-FEM-Wall//c.[list ${cgroupid}]//i.SetActive"
	    set active_or_not [::xmlutils::setXml $cxpath $cproperty]
	    if {$active_or_not=="No"} {
		continue
	    }
	    # Get the group node list
	    set nlist [list]
	    lassign [::wkcf::GetDemFemWallGroupNodes $cgroupid] fail msg nlist
	    if { $fail == 1 } {
		return [list 1 $msg]
	    }
	    #set nlist [::wkcf::GetDemFemWallGroupNodes $AppId $cgroupid]
	    if {[llength $nlist]} {
		# Write all nodes for this group in increasing orden
		GiD_File fprintf $demfemchannel "Begin Nodes // GUI DEM-FEM-Wall group identifier: $cgroupid"
		if {$ndime == "2D"} {
		    foreach nodeid $nlist {
		        lassign [GiD_Mesh get node $nodeid] layer x y z
		        GiD_File fprintf $demfemchannel "$nodeid [format {%#.6g %#.6g} $x $y] 0"
		    }
		} else {
		    foreach nodeid $nlist {
		        lassign [GiD_Mesh get node $nodeid] layer x y z
		        GiD_File fprintf $demfemchannel "$nodeid [format {%#.6g %#.6g %#.6g} $x $y $z]"
		    }
		}
		GiD_File fprintf $demfemchannel "End Nodes"
		GiD_File fprintf $demfemchannel ""
	    }
	}
    }
    # For debug
    if {!$::wkcf::pflag} {
	set endtime [clock seconds]
	set ttime [expr $endtime-$inittime]
	WarnWinText "Write nodal coordinates: [::KUtils::Duration $ttime]"
    }
    return [list 0 ""]
}
proc ::wkcf::WriteDSOLIDElementConnectivities {AppId} {
  # Write the element connectivities block
  global solid_element_props

  # For debug
  if {!$::wkcf::pflag} {
  set inittime [clock seconds]
  }

  if {$AppId eq "DEM"} {
	variable demfemchannel
	set rootid DSOLID
	#set basexpath "$rootid//c.DEM-Conditions//c.DEM-FEM-Solid-Wall"
	set basexpath "$rootid//c.Solid-Elements//c.Solid-Element"
	set gproplist [::xmlutils::setXmlContainerIds $basexpath]

	foreach cgroupid $gproplist {
	    set cproperty "dv"
	    set cxpath "${basexpath}//c.[list ${cgroupid}]//c.Properties//i.SetActive"
	    set active_or_not [::xmlutils::setXml $cxpath $cproperty]
	    if {$active_or_not=="No"} {
		  continue
	    }
	    # Get the group node list
	    set nlist [GiD_EntitiesGroups get $cgroupid elements]
	    if {[llength $nlist]} {
		set elemname [::wkcf::GetDSOLIDElementName [lindex $nlist 0] ]
		GiD_File fprintf $demfemchannel "Begin Elements $elemname // write DSOLID elements"
		foreach elemid $nlist {
		    set elemtype [lindex [GiD_Mesh get element $elemid] 1]
		    if { ($elemtype eq "Tetrahedra") || ($elemtype eq "Hexahedra") } {
		        GiD_File fprintf $demfemchannel "$elemid $solid_element_props($cgroupid) [lrange [GiD_Mesh get element $elemid] 3 end]"
		    }
		}
		GiD_File fprintf $demfemchannel "End Elements"
		GiD_File fprintf $demfemchannel ""
	  }
	  }
    }
    # For debug
    if {!$::wkcf::pflag} {
	    set endtime [clock seconds]
	    set ttime [expr $endtime-$inittime]
	    # WarnWinText "endtime:$endtime ttime:$ttime"
	    WarnWinText "Write element connectivities: [::KUtils::Duration $ttime]"
    }
}


proc ::wkcf::WriteDSOLIDContactElements {AppId} {
  variable demfemchannel
  variable solidcontact_element_props

  set rootid DSOLID
  set basexpath "$rootid//c.DEM-Conditions//c.DEM-FEM-Solid-Wall"
  set gproplist [::xmlutils::setXmlContainerIds $basexpath]
  set number_of_solidcontact 3
  foreach cgroupid $gproplist {
	set cproperty "dv"
	set cxpath "$rootid//c.DEM-Conditions//c.DEM-FEM-Solid-Wall//c.[list ${cgroupid}]//i.SetActive"
	set active_or_not [::xmlutils::setXml $cxpath $cproperty]
	if {$active_or_not=="No"} {
	    continue
	}
	# Get the group node list
	set nlist [GiD_EntitiesGroups get $cgroupid elements]
	if {[llength $nlist]} {
	    set elemname [::wkcf::GetDSOLIDContactElementName [lindex $nlist 0] ]
	    GiD_File fprintf $demfemchannel "Begin Conditions $elemname // write DSOLID Contact Conditions"
	    foreach elemid $nlist {
		  GiD_File fprintf $demfemchannel "$elemid $solidcontact_element_props($cgroupid) [lrange [GiD_Mesh get element $elemid] 3 end]"
	    }
	    GiD_File fprintf $demfemchannel "End Conditions"
	    GiD_File fprintf $demfemchannel ""
	  }
    }
}

proc ::wkcf::WriteElementConnectivities {AppId} {
    # Write the element connectivities block
    # Element block format
    # Begin Elements element_name
    # // id prop_id         n1        n2        n3        ...
    # End Elements
    # Arguments
    # AppId => Application identifier
    variable useqelem; variable dprops
    variable filechannel
    variable FluidApplication
    variable ActiveAppList

    # For debug
    if {!$::wkcf::pflag} {
	set inittime [clock seconds]
    }

    # General case for NON DEM application
    if {$AppId ne "DEM"} {
	# Check for all defined kratos elements
	if {([info exists dprops($AppId,AllKElemId)]) && ([llength $dprops($AppId,AllKElemId)])} {
	    set kwxpath "Applications/$AppId"
	    # For all defined kratos elements
	    foreach celemid $dprops($AppId,AllKElemId) {
		# Check for all defined group identifier for this element
		if {[info exists dprops($AppId,KElem,$celemid,AllGroupId)] && [llength $dprops($AppId,KElem,$celemid,AllGroupId)]} {
		    # For all defined group identifier for this element
		    foreach cgroupid $dprops($AppId,KElem,$celemid,AllGroupId) {
		        # Get the GiD entity type, element type and property identifier
		        lassign $dprops($AppId,KElem,$celemid,$cgroupid,GProps) GiDEntity GiDElemType PropertyId KEKWord nDim
		        if {[GiD_EntitiesGroups get $cgroupid elements -count -element_type $GiDElemType]} {
		            set etbf ""
		            set usennode [::xmlutils::getKKWord $kwxpath $celemid usennode]
		            if {$usennode eq "Yes"} {
		                set etbf [::wkcf::GetnDimnNode $GiDElemType $nDim]
		            }
		            set kelemtype [string trim ${KEKWord}${etbf}]

		            if {[dict get [::wkcf::GetFluidMaterialProperties $AppId] "NonNewtonianFluid"] eq "Yes"} {
		                set kelemtype "HerschelBulkleyVMS3D"
		            }

		            #set GlobalPId $dprops($AppId,KElem,$celemid,$cgroupid,GlobalPId)
		            set cproperty "dv"
		            set container [GetAnalysisDataContainer]
		            set cxpath "$AppId//c.$container//i.SolverType"
		            set MonolithicOption [::xmlutils::setXml $cxpath $cproperty]
		            if {$MonolithicOption eq "Monolithic"} {
		                if {[GetAnalysisDataContainer] eq "AnalysisData"} {
		                    GiD_File fprintf $filechannel "Begin Elements MonolithicDEMCoupled3D   // GUI group identifier: $cgroupid"
		                } else {
		                    GiD_File fprintf $filechannel "Begin Elements $kelemtype   // GUI group identifier: $cgroupid"
		                }
		            } else {
		                GiD_File fprintf $filechannel "Begin Elements $kelemtype   // GUI group identifier: $cgroupid"
		            }
		            foreach elem_id [GiD_EntitiesGroups get $cgroupid elements -element_type $GiDElemType] {
		                #GiD_File fprintf $filechannel "$elem_id $GlobalPId [lrange [GiD_Mesh get element $elem_id] 3 end]"
		                GiD_File fprintf $filechannel "$elem_id 0 [lrange [GiD_Mesh get element $elem_id] 3 end]"
		            }
		            GiD_File fprintf $filechannel "End Elements"
		            GiD_File fprintf $filechannel ""
		        }
		    }
		}
	    }
	    GiD_File fprintf $filechannel ""
	}
    }

    # Special case of DEM application
    if {$AppId eq "DEM"} {
	# DEM Channels
	variable deminletchannel
	variable demfemchannel
	variable dem_ref_to_props_number

	set basexpath "$AppId//c.DEM-Elements"
	#set gelemlist [::xmlutils::setXmlContainerIds $basexpath]

	set gelemlist [::wkcf::GetDEMActiveElements]

	foreach celemid $gelemlist {

	    set elembasexpath "$basexpath//c.DEM-Element"
	    # Get the group node list
	    set cgroupidlist [::xmlutils::setXmlContainerIds $elembasexpath]
	    foreach cgroupid $cgroupidlist {

		set cproperty "dv"
		set cxpath "${elembasexpath}//c.[list ${cgroupid}]//c.Properties//i.SetActive"
		set active_or_not [::xmlutils::setXml $cxpath $cproperty]
		if {$active_or_not=="No"} {
		    continue
		}
		# set props_number 0
		if {$cgroupid != ""} {
		    set elist [GiD_EntitiesGroups get $cgroupid elements]
		    if {[llength $elist]} {
		        # Write all nodes for this group in increasing orden
		        set realelemname [::xmlutils::getKKWord "Applications/$AppId" $celemid]
		        if {[::xmlutils::setXml "$AppId//c.DEM-Elements//c.DEM-Element//c.$cgroupid//c.Properties//i.Material" dv] eq "DEM-IceMaterial"} {
		            set realelemname IceContinuumParticle3D
		        }
		        GiD_File fprintf $filechannel "Begin Elements $realelemname   //  GUI group identifier: $cgroupid"
		        foreach elem_id [GiD_EntitiesGroups get $cgroupid elements] {
		            GiD_File fprintf $filechannel "$elem_id $dem_ref_to_props_number [lrange [GiD_Mesh get element $elem_id] 3 3]"
		        }
		        GiD_File fprintf $filechannel "End Elements"
		        GiD_File fprintf $filechannel ""
		    }
		}
		incr dem_ref_to_props_number
	    }
	}
	set rootid $AppId
	set basexpath "$rootid//c.DEM-Conditions//c.DEM-FEM-Wall"
	set gproplist [::xmlutils::setXmlContainerIds $basexpath]

  variable demwall_element_props
	foreach cgroupid $gproplist {

    set cproperty "dv"
    set cxpath "$rootid//c.DEM-Conditions//c.DEM-FEM-Wall//c.[list ${cgroupid}]//i.SetActive"
    set active_or_not [::xmlutils::setXml $cxpath $cproperty]
    if {$active_or_not=="No"} {
	  continue
    }

    # Get the group node list
    set nlist [GiD_EntitiesGroups get $cgroupid elements]
    if {[llength $nlist]} {
		set elemname [::wkcf::GetDEMFEMElementName [lindex $nlist 0]]
		#
		set cxpath "$rootid//c.DEM-Conditions//c.DEM-FEM-Wall//c.[list ${cgroupid}]//c.Options//i.AnalyticProps"
	    set compute_analytic_data [::xmlutils::setXml $cxpath dv]
	    if {$compute_analytic_data == "Yes"} {
	    set elemname AnalyticRigidFace3D3N
	    }
		#
		GiD_File fprintf $demfemchannel "Begin Conditions $elemname // GUI DEM-FEM-Wall group identifier: $cgroupid"
		foreach elemid $nlist {
		    GiD_File fprintf $demfemchannel "$elemid $demwall_element_props($cgroupid) [lrange [GiD_Mesh get element $elemid] 3 end]"
		}
		GiD_File fprintf $demfemchannel "End Conditions"
		GiD_File fprintf $demfemchannel ""
		set options_container Options
		if {"Fluid" in $ActiveAppList} {
		    set options_container [GetOptionsContainer]
		}
		set cxpath "${basexpath}//c.[list ${cgroupid}]//c.$options_container"
		set allgprop [::xmlutils::setXmlContainerPairs $cxpath "" "dv"]
		set findeb [lsearch -index 0 $allgprop "is-embedded-in-fluid"]
		set IsEmbedded "No"
		if {$findeb !="-1"} {
		    set IsEmbedded [lindex $allgprop $findeb 1]
		}
		if {$IsEmbedded eq "Yes"} {
		    GiD_File fprintf $demfemchannel "Begin Elements RigidShellElement // GUI Inlet group identifier: $cgroupid"
		    foreach elemid $nlist {
		        GiD_File fprintf $demfemchannel "$elemid  $number_of_demfem_wall  [lrange [GiD_Mesh get element $elemid] 3 end]"
		    }
		    GiD_File fprintf $demfemchannel "End Elements"
		    GiD_File fprintf $demfemchannel ""
		}
	    }
	    incr number_of_demfem_wall
	}
    }

    # For debug
    if {!$::wkcf::pflag} {
	set endtime [clock seconds]
	set ttime [expr $endtime-$inittime]
	# WarnWinText "endtime:$endtime ttime:$ttime"
	WarnWinText "Write element connectivities: [::KUtils::Duration $ttime]"
    }
}

proc ::wkcf::WriteBoundaryConditions {AppId} {
    # Write the boundary condition block
    variable dprops
    # Check for all defined condition type
    if {([info exists dprops($AppId,AllBCTypeId)]) && ([llength $dprops($AppId,AllBCTypeId)])} {
	# For debug
	if {!$::wkcf::pflag} {
	    set inittime [clock seconds]
	}
	set inletvelglist [list]; set noslipglist [list]
	set flagvariablelist [list]
	# For all defined condition identifier
	foreach ccondid $dprops($AppId,AllBCTypeId) {
	    # WarnWinText "ccondid:$ccondid"
	    # Check for all defined group identifier inside this condition type
	    if {([info exists dprops($AppId,BC,$ccondid,AllGroupId)]) && ([llength $dprops($AppId,BC,$ccondid,AllGroupId)])} {
		# Select the condition type
		switch -exact -- $ccondid {
		    "Displacements" {
		        set kwxpath "Applications/StructuralAnalysis"
		        set kwordlist [list [::xmlutils::getKKWord $kwxpath "Dx" "kkword"] [::xmlutils::getKKWord $kwxpath "Dy" "kkword"] [::xmlutils::getKKWord $kwxpath "Dz" "kkword"]]
		        # Process displacement properties
		        ::wkcf::WriteDispRotBC $AppId $ccondid $kwordlist
		    }
		    "Rotations" {
		        set usenbst2 "No"
		        set useshells2 "No"
		        set usebeams "No"
		        set beamlist [list "BeamElement"]
		        set shelllist2 [list "ShellIsotropic" "ShellAnisotropic" "EBST"]
		        if {([info exists dprops($AppId,AllKElemId)]) && ($dprops($AppId,AllKElemId)>0)} {
		            # For all defined kratos elements
		            foreach celemid $dprops($AppId,AllKElemId) {
		                if {$celemid in $shelllist2} {
		                    set useshells2 "Yes"
		                } elseif {$celemid in $beamlist} {
		                    set usebeams "Yes"
		                }
		            }
		        }
		        if {($useshells2 eq "Yes")||($usebeams eq "Yes")} {
		            set kwxpath "Applications/StructuralAnalysis"
		            set kwordlist [list [::xmlutils::getKKWord $kwxpath "Rx"] [::xmlutils::getKKWord $kwxpath "Ry"] [::xmlutils::getKKWord $kwxpath "Rz"]]
		            # Process displacement properties
		            ::wkcf::WriteDispRotBC $AppId $ccondid $kwordlist
		        }
		    }
		    "OutletPressure" {
		        if {$AppId=="Fluid"} {
		            set cproperty "dv"
		            set cxpath "$AppId//c.AnalysisData//i.FluidApproach"
		            set FluidApproach [::xmlutils::setXml $cxpath $cproperty]
		            if { $FluidApproach eq "Eulerian" } {
		                # Solver type
		                set container [GetAnalysisDataContainer]
		                set cxpath "$AppId//c.$container//i.SolverType"
		                set SolverType [::xmlutils::setXml $cxpath $cproperty]
		                # WarnWinText "SolverType:$SolverType"
		                set kwxpath "Applications/$AppId"

		                switch -exact -- $SolverType {
		                    "ElementBased" {
		                        set kwid "OutletPressureFractionalStep"
		                        set kwordlist [list "[::xmlutils::getKKWord $kwxpath $kwid]"]
		                    }
		                    "Monolithic" {
		                        set kwid "OutletPressureMonolithic"
		                        set kwordlist [list "[::xmlutils::getKKWord $kwxpath $kwid]"]
		                    }
		                }

		                # Process outlet pressure
		                ::wkcf::WriteOutLetPressureBC $AppId $ccondid $kwordlist
		            }
		        }
		    }
		    "InletVelocity" {
		        if {[llength $dprops($AppId,BC,$ccondid,AllGroupId)]} {
		            foreach _gid $dprops($AppId,BC,$ccondid,AllGroupId) {
		                lappend inletvelglist $_gid
		            }
		        }
		    }
		    "No-Slip" {
		        if {[llength $dprops($AppId,BC,$ccondid,AllGroupId)]} {
		            foreach _gid $dprops($AppId,BC,$ccondid,AllGroupId) {
		                lappend noslipglist $_gid
		            }
		        }
		    }
		    "Flag-Variable" {
		        if {[llength $dprops($AppId,BC,$ccondid,AllGroupId)]} {
		            foreach _gid $dprops($AppId,BC,$ccondid,AllGroupId) {
		                lappend flagvariablelist $_gid
		            }
		        }
		    }
		    "Is-Slip" - "WallLaw" {
		        # Write is-slip/walllaw boundary condition
		        set kwordlist [list "IS_STRUCTURE" "Y_WALL"]
		        ::wkcf::WriteFluidIsSlipWallLawBC $AppId $ccondid $kwordlist
		    }
		    "Distance" {
		        # Write distance boundary condition
		        set kwordlist [list "DISTANCE"]
		        ::wkcf::WriteFluidDistanceBC $AppId $ccondid $kwordlist
		    }
		    "PFEMWall--" {
		        # Commented by J. Garate on 17/12/2012
		        set kwordlist [list "LINEAR_VELOCITY_X" "LINEAR_VELOCITY_Y" "LINEAR_VELOCITY_Z" "ANGULAR_VELOCITY_X" "ANGULAR_VELOCITY_Y" "ANGULAR_VELOCITY_Z" ]
		        ::wkcf::WriteFluidPFEMWallBC $AppId $ccondid $kwordlist
		    }
		    "PFEMFluidInlet" {
		        # Write PFEM Fluid Velocity
		        set kwordlist [list "VELOCITY_X" "VELOCITY_Y" "VELOCITY_Z"]
		        ::wkcf::WriteFluidPFEMInletBC $AppId $ccondid $kwordlist
		    }
		    "PrescribedTemperature" {
		        # Write prescribed temperature condition
		        set kwordlist [list "TEMPERATURE"]
		        ::wkcf::WriteConvectionDiffusionPrescribedTemperatureBC $AppId $ccondid $kwordlist
		    }
		    "HeatFlux" {
		        # Write prescribed heat flux condition
		        set kwordlist [list "HEAT_FLUX"]
		        ::wkcf::WriteConvectionDiffusionPrescribedHeatFluxBC $AppId $ccondid $kwordlist
		    }
		    "FaceHeatFlux" {
		        # Write prescribed face heat flux condition
		        set kwordlist [list "FACE_HEAT_FLUX"]
		        ::wkcf::WriteConvectionDiffusionPrescribedFaceHeatFluxBC $AppId $ccondid $kwordlist
		    }
		    "PFEMFixedWall" {
		        # Write PFEM fixed wall boundary condition
		        set kwordlist [list "IS_STRUCTURE" "DISPLACEMENT_X" "DISPLACEMENT_Y" "DISPLACEMENT_Z"]
		        ::wkcf::WritePFEMLagrangianFluidFixedWallBC $AppId $ccondid $kwordlist
		    }
		}
	    }
	}

	# For fluid application
	if {$AppId=="Fluid"} {
	    set kwxpath "Applications/$AppId"
	    set kwordlist [list [::xmlutils::getKKWord $kwxpath "Vx"] [::xmlutils::getKKWord $kwxpath "Vy"] [::xmlutils::getKKWord $kwxpath "Vz"]]
	    ::wkcf::WriteFluidBC $AppId $inletvelglist $noslipglist $flagvariablelist $kwordlist
	}

	# For debug
	if {!$::wkcf::pflag} {
	    set endtime [clock seconds]
	    set ttime [expr $endtime-$inittime]
	    # WarnWinText "endtime:$endtime ttime:$ttime"
	    WarnWinText "Write boundary conditions: [::KUtils::Duration $ttime]"
	}
    }
}

proc ::wkcf::WriteInitialConditions {AppId} {
    # Write the initial condition block
    variable dprops

    # Check for all defined initial condition type
    if {([info exists dprops($AppId,AllICTypeId)]) && ([llength $dprops($AppId,AllICTypeId)])} {
	# For debug
	if {!$::wkcf::pflag} {
	    set inittime [clock seconds]
	}
	set inletvelglist [list]; set noslipglist [list]
	set flagvariablelist [list]
	# For all defined initial condition identifier
	foreach ccondid $dprops($AppId,AllICTypeId) {
	    # WarnWinText "ccondid:$ccondid"
	    # Check for all defined group identifier inside this initial condition type
	    if {([info exists dprops($AppId,IC,$ccondid,AllGroupId)]) && ([llength $dprops($AppId,IC,$ccondid,AllGroupId)])} {
		# Select the initial condition type
		switch -exact -- $ccondid {
		    "InitialTemperature" {
		        set kwxpath "Applications/$AppId"
		        set kwordlist [list [::xmlutils::getKKWord $kwxpath "$ccondid" "kkword"]]
		        # Process initial temperature properties
		        ::wkcf::WriteInitialTemperatureIC $AppId $ccondid $kwordlist
		    }
		}
	    }
	}

	# For debug
	if {!$::wkcf::pflag} {
	    set endtime [clock seconds]
	    set ttime [expr $endtime-$inittime]
	    # WarnWinText "endtime:$endtime ttime:$ttime"
	    WarnWinText "Write initial conditions: [::KUtils::Duration $ttime]"
	}
    }
}

proc ::wkcf::WriteConditions {AppId} {
    # ABSTRACT: Write condition properties
    variable dprops; variable ndime
    variable ctbclink

    variable filechannel

    # Check for all defined kratos elements
    if {([info exists dprops($AppId,AllKElemId)]) && ([llength $dprops($AppId,AllKElemId)])} {
	# For debug
	if {!$::wkcf::pflag} {
	    set inittime [clock seconds]
	}

	set fixval "0"
	# Write conditions
	# Select the condition identifier
	set ConditionId "Condition"
	set cproperty "dv"

	# Select the application type
	# Fluid
	if {$AppId == "Fluid"} {
	    set cxpath "$AppId//c.AnalysisData//i.FluidApproach"
	    set FluidApproach [::xmlutils::setXml $cxpath $cproperty]
	    if { $FluidApproach eq "Eulerian" } {
		# Solver type
		set container [GetAnalysisDataContainer]
		set cxpath "$AppId//c.$container//i.SolverType"
		set SolverType [::xmlutils::setXml $cxpath $cproperty]
		# WarnWinText "SolverType:$SolverType"
		switch -exact -- $SolverType {
		    "ElementBased" {
		        set ConditionId "WallConditionDiscontinuous${ndime}"
		    }
		    "Monolithic" {
		        set ConditionId "MonolithicWallCondition${ndime}"
		    }
		}
	    } elseif {$FluidApproach eq "PFEM-Lagrangian"} {
		set ConditionId "Condition${ndime}"
	    }
	} elseif {$AppId == "ConvectionDiffusion"} {

	    # Convection diffusion
	    set ConditionId "ThermalFace${ndime}"
	}

	if {$ndime =="2D"} {
	    # 2D
	    set cgroupid "-AKGSkinMesh2D"
	    set GiDElemType "Linear"
	    if {[GiD_Groups exists $cgroupid]} {
		if {[GiD_EntitiesGroups get $cgroupid elements -count -element_type $GiDElemType]} {
		    # Write all conditions for 2D case
		    GiD_File fprintf $filechannel "%s" "Begin Conditions $ConditionId"
		    set condid 0

		    foreach elem_id [GiD_EntitiesGroups get $cgroupid elements -element_type $GiDElemType] {
		        incr condid 1
		        set nodes [lrange [GiD_Mesh get element $elem_id] 3 end]
		        # nodes = [list nodei nodej]
		        GiD_File fprintf $filechannel "%10d %4d %10d %10d" $condid $fixval [lindex $nodes 0] [lindex $nodes 1]
		        # Update the link between the condition id. and the BC element id
		        dict set ctbclink $elem_id $condid
		    }

		    GiD_File fprintf $filechannel "%s" "End Conditions"
		    GiD_File fprintf $filechannel ""
		}
	    }
	} elseif {$ndime =="3D"} {
	    # 3D
	    set cgroupid "-AKGSkinMesh3D"
	    set GiDElemType "Triangle"
	    set usetriangle 0
	    if {[GiD_Groups exists $cgroupid]} {

		if {[GiD_EntitiesGroups get $cgroupid elements -count -element_type $GiDElemType]} {
		    # Write all conditions for 3D case
		    GiD_File fprintf $filechannel "%s" "Begin Conditions $ConditionId"
		    set condid 0

		    foreach elem_id [GiD_EntitiesGroups get $cgroupid elements -element_type $GiDElemType] {
		        incr condid 1
		        set nodes [lrange [GiD_Mesh get element $elem_id] 3 end]
		        set ni [lindex $nodes 0]
		        set nj [lindex $nodes 1]
		        set nk [lindex $nodes 2]
		        # msg "$ni $nj $nk $condid"
		        GiD_File fprintf $filechannel "%10d %4d %10d %10d %10d" $condid $fixval $ni $nj $nk
		        # Update the link between the condition id. and the BC element id
		        dict set ctbclink $elem_id $condid
		    }
		    GiD_File fprintf $filechannel "%s" "End Conditions"
		    GiD_File fprintf $filechannel ""
		    set usetriangle 1
		}

		set GiDElemType "Quadrilateral"

		if {[GiD_EntitiesGroups get $cgroupid elements -count -element_type $GiDElemType]} {
		    # Write all conditions for 3D case
		    GiD_File fprintf $filechannel "%s" "Begin Conditions $ConditionId"
		    if {!$usetriangle} {
		        set condid 0
		    }
		    foreach elem_id [GiD_EntitiesGroups get $cgroupid elements -element_type $GiDElemType] {
		        incr condid 1
		        set nodes [lrange [GiD_Mesh get element $elem_id] 3 end]
		        set ni [lindex $nodes 0]
		        set nj [lindex $nodes 1]
		        set nk [lindex $nodes 2]
		        set nl [lindex $nodes 3]
		        # msg "$ni $nj $nk $nl $condid"
		        GiD_File fprintf $filechannel "%10d %4d %10d %10d %10d %10d" $condid $fixval $ni $nj $nk $nl
		        # Update the link between the condition id. and the BC element id
		        dict set ctbclink $elem_id $condid
		    }
		    GiD_File fprintf $filechannel "%s" "End Conditions"
		    GiD_File fprintf $filechannel ""
		}
	    }
	}

	# For debug
	if {!$::wkcf::pflag} {
	    set endtime [clock seconds]
	    set ttime [expr $endtime-$inittime]
	    # WarnWinText "endtime:$endtime ttime:$ttime"
	    WarnWinText "Write conditions: [::KUtils::Duration $ttime]"
	}
    }
}

proc ::wkcf::WriteDSOLIDMaterialsFile {} {
    variable DSOLID
    set ppfilename "materials.py"
    set PDir [::KUtils::GetPaths "PDir"]

    set fullname [file native [file join $PDir $ppfilename]]
    # First delete the file
    if {[file exists $fullname]} {
	set res [file delete -force $fullname]
    }
    if { [catch { set fileid [open $fullname w+] }] } {
	WarnWin [= "Cannot write file %s. Permission denied" $fullname].
	return 0
    }

    # Write DSOLID materials parameter file
    if {$DSOLID =="Yes"} {
	::wkcf::WriteDSOLIDMaterials "DSOLID" $fileid $PDir
    }
    close $fileid
}
proc ::wkcf::WriteProjectParameters {} {
    # Write the project parameters file
    variable StructuralAnalysis;     variable FluidApplication
    variable ConvectionDiffusionApplication
    variable DSOLID

    set ppfilename "ProjectParameters.py"
    set PDir [::KUtils::GetPaths "PDir"]

    set fullname [file native [file join $PDir $ppfilename]]
    # First delete the file
    if {[file exists $fullname]} {
	set res [file delete -force $fullname]
    }

    if { [catch { set fileid [open $fullname w+] }] } {

	WarnWin [= "Cannot write file %s. Permission denied" $fullname].
	return 0
    }


    if {$StructuralAnalysis =="Yes"} {
	::wkcf::WriteStructuralProjectParameters "StructuralAnalysis" $fileid $PDir
    }

    # Write solid application project parameter file
    if {$DSOLID =="Yes"} {
	::wkcf::WriteSolidProjectParameters "StructuralAnalysis" $fileid $PDir
    }
    # Write fluid application project parameter file
    if {$FluidApplication =="Yes"} {
	::wkcf::WriteFluidProjectParameters "Fluid" $fileid $PDir
    }

    # Write convection diffusion application project parameter file
    if {$ConvectionDiffusionApplication =="Yes"} {
	::wkcf::WriteConvectionDiffusionProjectParameters "ConvectionDiffusion" $fileid $PDir
    }

    close $fileid
}

proc ::wkcf::WriteGiDPostMode {AppId fileid} {
    # Write the GiD post mode variables for each applications
    # kratos key word xpath
    set kwxpath "Applications/$AppId"
    # Gid results
    set gidrlist [list "GiDPostMode" "GiDWriteMeshFlag" "GiDWriteConditionsFlag" "GiDWriteParticlesFlag" "GiDMultiFileFlag"]
    foreach gidr $gidrlist {
	# Get the value
	set cxpath "$AppId//c.Results//c.GiDOptions//i.[list ${gidr}]"
	set cproperty "dv"
	set cvalue [::xmlutils::setXml $cxpath $cproperty]
	if {$cvalue != ""} {
	    # Get the kratos keyword
	    set gidrkw [::xmlutils::getKKWord $kwxpath $gidr]
	    if {($gidr == "GiDWriteMeshFlag") || ($gidr == "GiDWriteConditionsFlag") || ($gidr == "GiDWriteParticlesFlag")} {
		if {$cvalue == "Yes"} {
		    set cvalue True
		} else {
		    set cvalue False
		}
		puts $fileid "$gidrkw = $cvalue"
	    } else {
		puts $fileid "$gidrkw = \"$cvalue\""
	    }
	}
    }
}

proc ::wkcf::CloseChannels {AppId} {
    variable filechannel
    variable deminletchannel
    variable demfemchannel
    # Close the file
    GiD_File fclose $filechannel
    if {$AppId eq "DEM"} {
	GiD_File fclose $deminletchannel
	GiD_File fclose $demfemchannel
    }
}

proc ::wkcf::GetElementCenter {element_id} {
    set element_data [GiD_Mesh get element $element_id]
    set num_nodes [lindex $element_data 2]
    set node_ids [lrange $element_data 3 2+$num_nodes]
    set sum {0 0 0}
    foreach node_id $node_ids {
	set coordinates [lrange [GiD_Mesh get node $node_id] 1 end]
	set sum [MathUtils::VectorSum $coordinates $sum]
    }
    return [MathUtils::ScalarByVectorProd [expr {1.0/$num_nodes}] $sum]
}

proc ::wkcf::GetNodeHigherentities {node_id} {
    set node_data [GiD_Info list_entities nodes $node_id]
    if {![regexp {HigherEntity: ([0-9]+)} $node_data dummy higherentity]} {
	set higherentity 9999; #the node does not exist, return > 0 to not delete it
    }
    return $higherentity
}

proc ::wkcf::Compute_External_Elements {ndime cgroupid element_ids} {

    set mesh_elements [GiD_EntitiesGroups get $cgroupid all_mesh]
    set real_mesh_elements [lindex $mesh_elements 1]
    set list_of_faces [list]
    foreach mesh_element_id $real_mesh_elements {
	set line($mesh_element_id) [GiD_Mesh get element $mesh_element_id]
	set partial_list_of_faces [lrange [GiD_Mesh get element $mesh_element_id] 3 end]
	lappend list_of_faces {*}$partial_list_of_faces
    }
    set unrepeated_list [lsort -integer -unique $list_of_faces]
    if {$ndime == "3D"} {
	set elements_in_common 6 ; #TODO: Check this constant
    } else {
	set elements_in_common 3 ; #TODO: Check this constant
    }

    foreach list_elem $unrepeated_list {
	set result($list_elem) [lsearch -all $list_of_faces $list_elem]
	set length($list_elem) [llength $result($list_elem)]
	if {$length($list_elem)>$elements_in_common} {
	    set todelete($list_elem) 1
	} else {
	    set todelete($list_elem) 0
	}
    }
    foreach list_elem $unrepeated_list {
	if {$todelete($list_elem)==1} {
	    set list_of_faces [lsearch -all -inline -not -exact $list_of_faces $list_elem]
	}
    }
    set unrepeated_list_exterior_nodes [lsort -integer -unique $list_of_faces]

    foreach element_id $element_ids { ; # Here we loop on each of the elements by id
	if {$ndime == "3D"} {
	    set element_nodes [lrange [GiD_Mesh get element $element_id] 3 end] ; # We get the nodes of the element
	} else {
	    set element_nodes [lrange [GiD_Mesh get element $element_id] 3 end]
	}
	set is_external_element($element_id) 0
	foreach element_node $element_nodes {
	    if {[lsearch $unrepeated_list_exterior_nodes $element_node] != -1} {
		set is_external_element($element_id) 1
		break
	    }
	}
    }
    return [array get is_external_element]
}

proc ::wkcf::Elements_Substitution {} {
    global KPriv
    variable ndime
    set fail 0
    package require math::statistics
    set seed [expr srand(0)]
    set gelemlist [::wkcf::GetDEMActiveElements]
    set final_list_of_isolated_nodes [list]
    foreach celemid $gelemlist {
	set cgroupidlist [::xmlutils::setXmlContainerIds {DEM//c.DEM-Elements//c.DEM-Element}]
	set cohesive_groups_list [::xmlutils::setXmlContainerIds {DEM//c.DEM-Cohesivegroup}]

	if {![GiD_Groups exists SKIN_SPHERE_DO_NOT_DELETE]} {
	    GiD_Groups create SKIN_SPHERE_DO_NOT_DELETE
	}

	foreach cgroupid $cgroupidlist { ; #Loop on each of the groups, be it a point, a line, a triangle, etc
	    set list_of_elements_to_add_to_skin_sphere_group [list]
	    array set my_prop [join [::xmlutils::setXmlContainerPairs "DEM//c.DEM-Elements//c.DEM-Element//c.$cgroupid//c.Properties" {} dv]]
	    set use_advanced_meshing_features $my_prop(AdvancedMeshingFeatures)
	    if {([lsearch $cohesive_groups_list $cgroupid] == -1) && ($KPriv(what_dempack_package) eq "C-DEMPack") && ($use_advanced_meshing_features eq "No")} { ; #Non-cohesive, already meshed with cuban mesher
		set cgroupid_elements [GiD_EntitiesGroups get $cgroupid elements]

		foreach cgroupid_element $cgroupid_elements {
		    if {[lsearch [GiD_EntitiesGroups get SKIN_SPHERE_DO_NOT_DELETE elements] $cgroupid_element] == -1} {
		        lappend list_of_elements_to_add_to_skin_sphere_group $cgroupid_element
		    }
		}
		GiD_EntitiesGroups assign SKIN_SPHERE_DO_NOT_DELETE elements $list_of_elements_to_add_to_skin_sphere_group
		continue
	    }

	    if {$use_advanced_meshing_features eq "Yes"} {
		if {$my_prop(AdvancedMeshingFeaturesAlgorithmType) eq "FEMtoDEM"} {
		    set use_fem_to_dem $my_prop(FEMtoDEM)
		    set element_radius [expr {0.5*$my_prop(Diameter)}]
		    set standard_deviation $my_prop(StandardDeviation)
		    set probldistr $my_prop(ProbabilityDistribution)
		    set min_radius [expr {0.5*$element_radius}] ; # HARDCODED THE CRITERIA TO CHOOSE BOTH MINIMUM AND MAXIMUM RADIUS
		    set max_radius [expr {1.5*$element_radius}]
		    if {$use_fem_to_dem == "AttheCentroid"} {
		        set nodes_to_delete [list]
		        set element_ids [GiD_EntitiesGroups get $cgroupid elements] ; # We get the ids of all the elements in cgroupid
		        array set is_external_element [::wkcf::Compute_External_Elements $ndime $cgroupid $element_ids]

		        foreach element_id $element_ids { ; # Here we loop on each of the elements by id
		            if {$ndime == "3D"} {
		                set element_nodes [lrange [GiD_Mesh get element $element_id] 3 end] ; # We get the nodes of the element
		            } else {
		                set element_nodes [lrange [GiD_Mesh get element $element_id] 3 end]
		            }
		            lappend nodes_to_delete {*}$element_nodes ; # We add those nodes to the nodes_to_delete list
		            if {$probldistr == "NormalDistribution"} {
		                set final_elem_radius [::wkcf::NormalDistribution $element_radius $standard_deviation $min_radius $max_radius]
		            } else {
		                set final_elem_radius [::wkcf::LognormalDistribution $element_radius $standard_deviation $min_radius $max_radius]
		            }
		            set node_id [GiD_Mesh create node append [wkcf::GetElementCenter $element_id]] ; # We create a new node starting from the center of the given element

		            if {$ndime == "3D"} {
		                set new_element_id [GiD_Mesh create element append sphere 1 $node_id $final_elem_radius] ; # We create a new sphere element starting from the previous node and obtain its id
		                lappend list_of_elements_to_add_to_skin_sphere_group {*}$new_element_id
		                if {($is_external_element($element_id)==1) && ([lsearch $cohesive_groups_list $cgroupid] != -1)} {
		                    GiD_EntitiesGroups assign SKIN_SPHERE_DO_NOT_DELETE elements $new_element_id
		                }
		            } else {
		                set new_element_id [GiD_Mesh create element append circle 1 $node_id $final_elem_radius 0 0 1] ; # We assume the 2D problem is contained in the XY plane
		                lappend list_of_elements_to_add_to_skin_sphere_group {*}$new_element_id
		                if {($is_external_element($element_id)==1) && ([lsearch $cohesive_groups_list $cgroupid] != -1)} {
		                    GiD_EntitiesGroups assign SKIN_SPHERE_DO_NOT_DELETE elements $new_element_id
		                }
		            }
		            foreach container_group [GiD_EntitiesGroups entity_groups elements $element_id] { ; # We get the list of groups to which the element with id $element_id belongs
		                GiD_EntitiesGroups assign $container_group elements $new_element_id ; # We assign the element with id $new_element_id to each of the groups in the loop
		            }
		        }

		        if {[lsearch $cohesive_groups_list $cgroupid] == -1} {
		            GiD_EntitiesGroups assign SKIN_SPHERE_DO_NOT_DELETE elements $list_of_elements_to_add_to_skin_sphere_group
		        }

		        GiD_Mesh delete element [GiD_EntitiesGroups get $cgroupid elements -element_type hexahedra]
		        GiD_Mesh delete element [GiD_EntitiesGroups get $cgroupid elements -element_type tetrahedra]
		        set nodes_to_delete [lsort -integer -unique $nodes_to_delete] ; # We reorder the list and remove repeated nodes
		        foreach node_id $nodes_to_delete {
		            set gid_info [GiD_Info list_entities nodes $node_id]
		            if {![wkcf::GetNodeHigherentities $node_id]} { ; # if this node does not have higher entities
		                GiD_Mesh delete node $node_id ; #delete the nodes of the element as long as it does not have higher entities
		            }
		        }
		        set point_node_ids [GiD_EntitiesGroups get $cgroupid nodes] ; # This list exists only for groups made up of isolated points
		        foreach node_id $point_node_ids {
		            if {![wkcf::GetNodeHigherentities $node_id]} {
		                if {$probldistr == "NormalDistribution"} {
		                    set final_elem_radius [::wkcf::NormalDistribution $element_radius $standard_deviation $min_radius $max_radius]
		                } else {
		                    set final_elem_radius [::wkcf::LognormalDistribution $element_radius $standard_deviation $min_radius $max_radius]
		                }
		                if {$ndime == "3D"} {
		                    set new_element_id [GiD_Mesh create element append sphere 1 $node_id $final_elem_radius] ; # We create a new sphere element starting from the previous node and obtain its id
		                } else {
		                    set new_element_id [GiD_Mesh create element append circle 1 $node_id $final_elem_radius 0 0 1] ; # We assume the 2D problem is contained in the XY plane
		                }
		                set list_of_groups_containing_this_elem [GiD_EntitiesGroups entity_groups nodes $node_id]
		                foreach container_group $list_of_groups_containing_this_elem {
		                    GiD_EntitiesGroups assign $container_group elements $new_element_id
		                }
		            }
		        }
		        set extra_nodes [GiD_EntitiesGroups get $cgroupid nodes]
		        foreach node_id $extra_nodes {
		            if {![wkcf::GetNodeHigherentities $node_id]} {
		                GiD_Mesh delete node $node_id
		            }
		        }
		    } elseif {$use_fem_to_dem == "AttheNodes"} {
		        # We first delete the elements (lines, triangles, quadrilaterals, tetraedra or hexahedra) of this group,
		        # but not their nodes, which will be used for creating the new sheres
		        GiD_Mesh delete element [GiD_EntitiesGroups get $cgroupid elements]
		        foreach node_id [GiD_EntitiesGroups get $cgroupid nodes] {
		            if {$probldistr == "NormalDistribution"} {
		                set final_elem_radius [::wkcf::NormalDistribution $element_radius $standard_deviation $min_radius $max_radius]
		            } else {
		                set final_elem_radius [::wkcf::LognormalDistribution $element_radius $standard_deviation $min_radius $max_radius]
		            }
		            if {$ndime == "3D"} {
		                set new_element_id [GiD_Mesh create element append sphere 1 $node_id $final_elem_radius] ; # We create a new sphere element starting from the previous node and obtain its id
		                lappend list_of_elements_to_add_to_skin_sphere_group {*}$new_element_id
		            } else {
		                set new_element_id [GiD_Mesh create element append circle 1 $node_id $final_elem_radius 0 0 1] ; # We assume the 2D problem is contained in the XY plane
		                lappend list_of_elements_to_add_to_skin_sphere_group {*}$new_element_id
		            }

		            set list_of_groups_containing_this_elem [GiD_EntitiesGroups entity_groups nodes $node_id]
		            foreach container_group $list_of_groups_containing_this_elem {
		                GiD_EntitiesGroups assign $container_group elements $new_element_id
		            }
		        }

		        if {[lsearch $cohesive_groups_list $cgroupid] == -1} {
		            GiD_EntitiesGroups assign SKIN_SPHERE_DO_NOT_DELETE elements $list_of_elements_to_add_to_skin_sphere_group
		        }

		    } elseif {$use_fem_to_dem == "AtBothNodesAndCentroids"} {
		        set nodes_to_delete [list]
		        set element_ids [GiD_EntitiesGroups get $cgroupid elements] ; # We get the ids of all the elements in cgroupid

		        foreach element_id $element_ids { ; # Here we loop on each of the elements by id
		            if {$ndime == "3D"} {
		                set element_nodes [lrange [GiD_Mesh get element $element_id] 3 end] ; # We get the nodes of the element
		            } else {
		                set element_nodes [lrange [GiD_Mesh get element $element_id] 3 end]
		            }
		            lappend nodes_to_delete {*}$element_nodes ; # We add those nodes to the nodes_to_delete list
		            if {$probldistr == "NormalDistribution"} {
		                set final_elem_radius [::wkcf::NormalDistribution $element_radius $standard_deviation $min_radius $max_radius]
		            } else {
		                set final_elem_radius [::wkcf::LognormalDistribution $element_radius $standard_deviation $min_radius $max_radius]
		            }
		            set node_id [GiD_Mesh create node append [wkcf::GetElementCenter $element_id]] ; # We create a new node starting from the center of the given element

		            if {$ndime == "3D"} {
		                set new_element_id [GiD_Mesh create element append sphere 1 $node_id $final_elem_radius] ; # We create a new sphere element starting from the previous node and obtain its id
		                lappend list_of_elements_to_add_to_skin_sphere_group {*}$new_element_id
		            } else {
		                set new_element_id [GiD_Mesh create element append circle 1 $node_id $final_elem_radius 0 0 1] ; # We assume the 2D problem is contained in the XY plane
		                lappend list_of_elements_to_add_to_skin_sphere_group {*}$new_element_id
		            }

		            foreach container_group [GiD_EntitiesGroups entity_groups elements $element_id] { ; # We get the list of groups to which the element with id $element_id belongs
		                GiD_EntitiesGroups assign $container_group elements $new_element_id ; # We assign the element with id $new_element_id to each of the groups in the loop
		            }
		        }

		        GiD_Mesh delete element [GiD_EntitiesGroups get $cgroupid elements -element_type hexahedra] ; # TODO done again at the end?
		        GiD_Mesh delete element [GiD_EntitiesGroups get $cgroupid elements -element_type tetrahedra] ; # TODO done again at the end?

		        foreach node_id [GiD_EntitiesGroups get $cgroupid nodes] {
		            if {$probldistr == "NormalDistribution"} {
		                set final_elem_radius [::wkcf::NormalDistribution $element_radius $standard_deviation $min_radius $max_radius]
		            } else {
		                set final_elem_radius [::wkcf::LognormalDistribution $element_radius $standard_deviation $min_radius $max_radius]
		            }
		            if {$ndime == "3D"} {
		                set new_element_id [GiD_Mesh create element append sphere 1 $node_id $final_elem_radius] ; # We create a new sphere element starting from the previous node and obtain its id
		                lappend list_of_elements_to_add_to_skin_sphere_group {*}$new_element_id
		            } else {
		                set new_element_id [GiD_Mesh create element append circle 1 $node_id $final_elem_radius 0 0 1] ; # We assume the 2D problem is contained in the XY plane
		                lappend list_of_elements_to_add_to_skin_sphere_group {*}$new_element_id
		            }

		            set list_of_groups_containing_this_elem [GiD_EntitiesGroups entity_groups nodes $node_id]
		            foreach container_group $list_of_groups_containing_this_elem {
		                GiD_EntitiesGroups assign $container_group elements $new_element_id
		            }
		        }

		        if {[lsearch $cohesive_groups_list $cgroupid] == -1} {
		            GiD_EntitiesGroups assign SKIN_SPHERE_DO_NOT_DELETE elements $list_of_elements_to_add_to_skin_sphere_group
		        }

		    }
		} else {
		    # 2D to 3D algorithm
		    lassign [lindex [GiD_Info Mesh elements Circle -array] 0] type element_ids element_nodes element_materials element_radii
		    foreach element_id $element_ids element_node [lindex $element_nodes 0] element_radius $element_radii {
		        set element_info($element_id) [list $element_node $element_radius]
		    }
		    set element_list [GiD_EntitiesGroups get $cgroupid elements]
		    foreach element_id $element_list {
		        lassign $element_info($element_id) element_node element_radius
		        lappend group_nodes($cgroupid) $element_node
		        lappend group_radius($cgroupid) $element_radius
		    }

		    GiD_Mesh delete element [GiD_EntitiesGroups get $cgroupid elements]
		    foreach node_id $group_nodes($cgroupid) radius $group_radius($cgroupid) {
		        set final_elem_radius $radius
		        set new_element_id [GiD_Mesh create element append sphere 1 $node_id $final_elem_radius] ; # We create a new sphere element starting from the previous node and obtain its id
		        lappend list_of_elements_to_add_to_skin_sphere_group {*}$new_element_id
		        set list_of_groups_containing_this_elem [GiD_EntitiesGroups entity_groups nodes $node_id]
		        foreach container_group $list_of_groups_containing_this_elem {
		            GiD_EntitiesGroups assign $container_group elements $new_element_id
		        }
		    }

		    if {[lsearch $cohesive_groups_list $cgroupid] == -1} {
		        GiD_EntitiesGroups assign SKIN_SPHERE_DO_NOT_DELETE elements $list_of_elements_to_add_to_skin_sphere_group
		    }
		}
	    }
	    lappend final_list_of_isolated_nodes {*}[lindex [GiD_EntitiesGroups get $cgroupid all_mesh] 0]

	    ::wkcf::Delete_Unnecessary_Elements_From_Mesh $cgroupid
	}
    }

    ::wkcf::Cleaning_Up_Skin_And_Removing_Isolated_Nodes $final_list_of_isolated_nodes

    ::wkcf::Destroy_Skin_Sphere_Group $KPriv(what_dempack_package) ; # Getting rid of the SKIN_SPHERE_DO_NOT_DELETE group when in discontinuum or swimming

    return $fail
}

proc ::wkcf::Delete_Unnecessary_Elements_From_Mesh {cgroupid} {

    #GiD_Mesh delete element [GiD_EntitiesGroups get $cgroupid nodes]
    GiD_Mesh delete element [GiD_EntitiesGroups get $cgroupid elements -element_type linear]
    GiD_Mesh delete element [GiD_EntitiesGroups get $cgroupid elements -element_type triangle]
    GiD_Mesh delete element [GiD_EntitiesGroups get $cgroupid elements -element_type quadrilateral]
    GiD_Mesh delete element [GiD_EntitiesGroups get $cgroupid elements -element_type tetrahedra]
    GiD_Mesh delete element [GiD_EntitiesGroups get $cgroupid elements -element_type hexahedra]
}

proc ::wkcf::Cleaning_Up_Skin_And_Removing_Isolated_Nodes {final_list_of_isolated_nodes} {

    GiD_EntitiesGroups unassign SKIN_SPHERE_DO_NOT_DELETE nodes
    #GiD_Mesh delete element [GiD_EntitiesGroups get SKIN_SPHERE_DO_NOT_DELETE elements -element_type triangle]
    GiD_Mesh delete element [GiD_EntitiesGroups get SKIN_SPHERE_DO_NOT_DELETE elements -element_type quadrilateral]
    GiD_EntitiesGroups unassign SKIN_SPHERE_DO_NOT_DELETE elements [GiD_EntitiesGroups get SKIN_SPHERE_DO_NOT_DELETE elements -element_type linear]
    GiD_EntitiesGroups unassign SKIN_SPHERE_DO_NOT_DELETE elements [GiD_EntitiesGroups get SKIN_SPHERE_DO_NOT_DELETE elements -element_type triangle]
    GiD_EntitiesGroups unassign SKIN_SPHERE_DO_NOT_DELETE elements [GiD_EntitiesGroups get SKIN_SPHERE_DO_NOT_DELETE elements -element_type quadrilateral]

    foreach node_id [lsort -integer -unique $final_list_of_isolated_nodes] {
	if {![wkcf::GetNodeHigherentities $node_id]} {
	    GiD_Mesh delete node $node_id
	}
    }
}

proc ::wkcf::NormalDistribution {mean standard_deviation min_rad max_rad} {
    if {$standard_deviation} {
	set max_iterations 1000 ; #set a maximun number of iterations to avoid an infinite loop
	for {set i 0} {$i < $max_iterations} {incr i} {
	    set u1 [::tcl::mathfunc::rand]
	    set u2 [::tcl::mathfunc::rand]
	    #set distribution [expr {$mean + $standard_deviation * sqrt(-2.0 * log($u1)) * cos(6.28318530717958647692 * $u2)}]
	    set distribution [math::statistics::random-normal $mean $standard_deviation 1] ; # We use the math::statistics library instead
	    if {$distribution > $min_rad && $distribution < $max_rad} {
		return $distribution
	    }
	}
	error "wkcf::NormalDistribution failed after $max_iterations iterations. mean=$mean std_dev=$standard_deviation min_rad=$min_rad max_rad=$max_rad"
    }
    return $mean
}

proc ::wkcf::LognormalDistribution {mean standard_deviation min_rad max_rad} {
    if {$standard_deviation} {
	set log_min [expr log($min_rad)]
	set log_max [expr log($max_rad)]
	set NormalMean [expr {log($mean * $mean / sqrt($standard_deviation * $standard_deviation + $mean * $mean))}]
	set NormalStdDev [expr sqrt(log(1.0 + $standard_deviation * $standard_deviation / ($mean * $mean)))]
	return [expr exp([::wkcf::NormalDistribution $NormalMean $NormalStdDev $log_min $log_max])]
    }
    return $mean
}

proc ::wkcf::Destroy_Skin_Sphere_Group {what_dempack_package} {
    if {$what_dempack_package eq "G-DEMPack"} {
	if [GiD_Groups exists SKIN_SPHERE_DO_NOT_DELETE] {
	    GiD_Groups delete SKIN_SPHERE_DO_NOT_DELETE
	    GidUtils::EnableGraphics
	    GidUtils::UpdateWindow GROUPS
	    GidUtils::DisableGraphics
	}
    }
}

###############################################################################
#
#    NAME: wkcfdem.tcl
#
#    PURPOSE: Useful procedures to write the properties of DEM application to the Kratos *.mdpa
#
#    QUANTECH ATZ-DEVELOPMENT DEPARTMENT
#
#    AUTHOR : G. Socorro
#
#    CREATED AT: 01/10/13
#
#    HISTORY:
#
#     0.3- 03/11/13-G. Socorro, add the proc GetInletGroupNodes to get the inlet condition group node list
#     0.2- 02/11/13-G. Socorro, add the proc AssignSpecialBoundaries and GetBoundariesNodeList
#     0.1- 01/10/13-G. Socorro, create a base source code from wkcf.tcl
#
###############################################################################

proc ::wkcf::GetInletGroupMeshProperties {AppId} {
    # ABSTRACT: Get the inlet condition group properties mdpa file (only the nodes)
    variable dprops
   
    # For debug
    if {!$::wkcf::pflag} {
        set inittime [clock seconds]
    }

    # Write assigned group to the elements
    set meshgroupid 0

    # Init the global mesh identifier list
    set dprops($AppId,AllMeshId) [list]
    set dprops($AppId,AllMeshGroupId) [list]

    # Set the rootid
    set rootid "$AppId"
    set cproperty "dv"
    
    # Get the values
    set basexpath "$rootid//c.Conditions//c.DEM-Inlet"
    set gproplist [::xmlutils::setXmlContainerIds $basexpath]
    # wa "gproplist:$gproplist"
    foreach cgroupid $gproplist {
	# Get the group properties
	set cxpath "${basexpath}//c.[list ${cgroupid}]//c.MainProperties"
	set allgprop [::xmlutils::setXmlContainerPairs $cxpath "" "dv"]
	# wa "allgprop:$allgprop"
	if {[llength $allgprop]} {
	    if {[GiD_EntitiesGroups get $cgroupid nodes -count]} {
		incr meshgroupid 1
		# Create the meshid-group identifier mapping
		lappend dprops($AppId,AllMeshId) $meshgroupid
		# Set mesh-group id mapping
		if {$cgroupid ni $dprops($AppId,AllMeshGroupId)} {
		    lappend dprops($AppId,AllMeshGroupId) $cgroupid
		}
		set dprops($AppId,Mesh,$cgroupid,MeshIdGroup) $meshgroupid
		
		# Set the properties
		# Modify some properties
		set newprops [list]
	       	foreach propid $allgprop {
		    lassign $propid key value
		    if {($key eq "MatModel")||($key eq "DirectionVectorX") ||($key eq "DirectionVectorY") ||($key eq "DirectionVectorZ")} {
			continue
		    } else {
			if {$key eq "VelocityModulus"} {
			    set uservelmod $value
			    # Find velocity module 
			    set findvelmol [lsearch -index 0 $allgprop "VelocityModulus"]
			    # wa "findvelmol:$findvelmol"
			    set findvx [lsearch -index 0 $allgprop "DirectionVectorX"]
			    # wa "findvx:$findvx"
			    set findvy [lsearch -index 0 $allgprop "DirectionVectorY"]
			    # wa "findvy:$findvy"
			    set findvz [lsearch -index 0 $allgprop "DirectionVectorZ"]
			    # wa "findvz:$findvz"
			    if {($findvx !="-1") && ($findvy !="-1")&& ($findvz !="-1")} {
				# Get the velocity values
				set vxval [lindex $allgprop $findvx 1]
				set vyval [lindex $allgprop $findvy 1]
				set vzval [lindex $allgprop $findvz 1]
				set velvector [list $vxval $vyval $vzval] 
				# Calculate the modulus
				set velmod [::MathUtils::VectorModulus $velvector]
				# wa "velmod:$velmod"
				set endvelmodvalue "\[3\] (0.0,0.0,0.0)"
				if {$velmod !=0.0} {
				    set vxval [expr ($vxval*$uservelmod)/$velmod]
				    set vyval [expr ($vyval*$uservelmod)/$velmod]
				    set vzval [expr ($vzval*$uservelmod)/$velmod]
				    # wa "vxval:$vxval vyval:$vyval vzval:$vzval"
				    set endvelmodvalue "\[3\] ($vxval,$vyval,$vzval)"
				} 
				lappend newprops [list $key $endvelmodvalue]	
			    }
			} else {
			    lappend newprops [list $key $value]	
			}
		    }
		}
		# wa "newprops:$newprops"
		if {[llength $newprops]} {
		    set dprops($AppId,Mesh,$cgroupid,MeshIdGroupProp) $newprops
		}
	    }
	}
    }
    
    # For debug
    if {!$::wkcf::pflag} {
	set endtime [clock seconds]
	set ttime [expr $endtime-$inittime]
	# WarnWinText "endtime:$endtime ttime:$ttime"
	WarnWinText "Obtening DEM-Inlet group using the new mesh format: [::KUtils::Duration $ttime]"
    }   
}

proc ::wkcf::GetBoundariesNodeList {} {
    variable ndime

    set nlist [list]
    set cgroupid "-AKGDEMSkinMesh2D"
    if {$ndime =="3D"} {
	set cgroupid "-AKGDEMSkinMesh3D"
    }
    if {[GiD_EntitiesGroups get $cgroupid nodes -count]} {
	set nlist [GiD_EntitiesGroups get $cgroupid nodes]
    }
    return $nlist
}

proc ::wkcf::AssignSpecialBoundaries {ndime entitylist} {

    # wa "ndime:$ndime entitylist:$entitylist"
    set DEMApplication "No"
    set cproperty "dv"
    set cxpath "GeneralApplicationData//c.ApplicationTypes//i.DEM"
    catch { set DEMApplication [::xmlutils::setXml $cxpath $cproperty] }
    if {$DEMApplication eq "Yes"} {
	if { $ndime =="2D" } {        
	    # Automatic Kratos Group for all DEM boundary points
	    set groupid "-AKGDEMSkinMesh2D"
	    set entitytype "point" 
	    ::wkcf::CleanAutomaticConditionGroupGiD $entitytype $groupid        
	    # Get all end point list from the boundary lines
	    set endpointlist [list]
	    foreach lineid $entitylist {
		lappend endpointlist {*}[lrange [GiD_Geometry get line $lineid] 2 end]
	    }
	    set endpointlist [lsort -integer -unique $endpointlist]
	    # wa "endpointlist:$endpointlist"
	    # Assign the boundary condition
	    ::wkcf::AssignConditionToGroupGID $entitytype $endpointlist $groupid 

	} elseif { $ndime =="3D" } {

	    # Automatic Kratos Group for all DEM boundary lines
	    set groupid "-AKGDEMSkinMesh3D"
	    set entitytype "line" 
	    ::wkcf::CleanAutomaticConditionGroupGiD $entitytype $groupid
	    # Get all end line list from the boundary surfaces
	    set endlinelist [list]
	    foreach surfid $entitylist {
		set surfprop [GiD_Geometry get surface $surfid]
		# set surfacetype [lindex $surfprop 0]
		set nline [lindex $surfprop 2]
		# wa "surfacetype:$surfacetype\nnline:$nline\nsurfprop:$surfprop"
		set lineprop [list]
		#if {$surfacetype eq "nurbssurface"} {
		    set lineprop [lrange $surfprop 9 [expr 9+$nline-1]]
		#}
		foreach lprop $lineprop {
		    lassign $lprop lineid orientation 
		    lappend endlinelist $lineid
		}
	    }
	    # wa "before endlinelist:$endlinelist"
	    set endlinelist [lsort -integer -unique $endlinelist]
	    # wa "endlinelist:$endlinelist"
	    # Assign the boundary condition
	    ::wkcf::AssignConditionToGroupGID $entitytype $endlinelist $groupid
	}
    }
}


proc ::wkcf::GetInletGroupNodes {AppId cgroupid} {
    # ABSTRACT: Get the inlet condition group node list
    variable dprops
    
    set cprop [list]

    # For debug
    if {!$::wkcf::pflag} {
        set inittime [clock seconds]
    }
    
    # Set the rootid
    set rootid "$AppId"
    set cproperty "dv"
    
    # Get the values
    set basexpath "$rootid//c.Conditions//c.DEM-Inlet"
    # Get the group properties
    set cxpath "${basexpath}//c.[list ${cgroupid}]//c.MainProperties"
    set allgprop [::xmlutils::setXmlContainerPairs $cxpath "" "dv"]
    # wa "allgprop:$allgprop"
    if {[llength $allgprop]} {
	if {[GiD_EntitiesGroups get $cgroupid nodes -count]} {
	    # Check for the exclude boundaries option
	    set findeb [lsearch -index 0 $allgprop "ExcludeBoundaries"]
	    # wa "findeb:$findeb"
	    set ExcludeBoundaries "No"
	    if {$findeb !="-1"} {
		set ExcludeBoundaries [lindex $allgprop $findeb 1]
	    }
	    # wa "ExcludeBoundaries:$ExcludeBoundaries"
	    if {$ExcludeBoundaries eq "No"} {
		set cprop [GiD_EntitiesGroups get $cgroupid nodes]
	    } else {
		# Get the boundary node list 
		set nlist [::wkcf::GetBoundariesNodeList]
		foreach node_id [GiD_EntitiesGroups get $cgroupid nodes] {
		    set findnode [lsearch $nlist $node_id]
		    # wa "findnode:$findnode"
		    if {$findnode =="-1"} {
			lappend cprop $node_id
		    }
		}
	    }
	}
    }
    
    # For debug
    if {!$::wkcf::pflag} {
	set endtime [clock seconds]
	set ttime [expr $endtime-$inittime]
	WarnWinText "Get DEM-Inlet group nodes list: [::KUtils::Duration $ttime]"
    }
    return $cprop  
}

proc ::wkcf::WriteInletGroupMeshProperties {AppId} {
    # ABSTRACT: Write inlet condition group properties mdpa file (only the nodes)
    variable dprops
    variable filechannel
    
    # For debug
    if {!$::wkcf::pflag} {
        set inittime [clock seconds]
    }
    
    # Set the rootid
    set rootid "$AppId"
    
    # Get the values
    set basexpath "$rootid//c.Conditions//c.DEM-Inlet"
    set gproplist [::xmlutils::setXmlContainerIds $basexpath]
    # wa "gproplist:$gproplist"
    foreach cgroupid $gproplist {
	# Get the group node list
	set nlist [::wkcf::GetInletGroupNodes $AppId $cgroupid]
	if {[llength $nlist]} {
	    incr meshgroupid 1
	    # Get the meshid-group identifier mapping
	    set meshgroupid $dprops($AppId,Mesh,$cgroupid,MeshIdGroup) 
	    
	    # Write mesh properties for this group
	    GiD_File fprintf $filechannel "%s" "Begin Mesh $meshgroupid \/\/ GUI group identifier: $cgroupid"
	    # Write nodes
	    GiD_File fprintf $filechannel "%s" " "
	    GiD_File fprintf $filechannel "%s" " Begin MeshNodes"
	    foreach node_id $nlist {
		GiD_File fprintf $filechannel "%10i" $node_id
	    }
	    GiD_File fprintf $filechannel "%s" " End MeshNodes"
	    GiD_File fprintf $filechannel "%s" " "
	    GiD_File fprintf $filechannel "%s" "End Mesh"
	    GiD_File fprintf $filechannel "%s" ""
	}
    }
    
    # For debug
    if {!$::wkcf::pflag} {
	set endtime [clock seconds]
	set ttime [expr $endtime-$inittime]
	# WarnWinText "endtime:$endtime ttime:$ttime"
	WarnWinText "Write DEM-Inlet group using the new mesh format: [::KUtils::Duration $ttime]"
    }   
}

proc ::wkcf::WriteBoundingBoxDefaults {fileid} {

    puts $fileid "BoundingBoxMaxX                  = 3.00000e+00"  
    puts $fileid "BoundingBoxMaxY                  = 3.00000e+00"
    puts $fileid "BoundingBoxMaxZ                  = 2.00000e+01"
    puts $fileid "BoundingBoxMinX                  = -3.00000e+00"
    puts $fileid "BoundingBoxMinY                  = -3.00000e+00"
    puts $fileid "BoundingBoxMinZ                  = 0.00000e+00"

}

proc ::wkcf::WriteExplicitSolverVariables {} {
    # Write constitutive laws properties
    variable dprops;  variable ActiveAppList
    
    set AppId "DEM"
    set filename "DEM_explicit_solver_var.py"
    set PDir [::KUtils::GetPaths "PDir"]
    
    set fullname [file native [file join $PDir $filename]]
    
    # First delete the file
    if {[file exists $fullname]} {
	set res [file delete -force $fullname]
    }

    if { [catch { set fileid [open $fullname w+] }] } {
	WarnWin [= "Cannot write file %s. Permission denied" $fullname].
	return 0
    }

    # Write the group properties
    set cproperty "dv"
    set rootid $AppId

    puts $fileid ""

    puts $fileid "# Model Type"
    puts $fileid ""
    puts $fileid "ContinuumOption                  = \"OFF\""
    puts $fileid "RotationOption                   = \"ON\""
    puts $fileid "HomogeneousMaterialOption        = \"OFF\""
    puts $fileid "ElementType                      = \"SphericSwimmingParticle3D\""
    puts $fileid ""
    puts $fileid "# Meshing Settings"
    puts $fileid ""
    puts $fileid "CleanIndentationsOption          = \"OFF\""
    puts $fileid ""
    puts $fileid "# General Settings"
    puts $fileid ""
    puts $fileid "FinalTime                        = 1.50000e+02"
    puts $fileid "GravityX                         = 0.00000e+00"
    puts $fileid "GravityY                         = 0.00000e+00"
    puts $fileid "GravityZ                         = -9.81000e+00"

    # Get the use bounding box 
    set cxpath "$rootid//c.DEM-Options//c.Boundingbox//i.UseBoundingBox"
    set UseBoundingBox [::xmlutils::setXml $cxpath $cproperty]
    if {$UseBoundingBox eq "Active"} {
	puts $fileid "BoundingBoxOption                = \"ON\""

	# Get the bounding box type
	set cxpath "$rootid//c.DEM-Options//c.Boundingbox//i.BoundingBoxType"
	set BoundingBoxType [::xmlutils::setXml $cxpath $cproperty]
	if {$BoundingBoxType eq "Automatic"} {
	    puts $fileid "AutomaticBoundingBoxOption       = \"ON\""
	    # Get the enlargement factor
	    set cxpath "$rootid//c.DEM-Options//c.Boundingbox//i.EnlargementFactor"
	    set EnlargementFactor [::xmlutils::setXml $cxpath $cproperty]
	    puts $fileid "BoundingBoxEnlargementFactor     = $EnlargementFactor"
	
	    # Write default bounding box values
	    ::wkcf::WriteBoundingBoxDefaults $fileid
	   
	} elseif {$BoundingBoxType eq "Fixed"} {

	    puts $fileid "AutomaticBoundingBoxOption       = \"OFF\""
	    puts $fileid "BoundingBoxEnlargementFactor     = 1.0"

	    # Get the bounding limit
	    set varlist [list MaxX MaxY MaxZ MinX MinY MinZ]
	    foreach varid $varlist {
		set cxpath "$rootid//c.DEM-Options//c.Boundingbox//i.$varid"
		set $varid [::xmlutils::setXml $cxpath $cproperty]
		puts $fileid "BoundingBox$varid                  = [set $varid]"
	    }
	}
    } else {
	puts $fileid "BoundingBoxOption                = \"OFF\""
	puts $fileid "AutomaticBoundingBoxOption       = \"OFF\""
	puts $fileid "BoundingBoxEnlargementFactor     = 1.0"

	# Write default bounding box values
	::wkcf::WriteBoundingBoxDefaults $fileid
    }

  
    puts $fileid "Dimension                        = 3"
    puts $fileid "OutputFileType                   = \"Ascii\""
    puts $fileid "Multifile                        = \"single_file\""
    puts $fileid "PrintNeighbourLists              = \"OFF\""
    puts $fileid "ModelDataInfo                    = \"OFF\""
    puts $fileid ""
    puts $fileid "# Special features"
    puts $fileid ""
    puts $fileid "VirtualMassOption                = \"OFF\""
    puts $fileid "VirtualMassCoefficient           = 0.00000e+00"
    puts $fileid "MagicFactor                      = 1.00000e+00"
    puts $fileid "DeltaOption                      = \"OFF\""
    puts $fileid "SearchRadiusExtension            = 1.00000e-02"
    puts $fileid "AmplifiedSearchRadiusExtension   = 1.10000e+00"
    puts $fileid "FixVelocitiesOption              = \"OFF\""
    puts $fileid "TotalTimePercentageFixVelocities = 1.00000e+00"
    puts $fileid "TrihedronOption                  = \"OFF\""
    puts $fileid "LimitSurfaceOption               = 0"
    puts $fileid "SurfaceNormalDirX1               = 0.00000e+00"
    puts $fileid "SurfaceNormalDirY1               = 1.00000e+00"
    puts $fileid "SurfaceNormalDirZ1               = 0.00000e+00"
    puts $fileid "SurfacePointCoorX1               = 0.00000e+00"
    puts $fileid "SurfacePointCoorY1               = 0.00000e+00"
    puts $fileid "SurfacePointCoorZ1               = 0.00000e+00"
    puts $fileid "SurfaceFrictionAngle1            = 45"
    puts $fileid "SurfaceNormalDirX2               = 0.00000e+00"
    puts $fileid "SurfaceNormalDirY2               = 1.00000e+00"
    puts $fileid "SurfaceNormalDirZ2               = 0.00000e+00"
    puts $fileid "SurfacePointCoorX2               = 0.00000e+00"
    puts $fileid "SurfacePointCoorY2               = 0.00000e+00"
    puts $fileid "SurfacePointCoorZ2               = 0.00000e+00"
    puts $fileid "SurfaceFrictionAngle2            = 45"
    puts $fileid "SurfaceNormalDirX3               = 0.00000e+00"
    puts $fileid "SurfaceNormalDirY3               = 1.00000e+00"
    puts $fileid "SurfaceNormalDirZ3               = 0.00000e+00"
    puts $fileid "SurfacePointCoorX3               = 0.00000e+00"
    puts $fileid "SurfacePointCoorY3               = 0.00000e+00"
    puts $fileid "SurfacePointCoorZ3               = 0.00000e+00"
    puts $fileid "SurfaceFrictionAngle3            = 45"
    puts $fileid "SurfaceNormalDirX4               = 0.00000e+00"
    puts $fileid "SurfaceNormalDirY4               = 1.00000e+00"
    puts $fileid "SurfaceNormalDirZ4               = 0.00000e+00"
    puts $fileid "SurfacePointCoorX4               = 0.00000e+00"
    puts $fileid "SurfacePointCoorY4               = 0.00000e+00"
    puts $fileid "SurfacePointCoorZ4               = 0.00000e+00"
    puts $fileid "SurfaceFrictionAngle4            = 45"
    puts $fileid "SurfaceNormalDirX5               = 0.00000e+00"
    puts $fileid "SurfaceNormalDirY5               = 1.00000e+00"
    puts $fileid "SurfaceNormalDirZ5               = 0.00000e+00"
    puts $fileid "SurfacePointCoorX5               = 0.00000e+00"
    puts $fileid "SurfacePointCoorY5               = 0.00000e+00"
    puts $fileid "SurfacePointCoorZ5               = 0.00000e+00"
    puts $fileid "SurfaceFrictionAngle5            = 45"
    puts $fileid "LimitCylinderOption              = 0"
    puts $fileid "CylinderVelocity1                = 0.00000e+00"
    puts $fileid "CylinderAngularVelocity1         = 0.00000e+00"
    puts $fileid "CylinderInitialBaseCentreX1      = 0.00000e+00"
    puts $fileid "CylinderInitialBaseCentreY1      = 0.00000e+00"
    puts $fileid "CylinderInitialBaseCentreZ1      = 0.00000e+00"
    puts $fileid "CylinderAxisX1                   = 0.00000e+00"
    puts $fileid "CylinderAxisY1                   = 0.00000e+00"
    puts $fileid "CylinderAxisZ1                   = 1.00000e+00"
    puts $fileid "CylinderRadius1                  = 1.00000e+00"
    puts $fileid "CylinderFrictionAngle1           = 45"
    puts $fileid "CylinderVelocity2                = 0.00000e+00"
    puts $fileid "CylinderAngularVelocity2         = 0.00000e+00"
    puts $fileid "CylinderInitialBaseCentreX2      = 0.00000e+00"
    puts $fileid "CylinderInitialBaseCentreY2      = 0.00000e+00"
    puts $fileid "CylinderInitialBaseCentreZ2      = 0.00000e+00"
    puts $fileid "CylinderAxisX2                   = 0.00000e+00"
    puts $fileid "CylinderAxisY2                   = 0.00000e+00"
    puts $fileid "CylinderAxisZ2                   = 1.00000e+00"
    puts $fileid "CylinderRadius2                  = 2.50000e+00"
    puts $fileid "CylinderFrictionAngle2           = 45"
    puts $fileid "CylinderRadius3                  = 1.00000e+00"
    puts $fileid "CylinderFrictionAngle3           = 45"
    puts $fileid "CylinderRadius4                  = 1.00000e+00"
    puts $fileid "CylinderFrictionAngle4           = 45"
    puts $fileid "CylinderRadius5                  = 1.00000e+00"
    puts $fileid "CylinderFrictionAngle5           = 45"
    puts $fileid ""
    puts $fileid "# Time Discretization Settings" 
    puts $fileid ""
    puts $fileid "IntegrationScheme               = \"forward_euler\""
    puts $fileid "TimeStepsPerSearchStep           = 1.00000e+00"

    # Get automatic delta time type
    set cxpath "$rootid//c.DEM-Options//c.TimeStep//i.UseAutomaticDeltaTime"
    set UseAutomaticDeltaTime [::xmlutils::setXml $cxpath $cproperty]
    if {$UseAutomaticDeltaTime eq "Fixed"} {
	puts $fileid "AutoReductionOfTimeStepOption    = \"OFF\"" 
	
	# Get the value of DeltaTime
	set cxpath "$rootid//c.DEM-Options//c.TimeStep//i.DeltaTime"
	set DeltaTime [::xmlutils::setXml $cxpath $cproperty]
	puts $fileid "MaxTimeStep                      = $DeltaTime" 

	# Get the value of DeltaTimeSafetyFactor
	set cxpath "$rootid//c.DEM-Options//c.TimeStep//i.DeltaTimeSafetyFactor"
	set DeltaTimeSafetyFactor [::xmlutils::setXml $cxpath $cproperty]
	puts $fileid "DeltaTimeSafetyFactor            = $DeltaTimeSafetyFactor"

    } elseif {$UseAutomaticDeltaTime eq "Automatic"} {
	puts $fileid "AutoReductionOfTimeStepOption    = \"ON\""
    }
    
    puts $fileid "OutputTimeStep                   = 1.00000e-02"
    puts $fileid "ControlTime                      = 2.00000e+01"
    puts $fileid ""
    puts $fileid "# Material Model"
    puts $fileid ""
    puts $fileid "NormalForceCalculationType       = \"Linear\""
    puts $fileid "NormalDampingType                = \"ViscDamp\""
    puts $fileid "TangentialDampingType            = \"NoDamp\""
    puts $fileid "FailureCriterionType             = \"Uncoupled\""
    puts $fileid "TauZero                          = 5.26000e+00"
    puts $fileid "SigmaMax                         = 3.32000e+01"
    puts $fileid "SigmaMin                         = 3.32000e+00"
    puts $fileid "InternalFriction                 = 3.50000e+01"
    puts $fileid "NonLinearNormalElasticOption     = \"OFF\""
    puts $fileid "C1                               = 4.00000e-01"
    puts $fileid "N1                               = 2.00000e+01"
    puts $fileid "C2                               = 1.50000e+00"
    puts $fileid "N2                               = 1.50000e+01"
    puts $fileid "RotationalSpringOption           = \"OFF\""
    puts $fileid "RotaDampingType                  = \"NoDamp\""
    puts $fileid ""
    puts $fileid "# Global Material Parameters"
    puts $fileid ""
    puts $fileid "GeneralDensity                   = 1.00000e+02"
    puts $fileid "GeneralYoungModulus              = 2.10000e+07"
    puts $fileid "GeneralPoissonRatio              = 5.00000e-01"
    puts $fileid "GeneralCohesion                  = 4.16000e+06"
    puts $fileid "GeneralRollingFriction           = 0.00000e+00"
    puts $fileid "GeneralTension                   = 2.01000e+06"
    puts $fileid "GeneralRotaDampRatio             = 5.00000e-01"
    puts $fileid "GeneralStaticFrictionCoef        = 0.00000e+00"
    puts $fileid "GeneralDynamicFrictionCoef       = 0.00000e+00"
    puts $fileid "GeneralRestitutionCoef           = 5.00000e-01"
    puts $fileid "GeneralColour                    = 1.00000e+00"
    puts $fileid "GlobalVariablesOption            = \"OFF\""
    puts $fileid "GlobalKn                         = 3.00000e+03"
    puts $fileid "GlobalKt                         = 1.00000e+03"
    puts $fileid "GlobalKr                         = 1.00000e+03"
    puts $fileid "GlobalRn                         = 1.00000e+03"
    puts $fileid "GlobalRT                         = 1.00000e+03"
    puts $fileid "GlobalRr                         = 1.00000e+03"
    puts $fileid "GlobalFrictionAngle              = 4.00000e+01"
    puts $fileid ""
    puts $fileid "# Continuum Options"
    puts $fileid ""
    puts $fileid "StressStrainOperationsOption     = \"OFF\""
    puts $fileid "ContactMeshOption                = \"OFF\""
    puts $fileid "ConcreteTestOption               = \"OFF\""
    puts $fileid "RealTimeGraphOption              = \"OFF\""
    puts $fileid "TriaxialOption                   = \"OFF\""
    puts $fileid "ConfinementPressure              = 0.00000e+00"
    puts $fileid "InitialPressureAplicationTime    = 0.00000e+00"
    puts $fileid "TotalTimePercentAsForceAplTime   = 1.50000e+01"
    puts $fileid ""
    puts $fileid "#POSTPROCES"
    puts $fileid ""
    puts $fileid "PostVelocity                     = \"1\""
    puts $fileid "PostDisplacement                 = \"1\""
    puts $fileid "PostRadialDisplacement           = \"0\""
    puts $fileid "PostRHS                          = \"0\""
    puts $fileid "PostTotalForces                  = \"0\""
    puts $fileid "PostDampForces                   = \"0\""
    puts $fileid "PostAppliedForces                = \"0\""
    puts $fileid "PostRadius                       = \"1\""
    puts $fileid "PostParticleCohesion             = \"0\""
    puts $fileid "PostParticleTension              = \"0\""
    puts $fileid "PostGroupId                      = \"0\""
    puts $fileid "PostExportId                     = \"1\""
    puts $fileid "PostExportParticleFailureId      = \"0\""
    puts $fileid "PostExportSkinSphere             = \"0\""
    puts $fileid "PostLocalContactForceLow         = \"0\""
    puts $fileid "PostLocalContactForceHigh        = \"0\""
    puts $fileid "PostFailureCriterionState        = \"0\""
    puts $fileid "PostContactFailure               = \"0\""
    puts $fileid "PostContactTau                   = \"0\""
    puts $fileid "PostContactSigma                 = \"0\""
    puts $fileid "PostAngularVelocity              = \"0\""
    puts $fileid "PostParticleMoment               = \"0\""
    puts $fileid "PostEulerAngles                  = \"0\""
    puts $fileid "PostRepresentativeVolume         = \"0\""
    puts $fileid "PostMeanContactArea              = \"0\""
    puts $fileid "PostStressTensor                 = \"0\""
    puts $fileid ""
    puts $fileid "#FROM CND:"
    puts $fileid ""
    puts $fileid "PredefinedSkinOption             = \"OFF\""
    puts $fileid ""
    puts $fileid "TotalElementsVolume              = 1.13097e-01"
    puts $fileid ""
    puts $fileid "# Declare Python Variables"
    puts $fileid ""

    set PName [::KUtils::GetPaths "PName"]
    puts $fileid "problem_name=\"${PName}${AppId}\"" 
    puts $fileid "problem_path=\"[file join $PDir]\"" 
    
    # Get the kratos path 
    set cxpath "GeneralApplicationData//c.ProjectConfiguration//i.KratosPath"
    set cproperty "dv"
    set KratosPath [::xmlutils::setXml $cxpath $cproperty]
    set KratosPath [file native $KratosPath]
    
    # Write the kratos path
    puts $fileid "kratos_path=\"$[file join $KratosPath]\""

    close $fileid
}
















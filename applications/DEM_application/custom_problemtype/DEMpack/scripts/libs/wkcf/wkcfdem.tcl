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
#     0.6- 31/01/14-G. Socorro, correct a bug in the proc WriteDEMConditionProperties (write the value of IMPOSED_VELOCITY_X_VALUE, etc.)
#     0.5- 23/01/14-J. Garate, procs GetDEMGroupMeshProperties, ::wkcf::WriteDEMGroupMeshProperties, WriteDampingData,WriteMatTestData
#                              , WritePostProcessData, WriteMaterialModelData,WriteExplicitSolverVariables created or modified
#     0.4- 08/11/13-G. Socorro, modify some variable identifier (add DEM-*)
#     0.3- 03/11/13-G. Socorro, add the proc GetInletGroupNodes to get the inlet condition group node list
#     0.2- 02/11/13-G. Socorro, add the proc AssignSpecialBoundaries and GetBoundariesNodeList
#     0.1- 01/10/13-G. Socorro, create a base source code from wkcf.tcl
#
###############################################################################

proc ::wkcf::GetBoundariesNodeList {} {
    variable ndime

    set nlist [list]
    set cgroupid "-AKGDEMSkinMesh2D"
    if {$ndime =="3D"} {
	set cgroupid "-AKGDEMSkinMesh3D"
    }
    if { [GiD_Groups exists $cgroupid] } {
	if {[GiD_EntitiesGroups get $cgroupid nodes -count]} {
	    set nlist [GiD_EntitiesGroups get $cgroupid nodes]
	}
    }
    return $nlist
}

proc ::wkcf::GetElementType {} {

    variable ActiveAppList
    variable ndime
    global KPriv
    set gelemlist ""
    set elembasexpath "DEM//c.DEM-Elements//c.DEM-Element"
    set grouplist [::xmlutils::setXmlContainerIds $elembasexpath]

    if {$grouplist != ""} {
	foreach cgroup $grouplist {
	    if {[::xmlutils::setXml "$elembasexpath//c.[list ${cgroup}]//c.Properties//i.Material" dv] eq "DEM-IceMaterial"} {
		return [list "IceContPartDEMElement3D"]
	    }
	}
    }

    if {"Fluid" in $ActiveAppList} {
	set grouplist [::xmlutils::setXmlContainerIds $elembasexpath]
	if {$grouplist != ""} {
	    foreach cgroup $grouplist {
		set cproperty "dv"
		set cxpath "$elembasexpath//c.[list ${cgroup}]//c.Properties//i.Material"
		set material_name [::xmlutils::setXml $cxpath $cproperty]
		set cxpath "DEM//c.DEM-Elements//c.DEM-Element//c.[list ${cgroup}]//c.Properties//i.ElementType"
		set element_type [::xmlutils::setXml $cxpath "dv"]
		  if {$element_type eq "Nano"} {
		    set gelemlist [list "SwimmingNanoParticle"]
		    break
		} else {
		    set gelemlist [list "SwimmingDEMElement"]
		}
	    }
	}
    } else {
	if {$ndime eq "3D"} {
	    set gelemlist [list "SphericPartDEMElement3D"]
	    if {$KPriv(what_dempack_package) eq "C-DEMPack"} {
		set thermal_option [::xmlutils::setXml "DEM//c.DEM-Options//c.DEM-AdvancedOptions//i.DEM-Thermal" dv]
		if {$thermal_option eq "Yes"} {
		    set gelemlist [list "ThermalSphericContPartDEMElement3D"]
		} else {
		    set gelemlist [list "SphericContPartDEMElement3D"]
		}
	    }
	} else {
	    set gelemlist [list "CylinderPartDEMElement2D"]
	    if {$KPriv(what_dempack_package) eq "C-DEMPack"} {
		set cproperty "dv"
		set cxpath "DEM//c.DEM-Options//c.DEM-AdvancedOptions//i.DEM-Thermal"
		set thermal_option [::xmlutils::setXml $cxpath $cproperty]
		if {$thermal_option eq "Yes"} {
		    set gelemlist [list "ThermalCylinderContPartDEMElement2D"]
		} else {
		    set gelemlist [list "CylinderContPartDEMElement2D"]
		}
	    }
	}
    }
    return $gelemlist
}

proc ::wkcf::GetDEMActiveElements {} {

    variable ActiveAppList
    variable ndime
    global KPriv
    set gelemlist ""
    set elembasexpath "DEM//c.DEM-Elements//c.DEM-Element"

    if {"Fluid" in $ActiveAppList} {
	set grouplist [::xmlutils::setXmlContainerIds $elembasexpath]
	if {$grouplist != ""} {
	    foreach cgroup $grouplist {
		set cproperty "dv"
		set cxpath "$elembasexpath//c.[list ${cgroup}]//c.Properties//i.Material"
		set material_name [::xmlutils::setXml $cxpath $cproperty]
		set cxpath "DEM//c.DEM-Elements//c.DEM-Element//c.[list ${cgroup}]//c.Properties//i.ElementType"
		set element_type [::xmlutils::setXml $cxpath "dv"]
		  if {$element_type eq "Nano"} {
		    set gelemlist [list "SwimmingNanoParticle"]
		    break
		} else {
		    set gelemlist [list "SwimmingDEMElement"]
		}
	    }
	}
    } else {
	if {$ndime eq "3D"} {
	    set gelemlist [list "SphericPartDEMElement3D"]
	    if {$KPriv(what_dempack_package) eq "C-DEMPack"} {
		set thermal_option [::xmlutils::setXml "DEM//c.DEM-Options//c.DEM-AdvancedOptions//i.DEM-Thermal" dv]
		if {$thermal_option eq "Yes"} {
		    set gelemlist [list "ThermalSphericContPartDEMElement3D"]
		} else {
		    set gelemlist [list "SphericContPartDEMElement3D"]
		}
	    }
	} else {
	    set gelemlist [list "CylinderPartDEMElement2D"]
	    if {$KPriv(what_dempack_package) eq "C-DEMPack"} {
		set cproperty "dv"
		set cxpath "DEM//c.DEM-Options//c.DEM-AdvancedOptions//i.DEM-Thermal"
		set thermal_option [::xmlutils::setXml $cxpath $cproperty]
		if {$thermal_option eq "Yes"} {
		    set gelemlist [list "ThermalCylinderContPartDEMElement2D"]
		} else {
		    set gelemlist [list "CylinderContPartDEMElement2D"]
		}
	    }
	}
    }
    return $gelemlist
}

proc ::wkcf::AssignSpecialBoundaries {ndime entitylist} {
    set DEMApplication "No"
    catch {set DEMApplication [::xmlutils::setXml {GeneralApplicationData//c.ApplicationTypes//i.DEM} dv]}
    if {$DEMApplication eq "Yes"} {
	if {$ndime =="2D"} {
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
	    # Assign the boundary condition
	    ::wkcf::AssignConditionToGroupGID $entitytype $endpointlist $groupid
	} elseif {$ndime =="3D"} {
	    # Automatic Kratos Group for all DEM boundary lines
	    set groupid "-AKGDEMSkinMesh3D"
	    set entitytype "line"
	    ::wkcf::CleanAutomaticConditionGroupGiD $entitytype $groupid
	    # Get all end line list from the boundary surfaces
	    set endlinelist [list]
	    foreach surfid $entitylist {
		set surfprop [GiD_Geometry get surface $surfid]
		set surfacetype [lindex $surfprop 0]
		set nline [lindex $surfprop 2]
		set lineprop [list]
		if {$surfacetype eq "nurbssurface"} {
		    set lineprop [lrange $surfprop 9 [expr {9+$nline-1}]]
		} else {
			set lineprop [lrange $surfprop 3 [expr {3+$nline-1}]]
		}
		foreach lprop $lineprop {
		    lassign $lprop lineid orientation
		    lappend endlinelist $lineid
		}
	    }
	    set endlinelist [lsort -integer -unique $endlinelist]
		#W "endlinelist: $endlinelist"

	    # Assign the boundary condition
	    ::wkcf::AssignConditionToGroupGID $entitytype $endlinelist $groupid
	}
    }
}

proc ::wkcf::ForceTheMeshingOfDEMFEMWallGroups {} {
    foreach group_id [::xmlutils::setXmlContainerIds "DEM//c.DEM-Conditions//c.DEM-FEM-Wall"] {
		GiD_Process Mescape Meshing MeshCriteria Mesh Surfaces {*}[lindex [GiD_EntitiesGroups get $group_id all_geometry] 2] escape
    }
}

proc ::wkcf::ForceTheMeshingOfDEMInletGroups {} {
    foreach group_id [::xmlutils::setXmlContainerIds "DEM//c.DEM-Conditions//c.DEM-Inlet"] {
		GiD_Process Mescape Meshing MeshCriteria Mesh Surfaces {*}[lindex [GiD_EntitiesGroups get $group_id all_geometry] 2] escape
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
    set basexpath "$rootid//c.DEM-Conditions//c.DEM-Inlet"
    # Get the group properties
    set cxpath "${basexpath}//c.[list ${cgroupid}]//c.MainProperties"
    set allgprop [::xmlutils::setXmlContainerPairs $cxpath "" "dv"]
    if {[llength $allgprop]} {
	if {[GiD_EntitiesGroups get $cgroupid nodes -count]} {
	    # Check for the exclude boundaries option
	    set findeb [lsearch -index 0 $allgprop "ExcludeBoundaries"]
	    set ExcludeBoundaries "No"
	    if {$findeb !="-1"} {
		set ExcludeBoundaries [lindex $allgprop $findeb 1]
	    }
	    if {$ExcludeBoundaries eq "No"} {
		set cprop [GiD_EntitiesGroups get $cgroupid nodes]
	    } else {
		# Get the boundary node list
		set nlist [::wkcf::GetBoundariesNodeList]
		foreach node_id [GiD_EntitiesGroups get $cgroupid nodes] {
		    set findnode [lsearch $nlist $node_id]
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
	set ttime [expr {$endtime-$inittime}]
	WarnWinText "Get DEM-Inlet group nodes list: [::KUtils::Duration $ttime]"
    }
	#W "cprop: $cprop"
    return $cprop
}

proc ::wkcf::GetDemFemWallGroupNodes {cgroupid} {
    set list_of_nodes [list]
    if {[GiD_EntitiesGroups get $cgroupid elements -count]} {
	foreach element_id [GiD_EntitiesGroups get $cgroupid elements] {
	    # JG cambiar a GiD_Mesh
	    set info_from_gid [GiD_Mesh get element $element_id]
	    set element_type [lindex $info_from_gid 1]
	    if {$element_type eq "Sphere" || $element_type eq "Circle"} {
		GiD_EntitiesGroups unassign $cgroupid elements $element_id
		GiD_EntitiesGroups unassign $cgroupid nodes [lindex $info_from_gid 3]
		continue
	    }
	    if {$element_type != "Triangle" && $element_type != "Quadrilateral" && $element_type != "Linear"} {
		return [list 1 [= "Non-triangular, non-quadrilateral or non-linear elements found in a DEM-FEM Wall group. Please check the entities inside group '%s'." $cgroupid] {}]
	    }
	    # JG OJO REPETIDOS LSORT + FLAGS
	    lappend list_of_nodes {*}[lrange $info_from_gid 3 end]
	}
    }
    return [list 0 "" [lsort -increasing -integer -unique $list_of_nodes]]
}

proc ::wkcf::GetDSOLIDContactGroupNodes {AppId cgroupid} {
    variable dprops
    set cprop [list]

    set rootid DSOLID
    set cproperty "dv"
    # Get the values
    set basexpath "$rootid//c.DEM-Conditions//c.DEM-FEM-Solid-Wall"
    if {[GiD_EntitiesGroups get $cgroupid elements -count]} {
	set elems [GiD_EntitiesGroups get $cgroupid elements]
	foreach eid $elems {
	    set info_from_gid [GiD_Mesh get element $eid]
	    set elem_type [lrange $info_from_gid 1 1]
	    set nodes [lrange $info_from_gid 3 end]
	    lappend cprop {*}$nodes
	}
    }
    set nlist [lsort -increasing -integer -unique $cprop]
    return [list 0 "" $nlist]
}


proc ::wkcf::GetDSOLIDBoundaryConditionNodes {AppId cgroupid} {
    variable dprops
    set cprop [list]

    set rootid DSOLID
    set cproperty "dv"
    # Get the values
    set basexpath "$rootid//c.DEM-Conditions//c.Displacements"

    if {[GiD_EntitiesGroups get $cgroupid nodes -count]} {

      set nodes [GiD_EntitiesGroups get $cgroupid nodes]

    }
    set nlist [lsort -increasing -integer -unique $nodes]
    return [list 0 "" $nlist]
}


proc ::wkcf::GetDSOLIDGroupNodes {AppId cgroupid} {
    variable dprops

    set cprop [list]
    set rootid DSOLID
    set cproperty "dv"
    set basexpath "$rootid//c.Solid-Elements//c.Solid-Element"

    if {[GiD_EntitiesGroups get $cgroupid elements -count]} {
	set elems [GiD_EntitiesGroups get $cgroupid elements]
	foreach eid $elems {
	    set info_from_gid [GiD_Mesh get element $eid]
	    set nodes [lrange $info_from_gid 3 end]
	    lappend cprop {*}$nodes
	}
    }
    set nlist [lsort -increasing -integer -unique $cprop]
    return [list 0 "" $nlist]
}
###############################################################################

proc ::wkcf::GetDemGroupNodes {AppId cgroupid} {

    variable dprops

    set cprop [list]


    # Set the rootid
    set rootid "$AppId"
    set cproperty "dv"

    # Get the values
    set basexpath "$rootid//c.DEM-Conditions//c.DEM-FEM-Wall"
    # Get the group properties

    if {[GiD_EntitiesGroups get $cgroupid elements -count]} {
	set elems [GiD_EntitiesGroups get $cgroupid elements]
	foreach eid $elems {
	    set nodes [lrange [GiD_Mesh get element $eid] 3 end]
	    lappend cprop {*}$nodes

	}
    }

    return [lsort -increasing -integer -unique $cprop]
}

###############################################################################

proc ::wkcf::WriteDEMGroupMeshProperties {AppId} {
    # ABSTRACT: Write all DEM condition group properties mdpa file (only the nodes)
    variable filechannel
    variable ndime
    global KPriv
    # For debug
    if {!$::wkcf::pflag} {
	set inittime [clock seconds]
    }
    set basexpath "$AppId//c.DEM-Conditions//c.DEM-FEM-Wall"
    set gFEMproplist [::xmlutils::setXmlContainerIds $basexpath]
    variable dem_group_mesh_property_number
    set condlist [list "DEM-ForceIntegrationGroup" "DEM-TopLayerGroup" \
	    "DEM-BottomLayerGroup" "DEM-VelocityBC" "DEM-VelocityIC"]
    set pathlist [list "$AppId//c.DEM-Results//c.DEM-Graphs//c.DEM-ForceIntegrationGroup" \
	    "$AppId//c.DEM-MaterialTest//c.DEM-TopLayerGroup" \
	    "$AppId//c.DEM-MaterialTest//c.DEM-BottomLayerGroup" \
	    "$AppId//c.DEM-Conditions//c.DEM-VelocityBC" \
	    "$AppId//c.DEM-Conditions//c.DEM-VelocityIC"]
    foreach condid $condlist cpath $pathlist {
	set gproplist [::xmlutils::setXmlContainerIds $cpath]
	foreach cgroupid $gproplist {
	    if {$cgroupid ni $gFEMproplist} {
		if {[GiD_EntitiesGroups get $cgroupid elements -count]} {
		    set my_element [lindex [GiD_Mesh get element [lindex [GiD_EntitiesGroups get $cgroupid elements] 0] ] 1]
		    if {$my_element == "Sphere" || $my_element == "Circle"} {
		        #incr meshgroupid 1
		        set localpath "$cpath//c.$cgroupid"
		        set allgprop [::xmlutils::setXmlContainerPairs $localpath "" "dv"]
		        set motion_type ""
		        if {$condid eq "DEM-VelocityBC"} {
		            set cxpath "${localpath}//i.DEM-VelocityBCMotion"
		            set motion_type [::xmlutils::setXml $cxpath "dv"]
		        }
		        set is_condition_active "Yes"
		        if {$condid eq "DEM-VelocityBC"  || $condid eq "DEM-VelocityIC"} {
		            set cxpath "$localpath//i.SetActive"
		            set is_condition_active [::xmlutils::setXml $cxpath "dv"]
		        }
		        set dem_motion_local_path "$cpath//c.$cgroupid"
		        set cxpath "${dem_motion_local_path}//i.SetActive"
		        if {$condid eq "DEM-ForceIntegrationGroup"} {
		            set is_motion_activated "Yes"
		        } else {
		            set is_motion_activated [::xmlutils::setXml $cxpath "dv"]
		        }
		        set dem_elem_full_path "$AppId//c.DEM-Elements//c.DEM-Element//c.$cgroupid//c.Properties//i.SetActive"
		        set is_dem_element_active [::xmlutils::setXml $dem_elem_full_path dv]
		        if {$motion_type=="RigidBodyMotion" || $is_motion_activated=="No" || $is_condition_active=="No"|| $is_dem_element_active=="No"} {
		        #if {$motion_type=="RigidBodyMotion" || $is_motion_activated=="No" || $is_condition_active=="No"} {}
		            continue
		        }
		        incr dem_group_mesh_property_number
		        # Ai Vi FIG Top Bottom
		        GiD_File fprintf $filechannel "%s" "Begin SubModelPart $dem_group_mesh_property_number \/\/ GUI conditionid $condid group identifier: $cgroupid"
		        #::wkcf::WriteDEMConditionProperties $AppId $cgroupid $valuelist

		        GiD_File fprintf $filechannel "%s" "  Begin SubModelPartData"
		        ##############IMPOSED VELOCITIES########################
		        if {$condid eq "DEM-VelocityBC"} {
		            set active [::xmlutils::setXml "${localpath}//i.SetActive" "dv"]
		            if {$active eq "No"} { continue }
		            set motion_type [::xmlutils::setXml "${localpath}//i.DEM-VelocityBCMotion" "dv"]
		            if {$motion_type eq "RigidBodyMotion" } {
		        GiD_File fprintf $filechannel "  RIGID_BODY_MOTION 1"

		        GiD_File fprintf $filechannel "  LINEAR_VELOCITY [3] ([::xmlutils::setXml "${localpath}//c.LinearVelocity//i.LinearVelocityX" "dv"],[::xmlutils::setXml "${localpath}//c.LinearVelocity//i.LinearVelocityY" "dv"],[::xmlutils::setXml "${localpath}//c.LinearVelocity//i.LinearVelocityZ" "dv"])"

		        set linear_is_periodic [::xmlutils::setXml "${localpath}//c.LinearVelocity//i.LinearPeriodic" "dv"]
		        if {linear_is_periodic eq "Yes"} {
		            GiD_File fprintf $filechannel "  VELOCITY_PERIOD [::xmlutils::setXml "${localpath}//c.LinearVelocity//i.LinearPeriod" "dv"]"
		        } else {
		            GiD_File fprintf $filechannel "  VELOCITY_PERIOD 0.0"
		        }

		        GiD_File fprintf $filechannel "  VELOCITY_START_TIME [::xmlutils::setXml "${localpath}//c.LinearVelocity//i.LinearStartTime" "dv"]"
		        GiD_File fprintf $filechannel "  VELOCITY_STOP_TIME [::xmlutils::setXml "${localpath}//c.LinearVelocity//i.LinearEndTime" "dv"]"

		        GiD_File fprintf $filechannel "  ANGULAR_VELOCITY [3] ([::xmlutils::setXml "${localpath}//c.AngularVelocity//i.AngularVelocityX" "dv"],[::xmlutils::setXml "${localpath}//c.AngularVelocity//i.AngularVelocityY" "dv"],[::xmlutils::setXml "${localpath}//c.AngularVelocity//i.AngularVelocityZ" "dv"])"

		        GiD_File fprintf $filechannel "  ROTATION_CENTER [3] ([::xmlutils::setXml "${localpath}//c.AngularVelocity//i.CenterOfRotationX" "dv"],[::xmlutils::setXml "${localpath}//c.AngularVelocity//i.CenterOfRotationY" "dv"],[::xmlutils::setXml "${localpath}//c.AngularVelocity//i.CenterOfRotationZ" "dv"])"

		        set angular_is_periodic [::xmlutils::setXml "${localpath}//c.AngularVelocity//i.AngularPeriodic" "dv"]
		        if {angular_is_periodic eq "Yes"} {
		            GiD_File fprintf $filechannel "  ANGULAR_VELOCITY_PERIOD [::xmlutils::setXml "${localpath}//c.AngularVelocity//i.AngularPeriod" "dv"]"
		        } else {
		            GiD_File fprintf $filechannel "  ANGULAR_VELOCITY_PERIOD 0.0"
		        }

		        GiD_File fprintf $filechannel "  ANGULAR_VELOCITY_START_TIME [::xmlutils::setXml "${localpath}//c.AngularVelocity//i.AngularStartTime" "dv"]"
		        GiD_File fprintf $filechannel "  ANGULAR_VELOCITY_STOP_TIME [::xmlutils::setXml "${localpath}//c.AngularVelocity//i.AngularEndTime" "dv"]"
		    } else {
		        GiD_File fprintf $filechannel "  RIGID_BODY_MOTION 0"
		        GiD_File fprintf $filechannel "  VELOCITY_START_TIME [::xmlutils::setXml "${localpath}//i.VStart" "dv"]"
		        GiD_File fprintf $filechannel "  VELOCITY_STOP_TIME [::xmlutils::setXml "${localpath}//i.VEnd" "dv"]"

		        if { [::xmlutils::setXml "${localpath}//i.Ax" "dv"] eq "Yes" } {
		            GiD_File fprintf $filechannel "  IMPOSED_VELOCITY_X_VALUE [::xmlutils::setXml "${localpath}//i.Vx" "dv"]"
		        }
		        if { [::xmlutils::setXml "${localpath}//i.Ay" "dv"] eq "Yes" } {
		            GiD_File fprintf $filechannel "  IMPOSED_VELOCITY_Y_VALUE [::xmlutils::setXml "${localpath}//i.Vy" "dv"]"
		        }
		        if { [::xmlutils::setXml "${localpath}//i.Az" "dv"] eq "Yes" } {
		            GiD_File fprintf $filechannel "  IMPOSED_VELOCITY_Z_VALUE [::xmlutils::setXml "${localpath}//i.Vz" "dv"]"
		        }
		        if { [::xmlutils::setXml "${localpath}//i.Bx" "dv"] eq "Yes" } {
		            GiD_File fprintf $filechannel "  IMPOSED_ANGULAR_VELOCITY_X_VALUE [::xmlutils::setXml "${localpath}//i.AVx" "dv"]"
		        }
		        if { [::xmlutils::setXml "${localpath}//i.By" "dv"] eq "Yes" } {
		            GiD_File fprintf $filechannel "  IMPOSED_ANGULAR_VELOCITY_Y_VALUE [::xmlutils::setXml "${localpath}//i.AVy" "dv"]"
		        }
		        if { [::xmlutils::setXml "${localpath}//i.Bz" "dv"] eq "Yes" } {
		            GiD_File fprintf $filechannel "  IMPOSED_ANGULAR_VELOCITY_Z_VALUE [::xmlutils::setXml "${localpath}//i.AVz" "dv"]"
		        }
		    }

		##############INITIAL VELOCITIES########################
		        } elseif {$condid eq "DEM-VelocityIC"} {
		            if {[::xmlutils::setXml "${localpath}//i.SetActive" "dv"] eq "No"} { continue }
		    GiD_File fprintf $filechannel "  INITIAL_VELOCITY_X_VALUE [::xmlutils::setXml "${localpath}//c.Values//i.Vx" "dv"]"
		    GiD_File fprintf $filechannel "  INITIAL_VELOCITY_Y_VALUE [::xmlutils::setXml "${localpath}//c.Values//i.Vy" "dv"]"
		    GiD_File fprintf $filechannel "  INITIAL_VELOCITY_Z_VALUE [::xmlutils::setXml "${localpath}//c.Values//i.Vz" "dv"]"
		    GiD_File fprintf $filechannel "  INITIAL_ANGULAR_VELOCITY_X_VALUE [::xmlutils::setXml "${localpath}//c.Values//i.AVx" "dv"]"
		    GiD_File fprintf $filechannel "  INITIAL_ANGULAR_VELOCITY_Y_VALUE [::xmlutils::setXml "${localpath}//c.Values//i.AVy" "dv"]"
		    GiD_File fprintf $filechannel "  INITIAL_ANGULAR_VELOCITY_Z_VALUE [::xmlutils::setXml "${localpath}//c.Values//i.AVz" "dv"]"

		        } elseif {$condid eq "DEM-ForceIntegrationGroup"} {
		    set active [::xmlutils::setXml "${localpath}//i.SetActive" "dv"]
		            if {$active eq "No"} { continue }
		            GiD_File fprintf $filechannel "  FORCE_INTEGRATION_GROUP 1"
		            GiD_File fprintf $filechannel "  IDENTIFIER $cgroupid"
		        } elseif {$condid eq "DEM-TopLayerGroup"} {
		    set active [::xmlutils::setXml "${localpath}//i.SetActive" "dv"]
		            if {$active eq "No"} { continue }
		            GiD_File fprintf $filechannel "TOP 1"
		        } elseif {$condid eq "DEM-BottomLayerGroup"} {
		    set active [::xmlutils::setXml "${localpath}//i.SetActive" "dv"]
		            if {$active eq "No"} { continue }
		            GiD_File fprintf $filechannel "  BOTTOM 1"

		        }


		        GiD_File fprintf $filechannel "%s" "  End SubModelPartData"
		        # Write nodes
		        GiD_File fprintf $filechannel "%s" "  Begin SubModelPartNodes"
		        foreach eid [GiD_EntitiesGroups get $cgroupid elements] {
		            # JG cambiar a Mesh
		            if {$ndime eq "2D"} {
		                set node_id [lindex [GiD_Info Mesh Elements circle $eid $eid] 1]
		            } else {
		                set node_id [lindex [GiD_Info Mesh Elements sphere $eid $eid] 1]
		            }
		            GiD_File fprintf $filechannel "  $node_id"
		        }
		        GiD_File fprintf $filechannel "%s" "  End SubModelPartNodes"
		        GiD_File fprintf $filechannel "%s" "End SubModelPart"
		        GiD_File fprintf $filechannel "%s" ""
		    }
		}
	    }
	}
    }
    # Write Mesh Properties for DEM Element Groups
    ::wkcf::WriteDEMElementMeshProperties $AppId
    # Write Applied Loads condition
    ::wkcf::WriteAppliedLoadsData $AppId
    # Write Mesh Properties for DEM FEM Wall Conditions
    ::wkcf::WriteDEMFEMWallMeshProperties $AppId
    # For debug
    if {!$::wkcf::pflag} {
	set endtime [clock seconds]
	set ttime [expr {$endtime-$inittime}]
	WarnWinText "Write DEM-Mesh Properties group: [::KUtils::Duration $ttime]"
    }
}

proc ::wkcf::WriteAppliedLoadsData {AppId} {

    variable filechannel
    variable ndime
    variable dem_group_mesh_property_number

    set basepath "DEM//c.DEM-LoadConditions//c.BlastLoad"
    set cgrouplist [::xmlutils::setXmlContainerIds $basepath]
    set dem_applied_loads [ expr $dem_group_mesh_property_number ]

    foreach cgroupid $cgrouplist {

	incr dem_applied_loads
	GiD_File fprintf $filechannel "Begin SubModelPart $dem_applied_loads // GUI DEM-LoadConditions group identifier: $cgroupid"
	GiD_File fprintf $filechannel "  Begin SubModelPartData // GUI DEM-LoadConditions group identifier: $cgroupid"

	set basegrouppath "$basepath//c.$cgroupid//c.MainProperties"

	set BlR [::xmlutils::setXml "$basegrouppath//i.Radius" "dv"]
	GiD_File fprintf $filechannel "   BLAST_RADIUS            $BlR"

	set cxpath "$basegrouppath//i.Pressure_curve"
	set BlC [::xmlutils::setXml $cxpath "dv"]
	GiD_File fprintf $filechannel "   BLAST_CURVE             $BlC"

	set cxpath "$basegrouppath//i.Pressure_max"
	set BlP [::xmlutils::setXml $cxpath "dv"]
	GiD_File fprintf $filechannel "   BLAST_PRESSURE_MAX      $BlP"

	set cxpath "$basegrouppath//i.Time_max"
	set BlT [::xmlutils::setXml $cxpath "dv"]
	GiD_File fprintf $filechannel "   BLAST_TIME_PRESSURE_MAX $BlT"

	set cxpath "$basegrouppath//i.ShapeFactor"
	set BlSF [::xmlutils::setXml $cxpath "dv"]
	GiD_File fprintf $filechannel "   BLAST_SHAPE_FACTOR      $BlSF"

	set cxpath "$basegrouppath//i.Delay"
	set BlD [::xmlutils::setXml $cxpath "dv"]
	GiD_File fprintf $filechannel "   BLAST_TIME_DELAY        $BlD"

	set cxpath "$basegrouppath//i.BoreHole"
	set BlBH [::xmlutils::setXml $cxpath "dv"]
	GiD_File fprintf $filechannel "   BLAST_BOREHOLE          $BlBH"

	set icoord 0
	set boreholes [GiD_EntitiesGroups get $cgroupid point]
	set BlNBH [ llength $boreholes ]
	GiD_File fprintf $filechannel "   BLAST_NPOINTS           $BlNBH"

	foreach iborehole $boreholes {
	    incr icoord
	    set coordinates_string [lrange [GiD_Geometry get point $iborehole] 1 end ]
	    set coordX [lindex $coordinates_string 0]
	    set coordY [lindex $coordinates_string 1]
	    set coordZ [lindex $coordinates_string 2]
	    GiD_File fprintf $filechannel "   BLAST_COORDINATES_$icoord \[3\] ( $coordX, $coordY, $coordZ )"
	}

	GiD_File fprintf $filechannel "  End SubModelPartData"
	GiD_File fprintf $filechannel "  Begin SubModelPartNodes"
	GiD_File fprintf $filechannel "  End SubModelPartNodes"
	GiD_File fprintf $filechannel "End SubModelPart"
	GiD_File fprintf $filechannel ""

    }
}

proc ::wkcf::WriteBoundingBoxDefaults {fileid} {
	global KPriv

	if {$KPriv(what_dempack_package) eq "C-DEMPack"} {
		puts $fileid "BoundingBoxMaxX                  =  10.0"
		puts $fileid "BoundingBoxMaxY                  =  10.0"
		puts $fileid "BoundingBoxMaxZ                  =  10.0"
		puts $fileid "BoundingBoxMinX                  = -10.0"
		puts $fileid "BoundingBoxMinY                  = -10.0"
		puts $fileid "BoundingBoxMinZ                  = -10.0"
	} else {
		puts $fileid "BoundingBoxMaxX                  =  1.00000e+01"
		puts $fileid "BoundingBoxMaxY                  =  1.00000e+01"
		puts $fileid "BoundingBoxMaxZ                  =  1.00000e+01"
		puts $fileid "BoundingBoxMinX                  = -1.00000e+01"
		puts $fileid "BoundingBoxMinY                  = -1.00000e+01"
		puts $fileid "BoundingBoxMinZ                  = -1.00000e+01"
	}
}

proc ::wkcf::WriteBoundingBoxDefaultsInJsonFile {fileid} {

	global KPriv

	if {$KPriv(what_dempack_package) eq "C-DEMPack"} {
		puts $fileid "\"BoundingBoxMaxX\"                  : 10.0,"
		puts $fileid "\"BoundingBoxMaxY\"                  : 10.0,"
		puts $fileid "\"BoundingBoxMaxZ\"                  : 10.0,"
		puts $fileid "\"BoundingBoxMinX\"                  : -10.0,"
		puts $fileid "\"BoundingBoxMinY\"                  : -10.0,"
		puts $fileid "\"BoundingBoxMinZ\"                  : -10.0,"
	} else {
		puts $fileid "\"BoundingBoxMaxX\"                  :  1.00000e+01,"
		puts $fileid "\"BoundingBoxMaxY\"                  :  1.00000e+01,"
		puts $fileid "\"BoundingBoxMaxZ\"                  :  1.00000e+01,"
		puts $fileid "\"BoundingBoxMinX\"                  : -1.00000e+01,"
		puts $fileid "\"BoundingBoxMinY\"                  : -1.00000e+01,"
		puts $fileid "\"BoundingBoxMinZ\"                  : -1.00000e+01,"
	}
}

proc ::wkcf::WriteMatTestData {fileid} {
    global KPriv
    set cxpath "DEM//c.DEM-Options//c.DEM-AdvancedOptions//i.DEM-MaterialTestActivate"
    set TestTypeOn [::xmlutils::setXml $cxpath "dv"]
    set cxpath "DEM//c.DEM-MaterialTest//i.DEM-TestType"
    set TestType [::xmlutils::setXml $cxpath "dv"]
    puts $fileid ""
    if {$TestTypeOn eq "Yes" && $KPriv(what_dempack_package) eq "C-DEMPack"} {
	puts $fileid "\"TestType\"                       : \"$TestType\","
    } else {
	puts $fileid "\"TestType\"                       : \"None\","
    }
    set LVelt 0.0
    set LVelb 0.0
    set cxpath "DEM//c.DEM-MaterialTest//i.DEM-ConfinementPressure"
    set ConfPress [::xmlutils::setXml $cxpath "dv"]
    puts $fileid "\"ConfinementPressure\"              : $ConfPress,"

    set basexpath "DEM//c.DEM-MaterialTest//c.DEM-TopLayerGroup"
    set topgroup [::xmlutils::setXmlContainerIds $basexpath]
    if {[llength $topgroup]} {
    set basexpath "DEM//c.DEM-Conditions//c.DEM-FEM-Wall"
    set gproplist [::xmlutils::setXmlContainerIds $basexpath]
    if {[llength $gproplist]} {
        set cxpath "DEM//c.DEM-Conditions//c.DEM-FEM-Wall//c.[lindex $topgroup 0]//c.LinearVelocity//i.LinearVelocityY"
        set LVelt [::xmlutils::setXml $cxpath "dv"]
    }
    }

    set basexpath "DEM//c.DEM-MaterialTest//c.DEM-BotLayerGroup"
    set botgroup [::xmlutils::setXmlContainerIds $basexpath]
    if {[llength $botgroup]} {
    set basexpath "DEM//c.DEM-Conditions//c.DEM-FEM-Wall"
    set gproplist [::xmlutils::setXmlContainerIds $basexpath]
    if {[llength $gproplist]} {
        set cxpath "DEM//c.DEM-Conditions//c.DEM-FEM-Wall//c.[lindex $botgroup 0]//c.LinearVelocity//i.LinearVelocityY"
        set LVelb [::xmlutils::setXml $cxpath "dv"]
    }
    }

    if {$TestTypeOn eq "No"} {set LVelt 0.0}
    if {$TestTypeOn eq "No"} {set LVelb 0.0}
    puts $fileid "\"LoadingVelocity\"              : [expr ($LVelt-$LVelb)],"
    set cxpath "DEM//c.DEM-MaterialTest//i.DEM-MeshType"
    set mt [::xmlutils::setXml $cxpath "dv"]
    puts $fileid "\"MeshType\"                        : \"$mt\","
    set cxpath "DEM//c.DEM-MaterialTest//i.DEM-MeshPath"
    set mp [::xmlutils::setXml $cxpath "dv"]
    if {$mp eq ""} {set mp "0"}
    puts $fileid "\"MeshPath\"                        : \"$mp\","
    set cxpath "DEM//c.DEM-MaterialTest//i.DEM-ProbeLength"
    set pl [::xmlutils::setXml $cxpath "dv"]
    puts $fileid "\"SpecimenLength\"                  : $pl,"
    set cxpath "DEM//c.DEM-MaterialTest//i.DEM-ProbeDiameter"
    set pd [::xmlutils::setXml $cxpath "dv"]
    puts $fileid "\"SpecimenDiameter\"                : $pd,"
    set ms [expr {3.14159265359*$pd*$pd/4}]
    puts $fileid "\"MeasuringSurface\"                 : $ms,"
}

proc ::wkcf::TranslateToBinary {yes_no_var} {
    if {($yes_no_var eq "1") || ($yes_no_var eq "Yes") || ($yes_no_var eq "Active") || ($yes_no_var eq "Activate")} {
	return 1
    } else {
	return 0}
}

proc ::wkcf::TranslateToBinaryJson {yes_no_var} {
    if { [::wkcf::TranslateToBinary $yes_no_var] == 1 } {
	return true
    } else {
	return false
    }
}

proc ::wkcf::WritePostProcessData {fileid} {
    variable ActiveAppList
    global KPriv
    puts $fileid "# PostProcess Results"

    # GraphExportFreq
    set cxpath "DEM//c.DEM-Results//c.DEM-Graphs//i.DEM-GraphExportFreq"
    set GEF [::xmlutils::setXml $cxpath "dv"]
    puts $fileid "GraphExportFreq                  = $GEF"

    # VelocityTrapGraphExportFreq
    set cxpath "DEM//c.DEM-Results//c.DEM-VelocityTrap//i.DEM-VelTrapGraphExportFreq"
    set VTGEF [::xmlutils::setXml $cxpath "dv"]
    puts $fileid "VelTrapGraphExportFreq           = $VTGEF"

    # Output Time Step
    if {"Fluid" in $ActiveAppList} {
	set cxpath "GeneralApplicationData//c.SimulationOptions//i.OutputDeltaTime"
	set OTS [::xmlutils::setXml $cxpath "dv"]
	puts $fileid "OutputTimeStep                   = $OTS"
    } else {
	set cxpath "DEM//c.DEM-Results//i.DEM-OTimeStepType"
	set OTimeStepType [::xmlutils::setXml $cxpath "dv"]

	if {$OTimeStepType eq "Detail_priority"} {
	    set cxpath "DEM//c.DEM-Results//i.DEM-OTimeStepDetail"
	    set OTimeStepDetail [::xmlutils::setXml $cxpath "dv"]
	    puts $fileid "OutputTimeStep                   = $OTimeStepDetail"
	}

	if {$OTimeStepType eq "Storage_priority"} {

	    set cxpath "DEM//c.DEM-SolutionStrategy//c.DEM-TimeParameters//i.DEM-TotalTime"
	    set TTime [::xmlutils::setXml $cxpath "dv"]
	    set cxpath "DEM//c.DEM-Results//i.DEM-OTimeStepStorage"
	    set amount [::xmlutils::setXml $cxpath "dv"]
	    set OTimeStepStorage [expr {double($TTime)/$amount}]
	    set cxpath "DEM//c.DEM-SolutionStrategy//c.DEM-TimeParameters//i.DeltaTime"
	    set MaxTimeStep [::xmlutils::setXml $cxpath "dv"]
	    set maxamount [expr {$TTime/$MaxTimeStep}]
	    if {$amount < $maxamount} {
		puts $fileid "OutputTimeStep                   = $OTimeStepStorage"
	    } else {
		puts $fileid "OutputTimeStep                   = $MaxTimeStep" }
	}
    }

    set cxpath "DEM//c.DEM-Options//c.DEM-Boundingbox//i.PrintBoundingBox"
    set PrintOrNot [::xmlutils::setXml $cxpath "dv"]
    puts $fileid "PostBoundingBox                  = [::wkcf::TranslateToBinary $PrintOrNot]"

    # Size distribution curve
    if {$KPriv(what_dempack_package) eq "C-DEMPack"} {
    set cproperty "dv"
    set cxpath "DEM//c.DEM-Results//i.DEM-Granulometry"
    set granulometry_option [::xmlutils::setXml $cxpath $cproperty]
    if {$granulometry_option eq "Yes"} {
    puts $fileid "Granulometry                     = \"[::xmlutils::setXml {DEM//c.DEM-Results//i.DEM-Granulometry} dv]\""
    }
    }

    # PostDisplacement
    if {"Fluid" in $ActiveAppList} {
	set cxpathtoDEMresults "GeneralApplicationData//c.SDEM-Results//c.ResultsToPrint//c.DEMResults"
	set cxpathtoFLUIDresults "GeneralApplicationData//c.SDEM-Results//c.ResultsToPrint//c.FluidResults"
    } else {
	set cxpathtoDEMresults "DEM//c.DEM-Results//c.DEM-PartElem"
    }

    set cxpath "$cxpathtoDEMresults//i.DEM-Displacement"
    set PrintOrNot [::xmlutils::setXml $cxpath "dv"]
    puts $fileid "PostDisplacement                 = [::wkcf::TranslateToBinary $PrintOrNot]"

    # PostVelocity
    set cxpath "$cxpathtoDEMresults//i.DEM-PostVel"
    set PrintOrNot [::xmlutils::setXml $cxpath "dv"]
    puts $fileid "PostVelocity                     = [::wkcf::TranslateToBinary $PrintOrNot]"

    puts $fileid "# DEM only Results"

    # PostTotalForces
    set cxpath "$cxpathtoDEMresults//i.DEM-TotalForces"
    set PrintOrNot [::xmlutils::setXml $cxpath "dv"]
    puts $fileid "PostTotalForces                  = [::wkcf::TranslateToBinary $PrintOrNot]"

    # PostRigidElementForces
    set cxpath "$cxpathtoDEMresults//i.DEM-RigidElementForces"
    set PrintOrNot [::xmlutils::setXml $cxpath "dv"]
    puts $fileid "PostRigidElementForces           = [::wkcf::TranslateToBinary $PrintOrNot]"

    if {$KPriv(what_dempack_package) eq "C-DEMPack"} {
	# PostExportSkinSphere
	set PrintOrNot [::xmlutils::setXml "$cxpathtoDEMresults//i.DEM-SkinSphere" dv]
	puts $fileid "PostSkinSphere                   = [::wkcf::TranslateToBinary $PrintOrNot]"
	# PostPoissonRatio
	puts $fileid "PostPoissonRatio                 = [::wkcf::TranslateToBinary [::xmlutils::setXml $cxpathtoDEMresults//i.DEM-PoissonRatio dv]]"

	set cproperty "dv"
	set cxpath "DEM//c.DEM-Options//c.DEM-AdvancedOptions//i.DEM-Thermal"
	set thermal_option [::xmlutils::setXml $cxpath $cproperty]
	if {$thermal_option eq "Yes"} {
	    # PostThermalTemp
	    set cxpath "DEM//c.DEM-Results//c.DEM-ThermalPost//i.DEM-ParticleTemperature"
	    set PrintOrNot [::xmlutils::setXml $cxpath "dv"]
	    puts $fileid "PostTemperature                  = [::wkcf::TranslateToBinary $PrintOrNot]"
	    # PostThermalFlux
	    set cxpath "DEM//c.DEM-Results//c.DEM-ThermalPost//i.DEM-ParticleHeatFlux"
	    set PrintOrNot [::xmlutils::setXml $cxpath "dv"]
	    puts $fileid "PostHeatFlux                     = [::wkcf::TranslateToBinary $PrintOrNot]"
	}

	# PostPrintVirtualSeaSurface
	if {[::wkcf::TranslateToBinary [::xmlutils::setXml "$cxpathtoDEMresults//i.DEM-SeaSurface" dv]]} {
	    puts $fileid "PostVirtualSeaSurfaceX1          = [::xmlutils::setXml "$cxpathtoDEMresults//i.DEM-VirtualSeaSurfaceX1" dv]"
	    puts $fileid "PostVirtualSeaSurfaceY1          = [::xmlutils::setXml "$cxpathtoDEMresults//i.DEM-VirtualSeaSurfaceY1" dv]"
	    puts $fileid "PostVirtualSeaSurfaceX2          = [::xmlutils::setXml "$cxpathtoDEMresults//i.DEM-VirtualSeaSurfaceX2" dv]"
	    puts $fileid "PostVirtualSeaSurfaceY2          = [::xmlutils::setXml "$cxpathtoDEMresults//i.DEM-VirtualSeaSurfaceY2" dv]"
	    puts $fileid "PostVirtualSeaSurfaceX3          = [::xmlutils::setXml "$cxpathtoDEMresults//i.DEM-VirtualSeaSurfaceX3" dv]"
	    puts $fileid "PostVirtualSeaSurfaceY3          = [::xmlutils::setXml "$cxpathtoDEMresults//i.DEM-VirtualSeaSurfaceY3" dv]"
	    puts $fileid "PostVirtualSeaSurfaceX4          = [::xmlutils::setXml "$cxpathtoDEMresults//i.DEM-VirtualSeaSurfaceX4" dv]"
	    puts $fileid "PostVirtualSeaSurfaceY4          = [::xmlutils::setXml "$cxpathtoDEMresults//i.DEM-VirtualSeaSurfaceY4" dv]"
	}
    }

    # Radius
    set cxpath "$cxpathtoDEMresults//i.DEM-Radius"
    set PrintOrNot [::xmlutils::setXml $cxpath "dv"]
    puts $fileid "PostRadius                       = [::wkcf::TranslateToBinary $PrintOrNot]"

    # Calculate Rotations
    set cxpath "DEM//c.DEM-Options//c.DEM-AdvancedOptions//i.DEM-CalculateRotations"
    set useRotationOption [::xmlutils::setXml $cxpath "dv"]
    if {$useRotationOption == "Yes"} {
	# PostAngularVelocity
	set cxpath "$cxpathtoDEMresults//i.DEM-AngularVelocity"
	set PrintOrNot [::xmlutils::setXml $cxpath "dv"]
	puts $fileid "PostAngularVelocity              = [::wkcf::TranslateToBinary $PrintOrNot]"

	# PostParticleMoment
	set cxpath "$cxpathtoDEMresults//i.DEM-ParticleMoment"
	set PrintOrNot [::xmlutils::setXml $cxpath "dv"]
	puts $fileid "PostParticleMoment               = [::wkcf::TranslateToBinary $PrintOrNot]"

	# PostEulerAngles
	set cxpath "$cxpathtoDEMresults//i.DEM-EulerAngles"
	set PrintOrNot [::xmlutils::setXml $cxpath "dv"]
	puts $fileid "PostEulerAngles                  = [::wkcf::TranslateToBinary $PrintOrNot]"

	set cxpath "DEM//c.DEM-Options//c.DEM-AdvancedOptions//i.DEM-RollingFriction"
	set rf [::xmlutils::setXml $cxpath "dv"]
	if {$rf == "Yes"} {
	    # PostRollingResistanceMoment
	    set cxpath "$cxpathtoDEMresults//i.DEM-RollingResistanceMoment"
	    set PrintOrNot [::xmlutils::setXml $cxpath "dv"]
	    puts $fileid "PostRollingResistanceMoment      = [::wkcf::TranslateToBinary $PrintOrNot]"
	} else {
	    puts $fileid "PostRollingResistanceMoment      = 0"
	}
    } else {
	puts $fileid "PostAngularVelocity              = 0"
	puts $fileid "PostParticleMoment               = 0"
	puts $fileid "PostEulerAngles                  = 0"
	puts $fileid "PostRollingResistanceMoment      = 0"
    }

    puts $fileid "# FEM only Results"
    # PostElasticForces
    set cxpath "$cxpathtoDEMresults//i.DEM-ElasForces"
    set PrintOrNot [::xmlutils::setXml $cxpath "dv"]
    puts $fileid "PostElasticForces                = [::wkcf::TranslateToBinary $PrintOrNot]"

    # PostContactForces
    set cxpath "$cxpathtoDEMresults//i.DEM-ContactForces"
    set PrintOrNot [::xmlutils::setXml $cxpath "dv"]
    puts $fileid "PostContactForces                = [::wkcf::TranslateToBinary $PrintOrNot]"

    # PostTangentialElasticForces
    set cxpath "$cxpathtoDEMresults//i.DEM-TangElasForces"
    set PrintOrNot [::xmlutils::setXml $cxpath "dv"]
    puts $fileid "PostTangentialElasticForces      = [::wkcf::TranslateToBinary $PrintOrNot]"

    # PostShearStress
    set cxpath "$cxpathtoDEMresults//i.DEM-ShearStress"
    set PrintOrNot [::xmlutils::setXml $cxpath "dv"]
    puts $fileid "PostShearStress                  = [::wkcf::TranslateToBinary $PrintOrNot]"

    # PostPressure
    set cxpath "$cxpathtoDEMresults//i.DEM-Pressure"
    set PrintOrNot [::xmlutils::setXml $cxpath "dv"]
    puts $fileid "PostPressure                     = [::wkcf::TranslateToBinary $PrintOrNot]"

    puts $fileid "# FEM_wear only Results"
    # PostWear
    set cxpath "$cxpathtoDEMresults//i.DEM-Wear"
    set PrintOrNot [::xmlutils::setXml $cxpath "dv"]
    puts $fileid "PostNonDimensionalVolumeWear     = [::wkcf::TranslateToBinary $PrintOrNot]"

    # PostNodalArea
    set cxpath "$cxpathtoDEMresults//i.DEM-NodalArea"
    set PrintOrNot [::xmlutils::setXml $cxpath "dv"]
    puts $fileid "PostNodalArea                    = [::wkcf::TranslateToBinary $PrintOrNot]"

    puts $fileid "# Results on bond elements"
    if {$KPriv(what_dempack_package) eq "C-DEMPack"} {
		# PostStressStrainOnSpheres
		puts $fileid "PostStressStrainOption           = [::wkcf::TranslateToBinary [::xmlutils::setXml $cxpathtoDEMresults//i.DEM-Stresses dv]]"

		# Write all Dem Bond Elem Properties
		set basexpath "DEM//c.DEM-Results//c.DEM-BondElem"
		set ilist [::xmlutils::setXmlContainerIds $basexpath "Item"]
		set kxpath "Applications/DEM"

		set cxpath "DEM//c.DEM-Results//c.DEM-BondElem//i.DEM-ContactSigma"
		set PrintOrNot [::xmlutils::setXml $cxpath "dv"]
		puts $fileid "PostContactSigma                 = [::wkcf::TranslateToBinary $PrintOrNot]"

		set cxpath "DEM//c.DEM-Results//c.DEM-BondElem//i.DEM-ContactTau"
		set PrintOrNot [::xmlutils::setXml $cxpath "dv"]
		puts $fileid "PostContactTau                   = [::wkcf::TranslateToBinary $PrintOrNot]"

		set cxpath "DEM//c.DEM-Results//c.DEM-BondElem//i.DEM-LocalContactForce"
		set PrintOrNot [::xmlutils::setXml $cxpath "dv"]
		puts $fileid "PostLocalContactForce            = [::wkcf::TranslateToBinary $PrintOrNot]"

		set cxpath "DEM//c.DEM-Results//c.DEM-BondElem//i.DEM-FailureCrit"
		set PrintOrNot [::xmlutils::setXml $cxpath "dv"]
		puts $fileid "PostFailureCriterionState        = [::wkcf::TranslateToBinary $PrintOrNot]"

		set cxpath "DEM//c.DEM-Results//c.DEM-BondElem//i.DEM-Failureid"
		set PrintOrNot [::xmlutils::setXml $cxpath "dv"]
		puts $fileid "PostContactFailureId             = [::wkcf::TranslateToBinary $PrintOrNot]"

		set cxpath "DEM//c.DEM-Results//c.DEM-BondElem//i.DEM-MeanContactArea"
		set PrintOrNot [::xmlutils::setXml $cxpath "dv"]
		puts $fileid "PostMeanContactArea              = [::wkcf::TranslateToBinary $PrintOrNot]"
	}
    puts $fileid "# Under revision"
    # PostRHS
    set cxpath "$cxpathtoDEMresults//i.DEM-Rhs"
    set PrintOrNot [::xmlutils::setXml $cxpath "dv"]
    puts $fileid "PostRHS                          = [::wkcf::TranslateToBinary $PrintOrNot]"

    # DampForces
    set cxpath "$cxpathtoDEMresults//i.DEM-DampForces"
    set PrintOrNot [::xmlutils::setXml $cxpath "dv"]
    puts $fileid "PostDampForces                   = [::wkcf::TranslateToBinary $PrintOrNot]"

    # AppliedForces
    set cxpath "$cxpathtoDEMresults//i.DEM-AppliedForces"
    set PrintOrNot [::xmlutils::setXml $cxpath "dv"]
    puts $fileid "PostAppliedForces                = [::wkcf::TranslateToBinary $PrintOrNot]"

    # PostGroupId
    set cxpath "$cxpathtoDEMresults//i.DEM-GroupId"
    set PrintOrNot [::xmlutils::setXml $cxpath "dv"]
    puts $fileid "PostGroupId                      = [::wkcf::TranslateToBinary $PrintOrNot]"

    # PostExportId
    set cxpath "$cxpathtoDEMresults//i.DEM-ExportId"
    set PrintOrNot [::xmlutils::setXml $cxpath "dv"]
    puts $fileid "PostExportId                     = [::wkcf::TranslateToBinary $PrintOrNot]"



    # SWIMMING-SPECIFIC SECTION BEGINS ###########################################################################

    if {"Fluid" in $ActiveAppList} {
	puts $fileid "# Swimming DEM-specific section begins"
	puts $fileid "#-------------------------------------"

	# PostPressure
	# DEM variables
	set cxpath "$cxpathtoFLUIDresults//i.Fluid-Pressure"
	set PrintOrNot [::xmlutils::setXml $cxpath "dv"]
	puts $fileid "PostFluidPressure                          = [::wkcf::TranslateToBinary $PrintOrNot]"
	set cxpath "DEM//c.DEM-Options//c.DEM-Boundingbox//i.PrintBoundingBox"
	set PBB [::xmlutils::setXml $cxpath "dv"]
	set cxpath "$cxpathtoDEMresults//i.DEM-ReynoldsNumber"
	set PrintOrNot [::xmlutils::setXml $cxpath "dv"]
	puts $fileid "print_REYNOLDS_NUMBER_option               = [::wkcf::TranslateToBinary $PrintOrNot]"
	set cxpath "$cxpathtoDEMresults//i.DEM-PressureGradProj"
	set PrintOrNot [::xmlutils::setXml $cxpath "dv"]
	puts $fileid "print_PRESSURE_GRAD_PROJECTED_option       = [::wkcf::TranslateToBinary $PrintOrNot]"
	set cxpath "$cxpathtoDEMresults//i.DEM-FluidVelProj"
	set PrintOrNot [::xmlutils::setXml $cxpath "dv"]
	puts $fileid "print_FLUID_VEL_PROJECTED_option           = [::wkcf::TranslateToBinary $PrintOrNot]"
	set cxpath "$cxpathtoDEMresults//i.DEM-FluidAccelProj"
	set PrintOrNot [::xmlutils::setXml $cxpath "dv"]
	puts $fileid "print_FLUID_ACCEL_PROJECTED_option         = [::wkcf::TranslateToBinary $PrintOrNot]"
	set cxpath "$cxpathtoDEMresults//i.DEM-Buoyancy"
	set PrintOrNot [::xmlutils::setXml $cxpath "dv"]
	puts $fileid "print_BUOYANCY_option                      = [::wkcf::TranslateToBinary $PrintOrNot]"
	set cxpath "$cxpathtoDEMresults//i.DEM-Drag"
	set PrintOrNot [::xmlutils::setXml $cxpath "dv"]
	puts $fileid "print_DRAG_FORCE_option                    = [::wkcf::TranslateToBinary $PrintOrNot]"
	set cxpath "$cxpathtoDEMresults//i.DEM-VirtualMass"
	set PrintOrNot [::xmlutils::setXml $cxpath "dv"]
	puts $fileid "print_VIRTUAL_MASS_FORCE_option            = [::wkcf::TranslateToBinary $PrintOrNot]"
	set cxpath "$cxpathtoDEMresults//i.DEM-Basset"
	set PrintOrNot [::xmlutils::setXml $cxpath "dv"]
	puts $fileid "print_BASSET_FORCE_option                  = [::wkcf::TranslateToBinary $PrintOrNot]"
	set cxpath "$cxpathtoDEMresults//i.DEM-Lift"
	set PrintOrNot [::xmlutils::setXml $cxpath "dv"]
	puts $fileid "print_LIFT_FORCE_option                    = [::wkcf::TranslateToBinary $PrintOrNot]"
	set cxpath "$cxpathtoDEMresults//i.DEM-FluidVelRate"
	set PrintOrNot [::xmlutils::setXml $cxpath "dv"]
	puts $fileid "print_FLUID_VEL_PROJECTED_RATE_option      = [::wkcf::TranslateToBinary $PrintOrNot]"
	set cxpath "$cxpathtoDEMresults//i.DEM-FluidViscosity"
	set PrintOrNot [::xmlutils::setXml $cxpath "dv"]
	puts $fileid "print_FLUID_VISCOSITY_PROJECTED_option     = [::wkcf::TranslateToBinary $PrintOrNot]"
	set cxpath "$cxpathtoDEMresults//i.DEM-FluidFractionProj"
	set PrintOrNot [::xmlutils::setXml $cxpath "dv"]
	puts $fileid "print_FLUID_FRACTION_PROJECTED_option      = [::wkcf::TranslateToBinary $PrintOrNot]"
	set cxpath "$cxpathtoDEMresults//i.DEM-VelocityLaplacianProjected"
	set PrintOrNot [::xmlutils::setXml $cxpath "dv"]
	puts $fileid "print_FLUID_VEL_LAPL_PROJECTED_option      = [::wkcf::TranslateToBinary $PrintOrNot]"
	set cxpath "$cxpathtoDEMresults//i.DEM-VelocityLaplacianRateProjected"
	set PrintOrNot [::xmlutils::setXml $cxpath "dv"]
	puts $fileid "print_FLUID_VEL_LAPL_RATE_PROJECTED_option = [::wkcf::TranslateToBinary $PrintOrNot]"
	set cxpath "$cxpathtoDEMresults//i.DEM-HydroForce"
	set PrintOrNot [::xmlutils::setXml $cxpath "dv"]
	puts $fileid "print_HYDRODYNAMIC_FORCE_option            = [::wkcf::TranslateToBinary $PrintOrNot]"
	set cxpath "$cxpathtoDEMresults//i.DEM-HydroMoment"
	set PrintOrNot [::xmlutils::setXml $cxpath "dv"]
	puts $fileid "print_HYDRODYNAMIC_MOMENT_option           = [::wkcf::TranslateToBinary $PrintOrNot]"

	# Fluid variables
    set cxpath "$cxpathtoFLUIDresults//i.Fluid-AugmentedVel"
	set PrintOrNot [::xmlutils::setXml $cxpath "dv"]
	puts $fileid "print_MESH_VELOCITY1_option                = [::wkcf::TranslateToBinary $PrintOrNot]"
    set cxpath "$cxpathtoFLUIDresults//i.Fluid-BodyForce"
	set PrintOrNot [::xmlutils::setXml $cxpath "dv"]
	puts $fileid "print_BODY_FORCE_option                    = [::wkcf::TranslateToBinary $PrintOrNot]"
	set cxpath "$cxpathtoFLUIDresults//i.Fluid-FluidFraction"
	set PrintOrNot [::xmlutils::setXml $cxpath "dv"]
	puts $fileid "print_FLUID_FRACTION_option                = [::wkcf::TranslateToBinary $PrintOrNot]"
	set cxpath "$cxpathtoFLUIDresults//i.Fluid-FluidFractionGrad"
	set PrintOrNot [::xmlutils::setXml $cxpath "dv"]
	puts $fileid "print_FLUID_FRACTION_GRADIENT_option       = [::wkcf::TranslateToBinary $PrintOrNot]"
	set cxpath "$cxpathtoFLUIDresults//i.Fluid-HydroReaction"
	set PrintOrNot [::xmlutils::setXml $cxpath "dv"]
	puts $fileid "print_HYDRODYNAMIC_REACTION_option         = [::wkcf::TranslateToBinary $PrintOrNot]"
	set cxpath "$cxpathtoFLUIDresults//i.Fluid-Pressure"
	set PrintOrNot [::xmlutils::setXml $cxpath "dv"]
	puts $fileid "print_PRESSURE_option                      = [::wkcf::TranslateToBinary $PrintOrNot]"
	set cxpath "$cxpathtoFLUIDresults//i.Fluid-PressureGrad"
	set PrintOrNot [::xmlutils::setXml $cxpath "dv"]
	puts $fileid "print_PRESSURE_GRADIENT_option             = [::wkcf::TranslateToBinary $PrintOrNot]"
	set cxpath "$cxpathtoFLUIDresults//i.Fluid-SolidFraction"
	set PrintOrNot [::xmlutils::setXml $cxpath "dv"]
	puts $fileid "print_DISPERSE_FRACTION_option                = [::wkcf::TranslateToBinary $PrintOrNot]"
    set cxpath "$cxpathtoFLUIDresults//i.Fluid-MeanHydroReaction"
    set PrintOrNot [::xmlutils::setXml $cxpath "dv"]
    puts $fileid "print_MEAN_HYDRODYNAMIC_REACTION_option    = [::wkcf::TranslateToBinary $PrintOrNot]"
	set cxpath "$cxpathtoFLUIDresults//i.Fluid-VelocityLaplacian"
	set PrintOrNot [::xmlutils::setXml $cxpath "dv"]
	puts $fileid "print_VELOCITY_LAPLACIAN_option            = [::wkcf::TranslateToBinary $PrintOrNot]"
	set cxpath "$cxpathtoFLUIDresults//i.Fluid-VelocityLaplacianRate"
	set PrintOrNot [::xmlutils::setXml $cxpath "dv"]
	puts $fileid "print_VELOCITY_LAPLACIAN_RATE_option       = [::wkcf::TranslateToBinary $PrintOrNot]"
	puts $fileid ""
	puts $fileid "#-------------------------------------"
	puts $fileid "# Swimming DEM-specific section ends"
	puts $fileid ""
	# SWIMMING-SPECIFIC SECTION ENDS ###########################################################################
    }
}


proc ::wkcf::WritePostProcessDataForJson {fileid} {

    variable ActiveAppList
    global KPriv

    # GraphExportFreq
    set cxpath "DEM//c.DEM-Results//c.DEM-Graphs//i.DEM-GraphExportFreq"
    set GEF [::xmlutils::setXml $cxpath "dv"]
    puts $fileid "\"GraphExportFreq\"                  : $GEF,"

    # VelocityTrapGraphExportFreq
    set cxpath "DEM//c.DEM-Results//c.DEM-VelocityTrap//i.DEM-VelTrapGraphExportFreq"
    set VTGEF [::xmlutils::setXml $cxpath "dv"]
    puts $fileid "\"VelTrapGraphExportFreq\"           : $VTGEF,"

    # Output Time Step
    if {"Fluid" in $ActiveAppList} {
	set cxpath "GeneralApplicationData//c.SimulationOptions//i.OutputDeltaTime"
	set OTS [::xmlutils::setXml $cxpath "dv"]
	puts $fileid "\"OutputTimeStep\"                   : $OTS,"
    } else {
	set cxpath "DEM//c.DEM-Results//i.DEM-OTimeStepType"
	set OTimeStepType [::xmlutils::setXml $cxpath "dv"]

	if {$OTimeStepType eq "Detail_priority"} {
	    set cxpath "DEM//c.DEM-Results//i.DEM-OTimeStepDetail"
	    set OTimeStepDetail [::xmlutils::setXml $cxpath "dv"]
	    puts $fileid "\"OutputTimeStep\"                   : $OTimeStepDetail,"
	}

	if {$OTimeStepType eq "Storage_priority"} {

	    set cxpath "DEM//c.DEM-SolutionStrategy//c.DEM-TimeParameters//i.DEM-TotalTime"
	    set TTime [::xmlutils::setXml $cxpath "dv"]
	    set cxpath "DEM//c.DEM-Results//i.DEM-OTimeStepStorage"
	    set amount [::xmlutils::setXml $cxpath "dv"]
	    set OTimeStepStorage [expr {double($TTime)/$amount}]
	    set cxpath "DEM//c.DEM-SolutionStrategy//c.DEM-TimeParameters//i.DeltaTime"
	    set MaxTimeStep [::xmlutils::setXml $cxpath "dv"]
	    set maxamount [expr {$TTime/$MaxTimeStep}]
	    if {$amount < $maxamount} {
		puts $fileid "\"OutputTimeStep\"                   : $OTimeStepStorage,"
	    } else {
		puts $fileid "\"OutputTimeStep\"                   : $MaxTimeStep," }
	}
    }

    set cxpath "DEM//c.DEM-Options//c.DEM-Boundingbox//i.PrintBoundingBox"
    set PrintOrNot [::xmlutils::setXml $cxpath "dv"]
    puts $fileid "\"PostBoundingBox\"                  : [::wkcf::TranslateToBinaryJson $PrintOrNot],"

    # Size distribution curve
    if {$KPriv(what_dempack_package) eq "C-DEMPack"} {
    set cproperty "dv"
    set cxpath "DEM//c.DEM-Results//i.DEM-Granulometry"
    set granulometry_option [::xmlutils::setXml $cxpath $cproperty]
    if {$granulometry_option eq "Yes"} {
    puts $fileid "\"Granulometry\"                     : \"[::xmlutils::setXml {DEM//c.DEM-Results//i.DEM-Granulometry} dv]\","
    }
    }

    # PostDisplacement
    if {"Fluid" in $ActiveAppList} {
	set cxpathtoDEMresults "GeneralApplicationData//c.SDEM-Results//c.ResultsToPrint//c.DEMResults"
	set cxpathtoFLUIDresults "GeneralApplicationData//c.SDEM-Results//c.ResultsToPrint//c.FluidResults"
    } else {
	set cxpathtoDEMresults "DEM//c.DEM-Results//c.DEM-PartElem"
    }

    set cxpath "$cxpathtoDEMresults//i.DEM-Displacement"
    set PrintOrNot [::xmlutils::setXml $cxpath "dv"]
    puts $fileid "\"PostDisplacement\"                 : [::wkcf::TranslateToBinaryJson $PrintOrNot],"

    # PostVelocity
    set cxpath "$cxpathtoDEMresults//i.DEM-PostVel"
    set PrintOrNot [::xmlutils::setXml $cxpath "dv"]
    puts $fileid "\"PostVelocity\"                     : [::wkcf::TranslateToBinaryJson $PrintOrNot],"

    # PostTotalForces
    set cxpath "$cxpathtoDEMresults//i.DEM-TotalForces"
    set PrintOrNot [::xmlutils::setXml $cxpath "dv"]
    puts $fileid "\"PostTotalForces\"                  : [::wkcf::TranslateToBinaryJson $PrintOrNot],"

    # PostRigidElementForces
    set cxpath "$cxpathtoDEMresults//i.DEM-RigidElementForces"
    set PrintOrNot [::xmlutils::setXml $cxpath "dv"]
    puts $fileid "\"PostRigidElementForces\"           : [::wkcf::TranslateToBinaryJson $PrintOrNot],"

    if {$KPriv(what_dempack_package) eq "C-DEMPack"} {
	# PostExportSkinSphere
	set PrintOrNot [::xmlutils::setXml "$cxpathtoDEMresults//i.DEM-SkinSphere" dv]
	puts $fileid "\"PostSkinSphere\"                   : [::wkcf::TranslateToBinaryJson $PrintOrNot],"
	# PostPoissonRatio
	puts $fileid "\"PostPoissonRatio\"                 : [::wkcf::TranslateToBinaryJson [::xmlutils::setXml $cxpathtoDEMresults//i.DEM-PoissonRatio dv]],"

	set cproperty "dv"
	set cxpath "DEM//c.DEM-Options//c.DEM-AdvancedOptions//i.DEM-Thermal"
	set thermal_option [::xmlutils::setXml $cxpath $cproperty]
	if {$thermal_option eq "Yes"} {
	    # PostThermalTemp
	    set cxpath "DEM//c.DEM-Results//c.DEM-ThermalPost//i.DEM-ParticleTemperature"
	    set PrintOrNot [::xmlutils::setXml $cxpath "dv"]
	    puts $fileid "\"PostTemperature\"                  : [::wkcf::TranslateToBinaryJson $PrintOrNot],"
	    # PostThermalFlux
	    set cxpath "DEM//c.DEM-Results//c.DEM-ThermalPost//i.DEM-ParticleHeatFlux"
	    set PrintOrNot [::xmlutils::setXml $cxpath "dv"]
	    puts $fileid "\"PostHeatFlux\"                     : [::wkcf::TranslateToBinaryJson $PrintOrNot],"
	}

	# PostPrintVirtualSeaSurface
	if {[::wkcf::TranslateToBinaryJson [::xmlutils::setXml "$cxpathtoDEMresults//i.DEM-SeaSurface" dv]]} {
	    puts $fileid "\"virtual_sea_surface_settings\"     : {"
	    puts $fileid "\"print_sea_surface\"                : true,"
	    puts $fileid "\"PostVirtualSeaSurfaceX1\"          : [::xmlutils::setXml "$cxpathtoDEMresults//i.DEM-VirtualSeaSurfaceX1" dv],"
	    puts $fileid "\"PostVirtualSeaSurfaceY1\"          : [::xmlutils::setXml "$cxpathtoDEMresults//i.DEM-VirtualSeaSurfaceY1" dv],"
	    puts $fileid "\"PostVirtualSeaSurfaceX2\"          : [::xmlutils::setXml "$cxpathtoDEMresults//i.DEM-VirtualSeaSurfaceX2" dv],"
	    puts $fileid "\"PostVirtualSeaSurfaceY2\"          : [::xmlutils::setXml "$cxpathtoDEMresults//i.DEM-VirtualSeaSurfaceY2" dv],"
	    puts $fileid "\"PostVirtualSeaSurfaceX3\"          : [::xmlutils::setXml "$cxpathtoDEMresults//i.DEM-VirtualSeaSurfaceX3" dv],"
	    puts $fileid "\"PostVirtualSeaSurfaceY3\"          : [::xmlutils::setXml "$cxpathtoDEMresults//i.DEM-VirtualSeaSurfaceY3" dv],"
	    puts $fileid "\"PostVirtualSeaSurfaceX4\"          : [::xmlutils::setXml "$cxpathtoDEMresults//i.DEM-VirtualSeaSurfaceX4" dv],"
	    puts $fileid "\"PostVirtualSeaSurfaceY4\"          : [::xmlutils::setXml "$cxpathtoDEMresults//i.DEM-VirtualSeaSurfaceY4" dv]"
	    puts $fileid "},"
	}
    }

    # Radius
    set cxpath "$cxpathtoDEMresults//i.DEM-Radius"
    set PrintOrNot [::xmlutils::setXml $cxpath "dv"]
    puts $fileid "\"PostRadius\"                       : [::wkcf::TranslateToBinaryJson $PrintOrNot],"

    # Calculate Rotations
    set cxpath "DEM//c.DEM-Options//c.DEM-AdvancedOptions//i.DEM-CalculateRotations"
    set useRotationOption [::xmlutils::setXml $cxpath "dv"]
    if {$useRotationOption == "Yes"} {
	# PostAngularVelocity
	set cxpath "$cxpathtoDEMresults//i.DEM-AngularVelocity"
	set PrintOrNot [::xmlutils::setXml $cxpath "dv"]
	puts $fileid "\"PostAngularVelocity\"              : [::wkcf::TranslateToBinaryJson $PrintOrNot],"

	# PostParticleMoment
	set cxpath "$cxpathtoDEMresults//i.DEM-ParticleMoment"
	set PrintOrNot [::xmlutils::setXml $cxpath "dv"]
	puts $fileid "\"PostParticleMoment\"               : [::wkcf::TranslateToBinaryJson $PrintOrNot],"

	# PostEulerAngles
	set cxpath "$cxpathtoDEMresults//i.DEM-EulerAngles"
	set PrintOrNot [::xmlutils::setXml $cxpath "dv"]
	puts $fileid "\"PostEulerAngles\"                  : [::wkcf::TranslateToBinaryJson $PrintOrNot],"

	set cxpath "DEM//c.DEM-Options//c.DEM-AdvancedOptions//i.DEM-RollingFriction"
	set rf [::xmlutils::setXml $cxpath "dv"]
	if {$rf == "Yes"} {
	    # PostRollingResistanceMoment
	    set cxpath "$cxpathtoDEMresults//i.DEM-RollingResistanceMoment"
	    set PrintOrNot [::xmlutils::setXml $cxpath "dv"]
	    puts $fileid "\"PostRollingResistanceMoment\"      : [::wkcf::TranslateToBinaryJson $PrintOrNot],"
	} else {
	    puts $fileid "\"PostRollingResistanceMoment\"      : false,"
	}
    } else {
	puts $fileid "\"PostAngularVelocity\"              : false,"
	puts $fileid "\"PostParticleMoment\"               : false,"
	puts $fileid "\"PostEulerAngles\"                  : false,"
	puts $fileid "\"PostRollingResistanceMoment\"      : false,"
    }
    # PostCharacteristicLength
    if {$KPriv(what_dempack_package) ne "F-DEMPack"} {
	set PrintOrNot [::xmlutils::setXml "$cxpathtoDEMresults//i.DEM-CharacteristicLength" dv]
	puts $fileid "\"PostCharacteristicLength\"         : [::wkcf::TranslateToBinaryJson $PrintOrNot],"
    }
    puts $fileid ""
    # PostElasticForces
    set cxpath "$cxpathtoDEMresults//i.DEM-ElasForces"
    set PrintOrNot [::xmlutils::setXml $cxpath "dv"]
    puts $fileid "\"PostElasticForces\"                : [::wkcf::TranslateToBinaryJson $PrintOrNot],"

    # PostContactForces
    set cxpath "$cxpathtoDEMresults//i.DEM-ContactForces"
    set PrintOrNot [::xmlutils::setXml $cxpath "dv"]
    puts $fileid "\"PostContactForces\"                : [::wkcf::TranslateToBinaryJson $PrintOrNot],"

    # PostTangentialElasticForces
    set cxpath "$cxpathtoDEMresults//i.DEM-TangElasForces"
    set PrintOrNot [::xmlutils::setXml $cxpath "dv"]
    puts $fileid "\"PostTangentialElasticForces\"      : [::wkcf::TranslateToBinaryJson $PrintOrNot],"

    # PostShearStress
    set cxpath "$cxpathtoDEMresults//i.DEM-ShearStress"
    set PrintOrNot [::xmlutils::setXml $cxpath "dv"]
    puts $fileid "\"PostShearStress\"                  : [::wkcf::TranslateToBinaryJson $PrintOrNot],"

    # PostReactions
	if {$KPriv(what_dempack_package) eq "C-DEMpack"} {
		set cxpath "$cxpathtoDEMresults//i.DEM-Reactions"
		set PrintOrNot [::xmlutils::setXml $cxpath "dv"]
		puts $fileid "\"PostReactions\"                    : [::wkcf::TranslateToBinaryJson [::xmlutils::setXml "$cxpathtoDEMresults//i.DEM-Reactions" "dv"]],"
	}

    # PostPressure
    set cxpath "$cxpathtoDEMresults//i.DEM-Pressure"
    set PrintOrNot [::xmlutils::setXml $cxpath "dv"]
    puts $fileid "\"PostPressure\"                     : [::wkcf::TranslateToBinaryJson $PrintOrNot],"

    puts $fileid ""
    # PostWear
    set cxpath "$cxpathtoDEMresults//i.DEM-Wear"
    set PrintOrNot [::xmlutils::setXml $cxpath "dv"]
    puts $fileid "\"PostNonDimensionalVolumeWear\"     : [::wkcf::TranslateToBinaryJson $PrintOrNot],"

    # PostNodalArea
    set cxpath "$cxpathtoDEMresults//i.DEM-NodalArea"
    set PrintOrNot [::xmlutils::setXml $cxpath "dv"]
    puts $fileid "\"PostNodalArea\"                    : [::wkcf::TranslateToBinaryJson $PrintOrNot],"

    puts $fileid ""
    if {$KPriv(what_dempack_package) eq "C-DEMPack"} {
	# PostStressStrainOnSpheres
	puts $fileid "\"PostStressStrainOption\"           : [::wkcf::TranslateToBinaryJson [::xmlutils::setXml $cxpathtoDEMresults//i.DEM-Stresses dv]],"
   }

   if {$KPriv(what_dempack_package) eq "C-DEMPack"} {
       # Write all Dem Bond Elem Properties
       set basexpath "DEM//c.DEM-Results//c.DEM-BondElem"
       set ilist [::xmlutils::setXmlContainerIds $basexpath "Item"]
       set kxpath "Applications/DEM"

       set cxpath "DEM//c.DEM-Results//c.DEM-BondElem//i.DEM-ContactSigma"
       set PrintOrNot [::xmlutils::setXml $cxpath "dv"]
       puts $fileid "\"PostContactSigma\"                 : [::wkcf::TranslateToBinaryJson $PrintOrNot],"

       set cxpath "DEM//c.DEM-Results//c.DEM-BondElem//i.DEM-ContactTau"
       set PrintOrNot [::xmlutils::setXml $cxpath "dv"]
       puts $fileid "\"PostContactTau\"                   : [::wkcf::TranslateToBinaryJson $PrintOrNot],"

       set cxpath "DEM//c.DEM-Results//c.DEM-BondElem//i.DEM-LocalContactForce"
       set PrintOrNot [::xmlutils::setXml $cxpath "dv"]
       puts $fileid "\"PostLocalContactForce\"            : [::wkcf::TranslateToBinaryJson $PrintOrNot],"

       set cxpath "DEM//c.DEM-Results//c.DEM-BondElem//i.DEM-FailureCrit"
       set PrintOrNot [::xmlutils::setXml $cxpath "dv"]
       puts $fileid "\"PostFailureCriterionState\"        : [::wkcf::TranslateToBinaryJson $PrintOrNot],"

       set cxpath "DEM//c.DEM-Results//c.DEM-BondElem//i.DEM-Failureid"
       set PrintOrNot [::xmlutils::setXml $cxpath "dv"]
       puts $fileid "\"PostContactFailureId\"             : [::wkcf::TranslateToBinaryJson $PrintOrNot],"

       set cxpath "DEM//c.DEM-Results//c.DEM-BondElem//i.DEM-MeanContactArea"
       set PrintOrNot [::xmlutils::setXml $cxpath "dv"]
       puts $fileid "\"PostMeanContactArea\"              : [::wkcf::TranslateToBinaryJson $PrintOrNot],"
   }

    # PostRHS
    set cxpath "$cxpathtoDEMresults//i.DEM-Rhs"
    set PrintOrNot [::xmlutils::setXml $cxpath "dv"]
    puts $fileid "\"PostRHS\"                          : [::wkcf::TranslateToBinaryJson $PrintOrNot],"

    # DampForces
    set cxpath "$cxpathtoDEMresults//i.DEM-DampForces"
    set PrintOrNot [::xmlutils::setXml $cxpath "dv"]
    puts $fileid "\"PostDampForces\"                   : [::wkcf::TranslateToBinaryJson $PrintOrNot],"

    # AppliedForces
    set cxpath "$cxpathtoDEMresults//i.DEM-AppliedForces"
    set PrintOrNot [::xmlutils::setXml $cxpath "dv"]
    puts $fileid "\"PostAppliedForces\"                : [::wkcf::TranslateToBinaryJson $PrintOrNot],"

    # PostGroupId
    set cxpath "$cxpathtoDEMresults//i.DEM-GroupId"
    set PrintOrNot [::xmlutils::setXml $cxpath "dv"]
    puts $fileid "\"PostGroupId\"                      : [::wkcf::TranslateToBinaryJson $PrintOrNot],"

    # PostExportId
    set cxpath "$cxpathtoDEMresults//i.DEM-ExportId"
    set PrintOrNot [::xmlutils::setXml $cxpath "dv"]
    puts $fileid "\"PostExportId\"                     : [::wkcf::TranslateToBinaryJson $PrintOrNot],"


    # SWIMMING-SPECIFIC SECTION BEGINS ###########################################################################

    if {"Fluid" in $ActiveAppList} {
	puts $fileid ""

	# PostPressure
	# DEM variables
	set cxpath "$cxpathtoFLUIDresults//i.Fluid-Pressure"
	set PrintOrNot [::xmlutils::setXml $cxpath "dv"]
	puts $fileid "\"PostFluidPressure\"                          : [::wkcf::TranslateToBinaryJson $PrintOrNot],"
	set cxpath "DEM//c.DEM-Options//c.DEM-Boundingbox//i.PrintBoundingBox"
	set PBB [::xmlutils::setXml $cxpath "dv"]
	set cxpath "$cxpathtoDEMresults//i.DEM-ReynoldsNumber"
	set PrintOrNot [::xmlutils::setXml $cxpath "dv"]
	puts $fileid "\"print_REYNOLDS_NUMBER_option\"               : [::wkcf::TranslateToBinaryJson $PrintOrNot],"
	set cxpath "$cxpathtoDEMresults//i.DEM-PressureGradProj"
	set PrintOrNot [::xmlutils::setXml $cxpath "dv"]
	puts $fileid "\"print_PRESSURE_GRAD_PROJECTED_option\"       : [::wkcf::TranslateToBinaryJson $PrintOrNot],"
	set cxpath "$cxpathtoDEMresults//i.DEM-FluidVelProj"
	set PrintOrNot [::xmlutils::setXml $cxpath "dv"]
	puts $fileid "\"print_FLUID_VEL_PROJECTED_option\"           : [::wkcf::TranslateToBinaryJson $PrintOrNot],"
	set cxpath "$cxpathtoDEMresults//i.DEM-FluidAccelProj"
	set PrintOrNot [::xmlutils::setXml $cxpath "dv"]
	puts $fileid "\"print_FLUID_ACCEL_PROJECTED_option\"         : [::wkcf::TranslateToBinaryJson $PrintOrNot],"
	set cxpath "$cxpathtoDEMresults//i.DEM-Buoyancy"
	set PrintOrNot [::xmlutils::setXml $cxpath "dv"]
	puts $fileid "\"print_BUOYANCY_option\"                      : [::wkcf::TranslateToBinaryJson $PrintOrNot],"
	set cxpath "$cxpathtoDEMresults//i.DEM-Drag"
	set PrintOrNot [::xmlutils::setXml $cxpath "dv"]
	puts $fileid "\"print_DRAG_FORCE_option\"                    : [::wkcf::TranslateToBinaryJson $PrintOrNot],"
	set cxpath "$cxpathtoDEMresults//i.DEM-VirtualMass"
	set PrintOrNot [::xmlutils::setXml $cxpath "dv"]
	puts $fileid "\"print_VIRTUAL_MASS_FORCE_option\"            : [::wkcf::TranslateToBinaryJson $PrintOrNot],"
	set cxpath "$cxpathtoDEMresults//i.DEM-Basset"
	set PrintOrNot [::xmlutils::setXml $cxpath "dv"]
	puts $fileid "\"print_BASSET_FORCE_option\"                  : [::wkcf::TranslateToBinaryJson $PrintOrNot],"
	set cxpath "$cxpathtoDEMresults//i.DEM-Lift"
	set PrintOrNot [::xmlutils::setXml $cxpath "dv"]
	puts $fileid "\"print_LIFT_FORCE_option\"                    : [::wkcf::TranslateToBinaryJson $PrintOrNot],"
	set cxpath "$cxpathtoDEMresults//i.DEM-FluidVelRate"
	set PrintOrNot [::xmlutils::setXml $cxpath "dv"]
	puts $fileid "\"print_FLUID_VEL_PROJECTED_RATE_option\"      : [::wkcf::TranslateToBinaryJson $PrintOrNot],"
	set cxpath "$cxpathtoDEMresults//i.DEM-FluidViscosity"
	set PrintOrNot [::xmlutils::setXml $cxpath "dv"]
	puts $fileid "\"print_FLUID_VISCOSITY_PROJECTED_option\"     : [::wkcf::TranslateToBinaryJson $PrintOrNot],"
	set cxpath "$cxpathtoDEMresults//i.DEM-FluidFractionProj"
	set PrintOrNot [::xmlutils::setXml $cxpath "dv"]
	puts $fileid "\"print_FLUID_FRACTION_PROJECTED_option\"      : [::wkcf::TranslateToBinaryJson $PrintOrNot],"
	set cxpath "$cxpathtoDEMresults//i.DEM-VelocityLaplacianProjected"
	set PrintOrNot [::xmlutils::setXml $cxpath "dv"]
	puts $fileid "\"print_FLUID_VEL_LAPL_PROJECTED_option\"      : [::wkcf::TranslateToBinaryJson $PrintOrNot],"
	set cxpath "$cxpathtoDEMresults//i.DEM-VelocityLaplacianRateProjected"
	set PrintOrNot [::xmlutils::setXml $cxpath "dv"]
	puts $fileid "\"print_FLUID_VEL_LAPL_RATE_PROJECTED_option\" : [::wkcf::TranslateToBinaryJson $PrintOrNot],"
	set cxpath "$cxpathtoDEMresults//i.DEM-HydroForce"
	set PrintOrNot [::xmlutils::setXml $cxpath "dv"]
	puts $fileid "\"print_HYDRODYNAMIC_FORCE_option\"            : [::wkcf::TranslateToBinaryJson $PrintOrNot],"
	set cxpath "$cxpathtoDEMresults//i.DEM-HydroMoment"
	set PrintOrNot [::xmlutils::setXml $cxpath "dv"]
	puts $fileid "\"print_HYDRODYNAMIC_MOMENT_option\"           : [::wkcf::TranslateToBinaryJson $PrintOrNot],"

	# Fluid variables
    set cxpath "$cxpathtoFLUIDresults//i.Fluid-AugmentedVel"
	set PrintOrNot [::xmlutils::setXml $cxpath "dv"]
	puts $fileid "\"print_MESH_VELOCITY1_option\"                : [::wkcf::TranslateToBinaryJson $PrintOrNot],"
    set cxpath "$cxpathtoFLUIDresults//i.Fluid-BodyForce"
	set PrintOrNot [::xmlutils::setXml $cxpath "dv"]
	puts $fileid "\"print_BODY_FORCE_option\"                    : [::wkcf::TranslateToBinaryJson $PrintOrNot],"
	set cxpath "$cxpathtoFLUIDresults//i.Fluid-FluidFraction"
	set PrintOrNot [::xmlutils::setXml $cxpath "dv"]
	puts $fileid "\"print_FLUID_FRACTION_option\"                : [::wkcf::TranslateToBinaryJson $PrintOrNot],"
	set cxpath "$cxpathtoFLUIDresults//i.Fluid-FluidFractionGrad"
	set PrintOrNot [::xmlutils::setXml $cxpath "dv"]
	puts $fileid "\"print_FLUID_FRACTION_GRADIENT_option\"       : [::wkcf::TranslateToBinaryJson $PrintOrNot],"
	set cxpath "$cxpathtoFLUIDresults//i.Fluid-HydroReaction"
	set PrintOrNot [::xmlutils::setXml $cxpath "dv"]
	puts $fileid "\"print_HYDRODYNAMIC_REACTION_option\"         : [::wkcf::TranslateToBinaryJson $PrintOrNot],"
	set cxpath "$cxpathtoFLUIDresults//i.Fluid-Pressure"
	set PrintOrNot [::xmlutils::setXml $cxpath "dv"]
	puts $fileid "\"print_PRESSURE_option\"                      : [::wkcf::TranslateToBinaryJson $PrintOrNot],"
	set cxpath "$cxpathtoFLUIDresults//i.Fluid-PressureGrad"
	set PrintOrNot [::xmlutils::setXml $cxpath "dv"]
	puts $fileid "\"print_PRESSURE_GRADIENT_option\"             : [::wkcf::TranslateToBinaryJson $PrintOrNot],"
	set cxpath "$cxpathtoFLUIDresults//i.Fluid-SolidFraction"
	set PrintOrNot [::xmlutils::setXml $cxpath "dv"]
	puts $fileid "\"print_DISPERSE_FRACTION_option\"             : [::wkcf::TranslateToBinaryJson $PrintOrNot],"
	set cxpath "$cxpathtoFLUIDresults//i.Fluid-MeanHydroReaction"
	set PrintOrNot [::xmlutils::setXml $cxpath "dv"]
	puts $fileid "\"print_MEAN_HYDRODYNAMIC_REACTION_option\"    : [::wkcf::TranslateToBinaryJson $PrintOrNot],"
	set cxpath "$cxpathtoFLUIDresults//i.Fluid-VelocityLaplacian"
	set PrintOrNot [::xmlutils::setXml $cxpath "dv"]
	puts $fileid "\"print_VELOCITY_LAPLACIAN_option\"            : [::wkcf::TranslateToBinaryJson $PrintOrNot],"
	set cxpath "$cxpathtoFLUIDresults//i.Fluid-VelocityLaplacianRate"
	set PrintOrNot [::xmlutils::setXml $cxpath "dv"]
	puts $fileid "\"print_VELOCITY_LAPLACIAN_RATE_option\"       : [::wkcf::TranslateToBinaryJson $PrintOrNot],"

	# SWIMMING-SPECIFIC SECTION ENDS ###########################################################################
    }
}

proc ::wkcf::WriteExplicitSolverVariables {} {
    # Write constitutive laws properties
    variable dprops;  variable ActiveAppList
    global KPriv
    variable ndime
    set AppId "DEM"
    set filename "DEM_explicit_solver_var.py"
    set PDir [::KUtils::GetPaths "PDir"]
    set fullname [file native [file join $PDir $filename]]
    # First delete the file
    if {[file exists $fullname]} {
	set res [file delete -force -- $fullname]
    }
    if {[catch {set fileid [open $fullname w+]}]} {
	WarnWin [= "Cannot write file %s. Permission denied" $fullname].
	return 0
    }
    # Write the group properties
    set cproperty "dv"
    set rootid $AppId
    puts $fileid ""
    puts $fileid "# DEM General Options"

    # Dimension
    set cxpath "GeneralApplicationData//c.Domain//i.SpatialDimension"
    set ndim [string range [::xmlutils::setXml $cxpath $cproperty] 0 0]
    puts $fileid "Dimension                        = $ndim"
    #if {$KPriv(what_dempack_package) ne "C-DEMPack"} {}
    # Periodicity of the DEM domain
    set cxpath "$rootid//c.DEM-Options//c.DEM-Physical-opts//i.PeriodicDomainOption"
    set PeriodicDomainOption [::xmlutils::setXml $cxpath $cproperty]

    if {$PeriodicDomainOption == "Yes"} {
	puts $fileid "PeriodicDomainOption             = \"ON\""
    } else {
	puts $fileid "PeriodicDomainOption             = \"OFF\""
    }
    #
    # Get the use bounding box
    set cxpath "$rootid//c.DEM-Options//c.DEM-Boundingbox//i.UseBoundingBox"
    set UseBoundingBox [::xmlutils::setXml $cxpath $cproperty]
    if {$UseBoundingBox == "Yes"} {
	puts $fileid "BoundingBoxOption                = \"ON\""

	# Get the enlargement factor
	set enl_fct_path "$rootid//c.DEM-Options//c.DEM-Boundingbox//i.EnlargementFactor"
	set EnlargementFactor [::xmlutils::setXml $enl_fct_path $cproperty]

	# Get the bounding box type
	set cxpath "$rootid//c.DEM-Options//c.DEM-Boundingbox//i.BoundingBoxType"
	set BoundingBoxType [::xmlutils::setXml $cxpath $cproperty]
	set cxpath "$rootid//c.DEM-Options//c.DEM-Boundingbox//i.StartTime"
	set BoundingBoxStartTime [::xmlutils::setXml $cxpath $cproperty]
	set cxpath "$rootid//c.DEM-Options//c.DEM-Boundingbox//i.StopTime"
	set BoundingBoxStopTime [::xmlutils::setXml $cxpath $cproperty]
	if {$BoundingBoxType eq "Automatic"} {
	    puts $fileid "AutomaticBoundingBoxOption       = \"ON\""
	    puts $fileid "BoundingBoxEnlargementFactor     = $EnlargementFactor"
	    puts $fileid "BoundingBoxStartTime             = $BoundingBoxStartTime"
	    puts $fileid "BoundingBoxStopTime              = $BoundingBoxStopTime"
	    # Write default bounding box values
	    ::wkcf::WriteBoundingBoxDefaults $fileid

	} elseif {$BoundingBoxType eq "Fixed"} {
	    puts $fileid "AutomaticBoundingBoxOption       = \"OFF\""
	    puts $fileid "BoundingBoxEnlargementFactor     = 1.0"
	    puts $fileid "BoundingBoxStartTime             = $BoundingBoxStartTime"
	    puts $fileid "BoundingBoxStopTime              = $BoundingBoxStopTime"

	    # Get the bounding limit
	    set varlist [list MaxX MaxY MaxZ MinX MinY MinZ]
	    foreach varid $varlist {
		set cxpath "$rootid//c.DEM-Options//c.DEM-Boundingbox//i.$varid"
		set $varid [::xmlutils::setXml $cxpath $cproperty]
		puts $fileid "BoundingBox$varid                  = [set $varid]"
	    }
	}
    } else {
	puts $fileid "BoundingBoxOption                = \"OFF\""
	puts $fileid "BoundingBoxEnlargementFactor     = 1.0"
	puts $fileid "AutomaticBoundingBoxOption       = \"OFF\""

	# Write default bounding box values
	::wkcf::WriteBoundingBoxDefaults $fileid
    }
    puts $fileid ""

    #Inlet
    set cxpath "$rootid//c.DEM-Conditions//c.DEM-Inlet"
    set group [::xmlutils::setXmlContainerIds $cxpath]

    if {$group != ""} {
	puts $fileid "dem_inlet_option                 = 1"
    } else {
	puts $fileid "dem_inlet_option                 = 0"
    }

    # Get the gravity
    if {"Fluid" in $ActiveAppList} {
	set cxpathtogravity "GeneralApplicationData//c.SDEM-Physical-opts"
    } else {
	set cxpathtogravity "$rootid//c.DEM-Options//c.DEM-Physical-opts"
    }

    set cxpath "$cxpathtogravity//i.GravityValue"
    set usergravitymod [::xmlutils::setXml $cxpath $cproperty]
    set cxpath "$cxpathtogravity//i.Cx"
    set gravityCx [::xmlutils::setXml $cxpath $cproperty]
    set cxpath "$cxpathtogravity//i.Cy"
    set gravityCy [::xmlutils::setXml $cxpath $cproperty]
    set cxpath "$cxpathtogravity//i.Cz"
    if {$ndime eq "3D"} {
	set gravityCz [::xmlutils::setXml $cxpath $cproperty]
    } else {
	set gravityCz "0.0"
    }

    set gravityvector [list $gravityCx $gravityCy $gravityCz]
    set gravitymod [::MathUtils::VectorModulus $gravityvector]

    if {$gravitymod==0.0} {
	W "The null vector is not a valid direction for the gravity force!"
	W "Computations will proceed anyway...\nPlease make sure that value was intentional..."
    }

    if {$gravitymod !=0.0} {
	set endgravityCx [expr {($gravityCx*$usergravitymod)/$gravitymod}]
	set endgravityCy [expr {($gravityCy*$usergravitymod)/$gravitymod}]
	set endgravityCz [expr {($gravityCz*$usergravitymod)/$gravitymod}]
    } else {
	set endgravityCx $gravityCx
	set endgravityCy $gravityCy
	set endgravityCz $gravityCz
    }
    puts $fileid "GravityX                         = $endgravityCx"
    puts $fileid "GravityY                         = $endgravityCy"
    puts $fileid "GravityZ                         = $endgravityCz"
    puts $fileid ""

    # Compute energies
    set UseComputeEnergies [::xmlutils::setXml "$rootid//c.DEM-Results//c.DEM-Graphs//i.UseComputeEnergies" dv]
    if {$UseComputeEnergies eq "Yes"} {
	puts $fileid "EnergyCalculationOption          = 1"
	# Get the bounding limit
	set varlist [list PotentialEnergyReferencePointX PotentialEnergyReferencePointY PotentialEnergyReferencePointZ]
	foreach varid $varlist {
	    set cxpath "$rootid//c.DEM-Results//c.DEM-Graphs//i.$varid"
	    set $varid [::xmlutils::setXml $cxpath $cproperty]
	    puts $fileid "$varid   = [set $varid]"
	}
    } else {
	puts $fileid "EnergyCalculationOption          = 0"
    }

    # Get the use velocity trap
    set cxpath "$rootid//c.DEM-Results//c.DEM-VelocityTrap//i.UseVelocityTrap"
    set UseVelocityTrap [::xmlutils::setXml $cxpath $cproperty]
    if {$UseVelocityTrap eq "Yes"} {
	puts $fileid "VelocityTrapOption               = 1"

	# Get the bounding limit
	set varlist [list MaxX MaxY MaxZ MinX MinY MinZ]
	foreach varid $varlist {
	    set cxpath "$rootid//c.DEM-Results//c.DEM-VelocityTrap//i.$varid"
	    set $varid [::xmlutils::setXml $cxpath $cproperty]
	    puts $fileid "VelocityTrap$varid                  = [set $varid]"
	}

    } else {
	puts $fileid "VelocityTrapOption               = 0"
    }

    # Calculate Rotations
    set cxpath "$rootid//c.DEM-Options//c.DEM-AdvancedOptions//i.DEM-CalculateRotations"
    set useRotationOption [::xmlutils::setXml $cxpath $cproperty]
    if {$useRotationOption == "Yes"} {
	puts $fileid "RotationOption                   = \"ON\""
    } else {
	puts $fileid "RotationOption                   = \"OFF\""
    }

    # Clean IndentationsOption
    set cxpath "$rootid//c.DEM-Options//c.DEM-AdvancedOptions//i.DEM-CleanInitialIndentations"
    set CleanIndentationsOption [::xmlutils::setXml $cxpath $cproperty]
    if {$KPriv(what_dempack_package) eq "C-DEMPack"} {
	puts $fileid "CleanIndentationsOption          = \"OFF\""
    } elseif {$CleanIndentationsOption == "Yes"} {
	puts $fileid "CleanIndentationsOption          = \"ON\""
    } else {
	puts $fileid "CleanIndentationsOption          = \"OFF\""
    }

    if {$KPriv(what_dempack_package) ne "C-DEMPack"} {
	# RemoveBallsInEmbedded
	set cxpath "$rootid//c.DEM-Options//c.DEM-AdvancedOptions//i.DEM-RemoveBallsInEmbedded"
	set RemoveBallsInEmbedded [::xmlutils::setXml $cxpath $cproperty]
	if {$RemoveBallsInEmbedded == "Yes"} {
	    puts $fileid "RemoveBallsInEmbeddedOption      = 1"
	} else {
	    puts $fileid "RemoveBallsInEmbeddedOption      = 0"
	}
	puts $fileid ""
    }

    # Tangency Tolerance
    set cxpath "$rootid//c.DEM-Options//c.DEM-AdvancedOptions//i.DEM-TangencyTolerance"
    set TangencyToleranceType [::xmlutils::setXml $cxpath $cproperty]
    puts $fileid "DeltaOption                      = \"$TangencyToleranceType\""

    if {$TangencyToleranceType eq "Absolute"} {
	set cxpath "$rootid//c.DEM-Options//c.DEM-AdvancedOptions//i.DEM-TangencyToleranceValueAbsolute"
	set TangenTolerValue [::xmlutils::setXml $cxpath $cproperty]
	puts $fileid "SearchTolerance                  = $TangenTolerValue"
    }

    if {$TangencyToleranceType eq "Coordination_Number"} {
	set cxpath "$rootid//c.DEM-Options//c.DEM-AdvancedOptions//i.DEM-TangencyCoordinationNumber"
	set CoordNumber [::xmlutils::setXml $cxpath $cproperty]
	puts $fileid "CoordinationNumber               = $CoordNumber"
    }

    if {$KPriv(what_dempack_package) eq "C-DEMPack"} {
	# Amplified Search Radius Extension
	set cxpath "$rootid//c.DEM-Options//c.DEM-AdvancedOptions//i.DEM-NeighbourSearchAcceptedGap"
	set AcceptedGap [::xmlutils::setXml $cxpath $cproperty]
	puts $fileid "AmplifiedSearchRadiusExtension   = $AcceptedGap"

	set cxpath "$rootid//c.DEM-Options//c.DEM-AdvancedOptions//i.DEM-MaxAmplificationRatioOfSearchRadius"
	set MaxAmplificationRatioOfSearchRadius [::xmlutils::setXml $cxpath $cproperty]
	puts $fileid "MaxAmplificationRatioOfSearchRadius = $MaxAmplificationRatioOfSearchRadius"

    } else {
	puts $fileid "AmplifiedSearchRadiusExtension   = 0.0"
    }

    # Export Model Data Info
    set cxpath "$rootid//c.DEM-Options//c.DEM-AdvancedOptions//i.DEM-ExportModelDataInformation"
    set modatainfo [::xmlutils::setXml $cxpath $cproperty]
    if {$modatainfo eq "Yes"} {
	puts $fileid "ModelDataInfo                    = \"ON\""
    } else {
	puts $fileid "ModelDataInfo                    = \"OFF\""
    }

    # Virtual Mass Coefficient
    set cxpath "$rootid//c.DEM-Options//c.DEM-AdvancedOptions//i.DEM-VirtualMassCoef"
    set vm [::xmlutils::setXml $cxpath $cproperty]
    puts $fileid "VirtualMassCoefficient           = $vm"

    # Rolling Friction
    set cxpath "$rootid//c.DEM-Options//c.DEM-AdvancedOptions//i.DEM-RollingFriction"
    set rf [::xmlutils::setXml $cxpath $cproperty]
    if {$rf eq "Yes"} {
	puts $fileid "RollingFrictionOption            = \"ON\""
    } else {
	puts $fileid "RollingFrictionOption            = \"OFF\""
    }


    # Compute Stress Tensor
    if {$KPriv(what_dempack_package) eq "C-DEMPack"} {
	set tensor [::xmlutils::setXml "$rootid//c.DEM-Options//c.DEM-AdvancedOptions//i.DEM-ComputeStressTensorOption" "dv"]
	if {$tensor eq "Yes"} {
	    puts $fileid "ComputeStressTensorOption        = \"ON\""
	} else {
	    puts $fileid "ComputeStressTensorOption        = \"OFF\""
	}
    }

    # Poisson Effect Option
    if {$KPriv(what_dempack_package) eq "C-DEMPack"} {
	set pe [::xmlutils::setXml "$rootid//c.DEM-Options//c.DEM-AdvancedOptions//i.DEM-PoissonEffect" "dv"]
	if {$pe eq "Yes"} {
	    puts $fileid "PoissonEffectOption              = \"ON\""
	} else {
	    puts $fileid "PoissonEffectOption              = \"OFF\""
	}
    }

    # Shear Strains Parallel To Bonds Effect Option
    if {$KPriv(what_dempack_package) eq "C-DEMPack"} {
	set ssptb [::xmlutils::setXml "$rootid//c.DEM-Options//c.DEM-AdvancedOptions//i.DEM-ShearStrainParallelToBondEffect" "dv"]
	if {$ssptb eq "Yes"} {
	    puts $fileid "ShearStrainParallelToBondOption  = \"ON\""
	} else {
	    puts $fileid "ShearStrainParallelToBondOption  = \"OFF\""
	}
    }

    #set cxpath "$rootid//c.DEM-Options//c.DEM-AdvancedOptions//i.DEM-ComputeEnergies"
    #set ceo [::xmlutils::setXml $cxpath $cproperty]
    #if {$ceo eq "Activate"} {
	#    set ceo 1
	#} else {
	#    set ceo 0
	#}
    #puts $fileid "ComputeEnergiesOption            = $ceo"

    # Dont search until failure
    if {$KPriv(what_dempack_package) eq "C-DEMPack"} {
	set io [::xmlutils::setXml "$rootid//c.DEM-Options//c.DEM-AdvancedOptions//i.DEM-DontSearchUntilFail" "dv"]
	if {$io eq "Yes"} {
	    puts $fileid {DontSearchUntilFailure           = "ON"}
	} else {
	    puts $fileid {DontSearchUntilFailure           = "OFF"}
	}
    }

    #set cxpath "$rootid//c.DEM-Options//c.DEM-AdvancedOptions//i.DEM-DontSearchUntilFail"
    #set DSUF [::xmlutils::setXml $cxpath $cproperty]
    #puts $fileid "DontSearchUntilFailure           = \"$DSUF\""

    # Contact Mesh Option
    if {$KPriv(what_dempack_package) eq "C-DEMPack"} {
	set CMO [::xmlutils::setXml "$rootid//c.DEM-Results//i.DEM-ShowBondElements" dv]
	if {$CMO eq "Yes"} {
	    puts $fileid "ContactMeshOption                = \"ON\""
	} else {
	    puts $fileid "ContactMeshOption                = \"OFF\""
	}
    } else {
	puts $fileid "ContactMeshOption                = \"OFF\""
    }

    # Output file Type
    if {"Fluid" in $ActiveAppList} {
	set cxpath "GeneralApplicationData//c.SDEM-Results//c.GiDOptions//i.GiDPostMode"
    } else {
	set cxpath "$rootid//c.DEM-Results//c.GiDOptions//i.GiDPostMode"
    }
    set GiDMode [::xmlutils::setXml $cxpath $cproperty]
    puts $fileid "OutputFileType                   = \"$GiDMode\""

    # Multifile
    if {"Fluid" in $ActiveAppList} {
	set cxpath "GeneralApplicationData//c.SDEM-Results//c.GiDOptions//i.GiDMultiFileFlag"
    } else {
	set cxpath "$rootid//c.DEM-Results//c.GiDOptions//i.GiDMultiFileFlag"
    }
    set GiDMultiFileFlag [::xmlutils::setXml $cxpath $cproperty]
    if {$GiDMultiFileFlag eq "Single"} {
	puts $fileid "Multifile                        = \"single_file\""
    } else {
	puts $fileid "Multifile                        = \"multiple_files\""
    }

    puts $fileid "ElementType                      = \"[::wkcf::GetElementType]\""
    ######################################################################################

    puts $fileid ""
    puts $fileid "# Solution Strategy"

    # Translational Integration Scheme
    set cxpath "$rootid//c.DEM-SolutionStrategy//i.DEM-TimeTranslationalIntegrationScheme"
    set TransIntegScheme [::xmlutils::setXml $cxpath $cproperty]

    puts $fileid "TranslationalIntegrationScheme       = \"$TransIntegScheme\""

    # Rotational Integration Scheme
    set cxpath "$rootid//c.DEM-SolutionStrategy//i.DEM-TimeRotationalIntegrationScheme"
    set RotIntegScheme [::xmlutils::setXml $cxpath $cproperty]

    puts $fileid "RotationalIntegrationScheme          = \"$RotIntegScheme\""

    # Auto Reduction of Time Step Option
    set cxpath "$rootid//c.DEM-SolutionStrategy//c.DEM-TimeParameters//i.UseAutomaticDeltaTime"
    set DTimeSelection [::xmlutils::setXml $cxpath $cproperty]

    if {$DTimeSelection eq "Fixed"} {
	puts $fileid "AutomaticTimestep                = \"OFF\""
    } else {
	puts $fileid "AutomaticTimestep                = \"ON\""
    }
    set cxpath "$rootid//c.DEM-SolutionStrategy//c.DEM-TimeParameters//i.DeltaTimeSafetyFactor"
    set SafetyFact [::xmlutils::setXml $cxpath $cproperty]
    puts $fileid "DeltaTimeSafetyFactor            = $SafetyFact"

    set cxpath "$rootid//c.DEM-SolutionStrategy//c.DEM-TimeParameters//i.DeltaTime"
    set MaxTimeStep [::xmlutils::setXml $cxpath $cproperty]
    puts $fileid "MaxTimeStep                      = $MaxTimeStep"

    # Final Time
    if {"Fluid" in $ActiveAppList} {
	set cxpath "GeneralApplicationData//c.SimulationOptions//i.EndTime"
	set TTime [::xmlutils::setXml $cxpath $cproperty]
    } else {
	set cxpath "$rootid//c.DEM-SolutionStrategy//c.DEM-TimeParameters//i.DEM-TotalTime"
	set TTime [::xmlutils::setXml $cxpath $cproperty]
    }
    puts $fileid "FinalTime                        = $TTime"

    # Control Time
    set cxpath "$rootid//c.DEM-SolutionStrategy//c.DEM-TimeParameters//i.DEM-ScreenInfoOutput"
    set CTime [::xmlutils::setXml $cxpath $cproperty]
    puts $fileid "ControlTime                      = $CTime"

    # Neighbour Search Frequency
    set cxpath "$rootid//c.DEM-SolutionStrategy//c.DEM-TimeParameters//i.DEM-NeighbourSearchFrequency"
    set FrecTime [::xmlutils::setXml $cxpath $cproperty]
    puts $fileid "NeighbourSearchFrequency         = $FrecTime"

    # Material Test Data #########################################################################################
    if {$KPriv(what_dempack_package) eq "C-DEMPack"} {
		set cproperty "dv"
		set cxpath "DEM//c.DEM-Options//c.DEM-AdvancedOptions//i.DEM-MaterialTestActivate"
		set material_option [::xmlutils::setXml $cxpath $cproperty]
		if {$material_option eq "Yes"} {
		        ::wkcf::WriteMatTestData $fileid
		}
    }

    # SWIMMING-SPECIFIC SECTION BEGINS ###########################################################################
    if {"Fluid" in $ActiveAppList} {
	set cxpath "$rootid//c.DEM-Fluid-interaction//i.DEM-TwoWayCoupling"
	set cproperty "dv"
	set cxpath "GeneralApplicationData//c.CouplingParameters//i.CouplingLevel"
	set coupling_level_type [::xmlutils::setXml $cxpath $cproperty]
    set cproperty "dv"
    set cxpath "GeneralApplicationData//c.CouplingParameters//i.TimeAveragingType"
    set time_averaging_type [::xmlutils::setXml $cxpath $cproperty]

    set interaction_start_time [::xmlutils::setXml {GeneralApplicationData//c.CouplingParameters//i.InteractionStartTime} dv]

	set cxpath "GeneralApplicationData//c.CouplingParameters//c.HydrodynamicForceModel//i.PickIndividualForces"
	set pick_individual_forces_option [::wkcf::TranslateToBinary [::xmlutils::setXml $cxpath $cproperty]]

	set cxpath "GeneralApplicationData//c.CouplingParameters//c.HydrodynamicForceModel//i.BuoyancyForceType"
	set buoyancy_force_type [::xmlutils::setXml $cxpath $cproperty]

	set cxpath "GeneralApplicationData//c.CouplingParameters//c.HydrodynamicForceModel//i.DragForceType"
	set drag_force_type [::xmlutils::setXml $cxpath $cproperty]

	set cxpath "GeneralApplicationData//c.CouplingParameters//c.HydrodynamicForceModel//i.VirtualMassForceType"
	set virtual_mass_force_type [::xmlutils::setXml $cxpath $cproperty]

	set cxpath "GeneralApplicationData//c.CouplingParameters//c.HydrodynamicForceModel//i.SaffmanLiftForceType"
	set lift_force_type [::xmlutils::setXml $cxpath $cproperty]

	set cxpath "GeneralApplicationData//c.CouplingParameters//c.HydrodynamicForceModel//i.MagnusLiftForceType"
	set magnus_force_type [::xmlutils::setXml $cxpath $cproperty]

	set cxpath "GeneralApplicationData//c.CouplingParameters//c.HydrodynamicForceModel//i.HydrodynamicTorqueType"
	set hydro_torque_type [::xmlutils::setXml $cxpath $cproperty]

	set cxpath "GeneralApplicationData//c.CouplingParameters//c.HydrodynamicForceModel//i.DragModifier"
	set drag_modifier_type [::xmlutils::setXml $cxpath $cproperty]

	set cxpath "GeneralApplicationData//c.CouplingParameters//c.HydrodynamicForceModel//i.ViscosityModificationType"
	set viscosity_modification_type [::xmlutils::setXml $cxpath $cproperty]

	set cxpath "GeneralApplicationData//c.CouplingParameters//c.HydrodynamicForceModel//c.MaxeyRileyEquations//i.MRDragModifier"
	set MR_drag_modifier_type [::xmlutils::setXml $cxpath $cproperty]

	set cxpath "GeneralApplicationData//c.PostProcessingVariables//i.GradientCalculationType"
	set gradient_calculation_type [::xmlutils::setXml $cxpath $cproperty]

	set cxpath "GeneralApplicationData//c.PostProcessingVariables//i.VelocityLaplacianCalculationType"
	set velocity_laplacian_calculation_type [::xmlutils::setXml $cxpath $cproperty]

	if {$velocity_laplacian_calculation_type == "0"} {
	    set cxpath "GeneralApplicationData//c.CouplingParameters//c.HydrodynamicForceModel//i.IncludeFaxenTerms0"
	} elseif {$velocity_laplacian_calculation_type == "1"} {
	    set cxpath "GeneralApplicationData//c.CouplingParameters//c.HydrodynamicForceModel//i.IncludeFaxenTerms1"
	} elseif {$velocity_laplacian_calculation_type == "2"} {
	    set cxpath "GeneralApplicationData//c.CouplingParameters//c.HydrodynamicForceModel//i.IncludeFaxenTerms"
	}

	set include_faxen_terms_option [::wkcf::TranslateToBinary [::xmlutils::setXml $cxpath $cproperty]]

	if {!$pick_individual_forces_option} { # Maxey-Riley
	    if {!$include_faxen_terms_option} {
		set drag_force_type 2
		set virtual_mass_force_type 0
	    } else {
		set drag_force_type 11
		set virtual_mass_force_type 11
	    }
	    set lift_force_type 0
	    set magnus_force_type 0
	    set hydro_torque_type 0
	    set drag_modifier_type $MR_drag_modifier_type

		set buoyancy_force_type 2
	}

	puts $fileid ""
	puts $fileid "# Swimming DEM-specific section begins"
	puts $fileid "#-------------------------------------"
	puts $fileid "coupling_level_type                    = $coupling_level_type"
    puts $fileid "time_averaging_type                    = $time_averaging_type"
    puts $fileid "interaction_start_time                 = $interaction_start_time"
    puts $fileid "pick_individual_forces_option          = $pick_individual_forces_option"
    puts $fileid "include_faxen_terms_option             = $include_faxen_terms_option # (relevant if the Maxey Riley equation is used)"
    puts $fileid "gradient_calculation_type              = $gradient_calculation_type # (Not calculated (0), volume-weighed average(1), Superconvergent recovery(2))"
    puts $fileid "laplacian_calculation_type             = $velocity_laplacian_calculation_type # (Not calculated (0), Finite element projection (1), Superconvergent recovery(2))"
    puts $fileid "buoyancy_force_type                    = $buoyancy_force_type # null buoyancy (0), compute buoyancy (1)  if drag_force_type is 2 buoyancy is always parallel to gravity"
    puts $fileid "drag_force_type                        = $drag_force_type # null drag (0), Stokes (1), Weatherford (2), Ganser (3), Ishii (4), Newtonian Regime (5)"
    puts $fileid "virtual_mass_force_type                = $virtual_mass_force_type # null virtual mass force (0)"
    puts $fileid "lift_force_type                        = $lift_force_type # null lift force (0), Saffman (1)"
    puts $fileid "magnus_force_type                      = $magnus_force_type # null magnus force (0), Rubinow and Keller (1), Oesterle and Bui Dihn (2)"
    puts $fileid "hydro_torque_type                      = $hydro_torque_type # null hydrodynamic torque (0), Dennis (1)"
    puts $fileid "drag_modifier_type                     = $drag_modifier_type"
    puts $fileid "viscosity_modification_type            = $viscosity_modification_type"
    puts $fileid ""
    puts $fileid "# Parameters not yet settable from interface"
    puts $fileid ""
	puts $fileid "coupling_weighing_type                 = 2 # {fluid_to_DEM, DEM_to_fluid, fluid_fraction} = {lin, lin, imposed} (-1), {lin, const, const} (0), {lin, lin, const} (1), {lin, lin, lin} (2), averaging method (3)"
	puts $fileid "fluid_model_type                       = 1 # untouched, velocity incremented by 1/fluid_fraction (0), modified mass conservation only (1)"
	puts $fileid "coupling_scheme_type                   = \"UpdatedFluid\" # \"UpdatedFluid\", \"UpdatedDEM\""
	puts $fileid "print_particles_results_option         = 0"
	puts $fileid "add_each_hydro_force_option            = 1 # add each of the hydrodynamic forces (drag, lift and virtual mass)"
	puts $fileid "project_at_every_substep_option        = 1"
	puts $fileid "velocity_trap_option                   = 0"
	puts $fileid "inlet_option                           = 1"
	puts $fileid "manually_imposed_drag_law_option       = 0"
	puts $fileid "stationary_problem_option              = 0 # stationary, stop calculating the fluid after it reaches the stationary state"
	puts $fileid "flow_in_porous_medium_option           = 0 # the porosity is an imposed field"
	puts $fileid "flow_in_porous_DEM_medium_option       = 0 # the DEM part is kept static"
	puts $fileid "embedded_option                        = 1 # the embedded domain tools are to be used"
	puts $fileid "make_results_directories_option        = 1 # results are written into a folder (../results) inside the problem folder"
	puts $fileid "body_force_on_fluid_option             = 1"
	puts $fileid "print_debug_info_option                = 0 # print a summary of global physical measures"
	puts $fileid "print_particles_results_cycle          = 1 # number of 'ticks' per printing cycle"
	puts $fileid "debug_tool_cycle                       = 10 # number of 'ticks' per debug computations cycle"
	puts $fileid "similarity_transformation_type         = 0 # no transformation (0), Tsuji (1)"
	puts $fileid "dem_inlet_element_type                 = \"SphericSwimmingParticle3D\"  # \"SphericParticle3D\", \"SphericSwimmingParticle3D\""
	puts $fileid "drag_modifier_type                     = 2 # Hayder (2), Chien (3) # problemtype option"
	puts $fileid "drag_porosity_correction_type          = 0 # No correction (0), Richardson and Zaki (1)"
	puts $fileid "min_fluid_fraction                     = 0.2"
	puts $fileid "initial_drag_force                     = 0.0   # problemtype option"
	puts $fileid "drag_law_slope                         = 0.0   # problemtype option"
	puts $fileid "power_law_tol                          = 0.0"
	puts $fileid "model_over_real_diameter_factor        = 1.0 # not active if similarity_transformation_type = 0"
	puts $fileid "max_pressure_variation_rate_tol        = 1e-3 # for stationary problems, criterion to stop the fluid calculations"
	puts $fileid "time_steps_per_stationarity_step       = 15 # number of fluid time steps between consecutive assessment of stationarity steps"
	puts $fileid "meso_scale_length                      = 0.2 # the radius of the support of the averaging function for homogenization (<=0 for automatic calculation)"
	puts $fileid "shape_factor                           = 0.5 # the density function's maximum over its support's radius (only relevant if coupling_weighing_type == 3)"
	puts $fileid "#-------------------------------------"
	puts $fileid "# Swimming DEM-specific section ends"
	puts $fileid ""
    }
    # SWIMMING-SPECIFIC SECTION ENDS ###########################################################################

    puts $fileid ""
    ::wkcf::WritePostProcessData $fileid
    puts $fileid ""
    # Get the problem name
    set PName [::KUtils::GetPaths "PName"]
    puts $fileid "#"
    puts $fileid "problem_name=\"${PName}\""

    # Get the kratos path
    set KratosPath [::xmlutils::setXml {GeneralApplicationData//c.ProjectConfiguration//i.KratosPath} dv]
    set KratosPath [file native $KratosPath]

    close $fileid
}

proc ::wkcf::WriteExplicitSolverVariablesInJsonFile {} {
    # Write constitutive laws properties
    variable dprops;  variable ActiveAppList
    global KPriv
    variable ndime
    set AppId "DEM"
    set filename "ProjectParametersDEM.json"
    set PDir [::KUtils::GetPaths "PDir"]
    set fullname [file native [file join $PDir $filename]]
    # First delete the file
    if {[file exists $fullname]} {
	set res [file delete -force -- $fullname]
    }
    if {[catch {set fileid [open $fullname w+]}]} {
	WarnWin [= "Cannot write file %s. Permission denied" $fullname].
	return 0
    }
    # Write the group properties
    set cproperty "dv"
    set rootid $AppId
    puts $fileid "{"

    # Dimension
    set cxpath "GeneralApplicationData//c.Domain//i.SpatialDimension"
    set ndim [string range [::xmlutils::setXml $cxpath $cproperty] 0 0]
    puts $fileid "\"Dimension\"                        : $ndim,"
    #if {$KPriv(what_dempack_package) ne "C-DEMPack"} {}
    # Periodicity of the DEM domain
    set cxpath "$rootid//c.DEM-Options//c.DEM-Physical-opts//i.PeriodicDomainOption"
    set PeriodicDomainOption [::xmlutils::setXml $cxpath $cproperty]

    if {$PeriodicDomainOption == "Yes"} {
	puts $fileid "\"PeriodicDomainOption\"             : true,"
    } else {
	puts $fileid "\"PeriodicDomainOption\"             : false,"
    }
    #
    # Get the use bounding box
    set cxpath "$rootid//c.DEM-Options//c.DEM-Boundingbox//i.UseBoundingBox"
    set UseBoundingBox [::xmlutils::setXml $cxpath $cproperty]
    if {$UseBoundingBox == "Yes"} {
	puts $fileid "\"BoundingBoxOption\"                : true,"

	# Get the enlargement factor
	set enl_fct_path "$rootid//c.DEM-Options//c.DEM-Boundingbox//i.EnlargementFactor"
	set EnlargementFactor [::xmlutils::setXml $enl_fct_path $cproperty]

	# Get the bounding box type
	set cxpath "$rootid//c.DEM-Options//c.DEM-Boundingbox//i.BoundingBoxType"
	set BoundingBoxType [::xmlutils::setXml $cxpath $cproperty]
	set cxpath "$rootid//c.DEM-Options//c.DEM-Boundingbox//i.StartTime"
	set BoundingBoxStartTime [::xmlutils::setXml $cxpath $cproperty]
	set cxpath "$rootid//c.DEM-Options//c.DEM-Boundingbox//i.StopTime"
	set BoundingBoxStopTime [::xmlutils::setXml $cxpath $cproperty]
	if {$BoundingBoxType eq "Automatic"} {
	    puts $fileid "\"AutomaticBoundingBoxOption\"       : true,"
	    puts $fileid "\"BoundingBoxEnlargementFactor\"     : $EnlargementFactor,"
	    puts $fileid "\"BoundingBoxStartTime\"             : $BoundingBoxStartTime,"
	    puts $fileid "\"BoundingBoxStopTime\"              : $BoundingBoxStopTime,"
	    # Write default bounding box values
	    ::wkcf::WriteBoundingBoxDefaultsInJsonFile $fileid

	} elseif {$BoundingBoxType eq "Fixed"} {
	    puts $fileid "\"AutomaticBoundingBoxOption\"       : false,"
	    puts $fileid "\"BoundingBoxEnlargementFactor\"     : 1.0,"
	    puts $fileid "\"BoundingBoxStartTime\"             : $BoundingBoxStartTime,"
	    puts $fileid "\"BoundingBoxStopTime\"              : $BoundingBoxStopTime,"

	    # Get the bounding limit
	    set varlist [list MaxX MaxY MaxZ MinX MinY MinZ]
	    foreach varid $varlist {
		set cxpath "$rootid//c.DEM-Options//c.DEM-Boundingbox//i.$varid"
		set $varid [::xmlutils::setXml $cxpath $cproperty]
		puts $fileid "\"BoundingBox$varid\"                  : [set $varid],"
	    }
	}
    } else {
	puts $fileid "\"BoundingBoxOption\"                : false,"
	puts $fileid "\"BoundingBoxEnlargementFactor\"     : 1.0,"
	puts $fileid "\"AutomaticBoundingBoxOption\"       : false,"

	# Write default bounding box values
	::wkcf::WriteBoundingBoxDefaultsInJsonFile $fileid
    }
    puts $fileid ""

    #Inlet
    set cxpath "$rootid//c.DEM-Conditions//c.DEM-Inlet"
    set group [::xmlutils::setXmlContainerIds $cxpath]

    if {$group != ""} {
	puts $fileid "\"dem_inlet_option\"                 : true,"
    } else {
	puts $fileid "\"dem_inlet_option\"                 : false,"
    }

    # Get the gravity
    if {"Fluid" in $ActiveAppList} {
	set cxpathtogravity "GeneralApplicationData//c.SDEM-Physical-opts"
    } else {
	set cxpathtogravity "$rootid//c.DEM-Options//c.DEM-Physical-opts"
    }

    set cxpath "$cxpathtogravity//i.GravityValue"
    set usergravitymod [::xmlutils::setXml $cxpath $cproperty]
    set cxpath "$cxpathtogravity//i.Cx"
    set gravityCx [::xmlutils::setXml $cxpath $cproperty]
    set cxpath "$cxpathtogravity//i.Cy"
    set gravityCy [::xmlutils::setXml $cxpath $cproperty]
    set cxpath "$cxpathtogravity//i.Cz"
    if {$ndime eq "3D"} {
	set gravityCz [::xmlutils::setXml $cxpath $cproperty]
    } else {
	set gravityCz "0.0"
    }

    set gravityvector [list $gravityCx $gravityCy $gravityCz]
    set gravitymod [::MathUtils::VectorModulus $gravityvector]

    if {$gravitymod==0.0} {
	W "The null vector is not a valid direction for the gravity force!"
	W "Computations will proceed anyway...\nPlease make sure that value was intentional..."
    }

    if {$gravitymod !=0.0} {
	set endgravityCx [expr {($gravityCx*$usergravitymod)/$gravitymod}]
	set endgravityCy [expr {($gravityCy*$usergravitymod)/$gravitymod}]
	set endgravityCz [expr {($gravityCz*$usergravitymod)/$gravitymod}]
    } else {
	set endgravityCx $gravityCx
	set endgravityCy $gravityCy
	set endgravityCz $gravityCz
    }
    puts $fileid "\"GravityX\"                         : $endgravityCx,"
    puts $fileid "\"GravityY\"                         : $endgravityCy,"
    puts $fileid "\"GravityZ\"                         : $endgravityCz,"
    puts $fileid ""

    # Compute energies
    set UseComputeEnergies [::xmlutils::setXml "$rootid//c.DEM-Results//c.DEM-Graphs//i.UseComputeEnergies" dv]
    if {$UseComputeEnergies eq "Yes"} {
	puts $fileid "\"EnergyCalculationOption\"          : true,"
	# Get the bounding limit
	set varlist [list PotentialEnergyReferencePointX PotentialEnergyReferencePointY PotentialEnergyReferencePointZ]
	foreach varid $varlist {
	    set cxpath "$rootid//c.DEM-Results//c.DEM-Graphs//i.$varid"
	    set $varid [::xmlutils::setXml $cxpath $cproperty]
	    puts $fileid "\"$varid\"   : [set $varid],"
	}
    } else {
	puts $fileid "\"EnergyCalculationOption\"          : false,"
    }

    # Get the use velocity trap
    set cxpath "$rootid//c.DEM-Results//c.DEM-VelocityTrap//i.UseVelocityTrap"
    set UseVelocityTrap [::xmlutils::setXml $cxpath $cproperty]
    if {$UseVelocityTrap eq "Yes"} {
	puts $fileid "\"VelocityTrapOption\"               : true,"

	# Get the bounding limit
	set varlist [list MaxX MaxY MaxZ MinX MinY MinZ]
	foreach varid $varlist {
	    set cxpath "$rootid//c.DEM-Results//c.DEM-VelocityTrap//i.$varid"
	    set $varid [::xmlutils::setXml $cxpath $cproperty]
	    puts $fileid "\"VelocityTrap$varid\"                  : [set $varid],"
	}

    } else {
	puts $fileid "\"VelocityTrapOption\"               : false,"
    }

    # Calculate Rotations
    set cxpath "$rootid//c.DEM-Options//c.DEM-AdvancedOptions//i.DEM-CalculateRotations"
    set useRotationOption [::xmlutils::setXml $cxpath $cproperty]
    if {$useRotationOption == "Yes"} {
	puts $fileid "\"RotationOption\"                   : true,"
    } else {
	puts $fileid "\"RotationOption\"                   : false,"
    }

    # Clean IndentationsOption
    set cxpath "$rootid//c.DEM-Options//c.DEM-AdvancedOptions//i.DEM-CleanInitialIndentations"
    set CleanIndentationsOption [::xmlutils::setXml $cxpath $cproperty]
    if {$KPriv(what_dempack_package) eq "C-DEMPack"} {
	puts $fileid "\"CleanIndentationsOption\"          : false,"
    } elseif {$CleanIndentationsOption == "Yes"} {
	puts $fileid "\"CleanIndentationsOption\"          : true,"
    } else {
	puts $fileid "\"CleanIndentationsOption\"          : false,"
    }

	puts $fileid "\"strategy_parameters\" : {"

	# Remove initially indented balls with walls
	set cxpath "$rootid//c.DEM-Options//c.DEM-AdvancedOptions//i.DEM-RemoveBallsInitiallyTouchingWalls"
    set RemoveBallsInitiallyTouchingWallsOption [::xmlutils::setXml $cxpath $cproperty]
  	if {$RemoveBallsInitiallyTouchingWallsOption == "Yes"} {
		puts $fileid "    \"RemoveBallsInitiallyTouchingWalls\"          : true"
    } else {
		puts $fileid "    \"RemoveBallsInitiallyTouchingWalls\"          : false"
    }

	puts $fileid "},"

    if {$KPriv(what_dempack_package) ne "C-DEMPack"} {
	# RemoveBallsInEmbedded
	set cxpath "$rootid//c.DEM-Options//c.DEM-AdvancedOptions//i.DEM-RemoveBallsInEmbedded"
	set RemoveBallsInEmbedded [::xmlutils::setXml $cxpath $cproperty]
	if {$RemoveBallsInEmbedded == "Yes"} {
	    puts $fileid "\"RemoveBallsInEmbeddedOption\"      : true,"
	} else {
	    puts $fileid "\"RemoveBallsInEmbeddedOption\"      : false,"
	}
	puts $fileid ""
    }

    # Tangency Tolerance
    set cxpath "$rootid//c.DEM-Options//c.DEM-AdvancedOptions//i.DEM-TangencyTolerance"
    set TangencyToleranceType [::xmlutils::setXml $cxpath $cproperty]
    puts $fileid "\"DeltaOption\"                      : \"$TangencyToleranceType\","

    if {$TangencyToleranceType eq "Absolute"} {
	set cxpath "$rootid//c.DEM-Options//c.DEM-AdvancedOptions//i.DEM-TangencyToleranceValueAbsolute"
	set TangenTolerValue [::xmlutils::setXml $cxpath $cproperty]
	puts $fileid "\"SearchTolerance\"                  : $TangenTolerValue,"
    }

    if {$TangencyToleranceType eq "Coordination_Number"} {
	set cxpath "$rootid//c.DEM-Options//c.DEM-AdvancedOptions//i.DEM-TangencyCoordinationNumber"
	set CoordNumber [::xmlutils::setXml $cxpath $cproperty]
	puts $fileid "\"CoordinationNumber\"               : $CoordNumber,"
    }

    if {$KPriv(what_dempack_package) eq "C-DEMPack"} {
	# Amplified Search Radius Extension
	set cxpath "$rootid//c.DEM-Options//c.DEM-AdvancedOptions//i.DEM-NeighbourSearchAcceptedGap"
	set AcceptedGap [::xmlutils::setXml $cxpath $cproperty]
	puts $fileid "\"AmplifiedSearchRadiusExtension\"   : $AcceptedGap,"

	set cxpath "$rootid//c.DEM-Options//c.DEM-AdvancedOptions//i.DEM-MaxAmplificationRatioOfSearchRadius"
	set MaxAmplificationRatioOfSearchRadius [::xmlutils::setXml $cxpath $cproperty]
	puts $fileid "\"MaxAmplificationRatioOfSearchRadius\" : $MaxAmplificationRatioOfSearchRadius,"

    } else {
	puts $fileid "\"AmplifiedSearchRadiusExtension\"   : 0.0,"
    }

    # Export Model Data Info
    set cxpath "$rootid//c.DEM-Options//c.DEM-AdvancedOptions//i.DEM-ExportModelDataInformation"
    set modatainfo [::xmlutils::setXml $cxpath $cproperty]
    if {$modatainfo eq "Yes"} {
	puts $fileid "\"ModelDataInfo\"                    : true,"
    } else {
	puts $fileid "\"ModelDataInfo\"                    : false,"
    }

    # Virtual Mass Coefficient
    set cxpath "$rootid//c.DEM-Options//c.DEM-AdvancedOptions//i.DEM-VirtualMassCoef"
    set vm [::xmlutils::setXml $cxpath $cproperty]
    puts $fileid "\"VirtualMassCoefficient\"           : $vm,"

    # Rolling Friction
    set cxpath "$rootid//c.DEM-Options//c.DEM-AdvancedOptions//i.DEM-RollingFriction"
    set rf [::xmlutils::setXml $cxpath $cproperty]
    if {$rf eq "Yes"} {
	puts $fileid "\"RollingFrictionOption\"            : true,"
    } else {
	puts $fileid "\"RollingFrictionOption\"            : false,"
    }

    # Global Damping
    puts $fileid "\"GlobalDamping\"                    : [::xmlutils::setXml "$rootid//c.DEM-Options//c.DEM-AdvancedOptions//i.DEM-GlobalDamping" dv],"

    # Compute Stress Tensor
    if {$KPriv(what_dempack_package) eq "C-DEMPack"} {
	set tensor [::xmlutils::setXml "$rootid//c.DEM-Options//c.DEM-AdvancedOptions//i.DEM-ComputeStressTensorOption" "dv"]
	if {$tensor eq "Yes"} {
	    puts $fileid "\"ComputeStressTensorOption\"        : true,"
	} else {
	    puts $fileid "\"ComputeStressTensorOption\"        : false,"
	}
    }

    # Poisson Effect Option
    if {$KPriv(what_dempack_package) eq "C-DEMPack"} {
	set pe [::xmlutils::setXml "$rootid//c.DEM-Options//c.DEM-AdvancedOptions//i.DEM-PoissonEffect" "dv"]
	if {$pe eq "Yes"} {
	    puts $fileid "\"PoissonEffectOption\"              : true,"
	} else {
	    puts $fileid "\"PoissonEffectOption\"              : false,"
	}
    }

    # Shear Strains Parallel To Bonds Effect Option
    if {$KPriv(what_dempack_package) eq "C-DEMPack"} {
	set ssptb [::xmlutils::setXml "$rootid//c.DEM-Options//c.DEM-AdvancedOptions//i.DEM-ShearStrainParallelToBondEffect" "dv"]
	if {$ssptb eq "Yes"} {
	    puts $fileid "\"ShearStrainParallelToBondOption\"  : true,"
	} else {
	    puts $fileid "\"ShearStrainParallelToBondOption\"  : false,"
	}
    }

    #set cxpath "$rootid//c.DEM-Options//c.DEM-AdvancedOptions//i.DEM-ComputeEnergies"
    #set ceo [::xmlutils::setXml $cxpath $cproperty]
    #if {$ceo eq "Activate"} {
	#    set ceo 1
	#} else {
	#    set ceo 0
	#}
    #puts $fileid "\"ComputeEnergiesOption\"            : $ceo"

    # Dont search until failure
    if {$KPriv(what_dempack_package) eq "C-DEMPack"} {
	set io [::xmlutils::setXml "$rootid//c.DEM-Options//c.DEM-AdvancedOptions//i.DEM-DontSearchUntilFail" "dv"]
	if {$io eq "Yes"} {
	    puts $fileid "\"DontSearchUntilFailure\"           : true,"
	} else {
	    puts $fileid "\"DontSearchUntilFailure\"           : false,"
	}
    }

    #set cxpath "$rootid//c.DEM-Options//c.DEM-AdvancedOptions//i.DEM-DontSearchUntilFail"
    #set DSUF [::xmlutils::setXml $cxpath $cproperty]
    #puts $fileid "\"DontSearchUntilFailure\"           : \"$DSUF\""

    # Contact Mesh Option
    if {$KPriv(what_dempack_package) eq "C-DEMPack"} {
	set CMO [::xmlutils::setXml "$rootid//c.DEM-Results//i.DEM-ShowBondElements" dv]
	if {$CMO eq "Yes"} {
	    puts $fileid "\"ContactMeshOption\"                : true,"
	} else {
	    puts $fileid "\"ContactMeshOption\"                : false,"
	}
    } else {
	puts $fileid "\"ContactMeshOption\"                : false,"
    }

    # Output file Type
    if {"Fluid" in $ActiveAppList} {
	set cxpath "GeneralApplicationData//c.SDEM-Results//c.GiDOptions//i.GiDPostMode"
    } else {
	set cxpath "$rootid//c.DEM-Results//c.GiDOptions//i.GiDPostMode"
    }
    set GiDMode [::xmlutils::setXml $cxpath $cproperty]
    puts $fileid "\"OutputFileType\"                   : \"$GiDMode\","

    # Multifile
    if {"Fluid" in $ActiveAppList} {
	set cxpath "GeneralApplicationData//c.SDEM-Results//c.GiDOptions//i.GiDMultiFileFlag"
    } else {
	set cxpath "$rootid//c.DEM-Results//c.GiDOptions//i.GiDMultiFileFlag"
    }
    set GiDMultiFileFlag [::xmlutils::setXml $cxpath $cproperty]
    if {$GiDMultiFileFlag eq "Single"} {
	puts $fileid "\"Multifile\"                        : \"single_file\","
    } else {
	puts $fileid "\"Multifile\"                        : \"multiple_files\","
    }

    puts $fileid "\"ElementType\"                      : \"[::wkcf::GetElementType]\","
    ######################################################################################

    puts $fileid ""

    # Translational Integration Scheme
    set cxpath "$rootid//c.DEM-SolutionStrategy//i.DEM-TimeTranslationalIntegrationScheme"
    set TransIntegScheme [::xmlutils::setXml $cxpath $cproperty]

    puts $fileid "\"TranslationalIntegrationScheme\"   : \"$TransIntegScheme\","

    # Rotational Integration Scheme
    set cxpath "$rootid//c.DEM-SolutionStrategy//i.DEM-TimeRotationalIntegrationScheme"
    set RotIntegScheme [::xmlutils::setXml $cxpath $cproperty]

    puts $fileid "\"RotationalIntegrationScheme\"      : \"$RotIntegScheme\","

    # Auto Reduction of Time Step Option
    set cxpath "$rootid//c.DEM-SolutionStrategy//c.DEM-TimeParameters//i.UseAutomaticDeltaTime"
    set DTimeSelection [::xmlutils::setXml $cxpath $cproperty]

    if {$DTimeSelection eq "Fixed"} {
	puts $fileid "\"AutomaticTimestep\"                : false,"
    } else {
	puts $fileid "\"AutomaticTimestep\"                : true,"
    }
    set cxpath "$rootid//c.DEM-SolutionStrategy//c.DEM-TimeParameters//i.DeltaTimeSafetyFactor"
    set SafetyFact [::xmlutils::setXml $cxpath $cproperty]
    puts $fileid "\"DeltaTimeSafetyFactor\"            : $SafetyFact,"

    set cxpath "$rootid//c.DEM-SolutionStrategy//c.DEM-TimeParameters//i.DeltaTime"
    set MaxTimeStep [::xmlutils::setXml $cxpath $cproperty]
    puts $fileid "\"MaxTimeStep\"                      : $MaxTimeStep,"

    if {"Fluid" in $ActiveAppList} {
	set cxpath "GeneralApplicationData//c.SimulationOptions//i.EndTime"
	set TTime [::xmlutils::setXml $cxpath $cproperty]
    } else {
	set cxpath "$rootid//c.DEM-SolutionStrategy//c.DEM-TimeParameters//i.DEM-TotalTime"
	set TTime [::xmlutils::setXml $cxpath $cproperty]
    }
    puts $fileid "\"FinalTime\"                        : $TTime,"


    set cxpath "$rootid//c.DEM-SolutionStrategy//c.DEM-TimeParameters//i.DEM-ScreenInfoOutput"
    set CTime [::xmlutils::setXml $cxpath $cproperty]
    puts $fileid "\"ControlTime\"                      : $CTime,"


    set cxpath "$rootid//c.DEM-SolutionStrategy//c.DEM-TimeParameters//i.DEM-NeighbourSearchFrequency"
    set FrecTime [::xmlutils::setXml $cxpath $cproperty]
    puts $fileid "\"NeighbourSearchFrequency\"         : $FrecTime,"

    if {$KPriv(what_dempack_package) eq "C-DEMPack"} {
	set cproperty "dv"
	set cxpath "DEM//c.DEM-Options//c.DEM-AdvancedOptions//i.DEM-MaterialTestActivate"
	set material_option [::xmlutils::setXml $cxpath $cproperty]
	if {$material_option eq "Yes"} {
	    ::wkcf::WriteMatTestData $fileid
	}
    }

    if {"Fluid" in $ActiveAppList} {
	set cxpath "$rootid//c.DEM-Fluid-interaction//i.DEM-TwoWayCoupling"
	set cproperty "dv"
	set cxpath "GeneralApplicationData//c.CouplingParameters//i.CouplingLevel"
	set coupling_level_type [::xmlutils::setXml $cxpath $cproperty]
    set cproperty "dv"
    set cxpath "GeneralApplicationData//c.CouplingParameters//i.TimeAveragingType"
    set time_averaging_type [::xmlutils::setXml $cxpath $cproperty]

    set interaction_start_time [::xmlutils::setXml {GeneralApplicationData//c.CouplingParameters//i.InteractionStartTime} dv]

	set cxpath "GeneralApplicationData//c.CouplingParameters//c.HydrodynamicForceModel//i.PickIndividualForces"
	set pick_individual_forces_option [::wkcf::TranslateToBinary [::xmlutils::setXml $cxpath $cproperty]]

	set cxpath "GeneralApplicationData//c.CouplingParameters//c.HydrodynamicForceModel//i.BuoyancyForceType"
	set buoyancy_force_type [::xmlutils::setXml $cxpath $cproperty]

	set cxpath "GeneralApplicationData//c.CouplingParameters//c.HydrodynamicForceModel//i.DragForceType"
	set drag_force_type [::xmlutils::setXml $cxpath $cproperty]

	set cxpath "GeneralApplicationData//c.CouplingParameters//c.HydrodynamicForceModel//i.VirtualMassForceType"
	set virtual_mass_force_type [::xmlutils::setXml $cxpath $cproperty]

	set cxpath "GeneralApplicationData//c.CouplingParameters//c.HydrodynamicForceModel//i.SaffmanLiftForceType"
	set lift_force_type [::xmlutils::setXml $cxpath $cproperty]

	set cxpath "GeneralApplicationData//c.CouplingParameters//c.HydrodynamicForceModel//i.MagnusLiftForceType"
	set magnus_force_type [::xmlutils::setXml $cxpath $cproperty]

	set cxpath "GeneralApplicationData//c.CouplingParameters//c.HydrodynamicForceModel//i.HydrodynamicTorqueType"
	set hydro_torque_type [::xmlutils::setXml $cxpath $cproperty]

	set cxpath "GeneralApplicationData//c.CouplingParameters//c.HydrodynamicForceModel//i.DragModifier"
	set drag_modifier_type [::xmlutils::setXml $cxpath $cproperty]

	set cxpath "GeneralApplicationData//c.CouplingParameters//c.HydrodynamicForceModel//i.ViscosityModificationType"
	set viscosity_modification_type [::xmlutils::setXml $cxpath $cproperty]

	set cxpath "GeneralApplicationData//c.CouplingParameters//c.HydrodynamicForceModel//c.MaxeyRileyEquations//i.MRDragModifier"
	set MR_drag_modifier_type [::xmlutils::setXml $cxpath $cproperty]

	set cxpath "GeneralApplicationData//c.PostProcessingVariables//i.GradientCalculationType"
	set gradient_calculation_type [::xmlutils::setXml $cxpath $cproperty]

	set cxpath "GeneralApplicationData//c.PostProcessingVariables//i.VelocityLaplacianCalculationType"
	set velocity_laplacian_calculation_type [::xmlutils::setXml $cxpath $cproperty]

	if {$velocity_laplacian_calculation_type == "0"} {
	    set cxpath "GeneralApplicationData//c.CouplingParameters//c.HydrodynamicForceModel//i.IncludeFaxenTerms0"
	} elseif {$velocity_laplacian_calculation_type == "1"} {
	    set cxpath "GeneralApplicationData//c.CouplingParameters//c.HydrodynamicForceModel//i.IncludeFaxenTerms1"
	} elseif {$velocity_laplacian_calculation_type == "2"} {
	    set cxpath "GeneralApplicationData//c.CouplingParameters//c.HydrodynamicForceModel//i.IncludeFaxenTerms"
	}

	set include_faxen_terms_option [::wkcf::TranslateToBinary [::xmlutils::setXml $cxpath $cproperty]]

	if {!$pick_individual_forces_option} { # Maxey-Riley
	    if {!$include_faxen_terms_option} {
		set drag_force_type 2
		set virtual_mass_force_type 0
	    } else {
		set drag_force_type 11
		set virtual_mass_force_type 11
	    }
	    set lift_force_type 0
	    set magnus_force_type 0
	    set hydro_torque_type 0
	    set drag_modifier_type $MR_drag_modifier_type
	}

	puts $fileid ""
	puts $fileid "\"coupling_level_type\"                    : $coupling_level_type,"
    puts $fileid "\"time_averaging_type\"                    : $time_averaging_type,"
    puts $fileid "\"interaction_start_time\"                 : $interaction_start_time,"
    if {$pick_individual_forces_option} {
	puts $fileid "\"pick_individual_forces_option\"          : true,"
    } else {
	puts $fileid "\"pick_individual_forces_option\"          : false,"
    }
    set SearchNeighboursOption [::xmlutils::setXml "$rootid//c.DEM-Options//c.DEM-Physical-opts//i.SearchNeighboursOption" dv]
    if {$SearchNeighboursOption == "Yes"} {
	puts $fileid "\"do_search_neighbours\"                   : true,"
    } else {
	puts $fileid "\"do_search_neighbours\"                   : false,"
    }
    if {$include_faxen_terms_option} {
	puts $fileid "\"include_faxen_terms_option\"             : true,"
    } else {
	puts $fileid "\"include_faxen_terms_option\"             : false,"
    }
    puts $fileid "\"include_faxen_terms_option_comment\"     : \"(relevant if the Maxey Riley equation is used)\","

    puts $fileid "\"gradient_calculation_type\"              : $gradient_calculation_type,"
    puts $fileid "\"gradient_calculation_type_comment\"      : \"(Not calculated (0), volume-weighed average(1), Superconvergent recovery(2))\","

    puts $fileid "\"laplacian_calculation_type\"             : $velocity_laplacian_calculation_type,"
    puts $fileid "\"laplacian_calculation_type_comment\"     : \"(Not calculated (0), Finite element projection (1), Superconvergent recovery(2))\","
    puts $fileid "\"buoyancy_force_type\"                    : $buoyancy_force_type,"
    puts $fileid "\"buoyancy_force_type_comment\"            : \"null buoyancy (0), compute buoyancy (1)  if drag_force_type is 2 buoyancy is always parallel to gravity\","

    puts $fileid "\"drag_force_type\"                        : $drag_force_type,"
    puts $fileid "\"drag_force_type_comment\"                : \" null drag (0), Stokes (1), Weatherford (2), Ganser (3), Ishii (4), Newtonian Regime (5)\","

    puts $fileid "\"virtual_mass_force_type\"                : $virtual_mass_force_type,"
    puts $fileid "\"virtual_mass_force_type_comment\"        : \"null virtual mass force (0)\","

    puts $fileid "\"lift_force_type\"                        : $lift_force_type,"
    puts $fileid "\"lift_force_type_comment\"                : \"# null lift force (0), Saffman (1)\","

    puts $fileid "\"magnus_force_type\"                      : $magnus_force_type,"
    puts $fileid "\"magnus_force_type_comment\"              : \" null magnus force (0), Rubinow and Keller (1), Oesterle and Bui Dihn (2)\","

    puts $fileid "\"hydro_torque_type\"                      : $hydro_torque_type,"
    puts $fileid "\"hydro_torque_type_comment\"              : \"null hydrodynamic torque (0), Dennis (1)\","

    puts $fileid "\"drag_modifier_type\"                     : $drag_modifier_type,"
    puts $fileid "\"viscosity_modification_type\"            : $viscosity_modification_type,"
    puts $fileid ""

	puts $fileid "\"coupling_weighing_type\"                 : 2,"
	puts $fileid "\"coupling_weighing_type_comment\"         : \"{fluid_to_DEM, DEM_to_fluid, fluid_fraction} = {lin, lin, imposed} (-1), {lin, const, const} (0), {lin, lin, const} (1), {lin, lin, lin} (2), averaging method (3)\","

	puts $fileid "\"fluid_model_type\"                       : 1,"
	puts $fileid "\"fluid_model_type_comment\"               : \" untouched, velocity incremented by 1/fluid_fraction (0), modified mass conservation only (1)\","

	puts $fileid "\"coupling_scheme_type\"                   : \"UpdatedFluid\","
	puts $fileid "\"coupling_scheme_type_comment\"           : \" UpdatedFluid, UpdatedDEM\","

	puts $fileid "\"print_particles_results_option\"         : false,"

	puts $fileid "\"add_each_hydro_force_option\"            : true,"
	puts $fileid "\"add_each_hydro_force_option_comment\"    : \" add each of the hydrodynamic forces (drag, lift and virtual mass)\","

	puts $fileid "\"project_at_every_substep_option\"        : true,"
	puts $fileid "\"velocity_trap_option\"                   : false,"
	puts $fileid "\"inlet_option\"                           : true,"
	puts $fileid "\"manually_imposed_drag_law_option\"       : false,"

	puts $fileid "\"stationary_problem_option\"              : false,"
	puts $fileid "\"stationary_problem_option_comment\"      : \" stationary, stop calculating the fluid after it reaches the stationary state\","

	puts $fileid "\"flow_in_porous_medium_option\"           : false,"
	puts $fileid "\"flow_in_porous_medium_option_comment\"   : \" the porosity is an imposed field\","

	puts $fileid "\"flow_in_porous_DEM_medium_option\"       : false,"
	puts $fileid "\"flow_in_porous_DEM_medium_option_comment\" : \"the DEM part is kept static\","

	puts $fileid "\"embedded_option\"                        : true,"
	puts $fileid "\"embedded_option_comment\"                : \"the embedded domain tools are to be used\","

	puts $fileid "\"make_results_directories_option\"        : true,"
	puts $fileid "\"make_results_directories_option_comment\": \"results are written into a folder (../results) inside the problem folder\","

	puts $fileid "\"body_force_on_fluid_option\"             : true,"

	puts $fileid "\"print_debug_info_option\"                : false,"
	puts $fileid "\"print_debug_info_option_comment\"        : \" print a summary of global physical measures\","

	puts $fileid "\"print_particles_results_cycle\"          : 1,"
	puts $fileid "\"print_particles_results_cycle_comment\"  : \" number of 'ticks' per printing cycle\","

	puts $fileid "\"debug_tool_cycle\"                       : 10,"
	puts $fileid "\"debug_tool_cycle_comment\"               : \" number of 'ticks' per debug computations cycle\","

	puts $fileid "\"similarity_transformation_type\"         : 0,"
	puts $fileid "\"similarity_transformation_type_comment\" : \" no transformation (0), Tsuji (1)\","

	puts $fileid "\"dem_inlet_element_type\"                 : \"SphericSwimmingParticle3D\","
	puts $fileid "\"dem_inlet_element_type_comment\"         : \" SphericParticle3D, SphericSwimmingParticle3D\","

	puts $fileid "\"drag_modifier_type\"                     : 2,"
	puts $fileid "\"drag_modifier_type_comment\"             : \" Hayder (2), Chien (3) # problemtype option\","

	puts $fileid "\"drag_porosity_correction_type\"          : 0,"
	puts $fileid "\"drag_porosity_correction_type_comment\"  : \" No correction (0), Richardson and Zaki (1)\","

	puts $fileid "\"min_fluid_fraction\"                     : 0.2,"
	puts $fileid "\"initial_drag_force\"                     : 0.0,"
	puts $fileid "\"drag_law_slope\"                         : 0.0,"
	puts $fileid "\"power_law_tol\"                          : 0.0,"
	puts $fileid "\"model_over_real_diameter_factor\"        : 1.0,"
	puts $fileid "\"model_over_real_diameter_factor_comment\": \" not active if similarity_transformation_type = 0\","

	puts $fileid "\"max_pressure_variation_rate_tol\"        : 1e-3,"
	puts $fileid "\"max_pressure_variation_rate_tol_comment\": \" for stationary problems, criterion to stop the fluid calculations\","

	puts $fileid "\"time_steps_per_stationarity_step\"       : 15,"
	puts $fileid "\"time_steps_per_stationarity_step_comment\": \" number of fluid time steps between consecutive assessment of stationarity steps\","

	puts $fileid "\"meso_scale_length\"                      : 0.2,"
	puts $fileid "\"meso_scale_length_comment\"              : \" the radius of the support of the averaging function for homogenization (<=0 for automatic calculation)\","

	puts $fileid "\"shape_factor\"                           : 0.5,"
	if {[dict get [::wkcf::GetFluidMaterialProperties "Fluid"] "NonNewtonianFluid"] eq "No"} {
	    puts $fileid "\"non_newtonian_option\"                   : false,"
	} else {
	    puts $fileid "\"non_newtonian_option\"                   : true,"
	    puts $fileid "\"yield_stress\"                           : [dict get [::wkcf::GetFluidMaterialProperties "Fluid"] "YieldStress"],"
	    puts $fileid "\"regularization_coefficient\"             : [dict get [::wkcf::GetFluidMaterialProperties "Fluid"] "BinghamSmoother"],"
	    puts $fileid "\"power_law_k\"                            : [dict get [::wkcf::GetFluidMaterialProperties "Fluid"] "PowerLawK"],"
	    puts $fileid "\"power_law_n\"                            : [dict get [::wkcf::GetFluidMaterialProperties "Fluid"] "PowerLawN"],"
	}
    }
    # SWIMMING-SPECIFIC SECTION ENDS ###########################################################################

    puts $fileid ""
    ::wkcf::WritePostProcessDataForJson $fileid
    puts $fileid ""
    # Get the problem name
    set PName [::KUtils::GetPaths "PName"]
    puts $fileid "\"problem_name\" : \"${PName}\""

    puts $fileid "}"

    # Get the kratos path
    set KratosPath [::xmlutils::setXml {GeneralApplicationData//c.ProjectConfiguration//i.KratosPath} dv]
    set KratosPath [file native $KratosPath]

    close $fileid
}

proc ::wkcf::WriteDemNodalVariables { AppId } {
    variable filechannel
    variable dprops
    variable ndime

    set nodelistbuffer [dict create ]

    # Write RADIUS
    set nodvar [::xmlutils::getKKWord "Applications/$AppId" "Radius"]

    set basexpath "$AppId//c.DEM-Elements"

    set gelemlist [::wkcf::GetDEMActiveElements]
    set nodelistbuffer ""

    foreach celemid $gelemlist {

	set elembasexpath "$basexpath//c.DEM-Element"
	# Get the group node list
	set grouplist [::xmlutils::setXmlContainerIds $elembasexpath]
	if {$grouplist != ""} {

	    foreach cgroupid $grouplist {

		if {$cgroupid != ""} {
		    if {[GiD_EntitiesGroups get $cgroupid elements -count]} {

		        GiD_File fprintf $filechannel "%s" "Begin NodalData $nodvar  // GUI group identifier: $cgroupid Elementid $celemid"
		        # De cada elemento necesitamos su nodo y su radio
		        foreach elem [GiD_EntitiesGroups get $cgroupid elements] {
		            set r 0
		            set info [GiD_Info list_entities -more elements $elem]
		            if { ![regexp {Type=([^ ]+) Nnode=([0-9]+) Volume=([^ ]+)} $info dummy type nnode volume] } {
		                set type ?
		                set nnode 0
		                set volume 0
		            }
		            if { $nnode > 0} {
		                set nodeid [lindex [GiD_Mesh get element $elem] 3]
		                dict lappend nodelistbuffer $cgroupid $nodeid
		                if {$ndime eq "3D"} {
		                    set r [expr {pow(3.0/4.0*$volume/$::MathUtils::Pi,1.0/3.0)}]
		                } else {
		                    set r [expr {sqrt($volume/$::MathUtils::Pi)}]
		                }
		                set r [GiD_FormatReal "%10.5e" $r]
		                GiD_File fprintf $filechannel "%s" "$nodeid 0 $r"
		            }
		        }

		        GiD_File fprintf $filechannel "%s" "End NodalData"
		        GiD_File fprintf $filechannel "%s" ""
		    }

		}
	    }
	}
    }

    # Write Cohesive Groups
    set nodvar [::xmlutils::getKKWord "Applications/$AppId" "Cohesive_Group"]

    set basexpath "$AppId//c.DEM-Cohesivegroup"
    set grouplist [::xmlutils::setXmlContainerIds $basexpath]

    foreach cgroupid $grouplist {
	if {$cgroupid != ""} {
	    if {[GiD_EntitiesGroups get $cgroupid elements -count]} {
		GiD_File fprintf $filechannel "%s" "Begin NodalData $nodvar  // GUI group identifier: $cgroupid"
		set cproperty "dv"
		set cxpath "$elembasexpath//c.[list ${cgroupid}]//c.Properties//i.SetActive"
		set active_or_not [::xmlutils::setXml $cxpath $cproperty]
		if {$active_or_not=="No"} {
		    continue
		}
		set elembasexpath "$basexpath//c.$cgroupid//c.Properties//i.CohGroup"
		# Get the group cohesive group
		set cohgroup [::xmlutils::setXml $elembasexpath "dv" ]

		if {[dict exists $nodelistbuffer $cgroupid]} {
		    foreach nodeid [dict get $nodelistbuffer $cgroupid] {
		        GiD_File fprintf $filechannel "%s" "$nodeid 0 $cohgroup"
		    }
		} else {
		    foreach elem [GiD_EntitiesGroups get $cgroupid elements] {
		        # Get the node from the element
		        set nodeid [lindex [GiD_Mesh get element $elem] 3]
		        dict lappend nodelistbuffer $cgroupid $nodeid
		        GiD_File fprintf $filechannel "%s" "$nodeid 0 $cohgroup"
		    }
		}
		GiD_File fprintf $filechannel "%s" "End NodalData"
		GiD_File fprintf $filechannel "%s" ""
	    }
	}
    }

    # Write All Material Dependent Properties; pe. Density, young modulus...
    #     set proplist [list "Density" "YoungModulus" "PoissonRatio" "ParticleFrictionAngle" "CoefficientOfRestitution" "Color"]
    #     foreach cproperty $proplist {
	#         set nodvar [::xmlutils::getKKWord "Applications/$AppId" $cproperty]
	#
	#         set basexpath "$AppId//c.DEM-Elements"
	#
	#         foreach celemid $gelemlist {
	    #
	    #             set elembasexpath "$basexpath//c.DEM-Element"
	    #             # Get the group node list
	    #             set grouplist [::xmlutils::setXmlContainerIds $elembasexpath]
	    #             if {$grouplist != ""} {
		#                 foreach cgroupid $grouplist {
		    #                     if {$cgroupid != ""} {
		        #                         if {[GiD_EntitiesGroups get $cgroupid elements -count]} {
		            #                             set matxpath "$elembasexpath//c.$cgroupid//c.Properties//i.Material"
		            #                             set material [::xmlutils::setXml $matxpath "dv" ]
		            #
		            #                             set cxpath "DEMMaterial//m.$material//p.$cproperty"
		            #                             set propvalue [::xmlutils::setXml $cxpath "dv" "read" "" "mat"]
		            #
		            #                             if {$cproperty eq "CoefficientOfRestitution"} {
		                #                                 if {$propvalue == 0.0} {
		                    #                                     set propvalue 1.0
		                    #                                 } else {
		                    #                                     set propvalue [expr {log($propvalue)}]
		                    #                                 }
		                #                             }
		            #                             set propvalue [GiD_FormatReal "%10.5e" $propvalue]
		            #                             GiD_File fprintf $filechannel "%s" "Begin NodalData $nodvar  // GUI group identifier: $cgroupid"
		            #                             # De cada elemento necesitamos su nodo y su radio
		            #                             if {[dict exists $nodelistbuffer $cgroupid]} {
		                #                                 foreach nodeid [dict get $nodelistbuffer $cgroupid] {
		                    #                                     GiD_File fprintf $filechannel "%s" "$nodeid 0 $propvalue"
		                    #                                 }
		                #                             } else {
		                #                                 foreach elem [GiD_EntitiesGroups get $cgroupid elements] {
		                    #                                     # Get the node from the element
		                    #                                     set nodeid [lindex [GiD_Mesh get element $elem] 3]
		                    #                                     dict lappend nodelistbuffer $cgroupid $nodeid
		                    #                                     GiD_File fprintf $filechannel "%s" "$nodeid 0 $propvalue"
		                    #                                 }
		                #                             }
		            #
		            #                             GiD_File fprintf $filechannel "%s" "End NodalData"
		            #                             GiD_File fprintf $filechannel "%s" ""
		            #                         }
		        #                     }
		    #                 }
		#             }
	    #         }
	#     }
}

proc ::wkcf::WriteDemNodalVariables2 {AppId filechannel} {

    global KPriv
    # Write RADIUS
    variable ndime
    set nodvar [::xmlutils::getKKWord Applications/$AppId Radius]
    set basexpath $AppId//c.DEM-Elements
    set gelemlist [::wkcf::GetDEMActiveElements]
    set elembasexpath $basexpath//c.DEM-Element
    if {$ndime eq "2D"} {
	lassign [lindex [GiD_Info Mesh elements Circle -array] 0] type element_ids element_nodes element_materials element_radii
	foreach element_id $element_ids element_node [lindex $element_nodes 0] element_radius $element_radii {
	    set element_data($element_id) [list $element_node $element_radius]}
    } else {
	set a [GiD_Info Mesh elements Sphere -array]
	set b [lindex $a 0]
	lassign $b type element_ids element_nodes element_materials element_radii
	foreach element_id $element_ids element_node [lindex $element_nodes 0] element_radius $element_radii {
	    set element_data($element_id) [list $element_node $element_radius]}
    }
    #release memory
    foreach variable_name {element_ids element_nodes element_materials element_radii} {
	unset $variable_name
    }
    # Get the group node list
    set list_of_active_dem_elements ""
    set grouplist [::xmlutils::setXmlContainerIds $elembasexpath]
    if {[llength $grouplist]} {
	foreach celemid $gelemlist {
	    foreach cgroupid $grouplist {

		set cproperty "dv"
		set cxpath "$elembasexpath//c.[list ${cgroupid}]//c.Properties//i.SetActive"
		set active_or_not [::xmlutils::setXml $cxpath $cproperty]
		if {$active_or_not=="No"} {
		    continue
		}

		if {[GiD_EntitiesGroups get $cgroupid elements -count]} {
		    GiD_File fprintf $filechannel "Begin NodalData $nodvar  // GUI group identifier: $cgroupid Elementid $celemid"
		    set element_list [GiD_EntitiesGroups get $cgroupid elements]
		    lappend list_of_active_dem_elements $element_list
		    foreach element_id $element_list {
		        lassign $element_data($element_id) element_node element_radius
		        lappend group_nodes($cgroupid) $element_node
		        GiD_File fprintf $filechannel "$element_node 0 $element_radius"
		    }
		    GiD_File fprintf $filechannel "End NodalData"
		    GiD_File fprintf $filechannel ""
		}
	    }
	}
    }
    # Write Cohesive Groups
    set nodvar [::xmlutils::getKKWord Applications/$AppId Cohesive_Group]
    set basexpath $AppId//c.DEM-Cohesivegroup
    set grouplist [::xmlutils::setXmlContainerIds $basexpath]
    foreach cgroupid $grouplist {
	if {[GiD_EntitiesGroups get $cgroupid elements -count]} {
	    set cproperty "dv"
	    set cxpath "$AppId//c.DEM-Elements//c.DEM-Element//c.[list ${cgroupid}]//c.Properties//i.SetActive"
	    set active_or_not [::xmlutils::setXml $cxpath $cproperty]
	    if {$active_or_not=="No"} {
		continue
	    }
	    GiD_File fprintf $filechannel "Begin NodalData $nodvar  // GUI group identifier: $cgroupid"
	    set elembasexpath $basexpath//c.$cgroupid//c.Properties//i.CohGroup
	    set cohesive_group [::xmlutils::setXml $elembasexpath dv]
	    if { [info exists group_nodes($cgroupid)] } {
		foreach node_id $group_nodes($cgroupid) {
		    GiD_File fprintf $filechannel "$node_id 0 $cohesive_group"
		}
	    } else {
		# Get the node from the element
		foreach element_id [GiD_EntitiesGroups get $cgroupid elements] {
		    set node_id [lindex $element_data($element_id) 0]
		    lappend group_nodes($cgroupid) $node_id
		    GiD_File fprintf $filechannel "$node_id 0 $cohesive_group"
		}
	    }
	    GiD_File fprintf $filechannel "End NodalData"
	    GiD_File fprintf $filechannel ""
	}
    }
    if {$KPriv(what_dempack_package) eq "C-DEMPack"} {

	# Write Skin Sphere
	if {[GiD_Groups exists SKIN_SPHERE_DO_NOT_DELETE]} {
	if {$ndime eq "2D"} {
	    set skin_element_ids [GiD_EntitiesGroups get SKIN_SPHERE_DO_NOT_DELETE all_mesh -element_type circle] ; # Get the ids of elements in SKIN_SPHERE
	} else {
	    set skin_element_ids [GiD_EntitiesGroups get SKIN_SPHERE_DO_NOT_DELETE all_mesh -element_type sphere]
	}
	} else {
	set skin_element_ids [list]
	}

	set only_skin_elems_ids [lindex $skin_element_ids 1]
	set list_of_active_dem_elements [regsub -all {\{|\}} $list_of_active_dem_elements ""]

	foreach active_element $list_of_active_dem_elements {
	    set temporal_list($active_element) 1
	}

	set elements_in_common_list ""
	foreach skin_element $only_skin_elems_ids {
	    if {[info exists temporal_list($skin_element)]} {
		lappend elements_in_common_list $skin_element
	    }
	}

	foreach element_id $elements_in_common_list { ; # Here we loop on each of the elements by id
	    set element_nodes_id($element_id) [lindex [GiD_Mesh get element $element_id] 3] ; # We get the nodes of the element
	}
	GiD_File fprintf $filechannel "Begin NodalData SKIN_SPHERE"
	foreach element_id $elements_in_common_list {
	    GiD_File fprintf $filechannel "$element_nodes_id($element_id) 0 1"
	}
	GiD_File fprintf $filechannel "End NodalData"
	GiD_File fprintf $filechannel ""
    }
    #    # Write All Material Dependent Properties; pe. Density, young modulus...
    #     set basexpath $AppId//c.DEM-Elements
    #     set elembasexpath $basexpath//c.DEM-Element
    #     set grouplist [::xmlutils::setXmlContainerIds $elembasexpath]
    #     foreach cproperty {Density YoungModulus PoissonRatio ParticleFrictionAngle CoefficientOfRestitution Color} {
	#         set nodvar [::xmlutils::getKKWord Applications/$AppId $cproperty]
	#         foreach celemid $gelemlist {
	    #             foreach cgroupid $grouplist {
		#                 if {[GiD_EntitiesGroups get $cgroupid elements -count]} {
		    #                     set matxpath $elembasexpath//c.$cgroupid//c.Properties//i.Material
		    #                     set material [::xmlutils::setXml $matxpath dv]
		    #                     set cxpath DEMMaterial//m.$material//p.$cproperty
		    #                     set propvalue [::xmlutils::setXml $cxpath dv read "" mat]
		    #                     if {$cproperty eq "CoefficientOfRestitution"} {
		        #                         if {$propvalue == 0.0} {
		            #                             set propvalue 1.0
		            #                         } else {
		            #                             set propvalue [expr {log($propvalue)}]
		            #                         }
		        #                     }
		    #                     GiD_File fprintf $filechannel "Begin NodalData $nodvar  // GUI group identifier: $cgroupid"
		    #                     if { [info exists group_nodes($cgroupid)] } {
		        #                         foreach node_id $group_nodes($cgroupid) {
		            #                             GiD_File fprintf $filechannel "$node_id 0 $propvalue"
		            #                         }
		        #                     } else {
		        #                         # Get the node from the element
		        #                         foreach element_id [GiD_EntitiesGroups get $cgroupid elements] {
		            #                             set node_id [lindex $element_data($element_id) 0]
		            #                             lappend group_nodes($cgroupid) $node_id
		            #                             GiD_File fprintf $filechannel "$node_id 0 $propvalue"
		            #                         }
		        #                     }
		    #                     GiD_File fprintf $filechannel "End NodalData"
		    #                     GiD_File fprintf $filechannel ""
		    #                 }
		#             }
	    #         }
	#     }
    return 0
}

proc ::wkcf::WriteDEMConditionProperties {AppId cgroupid cpropvalue} {
    variable dprops; variable pflag
    variable filechannel

    # For debug
    if {!$pflag} {
	set inittime [clock seconds]
    }

    # Kratos key word xpath
    set kxpath "Applications/$AppId"

    set proplist [list "IMPOSED_VELOCITY_X" "IMPOSED_VELOCITY_X_VALUE" \
	    "IMPOSED_VELOCITY_Y" "IMPOSED_VELOCITY_Y_VALUE" \
	    "IMPOSED_VELOCITY_Z" "IMPOSED_VELOCITY_Z_VALUE" \
	    "IMPOSED_ANGULAR_VELOCITY_X" "IMPOSED_ANGULAR_VELOCITY_X_VALUE" \
	    "IMPOSED_ANGULAR_VELOCITY_Y" "IMPOSED_ANGULAR_VELOCITY_Y_VALUE" \
	    "IMPOSED_ANGULAR_VELOCITY_Z" "IMPOSED_ANGULAR_VELOCITY_Z_VALUE" \
	    "VELOCITY_START_TIME" "VELOCITY_STOP_TIME" \
	    "FORCE_INTEGRATION_GROUP" \
	    "TOP" "BOTTOM" \
	    "INITIAL_VELOCITY_X_VALUE" "INITIAL_VELOCITY_Y_VALUE" "INITIAL_VELOCITY_Z_VALUE" \
	    "INITIAL_ANGULAR_VELOCITY_X_VALUE" "INITIAL_ANGULAR_VELOCITY_Y_VALUE" \
	    "INITIAL_ANGULAR_VELOCITY_Z_VALUE" "IDENTIFIER"]

    GiD_File fprintf $filechannel "%s" "  Begin SubModelPartData // GUI Groupid: $cgroupid"

    foreach cprop $proplist cval $cpropvalue {
	set my_data($cprop) $cval
    }

    foreach cprop [array names my_data] {
	if {$cprop == "IMPOSED_VELOCITY_X_VALUE" && !$my_data(IMPOSED_VELOCITY_X)} {
	    continue
	}
	if {$cprop == "IMPOSED_VELOCITY_Y_VALUE" && !$my_data(IMPOSED_VELOCITY_Y)} {
	    continue
	}
	if {$cprop == "IMPOSED_VELOCITY_Z_VALUE" && !$my_data(IMPOSED_VELOCITY_Z)} {
	    continue
	}
	if {$cprop == "IMPOSED_ANGULAR_VELOCITY_X_VALUE" && !$my_data(IMPOSED_ANGULAR_VELOCITY_X)} {
	    continue
	}
	if {$cprop == "IMPOSED_ANGULAR_VELOCITY_Y_VALUE" && !$my_data(IMPOSED_ANGULAR_VELOCITY_Y)} {
	    continue
	}
	if {$cprop == "IMPOSED_ANGULAR_VELOCITY_Z_VALUE" && !$my_data(IMPOSED_ANGULAR_VELOCITY_Z)} {
	    continue
	}
	if {$cprop == "INITIAL_VELOCITY_X_VALUE" && !$my_data(INITIAL_VELOCITY_X_VALUE)} {
	    continue
	}
	if {$cprop == "INITIAL_VELOCITY_Y_VALUE" && !$my_data(INITIAL_VELOCITY_Y_VALUE)} {
	    continue
	}
	if {$cprop == "INITIAL_VELOCITY_Z_VALUE" && !$my_data(INITIAL_VELOCITY_Z_VALUE)} {
	    continue
	}
	if {$cprop == "INITIAL_ANGULAR_VELOCITY_X_VALUE" && !$my_data(INITIAL_ANGULAR_VELOCITY_X_VALUE)} {
	    continue
	}
	if {$cprop == "INITIAL_ANGULAR_VELOCITY_Y_VALUE" && !$my_data(INITIAL_ANGULAR_VELOCITY_Y_VALUE)} {
	    continue
	}
	if {$cprop == "INITIAL_ANGULAR_VELOCITY_Z_VALUE" && !$my_data(INITIAL_ANGULAR_VELOCITY_Z_VALUE)} {
	    continue
	}
	set fix_list [list "IMPOSED_VELOCITY_X" "IMPOSED_VELOCITY_Y" "IMPOSED_VELOCITY_Z" "IMPOSED_ANGULAR_VELOCITY_X" \
		"IMPOSED_ANGULAR_VELOCITY_Y" "IMPOSED_ANGULAR_VELOCITY_Z"]
	if {[lsearch -exact $fix_list $cprop] >= 0} {continue}

	GiD_File fprintf $filechannel "%s" "  $cprop $my_data($cprop)"
    }

    GiD_File fprintf $filechannel "%s" "  End SubModelPartData"

    # For debug
    if {!$::wkcf::pflag} {
	set endtime [clock seconds]
	set ttime [expr {$endtime-$inittime}]
	WarnWinText "Obtaining DEM-Inlet group using the new mesh format: [::KUtils::Duration $ttime]"
    }
}

proc ::wkcf::WriteDEMElementMeshProperties {AppId} {

    variable filechannel
    variable ndime
    set basexpath "$AppId//c.DEM-Conditions//c.DEM-VelocityBC"
    set gproplist [::xmlutils::setXmlContainerIds $basexpath]
    set pathtoDEMElems "$AppId//c.DEM-Elements//c.DEM-Element"

    foreach cgroupid $gproplist {
	set cproperty "dv"
	set cxpath "${basexpath}//c.[list ${cgroupid}]//i.SetActive"
	set active_or_not [::xmlutils::setXml $cxpath $cproperty]
	set cproperty "dv"
	set cxpath "${pathtoDEMElems}//c.[list ${cgroupid}]//c.Properties//i.SetActive"
	set dem_elem_active_or_not [::xmlutils::setXml $cxpath $cproperty]

	if {$active_or_not=="No" || $dem_elem_active_or_not=="No"} {
	    continue
	}

	set cproperty "dv"
	set cxpath "${basexpath}//c.[list ${cgroupid}]//i.DEM-VelocityBCMotion"
	set motion_type [::xmlutils::setXml $cxpath $cproperty]
    # Get the group node list
    set elist [GiD_EntitiesGroups get $cgroupid elements]
    if {[llength $elist]} {
    # Write all nodes for this group in increasing order
    set nodeslist [list]
    foreach eid $elist {
	if { $ndime eq "2D"} {
	    set nodeid [lindex [GiD_Info Mesh Elements circle $eid $eid] 1]
	} else {
	    set nodeid [lindex [GiD_Info Mesh Elements sphere $eid $eid] 1]
	}
	#if { $nodeid == "" } {
	    #    return [list 1 [= "Some elements in this group are not spheres! Check Group %s" $cgroupid]]
	    #}
	lappend nodeslist $nodeid
    }
    set nodeslist [lsort -integer -unique $nodeslist]
    set cproperty "dv"
    set cxpath "${basexpath}//c.[list ${cgroupid}]//c.LinearVelocity//i.LinearVelocityX"
    set LinearVelocityX [::xmlutils::setXml $cxpath $cproperty]
    set cxpath "${basexpath}//c.[list ${cgroupid}]//c.LinearVelocity//i.LinearVelocityY"
    set LinearVelocityY [::xmlutils::setXml $cxpath $cproperty]
    set cxpath "${basexpath}//c.[list ${cgroupid}]//c.LinearVelocity//i.LinearVelocityZ"
    set LinearVelocityZ [::xmlutils::setXml $cxpath $cproperty]
    set cxpath "${basexpath}//c.[list ${cgroupid}]//c.LinearVelocity//i.LinearPeriodic"
    set IsPeriodic [::xmlutils::setXml $cxpath $cproperty]
    if {$IsPeriodic=="Yes"} {
	set cxpath "${basexpath}//c.[list ${cgroupid}]//c.LinearVelocity//i.LinearPeriod"
	set Period [::xmlutils::setXml $cxpath $cproperty]
    } else {
	set Period "0.0"
    }
    set cxpath "${basexpath}//c.[list ${cgroupid}]//c.LinearVelocity//i.LinearStartTime"
    set LinearVelocityStartTime [::xmlutils::setXml $cxpath $cproperty]
    set cxpath "${basexpath}//c.[list ${cgroupid}]//c.LinearVelocity//i.LinearEndTime"
    set LinearVelocityEndTime [::xmlutils::setXml $cxpath $cproperty]
    set cxpath "${basexpath}//c.[list ${cgroupid}]//c.AngularVelocity//i.AngularVelocityX"
    set AngularVelocityX [::xmlutils::setXml $cxpath $cproperty]
    set cxpath "${basexpath}//c.[list ${cgroupid}]//c.AngularVelocity//i.AngularVelocityY"
    set AngularVelocityY [::xmlutils::setXml $cxpath $cproperty]
    set cxpath "${basexpath}//c.[list ${cgroupid}]//c.AngularVelocity//i.AngularVelocityZ"
    set AngularVelocityZ [::xmlutils::setXml $cxpath $cproperty]
    set cxpath "${basexpath}//c.[list ${cgroupid}]//c.AngularVelocity//i.CenterOfRotationX"
    set CenterOfRotationX [::xmlutils::setXml $cxpath $cproperty]
    set cxpath "${basexpath}//c.[list ${cgroupid}]//c.AngularVelocity//i.CenterOfRotationY"
    set CenterOfRotationY [::xmlutils::setXml $cxpath $cproperty]
    set cxpath "${basexpath}//c.[list ${cgroupid}]//c.AngularVelocity//i.CenterOfRotationZ"
    set CenterOfRotationZ [::xmlutils::setXml $cxpath $cproperty]
    set cxpath "${basexpath}//c.[list ${cgroupid}]//c.AngularVelocity//i.AngularPeriodic"
    set AngularIsPeriodic [::xmlutils::setXml $cxpath $cproperty]
    if {$AngularIsPeriodic=="Yes"} {
	set cxpath "${basexpath}//c.[list ${cgroupid}]//c.AngularVelocity//i.AngularPeriod"
	set AngularPeriod [::xmlutils::setXml $cxpath $cproperty]
    } else {
	set AngularPeriod "0.0"
    }
    set cxpath "${basexpath}//c.[list ${cgroupid}]//c.AngularVelocity//i.AngularStartTime"
    set AngularVelocityStartTime [::xmlutils::setXml $cxpath $cproperty]
    set cxpath "${basexpath}//c.[list ${cgroupid}]//c.AngularVelocity//i.AngularEndTime"
    set AngularVelocityEndTime [::xmlutils::setXml $cxpath $cproperty]
    set RigidBodyMotion 1
    set TableNumber 0
    set TableVelocityComponent 0

    variable dem_group_mesh_property_number
    incr dem_group_mesh_property_number

	if {$motion_type=="FromATable"} {
	set RigidBodyMotion 0
	set TableNumber $dem_group_mesh_property_number
	set TableVelocityComponent [::xmlutils::setXml "${basexpath}//c.[list ${cgroupid}]//i.TableVelocityComponent" dv]
	}
	if {$motion_type=="None"} {
	foreach {LinearVelocityX LinearVelocityY LinearVelocityZ AngularVelocityX AngularVelocityY AngularVelocityZ} {0.0 0.0 0.0 0.0 0.0 0.0} {}
    }
    if {$motion_type=="FixedDOFs"} {
	set RigidBodyMotion 0
    }

    GiD_File fprintf $filechannel "Begin SubModelPart $dem_group_mesh_property_number // DEM-Element-RigidBodyMotion. Group name: $cgroupid"
    GiD_File fprintf $filechannel "  Begin SubModelPartData // DEM-Element-RigidBodyMotion. Group name: $cgroupid"
    if {$motion_type=="None" || $motion_type=="LinearPeriodic"} {
	GiD_File fprintf $filechannel "  LINEAR_VELOCITY \[3\] ($LinearVelocityX,$LinearVelocityY,$LinearVelocityZ)"
	GiD_File fprintf $filechannel "  VELOCITY_PERIOD $Period"
	GiD_File fprintf $filechannel "  ANGULAR_VELOCITY \[3\] ($AngularVelocityX,$AngularVelocityY,$AngularVelocityZ)"
	GiD_File fprintf $filechannel "  ROTATION_CENTER \[3\] ($CenterOfRotationX,$CenterOfRotationY,$CenterOfRotationZ)"
	GiD_File fprintf $filechannel "  ANGULAR_VELOCITY_PERIOD $AngularPeriod"
	GiD_File fprintf $filechannel "  VELOCITY_START_TIME $LinearVelocityStartTime"
	GiD_File fprintf $filechannel "  VELOCITY_STOP_TIME $LinearVelocityEndTime"
	GiD_File fprintf $filechannel "  ANGULAR_VELOCITY_START_TIME $AngularVelocityStartTime"
	GiD_File fprintf $filechannel "  ANGULAR_VELOCITY_STOP_TIME $AngularVelocityEndTime"
    }
    GiD_File fprintf $filechannel "  RIGID_BODY_MOTION $RigidBodyMotion"
    GiD_File fprintf $filechannel "  IDENTIFIER $cgroupid"
    GiD_File fprintf $filechannel "  End SubModelPartData"
    GiD_File fprintf $filechannel "  Begin SubModelPartNodes"
    foreach nid $nodeslist {
	#lassign [GiD_Mesh get node $nid] layer x y z
	GiD_File fprintf $filechannel "  $nid"
    }
    GiD_File fprintf $filechannel "  End SubModelPartNodes"
    GiD_File fprintf $filechannel "End SubModelPart"
    GiD_File fprintf $filechannel ""

    if {$motion_type=="FromATable"} {
	set filename [::xmlutils::setXml "${basexpath}//c.[list ${cgroupid}]//i.VelocitiesFilename" dv]
	GiD_File fprintf $filechannel "Begin Table $TableNumber TIME VELOCITY"
	set file_open [open [file native [file join [::KUtils::GetPaths "PDir"] $filename]] r]
	set file_data [read $file_open]
	close $file_open
	GiD_File fprintf -nonewline $filechannel $file_data
	GiD_File fprintf $filechannel "End Table"
	GiD_File fprintf $filechannel ""
    }
	}
    }
    GiD_File fprintf $filechannel "Begin Table 0 TIME VELOCITY"
    GiD_File fprintf $filechannel "0.0  0.0"
    GiD_File fprintf $filechannel "1.0  0.0"
    GiD_File fprintf $filechannel "End Table"
    GiD_File fprintf $filechannel ""
}
proc ::wkcf::WriteDSOLIDContactKinematics {AppId} {
    variable demfemchannel
    set rootid DSOLID

    set basexpath "$rootid//c.DEM-Conditions//c.Displacements"
    set gproplist [::xmlutils::setXmlContainerIds $basexpath]
    foreach cgroupid $gproplist {

	# Get the group node list
  set nlist [list]
	lassign [::wkcf::GetDSOLIDBoundaryConditionNodes $AppId $cgroupid] fail msg nlist
	if {$fail == 1} {

	    return [list 1 $msg]
	}

	if {[llength $nlist]} {

	    set cproperty "dv"
	    set cxpath "${basexpath}//c.[list ${cgroupid}]//c.Values//i.Vx"
	    set LinearX [::xmlutils::setXml $cxpath $cproperty]
	    set cxpath "${basexpath}//c.[list ${cgroupid}]//c.Values//i.Vy"
	    set LinearY [::xmlutils::setXml $cxpath $cproperty]
	    set cxpath "${basexpath}//c.[list ${cgroupid}]//c.Values//i.Vz"
	    set LinearZ [::xmlutils::setXml $cxpath $cproperty]

	    set cxpath "${basexpath}//c.[list ${cgroupid}]//c.Activation//i.Ax"
	    set FixX [::xmlutils::setXml $cxpath $cproperty]
	    set cxpath "${basexpath}//c.[list ${cgroupid}]//c.Activation//i.Ay"
	    set FixY [::xmlutils::setXml $cxpath $cproperty]
	    set cxpath "${basexpath}//c.[list ${cgroupid}]//c.Activation//i.Az"
	    set FixZ [::xmlutils::setXml $cxpath $cproperty]

	    variable contact_property_number 3
	    #set contact_property_number [expr $contact_property_number + 1 ]


	    GiD_File fprintf $demfemchannel "Begin NodalData DISPLACEMENT_X"
	    foreach nodeid $nlist {
		GiD_File fprintf $demfemchannel "$nodeid $FixX $LinearX"
	    }
	    GiD_File fprintf $demfemchannel "End NodalData"
	    GiD_File fprintf $demfemchannel ""

	    GiD_File fprintf $demfemchannel "Begin NodalData DISPLACEMENT_Y"
	    foreach nodeid $nlist {
		GiD_File fprintf $demfemchannel "$nodeid $FixY $LinearY"
	    }
	    GiD_File fprintf $demfemchannel "End NodalData"
	    GiD_File fprintf $demfemchannel ""

	    GiD_File fprintf $demfemchannel "Begin NodalData DISPLACEMENT_Z"
	    foreach nodeid $nlist {
		GiD_File fprintf $demfemchannel "$nodeid $FixZ $LinearZ"
	    }
	    GiD_File fprintf $demfemchannel "End NodalData"
	    GiD_File fprintf $demfemchannel ""


	}
    }
}

proc ::wkcf::WriteDSOLIDVolumeAccelerationOnNodes {AppId} {
    variable demfemchannel
    set rootid DSOLID

    set basexpath "$rootid//c.Solid-Elements//c.Solid-Element"

    set gproplist [::xmlutils::setXmlContainerIds $basexpath]
    foreach cgroupid $gproplist {


   # Get the group node list
   set nlist [list]
   lassign [::wkcf::GetDSOLIDGroupNodes $AppId $cgroupid] fail msg nlist

    set cproperty "dv"
    set cxpathtogravity "DEM//c.DEM-Options//c.DEM-Physical-opts"

    set cxpath "$cxpathtogravity//i.GravityValue"
    set usergravitymod [::xmlutils::setXml $cxpath $cproperty]
    set cxpath "$cxpathtogravity//i.Cx"
    set gravityCx [::xmlutils::setXml $cxpath $cproperty]
    set cxpath "$cxpathtogravity//i.Cy"
    set gravityCy [::xmlutils::setXml $cxpath $cproperty]
    set cxpath "$cxpathtogravity//i.Cz"
    set gravityCz [::xmlutils::setXml $cxpath $cproperty]

    set gravityvector [list $gravityCx $gravityCy $gravityCz]
    set gravitymod [::MathUtils::VectorModulus $gravityvector]

    if {$gravitymod==0.0} {
	W "The null vector is not a valid direction for the gravity force!"
	W "Computations will proceed anyway...\nPlease make sure that value was intentional..."
    }

    if {$gravitymod !=0.0} {
  set endgravityCx [expr {($gravityCx*$usergravitymod)/$gravitymod}]
  set endgravityCy [expr {($gravityCy*$usergravitymod)/$gravitymod}]
  set endgravityCz [expr {($gravityCz*$usergravitymod)/$gravitymod}]
    } else {
  set endgravityCx $gravityCx
  set endgravityCy $gravityCy
  set endgravityCz $gravityCz
    }

  if {$fail == 1} {
      return [list 1 $msg]
  }
  if {[llength $nlist]} {

      GiD_File fprintf $demfemchannel "Begin NodalData VOLUME_ACCELERATION_X"
      foreach nodeid $nlist {
    GiD_File fprintf $demfemchannel "$nodeid 1 $endgravityCx"
      }
      GiD_File fprintf $demfemchannel "End NodalData"
      GiD_File fprintf $demfemchannel ""

      GiD_File fprintf $demfemchannel "Begin NodalData VOLUME_ACCELERATION_Y"
      foreach nodeid $nlist {
    GiD_File fprintf $demfemchannel "$nodeid 1 $endgravityCy"
      }
      GiD_File fprintf $demfemchannel "End NodalData"
      GiD_File fprintf $demfemchannel ""

      GiD_File fprintf $demfemchannel "Begin NodalData VOLUME_ACCELERATION_Z"
      foreach nodeid $nlist {
    GiD_File fprintf $demfemchannel "$nodeid 1 $endgravityCz"
      }
      GiD_File fprintf $demfemchannel "End NodalData"
      GiD_File fprintf $demfemchannel ""

}
}
}

proc ::wkcf::WriteDEMFEMWallMeshProperties {AppId} {
    variable demfemchannel
    variable demfem_ref_to_props_number
    set demfem_ref_to_props_number  0
	variable demfem_motion_table
    set demfem_motion_table  0
    global KPriv
    set rootid $AppId

    set basexpath "$rootid//c.DEM-Results//c.DEM-Graphs//c.DEM-ForceIntegrationGroup"
    set FIGgrlist [::xmlutils::setXmlContainerIds $basexpath]
    set basexpath "$AppId//c.DEM-MaterialTest//c.DEM-TopLayerGroup"
    set TOPgrlist [::xmlutils::setXmlContainerIds $basexpath]
    set basexpath "$AppId//c.DEM-MaterialTest//c.DEM-BottomLayerGroup"
    set BOTgrlist [::xmlutils::setXmlContainerIds $basexpath]

    set basexpath "$rootid//c.DEM-Conditions//c.DEM-FEM-Wall"
    set gproplist [::xmlutils::setXmlContainerIds $basexpath]

    foreach cgroupid $gproplist {
	set cproperty "dv"
	set cxpath "$rootid//c.DEM-Conditions//c.DEM-FEM-Wall//c.[list ${cgroupid}]//i.SetActive"
	set active_or_not [::xmlutils::setXml $cxpath $cproperty]
	if {$active_or_not=="No"} {
	    continue
	}
	# Get the group node list
	lassign [::wkcf::GetDemFemWallGroupNodes $cgroupid] fail msg nlist
	if {$fail == 1} {
	    return [list 1 $msg]
	}
	#set nlist [::wkcf::GetDemFemWallGroupNodes $AppId $cgroupid]

	if {[llength $nlist]} {

	    set cproperty "dv"
	    set cxpath "${basexpath}//c.[list ${cgroupid}]//c.LinearVelocity//i.LinearVelocityX"
	    set LinearVelocityX [::xmlutils::setXml $cxpath $cproperty]
	    set cxpath "${basexpath}//c.[list ${cgroupid}]//c.LinearVelocity//i.LinearVelocityY"
	    set LinearVelocityY [::xmlutils::setXml $cxpath $cproperty]
	    set cxpath "${basexpath}//c.[list ${cgroupid}]//c.LinearVelocity//i.LinearVelocityZ"
	    set LinearVelocityZ [::xmlutils::setXml $cxpath $cproperty]
	    set cxpath "${basexpath}//c.[list ${cgroupid}]//c.LinearVelocity//i.LinearPeriodic"
	    set IsPeriodic [::xmlutils::setXml $cxpath $cproperty]
	    if {$IsPeriodic=="Yes"} {
		set cxpath "${basexpath}//c.[list ${cgroupid}]//c.LinearVelocity//i.LinearPeriod"
		set Period [::xmlutils::setXml $cxpath $cproperty]
	    } else {
		set Period "0.0"
	    }
	    set cxpath "${basexpath}//c.[list ${cgroupid}]//c.LinearVelocity//i.LinearStartTime"
	    set LinearVelocityStartTime [::xmlutils::setXml $cxpath $cproperty]
	    set cxpath "${basexpath}//c.[list ${cgroupid}]//c.LinearVelocity//i.LinearEndTime"
	    set LinearVelocityEndTime [::xmlutils::setXml $cxpath $cproperty]

	    set cxpath "${basexpath}//c.[list ${cgroupid}]//c.AngularVelocity//i.AngularVelocityX"
	    set AngularVelocityX [::xmlutils::setXml $cxpath $cproperty]
	    set cxpath "${basexpath}//c.[list ${cgroupid}]//c.AngularVelocity//i.AngularVelocityY"
	    set AngularVelocityY [::xmlutils::setXml $cxpath $cproperty]
	    set cxpath "${basexpath}//c.[list ${cgroupid}]//c.AngularVelocity//i.AngularVelocityZ"
	    set AngularVelocityZ [::xmlutils::setXml $cxpath $cproperty]
	    set cxpath "${basexpath}//c.[list ${cgroupid}]//c.AngularVelocity//i.CenterOfRotationX"
	    set CenterOfRotationX [::xmlutils::setXml $cxpath $cproperty]
	    set cxpath "${basexpath}//c.[list ${cgroupid}]//c.AngularVelocity//i.CenterOfRotationY"
	    set CenterOfRotationY [::xmlutils::setXml $cxpath $cproperty]
	    set cxpath "${basexpath}//c.[list ${cgroupid}]//c.AngularVelocity//i.CenterOfRotationZ"
	    set CenterOfRotationZ [::xmlutils::setXml $cxpath $cproperty]
	    set cxpath "${basexpath}//c.[list ${cgroupid}]//c.AngularVelocity//i.AngularPeriodic"
	    set AngularIsPeriodic [::xmlutils::setXml $cxpath $cproperty]
	    if {$AngularIsPeriodic=="Yes"} {
		set cxpath "${basexpath}//c.[list ${cgroupid}]//c.AngularVelocity//i.AngularPeriod"
		set AngularPeriod [::xmlutils::setXml $cxpath $cproperty]
	    } else {
		set AngularPeriod "0.0"
	    }
	    set cxpath "${basexpath}//c.[list ${cgroupid}]//c.Options//i.fixed_wall"
            set fixed_wall [::xmlutils::setXml $cxpath $cproperty]
	    if {$fixed_wall=="Yes"} {
		set fixed_wall_value 1
	    } else {
		set fixed_wall_value 0
	    }
            set ghostpath "${basexpath}//c.[list ${cgroupid}]//c.Options//i.AnalyticProps"
            set ghost_wall [::xmlutils::setXml $ghostpath $cproperty]
	    if {$ghost_wall=="Yes"} {
		set ghost_wall_value 1
	    } else {
		set ghost_wall_value 0
	    }
	    set cxpath "${basexpath}//c.[list ${cgroupid}]//c.AngularVelocity//i.AngularStartTime"
	    set AngularVelocityStartTime [::xmlutils::setXml $cxpath $cproperty]
	    set cxpath "${basexpath}//c.[list ${cgroupid}]//c.AngularVelocity//i.AngularEndTime"
	    set AngularVelocityEndTime [::xmlutils::setXml $cxpath $cproperty]
	    set RigidBodyMotionOption 1
	    set TableNumber 0
	    set TableVelocityComponent 0
	    foreach {FreeBodyMotion RigidBodyMass CentroidX CentroidY CentroidZ InertiaX InertiaY InertiaZ Buoyancy} {0 0 0.0 0.0 0.0 1.0 1.0 1.0 "No"} {}
	    set type_of_motion [::xmlutils::setXml "${basexpath}//c.[list ${cgroupid}]//i.DEM-RBImposedMotion" dv]
	    if {$type_of_motion=="None"} {
		    foreach {LinearVelocityX LinearVelocityY LinearVelocityZ AngularVelocityX AngularVelocityY AngularVelocityZ RigidBodyMotionOption} {0.0 0.0 0.0 0.0 0.0 0.0 0} {}
	    }
	if {$type_of_motion=="FreeMotion"} {
	    set RigidBodyMotionOption 0
	    set FreeBodyMotion 1
	    set RigidBodyMass [::xmlutils::setXml "${basexpath}//c.[list ${cgroupid}]//i.Mass" dv]
	    set CentroidX [::xmlutils::setXml "${basexpath}//c.[list ${cgroupid}]//c.Centroid//i.CX" dv]
	    set CentroidY [::xmlutils::setXml "${basexpath}//c.[list ${cgroupid}]//c.Centroid//i.CY" dv]
	    set CentroidZ [::xmlutils::setXml "${basexpath}//c.[list ${cgroupid}]//c.Centroid//i.CZ" dv]
	    set InertiaX [::xmlutils::setXml "${basexpath}//c.[list ${cgroupid}]//c.DEM-Inertias//i.IX" dv]
	    set InertiaY [::xmlutils::setXml "${basexpath}//c.[list ${cgroupid}]//c.DEM-Inertias//i.IY" dv]
	    set InertiaZ [::xmlutils::setXml "${basexpath}//c.[list ${cgroupid}]//c.DEM-Inertias//i.IZ" dv]
	    #Imposed velocities
		set imposed_velocity_X [::xmlutils::setXml "${basexpath}//c.[list ${cgroupid}]//c.DEM-RBE-DOFS//i.Ax" dv]
		set imposed_velocity_Y [::xmlutils::setXml "${basexpath}//c.[list ${cgroupid}]//c.DEM-RBE-DOFS//i.Ay" dv]
		set imposed_velocity_Z [::xmlutils::setXml "${basexpath}//c.[list ${cgroupid}]//c.DEM-RBE-DOFS//i.Az" dv]
		set imposed_angular_velocity_X [::xmlutils::setXml "${basexpath}//c.[list ${cgroupid}]//c.DEM-RBE-DOFS//i.Bx" dv]
		set imposed_angular_velocity_Y [::xmlutils::setXml "${basexpath}//c.[list ${cgroupid}]//c.DEM-RBE-DOFS//i.By" dv]
		set imposed_angular_velocity_Z [::xmlutils::setXml "${basexpath}//c.[list ${cgroupid}]//c.DEM-RBE-DOFS//i.Bz" dv]
		foreach {TableNumberVX TableNumberVY TableNumberVZ TableNumberAVX TableNumberAVY TableNumberAVZ} {0 0 0 0 0 0} {}
	    if {$imposed_velocity_X=="Constant"} {
			set VelocityX [::xmlutils::setXml "${basexpath}//c.[list ${cgroupid}]//c.DEM-RBE-DOFS//i.Vx" dv]
	    }
	    if {$imposed_velocity_X=="FromATable"} {
			incr demfem_motion_table
	        set TableNumberVX $demfem_motion_table
	    }
	    if {$imposed_velocity_Y=="Constant"} {
			set VelocityY [::xmlutils::setXml "${basexpath}//c.[list ${cgroupid}]//c.DEM-RBE-DOFS//i.Vy" dv]
	    }
	    if {$imposed_velocity_Y=="FromATable"} {
			incr demfem_motion_table
	        set TableNumberVY $demfem_motion_table
	    }
	    if {$imposed_velocity_Z=="Constant"} {
			set VelocityZ [::xmlutils::setXml "${basexpath}//c.[list ${cgroupid}]//c.DEM-RBE-DOFS//i.Vz" dv]
	    }
	    if {$imposed_velocity_Z=="FromATable"} {
			incr demfem_motion_table
	        set TableNumberVZ $demfem_motion_table
	    }
	    if {$imposed_angular_velocity_X=="Constant"} {
			set AngularVelocityX [::xmlutils::setXml "${basexpath}//c.[list ${cgroupid}]//c.DEM-RBE-DOFS//i.AVx" dv]
	    }
	    if {$imposed_angular_velocity_X=="FromATable"} {
			incr demfem_motion_table
	        set TableNumberAVX $demfem_motion_table
	    }
	    if {$imposed_angular_velocity_Y=="Constant"} {
			set AngularVelocityY [::xmlutils::setXml "${basexpath}//c.[list ${cgroupid}]//c.DEM-RBE-DOFS//i.AVy" dv]
	    }
	    if {$imposed_angular_velocity_Y=="FromATable"} {
			incr demfem_motion_table
	        set TableNumberAVY $demfem_motion_table
	    }
	    if {$imposed_angular_velocity_Z=="Constant"} {
			set AngularVelocityZ [::xmlutils::setXml "${basexpath}//c.[list ${cgroupid}]//c.DEM-RBE-DOFS//i.AVz" dv]
	    }
	    if {$imposed_angular_velocity_Z=="FromATable"} {
			incr demfem_motion_table
	        set TableNumberAVZ $demfem_motion_table
	    }
		set velocity_start_time [::xmlutils::setXml "${basexpath}//c.[list ${cgroupid}]//c.DEM-RBE-DOFS//i.VStart" dv]
	    set velocity_stop_time [::xmlutils::setXml "${basexpath}//c.[list ${cgroupid}]//c.DEM-RBE-DOFS//i.VEnd" dv]
	    #Initial velocities
		set initial_velocity_X [::xmlutils::setXml "${basexpath}//c.[list ${cgroupid}]//c.DEM-RBE-InitialVelocities//i.Ax" dv]
		set initial_velocity_Y [::xmlutils::setXml "${basexpath}//c.[list ${cgroupid}]//c.DEM-RBE-InitialVelocities//i.Ay" dv]
		set initial_velocity_Z [::xmlutils::setXml "${basexpath}//c.[list ${cgroupid}]//c.DEM-RBE-InitialVelocities//i.Az" dv]
		set initial_angular_velocity_X [::xmlutils::setXml "${basexpath}//c.[list ${cgroupid}]//c.DEM-RBE-InitialVelocities//i.Bx" dv]
		set initial_angular_velocity_Y [::xmlutils::setXml "${basexpath}//c.[list ${cgroupid}]//c.DEM-RBE-InitialVelocities//i.By" dv]
		set initial_angular_velocity_Z [::xmlutils::setXml "${basexpath}//c.[list ${cgroupid}]//c.DEM-RBE-InitialVelocities//i.Bz" dv]
	    if {$initial_velocity_X=="Yes"} {
		    set InitialVelocityX [::xmlutils::setXml "${basexpath}//c.[list ${cgroupid}]//c.DEM-RBE-InitialVelocities//i.Vx" dv]
	    }
	    if {$initial_velocity_Y == "Yes" } {
		    set InitialVelocityY [::xmlutils::setXml "${basexpath}//c.[list ${cgroupid}]//c.DEM-RBE-InitialVelocities//i.Vy" dv]
	    }
	    if {$initial_velocity_Z == "Yes" } {
		    set InitialVelocityZ [::xmlutils::setXml "${basexpath}//c.[list ${cgroupid}]//c.DEM-RBE-InitialVelocities//i.Vz" dv]
	    }
	    if {$initial_angular_velocity_X == "Yes" } {
		    set InitialAngularVelocityX [::xmlutils::setXml "${basexpath}//c.[list ${cgroupid}]//c.DEM-RBE-InitialVelocities//i.AVx" dv]
	    }
	    if {$initial_angular_velocity_Y == "Yes" } {
		    set InitialAngularVelocityY [::xmlutils::setXml "${basexpath}//c.[list ${cgroupid}]//c.DEM-RBE-InitialVelocities//i.AVy" dv]
	    }
	    if {$initial_angular_velocity_Z == "Yes" } {
		    set InitialAngularVelocityZ [::xmlutils::setXml "${basexpath}//c.[list ${cgroupid}]//c.DEM-RBE-InitialVelocities//i.AVz" dv]
	    }
        #External forces and moments
	    set external_force_X [::xmlutils::setXml "${basexpath}//c.[list ${cgroupid}]//c.DEM-RBImposedForces//i.DEM-RBE-ExternalForceX" dv]
	    set external_force_Y [::xmlutils::setXml "${basexpath}//c.[list ${cgroupid}]//c.DEM-RBImposedForces//i.DEM-RBE-ExternalForceY" dv]
	    set external_force_Z [::xmlutils::setXml "${basexpath}//c.[list ${cgroupid}]//c.DEM-RBImposedForces//i.DEM-RBE-ExternalForceZ" dv]
	    set external_moment_X [::xmlutils::setXml "${basexpath}//c.[list ${cgroupid}]//c.DEM-RBImposedForces//i.DEM-RBE-ExternalForceX" dv]
	    set external_moment_Y [::xmlutils::setXml "${basexpath}//c.[list ${cgroupid}]//c.DEM-RBImposedForces//i.DEM-RBE-ExternalForceY" dv]
	    set external_moment_Z [::xmlutils::setXml "${basexpath}//c.[list ${cgroupid}]//c.DEM-RBImposedForces//i.DEM-RBE-ExternalForceZ" dv]
		foreach {TableNumberFX TableNumberFY TableNumberFZ TableNumberMX TableNumberMY TableNumberMZ} {0 0 0 0 0 0} {}
	    if {$external_force_X=="None"} {
		    set ExternalFX 0.0
	    }
	    if {$external_force_Y=="None"} {
		    set ExternalFY 0.0
	    }
	    if {$external_force_Z=="None"} {
		    set ExternalFZ 0.0
	    }
	    if {$external_force_X=="Constant"} {
	        set ExternalFX [::xmlutils::setXml "${basexpath}//c.[list ${cgroupid}]//c.DEM-RBImposedForces//i.FX" dv]
		}
	    if {$external_force_Y=="Constant"} {
	        set ExternalFY [::xmlutils::setXml "${basexpath}//c.[list ${cgroupid}]//c.DEM-RBImposedForces//i.FY" dv]
		}
	    if {$external_force_Z=="Constant"} {
	        set ExternalFZ [::xmlutils::setXml "${basexpath}//c.[list ${cgroupid}]//c.DEM-RBImposedForces//i.FZ" dv]
		}
	    if {$external_force_X=="FromATable"} {
			incr demfem_motion_table
	        set TableNumberFX $demfem_motion_table
		}
	    if {$external_force_Y=="FromATable"} {
			incr demfem_motion_table
	        set TableNumberFY $demfem_motion_table
		}
	    if {$external_force_Z=="FromATable"} {
			incr demfem_motion_table
	        set TableNumberFZ $demfem_motion_table
		}
	    if {$external_moment_X=="None"} {
		    set ExternalMX 0.0
	    }
	    if {$external_moment_Y=="None"} {
		    set ExternalMY 0.0
	    }
	    if {$external_moment_Z=="None"} {
		    set ExternalMZ 0.0
	    }
	    if {$external_moment_X=="Constant"} {
	        set ExternalMX [::xmlutils::setXml "${basexpath}//c.[list ${cgroupid}]//c.DEM-RBImposedMoments//i.MX" dv]
	    }
	    if {$external_moment_Y=="Constant"} {
	        set ExternalMY [::xmlutils::setXml "${basexpath}//c.[list ${cgroupid}]//c.DEM-RBImposedMoments//i.MY" dv]
	    }
	    if {$external_moment_Z=="Constant"} {
	        set ExternalMZ [::xmlutils::setXml "${basexpath}//c.[list ${cgroupid}]//c.DEM-RBImposedMoments//i.MZ" dv]
	    }
	    if {$external_moment_X=="FromATable"} {
			incr demfem_motion_table
	        set TableNumberMX $demfem_motion_table
	    }
	    if {$external_moment_Y=="FromATable"} {
			incr demfem_motion_table
	        set TableNumberMY $demfem_motion_table
	    }
	    if {$external_moment_Z=="FromATable"} {
			incr demfem_motion_table
	        set TableNumberMZ $demfem_motion_table
	    }
	    if {$KPriv(what_dempack_package) eq "C-DEMPack"} {
		    set cxpath "${basexpath}//c.[list ${cgroupid}]//c.Options//i.ShipElement"
		    set Buoyancy [::xmlutils::setXml $cxpath dv]
		    if {$Buoyancy=="Yes"} {
		        set Buoyancy 1
		    } else {
		        set Buoyancy 0
		    }
		    set cxpath "${basexpath}//c.[list ${cgroupid}]//c.Options//i.EnginePower"
		    set enginepower [::xmlutils::setXml $cxpath dv]
		    set cxpath "${basexpath}//c.[list ${cgroupid}]//c.Options//i.MaxEngineForce"
		    set maxengineforce [::xmlutils::setXml $cxpath dv]
		    set cxpath "${basexpath}//c.[list ${cgroupid}]//c.Options//i.ThresholdVelocity"
		    set thresholdvelocity [::xmlutils::setXml $cxpath dv]
		    set cxpath "${basexpath}//c.[list ${cgroupid}]//c.Options//i.EnginePerformance"
		    set engineperformance [::xmlutils::setXml $cxpath dv]
		    set cxpath "${basexpath}//c.[list ${cgroupid}]//c.Options//i.DragConstantX"
		    set dragconstantx [::xmlutils::setXml $cxpath dv]
		    set cxpath "${basexpath}//c.[list ${cgroupid}]//c.Options//i.DragConstantY"
		    set dragconstanty [::xmlutils::setXml $cxpath dv]
		    set cxpath "${basexpath}//c.[list ${cgroupid}]//c.Options//i.DragConstantZ"
		    set dragconstantz [::xmlutils::setXml $cxpath dv]
	    }
	}
        #Imposed velocity tables
	    if {$imposed_velocity_X=="FromATable"} {
			set filename [::xmlutils::setXml "${basexpath}//c.[list ${cgroupid}]//c.DEM-RBE-DOFS//i.FilenameVx" dv]
	        GiD_File fprintf $demfemchannel "Begin Table $TableNumberVX TIME VELOCITY_X"
	        set file_open [open [file native [file join [::KUtils::GetPaths "PDir"] $filename]] r]
	        set file_data [read $file_open]
	        close $file_open
	        GiD_File fprintf -nonewline $demfemchannel $file_data
	        GiD_File fprintf $demfemchannel "End Table"
	        GiD_File fprintf $demfemchannel ""
	    }
	    if {$imposed_velocity_Y=="FromATable"} {
			set filename [::xmlutils::setXml "${basexpath}//c.[list ${cgroupid}]//c.DEM-RBE-DOFS//i.FilenameVy" dv]
	        GiD_File fprintf $demfemchannel "Begin Table $TableNumberVY TIME VELOCITY_Y"
	        set file_open [open [file native [file join [::KUtils::GetPaths "PDir"] $filename]] r]
	        set file_data [read $file_open]
	        close $file_open
	        GiD_File fprintf -nonewline $demfemchannel $file_data
	        GiD_File fprintf $demfemchannel "End Table"
	        GiD_File fprintf $demfemchannel ""
	    }
	    if {$imposed_velocity_Z=="FromATable"} {
			set filename [::xmlutils::setXml "${basexpath}//c.[list ${cgroupid}]//c.DEM-RBE-DOFS//i.FilenameVz" dv]
	        GiD_File fprintf $demfemchannel "Begin Table $TableNumberVZ TIME VELOCITY_Z"
	        set file_open [open [file native [file join [::KUtils::GetPaths "PDir"] $filename]] r]
	        set file_data [read $file_open]
	        close $file_open
	        GiD_File fprintf -nonewline $demfemchannel $file_data
	        GiD_File fprintf $demfemchannel "End Table"
	        GiD_File fprintf $demfemchannel ""
	    }
	    if {$imposed_angular_velocity_X=="FromATable"} {
			set filename [::xmlutils::setXml "${basexpath}//c.[list ${cgroupid}]//c.DEM-RBE-DOFS//i.FilenameAVx" dv]
	        GiD_File fprintf $demfemchannel "Begin Table $TableNumberAVX TIME ANGULAR_VELOCITY_X"
	        set file_open [open [file native [file join [::KUtils::GetPaths "PDir"] $filename]] r]
	        set file_data [read $file_open]
	        close $file_open
	        GiD_File fprintf -nonewline $demfemchannel $file_data
	        GiD_File fprintf $demfemchannel "End Table"
	        GiD_File fprintf $demfemchannel ""
	    }
	    if {$imposed_angular_velocity_Y=="FromATable"} {
			set filename [::xmlutils::setXml "${basexpath}//c.[list ${cgroupid}]//c.DEM-RBE-DOFS//i.FilenameAVy" dv]
	        GiD_File fprintf $demfemchannel "Begin Table $TableNumberAVY TIME ANGULAR_VELOCITY_Y"
	        set file_open [open [file native [file join [::KUtils::GetPaths "PDir"] $filename]] r]
	        set file_data [read $file_open]
	        close $file_open
	        GiD_File fprintf -nonewline $demfemchannel $file_data
	        GiD_File fprintf $demfemchannel "End Table"
	        GiD_File fprintf $demfemchannel ""
	    }
	    if {$imposed_angular_velocity_Z=="FromATable"} {
			set filename [::xmlutils::setXml "${basexpath}//c.[list ${cgroupid}]//c.DEM-RBE-DOFS//i.FilenameAVz" dv]
	        GiD_File fprintf $demfemchannel "Begin Table $TableNumberAVZ TIME ANGULAR_VELOCITY_Z"
	        set file_open [open [file native [file join [::KUtils::GetPaths "PDir"] $filename]] r]
	        set file_data [read $file_open]
	        close $file_open
	        GiD_File fprintf -nonewline $demfemchannel $file_data
	        GiD_File fprintf $demfemchannel "End Table"
	        GiD_File fprintf $demfemchannel ""
	    }
        #Imposed external applied force tables
	    if {$external_force_X=="FromATable"} {
			set filename [::xmlutils::setXml "${basexpath}//c.[list ${cgroupid}]//c.DEM-RBImposedForces//i.FilenameFX" dv]
	        GiD_File fprintf $demfemchannel "Begin Table $TableNumberFX TIME EXTERNAL_APPLIED_FORCE_X"
	        set file_open [open [file native [file join [::KUtils::GetPaths "PDir"] $filename]] r]
	        set file_data [read $file_open]
	        close $file_open
	        GiD_File fprintf -nonewline $demfemchannel $file_data
	        GiD_File fprintf $demfemchannel "End Table"
	        GiD_File fprintf $demfemchannel ""
		}
	    if {$external_force_Y=="FromATable"} {
			set filename [::xmlutils::setXml "${basexpath}//c.[list ${cgroupid}]//c.DEM-RBImposedForces//i.FilenameFY" dv]
	        GiD_File fprintf $demfemchannel "Begin Table $TableNumberFY TIME EXTERNAL_APPLIED_FORCE_Y"
	        set file_open [open [file native [file join [::KUtils::GetPaths "PDir"] $filename]] r]
	        set file_data [read $file_open]
	        close $file_open
	        GiD_File fprintf -nonewline $demfemchannel $file_data
	        GiD_File fprintf $demfemchannel "End Table"
	        GiD_File fprintf $demfemchannel ""
		}
	    if {$external_force_Z=="FromATable"} {
			set filename [::xmlutils::setXml "${basexpath}//c.[list ${cgroupid}]//c.DEM-RBImposedForces//i.FilenameFZ" dv]
	        GiD_File fprintf $demfemchannel "Begin Table $TableNumberFZ TIME EXTERNAL_APPLIED_FORCE_Z"
	        set file_open [open [file native [file join [::KUtils::GetPaths "PDir"] $filename]] r]
	        set file_data [read $file_open]
	        close $file_open
	        GiD_File fprintf -nonewline $demfemchannel $file_data
	        GiD_File fprintf $demfemchannel "End Table"
	        GiD_File fprintf $demfemchannel ""
	    }
	    if {$external_moment_X=="FromATable"} {
			set filename [::xmlutils::setXml "${basexpath}//c.[list ${cgroupid}]//c.DEM-RBImposedMoments//i.FilenameMX" dv]
	        GiD_File fprintf $demfemchannel "Begin Table $TableNumberMX TIME EXTERNAL_APPLIED_MOMENT_X"
	        set file_open [open [file native [file join [::KUtils::GetPaths "PDir"] $filename]] r]
	        set file_data [read $file_open]
	        close $file_open
	        GiD_File fprintf -nonewline $demfemchannel $file_data
	        GiD_File fprintf $demfemchannel "End Table"
	        GiD_File fprintf $demfemchannel ""
		}
	    if {$external_moment_Y=="FromATable"} {
			set filename [::xmlutils::setXml "${basexpath}//c.[list ${cgroupid}]//c.DEM-RBImposedMoments//i.FilenameMY" dv]
	        GiD_File fprintf $demfemchannel "Begin Table $TableNumberMY TIME EXTERNAL_APPLIED_MOMENT_Y"
	        set file_open [open [file native [file join [::KUtils::GetPaths "PDir"] $filename]] r]
	        set file_data [read $file_open]
	        close $file_open
	        GiD_File fprintf -nonewline $demfemchannel $file_data
	        GiD_File fprintf $demfemchannel "End Table"
	        GiD_File fprintf $demfemchannel ""
		}
	    if {$external_moment_Z=="FromATable"} {
			set filename [::xmlutils::setXml "${basexpath}//c.[list ${cgroupid}]//c.DEM-RBImposedMoments//i.FilenameMZ" dv]
	        GiD_File fprintf $demfemchannel "Begin Table $TableNumberMZ TIME EXTERNAL_APPLIED_MOMENT_Z"
	        set file_open [open [file native [file join [::KUtils::GetPaths "PDir"] $filename]] r]
	        set file_data [read $file_open]
	        close $file_open
	        GiD_File fprintf -nonewline $demfemchannel $file_data
	        GiD_File fprintf $demfemchannel "End Table"
	        GiD_File fprintf $demfemchannel ""
	    }
	    # Write mesh properties for this group
	    GiD_File fprintf $demfemchannel "%s" "Begin SubModelPart $demfem_ref_to_props_number \/\/ DEM-FEM-Wall. Group name: $cgroupid"
	    GiD_File fprintf $demfemchannel "  Begin SubModelPartData // DEM-FEM-Wall. Group name: $cgroupid"
	    GiD_File fprintf $demfemchannel "  PROPERTIES_ID $demfem_ref_to_props_number"
	    GiD_File fprintf $demfemchannel "  RIGID_BODY_MOTION $RigidBodyMotionOption"
	    GiD_File fprintf $demfemchannel "  FREE_BODY_MOTION $FreeBodyMotion"
	    if {$type_of_motion=="LinearPeriodic"} {
	    GiD_File fprintf $demfemchannel "  FIXED_MESH_OPTION $fixed_wall_value"
	    GiD_File fprintf $demfemchannel "  LINEAR_VELOCITY \[3\] ($LinearVelocityX,$LinearVelocityY,$LinearVelocityZ)"
	    GiD_File fprintf $demfemchannel "  VELOCITY_PERIOD $Period"
	    GiD_File fprintf $demfemchannel "  ANGULAR_VELOCITY \[3\] ($AngularVelocityX,$AngularVelocityY,$AngularVelocityZ)"
	    GiD_File fprintf $demfemchannel "  ROTATION_CENTER \[3\] ($CenterOfRotationX,$CenterOfRotationY,$CenterOfRotationZ)"
	    GiD_File fprintf $demfemchannel "  ANGULAR_VELOCITY_PERIOD $AngularPeriod"
	    GiD_File fprintf $demfemchannel "  VELOCITY_START_TIME $LinearVelocityStartTime"
	    GiD_File fprintf $demfemchannel "  VELOCITY_STOP_TIME $LinearVelocityEndTime"
	    GiD_File fprintf $demfemchannel "  ANGULAR_VELOCITY_START_TIME $AngularVelocityStartTime"
	    GiD_File fprintf $demfemchannel "  ANGULAR_VELOCITY_STOP_TIME $AngularVelocityEndTime"
	    }
	    if {$type_of_motion=="FreeMotion"} {
	    GiD_File fprintf $demfemchannel "  RIGID_BODY_MASS $RigidBodyMass"
	    GiD_File fprintf $demfemchannel "  RIGID_BODY_CENTER_OF_MASS \[3\] ($CentroidX,$CentroidY,$CentroidZ)"
	    GiD_File fprintf $demfemchannel "  RIGID_BODY_INERTIAS \[3\] ($InertiaX,$InertiaY,$InertiaZ)"
	    #Imposed velocities
	    if {$imposed_velocity_X=="Constant"} {
			GiD_File fprintf $demfemchannel "  IMPOSED_VELOCITY_X_VALUE $VelocityX"
	    }
	    if {$imposed_velocity_Y=="Constant"} {
			GiD_File fprintf $demfemchannel "  IMPOSED_VELOCITY_Y_VALUE $VelocityY"
	    }
	    if {$imposed_velocity_Z=="Constant"} {
			GiD_File fprintf $demfemchannel "  IMPOSED_VELOCITY_Z_VALUE $VelocityZ"
	    }
		GiD_File fprintf $demfemchannel "  TABLE_NUMBER_VELOCITY \[3\] ($TableNumberVX,$TableNumberVY,$TableNumberVZ)"
	    if {$imposed_angular_velocity_X=="Constant"} {
			GiD_File fprintf $demfemchannel "  IMPOSED_ANGULAR_VELOCITY_X_VALUE $AngularVelocityX"
	    }
	    if {$imposed_angular_velocity_Y=="Constant"} {
			GiD_File fprintf $demfemchannel "  IMPOSED_ANGULAR_VELOCITY_Y_VALUE $AngularVelocityY"
	    }
	    if {$imposed_angular_velocity_Z=="Constant"} {
			GiD_File fprintf $demfemchannel "  IMPOSED_ANGULAR_VELOCITY_Z_VALUE $AngularVelocityZ"
	    }
	    GiD_File fprintf $demfemchannel "  TABLE_NUMBER_ANGULAR_VELOCITY \[3\] ($TableNumberAVX,$TableNumberAVY,$TableNumberAVZ)"
	    GiD_File fprintf $demfemchannel "  VELOCITY_START_TIME $velocity_start_time"
	    GiD_File fprintf $demfemchannel "  VELOCITY_STOP_TIME $velocity_stop_time"
	    #Initial velocities
	    if {$initial_velocity_X=="Yes"} {
		    GiD_File fprintf $demfemchannel "  INITIAL_VELOCITY_X_VALUE $InitialVelocityX"
	    }
	    if {$initial_velocity_Y=="Yes"} {
    		GiD_File fprintf $demfemchannel "  INITIAL_VELOCITY_Y_VALUE $InitialVelocityY"
	    }
	    if {$initial_velocity_Z=="Yes"} {
	    	GiD_File fprintf $demfemchannel "  INITIAL_VELOCITY_Z_VALUE $InitialVelocityZ"
	    }
	    if {$initial_angular_velocity_X=="Yes"} {
		    GiD_File fprintf $demfemchannel "  INITIAL_ANGULAR_VELOCITY_X_VALUE $InitialAngularVelocityX"
	    }
	    if {$initial_angular_velocity_Y=="Yes"} {
    		GiD_File fprintf $demfemchannel "  INITIAL_ANGULAR_VELOCITY_Y_VALUE $InitialAngularVelocityY"
	    }
	    if {$initial_angular_velocity_Z=="Yes"} {
	    	GiD_File fprintf $demfemchannel "  INITIAL_ANGULAR_VELOCITY_Z_VALUE $InitialAngularVelocityZ"
	    }
        #External forces
	    if {$external_force_X=="Constant"} {
    	    GiD_File fprintf $demfemchannel "  EXTERNAL_APPLIED_FORCE_X $ExternalFX"
		}
        if {$external_force_Y=="Constant"} {
	        GiD_File fprintf $demfemchannel "  EXTERNAL_APPLIED_FORCE_Y $ExternalFY"
		}
        if {$external_force_Z=="Constant"} {
	        GiD_File fprintf $demfemchannel "  EXTERNAL_APPLIED_FORCE_Z $ExternalFZ"
		}
	    GiD_File fprintf $demfemchannel "  TABLE_NUMBER_FORCE \[3\] ($TableNumberFX,$TableNumberFY,$TableNumberFZ)"
	    if {$external_moment_X=="Constant"} {
    	    GiD_File fprintf $demfemchannel "  EXTERNAL_APPLIED_MOMENT_X $ExternalMX"
		}
	    if {$external_moment_Y=="Constant"} {
    	    GiD_File fprintf $demfemchannel "  EXTERNAL_APPLIED_MOMENT_Y $ExternalMY"
		}
	    if {$external_moment_Z=="Constant"} {
    	    GiD_File fprintf $demfemchannel "  EXTERNAL_APPLIED_MOMENT_Z $ExternalMZ"
		}
	    GiD_File fprintf $demfemchannel "  TABLE_NUMBER_MOMENT \[3\] ($TableNumberMX,$TableNumberMY,$TableNumberMZ)"
	    if {$KPriv(what_dempack_package) eq "C-DEMPack"} {
		GiD_File fprintf $demfemchannel "  FLOATING_OPTION $Buoyancy"
		GiD_File fprintf $demfemchannel "  DEM_ENGINE_POWER $enginepower"
		GiD_File fprintf $demfemchannel "  DEM_MAX_ENGINE_FORCE $maxengineforce"
		GiD_File fprintf $demfemchannel "  DEM_THRESHOLD_VELOCITY $thresholdvelocity"
		GiD_File fprintf $demfemchannel "  DEM_ENGINE_PERFORMANCE $engineperformance"
		GiD_File fprintf $demfemchannel "  DEM_DRAG_CONSTANT_X $dragconstantx"
		GiD_File fprintf $demfemchannel "  DEM_DRAG_CONSTANT_Y $dragconstanty"
		GiD_File fprintf $demfemchannel "  DEM_DRAG_CONSTANT_Z $dragconstantz"
	    }
	    }
        GiD_File fprintf $demfemchannel "  IS_GHOST $ghost_wall_value"
	    GiD_File fprintf $demfemchannel "  IDENTIFIER $cgroupid"

	    set TOP 0
	    if {$cgroupid in $TOPgrlist} {set TOP 1}
	    GiD_File fprintf $demfemchannel "  [::xmlutils::getKKWord "Applications/$AppId" "TopGroup"] $TOP"
	    set BOT 0
	    if {$cgroupid in $BOTgrlist} {set BOT 1}
	    GiD_File fprintf $demfemchannel "  [::xmlutils::getKKWord "Applications/$AppId" "BottomGroup"] $BOT"
	    set FIG 0
	    if {$cgroupid in $FIGgrlist} {set FIG 1}
	    if {$BOT} {set FIG 1}
	    if {$TOP} {set FIG 1}
	    GiD_File fprintf $demfemchannel "  [::xmlutils::getKKWord "Applications/$AppId" "ForceIntegrationGroup"] $FIG"

	    GiD_File fprintf $demfemchannel "  End SubModelPartData"
	    GiD_File fprintf $demfemchannel "  Begin SubModelPartNodes"

	    foreach nodeid $nlist {
		GiD_File fprintf $demfemchannel "  $nodeid"
	    }
	    GiD_File fprintf $demfemchannel "  End SubModelPartNodes"
	    #
	    set nlist [GiD_EntitiesGroups get $cgroupid elements]
        if {[llength $nlist]} {
            GiD_File fprintf $demfemchannel "  Begin SubModelPartConditions"
            foreach elemid $nlist {
                GiD_File fprintf $demfemchannel "  $elemid"
            }
            GiD_File fprintf $demfemchannel "  End SubModelPartConditions"
        }
        GiD_File fprintf $demfemchannel "  Begin SubModelPartTables"
        if {$TableNumberVX>0} {
			GiD_File fprintf $demfemchannel "  $TableNumberVX"
		}
        if {$TableNumberVY>0} {
			GiD_File fprintf $demfemchannel "  $TableNumberVY"
		}
        if {$TableNumberVZ>0} {
			GiD_File fprintf $demfemchannel "  $TableNumberVZ"
		}
        if {$TableNumberAVX>0} {
			GiD_File fprintf $demfemchannel "  $TableNumberAVX"
		}
        if {$TableNumberAVY>0} {
			GiD_File fprintf $demfemchannel "  $TableNumberAVY"
		}
        if {$TableNumberAVZ>0} {
			GiD_File fprintf $demfemchannel "  $TableNumberAVZ"
		}
        if {$TableNumberFX>0} {
			GiD_File fprintf $demfemchannel "  $TableNumberFX"
		}
        if {$TableNumberFY>0} {
			GiD_File fprintf $demfemchannel "  $TableNumberFY"
		}
        if {$TableNumberFZ>0} {
			GiD_File fprintf $demfemchannel "  $TableNumberFZ"
		}
        if {$TableNumberMX>0} {
			GiD_File fprintf $demfemchannel "  $TableNumberMX"
		}
        if {$TableNumberMY>0} {
			GiD_File fprintf $demfemchannel "  $TableNumberMY"
		}
        if {$TableNumberMZ>0} {
			GiD_File fprintf $demfemchannel "  $TableNumberMZ"
		}
        GiD_File fprintf $demfemchannel "  End SubModelPartTables"
	    #
	    GiD_File fprintf $demfemchannel "End SubModelPart"
	    GiD_File fprintf $demfemchannel ""
	    incr demfem_ref_to_props_number
	}
    }
}


proc ::wkcf::WriteCustomSubModelParts {AppId} {
    variable deminletchannel
	variable filechannel
	variable demfemchannel
    global KPriv

    # Set the rootid
    set rootid "$AppId"
	set basexpath "$rootid//c.DEM-CustomSubModelParts"

    set gproplist [::xmlutils::setXmlContainerIds $basexpath]

    foreach cgroupid $gproplist {

		set properties_path "${basexpath}//c.[list ${cgroupid}]//c.MainProperties"
		set cproperty "dv"
		set destination_mdpa [::xmlutils::setXml "${properties_path}//i.WhatMdpa" $cproperty]

		if { $destination_mdpa=="DEM" } {
			set output_channel $filechannel
		} elseif { $destination_mdpa == "FEM" } {
			set output_channel $demfemchannel
		} elseif { $destination_mdpa == "DEM-Inlet" } {
			set output_channel $deminletchannel
		}

		GiD_File fprintf $output_channel "%s" "Begin SubModelPart $cgroupid \/\/ Custom SubModelPart. Group name: $cgroupid"
		GiD_File fprintf $output_channel "  Begin SubModelPartData // DEM-FEM-Wall. Group name: $cgroupid"
		GiD_File fprintf $output_channel "  End SubModelPartData"
		GiD_File fprintf $output_channel "  Begin SubModelPartNodes"

		set nlist [GiD_EntitiesGroups get $cgroupid nodes]
		foreach nodeid $nlist {
			GiD_File fprintf $output_channel "  $nodeid"
		}
		GiD_File fprintf $output_channel "  End SubModelPartNodes"

		if { $destination_mdpa != "DEM-Inlet" } {
			set nlist [GiD_EntitiesGroups get $cgroupid elements]
			if {[llength $nlist]} {
				if { $destination_mdpa == "FEM" } {
					GiD_File fprintf $output_channel "  Begin SubModelPartConditions"
				} else {
					GiD_File fprintf $output_channel "  Begin SubModelPartElements"
				}
				foreach elemid $nlist {
					GiD_File fprintf $output_channel "  $elemid"
				}
				if { $destination_mdpa == "FEM" } {
					GiD_File fprintf $output_channel "  End SubModelPartConditions"
				} else {
					GiD_File fprintf $output_channel "  End SubModelPartElements"
				}
			}
		}

		GiD_File fprintf $output_channel "End SubModelPart"
		GiD_File fprintf $output_channel ""
	}

}


proc ::wkcf::WriteInletGroupMeshProperties {AppId} {
    # ABSTRACT: Write inlet condition group properties mdpa file (only the nodes)
    variable dprops
    variable deminletchannel
    variable dem_ref_to_props_number
    global KPriv
    variable ActiveAppList
    variable ndime

    # For debug
    if {!$::wkcf::pflag} {
	set inittime [clock seconds]
    }

    # Set the rootid
    set rootid "$AppId"

    # Get the values
    set basexpath "$rootid//c.DEM-Conditions//c.DEM-Inlet"
    set gproplist [::xmlutils::setXmlContainerIds $basexpath]

    foreach cgroupid $gproplist {

	set properties_path "${basexpath}//c.[list ${cgroupid}]//c.MainProperties"
	set cproperty "dv"
	set active_or_not [::xmlutils::setXml "${properties_path}//i.SetActive" $cproperty]
	if {$active_or_not=="No"} {
	    continue
	}

    set TableNumber 0
    set TableVelocityComponent 0
    set type_of_motion [::xmlutils::setXml "${properties_path}//i.DEM-RBImposedMotion" dv]
    if {$type_of_motion=="FromATable"} {
	set TableNumber $dem_ref_to_props_number
	set TableVelocityComponent [::xmlutils::setXml "${properties_path}//i.TableVelocityComponent" dv]
    }

	set cproperty "dv"
	set LinearVelocityX [::xmlutils::setXml "${properties_path}//c.LinearVelocity//i.LinearVelocityX" $cproperty]
	set LinearVelocityY [::xmlutils::setXml "${properties_path}//c.LinearVelocity//i.LinearVelocityY" $cproperty]
	set LinearVelocityZ [::xmlutils::setXml "${properties_path}//c.LinearVelocity//i.LinearVelocityZ" $cproperty]
	set IsPeriodic [::xmlutils::setXml "${properties_path}//c.LinearVelocity//i.LinearPeriodic" $cproperty]
	if {$IsPeriodic=="Yes"} {
	    set Period [::xmlutils::setXml "${basexpath}//c.[list ${cgroupid}]//c.MainProperties//c.LinearVelocity//i.LinearPeriod" $cproperty]
	} else {
	    set Period "0.0"
	}
	set LinearVelocityStartTime [::xmlutils::setXml "${properties_path}//c.LinearVelocity//i.LinearStartTime" $cproperty]
	set LinearVelocityEndTime [::xmlutils::setXml "${properties_path}//c.LinearVelocity//i.LinearEndTime" $cproperty]

	set AngularVelocityX  [::xmlutils::setXml "${properties_path}//c.AngularVelocity//i.AngularVelocityX" $cproperty]
	set AngularVelocityY  [::xmlutils::setXml "${properties_path}//c.AngularVelocity//i.AngularVelocityY" $cproperty]
	set AngularVelocityZ  [::xmlutils::setXml "${properties_path}//c.AngularVelocity//i.AngularVelocityZ" $cproperty]
	set CenterOfRotationX [::xmlutils::setXml "${properties_path}//c.AngularVelocity//i.CenterOfRotationX" $cproperty]
	set CenterOfRotationY [::xmlutils::setXml "${properties_path}//c.AngularVelocity//i.CenterOfRotationY" $cproperty]
	set CenterOfRotationZ [::xmlutils::setXml "${properties_path}//c.AngularVelocity//i.CenterOfRotationZ" $cproperty]
	set AngularIsPeriodic [::xmlutils::setXml "${properties_path}//c.AngularVelocity//i.AngularPeriodic" $cproperty]
	if {$AngularIsPeriodic=="Yes"} {
	    set AngularPeriod [::xmlutils::setXml "${properties_path}//c.AngularVelocity//i.AngularPeriod" $cproperty]
	} else {
	    set AngularPeriod "0.0"
	}
	set AngularVelocityStartTime [::xmlutils::setXml "${properties_path}//c.AngularVelocity//i.AngularStartTime" $cproperty]
	set AngularVelocityEndTime [::xmlutils::setXml "${properties_path}//c.AngularVelocity//i.AngularEndTime" $cproperty]

	# Get the group node list
	set nlist [::wkcf::GetInletGroupNodes $AppId $cgroupid]
	if {[llength $nlist]} {
	    # Write mesh properties for this group
	    GiD_File fprintf $deminletchannel "%s" "Begin SubModelPart $dem_ref_to_props_number \/\/ Group name: $cgroupid"

	    GiD_File fprintf $deminletchannel "  Begin SubModelPartData // DEM-Inlet. Group name: $cgroupid"
	    GiD_File fprintf $deminletchannel "  PROPERTIES_ID $dem_ref_to_props_number"
	    if {$type_of_motion=="LinearPeriodic"} {
		GiD_File fprintf $deminletchannel "  LINEAR_VELOCITY \[3\] ($LinearVelocityX,$LinearVelocityY,$LinearVelocityZ)"
		GiD_File fprintf $deminletchannel "  VELOCITY_PERIOD $Period"
		GiD_File fprintf $deminletchannel "  ANGULAR_VELOCITY \[3\] ($AngularVelocityX,$AngularVelocityY,$AngularVelocityZ)"
		GiD_File fprintf $deminletchannel "  ROTATION_CENTER \[3\] ($CenterOfRotationX,$CenterOfRotationY,$CenterOfRotationZ)"
		GiD_File fprintf $deminletchannel "  ANGULAR_VELOCITY_PERIOD $AngularPeriod"
		GiD_File fprintf $deminletchannel "  VELOCITY_START_TIME $LinearVelocityStartTime"
		GiD_File fprintf $deminletchannel "  VELOCITY_STOP_TIME $LinearVelocityEndTime"
		GiD_File fprintf $deminletchannel "  ANGULAR_VELOCITY_START_TIME $AngularVelocityStartTime"
		GiD_File fprintf $deminletchannel "  ANGULAR_VELOCITY_STOP_TIME $AngularVelocityEndTime"
		GiD_File fprintf $deminletchannel "  RIGID_BODY_MOTION 1"
	    } else {
		GiD_File fprintf $deminletchannel "  RIGID_BODY_MOTION 0"
		GiD_File fprintf $deminletchannel "  //TABLE_VELOCITY_COMPONENT $TableVelocityComponent"
		}
	    GiD_File fprintf $deminletchannel "  IDENTIFIER $cgroupid"
	    set material [::xmlutils::setXml "${properties_path}//i.Material" "dv"]
	    set constitutive_law ""
	    set contains_clusters 0
	    set random_orientation 0
	    set using_dem_kdem 0
	    set active_or_not [::xmlutils::setXml "${properties_path}//i.SetActive" "dv"]

	    if {$KPriv(what_dempack_package) eq "G-DEMPack"} {
	    if {"Fluid" in $ActiveAppList} {
		set inlet_element_type SphericSwimmingParticle3D
	    } else {
		set inlet_element_type SphericParticle3D
	    }
	    set inlet_injector_element_type SphericParticle3D
	} elseif {$KPriv(what_dempack_package) eq "F-DEMPack"} {
		        set inlet_element_type SphericSwimmingParticle3D
		    set inlet_injector_element_type SphericParticle3D
		} else {
	    set inlet_element_type SphericContinuumParticle3D
	    set inlet_injector_element_type SphericContinuumParticle3D
	}
	if {"Fluid" ni $ActiveAppList} {
	    if {[::xmlutils::setXml "${properties_path}//i.InletElementType" "dv"] eq "Cluster3D"} {
		set inlet_element_type [::xmlutils::setXml "${properties_path}//i.Cluster3D" dv]
		set contains_clusters 1
		lassign [::wkcf::GetClusterFileNameAndReplaceInletElementType $inlet_element_type] inlet_element_type cluster_file_name
	    }
	}
	    if {$inlet_element_type eq "Cluster3D"} {
		lappend KPriv(list_of_cluster_files) $cluster_file_name
		GiD_File fprintf $deminletchannel "  CLUSTER_FILE_NAME $cluster_file_name"
	    }
	    GiD_File fprintf $deminletchannel "  INJECTOR_ELEMENT_TYPE $inlet_injector_element_type"
	    GiD_File fprintf $deminletchannel "  ELEMENT_TYPE $inlet_element_type"
	    GiD_File fprintf $deminletchannel "  CONTAINS_CLUSTERS $contains_clusters"
	    set velocity_modulus [::xmlutils::setXml "${properties_path}//i.VelocityModulus" "dv"]
	    set velocity_x [expr {$velocity_modulus * [::xmlutils::setXml "${properties_path}//i.DirectionVectorX" "dv"]}]
	    set velocity_y [expr {$velocity_modulus * [::xmlutils::setXml "${properties_path}//i.DirectionVectorY" "dv"]}]
	    set velocity_z [expr {$velocity_modulus * [::xmlutils::setXml "${properties_path}//i.DirectionVectorZ" "dv"]}]
	    GiD_File fprintf $deminletchannel "  VELOCITY \[3\] ($velocity_x, $velocity_y, $velocity_z)"
	    set max_deviation_angle [::xmlutils::setXml "${properties_path}//i.VelocityDeviation" "dv"]
	    GiD_File fprintf $deminletchannel "  MAX_RAND_DEVIATION_ANGLE $max_deviation_angle"

	    if {"Fluid" ni $ActiveAppList} {
	    if {[::xmlutils::setXml "${properties_path}//i.Cluster3D" dv] eq "SingleSphereCluster3D"} {
		#set excentricity [::xmlutils::setXml "${properties_path}//i.Excentricity" dv]
		GiD_File fprintf $deminletchannel "  EXCENTRICITY [::xmlutils::setXml "${properties_path}//i.Excentricity" dv]"
		#set probability_distribution [::xmlutils::setXml "${properties_path}//i.ProbabilityDistributionOfExcentricity" dv]
		GiD_File fprintf $deminletchannel "  EXCENTRICITY_PROBABILITY_DISTRIBUTION [::xmlutils::setXml "${properties_path}//i.ProbabilityDistributionOfExcentricity" dv]"
		#set standard_deviation [::xmlutils::setXml "${properties_path}//i.StandardDeviationOfExcentricity" dv]
		GiD_File fprintf $deminletchannel "  EXCENTRICITY_STANDARD_DEVIATION [::xmlutils::setXml "${properties_path}//i.StandardDeviationOfExcentricity" dv]"
	    }
	    }

	    set inlet_number_of_particles [::xmlutils::setXml "${properties_path}//i.NumberOfParticles" "dv"]
	    GiD_File fprintf $deminletchannel "  INLET_NUMBER_OF_PARTICLES $inlet_number_of_particles"

	    set type_of_measurement [::xmlutils::setXml "${properties_path}//i.TypeOfFlowMeasurement" "dv"]
	    if {$type_of_measurement eq "mass_flow"} {
	    set mass_flow_option 1
	} else {
	    set mass_flow_option 0
	}
	    GiD_File fprintf $deminletchannel "  IMPOSED_MASS_FLOW_OPTION $mass_flow_option"

	    set inlet_mass_flow [::xmlutils::setXml "${properties_path}//i.InletMassFlow" "dv"]
	    GiD_File fprintf $deminletchannel "  MASS_FLOW $inlet_mass_flow"
	    set inlet_start_time [::xmlutils::setXml "${properties_path}//i.InletStartTime" "dv"]
	    GiD_File fprintf $deminletchannel "  INLET_START_TIME $inlet_start_time"
	    set inlet_stop_time [::xmlutils::setXml "${properties_path}//i.InletStopTime" "dv"]
	    GiD_File fprintf $deminletchannel "  INLET_STOP_TIME $inlet_stop_time"
	    set particle_diameter [::xmlutils::setXml "${properties_path}//i.ParticleDiameter" "dv"]
	    GiD_File fprintf $deminletchannel "  RADIUS [expr {0.5 * $particle_diameter}]"
	    set probability_distribution [::xmlutils::setXml "${properties_path}//i.ProbabilityDistribution" "dv"]
	    GiD_File fprintf $deminletchannel "  PROBABILITY_DISTRIBUTION $probability_distribution"
	    set standard_deviation [::xmlutils::setXml "${properties_path}//i.StandardDeviation" "dv"]
	    GiD_File fprintf $deminletchannel "  STANDARD_DEVIATION $standard_deviation"
	    if {[::xmlutils::setXml "${properties_path}//i.RandomOrientation" "dv"] == "Yes"} {
	    set random_orientation 1
	    }

	    GiD_File fprintf $deminletchannel "  RANDOM_ORIENTATION $random_orientation"
	    set orientation_x [::xmlutils::setXml "${properties_path}//i.OrientationX" "dv"]
	    set orientation_y [::xmlutils::setXml "${properties_path}//i.OrientationY" "dv"]
	    set orientation_z [::xmlutils::setXml "${properties_path}//i.OrientationZ" "dv"]
	    set orientation_w [::xmlutils::setXml "${properties_path}//i.OrientationW" "dv"]
	    GiD_File fprintf $deminletchannel "  ORIENTATION \[4\] ($orientation_x, $orientation_y, $orientation_z, $orientation_w)"

	    GiD_File fprintf $deminletchannel "  End SubModelPartData"

	    # Write nodes
	    GiD_File fprintf $deminletchannel "%s" "  Begin SubModelPartNodes"
	    foreach node_id $nlist {
		GiD_File fprintf $deminletchannel "  $node_id"
	    }
	    GiD_File fprintf $deminletchannel "%s" "  End SubModelPartNodes"
	    GiD_File fprintf $deminletchannel "%s" "End SubModelPart"
	    GiD_File fprintf $deminletchannel "%s" ""
	    incr dem_ref_to_props_number

	}
	if {$type_of_motion=="FromATable"} {
	    set properties_path "${basexpath}//c.[list ${cgroupid}]//c.MainProperties"
	    set filename [::xmlutils::setXml "${properties_path}//i.VelocitiesFilename" dv]
	    GiD_File fprintf $deminletchannel "Begin Table $TableNumber TIME VELOCITY"
	    set file_open [open [file native [file join [::KUtils::GetPaths "PDir"] $filename]] r]
	    set file_data [read $file_open]
	    close $file_open
	    GiD_File fprintf -nonewline $deminletchannel $file_data
	    GiD_File fprintf $deminletchannel "End Table"
	    GiD_File fprintf $deminletchannel ""
	}
    }
    GiD_File fprintf $deminletchannel "Begin Table 0 TIME VELOCITY"
    GiD_File fprintf $deminletchannel "0.0  0.0"
    GiD_File fprintf $deminletchannel "1.0  0.0"
    GiD_File fprintf $deminletchannel "End Table"
    GiD_File fprintf $deminletchannel ""

    # For debug
    if {!$::wkcf::pflag} {
	set endtime [clock seconds]
	set ttime [expr $endtime-$inittime]
	# WarnWinText "endtime:$endtime ttime:$ttime"
	WarnWinText "Write DEM-Inlet group using the new mesh format: [::KUtils::Duration $ttime]"
    }
}

proc ::wkcf::GetClusterFileNameAndReplaceInletElementType {inlet_element_type} {
    if {$inlet_element_type eq "LineCluster3D"} {
	set inlet_element_type "Cluster3D"
	set cluster_file_name "linecluster3D.clu"
    } elseif {$inlet_element_type eq "RingCluster3D"} {
	set inlet_element_type "Cluster3D"
	set cluster_file_name "ringcluster3D.clu"
    } elseif {$inlet_element_type eq "Wheat5Cluster3D"} {
	set inlet_element_type "Cluster3D"
	set cluster_file_name "wheat5cluster3D.clu"
    } elseif {$inlet_element_type eq "SoyBeanCluster3D"} {
	set inlet_element_type "Cluster3D"
	set cluster_file_name "soybeancluster3D.clu"
    } elseif {$inlet_element_type eq "CornKernel3Cluster3D"} {
	set inlet_element_type "Cluster3D"
	set cluster_file_name "corn3cluster3D.clu"
    } elseif {$inlet_element_type eq "CornKernelCluster3D"} {
	set inlet_element_type "Cluster3D"
	set cluster_file_name "cornkernelcluster3D.clu"
    } elseif {$inlet_element_type eq "Rock1Cluster3D"} {
	set inlet_element_type "Cluster3D"
	set cluster_file_name "rock1cluster3D.clu"
    } elseif {$inlet_element_type eq "Rock2Cluster3D"} {
	set inlet_element_type "Cluster3D"
	set cluster_file_name "rock2cluster3D.clu"
    } elseif {$inlet_element_type eq "Ballast1Cluster3D"} {
	set inlet_element_type "Cluster3D"
	set cluster_file_name "ballast1cluster3D.clu"
    } elseif {$inlet_element_type eq "Ballast1Cluster3Dred"} {
	set inlet_element_type "Cluster3D"
	set cluster_file_name "ballast1cluster3Dred.clu"
    } elseif {$inlet_element_type eq "Ballast2Cluster3D"} {
	set inlet_element_type "Cluster3D"
	set cluster_file_name "ballast2cluster3D.clu"
    } elseif {$inlet_element_type eq "Ballast2Cluster3Dred"} {
	set inlet_element_type "Cluster3D"
	set cluster_file_name "ballast2cluster3Dred.clu"
    } elseif {$inlet_element_type eq "Ballast3Cluster3D"} {
	set inlet_element_type "Cluster3D"
	set cluster_file_name "ballast3cluster3D.clu"
    } elseif {$inlet_element_type eq "Ballast3Cluster3Dred"} {
	set inlet_element_type "Cluster3D"
	set cluster_file_name "ballast3cluster3Dred.clu"
    } elseif {$inlet_element_type eq "Ballast4Cluster3D"} {
	set inlet_element_type "Cluster3D"
	set cluster_file_name "ballast4cluster3D.clu"
    } elseif {$inlet_element_type eq "Ballast4Cluster3Dred"} {
	set inlet_element_type "Cluster3D"
	set cluster_file_name "ballast4cluster3Dred.clu"
    } elseif {$inlet_element_type eq "Ballast5Cluster3D"} {
	set inlet_element_type "Cluster3D"
	set cluster_file_name "ballast5cluster3D.clu"
    } elseif {$inlet_element_type eq "Ballast5Cluster3Dred"} {
	set inlet_element_type "Cluster3D"
	set cluster_file_name "ballast5cluster3Dred.clu"
    } elseif {$inlet_element_type eq "Ballast6Cluster3D"} {
	set inlet_element_type "Cluster3D"
	set cluster_file_name "ballast6cluster3D.clu"
    } elseif {$inlet_element_type eq "Ballast6Cluster3Dred"} {
	set inlet_element_type "Cluster3D"
	set cluster_file_name "ballast6cluster3Dred.clu"
    } elseif {$inlet_element_type eq "SoyBean3Cluster3D"} {
	set inlet_element_type "Cluster3D"
	set cluster_file_name "soybean3cluster3D.clu"
    } elseif {$inlet_element_type eq "CapsuleCluster3D"} {
	set inlet_element_type "Cluster3D"
	set cluster_file_name "capsulecluster3D.clu"
    } elseif {$inlet_element_type eq "SingleSphereCluster3D"} {
	set inlet_element_type "Cluster3D"
	set cluster_file_name "singlespherecluster3D.clu"
    } elseif {$inlet_element_type eq "Rock3RefinedCluster3D"} {
    set inlet_element_type "Cluster3D"
    set cluster_file_name "rock3refinedcluster3D.clu"
    }

    return [list $inlet_element_type $cluster_file_name]
}

proc ::wkcf::GetDEMFEMElementName {firstelementid} {
    variable ndime
    if {$ndime eq "2D"} {set elemname "RigidEdge3D"
    } else {set elemname "RigidFace3D"}
    #set elemname [append elemname $ndime]
    set elemname [append elemname [lindex [GiD_Mesh get element $firstelementid] 2] ]
    set elemname [append elemname "N"]
    return $elemname
}

proc ::wkcf::GetDSOLIDElementName {firstelementid} {
    set elemtype [lindex [GiD_Mesh get element $firstelementid] 1]
    if {$elemtype eq "Tetrahedra"} {set elemname "TotalLagrangianElement3D4N"
    } elseif {$elemtype eq "Hexahedra"} {set elemname "TotalLagrangianElement3D8N"
    } else {WarnWinText "Element type not supported. Some elements will not be considered"}

    return $elemname
}


proc ::wkcf::GetDSOLIDContactElementName {firstelementid} {
    set elemtype [lindex [GiD_Mesh get element $firstelementid] 1]
    if {$elemtype eq "Triangle"} {set elemname "SolidFace3D3N"
    } else {set elemname "SolidFace3D4N"}
    return $elemname
}

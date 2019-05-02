###############################################################################
#
#    NAME: wkcfgroups.tcl
#
#    PURPOSE: Useful procedures to work with groups
#
#    QUANTECH ATZ-DEVELOPMENT DEPARTMENT
#
#    AUTHOR : G. Socorro
#
#    CREATED AT: 02/04/12
#
#    HISTORY:
#
#     1.0- 18/06/13-G. Socorro, delete the proc WritePythonGroupProperties
#     0.9- 17/06/13-G. Socorro, delete wmethod variable and all related procedures (*_m0,*_m1,*_m2) => now we are using only the new GiD groups
#     0.8- 08/05/13-G. Socorro, correct a bug in the proc WriteGroupMeshProperties_m2 => write the element identifier only for the element type assigned to the group
#     0.7- 07/11/12-J. Garate,  Modification and adaptation on functions: WriteGroupMeshProperties, WriteGroupProperties
#                               Creation of functions using GiD_File fprintf $filechannel "%s" format
#     0.6- 24/09/12-G. Socorro, update the proc WriteGroupMeshProperties_m1 to write the Spalart-Allmaras turbulence model
#     0.5- 23/09/12-G. Socorro, update the proc WriteGroupMeshProperties_m1 to write the drag forces
#     0.4- 22/07/12-G. Socorro, write the group properties using the new mesh format
#     0.3- 05/05/12-G. Socorro, write the group properties using the fast method 
#     0.2- 03/04/12-G. Socorro, change Group by GroupNodes
#     0.1- 02/04/12-G. Socorro, create a base source code from wkcf.tcl
#
###############################################################################

proc ::wkcf::WriteGroupMeshProperties {AppId} {
    # ABSTRACT: Write the group properties to the end of the mdpa file
    variable dprops
    variable filechannel
    
    # For debug
    if {!$::wkcf::pflag} {
        set inittime [clock seconds]
    }
    # Write assigned group to the elements
    set meshgroupid 0

    # Init the global mesh identifier list
    set dprops($AppId,AllMeshId) [list]
    
    if {$AppId eq "DEM"} {
	
	# Write inlet group to the kratos mesh format 
	::wkcf::WriteDEMGroupMeshProperties $AppId 

    } else {
    
	# Get all WriteAllGiDGroups from the project settings 
	set WriteAllGiDGroups [::kps::getConfigValue "WriteAllGiDGroups"]
	# wa "WriteAllGiDGroups:$WriteAllGiDGroups"
	
	if {$WriteAllGiDGroups || $AppId ne "Fluid"} {
	# Write all defined groups properties
	
	    # Write all GiD group in the Kratos mesh format
	    ::wkcf::WriteAllGiDGroupMesh $AppId
	} else {
	    
	    # Write only for the Kratos defined elements
	    ::wkcf::WriteAssignedElementGroupMesh $AppId
	}
    }
}

proc ::wkcf::WriteAssignedElementGroupMesh {AppId} {
    # ABSTRACT: Write the group properties to the end of the mdpa file
    variable dprops
    variable filechannel
    
    # For debug
    if {!$::wkcf::pflag} {
        set inittime [clock seconds]
    }
    # Write assigned group to the elements
    set meshgroupid 0

    # Init the global mesh identifier list
    set dprops($AppId,AllMeshId) [list]
    
    # Write only for the Kratos defined elements
    if {([info exists dprops($AppId,AllKElemId)]) && ([llength $dprops($AppId,AllKElemId)])} {
	
	# For all defined kratos elements        
	foreach celemid $dprops($AppId,AllKElemId) {
	    # Check for all defined group identifier for this element
	    if {([info exists dprops($AppId,KElem,$celemid,AllGroupId)]) && ([llength $dprops($AppId,KElem,$celemid,AllGroupId)])} {
		# wa "celemid:$celemid"
		# For all defined group identifier for this element
		foreach cgroupid $dprops($AppId,KElem,$celemid,AllGroupId) {
		    # wa "cgroupid:$cgroupid"
		    # wa "GProps:$dprops($AppId,KElem,$celemid,$cgroupid,GProps)"
		    
		    # Get the GiD entity type, element type and property identifier
		    lassign $dprops($AppId,KElem,$celemid,$cgroupid,GProps) GiDEntity GiDElemType PropertyId KEKWord nDim
		    
		    # Write nodes
		    set wnodes 0
		    if {[GiD_EntitiesGroups get $cgroupid nodes -count]} {
			incr meshgroupid 1
			# Create the meshid-group identifier mapping
			lappend dprops($AppId,AllMeshId) $meshgroupid
			set dprops($AppId,Mesh,$cgroupid,MeshIdGroup) $meshgroupid
			
			# Write mesh properties for this group
			GiD_File fprintf $filechannel "%s" "Begin Mesh $meshgroupid \/\/ GUI group identifier: $cgroupid"
			# Write nodes
			GiD_File fprintf $filechannel "%s" " "
			GiD_File fprintf $filechannel "%s" " Begin MeshNodes"
			foreach node_id [GiD_EntitiesGroups get $cgroupid nodes] {
			    GiD_File fprintf $filechannel "%10i" $node_id
			}
			GiD_File fprintf $filechannel "%s" " End MeshNodes"
			set wnodes 1
		    }
		    
		    # Write elements                
		    set welements 0
		    if {[GiD_EntitiesGroups get $cgroupid elements -count]} {                
			GiD_File fprintf $filechannel "%s" " "
			GiD_File fprintf $filechannel "%s" " Begin MeshElements"
			foreach elem_id [GiD_EntitiesGroups get $cgroupid elements -element_type $GiDElemType] {                           
			    GiD_File fprintf $filechannel %10i $elem_id
			}
			GiD_File fprintf $filechannel "%s" " End MeshElements"
			GiD_File fprintf $filechannel "%s" " "
			set welements 1
		    }
		    
		    # Write end mesh
		    if {($wnodes)&&($welements)} {
			GiD_File fprintf $filechannel "%s" "End Mesh"
			GiD_File fprintf $filechannel "%s" ""
		    }
		}
	    }
	}
	
	# For all defined group in the drag calculation
	# Get the value
	set basexpath "$AppId//c.Results//c.DragOptions"
	set dragproplist [::xmlutils::setXmlContainerIds $basexpath]
	# wa "dragproplist:$dragproplist"
	foreach cgroupid $dragproplist {
	    if {[GiD_EntitiesGroups get $cgroupid nodes -count]} {
		incr meshgroupid 1
		# Create the meshid-group identifier mapping
		lappend dprops($AppId,AllMeshId) $meshgroupid
		set dprops($AppId,Mesh,$cgroupid,MeshIdGroup) $meshgroupid
		
		# Write mesh properties for this group
		GiD_File fprintf $filechannel "%s" "Begin Mesh $meshgroupid \/\/ GUI group identifier: $cgroupid"
		# Write nodes
		GiD_File fprintf $filechannel "%s" " "
		GiD_File fprintf $filechannel "%s" " Begin MeshNodes"
		foreach node_id [GiD_EntitiesGroups get $cgroupid nodes] {
		    GiD_File fprintf $filechannel "%10i" $node_id
		}
		GiD_File fprintf $filechannel "%s" " End MeshNodes"
		GiD_File fprintf $filechannel "%s" " "
		GiD_File fprintf $filechannel "%s" "End Mesh"
		GiD_File fprintf $filechannel "%s" ""
	    }
	}
	
	# Check for used Spalart-Allmaras turbulence model
	set rootid "$AppId"
	set cproperty "dv"
	# Get the turbulence properties
	set cxpath "$rootid//c.AnalysisData//i.TurbulenceModel"
	set TurbulenceModel [::xmlutils::setXml $cxpath $cproperty]
	# WarnWinText "TurbulenceModel:$TurbulenceModel"
	if {$TurbulenceModel eq "Spalart-Allmaras"} {
	    variable ndime
	    # Get the values
	    set basexpath "$rootid//c.AnalysisData//c.Spalart-AllmarasGroupId${ndime}"
	    set gproplist [::xmlutils::setXmlContainerIds $basexpath]
	    # wa "gproplist:$gproplist"
	    foreach cgroupid $gproplist {
		# Get the group properties
		set cxpath "${basexpath}//c.[list ${cgroupid}]//c.MainProperties"
		set allgprop [::xmlutils::setXmlContainerPairs $cxpath "" "dv"]
		# wa "allgprop:$allgprop"
		if {[llength $allgprop]} {
		    set Activate [lindex $allgprop 0 1]
		    # wa "Activate:$Activate"
		    if {$Activate} {
			if {[GiD_EntitiesGroups get $cgroupid nodes -count]} {
			    incr meshgroupid 1
			    # Create the meshid-group identifier mapping
			    lappend dprops($AppId,AllMeshId) $meshgroupid
			    set dprops($AppId,Mesh,$cgroupid,MeshIdGroup) $meshgroupid
			    
			    # Write mesh properties for this group
			    GiD_File fprintf $filechannel "%s" "Begin Mesh $meshgroupid \/\/ GUI group identifier: $cgroupid"
			    # Write nodes
			    GiD_File fprintf $filechannel "%s" " "
			    GiD_File fprintf $filechannel "%s" " Begin MeshNodes"
			    foreach node_id [GiD_EntitiesGroups get $cgroupid nodes] {
				GiD_File fprintf $filechannel "%10i" $node_id
			    }
			    GiD_File fprintf $filechannel "%s" " End MeshNodes"
			    GiD_File fprintf $filechannel "%s" " "
			    GiD_File fprintf $filechannel "%s" "End Mesh"
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
	WarnWinText "Write group using the new mesh format: [::KUtils::Duration $ttime]"
    }

}

proc ::wkcf::WriteAllGiDGroupMesh {AppId} {
    # ABSTRACT: Write all GiD group properties to the end of the mdpa file using the new kratos mesh format
    variable dprops
    variable filechannel
    
    # For debug
    if {!$::wkcf::pflag} {
        set inittime [clock seconds]
    }
    # Write assigned group to the elements
    set meshgroupid 0

    # Init the global mesh identifier list
    set dprops($AppId,AllMeshId) [list]

    foreach cgroupid [GiD_Groups list] {
	# Get the group state (write only normal or disabled groups)
	if {[GiD_Groups get state $cgroupid] ne "hidden"} {
	    # Write nodes
	    set wnodes 0
	    if {[GiD_EntitiesGroups get $cgroupid nodes -count]} {
		incr meshgroupid 1
		# Create the meshid-group identifier mapping
		lappend dprops($AppId,AllMeshId) $meshgroupid
		set dprops($AppId,Mesh,$cgroupid,MeshIdGroup) $meshgroupid
		
		# Write mesh properties for this group
		GiD_File fprintf $filechannel "%s" "Begin Mesh $meshgroupid \/\/ GUI group identifier: $cgroupid"
		# Write nodes
		GiD_File fprintf $filechannel "%s" " "
		GiD_File fprintf $filechannel "%s" " Begin MeshNodes"
		foreach node_id [GiD_EntitiesGroups get $cgroupid nodes] {
		    GiD_File fprintf $filechannel "%10i" $node_id
		}
		GiD_File fprintf $filechannel "%s" " End MeshNodes"
		set wnodes 1
	    }
	    
	    # Write elements                
	    set welements 0
	    if {[GiD_EntitiesGroups get $cgroupid elements -count]} {                
		GiD_File fprintf $filechannel "%s" " "
		GiD_File fprintf $filechannel "%s" " Begin MeshElements"
		foreach elem_id [GiD_EntitiesGroups get $cgroupid elements] {                           
		    GiD_File fprintf $filechannel %10i $elem_id
		}
		GiD_File fprintf $filechannel "%s" " End MeshElements"
		GiD_File fprintf $filechannel "%s" " "
		set welements 1
	    }
	    
	    # Write end mesh
	    if {($wnodes)||($welements)} {
		GiD_File fprintf $filechannel "%s" "End Mesh"
		GiD_File fprintf $filechannel "%s" ""
	    }
	}
    }

    # For debug
    if {!$::wkcf::pflag} {
	set endtime [clock seconds]
	set ttime [expr $endtime-$inittime]
	# WarnWinText "endtime:$endtime ttime:$ttime"
	WarnWinText "Write all GiD group using the new mesh format: [::KUtils::Duration $ttime]"
    }
}
















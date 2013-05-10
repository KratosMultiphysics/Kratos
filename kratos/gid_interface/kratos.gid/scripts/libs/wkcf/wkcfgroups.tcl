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
    # ABSTRACT: Write the group properties to the mdpa file
    variable wmethod
    
    # For debug
    if {!$::wkcf::pflag} {
        set inittime [clock seconds]
    }
    switch -exact -- $wmethod {
        "1" {
            ::wkcf::WriteGroupMeshProperties_m1 $AppId
        }
        "2" {
            ::wkcf::WriteGroupMeshProperties_m2 $AppId
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

proc ::wkcf::WriteGroupMeshProperties_m1 {AppId} {
    # ABSTRACT: Write the group properties to the end of the mdpa file
    variable dprops

    #
    # Write assigned group to the elements
    #
    set meshgroupid 0

    # Init the global mesh identifier list
    set dprops($AppId,AllMeshId) [list]
    
    if {([info exists dprops($AppId,AllKElemId)]) && ([llength $dprops($AppId,AllKElemId)])} {

	# For all defined kratos elements        
	foreach celemid $dprops($AppId,AllKElemId) {
	    # Check for all defined group identifier for this element
	    if {([info exists dprops($AppId,KElem,$celemid,AllGroupId)]) && ([llength $dprops($AppId,KElem,$celemid,AllGroupId)])} {
		# wa "celemid:$celemid"
		# For all defined group identifier for this element
		foreach cgroupid $dprops($AppId,KElem,$celemid,AllGroupId) {
		    # Group properties format 
		    set gprop [dict create]
		    set f "%10i\n"
		    set f [subst $f]
		    dict set gprop $cgroupid "$f"
		    if {[write_calc_data nodes -count $gprop]>0} {
			incr meshgroupid 1
			# Create the meshid-group identifier mapping
			lappend dprops($AppId,AllMeshId) $meshgroupid
			set dprops($AppId,Mesh,$cgroupid,MeshIdGroup) $meshgroupid

			# Write mesh properties for this group
			write_calc_data puts "Begin Mesh $meshgroupid \/\/ GUI group identifier: $cgroupid"
			# Write nodes
			write_calc_data puts " "
			write_calc_data puts " Begin MeshNodes"
   			write_calc_data nodes -sorted $gprop
			write_calc_data puts " End MeshNodes"
			# Write elements
			write_calc_data puts " "
			write_calc_data puts " Begin MeshElements"
   			write_calc_data elements -sorted $gprop
			write_calc_data puts " End MeshElements"
			write_calc_data puts " "
			write_calc_data puts "End Mesh"
		        write_calc_data puts ""
		    }
		    unset gprop
		}
	    }
	}

	# For all defined group in the drag calculation
	# Get the value
	set basexpath "$AppId//c.Results//c.DragOptions"
	set dragproplist [::xmlutils::setXmlContainerIds $basexpath]
	# wa "dragproplist:$dragproplist"
	foreach cgroupid $dragproplist {
	    # Group properties format 
	    set gprop [dict create]
	    set f "%10i\n"
	    set f [subst $f]
	    dict set gprop $cgroupid "$f"
	    if {[write_calc_data nodes -count $gprop]>0} {
		incr meshgroupid 1
		# Create the meshid-group identifier mapping
		lappend dprops($AppId,AllMeshId) $meshgroupid
		set dprops($AppId,Mesh,$cgroupid,MeshIdGroup) $meshgroupid
		
		# Write mesh properties for this group
		write_calc_data puts "Begin Mesh $meshgroupid \/\/ GUI group identifier: $cgroupid"
		# Write nodes
		write_calc_data puts " "
		write_calc_data puts " Begin MeshNodes"
		write_calc_data nodes -sorted $gprop
		write_calc_data puts " End MeshNodes"
		write_calc_data puts " "
		write_calc_data puts "End Mesh"
		write_calc_data puts ""
	    }
	    unset gprop
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
		set cxpath "${basexpath}//c.${cgroupid}//c.MainProperties"
		set allgprop [::xmlutils::setXmlContainerPairs $cxpath "" "dv"]
		# wa "allgprop:$allgprop"
		if {[llength $allgprop]} {
		    set Activate [lindex $allgprop 0 1]
		    # wa "Activate:$Activate"
		    if {$Activate} {
			# Group properties format 
			set gprop [dict create]
			set f "%10i\n"
			set f [subst $f]
			dict set gprop $cgroupid "$f"
			if {[write_calc_data nodes -count $gprop]>0} {
			    incr meshgroupid 1
			    # Create the meshid-group identifier mapping
			    lappend dprops($AppId,AllMeshId) $meshgroupid
			    set dprops($AppId,Mesh,$cgroupid,MeshIdGroup) $meshgroupid
			    
			    # Write mesh properties for this group
			    write_calc_data puts "Begin Mesh $meshgroupid \/\/ GUI group identifier: $cgroupid"
			    # Write nodes
			    write_calc_data puts " "
			    write_calc_data puts " Begin MeshNodes"
			    write_calc_data nodes -sorted $gprop
			    write_calc_data puts " End MeshNodes"
			    write_calc_data puts " "
			    write_calc_data puts "End Mesh"
			    write_calc_data puts ""
			}
			unset gprop
		    }
		}
	    }
	}
    }
 }

proc ::wkcf::WriteGroupMeshProperties_m2 {AppId} {
    # ABSTRACT: Write the group properties to the end of the mdpa file
    variable dprops
    variable filechannel

    #
    # Write assigned group to the elements
    #
    set meshgroupid 0

    # Init the global mesh identifier list
    set dprops($AppId,AllMeshId) [list]
    
    if {([info exists dprops($AppId,AllKElemId)]) && ([llength $dprops($AppId,AllKElemId)])} {

        # For all defined kratos elements        
        foreach celemid $dprops($AppId,AllKElemId) {
            # Check for all defined group identifier for this element
            if {([info exists dprops($AppId,KElem,$celemid,AllGroupId)]) && ([llength $dprops($AppId,KElem,$celemid,AllGroupId)])} {
		# wa "celemid:$celemid"
                # For all defined group identifier for this element
                foreach cgroupid $dprops($AppId,KElem,$celemid,AllGroupId) {
		    # wa "cgroupid:$cgroupid"
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
			    # wa "elem_id:$elem_id"
                            GiD_File fprintf $filechannel "%10i" $elem_id
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
                set cxpath "${basexpath}//c.${cgroupid}//c.MainProperties"
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
 }

proc ::wkcf::WriteGroupProperties {AppId} {
    # ABSTRACT: Write the group properties file
    variable wmethod
    
    # For debug
    if {!$::wkcf::pflag} {
        set inittime [clock seconds]
    }
    switch -exact -- $wmethod {
        "0" {
            ::wkcf::WriteGroupProperties_m0 $AppId
        }
        "1" {
            ::wkcf::WriteGroupProperties_m1 $AppId
        }
        "2" {
            ::wkcf::WriteGroupProperties_m2 $AppId
        }
    }
    # For debug
    if {!$::wkcf::pflag} {
        set endtime [clock seconds]
        set ttime [expr $endtime-$inittime]
        # WarnWinText "endtime:$endtime ttime:$ttime"
        WarnWinText "Write group in nodes: [::KUtils::Duration $ttime]"
    }
}

proc ::wkcf::WriteGroupProperties_m1 {AppId} {
    # ABSTRACT: Write the group properties file
    variable dprops

    set filename "GroupDefinition.txt"
    set PDir [::KUtils::GetPaths "PDir"]

    set fullname [file native [file join $PDir $filename]]
    
    # Use the write_calc_data procedure from the GiD kernel
    # Init
    write_calc_data init $fullname

    #
    # Write assigned group to the elements
    #
    set ggroupid 0
     
    if {([info exists dprops($AppId,AllKElemId)]) && ([llength $dprops($AppId,AllKElemId)])} {
	# For all defined kratos elements        
	foreach celemid $dprops($AppId,AllKElemId) {
	    # Check for all defined group identifier for this element
	    if {([info exists dprops($AppId,KElem,$celemid,AllGroupId)]) && ([llength $dprops($AppId,KElem,$celemid,AllGroupId)])} {
		# wa "celemid:$celemid"
		# For all defined group identifier for this element
		foreach cgroupid $dprops($AppId,KElem,$celemid,AllGroupId) {
		    # Group properties format 
		    set gprop [dict create]
		    set f "%10i\n"
		    set f [subst $f]
		    dict set gprop $cgroupid "$f"
		    if {[write_calc_data nodes -count $gprop]>0} {
			incr ggroupid 1
			write_calc_data puts "Begin GroupNodes $ggroupid $cgroupid  \/\/ GUI group identifier: $cgroupid"
			write_calc_data nodes -sorted $gprop
			write_calc_data puts "End GroupNodes"
		        write_calc_data puts ""
		    }
		    unset gprop
		}
	    }
	}
    }
   
    #
    # Write the group assigned to the boundary condition
    #
  
    # Check for all defined condition type
    if {([info exists dprops($AppId,AllBCTypeId)]) && ([llength $dprops($AppId,AllBCTypeId)])} {
	variable gidentitylist; variable useqelem; variable ndime
	# For all defined condition identifier
	foreach ccondid $dprops($AppId,AllBCTypeId) {
	    # wa "ccondid:$ccondid"
	    # Check for all defined group identifier inside this condition type
	    if {([info exists dprops($AppId,BC,$ccondid,AllGroupId)]) && ([llength $dprops($AppId,BC,$ccondid,AllGroupId)])} {
		# For each group of this BC 	
		foreach cgroupid $dprops($AppId,BC,$ccondid,AllGroupId) {
		    # Group properties format for the boundary condition groups
		    set gprop_bc [dict create]
		    set f "%10i\n"
		    set f [subst $f]
		    dict set gprop_bc $cgroupid "$f"
		    if {[write_calc_data nodes -count $gprop_bc]>0} {
			incr ggroupid 1
			write_calc_data puts "Begin GroupNodes $ggroupid $cgroupid  \/\/ GUI group identifier: $cgroupid"
			write_calc_data nodes -sorted $gprop_bc
			write_calc_data puts "End GroupNodes"
		        write_calc_data puts ""
		    }
		    # Unset the dictionaries
		    unset gprop_bc
		}
	    }
	}
    }
     
    # End
    write_calc_data end
}

proc ::wkcf::WriteGroupProperties_m2 {AppId} {
    # ABSTRACT: Write the group properties file
    variable dprops
    variable filechannel

    set filename "GroupDefinition.txt"
    set PDir [::KUtils::GetPaths "PDir"]

    set fullname [file native [file join $PDir $filename]]
    
    # Use the write_calc_data procedure from the GiD kernel
    # Init
    write_calc_data init $fullname

    #
    # Write assigned group to the elements
    #
    set ggroupid 0
     
    if {([info exists dprops($AppId,AllKElemId)]) && ([llength $dprops($AppId,AllKElemId)])} {
        # For all defined kratos elements        
        foreach celemid $dprops($AppId,AllKElemId) {
            # Check for all defined group identifier for this element
            if {([info exists dprops($AppId,KElem,$celemid,AllGroupId)]) && ([llength $dprops($AppId,KElem,$celemid,AllGroupId)])} {
                # wa "celemid:$celemid"
                # For all defined group identifier for this element
                foreach cgroupid $dprops($AppId,KElem,$celemid,AllGroupId) {
                    # Group properties format 
                    set gprop [dict create]
                    set f "%10i\n"
                    set f [subst $f]
                    dict set gprop $cgroupid "$f"
                    if {[write_calc_data nodes -count $gprop]>0} {
                        incr ggroupid 1
                        write_calc_data puts "Begin GroupNodes $ggroupid $cgroupid  \/\/ GUI group identifier: $cgroupid"
                        write_calc_data nodes -sorted $gprop
                        write_calc_data puts "End GroupNodes"
                        write_calc_data puts ""
                    }
                    unset gprop
                }
            }
        }
    }
   
    #
    # Write the group assigned to the boundary condition
    #
  
    # Check for all defined condition type
    if {([info exists dprops($AppId,AllBCTypeId)]) && ([llength $dprops($AppId,AllBCTypeId)])} {
        variable gidentitylist; variable useqelem; variable ndime
        # For all defined condition identifier
        foreach ccondid $dprops($AppId,AllBCTypeId) {
            # wa "ccondid:$ccondid"
            # Check for all defined group identifier inside this condition type
            if {([info exists dprops($AppId,BC,$ccondid,AllGroupId)]) && ([llength $dprops($AppId,BC,$ccondid,AllGroupId)])} {
                # For each group of this BC 	
                foreach cgroupid $dprops($AppId,BC,$ccondid,AllGroupId) {
                    # Group properties format for the boundary condition groups
                    set gprop_bc [dict create]
                    set f "%10i\n"
                    set f [subst $f]
                    dict set gprop_bc $cgroupid "$f"
                    if {[write_calc_data nodes -count $gprop_bc]>0} {
                        incr ggroupid 1
                        write_calc_data puts "Begin GroupNodes $ggroupid $cgroupid  \/\/ GUI group identifier: $cgroupid"
                        write_calc_data nodes -sorted $gprop_bc
                        write_calc_data puts "End GroupNodes"
                        write_calc_data puts ""
                    }
                    # Unset the dictionaries
                    unset gprop_bc
                }
            }
        }
    }
     
    # End
    write_calc_data end
}

proc ::wkcf::WriteGroupProperties_m0 {AppId} {
    # ABSTRACT: Write the group properties file
    variable dprops

    set filename "GroupDefinition.txt"
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
    puts $fileid ""
    
    #
    # Write assigned group to the elements
    #
    set ggroupid 0
    if {([info exists dprops($AppId,AllKElemId)]) && ([llength $dprops($AppId,AllKElemId)])} {
	# For all defined kratos elements        
	foreach celemid $dprops($AppId,AllKElemId) {
	    # Check for all defined group identifier for this element
	    if {([info exists dprops($AppId,KElem,$celemid,AllGroupId)]) && ([llength $dprops($AppId,KElem,$celemid,AllGroupId)])} {
		# wa "celemid:$celemid"
		# For all defined group identifier for this element
		foreach cgroupid $dprops($AppId,KElem,$celemid,AllGroupId) {
		    # Get the GiD entity type, element type and property identifier
		    lassign $dprops($AppId,KElem,$celemid,$cgroupid,GProps) GiDEntity GiDElemType PropertyId KEKWord nDim
		    # WarnWinText "GiDEntity:$GiDEntity GiDElemType:$GiDElemType PropertyId:$PropertyId KEKWord:$KEKWord nDim:$nDim"
		    # Get all defined entities for this group identifier
		    set allnlist [::KUtils::GetDefinedMeshGiDEntities $cgroupid $GiDEntity Nodes]
		    if {[llength $allnlist]} {
			set pgroup ""
			incr ggroupid 1
			puts $fileid "Begin GroupNodes $ggroupid $cgroupid  \/\/ GUI group identifier: $cgroupid"
			foreach nid [lsort -integer $allnlist] {
			    append pgroup "$nid\n" 
			}
			# Write the group nodes entities list
			puts $fileid "[string trimright $pgroup]"
			puts $fileid "End GroupNodes"
			puts $fileid ""
		    }
		    unset allnlist
		}
	    }
	}
    }
  
    #
    # Write the group assigned to the boundary condition
    #
 
    # Check for all defined condition type
    if {([info exists dprops($AppId,AllBCTypeId)]) && ([llength $dprops($AppId,AllBCTypeId)])} {
	variable gidentitylist; variable useqelem; variable ndime
	# For all defined condition identifier
	foreach ccondid $dprops($AppId,AllBCTypeId) {
	    # wa "ccondid:$ccondid"
	    # Check for all defined group identifier inside this condition type
	    if {([info exists dprops($AppId,BC,$ccondid,AllGroupId)]) && ([llength $dprops($AppId,BC,$ccondid,AllGroupId)])} {
		# For each group of this BC 	
		foreach cgroupid $dprops($AppId,BC,$ccondid,AllGroupId) {
		    set NodeList [list]
		    foreach GiDEntity $gidentitylist {
			# Get all defined entities for this group identifier
			switch $GiDEntity {
			    "point" {
				set callnlist [::KUtils::GetDefinedMeshGiDEntities $cgroupid $GiDEntity "Nodes" $useqelem]
				if {[llength $callnlist]} {
				    lappend NodeList {*}$callnlist
				}
			    }
			    "line" - "surface" {
				set callnlist [::KUtils::GetDefinedMeshGiDEntities $cgroupid $GiDEntity "Nodes" $useqelem]
				if {[llength $callnlist]} {
				    lappend NodeList {*}$callnlist
				}
			    }
			    "volume" {
				if {$ndime =="3D"} {
				    set callnlist [::KUtils::GetDefinedMeshGiDEntities $cgroupid $GiDEntity "Nodes" $useqelem]
				    if {[llength $callnlist]} {
					lappend NodeList {*}$callnlist
				    }
				}
			    }
			} 
		    }
		    # wa "groupid:$cgroupid NodeList:$NodeList"
		    if {[llength $NodeList]} {
		       	set pgroup ""
			incr ggroupid 1
			puts $fileid "Begin GroupNodes $ggroupid $cgroupid  \/\/ GUI group identifier: $cgroupid BC Type: $ccondid"
			foreach nid [lsort -integer $NodeList] {
			    append pgroup "$nid\n" 
			}
			# Write the group nodes entities list
			puts $fileid "[string trimright $pgroup]"
			puts $fileid "End GroupNodes"
			puts $fileid ""
		    }
		    unset NodeList
		}
	    }
	}
    }

    close $fileid
}

proc ::wkcf::WritePythonGroupProperties {} {
    # Write the python group properties file
    variable dprops;  variable ActiveAppList

    set filename "ProjectGroups.py"
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
    puts $fileid "from Kratos import *"
    puts $fileid ""
    puts $fileid ""
    puts $fileid "## This part to be configured for each problem"
    # Group identifier
    #project_groups_names=["support", "Load"]
     puts $fileid ""
    # Group nodes identifier
    #project_groups_nodes_ids={"support" : [1,3,6,11],"Load" : [314]}
     puts $fileid ""
    # Group element identifier
    #project_groups_elements_ids={"support" : [74,511,509,21],"Load":[]}

    # Check for all defined kratos elements
    set pgeids "project_groups_elements_ids=\{"
    # For each active application
    foreach AppId $ActiveAppList {
	if {([info exists dprops($AppId,AllKElemId)]) && ($dprops($AppId,AllKElemId)>0)} {
	    # For all defined kratos elements        
	    foreach celemid $dprops($AppId,AllKElemId) {
		# Check for all defined group identifier for this element
		if {([info exists dprops($AppId,KElem,$celemid,AllGroupId)]) && ($dprops($AppId,KElem,$celemid,AllGroupId)>0)} {
		    # WarnWinText "celemid:$celemid"
		    # For all defined group identifier for this element
		    foreach cgroupid $dprops($AppId,KElem,$celemid,AllGroupId) {
		        # Get the GiD entity type, element type and property identifier
		        lassign $dprops($AppId,KElem,$celemid,$cgroupid,GProps) GiDEntity GiDElemType PropertyId KEKWord nDim
		        # WarnWinText "GiDEntity:$GiDEntity GiDElemType:$GiDElemType PropertyId:$PropertyId KEKWord:$KEKWord nDim:$nDim"
		        
		        # Get all defined entities for this group identifier
		        set allelist [::KUtils::GetDefinedMeshGiDEntities $cgroupid $GiDEntity]
		        # WarnWinText "alllist:$allelist"
		        if {[llength $allelist]>0} {
		            append pgeids "\"$cgroupid\":\[[join $allelist ,]\],"
		        }
		    }
		}
	    }
	}
    }
    set findcomma [string index $pgeids end]
    # WarnWinText "findcomma:$findcomma"
    if {$findcomma ==","} {
	set pgeids [string range $pgeids 0 end-1]
    }
    append pgeids "\}"
    puts $fileid "$pgeids"
    
    puts $fileid ""
    # Group condition identifier
    #project_groups_conditions_ids={"support" : [],"Load":[]}

    # This part is always the same for all problems
    puts $fileid "## This part is always the same for all problems"
    puts $fileid ""
    puts $fileid "model_groups = \{\}"
    puts $fileid ""
    puts $fileid "def InitializeGroups(model_part):"
    puts $fileid "    for group_name in project_groups_names:"
    puts $fileid "        mesh = Mesh()"
    puts $fileid ""
    puts $fileid "        for node_id in project_groups_nodes_ids\[group_name\]:"
    puts $fileid "            mesh.Nodes\[node_id\]=model_part.Nodes\[node_id\]"
    puts $fileid ""
    puts $fileid "        for element_id in project_groups_elements_ids\[group_name\]:"
    puts $fileid "            mesh.Elements\[element_id\]=model_part.Elements\[element_id\]"
    puts $fileid ""
    puts $fileid "        for condition_id in project_groups_conditions_ids\[group_name\]:"
    puts $fileid "            mesh.Conditions\[condition_id\]=model_part.Conditions\[condition_id\]"
    puts $fileid ""
    puts $fileid ""
    puts $fileid "        print group_name, \"mesh:\", mesh"
    puts $fileid ""       
    puts $fileid "        model_groups\[group_name\]=mesh"
    puts $fileid ""

    close $fileid
}


















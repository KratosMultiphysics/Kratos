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
#     0.1- 12/09/13-G. Socorro, create a base source code from wkcf.tcl
#
###############################################################################

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
	
	# Write mesh velocity and velocity for each node identifier
	set MeshVelocity 0.0; set Velocity 0.0

	set kxpath "Applications/$AppId"
	set cpropid "0"        

	# Write the group nodal properties
	foreach celemid $dprops($AppId,AllKElemId) {
	    # Check for all defined group identifier for this element
	    if {([info exists dprops($AppId,KElem,$celemid,AllGroupId)]) && ([llength $dprops($AppId,KElem,$celemid,AllGroupId)])} {
		# For all defined group identifier for this element
		foreach cgroupid $dprops($AppId,KElem,$celemid,AllGroupId) {                    
		    # Write mesh velocity value for this group
		    if {[GiD_EntitiesGroups get $cgroupid nodes -count]} {
		        set vkword [::xmlutils::getKKWord $kxpath "MeshVelocity" "kkword"]
		        GiD_File fprintf $filechannel "%s" "Begin NodalData $vkword \/\/ GUI group identifier: $cgroupid"
		        foreach node_id [GiD_EntitiesGroups get $cgroupid nodes] {
		            GiD_File fprintf $filechannel "%10i %4i %4s" $node_id $cpropid $MeshVelocity
		        }
		        GiD_File fprintf $filechannel "%s " "End NodalData"
		        GiD_File fprintf $filechannel ""
		    }

		    # Write velocity value for this group
		    if {[GiD_EntitiesGroups get $cgroupid nodes -count]} {
		        set vkword [::xmlutils::getKKWord $kxpath "Velocity" "kkword"]
		        GiD_File fprintf $filechannel "%s" "Begin NodalData $vkword \/\/ GUI group identifier: $cgroupid"
		        foreach node_id [GiD_EntitiesGroups get $cgroupid nodes] {
		            GiD_File fprintf $filechannel "%10i %4i %4s" $node_id $cpropid $Velocity
		        }
		        GiD_File fprintf $filechannel "%s" "End NodalData"
		        GiD_File fprintf $filechannel "%s" ""
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


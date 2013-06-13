###############################################################################
#
#    NAME: wkcfstructuralanalysis.tcl
#
#    PURPOSE: Useful procedures to work with the fluid application
#
#    QUANTECH ATZ-DEVELOPMENT DEPARTMENT
#
#    AUTHOR : G. Socorro
#
#    CREATED AT: 02/04/12
#
#    HISTORY:
#   
#     1.5- 22/05/13-G. Socorro, correct a bug in the procs WritePuntualLoads and WriteBeamUniformlyDistributedLoads when write condition 
#                               identifier in the case of more than one load group
#     1.4- 17/05/13-G. Socorro, add the proc WriteBeamUniformlyDistributedLoads, update the proc WriteConstitutiveLawsProperties
#                               to write the cross section properties
#     1.3- 24/04/13-G. Socorro, write Rotational_Dofs = True for beam element type
#     1.2- 18/03/13-G. Socorro, correct a bug in the proc WritePressureLoads_m2 change $group by $cgroupid
#     1.1- 17/12/12-J. Garate,  Kratos Path is no longer written at ProjectParameters.py
#     1.0- 12/11/12-J. Garate,  Fixed some errors
#     0.9- 07/11/12-J. Garate,  Modification and adaptation on functions: WriteDispRotBC, WritePressureLoads, WritePuntualLoads, WriteBodyForceValues
#                               Creation of functions using GiD_File fprintf $filechannel "%s" format
#     0.8- 09/10/12-G. Socorro, correct a bug with the output format in the proc WritePressureLoads_m1 
#     0.7- 14/05/12-G. Socorro, modify the proc WriteDispRotBC_m1 to write many groups with displacement/rotation bc
#     0.6- 06/05/12-G. Socorro, update the proc WritePressureLoads to write using the fast method (write_calc_data)   
#     0.5- 05/05/12-G. Socorro, improve the proc WriteDispRotBC to write using the fast method (write_calc_data)
#     0.4- 20/04/12-J. Garate, Acabada ::wkcf::WritePuntualLoads_m1, empezada Pressure
#     0.3- 19/04/12-J. Garate, ::wkcf::WritePuntualLoads_m1
#     0.2- 16/04/12-G. Socorro, ::wkcf::WriteDispRotBC_m1
#     0.1- 02/04/12-G. Socorro, create a base source code from wkcf.tcl
#
###############################################################################

proc ::wkcf::WriteDispRotBC {AppId ccondid kwordlist} {
    # ABSTRACT: Write displacements or rotation boundary condition
    variable wmethod
       
    # For debug
    if {!$::wkcf::pflag} {
        set inittime [clock seconds]
    }
    if {$wmethod eq 0} {
        ::wkcf::WriteDispRotBC_m0 $AppId $ccondid $kwordlist
    } elseif {$wmethod eq 1} {
        ::wkcf::WriteDispRotBC_m1 $AppId $ccondid $kwordlist
    } elseif {$wmethod eq 2} {
        ::wkcf::WriteDispRotBC_m2 $AppId $ccondid $kwordlist
    }
    # For debug
    if {!$::wkcf::pflag} {
        set endtime [clock seconds]
        set ttime [expr $endtime-$inittime]
        # WarnWinText "endtime:$endtime ttime:$ttime"
        WarnWinText "Write structural analysis displacements or rotation boundary condition: [::KUtils::Duration $ttime]"
    }
}

proc ::wkcf::WriteDispRotBC_m0 {AppId ccondid kwordlist} {
    # Write displacements or rotation boundary condition
    variable ndime;    variable gidentitylist
    variable useqelem; variable dprops

    set dxnodeplist [list]
    set dynodeplist [list]
    set dznodeplist [list]
    
    # For all defined group identifier inside this condition type
    foreach cgroupid $dprops($AppId,BC,$ccondid,AllGroupId) {
        # Get the condition properties
        set GProps $dprops($AppId,BC,$ccondid,$cgroupid,GProps)
        # WarnWinText "GProps:$GProps"
        # Assign values
        lassign $GProps fix_x xval fix_y yval fix_z zval
        # WarnWinText "fix_x:$fix_x xval:$xval fix_y:$fix_y yval:$yval fix_z:$fix_z zval:$zval"
        set allnlist [list]
        foreach GiDEntity $gidentitylist {
            # WarnWin "GiDEntity:$GiDEntity cgroupid:$cgroupid"
            # Get all defined entities for this group identifier
            switch $GiDEntity {
                "point" {
                    set callnlist [::KUtils::GetDefinedMeshGiDEntities $cgroupid $GiDEntity "Nodes" $useqelem]
                    if {[llength $callnlist]} {
                        lappend allnlist $callnlist
                    }
                    # WarnWinText "$GiDEntity alllist:$allnlist"
                }
                "line" - "surface" {
                    set callnlist [::KUtils::GetDefinedMeshGiDEntities $cgroupid $GiDEntity "Nodes" $useqelem]
                    if {[llength $callnlist]} {
                        lappend allnlist $callnlist
                    }
                    # WarnWinText "$GiDEntity alllist:$allnlist"
                } 
                "volume" {
                    # Only for displacement BC
                    if {$ccondid =="Displacements"} {
                        set callnlist [::KUtils::GetDefinedMeshGiDEntities $cgroupid $GiDEntity "Nodes" $useqelem]
                        if {[llength $callnlist]} {
                            lappend allnlist $callnlist
                        }
                        # WarnWinText "$GiDEntity alllist:$allnlist"
                    }
                } 
            }
        }
        # WarnWinText "$GiDEntity alllist:$allnlist"
        foreach cprop $allnlist {
            set cprop [lsort -integer -unique $cprop]
            foreach nodeid $cprop {
                # Fix x
                if {$fix_x =="1"} {
                    lappend dxnodeplist "$nodeid 1 $xval"
                }
                # Fix y
                if {$fix_y =="1"} {
                    lappend dynodeplist "$nodeid 1 $yval"
                }
                # Check for 3D problem
                if {$ndime =="3D"} {
                    # Fix z
                    if {$fix_z =="1"} {
                        lappend dznodeplist "$nodeid 1 $zval"
                    }
                }
            }
        }
    }
    unset allnlist
    # WarnWinText "dxnodeplist:$dxnodeplist"
    # WarnWinText "dynodeplist:$dynodeplist"
    # WarnWinText "dznodeplist:$dznodeplist"
    
    # DISPLACEMENT_X or ROTATION_X
    if {[llength $dxnodeplist]>0} {
	set xitem [lindex $kwordlist 0]
	write_calc_data puts "Begin NodalData $xitem"
	foreach citem $dxnodeplist {
	    lassign $citem nodeid fix_x xval
	    set cf "[format "%4i%4i%10.5f" $nodeid $fix_x $xval]"
	    write_calc_data puts "$cf"
	}
	write_calc_data puts "End NodalData"
	write_calc_data puts ""
    }
    
    # DISPLACEMENT_Y or ROTATION_Y
    if {[llength $dynodeplist]>0} {
	set yitem [lindex $kwordlist 1]
	write_calc_data puts "Begin NodalData $yitem"
	foreach citem $dynodeplist {
	    lassign $citem nodeid fix_y yval
	    set cf "[format "%4i%4i%10.5f" $nodeid $fix_y $yval]"
	    write_calc_data puts "$cf"
	}
	write_calc_data puts "End NodalData"
	write_calc_data puts ""
    }
    
    # Check for 3D problem
    if {$ndime =="3D"} {
	# DISPLACEMENT_Z or ROTATION_Z
	if {[llength $dznodeplist]>0} {
	    set zitem [lindex $kwordlist 2]
	    write_calc_data puts "Begin NodalData $zitem"
	    foreach citem $dznodeplist {
		lassign $citem nodeid fix_z zval
		set cf "[format "%4i%4i%10.5f" $nodeid 1 $zval]"
		write_calc_data puts "$cf"
	    }
	    write_calc_data puts "End NodalData"
	    write_calc_data puts ""
	}
    }
    unset dxnodeplist
    unset dynodeplist
    unset dznodeplist
}

proc ::wkcf::WriteDispRotBC_m1 {AppId ccondid kwordlist} {
    # Write displacements or rotation boundary condition
    variable ndime; variable dprops

    # For all defined group identifier inside this condition type
    foreach cgroupid $dprops($AppId,BC,$ccondid,AllGroupId) {
	    
	# Get the condition properties
	lassign $dprops($AppId,BC,$ccondid,$cgroupid,GProps) fix_x xval fix_y yval fix_z zval
	# WarnWinText "fix_x:$fix_x xval:$xval fix_y:$fix_y yval:$yval fix_z:$fix_z zval:$zval

	# Create the dictionary formats
	set gpropdx [dict create]
	set gpropdy [dict create]
	if {$ndime eq "3D"} {
	    set gpropdz [dict create]
	}
		
	# For X
	set f "%10i [format "%4i" $fix_x] [format "%10.5f" $xval]\n"
	set f [subst $f]
	dict set gpropdx $cgroupid $f
	# For Y
	set f "%10i [format "%4i" $fix_y] [format "%10.5f" $yval]\n"
	set f [subst $f]
	dict set gpropdy $cgroupid $f
	if {$ndime eq "3D"} {
	    # For Z
	    set f "%10i [format "%4i" $fix_z] [format "%10.5f" $zval]\n"
	    set f [subst $f]
	    dict set gpropdz $cgroupid $f
	}
	
	# DISPLACEMENT_X or ROTATION_X
	if {[write_calc_data nodes -count $gpropdx]>0} {
	    set xitem [lindex $kwordlist 0]
	    write_calc_data puts "Begin NodalData $xitem"
	    write_calc_data nodes -sorted $gpropdx
	    write_calc_data puts "End NodalData"
	    write_calc_data puts ""        
	}
	
	# DISPLACEMENT_Y or ROTATION_Y
	if {[write_calc_data nodes -count $gpropdy]>0} {
	    set yitem [lindex $kwordlist 1]
	    write_calc_data puts "Begin NodalData $yitem"
	    write_calc_data nodes -sorted $gpropdy
	    write_calc_data puts "End NodalData"
	    write_calc_data puts ""        
	}

	# DISPLACEMENT_Z or ROTATION_Z
	if {($ndime eq "3D") && ([write_calc_data nodes -count $gpropdz]>0)} {
	    set zitem [lindex $kwordlist 2]
	    write_calc_data puts "Begin NodalData $zitem"
	    write_calc_data nodes -sorted $gpropdz
	    write_calc_data puts "End NodalData"
	    write_calc_data puts ""        
	}

	# Unset Local Variables
	unset gpropdx
	unset gpropdy
	if {$ndime ne "2D"} {
	    unset gpropdz
	}
    }  
   
}

proc ::wkcf::WriteDispRotBC_m2 {AppId ccondid kwordlist} {
    # Write displacements or rotation boundary condition
    variable ndime; variable dprops
    variable filechannel
    
    # For all defined group identifier inside this condition type
    foreach cgroupid $dprops($AppId,BC,$ccondid,AllGroupId) {
            
        # Get the condition properties
        lassign $dprops($AppId,BC,$ccondid,$cgroupid,GProps) fix_x xval fix_y yval fix_z zval
        # WarnWinText "fix_x:$fix_x xval:$xval fix_y:$fix_y yval:$yval fix_z:$fix_z zval:$zval"
       
        # DISPLACEMENT_X or ROTATION_X
        set nodes [GiD_EntitiesGroups get $cgroupid nodes]
        
        if {[llength $nodes]} {
            set xitem [lindex $kwordlist 0]
            GiD_File fprintf $filechannel "%s" "Begin NodalData $xitem"
            foreach nodeid $nodes {
                GiD_File fprintf $filechannel "%10i %4i %10f" $nodeid $fix_x $xval
            }
            GiD_File fprintf $filechannel "End NodalData"
            GiD_File fprintf $filechannel "" 
        }
        
        # DISPLACEMENT_Y or ROTATION_Y
        if {[llength $nodes]} {
            set yitem [lindex $kwordlist 1]
            GiD_File fprintf $filechannel "%s" "Begin NodalData $yitem"
            foreach nodeid $nodes {
                GiD_File fprintf $filechannel "%10i %4i %10f" $nodeid $fix_y $yval
            }
            GiD_File fprintf $filechannel "End NodalData"
            GiD_File fprintf $filechannel "" 
        }

        # DISPLACEMENT_Z or ROTATION_Z
        if {($ndime eq "3D") && ([llength $nodes])} {
            set zitem [lindex $kwordlist 2]
            GiD_File fprintf $filechannel "%s" "Begin NodalData $zitem"
            foreach nodeid $nodes {
                GiD_File fprintf $filechannel "%10i %4i %10f" $nodeid $fix_z $zval
            }
            GiD_File fprintf $filechannel "End NodalData"
            GiD_File fprintf $filechannel "" 
        }
    }  
   
}

proc ::wkcf::WriteLoads {AppId} {
    # Write the load properties block
    variable dprops     

    # Check for all defined load type
    if {([info exists dprops($AppId,AllLoadTypeId)]) && ([llength $dprops($AppId,AllLoadTypeId)]>0)} {
        set kwxpath "Applications/$AppId"
        # For all defined load identifier
        foreach cloadtid $dprops($AppId,AllLoadTypeId) {
	    # WarnWinText "cloadtid:$cloadtid"
            # Check for all defined group identifier inside this load type
            if {([info exists dprops($AppId,Loads,$cloadtid,AllGroupId)]) && ([llength $dprops($AppId,Loads,$cloadtid,AllGroupId)])} {
            # Select the load type
            switch -exact -- $cloadtid {
                "Puntual"
                {
                    # Concentrate or puntual loads (forces or moments)
                    set kwordlist [list [::xmlutils::getKKWord $kwxpath "Fx"] [::xmlutils::getKKWord $kwxpath "Fy"] [::xmlutils::getKKWord $kwxpath "Fz"]]
                    # Write puntual loads
                    ::wkcf::WritePuntualLoads $AppId $cloadtid $kwordlist
                }
		"BeamUniformlyDistributed"
                {
		    # Write beam uniformly distributed loads
                    ::wkcf::WriteBeamUniformlyDistributedLoads $AppId $cloadtid
                }
                "Pressure" {
                    # Pressure loads
                    ::wkcf::WritePressureLoads $AppId $cloadtid
                }
            }
            }
        }
    }
}

proc ::wkcf::WriteBeamUniformlyDistributedLoads {AppId cloadtid} {
    # ASBTRACT: Write beam uniformly distributed loads (forces)
    variable wmethod
    
    # For debug
    if {!$::wkcf::pflag} {
        set inittime [clock seconds]
    }
    if {$wmethod eq 2} {
	::wkcf::WriteBeamUniformlyDistributedLoads_m2 $AppId $cloadtid
    } 
    # For debug
    if {!$::wkcf::pflag} {
        set endtime [clock seconds]
        set ttime [expr $endtime-$inittime]
        # WarnWinText "endtime:$endtime ttime:$ttime"
        WarnWinText "Write structural analysis beam uniformly distributed loads (forces): [::KUtils::Duration $ttime]"
    }
}

proc ::wkcf::WriteBeamUniformlyDistributedLoads_m2 {AppId cloadtid} {
    # ABSTRACT: Write beam uniformly distributed loads (forces)
    variable ndime;  variable dprops
    variable filechannel;  variable sa_icondid 
    
    # Kratos key word xpath
    set kwxpath "Applications/$AppId"
    # Set the GiD element type
    set GiDElemType "Line"
    # Set the line force 3d-3n condition keyword
    set LineForce3D2N "LineForce3D2N"
    # Set the reference property id
    set RefPropId "1"
      
    # For all defined group identifier inside this load type
    foreach cgroupid $dprops($AppId,Loads,$cloadtid,AllGroupId) {
	# wa "cgroupid:$cgroupid"
        if {[GiD_EntitiesGroups get $cgroupid elements -count -element_type $GiDElemType]} {
            # Write LineForce3D2N condition
            GiD_File fprintf $filechannel "%s" "Begin Conditions $LineForce3D2N // GUI beam uniformly distributed load group identifier: $cgroupid"
	    foreach elem_id [GiD_EntitiesGroups get $cgroupid elements -element_type $GiDElemType] {
                set nodes [lrange [GiD_Mesh get element $elem_id] 3 end]
                set N1 [lindex $nodes 0]
                set N2 [lindex $nodes 1]
		incr sa_icondid
                set cf "[format "%4i%4i%8i%8i" $sa_icondid $RefPropId $N1 $N2]"
                GiD_File fprintf $filechannel "%s" "$cf"
            }
            GiD_File fprintf $filechannel "%s" "End Conditions"
            GiD_File fprintf $filechannel ""
        }
    }


    # Write uniformly distributed load values for all nodes inside a group identifier 
    foreach cgroupid $dprops($AppId,Loads,$cloadtid,AllGroupId) {
	
        # Get the load properties
        set GProps $dprops($AppId,Loads,$cloadtid,$cgroupid,GProps)
	# WarnWinText "cgroupid:$cgroupid GProps:$GProps"
        # Assign values
        foreach item $GProps {
            foreach {cvar cval} $item {
                set $cvar $cval
		# wa "cvar:$cvar cval:$cval"
		
		# Get the current keyword
		set ckword [::xmlutils::getKKWord $kwxpath "$cvar"]
		# WarnWinText "ckword:$ckword"
		if {([GiD_EntitiesGroups get $cgroupid nodes -count]>0)&&($ndime=="3D")} {
		    
		    if {$cval !="0.0"} {
			# Write beam uniformly properties values
			GiD_File fprintf $filechannel "%s" "Begin NodalData $ckword // GUI beam uniformly distributed load group identifier: $cgroupid"
		        foreach node_id [GiD_EntitiesGroups get $cgroupid nodes] {
		            GiD_File fprintf $filechannel "%10i %6i %20.10f" $node_id $RefPropId $cval
		        }
		        GiD_File fprintf $filechannel "%s" "End NodalData"
		        GiD_File fprintf $filechannel ""
		    }
		}
            }
        } 
    }
}

proc ::wkcf::WritePressureLoads {AppId cloadtid} {
    # ABSTRACT: Write pressure loads (positive or negative shell pressure)
    variable wmethod
       
    # For debug
    if {!$::wkcf::pflag} {
        set inittime [clock seconds]
    }
    if {$wmethod eq 1} {
        ::wkcf::WritePressureLoads_m1 $AppId $cloadtid
    } elseif {$wmethod eq 0} {
        ::wkcf::WritePressureLoads_m0 $AppId $cloadtid
    } elseif {$wmethod eq 2} {
        ::wkcf::WritePressureLoads_m2 $AppId $cloadtid
    }
    # For debug
    if {!$::wkcf::pflag} {
        set endtime [clock seconds]
        set ttime [expr $endtime-$inittime]
        # WarnWinText "endtime:$endtime ttime:$ttime"
        WarnWinText "Write structural analysis write pressure loads (positive or negative shell pressure): [::KUtils::Duration $ttime]"
    }
}

proc ::wkcf::WritePressureLoads_m1 {AppId cloadtid} {
    # ABSTRACT: Write pressure loads (positive or negative shell pressure)
    variable ndime;  variable dprops
    variable sa_icondid

    # Kratos key word xpath
    set kxpath "Applications/$AppId"
    # Set the GiD element type
    set GiDElemType "Triangle"
    # Set the face 3d-3n condition keyword
    set face3d3nkword "Face3D3N"
    # Set the reference property id
    set RefPropId "1"
       
    # For all defined group identifier inside this load type
    foreach cgroupid $dprops($AppId,Loads,$cloadtid,AllGroupId) {
    	set gpressface [dict create]
	# Diccionario gpressface
	set f "%8i%8i%8i%8i\n"
	set f [subst $f]
	dict set gpressface $cgroupid $f
	if {[write_calc_data has_elements -elemtype $GiDElemType $gpressface]} {
	    # Write Face3D3N condition
	    write_calc_data puts "Begin Conditions $face3d3nkword // GUI pressure load group identifier: $cgroupid"
	    foreach {elemid N1 N2 N3} [write_calc_data connectivities -return -elements_faces elements $gpressface] {
		incr sa_icondid
		set cf "[format "%4i%4i%8i%8i%8i" $sa_icondid $RefPropId $N1 $N2 $N3]"
		write_calc_data puts "$cf"
	    }
	    write_calc_data puts "End Conditions"
	    write_calc_data puts ""
	}
	# Unset the local dictionary
	unset gpressface
    }


    # Write pressure values for all nodes inside a group identifier 
    foreach cgroupid $dprops($AppId,Loads,$cloadtid,AllGroupId) {

	# Get the load properties
	set GProps $dprops($AppId,Loads,$cloadtid,$cgroupid,GProps)
	# WarnWinText "cgroupid:$cgroupid GProps:$GProps"
	# Assign values
	foreach item $GProps {
	    foreach {cvar cval} $item {
		set $cvar $cval
	    }
	} 

	# WarnWinText "FixPressure:$FixPressure PressureType:$PressureType PressureValue:$PressureValue"

	# Create the dictionary
	set gprop [dict create]

	# Node format
	set f "%10i"
	set f [subst $f]
	dict set gprop $cgroupid $f
	if {([write_calc_data nodes -count $gprop]>0)&&($ndime=="3D")} {
	    # Get the current pressure keyword 
	    set kwid "${PressureType}Pressure"
	    set ckword [::xmlutils::getKKWord $kxpath $kwid]
	    # WarnWinText "ckword:$ckword"
		
	    if {$PressureValue !="0.0"} {
		set f "%10i [format "%6i%20.10f" $FixPressure $PressureValue]\n"
		set f [subst $f]
		dict set gprop $cgroupid $f
		# Set the real pressure value
		# set PressureValue [expr double($PressureValue)/$ngroupnodes]
		# Write the pressure values
		# Pressure properties
		write_calc_data puts "Begin NodalData $ckword // GUI pressure load group identifier: $cgroupid"
		write_calc_data nodes -sorted $gprop
		write_calc_data puts "End NodalData"
		write_calc_data puts ""
	    }
	}
	# Unset the local dictionary
	unset gprop
    }
}

proc ::wkcf::WritePressureLoads_m2 {AppId cloadtid} {
    # ABSTRACT: Write pressure loads (positive or negative shell pressure)
    variable ndime;  variable dprops
    variable filechannel; variable sa_icondid

    # Kratos key word xpath
    set kxpath "Applications/$AppId"
    # Set the GiD element type
    set GiDElemType "Triangle"
    # Set the face 3d-3n condition keyword
    set face3d3nkword "Face3D3N"
    # Set the reference property id
    set RefPropId "1"
       
    # For all defined group identifier inside this load type
    foreach cgroupid $dprops($AppId,Loads,$cloadtid,AllGroupId) {
        if {[GiD_EntitiesGroups get $cgroupid elements -count -element_type $GiDElemType]} {
            # Write Face3D3N condition
            GiD_File fprintf $filechannel "%s" "Begin Conditions $face3d3nkword // GUI pressure load group identifier: $cgroupid"
            foreach elem_id [GiD_EntitiesGroups get $cgroupid elements -element_type $GiDElemType] {
                set nodes [lrange [GiD_Mesh get element $elem_id] 3 end]
                set N1 [lindex $nodes 0]
                set N2 [lindex $nodes 1]
                set N3 [lindex $nodes 2]
                incr sa_icondid
                set cf "[format "%4i%4i%8i%8i%8i" $sa_icondid $RefPropId $N1 $N2 $N3]"
                GiD_File fprintf $filechannel "%s" "$cf"
            }
            GiD_File fprintf $filechannel "%s" "End Conditions"
            GiD_File fprintf $filechannel ""
        }
    }


    # Write pressure values for all nodes inside a group identifier 
    foreach cgroupid $dprops($AppId,Loads,$cloadtid,AllGroupId) {

        # Get the load properties
        set GProps $dprops($AppId,Loads,$cloadtid,$cgroupid,GProps)
        # WarnWinText "cgroupid:$cgroupid GProps:$GProps"
        # Assign values
        foreach item $GProps {
            foreach {cvar cval} $item {
                set $cvar $cval
            }
        } 

        # WarnWinText "FixPressure:$FixPressure PressureType:$PressureType PressureValue:$PressureValue"
        if {([GiD_EntitiesGroups get $cgroupid nodes -count]>0)&&($ndime=="3D")} {
            # Get the current pressure keyword 
            set kwid "${PressureType}Pressure"
            set ckword [::xmlutils::getKKWord $kxpath $kwid]
            # WarnWinText "ckword:$ckword"
            
            if {$PressureValue !="0.0"} {
                # Set the real pressure value
                # set PressureValue [expr double($PressureValue)/$ngroupnodes]
                # Write the pressure values
                # Pressure properties
                GiD_File fprintf $filechannel "%s" "Begin NodalData $ckword // GUI pressure load group identifier: $cgroupid"
                foreach node_id [GiD_EntitiesGroups get $cgroupid nodes] {
                    GiD_File fprintf $filechannel "%10i %6i %20.10f" $node_id $FixPressure $PressureValue
                }
                GiD_File fprintf $filechannel "%s" "End NodalData"
                GiD_File fprintf $filechannel ""
            }
        }
    }
}

proc ::wkcf::WritePressureLoads_m0 {AppId cloadtid} {
    # Write pressure loads (positive or negative shell pressure)
    variable ndime;    variable gidentitylist
    variable useqelem; variable dprops
    variable sa_icondid

    # Kratos key word xpath
    set kxpath "Applications/$AppId"
    # Set GiD entity
    set GiDEntity "surface"
    # Set the GiD element type
    set GiDElemType "Triangle"
    # Set the face 3d-3n condition keyword
    set face3d3nkword "Face3D3N"
    # Set the reference property id
    set RefPropId "1"
    set cnodeglist [list]
    # Node dictionary
    set ndict [dict create]
    # Create an element dictionary
    set edict [dict create]

    # For all defined group identifier inside this load type
    foreach cgroupid $dprops($AppId,Loads,$cloadtid,AllGroupId) {
	# Get all defined entities for this group identifier
	set allelist [::KUtils::GetDefinedMeshGiDEntities $cgroupid $GiDEntity]
	# WarnWinText "cgroupid:$cgroupid alllist:$allelist"
	if {([llength $allelist])&&($ndime=="3D")} {
	    # Create a node group dictionary
	    set nc [dict create]
	    # Non repeated nodes identifier
	    # Get all defined nodes
	    foreach elemid $allelist {
		# Update the element group (delete repeated identifier)
		dict set edict $elemid $cgroupid
		# Get the element properties
		lassign [lrange [GiD_Info Mesh Elements $GiDElemType $elemid] 1 end-1] N1 N2 N3
		# Update node group dictionary variable (delete repeated nodes)
		dict set nc $N1 ""
		dict set nc $N2 ""
		dict set nc $N3 ""
	    }
	    # Update nodal group dictionary
	    dict set ndict $cgroupid $nc
	    unset nc
	}
    }


    # Write Face3D3N condition
    foreach cgroupid $dprops($AppId,Loads,$cloadtid,AllGroupId) {
	
	# Init the element group list
	set celemglist [list]
	dict for {elemid dgroupid} $edict {
	    if {$cgroupid ==$dgroupid} {
		lappend celemglist $elemid
	    }
	}
	set ngroupelems [llength $celemglist]
	# WarnWinText "ngroupelems:$ngroupelems"
	if {($ngroupelems) && ($ndime=="3D")} {
	    write_calc_data puts "Begin Conditions $face3d3nkword // GUI pressure load group identifier: $cgroupid"
	    # Get all defined nodes
	    foreach elemid $celemglist {
		incr sa_icondid 1
		# Get the element properties
		lassign [lrange [GiD_Info Mesh Elements $GiDElemType $elemid] 1 end-1] N1 N2 N3
		set cf "[format "%4i%4i%8i%8i%8i" $sa_icondid $RefPropId $N1 $N2 $N3]"
		write_calc_data puts "$cf"
	    }
	    write_calc_data puts "End Conditions"
	    write_calc_data puts ""
	}
    }
    unset edict

    # Write pressure values for all nodes inside a group identifier 
    foreach cgroupid $dprops($AppId,Loads,$cloadtid,AllGroupId) {

	# Get the load properties
	set GProps $dprops($AppId,Loads,$cloadtid,$cgroupid,GProps)
	# WarnWinText "cgroupid:$cgroupid GProps:$GProps"
	# Assign values
	foreach item $GProps {
	    foreach {cvar cval} $item {
		set $cvar $cval
	    }
	} 

	# WarnWinText "FixPressure:$FixPressure PressureType:$PressureType PressureValue:$PressureValue"

	# Init the node group list
	if {[llength [dict get $ndict $cgroupid]]} {
	    set cnodeglist [list]
	    foreach {nodeid empty} [dict get $ndict $cgroupid] {
		lappend cnodeglist $nodeid
	    }
	    
	    set ngroupnodes [llength $cnodeglist]
	    # WarnWinText "ngroupnodes:$ngroupnodes"
	    if {($ngroupnodes) && ($ndime=="3D")} {
		# Get the current pressure keyword 
		set kwid "${PressureType}Pressure"
		set ckword [::xmlutils::getKKWord $kxpath $kwid]
		# WarnWinText "ckword:$ckword"
		
		if {$PressureValue !="0.0"} {
		    # Set the real pressure value
		    # set PressureValue [expr double($PressureValue)/$ngroupnodes]
		    # Write the pressure values
		    # Pressure properties
		    write_calc_data puts "Begin NodalData $ckword // GUI pressure load group identifier: $cgroupid"
		    foreach nodeid [lsort -integer $cnodeglist] {
		        set cf "[format "%6i%4i%20.10f" $nodeid $FixPressure $PressureValue]"
		        write_calc_data puts "$cf"
		    }
		    write_calc_data puts "End NodalData"
		    write_calc_data puts ""
		}
		unset cnodeglist
	    }
	}
    }
    unset ndict
}

proc ::wkcf::WritePuntualLoads {AppId cloadtid kwordlist} {
    # ASBTRACT: Write concentrated loads (puntual force or moment)
    variable wmethod
    
    # For debug
    if {!$::wkcf::pflag} {
        set inittime [clock seconds]
    }
    if {$wmethod eq 1} {
        ::wkcf::WritePuntualLoads_m1 $AppId $cloadtid $kwordlist
    } elseif {$wmethod eq 2} {
        ::wkcf::WritePuntualLoads_m2 $AppId $cloadtid $kwordlist
    } elseif {$wmethod eq 0} {
        ::wkcf::WritePuntualLoads_m0 $AppId $cloadtid $kwordlist
    }
    # For debug
    if {!$::wkcf::pflag} {
        set endtime [clock seconds]
        set ttime [expr $endtime-$inittime]
        # WarnWinText "endtime:$endtime ttime:$ttime"
        WarnWinText "Write structural analysis concentrated loads (puntual force or moment): [::KUtils::Duration $ttime]"
    }
}

proc ::wkcf::WritePuntualLoads_m1 {AppId cloadtid kwordlist} {
    # ASBTRACT: Write concentrated loads (puntual force or moment)
    variable ndime; variable dprops
    variable sa_icondid

    # For all defined group identifier inside this load type
    foreach cgroupid $dprops($AppId,Loads,$cloadtid,AllGroupId) {
	# Get the load properties
	set GProps $dprops($AppId,Loads,$cloadtid,$cgroupid,GProps)
	# WarnWinText "cgroupid:$cgroupid GProps:$GProps"
	# Assign values
	foreach item $GProps {
	    foreach {cvar cval} $item {
		if {$cval eq "0.0"} {
		    set cval 0
		}
		set $cvar $cval
	    }
	} 
	# Declaración del diccionario
	set gpointforce [dict create]
	
	set usepointforce "No"
	# Diccionario gpointforce
	set f "%10i"
	set f [subst $f]
	dict set gpointforce $cgroupid $f
	
	switch -exact -- $ndime {
	    "2D" {
		set pointforcekword "PointForce2D"
		foreach item [list $Fx $Fy $Mx $My] {
		    if {$item !="0"} {
			set usepointforce "Yes"
			break
		    }
		}
	    }
	    "3D" {
		set pointforcekword "PointForce3D"
		foreach item [list $Fx $Fy $Fz $Mx $My $Mz] {
		    if {$item !="0"} {
			set usepointforce "Yes"
			break
		    }
		}
	    }
	}
	if {$usepointforce =="Yes"} {
	    # Active the PointForce condition
	    
	    write_calc_data puts "Begin Conditions $pointforcekword // GUI puntual load group identifier: $cgroupid"
	    foreach output [write_calc_data nodes -return -sorted $gpointforce] {
		incr sa_icondid
		set cf "[format "%4i %4i %10i" $sa_icondid 1 $output]"
		write_calc_data puts "$cf"
	    }
	    write_calc_data puts "End Conditions"
	    write_calc_data puts ""
	}
	
	
	set gforcex [dict create]
	set gforcey [dict create]
	if {$ndime eq "3D"} {
	    set gforcez [dict create]
	}
	
	# Assign values
	foreach item $GProps {
	    foreach {cvar cval} $item {
		if {$cval eq "0.0"} {
		    set cval 0
		}
		set $cvar $cval
	    }
	} 
	set IsFixed "0"
	set f "%10i [format "%2i %20.10f" $IsFixed $Fx]\n"
	set f [subst $f]
	dict set gforcex $cgroupid $f
	set f "%10i [format "%2i %20.10f" $IsFixed $Fy]\n"
	set f [subst $f]
	dict set gforcey $cgroupid $f
	if {$ndime eq "3D"} {
	    set f "%10i [format "%2i %20.10f" $IsFixed $Fz]\n"
	    set f [subst $f]
	    dict set gforcez $cgroupid $f
	}
	# WarnWinText "Fx:$Fx Fy:$Fy Fz:$Fz"
	
	# DISPLACEMENT_X or ROTATION_X
	if {$Fx ne 0} {
	    set xitem [lindex $kwordlist 0]
	    write_calc_data puts "Begin NodalData $xitem // GUI puntual load group identifier: $cgroupid"
	    write_calc_data nodes -sorted $gforcex
	    write_calc_data puts "End NodalData"
	    write_calc_data puts ""        
	}
	if {$Fy ne 0} {
	    set yitem [lindex $kwordlist 1]
	    write_calc_data puts "Begin NodalData $yitem // GUI puntual load group identifier: $cgroupid"
	    write_calc_data nodes -sorted $gforcey
	    write_calc_data puts "End NodalData"
	    write_calc_data puts ""        
	}
	if {$Fz ne 0} {
	    set zitem [lindex $kwordlist 2]
	    write_calc_data puts "Begin NodalData $zitem // GUI puntual load group identifier: $cgroupid"
	    write_calc_data nodes -sorted $gforcez
	    write_calc_data puts "End NodalData"
	    write_calc_data puts ""        
	}
	if { [info exists gforcex] } {
	    unset gforcex
	}
	if { [info exists gforcey] } {
	    unset gforcey
	}
	if { [info exists gforcez] } {
	    unset gforcez
	}
	if { [info exists gpointforce] } {
	    unset gpointforce
	}
    }
}

proc ::wkcf::WritePuntualLoads_m2 {AppId cloadtid kwordlist} {
    # ASBTRACT: Write concentrated loads (puntual force or moment)
    variable ndime; variable dprops
    variable filechannel; variable sa_icondid

    # For all defined group identifier inside this load type
    foreach cgroupid $dprops($AppId,Loads,$cloadtid,AllGroupId) {
        # Get the load properties
        set GProps $dprops($AppId,Loads,$cloadtid,$cgroupid,GProps)
        # WarnWinText "cgroupid:$cgroupid GProps:$GProps"
        # Assign values
        foreach item $GProps {
            foreach {cvar cval} $item {
                if {$cval eq "0.0"} {
                    set cval 0
                }
                set $cvar $cval
            }
        } 

        set usepointforce "No"
        
        switch -exact -- $ndime {
            "2D" {
                set pointforcekword "PointForce2D"
                foreach item [list $Fx $Fy $Mx $My] {
                    if {$item !="0"} {
                        set usepointforce "Yes"
                        break
                    }
                }
            }
            "3D" {
                set pointforcekword "PointForce3D"
                foreach item [list $Fx $Fy $Fz $Mx $My $Mz] {
                    if {$item !="0"} {
                        set usepointforce "Yes"
                        break
                    }
                }
            }
        }
        if {$usepointforce =="Yes"} {
            # Active the PointForce condition
            
            GiD_File fprintf $filechannel "%s" "Begin Conditions $pointforcekword // GUI puntual load group identifier: $cgroupid"
            foreach output [GiD_EntitiesGroups get $cgroupid nodes] {
                incr sa_icondid
                set cf "[format "%4i %4i %10i" $sa_icondid 1 $output]"
                GiD_File fprintf $filechannel "%s" "$cf"
            }
            GiD_File fprintf $filechannel "%s" "End Conditions"
            GiD_File fprintf $filechannel "%s" ""
        }

        # Assign values
        foreach item $GProps {
            foreach {cvar cval} $item {
                if {$cval eq "0.0"} {
                    set cval 0
                }
                set $cvar $cval
            }
        } 
        set IsFixed "0"
        # WarnWinText "Fx:$Fx Fy:$Fy Fz:$Fz"
        
        # DISPLACEMENT_X or ROTATION_X
        if {$Fx ne 0} {
            set xitem [lindex $kwordlist 0]
            GiD_File fprintf $filechannel "%s" "Begin NodalData $xitem // GUI puntual load group identifier: $cgroupid"
            foreach node_id [GiD_EntitiesGroups get $cgroupid nodes] {
                GiD_File fprintf $filechannel "%10i %2i %20.10f" $node_id $IsFixed $Fx
            }
            GiD_File fprintf $filechannel "%s" "End NodalData"
            GiD_File fprintf $filechannel ""        
        }
        if {$Fy ne 0} {
            set yitem [lindex $kwordlist 1]
            GiD_File fprintf $filechannel "%s" "Begin NodalData $yitem // GUI puntual load group identifier: $cgroupid"
            foreach node_id [GiD_EntitiesGroups get $cgroupid nodes] {
                GiD_File fprintf $filechannel "%10i %2i %20.10f" $node_id $IsFixed $Fy
            }
            GiD_File fprintf $filechannel "%s" "End NodalData"
            GiD_File fprintf $filechannel ""       
        }
        if {$Fz ne 0} {
            set zitem [lindex $kwordlist 2]
            GiD_File fprintf $filechannel "%s" "Begin NodalData $zitem // GUI puntual load group identifier: $cgroupid"
            foreach node_id [GiD_EntitiesGroups get $cgroupid nodes] {
                GiD_File fprintf $filechannel "%10i %2i %20.10f" $node_id $IsFixed $Fz
            }
            GiD_File fprintf $filechannel "%s" "End NodalData"
            GiD_File fprintf $filechannel "%s" ""        
        }
    }
}

proc ::wkcf::WritePuntualLoads_m0 {AppId cloadtid kwordlist} {
    # Write concentrated loads (puntual force or moment)
    variable ndime;    variable gidentitylist
    variable useqelem; variable dprops
    variable sa_icondid

    # For all defined group identifier inside this load type
    foreach cgroupid $dprops($AppId,Loads,$cloadtid,AllGroupId) {
	# Get the load properties
	set GProps $dprops($AppId,Loads,$cloadtid,$cgroupid,GProps)
	# WarnWinText "cgroupid:$cgroupid GProps:$GProps"
	# Assign values
	foreach item $GProps {
	    foreach {cvar cval} $item {
		set $cvar $cval
	    }
	} 
	# WarnWinText "Fx:$Fx Fy:$Fy Fz:$Fz"
	set allnlist [list]
	set GiDEntity "point"
	# Get all defined entities for this group identifier
	set callnlist [::KUtils::GetDefinedMeshGiDEntities $cgroupid $GiDEntity "Nodes" $useqelem]
	# WarnWinText "callnlist:$callnlist"
	if {[llength $callnlist]} {
	    lappend allnlist $callnlist
	} 

	# Sort the node list
	# WarnWinText "$GiDEntity alllist:$allnlist"
	
	if {[llength $allnlist]} {
	    set allnlist [lsort -integer -unique {*}$allnlist]         

	    set usepointforce "No"
	    switch -exact -- $ndime {
		"2D" {
		    set pointforcekword "PointForce2D"
		    foreach item [list $Fx $Fy $Mx $My] {
		        if {$item !="0.0"} {
		            set usepointforce "Yes"
		            break
		        }
		    }
		}
		"3D" {
		    set pointforcekword "PointForce3D"
		    foreach item [list $Fx $Fy $Fz $Mx $My $Mz] {
		        if {$item !="0.0"} {
		            set usepointforce "Yes"
		            break
		        }
		    }
		}
	    }


	    if {$usepointforce =="Yes"} {
		# Active the PointForce condition
		write_calc_data puts "Begin Conditions $pointforcekword // GUI puntual load group identifier: $cgroupid"
		foreach nodeid $allnlist {
		    incr sa_icondid
		    set cf "[format "%4i%4i%8i" $sa_icondid 1 $nodeid]"
		    write_calc_data puts "$cf"
		}
		write_calc_data puts "End Conditions"
		write_calc_data puts ""
	    }
	    
	    # WarnWinText "Fx:$Fx"
	    set nplist [list "0.0" ""]
	    set IsFixed "0"
	    # FORCE_X
	    if {$Fx ni $nplist} {
		set xitem [lindex $kwordlist 0]
		write_calc_data puts "Begin NodalData $xitem // GUI puntual load group identifier: $cgroupid"
		foreach nodeid $allnlist {
		    set cf "[format "%6i%4i%20.10f" $nodeid $IsFixed $Fx]"
		    write_calc_data puts "$cf"
		}
		write_calc_data puts "End NodalData"
		write_calc_data puts ""
	    }

	    # FORCE_Y
	    if {$Fy ni $nplist} {
		set xitem [lindex $kwordlist 1]
		write_calc_data puts "Begin NodalData $xitem // GUI puntual load group identifier: $cgroupid"
		foreach nodeid $allnlist {
		    set cf "[format "%6i%4i%20.10f" $nodeid $IsFixed $Fy]"
		    write_calc_data puts "$cf"
		}
		write_calc_data puts "End NodalData"
		write_calc_data puts ""
	    }
	    
	    if {$ndime =="3D"} {
		# FORCE_Z
		if {$Fz ni $nplist} {
		    set xitem [lindex $kwordlist 2]
		    write_calc_data puts "Begin NodalData $xitem // GUI puntual load group identifier: $cgroupid"
		    foreach nodeid $allnlist {
		        set cf "[format "%6i%4i%20.10f" $nodeid $IsFixed $Fz]"
		        write_calc_data puts "$cf"
		    }
		    write_calc_data puts "End NodalData"
		    write_calc_data puts ""
		}
	    }
	}
	# si allnlist no está vacío...        
	#else {
	#      msg "Check the Loads"
	#}
    }
}

proc ::wkcf::WriteBodyForceValues {props} {
    # Write the gravity properties to the kratos data file
    # Arguments
    # props => Body force properties
    variable ndime
    variable wmethod 
    variable filechannel 
    
    # WarnWinText "props:$props"
    set GravityValue [lindex $props 0 1]
    set Cx [expr [lindex $props 1 1] * $GravityValue]
    set Cy [expr [lindex $props 2 1] * $GravityValue]
    if {$ndime =="3D"} {
        set Cz [expr [lindex $props 3 1] * $GravityValue]        
        set vector "\($Cx, $Cy, $Cz\)"
    } else {
        set vector "\($Cx, $Cy, 0.0\)"
    }
    if {$wmethod eq 1} {
        write_calc_data puts " BODY_FORCE \[3\] $vector"
    } elseif {$wmethod eq 2} {
        GiD_File fprintf $filechannel "%s" " BODY_FORCE \[3\] $vector"
    }
}

proc ::wkcf::WriteStructuralProjectParameters {AppId fileid PDir} {
    # Write the project parameters file
    variable ndime; variable dprops

    # kratos key word xpath
    set kwxpath "Applications/$AppId"
    
    set domain_size 2
    if {$ndime =="2D"} {
	set domain_size 2
    } elseif {$ndime =="3D"} {
	set domain_size 3
    }
    
    # Domain size
    puts $fileid "domain_size = $domain_size"
    puts $fileid ""

    # Get others properties
    # Analysis type
    set cxpath "$AppId//c.AnalysisData//i.AnalysisType"
    set cproperty "dv"
    set AnalysisType [::xmlutils::setXml $cxpath $cproperty]
    # WarnWinText "AnalysisType:$AnalysisType"
    if {$AnalysisType =="Non-Linear"} {
	# Solution method
	set cxpath "$AppId//c.SolutionStrategy//c.Non-Linear//i.SolutionMethod"
	set SolutionMethod [::xmlutils::setXml $cxpath $cproperty]
	puts $fileid "Solution_method = \"$SolutionMethod\""
	# Convergence criteria
	set cxpath "$AppId//c.SolutionStrategy//c.Non-Linear//i.ConvergenceCriteria"
	set ConvergenceCriteria [::xmlutils::setXml $cxpath $cproperty]
	# Residual convergence tolerance
	set cxpath "$AppId//c.SolutionStrategy//c.Non-Linear//i.ResidualConvergenceTolerance"
	set ResidualConvergenceTolerance [::xmlutils::setXml $cxpath $cproperty]
	# Residual absolute tolerance
	set cxpath "$AppId//c.SolutionStrategy//c.Non-Linear//i.ResidualAbsoluteTolerance"
	set ResidualAbsoluteTolerance [::xmlutils::setXml $cxpath $cproperty]
	# Displacement convergence tolerance
	set cxpath "$AppId//c.SolutionStrategy//c.Non-Linear//i.DisplacementConvergenceTolerance"
	set DisplacementConvergenceTolerance [::xmlutils::setXml $cxpath $cproperty]
	# Displacement absolute tolerance
	set cxpath "$AppId//c.SolutionStrategy//c.Non-Linear//i.DisplacementAbsoluteTolerance"
	set DisplacementAbsoluteTolerance [::xmlutils::setXml $cxpath $cproperty]
	# Maximum iterations
	set cxpath "$AppId//c.SolutionStrategy//c.Non-Linear//i.MaximumIterations"
	set MaximumIterations [::xmlutils::setXml $cxpath $cproperty]

	switch -exact -- $ConvergenceCriteria {
	    "Displacement" {
		puts $fileid "Displacement_Convergence_Tolerance = $DisplacementConvergenceTolerance"
		puts $fileid "Displacement_Absolute_Tolerance = $DisplacementAbsoluteTolerance"
		puts $fileid "Residual_Convergence_Tolerance = $ResidualConvergenceTolerance"
		puts $fileid "Residual_Absolute_Tolerance = $ResidualAbsoluteTolerance"
	    }
	    "Residual" {
		puts $fileid "Residual_Convergence_Tolerance = $ResidualConvergenceTolerance"
		puts $fileid "Residual_Absolute_Tolerance = $ResidualAbsoluteTolerance"
		puts $fileid "Displacement_Convergence_Tolerance = $DisplacementConvergenceTolerance"
		puts $fileid "Displacement_Absolute_Tolerance = $DisplacementAbsoluteTolerance"
	    }
	    "DisplacementAndResidual" {
		puts $fileid "Residual_Convergence_Tolerance = $ResidualConvergenceTolerance"
		puts $fileid "Residual_Absolute_Tolerance = $ResidualAbsoluteTolerance"
		puts $fileid "Displacement_Convergence_Tolerance = $DisplacementConvergenceTolerance"
		puts $fileid "Displacement_Absolute_Tolerance = $DisplacementAbsoluteTolerance"
	    }
	    "DisplacementOrResidual" {
		puts $fileid "Residual_Convergence_Tolerance = $ResidualConvergenceTolerance"
		puts $fileid "Residual_Absolute_Tolerance = $ResidualAbsoluteTolerance"
		puts $fileid "Displacement_Convergence_Tolerance = $DisplacementConvergenceTolerance"
		puts $fileid "Displacement_Absolute_Tolerance = $DisplacementAbsoluteTolerance"
	    }
	}
	puts $fileid "Max_Iter = $MaximumIterations"
	set cConvergenceCriteria "Displacement_Criteria"
	set cConvergenceCriteria [::xmlutils::getKKWord $kwxpath $ConvergenceCriteria]
	puts $fileid "Convergence_Criteria = \"$cConvergenceCriteria\""
	# WarnWinText "SolutionMethod:$SolutionMethod ConvergenceCriteria:$ConvergenceCriteria ConvergenceTolerance:$ConvergenceTolerance AbsoluteTolerance:$AbsoluteTolerance MaximumIterations:$MaximumIterations"
	#END NON LINEAR
    } elseif {$AnalysisType =="Linear"} {

	puts $fileid "Convergence_Tolerance = 1.0"
	puts $fileid "Absolute_Tolerance = 1.0"
	puts $fileid "Max_Iter = 1"
	set cConvergenceCriteria "Displacement_Criteria"
	puts $fileid "Convergence_Criteria = \"$cConvergenceCriteria\""
	#Hay que ponerlos para que Kratos no se queje, pero luego no lo va a usar
	puts $fileid ""
	puts $fileid "Residual_Convergence_Tolerance = 1.0E-3"
	puts $fileid "Residual_Absolute_Tolerance = 1.0E-6"
	puts $fileid "Displacement_Convergence_Tolerance = 1.0E-6"
	puts $fileid "Displacement_Absolute_Tolerance = 1.0E-9"
    }
    # Solution method
    set cxpath "$AppId//c.SolutionStrategy//c.Non-Linear//i.SolutionMethod"
    set SolutionMethod [::xmlutils::setXml $cxpath $cproperty]
    puts $fileid "Solution_method = \"$SolutionMethod\""
    # Solution type
    set cxpath "$AppId//c.AnalysisData//i.SolutionType"
    set SolutionType [::xmlutils::setXml $cxpath $cproperty]
    # WarnWinText "SolutionType:$SolutionType"
    if {$SolutionType =="Dynamic"} {
	# Start time
	set cxpath "$AppId//c.SolutionStrategy//c.Dynamic//i.StartTime"
	set StartTime [::xmlutils::setXml $cxpath $cproperty]
	# End time
	set cxpath "$AppId//c.SolutionStrategy//c.Dynamic//i.EndTime"
	set EndTime [::xmlutils::setXml $cxpath $cproperty]
	# Delta time
	set cxpath "$AppId//c.SolutionStrategy//c.Dynamic//i.DeltaTime"
	set DeltaTime [::xmlutils::setXml $cxpath $cproperty]
	
	# WarnWinText "StartTime:$StartTime EndTime:$EndTime DeltaTime:$DeltaTime"
	puts $fileid "Dt = $StartTime"
	puts $fileid "max_time = $EndTime"
	set nsteps [expr int(double($EndTime-$StartTime)/double($DeltaTime))]
	puts $fileid "nsteps = $nsteps"

	# Output step 
	set cxpath "$AppId//c.Results//i.OutputDeltaTime"
	set OutputDeltaTime [::xmlutils::setXml $cxpath $cproperty]
	puts $fileid "output_time = $OutputDeltaTime"
	# WarnWinText "OutputDeltaTime:$OutputDeltaTime"
	set output_step [expr int($OutputDeltaTime/double($DeltaTime))]
	# WarnWinText "output_step:$output_step"
	puts $fileid "output_step = $output_step"
	
	puts $fileid "SolverType = \"DynamicSolver\""

    } elseif {$SolutionType =="RelaxedDynamic"} {

	# Start time
	set cxpath "$AppId//c.SolutionStrategy//c.Dynamic//i.StartTime"
	set StartTime [::xmlutils::setXml $cxpath $cproperty]
	# End time
	set cxpath "$AppId//c.SolutionStrategy//c.Dynamic//i.EndTime"
	set EndTime [::xmlutils::setXml $cxpath $cproperty]
	# Delta time
	set cxpath "$AppId//c.SolutionStrategy//c.Dynamic//i.DeltaTime"
	set DeltaTime [::xmlutils::setXml $cxpath $cproperty]

	# WarnWinText "StartTime:$StartTime EndTime:$EndTime DeltaTime:$DeltaTime"
	puts $fileid "Dt = $StartTime"
	puts $fileid "max_time = $EndTime"
	set nsteps [expr int(double($EndTime-$StartTime)/double($DeltaTime))]
	puts $fileid "nsteps = $nsteps"

	# Output step 
	set cxpath "$AppId//c.Results//i.OutputDeltaTime"
	set OutputDeltaTime [::xmlutils::setXml $cxpath $cproperty]
	puts $fileid "output_time = $OutputDeltaTime"
	# WarnWinText "OutputDeltaTime:$OutputDeltaTime"
	set output_step [expr int($OutputDeltaTime/double($DeltaTime))]
	# WarnWinText "output_step:$output_step"
	puts $fileid "output_step = $output_step"
	
	puts $fileid "SolverType = \"RelaxedDynamicSolver\""

    } elseif {$SolutionType =="Static"} {
	puts $fileid "SolverType = \"StaticSolver\""
    }


    # Solution strategy
    set cxpath "$AppId//c.SolutionStrategy//i.LinearSolverType"
    set LinearSolverType [::xmlutils::setXml $cxpath $cproperty]
    if {$LinearSolverType =="Direct"} {
	# Direct solver type
	set cxpath "$AppId//c.SolutionStrategy//i.DirectSolverType"
	set DirectSolverType [::xmlutils::setXml $cxpath $cproperty]
	# WarnWinText "DirectSolverType:$DirectSolverType"
	set cDirectSolverType "SkylineLUFactorization"
	set cDirectSolverType [::xmlutils::getKKWord $kwxpath $DirectSolverType]
	puts $fileid "LinearSolver = \"$cDirectSolverType\""

    } elseif {$LinearSolverType =="Iterative"} {

	# Iterative solver type        
	set cxpath "$AppId//c.SolutionStrategy//i.IterativeSolverType"
	set IterativeSolverType [::xmlutils::setXml $cxpath $cproperty]
	# Tolerance
	set cxpath "$AppId//c.SolutionStrategy//i.Tolerance"
	set Tolerance [::xmlutils::setXml $cxpath $cproperty]
	# Maximum iteration
	set cxpath "$AppId//c.SolutionStrategy//i.MaximumIteration"
	set MaximumIteration [::xmlutils::setXml $cxpath $cproperty]
	# preconditioner type
	set cxpath "$AppId//c.SolutionStrategy//i.PreconditionerType"
	set PreconditionerType [::xmlutils::setXml $cxpath $cproperty]
	# WarnWinText "IterativeSolverType:$IterativeSolverType Tolerance:$Tolerance MaximumIteration:$MaximumIteration PreconditionerType:$PreconditionerType"

	# Solver type
	set lsolver [::xmlutils::getKKWord $kwxpath $IterativeSolverType]
	set precond [::xmlutils::getKKWord $kwxpath $PreconditionerType]
	puts $fileid "LinearSolver = \"${lsolver}_${precond}\""
	puts $fileid "Linear_Solver_Tolerance = $Tolerance"
	puts $fileid "Linear_Solver_Max_Iteration = $MaximumIteration"

    }        

    puts $fileid "FindNodalNeighbours = \"False\""
    
    # Check for use shell elements
    set usenbst "No"
    set useshells "No"
    set usebeams "No"

    set shelllist [list "ShellIsotropic" "ShellAnisotropic" "EBST"]
    set beamlist [list "BeamElement"]

    if {([info exists dprops($AppId,AllKElemId)]) && ($dprops($AppId,AllKElemId)>0)} {
	# For all defined kratos elements        
	foreach celemid $dprops($AppId,AllKElemId) {
	    if {$celemid in $shelllist} {
		set useshells "Yes"
		if {$celemid == "EBST"} {
		    set usenbst "Yes"
		    break
		}
	    } elseif {$celemid in $beamlist} {
		set usebeams "Yes"
	    }
	}
    }
    if {$usenbst =="Yes"} {
	puts $fileid "FindElementalNeighbours = \"True\""
    } else {
	puts $fileid "FindElementalNeighbours = \"False\""
    }
    if {($useshells eq "Yes")||($usebeams eq "Yes")} {
	puts $fileid "Rotational_Dofs = \"True\""
    } else {
	puts $fileid "Rotational_Dofs = \"False\""
    }

    # For results
    puts $fileid ""
    # On nodes results
    set cnrlist [list "Displacements" "Forces" "Reactions"]
    set nodal_results "nodal_results=\["
    foreach cnr $cnrlist {
	set cxpath "$AppId//c.Results//c.OnNodes//i.${cnr}"
	set cproperty "dv"
	set cvalue [::xmlutils::setXml $cxpath $cproperty]
	if {$cvalue =="Yes"} {
	    set cnkr [::xmlutils::getKKWord $kwxpath $cnr]
	    append nodal_results "\"$cnkr\","
	}
    }
    set findcomma [string last "," $nodal_results]
    if {$findcomma !="-1"} {
	set nodal_results [string range $nodal_results 0 end-1]
	append nodal_results "\]" 
	puts $fileid "$nodal_results"
    }
    # WarnWinText "nodal_results:$nodal_results"

    # On Gauss point results
    #puts $fileid ""
    set cgrlist [list "GreenLagrangeStrainTensor" "Rotations" "PK2StressTensor" "BeamMoments" "BeamForces"]
    set gauss_points_results "gauss_points_results=\["
    foreach cgr $cgrlist {
	set cxpath "$AppId//c.Results//c.OnGaussPoints//i.${cgr}"
	set cproperty "dv"
	set cvalue [::xmlutils::setXml $cxpath $cproperty]
	if {$cvalue =="Yes"} {
	    set cgkr [::xmlutils::getKKWord $kwxpath $cgr]
	    append gauss_points_results "\"$cgkr\","
	}
    }
    set findcomma [string last "," $gauss_points_results]
    if {$findcomma !="-1"} {
	set gauss_points_results [string range $gauss_points_results 0 end-1]
	append gauss_points_results "\]" 
	puts $fileid "$gauss_points_results"
    }
    # WarnWinText "gauss_points_results:$gauss_points_results"

    # GiD post mode variables
    ::wkcf::WriteGiDPostMode $AppId $fileid 

    puts $fileid ""
    set PName [::KUtils::GetPaths "PName"]
    puts $fileid "problem_name=\"${PName}${AppId}\"" 
    puts $fileid "problem_path=\"$PDir\"" 
    
    
    # Commented by J. Garate on 17/12/2012
    # Get the kratos path 
    # set cxpath "GeneralApplicationData//c.ProjectConfiguration//i.KratosPath"
    # set cproperty "dv"
    # set KratosPath [::xmlutils::setXml $cxpath $cproperty]
    # set KratosPath [file native $KratosPath]

    # Write the kratos path
    # puts $fileid "kratos_path=\"${KratosPath}\"" 
	
}

proc ::wkcf::WriteConstitutiveLawsProperties {} {
    # Write constitutive laws properties
    variable dprops;  variable ActiveAppList
    
    set AppId "StructuralAnalysis"
    set filename "materials.py"
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
    # puts $fileid "import os, sys"
    puts $fileid ""
    puts $fileid "# Importing the Kratos Library"
    puts $fileid "from KratosMultiphysics import *"
    puts $fileid "from KratosMultiphysics.StructuralApplication import *"
    # Cross section properties
    puts $fileid "from apply_sections import SetProperties"
  

    puts $fileid "def AssignMaterial(Properties):"
    
    # For each properties
    set propid 0
    foreach PropertyId $dprops($AppId,GKProps,AllPropertyId) {
	incr propid 1
	puts $fileid "# GUI property identifier: $PropertyId"
	# Get the material identifier for this property 
	set MatId $dprops($AppId,Property,$PropertyId,MatId)
	puts $fileid "# GUI material identifier: $MatId"
	# Check for use fluency
	if {$dprops($AppId,Material,$MatId,UseFluency) =="Yes"} {
	    puts $fileid "    fluency = $dprops($AppId,Material,$MatId,Fluency)"
	}
	# Check for use behavior
	if {$dprops($AppId,Material,$MatId,UseBehavior) =="Yes"} {
	    puts $fileid "    behavior = $dprops($AppId,Material,$MatId,Behavior)"
	}
	# Write the property identifier
	puts $fileid "    prop_id = $propid\;"
	puts $fileid "    prop = Properties\[prop_id\]"
	# Write material model
	puts $fileid "    mat = $dprops($AppId,Material,$MatId,MatModel)\;"
	puts $fileid "    prop.SetValue(CONSTITUTIVE_LAW, mat.Clone())\;"

	# Write cross section properties
	# Get all properties
	set CrossSectionProps [list]
	if {[info exists dprops($AppId,CrossSection,$PropertyId,CProps)]} { 
	    set CrossSectionProps $dprops($AppId,CrossSection,$PropertyId,CProps)
	}
	
	# Select the section type
	if {([info exists dprops($AppId,CrossSection,$PropertyId,SectionType)]) && ([llength $CrossSectionProps])} {
	    set SectionType $dprops($AppId,CrossSection,$PropertyId,SectionType)
	    # Rectangular cross section
	    if {$SectionType eq "Rectangular"} {
		set csectiontype "\"Square\""
		set proplist [list 1 1]
		# Get the property list
		foreach cprop $CrossSectionProps {
		    lassign $cprop cspropid csproptypeid cvalue _SectionType
		    if {$cspropid eq "RectangularHeight"} {
			lset proplist 0 $cvalue
		    } elseif {$cspropid eq "RectangularWidth"} {
			lset proplist 1 $cvalue
		    }
		}
		# wa "proplist:$proplist"
		set proplist [join $proplist ","]
		set endprop "\[$proplist\]"
		# wa "endprop:$endprop"
		puts $fileid "    prop = SetProperties($csectiontype,$endprop,prop)\;" 

	    } elseif {$SectionType eq "Circular"} {
		set csectiontype "\"$SectionType\""
		set proplist [list 1]
		# Get the property list
		foreach cprop $CrossSectionProps {
		    lassign $cprop cspropid csproptypeid cvalue _SectionType
		    if {$cspropid eq "CircularDiameter"} {
			lset proplist 0 $cvalue
		    } 
		}
		# wa "proplist:$proplist"
		set proplist [join $proplist ","]
		set endprop "\[$proplist\]"
		# wa "endprop:$endprop"
		puts $fileid "    prop = SetProperties($csectiontype,$endprop,prop)\;" 
	    } else {
		# Profile sections
		set csectiontype "\"$SectionType\""
		set proplist [list 1]
		# Get the property list
		foreach cprop $CrossSectionProps {
		    lassign $cprop cspropid csproptypeid cvalue _SectionType
		    if {$cspropid eq "ProfileDB"} {
			lset proplist 0 $cvalue
		    } 
		}
		# wa "proplist:$proplist"
		set proplist [join $proplist ","]
		set endprop "\[\"$proplist\"\]"
		# wa "endprop:$endprop"
		puts $fileid "    prop = SetProperties($csectiontype,$endprop,prop)\;" 
	    }
	}
	
    }

    close $fileid

}
 
proc ::wkcf::GetBehaviorFluencyProperties {AppId MatId MatModel cptype ptype} {
    # Get the behavior and fluency material properties
    variable dprops
    variable ndime

    # WarnWinText "AppId:$AppId MatId:$MatId cptype:$cptype ptype:$ptype"
    # Xpath for constitutive laws
    set clxpath "CLawProperties"
    # Get all material properties
    set mpxpath "[::KMat::findMaterialParent $MatId]//m.${MatId}"
    # WarnWinText "mpxpath:$mpxpath"

    # Set fluency and behavior variables
    set dprops($AppId,Material,$MatId,UseFluency) "Yes"
    set dprops($AppId,Material,$MatId,UseBehavior) "Yes"
    # Get the softening behavior
    set mbehavior [::xmlutils::getKKWord $clxpath $ptype "mbehavior"]
    # Get the softening behavior xpath values
    set mbxpath [::xmlutils::getKKWord $clxpath $ptype "mbxpath"]
    # WarnWinText "mbehavior:$mbehavior mbxpath:$mbxpath"
    # Get the current behavior 
    set cbvalue [lindex [::KMat::getMaterialProperties "p" "$mpxpath//$mbxpath//p.$mbehavior"] 0 1]
    # Get the internal behavior properties
    set mbivalues [split [::xmlutils::getKKWord $clxpath $ptype "mbivalues"] ,]
    # Get the write behavior properties
    set mbwritev [split [::xmlutils::getKKWord $clxpath $ptype "mbwritev"] ,]
    # WarnWinText "mbwritev:$mbwritev mbivalues:$mbivalues\n$mpxpath//$mbxpath//p.$mbehavior cbvalue:$cbvalue"
    foreach mbiv $mbivalues mbwv $mbwritev {
        if {$mbiv ==$cbvalue} {
            set dprops($AppId,Material,$MatId,Behavior) "$mbwv"
            break
        }
    }
    # WarnWinText "dprops($AppId,Material,$MatId,Behavior):$dprops($AppId,Material,$MatId,Behavior)"
    if {$MatModel =="Damage"} {
        # Damage models
        # Get the energy yield function
        # Get the internal state properties
        set msivalues [split [::xmlutils::getKKWord $clxpath "MState" "msivalues"] ,]
        # Get the write behavior properties
        set mswritev [split [::xmlutils::getKKWord $clxpath "MState" "mswritev"] ,]
        # WarnWinText "mswritev:$mswritev msivalues:$msivalues"
        foreach msiv $msivalues mswv $mswritev {
            # WarnWinText "msiv:$msiv cptype:$cptype mswv:$mswv" 
            if {$msiv ==$cptype} {
                set dprops($AppId,Material,$MatId,Fluency) "EnergyYieldFunction(State.${mswv})"
                break
            }
        }
    } elseif {$MatModel == "Elasto-Plastic"} {
	# Elasto-plastic models
	# Get the internal state properties
	set msivalues [split [::xmlutils::getKKWord $clxpath "MState" "msivalues"] ,]
	# Get the write behavior properties
	set mswritev [split [::xmlutils::getKKWord $clxpath "MState" "mswritev"] ,]
	# WarnWinText "mswritev:$mswritev msivalues:$msivalues"
	set cstate ""
	foreach msiv $msivalues mswv $mswritev {
	    # WarnWinText "msiv:$msiv cptype:$cptype mswv:$mswv" 
	    if {$msiv ==$cptype} {
		set cstate "myState.${mswv}"
		break
	    }
	}
	# WarnWinText "cstate:$cstate"
	# Get the yield function properties
	# Get the yield criteria
	set yfid "YieldFunctions"
	set myieldcriteria [::xmlutils::getKKWord "$clxpath" $ptype "myieldcriteria"]
	# Get the yield criteria xpath values
	set mycxpath [::xmlutils::getKKWord "$clxpath" $ptype "mycxpath"]
	# Get the current yield criteria
	set cycvalue [lindex [::KMat::getMaterialProperties "p" "$mpxpath//$mycxpath//p.$myieldcriteria"] 0 1]
	# WarnWinText "myieldcriteria:$myieldcriteria mycxpath:$mycxpath cycvalue:$cycvalue"
	# Get the yield function options
	set yfivalues [split [::xmlutils::getKKWord "$clxpath//$yfid" "AvailableYieldFunction" "yfivalues"] ,]
	# Get the write yield function properties
	set yfwritev [split [::xmlutils::getKKWord "$clxpath//$yfid" "AvailableYieldFunction" "yfwritev"] ,]
	# WarnWinText "yfwritev:$yfwritev yfivalues:$yfivalues"
	set cyf ""
	foreach yfiv $yfivalues yfwv $yfwritev {
	    # WarnWinText "yfiv:$yfiv cycvalue:$cycvalue yfwv:$yfwv" 
	    if {$yfiv ==$cycvalue} {
            set cyf "$yfwv"
            break
	    }
	}
	# WarnWinText "cyf:$cyf"
	
	set dprops($AppId,Material,$MatId,Fluency) "${cyf}(${cstate},PotencialPlastic.Associated)"
    }
    # WarnWinText "dprops($AppId,Material,$MatId,Fluency):$dprops($AppId,Material,$MatId,Fluency)"
}

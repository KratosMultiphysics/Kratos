###############################################################################
#
#    NAME: wkcfstructuralanalysis.tcl
#
#    PURPOSE: Useful procedures to work with the structural application
#
#    QUANTECH ATZ-DEVELOPMENT DEPARTMENT
#
#    AUTHOR : G. Socorro
#
#    CREATED AT: 02/04/12
#
#    HISTORY:
#
#     2.0- 10/10/13-G. Socorro, modify the GetBehaviorFluencyProperties to get the new keyword from kratos_key_words.xml file
#     1.9- 20/09/13-G. Socorro, modify the proc WriteSurfaceLoad to correct a bug when write linear and quadratic quadrilateral elements
#     1.8- 19/09/13-G. Socorro, modify the proc WriteVolumeAcceleration to write VOLUME_ACCELERATION_* nodal properties
#     1.7- 25/06/13-A. Melendo, new proc WriteFaceLoads
#     1.6- 17/06/13-G. Socorro, delete wmethod variable and all relalated procedure (*_m0,*_m1,*_m2) => now we are using only the new GiD groups
#     1.5- 22/05/13-G. Socorro, correct a bug in the procs WritePointLoads and WriteBeamLoads when write condition
#                               identifier in the case of more than one load group
#     1.4- 17/05/13-G. Socorro, add the proc WriteBeamUniformlyDistributedLoads, update the proc WriteConstitutiveLawsProperties
#                               to write the cross section properties
#     1.3- 24/04/13-G. Socorro, write Rotational_Dofs = True for beam element type
#     1.2- 18/03/13-G. Socorro, correct a bug in the proc WritePressureLoads_m2 change $group by $cgroupid
#     1.1- 17/12/12-J. Garate,  Kratos Path is no longer written at ProjectParameters.py
#     1.0- 12/11/12-J. Garate,  Fixed some errors
#     0.9- 07/11/12-J. Garate,  Modification and adaptation on functions: WriteDispRotBC, WritePressureLoads, WritePointLoads, WriteBodyForceValues
#                               Creation of functions using GiD_File fprintf $filechannel "%s" format
#     0.8- 09/10/12-G. Socorro, correct a bug with the output format in the proc WritePressureLoads_m1
#     0.7- 14/05/12-G. Socorro, modify the proc WriteDispRotBC_m1 to write many groups with displacement/rotation bc
#     0.6- 06/05/12-G. Socorro, update the proc WritePressureLoads to write using the fast method (write_calc_data)
#     0.5- 05/05/12-G. Socorro, improve the proc WriteDispRotBC to write using the fast method (write_calc_data)
#     0.4- 20/04/12-J. Garate, Acabada ::wkcf::WritePointLoads_m1, empezada Pressure
#     0.3- 19/04/12-J. Garate, ::wkcf::WritePointLoads_m1
#     0.2- 16/04/12-G. Socorro, ::wkcf::WriteDispRotBC_m1
#     0.1- 02/04/12-G. Socorro, create a base source code from wkcf.tcl
#
###############################################################################

proc ::wkcf::WriteDispRotBC {AppId ccondid kwordlist} {
    # ABSTRACT: Write displacements or rotation boundary condition
    variable ndime; variable dprops
    variable filechannel

    # For debug
    if {!$::wkcf::pflag} {
	set inittime [clock seconds]
    }

    # For all defined group identifier inside this condition type
    foreach cgroupid $dprops($AppId,BC,$ccondid,AllGroupId) {
	# Get the condition properties
	lassign $dprops($AppId,BC,$ccondid,$cgroupid,GProps) fix_x xval fix_y yval fix_z zval
	# WarnWinText "fix_x:$fix_x xval:$xval fix_y:$fix_y yval:$yval fix_z:$fix_z zval:$zval"

	# DISPLACEMENT_X or ROTATION_X
	set nodes [GiD_EntitiesGroups get $cgroupid nodes]

	if {[llength $nodes]} {
	    set xitem [lindex $kwordlist 0]
      if {$fix_x != "0" || $xval != "0.000000"} {
	    GiD_File fprintf $filechannel "%s" "Begin NodalData $xitem"
	    foreach nodeid $nodes {
		GiD_File fprintf $filechannel "%10i %4i %10f" $nodeid $fix_x $xval
	    }
	    GiD_File fprintf $filechannel "End NodalData"
	    GiD_File fprintf $filechannel ""
	}
  }

	# DISPLACEMENT_Y or ROTATION_Y
	if {[llength $nodes]} {
	    set yitem [lindex $kwordlist 1]
      if {$fix_y != "0" || $yval != "0.000000"} {
	    GiD_File fprintf $filechannel "%s" "Begin NodalData $yitem"
	    foreach nodeid $nodes {
		GiD_File fprintf $filechannel "%10i %4i %10f" $nodeid $fix_y $yval
	    }
	    GiD_File fprintf $filechannel "End NodalData"
	    GiD_File fprintf $filechannel ""
	}
  }

	# DISPLACEMENT_Z or ROTATION_Z
	if {($ndime eq "3D") && ([llength $nodes])} {
	    set zitem [lindex $kwordlist 2]
      if {$fix_z != "0" || $zval != "0.000000"} {
	    GiD_File fprintf $filechannel "%s" "Begin NodalData $zitem"
	    foreach nodeid $nodes {
		GiD_File fprintf $filechannel "%10i %4i %10f" $nodeid $fix_z $zval
	    }
	    GiD_File fprintf $filechannel "End NodalData"
	    GiD_File fprintf $filechannel ""
	}
    }
    }

    # For debug
    if {!$::wkcf::pflag} {
	set endtime [clock seconds]
	set ttime [expr $endtime-$inittime]
	# WarnWinText "endtime:$endtime ttime:$ttime"
	WarnWinText "Write structural analysis displacements or rotation boundary condition: [::KUtils::Duration $ttime]"

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
	    # wa "cloadtid:$cloadtid"
	    # Check for all defined group identifier inside this load type
	    if {([info exists dprops($AppId,Loads,$cloadtid,AllGroupId)]) && ([llength $dprops($AppId,Loads,$cloadtid,AllGroupId)])} {
	    # Select the load type
	    switch -exact -- $cloadtid {
		"VolumeAcceleration"
		{
		    # Concentrate or puntual loads (forces or moments)
		    set kwordlist [list [::xmlutils::getKKWord $kwxpath "Gx"] [::xmlutils::getKKWord $kwxpath "Gy"] [::xmlutils::getKKWord $kwxpath "Gz"]]
		    # Write puntual loads
		    ::wkcf::WriteVolumeAcceleration $AppId $cloadtid $kwordlist
		}
		"PointLoad2D"
		{
		    # Concentrate or puntual loads (forces or moments)
		    set kwordlist [list [::xmlutils::getKKWord $kwxpath "Fx"] [::xmlutils::getKKWord $kwxpath "Fy"] [::xmlutils::getKKWord $kwxpath "Fz"]]
		    # Write puntual loads
		    ::wkcf::WritePointLoads $AppId $cloadtid $kwordlist
		}
		"AxisymPointLoad"
		{
		    # Concentrate or puntual loads (forces or moments)
		    set kwordlist [list [::xmlutils::getKKWord $kwxpath "Fx"] [::xmlutils::getKKWord $kwxpath "Fy"] [::xmlutils::getKKWord $kwxpath "Fz"]]
		    # Write puntual loads
		    ::wkcf::WritePointLoadsAxisym $AppId $cloadtid $kwordlist
		}
		"PointLoad3D"
		{
		    # Concentrate or puntual loads (forces or moments)
		    set kwordlist [list [::xmlutils::getKKWord $kwxpath "Fx"] [::xmlutils::getKKWord $kwxpath "Fy"] [::xmlutils::getKKWord $kwxpath "Fz"]]
		    # Write puntual loads
		    ::wkcf::WritePointLoads $AppId $cloadtid $kwordlist
		}
		"PointMoment"
		{
		    # Concentrate or puntual loads (forces or moments)
		    set kwordlist [list [::xmlutils::getKKWord $kwxpath "Mx"] [::xmlutils::getKKWord $kwxpath "My"] [::xmlutils::getKKWord $kwxpath "Mz"]]
		    # Write puntual loads
		    ::wkcf::WritePointMoments $AppId $cloadtid $kwordlist
		}
		"LineLoad"
		{
		    # Write line loads in 2D
		    ::wkcf::WriteLineLoad $AppId $cloadtid
		}
		"AxisymLineLoad"
		{
		    # Write line loads in 2D beams
		    ::wkcf::WriteLineLoadAxisym $AppId $cloadtid
		}
		"BeamLoad3D"
		{
		    # Write line loads in 3D beams
		    ::wkcf::WriteLineLoad $AppId $cloadtid
		}
		"LinePressureLoad" {
		    # Pressure loads
		    ::wkcf::WriteLinePressureLoad $AppId $cloadtid
		}
		"AxisymLinePressureLoad" {
		    # Pressure loads
		    ::wkcf::WriteLinePressureLoadAxisym $AppId $cloadtid
		}
		"BeamPressure3D" {
		    # Pressure loads
		    ::wkcf::WriteLinePressureLoad $AppId $cloadtid
		}
		"SurfaceLoad2D" {
		    # Surface Face loads
		    ::wkcf::WriteSurfaceLoad $AppId $cloadtid
		}
		"SurfaceLoad3D" {
		    # Surface Face loads
		    ::wkcf::WriteSurfaceLoad $AppId $cloadtid
		}
		"SurfacePressureLoad2D" {
		    # Pressure loads
		    ::wkcf::WriteSurfacePressure $AppId $cloadtid
		}
		"SurfacePressureLoad3D" {
		    # Pressure loads
		    ::wkcf::WriteSurfacePressure $AppId $cloadtid
		}
	    }
	    }
	}
    }
}


proc ::wkcf::WriteVolumeAcceleration {AppId cloadtid kwordlist} {
    # ABSTRACT: Write volume acceleration to kratos data file
    variable ndime; variable dprops
    variable filechannel

    # For debug
    if {!$::wkcf::pflag} {
	set inittime [clock seconds]
    }

    # Check for all defined group identifier inside this load type
    if {([info exists dprops($AppId,Loads,$cloadtid,AllGroupId)]) && ([llength $dprops($AppId,Loads,$cloadtid,AllGroupId)])} {

	set fix_x 1
	set fix_y 1
	set fix_z 1

	foreach cgroupid $dprops($AppId,Loads,$cloadtid,AllGroupId) {

	    # Get the group properties
	    set cprop $dprops($AppId,Loads,$cloadtid,$cgroupid,GProps)
	    # wa "cprop:$cprop"
	    set GravityValue [lindex $cprop 0 1]
	    set Cx [expr [lindex $cprop 1 1] * $GravityValue]
	    set Cy [expr [lindex $cprop 2 1] * $GravityValue]
	    if {$ndime =="3D"} {
		set Cz [expr [lindex $cprop 3 1] * $GravityValue]
	    } else {
		set Cz 0
	    }

	    # Get the group entities list
	    set nodes [GiD_EntitiesGroups get $cgroupid nodes]

	    # VOLUME_ACCELERATION_X
	    if {[llength $nodes]} {
		set xitem [lindex $kwordlist 0]
		GiD_File fprintf $filechannel "%s" "Begin NodalData $xitem"
		foreach nodeid $nodes {
		    GiD_File fprintf $filechannel "%10i %4i %10f" $nodeid $fix_x $Cx
		}
		GiD_File fprintf $filechannel "End NodalData"
		GiD_File fprintf $filechannel ""


		# VOLUME_ACCELERATION_Y
		set yitem [lindex $kwordlist 1]
		GiD_File fprintf $filechannel "%s" "Begin NodalData $yitem"
		foreach nodeid $nodes {
		    GiD_File fprintf $filechannel "%10i %4i %10f" $nodeid $fix_y $Cy
		}
		GiD_File fprintf $filechannel "End NodalData"
		GiD_File fprintf $filechannel ""


		# VOLUME_ACCELERATION_Z
		if {$ndime eq "3D"} {
		    set zitem [lindex $kwordlist 2]
		    GiD_File fprintf $filechannel "%s" "Begin NodalData $zitem"
		    foreach nodeid $nodes {
		        GiD_File fprintf $filechannel "%10i %4i %10f" $nodeid $fix_z $Cz
		    }
		    GiD_File fprintf $filechannel "End NodalData"
		    GiD_File fprintf $filechannel ""
		}
	    }
	}
    }
}

proc ::wkcf::WritePointLoads {AppId cloadtid kwordlist} {
    # ABSTRACT: Write concentrated loads (puntual force or moment)
    variable ndime; variable dprops
    variable filechannel; variable sa_icondid

    # For debug
    if {!$::wkcf::pflag} {
	set inittime [clock seconds]
    }

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
    set pointforcekword "PointLoadCondition2D1N"
		foreach item [list $Fx $Fy $Mz] {
		    if {$item !="0"} {
		        set usepointforce "Yes"
		        break
		    }
		}
	    }
	    "3D" {
    set pointforcekword "PointLoadCondition3D1N"
		foreach item [list $Fx $Fy $Fz $Mx $My $Mz] {
		    if {$item !="0"} {
		        set usepointforce "Yes"
		        break
		    }
		}
	    }
	}
	if {$usepointforce =="Yes"} {
	    # Active the PointLoad condition

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

    # For debug
    if {!$::wkcf::pflag} {
	set endtime [clock seconds]
	set ttime [expr $endtime-$inittime]
	# WarnWinText "endtime:$endtime ttime:$ttime"
	WarnWinText "Write structural analysis concentrated loads (puntual force or moment): [::KUtils::Duration $ttime]"
    }
}

proc ::wkcf::WritePointMoments {AppId cloadtid kwordlist} {
    # ABSTRACT: Write concentrated loads (puntual force or moment)
    variable ndime; variable dprops
    variable filechannel; variable sa_icondid

    # For debug
    if {!$::wkcf::pflag} {
	set inittime [clock seconds]
    }

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
    set pointforcekword "PointTorqueCondition2D1N"
		foreach item [list $Mz] {
		    if {$item !="0"} {
		        set usepointforce "Yes"
		        break
		    }
		}
	    }
	    "3D" {
    set pointforcekword "PointTorqueCondition3D1N"
		foreach item [list $Mx $My $Mz] {
		    if {$item !="0"} {
		        set usepointforce "Yes"
		        break
		    }
		}
	    }
	}
	if {$usepointforce =="Yes"} {
	    # Active the PointLoad condition

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
	if {$Mx ne 0} {
	    set xitem [lindex $kwordlist 0]
	    GiD_File fprintf $filechannel "%s" "Begin NodalData $xitem // GUI puntual load group identifier: $cgroupid"
	    foreach node_id [GiD_EntitiesGroups get $cgroupid nodes] {
    GiD_File fprintf $filechannel "%10i %2i %20.10f" $node_id $IsFixed $Mx
	    }
	    GiD_File fprintf $filechannel "%s" "End NodalData"
	    GiD_File fprintf $filechannel ""
	}
	if {$My ne 0} {
	    set yitem [lindex $kwordlist 1]
	    GiD_File fprintf $filechannel "%s" "Begin NodalData $yitem // GUI puntual load group identifier: $cgroupid"
	    foreach node_id [GiD_EntitiesGroups get $cgroupid nodes] {
    GiD_File fprintf $filechannel "%10i %2i %20.10f" $node_id $IsFixed $My
	    }
	    GiD_File fprintf $filechannel "%s" "End NodalData"
	    GiD_File fprintf $filechannel ""
	}
	if {$Mz ne 0} {
	    set zitem [lindex $kwordlist 2]
	    GiD_File fprintf $filechannel "%s" "Begin NodalData $zitem // GUI puntual load group identifier: $cgroupid"
	    foreach node_id [GiD_EntitiesGroups get $cgroupid nodes] {
    GiD_File fprintf $filechannel "%10i %2i %20.10f" $node_id $IsFixed $Mz
	    }
	    GiD_File fprintf $filechannel "%s" "End NodalData"
	    GiD_File fprintf $filechannel "%s" ""
	}
    }

    # For debug
    if {!$::wkcf::pflag} {
	set endtime [clock seconds]
	set ttime [expr $endtime-$inittime]
	# WarnWinText "endtime:$endtime ttime:$ttime"
	WarnWinText "Write structural analysis concentrated loads (puntual force or moment): [::KUtils::Duration $ttime]"
    }
}

proc ::wkcf::WriteLineLoad {AppId cloadtid} {
    # ASBTRACT: Write beam uniformly distributed loads (forces)
    variable ndime;  variable dprops
    variable filechannel;  variable sa_icondid
    variable useqelem

    # For debug
    if {!$::wkcf::pflag} {
	set inittime [clock seconds]
    }

    # Kratos key word xpath
    set kwxpath "Applications/$AppId"
    # Set the GiD element type
    set GiDElemType "Linear"
    # Set dimension and number of nodes
    set etbf ""
    set etbf [::wkcf::GetnDimnNode $GiDElemType $ndime]
    # Set the condition keyword
    set ConditionKeyWord "LineLoadCondition"
    set LineLoadCondition [string trim ${ConditionKeyWord}${etbf}]
    # Set the reference property id
    set RefPropId "1"

    # For all defined group identifier inside this load type
    foreach cgroupid $dprops($AppId,Loads,$cloadtid,AllGroupId) {
	# wa "cgroupid:$cgroupid"
	if {[GiD_EntitiesGroups get $cgroupid elements -count -element_type $GiDElemType]} {
	    # Write LineLoadCondition2N condition
	    GiD_File fprintf $filechannel "%s" "Begin Conditions $LineLoadCondition // GUI beam uniformly distributed load group identifier: $cgroupid"
	    foreach elem_id [GiD_EntitiesGroups get $cgroupid elements -element_type $GiDElemType] {
		if {$useqelem=="1"} {
		    set nodes [lrange [GiD_Mesh get element $elem_id] 3 end]
		    set N1 [lindex $nodes 0]
		    set N2 [lindex $nodes 1]
		    set N3 [lindex $nodes 2]
		    incr sa_icondid
		    set cf "[format "%4i%4i%8i%8i%8i" $sa_icondid $RefPropId $N1 $N2 $N3]"
		  } else {
		    set nodes [lrange [GiD_Mesh get element $elem_id] 3 end]
		    set N1 [lindex $nodes 0]
		    set N2 [lindex $nodes 1]
		    incr sa_icondid
		    set cf "[format "%4i%4i%8i%8i" $sa_icondid $RefPropId $N1 $N2]"
		}
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
		if {([GiD_EntitiesGroups get $cgroupid nodes -count]>0)} {

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

    # For debug
    if {!$::wkcf::pflag} {
	set endtime [clock seconds]
	set ttime [expr $endtime-$inittime]
	# WarnWinText "endtime:$endtime ttime:$ttime"
	WarnWinText "Write structural analysis beam uniformly distributed loads (forces): [::KUtils::Duration $ttime]"
    }
}

proc ::wkcf::WriteLinePressureLoad {AppId cloadtid} {
    # ABSTRACT: Write pressure loads (positive or negative shell pressure)
    variable ndime;  variable dprops
    variable filechannel; variable sa_icondid
    variable useqelem

    # For debug
    if {!$::wkcf::pflag} {
	set inittime [clock seconds]
    }

    # Kratos key word xpath
    set kwxpath "Applications/$AppId"
    # Set the GiD element type
    set GiDElemType "Linear"
    # Set dimension and number of nodes
    set etbf ""
    set etbf [::wkcf::GetnDimnNode $GiDElemType $ndime]
    # Set the condition keyword
    set ConditionKeyWord "LineLoadCondition"
    set LinePressureCondition [string trim ${ConditionKeyWord}${etbf}]

    # Set the reference property id
    set RefPropId "1"

    # For all defined group identifier inside this load type
    foreach cgroupid $dprops($AppId,Loads,$cloadtid,AllGroupId) {
	if {[GiD_EntitiesGroups get $cgroupid elements -count -element_type $GiDElemType]} {
	    # Write LineLoadCondition2N condition
	    GiD_File fprintf $filechannel "%s" "Begin Conditions $LinePressureCondition // GUI pressure load group identifier: $cgroupid"
	    foreach elem_id [GiD_EntitiesGroups get $cgroupid elements -element_type $GiDElemType] {
		if {$useqelem=="1"} {
		    set nodes [lrange [GiD_Mesh get element $elem_id] 3 end]
		    set N1 [lindex $nodes 0]
		    set N2 [lindex $nodes 1]
		    set N3 [lindex $nodes 2]
		    incr sa_icondid
		    set cf "[format "%4i%4i%8i%8i%8i" $sa_icondid $RefPropId $N1 $N2 $N3]"
		  } else {
		    set nodes [lrange [GiD_Mesh get element $elem_id] 3 end]
		    set N1 [lindex $nodes 0]
		    set N2 [lindex $nodes 1]
		    incr sa_icondid
		    set cf "[format "%4i%4i%8i%8i" $sa_icondid $RefPropId $N1 $N2]"
		}
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
	if {([GiD_EntitiesGroups get $cgroupid nodes -count]>0)} {
	    # Get the current pressure keyword
	    set kwid "${PressureType}Pressure"
	    set ckword [::xmlutils::getKKWord $kwxpath $kwid]
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

    # For debug
    if {!$::wkcf::pflag} {
	set endtime [clock seconds]
	set ttime [expr $endtime-$inittime]
	# WarnWinText "endtime:$endtime ttime:$ttime"
	WarnWinText "Write structural analysis write pressure loads (positive or negative shell pressure): [::KUtils::Duration $ttime]"
    }
}

# Axisymmetric loads:

proc ::wkcf::WritePointLoadsAxisym {AppId cloadtid kwordlist} {
    # ABSTRACT: Write concentrated loads (puntual force or moment)
    variable ndime; variable dprops
    variable filechannel; variable sa_icondid

    # For debug
    if {!$::wkcf::pflag} {
	set inittime [clock seconds]
    }

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
    set pointforcekword "AxsymPointLoadCondition2D1N"
		foreach item [list $Fx $Fy $Mz] {
		    if {$item !="0"} {
		        set usepointforce "Yes"
		        break
		    }
		}
	    }

	}
	if {$usepointforce =="Yes"} {
	    # Active the PointLoad condition

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

    # For debug
    if {!$::wkcf::pflag} {
	set endtime [clock seconds]
	set ttime [expr $endtime-$inittime]
	# WarnWinText "endtime:$endtime ttime:$ttime"
	WarnWinText "Write structural analysis concentrated loads (puntual force or moment): [::KUtils::Duration $ttime]"
    }
}

proc ::wkcf::WriteLineLoadAxisym {AppId cloadtid} {
    # ASBTRACT: Write beam uniformly distributed loads (forces)
    variable ndime;  variable dprops
    variable filechannel;  variable sa_icondid
    variable useqelem

    # For debug
    if {!$::wkcf::pflag} {
	set inittime [clock seconds]
    }

    # Kratos key word xpath
    set kwxpath "Applications/$AppId"
    # Set the GiD element type
    set GiDElemType "Linear"
    # Set dimension and number of nodes
    set etbf ""
    set etbf [::wkcf::GetnDimnNode $GiDElemType $ndime]
    # Set the condition keyword
    set ConditionKeyWord "LineLoadAxisymCondition"
    set LineLoadCondition [string trim ${ConditionKeyWord}${etbf}]
    # Set the reference property id
    set RefPropId "1"

    # For all defined group identifier inside this load type
    foreach cgroupid $dprops($AppId,Loads,$cloadtid,AllGroupId) {
	# wa "cgroupid:$cgroupid"
	if {[GiD_EntitiesGroups get $cgroupid elements -count -element_type $GiDElemType]} {
	    # Write LineLoadCondition2N condition
	    GiD_File fprintf $filechannel "%s" "Begin Conditions $LineLoadCondition // GUI beam uniformly distributed load group identifier: $cgroupid"
	    foreach elem_id [GiD_EntitiesGroups get $cgroupid elements -element_type $GiDElemType] {
		if {$useqelem=="1"} {
		    set nodes [lrange [GiD_Mesh get element $elem_id] 3 end]
		    set N1 [lindex $nodes 0]
		    set N2 [lindex $nodes 1]
		    set N3 [lindex $nodes 2]
		    incr sa_icondid
		    set cf "[format "%4i%4i%8i%8i%8i" $sa_icondid $RefPropId $N1 $N2 $N3]"
		  } else {
		    set nodes [lrange [GiD_Mesh get element $elem_id] 3 end]
		    set N1 [lindex $nodes 0]
		    set N2 [lindex $nodes 1]
		    incr sa_icondid
		    set cf "[format "%4i%4i%8i%8i" $sa_icondid $RefPropId $N1 $N2]"
		}
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
		if {([GiD_EntitiesGroups get $cgroupid nodes -count]>0)} {

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

    # For debug
    if {!$::wkcf::pflag} {
	set endtime [clock seconds]
	set ttime [expr $endtime-$inittime]
	# WarnWinText "endtime:$endtime ttime:$ttime"
	WarnWinText "Write structural analysis beam uniformly distributed loads (forces): [::KUtils::Duration $ttime]"
    }
}

proc ::wkcf::WriteLinePressureLoadAxisym {AppId cloadtid} {
    # ABSTRACT: Write pressure loads (positive or negative shell pressure)
    variable ndime;  variable dprops
    variable filechannel; variable sa_icondid
    variable useqelem

    # For debug
    if {!$::wkcf::pflag} {
	set inittime [clock seconds]
    }

    # Kratos key word xpath
    set kwxpath "Applications/$AppId"
    # Set the GiD element type
    set GiDElemType "Linear"
    # Set dimension and number of nodes
    set etbf ""
    set etbf [::wkcf::GetnDimnNode $GiDElemType $ndime]
    # Set the condition keyword
    set ConditionKeyWord "LineLoadAxisymCondition"
    set LinePressureCondition [string trim ${ConditionKeyWord}${etbf}]

    # Set the reference property id
    set RefPropId "1"

    # For all defined group identifier inside this load type
    foreach cgroupid $dprops($AppId,Loads,$cloadtid,AllGroupId) {
	if {[GiD_EntitiesGroups get $cgroupid elements -count -element_type $GiDElemType]} {
	    # Write LineLoadCondition2N condition
	    GiD_File fprintf $filechannel "%s" "Begin Conditions $LinePressureCondition // GUI pressure load group identifier: $cgroupid"
	    foreach elem_id [GiD_EntitiesGroups get $cgroupid elements -element_type $GiDElemType] {
		if {$useqelem=="1"} {
		    set nodes [lrange [GiD_Mesh get element $elem_id] 3 end]
		    set N1 [lindex $nodes 0]
		    set N2 [lindex $nodes 1]
		    set N3 [lindex $nodes 2]
		    incr sa_icondid
		    set cf "[format "%4i%4i%8i%8i%8i" $sa_icondid $RefPropId $N1 $N2 $N3]"
		  } else {
		    set nodes [lrange [GiD_Mesh get element $elem_id] 3 end]
		    set N1 [lindex $nodes 0]
		    set N2 [lindex $nodes 1]
		    incr sa_icondid
		    set cf "[format "%4i%4i%8i%8i" $sa_icondid $RefPropId $N1 $N2]"
		}
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
	if {([GiD_EntitiesGroups get $cgroupid nodes -count]>0)} {
	    # Get the current pressure keyword
	    set kwid "${PressureType}Pressure"
	    set ckword [::xmlutils::getKKWord $kwxpath $kwid]
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


    # Write faceforce values for all nodes inside a group identifier
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


  if {([GiD_EntitiesGroups get $cgroupid nodes -count]>0)} {

      foreach Component {LineDFx LineDFy LineDFz} {
        # Get the current face load component keyword
        set kwid "$Component"
        set ckword [::xmlutils::getKKWord $kwxpath $kwid]
        if {[set $Component] !="0.0"} {
    # Set the real pressure value
    # set PressureValue [expr double($PressureValue)/$ngroupnodes]
    # Write the pressure values
    # Pressure properties
    GiD_File fprintf $filechannel "%s" "Begin NodalData $ckword // GUI pressure load group identifier: $cgroupid"
    foreach node_id [GiD_EntitiesGroups get $cgroupid nodes] {
        GiD_File fprintf $filechannel "%10i %6i %20.10f" $node_id 0 [set $Component]
    }
    GiD_File fprintf $filechannel "%s" "End NodalData"
    GiD_File fprintf $filechannel ""
        }
	    }
	}
    }

    # For debug
    if {!$::wkcf::pflag} {
	set endtime [clock seconds]
	set ttime [expr $endtime-$inittime]
	# WarnWinText "endtime:$endtime ttime:$ttime"
	WarnWinText "Write structural analysis write pressure loads (positive or negative shell pressure): [::KUtils::Duration $ttime]"
    }
}


#Surface Loads:

proc ::wkcf::WriteSurfaceLoad {AppId cloadtid} {
    # ABSTRACT: Write surface loads
    variable ndime;  variable dprops
    variable filechannel; variable sa_icondid
    variable useqelem

    # For debug
    if {!$::wkcf::pflag} {
	set inittime [clock seconds]
    }

    # Kratos key word xpath
    set kwxpath "Applications/$AppId"

    # TRIANGULAR SURFACE:

    # Set the GiD element type
    set GiDElemType "Triangle"
    # Set dimension and number of nodes
    set etbf ""
    set etbf [::wkcf::GetnDimnNode $GiDElemType $ndime]
    # Set the condition keyword
    set ConditionKeyWord "SurfaceLoadCondition"
    set SurfaceLoadCondition [string trim ${ConditionKeyWord}${etbf}]

    # Set the reference property id
    set RefPropId "1"

    # For all defined group identifier inside this load type
    foreach cgroupid $dprops($AppId,Loads,$cloadtid,AllGroupId) {
	if {[GiD_EntitiesGroups get $cgroupid elements -count -element_type $GiDElemType]} {
	    # Write SurfaceLoad condition
	    GiD_File fprintf $filechannel "%s" "Begin Conditions $SurfaceLoadCondition // GUI face load group identifier: $cgroupid"
	    foreach elem_id [GiD_EntitiesGroups get $cgroupid elements -element_type $GiDElemType] {
		if {$useqelem=="1"} {
		    set nodes [lrange [GiD_Mesh get element $elem_id] 3 end]
		    set N1 [lindex $nodes 0]
		    set N2 [lindex $nodes 1]
		    set N3 [lindex $nodes 2]
		    set N4 [lindex $nodes 3]
		    set N5 [lindex $nodes 4]
		    set N6 [lindex $nodes 5]
		    incr sa_icondid
		    set cf "[format "%4i%4i%8i%8i%8i%8i%8i%8i" $sa_icondid $RefPropId $N1 $N2 $N3 $N4 $N5 $N6]"
		} else {
		     set nodes [lrange [GiD_Mesh get element $elem_id] 3 end]
		    set N1 [lindex $nodes 0]
		    set N2 [lindex $nodes 1]
		    set N3 [lindex $nodes 2]
		    incr sa_icondid
		    set cf "[format "%4i%4i%8i%8i%8i" $sa_icondid $RefPropId $N1 $N2 $N3]"
		}
		  GiD_File fprintf $filechannel "%s" "$cf"
	    }
	    GiD_File fprintf $filechannel "%s" "End Conditions"
	    GiD_File fprintf $filechannel ""
	}
    }

    # QUADRILATERAL SURFACE:

    # Set the GiD element type
    set GiDElemType "Quadrilateral"
    # Set dimension and number of nodes
    set etbf ""
    set etbf [::wkcf::GetnDimnNode $GiDElemType $ndime]
    # Set the condition keyword
    set ConditionKeyWord "SurfaceLoadCondition"
    set SurfaceLoadCondition [string trim ${ConditionKeyWord}${etbf}]
    # wa "SurfaceLoadCondition:$SurfaceLoadCondition"

    # For all defined group identifier inside this load type
    foreach cgroupid $dprops($AppId,Loads,$cloadtid,AllGroupId) {
	# wa "cgroupid:$cgroupid"
	if {[GiD_EntitiesGroups get $cgroupid elements -count -element_type $GiDElemType]} {
	    # Write SurfaceLoad condition
	    GiD_File fprintf $filechannel "%s" "Begin Conditions $SurfaceLoadCondition // GUI face load group identifier: $cgroupid"
	    foreach elem_id [GiD_EntitiesGroups get $cgroupid elements -element_type $GiDElemType] {
		if {$useqelem=="1"} {
		    set nodes [lrange [GiD_Mesh get element $elem_id] 3 end]
		    set N1 [lindex $nodes 0]
		    set N2 [lindex $nodes 1]
		    set N3 [lindex $nodes 2]
		    set N4 [lindex $nodes 3]
		    set N5 [lindex $nodes 4]
		    set N6 [lindex $nodes 5]
		    set N7 [lindex $nodes 6]
		    set N8 [lindex $nodes 7]
		    incr sa_icondid
		    set cf "[format "%4i%4i%8i%8i%8i%8i%8i%8i%8i%8i" $sa_icondid $RefPropId $N1 $N2 $N3 $N4 $N5 $N6 $N7 $N8]"

		} else {
		    set nodes [lrange [GiD_Mesh get element $elem_id] 3 end]
		    set N1 [lindex $nodes 0]
		    set N2 [lindex $nodes 1]
		    set N3 [lindex $nodes 2]
		    set N4 [lindex $nodes 3]
		    incr sa_icondid
		    set cf "[format "%4i%4i%8i%8i%8i%8i" $sa_icondid $RefPropId $N1 $N2 $N3 $N4]"
		}
		  GiD_File fprintf $filechannel "%s" "$cf"
	    }
	    GiD_File fprintf $filechannel "%s" "End Conditions"
	    GiD_File fprintf $filechannel ""
	}
    }

    # Write faceforce values for all nodes inside a group identifier
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

	#values:
	#<Container id="Values" pid="Values" help="Set the values">
	#        <Item id="Vx" pid="X" dv="0.0" help="X coordinate"/>
	#        <Item id="Vy" pid="Y" dv="0.0" help="Y coordinate"/>
	#        <Item id="Vz" pid="Z" dv="0.0" nDim="3D" help="Z coordinate"/>
	#</Container>
	#<Container id="Activation" pid="Fixed" help="Fix/release some degree of freedom">
	#        <Item id="Ax" pid="X active" dv="1" ivalues="1,0" values="1,0" help="Fix X degree of freedom"/>
	#        <Item id="Ay" pid="Y active" dv="1" ivalues="1,0" values="1,0" help="Fix Y degree of freedom"/>
	#        <Item id="Az" pid="Z active" dv="0" nDim="3D" ivalues="1,0" values="1,0" help="Fix Z degree of freedom"/>
	#</Container>

	# WarnWinText "FixPressure:$FixPressure PressureType:$PressureType PressureValue:$PressureValue"
	if {([GiD_EntitiesGroups get $cgroupid nodes -count]>0)} {

      foreach Component {SurfaceDFx SurfaceDFy SurfaceDFz} {
	      # Get the current face load component keyword
	      set kwid "$Component"
	      set ckword [::xmlutils::getKKWord $kwxpath $kwid]
	      if {[set $Component] !="0.0"} {
		# Set the real pressure value
		# set PressureValue [expr double($PressureValue)/$ngroupnodes]
		# Write the pressure values
		# Pressure properties
		GiD_File fprintf $filechannel "%s" "Begin NodalData $ckword // GUI pressure load group identifier: $cgroupid"
		foreach node_id [GiD_EntitiesGroups get $cgroupid nodes] {
		    GiD_File fprintf $filechannel "%10i %6i %20.10f" $node_id 0 [set $Component]
		}
		GiD_File fprintf $filechannel "%s" "End NodalData"
		GiD_File fprintf $filechannel ""
	      }
	    }
	}
    }

    # For debug
    if {!$::wkcf::pflag} {
	set endtime [clock seconds]
	set ttime [expr $endtime-$inittime]
	# WarnWinText "endtime:$endtime ttime:$ttime"
	WarnWinText "Write structural analysis write face loads: [::KUtils::Duration $ttime]"
    }
}



proc ::wkcf::WriteSurfacePressure {AppId cloadtid} {
    # ABSTRACT: Write face loads (positive or negative shell pressure)
    variable ndime;  variable dprops
    variable filechannel; variable sa_icondid
    variable useqelem

    # For debug
    if {!$::wkcf::pflag} {
	set inittime [clock seconds]
    }

    # Kratos key word xpath
    set kwxpath "Applications/$AppId"

    # TRIANGULAR SURFACE:

    # Set the GiD element type
    set GiDElemType "Triangle"
    # Set dimension and number of nodes
    set etbf ""
    set etbf [::wkcf::GetnDimnNode $GiDElemType $ndime]
    # Set the condition keyword
    set ConditionKeyWord "SurfaceLoadCondition"
    set SurfacePressureCondition [string trim ${ConditionKeyWord}${etbf}]

    # Set the reference property id
    set RefPropId "1"

    # For all defined group identifier inside this load type
    foreach cgroupid $dprops($AppId,Loads,$cloadtid,AllGroupId) {
	if {[GiD_EntitiesGroups get $cgroupid elements -count -element_type $GiDElemType]} {
	    # Write SurfaceLoad condition
	    GiD_File fprintf $filechannel "%s" "Begin Conditions $SurfacePressureCondition // GUI face load group identifier: $cgroupid"
	    foreach elem_id [GiD_EntitiesGroups get $cgroupid elements -element_type $GiDElemType] {

		if {$useqelem=="1"} {
        set nodes [lrange [GiD_Mesh get element $elem_id] 3 end]
		    set N1 [lindex $nodes 0]
		    set N2 [lindex $nodes 1]
		    set N3 [lindex $nodes 2]
		    set N4 [lindex $nodes 3]
		    set N5 [lindex $nodes 4]
		    set N6 [lindex $nodes 5]
		    incr sa_icondid
		    set cf "[format "%4i%4i%8i%8i%8i%8i%8i%8i" $sa_icondid $RefPropId $N1 $N2 $N3 $N4 $N5 $N6]"

		} else {
		     set nodes [lrange [GiD_Mesh get element $elem_id] 3 end]
		    set N1 [lindex $nodes 0]
		    set N2 [lindex $nodes 1]
		    set N3 [lindex $nodes 2]
		    incr sa_icondid
		    set cf "[format "%4i%4i%8i%8i%8i" $sa_icondid $RefPropId $N1 $N2 $N3]"
		}
		  GiD_File fprintf $filechannel "%s" "$cf"
	    }
	    GiD_File fprintf $filechannel "%s" "End Conditions"
	    GiD_File fprintf $filechannel ""
	}
    }


    # QUADRILATERAL SURFACE:

    # Set the GiD element type
    set GiDElemType "Quadrilateral"
    # Set dimension and number of nodes
    set etbf ""
    set etbf [::wkcf::GetnDimnNode $GiDElemType $ndime]
    # Set the condition keyword
    set ConditionKeyWord "SurfaceLoadCondition"
    set SurfacePressureCondition [string trim ${ConditionKeyWord}${etbf}]

    # For all defined group identifier inside this load type
    foreach cgroupid $dprops($AppId,Loads,$cloadtid,AllGroupId) {
	if {[GiD_EntitiesGroups get $cgroupid elements -count -element_type $GiDElemType]} {
	    # Write SurfaceLoad condition
	    GiD_File fprintf $filechannel "%s" "Begin Conditions $SurfacePressureCondition // GUI face load group identifier: $cgroupid"
	    foreach elem_id [GiD_EntitiesGroups get $cgroupid elements -element_type $GiDElemType] {

		if {$useqelem=="1"} {
        set nodes [lrange [GiD_Mesh get element $elem_id] 3 end]
		    set N1 [lindex $nodes 0]
		    set N2 [lindex $nodes 1]
		    set N3 [lindex $nodes 2]
		    set N4 [lindex $nodes 3]
		    set N5 [lindex $nodes 4]
		    set N6 [lindex $nodes 5]
		    set N7 [lindex $nodes 6]
		    set N8 [lindex $nodes 7]
		    incr sa_icondid
		    set cf "[format "%4i%4i%8i%8i%8i%8i%8i%8i%8i%8i" $sa_icondid $RefPropId $N1 $N2 $N3 $N4 $N5 $N6 $N7 $N8]"
		} else {
        set nodes [lrange [GiD_Mesh get element $elem_id] 3 end]
		    set N1 [lindex $nodes 0]
		    set N2 [lindex $nodes 1]
		    set N3 [lindex $nodes 2]
		    set N4 [lindex $nodes 3]
		    incr sa_icondid
		    set cf "[format "%4i%4i%8i%8i%8i%8i" $sa_icondid $RefPropId $N1 $N2 $N3 $N4]"
		}
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
	    set ckword [::xmlutils::getKKWord $kwxpath $kwid]
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

    # For debug
    if {!$::wkcf::pflag} {
	set endtime [clock seconds]
	set ttime [expr $endtime-$inittime]
	# WarnWinText "endtime:$endtime ttime:$ttime"
	WarnWinText "Write structural analysis write pressure loads (positive or negative shell pressure): [::KUtils::Duration $ttime]"
    }
}


proc ::wkcf::WriteBodyForceValues {props} {
    # Write the gravity properties to the kratos data file
    # Arguments
    # props => Body force properties
    variable ndime; variable filechannel

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
    GiD_File fprintf $filechannel "%s" " BODY_FORCE \[3\] $vector"
}

proc ::wkcf::WriteDSOLIDMaterials {AppId fileid PDir} {
    variable DSOLID
    global solid_element_props
    set rootid DSOLID
  	set basexpath "$rootid//c.Solid-Elements//c.Solid-Element"
  	set gproplist [::xmlutils::setXmlContainerIds $basexpath]
    puts $fileid "# Importing the Kratos Library"
    puts $fileid "from KratosMultiphysics import *"
    puts $fileid "from KratosMultiphysics.SolidMechanicsApplication import *"
    puts $fileid "def AssignMaterial(Properties):"
    foreach cgroupid $gproplist {
    puts $fileid "# GUI property identifier: Property $solid_element_props($cgroupid)"
    puts $fileid "# GUI material identifier: Steel_AISI1059"
    puts $fileid "    prop_id = $solid_element_props($cgroupid);"
    puts $fileid "    prop = Properties\[prop_id\]"
    set cproperty "dv"
    set cxpath "$rootid//c.Solid-Elements//c.Solid-Element//c.[list ${cgroupid}]//c.Properties//i.MatModel"
    set matmodel [::xmlutils::setXml $cxpath $cproperty]
    if {$matmodel eq "Elastic-Isotropic"} {    puts $fileid "    mat = LinearElastic3DLaw();"
    } elseif {$matmodel eq "HyperElastic-Isotropic"} {    puts $fileid "    mat = HyperElastic3DLaw();"
    } elseif {$matmodel eq "HyperElastic-Plastic"} {    puts $fileid "    mat = HyperElasticPlasticJ23DLaw();"
    }
    puts $fileid "    prop.SetValue(CONSTITUTIVE_LAW, mat.Clone());"
    }
}

proc ::wkcf::WriteSolidProjectParameters {AppId fileid PDir} {
    # Write the project parameters file
    variable ndime; variable dprops
    variable DSOLID

    # kratos key word xpath
    set kwxpath "Applications/$AppId"
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
    puts $fileid "#Problem Data"
    puts $fileid "#####WriteDSolidProjectParameters##########"
    puts $fileid ""
    puts $fileid "ProblemType = \"Mechanical\""

    # Number of threads
    set cxpath "$AppId//c.SolutionStrategy//c.ParallelType//i.ParallelSolutionType"
    set ParallelType [::xmlutils::setXml $cxpath $cproperty]
    if {$ParallelType =="OpenMP"} {
	# Number of Steps
	set cxpath "$AppId//c.SolutionStrategy//c.ParallelType//i.OpenMPNumberOfThreads"
	set NumberOfThreads [::xmlutils::setXml $cxpath $cproperty]
	puts $fileid "NumberofThreads = $NumberOfThreads"
    } else {
	puts $fileid "NumberofThreads = 1"
    }

    # Check for use shell elements
    set usenbst "No"
    set useshells "No"
    set usebeams "No"

    set shelllist [list "EBST" "ShellThick" "ShellThickCR" "ShellThin" "ShellThinCR"]
    set beamlist  [list "BeamElement"]

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


    # Solution method
    set cxpath "$AppId//c.AnalysisData//i.AnalysisType"
    set AnalysisType [::xmlutils::setXml $cxpath $cproperty]
    if {$AnalysisType =="Non-Linear"} {

	set cxpath "$AppId//c.SolutionStrategy//c.Non-Linear//i.SolutionMethod"
	set SolutionMethod [::xmlutils::setXml $cxpath $cproperty]
	puts $fileid "Solution_method = \"$SolutionMethod\""
      } else {
  	puts $fileid "Solution_method = \"Newton-Raphson\""
      }

    # Solution type
    set cxpath "$AppId//c.AnalysisData//i.SolutionType"
    set SolutionType [::xmlutils::setXml $cxpath $cproperty]
    # WarnWinText "SolverType:$SolutionType"
    
    puts $fileid "SolverType = \"DynamicSolver\""

  	# Delta time
  	set cxpath "$AppId//c.SolutionStrategy//c.Dynamic//i.DeltaTime"
  	set DeltaTime [::xmlutils::setXml $cxpath $cproperty]
  	puts $fileid "time_step = $DeltaTime"

  	# Number of Steps
  	set cxpath "$AppId//c.SolutionStrategy//c.Dynamic//i.EndTime"
  	set EndTime [::xmlutils::setXml $cxpath $cproperty]
    puts $fileid "end_time = $EndTime"

    puts $fileid ""
    puts $fileid "#Solver Data"
    puts $fileid "#################################################"

    set trailing_spaces  "    "
    # Structural solver configuration **********************
    puts $fileid ""
    puts $fileid "class SolverSettings:"
    puts $fileid "${trailing_spaces}solver_type = \"mechanical_solver\""
    puts $fileid "${trailing_spaces}domain_size = $domain_size"
    puts $fileid "${trailing_spaces}echo_level  = 0"

    # Solution Methods
    set cxpath "$AppId//c.AnalysisData//i.TimeIntegrationMethod"
    set cproperty "dv"
    set TimeIntegrationMethod [::xmlutils::setXml $cxpath $cproperty]
    set cxpath "$AppId//c.SolutionStrategy//c.Dynamic//i.TimeStepPredictionLevel"
    set TimeStepPredictionLevel [::xmlutils::setXml $cxpath $cproperty]
    set cxpath "$AppId//c.SolutionStrategy//c.Dynamic//i.UseRayleighDamping"
    set UseRayleighDamping [::xmlutils::setXml $cxpath $cproperty]

    puts $fileid ""
    puts $fileid "${trailing_spaces}max_delta_time  = time_step"
    puts $fileid "${trailing_spaces}time_integration_method = \"Explicit\""
    puts $fileid "${trailing_spaces}explicit_integration_scheme = \"CentralDifferences\""
    puts $fileid "${trailing_spaces}time_step_prediction_level  = \"$TimeStepPredictionLevel\""

    if {($UseRayleighDamping eq "Yes")} {
    puts $fileid "${trailing_spaces}rayleigh_damping = True"
    } else {
    puts $fileid "${trailing_spaces}rayleigh_damping = False"
    }
    puts $fileid ""

    if {($useshells eq "Yes")||($usebeams eq "Yes")} {
	puts $fileid "${trailing_spaces}RotationDofs = True"
    } else {
	puts $fileid "${trailing_spaces}RotationDofs = False"
    }
    puts $fileid "${trailing_spaces}PressureDofs = False"

    puts $fileid "${trailing_spaces}ReformDofSetAtEachStep = False"
    puts $fileid "${trailing_spaces}LineSearch = False"
    puts $fileid "${trailing_spaces}Implex = False"
    puts $fileid "${trailing_spaces}ComputeReactions = True"
    puts $fileid "${trailing_spaces}ComputeContactForces = False"

    # Solution type
    set cxpath "$AppId//c.AnalysisData//i.SolutionType"
    set cproperty "dv"
    set SolutionType [::xmlutils::setXml $cxpath $cproperty]
    # WarnWinText "SolverType:$SolutionType"

    puts $fileid "${trailing_spaces}scheme_type = \"DynamicSolver\""
   
    # Analysis type
    set cxpath "$AppId//c.AnalysisData//i.AnalysisType"
    set AnalysisType [::xmlutils::setXml $cxpath $cproperty]
    # WarnWinText "AnalysisType:$AnalysisType"
    if {$AnalysisType =="Non-Linear"} {
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
		puts $fileid "${trailing_spaces}convergence_criterion = \"Displacement_criteria\""
	    }
	    "Residual" {
		puts $fileid "${trailing_spaces}convergence_criterion = \"Residual_criteria\""
	    }
	    "DisplacementAndResidual" {
		puts $fileid "${trailing_spaces}convergence_criterion = \"And_criteria\""
	    }
	    "DisplacementOrResidual" {
		puts $fileid "${trailing_spaces}convergence_criterion = \"Or_criteria\""
	    }
	}
	puts $fileid "${trailing_spaces}displacement_relative_tolerance = $DisplacementConvergenceTolerance"
	puts $fileid "${trailing_spaces}displacement_absolute_tolerance = $DisplacementAbsoluteTolerance"
	puts $fileid "${trailing_spaces}residual_relative_tolerance = $ResidualConvergenceTolerance"
	puts $fileid "${trailing_spaces}residual_absolute_tolerance = $ResidualAbsoluteTolerance"
	puts $fileid "${trailing_spaces}max_iteration = $MaximumIterations"

    } elseif {$AnalysisType =="Linear"} {

	puts $fileid "${trailing_spaces}convergence_criterion = \"Residual_criteria\""
	puts $fileid "${trailing_spaces}displacement_relative_tolerance = 1.0E-4"
	puts $fileid "${trailing_spaces}displacement_absolute_tolerance = 1.0E-9"
	puts $fileid "${trailing_spaces}residual_relative_tolerance = 1.0E-4"
	puts $fileid "${trailing_spaces}residual_absolute_tolerance = 1.0E-9"
	puts $fileid "${trailing_spaces}max_iteration = 10"
    }

    set rootid "$AppId"
    ::wkcf::WriteLinearSolver $rootid $fileid $trailing_spaces "linear_solver_config"


    puts $fileid ""
    puts $fileid "#Constraints Data"
    puts $fileid "#################################################"
    puts $fileid ""
    # Incremental Load
    if {($SolutionType =="Dynamic")||($SolutionType =="Quasi-Static")||($SolutionType =="Pseudo-Dynamic")} {

	set cxpath "$AppId//c.Loads//i.IncrementalLoad"
	set incremental_load [::xmlutils::setXml $cxpath $cproperty]
	if {($incremental_load eq "Yes")} {
	    puts $fileid "Incremental_Load = \"True\""
	} else {
	    puts $fileid "Incremental_Load = \"False\""
	}
	# Incremental Displacement
	set cxpath "$AppId//c.Conditions//i.IncrementalMovement"
	set incremental_movement [::xmlutils::setXml $cxpath $cproperty]
	if {($incremental_movement eq "Yes")} {
	    puts $fileid "Incremental_Displacement = \"True\""
	} else {
	    puts $fileid "Incremental_Displacement = \"False\""
	}

    } else {

	puts $fileid "Incremental_Load = \"False\""
	puts $fileid "Incremental_Displacement = \"False\""

    }

    puts $fileid ""
    puts $fileid "#PostProcess Data"
    puts $fileid "#################################################"
    puts $fileid ""
    # For results

    # On nodes results
    set cnrlist [list "Displacements" "Velocities" ]
     #"Accelerations" "Rotations" "Reactions" "Forces"
    set nodal_results "nodal_results=\["
    foreach cnr $cnrlist {
    	set cxpath "DSOLID//c.Solid-Post//c.OnNodes//i.[list ${cnr}]"
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
    }
	  append nodal_results "\]"
	  puts $fileid "$nodal_results"
    # WarnWinText "nodal_results:$nodal_results"

    # On Gauss point results
    set cgrlist [list "StrainTensor" "StressTensor" "VonMises" ]
    #"PlasticStrain" "DeltaPlasticStrain" "BeamMoments" "BeamForces" "ShellForcesLocal" "ShellForcesGlobal" "ShellMomentsLocal" "ShellMomentsGlobal" "ShellStrainLocal" "ShellStrainGlobal" "ShellCurvatureLocal" "ShellCurvatureGlobal" "MaterialDirectionX" "MaterialDirectionY" "MaterialDirectionZ"]
    set gauss_points_results "gauss_points_results=\["
    foreach cgr $cgrlist {
  	set cxpath "DSOLID//c.Solid-Post//c.OnGaussPoints//i.[list ${cgr}]"
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
    }
	append gauss_points_results "\]"
	puts $fileid "$gauss_points_results"
    # WarnWinText "gauss_points_results:$gauss_points_results"

    # GiD post mode variables
    #::wkcf::WriteGiDPostMode $AppId $fileid
    ::wkcf::WriteGiDPostModeNew $AppId $fileid
    puts $fileid ""
    set cxpath "$AppId//c.SolutionStrategy//c.Dynamic//i.DeltaTime"
    set TimeStep [::xmlutils::setXml $cxpath $cproperty]

    # Number of Steps
    set cxpath "$AppId//c.Results//i.OutputDeltaTime"
    set TimeWriting [::xmlutils::setXml $cxpath $cproperty]
    set FrequencyWriting [expr int(double($TimeWriting)/double($TimeStep))]

    if {$FrequencyWriting != 0 && $FrequencyWriting > 1} {
  puts $fileid "GiDWriteFrequency = $TimeWriting"
    } else {
  puts $fileid "GiDWriteFrequency = $TimeStep"
    }
    puts $fileid "WriteResults = \"PreMeshing\""
    puts $fileid "echo_level = 0"
    puts $fileid ""
    puts $fileid "# graph_options"
    puts $fileid "PlotGraphs = \"False\""
    puts $fileid "PlotFrequency = 0 "
    puts $fileid ""
    puts $fileid "# list options"
    puts $fileid "PrintLists = \"True\""
    puts $fileid "file_list = \[\] "
    puts $fileid ""
    puts $fileid "# restart options"
    puts $fileid "SaveRestart = False"
    puts $fileid "RestartFrequency = 0"
    puts $fileid "LoadRestart = False"
    puts $fileid "Restart_Step = 0"


    puts $fileid ""
    set PName [::KUtils::GetPaths "PName"]
    puts $fileid "problem_name=\"${PName}${AppId}\""
    puts $fileid "problem_path=\"[file join $PDir]\""

}









proc ::wkcf::WriteStructuralProjectParameters {AppId fileid PDir} {
    # Write the project parameters file
    variable ndime; variable dprops

    # kratos key word xpath
    set kwxpath "Applications/$AppId"

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
    puts $fileid "#Problem Data"
    puts $fileid "#################################################"
    puts $fileid ""
    puts $fileid "ProblemType = \"Mechanical\""

    # Number of threads
    set cxpath "$AppId//c.SolutionStrategy//c.ParallelType//i.ParallelSolutionType"
    set ParallelType [::xmlutils::setXml $cxpath $cproperty]
    if {$ParallelType =="OpenMP"} {
    	# Number of Steps
    	set cxpath "$AppId//c.SolutionStrategy//c.ParallelType//i.OpenMPNumberOfThreads"
    	set NumberOfThreads [::xmlutils::setXml $cxpath $cproperty]
    	puts $fileid "NumberofThreads = $NumberOfThreads"
        } else {
    	puts $fileid "NumberofThreads = 1"
    }

    # Check for use shell elements
    set usenbst "No"
    set useshells "No"
    set usebeams "No"

    set shelllist [list "EBST" "ShellThick" "ShellThickCR" "ShellThin" "ShellThinCR"]
    set beamlist  [list "BeamElement"]

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

    # Solution method
    set cxpath "$AppId//c.AnalysisData//i.AnalysisType"
    set AnalysisType [::xmlutils::setXml $cxpath $cproperty]
    if {$AnalysisType =="Non-Linear"} {
  	set cxpath "$AppId//c.SolutionStrategy//c.Non-Linear//i.SolutionMethod"
  	set SolutionMethod [::xmlutils::setXml $cxpath $cproperty]
  	puts $fileid "Solution_method = \"$SolutionMethod\""
      } else {
  	puts $fileid "Solution_method = \"Newton-Raphson\""
      }

    # Solution type
    set cxpath "$AppId//c.AnalysisData//i.SolutionType"
    set SolutionType [::xmlutils::setXml $cxpath $cproperty]
    # WarnWinText "SolverType:$SolutionType"
    if {$SolutionType =="Dynamic"} {
	puts $fileid "SolverType = \"DynamicSolver\""

	# Delta time
	set cxpath "$AppId//c.SolutionStrategy//c.Dynamic//i.DeltaTime"
	set DeltaTime [::xmlutils::setXml $cxpath $cproperty]
	puts $fileid "time_step = $DeltaTime"

	# Number of Steps
	set cxpath "$AppId//c.SolutionStrategy//c.Dynamic//i.EndTime"
	set EndTime [::xmlutils::setXml $cxpath $cproperty]
    puts $fileid "end_time = $EndTime"

    } elseif {$SolutionType =="Quasi-Static"} {

	puts $fileid "SolverType = \"QuasiStaticSolver\""

	# Delta time
	set cxpath "$AppId//c.SolutionStrategy//c.Dynamic//i.DeltaTime"
	set DeltaTime [::xmlutils::setXml $cxpath $cproperty]
	puts $fileid "time_step = $DeltaTime"

	# Number of Steps
	set cxpath "$AppId//c.SolutionStrategy//c.Dynamic//i.EndTime"
	set EndTime [::xmlutils::setXml $cxpath $cproperty]
    puts $fileid "end_time    = $EndTime"

    } elseif {$SolutionType =="Pseudo-Dynamic"} {

	puts $fileid "SolverType = \"PseudoDynamicSolver\""

	# Delta time
	set cxpath "$AppId//c.SolutionStrategy//c.Dynamic//i.DeltaTime"
	set DeltaTime [::xmlutils::setXml $cxpath $cproperty]
	puts $fileid "time_step = $DeltaTime"

	# Number of Steps
	set cxpath "$AppId//c.SolutionStrategy//c.Dynamic//i.EndTime"
	set EndTime [::xmlutils::setXml $cxpath $cproperty]
    puts $fileid "end_time    = $EndTime"

    } elseif {$SolutionType =="Static"} {
	puts $fileid "SolverType = \"StaticSolver\""
	puts $fileid "time_step = 1"
    puts $fileid "end_time    = 1"

    }

    puts $fileid ""
    puts $fileid "#Solver Data"
    puts $fileid "#################################################"

    set trailing_spaces  "    "
    # Structural solver configuration **********************
    puts $fileid ""
    puts $fileid "class SolverSettings:"
    puts $fileid "${trailing_spaces}solver_type = \"mechanical_solver\""
    puts $fileid "${trailing_spaces}domain_size = $domain_size"
    puts $fileid "${trailing_spaces}echo_level  = 0"

    # Solution Methods
    set cxpath "$AppId//c.AnalysisData//i.TimeIntegrationMethod"
    set cproperty "dv"
    set TimeIntegrationMethod [::xmlutils::setXml $cxpath $cproperty]
    set cxpath "$AppId//c.SolutionStrategy//c.Dynamic//i.TimeStepPredictionLevel"
    set TimeStepPredictionLevel [::xmlutils::setXml $cxpath $cproperty]
    set cxpath "$AppId//c.SolutionStrategy//c.Dynamic//i.UseRayleighDamping"
    set UseRayleighDamping [::xmlutils::setXml $cxpath $cproperty]

    puts $fileid ""
    puts $fileid "${trailing_spaces}max_delta_time  = time_step"
    puts $fileid "${trailing_spaces}time_integration_method = \"$TimeIntegrationMethod\""
    puts $fileid "${trailing_spaces}explicit_integration_scheme = \"CentralDifferences\""
    puts $fileid "${trailing_spaces}time_step_prediction_level  = \"$TimeStepPredictionLevel\""

    if {($UseRayleighDamping eq "Yes")} {
    puts $fileid "${trailing_spaces}rayleigh_damping = True"
    } else {
    puts $fileid "${trailing_spaces}rayleigh_damping = False"
    }
    puts $fileid ""

    if {($useshells eq "Yes")||($usebeams eq "Yes")} {
	puts $fileid "${trailing_spaces}RotationDofs = True"
    } else {
	puts $fileid "${trailing_spaces}RotationDofs = False"
    }
    puts $fileid "${trailing_spaces}PressureDofs = False"

    puts $fileid "${trailing_spaces}ReformDofSetAtEachStep = False"
    puts $fileid "${trailing_spaces}LineSearch = False"
    puts $fileid "${trailing_spaces}Implex = False"
    puts $fileid "${trailing_spaces}ComputeReactions = True"
    puts $fileid "${trailing_spaces}ComputeContactForces = False"

    # Solution type
    set cxpath "$AppId//c.AnalysisData//i.SolutionType"
    set cproperty "dv"
    set SolutionType [::xmlutils::setXml $cxpath $cproperty]
    # WarnWinText "SolverType:$SolutionType"
    if {$SolutionType =="Dynamic"} {
	puts $fileid "${trailing_spaces}scheme_type = \"DynamicSolver\""
    } elseif {$SolutionType =="Quasi-Static"} {
	puts $fileid "${trailing_spaces}scheme_type = \"QuasiStaticSolver\""
    } elseif {$SolutionType =="Pseudo-Dynamic"} {
	puts $fileid "${trailing_spaces}scheme_type = \"PseudoDynamicSolver\""
    } elseif {$SolutionType =="Static"} {
	puts $fileid "${trailing_spaces}scheme_type = \"StaticSolver\""
    }
    # Analysis type
    set cxpath "$AppId//c.AnalysisData//i.AnalysisType"
    set AnalysisType [::xmlutils::setXml $cxpath $cproperty]
    # WarnWinText "AnalysisType:$AnalysisType"
    if {$AnalysisType =="Non-Linear"} {
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
		puts $fileid "${trailing_spaces}convergence_criterion = \"Displacement_criteria\""
	    }
	    "Residual" {
		puts $fileid "${trailing_spaces}convergence_criterion = \"Residual_criteria\""
	    }
	    "DisplacementAndResidual" {
		puts $fileid "${trailing_spaces}convergence_criterion = \"And_criteria\""
	    }
	    "DisplacementOrResidual" {
		puts $fileid "${trailing_spaces}convergence_criterion = \"Or_criteria\""
	    }
	}
	puts $fileid "${trailing_spaces}displacement_relative_tolerance = $DisplacementConvergenceTolerance"
	puts $fileid "${trailing_spaces}displacement_absolute_tolerance = $DisplacementAbsoluteTolerance"
	puts $fileid "${trailing_spaces}residual_relative_tolerance = $ResidualConvergenceTolerance"
	puts $fileid "${trailing_spaces}residual_absolute_tolerance = $ResidualAbsoluteTolerance"
	puts $fileid "${trailing_spaces}max_iteration = $MaximumIterations"

	# WarnWinText "SolutionMethod:$SolutionMethod ConvergenceCriteria:$ConvergenceCriteria ConvergenceTolerance:$ConvergenceTolerance AbsoluteTolerance:$AbsoluteTolerance MaximumIterations:$MaximumIterations"
	#END NON LINEAR
    } elseif {$AnalysisType =="Linear"} {

	puts $fileid "${trailing_spaces}convergence_criterion = \"Residual_criteria\""
	puts $fileid "${trailing_spaces}displacement_relative_tolerance = 1.0E-4"
	puts $fileid "${trailing_spaces}displacement_absolute_tolerance = 1.0E-9"
	puts $fileid "${trailing_spaces}residual_relative_tolerance = 1.0E-4"
	puts $fileid "${trailing_spaces}residual_absolute_tolerance = 1.0E-9"
	puts $fileid "${trailing_spaces}max_iteration = 10"
    }

    set rootid "$AppId"
    ::wkcf::WriteLinearSolver $rootid $fileid $trailing_spaces "linear_solver_config"

    # Structural solver configuration **********************


    # puts $fileid "LineSearch = \"False\""
    # #puts $fileid "FindNodalNeighbours = \"False\""

    # if {$usenbst =="Yes"} {
    # 	#puts $fileid "FindElementalNeighbours = \"True\""
    # } else {
    # 	#puts $fileid "FindElementalNeighbours = \"False\""
    # }
    # if {($useshells eq "Yes")||($usebeams eq "Yes")} {
    # 	puts $fileid "Rotational_Dofs = \"True\""
    # } else {
    # 	puts $fileid "Rotational_Dofs = \"False\""
    # }



    # puts $fileid ""
    # puts $fileid "#Solver Data"
    # puts $fileid "#################################################"

    # # Linear Solver type

    # set cxpath "$AppId//c.SolutionStrategy//c.LinearSolver//i.LinearSolverType"
    # set LinearSolverType [::xmlutils::setXml $cxpath $cproperty]
    # if {$LinearSolverType =="Direct"} {

    # 	# Direct solver type
    # 	set cxpath "$AppId//c.SolutionStrategy//c.LinearSolver//i.DirectSolverType"
    # 	set solvertype [::xmlutils::setXml $cxpath $cproperty]
    # 	set DirectSolverType [::xmlutils::getKKWord $kwxpath $solvertype]
    # 	puts $fileid "LinearSolver = \"$DirectSolverType\""

    # 	puts $fileid "Linear_Solver_Tolerance = 1.0E-6"
    # 	puts $fileid "Linear_Solver_Max_Iteration = 5000"


    # } elseif {$LinearSolverType =="Iterative"} {

    # 	# Iterative solver type
    # 	set cxpath "$AppId//c.SolutionStrategy//c.LinearSolver//i.IterativeSolverType"
    # 	set solvertype [::xmlutils::setXml $cxpath $cproperty]
    # 	set IterativeSolverType [::xmlutils::getKKWord $kwxpath $solvertype]
    # 	puts $fileid "LinearSolver = \"$IterativeSolverType\""

    # 	# Tolerance
    # 	set cxpath "$AppId//c.SolutionStrategy//c.LinearSolver//i.Tolerance"
    # 	set Tolerance [::xmlutils::setXml $cxpath $cproperty]
    # 	puts $fileid "Linear_Solver_Tolerance = $Tolerance"

    # 	# Maximum iteration
    # 	set cxpath "$AppId//c.SolutionStrategy//c.LinearSolver//i.MaximumIteration"
    # 	set MaximumIteration [::xmlutils::setXml $cxpath $cproperty]

    # 	puts $fileid "Linear_Solver_Max_Iteration = $MaximumIteration"

    # }



    # # Analysis type
    # set cxpath "$AppId//c.AnalysisData//i.AnalysisType"
    # set cproperty "dv"
    # set AnalysisType [::xmlutils::setXml $cxpath $cproperty]
    # # WarnWinText "AnalysisType:$AnalysisType"
    # if {$AnalysisType =="Non-Linear"} {
    # 	# Convergence criteria
    # 	set cxpath "$AppId//c.SolutionStrategy//c.Non-Linear//i.ConvergenceCriteria"
    # 	set ConvergenceCriteria [::xmlutils::setXml $cxpath $cproperty]
    # 	# Residual convergence tolerance
    # 	set cxpath "$AppId//c.SolutionStrategy//c.Non-Linear//i.ResidualConvergenceTolerance"
    # 	set ResidualConvergenceTolerance [::xmlutils::setXml $cxpath $cproperty]
    # 	# Residual absolute tolerance
    # 	set cxpath "$AppId//c.SolutionStrategy//c.Non-Linear//i.ResidualAbsoluteTolerance"
    # 	set ResidualAbsoluteTolerance [::xmlutils::setXml $cxpath $cproperty]
    # 	# Displacement convergence tolerance
    # 	set cxpath "$AppId//c.SolutionStrategy//c.Non-Linear//i.DisplacementConvergenceTolerance"
    # 	set DisplacementConvergenceTolerance [::xmlutils::setXml $cxpath $cproperty]
    # 	# Displacement absolute tolerance
    # 	set cxpath "$AppId//c.SolutionStrategy//c.Non-Linear//i.DisplacementAbsoluteTolerance"
    # 	set DisplacementAbsoluteTolerance [::xmlutils::setXml $cxpath $cproperty]
    # 	# Maximum iterations
    # 	set cxpath "$AppId//c.SolutionStrategy//c.Non-Linear//i.MaximumIterations"
    # 	set MaximumIterations [::xmlutils::setXml $cxpath $cproperty]

    # 	switch -exact -- $ConvergenceCriteria {
    # 	    "Displacement" {
    # 		puts $fileid "Convergence_Criteria = \"Displacement_criteria\""
    # 		puts $fileid "Convergence_Tolerance = $DisplacementConvergenceTolerance"
    # 		puts $fileid "Absolute_Tolerance = $DisplacementAbsoluteTolerance"
    # 		# puts $fileid "Displacement_Convergence_Tolerance = $DisplacementConvergenceTolerance"
    # 		# puts $fileid "Displacement_Absolute_Tolerance = $DisplacementAbsoluteTolerance"
    # 		# puts $fileid "Residual_Convergence_Tolerance = $ResidualConvergenceTolerance"
    # 		# puts $fileid "Residual_Absolute_Tolerance = $ResidualAbsoluteTolerance"
    # 	    }
    # 	    "Residual" {
    # 		puts $fileid "Convergence_Criteria = \"Residual_criteria\""
    # 		puts $fileid "Convergence_Tolerance = $ResidualConvergenceTolerance"
    # 		puts $fileid "Absolute_Tolerance = $ResidualAbsoluteTolerance"
    # 		# puts $fileid "Residual_Convergence_Tolerance = $ResidualConvergenceTolerance"
    # 		# puts $fileid "Residual_Absolute_Tolerance = $ResidualAbsoluteTolerance"
    # 		# puts $fileid "Displacement_Convergence_Tolerance = $DisplacementConvergenceTolerance"
    # 		# puts $fileid "Displacement_Absolute_Tolerance = $DisplacementAbsoluteTolerance"
    # 	    }
    # 	    "DisplacementAndResidual" {
    # 		puts $fileid "Convergence_Criteria = \"And_criteria\""
    # 		puts $fileid "Convergence_Tolerance = $ResidualConvergenceTolerance"
    # 		puts $fileid "Absolute_Tolerance = $ResidualAbsoluteTolerance"
    # 		# puts $fileid "Residual_Convergence_Tolerance = $ResidualConvergenceTolerance"
    # 		# puts $fileid "Residual_Absolute_Tolerance = $ResidualAbsoluteTolerance"
    # 		# puts $fileid "Displacement_Convergence_Tolerance = $DisplacementConvergenceTolerance"
    # 		# puts $fileid "Displacement_Absolute_Tolerance = $DisplacementAbsoluteTolerance"
    # 	    }
    # 	    "DisplacementOrResidual" {
    # 		puts $fileid "Convergence_Criteria = \"Or_criteria\""
    # 		puts $fileid "Convergence_Tolerance = $ResidualConvergenceTolerance"
    # 		puts $fileid "Absolute_Tolerance = $ResidualAbsoluteTolerance"
    # 		# puts $fileid "Residual_Convergence_Tolerance = $ResidualConvergenceTolerance"
    # 		# puts $fileid "Residual_Absolute_Tolerance = $ResidualAbsoluteTolerance"
    # 		# puts $fileid "Displacement_Convergence_Tolerance = $DisplacementConvergenceTolerance"
    # 		# puts $fileid "Displacement_Absolute_Tolerance = $DisplacementAbsoluteTolerance"
    # 	    }
    # 	}
    # 	puts $fileid "Max_Iter = $MaximumIterations"


    # 	# WarnWinText "SolutionMethod:$SolutionMethod ConvergenceCriteria:$ConvergenceCriteria ConvergenceTolerance:$ConvergenceTolerance AbsoluteTolerance:$AbsoluteTolerance MaximumIterations:$MaximumIterations"
    # 	#END NON LINEAR
    # } elseif {$AnalysisType =="Linear"} {

    # 	puts $fileid "Convergence_Criteria  = \"Residual_criteria\""
    # 	puts $fileid "Convergence_Tolerance = 1.0E-4"
    # 	puts $fileid "Absolute_Tolerance    = 1.0E-9"
    # 	puts $fileid "Max_Iter = 10"
    # }

    puts $fileid ""
    puts $fileid "#Constraints Data"
    puts $fileid "#################################################"
    puts $fileid ""
    # Incremental Load
    if {($SolutionType =="Dynamic")||($SolutionType =="Quasi-Static")||($SolutionType =="Pseudo-Dynamic")} {

	set cxpath "$AppId//c.Loads//i.IncrementalLoad"
	set incremental_load [::xmlutils::setXml $cxpath $cproperty]
	if {($incremental_load eq "Yes")} {
	    puts $fileid "Incremental_Load = \"True\""
	} else {
	    puts $fileid "Incremental_Load = \"False\""
	}
	# Incremental Displacement
	set cxpath "$AppId//c.Conditions//i.IncrementalMovement"
	set incremental_movement [::xmlutils::setXml $cxpath $cproperty]
	if {($incremental_movement eq "Yes")} {
	    puts $fileid "Incremental_Displacement = \"True\""
	} else {
	    puts $fileid "Incremental_Displacement = \"False\""
	}

    } else {

	puts $fileid "Incremental_Load = \"False\""
	puts $fileid "Incremental_Displacement = \"False\""

    }

    puts $fileid ""
    puts $fileid "#PostProcess Data"
    puts $fileid "#################################################"
    puts $fileid ""
    # For results

    # On nodes results
    set cnrlist [list "Displacements" "Velocities" "Accelerations" "Rotations" "Reactions" "Forces"]
    set nodal_results "nodal_results=\["
    foreach cnr $cnrlist {
	set cxpath "$AppId//c.Results//c.OnNodes//i.[list ${cnr}]"
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
    }
	append nodal_results "\]"
	  puts $fileid "$nodal_results"
    # WarnWinText "nodal_results:$nodal_results"

    # On Gauss point results
    set cgrlist [list "StrainTensor" "StressTensor" "VonMises" "PlasticStrain" "DeltaPlasticStrain" "BeamMoments" "BeamForces" "ShellForcesLocal" "ShellForcesGlobal" "ShellMomentsLocal" "ShellMomentsGlobal" "ShellStrainLocal" "ShellStrainGlobal" "ShellCurvatureLocal" "ShellCurvatureGlobal" "MaterialDirectionX" "MaterialDirectionY" "MaterialDirectionZ"]
    set gauss_points_results "gauss_points_results=\["
    foreach cgr $cgrlist {
	set cxpath "$AppId//c.Results//c.OnGaussPoints//i.[list ${cgr}]"
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
    }
	append gauss_points_results "\]"
	puts $fileid "$gauss_points_results"
    # WarnWinText "gauss_points_results:$gauss_points_results"

    # GiD post mode variables
    #::wkcf::WriteGiDPostMode $AppId $fileid
    ::wkcf::WriteGiDPostModeNew $AppId $fileid
    puts $fileid ""
    set cxpath "$AppId//c.SolutionStrategy//c.Dynamic//i.DeltaTime"
    set TimeStep [::xmlutils::setXml $cxpath $cproperty]

    # Number of Steps
    set cxpath "$AppId//c.Results//i.OutputDeltaTime"
    set TimeWriting [::xmlutils::setXml $cxpath $cproperty]
    set FrequencyWriting [expr int(double($TimeWriting)/double($TimeStep))]

    if {$FrequencyWriting != 0 && $FrequencyWriting > 1} {
  puts $fileid "GiDWriteFrequency = $TimeWriting"
    } else {
  puts $fileid "GiDWriteFrequency = $TimeStep"
    }
    puts $fileid "WriteResults = \"PreMeshing\""
    puts $fileid "echo_level = 0"
    puts $fileid ""
    puts $fileid "# graph_options"
    puts $fileid "PlotGraphs = \"False\""
    puts $fileid "PlotFrequency = 0 "
    puts $fileid ""
    puts $fileid "# list options"
    puts $fileid "PrintLists = \"True\""
    puts $fileid "file_list = \[\] "
    puts $fileid ""
    puts $fileid "# restart options"
    puts $fileid "SaveRestart = False"
    puts $fileid "RestartFrequency = 0"
    puts $fileid "LoadRestart = False"
    puts $fileid "Restart_Step = 0"


    puts $fileid ""
    set PName [::KUtils::GetPaths "PName"]
    puts $fileid "problem_name=\"${PName}${AppId}\""
    puts $fileid "problem_path=\"[file join $PDir]\""


    # Commented by J. Garate on 17/12/2012
    # Get the kratos path
    # set cxpath "GeneralApplicationData//c.ProjectConfiguration//i.KratosPath"
    # set cproperty "dv"
    # set KratosPath [::xmlutils::setXml $cxpath $cproperty]
    # set KratosPath [file native $KratosPath]

    # Write the kratos path
    # puts $fileid "kratos_path=\"${KratosPath}\""

}



proc ::wkcf::WriteGiDPostModeNew {AppId fileid} {
    # Write the GiD post mode variables for each applications

    # kratos key word xpath
    set kwxpath "Applications/$AppId"

    puts $fileid ""
    set trailing_spaces  "    "
    puts $fileid "# GiD output configuration"
    puts $fileid "class GidOutputConfiguration:"
     # Gid results
    set gidrlist [list "GiDPostMode" "GiDWriteMeshFlag" "GiDWriteConditionsFlag" "GiDWriteParticlesFlag" "GiDMultiFileFlag"]
    foreach gidr $gidrlist {
        # Get the value
        set cxpath "$AppId//c.Results//c.GiDOptions//i.[list ${gidr}]"
        set cproperty "dv"
        set cvalue [::xmlutils::setXml $cxpath $cproperty]
        # Get the kratos keyword
        set gidrkw [::xmlutils::getKKWord $kwxpath $gidr]
        # WarnWinText "gidr:$gidr cvalue:$cvalue gidrkw:$gidrkw"
        if {($gidr=="GiDWriteMeshFlag") || ($gidr=="GiDWriteConditionsFlag") || ($gidr=="GiDWriteParticlesFlag")} {
            if {$cvalue =="Yes"} {
              set cvalue True
            } else {
              set cvalue False
            }
            puts $fileid "${trailing_spaces}$gidrkw = $cvalue"
        } else {
            puts $fileid "${trailing_spaces}$gidrkw = \"$cvalue\""
        }
    }
}



proc ::wkcf::WriteLinearSolver {rootid fileid trailing_spaces config_name} {
    # Write linear solver with a config_name

    # Kratos key word xpath
    set kwxpath "Applications/$rootid"
    # Set default value xml variable
    set cproperty "dv"

    puts $fileid "    class ${config_name}:"

    # Linear Solver type
    set cxpath "$rootid//c.SolutionStrategy//c.LinearSolver//i.LinearSolverType"
    set LinearSolverType [::xmlutils::setXml $cxpath $cproperty]
    if {$LinearSolverType =="Direct"} {

	# Direct solver type
	set cxpath "$rootid//c.SolutionStrategy//c.LinearSolver//i.DirectSolverType"
	set solvertype [::xmlutils::setXml $cxpath $cproperty]
	set DirectSolverType [::xmlutils::getKKWord $kwxpath $solvertype]
	puts $fileid "${trailing_spaces}    solver_type = \"$DirectSolverType\""
	puts $fileid "${trailing_spaces}    scaling = False"

    } elseif {$LinearSolverType =="Iterative"} {

	# Iterative solver type
	set cxpath "$rootid//c.SolutionStrategy//c.LinearSolver//i.IterativeSolverType"
	set solvertype [::xmlutils::setXml $cxpath $cproperty]
	set IterativeSolverType [::xmlutils::getKKWord $kwxpath $solvertype]
	# Tolerance
	set cxpath "$rootid//c.SolutionStrategy//c.LinearSolver//i.Tolerance"
	set Tolerance [::xmlutils::setXml $cxpath $cproperty]

	# Maximum iteration
	set cxpath "$rootid//c.SolutionStrategy//c.LinearSolver//i.MaximumIteration"
	set MaximumIteration [::xmlutils::setXml $cxpath $cproperty]

	puts $fileid "${trailing_spaces}    solver_type = \"$IterativeSolverType\""
	puts $fileid "${trailing_spaces}    tolerance = $Tolerance"
	puts $fileid "${trailing_spaces}    max_iteration = $MaximumIteration"
	puts $fileid "${trailing_spaces}    scaling = False"
	puts $fileid "${trailing_spaces}    verbosity = 0"


        puts $fileid "${trailing_spaces}    #Pastix Iterative Solver:"
        puts $fileid "${trailing_spaces}    gmres_krylov_space_dimension = 100"
        puts $fileid "${trailing_spaces}    ilu_level_of_fill            = 3"

        puts $fileid "${trailing_spaces}    #GMRES or CG:"
	puts $fileid "${trailing_spaces}    preconditioner_type          = \"None\""

        puts $fileid "${trailing_spaces}    #Deflated CG:"
	puts $fileid "${trailing_spaces}    assume_constant_structure    = True"
        puts $fileid "${trailing_spaces}    max_reduced_size             = 1000"

        puts $fileid "${trailing_spaces}    #AMG: (requires block_builder)"
        puts $fileid "${trailing_spaces}    smoother_type  = \"ILU0\" #\"DAMPED_JACOBI\""
	puts $fileid "${trailing_spaces}    krylov_type    = \"GMRES\""

    }

    puts $fileid ""
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
    puts $fileid "from KratosMultiphysics.SolidMechanicsApplication import *"
    # Cross section properties
    puts $fileid "from beam_sections_python_utility import SetProperties"


    puts $fileid "def AssignMaterial(Properties):"

    # For each properties
    set propid 0
    foreach PropertyId $dprops($AppId,GKProps,AllPropertyId) {
	incr propid 1
	puts $fileid "# GUI property identifier: $PropertyId"
	# Get the material identifier for this property
	set MatId $dprops($AppId,Property,$PropertyId,MatId)
	puts $fileid "# GUI material identifier: $MatId"
	# Check for use of the flow rule
	if {$dprops($AppId,Material,$MatId,UseFlowRule) =="Yes"} {
	    #puts $fileid "    FlowRule = $dprops($AppId,Material,$MatId,FlowRule)"
	}
	# Check for use of the hardening law
	if {$dprops($AppId,Material,$MatId,UseHardeningLaw) =="Yes"} {
	    #puts $fileid "    HardeningLaw = $dprops($AppId,Material,$MatId,HardeningLaw)"
	}
	# Check for use of the yield criterion
	if {$dprops($AppId,Material,$MatId,UseYieldCriterion) =="Yes"} {
	    #puts $fileid "    YieldCriterion = $dprops($AppId,Material,$MatId,YieldCriterion)"
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
    set csectiontype "\"Rectangular\""
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
    set mpxpath "[::KMat::findMaterialParent $MatId]//m.[list ${MatId}]"
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
    set cbvalue [lindex [::KMat::getMaterialProperties "p" "$mpxpath//$mbxpath//p.[list $mbehavior]"] 0 1]
    # Get the internal behavior properties
    set mbivalues [split [::xmlutils::getKKWord $clxpath $ptype "mbivalues"] ,]
    # Get the write behavior properties
    set mbwritev [split [::xmlutils::getKKWord $clxpath $ptype "mbwritev"] ,]
    # WarnWinText "mbwritev:$mbwritev mbivalues:$mbivalues\n$mpxpath//$mbxpath//p.[list $mbehavior] cbvalue:$cbvalue"
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
		set dprops($AppId,Material,$MatId,Fluency) "EnergyYieldFunction(State.[list ${mswv}])"
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
		set cstate "myState.[list ${mswv}]"
		break
	    }
	}
	# WarnWinText "cstate:$cstate"
	# Get the yield function properties
	# Get the yield criteria
	set yfid "YieldCriteria"
	set myieldcriteria [::xmlutils::getKKWord "$clxpath" $ptype "myieldcriteria"]
	# Get the yield criteria xpath values
	set mycxpath [::xmlutils::getKKWord "$clxpath" $ptype "mycxpath"]
	# Get the current yield criteria
	set cycvalue [lindex [::KMat::getMaterialProperties "p" "$mpxpath//$mycxpath//p.[list $myieldcriteria]"] 0 1]
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


proc ::wkcf::GetPlasticityProperties {AppId MatId MatModel cptype ptype} {
    # Get the platicity  material properties
    variable dprops
    variable ndime

    # WarnWinText "AppId:$AppId MatId:$MatId cptype:$cptype ptype:$ptype"
    # Xpath for constitutive laws
    set clxpath "CLawProperties"
    # Get all material properties
    set mpxpath "[::KMat::findMaterialParent $MatId]//m.[list ${MatId}]"
    # WarnWinText "mpxpath:$mpxpath"

    # Set global FlowRule, HardeningLaw and YieldCriterion
    set dprops($AppId,Material,$MatId,UseFlowRule) "Yes"
    set dprops($AppId,Material,$MatId,UseYieldCriterion) "Yes"
    set dprops($AppId,Material,$MatId,UseHardeningLaw) "Yes"


    if {$MatModel == "HyperElastic-Plastic"} {

	# HyperElastic-plastic models

	# Get the FLOW RULE
	set mflowrule [::xmlutils::getKKWord $clxpath $ptype "mflowrule"]
	# Get the hardning law xpath values
	set mfrxpath [::xmlutils::getKKWord $clxpath $ptype "mfrxpath"]
	# WarnWinText "mflowrule:$mflowrule mfrxpath:$mfrxpath"
	# Get the current flow rule
	set cfrvalue [lindex [::KMat::getMaterialProperties "p" "$mpxpath//$mfrxpath//p.[list $mflowrule]"] 0 1]
	# Get the internal flow rule
	set mfrivalues [split [::xmlutils::getKKWord $clxpath $ptype "mfrivalues"] ,]
	# Get the writen internal flow rule properties
	set mfrwritev [split [::xmlutils::getKKWord $clxpath $ptype "mfrwritev"] ,]
	# WarnWinText "mfrwritev:$mfrwritev mfrivalues:$mfrivalues\n$mpxpath//$mfrxpath//p.[list $mflowrule] cfrvalue:$cfrvalue"
	foreach mfriv $mfrivalues mfrwv $mfrwritev {
	    if {$mfriv ==$cfrvalue} {
		set dprops($AppId,Material,$MatId,FlowRule) "$mfrwv"
		break
	    }
	}
	# WarnWinText "dprops($AppId,Material,$MatId,FlowRule):$dprops($AppId,Material,$MatId,FlowRule)"


	# Get the HARDENING LAW
	set hlid "HardeningLaw"
	set mhardeninglaw [::xmlutils::getKKWord $clxpath $ptype "mhardeninglaw"]
	# Get the hardning law xpath values
	set mhxpath [::xmlutils::getKKWord $clxpath $ptype "mhxpath"]
	# WarnWinText "mhardeninglaw:$mhardeninglaw mhxpath:$mhxpath"
	# Get the current hardening law
	set chvalue [lindex [::KMat::getMaterialProperties "p" "$mpxpath//$mhxpath//p.[list $mhardeninglaw]"] 0 1]

	# WarnWinText "mhardeninglaw:$mhardeninglaw mhxpath:$mhxpath chvalue:$chvalue"
	# Get the yield criterion options
	set mhivalues [split [::xmlutils::getKKWord "$clxpath//$hlid" "AvailableHardeningLaw" "mhivalues"] ,]
	# Get the write yield criterion properties
	set mhwritev [split [::xmlutils::getKKWord "$clxpath//$hlid" "AvailableHardeningLaw" "mhwritev"] ,]
	# WarnWinText "mhwritev:$mhwritev mhivalues:$mhivalues\n$mpxpath//$mhxpath//p.[list $mhardeninglaw] chvalue:$chvalue"
	set chl ""
	foreach mhiv $mhivalues mhwv $mhwritev {
	    # WarnWinText "mhiv:$mhiv mhivalue:$mhivalue mhwv:$mhwv"
	    if {$mhiv ==$chvalue} {
		set chl "$mhwv"
		break
	    }
	}
	# WarnWinText "chl:$chl"
	if {$chl !=""} {
	    set dprops($AppId,Material,$MatId,HardeningLaw) "$chl"

	    # Update material properties

	}
	# WarnWinText "dprops($AppId,Material,$MatId,HardeningLaw):$dprops($AppId,Material,$MatId,HardeningLaw)"

	# Get the SATURATION LAW
	set cslid "SaturationLaw"
	set msaturationlaw [::xmlutils::getKKWord $clxpath $ptype "msaturationlaw"]
	# Get the saturation law xpath values
	set msxpath [::xmlutils::getKKWord $clxpath $ptype "msxpath"]
	# WarnWinText "msaturationlaw:$msaturationlaw msxpath:$msxpath"
	# Get the current hardening law
	set csvalue [lindex [::KMat::getMaterialProperties "p" "$mpxpath//$msxpath//p.[list $msaturationlaw]"] 0 1]
	# WarnWinText "msaturationlaw:$msaturationlaw msxpath:$msxpath csvalue:$csvalue"
	# Get the yield criterion options
	set msivalues [split [::xmlutils::getKKWord "$clxpath//$cslid" "AvailableSaturationLaws" "msivalues"] ,]
	# Get the write yield criterion properties
	set mswritev [split [::xmlutils::getKKWord "$clxpath//$cslid" "AvailableSaturationLaws" "mswritev"] ,]
	# WarnWinText "mswritev:$mswritev mhivalues:$msivalues\n$mpxpath//$msxpath//p.[list $msaturationlaw] csvalue:$csvalue"
	set csl ""
	foreach msiv $msivalues mswv $mswritev {
	    # WarnWinText "msiv:$msiv msivalue:$msiv mswv:$mswv"
	    if {$msiv ==$csvalue} {
		set csl "$mswv"
		break
	    }
	}
	# WarnWinText "csl:$csl"

	set dprops($AppId,Material,$MatId,SaturationLaw) "$csl"

	# WarnWinText "dprops($AppId,Material,$MatId,SaturationLaw):$dprops($AppId,Material,$MatId,SaturationLaw)"

	# Get the yield function properties

	# Get the YIELD CRITERION
	set ycid "YieldCriteria"
	set myieldcriterion [::xmlutils::getKKWord "$clxpath" $ptype "myieldcriterion"]
	# Get the yield criteria xpath values
	set mycxpath [::xmlutils::getKKWord "$clxpath" $ptype "mycxpath"]
	# Get the current yield criterion
	set cycvalue [lindex [::KMat::getMaterialProperties "p" "$mpxpath//$mycxpath//p.[list $myieldcriterion]"] 0 1]
	# WarnWinText "myieldcriterion:$myieldcriterion mycxpath:$mycxpath cycvalue:$cycvalue"
	# Get the yield criterion options
	set ycivalues [split [::xmlutils::getKKWord "$clxpath//$ycid" "AvailableYieldCriteria" "ycivalues"] ,]
	# Get the write yield criterion properties
	set ycwritev [split [::xmlutils::getKKWord "$clxpath//$ycid" "AvailableYieldCriteria" "ycwritev"] ,]
	# WarnWinText "ycwritev:$ycwritev ycivalues:$ycivalues"
	set cyc ""
	foreach yciv $ycivalues ycwv $ycwritev {
	    # WarnWinText "yciv:$yciv cycvalue:$cycvalue ycwv:$ycwv"
	    if {$yciv ==$cycvalue} {
		set cyc "$ycwv"
		break
	    }
	}
	# WarnWinText "cyc:$cyc"

	set dprops($AppId,Material,$MatId,YieldCriterion) "$cyc"
    }
    # WarnWinText "dprops($AppId,Material,$MatId,YieldCriterion):$dprops($AppId,Material,$MatId,YieldCriterion)"
}

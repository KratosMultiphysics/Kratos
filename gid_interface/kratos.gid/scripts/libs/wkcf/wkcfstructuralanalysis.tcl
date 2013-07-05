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
#     1.7- 25/06/13-A. Melendo, new proc WriteFaceLoads
#     1.6- 17/06/13-G. Socorro, delete wmethod variable and all relalated procedure (*_m0,*_m1,*_m2) => now we are using only the new GiD groups
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
                "Face" {
                    # Face loads
                    ::wkcf::WriteFaceLoads $AppId $cloadtid
                }
            }
            }
        }
    }
}

proc ::wkcf::WriteBeamUniformlyDistributedLoads {AppId cloadtid} {
    # ASBTRACT: Write beam uniformly distributed loads (forces)
    variable ndime;  variable dprops
    variable filechannel;  variable sa_icondid 
    
    # For debug
    if {!$::wkcf::pflag} {
        set inittime [clock seconds]
    }
    
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
       
    # For debug
    if {!$::wkcf::pflag} {
        set endtime [clock seconds]
        set ttime [expr $endtime-$inittime]
        # WarnWinText "endtime:$endtime ttime:$ttime"
        WarnWinText "Write structural analysis beam uniformly distributed loads (forces): [::KUtils::Duration $ttime]"
    }
}

proc ::wkcf::WritePressureLoads {AppId cloadtid} {
    # ABSTRACT: Write pressure loads (positive or negative shell pressure)
    variable ndime;  variable dprops
    variable filechannel; variable sa_icondid
       
    # For debug
    if {!$::wkcf::pflag} {
        set inittime [clock seconds]
    }

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
    
    # For debug
    if {!$::wkcf::pflag} {
        set endtime [clock seconds]
        set ttime [expr $endtime-$inittime]
        # WarnWinText "endtime:$endtime ttime:$ttime"
        WarnWinText "Write structural analysis write pressure loads (positive or negative shell pressure): [::KUtils::Duration $ttime]"
    }
}

proc ::wkcf::WriteFaceLoads {AppId cloadtid} {
    # ABSTRACT: Write face loads (positive or negative shell pressure)
    variable ndime;  variable dprops
    variable filechannel; variable sa_icondid
       
    # For debug
    if {!$::wkcf::pflag} {
        set inittime [clock seconds]
    }

    # Kratos key word xpath
    set kxpath "Applications/$AppId"
    # Set the GiD element type
    set GiDElemType "Triangle"
    # Set the face 3d-3n condition keyword
    set faceforce3d3nkword "FaceForce3D3N"
    # Set the reference property id
    set RefPropId "1"
       
    # For all defined group identifier inside this load type
    foreach cgroupid $dprops($AppId,Loads,$cloadtid,AllGroupId) {
        if {[GiD_EntitiesGroups get $cgroupid elements -count -element_type $GiDElemType]} {
            # Write Face3D3N condition
            GiD_File fprintf $filechannel "%s" "Begin Conditions $faceforce3d3nkword // GUI face load group identifier: $cgroupid"
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
        if {([GiD_EntitiesGroups get $cgroupid nodes -count]>0)&&($ndime=="3D")} {
            
            foreach Component {BUDFx BUDFy BUDFz} {
              # Get the current face load component keyword                       
              set kwid "$Component"
              set ckword [::xmlutils::getKKWord $kxpath $kwid]              
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

proc ::wkcf::WritePuntualLoads {AppId cloadtid kwordlist} {
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

    # For debug
    if {!$::wkcf::pflag} {
        set endtime [clock seconds]
        set ttime [expr $endtime-$inittime]
        # WarnWinText "endtime:$endtime ttime:$ttime"
        WarnWinText "Write structural analysis concentrated loads (puntual force or moment): [::KUtils::Duration $ttime]"
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
	append nodal_results "\]" 
	puts $fileid "$nodal_results"
    }
    # WarnWinText "nodal_results:$nodal_results"

    # On Gauss point results
    #puts $fileid ""
    set cgrlist [list "GreenLagrangeStrainTensor" "Rotations" "PK2StressTensor" "BeamMoments" "BeamForces"]
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
	set yfid "YieldFunctions"
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

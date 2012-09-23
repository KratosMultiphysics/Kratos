###############################################################################
#
#    NAME: wkcf.tcl
#
#    PURPOSE: Create the base namespace to write all the Kratos calculation files
#             Get all the necesary data from the spd in xml format
#
#    QUANTECH ATZ-DEVELOPMENT DEPARTMENT
#
#    AUTHOR : G. Socorro
#
#    CREATED AT: 01/11/09
#
#    HISTORY:
#   
#     5.2- 23/09/12-G. Socorro, Write the GroupMeshProperties before  the cut and graph properties to get the mesh reference
#     5.1- 21/09/12-G. Socorro, write bat file only when ::tcl_platform(os) eq "Linux"
#     5.0- 07/06/12-G. Socorro, modify the proc WriteConditions to write WallCondition2D/WallCondition3D for fractional step solver and 
#                               MonolithicWallCondition2D/MonolithicWallCondition3D for monolithic solver
#     4.9- 04/06/12-J. Garate,  Correct "WALL_LAW_Y"
#     4.8- 13/05/12-G. Socorro, create a new variable to link the conditions with others BC applied are the boundary of the body
#     4.7- 10/05/12-G. Socorro, comment a warning message in the proc WriteConditions
#     4.6- 09/05/12-G. Socorro, update the procedure WriteElementConnectivities to write the format for each gid element type
#     4.5- 08/05/12-G. Socorro, change IS-SLIP keyword by IS_STRUCTURE, propose by Riccardo Rossi, update the proc WriteConditions write condid from 1 to n
#     4.4- 07/05/12-G. Socorro, Modify the proc WriteConditions, improve others procs ([write_calc_data) 
#     4.3- 16/04/12-G. Socorro, New procedure to Write Coordinates.
#     4.2- 12/04/12-G. Socorro, New procedure to Write Connectivities.
#     4.1- 03/04/12-G. Socorro, modify the proc WriteElementConnectivities to use the new variable "usennode" from the kratos keyword file
#     4.0- 02/04/12-G. Socorro, split the source code into many scripts
#     3.9- 30/04/12-G. Socorro, add the proc WriteCutAndGraph to write the cut and graph properties
#     3.8- 23/02/12 J. Gárate, Hay que revisarlo. (line 1440) $allnlist está vacía, que no trabaje con ella
#     3.7- 10/02/12 J. Gárate, Corregidos los bugs de ROTATIONS y de ANALYSIS TYPE LINEAR, Tabulacion. Cerrar $fileid.
#     3.6- 08/06/11-G. Socorro, write always the body force in the property block
#     3.5- 07/06/11-G. Socorro, correct a bug when write material and section properties
#     3.4- 06/06/11-G. Socorro, add "OutputDeltaTime" to the structural analysis application project parameter
#     3.3- 26/05/11-G. Socorro, add RelaxedDynamic options
#     3.2- 24/05/11-G. Socorro, add the FACE3D3N condition when use surface pressure load type
#     3.1- 23/05/11-G. Socorro, modify some procedure to use the membrane element 
#     3.0- 03/02/11-G. Socorro, comment the procedure WritePythonGroupProperties until the Python group design is finished
#     2.9- 02/02/11-G. Socorro, add a new variable to debug the source code (written process for big models)
#     2.8- 25/01/11-G. Socorro, correct a bug in the procedure WriteFluidBC (write duplicate node identifier and repeated the y values in the z values)
#     2.7- 02/11/10-G. Socorro, modify the procedure WritePropertyAtNodes to use dict option => fast write option for big models
#     2.6- 23/09/10-G. Socorro, add Start_time variable
#     2.5- 08/09/10-G. Socorro, add "Distance" variable to the fluid application results
#     2.4- 06/09/10-G. Socorro, reset xcomp, ycomp and zcomp at the end of the foreach group when write no-slip BC in the fluid application  
#     2.3- 06/09/10-G. Socorro, correct an error when defined body force for group of element with the same property
#     2.2- 03/09/10-G. Socorro, add monolithic and pressuresplitting solver option, change "True" by true
#     2.1- 05/06/10-G. Socorro, add quotes to the Solution_method variable
#     2.0- 30/06/10-G. Socorro, update fluid parameter file
#     1.9- 28/06/10-G. Socorro, add GiDOptions variable in the parameter file
#     1.8- 25/06/10-G. Socorro, improve the procedure WriteNodalCoordinates to write all node identifier without repeated id and in increasing order
#     1.7- 22/06/10-G. Socorro, update project parameters script
#     1.6- 11/06/10-G. Socorro, update some procedure to include the new constitutive load definition
#     1.5- 09/06/10-G. Socorro, add the procedure WritePythonGroupProperties and update some proc to use a FSI application
#     1.4- 08/06/10-G. Socorro, add FSI application
#     1.3- 31/05/10-G. Socorro, get all kratos key word from a XML file (kratos_key_words.xml)
#     1.2- 20/05/10-G. Socorro, modify some function to include the Fluid application
#     1.1- 18/05/10-G. Socorro, modify the procedure SelectPythonScript to select the script as a solution function
#     1.0- 12/05/10-G. Socorro, start to write the data file for the Fluid application
#     0.9- 11/05/10-G. Socorro, select the correct Python script and update the project parameters file
#     0.8- 10/05/10-G. Socorro, pass some utilities procedure to the script wkcfutils.tcl
#     0.7- 05/05/10-G. Socorro, start to write the Kratos global properties
#     0.6- 04/05/10-G. Socorro, write nodal displacement condition for structural analysis application
#     0.5- 03/05/10-G. Socorro, write nodal coordinates block
#     0.4- 30/04/10-G. Socorro, write element properties block
#     0.3- 27/04/10-G. Socorro, write the project parameters file and the run scripts
#     0.2- 20/04/10-G. Socorro, add sqlite3 database options
#     0.1- 01/11/09-G. Socorro, create a base source code
#
###############################################################################

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
    
    # 0 => Metodo antiguo (Poco eficiente)
    # 1 => Metodo nuevo (write_calc_data)
    variable wmethod 
}

proc ::wkcf::WriteCalculationFiles {filename} {
    global KPriv     
    variable FluidApplication; variable StructuralAnalysis
    variable ActiveAppList
    
    # Unset some local variables
    ::wkcf::UnsetLocalVariables
    
    # Init some namespace global variables
    ::wkcf::Preprocess
  
    # Write each block of the files *.mdpa
    # Rename the file name => Change .dat by .mdpa
    set basefilename "[string range $filename 0 end-4]"
    # Write a mdpa file for each application
    foreach AppId $ActiveAppList {
	set filename "${basefilename}${AppId}.mdpa"
	# WarnWinText "filename:$filename"
	# Use the write_calc_data procedure from the GiD kernel
	# Init
	write_calc_data init $filename

	# Write model part data
	::wkcf::WriteModelPartData $AppId

	# Write properties block
	::wkcf::WriteProperties $AppId
	
	# Write nodes block    
	::wkcf::WriteNodalCoordinates $AppId

	# Write elements block 
	::wkcf::WriteElementConnectivities $AppId

	# For fluid application
	if {$AppId == "Fluid"} {
	    # Write conditions (Condition2D and Condition3D)
	    ::wkcf::WriteConditions $AppId  
	}

	# Write boundary condition block
	::wkcf::WriteBoundaryConditions $AppId
	
	# For structural analysis application
	if {$AppId =="StructuralAnalysis"} {
	    # Write load properties block
	    ::wkcf::WriteLoads $AppId
	}

	# For fluid application
	if {$AppId == "Fluid"} {
	    # Write nodal data for density and viscosity            
	    ::wkcf::WritePropertyAtNodes $AppId
	    
	    # Write all group properties
	    # ::wkcf::WriteGroupProperties $AppId
	    ::wkcf::WriteGroupMeshProperties $AppId

	    # Write the cutting and point history properties
	    ::wkcf::WriteCutAndGraph $AppId
	}


	# End
	write_calc_data end
    }
    
    # Write the project parameters file
    ::wkcf::WriteProjectParameters
    
    # Write python group properties
    # ::wkcf::WritePythonGroupProperties
    
    # Select python scripts
    ::wkcf::SelectPythonScript
    
    # Write constitutive laws properties
    if {$StructuralAnalysis=="Yes"} {
	::wkcf::WriteConstitutiveLawsProperties
    }

    # Write bat file only for the Linux OS
    if {($::tcl_platform(os) eq "Linux")} {
	::wkcf::WriteBatFile $AppId
    }
    
    # Unset some local variables
    ::wkcf::UnsetLocalVariables
    
    return 1
}

proc ::wkcf::SelectPythonScript {} {
    # Select the correct Python script
    variable FluidApplication; variable StructuralAnalysis
    
    set endfilename "KratosOpenMP.py"
    set mpiendfilename "KratosMPI.py"
    
    # Get the application root identifier    
    set rootdataid [::wkcf::GetApplicationRootId]
    
    if {$StructuralAnalysis =="Yes"} {
	# Solution type
	set cproperty "dv"
	set cxpath "$rootdataid//c.AnalysisData//i.SolutionType"
	set SolutionType [::xmlutils::setXml $cxpath $cproperty]
	# WarnWinText "SolutionType:$SolutionType"
	if {($SolutionType =="Dynamic")||($SolutionType =="RelaxedDynamic")} {
	    set ppfilename "KratosOpenMPStructuralDynamic.py"
	    set mpifilename "KratosMPIStructuralDynamic.py"
	} elseif {$SolutionType =="Static"} {
	    set ppfilename "KratosOpenMPStructuralStatic.py"
	    set mpifilename "KratosMPIStructuralStatic.py"
	}
    }
    if {$FluidApplication =="Yes"} {
	set ppfilename "KratosOpenMPFluid.py"
	set mpifilename "KratosMPIFluid.py"
    }
    
    set PTDir [::KUtils::GetPaths "PTDir"]
    set fromfname [file native [file join "$PTDir/python" $ppfilename]]
    set mpifromfname [file native [file join "$PTDir/python" $mpifilename]]
    
    set PDir [file native [::KUtils::GetPaths "PDir"]]
    set tofname [file native [file join $PDir $endfilename]]
    set mpitofname [file native [file join $PDir $mpiendfilename]]
    
    # Copy the script file
    if {[catch {file copy -force "$fromfname" "$tofname"} error]} {
	WarnWin [= "Could not copy the Kratos Python script (%s) to (%s): Error (%)" $fromfname $tofname $error ]
	return ""
    }
    if {[catch {file copy -force "$mpifromfname" "$mpitofname"} error]} {
	WarnWin [= "Could not copy the Kratos Python script (%s) to (%s): Error (%)" $mpifromfname $mpitofname $error ]
	return ""
    }    
    
}

proc ::wkcf::WriteModelPartData {AppId} {
    # Write the model part data
    # Arguments
    # AppId => Application identifier
    
    write_calc_data puts "Begin ModelPartData"
    write_calc_data puts "//  VARIABLE_NAME value"
    write_calc_data puts "End ModelPartData"
    write_calc_data puts ""
    
}

proc ::wkcf::WriteProperties {AppId} {
    # Write the properties block
    # Arguments
    # AppId => Application identifier
    variable dprops;  variable ndime
    
    # For structural analysis application
    if {$AppId =="StructuralAnalysis"} {
	set propid 0
	foreach PropertyId $dprops($AppId,GKProps,AllPropertyId) {
	    incr propid 1
	    write_calc_data puts "Begin Properties $propid // GUI property identifier: $PropertyId"
	    # Get the material identifier for this property 
	    set MatId $dprops($AppId,Property,$PropertyId,MatId) 
	    # WarnWinText "PropertyId:$PropertyId MatId:$MatId"
	    write_calc_data puts "// GUI material identifier: $MatId"
	    # Write material properties
	    foreach cpropid $dprops($AppId,Material,$MatId,Props) {
		# WarnWinText "material propid:$cpropid"
		lassign $cpropid key value
		write_calc_data puts " $key $value"
	    } 
	    
	    # Write section (others) properties (thickness, etc)
	    if {[info exists dprops($AppId,Material,$PropertyId,CProps)]} {
		foreach cpropid $dprops($AppId,Material,$PropertyId,CProps) {
		    # WarnWinText "cpropid:$cpropid"
		    lassign $cpropid key value
		    # WarnWinText "key:$key value:$value"
		    write_calc_data puts " $key $value"
		} 
	    }
	    
	    if {$dprops($AppId,GKProps,$PropertyId,AddBF)=="Yes"} {
		# Get the body force properties
		set cloadtid "BodyForce"
		# Get the group identifier list for this property
		if {([info exists dprops($AppId,Property,$PropertyId,GroupId)]) && ([llength $dprops($AppId,Property,$PropertyId,GroupId)]>0)} {
		    foreach GroupId $dprops($AppId,Property,$PropertyId,GroupId) {
			# WarnWinText "GroupId:$GroupId"
			if {[info exists dprops($AppId,Loads,$cloadtid,$GroupId,GProps)]} {
			    set cprop $dprops($AppId,Loads,$cloadtid,$GroupId,GProps)
			    # WarnWinText "cprop:$cprop GroupId:$GroupId"
			    write_calc_data puts "// GUI body force group identifier: $GroupId"
			    ::wkcf::WriteBodyForceValues $cprop
			}
		    }
		}
	    } else {
		# Write the body forces with default values
		set cprop [list [list GravityValue 9.8] [list Cx 0.0] [list Cy 0.0] [list Cz 0.0]]
		::wkcf::WriteBodyForceValues $cprop
	    }
	    
	    write_calc_data puts "End Properties"
	}
	write_calc_data puts ""
    }
    
    # For fluid application
    if {$AppId =="Fluid"} {
	write_calc_data puts "Begin Properties 0"
	write_calc_data puts "End Properties"
	write_calc_data puts ""
    }
    
}

proc ::wkcf::WriteNodalCoordinates {AppId} {
    # Write the nodal coordinates block
    # Nodes block format
    # Begin Nodes
    # // id          X        Y        Z
    # End Nodes
    # Arguments
    # AppId  => Application identifier
    variable ndime; variable dprops; variable wmethod
    
    
    # Check for all defined kratos elements
    if {([info exists dprops($AppId,AllKElemId)]) && ([llength $dprops($AppId,AllKElemId)])} {
	# For debug
	if {!$::wkcf::pflag} {
	    set inittime [clock seconds]
	}
	if {$wmethod eq "1"} {
	    set cformat "%10d"
	    # For all defined kratos nodes        
	    foreach celemid $dprops($AppId,AllKElemId) {
		# Check for all defined group identifier for this element
		if {([info exists dprops($AppId,KElem,$celemid,AllGroupId)]) && ([llength $dprops($AppId,KElem,$celemid,AllGroupId)])} {
		    # For all defined group identifier for this element
		    foreach cgroupid $dprops($AppId,KElem,$celemid,AllGroupId) {
			# Get the GiD entity type, element type and property identifier
			lassign $dprops($AppId,KElem,$celemid,$cgroupid,GProps) GiDEntity GiDElemType PropertyId KEKWord nDim
			# WarnWinText "GiDEntity:$GiDEntity GiDElemType:$GiDElemType PropertyId:$PropertyId KEKWord:$KEKWord nDim:$nDim"
			# Create the dictionary used in the format
			set gprop [dict create]	
			dict set gprop $cgroupid "$cformat"
			
			if {[write_calc_data nodes -count $gprop]>0} {
			    # Write all nodes for this group in incresing orden
			    write_calc_data puts "Begin Nodes \/\/ GUI group identifier: $cgroupid"
			    set wbuff ""
			    if {$ndime =="2D"} {
				foreach nodeid [write_calc_data nodes -sorted -return $gprop] {
				    lassign  [lrange [lindex [GiD_Info Coordinates $nodeid mesh] 0] 0 2] xval yval
				    # WarnWinText "nodeid:$nodeid ncoord:$ncoord"
				    append wbuff "[format "%4i  %10.5f  %10.5f  %10.5f" $nodeid $xval $yval 0.0]\n"
				    
				} 
			    } else {
				foreach nodeid [write_calc_data nodes -sorted -return $gprop] {
				    lassign  [lrange [lindex [GiD_Info Coordinates $nodeid mesh] 0] 0 2] xval yval zval
				    # WarnWinText "nodeid:$nodeid ncoord:$ncoord"
				    append wbuff "[format "%4i  %10.5f  %10.5f  %10.5f" $nodeid $xval $yval $zval]\n"
				}
			    }
			    write_calc_data puts "[string trimright $wbuff]"
			}
			write_calc_data puts "End Nodes"
			write_calc_data puts ""
			# Unset the dictionary
			unset gprop
		    }
		}
	    }
	} else {
	    set cnodeglist [list]
	    # Create a dictionary
	    set nc [dict create 0 0]
	    # For all defined kratos nodes        
	    foreach celemid $dprops($AppId,AllKElemId) {
		# Check for all defined group identifier for this element
		if {([info exists dprops($AppId,KElem,$celemid,AllGroupId)]) && ([llength $dprops($AppId,KElem,$celemid,AllGroupId)])} {
		    # For all defined group identifier for this element
		    foreach cgroupid $dprops($AppId,KElem,$celemid,AllGroupId) {
			# Get the GiD entity type, element type and property identifier
			lassign $dprops($AppId,KElem,$celemid,$cgroupid,GProps) GiDEntity GiDElemType PropertyId KEKWord nDim
			# WarnWinText "GiDEntity:$GiDEntity GiDElemType:$GiDElemType PropertyId:$PropertyId KEKWord:$KEKWord nDim:$nDim"
			# Get all defined entities for this group identifier
			set allelist [::KUtils::GetDefinedMeshGiDEntities $cgroupid $GiDEntity]
			# WarnWinText "alllist:$allelist"
			if {[llength $allelist]} {
			    # Non repeated nodes identifier
			    # Get all defined nodes
			    foreach elemid $allelist {
				# Get the element properties
				foreach nodeid [lrange [GiD_Info Mesh Elements $GiDElemType $elemid] 1 end-1] {
				    dict set nc $nodeid $cgroupid 
				}
			    }
			}
		    }
		}
	    }
	    # Write the group nodal properties
	    foreach celemid $dprops($AppId,AllKElemId) {
		# Check for all defined group identifier for this element
		if {([info exists dprops($AppId,KElem,$celemid,AllGroupId)]) && ($dprops($AppId,KElem,$celemid,AllGroupId)>0)} {
		    # For all defined group identifier for this element
		    foreach cgroupid $dprops($AppId,KElem,$celemid,AllGroupId) {
			# Init the node group list
			set cnodeglist [list]
			dict for {nodeid dgroupid} $nc {
			    if {$cgroupid ==$dgroupid} {
				lappend cnodeglist $nodeid
			    }
			}
			if {[llength $cnodeglist]} {
			    # Write all nodes for this group in incresing orden
			    write_calc_data puts "Begin Nodes \/\/ GUI group identifier: $cgroupid"
			    set wbuff ""
			    foreach nodeid [lsort -integer $cnodeglist] {
				set ncoord [lrange [lindex [GiD_Info Coordinates $nodeid mesh] 0] 0 2]
				# WarnWinText "nodeid:$nodeid ncoord:$ncoord"
				if {$ndime =="2D"} {
				    append wbuff "[format "%4i  %10.5f  %10.5f  %10.5f" $nodeid [lindex $ncoord 0] [lindex $ncoord 1] 0.0]\n"
				} elseif {$ndime =="3D"} {
				    append wbuff "[format "%4i  %10.5f  %10.5f  %10.5f" $nodeid [lindex $ncoord 0] [lindex $ncoord 1] [lindex $ncoord 2]]\n"
				}
			    }
			    write_calc_data puts "[string trimright $wbuff]"
			}
			write_calc_data puts "End Nodes"
			write_calc_data puts ""
			unset cnodeglist 
		    }
		}
	    }
	    unset nc
	    write_calc_data puts ""
	}
	# For debug
	if {!$::wkcf::pflag} {
	    set endtime [clock seconds]
	    set ttime [expr $endtime-$inittime]
	    # WarnWinText "endtime:$endtime ttime:$ttime"
	    WarnWinText "Write nodal coordinates: [::KUtils::Duration $ttime]"
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
    variable useqelem; variable dprops; variable wmethod
    global KPriv
    
    # Check for all defined kratos elements
    if {([info exists dprops($AppId,AllKElemId)]) && ([llength $dprops($AppId,AllKElemId)])} {
	# For debug
	if {!$::wkcf::pflag} {
	    set inittime [clock seconds]
	}
	set kwxpath "Applications/$AppId"
	# For all defined kratos elements        
	foreach celemid $dprops($AppId,AllKElemId) {
	    # Check for all defined group identifier for this element
	    if {([info exists dprops($AppId,KElem,$celemid,AllGroupId)]) && ([llength $dprops($AppId,KElem,$celemid,AllGroupId)])} {
		# WarnWinText "celemid:$celemid"
		# For all defined group identifier for this element
		foreach cgroupid $dprops($AppId,KElem,$celemid,AllGroupId) {
		    # Get the GiD entity type, element type and property identifier
		    lassign $dprops($AppId,KElem,$celemid,$cgroupid,GProps) GiDEntity GiDElemType PropertyId KEKWord nDim
		    # WarnWinText "GiDEntity:$GiDEntity GiDElemType:$GiDElemType PropertyId:$PropertyId KEKWord:$KEKWord nDim:$nDim"
		    
		    # Get the element identifier => nDimNNode
		    set etbf ""
		    # Get use ndimNNode from the kratos keyword file
		    set usennode [::xmlutils::getKKWord $kwxpath "$celemid" "usennode"]
		    # wa "usennode:$usennode elemkword:$KEKWord"
		    if {$usennode eq "Yes"} {
			set etbf [::wkcf::GetnDimnNode $GiDElemType $nDim]
		    } 
		    if {$wmethod ne "1"} {
			# Get all defined entities for this group identifier
			set allelist [::KUtils::GetDefinedMeshGiDEntities $cgroupid $GiDEntity]
			# WarnWinText "alllist:$allelist"
			if {[llength $allelist]} {
			    set kelemtype "[string trim ${KEKWord}${etbf}]"
			    set GlobalPId $dprops($AppId,KElem,$celemid,$cgroupid,GlobalPId)
			    write_calc_data puts "Begin Elements $kelemtype   \/\/ GUI group identifier: $cgroupid"
			    foreach elemid $allelist {
				# WarnWinText "elemid:$elemid"
				# Get the element properties
				set eprop [GiD_Info Mesh Elements $GiDElemType $elemid]
				# WarnWinText "eprop:$eprop"
				if {[llength $eprop]} {
				    set elemf ""
				    append elemf [format "%8i%8i" $elemid $GlobalPId]
				    foreach citem [lrange $eprop 1 end-1] {
					append elemf [format "%8i" $citem] " "
				    }
				    write_calc_data puts "$elemf"
				}
			    }
			}
			write_calc_data puts "End Elements"
			write_calc_data puts ""
		    } else {
			# Create the dictionary
			set gprop [dict create]
			set f "%10d"
			set f [subst $f]
			dict set gprop $cgroupid "$f"
			if {[write_calc_data has_elements -elemtype $GiDElemType $gprop]} { 
			    set kelemtype "[string trim ${KEKWord}${etbf}]"
			    set GlobalPId $dprops($AppId,KElem,$celemid,$cgroupid,GlobalPId)
			    write_calc_data puts "Begin Elements $kelemtype   \/\/ GUI group identifier: $cgroupid"
			    switch -exact -- $GiDElemType {
				"Triangle" {
				    set f "%10d [format "%10d" $GlobalPId] %10d %10d %10d\n"
				}
				"Quadrilateral" {
				    set f "%10d [format "%10d" $GlobalPId] %10d %10d %10d %10d\n"
				}
				"Tetrahedra" {
				    set f "%10d [format "%10d" $GlobalPId] %10d %10d %10d %10d\n"
				}
				"Hexahedra" {
				    set f "%10d [format "%10d" $GlobalPId] %10d %10d %10d %10d %10d %10d %10d %10d\n"
				}
			    }
			    # wa "f:$f"
			    set f [subst $f]
			    dict set gprop $cgroupid "$f"
			    # wa "GiDElemType:$GiDElemType"
			    write_calc_data connectivities -elemtype "$GiDElemType" $gprop
			    write_calc_data puts "End Elements"
			    write_calc_data puts ""
			}
			# Unset the dictionary
			unset gprop
		    }
		}
	    }
	}
	write_calc_data puts ""
	
	# For debug
	if {!$::wkcf::pflag} {
	    set endtime [clock seconds]
	    set ttime [expr $endtime-$inittime]
	    # WarnWinText "endtime:$endtime ttime:$ttime"
	    WarnWinText "Write element connectivities: [::KUtils::Duration $ttime]"
	}
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
			set shelllist2 [list "ShellIsotropic" "ShellAnisotropic" "EBST"]
			if {([info exists dprops($AppId,AllKElemId)]) && ($dprops($AppId,AllKElemId)>0)} {
			    # For all defined kratos elements        
			    foreach celemid $dprops($AppId,AllKElemId) {
				if {$celemid in $shelllist2} {
				    set useshells2 "Yes"
				}
			    }
			}
			if {$useshells2 == "Yes"} {
			    set kwxpath "Applications/StructuralAnalysis"
			    set kwordlist [list [::xmlutils::getKKWord $kwxpath "Rx"] [::xmlutils::getKKWord $kwxpath "Ry"] [::xmlutils::getKKWord $kwxpath "Rz"]]
			    # Process displacement properties
			    ::wkcf::WriteDispRotBC $AppId $ccondid $kwordlist 
			}
		    }
		    "OutletPressure" {
			set kwordlist [list "PRESSURE"]
			# Process outlet pressure
			::wkcf::WriteOutLetPressureBC $AppId $ccondid $kwordlist 
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
		    "Is-Slip" {
			# Write is-slip boundary condition
			# set kwordlist [list "IS-SLIP"]
			set kwordlist [list "IS_STRUCTURE" "Y_WALL"]
			::wkcf::WriteFluidIsSlipBC $AppId $ccondid $kwordlist
		    }
		    "WallLaw" {
			# Write wall law boundary condition
			set kwordlist [list "WALL_LAW_Y"]
			::wkcf::WriteFluidWallLawBC $AppId $ccondid $kwordlist
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

proc ::wkcf::WriteConditions {AppId} {
    # ABSTRACT: Write condition properties
    variable dprops; variable ndime
    variable ctbclink

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
	# Solver type
	set cxpath "$AppId//c.AnalysisData//i.SolverType"
	set SolverType [::xmlutils::setXml $cxpath $cproperty]
	# WarnWinText "SolverType:$SolverType"
	switch -exact -- $SolverType {
	    "ElementBased" {
		set ConditionId "WallCondition${ndime}"     
	    }
	    "Monolithic" {
		set ConditionId "MonolithicWallCondition${ndime}"
	    }
	}

	if {$ndime =="2D"} {
	    set cgroupid "-@kratos@b2d"
	    set GiDElemType "Linear"
	    set gprop [dict create]
	    set f "%10i"
	    dict set gprop $cgroupid "$f"
	    if {[write_calc_data has_elements -elemtype $GiDElemType $gprop]} {
		set f "%10d [format "%4d" $fixval] %10d %10d\n"
		set f [subst $f]
		dict set gprop $cgroupid "$f"
		# Write the pressure value
		write_calc_data puts "Begin Conditions $ConditionId"
		# write_calc_data connectivities -elemtype "$GiDElemType" $gprop
		set condid 0
		foreach {elemid cfixval nodei nodej} [write_calc_data connectivities -return -elemtype "$GiDElemType" $gprop] {
		    incr condid 1 
		    # wa "elemid:$elemid cfixval:$cfixval nodei:$nodei nodej:$nodej"
		    write_calc_data puts "[format "%10d %4d %10d %10d" $condid $cfixval $nodei $nodej]"
		    
		    # Update the link between the condition id. and the BC element id
		    dict set ctbclink $elemid $condid
		}
		write_calc_data puts "End Conditions"
		write_calc_data puts ""
	    }
	    unset gprop

	} elseif {$ndime =="3D"} {

	    set cgroupid "-@kratos@b3d"
	    set GiDElemType "Triangle"
	    set gprop [dict create]
	    set f "%10i"
	    dict set gprop $cgroupid "$f"
	    if {[write_calc_data has_elements -elemtype $GiDElemType $gprop]} {
		set f "%10d [format "%4d" $fixval] %10d %10d %10d\n"
		set f [subst $f]
		dict set gprop $cgroupid "$f"
		# Write the condition3D
		write_calc_data puts "Begin Conditions $ConditionId"
		# write_calc_data connectivities -elemtype "$GiDElemType" $gprop
		set condid 0
		foreach {elemid cfixval nodei nodej nodek} [write_calc_data connectivities -return -elemtype "$GiDElemType" $gprop] {
		    incr condid 1 
		    # wa "elemid:$elemid cfixval:$cfixval nodei:$nodei nodej:$nodej"
		    write_calc_data puts "[format "%10d %4d %10d %10d %10d" $condid $cfixval $nodei $nodej $nodek]"

		    # Update the link between the condition id. and the BC element id
		    dict set ctbclink $elemid $condid
		}
		write_calc_data puts "End Conditions"
		write_calc_data puts ""
	    }
	    unset gprop
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

proc ::wkcf::WriteProjectParameters {} {
    # Write the project parameters file
    variable StructuralAnalysis;     variable FluidApplication
    
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
    
    
    # Write structural application parameter file
    if {$StructuralAnalysis =="Yes"} {
	::wkcf::WriteStructuralProjectParameters "StructuralAnalysis" $fileid $PDir
    }    
    
    if {$FluidApplication =="Yes"} {
	::wkcf::WriteFluidProjectParameters "Fluid" $fileid $PDir
    }
    
    close $fileid
    
}

proc ::wkcf::WriteGiDPostMode {AppId fileid} {
    # Write the GiD post mode variables for each applications
    
    # kratos key word xpath
    set kwxpath "Applications/$AppId"
    
    # Gid results
    set gidrlist [list "GiDPostMode" "GiDWriteMeshFlag" "GiDWriteConditionsFlag" "GiDMultiFileFlag"]
    foreach gidr $gidrlist {
	# Get the value
	set cxpath "$AppId//c.Results//c.GiDOptions//i.${gidr}"
	set cproperty "dv"
	set cvalue [::xmlutils::setXml $cxpath $cproperty]
	# Get the kratos keyword
	set gidrkw [::xmlutils::getKKWord $kwxpath $gidr]
	# WarnWinText "gidr:$gidr cvalue:$cvalue gidrkw:$gidrkw"
	if {($gidr=="GiDWriteMeshFlag") || ($gidr=="GiDWriteConditionsFlag")} {
	    if {$cvalue =="Yes"} {
		set cvalue "True"
	    } elseif {$cvalue =="Yes"} {
		set cvalue "False"
	    }
	    puts $fileid "$gidrkw = $cvalue"
	} else {
	    puts $fileid "$gidrkw = \"$cvalue\""
	}
    }
}

















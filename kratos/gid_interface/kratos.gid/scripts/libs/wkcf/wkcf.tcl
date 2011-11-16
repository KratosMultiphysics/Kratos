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
# package require sqlite3

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

    # Debug/Release variable
    variable pflag 

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

	# Write boundary condition block
	::wkcf::WriteBoundaryConditions $AppId

	# For structural analysis application
	if {$AppId =="StructuralAnalysis"} {
	    # Write load properties block
	    ::wkcf::WriteLoads $AppId
	}

	# Write nodal data for density and viscosity
	if {$AppId == "Fluid"} {
	    ::wkcf::WritePropertyAtNodes $AppId
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

    # Write bat file
    # ::wkcf::WriteBatFile

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

proc ::wkcf::WritePropertyAtNodes {AppId} {
    # Write some properties at the nodal level for Fluid application
    variable dprops

    set cproplist [list "Density" "Viscosity"]
    # Check for all defined kratos elements
    if {([info exists dprops($AppId,AllKElemId)]) && ($dprops($AppId,AllKElemId)>0)} {
	# For debug
	if {!$::wkcf::pflag} {
	    set inittime [clock seconds]
	}
	# Create a dictionary
	set nc [dict create 0 0]
	# For all defined kratos elements	
	foreach celemid $dprops($AppId,AllKElemId) {
	    # Check for all defined group identifier for this element
	    if {([info exists dprops($AppId,KElem,$celemid,AllGroupId)]) && ($dprops($AppId,KElem,$celemid,AllGroupId)>0)} {
		# For all defined group identifier for this element
		foreach cgroupid $dprops($AppId,KElem,$celemid,AllGroupId) {
		    # Get the GiD entity type, element type and property identifier
		    lassign $dprops($AppId,KElem,$celemid,$cgroupid,GProps) GiDEntity GiDElemType PropertyId KEKWord nDim
		    # WarnWinText "GiDEntity:$GiDEntity GiDElemType:$GiDElemType PropertyId:$PropertyId KEKWord:$KEKWord nDim:$nDim"
		    # Get all defined entities for this group identifier
 		    set allelist [::KUtils::GetDefinedMeshGiDEntities $cgroupid $GiDEntity]
 		    # WarnWinText "alllist:$allelist"
 		    if {[llength $allelist]>0} {
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

	# Write viscosity and density for each node identifier
	set Density 0.0; set Viscosity 0.0
	foreach PropertyId $dprops($AppId,GKProps,AllPropertyId) {
	    # Get the material identifier for this property 
	    set MatId $dprops($AppId,Property,$PropertyId,MatId) 
	    # Get the group identifier
	    set GroupId $dprops($AppId,Property,$PropertyId,GroupId)
	    # Get all material properties
	    set mpxpath "[::KMat::findMaterialParent $MatId]//m.${MatId}"
	    # WarnWinText "mpxpath:$mpxpath"
	    # Get the material properties
	    foreach pid $cproplist {
		if {$pid =="Density"} {
		    set xpath "c.General"
		    # Get the current value for this properties
		    set cvalue [lindex [::KMat::getMaterialProperties "p" "$mpxpath//$xpath//p.$pid"] 0 1]
		    set Density [GiD_FormatReal "%10.5e" $cvalue]
		} elseif {$pid =="Viscosity"} {
		    set xpath "c.Fluid"
		    # Get the current value for this properties
		    set cvalue [lindex [::KMat::getMaterialProperties "p" "$mpxpath//$xpath//p.$pid"] 0 1]
		    set Viscosity [GiD_FormatReal "%10.5e" $cvalue]
		}
	    }
	    # Only the first property
	    break 
	}
	# WarnWinText "Density:$Density Viscosity:$Viscosity"
	    
	set kxpath "Materials"
	set cpropid "0"	

	# Write the group nodal properties
	foreach celemid $dprops($AppId,AllKElemId) {
	    # Check for all defined group identifier for this element
	    if {([info exists dprops($AppId,KElem,$celemid,AllGroupId)]) && ($dprops($AppId,KElem,$celemid,AllGroupId)>0)} {
	 	# For all defined group identifier for this element
	 	foreach cgroupid $dprops($AppId,KElem,$celemid,AllGroupId) {
		    set cnodeglist [list]
		    dict for {nodeid dgroupid} $nc {
			if {$cgroupid ==$dgroupid} {
			    lappend cnodeglist $nodeid
			}
		    }
		    if {[llength $cnodeglist]} {
		        # Write all nodes for this group in incresing orden
			set viscobf ""; set densibf ""
		        foreach nodeid [lsort -integer $cnodeglist] {
			    append viscobf "[format "%4i  %4i" $nodeid $cpropid]    $Viscosity\n"
			    append densibf "[format "%4i  %4i" $nodeid $cpropid]    $Density\n"
			}
			# Write viscosity value for this group
			set vkword [::xmlutils::getKKWord $kxpath "Viscosity" "kkword"]
			write_calc_data puts "Begin NodalData $vkword \/\/ GUI group identifier: $cgroupid"
			write_calc_data puts "[string trimright $viscobf]"
			write_calc_data puts "End NodalData"
			write_calc_data puts ""
			
			# Write density value for this group 
			set dkword [::xmlutils::getKKWord $kxpath "Density" "kkword"]
			write_calc_data puts "Begin NodalData $dkword \/\/ GUI group identifier: $cgroupid"
			write_calc_data puts "[string trimright $densibf]"
			write_calc_data puts "End NodalData"
			write_calc_data puts ""
		    }
		}
		unset cnodeglist
	    }
	}
 	unset nc

	# For debug
	if {!$::wkcf::pflag} {
	    set endtime [clock seconds]
	    set ttime [expr $endtime-$inittime]
	    WarnWinText "endtime:$endtime ttime:$ttime"
	    WarnWinText "Property at nodes [::KUtils::Duration $ttime]"
	}
    }
}

proc ::wkcf::WriteNodalCoordinates {AppId} {
    # Write the nodal coordinates block
    # Nodes block format
    # Begin Nodes
    # // id	  X	Y	Z
    # End Nodes
    # Arguments
    # AppId  => Application identifier
    variable ndime; variable dprops
   
     # Check for all defined kratos elements
    if {([info exists dprops($AppId,AllKElemId)]) && ($dprops($AppId,AllKElemId)>0)} {
	# For debug
	if {!$::wkcf::pflag} {
	    set inittime [clock seconds]
	}
	set cnodeglist [list]
	# Create a dictionary
	set nc [dict create 0 0]
	# For all defined kratos nodes	
	foreach celemid $dprops($AppId,AllKElemId) {
	    # Check for all defined group identifier for this element
	    if {([info exists dprops($AppId,KElem,$celemid,AllGroupId)]) && ($dprops($AppId,KElem,$celemid,AllGroupId)>0)} {
		# For all defined group identifier for this element
		foreach cgroupid $dprops($AppId,KElem,$celemid,AllGroupId) {
		    # Get the GiD entity type, element type and property identifier
		    lassign $dprops($AppId,KElem,$celemid,$cgroupid,GProps) GiDEntity GiDElemType PropertyId KEKWord nDim
		    # WarnWinText "GiDEntity:$GiDEntity GiDElemType:$GiDElemType PropertyId:$PropertyId KEKWord:$KEKWord nDim:$nDim"
		    # Get all defined entities for this group identifier
 		    set allelist [::KUtils::GetDefinedMeshGiDEntities $cgroupid $GiDEntity]
		    # WarnWinText "alllist:$allelist"
 		    if {[llength $allelist]>0} {
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

	# For debug
	if {!$::wkcf::pflag} {
	    set endtime [clock seconds]
	    set ttime [expr $endtime-$inittime]
	    WarnWinText "endtime:$endtime ttime:$ttime"
	    WarnWinText "Coordinates [::KUtils::Duration $ttime]"
	}
    }
}

proc ::wkcf::WriteElementConnectivities {AppId} {
    # Write the element connectivities block
    # Element block format
    # Begin Elements element_name
    # // id prop_id	 n1	n2	n3	...
    # End Elements
    # Arguments
    # AppId => Application identifier
    variable useqelem; variable dprops
    global KPriv
    
    set shellidlist [list "ShellIsotropic" "ShellAnisotropic" "Membrane" "EBST" "ASGS2D" "ASGS3D" "Fluid2D" "Fluid3D"]
    # Check for all defined kratos elements
    if {([info exists dprops($AppId,AllKElemId)]) && ($dprops($AppId,AllKElemId)>0)} {
	# For debug
	if {!$::wkcf::pflag} {
	    set inittime [clock seconds]
	}

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

		    # Get the element identifier => nDimNNode
		    set etbf ""
		    if {$celemid ni $shellidlist} {
			set etbf [::wkcf::GetnDimnNode $GiDElemType $nDim]
		    } 
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
			    if {[llength $eprop]>0} {
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
		}
	    }
	}
	write_calc_data puts ""
	
	# For debug
	if {!$::wkcf::pflag} {
	    set endtime [clock seconds]
	    set ttime [expr $endtime-$inittime]
	    WarnWinText "endtime:$endtime ttime:$ttime"
	    WarnWinText "Connectivities [::KUtils::Duration $ttime]"
	}
    }
}

proc ::wkcf::WriteBoundaryConditions {AppId} {
    # Write the boundary condition block
    variable dprops
  
    # Check for all defined condition type
    if {([info exists dprops($AppId,AllBCTypeId)]) && ([llength $dprops($AppId,AllBCTypeId)]>0)} {
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
	    if {([info exists dprops($AppId,BC,$ccondid,AllGroupId)]) && ($dprops($AppId,BC,$ccondid,AllGroupId)>0)} {
		# Select the condition type
		switch -exact -- $ccondid {
		    "Displacements"
		    {
			set kwxpath "Applications/StructuralAnalysis"
			set kwordlist [list [::xmlutils::getKKWord $kwxpath "Dx" "kkword"] [::xmlutils::getKKWord $kwxpath "Dy" "kkword"] [::xmlutils::getKKWord $kwxpath "Dz" "kkword"]]
			# Process displacement properties
			::wkcf::WriteDispRotBC $AppId $ccondid $kwordlist 
		    }
		    "Rotations"
		    {
			set kwxpath "Applications/StructuralAnalysis"
			set kwordlist [list [::xmlutils::getKKWord $kwxpath "Rx"] [::xmlutils::getKKWord $kwxpath "Ry"] [::xmlutils::getKKWord $kwxpath "Rz"]]
			# Process displacement properties
			::wkcf::WriteDispRotBC $AppId $ccondid $kwordlist 
		    }
		    "OutletPressure"
		    {
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
	    WarnWinText "endtime:$endtime ttime:$ttime"
	    WarnWinText "BC [::KUtils::Duration $ttime]"
	}
    }
}

proc ::wkcf::WriteFluidBC {AppId inletvelglist noslipglist flagvariablelist kwordlist} {
    variable gidentitylist; variable ndime
    variable useqelem; variable dprops

    # WarnWinText "inletvelglist:$inletvelglist\nnoslipglist:$noslipglist\nflagvariablelist:$flagvariablelist\nkwordlist:$kwordlist"
    # Map Inlet-NoSlip => Use no-slip values at share nodes
    set icondid "InletVelocity"; set nscondid "No-Slip"
    set cpropid "1"
    set xitem [lindex $kwordlist 0]
    set yitem [lindex $kwordlist 1]
    set zitem [lindex $kwordlist 2]
   
    if {[llength $noslipglist]} {
	# No-slip
	foreach cgroupid $noslipglist {
	    set allnslip($cgroupid,NodeList) [list]
	    foreach GiDEntity $gidentitylist {
		# Get all defined entities for this group identifier
		switch $GiDEntity {
		    "point" {
			set callnlist [::KUtils::GetDefinedMeshGiDEntities $cgroupid $GiDEntity "Nodes" $useqelem]
			if {[llength $callnlist]} {
			    lappend allnslip($cgroupid,NodeList) {*}$callnlist
			}
		    }
		    "line" - "surface" {
			set callnlist [::KUtils::GetDefinedMeshGiDEntities $cgroupid $GiDEntity "Nodes" $useqelem]
			if {[llength $callnlist]} {
			    lappend allnslip($cgroupid,NodeList) {*}$callnlist
			}
		    }
		    "volume" {
			if {$ndime =="3D"} {
			    set callnlist [::KUtils::GetDefinedMeshGiDEntities $cgroupid $GiDEntity "Nodes" $useqelem]
			    if {[llength $callnlist]} {
				lappend allnslip($cgroupid,NodeList) {*}$callnlist
			    }
			}
		    }
		} 
	    }
	    # WarnWinText "groupid:$cgroupid NodeList:$allnslip($cgroupid,NodeList)"
	}

	# Check to match node identifier
	# Use first the inlet
	set xcomp ""; set ycomp ""; set zcomp ""
	if {[llength $inletvelglist]} {
	    set condmatch [dict create none 0]
	    # For each group in the no-slip condition
	    foreach nsgroupid $noslipglist {
		set nsGProps $dprops($AppId,BC,$nscondid,$nsgroupid,GProps)
		# WarnWinText "nsgroupid:$nsgroupid nsGProps:$nsGProps"
		foreach nsnodeid $allnslip($nsgroupid,NodeList) {
		    set clist [list 0 0 0]
		    if {[lindex $nsGProps 0]} {
			append xcomp "[format "%8i%8i" $nsnodeid $cpropid]   [GiD_FormatReal "%10.5e" [lindex $nsGProps 1]]\n"
			lset clist 0 1 
		    }
		    if {[lindex $nsGProps 2]} {
			append ycomp "[format "%8i%8i" $nsnodeid $cpropid]   [GiD_FormatReal "%10.5e" [lindex $nsGProps 3]]\n"
			lset clist 1 1
		    }
		    if {[lindex $nsGProps 4]} {
			append zcomp "[format "%8i%8i" $nsnodeid $cpropid]   [GiD_FormatReal "%10.5e" [lindex $nsGProps 5]]\n"
			lset clist 2 1
		    }
		    dict set condmatch $nsnodeid $clist 
		}
		
		# Write this group identifier
		if {[string length $xcomp]} {
		    write_calc_data puts "Begin NodalData $xitem \/\/ No-slip condition GUI group identifier: $nsgroupid"
		    write_calc_data puts "[string trimright ${xcomp}]"
		    write_calc_data puts "End NodalData"
		    write_calc_data puts ""
		}
		if {[string length $ycomp]} {
		    write_calc_data puts "Begin NodalData $yitem \/\/ No-slip condition GUI group identifier: $nsgroupid"
		    write_calc_data puts "[string trimright ${ycomp}]"
		    write_calc_data puts "End NodalData"
		    write_calc_data puts ""
		}
		if {$ndime =="3D"} {
		    if {[string length $zcomp]} {
			write_calc_data puts "Begin NodalData $zitem \/\/ No-slip condition GUI group identifier: $nsgroupid"
			write_calc_data puts "[string trimright ${zcomp}]"
			write_calc_data puts "End NodalData"
			write_calc_data puts ""
		    }
		}

		# Reset xcomp, ycomp and zcomp
		set xcomp ""; set ycomp ""; set zcomp ""
	    }

	    # Get the Inlet velocity entities
	    foreach cgroupid $inletvelglist {
		set allninlet($cgroupid,NodeList) [list]
		foreach GiDEntity $gidentitylist {
		    # Get all defined entities for this group identifier
		    switch $GiDEntity {
			"point" {
			    set callnlist [::KUtils::GetDefinedMeshGiDEntities $cgroupid $GiDEntity "Nodes" $useqelem]
			    if {[llength $callnlist]} {
				lappend allninlet($cgroupid,NodeList) {*}$callnlist
			    }
			}
			"line" - "surface" {
			    set callnlist [::KUtils::GetDefinedMeshGiDEntities $cgroupid $GiDEntity "Nodes" $useqelem]
			    if {[llength $callnlist]} {
				lappend allninlet($cgroupid,NodeList) {*}$callnlist
			    }
			}
			"volume" {
			    if {$ndime =="3D"} {
				set callnlist [::KUtils::GetDefinedMeshGiDEntities $cgroupid $GiDEntity "Nodes" $useqelem]
				if {[llength $callnlist]} {
				    lappend allninlet($cgroupid,NodeList) {*}$callnlist
				}
			    }
			}
		    } 
		}
		# WarnWinText "groupid:$cgroupid NodeList:$allninlet($cgroupid,NodeList)"
	    }
  
	    # For all inlet velocity group identifier
	    set ixcomp ""; set iycomp ""; set izcomp ""
	    foreach igroupid $inletvelglist {
		set iGProps $dprops($AppId,BC,$icondid,$igroupid,GProps)
		# WarnWinText "igroupid:$igroupid iGProps:$iGProps"
		foreach inodeid $allninlet($igroupid,NodeList) {
		    # WarnWinText "inodeid:$inodeid"
		    # Check that this node identifier exists in the dictionary
		    if {[dict exists $condmatch $inodeid]} {
			# Get the properties
			set nprop [dict get $condmatch $inodeid] 
			# WarnWinText "nprop:$nprop"
			# Check x flag
			if {[lindex $nprop 0]=="0"} {
			    # Write this node identifier
			    if {[lindex $iGProps 0]} {
				append ixcomp "[format "%8i%8i" $inodeid $cpropid]   [GiD_FormatReal "%10.5e" [lindex $iGProps 1]]\n"
			    }
			}
			# Check y flag
			if {[lindex $nprop 1]=="0"} {
			    if {[lindex $iGProps 2]} {
				append iycomp "[format "%8i%8i" $inodeid $cpropid]   [GiD_FormatReal "%10.5e" [lindex $iGProps 3]]\n"
			    }
			}
			# Check z flag
			if {[lindex $nprop 2]=="0"} {
			    if {[lindex $iGProps 4]} {
				append izcomp "[format "%8i%8i" $inodeid $cpropid]   [GiD_FormatReal "%10.5e" [lindex $iGProps 5]]\n"
			    }
			}
		    } else {
			# Write this node identifier
			if {[lindex $iGProps 0]} {
			    append ixcomp "[format "%8i%8i" $inodeid $cpropid]   [GiD_FormatReal "%10.5e" [lindex $iGProps 1]]\n"
			}
			if {[lindex $iGProps 2]} {
			    append iycomp "[format "%8i%8i" $inodeid $cpropid]   [GiD_FormatReal "%10.5e" [lindex $iGProps 3]]\n"
			}
			if {[lindex $iGProps 4]} {
			    append izcomp "[format "%8i%8i" $inodeid $cpropid]   [GiD_FormatReal "%10.5e" [lindex $iGProps 5]]\n"
			}
		    }
		}
		
		# Write this group identifier
		if {[string length $ixcomp]} {
		    write_calc_data puts "Begin NodalData $xitem \/\/ Inlet velocity condition GUI group identifier: $igroupid"
		    write_calc_data puts "[string trimright ${ixcomp}]"
		    write_calc_data puts "End NodalData"
		    write_calc_data puts ""
		}
		if {[string length $iycomp]} {
		    write_calc_data puts "Begin NodalData $yitem \/\/ Inlet velocity condition GUI group identifier: $igroupid"
		    write_calc_data puts "[string trimright ${iycomp}]"
		    write_calc_data puts "End NodalData"
		    write_calc_data puts ""
		}
		if {$ndime =="3D"} {
		    if {[string length $izcomp]} {
			write_calc_data puts "Begin NodalData $zitem \/\/ Inlet velocity condition GUI group identifier: $igroupid"
			write_calc_data puts "[string trimright ${izcomp}]"
			write_calc_data puts "End NodalData"
			write_calc_data puts ""
		    }
		}
		
		# Reset ixcomp, iycomp and zcomp
		set ixcomp ""; set iycomp ""; set izcomp ""
	    }
	    
	    # unset dictionary variable
	    unset condmatch
	    if {[info exists allninlet]} {
		unset allninlet
	    }
	    if {[info exists allnslip]} {
		unset allnslip
	    }
	
	} else {

	    # Write all no-slip condition
	    # For each group in the no-slip condition
	    foreach nsgroupid $noslipglist {
		set nsGProps $dprops($AppId,BC,$nscondid,$nsgroupid,GProps)
		foreach nsnodeid $allnslip($nsgroupid,NodeList) {
		    if {[lindex $nsGProps 0]} {
			append xcomp "[format "%8i%8i" $nsnodeid $cpropid]   [GiD_FormatReal "%10.5e" [lindex $nsGProps 1]]\n"
		    }
		    if {[lindex $nsGProps 2]} {
			append ycomp "[format "%8i%8i" $nsnodeid $cpropid]   [GiD_FormatReal "%10.5e" [lindex $nsGProps 3]]\n"
		    }
		    if {[lindex $nsGProps 4]} {
			append zcomp "[format "%8i%8i" $nsnodeid $cpropid]   [GiD_FormatReal "%10.5e" [lindex $nsGProps 5]]\n"
		    }
		}
		
		# Write this group identifier
		if {[string length $xcomp]} {
		    write_calc_data puts "Begin NodalData $xitem \/\/ No-slip condition GUI group identifier: $nsgroupid"
		    write_calc_data puts "[string trimright ${xcomp}]"
		    write_calc_data puts "End NodalData"
		    write_calc_data puts ""
		}
		if {[string length $ycomp]} {
		    write_calc_data puts "Begin NodalData $yitem \/\/ No-slip condition GUI group identifier: $nsgroupid"
		    write_calc_data puts "[string trimright ${ycomp}]"
		    write_calc_data puts "End NodalData"
		    write_calc_data puts ""
		}
		if {$ndime =="3D"} {
		    if {[string length $zcomp]} {
			write_calc_data puts "Begin NodalData $zitem \/\/ No-slip condition GUI group identifier: $nsgroupid"
			write_calc_data puts "[string trimright ${zcomp}]"
			write_calc_data puts "End NodalData"
			write_calc_data puts ""
		    }
		}
		
		# Reset xcomp, ycomp and zcomp
		set xcomp ""; set ycomp ""; set zcomp ""
	    }
	}
    }
    
    # Write Flag-variable and is_boundary nodal data conditions
    ::wkcf::WriteFluidFlagVariableBC $AppId $flagvariablelist

}

proc ::wkcf::WriteFluidFlagVariableBC {AppId flagvariablelist} {
    variable gidentitylist; variable ndime
    variable useqelem; variable dprops

    WarnWinText "flagvariablelist:$flagvariablelist"
    # For nodes with many flag variable defined flag of level two have the priority over flag of level one
    set flagvarcondid "Flag-Variable"
    set cpropid "0"
    set isbcpropid "1"

    if {[llength $flagvariablelist]} {
	# Flag-Variable
	foreach cgroupid $flagvariablelist {
	    set allnflagvar($cgroupid,NodeList) [list]
	    foreach GiDEntity $gidentitylist {
		# Get all defined entities for this group identifier
		switch $GiDEntity {
		    "point" {
			set callnlist [::KUtils::GetDefinedMeshGiDEntities $cgroupid $GiDEntity "Nodes" $useqelem]
			if {[llength $callnlist]} {
			    lappend allnflagvar($cgroupid,NodeList) {*}$callnlist
			}
		    }
		    "line" - "surface" {
			set callnlist [::KUtils::GetDefinedMeshGiDEntities $cgroupid $GiDEntity "Nodes" $useqelem]
			if {[llength $callnlist]} {
			    lappend allnflagvar($cgroupid,NodeList) {*}$callnlist
			}
		    }
		    "volume" {
			if {$ndime =="3D"} {
			    set callnlist [::KUtils::GetDefinedMeshGiDEntities $cgroupid $GiDEntity "Nodes" $useqelem]
			    if {[llength $callnlist]} {
				lappend allnflagvar($cgroupid,NodeList) {*}$callnlist
			    }
			}
		    }
		} 
	    }
	    # WarnWinText "groupid:$cgroupid NodeList:$allnflagvar($cgroupid,NodeList)"
	}

	# Write the flag condition
	set flag1 0
	set fvcomp ""; set isbcomp ""
	set fvitem "FLAG-VARIABLE"
	set isbitem "IS_BOUNDARY"
	set flagvar2 [dict create]
	# For each group in the flag-variable condition
	# Create a dict for all nodes with flag equal to 2
	foreach groupid $flagvariablelist {
	    set GProps $dprops($AppId,BC,$flagvarcondid,$groupid,GProps)
	    # WarnWinText "groupid:$groupid GProps:$GProps"
	    if {[lindex $GProps 0]=="2"} {
		foreach nodeid $allnflagvar($groupid,NodeList) {
		    append fvcomp "[format "%8i%8i%8i" $nodeid $cpropid [lindex $GProps 0]]\n"
		    append isbcomp "[format "%8i%8i%8i" $nodeid $cpropid $isbcpropid]\n"
		    dict set flagvar2 $nodeid $groupid
		    set flag1 1
		}
		# Write this group identifier
		# Flag-Variable
		if {[string length $fvcomp]} {
		    write_calc_data puts "Begin NodalData $fvitem \/\/ Flag-Variable condition GUI group identifier: $groupid"
		    write_calc_data puts "[string trimright ${fvcomp}]"
		    write_calc_data puts "End NodalData"
		    write_calc_data puts ""
		}
		# is_boundary
		if {[string length $isbcomp]} {
		    write_calc_data puts "Begin NodalData $isbitem \/\/ is_boundary associated with Flag-Variable condition GUI group identifier: $groupid"
		    write_calc_data puts "[string trimright ${isbcomp}]"
		    write_calc_data puts "End NodalData"
		    write_calc_data puts ""
		}
		# Reset components
		set fvcomp "";	set isbcomp ""
	    }
	}

	# Write all group with flag-variable equal to 1
	if {$flag1} {
	    set fvcomp ""; set isbcomp ""
	    # For each group in the flag-variable condition
	    foreach groupid $flagvariablelist {
		set GProps $dprops($AppId,BC,$flagvarcondid,$groupid,GProps)
		# WarnWinText "groupid:$groupid GProps:$GProps"
		if {[lindex $GProps 0]=="1"} {
		    foreach nodeid $allnflagvar($groupid,NodeList) {
			if {![dict exists $flagvar2 $nodeid]} {
			    append fvcomp "[format "%8i%8i%8i" $nodeid $cpropid [lindex $GProps 0]]\n"
			    append isbcomp "[format "%8i%8i%8i" $nodeid $cpropid $isbcpropid]\n"
			}
		    }
		    # Write this group identifier
		    # Flag-Variable
		    if {[string length $fvcomp]} {
			write_calc_data puts "Begin NodalData $fvitem \/\/ Flag-Variable condition GUI group identifier: $groupid"
			write_calc_data puts "[string trimright ${fvcomp}]"
			write_calc_data puts "End NodalData"
			write_calc_data puts ""
		    }
		    # is_boundary
		    if {[string length $isbcomp]} {
			write_calc_data puts "Begin NodalData $isbitem \/\/ is_boundary associated with Flag-Variable condition GUI group identifier: $groupid"
			write_calc_data puts "[string trimright ${isbcomp}]"
			write_calc_data puts "End NodalData"
			write_calc_data puts ""
		    }
		    # Reset components
		    set fvcomp ""; set isbcomp ""
		}
	    }
	}
    }
}

proc ::wkcf::WriteConditions {} {
    # ABSTRACT: Write condition properties
    variable dprops; variable ndime

    # Check for all defined kratos elements
    if {([info exists dprops(AllKElemId)]) && ($dprops(AllKElemId)>0)} {

	# Find boundaries
	if {$ndime =="2D"} {
	    set blinelist [::wkcf::FindBoundaries line]
	    # WarnWinText "belist:$blinelist"

	} elseif {$ndime =="3D"} {

	    set bsurfacelist [::wkcf::FindBoundaries surface]
	    # WarnWinText "bsurfacelist:$bsurfacelist"
	}

	# For all defined kratos elements	
	foreach celemid $dprops(AllKElemId) {
	    # Check for all defined group identifier for this element
	    if {([info exists dprops(KElem,$celemid,AllGroupId)]) && ($dprops(KElem,$celemid,AllGroupId)>0)} {
		# For all defined group identifier for this element
		foreach cgroupid $dprops(KElem,$celemid,AllGroupId) {
		    # Get the GiD entity type, element type and property identifier
		    lassign $dprops(KElem,$celemid,$cgroupid,GProps) GiDEntity GiDElemType PropertyId KEKWord nDim
		    # WarnWinText "GiDEntity:$GiDEntity GiDElemType:$GiDElemType PropertyId:$PropertyId KEKWord:$KEKWord nDim:$nDim"
		    
		    # Get all defined entities for this group identifier
		    set allelist [::KUtils::GetDefinedMeshGiDEntities $cgroupid $GiDEntity]
		    # WarnWinText "alllist:$allelist"
		    if {[llength $allelist]>0} {
		    }
		}
	    }
	}
    }
}

proc ::wkcf::WriteOutLetPressureBC {AppId ccondid kwordlist} {
    # Write outlet pressure boundary condition
    variable ndime;    variable gidentitylist
    variable useqelem; variable dprops

    set nodelist [list]
     
    # For all defined group identifier inside this condition type
    foreach cgroupid $dprops($AppId,BC,$ccondid,AllGroupId) {
	# Get the condition properties
	set GProps $dprops($AppId,BC,$ccondid,$cgroupid,GProps)
	# WarnWinText "GProps:$GProps"
	# Assign values
	lassign $GProps fixval pressureval
	# WarnWinText "fixval:$fixval pressureval:$pressureval"
	set allnlist [list]
	foreach GiDEntity $gidentitylist {
	    # Get all defined entities for this group identifier
	    switch $GiDEntity {
	 	"point" {
	 	    set callnlist [::KUtils::GetDefinedMeshGiDEntities $cgroupid $GiDEntity "Nodes" $useqelem]
	 	    if {[llength $callnlist]} {
	 		lappend allnlist $callnlist
	 	    }
		}
	 	"line" - "surface" {
	 	    set callnlist [::KUtils::GetDefinedMeshGiDEntities $cgroupid $GiDEntity "Nodes" $useqelem]
	 	    if {[llength $callnlist]} {
	 		lappend allnlist $callnlist
	 	    }
		} 
	    }
	}
	# WarnWinText "$GiDEntity alllist:$allnlist"
	foreach cprop $allnlist {
	    set cprop [lsort -integer -unique $cprop]
	    foreach nodeid $cprop {
	 	# Fix x
	 	if {$fixval =="1"} {
	 	    lappend nodelist "$nodeid 1 $pressureval"
	 	}
	    }
	}
    }
    # WarnWinText "nodelist:$nodelist"
       
    # PRESSURE
    if {[llength $nodelist]>0} {
    	set kitem [lindex $kwordlist 0]
    	write_calc_data puts "Begin NodalData $kitem"
    	foreach citem $nodelist {
    	    lassign $citem nodeid fix pval
    	    set cf "[format "%4i%4i%10.5f" $nodeid $fix $pval]"
    	    write_calc_data puts "$cf"
    	}
    	write_calc_data puts "End NodalData"
    	write_calc_data puts ""
	unset nodelist
    }
}

proc ::wkcf::WriteDispRotBC {AppId ccondid kwordlist} {
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
	    if {([info exists dprops($AppId,Loads,$cloadtid,AllGroupId)]) && ($dprops($AppId,Loads,$cloadtid,AllGroupId)>0)} {
		# Select the load type
		switch -exact -- $cloadtid {
		    "Puntual"
		    {
			# Concentrate or puntual loads (forces or moments)
			set kwordlist [list [::xmlutils::getKKWord $kwxpath "Fx"] [::xmlutils::getKKWord $kwxpath "Fy"] [::xmlutils::getKKWord $kwxpath "Fz"]]
			# Write puntual loads
			::wkcf::WritePuntualLoads $AppId $cloadtid $kwordlist
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

proc ::wkcf::WritePressureLoads {AppId cloadtid} {
    # Write pressure loads (positive or negative shell pressure)
    variable ndime;    variable gidentitylist
    variable useqelem; variable dprops

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
	    set icondid 0
	    # Get all defined nodes
	    foreach elemid $celemglist {
		incr icondid 1
		# Get the element properties
		lassign [lrange [GiD_Info Mesh Elements $GiDElemType $elemid] 1 end-1] N1 N2 N3
		set cf "[format "%4i%4i%8i%8i%8i" $icondid $RefPropId $N1 $N2 $N3]"
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
    # Write concentrated loads (puntual force or moment)
    variable ndime;    variable gidentitylist
    variable useqelem; variable dprops

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
	set allnlist [lsort -integer -unique {*}$allnlist]
	# WarnWinText "$GiDEntity alllist:$allnlist"

	set usepointforce "No"
	switch -exact -- $ndime {
	    "2D"
	    {
		set pointforcekword "PointForce2D"
		foreach item [list $Fx $Fy $Mx $My] {
		    if {$item !="0.0"} {
			set usepointforce "Yes"
			break
		    }
		}
	    }
	    "3D"
	    {
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
	    set jj 0
	    foreach nodeid $allnlist {
		incr jj
		set cf "[format "%4i%4i%8i" $jj 1 $nodeid]"
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
}

proc ::wkcf::WriteBodyForceValues {props} {
    # Write the gravity properties to the kratos data file
    # Arguments
    # props => Body force properties
    variable ndime

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
    write_calc_data puts " BODY_FORCE \[3\] $vector"
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

proc ::wkcf::WriteFluidProjectParameters {AppId fileid PDir} {
    variable ndime

    set rootid "$AppId"
    # Kratos key word xpath
    set kxpath "Applications/$rootid"

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
    set cproperty "dv"
    
    # Fluid type
    set cxpath "$rootid//c.AnalysisData//i.FluidType"
    set FluidType [::xmlutils::setXml $cxpath $cproperty]
    # WarnWinText "FluidType:$FluidType"
    if {$FluidType =="Compressible"} {

    } elseif {$FluidType =="Incompressible"} {

	# Free surface
	set cxpath "$rootid//c.AnalysisData//i.FreeSurface"
	set FreeSurface [::xmlutils::setXml $cxpath $cproperty]
	# WarnWinText "FreeSurface:$FreeSurface"
	if {$FreeSurface =="No"} {
	    # Solver type
	    set cxpath "$rootid//c.AnalysisData//i.SolverType"
	    set SolverType [::xmlutils::setXml $cxpath $cproperty]
	    # WarnWinText "SolverType:$SolverType"
	    # Get the kratos keyword
	    set ckword [::xmlutils::getKKWord $kxpath $SolverType]
	    # WarnWinText "ckword:$ckword"
	    puts $fileid "SolverType = \"$ckword\""
	   
	    # Monolithic,PressureSplitting,ElementBased,EdgeBased
	    puts $fileid ""
	    switch -exact -- $SolverType {
		"ElementBased" {
		    # Fractional step options => ElementBased

		    # Solution strategy
		    # Linear solvers
		    # Velocity
		    ::wkcf::WriteFluidSolvers $rootid $fileid "Velocity"
		    puts $fileid ""
		    # Pressure
		    ::wkcf::WriteFluidSolvers $rootid $fileid "Pressure"
		 
		}
		"PressureSplitting" {
		    # Pressure splitting
		    
		    # Solution strategy
		    # Linear solvers
		    # Velocity
		    ::wkcf::WriteFluidSolvers $rootid $fileid "Velocity"
		    puts $fileid ""
		    # Pressure
		    ::wkcf::WriteFluidSolvers $rootid $fileid "Pressure"

		}
		"Monolithic" {
		    # Monolithic

		    # Solution strategy
		    # Linear solvers
		    # Velocity
		    ::wkcf::WriteFluidSolvers $rootid $fileid "Monolithic"
		    puts $fileid ""

		}
	    }
	  
	    # Write relative and absolute tolerances
	    puts $fileid ""
	    set ctlist [list "RelativeVelocityTolerance" "AbsoluteVelocityTolerance" "RelativePressureTolerance" "AbsolutePressureTolerance"]
	    foreach cv $ctlist {
		set cxpath "$rootid//c.SolutionStrategy//c.Advanced//i.${cv}"
		set cvalue [::xmlutils::setXml $cxpath $cproperty]
		set ckword [::xmlutils::getKKWord $kxpath $cv]
		puts $fileid "$ckword = $cvalue"
	    }
	  	  
	    puts $fileid ""  
	    # Time order
	    set cxpath "$rootid//c.SolutionStrategy//c.Advanced//i.TimeOrder"
	    set TimeOrder [::xmlutils::setXml $cxpath $cproperty]
	    puts $fileid "time_order = $TimeOrder"

	    # Predictor corrector
	    set cxpath "$rootid//c.SolutionStrategy//c.Advanced//i.PredictorCorrector"
	    set PredictorCorrector [::xmlutils::setXml $cxpath $cproperty]
	    puts $fileid "predictor_corrector = $PredictorCorrector"
	    
	    if {$SolverType in [list "ElementBased" "EdgeBased"]} {
		# Maximum velocity iterations
		set cxpath "$rootid//c.SolutionStrategy//c.Advanced//i.MaximumVelocityIterations"
		set MaximumVelocityIterations [::xmlutils::setXml $cxpath $cproperty]
		puts $fileid "max_vel_its = $MaximumVelocityIterations"
		
		# Maximum pressure iterations
		set cxpath "$rootid//c.SolutionStrategy//c.Advanced//i.MaximumPressureIterations"
		set MaximumPressureIterations [::xmlutils::setXml $cxpath $cproperty]
		puts $fileid "max_press_its = $MaximumPressureIterations"

	    } elseif {$SolverType in [list "Monolithic" "PressureSplitting"]} {
		# Maximum iterations
		set cxpath "$rootid//c.SolutionStrategy//c.Advanced//i.MaximumIterations"
		set MaximumIterations [::xmlutils::setXml $cxpath $cproperty]
		puts $fileid "max_iterations = $MaximumIterations"

	    }

	    # Laplacian form
	    set cxpath "$rootid//c.SolutionStrategy//c.Advanced//i.LaplacianForm"
	    set LaplacianForm [::xmlutils::setXml $cxpath $cproperty]
	    # Get the kratos keyword
	    set ckword [::xmlutils::getKKWord $kxpath $LaplacianForm]
	    puts $fileid "laplacian_form = $ckword"
	    
	} else {
	    # Solver type for free surface
	    set cxpath "$rootid//c.AnalysisData//i.SolverTypeFreeSurf"
	    set SolverTypeFreeSurf [::xmlutils::setXml $cxpath $cproperty]
	    # WarnWinText "SolverTypeFreeSurf:$SolverTypeFreeSurf"
	}
    }

    puts $fileid ""
    # Start time
    set cxpath "$rootid//c.SolutionStrategy//i.StartTime"
    set StartTime [::xmlutils::setXml $cxpath $cproperty]
    # End time
    set cxpath "$rootid//c.SolutionStrategy//i.EndTime"
    set EndTime [::xmlutils::setXml $cxpath $cproperty]
    # Delta time
    set cxpath "$rootid//c.SolutionStrategy//i.DeltaTime"
    set DeltaTime [::xmlutils::setXml $cxpath $cproperty]
    
    # WarnWinText "StartTime:$StartTime EndTime:$EndTime DeltaTime:$DeltaTime"
    puts $fileid "Dt = $DeltaTime"
    puts $fileid "Start_time = $StartTime"
    puts $fileid "max_time = $EndTime"
    set nsteps [expr int(double($EndTime-$StartTime)/double($DeltaTime))]
    puts $fileid "nsteps = $nsteps"

    puts $fileid ""
    # Use dt in stabilization => DynamicTau
    set cxpath "$rootid//c.SolutionStrategy//c.Advanced//i.DynamicTau"
    set DynamicTau [::xmlutils::setXml $cxpath $cproperty]
    puts $fileid "use_dt_in_stabilization = $DynamicTau"
    
    # Use ortogonal subscales => OssSwitch
    set cxpath "$rootid//c.SolutionStrategy//c.Advanced//i.OssSwitch"
    set OssSwitch [::xmlutils::setXml $cxpath $cproperty]
    puts $fileid "use_orthogonal_subscales = $OssSwitch"
 
    # Calculate reactions
    set cxpath "$rootid//c.Results//c.OnNodes//i.Reactions"
    set Reactions [::xmlutils::setXml $cxpath $cproperty]
    if {$Reactions =="Yes"} {
	puts $fileid "Calculate_reactions = True"
    } else {
	puts $fileid "Calculate_reactions = False"
    }

    puts $fileid ""
    # Output step 
    set cxpath "$rootid//c.Results//i.OutputDeltaTime"
    set OutputDeltaTime [::xmlutils::setXml $cxpath $cproperty]
    puts $fileid "output_time = $OutputDeltaTime"
    # WarnWinText "OutputDeltaTime:$OutputDeltaTime"
    set output_step [expr int($OutputDeltaTime/double($DeltaTime))]
    # WarnWinText "output_step:$output_step"
    puts $fileid "output_step = $output_step"

    # For results
    puts $fileid ""
    # On nodes results
    set cnrlist [list "Velocity" "Pressure" "Reactions" "Distance"]
    # set cnrlist [list "Velocity" "Pressure" "Reactions"]
    set nodal_results "nodal_results=\["
    foreach cnr $cnrlist {
     	set cxpath "$rootid//c.Results//c.OnNodes//i.${cnr}"
     	set cproperty "dv"
     	set cvalue [::xmlutils::setXml $cxpath $cproperty]
     	if {$cvalue =="Yes"} {
	    set cnkr [::xmlutils::getKKWord $kxpath $cnr]
     	    append nodal_results "\"$cnkr\","
     	}
     }
     set findcomma [string last "," $nodal_results]
     if {$findcomma !="-1"} {
     	set nodal_results [string range $nodal_results 0 end-1]
     	append nodal_results "\]" 
     	puts $fileid "$nodal_results"
     }

    # Set gauss_points_results to empty
    puts $fileid "gauss_points_results=\[\]"

    # WarnWinText "nodal_results:$nodal_results"
 
    # GiD post mode variables
    ::wkcf::WriteGiDPostMode $AppId $fileid 

    puts $fileid ""
    set PName [::KUtils::GetPaths "PName"]
    puts $fileid "problem_name=\"${PName}${AppId}\"" 
    puts $fileid "problem_path=\"$PDir\"" 
    
    # Get the kratos path 
    set cxpath "GeneralApplicationData//c.ProjectConfiguration//i.KratosPath"
    set cproperty "dv"
    set KratosPath [::xmlutils::setXml $cxpath $cproperty]
    set KratosPath [file native $KratosPath]

    # Write the kratos path
    puts $fileid "kratos_path=\"${KratosPath}\"" 
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
	    "Displacement"
	    {
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

    } elseif {$AnalysisType =="Linear"} {

	puts $fileid "Convergence_Tolerance = 1.0"
	puts $fileid "Absolute_Tolerance = 1.0"
	puts $fileid "Max_Iter = 1"
	set cConvergenceCriteria "Displacement_Criteria"
	puts $fileid "Convergence_Criteria = \"$cConvergenceCriteria\""
    }

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

    set shelllist [list "ShellIsotropic" "ShellAnisotropic" "EBST"]

    if {([info exists dprops($AppId,AllKElemId)]) && ($dprops($AppId,AllKElemId)>0)} {
	# For all defined kratos elements	
	foreach celemid $dprops($AppId,AllKElemId) {
	    if {$celemid in $shelllist} {
		set useshells "Yes"
		if {$celemid == "EBST"} {
		    set usenbst "Yes"
		    break
		}
	    }
	}
    }
    if {$usenbst =="Yes"} {
	puts $fileid "FindElementalNeighbours = \"True\""
    } else {
	puts $fileid "FindElementalNeighbours = \"False\""
    }
    if {$useshells =="Yes"} {
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
    puts $fileid ""
    set cgrlist [list "GreenLagrangeStrainTensor" "Rotations" "PK2StressTensor"]
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
    
    # Get the kratos path 
    set cxpath "GeneralApplicationData//c.ProjectConfiguration//i.KratosPath"
    set cproperty "dv"
    set KratosPath [::xmlutils::setXml $cxpath $cproperty]
    set KratosPath [file native $KratosPath]

    # Write the kratos path
    puts $fileid "kratos_path=\"${KratosPath}\"" 
 
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
    puts $fileid "import os, sys"
    puts $fileid ""
    puts $fileid "# Importing the Kratos Library"
    puts $fileid "from Kratos import *"
    puts $fileid "from KratosStructuralApplication import *"

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
	puts $fileid "    prop_id = $propid;"
	puts $fileid "    prop = Properties(prop_id)"
	# Write material model
	puts $fileid "    mat = $dprops($AppId,Material,$MatId,MatModel);"
      	puts $fileid "    prop.SetValue(CONSTITUTIVE_LAW, mat.Clone());"
	
    }

    close $fileid

}

















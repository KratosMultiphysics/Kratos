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
#     7.2- 04/04/14-G. Socorro, modify the proc WriteConditions to check that some hided  group exists (–AKGSkinMesh2D and -AKGSkinMesh3D)
#     7.1- 20/12/13-G. Socorro, update to use DEM and convection-diffusion applications
#     7.0- 19/09/13-G. Socorro, modify the proc WriteConditions to write condition for Quadrilateral surface element
#                               write convection-diffusion conditions
#     6.9- 18/09/13-G. Socorro, start to add the convection-diffusion properties 
#                               (add the proc WriteConvectionDiffusionPropertyAtNodes and WriteConvectionDiffusionPrescribedTemperatureBC)
#     6.8- 16/07/13-G. Socorro, modify the proc WriteBoundaryConditions to write the OutletPressure BC as a function of the solver type
#     6.7- 14/07/13-G. Socorro, modify the proc WriteBoundaryConditions to write walllaw BC
#     6.6- 18/06/13-G. Socorro, delete the call to the proc WritePythonGroupProperties
#     6.5- 17/06/13-G. Socorro, delete wmethod variable and all related procedures (*_m0,*_m1,*_m2) => now we are using only the new GiD groups
#     6.4- 24/04/13-G. Socorro, write rotational dofs boundary condition for beam element type
#     6.3- 19/03/13-G. Socorro, update the proc SelectPythonScript to select the correct python script
#     6.2- 17/12/12-J. Garate,  PFEM Wall is disabled for .mdpa
#     6.1- 21/11/12-J. Garate,  Modified ::wkcf::WriteNodalCoordinates_m2 ::wkcf::WriteElementConnectivities_m2 ::wkcf::WriteBoundaryCondition for PFEM
#     6.0- 26/11/12-J. Garate,  Support to PFEM Application
#     5.9- 12/11/12-J. Garate,  Fixed some errors
#     5.8- 07/11/12-J. Garate,  Added to namespace variable filechannel, New GiD method to write on files will use this channel the .mdpa file
#                               Modification and adaptation on functions: WriteCalculationFiles , WriteModelPartData, ::wkcf::WriteProperties
#                               WriteNodalCoordinates, WriteElementConnectivities, WriteBoundaryConditions, WriteConditions.
#                               Creation of functions using GiD_File fprintf $filechannel "%s" format
#     5.7- 22/10/12-J. Garate,  Massive group erasing correction
#     5.6- 10/10/12-J. Garate,  Adaptation for New GiD Groups
#     5.5- 10/10/12-G. Socorro, update the proc SelectPythonScript to select the the active structural analysis python script
#     5.4- 03/10/12-G. Socorro, add a call to the proc WriteFluidDistanceBC to write the distance variable
#     5.3- 24/09/12-G. Socorro, create a new variable "wbatfile" to control when write the bat file for Kratos 
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
    # ConvectionDiffusion application
    variable ConvectionDiffusionApplication
    # DEM application
    variable DEMApplication
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
    

    # To write the bat file
    variable wbatfile
    
    # filechannel -> output file .mdpa
    variable filechannel
    variable demfilechannel

    # Structural analysis condition counter
    variable sa_icondid
    
     # Debug/Release variable [0 => Debug, 1 => Release] Timers
    variable pflag
    set pflag 1
}

proc ::wkcf::WriteCalculationFiles {filename} {
    variable FluidApplication; variable StructuralAnalysis
    variable ActiveAppList;    variable filechannel
    variable ConvectionDiffusionApplication; variable DEMApplication
    variable demfilechannel
    variable meshgroupid    
    set meshgroupid 1

    # For debug
    if {!$::wkcf::pflag} {
        set inittime [clock seconds]
    }
    
    # Unset some local variables
    ::wkcf::UnsetLocalVariables
    
    # Init some namespace global variables
    ::wkcf::Preprocess
    
    # Write each block of the files *.mdpa
    # Rename the file name => Change .dat by .mdpa
    set basefilename "[string range $filename 0 end-4]"
    # Write a mdpa file for each application
    foreach AppId $ActiveAppList {
	if {$AppId eq "DEM"} {
	    set filename "${basefilename}${AppId}_Inlet.mdpa"
	    set demfilename "${basefilename}${AppId}.mdpa"
	    # Open the file
	    set demfilechannel [GiD_File fopen $demfilename]
	} else {
	    set filename "${basefilename}${AppId}.mdpa"
	}
	# wa "filename:$filename"
        
	# Open the file
	set filechannel [GiD_File fopen $filename]
	
        # Write model part data
        ::wkcf::WriteModelPartData $AppId
	
        # Write properties block
	::wkcf::WriteProperties $AppId
	
        # Write nodes block    
        ::wkcf::WriteNodalCoordinates $AppId
	
        # Write elements block 
        ::wkcf::WriteElementConnectivities $AppId
	
        # For fluid or convection-diffusion applications
        if {($AppId == "Fluid")} {
            # Write conditions (Condition2D and Condition3D)
            ::wkcf::WriteConditions $AppId  

        } elseif {($AppId == "ConvectionDiffusion")} {
	    # For convection-diffusion application

            # Write conditions (ThermalFace2D and ThermalFace3D)
	    ::wkcf::WriteConditions $AppId  
	    

            # Write initial conditions
            ::wkcf::WriteInitialConditions $AppId  
        }

	# Write boundary condition block
	::wkcf::WriteBoundaryConditions $AppId
	
	# For structural analysis application
	if {$AppId =="StructuralAnalysis"} {
	    # Write load properties block
	    ::wkcf::WriteLoads $AppId

	} elseif {$AppId == "Fluid"} {
	    # For fluid application
	    # Write nodal data for density and viscosity            
	    ::wkcf::WritePropertyAtNodes $AppId
	    
	    # Write all group properties
	    ::wkcf::WriteGroupMeshProperties $AppId
	    
	    # Write the cutting and point history properties
	    ::wkcf::WriteCutAndGraph $AppId
	    
	} elseif {$AppId == "ConvectionDiffusion"} {
	    
	    # Write nodal data for mesh velocity and velocity
	    ::wkcf::WriteConvectionDiffusionPropertyAtNodes $AppId

	    # Write the cutting and point history properties
	    ::wkcf::WriteCutAndGraph $AppId

	} elseif {$AppId == "DEM"} {
	    
	    # Write all group properties
	    ::wkcf::WriteGroupMeshProperties $AppId
	}
	
	# Close the file
	GiD_File fclose $filechannel
	if {$AppId eq "DEM"} {
	    GiD_File fclose $demfilechannel
	}
    }
    
    # Write the project parameters file
    ::wkcf::WriteProjectParameters
    
    
    # Select python scripts
    ::wkcf::SelectPythonScript
   
    # Write constitutive laws properties
    if {$StructuralAnalysis=="Yes"} {
	::wkcf::WriteConstitutiveLawsProperties
    }
   
    #  Write DEM explicit solver variables
    if {$AppId == "DEM"} {
	
	::wkcf::WriteExplicitSolverVariables
    }

    variable wbatfile
    if {$wbatfile} {
        # Write bat file only for the Linux OS
        if {($::tcl_platform(os) eq "Linux")} {
            ::wkcf::WriteBatFile $AppId
        }
    }
    
    # Unset some local variables
    ::wkcf::UnsetLocalVariables
    
    # For debug
    if {!$::wkcf::pflag} {
        set endtime [clock seconds]
        set ttime [expr {$endtime-$inittime}]
        WarnWinText "Total time: [::KUtils::Duration $ttime]"
    }
    
    return 1
}

proc ::wkcf::SelectPythonScript {} {
    # Select the correct Python script
    variable FluidApplication; variable StructuralAnalysis
    variable ConvectionDiffusionApplication
    
    set endfilename "KratosOpenMP.py"
    set mpiendfilename "KratosMPI.py"
    
    # Get the application root identifier    
    set rootdataid [::wkcf::GetApplicationRootId]
    
    set cproperty "dv"

    set PDir [file native [::KUtils::GetPaths "PDir"]]
    set PTDir [::KUtils::GetPaths "PTDir"]

    if {$StructuralAnalysis =="Yes"} {
	# Check for use OpenMP
	set cxpath "$rootdataid//c.SolutionStrategy//c.ParallelType//i.ParallelSolutionType"
	set ParallelSolutionType [::xmlutils::setXml $cxpath $cproperty]
	# wa "ParallelSolutionType:$ParallelSolutionType"
	
	# Solution type
	set cxpath "$rootdataid//c.AnalysisData//i.SolutionType"
	set SolutionType [::xmlutils::setXml $cxpath $cproperty]
	# wa "SolutionType:$SolutionType"
	
	if {$ParallelSolutionType eq "OpenMP"} {
	    # OpenMP
	    if {($SolutionType =="Dynamic")||($SolutionType =="Quasi-Static")||($SolutionType =="Pseudo-Dynamic")} {
		set ppfilename "KratosStructuralOpenMP.py"
	    } elseif {$SolutionType =="Static"} {
		set ppfilename "KratosStructuralOpenMP.py"
	    }
	
	    set fromfname [file native [file join "$PTDir/python" $ppfilename]]
	    set tofname [file native [file join $PDir $ppfilename]]
	
	    # Copy the script file
	    if {[catch {file copy -force "$fromfname" "$tofname"} error]} {
		WarnWin [= "Could not copy the Kratos Python script (%s) to (%s): Error (%s)" $fromfname $tofname $error ]
		return ""
	    }
	
	} elseif {$ParallelSolutionType eq "MPI"} {
	    # MPI
	    if {($SolutionType =="Dynamic")||($SolutionType =="Quasi-Static")||($SolutionType =="Pseudo-Dynamic")} {
		set mpifilename "KratosStructuralOpenMP.py"
	    } elseif {$SolutionType =="Static"} {
		set mpifilename "KratosStructuralOpenMP.py"
	    }
	    
	    set mpifromfname [file native [file join "$PTDir/python" $mpifilename]]
	    set mpitofname [file native [file join $PDir $mpifilename]]
	    
	    # Copy the script file
	    if {[catch {file copy -force "$mpifromfname" "$mpitofname"} error]} {
		WarnWin [= "Could not copy the Kratos Python script (%s) to (%s): Error (%s)" $mpifromfname $mpitofname $error ]
		return ""
	    }
	}
    }

    # For fluid application
    if {$FluidApplication =="Yes"} {
	# Free surface
	set cxpath "$rootdataid//c.AnalysisData//i.FreeSurface"
	set FreeSurface [::xmlutils::setXml $cxpath $cproperty]
	# wa "FreeSurface:$FreeSurface"

	# Solver type for free surface
	set cxpath "$rootdataid//c.AnalysisData//i.SolverTypeFreeSurf"
	set SolverTypeFreeSurf [::xmlutils::setXml $cxpath $cproperty]
	# wa "SolverTypeFreeSurf:$SolverTypeFreeSurf"
	
	# Get the fluid approach
	set cxpath "$rootdataid//c.AnalysisData//i.FluidApproach"
	set FluidApproach [::xmlutils::setXml $cxpath $cproperty]

	# wa "FluidApproach:$FluidApproach"

	# Check for use OpenMP
	set cxpath "$rootdataid//c.SolutionStrategy//c.ParallelType//i.ParallelSolutionType"
	set ParallelSolutionType [::xmlutils::setXml $cxpath $cproperty]
	
	if {$ParallelSolutionType eq "OpenMP"} {
	    
	    if {$FluidApproach eq "Eulerian"} {
		# Eulerian fluid case

		set ppfilename "KratosOpenMPFluid.py"
		if {($FreeSurface eq "Yes") && ($SolverTypeFreeSurf eq "LevelSet")} {
		    set ppfilename "KratosOpenMPFluidLevelSet.py"
		}
		
		set tofname [file native [file join $PDir $endfilename]]
		set fromfname [file native [file join $PTDir python $ppfilename]]
		
		# Copy the script file
		if {[catch {file copy -force -- "$fromfname" "$tofname"} error]} {
		    WarnWin [= "Could not copy the Kratos Python script (%s) to (%s): Error (%s)" $fromfname $tofname $error ]
		    return ""
		}

	    } elseif {$FluidApproach eq "PFEM-Lagrangian"} {
		# PFEM case

		set ppfilename "KratosOpenMPPFEM.py"
		set tofname [file native [file join $PDir $endfilename]]
		set fromfname [file native [file join $PTDir python $ppfilename]]
		
		# Copy the script file
		if {[catch {file copy -force -- "$fromfname" "$tofname"} error]} {
		    WarnWin [= "Could not copy the Kratos Python script (%s) to (%s): Error (%s)" $fromfname $tofname $error ]
		    return ""
		}
	    }
	    
	} elseif {$ParallelSolutionType eq "MPI"} {
	    set mpifilename "KratosMPIFluid.py"
	    
	    set mpitofname [file native [file join $PDir $mpiendfilename]]
	    set mpifromfname [file native [file join $PTDir python $mpifilename]]
 
	    if {[catch {file copy -force -- "$mpifromfname" "$mpitofname"} error]} {
		WarnWin [= "Could not copy the Kratos Python script (%s) to (%s): Error (%s)" $mpifromfname $mpitofname $error ]
		return ""
	    }
	}
    }

    # For convection-diffusion application
    if {$ConvectionDiffusionApplication =="Yes"} {
	# Check for use OpenMP
	set cxpath "$rootdataid//c.SolutionStrategy//c.ParallelType//i.ParallelSolutionType"
	set ParallelSolutionType [::xmlutils::setXml $cxpath $cproperty]
	
	if {$ParallelSolutionType eq "OpenMP"} {
	    
	    set ppfilename "KratosOpenMPConvDiff.py"
	     
	    set tofname [file native [file join $PDir $endfilename]]
	    set fromfname [file native [file join $PTDir python $ppfilename]]
	    
	    # Copy the script file
	    if {[catch {file copy -force -- "$fromfname" "$tofname"} error]} {
		WarnWin [= "Could not copy the Kratos Python script (%s) to (%s): Error (%s)" $fromfname $tofname $error ]
		return ""
	    }
	} 
    }
}

proc ::wkcf::WriteModelPartData {AppId} {
    # Write the model part data
    # Arguments
    # AppId => Application identifier
    variable filechannel; variable demfilechannel
     
    if {$AppId eq "DEM"} {
    	GiD_File fprintf $demfilechannel "%s" "Begin ModelPartData"
	GiD_File fprintf $demfilechannel "%s" "//  VARIABLE_NAME value"
	GiD_File fprintf $demfilechannel "%s" "End ModelPartData"
	GiD_File fprintf $demfilechannel "%s" ""
    } else {
	GiD_File fprintf $filechannel "%s" "Begin ModelPartData"
	GiD_File fprintf $filechannel "%s" "//  VARIABLE_NAME value"
	GiD_File fprintf $filechannel "%s" "End ModelPartData"
	GiD_File fprintf $filechannel "%s" ""
    }
}   

proc ::wkcf::WriteProperties {AppId} {
    # Write the properties block
    # Arguments
    # AppId => Application identifier
    variable dprops;  variable ndime
    
    variable filechannel
    
    # For structural analysis application
    if {$AppId =="StructuralAnalysis"} {
        set propid 0
        foreach PropertyId $dprops($AppId,GKProps,AllPropertyId) {
            incr propid 1
            GiD_File fprintf $filechannel "%s" "Begin Properties $propid // GUI property identifier: $PropertyId"
            
            # Get the material identifier for this property 
            set MatId $dprops($AppId,Property,$PropertyId,MatId) 
            # wa "PropertyId:$PropertyId MatId:$MatId"
            GiD_File fprintf $filechannel "%s" "// GUI material identifier: $MatId"
            
            # Write material properties
            foreach cpropid $dprops($AppId,Material,$MatId,Props) {
                # WarnWinText "material propid:$cpropid"
                lassign $cpropid key value
                GiD_File fprintf $filechannel "%s" " $key $value"
            } 
            
            # Write section (others) properties (thickness, etc.)
            if {[info exists dprops($AppId,Material,$PropertyId,CProps)]} {
                foreach cpropid $dprops($AppId,Material,$PropertyId,CProps) {
                    # wa "cpropid:$cpropid"
                    lassign $cpropid cvaluetype ckeyword cvalue
                    # wa "cvaluetype:$cvaluetype ckeyword:$ckeyword cvalue:$cvalue"
                    #  Check the data type
                    if {$cvaluetype eq "Scalar"} {
                        # Scalar
                        GiD_File fprintf $filechannel "%s" " $ckeyword $cvalue"
                    } elseif {$cvaluetype eq "Matrix"} {
                        # Matrix
                        lassign $cvalue matdim matvalue 
                        GiD_File fprintf $filechannel "%s" " $ckeyword ${matdim} $matvalue"
                    }
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
                            GiD_File fprintf $filechannel "%s" "// GUI body force group identifier: $GroupId"
                            ::wkcf::WriteBodyForceValues $cprop
                        }
                    }
                }
            } else {
                # Write the body forces with default values
                set cprop [list [list GravityValue 9.8] [list Cx 0.0] [list Cy 0.0] [list Cz 0.0]]
                ::wkcf::WriteBodyForceValues $cprop
            }
            
	    GiD_File fprintf $filechannel "%s" "End Properties"
        }
       
        GiD_File fprintf $filechannel ""

    } elseif {$AppId =="Fluid"} {
	# For fluid application
   
        GiD_File fprintf $filechannel "%s" "Begin Properties 0"
        GiD_File fprintf $filechannel "%s" "End Properties"
        GiD_File fprintf $filechannel ""

    } elseif {$AppId =="ConvectionDiffusion"} {
	# For convection-diffusion application

	# Kratos key word xpath
	set kxpath "Applications/$AppId"
	set cproperty "dv"
	
	# Only in the case of linear analysis
	# Analysis type
	set cxpath "$AppId//c.AnalysisData//i.AnalysisType"
	set AnalysisType [::xmlutils::setXml $cxpath $cproperty]

	# Check the face flux case
	set facefluxflag [expr {([info exists dprops($AppId,GBCKProps,AllPropertyId)]) && ([llength $dprops($AppId,GBCKProps,AllPropertyId)])}]
	if {$facefluxflag} {

	    if {$AnalysisType eq "Linear"} {
		
		# Write the base properties => Properties identifier =0
		::wkcf::CDWriteDefaultProperties $AppId
	    }
	    
	    # Write boundary condition properties 
	    foreach PropertyId $dprops($AppId,GBCKProps,AllPropertyId) {
		# Get the propertyid-group identifier
		set PropertyGroupId $dprops($AppId,GBCKProps,$PropertyId,PropertyGroupId)
		GiD_File fprintf $filechannel "%s" "Begin Properties $PropertyId // GUI face heat flux BC group identifier: $PropertyGroupId"
		
		# Write face heat flux BC properties 
		foreach cpropid $dprops($AppId,GBCKProps,$PropertyId,PropertyList) {
		    # wa "Properties propid:$cpropid"
		    lassign $cpropid key value
		    if {$value ne ""} {
			# Get the kratos keyword
			set ckword [::xmlutils::getKKWord $kxpath $key]
			# wa "ckword:$ckword"
			if {$ckword !=""} {
			    GiD_File fprintf $filechannel "%s" " $ckword $value"
			}
		    }
		}
		
		GiD_File fprintf $filechannel "%s" "End Properties"
		GiD_File fprintf $filechannel ""
	    }
	} else {
	    
	     # Write the base properties => Properties identifier =0
	    ::wkcf::CDWriteDefaultProperties $AppId	   
	}

    } elseif {$AppId =="DEM"} {
	
	# Kratos key word xpath
	set kxpath "Applications/$AppId"

        foreach cgroupid $dprops($AppId,AllMeshGroupId) {
	    # Get the mesh-group identifier
	    set meshgroupid $dprops($AppId,Mesh,$cgroupid,MeshIdGroup) 
            GiD_File fprintf $filechannel "%s" "Begin Properties $meshgroupid // GUI property identifier: $cgroupid"
            
	    # Write inlet properties 
            foreach cpropid $dprops($AppId,Mesh,$cgroupid,MeshIdGroupProp) {
		# WarnWinText "material propid:$cpropid"
                lassign $cpropid key value
		if {$value ne ""} {
		    if {$key eq "Material"} {
			GiD_File fprintf $filechannel "%s" "// GUI inlet material identifier: $value"
			# Get the material properties
			# Get all material properties
			set mpxpath "[::KMat::findMaterialParent $value]//m.[list ${value}]"
			# WarnWinText "mpxpath:$mpxpath"
			# Get all the properties
			set allmatprops [::xmlutils::setXmlContainerPairs $mpxpath "" "value" "Property" "mat"]
			# wa "allmatprops:$allmatprops"
			foreach propid $allmatprops {
			    lassign $propid key value
			    # wa "key:$key value:$value"
			    
			    # Get the kratos keyword
			    set ckword [::xmlutils::getKKWord $kxpath $key]
			    # wa "ckword:$ckword"
			    if {$ckword !=""} {
				GiD_File fprintf $filechannel "%s" " $ckword $value"
			    }
			}
		

		    } else { 
			# Get the kratos keyword
			set ckword [::xmlutils::getKKWord $kxpath $key]
			# wa "ckword:$ckword"
			if {$ckword !=""} {
			    GiD_File fprintf $filechannel "%s" " $ckword $value"
			}
		    }
		}
            }

	    GiD_File fprintf $filechannel "%s" "End Properties"
	    GiD_File fprintf $filechannel ""
	}
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
    variable ndime; variable dprops
    variable filechannel

    # For debug
    if {!$::wkcf::pflag} {
        set inittime [clock seconds]
    }
    
    # Check for all defined kratos elements
    if {([info exists dprops($AppId,AllKElemId)]) && ([llength $dprops($AppId,AllKElemId)])} {
       
	# For all defined kratos nodes        
	foreach celemid $dprops($AppId,AllKElemId) {
            # Check for all defined group identifier for this element
            if {[info exists dprops($AppId,KElem,$celemid,AllGroupId)] && [llength $dprops($AppId,KElem,$celemid,AllGroupId)]} {
                # For all defined group identifier for this element
                foreach cgroupid $dprops($AppId,KElem,$celemid,AllGroupId) {
                    # Get the GiD entity type, element type and property identifier
                    lassign $dprops($AppId,KElem,$celemid,$cgroupid,GProps) GiDEntity GiDElemType PropertyId KEKWord nDim
                    if { [GiD_EntitiesGroups get $cgroupid nodes -count] } {
                        # Write all nodes for this group in increasing orden
                        GiD_File fprintf $filechannel "Begin Nodes // GUI group identifier: $cgroupid"			
                        if {$ndime == "2D"} {			    
                            foreach nodeid [GiD_EntitiesGroups get $cgroupid nodes] {                               
                                lassign [GiD_Mesh get node $nodeid] layer x y z
                                # msg "nodeid:$nodeid layer:$layer x:$x y:$y z:$z"
                                GiD_File fprintf $filechannel "$nodeid [format {%.5f %.5f} $x $y] 0"
                            }
                        } else {			            
                            foreach nodeid [GiD_EntitiesGroups get $cgroupid nodes] {                                
                                lassign [GiD_Mesh get node $nodeid] layer x y z
                                # msg "nodeid:$nodeid layer:$layer x:$x y:$y z:$z"
                                GiD_File fprintf $filechannel "$nodeid [format {%.5f %.5f %.5f} $x $y $z]"                 
                            }
                        }
			GiD_File fprintf $filechannel "End Nodes"
			GiD_File fprintf $filechannel ""
                    }
		}
            }
        }
    }
    
    # Special case of DEM application
    if {$AppId eq "DEM"} {
	# Get the values
	set rootid $AppId
	set basexpath "$rootid//c.Conditions//c.DEM-Inlet"
	set gproplist [::xmlutils::setXmlContainerIds $basexpath]
	# wa "gproplist:$gproplist"
	foreach cgroupid $gproplist {
	    # Get the group node list
	    set nlist [::wkcf::GetInletGroupNodes $AppId $cgroupid]
	    if {[llength $nlist]} {
		# Write all nodes for this group in increasing orden
		GiD_File fprintf $filechannel "Begin Nodes // GUI Inlet group identifier: $cgroupid"			
		if {$ndime == "2D"} {			    
		    foreach nodeid $nlist {                               
			lassign [GiD_Mesh get node $nodeid] layer x y z
			# msg "nodeid:$nodeid layer:$layer x:$x y:$y z:$z"
			GiD_File fprintf $filechannel "$nodeid [format {%.5f %.5f} $x $y] 0"
		    }
		} else {			            
		    foreach nodeid $nlist {                                
			lassign [GiD_Mesh get node $nodeid] layer x y z
			# msg "nodeid:$nodeid layer:$layer x:$x y:$y z:$z"
			GiD_File fprintf $filechannel "$nodeid [format {%.5f %.5f %.5f} $x $y $z]"                 
		    }
		}
		GiD_File fprintf $filechannel "End Nodes"
		GiD_File fprintf $filechannel ""
	    }
	}

	# Update the standard mdpa file
	variable demfilechannel
	
	GiD_File fprintf $demfilechannel "Begin Nodes"
	GiD_File fprintf $demfilechannel "End Nodes"

    }

    # For debug
    if {!$::wkcf::pflag} {
        set endtime [clock seconds]
        set ttime [expr $endtime-$inittime]
        # WarnWinText "endtime:$endtime ttime:$ttime"
        WarnWinText "Write nodal coordinates: [::KUtils::Duration $ttime]"
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
    variable useqelem; variable dprops
    variable filechannel
    
    # For debug
    if {!$::wkcf::pflag} {
        set inittime [clock seconds]
    }
    
    # Check for all defined kratos elements
    if {([info exists dprops($AppId,AllKElemId)]) && ([llength $dprops($AppId,AllKElemId)])} {
        set kwxpath "Applications/$AppId"        
        # For all defined kratos elements        
        foreach celemid $dprops($AppId,AllKElemId) {
	    # wa "celemid:$celemid"
            # Check for all defined group identifier for this element
            if {[info exists dprops($AppId,KElem,$celemid,AllGroupId)] && [llength $dprops($AppId,KElem,$celemid,AllGroupId)]} {
                # For all defined group identifier for this element
                foreach cgroupid $dprops($AppId,KElem,$celemid,AllGroupId) {
                    # Get the GiD entity type, element type and property identifier
                    lassign $dprops($AppId,KElem,$celemid,$cgroupid,GProps) GiDEntity GiDElemType PropertyId KEKWord nDim
		    # wa "cgroupid:$cgroupid GiDEntity:$GiDEntity GiDElemType:$GiDElemType PropertyId:$PropertyId KEKWord:$KEKWord nDim:$nDim"
                    if {[GiD_EntitiesGroups get $cgroupid elements -count -element_type $GiDElemType]} {
                        set etbf ""
                        set usennode [::xmlutils::getKKWord $kwxpath $celemid usennode]                                       
                        if {$usennode eq "Yes"} {
                            set etbf [::wkcf::GetnDimnNode $GiDElemType $nDim]
                        }   
                        set kelemtype [string trim ${KEKWord}${etbf}]
                        set GlobalPId $dprops($AppId,KElem,$celemid,$cgroupid,GlobalPId)
                        GiD_File fprintf $filechannel "Begin Elements $kelemtype   // GUI group identifier: $cgroupid"
                        foreach elem_id [GiD_EntitiesGroups get $cgroupid elements -element_type $GiDElemType] {
                            GiD_File fprintf $filechannel "$elem_id $GlobalPId [lrange [GiD_Mesh get element $elem_id] 3 end]"
                        }
                        GiD_File fprintf $filechannel "End Elements"
                        GiD_File fprintf $filechannel ""
                    }
                }
            }
        }
        GiD_File fprintf $filechannel ""
    }
    
    # For debug
    if {!$::wkcf::pflag} {
        set endtime [clock seconds]
        set ttime [expr $endtime-$inittime]
        # WarnWinText "endtime:$endtime ttime:$ttime"
        WarnWinText "Write element connectivities: [::KUtils::Duration $ttime]"
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
			set usebeams "No"
			set beamlist [list "BeamElement"]
                        set shelllist2 [list "EBST" "ShellThick" "ShellThickCR" "ShellThin" "ShellThinCR"]
                        if {([info exists dprops($AppId,AllKElemId)]) && ($dprops($AppId,AllKElemId)>0)} {
                            # For all defined kratos elements        
                            foreach celemid $dprops($AppId,AllKElemId) {
                                if {$celemid in $shelllist2} {
                                    set useshells2 "Yes"
                                } elseif {$celemid in $beamlist} {
                                    set usebeams "Yes"
				}
                            }
                        }
                        if {($useshells2 eq "Yes")||($usebeams eq "Yes")} {
                            set kwxpath "Applications/StructuralAnalysis"
                            set kwordlist [list [::xmlutils::getKKWord $kwxpath "Rx"] [::xmlutils::getKKWord $kwxpath "Ry"] [::xmlutils::getKKWord $kwxpath "Rz"]]
                            # Process displacement properties
                            ::wkcf::WriteDispRotBC $AppId $ccondid $kwordlist 
                        }
                    }
                    "OutletPressure" {
			if {$AppId=="Fluid"} {
			    set cproperty "dv"
			    set cxpath "$AppId//c.AnalysisData//i.FluidApproach"
			    set FluidApproach [::xmlutils::setXml $cxpath $cproperty]
			    if { $FluidApproach eq "Eulerian" } {
				# Solver type
				set cxpath "$AppId//c.AnalysisData//i.SolverType"
				set SolverType [::xmlutils::setXml $cxpath $cproperty]
				# WarnWinText "SolverType:$SolverType"
				set kwxpath "Applications/$AppId"

				switch -exact -- $SolverType {
				    "ElementBased" {
					set kwid "OutletPressureFractionalStep"
					set kwordlist [list "[::xmlutils::getKKWord $kwxpath $kwid]"]
				    }
				    "Monolithic" {
					set kwid "OutletPressureMonolithic"
					set kwordlist [list "[::xmlutils::getKKWord $kwxpath $kwid]"]
				    }
				}
							                       
				# Process outlet pressure
				::wkcf::WriteOutLetPressureBC $AppId $ccondid $kwordlist 
			    }
			}
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
                    "Is-Slip" - "WallLaw" {
                        # Write is-slip/walllaw boundary condition
                        set kwordlist [list "IS_STRUCTURE" "Y_WALL"]
                        ::wkcf::WriteFluidIsSlipWallLawBC $AppId $ccondid $kwordlist
                    }
		    "Distance" {
                        # Write distance boundary condition
                        set kwordlist [list "DISTANCE"]
                        ::wkcf::WriteFluidDistanceBC $AppId $ccondid $kwordlist
                    }
                    "PFEMWall--" {
                        # Commented by J. Garate on 17/12/2012
                        set kwordlist [list "LINEAR_VELOCITY_X" "LINEAR_VELOCITY_Y" "LINEAR_VELOCITY_Z" "ANGULAR_VELOCITY_X" "ANGULAR_VELOCITY_Y" "ANGULAR_VELOCITY_Z" ]
                        ::wkcf::WriteFluidPFEMWallBC $AppId $ccondid $kwordlist
                    }
                    "PFEMFluidInlet" {
                        # Write PFEM Fluid Velocity
                        set kwordlist [list "VELOCITY_X" "VELOCITY_Y" "VELOCITY_Z"]
                        ::wkcf::WriteFluidPFEMInletBC $AppId $ccondid $kwordlist
                    }
		    "PrescribedTemperature" {
			# Write prescribed temperature condition
			set kwordlist [list "TEMPERATURE"]
                        ::wkcf::WriteConvectionDiffusionPrescribedTemperatureBC $AppId $ccondid $kwordlist
		    }
		    "HeatFlux" {
			# Write prescribed heat flux condition
			set kwordlist [list "HEAT_FLUX"]
                        ::wkcf::WriteConvectionDiffusionPrescribedHeatFluxBC $AppId $ccondid $kwordlist
		    }
		    "FaceHeatFlux" {
			# Write prescribed face heat flux condition
			set kwordlist [list "FACE_HEAT_FLUX"]
                        ::wkcf::WriteConvectionDiffusionPrescribedFaceHeatFluxBC $AppId $ccondid $kwordlist
		    }
		    "PFEMFixedWall" {
                        # Write PFEM fixed wall boundary condition
                        set kwordlist [list "IS_STRUCTURE" "DISPLACEMENT_X" "DISPLACEMENT_Y" "DISPLACEMENT_Z"]
                        ::wkcf::WritePFEMLagrangianFluidFixedWallBC $AppId $ccondid $kwordlist
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

proc ::wkcf::WriteInitialConditions {AppId} {
    # Write the initial condition block
    variable dprops
    
    # Check for all defined initial condition type
    if {([info exists dprops($AppId,AllICTypeId)]) && ([llength $dprops($AppId,AllICTypeId)])} {
        # For debug
        if {!$::wkcf::pflag} {
            set inittime [clock seconds]
        }
        set inletvelglist [list]; set noslipglist [list]
        set flagvariablelist [list] 
        # For all defined initial condition identifier
        foreach ccondid $dprops($AppId,AllICTypeId) {
	    # WarnWinText "ccondid:$ccondid"
            # Check for all defined group identifier inside this initial condition type
            if {([info exists dprops($AppId,IC,$ccondid,AllGroupId)]) && ([llength $dprops($AppId,IC,$ccondid,AllGroupId)])} {
                # Select the initial condition type
                switch -exact -- $ccondid {
                    "InitialTemperature" {
                        set kwxpath "Applications/$AppId"
                        set kwordlist [list [::xmlutils::getKKWord $kwxpath "$ccondid" "kkword"]]
                        # Process initial temperature properties
                        ::wkcf::WriteInitialTemperatureIC $AppId $ccondid $kwordlist 
                    }
		}
            }
        }
       
        # For debug
        if {!$::wkcf::pflag} {
            set endtime [clock seconds]
            set ttime [expr $endtime-$inittime]
            # WarnWinText "endtime:$endtime ttime:$ttime"
            WarnWinText "Write initial conditions: [::KUtils::Duration $ttime]"
        }
    }
}

proc ::wkcf::WriteConditions {AppId} {
    # ABSTRACT: Write condition properties
    variable dprops; variable ndime
    variable ctbclink

    variable filechannel

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

	# Select the application type
	# Fluid
	if {$AppId == "Fluid"} {
	    set cxpath "$AppId//c.AnalysisData//i.FluidApproach"
	    set FluidApproach [::xmlutils::setXml $cxpath $cproperty]
	    if { $FluidApproach eq "Eulerian" } {
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
	    } elseif {$FluidApproach eq "PFEM-Lagrangian"} {
		set ConditionId "Condition${ndime}"
	    }
	} elseif {$AppId == "ConvectionDiffusion"} {
	    
	    # Convection diffusion
	    set ConditionId "ThermalFace${ndime}"  
	}
        
        if {$ndime =="2D"} {
	    # 2D
	    set cgroupid "-AKGSkinMesh2D"
            set GiDElemType "Linear"

	    if {[GiD_Groups exists $cgroupid]} {
    
	      if {[GiD_EntitiesGroups get $cgroupid elements -count -element_type $GiDElemType]} {
		  # Write all conditions for 2D case
		  GiD_File fprintf $filechannel "%s" "Begin Conditions $ConditionId"
		  set condid 0
		  
		  foreach elem_id [GiD_EntitiesGroups get $cgroupid elements -element_type $GiDElemType] {
		      incr condid 1 
		      set nodes [lrange [GiD_Mesh get element $elem_id] 3 end]
		      # nodes = [list nodei nodej]
		      GiD_File fprintf $filechannel "%10d %4d %10d %10d" $condid $fixval [lindex $nodes 0] [lindex $nodes 1]
		      # Update the link between the condition id. and the BC element id
		      dict set ctbclink $elem_id $condid
		  }
		  
		  GiD_File fprintf $filechannel "%s" "End Conditions"
		  GiD_File fprintf $filechannel ""
	      }
	    }
        } elseif {$ndime =="3D"} {

	    # 3D
	    set cgroupid "-AKGSkinMesh3D"
            set GiDElemType "Triangle"
	    set usetriangle 0

    	    if {[GiD_Groups exists $cgroupid]} {
 
		if {[GiD_EntitiesGroups get $cgroupid elements -count -element_type $GiDElemType]} {
		    # Write all conditions for 3D case
		    GiD_File fprintf $filechannel "%s" "Begin Conditions $ConditionId"
		    set condid 0
		    
		    foreach elem_id [GiD_EntitiesGroups get $cgroupid elements -element_type $GiDElemType] {
			incr condid 1 
			set nodes [lrange [GiD_Mesh get element $elem_id] 3 end]
			set ni [lindex $nodes 0]
			set nj [lindex $nodes 1]
			set nk [lindex $nodes 2]
			# msg "$ni $nj $nk $condid"
			GiD_File fprintf $filechannel "%10d %4d %10d %10d %10d" $condid $fixval $ni $nj $nk
			# Update the link between the condition id. and the BC element id
			dict set ctbclink $elem_id $condid
		    }
		    GiD_File fprintf $filechannel "%s" "End Conditions"
		    GiD_File fprintf $filechannel ""
		    set usetriangle 1
		}
    
		set GiDElemType "Quadrilateral"
    
		if {[GiD_EntitiesGroups get $cgroupid elements -count -element_type $GiDElemType]} {
		    # Write all conditions for 3D case
		    GiD_File fprintf $filechannel "%s" "Begin Conditions $ConditionId"
		    if {!$usetriangle} {
			set condid 0
		    }
		    foreach elem_id [GiD_EntitiesGroups get $cgroupid elements -element_type $GiDElemType] {
			incr condid 1 
			set nodes [lrange [GiD_Mesh get element $elem_id] 3 end]
			set ni [lindex $nodes 0]
			set nj [lindex $nodes 1]
			set nk [lindex $nodes 2]
			set nl [lindex $nodes 3]
			# msg "$ni $nj $nk $nl $condid"
			GiD_File fprintf $filechannel "%10d %4d %10d %10d %10d %10d" $condid $fixval $ni $nj $nk $nl
			# Update the link between the condition id. and the BC element id
			dict set ctbclink $elem_id $condid
		    }
		    GiD_File fprintf $filechannel "%s" "End Conditions"
		    GiD_File fprintf $filechannel ""
		}
	    }
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
    variable ConvectionDiffusionApplication

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
    
    # Write fluid application project parameter file
    if {$FluidApplication =="Yes"} {
	::wkcf::WriteFluidProjectParameters "Fluid" $fileid $PDir
    }
    
    # Write convection diffusion application project parameter file
    if {$ConvectionDiffusionApplication =="Yes"} {
	::wkcf::WriteConvectionDiffusionProjectParameters "ConvectionDiffusion" $fileid $PDir
    }

    close $fileid
    
}

proc ::wkcf::WriteGiDPostMode {AppId fileid} {
    # Write the GiD post mode variables for each applications
    
    # kratos key word xpath
    set kwxpath "Applications/$AppId"
    
    # Gid results
    set gidrlist [list "GiDPostMode" "GiDWriteMeshFlag" "GiDWriteConditionsFlag" "GiDWriteParticlesFlag" "GiDMultiFileFlag"]
    foreach gidr $gidrlist {
        # Get the value
        set cxpath "$AppId//c.Results//c.GiDOptions//i.[list ${gidr}]"
        set cproperty "dv"
        set cvalue [::xmlutils::setXml $cxpath $cproperty]
	if {$cvalue !=""} {
	    # Get the kratos keyword
	    set gidrkw [::xmlutils::getKKWord $kwxpath $gidr]
	    # WarnWinText "gidr:$gidr cvalue:$cvalue gidrkw:$gidrkw"
	    if {($gidr=="GiDWriteMeshFlag") || ($gidr=="GiDWriteConditionsFlag") || ($gidr=="GiDWriteParticlesFlag")} {
		if {$cvalue =="Yes"} {
		    set cvalue True
		} else {
		    set cvalue False
		}
		# wa "cvalue:$cvalue"
		puts $fileid "$gidrkw = $cvalue"
	    } else {
		puts $fileid "$gidrkw = \"$cvalue\""
	    }
	}
    }
}

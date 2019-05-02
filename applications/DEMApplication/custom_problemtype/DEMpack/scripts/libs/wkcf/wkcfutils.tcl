###############################################################################
#
#    NAME: wkcfutils.tcl
#
#    PURPOSE: Some utilities procedures to write all the Kratos calculation files
#
#    QUANTECH ATZ-DEVELOPMENT DEPARTMENT
#
#    AUTHOR : G. Socorro
#
#    CREATED AT: 10/05/10
#
#    HISTORY:
#	  4.1- 08/05/18-M. Lopez,   create a new proc FindBoundariesOfNonSphericElements to substitute the old FindBoundaries_no_spheric_elems
#     4.0- 10/10/13-G. Socorro, create a new proc GetMaterialPropertiesFromAttributes to get some material properties using
#                                the container CLawProperties in the kratos_key_words.xml file
#                               - Modify the proc GetMaterialProperties to write the the properties of the material model "HyperElastic-Plastic"
#                               - Add some utilities proc and variable for Convection-Diffusion and DEM applications
#     3.9- 03/10/13-G. Socorro, correct a bug in the proc GetSurfaceTypeList (special case of default element type for volumes)
#     3.8- 22/09/13-G. Socorro, add the proc GetElementIdFromPropertyId to get the list of element identifier from a property identifier
#     3.7- 19/09/13-G. Socorro, add the proc GetSurfaceTypeList to classify the surface type
#     3.6- 18/09/13-G. Socorro, add some proc and variable for the convection-diffusion application
#     3.5- 14/07/13-G. Socorro, modify the proc GetBoundaryConditionProperties to enable the condition WallLaw
#     3.4- 17/06/13-G. Socorro, delete wmethod variable => now we are using only the new GiD groups
#     3.3- 17/05/13-G. Socorro, add the proc SumInertia, modify the proc GetCrossSectionProperties to include the new cross section properties (database)
#     3.2- 11/02/13-G. Socorro, correct a bug in the proc GetCrossSectionProperties when write the matrix properties (delete a parenthesis was left)
#     3.1- 17/12/12-J. Garate,  Disabled the Wall PFEM .mdpa Write. Add Beam 2D
#     3.0- 05/12/12-J. Garate,  PFEM Slip velocity format correction
#     2.9- 03/12/12-J. Garate,  Corrected a Bug on Angular Velocity for PFEM
#     2.8- 28/11/12-J. Garate,  Update ::wkcf::GetBoundaryConditionProperties for PFEM
#     2.7- 07/11/12-J. Garate,  Preprocess now selects the wmethod value { 0, 1, 2 } depending on GiD_Version.
#                               Added the CleanAutomaticConditionGroupGiD and AssignConditionToGroupGID functions for "Before Mesh Generation" Event Procedure
#     2.6- 10/10/12-G. Socorro, update the proc GetPropertiesData to assign the cross section properties for all constitutive laws
#     2.5- 09/10/12-G. Socorro, add the proc GetCrossSectionProperties, modify other procs to include the structural
#                               analysis cross section properties options
#     2.4- 04/10/12-G. Socorro, add the proc GetPorousZonesProperties to get the properties of the porous zones
#     2.3- 02/10/12-G. Socorro, update the proc GetBoundaryConditionProperties to include the Distance BC
#     2.2- 24/09/12-G. Socorro, init the new variable "wbatfile" to 0
#     2.1- 21/09/12-G. Socorro, update the proc WriteBatFile to write the bat file for Linux OS
#     2.0- 23/07/12-G. Socorro, add the "ConstantValue" to the Is-Slip condition
#     1.9- 13/05/12-G. Socorro, set/unset the local variable ctbclink (condition to bc linker)
#     1.8- 07/05/12-G. Socorro, add the CleanAutomatic, CleanAutomaticConditionGroup and AssignConditionToGroup
#     1.7- 20/04/12-G. Socorro, update the proc GetBoundaryConditionProperties to include is-slip and walllaw boundary conditions
#     1.6- 10/03/12-G. Socorro, pass some proc to fluid or structural analysis script (WriteFluidSolvers,etc.)
#     1.5- 27/03/12-G. Socorro, modify some structural analysis application properties (constitutive modeling)
#     1.4- 20/02/12-J. Garate,  añadida una opcion para corregir un bug con Solid 2D
#     1.3- 30/06/11-G. Socorro, add the new condition Flag-Variable to the fluid application
#     1.2- 07/06/11-G. Socorro, correct a bug when defined many properties which different thickness => proc GetMaterialProperties
#     1.1- 24/05/11-G. Socorro, update some procedure to use the membrane element => Use thickness property
#     1.0- 09/12/10-G. Socorro, add thickness for "Elastic-Isotropic" material model and shell element
#     0.9- 07/09/10-G. Socorro, correct an error when defined body force for group of element with the same property
#     0.8- 06/09/10-G. Socorro, check for active group variable when get the element properties
#     0.7- 03/09/10-G. Socorro, correct an error with the BC for inlet and no-slip
#     0.6- 16/06/10-G. Socorro, add the elasto-plastic material model to the constitutive laws
#     0.5- 15/06/10-G. Socorro, add the damage material model to the constitutive laws
#     0.4- 11/06/10-G. Socorro, update material properties using the new constitutive equation
#     0.3- 09/06/10-G. Socorro, add FSI application
#     0.3- 31/05/10-G. Socorro, get all kratos key word from a XML file (kratos_key_words.xml)
#     0.2- 12/05/10-G. Socorro, start to add the fluid application options
#     0.1- 10/05/10-G. Socorro, create the base source code from wkcf.tcl
#
###############################################################################

proc ::wkcf::Preprocess {} {
    # Create some global variables used when write data to the file
    # Get the problem dimension => spatial dimension
    variable ndime;    variable StructuralAnalysis
    variable gidetype; variable useqelem
    variable dprops;   variable FluidApplication
    variable FSIApplication; variable ActiveAppList
  variable DSOLID

    # Check for use quadratic elements
    set useqelem [GiD_Info Project Quadratic]
    # WarnWinText "useqelem:$useqelem"

    # Init the structural analysis condition counter
    variable sa_icondid
    set sa_icondid 0

    # Get the spatial dimension
    set ndime [::xmlutils::GetSpatialDimension]

    # Get application type
    # Structural analysis
    set cxpath "GeneralApplicationData//c.ApplicationTypes//i.StructuralAnalysis"
    set cproperty "dv"
    set StructuralAnalysis [::xmlutils::setXml $cxpath $cproperty]

    # WarnWinText "StructuralAnalysis:$StructuralAnalysis ndime:$ndime"

    # Fluid application
    set cxpath "GeneralApplicationData//c.ApplicationTypes//i.Fluid"
    set FluidApplication [::xmlutils::setXml $cxpath $cproperty]


    # WarnWinText "FluidApplication:$FluidApplication"

    # FSI application
    set cxpath "GeneralApplicationData//c.ApplicationTypes//i.FluidStructureInteraction"
    set FSIApplication [::xmlutils::setXml $cxpath $cproperty]

    # WarnWinText "FSIApplication:$FSIApplication"

    # DEM application
    variable DEMApplication
    set cxpath "GeneralApplicationData//c.ApplicationTypes//i.DEM"
    set DEMApplication [::xmlutils::setXml $cxpath $cproperty]

    # WarnWinText "DEMApplication:$DEMApplication"

  # DEMSOLID application
  variable DSOLID
  set cxpath "GeneralApplicationData//c.ApplicationTypes//i.DSOLID"
  set DSOLID [::xmlutils::setXml $cxpath $cproperty]
    # Convection diffusion
  variable ConvectionDiffusionApplication
  set cxpath "GeneralApplicationData//c.ApplicationTypes//i.ConvectionDiffusion"
  set ConvectionDiffusionApplication [::xmlutils::setXml $cxpath $cproperty]

    # WarnWinText "ConvectionDiffusionApplication:$ConvectionDiffusionApplication"

    # Update active application list

  if {$FSIApplication =="Yes"} {
	set ActiveAppList [list "StructuralAnalysis" "Fluid"]
	# set structural analysis and fluid application to Yes
	set StructuralAnalysis "Yes"
	set FluidApplication "Yes"
    } else {
	# Only one application is active at the same time
	# Structural analysis
	if {$StructuralAnalysis =="Yes"} {
	    set ActiveAppList [list "StructuralAnalysis"]
	}
	# Fluid

	if {$FluidApplication =="Yes"} {
	    set ActiveAppList [list "Fluid"]
	}
	# Convection diffusion
	if {$ConvectionDiffusionApplication eq "Yes"} {
	    set ActiveAppList [list "ConvectionDiffusion"]
	}
	# DEM application
	if {($DEMApplication eq "Yes")&& ($FluidApplication =="Yes")} {
	    set ActiveAppList [list "Fluid" "DEM"]
	} elseif {($DEMApplication eq "Yes")&& ($FluidApplication =="No")} {
	    set ActiveAppList [list "DEM"]
	}
  # DEM+SOLID application
	if {($DEMApplication eq "Yes")&& ($DSOLID =="Yes")&& ($FluidApplication =="No")} {
	    set ActiveAppList [list "DSOLID" "DEM"]
	} elseif {($DEMApplication eq "Yes")&& ($DSOLID =="No")&& ($FluidApplication =="No")} {
	    set ActiveAppList [list "DEM"]
	}
    }

    # WarnWinText "ActiveAppList:$ActiveAppList"

    # Debug/Release variable [0 => Debug, 1 => Release] Timers
    variable pflag
    set pflag 1

    # Get the element properties
    ::wkcf::GetElementProperties

    # Get properties data
    ::wkcf::GetPropertiesData

    # Convection diffusion
    if {$ConvectionDiffusionApplication eq "Yes"} {
	# Get the initial conditions
	::wkcf::GetInitialConditionProperties

	# Get the mapping between face heat flux and the property data
	::wkcf::GetPropertyDataFromFaceHeatFluxBC "ConvectionDiffusion"
    }

    # Get boundary condition properties
    ::wkcf::GetBoundaryConditionProperties

    if {$FluidApplication =="Yes"} {
	set AppId "Fluid"
	set cproperty "dv"
	# Free surface
	set cxpath "$AppId//c.AnalysisData//i.FreeSurface"
	set FreeSurface [::xmlutils::setXml $cxpath $cproperty]
	# WarnWinText "FreeSurface:$FreeSurface"
	if {$FreeSurface =="Yes"} {
	    # Get porous zones properties
	    ::wkcf::GetPorousZonesProperties $AppId
	}
    }

    if {$StructuralAnalysis =="Yes"} {
	# Get load properties
	::wkcf::GetLoadProperties
    }

    # Create the kratos global properties identifier
    ::wkcf::CreateKratosPropertiesIdentifier


    # To write the bat file
    variable wbatfile
    set wbatfile 0

    # Conditions to BC link
    variable ctbclink
    set ctbclink [dict create]

    # ::WinUtils::PrintArray dprops
}

proc ::wkcf::CreateKratosPropertiesIdentifier {} {
    # Create the kratos global properties identifier
    variable dprops; variable StructuralAnalysis
    variable FluidApplication

    # Check for used body forces
    set usebforce "No"
    if {$StructuralAnalysis =="Yes"} {
	set AppId "StructuralAnalysis"
	# List with all group that use body forces
	set dprops($AppId,AllBodyForceGroupId) [list]
	foreach cloadtid $dprops($AppId,AllLoadTypeId) {
	    if {$cloadtid =="BodyForce"} {
		# Check the group identifier
		if {([info exists dprops($AppId,Loads,$cloadtid,AllGroupId)]) && ([llength $dprops($AppId,Loads,$cloadtid,AllGroupId)]>0)} {
		    foreach cgroupid $dprops($AppId,Loads,$cloadtid,AllGroupId) {
		        if {$cgroupid ni $dprops($AppId,AllKEGroupId)} {
		            # Error => In this version all body force must belong to the same element group
		        } else {
		            set usebforce "Yes"
		            if {$cgroupid ni $dprops($AppId,AllBodyForceGroupId)} {
		                lappend dprops($AppId,AllBodyForceGroupId) $cgroupid
		            }
		        }
		    }
		} else {
		    # Error: First define some group identifier for this body force
		}
		break
	    }
	}

	# WarnWinText "usebforce:$usebforce"

	# Create the global kratos properties list
	set dprops($AppId,GKProps,AllPropertyId) $dprops($AppId,AllKPropertyId)
	# For all defined kratos elements
	foreach celemid $dprops($AppId,AllKElemId) {
	    # WarnWinText "celemid:$celemid"
	    # For all defined group identifier for this elements
	    foreach cgroupid $dprops($AppId,KElem,$celemid,AllGroupId) {
		# Get the group properties
		lassign $dprops($AppId,KElem,$celemid,$cgroupid,GProps) GiDEntity GiDElemType PropertyId KEKWord nDim
		# WarnWinText "cgroupid:$cgroupid PropertyId:$PropertyId"
		if {![info exists dprops($AppId,GKProps,$PropertyId,AddBF)]} {
		    set dprops($AppId,GKProps,$PropertyId,AddBF) "No"
		}
		if {$usebforce =="Yes"} {
		    # Check to add body force properties
		    if {$cgroupid in $dprops($AppId,AllBodyForceGroupId)} {
		        set dprops($AppId,GKProps,$PropertyId,AddBF) "Yes"
		    }
		}

		# Update the global kratos property identifier
		set GlobalPId [expr [lsearch $dprops($AppId,GKProps,AllPropertyId) $PropertyId]+1]
		set dprops($AppId,KElem,$celemid,$cgroupid,GlobalPId) $GlobalPId
	    }
	}
    }

    # For fluid application
    set AppId ""
    if {$FluidApplication eq "Yes"} {
	set AppId "Fluid"
    }

    # For convection-diffusion application
    variable ConvectionDiffusionApplication
    if {$ConvectionDiffusionApplication eq "Yes"} {
	set AppId "ConvectionDiffusion"
    }
    variable DEMApplication
    if {$DEMApplication eq "Yes"} {
	set AppId "DEM"
    }
    variable DSOLID
    if {$DSOLID eq "Yes"} {
  set AppId "DSOLID"
    }

    if {$AppId !=""} {
	# Create the global kratos properties list
	set dprops($AppId,GKProps,AllPropertyId) $dprops($AppId,AllKPropertyId)
	# For all defined kratos elements
	foreach celemid $dprops($AppId,AllKElemId) {
	    # For all defined group identifier for this elements
	    foreach cgroupid $dprops($AppId,KElem,$celemid,AllGroupId) {
		# Get the group properties
		# lassign $dprops($AppId,KElem,$celemid,$cgroupid,GProps) GiDEntity GiDElemType PropertyId KEKWord nDim

		# Update the global kratos property identifier
		set GlobalPId 0
		set dprops($AppId,KElem,$celemid,$cgroupid,GlobalPId) $GlobalPId
	    }
	}
    }
}

proc ::wkcf::GetLoadProperties {} {
    # Get all load properties
    variable dprops

    # Set the application root identifier
    set rootdataid "StructuralAnalysis"

    # Get all load properties
    set cxpath "$rootdataid//c.Loads"
    set clproplist_temp [::xmlutils::setXmlContainerIds $cxpath]

    #check that are containers class="Groups" , if are class="SameTemplateGroups" , add childs
    set clproplist [list]
    foreach cloadtid $clproplist_temp {
	set fullname "$rootdataid//c.Loads//c.[list $cloadtid]"
	set class [::xmlutils::setXml $fullname class]
	if { $class == "Groups" } {
	    lappend clproplist $cloadtid
	    set cloadtidpath "$cxpath//c.[list ${cloadtid}]"
	    lappend clproplist $cloadtidpath
	} elseif { $class == "SameTemplateGroups" } {
	    foreach childcloadtid [::xmlutils::setXmlContainerIds $fullname] {
		lappend clproplist $childcloadtid
		set cloadtidpath "$cxpath//c.[list ${cloadtid}]//c.[list ${childcloadtid}]"
		lappend clproplist $cloadtidpath
	    }
	}
    }

    # WarnWinText "clproplist:$clproplist"
    # Load type list
    set dprops($rootdataid,AllLoadTypeId) [list]
    foreach {cloadtid cxpath} $clproplist {
	# WarnWinText "cloadtid:$cloadtid"
	# Get the group identifier defined for this load type

	set cgrouplist [::xmlutils::setXmlContainerIds $cxpath]
	# WarnWinText "current cgrouplist:$cgrouplist"
	if {[llength $cgrouplist]} {
	    # Update load type identifier
	    lappend dprops($rootdataid,AllLoadTypeId) $cloadtid
	    # WarnWinText "inside cgrouplist:$cgrouplist"
	    foreach cgroupid $cgrouplist {
		# WarnWinText "cgroupid:$cgroupid"
		# Get the main properties
		set cgxpath "$cxpath//c.[list ${cgroupid}]"
		# WarnWinText "cgxpath:$cgxpath"
		set allmprop [::xmlutils::setXmlContainerPairs $cgxpath "" "dv"]
		# WarnWinText "allmprop:$allmprop"

		# Kratos load to group link
		# Group list
		if {![info exists dprops($rootdataid,Loads,$cloadtid,AllGroupId)]} {
		    set dprops($rootdataid,Loads,$cloadtid,AllGroupId) [list]
		}
		if {$cgroupid ni $dprops($rootdataid,Loads,$cloadtid,AllGroupId)} {
		    lappend dprops($rootdataid,Loads,$cloadtid,AllGroupId) $cgroupid
		}
		# Group properties
		if {![info exists dprops($rootdataid,Loads,$cloadtid,$cgroupid,GProps)]} {
		    set dprops($rootdataid,Loads,$cloadtid,$cgroupid,GProps) [list]
		}
		set dprops($rootdataid,Loads,$cloadtid,$cgroupid,GProps) $allmprop

	    }
	}
	# Reset the path
	set cxpath "$rootdataid//c.Loads"
    }
}

proc ::wkcf::GetInitialConditionProperties {} {
    # Get all initial condition properties
    variable dprops; variable ActiveAppList

    # For each active application
    foreach AppId $ActiveAppList {
	# Get the application root identifier
	set rootdataid $AppId
	# Get all defined initial condition groups
	set cxpath "$rootdataid//c.InitialConditions"
	set cicproplist [::xmlutils::setXmlContainerIds $cxpath]
	# WarnWinText "cicproplist:$cicproplist"
	# Initial condition type list
	set dprops($AppId,AllICTypeId) [list]
	foreach cictid $cicproplist {
	    # WarnWinText "cictid:$cictid"
	    # Get the group identifier defined for this condition
	    set cxpath "${cxpath}//c.[list ${cictid}]"
	    set cicgrouplist [::xmlutils::setXmlContainerIds $cxpath]
	    # WarnWinText "cicgrouplist:$cicgrouplist"
	    if {[llength $cicgrouplist]} {
		# Update initial condition type identifier
		lappend dprops($AppId,AllICTypeId) $cictid
		# WarnWinText "inside cicgrouplist:$cicgrouplist"
		foreach cgroupid $cicgrouplist {
		    set proplist [list]
		    switch -exact -- $cictid {
		        "InitialTemperature" {
		            # Get properties
		            foreach citem [list "InitialTemperature"] {
		                # set xpath
		                set pcxpath "$cxpath//c.[list ${cgroupid}]//c.MainProperties//i.[list ${citem}]"
		                set cproperty "dv"
		                set CValue [::xmlutils::setXml $pcxpath $cproperty]
		                lappend proplist $CValue
		            }
		        }
		    }

		    # WarnWinText "proplist:$proplist"
		    # Kratos IC to group link
		    # Group list
		    if {![info exists dprops($AppId,IC,$cictid,AllGroupId)]} {
		        set dprops($AppId,IC,$cictid,AllGroupId) [list]
		    }
		    if {$cgroupid ni $dprops($AppId,IC,$cictid,AllGroupId)} {
		        lappend dprops($AppId,IC,$cictid,AllGroupId) $cgroupid
		    }
		    # Group properties
		    if {![info exists dprops($AppId,IC,$cictid,$cgroupid,GProps)]} {
		        set dprops($AppId,IC,$cictid,$cgroupid,GProps) [list]
		    }
		    set dprops($AppId,IC,$cictid,$cgroupid,GProps) $proplist
		}
	    }
	    # Reset the path
	    set cxpath "$rootdataid//c.InitialConditions"
	}
    }
}

proc ::wkcf::GetBoundaryConditionProperties {} {
    # Get all boundary condition properties
    variable dprops; variable ActiveAppList

    # For each active application
    foreach AppId $ActiveAppList {
	# Get the application root identifier
	set rootdataid $AppId
	    # Get all defined condition groups
	set cxpath "$rootdataid//c.Conditions"
	set cbcproplist [::xmlutils::setXmlContainerIds $cxpath]
	# WarnWinText "cbcproplist:$cbcproplist"
	# Boundary condition type list
	set dprops($AppId,AllBCTypeId) [list]
	foreach cbctid $cbcproplist {
	    # WarnWinText "cbctid:$cbctid"
	    # Get the group identifier defined for this condition
	    set cxpath "${cxpath}//c.[list ${cbctid}]"
	    set cbcgrouplist [::xmlutils::setXmlContainerIds $cxpath]
	    # WarnWinText "cbcgrouplist:$cbcgrouplist"
	    if {[llength $cbcgrouplist]} {
		# Update boundary condition type identifier
		lappend dprops($AppId,AllBCTypeId) $cbctid
		# WarnWinText "inside cbcgrouplist:$cbcgrouplist"
		foreach cgroupid $cbcgrouplist {
		    set proplist [list]
		    switch -exact -- $cbctid {
		        "Displacements" - "Rotations" - "InletVelocity" - "No-Slip" {
		            # Get activation properties
		            foreach ca [list Ax Ay Az] cv [list Vx Vy Vz] {
		                # Activation
		                set acxpath "$cxpath//c.[list ${cgroupid}]//c.Activation//i.[list ${ca}]"
		                # WarnWinText "Activation :cxpath:$cxpath"
		                set cproperty "dv"
		                set CActive [::xmlutils::setXml $acxpath $cproperty]
		                # Values
		                set vcxpath "$cxpath//c.[list ${cgroupid}]//c.Values//i.[list ${cv}]"
		                # WarnWinText "Values :cxpath:$cxpath"
		                set cproperty "dv"
		                set CValue [::xmlutils::setXml $vcxpath $cproperty]
		                lappend proplist $CActive $CValue
		            }
		        }
		        "OutletPressure" {
		            # Get properties
		            foreach citem [list "FixPressure" "PressureValue"] {
		                # set xpath
		                set pcxpath "$cxpath//c.[list ${cgroupid}]//c.MainProperties//i.[list ${citem}]"
		                set cproperty "dv"
		                set CValue [::xmlutils::setXml $pcxpath $cproperty]
		                lappend proplist $CValue
		            }
		        }
		        "Flag-Variable" {
		            # Get properties
		            foreach citem [list "Flag"] {
		                # set xpath
		                set pcxpath "$cxpath//c.[list ${cgroupid}]//c.MainProperties//i.[list ${citem}]"
		                set cproperty "dv"
		                set CValue [::xmlutils::setXml $pcxpath $cproperty]
		                lappend proplist $CValue
		            }
		        }
		        "Is-Slip" {
		            # Get properties
		            foreach citem [list "Activate" "ConstantValue"] {
		                # set xpath
		                set pcxpath "$cxpath//c.[list ${cgroupid}]//c.MainProperties//i.[list ${citem}]"
		                set cproperty "dv"
		                set CValue [::xmlutils::setXml $pcxpath $cproperty]
		                lappend proplist $CValue
		            }
		        }
		        "WallLaw" {
		            # Get properties
		            foreach citem [list "Activate" "ConstantValue"] {
		                # set xpath
		                set pcxpath "$cxpath//c.[list ${cgroupid}]//c.MainProperties//i.[list ${citem}]"
		                set cproperty "dv"
		                set CValue [::xmlutils::setXml $pcxpath $cproperty]
		                lappend proplist $CValue
		            }
		        }
		        "Distance" {
		            # Get properties
		            foreach citem [list "DistanceValue"] {
		                # set xpath
		                set pcxpath "$cxpath//c.[list ${cgroupid}]//c.MainProperties//i.[list ${citem}]"
		                set cproperty "dv"
		                set CValue [::xmlutils::setXml $pcxpath $cproperty]
		                lappend proplist $CValue
		            }
		        }
		        "PFEMWall--" {
		        # Commented by J. Garate on 17/12/2012
		            # Get properties
		            set citem "LinearVelocity"
		            # set xpath
		            foreach cv [list LVx LVy LVz] {
		                set pcxpath "$cxpath//c.[list ${cgroupid}]//c.[list $citem]//i.[list ${cv}]"
		                set cproperty "dv"
		                # msg "pcxpath : $pcxpath"
		                # msg "cproperty : $cproperty"
		                set CValue [::xmlutils::setXml $pcxpath $cproperty]
		                # msg "CValue : $CValue"
		                lappend proplist $CValue
		            }

		            set citem "AngularVelocity"
		            # set xpath
		            foreach cv [list AVx AVy AVz] {
		                set pcxpath "$cxpath//c.[list ${cgroupid}]//c.[list $citem]//i.[list ${cv}]"
		                set cproperty "dv"
		                # msg "pcxpath : $pcxpath"
		                # msg "cproperty : $cproperty"
		                set CValue [::xmlutils::setXml $pcxpath $cproperty]
		                # msg "CValue : $CValue"
		                lappend proplist $CValue
		            }

		            set citem "RotationCenter"
		            foreach cv [list Gx Gy Gz] {
		                set pcxpath "$cxpath//c.[list ${cgroupid}]//c.[list $citem]//i.[list ${cv}]"
		                set cproperty "dv"
		                set CValue [::xmlutils::setXml $pcxpath $cproperty]
		                lappend proplist $CValue
		            }
		        }
		        "PFEMFluidInlet" {
		            # Get properties
		            foreach ca [list Ax Ay Az] cv [list Vx Vy Vz] {
		                # Activation
		                set acxpath "$cxpath//c.[list ${cgroupid}]//c.Activation//i.[list ${ca}]"
		                # WarnWinText "Activation :acxpath:$acxpath"
		                set cproperty "dv"
		                set CActive [::xmlutils::setXml $acxpath $cproperty]
		                # Values
		                set vcxpath "$cxpath//c.[list ${cgroupid}]//c.Values//i.[list ${cv}]"
		                # WarnWinText "Values :vcxpath:$vcxpath"
		                set cproperty "dv"
		                set CValue [::xmlutils::setXml $vcxpath $cproperty]
		                lappend proplist $CActive $CValue
		                # msg "$CActive $CValue"
		            }
		        }
		        "PrescribedTemperature" {
		            # Get properties
		            set cproperty "dv"
		            foreach ca [list APTemperature] cv [list VPTemperature] {
		                # Activation
		                set acxpath "$cxpath//c.[list ${cgroupid}]//c.Activation//i.[list ${ca}]"
		                # WarnWinText "Activation :acxpath:$acxpath"
		                set CActive [::xmlutils::setXml $acxpath $cproperty]
		                # Values
		                set vcxpath "$cxpath//c.[list ${cgroupid}]//c.Values//i.[list ${cv}]"
		                # WarnWinText "Values :vcxpath:$vcxpath"
		                set CValue [::xmlutils::setXml $vcxpath $cproperty]
		                lappend proplist $CActive $CValue
		                # msg "$CActive $CValue"
		            }
		        }
		        "HeatFlux" {
		            # Get properties
		            foreach citem [list "VHeatFlux"] {
		                # set xpath
		                set pcxpath "$cxpath//c.[list ${cgroupid}]//c.MainProperties//i.[list ${citem}]"
		                set cproperty "dv"
		                set CValue [::xmlutils::setXml $pcxpath $cproperty]
		                lappend proplist $CValue
		            }
		        }
		        "FaceHeatFlux" {
		            # Get properties
		            foreach citem [list "VFaceHeatFlux"] {
		                # set xpath
		                set pcxpath "$cxpath//c.[list ${cgroupid}]//c.MainProperties//i.[list ${citem}]"
		                set cproperty "dv"
		                set CValue [::xmlutils::setXml $pcxpath $cproperty]
		                lappend proplist $CValue
		            }
		        }
		        "PFEMFixedWall" {
		            # PFEM fixed wall boundary condition
		            # Get properties
		            foreach citem [list "Dx" "Dy" "Dz"] {
		                # set xpath
		                set pcxpath "$cxpath//c.[list ${cgroupid}]//c.Displacement//i.[list ${citem}]"
		                set cproperty "dv"
		                set CValue [::xmlutils::setXml $pcxpath $cproperty]
		                lappend proplist $CValue
		            }
		        }
		    }
		    # WarnWinText "proplist:$proplist"
		    # Kratos BC to group link
		    # Group list
		    if {![info exists dprops($AppId,BC,$cbctid,AllGroupId)]} {
		        set dprops($AppId,BC,$cbctid,AllGroupId) [list]
		    }
		    if {$cgroupid ni $dprops($AppId,BC,$cbctid,AllGroupId)} {
		        lappend dprops($AppId,BC,$cbctid,AllGroupId) $cgroupid
		    }
		    # Group properties
		    if {![info exists dprops($AppId,BC,$cbctid,$cgroupid,GProps)]} {
		        set dprops($AppId,BC,$cbctid,$cgroupid,GProps) [list]
		    }
		    set dprops($AppId,BC,$cbctid,$cgroupid,GProps) $proplist
		}
	    }
	    # Reset the path
	    set cxpath "$rootdataid//c.Conditions"
	}
    }
}

proc ::wkcf::GetPorousZonesProperties {AppId} {
    # ABSTRACT: Get all porous zones properties
    variable dprops

    # Get the application root identifier
    set rootdataid $AppId

    # Get use Ergun equation value
    set cproperty "dv"
    set contid "PorousZones"
    set cxpath "$rootdataid//c.SolutionStrategy//c.[list ${contid}]//i.UseErgunEquation"
    set UseErgunEquation [::xmlutils::setXml $cxpath $cproperty]

    # wa "UseErgunEquation:$UseErgunEquation"
    # Get all defined porous zone groups
    set cxpath "$rootdataid//c.SolutionStrategy//c.[list ${contid}]"
    set cbcproplist [::xmlutils::setXmlContainerIds $cxpath]
    # wa "cbcproplist:$cbcproplist"
    # Porous zones type list
    set dprops($AppId,AllPorousZonesTypeId) [list]
    foreach cbctid $cbcproplist {
	# wa "cbctid:$cbctid"
	set cproplist [list]
	if {(($cbctid eq "ErgunEquationNo") && ($UseErgunEquation eq "No"))} {
	    set cproplist [list "PorosityValue" "LinearDarcyCoefficient" "NonLinearDarcyCoefficient"]
	} elseif {(($cbctid eq "ErgunEquationYes") && ($UseErgunEquation eq "Yes"))} {
	    set cproplist [list "PorosityValue" "DiameterValue"]
	}
	if {![llength $cproplist]} {
	    continue
	}
	# Get the group identifier defined for this condition
	set cxpath "${cxpath}//c.[list ${cbctid}]"
	set cbcgrouplist [::xmlutils::setXmlContainerIds $cxpath]
	# wa "cbcgrouplist:$cbcgrouplist cxpath:$cxpath"
	if {[llength $cbcgrouplist]} {
	    # Update zones type identifier
	    lappend dprops($AppId,AllPorousZonesTypeId) $cbctid
	    # wa "inside cbcgrouplist:$cbcgrouplist cxpath:$cxpath"
	    foreach cgroupid $cbcgrouplist {
		set proplist [list]
		# Get properties
		foreach cprop $cproplist {
		    set cvxpath "$cxpath//c.[list ${cgroupid}]//c.MainProperties//i.[list ${cprop}]"
		    # wa "cvxpath:$cvxpath cprop:$cprop"
		    set cproperty "dv"
		    set cvalue [::xmlutils::setXml $cvxpath $cproperty]
		    lappend proplist $cvalue
		}
		# wa "proplist:$proplist"
		# Kratos porous zones to group link
		# Group list
		if {![info exists dprops($AppId,$contid,$cbctid,AllGroupId)]} {
		    set dprops($AppId,$contid,$cbctid,AllGroupId) [list]
		}
		if {$cgroupid ni $dprops($AppId,$contid,$cbctid,AllGroupId)} {
		    lappend dprops($AppId,$contid,$cbctid,AllGroupId) $cgroupid
		}
		# Group properties
		if {![info exists dprops($AppId,$contid,$cbctid,$cgroupid,GProps)]} {
		    set dprops($AppId,$contid,$cbctid,$cgroupid,GProps) [list]
		}
		set dprops($AppId,$contid,$cbctid,$cgroupid,GProps) $proplist

	    }
	}
	# Reset the path
	set cxpath "//c.SolutionStrategy//c.[list ${contid}]"
    }
}

proc ::wkcf::GetApplicationRootId {} {
    # Get the application root identifier
    variable FluidApplication;  variable StructuralAnalysis
    variable ConvectionDiffusionApplication; variable DEMApplication

    # Select the current active application
    if {$StructuralAnalysis =="Yes"} {
	set rootdataid "StructuralAnalysis"
    } elseif {$FluidApplication =="Yes"} {
	    set rootdataid "Fluid"
    } elseif {$ConvectionDiffusionApplication eq "Yes"} {
	set rootdataid "ConvectionDiffusion"
    } elseif {$DEMApplication eq "Yes" && $FluidApplication =="No"} {
	set rootdataid "DEM"
    }

    return $rootdataid
}

proc ::wkcf::GetPropertiesData {} {
    # Process all properties
    variable dprops; variable ActiveAppList
    variable ndime

    # Xpath for element constitutive laws
    set cexpath "ElementCLaws"

    # For each active application
    foreach AppId $ActiveAppList {
	# Get the application root identifier
	set rootdataid $AppId
	# Get the properties identifier
	set cxpath "$rootdataid//c.Properties"
	set cproplist [::xmlutils::setXmlContainerIds $cxpath]
	# WarnWinText "cproplist:$cproplist"
	# All material list
	set dprops($AppId,AllMatId) [list]
	# Get the properties
	foreach propid $cproplist {
	    # Material identifier
	    set mxpath "$cxpath//c.[list ${propid}]//c.MainProperties//i.Material"
	    set cproperty "dv"
	    set MatId [::xmlutils::setXml $mxpath $cproperty]
	    # WarnWinText "MatId:$MatId"
	    set dprops($AppId,Property,$propid,MatId) "$MatId"
	    # Get the material properties
	    if {$MatId ni $dprops($AppId,AllMatId)} {
		lappend dprops($AppId,AllMatId) $MatId
	    }

	    # Property type => Base element type
	    set ptypexpath "$cxpath//c.[list ${propid}]//c.MainProperties//i.ElemType"
	    set cproperty "dv"
	    set ptype [::xmlutils::setXml $ptypexpath $cproperty]
	    # WarnWinText "ptype:$ptype"
	    set dprops($AppId,Property,$propid,BaseElemType) $ptype

	    # Material model
	    set xpath "$cxpath//c.[list ${propid}]//c.MainProperties//i.MatModel"
	    set cproperty "dv"
	    set MatModel [::xmlutils::setXml $xpath $cproperty]
	    # WarnWinText "MatModel:$MatModel"
	    set dprops($AppId,Property,$propid,MatModel) $MatModel

	    if {$AppId ne "StructuralAnalysis"} {
		continue
	    }

	    # Only for structural analysis
	    # For cross section properties
	    # Get the property list
	    set pid "propertylist"
	    set id "CSProperty"
	    set CSProperty [split [::xmlutils::getKKWord $cexpath "$id" "$pid"] ","]
	    # wa "CSProperty:$CSProperty"
	    foreach cspropid $CSProperty {
		# Get the current value
		set txpath "$cxpath//c.[list ${propid}]//c.MainProperties//i.[list ${cspropid}]"
		set cproperty "dv"
		set cvalue [::xmlutils::setXml $txpath $cproperty]
		# wa "cspropid:$cspropid cvalue:$cvalue"
		set dprops($AppId,Property,$propid,${cspropid}) $cvalue
	    }

	    # Get the list of element that have this propertyid assigned
	    set ElementIdList [::wkcf::GetElementIdFromPropertyId $AppId $propid]
	    # Get the first element of the list
	    set etype [lindex $ElementIdList 0]
	    # wa "ElementIdList:$ElementIdList"

	    # Set fluency and behavior variables
	    set dprops($AppId,Material,$MatId,UseFluency) "No"
	    set dprops($AppId,Material,$MatId,Fluency) ""
	    set dprops($AppId,Material,$MatId,UseBehavior) "No"
	    set dprops($AppId,Material,$MatId,Behavior) ""

	    set dprops($AppId,Material,$MatId,UseFlowRule) "No"
	    set dprops($AppId,Material,$MatId,FlowRule) ""
	    set dprops($AppId,Material,$MatId,UseYieldCriterion) "No"
	    set dprops($AppId,Material,$MatId,YieldCriterion) ""
	    set dprops($AppId,Material,$MatId,UseHardeningLaw) "No"
	    set dprops($AppId,Material,$MatId,HardeningLaw) ""

	    # Get material properties
	    switch -exact -- $MatModel {
		"Elastic-Isotropic" {
		    if {$ndime eq "2D"} {
		        # 2D case
		        if {$ptype=="Solid"} {
		            # I need to check if the element is PlaneStrain2D or PlaneStress3D or
		            # Axisymmetric2D
		            if {($etype eq "PlaneStrain2D") || ($etype eq "TotalLagrangian2DPlaneStrain") || ($etype eq "UpdatedLagrangian2DPlaneStrain") || ($etype eq "SpatialLagrangian2DPlaneStrain")} {
		                set cptype "LinearElasticPlaneStrain2D"
		                # Get the material properties
		                ::wkcf::GetMaterialProperties $AppId $propid $MatId $cptype $MatModel

		                # Get the cross section properties
		                ::wkcf::GetCrossSectionProperties $AppId $propid $ptype

		            } elseif {($etype eq "PlaneStress2D") || ($etype eq "TotalLagrangian2DPlaneStress") || ($etype eq "UpdatedLagrangian2DPlaneStress") || ($etype eq "SpatialLagrangian2DPlaneStress")} {
		                set cptype "LinearElasticPlaneStress2D"
		                # Get the material properties
		                ::wkcf::GetMaterialProperties $AppId $propid $MatId $cptype $MatModel

		                # Get the cross section properties
		                ::wkcf::GetCrossSectionProperties $AppId $propid $ptype

		            } elseif {($etype eq "Axisymmetric2D") || ($etype eq "TotalLagrangian2DAxisymmetric") || ($etype eq "UpdatedLagrangian2DAxisymmetric") || ($etype eq "SpatialLagrangian2DAxisymmetric")} {
		                set cptype "LinearElasticAxisymmetric2D"
		                # Get the material properties
		                ::wkcf::GetMaterialProperties $AppId $propid $MatId $cptype $MatModel

		                # Get the cross section properties
		                ::wkcf::GetCrossSectionProperties $AppId $propid $ptype


		            } else { # added here because I need another solution
		                # Confirmar la veracidad!
		                set cptype "LinearElasticPlaneStrain2D"
		                # Get the material properties
		                ::wkcf::GetMaterialProperties $AppId $propid $MatId $cptype $MatModel

		                # Get the cross section properties
		                ::wkcf::GetCrossSectionProperties $AppId $propid $ptype

		            }
		        } elseif {$ptype=="BeamElement"} {
		            set cptype "LinearElasticPlaneStrain2D"
		            # Get the material properties
		            ::wkcf::GetMaterialProperties $AppId $propid $MatId $cptype $MatModel

		            # Get the cross section properties
		            ::wkcf::GetCrossSectionProperties $AppId $propid $ptype
		        }
		    } elseif {$ndime =="3D"} {
		        # 3D case
		        set ltypelist [list "Solid" "Shell" "Membrane" "Beam" "Truss" "EBST"]
		        # Check that this ptype is in the list
		        if {$ptype in $ltypelist} {
		            if {($ptype eq "Shell") || ($ptype eq "Membrane") || ($ptype eq "EBST")} {
		                set cptype "LinearElasticPlaneStress2D"
		            } elseif {($ptype eq "Beam") || ($ptype eq "Truss")} {
		                set cptype "Isotropic2D"
		            } else { #SolidElement3D, TotalLagrangian3D, UpdatedLagrangian3D, SpatialLagrangian3D
		                set cptype "LinearElastic3D"
		            }

		            # Get the material properties
		            ::wkcf::GetMaterialProperties $AppId $propid $MatId $cptype $MatModel

		            # Get the cross section properties
		            ::wkcf::GetCrossSectionProperties $AppId $propid $ptype
		        }
		    }
		}
		"HyperElastic-Isotropic" {
		    if {$ndime eq "2D"} {
		        # 2D case
		        if {$ptype=="Solid"} {
		            # I need to check if the element is PlaneStrain2D or PlaneStress3D or
		            # Axisymmetric2D
		            if {($etype eq "PlaneStrain2D") || ($etype eq "TotalLagrangian2DPlaneStrain") || ($etype eq "UpdatedLagrangian2DPlaneStrain") || ($etype eq "SpatialLagrangian2DPlaneStrain")} {
		                set cptype "HyperElasticPlaneStrain2D"

		                # Get the material properties
		                ::wkcf::GetMaterialProperties $AppId $propid $MatId $cptype $MatModel

		                # Get the cross section properties
		                ::wkcf::GetCrossSectionProperties $AppId $propid $ptype

		            } elseif {($etype eq "PlaneStress2D") || ($etype eq "TotalLagrangian2DPlaneStress") || ($etype eq "UpdatedLagrangian2DPlaneStress") || ($etype eq "SpatialLagrangian2DPlaneStress")} {
		                set cptype "HyperElasticPlaneStress2D"

		                # Get the material properties
		                ::wkcf::GetMaterialProperties $AppId $propid $MatId $cptype $MatModel

		                # Get the cross section properties
		                ::wkcf::GetCrossSectionProperties $AppId $propid $ptype


		            } elseif {($etype eq "Axisymmetric2D") || ($etype eq "TotalLagrangian2DAxisymmetric") || ($etype eq "UpdatedLagrangian2DAxisymmetric") || ($etype eq "SpatialLagrangian2DAxisymmetric")} {
		                set cptype "HyperElasticAxisymmetric2D"
		                # Get the material properties
		                ::wkcf::GetMaterialProperties $AppId $propid $MatId $cptype $MatModel

		                # Get the cross section properties
		                ::wkcf::GetCrossSectionProperties $AppId $propid $ptype


		            } else { # added here because I need another solution
		                # Confirmar la veracidad!
		                set cptype "HyperElasticPlaneStrain2D"
		                # Get the material properties
		                ::wkcf::GetMaterialProperties $AppId $propid $MatId $cptype $MatModel

		                # Get the cross section properties
		                ::wkcf::GetCrossSectionProperties $AppId $propid $ptype

		            }
		        } elseif {$ptype=="Beam"} {
		            set cptype "Isotropic2D"
		            # Get the material properties
		            ::wkcf::GetMaterialProperties $AppId $propid $MatId $cptype $MatModel

		            # Get the cross section properties
		            ::wkcf::GetCrossSectionProperties $AppId $propid $ptype
		        }
		    } elseif {$ndime =="3D"} {
		        # 3D case
		        set ltypelist [list "Solid" "Shell" "Membrane" "Beam" "Truss" "EBST"]
		        # Check that this ptype is in the list
		        if {$ptype in $ltypelist} {
		            if {($ptype eq "Shell") || ($ptype eq "Membrane") || ($ptype eq "EBST")} {
		                set cptype "LinearElasticPlaneStress2D"
		            } elseif {($ptype eq "Beam") || ($ptype eq "Truss")} {
		                set cptype "Isotropic2D"
		            } else { #SolidElement3D, TotalLagrangian3D, UpdatedLagrangian3D, SpatialLagrangian3D
		                set cptype "HyperElastic3D"
		            }

		            # Get the material properties
		            ::wkcf::GetMaterialProperties $AppId $propid $MatId $cptype $MatModel

		            # Get the cross section properties
		            ::wkcf::GetCrossSectionProperties $AppId $propid $ptype
		        }
		    }
		}
		"HyperElastic-Plastic" {
		    if {$ndime eq "2D"} {
		        # 2D case
		        if {$ptype=="Solid"} {
		            # I need to check if the element is PlaneStrain2D or PlaneStress3D or
		            # Axisymmetric2D
		            if {($etype eq "PlaneStrain2D") || ($etype eq "TotalLagrangian2DPlaneStrain") || ($etype eq "UpdatedLagrangian2DPlaneStrain") || ($etype eq "SpatialLagrangian2DPlaneStrain")} {
		                set cptype "HyperElasticPlasticPlaneStrain2D"
		                # Get the material properties
		                ::wkcf::GetMaterialProperties $AppId $propid $MatId $cptype $MatModel

		                # Get the material behavior and fluency properties
		                ::wkcf::GetPlasticityProperties $AppId $MatId $MatModel $ptype $cptype

		                # Get the cross section properties
		                ::wkcf::GetCrossSectionProperties $AppId $propid $ptype

		            } elseif {($etype eq "PlaneStress2D") || ($etype eq "TotalLagrangian2DPlaneStress") || ($etype eq "UpdatedLagrangian2DPlaneStress") || ($etype eq "SpatialLagrangian2DPlaneStress")} {
		                set cptype "HyperElasticPlasticPlaneStress2D"
		                # Get the material properties
		                ::wkcf::GetMaterialProperties $AppId $propid $MatId $cptype $MatModel

		                # Get the material behavior and fluency properties
		                ::wkcf::GetPlasticityProperties $AppId $MatId $MatModel $ptype $cptype

		                # Get the cross section properties
		                ::wkcf::GetCrossSectionProperties $AppId $propid $ptype

		            } elseif {($etype eq "Axisymmetric2D") || ($etype eq "TotalLagrangian2DAxisymmetric") || ($etype eq "UpdatedLagrangian2DAxisymmetric") || ($etype eq "SpatialLagrangian2DAxisymmetric")} {
		                set cptype "HyperElasticPlasticAxisymmetric2D"
		                # Get the material properties
		                ::wkcf::GetMaterialProperties $AppId $propid $MatId $cptype $MatModel

		                # Get the material behavior and fluency properties
		                ::wkcf::GetPlasticityProperties $AppId $MatId $MatModel $ptype $cptype

		                # Get the cross section properties
		                ::wkcf::GetCrossSectionProperties $AppId $propid $ptype

		            } else { # added here because I need another solution
		                # Confirmar la veracidad!
		                set cptype "HyperElasticPlasticPlaneStrain2D"
		                # Get the material properties
		                ::wkcf::GetMaterialProperties $AppId $propid $MatId $cptype $MatModel

		                # Get the material behavior and fluency properties
		                ::wkcf::GetPlasticityProperties $AppId $MatId $MatModel $ptype $cptype

		                # Get the cross section properties
		                ::wkcf::GetCrossSectionProperties $AppId $propid $ptype
		            }
		        } elseif {$ptype=="Beam"} {
		            set cptype "Isotropic2D"
		            # Get the material properties
		            ::wkcf::GetMaterialProperties $AppId $propid $MatId $cptype $MatModel

		            # Get the cross section properties
		            ::wkcf::GetCrossSectionProperties $AppId $propid $ptype
		        }
		    } elseif {$ndime =="3D"} {
		        # 3D case
		        set ltypelist [list "Solid" "Shell" "Membrane" "Beam" "Truss" "EBST"]
		        # Check that this ptype is in the list
		        if {$ptype in $ltypelist} {
		            if {($ptype eq "Shell") || ($ptype eq "Membrane") || ($ptype eq "EBST")} {
		                set cptype "Isotropic2D"
		            } elseif {($ptype eq "Beam") || ($ptype eq "Truss")} {
		                set cptype "Isotropic2D"
		            } else { #SolidElement3D, TotalLagrangian3D, UpdatedLagrangian3D, SpatialLagrangian3D
		                set cptype "HyperElasticPlastic3D"
		            }

		            # Get the material properties
		            ::wkcf::GetMaterialProperties $AppId $propid $MatId $cptype $MatModel

		            # Get the material behavior and fluency properties
		            ::wkcf::GetPlasticityProperties $AppId $MatId $MatModel $ptype $cptype
		            # Get the cross section properties
		            ::wkcf::GetCrossSectionProperties $AppId $propid $ptype
		        }
		    }
		}
		"Elastic-Orthotropic" {
		    msgS "Elastic-Orthotropic is not implemented yet. Contact the Kratos team."
		}
		"Elasto-Plastic" {
		    if {($ptype=="PlaneStrain") && ($ndime =="2D")} {
		        set cptype "Plasticity2D"
		        # Get the material properties
		        ::wkcf::GetMaterialProperties $AppId $propid $MatId $cptype $MatModel

		        # Get the material behavior and fluency properties
		        ::wkcf::GetBehaviorFluencyProperties $AppId $MatId $MatModel $ptype $cptype

		        # Get the cross section properties
		        ::wkcf::GetCrossSectionProperties $AppId $propid $ptype

		    } elseif {($ptype=="PlaneStress") && ($ndime =="2D")} {
		        set cptype "Plasticity2D"
		        # Get the material properties
		        ::wkcf::GetMaterialProperties $AppId $propid $MatId $cptype $MatModel

		        # Get the material behavior and fluency properties
		        ::wkcf::GetBehaviorFluencyProperties $AppId $MatId $MatModel $ptype $cptype

		        # Get the cross section properties
		        ::wkcf::GetCrossSectionProperties $AppId $propid $ptype
		    } elseif {($ptype=="Solid") && ($ndime =="2D")} {
		        set cptype "Plasticity2D"
		        # Get the material properties
		        ::wkcf::GetMaterialProperties $AppId $propid $MatId $cptype $MatModel

		        # Get the material behavior and fluency properties
		        ::wkcf::GetBehaviorFluencyProperties $AppId $MatId $MatModel $ptype $cptype

		        # Get the cross section properties
		        ::wkcf::GetCrossSectionProperties $AppId $propid $ptype

		    } elseif {($ptype=="Solid") && ($ndime =="3D")} {
		        set cptype "Plasticity3D"
		        # Get the material properties
		        ::wkcf::GetMaterialProperties $AppId $propid $MatId $cptype $MatModel

		        # Get the material behavior and fluency properties
		        ::wkcf::GetBehaviorFluencyProperties $AppId $MatId $MatModel $ptype $cptype

		        # Get the cross section properties
		        ::wkcf::GetCrossSectionProperties $AppId $propid $ptype
		    }
		}
		"Damage" {
		    if {($ptype=="PlaneStrain") && ($ndime =="2D")} {
		        set cptype "IsotropicDamage"
		        # Get the material properties
		        ::wkcf::GetMaterialProperties $AppId $propid $MatId $cptype $MatModel

		        # Get the material behavior and fluency properties
		        ::wkcf::GetBehaviorFluencyProperties $AppId $MatId $MatModel $ptype $cptype

		        # Get the cross section properties
		        ::wkcf::GetCrossSectionProperties $AppId $propid $ptype

		    } elseif {($ptype=="PlaneStress") && ($ndime =="2D")} {
		        set cptype "IsotropicDamage"
		        # Get the material properties
		        ::wkcf::GetMaterialProperties $AppId $propid $MatId $cptype $MatModel

		        # Get the material behavior and fluency properties
		        ::wkcf::GetBehaviorFluencyProperties $AppId $MatId $MatModel $ptype $cptype

		        # Get the cross section properties
		        ::wkcf::GetCrossSectionProperties $AppId $propid $ptype

		    } elseif {($ptype=="Solid") && ($ndime =="3D")} {
		        set cptype "IsotropicDamage3D"
		        # Get the material properties
		        ::wkcf::GetMaterialProperties $AppId $propid $MatId $cptype $MatModel

		        # Get the material behavior and fluency properties
		        ::wkcf::GetBehaviorFluencyProperties $AppId $MatId $MatModel $ptype $cptype

		        # Get the cross section properties
		        ::wkcf::GetCrossSectionProperties $AppId $propid $ptype

		    } elseif {($ptype=="Solid") && ($ndime =="2D")} {
		        set cptype "IsotropicDamage2D"
		        # Get the material properties
		        ::wkcf::GetMaterialProperties $AppId $propid $MatId $cptype $MatModel

		        # Get the material behavior and fluency properties
		        ::wkcf::GetBehaviorFluencyProperties $AppId $MatId $MatModel $ptype $cptype

		        # Get the cross section properties
		        ::wkcf::GetCrossSectionProperties $AppId $propid $ptype
		    }
		}
	    }
	}
    }
}

proc ::wkcf::GetMaterialProperties {AppId propid MatId ptype CMatModel} {
    # Get material properties
    variable dprops
    variable ndime

    # wa "MaterialProperties =>AppId:$AppId propid:$propid MatId:$MatId ptype:$ptype CMatModel:$CMatModel"
    # Xpath for constitutive laws
    set clxpath "CLawProperties"
    # Xpath for materials
    set matxpath "Materials"

    # Get all material properties
    set mpxpath "[::KMat::findMaterialParent $MatId]//m.[list ${MatId}]"
    # wa "mpxpath:$mpxpath"

    # Set the kratos model keyword base xpath
    set kmxpath "Applications//$AppId"

    # Material model
    set MatModel [::xmlutils::getKKWord $clxpath $ptype "matm"]
    set dprops($AppId,Material,$MatId,MatModel) "$MatModel"
    # Get the used material properties
    set mprops [split [::xmlutils::getKKWord $clxpath $ptype "mprops"] ,]
    # Get xpath values
    set mxpath [split [::xmlutils::getKKWord $clxpath $ptype "mxpath"] ,]
    # wa "MatModel:$MatModel mprops:$mprops mxpath:$mxpath"
    set matplist [list]
    foreach pid $mprops xpath $mxpath {
	# Get the kratos key word
	set kkword [::xmlutils::getKKWord $matxpath $pid "kkword"]
	# wa "pid:$pid kkword:$kkword"
	# Get the current value for this properties
	# WarnWinText "xpath:$mpxpath//$xpath//p.[list $pid]"
	set cvalue [lindex [::KMat::getMaterialProperties "p" "$mpxpath//$xpath//p.[list $pid]"] 0 1]
	# wa "cvalue:$cvalue"
	if {($kkword !="") && ($cvalue !="")} {
	    lappend matplist [list $kkword $cvalue]
	}
    }

    # Add others properties for specific constitutive models
    if {$CMatModel == "HyperElastic-Plastic"} {
	# Get the yield function properties
	set yfid "YieldCriteria"
	# Get the yield criteria
	set myieldcriteria [::xmlutils::getKKWord $clxpath $ptype "myieldcriterion"]
	# wa "clxpath:$clxpath ptype:$ptype"
	# Get the yield criteria xpath values
	set mycxpath [::xmlutils::getKKWord $clxpath $ptype "mycxpath"]
	# Get the current yield criteria
	set cycvalue [lindex [::KMat::getMaterialProperties "p" "$mpxpath//$mycxpath//p.[list $myieldcriteria]"] 0 1]
	# WarnWinText "myieldcriteria:$myieldcriteria mycxpath:$mycxpath cycvalue:$cycvalue"
	# Get the yield function options
	set cyf ""
	set yfivalues [split [::xmlutils::getKKWord "$clxpath//$yfid" "AvailableYieldCriteria" "ycivalues"] ,]
	# wa "yfivalues:$yfivalues"
	foreach yfiv $yfivalues {
	    if {$yfiv ==$cycvalue} {
		set cyf "$yfiv"
		break
	    }
	}
	# wa "cyf:$cyf"
	# Get other properties for the specific yield function
	if {$cyf !=""} {
	    set cprop [::wkcf::GetMaterialPropertiesFromAttributes $AppId $mpxpath $yfid $cyf]
	    if {[llength $cprop]} {
		set matplist [concat $matplist $cprop]
	    }
	}

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
	set chl ""
	foreach mhiv $mhivalues {
	    # WarnWinText "mhiv:$mhiv"
	    if {$mhiv ==$chvalue} {
		set chl "$mhiv"
		break
	    }
	}
	# WarnWinText "chl:$chl"
	if {$chl !=""} {
	    set cprop [::wkcf::GetMaterialPropertiesFromAttributes $AppId $mpxpath $hlid $chl]
	    if {[llength $cprop]} {
		set matplist [concat $matplist $cprop]
	    }
	}

	# Get the SATURATION LAW
	set cslid "SaturationLaw"
	set msaturationlaw [::xmlutils::getKKWord $clxpath $ptype "msaturationlaw"]
	# Get the saturation law xpath values
	set msxpath [::xmlutils::getKKWord $clxpath $ptype "msxpath"]
	# WarnWinText "msaturationlaw:$msaturationlaw msxpath:$msxpath"
	# Get the current saturation law
	set csvalue [lindex [::KMat::getMaterialProperties "p" "$mpxpath//$msxpath//p.[list $msaturationlaw]"] 0 1]
	# WarnWinText "msaturationlaw:$msaturationlaw msxpath:$msxpath csvalue:$csvalue"
	# Get the saturation lawoptions
	set msivalues [split [::xmlutils::getKKWord "$clxpath//$cslid" "AvailableSaturationLaws" "msivalues"] ,]
	set csl ""
	foreach msiv $msivalues {
	    # wa "msiv:$msiv "
	    if {$msiv ==$csvalue} {
		set csl "$msiv"
		break
	    }
	}
	# wa "csl:$csl"
	if {$chl !=""} {
	    set cprop [::wkcf::GetMaterialPropertiesFromAttributes $AppId $mpxpath $cslid $csl]
	    if {[llength $cprop]} {
		set matplist [concat $matplist $cprop]
	    }
	}

	# Get the FLOW RULE
	set mflowrule [::xmlutils::getKKWord $clxpath $ptype "mflowrule"]
	# Get the flow rule xpath values
	set mfrxpath [::xmlutils::getKKWord $clxpath $ptype "mfrxpath"]
	# WarnWinText "mflowrule:$mflowrule mfrxpath:$mfrxpath"
	# Get the current flow rule
	set cfrvalue [lindex [::KMat::getMaterialProperties "p" "$mpxpath//$mfrxpath//p.[list $mflowrule]"] 0 1]
	# Get the internal flow rule
	set msivalues [split [::xmlutils::getKKWord "$clxpath//$mflowrule" "AvailableFlowRules" "msivalues"] ,]
	# WarnWinText "msivalues:$msivalues\n$mpxpath//$mfrxpath//p.[list $mflowrule] cfrvalue:$cfrvalue"
	set frl ""
	foreach mfriv $msivalues {
	    if {$mfriv ==$cfrvalue} {
		set frl "$mfriv"
		break
	    }
	}
	# wa "frl:$frl"
	if {$frl !=""} {
	    set cprop [::wkcf::GetMaterialPropertiesFromAttributes $AppId $mpxpath $mflowrule $frl]
	    if {[llength $cprop]} {
		set matplist [concat $matplist $cprop]
	    }
	}
    }

    # wa "matplist:$matplist"
    set dprops($AppId,Material,$MatId,Props) $matplist

}

proc ::wkcf::GetMaterialPropertiesFromAttributes {AppId mpxpath containerid attributeid} {
    # Get other material properties (constitutive laws)

    set matplist [list]

    # Xpath for constitutive laws
    set clxpath "CLawProperties"

    # Xpath for materials
    set matxpath "Materials"

    # Set the kratos model keyword base xpath
    set kmxpath "Applications//$AppId"

    # Get the property internal values
    set privalues [split [::xmlutils::getKKWord "$clxpath//$containerid" "$attributeid" "privalues"] ,]
    # Get the write yield function parameters
    set prxpath [split [::xmlutils::getKKWord "$clxpath//$containerid" "$attributeid" "prxpath"] ,]
    # WarnWinText "$clxpath//$containerid prxpath:$prxpath privalues:$privalues"
    foreach pid $privalues xpath $prxpath {
	# Get the kratos key word
	set kkword [::xmlutils::getKKWord $matxpath $pid "kkword"]
	# WarnWinText "pid:$pid kkword:$kkword"
	# Get the current value for this properties
	# WarnWinText "xpath:$mpxpath//$xpath//p.[list $pid]"
	set cvalue [lindex [::KMat::getMaterialProperties "p" "$mpxpath//$xpath//p.[list $pid]"] 0 1]
	# wa "cvalue:$cvalue"
	if {($kkword !="") && ($cvalue !="")} {
	    lappend matplist [list $kkword $cvalue]
	}
    }
    return $matplist
}

proc ::wkcf::SumInertia {InertiaIx InertiaIy} {

    return [expr $InertiaIx + $InertiaIy]
}

proc ::wkcf::GetCrossSectionProperties {AppId propid belemtype} {
    # Get the cross section properties
    # Arguments
    # AppId      => Application identifier
    # propid     => Current property identifier
    # belemtype  => Base element type
    variable dprops
    variable ndime

    # wa "MaterialProperties =>AppId:$AppId propid:$propid belemtype:$belemtype"
    # Xpath for constitutive laws
    set clxpath "CLawProperties"

    # Xpath for element constitutive laws
    set cexpath "ElementCLaws"

    # Set the kratos model keyword base xpath
    set kmxpath "Applications//$AppId"

    # Get others properties
    # wa "clxpath:$clxpath belemtype:$belemtype"
    # Get the property list
    set pid "propertylist"
    set id "CSProperty"
    set CSProperty [split [::xmlutils::getKKWord $cexpath "$id" "$pid"] ","]
    # wa "CSProperty:$CSProperty"
    # Get the property type list
    set pid "propertytypelist"
    set CSPropertyType [split [::xmlutils::getKKWord $cexpath "$id" "$pid"] ","]
    # wa "CSPropertyType:$CSPropertyType"

    # Get the section type list
    set id "SectionTypes"
    set pid "ivalues"
    set SectionTypes [split [::xmlutils::getKKWord $cexpath "$id" "$pid"] ","]
    # wa "SectionTypes:$SectionTypes"

    # Get the active property list and the cross section property list for this property identifier
    set CPropsList [list]
    set CrossSectionPropsList [list]
    foreach cspropid $CSProperty csproptypeid $CSPropertyType {
	# wa "cspropid:$cspropid csproptypeid:$csproptypeid"
	if {$cspropid eq "Thickness"} {
	    set endcspropid "${cspropid}${ndime}"
	} else {
	    set endcspropid $cspropid
	}
	# Get the base element type
	set pid "elementType"
	set ElementTypeList [split [::xmlutils::getKKWord $cexpath $endcspropid "$pid"] ","]
	# wa "ElementTypeList:$ElementTypeList"

	# Get the section type
	set pid "sectionType"
	set SectionTypeList [split [::xmlutils::getKKWord $cexpath $endcspropid "$pid"] ","]
	# wa "SectionTypeList:$SectionTypeList"
	set flag1 [expr {($cspropid eq "SectionType")||($cspropid eq "ProfileDB")}]
	set flag2 [expr {($SectionTypeList in $SectionTypes) && ($SectionTypeList !="")}]
	if {$flag1 || $flag2} {
	    if {$SectionTypeList eq "UserDefined"} {
		# Write to the material properties
		if {$belemtype in $ElementTypeList} {
		    lappend CPropsList [list $cspropid $csproptypeid]
		}
	    } elseif {($SectionTypeList eq "Rectangular")||($SectionTypeList eq "Circular")} {
		# Write to the cross section properties
		lappend CrossSectionPropsList [list $cspropid $csproptypeid $SectionTypeList]
	    } else {
		# For "SectionType" and "ProfileDB"
		# Write to the cross section properties
		lappend CrossSectionPropsList [list $cspropid $csproptypeid $SectionTypeList]
	    }
	} else {
	if {$belemtype in $ElementTypeList} {
	    lappend CPropsList [list $cspropid $csproptypeid]
	}
    }
    }
    # wa "CPropsList:$CPropsList"
    # wa "CrossSectionPropsList:$CrossSectionPropsList"

    # Init cross section properties list
    set dprops($AppId,CrossSection,$propid,CProps) [list]

    # Process all scalar and list data type only for property type that are different of "UserDefined"
    set UseCrossSectionProps 0

    foreach cpropid $CrossSectionPropsList {
	lassign $cpropid cspropid csproptypeid SectionType
	# wa "cspropid:$cspropid csproptypeid:$csproptypeid SectionType:$SectionType"
	# Get the current value
	set cvalue ""
	if {[info exists dprops($AppId,Property,$propid,$cspropid)]} {
	    set cvalue $dprops($AppId,Property,$propid,$cspropid)
	}
	# Check the section type
	if {($cspropid eq "SectionType") && ($cvalue !="UserDefined")} {
	    set UseCrossSectionProps 1

	    # Update the section type
	    set dprops($AppId,CrossSection,$propid,SectionType) $cvalue
	}
    }

    # Check for use cross section properties
    if {$UseCrossSectionProps} {
	foreach cpropid $CrossSectionPropsList {
	    lassign $cpropid cspropid csproptypeid SectionType
	    # wa "cspropid:$cspropid csproptypeid:$csproptypeid SectionType:$SectionType"
	    # Get the current value
	    set cvalue ""
	    if {[info exists dprops($AppId,Property,$propid,$cspropid)]} {
		set cvalue $dprops($AppId,Property,$propid,$cspropid)
	    }
	    # wa "cvalue:$cvalue"
	    # Check the section type
	    if {($SectionType !="UserDefined") && ($cvalue !="")} {
		if {$csproptypeid eq "Scalar"} {
		    # Update cross section properties (property type identifier, current value, section type)
		    lappend dprops($AppId,CrossSection,$propid,CProps) [list $cspropid $csproptypeid $cvalue $SectionType]
		} elseif {$csproptypeid eq "List"} {
		    # Update cross section properties (property type identifier, current value, section type)
		    lappend dprops($AppId,CrossSection,$propid,CProps) [list $cspropid $csproptypeid $cvalue $SectionType]
		}
		# wa "dprops($AppId,CrossSection,$propid,CProps):$dprops($AppId,CrossSection,$propid,CProps)"
	    }
	}
    }

    # Init CProps list
    set dprops($AppId,Material,$propid,CProps) [list]

    # Process all scalar data type
    foreach cpropid $CPropsList {
	lassign $cpropid cspropid csproptypeid
	# wa "cspropid:$cspropid csproptypeid:$csproptypeid"
	if {$csproptypeid eq "Scalar"} {
	    # Get the current value
	    set cvalue $dprops($AppId,Property,$propid,$cspropid)
	    set kword [::xmlutils::getKKWord $kmxpath $cspropid "kkword"]
	    if {$kword ne ""} {
	    # Update section properties (property type identifier, kratos keyword, current value)
	    lappend dprops($AppId,Material,$propid,CProps) [list $csproptypeid $kword $cvalue]
	    # wa "dprops($AppId,Material,$propid,CProps):$dprops($AppId,Material,$propid,CProps)"
	}
    }
    }

    # Process all matrix data type
    # a) First get all defined matrix identifier
    set AllMatrixIdList [list]
    foreach cpropid $CPropsList {
	lassign $cpropid cspropid csproptypeid
	if {$csproptypeid eq "Matrix"} {
	    # Get the matrix identifier for this property
	    set id "$cspropid"
	    set pid "matrixid"
	    set MatrixId [::xmlutils::getKKWord $cexpath "$id" "$pid"]
	    # wa "MatrixId:$MatrixId"
	    if {$MatrixId ni $AllMatrixIdList} {
		lappend AllMatrixIdList $MatrixId
	    }
	}
    }
    # wa "AllMatrixIdList:$AllMatrixIdList"

    # b) For each matrix identifier create the matrix properties
    set CMatrixProp [list]
    foreach MatrixId $AllMatrixIdList {

	# Get the keyword for this matrix identifier
	set id "$MatrixId"
	set pid "matrixkwordid"
	set MatrixKWordId [::xmlutils::getKKWord $cexpath "$id" "$pid"]
	# wa "MatrixKWordId:$MatrixKWordId"
	# Get the kratos keyword for this matrix
	set kword [::xmlutils::getKKWord $kmxpath $MatrixKWordId "kkword"]

	# Get the matrix dimension for this matrix identifier
	set id "$MatrixId"
	set pid "matrixndim"
	set MatrixNDim [split [::xmlutils::getKKWord $cexpath "$id" "$pid"] ","]
	lassign $MatrixNDim nrows ncols
	# wa "MatrixNDim:$MatrixNDim nrows:$nrows ncols:$ncols"

	# Get the matrix components for this matrix identifier
	set id "$MatrixId"
	set pid "matrixcomp"
	set MatrixComp [split [::xmlutils::getKKWord $cexpath "$id" "$pid"] ","]
	# wa "MatrixComp:$MatrixComp"

	# Get the current value for each component
	set cvaluelist [list]
	foreach ccompid $MatrixComp {
	    set cvalue "0.0"
	    if {[info exists dprops($AppId,Property,$propid,$ccompid)]} {
		set cvalue $dprops($AppId,Property,$propid,$ccompid)
	    }
	    lappend cvaluelist $cvalue
	}
	# wa "cvaluelist:$cvaluelist"
	# Get the matrix components position for this matrix identifier
	set id "$MatrixId"
	set pid "matrixcomppos"
	set MatrixCompPos [split [::xmlutils::getKKWord $cexpath "$id" "$pid"] ","]
	# wa "MatrixCompPos:$MatrixCompPos"

	# Get the matrix default values for this matrix identifier
	set id "$MatrixId"
	set pid "matrixdefvals"
	set MatrixDefVals [split [::xmlutils::getKKWord $cexpath "$id" "$pid"] ","]
	# wa "MatrixDefVals:$MatrixDefVals"

	# Set the values for the matrix properties
	set matbf ""
	set cmatrixdim "[join $MatrixNDim ","]"
	append matbf "\[$cmatrixdim\]" " "

	# Init the matrix with the default values
	set MatrixValues $MatrixDefVals
	# For each component assign the correct value
	set cpos 0
	foreach mcpos $MatrixCompPos {
	    # Get the value
	    set cvalue [lindex $cvaluelist $cpos]
	    incr cpos 1
	    # Modify the matrix value
	    lset MatrixValues [expr $mcpos-1] $cvalue
	}
	# wa "MatrixValues:$MatrixValues"

	# Resolve the special case that the properties are calculated
	# Get the matrix special calculation flag
	set id "$MatrixId"
	set pid "matrixspecialcal"
	set MatrixSpecialCal [split [::xmlutils::getKKWord $cexpath "$id" "$pid"] ","]
	# wa "MatrixSpecialCal:$MatrixSpecialCal"

	# Get the matrix special procedure identifiers
	set id "$MatrixId"
	set pid "matrixspecialpid"
	set MatrixSpecialProcId [split [::xmlutils::getKKWord $cexpath "$id" "$pid"] ","]
	# wa "MatrixSpecialProcId:$MatrixSpecialProcId"
	set cpos -1
	foreach spcal $MatrixSpecialCal spprocid $MatrixSpecialProcId {
	    # wa "spcal:$spcal spprocid:$spprocid"
	    incr cpos 1
	    if {$spcal} {
		# Get the procedure args identifier
		set existsproc [info procs $spprocid]
		if {$existsproc !=""} {
		    set argslist [info args $spprocid]
		    # wa "existsproc:$existsproc argslist:$argslist"
		    set pid "$spprocid"
		    # Create the procedure
		    foreach args $argslist {
		        # Get this property
		        set cvalue "0.0"
		        if {[info exists dprops($AppId,Property,$propid,$args)]} {
		            set cvalue $dprops($AppId,Property,$propid,$args)
		        }
		        # wa "cvalue:$cvalue"
		        append pid " " $cvalue
		    }
		    # wa "pid:$pid"
		    set currentvalue [eval $pid]
		    # wa "currentvalue:$currentvalue"
		    if {$currentvalue !=""} {
		        # Modify the matrix value
		        lset MatrixValues $cpos $currentvalue
		    }
		}
	    }
	}

	# Print the value into the buffer
	# Open the matrix brackets
	append matbf "(("
	set cpos 1
	set vallen [llength $MatrixValues]
	set valcount 0
	foreach cval $MatrixValues {
	    incr valcount 1
	    if {$cpos<$ncols} {
		append matbf "$cval,"
		incr cpos 1
	    } elseif {$cpos==$ncols} {
		if {$valcount == $vallen} {
		    append matbf "$cval)"
		    break
		} else {
		    append matbf "$cval),("
		}
		set cpos 1
	    }
	}
	# End the matrix brackets
	append matbf ")"
	# wa "matbf:$matbf"
	# Set the matrix properties
	lappend CMatrixProp [list "Matrix" "$kword" "$matbf"]
    }
    # wa "CMatrixProp:$CMatrixProp"

    # c) Update the cross section properties array
    if {[llength $CMatrixProp]} {
	# Update section properties (property type identifier, kratos keyword, current value)
	foreach MProp $CMatrixProp {
	    # wa "MProp:$MProp"
	    lappend dprops($AppId,Material,$propid,CProps) $MProp
	}
    }
    # wa "dprops($AppId,Material,$propid,CProps):$dprops($AppId,Material,$propid,CProps)"

}

proc ::wkcf::GetElementProperties {} {
    variable dprops
    variable ActiveAppList

    # For each active application
    foreach AppId $ActiveAppList {
	# Get the application root identifier
	set rootdataid $AppId
	# Get the properties element links
	# All kratos element identifier
	set dprops($AppId,AllKElemId) [list]
	# All kratos property identifier
	set dprops($AppId,AllKPropertyId) [list]
	# All Kratos element group identifier
	set dprops($AppId,AllKEGroupId) [list]
	# Get all defined element types
	set cbasexpath "$rootdataid//c.Elements"
	if {$AppId eq "DEM"} {
	    set cbasexpath "$rootdataid//c.DEM-Elements"
	}
	set glist [::xmlutils::setXmlContainerIds $cbasexpath]
	set kwxpath "Applications/$rootdataid"
	# WarnWinText "AppId:$AppId"
	# WarnWinText "glist:$glist"

	foreach celemid $glist {
	    # Get the group identifier defined for this element
	    set cxpath "$cbasexpath//c.[list ${celemid}]"
	    set cgrouplist [::xmlutils::setXmlContainerIds $cxpath]
	    if {[llength $cgrouplist]>0} {
		# WarnWinText "cgrouplist:$cgrouplist"
		foreach cgroupid $cgrouplist {
		    # Check if this group is in active state
		    # Get active property
		    set cxpath "$cbasexpath//c.[list ${celemid}]//c.[list ${cgroupid}]"
		    set cproperty "active"
		    set ActiveGroup [::xmlutils::setXml $cxpath $cproperty]
		    # WarnWinText "ActiveGroup:$ActiveGroup"
		    if {$ActiveGroup =="0"} {
		        continue
		    }
		    # Update all kratos element identifier
		    if {$celemid ni $dprops($AppId,AllKElemId)} {
		        lappend dprops($AppId,AllKElemId) $celemid
		    }
		    # Update all group identifier
		    if {$cgroupid ni $dprops($AppId,AllKEGroupId)} {
		        lappend dprops($AppId,AllKEGroupId) $cgroupid
		    }
		    # Get the GiD entity type
		    set cxpath "$cbasexpath//c.[list ${celemid}]"
		    set cproperty "GiDEntity"
		    set GiDEntity [::xmlutils::setXml $cxpath $cproperty]
		    # WarnWinText "GiDEntity:$GiDEntity"
		    # Get the GiD element type for Kratos elements type
		    set cxpath "$cbasexpath//c.[list ${celemid}]//c.[list ${cgroupid}]//c.Properties//i.ElementType"
		    set cproperty "dv"
		    set GiDElemType [::xmlutils::setXml $cxpath $cproperty]
		    # WarnWinText "GiDElemType:$GiDElemType"
		    # Get the property identifier
		    set cxpath "$cbasexpath//c.[list ${celemid}]//c.[list ${cgroupid}]//c.Properties//i.Property"
		    set cproperty "dv"
		    set PropertyId [::xmlutils::setXml $cxpath $cproperty]
		    # WarnWinText "PropertyId:$PropertyId"
		    # Get the Key word for the element type
		    set kelemtype [::xmlutils::getKKWord $kwxpath $celemid "kkword"]
		    # set cxpath "$rootdataid//c.Elements//c.[list ${celemid}]"
		    # set cproperty "kkword"
		    # set kelemtype [::xmlutils::setXml $cxpath $cproperty]
		    # WarnWinText "kelemtype:$kelemtype"
		    # Get ndim
		    set cxpath "$cbasexpath//c.[list ${celemid}]"
		    set cproperty "nDim"
		    set nDim ""; #[::xmlutils::setXml $cxpath $cproperty]
		    #W "nDim:$nDim"

		    set GProps [list $GiDEntity $GiDElemType $PropertyId $kelemtype $nDim]

		    # Kratos element to group link
		    # Group list
		    if {![info exists dprops($AppId,KElem,$celemid,AllGroupId)]} {
		        set dprops($AppId,KElem,$celemid,AllGroupId) [list]
		    }
		    if {$cgroupid ni $dprops($AppId,KElem,$celemid,AllGroupId)} {
		        lappend dprops($AppId,KElem,$celemid,AllGroupId) $cgroupid
		    }
		    # Group properties
		    if {![info exists dprops($AppId,KElem,$celemid,$cgroupid,GProps)]} {
		        set dprops($AppId,KElem,$celemid,$cgroupid,GProps) [list]
		    }
		    set dprops($AppId,KElem,$celemid,$cgroupid,GProps) $GProps

		    # Update AllKPropertyId list
		    if {$PropertyId ni $dprops($AppId,AllKPropertyId)} {
		        lappend dprops($AppId,AllKPropertyId) $PropertyId
		    }
		    # Property to group identifier link
		    if {![info exists dprops($AppId,Property,$PropertyId,GroupId)]} {
		        set dprops($AppId,Property,$PropertyId,GroupId) [list]
		    }
		    # Update the link between group and property => if the property exist add to the list
		    if {[info exists dprops($AppId,Property,$PropertyId,GroupId)]} {
		        lappend dprops($AppId,Property,$PropertyId,GroupId) $cgroupid
		    } else {
		        set dprops($AppId,Property,$PropertyId,GroupId) $cgroupid
		    }

		    # Update the link between the element type and property => if the property exist add to the list
		    if {[info exists dprops($AppId,Property,$PropertyId,ElementId)]} {
		        lappend dprops($AppId,Property,$PropertyId,ElementId) $celemid
		    } else {
		        set dprops($AppId,Property,$PropertyId,ElementId) $celemid
		    }
		}
	    }
	}
    }
}

proc ::wkcf::GetElementIdFromPropertyId {AppId PropertyId} {
    variable dprops

    set ElementId ""
    if {[info exists dprops($AppId,Property,$PropertyId,ElementId)]} {
       set ElementId $dprops($AppId,Property,$PropertyId,ElementId)
    }
    return $ElementId
}

proc ::wkcf::UnsetLocalVariables {} {
    variable dprops; variable AppId
    variable ctbclink

    if {[info exists AppId]} {
	unset AppId
    }
    if {[info exists dprops]} {
	unset dprops
    }

    if {[info exists ctbclink]} {
	unset ctbclink
    }
    variable ActiveAppList
    if {[info exists ActiveAppList]} {
	unset ActiveAppList
    }
}

proc ::wkcf::GetnDimnNode {GiDElemType nDim} {
    # ABSTRACT: Get the kratos element identifier as a function of GiD element type
    variable useqelem

    # wa "GiDElemType:$GiDElemType nDim:$nDim"
    set etbf ""
    set nid "3N"
    if {$nDim =="2D"} {
	if {$GiDElemType =="Triangle"} {
	    if {$useqelem=="1"} {
		set nid "6N"
	    } else {
		set nid "3N"
	    }
	} elseif {$GiDElemType =="Quadrilateral"} {
	    if {$useqelem=="1"} {
		set nid "8N"
	    } else {
		set nid "4N"
	    }
	} elseif {$GiDElemType =="Linear"} {
	    if {$useqelem=="1"} {
		set nid "3N"
	    } else {
		set nid "2N"
	    }
	}
    } elseif {$nDim =="3D"} {
	if {$GiDElemType =="Tetrahedra"} {
	    if {$useqelem=="1"} {
		set nid "10N"
	    } else {
		set nid "4N"
	    }
	} elseif {$GiDElemType =="Hexahedra"} {
	    if {$useqelem=="1"} {
		set nid "20N"
	    } elseif {$useqelem=="2"} {
		set nid "27N"
	    } else {
		set nid "8N"
	    }
	} elseif {$GiDElemType =="Triangle"} {
	    if {$useqelem=="1"} {
		set nid "6N"
	    } else {
		set nid "3N"
	    }
	} elseif {$GiDElemType =="Quadrilateral"} {
	    if {$useqelem=="1"} {
		set nid "8N"
	    } else {
		set nid "4N"
	    }
	} elseif {$GiDElemType =="Linear"} {
	    if {$useqelem=="1"} {
		set nid "3N"
	    } else {
		set nid "2N"
	    }
	}
    }
    set etbf "${nDim}${nid}"
    return $etbf
}

proc ::wkcf::GetSurfaceTypeList {surfacelist} {

    # set vollist [GiD_Geometry list volume 1:end]
    # if {[llength $vollist]} {
    #         # Get the element type
    #         foreach volid $vollist {
    #             set cvprop [GiD_Info list_entities volumes $volid]
    #             # wa "volid:$volid cvprop:$cvprop"
    #             regexp -nocase {Elemtype=([0-9]*)} $cvprop none voltype
    #             # wa "voltype:$voltype"
    #         }
    # }
    set tetrasurf [list]
    set hexasurf [list]
    foreach surfid $surfacelist {
	# Check for higher entity
	set cprop [GiD_Info list_entities -More surfaces $surfid]
	# wa "surfid:$surfid cprop:$cprop"
	set isve 0
	regexp -nocase {higherentity: ([0-9]+)} $cprop none ivhe
	# wa "ivhe:$ivhe"
	if {$ivhe} {
	    set he [regexp -nocase {Higher entities volumes: (.)*} $cprop vol]
	    # wa "he:$he vol:$vol vlist:[lrange $vol 3 end-2]"
	    if {$he && $vol !=""} {
	    set voltype ""
	    set vlist [lindex [lrange $vol 3 end-2] 0]
	    # wa "vlist:$vlist"
	    set cvprop [GiD_Info list_entities volumes $vlist]
	    # wa "cvprop:$cvprop"
	    regexp -nocase {Elemtype=([0-9]*)} $cvprop none voltype
	    # wa "voltype:$voltype"
	    if {($voltype == 4) || ($voltype == 0) || ($voltype=="")} {
		lappend tetrasurf $surfid
	    } elseif {$voltype == 5} {
		lappend hexasurf $surfid
	    } else {
		lappend tetrasurf $surfid
	    }
	    }
	}
    }
    return [list $tetrasurf $hexasurf]
}

proc ::wkcf::FindBoundaries {entity} {
    # ABSTRACT: Return a list containing all boundaries entities
    # Arguments
    # entity => Entity to be processed
    #  * entity=line for models made of surfaces
    #  * entity=surface for models made of volumes
    # Note: This procedure in the same used in the fluid_only problem type

    set boundarylist [list]

    # Generate some names
    set Entity [string toupper $entity 0 0]
    set entities [format "%ss" $entity]

    # Get the number of the last entity
    set instruction [format "MaxNum%ss" $Entity]
    set Max [GiD_Info Geometry $instruction]
    # Generate a list containing all entities and record their ids and number of HigherEntities
    set EntityList [GiD_Info list_entities $entities 1:$Max]
    set candidates [regexp -all -inline {Num: ([0-9]*) HigherEntity: ([0-9]*)} $EntityList]
    # Find ids of entities with exactly 1 HigherEntity (this means they are in the boundary)
    for {set i 1} {$i < [llength $candidates]} {incr i 3} {
	set j [expr {$i + 1}]
	if {[lindex $candidates $j] == 1} {lappend boundarylist [lindex $candidates $i]}
    }
    return $boundarylist
}

proc ::wkcf::FindBoundariesOfNonSphericElements {entity} {
    # ABSTRACT: Return a list containing all boundaries entities
    # Arguments
    # entity => Entity to be processed
    #  * entity=line for models made of surfaces
    #  * entity=surface for models made of volumes
    # Note: This procedure in the same used in the fluid_only problem type

    set groups_to_spherize_list [::xmlutils::setXmlContainerIds {DEM//c.DEM-Elements//c.DEM-Element}]
    foreach volume_id [GiD_Geometry list volume 1:end] { ; #list of volume identifiers in the whole range
	set volume_info [GiD_Info list_entities volume $volume_id] ; #info about those volumes
	set is_spheric [regexp {Elemtype=9} $volume_info] ; #finding out if the element type is spheric

	foreach group_that_includes_this_volume [GiD_EntitiesGroups entity_groups volumes $volume_id] { ; #list of groups to which to $volume_id belongs
	#next we search $group_that_includes_this_volume among $groups_to_spherize_list:
	    if {[lsearch $groups_to_spherize_list $group_that_includes_this_volume] >= 0} {
		set is_spheric 1
	    }
	}

	if {$is_spheric==0} {
	    foreach item [lrange [GiD_Geometry get volume $volume_id] 2 end] {
		set surface_id [lindex $item 0]
		incr surfaces_higher_entities_list($surface_id)
	    }
	}
    }

    set boundarylist [list]
    foreach surface_id [lsort -integer [array names surfaces_higher_entities_list]] {
	if {$surfaces_higher_entities_list($surface_id) == 1} {
	    lappend boundarylist $surface_id
	}
    }
    return $boundarylist
}

proc ::wkcf::FindAllSurfacesOfNonSphericElements {entity} {
    # ABSTRACT: Return a list containing all boundaries entities
    # Arguments
    # entity => surface

	set surf_high_entities [list]
	set surf_no_high_entities [list]
	set boundarylist [list]

	# Boundary surfaces of all the volumes in the domain
    foreach volume_id [GiD_Geometry list volume 1:end] { ; #list of volume identifiers in the whole range
		set volume_info [GiD_Info list_entities volume $volume_id] ; #info about those volumes
		set is_spheric [regexp {Elemtype=9} $volume_info] ; #finding out if the element type is spheric

		# Sphere volumes are excluded
		if {$is_spheric==0} {
		    foreach item [lrange [GiD_Geometry get volume $volume_id] 2 end] {
		        lappend surf_high_entities [lindex $item 0]
		   }
		}
    }

	# Surfaces with no higher entities (not belonging to a volume)
	set layers [GiD_Info layers]
	foreach layer $layers {
		lappend surf_no_high_entities [GiD_Info layers -entities surfaces -higherentity 0 $layer]
	}

	set boundarylist [concat {*}$surf_high_entities {*}$surf_no_high_entities]
	#W "boundarylist surfaces: $boundarylist"

    return $boundarylist
}

proc ::wkcf::FindBoundariesOfSphericElements {entity} {

    set groups_to_spherize_list [::xmlutils::setXmlContainerIds {DEM//c.DEM-Elements//c.DEM-Element}]
    foreach volume_id [GiD_Geometry list volume 1:end] { ; #list of volume identifiers in the whole range
	set volume_info [GiD_Info list_entities volume $volume_id] ; #info about those volumes
	set is_spheric [regexp {Elemtype=9} $volume_info] ; #finding out if the element type is spheric

	foreach group_that_includes_this_volume [GiD_EntitiesGroups entity_groups volumes $volume_id] { ; #list of groups to which to $volume_id belongs
	#next we search $group_that_includes_this_volume among $groups_to_spherize_list:
	    if {[lsearch $groups_to_spherize_list $group_that_includes_this_volume] >= 0} {
		set is_spheric 1
	    }
	}

	if {$is_spheric==1} {
	    foreach item [lrange [GiD_Geometry get volume $volume_id] 2 end] {
		set surface_id [lindex $item 0]
		incr surfaces_higher_entities_list($surface_id)
	    }
	}
    }

    set boundarylist [list]
    foreach surface_id [lsort -integer [array names surfaces_higher_entities_list]] {
	if {$surfaces_higher_entities_list($surface_id) == 1} {
	    lappend boundarylist $surface_id
	}
    }
    return $boundarylist
}

proc ::wkcf::FindBoundariesOfCircularElements {entity} {

    set groups_to_circularize_list [::xmlutils::setXmlContainerIds {DEM//c.DEM-Elements//c.DEM-Element}]
    foreach surface_id [GiD_Geometry list surface 1:end] { ; #list of surface identifiers in the whole range
	set surface_info [GiD_Info list_entities surface $surface_id] ; #info about those surfaces
	set is_circular [regexp {Elemtype=10} $surface_info] ; #finding out if the element type is circular

	foreach group_that_includes_this_surface [GiD_EntitiesGroups entity_groups surfaces $surface_id] {
	#next we search $group_that_includes_this_surface among $groups_to_circularize_list:
	    if {[lsearch $groups_to_circularize_list $group_that_includes_this_surface] >= 0} {
		set is_circular 1
	    }
	}

	set number_of_lines_in_the_surface [lindex [GiD_Geometry get surface $surface_id] 2]

	if {$is_circular==1} {
	    foreach item [lrange [GiD_Geometry get surface $surface_id] 9 [expr {8 + $number_of_lines_in_the_surface}]] {
		set line_id [lindex $item 0]
		incr lines_higher_entities_list($line_id)
	    }
	}
    }

    set boundarylist [list]
    foreach line_id [lsort -integer [array names lines_higher_entities_list]] {
	if {$lines_higher_entities_list($line_id) == 1} {
	    lappend boundarylist $line_id
	}
    }
    return $boundarylist
}

proc ::wkcf::AssignGeometricalEntitiesToSkinSphere2D {} {

    set list_of_points [GiD_Geometry list point 1:end]
    set list_of_lines [GiD_Geometry list line 1:end]
    if {![GiD_Groups exists SKIN_SPHERE_DO_NOT_DELETE]} {
	GiD_Groups create SKIN_SPHERE_DO_NOT_DELETE
    }

    set points_to_add_to_skin_circles [list]
    set lines_to_add_to_skin_circles [list]
    set boundary_circle_line_list [::wkcf::FindBoundariesOfCircularElements line]

    set total_skin_line_circle_list [concat $lines_to_add_to_skin_circles $boundary_circle_line_list]
    set total_skin_circle_list [list $points_to_add_to_skin_circles $total_skin_line_circle_list {} {}]
    GiD_EntitiesGroups assign SKIN_SPHERE_DO_NOT_DELETE all_geometry $total_skin_circle_list
}

proc ::wkcf::AssignGeometricalEntitiesToSkinSphere3D {} {

    set list_of_points [GiD_Geometry list point 1:end]
    set list_of_lines [GiD_Geometry list line 1:end]
    set list_of_surfaces [GiD_Geometry list surface 1:end]
    if {![GiD_Groups exists SKIN_SPHERE_DO_NOT_DELETE]} {
	GiD_Groups create SKIN_SPHERE_DO_NOT_DELETE
    }

    set points_to_add_to_skin_spheres [list]
    set lines_to_add_to_skin_spheres [list]
    set surfaces_to_add_to_skin_spheres [list]
    set bound_sphere_surface_list [::wkcf::FindBoundariesOfSphericElements surface]

    foreach point_id $list_of_points line_id $list_of_lines surface_id $list_of_surfaces {
	set point_info [GiD_Info list_entities point $point_id]
	set line_info [GiD_Info list_entities line $line_id]
	set surface_info [GiD_Info list_entities surface $surface_id]
	set point_has_no_higher_entities [regexp {HigherEntity: 0} $point_info]
	set line_has_no_higher_entities [regexp {HigherEntity: 0} $line_info]
	set surface_has_no_higher_entities [regexp {HigherEntity: 0} $surface_info]
	if {$point_has_no_higher_entities == 1} {
	    lappend points_to_add_to_skin_spheres $point_id
	}
	if {$line_has_no_higher_entities == 1} {
	    lappend lines_to_add_to_skin_spheres $line_id
	}
	if {$surface_has_no_higher_entities == 1} {
	    lappend surfaces_to_add_to_skin_spheres $surface_id
	}
    }

    set total_skin_surface_sphere_list [concat $surfaces_to_add_to_skin_spheres $bound_sphere_surface_list]
    set total_skin_sphere_list [list $points_to_add_to_skin_spheres $lines_to_add_to_skin_spheres $total_skin_surface_sphere_list {}]
    GiD_EntitiesGroups assign SKIN_SPHERE_DO_NOT_DELETE all_geometry $total_skin_sphere_list
}

proc ::wkcf::AlignLineNormals {direction} {
    # ABSTRACT: Makes all of boundary lines' normals point inwards or outwards
    # Arguments
    # direction => Direction option ["Inwards"|"Outwards"]
    # Note: This procedure in the same used in the fluid_only problem type

    switch $direction {
	Inwards {
	    set wrong_way "DIFF1ST"
	}
	Outwards {
	    set wrong_way "SAME1ST"
	}
	default {puts "Unknown direction, line normals not aligned"}
    }

    set surfacelist [GiD_Geometry list surface 1:]

    # For each surface, we look for boundary lines oriented in the wrong direction
    set linelist [list]
    foreach surface $surfacelist {
	set surfaceinfo [GiD_Info list_entities surfaces $surface]
	set numpos [lsearch $surfaceinfo "NumLines:"]
	set numlines [lindex $surfaceinfo [expr {$numpos +1}]]
	for {set i 0} {$i < $numlines} {incr i} {
	    set orient [lindex $surfaceinfo [expr {$numpos+5+4*$i}]]
	    if {[string compare $orient $wrong_way]==0} {
		# If the normal is pointing in the wrong direction,
		# Check if it's a contour line
		set linenum [lindex $surfaceinfo [expr {$numpos+3+4*$i}]]
		set lineinfo [GiD_Info list_entities lines $linenum]
		#set highpos [lsearch $surfinfo "HigherEntity:"]
		set higherentities [lindex $lineinfo 4]
		if {$higherentities==1} {
		    lappend linelist $linenum
		}
	    }
	}
    }

    if {[llength $linelist]} {
	# If its in the contour, switch its normal
	eval GiD_Process Mescape Utilities SwapNormals Lines Select $linelist
    }

}

proc ::wkcf::AlignSurfNormals {direction} {
    # ABSTRACT: Makes all of boundary surfaces' normals point inwards or outwards
    # Arguments
    # direction => Direction option ["Inwards"|"Outwards"]
    # Note: This procedure in the same used in the fluid_only problem type

    switch $direction {
	Inwards {
	    set wrong_way "DIFF1ST"
	}
	Outwards {
	    set wrong_way "SAME1ST"
	}
	default {puts "Unknown Direction, surface normals not aligned"}
    }

    set volumelist [GiD_Geometry list volume 1:]

    set surfacelist [list]
    # For each volume, we look for face surfaces with oriented in the wrong direction
    foreach volume $volumelist {
	set volumeinfo [GiD_Info list_entities volumes $volume]
	set numpos [lsearch $volumeinfo "NumSurfaces:"]
	set numsurf [lindex $volumeinfo [expr {$numpos +1 }]]
	for {set i 0} {$i < $numsurf} {incr i} {
	    set orient [lindex $volumeinfo [expr {$numpos+5+4*$i}]]
	    if {[string compare $orient $wrong_way]==0} {
		# If the normal is pointing in the wrong direction,
		# Check if it's a contour surface
		set surfnum [lindex $volumeinfo [expr {$numpos+3+4*$i}]]
		set surfinfo [GiD_Info list_entities surfaces $surfnum]
		set higherentities [lindex $surfinfo 4]
		if {$higherentities==1} {
		 lappend surfacelist $surfnum
		}
	    }
	}
    }

    if {[llength $surfacelist]} {
	# If its in the contour, switch its normal
	eval GiD_Process Mescape Utilities SwapNormals Surfaces Select $surfacelist
    }
}

proc ::wkcf::WriteBatFile {AppId} {
    # ABSTRACT: Write the Kratos bat files

    # Get the parallel solution type
    set rootid "$AppId"
    # Kratos key word xpath
    set kxpath "Applications/$rootid"
    set cproperty "dv"
    set cxpath "$rootid//c.SolutionStrategy//c.ParallelType//i.ParallelSolutionType"
    set ParallelSolutionType [::xmlutils::setXml $cxpath $cproperty]

    if {$ParallelSolutionType eq "MPI"} {
	#  Get the number of processors
	set cxpath "$rootid//c.SolutionStrategy//c.ParallelType//i.MPINumberOfProcessors"
	set MPINumberOfProcessors [::xmlutils::setXml $cxpath $cproperty]
	# wa "MPINumberOfProcessors:$MPINumberOfProcessors"

	if {($::tcl_platform(os) eq "Linux")} {
	    # Linux
	    set batfilename "kratos-mpi.unix.bat"
	    set ProblemTypePath [::KUtils::GetPaths "PTDir"]
	    set batfullname [file native [file join $ProblemTypePath $batfilename]]

	    # First delete the file
	    set res ""
	    catch { set res [file delete -force $batfullname] }

	    # Create the new file
	    set f [open $batfullname w]
	    # WarnWinText "batfullname:$batfullname res:$res f:$f"

	    puts $f "\#\!\/bin\/bash -i"
	    puts $f "#  echo hola"
	    puts $f "#  echo par1 : $1"
	    puts $f "#  echo \"1. param: $1\" > /tmp/kk.txt"
	    puts $f "#  echo \"2. param: $2\" >> /tmp/kk.txt"
	    puts $f "#  echo \"3. param: $3\" >> /tmp/kk.txt"
	    puts $f "#  echo \"4. param: $4\" >> /tmp/kk.txt"
	    puts $f "#  echo \"5. param: $5\" >> /tmp/kk.txt"
	    puts $f "#    OutputFile: \"$2/$1.info\""
	    puts $f "#    ErrorFile: \"$2/$1.err\""
	    puts $f "# Delete previous result file"
	    puts $f "rm -f -r \"$2/$1*.post.bin\""
	    puts $f "rm -f -r \"$2/$1*.post.res\""
	    puts $f "rm -f -r \"$2/$1*.post.msh\""
	    puts $f "rm -f -r \"$2/$1*.info\""
	    puts $f "rm -f \"$2/$1.err\""
	    puts $f "rm -f \"$2/$1.flavia.dat\""
	    puts $f ""
	    puts $f "export OMP_NUM_THREADS=1"
	    puts $f "mpirun -np $MPINumberOfProcessors /usr/bin/python KratosMPI.py >\"$2/$1.info\" 2>\"$2/$1.err\""

	} else {

	    # Windows

	    set batfilename "kratos-mpi.win.bat"
	    set ProblemTypePath [::KUtils::GetPaths "PTDir"]
	    set batfullname [file native [file join $ProblemTypePath $batfilename]]

	    # First delete the file
	    set res ""
	    catch { set res [file delete -force $batfullname] }

	    # Create the new file
	    set f [open $batfullname w]

	    puts $f "REM @ECHO OFF"
	    puts $f "REM Identification for arguments"
	    puts $f "REM basename                          = %1"
	    puts $f "REM Project directory                 = %2"
	    puts $f "REM Problem directory                 = %3"
	    puts $f " "
	    puts $f "REM OutputFile: %2\\%1.info"
	    puts $f "REM ErrorFile: %2\\%1.err"
	    puts $f " "
	    puts $f "DEL %2\\%1.info"
	    puts $f "DEL %2\\%1*.post.bin"
	    puts $f "DEL %2\\%1*.post.res"
	    puts $f "DEL %2\\%1*.post.msh"
	    puts $f "DEL %2\\%1.flavia.dat"
	    puts $f "DEL %2\\%1.err"
	    puts $f "REM Run the python script"

	    if {$MPINumberOfProcessors>0} {
		puts $f "mpirun -np $MPINumberOfProcessors python kratosMPI.py > %2\\%1.info 2> %2\\%1.err"
	    }
	}

	close $f
    }

}

# From Fluid_only problem type: Unassigns automatically assigned GiD Conditions (Model Parts, Conditions, Elements) from given entity types
proc ::wkcf::CleanAutomatic {Condition args} {
    foreach entity $args {
	set autolist {}
	set infolist [GiD_Info conditions ${entity}_${Condition} geometry]
	foreach item $infolist {
	    set id [regexp -inline {^E ([0-9]*) - ([0-9]*)} $item]
	    # If it's "automatic" value is >0 (its always 0 for user-assigned data), store its id
	    if {[lindex $id 2] > 0} {lappend autolist [lindex $id 1]}
	}
	GiD_UnAssignData Condition ${entity}_${Condition} ${entity}s $autolist
    }
}

proc ::wkcf::CleanAutomaticConditionGroup {what args {fieldname ""} {fieldvalue ""}} {
    # Need New GiD_group adaptation
    set Condition "groups"

    switch -- -exact $what {
	"SelectEntity" {
	    foreach entity $args {
		set autolist {}
		set infolist [GiD_Info conditions ${entity}_${Condition} geometry]
		foreach item $infolist {
		    set id [regexp -inline {^E ([0-9]*) - ([0-9]*)} $item]
		    # If it's "automatic" value is >0 (its always 0 for user-assigned data), store its id
		    if {[lindex $id 2] > 0} {lappend autolist [lindex $id 1]}
		}
		# wa "autolist:$autolist"
		GiD_UnAssignData Condition ${entity}_${Condition} ${entity}s $autolist
	    }
	}
	"UseAllWhereField" {
	    foreach entity $args {
	    GiD_UnAssignData Condition ${entity}_${Condition} ${entity}s all $fieldname $fieldvalue
	    }
	}
    }
}

proc ::wkcf::CleanAutomaticConditionGroupGiD {args {fieldvalue ""}} {
    if {![GiD_Groups exists $fieldvalue]} {
	GiD_Groups create $fieldvalue
    }
    GiD_Groups edit state $fieldvalue hidden
    # msg [GiD_Groups get state $fieldvalue]
    # msg "$fieldvalue [GiD_EntitiesGroups get $fieldvalue elements]"
    foreach entity $args {
	GiD_EntitiesGroups unassign $fieldvalue $entity
    }
    GidUtils::UpdateWindow GROUPS
}

proc ::wkcf::AssignConditionToGroup {entity elist groupid} {
    # Need New GiD_group adaptation
    set Condition "groups"

    GiD_AssignData Condition ${entity}_${Condition} ${entity}s $groupid $elist

}
proc ::wkcf::AssignConditionToGroupGID {entity elist groupid} {
    # Need New GiD_group adaptation
    if {![GiD_Groups exists $groupid]} {
	GiD_Groups create $groupid
    }
    GiD_Groups edit state $groupid hidden
    # msg [GiD_Groups get state $groupid]
    # msg "$groupid [GiD_EntitiesGroups get $groupid elements]"
    GiD_EntitiesGroups assign $groupid $entity $elist
    GidUtils::UpdateWindow GROUPS
}

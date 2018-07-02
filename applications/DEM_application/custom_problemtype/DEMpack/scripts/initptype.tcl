
###############################################################################
#
#    NAME: kinitptype.tcl
#
#    PURPOSE: Init Kratos problem type
#
#    QUANTECH ATZ-DEVELOPMENT DEPARTMENT
#
#    AUTHOR : G. Socorro
#
#    CREATED AT: 01/02/10
#
#    HISTORY:
#
#     1.4-09/10/13-G. Socorro, add two new script file wkcfconvectiondiffusion.tcl wkcfdem.tcl
#     1.3-19/09/13-G. Socorro, modify the proc BeforeMeshGeneration to assign automatically triangle or quadrilateral element to the skin surfaces
#     1.2-18/06/13-G. Socorro, delete the proc kipt::NewGiDGroups
#     1.1-22/10/12-J. Garate, Support for new GiD Groups
#     1.0-08/10/12-J. Garate, Enable/disable kipt::NewGiDGroups
#     0.9-01/10/12-J. Garate, Enable/disable Curves Module
#     0.8-20/09/12-J. Garate, add Curves, Tables and Plotgraph source files
#     0.7-04/05/12-G. Socorro, add a new variable to control the group deletion (when exists from the problem type)
#     0.6-03/05/12-G. Socorro, Delete all group identifier using ::KUtils::DeleteAllGroupIdentifier and
#                              close all group window Cond_Groups window close
#     0.5-10/04/12-G. Socorro, load new script in the (wkcffluid.tcl, etc.)
#     0.4-02/04/12-J. Garate, icon path change to adapt to GiD Themes
#     0.3-29/03/12-G. Socorro, load the new kmprops scripts
#     0.2-22/06/11-G. Socorro, delete KPriv(release) and create KPriv(RDConfig) in the Kratos.tcl
#     0.1-01/02/10-G. Socorro, create the base source code
#
###############################################################################

#short alias
proc msg {message} {
    WarnWinText $message
}

proc msgS {message} {
    WarnWin $message
}

proc wa {message} {
    WarnWinText $message
}

# kipt => kratos Init Problem Type

#######################################################
### namespace Kratos ###
#######################################################

namespace eval kipt {} {
    variable ProgramName kratos
    variable VersionNumber; #interface version, get it from xml to avoid duplication
    variable Web http://www.cimne.com/kratos
}

proc kipt::Splash {} {
    global KPriv
    variable ProgramName
    variable VersionNumber

    set prev_splash_state [GiD_Set SplashWindow]
    GiD_Set SplashWindow 1 ;#set temporary to 1 to force show splash without take care of the GiD splash preference
    set off_x 150
    set fnt "Sans 10"
    if { $::tcl_platform(platform) == "windows" } {
	set fnt "verdana 10"
	set off_x 130
    }

    if {$KPriv(what_dempack_package) eq "C-DEMPack"} {
	set imagename splash_C.png
    } elseif { $KPriv(what_dempack_package) eq "G-DEMPack"} {
	set imagename splash_D.png
    } elseif { $KPriv(what_dempack_package) eq "S-DEMPack"} {
	set imagename splash_S.png
    } else {
	set imagename splash_F.png
    }

    ::GidUtils::Splash [file join $::KPriv(dir) images Classic $imagename] .splash 1 \
	[list "$::kipt::ProgramName Version $::kipt::VersionNumber \n$::kipt::Web" $off_x 230]
    # next time to be removed...
    # new_gid = -1 for gid v. 11.1.5d or below
    # dev_kit's splash has been modified and some changes are already there...
    set new_splash [ ::GidUtils::VersionCmp 11.1.6d]
    if {[ winfo exists .splash.lv] && ( $new_splash < 0)} {
	.splash.lv configure -font $fnt -background white -foreground black \
	    -relief solid -borderwidth 1 -padx 12 -pady 3
	update
    }
    GiD_Set SplashWindow $prev_splash_state
}

proc kipt::About {} {
    variable KratosPriv
    variable ProgramName
    variable VersionNumber
    set prev_splash_state [GiD_Set SplashWindow]
    GiD_Set SplashWindow 1
    set off_x 150
    set fnt "Sans 10"
    if { $::tcl_platform(platform) == "windows" } {
	set fnt "verdana 10"
	set off_x 130
    }
    ::GidUtils::Splash [file join $::KPriv(dir) images Classic splash.png] .splash 0 \
	[list "$::kipt::ProgramName Version $::kipt::VersionNumber \n$::kipt::Web" $off_x 230]
    # next time to be removed...
    # new_gid = -1 for gid v. 11.1.5d or below
    # dev_kit's splash has been modified and some changes are already there...
    set new_splash [::GidUtils::VersionCmp 11.1.6d]
    if {[winfo exists .splash.lv] && ($new_splash < 0)} {
	.splash.lv configure -font $fnt -background white -foreground black \
	    -relief solid -borderwidth 1 -padx 12 -pady 3
	update
    }
    GiD_Set SplashWindow $prev_splash_state
}

proc kipt::BeforeMeshGeneration {elementsize} {

    set ndime 3D
    # Get the spatial dimension
    catch {set ndime [::xmlutils::setXml {GeneralApplicationData//c.Domain//i.SpatialDimension} dv]}

    # Reset Automatic Conditions from previous executions
    #GiD_Process Mescape Meshing MeshCriteria DefaultMesh Lines 1:end
    #GiD_Process Mescape Meshing MeshCriteria DefaultMesh Surfaces 1:end escape
    #end common

    if {$ndime =="2D"} {
	# Align the normal
	::wkcf::AlignLineNormals Outwards
	# Reset Automatic Conditions from previous executions
	set entitytype "line"
	# Automatic Kratos Group for Boundary Condition
	set groupid "-AKGSkinMesh2D"
	::wkcf::CleanAutomaticConditionGroupGiD $entitytype $groupid
	# Find boundaries
	set blinelist [::wkcf::FindBoundaries $entitytype]
	# Automatically meshing all the boundary lines
	#GiD_Process Mescape Meshing MeshCriteria Mesh Lines {*}$blinelist escape
	# Assign the boundary condition
	::wkcf::AssignConditionToGroupGID $entitytype $blinelist $groupid
	#
	::wkcf::AssignGeometricalEntitiesToSkinSphere2D
	# Special case of DEM
	::wkcf::AssignSpecialBoundaries $ndime $blinelist

    } elseif {$ndime=="3D"} {

	# Align the normal
	::wkcf::AlignSurfNormals Outwards
	# Reset Automatic Conditions from previous executions
	set entitytype "surface"
	# Automatic Kratos Group for Boundary Condition
	set groupid "-AKGSkinMesh3D"
	::wkcf::CleanAutomaticConditionGroupGiD $entitytype $groupid
	# Find boundaries
	set bsurfacelist [::wkcf::FindBoundariesOfNonSphericElements $entitytype]
	set allsurfacelist [::wkcf::FindAllSurfacesOfNonSphericElements $entitytype]
	#
	::wkcf::AssignGeometricalEntitiesToSkinSphere3D
	# Get the surface type list
	lassign [::wkcf::GetSurfaceTypeList $bsurfacelist] tetrasurf hexasurf
	# Triangle
	if {[llength $tetrasurf]} {
	    # Assign the triangle element type
	    GiD_Process Mescape Meshing ElemType Triangle $tetrasurf escape
	    # Automatically meshing all the boundary surfaces
	    GiD_Process Mescape Meshing MeshCriteria Mesh Surfaces {*}$tetrasurf escape
	}
	# Quadrilateral
	if {[llength $hexasurf]} {
	    # Assign the quadrilateral element type
	    GiD_Process Mescape Meshing ElemType Quadrilateral $hexasurf escape
	    # Automatically meshing all the boundary surfaces
	    GiD_Process Mescape Meshing MeshCriteria Mesh Surfaces {*}$hexasurf escape
	}
	::wkcf::AssignConditionToGroupGID $entitytype $bsurfacelist $groupid

	# Special case of DEM
	::wkcf::AssignSpecialBoundaries $ndime $allsurfacelist
	::wkcf::ForceTheMeshingOfDEMFEMWallGroups
	::wkcf::ForceTheMeshingOfDEMInletGroups
    }
}

proc kipt::InitGIDProject {dir scripts_dir} {

    global KPriv
    global GidPriv
    # Set dir to a global variable
    set KPriv(dir) $dir
    set KPriv(problemTypeDir) $dir

    set KPriv(CurvesModule) 0  ;# For activating the Curves Module [Disabled -> 0 | Enabled -> 1]
    set KPriv(RDConfig) 0 ;# For release/debug options [Release =>1|Debug => 0]

    # Read kratos.xml file
    if {[info procs ::ReadProblemtypeXml] == ""} {
	WarnWin [=  "This GiD version is too old, get the latest version available"]
	return 1
    }

    set xmlfile [file join $dir kratos.xml]
    set data [ReadProblemtypeXml $xmlfile Infoproblemtype {Version MinimumGiDVersion}]
    if {$data == ""} {
	if {[file exists $xmlfile]} {
	    WarnWinText [= "Could not read file %s" $xmlfile]
	} else {
	    WarnWinText [= "Configuration file %s not found" $xmlfile]
	}
	return 1
    }
    array set problemtype_local $data
    set kipt::VersionNumber $problemtype_local(Version)

    # Check the required GiD version
    if { [::GidUtils::VersionCmp $problemtype_local(MinimumGiDVersion)]  < 0 } {
	WarnWinText [= "Error: This interface requires GiD %s or later" $problemtype_local(MinimumGiDVersion)].
    }

    kipt::Splash

    # Init some xml global variables
    ::kipt::InitGlobalXMLVariables
    ::kipt::LoadSourceFiles $dir $scripts_dir

    # Variable to control the group deletion (when exists from the problem type)
    set KPriv(Groups,DeleteGroup) 1

    set ptypeName [lindex [split $KPriv(problemTypeDir) "/"] end]
    set KPriv(pTypeName) [string map {".gid" ""} $ptypeName]

    # Set images directory
    if {[gid_themes::GetCurrentTheme] == "GiD_black"} {
	set KPriv(imagesdir) [file join images Dark]
    } else {
	set KPriv(imagesdir) [file join images Classic]
    }

    # Change system menu
    # Preprocess
    # scripts/menus.tcl
    ::kmtb::ChangePreprocessMenu $dir

    if {$GidPriv(Language) == "es"} {
	WarnWin [= "You are currently working with the Spanish version of GiD. The DEMPack Problem\
	Types need the GiD English version in order to fully work. GiD will now automatically shift\
	to English. In the Preferences Window that will open next, click on 'Close' if you agree with\
	this change or go to 'Language' to set GiD back into Spanish. In the latter case, DEMPack will\
	remain unloaded."]
	MsgcatSetLocale en
	UpdateWidgetsLanguage
	set GidPriv(Language) en
    }

    # Create the process toolbar
    ::kmtb::CreatePreprocessModelTBar

    ChangeWindowTitle

    #::kwiz::CreateDEMWizard

    # Maintain the problem type
    GiD_Set MaintainProblemTypeInNew 1

    return 0
}

proc kipt::InitGlobalXMLVariables {} {
    global KPriv

    # List of node Id's
    set KPriv(groupsId) {}

    # List of material Id´s
    set KPriv(materialsId) {}

    # Xml root
    set KPriv(xml) ""
    set KPriv(encrXml) ""
    set KPriv(xmlDoc) ""

    # Xml materials root
    set KPriv(xmlMat) ""
    set KPriv(encrXmlMat) ""
    set KPriv(xmlDocMat) ""

    # Xml functions
    set KPriv(xmlFun) ""
    set KPriv(encrXmlFun) ""
    set KPriv(xmlDocFun) ""

    # kratos_key_words.xml
    set KPriv(xmlKKW) ""
    set KPriv(xmlDocKKW) ""

    # kratos.ini
    set KPriv(xmlIni) ""
    set KPriv(xmlDocIni) ""
}

proc kipt::EndGIDProject {} {
    global KPriv
    # Destroy all pre/post open window

    # For group editor
    set w .gid.kegroups
    if {[winfo exists $w]} {
	destroy $w
    }

    set w .gid.kmprops
    if {[winfo exists $w]} {
	::KMProps::CloseWindowInside $w
	if {[winfo exists $w]} { destroy $w }
    }

    # Validation window
    set w .gid.modelvalidation
    if {[winfo exists $w]} {
	::KMValid::CreateReportWindowbClose $w
    }

    # Close Project Settings Window if it exists
    set w .gid.settingWin
    if {[winfo exists $w]} {
	::kps::WindowbClose $w
    }

    # ********************************
    #     End the bitmaps
    # ********************************
    # Preprocess
    ::kmtb::EndCreatePreprocessTBar

    # Postprocess


    # clean tDom object if exists
    catch { [$KPriv(xml) delete] }
    catch { [$KPriv(xmlMat) delete] }
    if { [info exists KPriv(xmlDoc) ] } {
	$KPriv(xmlDoc) delete
    }
    if { [info exists KPriv(xml) ] } {
	unset KPriv(xml)
    }

    # Unset the problem type global variables
    kipt::UnsetGlobalVars

    ChangeWindowTitle
}

proc kipt::LoadSourceFiles {dir scripts_dir} {

    # Load the application scripts for Kratos applications
    global KPriv

    # Load some packages
    set lib_paths [list]
    set lib_filenames [list]

    # For scripts directory
    lappend lib_paths [file join $scripts_dir scripts]
    lappend lib_filenames {files.tcl winutils.tcl menus.tcl utils.tcl stringutils.tcl \
	modelvalidation.tcl projectSettings.tcl}

    # For xml libs
    lappend lib_paths [file join $scripts_dir scripts libs xml]
    lappend lib_filenames {xmlutils.tcl xpathq.tcl}

    if { [kipt::CurvesModule] } {

		# For Curves, graphics and tables
		lappend lib_paths [file join $scripts_dir scripts libs curves]
		lappend lib_filenames {curves.tcl tables.tcl}
		lappend lib_paths [file join $scripts_dir scripts libs graphics]
		lappend lib_filenames {plotgraph.tcl}
    }

    lappend lib_paths [file join $scripts_dir scripts libs graphics]
    lappend lib_filenames {plotgraph.tcl}

    # For Wizard library
    lappend lib_paths [file join $scripts_dir scripts libs wizard]
    lappend lib_filenames {snitwiz.tcl}

    # For wcb library
    lappend lib_paths [file join $scripts_dir scripts libs wcb]
    lappend lib_filenames {wcb.tcl}

    # For write calculation file
    lappend lib_paths [file join $scripts_dir scripts libs wkcf]
    lappend lib_filenames {wkcf.tcl wkcfutils.tcl wkcffluid.tcl wkcfstructuralanalysis.tcl wkcfgroups.tcl wkcfconvectiondiffusion.tcl wkcfdem.tcl}

    # Load kegroups
    lappend lib_paths [file join $scripts_dir scripts kegroups]
    lappend lib_filenames {kegroups.tcl kGroupEntities.tcl}

    # Load KMProps
    lappend lib_paths [file join $scripts_dir scripts kmprops]
    lappend lib_filenames {kmprops.tcl kmpropswin.tcl kmpropsfwg.tcl kmpropstree.tcl \
	kmpropsgroups.tcl kmpropscbwd.tcl kmaterials.tcl kFunctions.tcl}

    # Source them
    foreach lib_path $lib_paths filenames $lib_filenames {
		foreach filename $filenames {
			set full_filename [file join $lib_path $filename]
			# msg "sourcing $full_filename"
			if { [catch {source $full_filename} msg ] } {
			WarnWin [= "Error reading file '%s': %s" $full_filename $msg].
			return 1
			}
		}
    }

    # For Wizards
    set lib_paths [list ]
    set lib_filenames [list ]

	if {$KPriv(what_dempack_package) eq "C-DEMPack"} {
		# Compressiontest wizard
		lappend lib_paths [file join $dir wizards compressiontest]
		lappend lib_filenames {compressiontest_defvals.tcl compressiontest_body.tcl}
	}

    # source them
    foreach lib_path $lib_paths filenames $lib_filenames {
		foreach filename $filenames {
			set full_filename [file join $lib_path $filename]
			if { [catch {source $full_filename} msg ] } {
			WarnWin [= "Error reading file '%s': %s" $full_filename $msg].
			return 1
			}
		}
    }

    package require snitwiz
    package require wcb
    package require KEGroups

    #kike: en kGroupEntities.tcl habia un "package provide KEGroups 1.0" !!
    #      con lo cual este package require KEGroups no hara un source de kegroups.tcl
    #      imagino que estaba mal kGroupEntities.tcl, lo he cambiado a "package require KEGroups"

    return 0
}

proc kipt::InitGIDPostProcess {} {
    set ::KMProps::RestoreWinFromPost 0
    if {[info exists ::KMProps::Layout]} {
	if {($::KMProps::Layout eq "INSIDE_LEFT") ||($::KMProps::Layout eq "INSIDE_RIGHT")} {
	    set w .gid.kmprops
	    if {[winfo exists $w]} {
		destroy $w
		set ::KMProps::RestoreWinFromPost 1
	    }
	}
    }

    # Get application type
    # Structural analysis
    set cxpath "GeneralApplicationData//c.ApplicationTypes//i.StructuralAnalysis"
    set cproperty dv
    set StructuralAnalysis [::xmlutils::setXml $cxpath $cproperty]

    # Fuild application
    set cxpath "GeneralApplicationData//c.ApplicationTypes//i.Fluid"
    set cproperty "dv"
    set FluidApplication [::xmlutils::setXml $cxpath $cproperty]

    set appid ""
    if {$FluidApplication =="Yes"} {
	set appid Fluid
    } elseif {$StructuralAnalysis=="Yes"} {
	set appid StructuralAnalysis
    }

    if {$appid !=""} {
	# Get the result type
	set cprop GiDMultiFileFlag
	set cxpath "$appid//c.Results//c.GiDOptions//i.[list ${cprop}]"
	set cproperty dv
	set rtype [::xmlutils::setXml $cxpath $cproperty]

	# Get the GiD post mode
	set cprop GiDPostMode
	set cxpath "$appid//c.Results//c.GiDOptions//i.[list ${cprop}]"
	set cproperty dv
	set pmode [::xmlutils::setXml $cxpath $cproperty]

	set existfiles [::KUtils::ReadResultsFromFiles $appid $rtype $pmode CheckRFiles]
	if {!$existfiles} {
	    WarnWin [= "The simulation is not calculated yet or is currently being calculated"].
	    return ""
	} else {
	    # Try to read the result files
	    set ok [::KUtils::ReadResultsFromFiles $appid $rtype $pmode ReadRFiles]
	}
    }
}

proc kipt::SelectGIDBatFile {directory basename} {
    set batfilename ""
    set args ""

    # Get application type
    # Structural analysis
    set cxpath "GeneralApplicationData//c.ApplicationTypes//i.StructuralAnalysis"
    set cproperty "dv"
    set StructuralAnalysis [::xmlutils::setXml $cxpath $cproperty]

    # WarnWinText "StructuralAnalysis:$StructuralAnalysis"

    # Fluid application
    set cxpath "GeneralApplicationData//c.ApplicationTypes//i.Fluid"
    set FluidApplication [::xmlutils::setXml $cxpath $cproperty]

    # DEM application
    set cxpath "GeneralApplicationData//c.ApplicationTypes//i.DEM"
    set DEMApplication [::xmlutils::setXml $cxpath $cproperty]

    # WarnWinText "FluidApplication:$FluidApplication"

    # Structural analyis
    if {$StructuralAnalysis eq "Yes"} {
	set rootid "StructuralAnalysis"
	# Kratos key word xpath
	set kxpath "Applications/$rootid"
	# Get the parallel solution type
	set cxpath "$rootid//c.SolutionStrategy//c.ParallelType//i.ParallelSolutionType"
	set ParallelSolutionType [::xmlutils::setXml $cxpath $cproperty]

	# Solution type
	set cxpath "$rootid//c.AnalysisData//i.SolutionType"
	set SolutionType [::xmlutils::setXml $cxpath $cproperty]

	if {$ParallelSolutionType eq "MPI"} {
	    if {($SolutionType =="Dynamic")||($SolutionType =="RelaxedDynamic")} {
		if {($::tcl_platform(os) eq "Linux")} {
		    set batfilename "kratos-structural-mpi.unix.bat"
		} else {
		    # set batfilename "kratos-structural-mpi.win.bat"
		}
	    } elseif {$SolutionType =="Static"} {
		if {($::tcl_platform(os) eq "Linux")} {
		    set batfilename "kratos-structural-mpi.unix.bat"
		} else {
		    # set batfilename "kratos-structural-mpi.win.bat"
		}
	    }

	    #  Get the number of processors
	    set cxpath "$rootid//c.SolutionStrategy//c.ParallelType//i.MPINumberOfProcessors"
	    set MPINumberOfProcessors [::xmlutils::setXml $cxpath $cproperty]
	    if {$MPINumberOfProcessors>0} {
		# Calculate arguments
		set args "$MPINumberOfProcessors"
	    }
	} else {
	    # OpenMP
	    #  Get the number of threads
	    set cxpath "$rootid//c.SolutionStrategy//c.ParallelType//i.OpenMPNumberOfThreads"
	    set OpenMPNumberOfThreads [::xmlutils::setXml $cxpath $cproperty]
	    if {$OpenMPNumberOfThreads>0} {
		# Calculate arguments
		set args "$OpenMPNumberOfThreads"
	    }
	    if {($SolutionType =="Dynamic")||($SolutionType =="Quasi-Static")||($SolutionType =="Pseudo-Dynamic")} {
		if {($::tcl_platform(os) eq "Linux")} {
		    set batfilename "kratos-structural-openmp.unix.bat"
		} else {
		    set batfilename "kratos-structural-openmp.win.bat"
		}
	    } elseif {$SolutionType =="Static"} {
		if {($::tcl_platform(os) eq "Linux")} {
		    set batfilename "kratos-structural-openmp.unix.bat"
		} else {
		    set batfilename "kratos-structural-openmp.win.bat"
		}
	    }
	}
    }

    # Fluid application
    if {$FluidApplication eq "Yes" && $DEMApplication eq "Yes"} {
	set rootid "Fluid"
	# Kratos key word xpath
	set kxpath "Applications/$rootid"
	# Get the parallel solution type
	set cxpath "$rootid//c.SolutionStrategy//c.ParallelType//i.ParallelSolutionType"
	set ParallelSolutionType [::xmlutils::setXml $cxpath $cproperty]

	# Free surface
	set cxpath "$rootid//c.AnalysisData//i.FreeSurface"
	set FreeSurface [::xmlutils::setXml $cxpath $cproperty]

	# Get the fluid approach
	set cxpath "$rootid//c.AnalysisData//i.FluidApproach"
	set FluidApproach [::xmlutils::setXml $cxpath $cproperty]

	# Solver type for free surface
	set cxpath "$rootid//c.AnalysisData//i.SolverTypeFreeSurf"
	set SolverTypeFreeSurf [::xmlutils::setXml $cxpath $cproperty]
	# WarnWinText "SolverTypeFreeSurf:$SolverTypeFreeSurf"

	if {$ParallelSolutionType eq "MPI"} {
	    if {($::tcl_platform(os) eq "Linux")} {
		set batfilename "kratos-mpi.unix.bat"
		#  Get the number of processors
		set cxpath "$rootid//c.SolutionStrategy//c.ParallelType//i.MPINumberOfProcessors"
		set MPINumberOfProcessors [::xmlutils::setXml $cxpath $cproperty]
		if {$MPINumberOfProcessors>0} {
		    # Calculate arguments
		    set args "$MPINumberOfProcessors"
		}
	    }
	} else {
	    # OpenMP
	    #  Get the number of threads
	    set cxpath "GeneralApplicationData//c.SimulationOptions//i.OpenMPNumberOfThreads"
	    set OpenMPNumberOfThreads [::xmlutils::setXml $cxpath $cproperty]
	    if {$OpenMPNumberOfThreads>0} {
		# Calculate arguments
		set args "$OpenMPNumberOfThreads"
	    }
	    if {($FreeSurface eq "Yes") && ($SolverTypeFreeSurf eq "LevelSet")} {
		if {($::tcl_platform(os) eq "Linux")} {
		    set batfilename "kratos-openmplevelset.unix.bat"
		} else {
		    set batfilename "kratos-openmplevelset.win.bat"
		}
	    } else {
		if {($::tcl_platform(os) eq "Linux")} {
		    set batfilename "kratos.unix.bat"
		} else {
		    set batfilename "kratos.win.bat"
		}
	    }
	}
    }


    # DEM application
    if {$DEMApplication eq "Yes" && $FluidApplication eq "No"} {
	set rootid "DEM"
	# Kratos key word xpath
	set kxpath "Applications/$rootid"
	# Get the parallel solution type
	set cxpath "$rootid//c.DEM-SolutionStrategy//c.DEM-ParallelType//i.ParallelType"
	set ParallelSolutionType [::xmlutils::setXml $cxpath $cproperty]

	if {$ParallelSolutionType eq "MPI"} {
	    if {($::tcl_platform(os) eq "Linux")} {
		set batfilename "kratos-mpi.unix.bat"
		#  Get the number of processors
		set cxpath "$rootid//c.DEM-SolutionStrategy//c.DEM-ParallelType//i.NumberOfThreads"
		set MPINumberOfProcessors [::xmlutils::setXml $cxpath $cproperty]
		if {$MPINumberOfProcessors>0} {
		    # Calculate arguments
		    set args "$MPINumberOfProcessors"
		}
	    } else {
		#Windows or Mac
		WarnWinText [= "MPI is not available for Windows yet. You can use OpenMP instead."]
	    }
	} else {
	    # OpenMP
	    #  Get the number of threads
	    set cxpath "$rootid//c.DEM-SolutionStrategy//c.DEM-ParallelType//i.NumberOfThreads"
	    set OpenMPNumberOfThreads [::xmlutils::setXml $cxpath $cproperty]
	    if {$OpenMPNumberOfThreads>0} {
		# Calculate arguments
		set args "$OpenMPNumberOfThreads"
	    }
	    if {($::tcl_platform(os) eq "Linux")} {
		set batfilename "kratos.unix.bat"
	    } else {
		set batfilename "kratos.win.bat"
	    }
	}
    }

    set ret "$batfilename $args"
    if {$batfilename != ""} {
		return $ret
    } else {
		return ""
    }
}

proc kipt::ChangedLanguage {language} {
    global KPriv
    #set dir [GiD_Info problemtypepath]
    set dir $KPriv(dir)
    # Preprocess
    ::kmtb::ChangePreprocessMenu $dir
    # Postprocess
}

proc kipt::CurvesModule {} {
    global KPriv
    if {[info exists KPriv(CurvesModule)]} {
	return $KPriv(CurvesModule)
    }
    return 0
}

proc kipt::UnsetGlobalVars {} {
    global KData KPriv
    foreach arrid [list KData KPriv] {
		if {[info exists $arrid]} {
		        unset $arrid
		}
    }
}

########################################################################
#    KRATOS: GiD interface for the Kratos problem type
########################################################################
#
#    NAME: Kratos.tcl
#
#    PURPOSE: Init script for Kratos problem type
#
#    AUTHOR : G. Socorro
#
#    CREATED AT: 01/11/09
#
#    HISTORY: 
# 
#     1,6- 22/07/12- G. Socorro, modify the BeforeMeshGeneration to automatically mesh de boundary lines/surfaces when use Is-Slip BC
#     1.5- 09/05/12- G. Socorro, use conditions only in the fluid application
#     1.4- 07/05/12- G. Socorro, update the proc BeforeMeshGeneration to write the Condition2D and Condition3D properties
#     1.3- 04/05/12- G. Socorro, update the proc BeforeDeleteGroup
#     1.2- 03/05/12- J. Garate, proc EndGIDPostProces
#     1.1- 03/05/12- G. Socorro, load the new gid_groups package
#     1.0- 23/03/12- J.G�rate, No se puede cancelar un renombre
#     0.9- 20/03/12- J.G�rate, Group Event Procedures
#     0.8- 16/11/11- G. Socorro, add the global variable KPriv(SRCConfig) to use TCL or TBE file distribution
#     0.7- 22/06/11- G. Socorro, add the proc KLoadTBEFiles to load the sources from TBE files
#     0.6- 03/02/11- G. Socorro, update the procedure UnsetGlobalVars
#     0.5- 08/09/10- G. Socorro, add the event InitGIDPostProcess to read Kratos result files when pass from preproces to postprocess
#     0.4- 03/09/10- G. Socorro, add BeforeMeshGeneration option to modify the normal to the line and surfaces
#     0.3- 01/02/10- G. Socorro, add a new global procedure msg to call WarnWinText
#     0.2- 24/12/09- G. Socorro, add new xmlutils, xpathq and wkcf namespaces
#     0.1- 01/11/09- G. Socorro, create the base source file
#
########################################################################

proc ReadSomePTData {dir} {
	
    global ProgramName VersionNumber MinimumGiDVersion

    dom parse [tDOM::xmlReadFile [file join $dir kratos.xml]] doc
    
    set ProgramName [$doc selectNodes string(Infoproblemtype/Program/Name)]
    set VersionNumber [$doc selectNodes string(Infoproblemtype/Program/Version)]
    set MinimumGiDVersion [$doc selectNodes string(Infoproblemtype/Program/MinimumGiDVersion)]

    $doc delete
}

proc CheckRequiredGiDVersion {VersionRequired} {
    
    set comp -1
    catch { 
	set comp [::GidUtils::VersionCmp $VersionRequired]
    }
    if { $comp < 0 } {
	msg "Error: This interface requires GiD $VersionRequired or later"
    }
}

proc UnsetGlobalVars {} {
	
    global KData KPriv MinimumGiDVersion
    global VersionNumber ProgramName

    foreach arrid [list KData KPriv VersionNumber ProgramName MinimumGiDVersion] {
	if {[info exists $arrid]} {
	    unset $arrid
	}
    }
}

proc LoadGIDProject {filename} {

    ::kfiles::LoadSPD $filename
}

proc SaveGIDProject {filename} {
    
    ::kfiles::SaveSPD $filename
}

proc AfterTransformProblemType { file oldproblemtype newproblemtype } {
	
	set name [lindex [split $file "/"] end]
	#msg "${file}/${name}.spd"
	LoadGIDProject "${file}.gid/${name}.spd"
	
	return 0
}

#proc BeforeTransformProblemType { file oldproblemtype newproblemtype } {
	
	#msg "before transform"
	#global KPriv
	#msg "$KPriv(dir)"
	#set path [GiD_Info problemtypepath]
	#set name [lindex [split $file "/"] end]
	#
	#msg "after transform path:$path name:$name\n$file $oldproblemtype $newproblemtype"
	#
	##LoadGIDProject "$file/kratos_default.spd"
	#
	#
	##
	###Transforma el spd si son versiones distintas
	 #::xmlutils::checkSpdVersion
	 #
	 #return -cancel-
#}

proc InitGIDProject { dir } {
	
    global KData KPriv
    global VersionNumber ProgramName MinimumGiDVersion

    # Unset global variables
    UnsetGlobalVars
    
    # Read kratos.xml file
    ReadSomePTData $dir
    
    # WarnWinText "VersionNumber:$VersionNumber ProgramName:$ProgramName MinimumGiDVersion:$MinimumGiDVersion"
    # Check the required GiD version
    set VersionRequired "$MinimumGiDVersion"
    CheckRequiredGiDVersion $VersionRequired
    
    # Init packages
    SRC gid_groups_public.tcl
    gid_groups_conds::init_package
    
    # For release/debug options [Release =>1|Debug => 0]
    set KPriv(RDConfig) 0
    # For distribution srctcl/srctbe options [srctbe =>1|srctcl => 0]
    set KPriv(SRCConfig) 0

    # Load the application scripts 
    if {!$KPriv(SRCConfig)} {
	# For scripts directory
	set scriptspath "$dir/scripts/"
	if { [catch {source $scriptspath/initptype.tcl}] } {
	    return 0
	} else {
	    # Init some xml global variables
	    ::kipt::InitGlobalXMLVariables
	    
	    ::kipt::LoadSourceFiles $dir
	}
    } else {
	# Load tbe files
	KLoadTBEFiles $dir

	# Init some xml global variables
	::kipt::InitGlobalXMLVariables
    }

    # Init problem type
    ::kipt::InitPType $dir
    
    
}

proc EndGIDProject {} {

    # Free problem type
    ::kipt::FreePType
}

proc BeforeWriteCalcFileGIDProject { file } {
 
    #return 1 ;
    return [::KMValid::ValidateModel]

}

proc AfterWriteCalcFileGIDProject {filename errorflag } {
    # WarnWinText "AfterWriteCalcFileGIDProject\n filename:$filename\n errorflag:$errorflag"
    
    set ret 1

    # Try to write the Kratos input data file
    set err [catch { ::wkcf::WriteCalculationFiles $filename} ret]
    if { $err } {
	snit_messageBox -parent .gid -message \
	    [= "Error when preparing data for analysis (%s)" $ret]
	return "-cancel-"
    }
    return $ret
}

proc msg {mesage} {
    
    WarnWinText $mesage
}

proc wa {mesage} {
    
    WarnWinText $mesage
}

proc BeforeMeshGeneration {elementsize} { 

    set ndime "3D"
    # Get the spatial dimension
    set cxpath "GeneralApplicationData//c.Domain//i.SpatialDimension"
    set cproperty "dv"
    catch { set ndime [::xmlutils::setXml $cxpath $cproperty] }
    
    # Get application type
  
    # Fuild application
    set cxpath "GeneralApplicationData//c.ApplicationTypes//i.Fluid"
    set cproperty "dv"
    set FluidApplication [::xmlutils::setXml $cxpath $cproperty]

    if {$FluidApplication eq "Yes"} {
	if {$ndime =="2D"} {
	    
	    # Align the normal
	    ::wkcf::AlignLineNormals Outwards 
	    
	    # Reset Automatic Conditions from previous executions 
	    set what "UseAllWhereField"
	    set entitytype "line"
	    set groupid "-@kratos@b2d"
	    # Old groups
	    set fieldname "groupid"
	    ::wkcf::CleanAutomaticConditionGroup $what $entitytype $fieldname $groupid
	    # New groups
	    set fieldname "name"
	    ::wkcf::CleanAutomaticConditionGroup $what $entitytype $fieldname $groupid
	    
	    # Check for use Is-Slip BC
	    set issliplist [::KMValid::SlipNoSlipList "slip"]
	    if {[llength $issliplist]} {
		# Find boundaries
		set blinelist [::wkcf::FindBoundaries $entitytype]
		# wa "belist:$blinelist"
		
		# Automatically meshing all the boundary lines
		GiD_Process Mescape Meshing MeshCriteria Mesh Lines $blinelist escape 
		
		# Assign the boundary condition
		::wkcf::AssignConditionToGroup $entitytype $blinelist $groupid
	    }

	} elseif {$ndime =="3D"} {
	    
	    # Align the normal
	    ::wkcf::AlignSurfNormals Outwards
	    
	    # Reset Automatic Conditions from previous executions 
	    set what "UseAllWhereField"
	    set entitytype "surface"
	    set groupid "-@kratos@b3d"
	    # Old groups
	    set fieldname "groupid"
	    ::wkcf::CleanAutomaticConditionGroup $what $entitytype $fieldname $groupid
	    # New groups
	    set fieldname "name"
	    ::wkcf::CleanAutomaticConditionGroup $what $entitytype $fieldname $groupid
	    
	    # Check for use Is-Slip BC
	    set issliplist [::KMValid::SlipNoSlipList "slip"]
	    if {[llength $issliplist]} {
		# Find boundaries
		set bsurfacelist [::wkcf::FindBoundaries $entitytype]
		# WarnWinText "bsurfacelist:$bsurfacelist"
		
		# Automatically meshing all the boundary surfaces
		GiD_Process Mescape Meshing MeshCriteria Mesh Surfaces $bsurfacelist escape 
		
		# Assign the boundary condition
		::wkcf::AssignConditionToGroup $entitytype $bsurfacelist $groupid
	    }
	}
    }

}

proc InitGIDPostProcess {} { 

    set ::KMProps::RestoreWinFromPost 0
    if {[info exists ::KMProps::Layout]} {
	if {($::KMProps::Layout eq "INSIDE_LEFT") ||($::KMProps::Layout eq "INSIDE_RIGHT")} {
	    set w ".gid.kmprops" 
	    if {[winfo exists $w]} {
		destroy $w
		set ::KMProps::RestoreWinFromPost 1
	    }
    	}
    }

    # Get application type
    # Structural analysis
    set cxpath "GeneralApplicationData//c.ApplicationTypes//i.StructuralAnalysis"
    set cproperty "dv"
    set StructuralAnalysis [::xmlutils::setXml $cxpath $cproperty]
    
    # WarnWinText "StructuralAnalysis:$StructuralAnalysis"

    # Fuild application
    set cxpath "GeneralApplicationData//c.ApplicationTypes//i.Fluid"
    set cproperty "dv"
    set FluidApplication [::xmlutils::setXml $cxpath $cproperty]

    # WarnWinText "FluidApplication:$FluidApplication"
    set appid ""
    if {$FluidApplication =="Yes"} {
	set appid "Fluid"
    } elseif {$StructuralAnalysis=="Yes"} {
	set appid "StructuralAnalysis"
    }

    if {$appid !=""} {
	# Get the result type
	set cprop "GiDMultiFileFlag"
	set cxpath "$appid//c.Results//c.GiDOptions//i.${cprop}"
     	set cproperty "dv"
     	set rtype [::xmlutils::setXml $cxpath $cproperty]
	# WarnWinText "rtype:$rtype"
	
	# Get the GiD post mode
	set cprop "GiDPostMode"
	set cxpath "$appid//c.Results//c.GiDOptions//i.${cprop}"
     	set cproperty "dv"
     	set pmode [::xmlutils::setXml $cxpath $cproperty]
	# WarnWinText "pmode:$pmode"

	set existfiles [::KUtils::ReadResultsFromFiles $appid $rtype $pmode "CheckRFiles"]
	if {!$existfiles} {
	    WarnWin [= "The simulation is not calculated yet or is currently being calculated"].
	    return ""
	} else {
	    # Try to read the result files
	    set ok [::KUtils::ReadResultsFromFiles $appid $rtype $pmode "ReadRFiles"]
	}
    } 
}

proc KLoadTBEFiles {dir} {
    
    # For scripts directory
    set scriptspath "$dir/scripts/"
    cd $scriptspath
      
    set dirlist [glob *]
    # WarnWinText "dirlist:$dirlist\n\n"
    foreach cdir $dirlist {
	if {[file isdirectory $cdir]} {
	    cd $cdir
	    set tbelist ""
	    catch { set tbelist [glob *] }
	    # WarnWinText "cdir:$cdir => tbelist:$tbelist\n"
	    if {[llength $tbelist]} {
		foreach tbe_level1 $tbelist {
		    # WarnWinText "Current tbe_level1:$tbe_level1"
		    if {[file isdirectory $tbe_level1]} {
			cd $tbe_level1
			set tbe_l2list ""
			catch { set tbe_l2list [glob *] }
			# WarnWinText "tbe_l2list:$tbe_l2list"
			if {[llength $tbe_l2list]} {
			    foreach tbe_level2 $tbe_l2list {
				if {[file extension $tbe_level2]==".tbe"} {
				    # WarnWinText "current l2 tbe_l2list:$tbe_level2"
				    loadtbefile $tbe_level2		    
				}
			    }
			}
			cd ..
		    } else {
			if {[file extension $tbe_level1]==".tbe"} {
			    # WarnWinText "current level1:$tbe_level1\n"
			    loadtbefile $tbe_level1
			}
		    }
		}
	    }
	    cd ..
	} else {
	    if {[file extension $cdir]==".tbe"} {
		# WarnWinText "current cdir:$cdir\n"
		loadtbefile $cdir
	    }
	}
    }
} 

proc BeforeDeleteGroup { name } {
    # wa "delete name:$name"
    
    set DeleteGroup "Delete" 
    if {[info exists ::KPriv(Groups,DeleteGroup)]} {
	if {$::KPriv(Groups,DeleteGroup)} {
	    set DeleteGroup [::KEGroups::BorraGrupo $name]
	} 
    }
    if { $DeleteGroup eq "-cancel-" } {
	return $DeleteGroup
    }
}

proc AfterCreateGroup { name } {
     #wa "name:$name"
}

proc AfterRenameGroup { oldname newname } {
    #wa "oldname:$oldname newname:$newname"
    ::KEGroups::RenombraGrupo $oldname $newname 1
    #Si se renombra un grupo, no nos queda otra... no se puede impedir.
    #return
}

proc EndGIDPostProcess { } {
    
    # Try to restore the properties window
    if {[info exists ::KMProps::RestoreWinFromPost]} {
	if {$::KMProps::RestoreWinFromPost} {
	    ::KMProps::StartBaseWindow
	    
    	}
    }
}
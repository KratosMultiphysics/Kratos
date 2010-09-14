########################################################################
#    KRATOS: GiD interface for the Kratos problem type
########################################################################
#
#    NAME: Blade-Pre.tcl
#
#    PURPOSE: Init script for Kratos problem type
#
#    AUTHOR : G. Socorro
#
#    CREATED AT: 01/11/09
#
#    LAST MODIFICATION : add the event InitGIDPostProcess to read Kratos result files when pass from preproces to postprocess
#
#    VERSION : 0.5
#
#    HISTORY: 
#
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
	WarnWin [= "Error: This interface requires GiD %s or later" $VersionRequired].
    }
}

proc UnsetGlobalVars {} {
	
    global KData KPriv MinimumGiDVersion
    global VersionNumber ProgramName

    catch { unset KData }
    catch { unset KPriv }
    catch { unset VersionNumber }
    catch { unset ProgramName }
    catch { unset MinimumGiDVersion }

}
proc LoadGIDProject {filename} {
    
    ::kfiles::LoadSPD $filename
}

proc SaveGIDProject {filename} {
	
	::kfiles::SaveSPD $filename
}

proc AfterTransformProblemType { file oldproblemtype newproblemtype } {
	
	set name [lindex [split $file "/"] end]
	msg "${file}/${name}.spd"
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
    
    # Load init problem type script
    source $dir/scripts/initptype.tcl  

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

proc BeforeMeshGeneration {elementsize} { 

    set ndime "3D"
    # Get the spatial dimension
    set cxpath "GeneralApplicationData//c.Domain//i.SpatialDimension"
    set cproperty "dv"
    catch { set ndime [::xmlutils::setXml $cxpath $cproperty] }
    
    if {$ndime =="2D"} {
	::wkcf::AlignLineNormals Outwards 	
    } elseif {$ndime =="3D"} {
	::wkcf::AlignSurfNormals Outwards
    }

}

proc InitGIDPostProcess {} { 
 
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
 
 
 
 
 
 
 
 
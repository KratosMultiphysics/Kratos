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
#    LAST MODIFICATION : add a new global procedure msg to call WarnWinText
#
#    VERSION : 0.3
#
#    HISTORY: 
#
#     0.3- 01/02/10-G. Socorro, add a new global procedure msg to call WarnWinText
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

 
 
 
 
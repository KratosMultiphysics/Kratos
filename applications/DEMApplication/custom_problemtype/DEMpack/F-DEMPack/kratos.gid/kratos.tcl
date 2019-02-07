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
#########################################################
### START GiD-Tcl Events (procedures raised from GiD) ###
#########################################################

proc LoadGIDProject {filename} {
    ::SetDempacktype
    ::kfiles::LoadSPD $filename
}

proc SaveGIDProject {filename} {
    ::kfiles::SaveSPD $filename
}

proc AfterTransformProblemType {filename oldproblemtype newproblemtype} {
    set file_tail [file tail $filename]
    set spd_filename [file join ${filename}.gid ${file_tail}.spd]
    ::kfiles::LoadSPD $spd_filename
    return 0
}

proc BeforeRunCalculation {batfilename basename dir problemtypedir gidexe args} {

    # WarnWinText  " BeforeRunCalculation => batfilename:$batfilename basename:$basename dir:$dir problemtypedir:$problemtypedir gidexe:$gidexe args:$args"

    # WarnWin "BeforeRunCalculation"

    # Start the check running process event
    ::KUtils::IsProcessRunningVar Set 1

    # set err [catch { GiD_Info Mesh Quadrilateral }]

    # if { $err } {
    #    WarnWinText [= "Cannot run the simulation. This version of the problem type does not support quadrilateral elements."]
    #    return -cancel-
    # }

    set value ok

    return $value
}

proc AfterRunCalculation {basename dir problemtypedir where error errorfilename} {

    # WarnWinText "AfterRunCalculation => basename:$basename dir:$dir problemtypedir:$problemtypedir where:$where error:$error errorfilename:$errorfilename"

    # Reset the check running process event
    ::KUtils::IsProcessRunningVar Set 0

    # return nowindow
    return 0
}

proc BeforeTransformProblemType {file oldproblemtype newproblemtype} {
    #set path [GiD_Info problemtypepath]
    #set name [lindex [split $file "/"] end]
    ##LoadGIDProject "$file/kratos_default.spd"
    ###Transforma el spd si son versiones distintas
    #::xmlutils::checkSpdVersion
    #return -cancel-
}

proc SetDempacktype {} {

global KPriv
set KPriv(what_dempack_package) F-DEMPack
}

proc InitGIDProject {dir} {
    global KPriv
    SetDempacktype
    # Load the application scripts
    set scripts_dir [file join $dir .. .. ]
    set tcl_filename [file join $scripts_dir scripts initptype.tcl]
    if {[catch {source $tcl_filename} msg]} {
        WarnWinText $msg
        return 1
    }
    GiD_Set CalcWithoutMesh 1
    GiDMenu::RemoveOption "Mesh" [list "Element type" "Quadrilateral"] PRE
    GiDMenu::RemoveOption Mesh [list {Quadratic type}] PRE

    
    return [::kipt::InitGIDProject $dir $scripts_dir]
}

proc EndGIDProject {} {
    return [::kipt::EndGIDProject]
}

proc BeforeWriteCalcFileGIDProject {file} {
    #return 1
    return [::KMValid::ValidateModel]
}

proc AfterWriteCalcFileGIDProject {filename errorflag} {
    set ret 1
    # Try to write the Kratos input data file
    set err [catch {lassign [::wkcf::WriteCalculationFiles $filename] fail msg} ret]
    if {$err} {
	WarnWinText [= "Error when preparing data for analysis (%s)" $ret]
	return -cancel-
    } else {
	if {$fail} {
	    WarnWinText $msg
	    return -cancel-
	}
    }
    return $ret
}

proc BeforeMeshGeneration {elementsize} {
    # ::GidUtils::DisableGraphics
    #catch to not let disabledGraphics never
    wkcf::Preprocess; # DO THIS PROPERLY
    if {[catch {kipt::BeforeMeshGeneration $elementsize} err]} {
	WarnWinText $err
    }
    # ::GidUtils::EnableGraphics
}

proc AfterMeshGeneration {fail} {
    # After Mesh Generation
    set without_window [GidUtils::AreWindowsDisabled]; #batch mode
    if {!$without_window} {
	GidUtils::DisableGraphics
    }
    if {[catch {wkcf::Elements_Substitution} msg]} {
      W "wkcf::Elements_Substitution. $msg"
    }
    #wkcf::Elements_Elimination
    #GiD_Process Mescape Meshing EditMesh DelLonelyNods Yes
    #GiD_Process Mescape Meshing MeshView
    if {!$without_window} {
	GidUtils::EnableGraphics
    }
}

proc InitGIDPostProcess {} {
    return [kipt::InitGIDPostProcess]
}

proc BeforeDeleteGroup {name} {
    global KPriv
    set DeleteGroup "Delete"
    update
    if {[info exists ::KPriv(Groups,DeleteGroup)]} {
	if {$::KPriv(Groups,DeleteGroup)} {
	    set DeleteGroup [::KEGroups::BorraGrupo $name]
	}
    }
    update
    if {$DeleteGroup eq "-cancel-"} {
	after 100 {ChangeLayers}
	return $DeleteGroup
    }
}

proc AfterRenameGroup {oldname newname} {
    # Valida para los grupos de GiD 11.1.1d
    #(work in progress, falta poco para poderlo activar) modificado kratos para que acepte cualquier nombre grupo
    #       ::KEGroups::RenombraGrupo $oldname $newname 0

    ::KEGroups::RenombraGrupo $oldname $newname 1
    #Si se renombra un grupo, no nos queda otra... no se puede impedir.
}

proc EndGIDPostProcess {} {
    # Try to restore the properties window
    if {[info exists ::KMProps::RestoreWinFromPost]} {
	if {$::KMProps::RestoreWinFromPost} {
	    ::KMProps::StartBaseWindow
	}
    }
}

proc SelectGIDBatFile {directory basename} {
    return [kipt::SelectGIDBatFile $directory $basename]
}

proc ChangedLanguage {language} {
    return [kipt::ChangedLanguage $language]
}

#####################################################
### END GID-TCL Events (procedures raised by GiD) ###
#####################################################

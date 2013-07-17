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
#     3.9- 04/07/13- A. Melendo, group name without validation
#     3.8- 25/06/13- A. Melendo, to not let disabledGraphics never
#     3.7- 19/06/13- G. Socorro, delete the event AfterRenameCondGroup, AfterCreateCondGroup and BeforeDeleteCondGroup used in the Compass library
#     3.6- 18/06/13- G. Socorro, delete the global variable KPriv(NewGiDGroups)
#     3.5- 17/06/13- G. Socorro, modify the proc BeforeMeshGeneration to delete the old Compass condition
#     3.4- 14/06/13- G. Socorro, modify the procs CheckRequiredGiDVersion
#     3.3- 18/06/13- A. Melendo, simplify the proc BeforeMeshGeneration to write conditional data condition in all applications
#     3.2- 19/04/13- G. Socorro, modify the proc BeforeMeshGeneration to write conditional data condition in the fluid application
#     3.1- 26/11/12- J. Garate,  BeforeMeshGeneration modified, support PFEM Application
#     3.0- 12/11/12- J. Garate, Minor Fixing
#     2.9- 07/11/12- J. Garate, GiD 11.1.2d is the minimum Required version for New GiD Groups
#     2.8- 29/10/12- G. Socorro, Add the proc AfterCreateGroup
#     2.7- 22/10/12- J. Garate, Function correction 
#     2.6- 17/10/12- J. Garate, Correction when transferring old Cond_Groups to GiD_Groups
#     2.5- 10/10/12- J. Garate, Adapatation for New GiD_Groups
#     2.4- 10/10/12- G. Socorro, update the proc SelectGIDBatFile to include the structural analysis bat files for MPI and OpenMP
#     2.3- 08/10/12- J. Gárate, enable/disable New GiD Groups module KPriv(NewGiDGroups) [0 -> Old Groups  | 1 -> New Groups]
#     2.2- 04/10/12- G. Socorro, update the proc SelectGIDBatFile to include the LevelSet bat file 
#     2.1- 04/10/12- G. Socorro, add the proc BeforeDeleteCondGroup,AfterCreateCondGroup and AfterRenameCondGroup
#     2.0- 01/10/12- J. Gárate, enable/disable curves module KPriv(CurvesModule) [0 -> Disable  | 1 -> Enable]
#     1.9- 28/09/12- G. Socorro, modify the event SelectGIDBatFile tp include the number of threads for the OpenMP case
#     1.8- 24/09/12- G. Socorro, modify the event SelectGIDBatFile to include the number of processors in the command line
#     1.7- 21/09/12- G. Socorro, correct a bug in the proc BeforeMeshGeneration use {*} to get all the surface identifier
#     1,6- 22/07/12- G. Socorro, modify the BeforeMeshGeneration to automatically mesh de boundary lines/surfaces when use Is-Slip BC
#     1.5- 09/05/12- G. Socorro, use conditions only in the fluid application
#     1.4- 07/05/12- G. Socorro, update the proc BeforeMeshGeneration to write the Condition2D and Condition3D properties
#     1.3- 04/05/12- G. Socorro, update the proc BeforeDeleteGroup
#     1.2- 03/05/12- J. Garate, proc EndGIDPostProces
#     1.1- 03/05/12- G. Socorro, load the new gid_groups package
#     1.0- 23/03/12- J.Gárate, No se puede cancelar un renombre
#     0.9- 20/03/12- J.Gárate, Group Event Procedures
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


#########################################################
### START GiD-Tcl Events (procedures raised from GiD) ###
#########################################################

proc LoadGIDProject {filename} {
    ::kfiles::LoadSPD $filename
}

proc SaveGIDProject {filename} {
    ::kfiles::SaveSPD $filename
}

proc AfterTransformProblemType { filename oldproblemtype newproblemtype } {
    set file_tail [file tail $filename]
    set spd_filename [file join ${filename}.gid ${file_tail}.spd]
    ::kfiles::LoadSPD $spd_filename
    return 0
}

proc BeforeTransformProblemType { file oldproblemtype newproblemtype } {      
    #set path [GiD_Info problemtypepath]
    #set name [lindex [split $file "/"] end]    
    ##LoadGIDProject "$file/kratos_default.spd"   
    ###Transforma el spd si son versiones distintas
    #::xmlutils::checkSpdVersion    
    return -cancel-
}

proc InitGIDProject { dir } {                 
    # Load the application scripts     
    set tcl_filename [file join $dir scripts initptype.tcl]
    if { [catch {source $tcl_filename} msg] } {
        WarnWinText $msg
        return 1
    }
    return [::kipt::InitGIDProject $dir]
}

proc EndGIDProject {} {
    return [::kipt::EndGIDProject]
}

proc BeforeWriteCalcFileGIDProject { file } {
    #return 1
    return [::KMValid::ValidateModel]    
}

proc AfterWriteCalcFileGIDProject {filename errorflag } {   
    set ret 1
    # Try to write the Kratos input data file
    set err [catch { ::wkcf::WriteCalculationFiles $filename } ret]
    if { $err } {
        WarnWinText [= "Error when preparing data for analysis (%s)" $ret]
        return -cancel-
    }
    return $ret
}

proc BeforeMeshGeneration { elementsize } {    
    ::GidUtils::DisableGraphics
    #catch to not let disabledGraphics never
    if { [catch { kipt::BeforeMeshGeneration $elementsize } err] } {
        WarnWinText $err
    }
    ::GidUtils::EnableGraphics       
}

proc InitGIDPostProcess {} { 
    return [kipt::InitGIDPostProcess]
}

proc BeforeDeleteGroup { name } {
    global KPriv
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

proc AfterRenameGroup { oldname newname } {
    # Valida para los grupos de GiD 11.1.1d    
    #(work in progress, falta poco para poderlo activar) modificado kratos para que acepte cualquier nombre grupo 
    #       ::KEGroups::RenombraGrupo $oldname $newname 0
    
    ::KEGroups::RenombraGrupo $oldname $newname 1
    #Si se renombra un grupo, no nos queda otra... no se puede impedir.
}


proc EndGIDPostProcess { } {    
    # Try to restore the properties window
    if {[info exists ::KMProps::RestoreWinFromPost]} {
        if {$::KMProps::RestoreWinFromPost} {
            ::KMProps::StartBaseWindow
        }
    }
}

proc SelectGIDBatFile {directory basename } {
    return [kipt::SelectGIDBatFile $directory $basename]
    
}  

proc ChangedLanguage { language } {
    return [kipt::ChangedLanguage $language]
}

#######################################################
### END GID-TCL Events (procedures raised from GiD) ###
#######################################################



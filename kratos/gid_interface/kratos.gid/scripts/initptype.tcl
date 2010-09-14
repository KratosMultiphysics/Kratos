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
#    LAST MODIFICATION: create the base source code
#
#    VERSION : 0.1
#
#    HISTORY:
#
#     0.1-01/02/10-G. Socorro, create the base source code
#
###############################################################################

# ikpt => kratos Init Problem Type 

namespace eval kipt {

}

proc kipt::InitPType { dir } {
    
    global KPriv ProgramName
    
    #Switch between RELEASE and DEBUG mode
    set KPriv(release) 1
    
    # Set dir to a global variable
    set KPriv(dir) $dir
    set KPriv(problemTypeDir) $dir
    
    set ptypeName [lindex [split $KPriv(problemTypeDir) "/"] end]
    set KPriv(pTypeName) [string map {".gid" ""} $ptypeName]
    
    #List of node Id's
    set KPriv(groupsId) {}
    
    #List of material Id´s
    set KPriv(materialsId) {}
    set KPriv(materialsList) {}

    
    #Xml root
    set KPriv(xml) ""
    set KPriv(encrXml) ""
    set KPriv(xmlDoc) ""
    
    #Xml Materials root
    set KPriv(xmlMat) ""
    set KPriv(encrXmlMat) ""
    set KPriv(xmlDocMat) ""
    
    #Xml Functions
    set KPriv(xmlFun) ""
    set KPriv(encrXmlFun) ""
    set KPriv(xmlDocFun) ""
    
    #kratos_key_words.xml 
    set KPriv(xmlKKW) ""
    set KPriv(xmlDocKKW) ""
    
    #kratos.ini
    set KPriv(xmlIni) ""
    set KPriv(xmlDocIni) ""
    
    #No es necesario porque solo lo necesitamos para leer
    #set KPriv(encrXmlKKW) [lindex $xmlArray 1]
    
    
    # Set images directory
    set imagespath "$dir/images"
    set KPriv(imagesdir) $imagespath
    
    # Load all sources scripts
    ::kipt::LoadSourceFiles $dir

    # Change system menu
    # Preprocess
    ::kmtb::ChangePreprocessMenu $dir
    
    # Postprocess
    
    set GidPriv(ProgName) $ProgramName
    ChangeWindowTitle
    
    # Maintain the problem type
    GiD_Set MaintainProblemTypeInNew 1
}

proc kipt::FreePType {} {
    
    global KPriv

    # Destroy all pre/post open window
    
    # For group editor
    set w ".gid.kegroups" 
    if {[winfo exists $w]} {
    destroy $w
    }
    
    set w ".gid.kmprops" 
    if {[winfo exists $w]} {
    destroy $w
    }

    # Validation window
    set w ".gid.modelvalidation" 
    if {[winfo exists $w]} {
        ::KMValid::CreateReportWindowbClose $w 
    }
    
    # Close Project Settings Window if it exists
    set w ".gid.settingWin" 
    if {[winfo exists $w]} {
        ::kps::WindowbClose $w 
    }



    # ********************************
    #     End the bitmaps
    # ********************************
    # Preprocess 
    ::kmtb::EndCreatePreprocessTBar
    
    # Postprocess 
    
    #Limpiar los objetos tDom si aun existían
    catch { [$KPriv(xml) delete] }
    
    catch { [$KPriv(xmlMat) delete] }
    
    # Unset the problem type global variables
    UnsetGlobalVars

    # Reset namespaces

    global GidPriv
    set GidPriv(ProgName) "GiD"
    ChangeWindowTitle ""
}

proc kipt::LoadSourceMessage {tclfname} {
    global VersionNumber ProgramName
    
    WarnWin [= "Error reading file %s %s %s problem type has been not correctly installed" $tclfname $ProgramName $VersionNumber].
}

proc kipt::LoadSourceFiles {dir} {
    
    # Load the application scripts for Kratos applications
    global KPriv
    
    #if { [lsearch -exact $::auto_path [file join $dir scripts]] == -1 } {
    #        lappend ::auto_path [file join $dir scripts]
    #}
    # WarnWinText "::auto_path:$::auto_path"
    
    # Load some packages
    
    # WarnWinText "dir:$dir"
    # For scripts directory
    set scriptspath "$dir/scripts"
    
    if { [catch {source $scriptspath/files.tcl}] } {
    ::kipt::LoadSourceMessage files.tcl        
    return ""
    }
    if { [catch {source $scriptspath/winutils.tcl}] } {
    ::kipt::LoadSourceMessage winutils.tcl.tcl        
    return ""
    }
    if { [catch {source $scriptspath/menus.tcl}] } {
    ::kipt::LoadSourceMessage menus.tcl        
    return ""
    } 
    if { [catch {source $scriptspath/utils.tcl}] } {
    ::kipt::LoadSourceMessage utils.tcl        
    return ""
    }
    if { [catch {source $scriptspath/stringutils.tcl}] } {
    ::kipt::LoadSourceMessage stringutils.tcl        
    return ""
    } 
    if { [catch {source $scriptspath/modelvalidation.tcl}] } {
    ::kipt::LoadSourceMessage modelvalidation.tcl        
    return ""
    }
    if { [catch {source $scriptspath/projectSettings.tcl}] } {
    ::kipt::LoadSourceMessage projectSettings.tcl        
    return ""
    }
    
    # For xml libs
    set xmlpath "$dir/scripts/libs/xml"
    if { [catch {source $xmlpath/xmlutils.tcl}] } {
    ::kipt::LoadSourceMessage xmlutils.tcl
    return ""
    }
    if { [catch {source $xmlpath/xpathq.tcl}] } {
    ::kipt::LoadSourceMessage xpathq.tcl
    return ""
    }
    
    # For write calculation file
    set wkcfpath "$dir/scripts/libs/wkcf"        
    if { [catch {source $wkcfpath/wkcf.tcl}] } {
    ::kipt::LoadSourceMessage wkcf.tcl
    return ""
    }
    if { [catch {source $wkcfpath/wkcfutils.tcl}] } {
    ::kipt::LoadSourceMessage wkcfutils.tcl
    return ""
    }

    # Load kegroups
    set kegrouppath "$dir/scripts/kegroups"
    if { [catch {source $kegrouppath/kegroups.tcl}] } {
    ::kipt::LoadSourceMessage kegroups.tcl
    return ""
    }
    if { [catch {source $kegrouppath/kGroupEntities.tcl}] } {
    ::kipt::LoadSourceMessage kGroupEntities.tcl
    return ""
    }
    
    package require KEGroups
    
    # Load KMProps
    set kPropsPath "$dir/scripts/kmprops"
    
    if { [catch {source $kPropsPath/kmprops.tcl}] } {
        ::kipt::LoadSourceMessage kmprops.tcl
        return ""
    }
     if { [catch {source $kPropsPath/kmaterials.tcl} er] } {
        ::kipt::LoadSourceMessage kmaterials.tcl
        return ""
    }
    if { [catch {source $kPropsPath/kFunctions.tcl}] } {
        ::kipt::LoadSourceMessage kFunctions.tcl
        return ""
    }
}

proc kipt::InitPostProcess { } {
    
}

proc kipt::FreePostProcess { } {
    
    # Destroy all postprocess open window
    return "" 
    
}

proc kipt::LoadResultsToPostProcess {filename} {
    return ""
}

proc kipt::UpdateLanguage {language} {
    
    # WarnWinText "language:$language"
    set dir [GiD_Info problemtypepath]
    # WarnWinText "ProblemType: [GiD_Info problemtypepath]" 
    
    # Preprocess
    ::kmtb::ChangePreprocessMenu $dir
    
    # Postprocess
    
}

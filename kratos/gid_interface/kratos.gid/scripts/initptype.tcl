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
#     1.2-18/06/13-G. Socorro, delete the proc kipt::NewGiDGroups
#     1.1-22/10/12-J. Garate, Support for new GiD Groups
#     1.0-08/10/12-J. Garate, Enable/disable kipt::NewGiDGroups
#     0.9-01/10/12-J. Garate, Enable/disable Curves Module
#     0.8-20/09/12-J. Garate, add Curves, Tables and Plotgraph source files
#     0.7-04/05/12-G. Socorro, add a new variable to control the group deletion (when exists from the problem type)  
#     0.6-03/05/12-G. Socorro, Delete all group identifier using ::KUtils::DeleteAllGroupIdentifier and 
#                              close all group window Cond_Groups window close
#     0.5-10/04/12-G. Socorro, load new script in the (wkcffluid.tcl, etc.)
#     0.4-02/04/12-J.Garate, icon path change to adapt to GiD Themes
#     0.3-29/03/12-G. Socorro, load the new kmprops scripts
#     0.2-22/06/11-G. Socorro, delete KPriv(release) and create KPriv(RDConfig) in the Kratos.tcl
#     0.1-01/02/10-G. Socorro, create the base source code
#
###############################################################################

# ikpt => kratos Init Problem Type 

namespace eval kipt {

}

proc kipt::InitPType { dir } {
    
    global KPriv ProgramName
    
    # kipt::CheckLicense
    
    # Set dir to a global variable
    set KPriv(dir) $dir
    set KPriv(problemTypeDir) $dir
    
    # Variable to control the group deletion (when exists from the problem type)  
    set KPriv(Groups,DeleteGroup) 1

    set ptypeName [lindex [split $KPriv(problemTypeDir) "/"] end]
    set KPriv(pTypeName) [string map {".gid" ""} $ptypeName]
    
    # Set images directory
    ###########Aqui hay que separar entre GiD Classic y GiD Dark #########
    set imagespath "images/Classic"
    if {[gid_themes::GetCurrentTheme] == "GiD_black"} {
        set imagespath "images/Dark"
    }
    set KPriv(imagesdir) $imagespath
    
    # Change system menu
    # Preprocess
    # scripts/menus.tcl
    ::kmtb::ChangePreprocessMenu $dir
    
    # Postprocess
    
    set GidPriv(ProgName) $ProgramName
    ChangeWindowTitle
    
    # Maintain the problem type
    GiD_Set MaintainProblemTypeInNew 1
}

proc kipt::CheckLicense { } {
    package require verifp   
    # get list of all sysinfos: local and usb's
    WarnWin "all devices sysinfos: [vp_getsysinfo]"

    #try for a valid password
    set res [vp_getauthorization myprogname 1.1 * my_secret_key]
    set status [lindex $res 0]    
    if { $status != "VERSION_PRO" } {
        set msg [lindex $res 1]
        WarnWin "unregistered version. msg:$msg"
    } else {
        WarnWin "professional version."
    }
    #release password (specially if password is floating)
    vp_releaseauthorization myprogname 1.1 *
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
    

    # Limpiar los objetos tDom si aun existían
    catch { [$KPriv(xml) delete] }
    
    catch { [$KPriv(xmlMat) delete] }

    if { [info exists KPriv(xmlDoc) ] } { 
        $KPriv(xmlDoc) delete
    }
    
    
    if { [info exists KPriv(xml) ] } { 
        unset KPriv(xml)
    }
    
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
    if { [kipt::CurvesModule ] } {
        # For Curves, graphics and tables
        set curvepath "$dir/scripts/libs/curves"
        if { [catch {source $curvepath/curves.tcl}] } {
        ::kipt::LoadSourceMessage curves.tcl
        return ""
        }
        if { [catch {source $curvepath/tables.tcl}] } {
        ::kipt::LoadSourceMessage tables.tcl
        return ""
        }
        set curvepath "$dir/scripts/libs/graphics"
        if { [catch {source $curvepath/plotgraph.tcl}] } {
        ::kipt::LoadSourceMessage plotgraph.tcl
        return ""
        }
    }
    
    # For write calculation file
    set wkcfpath "$dir/scripts/libs/wkcf"
    if { [catch {source $wkcfpath/wkcf.tcl} cerror] } {
    # WarnWin $cerror
	::kipt::LoadSourceMessage wkcf.tcl
	return ""
    }
    if { [catch {source $wkcfpath/wkcfutils.tcl}] } {
	::kipt::LoadSourceMessage wkcfutils.tcl
	return ""
    }
    if { [catch {source $wkcfpath/wkcffluid.tcl}] } {
	::kipt::LoadSourceMessage wkcffluid.tcl
	return ""
    }
    if { [catch {source $wkcfpath/wkcfstructuralanalysis.tcl}] } {
	::kipt::LoadSourceMessage wkcfstructuralanalysis.tcl
	return ""
    }
    if { [catch {source $wkcfpath/wkcfgroups.tcl}] } {
	::kipt::LoadSourceMessage wkcfgroups.tcl
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
    if { [catch {source $kPropsPath/kmpropswin.tcl}] } {
	::kipt::LoadSourceMessage kmpropswin.tcl
	return ""
    }
    if { [catch {source $kPropsPath/kmpropsfwg.tcl}] } {
	::kipt::LoadSourceMessage kmpropsfwg.tcl
	return ""
    }
    if { [catch {source $kPropsPath/kmpropstree.tcl}] } {
	::kipt::LoadSourceMessage kmpropstree.tcl
	return ""
    }
    if { [catch {source $kPropsPath/kmpropsgroups.tcl}] } {
	::kipt::LoadSourceMessage kmpropsgroups.tcl
	return ""
    }
    if { [catch {source $kPropsPath/kmpropscbwd.tcl}] } {
	::kipt::LoadSourceMessage kmpropscbwd.tcl
	return ""
    }
    if { [catch {source $kPropsPath/kmaterials.tcl} er] } {
	msg $er
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

proc kipt::CurvesModule { } {
    global KPriv
    if { [info exists KPriv(CurvesModule)] } {
        return $KPriv(CurvesModule)
    } 
    return 0
}


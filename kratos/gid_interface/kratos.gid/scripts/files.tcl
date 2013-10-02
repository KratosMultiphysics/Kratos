###############################################################################
#
#    NAME: kfiles.tcl
#
#    PURPOSE: TCL script to work with the Kratos problem type files
#
#    QUANTECH ATZ-DEVELOPMENT DEPARTMENT
#
#    AUTHOR : G. Socorro
#
#    CREATED AT: 01/11/09
#
#    HISTORY:
#
#     1.6- 15/07/13- G. Socorro, modify the file giveConfigFile to add a dot before the name in linux platform
#     1.5- 17/06/13- G. Socorro, delete the procs RecursiveChildTransfer and TransferOldGroupstoGID to use only the new GiD groups
#     1.4- 12/02/13- J. Garate, ::kfiles::TransferOldGroupstoGID FIXED
#     1.3- 12/11/12- J. Garate, ::kfiles::TransferOldGroupstoGID modifications
#     1.2- 07/11/12- J. Garate, Adaptation for New GiD Groups on Transfer Function
#     1.1- 17/10/12- J. Garate, Correction when transferring old Cond_Groups to GiD_Groups
#     1.0- 08/10/12- J. Garate, Adaptation for New GiD Groups
#     0.9- 03/10/12- J. Garate, update transform spd proc
#     0.8- 27/05/12- J. Garate, Preparacion para actualizar la base de datos de materiales.
#     0.7- 03/05/12- J. Garate, state/visibility while transfering groups from .spd 
#     0.6- 26/04/12- G. Socorro, change GiD_Groups by Cond_Groups
#     0.5- 22/03/12- J. Garate,  Cambio de funciones a funciones publicas de GiD
#     0.4- 08/03/12- J. Garate,  ::kfiles::TransferOldGroupstoGID Mantiene la jerarquía Padre-Hijos en los grupos. 
#                                 El color identificativo de cada grupo se pasa de forma correcta.
#     0.3- 05/03/12- J. Garate,  ::kfiles::TransferOldGroupstoGID Transfiere los grupos del .spd a los grupos de GiD
#     0.2- 18/06/10- G. Socorro, Set ::KMat::xml path variable
#     0.1- 01/11/09- G. Socorro, create a base source code
#
###############################################################################

namespace eval ::kfiles {
    
}

proc ::kfiles::LoadSPD {filename} {
    
    global KPriv
    
    # wa "LoadSPD => filename:$filename"
    
    set KPriv(problemTypeDir) [file dirname $filename]
    
    # PROPERTIES
    set xmlNameFile "kratos_default.spd"
    
    if {![file exists $filename] || [file size $filename] < 1} {
        
        set filename "$KPriv(dir)/$xmlNameFile"
        
    } else {        
        #Se guarda una copia del archivo original antes de modificarlo
        
        if { [file tail $filename] == "kratos_default.spd" } {
            # kike: commented, this backup copy is unneeded 
            # (the 'kratos_default.spd' file must be opened for read only and then can't be corrupted)
            # and in general will be forbidden to write in the problemtype folder !!
            # it is only allowed at runtime to write in user folders (preferences, tmp, or model folder)
            # maybe this backup is interesting when reading the model .spd file to avoid possible corruption during transform??
            
            #::kfiles::MakeBackupCopyOfSPDFile $filename
        } else {
            ::kfiles::MakeBackupCopyOfSPDFile $filename
        }       
    }
    
    #::KEGroups::Init
    set xmlArray [::xmlutils::openFile "." "$filename"]
    set KPriv(xml) [lindex $xmlArray 0]
    # wa "[$KPriv(xml) asXML]"
    set KPriv(encrXml) [lindex $xmlArray 1]
    set KPriv(xmlDoc) [lindex $xmlArray 2]
    
    # Crea el archivo Kratos.ini, si aun no existe e inicializa
    #  las variables globales KPriv(xmlIni)...
    ::kps::updateKratosIniFile
    
    #Si estamos cargando el default se tendrán que comprobar las versiones
    if {$filename != "$KPriv(dir)/$xmlNameFile"} {
        
        #Transforma el spd si son versiones distintas
        if { [::xmlutils::checkSpdVersion $filename] } {
            ::xmlutils::UpdateSpd $filename
        }
        
    }
    # wa "after [$KPriv(xml) asXML]"
    #
    # MATERIALS
    #
    set filename_mat "[string range $filename 0 [expr [string length $filename] - 5]].kmdb"
    set xmlFile_mat "kratos_default.kmdb"
    #msg "filename_mat:$filename_mat xmlFile_mat:$xmlFile_mat"
    if {![file exists $filename_mat] || [file size $filename_mat] < 1} {
        
        set filename_mat "$KPriv(dir)/$xmlFile_mat"
        
    } else {
        #Se guarda una copia del archivo original antes de modificarlo
        
        if { [file tail $filename_mat] == "kratos_default.kmdb" } {
            # kike: commented, this backup copy is unneeded 
            # (the 'kratos_default.kmdb' file must be opened for read only and then can't be corrupted)
            # and in general will be forbidden to write in the problemtype folder !!
            # it is only allowed at runtime to write in user folders (preferences, tmp, or model folder)
            # maybe this backup is interesting when reading the model .kmdb file to avoid possible corruption during transform??
            
            #::kfiles::MakeBackupCopyOfSPDFile $filename_mat .kmdb
        } else {
            ::kfiles::MakeBackupCopyOfSPDFile $filename_mat .kmdb
        }
    }
    
    set xmlArray [::xmlutils::openFile "." "$filename_mat"]
    
    set KPriv(xmlMat) [lindex $xmlArray 0]
    set KPriv(encrXmlMat) [lindex $xmlArray 1]
    set KPriv(xmlDocMat) [lindex $xmlArray 2]
    
    ::xmlutils::checkMatVersion $filename_mat
    
    #
    # KKWORDS (KRATOS KEY WORDS)
    #
    set filePath "$KPriv(dir)/kratos_key_words.xml"
    set xmlArray [::xmlutils::openFile "." "$filePath" 0]
    
    set KPriv(xmlKKW) [lindex $xmlArray 0]
    set KPriv(xmlDocKKW) [lindex $xmlArray 2]
    #No es necesario porque solo lo necesitamos para leer
    #set KPriv(encrXml) [lindex $xmlArray 1]
    
    # IDIOMA: Lee todas las palabras del spd por si se quieren    
    # kike: commented, now ramtranslator scans also the kratos_default.spd and kratos_default.kmdb
    # and then is not necessary the trick of invoke ::xmlutils::getLanguageWords to parse them creating an auxiliary msgs\words.tcl file
    # and must not invoke this parsing every time a model is load, it is a kratos problemtype developer task only!!
    #::xmlutils::getLanguageWords
    
    
    
    # Set KMat xml path
    if {[info exists ::KMat::xml]} {
        set ::KMat::xml $KPriv(xmlMat)
    }
    
    # Nos guardamos el nombre del problemType cargado
    set ptypeName [lindex [split $KPriv(problemTypeDir) "/"] end]
    set KPriv(pTypeName) [string map {".gid" ""} $ptypeName]    
}

proc ::kfiles::SaveSPD {filename} {
    global KPriv;
    
    #Actualizamos los posibles cambios que haya habido en el ".spd"
    if {[info exists ::KMProps::WinPath]} {
        if {[winfo exists $::KMProps::WinPath]} {
            ::KMProps::RefreshTree "" 1
        }
    }
    
    # Coger de la ventana de grupos y guardar en $KPriv(xmlDoc)
    ::KEGroups::GroupsToXml
    # msg $filename
    # Escribimos en el fichero .spd el xml almacenado en memoria
    ::xmlutils::writeFile "${filename}" $KPriv(dir) $KPriv(encrXml) $KPriv(xmlDoc) $KPriv(RDConfig) 0
    
    # Escribimos en el fichero .kmdb de materiales el xml almacenado en memoria
    set materialFile "[string range $filename 0 [expr [string length $filename] - 5]].kmdb"
    ::xmlutils::writeFile "${materialFile}" $KPriv(dir) $KPriv(encrXmlMat) $KPriv(xmlDocMat) $KPriv(RDConfig) 0
    
    if {$KPriv(xmlDocFun) != ""} {
        set encryptFile 0
        ::xmlutils::writeFile "$KPriv(dir)/python_functions.xml" $KPriv(dir) $KPriv(encrXmlFun) $KPriv(xmlDocFun) $KPriv(RDConfig) $encryptFile
    }
    
    if {$KPriv(RDConfig) == 0} {
        # En modo debug tenemos que encriptar los defaults para pasarlos con la versión
        ::xmlutils::writeFile "${filename}_encrypt" $KPriv(dir) $KPriv(encrXml) $KPriv(xmlDoc) 1
        ::xmlutils::writeFile "${materialFile}_encrypt" $KPriv(dir) $KPriv(encrXmlMat) $KPriv(xmlDocMat) 1
    }
}

proc ::kfiles::MakeBackupCopyOfSPDFile {filename {extension ".spd"}} {
    
    # ABSTRACT:
    # Make a backup copy of file
    # ARGUMENTS:
    # Name -> Path and name of the file
    
    set basename [file tail $filename]
    # Rename the file
    # Get position of spd extension
    set word "$extension"
    #msg "filename$filename \n word:$word"
    set found [string first $word $basename]
    if {$found !="-1"} {
        set name [string range $basename 0 [expr $found-1]]
    } else {
        return ""        
    }
    set SPDBackup ${name}_backup$extension
    
    # Make backup copy at the problemtype dir
    set dir [file dirname $filename]
    #WarnWinText "SPDBackup:$dirGid/$SPDBackup"
    
    if {[catch {file copy -force "$filename" "$dir/$SPDBackup"} error]} {
        WarnWin [= "Could not make a backup copy of Kratos interface data file (%s)" $error]
        return ""
    }
}

#------------------------------------------------------------------------------
#
# Funciones para manejar el PROJECT SETTINGS
#
#------------------------------------------------------------------------------

proc ::kfiles::giveConfigFile { {ptypeName "kratos"} } {
    
    global KPriv
    
    set filename [GiveGidDefaultsFile]
    set dirname [file dirname $filename]
    set extname [file extension $filename]
    # wa "filename:$filename extname:$extname dirname:$dirname"
    if {$extname ne ".ini"} {
        set extname .ini
    }
    if {($::tcl_platform(os) eq "Linux")} {
        set kname ".${ptypeName}${extname}"
    } else {
        set kname "${ptypeName}${extname}"
    }
    
    set endpath [file native [file join $dirname $kname]]
    return $endpath
}

proc ::kfiles::varOnConfigFile { var {def 0} } {
    
    # Check the value of variable in the configuration file
    set file [::kfiles::giveConfigFile]
    if { [catch { set fileid [open $file r] }] } {
        return $def
    }
    while { ![eof $fileid] } {
        set aa [gets $fileid]
        if { [regexp "$var\[ ]*(\[01])" $aa {} val] } {
            close $fileid
            return $val
        }
    }
    close $fileid
    return $def
}

proc ::kfiles::writeVarConfigFile { varArray valArray } {
    
    # Write the value of variable to the configuration file
    
    set fileid ""
    
    set file [::kfiles::giveConfigFile]
    if { [file readable $file] } {
        if { [catch { set fileid [open $file r] }] } {
            return 0
        }
        set fileread [read $fileid]
        close $fileid
    } else {
        set fileread ""
    }
    
    foreach var $varArray val $valArray {
        
        if { [regexp "$var" $fileread {}] } {
            regsub "$var\[ ]*(\[01])" $fileread "$var $val"  filewrite
        } else {
            set filewrite ""
            append filewrite $fileread "\n$var $val"
        }
        if { [catch { set fileid [open $file w+] }] } {
            WarnWin [= "Cannot write file %s. Permission denied" $file].
            return 0
        }
        puts $fileid $filewrite
    }
    if { $fileid != "" } {
        close $fileid
    }
    return 1
}

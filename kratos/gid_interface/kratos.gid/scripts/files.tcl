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
#    LAST MODIFICATION : Set ::KMat::xml path variable
#
#    VERSION : 0.2
#
#    HISTORY:
# 
#     0.2- 18/06/10-G. Socorro, Set ::KMat::xml path variable
#     0.1- 01/11/09-G. Socorro, create a base source code
#
###############################################################################

namespace eval ::kfiles:: {

}
 
proc ::kfiles::LoadSPD {filename} {
    
    global KPriv
    
    set KPriv(problemTypeDir) [file dirname $filename]
    
    #
   	# PROPERTIES
   	#
    set xmlNameFile "kratos_default.spd"
    
    if {![file exists $filename] || [file size $filename] < 1} {
		
		set filename "$KPriv(dir)/$xmlNameFile"
			
   	} else {
   		
		#Se guarda una copia del archivo original antes de modificarlo
    	::kfiles::MakeBackupCopyOfSPDFile $filename
   	}
   	
   	#::KEGroups::Init
   	set xmlArray [::xmlutils::openFile "." "$filename"]
   	
   	set KPriv(xml) [lindex $xmlArray 0]
   	set KPriv(encrXml) [lindex $xmlArray 1]
   	set KPriv(xmlDoc) [lindex $xmlArray 2]
   	
   	# Crea el archivo Kratos.ini, si aun no existe e inicializa
	#  las variables globales KPriv(xmlIni)...
	::kps::updateKratosIniFile
    
   	#Si estamos cargando el default se tendrán que comprobar las versiones
   	if {$filename != "$KPriv(dir)/$xmlNameFile"} {
   		
   		#Transforma el spd si son versiones distintas
 		::xmlutils::checkSpdVersion $filename
 	}
 	
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
    	::kfiles::MakeBackupCopyOfSPDFile $filename_mat ".kmdb"
	}
	
   	set xmlArray [::xmlutils::openFile "." "$filename_mat"]
   	
   	set KPriv(xmlMat) [lindex $xmlArray 0]
   	set KPriv(encrXmlMat) [lindex $xmlArray 1]
   	set KPriv(xmlDocMat) [lindex $xmlArray 2]
    
    #
    # KKWORDS (KRATOS KEY WORDS)
    #
    set filePath "$KPriv(dir)/kratos_key_words.xml"
	set xmlArray [::xmlutils::openFile "." "$filePath" 0]
	
	set KPriv(xmlKKW) [lindex $xmlArray 0]
	set KPriv(xmlDocKKW) [lindex $xmlArray 2]
	#No es necesario porque solo lo necesitamos para leer
	#set KPriv(encrXml) [lindex $xmlArray 1]
    
    
	#
	# IDIOMA: Lee todas las palabras del spd por si se quieren 
	#
	::xmlutils::getLanguageWords
	
	
	

    # Set KMat xml path
    if {[info exists ::KMat::xml]} {
		set ::KMat::xml $KPriv(xmlMat)
    }
    
    # Nos guardamos el nombre del problemType cargado
    set ptypeName [lindex [split $KPriv(problemTypeDir) "/"] end]
    set KPriv(pTypeName) [string map {".gid" ""} $ptypeName]
    
    
    #PRUEBAS: Inicios rápidos para probar
    #::kps::InitSettingWindow
	# ::KMValid::ValidateModel
	#::KMProps::InitBaseWindow Materials
	#::KMProps::InitBaseWindow
    #::KFun::InitBaseWindow ".gid.kmprops.nb.fProp.fBottom.nb.fMainProperties" "Fx"
    
    
}

proc ::kfiles::SaveSPD {filename} {
    global KPriv
    
    #Actualizamos los posibles cambios que haya habido en el ".spd"
    if {[info exists ::KMProps::WinPath]} {
    	if {[winfo exists $::KMProps::WinPath]} {
	    ::KMProps::refreshTree "" 1
    	}
    }
	
    # Escribimos en el fichero .spd el xml almacenado en memoria
    ::xmlutils::writeFile "${filename}" $KPriv(dir) $KPriv(encrXml) $KPriv(xmlDoc) $KPriv(RDConfig)
    
    # Escribimos en el fichero .spd de materiales el xml almacenado en memoria
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
    set kname "${ptypeName}${extname}"
    
    return [file join $dirname $kname]
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
		if { [catch { set fileid [open $file r] }] } { return 0 }
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

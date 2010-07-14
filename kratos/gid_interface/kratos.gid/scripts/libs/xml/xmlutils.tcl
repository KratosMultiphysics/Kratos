###############################################################################
#
#    NAME: xmlutils.tcl
#
#    PURPOSE: Utilities procedures to work with XML files
#
#    QUANTECH ATZ-DEVELOPMENT DEPARTMENT
#
#    AUTHOR : G. Socorro
#
#    CREATED AT: 01/11/09
#
#    LAST MODIFICATION : add some new utilities procedures from the wiki http://wiki.tcl.tk/4193
#
#    VERSION : 0.2
#
#    HISTORY:
#
#     0.2- 24/12/09-G. Socorro, add some new utilities procedures from the wiki http://wiki.tcl.tk/4193
#     0.1- 01/11/09-G. Socorro, create a base source code
#
###############################################################################

package require tdom

# Create a base namespace xmlutils
namespace eval ::xmlutils:: {

}

proc ::xmlutils::initKKWord { } {
	
	#global KPriv
    #
    #if { $KPriv(xmlKKW) == "" } {
	#
	#        set filePath "$KPriv(dir)/kratos_key_words.xml"
	#        set xmlArray [::xmlutils::openFile "." "$filePath" 0]
	#        
	#        set KPriv(xmlKKW) [lindex $xmlArray 0]
	#        set KPriv(xmlDocKKW) [lindex $xmlArray 2]
	#        #No es necesario porque solo lo necesitamos para leer
	#        #set KPriv(encrXml) [lindex $xmlArray 1]
	#
    #}
}

proc ::xmlutils::getKKWord {xpath id {cattr kkword}} {
    
    global KPriv
    
    set node [$KPriv(xmlDocKKW) selectNodes "Kratos_KWords/$xpath/Item\[@id='$id'\]"]
    
    if { $node != "" } {
	return [$node getAttribute $cattr ""]
    }
}

#******************************************************************************
#
#################################################
proc ::xmlutils::AsXml {content {addtag No} {tag document}} {
    # From wiki: http://wiki.tcl.tk/1740

    set XML_MAP [list < "&lt;" > "&gt;" & "&amp;" \" "&quot;" ' "&apos;"]
    if {$addtag=="Yes"} {
	return <$tag>[string map $XML_MAP $content]</$tag>
    } else {
	return [string map $XML_MAP $content]
    }
}

#-------------------------------------------------------------------------------------------------
# Abre el fichero "inputfile" del directorio "iniDir" y extrae el código xml 
#-------------------------------------------------------------------------------------------------
proc ::xmlutils::openFile {w fullname {encodeFile 1}} {
    
    #set fullname [file native [file join $iniDir $inputfile]]
    #set dirpath [file dirname $fullname]
    
    # Read the input file 
    if {![file exists $fullname]} {
	set msg "The file ($fullname) does not exists. Check that this file exists"
	WarnWin "$msg." 
	return 0
    }
    
    
    set xml [tDOM::xmlReadFile $fullname encriptStr]

    #set xmlDoc [dom parse $xml]
    #set fullname "faslkdfjlaskdgjkratos_default.spd"
    #msg $fullname
    #msg [string first "_default" $fullname]
    
    set xmlDoc ""
    catch {
    	set xmlDoc [dom parse [::xmlutils::DecodeK [string trim $xml]]]
    }
    
    catch {
	    if { $xmlDoc == "" } {
	    	set xmlDoc  [dom parse [string trim $xml]]
	    }
	}
   	if { $xmlDoc == "" } {
   		WarnWin "Format error parsing the xml document\n'$fullname'"
   	}
    
    #if { $encodeFile && [string first "_default" $fullname] == -1 } {
	#    set xmlDoc  [dom parse [::xmlutils::DecodeK [string trim $xml]]]
    #} else {
	#    set xmlDoc  [dom parse [string trim $xml]]
    #}
    
    #Controlar un posible error al hacer el dom parse
    #msg "error?:[catch { dom parse $xml }]"
    #if { [catch { [set xmlDoc [dom parse $xml]] }] == 0 } {
    #        WarnWin "Format error parsing the xml document\n'$fullname'"
    #        return -1
    #}
    
    set xmlList [$xmlDoc documentElement]
    set xmlList [lappend xmlList $encriptStr]
    set xmlList [lappend xmlList $xmlDoc]
    
    #msg "xmlList:$xmlList"
    return $xmlList

}

#-------------------------------------------------------------------------------------------------
# Crea o Abre el fichero "outputfile" en el directorio "iniDir" y 
#   le inserta el código xml q tenemos en memoria 
#-------------------------------------------------------------------------------------------------
proc ::xmlutils::writeFile {outputfile iniDir encStr xmlDoc {release 1} {encryptFile 1}} {
    
    #Si este archivo no hay que encriptarlo no comprobamos nada aquí
    if { $encryptFile } {
    	
	    #Si estamos en modo DEBUG, guardamos el archivo de las dos formas
	    if { !$release } {
	    	
	    	::xmlutils::writeFile "$outputfile" $iniDir $encStr $xmlDoc 1
	    	set outputfile "${outputfile}_debug"
	    }
	}
    
    set fullname [file native [file join $iniDir $outputfile]]
    set dirpath [file dirname $fullname]
    
    #msg "encStr: $encStr xmlDoc:$xmlDoc"
    #Extrae la cabecera
    set outfd [open $fullname w+]
    fconfigure $outfd -encoding [::tDOM::IANAEncoding2TclEncoding $encStr]
    
    #La escribe en el fichero de salida
    puts $outfd "<?xml version='1.0' encoding='$encStr'?>"
    
    #Comprobamos si es necesario encriptarlo
    if { $release && $encryptFile} {
    	
	    puts $outfd [::xmlutils::EncodeK [$xmlDoc asXML]]
    } else {
	    puts $outfd [$xmlDoc asXML]
    }
    
    close $outfd
    #$xmlDoc delete
    
    # Para tener el código visible por ejemplo para pruebas
    #set outfd [open ${fullname}pruebas.spd w+]
    #fconfigure $outfd -encoding [::tDOM::IANAEncoding2TclEncoding $encStr]
    
    #La escribe en el fichero de salida
    #puts $outfd "<?xml version='1.0' encoding='$encStr'?>"
    #puts $outfd [$xmlDoc asXML]
    #close $outfd
}

#Si el problemtype ya estaba guardado, guardamos el xml y lo volvemos a cargar
proc ::xmlutils::reloadFile { xmlDocument encriptXml } {
    
    set dirGid [GiD_Info problemtypepath]
    #msg "dirGid:$dirGid"
}

proc ::xmlutils::EncodeK { x } {
    return [encrypter BASE64 encode [encoding convertto utf-8 $x]]
}

proc ::xmlutils::DecodeK { x } {
    return [encoding convertfrom utf-8 [encrypter BASE64 decode $x]]
}

#-------------------------------------------------------------------------------------------------
#
# Crea o Abre el fichero "/msgs/words.tcl" y le añade todas las palabras para traducir
#
#-------------------------------------------------------------------------------------------------
proc ::xmlutils::getLanguageWords { } {
    
    global KPriv
    set xml $KPriv(xml)
    
    set xmlMat $KPriv(xmlMat)
    
    set nodes {}
    lappend nodes [$xml getElementsByTagName "RootData"]
    lappend nodes [$xml getElementsByTagName "Container"]
    lappend nodes [$xml getElementsByTagName "Item"]
    lappend nodes [$xml getElementsByTagName "Template"]
    lappend nodes [$xmlMat getElementsByTagName "Material"]
    lappend nodes [$xmlMat getElementsByTagName "Property"]
    lappend nodes [$xmlMat getElementsByTagName "MaterialGroup"]
    
    
    dict create word val
    
    #Recorremos todos los nodos del xml
    foreach listNodes $nodes {
	foreach node $listNodes {
	    
	    #Palabras clave de usuario PID
	    set pid [$node getAttribute pid ""]
	    if { $pid != "" } {
		
		dict set word $pid val ""
	    }
	    
	    #Palabras clave de usuario PID
	    set help [$node getAttribute help ""]
	    if { $help != "" } {
		
		dict set word $help val ""
	    }
	    
	    #Contenido de los combos
	    set values [$node getAttribute values ""]
	    set Lvalues [split $values ","]
	    foreach cmbVal $Lvalues {
		dict set word $cmbVal val ""
	    }
	}}
    
    #Crea un archivo donde poner las palabras del diccionario
    set fullname "$KPriv(dir)/msgs/words.tcl"
    if { [catch { set file [open $fullname w] }] } {
	return 0
    }
    
    set lista {}
    dict for {id nada} $word {
	#puts $file "[= $id ]"
	#msg "id:$id"
	lappend lista $id
    }
    
    set lista [lsort $lista]
    
    foreach id $lista {
	puts $file "\[\= $id \]"
    }
    
    close $file
}


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Compara la versión del problemTypeName.spd abierto con la versión de kratos_default.spd
# Y si son distintas añade al spd los nodos y atributos nuevos del default
#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
proc ::xmlutils::checkSpdVersion { filename } {
    
    global KPriv
    
    #Abrimos el spd default
    set xmlFileDefault "$KPriv(dir)/kratos_default.spd"
    set xmlArray [::xmlutils::openFile "." "$xmlFileDefault"]
    set xmlDef [lindex $xmlArray 0]
    
    set xml $KPriv(xml)
    
    #Comprobamos las versiones
    set pTypeVersion [::xmlutils::xmlVersion]
    set defaultVersion [::xmlutils::xmlVersion $xmlDef]
    
    #msg "\n\n$pTypeVersion != $defaultVersion"
    
    if {$pTypeVersion != $defaultVersion } {
	
	variable logChanges {}
	
	set path [GiD_Info problemtypepath]
	set name [lindex [split $path "/"] end]
	set aviso [= "This model version it is older than the problem type one. The file: \n '$filename'\n it is going to be updated.\n A back-up of the original file will be generated."]
	set answer [WarnWin "$aviso" "."]
	
	#----------------------------------------------------------------
	# BLOQUE DE BORRADO
	#----------------------------------------------------------------
	set baseNode [$xmlDef selectNodes "/Kratos_Data"]
	
	set nodeNames {RootData Container Item Template}
	
	foreach nodeName $nodeNames {
	    
	    set nodes [$xml getElementsByTagName "$nodeName"]
	    foreach node $nodes {
		
			#Bloqueamos los posibles errores intentando eliminar nodos ya eliminados
			catch {
			    set class [$node getAttribute class ""]
			    if { $nodeName != "Container" || ($class != "Group" && $class != "Property") } {
			        
			        set id [$node getAttribute id ""]
			        set foundNode [$baseNode find id $id]
			        if { $foundNode == "" } {
			            
			    		lappend ::xmlutils::logChanges "DELETE NODE [$node nodeName] Id: [$node getAttribute id 0]"
			    		
			            #Si el nodo ya no está lo borramos
			            $node delete
			        }
			    }
			}
	    }        
	}
	
	
	#----------------------------------------------------------------
	# INSERCIÓN NODOS / AÑADIDO DE ATRIBUTOS
	#----------------------------------------------------------------
	set baseNode [$xml selectNodes "/Kratos_Data"]
	
	#Tipos de nodos a actualizar
	set nodeNames {RootData Container Item Template}
	
	foreach nodeName $nodeNames {
	    
	    set nodes [$xmlDef getElementsByTagName "$nodeName"]
	    
	    foreach node $nodes {
			
			set id [$node getAttribute id ""]
			
			set foundNode [$baseNode find id $id]
	
			if { $foundNode != "" } {
			    
			    if { [$foundNode nodeName] == [$node nodeName] } {
			    	
			    	::xmlutils::updateAtributesRecursive $baseNode $foundNode $node $id
			    }
			    
			} else {
			    
			    #Copia el node de Kratos_default q no aparece en el problemtype
			    ::xmlutils::copyNode $node $baseNode
			}
	    }
	}
	
	#Crea un archivo de log con los cambios realizados
	if { [llength $logChanges] } {
		
		set fullname [file native [file join $KPriv(problemTypeDir) updateSPD.log]]
		set outfd [open $fullname w+]
	    
	    puts "Changes in configuration file .spd to the new version $defaultVersion"
	    puts ""
	    foreach line $logChanges {
	    	#La escribe en el fichero de salida
	    	puts $outfd $line
	    }
	    
	    close $outfd
	 }
    
	
	::xmlutils::xmlVersion $xml $defaultVersion
	
    }
}

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Compara los atributos y añade los q falten y elimina los que sobren
#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
proc ::xmlutils::updateAtributes { targetNode sourceNode } {
    
    #set Tpath [::xmlutils::getPathFromNode $targetNode]
    #set Spath [::xmlutils::getPathFromNode $sourceNode]
    #
    #if {$Tpath != $Spath } {
    #        msg "Be carefull, trying to update attributes with different paths:\n$Tpath\n$Spath\n"
    #        return
    #}
    
    #Si el nodo es un grupo o una propiedad, no estaba en default
    set class [$targetNode getAttribute class ""]
    if { $class != "Group" && $class != "Property" } {

	set atributes [$sourceNode attributes]
	foreach atr $atributes {
	    
	    if { ![$targetNode hasAttribute $atr] } {
		
			set value [$sourceNode @$atr]                
			$targetNode setAttribute $atr $value
			
			lappend ::xmlutils::logChanges "INSERT ATRIBUTE [$targetNode nodeName] Id: [$targetNode getAttribute id ""]:$atr=\"$value\""
	    }
	}
	
	set atributes [$targetNode attributes]
	foreach atr $atributes {
	    
	    if { ![$sourceNode hasAttribute $atr] } {
		
		#Si el nodo original ya no lo tiene, es que se ha borrado
		$targetNode removeAttribute $atr
		
		lappend ::xmlutils::logChanges "DELETE ATTRIBUTE [$targetNode nodeName] Id: [$targetNode getAttribute id ""]->$atr=\"$value\""
	    }
	}
    }
}
proc ::xmlutils::updateAtributesRecursive { baseNode targetNode sourceNode id } {
    
    set Tpath [::xmlutils::getPathFromNode $targetNode]
    set Spath [::xmlutils::getPathFromNode $sourceNode]
    
    if {$Tpath != $Spath } {
	
		set nodes [$baseNode descendant all]
		foreach node $nodes {
		    set fNode [$node find id $id]
		    if { $fNode != "" } {
			::xmlutils::updateAtributesRecursive $node $fNode $sourceNode $id
		    }
		}
		return
    }
    
    #Si el nodo es un grupo o una propiedad, no estaba en default
    set class [$targetNode getAttribute class ""]
    if { $class != "Group" && $class != "Property" } {
	
	set atributes [$sourceNode attributes]
	foreach atr $atributes {
	    
	    if { ![$targetNode hasAttribute $atr] } {
		
			set value [$sourceNode @$atr]                
			$targetNode setAttribute $atr $value
			
			lappend logChanges "INSERT ATTRIBUTE [$targetNode nodeName] Id: [$targetNode getAttribute id ""] -> $atr=\"$value\""
	    }
	}
	
	set atributes [$targetNode attributes]
	foreach atr $atributes {
	    
	    if { ![$sourceNode hasAttribute $atr] } {
		
		#Si el nodo original ya no lo tiene, es que se ha borrado
		$targetNode removeAttribute $atr
		
		lappend ::xmlutils::logChanges "DELETE ATTRIBUTE [$targetNode nodeName] Id: [$targetNode getAttribute id ""] -> atr:$atr"
	    }
	}
    }
}

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Dado un determinado nodo de TDom obtiene su path en las formas: 
# 1:  RootDataNode//..c.nodoN..//c.nodoPadre//c.finalNode
# 2:  RootDataNode//..c.nodoN..//c.nodoPadre
#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
proc ::xmlutils::getPathFromNode { finalNode {includeFinal 1}} {
    
    set path ""
    
    set nodes [$finalNode ancestor all]
    #msg "$finalNode $nodes"
    set i 0
    foreach node $nodes {
	
	set id [$node getAttribute id ""]
	if { $id != "" } {
	    set path "c.$id//$path"
	}
    }
    
    if {$path != ""} {
	set path [string range $path 2 end]
	if { $includeFinal } {
	    set path "${path}c.[$finalNode getAttribute id 0]"
	} else {
	    set path [string range $path 0 end-2]
	}
	
	return [::xmlutils::setXPath $path]
    }
    return ""
    
}

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# SourceNode es el nodo del default.spd que queremos copiar en el spd del problemType
# BaseNode es el nodo Kratos_data del problemType.sdp (donde vamos a copiar el nodo source)
#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
proc ::xmlutils::copyNode { sourceNode baseNode } {
    
    global KPriv
    
    set parentNode [$sourceNode parentNode]
    
    if {$parentNode == "" } { return }
    
    set id [$parentNode getAttribute id ""]
    
    set foundParentNode [$baseNode find id $id]
    
    #Caso especial ROOT_DATA
    if { [$parentNode nodeName] == "Kratos_Data" } {
	
	set foundParentNode $baseNode
	$foundParentNode appendChild $sourceNode
	
	lappend ::xmlutils::logChanges "INSERT NODE [$sourceNode nodeName] Id: [$sourceNode getAttribute id 0]"
	#msg "INSERT [$sourceNode nodeName] Id:[$sourceNode getAttribute id 0] xml:[$sourceNode asXML]"
	return ""
	
	#Caso especial TEMPLATES
    } elseif { [$parentNode nodeName] == "Templates" } {
	
	$parentNode appendChild $sourceNode
	lappend ::xmlutils::logChanges "INSERT NODE [$sourceNode nodeName] Id: [$sourceNode getAttribute id 0]"
	#msg "INSERT [$sourceNode nodeName] Id:[$sourceNode getAttribute id 0] xml:[$sourceNode asXML]"
	
	#El resto de nodos
    } else {
	
	#Si no es vacío ya sabemos donde meter "sourceNode"
	if { $foundParentNode != "" } {
	    
	    set nodoSiguiente [$sourceNode nextSibling]
	    
	    if {$nodoSiguiente != "" } {
		
		set idSiguiente [$nodoSiguiente getAttribute id ""]
		
		#Buscamos si también está este nodo donde lo vamos a insertar,
		# y si no está, lo pondremos el último
		set foundBro [$foundParentNode find id $idSiguiente]
		
		if { $foundBro != "" } {
		    $foundParentNode insertBefore $sourceNode $foundBro
		} else {
		    #Si no tiene nodo siguiente lo ponemos el último
		    $foundParentNode appendChild $sourceNode
		}
		lappend ::xmlutils::logChanges "INSERT NODE [$sourceNode nodeName] Id: [$sourceNode getAttribute id 0]"
		#msg "INSERT [$sourceNode nodeName] Id:[$sourceNode getAttribute id 0] xml:[$sourceNode asXML]"
		
	    } else {
		#Si no tiene nodo siguiente lo ponemos el último
		$foundParentNode appendChild $sourceNode
		lappend ::xmlutils::logChanges "INSERT NODE [$sourceNode nodeName] Id: [$sourceNode getAttribute id 0]"
		#msg "INSERT [$sourceNode nodeName] Id:[$sourceNode getAttribute id 0] xml:[$sourceNode asXML]"
	    }
	}
    }
    return ""
}

#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#           #   #     #     #     #
#            # #      # # # #     #
#             #       #  #  #     #
#            # #      #     #     #
#           #   #     #     #     ######
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#

proc ::xmlutils::xmlVersion { {xml ""} {change ""} } {
    
    if { $xml == "" } {
	global KPriv
	set xml $KPriv(xml)
    }
    if { $change == "" } {
	set version [$xml set "/Kratos_Data/@version"]
	return $version
    } else {
	$xml set "/Kratos_Data/@version" $change
    }
}

# Prepara la query para utilizar las funciones de domNOde
proc ::xmlutils::setXPath { path { type "props"} } {
    
    set splitted [::KMProps::split2 $path //]
    
    if { $type != "props"} {
	set i 0
	foreach itemId $splitted {
	    if { $i == 0 } {                
		#                                                        set xpath "/Kratos_KMat_DB/Materials/Material\[@id='$itemId'\]"                                                  
		set xpath "/Kratos_KMat_DB/Materials/MaterialGroup\[@id='$itemId'\]"                                                  
		
	    } else {
		if { [string index $itemId 0] == "p" } {
		    set xpath "$xpath/Property\[@id='[string range $itemId 2 end]'\]"
		}
		if { [string index $itemId 0] == "m" } {
		    set xpath "$xpath/Material\[@id='[string range $itemId 2 end]'\]"
		}
		if { [string index $itemId 0] == "c" } {
		    set xpath "$xpath/Container\[@id='[string range $itemId 2 end]'\]"
		}
	    }        
	    incr i
	}
	return $xpath
    } else {
	set i 0
	#msg "s:$splitted"
	foreach itemId $splitted {
	    
	    if { $i == 0 } {                
		#El primer elemento será siempre del nivel 'RootData'
		set xpath "/Kratos_Data/RootData\[@id='$itemId'\]"
	    } else {
		
		if { [string index $itemId 0] == "c" } {
		    set xpath "$xpath/Container\[@id='[string range $itemId 2 end]'\]"
		} else {
		    set xpath "$xpath/Item\[@id='[string range $itemId 2 end]'\]"
		}
	    }
	    incr i
	}
    }
    return $xpath
}

#
# Editar o extraer propiedades del xml en memoria
#
proc ::xmlutils::setXml { path property {command "read"} {value ""} } {
    
    global KPriv
    
    set xpath "[::xmlutils::setXPath $path]"
    
    if { $command == "read" } {
	
	set value [$KPriv(xml) set "$xpath/@$property" ]
	
	if { $property == "dvText" } {
	    
	    #En vez de devolver el valor de dv, devuelve su equivalente traducible
	    set value [$KPriv(xml) set "$xpath/@dv" ]
	    
	    if { [lindex $value 0] != "" } {
		
		set values [::xmlutils::getXMLValues $path]
		set ivalues [::xmlutils::getXMLValues $path "" "iValues"]
		
		if { [llength $ivalues] > 0 } { 
		    
		    set index [::xmlutils::getSelected $value $ivalues]
		    set value [lindex $values $index]
		}
	    }
	}
	
	#Cuando hay espacios el xml devuelve una lista y si la imprimes tal cual aparecen corchetes
	if { [llength $value] == 1 } {
	    set value [lindex $value 0]
	}
	return $value
	
    } else {
	
	$KPriv(xml) set "$xpath/@$property" "$value"
	return "1"                        
    }
}

#
# Accede a la ruta $path y devuelve una lista con todos los ID de primer nivel
#
proc ::xmlutils::getXmlNodeName { path {type "props"}} {
    
    global KPriv
    
    set xpath "[::xmlutils::setXPath $path $type]"
    
    if { $type == "props" } {
	set nodes [$KPriv(xml) selectNodes "${xpath}"]
    } else {
	set nodes [$KPriv(xmlMat) selectNodes "${xpath}"]
    }
    
    if {[llength $nodes]} {
	return [[lindex $nodes 0] nodeName]
    } else {
	return "noNodes"
    }
}

#
# Accede a la ruta $path y devuelve una lista con todos los ID de primer nivel
#
proc ::xmlutils::setXmlContainerIds { path {nodeType "Container"} { type "props" } } {
    
    global KPriv
    
    set listIds {}
    
    set xpath "[::xmlutils::setXPath $path $type]"
    
    set nodes [$KPriv(xml) selectNodes "${xpath}/$nodeType"]
    
    foreach node $nodes {
	
		lappend listIds [$node getAttribute id ""]
    }
    
    return $listIds
}

#
# Accede a la ruta $path y devuelve una lista de parejas { propiedad valor } de sus items
# idNode: Si nos llega un "id" solo devolvemos el valor de su "atributo"
#
proc ::xmlutils::setXmlContainerPairs { path {id ""} {atributo "dv"} {tagname "Item"} {type "props"}} {
    
    global KPriv
    
    set listPares {}
    
    set xpath "[::xmlutils::setXPath $path $type]"
    
    if { $type == "mat"} {
	set xml $KPriv(xmlMat)
    } else {
	set xml $KPriv(xml)
    }

    set node [$xml selectNodes "${xpath}"]
    
    #        set items [$node descendant all "Item"]
    set items [$node descendant all $tagname]
    
    foreach node $items {
	
	if { $id != "" && $id == [$node getAttribute id ""] } {
	    return  [$node getAttribute $atributo ""]
	}
	lappend listPares [list [$node getAttribute id ""] [$node getAttribute $atributo ""]]
    }
    
    return $listPares
}

proc ::xmlutils::GetMatXmlContainerId {path} {
    global KPriv
    
    set xpath "[::xmlutils::setXPath $path mat]"
    set xml $KPriv(xmlMat)
    set node [$xml selectNodes "${xpath}"]
    WarnWinText "xpath:$xpath node:$node"
    set cidlist [list]
    foreach node [$node descendant all Container] {
	lappend cidlist [$node getAttribute id ""]
    }
    WarnWinText "cidlist:$cidlist"
    return $cidlist
}

proc ::xmlutils::unsetXml { path { type "props" } } {
    
    global KPriv                
    if { $type == "props" } {
	set xpath "[::xmlutils::setXPath $path]"
	$KPriv(xml) unset $xpath
	
    } else {
	set xpath "[::xmlutils::setXPath $path $type]"
	$KPriv(xmlMat) unset $xpath        
    }
}

proc ::xmlutils::insertXml { path nodeName properties {xml ""} } {
    
    if { $xml == "" } {
	global KPriv
	$KPriv(xml)
    }
    
    #RootData\[@id='$itemId'\]
    if { $path == "root" } { 
	set xpath "/Kratos_Data/$nodeName"
    } else {
	set xpath "[::xmlutils::setXPath $path]/$nodeName"
    }
    
    foreach prop value $properties {
	set xpath "$xpath $prop=\"$value\""
    }
    
    $xml lappend "$xpath"
}

#
#*************************************************************
#*        #####    #####    #     #    ####    #####          
#*        #        #   #    # # # #    #   #   #   #          
#*        #        #   #    #  #  #    ####    #   #          
#*        #        #   #    #     #    #   #   #   #          
#*        #####    #####    #     #    ####    #####          
#*************************************************************
#
# Accede al XML al nodo con path 'fullname' y coge el 
# atributo "values" o la lista especial correspondiente
#
proc ::xmlutils::getXMLValues { fullname {idTemplate ""} {iValues ""} {idTemplateFull ""} {specialFilter ""}} {
    
    #msg "$application --> $comboList\nargs:1$fullname 2$idTemplate 3$iValues 4$idTemplateFull 5$specialFilter"
    
    global KPriv
    
    set comboList {}
    
    if { $iValues != "" } {
	    set atrValues "ivalues"
	} else {
	    set atrValues "values"
	}
		
    if { $idTemplate == "" } {
	#Se utiliza el fullname normalmente
	set specialList [::xmlutils::setXml $fullname GCV]
    } else {
	#Se tiene que consultar la propiedad en el template
	set specialList [::KMProps::getPropTemplate $idTemplate GCV "$fullname"]
    }
    #Si no hay lista especial, cargamos el contenido de "values="
    if { $specialList == "" } {
		
		if { $idTemplate == "" } {
		    set comboList [split [::xmlutils::setXml $fullname $atrValues] ","]
		} else {
		    set comboList [split [::KMProps::getPropTemplate $idTemplate $atrValues "$fullname"] ","]
		}
	
    } elseif { $specialList == "Materials" } {
		
		if {$idTemplateFull == "" } {
			set application [::KMProps::getApplication $fullname]
		} else {
			set application [::KMProps::getApplication $idTemplateFull]
		}
		
		set comboList [::KMat::getMaterials $application]		
		
    } elseif { $specialList == "Properties" } {
		
		if {$idTemplateFull != "" } {
	
		        set comboList [::KMProps::getProps $idTemplateFull]
		} else {
		        
		        set comboList [::KMProps::getProps $fullname]
		}
		
    } elseif { $specialList == "ElemType" } {
    
	    #::xmlutils::initKKWord
	    
	    set node [$KPriv(xmlDocKKW) selectNodes "Kratos_KWords/ElementCLaws/Item\[@id='ElementTypes'\]"]
	    
	    if { $node != "" } {
		    set comboList [split [$node getAttribute $atrValues ""] ","]
		    $node setAttribute dv "$specialFilter"
	    }
	    
    } elseif { $specialList == "MatModel" } {
    
    	set dvElemFilter ""
    	if { $specialFilter == "" } {
    		#Caso especial para los combos dinámicos:
	    	#Mira a ver si el elemento anterior era ElemType para utilizar su valor
	    	set elemtypeFullname [string map {"MatModel" "ElemType"} $fullname]
	    	set dvElemFilter [::xmlutils::setXml $elemtypeFullname dv]
	    }
	    #Caso estándard: accede al valor dv de element type para coger el matModel filtrando por este elementType
	    if { $dvElemFilter == "" } {
	    	#Si el valor no existe, miramos el que hay en el xml KKWORDS
	    	set node [$KPriv(xmlDocKKW) selectNodes "Kratos_KWords/ElementCLaws/Item\[@id='ElementTypes'\]"]
		    if { $node != "" } {
			    set dvElemFilter [$node getAttribute dv ""]
			}
		}
		#Cogemos la lista correspondiente
	    set node [$KPriv(xmlDocKKW) selectNodes "Kratos_KWords/ElementCLaws/Item\[@id='$dvElemFilter'\]"]
	    if { $node != "" } {
	            set comboList [split [$node getAttribute $atrValues ""] ","]
	    }
	}
    #Si sacamos los ids no los traducimos
    if { $iValues == "" } {
		set resultList {}
		foreach val $comboList {
		    set traduction "[= $val ]"
		    #set traduction "$val"
		    lappend resultList $traduction
		}
		return $resultList
	} else {
		
		return $comboList
    }
    
    
}

#
# Dado el path del combo, y el fullname para saber donde encontrar los values
# devolvemos $dv (id $value = "id") o su equivalente traducible
#
proc ::xmlutils::getComboDv { fcmb fullname {value "id"} {idTemplate ""}} {
    
    set index [$fcmb current]
    
    set comboList [::xmlutils::getXMLValues $fullname "$idTemplate"]
    set icomboList [::xmlutils::getXMLValues $fullname "$idTemplate" "iValues"]
    
    if { $index < [llength $icomboList] } {
	
	if { $value == "id" } {
	    return [lindex $icomboList $index]
	} else {
	    return [lindex $comboList $index]
	}
    } else {
    	#El valor
    	if {[winfo exists $fcmb] } {
    		return [$fcmb get]
    	}
		msg "Fuera de índice $index (length: [llength $comboList])"
    }
}

proc ::xmlutils::setComboDv { fcmb fullname dv {idTemplate ""} } {
    
    set icomboList [::xmlutils::getXMLValues $fullname "$idTemplate" "iValues"]
    set selected [::xmlutils::getSelected $dv $icomboList]
    $fcmb current $selected
}

proc ::xmlutils::getSelected { comboItem comboList} {
	
	set i 0
	foreach iElem $comboList {
		if {$iElem == $comboItem} {
		        return $i
		}
		set i [expr $i + 1]
	}
	return 0
}


proc ::xmlutils::getComboBoxState { fullname } {
    
    set values [::xmlutils::setXml $fullname values]
    set GCV [::xmlutils::setXml $fullname GCV]
    set CBState [::xmlutils::setXml $fullname CBState]
	
    if { ($GCV == "" && [llength $values] == 0 ) || $CBState == "normal"} {
		return normal
    } else {
		return readonly
    } 
}

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++

proc ::xmlutils::AsXml {content {action toXml}} {
	
    # From wiki: http://wiki.tcl.tk/1740
	if { $action == "toXml" } {
	     set XML_MAP [list < "&lt;" > "&gt;" & "&amp;" \" "&quot;" ' "&apos;"]
    } else {
	     set XML_MAP [list "&lt;" <  "&gt;" > "&amp;" & "&quot;" \" "&apos;" ']
    }
    
    return [string map $XML_MAP [string trim $content]]
}

#
# TODO: return the value when setting a value and return a list of values when multiple values are set
#
# From wiki: http://wiki.tcl.tk/4193
# The basic idea is that values from an xml document can be read and modified using XPath queries. 
# This really just amounts to a simpler interface to the domNode object than what tDOM provides.
# Creation command:
#    * xmlstruct::create xml - Returns an extended domNode object command
# The following methods are supported in addition to those the tDOM domNode object already provides:
#    * $node set xpathQuery ?value? - retrieves or modifies portions of the xml document that match the given xpathQuery
#    * $node unset xpathQuery - deletes the portions of the xml document that match the given xpathQuery
#    * $node lappend xpathQuery value ?value? ?value? ... - appends the given values to the node(s) that match the given xpathQuery

# By placing these procs in the ::dom::domNode namespace, they automatically become add-on domNode methods
proc ::dom::domNode::unset {node query} {
    ::set resultNodes [$node selectNodes $query type]
    switch $type {
	attrnodes {::xmlstruct::unsetattrs $node $query}
	nodes {::xmlstruct::unsetnodes $resultNodes}
	empty {error "No results found for '$query'"}
	default {error "$type is an unsupported query result type"}
    }
}

proc ::dom::domNode::set {node query args} {
    switch [llength $args] {
	0 {return [::xmlstruct::getvalue $node $query]}
	1 {return [::xmlstruct::setvalue $node $query [lindex $args 0]]}
	default {error "wrong # args: should be \"set xpathQuery ?newValue?\""}
    }
}
proc ::dom::domNode::nodeType {node} {

}

proc ::dom::domNode::lappend {node query args} {
    foreach arg $args {
	::xmlstruct::setnew $node $query $arg
    }
}

# Create the name space xmlstruct 
namespace eval ::xmlstruct:: {}

proc ::xmlstruct::create {xml} {
    # Convenience function for creating an xml doc and returning the root
    ::set doc [dom parse $xml]
    return [$doc documentElement]
}

proc ::xmlstruct::getvalue {node query} {
	
    # For '$node set query' calls
    ::set resultNodes [$node selectNodes $query type]
    switch $type {
	attrnodes {
	    ::set retVal {}
	    foreach attrVal $resultNodes {
		lappend retVal [lindex $attrVal 1]
	    }
	    return $retVal
	}
	nodes {
	    ::set retVal {}
	    foreach node $resultNodes {
		::set xml ""
		foreach child [$node childNodes] {
		    append xml [$child asXML]
		}
		lappend retVal $xml
	    }
	    # This is so the curly braces are not there due to the above lappend
	    if {[llength $resultNodes] == 1} {::set retVal [lindex $retVal 0]}
	    return $retVal
	}
	empty {return ""}
	default {error "$type is an unsupported query result type"}
    }
}


proc ::xmlstruct::setvalue {node query value} {
    # For '$node set query value' calls
    ::set targetNodes [$node selectNodes $query type]
    switch $type {
	nodes {::xmlstruct::setnodes $targetNodes $query $value}
	attrnodes {::xmlstruct::setattrs $node $query $value}
	empty {::xmlstruct::setnew $node $query $value}
	default {error "$type is an unsupported query result type"}
    }
}

proc ::xmlstruct::setnew {node query value} {
    # Creates a new attribute/element for an xpath query in which all
    # the elements of the query up to the last exist
    set possibleMatch [split $query /]
    set unmatched [lindex $possibleMatch end]
    set possibleMatch [lreplace $possibleMatch end end]
    if {[llength $possibleMatch] == 0} {
	set possibleMatch .
    }
    
    set nodes [$node selectNodes [join $possibleMatch /] type]
    switch $type {
	nodes {
	    if {[string index $unmatched 0] == "@"} {
		foreach node $nodes {
		    $node setAttribute [string range $unmatched 1 end] $value
		}
	    } else {
		foreach node $nodes {
		    $node appendXML "<$unmatched/>"
		    set newNode [$node lastChild]
		    $newNode set . $value
		}
	    }
	}
	attrnodes {error "Can't add children to attributes ($possibleMatch)"}
	empty {error "Create elements matching $possibleMatch first"}
    }
}

proc ::xmlstruct::unsetattrs {node query} {
    # For i.e. '$node unset {/employees/employee[1]/@age}' calls
    ::set nodeQuery [join [lrange [split $query /] 0 end-1] /]
    ::set attribute [string range [lindex [split $query /] end] 1 end]
    foreach matchingNode [$node selectNodes $nodeQuery] {
	$matchingNode removeAttribute $attribute
    }
}

proc ::xmlstruct::setattrs {node query value} {
    # For i.e. '$node set {/employees/employee[1]/@age} 25' calls
    ::set nodeQuery [join [lrange [split $query /] 0 end-1] /]
    ::set attribute [string range [lindex [split $query /] end] 1 end]
    foreach matchingNode [$node selectNodes $nodeQuery] {
	$matchingNode setAttribute $attribute $value
    }
    return $value
}

proc ::xmlstruct::unsetnodes {nodes} {
    # For i.e. '$node unset {/employees/employee[1]}' calls
    # This probably breaks if some nodes are descendents of each other and
    # they don't get deleted in the right order
    foreach node $nodes {
	$node delete
    }
}

proc ::xmlstruct::isXml {string} {
    # Determines if the given string is intended to be valid xml
    ::set string [string trim $string]
    if {([string index $string 0] == "<") && [string index $string end] == ">"} {
	return 1
    } else {
	return 0
    }
}

proc ::xmlstruct::setnodes {targetNodes query value} {
    # For i.e. '$node set {/employees/employee[1]} value' calls
    if {[::xmlstruct::isXml $value]} {
	foreach target $targetNodes {::xmlstruct::setxml $target $value}
    } else {
	foreach target $targetNodes {::xmlstruct::settext $target $value}
    }
}

proc ::xmlstruct::settext {node text} {
    # TODO: don't allow this to be called for the documentElement node
    # (i.e. $obj set / "some text"  should not be allowed)
    # For i.e. '$node set {/employees/employee/name} Bill' calls
    ::set doc [$node ownerDocument]
    foreach child [$node childNodes] {$child delete}
    if {[string length $text] > 0} {
	::set textNode [$doc createTextNode $text]
	$node appendChild $textNode
    }
    return $text
}

proc ::xmlstruct::setxml {node xml} {
    # For i.e. '$node set {/employees/employee} <name>Bill</name>' calls
    foreach child [$node childNodes] {$child delete}
    $node appendXML $xml
    return $xml
}


#
#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#           #   #     #     #     #
#            # #      # # # #     #
#             #       #  #  #     #
#            # #      #     #     #
#           #   #     #     #     ######
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
#

#
# Si no hay atributo devuelve uno o más nodos, y si lo hay, uno o más atributos.
# Además, si value no es nulo, hace un set en el xml
#
proc ::xmlutils::getAttribute { xml xPath {attribute ""} {value ""} } {
	
	set nodes [$xml selectNodes "$xPath"]
    
    if { $nodes != "" } {
		 
		 if { $attribute == "" } {
		         return $nodes
		 } else {
		         if { [llength $nodes] == 1 } {
		                 
		                 if { $value != "" } {
		                         return [$nodes setAttribute $attribute $value]
		                 } else {        
		                         return [$nodes getAttribute $attribute ""]
		                 }
		         } else {
		                 set listAttr {}
		                 foreach node $nodes {
		                         if { $value != "" } {
		                                 $nodes setAttribute $attribute $value
		                         } else {        
		                                 lappend listAttr [$node getAttribute $attribute ""]
		                         }
		                         
		                 }
		                 return $listAttr
		         }
		 }        
    }
}

#
# Accede al XML al nodo con xpath y coge el 
# atributo "values" o la lista especial correspondiente
#
proc ::xmlutils::getValues { xml xpath {iValues ""} } {
      
    set comboList {}
    
    if { $iValues != "" } {
	    set atrValues "ivalues"
	} else {
	    set atrValues "values"
	}
		
	set specialList [::xmlutils::getAttribute $xml $xpath GCV]
    
    #Si no hay lista especial, cargamos el contenido de "values="
    if { $specialList == "" } {
		
		set comboList [split [::xmlutils::getAttribute $xml $xpath $atrValues] ","]
    
    } elseif { $specialList == "OtrasAcciones" } {
		#Aquí se extraerían los values de otros lugares ...
	}
	
    #Si sacamos los ids no los traducimos
    if { $iValues == "" } {
		set resultList {}
		foreach val $comboList {
		    set traduction "[= $val ]"
		    #set traduction "$val"
		    lappend resultList $traduction
		}
		return $resultList
	} else {
		
		return $comboList
    }
    
    
}

#
# Consulta el valor en dv, lo compara con la lista de valores 
# y devuelve el string traducible correspondiente
#
proc ::xmlutils::getValueText { xml xpath {attribute "dv"} } {
	    
	    #En vez de devolver el valor de dv, devuelve su equivalente traducible
	    set value [::xmlutils::getAttribute $xml $xpath $attribute]
	    
	    if { [lindex $value 0] != "" } {
		
			set values [::xmlutils::getValues $xml $xpath]
    		set ivalues [::xmlutils::getValues $xml $xpath "iValues"]
			
			if { [llength $ivalues] > 0 } { 
			    
			    set index [::xmlutils::getSelected $value $ivalues]
			    set value [lindex $values $index]
			}
	    }
	    #Retorna el valor traducible equivalente o simplemente "dv" si no hay ivalues
	    return $value
}

#
# Dado el path del combo, y el xpath para saber donde encontrar los values
# devolvemos $dv (id $value = "id") o su equivalente traducible
#
proc ::xmlutils::getComboValue { xml xpath fcmb {value "id"} } {
    
    if { ![winfo exists $fcmb] } {
    	return "error"
    }
    
    set index [$fcmb current]
    
    set comboList [::xmlutils::getValues $xml $xpath]
    set icomboList [::xmlutils::getValues $xml $xpath "iValues"]
    
    if { $index < [llength $icomboList] } {
	
	if { $value == "id" } {
	    return [lindex $icomboList $index]
	} else {
	    return [lindex $comboList $index]
	}
    } else {
		msg "Fuera de índice $index (length: [llength $comboList])"
    }
}

proc ::xmlutils::setComboValue { xml xpath fcmb valor } {
    
    set icomboList [::xmlutils::getValues $xml $xpath "iValues"]
    set selected [::xmlutils::getSelected $valor $icomboList]
    $fcmb current $selected
}

proc ::xmlutils::getComboState { xml xpath } {
    
    set values [::xmlutils::getAttribute $xml $xpath values]
    set GCV [::xmlutils::getAttribute $xml $xpath GCV]
    
    if { $GCV == "" && [llength $values] == 0 } {
		return normal
    } else {
		return readonly
    }
}

proc ::xmlutils::copyTemplate { xml xpath templatePath idTemplate nodeName attributesArray} {
    
    #set idTemplate "Forces"
#    set id [::xmlutils::getPropTemplate "$idTemplate" id]
#    set pid [::xmlutils::getPropTemplate "$idTemplate" pid]
#    set help [::xmlutils::getPropTemplate "$idTemplate" help]
#    set icon [::xmlutils::getPropTemplate "$idTemplate" icon]
    
#    set xpath "[::xmlutils::setXPath $fullname]"
    
#    set node [$xml selectNodes "${templatePath}\[@id='noPath'\]"]
    set template [$xml set "${templatePath}\[@id='$idTemplate'\]"]
    
    #No se puede insertar en el xml un fragmento con mas de un nodo, por eso utilizamos 
    #las etiquetas auxiliares "<gouptemplate>$template</gouptemplate>"
    set template "<groupTemplate>$template</groupTemplate>"

	set textAttr ""
    foreach attrValue $attributesArray {
	set textAttr "$textAttr $attrValue"
    }
    
    $xml lappend "$xpath/$nodeName $textAttr" $template

    #Eliminamos ahora las etiquetas auxiliares "<gouptemplate></gouptemplate>"
#    ::xmlutils::replaceTemplate
    
    return $template
}

proc ::xmlutils::replaceTemplate { xml {xpath ""} } {
    
    global KPriv
    
    set xmlText [$xml asXML]
    
    set xmlText [string map {"<groupTemplate>" "" "</groupTemplate>" ""} $xmlText]
    
    set newXmlDoc [dom parse $xmlText]
    
    set newXml [$newXmlDoc documentElement]

    return [list $newXmlDoc $newXml]
}

#
# Accede a la ruta $path y devuelve una lista con todos los ID de primer nivel
#
proc ::xmlutils::getXmlChildIds { xml xpath } {
    
    set listIds {}
    
    set nodes [$xml selectNodes "$xpath"]
    
    foreach node $nodes {
		
		set id [$node getAttribute id ""]
		if { $id != "" } {
			lappend listIds $id
		}
    }
    
    return $listIds
}

package provide xmlstruct 1.0
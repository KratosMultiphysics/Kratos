###############################################################################
#
#        NAME: xmlutils.tcl
#
#        PURPOSE: Utilities procedures to work with XML files
#
#        QUANTECH ATZ-DEVELOPMENT DEPARTMENT
#
#        AUTHOR : G. Socorro
#
#        CREATED AT: 01/11/09
#
#        HISTORY:
#
#       2.1- 16/07/13- G. Socorro, prepare the proc openFile to use the problemtype translation
#       2.0- 15/07/13- G. Socorro, correct a bug in the proc copyTemplate, add some comment to the proc ::xmlstruct::setnew
#       1.9- 04/07/13- A.Melendo, modify ::xmlstruct::setnew to accept groups with any name
#       1.8- 06/05/13- G. Socorro, add the proc GetPropertySectionType, indent the source code
#       1.7- 12/02/13- J. Garate,   Modification on UpdateSpd. modeltype for Groups.
#       1.6- 13/12/12- J. Garate,   Modification on UpdateSpd. If an Item or Container is hidden at Default.spd, keep it hidden.
#       1.5- 07/11/12- J. Garate,   Modification on ::xmlutils::setXml to accept xpath or path as parameter
#       1.4- 09/10/12- G. Socorro,  add the proc GetPropertyElemType 
#       1.3- 03/10/12- J. Garate,   Update ::xmlutils::UpdateSpd
#       1.2- 01/10/12- J. Garate,   Enable/disable Curves Module
#       1.1- 20/09/12- J. Garate,   Adaptacion de más funciones para las creacion / edicion de curvas 
#       1.0- 20/07/12- J. Garate,   Adaptacion de más funciones para los Tabs de Materiales 
#       0.9- 19/07/12- J. Garate,   Adaptacion de funciones para los Tabs de Materiales (MyPathFromNode, parentNodePath, setXML, getXMLvalues)
#       0.8- 04/06/12- J. Garate,   Template select when transferring user materials
#       0.7- 27/05/12- J. Garate,   ::xmlutils::checkMatVersion y funciones auxiliares. Actualiza la base de datos de materiales
#       0.6- 04/05/12- J. Garate,   cambio en el log
#       0.5- 26/04/12- J. Garate,   ::xmlutils::checkSpdVersion  ::xmlutils::myPathFromNode   ::xmlutils::GetOldDvFromNewNode
#       0.4- 27/02/12- J. Garate,   Correccion de la condicion de entrada a la validacion del .spd. Edicion de la misma funcion.
#       0.3- 03/09/10- G. Socorro,  correct an error with outfd
#       0.2- 24/12/09- G. Socorro,  add some new utilities procedures from the wiki http://wiki.tcl.tk/4193
#       0.1- 01/11/09- G. Socorro,  create a base source code
#
###############################################################################

package require tdom

# Create a base namespace xmlutils
namespace eval ::xmlutils:: {
    
    variable logchanges {}
}

proc ::xmlutils::initKKWord { } {
    
}

proc ::xmlutils::GetPropertyElemType {baseappprop propId} {
    # ABSTRACT: Get the property element type for a specific property identifier
    set cxpath "$baseappprop//c.[list ${propId}]//c.MainProperties//i.ElemType"
    set cproperty "dv"
    set PropertyElemType ""
    set checkitem [::xmlutils::setXml $cxpath $cproperty check]
    if {$checkitem !=""} {
        set PropertyElemType [::xmlutils::setXml $cxpath $cproperty]
    } 
    return $PropertyElemType
}

proc ::xmlutils::GetPropertySectionType {baseappprop propId} {
    # ABSTRACT: Get the property section type for a specific property identifier
    
    set cxpath "$baseappprop//c.[list ${propId}]//c.MainProperties//i.SectionType"
    set cproperty "dv"
    set PropertySectionType ""
    set checkitem [::xmlutils::setXml $cxpath $cproperty check]
    # wa "checkitem:$checkitem"
    if {$checkitem !=""} {
        set PropertySectionType [::xmlutils::setXml $cxpath $cproperty]
    }
    return $PropertySectionType
}

proc ::xmlutils::GetSpatialDimension {} {
    # ABSTRACT: Get the spatial dimension
    
    set cxpath "GeneralApplicationData//c.Domain//i.SpatialDimension"
    set cproperty "dv"
    set ndime [::xmlutils::setXml $cxpath $cproperty]
    
    return $ndime
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

#
#-------------------------------------------------------------------------------------------------
# Abre el fichero "inputfile" del directorio "iniDir" y extrae el código xml 
#-------------------------------------------------------------------------------------------------
#
proc ::xmlutils::openFile {w fullname {encodeFile 1}} {
    
    #set fullname [file native [file join $iniDir $inputfile]]
    #set dirpath [file dirname $fullname]
    
    # Read the input file 
    if {![file exists $fullname]} {
        set msg [= "The file (%s) does not exist. Check that this file exists" $fullname]
        WarnWin "${msg}." 
        return 0
    }
    
    set xml [tDOM::xmlReadFile $fullname encriptStr]
    
    #Inicializamos la variable que contendrá el documento parseado
    set xmlDoc ""
    catch {
        #Intentamos leer y parsear el xml encriptado
        set xmlDoc [dom parse [::xmlutils::DecodeK [string trim $xml]]]
    }
    
    catch {
        #Si está vacío es que falló porque no estaba encriptado, lo intentamos leer normal
        if { $xmlDoc == "" } {
            set xmlDoc  [dom parse [string trim $xml]]
        }
    }
    if { $xmlDoc == "" } {
        #Si sigue fallando el xml tiene que ser erróneo
        WarnWin [= "Format error parsing the xml document '%s'" $fullname]
    }
    
    #Devolvemos las 3 variables la primera "xml" para utilizar el documento, 
    # y las dos siguientes "encr" y "xmlDoc" para cuando lo guardemos        
    
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
            
            ::xmlutils::writeFile "${outputfile}_debug" $iniDir $encStr $xmlDoc 1
            #set outputfile "${outputfile}"
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
        }
    }
    
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

proc ::xmlutils::getMatpTypeVersion { xml } {
    set raw [$xml set "/Kratos_KMat_DB/@version"]
    if {[string range $raw 0 0] == "\{"} {
        set raw [string range [$xml set "/Kratos_KMat_DB/@version"] 1 end-1]
    }
    return $raw
}

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#                             ---- UPDATE KMDB -----
# Compara la versión del problemTypeName.kmdb abierto con la versión de kratos_default.kmdb
# Y si son distintas crea un nuevo problemTypeName.kmdb copia del default, y luego le añade
# las diferencias como grupos, properties y valores de combos (atributo "value")
#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
proc ::xmlutils::checkMatVersion { filename } {
    
    global KPriv
    
    #Abrimos el kmdb default
    set xmlFileDefault "$KPriv(dir)/kratos_default.kmdb"
    set xmlArray [::xmlutils::openFile "." "$xmlFileDefault"]
    set xmlDef [lindex $xmlArray 0]
    set encrXmlDef [lindex $xmlArray 1]
    
    #Este es el xml del modelo actual
    set xmlOld $KPriv(xmlMat)
    
    #Comprobamos las versiones
    set pTypeVersion [::xmlutils::getMatpTypeVersion $xmlOld]
    set defaultVersion [::xmlutils::getMatpTypeVersion $xmlDef]
    
    if {$pTypeVersion != $defaultVersion} { 
        
        #msg [$xmlOld asXML]
        variable logChanges {}
        
        set xmlDocNew [dom parse [$xmlDef asXML]]
        set xmlNew [$xmlDocNew documentElement]
        
        set path [GiD_Info problemtypepath]
        set name [lindex [split $path "/"] end]
        msg "You are working on the $pTypeVersion Kratos Material DataBase version and the current version is $defaultVersion ."
        set warning [= "This model version it is older than the problem type one. The file: \n '$filename'\n it is going to be updated.\n A back-up of the original file will be generated."]
        msg $warning
        
        #--------------------------------------------------------------------
        # RECORRER TODOS LOS NODOS COMPROBANDO EL VALOR VALUE (tb open y state)
        #--------------------------------------------------------------------
        
        set baseNodePaths {}
        lappend baseNodePaths "/Kratos_KMat_DB/Materials/MaterialGroup\[@id='Metal'\]"
        lappend baseNodePaths "/Kratos_KMat_DB/Materials/MaterialGroup\[@id='Fluid'\]"
        lappend baseNodePaths "/Kratos_KMat_DB/Materials/MaterialGroup\[@id='Plastic'\]"
        lappend baseNodePaths "/Kratos_KMat_DB/Materials/MaterialGroup\[@id='Composite'\]"
        lappend baseNodePaths "/Kratos_KMat_DB/Materials/MaterialGroup\[@id='DEMMaterial'\]"
        
        # Añadir todos los tipos de material que haya en MaterialGroup
        foreach baseNodePath $baseNodePaths {
            # wa "baseNodePath:$baseNodePath"
            # Para cada item, seleccionamos los nodos en el xml antiguo
            set baseNodeOld [$xmlOld selectNodes $baseNodePath]
            set materialsOld [list]
            # Avoid problem with not defined material group in the old material database
            if {$baseNodeOld !=""} {
                set xmlOldPartial [[dom parse [$baseNodeOld asXML]] documentElement]
                set materialsOld [$xmlOldPartial getElementsByTagName "Material"]
            }
            
            # Recorremos el nuevo xml también por partes
            set baseNodeNew [$xmlNew selectNodes $baseNodePath]
            set xmlNewPartial [[dom parse [$baseNodeNew asXML]] documentElement]
            # Creamos listas de materiales para verficar si el usuario tiene materiales suyos
            set materialsNew [$xmlNewPartial getElementsByTagName "Material"]
            
            set matnewlist {}
            set matoldlist {}
            
            foreach matNew $materialsNew {
                set matname [$matNew getAttribute id ""]
                lappend matnewlist $matname $matNew
            }
            foreach matOld $materialsOld {
                set matname [$matOld getAttribute id ""]
                lappend matoldlist $matname $matOld
            }
            
            # wa "matnewlist:$matnewlist\nmatoldlist:$matoldlist"
            # Organizamos los materiales en listas. 
            # changelist es la lista de materiales nuevos, añadidos por el usuario. Hay que copiar el template y copiar los valores
            # checklist es la lista de materiales que ya estaban, hay que comprobar si hay cambios en el default
            set changelist {}
            set checklist {}
            foreach { matOldname matOldnode} $matoldlist {
                if { $matOldname ni $matnewlist } {
                    lappend changelist $matOldname $matOldnode
                    #msg [$matOldnode asXML]
                } else { 
                    lappend checklist $matOldname $matOldnode
                }
            }
            # wa "changelist:$changelist"
            if {[llength $changelist]} {
                ::xmlutils::CopyUserMaterialtoxmlNew $xmlOld $xmlNew $changelist $baseNodeNew
            }
            
            # En este punto, tenemos en xmlNew una copia del Default, con los materiales que el usuario haya añadido.
            # Falta mirar si el usuario habia tocado los materiales standard, para quedarnos con sus valores.
            # wa "checklist:$checklist"
            if {[llength $checklist]} {
                ::xmlutils::CheckMateriallist $xmlOld $xmlNew $checklist
            }
            
            # En este punto, tenemos en xmlNew una copia del Default, con los materiales que el usuario haya añadido
            # y con los valores definidos por el usuario.
            # Falta definir el xmlNew como xml a usar a partir de ahora.
            
            set KPriv(xmlMat) $xmlNew
            #msg [$KPriv(xmlMat) asXML]
            set KPriv(xmlDocMat) $xmlDocNew
            set KPriv(encrXmlMat) $encrXmlDef
            
            # Finalizada la transferencia de materiales, no olvide guardar para no perder los cambios
        }
    }
}

proc ::xmlutils::findMaterialParent { xml matid } {
    
    
    set nodes [$xml selectNodes "/Kratos_KMat_DB/Materials/MaterialGroup\[@id\]"]
    
    foreach node $nodes {                          
        set nodes2 [$node childNodes]
        foreach node2 $nodes2 {  
            set aux [$node2 getAttribute id ""]
            if { $aux == $matid} {
                set parent [$node getAttribute id ""]
                return $parent
            }
        }                
    }
}

proc ::xmlutils::CopyUserMaterialtoxmlNew {xmlOld xmlNew changelist baseNodeNew} {
    # Copia los materiales que hay en la changelist, de xmlOld a xmlNew
    # copiando el template y pasando las caracteristicas una a una.
    
    foreach {UserMaterialName UserMaterialNode} $changelist  {
        
        set mattype [::xmlutils::findMaterialParent $xmlOld $UserMaterialName]
        set matTemplateid "NewMaterial"
        if {$mattype eq "DEMMaterial"} {set matTemplateid "NewDEMMaterial"}
        
        set nodeTempl [[$xmlNew find id $matTemplateid] asList]
        set nodeTempl [lreplace $nodeTempl [lsearch $nodeTempl Template] [lsearch $nodeTempl Template] Material]
        set aux [lindex $nodeTempl 1]
        set newaux [lreplace $aux 1 1 AuxMaterial]
        set nodeTempl [lreplace $nodeTempl [lsearch $nodeTempl $aux] [lsearch $nodeTempl $aux] $newaux]
        $baseNodeNew appendFromList $nodeTempl
        
        set nodeTemplate [$xmlNew find id "AuxMaterial"]
        set NewNode [$nodeTemplate cloneNode -deep]
        
        set vOld [$UserMaterialNode getAttribute id ""]
        $NewNode setAttribute id $vOld
        
        set vOld [$UserMaterialNode getAttribute pid ""]
        $NewNode setAttribute pid $vOld
        
        ::xmlutils::RecursiveValueTrans $NewNode $UserMaterialNode
        
        $baseNodeNew appendChild $NewNode 
        set auxnode [$xmlNew find id "AuxMaterial"]
        $auxnode delete
    }
}

proc ::xmlutils::RecursiveValueTrans { nodeNew nodeOld } {
    
    if { [$nodeNew hasAttribute value] } {
        set vOld [$nodeOld getAttribute value ""]
        $nodeNew setAttribute value $vOld
    } 
    if { [$nodeNew hasAttribute unit] } {
        set vOld [$nodeOld getAttribute unit ""]
        $nodeNew setAttribute unit $vOld
    }
    if { [$nodeNew hasChildNodes] } {
        foreach childNew [$nodeNew childNodes] childOld [$nodeOld childNodes] {
            if {$childOld != "" } {
                if {$childNew != "" } {
                    ::xmlutils::RecursiveValueTrans $childNew  $childOld
                }
            }
        }
    }
}

proc ::xmlutils::CheckMateriallist { xmlOld xmlNew checklist } {
    # Comprueba las propiedades de los materiales que hay en la checklist.
    # si hay diferencias entre el xmlOld y el xmlNew, copia de Old a New
    
    foreach {UserMaterialName UserMaterialNode} $checklist  {
        set NewNode [$xmlNew find id $UserMaterialName]
        #set NewNode [$nodeTemplate cloneNode -deep]
        
        set vOld [$UserMaterialNode getAttribute id ""]
        $NewNode setAttribute id $vOld
        
        set vOld [$UserMaterialNode getAttribute pid ""]
        $NewNode setAttribute pid $vOld
        
        ::xmlutils::RecursiveValueTrans $NewNode $UserMaterialNode
        
        #$baseNodeNew appendChild $NewNode 
        #msg [$baseNodeNew asXML]
    }
}

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#                                                  ---- UPDATE SPD -----
# Compara la versión del problemTypeName.spd abierto con la versión de kratos_default.spd
# Y si son distintas crea un nuevo problemTypeName.spd copia del default, y luego le añade
# las diferencias como grupos, properties y valores de combos (atributo "dv")
#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
proc ::xmlutils::checkSpdVersion { filename } {
    
    global KPriv
    
    #Abrimos el spd default
    set xmlFileDefault "$KPriv(dir)/kratos_default.spd"
    set xmlArray [::xmlutils::openFile "." "$xmlFileDefault"]
    set xmlDef [lindex $xmlArray 0]
    set encrXmlDef [lindex $xmlArray 1]
    
    #Este es el xml del modelo actual
    set xmlOld $KPriv(xml)
    
    #Comprobamos en los settings del proyecto si es necesario validar la versión
    set validateSPD [::kps::getConfigValue "CheckSpdVersion"]
    
    # msg "validateSPD: $validateSPD"
    if { !$validateSPD } {
        
        return 0
    }
    
    #Comprobamos las versiones
    set pTypeVersion [::xmlutils::xmlVersion]
    set defaultVersion [::xmlutils::xmlVersion $xmlDef]
    
    # msg "pTypeVersion $pTypeVersion"
    #msg "defaultVersion $defaultVersion"
    
    if {$pTypeVersion != $defaultVersion} { 
        return 1
    }
    return 0
}

proc ::xmlutils::UpdateSpd {filename {outputDisplay 1} {outputLog 1}} {
    global KPriv
    
    # Abrimos el spd default
    set xmlFileDefault "$KPriv(dir)/kratos_default.spd"
    set xmlArray [::xmlutils::openFile "." "$xmlFileDefault"]
    set xmlDef [lindex $xmlArray 0]
    set encrXmlDef [lindex $xmlArray 1]
    
    # Este es el xml del modelo actual
    set xmlOld $KPriv(xml)
    
    # wa "Old xml:[$xmlOld asXML]"
    variable logChanges {}
    
    set xmlDocNew [dom parse [$xmlDef asXML]]
    set xmlNew [$xmlDocNew documentElement]
    
    set pTypeVersion [::xmlutils::xmlVersion]
    set defaultVersion [::xmlutils::xmlVersion $xmlDef]
    
    set path [GiD_Info problemtypepath]
    set name [lindex [split $path "/"] end]
    if {$outputDisplay} {
        msg [= "You are working on the %s Kratos version and the current version is %s." $pTypeVersion $defaultVersion]
        set warning [= "This model version it is older than the problem type one. The file: \n '%s'\n it is going to be updated.\n A back-up of the original file will be generated." $filename]
        msg $warning
    }
    
    #--------------------------------------------------------------------
    # RECORRER TODOS LOS NODOS COMPROBANDO EL VALOR DV (tb open y state)
    #--------------------------------------------------------------------
    
    # All rootdata
    set baseNodePaths {}
    lappend baseNodePaths "/Kratos_Data/RootData\[@id='GeneralApplicationData'\]"
    lappend baseNodePaths "/Kratos_Data/RootData\[@id='StructuralAnalysis'\]"
    #lappend baseNodePaths "/Kratos_Data/RootData\[@id='Fluid'\]"
    #lappend baseNodePaths "/Kratos_Data/RootData\[@id='PFEM'\]"
    #lappend baseNodePaths "/Kratos_Data/RootData\[@id='FluidStructureInteraction'\]" 
    #lappend baseNodePaths "/Kratos_Data/RootData\[@id='ConvectionDiffusion'\]" 
    lappend baseNodePaths "/Kratos_Data/RootData\[@id='DEM'\]"
    
    if { [kipt::CurvesModule ] } {
        lappend baseNodePaths "/Kratos_Data/RootData\[@id='Curves'\]"
    }
    #msg $baseNodePaths
    
    foreach baseNodePath $baseNodePaths {
        set numLogs 0
        # lappend ::xmlutils::logChanges "Sección $baseNodePath:"
        
        # Para cada item, buscamos su correspondiente por partes 
        set baseNodeOld [$xmlOld selectNodes $baseNodePath]
        # Si no existe el nodo antiguo, salta al siguiente
        if { $baseNodeOld == ""} { continue }
        
        # Recorremos el nuevo xml también por partes
        set baseNodeNew [$xmlNew selectNodes $baseNodePath]
        set xmlNewPartial [[dom parse [$baseNodeNew asXML]] documentElement]
        
        # Para actualizar los rootData
        $baseNodeNew setAttribute state [$baseNodeOld getAttribute state ""]
        $baseNodeNew setAttribute open [$baseNodeOld getAttribute open ""]
        
        set nodes [$xmlNewPartial getElementsByTagName "Item"]
        set nodes [concat $nodes [$xmlNewPartial getElementsByTagName "Container"]]
        
        # Cada node aux es un nodo del documento auxiliar que utilizamos para recorrer por partes
        foreach nodeAux $nodes {
            
            #set attributes [list state dv open]
            # If an Item or Container is hidden at Default.spd, keep it hidden
            set attributes [list dv open]
            foreach att $attributes {
                
                #Si tiene dv actualizamos su valor
                if { [$nodeAux hasAttribute $att] } {
                    
                    #Buscamos el valor en el antiguo spd para actualizarlo
                    set id [$nodeAux getAttribute id ""]
                    set foundNode [$baseNodeOld find id $id]
                    if { $foundNode != "" &&  [$foundNode nodeName] == [$nodeAux nodeName]  } {
                        
                        #Hemos encontrado el mismo nodo en el modelo del problemtype
                        
                        set dvOld [$foundNode getAttribute $att ""]
                        set dvNew [$nodeAux getAttribute $att ""]
                        
                        if { $dvOld != $dvNew } {
                            
                            #Estamos en condiciones de actualizar el valor DV de xmlOld en xmlNew
                            set xPath [$nodeAux toXPath]
                            set xPath [string map [list "/RootData" $baseNodePath] $xPath]
                            set nodeNew [$xmlNew selectNodes $xPath]
                            
                            if {$nodeNew != "" } {
                                #msg "  .. . .  newid: [$nodeNew getAttribute id 0]"
                                #Actualizamos el valor de DV validando si es existen iValores
                                set ivalues [split [$nodeNew getAttribute ivalues ""] ","]
                                
                                if { [llength $ivalues] } {
                                    
                                    if { ($dvOld in $ivalues) } {
                                        $nodeNew setAttribute $att $dvOld
                                        if {$outputLog} {
                                            lappend ::xmlutils::logChanges "UPDATE $att (COMBO): $id  default $att: $dvNew  new $att: $dvOld"
                                            lappend ::xmlutils::logChanges "          path:[::xmlutils::printPathFromNode $nodeNew]\n"
                                        }
                                    }
                                } else {
                                    if {$outputLog} {
                                        lappend ::xmlutils::logChanges "UPDATE $att: $id  default $att: $dvNew  new $att: $dvOld"
                                        lappend ::xmlutils::logChanges "         path:[::xmlutils::printPathFromNode $nodeNew]\n"
                                    }
                                    $nodeNew setAttribute $att $dvOld
                                }
                            }
                            incr numLogs
                        }
                    }
                }
            }
        }
    }
    
    #------------------------------------------------------------------------
    # Insertar los grupos definidos y la asignación de grupos y propiedades
    #------------------------------------------------------------------------
    #--- GRUPOS
    set nodeGroups [$xmlOld selectNodes "/Kratos_Data/Groups/Group"]
    set newNodeGroups [$xmlNew selectNodes "/Kratos_Data/Groups"]
    
    foreach node $nodeGroups {
        $newNodeGroups appendChild $node
    }
    set modeltypenode [$xmlOld selectNodes "/Kratos_Data/Groups"]
    set modeltype [lindex [lindex [$modeltypenode asList] 1] 1]
    #::xmlutils::setXml [$xmlNew selectNodes "/Kratos_Data/Groups"] "modeltype" "write" $modeltype
    $xmlNew set "/Kratos_Data/Groups/@modeltype" "$modeltype"
    
    
    #--- ASIGNACION DE GRUPOS Y PROPIEDADES
    set baseNodePaths {}
    if {$KPriv(what_dempack_package) eq "C-DEMPack"} {
        lappend baseNodePaths "/Kratos_Data/RootData\[@id='DEM'\]"        
    } else {
        lappend baseNodePaths "/Kratos_Data/RootData\[@id='StructuralAnalysis'\]"
        #lappend baseNodePaths "/Kratos_Data/RootData\[@id='Fluid'\]"
        #lappend baseNodePaths "/Kratos_Data/RootData\[@id='PFEM'\]"
        #lappend baseNodePaths "/Kratos_Data/RootData\[@id='FluidStructureInteraction'\]" 
        #lappend baseNodePaths "/Kratos_Data/RootData\[@id='ConvectionDiffusion'\]" 
        lappend baseNodePaths "/Kratos_Data/RootData\[@id='DEM'\]"
    }
    
    foreach baseNodePath $baseNodePaths {
        
        lappend ::xmlutils::logChanges "\nGroups and properties ($baseNodePath):\n"
        set baseNodeOld [$xmlOld selectNodes $baseNodePath]
        
        #Recorremos el nuevo xml también por partes
        set baseNodeNew [$xmlNew selectNodes $baseNodePath]
        set xmlNewPartial [[dom parse [$baseNodeNew asXML]] documentElement]
        
        #Para cada container de clase groups, comprobamos si tiene grupos
        set nodes [$xmlNewPartial getElementsByTagName "Container"]
        #Cada node aux es un nodo del documento auxiliar que utilizamos para recorrer por partes
        foreach nodeAux $nodes {
            set class [$nodeAux getAttribute class ""]
            if { $class == "Groups" || $class == "Properties"} {        
                
                set idTemplate [$nodeAux getAttribute idTemplate ""]
                set xpath [::xmlutils::getPathFromNode $nodeAux]
                set nodes [$xmlOld selectNodes "${xpath}/Container"]
                foreach groupNode $nodes {
                    #msg "copy $xmlNew $groupNode $xpath $idTemplate"
                    ::xmlutils::copyGroupNode $xmlNew $groupNode $xpath $idTemplate $nodeAux [$nodeAux getAttribute class ""] 
                }
            }
        }
    }
    
    set KPriv(xml) $xmlNew
    # wa "xmlNew Final"
    # wa [$KPriv(xml) asXML]
    set KPriv(xmlDoc) $xmlDocNew
    set KPriv(encrXml) $encrXmlDef
    if {$outputLog} {
        #Crea un archivo de log con los cambios realizados
        set ptypeName "updateSPD.log"
        catch {
            set ptypeName [string map {".gid" ""} [lindex [split $KPriv(problemTypeDir) "/"] end] ]
        }
        
        # wa "llength logChanges: [llength $logChanges]\nptypename: $ptypeName"
        set fullname [file native [file join $KPriv(problemTypeDir) "${ptypeName}.log"]]
        set outfd [open $fullname w+]
        
        puts $outfd "#----------------------------------------------------------------"
        puts $outfd "# Changes in configuration file '.spd'"
        puts $outfd "# converting $pTypeVersion to $defaultVersion version."
        puts $outfd "#----------------------------------------------------------------"
        puts $outfd ""
        foreach line $logChanges {
            #La escribe en el fichero de salida
            puts $outfd $line
        }
        
        close $outfd
        
    }
}

proc ::xmlutils::copyGroupNode { xmlNew groupNode nodexPath idTemplate oldXmlNode {class "Groups"}} {
    
    
    set idGroup [$groupNode getAttribute id ""]
    
    set templxPath "/Kratos_Data/Templates/Template\[@id='$idTemplate'\]"
    
    set template_node [$xmlNew selectNodes $templxPath]
    
    ::KMProps::addSubtemplate $xmlNew $template_node
    
    set targetNode [$xmlNew selectNodes $nodexPath]
    
    $targetNode appendXML [$groupNode asXML]
    
    #msg [$groupNode asXML]
    
    set newNode [$xmlNew selectNodes "${nodexPath}/Container\[@id='$idGroup'\]"]
    #msg "New Nodes"
    foreach node [$newNode childNodes] {
        #msg [$node asXML]
        $node delete
    }
    
    foreach node [$template_node childNodes] {
        $newNode appendXML [$node asXML]
    }
    
    #msg [$newNode asXML]
    
    foreach node [$newNode childNodes] {
        ::xmlutils::GetOldDvFromNewNode $node
        foreach node2 [$node childNodes] {
            ::xmlutils::GetOldDvFromNewNode $node2
            foreach node3 [$node2 childNodes] {
                ::xmlutils::GetOldDvFromNewNode $node3
            }
        }
    }
    #msg [$newNode asXML]
    
    lappend ::xmlutils::logChanges "[$targetNode getAttribute id 0] --> $idGroup"
}


proc ::xmlutils::GetOldDvFromNewNode { node } {
    global KPriv
    
    if {[$node hasAttribute dv]} {
        
        set path [::xmlutils::myPathFromNode $node]
        set oldNodeDv [$KPriv(xml) selectNodes $path]
        
        if { $oldNodeDv != "" } {
            set dvOld [$oldNodeDv getAttribute dv ""]
            $node setAttribute dv $dvOld
        }
    }
}

proc ::xmlutils::myPathFromNode { finalNode {type "props"} } {
    
    set path ""
    set includeFinal 1
    set nodes [$finalNode ancestor all]
    #msg "$finalNode $nodes"
    
    set i 0
    foreach node $nodes {
        
        set id [$node getAttribute id ""]
        if { $id != "" } {
            set path "c.[list $id]//$path"
        }
    }
    
    if {$path != ""} {
        set path [string range $path 2 end]
        if { $includeFinal } {
            set path "${path}i.[$finalNode getAttribute id 0]"
        } else {
            set path [string range $path 0 end-2]
        }
        
        return [::xmlutils::setXPath $path $type]
    }
    return ""
}

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#                         ---- UPDATE SPD -----
# Compara la versión del problemTypeName.spd abierto con la versión de kratos_default.spd
# Y si son distintas añade al spd los nodos y atributos nuevos del default
#  Compara los atributos, añade los q falten y elimina los que sobren
#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
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
    
    # Si el nodo es un grupo o una propiedad, no estaba en default
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
                
                # Si el nodo original ya no lo tiene, es que se ha borrado
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
proc ::xmlutils::getPathFromNode { finalNode {includeFinal 1} {type "props"}} {
    
    set path ""
    
    set nodes [$finalNode ancestor all]
    set i 0
    foreach node $nodes {
        
        set id [$node getAttribute id ""]
        if { $id != "" } {
            set path "c.[list $id]//$path"
        }
    }
    
    if {$path != ""} {
        set path [string range $path 2 end]
        if { $includeFinal } {
            set path "${path}c.[$finalNode getAttribute id 0]"
        } else {
            set path [string range $path 0 end-2]
        }
        return [::xmlutils::setXPath $path $type]
    }
    return ""
    
}

proc ::xmlutils::getFullnameFromNode { finalNode {includeFinal 1} {type "props"}} {
    
    set path ""
    
    set nodes [$finalNode ancestor all]
    set i 0
    foreach node $nodes {
        
        set id [$node getAttribute id ""]
        if { $id != "" } {
            set path "c.[list $id]//$path"
        }
    }
    
    if {$path != ""} {
        set path [string range $path 2 end]
        if { $includeFinal } {
            set path "${path}c.[$finalNode getAttribute id 0]"
        } else {
            set path [string range $path 0 end-2]
        }
        return $path
    }
    return ""
    
}

proc ::xmlutils::getPathFromXPath { xpath {type "props"} } {
    
    set path ""
    
    set nodes [$finalNode ancestor all]
    
    set i 0
    foreach node $nodes {
        set id [$node getAttribute id ""]
        if { $id != "" } {
            set path "c.[list $id]//$path"
        }
    }
    
    if {$path != ""} {
        set path [string range $path 2 end]
        if { $includeFinal } {
            set path "${path}c.[$finalNode getAttribute id 0]"
        } else {
            set path [string range $path 0 end-2]
        }
        return [::xmlutils::setXPath $path $type]
    }
    return ""
    
}

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Dado un determinado nodo de TDom IMPRIME su path en las formas: 
# 1:  RootDataID//..IDnodoN..//IDnodoPadre//IDfinalNode
# 2:  RootDataID//..IDnodoN..//IDnodoPadre
#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

proc ::xmlutils::printPathFromNode { finalNode {includeFinal 1}} {
    
    set path ""
    
    set nodes [$finalNode ancestor all]
    
    set i 0
    foreach node $nodes {
        
        set id [$node getAttribute id ""]
        if { $id != "" } {
            set path "$id/$path"
        }
    }
    if {$path != ""} {
        
        if { $includeFinal } {
            set path "${path}[$finalNode getAttribute id 0]"
        } else {
            set path [string range $path 0 end-1]
        }
        return $path
    }
    return ""
    
}

proc ::xmlutils::parentNodePath { nodePath } {
    set parent [::KMProps::split2 $nodePath "/" ]
    set parent [lrange $parent 0 end-1]
    set ret ""
    foreach pc $parent {
        append ret $pc
        append ret "/"
    }
    set ret [string range $ret 0 end-2]
    return $ret
    
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
    
    # Caso especial ROOT_DATA
    if { [$parentNode nodeName] == "Kratos_Data" } {
        
        set foundParentNode $baseNode
        $foundParentNode appendChild $sourceNode
        
        lappend ::xmlutils::logChanges "INSERT NODE [$sourceNode nodeName] Id: [$sourceNode getAttribute id 0]"
        #msg "INSERT [$sourceNode nodeName] Id:[$sourceNode getAttribute id 0] xml:[$sourceNode asXML]"
        return ""
        
        # Caso especial TEMPLATES
    } elseif { [$parentNode nodeName] == "Templates" } {
        
        $parentNode appendChild $sourceNode
        lappend ::xmlutils::logChanges "INSERT NODE [$sourceNode nodeName] Id: [$sourceNode getAttribute id 0]"
        #msg "INSERT [$sourceNode nodeName] Id:[$sourceNode getAttribute id 0] xml:[$sourceNode asXML]"
        
        # El resto de nodos
    } else {
        
        # Si no es vacío ya sabemos donde meter "sourceNode"
        if { $foundParentNode != "" } {
            
            set nodoSiguiente [$sourceNode nextSibling]
            
            if {$nodoSiguiente != "" } {
                
                set idSiguiente [$nodoSiguiente getAttribute id ""]
                
                # Buscamos si también está este nodo donde lo vamos a insertar,
                # y si no está, lo pondremos el último
                set foundBro [$foundParentNode find id $idSiguiente]
                
                set insertOK 0
                catch {
                    # Intentamos insertarlo en su sitio
                    $foundParentNode insertBefore $sourceNode $foundBro
                    set insertOK 1
                }
                # Si ha fallado insertándolo en su lugar, lo ponemos el último
                if { !$insertOK } {
                    $foundParentNode appendChild $sourceNode
                }
                
                lappend ::xmlutils::logChanges "INSERT NODE [$sourceNode nodeName] Id: [$sourceNode getAttribute id 0]"
                
            } else {
                # Si no tiene nodo siguiente lo ponemos el último
                $foundParentNode appendChild $sourceNode
                lappend ::xmlutils::logChanges "INSERT NODE [$sourceNode nodeName] Id: [$sourceNode getAttribute id 0]"
                
            }
        }
    }
    return ""
}

#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#              #   #    #     #     #
#               # #     # # # #     #
#                #      #  #  #     #
#               # #     #     #     #
#              #   #    #     #     ######
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
    
    if {$path == ""} {
        return ""
    }
    set splitted [::KMProps::split2 $path //]
    
    if { $type != "props"} {
        set i 0
        foreach itemId $splitted {
            #wa "itemid splitted $i $itemId"
            if { $i == 0 } {
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
        foreach itemId $splitted {
            
            if { $i == 0 } {
                # El primer elemento será siempre del nivel 'RootData'
                set xpath "/Kratos_Data/RootData\[@id='$itemId'\]"
            } else {
                
                if { [string index $itemId 0] == "c" } {
                    set xpath "$xpath/Container\[@id='[string range $itemId 2 end]'\]"
                } elseif { [string index $itemId 0] == "i" } {
                    set xpath "$xpath/Item\[@id='[string range $itemId 2 end]'\]"
                } elseif { [string index $itemId 0] == "t" } {
                    set xpath "$xpath/ContainerTable\[@id='[string range $itemId 2 end]'\]"
                } elseif { [string index $itemId 0] == "T" } {
                    set xpath "$xpath/TItem\[@id='[string range $itemId 2 end]'\]"
                }
            }
            incr i
        }
    }
    return $xpath
}

proc ::xmlutils::getBackTrace backtraceref {
    upvar $backtraceref backTrace
    
    set startlevel [expr {[info level] - 2}]
    for {set level 1} {$level <= $startlevel} {incr level} {
        lappend backTrace [lindex [info level $level] 0]
    }
}

proc ::xmlutils::EvaluatePreProcessState {} {
    
    set backtrace {}
    ::xmlutils::getBackTrace backtrace
    set stackinfo $backtrace
    set after_executing False
    
    if {([lsearch $stackinfo ::AfterWriteCalcFileGIDProject*] != -1) || ([lsearch $stackinfo ::SelectGIDBatFile*] != -1)} {
        set after_executing True
    }
    
    return $after_executing
}

#
# Editar o extraer propiedades del xml en memoria
#

proc ::xmlutils::setXml {path property {command "read"} {value ""} {type "props"} {xpathvar 0}} {
    
    global KPriv
    
    set after_executing [::xmlutils::EvaluatePreProcessState]
    
    if {$type == "props"} {
        
        if {$xpathvar} {
            set xpath $path
        } else {
            set xpath "[::xmlutils::setXPath $path]"
        }
        
        if {$KPriv(what_dempack_package) ne "C-DEMPack"} {
            set check [$KPriv(xml) selectNodes $xpath]
            if {$command eq "check"} {
                return $check
            }
            if {([llength $check]==0) && ($command eq "read")} {
                W "props path $xpath does not exist"
                PrintStack
            }
        }
        
        if {$command == "read"} {
            
            set value [$KPriv(xml) set "$xpath/@$property"]
            
            if {$property == "dvText"} {
                
                # En vez de devolver el valor de dv, devuelve su equivalente traducible
                set value [$KPriv(xml) set "$xpath/@dv" ]
                if {[lindex $value 0] != ""} {
                    
                    set values [::xmlutils::getXMLValues $path]
                    set ivalues [::xmlutils::getXMLValues $path "" "iValues"]
                    
                    if {[llength $ivalues] > 0} { 
                        
                        set index [::xmlutils::getSelected $value $ivalues]
                        set value [lindex $values $index]
                    }
                }
            }
            
            # Cuando hay espacios el xml devuelve una lista y si la imprimes tal cual aparecen corchetes
            if {[llength $value] == 1} {
                set value [lindex $value 0]
            }
            
            if {($value == "") && ([string match "*\[@id='SetActive'\]*" $xpath])} {
                return $value
            }
            
            #to check if everything is OK and all fields are found when calculation is launched
            if {($value == "") && ($after_executing == True)} {
                W "path $path returns an empty value!"
            }
            
            return $value
            
        } else {
            $KPriv(xml) set "$xpath/@$property" "$value"
            #to check if everything is OK and all fields are found when calculation is launched
            if {$value == ""} {
                W "path $path returns an empty value $value for attribute $property!"
            }
            #return $value   
            return "1"
        }
        # Here we read the Material section in the tree
    } elseif { $type == "mat" } {
        set xpath "[::xmlutils::setXPath $path $type]"
        if {$KPriv(what_dempack_package) ne "C-DEMPack"} {
            set check [$KPriv(xmlMat) selectNodes $xpath]
            if { [llength $check]==0 } {
                W "material path $xpath does not exist" 
            }        
        }     
        
        if {$command == "read"} {
            set value [$KPriv(xmlMat) set "$xpath/@$property" ]
            
            if {$property == "dv"} {
                
                # En vez de devolver el valor de dv, devuelve su equivalente traducible
                set value [$KPriv(xmlMat) set "$xpath/@value" ]
                
                if {[lindex $value 0] != ""} {
                    
                    set values [::xmlutils::getXMLValues $path "" "" "" "" $type]
                    set ivalues [::xmlutils::getXMLValues $path "" "iValues" "" "" $type]
                    
                    if { [llength $ivalues] > 0 } { 
                        
                        set index [::xmlutils::getSelected $value $ivalues]
                        set value [lindex $values $index]
                    }
                }
            }
            
            # Cuando hay espacios el xml devuelve una lista y si la imprimes tal cual aparecen corchetes
            if { [llength $value] == 1 } {
                set value [lindex $value 0]
            }  
            
            #to check if everything is OK and all fields are found when calculation is launched
            #if...
            #    W "path $path returns empty value for attribute $property!"
            
            if {($value == "") && ($property == "dv")} {
                W "An empty value was found in path $path."
                W "There might be errors in your simulation."
            }
            return $value
        } else {
            $KPriv(xmlMat) set "$xpath/@$property" "$value"
            
            #to check if everything is OK and all fields are found when calculation is launched
            #if ...
            #    W "path $path returns empty value for attribute $property!
            return "1"
        }
    } elseif {$type == "curve" } {
        if { [kipt::CurvesModule ] } {
            set xpath "[::xmlutils::setXPath $path]"
            if { $command == "read" } {
                
                set value [$KPriv(xml) set "$xpath/@$property" ]
                
                if { $property == "dvText" } {
                    
                    # En vez de devolver el valor de dv, devuelve su equivalente traducible
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
                
                # Cuando hay espacios el xml devuelve una lista y si la imprimes tal cual aparecen corchetes
                if { [llength $value] == 1 } {
                    set value [lindex $value 0]
                }
                
                return $value
                
            } else {
                $KPriv(xml) set "$xpath/@$property" "$value"                            
                return "1"
            }
        }
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
    } elseif { $type == "mat" } {
        set nodes [$KPriv(xmlMat) selectNodes "${xpath}"]
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
    if {$type eq "props"} {
        set nodes [$KPriv(xml) selectNodes "${xpath}/$nodeType"]
    } elseif {$type eq "mats"} {
        set nodes [$KPriv(xmlMat) selectNodes "${xpath}/$nodeType"]
    }
    # wa "path:$path nodeType:$nodeType nodes:$nodes"
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
    
    # wa "path:$path id:$id atributo:$atributo tagname:$tagname type:$type"
    set listPares {}
    
    set xpath "[::xmlutils::setXPath $path $type]"
    
    if { $type == "mat"} {
        set xml $KPriv(xmlMat)
    } else {
        set xml $KPriv(xml)
    }
    
    set node [$xml selectNodes "${xpath}"]
    
    set items [$node descendant all $tagname]
    
    foreach node $items {
        # wa "node:$node"        
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
    #WarnWinText "xpath:$xpath node:$node"
    set cidlist [list]
    foreach node [$node descendant all Container] {
        lappend cidlist [$node getAttribute id ""]
    }
    #WarnWinText "cidlist:$cidlist"
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
#*          #####         #####  #          #  ####   #####
#*          #      #   #  # # # #  #   #  #   #
#*          #      #   #  #  #  #  ####   #   #
#*          #      #   #  #     #  #   #  #   #
#*          #####         #####  #     #  ####   #####
#*************************************************************
#
# Accede al XML al nodo con path 'fullname' y coge el 
# atributo "values" o la lista especial correspondiente
#

proc ::xmlutils::getXMLValues { fullname {idTemplate ""} {iValues ""} {idTemplateFull ""} {specialFilter ""} {type "props"}} {
    
    # wa "$application --> $comboList\nargs:1$fullname 2$idTemplate 3$iValues 4$idTemplateFull 5$specialFilter"
    
    
    global KPriv
    
    set comboList {}
    
    if { $iValues != "" } {
        set atrValues "ivalues"
    } else {
        set atrValues "values"
    }
    
    if { $idTemplate == "" } {
        # Se utiliza el fullname normalmente
        set specialList [::xmlutils::setXml $fullname GCV "read" "" $type]
    } else {
        # Se tiene que consultar la propiedad en el template
        set specialList [::KMProps::getPropTemplate $idTemplate GCV "$fullname"]
    }
    # Si no hay lista especial, cargamos el contenido de "values="
    if { $specialList == "" } {
        
        if { $idTemplate == "" } {
            set comboList [split [::xmlutils::setXml $fullname $atrValues "read" "" $type] ","]
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
        
        set node [$KPriv(xmlDocKKW) selectNodes "Kratos_KWords/ElementCLaws/Item\[@id='ElementTypes$::KMProps::nDim'\]"]
        #set node [$KPriv(xmlDocKKW) selectNodes "Kratos_KWords/ElementCLaws/Item\[@id='ElementTypes'\]"]
        
        if { $node != "" } {
            set comboList [split [$node getAttribute $atrValues ""] ","]
            # Set the dv value
            if {($idTemplateFull ne "")&&($specialFilter ne "")} {
                set ok [::xmlutils::setXml $idTemplateFull dv "write" "$specialFilter" $type]
            }
        }
        
    } elseif { $specialList == "SectType" } {
        
        set node [$KPriv(xmlDocKKW) selectNodes "Kratos_KWords/ElementCLaws/Item\[@id='SectionTypes'\]"]
        
        if { $node != "" } {
            set comboList [split [$node getAttribute $atrValues ""] ","]
            # wa "comboList:$comboList specialFilter:$specialFilter atrValues:$atrValues idTemplateFull:$idTemplateFull"
            if {($idTemplateFull ne "")&&($specialFilter ne "")} {
                set ok [::xmlutils::setXml $idTemplateFull dv "write" "$specialFilter" $type]
            }
        }
        
    } elseif { $specialList == "MatModel" } {
        
        set dvElemFilter ""
        if { $specialFilter == "" } {
            # Caso especial para los combos dinámicos:
            # Mira a ver si el elemento anterior era ElemType para utilizar su valor
            set elemtypeFullname [string map {"MatModel" "ElemType"} $fullname]
            set dvElemFilter [::xmlutils::setXml $elemtypeFullname dv]
        }
        # Caso estándard: accede al valor dv de element type para coger el matModel filtrando por este elementType
        if { $dvElemFilter == "" } {
            # Si el valor no existe, miramos el que hay en el xml KKWORDS
            set node [$KPriv(xmlDocKKW) selectNodes "Kratos_KWords/ElementCLaws/Item\[@id='ElementTypes$::KMProps::nDim'\]"]
            #set node [$KPriv(xmlDocKKW) selectNodes "Kratos_KWords/ElementCLaws/Item\[@id='ElementTypes'\]"]
            if { $node != "" } {
                set dvElemFilter [$node getAttribute dv ""]
            }
        }
        # Cogemos la lista correspondiente
        set node [$KPriv(xmlDocKKW) selectNodes "Kratos_KWords/ElementCLaws/Item\[@id='$dvElemFilter'\]"]
        if { $node != "" } {
            set comboList [split [$node getAttribute $atrValues ""] ","]
        }
    } elseif { $specialList == "ProfileType" } {
        
        set dvSectFilter ""
        if { $specialFilter == "" } {
            # Look if the previous section type to use their value
            set sectiontypeFullname [string map {"ProfileDB" "SectionType"} $fullname]
            set dvSectFilter [::xmlutils::setXml $sectiontypeFullname dv]
            # wa "sectiontypeFullname:$sectiontypeFullname dvSectFilter:$dvSectFilter"
        }
        # Standard case: access to the value of dv in the section type to get the profiledb filter by this section type
        if {$dvSectFilter == "" } {
            # If the values does not exist, get it from xml KKWORDS
            set node [$KPriv(xmlDocKKW) selectNodes "Kratos_KWords/ElementCLaws/Item\[@id='SectionTypes'\]"]
            if {$node != "" } {
                # Get dv value
                set dvSectFilter [$node getAttribute dv ""]
            }
        }
        # Get the correct list
        set node [$KPriv(xmlDocKKW) selectNodes "Kratos_KWords/ElementCLaws/Item\[@id='$dvSectFilter'\]"]
        if { $node != "" } {
            set comboList [split [$node getAttribute $atrValues ""] ","]
        }
        # wa "comboList:$comboList"
    }
    # Si sacamos los ids no los traducimos
    if { $iValues == "" } {
        set resultList {}
        foreach val $comboList {
            set traduction "[= $val ]"
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
        if {[winfo exists $fcmb] } {
            return [$fcmb get]
        }
        msg "Fuera de índice $index (length: [llength $comboList])"
    }
}

proc ::xmlutils::setComboDv { fcmb fullname dv {idTemplate ""} } {
    
    set icomboList [::xmlutils::getXMLValues $fullname "$idTemplate" "iValues"]
    set selected [::xmlutils::getSelected $dv $icomboList]
    # wa "icomboList:$icomboList selected:$selected"
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
#        * xmlstruct::create xml - Returns an extended domNode object command
# The following methods are supported in addition to those the tDOM domNode object already provides:
#        * $node set xpathQuery ?value? - retrieves or modifies portions of the xml document that match the given xpathQuery
#        * $node unset xpathQuery - deletes the portions of the xml document that match the given xpathQuery
#        * $node lappend xpathQuery value ?value? ?value? ... - appends the given values to the node(s) that match the given xpathQuery

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
    set aux [$node selectNodes $query type]
    ::set resultNodes $aux
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
    
    # wa "node:$node query:$query value:$value" 
    
    #hacer algo con xpath
    #al hacer el split se destruye todo ya que query contiene "//" (carpeta grupos)
    
    #set possibleMatch [split $query /]
    set list_splited [split $query \]]
    # wa "list_splited:$list_splited"
    if { [lindex $list_splited end] == "" } {
        #last element ends with \]
        set unmatched [lindex $list_splited end-1]
        set unmatched "[string range $unmatched 1 end]\]"
    } else {
        set unmatched [string range [lindex $list_splited end] 1 end]
    }
    set possibleMatch [string range $query 0 end-[expr [string length $unmatched]+1]]
    if {[llength $possibleMatch] == 0} {
        set possibleMatch .
    }
    # wa "unmatched:$unmatched possibleMatch:$possibleMatch"
    set nodes [$node selectNodes $possibleMatch type]
    # wa "nodes:$nodes"
    switch $type {
        nodes {            
            if {[string index $unmatched 0] == "@"} {
                foreach node $nodes {
                    $node setAttribute [string range $unmatched 1 end] $value
                }
            } else {
                if { [regexp -start 0 -indices {\[} $unmatched idxs] } {
                    set nametag [string range $unmatched 0 [expr [lindex $idxs 0]-1]]
                    set unmatched [string range $unmatched [expr [lindex $idxs 0]+1] end-1]
                } else {
                    set nametag $unmatched
                    set unmatched ""
                }
                # wa "nametag:$nametag unmatched:$unmatched"
                foreach node $nodes {
                    if {[string index $unmatched 0] == "@"} {
                        set child [[$node ownerDocument] createElement "$nametag"]
                        $child setAttribute {*}[regsub -all ' [regsub -all =' [string range $unmatched 1 end] " \{"] \}]
                    }                    
                    $node appendChild $child
                    #appendXML "<$unmatched/>"
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
#          #   #       #     #       #
#           # #        # # # #       #
#            #         #  #  #       #
#           # #        #     #       #
#          #   #       #     #       ######
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
#

#
# Si no hay atributo devuelve uno o más nodos, y si lo hay, uno o más atributos.
# Además, si value no es nulo, hace un set en el xml
#
proc ::xmlutils::getAttribute { xml xPath {attribute ""} {value ""} } {
    set nodes ""
    if { $xml != "" } {
        set nodes [$xml selectNodes $xPath]
    }
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
    } else {
        return "No_node_in_xml"
    }
}

#
# Accede al XML al nodo con xpath y coge el 
# atributo "values" o la lista especial correspondiente
#
proc ::xmlutils::getValues { xml xpath {iValues ""} } {
    #msg "Getting values from /n [$xml asXML] /n on $xpath with $iValues"
    set comboList {}
    
    if { $iValues != "" } {
        set atrValues "ivalues"
    } else {
        set atrValues "values"
    }
    
    set specialList [::xmlutils::getAttribute $xml $xpath GCV]
    #msg "Special list $specialList"
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
            #set traduction "[= $val ]"
            set traduction "$val"
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
    if {$attribute ni {pid help}} {
        if { [lindex $value 0] != "" } {
            
            set values [::xmlutils::getValues $xml $xpath]
            set ivalues [::xmlutils::getValues $xml $xpath "iValues"]
            
            if { [llength $ivalues] > 0 } { 
                
                set index [::xmlutils::getSelected $value $ivalues]
                set value [lindex $values $index]
            }
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
    
    set template [$xml set "${templatePath}\[@id=\'$idTemplate\'\]"]
    
    #No se puede insertar en el xml un fragmento con mas de un nodo, por eso utilizamos 
    #las etiquetas auxiliares "<grouptemplate>$template</grouptemplate>"
    set template "<groupTemplate>$template</groupTemplate>"
    # wa "template:$template"
    set textAttr ""
    foreach attrValue $attributesArray {
        set textAttr "$textAttr $attrValue"
    }
    # wa "xpath/nodeName:$xpath/$nodeName textAttr:$textAttr"
    set firstattr [lindex $textAttr 0]
    # wa "firstattr:$firstattr"
    $xml lappend "$xpath/$nodeName\[@id=\'$firstattr\' $textAttr\]" $template
    
    return $template
}

proc ::xmlutils::replaceTemplate { xml {xpath ""} } {
    
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

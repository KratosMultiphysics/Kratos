###############################################################################
#
#  NAME: kmprops.tcl
#
#  PURPOSE: Main window to manage model properties
#
#  QUANTECH ATZ-DEVELOPMENT DEPARTMENT
#
#  AUTHORS : L. Calvo, G. Socorro
#
#  CREATED AT: 25/02/2010
#
#  HISTORY:
# 
#   1.9- 16/07/13-G. Socorro, modify the proc stateNode to show the thickness as a function of the selected element type
#   1.8- 06/05/13-G. Socorro, add the option to work with the cross section properties (simple property or database)
#   1.7- 09/10/12-G. Socorro, update others procs to include the cross property functionality
#   1.6- 20/09/12-J.Garate, Adaptation for New Kratos Interface Version, including Curves support
#   1.5- 04/04/12-J.Garate, finalizada la implementacion de la funcion recursiva de llenado del arbol desde xml, ya en uso
#   1.4- 02/04/12-J.Garate, add ::KMProps::FillTreePropsRecursive , not in use yet until validation.
#   1.3- 28/03/12-G. Socorro, pass some procs to the tcl file kmpropswin.tcl
#   1.2- 27/03/12-G. Socorro, update the procs GetAutomaticPropertyName and buildPropertyFrame to be used in the fluid application
#   1.1- 22/03/12-J. Gárate, Cambio a funciones públicas de los grupos de GiD
#   1.1- 20/03/12-J. Garate, Renombrado de grupos y actualizacion del arbol tras el evento de renombrado
#   1.0- 19/03/12- J. Garate, Borrado de grupos y actualizacion del arbol tras el evento de borrado
#   0.9- 15/03/12 J.Garate: Reparado el boton Auto Group, que aparecia fuera del frame cuando la ventana estaba INSIDE 
#   0.8- 15/03/12 J.Garate: Solucionado el bug de la ventana INSIDE LEFT, que aparecia minimizada. Reestructurado de una funcion ocndicional.
#   0.7- 12/03/12 J.Garate: Total adaptacion de los grupos antiguos a GiD Groups. Arreglado el bug de los botones del AutoGroup
#   0.7- 31/01/12 J.Garate: Arreglado el bug del "Doble Click"
#   0.6- 24/05/11-G. Socorro, correct a bug with the proc InsertCaption when the GiD layer window is open
#   0.5- 12/08/10 Miguel Pasenau: se puede integrar la ventana dentro de la ventana 
#                 principal de GiD;
#                 mediante una barrita (caption) debajo del titulo de la ventana se puede
#                 seleccionar si se quiere integrar la ventana o no;
#                 tambien se ha cambiado el empaquetamiendo de 'pack' a 'grid'
#                 de momento esta puesto OUTSIDE / INSIDE_LEFT, 
#                 nueva variable Location para indicar donde esta
#   0.4- 05/05/10 LCA: Se han habilitado combos con idioma añadiendo al spd la propiedad "ivalues"
#   0.3- 26/04/10 LCA: Estructura correcta del árbol, con filtros por dimensión(2D,3D) 
#                y por Types (Structural, Solution y Analysis). Funciona link entre propiedades y materiales.
#   0.2- 29/03/10 LCA: Hasta este momento tenemos funcionando el tab Properties
#                (el de Materials está en proceso). Se pueden editar items y asignar grupos a condiciones.
#   0.1- 25/02/10 LCA, create a base source code from the kegroups.tcl script
#
###############################################################################

# Package requires
package require treectrl
package require tooltip
package require tdom
package require xmlstruct 1.0

package provide KMprops 1.0

# Create a base namespace KMProps
namespace eval ::KMProps:: {
    
    # Properties array
    variable Props
    # Path of the base window 
    variable WinPath ".gid.kmprops"
    # Window layout ["OUTSIDE"|"LEFT"|"RIGHT"] 
    variable WinLayout 
    # Location: OUTSIDE | INSIDE_LEFT | INSIDE_RIGHT
    variable Location         
    variable SystemHighlight
    variable SystemHighlightText

    # Tree properties path
    variable TreePropsPath
    # Notebook properties path
    variable NbPropsPath
    
    variable ngroups
    # Last selected item
    variable lastSelected {}
    variable selectedEntity ""
    
    variable application "StructuralAnalysis"
    
    # Se inicializan las clases dinámicamente leyendo del xml
    variable visibilityVars {}
    # Maximum number of iterations to find a new property
    variable MaxIdIter

    # Window path name
    variable winpath ".gid.kmprops"
    
    # Window position
    variable StartLayout "INSIDE_LEFT"
   
    # Current window position ["INSIDE_LEFT" "INSIDE_RIGTH" "OUTSIDE"
    variable Layout $StartLayout
    
    # 1 If need to restore when switch from/to postProcess
    variable RestoreWinFromPost

}

proc ::KMProps::Init {} {
	
    variable WinLayout
    variable SystemHighlight
    variable SystemHighlightText
    variable ngroups
    variable Props
   
    # Get default colors
    set w [listbox .listbox]
    set SystemHighlight [$w cget -selectbackground]
    set SystemHighlightText [$w cget -selectforeground]
    destroy $w
    if { $::tcl_platform(platform) eq "unix" } {
	# I hate that gray selection color
	set SystemHighlight #316ac5
	set SystemHighlightText White
    }
    
    set WinLayout "OUTSIDE"
    set ngroups 10000
    # Maximum number of iterations to find a new property
    variable MaxIdIter
    set MaxIdIter "1000"
    variable RestoreWinFromPost
    set RestoreWinFromPost 0
}

proc ::KMProps::initVisibilityClass { } {
    variable visibilityVars
    global KPriv

    set visibilityVars {}

    set classes [$KPriv(xml) set "/Kratos_Data/ClassConfiguration/Class" ]
    
    set classes [split $classes ","]
    
    foreach classid $classes {
	lappend visibilityVars $classid
	set ::KMProps::${classid} ""
    }
    
    # wa "visibilityVars:$visibilityVars\nkinemtype:[info exists ::KMProps::kinemType]"
}

#---------------------------------------------------------------------------------------------- 
# Lee el xml y carga el árbol de propiedades de forma recursiva, sin limite de niveles
#----------------------------------------------------------------------------------------------
proc ::KMProps::FillTreeProps { } {
	
	variable TreePropsPath
	
	set T $TreePropsPath

	# Seleccionamos todos los nodos del primer nivel
	set nodes [$::KPriv(xml) selectNodes "/Kratos_Data/RootData\[@id\]"]
	
	foreach node $nodes {
	    # Insertamos cada RootData de 1er nivel en el árbol
	    set item [::KMProps::InsertNewProp $node [$node getAttribute id ""] $T "" "root" [$node hasChildNodes] [::KMProps::stateNode $node] [$node getAttribute open "0"]]
	    # wa "item:$item"
	    if {$item != -1} {
		set path "[$node getAttribute id 0]//"
		::KMProps::FillRecursiveChilds $T $path $node $item
	    }
	}
    }

proc ::KMProps::FillRecursiveChilds { T path node item} {
	
	set nodes2 [$node childNodes]
	set pathcp $path
	foreach node2 $nodes2 {
	    # msg $path
	    # msg $pathcp
	    # msg "id:[$node2 getAttribute id ""]"
	    #set item3 [::KMProps::InsertNewProp $node3 [::KMProps::splitNode $node3] $T "[$node getAttribute id 0]//[::KMProps::splitNode $node2]//" "$item2" [$node3 hasChildNodes] [::KMProps::stateNode $node3] [$node3 getAttribute open "0"]]
	    set item2 [::KMProps::InsertNewProp $node2 [::KMProps::splitNode $node2] $T $path "$item" [$node2 hasChildNodes] [::KMProps::stateNode $node2] [$node2 getAttribute open "0"]]
	    if {$item2 != -1} {
		
		#set pathcp [join [concat $path "[::KMProps::splitNode $node2]//"] ""]
		set pathcp [regsub -all "// " [concat $path "[::KMProps::splitNode $node2]//"] "//"]
		
		#msg "$pathcp \n"
		::KMProps::FillRecursiveChilds $T $pathcp $node2 $item2
		#set pathcp [join [concat $path "[::KMProps::splitNode $node2]//"] ""]
		set pathcp [regsub -all "// " [concat $path "[::KMProps::splitNode $node2]//"] "//"]
	    }
	}
    }

#
# Valida varias cosas para cada nodo
#  DIMENSION (2D / 3D)
#  STATE (normal, hidden, disabled)
#
proc ::KMProps::stateNode { node } {
	
    # Validate for each node if it has to be visible
    # (depending on the chosen values in some comboboxes)
    if { [$node nodeName] == "Item" } {
	# Get the cross section property list
	set PropertyList [::KMProps::GetCrossSectionPropertyList]
	# wa "PropertyList:$PropertyList"
	
	# Special case for the properties (cross section properties)
	set id [$node getAttribute id ""]
	# wa "\nid:$id\n"
	if { $id == "ElemType" } {
	    #  Element type property
	    
	    set ::KMProps::ElemTypeProperty [$node getAttribute dv ""]
	    
	} elseif {$id in $PropertyList} {
	    # Get the parent node
	    set parentNode [[$node parentNode] parentNode]
	    set CurrentPropertyId [$parentNode getAttribute id ""]
	    # wa "CurrentPropertyId:$CurrentPropertyId parentNode:$parentNode"
	    if {$CurrentPropertyId !=""} {
		set path [::xmlutils::getFullnameFromNode $parentNode 0]
		# wa "path:$path"
		# Get the select element base type
		set ::KMProps::ElemTypeProperty [::xmlutils::GetPropertyElemType $path $CurrentPropertyId]
		# wa "ElemTypeProperty:$::KMProps::ElemTypeProperty"
	       
		# Get the select section base type
		set ::KMProps::SectionTypeProperty [::xmlutils::GetPropertySectionType $path $CurrentPropertyId]
		# wa "SectionTypeProperty:$::KMProps::SectionTypeProperty"
	    }
	    # Special case of thickness
	    if {$id eq "Thickness"} {
		# Check that this cross section property is active
		set ShowProperty [::KMProps::ShowPropertyByElementType $id]
		# wa "ShowProperty:$ShowProperty"
		if {!$ShowProperty} {
		    return "hidden"
		} 
	    } else {
		# Check that this cross section property is active => Checking the section type
		set SectionTypeShowProperty [::KMProps::ShowPropertyBySectionType $id RTree]
		# wa "SectionTypeShowProperty:$SectionTypeShowProperty"
		if {!$SectionTypeShowProperty} {
		    return "hidden"
		}
	    } 
	}
    }
    
    set state [$node getAttribute state "normal"]
    # wa "\nid:[$node getAttribute id ""]        state:$state"
    
    #Si el estado es hiddenAll se oculta el nodo y toda su descendencia
    if {$state == "hiddenAll"} {
	return "-1"
    }
    
    #if { [$node nodeName] == "Item" } {
    #        #Salvedad para no mostrar la propiedad Thickness en algunos casos
    #        set id [$node getAttribute id ""]
    #        if { $id == "PressureValue" } {
    #                set id $id
    #        }
    #}
    
    # wa "::KMProps::visibilityVars:$::KMProps::visibilityVars"
    foreach var $::KMProps::visibilityVars {
	
	set globalVar [set ::KMProps::$var]
	
	set nodeValuesVar [split [$node getAttribute $var ""] ","]
	
	# wa "id:[$node getAttribute id ""] var:$var globalVar:$globalVar nodeValuesVar:$nodeValuesVar"
	#Si el nodo tiene alguna restriccion de clase (p.ej. del tipo strucType=Shell)
	# y no coincide con el valor seleccionado, ocultamos el nodo        
	 if { $nodeValuesVar != "" } {
	    
	    if { !($globalVar in $nodeValuesVar) } {
		
		# msg "$state:  --------- > g: $globalVar in nodeVals: $nodeValuesVar"
		if {$var == "strucType" } {
		    if { $globalVar != "Generic"} { 
		        # msg "aqui $nodeValuesVar \"\" !=  && $globalVar in $nodeValuesVar"
		        ::KMProps::setNoActiveGroups $node
		        return -1
		    }
		} else {
		    # msg "este caso $nodeValuesVar \"\" !=  && $globalVar in $nodeValuesVar"
		    ::KMProps::setNoActiveGroups $node
		    return -1
		}
	    }
	}
    }
    
    #Leemos la class del nodo para ver si requiere de acciones especiales (ocultar nodos)
    #Si no se llega aquí es porque el nodo no era visible, y en ese caso no se tiene que hacer
    set value [$node getAttribute dv ""]
    
    set class [$node getAttribute class ""]
    # wa "value:$value class:$class"
    #Equivalente a Switch $class
    foreach var $::KMProps::visibilityVars {
	
	if {$var == $class} {
	    
	    #Caso especial para el solver de fluidos
	    if { $var == "fluidSolvTyp" } {
		
		#El solver de fluidos está duplicado dependiendo de una variable previa
		set freeYesOrNo [$node getAttribute freeSurf ""]
		
		if { $freeYesOrNo == "" || $freeYesOrNo == $::KMProps::freeSurf } {
		    set ::KMProps::fluidSolvTyp "$value"
		}
	    } else {
		# General case
		set ::KMProps::$var $value
	    }
	}
    }
    
    if { $class == "Group" } {
	# If you come here to say that the node is visible
	if { $state == "normal" } {        
	    $node setAttribute active 1
	}
    }
    # wa "state:$state"
    
    # Back node status ("normal" by default)
    return $state
}


#
# Separa cada node en "inicialNombreTag.idNodo"
#
proc ::KMProps::splitNode { node } {
	
	set id [$node getAttribute id ""]
	
	if { [$node tagName] == "Container"} {
		return "c.[list $id]"
	} elseif { [$node tagName] == "Item"} {
		return "i.[list $id]"
	} elseif { [$node tagName] == "Property"} {
		return "p.[list [$node getAttribute id ""]]"
	} elseif { [$node tagName] == "Material"} {
		return "m.[list [$node getAttribute id ""]]"
	} elseif { [$node tagName] == "ContainerTable"} {
		return "t.[list $id]"
	} elseif { [$node tagName] == "TItem"} {
		return "T.[list $id]"
	} else {
		return "NoTree"
	}
}


#
###################################################################################################
#--------------------------------------------------------------------------------------------------
# Funciones para manejar el XML de propiedades (las genéricas están en /lib/xml/xmlutils.tcl
#--------------------------------------------------------------------------------------------------
###################################################################################################
#
# Prepara la query para utilizar las funciones de domNOde
proc ::KMProps::setXPath { path  } {
	
	set splitted [::KMProps::split2 $path //]
	
	set i 0
	
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
	
	return $xpath
}

proc ::KMProps::getTemplateStructure { id {type "containers"}} {
	
	global KPriv
	
	set T $::KMProps::TreePropsPath
	
	set listTemplate {}
	
	if { $type == "containers" } {
		
		#Seleccionamos todos los containers del template y sus respectivos items
		set nodes [$KPriv(xml) selectNodes "/Kratos_Data/Templates/Template\[@id='$id'\]/Container"]
		
		foreach node $nodes {
		        set listChilds [$node getAttribute id ""]
		        foreach item [$node childNodes] {
		        
		            set listChilds [lappend listChilds [$item getAttribute id ""]]
		        }
		        set listTemplate [lappend listTemplate $listChilds]
		}
	} else { # "Items"
		set listTemplate [$KPriv(xml) selectNodes "/Kratos_Data/Templates/Template\[@id='$id'\]/Items"]
	}
	
	return $listTemplate
}

proc ::KMProps::copyTemplate { idTemplate fullname groupId clase } {
	
	global KPriv
	
	#set idTemplate "Forces"
	set id [::KMProps::getPropTemplate "$idTemplate" id]
	set pid [::KMProps::getPropTemplate "$idTemplate" pid]
	set help [::KMProps::getPropTemplate "$idTemplate" help]
	set icon [::KMProps::getPropTemplate "$idTemplate" icon]
	
	set xpath "[::xmlutils::setXPath $fullname]"
		
	set template_node [$KPriv(xml) selectNodes "/Kratos_Data/Templates/Template\[@id='$idTemplate'\]"]
	
	::KMProps::addSubtemplate $KPriv(xml) $template_node
	
	set template [$KPriv(xml) set "/Kratos_Data/Templates/Template\[@id='$idTemplate'\]"]
	
	#No se puede insertar en el xml un fragmento con mas de un nodo, por eso utilizamos 
	#las etiquetas auxiliares "<grouptemplate>$template</grouptemplate>"
	    set template "<groupTemplate>$template</groupTemplate>"
	
	if { $clase == "OnlyGetText" } {
		return $template
	}
	
	    $KPriv(xml) lappend "$xpath/Container\[@id='$groupId' pid='$groupId' class='$clase' icon='$icon' help='$help' open='1'\]" $template
	
	    #Eliminamos ahora las etiquetas auxiliares "<gouptemplate></gouptemplate>"
	    ::KMProps::replaceTemplate
	
	return $template
}

proc ::KMProps::addSubtemplate { xml_to_edit template_node } {

        foreach subnode [$template_node childNodes] {            
            if {[$subnode hasAttribute idTemplate]} {
                set subtemplate_name [$subnode getAttribute idTemplate] 
                set subtemplate_node [$xml_to_edit selectNodes "/Kratos_Data/Templates/Template\[@id='$subtemplate_name'\]"]
                set copy_of_subtemplate_node [$subtemplate_node cloneNode -deep]
                foreach container_in_subtemplate_element [$copy_of_subtemplate_node childNodes] {
                    $template_node insertBefore $container_in_subtemplate_element $subnode                                                                     
                 }
                $template_node removeChild $subnode
            } elseif {[$subnode nodeName] == "Container"} {
                set subnode_name [$subnode getAttribute id]
                foreach subsubnode [$subnode childNodes] {
                    if {[$subsubnode hasAttribute idTemplate]} {
                        set subtemplate_name [$subsubnode getAttribute idTemplate] 
                        set subtemplate_node [$xml_to_edit selectNodes "/Kratos_Data/Templates/Template\[@id='$subtemplate_name'\]"]                                                                                 
                        set copy_of_subtemplate_node [$subtemplate_node cloneNode -deep]
                        foreach container_in_subtemplate_element [$copy_of_subtemplate_node childNodes] {
                            $subnode insertBefore $container_in_subtemplate_element $subsubnode                                                                     
                        }
                        $subnode removeChild $subsubnode
                    }
                }
            }
        }
}

proc ::KMProps::replaceTemplate { } {
	
	global KPriv
	
	set xmlText [$KPriv(xml) asXML]
	
	set xmlText [string map {"<groupTemplate>" "" "</groupTemplate>" ""} $xmlText]
	
	set KPriv(xmlDoc) [dom parse $xmlText]
	
	set KPriv(xml) [$KPriv(xmlDoc) documentElement]
}

proc ::KMProps::getPropTemplate {id property {templatePath ""}} {
	
	global KPriv
	
	set xpath "/Kratos_Data/Templates/Template\[@id='$id'\]"
	
	set splitted [::KMProps::split2 $templatePath //]
	if {[llength $splitted] >= 1} {

		set idContainer [lindex $splitted 0]
		if { $idContainer != "OnlyItems" } {
		        set xpath "$xpath/Container\[@id='$idContainer'\]"
		}
	}
	if { [llength $splitted] >= 2 } {
		set idItem [lindex $splitted 1]
		set xpath "$xpath/Item\[@id='$idItem'\]"
	}

	set value [$KPriv(xml) set "$xpath/@$property" ]
	
	return [lindex $value 0]
	
}

proc ::KMProps::itemType {fullname} {
	
	set splitted [::KMProps::split2 $fullname //]
	
	if { [llength $splitted] < 2 } {
		return "r"
	} else {
		set letra [string index [lindex $splitted end] 0]
		#Devuelve identificador de item (en esta caso "c" o "i" de container o item)
		return $letra
	}
}


# - END XML FUNCTIONS

proc ::KMProps::getApplication { fullname } {
	
	set application [lindex [split $fullname "//"] 0]
	
	return $application
}

proc ::KMProps::getProps { fullname } {
	
	global KPriv
	
	set application [::KMProps::getApplication $fullname]
	
	#Miramos si hay alguna propiedad dada de alta (si la hay,obligatoriamente estará seleccionada)
	set props [::xmlutils::setXmlContainerIds "${application}//c.Properties"]
	
	#Miramos si es necesario filtrar por Constitutive Laws
	
	#Extraemos del fullname, el id del elemento pulsado
	set lSplit [split $fullname "//"]
	set lRange [lindex $lSplit 4] 
	set id [string range $lRange 2 end]
	
	#msg "$fullname \n $lSplit \n $lRange \n $id \n\n"
	
	set cLawsList {}
	set xPath "Kratos_KWords/ElementValidCLaws/Item\[@id='$id'\]"
	set active [::xmlutils::getAttribute $KPriv(xmlDocKKW) $xPath "active"]
	if { $active == 1 } {
		set cLawsList [split [::xmlutils::getAttribute $KPriv(xmlDocKKW) $xPath claws] ","]   
		
		if { $cLawsList != "" } {
		        
		        set FilterProps {} 
		        foreach idProp $props {
		        
		        set dv [::xmlutils::setXml "${application}//c.Properties//c.[list ${idProp}]//c.MainProperties//i.MatModel" dv]
		        
		        if {$dv in $cLawsList } {
		                
		                lappend FilterProps $idProp
		        }
		        }
		        return $FilterProps
		}
	}
	
	return $props
}

#separator can be a multicharacter, like //
proc ::KMProps::split2 { x separator } {
	set splitchar \uFFFF ;#forbidden unicode character, x must never contain it
	return [split [string map "$separator $splitchar" $x] $splitchar]
}

proc ::KMProps::randomColor { } {

	set r [expr {int (255 * rand())}]
	set g [expr {int (255 * rand())}]
	set b [expr {int (255 * rand())}]                
	
	return [format "\#%02x%02x%02x" $r $g $b]
}

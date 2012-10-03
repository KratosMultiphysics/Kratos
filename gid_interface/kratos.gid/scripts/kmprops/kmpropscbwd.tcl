#####################################################################################
#
#  NAME: kmpropscbwd.tcl
#
#  PURPOSE: Combo box widget option (edit,delete, etc.) from the kratos main window 
#           properties
#
#  QUANTECH ATZ-DEVELOPMENT DEPARTMENT
#
#  AUTHORS : G. Socorro => GSM
#
#  CREATED AT: 29/03/2012
#
#  HISTORY:
# 
#   0.5- 03/10/12- GSM, add a message for the "Compressible" fluid case
#   0.4- 27/09/12- J.Garate, Change combo's size
#   0.3- 23/09/12- GSM, update the proc specialComboAction to disabled FSI application
#   0.2- 23/07/12- GSM, modify some proc to use PFEM options 
#   0.1- 29/03/2012 G. Socorro, create a base source code from the kmprops.tcl script
#
######################################################################################
#                      Procedures that belong to this file
###############################################################################
#         Name                      |        Functionality
#------------------------------------------------------------------------------
# 1.            | 

proc ::KMProps::changeCmbValues { f path {idTemplate ""} {elemTypeDv ""} {noTemplateFullname ""} } { 
    
    set values [::xmlutils::getXMLValues $path $idTemplate "" "$noTemplateFullname" "$elemTypeDv"]
    
    $f configure -values $values
    
    if { [llength $values] > 0 && [$f current] == -1 } {
	$f current 0        
    }
}

proc ::KMProps::cmbElemTypeChange { f fullname {idTemplate ""} } {
    
    global KPriv
    
    if { [info exists ::KMProps::cmbElemType] } {
	
	set dv [::xmlutils::getComboDv $f $fullname "id" $idTemplate]
	
	set xpath "Kratos_KWords/ElementCLaws/Item\[@id='ElementTypes'\]"
	::xmlutils::getAttribute $KPriv(xmlDocKKW) $xpath dv $dv
	
	#Actualizamos los valores de Material Mode en función del nuevo dv
	set fMat [string map {"ElemType" "MatModel"} $f]
	
	if { $idTemplate == "" } {
	    
	    set fullnameMat [string map {"ElemType" "MatModel"} $fullname]
	    ::KMProps::changeCmbValues $fMat "$fullnameMat" $idTemplate "NoSearchFullname"
	} else {
	    
	    ::KMProps::changeCmbValues $fMat "MainProperties//MatModel" $idTemplate
	}
	
	set fCmbThick [string map {"ElemType" "Thickness"} $f]
	set fLblThick [string map {"cmbThickness" "lblThickness"} $fCmbThick]
	if { [::KMProps::showThickness] } {
	    
	    if { [winfo exists $fCmbThick] } { 
		grid $fCmbThick 
	    }
	    if { [winfo exists $fLblThick] } { 
		grid $fLblThick 
	    }
	} else {
	    if { [winfo exists $fCmbThick] } { 
		grid remove $fCmbThick 
	    }
	    if { [winfo exists $fLblThick] } { 
		grid remove $fLblThick 
	    }
	    
	}
    }
}

#
# Deshabilita el widget de values si el flag activation está a 0
#
proc ::KMProps::cmbDisable { fullname f {id "" }} {
    
    if {$id == "" } {
	set iniIdEmpty ""
	set id [::xmlutils::setXml $fullname id]
    } else {
	set iniIdEmpty "$id"
    }
    
    set selCombo [set ::KMProps::cmb$id]
    
    set whatDisable [string map {A V} $id]
    
    set fullname [string map [list "Activation//i.$id" "Values//i.$whatDisable"] $fullname]
    
    #msg "state?[$f.fValues.cmb$whatDisable  "
    if {$selCombo == 0 } {
	
	if {[winfo exists $f.fValues.cmb$whatDisable]} {
	    $f.fValues.cmb$whatDisable configure -state disabled
	}
	if { $iniIdEmpty == "" } {
	    #set ::KMProps::cmb$whatDisable "Disabled"
	    ::xmlutils::setXml $fullname state "write" disabled
	    
	}
    } else {
	if {[winfo exists $f.fValues.cmb$whatDisable]} {
	    $f.fValues.cmb$whatDisable configure -state normal
	}
	if { $iniIdEmpty == "" } {
	    set dvChild [::xmlutils::setXml $fullname dv]
	    set ::KMProps::cmb$whatDisable $dvChild                
	    ::xmlutils::setXml $fullname state "write" normal
	}
    }
}

proc ::KMProps::getCmbWidth { comboList } {
    
    set width 10
    
    # Validamos el tamaño de los string del combo para ponerle uno o otro tamaño
    foreach c $comboList {
    if { [string length $c] > 5 } {
	    set width 10
	}
	if { [string length $c] > 10 } {
	    set width 15
	}
	if { [string length $c] > 15 } {
	    set width 20
	}
    if { [string length $c] > 20 } {
	    set width 25
	}
    }
    return $width
}   

proc ::KMProps::cmbCancel { item T  } {
    
    set fullname [DecodeName [$T item tag names $item]]
    set idFull [string map { "." "" "//" ""} $fullname]
    
    set pid [::xmlutils::setXml $fullname pid]
    set id [::xmlutils::setXml $fullname id]
    set dv [::xmlutils::setXml $fullname dv]
    set dvText [::xmlutils::setXml $fullname dvText]
    
    set ::KMProps::cmb$idFull "$dv"
    
    set f "$T.f$idFull"
    if {[winfo exists $f]} {
	destroy $f
    }
    
    # Miramos si tiene algun estilo especial
    $T item style set $item C0 styAnyRead
    if { [::xmlutils::setXml $fullname style] == "*" } {
	$T item element configure $item C0 elemTxtRead -text "$pid* : $dvText"
    } else {
	$T item element configure $item C0 elemTxtRead -text "$pid: $dvText"
    }
    
    set ::KMProps::lastSelected {}
}

proc ::KMProps::cmbSelectChange { item T {remove 1} {selectVal "current"} } {
    
    #msg "cmbselectChange $remove"
    set refresh 0
    set fullname [DecodeName [$T item tag names $item]]
    set idFull [string map { "." "" "//" ""} $fullname]
    
    set id [::xmlutils::setXml $fullname id]
    
    set comboState [::xmlutils::getComboBoxState $fullname]
    if { $comboState == "normal" } {
	
	# Si han pulsado el combo en modo editar no tiene que hacer nada (antes desaparecía)
	set f "${T}.f${idFull}.cmb"
	if { [winfo exists $f] && [$f state] == "pressed" } {
	    return
	}
	
	set selCombo [set ::KMProps::cmb$idFull]
	set selComboText $selCombo
	if {$selectVal != "current"} {
	    set ::KMProps::lastSelected {}
	}
    } else {
	
	if {$selectVal == "current"} {
	    
	    set f "$T.f$idFull.cmb"
	    set selCombo [::xmlutils::getComboDv $f $fullname]
	    set selComboText [::xmlutils::getComboDv $f $fullname "text"]
	    set ::KMProps::lastSelected [list $item $selCombo $selComboText]
	    
	} else {  
	    
	    set selCombo [lindex $::KMProps::lastSelected 1]
	    set selComboText [lindex $::KMProps::lastSelected 2]
	    set ::KMProps::lastSelected {}
	}
    }
    
    if {$id == "Ax" || $id == "Ay" || $id == "Az"} { 
	
	set whatDisable [string map {A V} $id]
	set abuelo [$T item parent [$T item parent $item]]
	foreach itemChild [$T item descendants $abuelo] {
	    
	    set fullnameChild [DecodeName [$T item tag names $itemChild]]
	    set idChild [::xmlutils::setXml $fullnameChild id]
	    if { $idChild == $whatDisable } {
		if {$selCombo == 0 } {
		    $T item enabled $itemChild 0
		    $T item element configure $itemChild C0 elemTxtRead -fill { gray }
		    #$T item element configure $item C0 elemTxtRead -text "[$T item text $itemChild 0]: Disabled" -fill { gray }
		    ::xmlutils::setXml $fullnameChild state "write" disabled
		} else {
		    $T item enabled $itemChild 1
		    set dvChild [::xmlutils::setXml $fullnameChild dv]
		    $T item element configure $itemChild C0 elemTxtRead -fill { black }
		    ::xmlutils::setXml $fullnameChild state "write" normal
		}
	    }
	}
    } elseif { $id == "ElemType" } {
	
	global KPriv        
	set matFullname [string map {"ElemType" "MatModel"} $fullname]
	
	::xmlutils::setXml $matFullname dv
	
	##Si se cambia el ElementType hay que cambiar el dv de Material Model
	set xpath "Kratos_KWords/ElementCLaws/Item\[@id='ElementTypes'\]"
	::xmlutils::getAttribute $KPriv(xmlDocKKW) $xpath dv "$selCombo"
	set values [::xmlutils::getXMLValues "$matFullname" "" "" "" "NoElementFilter"]
	
	set dvMat [::xmlutils::setXml $matFullname dv]
	
	if { $dvMat ni $values } {
	    ::xmlutils::setXml $matFullname dv "write" [lindex $values 0]
	    #msg "$$matFullname \n $selCombo --> $values dvMat:$dvMat dv:[::xmlutils::setXml $matFullname dv] values0:[lindex $values 0]"
	    
	}
	#Hacemos un refresh para refrescar los cambios, tanto en matModel com en Thickness
	set refresh 1
	
	#
	#
	#
	#set values [::xmlutils::getXMLValues "$matFullname" "" "" "" "NoElementFilter"]
	#
	#if { $dvMat ni $values } {
	#        ::xmlutils::setXml $matFullname dv [lindex $values 0]
	#        ::KMProps::RefreshTree $T
	#}
    }
    
    set pid [::xmlutils::setXml $fullname pid]
    
    #msg "id:$id rem:$remove sel:$selCombo"
    if { $remove } {
	
	#Eliminamos el frame
	set f "$T.f$idFull"
	if {[winfo exists $f]} {
	    destroy $f
	}
	
	#Miramos si tiene algun estilo especial
	$T item style set $item C0 styAnyRead
	if { [::xmlutils::setXml $fullname style] == "*" } {
	    $T item element configure $item C0 elemTxtRead -text "$pid* : $selComboText"
	} else {
	    $T item element configure $item C0 elemTxtRead -text "$pid: $selComboText"
	}
    }
    
    #Vuelve a poner la imagen por si no recargamos el árbol
    set icon [::xmlutils::setXml $fullname icon]
    set imagen [::WinUtils::GetImage $icon]
    if { $imagen != -1 } {
	$T item image $item C0 $imagen
    }
    
    #Miramos si ha cambiado para actualizar en consecuencia
    if { [::xmlutils::setXml $fullname dv] != $selCombo } {
	
	#Guarda el nuevo valor en el xml
	::xmlutils::setXml $fullname dv "write" $selCombo
	
	#Si el combo tiene un class especial habrá q reconstruir el árbol
	set clase [::xmlutils::setXml $fullname class]
	if { $clase != "" } {
	    
	    #Comprueba si el cambio va a afectar a la estructura del arbol
	    ::KMProps::specialComboAction $T $clase $selCombo $item $id
	    
	    #Volvemos a cargar el árbol de propiedades
	    set refresh 1
	}
    }
    
    if {$refresh} {
	::KMProps::RefreshTree $T
    }
}

#
# Si el atributo class de un nodo requiere algún cambio en el xml
# se prepara todo para q el refresh tree cambie el treeCtl
#
proc ::KMProps::specialComboAction { T clase selCombo item id } {
    
    global KPriv 

    # Variable used to disabked the FSI
    set disabledFSI 1
    # wa "clase:$clase selCombo:$selCombo item:$item id:$id"
    if { $clase == "application" } {
	
	if {$selCombo == "Yes"} {
	    
	    # Este combo solo permite una aplicación activa al mismo tiempo
	    # Así que desactivamos el resto
	    set padre [$T item parent $item]
	    set fullParent [DecodeName [$T item tag names $padre]]
	    set parentId [::xmlutils::setXml $fullParent id]
	    set applications [::xmlutils::setXmlContainerIds $fullParent "Item"]
	    foreach aplicId $applications {
		if { $aplicId != $id } {
		    ::xmlutils::setXml "${fullParent}//i.${aplicId}" dv "write" "No"
		    ::xmlutils::setXml "$aplicId" state "write" "hiddenAll"
		}
	    }
	    
	    # Activamos la que corresponda
	    if { $id == "FluidStructureInteraction" } {
	
		if {$disabledFSI} {
		    ::xmlutils::setXml "FluidStructureInteraction" state "write" "hiddenAll"
		    ::xmlutils::setXml "FluidStructureInteraction" open "write" 0
		} else {
		    # Activamos Fluidos y Estructuras
		    ::xmlutils::setXml "StructuralAnalysis" state "write" "normal"
		    ::xmlutils::setXml "StructuralAnalysis" open "write" 1
		    ::xmlutils::setXml "Fluid" state "write" "normal"
		    ::xmlutils::setXml "Fluid" open "write" 1
		    
		    ::xmlutils::setXml "PFEM" state "write" "hiddenAll"
		    ::xmlutils::setXml "PFEM" open "write" 0
		    
		    # Forzamos el tipo de fluido a "Incompressible"
		    set fluidType [::xmlutils::setXml "Fluid//c.AnalysisData//i.FluidType" dv]
		    if { $fluidType == "Compressible"} {
			::xmlutils::setXml "Fluid//c.AnalysisData//i.FluidType" dv "write" "Incompressible"
		    }
		}
	    } else {
		# Activamos la aplicación seleccionada
		::xmlutils::setXml "$id" state "write" "normal"
		# La desplegamos
		::xmlutils::setXml "$id" open "write" 1

		if {$disabledFSI} {
		    ::xmlutils::setXml "FluidStructureInteraction" state "write" "hiddenAll"
		    ::xmlutils::setXml "FluidStructureInteraction" open "write" 0
		}
	    }
	    
	} else {
	    
	    if { $id == "FluidStructureInteraction" } {
		if {$disabledFSI} {
		    ::xmlutils::setXml "FluidStructureInteraction" state "write" "hiddenAll"
		    ::xmlutils::setXml "FluidStructureInteraction" open "write" 0
		} else {
		    #if{ [::xmlutils::setXml "${fullParent}//i.${aplicId}" dv] } {}
		    ::xmlutils::setXml "StructuralAnalysis" state "write" "hiddenAll"
		    ::xmlutils::setXml "Fluid" state "write" "hiddenAll"
		    ::xmlutils::setXml "PFEM" state "write" "hiddenAll"
		}
	    } else {
		::xmlutils::setXml "$id" state "write" "hiddenAll"
	    }
	    
	}
    }
    
    foreach var $::KMProps::visibilityVars {
	
	if {$var == $clase} {
	    
	    # General case
	    if { [set ::KMProps::$var] != "$selCombo" } {
		set ::KMProps::$var $selCombo
		#msg "var$var"
	    }
	    
	    # Special cases
	    if { $var == "fluidType" } {
		if {$selCombo == "Compressible" } {
		    set txt [= "Compressible fluids are not available in this version"]   
		    WarnWin "$txt."
		    set fullname [DecodeName [$T item tag names $item]]
		    ::xmlutils::setXml "$fullname" dv "write" "Incompressible"
		} 
	    } 
	}
    }
    return ""
}

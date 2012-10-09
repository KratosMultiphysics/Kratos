#####################################################################################
#
#  NAME: kmpropsfwg.tcl
#
#  PURPOSE: Frame and widget used in the kratos main window to manage model/material
#           properties
#
#  QUANTECH ATZ-DEVELOPMENT DEPARTMENT
#
#  AUTHORS : G. Socorro
#
#  CREATED AT: 29/03/2012
#
#  HISTORY:
# 
#   0.8- 09/10/12- G. Socorro, correct a bug when editing the tree in Linux OS => Modify bind $f.cmb <KeyPress> by bind $f.cmb <Return>
#   0.7- 27/09/12-J. Garate, Pick Coordinates Button
#   0.6- 26/04/12-G. Socorro, change GiD_Groups by Cond_Groups
#   0.5- 12/04/2012 JGarate, Ahora selecciona por defecto el primer elemento en Geometry Auto Group
#   0.4- 10/04/2012 JGarate, Cambio de iconos para AutoGroup. Corregidos los problemas en mesh con algunos elementos
#   0.3- 04/04/2012 JGarate, Corregido el error en mallado, por el cual no se permitia usar el Auto New Group, porque no mostraba los botones de seleccion
#   0.2- 02/04/2012 JGarate, labels and combobox grouped left. Combobox width edited. Autogroup button changed
#   0.1- 29/03/2012 G. Socorro, create a base source code from the kmprops.tcl script
#
######################################################################################
#                      Procedures that belong to this file
###############################################################################
#         Name                      |        Functionality
#------------------------------------------------------------------------------
# 1. CreateTreeAndToolbar           | Create the treectrl properties
# 2. CreateBottomFrame              | Create the botton frame
# 3. DestroyBottomFrame             | Destroy the botton frame

proc ::KMProps::CreateTreeAndToolbar { w } {
    # ABSTRACT: Create the treectrl properties 
    # Arguments
    # w => Frame path
    # Return
    # T => The tree path

    set mdf [ttk::frame $w.middle]
    set T [::KMProps::CreateTreeProperties $w]
    
    grid $w.middle -sticky wens
    
    focus $T
    
    return $T
}

proc ::KMProps::CreateBottomFrame { } {
    # ABSTRACT: Create the botton frame
    
    # If exists the bottom frame destroy it
    set f [::KMProps::DestroyBottomFrame]
    
    # Create the frame where set the properties
    ttk::frame $f -borderwidth 0
    
    # Grid for toolbar
    grid $f -row 2 -column 0 -sticky wes
    
    return $f
    
}

proc ::KMProps::DestroyBottomFrame { } {
    # ABSTRACT: Destroy the botton frame
    variable NbPropsPath
    
    set f ${NbPropsPath}.fBottom    
    if {[winfo exists $f]} {
        foreach w [winfo children $f] {
            destroy $w
        }
    	destroy $f
    }
    return $f
}

#
# Construye un frame con un tab por cada container que cuelgue de "item"
# y dentro de cada tab, etiquetas y combos para cada item (si los hay)
#
proc ::KMProps::buildGroupsFrame { T idTemplate item fullname} {
    
    global KPriv
    
    set f [::KMProps::CreateBottomFrame]
    
    # El combo de propiedades es importante resetearlo
    # (posteriormente se utiliza para validar si hay alguna seleccionada)
    if {[info exists ::KMProps::cmbProperty]} {
	unset ::KMProps::cmbProperty
    }
    #Parametro para utilizar o no los combos "Activation"
    set activation 0
    set listT [::KMProps::getTemplateStructure $idTemplate]
    
    if {[llength $listT] >= 1 } {
	set nb ${f}.nb
	grid [ttk::notebook $nb ] -row 0 -column 0 -columnspan 2 -padx 0 -sticky nw -in $f
	
	
	#Lista de listas con formato {idContainer idItem1 idItem2...}
	foreach listContainer $listT {
	    
	    #Si tiene como mínimo el container y un item ponemos tab y dentro label-combos
	    if {[llength $listContainer] >= 2} {
		
            set idContainer [lindex $listContainer 0]
            
            #Si solo hay un container le damos el nombre del item pulsado, no del template
            if {[llength $listT] == 1} {
                set pid "[$T item text $item 0]"
            } else {
                set pid [= [::KMProps::getPropTemplate $idTemplate pid $idContainer]]
                #Lo limitamos por si es demasiado largo
                set pid "[string range $pid 0 20]"
            }
            
            # Para cada container declaramos un tab
            set fTab ${nb}.f$idContainer
            $nb add [ttk::labelframe $fTab -text "[= Properties]" -padding {10 0 10 10}] \
                -text "$pid"
            
            for {set i 1} { $i < [llength $listContainer] } {incr i} {
                
                set id [lindex $listContainer $i]
                
                #Si no coincide la dimensión 2D/3D no lo ponemos
                set nDim [::KMProps::getPropTemplate $idTemplate nDim "$idContainer//$id"]
                if { $nDim == "" || $nDim == $::KMProps::nDim } {
                    
                    #Comprobamos el estado para ver si está activo o no
                    set state [::KMProps::getPropTemplate $idTemplate state "$idContainer//$id"]
                    
                    if {$state != "hidden"} {
                        
                        set pid [::KMProps::getPropTemplate $idTemplate pid "$idContainer//$id"]
                        set tooltip [::KMProps::getPropTemplate $idTemplate tooltip "$idContainer//$id"]
                        set dv [::KMProps::getPropTemplate $idTemplate dv "$idContainer//$id"]
                        set function [::KMProps::getPropTemplate $idTemplate function "$idContainer//$id"]
                        
                        #Obtenemos la lista de valores para el combo si existe
                        set comboList [::xmlutils::getXMLValues "$idContainer//$id" $idTemplate "" $fullname]
                        set icomboList [::xmlutils::getXMLValues "$idContainer//$id" $idTemplate "ivalues" $fullname]
                        
                        set CBState [::KMProps::getPropTemplate $idTemplate CBState "$idContainer//$id"]
                        if { $CBState == "normal" } {
                        set values $comboList
                        set comboList {}
                        } else {
                        set values {}
                        }
                        
                        #Para cada item añadimos label y combo
                        grid [ttk::label $fTab.lbl$id -text "$pid:" ] \
                        -row $i -column 0 -pady 2 -sticky nw -in $fTab
                        
                        if { [llength $comboList] > 0 } {
                            
                            if {$state != "disable"} {
                                set state "readonly"
                            }
                            if {[string length $id] == 2 && ([string index $id end] == "x" || [string index $id end] == "y" || [string index $id end] == "z") } { 
                                set width 15
                            } else {
                                set width 20
                            }
                            grid [ttk::combobox $fTab.cmb$id -values $comboList -state $state -width [::KMProps::getCmbWidth $comboList] -textvariable "::KMProps::cmb$id"] \
                                -row $i -column 1 -padx 3 -pady 2 -sticky nw -in $fTab
                            tooltip::tooltip $fTab.cmb$id [= "%s" $tooltip]
                            
                            if {$id == "Ax" || $id == "Ay" || $id == "Az"} {
                                bind $fTab.cmb$id <<ComboboxSelected>> "::KMProps::cmbDisable $fullname $f.nb $id"
                                set activation 1
                            }
                            ::xmlutils::setComboDv $fTab.cmb$id $fullname $dv $idTemplate 
                        
                        } else {
                        
                            if {$id == "Vx" || $id == "Vy" || $id == "Vz"} {
                                if { $activation } {
                                set activeId "A[string range $id 1 1]"
                                set value [set "::KMProps::cmb$activeId"]
                                if { $value == 0 } {
                                    set state "disabled"
                                }
                                }
                            }
                            
                            grid [ttk::combobox $fTab.cmb$id -state $state -values $values -textvariable "::KMProps::cmb$id" -width [::KMProps::getCmbWidth $comboList]] \
                                -row $i -column 1 -padx 3 -pady 2 -sticky nw -in $fTab
                            tooltip::tooltip $fTab.cmb$id [= "%s" $tooltip]
                            
                            #set dv [::KMProps::getPropTemplate $idTemplate dvText "$idContainer//$id"]
                            set ::KMProps::cmb$id $dv
                        
                        }
                        
                
                        if {$function != "" } {
                            grid [ttk::button $fTab.funct$id -text "functions" -command "KFun::InitBaseWindow $fTab $id" -style TMenubutton.Toolbutton -width [::KMProps::getCmbWidth $comboList]]  \
                                -row $i -column 2 -sticky nw  -pady 0 -padx 3 -in $fTab
                            tooltip::tooltip $fTab.funct$id [= "Function manager"]
                            set img [::WinUtils::GetImage "functions.gif"]
                            if { $img != -1 } { $fTab.funct$id configure -image $img }
                            
                            grid [ttk::button $fTab.deleteFunct$id -text "delete" -command "::KMProps::unassignFunction $fTab $id $fullname" -style TMenubutton.Toolbutton -width [::KMProps::getCmbWidth $comboList]] \
                                -row $i -column 3 -sticky nw  -pady 0 -padx 2 -in $fTab
                            tooltip::tooltip $fTab.deleteFunct$id [= "Unassign function"]
                            set img [::WinUtils::GetImage "delete_icon.gif"]
                            if { $img != -1 } { $fTab.deleteFunct$id configure -image $img }
                        
                        }
                    }
                }
            }
	    }
	}
	
	#
	# PUNTOS, LINEAS, SUPERFICIES Y VOLUMENES
	#
	set whatuse [GiD_Info Project ViewMode]
	set geomlist [GetGeometryImageFiles]
	set meshlist {node.gif element_line.gif element_surface.gif element_volume.gif}
	set col 0
	
	#Miramos si hay restricción de entidades
	set entityList [split [::xmlutils::setXml $fullname GiDEntity] "," ]
	
	#set entityList [split [::KMProps::getPropTemplate $idTemplate GiDEntity] "," ]
	
	switch $whatuse {
	    
	    GEOMETRYUSE {
		foreach i $geomlist {
		    #::KMProps::CorrectFileExtensionTo ".png" $i
		    
		    set command [file rootname $i]
		    if {$command in $entityList || $entityList == "" } {
			
			
			set fb "${f}.b$command"
			grid [ttk::button $fb -text "$i" -command "::KMProps::changeImage $command $i $f" -style TMenubutton.Toolbutton -width 15] \
			    -row 1 -column 0 -sticky nw  -pady 3 -padx [expr (50 * $col) + 15] -in $f
			tooltip::tooltip $fb [= "Entity %s" $command]
			set im {}
			append im $command ".gif"
			set img [::WinUtils::GetImage $im]
			if { $img != -1 } {
			    $fb configure -image $img
			}
			if {$col == 0} {
			    #Por defecto dejamos marcada la primera
			    set i "[string range $i 0 [expr [string length $i] - 5]]_sel.gif"
			    set ::KMProps::selectedEntity $command
			    ::KMProps::changeImage $::KMProps::selectedEntity $i $f
			}
			incr col
		    }
		}
	    }
	    MESHUSE {
		    set meshentityList ""
		    if { "point" in $entityList } {
			lappend meshentityList "node"
		    } 
		    if { "line" in $entityList } {
			lappend meshentityList "element_line"
		    }
		    if { "surface" in $entityList } {
			lappend meshentityList "element_surface"
		    }
		    if { "volume" in $entityList } {
			lappend meshentityList "element_volume"
		    }
		    
		foreach i $meshlist {
		    set command [file rootname $i]
		    if {$meshentityList == "" || $command in $meshentityList} {
			if {$col == 0} {
			    #Por defecto dejamos marcada la primera
			    set i "[string range $i 0 [expr [string length $i] - 5]]_sel.gif"
			    set ::KMProps::selectedEntity $command
			}
			
			set fb "${f}.b$command"
			grid [ttk::button $fb -text "$i" -command "::KMProps::changeImage $command $i $f" -style TMenubutton.Toolbutton -width 15] \
			    -row 1 -column 0 -sticky nw  -pady 3 -padx [expr (50 * $col) + 15] -in $f
			if { [string range $command 0 6] == "element" } {
				set command [string range $command 8 end]
				#append i "element_" $i
			}
			
			tooltip::tooltip $fb [= "Entity $command"]
			set img [::WinUtils::GetImage $i]
			
			if { $img != -1 } {
			    $fb configure -image $img
			}
			incr col
		    }
		}
	    }
	}
	
	grid [ttk::label $f.lGroups -text "Group:" ] \
	    -row 2 -column 0 -pady 3 -padx 5 -sticky nw -in $f
	
	#COMBO DE GRUPOS
	set fGroups $f.cGroups
    
	set filterGroups [Cond_Groups list]
	grid [ttk::combobox $fGroups -state readonly -values "$filterGroups" -textvariable "::KMProps::selGroup" \
		-postcommand "::KMProps::changeGroups [list $entityList] $fGroups" -width 15] \
		-row 2 -column 0 -pady 3 -padx 55 -sticky nw -in $f
	
	set ::KMProps::selGroup ""
	if { [llength $filterGroups] > 0 } {
	    $f.cGroups current 0
	}
	bind $f.cGroups <<ComboboxSelected>> "::KMProps::cmbChangeCheckGroups $f"
	
	# BOTON A LA DERECHA DE LOS GRUPOS (CREAR GRUPO AUTOMATICAMENTE)
	grid [ttk::button $f.iGroups -text [= "newGroup"] -command "::KMProps::autoNewGroup $id" ] \
	    -row 2 -column 0 -sticky nw  -pady 0 -padx 180 -in $f
	tooltip::tooltip $f.iGroups [= "Create automatic new group"]
	$f.iGroups configure -image [::WinUtils::GetImage "newAutoGroup.gif" ]
	
	grid [ttk::button $f.bPropOk -text [= "Ok"]  -command "::KMProps::acceptGroups $T $idTemplate $fullname $item {$listT} {$entityList} $fGroups" ] \
	    -row 3 -column 0 -sticky nw  -pady 5 -padx 5 -in $f
	tooltip::tooltip $f.bPropOk [= "Assign condition to the selected group"]
	
	grid [ttk::button $f.bPropCancel -text [= "Cancel"]  -command "::KMProps::DestroyBottomFrame" ] \
	    -row 3 -column 0 -sticky nw  -pady 5 -padx 130  -in $f
	tooltip::tooltip $f.bPropCancel [= "Cancel assignation"]
	
    }
    
    bind $T <KeyPress> "if { %k == 27   } { ::KMProps::DestroyBottomFrame }"
}

proc ::KMProps::unassignFunction { f id fullname } {
    
    set fcmb "$f.cmb$id"
    
    if {[winfo exists $fcmb]} {
	$fcmb configure -state normal
    }
    
    set ::KMProps::cmb$id "0.0"    
}

proc ::KMProps::assignFunction { f id idFunction } {
    
    set fcmb "$f.cmb$id"
    
    if {[winfo exists $fcmb]} {
	$fcmb configure -state disabled
    }
    set ::KMProps::cmb$id "$idFunction"
}

#
#Consulta el xml KKWord para ver si se tiene que mostrar el elemento del Tickness
#
proc ::KMProps::showThickness { } {
    
    #return 1
    
    global KPriv
    
    #Leemos el xml si aun no se ha leído
    #::xmlutils::initKKWord
    
    if { [info exists ::KMProps::ElemTypeThickness] } {
	set dv $::KMProps::ElemTypeThickness
	unset ::KMProps::ElemTypeThickness
    } else {
	set xpath "Kratos_KWords/ElementCLaws/Item\[@id='ElementTypes'\]"
	set dv [::xmlutils::getAttribute $KPriv(xmlDocKKW) $xpath dv]
    }
    
    set xpath "Kratos_KWords/ElementCLaws/Item\[@id='Tickness$::KMProps::nDim'\]"
    set ListElementType [split [::xmlutils::getAttribute $KPriv(xmlDocKKW) $xpath elementType] ","]
    
    if { $dv in $ListElementType } {
	return 1
    } else {
	return 0
    }
}


#
# Construye un frame consultando el template de propiedades correspondiente,
#  con un tab por cada container y dentro de cada tab, 
#  etiquetas y combos para cada item (si los hay)
#
proc ::KMProps::buildPropertyFrame { T idTemplate item fullname } { 
    
    global KPriv
    # WarnWinText "T:$T idTemplate:$idTemplate item:$item fullname:$fullname"
    set f [::KMProps::CreateBottomFrame]
    
    set listT [::KMProps::getTemplateStructure $idTemplate]
    # WarnWinText "listT:$listT"
    if {[llength $listT] >= 1 } {
        set nb ${f}.nb
        grid [ttk::notebook $nb ] -row 0 -column 0 -columnspan 2 -padx 0 -sticky nw -in $f
        
        #Lista de listas con formato {idContainer idItem1 idItem2...}
        foreach listContainer $listT {
            
            #Si tiene como mínimo el container(1er elemento) y un item ponemos tab y dentro label-combos
            if {[llength $listContainer] > 1} {
            
                set idContainer [lindex $listContainer 0]
                
                set pid [::KMProps::getPropTemplate $idTemplate pid $idContainer]
                # wa "idContainer:$idContainer pid:$pid cero:[lindex $listT 0 0]"
                # Para cada container declaramos un tab
                set fTab ${nb}.f$idContainer
                set ptxt "[= Properties]"
                $nb add [ttk::labelframe $fTab -text "$ptxt" -padding {10 0 10 10}] \
                    -text "[string range $pid 0 20]"
                
                #En el caso del primer tab forzamos el item de "Nombre de propiedad" y 2campos mas
                if { $idContainer == [lindex $listT 0 0] } {
                    
                    # Property identifier
                    if {$idTemplate eq "CutProperties"} {
                    set cptxt [= "Cut Id"]
                    set cbhelp [= "Choose a new cut identifier"]
                    set whatoption "Cut"
                    } elseif {$idTemplate eq "HistoryOutputOnPoints"} {
                    set cptxt [= "Graph Id"]
                    set cbhelp [= "Choose a new graph identifier"]
                    set whatoption "Graph"
                    } else {
                    set cptxt [= "Property Id"]
                    set cbhelp [= "Choose a new property identifier"]
                    set whatoption "Property"
                    }

                    grid [ttk::label $fTab.lblName -text "$cptxt:" ] \
                    -row 0 -column 0 -pady 5 -sticky nw -in $fTab
                    
                    grid [ttk::combobox $fTab.cmbPropertyName -state normal -textvariable "::KMProps::propertyName" -width 10	 ] \
                    -row 0 -column 1 -padx 3 -pady 5 -sticky nw -in $fTab

                    tooltip::tooltip $fTab.cmbPropertyName $cbhelp
                    
                    set ::KMProps::propertyName "[::KEGroups::GetAutomaticPropertyName $fullname "$whatoption"]"
                    focus $fTab.cmbPropertyName
                }
                for {set i 1} { $i < [llength $listContainer] } {incr i} {
                    
                    set id [lindex $listContainer $i]
                    
                    #Los nodos ocultos no se deben mostrar
                    set state [::KMProps::getPropTemplate $idTemplate state "$idContainer//$id"]
                    if {$state != "hidden" } {
                    
                        #Si no coincide la dimensión 2D/3D no lo ponemos
                        set nDim [::KMProps::getPropTemplate $idTemplate nDim "$idContainer//$id"]
                        if { $nDim == "" || $nDim == $::KMProps::nDim } {
                            
                            set pid [::KMProps::getPropTemplate $idTemplate pid "$idContainer//$id"]
                            set tooltip [::KMProps::getPropTemplate $idTemplate tooltip "$idContainer//$id"]
                            set dv [::KMProps::getPropTemplate $idTemplate dv "$idContainer//$id"]
                            
                            set CBState [::KMProps::getPropTemplate $idTemplate CBState "$idContainer//$id"]
                            if { $CBState == "normal" } {
                                set values $comboList
                                set comboList {}
                            } else {
                                set values {}
                            }
                            
                            # Es importante que ElemType sea el primer combo para actualizar el "dv" en el xml
                            if { $id == "ElemType" } {
                                set comboList [::xmlutils::getXMLValues "$idContainer//$id" $idTemplate "" "" $dv]
                            } else {
                                set comboList [::xmlutils::getXMLValues "$idContainer//$id" $idTemplate "" "$fullname"]
                            }
                            
                            #Para cada item añadimos label y combo
                            grid [ttk::label $fTab.lbl$id -text "$pid:" ] \
                            -row $i -column 0 -pady 5 -sticky nw -in $fTab
                            
                            
                            if { [llength $comboList] } {
                            
                            grid [ttk::combobox $fTab.cmb$id -values $comboList -state readonly -width [::KMProps::getCmbWidth $comboList] \
                                  -textvariable "::KMProps::cmb$id" \
                                  -postcommand [list ::KMProps::changeCmbValues "$fTab.cmb$id" "$idContainer//$id" "$idTemplate" "" "$fullname"] ] \
                                -row $i -column 1 -padx 5 -pady 2 -sticky nw -in $fTab
                            
                            ::xmlutils::setComboDv $fTab.cmb$id $fullname $dv $idTemplate
                            tooltip::tooltip $fTab.cmb$id $tooltip
                            
                            set psb [::KMProps::getPropTemplate $idTemplate sbi "$idContainer//$id"]
                            
                            
                            #En este caso se tendrá qué recargar si existe el combo de "Material Model"
                            if { $id == "ElemType" } {
                                bind $fTab.cmb$id <<ComboboxSelected>> "::KMProps::cmbElemTypeChange $fTab.cmb$id $idContainer//$id $idTemplate"
                            }
                            
                            } else {
                                
                                grid [ttk::combobox $fTab.cmb$id -state normal -width [::KMProps::getCmbWidth $comboList] -values $values -textvariable "::KMProps::cmb$id"] \
                                    -row $i -column 1 -padx 5 -pady 2 -sticky nw -in $fTab
                                set ::KMProps::cmb$id $dv
                                tooltip::tooltip $fTab.cmb$id $tooltip
                                
                                # Pick Coordinates Button
                                set psb [::KMProps::getPropTemplate $idTemplate sbi "$idContainer//$id"]
                                if { $psb == "PickCoordinates" } {
                                    set sbxp [::KMProps::getPropTemplate $idTemplate sbxp "$idContainer//$id"]
                                    set sbp [::KMProps::getPropTemplate $idTemplate sbp "$idContainer//$id"]
                                    grid [ttk::button $fTab.btn$id -text $sbp -command [list ::KMProps::selectionButton $sbxp] -width 4 ] \
                                    -row $i -column 2 -padx 5 -pady 2 -sticky ew -in $fTab
                                    tooltip::tooltip $fTab.btn$id "Pick Coordinates"
                                }
                            }
                            
                            
                            #En el caso de ElemType actualizamos la variable de filtrado del Thickness
                            if { $id == "ElemType" } {
                                set ::KMProps::ElemTypeThickness $dv
                                
                            #En el caso del Thickness, no siempre se mostrará
                            } elseif { $id == "Thickness" && ![::KMProps::showThickness]} {
                            
                                grid remove $fTab.lbl$id
                                grid remove $fTab.cmb$id
                            
                            }
                        }
                    }
                }
            }                                                
        }
    }
    
    grid [ttk::button $f.bPropOk -text [= "Ok"]  -command "::KMProps::acceptProperty $T $idTemplate $fullname $item {$listT}" ] \
	-row 3 -column 0 -sticky nw  -pady 3 -padx 20  -in $f
    
    tooltip::tooltip $f.bPropOk [= "Assign the defined properties"]
    
    grid [ttk::button $f.bPropCancel -text [= "Cancel"]  -command "::KMProps::DestroyBottomFrame" ] \
	-row 3 -column 0 -sticky nw  -pady 3 -padx 100  -in $f
    tooltip::tooltip $f.bPropCancel [= "Cancel the defined properties"]
    
    bind $T <KeyPress> "if { %k == 27   } { ::KMProps::DestroyBottomFrame }"
    
}

proc ::KMProps::selectionButton { psb } {
    set xyz [GidUtils::GetCoordinates]
    set psb [split $psb ","]
    
    foreach var $psb val $xyz {
        set ::KMProps::cmb$var $val
    }
    
}

proc ::KEGroups::GetAutomaticPropertyName { fullname {baseid "Property"}} {
    
    set name ""
    if {$baseid eq "Property"} {
	set fpath "[::KMProps::getApplication $fullname]//c.Properties"
    } elseif {$baseid eq "Cut"} {
	set fpath "$fullname"
    } elseif {$baseid eq "Graph"} {
	set fpath "$fullname"
    }
    # wa "fpath:$fpath"
    set propIds [::xmlutils::setXmlContainerIds $fpath]
    set i 0
    # wa "baseid:$baseid fullname:$fullname propIds:$propIds"
    if {[llength propIds]} {
	for {set i 1} {$i<=$::KMProps::MaxIdIter} {incr i} {
	    set name "${baseid}${i}"
	    if { [lsearch -exact $propIds $name] == -1 } { break }
	}
    } else {
	set name "${baseid}"
    }
    return $name
}

#
# Click en botón OK de Property
#
proc ::KMProps::acceptProperty { T idTemplate fullname item listT} {
    
    set property $::KMProps::propertyName
    
    #Validamos que la propiedad no tenga carácteres extraños
    set property [::KUtils::parseTreeStr $property]
    if { $property == -1 } {
	
	WarnWin [= "You can't use some reservate chars like:\n  :   /   $   .   \\  %  "]
	set ::KMProps::propertyName ""
	return ""
    }
    
    #Comprobamos q el nombre no sea vacío
    if { $property == "" } {
	WarnWin [= "The property name can not be empty"]
	return ""
    }
    
    #Comprobamos q la propiedad aun no exista
    set id [::xmlutils::setXml "$fullname//c.$property" id]

    if { $id != "" } {
	
	WarnWin [= "This name property it is already assigned."]
	return
    }

    ::KMProps::copyTemplate ${idTemplate} $fullname "$property" "Property"
    
    #
    # Ahora debemos actualizar todos los valores en el xml
    #
    if {[llength $listT] >= 1 } {
	
	#Lista de listas con formato {idContainer idItem1 idItem2...}
	foreach listContainer $listT {
	    
	    #Si tiene como mínimo el container y un item entramos
	    if {[llength $listContainer] >= 2} {
		
		set idContainer [lindex $listContainer 0]
		
		#Recorremos los items
		for {set i 1} { $i < [llength $listContainer] } {incr i} {
		    
		    set id [lindex $listContainer $i]
		    
		    #if { $id != "Thickness" || ($id == "Thickness" && [::KMProps::showThickness])} {
		    
		    set fullNombre "$fullname//c.$property//c.$idContainer//i.$id"
		    
		    #Los nodos ocultos no existían en el formulario
		    set state [::xmlutils::setXml $fullNombre state]
		    if {$state != "hidden" } {
			
			#Validamos la dimensión de cada elemento
			set nDim [::xmlutils::setXml $fullNombre nDim]
			if { $nDim == "" || $nDim == $::KMProps::nDim } {
			    
			    set value [set ::KMProps::cmb$id]
			    
			    ::xmlutils::setXml $fullNombre dv "write" $value
			}
		    }
		    #}
		}
	    }
	}
    }

    #Destruimos el frame inferior
    ::KMProps::DestroyBottomFrame
    
    #Recargar tree
    ::KMProps::RefreshTree $T
}

proc ::KMProps::changeImage {entity img path} {
    
    variable selectedEntity
    
    if {$selectedEntity != "" } {
	
	set f "$path.b$selectedEntity"
	
	$f configure -image [::WinUtils::GetImage "${selectedEntity}.gif"]
    }
    
    "$path.b$entity" configure -image [::WinUtils::GetImage "${entity}_sel.gif"]

    set selectedEntity $entity
}

#
# Construye un frame con un tab por cada container que cuelgue de "item"
# y dentro de cada tab, etiquetas y combos para cada item (si los hay)
# 
proc ::KMProps::buildTabFrame { T item {class "Tab"} } {
    
    set f [::KMProps::CreateBottomFrame]

    # Miramos los descendientes directos y si son container ponemos un tab por cada uno q tenga items
    set children [$T item children $item]

    set listTabs {}
    set listItems {}
    set acceptItems {}
    foreach itemChild $children {
	
	set fullname [DecodeName [$T item tag names $itemChild]]
	
	set nodeName [::xmlutils::getXmlNodeName $fullname]
	
	#Miramos si cada hijo es container o item
	if { $nodeName == "Container" } {
	    #Si no tiene items no agregamos el tab
	    if { [$T item numchildren $itemChild] > 0 } {
		lappend listTabs $fullname
	    }
	    #msg "numchild:[$T item numchildren $itemChild]"
	} elseif { $nodeName == "Item" } {
	    lappend listItems $itemChild
	}
    }
    
    #Reseteamos la variable que nos indica si estamos en una propiedad
    set ::KMProps::propertyName ""
    
    # Si no tiene containers pero tiene items, utilizamos como tab el elemento padre seleccionado
    if { [llength $listTabs] == 0 && [llength $listItems] >= 1 } {
	
	set fullname [DecodeName [$T item tag names $item]]
	lappend listTabs $fullname
    } else {
	set listItems {}
    }
    
    if { [llength $listTabs] >= 1 } {
	
	set nb ${f}.nb
	grid [ttk::notebook $nb ] -row 0 -sticky ewn
	
	set i 0
	foreach fullname $listTabs {
	    
	    set id [::xmlutils::setXml $fullname id]
	    set pid [::xmlutils::setXml $fullname pid]                                
	    set help [::xmlutils::setXml $fullname help]
	    
	    #Para cada container declaramos un tab
	    set fTab ${nb}.f$id
	    
	    $nb add [ttk::labelframe $fTab -text "[= Properties]" -padding {10 10 10 10}] -text "[string range [= $pid] 0 20]"
	    
	    if {[llength $listItems] > 0 } {
		#En este caso en realidad los nietos son los propios hijos del item pulsado
		set nietos [$T item children $item]
	    } else {
		#Miramos los hijos de cada container (nietos del pulsado originalmente)
		set itemChild [lindex $children $i]
		set nietos [$T item children $itemChild]
	    }
	    
	    ######################
	    # NOMBRE PROPIEDAD - Añadimos "Nombre de propiedad"
	    if {$class == "Property" && $i == 0 } {
		
		#Nombre de propiedad
		grid [ttk::label $fTab.lblName -text "[= Property Name:]" ] \
		    -row 0 -column 0 -pady 5 -sticky nw -in $fTab
		
		grid [ttk::combobox $fTab.cmbPropertyName -state normal -textvariable "::KMProps::propertyName" -width 15] \
		    -row 0 -column 1 -padx 3 -pady 5 -sticky nw -in $fTab
		
		set ::KMProps::propertyName "[$T item text $item 0]"
		
	    } elseif { $class == "Group" && $i == 0 } {
		
		######################
		#COMBO DE GRUPOS
		set fullname [DecodeName [$T item tag names [$T item parent $item]]]
		
		set ::KMProps::selGroup "[$T item text $item 0]"
		
		set entityList [split [::xmlutils::setXml $fullname GiDEntity] ","]
		#set filterGroups [::KMProps::getGroups $entityList $fullname]
		set filterGroups [Cond_Groups list]
		#msg "lleggo a  filterGroups = $filterGroups"
		grid [ttk::label $fTab.lblName -text "[= Group:]" ] \
		    -row 0 -column 0 -pady 5 -sticky nw -in $fTab
		
		set fGroups $fTab.cGroups
		grid [ttk::combobox $fGroups -state readonly -values "$filterGroups" -textvariable "::KMProps::selGroup"  -postcommand "::KMProps::changeGroups [list $entityList] $fGroups $fullname" -width 15] \
		    -row 0 -column 1 -pady 5 -sticky nw -in $fTab
		
		
		
		bind $fTab.cGroups <<ComboboxSelected>> "::KMProps::cmbChangeCheckGroups $fTab"
		######################
	    }
	    
	    set row 1
	    foreach nieto $nietos {
		
		set fullname [DecodeName [$T item tag names $nieto]]
		set dv [::xmlutils::setXml $fullname dv]
		#Comprobamos q sea un item
		if { [::xmlutils::getXmlNodeName $fullname]  == "Item" } {
		    
		    set state [::xmlutils::setXml $fullname state]
		    
		    if {$state != "hidden"} {
			
			set id [::xmlutils::setXml $fullname id]
			set pid [::xmlutils::setXml $fullname pid]
			set tooltip [::xmlutils::setXml $fullname help]
			set function [::xmlutils::setXml $fullname function]
			
			#Nos guardamos cada item para actualizar sus valores al final
			lappend acceptItems $nieto $fTab

			set comboList [::xmlutils::getXMLValues "$fullname"]
			
			set CBState [::xmlutils::setXml $fullname CBState]
			if { $CBState == "normal" } {
			    set values $comboList
			    set comboList {}
			} else {
			    set values {}
			}
			
			#Para cada item añadimos label y combo
			grid [ttk::label $fTab.lbl$id -text "${pid}:" ] \
			    -row $row -column 0 -padx 3 -pady 5 -sticky nw -in $fTab
			
			if { [llength $comboList] > 0 } {
			    
			    if {$state != "disabled"} {
				set state "readonly"
			    }
			    grid [ttk::combobox $fTab.cmb$id -values $comboList -state $state \
				      -textvariable "::KMProps::cmb$id" -width [::KMProps::getCmbWidth $comboList]]\
				-row $row -column 1 -padx 3 -pady 5 -sticky nw -in $fTab 
			    
			    ::xmlutils::setComboDv $fTab.cmb$id $fullname $dv
			    #set selected [::xmlutils::getSelected $dv $comboList]
			    #$fTab.cmb$id current $selected
			    
			    if {$id == "Ax" || $id == "Ay" || $id == "Az"} {
				::KMProps::cmbDisable $fullname $f.nb
				bind $fTab.cmb$id <<ComboboxSelected>> "::KMProps::cmbDisable $fullname $f.nb"
			    } elseif { $id == "ElemType" } {
				#En este caso se tendrá qué recargar si existe el combo de "Material Model"
				bind $fTab.cmb$id <<ComboboxSelected>> "::KMProps::cmbElemTypeChange $fTab.cmb$id $fullname"
			    }
			    
			} else {
			    if {$state != "disabled"} {
				set state "normal"
			    }
			    
			    grid [ttk::combobox $fTab.cmb$id -state $state -values $values -textvariable "::KMProps::cmb$id" -width [::KMProps::getCmbWidth $comboList]] \
				-row $row -column 1 -padx 3 -pady 5 -sticky nwe -in $fTab 
			    set ::KMProps::cmb$id $dv
			}
			
			if {$function != "" } {
			    grid [ttk::button $fTab.funct$id -text "functions" -command "KFun::InitBaseWindow $fTab $id" -style TMenubutton.Toolbutton -width [::KMProps::getCmbWidth $comboList]] \
				-row $row -column 2 -sticky nw  -pady 0 -padx 3 -in $fTab
			    tooltip::tooltip $fTab.funct$id [= "Function manager"]
			    set img [::WinUtils::GetImage "functions.gif"]
			    if { $img != -1 } { $fTab.funct$id configure -image $img }
			    
			    grid [ttk::button $fTab.deleteFunct$id -text "delete" -command "::KMProps::unassignFunction $fTab $id $fullname" -style TMenubutton.Toolbutton -width [::KMProps::getCmbWidth $comboList]] \
				-row $row -column 3 -sticky nw  -pady 0 -padx 2 -in $fTab
			    tooltip::tooltip $fTab.deleteFunct$id [= "Unassign function"]
			    set img [::WinUtils::GetImage "delete_icon.gif"]
			    if { $img != -1 } { $fTab.deleteFunct$id configure -image $img }
			    
			}
			
			#En el caso de ElemType actualizamos la variable de filtrado del Thickness
			if { $id == "ElemType" } {
			    set ::KMProps::ElemTypeThickness $dv
			    
			    #En el caso del Thickness, no siempre se mostrará
			} elseif { $id == "Thickness" && ![::KMProps::showThickness]} {
			    
			    grid remove $fTab.lbl$id
			    grid remove $fTab.cmb$id
			}
			
			incr row
		    }
		}
	    }
	    incr i
	}
	# Si pulsan Esc también forzamos la salida del Tab
	bind $T <KeyPress> "if { %k == 27   } { ::KMProps::DestroyBottomFrame }"
	#bind $T <KeyPress> "if { %k == 13   } {  [list ::KMProps::acceptTabFrame $T $acceptItems $class $item] }"
	#bind $::KMProps::WinPath <KeyPress> "if { %k == 13   } { [list ::KMProps::acceptTabFrame $T $acceptItems $class $item] }"
	
	if { [llength $acceptItems] } {
	    grid [ttk::button $f.bPropOk -text "Ok"  -command "[list ::KMProps::acceptTabFrame $T $acceptItems $class $item]" ] \
		-row 1 -column 0 -sticky sw  -pady 3 -padx 20  -in $f
	    tooltip::tooltip $f.bPropOk [= "Confirm values"]
	    
	    grid [ttk::button $f.bPropCancel -text "Cancel"  -command "::KMProps::DestroyBottomFrame" ] \
		-row 1 -column 0 -sticky sw  -pady 3 -padx 100  -in $f
	    tooltip::tooltip $f.bPropCancel [= "Cancel assignation"]
	}
    }
}

#
# Contruye un combo dinámico en el item pulsado del árbol y posteriormente se elimina
#
proc ::KMProps::buildFrame { T item { type "props" } } {
    
    set fullname [DecodeName [$T item tag names $item]]
    set idFull [string map { "." "" "//" ""} $fullname]
    
    #Comprobamos que sea un item
    if { [::xmlutils::getXmlNodeName $fullname]  == "Item" } {

	set id [::xmlutils::setXml $fullname id]
	set pid [::xmlutils::setXml $fullname pid]
	set dv [::xmlutils::setXml $fullname dv]
	set tooltip [::xmlutils::setXml $fullname help]
	set comboList [::xmlutils::getXMLValues $fullname]
	
	set CBState [::xmlutils::setXml $fullname CBState]
	if { $CBState == "normal" } {
	    set values $comboList
	    set comboList {}
	} else {
	    set values {}
	}
	
	#---------------------------#---------------------------#
	# Configurar frame en función del XML
	#---------------------------#---------------------------#
	set bg "#F8F8F8"
	set f "$T.f$idFull"
	
	if { [winfo exists $f] } {
	    return
	}
	set f [frame "$T.f$idFull" -borderwidth 0 -background $bg]
	
	if { [llength $comboList] > 0 } {
	    
	    grid [ttk::combobox $f.cmb -values $comboList -state readonly -width [::KMProps::getCmbWidth $comboList] -textvariable "::KMProps::cmb$idFull"] \
		-row 0 -column 0 -padx 3 -sticky nw -in $f 
	    
	    ::xmlutils::setComboDv $f.cmb $fullname $dv
	    #set selected [::xmlutils::getSelected $dv $comboList]
	    #$f.cmb current $selected
	    
	    bind $f.cmb <<ComboboxSelected>> "::KMProps::cmbSelectChange $item $T 1 current"
	} else {
	    
	    grid [ttk::combobox $f.cmb -state normal -values $values -textvariable "::KMProps::cmb$idFull" -width [::KMProps::getCmbWidth $comboList]] \
		-row 0 -column 0 -padx 3 -sticky nw -in $f
	    
	    set ::KMProps::cmb$idFull $dv
	    
	    bind $f.cmb <Leave> "::KMProps::cmbSelectChange $item $T 0 current"
	    bind $f.cmb <FocusOut> "::KMProps::cmbSelectChange $item $T 1 current"
	    bind $f.cmb <Escape> "::KMProps::cmbCancel $item $T"
	}
	# Si pulsan intro o Esc también forzamos la salida del combo
	bind $f.cmb <Return> [list ::KMProps::cmbSelectChange $item $T 1 current]
	bind $T <KeyPress> "if { %k == 27   } { ::KMProps::cmbCancel $item $T }"
	
	tooltip::tooltip $f.cmb "$tooltip"
	
	return $f
	
    } else {
	
	return ""
    }
}

proc ::KMProps::acceptTabFrame { T listItems class {itemSel ""}} {
    
    #
    # Actualiza el resto de valores del xml con el contenido de cada combo
    #
    foreach {item f} $listItems {
	
	set fullname [DecodeName [$T item tag names $item]]
	set idFull [string map { "." "" "//" ""} $fullname]
	set id [::xmlutils::setXml $fullname id]
	set pid [::xmlutils::setXml $fullname pid]
	set estado [::xmlutils::setXml $fullname state]
	
	
	#msg "id: $id estado:[$f.cmb$id cget -state]"
	
	set cmbState [$f.cmb$id cget -state]
	
	if { $cmbState == "readonly" } {
	    
	    set selCombo [::xmlutils::getComboDv "$f.cmb$id" $fullname]
	    set selComboText [::xmlutils::getComboDv "$f.cmb$id" $fullname "text"]
	} else {
	    
	    set selCombo [set ::KMProps::cmb$id]
	    set selComboText $selCombo 
	    
	}
	#$cmbState == "readonly"
	
	
	set function [::xmlutils::setXml $fullname function]
	if { $function != "" } {
	    if { $cmbState == "disabled" } {
		::xmlutils::setXml $fullname state "write" "disabled"
		::xmlutils::setXml $fullname function "write" 1
	    } else {
		::xmlutils::setXml $fullname state "write" "normal"
	    }
	}
	
	#Guarda el nuevo valor en el xml
	::xmlutils::setXml $fullname dv "write" $selCombo
	$T item style set $item C0 styAnyRead
	
	#Si el combo era especial se tendrá q reconstruir el árbol
	set clase [::xmlutils::setXml $fullname class]
	
	#Comprueba si el cambio va a afectar a la estructura del arbol
	::KMProps::specialComboAction $T $clase $selCombo $item $id
    }
    
    #
    # Caso especial para las propiedades
    #
    if {$class == "Property" && $itemSel != "" } {
	
	set newPropertyName $::KMProps::propertyName
	
	#Validamos que la propiedad no tenga carácteres extraños
	set newPropertyName [::KUtils::parseTreeStr $newPropertyName]
	if { $newPropertyName == -1 } {
	    
	    WarnWin [= "You can't use some reservate chars like:\n  :   /   $   .   \\  %  "]
	    set ::KMProps::propertyName ""
	    return
	}
	
	#Comprobamos q el nombre no sea vacío
	if { $newPropertyName == "" } {
	    WarnWin [= "The property name can not be empty"]
	    return
	}
	
	#Comprobamos q si la propiedad ha cambiado no exista el nuevo nombre
	set fullname [DecodeName [$T item tag names $itemSel]]
	set oldId [::xmlutils::setXml "$fullname" id]
	set repeatId [::xmlutils::setXml "[::KMProps::getApplication $fullname]//c.Properties//c.$newPropertyName" id]
	
	if { $repeatId != "" && $oldId != $newPropertyName } {
	    WarnWin [= "This property name already exist."]
	    return
	}
	
	#Si todo es correcto, renombramos la propiedad en todo el árbol y luego lo recargamos
	#(deleteProps hace un rename si llega el 3er argumento "newPropertyName")
	::KMProps::deleteProps $T $itemSel $newPropertyName
	
	::xmlutils::setXml $fullname pid "write" $newPropertyName
	::xmlutils::setXml $fullname id "write" $newPropertyName
	
    } elseif { $class == "Group" } {

	#
	# Caso especial para los Grupos
	#
	set fullname [DecodeName [$T item tag names $itemSel]]
	
	::xmlutils::setXml $fullname pid "write" $::KMProps::selGroup
	::xmlutils::setXml $fullname id "write" $::KMProps::selGroup
    }
    
    #Volvemos a cargar el árbol para q el path de los items sea correcto
    ::KMProps::RefreshTree $T
    
    ::KMProps::DestroyBottomFrame
}

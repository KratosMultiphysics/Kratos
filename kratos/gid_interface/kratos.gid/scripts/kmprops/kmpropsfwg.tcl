#####################################################################################
#
#  NAME: kmpropsfwg.tcl
#
#  PURPOSE: Frame and widget used in the kratos main window to manage model/material
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
#   1.8- 09/11/13- GSM, change Ecape by Escape
#   1.7- 16/07/13- GSM, modify the proc ShowPropertyBySectionType to enable the properties as a function of the selected element type and section type
#                       - modify the proc buildTabFrame to show the thickness as a function of the selected element type
#   1.6- 15/07/13- GSM, set filterGroups to [::KMProps::GetAvailableGiDGroups] to use only normal or disabled GiD group 
#   1.5- 04/07/13- A.Melendo, created short_name_autogroup and long_name_description
#   1.4- 18/06/13- GSM, delete the use of the proc kipt::NewGiDGroups 
#   1.3- 13/05/13- G. Socorro, add the ShowPropertyBySectionType, update the proc buildGroupsFrame
#                               - modify the proc buildPropertyFrame to show the help, add translation option
#                               - add cross section database options
#   1.2- 05/12/12-  J. Garate,  Add Pick Coordinates Button to Tab frame
#   1.1- 28/11/12-  J. Garate,  Add Pick Coordinates Button to Groups Template frame
#   1.0- 10/10/12-  J. Garate,  Adaptation for New GiD Groups, Autogroup Frame Bug Corrected
#   0.9- 10/10/12-  G. Socorro, add the procs GetCrossSectionPropertyList and ShowPropertyByElementType
#                               update others procs to include the cross property functionality
#   0.8- 09/10/12-  G. Socorro, correct a bug when editing the tree in Linux OS => Modify bind $f.cmb <KeyPress> by bind $f.cmb <Return>
#   0.7- 27/09/12-  J. Garate,  Pick Coordinates Button
#   0.6- 26/04/12-  G. Socorro, change GiD_Groups by Cond_Groups
#   0.5- 12/04/2012 JGarate,    Ahora selecciona por defecto el primer elemento en Geometry Auto Group
#   0.4- 10/04/2012 JGarate,    Cambio de iconos para AutoGroup. Corregidos los problemas en mesh con algunos elementos
#   0.3- 04/04/2012 JGarate,    Corregido el error en mallado, por el cual no se permitia usar el Auto New Group, porque no mostraba los botones de seleccion
#   0.2- 02/04/2012 JGarate,    labels and combobox grouped left. Combobox width edited. Autogroup button changed
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
    grid $f -row 2 -column 0 -sticky wes -columnspan 2
    
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
    # wa "buildGroupsFrame fullname $fullname"
    set f [::KMProps::CreateBottomFrame]
    
    # El combo de propiedades es importante resetearlo
    # (posteriormente se utiliza para validar si hay alguna seleccionada)
    if {[info exists ::KMProps::cmbProperty]} {
	unset ::KMProps::cmbProperty
    }
    
    # Parametro para utilizar o no los combos "Activation"
    set activation 0
    set listT [::KMProps::getTemplateStructure $idTemplate]
    
    if {[llength $listT] >= 1 } {
	set nb $f.nb
	if {[winfo exists $nb]} {
	    destroy $nb
	}
	
	set long_name_description [$T item text $item 0]
	set parent_item [$T item parent $item]
	set parent_fullname [DecodeName [$T item tag names $parent_item]]
	set parent_class [::xmlutils::setXml $parent_fullname class]
	
	if { $parent_class=="SameTemplateGroups" } {
	    set long_name_description "[$T item text $parent_item 0] $long_name_description"      
	}  
	set short_name_autogroup ""
	if { [string length $long_name_description] < 10 } {
	    set short_name_autogroup $long_name_description
	} else {
	    foreach word $long_name_description {
		append short_name_autogroup [string range $word 0 3]      
	    }
	    if { [string length $short_name_autogroup] > 20 } {
		set short_name_autogroup [string range $short_name_autogroup 0 20] 
	    }
	}
	
	grid [ ttk::notebook $nb ] -sticky nwe
	
	
	# Lista de listas con formato {idContainer idItem1 idItem2...}
	foreach listContainer $listT {
	    
	    # Si tiene como m�nimo el container y un item ponemos tab y dentro label-combos
	    if {[llength $listContainer] >= 2} {
		
		set idContainer [lindex $listContainer 0]
		
		# Si solo hay un container le damos el nombre del item pulsado, no del template
		if {[llength $listT] == 1} {
		    set pid $long_name_description
		} else {
		    set pid [= [::KMProps::getPropTemplate $idTemplate pid $idContainer]]
		    # Lo limitamos por si es demasiado largo
		    set pid "[string range $pid 0 20]"
		}
		
		# Para cada container declaramos un tab
		set fTab ${nb}.f$idContainer
		$nb add [ttk::labelframe $fTab -text "[= Properties]" -padding {10 0 10 10}] \
		    -text "$pid"
		
		for {set i 1} { $i < [llength $listContainer] } {incr i} {
		    
		    set id [lindex $listContainer $i]
		    
		    # Si no coincide la dimensi�n 2D/3D no lo ponemos
		    set nDim [::KMProps::getPropTemplate $idTemplate nDim "$idContainer//$id"]
		    if { $nDim == "" || $nDim == $::KMProps::nDim } {
		        
		        # Comprobamos el estado para ver si est� activo o no
		        set state [::KMProps::getPropTemplate $idTemplate state "$idContainer//$id"]
		        
		        if {$state != "hidden"} {
		            
		            set pid [::KMProps::getPropTemplate $idTemplate pid "$idContainer//$id"]
		            set tooltip [::KMProps::getPropTemplate $idTemplate tooltip "$idContainer//$id"]
		            set dv [::KMProps::getPropTemplate $idTemplate dv "$idContainer//$id"]
		            set function [::KMProps::getPropTemplate $idTemplate function "$idContainer//$id"]
		            
		            # Obtenemos la lista de valores para el combo si existe
		            set comboList [::xmlutils::getXMLValues "$idContainer//$id" $idTemplate "" $fullname]
		            set icomboList [::xmlutils::getXMLValues "$idContainer//$id" $idTemplate "ivalues" $fullname]
		            
		            set CBState [::KMProps::getPropTemplate $idTemplate CBState "$idContainer//$id"]
		            if { $CBState == "normal" } {
		                set values $comboList
		                set comboList {}
		            } else {
		                set values {}
		            }
		            
		            # Para cada item a�adimos label y combo
		            grid [ttk::label $fTab.lbl$id -text "$pid:" ] -row $i -column 0 -pady 2 -sticky nw
		            
		            if { [llength $comboList] > 0 } {
		                
		                if {$state != "disable"} {
		                    set state "readonly"
		                }
		                if {[string length $id] == 2 && \
					([string index $id end] == "x" || [string index $id end] == "y" ||\
					     [string index $id end] == "z") } { 
		                    set width 15
		                } else {
		                    set width 20
		                }
		                grid [ttk::combobox $fTab.cmb$id -values $comboList \
					  -state $state -width [::KMProps::getCmbWidth $comboList] \
					  -textvariable "::KMProps::cmb$id"] \
		                    -row $i -column 1 -padx 3 -pady 2 -sticky nw
		                tooltip::tooltip $fTab.cmb$id [= $tooltip]
		                
		                if {$id == "Ax" || $id == "Ay" || $id == "Az"} {
		                    bind $fTab.cmb$id <<ComboboxSelected>> [list ::KMProps::cmbDisable $fullname $f.nb $id]
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
		                    -row $i -column 1 -padx 3 -pady 2 -sticky nw
		                tooltip::tooltip $fTab.cmb$id [= $tooltip]
		                
		                #set dv [::KMProps::getPropTemplate $idTemplate dvText "$idContainer//$id"]
		                set ::KMProps::cmb$id $dv
		                
		                # Pick Coordinates Button
		                set psb [::KMProps::getPropTemplate $idTemplate sbi "$idContainer//$id"]
		                if { $psb == "PickCoordinates" } {
		                    # msg "pick for $id"
		                    set sbxp [::KMProps::getPropTemplate $idTemplate sbxp "$idContainer//$id"]
		                    set sbp [::KMProps::getPropTemplate $idTemplate sbp "$idContainer//$id"]
		                    grid [ttk::button $fTab.btn$id -text $sbp -command [list ::KMProps::selectionButton $sbxp] -width 4 ] \
					-row $i -column 2 -padx 5 -pady 2 -sticky ew
		                    tooltip::tooltip $fTab.btn$id [= "Pick Coordinates"]
		                }
				
		            }
		            
			    
		            if {$function != "" } {
		                grid [ttk::button $fTab.funct$id -text "functions" -command "KFun::InitBaseWindow $fTab $id" -style TMenubutton.Toolbutton -width [::KMProps::getCmbWidth $comboList]]  \
		                    -row $i -column 2 -sticky nw  -pady 0 -padx 3
		                tooltip::tooltip $fTab.funct$id [= "Function manager"]
		                set img [::WinUtils::GetImage "functions.gif"]
		                if { $img != -1 } { $fTab.funct$id configure -image $img }
		                
		                grid [ttk::button $fTab.deleteFunct$id -text "delete" -command "::KMProps::unassignFunction $fTab $id $fullname" -style TMenubutton.Toolbutton -width [::KMProps::getCmbWidth $comboList]] \
		                    -row $i -column 3 -sticky nw  -pady 0 -padx 2
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
	
	#Miramos si hay restricci�n de entidades
	set entityList [split [::xmlutils::setXml $fullname GiDEntity] "," ]
	
	#set entityList [split [::KMProps::getPropTemplate $idTemplate GiDEntity] "," ]

	set fet [ ttk::frame $f.fEntityTypes]

	set idx_col 0

	switch $whatuse {
	    
	    GEOMETRYUSE {
		foreach i $geomlist {
		    #::KMProps::CorrectFileExtensionTo ".png" $i
		    
		    set command [ file rootname $i]
		    if {$command in $entityList || $entityList == "" } {
		        set fb $fet.b$command
		        grid [ ttk::button $fb -text $i \
				   -command [ list ::KMProps::changeImage $command $i $fet] \
				   -style TMenubutton.Toolbutton -width 15] \
		            -sticky nw -pady 3 -padx 5 \
			    -row 0 -column $idx_col
			incr idx_col
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
		            ::KMProps::changeImage $::KMProps::selectedEntity $i $fet
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
		        
		        set fb $fet.b$command
		        grid [ ttk::button $fb -text $i \
				   -command [ list ::KMProps::changeImage $command $i $fet] \
				   -style TMenubutton.Toolbutton -width 15] \
		            -sticky nw -pady 3 -padx 5 \
			    -row 0 -column $idx_col
			incr idx_col
		        if { [string range $command 0 6] == "element" } {
		            set command [string range $command 8 end]
		            #append i "element_" $i
		        }
		        
		        tooltip::tooltip $fb [= "Entity %s" $command]
		        set img [::WinUtils::GetImage $i]
		        
		        if { $img != -1 } {
		            $fb configure -image $img
		        }
		        incr col
		    }
		}
	    }
	}

	grid $fet -sticky new

	set fg [ ttk::frame $f.fGroups]

	ttk::label $fg.lGroups -text [= "Group"]:
	
	#COMBO DE GRUPOS
	set fGroups $fg.cGroups
	
	# Get the group list
	set filterGroups [::KMProps::GetAvailableGiDGroups]
	
	ttk::combobox $fGroups -state readonly -values "$filterGroups" \
	    -textvariable "::KMProps::selGroup" \
	    -postcommand "::KMProps::changeGroups [list $entityList] $fGroups" \
	    -width 15
	# -pady 3 -padx 55

	set ::KMProps::selGroup ""
	#Better always start frame with empty combo
	#if { [llength $filterGroups] > 0 } {
	#    $fg.cGroups current 0
	#}
	bind $fg.cGroups <<ComboboxSelected>> [list ::KMProps::cmbChangeCheckGroups $fg]
	
	# BOTON A LA DERECHA DE LOS GRUPOS (CREAR GRUPO AUTOMATICAMENTE)

	# just define button frame ( ok, cancel) so that we can pass it to the autogroup command
	set fbut [ ttk::frame $f.fbuttons]

	ttk::button $fg.iGroups -text [= "newGroup"] \
	    -command [list ::KMProps::autoNewGroup [list $short_name_autogroup] $fbut] 
	#  -pady 0 -padx 180
	tooltip::tooltip $fg.iGroups [= "Create automatic new group"]
	$fg.iGroups configure -image [::WinUtils::GetImage "newAutoGroup.gif" ]

	grid $fg.lGroups -row 0 -column 0 -sticky ne -padx 5 -pady 3
	grid $fGroups -row 0 -column 1 -sticky wne  -padx 5 -pady 3
	grid $fg.iGroups -row 0 -column 2 -sticky nw -padx 5 -pady 3
	grid rowconfigure  $fg 0 -weight 1
	grid columnconfigure  $fg 1 -weight 1

	grid $fg -sticky new

	#$fbut.bBottomOk and $fbut.bBottomCancel names need to popup click ok cancel before change item
	ttk::button $fbut.bBottomOk -text [= "Ok"] \
	    -command "::KMProps::acceptGroups $T $idTemplate $fullname $item {$listT} {$entityList} $fGroups"
	tooltip::tooltip $fbut.bBottomOk [= "Assign condition to the selected group"]
	
	ttk::button $fbut.bBottomCancel -text [= "Cancel"] \
	    -command "::KMProps::DestroyBottomFrame"
	tooltip::tooltip $fbut.bBottomCancel [= "Cancel assignation"]

	grid $fbut.bBottomOk  $fbut.bBottomCancel -sticky n -padx 5 -pady 3

	grid $fbut -sticky news
	if { $::tcl_version >= 8.5 } { grid anchor $fbut center }

	# $nb
	# $fet
	# $fg
	# $fbut
	grid columnconfigure $f 0 -weight 1
	grid rowconfigure $f {0 1 2 3} -weight 1
    }
    
    bind $T <Escape> [list ::KMProps::DestroyBottomFrame]
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
    
    #Leemos el xml si aun no se ha le�do
    #::xmlutils::initKKWord
    
    if { [info exists ::KMProps::ElemTypeThickness] } {
	set dv $::KMProps::ElemTypeThickness
	unset ::KMProps::ElemTypeThickness
    } else {
	set xpath "Kratos_KWords/ElementCLaws/Item\[@id='ElementTypes$::KMProps::nDim'\]"
	#set xpath "Kratos_KWords/ElementCLaws/Item\[@id='ElementTypes'\]"
	set dv [::xmlutils::getAttribute $KPriv(xmlDocKKW) $xpath dv]
    }
    
    set xpath "Kratos_KWords/ElementCLaws/Item\[@id='Tickness$::KMProps::nDim'\]"
    set ListElementType [split [::xmlutils::getAttribute $KPriv(xmlDocKKW) $xpath elementType] ","]
    # wa "dv:$dv ListElementType:$ListElementType"
    if { $dv in $ListElementType } {
	return 1
    } else {
	return 0
    }
}

proc ::KMProps::GetCrossSectionPropertyList {} {
    # ABSTRACT: Return the cross section property list
    # Result
    # Crosss section property list
    global KPriv
    
    # Xpath to all defined the crosss section property list
    set xpath "Kratos_KWords/ElementCLaws/Item\[@id='CSProperty'\]"
    set PropertyList [split [::xmlutils::getAttribute $KPriv(xmlDocKKW) $xpath "propertylist"] ","]

    return $PropertyList
} 

proc ::KMProps::ShowPropertyByElementType {propertyid} {
    # ABSTRACT: Show some properties as a function of the select element type    
    # Arguments
    # propertyid  -> Property Id in the xml file
    # Result
    # 1 -> Show this property
    # 0 -> Hide this property
    global KPriv
    
    set ndime "2D"
    if {[info exists ::KMProps::nDim]} {
	set ndime $::KMProps::nDim
    } else {
	# Get the spatial dimension
	set ndime [::xmlutils::GetSpatialDimension]
	# wa "ndime:$ndime"
    }
    
    # Current edited properties
    if {[info exists ::KMProps::ElemTypeProperty] } {
	set dv $::KMProps::ElemTypeProperty
	unset ::KMProps::ElemTypeProperty
    } else {
	# Default case
	# Xpath to all defined element type
	set xpath "Kratos_KWords/ElementCLaws/Item\[@id='ElementTypes$::KMProps::nDim'\]"
	#set xpath "Kratos_KWords/ElementCLaws/Item\[@id='ElementTypes'\]"
	set dv [::xmlutils::getAttribute $KPriv(xmlDocKKW) $xpath dv]
    } 
    # wa "dv:$dv"
    if {$dv !=""} {
	# Get the property from the Kratos keyword mapping file
	# Process the special case of thickness
	if {$propertyid eq "Thickness"} {
	    set propertyid "${propertyid}${ndime}"
	}
	set xpath "Kratos_KWords/ElementCLaws/Item\[@id='${propertyid}'\]"
	
	# Get the element type list
	set ListElementType [split [::xmlutils::getAttribute $KPriv(xmlDocKKW) $xpath elementType] ","]
	# wa "propertyid:$propertyid ListElementType:$ListElementType"
	if { $dv in $ListElementType } {
	    return 1
	} else {
	    return 0
	}
    }
    return 0
}


proc ::KMProps::ShowPropertyBySectionType {propertyid {from CreateFrame}} {
    # ABSTRACT: Show some properties as a function of the select section type    
    # Arguments
    # propertyid  -> Property Id in the xml file
    # Result
    # 1 -> Show this property
    # 0 -> Hide this property
    global KPriv
    
    set ndime "2D"
    if {[info exists ::KMProps::nDim]} {
	set ndime $::KMProps::nDim
    } else {
	# Get the spatial dimension
	set ndime [::xmlutils::GetSpatialDimension]
	# wa "ndime:$ndime"
    }
    
    # Current edited properties
    if {[info exists ::KMProps::SectionTypeProperty] } {
	set dv $::KMProps::SectionTypeProperty
	unset ::KMProps::SectionTypeProperty
    } else {
	# Default case
	# Xpath to all defined section type
	set xpath "Kratos_KWords/ElementCLaws/Item\[@id='SectionTypes'\]"
	set dv [::xmlutils::getAttribute $KPriv(xmlDocKKW) $xpath dv]
    } 
    # wa "dv:$dv"
    if {$dv !=""} {
	
	if {$from eq "RTree"} {
	    # Get the element type properties
	    # Current edited properties
	    if {[info exists ::KMProps::ElemTypeProperty] } {
		set etypedv $::KMProps::ElemTypeProperty
	    } else {
		# Default case
		# Xpath to all defined element type
		set xpath "Kratos_KWords/ElementCLaws/Item\[@id='ElementTypes$::KMProps::nDim'\]"
		#set xpath "Kratos_KWords/ElementCLaws/Item\[@id='ElementTypes'\]"
		set etypedv [::xmlutils::getAttribute $KPriv(xmlDocKKW) $xpath dv]
	    } 
	    
	    # Get the property from the Kratos keyword mapping file
	    set xpath "Kratos_KWords/ElementCLaws/Item\[@id='${propertyid}'\]"
	    
	    set ListElementType [split [::xmlutils::getAttribute $KPriv(xmlDocKKW) $xpath elementType] ","]
	    # wa "etypedv:$etypedv ListElementType:$ListElementType"
	    
	    
	    # Get the section type list
	    set ListSectionType [split [::xmlutils::getAttribute $KPriv(xmlDocKKW) $xpath sectionType] ","]
	    
	    # wa "propertyid:$propertyid ListSectionType:$ListSectionType"
	    if {($dv in $ListSectionType) && ($etypedv in $ListElementType) } {
		return 1
	    } else {
		return 0
	    }
	} else {

	    # Get the property from the Kratos keyword mapping file
	    set xpath "Kratos_KWords/ElementCLaws/Item\[@id='${propertyid}'\]"
	    
	    # Get the section type list
	    set ListSectionType [split [::xmlutils::getAttribute $KPriv(xmlDocKKW) $xpath sectionType] ","]
	    
	    # wa "propertyid:$propertyid ListSectionType:$ListSectionType"
	    if {$dv in $ListSectionType} {
		return 1
	    } else {
		return 0
	    }
	}
    }
    return 0
}

#
# Construye un frame consultando el template de propiedades correspondiente,
#  con un tab por cada container y dentro de cada tab, 
#  etiquetas y combos para cada item (si los hay)
#
proc ::KMProps::buildPropertyFrame { T idTemplate item fullname } { 

    global KPriv
    
    # wa "T:$T idTemplate:$idTemplate item:$item fullname:$fullname"
    set f [::KMProps::CreateBottomFrame]

    set listT [::KMProps::getTemplateStructure $idTemplate]
    # wa "listT:$listT"
    if {[llength $listT] >= 1 } {
	set nb ${f}.nb
	if {[winfo exists $nb]} {
	    destroy $nb
	}
	grid [ttk::notebook $nb ] -row 0 -column 0 -columnspan 2 -padx 0 -sticky nw -in $f

	# Lista de listas con formato {idContainer idItem1 idItem2...}
	foreach listContainer $listT {

	    # Si tiene como m�nimo el container(1er elemento) y un item ponemos tab y dentro label-combos
	    if {[llength $listContainer] > 1} {

		set idContainer [lindex $listContainer 0]

		set pid [::KMProps::getPropTemplate $idTemplate pid $idContainer]
		# wa "idContainer:$idContainer pid:$pid cero:[lindex $listT 0 0]"
		# Para cada container declaramos un tab
		set fTab ${nb}.f$idContainer
		set ptxt "[= Properties]"
		set tabtxt "[= "%s" $pid]" 
		$nb add [ttk::labelframe $fTab -text "$ptxt" -padding {10 0 10 10}] \
		    -text "$tabtxt"

		# En el caso del primer tab forzamos el item de "Nombre de propiedad" y 2 campos mas
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

		    grid [ttk::combobox $fTab.cmbPropertyName -state normal -textvariable "::KMProps::propertyName" -width 10] \
			-row 0 -column 1 -padx 3 -pady 5 -sticky nw -in $fTab

		    tooltip::tooltip $fTab.cmbPropertyName $cbhelp

		    set ::KMProps::propertyName "[::KEGroups::GetAutomaticPropertyName $fullname "$whatoption"]"

		    focus $fTab.cmbPropertyName
		}
		for {set i 1} { $i < [llength $listContainer] } {incr i} {

		    set id [lindex $listContainer $i]
		    # wa "id:$id"
		    # Los nodos ocultos no se deben mostrar
		    set state [::KMProps::getPropTemplate $idTemplate state "$idContainer//$id"]
		    # wa "state:$state"
		    
		    if {$state != "hidden" } {

		        # Si no coincide la dimensi�n 2D/3D no lo ponemos
		        set nDim [::KMProps::getPropTemplate $idTemplate nDim "$idContainer//$id"]
			# wa "nDim:$nDim"
		        if { $nDim == "" || $nDim == $::KMProps::nDim } {

		            set pid [::KMProps::getPropTemplate $idTemplate pid "$idContainer//$id"]
		            set tooltip [::KMProps::getPropTemplate $idTemplate help "$idContainer//$id"]
		            set dv [::KMProps::getPropTemplate $idTemplate dv "$idContainer//$id"]

		            set CBState [::KMProps::getPropTemplate $idTemplate CBState "$idContainer//$id"]
			    # wa "pid:$pid tooltip:$tooltip dv:$dv CBState:$CBState" 
		            if { $CBState == "normal" } {
		                set values $comboList
		                set comboList {}
		            } else {
		                set values {}
		            }

		            # Es importante que ElemType sea el primer combo para actualizar el "dv" en el xml
		            if {($id eq "ElemType")||($id eq "SectionType")} {
		                set comboList [::xmlutils::getXMLValues "$idContainer//$id" $idTemplate "" "$fullname" $dv]
		            } else {
		                set comboList [::xmlutils::getXMLValues "$idContainer//$id" $idTemplate "" "$fullname"]
		            }
		            set cdv [::KMProps::getPropTemplate $idTemplate dv "$idContainer//$id"]
			    # wa "after cdv:$cdv"

		            # Para cada item a�adimos label y combo
		            set lpath "$fTab.lbl$id"
		            if {[winfo exists $lpath]} {
		                destroy $lpath
		            }
		            grid [ttk::label $lpath -text "$pid:" ] \
				-row $i -column 0 -pady 5 -sticky nw -in $fTab


		            # Init the global combobox variable
		            set varid "::KMProps::cmb$id"
		            if {![info exists $varid]} {
		                set $varid ""
		            }

		            # Destroy the current combobox
			    set cbpath "$fTab.cmb$id"
			    if {[winfo exists $cbpath]} {
		                destroy $cbpath
			    }

			    # wa "comboList:$comboList"
		            if {[llength $comboList]} {

		                grid [ttk::combobox $cbpath -values $comboList -state readonly -width [::KMProps::getCmbWidth $comboList] \
					  -textvariable "$varid" \
					  -postcommand [list ::KMProps::changeCmbValues "$cbpath" "$idContainer//$id" "$idTemplate" "" "$fullname"] ] \
				    -row $i -column 1 -padx 5 -pady 2 -sticky nw -in $fTab
		                
		                # Set the combobox values

		                ::xmlutils::setComboDv $cbpath $fullname $dv $idTemplate
		                
		                # Set the help text
		                tooltip::tooltip $cbpath $tooltip

		                set psb [::KMProps::getPropTemplate $idTemplate sbi "$idContainer//$id"]
				

		                # En este caso se tendr� qu� recargar si existe el combo de "Material Model"
		                if { $id == "ElemType" } {
		                    bind $cbpath <<ComboboxSelected>> [list ::KMProps::cmbElemTypeChange $cbpath $idContainer//$id $idTemplate]
		                }

		                # Get the section type
		                if {$id eq "SectionType"} {
		                    bind $cbpath <<ComboboxSelected>> [list ::KMProps::cmbSectionTypeChange $cbpath $idContainer//$id $idTemplate]
		                }

		            } else {

		                # Create the combobox
		                grid [ttk::combobox $cbpath -state normal -width [::KMProps::getCmbWidth $comboList] -values $values -textvariable "$varid"] \
				    -row $i -column 1 -padx 5 -pady 2 -sticky nw -in $fTab

		                # Set the current value
		                set $varid $dv
		                
		                # Set the help text 
		                tooltip::tooltip $fTab.cmb$id $tooltip

		                
		            }
			    
			    # Pick Coordinates Button
			    set psb [::KMProps::getPropTemplate $idTemplate sbi "$idContainer//$id"]
			    if { $psb == "PickCoordinates" } {
				# msg "pick for $id"
				set sbxp [::KMProps::getPropTemplate $idTemplate sbxp "$idContainer//$id"]
				set sbp [::KMProps::getPropTemplate $idTemplate sbp "$idContainer//$id"]
				grid [ttk::button $fTab.btn$id -text $sbp -command [list ::KMProps::selectionButton $sbxp] -width 4 ] \
		                    -row $i -column 2 -padx 5 -pady 2 -sticky ew -in $fTab
				tooltip::tooltip $fTab.btn$id [= "Pick Coordinates"]
			    }

		            # In the case of ElemType update the property filter
		            # Get the cross section property list
		            set PropertyList [::KMProps::GetCrossSectionPropertyList]
			    # wa "buildPropertyFrame =>PropertyList:$PropertyList"
		            if { $id == "ElemType" } {
		                set ::KMProps::ElemTypeProperty $dv
				
		            } elseif { $id == "SectionType" } {
		                set ::KMProps::SectionTypeProperty $dv

		            } elseif {$id in $PropertyList} {
		                
		                # Get the current element type
		                set cdv [::KMProps::getPropTemplate $idTemplate dv "$idContainer//ElemType"]
		                # wa "cdv:$cdv"
		                set ::KMProps::ElemTypeProperty $cdv
		                if {![::KMProps::ShowPropertyByElementType $id]} {
		                    # Remove some properties
		                    grid remove $fTab.lbl$id
		                    grid remove $fTab.cmb$id
		                }
				
		                # Get the current section type
		                set cdv [::KMProps::getPropTemplate $idTemplate dv "$idContainer//SectionType"]
				# wa "cdv:$cdv"
		                set ::KMProps::SectionTypeProperty $cdv
		                if {![::KMProps::ShowPropertyBySectionType $id]} {
		                    # Remove some properties
		                    grid remove $fTab.lbl$id
		                    grid remove $fTab.cmb$id
		                }
		            }
		        }
		    }
		}
	    }
	}
    }                                                

    # First delete buttons    
    set wdlist [list $f.bBottomOk $f.bBottomCancel]
    foreach wd $wdlist {
	if {[winfo exists $wd]} {
	    destroy $wd
	}
    }
    #$f.bBottomOk and $f.bBottomCancel names need to popup click ok cancel before change item
    grid [ttk::button $f.bBottomOk -text [= "Ok"]  -command "::KMProps::acceptProperty $T $idTemplate $fullname $item {$listT}" ] \
	-row 3 -column 0 -sticky nw  -pady 3 -padx 20  -in $f

    tooltip::tooltip $f.bBottomOk [= "Assign the defined properties"]

    grid [ttk::button $f.bBottomCancel -text [= "Cancel"]  -command "::KMProps::DestroyBottomFrame" ] \
	-row 3 -column 0 -sticky nw  -pady 3 -padx 100  -in $f
    tooltip::tooltip $f.bBottomCancel [= "Cancel the defined properties"]

    bind $T <Escape> [list ::KMProps::DestroyBottomFrame]

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
# Click en bot�n OK de Property
#
proc ::KMProps::acceptProperty { T idTemplate fullname item listT} {
    
    # wa "T:$T idTemplate:$idTemplate fullname:$fullname item:$item listT:$listT"

    # Get the property identifier
    
    set property $::KMProps::propertyName
    
    # Validamos que la propiedad no tenga car�cteres extra�os
    set property [::KUtils::parseTreeStr $property]
    if { $property == -1 } {
	WarnWin [= "You can't use some reservate chars like:\n  :   /   $   .   \\  %  "]
	set ::KMProps::propertyName ""
	return ""
    }
    
    # Comprobamos q el nombre no sea vac�o
    if { $property == "" } {
	WarnWin [= "The property name can not be empty"]
	return ""
    }
    
    # Comprobamos q la propiedad aun no exista
    set id [::xmlutils::setXml "$fullname//c.[list $property]" id]

    if { $id != "" } {
	WarnWin [= "This name property it is already assigned."]
	return
    }

    # Copy the template 

    ::KMProps::copyTemplate ${idTemplate} $fullname "$property" "Property"
    
    #
    # Ahora debemos actualizar todos los valores en el xml
    #
    if {[llength $listT] >= 1 } {
	
	# Lista de listas con formato {idContainer idItem1 idItem2...}
	foreach listContainer $listT {
	    
	    # Si tiene como m�nimo el container y un item entramos
	    if {[llength $listContainer] >= 2} {
		
		set idContainer [lindex $listContainer 0]
		# wa "idContainer:$idContainer"
		
		# Recorremos los items
		for {set i 1} { $i < [llength $listContainer] } {incr i} {
		    
		    set id [lindex $listContainer $i]
		    # wa "id:$id"
		    set fullNombre "$fullname//c.[list $property]//c.[list $idContainer]//i.[list $id]"
		    # wa "fullNombre:$fullNombre"
		    
		    # Los nodos ocultos no exist�an en el formulario
		    set state [::xmlutils::setXml $fullNombre state]
		    # wa "state:$state"
		    if {$state != "hidden" } {
		        
		        # Validamos la dimensi�n de cada elemento
		        set nDim [::xmlutils::setXml $fullNombre nDim]
			# wa "nDim:$nDim ::KMProps::nDim:$::KMProps::nDim"
		        if { $nDim == "" || $nDim == $::KMProps::nDim } {
		            
		            # Get the value from the internal variable
		            
		            set value [set ::KMProps::cmb$id]
			    # wa "value:$value"
		            # Update the value in the xml file
		            
		            ::xmlutils::setXml $fullNombre dv "write" $value
		        }
		    }
		}
	    }
	}
    }

    # Destroy the buttom frame
    ::KMProps::DestroyBottomFrame
    
    # Reload the tree
    ::KMProps::RefreshTree $T
}

proc ::KMProps::changeImage {entity img path} {
    
    variable selectedEntity
    
    if {$selectedEntity != "" } {
	
	set f $path.b$selectedEntity
	
	$f configure -image [::WinUtils::GetImage $selectedEntity.gif]
    }
    
    $path.b$entity configure -image [::WinUtils::GetImage ${entity}_sel.gif]

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
    # wa "children:$children item:$item class:$class"

    set listTabs {}
    set listItems {}
    set acceptItems {}
    foreach itemChild $children {

	set fullname [DecodeName [$T item tag names $itemChild]]

	set nodeName [::xmlutils::getXmlNodeName $fullname]

	# Miramos si cada hijo es container o item
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

    # Reset the variable used for this property
    set ::KMProps::propertyName ""

    # Si no tiene containers pero tiene items, utilizamos como tab el elemento padre seleccionado
    if { [llength $listTabs] == 0 && [llength $listItems] >= 1 } {

	set fullname [DecodeName [$T item tag names $item]]
	lappend listTabs $fullname
    } else {
	set listItems {}
    }
    # wa "listTabs:$listTabs listItems:$listItems acceptItems:$acceptItems"

    if { [llength $listTabs] >= 1 } {

	# Delete the notebook widget

	set nb ${f}.nb
	if {[winfo exists $nb]} {
	    destroy $nb
	}
	grid [ttk::notebook $nb ] -row 0 -sticky ewn

	set i 0
	foreach fullname $listTabs {

	    set id [::xmlutils::setXml $fullname id]
	    set pid [::xmlutils::setXml $fullname pid]                                
	    set help [::xmlutils::setXml $fullname help]

	    # wa "id:$id pid:$pid help:$help"

	    # For each container declare a new tab
	    set fTab ${nb}.f$id

	    set tabtxt "[= "%s" $pid]"
	    set ftxt [= "Properties"]
	    $nb add [ttk::labelframe $fTab -text "$ftxt" -padding {10 10 10 10}] -text "$tabtxt"

	    if {[llength $listItems] > 0 } {
		#En este caso en realidad los nietos son los propios hijos del item pulsado
		set nietos [$T item children $item]
	    } else {
		#Miramos los hijos de cada container (nietos del pulsado originalmente)
		set itemChild [lindex $children $i]
		set nietos [$T item children $itemChild]
	    }

	    # Property name => Add the property identifier
	    if {($class == "Property") && ($i == 0) } {

		# Property name/identifier
		set txt "[= "Property Id"]:"
		grid [ttk::label $fTab.lblName -text "$txt" ] \
		    -row 0 -column 0 -pady 5 -sticky nw -in $fTab

		grid [ttk::combobox $fTab.cmbPropertyName -state normal -textvariable "::KMProps::propertyName" -width 15] \
		    -row 0 -column 1 -padx 3 -pady 5 -sticky nw -in $fTab

		set ::KMProps::propertyName "[$T item text $item 0]"
		

	    } elseif { ($class == "Group") && ($i == 0) } {

		######################
		#COMBO DE GRUPOS
		set fullname [DecodeName [$T item tag names [$T item parent $item]]]

		set ::KMProps::selGroup "[$T item text $item 0]"

		set entityList [split [::xmlutils::setXml $fullname GiDEntity] ","]
		
		# Get the group list               
		set filterGroups [::KMProps::GetAvailableGiDGroups]
		grid [ttk::label $fTab.lblName -text "[= Group:]" ] \
		    -row 0 -column 0 -pady 5 -sticky nw -in $fTab

		set fGroups $fTab.cGroups
		grid [ttk::combobox $fGroups -state readonly -values "$filterGroups" -textvariable "::KMProps::selGroup"  -postcommand "::KMProps::changeGroups [list $entityList] $fGroups $fullname" -width 15] \
		    -row 0 -column 1 -pady 5 -sticky nw -in $fTab



		bind $fTab.cGroups <<ComboboxSelected>> [list ::KMProps::cmbChangeCheckGroups $fTab]
		######################
	    }

	    set row 1
	    foreach nieto $nietos {

		set fullname [DecodeName [$T item tag names $nieto]]
		set dv [::xmlutils::setXml $fullname dv]
		
		# wa "fullname:$fullname nieto:$nieto dv:$dv"

		# Check for item
		set NodeName [::xmlutils::getXmlNodeName $fullname] 
		if {$NodeName eq "Item" } {

		    set state [::xmlutils::setXml $fullname state]
		    # wa "state:$state"

		    if {$state != "hidden"} {

		        set id [::xmlutils::setXml $fullname id]
		        set pid [::xmlutils::setXml $fullname pid]
		        set tooltip [::xmlutils::setXml $fullname help]
		        set function [::xmlutils::setXml $fullname function]
		        # wa "id:$id pid:$pid tooltip:$tooltip function:$function"
		        
		        # Store each item to update the values when finish
		        lappend acceptItems $nieto $fTab

		        set comboList [::xmlutils::getXMLValues "$fullname"]
		        # wa "comboList:$comboList"

		        set CBState [::xmlutils::setXml $fullname CBState]
		        if { $CBState == "normal" } {
		            set values $comboList
		            set comboList {}
		        } else {
		            set values {}
		        }

		        # For each item add the label and combobox
		        grid [ttk::label $fTab.lbl$id -text "${pid}:" ] \
			    -row $row -column 0 -padx 3 -pady 5 -sticky nw -in $fTab

		        # Set the combobox path
		        set cbpath $fTab.cmb$id
		        if {[winfo exists $cbpath]} {
		            destroy $cbpath
		        }

		        # Create the combobox variable and init it
		        set varid "::KMProps::cmb$id"
		        if {![info exists $varid]} {
		            set $varid "" 
		        }

		        if { [llength $comboList] > 0 } { 

		            if {$state != "disabled"} {
		                set state "readonly"
		            }
			    
		            # Create the combobox
		            grid [ttk::combobox $cbpath -values $comboList -state $state \
				      -textvariable "$varid" -width [::KMProps::getCmbWidth $comboList]]\
				-row $row -column 1 -padx 3 -pady 5 -sticky nw -in $fTab 

		            # Set the help 
		            ::tooltip::tooltip $cbpath [= "%s" $tooltip]

		            # Set the combobox values
		            ::xmlutils::setComboDv $cbpath $fullname $dv

		            if {$id == "Ax" || $id == "Ay" || $id == "Az"} {
		                ::KMProps::cmbDisable $fullname $f.nb
		                bind $fTab.cmb$id <<ComboboxSelected>> [list ::KMProps::cmbDisable $fullname $f.nb]
		            } elseif {$id eq "ElemType" } {
		                # En este caso se tendr� qu� recargar si existe el combo de "Material Model"
		                bind $cbpath <<ComboboxSelected>> [list ::KMProps::cmbElemTypeChange $cbpath $fullname]
		            } elseif {$id eq "SectionType" } {
		                # Update the section type combobox
		                bind $cbpath <<ComboboxSelected>> [list ::KMProps::cmbSectionTypeChange $cbpath $fullname]
		            }

		        } else {
		            if {$state != "disabled"} {
		                set state "normal"
		            }

		            # Create the combobox
		            grid [ttk::combobox $cbpath -state $state -values $values -textvariable "$varid" -width [::KMProps::getCmbWidth $comboList]] \
				-row $row -column 1 -padx 3 -pady 5 -sticky nwe -in $fTab 

		            # Set the help 
		            ::tooltip::tooltip $cbpath [= "%s" $tooltip]

		            # Update variable value from the current value of dv 
		            set $varid $dv
		            
		        }
			
			# Pick Coordinates Button
		        set psb [::xmlutils::setXml $fullname sbi]
		        if { $psb eq "PickCoordinates" } {
		            # msg "pick for $id"
		            #[::xmlutils::setXml $fullname sbp]
		            set sbxp [::xmlutils::setXml $fullname sbxp]
		            set sbp [::xmlutils::setXml $fullname sbp]
		            grid [ttk::button $fTab.btn$id -text $sbp -command [list ::KMProps::selectionButton $sbxp] -width 4 ] \
				-row $row -column 2 -padx 5 -pady 2 -sticky ew -in $fTab
		            tooltip::tooltip $fTab.btn$id [= "Pick Coordinates"]
		        }

		        # Create the python function

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

		        # In the case of ElemType update the property filter
		        set PropertyList [::KMProps::GetCrossSectionPropertyList]
		        # wa "buildTabFrame => PropertyList:$PropertyList"
		        if { $id == "ElemType" } {
		            set ::KMProps::ElemTypeProperty $dv

			} elseif {$id in $PropertyList} {
		            # Get the parent node
		            set ParentNode [$T item parent $nieto]
		            set fpath [DecodeName [$T item tag names $ParentNode]]
		            set CurrentPropertyId [::xmlutils::setXml $fpath id]
		            # wa "CurrentPropertyId:$CurrentPropertyId"

		            # Get the select element base type
		            set ::KMProps::ElemTypeProperty [::xmlutils::GetPropertyElemType $CurrentPropertyId]
		            # wa "ElemTypeProperty:$::KMProps::ElemTypeProperty"

			    # Get the select section base type
			    set ::KMProps::SectionTypeProperty [::xmlutils::GetPropertySectionType $CurrentPropertyId]
			    # wa "SectionTypeProperty:$::KMProps::SectionTypeProperty"

			    # Special case of thickness
			    if {$id eq "Thickness"} {
				set ElemShowProperty [::KMProps::ShowPropertyByElementType $id]
				# wa "ElemShowProperty:$ElemShowProperty"
				if {!$ElemShowProperty} {
				    # Remove some properties
				    grid remove $fTab.lbl$id
				    grid remove $fTab.cmb$id
				}
			    } else {
				set SectionShowProperty [::KMProps::ShowPropertyBySectionType $id]
				# wa "SectionShowProperty:$SectionShowProperty"
				if {!$SectionShowProperty} {
				    # Remove some properties
				    grid remove $fTab.lbl$id
				    grid remove $fTab.cmb$id
				}
			    }
		        }
		        
		        incr row
		    }
		}
	    }
	    incr i
	}
	
	# Si pulsan Esc tambi�n forzamos la salida del Tab
	bind $T <Escape> [list ::KMProps::DestroyBottomFrame]
	#bind $T <Return> [list ::KMProps::acceptTabFrame $T $acceptItems $class $item]
	#bind $::KMProps::WinPath <Return> [list ::KMProps::acceptTabFrame $T $acceptItems $class $item]

	if { [llength $acceptItems] } {
	    # First remove the button widget if it exists
	    set wdlist [list $f.bBottomOk $f.bBottomCancel]
	    foreach wd $wdlist {
		if {[winfo exists $wd]} {
		    destroy $wd
		}
	    }
	    #$f.bBottomOk and $f.bBottomCancel names need to popup click ok cancel before change item
	    grid [ttk::button $f.bBottomOk -text "Ok"  -command "[list ::KMProps::acceptTabFrame $T $acceptItems $class $item]" ] \
		-row 1 -column 0 -sticky sw  -pady 3 -padx 20  -in $f
	    tooltip::tooltip $f.bBottomOk [= "Confirm values"]

	    grid [ttk::button $f.bBottomCancel -text "Cancel"  -command "::KMProps::DestroyBottomFrame" ] \
		-row 1 -column 0 -sticky sw  -pady 3 -padx 100  -in $f
	    tooltip::tooltip $f.bBottomCancel [= "Cancel assignation"]
	}
    }
}

#
# Construye un combo din�mico en el item pulsado del �rbol y posteriormente se elimina
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
	# Configurar frame en funci�n del XML
	#---------------------------#---------------------------#
	set f "$T.f$idFull"
	
	if { [winfo exists $f] } {
	    return
	}
	# Remove the backgroud color to be used with the new GiD dark theme
	set bg "#F8F8F8"
	# set f [frame "$T.f$idFull" -borderwidth 0 -background $bg]
	set f [frame "$T.f$idFull" -borderwidth 0]
	
	if { [llength $comboList] > 0 } {
	    
	    grid [ttk::combobox $f.cmb -values $comboList -state readonly -width [::KMProps::getCmbWidth $comboList] -textvariable "::KMProps::cmb$idFull"] \
		-row 0 -column 0 -padx 3 -sticky nw -in $f 
	    
	    ::xmlutils::setComboDv $f.cmb $fullname $dv
	    #set selected [::xmlutils::getSelected $dv $comboList]
	    #$f.cmb current $selected
	    
	    bind $f.cmb <<ComboboxSelected>> [list ::KMProps::cmbSelectChange $item $T 1 current]
	} else {
	    
	    grid [ttk::combobox $f.cmb -state normal -values $values -textvariable "::KMProps::cmb$idFull" -width [::KMProps::getCmbWidth $comboList]] \
		-row 0 -column 0 -padx 3 -sticky nw -in $f
	    
	    set ::KMProps::cmb$idFull $dv
	    
	    bind $f.cmb <Leave> [list ::KMProps::cmbSelectChange $item $T 0 current]
	    bind $f.cmb <FocusOut> [list ::KMProps::cmbSelectChange $item $T 1 current]
	    bind $f.cmb <Escape> [list ::KMProps::cmbCancel $item $T]
	}
	# Si pulsan intro o Esc tambi�n forzamos la salida del combo
	bind $f.cmb <Return> [list ::KMProps::cmbSelectChange $item $T 1 current]
	bind $T <Escape> [list ::KMProps::cmbCancel $item $T]
	
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
	
	#Si el combo era especial se tendr� q reconstruir el �rbol
	set clase [::xmlutils::setXml $fullname class]
	
	#Comprueba si el cambio va a afectar a la estructura del arbol
	::KMProps::specialComboAction $T $clase $selCombo $item $id
    }
    
    #
    # Caso especial para las propiedades
    #
    if {$class == "Property" && $itemSel != "" } {
	
	set newPropertyName $::KMProps::propertyName
	
	#Validamos que la propiedad no tenga car�cteres extra�os
	set newPropertyName [::KUtils::parseTreeStr $newPropertyName]
	if { $newPropertyName == -1 } {
	    
	    WarnWin [= "You can't use some reservate chars like:\n  :   /   $   .   \\  %  "]
	    set ::KMProps::propertyName ""
	    return
	}
	
	#Comprobamos q el nombre no sea vac�o
	if { $newPropertyName == "" } {
	    WarnWin [= "The property name can not be empty"]
	    return
	}
	
	#Comprobamos q si la propiedad ha cambiado no exista el nuevo nombre
	set fullname [DecodeName [$T item tag names $itemSel]]
	set oldId [::xmlutils::setXml "$fullname" id]
	set repeatId [::xmlutils::setXml "[::KMProps::getApplication $fullname]//c.Properties//c.[list $newPropertyName]" id]
	
	if { $repeatId != "" && $oldId != $newPropertyName } {
	    WarnWin [= "This property name already exist."]
	    return
	}
	
	#Si todo es correcto, renombramos la propiedad en todo el �rbol y luego lo recargamos
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
    
    #Volvemos a cargar el �rbol para q el path de los items sea correcto
    ::KMProps::RefreshTree $T
    
    ::KMProps::DestroyBottomFrame
}

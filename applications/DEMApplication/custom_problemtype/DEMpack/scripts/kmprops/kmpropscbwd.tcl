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
#   1.0- 03/03/14-GSM, correct a bug in the proc DemMaterialSpecialInteraction change Active by Activate to hidden/normal the RollingFriction property
#   0.9- 18/06/13-GSM, delete the use of the proc kipt::NewGiDGroups 
#   0.8- 06/05/13-G. Socorro, add the option to work with the cross section properties (simple property or database)
#                             - modify and update some procedures
#   0.7- 13/12/12- J. Garate, add a message for the "PFEM" fluid case on old versions (11.1.2d)
#   0.6- 09/10/12- G. Socorro, update others procs to include the cross property functionality
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
    
    if {[winfo exists $f]} {
    
    set values [::xmlutils::getXMLValues $path $idTemplate "" "$noTemplateFullname" "$elemTypeDv"]
    
    $f configure -values $values
    
	if {([llength $values]) && ([$f current] == -1)} {
	$f current 0        
    }
}
}

proc ::KMProps::cmbElemTypeChange { f fullname {idTemplate ""} } {
    
    global KPriv
    
    # wa "ElemType f:$f fullname:$fullname idTemplate:$idTemplate"
    if { [info exists ::KMProps::cmbElemType] } {
	
	set dv [::xmlutils::getComboDv $f $fullname "id" $idTemplate]
	# wa "dv:$dv"

	set xpath "Kratos_KWords/ElementCLaws/Item\[@id='ElementTypes$::KMProps::nDim'\]"
	#set xpath "Kratos_KWords/ElementCLaws/Item\[@id='ElementTypes'\]"
	::xmlutils::getAttribute $KPriv(xmlDocKKW) $xpath dv $dv
	
	# We update node Material values based on the new dv
	set fMat [string map {"ElemType" "MatModel"} $f]
	
	if { $idTemplate == "" } {
	    
	    set fullnameMat [string map {"ElemType" "MatModel"} $fullname]
	    ::KMProps::changeCmbValues $fMat "$fullnameMat" $idTemplate "NoSearchFullname"
	} else {
	    
	    ::KMProps::changeCmbValues $fMat "MainProperties//MatModel" $idTemplate
	}
	
	# Get the cross section property list
	set PropertyList [::KMProps::GetCrossSectionPropertyList]
	# wa "PropertyList:$PropertyList"

	# For each property
	foreach propid $PropertyList {

	    # For element type
	    set fCmbPropId [string map [list "ElemType" "$propid"] $f]
	    set fLblPropId [string map [list "cmb${propid}" "lbl${propid}"] $fCmbPropId]
	    # wa "fCmbPropId:$fCmbPropId fLblPropId:$fLblPropId"
	    # Get the show property state for the current property identifier
	    set ShowProperty [::KMProps::ShowPropertyByElementType $propid]
	    # wa "ShowProperty by element type:$ShowProperty"
	    if {$ShowProperty} {
		# Show this property	
		if { [winfo exists $fCmbPropId] } { 
		    grid $fCmbPropId 
		}
		if { [winfo exists $fLblPropId] } { 
		    grid $fLblPropId 
		}
	    } else {
		# Hide this property
		if { [winfo exists $fCmbPropId] } { 
		    grid remove $fCmbPropId 
		}
		if { [winfo exists $fLblPropId] } { 
		    grid remove $fLblPropId 
		}
	    }
	    
	    set xpath "Kratos_KWords/ElementCLaws/Item\[@id='${propid}'\]"
	    # Get the element type list
	    set ListElementType [split [::xmlutils::getAttribute $KPriv(xmlDocKKW) $xpath elementType] ","]
	    # wa "propid:$propid ListElementType:$ListElementType"
	    if {($dv in $ListElementType) && ($dv ne "Truss") } {
		# Note: Special case of truss element => Without using section type
		
		# For section type
		set fCmbPropId [string map [list "ElemType" "SectionType"] $f]
		set fCmbPropId [string map [list "SectionType" "$propid"] $fCmbPropId]
		set fLblPropId [string map [list "cmb${propid}" "lbl${propid}"] $fCmbPropId]
		# wa "fCmbPropId:$fCmbPropId fLblPropId:$fLblPropId"
		
		# Get the show property state for the current property identifier => Using the section type
		set ShowProperty [::KMProps::ShowPropertyBySectionType $propid]
		# wa "propid:$propid ShowProperty:$ShowProperty"
		if {$ShowProperty} {
		    # Show this property	
		    if { [winfo exists $fCmbPropId] } { 
			grid $fCmbPropId 
		    }
		    if { [winfo exists $fLblPropId] } { 
			grid $fLblPropId 
		    }
		} else {
		    # Hide this property
		    if { [winfo exists $fCmbPropId] } { 
			grid remove $fCmbPropId 
		    }
		    if { [winfo exists $fLblPropId] } { 
			grid remove $fLblPropId 
		    }
		}
	    }
	}
    }
}

proc ::KMProps::cmbSectionTypeChange { f fullname {idTemplate ""} } {
    
    global KPriv
    
    # wa "SectionType f:$f fullname:$fullname idTemplate:$idTemplate"
    # Check that exists this combobox variable
    if {[info exists ::KMProps::cmbSectionType] } {
	
	set dv [::xmlutils::getComboDv $f $fullname "id" $idTemplate]
	# wa "dv:$dv"
	set xpath "Kratos_KWords/ElementCLaws/Item\[@id='SectionTypes'\]"
	::xmlutils::getAttribute $KPriv(xmlDocKKW) $xpath dv $dv
	
	# We update node ProfileDB values based on the new dv
	set fMat [string map {"SectionType" "ProfileDB"} $f]
	# wa "fMat:$fMat"
	if { $idTemplate == "" } {
	    
	    set fullnameMat [string map {"SectionType" "ProfileDB"} $fullname]
	    ::KMProps::changeCmbValues $fMat "$fullnameMat" $idTemplate "NoSearchFullname"
	} else {
	    
	    ::KMProps::changeCmbValues $fMat "MainProperties//ProfileDB" $idTemplate
	}
	
	# Get the cross section property list
	set PropertyList [::KMProps::GetCrossSectionPropertyList]
	# For each property
	 foreach propid $PropertyList {
	     set fCmbPropId [string map [list "SectionType" "$propid"] $f]
	     set fLblPropId [string map [list "cmb${propid}" "lbl${propid}"] $fCmbPropId]
	     # wa "fCmbPropId:$fCmbPropId fLblPropId:$fLblPropId"

	     # Get the show property state for the current property identifier => Using the section type
	     set ShowProperty [::KMProps::ShowPropertyBySectionType $propid]
	     # wa "propid:$propid ShowProperty:$ShowProperty"
	     if {$ShowProperty} {
	 	# Show this property	
	 	if { [winfo exists $fCmbPropId] } { 
	 	    grid $fCmbPropId 
	 	}
	 	if { [winfo exists $fLblPropId] } { 
	 	    grid $fLblPropId 
	 	}
	     } else {
	 	# Hide this property
	 	if { [winfo exists $fCmbPropId] } { 
	 	    grid remove $fCmbPropId 
	 	}
	 	if { [winfo exists $fLblPropId] } { 
	 	    grid remove $fLblPropId 
	 	}
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
    
    set fullname [string map [list "Activation//i.[list $id]" "Values//i.[list $whatDisable]"] $fullname]
    
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
    foreach Item $comboList {
	set strlen [string length $Item] 
	if {$strlen  > 5 } {
	    set width 10
	}
	if {$strlen > 10 } {
	    set width 15
	}
	if {$strlen > 15 } {
	    set width 20
	}
	if {$strlen > 20 } {
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
    
    # WarnWin "cmbselectChange remove:$remove selectVal:$selectVal"
    set refresh 0
    set fullname [DecodeName [$T item tag names $item]]
    set idFull [string map { "." "" "//" ""} $fullname]
    
    set id [::xmlutils::setXml $fullname id]
    
    set comboState [::xmlutils::getComboBoxState $fullname]
    # wa "comboState:$comboState fullname:$fullname idFull:$idFull id:$id"
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
    # wa "selCombo:$selCombo selComboText:$selComboText"
    if {$id == "Ax" || $id == "Ay" || $id == "Az"} {
    #if {0} {}
	set whatDisable [string map {A V} $id]
	set abuelo [$T item parent [$T item parent $item]]
	foreach itemChild [$T item descendants $abuelo] {
	    
	    set fullnameChild [DecodeName [$T item tag names $itemChild]]
	    set idChild [::xmlutils::setXml $fullnameChild id]
	    if {$idChild == $whatDisable} {
		if {$selCombo == 0} {
		    $T item enabled $itemChild 0
		    #$T item element configure $itemChild C0 elemTxtRead -fill { gray }
		    $T item element configure $itemChild C0 elemTxtRead -fill {black}
		    #$T item element configure $item C0 elemTxtRead -text "[$T item text $itemChild 0]: Disabled" -fill { gray }
		    #::xmlutils::setXml $fullnameChild state "write" disabled
		    ::xmlutils::setXml $fullnameChild state "write" normal
		} else {
		    $T item enabled $itemChild 1
		    set dvChild [::xmlutils::setXml $fullnameChild dv]
		    $T item element configure $itemChild C0 elemTxtRead -fill {black}
		    ::xmlutils::setXml $fullnameChild state "write" normal
		}
	    }
	}
    } elseif {$id == "ElemType"} {
	
	set matFullname [string map {"ElemType" "MatModel"} $fullname]
	
	::xmlutils::setXml $matFullname dv
	
	##Si se cambia el ElementType hay que cambiar el dv de Material Model
	set xpath "Kratos_KWords/ElementCLaws/Item\[@id='ElementTypes$::KMProps::nDim'\]"
	#set xpath "Kratos_KWords/ElementCLaws/Item\[@id='ElementTypes'\]"
	::xmlutils::getAttribute $::KPriv(xmlDocKKW) $xpath dv "$selCombo"
	set values [::xmlutils::getXMLValues "$matFullname" "" "" "" "NoElementFilter"]
	
	set dvMat [::xmlutils::setXml $matFullname dv]
	
	if {$dvMat ni $values} {
	    ::xmlutils::setXml $matFullname dv "write" [lindex $values 0]
	    #msg "$$matFullname \n $selCombo --> $values dvMat:$dvMat dv:[::xmlutils::setXml $matFullname dv] values0:[lindex $values 0]"
	    
	}
	#Hacemos un refresh para refrescar los cambios, tanto en matModel com en Thickness
	set refresh 1
	
    } elseif {$id eq "SectionType"} {
	
	set matFullname [string map {"SectionType" "ProfileDB"} $fullname]
	# wa "matFullname:$matFullname"
	::xmlutils::setXml $matFullname dv
	
	# When change the section type => Change the profile database and others properties
	set xpath "Kratos_KWords/ElementCLaws/Item\[@id='SectionTypes'\]"
	::xmlutils::getAttribute $::KPriv(xmlDocKKW) $xpath dv "$selCombo"
	set values [::xmlutils::getXMLValues "$matFullname" "" "" "" "NoElementFilter"]
	# wa "values:$values"
	set dvMat [::xmlutils::setXml $matFullname dv]
	
	if {$dvMat ni $values} {
	    ::xmlutils::setXml $matFullname dv "write" [lindex $values 0]
	}

	# Update the changes
	set refresh 1
	
    }
    
    set pid [::xmlutils::setXml $fullname pid]
    
    #msg "id:$id rem:$remove sel:$selCombo"
    if {$remove} {
	
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
    if {$imagen != -1} {
	$T item image $item C0 $imagen
    }
    
    #Miramos si ha cambiado para actualizar en consecuencia
    if {[::xmlutils::setXml $fullname dv] != $selCombo} {
	
	#Guarda el nuevo valor en el xml
	::xmlutils::setXml $fullname dv "write" $selCombo
	
	#Si el combo tiene un class especial habrá q reconstruir el árbol
	set clase [::xmlutils::setXml $fullname class]
	if {$clase != ""} {
	    
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
# se prepara todo para que el refresh tree cambie el treeCtl
#
proc ::KMProps::specialComboAction { T clase selCombo item id } {
    
    global KPriv 

    # Variable used to disabled the FSI
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
		 #wa "aplicId:$aplicId id.$id"
		if { $aplicId != $id } {
		    # Check the spacial case of DEM
		    if {$id eq "DEM"} {
			# Check for the case of fluid
			#msg "id $id appid $aplicId"
			if {$aplicId ne "Fluid"} {
			    ::xmlutils::setXml "${fullParent}//i.[list ${aplicId}]" dv "write" "No"
			    ::xmlutils::setXml "$aplicId" state "write" "hiddenAll"
			} else {
			    
			}
		    } elseif {$id eq "Fluid"} {
			# Check for the case of DEM
			if {$aplicId ne "DEM"} {
			    ::xmlutils::setXml "${fullParent}//i.[list ${aplicId}]" dv "write" "No"
			    ::xmlutils::setXml "$aplicId" state "write" "hiddenAll"
			} else {
			    # Seleccion Fluid / Si DEM -> show
			    ::KMProps::DemFluidSpecialInteraction "Yes"
			}
		    } else {
			::xmlutils::setXml "${fullParent}//i.[list ${aplicId}]" dv "write" "No"
			::xmlutils::setXml "$aplicId" state "write" "hiddenAll"
		    }
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
		    
		    ::xmlutils::setXml "ConvectionDiffusion" state "write" "hiddenAll"
		    ::xmlutils::setXml "ConvectionDiffusion" open "write" 0

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
	    # Se pulsa  NO
	    if { $id == "FluidStructureInteraction" } {
		if {$disabledFSI} {
		    ::xmlutils::setXml "FluidStructureInteraction" state "write" "hiddenAll"
		    ::xmlutils::setXml "FluidStructureInteraction" open "write" 0
		} else {
		    #if{ [::xmlutils::setXml "[list ${fullParent}]//i.[list ${aplicId}]" dv] } {}
		    ::xmlutils::setXml "StructuralAnalysis" state "write" "hiddenAll"
		    ::xmlutils::setXml "Fluid" state "write" "hiddenAll"
		    ::xmlutils::setXml "PFEM" state "write" "hiddenAll"
		    ::xmlutils::setXml "ConvectionDiffusion" state "write" "hiddenAll"
		}
	    } else {
		::xmlutils::setXml "$id" state "write" "hiddenAll"
	    }
	    if { $id == "Fluid" } {
		
		::KMProps::DemFluidSpecialInteraction "No"
		
	    }
	}
    } elseif {$clase eq "DEM-Rolling"} {
	::KMProps::DemMaterialSpecialInteraction $id
	
    }
    
    # WarnWin "selCombo:$selCombo"
    foreach var $::KMProps::visibilityVars {
	# wa "var $var"
        if {$var == $clase} {
            
            # General case
            if { [set ::KMProps::$var] != "$selCombo" } {
                set ::KMProps::$var $selCombo
                #msg "var$var"
            }
            
            # Special cases
            # Compressible Fluids
            if { $var == "fluidType" } {
                if {$selCombo == "Compressible" } {
                    set txt [= "Compressible fluids are not available in this version"]   
                    WarnWin "$txt."
                    set fullname [DecodeName [$T item tag names $item]]
                    ::xmlutils::setXml "$fullname" dv "write" "Incompressible"
                } 
            } elseif { $var == "strucType" } {
		if { $selCombo == "Beam" } {
		    # Check the special case of kinematic type
		    set xpath "StructuralAnalysis//c.AnalysisData//i.KinematicType"
		    set KinematicType [::xmlutils::setXml "$xpath" dv]
		    # wa "KinematicType:$KinematicType"
		    if { $KinematicType == "LargeDisplacements" } {
			set txt [= "Beam structural type can not be used with the large displacements formulation in this version"]  
			WarnWin "$txt."
			::xmlutils::setXml "$xpath" dv "write" "SmallDisplacements"
		    }
                } elseif { $selCombo == "Shell" } {
		    # Check the special case of kinematic type
		    set xpath "StructuralAnalysis//c.AnalysisData//i.KinematicType"
		    set KinematicType [::xmlutils::setXml "$xpath" dv]
		    # wa "KinematicType:$KinematicType"
		    if { $KinematicType == "LargeDisplacements" } {
			set txt [= "Shell structural type can not be used with the large displacements formulation in this version"]  
			WarnWin "$txt."
			::xmlutils::setXml "$xpath" dv "write" "SmallDisplacements"
		    }
                }
	    } elseif { $var == "kinemType" } {
		if { $selCombo == "LargeDisplacements" } {
		    # Check the special case of kinematic type
		    set xpath "StructuralAnalysis//c.AnalysisData//i.StructuralType"
		    set StructuralType [::xmlutils::setXml "$xpath" dv]
		    # wa "StructuralType:$StructuralType"
		    if { $StructuralType == "Beam" } {
			set txt [= "The large displacements formulation can not be used with the beam structural type in this version"]  
			WarnWin "$txt."
			set xpath "StructuralAnalysis//c.AnalysisData//i.KinematicType"
			::xmlutils::setXml "$xpath" dv "write" "SmallDisplacements"
		    } elseif { $StructuralType == "Shell" } {
			set txt [= "The large displacements formulation can not be used with the shell structural type in this version"]  
			WarnWin "$txt."
			set xpath "StructuralAnalysis//c.AnalysisData//i.KinematicType"
			::xmlutils::setXml "$xpath" dv "write" "SmallDisplacements"
		    }
                }
	    } 
        }
    }
    return ""
}

# Sirve para modificar el arbol de materiales en funcion de combos especiales de DEM
proc ::KMProps::DemMaterialSpecialInteraction { id } {
   
    if {$id eq "DEM-RollingFriction"} {
	set gproplist [::xmlutils::setXmlContainerIds "DEMMaterial" "Material" "mats"]
	foreach demMaterial $gproplist {
	    set path "DEM//c.DEM-Options//c.DEM-AdvancedOptions//i.DEM-RollingFriction"
	    set rollstate [::xmlutils::setXml $path "dv"]
	    # Caso especial para ocultar propiedad de materiales RollingFriction
	    if {$rollstate eq "Yes"} {
		::xmlutils::setXml "DEMMaterial//m.$demMaterial//p.RollingFriction" state "write" "normal" "mat"
		::xmlutils::setXml "DEMMaterial//m.$demMaterial//p.RollingFrictionWithWalls" state "write" "normal" "mat"
	    } else {
		::xmlutils::setXml "DEMMaterial//m.$demMaterial//p.RollingFriction" state "write" "hidden" "mat"
		::xmlutils::setXml "DEMMaterial//m.$demMaterial//p.RollingFrictionWithWalls" state "write" "hidden" "mat"
	    }  
	}
    } elseif {$id eq "Fluid"} {
	set gproplist [::xmlutils::setXmlContainerIds "DEMMaterial" "Material" "mats"]
	foreach demMaterial $gproplist {
	    
	    set path "GeneralApplicationData//c.ApplicationTypes//i.Fluid"
	    set Fluidstate [::xmlutils::setXml $path "dv"]
	    if {$Fluidstate eq "Yes"} {
		# Caso especial para ocultar propiedad de materiales ParticleSphericity
		::xmlutils::setXml "DEMMaterial//m.$demMaterial//p.ParticleSphericity" state "write" "normal" "mat"
		::xmlutils::setXml "DEM//c.DEM-Fluid-interaction" state "write" "normal"
	    } else {
		# Caso especial para ocultar propiedad de materiales ParticleSphericity
		::xmlutils::setXml "DEMMaterial//m.$demMaterial//p.ParticleSphericity" state "write" "hidden" "mat"
		::xmlutils::setXml "DEM//c.DEM-Fluid-interaction" state "write" "hiddenAll"
	    }
	}  
    }
    ::KMat::refreshTree
}

proc ::KMProps::DemFluidSpecialInteraction { Fluidstate } {
    #msg "$Fluidstate"
#    if {$Fluidstate eq "No"} {
#	::xmlutils::setXml "DEM//c.DEM-Fluid-interaction" state "write" "hiddenAll"
#	::xmlutils::setXml "DEM//c.DEM-Elements//c.SwimmingDEMElement" state "write" "hiddenAll"
#	::xmlutils::setXml "DEM//c.DEM-Elements//c.SphericContPartDEMElement3D" state "write" "normal"
#	::xmlutils::setXml "DEM//c.DEM-Elements//c.SphericPartDEMElement3D" state "write" "normal"
#	::xmlutils::setXml "DEM//c.DEM-Elements//c.CilContPartDEMElement2D" state "write" "normal"
#	::xmlutils::setXml "DEM//c.DEM-Elements//c.CilPartDEMElement2D" state "write" "normal"
#	} else {
#	::xmlutils::setXml "DEM//c.DEM-Fluid-interaction" state "write" "normal"
#	::xmlutils::setXml "DEM//c.DEM-Elements//c.SwimmingDEMElement" state "write" "normal"
#	::xmlutils::setXml "DEM//c.DEM-Elements//c.SphericContPartDEMElement3D" state "write" "hiddenAll"
#	::xmlutils::setXml "DEM//c.DEM-Elements//c.SphericPartDEMElement3D" state "write" "hiddenAll"
#	::xmlutils::setXml "DEM//c.DEM-Elements//c.CilContPartDEMElement2D" state "write" "hiddenAll"
#	::xmlutils::setXml "DEM//c.DEM-Elements//c.CilPartDEMElement2D" state "write" "hiddenAll"
#    }
    ::KMProps::DemMaterialSpecialInteraction "Fluid"
}

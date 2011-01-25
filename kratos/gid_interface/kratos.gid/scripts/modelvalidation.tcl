###############################################################################
#
#	NAME: modelvalidation.tcl
#
#	PURPOSE: Validation of interface data
#
#	QUANTECH ATZ-DEVELOPMENT DEPARTMENT
#
#	AUTHOR : G. Socorro => GSM
#
#	CREATED AT: 09/06/2010
#
#	HISTORY:
#
#       0.5- 25/01/11-GSM, show the warning/error message only when find some error
#       0.4- 01/20/10 LC, Corregido bug con la gestión de initalConditions al mostrar errores
#	0.3- 27/09/10 LC, Se ha pasado InitialConditions de errores a warnings, 
#		se valida que los elements, conditions e InitialConditions tengan algún grupo activo, 
#		y se valida que el path de Kratos exista en la computadora.
#	0.2- 11/06/10 GS, Update CreateReportWindow to use InitWindow and add a close button
#	0.1- 25/02/10 KS, create a base source code
#
###############################################################################

package provide KMat 1.0 


# Create a base namespace KMat
namespace eval ::KMValid:: {
	
	# Window path
	variable winpath
	variable Errors 0
	variable Warnings 0
	variable initialConditions ""
}

proc ::KMValid::Init {} {
	
	variable winpath

	set winpath ".gid.modelvalidation"
	
}

proc ::KMValid::ValidateModel {{w .gid.modelvalidation}} {
	
	# Init some properties
	::KMValid::Init

	# create report window
	set result [ ::KMValid::CreateReportWindow $w]
	
	return $result
	
}

proc ::KMValid::CreateReportWindow {w} {
	
	global KPriv

	#Comprobamos en los settings del proyecto si es necesario validar la versión
	set ValidateModel [::kps::getConfigValue "ValidateModel"]
	
	if { !$ValidateModel } {
		
		return 0
	}
	
	# Init the window
	set title [= "Model validation"]
	InitWindow $w $title ::KMValid::CreateReportWindowWindowGeom ::KMValid::CreateReportWindow

	# Text with scrolls only if required
	frame $w.fr
	scrollbar $w.fr.scrolly -command "$w.fr.t yview" \
		-orient vertical -relief sunken
	scrollbar $w.fr.scrollx -command "$w.fr.t xview"\
		-orient horizontal -relief sunken
	text $w.fr.t -font SmallFont -yscrollcommand "$w.fr.scrolly set" \
		-xscrollcommand "$w.fr.scrollx set" -wrap none
	
	grid $w.fr.t -row 1 -column 1 -sticky nsew
	grid $w.fr.scrolly -row 1 -column 2 -sticky ns
	grid $w.fr.scrollx -row 2 -column 1 -sticky ew
	grid rowconf $w.fr 1 -weight 1
	grid columnconf $w.fr 1 -weight 1
	grid $w.fr -sticky nsew -padx 5 -pady 5
	
	grid remove $w.fr.scrolly
	grid remove $w.fr.scrollx
	bind $w.fr.t <Configure> "::WinUtils::ConfigureListScrollbars $w.fr.t $w.fr.scrollx $w.fr.scrolly"
	
	# Lower buttons
	set def_back [$w cget -background]
	frame $w.frmButtons -bg [CCColorActivo $def_back]
	
	# Close button
	Button $w.frmButtons.btnclose -text [= "Close"] \
		-command "::KMValid::CreateReportWindowbClose $w" \
		-helptext [= "Close the report window"] \
		-underline 0 -width 6
	
	# Print button
	Button $w.frmButtons.btnprint -text [= "Print"] \
		-command "::WinUtils::Print $w.fr.t" \
		-helptext [= "Print the report via a temporal file using notepad"] \
		-underline 0 -width 6
	
	SetWidgetsWidthFromText $w.frmButtons.btnclose $w.frmButtons.btnprint
	
	# Geometry manager
	grid $w.frmButtons \
		-sticky ews \
		-columnspan 7
	
	grid $w.frmButtons.btnprint \
		-in $w.frmButtons \
		-row 0 -column 0 \
		-sticky wn \
		-padx 4m -pady 2m
	
	grid $w.frmButtons.btnclose \
		-in $w.frmButtons \
		-row 0 -column 7 \
		-sticky en \
		-padx 4m -pady 2m

	
	# For w
	wm protocol $w WM_DELETE_WINDOW "::KMValid::CreateReportWindowbClose $w"

	grid rowconfigure $w 0 -weight 1
	grid columnconf $w 0 -weight 1
	
	# For button
	grid rowconfigure $w.frmButtons 0 -weight 1
	grid columnconf $w.frmButtons {0 1 2 3 4 5 6 7} -weight 1
	
	# Set the focus to the close button   
	focus $w.frmButtons.btnclose
	
	# Binding
	bind $w <Alt-c> "tkButtonInvoke $w.frmButtons.btnclose"
	bind $w <Escape> "tkButtonInvoke $w.frmButtons.btnclose"
	
	# fill validation report list
	set allreportlist [list]
	
	set allreportlist [::KMValid::ValidateProjectConfiguration $allreportlist]
	
	set allreportlist [::KMValid::ValidateGroups $allreportlist]
	
	#lappend allreportlist "ErrorAplication:ApplicationTitle Structural Analysis:"
	
	set xml $KPriv(xml)	
	set apliXpath "/Kratos_Data/RootData\[@id='GeneralApplicationData'\]/Container\[@id='ApplicationTypes'\]/Item"
	
	
	set FSI [::xmlutils::getAttribute $xml "${apliXpath}\[@id='FluidStructureInteraction'\]" dv]
	
	if { $FSI == "Yes" } {
		set structural "Yes"
		set fluid "Yes"
	} else {
		set structural [::xmlutils::getAttribute $xml "${apliXpath}\[@id='StructuralAnalysis'\]" dv]	
		set fluid [::xmlutils::getAttribute $xml "${apliXpath}\[@id='Fluid'\]" dv]
	}
	
	if { $structural == "Yes" } {
		
		set reportAux {}
		set reportAux [::KMValid::ValidateAssignedGroupsInModel $reportAux StructuralAnalysis]
		
		if {[llength $reportAux]} {
			lappend allreportlist "Error: STRUCTURAL ANALYSIS APPLICATION:\n"
			lappend allreportlist "Error:\n"
			set allreportlist [concat $allreportlist $reportAux]
		}		
	}
	
	if {$fluid == "Yes" } {
		
		set reportAux {}
		set reportAux [::KMValid::ValidateAssignedGroupsInModel $reportAux Fluid]
		if {[llength $reportAux]} {
			if { $structural == "Yes" } {
				lappend allreportlist "Error:\n"
			}
			lappend allreportlist "Error: FLUID  APPLICATION:"
			lappend allreportlist "Error:\n"
			set allreportlist [concat $allreportlist $reportAux]
		}
		
		#Se comprueba a parte si hay que añadir el Warnning de InitialConditions
		if { $::KMValid::initialConditions != "" } {
			
			lappend allreportlist "$::KMValid::initialConditions"
			lappend allreportlist "Warning:"
		}
	}
	
	set allreportlist [::KMValid::ValidateProjectInformation $allreportlist]
	
	# process messages
	::KMValid::ProcessMessages $w $allreportlist
		
	set aviso ""
	if { $::KMValid::Errors } {
	    set aviso "There are Errors in the model, do you want to continue anyway?"
	} elseif { $::KMValid::Warnings } {
	   # set aviso "There are warnings in the model, do you want to continue anyway?"
	} else {
	    destroy $w
	}
	
	if { $aviso != "" } {
		set answer [::WinUtils::confirmBox "." "$aviso" "okcancel"]
		if { $answer == "ok" } {
			return 1
		} else {
			return "-cancel-"
		}
	}
}

proc ::KMValid::CreateReportWindowbClose {{w .gid.modelvalidation}} {

	# Destroy the window widget	
	if {[winfo exists $w]} { 
		destroy $w
	}
}

proc ::KMValid::ProcessMessages {w allreportlist} {
	
	# Configure the result tag
	set tpath $w.fr.t
	$tpath tag configure bold -font {Courier 11 bold italic}
	$tpath tag configure big -font {Courier 13 bold}
	$tpath tag configure bigblue -font {Courier 13 bold} -foreground blue
	$tpath tag configure bigred -font {Courier 13 bold} -foreground red
	$tpath tag configure bigblack -font {Courier 11 bold} -foreground black
	$tpath tag configure verybig -font {Helvetica 24 bold}
	$tpath tag configure margins1 -lmargin1 10m -lmargin2 6m -rmargin 10m -underline on
	$tpath tag configure margins1a -lmargin1 14m -lmargin2 6m -rmargin 10m 
	$tpath tag configure margins2 -lmargin1 14m -lmargin2 6m -rmargin 10m
	$tpath tag configure redcolor -foreground red
	$tpath tag configure bluecolor -foreground blue
	$tpath tag configure bluemargins2 -foreground blue -lmargin1 14m -lmargin2 6m -rmargin 10m

	if {[llength $allreportlist]>0} {
		
		set nrlinelist [list]
		set errorlist [list]
		set warninglist [list]
		# Add all the text with your format
		foreach items $allreportlist {
			set keyword [lindex $items 0]
			# WarnWinText "keyword:$keyword"
			switch $keyword {
				"NRLine:" {
					set cstr [lrange $items 1 end]
					lappend nrlinelist "${cstr}"
				}
				"Error:" {
					lappend errorlist [lrange $items 1 end]
					
				}
				"Warning:" {
					lappend warninglist [lrange $items 1 end]
					
				}
			}
		}

		set bcount -1
		# Error
		#set applications [list "Structural Annalysis" "Fluid" ]
		
		if {[llength $errorlist]>0} {
			set ::KMValid::Errors 1
			$tpath insert end "***** Begin error *****\n\n"  bigred
			
			foreach items $errorlist {
				
				#Opciones especiales para resaltar la aplicación y solo imprimirla si tiene errores
				if { $items == "STRUCTURAL ANALYSIS APPLICATION:" } {
					
					#Comprobamos que haya algún error
					$tpath insert end "${items}\n" bigblack
					
				} elseif { $items == "FLUID APPLICATION:" } {
						
						$tpath insert end "${items}\n" bigblack
						
				#Caso normal: imprimir los errores
				} else {
					# Check the end key
					set bkey [lindex $items end]
					if {$bkey =="BPSetting"} {
						incr bcount 1
						set items [lrange $items 0 end-1]
						::ImportBPM::ProcessTags $tpath $bcount $items $bkey
					} elseif {$bkey =="BPStation"} {
						incr bcount 1
						set items [lrange $items 0 end-1]
						::ImportBPM::ProcessTags $tpath $bcount $items $bkey
					} else {

						$tpath insert end "${items}\n"
					}
				}
			}
			$tpath insert end "\n***** End error *****\n\n"  bigred
		}

		# Warning
		if {[llength $warninglist]>0} {
			
			set ::KMValid::Warnings 1
			
			$tpath insert end "***** Begin warning *****\n\n"  bigblue
			foreach items $warninglist {
				# Check the end key
				set bkey [lindex $items end]
				if {$bkey =="BPSetting"} {
					incr bcount 1
					set items [lrange $items 0 end-1]
					::ImportBPM::ProcessTags $tpath $bcount $items $bkey
				} elseif {$bkey =="BPStation"} {
					incr bcount 1
					set items [lrange $items 0 end-1]
					::ImportBPM::ProcessTags $tpath $bcount $items $bkey
				} else {
					$tpath insert end "${items}\n" 
				}
			}
			$tpath insert end "\n***** End warning *****\n\n"  bigblue
		}

		# Not decoded lines
		if {[llength $nrlinelist]>0} {
			$tpath insert end "***** Begin lines without decoding *****\n\n" big
			foreach items $nrlinelist {
				# WarnWinText "antes items:$items"
				set items [::ImportBPM::ChangeDollarString $items]
				# WarnWinText "items:$items"
				$tpath insert end "${items}\n" 
			}
			$tpath insert end "\n***** End lines without decoding *****\n\n" big
		}
	}
	
	$w.fr.t see end
	update

	::WinUtils::ConfigureListScrollbars $w.fr.t $w.fr.scrollx $w.fr.scrolly
	$w.fr.t configure -state disabled
	
	# Activate bindings
	wm protocol $w WM_DELETE_WINDOW "destroy $w" 
	
	update
}

#
# Error: Como mínimo hay un grupo
# Warning: si hay grupos vacíos (sin entidades)
# 
proc ::KMValid::ValidateGroups { allreportlist } {
		
	global KPriv
	set xml $KPriv(xml)
	
	set xpath "/Kratos_Data/Groups/Group"
	
	set groups [::xmlutils::getXmlChildIds $xml $xpath ] 
	if { $groups == 0 } {
		
		lappend allreportlist "Error: There are no groups defined."
		lappend allreportlist "Error: \n"
		
		set line1 ""
        set line2 ""
        lappend allreportlist ${line1}
        lappend allreportlist ${line2}
        
	} else {
		
		set noEntitiesGroups [list ]
		#set numnoentity 0
		foreach grup $groups {
			
			if {![::KEGroups::getGroupGiDEntities $grup ALL hasEntities]} {
				
				lappend noEntitiesGroups $grup
				#incr numNoEntity
			}
		}
		
		if { [llength $noEntitiesGroups] } {
		
			if { [llength $noEntitiesGroups] == [llength $groups] } {
				
				lappend allreportlist "Error: There are no groups with entities defined."
				lappend allreportlist "Error: \n"
				
			} else  {
	
				lappend allreportlist "Warning: There are groups with no entities: "
				foreach group $noEntitiesGroups {
					lappend allreportlist "Warning: $group"
				}
				lappend allreportlist "Warning:"
			}
			set line1 ""
	        set line2 ""
	        lappend allreportlist ${line1}
	        lappend allreportlist ${line2}
		}
	}

	return $allreportlist
}


#
# Error: Por lo menos 1 PROPERTY (con id)
# Error: Que la property tenga un material ID
# Error: Que las propiedades activas no sean nulas
# Error: Que la property tenga Thickness no nulo (donde sea necesario)
#
# Error: Como mínimo hay un grupo activo con entidades en ELEMENTS
# Error: Como mínimo hay un grupo activo con entidades en LOADS 
# Error: BODYFORCE y PUNTUAL (Que las propiedades no sean nulas)
# Error: Como mínimo hay un grupo con entidades en CONDITIONS
# 
proc ::KMValid::ValidateAssignedGroupsInModel { allreportlist {application "StructuralAnalysis"} } {
	
	set appPath "/Kratos_Data/RootData\[@id='$application'\]"
	
	#Avisaremos de todos los errores al final
	set elementsGroups 0
	set conditionsGroups 0
	
	set hasMaterialId 0
	set hasThickness 0
	set noIdProps 0
	set noMaterialProps {}
	set nullThickness {}
	
	## StructuralAnalysis
	global KPriv
	set xml $KPriv(xml)
	
	#ELEMENTS

	#Buscamos los grupos asignados a algún element
	set xpath "$appPath/Container\[@id='Elements'\]/Container"
	set nodes [$xml selectNodes $xpath]
	
	foreach node $nodes {
		
		set idElem [$node getAttribute id ""]
		
		set groupsXPath "${xpath}\[@id='$idElem'\]/Container"
		set groups [::xmlutils::getXmlChildIds $xml "${groupsXPath}" ] 
		#msg "idElem$idElem   groups:$groups\n"
		foreach group $groups {
			
			set grNodeXPath "${groupsXPath}\[@id='$group'\]"
			set active [::xmlutils::getAttribute $xml $grNodeXPath active]
			
			if {$active} {
			
			
				if {[::KEGroups::getGroupGiDEntities $group ALL hasEntities]} {
					#Si uno de los grupos tiene entidades ya no sacaremos ese error
					
					#Pero además tiene que estar activo para esa convinación de filtros (visible en el árbol)
					#set isActiveGroup $groupsXPath...
					
					
					set elementsGroups 1
				}
				
				#Ahora comprobamos la propiedad asignada a este grupo
				set grNodeXPath "${groupsXPath}\[@id='$group'\]/Container"
				set grNodeConainers [$xml selectNodes "$grNodeXPath"]
				foreach nodeContainer $grNodeConainers {
					
					set nodeItems [$nodeContainer childNodes]
					foreach item $nodeItems {
						
						set idItem [$item getAttribute id ""]
						if { $idItem == "Property" } {
							set propId [$item getAttribute dv ""]
							
							#Ahora buscamos esa propiedad en 'Properties'
							set xpath "$appPath/Container\[@id='Properties'\]/Container\[@id='$propId'\]"
							set node [$xml selectNodes $xpath]
							if { $node != "" } {
								
								foreach nodeCont [$node childNodes] {
									foreach nodePropIt [$nodeCont childNodes] {
										
										#Aquí tenemos los items de una propiedad
										set idItem [$nodePropIt getAttribute id ""]
										#Comprobamos que el material no sea nulo
										if { $idItem == "Material" } {
											if {  [$nodePropIt getAttribute dv ""] == "" } {
												
												#Si la propiedad no tiene material asignado, 
												#la añadimos a la lista de propiedades sin material
												if { !($propId in $noMaterialProps) } {
												
													lappend noMaterialProps $propId
												}
											}
										} elseif { $idItem == "ElemType" } {
											
											set ::KMProps::ElemTypeThickness [$nodePropIt getAttribute dv ""]
											
										} elseif { $idItem == "Thickness" } {
											
											#Comprobamos si necesita Thicness, y en ese caso que no sea nulo 
											#set ::KMProps::ElemTypeThickness $idElem
											set xpath "/Kratos_Data/RootData\[@id='GeneralApplicationData'\]/Container\[@id='Domain'\]/Item\[@id='SpatialDimension'\]"
											set ::KMProps::nDim  [::xmlutils::getAttribute $xml $xpath dv]
											if { [::KMProps::showThickness] && [$nodePropIt getAttribute dv ""] == "" } {
												
												#Si la propiedad necesita tickness y no lo tiene, error
												if { !($propId in $nullThickness) } {
												
													lappend nullThickness $propId
												}
											}
										}
									}
								}
							}
						}	
					}
				}
			}
		}
	}
	
	if { ! $elementsGroups } {
		
		lappend allreportlist "Error: There are no groups assigned to Elements."
		lappend allreportlist "\n"
	}
	
	#
	# PROPERTIES
	#
	set propsMessages {}
	
	set xpath "$appPath/Container\[@id='Properties'\]/Container"
	set nodes [$xml selectNodes $xpath]
	if { [llength $nodes] == 0 } {
		
		lappend propsMessages "Error: There are no properties defined."
	} else {
		foreach node $nodes {
			if { [$node getAttribute id ""] == "" } { 
				set noIdProps 1
			}
		}
		if { $noIdProps } {
			lappend propsMessages "Error: There are properties without a correct ID."
		}
	}
	
	if { [llength $noMaterialProps] } {
		
		lappend propsMessages "Error: There are properties without an associated material:\n$noMaterialProps"
	}
	if { [llength $nullThickness] } {
		
		lappend propsMessages "Error: There are properties without a correct Thickness:\n$nullThickness"
	}
	
	if { [llength $propsMessages] } {
		
		#lappend allreportlist "Error:PROPERTIES:"	
		foreach prop $propsMessages {
	
			lappend allreportlist $prop
		}
	}
	
	set allreportlist [::KMValid::checkGroups $xml $appPath "Conditions" $allreportlist]
	
	if { $application == "StructuralAnalysis" } {
		
		set allreportlist [::KMValid::checkGroups $xml $appPath "Loads" $allreportlist]
		
	} elseif { $application == "Fluid" } {
		
		#Ahora ya no es un error, es un warning, por lo que se gestiona de forma distinta
		set allreportlist [::KMValid::checkGroups $xml $appPath "InitialConditions" $allreportlist]
	}
	
	#set xpath "$appPath/Container\[@id='Conditions'\]/Container"
	#if { [llength [::KMValid::groupsWithEntities $xml $xpath]] == 0 } {	
	#	lappend allreportlist "Error: There are no groups assigned to Conditions."
	#}
	
	#set line1 ""
	#set line2 ""
	#lappend allreportlist ${line1}
	#lappend allreportlist ${line2}

	return $allreportlist
}

#
# SEPARAR POR APLICACIONES (con color)
# Poner el path en los errores
#
proc ::KMValid::groupsWithEntities { xml xpath {returnNodes 0}} {
	
	set groupEntities {}
	
	#Buscamos los grupos asignados a algún element
	set nodes [$xml selectNodes $xpath]
	
	foreach node $nodes {
		
		set idElem [$node getAttribute id ""]
		
		set groupsXPath "${xpath}\[@id='$idElem'\]/Container"
		set groups [$xml selectNodes $groupsXPath]
		#set groups [::xmlutils::getXmlChildIds $xml "${groupsXPath}" ]
		
		foreach nodeGroup $groups {
			
			set idGroup [$nodeGroup getAttribute id ""]
			
			#Se tiene que cumplir que el nodo esté activo (visible en árbol según los filtros seleccionados)
			if {[$nodeGroup getAttribute active 0] != 0} {
				
				if {[::KEGroups::getGroupGiDEntities $idGroup ALL hasEntities]} {
					#Si uno de los grupos tiene entidades ya no sacaremos ese error
					if { $returnNodes } {
						lappend groupEntities $nodeGroup
					} else {
						lappend groupEntities $idGroup
					}
				}
			}
		}
	}
	return $groupEntities
}

proc getPropertyItems { } {

	set xpath "$appPath/Container\[@id='Properties'\]/Container"
	set nodes [$xml selectNodes $xpath]
	foreach nodeCont $nodes {
		foreach node [$nedeCont childNodes] {
			#Aquí tenemos los items de una propiedad
		}
	}
}

#
# Campos obligatorios con * error y los que no Warning
#
proc ::KMValid::ValidateProjectInformation { allreportlist } {
	
	global KPriv
	set xml $KPriv(xml)
	
	set importantFields ""
	
	set xpath "/Kratos_Data/RootData\[@id='GeneralApplicationData'\]/Container\[@id='ProjectInfo'\]/Item"
	set nodes [$xml selectNodes $xpath]
	foreach node $nodes {
		
		set style [$node getAttribute style ""]
		if { $style == "*" } {
			if { [$node getAttribute dv ""]	== "" } {
				set id [$node getAttribute id ""]
				lappend importantFields $id
			}		
		}
	}
	if { [llength $importantFields] } {
		lappend allreportlist "Warning: There are mandatory fields (*) with no info: "
		foreach field $importantFields {
			lappend allreportlist "Warning: $field"
		}
		lappend allreportlist "Warning:\n"
	}
	
	set xpath "/Kratos_Data/RootData\[@id='GeneralApplicationData'\]/Container\[@id='ProjectConfiguration'\]/Item\[@id='KratosPath']"
	set node [$xml selectNodes $xpath]
	
	if {$node != "" } {
		set dv [$node getAttribute dv ""]
		
		if { ![file isdirectory $dv] } {
			lappend allreportlist ""
			lappend allreportlist "Error: The kratos path '$dv' does not exists."
			lappend allreportlist "Error:"
		}
	}
	
	set line1 ""
	set line2 ""
	lappend allreportlist ${line1}
	lappend allreportlist ${line2}

	return $allreportlist
	
}

#
#
#
proc ::KMValid::ValidateProjectConfiguration { allreportlist } {

	set line1 ""
	set line2 ""
	lappend allreportlist ${line1}
	lappend allreportlist ${line2}

	return $allreportlist
}

#
# Accede a un nodo principal como LOADS o CONDITIONS y mira si tiene algún grupo asignado
#
proc ::KMValid::checkGroups { xml appPath rootContainer allreportlist} {
	
	set nullProps {}
	set xpath "$appPath/Container\[@id='$rootContainer'\]/Container"
	set returnNodes 1
	set assignedGroups [::KMValid::groupsWithEntities $xml $xpath $returnNodes]
	
	if { [llength $assignedGroups] } {
		
		foreach nodeGroup $assignedGroups {
			
			foreach nodeContainer [$nodeGroup childNodes] {
				
				set nodeItems [$nodeContainer childNodes]
				foreach item $nodeItems {
					set idItem [$item getAttribute id ""]
					set dv [$item getAttribute dv ""]
					if { $dv == "" } {
						lappend nullProps [$item getAttribute id ""]
					}
				}
			}
			if { [llength $nullProps] } {
				
				set parentPID [[$nodeGroup parentNode] getAttribute pid ""]
				set groupId [$nodeGroup getAttribute id ""]
				
				#if { $firstEmptyProp } {
					#set firstEmptyProp 0
					#lappend allreportlist "Error: LOADS:"
				#}
				
				lappend allreportlist "Error: $rootContainer-> $parentPID-> group '$groupId': Empty properties '$nullProps'"
				
				
			}
		}
		#if { ! $firstEmptyProp } {
			#lappend allreportlist "Error:\n"
		#}
	} else {
		
		#En el caso de Initial conditions se ha cambiado de error a warning, por lo que se ha creado
		# una variable global para tratar este caso particular en vez de cambiar el diseño
		if { $rootContainer == "InitialConditions" } {

			set ::KMValid::initialConditions "Warning: There are no groups assigned to $rootContainer."
		} else {
			
			lappend allreportlist "Error: There are no groups assigned to $rootContainer."
		}
	}
	
	return $allreportlist
}
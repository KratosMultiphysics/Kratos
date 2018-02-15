###############################################################################
#
#        NAME: modelvalidation.tcl
#
#        PURPOSE: Validation of interface data
#
#        QUANTECH ATZ-DEVELOPMENT DEPARTMENT
#
#        AUTHOR : G. Socorro => GSM
#
#        CREATED AT: 09/06/2010
#
#        HISTORY:
#
#   1.7- 12/07/13- GSM, correct a bug in the proc CheckRepeatedNodes (change write_calc_data by GiD_EntitiesGroups)
#   1.6- 17/06/13- GSM, modify the procs ValidateGroups, groupsWithEntities and ValidateAssignedGroupsInModel to use only the new GiD group,
#   1.5- 17/05/13- GSM, add the validation of the cross section properties for the structural analysis application
#   1.4- 26/11/12- J. Gárate, PFEM correction on ::KMValid::isWall, join 3D and 2D
#   1.3- 26/11/12- J. Gárate, PFEM support, ::KMValid::isWall
#   1.2- 22/10/12- J. Gárate, Validation function for drag files field. Support fot new GiD_Groups
#   1.1- 10/10/12- GSM, improve the model validation to element type and cross section properties
#   1.0- 03/10/12- GSM, correct a bug with the namespace variable Errors
#   0.9- 20/07/12- GSM, update some proc to delete old source code
#   0.8- 19/07/12- J. Gárate, Check if any group or node is shared between Slip and NoSlip
#   0.7- 09/02/12- J. Gárate, Deshabilitada la comprobacion del Kratos Path
#   0.6- 08/02/12- J. Gárate, Si no hay ni error ni warning,  cierra la ventana de Model Validation
#   0.5- 25/01/11-GSM, show the warning/error message only when find some error
#   0.4- 01/20/10 LC, Corregido bug con la gestión de initalConditions al mostrar errores
#   0.3- 27/09/10 LC, Se ha pasado InitialConditions de errores a warnings, 
#                     se valida que los elements, conditions e InitialConditions tengan algún grupo activo, 
#                     y se valida que el path de Kratos exista en la computadora.
#   0.2- 11/06/10 GS, Update CreateReportWindow to use InitWindow and add a close button
#   0.1- 25/02/10 KS, create a base source code
#
###############################################################################



# Create a base namespace KMat
namespace eval ::KMValid:: {
    
    # Window path
    variable winpath
    variable Errors 
    variable Warnings 
    variable initialConditions ""
}

proc ::KMValid::Init {} {
    
    variable winpath

    set winpath ".gid.modelvalidation"

    variable Errors 
    set Errors 0

    variable Warnings 
    set Warnings 0

    variable initialConditions ""
    set initialConditions ""
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

    # Comprobamos en los settings del proyecto si es necesario validar la versión
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
    
    # Get the time
    set systemTime [clock seconds]
    set monthNumber [string trimleft [clock format [clock seconds] -format %m] 0]
    set yearNumber [string trimleft [clock format [clock seconds] -format %Y] 0]
    #WarnWin $monthNumber
    #WarnWin $yearNumber
    
    # fill validation report list
    set allreportlist [list]
    
    #set limityear 2015
    #set limitmonth 7
    
    #if { $yearNumber == 2015 && $monthNumber == 6 } {
    #    lappend allreportlist "Warning: The License for this code will expire at the beginning of the next month.\n This is a code under development. Contact the CIMNE DEMPack Team for a new version of the code before it expires.\nSorry about the inconvenience.\n"
    #}
    #
    #if { ($yearNumber >= $limityear && $monthNumber >= $limitmonth) || $yearNumber > $limityear } {
    #    lappend allreportlist "Error: LICENSE EXPIRED!!\nThis was a code under development. Contact the CIMNE DEMPack Team for a new version of the code.\nSorry about the inconvenience."
        # process messages
    #    ::KMValid::ProcessMessages $w $allreportlist   
    #    return "-cancel-"
    #}                        
    
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
    set allreportlist [::KMValid::ValidateDragFiles $allreportlist]
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
	# Validate BC  Que no coincidan grupos con condiciones Slip y No Slip
	set msg [::KMValid::ValidateGroupRepeat "Is-Slip" "No-Slip"]
	if {[llength $msg]} {
	    lappend allreportlist $msg
	    lappend allreportlist "Warning:"
	}
    }
    
    set allreportlist [::KMValid::ValidateProjectInformation $allreportlist]
    
    # WarnWin "allreportlist:$allreportlist ::KMValid::Errors:$::KMValid::Errors"

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
    $tpath tag configure big -font {Courier 12 bold}
    $tpath tag configure bigblue -font {Courier 12 bold} -foreground blue
    $tpath tag configure bigred -font {Courier 12 bold} -foreground red
    $tpath tag configure bigblack -font {Courier 11 bold} -foreground black
    $tpath tag configure verybig -font {Helvetica 24 bold}
    $tpath tag configure margins1 -lmargin1 10m -lmargin2 6m -rmargin 10m -underline on
    $tpath tag configure margins1a -lmargin1 14m -lmargin2 6m -rmargin 10m 
    $tpath tag configure margins2 -lmargin1 14m -lmargin2 6m -rmargin 10m
    $tpath tag configure redcolor -foreground red
    $tpath tag configure bluecolor -foreground blue
    $tpath tag configure bluemargins2 -foreground blue -lmargin1 14m -lmargin2 6m -rmargin 10m
    
    if {[llength $allreportlist]>0} {
	
	set errorlist [list]
	set warninglist [list]
	# Add all the text with your format
	foreach items $allreportlist {
	    set keyword [lindex $items 0]
	    # WarnWinText "keyword:$keyword"
	    switch $keyword {
		"Error:" {
		    lappend errorlist [lrange $items 1 end]
		    
		}
		"Warning:" {
		    lappend warninglist [lrange $items 1 end]
		    
		}
	    }
	}
	# wa "errorlist:$errorlist"
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
		    $tpath insert end "${items}\n"
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
		$tpath insert end "${items}\n" 
		
	    }
	    $tpath insert end "\n***** End warning *****\n\n"  bigblue
	}
	
    } else {
	
	::KMValid::CreateReportWindowbClose $w
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
    
	foreach group $groups {
	    
	    set hasent [::KEGroups::getGroupGiDEntitiesNew $group "hasEntities"]
	    # msg "Hasent for $group $hasent"
	    if {$hasent eq "0"} {
		
		lappend noEntitiesGroups $group
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
    
    # Avisaremos de todos los errores al final
    set elementsGroups 0
    set groupwithoutEntities 0
    set groupwithoutEntitiesList [list]
    set conditionsGroups 0
    
    set hasMaterialId 0
    set noIdProps 0
    set noMaterialProps {}
    set nullProperty {}
    set nullPropertyDict [dict create]

    # StructuralAnalysis
    global KPriv
    
    # Set a local xml variable
    set xml $KPriv(xml)
    
    # Get the spatial dimension
    set ndime [::xmlutils::GetSpatialDimension]
    # wa "ndime:$ndime"
    
    # For elements

    # We look for the groups assigned to any element
    set xpath "$appPath/Container\[@id='Elements'\]/Container"
    set nodes [$xml selectNodes $xpath]
    
    foreach node $nodes {
	
	set idElem [$node getAttribute id ""]
	
	set groupsXPath "${xpath}\[@id='$idElem'\]/Container"
	set groups [::xmlutils::getXmlChildIds $xml "${groupsXPath}" ] 
	# msg "idElem$idElem   groups:$groups\n"
	foreach group $groups {
	    
	    set grNodeXPath "${groupsXPath}\[@id='$group'\]"
	    set active [::xmlutils::getAttribute $xml $grNodeXPath active]
	    # wa "grNodeXPath:$grNodeXPath"
	    if {$active} {
		        set hasEntities [::KEGroups::getGroupGiDEntitiesNew $group "hasEntities"]
		    # wa "hasEntities:$hasEntities Existe active en grNodeXPath $group"
		if {$hasEntities} {
		    #Si uno de los grupos tiene entidades ya no sacaremos ese error
		    
		    #Pero además tiene que estar activo para esa combinación de filtros (visible en el árbol)
		    #set isActiveGroup $groupsXPath...
		    set elementsGroups 1
		    set groupwithoutEntities 1
		    
		} else {
		    lappend groupwithoutEntitiesList [list $idElem $group]
		}
		
		# Now check the property assigned to this group
		set grNodeXPath "${groupsXPath}\[@id='$group'\]/Container"
		set grNodeConainers [$xml selectNodes "$grNodeXPath"]
		foreach nodeContainer $grNodeConainers {
		    
		    set nodeItems [$nodeContainer childNodes]
		    foreach item $nodeItems {
		    
		        set idItem [$item getAttribute id ""]
		        if { $idItem == "Property" } {
		            set propId [$item getAttribute dv ""]
		            
		            # Now look that property in 'Properties'
		            set xpath "$appPath/Container\[@id='Properties'\]/Container\[@id='$propId'\]"
		            # wa "propId:$propId xpath:$xpath"
		            set node [$xml selectNodes $xpath]
		            if { $node != "" } {
		                # Get the cross section property list
		                set PropertyList [::KMProps::GetCrossSectionPropertyList]
		                # wa "propId:$propId PropertyList:$PropertyList"
		                foreach nodeCont [$node childNodes] {
		                    set plist [list]
		                    foreach nodePropIt [$nodeCont childNodes] {
		                    
		                        # Here are the items of property
		                        set idItem [$nodePropIt getAttribute id ""]
		                        # wa "idItem:$idItem"
		                        # We found that the material is not null
		                        if { $idItem == "Material" } {
		                            if {  [$nodePropIt getAttribute dv ""] == "" } {
		                            
		                                #Si la propiedad no tiene material asignado, 
		                                #la añadimos a la lista de propiedades sin material
		                                if { !($propId in $noMaterialProps) } {
		                                    
		                                    lappend noMaterialProps $propId
		                                }
		                            }
		                        } elseif { $idItem == "ElemType" } {
		                            
		                            set ::KMProps::ElemTypeProperty [$nodePropIt getAttribute dv ""]
		                            
		                        } elseif { $idItem == "SectionType" } {
		                            
		                            set ::KMProps::SectionTypeProperty [$nodePropIt getAttribute dv ""]
		                            
		                        } elseif { $idItem in $PropertyList } {
		                            # Check only in the structural analysis application
		                            if {$application eq "StructuralAnalysis"} {
		                                # Get the element type
		                            set cxpath "$xpath/Container\[@id='MainProperties'\]/Item\[@id='ElemType'\]"
		                            # wa "cxpath:$cxpath"
		                            set cnode [$xml selectNodes $cxpath]
		                            # wa "cnode:$cnode"
		                                # Set the element type
		                            set ::KMProps::ElemTypeProperty [$cnode getAttribute dv ""]
		                            # wa "ElemTypeProperty:$::KMProps::ElemTypeProperty"
		                                
		                                # Set the section type
		                                set cxpath "$xpath/Container\[@id='MainProperties'\]/Item\[@id='SectionType'\]"
		                                # wa "cxpath:$cxpath"
		                                set cnode [$xml selectNodes $cxpath]
		                                # wa "cnode:$cnode"
		                                # Set the section type
		                                set ::KMProps::SectionTypeProperty [$cnode getAttribute dv ""]
		                                # wa "SectionTypeProperty:$::KMProps::SectionTypeProperty"
		                                
		                            # Check that this cross section property is necessary and this is not empty
		                            set ::KMProps::nDim $ndime
		                            set ShowProperty [::KMProps::ShowPropertyByElementType $idItem]
		                            set CurrentValue [$nodePropIt getAttribute dv ""]
		                            # wa "ShowProperty:$ShowProperty CurrentValue:$CurrentValue "
		                            if {(($ShowProperty) && ($CurrentValue == ""))} {
		                                # wa "inside propId:$propId"
		                                    # Check the section type
		                                    set ShowSectionProperty [::KMProps::ShowPropertyBySectionType $idItem]
		                                    # wa "ShowSectionProperty:$ShowSectionProperty idItem:$idItem"
		                                    if {$ShowSectionProperty} {
		                                # If the property needs tickness and do not have it, error
		                                if { !($propId in $nullProperty) } {
		                                    
		                                    lappend nullProperty $propId
		                                }
		                                # Update the property list
		                                lappend plist $idItem
		                                # Update the dictionary
		                                dict set nullPropertyDict $propId $plist
		                                # wa "nullProperty:$nullProperty"
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
    }
    
    if { ! $elementsGroups && !$groupwithoutEntities} {
	if {[llength $groupwithoutEntitiesList]} {
	    set bf ""
	    foreach cprop $groupwithoutEntitiesList {
	    lassign $cprop ElemId GroupId
	    append bf "Element:$ElemId -> Group:$GroupId "
	    }
	    lappend allreportlist "Error: There are no entities assigned to this element and group ($bf)."
	    lappend allreportlist "\n"
	}
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
    if { [llength $nullProperty] } {
	set bf ""
	foreach propid $nullProperty {
	    if {[dict exists $nullPropertyDict $propid]} {
		set plist [dict get $nullPropertyDict $propid]
		append bf "($propid -> [join $plist ","]) "
	    }
	}
	if {[string length $bf]} {
	    lappend propsMessages "Error: There are properties without a correct value:\n$bf"
	} else {
	    lappend propsMessages "Error: There are properties without a correct value:\n$nullProperty"
	}
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
	
	# Ahora ya no es un error, es un warning, por lo que se gestiona de forma distinta
	set allreportlist [::KMValid::checkGroups $xml $appPath "InitialConditions" $allreportlist]
    }

    #set xpath "$appPath/Container\[@id='Conditions'\]/Container"
    #if { [llength [::KMValid::groupsWithEntities $xml $xpath]] == 0 } {        
    #        lappend allreportlist "Error: There are no groups assigned to Conditions."
    #}
    
    #set line1 ""
    #set line2 ""
    #lappend allreportlist ${line1}
    #lappend allreportlist ${line2}
    # WarnWin "allreportlist:$allreportlist"
    return $allreportlist
}

#
# SEPARAR POR APLICACIONES (con color)
# Poner el path en los errores
#
proc ::KMValid::groupsWithEntities { xml xpath {returnNodes 0}} {
    
    set groupEntities {}
    
    #Buscamos los grupos asignados a algún element
    set nodes_inside_load [$xml selectNodes $xpath]
    set nodes [list]
    foreach node $nodes_inside_load {
	set idElem [$node getAttribute id ""]
	if { [$node getAttribute class]=="Groups" } {            
	    set groupsXPath "${xpath}\[@id='[list $idElem]'\]/Container"            
	    lappend nodes {*}[list $node $groupsXPath]
	} elseif { [$node getAttribute class]=="SameTemplateGroups" } {
	    foreach childnode [$node childNodes] {
		set childidElem [$childnode getAttribute id ""]
		set groupsXPath "${xpath}\[@id='[list $idElem]'\]/Container\[@id='[list $childidElem]'\]/Container"
		lappend nodes {*}[list $childnode $groupsXPath]
	    }            
	}
    }
	
    
    
    foreach {node groupsXPath} $nodes {
	
	set groups [$xml selectNodes $groupsXPath]
	#set groups [::xmlutils::getXmlChildIds $xml "${groupsXPath}" ]
	
	foreach nodeGroup $groups {
	    
	    set idGroup [$nodeGroup getAttribute id ""]
	    
	    #Se tiene que cumplir que el nodo esté activo (visible en árbol según los filtros seleccionados)
	    if {[$nodeGroup getAttribute active 0] != 0} {
		    set hasEntities [::KEGroups::getGroupGiDEntitiesNew $idGroup "hasEntities"]
		if {$hasEntities} {
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
	    if { [$node getAttribute dv ""]        == "" } {
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

		        lappend allreportlist "Error: $rootContainer-> $parentPID-> group '$groupId': Empty properties '$nullProps'"
		        
		        
		    }
		}

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

proc ::KMValid::SlipNoSlipList { {type "slip"} } {
    # De momento solo funciona con Slip y No Slip
	# Busca busca los grupos que pertenecen a $type
    global KPriv
    
    if { $type == "slip"} {
	set sliplist ""
	
	set root "/Kratos_Data/RootData\[@id='Fluid'\]/Container\[@id='Conditions'\]/Container\[@id='Is-Slip'\]"
	set slipnode [$KPriv(xml) selectNodes "$root"]
	set nodes [$slipnode childNodes]
	foreach node $nodes {
	    lappend sliplist [$node getAttribute id ""]
	}
	return $sliplist
    } else {
	
	set nosliplist ""
	set root "/Kratos_Data/RootData\[@id='Fluid'\]/Container\[@id='Conditions'\]/Container\[@id='No-Slip'\]"
	set noslipnode [$KPriv(xml) selectNodes "$root"]
	set nodes [$noslipnode childNodes]
	foreach node $nodes {
	    lappend nosliplist [$node getAttribute id ""]
	}
	
	return $nosliplist
    }
}


proc ::KMValid::isWall { nDim } {
    # De momento solo funciona con Slip y No Slip
	# Busca busca los grupos que pertenecen a $type
    global KPriv
    # msg "entro en wall"
    set wallList ""
    set root "/Kratos_Data/RootData\[@id='Fluid'\]/Container\[@id='Conditions'\]/Container\[@id='PFEMWall'\]"
    
    set Wallnode [$KPriv(xml) selectNodes "$root"]
    # msg "List [$Wallnode asXML]"
    set nodes [$Wallnode childNodes]
    foreach node $nodes {
	lappend wallList [$node getAttribute id ""]
    }
    return $wallList

}


proc ::KMValid::ValidateGroupRepeat { firstcomp secondcomp } {
	# De momento solo funciona con Slip y No Slip
	# Busca si hay grupos o nodos en comun entre firstcomp y secondcomp

    set msg [list]
    # Primero buscamos los grupos que pertenecen a Slip y a No Slip
    set listcomp1 [::KMValid::SlipNoSlipList "slip"]
    set listcomp2 [::KMValid::SlipNoSlipList "noslip"]
    # wa "listcomp1:$listcomp1 listcomp2:$listcomp2"

    # Check that the to list have some defined groups
    if {[llength $listcomp1] && [llength $listcomp2]} {
	set retval [::KUtils::TwoListRepeatedItems $listcomp1 $listcomp2 "groups"]
	if { $retval != ""} {
	    # Si algun grupo está en las 2 condiciones, para la ejecución
	    set txt [= "Warning: Repeated groups (%s) between %s and %s" $retval $firstcomp $secondcomp]
	    return "$txt" 
	} 
	
	# Ahora vamos a comporbar que no hayan nodos repetidos
	set msn [::KMValid::CheckRepeatedNodes $listcomp1 $listcomp2 $firstcomp $secondcomp]
	return $msn
    }
    return $msg
}


proc ::KMValid::CheckRepeatedNodes { grouplist1 grouplist2 Item1 Item2 } {
    # Comprueba que ningun nodo esté en $grouplist1 y en $grouplist2.
    # Need New GiD_Group Validate
    
    set l1nlist [list]
    set l2nlist [list]
    
    foreach cgroupid $grouplist1 {
	if {[GiD_EntitiesGroups get $cgroupid nodes -count]} {
	    lappend l1nlist {*}[GiD_EntitiesGroups get $cgroupid nodes]
	}
    }
    # wa "l1nlist:$l1nlist"
    
    # Check for repeated nodes only if the first list is not empty
    if {[llength $l1nlist]} { 
	foreach cgroupid $grouplist2 {
	    if {[GiD_EntitiesGroups get $cgroupid nodes -count]} {
		lappend l2nlist {*}[GiD_EntitiesGroups get $cgroupid nodes]
	    }
	}
	# wa "l2nlist:$l2nlist"
	
	# Check only if l2nlist is not empty
	if {[llength $l2nlist]} {
	    set rnodes [list]
	    set rnodes [::KMValid::FindRepeated $l1nlist $l2nlist]
	    if {[llength $rnodes]} {
		set txt [= "Warning: Repeated node (%s) between the %s and %s boundary condition. Check the group properties assigned to this boundary conditions" $rnodes $Item1 $Item2]
		return "${txt}"
	    } 
	}
    }

    return ""
}

proc ::KMValid::FindRepeated {nlist1 nlist2} {
    # Busca si hay nodos repetidos entre nlist1 y nlist2
    # Necesita eficiencia. En el peor caso puede recorrer cientos de miles de nodos
    set total [append nlist1 $nlist2]
    set total [lsort $total]
    set rec [lsort -unique $total]
    
    set a [llength $total]
    set b [llength $rec]
    if { $a != $b } {
	set i 1
	
	foreach item $total {
	    
	    if { $item == [lindex $total $i] } {
		return $item
		break
	    }
	    incr i 1
	}
    }
    return ""
}


proc ::KMValid::ValidateDragFiles { allreportlist } {
    # Validates if Drag Files field is empty
    
    global KPriv
    
    set root "/Kratos_Data/RootData\[@id='Fluid'\]/Container\[@id='Results'\]/Container\[@id='DragOptions'\]"
	set noded [$KPriv(xml) selectNodes "$root"]
    #msg [$noded asXML]
    set nodes [$noded childNodes]
    foreach node $nodes {
	
	set nitem [lindex [$node childNodes] 0]
	#msg [$nitem asXML]
	set item [$nitem childNodes]
	#msg [$item asXML]
	set dv [$item getAttribute dv ""]
	#msg [[lindex [$node childNodes] 0] getAttribute id ""]
	if {$dv eq ""} {
	    lappend allreportlist "Error: Empty Field: Drag Calculation Output File id."
	}
    }
    return $allreportlist
}

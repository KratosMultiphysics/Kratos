#####################################################################################
#
#  NAME: kmpropsgroups.tcl
#
#  PURPOSE: Manage the group options in the kratos main model window 
#
#  QUANTECH ATZ-DEVELOPMENT DEPARTMENT
#
#  AUTHORS : G. Socorro and L. Calvo
#
#  CREATED AT: 29/03/2012
#
#  HISTORY:
# 
#   0.9- 13/07/13- G. Socorro, add the proc GetAvailableGiDGroups, select all the group with normal or disabled state
#   0.8- 25/06/13- A. Melendo, new List and Draw procs
#   0.7- 18/06/13- G. Socorro, delete the use of the proc kipt::NewGiDGroups (delete the call to the compass groups => Cond_Groups)
#   0.6- 12/02/12- G. Socorro, modify the link to the GiD group window to use ::WinUtils::OpenGiDGroupTab
#   0.5- 10/10/12- J.  Garate, adaptation for New GiD Groups, Autogroup Frame Bug Corrected when assigning entities
#   0.4- 01/10/12- J.  Garate, deleting group function corrected
#   0.3- 26/04/12- G. Socorro, change GiD_Groups by Cond_Groups
#   0.2- 02/04/12- G. Socorro, correct a bug with the combobox path (update the autoNewGroup proc)
#   0.1- 29/03/12- G. Socorro, create a base source code from the kmprops.tcl script
#
######################################################################################
#                      Procedures that belong to this file
###############################################################################
#         Name                      |        Functionality
#------------------------------------------------------------------------------
# 1.            | 


proc ::KMProps::GetAvailableGiDGroups {} {
    set allgrouplist [list] 
    foreach groupid [GiD_Groups list] {
	if {[GiD_Groups get state $groupid] != "hidden"} {
	    lappend allgrouplist $groupid
	}
    }
    return $allgrouplist
}

proc ::KMProps::changeGroups { entityList f {fullname ""} } {
    variable selGroup

    # Get the GiD group list
    set allgrouplist [::KMProps::GetAvailableGiDGroups] 
    if {[llength $allgrouplist]} {
	$f configure -values $allgrouplist
    }
    
    if {$selGroup !=""} {
	if { !($selGroup in $allgrouplist) } {
	    WarnWin [= "The new group '%s' has not any usefull entity assigned." $selGroup]
	    set selGroup ""
	}
    }
 }
 
 proc ::KMProps::cmbGetGidGroups {cbpath fullname} {
   if { [winfo exists $cbpath] } {
	
	# Get the GiD group list
	set allgrouplist [::KMProps::GetAvailableGiDGroups] 
	# wa "$allgrouplist"
	if {[llength $allgrouplist] } {
	    $cbpath configure -values $allgrouplist
	    if {[$cbpath get] == ""} {
		$cbpath set [lindex $allgrouplist 0]
	    }
	} else {
	    $cbpath configure -values {}
	}
    }
 }

proc ::KMProps::cmbChangeCheckGroups { f } {
    variable selGroup
    
    if { [winfo exists $f.cGroups] } {
	
	# Get the GiD group list
	set allgrouplist [::KMProps::GetAvailableGiDGroups] 
	
	if {[llength $allgrouplist] } {
	    
	    $f.cGroups configure -values $allgrouplist
	    if {$selGroup ni $allgrouplist} {
		set selGroup [lindex $allgrouplist 0]
	    }
	} else {
	    $f.cGroups configure -values {}
	    set selGroup ""
	}
    }
}


proc ::KMProps::setNoActiveGroups { node } {
    
    if {[$node getAttribute class ""] == "Groups" } {
	
	foreach nod [$node childNodes] {
	    if {[$nod getAttribute class ""] == "Group" } {
		$nod setAttribute active 0
	    }
	}
    }
}

proc ::KMProps::getGroups { entityList {fullname ""}} {
    
    global KPriv
    
    # Switch state
    if { $entityList == "" } {
	set PState [GiD_Info Project ViewMode]
	if {($PState == "GEOMETRYUSE")} {
	    set entityList {point line surface volume}
	} else {
	    set entityList {nodes elements}
	}
    }
    
    set grupos {}
	set grw "GiD_Groups"
    foreach groupId [$grw list] {
	foreach entity $entityList {
	    if { [::KEGroups::getGroupGiDEntities $groupId $entity "hasEntities"] } {
	    if { !( $groupId in $grupos) } {
		lappend grupos $groupId
	    }
	    }
	}
    }
    
    
    if {$fullname != ""} {
	
	#Eliminamos de la lista los grupos ya asignados a esta propiedad
	set assignedGroups [::xmlutils::setXmlContainerIds $fullname]
	foreach g $assignedGroups {
	    if { $g != $::KMProps::selGroup } {
	    set grupos [::KEGroups::listReplace $grupos $g]
	    }
	}
    }
    return $grupos
}

proc ::KMProps::autoNewGroup { group_name fpath } {
    variable selGroup
    global KPriv
    #WarnWinText "group_name: $group_name       fpath: $fpath"
    set GroupId [::KEGroups::GetAutomaticGroupName "$group_name"]
    
    # Assign the selected group identifier
    set selGroup $GroupId 
    
    # Create the new group
    GiD_Groups create $GroupId
    
    # Selection the entities to be assigned
    ::KEGroups::SelectionAssign $::KMProps::selectedEntity $GroupId $fpath
    
    # Ponemos el foco en la ventana de propiedades
    #focus $winpath
    GidUtils::UpdateWindow GROUPS
}

proc ::KMProps::acceptGroups { T idTemplate fullname item listT entityList fGroups} {
    variable selGroup
    variable NbPropsPath
    global KPriv

    set grupo $selGroup
    
    if { $grupo == "" } {
	WarnWin [= "You have to choose one group\n (you can create a new one pushing the button on the right)"]
    } else {
	
	#Primero comprobamos q el grupo aun no exista
	#set id [::xmlutils::setXml "$fullname//c.[list $grupo]" id]
	set path "$fullname//c.[list $grupo]"
	set xpath "[::xmlutils::setXPath $path]"
	set check [$KPriv(xml) selectNodes $xpath]        
	
	if { [llength $check]!=0 } {	    
	    WarnWin [= "This group is already assigned to this condition."]
	} else {
	    
	    #Validamos que haya alguna propiedad seleccionada
	    if {[info exists ::KMProps::cmbProperty]} {
	    if {$::KMProps::cmbProperty == "" } {
		WarnWin [= "You must define a Property before!"]
		return ""
	    }
	    }
	    
	    #Comprobamos que el grupo no sea un AutoGroup sin entidades
	    ::KMProps::changeGroups $entityList $fGroups
	    if { $::KMProps::selGroup == "" } {
		return ""
	    }
	    
	    set template [::KMProps::copyTemplate ${idTemplate} $fullname "$grupo" "Group"]
	    #Validamos que haya algún combo de "properties" en el template
	    if {[string match "*GCV=\"Properties*" $template]} {
	    
	    #Miramos si hay alguna propiedad dada de alta (si la hay,obligatoriamente estará seleccionada)
	    set props [::xmlutils::setXmlContainerIds "[::KMProps::getApplication $fullname]//c.Properties"]
	    if { [llength $props] < 1 } {
		WarnWin [= "You must define a Property before."]
		#::KMProps::deleteProps $T $itemSel $newPropertyName
		return  ""
	    }
	    }
	    
	    
	    #
	    # Ahora debemos actualizar todos los valores en el xml
	    #
	    set fBottom ${NbPropsPath}.fBottom
	    
	    if {[llength $listT] >= 1 } {
	    
	    #Lista de listas con formato {idContainer idItem1 idItem2...}
	    foreach listContainer $listT {
		
		#Si tiene como mínimo el container y un item entramos
		if {[llength $listContainer] >= 2} {
		    
		    set idContainer [lindex $listContainer 0]
		    
		    #set id [::KMProps::getPropTemplate $idTemplate id $idContainer]
		    
		    #Recorremos los items
		    for {set i 1} { $i < [llength $listContainer] } {incr i} {
		        
		        set id [lindex $listContainer $i]
		        
		        #Cuando tengamos nodos ocultos o incompatibles la variable no existirá
		        if {[info exists ::KMProps::cmb$id]} {
		            
		            set fullNombre "$fullname//c.[list $grupo]//c.[list $idContainer]//i.[list $id]"
		            
		            set f "${fBottom}.nb.f${idContainer}.fc.c.fsc.cmb${id}"
		            
		            if { [winfo exists $f] } {
		                
		                if { [$f cget -state] == "readonly" } {
		                    
		                    set value [::xmlutils::getComboDv $f $fullNombre]
		                } else {
		                    
		                    set value [set ::KMProps::cmb$id]                                                                          
		                }
		                
		                if {$id == "Vx" || $id == "Vy" || $id == "Vz"} {
		                    
		                    set activeId "A[string range $id 1 1]"
		                    set fullActive "$fullname//c.[list $grupo]//c.Activation//i.[list $activeId]"
		                    set active [::xmlutils::setXml $fullActive dv "read"]
		                    
		                    if { $active == 0 } {
		                        ::xmlutils::setXml $fullNombre state "write" "disabled"
		                    } else {
		                        ::xmlutils::setXml $fullNombre state "write" "normal"
		                    }
		                }
		                #Comprobamos si el combo tiene una función asignada
		                set function [::xmlutils::setXml $fullNombre function]
		                if { $function != "" && [$f cget -state] == "disabled"} {
		                    
		                    ::xmlutils::setXml $fullNombre function "write" 1
		                    ::xmlutils::setXml $fullNombre state "write" "disabled"
		                } 
		                
		                #Guarda el nuevo valor en el xml
		                ::xmlutils::setXml $fullNombre dv "write" $value
		            }
		        }
		    }
		}
	    }
	    }

	    #Destruimos el frame inferior
	    ::KMProps::DestroyBottomFrame
	    
	    ::KMProps::RefreshTree $T
	    
	    $T selection add $item
	    $T item expand $item
		# Open the new group tab
		::WinUtils::OpenGiDGroupTab
	}
    }
}

#
# Borrar el grupo de la condición preguntando si lo queremos borrar completamente
#
proc ::KMProps::deleteGroupCondition { T item } {
    set fullname [DecodeName [$T item tag names $item]]
    set GroupId [$T item text $item 0]
    
    set aviso "Removing group $GroupId. Please, choose the properly option:\n\n\n Yes: Desassign this group (recomended)\n\n No: Complete group removing (with all its descendants).\n\n Cancel: Keep the group assigned."
    set answer [::WinUtils::confirmBox "." "$aviso" "yesnocancel"]
    if { $answer == "yes" } {
	
	#Desasigna de la gemoetría el item seleccionado
	#::KEGroups::UnAssignCondition $GroupId
	
	#Elimina el grupo del xml
	::xmlutils::unsetXml [DecodeName [$T item tag names $item]]
	
	#Elimina el grupo del árbol
	::KMProps::deleteItem $T $item
	
    } elseif {$answer == "no" } {
    
	::KEGroups::BorraGrupo [lindex [split [DecodeName [$T item tag names $item]] "."] end]
	::xmlutils::unsetXml [DecodeName [$T item tag names $item]]
	::KMProps::deleteItem $T $item
	
	
    } else {
    }
    
}

#
# Dibujado el grupo de la condición
#
proc ::KMProps::drawGroupCondition { T item } {
    #msg "$T $item"
    set fullname [DecodeName [$T item tag names $item]]
    set GroupId [$T item text $item 0]
    
    #GiD_Groups draw $GroupId
    #GiD_Redraw
    
    set parent [winfo parent [winfo parent [winfo parent $T]]]    
    #FinishButton [winfo parent [winfo parent $parent]] $parent.caption [_ "Press 'Finish' to end selection"] {GiD_Groups end_draw; GiD_Redraw} disableall [GiD_Set SmallWinSelecting]
    GiD_Process 'Groups Draw {*}$GroupId escape
    
    FinishButton $parent $parent.caption [_ "Press 'Finish' to end selection"] "" disableall [GiD_Set SmallWinSelecting]
    
}
#
# Dibujado los grupos de una condición
#
proc ::KMProps::drawGroupsCondition { T item } {
    set ListGroupId ""
    foreach childitem [$T item children $item] {
      set fullname [DecodeName [$T item tag names $childitem]]
      set class [::xmlutils::setXml $fullname class]
      if {$class == "Group" } {
	lappend ListGroupId [$T item text $childitem 0]
      }
    }
    if { [llength $ListGroupId] > 0 } {
      set parent [winfo parent [winfo parent [winfo parent $T]]]
      GiD_Process 'Groups Draw {*}$ListGroupId escape
    
      FinishButton $parent $parent.caption [_ "Press 'Finish' to end selection"] "" disableall [GiD_Set SmallWinSelecting]
    }
}

proc ::KMProps::auxiliarfunction_ReturnEntitiesInsideGroups { ListGroupId } {

  set ReturnList ""

  set types [list [_ "Points"] [_ "Lines"] [_ "Surfaces"] [_ "Volumes"] [_ "Nodes"] [_ "Elements"] [_ "Faces"]]
  
  set totallistentities [list {*}[GiD_EntitiesGroups get [lindex $ListGroupId 0] all_geometry] {*}[GiD_EntitiesGroups get [lindex $ListGroupId 0] all_mesh]]
  #convert faces format
  if { [llength [lindex $totallistentities end]] > 0} {
    set txt ""
    foreach element [lindex [lindex $totallistentities end] 0] face [lindex [lindex $totallistentities end] 1] {
      #convert to real format to sort
      lappend txt $element.[list $face]
    }
    set totallistentities [lreplace $totallistentities end end $txt]
  }
  
  foreach GroupId [lrange $ListGroupId 1 end] {   
    set listentities [list {*}[GiD_EntitiesGroups get [lindex $GroupId 0] all_geometry] {*}[GiD_EntitiesGroups get [lindex $GroupId 0] all_mesh]]
    if { [llength [lindex $listentities end]] > 0} {
      set txt ""
      foreach element [lindex [lindex $listentities end] 0] face [lindex [lindex $listentities end] 1] {
	#convert to real format to sort
	lappend txt $element.[list $face]
      }
      set listentities [lreplace $listentities end end $txt]
    }
    set count 0
    foreach entititype $listentities {    
      set totallistentities [lreplace $totallistentities $count $count [ lsort -increasing -unique -real [list {*}[lindex $totallistentities $count] {*}$entititype ]]]
      incr count
    }
  }
  set count 0
  foreach type $types {
    if {[llength [lindex $totallistentities $count]]>0} {
      #convert from real format to face format
      append ReturnList "  " $type ": " [regsub -all {\.} [lindex $totallistentities $count] :] "\n"
    }
    incr count
  }
  return $ReturnList
 
}


proc ::KMProps::listGroupCondition { T item } {
    #msg "$T $item"
    #set fullname [DecodeName [$T item tag names $item]]
    set GroupId [$T item text $item 0]
    
    set listentities [::KMProps::auxiliarfunction_ReturnEntitiesInsideGroups $GroupId]
    
    # wa "${GroupId}: "
    # wa $listentities
}


proc ::KMProps::listGroupsCondition { T item } {

    set ListGroupId ""
    
    foreach childitem [$T item children $item] {
      set fullname [DecodeName [$T item tag names $childitem]]
      set class [::xmlutils::setXml $fullname class]
	    if {$class == "Group" } {
	lappend ListGroupId [$T item text $childitem 0]
	    }
    }
    if { [llength $ListGroupId] > 0 } {
      set listentities [::KMProps::auxiliarfunction_ReturnEntitiesInsideGroups $ListGroupId]
      set fullname [DecodeName [$T item tag names $item]]
      set name [::xmlutils::setXml $fullname pid]
    
      WarnWinText "$name [_ "has been applied to"]:"
      WarnWinText "  [_ "Groups"]: $ListGroupId"
      WarnWinText $listentities
    }
}
proc ::KMProps::auxiliarfunction_GetItemTextFromTreeItem { T item} {
  set text ""
  if { [$T item id $item] != "" } {  
      set fullname [DecodeName [$T item tag names $item]]
      set dv [::xmlutils::setXml $fullname dv]
      set pid [::xmlutils::setXml $fullname pid]
      if { $dv != "" } {        
	append text $pid ": " $dv
      }    
  }
  return $text
}
proc ::KMProps::auxiliarfunction_GetContainerTextFromTreeItem { T item} {
  set text ""
  if { [$T item id $item] != "" } {
    set fullname [DecodeName [$T item tag names $item]]
    set pid [::xmlutils::setXml $fullname pid]   
    if { $pid=="" } {
      set pid [_ "Root"]  
    }
    append text $pid    
  }
  return $text
}
proc ::KMProps::auxiliarfunction_WritteChildren { T item level} {
  incr level
  set numspaces "  "
  set separatoritem ", "
  
  set fullname [DecodeName [$T item tag names $item]]
  
  if { $fullname!="" && [::xmlutils::getXmlNodeName $fullname] == "Item" } {
    set information [::KMProps::auxiliarfunction_GetItemTextFromTreeItem $T $item]
    return $information
  } else {
    #like Container
    set tab ""
    for { set i 1 } { $i < $level } { incr i } {        
      append tab $numspaces
    }
    set text "$tab"    
    set importantinformation 0
    append text [::KMProps::auxiliarfunction_GetContainerTextFromTreeItem $T $item] ":"
    set listchildrens [$T item children $item]
    foreach children $listchildrens { 
      set information [::KMProps::auxiliarfunction_WritteChildren $T $children $level]
      if { $information != "" } {        
	set fullnamechildren [DecodeName [$T item tag names $children]]
      
	  if { [::xmlutils::getXmlNodeName $fullnamechildren] == "Item"} {
	    if { $importantinformation==0 } {
	      #tab de level + 1
	      append text "\n$tab" 
	      append text $numspaces
	    } else {
	      append text $separatoritem
	    }
	  } else {
	    append text "\n"           
	  }
		  
	set importantinformation 1
	append text $information
      }
    }
    if {$importantinformation==0} {
      set text ""
    }
  }
  return $text
}

proc ::KMProps::listSubtree { T item } {
  set listId ""

  set listsubitems [::KMProps::auxiliarfunction_WritteChildren $T $item 0]
	
  # wa WarnWinText $listsubitems

}


#
# Busca en el arbol de propiedades, todos los grupos que se llamen $GroupId
#
proc ::KMProps::findGroups { GroupId } {
    
    global KPriv
    
	set ElemList { }
    set bool 0
    #Primero hay que asegurarse de que exista el árbol
	
	#Si está el frameBottom activo nos lo cargamos
	::KMProps::DestroyBottomFrame
	
	#Si no está abierta la ventana, deberemos recorrer el xml a mano
	set nodes [$KPriv(xml) getElementsByTagName "Container"]
	
	foreach container $nodes {
	    
		#if {[$container getAttribute idTemplate ""] == ""} {continue}
		foreach node [$container childNodes] {
		        catch {
		                #Si encuentra el grupo intenta borrar el nodo y su descendencia
		                if { [$node getAttribute id ""] == "$GroupId" } {
		                        if { [$node getAttribute class ""] == "Group" } {
		                                lappend ElemList $node
		                                set bool 1
		                        }
		                }
		        }
		}
	}
	
	if { $bool == 0 } {
		set ElemList 0
    }
    return $ElemList
}


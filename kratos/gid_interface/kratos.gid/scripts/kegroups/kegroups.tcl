##############################################################################################
#
#	NAME: kegroups.tcl
#
#	PURPOSE: Utilities procedures to work with the Kratos entities groups editor
#
#	QUANTECH ATZ-DEVELOPMENT DEPARTMENT
#
#	AUTHOR : G. Socorro
#
#	CREATED AT: 01/11/09
#
#	LAST MODIFICATION : add the procedure ::KEGroups::RenameGroupIdGiDCond to rename the group identifier in the GiD condition database
#
#	VERSION : 0.3
#
#	HISTORY:
#
#        0.3- 13/05/10-G. Socorro, add the procedure ::KEGroups::RenameGroupIdGiDCond to rename the group identifier in the GiD condition database
#	 0.2- 19/03/10-LCA, Reparar acciones en el arbol (delete masivo, rename, ...) y aumentar
#												 el número de niveles a 5
#	 0.1- 01/11/09-G. Socorro, create a base source code from the GiD layer.tcl script
#	 
#
##############################################################################################

package require treectrl
package require tooltip
package provide KEGroups 1.0 
package require snit
package require tdom
package require xmlstruct 1.0

# Create a snit type to store the entity group properties
#snit::type EGroup {}

# Create a base namespace KEGroups
namespace eval ::KEGroups:: {
 
	# Path of the base window 
	variable WinPath ".gid.kegroups"
	variable WinLayout 
	variable SystemHighlight
	variable SystemHighlightText
	variable TreePath   ;# The tree path
}

proc ::KEGroups::Init {} {
	
	variable WinLayout;		   variable SystemHighlight
	variable SystemHighlightText; 
	variable Props
	
	
	# Get default colors
	set w [listbox .listbox]
	set SystemHighlight [$w cget -selectbackground]
	set SystemHighlightText [$w cget -selectforeground]
	destroy $w
	if { $::tcl_platform(platform) == "unix" } {
	 # I hate that gray selection color
	 set SystemHighlight #316ac5
	 set SystemHighlightText White
	}

	set WinLayout "OUTSIDE"
}

proc ::KEGroups::InitBaseWindow {{what "OUTSIDE"} } {
	
	
	global KPriv
	variable WinLayout
	
	#msg [$KPriv(xml) asXML]
	
	# Init KEGroups namaspace variables
	::KEGroups::Init
	
	#Si el problemtype ya estaba guardado, guardamos el xml y lo volvemos a cargar
	#set dirGid [GiD_Info problemtypepath]
	#msg "dirGid:$dirGid"
	#::kfiles::SaveSPD "$KPriv(dir)"
	
	set w "$::KEGroups::WinPath"

	if {$what == "OUTSIDE"} {
		
		# Open the window outside
		::KEGroups::OpenWindowOutside $w
	
		if {![winfo exists $w]} return 
	
		::KEGroups::CreateTreeAndToolbar $w
		
		::KEGroups::FillTree
		
		## For lower buttons
		#set tf [ttk::frame $w.buttons] 
		#grid $tf -sticky ews
	#
		#grid anchor $tf center
		#
		#grid [ttk::button $tf.bClose -text [= "Close"] -command [list destroy $w]]  -sticky ew -padx 5 -pady 6
		
		# Binding
		bind $w <Alt-c> "destroy $w"
		bind $w <Escape> "destroy $w"
		
	}
}

proc ::KEGroups::CreateTreeAndToolbar { w } {
	
	variable TreePath   


	 # Create the treectrl properties
	set mdf [ttk::frame $w.middle]
	set T [::KEGroups::CreateTreeProperties $w]
	# Set the tree path to a namespace variable
	set TreePath $T
	
	grid $w.middle -sticky wes
  
	# Create the toolbar 
	# Frame
	set tbf [ttk::frame $w.tbar -borderwidth 0]

	ttk::button $tbf.newgroupid -image [::WinUtils::GetImage new_tree.gif]  -command [list ::KEGroups::CreateNewGroupId $T] -style Toolbutton
	tooltip::tooltip $tbf.newgroupid [= "Add a new group identifier"]
	
	ttk::button $tbf.deleteGroupsId -image [::WinUtils::GetImage delete_tree.gif]  -command [list ::KEGroups::DeleteGroupsId $T] -style Toolbutton
	tooltip::tooltip $tbf.deleteGroupsId [= "Delete the selected group"]
	  
	ttk::menubutton $tbf.assign -image [::WinUtils::GetImage assign.gif] -menu $tbf.assign.m 
	
	menu $tbf.assign.m -postcommand "::KEGroups::rebuildMenu $tbf.assign.m"
	#$tbf.assign.m add command -command [::KEGroups::assignEntities "$tbf.assign.m"]
	
	ttk::button $tbf.unassign -image [::WinUtils::GetImage unassign.gif] -command [list ::KEGroups::unAssignEntities ] -style Toolbutton
	tooltip::tooltip $tbf.unassign [= "Unassing Entities from Groups"]
	
	ttk::button $tbf.listEntities -image [::WinUtils::GetImage list_entities.gif] -command [list ::LEntities::InitBaseWindow] -style Toolbutton
	tooltip::tooltip $tbf.listEntities [= "List Entities from Groups"]
	
	
	# Grid for toolbar 
	grid $w.tbar -row 2 -column 0 -sticky wes
	grid $tbf.newgroupid -sticky we -row 0 -column 0
	grid $tbf.deleteGroupsId -sticky we -row 0 -column 1
	grid $tbf.assign -sticky we -row 0 -column 2
	grid $tbf.unassign -sticky we -row 0 -column 3
	grid $tbf.listEntities -sticky w -row 0 -column 4
	
	grid [ttk::button $tbf.bClose -text [= "Close"] -command [list destroy $w]]  -row 1 -column 3 -columnspan 3 -sticky e -padx 5 -pady 3
	grid $tbf -sticky ews
	grid anchor $tbf.bClose center
	
	focus $T
}

proc ::KEGroups::rebuildMenu { w } {
	
	$w delete 0 10
	$w add command -command [::KEGroups::assignEntities "$w"]
	return ""
	
	destroy $w
	menu $w -postcommand "::KEGroups::rebuildMenu $w"
	$w add command -command [::KEGroups::assignEntities $w]
	#"::KEGroups::rebuildMenu $whatuse $w $command" 
	return ""
	set whatuse [GiD_Info Project ViewMode]
	
	if { $lastUse == $whatuse } {
		#Asignar entidades normalmente
		::KEGroups::GroupsSelectionAssign $command
	} else {
		WarnWin [= "The GiD ViewMode has changed to $whatuse.\nPush the menu button again."]
		#No dejamos asignar una entidad incorrecta, reconstuimos el menu
		destroy $w
		#menu $w -postcommand [::KEGroups::assignEntities "$w"]
		#$w add command -command [::KEGroups::assignEntities "$w"]
		#$w configure -state active
		#${::KEGroups::WinPath}.tbar.assign configure invoke
	}
	
}

#   
###################################################################################################
#--------------------------------------------------------------------------------------------------
# ASIGNACIÓN DE GRUPOS
#--------------------------------------------------------------------------------------------------
###################################################################################################
#
proc ::KEGroups::assignEntities { w } {
	 
	#Create the bitmaps menu
	set whatuse [GiD_Info Project ViewMode]
	
	set geomlist [GetGeometryImageFiles]
	set meshlist {node.gif element.gif}
	set postlist {node.gif element.gif resultonview.gif}
	set graflist {node.gif resultonview.gif}
	
	#$w delete 0 end
	#msg "-->What:$whatuse"
	switch $whatuse {
	GEOMETRYUSE {
		foreach i $geomlist {
		set command [file rootname $i]
		::KEGroups::AddImageCommandToMenu $w [GetImage $i] "::KEGroups::GroupsSelectionAssign $command"
		
		# $w add command -image [GetImage $i] -command "$basecommand" "$command" -hidemargin 1
		}
	}
	MESHUSE {
		foreach i $meshlist {
		set command [file rootname $i]
		::KEGroups::AddImageCommandToMenu $w [GetImage $i] "::KEGroups::GroupsSelectionAssign $command"
		}
	}
	POSTUSE {
		if { $args == "list"} {
		set postlist $meshlist
		}
		foreach i $postlist {
		set command [file rootname $i]
		
		#::KEGroups::AddImageCommandToMenu $w [GetImage $i] "$basecommand" "$command"
		# $w add command -image [GetImage $i] -command "$basecommand $command" -hidemargin 1
		}
	}
	GRAFUSE {
		foreach i $graflist {
		set command [file rootname $i]
		#::KEGroups::AddImageCommandToMenu $w [GetImage $i] "$basecommand" "$command"
		# $w add command -image [GetImage $i] -command "$basecommand $command" -hidemargin 1
		}
	}
	}
}

proc ::KEGroups::AddImageCommandToMenu { menu img cmd { st normal}} {
	
	if { $::tcl_platform(os) != "Darwin" } {
	if { $cmd != ""} {
		$menu add command -image $img -hidemargin 1 -state $st -command "$cmd"
	} else {
		$menu add command -image $img -hidemargin 1 -state $st
	}
	} else {
	if { $cmd != ""} {
		$menu add command -label [ file root [ file tail [ $img cget -file]]] -hidemargin 1 -state $st -command "$cmd"
	} else {
		$menu add command -label [ file root [ file tail [ $img cget -file]]] -hidemargin 1 -state $st
	}
	}
	return cm
}

proc ::KEGroups::GroupsSelectionAssign {entity} {
	
	variable TreePath
	
	set item [$TreePath selection get 0]
	
	if {$item != ""} {
		
		set GroupId [$TreePath item text $item 0]
			
		::KEGroups::SelectionAssign $entity $GroupId $::KEGroups::WinPath
		
	} else {
		
		WarnWin [= "Any group selected"]
	}
}

proc ::KEGroups::SelectionAssign { entity GroupId WinPath } {
	
	# Get the condition identifier
	#set condname [::groupProp::EditSelectionGetConditionId $groupid "Assign"]
	
	if { $entity == "element" } {
		WarnWin [= "Element assignation unavailable"]
		return ""
	}
	
	if { $entity == "node" } {
		set condname "point_groups"
	} else {
		set condname "${entity}_groups"
	}
	
	set OldSmallWinSelecting [GiD_Info variable SmallWinSelecting]
	
	#set OldSmallWinSelecting $entity
	if {$OldSmallWinSelecting == 0 } {
		set SmallWinSelecting 1
		GiD_Set SmallWinSelecting $SmallWinSelecting
	} else {
		set SmallWinSelecting 1
	}
	
	# ::GidUtils::DisableWarnLine
	FinishButton $WinPath "" [= "Press 'Finish' to stop the entities selection"] "::GidUtils::EnableWarnLine" disableall $SmallWinSelecting
	# Try to assign conditions
	
	GiD_Process MEscape
	GiD_Process Data Conditions AssignCond $condname NoRepeatField groupid
	GiD_Process change $GroupId
	#msg "GroupId:  $GroupId   .  condname:$condname \nentity:$entity WinPath:$WinPath"
		
	# ::GidUtils::EnableWarnLine

}

proc ::KEGroups::unAssignEntities { } {
	
	variable TreePath
	set item [$TreePath selection get]
	set grupos ""
	
	#Tiene que haber algún elemento seleccionado
	if {[llength $item] > 0 } {
		
		for {set i 0} { $i < [llength $item] } {incr i} {
				
			set GroupId [$TreePath  item text [lindex $item $i] 0]
				::KEGroups::UnAssignCondition $GroupId
				set grupos "$grupos '$GroupId'"
		}
		WarnWin [= "Succesfull unassigned!\nGroups: %s" $grupos]
	}
}

proc ::KEGroups::UnAssignCondition {groupId} {

	foreach entity [list point line surface volume] {
		set entityList [::KEGroups::getGroupGiDEntities $groupId $entity]
		if { [llength $entityList] } {
			eval [list GiD_Process MEscape Data Conditions AssignCond \
				  ${entity}_groups UnAssign Field groupid $groupId] $entityList \
			escape
		}
	}
	GiD_Process MEscape
}

proc ::KEGroups::FreezeLayers {frozen_layers} {
	foreach layer $frozen_layers {
	GiD_Process Layers Freeze $layer escape
	}
}

proc ::KEGroups::UnFreezeLayers {} {
	
	set frozen_layers ""
	foreach layer [GiD_Info Layers] {
	set frozen [lindex [GiD_Info Layers $layer] 0 1]
	if { $frozen } {
		GiD_Process Layers UnFreeze $layer escape
		lappend frozen_layers $layer
	}
	}
	return $frozen_layers
}


#
# Devuelve los id's de entidades que estén a asignadas a alguno de los grupos en "groupsIds"
##OLD CODE
proc ::KEGroups::getAssignedGiDEntities { groupIds } { 
	
	#OLD CODE
	::GidUtils::DisableGraphics
	set frozen_layers [::KEGroups::UnFreezeLayers]
	set infoProj [GiD_Info Project ViewMode]
	
	set entities {}
	foreach geom_mesh [list geometry mesh] {
	if {($infoProj == "GEOMETRYUSE") && ($geom_mesh == "mesh")} {
		GiD_Process MEscape Meshing MeshView
	} elseif {($infoProj == "MESHUSE") && ($geom_mesh == "geometry")} {
		GiD_Process MEscape Geometry
	}
	
	set points {}
	set lines {}
	set surfaces {}
	set volumes {}
	foreach i [list point line surface volume] {
		
		foreach j [GiD_Info conditions ${i}_groups $geom_mesh] {
			foreach "- num - group" $j break 
		
				if { $group in $groupIds } {
					
					switch $i {
						point {
							lappend points $num
						} 
						line {
							lappend lines $num
						}
						surface {
							lappend surfaces $num
						}
						volume {
							lappend volumes $num
						}
					} 
				}
		
		}
	}
	lappend entities [list $points $lines $surfaces $volumes]
	
	if {($infoProj == "GEOMETRYUSE") && ($geom_mesh == "mesh")} {
		GiD_Process MEscape Geometry
	} elseif { ($infoProj == "MESHUSE") && ($geom_mesh == "geometry") } {
		GiD_Process MEscape Meshing MeshView				
	}
	}
	   
	::KEGroups::FreezeLayers $frozen_layers
	::GidUtils::EnableGraphics
	
	return $entities
}

proc ::KEGroups::getGroupGiDEntities {groupId {givenEntity "point"} {action ""}} {
	
    # Update the link to the GiD condition properties for this group identifier

    # Disable graphics
    ::GidUtils::DisableGraphics
    # Store the freeze layers
    set flayerslist [::KEGroups::UnFreezeLayers]
    # Get the project view mode
    set PState [GiD_Info Project ViewMode]
       
	
	# Switch state
	if {($PState == "GEOMETRYUSE")} {
		
		set gmid "geometry"
		
	} elseif {($PState == "MESHUSE")} {
		
		set gmid "mesh"
	}
	
	# For each GiD group entities
	set EntityIdList {}
	
	
	if { $givenEntity == "ALL" } {
		set entities [list point line surface volume]
		
    } elseif { $givenEntity == "elements" } {
    	
    	set entities [list line surface volume]
    } elseif { $givenEntity == "nodes" } {
    	
    	set entities "point"
    } else {
    	
    	set entities $givenEntity
    }
    
    foreach entity $entities {
	    foreach CondProp [GiD_Info conditions ${entity}_groups $gmid] {
		lassign $CondProp - CEId - CGroupId
			# Update the entities list
			if {$CGroupId == $groupId } {
				
				#Si únicamente nos interesaba saber si el grupo tenía entidades acabamos aquí
				if {$action == "hasEntities" } {
					# Restore the layer state
				    ::KEGroups::FreezeLayers $flayerslist
				    # Enable graphics
				    ::GidUtils::EnableGraphics 
					return 1
				}
					
		    	lappend EntityIdList $CEId
		    }
	    }
    }
    #WarnWinText "EntityList:$EntityList"
	
    # Restore the layer state
    ::KEGroups::FreezeLayers $flayerslist
    # Enable graphics
    ::GidUtils::EnableGraphics 
    
    if {$action == "hasEntities" } {
    	return 0
	}
	
    return $EntityIdList
}


#
###################################################################################################
#--------------------------------------------------------------------------------------------------
# Funciones para manejar el XML
#--------------------------------------------------------------------------------------------------
###################################################################################################
#
# Prepara la query para utilizar las funciones de domNOde
proc ::KEGroups::setGroupsXPath { path } {
	
	set splitted [::KEGroups::split2 $path //]
	
	if { [llength $splitted] >= 1 } {
		
		set xpath "/Kratos_Data/Groups/Group\[@id='[lindex $splitted 0]'\]"
	} 
	if { [llength $splitted] >= 2 } {
		
		set xpath "$xpath/Group\[@id='[lindex $splitted 1]'\]"
	
	} 
	if { [llength $splitted] >= 3 } {
		
		set xpath "$xpath/Group\[@id='[lindex $splitted 2]'\]"
	}
	if { [llength $splitted] >= 4 } {
		
		set xpath "$xpath/Group\[@id='[lindex $splitted 3]'\]"
	}
	if { [llength $splitted] >= 5 } {
		
		set xpath "$xpath/Group\[@id='[lindex $splitted 4]'\]"
	}

	return $xpath
}
#
# Extraer los grupos del xml
#
proc ::KEGroups::getGroupsType { } {
	
	global KPriv
	
	set xpath ""
	
	set grupos [$KPriv(xml) set "/Kratos_Data/GroupTypes/Item/@pid" ]
	#msg "Grupos:$grupos"
	
	return $grupos
	#return {"BC" "BC2" "load" "Default"}
}

#
# Editar o extraer propiedades del xml en memoria
#
proc ::KEGroups::setXml { path property {value ""} } {
	
	global KPriv
	
	set xpath "[::KEGroups::setGroupsXPath $path]"
	
	if { $value == "" } {
		
		return [$KPriv(xml) set "$xpath/@$property" ]
	} else {
		
		$KPriv(xml) set "$xpath/@$property" "$value"
		#msg "$KPriv(xml) set $xpath/@$property $value"
		return "1"
	}
}

proc ::KEGroups::unsetXml { path } {
	
	global KPriv
	
	set xpath "[::KEGroups::setGroupsXPath $path]"
	
	$KPriv(xml) unset $xpath
}

proc ::KEGroups::insertXml { path id color state type } {
	
	global KPriv
	
	if { $path == "root" } { 
		set xpath "/Kratos_Data/Groups"
	} else {
		set xpath "[::KEGroups::setGroupsXPath $path]"
	}
	
	$KPriv(xml) lappend "$xpath/Group id=\"$id\" color=\"$color\" state=\"$state\" type=\"$type\"" ""
	
	lappend KPriv(groupsId) $id
	
}

#---------------------------------------------------------------------------------------------- 
# Lee el xml y carga el árbol de grupos con todas sus propiedades
#----------------------------------------------------------------------------------------------
proc ::KEGroups::FillTree { } {
	
	global KPriv
	
	#Reseteamos la lista de Id's por si contenía un estado anterior
	set KPriv(groupsId) {}
	
	set T $::KEGroups::TreePath
	
	#Seleccionamos todos los nodos del primer nivel
	set nodes [$KPriv(xml) selectNodes "/Kratos_Data/Groups/Group\[@id\]"]
	
	foreach node $nodes {
		
		if { [$node nodeType] == "ELEMENT_NODE" } {
				
		#Nos guardamos todos los Id
		lappend KPriv(groupsId) [$node getAttribute id ""]
		
		#Añadimos 
		
		# Insertamos cada grupo de 1er nivel en el árbol
		set item [::KEGroups::InsertNewGroup [$node getAttribute id ""] $T [$node getAttribute color ""] [$node getAttribute state ""] [$node getAttribute type ""] "" "root" [$node hasChildNodes] ]
		
		#Seleccionamos sus hijos (2º nivel)
		set nodes2 [$node childNodes]
		foreach node2 $nodes2 {
				
				if { [$node2 nodeType] == "ELEMENT_NODE" } {
				#Agregamos cada elemento de 2º nivel
				lappend KPriv(groupsId) [$node2 getAttribute id ""]		
				set item2 [::KEGroups::InsertNewGroup [$node2 getAttribute id ""] $T [$node2 getAttribute color ""] [$node2 getAttribute state ""] [$node2 getAttribute type ""] "[$node getAttribute id 0]//" "$item" [$node2 hasChildNodes] ]
				
				
				#Seleccionamos los hijos (3º nivel)
				set nodes3 [$node2 childNodes]
				foreach node3 $nodes3 {
						
						if { [$node3 nodeType] == "ELEMENT_NODE" } {
								#msg "..node3: $node3"
						lappend KPriv(groupsId) [$node3 getAttribute id ""]
						set item3 [::KEGroups::InsertNewGroup [$node3 getAttribute id ""] $T [$node3 getAttribute color ""] [$node3 getAttribute state ""] [$node3 getAttribute type ""] "[$node getAttribute id 0]//[$node2 getAttribute id 0]//" "$item2" [$node3 hasChildNodes]]
						#Seleccionamos los hijos (4º nivel)
								set nodes4 [$node3 childNodes]
								
								foreach node4 $nodes4 {
										if { [$node4 nodeType] == "ELEMENT_NODE" } {		
										#msg "...node4: [$node4 getAttribute id 0]"		
										lappend KPriv(groupsId) [$node4 getAttribute id ""]
										set item4 [::KEGroups::InsertNewGroup [$node4 getAttribute id ""] $T [$node4 getAttribute color ""] [$node4 getAttribute state ""] [$node3 getAttribute type ""] "[$node getAttribute id 0]//[$node2 getAttribute id 0]//[$node3 getAttribute id 0]//" "$item3" [$node4 hasChildNodes]]
									   
										#Seleccionamos los hijos (5º nivel)
												set nodes5 [$node4 childNodes]
												foreach node5 $nodes5 {
														if { [$node5 nodeType] == "ELEMENT_NODE" } {			
														lappend KPriv(groupsId) [$node5 getAttribute id ""]
														::KEGroups::InsertNewGroup [$node5 getAttribute id ""] $T [$node5 getAttribute color ""] [$node5 getAttribute state ""] [$node5 getAttribute type ""] "[$node getAttribute id 0]//[$node2 getAttribute id 0]//[$node3 getAttribute id 0]//[$node4 getAttribute id 0]//" "$item4" [$node5 hasChildNodes]
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
	
	return ""
}

proc ::KEGroups::GetAutomaticGroupName { {auto ""} } {
	
	set name ""
	global KPriv
	set i 0
	#foreach grup $KPriv(groupsId) {
	#		#msg "group$i:KPriv $grup"
	#		incr $i
	#}
	if { [llength $KPriv(groupsId)] > 0 } {
		
		for {set i 1} {$i<10000} {incr i} {
		
				set name ${auto}Group${i}
			if { [lsearch -exact $KPriv(groupsId) $name] == -1 } { break }
		}
	} else {
		if { $auto == "" } {
				set name "Group1"
		} else {
			set name "${auto}Group1"
		}
	}
	return $name
}

proc ::KEGroups::CreateNewGroupId {{T ""} {name ""}} {
	
	if { $name == "" } {
		 set name [::KEGroups::GetAutomaticGroupName]
	} else {
		 if { ![::KEGroups::isValidGroupName $name] } {
			 WarnWin [= "Bad group name, start or end by '//' is not allowed"]
			 return ""
		 }
	}
	
	#Obtenemos el item del nodo donde vamos a insertar uno nuevo (si hay varios el primero
	set item [$T selection get 0]
	
	if {$item == "" } {

		set item "root"
	}
	#Obtenemos el path del nodo
	set path [DecodeName [$T item tag names $item]]
	
	#msg "path:$path Num levels:[llength [::KEGroups::split2 $path //] ]"
	
	if { [llength [::KEGroups::split2 $path //] ] >= 5 } {
	
		#Si estamos en el último nivel, lo insertamos en el padre de ese nivel
		set item [$T item parent $item]
		#Y le ponemos su path correcto
		set path [DecodeName [$T item tag names $item]]
		#msg "--->path:$path"
		
	} else {
		#Forzamos el botón de plegar-desplegar por si no lo tenía
		$T item configure $item -button yes
	}
	
	#Preparamos un color aleatorio
	set color [::KEGroups::randomColor]
	
	if { $item == "root" } {
		
		::KEGroups::InsertNewGroup $name $T $color 1 Generic "" "$item" 0
		set path "root"
		
	} else {
		
		::KEGroups::InsertNewGroup $name $T $color 1 Generic "$path//" "$item" 0
	}
		
	::KEGroups::insertXml "$path" $name $color 1 Generic
	
	#Acceder a un elemento por Id's
	#::KEGroups::setXml $fullname color red
	#OLD WAY: $y set {/Kratos_Data/Groups/Group[@id='G1']/Group[@id='G2']/Group[@id='G3']/@color} red
	
	
	
	#GiD_Process Layers New $name escape
	#Bajar al nivel del árbol
	return $name
}
proc ::KEGroups::listReplace {listVariable value {newVal ""}} {
	
	set idx [lsearch -exact $listVariable $value]
	
	if { $idx != -1 } {
		
		if { $newVal == "" } { 
			return [lreplace $listVariable $idx $idx]
		} else {
			return [lreplace $listVariable $idx $idx $newVal]
		}
	} else {
		return $listVariable
	}
}

proc ::KEGroups::DeleteGroupsId { T } {
	
	global KPriv
	
	set items [$T selection get]
	
	if {[llength $items] > 0 } {
		
		# Para avisar de qué items se van a borrar
		set tuttoItem {}
		foreach it $items {
				lappend tuttoItem "[lindex [$T item text $it] 0]"
		}
		set aviso "The selected items and all their descendants are going to be removed:\n\n$tuttoItem"
		set confirmado [::WinUtils::confirmBox "." "$aviso"]
		if { $confirmado == "ok" } {
				 
				 #Buscamos todos los descendientes del item a eliminar
				set completList {}
				foreach item $items {
				
						lappend completList $item
						foreach i [$T item descendants $item] {
								set idx [lsearch -exact $items $i]
								if { $idx == -1 } {
										lappend completList $i		
								}
						}
				}
				
				for {set i 0} { $i < [llength $completList] } {incr i} {
						
						#Desasigna de la gemoetría cada item seleccionado
					set GroupId [$T item text [lindex $completList $i] 0]
						::KEGroups::UnAssignCondition $GroupId
						
						#Eliminamos el elemento de la lista de ID's
						set KPriv(groupsId) [::KEGroups::listReplace $KPriv(groupsId) $GroupId]
						
						#Elimina el frame asociado
						set GroupId [$T  item text [lindex $completList $i] 0]
						if { [winfo exists ${T}.f$GroupId] } {
								destroy ${T}.f$GroupId
						}
					   
							::KMProps::deleteGroups $GroupId
					   
				}
				
				for {set i 0} { $i < [llength $items] } {incr i} {
						
						#Eliminamos cada item seleccionado, con un catch por si se intenta borrar alguno que ya no existe				
						catch { ::KEGroups::unsetXml [DecodeName [$T item tag names [lindex $items $i]]]}
						
						
						#Elimina el item del árbol
						::KMProps::deleteItem $T [lindex $items $i]
					
				}
				#Si está la ventana de listEntities abierta la refrescamos
				if { [winfo exists ${::KEGroups::WinPath}.listEntities] } {
					::LEntities::InitBaseWindow
				}
		}		
	} else {
		WarnWin [= "No group selected"]
	}
}
	
proc ::KEGroups::isValidGroupName {name} {
	
	if { [string range $name end-1 end] == "//" } {
	return 0
	} elseif { [string range $name 0 1] == "//" } {
	return 0
	} elseif { $name == "" } {
	return 0
	}
	return 1
}

proc ::KEGroups::ClickTree { x y T } {
	
	set info [$T identify $x $y]
	
	if { [lindex $info 0] == "item" && [llength $info] >= 4 } {
	
		set item [lindex $info 1]
		set col [lindex $info 3]
		
		set fullname [DecodeName [$T item tag names $item]]
		#msg "Fullname:$fullname Item:$item Col:$col"
		
	} elseif { [lindex $info 0] == "header" && [lindex $info 1] == "0" } {
	
		if { [$T column cget C0 -arrow] == "up" } {
			$T column configure C0 -arrow down
			$T item sort 0 -dictionary -increasing
		} else {
			$T column configure C0 -arrow up
			$T item sort 0 -dictionary -decreasing
		}
		#msg "Info:$info"
		return ""
	} else {
		return ""
	}
	
	if { $col != 0 } {
		if { ![$T selection includes $item] } {
			$T selection clear
			$T selection add $item
		}
	}
	if { $col == 0 } {
		#SetLayersTo TOUSE $T
	} elseif { $col == 1 } {
					
		set parent [winfo parent $T]
   
		#Obtenemos el color del item seleccionado
		set cur_col [lindex [::KEGroups::setXml $fullname color] 0]
		
		set color [GIDChooseColor $parent.selectcolor -title [= "Select color"] -color $cur_col]
		if { $color != "" } {
			#Cambiamos el color de todos los item seleccionados
			set items [$T selection get]
			for {set i 0} { $i < [llength $items] } {incr i} {
						
					$T item element configure [lindex $items $i] C1 elemRectColor -fill $color
					set fullname [DecodeName [$T item tag names [lindex $items $i]]]		
					::KEGroups::setXml $fullname color $color
				}
		}	

	} elseif { $col == 2 } {
		
		if { [$T item element cget $item C$col elemImgAny -image] == [GetImage layer_off.gif] } {
			set estado 1
		} else {
		set estado 0
		}
		
		#Cambiamos la imagen (bombilla) de todos los item seleccionados
		set items [$T selection get]
		set grupos {}
		for {set i 0} { $i < [llength $items] } {incr i} {
		
		set item [lindex $items $i]
			set fullname [DecodeName [$T item tag names $item]]
				
				#Nos guardamos los idGrupo de la selección
				lappend grupos [lindex [$T item text $item] 0]
				
				if { $estado == 1 } {
					$T item element configure $item C$col elemImgAny -image [GetImage layer_on.gif]
					::KEGroups::setXml $fullname state "1"
					
				} else {
				$T item element configure $item C$col elemImgAny -image [GetImage layer_off.gif]
				::KEGroups::setXml $fullname state "0"
				
				}
		}
		#Escondemos o mostramos las entidades de los grupos seleccionados
		if { [llength $grupos] } {
			
			if {$estado == 1 } {
				
				# Para activar, es necesario volver a mostrar todas y desactivar todas las ocultas
				GiD_Process Layers BringToFrontAll
				#msg "GiD_Process Layers BringToFrontAll"
				#Para la activación, lo q necesitamos son los grupos no activos
				set grupos [::KEGroups::getDesactivatedGroups $T]
			}
			
			 # Get the project view mode
    		set PState [GiD_Info Project ViewMode]
			
			foreach groupId $grupos {
				# Switch state
				if {($PState == "GEOMETRYUSE") } {
					
					foreach entity [list point line surface volume] {
						set entityList [::KEGroups::getGroupGiDEntities $groupId $entity]
						if { [llength $entityList]} {
							GiD_Process Layers SendToBack $entity {*}$entityList escape
						}
					}
						
				} elseif {($PState == "MESHUSE") } {
					
					foreach entity [list nodes elements] {
						set entityList [::KEGroups::getGroupGiDEntities $groupId $entity]
						if { [llength $entityList]} {
							GiD_Process Layers SendToBack $entity {*}$entityList escape
						}
					}
					
					
				}
			}
		#Redibuja las entidades en pantalla	
		GiD_Redraw
		}
				
	} elseif { $col == 3 } {
		#No llega el evento porque lo absorve el combo
	}
	
	if { $col != 0 } {
	return -code break
	}
	return ""
	
}


proc ::KEGroups::SetGroupsToRename { T item newtext } {
	
	global KPriv
	
	if { $newtext == ""} {  
		WarnWin [= "You can not choose an empty group name."]
		return  
	}
	if { $item == 0 } {
		WarnWin [= "Root folder can't be edited"]
		return
	}
	#Validamos q el nombre no tenga carácteres que vulneran la seguridad y quitamos espacios
	set newtext [::KUtils::parseTreeStr $newtext]
	if { $newtext == -1 } {
		WarnWin [= "You can't use some reservate chars like:\n  :   /   $   .   \\  %  "]
		return
	}
	
	#Si han elegido el mismo nombre ya no es necesario hacer nada más
	set oldgroupid [lindex [$T item text $item] 0]
	if {$oldgroupid == $newtext } {
	    return
	}
	
	#Controlamos q el nombre no esté ya en el árbol (a no ser q no se haya cambiado)
	if { $newtext in $KPriv(groupsId) } {
		WarnWin [= "The group name '%s' already exist.\nChoose another, please." $newtext]
		return
	}
  
	# Si el nombre es correcto procedemos a renombrar las asignaciones en el arbol Model 
	# si es necesario (pasamos primero el nombre antiguo y luego el nuevo
	::KMProps::deleteGroups [$T item text $item 0] $newtext
	
	#Guardamos el nivel del item a renombrar
	set fullname [DecodeName [$T item tag names $item]]
	set splitted [::KEGroups::split2 $fullname //]
	set whereRename [llength $splitted]
	
	#Recorremos toda su familia		cambiando cada path
	set items [$T item descendants $item]
	
	foreach i $items {
		
		set fullname [DecodeName [$T item tag names $i]]
		set splitted [::KEGroups::split2 $fullname //]
		
		#Sustituimos la posición apropiada del path por el nuevo nombre
		set splitted [::KEGroups::listReplace $splitted [lindex $splitted [expr $whereRename - 1]] $newtext]
		
		#Reconstruimos el nuevo path
		set fullNewName ""
		foreach iSplit $splitted {
				
				set fullNewName "$fullNewName$iSplit//"
		}
		set fullNewName [string range $fullNewName 0 end-2]
		
		#Cambiamos la etiqueta (el path) en el arbol
		$T item tag remove $i [list names [$T item tag names $i]]
		$T item tag add $i [EncodeName $fullNewName]
	}
	
	#Cambiamos el item a renombrar 
	set fullname [DecodeName [$T item tag names $item]]
	
	::KEGroups::editTag $T $item $fullname $newtext
	
	# Update the link to the GiD condition properties for this group identifier
	::KEGroups::RenameGroupIdGiDCond $oldgroupid $newtext
	
	#Si está la ventana de listEntities abierta la refrescamos
	if { [winfo exists ${::KEGroups::WinPath}.listEntities] } {
		::LEntities::InitBaseWindow
	}

	return ""
}

proc ::KEGroups::RenameGroupIdGiDCond {oldgroupid newgroupid} {
    # Update the link to the GiD condition properties for this group identifier

    # Disable graphics
    ::GidUtils::DisableGraphics
    # Store the freeze layers
    set flayerslist [::KEGroups::UnFreezeLayers]
    # Get the project view mode
    set PState [GiD_Info Project ViewMode]
       
    # For each geometry o mesh mode
    foreach gmid [list geometry mesh] {
	
	# Switch state
	if {($PState == "GEOMETRYUSE") && ($gmid == "mesh")} {
		GiD_Process MEscape Meshing MeshView
	} elseif {($PState == "MESHUSE") && ($gmid == "geometry")} {
		GiD_Process MEscape Geometry
	}

	# For each GiD group entities
	foreach EntityId [list point line surface volume] {
	    set EntityList [list]
	    foreach CondProp [GiD_Info conditions ${EntityId}_groups $gmid] {
		lassign $CondProp - CEId - CGroupId
		# Update the entities list
		if { $CGroupId eq $oldgroupid } {
		    lappend EntityList $CEId
		}
	    }
	    #WarnWinText "EntityList:$EntityList"
	    if {[llength $EntityList]} {
		# Rename the GiD conditions
		# UnAssign condition
		eval [list GiD_Process MEscape Data Conditions AssignCond \
			  ${EntityId}_groups UnAssign Field groupid $oldgroupid] $EntityList escape
		# Reassign condition with the new name
		GiD_Process MEscape Data Conditions AssignCond ${EntityId}_groups NoRepeatField groupid
		eval [list GiD_Process change $newgroupid] $EntityList escape
	    }
	}
	# Switch state
	if {($PState == "GEOMETRYUSE") && ($gmid == "mesh")} {
	    GiD_Process MEscape Geometry
	} elseif { ($PState == "MESHUSE") && ($gmid == "geometry") } {
	    GiD_Process MEscape Meshing MeshView				
	}
    }
    # Restore the layer state
    ::KEGroups::FreezeLayers $flayerslist
    # Enable graphics
    ::GidUtils::EnableGraphics 
}

#
# Renombra un item, modificando su path y reconstruyendo el combo
#
proc ::KEGroups::editTag { T item fullname newtext } {
	
	global KPriv		
	
	set parts [::KEGroups::split2 $fullname //]
	lset parts end $newtext
	set newPath [join $parts //]

	#Renombra en la lista de ID's
	set idItem [lindex [$T item text $item] 0]
	set KPriv(groupsId) [::KEGroups::listReplace $KPriv(groupsId) $idItem $newtext]
	
	#Reconstruye el combo
	set type [set ::KEGroups::type$idItem]
	
	if { [winfo exists $T.f$idItem] } {
		destroy $T.f$idItem
	}
	set f [::KEGroups::buildCombo $T $item $newtext $type $newPath]	
	$T item style set $item C3 styFrame
	$T item element configure $item C3 eWindow -window $f	
		
		
	#Cambiar nombre en el árbol
	$T item tag remove $item [list names [$T item tag names $item]]
	
	$T item tag add $item [EncodeName $newPath]
	$T item element configure $item C0 elemTxtRead -text $newtext
	
	#Cambiar nombre en el XML
	::KEGroups::setXml $fullname id $newtext
	
	return $newPath
}					

#separator can be a multicharacter, like //
proc ::KEGroups::split2 { x separator } {
	
	set splitchar \uFFFF ;#forbidden unicode character, x must never contain it
	return [split [string map "$separator $splitchar" $x] $splitchar]
}

proc ::KEGroups::CreateTreeProperties {w} {
	
	variable SystemHighlight
	variable SystemHighlightText

	# Scrollbars
	set vsb $w.vsb1
	set hsb $w.hsb1

	# Create the treectrl and set the scrollbars
	set T [treectrl $w.t -xscrollcommand [list $hsb set] -yscrollcommand [list $vsb set]]
	ttk::scrollbar $vsb -orient vertical   -command [list $T yview]
	ttk::scrollbar $hsb -orient horizontal -command [list $T xview]
	
	# Set the height
	set height [font metrics [$T cget -font] -linespace]
	if {$height < 22} {
	 set height 22
	}

	# Configure the treectrl
	$T configure -indent 15 -itemheight $height -selectmode extended \
	 -showroot 0 -showrootbutton 0 -showbuttons 1 -showlines 1 \
	 -highlightthickness 0 -borderwidth 0 \
	 -xscrollincrement 20 -yscrollincrement 20
	
	# Create the column identifier list
	set collistid [list [= "Group Id"] [= "Color"] [= "On/Off"] [= "Type"]]
	
	set i 0
	foreach cid $collistid {
	 $T column create -text $cid -tags C$i -weight 0
	 incr i
	}
	# Configure the column weight and arrow
	$T column configure C0 -weight 1 -arrow up

	# Configure the column that have the tree
	$T configure -treecolumn C0

	# Create elements
	$T element create elemImgAny image
	$T element create elemTxtRead text -fill [list $SystemHighlightText {selected focus}] -lines 1
	$T element create elemRectSel rect -fill [list $SystemHighlight {selected focus} gray {selected !focus}] -showfocus yes
	$T element create elemRectColor rect -width 30
	$T element create eWindow window
	
	
	# Create styles using the elements
	set S [$T style create styAnyRead]
	$T style elements $S {elemRectSel elemTxtRead elemImgAny}
	$T style layout $S elemImgAny -expand ns
	$T style layout $S elemTxtRead -padx 4 -expand ns -squeeze x
	$T style layout $S elemRectSel -union [list elemTxtRead] -iexpand ns -ipadx 2
	
	set S [$T style create styRectColor]
	$T style elements $S {elemRectColor elemTxtRead}
	$T style layout $S elemRectColor -union [list elemTxtRead] -iexpand ns -ipadx 2 -pady 1

	set S [$T style create styAnyImage]
	$T style elements $S {elemImgAny}
	$T style layout $S elemImgAny -expand ns -pady 1
	
	 set S [$T style create styFrame -orient horizontal]
	$T style elements $S {eWindow}
	$T style layout $S eWindow -iexpand x -squeeze x -padx {0 1} -pady {2 2}
	
	# Items
	set item root
	$T item configure $item -button yes
	$T item style set $item C0 styAnyRead
	$T item element configure $item C0 elemTxtRead -text [= "Root"]

	# List of lists: {column style element ...} specifying text elements
	# the user can edit
	TreeCtrl::SetEditable $T {
	 {C0 styAnyRead elemTxtRead}
	}
	
	# List of lists: {column style element ...} specifying elements
	# the user can click on or select with the selection rectangle
	TreeCtrl::SetSensitive $T {
	 {C0 styAnyRead elemTxtRead}
	}

	# List of lists: {column style element ...} specifying elements
	# added to the drag image when dragging selected items
	TreeCtrl::SetDragImage $T {
	}

	# Some notify install
	#$T notify install <Drag-receive>
	$T notify install <Edit-accept>

	# Notify bind
	# TODO
	  #$T notify bind DragTag <Drag-receive> { ::KEGroups::ReceiveDragGroups  %T %l %I }
	  $T notify bind EditTag <Edit-accept> { ::KEGroups::SetGroupsToRename %T %I %t }

	  bind $T <Button-1> [list ::KEGroups::ClickTree %x %y $T]
#	 bind $T <Double-Button-1> [list DoubleClickTableList_New %x %y $T]
#	 bind $T <Return> [list SetLayersTo TOUSE $T]
#	 bind $T <Key-Delete> [list SetLayersToDelete $T]
#	 bind $T <Alt_L> [list InvertSelectionTableList $T]
#	 bind $T <Alt_R> [list InvertSelectionTableList $T]
#	 bind $T <Meta_L> [list InvertSelectionTableList $T]
#	 bind $T <Meta_R> [list InvertSelectionTableList $T]
	  bind $T <F2> [list ::KEGroups::BeginEditGroups $T]

	bind $T <Button-3> "[list ::KEGroups::MenuContextualGroup %W %x %y] ; break"

	# Grid the tree
	grid $T $vsb -sticky nsew
	grid configure $vsb -sticky ns
	grid $hsb -sticky ew
	grid remove $vsb $hsb
	bind $T <Configure> [list ConfigureListScrollbars $T $hsb $vsb]
	
	grid columnconfigure $w 0 -weight 1
	grid rowconfigure $w 0 -weight 1

	bindtags $T [list $T TreeCtrlFileList TreeCtrl [winfo toplevel $T] all]

	return $T
}

proc ::KEGroups::BeginEditGroups { T } {
	
	set I [$T item id active]
	set C 0
	set E elemTxtRead
	::TreeCtrl::FileListEdit $T $I $C $E
}

proc ::KEGroups::ReceiveDragGroups { T dragged_list dst } {
	
	#set dstname [DecodeName [$T item tag names $dst]]
	#if { [IsItemFolderOfLayers $T $dst] } {
	##folder of layers
	#set cmd [GetOldAndNewLayernamesRecursive $T $dragged_list $dstname ""]
	#GiD_Process Layers {*}$cmd escape
	#} else {
	#WarnWin "$dstname is not a folder"
	#}
}

proc ::KEGroups::WriteGeomToVar { w what geomname {InitComm ""}} {
	
	variable WinLayout
	global GidPriv

	set trans 1
	update idletasks
	if { $WinLayout == "RIGHT" } {
	# maintain all dimensions except the new width
	set prevgeom [lindex $GidPriv($geomname) 1]
	foreach {width height x y} [split $prevgeom x+] break
	set width [winfo width $w]
	} else {
	# "OUTSIDE"
	set width [winfo width $w]
	set height [winfo height $w]
	set x [winfo x $w]
	set y [winfo y $w]
	}
	set geom ${width}x${height}+${x}+${y}
	set GidPriv($geomname) "$what $geom $trans $InitComm"
	#  msg $GidPriv($geomname)
}

proc ::KEGroups::OpenWindowOutside { w } {
 
	if { [winfo exists $w] } {
	::KEGroups::CloseWindowOutside $w
	}

	# Init the window
	set title [= "Entities group editor"]
	InitWindow $w $title KEGroupsWindowGeom ::KEGroups::InitBaseWindow
	
	return $w
}

proc ::KEGroups::CloseWindowOutside { w } {
	destroy $w
}

proc ::KEGroups::IsPreprocessMode {} {
	# Check for preprocess mode
	if { [lsearch -exact {GEOMETRYUSE MESHUSE} [GiD_Info Project ViewMode]] == -1 } {
	return 0 
	} 
	return 1
}

proc ::KEGroups::InsertNewGroup { groupName T color {state 1} {type "Generic"} {parent ""} {parentitem root} {childs true} } {
	
	if { $parent != "" } {
		set fullname $parent$groupName
	} else {
		set fullname $groupName
	}
	
	if { $childs } {
		
		set item [$T item create -button yes -tags [EncodeName $fullname]]
		$T item lastchild $parentitem $item
		$T item style set $item C0 styAnyRead
		$T item element configure $item C0 elemTxtRead -text "$groupName"
		
	} else {
		
		set item [$T item create -button no -tags [EncodeName $fullname]]
		$T item lastchild $parentitem $item
		$T item style set $item C0 styAnyRead
		$T item element configure $item C0 elemTxtRead -text $groupName
	}

	
	  
	$T item style set $item C1 styRectColor
	$T item element configure $item C1 elemRectColor -fill $color

	$T item style set $item C2 styAnyImage
	if { $state == 1 } {
		$T item element configure $item C2 elemImgAny -image [GetImage layer_on.gif]
	} else {
		$T item element configure $item C2 elemImgAny -image [GetImage layer_off.gif]
	}
	
	#$T item style set $item C3 styAnyRead
	#$T item element configure $item C3 elemTxtRead -text "$type"
	
	#---------------------------#---------------------------#
	# Configurar frame para comboBox
	#---------------------------#---------------------------#
	
	#set id [::KEGroups::setXml $fullname id]
	#msg "id:$id  -  item$item: $T.f$id"
	
	set f [::KEGroups::buildCombo $T $item $groupName $type $fullname]
	
	$T item style set $item C3 styFrame
	$T item element configure $item C3 eWindow -window $f
   
	return $item
}

proc ::KEGroups::buildCombo {T item id type fullname } {
	
	set bg "#FFFFCC"
	
	set f [frame $T.f$id -borderwidth 0 -background $bg]
	   
	ttk::combobox $f.mb1 -values [::KEGroups::getGroupsType] -state readonly -width 9 -height 10 -textvariable "::KEGroups::type$id"
	
	set selected [::KEGroups::getSelected $type]
	
	$f.mb1 current $selected
	
	pack $f.mb1 -side left -padx {1 1}
	
	bind $f.mb1 <<ComboboxSelected>> "::KEGroups::changeTypeEvent $T $fullname $item"
	
	return $f
}

proc ::KEGroups::changeTypeEvent {T fullname item {changeTo ""} } {
	
	if { $changeTo == "" } {
		
		set id [lindex [$T item text $item] 0]
		set tipo [ set ::KEGroups::type${id}]
	} else {
		set tipo "$changeTo"
	}
	
	::KEGroups::setXml $fullname type $tipo
	
	#Cambiamos también el tipo de todos los item seleccionados
	set items [$T selection get]
	#msg "items:$items"
	for {set i 0} { $i < [llength $items] } {incr i} {
		
		set item [lindex $items $i]
		
		set fullname [DecodeName [$T item tag names $item]]
		set id [::KEGroups::setXml $fullname id]
		
		#Cambiamos el valor en el combo
		set ::KEGroups::type$id $tipo
		
		#Actualizamos el valor en el xml
		::KEGroups::setXml $fullname type $tipo
	}
}

proc ::KEGroups::getSelected { type } {

	set typeList [::KEGroups::getGroupsType]
	
	set i 0
	foreach iTipo $typeList {
		#msg "itip$i:$iTipo"
		if {$iTipo == $type} {
				return $i
		}
		set i [expr $i + 1]
	}
	return [expr [llength $typeList] - 1]
}

proc ::KEGroups::randomColor { } {

	set r [expr {int (255 * rand())}]
	set g [expr {int (255 * rand())}]
	set b [expr {int (255 * rand())}]		
   
	return [format "\#%02x%02x%02x" $r $g $b]
}

proc ::KEGroups::MenuContextualGroup { T x y } {
	set w $T.menucontextualgroup
	if { [winfo exists $w] } {
	destroy $w
	}

	menu $w

	$w add command -label [= "New Group#C#layer"] -command [list ::KEGroups::CreateNewGroupId $T]

	set nItems [llength [$T selection get]]
	
	if { $nItems == 1 } {
	$w add command -label [= "Rename"] -command [list ::KEGroups::BeginEditGroups $T] -state normal
	} else {
	$w add command -label [= "Rename"] -state disabled
	}
	if { $nItems > 0 } {
	$w add command -label [= "Delete"] -command [list ::KEGroups::DeleteGroupsId $T] -state normal
	} else {
	$w add command -label [= "Delete"] -state disabled
	}
	$w add separator

	if { $nItems > 0 } {
	set state normal
	} else {
	set state disabled
	}
	$w add command -label [= "Color"]... -command [list ::KEGroups::actionTree $T COLOR] -state $state
	$w add command -label [= "On#C#verb"] -command [list ::KEGroups::actionTree $T ON] -state $state
	$w add command -label [= "Off#C#verb"] -command [list ::KEGroups::actionTree $T OFF] -state $state
	
	  
	#Para añadir los tipos al menú contextual (no está probado)
	$w add separator
	if { $nItems > 0 } {
		$w add cascade -label [= "Type"] -menu $w.type -state normal
	} else {
		$w add cascade -label [= "Type"] -menu $w.type -state disabled
	}
	set typesList [::KEGroups::getGroupsType]
	
	menu $w.type
	foreach i $typesList {
		
		$w.type add command -label [= "%s" $i] -hidemargin 1 -state $state -command "::KEGroups::actionTree $T TYPE $i"
		
	}
	
	
	# ASIGNAR ENTIDADES A GRUPO
	$w add separator
	if { $nItems > 0 } {
		$w add cascade -label [= "Assign Entities"] -menu $w.assign -state normal
	} else {
		$w add cascade -label [= "Assign Entities"] -menu $w.assign -state disabled
	}
	menu $w.assign -postcommand "::KEGroups::rebuildMenu $w.assign"
	
	#DESASIGNAR
	$w add command -label [= "Unassign"] -command [list ::KEGroups::unAssignEntities] -state $state
	
	set x [expr [winfo rootx $T]+$x+2]
	set y [expr [winfo rooty $T]+$y]
	GiD_PopupMenu $w $x $y
}

proc ::KEGroups::actionTree { T action {tipo "Generic"} } {
	
	set items [$T selection get]
	
	if { $action == "ON" || $action == "OFF"} {
		
		if { [$T item element cget [lindex $items 0] C2 elemImgAny -image] == [GetImage layer_off.gif] } {
			set estado 1
		} else {
		set estado 0
		}		
	} elseif { $action == "COLOR" } {
		set fullname [DecodeName [$T item tag names [lindex $items 0]]]
		set cur_col [lindex [::KEGroups::setXml $fullname color] 0]
		set color [GIDChooseColor [winfo parent $T].selectcolor -title [= "Select color"] -color $cur_col]				
	
	}		
	
	#Obtenemos todos los items seleccionados		
	
	for {set i 0} { $i < [llength $items] } {incr i} {
		
		set fullname [DecodeName [$T item tag names [lindex $items $i]]]   
		
		if { $action == "ON" } {
				
				$T item element configure [lindex $items $i] C2 elemImgAny -image [GetImage layer_on.gif]
				::KEGroups::setXml $fullname state "1"
				
		} elseif { $action == "OFF" } {
				
				$T item element configure [lindex $items $i] C2 elemImgAny -image [GetImage layer_off.gif]
				::KEGroups::setXml $fullname state "0"
				
		} elseif { $action == "COLOR" } {
				
				$T item element configure [lindex $items $i] C1 elemRectColor -fill $color
			::KEGroups::setXml $fullname color $color
		} elseif { $action == "TYPE" } { 
				
				::KEGroups::changeTypeEvent $T $fullname [lindex $items $i] $tipo
		}
	}		
	
	return ""
	
}

#
# ACTUALMENTE NO SE UTILIZA
# Obtiene todos los GroupId's sin insertar nada en el arbol
# (si llega un groupId, devuelve toda su descendencia)
#
proc ::KEGroups::getXmlGroupsId { {groupId ""}} {
	
	global KPriv
	
	#Reseteamos la lista de Id's por si contenía un estado anterior
	set KPriv(groupsId) {}
	# Si llega un groupId
	set childNodes {}
	
	#Seleccionamos todos los nodos del primer nivel
	set nodes [$KPriv(xml) selectNodes "/Kratos_Data/Groups/Group\[@id\]"]
	
	foreach node $nodes {
		
		#Nos guardamos todos los Id
		#if { [llength [set childNodes [::KEGroups::addGroupId $node $groupId]]] == 0 } {
		#		return $childNodes
		#}
		lappend KPriv(groupsId) [$node getAttribute id ""]		
		#Seleccionamos sus hijos (2º nivel)
		set nodes2 [$node childNodes]
		foreach node2 $nodes2 {
		#Agregamos cada elemento de 2º nivel
		lappend KPriv(groupsId) [$node2 getAttribute id ""]		
		#Seleccionamos los hijos (3º nivel)
		set nodes3 [$node2 childNodes]
		foreach node3 $nodes3 {
			lappend KPriv(groupsId) [$node3 getAttribute id ""]
			#Seleccionamos los hijos (4º nivel)
				set nodes4 [$node3 childNodes]
				foreach node4 $nodes4 {
				lappend KPriv(groupsId) [$node4 getAttribute id ""]
				#Seleccionamos los hijos (5º nivel)
						set nodes5 [$node4 childNodes]
						foreach node5 $nodes5 {
								lappend KPriv(groupsId) [$node5 getAttribute id ""]
						}
						}
		}
		}
	}
	return ""
}

proc ::KEGroups::addGroupId { node groupId } {
	
	global KPriv
	
	set id [$node getAttribute id ""]
	lappend KPriv(groupsId) $id
	
	if { $groupId == $id } {
		
		set childs {$groupId}
		set descendants [$node descendant all]
		
		foreach des $descendants {
				lappend childs [$des getAttribute id ""]
				#msg [$des getAttribute id ""]
		}
		return $childs
	}
	return {}
}

proc ::KEGroups::getDesactivatedGroups { T } {
	
	set grupos {}
	set descendants [$T item descendants root]
		
	foreach item $descendants {
		
		if { [$T item element cget $item C2 elemImgAny -image] == [GetImage layer_off.gif] } {
			lappend grupos [lindex [$T item text $item] 0]   
		}   
	}
	return $grupos
}
##############################################################################################
#
#        NAME: kegroups.tcl
#
#        PURPOSE: Utilities procedures to work with the Kratos entities groups editor
#
#        QUANTECH ATZ-DEVELOPMENT DEPARTMENT
#
#        AUTHOR : G. Socorro
#
#        CREATED AT: 01/11/09
#
#        HISTORY:
#
#        2.6- 19/06/13- G. Socorro, modify/update some procedures (BorraGrupo, CheckGroup, RenombraGrupo, ValidateNewName)
#        2.5- 18/06/13- G. Socorro, delete the use of the proc kipt::NewGiDGroups (delete the call to the compass groups => Cond_Groups)
#        2.4- 12/04/13- G. Socorro, correct a bug in the proc getGroupGiDEntities for the special case of faces, add -count in the proc getGroupGiDEntitiesNew
#        2.3- 13/12/12- J. Garate,  Corrected ::KEGroups::getGroupGiDEntities bug for Model Validation and Transfrer CondGroups to GiD Groups
#        2.2- 28/11/12- J. Garate,  Corrected bug when transferring old groups to new gid groups, erasing old Cond
#        2.1- 07/11/12- J. Garate,  ::KEGroups::GroupsToXml is ready to accept New GiD Groups, and modificate the .spd field "Groups modeltype"
#        2.0- 22/10/12- J. Garate,  Corrected bugs on New GiD_Groups management
#        1.9- 10/10/12- J. Garate,  Full Adaptation to New GiD_Groups, including corrections and new functions
#        1.8- 27/05/12- J. Garate,  Creation of ::KEGroups::listReplace function
#        1.7- 07/05/12- J. Garate,  update renaming groups restrictions
#        1.6- 04/05/12- G. Socorro, update some procedures
#        1.5- 03/05/12- J. Garate,  GiD Groups transfer to .spd // W Child
#        1.4- 03/05/12- J. Garate,  GiD Groups transfer to .spd
#        1.3- 26/04/12- G. Socorro, change GiD_Groups by Cond_Groups
#        1.2- 26/03/12- J. Gárate,  Cambio de iconos para AutoGroup
#        1.1- 26/03/12- J. Gárate,  Renombrado de grupos. Ventana de error en el renombrado.
#        1.0- 22/03/12- J. Gárate,  Cambio a funciones públicas de los grupos de GiD. Borrado de funciones de Grupos antiguas, ahora están en OLDKEGROUPS.TCL
#        0.9- 20/03/12- J. Garate,  Renombrado de grupos y actualizacion del arbol tras el evento de renombrado
#        0.8- 19/03/12- J. Garate,  Borrado de grupos y actualizacion del arbol tras el evento de borrado
#        0.7- 12/03/12- J. Garate   AutoGroup arreglado para INSIDE Window
#        0.6- 12/03/12- J. Garate   Adaptacion a los nuevos Grupos de GiD. Pendiente: Eliminar funciones antiguas.
#        0.5- 07/02/12- J. Garate   Actualizada la funcion ::KEGroups::getGroupGiDEntities
#        0.4- 22/06/11- G. Socorro, delete snit, tdom and xmlstruct from the package require
#        0.3- 13/05/10- G. Socorro, add the procedure ::KEGroups::RenameGroupIdGiDCond to rename the group identifier in the GiD condition database
#        0.2- 19/03/10- Luis CA,    Reparar acciones en el arbol (delete masivo, rename, ...) y aumentar el número de niveles a 5
#        0.1- 01/11/09- G. Socorro, create a base source code from the GiD layer.tcl script
#         
#
##############################################################################################

package require treectrl
package require tooltip
package provide KEGroups 1.0 

# Create a base namespace KEGroups
namespace eval ::KEGroups:: {
    
    # Path of the base window 
    variable WinPath ".gid.kegroups"
    variable WinLayout 
    variable SystemHighlight
    variable SystemHighlightText
    # The tree path
    #variable TreePath
}

proc ::KEGroups::GroupsSelectionAssign {entity} {
    
	variable TreePath
    
	set item [$TreePath selection get 0]
    
	if {$item != ""} {
	
		set GroupId [$TreePath item text $item 0]
	
		::KEGroups::SelectionAssign $entity $GroupId $::KEGroups::WinPath
	
	} else {
	
		msg "No group selected"
	}
}

proc ::KEGroups::SelectionAssign { entity GroupId WinPath } {

    # Get the condition identifier
    
    # Mapeamos el tipo de entity
    if {$entity eq "point"} {
	set entity "Points"
    } elseif {$entity eq "line"} {
	set entity "Lines"
    } elseif {$entity eq "surface"} {
	set entity "Surfaces"
    } elseif {$entity eq "volume"} {
	set entity "Volumes"
    } elseif {$entity eq "node"} {
	set entity "Nodes"
    } elseif {$entity eq "all"} {
	if {[GiD_Info Project ViewMode] eq "MESHUSE"} {
	    set entity "All_Mesh_Types"
	} else {
	    set entity "All_Geom_Types"
	}
    } else {
	set entity "Elements"
    } 

    #set OldSmallWinSelecting [GiD_Info variable SmallWinSelecting]
    #if {$OldSmallWinSelecting == 0 } {
    #set SmallWinSelecting 1
    #GiD_Set SmallWinSelecting $SmallWinSelecting
    #} else {
    #set SmallWinSelecting 1
    #}
    
    # Try to assign the entities
    
    GiD_Process MEscape
    
    GiD_Process Utilities EntitiesGroups Assign $GroupId $entity
    FinishButton $WinPath $WinPath.bBottomOk [= "Press 'Finish' to stop the entities selection"] "" disableall [GiD_Set SmallWinSelecting]        
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

proc ::KEGroups::GetStateGeoMesh { entity } {
    switch $entity {   
	"point" { return "geometry" }
	"line" { return "geometry" }
	"surface" { return "geometry" }
	"volume" { return "geometry" }
	"nodes" { return "mesh" }
	"element" { return "mesh" }
	"faces" { return "mesh" }
    }
}

proc ::KEGroups::getGroupGiDEntities {groupId {givenEntity "point"} {action ""}} {
    
    # Update the link to the GiD condition properties for this group identifier
    
    # wa "groupId:$groupId givenEntity:$givenEntity action:$action"
    # Disable graphics
    ::GidUtils::DisableGraphics
    # Store the freeze layers
    set flayerslist [::KEGroups::UnFreezeLayers]
    
    # For each GiD group entities
    set EntityIdList {}
    
    
    if { $givenEntity eq "ALL" } {
	
	set entities [list point line surface volume node element]
	
    } elseif { $givenEntity == "element" } {
	
	set entities [list line surface volume]
	
    } elseif { $givenEntity == "nodes" } {
	
	set entities "point"
	
    } elseif { $givenEntity == "faces" } {
	
	set entities [list line surface]
	
    } else {
	
	set entities $givenEntity
    }
    
    # Special case of faces
    set elemIdList ""
    set faceIdList ""

    foreach entity $entities {
	if {$givenEntity eq "ALL" } {
	    set gmid [::KEGroups::GetStateGeoMesh $entity]
	} else {
	    set gmid [::KEGroups::GetStateGeoMesh $givenEntity]
	}
	foreach CondProp [GiD_Info conditions ${entity}_groups $gmid] {
	    lassign $CondProp CId CEId - CGroupId
	    # msg "Entity = $entity"
	    # msg "CId = $CId"
	    # msg "CEId = $CEId"
	    # msg "entity : $entity // Cprop: $CondProp // gmid: $gmid"
	    
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
		if {$givenEntity eq "faces"} {
		    if {($CId ne "E")&&($CId ne "N")} {
		               lappend elemIdList $EID
		        lappend faceIdList $FID
		    }
		} else {
		    if {($CId == "E")||($CId == "N")} {
		        lappend EntityIdList $CEId
		    } 
		}
	    }
	}
    }
    
    # Join element and face list for the special case of faces
    if {$givenEntity eq "faces"} {
	if {([llength $elemIdList]) && ([llength $faceIdList])} {
	    set EntityIdList [list $elemIdList $faceIdList]
	}
    }

    # Restore the layer state
    ::KEGroups::FreezeLayers $flayerslist
    
    # Enable graphics
    ::GidUtils::EnableGraphics 
    
    # wa "EntityIdList:$EntityIdList"
    if {$action == "hasEntities" } {
	return 0
    }
    return $EntityIdList
}

proc ::KEGroups::getGroupGiDEntitiesNew {groupId {action ""}} {
    
    # Update the link to the GiD condition properties for this group identifier

    # Get the project view mode
    set PState [GiD_Info Project ViewMode]

    # Switch state
    if {($PState == "GEOMETRYUSE")} {
	
	set gmid "all_geometry"
	
    } elseif {($PState == "MESHUSE")} {
	
	set gmid "all_mesh"
    }
    
    # For each GiD group entities
    set EntityList [GiD_EntitiesGroups get $groupId $gmid -count]
    #msg " Groups $groupId -> $EntityList"
    if {$action != ""} {
	if {[llength $EntityList]} { 
	    return 1 
	} else { 
	    return 0 
	}
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
#
# Extraer los grupos del xml
#

proc ::KEGroups::split2 { x separator } {
    
    set splitchar \uFFFF ;#forbidden unicode character, x must never contain it
    return [split [string map "$separator $splitchar" $x] $splitchar]
}

proc ::KEGroups::GetAutomaticGroupName { {auto ""} } {
    
    set name ""
	set groups [GiD_Groups list]
    set i 0

    if { $auto == "" } {
      set name "Group"
    }
    
    if { [llength $groups] > 0 } {
	
	for {set i 1} {$i<10000} {incr i} {
	    
	    set name ${auto}${i}
	    if { [lsearch -exact $groups $name] == -1 } { break }
	}
    } else {
	if { $auto == "" } {
	    set name "Group0"
	} else {
	    set name "${auto}0"
	}
    }
    return $name
}

proc ::KEGroups::isValidGroupName {name} {
    set checklist {0 0 end end}
  
    foreach {check1 check2} $checklist {
	if { [string range $name $check1 $check2] == "/" } {
	    return 0
	} elseif { [string range $name $check1 $check2] == "/" } {
	    return 0
	} elseif { $name == "" } {
	    return 0
	} elseif { [string range $name $check1 $check2] == "." } {
	    return 0
	    
	} elseif { [string range $name $check1 $check2] == "-" } {
	    return 0
	    
	} elseif { [string range $name $check1 $check2] == "@" } {
	    return 0
	    
	} elseif { [string range $name $check1 $check2] == "<" } {
	    return 0
	    
	} elseif { [string range $name $check1 $check2] == "%" } {
	    return 0
	}
    }
    return 1
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

proc ::KEGroups::randomColor { } {

    set r [expr {int (255 * rand())}]
    set g [expr {int (255 * rand())}]
    set b [expr {int (255 * rand())}]                
    
    return [format "\#%02x%02x%02x" $r $g $b]
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

## Transfer GROUPS to .spd
proc ::KEGroups::GroupsToXml { } {

    global KPriv
    
    set xpath "/Kratos_Data/Groups"
    set basenode [$KPriv(xmlDoc) selectNodes $xpath]
    
    foreach child [$basenode childNodes] {
	$child delete
    }
    set grw "GiD_Groups"
    ::xmlutils::setXml $xpath "modeltype" "write" "GiD" "props" 1
    
    set GiDGroups [$grw list]
    foreach group $GiDGroups {
	set lgroup [split $group //]
	set path [::KEGroups::nodePath $group]
	set color [::KEGroups::randomColor]
	set state [$grw get visible $group]
	::KEGroups::insertgroupXml $path $group $color $state 
    }
}

proc ::KEGroups::insertgroupXml { path id color state } {
    
    global KPriv
    
    set basenode [$KPriv(xmlDoc) selectNodes $path]
    if { $basenode != "" } {
	    set id [lindex [split $id //] end]
      #id to xml
	    #$basenode appendXML "<Group id=\"$id\" color=\"$color\" state=\"$state\" type=\"Generic\"/>"
      set child [$KPriv(xmlDoc) createElement "Group"]
      $child setAttribute id [list $id] color [list $color] state [list $state] type Generic
      $basenode appendChild $child
      
    } else {
    
	}
}

proc ::KEGroups::nodePath { group } {
	
	set lgroups [split $group //]
	if { [llength $lgroups] == 1 } {
		set path "/Kratos_Data/Groups"
	} else {
		set path "/Kratos_Data/Groups"
		set lgroups [lrange $lgroups 0 end-1]
		foreach gr $lgroups {
		        if { $gr != "" } {
		        set path "$path/Group\[@id='$gr'\]"
		    }
		}
	}
	
	return $path
}

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

proc ::KEGroups::BorraGrupo { name } {
    global KPriv
    
    # Cuando el usuario borra un grupo de la ventana de grupos de GiD
    set deletinglist [::KMProps::findGroups $name]
    
	
    if { $deletinglist != 0 } {
		
		if { [::KEGroups::CheckGroup $name "delete"] != "-cancel-" } {
	    
		        foreach node $deletinglist {
		                $node delete
		        }
	    ::KMProps::RefreshTree
	} else {
	    return "-cancel-"
	}
    }
}

proc ::KEGroups::CheckGroup {name action} {
    
    set txt [= "Note: This action will affect all the Project Properties"]
    set msg [= "There are properties assigned to %s, do you want to %s it?\n %s" $name $action $txt]
    
    set answer [::WinUtils::confirmBox ".gid" "$msg" "okcancel"]
    if { $answer == "ok" } {
	return 1
    } else {
	return "-cancel-"
    }
}

proc ::KEGroups::RenombraGrupo { oldname newname validation } {
    # if validation == 1 -> newname needs the validation process
    set valid 1
    if {$validation == 1} {
      #future: only call of this function ::KEGroups::RenombraGrupo is with parameter validation==0
      #code part not in use
      set valid [::KEGroups::ValidateNewName $oldname $newname]
    }
    if {$valid != 0 } {
	set editinglist [::KMProps::findGroups $oldname]
	#findGroups busca si el grupo que vamos a editar está en el arbol, y devuelve a editinglist, la lista de propiedades que estaban relacionadas con ese grupo.
	global KPriv
	#msg "editing $editinglist"
	if { $editinglist != 0 } {
	    if { [winfo exists $::KMProps::WinPath] } {
		foreach node $editinglist {
                    $node setAttribute pid $newname
                    $node setAttribute id $newname
		}
		::KMProps::RefreshTree 
	    } else {
		foreach { node } $editinglist {
		    $node setAttribute pid $newname
		    $node setAttribute id $newname
		}
	    }
	} 
    }
}

#near future: proc not in use!!
proc ::KEGroups::ValidateNewName { oldname newname } {
    # POST: Returns 0 if name is invalid, returns 1 if name is valid,
    
    set ret 1
    set errlist {}
    set newtext [::KUtils::parseTreeStr $newname]

    set txt [= "You can't use some reservate chars like"]    
    if { $newname == ""} {  
	append errlist [= "You can not choose an empty group name"].
	set ret 0
    } elseif { $newtext == -1 } {

	append errlist "$txt:\n  :   /   $   .   \\  %  (space)"
	set ret 0
    } elseif { [::KEGroups::isValidGroupName $newname] != 1} {
	append errlist "$txt:\n  :   /   $   .   \\  %  (space)"
	set ret 0
    }
    if { $ret == 0 } {
	::KEGroups::ErrorRenameWindow $oldname $newname $errlist
    }
    return $ret
}

#near future: proc not in use!!
proc ::KEGroups::ErrorRenameWindow { oldname newname errlist} {
    global KPriv
    
    set qq .errwin
    if {[winfo exists $qq] } {
	destroy $qq
    }
    set newnewname ""
    # Init the window
    set title [= "Rename Error Window"]
    InitWindow $qq $title "" ""
    focus $qq
    
    if {[info exists KPriv(newnewname)] } {
	unset KPriv(newnewname)
    }
    set KPriv(newnewname) ""
    set KPriv(grnamerror) ""
    
    grid [frame $qq.main] -sticky nswe -padx 4 -pady 4
    grid [ttk::label $qq.main.msg -text $errlist] -sticky nswe -pady 2 -padx 2 -column 0 -row 0
    grid [ttk::label $qq.main.ask -text [= "Please choose a new name"]:] -sticky nswe -pady 2 -padx 2 -column 0 -row 1
    grid [ttk::entry $qq.main.ent -textvariable KPriv(newnewname) -validate key] -sticky nswe -pady 2 -padx 2 -column 1 -row 1
    grid [ttk::button $qq.main.ok -text [= "Ok"] -command [list ::KEGroups::AuxRenameGroup $oldname $newname $qq]] -sticky nswe -pady 2 -padx 2 -ipady 2 -column 1 -row 2
 
    bind $qq <Return> [list ::KEGroups::AuxRenameGroup $oldname $newname $qq]
}


#near future: proc not in use!!
proc ::KEGroups::AuxRenameGroup { oldname newname qq} {
	global KPriv

	 set ret 1         
	  set newtext [::KUtils::parseTreeStr $KPriv(newnewname)]
	  if { $KPriv(newnewname) == ""} {              
	    set ret 0
	  } elseif { $newtext == -1 } {
	    set ret 0
	  } elseif { [::KEGroups::isValidGroupName $KPriv(newnewname)] != 1} {
	    set ret 0
	  }
    if { $ret==0 } {
      set KPriv(grnamerror) [= "Choose a valid name"]
    } else {
	 set valid 0
	 
	 GiD_Groups edit rename $newname $KPriv(newnewname)
   GidUtils::UpdateWindow GROUPS
   
	 ::KEGroups::RenombraGrupo $oldname $KPriv(newnewname) $valid
	 
	 if {[winfo exists $qq]} {
	     destroy $qq
	 }
	 if {[info exists KPriv(newnewname)] } {
	     unset KPriv(newnewname)
	 }
	 return
     }
     
     if { [winfo exists  $qq.main.warn] } {
	 destroy  $qq.main.warn
     }
     grid [ttk::label $qq.main.warn -text $KPriv(grnamerror)] -sticky nswe -pady 2 -padx 2 -column 0 -row 3
     $qq.main.warn configure -foreground red
}

proc ::KEGroups::TransferCondGroupstoGiDGroups { } {
    # Transfers from Cond_Groups -> GiD_Groups (&entities), deleting the old ones
    global KPriv
    set oldGroupList [Cond_Groups list]
    if {[info exists ::KPriv(Groups,DeleteGroup)]} {
	set oldvar $::KPriv(Groups,DeleteGroup)
	set ::KPriv(Groups,DeleteGroup) 0
    }
    
    foreach oldGr $oldGroupList {
	set GiDGroups [GiD_Groups list]
	if {$oldGr in $GiDGroups} {
	    # First, create the new groups with the old data
	    GiD_Groups delete $oldGr
	}
	#set clr [Cond_Groups get color $oldGr]
	set state [Cond_Groups get visible $oldGr]
	GiD_Groups create $oldGr
	#GiD_Groups edit color $oldGr $clr
	GiD_Groups edit visible $oldGr $state
	
	# Then we fill the new group with the old group's conditions
	foreach entity [list point line surface volume nodes element faces] {
	    set entityList [::KEGroups::getGroupGiDEntities $oldGr $entity]
	    # msg "Group $oldGr"
	    # msg "Entity $entity"
	    # msg "EntityList $entityList"
	    # msg "\n"
	   
	    if {[llength $entityList]} {
		GiD_EntitiesGroups assign $oldGr $entity $entityList
	    }
	}
	
	# Finally we can delete from the Old Group's System
	
	#msg $oldGr
	Cond_Groups delete $oldGr
	#msg $oldGr
	}
    
    if {[info exists ::KPriv(Groups,DeleteGroup)]} {
     set ::KPriv(Groups,DeleteGroup) $oldvar
    }
    GidUtils::UpdateWindow GROUPS
}

proc ::KEGroups::TransferGiDGroupstoCondGroups { { what "all" } } {
    # Transfers from GiD_Groups -> Cond_Groups (&entities)

    global KPriv
    
    set oldGroupList [GiD_Groups list]
    set condgrouplist [Cond_Groups list]
    
    if {[info exists ::KPriv(Groups,DeleteGroup)]} {
	set oldvar $::KPriv(Groups,DeleteGroup)
	set ::KPriv(Groups,DeleteGroup) 0
    }
    
    foreach oldGr $oldGroupList {
    
	if {$oldGr in $condgrouplist} {
	    Cond_Groups delete $oldGr
	}
	
		#set clr [Cond_Groups get color $oldGr]
		set state [GiD_Groups get visible $oldGr]
	Cond_Groups create $oldGr
	#GiD_Groups edit color $oldGr $clr
	Cond_Groups edit visible $oldGr $state
	
	# Then we fill the new group with the old group's conditions
	foreach entity [list elements nodes] {
		    set entityList [GiD_EntitiesGroups get $oldGr $entity]
	    if {$entityList != ""} {
		set entidad [append entity "_groups"]
##                 msg "entidad $entidad , grupo $oldGr , lista $entityList"
 #                 msg "Group $oldGr"
 #                 msg "Entity $entidad"
 #                 msg "EntityList $entityList"
 #                 msg " "
 ##
		eval [GiD_Process MEscape Data Conditions AssignCond $entidad NoRepeatField name Change $oldGr $entityList MEscape]
		#msg [::KEGroups::getGroupGiDEntities $oldGr $entity]
	    }
	}
    }
    if {[info exists ::KPriv(Groups,DeleteGroup)]} {
	set ::KPriv(Groups,DeleteGroup) $oldvar
    }
}

proc ::KEGroups::DestroyOldCondGroups { } {
    # Destruye todos los grupos antiguos
    
    global KPriv
    set oldGroupList [Cond_Groups list]
    if {[info exists ::KPriv(Groups,DeleteGroup)]} {
	set oldvar $::KPriv(Groups,DeleteGroup)
	set ::KPriv(Groups,DeleteGroup) 0
    }
    foreach oldGr $oldGroupList {
	Cond_Groups delete $oldGr
    }
    if {[info exists ::KPriv(Groups,DeleteGroup)]} {
	set ::KPriv(Groups,DeleteGroup) $oldvar
    }
}

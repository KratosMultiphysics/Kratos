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
#        1.8- 27/05/12- J. Garate, ::KEGroups::listReplace
#        1.7- 07/05/12- J. Garate, update renaming groups restrictions
#        1.6- 04/05/12-G. Socorro, update some procedures
#        1.5- 03/05/12-J. Garate, GiD Groups transfer to .spd // W Child
#        1.4- 03/05/12-J. Garate, GiD Groups transfer to .spd
#        1.3- 26/04/12-G. Socorro, change GiD_Groups by Cond_Groups
#        1.2- 26/03/12- J. Gárate, Cambio de iconos para AutoGroup
#        1.1- 26/03/12- J. Gárate, Renombrado de grupos. Ventana de error en el renombrado.
#        1.0- 22/03/12- J. Gárate, Cambio a funciones públicas de los grupos de GiD. Borrado de funciones de Grupos antiguas, ahora están en OLDKEGROUPS.TCL
#        0.9- 20/03/12- J. Garate, Renombrado de grupos y actualizacion del arbol tras el evento de renombrado
#        0.8- 19/03/12- J. Garate, Borrado de grupos y actualizacion del arbol tras el evento de borrado
#        0.7- 12/03/12- J. Garate AutoGroup arreglado para INSIDE Window
#        0.6- 12/03/12- J. Garate Adaptacion a los nuevos Grupos de GiD. Pendiente: Eliminar funciones antiguas.
#        0.5- 7/02/12- J. Garate Actualizada la funcion ::KEGroups::getGroupGiDEntities
#        0.4- 22/06/11-G. Socorro, delete snit, tdom and xmlstruct from the package require
#        0.3- 13/05/10-G. Socorro, add the procedure ::KEGroups::RenameGroupIdGiDCond to rename the group identifier in the GiD condition database
#        0.2- 19/03/10-LCA, Reparar acciones en el arbol (delete masivo, rename, ...) y aumentar el número de niveles a 5
#        0.1- 01/11/09-G. Socorro, create a base source code from the GiD layer.tcl script
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
	
		msg "Any group selected"
	}
}

proc ::KEGroups::SelectionAssign { entity GroupId WinPath } {
    
    # Get the condition identifier
    #set condname [::groupProp::EditSelectionGetConditionId $groupid "Assign"]
    
	if { [string range $entity 0 6] == "element" } {
		#WarnWin [= "Element assignation unavailable"]
		#return ""
		set condname [string range $entity 8 end]
		append condname "_groups"
	} else {
    
		if { $entity == "node" } {
		        set condname "point_groups"
		} else {
		        set condname "${entity}_groups"
		}
	}
    
	set OldSmallWinSelecting [GiD_Info variable SmallWinSelecting]
	if {$OldSmallWinSelecting == 0 } {
		set SmallWinSelecting 1
		GiD_Set SmallWinSelecting $SmallWinSelecting
	} else {
		set SmallWinSelecting 1
	}
    
	#set OldSmallWinSelecting $entity

	# set Location $::KMProps::Location
	#msg $Location
	#msg $WinPath
	# ::GidUtils::DisableWarnLine
	FinishButton $WinPath $WinPath [= "Press 'Finish' to stop the entities selection"] "::GidUtils::EnableWarnLine" disableall $SmallWinSelecting
	# Try to assign conditions
	#::KEGroups::RestoreWindow $WinPath
	#msg $condname
	GiD_Process MEscape
	GiD_Process Data Conditions AssignCond $condname NoRepeatField name
	GiD_Process change $GroupId
	#msg "GroupId:  $GroupId   .  condname:$condname \nentity:$entity WinPath:$WinPath"
    
	# ::GidUtils::EnableWarnLine

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
		    lassign $CondProp CId CEId - CGroupId
		    #WarnWinText "entity : $entity // Cprop: $CondProp // gmid: $gmid"
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
                if {$CId == "E" } {
                lappend EntityIdList $CEId

                } else {
                    if {$CId == "N" } {
                         lappend EntityIdList $CEId
                    } else {
                            lappend EntityIdList $CId
                    }
                }
		    }
		}
    }
    
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
#
# Extraer los grupos del xml
#

proc ::KEGroups::split2 { x separator } {
    
    set splitchar \uFFFF ;#forbidden unicode character, x must never contain it
    return [split [string map "$separator $splitchar" $x] $splitchar]
}

proc ::KEGroups::GetAutomaticGroupName { {auto ""} } {
    
    set name ""
    set groups [Cond_Groups list]
    set i 0
    #foreach grup $KPriv(groupsId) {
    #msg "group$i:KPriv $grup"
    #incr $i
    #}
    if { [llength $groups] > 0 } {
	
	for {set i 1} {$i<10000} {incr i} {
	    
	    set name ${auto}Group${i}
	    if { [lsearch -exact $groups $name] == -1 } { break }
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
	set GiDGroups [Cond_Groups list]
	foreach group $GiDGroups {
		set lgroup [split $group //]
		set path [::KEGroups::nodePath $group]
		set color [Cond_Groups get color $group]
		set state [Cond_Groups get visible $group]
		::KEGroups::insertgroupXml $path $group $color $state 
	}
}

proc ::KEGroups::insertgroupXml { path id color state } {
    
    global KPriv
    
    #$KPriv(xmlDoc) lappend "$xpath/Group id=\"$id\" color=\"$color\" state=\"$state\" " ""
    set basenode [$KPriv(xmlDoc) selectNodes $path]
    if { $basenode != "" } {
	    set id [lindex [split $id //] end]
	    $basenode appendXML "<Group id=\"$id\" color=\"$color\" state=\"$state\" type=\"Generic\"/>"
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
    #msg "deleting $deletinglist"
    if { $deletinglist != 0 } {
	#wa "si-2"
	if { [::KEGroups::CheckGroup $name "delete"] != "-cancel-" } {
	    if { [winfo exists $::KMProps::WinPath] } {
		foreach { fullname T item } $deletinglist {
		    
		    #Elimina el grupo del xml
		    ::xmlutils::unsetXml $fullname
		    #Elimina el grupo del árbol
		    ::KMProps::deleteItem $T $item
		    
		}
		::KMProps::RefreshTree $T
	    } else {
		foreach { node } $deletinglist {
		    $node delete
		    
		}
	    }
	} else {
	    return "-cancel-"
	}
    }
}

proc ::KEGroups::CheckGroup {name action} {
	
    set aviso "There are properties assigned to $name , Do you want to $action it? \n Note: This action will affect all the Project Properties"
    
    set answer [::WinUtils::confirmBox "." "$aviso" "okcancel"]
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
		set valid [::KEGroups::ValidateNewName $oldname $newname]
	}
	if {$valid != 0 } {
		set editinglist [::KMProps::findGroups $oldname]
		#findGroups busca si el grupo que vamos a editar está en el arbol, y devuelve a editinglist, la lista de propiedades que estaban relacionadas con ese grupo.
		global KPriv
		#msg "editing $editinglist"
		if { $editinglist != 0 } {
		        if { [winfo exists $::KMProps::WinPath] } {
		            foreach { fullname T item } $editinglist {

		                ::xmlutils::setXml $fullname pid "write" $newname
		                ::xmlutils::setXml $fullname id "write" $newname
		            }
		                ::KMProps::RefreshTree $T
		        } else {
		            foreach { node } $editinglist {
		                $node setAttribute pid $newname
		                $node setAttribute id $newname
		            }
		        }
		} 
	}
}

proc ::KEGroups::ValidateNewName { oldname newname } {
    # POST: Returns 0 if name is invalid, returns 1 if name is valid,
    
    set ret 1
    set errlist {}
    set newtext [::KUtils::parseTreeStr $newname]
    
    if { $newname == ""} {  
	append errlist "You can not choose an empty group name."
	set ret 0
    } elseif { $newtext == -1 } {
	append errlist "You can't use some reservate chars like:\n  :   /   $   .   \\  %  "
	set ret 0
    } elseif { [::KEGroups::isValidGroupName $newname] != 1} {
	append errlist "You can't use some reservate chars like:\n  :   /   $   .   \\  %  "
	set ret 0
    }
    if { $ret == 0 } {
	::KEGroups::ErrorRenameWindow $oldname $newname $errlist
    }
    return $ret
}

proc ::KEGroups::ErrorRenameWindow { oldname newname errlist} {
	set qq .errwin
	global KPriv
	if { [winfo exists $qq] } {
		catch {destroy $qq}
	}
	set newnewname ""
	# Init the window
	set title [= "Rename Error Window"]
	InitWindow $qq $title "" ""
	focus $qq
	#toplevel $qq
	if {[info exists Kpriv(newnewname)] } {
		unset Kpriv(newnewname)
	}
	set Kpriv(newnewname) ""
	set Kpriv(grnamerror) ""
	
	grid [frame $qq.main] -sticky nswe -padx 4 -pady 4
	grid [ttk::label $qq.main.msg -text $errlist] -sticky nswe -pady 2 -padx 2 -column 0 -row 0
	grid [ttk::label $qq.main.ask -text "Please choose a new name: "] -sticky nswe -pady 2 -padx 2 -column 0 -row 1
	grid [ttk::entry $qq.main.ent -textvariable KPriv(newnewname) -validate key] -sticky nswe -pady 2 -padx 2 -column 1 -row 1
	grid [ttk::button $qq.main.ok -text "Ok" -command [list ::KEGroups::AuxRenameGroup $oldname $newname $qq] ] -sticky nswe -pady 2 -padx 2 -ipady 2 -column 1 -row 2
	bind $qq <Return> [list ::KEGroups::AuxRenameGroup $oldname $newname $qq]        

}

 proc ::KEGroups::AuxRenameGroup { oldname newname qq} {
	global KPriv
	if { $KPriv(newnewname) == ""} {  
		set Kpriv(grnamerror) "Choose a valid name"
		
	} elseif { [::KEGroups::isValidGroupName $KPriv(newnewname)] != 1} {
		set Kpriv(grnamerror) "Choose a valid name"
		
	} elseif {$KPriv(newnewname) in [Cond_Groups list]} {
		set Kpriv(grnamerror) "Group name already exists"
		
	} else {
		set valid 0
		
		Cond_Groups edit rename $newname $KPriv(newnewname)
		::KEGroups::RenombraGrupo $oldname $KPriv(newnewname) $valid
		
		catch {destroy $qq}
		catch {destroy $KPriv(newnewname)}
		return
	}
	if { [winfo exists  $qq.main.warn] } {
		catch {destroy  $qq.main.warn}
	}
		grid [ttk::label $qq.main.warn -text $Kpriv(grnamerror)] -sticky nswe -pady 2 -padx 2 -column 0 -row 3
		$qq.main.warn configure -foreground red
}
namespace eval Pfem::write {
}

proc Pfem::write::Init { } {
}

proc Pfem::write::writeParametersEvent { } {
    write::WriteJSON [getParametersDict]
}

# Model Part Blocks
proc Pfem::write::writeModelPartEvent { } {
    set parts_un_list [GetPartsUN]
    foreach part_un $parts_un_list {
        write::initWriteData $part_un "PFEM_Materials"
    }
    
    write::writeModelPartData
    write::WriteString "Begin Properties 0"
    write::WriteString "End Properties"
    write::writeMaterials "Pfem"
    
    write::writeNodalCoordinates
    foreach part_un $parts_un_list {
        write::initWriteData $part_un "PFEM_Materials"
        write::writeElementConnectivities
    }
    Solid::write::writeConditions
    Pfem::write::writeMeshes
}

proc Pfem::write::writeMeshes { } {
    
    foreach part_un [GetPartsUN] {
        write::initWriteData $part_un "PFEM_Materials"
        write::writePartMeshes
    }
    # Solo Malla , no en conditions
    write::writeNodalConditions "PFEM_NodalConditions"
    
    # A Condition y a meshes-> salvo lo que no tenga topologia
    Solid::write::writeLoads
}

proc Pfem::write::GetPartsUN { } {
    customlib::UpdateDocument
    set lista [list ]
    set root [customlib::GetBaseRoot]
    set xp1 "[spdAux::getRoute "PFEM_Bodies"]/blockdata/condition"
    set i 0
    foreach part_node [$root selectNodes $xp1] {
        if {![$part_node hasAttribute "un"]} {
            set un "PFEM_Part$i"
            while {[spdAux::getRoute $un] ne ""} {
                incr i
                set un "PFEM_Part$i"
            }
            $part_node setAttribute un $un
            spdAux::setRoute $un [$part_node toXPath]
        }
        lappend lista [get_domnode_attribute $part_node un]
    }
    customlib::UpdateDocument
    return $lista
}
# Custom files (Copy python scripts, write materials file...)
proc Pfem::write::writeCustomFilesEvent { } {

}


Pfem::write::Init

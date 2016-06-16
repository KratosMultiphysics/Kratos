namespace eval Pfem::write {
}

proc Pfem::write::Init { } {
    # Namespace variables inicialization
}

# Project Parameters
proc Pfem::write::getParametersDict { } {
    set projectParametersDict [dict create]
    
    # Problem data
    # Create section
    set problemDataDict [dict create]
    
    # Add items to section
    set model_name [file tail [GiD_Info Project ModelName]]
    dict set problemDataDict problem_name $model_name
    
    set nDim [write::getValue nDim]
    dict set problemDataDict domain_size $nDim
    
    ## Parallelization
    #set paralleltype [write::getValue ParallelType]
    #if {$paralleltype eq "OpenMP"} {
    #    #set nthreads [write::getValue Parallelization OpenMPNumberOfThreads]
    #    #dict set problemDataDict NumberofThreads $nthreads
    #} else {
    #    #set nthreads [write::getValue Parallelization MPINumberOfProcessors]
    #    #dict set problemDataDict NumberofProcessors $nthreads
    #}
    
    # Add section to document
    dict set projectParametersDict problem_data $problemDataDict
    
    # Nodal data processes
    set nodalDataDict [GetNodalDataDict]
    dict set projectParametersDict nodal_data_process_list $nodalDataDict
    
    return $projectParametersDict
}


proc Pfem::write::GetNodalDataDict { } {
    set doc $gid_groups_conds::doc
    set root [$doc documentElement]
    set NodalData [list ]
    set parts [list "PFEM_Rigid2DParts" "PFEM_Rigid3DParts" "PFEM_Deformable2DParts" "PFEM_Deformable3DParts" "PFEM_Fluid2DParts" "PFEM_Fluid3DParts"]
    
    foreach part $parts {
        set xp1 "[spdAux::getRoute $part]/group"
        set groups [$root selectNodes $xp1]
        foreach group $groups {
            set partid [[$group parent] @n]
            set groupid [$group @n]
            set processDict [dict create]
            dict set processDict process_name "ApplyValuesToNodes"
            dict set processDict implemented_in_module "KratosMultiphysics.PFEMBaseApplication"
            
            set params [dict create]
            set xp2 "./value"
            set atts [$group selectNodes $xp2]
            #W "$group $groupid $atts"
            foreach att $atts {
                set state [get_domnode_attribute $att state]
                if {$state ne "hidden"} {
                    set paramName [$att @n]
                    set paramValue [get_domnode_attribute $att v]
                    if {[write::isBoolean $paramValue]} {set paramValue [expr $paramValue]}
                    dict set params $paramName $paramValue
                }
            }
            dict set params "model_part_name" [::write::getMeshId $partid $groupid]
            dict set processDict "Parameters" $params
            lappend NodalData $processDict
        }
    }
    
    return $NodalData
}
proc Pfem::write::writeParametersEvent { } {
    write::WriteJSON [getParametersDict]
}

# Model Part Blocks
proc Pfem::write::writeModelPartEvent { } {
    write::WriteString "Begin ModelPartData"
    write::WriteString "End ModelPartData"
    write::WriteString "Begin Properties 0"
    write::WriteString "End Properties"
    
    
    write::writeNodalCoordinates
    
     # Solo Solidos
    write::writeNodalConditions "PFEM_Rigid2DParts"
    write::writeNodalConditions "PFEM_Rigid3DParts"
    
    # Solo Fluidos
    write::writeNodalConditions "PFEM_Fluid2DParts"
    write::writeNodalConditions "PFEM_Fluid3DParts"
    
    writeMeshingDomains
}

proc Pfem::write::writeMeshingDomains { } {
    set doc $gid_groups_conds::doc
    set root [$doc documentElement]
    
    set meshing_domains(dummy) [list ]
    set parts [list "PFEM_Rigid2DParts" "PFEM_Rigid3DParts" "PFEM_Deformable2DParts" "PFEM_Deformable3DParts" "PFEM_Fluid2DParts" "PFEM_Fluid3DParts"]
    
    foreach part $parts {
        set xp1 "[spdAux::getRoute $part]/group"
        set groups [$root selectNodes $xp1]
        foreach group $groups {
            set groupid [$group @n]
            set md [get_domnode_attribute [$group selectNodes ".//value\[@n='MeshingDomain']"] v]
            lappend meshing_domains($md) $groupid
        }
    }
    array unset meshing_domains dummy
    #Kratos::PrintArray meshing_domains
    
    foreach domain [array names meshing_domains] {
        set gname "_TEMP$domain"
        GiD_Groups create $gname
        foreach group $meshing_domains($domain) {
            GiD_EntitiesGroups assign $gname nodes [GiD_EntitiesGroups get $group nodes]
        }
        set domname [string map {" " "_"} "$domain"]
        write::writeGroupMesh "$domname" $gname nodes
        GiD_Groups delete $gname
    }
}

# Custom files (Copy python scripts, write materials file...)
proc Pfem::write::writeCustomFilesEvent { } {

}


Pfem::write::Init

proc ExtractSurfaceTriangles_ext { basename dir problemtypedir } {

    ## Source auxiliar procedures
    source [file join $problemtypedir OBJFileAuxProcs.tcl]

    ## Start OBJ file
    #set filename [file join $dir ${basename}.obj]
    set filename [file join $dir generic.obj]
    set FileVar [open $filename w]

    # set all_mesh [GiD_EntitiesLayers get Layer0 all_mesh]
    # If <over> is all_mesh then is obtained a list with 2 sublists: node id's and element id's
 
    #set Entities [GiD_EntitiesGroups get [lindex $Group 1] elements -element_type $ElemType]
    set triangles [GiD_Mesh list -element_type {triangle} element]
    set tetrahedra [GiD_Mesh list -element_type {tetrahedra} element]
    
    W "list of triangles"
    W $triangles

    W "list of tetrahedra"
    W $tetrahedra
    
    # puts $FileVar "Face Elements iD - not required in obj file"
    # foreach element_id $triangles 
    #   puts $FileVar " f $element_id"
    # 
    # puts $FileVar "End Face Elements"
    # puts $FileVar ""

    
    
    # set element_ids [GiD_EntitiesGroups get $groupid elements] ;               # get ids of all elements in cgroupid
    # #array set is_external_element [DEM::write::Compute_External_Elements 3 $groupid $element_ids]

    # aux sumar a cada node de cada triangle la normal del triangle
    # aux nombre de normals de cada node

    set triangle_nodes [list]
    foreach element_id $triangles { ;
        set nodes [lrange [GiD_Mesh get element $element_id] 3 end] ;   # get the nodes of the element
        lappend triangle_nodes {*}$nodes ;                              # add those nodes to the nodes_to_delete list
    }


    ## Nodes only if belong to any triangle
    set Nodes [GiD_Info Mesh Nodes]
    #puts $FileVar "Ordered vertex list only if associated with any face triangle"
    for {set i 0} {$i < [llength $Nodes]} {incr i 4} {
        if {[lindex $Nodes $i] in $triangle_nodes} {
        puts -nonewline $FileVar "v "
        puts -nonewline $FileVar [format  "%.10f" [lindex $Nodes [expr { $i+1 }]]]
        puts -nonewline $FileVar " "
        puts -nonewline $FileVar [format  "%.10f" [lindex $Nodes [expr { $i+2 }]]]
        puts -nonewline $FileVar " "
        puts $FileVar [format  "%.10f" [lindex $Nodes [expr { $i+3 }]]]
        }
    }

    # Lo que tienes que hacer es recorrer los elementos que son los que apuntan a
    # sus nodos, y sumar a cada nodo la aportación de su normal.
    # Aun así, cualquier función script que haga cosas como esta en malla puede
    # ser lenta en mallas muy grandes, por crear arrays o listas grandes, pero si
    # la haces mal no será lenta sino muy-muy-lenta
    
    set aux 1
    dict set node_db_x $aux 0.0
    dict unset node_db_x $aux

    dict set node_db_y $aux 0.0
    dict unset node_db_y $aux

    dict set node_db_z $aux 0.0
    dict unset node_db_z $aux

    dict set node_db_counter $aux 0
    dict unset node_db_counter $aux

    foreach element_id $triangles { ;
        W "___________next element of the list"
        set nodes [lrange [GiD_Mesh get element $element_id] 3 end] ;
        W $nodes
        set ElementNormal_x [lrange [GiD_Mesh get element $element_id normal] 0 0] ; 
        set ElementNormal_y [lrange [GiD_Mesh get element $element_id normal] 1 1] ; 
        set ElementNormal_z [lrange [GiD_Mesh get element $element_id normal] 2 2] ; 
        
        foreach node $nodes {
            if {[dict exists $node_db_x $node]} {
                set old_value [dict get $node_db_x $node]
                set new_value [expr { $old_value+$ElementNormal_x}]
                # set new_value [expr { $old_value[0]+$ElementNormal_x }] [expr { $old_value[1]+$ElementNormal_y }]
                dict unset node_db_x $node
                dict set node_db_x $node $new_value

                dict incr node_db_counter $node 1
            } else {
                dict set node_db_x $node $ElementNormal_x
                dict set node_db_counter $node 1
            }
        }
        W "foreach element x "
        W $node_db_x

        foreach node $nodes {
            if {[dict exists $node_db_y $node]} {
                set old_value [dict get $node_db_y $node]
                set new_value [expr { $old_value+$ElementNormal_y}]
                dict unset node_db_y $node
                dict set node_db_y $node $new_value
            } else {
                dict set node_db_y $node $ElementNormal_y
            }
        }
        W "foreach element y "
        W $node_db_y

        foreach node $nodes {
            if {[dict exists $node_db_z $node]} {
                set old_value [dict get $node_db_z $node]
                set new_value [expr { $old_value+$ElementNormal_z}]
                dict unset node_db_z $node
                dict set node_db_z $node $new_value
            } else {
                dict set node_db_z $node $ElementNormal_z
            }
        }
        W "foreach element z"
        W $node_db_z
        W $node_db_counter

    }


    # Get a list of the keys and sort them
    W "BEFORE ORDERING KEYS"
    W $node_db_x
    set node_db_x [lsort -stride 2 $node_db_x]
    set node_db_y [lsort -stride 2 $node_db_y]
    set node_db_z [lsort -stride 2 $node_db_z]
    W "AFTER ORDERING KEYS"
    W $node_db_x        

    foreach item [dict keys $node_db_x] {
            set val [dict get $node_db_x $item]
            # set counter [dict get $node_db_counter $item]
            # set avg_val_normal [expr { $val/$counter }]
            # lappend normals_x $avg_val_normal ;
            lappend normals_x $val ;        
    }

    foreach item [dict keys $node_db_y] {
            set val [dict get $node_db_y $item]
            # set counter [dict get $node_db_counter $item]
            # set avg_val_normal [expr { $val/$counter }]
            # lappend normals_y $avg_val_normal ; 
            lappend normals_y $val ;          
    }

    foreach item [dict keys $node_db_z] {
            set val [dict get $node_db_z $item]
            # set counter [dict get $node_db_counter $item]
            # set avg_val_normal [expr { $val/$counter }]
            # lappend normals_z $avg_val_normal ; 
            lappend normals_z $val ;         
    }

    W $normals_x
    W $normals_y
    W $normals_z

    # foreach component1 $normals_x component2 $normals_y component3 $normals_z {
    #     lappend normals $component1 $component2 $component3 ;
    # }
    # # W $normals

    # normalization is required.
    foreach component1 $normals_x component2 $normals_y component3 $normals_z {
        set magnitude [expr { sqrt( $component1 * $component1 + $component2 * $component2 + $component3 * $component3) }]
        set component1_normalized [expr { $component1/$magnitude }]
        set component2_normalized [expr { $component2/$magnitude }]
        set component3_normalized [expr { $component3/$magnitude }]
        lappend normals $component1_normalized $component2_normalized $component3_normalized ;
    }
    W $normals

    

    

    ## vertex normals only if belong to any triangle
    puts $FileVar " "
    #puts $FileVar "Ordered vertex normals list only if associated with any face triangle"
    for {set i 0} {$i < [llength $normals]} {incr i 3} {
        puts -nonewline $FileVar "vn "
        puts -nonewline $FileVar [format  "%.10f" [lindex $normals [expr { $i }]]]
        puts -nonewline $FileVar " "
        puts -nonewline $FileVar [format  "%.10f" [lindex $normals [expr { $i+1 }]]]
        puts -nonewline $FileVar " "
        puts $FileVar [format  "%.10f" [lindex $normals [expr { $i+2 }]]]
    }


    # required format is
    # v 987.823009 -583.341002 122.360997    
    # vn 0.329248 0.150250 0.932213
    # v 979.430974 -499.442995 88.674001
    # vn 0.689488 0.597651 0.409169

    # f 18//18 10//10 9//9     f v1//vn1 v2//vn2 v3//vn3 ...
    # f 18//18 9//9 17//17 
    # f 19//19 11//11 10//10 


    #puts $FileVar "Triangle faces defined by its associated vertex"
    #f v1//vn1 v2//vn2 v3//vn3
    puts $FileVar " "
    foreach element_id $triangles { ;
        set nodes [lrange [GiD_Mesh get element $element_id] 3 end] ;
        puts $FileVar "f [lindex $nodes 0]//[lindex $nodes 0] [lindex $nodes 1]//[lindex $nodes 1] [lindex $nodes 2]//[lindex $nodes 2]" ;
    }






























    # ## Conditions
    # set ConditionId 0
    # set ConditionDict [dict create]
    # set Dim [GiD_AccessValue get gendata Domain_Size]
    # # Force
    # set Groups [GiD_Info conditions Force groups]
    # if {$Dim eq 2} {
    #     # UPwForceCondition2D1N
    #     WriteNodalConditions FileVar ConditionId ConditionDict $Groups UPwForceCondition2D1N $BodyElemsProp
    # } else {
    #     # UPwForceCondition3D1N
    #     WriteNodalConditions FileVar ConditionId ConditionDict $Groups UPwForceCondition3D1N $BodyElemsProp
    # }
    # # Face_Load
    # set Groups [GiD_Info conditions Face_Load groups]
    # if {$Dim eq 2} {
    #     if {$IsQuadratic eq 0} {
    #         # UPwFaceLoadCondition2D2N
    #         WriteFaceConditions FileVar ConditionId ConditionDict $Groups UPwFaceLoadCondition2D2N $PropertyDict
    #     } else {
    #         # LineLoadDiffOrderCondition2D3N
    #         WriteFaceConditions FileVar ConditionId ConditionDict $Groups LineLoadDiffOrderCondition2D3N $PropertyDict
    #     }
    # } else {
    #     if {$IsQuadratic eq 0} {
    #         for {set i 0} {$i < [llength $Groups]} {incr i} {
    #             set MyConditionList [list]
    #             # UPwFaceLoadCondition3D3N
    #             WriteTypeFaceConditions FileVar ConditionId MyConditionList [lindex $Groups $i] tetrahedra UPwFaceLoadCondition3D3N $PropertyDict
    #             WriteTypeFaceConditions FileVar ConditionId MyConditionList [lindex $Groups $i] prism UPwFaceLoadCondition3D3N $PropertyDict
    #             # UPwFaceLoadCondition3D4N
    #             WriteTypeFaceConditions FileVar ConditionId MyConditionList [lindex $Groups $i] hexahedra UPwFaceLoadCondition3D4N $PropertyDict
    #             dict set ConditionDict [lindex [lindex $Groups $i] 1] $MyConditionList
    #         }
    #     } elseif {$IsQuadratic eq 1} {
    #         for {set i 0} {$i < [llength $Groups]} {incr i} {
    #             set MyConditionList [list]
    #             # SurfaceLoadDiffOrderCondition3D6N
    #             WriteTypeFaceConditions FileVar ConditionId MyConditionList [lindex $Groups $i] tetrahedra SurfaceLoadDiffOrderCondition3D6N $PropertyDict
    #             # SurfaceLoadDiffOrderCondition3D8N
    #             WriteTypeFaceConditions FileVar ConditionId MyConditionList [lindex $Groups $i] hexahedra SurfaceLoadDiffOrderCondition3D8N $PropertyDict
    #             dict set ConditionDict [lindex [lindex $Groups $i] 1] $MyConditionList
    #         }
    #     } else {
    #         for {set i 0} {$i < [llength $Groups]} {incr i} {
    #             set MyConditionList [list]
    #             # SurfaceLoadDiffOrderCondition3D6N
    #             WriteTypeFaceConditions FileVar ConditionId MyConditionList [lindex $Groups $i] tetrahedra SurfaceLoadDiffOrderCondition3D6N $PropertyDict
    #             # SurfaceLoadDiffOrderCondition3D9N
    #             WriteTypeFaceConditions FileVar ConditionId MyConditionList [lindex $Groups $i] hexahedra SurfaceLoadDiffOrderCondition3D9N $PropertyDict
    #             dict set ConditionDict [lindex [lindex $Groups $i] 1] $MyConditionList
    #         }
    #     }
    # }
    # # Normal_Load
    # set Groups [GiD_Info conditions Normal_Load groups]
    # if {$Dim eq 2} {
    #     if {$IsQuadratic eq 0} {
    #         # UPwNormalFaceLoadCondition2D2N
    #         WriteFaceConditions FileVar ConditionId ConditionDict $Groups UPwNormalFaceLoadCondition2D2N $PropertyDict
    #     } else {
    #         # LineNormalLoadDiffOrderCondition2D3N
    #         WriteFaceConditions FileVar ConditionId ConditionDict $Groups LineNormalLoadDiffOrderCondition2D3N $PropertyDict
    #     }
    # } else {
    #     if {$IsQuadratic eq 0} {
    #         for {set i 0} {$i < [llength $Groups]} {incr i} {
    #             set MyConditionList [list]
    #             # UPwNormalFaceLoadCondition3D3N
    #             WriteTypeFaceConditions FileVar ConditionId MyConditionList [lindex $Groups $i] tetrahedra UPwNormalFaceLoadCondition3D3N $PropertyDict
    #             WriteTypeFaceConditions FileVar ConditionId MyConditionList [lindex $Groups $i] prism UPwNormalFaceLoadCondition3D3N $PropertyDict
    #             # UpwNormalFaceLoadCondition3D4N
    #             WriteTypeFaceConditions FileVar ConditionId MyConditionList [lindex $Groups $i] hexahedra UpwNormalFaceLoadCondition3D4N $PropertyDict
    #             dict set ConditionDict [lindex [lindex $Groups $i] 1] $MyConditionList
    #         }
    #     } elseif {$IsQuadratic eq 1} {
    #         for {set i 0} {$i < [llength $Groups]} {incr i} {
    #             set MyConditionList [list]
    #             # SurfaceNormalLoadDiffOrderCondition3D6N
    #             WriteTypeFaceConditions FileVar ConditionId MyConditionList [lindex $Groups $i] tetrahedra SurfaceNormalLoadDiffOrderCondition3D6N $PropertyDict
    #             # SurfaceNormalLoadDiffOrderCondition3D8N
    #             WriteTypeFaceConditions FileVar ConditionId MyConditionList [lindex $Groups $i] hexahedra SurfaceNormalLoadDiffOrderCondition3D8N $PropertyDict
    #             dict set ConditionDict [lindex [lindex $Groups $i] 1] $MyConditionList
    #         }
    #     } else {
    #         for {set i 0} {$i < [llength $Groups]} {incr i} {
    #             set MyConditionList [list]
    #             # SurfaceNormalLoadDiffOrderCondition3D6N
    #             WriteTypeFaceConditions FileVar ConditionId MyConditionList [lindex $Groups $i] tetrahedra SurfaceNormalLoadDiffOrderCondition3D6N $PropertyDict
    #             # SurfaceNormalLoadDiffOrderCondition3D9N
    #             WriteTypeFaceConditions FileVar ConditionId MyConditionList [lindex $Groups $i] hexahedra SurfaceNormalLoadDiffOrderCondition3D9N $PropertyDict
    #             dict set ConditionDict [lindex [lindex $Groups $i] 1] $MyConditionList
    #         }
    #     }
    # }
    # # Normal_Fluid_Flux
    # set Groups [GiD_Info conditions Normal_Fluid_Flux groups]
    # if {$Dim eq 2} {
    #     if {$IsQuadratic eq 0} {
    #         if {$FIC eq false} {
    #             # UPwNormalFluxCondition2D2N
    #             WriteFaceConditions FileVar ConditionId ConditionDict $Groups UPwNormalFluxCondition2D2N $PropertyDict
    #         } else {
    #             # UPwNormalFluxFICCondition2D2N
    #             WriteFaceConditions FileVar ConditionId ConditionDict $Groups UPwNormalFluxFICCondition2D2N $PropertyDict
    #         }
    #     } else {
    #         # LineNormalFluidFluxDiffOrderCondition2D3N
    #         WriteFaceConditions FileVar ConditionId ConditionDict $Groups LineNormalFluidFluxDiffOrderCondition2D3N $PropertyDict
    #     }
    # } else {
    #     if {$IsQuadratic eq 0} {
    #         if {$FIC eq false} {
    #             for {set i 0} {$i < [llength $Groups]} {incr i} {
    #                 set MyConditionList [list]
    #                 # UPwNormalFluxCondition3D3N
    #                 WriteTypeFaceConditions FileVar ConditionId MyConditionList [lindex $Groups $i] tetrahedra UPwNormalFluxCondition3D3N $PropertyDict
    #                 WriteTypeFaceConditions FileVar ConditionId MyConditionList [lindex $Groups $i] prism UPwNormalFluxCondition3D3N $PropertyDict
    #                 # UPwNormalFluxCondition3D4N
    #                 WriteTypeFaceConditions FileVar ConditionId MyConditionList [lindex $Groups $i] hexahedra UPwNormalFluxCondition3D4N $PropertyDict
    #                 dict set ConditionDict [lindex [lindex $Groups $i] 1] $MyConditionList
    #             }
    #         } else {
    #             for {set i 0} {$i < [llength $Groups]} {incr i} {
    #                 set MyConditionList [list]
    #                 # UPwNormalFluxFICCondition3D3N
    #                 WriteTypeFaceConditions FileVar ConditionId MyConditionList [lindex $Groups $i] tetrahedra UPwNormalFluxFICCondition3D3N $PropertyDict
    #                 WriteTypeFaceConditions FileVar ConditionId MyConditionList [lindex $Groups $i] prism UPwNormalFluxFICCondition3D3N $PropertyDict
    #                 # UPwNormalFluxFICCondition3D4N
    #                 WriteTypeFaceConditions FileVar ConditionId MyConditionList [lindex $Groups $i] hexahedra UPwNormalFluxFICCondition3D4N $PropertyDict
    #                 dict set ConditionDict [lindex [lindex $Groups $i] 1] $MyConditionList
    #             }
    #         }
    #     } elseif {$IsQuadratic eq 1} {
    #         for {set i 0} {$i < [llength $Groups]} {incr i} {
    #             set MyConditionList [list]
    #             # SurfaceNormalFluidFluxDiffOrderCondition3D6N
    #             WriteTypeFaceConditions FileVar ConditionId MyConditionList [lindex $Groups $i] tetrahedra SurfaceNormalFluidFluxDiffOrderCondition3D6N $PropertyDict
    #             # SurfaceNormalFluidFluxDiffOrderCondition3D8N
    #             WriteTypeFaceConditions FileVar ConditionId MyConditionList [lindex $Groups $i] hexahedra SurfaceNormalFluidFluxDiffOrderCondition3D8N $PropertyDict
    #             dict set ConditionDict [lindex [lindex $Groups $i] 1] $MyConditionList
    #         }
    #     } else {
    #         for {set i 0} {$i < [llength $Groups]} {incr i} {
    #             set MyConditionList [list]
    #             # SurfaceNormalFluidFluxDiffOrderCondition3D6N
    #             WriteTypeFaceConditions FileVar ConditionId MyConditionList [lindex $Groups $i] tetrahedra SurfaceNormalFluidFluxDiffOrderCondition3D6N $PropertyDict
    #             # SurfaceNormalFluidFluxDiffOrderCondition3D9N
    #             WriteTypeFaceConditions FileVar ConditionId MyConditionList [lindex $Groups $i] hexahedra SurfaceNormalFluidFluxDiffOrderCondition3D9N $PropertyDict
    #             dict set ConditionDict [lindex [lindex $Groups $i] 1] $MyConditionList
    #         }
    #     }
    # }
    # # Interface_Face_Load
    # set Groups [GiD_Info conditions Interface_Face_Load groups]
    # for {set i 0} {$i < [llength $Groups]} {incr i} {
    #     set MyConditionList [list]
    #     # UPwFaceLoadInterfaceCondition2D2N
    #     WriteInterfaceConditions FileVar ConditionId MyConditionList [lindex $Groups $i] linear UPwFaceLoadInterfaceCondition2D2N $InterfaceElemsProp Line2D2Connectivities
    #     # UPwFaceLoadInterfaceCondition3D4N
    #     WriteInterfaceConditions FileVar ConditionId MyConditionList [lindex $Groups $i] triangle UPwFaceLoadInterfaceCondition3D4N $InterfaceElemsProp TriangleInterface3D4Connectivities
    #     WriteInterfaceConditions FileVar ConditionId MyConditionList [lindex $Groups $i] quadrilateral UPwFaceLoadInterfaceCondition3D4N $InterfaceElemsProp QuadrilateralInterface3D4Connectivities
    #     dict set ConditionDict [lindex [lindex $Groups $i] 1] $MyConditionList
    # }
    # # Interface_Normal_Fluid_Flux
    # set Groups [GiD_Info conditions Interface_Normal_Fluid_Flux groups]
    # for {set i 0} {$i < [llength $Groups]} {incr i} {
    #     set MyConditionList [list]
    #     # UPwNormalFluxInterfaceCondition2D2N
    #     WriteInterfaceConditions FileVar ConditionId MyConditionList [lindex $Groups $i] linear UPwNormalFluxInterfaceCondition2D2N $InterfaceElemsProp Line2D2Connectivities
    #     # UPwNormalFluxInterfaceCondition3D4N
    #     WriteInterfaceConditions FileVar ConditionId MyConditionList [lindex $Groups $i] triangle UPwNormalFluxInterfaceCondition3D4N $InterfaceElemsProp TriangleInterface3D4Connectivities
    #     WriteInterfaceConditions FileVar ConditionId MyConditionList [lindex $Groups $i] quadrilateral UPwNormalFluxInterfaceCondition3D4N $InterfaceElemsProp QuadrilateralInterface3D4Connectivities
    #     dict set ConditionDict [lindex [lindex $Groups $i] 1] $MyConditionList
    # }

    # # Periodic_Bars
    # set IsPeriodic [GiD_AccessValue get gendata Periodic_Interface_Conditions]
    # if {$IsPeriodic eq true} {
    #     set PeriodicBarsDict [dict create]
    #     set Groups [GiD_Info conditions Interface_Part groups]
    #     for {set i 0} {$i < [llength $Groups]} {incr i} {
    #         if {[lindex [lindex $Groups $i] 20] eq true} {
    #             # Elements Property
    #             set InterfaceElemsProp [dict get $PropertyDict [lindex [lindex $Groups $i] 1]]
    #             set ConditionList [list]
    #             # InterfaceElement2D4N
    #             SavePeriodicBarsFromIE2D4N PeriodicBarsDict ConditionId ConditionList [lindex $Groups $i] $InterfaceElemsProp
    #             # InterfaceElement3D6N
    #             SavePeriodicBarsFromIE3D6N PeriodicBarsDict ConditionId ConditionList [lindex $Groups $i] $InterfaceElemsProp
    #             # InterfaceElement3D8N
    #             SavePeriodicBarsFromIE3D8N PeriodicBarsDict ConditionId ConditionList [lindex $Groups $i] $InterfaceElemsProp

    #             dict set ConditionDict Periodic_Bars_[lindex [lindex $Groups $i] 1] $ConditionList
    #         }
    #     }

    #     if {[dict size $PeriodicBarsDict] > 0} {
    #         puts $FileVar "Begin Conditions PeriodicCondition"
    #         dict for {Name PeriodicBar} $PeriodicBarsDict {
    #             puts $FileVar "  [dict get $PeriodicBar Id]  [dict get $PeriodicBar PropertyId]  [dict get $PeriodicBar Connectivities]"
    #         }
    #         puts $FileVar "End Conditions"
    #         puts $FileVar ""
    #     }
    # }

    # puts $FileVar ""

    # ## SubModelParts
    # # Body_Part
    # WriteElementSubmodelPart FileVar Body_Part
    # # Interface_Part
    # WriteElementSubmodelPart FileVar Interface_Part
    # # PropagationUnion (InterfaceElement)
    # if {[GiD_Groups exists PropagationUnion_3d_6] eq 1} {
    #     WritePropUnionElementSubmodelPart FileVar $PropUnionElementList
    # }
    # # Solid_Displacement
    # WriteConstraintSubmodelPart FileVar Solid_Displacement $TableDict
    # # Fluid_Pressure
    # WriteConstraintSubmodelPart FileVar Fluid_Pressure $TableDict
    # # Force
    # WriteLoadSubmodelPart FileVar Force $TableDict $ConditionDict
    # # Face_Load
    # WriteLoadSubmodelPart FileVar Face_Load $TableDict $ConditionDict
    # # Normal_Load
    # WriteLoadSubmodelPart FileVar Normal_Load $TableDict $ConditionDict
    # # Normal_Fluid_Flux
    # WriteLoadSubmodelPart FileVar Normal_Fluid_Flux $TableDict $ConditionDict
    # # Interface_Face_Load
    # WriteLoadSubmodelPart FileVar Interface_Face_Load $TableDict $ConditionDict
    # # Interface_Normal_Fluid_Flux
    # WriteLoadSubmodelPart FileVar Interface_Normal_Fluid_Flux $TableDict $ConditionDict
    # # Body_Acceleration
    # WriteConstraintSubmodelPart FileVar Body_Acceleration $TableDict

    # # Periodic_Bars
    # if {$IsPeriodic eq true} {
    #     WritePeriodicBarsSubmodelPart FileVar Interface_Part $ConditionDict
    # }

    close $FileVar

    set OBJOutput [list $triangles $tetrahedra]

    return $OBJOutput
}
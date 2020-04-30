proc ExtractSurfaceTriangles_ext { basename dir problemtypedir } {

    ## Source auxiliar procedures
    source [file join $problemtypedir OBJFileAuxProcs.tcl]

    ## Start OBJ file
    set filename [file join $dir ${basename}.obj]
    set FileVar [open $filename w]

    puts $FileVar ""

    # set all_mesh [GiD_EntitiesLayers get Layer0 all_mesh]
    # W "all_mesh"
    # W $all_mesh
    # If <over> is all_mesh then is obtained a list with 2 sublists: node id's and element id's
    # List of triangles defined by: its vertex and vertex normals
    #                               faces: ordered vertex ids


    #set Entities [GiD_EntitiesGroups get [lindex $Group 1] elements -element_type $ElemType]
    set triangles [GiD_Mesh list -element_type {triangle} element]
    set tetrahedra [GiD_Mesh list -element_type {tetrahedra} element]
    W $triangles
    W $tetrahedra
    
    puts $FileVar "Begin Elements"
    # for {set j 0} {$j < [llength $triangles]} {incr j} {
    #     puts $FileVar "  [lindex $triangles $j]   [[lindex $triangles $j]]"
    # }

    foreach element_id $triangles { ;
      puts $FileVar " f $element_id"
    }
    puts $FileVar "End Elements"
    puts $FileVar ""

    

    set triangle_nodes [list]
    # set element_ids [GiD_EntitiesGroups get $groupid elements] ;               # get ids of all elements in cgroupid
    # #array set is_external_element [DEM::write::Compute_External_Elements 3 $groupid $element_ids]

    foreach element_id $triangles { ;
        set element_nodes [lrange [GiD_Mesh get element $element_id] 3 end] ;   # get the nodes of the element
        lappend triangle_nodes {*}$element_nodes ;                              # add those nodes to the nodes_to_delete list
    }


    W "triangle_nodes"
    W $triangle_nodes
    #     triangle_nodes
    # 2 6 7 6 9 7 2 6 1 6 4 1 6 9 4 9 8 4 9 7 8 7 3 8 7 2 3 2 1 3 1 4 3 4 8 3



    # set ElementInfo [GiD_Mesh get element $ElemId]
    # #ElementInfo: <layer> <elemtype> <NumNodes> <N1> <N2> ...
    # return "[lindex $ElementInfo 3] [lindex $ElementInfo 4] [lindex $ElementInfo 5]"


    # TODO: Call GenerateOBJFile with the information extracted from the surface mesh
    # GenerateOBJFile


    # required format is
    # v 987.823009 -583.341002 122.360997
    # vn 0.329248 0.150250 0.932213
    # v 979.430974 -499.442995 88.674001
    # vn 0.689488 0.597651 0.409169

    # f 18//18 10//10 9//9 
    # f 18//18 9//9 17//17 
    # f 19//19 11//11 10//10 






    ## Nodes
    set Nodes [GiD_Info Mesh Nodes]
    puts $FileVar "Begin Nodes"
    for {set i 0} {$i < [llength $Nodes]} {incr i 4} {
        puts -nonewline $FileVar "  [lindex $Nodes $i]  "
        puts -nonewline $FileVar [format  "%.10f" [lindex $Nodes [expr { $i+1 }]]]
        puts -nonewline $FileVar " "
        puts -nonewline $FileVar [format  "%.10f" [lindex $Nodes [expr { $i+2 }]]]
        puts -nonewline $FileVar " "
        puts $FileVar [format  "%.10f" [lindex $Nodes [expr { $i+3 }]]]
    }
    puts $FileVar "End Nodes"
    puts $FileVar ""
    puts $FileVar ""


 





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
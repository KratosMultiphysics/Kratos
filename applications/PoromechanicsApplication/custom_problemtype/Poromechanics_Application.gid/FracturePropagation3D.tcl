proc WriteInitialFracturesData { dir problemtypedir gidpath } {
    
    ## Source auxiliar procedures
    source [file join $problemtypedir FracturePropagation3DAuxProcs.tcl]

    ## Set BodyVolumesDict
    set BodyGroups [GiD_Info conditions Body_Part groups]
    set BodyVolumesDict [dict create]
    for {set i 0} {$i < [llength $BodyGroups]} {incr i} {
        set BodyEntities [GiD_EntitiesGroups get [lindex [lindex $BodyGroups $i] 1] volumes]
        for {set j 0} {$j < [llength $BodyEntities]} {incr j} {
            set BodyVolume [GiD_Geometry get volume [lindex $BodyEntities $j]]
            set Groups [GiD_EntitiesGroups entity_groups volumes [lindex $BodyEntities $j]]
            dict set BodyVolumesDict [lindex $BodyEntities $j] Groups $Groups
            set Surfaces [list]
            for {set k 0} {$k < [lindex $BodyVolume 1]} {incr k} {
                lappend Surfaces [lindex [lindex $BodyVolume [expr {2+$k}]] 0]
            }
            dict set BodyVolumesDict [lindex $BodyEntities $j] Surfaces $Surfaces
            if {[lindex [GiD_Info list_entities Volumes [lindex $BodyEntities $j]] 11] eq "Meshing"} {
                set MeshSize [lindex [GiD_Info list_entities Volumes [lindex $BodyEntities $j]] 17]
                set MeshSize [string trimleft $MeshSize "size="]
                dict set BodyVolumesDict [lindex $BodyEntities $j] MeshSize $MeshSize
            } else {
                dict set BodyVolumesDict [lindex $BodyEntities $j] MeshSize 0
            }
        }
    }
    
    ## Set FracturesDict
    set InterfaceGroups [GiD_Info conditions Interface_Part groups]
    if {[llength $InterfaceGroups] < 1} {
        WarnWin "ERROR: Fracture Propagation needs at least 1 pre-defined fracture tip"
    }
    # Auxiliar Dictionary of points to determine the TipPoint of each crack
    set PointsDict [dict create]
    for {set i 0} {$i < [llength $InterfaceGroups]} {incr i} {
        # Notes: eq -> for string comparison (and int), == -> for numerical comparison
        #        ne -> for string comparison (and int), != -> for numerical comparison
        if {[lindex [lindex $InterfaceGroups $i] 3] eq false} {
            set InterfaceEntities [GiD_EntitiesGroups get [lindex [lindex $InterfaceGroups $i] 1] volumes]
            for {set j 0} {$j < [llength $InterfaceEntities]} {incr j} {
                # Note: The Ids of the two contact volumes of a crack tip must be, respectively, an even and an odd number.
                #       This is always the case if one creates first one volume of a crack tip and then the other volume of the same crack tip.
                if {$j%2} {
                    # Odd number
                    set FractureHalf "Right"
                } else {
                    # Even number
                    set FractureHalf "Left"
                }
                set InterfaceVolume [GiD_Geometry get volume [lindex $InterfaceEntities $j]]
                set Surface1 [GiD_Geometry get surface [lindex [lindex $InterfaceVolume 2] 0]]
                set Surface2 [GiD_Geometry get surface [lindex [lindex $InterfaceVolume 3] 0]]
                for {set k 0} {$k < [lindex $Surface1 2]} {incr k} {
                    set Line1Id [lindex [lindex $Surface1 [expr { 9+$k }]] 0]
                    for {set l 0} {$l < [lindex $Surface2 2]} {incr l} {
                        set Line2Id [lindex [lindex $Surface2 [expr { 9+$l }]] 0]
                        if {$Line1Id eq $Line2Id} {
                            set AxisLine [GiD_Geometry get line $Line1Id]
                            dict set PointsDict [lindex $AxisLine 2] ${FractureHalf}InterfaceVolume [lindex $InterfaceEntities $j]
                            dict set PointsDict [lindex $AxisLine 2] ${FractureHalf}Line $Line1Id
                            dict set PointsDict [lindex $AxisLine 3] ${FractureHalf}InterfaceVolume [lindex $InterfaceEntities $j]
                            dict set PointsDict [lindex $AxisLine 3] ${FractureHalf}Line $Line1Id
                        }
                    }
                }
            }
        }
    }
    set FracturesDict [dict create]
    set FractureId -1
    dict for {PointId Info} $PointsDict {
        if {[dict exists $Info RightInterfaceVolume] && [dict exists $Info LeftInterfaceVolume]} {
            # Define new fracture
            incr FractureId
            # Set TipPoint
            set TipPoint [GiD_Geometry get point $PointId]
            AddPointToFracturesDict FracturesDict $FractureId $PointId $TipPoint TipPoint
            # Set BodyVolume
            AddBodyVolumeToFracturesDict FracturesDict $FractureId $BodyVolumesDict
            # Set LeftInterfaceVolume and RightInterfaceVolume
            set LeftInterfaceVolume [GiD_Geometry get volume [dict get $Info LeftInterfaceVolume]]
            AddInterfaceVolumeToFracturesDict FracturesDict $FractureId [dict get $Info LeftInterfaceVolume] $LeftInterfaceVolume LeftInterfaceVolume
            set RightInterfaceVolume [GiD_Geometry get volume [dict get $Info RightInterfaceVolume]]
            AddInterfaceVolumeToFracturesDict FracturesDict $FractureId [dict get $Info RightInterfaceVolume] $RightInterfaceVolume RightInterfaceVolume
            # Set LeftLine and RightLine
            set LeftLine [GiD_Geometry get line [dict get $Info LeftLine]]
            dict set FracturesDict $FractureId LeftLine Id [dict get $Info LeftLine]
            set RightLine [GiD_Geometry get line [dict get $Info RightLine]]
            dict set FracturesDict $FractureId RightLine Id [dict get $Info RightLine]
            # Set LeftPoint and RightPoint
            if {[lindex $LeftLine 2] eq [lindex $RightLine 2]} {
                set LeftPoint [GiD_Geometry get point [lindex $LeftLine 3]]
                AddPointToFracturesDict FracturesDict $FractureId [lindex $LeftLine 3] $LeftPoint LeftPoint
                set RightPoint [GiD_Geometry get point [lindex $RightLine 3]]
                AddPointToFracturesDict FracturesDict $FractureId [lindex $RightLine 3] $RightPoint RightPoint
            } elseif {[lindex $LeftLine 3] eq [lindex $RightLine 2]} {
                set LeftPoint [GiD_Geometry get point [lindex $LeftLine 2]]
                AddPointToFracturesDict FracturesDict $FractureId [lindex $LeftLine 2] $LeftPoint LeftPoint
                set RightPoint [GiD_Geometry get point [lindex $RightLine 3]]
                AddPointToFracturesDict FracturesDict $FractureId [lindex $RightLine 3] $RightPoint RightPoint
            } elseif {[lindex $LeftLine 2] eq [lindex $RightLine 3]} {
                set LeftPoint [GiD_Geometry get point [lindex $LeftLine 3]]
                AddPointToFracturesDict FracturesDict $FractureId [lindex $LeftLine 3] $LeftPoint LeftPoint
                set RightPoint [GiD_Geometry get point [lindex $RightLine 2]]
                AddPointToFracturesDict FracturesDict $FractureId [lindex $RightLine 2] $RightPoint RightPoint
            } else {
                set LeftPoint [GiD_Geometry get point [lindex $LeftLine 2]]
                AddPointToFracturesDict FracturesDict $FractureId [lindex $LeftLine 2] $LeftPoint LeftPoint
                set RightPoint [GiD_Geometry get point [lindex $RightLine 2]]
                AddPointToFracturesDict FracturesDict $FractureId [lindex $RightLine 2] $RightPoint RightPoint
            }
            # Set TopPoint, BotPoint, TopLine, BotLine, TopLeftSurface, TopRightSurface, BotLeftSurface and BotRightSurface
            set LeftSurface1 [GiD_Geometry get surface [lindex [lindex $LeftInterfaceVolume 2] 0]]
            set LeftSurface2 [GiD_Geometry get surface [lindex [lindex $LeftInterfaceVolume 3] 0]]
            set RightSurface1 [GiD_Geometry get surface [lindex [lindex $RightInterfaceVolume 2] 0]]
            set RightSurface2 [GiD_Geometry get surface [lindex [lindex $RightInterfaceVolume 3] 0]]
            set IsRightSurface1 0
            set IsTop 0
            for {set i 0} {$i < [lindex $LeftSurface1 2]} {incr i} {
                set LeftLineId [lindex [lindex $LeftSurface1 [expr { 9+$i }]] 0]
                for {set j 0} {$j < [lindex $RightSurface1 2]} {incr j} {
                    set RightLineId [lindex [lindex $RightSurface1 [expr { 9+$j }]] 0]
                    if {$LeftLineId eq $RightLineId} {
                        set IsRightSurface1 1
                        set AxisLineId $LeftLineId
                        set AxisLine [GiD_Geometry get line $AxisLineId]
                        if {[lindex $AxisLine 2] ne $PointId} {
                            set TopBotPointId [lindex $AxisLine 2]
                            set TopBotPoint [GiD_Geometry get point $TopBotPointId]
                        } else {
                            set TopBotPointId [lindex $AxisLine 3]
                            set TopBotPoint [GiD_Geometry get point $TopBotPointId]
                        }
                        set IsTop [IsTopPoint $TopBotPoint $LeftPoint $RightPoint $TipPoint]
                    }
                }
            }
            if {$IsRightSurface1 eq 0} {
                for {set i 0} {$i < [lindex $LeftSurface1 2]} {incr i} {
                    set LeftLineId [lindex [lindex $LeftSurface1 [expr { 9+$i }]] 0]
                    for {set j 0} {$j < [lindex $RightSurface2 2]} {incr j} {
                        set RightLineId [lindex [lindex $RightSurface2 [expr { 9+$j }]] 0]
                        if {$LeftLineId eq $RightLineId} {
                            set AxisLineId $LeftLineId
                            set AxisLine [GiD_Geometry get line $AxisLineId]
                            if {[lindex $AxisLine 2] ne $PointId} {
                                set TopBotPointId [lindex $AxisLine 2]
                                set TopBotPoint [GiD_Geometry get point $TopBotPointId]
                            } else {
                                set TopBotPointId [lindex $AxisLine 3]
                                set TopBotPoint [GiD_Geometry get point $TopBotPointId]
                            }
                            set IsTop [IsTopPoint $TopBotPoint $LeftPoint $RightPoint $TipPoint]
                        }
                    }
                }
            }
            if {$IsTop eq 1} {
                AddSurfaceToFracturesDict FracturesDict $FractureId [lindex [lindex $LeftInterfaceVolume 2] 0] $LeftSurface1 TopLeftSurface
                dict set FracturesDict $FractureId TopLine Id $AxisLineId
                AddPointToFracturesDict FracturesDict $FractureId $TopBotPointId $TopBotPoint TopPoint
                if {$IsRightSurface1 eq 1} {
                    AddSurfaceToFracturesDict FracturesDict $FractureId [lindex [lindex $RightInterfaceVolume 2] 0] $RightSurface1 TopRightSurface
                    SearchAxis AxisLineId AxisLine TopBotPointId TopBotPoint $LeftSurface2 $RightSurface2 $PointId
                    AddSurfaceToFracturesDict FracturesDict $FractureId [lindex [lindex $RightInterfaceVolume 3] 0] $RightSurface2 BotRightSurface
                } else {
                    AddSurfaceToFracturesDict FracturesDict $FractureId [lindex [lindex $RightInterfaceVolume 3] 0] $RightSurface2 TopRightSurface
                    SearchAxis AxisLineId AxisLine TopBotPointId TopBotPoint $LeftSurface2 $RightSurface1 $PointId
                    AddSurfaceToFracturesDict FracturesDict $FractureId [lindex [lindex $RightInterfaceVolume 2] 0] $RightSurface1 BotRightSurface
                }
                AddSurfaceToFracturesDict FracturesDict $FractureId [lindex [lindex $LeftInterfaceVolume 3] 0] $LeftSurface2 BotLeftSurface
                dict set FracturesDict $FractureId BotLine Id $AxisLineId
                AddPointToFracturesDict FracturesDict $FractureId $TopBotPointId $TopBotPoint BotPoint
            } else {
                AddSurfaceToFracturesDict FracturesDict $FractureId [lindex [lindex $LeftInterfaceVolume 2] 0] $LeftSurface1 BotLeftSurface
                dict set FracturesDict $FractureId BotLine Id $AxisLineId
                AddPointToFracturesDict FracturesDict $FractureId $TopBotPointId $TopBotPoint BotPoint
                if {$IsRightSurface1 eq 1} {
                    AddSurfaceToFracturesDict FracturesDict $FractureId [lindex [lindex $RightInterfaceVolume 2] 0] $RightSurface1 BotRightSurface
                    SearchAxis AxisLineId AxisLine TopBotPointId TopBotPoint $LeftSurface2 $RightSurface2 $PointId
                    AddSurfaceToFracturesDict FracturesDict $FractureId [lindex [lindex $RightInterfaceVolume 3] 0] $RightSurface2 TopRightSurface
                } else {
                    AddSurfaceToFracturesDict FracturesDict $FractureId [lindex [lindex $RightInterfaceVolume 3] 0] $RightSurface2 BotRightSurface
                    SearchAxis AxisLineId AxisLine TopBotPointId TopBotPoint $LeftSurface2 $RightSurface1 $PointId
                    AddSurfaceToFracturesDict FracturesDict $FractureId [lindex [lindex $RightInterfaceVolume 2] 0] $RightSurface1 TopRightSurface
                }
                AddSurfaceToFracturesDict FracturesDict $FractureId [lindex [lindex $LeftInterfaceVolume 3] 0] $LeftSurface2 TopLeftSurface
                dict set FracturesDict $FractureId TopLine Id $AxisLineId
                AddPointToFracturesDict FracturesDict $FractureId $TopBotPointId $TopBotPoint TopPoint
            }
        }
    }
    
    ## Start FracturesData.json file
    set filename [file join $dir FracturesData.json]
    set FileVar1 [open $filename w]
        
    puts $FileVar1 "\{"

    ## fracture_data
    puts $FileVar1 "    \"fracture_data\": \{"
    puts $FileVar1 "        \"gid_path\":                             \"$gidpath\","
    # propagation parameters and body_domain_sub_model_part_list
    WriteFractureData FileVar1
    # interface_domain_sub_model_part_list
    set PutStrings \[
    set Groups [GiD_Info conditions Interface_Part groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
        append PutStrings \" [lindex [lindex $Groups $i] 1] \" ,
    }
    if {[GiD_Groups exists PropagationUnion_3d_6] eq 1} {
        append PutStrings \" PropagationUnion_3d_6 \" \]
    } else {
        set PutStrings [string trimright $PutStrings ,]
        append PutStrings \]
    }
    puts $FileVar1 "        \"interface_domain_sub_model_part_list\": $PutStrings"
    puts $FileVar1 "    \},"
    
    ## body_volumes_list
    WriteBodyVolumesList FileVar1 $BodyVolumesDict
        
    ## fractures_list
    WriteFracturesList FileVar1 $FracturesDict
    
    puts $FileVar1 "\}"
    
    close $FileVar1
}

proc GenerateNewFractures { dir problemtypedir PropagationData } {

    ## Source auxiliar procedures
    source [file join $problemtypedir FracturePropagation3DAuxProcs.tcl]

    # Previous Fractures dictionaries
    set BodyVolumesDict [lindex $PropagationData 1]
    set FracturesDict [lindex $PropagationData 2]
    
    # Propagation and Bifurcation dictionaries
    set PropagationDict [lindex $PropagationData 3]
    set BifurcationDict [lindex $PropagationData 4]

    # interface_domain_sub_model_part_old_list
    set InterfaceGroupsOld \[
    set Groups [GiD_Info conditions Interface_Part groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
        append InterfaceGroupsOld \" [lindex [lindex $Groups $i] 1] \" ,
    }
    if {[GiD_Groups exists PropagationUnion_3d_6] eq 1} {
        append InterfaceGroupsOld \" PropagationUnion_3d_6 \" \]
    } else {
        set InterfaceGroupsOld [string trimright $InterfaceGroupsOld ,]
        append InterfaceGroupsOld \]
        # Create PropagationUnion group if it does not exist
        GiD_Groups create PropagationUnion_3d_6
    }    
    set NumPropUnionGroups [llength [GiD_Groups list PropagationUnion_3d_6]]
    
    # create link interface group if bifurcation occurs and it does not exist
    if {[dict size $BifurcationDict] > 0} {
        set LinkInterfaceGroup 0
        #set Groups [GiD_Info conditions Interface_Part groups]
        for {set i 0} {$i < [llength $Groups]} {incr i} {
            if {[lindex [lindex $Groups $i] 3] eq true} {
                set LinkInterfaceGroup [lindex [lindex $Groups $i] 1]
                break
            }
        }
        if {$LinkInterfaceGroup eq 0} {
            GiD_Groups create Bifurcation_Link_Part_9
            set LinkInterfaceGroup Bifurcation_Link_Part_9
            set ConditionValues "true [lindex [lindex $Groups 0] 4] [lindex [lindex $Groups 0] 5] [lindex [lindex $Groups 0] 6] [lindex [lindex $Groups 0] 7] \
            [lindex [lindex $Groups 0] 8] [lindex [lindex $Groups 0] 9] [lindex [lindex $Groups 0] 10] [lindex [lindex $Groups 0] 11] [lindex [lindex $Groups 0] 12]\
            [lindex [lindex $Groups 0] 13] [lindex [lindex $Groups 0] 14] [lindex [lindex $Groups 0] 15] [lindex [lindex $Groups 0] 16] [lindex [lindex $Groups 0] 17]\
            0.0 [lindex [lindex $Groups 0] 19] [lindex [lindex $Groups 0] 20] [lindex [lindex $Groups 0] 21]"
            GiD_AssignData condition Interface_Part groups $ConditionValues $LinkInterfaceGroup
        }
    }
    
    # Delete all BodyVolumes
    dict for {BodyId BodyVolume} $BodyVolumesDict {
        GiD_Process Mescape Geometry Delete Volumes $BodyId escape
    }

    # Loop through all Propagation Fractures
    dict for {PropId Propagation} $PropagationDict {
        set MotherFractureId [dict get $Propagation MotherFractureId]

        set BodyVolumeId [lindex [dict get $FracturesDict $MotherFractureId BodyVolumes] 0]
        set BodyVolumeSurfaces [dict get $BodyVolumesDict $BodyVolumeId Surfaces]
        set OldTopLeftSurfaceLines [dict get $FracturesDict $MotherFractureId TopLeftSurface Lines]
        set OldTopRightSurfaceLines [dict get $FracturesDict $MotherFractureId TopRightSurface Lines]
        set OldBotLeftSurfaceLines [dict get $FracturesDict $MotherFractureId BotLeftSurface Lines]
        set OldBotRightSurfaceLines [dict get $FracturesDict $MotherFractureId BotRightSurface Lines]
        
        ## Modify Geometry
        # Delete InterfaceVolume
        GiD_Process Mescape Geometry Delete Volumes [dict get $FracturesDict $MotherFractureId LeftInterfaceVolume Id] [dict get $FracturesDict $MotherFractureId RightInterfaceVolume Id] escape
        # Delete old crack surfaces
        GiD_Process Mescape Geometry Delete Surfaces [dict get $FracturesDict $MotherFractureId TopLeftSurface Id] [dict get $FracturesDict $MotherFractureId TopRightSurface Id] \
            [dict get $FracturesDict $MotherFractureId BotLeftSurface Id] [dict get $FracturesDict $MotherFractureId BotRightSurface Id] escape
        set Index [lsearch $BodyVolumeSurfaces [dict get $FracturesDict $MotherFractureId TopLeftSurface Id]]
        set BodyVolumeSurfaces [lreplace $BodyVolumeSurfaces $Index $Index]
        set Index [lsearch $BodyVolumeSurfaces [dict get $FracturesDict $MotherFractureId TopRightSurface Id]]
        set BodyVolumeSurfaces [lreplace $BodyVolumeSurfaces $Index $Index]
        set Index [lsearch $BodyVolumeSurfaces [dict get $FracturesDict $MotherFractureId BotLeftSurface Id]]
        set BodyVolumeSurfaces [lreplace $BodyVolumeSurfaces $Index $Index]
        set Index [lsearch $BodyVolumeSurfaces [dict get $FracturesDict $MotherFractureId BotRightSurface Id]]
        set BodyVolumeSurfaces [lreplace $BodyVolumeSurfaces $Index $Index]
        # Delete old crack lines
        GiD_Process Mescape Geometry Delete Lines [dict get $FracturesDict $MotherFractureId TopLine Id] [dict get $FracturesDict $MotherFractureId BotLine Id] \
            [dict get $FracturesDict $MotherFractureId LeftLine Id] [dict get $FracturesDict $MotherFractureId RightLine Id] escape
        set Index [lsearch $OldTopLeftSurfaceLines [dict get $FracturesDict $MotherFractureId TopLine Id]]
        set OldTopLeftSurfaceLines [lreplace $OldTopLeftSurfaceLines $Index $Index]
        set Index [lsearch $OldTopLeftSurfaceLines [dict get $FracturesDict $MotherFractureId LeftLine Id]]
        set OldTopLeftSurfaceLines [lreplace $OldTopLeftSurfaceLines $Index $Index]
        set Index [lsearch $OldTopRightSurfaceLines [dict get $FracturesDict $MotherFractureId TopLine Id]]
        set OldTopRightSurfaceLines [lreplace $OldTopRightSurfaceLines $Index $Index]
        set Index [lsearch $OldTopRightSurfaceLines [dict get $FracturesDict $MotherFractureId RightLine Id]]
        set OldTopRightSurfaceLines [lreplace $OldTopRightSurfaceLines $Index $Index]
        set Index [lsearch $OldBotLeftSurfaceLines [dict get $FracturesDict $MotherFractureId BotLine Id]]
        set OldBotLeftSurfaceLines [lreplace $OldBotLeftSurfaceLines $Index $Index]
        set Index [lsearch $OldBotLeftSurfaceLines [dict get $FracturesDict $MotherFractureId LeftLine Id]]
        set OldBotLeftSurfaceLines [lreplace $OldBotLeftSurfaceLines $Index $Index]
        set Index [lsearch $OldBotRightSurfaceLines [dict get $FracturesDict $MotherFractureId BotLine Id]]
        set OldBotRightSurfaceLines [lreplace $OldBotRightSurfaceLines $Index $Index]
        set Index [lsearch $OldBotRightSurfaceLines [dict get $FracturesDict $MotherFractureId RightLine Id]]
        set OldBotRightSurfaceLines [lreplace $OldBotRightSurfaceLines $Index $Index]
        # Delete old TipPoint
        GiD_Process Mescape Geometry Delete Points [dict get $FracturesDict $MotherFractureId TipPoint Id] escape
        
        # Create new points
        GiD_Process Mescape Geometry Create Point
        # Create new point in TopInitCoordinates location
        GiD_Process [lindex [dict get $Propagation TopInitCoordinates] 0] [lindex [dict get $Propagation TopInitCoordinates] 1] [lindex [dict get $Propagation TopInitCoordinates] 2]
        # Create new point in BotInitCoordinates location
        GiD_Process [lindex [dict get $Propagation BotInitCoordinates] 0] [lindex [dict get $Propagation BotInitCoordinates] 1] [lindex [dict get $Propagation BotInitCoordinates] 2]
        # Create new point in LeftInitCoordinates location
        GiD_Process [lindex [dict get $Propagation LeftInitCoordinates] 0] [lindex [dict get $Propagation LeftInitCoordinates] 1] [lindex [dict get $Propagation LeftInitCoordinates] 2]
        # Create new point in RightInitCoordinates location
        GiD_Process [lindex [dict get $Propagation RightInitCoordinates] 0] [lindex [dict get $Propagation RightInitCoordinates] 1] [lindex [dict get $Propagation RightInitCoordinates] 2]
        # Create new point in TopEndCoordinates location
        GiD_Process [lindex [dict get $Propagation TopEndCoordinates] 0] [lindex [dict get $Propagation TopEndCoordinates] 1] [lindex [dict get $Propagation TopEndCoordinates] 2]
        # Create new point in BotEndCoordinates location
        GiD_Process [lindex [dict get $Propagation BotEndCoordinates] 0] [lindex [dict get $Propagation BotEndCoordinates] 1] [lindex [dict get $Propagation BotEndCoordinates] 2]
        # Create new point in LeftEndCoordinates location
        GiD_Process [lindex [dict get $Propagation LeftEndCoordinates] 0] [lindex [dict get $Propagation LeftEndCoordinates] 1] [lindex [dict get $Propagation LeftEndCoordinates] 2]
        # Create new point in RightEndCoordinates location
        GiD_Process [lindex [dict get $Propagation RightEndCoordinates] 0] [lindex [dict get $Propagation RightEndCoordinates] 1] [lindex [dict get $Propagation RightEndCoordinates] 2]
        # Create new point in TipCoordinates location
        GiD_Process [lindex [dict get $Propagation TipCoordinates] 0] [lindex [dict get $Propagation TipCoordinates] 1] [lindex [dict get $Propagation TipCoordinates] 2]
        GiD_Process escape
        
        # Regenerate lines for old BotLeftSurface
        GiD_Process Mescape Geometry Create Line Join [dict get $FracturesDict $MotherFractureId BotPoint Id] [expr { [GiD_Info Geometry MaxNumPoints]-7 }] \
            [expr { [GiD_Info Geometry MaxNumPoints]-6 }] [dict get $FracturesDict $MotherFractureId LeftPoint Id] escape escape
        # Regenerate old BotLeftSurface with its normal pointing outside the contact volume
        GiD_Process Mescape Geometry Create NurbsSurface [GiD_Info Geometry MaxNumLines] [expr {[GiD_Info Geometry MaxNumLines]-1}] [expr {[GiD_Info Geometry MaxNumLines]-2}]
        for {set i 0} {$i < [llength $OldBotLeftSurfaceLines]} {incr i} {
            GiD_Process [lindex $OldBotLeftSurfaceLines $i]
        }
        GiD_Process escape escape
        GiD_Process Mescape utilities SwapNormals Surfaces SelByNormal [GiD_Info Geometry MaxNumSurfaces] \
            escape [ComputeNormal [dict get $FracturesDict $MotherFractureId BotPoint Coordinates] [dict get $FracturesDict $MotherFractureId LeftPoint Coordinates] [dict get $Propagation BotInitCoordinates]] Yes escape
        lappend BodyVolumeSurfaces [GiD_Info Geometry MaxNumSurfaces]
        # Regenerate lines for old BotRightSurface
        GiD_Process Mescape Geometry Create Line Join [dict get $FracturesDict $MotherFractureId RightPoint Id] [expr { [GiD_Info Geometry MaxNumPoints]-5 }] \
            [expr { [GiD_Info Geometry MaxNumPoints]-7 }] escape escape
        # Regenerate old BotRightSurface with its normal pointing outside the contact volume
        GiD_Process Mescape Geometry Create NurbsSurface [GiD_Info Geometry MaxNumLines] [expr {[GiD_Info Geometry MaxNumLines]-1}] [expr {[GiD_Info Geometry MaxNumLines]-4}]
        for {set i 0} {$i < [llength $OldBotRightSurfaceLines]} {incr i} {
            GiD_Process [lindex $OldBotRightSurfaceLines $i]
        }
        GiD_Process escape escape
        GiD_Process Mescape utilities SwapNormals Surfaces SelByNormal [GiD_Info Geometry MaxNumSurfaces] \
            escape [ComputeNormal [dict get $FracturesDict $MotherFractureId BotPoint Coordinates] [dict get $Propagation BotInitCoordinates] [dict get $FracturesDict $MotherFractureId RightPoint Coordinates]] Yes escape
        lappend BodyVolumeSurfaces [GiD_Info Geometry MaxNumSurfaces]
        # Regenerate lines for old TopLeftSurface
        GiD_Process Mescape Geometry Create Line Join [dict get $FracturesDict $MotherFractureId TopPoint Id] [expr { [GiD_Info Geometry MaxNumPoints]-8 }] \
            [expr { [GiD_Info Geometry MaxNumPoints]-6 }] escape escape
        # Regenerate old TopLeftSurface with its normal pointing outside the contact volume
        GiD_Process Mescape Geometry Create NurbsSurface [GiD_Info Geometry MaxNumLines] [expr {[GiD_Info Geometry MaxNumLines]-1}] [expr {[GiD_Info Geometry MaxNumLines]-4}]
        for {set i 0} {$i < [llength $OldTopLeftSurfaceLines]} {incr i} {
            GiD_Process [lindex $OldTopLeftSurfaceLines $i]
        }
        GiD_Process escape escape
        GiD_Process Mescape utilities SwapNormals Surfaces SelByNormal [GiD_Info Geometry MaxNumSurfaces] \
            escape [ComputeNormal [dict get $FracturesDict $MotherFractureId TopPoint Coordinates] [dict get $Propagation TopInitCoordinates] [dict get $FracturesDict $MotherFractureId LeftPoint Coordinates]] Yes escape
        lappend BodyVolumeSurfaces [GiD_Info Geometry MaxNumSurfaces]
        # Regenerate lines for old TopRightSurface
        GiD_Process Mescape Geometry Create Line Join [expr { [GiD_Info Geometry MaxNumPoints]-8 }] [expr { [GiD_Info Geometry MaxNumPoints]-5 }] escape escape
        # Regenerate old TopRightSurface with its normal pointing outside the contact volume
        GiD_Process Mescape Geometry Create NurbsSurface [GiD_Info Geometry MaxNumLines] [expr {[GiD_Info Geometry MaxNumLines]-2}] [expr {[GiD_Info Geometry MaxNumLines]-4}]
        for {set i 0} {$i < [llength $OldTopRightSurfaceLines]} {incr i} {
            GiD_Process [lindex $OldTopRightSurfaceLines $i]
        }
        GiD_Process escape escape
        GiD_Process Mescape utilities SwapNormals Surfaces SelByNormal [GiD_Info Geometry MaxNumSurfaces] \
            escape [ComputeNormal [dict get $FracturesDict $MotherFractureId TopPoint Coordinates] [dict get $FracturesDict $MotherFractureId RightPoint Coordinates] [dict get $Propagation TopInitCoordinates]] Yes escape
        lappend BodyVolumeSurfaces [GiD_Info Geometry MaxNumSurfaces]
        # Regenerate old LeftInterfaceVolume
        set Surf1 [list [expr {[GiD_Info Geometry MaxNumSurfaces]-3}] 1]
        set Surf2 [list [expr {[GiD_Info Geometry MaxNumSurfaces]-1}] 0]
        # Note: the orientation 1 for the first surface means that its normal points outside the contact volume
        # Note: the orientation 0 for the second surface means that its normal points outside the contact volume
        set TransformMatrix [ComputeTransformMatrix [dict get $FracturesDict $MotherFractureId LeftPoint Coordinates] [dict get $FracturesDict $MotherFractureId BotPoint Coordinates] \
            [dict get $FracturesDict $MotherFractureId TopPoint Coordinates] [dict get $Propagation LeftInitCoordinates] [dict get $FracturesDict $MotherFractureId LeftPoint Coordinates]]
        GiD_Geometry create volume [dict get $FracturesDict $MotherFractureId LeftInterfaceVolume Id] [dict get $FracturesDict $MotherFractureId LeftInterfaceVolume Layer] 2 \
            $Surf1 $Surf2 contactvolume $TransformMatrix
        # Regenerate old RightInterfaceVolume
        set Surf1 [list [expr {[GiD_Info Geometry MaxNumSurfaces]-2}] 1]
        set Surf2 [list [GiD_Info Geometry MaxNumSurfaces] 0]
        set TransformMatrix [ComputeTransformMatrix [dict get $FracturesDict $MotherFractureId RightPoint Coordinates] [dict get $FracturesDict $MotherFractureId BotPoint Coordinates] \
            [dict get $FracturesDict $MotherFractureId TopPoint Coordinates] [dict get $FracturesDict $MotherFractureId RightPoint Coordinates] [dict get $Propagation RightInitCoordinates]]
        GiD_Geometry create volume [dict get $FracturesDict $MotherFractureId RightInterfaceVolume Id] [dict get $FracturesDict $MotherFractureId RightInterfaceVolume Layer] 2 \
            $Surf1 $Surf2 contactvolume $TransformMatrix
        
        # Create lines of PropagationUnion
        GiD_Process Mescape Geometry Create Line Join [expr {[GiD_Info Geometry MaxNumPoints]-7}] [expr {[GiD_Info Geometry MaxNumPoints]-3}] escape \
            [expr {[GiD_Info Geometry MaxNumPoints]-6}] [expr {[GiD_Info Geometry MaxNumPoints]-2}] escape \
            [expr {[GiD_Info Geometry MaxNumPoints]-8}] [expr {[GiD_Info Geometry MaxNumPoints]-4}] escape \
            [expr {[GiD_Info Geometry MaxNumPoints]-5}] [expr {[GiD_Info Geometry MaxNumPoints]-1}] escape \
            [expr {[GiD_Info Geometry MaxNumPoints]-7}] [expr {[GiD_Info Geometry MaxNumPoints]-1}] [expr {[GiD_Info Geometry MaxNumPoints]-8}] \
            [expr {[GiD_Info Geometry MaxNumPoints]-2}] [expr {[GiD_Info Geometry MaxNumPoints]-7}] escape \
            [expr {[GiD_Info Geometry MaxNumPoints]-3}] [expr {[GiD_Info Geometry MaxNumPoints]-1}] \
            [expr {[GiD_Info Geometry MaxNumPoints]-4}] [expr {[GiD_Info Geometry MaxNumPoints]-2}] [expr {[GiD_Info Geometry MaxNumPoints]-3}] escape escape
        # Create surfaces of PropagationUnion
        GiD_Process Mescape Geometry Create NurbsSurface [expr {[GiD_Info Geometry MaxNumLines]-15}] [expr {[GiD_Info Geometry MaxNumLines]-8}] [expr {[GiD_Info Geometry MaxNumLines]-7}] escape \
            [expr {[GiD_Info Geometry MaxNumLines]-12}] [expr {[GiD_Info Geometry MaxNumLines]-8}] [expr {[GiD_Info Geometry MaxNumLines]-6}] escape \
            [expr {[GiD_Info Geometry MaxNumLines]-11}] [expr {[GiD_Info Geometry MaxNumLines]-7}] [expr {[GiD_Info Geometry MaxNumLines]-3}] escape \
            [expr {[GiD_Info Geometry MaxNumLines]-9}] [expr {[GiD_Info Geometry MaxNumLines]-6}] [expr {[GiD_Info Geometry MaxNumLines]-2}] escape \
            [expr {[GiD_Info Geometry MaxNumLines]-18}] [expr {[GiD_Info Geometry MaxNumLines]-10}] [expr {[GiD_Info Geometry MaxNumLines]-4}] escape \
            [expr {[GiD_Info Geometry MaxNumLines]-13}] [expr {[GiD_Info Geometry MaxNumLines]-10}] [expr {[GiD_Info Geometry MaxNumLines]-5}] escape \
            [expr {[GiD_Info Geometry MaxNumLines]-11}] [expr {[GiD_Info Geometry MaxNumLines]-4}] [GiD_Info Geometry MaxNumLines] escape \
            [expr {[GiD_Info Geometry MaxNumLines]-9}] [expr {[GiD_Info Geometry MaxNumLines]-5}] [expr {[GiD_Info Geometry MaxNumLines]-1}] escape escape
        lappend BodyVolumeSurfaces [GiD_Info Geometry MaxNumSurfaces]
        lappend BodyVolumeSurfaces [expr {[GiD_Info Geometry MaxNumSurfaces]-1}]
        lappend BodyVolumeSurfaces [expr {[GiD_Info Geometry MaxNumSurfaces]-2}]
        lappend BodyVolumeSurfaces [expr {[GiD_Info Geometry MaxNumSurfaces]-3}]
        lappend BodyVolumeSurfaces [expr {[GiD_Info Geometry MaxNumSurfaces]-4}]
        lappend BodyVolumeSurfaces [expr {[GiD_Info Geometry MaxNumSurfaces]-5}]
        lappend BodyVolumeSurfaces [expr {[GiD_Info Geometry MaxNumSurfaces]-6}]
        lappend BodyVolumeSurfaces [expr {[GiD_Info Geometry MaxNumSurfaces]-7}]
        # Create PropagationUnion subgroups and assign points
        GiD_EntitiesGroups assign PropagationUnion_3d_6 points [expr {[GiD_Info Geometry MaxNumPoints]-7}]
        GiD_EntitiesGroups assign PropagationUnion_3d_6 points [expr {[GiD_Info Geometry MaxNumPoints]-5}]
        GiD_EntitiesGroups assign PropagationUnion_3d_6 points [expr {[GiD_Info Geometry MaxNumPoints]-1}]
        GiD_EntitiesGroups assign PropagationUnion_3d_6 points [expr {[GiD_Info Geometry MaxNumPoints]-8}]
        GiD_EntitiesGroups assign PropagationUnion_3d_6 points [expr {[GiD_Info Geometry MaxNumPoints]-3}]
        GiD_EntitiesGroups assign PropagationUnion_3d_6 points [expr {[GiD_Info Geometry MaxNumPoints]-4}]
        GiD_EntitiesGroups assign PropagationUnion_3d_6 points [expr {[GiD_Info Geometry MaxNumPoints]-6}]
        GiD_EntitiesGroups assign PropagationUnion_3d_6 points [expr {[GiD_Info Geometry MaxNumPoints]-2}]
        # Prism 1
        AddPropagationUnionPoint NumPropUnionGroups [expr {[GiD_Info Geometry MaxNumPoints]-7}]
        AddPropagationUnionPoint NumPropUnionGroups [expr {[GiD_Info Geometry MaxNumPoints]-5}]
        AddPropagationUnionPoint NumPropUnionGroups [expr {[GiD_Info Geometry MaxNumPoints]-1}]
        AddPropagationUnionPoint NumPropUnionGroups [expr {[GiD_Info Geometry MaxNumPoints]-8}]
        AddPropagationUnionPoint NumPropUnionGroups [expr {[GiD_Info Geometry MaxNumPoints]-5}]
        AddPropagationUnionPoint NumPropUnionGroups [expr {[GiD_Info Geometry MaxNumPoints]-1}]
        # Prism 2
        AddPropagationUnionPoint NumPropUnionGroups [expr {[GiD_Info Geometry MaxNumPoints]-7}]
        AddPropagationUnionPoint NumPropUnionGroups [expr {[GiD_Info Geometry MaxNumPoints]-1}]
        AddPropagationUnionPoint NumPropUnionGroups [expr {[GiD_Info Geometry MaxNumPoints]-3}]
        AddPropagationUnionPoint NumPropUnionGroups [expr {[GiD_Info Geometry MaxNumPoints]-8}]
        AddPropagationUnionPoint NumPropUnionGroups [expr {[GiD_Info Geometry MaxNumPoints]-1}]
        AddPropagationUnionPoint NumPropUnionGroups [expr {[GiD_Info Geometry MaxNumPoints]-4}]
        # Prism 3
        AddPropagationUnionPoint NumPropUnionGroups [expr {[GiD_Info Geometry MaxNumPoints]-6}]
        AddPropagationUnionPoint NumPropUnionGroups [expr {[GiD_Info Geometry MaxNumPoints]-7}]
        AddPropagationUnionPoint NumPropUnionGroups [expr {[GiD_Info Geometry MaxNumPoints]-2}]
        AddPropagationUnionPoint NumPropUnionGroups [expr {[GiD_Info Geometry MaxNumPoints]-6}]
        AddPropagationUnionPoint NumPropUnionGroups [expr {[GiD_Info Geometry MaxNumPoints]-8}]
        AddPropagationUnionPoint NumPropUnionGroups [expr {[GiD_Info Geometry MaxNumPoints]-2}]
        # Prism 4
        AddPropagationUnionPoint NumPropUnionGroups [expr {[GiD_Info Geometry MaxNumPoints]-2}]
        AddPropagationUnionPoint NumPropUnionGroups [expr {[GiD_Info Geometry MaxNumPoints]-7}]
        AddPropagationUnionPoint NumPropUnionGroups [expr {[GiD_Info Geometry MaxNumPoints]-3}]
        AddPropagationUnionPoint NumPropUnionGroups [expr {[GiD_Info Geometry MaxNumPoints]-2}]
        AddPropagationUnionPoint NumPropUnionGroups [expr {[GiD_Info Geometry MaxNumPoints]-8}]
        AddPropagationUnionPoint NumPropUnionGroups [expr {[GiD_Info Geometry MaxNumPoints]-4}]
        
        # Create lines for new crack tip
        GiD_Process Mescape Geometry Create Line Join [expr {[GiD_Info Geometry MaxNumPoints]-4}] [GiD_Info Geometry MaxNumPoints] escape \
            [expr {[GiD_Info Geometry MaxNumPoints]-3}] [GiD_Info Geometry MaxNumPoints] escape \
            [expr {[GiD_Info Geometry MaxNumPoints]-2}] [GiD_Info Geometry MaxNumPoints] escape \
            [expr {[GiD_Info Geometry MaxNumPoints]-1}] [GiD_Info Geometry MaxNumPoints] escape escape
        # Create surfaces for new crack tip
        # New BotLeftSurface
        GiD_Process Mescape Geometry Create NurbsSurface [expr {[GiD_Info Geometry MaxNumLines]-4}] [expr {[GiD_Info Geometry MaxNumLines]-2}] [expr {[GiD_Info Geometry MaxNumLines]-1}] escape escape
        GiD_Process Mescape utilities SwapNormals Surfaces SelByNormal [GiD_Info Geometry MaxNumSurfaces] \
            escape [ComputeNormal [dict get $Propagation LeftEndCoordinates] [dict get $Propagation TipCoordinates] [dict get $Propagation BotEndCoordinates]] Yes escape
        lappend BodyVolumeSurfaces [GiD_Info Geometry MaxNumSurfaces]
        # New TopLeftSurface
        GiD_Process Mescape Geometry Create NurbsSurface [expr {[GiD_Info Geometry MaxNumLines]-5}] [expr {[GiD_Info Geometry MaxNumLines]-1}] [expr {[GiD_Info Geometry MaxNumLines]-3}] escape escape
        GiD_Process Mescape utilities SwapNormals Surfaces SelByNormal [GiD_Info Geometry MaxNumSurfaces] \
            escape [ComputeNormal [dict get $Propagation LeftEndCoordinates] [dict get $Propagation TopEndCoordinates] [dict get $Propagation TipCoordinates]] Yes escape
        lappend BodyVolumeSurfaces [GiD_Info Geometry MaxNumSurfaces]
        # New BotRightSurface
        GiD_Process Mescape Geometry Create NurbsSurface [expr {[GiD_Info Geometry MaxNumLines]-7}] [expr {[GiD_Info Geometry MaxNumLines]-2}] [GiD_Info Geometry MaxNumLines] escape escape
        GiD_Process Mescape utilities SwapNormals Surfaces SelByNormal [GiD_Info Geometry MaxNumSurfaces] \
            escape [ComputeNormal [dict get $Propagation RightEndCoordinates] [dict get $Propagation BotEndCoordinates] [dict get $Propagation TipCoordinates]] Yes escape
        lappend BodyVolumeSurfaces [GiD_Info Geometry MaxNumSurfaces]
        # New TopRightSurface
        GiD_Process Mescape Geometry Create NurbsSurface [expr {[GiD_Info Geometry MaxNumLines]-6}] [expr {[GiD_Info Geometry MaxNumLines]-3}] [GiD_Info Geometry MaxNumLines] escape escape
        GiD_Process Mescape utilities SwapNormals Surfaces SelByNormal [GiD_Info Geometry MaxNumSurfaces] \
            escape [ComputeNormal [dict get $Propagation RightEndCoordinates] [dict get $Propagation TipCoordinates] [dict get $Propagation TopEndCoordinates]] Yes escape
        lappend BodyVolumeSurfaces [GiD_Info Geometry MaxNumSurfaces]
        # Create new LeftInterfaceVolume
        set Surf1 [list [expr {[GiD_Info Geometry MaxNumSurfaces]-3}] 1]
        set Surf2 [list [expr {[GiD_Info Geometry MaxNumSurfaces]-2}] 0]
        set TransformMatrix [ComputeTransformMatrix [dict get $Propagation LeftEndCoordinates] [dict get $Propagation BotEndCoordinates] \
            [dict get $Propagation TopEndCoordinates] [dict get $Propagation TipCoordinates] [dict get $Propagation LeftEndCoordinates]]
        GiD_Geometry create volume append [dict get $FracturesDict $MotherFractureId LeftInterfaceVolume Layer] 2 \
            $Surf1 $Surf2 contactvolume $TransformMatrix
        # Create new RightInterfaceVolume
        set Surf1 [list [expr {[GiD_Info Geometry MaxNumSurfaces]-1}] 1]
        set Surf2 [list [GiD_Info Geometry MaxNumSurfaces] 0]
        set TransformMatrix [ComputeTransformMatrix [dict get $Propagation RightEndCoordinates] [dict get $Propagation BotEndCoordinates] \
            [dict get $Propagation TopEndCoordinates] [dict get $Propagation RightEndCoordinates] [dict get $Propagation TipCoordinates]]
        GiD_Geometry create volume append [dict get $FracturesDict $MotherFractureId RightInterfaceVolume Layer] 2 \
            $Surf1 $Surf2 contactvolume $TransformMatrix
        
        
        ## Set Conditions
        for {set i 0} {$i < [llength [dict get $FracturesDict $MotherFractureId LeftInterfaceVolume Groups]]} {incr i} {
            GiD_EntitiesGroups assign [lindex [dict get $FracturesDict $MotherFractureId LeftInterfaceVolume Groups] $i] volumes [dict get $FracturesDict $MotherFractureId LeftInterfaceVolume Id]
            GiD_EntitiesGroups assign [lindex [dict get $FracturesDict $MotherFractureId LeftInterfaceVolume Groups] $i] volumes [expr {[GiD_Info Geometry MaxNumVolumes]-1}]
        }
        for {set i 0} {$i < [llength [dict get $FracturesDict $MotherFractureId RightInterfaceVolume Groups]]} {incr i} {
            GiD_EntitiesGroups assign [lindex [dict get $FracturesDict $MotherFractureId RightInterfaceVolume Groups] $i] volumes [dict get $FracturesDict $MotherFractureId RightInterfaceVolume Id]
            GiD_EntitiesGroups assign [lindex [dict get $FracturesDict $MotherFractureId RightInterfaceVolume Groups] $i] volumes [GiD_Info Geometry MaxNumVolumes]
        }
        
        
        ## Set Mesh options
        # GiD_Process Mescape Meshing Structured Lines 2 [expr {[GiD_Info Geometry MaxNumLines]-23}] \
        #     [expr {[GiD_Info Geometry MaxNumLines]-21}] [expr {[GiD_Info Geometry MaxNumLines]-20}] [expr {[GiD_Info Geometry MaxNumLines]-18}] \
        #     [expr {[GiD_Info Geometry MaxNumLines]-3}] [expr {[GiD_Info Geometry MaxNumLines]-2}] [expr {[GiD_Info Geometry MaxNumLines]-1}] \
        #     [GiD_Info Geometry MaxNumLines] escape escape
        GiD_Process Mescape Meshing Structured Surfaces [expr {[GiD_Info Geometry MaxNumSurfaces]-11}] \
            [expr {[GiD_Info Geometry MaxNumSurfaces]-10}] [expr {[GiD_Info Geometry MaxNumSurfaces]-9}] [expr {[GiD_Info Geometry MaxNumSurfaces]-8}] \
            [expr {[GiD_Info Geometry MaxNumSurfaces]-7}] [expr {[GiD_Info Geometry MaxNumSurfaces]-6}] [expr {[GiD_Info Geometry MaxNumSurfaces]-5}] \
            [expr {[GiD_Info Geometry MaxNumSurfaces]-4}] escape 1 [expr {[GiD_Info Geometry MaxNumLines]-12}] escape escape


        ## Update dictionaries
        # FracturesDict
        dict set FracturesDict $MotherFractureId TipPoint Id [GiD_Info Geometry MaxNumPoints]
        dict set FracturesDict $MotherFractureId TipPoint Coordinates [dict get $Propagation TipCoordinates]
        dict set FracturesDict $MotherFractureId TopPoint Id [expr {[GiD_Info Geometry MaxNumPoints]-4}]
        dict set FracturesDict $MotherFractureId TopPoint Coordinates [dict get $Propagation TopEndCoordinates]
        dict set FracturesDict $MotherFractureId BotPoint Id [expr {[GiD_Info Geometry MaxNumPoints]-3}]
        dict set FracturesDict $MotherFractureId BotPoint Coordinates [dict get $Propagation BotEndCoordinates]
        dict set FracturesDict $MotherFractureId LeftPoint Id [expr {[GiD_Info Geometry MaxNumPoints]-2}]
        dict set FracturesDict $MotherFractureId LeftPoint Coordinates [dict get $Propagation LeftEndCoordinates]
        dict set FracturesDict $MotherFractureId RightPoint Id [expr {[GiD_Info Geometry MaxNumPoints]-1}]
        dict set FracturesDict $MotherFractureId RightPoint Coordinates [dict get $Propagation RightEndCoordinates]
        dict set FracturesDict $MotherFractureId TopLine Id [expr {[GiD_Info Geometry MaxNumLines]-3}]
        dict set FracturesDict $MotherFractureId BotLine Id [expr {[GiD_Info Geometry MaxNumLines]-2}]
        dict set FracturesDict $MotherFractureId LeftLine Id [expr {[GiD_Info Geometry MaxNumLines]-1}]
        dict set FracturesDict $MotherFractureId RightLine Id [GiD_Info Geometry MaxNumLines]
        dict set FracturesDict $MotherFractureId TopLeftSurface Id [expr {[GiD_Info Geometry MaxNumSurfaces]-2}]
        set Lines [list [expr {[GiD_Info Geometry MaxNumLines]-1}] [expr {[GiD_Info Geometry MaxNumLines]-3}] [expr {[GiD_Info Geometry MaxNumLines]-5}]]
        dict set FracturesDict $MotherFractureId TopLeftSurface Lines $Lines
        dict set FracturesDict $MotherFractureId TopRightSurface Id [GiD_Info Geometry MaxNumSurfaces]
        set Lines [list [GiD_Info Geometry MaxNumLines] [expr {[GiD_Info Geometry MaxNumLines]-3}] [expr {[GiD_Info Geometry MaxNumLines]-6}]]
        dict set FracturesDict $MotherFractureId TopRightSurface Lines $Lines
        dict set FracturesDict $MotherFractureId BotLeftSurface Id [expr {[GiD_Info Geometry MaxNumSurfaces]-3}]
        set Lines [list [expr {[GiD_Info Geometry MaxNumLines]-1}] [expr {[GiD_Info Geometry MaxNumLines]-2}] [expr {[GiD_Info Geometry MaxNumLines]-4}]]
        dict set FracturesDict $MotherFractureId BotLeftSurface Lines $Lines
        dict set FracturesDict $MotherFractureId BotRightSurface Id [expr {[GiD_Info Geometry MaxNumSurfaces]-1}]
        set Lines [list [GiD_Info Geometry MaxNumLines] [expr {[GiD_Info Geometry MaxNumLines]-2}] [expr {[GiD_Info Geometry MaxNumLines]-7}]]
        dict set FracturesDict $MotherFractureId BotRightSurface Lines $Lines
        dict set FracturesDict $MotherFractureId LeftInterfaceVolume Id [expr {[GiD_Info Geometry MaxNumVolumes]-1}]
        dict set FracturesDict $MotherFractureId RightInterfaceVolume Id [GiD_Info Geometry MaxNumVolumes]
        # BodyVolumesDict
        dict set BodyVolumesDict $BodyVolumeId Surfaces $BodyVolumeSurfaces
    }

    # Loop through all Bifurcation Fractures
    dict for {BifId Bifurcation} $BifurcationDict {
        set MotherFractureId [dict get $Bifurcation MotherFractureId]

        set BodyVolumeId [lindex [dict get $FracturesDict $MotherFractureId BodyVolumes] 0]
        set BodyVolumeSurfaces [dict get $BodyVolumesDict $BodyVolumeId Surfaces]
        set OldTopLeftSurfaceLines [dict get $FracturesDict $MotherFractureId TopLeftSurface Lines]
        set OldTopRightSurfaceLines [dict get $FracturesDict $MotherFractureId TopRightSurface Lines]
        set OldBotLeftSurfaceLines [dict get $FracturesDict $MotherFractureId BotLeftSurface Lines]
        set OldBotRightSurfaceLines [dict get $FracturesDict $MotherFractureId BotRightSurface Lines]

        ## Modify Geometry
        # Delete InterfaceVolume
        GiD_Process Mescape Geometry Delete Volumes [dict get $FracturesDict $MotherFractureId LeftInterfaceVolume Id] [dict get $FracturesDict $MotherFractureId RightInterfaceVolume Id] escape
        # Delete old crack surfaces
        GiD_Process Mescape Geometry Delete Surfaces [dict get $FracturesDict $MotherFractureId TopLeftSurface Id] [dict get $FracturesDict $MotherFractureId TopRightSurface Id] \
            [dict get $FracturesDict $MotherFractureId BotLeftSurface Id] [dict get $FracturesDict $MotherFractureId BotRightSurface Id] escape
        set Index [lsearch $BodyVolumeSurfaces [dict get $FracturesDict $MotherFractureId TopLeftSurface Id]]
        set BodyVolumeSurfaces [lreplace $BodyVolumeSurfaces $Index $Index]
        set Index [lsearch $BodyVolumeSurfaces [dict get $FracturesDict $MotherFractureId TopRightSurface Id]]
        set BodyVolumeSurfaces [lreplace $BodyVolumeSurfaces $Index $Index]
        set Index [lsearch $BodyVolumeSurfaces [dict get $FracturesDict $MotherFractureId BotLeftSurface Id]]
        set BodyVolumeSurfaces [lreplace $BodyVolumeSurfaces $Index $Index]
        set Index [lsearch $BodyVolumeSurfaces [dict get $FracturesDict $MotherFractureId BotRightSurface Id]]
        set BodyVolumeSurfaces [lreplace $BodyVolumeSurfaces $Index $Index]
        # Delete old crack lines
        GiD_Process Mescape Geometry Delete Lines [dict get $FracturesDict $MotherFractureId TopLine Id] [dict get $FracturesDict $MotherFractureId BotLine Id] \
            [dict get $FracturesDict $MotherFractureId LeftLine Id] [dict get $FracturesDict $MotherFractureId RightLine Id] escape
        set Index [lsearch $OldTopLeftSurfaceLines [dict get $FracturesDict $MotherFractureId TopLine Id]]
        set OldTopLeftSurfaceLines [lreplace $OldTopLeftSurfaceLines $Index $Index]
        set Index [lsearch $OldTopLeftSurfaceLines [dict get $FracturesDict $MotherFractureId LeftLine Id]]
        set OldTopLeftSurfaceLines [lreplace $OldTopLeftSurfaceLines $Index $Index]
        set Index [lsearch $OldTopRightSurfaceLines [dict get $FracturesDict $MotherFractureId TopLine Id]]
        set OldTopRightSurfaceLines [lreplace $OldTopRightSurfaceLines $Index $Index]
        set Index [lsearch $OldTopRightSurfaceLines [dict get $FracturesDict $MotherFractureId RightLine Id]]
        set OldTopRightSurfaceLines [lreplace $OldTopRightSurfaceLines $Index $Index]
        set Index [lsearch $OldBotLeftSurfaceLines [dict get $FracturesDict $MotherFractureId BotLine Id]]
        set OldBotLeftSurfaceLines [lreplace $OldBotLeftSurfaceLines $Index $Index]
        set Index [lsearch $OldBotLeftSurfaceLines [dict get $FracturesDict $MotherFractureId LeftLine Id]]
        set OldBotLeftSurfaceLines [lreplace $OldBotLeftSurfaceLines $Index $Index]
        set Index [lsearch $OldBotRightSurfaceLines [dict get $FracturesDict $MotherFractureId BotLine Id]]
        set OldBotRightSurfaceLines [lreplace $OldBotRightSurfaceLines $Index $Index]
        set Index [lsearch $OldBotRightSurfaceLines [dict get $FracturesDict $MotherFractureId RightLine Id]]
        set OldBotRightSurfaceLines [lreplace $OldBotRightSurfaceLines $Index $Index]

        # Create new points
        GiD_Process Mescape Geometry Create Point
        # Create new point in TopInitCoordinates location
        GiD_Process [lindex [dict get $Bifurcation TopInitCoordinates] 0] [lindex [dict get $Bifurcation TopInitCoordinates] 1] [lindex [dict get $Bifurcation TopInitCoordinates] 2]
        # Create new point in BotInitCoordinates location
        GiD_Process [lindex [dict get $Bifurcation BotInitCoordinates] 0] [lindex [dict get $Bifurcation BotInitCoordinates] 1] [lindex [dict get $Bifurcation BotInitCoordinates] 2]
        # Create new point in LeftInitCoordinates location
        GiD_Process [lindex [dict get $Bifurcation LeftInitCoordinates] 0] [lindex [dict get $Bifurcation LeftInitCoordinates] 1] [lindex [dict get $Bifurcation LeftInitCoordinates] 2]
        # Create new point in RightInitCoordinates location
        GiD_Process [lindex [dict get $Bifurcation RightInitCoordinates] 0] [lindex [dict get $Bifurcation RightInitCoordinates] 1] [lindex [dict get $Bifurcation RightInitCoordinates] 2]
        # Create new point in TopTopEndCoordinates location
        GiD_Process [lindex [dict get $Bifurcation TopTopEndCoordinates] 0] [lindex [dict get $Bifurcation TopTopEndCoordinates] 1] [lindex [dict get $Bifurcation TopTopEndCoordinates] 2]
        # Create new point in TopBotEndCoordinates location
        GiD_Process [lindex [dict get $Bifurcation TopBotEndCoordinates] 0] [lindex [dict get $Bifurcation TopBotEndCoordinates] 1] [lindex [dict get $Bifurcation TopBotEndCoordinates] 2]
        # Create new point in TopLeftEndCoordinates location
        GiD_Process [lindex [dict get $Bifurcation TopLeftEndCoordinates] 0] [lindex [dict get $Bifurcation TopLeftEndCoordinates] 1] [lindex [dict get $Bifurcation TopLeftEndCoordinates] 2]
        # Create new point in TopRightEndCoordinates location
        GiD_Process [lindex [dict get $Bifurcation TopRightEndCoordinates] 0] [lindex [dict get $Bifurcation TopRightEndCoordinates] 1] [lindex [dict get $Bifurcation TopRightEndCoordinates] 2]
        # Create new point in TopTipCoordinates location
        GiD_Process [lindex [dict get $Bifurcation TopTipCoordinates] 0] [lindex [dict get $Bifurcation TopTipCoordinates] 1] [lindex [dict get $Bifurcation TopTipCoordinates] 2]
        # Create new point in BotTopEndCoordinates location
        GiD_Process [lindex [dict get $Bifurcation BotTopEndCoordinates] 0] [lindex [dict get $Bifurcation BotTopEndCoordinates] 1] [lindex [dict get $Bifurcation BotTopEndCoordinates] 2]
        # Create new point in BotBotEndCoordinates location
        GiD_Process [lindex [dict get $Bifurcation BotBotEndCoordinates] 0] [lindex [dict get $Bifurcation BotBotEndCoordinates] 1] [lindex [dict get $Bifurcation BotBotEndCoordinates] 2]
        # Create new point in BotLeftEndCoordinates location
        GiD_Process [lindex [dict get $Bifurcation BotLeftEndCoordinates] 0] [lindex [dict get $Bifurcation BotLeftEndCoordinates] 1] [lindex [dict get $Bifurcation BotLeftEndCoordinates] 2]
        # Create new point in BotRightEndCoordinates location
        GiD_Process [lindex [dict get $Bifurcation BotRightEndCoordinates] 0] [lindex [dict get $Bifurcation BotRightEndCoordinates] 1] [lindex [dict get $Bifurcation BotRightEndCoordinates] 2]
        # Create new point in BotTipCoordinates location
        GiD_Process [lindex [dict get $Bifurcation BotTipCoordinates] 0] [lindex [dict get $Bifurcation BotTipCoordinates] 1] [lindex [dict get $Bifurcation BotTipCoordinates] 2]
        GiD_Process escape
        
        # Create new lines
        GiD_Process Mescape Geometry Create Line Join [expr {[GiD_Info Geometry MaxNumPoints]-13}] [expr {[GiD_Info Geometry MaxNumPoints]-12}] escape \
            [dict get $FracturesDict $MotherFractureId TipPoint Id] [expr {[GiD_Info Geometry MaxNumPoints]-8}] escape \
            [dict get $FracturesDict $MotherFractureId TipPoint Id] [expr {[GiD_Info Geometry MaxNumPoints]-7}] escape \
            [dict get $FracturesDict $MotherFractureId TipPoint Id] [expr {[GiD_Info Geometry MaxNumPoints]-6}] escape \
            [dict get $FracturesDict $MotherFractureId TipPoint Id] [expr {[GiD_Info Geometry MaxNumPoints]-4}] escape \
            [dict get $FracturesDict $MotherFractureId TipPoint Id] [expr {[GiD_Info Geometry MaxNumPoints]-2}] escape \
            [dict get $FracturesDict $MotherFractureId TipPoint Id] [expr {[GiD_Info Geometry MaxNumPoints]-1}] escape \
            [expr {[GiD_Info Geometry MaxNumPoints]-11}] [dict get $FracturesDict $MotherFractureId TipPoint Id] escape \
            [expr {[GiD_Info Geometry MaxNumPoints]-11}] [expr {[GiD_Info Geometry MaxNumPoints]-7}] escape \
            [expr {[GiD_Info Geometry MaxNumPoints]-11}] [expr {[GiD_Info Geometry MaxNumPoints]-2}] escape \
            [expr {[GiD_Info Geometry MaxNumPoints]-10}] [dict get $FracturesDict $MotherFractureId TipPoint Id] escape \
            [expr {[GiD_Info Geometry MaxNumPoints]-10}] [expr {[GiD_Info Geometry MaxNumPoints]-6}] escape \
            [expr {[GiD_Info Geometry MaxNumPoints]-10}] [expr {[GiD_Info Geometry MaxNumPoints]-1}] escape \
            [expr {[GiD_Info Geometry MaxNumPoints]-13}] [expr {[GiD_Info Geometry MaxNumPoints]-6}] escape \
            [expr {[GiD_Info Geometry MaxNumPoints]-13}] [expr {[GiD_Info Geometry MaxNumPoints]-7}] escape \
            [expr {[GiD_Info Geometry MaxNumPoints]-13}] [dict get $FracturesDict $MotherFractureId TipPoint Id] escape \
            [expr {[GiD_Info Geometry MaxNumPoints]-13}] [expr {[GiD_Info Geometry MaxNumPoints]-9}] escape \
            [expr {[GiD_Info Geometry MaxNumPoints]-12}] [expr {[GiD_Info Geometry MaxNumPoints]-2}] escape \
            [expr {[GiD_Info Geometry MaxNumPoints]-12}] [expr {[GiD_Info Geometry MaxNumPoints]-1}] escape \
            [expr {[GiD_Info Geometry MaxNumPoints]-12}] [dict get $FracturesDict $MotherFractureId TipPoint Id] escape \
            [expr {[GiD_Info Geometry MaxNumPoints]-12}] [expr {[GiD_Info Geometry MaxNumPoints]-3}] escape \
            [expr {[GiD_Info Geometry MaxNumPoints]-11}] [expr {[GiD_Info Geometry MaxNumPoints]-13}] [expr {[GiD_Info Geometry MaxNumPoints]-10}] escape \
            [expr {[GiD_Info Geometry MaxNumPoints]-11}] [expr {[GiD_Info Geometry MaxNumPoints]-12}] [expr {[GiD_Info Geometry MaxNumPoints]-10}] escape \
            [expr {[GiD_Info Geometry MaxNumPoints]-7}] [expr {[GiD_Info Geometry MaxNumPoints]-9}] [expr {[GiD_Info Geometry MaxNumPoints]-6}] escape \
            [expr {[GiD_Info Geometry MaxNumPoints]-7}] [expr {[GiD_Info Geometry MaxNumPoints]-8}] [expr {[GiD_Info Geometry MaxNumPoints]-6}] escape \
            [expr {[GiD_Info Geometry MaxNumPoints]-2}] [expr {[GiD_Info Geometry MaxNumPoints]-4}] [expr {[GiD_Info Geometry MaxNumPoints]-1}] escape \
            [expr {[GiD_Info Geometry MaxNumPoints]-2}] [expr {[GiD_Info Geometry MaxNumPoints]-3}] [expr {[GiD_Info Geometry MaxNumPoints]-1}] escape \
            [dict get $FracturesDict $MotherFractureId RightPoint Id] [expr {[GiD_Info Geometry MaxNumPoints]-10}] escape \
            [dict get $FracturesDict $MotherFractureId LeftPoint Id] [expr {[GiD_Info Geometry MaxNumPoints]-11}] escape \
            [dict get $FracturesDict $MotherFractureId TopPoint Id] [expr {[GiD_Info Geometry MaxNumPoints]-13}] escape \
            [dict get $FracturesDict $MotherFractureId BotPoint Id] [expr {[GiD_Info Geometry MaxNumPoints]-12}] escape \
            [expr {[GiD_Info Geometry MaxNumPoints]-6}] [expr {[GiD_Info Geometry MaxNumPoints]-5}] [expr {[GiD_Info Geometry MaxNumPoints]-7}] escape \
            [expr {[GiD_Info Geometry MaxNumPoints]-9}] [expr {[GiD_Info Geometry MaxNumPoints]-5}] [expr {[GiD_Info Geometry MaxNumPoints]-8}] escape \
            [expr {[GiD_Info Geometry MaxNumPoints]-1}] [GiD_Info Geometry MaxNumPoints] [expr {[GiD_Info Geometry MaxNumPoints]-2}] escape \
            [expr {[GiD_Info Geometry MaxNumPoints]-4}] [GiD_Info Geometry MaxNumPoints] [expr {[GiD_Info Geometry MaxNumPoints]-3}] escape escape

        # Create new surfaces
        GiD_Process Mescape Geometry Create NurbsSurface [expr {[GiD_Info Geometry MaxNumLines]-22}] [expr {[GiD_Info Geometry MaxNumLines]-11}] [expr {[GiD_Info Geometry MaxNumLines]-9}]
        for {set i 0} {$i < [llength $OldTopRightSurfaceLines]} {incr i} {
            GiD_Process [lindex $OldTopRightSurfaceLines $i]
        }
        GiD_Process escape escape
        GiD_Process Mescape Geometry Create NurbsSurface [expr {[GiD_Info Geometry MaxNumLines]-23}] [expr {[GiD_Info Geometry MaxNumLines]-10}] [expr {[GiD_Info Geometry MaxNumLines]-9}]
        for {set i 0} {$i < [llength $OldTopLeftSurfaceLines]} {incr i} {
            GiD_Process [lindex $OldTopLeftSurfaceLines $i]
        }
        GiD_Process escape escape
        GiD_Process Mescape Geometry Create NurbsSurface [expr {[GiD_Info Geometry MaxNumLines]-20}] [expr {[GiD_Info Geometry MaxNumLines]-11}] [expr {[GiD_Info Geometry MaxNumLines]-8}]
        for {set i 0} {$i < [llength $OldBotRightSurfaceLines]} {incr i} {
            GiD_Process [lindex $OldBotRightSurfaceLines $i]
        }
        GiD_Process escape escape
        GiD_Process Mescape Geometry Create NurbsSurface [expr {[GiD_Info Geometry MaxNumLines]-21}] [expr {[GiD_Info Geometry MaxNumLines]-10}] [expr {[GiD_Info Geometry MaxNumLines]-8}]
        for {set i 0} {$i < [llength $OldBotLeftSurfaceLines]} {incr i} {
            GiD_Process [lindex $OldBotLeftSurfaceLines $i]
        }
        GiD_Process escape escape
        GiD_Process Mescape Geometry Create NurbsSurface [expr {[GiD_Info Geometry MaxNumLines]-20}] [expr {[GiD_Info Geometry MaxNumLines]-22}] [expr {[GiD_Info Geometry MaxNumLines]-44}] escape \
            [expr {[GiD_Info Geometry MaxNumLines]-21}] [expr {[GiD_Info Geometry MaxNumLines]-23}] [expr {[GiD_Info Geometry MaxNumLines]-44}] escape \
            [expr {[GiD_Info Geometry MaxNumLines]-44}] [expr {[GiD_Info Geometry MaxNumLines]-25}] [expr {[GiD_Info Geometry MaxNumLines]-29}] escape \
            [expr {[GiD_Info Geometry MaxNumLines]-20}] [expr {[GiD_Info Geometry MaxNumLines]-25}] [expr {[GiD_Info Geometry MaxNumLines]-34}] escape \
            [expr {[GiD_Info Geometry MaxNumLines]-21}] [expr {[GiD_Info Geometry MaxNumLines]-25}] [expr {[GiD_Info Geometry MaxNumLines]-37}] escape \
            [expr {[GiD_Info Geometry MaxNumLines]-22}] [expr {[GiD_Info Geometry MaxNumLines]-29}] [expr {[GiD_Info Geometry MaxNumLines]-34}] escape \
            [expr {[GiD_Info Geometry MaxNumLines]-23}] [expr {[GiD_Info Geometry MaxNumLines]-29}] [expr {[GiD_Info Geometry MaxNumLines]-37}] escape \
            [expr {[GiD_Info Geometry MaxNumLines]-21}] [expr {[GiD_Info Geometry MaxNumLines]-35}] [expr {[GiD_Info Geometry MaxNumLines]-27}] escape \
            [expr {[GiD_Info Geometry MaxNumLines]-20}] [expr {[GiD_Info Geometry MaxNumLines]-32}] [expr {[GiD_Info Geometry MaxNumLines]-26}] escape \
            [expr {[GiD_Info Geometry MaxNumLines]-32}] [expr {[GiD_Info Geometry MaxNumLines]-34}] [expr {[GiD_Info Geometry MaxNumLines]-38}] escape \
            [expr {[GiD_Info Geometry MaxNumLines]-37}] [expr {[GiD_Info Geometry MaxNumLines]-39}] [expr {[GiD_Info Geometry MaxNumLines]-35}] escape \
            [expr {[GiD_Info Geometry MaxNumLines]-26}] [expr {[GiD_Info Geometry MaxNumLines]-24}] [expr {[GiD_Info Geometry MaxNumLines]-12}] escape \
            [expr {[GiD_Info Geometry MaxNumLines]-27}] [expr {[GiD_Info Geometry MaxNumLines]-24}] [expr {[GiD_Info Geometry MaxNumLines]-13}] escape \
            [expr {[GiD_Info Geometry MaxNumLines]-38}] [expr {[GiD_Info Geometry MaxNumLines]-40}] [expr {[GiD_Info Geometry MaxNumLines]-14}] escape \
            [expr {[GiD_Info Geometry MaxNumLines]-39}] [expr {[GiD_Info Geometry MaxNumLines]-40}] [expr {[GiD_Info Geometry MaxNumLines]-15}] escape \
            [expr {[GiD_Info Geometry MaxNumLines]-36}] [expr {[GiD_Info Geometry MaxNumLines]-37}] [expr {[GiD_Info Geometry MaxNumLines]-42}] escape \
            [expr {[GiD_Info Geometry MaxNumLines]-33}] [expr {[GiD_Info Geometry MaxNumLines]-34}] [expr {[GiD_Info Geometry MaxNumLines]-41}] escape \
            [expr {[GiD_Info Geometry MaxNumLines]-41}] [expr {[GiD_Info Geometry MaxNumLines]-43}] [expr {[GiD_Info Geometry MaxNumLines]-16}] escape \
            [expr {[GiD_Info Geometry MaxNumLines]-42}] [expr {[GiD_Info Geometry MaxNumLines]-43}] [expr {[GiD_Info Geometry MaxNumLines]-17}] escape \
            [expr {[GiD_Info Geometry MaxNumLines]-30}] [expr {[GiD_Info Geometry MaxNumLines]-36}] [expr {[GiD_Info Geometry MaxNumLines]-23}] escape \
            [expr {[GiD_Info Geometry MaxNumLines]-31}] [expr {[GiD_Info Geometry MaxNumLines]-33}] [expr {[GiD_Info Geometry MaxNumLines]-22}] escape \
            [expr {[GiD_Info Geometry MaxNumLines]-28}] [expr {[GiD_Info Geometry MaxNumLines]-30}] [expr {[GiD_Info Geometry MaxNumLines]-19}] escape \
            [expr {[GiD_Info Geometry MaxNumLines]-28}] [expr {[GiD_Info Geometry MaxNumLines]-31}] [expr {[GiD_Info Geometry MaxNumLines]-18}] escape \
            [expr {[GiD_Info Geometry MaxNumLines]-12}] [expr {[GiD_Info Geometry MaxNumLines]-3}] [GiD_Info Geometry MaxNumLines] escape \
            [expr {[GiD_Info Geometry MaxNumLines]-13}] [expr {[GiD_Info Geometry MaxNumLines]-2}] [GiD_Info Geometry MaxNumLines] escape \
            [expr {[GiD_Info Geometry MaxNumLines]-14}] [expr {[GiD_Info Geometry MaxNumLines]-1}] [expr {[GiD_Info Geometry MaxNumLines]-3}] escape \
            [expr {[GiD_Info Geometry MaxNumLines]-15}] [expr {[GiD_Info Geometry MaxNumLines]-1}] [expr {[GiD_Info Geometry MaxNumLines]-2}] escape \
            [expr {[GiD_Info Geometry MaxNumLines]-16}] [expr {[GiD_Info Geometry MaxNumLines]-4}] [expr {[GiD_Info Geometry MaxNumLines]-7}] escape \
            [expr {[GiD_Info Geometry MaxNumLines]-17}] [expr {[GiD_Info Geometry MaxNumLines]-4}] [expr {[GiD_Info Geometry MaxNumLines]-6}] escape \
            [expr {[GiD_Info Geometry MaxNumLines]-18}] [expr {[GiD_Info Geometry MaxNumLines]-5}] [expr {[GiD_Info Geometry MaxNumLines]-7}] escape \
            [expr {[GiD_Info Geometry MaxNumLines]-19}] [expr {[GiD_Info Geometry MaxNumLines]-5}] [expr {[GiD_Info Geometry MaxNumLines]-6}] escape escape
        lappend BodyVolumeSurfaces [GiD_Info Geometry MaxNumSurfaces]
        lappend BodyVolumeSurfaces [expr {[GiD_Info Geometry MaxNumSurfaces]-1}]
        lappend BodyVolumeSurfaces [expr {[GiD_Info Geometry MaxNumSurfaces]-2}]
        lappend BodyVolumeSurfaces [expr {[GiD_Info Geometry MaxNumSurfaces]-3}]
        lappend BodyVolumeSurfaces [expr {[GiD_Info Geometry MaxNumSurfaces]-4}]
        lappend BodyVolumeSurfaces [expr {[GiD_Info Geometry MaxNumSurfaces]-5}]
        lappend BodyVolumeSurfaces [expr {[GiD_Info Geometry MaxNumSurfaces]-6}]
        lappend BodyVolumeSurfaces [expr {[GiD_Info Geometry MaxNumSurfaces]-7}]
        lappend BodyVolumeSurfaces [expr {[GiD_Info Geometry MaxNumSurfaces]-8}]
        lappend BodyVolumeSurfaces [expr {[GiD_Info Geometry MaxNumSurfaces]-9}]
        lappend BodyVolumeSurfaces [expr {[GiD_Info Geometry MaxNumSurfaces]-10}]
        lappend BodyVolumeSurfaces [expr {[GiD_Info Geometry MaxNumSurfaces]-11}]
        lappend BodyVolumeSurfaces [expr {[GiD_Info Geometry MaxNumSurfaces]-12}]
        lappend BodyVolumeSurfaces [expr {[GiD_Info Geometry MaxNumSurfaces]-13}]
        lappend BodyVolumeSurfaces [expr {[GiD_Info Geometry MaxNumSurfaces]-14}]
        lappend BodyVolumeSurfaces [expr {[GiD_Info Geometry MaxNumSurfaces]-15}]
        lappend BodyVolumeSurfaces [expr {[GiD_Info Geometry MaxNumSurfaces]-16}]
        lappend BodyVolumeSurfaces [expr {[GiD_Info Geometry MaxNumSurfaces]-17}]
        lappend BodyVolumeSurfaces [expr {[GiD_Info Geometry MaxNumSurfaces]-18}]
        lappend BodyVolumeSurfaces [expr {[GiD_Info Geometry MaxNumSurfaces]-19}]
        lappend BodyVolumeSurfaces [expr {[GiD_Info Geometry MaxNumSurfaces]-20}]
        lappend BodyVolumeSurfaces [expr {[GiD_Info Geometry MaxNumSurfaces]-21}]
        lappend BodyVolumeSurfaces [expr {[GiD_Info Geometry MaxNumSurfaces]-22}]
        lappend BodyVolumeSurfaces [expr {[GiD_Info Geometry MaxNumSurfaces]-23}]
        lappend BodyVolumeSurfaces [expr {[GiD_Info Geometry MaxNumSurfaces]-31}]
        lappend BodyVolumeSurfaces [expr {[GiD_Info Geometry MaxNumSurfaces]-32}]
        lappend BodyVolumeSurfaces [expr {[GiD_Info Geometry MaxNumSurfaces]-33}]
        lappend BodyVolumeSurfaces [expr {[GiD_Info Geometry MaxNumSurfaces]-34}]
        # Swap normals of new surfaces
        GiD_Process Mescape utilities SwapNormals Surfaces SelByNormal [GiD_Info Geometry MaxNumSurfaces] \
            escape [ComputeNormal [dict get $Bifurcation TopTopEndCoordinates] [dict get $Bifurcation TopTipCoordinates] [dict get $Bifurcation TopLeftEndCoordinates]] Yes escape
        GiD_Process Mescape utilities SwapNormals Surfaces SelByNormal [expr {[GiD_Info Geometry MaxNumSurfaces]-1}] \
            escape [ComputeNormal [dict get $Bifurcation TopTopEndCoordinates] [dict get $Bifurcation TopRightEndCoordinates] [dict get $Bifurcation TopTipCoordinates]] Yes escape
        GiD_Process Mescape utilities SwapNormals Surfaces SelByNormal [expr {[GiD_Info Geometry MaxNumSurfaces]-2}] \
            escape [ComputeNormal [dict get $Bifurcation TopBotEndCoordinates] [dict get $Bifurcation TopLeftEndCoordinates] [dict get $Bifurcation TopTipCoordinates]] Yes escape
        GiD_Process Mescape utilities SwapNormals Surfaces SelByNormal [expr {[GiD_Info Geometry MaxNumSurfaces]-3}] \
            escape [ComputeNormal [dict get $Bifurcation TopBotEndCoordinates] [dict get $Bifurcation TopTipCoordinates] [dict get $Bifurcation TopRightEndCoordinates]] Yes escape
        GiD_Process Mescape utilities SwapNormals Surfaces SelByNormal [expr {[GiD_Info Geometry MaxNumSurfaces]-4}] \
            escape [ComputeNormal [dict get $Bifurcation BotTopEndCoordinates] [dict get $Bifurcation BotTipCoordinates] [dict get $Bifurcation BotLeftEndCoordinates]] Yes escape
        GiD_Process Mescape utilities SwapNormals Surfaces SelByNormal [expr {[GiD_Info Geometry MaxNumSurfaces]-5}] \
            escape [ComputeNormal [dict get $Bifurcation BotTopEndCoordinates] [dict get $Bifurcation BotRightEndCoordinates] [dict get $Bifurcation BotTipCoordinates]] Yes escape
        GiD_Process Mescape utilities SwapNormals Surfaces SelByNormal [expr {[GiD_Info Geometry MaxNumSurfaces]-6}] \
            escape [ComputeNormal [dict get $Bifurcation BotBotEndCoordinates] [dict get $Bifurcation BotLeftEndCoordinates] [dict get $Bifurcation BotTipCoordinates]] Yes escape
        GiD_Process Mescape utilities SwapNormals Surfaces SelByNormal [expr {[GiD_Info Geometry MaxNumSurfaces]-7}] \
            escape [ComputeNormal [dict get $Bifurcation BotBotEndCoordinates] [dict get $Bifurcation BotTipCoordinates] [dict get $Bifurcation BotRightEndCoordinates]] Yes escape
        GiD_Process Mescape utilities SwapNormals Surfaces SelByNormal [expr {[GiD_Info Geometry MaxNumSurfaces]-31}] \
            escape [ComputeNormal [dict get $FracturesDict $MotherFractureId BotPoint Coordinates] [dict get $FracturesDict $MotherFractureId LeftPoint Coordinates] [dict get $Bifurcation BotInitCoordinates]] Yes escape
        GiD_Process Mescape utilities SwapNormals Surfaces SelByNormal [expr {[GiD_Info Geometry MaxNumSurfaces]-32}] \
            escape [ComputeNormal [dict get $FracturesDict $MotherFractureId BotPoint Coordinates] [dict get $Bifurcation BotInitCoordinates] [dict get $FracturesDict $MotherFractureId RightPoint Coordinates]] Yes escape
        GiD_Process Mescape utilities SwapNormals Surfaces SelByNormal [expr {[GiD_Info Geometry MaxNumSurfaces]-33}] \
            escape [ComputeNormal [dict get $FracturesDict $MotherFractureId TopPoint Coordinates] [dict get $Bifurcation TopInitCoordinates] [dict get $FracturesDict $MotherFractureId LeftPoint Coordinates]] Yes escape
        GiD_Process Mescape utilities SwapNormals Surfaces SelByNormal [expr {[GiD_Info Geometry MaxNumSurfaces]-34}] \
            escape [ComputeNormal [dict get $FracturesDict $MotherFractureId TopPoint Coordinates] [dict get $FracturesDict $MotherFractureId RightPoint Coordinates] [dict get $Bifurcation TopInitCoordinates]] Yes escape
        
        # Create new volumes
        # Old Contact volumes
        set Surf1 [list [expr {[GiD_Info Geometry MaxNumSurfaces]-31}] 1]
        set Surf2 [list [expr {[GiD_Info Geometry MaxNumSurfaces]-33}] 0]
        set TransformMatrix [ComputeTransformMatrix [dict get $FracturesDict $MotherFractureId LeftPoint Coordinates] [dict get $FracturesDict $MotherFractureId BotPoint Coordinates] \
            [dict get $FracturesDict $MotherFractureId TopPoint Coordinates] [dict get $Bifurcation LeftInitCoordinates] [dict get $FracturesDict $MotherFractureId LeftPoint Coordinates]]
        GiD_Geometry create volume [dict get $FracturesDict $MotherFractureId LeftInterfaceVolume Id] [dict get $FracturesDict $MotherFractureId LeftInterfaceVolume Layer] 2 \
            $Surf1 $Surf2 contactvolume $TransformMatrix
        set Surf1 [list [expr {[GiD_Info Geometry MaxNumSurfaces]-32}] 1]
        set Surf2 [list [expr {[GiD_Info Geometry MaxNumSurfaces]-34}] 0]
        set TransformMatrix [ComputeTransformMatrix [dict get $FracturesDict $MotherFractureId RightPoint Coordinates] [dict get $FracturesDict $MotherFractureId BotPoint Coordinates] \
            [dict get $FracturesDict $MotherFractureId TopPoint Coordinates] [dict get $FracturesDict $MotherFractureId RightPoint Coordinates] [dict get $Bifurcation RightInitCoordinates]]
        GiD_Geometry create volume [dict get $FracturesDict $MotherFractureId RightInterfaceVolume Id] [dict get $FracturesDict $MotherFractureId RightInterfaceVolume Layer] 2 \
            $Surf1 $Surf2 contactvolume $TransformMatrix
        # Link volumes
        GiD_Process Mescape Geometry Create volume [expr {[GiD_Info Geometry MaxNumSurfaces]-24}] [expr {[GiD_Info Geometry MaxNumSurfaces]-26}] [expr {[GiD_Info Geometry MaxNumSurfaces]-28}] [expr {[GiD_Info Geometry MaxNumSurfaces]-29}] escape \
            [expr {[GiD_Info Geometry MaxNumSurfaces]-25}] [expr {[GiD_Info Geometry MaxNumSurfaces]-27}] [expr {[GiD_Info Geometry MaxNumSurfaces]-28}] [expr {[GiD_Info Geometry MaxNumSurfaces]-30}] escape escape
        # New Contact volumes
        set Surf1 [list [expr {[GiD_Info Geometry MaxNumSurfaces]-6}] 1]
        set Surf2 [list [expr {[GiD_Info Geometry MaxNumSurfaces]-4}] 0]
        set TransformMatrix [ComputeTransformMatrix [dict get $Bifurcation BotLeftEndCoordinates] [dict get $Bifurcation BotBotEndCoordinates] \
            [dict get $Bifurcation BotTopEndCoordinates] [dict get $Bifurcation BotTipCoordinates] [dict get $Bifurcation BotLeftEndCoordinates]]
        GiD_Geometry create volume append [dict get $FracturesDict $MotherFractureId LeftInterfaceVolume Layer] 2 \
            $Surf1 $Surf2 contactvolume $TransformMatrix
        set Surf1 [list [expr {[GiD_Info Geometry MaxNumSurfaces]-7}] 1]
        set Surf2 [list [expr {[GiD_Info Geometry MaxNumSurfaces]-5}] 0]
        set TransformMatrix [ComputeTransformMatrix [dict get $Bifurcation BotRightEndCoordinates] [dict get $Bifurcation BotBotEndCoordinates] \
            [dict get $Bifurcation BotTopEndCoordinates] [dict get $Bifurcation BotRightEndCoordinates] [dict get $Bifurcation BotTipCoordinates]]
        GiD_Geometry create volume append [dict get $FracturesDict $MotherFractureId RightInterfaceVolume Layer] 2 \
            $Surf1 $Surf2 contactvolume $TransformMatrix
        set Surf1 [list [expr {[GiD_Info Geometry MaxNumSurfaces]-2}] 1]
        set Surf2 [list [GiD_Info Geometry MaxNumSurfaces] 0]
        set TransformMatrix [ComputeTransformMatrix [dict get $Bifurcation TopLeftEndCoordinates] [dict get $Bifurcation TopBotEndCoordinates] \
            [dict get $Bifurcation TopTopEndCoordinates] [dict get $Bifurcation TopTipCoordinates] [dict get $Bifurcation TopLeftEndCoordinates]]
        GiD_Geometry create volume append [dict get $FracturesDict $MotherFractureId LeftInterfaceVolume Layer] 2 \
            $Surf1 $Surf2 contactvolume $TransformMatrix
        set Surf1 [list [expr {[GiD_Info Geometry MaxNumSurfaces]-3}] 1]
        set Surf2 [list [expr {[GiD_Info Geometry MaxNumSurfaces]-1}] 0]
        set TransformMatrix [ComputeTransformMatrix [dict get $Bifurcation TopRightEndCoordinates] [dict get $Bifurcation TopBotEndCoordinates] \
            [dict get $Bifurcation TopTopEndCoordinates] [dict get $Bifurcation TopRightEndCoordinates] [dict get $Bifurcation TopTipCoordinates]]
        GiD_Geometry create volume append [dict get $FracturesDict $MotherFractureId RightInterfaceVolume Layer] 2 \
            $Surf1 $Surf2 contactvolume $TransformMatrix

        # Create PropagationUnion subgroups and assign points
        GiD_EntitiesGroups assign PropagationUnion_3d_6 points [expr {[GiD_Info Geometry MaxNumPoints]-1}]
        GiD_EntitiesGroups assign PropagationUnion_3d_6 points [expr {[GiD_Info Geometry MaxNumPoints]-2}]
        GiD_EntitiesGroups assign PropagationUnion_3d_6 points [expr {[GiD_Info Geometry MaxNumPoints]-3}]
        GiD_EntitiesGroups assign PropagationUnion_3d_6 points [expr {[GiD_Info Geometry MaxNumPoints]-4}]
        GiD_EntitiesGroups assign PropagationUnion_3d_6 points [expr {[GiD_Info Geometry MaxNumPoints]-6}]
        GiD_EntitiesGroups assign PropagationUnion_3d_6 points [expr {[GiD_Info Geometry MaxNumPoints]-7}]
        GiD_EntitiesGroups assign PropagationUnion_3d_6 points [expr {[GiD_Info Geometry MaxNumPoints]-8}]
        GiD_EntitiesGroups assign PropagationUnion_3d_6 points [expr {[GiD_Info Geometry MaxNumPoints]-9}]
        GiD_EntitiesGroups assign PropagationUnion_3d_6 points [expr {[GiD_Info Geometry MaxNumPoints]-10}]
        GiD_EntitiesGroups assign PropagationUnion_3d_6 points [expr {[GiD_Info Geometry MaxNumPoints]-11}]
        GiD_EntitiesGroups assign PropagationUnion_3d_6 points [expr {[GiD_Info Geometry MaxNumPoints]-12}]
        GiD_EntitiesGroups assign PropagationUnion_3d_6 points [expr {[GiD_Info Geometry MaxNumPoints]-13}]
        GiD_EntitiesGroups assign PropagationUnion_3d_6 points [dict get $FracturesDict $MotherFractureId TipPoint Id]
        # Bot Prism 1
        AddPropagationUnionPoint NumPropUnionGroups [expr {[GiD_Info Geometry MaxNumPoints]-12}]
        AddPropagationUnionPoint NumPropUnionGroups [expr {[GiD_Info Geometry MaxNumPoints]-10}]
        AddPropagationUnionPoint NumPropUnionGroups [expr {[GiD_Info Geometry MaxNumPoints]-1}]
        AddPropagationUnionPoint NumPropUnionGroups [dict get $FracturesDict $MotherFractureId TipPoint Id]
        AddPropagationUnionPoint NumPropUnionGroups [expr {[GiD_Info Geometry MaxNumPoints]-10}]
        AddPropagationUnionPoint NumPropUnionGroups [expr {[GiD_Info Geometry MaxNumPoints]-1}]
        # Bot Prism 2
        AddPropagationUnionPoint NumPropUnionGroups [expr {[GiD_Info Geometry MaxNumPoints]-12}]
        AddPropagationUnionPoint NumPropUnionGroups [expr {[GiD_Info Geometry MaxNumPoints]-1}]
        AddPropagationUnionPoint NumPropUnionGroups [expr {[GiD_Info Geometry MaxNumPoints]-3}]
        AddPropagationUnionPoint NumPropUnionGroups [dict get $FracturesDict $MotherFractureId TipPoint Id]
        AddPropagationUnionPoint NumPropUnionGroups [expr {[GiD_Info Geometry MaxNumPoints]-1}]
        AddPropagationUnionPoint NumPropUnionGroups [expr {[GiD_Info Geometry MaxNumPoints]-4}]
        # Bot Prism 3
        AddPropagationUnionPoint NumPropUnionGroups [expr {[GiD_Info Geometry MaxNumPoints]-12}]
        AddPropagationUnionPoint NumPropUnionGroups [expr {[GiD_Info Geometry MaxNumPoints]-2}]
        AddPropagationUnionPoint NumPropUnionGroups [expr {[GiD_Info Geometry MaxNumPoints]-11}]
        AddPropagationUnionPoint NumPropUnionGroups [dict get $FracturesDict $MotherFractureId TipPoint Id]
        AddPropagationUnionPoint NumPropUnionGroups [expr {[GiD_Info Geometry MaxNumPoints]-2}]
        AddPropagationUnionPoint NumPropUnionGroups [expr {[GiD_Info Geometry MaxNumPoints]-11}]
        # Bot Prism 4
        AddPropagationUnionPoint NumPropUnionGroups [expr {[GiD_Info Geometry MaxNumPoints]-12}]
        AddPropagationUnionPoint NumPropUnionGroups [expr {[GiD_Info Geometry MaxNumPoints]-3}]
        AddPropagationUnionPoint NumPropUnionGroups [expr {[GiD_Info Geometry MaxNumPoints]-2}]
        AddPropagationUnionPoint NumPropUnionGroups [dict get $FracturesDict $MotherFractureId TipPoint Id]
        AddPropagationUnionPoint NumPropUnionGroups [expr {[GiD_Info Geometry MaxNumPoints]-4}]
        AddPropagationUnionPoint NumPropUnionGroups [expr {[GiD_Info Geometry MaxNumPoints]-2}]
        # Top Prism 1
        AddPropagationUnionPoint NumPropUnionGroups [dict get $FracturesDict $MotherFractureId TipPoint Id]
        AddPropagationUnionPoint NumPropUnionGroups [expr {[GiD_Info Geometry MaxNumPoints]-10}]
        AddPropagationUnionPoint NumPropUnionGroups [expr {[GiD_Info Geometry MaxNumPoints]-6}]
        AddPropagationUnionPoint NumPropUnionGroups [expr {[GiD_Info Geometry MaxNumPoints]-13}]
        AddPropagationUnionPoint NumPropUnionGroups [expr {[GiD_Info Geometry MaxNumPoints]-10}]
        AddPropagationUnionPoint NumPropUnionGroups [expr {[GiD_Info Geometry MaxNumPoints]-6}]
        # Top Prism 2
        AddPropagationUnionPoint NumPropUnionGroups [dict get $FracturesDict $MotherFractureId TipPoint Id]
        AddPropagationUnionPoint NumPropUnionGroups [expr {[GiD_Info Geometry MaxNumPoints]-6}]
        AddPropagationUnionPoint NumPropUnionGroups [expr {[GiD_Info Geometry MaxNumPoints]-8}]
        AddPropagationUnionPoint NumPropUnionGroups [expr {[GiD_Info Geometry MaxNumPoints]-13}]
        AddPropagationUnionPoint NumPropUnionGroups [expr {[GiD_Info Geometry MaxNumPoints]-6}]
        AddPropagationUnionPoint NumPropUnionGroups [expr {[GiD_Info Geometry MaxNumPoints]-9}]
        # Top Prism 3
        AddPropagationUnionPoint NumPropUnionGroups [dict get $FracturesDict $MotherFractureId TipPoint Id]
        AddPropagationUnionPoint NumPropUnionGroups [expr {[GiD_Info Geometry MaxNumPoints]-7}]
        AddPropagationUnionPoint NumPropUnionGroups [expr {[GiD_Info Geometry MaxNumPoints]-11}]
        AddPropagationUnionPoint NumPropUnionGroups [expr {[GiD_Info Geometry MaxNumPoints]-13}]
        AddPropagationUnionPoint NumPropUnionGroups [expr {[GiD_Info Geometry MaxNumPoints]-7}]
        AddPropagationUnionPoint NumPropUnionGroups [expr {[GiD_Info Geometry MaxNumPoints]-11}]
        # Top Prism 4
        AddPropagationUnionPoint NumPropUnionGroups [dict get $FracturesDict $MotherFractureId TipPoint Id]
        AddPropagationUnionPoint NumPropUnionGroups [expr {[GiD_Info Geometry MaxNumPoints]-8}]
        AddPropagationUnionPoint NumPropUnionGroups [expr {[GiD_Info Geometry MaxNumPoints]-7}]
        AddPropagationUnionPoint NumPropUnionGroups [expr {[GiD_Info Geometry MaxNumPoints]-13}]
        AddPropagationUnionPoint NumPropUnionGroups [expr {[GiD_Info Geometry MaxNumPoints]-9}]
        AddPropagationUnionPoint NumPropUnionGroups [expr {[GiD_Info Geometry MaxNumPoints]-7}]


        ## Set Conditions
        for {set i 0} {$i < [llength [dict get $FracturesDict $MotherFractureId LeftInterfaceVolume Groups]]} {incr i} {
            GiD_EntitiesGroups assign [lindex [dict get $FracturesDict $MotherFractureId LeftInterfaceVolume Groups] $i] volumes [dict get $FracturesDict $MotherFractureId LeftInterfaceVolume Id]
            GiD_EntitiesGroups assign [lindex [dict get $FracturesDict $MotherFractureId LeftInterfaceVolume Groups] $i] volumes [expr {[GiD_Info Geometry MaxNumVolumes]-1}]
            GiD_EntitiesGroups assign [lindex [dict get $FracturesDict $MotherFractureId LeftInterfaceVolume Groups] $i] volumes [expr {[GiD_Info Geometry MaxNumVolumes]-3}]
        }
        for {set i 0} {$i < [llength [dict get $FracturesDict $MotherFractureId RightInterfaceVolume Groups]]} {incr i} {
            GiD_EntitiesGroups assign [lindex [dict get $FracturesDict $MotherFractureId RightInterfaceVolume Groups] $i] volumes [dict get $FracturesDict $MotherFractureId RightInterfaceVolume Id]
            GiD_EntitiesGroups assign [lindex [dict get $FracturesDict $MotherFractureId RightInterfaceVolume Groups] $i] volumes [GiD_Info Geometry MaxNumVolumes]
            GiD_EntitiesGroups assign [lindex [dict get $FracturesDict $MotherFractureId LeftInterfaceVolume Groups] $i] volumes [expr {[GiD_Info Geometry MaxNumVolumes]-2}]
        }
        GiD_EntitiesGroups assign $LinkInterfaceGroup volumes [expr {[GiD_Info Geometry MaxNumVolumes]-4}]
        GiD_EntitiesGroups assign $LinkInterfaceGroup volumes [expr {[GiD_Info Geometry MaxNumVolumes]-5}]


        ## Set Mesh options
        # GiD_Process Mescape Meshing Structured Lines 2 [expr {[GiD_Info Geometry MaxNumLines]-11}] \
        #     [expr {[GiD_Info Geometry MaxNumLines]-10}] [expr {[GiD_Info Geometry MaxNumLines]-9}] [expr {[GiD_Info Geometry MaxNumLines]-8}] \
        #     [expr {[GiD_Info Geometry MaxNumLines]-7}] [expr {[GiD_Info Geometry MaxNumLines]-6}] [expr {[GiD_Info Geometry MaxNumLines]-5}] \
        #     [expr {[GiD_Info Geometry MaxNumLines]-4}] [expr {[GiD_Info Geometry MaxNumLines]-3}] [expr {[GiD_Info Geometry MaxNumLines]-2}] \
        #     [expr {[GiD_Info Geometry MaxNumLines]-1}] [GiD_Info Geometry MaxNumLines] escape escape
        GiD_Process Mescape Meshing Structured Surfaces [expr {[GiD_Info Geometry MaxNumSurfaces]-30}] \
            [expr {[GiD_Info Geometry MaxNumSurfaces]-29}] [expr {[GiD_Info Geometry MaxNumSurfaces]-28}] [expr {[GiD_Info Geometry MaxNumSurfaces]-27}] \
            [expr {[GiD_Info Geometry MaxNumSurfaces]-26}] [expr {[GiD_Info Geometry MaxNumSurfaces]-25}] [expr {[GiD_Info Geometry MaxNumSurfaces]-24}] \
            [expr {[GiD_Info Geometry MaxNumSurfaces]-23}] [expr {[GiD_Info Geometry MaxNumSurfaces]-22}] [expr {[GiD_Info Geometry MaxNumSurfaces]-21}] \
            [expr {[GiD_Info Geometry MaxNumSurfaces]-20}] [expr {[GiD_Info Geometry MaxNumSurfaces]-19}] [expr {[GiD_Info Geometry MaxNumSurfaces]-18}] \
            [expr {[GiD_Info Geometry MaxNumSurfaces]-17}] [expr {[GiD_Info Geometry MaxNumSurfaces]-16}] [expr {[GiD_Info Geometry MaxNumSurfaces]-15}] \
            [expr {[GiD_Info Geometry MaxNumSurfaces]-14}] [expr {[GiD_Info Geometry MaxNumSurfaces]-13}] [expr {[GiD_Info Geometry MaxNumSurfaces]-12}] \
            [expr {[GiD_Info Geometry MaxNumSurfaces]-11}] [expr {[GiD_Info Geometry MaxNumSurfaces]-10}] [expr {[GiD_Info Geometry MaxNumSurfaces]-9}] \
            [expr {[GiD_Info Geometry MaxNumSurfaces]-8}] escape 1 [expr {[GiD_Info Geometry MaxNumLines]-44}] escape escape


        ## Update dictionaries
        # MotherFractureId is the new TopTip fracture
        dict set FracturesDict $MotherFractureId TipPoint Id [expr {[GiD_Info Geometry MaxNumPoints]-5}]
        dict set FracturesDict $MotherFractureId TipPoint Coordinates [dict get $Bifurcation TopTipCoordinates]
        dict set FracturesDict $MotherFractureId TopPoint Id [expr {[GiD_Info Geometry MaxNumPoints]-9}]
        dict set FracturesDict $MotherFractureId TopPoint Coordinates [dict get $Bifurcation TopTopEndCoordinates]
        dict set FracturesDict $MotherFractureId BotPoint Id [expr {[GiD_Info Geometry MaxNumPoints]-8}]
        dict set FracturesDict $MotherFractureId BotPoint Coordinates [dict get $Bifurcation TopBotEndCoordinates]
        dict set FracturesDict $MotherFractureId LeftPoint Id [expr {[GiD_Info Geometry MaxNumPoints]-7}]
        dict set FracturesDict $MotherFractureId LeftPoint Coordinates [dict get $Bifurcation TopLeftEndCoordinates]
        dict set FracturesDict $MotherFractureId RightPoint Id [expr {[GiD_Info Geometry MaxNumPoints]-6}]
        dict set FracturesDict $MotherFractureId RightPoint Coordinates [dict get $Bifurcation TopRightEndCoordinates]
        dict set FracturesDict $MotherFractureId TopLine Id [expr {[GiD_Info Geometry MaxNumLines]-5}]
        dict set FracturesDict $MotherFractureId BotLine Id [expr {[GiD_Info Geometry MaxNumLines]-4}]
        dict set FracturesDict $MotherFractureId LeftLine Id [expr {[GiD_Info Geometry MaxNumLines]-6}]
        dict set FracturesDict $MotherFractureId RightLine Id [expr {[GiD_Info Geometry MaxNumLines]-7}]
        dict set FracturesDict $MotherFractureId TopLeftSurface Id [GiD_Info Geometry MaxNumSurfaces]
        set Lines [list [expr {[GiD_Info Geometry MaxNumLines]-5}] [expr {[GiD_Info Geometry MaxNumLines]-6}] [expr {[GiD_Info Geometry MaxNumLines]-19}]]
        dict set FracturesDict $MotherFractureId TopLeftSurface Lines $Lines
        dict set FracturesDict $MotherFractureId TopRightSurface Id [expr {[GiD_Info Geometry MaxNumSurfaces]-1}]
        set Lines [list [expr {[GiD_Info Geometry MaxNumLines]-5}] [expr {[GiD_Info Geometry MaxNumLines]-7}] [expr {[GiD_Info Geometry MaxNumLines]-18}]]
        dict set FracturesDict $MotherFractureId TopRightSurface Lines $Lines
        dict set FracturesDict $MotherFractureId BotLeftSurface Id [expr {[GiD_Info Geometry MaxNumSurfaces]-2}]
        set Lines [list [expr {[GiD_Info Geometry MaxNumLines]-4}] [expr {[GiD_Info Geometry MaxNumLines]-6}] [expr {[GiD_Info Geometry MaxNumLines]-17}]]
        dict set FracturesDict $MotherFractureId BotLeftSurface Lines $Lines
        dict set FracturesDict $MotherFractureId BotRightSurface Id [expr {[GiD_Info Geometry MaxNumSurfaces]-3}]
        set Lines [list [expr {[GiD_Info Geometry MaxNumLines]-4}] [expr {[GiD_Info Geometry MaxNumLines]-7}] [expr {[GiD_Info Geometry MaxNumLines]-16}]]
        dict set FracturesDict $MotherFractureId BotRightSurface Lines $Lines
        dict set FracturesDict $MotherFractureId LeftInterfaceVolume Id [expr {[GiD_Info Geometry MaxNumVolumes]-1}]
        dict set FracturesDict $MotherFractureId RightInterfaceVolume Id [GiD_Info Geometry MaxNumVolumes]
        # New BotTip fracture
        set NewFractureId [dict size $FracturesDict]
        dict set FracturesDict $NewFractureId TipPoint Id [GiD_Info Geometry MaxNumPoints]
        dict set FracturesDict $NewFractureId TipPoint Coordinates [dict get $Bifurcation BotTipCoordinates]
        dict set FracturesDict $NewFractureId TopPoint Id [expr {[GiD_Info Geometry MaxNumPoints]-4}]
        dict set FracturesDict $NewFractureId TopPoint Coordinates [dict get $Bifurcation BotTopEndCoordinates]
        dict set FracturesDict $NewFractureId BotPoint Id [expr {[GiD_Info Geometry MaxNumPoints]-3}]
        dict set FracturesDict $NewFractureId BotPoint Coordinates [dict get $Bifurcation BotBotEndCoordinates]
        dict set FracturesDict $NewFractureId LeftPoint Id [expr {[GiD_Info Geometry MaxNumPoints]-2}]
        dict set FracturesDict $NewFractureId LeftPoint Coordinates [dict get $Bifurcation BotLeftEndCoordinates]
        dict set FracturesDict $NewFractureId RightPoint Id [expr {[GiD_Info Geometry MaxNumPoints]-1}]
        dict set FracturesDict $NewFractureId RightPoint Coordinates [dict get $Bifurcation BotRightEndCoordinates]
        dict set FracturesDict $NewFractureId TopLine Id [expr {[GiD_Info Geometry MaxNumLines]-1}]
        dict set FracturesDict $NewFractureId BotLine Id [GiD_Info Geometry MaxNumLines]
        dict set FracturesDict $NewFractureId LeftLine Id [expr {[GiD_Info Geometry MaxNumLines]-2}]
        dict set FracturesDict $NewFractureId RightLine Id [expr {[GiD_Info Geometry MaxNumLines]-3}]
        dict set FracturesDict $NewFractureId TopLeftSurface Id [expr {[GiD_Info Geometry MaxNumSurfaces]-4}]
        set Lines [list [expr {[GiD_Info Geometry MaxNumLines]-1}] [expr {[GiD_Info Geometry MaxNumLines]-2}] [expr {[GiD_Info Geometry MaxNumLines]-15}]]
        dict set FracturesDict $NewFractureId TopLeftSurface Lines $Lines
        dict set FracturesDict $NewFractureId TopRightSurface Id [expr {[GiD_Info Geometry MaxNumSurfaces]-5}]
        set Lines [list [expr {[GiD_Info Geometry MaxNumLines]-1}] [expr {[GiD_Info Geometry MaxNumLines]-3}] [expr {[GiD_Info Geometry MaxNumLines]-14}]]
        dict set FracturesDict $NewFractureId TopRightSurface Lines $Lines
        dict set FracturesDict $NewFractureId BotLeftSurface Id [expr {[GiD_Info Geometry MaxNumSurfaces]-6}]
        set Lines [list [expr {[GiD_Info Geometry MaxNumLines]-2}] [GiD_Info Geometry MaxNumLines] [expr {[GiD_Info Geometry MaxNumLines]-13}]]
        dict set FracturesDict $NewFractureId BotLeftSurface Lines $Lines
        dict set FracturesDict $NewFractureId BotRightSurface Id [expr {[GiD_Info Geometry MaxNumSurfaces]-7}]
        set Lines [list [expr {[GiD_Info Geometry MaxNumLines]-3}] [GiD_Info Geometry MaxNumLines] [expr {[GiD_Info Geometry MaxNumLines]-12}]]
        dict set FracturesDict $NewFractureId BotRightSurface Lines $Lines
        dict set FracturesDict $NewFractureId LeftInterfaceVolume Id [expr {[GiD_Info Geometry MaxNumVolumes]-3}]
        dict set FracturesDict $NewFractureId LeftInterfaceVolume Layer [dict get $FracturesDict $MotherFractureId LeftInterfaceVolume Layer]
        dict set FracturesDict $NewFractureId LeftInterfaceVolume Groups [dict get $FracturesDict $MotherFractureId LeftInterfaceVolume Groups]
        dict set FracturesDict $NewFractureId RightInterfaceVolume Id [expr {[GiD_Info Geometry MaxNumVolumes]-2}]
        dict set FracturesDict $NewFractureId RightInterfaceVolume Layer [dict get $FracturesDict $MotherFractureId RightInterfaceVolume Layer]
        dict set FracturesDict $NewFractureId RightInterfaceVolume Groups [dict get $FracturesDict $MotherFractureId RightInterfaceVolume Groups]
        # BodyVolumesDict
        dict set BodyVolumesDict $BodyVolumeId Surfaces $BodyVolumeSurfaces
    }

    # Create BodyVolumes again assign conditions and replace BodyVolumesDict and FracturesDict with the new values
    set BodyVolumesAuxDict $BodyVolumesDict
    unset BodyVolumesDict
    dict for {BodyId BodyVolume} $BodyVolumesAuxDict {
        set BodyVolumeSurfaces [dict get $BodyVolume Surfaces]
        GiD_Process Mescape Geometry Create volume
        for {set i 0} {$i < [llength $BodyVolumeSurfaces]} {incr i} {
            GiD_Process [lindex $BodyVolumeSurfaces $i]
        }
        GiD_Process escape escape
        
        set NewBodyVolumeId [GiD_Info Geometry MaxNumVolumes]
        
        set BodyVolumeGroups [dict get $BodyVolume Groups]
        for {set i 0} {$i < [llength $BodyVolumeGroups]} {incr i} {
            GiD_EntitiesGroups assign [lindex $BodyVolumeGroups $i] volumes $NewBodyVolumeId
        }

        if {[dict get $BodyVolume MeshSize] > 0.0} {
            GiD_Process Mescape Meshing AssignSizes Volumes [dict get $BodyVolume MeshSize] $NewBodyVolumeId escape escape
        }

        dict set BodyVolumesDict $NewBodyVolumeId Groups [dict get $BodyVolume Groups]
        dict set BodyVolumesDict $NewBodyVolumeId Surfaces [dict get $BodyVolume Surfaces]
        dict set BodyVolumesDict $NewBodyVolumeId MeshSize [dict get $BodyVolume MeshSize]
    }

    dict for {Id Fracture} $FracturesDict {
        set BodyVolumes [list]
        dict for {BodyId BodyVolume} $BodyVolumesDict {
            if {[GiD_Info IsPointInside Volume $BodyId [dict get $Fracture TipPoint Coordinates]] eq 1} {
                lappend BodyVolumes $BodyId
            }
        }
        dict set FracturesDict $Id BodyVolumes $BodyVolumes
    }

    # Generate New Mesh
    GiD_Process Mescape Meshing Generate Yes [GiD_Info Project LastElementSize] MeshingParametersFrom=Model Mescape

    ## Update FracturesData.json file
    set filename [file join $dir FracturesData.json]
    set FileVar [open $filename w]
        
    puts $FileVar "\{"

    ## fracture_data
    puts $FileVar "    \"fracture_data\": \{"
    puts $FileVar "        \"gid_path\":                             \"[lindex $PropagationData 0]\","
    # propagation parameters and body_domain_sub_model_part_list
    WriteFractureData FileVar
    # interface_domain_sub_model_part_list
    set PutStrings \[
    set Groups [GiD_Info conditions Interface_Part groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
        append PutStrings \" [lindex [lindex $Groups $i] 1] \" ,
    }
    append PutStrings \" PropagationUnion_3d_6 \" \]
    puts $FileVar "        \"interface_domain_sub_model_part_list\": $PutStrings,"
    # interface_domain_sub_model_part_old_list
    puts $FileVar "        \"interface_domain_sub_model_part_old_list\": $InterfaceGroupsOld"
    puts $FileVar "    \},"

    ## body_volumes_list
    WriteBodyVolumesList FileVar $BodyVolumesDict

    ## fractures_list
    WriteFracturesList FileVar $FracturesDict

    puts $FileVar "\}"

    close $FileVar
}

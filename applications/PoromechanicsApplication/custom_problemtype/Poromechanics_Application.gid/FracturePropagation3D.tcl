proc WriteInitialFracturesData { dir problemtypedir gidpath } {
    
    ## Source auxiliar procedures
    source [file join $problemtypedir Fracture3DAuxProcs.tcl]

    ## Set BodyVolumesDict
    set BodyGroups [GiD_Info conditions Body_Part groups]
    set BodyVolumesDict [dict create]
    for {set i 0} {$i < [llength $BodyGroups]} {incr i} {
        set BodyEntities [GiD_EntitiesGroups get [lindex [lindex $BodyGroups $i] 1] volumes]
        for {set j 0} {$j < [llength $BodyEntities]} {incr j} {
            set BodyVolume [GiD_Geometry get volume [lindex $BodyEntities $j]]
            set Groups [GiD_EntitiesGroups entity_groups volumes [lindex $BodyEntities $j]]
            dict set BodyVolumesDict [lindex $BodyEntities $j] Layer [lindex $BodyVolume 0]
            dict set BodyVolumesDict [lindex $BodyEntities $j] Groups $Groups
            set Surfaces [list]
            for {set k 0} {$k < [lindex $BodyVolume 1]} {incr k} {
                lappend Surfaces [lindex [lindex $BodyVolume [expr {2+$k}]] 0]
            }
            dict set BodyVolumesDict [lindex $BodyEntities $j] Surfaces $Surfaces
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
        # Notes: eq-> for string comparison (and int), ==-> for numerical comparison
        #        ne-> for string comparison (and int), !=-> for numerical comparison
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
                            dict set PointsDict [lindex $AxisLine 2] ${FractureHalf}Volume [lindex $InterfaceEntities $j]
                            dict set PointsDict [lindex $AxisLine 2] ${FractureHalf}Line $Line1Id
                            dict set PointsDict [lindex $AxisLine 3] ${FractureHalf}Volume [lindex $InterfaceEntities $j]
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
        if {[dict exists $Info RightVolume] && [dict exists $Info LeftVolume]} {
            # Define new fracture
            incr FractureId
            # Set TipPoint
            set TipPoint [GiD_Geometry get point $PointId]
            AddPointToFracturesDict FracturesDict $FractureId $PointId $TipPoint TipPoint
            # Set BodyVolume
            AddBodyVolumeToFracturesDict FracturesDict $FractureId $BodyVolumesDict
            # Set LeftVolume and RightVolume
            set LeftVolume [GiD_Geometry get volume [dict get $Info LeftVolume]]
            AddInterfaceVolumeToFracturesDict FracturesDict $FractureId [dict get $Info LeftVolume] $LeftVolume LeftVolume
            set RightVolume [GiD_Geometry get volume [dict get $Info RightVolume]]
            AddInterfaceVolumeToFracturesDict FracturesDict $FractureId [dict get $Info RightVolume] $RightVolume RightVolume
            # Set LeftLine and RightLine
            set LeftLine [GiD_Geometry get line [dict get $Info LeftLine]]
            AddLineToFracturesDict FracturesDict $FractureId [dict get $Info LeftLine] $LeftLine LeftLine
            set RightLine [GiD_Geometry get line [dict get $Info RightLine]]
            AddLineToFracturesDict FracturesDict $FractureId [dict get $Info RightLine] $RightLine RightLine
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
            set LeftSurface1 [GiD_Geometry get surface [lindex [lindex $LeftVolume 2] 0]]
            set LeftSurface2 [GiD_Geometry get surface [lindex [lindex $LeftVolume 3] 0]]
            set RightSurface1 [GiD_Geometry get surface [lindex [lindex $RightVolume 2] 0]]
            set RightSurface2 [GiD_Geometry get surface [lindex [lindex $RightVolume 3] 0]]
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
                AddSurfaceToFracturesDict FracturesDict $FractureId [lindex [lindex $LeftVolume 2] 0] $LeftSurface1 TopLeftSurface
                AddLineToFracturesDict FracturesDict $FractureId $AxisLineId $AxisLine TopLine
                AddPointToFracturesDict FracturesDict $FractureId $TopBotPointId $TopBotPoint TopPoint
                if {$IsRightSurface1 eq 1} {                
                    AddSurfaceToFracturesDict FracturesDict $FractureId [lindex [lindex $RightVolume 2] 0] $RightSurface1 TopRightSurface
                    SearchAxis AxisLineId AxisLine TopBotPointId TopBotPoint $LeftSurface2 $RightSurface2 $PointId
                    AddSurfaceToFracturesDict FracturesDict $FractureId [lindex [lindex $RightVolume 3] 0] $RightSurface2 BotRightSurface
                } else {
                    AddSurfaceToFracturesDict FracturesDict $FractureId [lindex [lindex $RightVolume 3] 0] $RightSurface2 TopRightSurface
                    SearchAxis AxisLineId AxisLine TopBotPointId TopBotPoint $LeftSurface2 $RightSurface1 $PointId
                    AddSurfaceToFracturesDict FracturesDict $FractureId [lindex [lindex $RightVolume 2] 0] $RightSurface1 BotRightSurface
                }
                AddSurfaceToFracturesDict FracturesDict $FractureId [lindex [lindex $LeftVolume 3] 0] $LeftSurface2 BotLeftSurface
                AddLineToFracturesDict FracturesDict $FractureId $AxisLineId $AxisLine BotLine
                AddPointToFracturesDict FracturesDict $FractureId $TopBotPointId $TopBotPoint BotPoint
            } else {
                AddSurfaceToFracturesDict FracturesDict $FractureId [lindex [lindex $LeftVolume 2] 0] $LeftSurface1 BotLeftSurface
                AddLineToFracturesDict FracturesDict $FractureId $AxisLineId $AxisLine BotLine
                AddPointToFracturesDict FracturesDict $FractureId $TopBotPointId $TopBotPoint BotPoint
                if {$IsRightSurface1 eq 1} {
                    AddSurfaceToFracturesDict FracturesDict $FractureId [lindex [lindex $RightVolume 2] 0] $RightSurface1 BotRightSurface
                    SearchAxis AxisLineId AxisLine TopBotPointId TopBotPoint $LeftSurface2 $RightSurface2 $PointId
                    AddSurfaceToFracturesDict FracturesDict $FractureId [lindex [lindex $RightVolume 3] 0] $RightSurface2 TopRightSurface
                } else {
                    AddSurfaceToFracturesDict FracturesDict $FractureId [lindex [lindex $RightVolume 3] 0] $RightSurface2 BotRightSurface
                    SearchAxis AxisLineId AxisLine TopBotPointId TopBotPoint $LeftSurface2 $RightSurface1 $PointId
                    AddSurfaceToFracturesDict FracturesDict $FractureId [lindex [lindex $RightVolume 2] 0] $RightSurface1 TopRightSurface
                }
                AddSurfaceToFracturesDict FracturesDict $FractureId [lindex [lindex $LeftVolume 3] 0] $LeftSurface2 TopLeftSurface
                AddLineToFracturesDict FracturesDict $FractureId $AxisLineId $AxisLine TopLine
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
    puts $FileVar1 "        \"gid_path\":                                 \"$gidpath\","
    # propagation parameters and body_domain_sub_model_part_list
    WriteFractureData FileVar1
    # interface_domain_sub_model_part_list
    set PutStrings \[
    set Groups [GiD_Info conditions Interface_Part groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
        append PutStrings \" [lindex [lindex $Groups $i] 1] \" ,
    }
    set PutStrings [string trimright $PutStrings ,]
    append PutStrings \]
    puts $FileVar1 "        \"interface_domain_sub_model_part_list\":     $PutStrings"
    puts $FileVar1 "    \},"
    
    ## body_volumes_list
    WriteBodyVolumesList FileVar1 $BodyVolumesDict
        
    ## fractures_list
    WriteFracturesList FileVar1 $FracturesDict
    
    puts $FileVar1 "\}"
    
    close $FileVar1

    ## Start gid_preferences.ini file
    set filename [file join $dir gid_preferences.ini]
    set FileVar2 [open $filename w]
        
    puts $FileVar2 "GID_OMP_NUM_THREADS [GiD_Set GID_OMP_NUM_THREADS]"
    puts $FileVar2 "AutomaticCorrectSizes [GiD_Set -meshing_parameters_model AutomaticCorrectSizes]"
    puts $FileVar2 "SplashWindow 0"
    puts $FileVar2 "SizeTransitionsFactor [GiD_Set -meshing_parameters_model SizeTransitionsFactor]"
    puts $FileVar2 "BoundaryWeightedTransition [GiD_Set -meshing_parameters_model BoundaryWeightedTransition]"
    puts $FileVar2 "SurfaceMesher [GiD_Set -meshing_parameters_model SurfaceMesher]"
    puts $FileVar2 "VolumeMesher [GiD_Set -meshing_parameters_model VolumeMesher]"
    
    close $FileVar2
}

proc GenerateNewFractures { dir problemtypedir PropagationData } {

    ## Source auxiliar procedures
    source [file join $problemtypedir Fracture3DAuxProcs.tcl]
    
    ###TODO
    ## Grups fills d'un grup especial (Prova) amb els nodes dels elements de colze
    #~ set Groups [GiD_Groups list Prova]
    #~ WarnWin $Groups
    #~ set NomGroup1 "Prova//[lindex $Groups 0]"
    #~ set NomGroup2 "Prova//[lindex $Groups 1]"
    #~ set Entities [GiD_EntitiesGroups get $NomGroup1 nodes]
    #~ WarnWin "Node Group1: $Entities"
    #~ set Entities [GiD_EntitiesGroups get $NomGroup2 nodes]
    #~ WarnWin "Node Group2: $Entities"
    ## Auxiliar Condition to save groups in it. Ja no cal
    #~ GiD_CreateData create condition AuxCond ovgroup ovnode {{Zero} {0.0}}
    ###TODO
}
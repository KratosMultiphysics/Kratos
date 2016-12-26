proc AddPointToFracturesDict {FracturesDict FractureId PointId Point PointType} {
    upvar $FracturesDict MyFracturesDict
    
    dict set MyFracturesDict $FractureId $PointType Id $PointId
    dict set MyFracturesDict $FractureId $PointType Layer [lindex $Point 0]
    dict set MyFracturesDict $FractureId $PointType Groups [GiD_EntitiesGroups entity_groups points $PointId]
    dict set MyFracturesDict $FractureId $PointType Coordinates "[lindex $Point 1] [lindex $Point 2] [lindex $Point 3]"
}

#-------------------------------------------------------------------------------

proc AddLineToFracturesDict {FracturesDict FractureId LineId Line LineType} {
    upvar $FracturesDict MyFracturesDict
    
    dict set MyFracturesDict $FractureId $LineType Id $LineId
    dict set MyFracturesDict $FractureId $LineType Layer [lindex $Line 1]
    dict set MyFracturesDict $FractureId $LineType Groups [GiD_EntitiesGroups entity_groups lines $LineId]
}

#-------------------------------------------------------------------------------

proc AddSurfaceToFracturesDict {FracturesDict FractureId SurfaceId Surface SurfaceType} {
    upvar $FracturesDict MyFracturesDict
    
    dict set MyFracturesDict $FractureId $SurfaceType Id $SurfaceId
    dict set MyFracturesDict $FractureId $SurfaceType Layer [lindex $Surface 1]
    dict set MyFracturesDict $FractureId $SurfaceType Groups [GiD_EntitiesGroups entity_groups surfaces $SurfaceId]
}

#-------------------------------------------------------------------------------

proc AddInterfaceVolumeToFracturesDict {FracturesDict FractureId VolumeId Volume VolumeType} {
    upvar $FracturesDict MyFracturesDict
    
    dict set MyFracturesDict $FractureId $VolumeType Id $VolumeId
    dict set MyFracturesDict $FractureId $VolumeType Layer [lindex $Volume 0]
    dict set MyFracturesDict $FractureId $VolumeType Groups [GiD_EntitiesGroups entity_groups volumes $VolumeId]
}

#-------------------------------------------------------------------------------

proc AddBodyVolumeToFracturesDict {FracturesDict FractureId BodyVolumesDict} {
    upvar $FracturesDict MyFracturesDict

    set TipCoordinates [dict get $MyFracturesDict $FractureId TipPoint Coordinates]
    set BodyVolumes [list]
    dict for {Id BodyVolume} $BodyVolumesDict {
        if {[GiD_Info IsPointInside Volume $Id $TipCoordinates] == 1} {
            lappend BodyVolumes $Id
        }
    }
    dict set MyFracturesDict $FractureId BodyVolumes $BodyVolumes
}

#-------------------------------------------------------------------------------

proc IsTopPoint {Point LeftPoint RightPoint TipPoint} {
        
    # LeftPoint Coordinates
    set Lp(0) [lindex $LeftPoint 1]
    set Lp(1) [lindex $LeftPoint 2]
    set Lp(2) [lindex $LeftPoint 3]
    
    # RightPoint Coordinates
    set Rp(0) [lindex $RightPoint 1]
    set Rp(1) [lindex $RightPoint 2]
    set Rp(2) [lindex $RightPoint 3]

    # TipPoint Coordinates
    set Tp(0) [lindex $TipPoint 1]
    set Tp(1) [lindex $TipPoint 2]
    set Tp(2) [lindex $TipPoint 3]
    
    # Unitary vector in local x direction
    set Vx(0) [expr {$Tp(0)-$Rp(0)}]
    set Vx(1) [expr {$Tp(1)-$Rp(1)}]
    set Vx(2) [expr {$Tp(2)-$Rp(2)}]
    set InvNorm [expr {1.0/sqrt($Vx(0)*$Vx(0)+$Vx(1)*$Vx(1)+$Vx(2)*$Vx(2))}]
    set Vx(0) [expr {$Vx(0)*$InvNorm}]
    set Vx(1) [expr {$Vx(1)*$InvNorm}]
    set Vx(2) [expr {$Vx(2)*$InvNorm}]

    # Unitary vector in local y direction
    set Vy(0) [expr {$Lp(0)-$Rp(0)}]
    set Vy(1) [expr {$Lp(1)-$Rp(1)}]
    set Vy(2) [expr {$Lp(2)-$Rp(2)}]
    set InvNorm [expr {1.0/sqrt($Vy(0)*$Vy(0)+$Vy(1)*$Vy(1)+$Vy(2)*$Vy(2))}]
    set Vy(0) [expr {$Vy(0)*$InvNorm}]
    set Vy(1) [expr {$Vy(1)*$InvNorm}]
    set Vy(2) [expr {$Vy(2)*$InvNorm}]
    
    # Unitary vector in local z direction (Cross product between Vx and Vy)
    set Vz(0) [expr {$Vx(1)*$Vy(2)-$Vx(2)*$Vy(1)}]
    set Vz(1) [expr {$Vx(2)*$Vy(0)-$Vx(0)*$Vy(2)}]
    set Vz(2) [expr {$Vx(0)*$Vy(1)-$Vx(1)*$Vy(0)}]
    
    # Rotation matrix (only z component)
    set R(20) $Vz(0)
    set R(21) $Vz(1)
    set R(22) $Vz(2)
    
    # Local z component of the TipPoint
    set Tz [expr {$R(20)*$Tp(0)+$R(21)*$Tp(1)+$R(22)*$Tp(2)}]
    
    # Global coordinates of the point
    set P(0) [lindex $Point 1]
    set P(1) [lindex $Point 2]
    set P(2) [lindex $Point 3]
    
    # Local z component of the Point
    set Pz [expr {$R(20)*$P(0)+$R(21)*$P(1)+$R(22)*$P(2)}]

    if {$Pz > $Tz} {
        return 1
    } else {
        return 0
    }
}

#-------------------------------------------------------------------------------

proc SearchAxis {AxisLineId AxisLine TopBotPointId TopBotPoint LeftSurface RightSurface TipPointId} {
    upvar $AxisLineId MyAxisLineId
    upvar $AxisLine MyAxisLine
    upvar $TopBotPointId MyTopBotPointId
    upvar $TopBotPoint MyTopBotPoint
    
    for {set i 0} {$i < [lindex $LeftSurface 2]} {incr i} {
        set LeftLineId [lindex [lindex $LeftSurface [expr { 9+$i }]] 0]
        for {set j 0} {$j < [lindex $RightSurface 2]} {incr j} {
            set RightLineId [lindex [lindex $RightSurface [expr { 9+$j }]] 0]
            if {$LeftLineId eq $RightLineId} {
                set MyAxisLineId $LeftLineId
                set MyAxisLine [GiD_Geometry get line $MyAxisLineId]
                if {[lindex $MyAxisLine 2] ne $TipPointId} {
                    set MyTopBotPointId [lindex $MyAxisLine 2]
                    set MyTopBotPoint [GiD_Geometry get point $MyTopBotPointId]
                } else {
                    set MyTopBotPointId [lindex $MyAxisLine 3]
                    set MyTopBotPoint [GiD_Geometry get point $MyTopBotPointId]
                }
            }
        }
    }
}

#--------------------------------------------------------------------------------------------------------------------------------------------------------------

proc WriteFractureData {FileVar} {
    upvar $FileVar MyFileVar

    puts $MyFileVar "        \"propagation_damage\":                       [GiD_AccessValue get gendata Propagation_Damage],"
    puts $MyFileVar "        \"propagation_length\":                       [GiD_AccessValue get gendata Propagation_Length],"
    puts $MyFileVar "        \"propagation_width\":                        [GiD_AccessValue get gendata Propagation_Width],"
    puts $MyFileVar "        \"propagation_height\":                       [GiD_AccessValue get gendata Propagation_Height],"
    puts $MyFileVar "        \"propagation_cosangle\":                     [GiD_AccessValue get gendata Propagation_CosAngle],"
    puts $MyFileVar "        \"propagation_frequency\":                    [GiD_AccessValue get gendata Propagation_Frequency],"
    # body_domain_sub_model_part_list
    set PutStrings \[
    set Groups [GiD_Info conditions Body_Part groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
        append PutStrings \" [lindex [lindex $Groups $i] 1] \" ,
    }
    set PutStrings [string trimright $PutStrings ,]
    append PutStrings \]
    puts $MyFileVar "        \"body_domain_sub_model_part_list\":          $PutStrings,"
}

#-------------------------------------------------------------------------------

proc WriteBodyVolumesList {FileVar BodyVolumesDict} {
    upvar $FileVar MyFileVar

    set iter 0
    puts $MyFileVar "    \"body_volumes_list\": \[\{"
    dict for {Id BodyVolume} $BodyVolumesDict {
        incr iter
        puts $MyFileVar "        \"id\":     $Id,"
        puts $MyFileVar "        \"layer\":  \"[dict get $BodyVolume Layer]\","
        if {[llength [dict get $BodyVolume Groups]]==0} {
            puts $MyFileVar "        \"groups\": \[\],"
        } else {
            set PutStrings \[
            for {set i 0} {$i < [llength [dict get $BodyVolume Groups]]} {incr i} {
                append PutStrings \" [lindex [dict get $BodyVolume Groups] $i] \" ,
            }
            set PutStrings [string trimright $PutStrings ,]
            append PutStrings \]
            puts $MyFileVar "        \"groups\": $PutStrings,"
        }
        set PutStrings \[
        for {set i 0} {$i < [llength [dict get $BodyVolume Surfaces]]} {incr i} {
            append PutStrings [lindex [dict get $BodyVolume Surfaces] $i] ,
        }
        set PutStrings [string trimright $PutStrings ,]
        append PutStrings \]
        puts $MyFileVar "        \"surfaces\":  $PutStrings"
        if {$iter < [dict size $BodyVolumesDict]} {
            puts $MyFileVar "    \},\{"
        } else {
            puts $MyFileVar "    \}\],"
        }
    }
}

#-------------------------------------------------------------------------------

proc WriteFracturesList {FileVar FracturesDict} {
    upvar $FileVar MyFileVar

    set iter 0
    puts $MyFileVar "    \"fractures_list\": \[\{"
    dict for {Id Fracture} $FracturesDict {
        incr iter
        puts $MyFileVar "        \"id\":                $Id,"
        # TipPoint
        puts $MyFileVar "        \"tip_point\":         \{"
        puts $MyFileVar "            \"id\":          [dict get $Fracture TipPoint Id],"
        puts $MyFileVar "            \"layer\":       \"[dict get $Fracture TipPoint Layer]\","
        if {[llength [dict get $Fracture TipPoint Groups]]==0} {
            puts $MyFileVar "            \"groups\":      \[\],"
        } else {
            set PutStrings \[
            for {set i 0} {$i < [llength [dict get $Fracture TipPoint Groups]]} {incr i} {
                append PutStrings \" [lindex [dict get $Fracture TipPoint Groups] $i] \" ,
            }
            set PutStrings [string trimright $PutStrings ,]
            append PutStrings \]
            puts $MyFileVar "            \"groups\":      $PutStrings,"
        }
        puts $MyFileVar "            \"coordinates\": \[[lindex [dict get $Fracture TipPoint Coordinates] 0],[lindex [dict get $Fracture TipPoint Coordinates] 1],[lindex [dict get $Fracture TipPoint Coordinates] 2]\]"
        puts $MyFileVar "        \},"
        # TopPoint
        puts $MyFileVar "        \"top_point\":         \{"
        puts $MyFileVar "            \"id\":          [dict get $Fracture TopPoint Id],"
        puts $MyFileVar "            \"layer\":       \"[dict get $Fracture TopPoint Layer]\","
        if {[llength [dict get $Fracture TopPoint Groups]]==0} {
            puts $MyFileVar "            \"groups\":      \[\],"
        } else {
            set PutStrings \[
            for {set i 0} {$i < [llength [dict get $Fracture TopPoint Groups]]} {incr i} {
                append PutStrings \" [lindex [dict get $Fracture TopPoint Groups] $i] \" ,
            }
            set PutStrings [string trimright $PutStrings ,]
            append PutStrings \]
            puts $MyFileVar "            \"groups\":      $PutStrings,"
        }
        puts $MyFileVar "            \"coordinates\": \[[lindex [dict get $Fracture TopPoint Coordinates] 0],[lindex [dict get $Fracture TopPoint Coordinates] 1],[lindex [dict get $Fracture TopPoint Coordinates] 2]\]"
        puts $MyFileVar "        \},"
        # BotPoint
        puts $MyFileVar "        \"bot_point\":         \{"
        puts $MyFileVar "            \"id\":          [dict get $Fracture BotPoint Id],"
        puts $MyFileVar "            \"layer\":       \"[dict get $Fracture BotPoint Layer]\","
        if {[llength [dict get $Fracture BotPoint Groups]]==0} {
            puts $MyFileVar "            \"groups\":      \[\],"
        } else {
            set PutStrings \[
            for {set i 0} {$i < [llength [dict get $Fracture BotPoint Groups]]} {incr i} {
                append PutStrings \" [lindex [dict get $Fracture BotPoint Groups] $i] \" ,
            }
            set PutStrings [string trimright $PutStrings ,]
            append PutStrings \]
            puts $MyFileVar "            \"groups\":      $PutStrings,"
        }
        puts $MyFileVar "            \"coordinates\": \[[lindex [dict get $Fracture BotPoint Coordinates] 0],[lindex [dict get $Fracture BotPoint Coordinates] 1],[lindex [dict get $Fracture BotPoint Coordinates] 2]\]"
        puts $MyFileVar "        \},"
        # LeftPoint
        puts $MyFileVar "        \"left_point\":         \{"
        puts $MyFileVar "            \"id\":          [dict get $Fracture LeftPoint Id],"
        puts $MyFileVar "            \"layer\":       \"[dict get $Fracture LeftPoint Layer]\","
        if {[llength [dict get $Fracture LeftPoint Groups]]==0} {
            puts $MyFileVar "            \"groups\":      \[\],"
        } else {
            set PutStrings \[
            for {set i 0} {$i < [llength [dict get $Fracture LeftPoint Groups]]} {incr i} {
                append PutStrings \" [lindex [dict get $Fracture LeftPoint Groups] $i] \" ,
            }
            set PutStrings [string trimright $PutStrings ,]
            append PutStrings \]
            puts $MyFileVar "            \"groups\":      $PutStrings,"
        }
        puts $MyFileVar "            \"coordinates\": \[[lindex [dict get $Fracture LeftPoint Coordinates] 0],[lindex [dict get $Fracture LeftPoint Coordinates] 1],[lindex [dict get $Fracture LeftPoint Coordinates] 2]\]"
        puts $MyFileVar "        \},"
        # RightPoint
        puts $MyFileVar "        \"right_point\":         \{"
        puts $MyFileVar "            \"id\":          [dict get $Fracture RightPoint Id],"
        puts $MyFileVar "            \"layer\":       \"[dict get $Fracture RightPoint Layer]\","
        if {[llength [dict get $Fracture RightPoint Groups]]==0} {
            puts $MyFileVar "            \"groups\":      \[\],"
        } else {
            set PutStrings \[
            for {set i 0} {$i < [llength [dict get $Fracture RightPoint Groups]]} {incr i} {
                append PutStrings \" [lindex [dict get $Fracture RightPoint Groups] $i] \" ,
            }
            set PutStrings [string trimright $PutStrings ,]
            append PutStrings \]
            puts $MyFileVar "            \"groups\":      $PutStrings,"
        }
        puts $MyFileVar "            \"coordinates\": \[[lindex [dict get $Fracture RightPoint Coordinates] 0],[lindex [dict get $Fracture RightPoint Coordinates] 1],[lindex [dict get $Fracture RightPoint Coordinates] 2]\]"
        puts $MyFileVar "        \},"
        # TopLine
        puts $MyFileVar "        \"top_line\":          \{"
        puts $MyFileVar "            \"id\":         [dict get $Fracture TopLine Id],"
        puts $MyFileVar "            \"layer\":      \"[dict get $Fracture TopLine Layer]\","
        if {[llength [dict get $Fracture TopLine Groups]]==0} {
            puts $MyFileVar "            \"groups\":     \[\]"
        } else {
            set PutStrings \[
            for {set i 0} {$i < [llength [dict get $Fracture TopLine Groups]]} {incr i} {
                append PutStrings \" [lindex [dict get $Fracture TopLine Groups] $i] \" ,
            }
            set PutStrings [string trimright $PutStrings ,]
            append PutStrings \]
            puts $MyFileVar "            \"groups\":     $PutStrings"
        }
        puts $MyFileVar "        \},"
        # BotLine
        puts $MyFileVar "        \"bot_line\":          \{"
        puts $MyFileVar "            \"id\":         [dict get $Fracture BotLine Id],"
        puts $MyFileVar "            \"layer\":      \"[dict get $Fracture BotLine Layer]\","
        if {[llength [dict get $Fracture BotLine Groups]]==0} {
            puts $MyFileVar "            \"groups\":     \[\]"
        } else {
            set PutStrings \[
            for {set i 0} {$i < [llength [dict get $Fracture BotLine Groups]]} {incr i} {
                append PutStrings \" [lindex [dict get $Fracture BotLine Groups] $i] \" ,
            }
            set PutStrings [string trimright $PutStrings ,]
            append PutStrings \]
            puts $MyFileVar "            \"groups\":     $PutStrings"
        }
        puts $MyFileVar "        \},"
        # LeftLine
        puts $MyFileVar "        \"left_line\":          \{"
        puts $MyFileVar "            \"id\":         [dict get $Fracture LeftLine Id],"
        puts $MyFileVar "            \"layer\":      \"[dict get $Fracture LeftLine Layer]\","
        if {[llength [dict get $Fracture LeftLine Groups]]==0} {
            puts $MyFileVar "            \"groups\":     \[\]"
        } else {
            set PutStrings \[
            for {set i 0} {$i < [llength [dict get $Fracture LeftLine Groups]]} {incr i} {
                append PutStrings \" [lindex [dict get $Fracture LeftLine Groups] $i] \" ,
            }
            set PutStrings [string trimright $PutStrings ,]
            append PutStrings \]
            puts $MyFileVar "            \"groups\":     $PutStrings"
        }
        puts $MyFileVar "        \},"
        # RightLine
        puts $MyFileVar "        \"right_line\":          \{"
        puts $MyFileVar "            \"id\":         [dict get $Fracture RightLine Id],"
        puts $MyFileVar "            \"layer\":      \"[dict get $Fracture RightLine Layer]\","
        if {[llength [dict get $Fracture RightLine Groups]]==0} {
            puts $MyFileVar "            \"groups\":     \[\]"
        } else {
            set PutStrings \[
            for {set i 0} {$i < [llength [dict get $Fracture RightLine Groups]]} {incr i} {
                append PutStrings \" [lindex [dict get $Fracture RightLine Groups] $i] \" ,
            }
            set PutStrings [string trimright $PutStrings ,]
            append PutStrings \]
            puts $MyFileVar "            \"groups\":     $PutStrings"
        }
        puts $MyFileVar "        \},"
        # TopLeftSurface
        puts $MyFileVar "        \"top_left_surface\": \{"
        puts $MyFileVar "            \"id\":       [dict get $Fracture TopLeftSurface Id],"
        puts $MyFileVar "            \"layer\":    \"[dict get $Fracture TopLeftSurface Layer]\","
        if {[llength [dict get $Fracture TopLeftSurface Groups]]==0} {
            puts $MyFileVar "            \"groups\":   \[\]"
        } else {
            set PutStrings \[
            for {set i 0} {$i < [llength [dict get $Fracture TopLeftSurface Groups]]} {incr i} {
                append PutStrings \" [lindex [dict get $Fracture TopLeftSurface Groups] $i] \" ,
            }
            set PutStrings [string trimright $PutStrings ,]
            append PutStrings \]
            puts $MyFileVar "            \"groups\":   $PutStrings"
        }
        puts $MyFileVar "        \},"
        # TopRightSurface
        puts $MyFileVar "        \"top_right_surface\": \{"
        puts $MyFileVar "            \"id\":       [dict get $Fracture TopRightSurface Id],"
        puts $MyFileVar "            \"layer\":    \"[dict get $Fracture TopRightSurface Layer]\","
        if {[llength [dict get $Fracture TopRightSurface Groups]]==0} {
            puts $MyFileVar "            \"groups\":   \[\]"
        } else {
            set PutStrings \[
            for {set i 0} {$i < [llength [dict get $Fracture TopRightSurface Groups]]} {incr i} {
                append PutStrings \" [lindex [dict get $Fracture TopRightSurface Groups] $i] \" ,
            }
            set PutStrings [string trimright $PutStrings ,]
            append PutStrings \]
            puts $MyFileVar "            \"groups\":   $PutStrings"
        }
        puts $MyFileVar "        \},"
        # BotLeftSurface
        puts $MyFileVar "        \"bot_left_surface\": \{"
        puts $MyFileVar "            \"id\":       [dict get $Fracture BotLeftSurface Id],"
        puts $MyFileVar "            \"layer\":    \"[dict get $Fracture BotLeftSurface Layer]\","
        if {[llength [dict get $Fracture BotLeftSurface Groups]]==0} {
            puts $MyFileVar "            \"groups\":   \[\]"
        } else {
            set PutStrings \[
            for {set i 0} {$i < [llength [dict get $Fracture BotLeftSurface Groups]]} {incr i} {
                append PutStrings \" [lindex [dict get $Fracture BotLeftSurface Groups] $i] \" ,
            }
            set PutStrings [string trimright $PutStrings ,]
            append PutStrings \]
            puts $MyFileVar "            \"groups\":   $PutStrings"
        }
        puts $MyFileVar "        \},"
        # BotRightSurface
        puts $MyFileVar "        \"bot_right_surface\": \{"
        puts $MyFileVar "            \"id\":       [dict get $Fracture BotRightSurface Id],"
        puts $MyFileVar "            \"layer\":    \"[dict get $Fracture BotRightSurface Layer]\","
        if {[llength [dict get $Fracture BotRightSurface Groups]]==0} {
            puts $MyFileVar "            \"groups\":   \[\]"
        } else {
            set PutStrings \[
            for {set i 0} {$i < [llength [dict get $Fracture BotRightSurface Groups]]} {incr i} {
                append PutStrings \" [lindex [dict get $Fracture BotRightSurface Groups] $i] \" ,
            }
            set PutStrings [string trimright $PutStrings ,]
            append PutStrings \]
            puts $MyFileVar "            \"groups\":   $PutStrings"
        }
        puts $MyFileVar "        \},"
        # LeftInterfaceVolume
        puts $MyFileVar "        \"left_interface_volume\": \{"
        puts $MyFileVar "            \"id\":       [dict get $Fracture LeftVolume Id],"
        puts $MyFileVar "            \"layer\":    \"[dict get $Fracture LeftVolume Layer]\","
        if {[llength [dict get $Fracture LeftVolume Groups]]==0} {
            puts $MyFileVar "            \"groups\":   \[\]"
        } else {
            set PutStrings \[
            for {set i 0} {$i < [llength [dict get $Fracture LeftVolume Groups]]} {incr i} {
                append PutStrings \" [lindex [dict get $Fracture LeftVolume Groups] $i] \" ,
            }
            set PutStrings [string trimright $PutStrings ,]
            append PutStrings \]
            puts $MyFileVar "            \"groups\":   $PutStrings"
        }
        puts $MyFileVar "        \},"
        # RightInterfaceVolume
        puts $MyFileVar "        \"right_interface_volume\": \{"
        puts $MyFileVar "            \"id\":       [dict get $Fracture RightVolume Id],"
        puts $MyFileVar "            \"layer\":    \"[dict get $Fracture RightVolume Layer]\","
        if {[llength [dict get $Fracture RightVolume Groups]]==0} {
            puts $MyFileVar "            \"groups\":   \[\]"
        } else {
            set PutStrings \[
            for {set i 0} {$i < [llength [dict get $Fracture RightVolume Groups]]} {incr i} {
                append PutStrings \" [lindex [dict get $Fracture RightVolume Groups] $i] \" ,
            }
            set PutStrings [string trimright $PutStrings ,]
            append PutStrings \]
            puts $MyFileVar "            \"groups\":   $PutStrings"
        }
        puts $MyFileVar "        \},"
        # BodyVolumes
        set PutStrings \[
        for {set i 0} {$i < [llength [dict get $Fracture BodyVolumes]]} {incr i} {
            append PutStrings [lindex [dict get $Fracture BodyVolumes] $i] ,
        }
        set PutStrings [string trimright $PutStrings ,]
        append PutStrings \]
        puts $MyFileVar "        \"body_volumes\":     $PutStrings"
        if {$iter < [dict size $FracturesDict]} {
            puts $MyFileVar "    \},\{"
        } else {
            puts $MyFileVar "    \}\]"
        }
    }
}
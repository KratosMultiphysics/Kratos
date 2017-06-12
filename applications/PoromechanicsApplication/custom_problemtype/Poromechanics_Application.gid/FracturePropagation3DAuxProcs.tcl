proc AddPointToFracturesDict {FracturesDict FractureId PointId Point PointType} {
    upvar $FracturesDict MyFracturesDict
    
    dict set MyFracturesDict $FractureId $PointType Id $PointId
    dict set MyFracturesDict $FractureId $PointType Coordinates "[lindex $Point 1] [lindex $Point 2] [lindex $Point 3]"
}

#-------------------------------------------------------------------------------

#~ proc AddLineToFracturesDict {FracturesDict FractureId LineId LineType} {
    #~ upvar $FracturesDict MyFracturesDict
    
    #~ dict set MyFracturesDict $FractureId $LineType Id $LineId
    #~ dict set MyFracturesDict $FractureId $LineType Groups [GiD_EntitiesGroups entity_groups lines $LineId]
#~ }

#-------------------------------------------------------------------------------

proc AddSurfaceToFracturesDict {FracturesDict FractureId SurfaceId Surface SurfaceType} {
    upvar $FracturesDict MyFracturesDict
    
    dict set MyFracturesDict $FractureId $SurfaceType Id $SurfaceId
    set Lines [list]
    for {set i 0} {$i<[lindex $Surface 2]} {incr i} {
        lappend Lines [lindex [lindex $Surface [expr { 9+$i }]] 0]
    }
    dict set MyFracturesDict $FractureId $SurfaceType Lines $Lines
    #dict set MyFracturesDict $FractureId $SurfaceType Groups [GiD_EntitiesGroups entity_groups surfaces $SurfaceId]
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
        if {[GiD_Info IsPointInside Volume $Id $TipCoordinates] eq 1} {
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
    
    # Vector in local x direction
    set Vx(0) [expr {$Tp(0)-$Rp(0)}]
    set Vx(1) [expr {$Tp(1)-$Rp(1)}]
    set Vx(2) [expr {$Tp(2)-$Rp(2)}]

    # Vector in local y direction
    set Vy(0) [expr {$Lp(0)-$Rp(0)}]
    set Vy(1) [expr {$Lp(1)-$Rp(1)}]
    set Vy(2) [expr {$Lp(2)-$Rp(2)}]
    
    # Unitary vector in local z direction (Cross product between Vx and Vy)    
    set Vz(0) [expr {$Vx(1)*$Vy(2)-$Vx(2)*$Vy(1)}]
    set Vz(1) [expr {$Vx(2)*$Vy(0)-$Vx(0)*$Vy(2)}]
    set Vz(2) [expr {$Vx(0)*$Vy(1)-$Vx(1)*$Vy(0)}]
    set InvNorm [expr {1.0/sqrt($Vz(0)*$Vz(0)+$Vz(1)*$Vz(1)+$Vz(2)*$Vz(2))}]
    set Vz(0) [expr {$Vz(0)*$InvNorm}]
    set Vz(1) [expr {$Vz(1)*$InvNorm}]
    set Vz(2) [expr {$Vz(2)*$InvNorm}]
    
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

    puts $MyFileVar "        \"propagation_damage\":                   [GiD_AccessValue get gendata Propagation_Damage],"
    puts $MyFileVar "        \"propagation_length\":                   [GiD_AccessValue get gendata Propagation_Length],"
    puts $MyFileVar "        \"propagation_width\":                    [GiD_AccessValue get gendata Propagation_Width],"
    puts $MyFileVar "        \"propagation_height\":                   [GiD_AccessValue get gendata Propagation_Height],"
    puts $MyFileVar "        \"correction_tolerance\":                 [GiD_AccessValue get gendata Correction_Tolerance],"
    puts $MyFileVar "        \"propagation_frequency\":                [GiD_AccessValue get gendata Propagation_Frequency],"
    # body_domain_sub_model_part_list
    set PutStrings \[
    set Groups [GiD_Info conditions Body_Part groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
        append PutStrings \" [lindex [lindex $Groups $i] 1] \" ,
    }
    set PutStrings [string trimright $PutStrings ,]
    append PutStrings \]
    puts $MyFileVar "        \"body_domain_sub_model_part_list\":      $PutStrings,"
}

#-------------------------------------------------------------------------------

proc WriteBodyVolumesList {FileVar BodyVolumesDict} {
    upvar $FileVar MyFileVar

    set iter 0
    puts $MyFileVar "    \"body_volumes_list\": \[\{"
    dict for {Id BodyVolume} $BodyVolumesDict {
        incr iter
        puts $MyFileVar "        \"id\":        $Id,"
        if {[llength [dict get $BodyVolume Groups]] eq 0} {
            puts $MyFileVar "        \"groups\":    \[\],"
        } else {
            set PutStrings \[
            for {set i 0} {$i < [llength [dict get $BodyVolume Groups]]} {incr i} {
                append PutStrings \" [lindex [dict get $BodyVolume Groups] $i] \" ,
            }
            set PutStrings [string trimright $PutStrings ,]
            append PutStrings \]
            puts $MyFileVar "        \"groups\":    $PutStrings,"
        }
        set PutStrings \[
        for {set i 0} {$i < [llength [dict get $BodyVolume Surfaces]]} {incr i} {
            append PutStrings [lindex [dict get $BodyVolume Surfaces] $i] ,
        }
        set PutStrings [string trimright $PutStrings ,]
        append PutStrings \]
        puts $MyFileVar "        \"surfaces\":  $PutStrings,"
        puts $MyFileVar "        \"mesh_size\": [dict get $BodyVolume MeshSize]"
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
        puts $MyFileVar "        \"id\":                     $Id,"
        # TipPoint
        puts $MyFileVar "        \"tip_point\":              \{"
        puts $MyFileVar "            \"id\":          [dict get $Fracture TipPoint Id],"
        puts $MyFileVar "            \"coordinates\": \[[lindex [dict get $Fracture TipPoint Coordinates] 0],[lindex [dict get $Fracture TipPoint Coordinates] 1],[lindex [dict get $Fracture TipPoint Coordinates] 2]\]"
        puts $MyFileVar "        \},"
        # TopPoint
        puts $MyFileVar "        \"top_point\":              \{"
        puts $MyFileVar "            \"id\":          [dict get $Fracture TopPoint Id],"
        puts $MyFileVar "            \"coordinates\": \[[lindex [dict get $Fracture TopPoint Coordinates] 0],[lindex [dict get $Fracture TopPoint Coordinates] 1],[lindex [dict get $Fracture TopPoint Coordinates] 2]\]"
        puts $MyFileVar "        \},"
        # BotPoint
        puts $MyFileVar "        \"bot_point\":              \{"
        puts $MyFileVar "            \"id\":          [dict get $Fracture BotPoint Id],"
        puts $MyFileVar "            \"coordinates\": \[[lindex [dict get $Fracture BotPoint Coordinates] 0],[lindex [dict get $Fracture BotPoint Coordinates] 1],[lindex [dict get $Fracture BotPoint Coordinates] 2]\]"
        puts $MyFileVar "        \},"
        # LeftPoint
        puts $MyFileVar "        \"left_point\":             \{"
        puts $MyFileVar "            \"id\":          [dict get $Fracture LeftPoint Id],"
        puts $MyFileVar "            \"coordinates\": \[[lindex [dict get $Fracture LeftPoint Coordinates] 0],[lindex [dict get $Fracture LeftPoint Coordinates] 1],[lindex [dict get $Fracture LeftPoint Coordinates] 2]\]"
        puts $MyFileVar "        \},"
        # RightPoint
        puts $MyFileVar "        \"right_point\":            \{"
        puts $MyFileVar "            \"id\":          [dict get $Fracture RightPoint Id],"
        puts $MyFileVar "            \"coordinates\": \[[lindex [dict get $Fracture RightPoint Coordinates] 0],[lindex [dict get $Fracture RightPoint Coordinates] 1],[lindex [dict get $Fracture RightPoint Coordinates] 2]\]"
        puts $MyFileVar "        \},"
        # TopLine
        puts $MyFileVar "        \"top_line\":               \{"
        puts $MyFileVar "            \"id\": [dict get $Fracture TopLine Id]"
        puts $MyFileVar "        \},"
        # BotLine
        puts $MyFileVar "        \"bot_line\":               \{"
        puts $MyFileVar "            \"id\": [dict get $Fracture BotLine Id]"
        puts $MyFileVar "        \},"
        # LeftLine
        puts $MyFileVar "        \"left_line\":              \{"
        puts $MyFileVar "            \"id\": [dict get $Fracture LeftLine Id]"
        puts $MyFileVar "        \},"
        # RightLine
        puts $MyFileVar "        \"right_line\":             \{"
        puts $MyFileVar "            \"id\": [dict get $Fracture RightLine Id]"
        puts $MyFileVar "        \},"
        # TopLeftSurface
        puts $MyFileVar "        \"top_left_surface\":       \{"
        puts $MyFileVar "            \"id\":    [dict get $Fracture TopLeftSurface Id],"
        set PutStrings \[
        for {set i 0} {$i < [llength [dict get $Fracture TopLeftSurface Lines]]} {incr i} {
            append PutStrings [lindex [dict get $Fracture TopLeftSurface Lines] $i] ,
        }
        set PutStrings [string trimright $PutStrings ,]
        append PutStrings \]
        puts $MyFileVar "            \"lines\": $PutStrings"
        puts $MyFileVar "        \},"
        # TopRightSurface
        puts $MyFileVar "        \"top_right_surface\":      \{"
        puts $MyFileVar "            \"id\":    [dict get $Fracture TopRightSurface Id],"
        set PutStrings \[
        for {set i 0} {$i < [llength [dict get $Fracture TopRightSurface Lines]]} {incr i} {
            append PutStrings [lindex [dict get $Fracture TopRightSurface Lines] $i] ,
        }
        set PutStrings [string trimright $PutStrings ,]
        append PutStrings \]
        puts $MyFileVar "            \"lines\": $PutStrings"
        puts $MyFileVar "        \},"
        # BotLeftSurface
        puts $MyFileVar "        \"bot_left_surface\":       \{"
        puts $MyFileVar "            \"id\":    [dict get $Fracture BotLeftSurface Id],"
        set PutStrings \[
        for {set i 0} {$i < [llength [dict get $Fracture BotLeftSurface Lines]]} {incr i} {
            append PutStrings [lindex [dict get $Fracture BotLeftSurface Lines] $i] ,
        }
        set PutStrings [string trimright $PutStrings ,]
        append PutStrings \]
        puts $MyFileVar "            \"lines\": $PutStrings"
        puts $MyFileVar "        \},"
        # BotRightSurface
        puts $MyFileVar "        \"bot_right_surface\":      \{"
        puts $MyFileVar "            \"id\":    [dict get $Fracture BotRightSurface Id],"
        set PutStrings \[
        for {set i 0} {$i < [llength [dict get $Fracture BotRightSurface Lines]]} {incr i} {
            append PutStrings [lindex [dict get $Fracture BotRightSurface Lines] $i] ,
        }
        set PutStrings [string trimright $PutStrings ,]
        append PutStrings \]
        puts $MyFileVar "            \"lines\": $PutStrings"
        puts $MyFileVar "        \},"
        # LeftInterfaceVolume
        puts $MyFileVar "        \"left_interface_volume\":  \{"
        puts $MyFileVar "            \"id\":     [dict get $Fracture LeftInterfaceVolume Id],"
        puts $MyFileVar "            \"layer\":  \"[dict get $Fracture LeftInterfaceVolume Layer]\","
        if {[llength [dict get $Fracture LeftInterfaceVolume Groups]] eq 0} {
            puts $MyFileVar "            \"groups\": \[\]"
        } else {
            set PutStrings \[
            for {set i 0} {$i < [llength [dict get $Fracture LeftInterfaceVolume Groups]]} {incr i} {
                append PutStrings \" [lindex [dict get $Fracture LeftInterfaceVolume Groups] $i] \" ,
            }
            set PutStrings [string trimright $PutStrings ,]
            append PutStrings \]
            puts $MyFileVar "            \"groups\": $PutStrings"
        }
        puts $MyFileVar "        \},"
        # RightInterfaceVolume
        puts $MyFileVar "        \"right_interface_volume\": \{"
        puts $MyFileVar "            \"id\":     [dict get $Fracture RightInterfaceVolume Id],"
        puts $MyFileVar "            \"layer\":  \"[dict get $Fracture RightInterfaceVolume Layer]\","
        if {[llength [dict get $Fracture RightInterfaceVolume Groups]] eq 0} {
            puts $MyFileVar "            \"groups\": \[\]"
        } else {
            set PutStrings \[
            for {set i 0} {$i < [llength [dict get $Fracture RightInterfaceVolume Groups]]} {incr i} {
                append PutStrings \" [lindex [dict get $Fracture RightInterfaceVolume Groups] $i] \" ,
            }
            set PutStrings [string trimright $PutStrings ,]
            append PutStrings \]
            puts $MyFileVar "            \"groups\": $PutStrings"
        }
        puts $MyFileVar "        \},"
        # BodyVolumes
        set PutStrings \[
        for {set i 0} {$i < [llength [dict get $Fracture BodyVolumes]]} {incr i} {
            append PutStrings [lindex [dict get $Fracture BodyVolumes] $i] ,
        }
        set PutStrings [string trimright $PutStrings ,]
        append PutStrings \]
        puts $MyFileVar "        \"body_volumes\":           $PutStrings"
        if {$iter < [dict size $FracturesDict]} {
            puts $MyFileVar "    \},\{"
        } else {
            puts $MyFileVar "    \}\]"
        }
    }
}

#--------------------------------------------------------------------------------------------------------------------------------------------------------------

proc ComputeNormal {VertexCoord Point1Coord Point2Coord} {
    # Vector in local x direction
    set Vx(0) [expr {[lindex $Point1Coord 0]-[lindex $VertexCoord 0]}]
    set Vx(1) [expr {[lindex $Point1Coord 1]-[lindex $VertexCoord 1]}]
    set Vx(2) [expr {[lindex $Point1Coord 2]-[lindex $VertexCoord 2]}]

    # Vector in local y direction
    set Vy(0) [expr {[lindex $Point2Coord 0]-[lindex $VertexCoord 0]}]
    set Vy(1) [expr {[lindex $Point2Coord 1]-[lindex $VertexCoord 1]}]
    set Vy(2) [expr {[lindex $Point2Coord 2]-[lindex $VertexCoord 2]}]
    
    # Vector in local z direction (Cross product between Vx and Vy)    
    set Vz(0) [expr {$Vx(1)*$Vy(2)-$Vx(2)*$Vy(1)}]
    set Vz(1) [expr {$Vx(2)*$Vy(0)-$Vx(0)*$Vy(2)}]
    set Vz(2) [expr {$Vx(0)*$Vy(1)-$Vx(1)*$Vy(0)}]
    
    return "$Vz(0) $Vz(1) $Vz(2)"
}

#-------------------------------------------------------------------------------

proc ComputeTransformMatrix {VertexCoord Point1Coord Point2Coord InitAxisCoord FinalAxisCoord} {
    # Unitary vector at the initial position of the Rotation
    set Ri(0) [expr {[lindex $Point1Coord 0]-[lindex $VertexCoord 0]}]
    set Ri(1) [expr {[lindex $Point1Coord 1]-[lindex $VertexCoord 1]}]
    set Ri(2) [expr {[lindex $Point1Coord 2]-[lindex $VertexCoord 2]}]
    set InvNorm [expr {1.0/sqrt($Ri(0)*$Ri(0)+$Ri(1)*$Ri(1)+$Ri(2)*$Ri(2))}]
    set Ri(0) [expr {$Ri(0)*$InvNorm}]
    set Ri(1) [expr {$Ri(1)*$InvNorm}]
    set Ri(2) [expr {$Ri(2)*$InvNorm}]
    # Unitary vector at the final position of the Rotation
    set Rf(0) [expr {[lindex $Point2Coord 0]-[lindex $VertexCoord 0]}]
    set Rf(1) [expr {[lindex $Point2Coord 1]-[lindex $VertexCoord 1]}]
    set Rf(2) [expr {[lindex $Point2Coord 2]-[lindex $VertexCoord 2]}]
    set InvNorm [expr {1.0/sqrt($Rf(0)*$Rf(0)+$Rf(1)*$Rf(1)+$Rf(2)*$Rf(2))}]
    set Rf(0) [expr {$Rf(0)*$InvNorm}]
    set Rf(1) [expr {$Rf(1)*$InvNorm}]
    set Rf(2) [expr {$Rf(2)*$InvNorm}]
    # Unitary rotation Axis
    set A(0) [expr {[lindex $FinalAxisCoord 0]-[lindex $InitAxisCoord 0]}]
    set A(1) [expr {[lindex $FinalAxisCoord 1]-[lindex $InitAxisCoord 1]}]
    set A(2) [expr {[lindex $FinalAxisCoord 2]-[lindex $InitAxisCoord 2]}]
    set InvNorm [expr {1.0/sqrt($A(0)*$A(0)+$A(1)*$A(1)+$A(2)*$A(2))}]
    set A(0) [expr {$A(0)*$InvNorm}]
    set A(1) [expr {$A(1)*$InvNorm}]
    set A(2) [expr {$A(2)*$InvNorm}]
    
    # Cosine of the angle of rotation
    set CosAngle [expr {$Ri(0)*$Rf(0)+$Ri(1)*$Rf(1)+$Ri(2)*$Rf(2)}]
    # Cross product between vectors Ri and Rf
    set n(0) [expr {$Ri(1)*$Rf(2)-$Ri(2)*$Rf(1)}]
    set n(1) [expr {$Ri(2)*$Rf(0)-$Ri(0)*$Rf(2)}]
    set n(2) [expr {$Ri(0)*$Rf(1)-$Ri(1)*$Rf(0)}]
    # Sine of the angle of rotation (positive angle between 0º and 90º)
    set SinAngle [expr {sqrt($n(0)*$n(0)+$n(1)*$n(1)+$n(2)*$n(2))}]
    
    # Transformation Matrix
    set Rxx [expr {$CosAngle+$A(0)*$A(0)*(1.0-$CosAngle)}]
    set Rxy [expr {$A(0)*$A(1)*(1.0-$CosAngle)-$A(2)*$SinAngle}]
    set Rxz [expr {$A(0)*$A(2)*(1.0-$CosAngle)+$A(1)*$SinAngle}]
    set Ryx [expr {$A(0)*$A(1)*(1.0-$CosAngle)+$A(2)*$SinAngle}]
    set Ryy [expr {$CosAngle+$A(1)*$A(1)*(1.0-$CosAngle)}]
    set Ryz [expr {$A(1)*$A(2)*(1.0-$CosAngle)-$A(0)*$SinAngle}]
    set Rzx [expr {$A(0)*$A(2)*(1.0-$CosAngle)-$A(1)*$SinAngle}]
    set Rzy [expr {$A(1)*$A(2)*(1.0-$CosAngle)+$A(0)*$SinAngle}]
    set Rzz [expr {$CosAngle+$A(2)*$A(2)*(1.0-$CosAngle)}]
    set Tx [expr {[lindex $InitAxisCoord 0]-($Rxx*[lindex $InitAxisCoord 0]+$Rxy*[lindex $InitAxisCoord 1]+$Rxz*[lindex $InitAxisCoord 2])}]
    set Ty [expr {[lindex $InitAxisCoord 1]-($Ryx*[lindex $InitAxisCoord 0]+$Ryy*[lindex $InitAxisCoord 1]+$Ryz*[lindex $InitAxisCoord 2])}]
    set Tz [expr {[lindex $InitAxisCoord 2]-($Rzx*[lindex $InitAxisCoord 0]+$Rzy*[lindex $InitAxisCoord 1]+$Rzz*[lindex $InitAxisCoord 2])}]
    
    return [list $Rxx $Rxy $Rxz $Tx \
                 $Ryx $Ryy $Ryz $Ty \
                 $Rzx $Rzy $Rzz $Tz \
                 0.0 0.0 0.0 1.0]
}

#-------------------------------------------------------------------------------

proc AddPropagationUnionPoint {NumPropUnionGroups PointId} {
    
    upvar $NumPropUnionGroups MyNumPropUnionGroups
    
    GiD_Groups create PropagationUnion_3d_6//SG$MyNumPropUnionGroups
    GiD_EntitiesGroups assign PropagationUnion_3d_6//SG$MyNumPropUnionGroups points $PointId
    incr MyNumPropUnionGroups
}
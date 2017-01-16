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

proc AddInterfaceSurfaceToFracturesDict {FracturesDict FractureId SurfaceId InterfaceSurface} {
    upvar $FracturesDict MyFracturesDict
    
    dict set MyFracturesDict $FractureId InterfaceSurface Id $SurfaceId
    dict set MyFracturesDict $FractureId InterfaceSurface Layer [lindex $InterfaceSurface 1]
    dict set MyFracturesDict $FractureId InterfaceSurface Groups [GiD_EntitiesGroups entity_groups surfaces $SurfaceId]
}

#-------------------------------------------------------------------------------

proc AddBodySurfaceToFracturesDict {FracturesDict FractureId BodySurfacesDict} {
    upvar $FracturesDict MyFracturesDict

    set TipCoordinates [dict get $MyFracturesDict $FractureId TipPoint Coordinates]
    set BodySurfaces [list]
    dict for {Id BodySurface} $BodySurfacesDict {
        if {[GiD_Info IsPointInside Surface $Id $TipCoordinates] eq 1} {
            lappend BodySurfaces $Id
        }
    }
    dict set MyFracturesDict $FractureId BodySurfaces $BodySurfaces
}


#--------------------------------------------------------------------------------------------------------------------------------------------------------------


proc WriteFractureData {FileVar} {
    upvar $FileVar MyFileVar
    
    puts $MyFileVar "        \"propagation_damage\":                   [GiD_AccessValue get gendata Propagation_Damage],"
    puts $MyFileVar "        \"propagation_length\":                   [GiD_AccessValue get gendata Propagation_Length],"
    puts $MyFileVar "        \"propagation_width\":                    [GiD_AccessValue get gendata Propagation_Width],"
    puts $MyFileVar "        \"propagation_cosangle\":                 [GiD_AccessValue get gendata Propagation_CosAngle],"
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

proc WriteBodySurfacesList {FileVar BodySurfacesDict} {
    upvar $FileVar MyFileVar

    set iter 0
    puts $MyFileVar "    \"body_surfaces_list\": \[\{"
    dict for {Id BodySurface} $BodySurfacesDict {
        incr iter
        puts $MyFileVar "        \"id\":     $Id,"
        if {[llength [dict get $BodySurface Groups]] eq 0} {
            puts $MyFileVar "        \"groups\": \[\],"
        } else {
            set PutStrings \[
            for {set i 0} {$i < [llength [dict get $BodySurface Groups]]} {incr i} {
                append PutStrings \" [lindex [dict get $BodySurface Groups] $i] \" ,
            }
            set PutStrings [string trimright $PutStrings ,]
            append PutStrings \]
            puts $MyFileVar "        \"groups\": $PutStrings,"
        }
        set PutStrings \[
        for {set i 0} {$i < [llength [dict get $BodySurface Lines]]} {incr i} {
            append PutStrings [lindex [dict get $BodySurface Lines] $i] ,
        }
        set PutStrings [string trimright $PutStrings ,]
        append PutStrings \]
        puts $MyFileVar "        \"lines\":  $PutStrings"
        if {$iter < [dict size $BodySurfacesDict]} {
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
        puts $MyFileVar "            \"coordinates\": \[[lindex [dict get $Fracture TipPoint Coordinates] 0],[lindex [dict get $Fracture TipPoint Coordinates] 1],[lindex [dict get $Fracture TipPoint Coordinates] 2]\]"
        puts $MyFileVar "        \},"
        # TopPoint
        puts $MyFileVar "        \"top_point\":         \{"
        puts $MyFileVar "            \"id\":          [dict get $Fracture TopPoint Id],"
        puts $MyFileVar "            \"coordinates\": \[[lindex [dict get $Fracture TopPoint Coordinates] 0],[lindex [dict get $Fracture TopPoint Coordinates] 1],[lindex [dict get $Fracture TopPoint Coordinates] 2]\]"
        puts $MyFileVar "        \},"
        # BotPoint
        puts $MyFileVar "        \"bot_point\":         \{"
        puts $MyFileVar "            \"id\":          [dict get $Fracture BotPoint Id],"
        puts $MyFileVar "            \"coordinates\": \[[lindex [dict get $Fracture BotPoint Coordinates] 0],[lindex [dict get $Fracture BotPoint Coordinates] 1],[lindex [dict get $Fracture BotPoint Coordinates] 2]\]"
        puts $MyFileVar "        \},"
        # TopLine
        puts $MyFileVar "        \"top_line\":          \{"
        puts $MyFileVar "            \"id\": [dict get $Fracture TopLine Id]"
        puts $MyFileVar "        \},"
        # BotLine
        puts $MyFileVar "        \"bot_line\":          \{"
        puts $MyFileVar "            \"id\": [dict get $Fracture BotLine Id]"
        puts $MyFileVar "        \},"
        # InterfaceSurface
        puts $MyFileVar "        \"interface_surface\": \{"
        puts $MyFileVar "            \"id\":     [dict get $Fracture InterfaceSurface Id],"
        puts $MyFileVar "            \"layer\":  \"[dict get $Fracture InterfaceSurface Layer]\","
        if {[llength [dict get $Fracture InterfaceSurface Groups]] eq 0} {
            puts $MyFileVar "            \"groups\": \[\]"
        } else {
            set PutStrings \[
            for {set i 0} {$i < [llength [dict get $Fracture InterfaceSurface Groups]]} {incr i} {
                append PutStrings \" [lindex [dict get $Fracture InterfaceSurface Groups] $i] \" ,
            }
            set PutStrings [string trimright $PutStrings ,]
            append PutStrings \]
            puts $MyFileVar "            \"groups\": $PutStrings"
        }
        puts $MyFileVar "        \},"
        # BodySurfaces
        set PutStrings \[
        for {set i 0} {$i < [llength [dict get $Fracture BodySurfaces]]} {incr i} {
            append PutStrings [lindex [dict get $Fracture BodySurfaces] $i] ,
        }
        set PutStrings [string trimright $PutStrings ,]
        append PutStrings \]
        puts $MyFileVar "        \"body_surfaces\":     $PutStrings"
        if {$iter < [dict size $FracturesDict]} {
            puts $MyFileVar "    \},\{"
        } else {
            puts $MyFileVar "    \}\]"
        }
    }
}
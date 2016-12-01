proc AddPointToFracturesDict {FracturesDict FractureId Line NumPoint PointType} {
    upvar $FracturesDict MyFracturesDict
        
    set Point [GiD_Geometry get point [lindex $Line $NumPoint]]
    set Groups [GiD_EntitiesGroups entity_groups points [lindex $Line $NumPoint]]
    set Coordinates "[lindex $Point 1] [lindex $Point 2] [lindex $Point 3]"
    dict set MyFracturesDict $FractureId $PointType Id [lindex $Line $NumPoint]
    dict set MyFracturesDict $FractureId $PointType Layer [lindex $Point 0]
    dict set MyFracturesDict $FractureId $PointType Groups $Groups
    dict set MyFracturesDict $FractureId $PointType Coordinates $Coordinates
}

#-------------------------------------------------------------------------------

proc AddLineToFracturesDict {FracturesDict FractureId InterfaceSurface NumLine Line LineType} {
    upvar $FracturesDict MyFracturesDict

    set Groups [GiD_EntitiesGroups entity_groups lines [lindex [lindex $InterfaceSurface $NumLine] 0]]
    dict set MyFracturesDict $FractureId $LineType Id [lindex [lindex $InterfaceSurface $NumLine] 0]
    dict set MyFracturesDict $FractureId $LineType Layer [lindex $Line 1]
    dict set MyFracturesDict $FractureId $LineType Groups $Groups
    dict set MyFracturesDict $FractureId $LineType InitPoint [lindex $Line 2]
    dict set MyFracturesDict $FractureId $LineType EndPoint [lindex $Line 3]
}

#-------------------------------------------------------------------------------

proc AddInterfaceSurfaceToFracturesDict {FracturesDict FractureId InterfaceEntity InterfaceSurface} {
    upvar $FracturesDict MyFracturesDict

    set Groups [GiD_EntitiesGroups entity_groups surfaces $InterfaceEntity]
    dict set MyFracturesDict $FractureId InterfaceSurface Id $InterfaceEntity
    dict set MyFracturesDict $FractureId InterfaceSurface Layer [lindex $InterfaceSurface 1]
    dict set MyFracturesDict $FractureId InterfaceSurface Groups $Groups
    dict set MyFracturesDict $FractureId InterfaceSurface TopLine [lindex $InterfaceSurface 4]
    dict set MyFracturesDict $FractureId InterfaceSurface BotLine [lindex $InterfaceSurface 3]
}

#-------------------------------------------------------------------------------

proc AddBodySurfaceToFracturesDict {FracturesDict FractureId BodySurfacesDict} {
    upvar $FracturesDict MyFracturesDict

    set TipCoordinates [dict get $MyFracturesDict $FractureId TipPoint Coordinates]
    set BodySurfaces [list]
    dict for {Id BodySurface} $BodySurfacesDict {
        if {[GiD_Info IsPointInside Surface $Id $TipCoordinates] == 1} {
            lappend BodySurfaces $Id
        }
    }
    dict set MyFracturesDict $FractureId BodySurfaces $BodySurfaces
}


#--------------------------------------------------------------------------------------------------------------------------------------------------------------


proc WriteFractureData {FileVar} {
    upvar $FileVar MyFileVar

    puts $MyFileVar "        \"propagation_length\":                       [GiD_AccessValue get gendata Propagation_Length],"
    puts $MyFileVar "        \"propagation_damage\":                       [GiD_AccessValue get gendata Propagation_Damage],"
    puts $MyFileVar "        \"propagation_width\":                        [GiD_AccessValue get gendata Propagation_Width],"
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

proc WriteBodySurfacesList {FileVar BodySurfacesDict} {
    upvar $FileVar MyFileVar

    set iter 0
    puts $MyFileVar "    \"body_surfaces_list\": \[\{"
    dict for {Id BodySurface} $BodySurfacesDict {
        incr iter
        puts $MyFileVar "        \"id\":     $Id,"
        puts $MyFileVar "        \"layer\":  \"[dict get $BodySurface Layer]\","
        if {[llength [dict get $BodySurface Groups]]==0} {
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
        # TopLine
        puts $MyFileVar "        \"top_line\":          \{"
        puts $MyFileVar "            \"id\":         [dict get $Fracture TopLine Id],"
        puts $MyFileVar "            \"layer\":      \"[dict get $Fracture TopLine Layer]\","
        if {[llength [dict get $Fracture TopLine Groups]]==0} {
            puts $MyFileVar "            \"groups\":     \[\],"
        } else {
            set PutStrings \[
            for {set i 0} {$i < [llength [dict get $Fracture TopLine Groups]]} {incr i} {
                append PutStrings \" [lindex [dict get $Fracture TopLine Groups] $i] \" ,
            }
            set PutStrings [string trimright $PutStrings ,]
            append PutStrings \]
            puts $MyFileVar "            \"groups\":     $PutStrings,"
        }
        puts $MyFileVar "            \"init_point\": [dict get $Fracture TopLine InitPoint],"
        puts $MyFileVar "            \"end_point\":  [dict get $Fracture TopLine EndPoint]"
        puts $MyFileVar "        \},"
        # BotLine
        puts $MyFileVar "        \"bot_line\":          \{"
        puts $MyFileVar "            \"id\":         [dict get $Fracture BotLine Id],"
        puts $MyFileVar "            \"layer\":      \"[dict get $Fracture BotLine Layer]\","
        if {[llength [dict get $Fracture BotLine Groups]]==0} {
            puts $MyFileVar "            \"groups\":     \[\],"
        } else {
            set PutStrings \[
            for {set i 0} {$i < [llength [dict get $Fracture BotLine Groups]]} {incr i} {
                append PutStrings \" [lindex [dict get $Fracture BotLine Groups] $i] \" ,
            }
            set PutStrings [string trimright $PutStrings ,]
            append PutStrings \]
            puts $MyFileVar "            \"groups\":     $PutStrings,"
        }
        puts $MyFileVar "            \"init_point\": [dict get $Fracture BotLine InitPoint],"
        puts $MyFileVar "            \"end_point\":  [dict get $Fracture BotLine EndPoint]"
        puts $MyFileVar "        \},"
        # InterfaceSurface
        puts $MyFileVar "        \"interface_surface\": \{"
        puts $MyFileVar "            \"id\":       [dict get $Fracture InterfaceSurface Id],"
        puts $MyFileVar "            \"layer\":    \"[dict get $Fracture InterfaceSurface Layer]\","
        if {[llength [dict get $Fracture InterfaceSurface Groups]]==0} {
            puts $MyFileVar "            \"groups\":   \[\],"
        } else {
            set PutStrings \[
            for {set i 0} {$i < [llength [dict get $Fracture InterfaceSurface Groups]]} {incr i} {
                append PutStrings \" [lindex [dict get $Fracture InterfaceSurface Groups] $i] \" ,
            }
            set PutStrings [string trimright $PutStrings ,]
            append PutStrings \]
            puts $MyFileVar "            \"groups\":   $PutStrings,"
        }
        puts $MyFileVar "            \"top_line\": \[[lindex [dict get $Fracture InterfaceSurface TopLine] 0],[lindex [dict get $Fracture InterfaceSurface TopLine] 1]\],"
        puts $MyFileVar "            \"bot_line\": \[[lindex [dict get $Fracture InterfaceSurface BotLine] 0],[lindex [dict get $Fracture InterfaceSurface BotLine] 1]\]"
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

#--------------------------------------------------------------------------------------------------------------------------------------------------------------

#~ proc OpenOldCrackTip {FracturesDict MotherFractureId PropagationDict BodySurfaceLines} {


    #~ # Delete InterfaceSurface
    #~ GiD_Process Mescape Geometry Delete Surfaces [dict get $FracturesDict $MotherFractureId InterfaceSurface Id] escape

    #~ WarnWin "Despres de borrar interface surface"

    #~ # Create new poin in TopInitCoordinates location
    #~ GiD_Process Mescape Geometry Create Point \
        #~ [lindex [dict get $PropagationDict TopInitCoordinates] 0] [lindex [dict get $PropagationDict TopInitCoordinates] 1] [lindex [dict get $PropagationDict TopInitCoordinates] 2] escape
    
    #~ WarnWin "Despres de crear nou topinitpoint"
    #~ # Delete old TopLine
    #~ GiD_Process Mescape Geometry Delete Lines [dict get $FracturesDict $MotherFractureId TopLine Id] escape
    #~ set Index [lsearch $BodySurfaceLines [dict get $FracturesDict $MotherFractureId TopLine Id]]
    #~ set BodySurfaceLines [lreplace $BodySurfaceLines $Index $Index]
    #~ # Create new line replacing old TopLine
    #~ GiD_Process Mescape Geometry Create Line Join \
        #~ [GiD_Info Geometry MaxNumPoints] [dict get $FracturesDict $MotherFractureId TopPoint Id] escape escape
    #~ lappend BodySurfaceLines [GiD_Info Geometry MaxNumLines]

    #~ # Create new poin in BotInitCoordinates location
    #~ GiD_Process Mescape Geometry Create Point \
        #~ [lindex [dict get $PropagationDict BotInitCoordinates] 0] [lindex [dict get $PropagationDict BotInitCoordinates] 1] [lindex [dict get $PropagationDict BotInitCoordinates] 2] escape
    #~ # Delete old BotLine
    #~ GiD_Process Mescape Geometry Delete Lines [dict get $FracturesDict $MotherFractureId BotLine Id] escape
    #~ set Index [lsearch $BodySurfaceLines [dict get $FracturesDict $MotherFractureId BotLine Id]]
    #~ set BodySurfaceLines [lreplace $BodySurfaceLines $Index $Index]
    #~ # Create new line replacing old BotLine
    #~ GiD_Process Mescape Geometry Create Line Join \
        #~ [dict get $FracturesDict $MotherFractureId BotPoint Id] [GiD_Info Geometry MaxNumPoints] escape escape
    #~ lappend BodySurfaceLines [GiD_Info Geometry MaxNumLines]

    #~ # Create ContactSurface for the old crack
    #~ set TopLine [expr { [GiD_Info Geometry MaxNumLines]-1 }]
    #~ lappend TopLine [lindex [dict get $FracturesDict $MotherFractureId InterfaceSurface TopLine] 1]
    #~ set BotLine [GiD_Info Geometry MaxNumLines]
    #~ lappend BotLine [lindex [dict get $FracturesDict $MotherFractureId InterfaceSurface BotLine] 1]
    #~ GiD_Geometry create surface [dict get $FracturesDict $MotherFractureId InterfaceSurface Id] contactsurface \
        #~ [dict get $FracturesDict $MotherFractureId InterfaceSurface Layer] 2 $BotLine $TopLine
#~ }

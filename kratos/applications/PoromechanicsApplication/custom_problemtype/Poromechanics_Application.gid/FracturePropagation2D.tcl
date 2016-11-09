proc WriteInitialFracturesData { dir gidpath } {
        
    ## Set Dictionaries
    set FractureId -1
    set BodySurfacesDict [dict create]
    set BodyGroups [GiD_Info conditions Body_Part groups]
    set InterfaceGroups [GiD_Info conditions Interface_Part groups]

    if {[llength $InterfaceGroups] < 1} {
        WarnWin "ERROR: For the moment Fracture Propagation needs at least 1 pre-existing fracture"
    }
    
    for {set i 0} {$i < [llength $InterfaceGroups]} {incr i} {
        if {[lindex [lindex $InterfaceGroups $i] 3] == false} {
            set InterfaceEntities [GiD_EntitiesGroups get [lindex [lindex $InterfaceGroups $i] 1] surfaces]
            for {set j 0} {$j < [llength $InterfaceEntities]} {incr j} {
                set InterfaceSurface [GiD_Geometry get surface [lindex $InterfaceEntities $j]]
                set BotLine [GiD_Geometry get line [lindex [lindex $InterfaceSurface 3] 0]]
                set TopLine [GiD_Geometry get line [lindex [lindex $InterfaceSurface 4] 0]]
                if {([lindex $BotLine 2] == [lindex $TopLine 2]) || ([lindex $BotLine 3] == [lindex $TopLine 2])} {
                    # Define new fracture
                    incr FractureId
                    
                    # Set TipPoint
                    set Point [GiD_Geometry get point [lindex $TopLine 2]]
                    set Groups [GiD_EntitiesGroups entity_groups points [lindex $TopLine 2]]
                    set Coordinates "[lindex $Point 1] [lindex $Point 2] [lindex $Point 3]"
                    dict set FracturesDict $FractureId TipPoint Id [lindex $TopLine 2]
                    dict set FracturesDict $FractureId TipPoint Layer [lindex $Point 0]
                    dict set FracturesDict $FractureId TipPoint Groups $Groups
                    dict set FracturesDict $FractureId TipPoint Coordinates $Coordinates
                    
                    # Set TopPoint
                    set Point [GiD_Geometry get point [lindex $TopLine 3]]
                    set Groups [GiD_EntitiesGroups entity_groups points [lindex $TopLine 3]]
                    set Coordinates "[lindex $Point 1] [lindex $Point 2] [lindex $Point 3]"
                    dict set FracturesDict $FractureId TopPoint Id [lindex $TopLine 3]
                    dict set FracturesDict $FractureId TopPoint Layer [lindex $Point 0]
                    dict set FracturesDict $FractureId TopPoint Groups $Groups
                    dict set FracturesDict $FractureId TopPoint Coordinates $Coordinates
                    
                    # Set BotPoint
                    if {[lindex $BotLine 2] != [lindex $TopLine 2]} {
                        set Point [GiD_Geometry get point [lindex $BotLine 2]]
                        set Groups [GiD_EntitiesGroups entity_groups points [lindex $BotLine 2]]
                        set Coordinates "[lindex $Point 1] [lindex $Point 2] [lindex $Point 3]"
                        dict set FracturesDict $FractureId BotPoint Id [lindex $BotLine 2]
                        
                    } else {
                        set Point [GiD_Geometry get point [lindex $BotLine 3]]
                        set Groups [GiD_EntitiesGroups entity_groups points [lindex $BotLine 3]]
                        set Coordinates "[lindex $Point 1] [lindex $Point 2] [lindex $Point 3]"
                        dict set FracturesDict $FractureId BotPoint Id [lindex $BotLine 3]
                    }
                    dict set FracturesDict $FractureId BotPoint Layer [lindex $Point 0]
                    dict set FracturesDict $FractureId BotPoint Groups $Groups
                    dict set FracturesDict $FractureId BotPoint Coordinates $Coordinates
                    
                    # Set TopLine
                    set Groups [GiD_EntitiesGroups entity_groups lines [lindex [lindex $InterfaceSurface 4] 0]]
                    dict set FracturesDict $FractureId TopLine Id [lindex [lindex $InterfaceSurface 4] 0]
                    dict set FracturesDict $FractureId TopLine Layer [lindex $TopLine 1]
                    dict set FracturesDict $FractureId TopLine Groups $Groups
                    dict set FracturesDict $FractureId TopLine InitPoint [lindex $TopLine 2]
                    dict set FracturesDict $FractureId TopLine EndPoint [lindex $TopLine 3]
                    
                    # Set BotLine
                    set Groups [GiD_EntitiesGroups entity_groups lines [lindex [lindex $InterfaceSurface 3] 0]]
                    dict set FracturesDict $FractureId BotLine Id [lindex [lindex $InterfaceSurface 3] 0]
                    dict set FracturesDict $FractureId BotLine Layer [lindex $BotLine 1]
                    dict set FracturesDict $FractureId BotLine Groups $Groups
                    dict set FracturesDict $FractureId BotLine InitPoint [lindex $BotLine 2]
                    dict set FracturesDict $FractureId BotLine EndPoint [lindex $BotLine 3]
                    
                    # Set InterfaceSurface
                    set Groups [GiD_EntitiesGroups entity_groups surfaces [lindex $InterfaceEntities $j]]
                    dict set FracturesDict $FractureId InterfaceSurface Id [lindex $InterfaceEntities $j]
                    dict set FracturesDict $FractureId InterfaceSurface Layer [lindex $InterfaceSurface 1]
                    dict set FracturesDict $FractureId InterfaceSurface Groups $Groups
                    dict set FracturesDict $FractureId InterfaceSurface TopLine [lindex $InterfaceSurface 4]
                    dict set FracturesDict $FractureId InterfaceSurface BotLine [lindex $InterfaceSurface 3]
                    
                    # Set BodySurface
                    set TipCoordinates [dict get $FracturesDict $FractureId TipPoint Coordinates]
                    set BodySurfaces [list]
                    for {set k 0} {$k < [llength $BodyGroups]} {incr k} {
                        set BodyEntities [GiD_EntitiesGroups get [lindex [lindex $BodyGroups $k] 1] surfaces]
                        for {set l 0} {$l < [llength $BodyEntities]} {incr l} {
                            if {[GiD_Info IsPointInside Surface [lindex $BodyEntities $l] $TipCoordinates] == 1} {
                                lappend BodySurfaces [lindex $BodyEntities $l]
                                #dict set FracturesDict $FractureId BodySurface [lindex $BodyEntities $l]
                                if {[dict exists $BodySurfacesDict [lindex $BodyEntities $l]]==0} {
                                    set BodySurface [GiD_Geometry get surface [lindex $BodyEntities $l]]
                                    set Groups [GiD_EntitiesGroups entity_groups surfaces [lindex $BodyEntities $l]]
                                    dict set BodySurfacesDict [lindex $BodyEntities $l] Layer [lindex $BodySurface 1]
                                    dict set BodySurfacesDict [lindex $BodyEntities $l] Groups $Groups
                                    set Lines [list]
                                    for {set m 0} {$m < [lindex $BodySurface 2]} {incr m} {
                                        lappend Lines [lindex [lindex $BodySurface [expr { 9+$m }]] 0]
                                    }
                                    dict set BodySurfacesDict [lindex $BodyEntities $l] Lines $Lines
                                }
                            }
                        }
                    }
                    dict set FracturesDict $FractureId BodySurfaces $BodySurfaces
                    
                } elseif {([lindex $BotLine 2] == [lindex $TopLine 3]) || ([lindex $BotLine 3] == [lindex $TopLine 3])} {
                    
                    # Define new fracture
                    incr FractureId
                    
                    # Set TipPoint
                    set Point [GiD_Geometry get point [lindex $TopLine 3]]
                    set Groups [GiD_EntitiesGroups entity_groups points [lindex $TopLine 3]]
                    set Coordinates "[lindex $Point 1] [lindex $Point 2] [lindex $Point 3]"
                    dict set FracturesDict $FractureId TipPoint Id [lindex $TopLine 3]
                    dict set FracturesDict $FractureId TipPoint Layer [lindex $Point 0]
                    dict set FracturesDict $FractureId TipPoint Groups $Groups
                    dict set FracturesDict $FractureId TipPoint Coordinates $Coordinates
                    
                    # Set TopPoint
                    set Point [GiD_Geometry get point [lindex $TopLine 2]]
                    set Groups [GiD_EntitiesGroups entity_groups points [lindex $TopLine 2]]
                    set Coordinates "[lindex $Point 1] [lindex $Point 2] [lindex $Point 3]"
                    dict set FracturesDict $FractureId TopPoint Id [lindex $TopLine 2]
                    dict set FracturesDict $FractureId TopPoint Layer [lindex $Point 0]
                    dict set FracturesDict $FractureId TopPoint Groups $Groups
                    dict set FracturesDict $FractureId TopPoint Coordinates $Coordinates
                    
                    # Set BotPoint
                    if {[lindex $BotLine 2] != [lindex $TopLine 3]} {
                        set Point [GiD_Geometry get point [lindex $BotLine 2]]
                        set Groups [GiD_EntitiesGroups entity_groups points [lindex $BotLine 2]]
                        set Coordinates "[lindex $Point 1] [lindex $Point 2] [lindex $Point 3]"
                        dict set FracturesDict $FractureId BotPoint Id [lindex $BotLine 2]
                    } else {
                        set Point [GiD_Geometry get point [lindex $BotLine 3]]
                        set Groups [GiD_EntitiesGroups entity_groups points [lindex $BotLine 3]]
                        set Coordinates "[lindex $Point 1] [lindex $Point 2] [lindex $Point 3]"
                        dict set FracturesDict $FractureId BotPoint Id [lindex $BotLine 3]
                    }
                    dict set FracturesDict $FractureId BotPoint Layer [lindex $Point 0]
                    dict set FracturesDict $FractureId BotPoint Groups $Groups
                    dict set FracturesDict $FractureId BotPoint Coordinates $Coordinates
                    
                    # Set TopLine
                    set Groups [GiD_EntitiesGroups entity_groups lines [lindex [lindex $InterfaceSurface 4] 0]]
                    dict set FracturesDict $FractureId TopLine Id [lindex [lindex $InterfaceSurface 4] 0]
                    dict set FracturesDict $FractureId TopLine Layer [lindex $TopLine 1]
                    dict set FracturesDict $FractureId TopLine Groups $Groups
                    dict set FracturesDict $FractureId TopLine InitPoint [lindex $TopLine 2]
                    dict set FracturesDict $FractureId TopLine EndPoint [lindex $TopLine 3]
                    
                    # Set BotLine
                    set Groups [GiD_EntitiesGroups entity_groups lines [lindex [lindex $InterfaceSurface 3] 0]]
                    dict set FracturesDict $FractureId BotLine Id [lindex [lindex $InterfaceSurface 3] 0]
                    dict set FracturesDict $FractureId BotLine Layer [lindex $BotLine 1]
                    dict set FracturesDict $FractureId BotLine Groups $Groups
                    dict set FracturesDict $FractureId BotLine InitPoint [lindex $BotLine 2]
                    dict set FracturesDict $FractureId BotLine EndPoint [lindex $BotLine 3]
                    
                    # Set InterfaceSurface
                    set Groups [GiD_EntitiesGroups entity_groups surfaces [lindex $InterfaceEntities $j]]
                    dict set FracturesDict $FractureId InterfaceSurface Id [lindex $InterfaceEntities $j]
                    dict set FracturesDict $FractureId InterfaceSurface Layer [lindex $InterfaceSurface 1]
                    dict set FracturesDict $FractureId InterfaceSurface Groups $Groups
                    dict set FracturesDict $FractureId InterfaceSurface TopLine [lindex $InterfaceSurface 4]
                    dict set FracturesDict $FractureId InterfaceSurface BotLine [lindex $InterfaceSurface 3]
                    
                    # Set BodySurface
                    set TipCoordinates [dict get $FracturesDict $FractureId TipPoint Coordinates]
                    set BodySurfaces [list]
                    for {set k 0} {$k < [llength $BodyGroups]} {incr k} {
                        set BodyEntities [GiD_EntitiesGroups get [lindex [lindex $BodyGroups $k] 1] surfaces]
                        for {set l 0} {$l < [llength $BodyEntities]} {incr l} {
                            if {[GiD_Info IsPointInside Surface [lindex $BodyEntities $l] $TipCoordinates] == 1} {
                                lappend BodySurfaces [lindex $BodyEntities $l]
                                #dict set FracturesDict $FractureId BodySurface [lindex $BodyEntities $l]
                                if {[dict exists $BodySurfacesDict [lindex $BodyEntities $l]]==0} {
                                    set BodySurface [GiD_Geometry get surface [lindex $BodyEntities $l]]
                                    set Groups [GiD_EntitiesGroups entity_groups surfaces [lindex $BodyEntities $l]]
                                    dict set BodySurfacesDict [lindex $BodyEntities $l] Layer [lindex $BodySurface 1]
                                    dict set BodySurfacesDict [lindex $BodyEntities $l] Groups $Groups
                                    set Lines [list]
                                    for {set m 0} {$m < [lindex $BodySurface 2]} {incr m} {
                                        lappend Lines [lindex [lindex $BodySurface [expr { 9+$m }]] 0]
                                    }
                                    dict set BodySurfacesDict [lindex $BodyEntities $l] Lines $Lines
                                }
                            }
                        }
                    }
                    dict set FracturesDict $FractureId BodySurfaces $BodySurfaces
                }
            }
        }
    }
        
    ## Start FracturesData.json file
    set filename [file join $dir FracturesData.json]
    set varfile [open $filename w]
        
    puts $varfile "\{"

    ## fracture_data
    puts $varfile "    \"fracture_data\": \{"
    puts $varfile "        \"gid_path\":                             \"$gidpath\","
    puts $varfile "        \"propagation_length\":                   [GiD_AccessValue get gendata Propagation_Length],"
    puts $varfile "        \"propagation_damage\":                   [GiD_AccessValue get gendata Propagation_Damage],"
    puts $varfile "        \"propagation_width\":                    [GiD_AccessValue get gendata Propagation_Width],"
    puts $varfile "        \"propagation_cosangle\":                 [GiD_AccessValue get gendata Propagation_CosAngle],"
    puts $varfile "        \"propagation_frequency\":                [GiD_AccessValue get gendata Propagation_Frequency],"
    # body_domain_sub_model_part_list
    set PutStrings \[
    set Groups [GiD_Info conditions Body_Part groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
        append PutStrings \" [lindex [lindex $Groups $i] 1] \" ,
    }
    set PutStrings [string trimright $PutStrings ,]
    append PutStrings \]
    puts $varfile "        \"body_domain_sub_model_part_list\":      $PutStrings,"
    # interface_domain_sub_model_part_list
    set PutStrings \[
    set Groups [GiD_Info conditions Interface_Part groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
        append PutStrings \" [lindex [lindex $Groups $i] 1] \" ,
    }
    set PutStrings [string trimright $PutStrings ,]
    append PutStrings \]
    puts $varfile "        \"interface_domain_sub_model_part_list\": $PutStrings"
    puts $varfile "    \},"
    
    ## fractures_list
    set iter 0
    puts $varfile "    \"fractures_list\": \[\{"
    dict for {Id Fracture} $FracturesDict {
        incr iter
        puts $varfile "        \"id\":                $Id,"
        # TipPoint
        puts $varfile "        \"tip_point\":         \{"
        puts $varfile "            \"id\":          [dict get $Fracture TipPoint Id],"
        puts $varfile "            \"layer\":       \"[dict get $Fracture TipPoint Layer]\","
        if {[llength [dict get $Fracture TipPoint Groups]]==0} {
            puts $varfile "            \"groups\":      \[\],"
        } else {
            set PutStrings \[
            for {set i 0} {$i < [llength [dict get $Fracture TipPoint Groups]]} {incr i} {
                append PutStrings \" [lindex [dict get $Fracture TipPoint Groups] $i] \" ,
            }
            set PutStrings [string trimright $PutStrings ,]
            append PutStrings \]
            puts $varfile "            \"groups\":      $PutStrings,"
        }
        puts $varfile "            \"coordinates\": \[[lindex [dict get $Fracture TipPoint Coordinates] 0],[lindex [dict get $Fracture TipPoint Coordinates] 1],[lindex [dict get $Fracture TipPoint Coordinates] 2]\]"
        puts $varfile "        \},"
        # TopPoint
        puts $varfile "        \"top_point\":         \{"
        puts $varfile "            \"id\":          [dict get $Fracture TopPoint Id],"
        puts $varfile "            \"layer\":       \"[dict get $Fracture TopPoint Layer]\","
        if {[llength [dict get $Fracture TopPoint Groups]]==0} {
            puts $varfile "            \"groups\":      \[\],"
        } else {
            set PutStrings \[
            for {set i 0} {$i < [llength [dict get $Fracture TopPoint Groups]]} {incr i} {
                append PutStrings \" [lindex [dict get $Fracture TopPoint Groups] $i] \" ,
            }
            set PutStrings [string trimright $PutStrings ,]
            append PutStrings \]
            puts $varfile "            \"groups\":      $PutStrings,"
        }
        puts $varfile "            \"coordinates\": \[[lindex [dict get $Fracture TopPoint Coordinates] 0],[lindex [dict get $Fracture TopPoint Coordinates] 1],[lindex [dict get $Fracture TopPoint Coordinates] 2]\]"
        puts $varfile "        \},"
        # BotPoint
        puts $varfile "        \"bot_point\":         \{"
        puts $varfile "            \"id\":          [dict get $Fracture BotPoint Id],"
        puts $varfile "            \"layer\":       \"[dict get $Fracture BotPoint Layer]\","
        if {[llength [dict get $Fracture BotPoint Groups]]==0} {
            puts $varfile "            \"groups\":      \[\],"
        } else {
            set PutStrings \[
            for {set i 0} {$i < [llength [dict get $Fracture BotPoint Groups]]} {incr i} {
                append PutStrings \" [lindex [dict get $Fracture BotPoint Groups] $i] \" ,
            }
            set PutStrings [string trimright $PutStrings ,]
            append PutStrings \]
            puts $varfile "            \"groups\":      $PutStrings,"
        }
        puts $varfile "            \"coordinates\": \[[lindex [dict get $Fracture BotPoint Coordinates] 0],[lindex [dict get $Fracture BotPoint Coordinates] 1],[lindex [dict get $Fracture BotPoint Coordinates] 2]\]"
        puts $varfile "        \},"
        # TopLine
        puts $varfile "        \"top_line\":          \{"
        puts $varfile "            \"id\":         [dict get $Fracture TopLine Id],"
        puts $varfile "            \"layer\":      \"[dict get $Fracture TopLine Layer]\","
        if {[llength [dict get $Fracture TopLine Groups]]==0} {
            puts $varfile "            \"groups\":     \[\],"
        } else {
            set PutStrings \[
            for {set i 0} {$i < [llength [dict get $Fracture TopLine Groups]]} {incr i} {
                append PutStrings \" [lindex [dict get $Fracture TopLine Groups] $i] \" ,
            }
            set PutStrings [string trimright $PutStrings ,]
            append PutStrings \]
            puts $varfile "            \"groups\":     $PutStrings,"
        }
        puts $varfile "            \"init_point\": [dict get $Fracture TopLine InitPoint],"
        puts $varfile "            \"end_point\":  [dict get $Fracture TopLine EndPoint]"
        puts $varfile "        \},"
        # BotLine
        puts $varfile "        \"bot_line\":          \{"
        puts $varfile "            \"id\":         [dict get $Fracture BotLine Id],"
        puts $varfile "            \"layer\":      \"[dict get $Fracture BotLine Layer]\","
        if {[llength [dict get $Fracture BotLine Groups]]==0} {
            puts $varfile "            \"groups\":     \[\],"
        } else {
            set PutStrings \[
            for {set i 0} {$i < [llength [dict get $Fracture BotLine Groups]]} {incr i} {
                append PutStrings \" [lindex [dict get $Fracture BotLine Groups] $i] \" ,
            }
            set PutStrings [string trimright $PutStrings ,]
            append PutStrings \]
            puts $varfile "            \"groups\":     $PutStrings,"
        }
        puts $varfile "            \"init_point\": [dict get $Fracture BotLine InitPoint],"
        puts $varfile "            \"end_point\":  [dict get $Fracture BotLine EndPoint]"
        puts $varfile "        \},"
        # InterfaceSurface
        puts $varfile "        \"interface_surface\": \{"
        puts $varfile "            \"id\":       [dict get $Fracture InterfaceSurface Id],"
        puts $varfile "            \"layer\":    \"[dict get $Fracture InterfaceSurface Layer]\","
        if {[llength [dict get $Fracture InterfaceSurface Groups]]==0} {
            puts $varfile "            \"groups\":   \[\],"
        } else {
            set PutStrings \[
            for {set i 0} {$i < [llength [dict get $Fracture InterfaceSurface Groups]]} {incr i} {
                append PutStrings \" [lindex [dict get $Fracture InterfaceSurface Groups] $i] \" ,
            }
            set PutStrings [string trimright $PutStrings ,]
            append PutStrings \]
            puts $varfile "            \"groups\":   $PutStrings,"
        }
        puts $varfile "            \"top_line\": \[[lindex [dict get $Fracture InterfaceSurface TopLine] 0],[lindex [dict get $Fracture InterfaceSurface TopLine] 1]\],"
        puts $varfile "            \"bot_line\": \[[lindex [dict get $Fracture InterfaceSurface BotLine] 0],[lindex [dict get $Fracture InterfaceSurface BotLine] 1]\]"
        puts $varfile "        \},"
        # BodySurfaces
        set PutStrings \[
        for {set i 0} {$i < [llength [dict get $Fracture BodySurfaces]]} {incr i} {
            append PutStrings [lindex [dict get $Fracture BodySurfaces] $i] ,
        }
        set PutStrings [string trimright $PutStrings ,]
        append PutStrings \]
        puts $varfile "        \"body_surfaces\":      $PutStrings"
        if {$iter < [dict size $FracturesDict]} {
            puts $varfile "    \},\{"
        } else {
            puts $varfile "    \}\],"
        }
    }
    
    ## body_surfaces_list
    set iter 0
    puts $varfile "    \"body_surfaces_list\": \[\{"
    dict for {Id BodySurface} $BodySurfacesDict {
        incr iter
        puts $varfile "        \"id\":     $Id,"
        puts $varfile "        \"layer\":  \"[dict get $BodySurface Layer]\","
        if {[llength [dict get $BodySurface Groups]]==0} {
            puts $varfile "        \"groups\": \[\],"
        } else {
            set PutStrings \[
            for {set i 0} {$i < [llength [dict get $BodySurface Groups]]} {incr i} {
                append PutStrings \" [lindex [dict get $BodySurface Groups] $i] \" ,
            }
            set PutStrings [string trimright $PutStrings ,]
            append PutStrings \]
            puts $varfile "        \"groups\": $PutStrings,"
        }
        set PutStrings \[
        for {set i 0} {$i < [llength [dict get $BodySurface Lines]]} {incr i} {
            append PutStrings [lindex [dict get $BodySurface Lines] $i] ,
        }
        set PutStrings [string trimright $PutStrings ,]
        append PutStrings \]
        puts $varfile "        \"lines\":  $PutStrings"
        if {$iter < [dict size $BodySurfacesDict]} {
            puts $varfile "    \},\{"
        } else {
            puts $varfile "    \}\]"
        }
    }
    
    puts $varfile "\}"
    
    close $varfile
}

proc GenerateNewFractures { dir PropagationData } {
    
    # Previous Fractures dictionaries
    set FracturesDict [lindex $PropagationData 1]
    set BodySurfacesDict [lindex $PropagationData 2]
    
    # Propagation and Bifurcation dictionaries
    set PropagationDict [lindex $PropagationData 3]
    set BifurcationDict [lindex $PropagationData 4]
    
    # Search TopBodySurface and BotBodySurface for Bifurcations with more than 1 BodySurface
    if {[dict size $BifurcationDict] > 0} {
        dict for {BifId Bifurcation} $BifurcationDict {
            set MotherFractureId [dict get $Bifurcation MotherFractureId]
            set BodySurfaces [dict get $FracturesDict $MotherFractureId BodySurfaces]
            if {[llength $BodySurfaces]>1} {
                set TopTipCoordinates "[lindex [dict get $Bifurcation TopTipCoordinates] 0] [lindex [dict get $Bifurcation TopTipCoordinates] 1] [lindex [dict get $Bifurcation TopTipCoordinates] 2]"
                set BotTipCoordinates "[lindex [dict get $Bifurcation BotTipCoordinates] 0] [lindex [dict get $Bifurcation BotTipCoordinates] 1] [lindex [dict get $Bifurcation BotTipCoordinates] 2]"
                for {set i 0} {$i < [llength $BodySurfaces]} {incr i} {
                    if {[GiD_Info IsPointInside Surface [lindex $BodySurfaces $i] $TopTipCoordinates] == 1} {
                        dict set BifurcationDict $BifId TopBodySurface [lindex $BodySurfaces $i]
                    } elseif {[GiD_Info IsPointInside Surface [lindex $BodySurfaces $i] $BotTipCoordinates] == 1} {
                        dict set BifurcationDict $BifId BotBodySurface [lindex $BodySurfaces $i]
                    }
                }
            }
        }
    }
        
    # Delete all BodySurfaces
    dict for {BodyId BodySurface} $BodySurfacesDict {
        GiD_Process Mescape Geometry Delete Surfaces $BodyId escape
    }
    
    # Loop through all Propagation Fractures
    dict for {PropId Propagation} $PropagationDict {
        set MotherFractureId [dict get $Propagation MotherFractureId]
        
        set BodySurfaceId [lindex [dict get $FracturesDict $MotherFractureId BodySurfaces] 0]
        set BodySurfaceLines [dict get $BodySurfacesDict $BodySurfaceId Lines]
        
        ## Modify Geometry
        # Delete InterfaceSurface
        GiD_Process Mescape Geometry Delete Surfaces [dict get $FracturesDict $MotherFractureId InterfaceSurface Id] escape
        # Move old TipPoint to TopInitCoordinates location
        GiD_Process Mescape Geometry Edit MovePoint [dict get $FracturesDict $MotherFractureId TipPoint Id] \
            [lindex [dict get $Propagation TopInitCoordinates] 0] [lindex [dict get $Propagation TopInitCoordinates] 1] [lindex [dict get $Propagation TopInitCoordinates] 2] escape
        # Create new poin in BotInitCoordinates location
        GiD_Process Mescape Geometry Create Point \
            [lindex [dict get $Propagation BotInitCoordinates] 0] [lindex [dict get $Propagation BotInitCoordinates] 1] [lindex [dict get $Propagation BotInitCoordinates] 2] escape
        # Delete old BotLine
        GiD_Process Mescape Geometry Delete Lines [dict get $FracturesDict $MotherFractureId BotLine Id] escape
        set Index [lsearch $BodySurfaceLines [dict get $FracturesDict $MotherFractureId BotLine Id]]
        set BodySurfaceLines [lreplace $BodySurfaceLines $Index $Index]
        # Create new line replacing old BotLine
        GiD_Process Mescape Geometry Create Line Join \
            [dict get $FracturesDict $MotherFractureId BotPoint Id] [GiD_Info Geometry MaxNumPoints] escape escape
        lappend BodySurfaceLines [GiD_Info Geometry MaxNumLines]
        # Create ContactSurface for the old crack
        set BotLine [GiD_Info Geometry MaxNumLines]
        lappend BotLine [lindex [dict get $FracturesDict $MotherFractureId InterfaceSurface BotLine] 1]
        GiD_Geometry create surface [dict get $FracturesDict $MotherFractureId InterfaceSurface Id] contactsurface \
            [dict get $FracturesDict $MotherFractureId InterfaceSurface Layer] 2 $BotLine [dict get $FracturesDict $MotherFractureId InterfaceSurface TopLine]
        # Create new TipPoint
        GiD_Process Mescape Geometry Create Point \
            [lindex [dict get $Propagation TipCoordinates] 0] [lindex [dict get $Propagation TipCoordinates] 1] [lindex [dict get $Propagation TipCoordinates] 2] escape
        # Create new TopLine
        GiD_Process Mescape Geometry Create Line Join \
            [GiD_Info Geometry MaxNumPoints] [dict get $FracturesDict $MotherFractureId TipPoint Id] escape escape
        lappend BodySurfaceLines [GiD_Info Geometry MaxNumLines]
        # Create new BotLine
        GiD_Process Mescape Geometry Create Line Join \
            [expr { [GiD_Info Geometry MaxNumPoints]-1 }] [GiD_Info Geometry MaxNumPoints] escape escape
        lappend BodySurfaceLines [GiD_Info Geometry MaxNumLines]
        # Create new InterfaceContactSurface
        set TopLine [expr { [GiD_Info Geometry MaxNumLines]-1 }]
        lappend TopLine [lindex [dict get $FracturesDict $MotherFractureId InterfaceSurface TopLine] 1]
        set BotLine [GiD_Info Geometry MaxNumLines]
        lappend BotLine [lindex [dict get $FracturesDict $MotherFractureId InterfaceSurface BotLine] 1]
        GiD_Geometry create surface append contactsurface \
            [dict get $FracturesDict $MotherFractureId InterfaceSurface Layer] 2 $BotLine $TopLine
        
        ## Set Conditions
        for {set i 0} {$i < [llength [dict get $FracturesDict $MotherFractureId InterfaceSurface Groups]]} {incr i} {
            GiD_EntitiesGroups assign [lindex [dict get $FracturesDict $MotherFractureId InterfaceSurface Groups] $i] surfaces [dict get $FracturesDict $MotherFractureId InterfaceSurface Id]
            GiD_EntitiesGroups assign [lindex [dict get $FracturesDict $MotherFractureId InterfaceSurface Groups] $i] surfaces [GiD_Info Geometry MaxNumSurfaces]
        }
        
        ## Set Mesh options
        set TopCoord [dict get $FracturesDict $MotherFractureId TopPoint Coordinates]
        set BotCoord [dict get $FracturesDict $MotherFractureId BotPoint Coordinates]
        set Distance [expr { sqrt(([lindex $TopCoord 0]-[lindex $BotCoord 0])**2 + ([lindex $TopCoord 1]-[lindex $BotCoord 1])**2) }]
        GiD_Process Mescape Meshing AssignSizes Lines [expr { $Distance*1.25 }] \
            [expr { [GiD_Info Geometry MaxNumLines]-2 }] [dict get $FracturesDict $MotherFractureId TopLine Id] escape escape
        GiD_Process Mescape Meshing AssignSizes Lines [expr { [GiD_AccessValue get gendata Propagation_Width]*1.25 }] \
            [expr { [GiD_Info Geometry MaxNumLines]-1 }] [GiD_Info Geometry MaxNumLines] escape escape
        GiD_Process Mescape Meshing AssignSizes Points 0.0 [dict get $FracturesDict $MotherFractureId TipPoint Id] escape escape
        #GiD_Process Mescape Meshing AssignSizes Points [expr { [GiD_AccessValue get gendata Propagation_Width]*0.25 }] [GiD_Info Geometry MaxNumPoints] escape escape
        GiD_Process Mescape Meshing AssignSizes Points [GiD_AccessValue get gendata Propagation_Width] [GiD_Info Geometry MaxNumPoints] escape escape
        
        ## Replace $FracturesDict $MotherFractureId by the new fracture
        dict set FracturesDict $MotherFractureId TopPoint Id [dict get $FracturesDict $MotherFractureId TipPoint Id]
        dict set FracturesDict $MotherFractureId TopPoint Coordinates [dict get $Propagation TopInitCoordinates]
        dict set FracturesDict $MotherFractureId BotPoint Id [expr { [GiD_Info Geometry MaxNumPoints]-1 }]
        dict set FracturesDict $MotherFractureId BotPoint Coordinates [dict get $Propagation BotInitCoordinates]
        dict set FracturesDict $MotherFractureId TipPoint Id [GiD_Info Geometry MaxNumPoints]
        dict set FracturesDict $MotherFractureId TipPoint Coordinates [dict get $Propagation TipCoordinates]
        dict set FracturesDict $MotherFractureId TopLine Id [expr { [GiD_Info Geometry MaxNumLines]-1 }]
        dict set FracturesDict $MotherFractureId TopLine InitPoint [dict get $FracturesDict $MotherFractureId TipPoint Id]
        dict set FracturesDict $MotherFractureId TopLine EndPoint [dict get $FracturesDict $MotherFractureId TopPoint Id]
        dict set FracturesDict $MotherFractureId BotLine Id [GiD_Info Geometry MaxNumLines]
        dict set FracturesDict $MotherFractureId BotLine InitPoint [dict get $FracturesDict $MotherFractureId BotPoint Id]
        dict set FracturesDict $MotherFractureId BotLine EndPoint [dict get $FracturesDict $MotherFractureId TipPoint Id]
        dict set FracturesDict $MotherFractureId InterfaceSurface Id [GiD_Info Geometry MaxNumSurfaces]
        dict set FracturesDict $MotherFractureId InterfaceSurface TopLine $TopLine
        dict set FracturesDict $MotherFractureId InterfaceSurface BotLine $BotLine
        
        dict set BodySurfacesDict $BodySurfaceId Lines $BodySurfaceLines
    }
    
    # interface_domain_sub_model_part_old_list
    set InterfaceGroupsOld \[
    set Groups [GiD_Info conditions Interface_Part groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
        append InterfaceGroupsOld \" [lindex [lindex $Groups $i] 1] \" ,
    }
    set InterfaceGroupsOld [string trimright $InterfaceGroupsOld ,]
    append InterfaceGroupsOld \]

    set LinkInterfaceGroup 0
    if {[dict size $BifurcationDict] > 0} {
        set Groups [GiD_Info conditions Interface_Part groups]
        for {set i 0} {$i < [llength $Groups]} {incr i} {
            if {[lindex [lindex $Groups $i] 3]==true} {
                set LinkInterfaceGroup [lindex [lindex $Groups $i] 1]
                break
            }
        }
        if {$LinkInterfaceGroup == 0} {
            GiD_Groups create Bifurcation_Link_Part_9
            set LinkInterfaceGroup Bifurcation_Link_Part_9
            set ConditionValues "true [lindex [lindex $Groups 0] 4] [lindex [lindex $Groups 0] 5] [lindex [lindex $Groups 0] 6] [lindex [lindex $Groups 0] 7] \
            [lindex [lindex $Groups 0] 8] [lindex [lindex $Groups 0] 9] [lindex [lindex $Groups 0] 10] [lindex [lindex $Groups 0] 11] [lindex [lindex $Groups 0] 12]\
            [lindex [lindex $Groups 0] 13] [lindex [lindex $Groups 0] 14] [lindex [lindex $Groups 0] 15] [lindex [lindex $Groups 0] 16] 0.0\
            [lindex [lindex $Groups 0] 18] [lindex [lindex $Groups 0] 19]"
            GiD_AssignData condition Interface_Part groups $ConditionValues $LinkInterfaceGroup
        }
    }
    
    # Loop through all Bifurcation Fractures
    dict for {BifId Bifurcation} $BifurcationDict {
        set MotherFractureId [dict get $Bifurcation MotherFractureId]
        set BodySurfaces [dict get $FracturesDict $MotherFractureId BodySurfaces]
        
        if {[llength $BodySurfaces]==1} {
            set BodySurfaceId [lindex $BodySurfaces 0]
            set BodySurfaceLines [dict get $BodySurfacesDict $BodySurfaceId Lines]
            
            ## Modify Geometry
            # Delete InterfaceSurface
            GiD_Process Mescape Geometry Delete Surfaces [dict get $FracturesDict $MotherFractureId InterfaceSurface Id] escape
            # Create new poin in TopInitCoordinates location
            GiD_Process Mescape Geometry Create Point \
                [lindex [dict get $Bifurcation TopInitCoordinates] 0] [lindex [dict get $Bifurcation TopInitCoordinates] 1] [lindex [dict get $Bifurcation TopInitCoordinates] 2] escape
            # Delete old TopLine
            GiD_Process Mescape Geometry Delete Lines [dict get $FracturesDict $MotherFractureId TopLine Id] escape
            set Index [lsearch $BodySurfaceLines [dict get $FracturesDict $MotherFractureId TopLine Id]]
            set BodySurfaceLines [lreplace $BodySurfaceLines $Index $Index]
            # Create new line replacing old TopLine
            GiD_Process Mescape Geometry Create Line Join \
                [GiD_Info Geometry MaxNumPoints] [dict get $FracturesDict $MotherFractureId TopPoint Id] escape escape
            lappend BodySurfaceLines [GiD_Info Geometry MaxNumLines]
            # Create new poin in BotInitCoordinates location
            GiD_Process Mescape Geometry Create Point \
                [lindex [dict get $Bifurcation BotInitCoordinates] 0] [lindex [dict get $Bifurcation BotInitCoordinates] 1] [lindex [dict get $Bifurcation BotInitCoordinates] 2] escape
            # Delete old BotLine
            GiD_Process Mescape Geometry Delete Lines [dict get $FracturesDict $MotherFractureId BotLine Id] escape
            set Index [lsearch $BodySurfaceLines [dict get $FracturesDict $MotherFractureId BotLine Id]]
            set BodySurfaceLines [lreplace $BodySurfaceLines $Index $Index]
            # Create new line replacing old BotLine
            GiD_Process Mescape Geometry Create Line Join \
                [dict get $FracturesDict $MotherFractureId BotPoint Id] [GiD_Info Geometry MaxNumPoints] escape escape
            lappend BodySurfaceLines [GiD_Info Geometry MaxNumLines]
            # Create ContactSurface for the old crack
            set TopLine [expr { [GiD_Info Geometry MaxNumLines]-1 }]
            lappend TopLine [lindex [dict get $FracturesDict $MotherFractureId InterfaceSurface TopLine] 1]
            set BotLine [GiD_Info Geometry MaxNumLines]
            lappend BotLine [lindex [dict get $FracturesDict $MotherFractureId InterfaceSurface BotLine] 1]
            GiD_Geometry create surface [dict get $FracturesDict $MotherFractureId InterfaceSurface Id] contactsurface \
                [dict get $FracturesDict $MotherFractureId InterfaceSurface Layer] 2 $BotLine $TopLine
            ## New TopTip
            # Create new TipPoint
            GiD_Process Mescape Geometry Create Point \
                [lindex [dict get $Bifurcation TopTipCoordinates] 0] [lindex [dict get $Bifurcation TopTipCoordinates] 1] [lindex [dict get $Bifurcation TopTipCoordinates] 2] escape
            # Create new TopLine
            GiD_Process Mescape Geometry Create Line Join \
                [GiD_Info Geometry MaxNumPoints] [expr { [GiD_Info Geometry MaxNumPoints]-2 }] escape escape
            lappend BodySurfaceLines [GiD_Info Geometry MaxNumLines]
            # Create new BotLine
            GiD_Process Mescape Geometry Create Line Join \
                [dict get $FracturesDict $MotherFractureId TipPoint Id] [GiD_Info Geometry MaxNumPoints] escape escape
            lappend BodySurfaceLines [GiD_Info Geometry MaxNumLines]
            # Create new InterfaceContactSurface
            set TopTopLine [expr { [GiD_Info Geometry MaxNumLines]-1 }]
            lappend TopTopLine [lindex [dict get $FracturesDict $MotherFractureId InterfaceSurface TopLine] 1]
            set TopBotLine [GiD_Info Geometry MaxNumLines]
            lappend TopBotLine [lindex [dict get $FracturesDict $MotherFractureId InterfaceSurface BotLine] 1]
            GiD_Geometry create surface append contactsurface \
                [dict get $FracturesDict $MotherFractureId InterfaceSurface Layer] 2 $TopBotLine $TopTopLine
            ## New BotTip
            # Create new TipPoint
            GiD_Process Mescape Geometry Create Point \
                [lindex [dict get $Bifurcation BotTipCoordinates] 0] [lindex [dict get $Bifurcation BotTipCoordinates] 1] [lindex [dict get $Bifurcation BotTipCoordinates] 2] escape
            # Create new TopLine
            GiD_Process Mescape Geometry Create Line Join \
                [GiD_Info Geometry MaxNumPoints] [dict get $FracturesDict $MotherFractureId TipPoint Id] escape escape
            lappend BodySurfaceLines [GiD_Info Geometry MaxNumLines]
            # Create new BotLine
            GiD_Process Mescape Geometry Create Line Join \
                [expr { [GiD_Info Geometry MaxNumPoints]-2 }] [GiD_Info Geometry MaxNumPoints] escape escape
            lappend BodySurfaceLines [GiD_Info Geometry MaxNumLines]
            # Create new InterfaceContactSurface
            set BotTopLine [expr { [GiD_Info Geometry MaxNumLines]-1 }]
            lappend BotTopLine [lindex [dict get $FracturesDict $MotherFractureId InterfaceSurface TopLine] 1]
            set BotBotLine [GiD_Info Geometry MaxNumLines]
            lappend BotBotLine [lindex [dict get $FracturesDict $MotherFractureId InterfaceSurface BotLine] 1]
            GiD_Geometry create surface append contactsurface \
                [dict get $FracturesDict $MotherFractureId InterfaceSurface Layer] 2 $BotBotLine $BotTopLine
            ## New Link Surface
            GiD_Process Mescape Geometry Create Line Join \
                [dict get $FracturesDict $MotherFractureId TipPoint Id] [expr { [GiD_Info Geometry MaxNumPoints]-3 }] \
                [expr { [GiD_Info Geometry MaxNumPoints]-2 }] [dict get $FracturesDict $MotherFractureId TipPoint Id] escape escape
            GiD_Process Mescape Geometry Create NurbsSurface [expr { [GiD_Info Geometry MaxNumLines]-2 }] \
                [expr { [GiD_Info Geometry MaxNumLines]-1 }] [GiD_Info Geometry MaxNumLines] escape escape
            
            ## Set Conditions
            for {set i 0} {$i < [llength [dict get $FracturesDict $MotherFractureId InterfaceSurface Groups]]} {incr i} {
                GiD_EntitiesGroups assign [lindex [dict get $FracturesDict $MotherFractureId InterfaceSurface Groups] $i] surfaces [dict get $FracturesDict $MotherFractureId InterfaceSurface Id]
                GiD_EntitiesGroups assign [lindex [dict get $FracturesDict $MotherFractureId InterfaceSurface Groups] $i] surfaces [expr { [GiD_Info Geometry MaxNumSurfaces]-1 }]
                GiD_EntitiesGroups assign [lindex [dict get $FracturesDict $MotherFractureId InterfaceSurface Groups] $i] surfaces [expr { [GiD_Info Geometry MaxNumSurfaces]-2 }]
            }
            GiD_EntitiesGroups assign $LinkInterfaceGroup surfaces [GiD_Info Geometry MaxNumSurfaces]
            
            ## Set Mesh options
            set TopCoord [dict get $FracturesDict $MotherFractureId TopPoint Coordinates]
            set BotCoord [dict get $FracturesDict $MotherFractureId BotPoint Coordinates]
            set Distance [expr { sqrt(([lindex $TopCoord 0]-[lindex $BotCoord 0])**2 + ([lindex $TopCoord 1]-[lindex $BotCoord 1])**2) }]
            GiD_Process Mescape Meshing AssignSizes Lines [expr { $Distance*1.25 }] \
                [expr { [GiD_Info Geometry MaxNumLines]-7 }] [expr { [GiD_Info Geometry MaxNumLines]-8 }] escape escape
            GiD_Process Mescape Meshing AssignSizes Lines [expr { [GiD_AccessValue get gendata Propagation_Width]*1.25 }] \
                [expr { [GiD_Info Geometry MaxNumLines]-3 }] [expr { [GiD_Info Geometry MaxNumLines]-4 }] \
                [expr { [GiD_Info Geometry MaxNumLines]-5 }] [expr { [GiD_Info Geometry MaxNumLines]-6 }] escape escape
            GiD_Process Mescape Meshing AssignSizes Points 0.0 [dict get $FracturesDict $MotherFractureId TipPoint Id] escape escape
            #GiD_Process Mescape Meshing AssignSizes Points [expr { [GiD_AccessValue get gendata Propagation_Width]*0.25 }] \
            #    [expr { [GiD_Info Geometry MaxNumPoints]-1 }] [GiD_Info Geometry MaxNumPoints] escape escape
            GiD_Process Mescape Meshing AssignSizes Points [GiD_AccessValue get gendata Propagation_Width] \
                [expr { [GiD_Info Geometry MaxNumPoints]-1 }] [GiD_Info Geometry MaxNumPoints] escape escape
            GiD_Process Mescape Meshing Structured Surfaces [GiD_Info Geometry MaxNumSurfaces] escape 1 [GiD_Info Geometry MaxNumLines] escape escape
            
            ## Replace $FracturesDict $MotherFractureId by the new TopTip fracture
            dict set FracturesDict $MotherFractureId TopPoint Id [expr { [GiD_Info Geometry MaxNumPoints]-3 }]
            dict set FracturesDict $MotherFractureId TopPoint Coordinates [dict get $Bifurcation TopInitCoordinates]
            dict set FracturesDict $MotherFractureId BotPoint Id [dict get $FracturesDict $MotherFractureId TipPoint Id]
            dict set FracturesDict $MotherFractureId BotPoint Coordinates [dict get $FracturesDict $MotherFractureId TipPoint Coordinates]
            dict set FracturesDict $MotherFractureId TipPoint Id [expr { [GiD_Info Geometry MaxNumPoints]-1 }]
            dict set FracturesDict $MotherFractureId TipPoint Coordinates [dict get $Bifurcation TopTipCoordinates]
            dict set FracturesDict $MotherFractureId TopLine Id [expr { [GiD_Info Geometry MaxNumLines]-6 }]
            dict set FracturesDict $MotherFractureId TopLine InitPoint [dict get $FracturesDict $MotherFractureId TipPoint Id]
            dict set FracturesDict $MotherFractureId TopLine EndPoint [dict get $FracturesDict $MotherFractureId TopPoint Id]
            dict set FracturesDict $MotherFractureId BotLine Id [expr { [GiD_Info Geometry MaxNumLines]-5 }]
            dict set FracturesDict $MotherFractureId BotLine InitPoint [dict get $FracturesDict $MotherFractureId BotPoint Id]
            dict set FracturesDict $MotherFractureId BotLine EndPoint [dict get $FracturesDict $MotherFractureId TipPoint Id]
            dict set FracturesDict $MotherFractureId InterfaceSurface Id [expr { [GiD_Info Geometry MaxNumSurfaces]-2 }]
            dict set FracturesDict $MotherFractureId InterfaceSurface TopLine $TopTopLine
            dict set FracturesDict $MotherFractureId InterfaceSurface BotLine $TopBotLine
            ## Add an additional value in $FracturesDict for the new BotTip fracture
            set NewFractureId [dict size $FracturesDict]
            dict set FracturesDict $NewFractureId TipPoint Id [GiD_Info Geometry MaxNumPoints]
            dict set FracturesDict $NewFractureId TipPoint Layer [dict get $FracturesDict $MotherFractureId TipPoint Layer]
            dict set FracturesDict $NewFractureId TipPoint Groups [dict get $FracturesDict $MotherFractureId TipPoint Groups]
            dict set FracturesDict $NewFractureId TipPoint Coordinates [dict get $Bifurcation BotTipCoordinates]
            dict set FracturesDict $NewFractureId TopPoint Id [dict get $FracturesDict $MotherFractureId BotPoint Id]
            dict set FracturesDict $NewFractureId TopPoint Layer [dict get $FracturesDict $MotherFractureId TopPoint Layer]
            dict set FracturesDict $NewFractureId TopPoint Groups [dict get $FracturesDict $MotherFractureId TopPoint Groups]
            dict set FracturesDict $NewFractureId TopPoint Coordinates [dict get $FracturesDict $MotherFractureId BotPoint Coordinates]
            dict set FracturesDict $NewFractureId BotPoint Id [expr { [GiD_Info Geometry MaxNumPoints]-2 }]
            dict set FracturesDict $NewFractureId BotPoint Layer [dict get $FracturesDict $MotherFractureId BotPoint Layer]
            dict set FracturesDict $NewFractureId BotPoint Groups [dict get $FracturesDict $MotherFractureId BotPoint Groups]
            dict set FracturesDict $NewFractureId BotPoint Coordinates [dict get $Bifurcation BotInitCoordinates]
            dict set FracturesDict $NewFractureId TopLine Id [expr { [GiD_Info Geometry MaxNumLines]-4 }]
            dict set FracturesDict $NewFractureId TopLine Layer [dict get $FracturesDict $MotherFractureId TopLine Layer]
            dict set FracturesDict $NewFractureId TopLine Groups [dict get $FracturesDict $MotherFractureId TopLine Groups]
            dict set FracturesDict $NewFractureId TopLine InitPoint [dict get $FracturesDict $NewFractureId TipPoint Id]
            dict set FracturesDict $NewFractureId TopLine EndPoint [dict get $FracturesDict $NewFractureId TopPoint Id]
            dict set FracturesDict $NewFractureId BotLine Id [expr { [GiD_Info Geometry MaxNumLines]-3 }]
            dict set FracturesDict $NewFractureId BotLine Layer [dict get $FracturesDict $MotherFractureId BotLine Layer]
            dict set FracturesDict $NewFractureId BotLine Groups [dict get $FracturesDict $MotherFractureId BotLine Groups]
            dict set FracturesDict $NewFractureId BotLine InitPoint [dict get $FracturesDict $NewFractureId BotPoint Id]
            dict set FracturesDict $NewFractureId BotLine EndPoint [dict get $FracturesDict $NewFractureId TipPoint Id]
            dict set FracturesDict $NewFractureId InterfaceSurface Id [expr { [GiD_Info Geometry MaxNumSurfaces]-1 }]
            dict set FracturesDict $NewFractureId InterfaceSurface Layer [dict get $FracturesDict $MotherFractureId InterfaceSurface Layer]
            dict set FracturesDict $NewFractureId InterfaceSurface Groups [dict get $FracturesDict $MotherFractureId InterfaceSurface Groups]
            dict set FracturesDict $NewFractureId InterfaceSurface TopLine $BotTopLine
            dict set FracturesDict $NewFractureId InterfaceSurface BotLine $BotBotLine
            
            dict set BodySurfacesDict $BodySurfaceId Lines $BodySurfaceLines
        } else {
            set TopBodySurfaceId [dict get $Bifurcation TopBodySurface]
            set TopBodySurfaceLines [dict get $BodySurfacesDict $TopBodySurfaceId Lines]
            set BotBodySurfaceId [dict get $Bifurcation BotBodySurface]
            set BotBodySurfaceLines [dict get $BodySurfacesDict $BotBodySurfaceId Lines]
            
            ## Modify Geometry
            # Delete InterfaceSurface
            GiD_Process Mescape Geometry Delete Surfaces [dict get $FracturesDict $MotherFractureId InterfaceSurface Id] escape
            # Create new poin in TopInitCoordinates location
            GiD_Process Mescape Geometry Create Point \
                [lindex [dict get $Bifurcation TopInitCoordinates] 0] [lindex [dict get $Bifurcation TopInitCoordinates] 1] [lindex [dict get $Bifurcation TopInitCoordinates] 2] escape
            # Delete old TopLine
            GiD_Process Mescape Geometry Delete Lines [dict get $FracturesDict $MotherFractureId TopLine Id] escape
            set Index [lsearch $TopBodySurfaceLines [dict get $FracturesDict $MotherFractureId TopLine Id]]
            set TopBodySurfaceLines [lreplace $TopBodySurfaceLines $Index $Index]
            # Create new line replacing old TopLine
            GiD_Process Mescape Geometry Create Line Join \
                [GiD_Info Geometry MaxNumPoints] [dict get $FracturesDict $MotherFractureId TopPoint Id] escape escape
            lappend TopBodySurfaceLines [GiD_Info Geometry MaxNumLines]
            # Create new poin in BotInitCoordinates location
            GiD_Process Mescape Geometry Create Point \
                [lindex [dict get $Bifurcation BotInitCoordinates] 0] [lindex [dict get $Bifurcation BotInitCoordinates] 1] [lindex [dict get $Bifurcation BotInitCoordinates] 2] escape
            # Delete old BotLine
            GiD_Process Mescape Geometry Delete Lines [dict get $FracturesDict $MotherFractureId BotLine Id] escape
            set Index [lsearch $BotBodySurfaceLines [dict get $FracturesDict $MotherFractureId BotLine Id]]
            set BotBodySurfaceLines [lreplace $BotBodySurfaceLines $Index $Index]
            # Create new line replacing old BotLine
            GiD_Process Mescape Geometry Create Line Join \
                [dict get $FracturesDict $MotherFractureId BotPoint Id] [GiD_Info Geometry MaxNumPoints] escape escape
            lappend BotBodySurfaceLines [GiD_Info Geometry MaxNumLines]
            # Create ContactSurface for the old crack
            set TopLine [expr { [GiD_Info Geometry MaxNumLines]-1 }]
            lappend TopLine [lindex [dict get $FracturesDict $MotherFractureId InterfaceSurface TopLine] 1]
            set BotLine [GiD_Info Geometry MaxNumLines]
            lappend BotLine [lindex [dict get $FracturesDict $MotherFractureId InterfaceSurface BotLine] 1]
            GiD_Geometry create surface [dict get $FracturesDict $MotherFractureId InterfaceSurface Id] contactsurface \
                [dict get $FracturesDict $MotherFractureId InterfaceSurface Layer] 2 $BotLine $TopLine
            ## New TopTip
            # Create new TipPoint
            GiD_Process Mescape Geometry Create Point \
                [lindex [dict get $Bifurcation TopTipCoordinates] 0] [lindex [dict get $Bifurcation TopTipCoordinates] 1] [lindex [dict get $Bifurcation TopTipCoordinates] 2] escape
            # Create new TopLine
            GiD_Process Mescape Geometry Create Line Join \
                [GiD_Info Geometry MaxNumPoints] [expr { [GiD_Info Geometry MaxNumPoints]-2 }] escape escape
            lappend TopBodySurfaceLines [GiD_Info Geometry MaxNumLines]
            # Create new BotLine
            GiD_Process Mescape Geometry Create Line Join \
                [dict get $FracturesDict $MotherFractureId TipPoint Id] [GiD_Info Geometry MaxNumPoints] escape escape
            lappend TopBodySurfaceLines [GiD_Info Geometry MaxNumLines]
            # Create new InterfaceContactSurface
            set TopTopLine [expr { [GiD_Info Geometry MaxNumLines]-1 }]
            lappend TopTopLine [lindex [dict get $FracturesDict $MotherFractureId InterfaceSurface TopLine] 1]
            set TopBotLine [GiD_Info Geometry MaxNumLines]
            lappend TopBotLine [lindex [dict get $FracturesDict $MotherFractureId InterfaceSurface BotLine] 1]
            GiD_Geometry create surface append contactsurface \
                [dict get $FracturesDict $MotherFractureId InterfaceSurface Layer] 2 $TopBotLine $TopTopLine
            ## New BotTip
            # Create new TipPoint
            GiD_Process Mescape Geometry Create Point \
                [lindex [dict get $Bifurcation BotTipCoordinates] 0] [lindex [dict get $Bifurcation BotTipCoordinates] 1] [lindex [dict get $Bifurcation BotTipCoordinates] 2] escape
            # Create new TopLine
            GiD_Process Mescape Geometry Create Line Join \
                [GiD_Info Geometry MaxNumPoints] [dict get $FracturesDict $MotherFractureId TipPoint Id] escape escape
            lappend BotBodySurfaceLines [GiD_Info Geometry MaxNumLines]
            # Create new BotLine
            GiD_Process Mescape Geometry Create Line Join \
                [expr { [GiD_Info Geometry MaxNumPoints]-2 }] [GiD_Info Geometry MaxNumPoints] escape escape
            lappend BotBodySurfaceLines [GiD_Info Geometry MaxNumLines]
            # Create new InterfaceContactSurface
            set BotTopLine [expr { [GiD_Info Geometry MaxNumLines]-1 }]
            lappend BotTopLine [lindex [dict get $FracturesDict $MotherFractureId InterfaceSurface TopLine] 1]
            set BotBotLine [GiD_Info Geometry MaxNumLines]
            lappend BotBotLine [lindex [dict get $FracturesDict $MotherFractureId InterfaceSurface BotLine] 1]
            GiD_Geometry create surface append contactsurface \
                [dict get $FracturesDict $MotherFractureId InterfaceSurface Layer] 2 $BotBotLine $BotTopLine
            ## New Link Surface
            GiD_Process Mescape Geometry Create Line Join \
                [dict get $FracturesDict $MotherFractureId TipPoint Id] [expr { [GiD_Info Geometry MaxNumPoints]-3 }] \
                [expr { [GiD_Info Geometry MaxNumPoints]-2 }] [dict get $FracturesDict $MotherFractureId TipPoint Id] escape escape
            GiD_Process Mescape Geometry Create NurbsSurface [expr { [GiD_Info Geometry MaxNumLines]-2 }] \
                [expr { [GiD_Info Geometry MaxNumLines]-1 }] [GiD_Info Geometry MaxNumLines] escape escape
            
            ## Set Conditions
            for {set i 0} {$i < [llength [dict get $FracturesDict $MotherFractureId InterfaceSurface Groups]]} {incr i} {
                GiD_EntitiesGroups assign [lindex [dict get $FracturesDict $MotherFractureId InterfaceSurface Groups] $i] surfaces [dict get $FracturesDict $MotherFractureId InterfaceSurface Id]
                GiD_EntitiesGroups assign [lindex [dict get $FracturesDict $MotherFractureId InterfaceSurface Groups] $i] surfaces [expr { [GiD_Info Geometry MaxNumSurfaces]-1 }]
                GiD_EntitiesGroups assign [lindex [dict get $FracturesDict $MotherFractureId InterfaceSurface Groups] $i] surfaces [expr { [GiD_Info Geometry MaxNumSurfaces]-2 }]
            }
            GiD_EntitiesGroups assign $LinkInterfaceGroup surfaces [GiD_Info Geometry MaxNumSurfaces]
            
            ## Set Mesh options
            set TopCoord [dict get $FracturesDict $MotherFractureId TopPoint Coordinates]
            set BotCoord [dict get $FracturesDict $MotherFractureId BotPoint Coordinates]
            set Distance [expr { sqrt(([lindex $TopCoord 0]-[lindex $BotCoord 0])**2 + ([lindex $TopCoord 1]-[lindex $BotCoord 1])**2) }]
            GiD_Process Mescape Meshing AssignSizes Lines [expr { $Distance*1.25 }] \
                [expr { [GiD_Info Geometry MaxNumLines]-7 }] [expr { [GiD_Info Geometry MaxNumLines]-8 }] escape escape
            GiD_Process Mescape Meshing AssignSizes Lines [expr { [GiD_AccessValue get gendata Propagation_Width]*1.25 }] \
                [expr { [GiD_Info Geometry MaxNumLines]-3 }] [expr { [GiD_Info Geometry MaxNumLines]-4 }] \
                [expr { [GiD_Info Geometry MaxNumLines]-5 }] [expr { [GiD_Info Geometry MaxNumLines]-6 }] escape escape
            GiD_Process Mescape Meshing AssignSizes Points 0.0 [dict get $FracturesDict $MotherFractureId TipPoint Id] escape escape
            #GiD_Process Mescape Meshing AssignSizes Points [expr { [GiD_AccessValue get gendata Propagation_Width]*0.25 }] \
            #    [expr { [GiD_Info Geometry MaxNumPoints]-1 }] [GiD_Info Geometry MaxNumPoints] escape escape
            GiD_Process Mescape Meshing AssignSizes Points [GiD_AccessValue get gendata Propagation_Width] \
                [expr { [GiD_Info Geometry MaxNumPoints]-1 }] [GiD_Info Geometry MaxNumPoints] escape escape
            GiD_Process Mescape Meshing Structured Surfaces [GiD_Info Geometry MaxNumSurfaces] escape 1 [GiD_Info Geometry MaxNumLines] escape escape
            
            ## Replace $FracturesDict $MotherFractureId by the new TopTip fracture
            dict set FracturesDict $MotherFractureId TopPoint Id [expr { [GiD_Info Geometry MaxNumPoints]-3 }]
            dict set FracturesDict $MotherFractureId TopPoint Coordinates [dict get $Bifurcation TopInitCoordinates]
            dict set FracturesDict $MotherFractureId BotPoint Id [dict get $FracturesDict $MotherFractureId TipPoint Id]
            dict set FracturesDict $MotherFractureId BotPoint Coordinates [dict get $FracturesDict $MotherFractureId TipPoint Coordinates]
            dict set FracturesDict $MotherFractureId TipPoint Id [expr { [GiD_Info Geometry MaxNumPoints]-1 }]
            dict set FracturesDict $MotherFractureId TipPoint Coordinates [dict get $Bifurcation TopTipCoordinates]
            dict set FracturesDict $MotherFractureId TopLine Id [expr { [GiD_Info Geometry MaxNumLines]-6 }]
            dict set FracturesDict $MotherFractureId TopLine InitPoint [dict get $FracturesDict $MotherFractureId TipPoint Id]
            dict set FracturesDict $MotherFractureId TopLine EndPoint [dict get $FracturesDict $MotherFractureId TopPoint Id]
            dict set FracturesDict $MotherFractureId BotLine Id [expr { [GiD_Info Geometry MaxNumLines]-5 }]
            dict set FracturesDict $MotherFractureId BotLine InitPoint [dict get $FracturesDict $MotherFractureId BotPoint Id]
            dict set FracturesDict $MotherFractureId BotLine EndPoint [dict get $FracturesDict $MotherFractureId TipPoint Id]
            dict set FracturesDict $MotherFractureId InterfaceSurface Id [expr { [GiD_Info Geometry MaxNumSurfaces]-2 }]
            dict set FracturesDict $MotherFractureId InterfaceSurface TopLine $TopTopLine
            dict set FracturesDict $MotherFractureId InterfaceSurface BotLine $TopBotLine
            ## Add an additional value in $FracturesDict for the new BotTip fracture
            set NewFractureId [dict size $FracturesDict]
            dict set FracturesDict $NewFractureId TipPoint Id [GiD_Info Geometry MaxNumPoints]
            dict set FracturesDict $NewFractureId TipPoint Layer [dict get $FracturesDict $MotherFractureId TipPoint Layer]
            dict set FracturesDict $NewFractureId TipPoint Groups [dict get $FracturesDict $MotherFractureId TipPoint Groups]
            dict set FracturesDict $NewFractureId TipPoint Coordinates [dict get $Bifurcation BotTipCoordinates]
            dict set FracturesDict $NewFractureId TopPoint Id [dict get $FracturesDict $MotherFractureId BotPoint Id]
            dict set FracturesDict $NewFractureId TopPoint Layer [dict get $FracturesDict $MotherFractureId TopPoint Layer]
            dict set FracturesDict $NewFractureId TopPoint Groups [dict get $FracturesDict $MotherFractureId TopPoint Groups]
            dict set FracturesDict $NewFractureId TopPoint Coordinates [dict get $FracturesDict $MotherFractureId BotPoint Coordinates]
            dict set FracturesDict $NewFractureId BotPoint Id [expr { [GiD_Info Geometry MaxNumPoints]-2 }]
            dict set FracturesDict $NewFractureId BotPoint Layer [dict get $FracturesDict $MotherFractureId BotPoint Layer]
            dict set FracturesDict $NewFractureId BotPoint Groups [dict get $FracturesDict $MotherFractureId BotPoint Groups]
            dict set FracturesDict $NewFractureId BotPoint Coordinates [dict get $Bifurcation BotInitCoordinates]
            dict set FracturesDict $NewFractureId TopLine Id [expr { [GiD_Info Geometry MaxNumLines]-4 }]
            dict set FracturesDict $NewFractureId TopLine Layer [dict get $FracturesDict $MotherFractureId TopLine Layer]
            dict set FracturesDict $NewFractureId TopLine Groups [dict get $FracturesDict $MotherFractureId TopLine Groups]
            dict set FracturesDict $NewFractureId TopLine InitPoint [dict get $FracturesDict $NewFractureId TipPoint Id]
            dict set FracturesDict $NewFractureId TopLine EndPoint [dict get $FracturesDict $NewFractureId TopPoint Id]
            dict set FracturesDict $NewFractureId BotLine Id [expr { [GiD_Info Geometry MaxNumLines]-3 }]
            dict set FracturesDict $NewFractureId BotLine Layer [dict get $FracturesDict $MotherFractureId BotLine Layer]
            dict set FracturesDict $NewFractureId BotLine Groups [dict get $FracturesDict $MotherFractureId BotLine Groups]
            dict set FracturesDict $NewFractureId BotLine InitPoint [dict get $FracturesDict $NewFractureId BotPoint Id]
            dict set FracturesDict $NewFractureId BotLine EndPoint [dict get $FracturesDict $NewFractureId TipPoint Id]
            dict set FracturesDict $NewFractureId InterfaceSurface Id [expr { [GiD_Info Geometry MaxNumSurfaces]-1 }]
            dict set FracturesDict $NewFractureId InterfaceSurface Layer [dict get $FracturesDict $MotherFractureId InterfaceSurface Layer]
            dict set FracturesDict $NewFractureId InterfaceSurface Groups [dict get $FracturesDict $MotherFractureId InterfaceSurface Groups]
            dict set FracturesDict $NewFractureId InterfaceSurface TopLine $BotTopLine
            dict set FracturesDict $NewFractureId InterfaceSurface BotLine $BotBotLine
            
            dict set BodySurfacesDict $TopBodySurfaceId Lines $TopBodySurfaceLines
            dict set BodySurfacesDict $BotBodySurfaceId Lines $BotBodySurfaceLines
        }
    }
    
    # Create BodySurfaces again assign conditions and replace BodySurfacesDict and FracturesDict with the new values
    set BodySurfacesAuxDict $BodySurfacesDict
    unset BodySurfacesDict
    dict for {BodyId BodySurface} $BodySurfacesAuxDict {
        set BodySurfaceLines [dict get $BodySurface Lines]
        GiD_Process Mescape Geometry Create NurbsSurface
        for {set i 0} {$i < [llength $BodySurfaceLines]} {incr i} {
            GiD_Process [lindex $BodySurfaceLines $i]
        }
        GiD_Process escape escape
        
        set NewBodySurfaceId [GiD_Info Geometry MaxNumSurfaces]
        
        set BodySurfaceGroups [dict get $BodySurface Groups]
        for {set i 0} {$i < [llength $BodySurfaceGroups]} {incr i} {
            GiD_EntitiesGroups assign [lindex $BodySurfaceGroups $i] surfaces $NewBodySurfaceId
        }
                        
        dict set BodySurfacesDict $NewBodySurfaceId Layer [dict get $BodySurface Layer]
        dict set BodySurfacesDict $NewBodySurfaceId Groups [dict get $BodySurface Groups]
        dict set BodySurfacesDict $NewBodySurfaceId Lines [dict get $BodySurface Lines]
    }
        
    dict for {Id Fracture} $FracturesDict {
        set BodySurfaces [list]
        dict for {BodyId BodySurface} $BodySurfacesDict {
            if {[GiD_Info IsPointInside Surface $BodyId [dict get $Fracture TipPoint Coordinates]] == 1} {
                lappend BodySurfaces $BodyId
            }
        }
        dict set FracturesDict $Id BodySurfaces $BodySurfaces
    }
    
    # Generate New Mesh
    GiD_Process Mescape Meshing Generate Yes [GiD_Info Project LastElementSize] MeshingParametersFrom=Preferences
    GiD_Process Mescape Meshing MeshView
    
    # Update FracturesData.json file
    set filename [file join $dir FracturesData.json]
    set varfile [open $filename w]
        
    puts $varfile "\{"

    ## fracture_data
    puts $varfile "    \"fracture_data\": \{"
    puts $varfile "        \"gid_path\":                                 \"[lindex $PropagationData 0]\","
    puts $varfile "        \"propagation_length\":                       [GiD_AccessValue get gendata Propagation_Length],"
    puts $varfile "        \"propagation_damage\":                       [GiD_AccessValue get gendata Propagation_Damage],"
    puts $varfile "        \"propagation_width\":                        [GiD_AccessValue get gendata Propagation_Width],"
    puts $varfile "        \"propagation_cosangle\":                     [GiD_AccessValue get gendata Propagation_CosAngle],"
    puts $varfile "        \"propagation_frequency\":                    [GiD_AccessValue get gendata Propagation_Frequency],"
    # body_domain_sub_model_part_list
    set PutStrings \[
    set Groups [GiD_Info conditions Body_Part groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
        append PutStrings \" [lindex [lindex $Groups $i] 1] \" ,
    }
    set PutStrings [string trimright $PutStrings ,]
    append PutStrings \]
    puts $varfile "        \"body_domain_sub_model_part_list\":          $PutStrings,"
    # interface_domain_sub_model_part_old_list
    puts $varfile "        \"interface_domain_sub_model_part_old_list\": $InterfaceGroupsOld,"
    # interface_domain_sub_model_part_list
    set PutStrings \[
    set Groups [GiD_Info conditions Interface_Part groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
        append PutStrings \" [lindex [lindex $Groups $i] 1] \" ,
    }
    set PutStrings [string trimright $PutStrings ,]
    append PutStrings \]
    puts $varfile "        \"interface_domain_sub_model_part_list\":     $PutStrings"
    puts $varfile "    \},"
    
    ## fractures_list
    set iter 0
    puts $varfile "    \"fractures_list\": \[\{"
    dict for {Id Fracture} $FracturesDict {
        incr iter
        puts $varfile "        \"id\":                $Id,"
        # TipPoint
        puts $varfile "        \"tip_point\":         \{"
        puts $varfile "            \"id\":          [dict get $Fracture TipPoint Id],"
        puts $varfile "            \"layer\":       \"[dict get $Fracture TipPoint Layer]\","
        if {[llength [dict get $Fracture TipPoint Groups]]==0} {
            puts $varfile "            \"groups\":      \[\],"
        } else {
            set PutStrings \[
            for {set i 0} {$i < [llength [dict get $Fracture TipPoint Groups]]} {incr i} {
                append PutStrings \" [lindex [dict get $Fracture TipPoint Groups] $i] \" ,
            }
            set PutStrings [string trimright $PutStrings ,]
            append PutStrings \]
            puts $varfile "            \"groups\":      $PutStrings,"
        }
        puts $varfile "            \"coordinates\": \[[lindex [dict get $Fracture TipPoint Coordinates] 0],[lindex [dict get $Fracture TipPoint Coordinates] 1],[lindex [dict get $Fracture TipPoint Coordinates] 2]\]"
        puts $varfile "        \},"
        # TopPoint
        puts $varfile "        \"top_point\":         \{"
        puts $varfile "            \"id\":          [dict get $Fracture TopPoint Id],"
        puts $varfile "            \"layer\":       \"[dict get $Fracture TopPoint Layer]\","
        if {[llength [dict get $Fracture TopPoint Groups]]==0} {
            puts $varfile "            \"groups\":      \[\],"
        } else {
            set PutStrings \[
            for {set i 0} {$i < [llength [dict get $Fracture TopPoint Groups]]} {incr i} {
                append PutStrings \" [lindex [dict get $Fracture TopPoint Groups] $i] \" ,
            }
            set PutStrings [string trimright $PutStrings ,]
            append PutStrings \]
            puts $varfile "            \"groups\":      $PutStrings,"
        }
        puts $varfile "            \"coordinates\": \[[lindex [dict get $Fracture TopPoint Coordinates] 0],[lindex [dict get $Fracture TopPoint Coordinates] 1],[lindex [dict get $Fracture TopPoint Coordinates] 2]\]"
        puts $varfile "        \},"
        # BotPoint
        puts $varfile "        \"bot_point\":         \{"
        puts $varfile "            \"id\":          [dict get $Fracture BotPoint Id],"
        puts $varfile "            \"layer\":       \"[dict get $Fracture BotPoint Layer]\","
        if {[llength [dict get $Fracture BotPoint Groups]]==0} {
            puts $varfile "            \"groups\":      \[\],"
        } else {
            set PutStrings \[
            for {set i 0} {$i < [llength [dict get $Fracture BotPoint Groups]]} {incr i} {
                append PutStrings \" [lindex [dict get $Fracture BotPoint Groups] $i] \" ,
            }
            set PutStrings [string trimright $PutStrings ,]
            append PutStrings \]
            puts $varfile "            \"groups\":      $PutStrings,"
        }
        puts $varfile "            \"coordinates\": \[[lindex [dict get $Fracture BotPoint Coordinates] 0],[lindex [dict get $Fracture BotPoint Coordinates] 1],[lindex [dict get $Fracture BotPoint Coordinates] 2]\]"
        puts $varfile "        \},"
        # TopLine
        puts $varfile "        \"top_line\":          \{"
        puts $varfile "            \"id\":         [dict get $Fracture TopLine Id],"
        puts $varfile "            \"layer\":      \"[dict get $Fracture TopLine Layer]\","
        if {[llength [dict get $Fracture TopLine Groups]]==0} {
            puts $varfile "            \"groups\":     \[\],"
        } else {
            set PutStrings \[
            for {set i 0} {$i < [llength [dict get $Fracture TopLine Groups]]} {incr i} {
                append PutStrings \" [lindex [dict get $Fracture TopLine Groups] $i] \" ,
            }
            set PutStrings [string trimright $PutStrings ,]
            append PutStrings \]
            puts $varfile "            \"groups\":     $PutStrings,"
        }
        puts $varfile "            \"init_point\": [dict get $Fracture TopLine InitPoint],"
        puts $varfile "            \"end_point\":  [dict get $Fracture TopLine EndPoint]"
        puts $varfile "        \},"
        # BotLine
        puts $varfile "        \"bot_line\":          \{"
        puts $varfile "            \"id\":         [dict get $Fracture BotLine Id],"
        puts $varfile "            \"layer\":      \"[dict get $Fracture BotLine Layer]\","
        if {[llength [dict get $Fracture BotLine Groups]]==0} {
            puts $varfile "            \"groups\":     \[\],"
        } else {
            set PutStrings \[
            for {set i 0} {$i < [llength [dict get $Fracture BotLine Groups]]} {incr i} {
                append PutStrings \" [lindex [dict get $Fracture BotLine Groups] $i] \" ,
            }
            set PutStrings [string trimright $PutStrings ,]
            append PutStrings \]
            puts $varfile "            \"groups\":     $PutStrings,"
        }
        puts $varfile "            \"init_point\": [dict get $Fracture BotLine InitPoint],"
        puts $varfile "            \"end_point\":  [dict get $Fracture BotLine EndPoint]"
        puts $varfile "        \},"
        # InterfaceSurface
        puts $varfile "        \"interface_surface\": \{"
        puts $varfile "            \"id\":       [dict get $Fracture InterfaceSurface Id],"
        puts $varfile "            \"layer\":    \"[dict get $Fracture InterfaceSurface Layer]\","
        if {[llength [dict get $Fracture InterfaceSurface Groups]]==0} {
            puts $varfile "            \"groups\":   \[\],"
        } else {
            set PutStrings \[
            for {set i 0} {$i < [llength [dict get $Fracture InterfaceSurface Groups]]} {incr i} {
                append PutStrings \" [lindex [dict get $Fracture InterfaceSurface Groups] $i] \" ,
            }
            set PutStrings [string trimright $PutStrings ,]
            append PutStrings \]
            puts $varfile "            \"groups\":   $PutStrings,"
        }
        puts $varfile "            \"top_line\": \[[lindex [dict get $Fracture InterfaceSurface TopLine] 0],[lindex [dict get $Fracture InterfaceSurface TopLine] 1]\],"
        puts $varfile "            \"bot_line\": \[[lindex [dict get $Fracture InterfaceSurface BotLine] 0],[lindex [dict get $Fracture InterfaceSurface BotLine] 1]\]"
        puts $varfile "        \},"
        # BodySurfaces
        set PutStrings \[
        for {set i 0} {$i < [llength [dict get $Fracture BodySurfaces]]} {incr i} {
            append PutStrings [lindex [dict get $Fracture BodySurfaces] $i] ,
        }
        set PutStrings [string trimright $PutStrings ,]
        append PutStrings \]
        puts $varfile "        \"body_surfaces\":      $PutStrings"
        if {$iter < [dict size $FracturesDict]} {
            puts $varfile "    \},\{"
        } else {
            puts $varfile "    \}\],"
        }
    }
    
    ## body_surfaces_list
    set iter 0
    puts $varfile "    \"body_surfaces_list\": \[\{"
    dict for {Id BodySurface} $BodySurfacesDict {
        incr iter
        puts $varfile "        \"id\":     $Id,"
        puts $varfile "        \"layer\":  \"[dict get $BodySurface Layer]\","
        if {[llength [dict get $BodySurface Groups]]==0} {
            puts $varfile "        \"groups\": \[\],"
        } else {
            set PutStrings \[
            for {set i 0} {$i < [llength [dict get $BodySurface Groups]]} {incr i} {
                append PutStrings \" [lindex [dict get $BodySurface Groups] $i] \" ,
            }
            set PutStrings [string trimright $PutStrings ,]
            append PutStrings \]
            puts $varfile "        \"groups\": $PutStrings,"
        }
        set PutStrings \[
        for {set i 0} {$i < [llength [dict get $BodySurface Lines]]} {incr i} {
            append PutStrings [lindex [dict get $BodySurface Lines] $i] ,
        }
        set PutStrings [string trimright $PutStrings ,]
        append PutStrings \]
        puts $varfile "        \"lines\":  $PutStrings"
        if {$iter < [dict size $BodySurfacesDict]} {
            puts $varfile "    \},\{"
        } else {
            puts $varfile "    \}\]"
        }
    }
        
    puts $varfile "\}"
    
    close $varfile
}
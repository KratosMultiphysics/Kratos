proc WriteInitialFracturesData { dir problemtypedir gidpath } {
    
    ## Source auxiliar procedures
    source [file join $problemtypedir FracturePropagation2DAuxProcs.tcl]

    ## Set BodySurfacesDict
    set BodySurfacesDict [dict create]
    set BodyGroups [GiD_Info conditions Body_Part groups]
    for {set k 0} {$k < [llength $BodyGroups]} {incr k} {
        set BodyEntities [GiD_EntitiesGroups get [lindex [lindex $BodyGroups $k] 1] surfaces]
        for {set l 0} {$l < [llength $BodyEntities]} {incr l} {
            set BodySurface [GiD_Geometry get surface [lindex $BodyEntities $l]]
            set Groups [GiD_EntitiesGroups entity_groups surfaces [lindex $BodyEntities $l]]
            dict set BodySurfacesDict [lindex $BodyEntities $l] Groups $Groups
            set Lines [list]
            for {set m 0} {$m < [lindex $BodySurface 2]} {incr m} {
                lappend Lines [lindex [lindex $BodySurface [expr { 9+$m }]] 0]
            }
            dict set BodySurfacesDict [lindex $BodyEntities $l] Lines $Lines
            if {[lindex [GiD_Info list_entities Surfaces [lindex $BodyEntities $l]] 11] eq "Meshing"} {
                set ElemType [lindex [GiD_Info list_entities Surfaces [lindex $BodyEntities $l]] 14]
                if {$ElemType eq "Elemtype=0" || $ElemType eq "Elemtype=2"} {
                    dict set BodySurfacesDict [lindex $BodyEntities $l] ElemType Triangle
                } elseif {$ElemType eq "Elemtype=3"} {
                    dict set BodySurfacesDict [lindex $BodyEntities $l] ElemType Quadrilateral
                } else {
                    dict set BodySurfacesDict [lindex $BodyEntities $l] ElemType Triangle
                }
                set MeshSize [lindex [GiD_Info list_entities Surfaces [lindex $BodyEntities $l]] 17]
                set MeshSize [string trimleft $MeshSize "size="]
                dict set BodySurfacesDict [lindex $BodyEntities $l] MeshSize $MeshSize
            } else {
                dict set BodySurfacesDict [lindex $BodyEntities $l] ElemType Triangle
                dict set BodySurfacesDict [lindex $BodyEntities $l] MeshSize 0
            }
        }
    }

    ## Set FracturesDict
    set FracturesDict [dict create]
    set FractureId -1
    set InterfaceGroups [GiD_Info conditions Interface_Part groups]
    if {[llength $InterfaceGroups] < 1} {
        WarnWin "ERROR: Fracture Propagation needs at least 1 pre-defined fracture tip"
    }
    for {set i 0} {$i < [llength $InterfaceGroups]} {incr i} {
        if {[lindex [lindex $InterfaceGroups $i] 3] eq false} {
            set InterfaceEntities [GiD_EntitiesGroups get [lindex [lindex $InterfaceGroups $i] 1] surfaces]
            for {set j 0} {$j < [llength $InterfaceEntities]} {incr j} {
                set InterfaceSurface [GiD_Geometry get surface [lindex $InterfaceEntities $j]]
                set BotLine [GiD_Geometry get line [lindex [lindex $InterfaceSurface 3] 0]]
                set TopLine [GiD_Geometry get line [lindex [lindex $InterfaceSurface 4] 0]]
                if {([lindex $BotLine 2] eq [lindex $TopLine 2]) || ([lindex $BotLine 3] eq [lindex $TopLine 2])} {
                    # Define new fracture
                    incr FractureId
                    # Set TipPoint
                    set PointId [lindex $TopLine 2]
                    AddPointToFracturesDict FracturesDict $FractureId $PointId [GiD_Geometry get point $PointId] TipPoint
                    # Set TopPoint
                    set PointId [lindex $TopLine 3]
                    AddPointToFracturesDict FracturesDict $FractureId $PointId [GiD_Geometry get point $PointId] TopPoint
                    # Set BotPoint
                    if {[lindex $BotLine 2] != [lindex $TopLine 2]} {
                        set PointId [lindex $BotLine 2]
                        AddPointToFracturesDict FracturesDict $FractureId $PointId [GiD_Geometry get point $PointId] BotPoint
                    } else {
                        set PointId [lindex $BotLine 3]
                        AddPointToFracturesDict FracturesDict $FractureId $PointId [GiD_Geometry get point $PointId] BotPoint
                    }
                    # Set TopLine
                    dict set FracturesDict $FractureId TopLine Id [lindex [lindex $InterfaceSurface 4] 0]
                    # Set BotLine
                    dict set FracturesDict $FractureId BotLine Id [lindex [lindex $InterfaceSurface 3] 0]
                    # Set InterfaceSurface
                    AddInterfaceSurfaceToFracturesDict FracturesDict $FractureId [lindex $InterfaceEntities $j] $InterfaceSurface
                    # Set BodySurface
                    AddBodySurfaceToFracturesDict FracturesDict $FractureId $BodySurfacesDict

                    # We need to check that the cross product between the vector in local x and the vector in local y points towards (0,0,1)
                    CheckJonintOrientation FracturesDict $FractureId
                    
                } elseif {([lindex $BotLine 2] eq [lindex $TopLine 3]) || ([lindex $BotLine 3] eq [lindex $TopLine 3])} {
                    # Define new fracture
                    incr FractureId
                    # Set TipPoint
                    set PointId [lindex $TopLine 3]
                    AddPointToFracturesDict FracturesDict $FractureId $PointId [GiD_Geometry get point $PointId] TipPoint
                    # Set TopPoint
                    set PointId [lindex $TopLine 2]
                    AddPointToFracturesDict FracturesDict $FractureId $PointId [GiD_Geometry get point $PointId] TopPoint
                    # Set BotPoint
                    if {[lindex $BotLine 2] != [lindex $TopLine 3]} {
                        set PointId [lindex $BotLine 2]
                        AddPointToFracturesDict FracturesDict $FractureId $PointId [GiD_Geometry get point $PointId] BotPoint
                    } else {
                        set PointId [lindex $BotLine 3]
                        AddPointToFracturesDict FracturesDict $FractureId $PointId [GiD_Geometry get point $PointId] BotPoint
                    }
                    # Set TopLine
                    dict set FracturesDict $FractureId TopLine Id [lindex [lindex $InterfaceSurface 4] 0]
                    # Set BotLine
                    dict set FracturesDict $FractureId BotLine Id [lindex [lindex $InterfaceSurface 3] 0]
                    # Set InterfaceSurface
                    AddInterfaceSurfaceToFracturesDict FracturesDict $FractureId [lindex $InterfaceEntities $j] $InterfaceSurface
                    # Set BodySurface
                    AddBodySurfaceToFracturesDict FracturesDict $FractureId $BodySurfacesDict

                    # We need to check that the cross product between the vector in local x and the vector in local y points towards (0,0,1)
                    CheckJonintOrientation FracturesDict $FractureId
                }
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
    set PutStrings [string trimright $PutStrings ,]
    append PutStrings \]
    puts $FileVar1 "        \"interface_domain_sub_model_part_list\": $PutStrings"
    puts $FileVar1 "    \},"

    ## body_surfaces_list
    WriteBodySurfacesList FileVar1 $BodySurfacesDict

    ## fractures_list
    WriteFracturesList FileVar1 $FracturesDict
    
    puts $FileVar1 "\}"
    
    close $FileVar1
}

proc GenerateNewFractures { dir problemtypedir PropagationData } {
    
    ## Source auxiliar procedures
    source [file join $problemtypedir FracturePropagation2DAuxProcs.tcl]
    
    # Previous Fractures dictionaries
    set BodySurfacesDict [lindex $PropagationData 1]
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
    set InterfaceGroupsOld [string trimright $InterfaceGroupsOld ,]
    append InterfaceGroupsOld \]
    
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
        # Create new poin in TopInitCoordinates location
        GiD_Process Mescape Geometry Create Point \
            [lindex [dict get $Propagation TopInitCoordinates] 0] [lindex [dict get $Propagation TopInitCoordinates] 1] [lindex [dict get $Propagation TopInitCoordinates] 2] escape
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
            [lindex [dict get $Propagation BotInitCoordinates] 0] [lindex [dict get $Propagation BotInitCoordinates] 1] [lindex [dict get $Propagation BotInitCoordinates] 2] escape
        # Delete old BotLine
        GiD_Process Mescape Geometry Delete Lines [dict get $FracturesDict $MotherFractureId BotLine Id] escape
        set Index [lsearch $BodySurfaceLines [dict get $FracturesDict $MotherFractureId BotLine Id]]
        set BodySurfaceLines [lreplace $BodySurfaceLines $Index $Index]
        # Create new line replacing old BotLine
        GiD_Process Mescape Geometry Create Line Join \
            [dict get $FracturesDict $MotherFractureId BotPoint Id] [GiD_Info Geometry MaxNumPoints] escape escape
        lappend BodySurfaceLines [GiD_Info Geometry MaxNumLines]
        # Delete old TipPoint
        GiD_Process Mescape Geometry Delete Points [dict get $FracturesDict $MotherFractureId TipPoint Id] escape
        # Create ContactSurface for the old crack
        set BotLine [list [GiD_Info Geometry MaxNumLines] 1]
        set TopLine [list [expr { [GiD_Info Geometry MaxNumLines]-1 }] 0]
        # Note: the orientation 1 for the first line means that the surface is at the left side of the line in the advancing sense
        # Note: the orientation 0 for the second line means that the surface is at the left side of the line in the advancing sense
        GiD_Geometry create surface [dict get $FracturesDict $MotherFractureId InterfaceSurface Id] contactsurface \
            [dict get $FracturesDict $MotherFractureId InterfaceSurface Layer] 2 $BotLine $TopLine
        
        # Create new TipPoint
        GiD_Process Mescape Geometry Create Point \
            [lindex [dict get $Propagation TipCoordinates] 0] [lindex [dict get $Propagation TipCoordinates] 1] [lindex [dict get $Propagation TipCoordinates] 2] escape
        # Create new TopLine
        GiD_Process Mescape Geometry Create Line Join \
            [GiD_Info Geometry MaxNumPoints] [expr { [GiD_Info Geometry MaxNumPoints]-2 }] escape escape
        lappend BodySurfaceLines [GiD_Info Geometry MaxNumLines]
        # Create new BotLine
        GiD_Process Mescape Geometry Create Line Join \
            [expr { [GiD_Info Geometry MaxNumPoints]-1 }] [GiD_Info Geometry MaxNumPoints] escape escape
        lappend BodySurfaceLines [GiD_Info Geometry MaxNumLines]
        # Create new InterfaceContactSurface
        set BotLine [list [GiD_Info Geometry MaxNumLines] 1]
        set TopLine [list [expr { [GiD_Info Geometry MaxNumLines]-1 }] 0]
        GiD_Geometry create surface append contactsurface \
            [dict get $FracturesDict $MotherFractureId InterfaceSurface Layer] 2 $BotLine $TopLine

        ## Set Conditions
        for {set i 0} {$i < [llength [dict get $FracturesDict $MotherFractureId InterfaceSurface Groups]]} {incr i} {
            GiD_EntitiesGroups assign [lindex [dict get $FracturesDict $MotherFractureId InterfaceSurface Groups] $i] surfaces [dict get $FracturesDict $MotherFractureId InterfaceSurface Id]
            GiD_EntitiesGroups assign [lindex [dict get $FracturesDict $MotherFractureId InterfaceSurface Groups] $i] surfaces [GiD_Info Geometry MaxNumSurfaces]
        }
        
        ## Set Mesh options
        # GiD_Process Mescape Meshing Structured Lines 2 \
        #     [expr { [GiD_Info Geometry MaxNumLines]-3 }] [expr { [GiD_Info Geometry MaxNumLines]-2 }] \
        #     [expr { [GiD_Info Geometry MaxNumLines]-1 }] [GiD_Info Geometry MaxNumLines] escape escape

        ## Replace $FracturesDict $MotherFractureId by the new fracture
        dict set FracturesDict $MotherFractureId TopPoint Id [expr { [GiD_Info Geometry MaxNumPoints]-2 }]
        dict set FracturesDict $MotherFractureId TopPoint Coordinates [dict get $Propagation TopInitCoordinates]
        dict set FracturesDict $MotherFractureId BotPoint Id [expr { [GiD_Info Geometry MaxNumPoints]-1 }]
        dict set FracturesDict $MotherFractureId BotPoint Coordinates [dict get $Propagation BotInitCoordinates]
        dict set FracturesDict $MotherFractureId TipPoint Id [GiD_Info Geometry MaxNumPoints]
        dict set FracturesDict $MotherFractureId TipPoint Coordinates [dict get $Propagation TipCoordinates]
        dict set FracturesDict $MotherFractureId TopLine Id [expr { [GiD_Info Geometry MaxNumLines]-1 }]
        dict set FracturesDict $MotherFractureId BotLine Id [GiD_Info Geometry MaxNumLines]
        dict set FracturesDict $MotherFractureId InterfaceSurface Id [GiD_Info Geometry MaxNumSurfaces]

        dict set BodySurfacesDict $BodySurfaceId Lines $BodySurfaceLines
    }

    # Loop through all Bifurcation Fractures
    dict for {BifId Bifurcation} $BifurcationDict {
        set MotherFractureId [dict get $Bifurcation MotherFractureId]

        set BodySurfaceId [lindex [dict get $FracturesDict $MotherFractureId BodySurfaces] 0]
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
        set BotLine [list [GiD_Info Geometry MaxNumLines] 1]
        set TopLine [list [expr { [GiD_Info Geometry MaxNumLines]-1 }] 0]
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
        set TopBotLine [list [GiD_Info Geometry MaxNumLines] 1]
        set TopTopLine [list [expr { [GiD_Info Geometry MaxNumLines]-1 }] 0]
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
        set BotBotLine [list [GiD_Info Geometry MaxNumLines] 1]
        set BotTopLine [list [expr { [GiD_Info Geometry MaxNumLines]-1 }] 0]
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
        # GiD_Process Mescape Meshing Structured Lines 2 \
        #     [expr { [GiD_Info Geometry MaxNumLines]-7 }] [expr { [GiD_Info Geometry MaxNumLines]-8 }] \
        #     [expr { [GiD_Info Geometry MaxNumLines]-3 }] [expr { [GiD_Info Geometry MaxNumLines]-4 }] \
        #     [expr { [GiD_Info Geometry MaxNumLines]-5 }] [expr { [GiD_Info Geometry MaxNumLines]-6 }] escape escape
        GiD_Process Mescape Meshing AssignSizes Points 0.0 [dict get $FracturesDict $MotherFractureId TipPoint Id] escape escape
        GiD_Process Mescape Meshing Structured Surfaces [GiD_Info Geometry MaxNumSurfaces] escape 1 [GiD_Info Geometry MaxNumLines] escape escape
        
        ## Replace $FracturesDict $MotherFractureId by the new TopTip fracture
        dict set FracturesDict $MotherFractureId TopPoint Id [expr { [GiD_Info Geometry MaxNumPoints]-3 }]
        dict set FracturesDict $MotherFractureId TopPoint Coordinates [dict get $Bifurcation TopInitCoordinates]
        dict set FracturesDict $MotherFractureId BotPoint Id [dict get $FracturesDict $MotherFractureId TipPoint Id]
        dict set FracturesDict $MotherFractureId BotPoint Coordinates [dict get $FracturesDict $MotherFractureId TipPoint Coordinates]
        dict set FracturesDict $MotherFractureId TipPoint Id [expr { [GiD_Info Geometry MaxNumPoints]-1 }]
        dict set FracturesDict $MotherFractureId TipPoint Coordinates [dict get $Bifurcation TopTipCoordinates]
        dict set FracturesDict $MotherFractureId TopLine Id [expr { [GiD_Info Geometry MaxNumLines]-6 }]
        dict set FracturesDict $MotherFractureId BotLine Id [expr { [GiD_Info Geometry MaxNumLines]-5 }]
        dict set FracturesDict $MotherFractureId InterfaceSurface Id [expr { [GiD_Info Geometry MaxNumSurfaces]-2 }]
        
        ## Add an additional value in $FracturesDict for the new BotTip fracture
        set NewFractureId [dict size $FracturesDict]
        dict set FracturesDict $NewFractureId TipPoint Id [GiD_Info Geometry MaxNumPoints]
        dict set FracturesDict $NewFractureId TipPoint Coordinates [dict get $Bifurcation BotTipCoordinates]
        dict set FracturesDict $NewFractureId TopPoint Id [dict get $FracturesDict $MotherFractureId BotPoint Id]
        dict set FracturesDict $NewFractureId TopPoint Coordinates [dict get $FracturesDict $MotherFractureId BotPoint Coordinates]
        dict set FracturesDict $NewFractureId BotPoint Id [expr { [GiD_Info Geometry MaxNumPoints]-2 }]
        dict set FracturesDict $NewFractureId BotPoint Coordinates [dict get $Bifurcation BotInitCoordinates]
        dict set FracturesDict $NewFractureId TopLine Id [expr { [GiD_Info Geometry MaxNumLines]-4 }]
        dict set FracturesDict $NewFractureId BotLine Id [expr { [GiD_Info Geometry MaxNumLines]-3 }]
        dict set FracturesDict $NewFractureId InterfaceSurface Id [expr { [GiD_Info Geometry MaxNumSurfaces]-1 }]
        dict set FracturesDict $NewFractureId InterfaceSurface Layer [dict get $FracturesDict $MotherFractureId InterfaceSurface Layer]
        dict set FracturesDict $NewFractureId InterfaceSurface Groups [dict get $FracturesDict $MotherFractureId InterfaceSurface Groups]
        
        dict set BodySurfacesDict $BodySurfaceId Lines $BodySurfaceLines
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
        
        GiD_Process Mescape Meshing ElemType [dict get $BodySurface ElemType] $NewBodySurfaceId escape escape
        if {[dict get $BodySurface MeshSize] > 0.0} {
            GiD_Process Mescape Meshing AssignSizes Surfaces [dict get $BodySurface MeshSize] $NewBodySurfaceId escape escape
        }

        dict set BodySurfacesDict $NewBodySurfaceId Groups [dict get $BodySurface Groups]
        dict set BodySurfacesDict $NewBodySurfaceId Lines [dict get $BodySurface Lines]
        dict set BodySurfacesDict $NewBodySurfaceId ElemType [dict get $BodySurface ElemType]
        dict set BodySurfacesDict $NewBodySurfaceId MeshSize [dict get $BodySurface MeshSize]
    }

    dict for {Id Fracture} $FracturesDict {
        set BodySurfaces [list]
        dict for {BodyId BodySurface} $BodySurfacesDict {
            if {[GiD_Info IsPointInside Surface $BodyId [dict get $Fracture TipPoint Coordinates]] eq 1} {
                lappend BodySurfaces $BodyId
            }
        }
        dict set FracturesDict $Id BodySurfaces $BodySurfaces
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
    set PutStrings [string trimright $PutStrings ,]
    append PutStrings \]
    puts $FileVar "        \"interface_domain_sub_model_part_list\": $PutStrings,"
    # interface_domain_sub_model_part_old_list
    puts $FileVar "        \"interface_domain_sub_model_part_old_list\": $InterfaceGroupsOld"
    puts $FileVar "    \},"

    ## body_surfaces_list
    WriteBodySurfacesList FileVar $BodySurfacesDict

    ## fractures_list
    WriteFracturesList FileVar $FracturesDict

    puts $FileVar "\}"

    close $FileVar
}

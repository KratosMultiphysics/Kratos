proc WriteMdpa { basename dir problemtypedir } {

    ## Source auxiliar procedures
    source [file join $problemtypedir MdpaAuxProcs.tcl]

    ## Start MDPA file
    set filename [file join $dir ${basename}.mdpa]
    set FileVar [open $filename w]

    ## Properties
    set PropertyId 0
    set PropertyDict [dict create]
    # Source_terms
    set Groups [GiD_Info conditions Computing_domain groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
        dict set PropertyDict [lindex [lindex $Groups $i] 1] $PropertyId
        puts $FileVar "Begin Properties $PropertyId"
        puts $FileVar "End Properties"
        puts $FileVar ""
        incr PropertyId
    }

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

    ## Elements
    # Computing_domain
    set Groups [GiD_Info conditions Computing_domain groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
        # Elements Property
        set BodyElemsProp [dict get $PropertyDict [lindex [lindex $Groups $i] 1]]
        # Trinagles: Element2D3N
        WriteElements FileVar [lindex $Groups $i] triangle Element2D3N $BodyElemsProp Triangle2D3Connectivities
        # Quadrilaterals: Element2D4N
        WriteElements FileVar [lindex $Groups $i] quadrilateral Element2D4N $BodyElemsProp Quadrilateral2D4Connectivities
    }


    ## Conditions
    set ConditionId 0
    set ConditionDict [dict create]
    set Dim 2

    # Slip_condition
    set Groups [GiD_Info conditions Slip_condition groups]
    if {$Dim eq 2} {
        WriteFaceConditions FileVar ConditionId ConditionDict $Groups LineCondition2D2N $PropertyDict
    }

    # Water_height
    set Groups [GiD_Info conditions Water_height groups]
    if {$Dim eq 2} {
        WriteFaceConditions FileVar ConditionId ConditionDict $Groups LineCondition2D2N $PropertyDict
    }

    # Imposed_flow_rate
    set Groups [GiD_Info conditions Imposed_flow_rate groups]
    if {$Dim eq 2} {
        WriteFaceConditions FileVar ConditionId ConditionDict $Groups LineCondition2D2N $PropertyDict
    }


    puts $FileVar ""

    ## SubModelParts
    # Topographic data
    WriteElementSubmodelPart FileVar Topography
    WriteElementSubmodelPart FileVar Bottom_friction
    # Initial_water_level
    WriteElementSubmodelPart FileVar Initial_water_level
    # Conditions
    WriteLoadSubmodelPart FileVar Slip_condition $ConditionDict
    WriteLoadSubmodelPart FileVar Water_height $ConditionDict
    WriteLoadSubmodelPart FileVar Imposed_flow_rate $ConditionDict


    close $FileVar
}

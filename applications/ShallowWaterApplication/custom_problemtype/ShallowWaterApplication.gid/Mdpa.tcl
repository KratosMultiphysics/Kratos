proc WriteMdpa { basename dir problemtypedir } {

    ## Source auxiliar procedures
    source [file join $problemtypedir MdpaAuxProcs.tcl]

    ## Start MDPA file
    set filename [file join $dir ${basename}.mdpa]
    set FileVar [open $filename w]

    ## Properties
    set PropertyId 0
    set PropertyDict [dict create]
    puts $FileVar "Begin Properties 0"
    puts $FileVar "End Properties"
    # Source_terms
    set Groups [GiD_Info conditions Bottom_friction groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
        incr PropertyId
        dict set PropertyDict [lindex [lindex $Groups $i] 1] $PropertyId
        puts $FileVar "Begin Properties $PropertyId"
        puts $FileVar "  MANNING [lindex [lindex $Groups $i] 3]"
        puts $FileVar "End Properties"
        puts $FileVar ""
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
    set VariablesType [GiD_AccessValue get gendata Variables]
    set FrameworkType [GiD_AccessValue get gendata Framework]
    # Body_Part
    set Groups [GiD_Info conditions Bottom_friction groups]
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
        WriteFaceConditions FileVar ConditionId ConditionDict $Groups Condition2D2N $PropertyDict
    }

    # Water_height
    set Groups [GiD_Info conditions Water_height groups]
    if {$Dim eq 2} {
        WriteFaceConditions FileVar ConditionId ConditionDict $Groups Condition2D2N $PropertyDict
    }

    # Imposed_flux
    set Groups [GiD_Info conditions Imposed_flux groups]
    if {$Dim eq 2} {
        WriteFaceConditions FileVar ConditionId ConditionDict $Groups Condition2D2N $PropertyDict
    }


    puts $FileVar ""

    ## SubModelParts
    # Body_Part
    WriteElementSubmodelPart FileVar Body_Part
    # Initial_water_level
    WriteConstraintSubmodelPart FileVar Initial_water_level
    # Slip_condition
    WriteLoadSubmodelPart FileVar Slip_condition $ConditionDict
    # Water_height
    WriteLoadSubmodelPart FileVar Water_height $ConditionDict
    # Imposed_flux
    WriteLoadSubmodelPart FileVar Imposed_flux $ConditionDict


    close $FileVar
}

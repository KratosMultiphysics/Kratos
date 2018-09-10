
# TODO: it may be dangerous to write Tables without format (puts -nonewline $FileVar [format  "%.10f" [lindex $Table $j]])

#-------------------------------------------------------------------------------

proc PressureTable {FileVar TableId TableDict CondName VarName} {
    set Groups [GiD_Info conditions $CondName groups]
    if {[llength $Groups]>0} {
        upvar $FileVar MyFileVar
        upvar $TableId MyTableId
        upvar $TableDict MyTableDict
        
        for {set i 0} {$i < [llength $Groups]} {incr i} {
            set AuxList [list]
            if {[lindex [lindex $Groups $i] 9] eq "Table_Interpolation"} {
                incr MyTableId
                dict set MyTableDict [lindex [lindex $Groups $i] 1] Table0 $MyTableId
                lappend AuxList $MyTableId
                puts $MyFileVar "Begin Table $MyTableId TIME $VarName"
                set Table [lindex [lindex $Groups $i] 10]
                for {set j 2} {$j <= [lindex $Table 1]} {incr j 2} {
                    puts $MyFileVar "  [lindex $Table $j] [lindex $Table [expr { $j+1 }]]"
                }
                puts $MyFileVar "End Table"
                puts $MyFileVar ""
            } else {
                dict set MyTableDict [lindex [lindex $Groups $i] 1] Table0 0
            }
            dict set MyTableDict [lindex [lindex $Groups $i] 1] TableList $AuxList
        }
    }
}

#-------------------------------------------------------------------------------

proc ScalarTable {FileVar TableId TableDict CondName VarName} {
    set Groups [GiD_Info conditions $CondName groups]
    if {[llength $Groups]>0} {
        upvar $FileVar MyFileVar
        upvar $TableId MyTableId
        upvar $TableDict MyTableDict
        
        for {set i 0} {$i < [llength $Groups]} {incr i} {
            set AuxList [list]
            if {[lindex [lindex $Groups $i] 4] eq "Table_Interpolation"} {
                incr MyTableId
                dict set MyTableDict [lindex [lindex $Groups $i] 1] Table0 $MyTableId
                lappend AuxList $MyTableId
                puts $MyFileVar "Begin Table $MyTableId TIME $VarName"
                set Table [lindex [lindex $Groups $i] 5]
                for {set j 2} {$j <= [lindex $Table 1]} {incr j 2} {
                    puts $MyFileVar "  [lindex $Table $j] [lindex $Table [expr { $j+1 }]]"
                }
                puts $MyFileVar "End Table"
                puts $MyFileVar ""
            } else {
                dict set MyTableDict [lindex [lindex $Groups $i] 1] Table0 0
            }
            dict set MyTableDict [lindex [lindex $Groups $i] 1] TableList $AuxList
        }
    }
}


#--------------------------------------------------------------------------------------------------------------------------------------------------------------


proc WriteElements {FileVar Group ElemType ElemName PropertyId ConnectivityType} {
    set Entities [GiD_EntitiesGroups get [lindex $Group 1] elements -element_type $ElemType]
    if {[llength $Entities] > 0} {
        upvar $FileVar MyFileVar
        
        puts $MyFileVar "Begin Elements $ElemName"
        for {set j 0} {$j < [llength $Entities]} {incr j} {
            puts $MyFileVar "  [lindex $Entities $j]  $PropertyId  [$ConnectivityType [lindex $Entities $j]]"
        }
        puts $MyFileVar "End Elements"
        puts $MyFileVar ""
    }
}

#-------------------------------------------------------------------------------

proc WriteNodalConditions {FileVar ConditionId ConditionDict Groups CondName PropertyId} {
    if {[llength $Groups]>0} {
        upvar $FileVar MyFileVar
        upvar $ConditionId MyConditionId
        upvar $ConditionDict MyConditionDict
                
        for {set i 0} {$i < [llength $Groups]} {incr i} {
            set MyConditionList [list]
            set Entities [GiD_EntitiesGroups get [lindex [lindex $Groups $i] 1] nodes]
            puts $MyFileVar "Begin Conditions $CondName"
            for {set j 0} {$j < [llength $Entities]} {incr j} {
                incr MyConditionId
                lappend MyConditionList $MyConditionId
                puts $MyFileVar "  $MyConditionId  $PropertyId  [lindex $Entities $j]"
            }
            puts $MyFileVar "End Conditions"
            puts $MyFileVar ""
            dict set MyConditionDict [lindex [lindex $Groups $i] 1] $MyConditionList
        }
    }
}

#-------------------------------------------------------------------------------

proc WriteFaceConditions {FileVar ConditionId ConditionDict Groups CondName PropertyDict} {
    if {[llength $Groups]>0} {
        upvar $FileVar MyFileVar
        upvar $ConditionId MyConditionId
        upvar $ConditionDict MyConditionDict
        
        for {set i 0} {$i < [llength $Groups]} {incr i} {
            set MyConditionList [list]
            set Entities [GiD_EntitiesGroups get [lindex [lindex $Groups $i] 1] faces]
            puts $MyFileVar "Begin Conditions $CondName"
            for {set j 0} {$j < [llength [lindex $Entities 1]]} {incr j} {
                incr MyConditionId
                lappend MyConditionList $MyConditionId
                #~ set ElementGroup [GiD_EntitiesGroups entity_groups element [lindex [lindex $Entities 0] $j]]
                #~ for {set k 0} {$k < [llength $ElementGroup]} {incr k} {
                    #~ if {[dict exists $PropertyDict [lindex $ElementGroup $k]] eq 1} {
                        #~ set PropertyId [dict get $PropertyDict [lindex $ElementGroup $k]]
                        #~ break
                    #~ }
                #~ }
                set Connectivities [GiD_Mesh get element [lindex [lindex $Entities 0] $j] face [lindex [lindex $Entities 1] $j]]
                puts $MyFileVar "  $MyConditionId  0  $Connectivities"
            }
            puts $MyFileVar "End Conditions"
            puts $MyFileVar ""
            dict set MyConditionDict [lindex [lindex $Groups $i] 1] $MyConditionList
        }
    }
}

#-------------------------------------------------------------------------------

proc WriteTypeFaceConditions {FileVar ConditionId ConditionList Group ElemType CondName PropertyDict} {
    set Entities [GiD_EntitiesGroups get [lindex $Group 1] faces -element_type $ElemType]
    if {[llength [lindex $Entities 1]] > 0} {
        upvar $FileVar MyFileVar
        upvar $ConditionId MyConditionId
        upvar $ConditionList MyConditionList

        puts $MyFileVar "Begin Conditions $CondName"
        for {set j 0} {$j < [llength [lindex $Entities 1]]} {incr j} {
            incr MyConditionId
            lappend MyConditionList $MyConditionId
            set ElementGroup [GiD_EntitiesGroups entity_groups element [lindex [lindex $Entities 0] $j]]
            for {set k 0} {$k < [llength ElementGroup]} {incr k} {
                if {[dict exists $PropertyDict [lindex $ElementGroup $k]] eq 1} {
                    set PropertyId [dict get $PropertyDict [lindex $ElementGroup $k]]
                    break
                }
            }
            set Connectivities [GiD_Mesh get element [lindex [lindex $Entities 0] $j] face [lindex [lindex $Entities 1] $j]]
            puts $MyFileVar "  $MyConditionId  $PropertyId  $Connectivities"
        }
        puts $MyFileVar "End Conditions"
        puts $MyFileVar ""
    }
}


#-------------------------------------------------------------------------------

proc WriteInterfaceConditions {FileVar ConditionId ConditionList Group ElemType CondName PropertyId ConnectivityType} {
    set Entities [GiD_EntitiesGroups get [lindex $Group 1] elements -element_type $ElemType]
    if {[llength $Entities] > 0} {
        upvar $FileVar MyFileVar
        upvar $ConditionId MyConditionId
        upvar $ConditionList MyConditionList

        puts $MyFileVar "Begin Conditions $CondName"
        for {set j 0} {$j < [llength $Entities]} {incr j} {
            incr MyConditionId
            lappend MyConditionList $MyConditionId
            puts $MyFileVar "  $MyConditionId  $PropertyId  [$ConnectivityType [lindex $Entities $j]]"
        }
        puts $MyFileVar "End Conditions"
        puts $MyFileVar ""
    }
}

#-------------------------------------------------------------------------------

proc Triangle2D3Connectivities { ElemId } {
    
    set ElementInfo [GiD_Mesh get element $ElemId]
    #ElementInfo: <layer> <elemtype> <NumNodes> <N1> <N2> ...
    return "[lindex $ElementInfo 3] [lindex $ElementInfo 4] [lindex $ElementInfo 5]"
}

#-------------------------------------------------------------------------------

proc Quadrilateral2D4Connectivities { ElemId } {
    
    #Note: It is the same for the Tethrahedron3D4
    
    set ElementInfo [GiD_Mesh get element $ElemId]
    #ElementInfo: <layer> <elemtype> <NumNodes> <N1> <N2> ...
    return "[lindex $ElementInfo 3] [lindex $ElementInfo 4] [lindex $ElementInfo 5]\
    [lindex $ElementInfo 6]"
}


#-------------------------------------------------------------------------------

proc Triangle2D6Connectivities { ElemId } {
    
    #It is the same for the Prism3D6
    
    set ElementInfo [GiD_Mesh get element $ElemId]
    #ElementInfo: <layer> <elemtype> <NumNodes> <N1> <N2> ...
    return "[lindex $ElementInfo 3] [lindex $ElementInfo 4] [lindex $ElementInfo 5]\
    [lindex $ElementInfo 6] [lindex $ElementInfo 7] [lindex $ElementInfo 8]"
}

#-------------------------------------------------------------------------------

proc Hexahedron3D8Connectivities { ElemId } {
    
    #It is the same for Quadrilateral2D8
    
    set ElementInfo [GiD_Mesh get element $ElemId]
    #ElementInfo: <layer> <elemtype> <NumNodes> <N1> <N2> ...
    return "[lindex $ElementInfo 3] [lindex $ElementInfo 4] [lindex $ElementInfo 5]\
    [lindex $ElementInfo 6] [lindex $ElementInfo 7] [lindex $ElementInfo 8]\
    [lindex $ElementInfo 9] [lindex $ElementInfo 10]"
}

#-------------------------------------------------------------------------------

proc Quadrilateral2D9Connectivities { ElemId } {
    
    set ElementInfo [GiD_Mesh get element $ElemId]
    #ElementInfo: <layer> <elemtype> <NumNodes> <N1> <N2> ...
    return "[lindex $ElementInfo 3] [lindex $ElementInfo 4] [lindex $ElementInfo 5]\
    [lindex $ElementInfo 6] [lindex $ElementInfo 7] [lindex $ElementInfo 8]\
    [lindex $ElementInfo 9] [lindex $ElementInfo 10] [lindex $ElementInfo 11]"
}

#-------------------------------------------------------------------------------

proc Tetrahedron3D10Connectivities { ElemId } {
    
    set ElementInfo [GiD_Mesh get element $ElemId]
    #ElementInfo: <layer> <elemtype> <NumNodes> <N1> <N2> ...
    return "[lindex $ElementInfo 3] [lindex $ElementInfo 4] [lindex $ElementInfo 5]\
    [lindex $ElementInfo 6] [lindex $ElementInfo 7] [lindex $ElementInfo 8]\
    [lindex $ElementInfo 9] [lindex $ElementInfo 10] [lindex $ElementInfo 11]\
    [lindex $ElementInfo 12]"
}

#-------------------------------------------------------------------------------

proc Hexahedron3D20Connectivities { ElemId } {
    
    set ElementInfo [GiD_Mesh get element $ElemId]
    #ElementInfo: <layer> <elemtype> <NumNodes> <N1> <N2> ...
    return "[lindex $ElementInfo 3] [lindex $ElementInfo 4] [lindex $ElementInfo 5]\
    [lindex $ElementInfo 6] [lindex $ElementInfo 7] [lindex $ElementInfo 8]\
    [lindex $ElementInfo 9] [lindex $ElementInfo 10] [lindex $ElementInfo 11]\
    [lindex $ElementInfo 12] [lindex $ElementInfo 13] [lindex $ElementInfo 14]\
    [lindex $ElementInfo 15] [lindex $ElementInfo 16] [lindex $ElementInfo 17]\
    [lindex $ElementInfo 18] [lindex $ElementInfo 19] [lindex $ElementInfo 20]\
    [lindex $ElementInfo 21] [lindex $ElementInfo 22]"
}

#-------------------------------------------------------------------------------

proc Hexahedron3D27Connectivities { ElemId } {
        
    set ElementInfo [GiD_Mesh get element $ElemId]
    #ElementInfo: <layer> <elemtype> <NumNodes> <N1> <N2> ...
    return "[lindex $ElementInfo 3] [lindex $ElementInfo 4] [lindex $ElementInfo 5]\
    [lindex $ElementInfo 6] [lindex $ElementInfo 7] [lindex $ElementInfo 8]\
    [lindex $ElementInfo 9] [lindex $ElementInfo 10] [lindex $ElementInfo 11]\
    [lindex $ElementInfo 12] [lindex $ElementInfo 13] [lindex $ElementInfo 14]\
    [lindex $ElementInfo 15] [lindex $ElementInfo 16] [lindex $ElementInfo 17]\
    [lindex $ElementInfo 18] [lindex $ElementInfo 19] [lindex $ElementInfo 20]\
    [lindex $ElementInfo 21] [lindex $ElementInfo 22] [lindex $ElementInfo 23]\
    [lindex $ElementInfo 24] [lindex $ElementInfo 25] [lindex $ElementInfo 26]\
    [lindex $ElementInfo 27] [lindex $ElementInfo 28] [lindex $ElementInfo 29]"
}

#-------------------------------------------------------------------------------

proc TriangleInterface2D4Connectivities { ElemId } {
    
    # Obtaining element nodes
    set ElementInfo [GiD_Mesh get element $ElemId]
    #ElementInfo: <layer> <elemtype> <NumNodes> <N1> <N2> ...
    return "[lindex $ElementInfo 3] [lindex $ElementInfo 4] [lindex $ElementInfo 4]\
    [lindex $ElementInfo 5]"
}

#-------------------------------------------------------------------------------

proc QuadrilateralInterface2D4Connectivities { ElemId } {
        
    # Obtaining element nodes
    set ElementInfo [GiD_Mesh get element $ElemId]
    #ElementInfo: <layer> <elemtype> <NumNodes> <N1> <N2> ...
    return "[lindex $ElementInfo 3] [lindex $ElementInfo 4] [lindex $ElementInfo 5]\
    [lindex $ElementInfo 6]"
}

#-------------------------------------------------------------------------------

proc TetrahedronInterface3D6Connectivities { ElemId } {
        
    set ElementInfo [GiD_Mesh get element $ElemId]
    #ElementInfo: <layer> <elemtype> <NumNodes> <N1> <N2> ...
    return "[lindex $ElementInfo 3] [lindex $ElementInfo 4] [lindex $ElementInfo 5]\
    [lindex $ElementInfo 6] [lindex $ElementInfo 4] [lindex $ElementInfo 5]"
}

#-------------------------------------------------------------------------------

proc PrismInterface3D6Connectivities { ElemId } {

    # Obtaining element nodes
    set ElementInfo [GiD_Mesh get element $ElemId]
    #ElementInfo: <layer> <elemtype> <NumNodes> <N1> <N2> ...
    return "[lindex $ElementInfo 3] [lindex $ElementInfo 4] [lindex $ElementInfo 5]\
    [lindex $ElementInfo 6] [lindex $ElementInfo 7] [lindex $ElementInfo 8]"
}

#-------------------------------------------------------------------------------

proc HexahedronInterface3D8Connectivities { ElemId } {

    # Obtaining element nodes
    set ElementInfo [GiD_Mesh get element $ElemId]
    #ElementInfo: <layer> <elemtype> <NumNodes> <N1> <N2> ...
    return "[lindex $ElementInfo 3] [lindex $ElementInfo 4] [lindex $ElementInfo 5]\
    [lindex $ElementInfo 6] [lindex $ElementInfo 7] [lindex $ElementInfo 8]\
    [lindex $ElementInfo 9] [lindex $ElementInfo 10]"
}

#-------------------------------------------------------------------------------

proc Line2D2Connectivities { ElemId } {
    
    set ElementInfo [GiD_Mesh get element $ElemId]
    #ElementInfo: <layer> <elemtype> <NumNodes> <N1> <N2> ...
    return "[lindex $ElementInfo 3] [lindex $ElementInfo 4]"
}

#-------------------------------------------------------------------------------

proc TriangleInterface3D4Connectivities { ElemId } {
    
    # Obtaining element nodes
    set ElementInfo [GiD_Mesh get element $ElemId]
    #ElementInfo: <layer> <elemtype> <NumNodes> <N1> <N2> ...
    return "[lindex $ElementInfo 3] [lindex $ElementInfo 4] [lindex $ElementInfo 4]\
    [lindex $ElementInfo 5]"
}

#-------------------------------------------------------------------------------

proc QuadrilateralInterface3D4Connectivities { ElemId } {
    
    # Obtaining element nodes
    set ElementInfo [GiD_Mesh get element $ElemId]
    #ElementInfo: <layer> <elemtype> <NumNodes> <N1> <N2> ...
    return "[lindex $ElementInfo 3] [lindex $ElementInfo 4] [lindex $ElementInfo 5]\
    [lindex $ElementInfo 6]"
}

#--------------------------------------------------------------------------------------------------------------------------------------------------------------

proc WriteElementSubmodelPart {FileVar CondName} {
    set Groups [GiD_Info conditions $CondName groups]
    if {[llength $Groups]>0} {
        upvar $FileVar MyFileVar
        
        for {set i 0} {$i < [llength $Groups]} {incr i} {
            puts $MyFileVar "Begin SubModelPart [lindex [lindex $Groups $i] 1]"
            # Nodes
            set Entities [GiD_EntitiesGroups get [lindex [lindex $Groups $i] 1] nodes]
            puts $MyFileVar "  Begin SubModelPartNodes"
            for {set j 0} {$j < [llength $Entities]} {incr j} {
                puts $MyFileVar "    [lindex $Entities $j]"
            }
            puts $MyFileVar "  End SubModelPartNodes"
            # Elements
            set Entities [GiD_EntitiesGroups get [lindex [lindex $Groups $i] 1] elements]
            puts $MyFileVar "  Begin SubModelPartElements"
            for {set j 0} {$j < [llength $Entities]} {incr j} {
                puts $MyFileVar "    [lindex $Entities $j]"
            }
            puts $MyFileVar "  End SubModelPartElements"
            # Conditions
            puts $MyFileVar "  Begin SubModelPartConditions"
            puts $MyFileVar "  End SubModelPartConditions"
            puts $MyFileVar "End SubModelPart"
            puts $MyFileVar ""
        }
    }    
}

#-------------------------------------------------------------------------------

proc WriteConstraintSubmodelPart {FileVar CondName} {
    set Groups [GiD_Info conditions $CondName groups]
    if {[llength $Groups]>0} {
        upvar $FileVar MyFileVar
        
        for {set i 0} {$i < [llength $Groups]} {incr i} {
            puts $MyFileVar "Begin SubModelPart [lindex [lindex $Groups $i] 1]"
            # Nodes
            set Entities [GiD_EntitiesGroups get [lindex [lindex $Groups $i] 1] nodes]
            puts $MyFileVar "  Begin SubModelPartNodes"
            for {set j 0} {$j < [llength $Entities]} {incr j} {
                puts $MyFileVar "    [lindex $Entities $j]"
            }
            puts $MyFileVar "  End SubModelPartNodes"
            # Elements
            puts $MyFileVar "  Begin SubModelPartElements"
            puts $MyFileVar "  End SubModelPartElements"
            # Conditions
            puts $MyFileVar "  Begin SubModelPartConditions"
            puts $MyFileVar "  End SubModelPartConditions"
            puts $MyFileVar "End SubModelPart"
            puts $MyFileVar ""
        }
    }
}

#-------------------------------------------------------------------------------

proc WriteLoadSubmodelPart {FileVar CondName ConditionDict} {
    set Groups [GiD_Info conditions $CondName groups]
    if {[llength $Groups]>0} {
        upvar $FileVar MyFileVar
        
        for {set i 0} {$i < [llength $Groups]} {incr i} {
            puts $MyFileVar "Begin SubModelPart [lindex [lindex $Groups $i] 1]"
            # Nodes
            set Entities [GiD_EntitiesGroups get [lindex [lindex $Groups $i] 1] nodes]
            puts $MyFileVar "  Begin SubModelPartNodes"
            for {set j 0} {$j < [llength $Entities]} {incr j} {
                puts $MyFileVar "    [lindex $Entities $j]"
            }
            puts $MyFileVar "  End SubModelPartNodes"
            # Elements
            puts $MyFileVar "  Begin SubModelPartElements"
            puts $MyFileVar "  End SubModelPartElements"
            # Conditions
            set ConditionList [dict get $ConditionDict [lindex [lindex $Groups $i] 1]]
            puts $MyFileVar "  Begin SubModelPartConditions"
            for {set j 0} {$j < [llength $ConditionList]} {incr j} {
                puts $MyFileVar "    [lindex $ConditionList $j]"
            }
            puts $MyFileVar "  End SubModelPartConditions"
            puts $MyFileVar "End SubModelPart"
            puts $MyFileVar ""
        }
    }
}

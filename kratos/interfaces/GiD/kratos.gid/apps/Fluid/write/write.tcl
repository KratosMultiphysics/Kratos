namespace eval Fluid::write {
    # Namespace variables declaration
    variable FluidConditions

    variable PartsUN
    variable BCUN
    variable writeCoordinatesByGroups

    variable validApps
}

proc Fluid::write::Init { } {
    # Namespace variables inicialization
    variable FluidConditions
    set FluidConditions(temp) 0
    unset FluidConditions(temp)

    variable PartsUN
    set PartsUN "FLParts"
    variable BCUN
    set BCUN "FLBC"

    variable writeCoordinatesByGroups
    set writeCoordinatesByGroups 0

    variable validApps
    set validApps [list "Fluid"]
}

proc Fluid::write::AddValidApps {appid} {
    variable validApps
    if {$appid ni $validApps} {lappend validApps $appid}
}

proc Fluid::write::SetCoordinatesByGroups {value} {
    variable writeCoordinatesByGroups
    set writeCoordinatesByGroups $value
}

# Events
proc Fluid::write::writeModelPartEvent { } {
    variable PartsUN
    variable writeCoordinatesByGroups
    variable validApps
    set err [Validate]
    if {$err ne ""} {error $err}
    write::initWriteData $PartsUN "FLMaterials"
    write::writeModelPartData
    writeProperties
    write::writeMaterials $validApps
    if {$writeCoordinatesByGroups} {write::writeNodalCoordinatesOnParts} {write::writeNodalCoordinates}
    write::writeElementConnectivities
    writeConditions
    writeMeshes
}
proc Fluid::write::writeCustomFilesEvent { } {
    write::CopyFileIntoModel "python/KratosFluid.py"
    write::RenameFileInModel "KratosFluid.py" "MainKratos.py"
}

proc Fluid::write::Validate {} {
    variable PartsUN
    set err ""
    set doc $gid_groups_conds::doc
    set root [$doc documentElement]

    # Check only 1 part in Parts
    set xp1 "[spdAux::getRoute $PartsUN]/group"
    if {[llength [$root selectNodes $xp1]] ne 1} {
        set err "You must set one part.\n"
    }
    # Check closed volume
    #if {[CheckClosedVolume] ne 1} {
    #    append err "Check boundary conditions."
    #}
    return $err
}

# MDPA Blocks
proc Fluid::write::writeProperties { } {
    # Begin Properties
    write::WriteString "Begin Properties 0"
    write::WriteString "End Properties"
    write::WriteString ""
}

proc Fluid::write::writeConditions { } {
    writeBoundaryConditions
    writeDrags
}

proc Fluid::write::writeBoundaryConditions { } {
    variable FluidConditions
    variable BCUN

    # Write the conditions
    set dict_group_intervals [write::writeConditions $BCUN]

    # Vamos a construir el array que nos permite escribir submodel parts y la malla de condiciones de contorno
    set doc $gid_groups_conds::doc
    set root [$doc documentElement]
    set xp1 "[spdAux::getRoute $BCUN]/condition/group"
    set iter 1
    foreach group [$root selectNodes $xp1] {
        set condid [[$group parent] @n]
        set groupid [get_domnode_attribute $group n]
        set cond [::Model::getCondition $condid]
        if {[$cond getAttribute SkinConditions]} {
            lassign [dict get $dict_group_intervals $groupid] ini fin
            set FluidConditions($groupid,initial) $ini
            set FluidConditions($groupid,final) $fin
            set FluidConditions($groupid,SkinCondition) 1
            #W "ARRAY [array get FluidConditions]"
        } else {
            set FluidConditions($groupid,initial) -1
            set FluidConditions($groupid,final) -1
            set FluidConditions($groupid,SkinCondition) 0
        }
    }
}

proc Fluid::write::writeDrags { } {

}

proc Fluid::write::writeMeshes { } {
    write::writePartMeshes
    write::writeNodalConditions "FLNodalConditions"
    writeConditionsMesh
    #writeSkinMesh
}

proc Fluid::write::writeConditionsMesh { } {
    variable FluidConditions
    variable BCUN
    set doc $gid_groups_conds::doc
    set root [$doc documentElement]
    set xp1 "[spdAux::getRoute $BCUN]/condition/group"
    #W "Conditions $xp1 [$root selectNodes $xp1]"
    foreach group [$root selectNodes $xp1] {
        set groupid [$group @n]
        set condid [[$group parent] @n]
        set ini $FluidConditions($groupid,initial)
        set end $FluidConditions($groupid,final)
        #W "$groupid $ini $end"
        if {$ini == -1} {
            ::write::writeGroupMesh $condid $groupid "Nodes"
        } else {
            ::write::writeGroupMesh $condid $groupid "Conditions" [list $ini $end]
        }
    }
}

proc Fluid::write::writeSkinMesh { } {
    variable FluidConditions
    variable BCUN
    set doc $gid_groups_conds::doc
    set root [$doc documentElement]
    set xp1 "[spdAux::getRoute $BCUN]/condition/group"
    #W "Conditions $xp1 [$root selectNodes $xp1]"
    set listiniend [list ]
    set listgroups [list ]
    foreach group [$root selectNodes $xp1] {
        set groupid [$group @n]
        set ini $FluidConditions($groupid,initial)
        set end $FluidConditions($groupid,final)
        lappend listiniend $ini $end
        lappend listgroups $groupid
    }
    set skinconfgroup "SKINCONDITIONS"
    catch {GiD_Groups delete $skinconfgroup}
    GiD_Groups create $skinconfgroup
    GiD_Groups edit state $skinconfgroup hidden
    foreach group $listgroups {
        GiD_EntitiesGroups assign $skinconfgroup nodes [GiD_EntitiesGroups get $group nodes]
    }
    ::write::writeGroupMesh EXTRA $skinconfgroup "Conditions" $listiniend
}

proc Fluid::write::CheckClosedVolume {} {
    variable BCUN
    set isclosed 1

    set doc $gid_groups_conds::doc
    set root [$doc documentElement]
    set xp1 "[spdAux::getRoute $BCUN]/condition/group"

    set listgroups [list ]
    foreach group [$root selectNodes $xp1] {
        set groupid [$group @n]
        set conditionName [[$group parent] @n]
        set cond [::Model::getCondition $conditionName]
        if {[$cond getAttribute "SkinConditions"] eq "True"} {
            set surfaces [GiD_EntitiesGroups get $groupid surfaces]
            foreach surf $surfaces {
                set linesraw [GiD_Geometry get surface $surf]
                set nlines [lindex $linesraw 2]
                set linespairs [lrange $linesraw 9 [expr 8 + $nlines]]
                foreach pair $linespairs {
                    set lid [lindex $pair 0]
                    incr usedsurfaceslines($lid)
                }
            }
        }
    }
    foreach lid [array names usedsurfaceslines] {
        if {$usedsurfaceslines($lid) ne "2"} {set isclosed 0;}
    }
    return $isclosed
}

Fluid::write::Init

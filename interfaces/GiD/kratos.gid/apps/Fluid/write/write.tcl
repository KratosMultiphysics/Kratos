namespace eval Fluid::write {
    # Namespace variables declaration
    variable FluidConditions
}

proc Fluid::write::Init { } {
    # Namespace variables inicialization
    variable FluidConditions
    set FluidConditions(temp) 0
    unset FluidConditions(temp)
}

# Events
proc Fluid::write::writeModelPartEvent { } {
    write::initWriteData "FLParts" "FLMaterials"
    write::writeModelPartData
    writeProperties
    write::writeMaterials
    write::writeNodalCoordinates
    write::writeElementConnectivities
    writeConditions
    writeMeshes
}
proc Fluid::write::writeCustomFilesEvent { } {

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
    
    # Write the conditions
    set dict_group_intervals [write::writeConditions "FLBC"]
    
    # Vamos a construir el array que nos permite escribir submodel parts y la malla de condiciones de contorno
    set doc $gid_groups_conds::doc
    set root [$doc documentElement]
    set xp1 "[apps::getRoute "FLBC"]/condition/group"
    set iter 1
    foreach group [$root selectNodes $xp1] {
        set condid [[$group parent] @n]
        set groupid [get_domnode_attribute $group n]
        set cond [::Model::getCondition $condid]
        
        lassign [dict get $dict_group_intervals $groupid] ini fin
        set FluidConditions($groupid,initial) $ini
        set FluidConditions($groupid,final) $fin
        set bc 0
        if {[$cond getAttribute SkinConditions]} {set bc 1}
        set FluidConditions($groupid,SkinCondition) $bc
        #W "ARRAY [array get FluidConditions]"
    }
}

proc Fluid::write::writeDrags { } {
    
}

proc Fluid::write::writeMeshes { } {
    write::writePartMeshes
    write::writeNodalConditions "FLNodalConditions"
    writeConditionsMesh
    writeSkinMesh
}

proc Fluid::write::writeConditionsMesh { } {
    variable FluidConditions
    set doc $gid_groups_conds::doc
    set root [$doc documentElement]
    set xp1 "[apps::getRoute "FLBC"]/condition/group"
    #W "Conditions $xp1 [$root selectNodes $xp1]"
    foreach group [$root selectNodes $xp1] {
        set groupid [$group @n]
        set ini $FluidConditions($groupid,initial)
        set end $FluidConditions($groupid,final)
        #W "$groupid $ini $end"
        ::write::writeGroupMesh [[$group parent] @n] $groupid "Conditions" [list $ini $end]
    }
}

proc Fluid::write::writeSkinMesh { } {
    variable FluidConditions
    set doc $gid_groups_conds::doc
    set root [$doc documentElement]
    set xp1 "[apps::getRoute "FLBC"]/condition/group"
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
    set isclosed 1
    
    set doc $gid_groups_conds::doc
    set root [$doc documentElement]
    set xp1 "[apps::getRoute "FLBC"]/condition/group"

    set listgroups [list ]
    foreach group [$root selectNodes $xp1] {
        set groupid [$group @n]
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
    foreach lid [array names usedsurfaceslines] {
        if {$usedsurfaceslines($lid) ne "2"} {set isclosed 0;}
    }
    return $isclosed
}

Fluid::write::Init

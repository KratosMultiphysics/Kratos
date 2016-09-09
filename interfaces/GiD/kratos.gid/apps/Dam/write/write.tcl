namespace eval Dam::write {
    variable nodalwrite
    
    variable ConditionsDictGroupIterators
    variable NodalConditionsGroup
}

proc Dam::write::Init { } {
    # Namespace variables inicialization
    variable nodalwrite
    set nodalwrite 0
    
    variable ConditionsDictGroupIterators
    variable NodalConditionsGroup
    set ConditionsDictGroupIterators [dict create]
    set NodalConditionsGroup [list ]
}


proc Dam::write::writeCustomFilesEvent { } {
    
    write::CopyFileIntoModel "python/dam_thermo_mechanic_script.py"
    write::RenameFileInModel "dam_thermo_mechanic_script.py" "MainKratos.py"
    
    #write::RenameFileInModel "ProjectParameters.json" "ProjectParameters.py"
}

# MDPA Blocks

proc Dam::write::writeModelPartEvent { } {
    write::initWriteData "DamParts" "DamMaterials"
    
    write::writeModelPartData
    write::WriteString "Begin Properties 0"
    write::WriteString "End Properties"
    
    UpdateMaterials
    write::writeMaterials
    #write::writeTables
    write::writeNodalCoordinates
    write::writeElementConnectivities
    writeConditions
    writeMeshes
    #writeCustomBlock
}

proc Dam::write::UpdateMaterials { } {
    set matdict [write::getMatDict]
    foreach {mat props} $matdict {
        set constlaw [dict get $props ConstitutiveLaw]
        # Modificar la ley constitutiva
        set newconstlaw $constlaw
        if {$constlaw eq "BilinearCohesive2DPlaneStress"} {set newconstlaw "BilinearCohesive2DLaw"}
        if {$constlaw eq "BilinearCohesive2DPlaneStrain"} {
            dict set matdict $mat THICKNESS  1.0000E+00
            set newconstlaw "BilinearCohesive2DLaw"
            }
        dict set matdict $mat CONSTITUTIVE_LAW_NAME $newconstlaw
    }
    write::setMatDict $matdict
}

proc Dam::write::writeConditions { } {
    variable ConditionsDictGroupIterators
    set ConditionsDictGroupIterators [write::writeConditions "DamLoads"]
}

proc Dam::write::writeMeshes { } {
    
    write::writePartMeshes
    
    # Solo Malla , no en conditions
    write::writeNodalConditions "DamNodalConditions"
    
    # A Condition y a meshes-> salvo lo que no tenga topologia
    writeLoads
}

proc Dam::write::writeLoads { } {
    #writeGravity "SLGravity"
    variable nodalwrite
    variable ConditionsDictGroupIterators
    set doc $gid_groups_conds::doc
    set root [$doc documentElement]
    set xp1 "[spdAux::getRoute "DamLoads"]/condition/group"
    foreach group [$root selectNodes $xp1] {
        set groupid [$group @n]
        #W "Writing mesh of Load $groupid"
        if {$groupid in [dict keys $ConditionsDictGroupIterators]} {
            ::write::writeGroupMesh [[$group parent] @n] $groupid "Conditions" [dict get $ConditionsDictGroupIterators $groupid]
        } else {
            ::write::writeGroupMesh [[$group parent] @n] $groupid "nodal"
        }
    }
}

Dam::write::Init

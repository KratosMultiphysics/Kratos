namespace eval Solid::write {
    variable validApps
    variable ConditionsDictGroupIterators
    variable NodalConditionsGroup
    variable writeCoordinatesByGroups
}

proc Solid::write::Init { } {
    # Namespace variables inicialization
    
    variable ConditionsDictGroupIterators
    variable NodalConditionsGroup
    set ConditionsDictGroupIterators [dict create]
    set NodalConditionsGroup [list ]
    
    variable validApps
    set validApps [list "Solid"]
    
    variable writeCoordinatesByGroups
    set writeCoordinatesByGroups 0
}

proc Solid::write::AddValidApps {appList} {
    variable validApps
    set validApps [list "Solid"]
    lappend validApps $appList
}

proc Solid::write::writeCustomFilesEvent { } {
    WriteMaterialsFile
    
    write::CopyFileIntoModel "python/KratosSolid.py"
    write::RenameFileInModel "KratosSolid.py" "MainKratos.py"
    
    #write::RenameFileInModel "ProjectParameters.json" "ProjectParameters.py"
}

proc Solid::write::SetCoordinatesByGroups {value} {
    variable writeCoordinatesByGroups
    set writeCoordinatesByGroups $value
}

# MDPA Blocks

proc Solid::write::writeModelPartEvent { } {
    variable writeCoordinatesByGroups
    variable validApps
    write::initWriteData "SLParts" "SLMaterials"
    
    write::writeModelPartData
    write::WriteString "Begin Properties 0"
    write::WriteString "End Properties"
    write::writeMaterials $validApps
    #write::writeTables
    if {$writeCoordinatesByGroups} {write::writeNodalCoordinatesOnParts} {write::writeNodalCoordinates}
    write::writeElementConnectivities
    writeConditions
    writeMeshes
    #writeCustomBlock
}


proc Solid::write::writeConditions { } {
    variable ConditionsDictGroupIterators
    set ConditionsDictGroupIterators [write::writeConditions "SLLoads"]
}

proc Solid::write::writeMeshes { } {
    
    write::writePartMeshes
    
    # Solo Malla , no en conditions
    write::writeNodalConditions "SLNodalConditions"
    
    # A Condition y a meshes-> salvo lo que no tenga topologia
    writeLoads
}

proc Solid::write::writeLoads { } {
    #writeGravity "SLGravity"
    variable ConditionsDictGroupIterators
    set doc $gid_groups_conds::doc
    set root [$doc documentElement]
    set xp1 "[spdAux::getRoute "SLLoads"]/condition/group"
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


proc Solid::write::writeCustomBlock { } {
    write::WriteString "Begin Custom"
    write::WriteString "Custom write for Solid, any app can call me, so be careful!"
    write::WriteString "End Custom"
    write::WriteString ""
}

# Custom files
proc Solid::write::WriteMaterialsFile { } {
    variable validApps
    
    write::OpenFile "materials.py"
    
    set str "
from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# Importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.SolidMechanicsApplication import *
#from beam_sections_python_utility import SetProperties
#from beam_sections_python_utility import SetMaterialProperties

def AssignMaterial(Properties):
    # material for solid material
"
    foreach {part mat} [write::getMatDict] {
        if {[dict get $mat APPID] in $validApps} {
            append str "
    prop_id = [dict get $mat MID];
    prop = Properties\[prop_id\]
    mat = [dict get $mat ConstitutiveLaw]()
    prop.SetValue(CONSTITUTIVE_LAW, mat.Clone())
        "
            write::WriteString $str
        }
    }
    write::CloseFile
    
}

Solid::write::Init

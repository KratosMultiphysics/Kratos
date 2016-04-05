namespace eval Solid::write {
    variable nodalwrite
    
    variable ConditionsDictGroupIterators
    variable NodalConditionsGroup
}

proc Solid::write::Init { } {
    # Namespace variables inicialization
    variable nodalwrite
    set nodalwrite 0
    
    variable ConditionsDictGroupIterators
    variable NodalConditionsGroup
    set ConditionsDictGroupIterators [dict create]
    set NodalConditionsGroup [list ]
}


proc Solid::write::writeCustomFilesEvent { } {
    WriteMaterialsFile
    
    write::CopyFileIntoModel "python/KratosSolid.py"
    write::RenameFileInModel "KratosSolid.py" "MainKratos.py"
    
    #write::RenameFileInModel "ProjectParameters.json" "ProjectParameters.py"
}

# MDPA Blocks

proc Solid::write::writeModelPartEvent { } {
    write::initWriteData "SLParts" "SLMaterials"
    
    write::writeModelPartData
    write::WriteString "Begin Properties 0"
    write::WriteString "End Properties"
    write::writeMaterials
    #write::writeTables
    write::writeNodalCoordinates
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
    variable nodalwrite
    variable ConditionsDictGroupIterators
    set doc $gid_groups_conds::doc
    set root [$doc documentElement]
    set xp1 "[apps::getRoute "SLLoads"]/condition/group"
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
    # Materials.py
    
    write::OpenFile "materials.py"
    
    set str "
from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# Importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.SolidApplication import *
#from beam_sections_python_utility import SetProperties
#from beam_sections_python_utility import SetMaterialProperties

def AssignMaterial(Properties):
    # material for solid material
"
    foreach {part mat} [write::getMatDict] {
        append str "
    prop_id = [dict get $mat MID];
    prop = Properties\[prop_id\]
    mat = [dict get $mat ConstitutiveLaw]()
    prop.SetValue(CONSTITUTIVE_LAW, mat.Clone())
        "
    }
    write::WriteString $str
    
    write::CloseFile
    
}

Solid::write::Init

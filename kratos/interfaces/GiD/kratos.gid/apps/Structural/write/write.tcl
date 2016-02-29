namespace eval Structural::write {
    variable nodalwrite
    
    variable ConditionsDictGroupIterators
    variable NodalConditionsGroup
}

proc Structural::write::Init { } {
    # Namespace variables inicialization
    variable nodalwrite
    set nodalwrite 0
    
    variable ConditionsDictGroupIterators
    variable NodalConditionsGroup
    set ConditionsDictGroupIterators [dict create]
    set NodalConditionsGroup [list ]
}

# Project Parameters
proc Structural::write::writeParametersEvent { } {
    #write::WriteString "Project Parameters"
    #cortar a nombre
    write::WriteString "ProblemName =\"[file tail [GiD_Info Project ModelName]]\""
    write::WriteString ""
    write::WriteString "#Problem Data"
    write::WriteString "#################################################"
    write::WriteString "ModelPartName = \"Structure\""
    write::WriteString "DomainSize = [string range [write::getValue nDim] 0 0]"
    
    # Parallelization
    set paralleltype [write::getValue SMParallelization ParallelSolutionType]
    if {$paralleltype eq "OpenMP"} {
        set nthreads [write::getValue SMParallelization OpenMPNumberOfThreads]
        write::WriteString "NumberofThreads = $nthreads"
    } else {
        set nthreads [write::getValue SMParallelization MPINumberOfProcessors]
        write::WriteString "NumberofProcessors = $nthreads"
    }
    
    # Time Parameters
    write::WriteString "time_step = [write::getValue SMTimeParameters DeltaTime]"
    write::WriteString "end_time = [write::getValue SMTimeParameters EndTime]"
    
    write::WriteString "EchoLevel = [write::getValue SMResults EchoLevel]"
    
    # Solution strategy
    write::WriteString "#Solver Data"
    write::WriteString "#################################################"
    write::WriteString "class SolverSettings:"
    write::WriteString "    solver_type = \"solid_mechanics_python_solver\""
    write::WriteString "    domain_size = DomainSize"
    write::WriteString "    echo_level = EchoLevel"
    write::WriteString "    solution_type  = \"[write::getValue SMSoluType]\""
    write::WriteString "    "
    
    printSolutionStrategy 4
    printSolvers 4
    
    printConstraints 0
    printLoads 0
    
    printResults 0
    
    # GiD output configuration
    write::WriteString "class GidOutputConfiguration:"
    write::WriteString "    GiDPostMode = \"[write::getValue SMResults GiDPostMode]\""
    write::WriteString "    GiDWriteMeshFlag = [write::getStringBinaryValue SMResults GiDWriteMeshFlag]"
    write::WriteString "    GiDWriteConditionsFlag = [write::getStringBinaryValue SMResults GiDWriteConditionsFlag]"
    write::WriteString "    GiDWriteParticlesFlag = [write::getStringBinaryValue SMResults GiDWriteParticlesFlag]"
    write::WriteString "    GiDMultiFileFlag = \"[write::getValue SMResults GiDMultiFileFlag]\""

    write::WriteString ""   
    write::WriteString "GiDWriteFrequency = [write::getValue SMResults OutputDeltaTime]"
    write::WriteString ""
    write::WriteString "# graph_options"
    write::WriteString "PlotGraphs = \"False\""
    write::WriteString "PlotFrequency = 0"
    write::WriteString ""
    write::WriteString "# list options"
    write::WriteString "PrintLists = \"True\""
    write::WriteString "file_list = \[\] "
    write::WriteString ""
    write::WriteString "# restart options"
    write::WriteString "SaveRestart = False"
    write::WriteString "RestartFrequency = 0"
    write::WriteString "LoadRestart = False"
    write::WriteString "Restart_Step = 0"
    write::WriteString ""

    
}

proc Structural::write::printResults {spacing} {
    set doc $gid_groups_conds::doc
    set root [$doc documentElement]
    set s [write::getSpacing $spacing]

    write::WriteString "#PostProcess Data"
    write::WriteString "#################################################"
    
    set xp1 "[apps::getRoute "SMNodalResults"]/value"
    set results [$root selectNodes $xp1]
    set min 0
    set str "nodal_results=\["
    foreach res $results {
        if {[get_domnode_attribute $res v] eq "Yes"} {
            set min 1
            set name [get_domnode_attribute $res n]
            append str "\"$name\","
        }
    }
    
    if {$min} {set str [string range $str 0 end-1]}
    append str "]"
    write::WriteString $str
    
    set xp1 "[apps::getRoute "SMElementResults"]/value"
    set results [$root selectNodes $xp1]
    set min 0
    set str "gauss_points_results=\["
    foreach res $results {
        if {[get_domnode_attribute $res v] in [list "Yes" "True"] && [get_domnode_attribute $res state] eq "normal"} {
            set min 1
            set name [get_domnode_attribute $res n]
            append str "\"$name\","
        }
    }
    
    if {$min} {set str [string range $str 0 end-1]}
    append str "]"
    write::WriteString $str
    write::WriteString ""
}

proc Structural::write::printConstraints {spacing} {
    set doc $gid_groups_conds::doc
    set root [$doc documentElement]
    set s [write::getSpacing $spacing]

    write::WriteString "#Constraints Data"
    write::WriteString "#################################################"
    write::WriteString "${s}constraints_process_list = \["
    set xp1 "[apps::getRoute "SMDoFs"]/condition/group"
    set groups [$root selectNodes $xp1]
    set str ""
    foreach group $groups {
        set groupName [$group @n]
        set cid [[$group parent] @n]
        set groupId [::write::getMeshId $cid $groupName]
        #W [[$group parent] @type]
        if {[[$group parent] @type] eq "vector"} {
            set FixX [expr [get_domnode_attribute [$group find n FixX] v] ? True : False]
            set FixY [expr [get_domnode_attribute [$group find n FixY] v] ? True : False]
            set FixZ [expr [get_domnode_attribute [$group find n FixZ] v] ? True : False]
            set ValX [get_domnode_attribute [$group find n ValX] v] 
            set ValY [get_domnode_attribute [$group find n ValY] v] 
            set ValZ [get_domnode_attribute [$group find n ValZ] v]
            set factornode [$group find n FACTOR]
            set Factor 1
            if {$factornode ne ""} {
                set Factor [get_domnode_attribute $factornode v]
            }
	    set arguments [list "a" "b"]
            write::WriteProcess "ApplyConstantVectorValueProcess" $arguments
            append str "{ \"process_name\" : \"ApplyConstantVectorValueProcess\",
              \"implemented_in_module\" : \"KratosMultiphysics\",
              \"implemented_in_file\": \"process_factory\",
              
              \"parameters\" : { 
                  \"mesh_id\": $groupId,
                  \"model_part_name\" : ModelPartName,
                  \"variable_name\": \"[[$group parent] @n]\",
                  \"factor\": $Factor,
                  \"value\": \[$ValX,$ValY,$ValZ\],
                  \"is_fixed_x\" : $FixX,
                  \"is_fixed_y\" : $FixY,
                  \"is_fixed_z\" : $FixZ,
                  }
            }, "
        } {
            set Value [get_domnode_attribute [$group firstChild] v] 
            set Fixed True
            
            append str "{ \"process_name\" : \"ApplyConstantScalarValueProcess\",
            \"implemented_in_python\" : True,
            \"implemented_in_module\" : \"KratosMultiphysics\",
            \"implemented_in_file\": \"process_factory\",
                        
            \"parameters\" : { 
                \"mesh_id\": $groupId,
                \"model_part_name\" : ModelPartName,
                \"variable_name\": \"$processName\",
                \"value\": $Value,
                \"is_fixed\": $Fixed
                }
          }, "
        }
    }
    # Quitar la coma
    if {[llength $groups]} {set str [string range $str 0 end-1]}
    write::WriteString $str
    write::WriteString "]"
    
}

proc Structural::write::printLoads {spacing} {
    set doc $gid_groups_conds::doc
    set root [$doc documentElement]
    set s [write::getSpacing $spacing]

    write::WriteString "#Loads Data"
    write::WriteString "#################################################"
    write::WriteString "${s}loads_process_list = \["
    set xp1 "[apps::getRoute "SMLoads"]/condition/group"
    set groups [$root selectNodes $xp1]
    set str ""
    foreach group $groups {
        set groupName [$group @n]
        set cid [[$group parent] @n]
        set groupId [::write::getMeshId $cid $groupName]
        set condId [[$group parent] @n]
        set condition [::Model::getCondition $condId]
        set type [$condition getProcessFormat]
        set processName [$condition getProcessName]
        if {$type eq "vector"} {
            set FixX False
            set FixY False
            set FixZ False
            set ValX [get_domnode_attribute [$group find n ValX] v] 
            set ValY [get_domnode_attribute [$group find n ValY] v] 
            set ValZ [get_domnode_attribute [$group find n ValZ] v] 
            set factornode [$group find n FACTOR]
            set Factor 1
            if {$factornode ne ""} {
                set Factor [get_domnode_attribute $factornode v]
            }
            
            append str "{ \"process_name\" : \"ApplyConstantVectorValueProcess\",
              \"implemented_in_module\" : \"KratosMultiphysics\",
              \"implemented_in_file\": \"process_factory\",
              
              \"parameters\" : { 
                \"mesh_id\": $groupId,
                  \"model_part_name\" : ModelPartName,
                  \"variable_name\": \"$processName\",
                  \"factor\": $Factor,
                  \"value\": \[$ValX,$ValY,$ValZ\],
                  \"is_fixed_x\" : $FixX,
                  \"is_fixed_y\" : $FixY,
                  \"is_fixed_z\" : $FixZ,
                  }
            }, "
        } {
            set Value [get_domnode_attribute [$group firstChild] v] 
            set Fixed True
            
            append str "{ \"process_name\" : \"ApplyConstantScalarValueProcess\",
            \"implemented_in_python\" : True,
            \"implemented_in_module\" : \"KratosMultiphysics\",
            \"implemented_in_file\": \"process_factory\",
                        
            \"parameters\" : { 
				\"mesh_id\": $groupId,
                \"model_part_name\" : ModelPartName,
                \"variable_name\": \"$processName\",
                \"value\": $Value,
                \"is_fixed\": $Fixed
                }
          }, "
        }
    }
    # Quitar la coma
    if {[llength $groups]} {set str [string range $str 0 end-1]}
    write::WriteString $str
    write::WriteString "]"
}

proc Structural::write::printSolutionStrategy {spacing} {
    set solstratName [write::getValue SMSolStrat]
    set schemeName [write::getValue SMScheme]
    set sol [::Model::GetSolutionStrategy $solstratName]
    set sch [$sol getScheme $schemeName]
    
    set spaces [write::getSpacing $spacing]
    write::WriteString "${spaces}scheme_type = \"[$sch getName]\""
    write::WriteString "${spaces}RotationDofs = False"
    write::WriteString "${spaces}PressureDofs = False"
    
    write::WriteString "${spaces}time_integration_method = \"[$sol getName]\""
    foreach {n in} [$sol getInputs] {
        write::WriteString "$spaces$n = [write::getValue SMStratParams $n ]"
    }
    foreach {n in} [$sch getInputs] {
        write::WriteString "$spaces$n = [write::getValue SMStratParams $n ]"
    }
}
proc Structural::write::printSolvers {spacing} {
    set solstratName [write::getValue SMSolStrat]
    set sol [::Model::GetSolutionStrategy $solstratName]
    set spaces [write::getSpacing $spacing]
    set spaces2 [write::getSpacing [expr $spacing +4]]
    foreach se [$sol getSolversEntries] {
        set un "SM$solstratName[$se getName]"
        set solverName [write::getValue $un Solver]
        write::WriteString "${spaces}class [$se getName]:"
        write::WriteString "${spaces2}solver_type = \"$solverName\""
          
        foreach {n in} [[::Model::GetSolver $solverName] getInputs] {
            write::WriteString "$spaces2$n = [write::getValue $un $n ]"
        }
    }
}

proc Structural::write::writeCustomFilesEvent { } {
    WriteMaterialsFile
    
    write::CopyFileIntoModel "python/KratosStructural.py"
}



# MDPA Blocks

proc Structural::write::writeModelPartEvent { } {
    write::initWriteData "SMParts" "SMMaterials"
    write::setGroupsTypeName "Mesh"
    
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


proc Structural::write::writeConditions { } {
    variable ConditionsDictGroupIterators
    set ConditionsDictGroupIterators [write::writeConditions "SMLoads"]
}

proc Structural::write::writeMeshes { } {
    
    write::writePartMeshes
    
    # Solo Malla , no en conditions
    write::writeNodalConditions "SMDoFs"
    
    # A Condition y a meshes-> salvo lo que no tenga topologia
    writeLoads
}

proc Structural::write::writeLoads { } {
    variable ConditionsDictGroupIterators
    set doc $gid_groups_conds::doc
    set root [$doc documentElement]
    set xp1 "[apps::getRoute "SMLoads"]/condition/group"
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

proc Structural::write::writeCustomBlock { } {
    write::WriteString "Begin Custom"
    write::WriteString "Custom write for Structural, any app can call me, so be careful!"
    write::WriteString "End Custom"
    write::WriteString ""
}

# Custom files
proc Structural::write::WriteMaterialsFile { } {
    # Materials.py
    
    write::OpenFile "materials.py"
    
    set str "
from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# Importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.StructuralApplication import *
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

Structural::write::Init

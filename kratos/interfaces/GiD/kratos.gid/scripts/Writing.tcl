##################################################################################
#   This file is common for all Kratos Applications.
#   Do not change anything here unless it's strictly necessary.
##################################################################################

namespace eval write {
    variable mat_dict
    variable dir
    variable parts
    variable matun
    variable meshes
    variable groups_type_name
}

proc write::Init { } {
    variable mat_dict
    variable dir
    variable parts
    variable meshes
    variable groups_type_name
    
    set mat_dict ""
    set dir ""
    set parts ""
    set meshes [dict create]
    set groups_type_name "SubModelPart"
}

proc write::initWriteData {partes mats} {
    variable parts
    variable matun
    variable meshes
    set parts $partes
    set matun $mats
    
    set meshes [dict create]
    processMaterials
}

proc write::setGroupsTypeName {name} {
    variable groups_type_name
    set groups_type_name $name
}

# Write Events
proc write::writeEvent { filename } {
    #W $filename
    variable dir
    set dir [file dirname $filename]
    set errcode 0
    #set inittime [clock seconds]
    set activeapp [::apps::getActiveApp]
    set appid [::apps::getActiveAppId]
    
    #### MDPA Write ####
    set errcode [writeAppMDPA $appid]
	
    #### Project Parameters Write ####
    set wevent [$activeapp getWriteParametersEvent]
    set filename "ProjectParameters.json"
    
    catch {CloseFile}
    OpenFile $filename
    eval $wevent
    #if {$errcode eq 0 && [catch {eval $wevent} fid] } {
    #    W "Problem Writing Project Parameters block:\n$fid\nEnd problems"
    #    set errcode 1
    #}
    catch {CloseFile}
	
    #### Custom File Write ####
    set wevent [$activeapp getWriteCustomEvent]
    
    catch {CloseFile}
    if { [catch {eval $wevent} fid] } {
	W "Problem Writing Custom block:\n$fid\nEnd problems"
	set errcode 1
    }
    catch {CloseFile}
	
    return $errcode
    # set endtime [clock seconds]
    # set ttime [expr {$endtime-$inittime}]
    # W "Total time: [Duration $ttime]"
}

proc write::writeAppMDPA {appid} {
    set errcode 0
    set activeapp [::apps::getAppById $appid]
    
    #### MDPA Write ####
    set wevent [$activeapp getWriteModelPartEvent]
    set filename "[file tail [GiD_Info project ModelName]].mdpa"
    
    catch {CloseFile}
    OpenFile $filename
    eval $wevent
    # Delegate in app
    #if { [catch {eval $wevent} fid] } {
    #    W "Problem Writing MDPA block:\n$fid\nEnd problems"
    #    set errcode 1
    #}
    catch {CloseFile}
    return $errcode
}

proc write::writeModelPartData { } {
    # Write the model part data
     
	WriteString "Begin ModelPartData"
	WriteString "//  VARIABLE_NAME value"
	WriteString "End ModelPartData"
	WriteString ""
}

proc write::writeTables { } {
    # Write the model part data
     
	WriteString "Begin Table"
	WriteString "Table content"
	WriteString "End Tablee"
	WriteString ""
}

proc write::writeMaterials { } {
    variable mat_dict
    W "Hey"
    W $mat_dict
    set exclusionList [list "MID" "ConstitutiveLaw" "Material"]
    # We print all the material data directly from the saved dictionary
    foreach material [dict keys $mat_dict] {
	WriteString "Begin Properties [dict get $mat_dict $material MID]"
	foreach prop [dict keys [dict get $mat_dict $material] ] {
	    if {$prop ni $exclusionList} {
		WriteString "    $prop [dict get $mat_dict $material $prop] "
	    }
	}
	WriteString "End Properties"
	WriteString ""
    }

}

proc write::writeNodalCoordinates { } {
    # Write the nodal coordinates block
    # Nodes block format
    # Begin Nodes
    # // id          X        Y        Z 
    # End Nodes
    
    WriteString "Begin Nodes"
    write_calc_data coordinates "%5d %14.5f %14.5f %14.5f%.0s\n"
    WriteString "End Nodes"
    WriteString "\n"
}

proc write::processMaterials {  } {
    variable parts
    variable matun
    set doc $gid_groups_conds::doc
    set root [$doc documentElement]
    
    variable mat_dict
    set xp1 "[spdAux::getRoute $parts]/group"
    set xp2 ".//value\[@n='Material']"
    
    set mat_dict ""
    set material_number 0 
    
    foreach gNode [$root selectNodes $xp1] {
        set group [$gNode getAttribute n]
        #set valueNode [$gNode selectNodes $xp2]
        #set material_name [get_domnode_attribute $valueNode v] 
        set material_name "material $material_number"
        if { ![dict exists $mat_dict $group] } {            
            incr material_number
            set mid $material_number
            
            set xp3 [spdAux::getRoute $matun]
            append xp3 [format_xpath {/blockdata[@n="material" and @name=%s]/value} $material_name]
        
            dict set mat_dict $group MID $material_number 
            
            set s1 [$gNode selectNodes ".//value"]
            set s2 [$root selectNodes $xp3]
            set us [join [list $s1 $s2]]
            
            foreach valueNode $us {
                set name [$valueNode getAttribute n]
                set state [get_domnode_attribute $valueNode state]
                if {$state eq "normal" && $name ne "Element"} {
                    # All the introduced values are translated to 'm' and 'kg' with the help of this function
                    set value [gid_groups_conds::convert_value_to_default $valueNode]
                    
                    if {[string is double $value]} {
                        set value [format "%13.5E" $value]
                    }
                    dict set mat_dict $group $name $value
                }
            }
        } 
    }
}


proc write::writeElementConnectivities { } {
    variable parts
    set doc $gid_groups_conds::doc
    set root [$doc documentElement]
    variable mat_dict
    
    set xp1 "[spdAux::getRoute $parts]/group"
    set material_number 0
    foreach gNode [$root selectNodes $xp1] {
        set formats ""
        set group [get_domnode_attribute $gNode n]
        if { [dict exists $mat_dict $group] } {          
            set mid [dict get $mat_dict $group MID]
            if {[$gNode hasAttribute ov]} {set ov [get_domnode_attribute $gNode ov] } {set ov [get_domnode_attribute [$gNode parent] ov] }
            #W $ov
            lassign [getEtype $ov] etype nnodes
            #W "$group $ov -> $etype $nnodes"
            if {$nnodes ne ""} {
                set formats [GetFormatDict $group $mid $nnodes]
                if {$etype ne "none"} {
                    set kelemtype [get_domnode_attribute [$gNode selectNodes ".//value\[@n='Element']"] v]
                    set elem [::Model::getElement $kelemtype]
                    #W $kelemtype
                    set top [$elem getTopologyFeature $etype $nnodes]
                    if {$top eq ""} {W "Element $kelemtype not available for $ov entities on group $group"; continue}
                    set kratosElemName [$top getKratosName]
                    WriteString "Begin Elements $kratosElemName// GUI group identifier: $group" 
                    write_calc_data connectivities $formats
                    WriteString "End Elements"
                    WriteString ""     
                } 
            } else {
                error [= "Error on $group -  no known element type"]
            } 
        }
    } 
}

proc write::writeConditions { baseUN } {
    set dictGroupsIterators [dict create]
    set doc $gid_groups_conds::doc
    set root [$doc documentElement]
    
    set xp1 "[spdAux::getRoute $baseUN]/condition/group"
    set iter 1
    foreach group [$root selectNodes $xp1] {
	set condid [[$group parent] @n]
	set groupid [get_domnode_attribute $group n]
	set ov [[$group parent] @ov]
	set cond [::Model::getCondition $condid]
	lassign [write::getEtype $ov] etype nnodes
	set kname [$cond getTopologyKratosName $etype $nnodes]
	if {$kname ne ""} {
	    set obj [list ]
	    WriteString "Begin Conditions $kname// GUI group identifier: $groupid"
	    if {$etype eq "Point"} {
		set formats [dict create $groupid "%10d \n"]
		
		#set obj [write_calc_data nodes -return $formats]
		#set obj [list 3 5]
		set obj [GiD_EntitiesGroups get $groupid nodes]
	    } else {
		set formats [write::GetFormatDict $groupid 0 $nnodes]
		set elems [write_calc_data connectivities -return $formats]
		#W $elems
		
		set obj [GetListsOfNodes $elems $nnodes 2]
	    }
	    set initial $iter
	    for {set i 0} {$i <[llength $obj]} {incr iter; incr i} {
	    set nids [lindex $obj $i]
	    WriteString "$iter 0 $nids"
	    }
	    set final [expr $iter -1]
	    WriteString "End Conditions"
	    WriteString ""
	    dict set dictGroupsIterators $groupid [list $initial $final]
	}
    }
    return $dictGroupsIterators
}

proc write::GetListsOfNodes {elems nnodes {ignore 0} } {
    set obj [list ]
    set imax [llength $elems]
    if {$nnodes eq 0} {return $obj}
    set i 0
    while {$i < $imax} {
	for {set j 0} {$j < $ignore} {incr j} {incr i; if {$i >= $imax} {return $obj}}
	set tmp [list ]
	for {set j 0} {$j < $nnodes} {incr j} {
	    if {$i >= $imax} {break}
	    lappend tmp [lindex $elems $i]
	    incr i
	}
	lappend obj $tmp
    }
    return $obj
}

proc write::getMeshId {cid group} {
    variable meshes
    
    set find [list $cid ${group}]
    if {[dict exists $meshes $find]} {
        return [dict get $meshes [list $cid ${group}]]
    } {
        return 0
    }
}

proc write::writeGroupMesh { cid group {what "Elements"} {iniend ""} } {
    variable meshes
    variable groups_type_name
    W "$cid $group $what"
    set gtn $groups_type_name
    if {![dict exists $meshes [list $cid ${group}]]} {
	set mid [expr [llength [dict keys $meshes]] +1]
	if {$gtn ne "Mesh"} {
	    set underscoredGroup [string map {" " " "} $group]
	    regsub -all { +} $group "_" underscoredGroup
	    set mid "${cid}_${underscoredGroup}"
	}
	dict set meshes [list $cid ${group}] $mid
	set gdict [dict create]
	set f "%10i\n"
	set f [subst $f]
	dict set gdict $group $f
	WriteString "Begin $gtn $mid // Group $group // Subtree $cid"
	WriteString "    Begin ${gtn}Nodes"
	write_calc_data nodes -sorted $gdict
	WriteString "    End ${gtn}Nodes"
	WriteString "    Begin ${gtn}Elements"
	if {$what eq "Elements"} {
	    write_calc_data elements -sorted $gdict
	}
	WriteString "    End ${gtn}Elements"
	WriteString "    Begin ${gtn}Conditions"
	if {$what eq "Conditions"} {
	    #write_calc_data elements -sorted $gdict
	    if {$iniend ne ""} {
		#W $iniend
		foreach {ini end} $iniend {
		    for {set i $ini} {$i<=$end} {incr i} {
		        WriteString [format %10d $i]
		    }
		}
	    }
	}
	WriteString "    End ${gtn}Conditions"
	WriteString "End $gtn"
    }
}

proc write::writeNodalConditions { keyword } {
    set doc $gid_groups_conds::doc
    set root [$doc documentElement]
    set xp1 "[spdAux::getRoute $keyword]/condition/group"
    set groups [$root selectNodes $xp1]
    if {$groups eq ""} {
        set xp1 "[spdAux::getRoute $keyword]/group"
        set groups [$root selectNodes $xp1]
    }
    foreach group $groups {
        set cid [[$group parent] @n]
        set groupid [$group @n]
        ::write::writeGroupMesh $cid $groupid "nodal"
    }
}

proc write::GetFormatDict { groupid n num} {
    set f "%10d [format "%10d" $n] [string repeat "%10d " $num]\n"
    #set f [subst $f]
    return [dict create $groupid $f]
}

proc write::getEtype {ov} {
    set isquadratic [isquadratic]
    set ret [list "" ""]
    
    if {$ov eq "point"} {
	set ret [list "Point" 1]
    }
    
    if {$ov eq "line"} {
	switch $isquadratic {
	    0 { set ret [list "Linear" 2] }
	    default { set ret [list "Linear" 2] }                                         
	} 
    }
    
    if {$ov eq "surface"} {
	foreach ielem [lrange [GiD_Info Mesh] 1 end] {
	    switch $ielem {
		Triangle {          
		    switch $isquadratic {
		        0 { set ret [list "Triangle" 3]  }
		        default { set ret [list "Triangle" 6]  }
		    }
		}
		Tetrahedra {          
		    switch $isquadratic {
		        0 { set ret [list "Triangle" 3]  }
		        default { set ret [list "Triangle" 6]  }
		    }
		}           
		Quadrilateral {          
		    switch $isquadratic {
		        0 { set ret [list "Quadrilateral" 4]  }                
		        1 { set ret [list "Quadrilateral" 8]  }                
		        2 { set ret [list "Quadrilateral" 9]  }                
		    }
		}
	    }
	}
    }
    
    if {$ov eq "volume"} {
	foreach ielem [lrange [GiD_Info Mesh] 1 end] {
	    switch $ielem {
		Tetrahedra {          
		    switch $isquadratic {
		        0 { set ret [list "Tetrahedra" 4]  }               
		        1 { set ret [list "Tetrahedra" 10] }                
		        2 { set ret [list "Tetrahedra" 10] }  
		    }
		}           
		Hexahedra {          
		    switch $isquadratic {
		        0 { set ret [list "Hexahedra" 8]  }                
		        1 { set ret [list "Hexahedra" 20]  }                
		        2 { set ret [list "Hexahedra" 27]  }                
		    }
		}
	    }
	}
    }
    
    return $ret
}
proc write::isquadratic {} {  
    set err [catch { GiD_Set Model(QuadraticType) } isquadratic]
    if { $err } {
	set isquadratic [lindex [GiD_Info Project] 5]
    }
    return $isquadratic
}

# GiD_Mesh get element $elem_id face $face_id
proc write::GetNodesFromElementFace {elem_id face_id} {
    set inf [GiD_Mesh get element $elem_id]
    set elem_type [lindex $inf 1]
    set nnodes [lindex $inf 2]
    set nodes [list ]
    switch $elem_type {
	Tetrahedra {
	    set matrix {{1 2 3 5 6 7} {2 4 3 9 10 6} {3 4 1 10 8 7} {4 2 1 9 5 8}}
	}
	Triangle {
	    set matrix {{1 2 4} {2 3 5} {3 1 6}}
	}
    }
    # Decrementamos porque la cara con id 1 corresponde a la posicion 0 de la matriz
    incr face_id -1
    set face_matrix [lindex $matrix $face_id]
    foreach node_index $face_matrix {
	set node [lindex $inf [expr $node_index +2]]
	if {$node ne ""} {lappend nodes $node}
    }
    #W "eid $elem_id fid $face_id nds $nodes"
    return $nodes
}


proc write::getPartsGroupsId {} {
    variable parts
    set doc $gid_groups_conds::doc
    set root [$doc documentElement]
    
    set listOfGroups [list ]
    set xp1 "[spdAux::getRoute $parts]/group"
    set groups [$root selectNodes $xp1]
    
    foreach group $groups {
	set groupName [get_domnode_attribute $group n]
	lappend listOfGroups $groupName
    }
    return $listOfGroups
}
proc write::getPartsMeshId {} {
    variable parts
    set doc $gid_groups_conds::doc
    set root [$doc documentElement]
    
    set listOfGroups [list ]

    foreach group [getPartsGroupsId] {
	lappend listOfGroups [getMeshId Parts $group]
    }
    return $listOfGroups
}

proc write::writePartMeshes { } {
    foreach group [getPartsGroupsId] {
	writeGroupMesh Parts $group "Elements"
    } 
}

proc write::dict2json {dictVal} {
    # XXX: Currently this API isn't symmetrical, as to create proper
    # XXX: JSON text requires type knowledge of the input data
    set json ""
    dict for {key val} $dictVal {
	# key must always be a string, val may be a number, string or
	# bare word (true|false|null)
	if {0 && ![string is double -strict $val] && ![regexp {^(?:true|false|null)$} $val]} {
	    set val "\"$val\""
	}
	if {[isDict $val]} {
	    set val [dict2json $val]
	    set val "\[${val}\]"
	} else {
	    set val \"$val\"
	}
	append json "\"$key\": $val," \n
    }
    if {[string range $json end-1 end] eq ",\n"} {set json [string range $json 0 end-2]}
    return "\{${json}\}"
}
proc write::json2dict {JSONtext} {
    string range [
    string trim [
	string trimleft [
	    string map {\t {} \n {} \r {} , { } : { } \[ \{ \] \}} $JSONtext
	    ] {\uFEFF}
	]
    ] 1 end-1
}
proc write::tcl2json { value } {
    # Guess the type of the value; deep *UNSUPPORTED* magic!
    # display the representation of a Tcl_Obj for debugging purposes. Do not base the behavior of any command on the results of this one; it does not conform to Tcl's value semantics!
    regexp {^value is a (.*?) with a refcount} [::tcl::unsupported::representation $value] -> type
    if {$value eq ""} {return [json::write array {*}[lmap v $value {tcl2json $v}]]}
    switch $type {
	string {
        if {$value eq "false"} {return [expr "false"]}
        if {$value eq "true"} {return [expr "true"]}
	    return [json::write string $value]
	}
	dict {
	    return [json::write object {*}[
	    dict map {k v} $value {tcl2json $v}]]
	}
	list {
	    return [json::write array {*}[lmap v $value {tcl2json $v}]]
	}
	int - double {
	    return [expr {$value}]
	}
	booleanString {
	    return [expr {$value ? "true" : "false"}]
	}
	default {
	    # Some other type; do some guessing...
	    if {$value eq "null"} {
	    # Tcl has *no* null value at all; empty strings are semantically
	    # different and absent variables aren't values. So cheat!
	    return $value
	    } elseif {[string is integer -strict $value]} {
	    return [expr {$value}]
	    } elseif {[string is double -strict $value]} {
	    return [expr {$value}]
	    } elseif {[string is boolean -strict $value]} {
	    return [expr {$value ? "true" : "false"}]
	    }
	    return [json::write string $value]
	}
    }
}

proc write::WriteJSON {processDict} {
    #W [dict2json $processDict]
    WriteString [write::tcl2json $processDict]
}

proc write::GetDefaultOutputDict {} {
    set outputDict [dict create]
    set resultDict [dict create]
    
    set GiDPostDict [dict create]
    dict set GiDPostDict GiDPostMode                [getValue Results GiDPostMode]
    dict set GiDPostDict WriteDeformedMeshFlag      [getValue Results GiDWriteMeshFlag]
    dict set GiDPostDict WriteConditionsFlag        [getValue Results GiDWriteConditionsFlag]
    #dict set GiDPostDict WriteParticlesFlag        [getValue Results GiDWriteParticlesFlag]
    dict set GiDPostDict MultiFileFlag              [getValue Results GiDMultiFileFlag]
    dict set resultDict gidpost_flags $GiDPostDict
    
    dict set resultDict file_label                 [getValue Results FileLabel]
    set outputCT [getValue Results OutputControlType]
    dict set resultDict output_control_type $outputCT
    if {$outputCT eq "time"} {set frequency [getValue Results OutputDeltaTime]} {set frequency [getValue Results OutputDeltaStep]}
    dict set resultDict output_frequency $frequency     
    
    dict set resultDict body_output           [getValue Results BodyOutput]
    dict set resultDict node_output           [getValue Results NodeOutput]
    dict set resultDict skin_output           [getValue Results SkinOutput]
    
    dict set resultDict plane_output [GetCutPlanesList]
    
    #dict set outputDict write_results "PreMeshing"
    #dict set outputDict plot_graphs false
    #dict set outputDict plot_frequency 0
    #dict set outputDict print_lists true
    #dict set outputDict file_list [list ]
    #dict set outputDict output_time 0.01
    #dict set outputDict volume_output true
    #dict set outputDict add_skin true
    
    dict set resultDict nodal_results [GetResultsList NodalResults]
    dict set resultDict gauss_point_results [GetResultsList ElementResults]
    
    dict set outputDict "result_file_configuration" $resultDict
    dict set outputDict "point_data_configuration" [GetEmptyList]
    return $outputDict
}
proc write::GetEmptyList { } {
    # This is a gipsy code
    set a [list ]
    return $a
}
proc write::GetCutPlanesList { } {
    set doc $gid_groups_conds::doc
    set root [$doc documentElement]

    set list_of_planes [list ]
    
    set xp1 "[spdAux::getRoute Results]/container\[@n='CutPlanes'\]/blockdata"
    set planes [$root selectNodes $xp1]
    
    foreach plane $planes {
	set pdict [dict create]
	set points [split [get_domnode_attribute [$plane firstChild] v] ","]
	set normals [split [get_domnode_attribute [$plane lastChild ] v] ","]
	dict set pdict point $points
	dict set pdict normal $normals
	if {![isVectorNull $normals]} {lappend list_of_planes $pdict}
	unset pdict 
    }
    return $list_of_planes
}

proc write::isVectorNull {vector} {
    set null 1
    foreach component $vector {
	if {$component != 0} {
	    set null 0
	    break
	}
    }
    return $null
}

proc write::GetDataType {value} {
    regexp {^value is a (.*?) with a refcount} [::tcl::unsupported::representation $value] -> type
    return $type
}

proc write::getSolutionStrategyParametersDict {} {
    set solStratUN [apps::getCurrentUniqueName SolStrat]
    set schemeUN [apps::getCurrentUniqueName Scheme]
    
    set solstratName [write::getValue $solStratUN]
    set schemeName [write::getValue $schemeUN]
    set sol [::Model::GetSolutionStrategy $solstratName]
    set sch [$sol getScheme $schemeName]
    
    set paramsPath [apps::getCurrentUniqueName StratParams]
    
    foreach {n in} [$sol getInputs] {
	dict set solverSettingsDict $n [write::getValue $paramsPath $n ]
    }
    foreach {n in} [$sch getInputs] {
	dict set solverSettingsDict $n [write::getValue $paramsPath $n ]
    }
    return $solverSettingsDict
}



proc write::getSolversParametersDict {} {
    set solStratUN  [apps::getCurrentUniqueName SolStrat]
    set solstratName [write::getValue $solStratUN]
    set sol [::Model::GetSolutionStrategy $solstratName]
    set solverSettingsDict [dict create]
    foreach se [$sol getSolversEntries] {
	set solverEntryDict [dict create]
	set un [apps::getCurrentUniqueName "$solstratName[$se getName]"] 
	set solverName [write::getValue $un Solver]
	if {$solverName ne "Default"} {
	    dict set solverEntryDict solver_type $solverName
	    #W "Solver Name $solverName"
	    foreach {n in} [[::Model::GetSolver $solverName] getInputs] {
		# JG temporal, para la precarga de combos
		if {[$in getType] ni [list "bool" "integer" "double"]} {
		    set v ""
		    catch {set v [write::getValue $un $n]}
		    if {$v eq ""} {set v [write::getValue $un $n]}
		    dict set solverEntryDict $n $v
		} {
		    dict set solverEntryDict $n [write::getValue $un $n]
		}
	    }
	    dict set solverSettingsDict [$se getName] $solverEntryDict
	}
	unset solverEntryDict
    }
    return $solverSettingsDict
}


proc ::write::getConditionsParametersDict {un {condition_type "Condition"}} {
    set doc $gid_groups_conds::doc
    set root [$doc documentElement]

    set bcCondsDict [list ]
    
    set xp1 "[spdAux::getRoute $un]/condition/group"
    set groups [$root selectNodes $xp1]
    if {$groups eq ""} {
        set xp1 "[spdAux::getRoute $un]/group"
        set groups [$root selectNodes $xp1]
    }
    foreach group $groups {
	set groupName [$group @n]
	#W "GROUP $groupName"
	set cid [[$group parent] @n]
	set groupId [::write::getMeshId $cid $groupName]
	set condId [[$group parent] @n]
	if {$condition_type eq "Condition"} {
	    set condition [::Model::getCondition $condId]
	} {
	    set condition [::Model::getNodalConditionbyId $condId]
	}
	#W "Condition = $condition"
	set processName [$condition getProcessName]
	set process [::Model::GetProcess $processName]
	set processDict [dict create]
	set paramDict [dict create]
	dict set paramDict mesh_id 0
	dict set paramDict model_part_name $groupId
	
	set process_attributes [$process getAttributes]
	set process_parameters [$process getInputs]
	
	dict set process_attributes process_name [dict get $process_attributes n]
	dict unset process_attributes n
	dict unset process_attributes pn
	
	set processDict [dict merge $processDict $process_attributes]
	foreach {inputName in_obj} $process_parameters {
	    set in_type [$in_obj getType]
	    catch {
		set variable_name [$condition getAttribute VariableName]
		#~ W $variable_name
		# "lindex" is a rough solution. Look for a better one.
		dict set paramDict variable_name [lindex $variable_name 0]
		
	    }
	    if {$in_type eq "vector"} {
		#W "input [$in_obj getName] fix: [$in_obj getFixity]"
		if {[$in_obj getFixity] eq "Imposed"} {
		    set is_fixed_x [expr True]
		    set is_fixed_y [expr True]
		    set is_fixed_z [expr True]
		    if {[$group find n FixX] ne ""} {
		        set is_fixed_x [expr [get_domnode_attribute [$group find n FixX] v] ? True : False]
		    }
		    if {[$group find n FixY] ne ""} {
		        set is_fixed_y [expr [get_domnode_attribute [$group find n FixY] v] ? True : False]
		    }
		    if {[$group find n FixZ] ne ""} {
		        set is_fixed_z [expr [get_domnode_attribute [$group find n FixZ] v] ? True : False]
		    }
		    dict set paramDict is_fixed_x $is_fixed_x
		    dict set paramDict is_fixed_y $is_fixed_y
		    dict set paramDict is_fixed_z $is_fixed_z    
		}
		
		set ValX [expr [get_domnode_attribute [$group find n ${inputName}X] v] ]
		set ValY [expr [get_domnode_attribute [$group find n ${inputName}Y] v] ] 
		set ValZ [expr 0.0]
		catch {set ValZ [expr [get_domnode_attribute [$group find n ${inputName}Z] v]]}
		dict set paramDict $inputName [list $ValX $ValY $ValZ]
	    } elseif {$in_type eq "double"} {
		set value [get_domnode_attribute [$group find n $inputName] v] 
		if {[$group find n Fix] ne ""} {
		    set is_fixed [expr [get_domnode_attribute [$group find n Fix] v] ? True : False]
		    dict set paramDict is_fixed $is_fixed
		}
		dict set paramDict $inputName [expr $value]
	    }
	}
	
	dict set processDict Parameters $paramDict
	lappend bcCondsDict $processDict
    }
    return $bcCondsDict
}


proc write::GetResultsList { un } {
    set doc $gid_groups_conds::doc
    set root [$doc documentElement]
    
    set result [list ]
    set xp1 "[spdAux::getRoute $un]/value"
    set resultxml [$root selectNodes $xp1]
    foreach res $resultxml {
        if {[get_domnode_attribute $res v] in [list "Yes" "True" "1"] && [get_domnode_attribute $res state] ne "hidden"} {
            set name [get_domnode_attribute $res n]
            lappend result $name
        }
    }
    return $result
}

proc write::getAllMaterialParametersDict {matname} {
    variable matun
    set doc $gid_groups_conds::doc
    set root [$doc documentElement]
    
    set md [dict create]
    
    set xp3 [spdAux::getRoute $matun]
    append xp3 [format_xpath {/blockdata[@n="material" and @name=%s]/value} $matname]

    set props [$root selectNodes $xp3]
    foreach prop $props {
        dict set md [$prop @n] [get_domnode_attribute $prop v]
    }
    return $md
}

# Auxiliar
proc write::Duration { int_time } {
    # W "entro con $int_time"
    set timeList [list]
    foreach div {86400 3600 60 1} mod {0 24 60 60} name {day hr min sec} {
	set n [expr {$int_time / $div}]
	if {$mod > 0} {set n [expr {$n % $mod}]}
	if {$n > 1} {
	    lappend timeList "$n ${name}s"
	} elseif {$n == 1} {
	    lappend timeList "$n $name"
	}
    }
    return [join $timeList]
}
 
proc write::getValue { name { it "" } } {
    # quitar de aqui!!!! Cuello de botella
    set doc $gid_groups_conds::doc
    set root [$doc documentElement]
    ##
    set xp [spdAux::getRoute $name]
    set node [$root selectNodes $xp]
    if {$it ne ""} {set node [$node find n $it]}
    
    if {[get_domnode_attribute $node v] eq ""} {
	catch {get_domnode_attribute $node values}
    }
    set v ""
    catch {set v [expr [get_domnode_attribute $node v]]}
    if {$v eq "" } {set v [get_domnode_attribute $node v]}
    # W "name: $name it: $it v: $v"
    return $v
 }
 
 proc write::isBoolean {value} {
    set goodList [list "Yes" "1" "yes" "ok" "YES" "Ok" "OK" "True" "TRUE" "true" "No" "0" "no" "NO" "False" "FALSE" "false"]
    if {$value in $goodList} {return 1} {return 0}
 }

proc write::getStringBinaryValue { name { it "" } } {
    set v [getValue $name $it]
    set goodList [list "Yes" "1" "yes" "ok" "YES" "Ok" "OK" "True" "TRUE" "true"]
    if {$v in $goodList} {return "True" } {return "False"}
}
 
proc write::OpenFile { fn } {
    variable dir
    set filename [file join $dir $fn]
    catch {CloseFile}
    customlib::InitWriteFile $filename
}

proc write::CloseFile { } {
    catch {customlib::EndWriteFile}
}

proc write::WriteString {str} {
    #W [format "%s" $str]
    write_calc_data puts [format "%s" $str]
}

proc write::getMatDict {} {
    variable mat_dict
    return $mat_dict
}

proc write::isDict {value} {
    return [expr {[string is list $value] && ([llength $value]&1) == 0}]
}

proc write::getSpacing {number} {
    set r ""
    for {set i 0} {$i<$number} {incr i} {append r " "}
    return $r
}

proc write::CopyFileIntoModel { filepath } {
    variable dir
    
    set activeapp [::apps::getActiveApp]
    set inidir [apps::getMyDir [$activeapp getName]]
    set totalpath [file join $inidir $filepath]
    file copy -force $totalpath $dir
}
proc write::RenameFileInModel { src target } {
    variable dir
    set srcfile [file join $dir $src]
    set tgtfile [file join $dir $target]
    file rename -force $srcfile $tgtfile
}

write::Init

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
    variable MDPA_loop_control
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
    
    variable MDPA_loop_control
    set MDPA_loop_control 0
}

proc write::initWriteData {partes mats} {
    variable parts
    variable matun
    variable meshes
    set parts $partes
    set matun $mats
    
    variable MDPA_loop_control
    set MDPA_loop_control 0
    
    #set meshes [dict create]
    processMaterials
}

proc write::setGroupsTypeName {name} {
    variable groups_type_name
    set groups_type_name $name
}

# Write Events
proc write::writeEvent { filename } {
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
    if {$::kratos_debug} {
        eval $wevent
    } else {
        if {$errcode eq 0 && [catch {eval $wevent} fid] } {
            W "Problem Writing Project Parameters block:\n$fid\nEnd problems"
            set errcode 1
        }
    }
    
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
    variable MDPA_loop_control
    incr MDPA_loop_control
    if {$MDPA_loop_control > 10} {error [= "Infinite loop on MDPA - Check recursive or cyclic calls"]}
    
    set errcode 0
    set activeapp [::apps::getAppById $appid]
    
    #### MDPA Write ####
    set wevent [$activeapp getWriteModelPartEvent]
    set filename "[file tail [GiD_Info project ModelName]].mdpa"
    
    catch {CloseFile}
    OpenFile $filename
    
    if {$::kratos_debug} {
        eval $wevent
    } else {
        if { [catch {eval $wevent} fid] } {
            W "Problem Writing MDPA block:\n$fid\nEnd problems"
            set errcode 1
        }
    }
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

proc write::writeMaterials { {appid ""}} {
    variable mat_dict
    
    set exclusionList [list "MID" "APPID" "ConstitutiveLaw" "Material" "Element"]
    # We print all the material data directly from the saved dictionary
    foreach material [dict keys $mat_dict] {
        set matapp [dict get $mat_dict $material APPID]
        if {$appid eq "" || $matapp in $appid} {
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
}

proc write::writeNodalCoordinatesOnParts { } {
    set formats [dict create]
    foreach group [getPartsGroupsId] {
        dict set formats $group "%5d %14.5f %14.5f %14.5f\n"
    }
    WriteString "Begin Nodes"
    GiD_WriteCalculationFile nodes $formats
    WriteString "End Nodes"
    WriteString "\n"
}
proc write::writeNodalCoordinates { } {
    # Write the nodal coordinates block
    # Nodes block format
    # Begin Nodes
    # // id          X        Y        Z 
    # End Nodes
    
    WriteString "Begin Nodes"
    customlib::WriteCoordinates "%5d %14.10f %14.10f %14.10f\n"
    WriteString "End Nodes"
    WriteString "\n"
}

proc write::processMaterials { } {
    variable parts
    variable matun
    variable mat_dict
    set doc $gid_groups_conds::doc
    set root [$doc documentElement]
    
    set xp1 "[spdAux::getRoute $parts]/group"
    set xp2 ".//value\[@n='Material']"
    
    set material_number [llength [dict keys $mat_dict] ] 

    foreach gNode [$root selectNodes $xp1] {
        set nodeApp [spdAux::GetAppIdFromNode $gNode]
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
            dict set mat_dict $group APPID $nodeApp
            
            set s1 [$gNode selectNodes ".//value"]
            set s2 [$root selectNodes $xp3]
            set us [join [list $s1 $s2]]
            
            foreach valueNode $us {
                set name [$valueNode getAttribute n]
                set state [get_domnode_attribute $valueNode state]
                if {$state ne "hidden"} {
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
    
    set xp1 "[spdAux::getRoute $parts]/group"
    foreach gNode [$root selectNodes $xp1] {
        write::writeGroupElementConnectivities $gNode
    } 
}

# gNode must be a tree group, have a value n = Element
proc write::writeGroupElementConnectivities { gNode } {
    variable mat_dict
    set formats ""
    set group [get_domnode_attribute $gNode n]
    if { [dict exists $mat_dict $group] } {          
        set mid [dict get $mat_dict $group MID]
        if {[$gNode hasAttribute ov]} {set ov [get_domnode_attribute $gNode ov] } {set ov [get_domnode_attribute [$gNode parent] ov] }
        lassign [getEtype $ov $group] etype nnodes
        if {$nnodes ne ""} {
            set formats [GetFormatDict $group $mid $nnodes]
            if {$etype ne "none"} {
                set kelemtype [get_domnode_attribute [$gNode selectNodes ".//value\[@n='Element']"] v]
                set elem [::Model::getElement $kelemtype]
                #W $kelemtype
                set top [$elem getTopologyFeature $etype $nnodes]
                if {$top ne ""} {
                    set kratosElemName [$top getKratosName]
                    WriteString "Begin Elements $kratosElemName// GUI group identifier: $group"
                    # formats -> tcl dict -> key : group id | value : format of printable information
                    # "%10d [format "%10d" $mid] [string repeat "%10d " $num]\n"
                    # "%10d 1 %10d %10d %10d %10d" -> Tetrahedra with property 1 assigned
                    GiD_WriteCalculationFile connectivities $formats
                    WriteString "End Elements"
                    WriteString ""
                } else {
                    error [= "Element $kelemtype not available for $ov entities on group $group"]
                }
            } else {
                error [= "You have not assigned a proper entity to group $group"]
            }
        } else {
            error [= "You have not assigned a proper entity to group $group"]
        } 
    } else {
        error [= "Group $group not in the material database. Call write::initWriteData first"]
    }
}

proc write::writeConditions { baseUN } {
    set dictGroupsIterators [dict create]
    set doc $gid_groups_conds::doc
    set root [$doc documentElement]
    
    set xp1 "[spdAux::getRoute $baseUN]/condition/group"
    set iter 1
    foreach groupNode [$root selectNodes $xp1] {
        set condid [[$groupNode parent] @n]
        set groupid [get_domnode_attribute $groupNode n]
        if {[$groupNode hasAttribute ov]} {set ov [$groupNode getAttribute ov]} {set ov [[$groupNode parent ] getAttribute ov]}
        set cond [::Model::getCondition $condid]
        lassign [write::getEtype $ov $groupid] etype nnodes
        set kname [$cond getTopologyKratosName $etype $nnodes]
        if {$kname ne ""} {
            set obj [list ]
            WriteString "Begin Conditions $kname// GUI group identifier: $groupid"
            if {$etype eq "Point"} {
                set formats [dict create $groupid "%10d \n"]
                
                #set obj [GiD_WriteCalculationFile nodes -return $formats]
                #set obj [list 3 5]
                set obj [GiD_EntitiesGroups get $groupid nodes]
            } else {
                set formats [write::GetFormatDict $groupid 0 $nnodes]
                set elems [GiD_WriteCalculationFile connectivities -return $formats]
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

proc write::transformGroupName {groupid} {
    set new_parts [list ]
    foreach part [GidUtils::Split $groupid "//"]  {
        foreach {bad good} {" " "_"} {
            lappend new_parts [string map [list $bad $good] $part]
        }
    }
    return [join $new_parts /]
}

proc write::writeGroupMesh { cid group {what "Elements"} {iniend ""} {tableid_list ""} } {
    variable meshes
    variable groups_type_name
    
    set gtn $groups_type_name
    if {![dict exists $meshes [list $cid ${group}]]} {
        set mid [expr [llength [dict keys $meshes]] +1]
        if {$gtn ne "Mesh"} {
            set good_name [write::transformGroupName $group]
            set mid "${cid}_${good_name}"
        }
        dict set meshes [list $cid ${group}] $mid
        set gdict [dict create]
        set f "%10i\n"
        set f [subst $f]
        dict set gdict $group $f
        WriteString "Begin $gtn $mid // Group $group // Subtree $cid"
        if {$tableid_list ne ""} {
            WriteString "    Begin SubModelPartTables"
            foreach tableid $tableid_list {
                WriteString "    $tableid"
            }
            WriteString "    End SubModelPartTables"
        }
        WriteString "    Begin ${gtn}Nodes"
        GiD_WriteCalculationFile nodes -sorted $gdict
        WriteString "    End ${gtn}Nodes"
        WriteString "    Begin ${gtn}Elements"
        if {$what eq "Elements"} {
            GiD_WriteCalculationFile elements -sorted $gdict
        }
        WriteString "    End ${gtn}Elements"
        WriteString "    Begin ${gtn}Conditions"
        if {$what eq "Conditions"} {
            #GiD_WriteCalculationFile elements -sorted $gdict
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

proc write::GetFormatDict { groupid mid num} {
    set f "%10d [format "%10d" $mid] [string repeat "%10d " $num]\n"
    return [dict create $groupid $f]
}

proc write::getEtype {ov group} {
    set isquadratic [isquadratic]
    set ret [list "" ""]
    set b 0
    if {$ov eq "point"} {
        if {$b} {error "Multiple element types in $group over $ov"}
        set ret [list "Point" 1]
        set b 1
    }
    
    if {$ov eq "line"} {
        if {$b} {error "Multiple element types in $group over $ov"}
        switch $isquadratic {
            0 { set ret [list "Linear" 2] }
            default { set ret [list "Linear" 2] }                                         
        } 
    }
    
    if {$ov eq "surface"} {
        if {[GiD_EntitiesGroups get $group elements -count -element_type Triangle]} {
            if {$b} {error "Multiple element types in $group over $ov"}
		    switch $isquadratic {
		        0 { set ret [list "Triangle" 3]  }
		        default { set ret [list "Triangle" 6]  }
		    }
            set b 1
		}
        if {[GiD_EntitiesGroups get $group elements -count -element_type Quadrilateral]} {
            if {$b} {error "Multiple element types in $group over $ov"}
		    switch $isquadratic {
		        0 { set ret [list "Quadrilateral" 4]  }                
		        1 { set ret [list "Quadrilateral" 8]  }                
		        2 { set ret [list "Quadrilateral" 9]  }                
		    }
            set b 1
		}
        if {[GiD_EntitiesGroups get $group elements -count -element_type Circle]} {
            if {$b} {error "Multiple element types in $group over $ov"}
		    switch $isquadratic {
		        0 { set ret [list "Circle" 1]  }                 
		    }
            set b 1
		}
    }
    
    if {$ov eq "volume"} {
        if {[GiD_EntitiesGroups get $group elements -count -element_type Tetrahedra]} {
            if {$b} {error "Multiple element types in $group over $ov"}
		    switch $isquadratic {
		        0 { set ret [list "Tetrahedra" 4]  }               
		        1 { set ret [list "Tetrahedra" 10] }                
		        2 { set ret [list "Tetrahedra" 10] }  
		    }
            set b 1
		}
        if {[GiD_EntitiesGroups get $group elements -count -element_type Hexahedra]} {
            if {$b} {error "Multiple element types in $group over $ov"}
		    switch $isquadratic {
		        0 { set ret [list "Hexahedra" 8]  }                
		        1 { set ret [list "Hexahedra" 20]  }                
		        2 { set ret [list "Hexahedra" 27]  }                
		    }
            set b 1
	    }
        if {[GiD_EntitiesGroups get $group elements -count -element_type Prism]} {
            if {$b} {error "Multiple element types in $group over $ov"}
		    switch $isquadratic {
		        0 { set ret [list "Hexahedra" 6]  }                
		        1 { set ret [list "Hexahedra" 1]  }                
		        2 { set ret [list "Hexahedra" 27]  }                
		    }
            set b 1
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
    # Decrementamos porque la cara con id = i corresponde a la posicion i-1 de la matriz
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
        if {$value eq "null"} {return null}
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
        if {[isBooleanFalse $value]} {return [expr "false"]}
        if {[isBooleanTrue $value]} {return [expr "true"]}
	    return [json::write string $value]
	    #return [expr {$value ? "true" : "false"}]
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
    WriteString [write::tcl2json $processDict]
}

proc write::GetDefaultOutputDict { {appid ""} } {
    set outputDict [dict create]
    set resultDict [dict create]
    
    if {$appid eq ""} {set results_UN Results } {set results_UN [apps::getAppUniqueName $appid Results]}
    set GiDPostDict [dict create]
    dict set GiDPostDict GiDPostMode                [getValue $results_UN GiDPostMode]
    dict set GiDPostDict WriteDeformedMeshFlag      [getValue $results_UN GiDWriteMeshFlag]
    dict set GiDPostDict WriteConditionsFlag        [getValue $results_UN GiDWriteConditionsFlag]
    dict set GiDPostDict MultiFileFlag              [getValue $results_UN GiDMultiFileFlag]
    dict set resultDict gidpost_flags $GiDPostDict
    
    dict set resultDict file_label                 [getValue $results_UN FileLabel]
    set outputCT [getValue $results_UN OutputControlType]
    dict set resultDict output_control_type $outputCT
    if {$outputCT eq "time"} {set frequency [getValue $results_UN OutputDeltaTime]} {set frequency [getValue $results_UN OutputDeltaStep]}
    dict set resultDict output_frequency $frequency     
    
    dict set resultDict body_output           [getValue $results_UN BodyOutput]
    dict set resultDict node_output           [getValue $results_UN NodeOutput]
    dict set resultDict skin_output           [getValue $results_UN SkinOutput]
    
    dict set resultDict plane_output [GetCutPlanesList $results_UN]

    dict set resultDict nodal_results [GetResultsList $results_UN OnNodes]
    dict set resultDict gauss_point_results [GetResultsList $results_UN OnElement]
    
    dict set outputDict "result_file_configuration" $resultDict
    dict set outputDict "point_data_configuration" [GetEmptyList]
    return $outputDict
}
proc write::GetEmptyList { } {
    # This is a gipsy code
    set a [list ]
    return $a
}
proc write::GetCutPlanesList { {results_UN Results} } {
    set doc $gid_groups_conds::doc
    set root [$doc documentElement]

    set list_of_planes [list ]
    
    set xp1 "[spdAux::getRoute $results_UN]/container\[@n='CutPlanes'\]/blockdata"
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


proc write::getSubModelPartNames { args } {
    set doc $gid_groups_conds::doc
    set root [$doc documentElement]
    
    set listOfProcessedGroups [list ]
    set groups [list ]
    foreach un $args {
        set xp1 "[spdAux::getRoute $un]/condition/group"
        set xp2 "[spdAux::getRoute $un]/group"
        set grs [$root selectNodes $xp1]
        if {$grs ne ""} {lappend groups {*}$grs}
        set grs [$root selectNodes $xp2]
        if {$grs ne ""} {lappend groups {*}$grs}
    }
    foreach group $groups {
        set groupName [$group @n]
        set cid [[$group parent] @n]
        set gname [::write::getMeshId $cid $groupName]
        if {$gname ni $listOfProcessedGroups} {lappend listOfProcessedGroups $gname}
    }
    
    return $listOfProcessedGroups
}


proc write::getSolversParametersDict { {appid ""} } {
    if {$appid eq ""} {
        set appid [apps::getActiveAppId]
    }
    set solStratUN [apps::getAppUniqueName $appid SolStrat]
    set solstratName [write::getValue $solStratUN]
    set sol [::Model::GetSolutionStrategy $solstratName]
    set solverSettingsDict [dict create]
    foreach se [$sol getSolversEntries] {
        set solverEntryDict [dict create]
        set un [apps::getAppUniqueName $appid "$solstratName[$se getName]"]
        if {[spdAux::getRoute $un] ne ""} {
            set solverName [write::getValue $un Solver]
            if {$solverName ne "Default"} {
                dict set solverEntryDict solver_type $solverName
                set solver [::Model::GetSolver $solverName]
                foreach {n in} [$solver getInputs] {
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
        set cid [[$group parent] @n]
        set groupId [::write::getMeshId $cid $groupName]
        set condId [[$group parent] @n]
        if {$condition_type eq "Condition"} {
            set condition [::Model::getCondition $condId]
        } {
            set condition [::Model::getNodalConditionbyId $condId]
        }
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
        catch {
            set variable_name [$condition getAttribute VariableName]
            # "lindex" is a rough solution. Look for a better one.
            if {$variable_name ne ""} {dict set paramDict variable_name [lindex $variable_name 0]}
        }
        foreach {inputName in_obj} $process_parameters {
            set in_type [$in_obj getType]
            if {$in_type eq "vector"} {
                if {[$in_obj getAttribute vectorType] eq "bool"} {
                    set ValX [expr [get_domnode_attribute [$group find n ${inputName}X] v] ? True : False]
                    set ValY [expr [get_domnode_attribute [$group find n ${inputName}Y] v] ? True : False]
                    set ValZ [expr False]
                    catch {set ValZ [expr [get_domnode_attribute [$group find n ${inputName}Z] v] ? True : False]}
                    dict set paramDict $inputName [list $ValX $ValY $ValZ]
                } {
                    if {[$in_obj getAttribute "enabled"] in [list "1" "0"]} {
                        foreach i [list "X" "Y" "Z"] {
                            if {[expr [get_domnode_attribute [$group find n Enabled_$i] v] ] ne "Yes"} {
                                set Val$i null
                            } else {
                                set printed 0
                                if {[$in_obj getAttribute "function"] eq "1"} {
                                    if {[get_domnode_attribute [$group find n "ByFunction$i"] v]  eq "Yes"} {
                                        set funcinputName "${i}function_$inputName"
                                        set value [get_domnode_attribute [$group find n $funcinputName] v]
                                        set Val$i $value
                                        set printed 1
                                    }
                                } 
                                if {!$printed} {
                                    set value [expr [get_domnode_attribute [$group find n ${inputName}$i] v] ]
                                    set Val$i $value
                                }
                            }
                        } 
                    } else {
                        set ValX [expr [get_domnode_attribute [$group find n ${inputName}X] v] ]
                        set ValY [expr [get_domnode_attribute [$group find n ${inputName}Y] v] ] 
                        set ValZ [expr 0.0]
                        catch {set ValZ [expr [get_domnode_attribute [$group find n ${inputName}Z] v]]}
                    }
                    dict set paramDict $inputName [list $ValX $ValY $ValZ]
                } 
            } elseif {$in_type eq "double" || $in_type eq "integer"} {
                set printed 0
                if {[$in_obj getAttribute "function"] eq "1"} {
                    if {[get_domnode_attribute [$group find n "ByFunction"] v]  eq "Yes"} {
                        set funcinputName "function_$inputName"
                        set value [get_domnode_attribute [$group find n $funcinputName] v]
                        dict set paramDict $inputName $value
                        set printed 1
                    }
                } 
                if {!$printed} {
                    set value [get_domnode_attribute [$group find n $inputName] v] 
                    dict set paramDict $inputName [expr $value]
                }
            } elseif {$in_type eq "bool"} {
                set value [get_domnode_attribute [$group find n $inputName] v] 
                set value [expr $value ? True : False]
                dict set paramDict $inputName [expr $value]
            } elseif {$in_type eq "tablefile"} {
                set value [get_domnode_attribute [$group find n $inputName] v] 
                dict set paramDict $inputName $value
            } else {
                if {[get_domnode_attribute [$group find n $inputName] state] ne "hidden" } {
                    set value [get_domnode_attribute [$group find n $inputName] v] 
                    dict set paramDict $inputName $value
                }
            }
        }
        if {[$group find n Interval] ne ""} {dict set paramDict interval [write::getInterval  [get_domnode_attribute [$group find n Interval] v]] }
        dict set processDict Parameters $paramDict
        lappend bcCondsDict $processDict
    }
    return $bcCondsDict
}

proc write::GetResultsList { un {cnd ""} } {
    set doc $gid_groups_conds::doc
    set root [$doc documentElement]
    
    set result [list ]
    if {$cnd eq ""} {set xp1 "[spdAux::getRoute $un]/value"} {set xp1 "[spdAux::getRoute $un]/container\[@n = '$cnd'\]/value"}
    set resultxml [$root selectNodes $xp1]
    foreach res $resultxml {
        if {[get_domnode_attribute $res v] in [list "Yes" "True" "1"] && [get_domnode_attribute $res state] ne "hidden"} {
            set name [get_domnode_attribute $res n]
            lappend result $name
        }
    }
    return $result
}
    
    #"output_process_list" : [{
    #    "python_module"   : "restart_process",
    #    "kratos_module"   : "KratosMultiphysics.SolidMechanicsApplication",
    #    "help"            : "This process writes restart files",
    #    "process_name"    : "RestartProcess",
    #    "Parameters"      : {
    #        "model_part_name"     : "Main Model",
    #        "save_restart"        : true,
    #        "restart_file_name"   : "cutting_test_rigid_2D_json",
    #        "restart_file_label"  : "step",
    #        "output_control_type" : "time",
    #        "output_frequency"    : 0,
    #        "json_output"         : false
    #    }
    #}],

proc write::GetRestartProcess { {un ""} {name "" } } {
    set doc $gid_groups_conds::doc
    set root [$doc documentElement]
    
    set resultDict [dict create ]
    if {$un eq ""} {set un "Restart"}
    if {$name eq ""} {set name "RestartOptions"}
    
    dict set resultDict "python_module" "restart_process"
    dict set resultDict "kratos_module" "KratosMultiphysics.SolidMechanicsApplication"
    dict set resultDict "help" "This process writes restart files"
    dict set resultDict "process_name" "RestartProcess"
    
    set params [dict create]
    set saveValue [write::getStringBinaryValue $un SaveRestart]

    dict set resultDict "process_name" "RestartProcess"
    dict set params "model_part_name" "Main Domain"
    dict set params "save_restart" $saveValue
    dict set params "restart_file_name" [file tail [GiD_Info Project ModelName]]
    set xp1 "[spdAux::getRoute $un]/container\[@n = '$name'\]/value"
    set file_label [getValue $un RestartFileLabel]
    dict set params "restart_file_label" $file_label
    set output_control [getValue $un RestartControlType]
    dict set params "output_control_type" $output_control
    if {$output_control eq "time"} {dict set params "output_frequency" [getValue $un RestartDeltaTime]} {dict set params "output_frequency" [getValue $un RestartDeltaStep]}
    set jsonoutput [write::getStringBinaryValue $un json_output]
    dict set params "json_output" $jsonoutput
	
    
    dict set resultDict "Parameters" $params
    return $resultDict
}

proc write::GetMeshFromCondition { base_UN condition_id } {
    set doc $gid_groups_conds::doc
    set root [$doc documentElement]
    
    set xp1 "[spdAux::getRoute $base_UN]/condition\[@n='$condition_id'\]/group"
    set groups [$root selectNodes $xp1]
    
    set meshes [list ]
    foreach gNode $groups {
        set group [$gNode @n]
        set meshid [getMeshId $condition_id $group]
        lappend meshes $meshid
    }
    return $meshes
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

proc write::getIntervalsDict { { un "Intervals" } {appid "" } } {
    set doc $gid_groups_conds::doc
    set root [$doc documentElement]
    
    set intervalsDict [dict create]
    set xp3 "[spdAux::getRoute $un]/blockdata\[@n='Interval'\]"
    if {$xp3 ne ""} {
        set intervals [$root selectNodes $xp3]
        foreach intNode $intervals {
            set name [get_domnode_attribute $intNode name]
            set xpini "value\[@n='IniTime'\]"
            set xpend "value\[@n='EndTime'\]"
            set ininode [$intNode selectNodes $xpini]
            set endnode [$intNode selectNodes $xpend]
            set ini ""
            set end ""
            catch {set ini [expr [get_domnode_attribute $ininode v]]}
            catch {set end [expr [get_domnode_attribute $endnode v]]}
            if {$ini eq ""} {set ini [get_domnode_attribute $ininode v]}
            if {$end eq ""} {set end [get_domnode_attribute $endnode v]}
            dict set intervalsDict $name [list $ini $end]
        }
    }
    return $intervalsDict
}
proc write::getInterval { interval {un "Intervals"} {appid "" }  } {
    set ini 0.0
    set end 0.0
    set intervals [write::getIntervalsDict $un]
    foreach int [dict keys $intervals] {
        if {$int eq $interval} {lassign [dict get $intervals $int] ini end}
    }
    return [list $ini $end]
}

proc write::SetParallelismConfiguration {{un "Parallelization"} {n "OpenMPNumberOfThreads"}} {
    set nt 0
    set paralleltype [write::getValue ParallelType]
    if {$paralleltype eq "OpenMP"} {
        catch {set nt [write::getValue $un $n]}
        if {$nt} {write::SetEnvironmentVariable OMP_NUM_THREADS $nt} {return 0}
    } else {
        #write::SetEnvironmentVariable OMP_NUM_THREADS 1
    }
}

proc write::SetEnvironmentVariable {name value} {
    set ::env($name) $value
}

proc write::Duration { int_time } {
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
 
proc write::getValueByNode { node } {
    if {[get_domnode_attribute $node v] eq ""} {
        catch {get_domnode_attribute $node dict}
        catch {get_domnode_attribute $node value}
    }
    set v ""
    catch {set v [expr [get_domnode_attribute $node v]]}
    if {$v eq "" } {set v [get_domnode_attribute $node v]}
    return $v
}
proc write::getValue { name { it "" } } {
    set doc $gid_groups_conds::doc
    set root [$doc documentElement]
    ##
    set xp [spdAux::getRoute $name]
    set node [$root selectNodes $xp]
    if {$it ne ""} {set node [$node find n $it]}
    
    return [getValueByNode $node]
 }
 
proc write::isBoolean {value} {
   set goodList [list "Yes" "1" "yes" "ok" "YES" "Ok" "OK" "True" "TRUE" "true" "No" "0" "no" "NO" "False" "FALSE" "false"]
   if {$value in $goodList} {return 1} {return 0}
}
proc write::isBooleanTrue {value} {
   set goodList [list "Yes" "1" "yes" "ok" "YES" "Ok" "OK" "True" "TRUE" "true"]
   if {$value in $goodList} {return 1} {return 0}
}
proc write::isBooleanFalse {value} {
   set goodList [list "No" "0" "no" "NO" "False" "FALSE" "false"]
   if {$value in $goodList} {return 1} {return 0}
}

proc write::getStringBinaryValue { name { it "" } } {
    set v [getValue $name $it]
    return [write::getStringBinaryFromValue $v]
}
proc write::getStringBinaryFromValue {v} {
    set goodList [list "Yes" "1" "yes" "ok" "YES" "Ok" "OK" "True" "TRUE" "true"]
    if {$v in $goodList} {return "true" } {return "false"}
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
    GiD_WriteCalculationFile puts [format "%s" $str]
}

proc write::getMatDict {} {
    variable mat_dict
    return $mat_dict
}
proc write::setMatDict {matdict} {
    variable mat_dict
    set mat_dict $matdict
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


proc write::WriteAssignedValues {condNode} {
    set assignedVector [list 1 0 1]
    set valuesVector [list 0.0 null 0.0]
    
    for {set i 0} {$i<3} {incr i} {
        set assigned [lindex $assignedVector $i]
        if {!$assigned} {set assignedVector [lreplace $assignedVector $i $i null]}
    }
    set ret [dict create value $valuesVector]
    return $ret
}

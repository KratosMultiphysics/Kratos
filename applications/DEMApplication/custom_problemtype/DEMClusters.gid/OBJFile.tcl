
proc GenerateOBJFile { basename dir problemtypedir } {

    ## Source auxiliar procedures
    ## source [file join $problemtypedir OBJFileAuxProcs.tcl]

    ## Start OBJ file
    #set filename [file join $dir ${basename}.obj]
    set filename [file join $dir generic.obj]
    set FileVar [open $filename w]

    set triangles [GiD_Mesh list -element_type {triangle} element]
    set tetrahedra [GiD_Mesh list -element_type {tetrahedra} element]

    set triangle_nodes [list]
    foreach element_id $triangles { ;
        set nodes [lrange [GiD_Mesh get element $element_id] 3 end] ;
        lappend triangle_nodes {*}$nodes ;
    }

    ## Nodes only if belong to any triangle
    set Nodes [GiD_Info Mesh Nodes]

    #puts $FileVar "Ordered vertex list only if associated with any face triangle"
    for {set i 0} {$i < [llength $Nodes]} {incr i 4} {
        if {[lindex $Nodes $i] in $triangle_nodes} {
        incr position
        puts -nonewline $FileVar "v "
        puts -nonewline $FileVar [format  "%.10f" [lindex $Nodes [expr { $i+1 }]]]
        puts -nonewline $FileVar " "
        puts -nonewline $FileVar [format  "%.10f" [lindex $Nodes [expr { $i+2 }]]]
        puts -nonewline $FileVar " "
        puts -nonewline $FileVar [format  "%.10f" [lindex $Nodes [expr { $i+3 }]]]
        puts $FileVar " "
        dict set position_list_dict [lindex $Nodes $i] $position
        }

    }

    set aux 1
    dict set node_db_x $aux 0.0
    dict unset node_db_x $aux

    dict set node_db_y $aux 0.0
    dict unset node_db_y $aux

    dict set node_db_z $aux 0.0
    dict unset node_db_z $aux

    dict set node_db_counter $aux 0
    dict unset node_db_counter $aux

    foreach element_id $triangles { ;
        set nodes [lrange [GiD_Mesh get element $element_id] 3 end] ;
        set ElementNormal_x [lrange [GiD_Mesh get element $element_id normal] 0 0] ;
        set ElementNormal_y [lrange [GiD_Mesh get element $element_id normal] 1 1] ;
        set ElementNormal_z [lrange [GiD_Mesh get element $element_id normal] 2 2] ;

        foreach node $nodes {
            if {[dict exists $node_db_x $node]} {
                set old_value [dict get $node_db_x $node]
                set new_value [expr { $old_value+$ElementNormal_x}]
                # set new_value [expr { $old_value[0]+$ElementNormal_x }] [expr { $old_value[1]+$ElementNormal_y }]
                dict unset node_db_x $node
                dict set node_db_x $node $new_value

                dict incr node_db_counter $node 1
            } else {
                dict set node_db_x $node $ElementNormal_x
                dict set node_db_counter $node 1
            }
        }

        foreach node $nodes {
            if {[dict exists $node_db_y $node]} {
                set old_value [dict get $node_db_y $node]
                set new_value [expr { $old_value+$ElementNormal_y}]
                dict unset node_db_y $node
                dict set node_db_y $node $new_value
            } else {
                dict set node_db_y $node $ElementNormal_y
            }
        }

        foreach node $nodes {
            if {[dict exists $node_db_z $node]} {
                set old_value [dict get $node_db_z $node]
                set new_value [expr { $old_value+$ElementNormal_z}]
                dict unset node_db_z $node
                dict set node_db_z $node $new_value
            } else {
                dict set node_db_z $node $ElementNormal_z
            }
        }
    }


    # Get a list of the keys and sort them
    set node_db_x [lsort -stride 2 $node_db_x]
    set node_db_y [lsort -stride 2 $node_db_y]
    set node_db_z [lsort -stride 2 $node_db_z]

    foreach item [dict keys $node_db_x] {
            set val [dict get $node_db_x $item]
            lappend normals_x $val ;
    }

    foreach item [dict keys $node_db_y] {
            set val [dict get $node_db_y $item]
            lappend normals_y $val ;
    }

    foreach item [dict keys $node_db_z] {
            set val [dict get $node_db_z $item]
            lappend normals_z $val ;
    }


    # normalization is not required according to documentation.
    # List of vertex normals in (x,y,z) form; normals might not be unit vectors.
    # verify that your geometry have all the normals pointing coherent.
    foreach component1 $normals_x component2 $normals_y component3 $normals_z {
        #set magnitude [expr { sqrt( $component1 * $component1 + $component2 * $component2 + $component3 * $component3) }]
        #set component1_normalized [expr { $component1/$magnitude }]
        #set component2_normalized [expr { $component2/$magnitude }]
        #set component3_normalized [expr { $component3/$magnitude }]
        #lappend normals $component1_normalized $component2_normalized $component3_normalized ;
        lappend normals $component1 $component2 $component3 ;
    }

    puts $FileVar " "
    #puts $FileVar "Ordered vertex normals list only if associated with any face triangle"
    for {set i 0} {$i < [llength $normals]} {incr i 3} {
        puts -nonewline $FileVar "vn "
        puts -nonewline $FileVar [format  "%.10f" [lindex $normals [expr { $i }]]]
        puts -nonewline $FileVar " "
        puts -nonewline $FileVar [format  "%.10f" [lindex $normals [expr { $i+1 }]]]
        puts -nonewline $FileVar " "
        puts $FileVar [format  "%.10f" [lindex $normals [expr { $i+2 }]]]
    }

    #puts $FileVar "Triangle faces defined by its associated vertex"
    puts $FileVar " "
    foreach element_id $triangles { ;
        set nodes [lrange [GiD_Mesh get element $element_id] 3 end] ;
        puts $FileVar "f [dict get $position_list_dict [lindex $nodes 0]]//[dict get $position_list_dict [lindex $nodes 0]] [dict get $position_list_dict [lindex $nodes 1]]//[dict get $position_list_dict [lindex $nodes 1]] [dict get $position_list_dict [lindex $nodes 2]]//[dict get $position_list_dict [lindex $nodes 2]]" ;
    }

    close $FileVar
    set OBJOutput [list $triangles $tetrahedra]

    return $OBJOutput
}


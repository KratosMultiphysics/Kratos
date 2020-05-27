proc GenerateOBJFile { basename dir problemtypedir } {

    ## Source auxiliar procedures
    ## source [file join $problemtypedir OBJFileAuxProcs.tcl]

    ## Start OBJ file
    #set filename [file join $dir ${basename}.obj]
    set filename [file join $dir generic.obj]
    set FileVar [open $filename w]

    set triangles [GiD_Mesh list -element_type {triangle} element]
    set tetrahedra [GiD_Mesh list -element_type {tetrahedra} element]
    
    W "list of triangles and tetrahedra"
    W $triangles
    W $tetrahedra

    set triangle_nodes [list]
    foreach element_id $triangles { ;
        set nodes [lrange [GiD_Mesh get element $element_id] 3 end] ;
        lappend triangle_nodes {*}$nodes ;                              
    }

    
    ## Nodes only if belong to any triangle
    set Nodes [GiD_Info Mesh Nodes]
    #set position_list [list]
    W $Nodes
    W $triangle_nodes
    W "________________________________________"

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
        puts -nonewline $FileVar " "
        #puts  $FileVar  $position
        #lappend position_list [lindex $Nodes $i] $position ;  
        dict set position_list_dict [lindex $Nodes $i] $position
        }
        
    }

    # cube.gid position_list
    # 1 1 2 2 3 3 4 4 6 5 7 6 8 7 9 8 10 9 11 10 13 11 14 12 15 13 28 14 29 15 30 16 31 17 32 18 33 19 38 20 39 21 40 22 42 23 43 24 44 25 45 26 46 27 47 28 55 29 56 30 57 31 58 32 59 33 60 34 62 35 63 36 64 37 65 38 66 39 67 40 71 41 72 42 73 43 74 44 75 45 76 46 77 47 78 48 79 49 81 50 82 51 83 52 85 53 86 54 87 55 88 56


    # REF:
    # Lo que tienes que hacer es recorrer los elementos que son los que apuntan a
    # sus nodos, y sumar a cada nodo la aportación de su normal.
    # Aun así, cualquier función script que haga cosas como esta en malla puede
    # ser lenta en mallas muy grandes, por crear arrays o listas grandes, pero si
    # la haces mal no será lenta sino muy-muy-lenta
    
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
        #W "___________next element of the list"   
        #W $element_id 
        set nodes [lrange [GiD_Mesh get element $element_id] 3 end] ;
        #W $nodes
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
        #W "foreach element x "
        #W $node_db_x

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
        #W "foreach element y "
        #W $node_db_y

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
        #W "foreach element z"
        #W $node_db_z
        #W $node_db_counter

    }


    # Get a list of the keys and sort them
    #W "BEFORE ORDERING KEYS"
    #W $node_db_x
    set node_db_x [lsort -stride 2 $node_db_x]
    set node_db_y [lsort -stride 2 $node_db_y]
    set node_db_z [lsort -stride 2 $node_db_z]
    #W "AFTER ORDERING KEYS"
    #W $node_db_x        

    foreach item [dict keys $node_db_x] {
            set val [dict get $node_db_x $item]
            # set counter [dict get $node_db_counter $item]
            # set avg_val_normal [expr { $val/$counter }]
            # lappend normals_x $avg_val_normal ;
            lappend normals_x $val ;        
    }

    foreach item [dict keys $node_db_y] {
            set val [dict get $node_db_y $item]
            # set counter [dict get $node_db_counter $item]
            # set avg_val_normal [expr { $val/$counter }]
            # lappend normals_y $avg_val_normal ; 
            lappend normals_y $val ;          
    }

    foreach item [dict keys $node_db_z] {
            set val [dict get $node_db_z $item]
            # set counter [dict get $node_db_counter $item]
            # set avg_val_normal [expr { $val/$counter }]
            # lappend normals_z $avg_val_normal ; 
            lappend normals_z $val ;         
    }

    #W $normals_x
    #W $normals_y
    #W $normals_z

    # normalization is not required according to documentation. 
    # List of vertex normals in (x,y,z) form; normals might not be unit vectors.
    # verify that your geometry have all the normals pointing inside the volume.
    foreach component1 $normals_x component2 $normals_y component3 $normals_z {
        set magnitude [expr { sqrt( $component1 * $component1 + $component2 * $component2 + $component3 * $component3) }]
        set component1_normalized [expr { $component1/$magnitude }]
        set component2_normalized [expr { $component2/$magnitude }]
        set component3_normalized [expr { $component3/$magnitude }]
        #lappend normals $component1_normalized $component2_normalized $component3_normalized ;
        lappend normals $component1 $component2 $component3 ;
    }
    #W "vector con todas las normales x y z,..."
    #W $normals

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
    #f v1//vn1 v2//vn2 v3//vn3
    puts $FileVar " "
    foreach element_id $triangles { ;
        set nodes [lrange [GiD_Mesh get element $element_id] 3 end] ;
        #puts $FileVar "f [lindex $nodes 0]//[lindex $nodes 0] [lindex $nodes 1]//[lindex $nodes 1] [lindex $nodes 2]//[lindex $nodes 2]" ;
        puts $FileVar "f [dict get $position_list_dict [lindex $nodes 0]]//[dict get $position_list_dict [lindex $nodes 0]] [dict get $position_list_dict [lindex $nodes 1]]//[dict get $position_list_dict [lindex $nodes 1]] [dict get $position_list_dict [lindex $nodes 2]]//[dict get $position_list_dict [lindex $nodes 2]]" ;

    }


    # required format is
    # v 987.823009 -583.341002 122.360997    
    # vn 0.329248 0.150250 0.932213
    # v 979.430974 -499.442995 88.674001
    # vn 0.689488 0.597651 0.409169

    # f 18//18 10//10 9//9     f v1//vn1 v2//vn2 v3//vn3 ...
    # f 18//18 9//9 17//17 
    # f 19//19 11//11 10//10 

    close $FileVar

    set OBJOutput [list $triangles $tetrahedra]

    return $OBJOutput
}


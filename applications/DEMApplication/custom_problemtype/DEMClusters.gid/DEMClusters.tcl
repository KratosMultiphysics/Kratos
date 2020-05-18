
### Brief methodology to generate the cluster file ###
# 1. extract from mesh: OBJ from surface mesh, MSH from the tetrahedra mesh
# 2. Generate SPH via sphereTree external executable passing OBJ file as argument (current issue accessing gid values)
# 3. Clean the SPH file
# 4. Generate Cluster.clu via mesh_to_cluster_converter.cpp passing SPH + MSH as arguments
#    modified and precompiled executable will be required to do this


# current order and links
# 1AfterMeshGeneration
# -ExtractSurfaceTriangles - aux run on calculate, folder and file save
# -GenerateOBJFile

# 2GenerateSPHFileFromOBJFile
# -call_SphereTree


# 3GenerateClusterFile

##  --------------------------------------------------------------------------------------------------------------------------------------------------
##  Current Blockers ## 

##- Data extraction and writing into separate files follwing a specific format ExtractSurfaceTriangles
##- Error with GiD_AccessValue get condition
##- Calling external precompiled cpp (already modified)


##  --------------------------------------------------------------------------------------------------------------------------------------------------


## GiD events --------------------------------------------------------------------------------------------------------------------------------------------------


proc InitGIDProject { dir } {

    # Create buttons to insert algorithm values
    if { [GidUtils::IsTkDisabled] eq 0} {
        GiDMenu::Create "Sphere Cluster Creation" PRE
        #GiDMenu::InsertOption "Sphere Cluster Creation" [list "SphereTree"] 0 PRE "GidOpenConditions \"SphereTree\"" "" ""
        GiDMenu::InsertOption "Sphere Cluster Creation" [list "Define Options"] 0 PRE "GidOpenProblemData" "" ""
        GiDMenu::InsertOption "Sphere Cluster Creation" [list "GenerateSPH" ] 1 PRE [list GenerateSPHFileFromOBJFile] "" ""
        GiDMenu::InsertOption "Sphere Cluster Creation" [list "GenerateCLU" ] 2 PRE [list GenerateClusterFile] "" ""
        GiDMenu::UpdateMenus
    }

    # Load the application scripts
    set scripts_dir [file join $dir .. .. ]
    # set tcl_filename [file join $scripts_dir scripts initptype.tcl]
    # if {[catch {source $tcl_filename} msg]} {
    #     WarnWinText $msg
    #     return 1
    # }

    # Save ProblemTypePath
    set ::DEMClusters::ProblemTypePath $dir
}


proc AfterReadGIDProject { filename } {

    # Save ProblemPath
    set projectpath $filename
    append projectpath .gid
    set ::DEMClusters::ProblemPath $projectpath

    # Save ProblemName
    if {[regexp -all {\\} $filename] > 0} {
        # Windows
        regsub -all {\\} $filename { } filename
    } else {
        # Unix
        regsub -all {/} $filename { } filename
    }
    set filename [lreplace $filename 0 [expr { [llength $filename]-2 }]]
    set ::DEMClusters::ProblemName $filename
}


proc BeforeMeshGeneration {elementsize} {

}

proc AfterMeshGeneration {fail} {

}


proc BeforeRunCalculation { batfilename basename dir problemtypedir gidexe args } {

    source [file join $problemtypedir OBJFile.tcl]
    set OBJOutput [ExtractSurfaceTriangles_ext $basename $dir $problemtypedir]
}


proc GenerateOBJFile { } {
    
    # TODO: Extract required mesh data and format it into a file:
    # Analyze the format of the OBJ and generate the file in GID
    # The format of the OBJ file is as follows:

    # block 1:surface vertex and vertex normal 
    # block 2: faces (with ordered vertex):

    # v 987.823009 -583.341002 122.360997
    # vn 0.329248 0.150250 0.932213
    # v 979.430974 -499.442995 88.674001
    # vn 0.689488 0.597651 0.409169


    # f 18//18 10//10 9//9 
    # f 18//18 9//9 17//17 
    # f 19//19 11//11 10//10 


    # SUR custom files can also be used with the following format:
    # The format of the SUR file is as follows:
    # -The number of vertices
    # -List of vertices with XYZ for position and XYZ for normal
    # -The number of triangles
    # -List of triangles with i0, i1, i2 (indices of the vertices - zero base)
    #                        n0, n1, n2 (neighbouring triangles, n0 shares edge
    #                                    from i0 to i1 etc.)

}

proc GenerateSPHFileFromOBJFile { } {

    # TODO: linked to GenerateSPHFileFromOBJFile button. On pressing Button, execute GenerateSPHFileFromOBJFile.
    call_SphereTree
}


namespace eval DEMClusters {
    variable ProblemName ""
    variable ProblemPath ""
    variable ProblemTypePath ""
}


proc call_SphereTree { } {

    W "executing call_SphereTree"
    set Algorithm [GiD_AccessValue get gendata Algorithm]
    if {$Algorithm == "MakeTreeMedial"} {
        ::DEMClusters::call_TreeMedial
    } elseif {$Algorithm == "MakeTreeGrid"} {
        ::DEMClusters::call_makeTreeGrid
    } elseif {$Algorithm == "MakeTreeSpawn"} {
        ::DEMClusters::call_makeTreeSpawn
    } elseif {$Algorithm == "MakeTreeOctree"} {
        ::DEMClusters::call_makeTreeOctree        
    } elseif {$Algorithm == "MakeTreeHubbard"} {
        ::DEMClusters::call_makeTreeHubbard
    } else {
        W "Select a valid algorithm"
    }

    # set Young_Modulus [GiD_AccessValue get condition Body_Part Young_Modulus]
    
}

proc DEMClusters::call_TreeMedial { } {

    W "executing call_TreeMedial"

    set Algorithm [GiD_AccessValue get gendata Algorithm]
    set branch [GiD_AccessValue get gendata branch]
    set depth [GiD_AccessValue get gendata depth]
    set testerLevels [GiD_AccessValue get gendata testerLevels]
    set numCover [GiD_AccessValue get gendata numCover]
    set minCover [GiD_AccessValue get gendata minCover]
    set initSpheres [GiD_AccessValue get gendata initSpheres]
    set minSpheres [GiD_AccessValue get gendata minSpheres]
    set erFact [GiD_AccessValue get gendata erFact]
    set numSamples [GiD_AccessValue get gendata numSamples]
    set minSamples [GiD_AccessValue get gendata minSamples]
    set genericOBJFilename [file join $::DEMClusters::ProblemPath generic.obj]
    #set filename_obj {$::DEMClusters::ProblemName}.obj ## custom names

    #TODO: now define the arguments and call the external script sphereTree:
    set argv "-depth $depth -branch $branch -numCover $numCover -minCover $minCover -initSpheres $initSpheres -minSpheres $minSpheres -erFact $erFact -testerLevels $testerLevels -verify -nopause -eval -expand -merge -burst -optimise balance -balExcess 0.001 -maxOptLevel 100 $genericOBJFilename"
    set program [file join $::DEMClusters::ProblemTypePath exec $Algorithm]
    exec $program {*}$argv

    # MakeTreeMedial -depth 3 -branch 8 -numCover 10000 -minCover 5 -initSpheres 1000 -minSpheres 200 -erFact 2 -testerLevels 2 -verify -nopause -eval -expand -merge -burst -optimise balance -balExcess 0.001 -maxOptLevel 100 generic.obj
    #set program [lindex $argv 0]
    # set program [file join $::DEMClusters::ProblemTypePath exec MakeTreeMedial.exe]
    # set arguments [lrange $argv 1 end]
    # exec $program {*}$arguments
   
    # TreeMedial ValidArgs:
    # -depth              Depth of the sphere-tree
    # -branch             Branching factor of sphere-tree
    # -numCover           Number of sample points to cover object with
    # -minCover           Minimum number of sample points per triangle
    # -initSpheres        Initial number of spheres in medial axis approx.
    # -minSpheres         Minimum number of spheres to create for each sub
    #                         region of the medial axis approximation.
    # -erFact             Amount by which to reduce the error when refining
    #                         the medial axis approximation.
    # -testerLevels       Controls the number of points to use to represent a
    #                         sphere when evaluating fit.  Use -1 for CONVEX
    #                         objects, 1 will generate 42 points and 2 will
    #                         generate 168 points.
    # -optimise           Which optimisation algorithm to use, SIMPLEX just
    #                         rearranges the spheres to try improve fit, BALANCE
    #                         tries to throw away spheres that don't improve the
    #                         approximation.
    # -maxOptLevel        Maximum level of the sphere-tree to apply the optimiser.
    #                         0 does first set only - i.e. children of level 0.
    # -balExcess          The amount of extra error the BALANCE algorithm is
    #                         allowed to introduce when throwing away error,
    #                         e.g. 0.05 allows a 5 percent increase in the error.
    # -verify             Verify the model is suitable for use
    # -nopause            Don't pause when processing, i.e. batch mode
    # -eval               Evaluate the fit of the sphere-tree and append the info
    #                         to the end of the output file.
    # -merge              Try the MERGE, BURST and EXPAND algorithms.  You can
    # -burst              specify any number of these that you wish.
    # -expand

}

proc DEMClusters::call_makeTreeGrid { } {

    set Algorithm [GiD_AccessValue get gendata Algorithm]
    set branch [GiD_AccessValue get gendata branch]
    set depth [GiD_AccessValue get gendata depth]
    set numCover [GiD_AccessValue get gendata numCover]
    set minCover [GiD_AccessValue get gendata minCover]
    set testerLevels [GiD_AccessValue get gendata testerLevels]
    set genericOBJFilename [file join $::DEMClusters::ProblemPath generic.obj]

    #TODO: now define the arguments and call the external script sphereTree:
    set argv "-depth $depth -branch $branch -numCover $numCover -minCover $minCover -testerLevels $testerLevels -verify -nopause -eval $genericOBJFilename"
    set program [file join $::DEMClusters::ProblemTypePath exec $Algorithm]
    exec $program {*}$argv

    # makeTreeGrid ValidArgs:
    # -depth              Depth of the sphere-tree
    # -branch             Branching factor of sphere-tree
    # -numCover           Number of sample points to cover object with
    # -minCover           Minimum number of sample points per triangle
    # -testerLevels       Controls the number of points to use to represent a
    #                         sphere when evaluating fit.  Use -1 for CONVEX
    #                         objects, 1 will generate 42 points and 2 will
    #                         generate 168 points.
    # -verify             Verify the model is suitable for use
    # -nopause            Don't pause when processing, i.e. batch mode
    # -eval               Evaluate the fit of the sphere-tree and append the info
    #                         to the end of the output file.
}

proc DEMClusters::call_makeTreeSpawn { } {

    set Algorithm [GiD_AccessValue get gendata Algorithm]
    set branch [GiD_AccessValue get gendata branch]
    set depth [GiD_AccessValue get gendata depth]
    set numCover [GiD_AccessValue get gendata numCover]
    set minCover [GiD_AccessValue get gendata minCover]
    set testerLevels [GiD_AccessValue get gendata testerLevels]
    set genericOBJFilename [file join $::DEMClusters::ProblemPath generic.obj]

    #TODO: now define the arguments and call the external script sphereTree:
    set argv "-depth $depth -branch $branch -numCover $numCover -minCover $minCover -testerLevels $testerLevels -verify -nopause -eval $genericOBJFilename"
    set program [file join $::DEMClusters::ProblemTypePath exec $Algorithm]
    exec $program {*}$argv

    # makeTreeSpawn ValidArgs:
    # -depth              Depth of the sphere-tree
    # -branch             Branching factor of sphere-tree
    # -numCover           Number of sample points to cover object with
    # -minCover           Minimum number of sample points per triangle
    # -testerLevels       Controls the number of points to use to represent a
    #                         sphere when evaluating fit.  Use -1 for CONVEX
    #                         objects, 1 will generate 42 points and 2 will
    #                         generate 168 points.
    # -verify             Verify the model is suitable for use
    # -nopause            Don't pause when processing, i.e. batch mode
    # -eval               Evaluate the fit of the sphere-tree and append the info
    #                         to the end of the output file.
}

proc DEMClusters::call_makeTreeOctree { } {

    set Algorithm [GiD_AccessValue get gendata Algorithm]
    set depth [GiD_AccessValue get gendata depth]
    set genericOBJFilename [file join $::DEMClusters::ProblemPath generic.obj]

    #TODO: now define the arguments and call the external script sphereTree:
    set argv "-depth $depth -nopause $genericOBJFilename"
    set program [file join $::DEMClusters::ProblemTypePath exec $Algorithm]
    exec $program {*}$argv

    # makeTreeOctree ValidArgs:
    # -depth              Depth of the sphere-tree
    # -nopause            Don't pause when processing, i.e. batch mode
}

proc DEMClusters::call_makeTreeHubbard { } {

    set Algorithm [GiD_AccessValue get gendata Algorithm]
    set branch [GiD_AccessValue get gendata branch]
    set depth [GiD_AccessValue get gendata depth]
    set numSamples [GiD_AccessValue get gendata numSamples]
    set minSamples [GiD_AccessValue get gendata minSamples]
    set genericOBJFilename [file join $::DEMClusters::ProblemPath generic.obj]

    #TODO: now define the arguments and call the external script sphereTree:
    set argv "-depth $depth -branch $branch -numSamples $numSamples -minSamples $minSamples -nopause $genericOBJFilename"
    set program [file join $::DEMClusters::ProblemTypePath exec $Algorithm]
    exec $program {*}$argv

    # makeTreeHubbard ValidArgs:
    # -depth              Depth of the sphere-tree
    # -branch             Branching factor of sphere-tree
    # -numSamples         Number of sample points to cover object with
    # -minSamples         Minimum number of sample points per triangle
    # makeTreeHubbard  -branch 8 -depth 3 -numSamples 500 -minSamples 1
    #                  -nopause  bunny-1500.obj
}


proc CleanSPHFile { sphfilename } {
    # TODO: Moved inside convert_to_cluster cpp. will Access generic_sph.sph file and execute the partial removal of some lines and columns as specified in the reference.
}

proc GenerateClusterFile { } {
    # TODO: linked to GenerateClusterFile button. Generate the cluster file from the joint information of generic_sph.sph and generic_msh.msh

    # Look for a way to edit, compile and execute mesh_to_clu_converter.cpp with the specified filenames
    # Create cpp executable able to generate cluster file directly from inputs.
    # Or pass always the same generic sph and msh names. Generate generic cluster filename.

    # Execute external compiled exe
    set program mesh_to_clu_converter
    exec $program

}






# proc Dam::write::Triangle2D3Connectivities { ElemId } {

#     set ElementInfo [GiD_Mesh get element $ElemId]
#     #ElementInfo: <layer> <elemtype> <NumNodes> <N1> <N2> ...
#     return "[lindex $ElementInfo 3] [lindex $ElementInfo 4] [lindex $ElementInfo 5]"
# }


#   if {$FEMtoDEM == "AttheCentroid"} {
#                     set nodes_to_delete [list]
#                     set element_ids [GiD_EntitiesGroups get $groupid elements] ;               # get ids of all elements in cgroupid
#                     #array set is_external_element [DEM::write::Compute_External_Elements 3 $groupid $element_ids]

#                     foreach element_id $element_ids { ;                                         # loop on each of the elements by id
#                         set element_nodes [lrange [GiD_Mesh get element $element_id] 3 end] ;   # get the nodes of the element
#                         lappend nodes_to_delete {*}$element_nodes ;                             # add those nodes to the nodes_to_delete list


# GiD_EntitiesLayers assign|get|entity_layer
# To handle the entities that belong to layers


# GiD_EntitiesLayers assign|assign_back_layer|assign_front_layer|get <layer> ?-also_lower_entities? <over> <selection>


# To add or know entities of a layer



# GiD_EntitiesLayers assign <layer> ?-also_lower_entities? ?-also_higher_entities? <over> <selection>



# To assing the selection of entities of kind over to the layer




# <layer> is the full name of the layer




# <-also_lower_entities> is an optional flag, to select also all lower entities of the selected ones (e.g. curves and points of the selected surfaces)




# <-also_higher_entities> is an optional flag, to select also all higher entities of the selected ones (e.g. volumes of the selected surfaces)




# <over> could be points, lines, surfaces, volumes, nodes, elements, all_geometry, all_mesh




# <selection> is a list of integer entity id's starting from 1.



# In case of all_geometry is expected a list with 4 items with the list of ids of points, lines, surfaces and volumes.



# In case of all_mesh is expected a list with 2 itemos with the list of ids of nodes and elements respectivelly















# GiD_Mesh get element <num|from_face|from_edge> ?face|face_linear|num_faces|edge_linear|num_edges|normal|tangent|center|connectivities ?<face_id>|<edge_id>??


# <num> is the identifier (integer > 0) for the element to be asked



# face optional, instead of the element nodes it returns the nodes of the face, first the linear corner nodes and then the quadratic nodes



# face_linear optional, instead of the element nodes it returns only the linear corner nodes, also is the element is quadratic



# num_faces returns the amount of faces of the element (for surface elements its edges act as faces)



# <face_id> is the local face index from 1 to the number of faces of the element. If <face_id> is missing then a list with all faces is returned



# edge_linear optional, instead of the element nodes it returns only the linear edge nodes, also is the element is quadratic



# <edge_id> is the local edge index from 1 to the number of edges of the element. If <edge_id> is missing then a list with all edges is returned



# num_edges returns the amount of edges of the element


# get element return the list: <element_layer> <elemtype> <nnode> <N1> ... <Nnnode>

# get element face|face_linear: <N1_face> ... <Nnnode_face>

# get element edge_linear: <N1_edge> <N2_edge>



# If from_face is specified then the command has this syntax


# GiD_Mesh get element from_face <face_nodes> ?-ordered?

# it find and return the list of element ids that have a face with these nodes


# <face_nodes> is the list of integer ids of the face nodes {<face_node_1> ... <face_node_n>} (only corner lineal nodes must be specified in the list)


# if -ordered is specified then only faces with the same orientation of the nodes will be taken into account (else the ordenation of the face nodes doesn't matter)


# If from_edge is specified then the command has this syntax


# GiD_Mesh get element from_edge <edge_nodes>

# it find and return the list of element ids that have an edge with these nodes


# <edge_nodes> is the list of integer ids of the edge nodes {<edge_node_1> <edge_node_2>} (only corner lineal nodes must be specified in the list)
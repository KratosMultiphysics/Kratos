
### Brief methodology to generate the cluster file ###
# 1. extract from mesh: OBJ from surface mesh, MSH from the tetrahedra mesh
# 2. Generate SPH via sphereTree external executable passing OBJ file as argument (current issue accessing gid values)
# 3. Clean the SPH file
# 4. Generate Cluster.clu via mesh_to_cluster_converter.cpp passing SPH + MSH as arguments
#    modified and precompiled executable will be required to do this


# current order and links
# 1AfterMeshGeneration
# -ExtractSurfaceTriangles
# -GenerateOBJFile

# 2GenerateSPHFileFromOBJFile
# -call_SphereTree
# -CleanSPHFile

# 3GenerateClusterFile

##  --------------------------------------------------------------------------------------------------------------------------------------------------


## GiD events --------------------------------------------------------------------------------------------------------------------------------------------------


proc InitGIDProject { dir } {

    # Create buttons to insert algorithm values
    if { [GidUtils::IsTkDisabled] eq 0} {
        GiDMenu::Create "Sphere Cluster Creation" PRE
        #GiDMenu::InsertOption "Sphere Cluster Creation" [list "SphereTree"] 0 PRE "GidOpenConditions \"SphereTree\"" "" ""
        GiDMenu::InsertOption "Sphere Cluster Creation" [list "Define SphereTree"] 0 PRE "GidOpenProblemData" "" ""
        GiDMenu::InsertOption "Sphere Cluster Creation" [list "GenerateSPH" ] 1 PRE [list GenerateSPHFileFromOBJFile] "" ""
        GiDMenu::Create "GenerateCLU" PRE
        GiDMenu::InsertOption "GenerateCLU" [list "GenerateCLU" ] 0 PRE [list GenerateClusterFile] "" ""
        GiDMenu::UpdateMenus
    }

    # Load the application scripts
    set scripts_dir [file join $dir .. .. ]
    # set tcl_filename [file join $scripts_dir scripts initptype.tcl]
    # if {[catch {source $tcl_filename} msg]} {
    #     WarnWinText $msg
    #     return 1
    # }
}




proc BeforeMeshGeneration {elementsize} {
    W "execute BeforeMeshGeneration"
}

proc AfterMeshGeneration {fail} {
    W "execute AfterMeshGeneration"
    ExtractSurfaceTriangles 
}

proc ExtractSurfaceTriangles { } {
    set all_mesh [GiD_EntitiesLayers get Layer0 all_mesh]
    W $all_mesh
    # If <over> is all_mesh then is obtained a list with 2 sublists: node id's and element id's
    # List of triangles defined by: its vertex and vertex normals
    #                               faces: ordered vertex ids

    # TODO: Call GenerateOBJFile with the information extracted from the surface mesh
    # GenerateOBJFile
}

proc BeforeRunCalculation { batfilename basename dir problemtypedir gidexe args } {

    # TODO:  On running calculate it will just to to locate and run a DEMClusters.bat

    # Write file
    #source [file join $problemtypedir file.tcl]

    # TODO: move out and create a button to call this. Do not depend on proc BeforeRunCalculation
    # WriteSphereTreeParameters.json
    # set TableDict [lindex $MDPAOutput 1]
    # source [file join $problemtypedir SphereTreeParameters.tcl]
    # WriteSphereTreeParameters $basename $dir $problemtypedir $TableDict

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

    W "executing GenerateSPHFileFromOBJFile"
    # TODO: linked to GenerateSPHFileFromOBJFile button. On pressing Button, execute GenerateSPHFileFromOBJFile.
    call_SphereTree
}

proc call_SphereTree { } {

    set $Young_Modulus [GiD_AccessValue get condition Body_Part Young_Modulus]
    W $Young_Modulus

    # TODO: error: seems the variables in the prb cannot be located.
    set $Algorithm [GiD_AccessValue get gendata Algorithm]
    W $Algorithm
    if {$Algorithm == "makeTreeMedial"} {
        call_TreeMedial
    } elseif {$Algorithm == "makeTreeGrid"} {
        call_makeTreeGrid
    } elseif {$Algorithm == "makeTreeSpawn"} {
        call_makeTreeSpawn
    } elseif {$Algorithm == "makeTreeOctree"} {
        call_makeTreeOctree        
    } elseif {$Algorithm == "makeTreeHubbard"} {
        call_makeTreeHubbard
    } else {
        W "Select a valid algorithm"
    }
    
    CleanSPHFile
}

proc call_TreeMedial { } {

    set $Algorithm [GiD_AccessValue get gendata Algorithm]
    set $branch [GiD_AccessValue get gendata branch]
    set $depth [GiD_AccessValue get gendata depth]
    set $testerLevels [GiD_AccessValue get gendata testerLevels]
    set $numCover [GiD_AccessValue get gendata numCover]
    set $minCover [GiD_AccessValue get gendata minCover]
    set $initSpheres [GiD_AccessValue get gendata initSpheres]
    set $minSpheres [GiD_AccessValue get gendata minSpheres]
    set $erFact [GiD_AccessValue get gendata erFact]
    set $numSamples [GiD_AccessValue get gendata numSamples]
    set $minSamples [GiD_AccessValue get gendata minSamples]
    set $genericOBJFilename  generic_obj.obj

    #TODO: now define the arguments and call the external script sphereTree:
    set argv {$Algorithm -depth $depth -branch $branch -numCover $numCover -minCover $minCover -initSpheres $initSpheres -minSpheres $minSpheres -erFact $erFact -testerLevels $testerLevels -verify -nopause -eval -expand -merge -burst -optimise balance -balExcess 0.001 -maxOptLevel 100 $genericOBJFilename}

    set program [lindex $argv 0]
    set arguments [lrange $argv 1 end]
    # Execute external script
    exec $program {*}$arguments


    #  Other References ##################
    # If $argv is foo bar baz, then
    # spawn [lindex $argv 0] [lrange $argv 1 end]
    # will invoke foo with 1 argument: "bar baz"
    # spawn [lindex $argv 0] {*}[lrange $argv 1 end]
    # will invoke foo with 2 arguments: "bar" and "baz"

    # set output [exec makeTreeMedial file_name.obj]
    # set output_full [exec makeTreeMedial -branch NS -depth 1 -testerLevels 2 -numCover 10000 -minCover 5 -initSpheres 1000 -minSpheres 200 -erFact 2 -verify -nopause -eval -expand -merge -burst -optimise balance -balExcess 0.001 -maxOptLevel 100 file_name.obj]

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

proc call_makeTreeGrid { } {

    set $Algorithm [GiD_AccessValue get gendata Algorithm]
    set $branch [GiD_AccessValue get gendata branch]
    set $depth [GiD_AccessValue get gendata depth]
    set $numCover [GiD_AccessValue get gendata numCover]
    set $minCover [GiD_AccessValue get gendata minCover]
    set $testerLevels [GiD_AccessValue get gendata testerLevels]
    set $genericOBJFilename  generic_obj.obj

    #TODO: now define the arguments and call the external script sphereTree:
    set argv {$Algorithm -depth $depth -branch $branch -numCover $numCover -minCover $minCover -testerLevels $testerLevels -verify -nopause -eval $genericOBJFilename}

    set program [lindex $argv 0]
    set arguments [lrange $argv 1 end]
    # Execute external script
    exec $program {*}$arguments

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

proc call_makeTreeSpawn { } {

    set $Algorithm [GiD_AccessValue get gendata Algorithm]
    set $branch [GiD_AccessValue get gendata branch]
    set $depth [GiD_AccessValue get gendata depth]
    set $numCover [GiD_AccessValue get gendata numCover]
    set $minCover [GiD_AccessValue get gendata minCover]
    set $testerLevels [GiD_AccessValue get gendata testerLevels]
    set $genericOBJFilename  generic_obj.obj

    #TODO: now define the arguments and call the external script sphereTree:
    set argv {$Algorithm -depth $depth -branch $branch -numCover $numCover -minCover $minCover -testerLevels $testerLevels -verify -nopause -eval $genericOBJFilename}

    set program [lindex $argv 0]
    set arguments [lrange $argv 1 end]
    # Execute external script
    exec $program {*}$arguments

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

proc call_makeTreeOctree { } {

    set $Algorithm [GiD_AccessValue get gendata Algorithm]
    set $depth [GiD_AccessValue get gendata depth]
    set $genericOBJFilename  generic_obj.obj

    #TODO: now define the arguments and call the external script sphereTree:
    set argv {$Algorithm -depth $depth -nopause $genericOBJFilename}

    set program [lindex $argv 0]
    set arguments [lrange $argv 1 end]
    # Execute external script
    exec $program {*}$arguments

    # makeTreeOctree ValidArgs:
    # -depth              Depth of the sphere-tree
    # -nopause            Don't pause when processing, i.e. batch mode
}

proc call_makeTreeOctree { } {

    set $Algorithm [GiD_AccessValue get gendata Algorithm]
    set $branch [GiD_AccessValue get gendata branch]
    set $depth [GiD_AccessValue get gendata depth]
    set $numSamples [GiD_AccessValue get gendata numSamples]
    set $minSamples [GiD_AccessValue get gendata minSamples]
    set $genericOBJFilename  generic_obj.obj

    #TODO: now define the arguments and call the external script sphereTree:
    set argv {$Algorithm -depth $depth -branch $branch -numSamples $numSamples -minSamples $minSamples $genericOBJFilename}

    set program [lindex $argv 0]
    set arguments [lrange $argv 1 end]
    # Execute external script
    exec $program {*}$arguments

    # makeTreeOctree ValidArgs:
    # -depth              Depth of the sphere-tree
    # -branch             Branching factor of sphere-tree
    # -numSamples         Number of sample points to cover object with
    # -minSamples         Minimum number of sample points per triangle
    # makeTreeHubbard  -branch 8 -depth 3 -numSamples 500 -minSamples 1
    #                  -nopause  bunny-1500.obj
}


proc CleanSPHFile { sphfilename } {
    # TODO: Access generic_sph.sph file and execute the partial removal of some lines and columns as specified in the reference.

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


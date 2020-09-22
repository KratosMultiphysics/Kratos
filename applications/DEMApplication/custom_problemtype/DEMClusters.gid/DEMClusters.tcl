
### Brief methodology to generate the cluster file ###
# 1. extract from mesh: OBJ from surface mesh, MSH from the tetrahedra mesh
# 2. Generate SPH via sphereTree external executable passing OBJ file as argument (current issue accessing gid values)
# 3. Clean the SPH file
# 4. Generate Cluster.clu via mesh_to_cluster_converter.cpp passing SPH + MSH as arguments
#    modified and precompiled executable will be required to do this


## Extra features:
## - Lateral icons
## - Show progress onscreen while generating SPH via spheretree
## - Remove finished calculate splash.( or modify with OBJ created succesfully)


##  --------------------------------------------------------------------------------------------------------------------------------------------------
##  Current Issues ##
## - Created option to capture and cancel external exec running from GiD, but GiD freezes during execution so cancel cannot be pressed.

##  --------------------------------------------------------------------------------------------------------------------------------------------------
##  Fixed Issues ##

## - Are ALL calculated Vertex normals correct? verified manually with testcubev4
## - Sphere tree paths with spaces
## - Dependencies in prb.  each line for an algorithm, one parenthesis
## - verify that your geometry have all the normals coherent. (should be as align normals has been added before meshing)
## - add some info on the help based on the manual
## - Calling external precompiled cpp, add criteria for negative radius when deleting line in sph
## - plan on msh to clu functionaly design of precompiled executable. args, paths and location,
## -     define arguments path and call exe from inside exec folder.
## - add dummy .bat to avoid error showing both unix.bat and win.bat
## - cluster visualizaton in GID pre
## - add export gidmesh as generic.msh on calculate
## - RECOMMENEDED MEDIAL - spheretree throws error if algorithm parameters are not quite good. example: small geom with high numsamples
## - Cluster visualization for both sph and principal axis of inertia.
## - Cancel button to stop external execs
##  --------------------------------------------------------------------------------------------------------------------------------------------------

## GiD events --------------------------------------------------------------------------------------------------------------------------------------------------


proc InitGIDProject { dir } {

    # Create buttons to insert algorithm values
    if { [GidUtils::IsTkDisabled] eq 0} {
        GiDMenu::Create "Sphere Cluster Creation" PRE
        #GiDMenu::InsertOption "Sphere Cluster Creation" [list "SphereTree"] 0 PRE "GidOpenConditions \"SphereTree\"" "" ""
        GiDMenu::InsertOption "Sphere Cluster Creation" [list "Select SPH options"] 0 PRE "GidOpenProblemData" "" ""
        GiDMenu::InsertOption "Sphere Cluster Creation" "---" 1 PRE "" "" "" replace =
        GiDMenu::InsertOption "Sphere Cluster Creation" [list "0.Generate OBJ, SPH, CLU and visualize" ] 2 PRE [list OneClickGo] "" ""
        GiDMenu::InsertOption "Sphere Cluster Creation" [list "1.Generate OBJ and export MSH" ] 3 PRE [list GenerateOBJMSH] "" ""
        GiDMenu::InsertOption "Sphere Cluster Creation" [list "2.Generate SPH file" ] 4 PRE [list GenerateSPHFileFromOBJFile] "" ""
        GiDMenu::InsertOption "Sphere Cluster Creation" [list "Cancel SPH generation" ] 5 PRE [list CancelSphereTree] "" ""
        GiDMenu::InsertOption "Sphere Cluster Creation" [list "3.Generate CLU file" ] 6 PRE [list GenerateClusterFile] "" ""
        #GiDMenu::InsertOption "Sphere Cluster Creation" [list "Visualize cluster over geometry" ] 4 PRE [list ReadClusterFileintoGeometry] "" ""
        GiDMenu::InsertOption "Sphere Cluster Creation" [list "Visualize cluster" ] 7 PRE [list ReadSPHFileintoMesh] "" ""
        GiDMenu::InsertOption "Sphere Cluster Creation" [list "Visualize cluster in principal axis" ] 8 PRE [list ReadClusterFileintoMesh] "" ""
        #GiDMenu::InsertOption "Sphere Cluster Creation" [list "Delete cluster over geometry" ] 6 PRE [list DeleteSpheresGeometry] "" ""
        GiDMenu::InsertOption "Sphere Cluster Creation" [list "Remove cluster visualization" ] 9 PRE [list DeleteSpheresMesh] "" ""

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

    Splash
}


proc Splash {} {
    set prev_splash_state [GiD_Set SplashWindow]
    GiD_Set SplashWindow 1 ;#set temporary to 1 to force show splash without take care of the GiD splash preference
    set off_x 150
    set fnt "Sans 10"
    if { $::tcl_platform(platform) == "windows" } {
	set fnt "verdana 10"
	set off_x 130
    }

    ::GidUtils::Splash [file join $::DEMClusters::ProblemTypePath images splash.png] .splash 1

    set new_splash [ ::GidUtils::VersionCmp 11.1.6d]
    if {[ winfo exists .splash.lv] && ( $new_splash < 0)} {
	.splash.lv configure -font $fnt -background white -foreground black \
	    -relief solid -borderwidth 1 -padx 12 -pady 3
	update
    }
    GiD_Set SplashWindow $prev_splash_state
}





proc DeleteSpheresGeometry { } {
    set volume_id_list [GiD_Geometry list volume 2:end]
    GiD_Layers create spheres_to_delete
    foreach id $volume_id_list { ;
        GiD_EntitiesLayers assign spheres_to_delete -also_lower_entities volume $id
    }
    GiD_Layers delete spheres_to_delete
    GiD_Process 'Render Normal

}

proc DeleteSpheresMesh { } {

    set spheres [GiD_Mesh list -element_type {sphere} element]
    # GiD_Mesh delete element $spheres
    foreach id $spheres { ;
        GiD_Process Mescape Meshing EditMesh DeleteElems LowerEntities $id
    }
    GiD_Process Mescape
    GiD_Process 'Render Normal

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

proc GiD_Event_BeforeMeshGeneration {elementsize} {

    # Align the normal
    AlignSurfNormals Outwards
    GiD_MeshData mesh_criteria to_be_meshed 2 surfaces [GiD_Geometry list surface]
}

proc AfterMeshGeneration {fail} {

}

proc OneClickGo {} {

    # Generate OBJ and MSH files on Calculate
    GiD_Process Mescape Utilities Calculate MEscape

    # Generate SPH file from OBJ
    GenerateSPHFileFromOBJFile

    # Generate CLU file from MSH and SPH
    GenerateClusterFile

    # Visualize cluster
    ReadSPHFileintoMesh
}




proc GenerateOBJMSH {} {

    # Generate OBJ and MSH files on Calculate
    GiD_Process Mescape Utilities Calculate MEscape

}

proc BeforeRunCalculation { batfilename basename dir problemtypedir gidexe args } {

    source [file join $problemtypedir OBJFile.tcl]

    # # Already centered when meshing. be sure to center the geometry before generating the OBJ or MSH
    # CenterGeometry

    set OBJOutput [GenerateOBJFile $basename $dir $problemtypedir]

    set modelname [GiD_Info Project ModelName]
    set export_msh [file join ${modelname}.gid generic.msh]

    set projectpath $modelname
    append projectpath .gid
    set ::DEMClusters::ProblemPath $projectpath

    GiD_Process Mescape Files WriteMesh $export_msh
}

proc AfterRunCalculation {basename dir problemtypedir where error errorfilename} {
    #return 0

    return -cancel- ; # To avoid the window
}


proc GenerateSPHFileFromOBJFile { } {
    call_SphereTree
}


namespace eval DEMClusters {
    variable ProblemName ""
    variable ProblemPath ""
    variable ProblemTypePath ""
}

proc CancelSphereTree { } {

    package require gid_cross_platform
    gid_cross_platform::end_process [pid $::DEMClusters::pid]
}


proc call_SphereTree { } {

    set Algorithm [GiD_AccessValue get gendata Algorithm]
    if {$Algorithm == "MakeTreeMedial"} {
        ::DEMClusters::call_makeTreeMedial
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
}

    proc my_finish_proc { pid status args } { ... }

    # returned user_stop must be 0 or 1 (1 to kill the process)

    proc my_progress_proc { args }  { ... return $user_stop }

proc DEMClusters::call_makeTreeMedial { } {

    set genericSPHFilename [file join $::DEMClusters::ProblemPath generic-medial.sph]
    if { [file exists $genericSPHFilename] } {
        W "Previous SPH has been removed. Generating new SPH file."
        file delete $genericSPHFilename
    }

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
    #set genericOBJFilename [file join $::DEMClusters::ProblemPath generic.obj]
    #set genericOBJFilename "\"$genericOBJFilename\""
    set modelname [GiD_Info Project ModelName]
    set genericOBJFilename [file join ${modelname}.gid generic.obj]


    #set filename_obj $::DEMClusters::ProblemName ## custom names
    #append filename_obj .obj
    # set Young_Modulus [GiD_AccessValue get condition Body_Part Young_Modulus]
    set argv "-depth $depth -branch $branch -numCover $numCover -minCover $minCover -initSpheres $initSpheres -minSpheres $minSpheres -erFact $erFact -testerLevels $testerLevels -nopause -eval -expand -merge -burst -optimise balance -balExcess 0.001 -maxOptLevel 100 $genericOBJFilename"

    package require platform
    set tcl_platform [platform::generic]
    if { $tcl_platform == "linux-x86_64" } {
        # SphereTree should be installed the linux system (compiled and installed)
        # set program makeTreeMedial
        set program [file join $::DEMClusters::ProblemTypePath exec makeTreeMedial]
    } else {
        set program [file join $::DEMClusters::ProblemTypePath exec $Algorithm]
    }

    # exec $program {*}$argv
    # catch { set ::DEMClusters::pid [exec $program {*}$argv] &} msg

    # proces can be canceled but not ouput on process info screen
    catch {set ::DEMClusters::pid [open "|$program $argv"]} msg

    # set ouputpath [file join $::DEMClusters::ProblemPath $::DEMClusters::ProblemName.info]
    # set outfl [open $ouputpath w]

    #puts $outfl [read $::DEMClusters::pid]
    #flush stdout
    # close $outfl

    #W [read $::DEMClusters::pid]
    #W "SPH generation finished"


    # view output using view process info but process cannot be canceled. requires tk_exec modified version of exec proc:
    # set ouputpath [file join $::DEMClusters::ProblemPath $::DEMClusters::ProblemName.info]
    # set outfl [open $ouputpath w]
    # catch {tk_exec $program {*}$argv > $ouputpath} msg

    # puts $outfl [set ::DEMClusters::pid [open "|$program $argv"]]
    # close $outfl


    # status will be= ok|user_stop|timeout
    # proc my_finish_proc { pid status args } { ... }
    # returned user_stop must be 0 or 1 (1 to kill the process)
    # proc my_progress_proc { args }  { ... return $user_stop }
    # set timeout 20000 ;#seconds, 0 to no limit
    # set maxmemory 0 ;#bytes 0 to no limit
    # gid_cross_platform::track_process $::DEMClusters::pid 1000 [clock seconds] $timeout $maxmemory [list my_finish_proc $args] [list my_progress_proc a b]




    # makeTreeMedial -depth 1 -branch 100 -numCover 10000 -minCover 5 -initSpheres 1000 -minSpheres 200 -erFact 2 -testerLevels 2 -nopause -eval -expand -merge -burst -optimise balance -balExcess 0.001 -maxOptLevel 100 generic.obj
    # set program [lindex $argv 0]
    # set program [file join $::DEMClusters::ProblemTypePath exec MakeTreeMedial.exe]
    # set arguments [lrange $argv 1 end]
    # exec $program {*}$arguments

    # TreeMedial ValidArgs:
    # -depth              Depth of the sphere-tree
    # -branch             Branching factor of sphere-tree - for depth 1, branch=num of spheres
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
    set depth [GiD_AccessValue get gendata depth]
    set numCover [GiD_AccessValue get gendata numCover]
    set minCover [GiD_AccessValue get gendata minCover]
    set testerLevels [GiD_AccessValue get gendata testerLevels]
    #set genericOBJFilename [file join $::DEMClusters::ProblemPath generic.obj]
    #set genericOBJFilename "\"$genericOBJFilename\""
    set modelname [GiD_Info Project ModelName]
    set genericOBJFilename [file join ${modelname}.gid generic.obj]

    set argv "-depth $depth -numCover $numCover -minCover $minCover -testerLevels $testerLevels -nopause -eval $genericOBJFilename"
    package require platform
    set tcl_platform [platform::generic]
    if { $tcl_platform == "linux-x86_64" } {
        # SphereTree should be installed the linux system (compiled and installed)
        # set program makeTreeGrid
        set program [file join $::DEMClusters::ProblemTypePath exec makeTreeGrid]
    } else {
        set program [file join $::DEMClusters::ProblemTypePath exec $Algorithm]
    }

    #exec $program {*}$argv
    catch {set ::DEMClusters::pid [open "|$program $argv"]} res

    # makeTreeGrid ValidArgs:
    # -depth              Depth of the sphere-tree
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
    set depth [GiD_AccessValue get gendata depth]
    set numCover [GiD_AccessValue get gendata numCover]
    set minCover [GiD_AccessValue get gendata minCover]
    set testerLevels [GiD_AccessValue get gendata testerLevels]
    #set genericOBJFilename [file join $::DEMClusters::ProblemPath generic.obj]
    #set genericOBJFilename "\"$genericOBJFilename\""
    set modelname [GiD_Info Project ModelName]
    set genericOBJFilename [file join ${modelname}.gid generic.obj]

    set argv "-depth $depth -numCover $numCover -minCover $minCover -testerLevels $testerLevels -nopause -eval $genericOBJFilename"
    package require platform
    set tcl_platform [platform::generic]
    if { $tcl_platform == "linux-x86_64" } {
        # SphereTree should be installed the linux system (compiled and installed)
        # set program makeTreeSpawn
        set program [file join $::DEMClusters::ProblemTypePath exec makeTreeSpawn]
    } else {
        set program [file join $::DEMClusters::ProblemTypePath exec $Algorithm]
    }

    #exec $program {*}$argv
    catch {set ::DEMClusters::pid [open "|$program $argv"]} res

    # makeTreeSpawn ValidArgs:
    # -depth              Depth of the sphere-tree
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



proc DEMClusters::call_makeTreeHubbard { } {

    set Algorithm [GiD_AccessValue get gendata Algorithm]
    set branch [GiD_AccessValue get gendata branch]
    set depth [GiD_AccessValue get gendata depth]
    set numSamples [GiD_AccessValue get gendata numSamples]
    set minSamples [GiD_AccessValue get gendata minSamples]
    #set genericOBJFilename [file join $::DEMClusters::ProblemPath generic.obj]
    #set genericOBJFilename "\"$genericOBJFilename\""
    set modelname [GiD_Info Project ModelName]
    set genericOBJFilename [file join ${modelname}.gid generic.obj]

    set argv "-depth $depth -branch $branch -numSamples $numSamples -minSamples $minSamples -nopause $genericOBJFilename"
    package require platform
    set tcl_platform [platform::generic]
    if { $tcl_platform == "linux-x86_64" } {
        # SphereTree should be installed the linux system (compiled and installed)
        # set program makeTreeHubbard
        set program [file join $::DEMClusters::ProblemTypePath exec makeTreeHubbard]
    } else {
        set program [file join $::DEMClusters::ProblemTypePath exec $Algorithm]
    }

    #exec $program {*}$argv
    catch {set ::DEMClusters::pid [open "|$program $argv"]} res

    # makeTreeHubbard ValidArgs:
    # -depth              Depth of the sphere-tree
    # -branch             Branching factor of sphere-tree
    # -numSamples         Number of sample points to cover object with
    # -minSamples         Minimum number of sample points per triangle
    # makeTreeHubbard  -branch 8 -depth 3 -numSamples 500 -minSamples 1 -nopause generic.obj
}

proc DEMClusters::call_makeTreeOctree { } {

    set Algorithm [GiD_AccessValue get gendata Algorithm]
    set depth [GiD_AccessValue get gendata depth]
    #set genericOBJFilename [file join $::DEMClusters::ProblemPath generic.obj]
    #set genericOBJFilename "\"$genericOBJFilename\""
    set modelname [GiD_Info Project ModelName]
    set genericOBJFilename [file join ${modelname}.gid generic.obj]

    set argv "-depth $depth -nopause $genericOBJFilename"
    package require platform
    set tcl_platform [platform::generic]
    if { $tcl_platform == "linux-x86_64" } {
        # SphereTree should be installed the linux system (compiled and installed)
        # set program makeTreeOctree
        set program [file join $::DEMClusters::ProblemTypePath exec makeTreeOctree]
    } else {
        set program [file join $::DEMClusters::ProblemTypePath exec $Algorithm]
    }

    #exec $program {*}$argv
    catch {set ::DEMClusters::pid [open "|$program $argv"]} res

    # makeTreeOctree ValidArgs:
    # -depth              Depth of the sphere-tree
    # -nopause            Don't pause when processing, i.e. batch mode
}


proc CenterGeometry { } {

    # last component of the list
    set last_volume [lrange [GiD_Geometry list volume 1:] end end]

    set all [GiD_Tools geometry mass_properties $last_volume]
    # Only GID 14.9d+

    # set mass [lrange $all 0 0]
    set a [join [lrange $all 1 1]]
    set centerX [lindex $a 0]
    set centerY [lindex $a 1]
    set centerZ [lindex $a 2]

    set nodeslist [GiD_Geometry list point 1:]

    W [GiD_Info listmassproperties Points $nodeslist]

    GiD_Process Mescape Utilities Move Volumes MaintainLayers Translation FNoJoin $centerX,$centerY,$centerZ FNoJoin 0.0,0.0,0.0 $last_volume escape Mescape
}



proc CenterMesh { } {

    set tetrahedra [GiD_Mesh list -element_type {tetrahedra} element]
    set triangles [GiD_Mesh list -element_type {triangle} element]
    # set all [GiD_Tools mesh mass_properties $tetrahedra]
    set all [GiD_Tools mesh mass_properties -boundary_elements $triangles]
    # Only GID 14.9d+

    set cdg [join [lrange $all 1 1]]
    set centerX [lindex $cdg 0]
    set centerY [lindex $cdg 1]
    set centerZ [lindex $cdg 2]

    GiD_Process Mescape Utilities Move Elements MaintainLayers Translation FNoJoin $centerX,$centerY,$centerZ FNoJoin 0.0,0.0,0.0 {*}$tetrahedra escape Mescape
    GiD_Process Mescape Utilities Move Elements MaintainLayers Translation FNoJoin $centerX,$centerY,$centerZ FNoJoin 0.0,0.0,0.0 {*}$triangles escape Mescape

}


proc GenerateClusterFile { } {

    # set all [GiD_Tools geometry mass_properties 1]
    # Only GID 14.9d+

    # set mass [lrange $all 0 0]

    # set a [join [lrange $all 1 1]]
    # set centerX [lindex $a 0]
    # set centerY [lindex $a 1]
    # set centerZ [lindex $a 2]

    # set a [join [lrange $all 2 2]]
    # set Ixx [expr {[lindex $a 0]/$mass}]
    # set Iyy [expr {[lindex $a 1]/$mass}]
    # set Izz [expr {[lindex $a 2]/$mass}]


    set Algorithm [GiD_AccessValue get gendata Algorithm]
    if {$Algorithm == "MakeTreeMedial"} {
        set genericSPHFilename [file join $::DEMClusters::ProblemPath generic-medial.sph]
        #set genericSPHFilename "\"$genericSPHFilename\""
    } elseif {$Algorithm == "MakeTreeGrid"} {
        set genericSPHFilename [file join $::DEMClusters::ProblemPath generic-grid.sph]
        #set genericSPHFilename "\"$genericSPHFilename\""
    } elseif {$Algorithm == "MakeTreeSpawn"} {
        set genericSPHFilename [file join $::DEMClusters::ProblemPath generic-spawn.sph]
        #set genericSPHFilename "\"$genericSPHFilename\""
    } elseif {$Algorithm == "MakeTreeOctree"} {
        set genericSPHFilename [file join $::DEMClusters::ProblemPath generic-octree.sph]
        #set genericSPHFilename "\"$genericSPHFilename\""
    } elseif {$Algorithm == "MakeTreeHubbard"} {
        #set genericSPHFilename generic-hubbard.sph
        set genericSPHFilename [file join $::DEMClusters::ProblemPath generic-hubbard.sph]
        #set genericSPHFilename "\"$genericSPHFilename\""
    } else {
        W "Select a valid algorithm"
    }

    # warning: this never ends. Only works with one click go.
    set i 0
    while { [file exists $genericSPHFilename] == 0 } {
        W "waiting for SPH to be generated... time: $i seconds"
        after 30000
        set i [expr {$i + 30}]
        # wait 30 seconds and recheck if sph has been generated.
    }

    W "SPH file detected in folder. Generating CLU file."
    set genericMSHFilename [file join $::DEMClusters::ProblemPath generic.msh]
    set ouputpath [file join $::DEMClusters::ProblemPath generic_cluster.clu]
    W "Cluster file generated in the following location:"
    W $ouputpath
    # set genericMSHFilename "\"$genericMSHFilename\""
    set argv_number 3
    package require platform
    set tcl_platform [platform::generic]
    if { $tcl_platform == "linux-x86_64" } {
        set cluster_exec create_cluster_unix
        set program [file join $::DEMClusters::ProblemTypePath exec $cluster_exec]
    } else {
        set cluster_exec create_cluster
        set program [file join $::DEMClusters::ProblemTypePath exec $cluster_exec]
    }

    exec $program 3 $genericMSHFilename $genericSPHFilename $ouputpath
    # DEBUG : exec $program $argv_number $genericMSHFilename $genericSPHFilename > [file join $::DEMClusters::ProblemPath output.txt]
    # REGENERATE EXEC FOR WINDOWS including required libraries: exec> g++ -static-libgcc -static-libstdc++ -o create_cluster.exe mesh_to_clu_converter.cpp
    # REGENERATE EXEC FOR unix including required libraries: exec> g++ -static-libgcc -static-libstdc++ -o create_cluster_unix mesh_to_clu_converter.cpp
}


proc AlignSurfNormals {direction} {
    # ABSTRACT: Makes all of boundary surfaces' normals point inwards or outwards
    # Arguments
    # direction => Direction option ["Inwards"|"Outwards"]
    # Note: This procedure in the same used in the fluid_only problem type

    switch $direction {
        Inwards {
            set wrong_way "DIFF1ST"
        }
        Outwards {
            set wrong_way "SAME1ST"
        }
        default {puts "Unknown Direction, surface normals not aligned"}
    }
    set volumelist [GiD_Geometry list volume 1:]
    set surfacelist [list]

    # For each volume, we look for face surfaces with oriented in the wrong direction
    foreach volume $volumelist {
        set volumeinfo [GiD_Info list_entities volumes $volume]
        set numpos [lsearch $volumeinfo "NumSurfaces:"]
        set numsurf [lindex $volumeinfo [expr {$numpos +1 }]]
        for {set i 0} {$i < $numsurf} {incr i} {
            set orient [lindex $volumeinfo [expr {$numpos+5+4*$i}]]
            if {[string compare $orient $wrong_way]==0} {
                # If the normal is pointing in the wrong direction,
                # Check if it's a contour surface
                set surfnum [lindex $volumeinfo [expr {$numpos+3+4*$i}]]
                set surfinfo [GiD_Info list_entities surfaces $surfnum]
                set higherentities [lindex $surfinfo 4]
                if {$higherentities==1} {
                    lappend surfacelist $surfnum
                }
            }
        }
    }

    if {[llength $surfacelist]} {
        # If its in the contour, switch its normal
        eval GiD_Process Mescape Utilities SwapNormals Surfaces Select $surfacelist
    }
}



proc ReadClusterFileintoGeometry { } {

    set modelname [GiD_Info Project ModelName]
    set clustername [file join ${modelname}.gid generic_cluster.clu]
    set FileVar [open $clustername "r+"]
    set file_data [read $FileVar]
    lreplace $file_data end end
    close $FileVar

    #  Process data file
    set data [split $file_data "\n"]
    set data [lreplace $data end-14 end]
    foreach line $data {
        if {[llength $line] >3} {
            GiD_Process Mescape Geometry Create Object Sphere [lrange $line 0 0] [lrange $line 1 1] [lrange $line 2 2] [lrange $line 3 3] escape escape
        }
    }
}






proc ReadClusterFileintoMesh { } {

    set modelname [GiD_Info Project ModelName]
    set clustername [file join ${modelname}.gid generic_cluster.clu]
    set FileVar [open $clustername "r+"]
    set file_data [read $FileVar]
    lreplace $file_data end end
    close $FileVar

    #  Process data file
    set data [split $file_data "\n"]
    set data [lreplace $data end-14 end]
    set sphere_nodes [list]
    set i 0
    foreach line $data {
        if {[llength $line] >3} {
       		set x [lindex $line 0]
            set y [lindex $line 1]
            set z [lindex $line 2]
            set node_id [GiD_Mesh create node append  [list $x $y $z]] ;
            lappend sphere_nodes {*}$node_id ;
            incr i 1
        }
    }
    set count 0
    foreach line $data {
        if {[llength $line] >3} {
            GiD_Process Mescape Meshing EditMesh CreateElement Sphere [lindex $line 3] [lindex $sphere_nodes $count] escape escape
            incr count 1
        }
    }
}


proc ReadSPHFileintoMesh { } {

    set Algorithm [GiD_AccessValue get gendata Algorithm]
    if {$Algorithm == "MakeTreeMedial"} {
        set genericSPHFilename [file join $::DEMClusters::ProblemPath generic-medial.sph]
        #set genericSPHFilename "\"$genericSPHFilename\""
    } elseif {$Algorithm == "MakeTreeGrid"} {
        set genericSPHFilename [file join $::DEMClusters::ProblemPath generic-grid.sph]
        #set genericSPHFilename "\"$genericSPHFilename\""
    } elseif {$Algorithm == "MakeTreeSpawn"} {
        set genericSPHFilename [file join $::DEMClusters::ProblemPath generic-spawn.sph]
        #set genericSPHFilename "\"$genericSPHFilename\""
    } elseif {$Algorithm == "MakeTreeOctree"} {
        set genericSPHFilename [file join $::DEMClusters::ProblemPath generic-octree.sph]
        #set genericSPHFilename "\"$genericSPHFilename\""
    } elseif {$Algorithm == "MakeTreeHubbard"} {
        #set genericSPHFilename generic-hubbard.sph
        set genericSPHFilename [file join $::DEMClusters::ProblemPath generic-hubbard.sph]
        #set genericSPHFilename "\"$genericSPHFilename\""
    } else {
        W "Select a valid algorithm"
    }

    #set modelname [GiD_Info Project ModelName]
    #set clustername [file join ${modelname}.gid generic_cluster.clu]
    set FileVar [open $genericSPHFilename "r+"]
    set file_data [read $FileVar]
    lreplace $file_data end end
    close $FileVar

    #  Process data file
    set data [split $file_data "\n"]
    set data [lreplace $data 0 1]
    set data [lreplace $data end-36 end]
    set sphere_nodes [list]
    set i 0
    foreach line $data {
        if {[llength $line] >3} {
       		set x [lindex $line 0]
            set y [lindex $line 1]
            set z [lindex $line 2]
            if {[lindex $line 0] != 0.0} {
                set node_id [GiD_Mesh create node append  [list $x $y $z]] ;
                lappend sphere_nodes {*}$node_id ;
                incr i 1
            }
        }
    }
    set count 0
    foreach line $data {
        if {[llength $line] >3} {
            if {[lindex $line 0] != 0.0} {
                GiD_Process Mescape Meshing EditMesh CreateElement Sphere [lindex $line 3] [lindex $sphere_nodes $count] escape escape
                incr count 1
            }
        }
    }
}



 proc tk_exec_fileevent {id} {
    global tk_exec_data
    global tk_exec_cond
    global tk_exec_pipe

    if {[eof $tk_exec_pipe($id)]} {
        fileevent $tk_exec_pipe($id) readable ""
        set tk_exec_cond($id) 1
        return
    }

    append tk_exec_data($id) [read $tk_exec_pipe($id) 1024]
  }

  proc tk_exec {args} {
    global tk_exec_id
    global tk_exec_data
    global tk_exec_cond
    global tk_exec_pipe
    global tcl_platform
    global env

    if {![info exists tk_exec_id]} {
        set tk_exec_id 0
    } else {
        incr tk_exec_id
    }

    set keepnewline 0

    for {set i 0} {$i < [llength $args]} {incr i} {
        set arg [lindex $args $i]
        switch -glob -- $arg {
            -keepnewline {
                set keepnewline 1
            }
            -- {
                incr i
                break
            }
            -* {
                error "unknown option: $arg"
            }
            * {
                break
            }
        }
    }

    if {$i > 0} {
        set args [lrange $args $i end]
    }

    if {$tcl_platform(platform) == "windows" && \
        [info exists env(COMSPEC)]} {
        set args [linsert $args 0 $env(COMSPEC) "/c"]
    }

    set pipe [open "|$args" r]

    set tk_exec_pipe($tk_exec_id) $pipe
    set tk_exec_data($tk_exec_id) ""
    set tk_exec_cond($tk_exec_id) 0

    fconfigure $pipe -blocking 0
    fileevent $pipe readable "tk_exec_fileevent $tk_exec_id"

    vwait tk_exec_cond($tk_exec_id)

    if {$keepnewline} {
        set data $tk_exec_data($tk_exec_id)
    } else {
        set data [string trimright $tk_exec_data($tk_exec_id) \n]
    }

    unset tk_exec_pipe($tk_exec_id)
    unset tk_exec_data($tk_exec_id)
    unset tk_exec_cond($tk_exec_id)

    if {[catch {close $pipe} err]} {
        error "pipe error: $err"
    }

    return $data
  }
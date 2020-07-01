
### Brief methodology to generate the cluster file ###
# 1. extract from mesh: OBJ from surface mesh, MSH from the tetrahedra mesh
# 2. Generate SPH via sphereTree external executable passing OBJ file as argument (current issue accessing gid values)
# 3. Clean the SPH file
# 4. Generate Cluster.clu via mesh_to_cluster_converter.cpp passing SPH + MSH as arguments
#    modified and precompiled executable will be required to do this


## Extra features:
## - Lateral icons
## - Show progress onscreen while generating SPH via spheretree
## - Cancel button to stop external execs
## - Remove finished calculate splash.( or modify with OBJ created succesfully)


##  --------------------------------------------------------------------------------------------------------------------------------------------------
##  Current Issues ##


##  --------------------------------------------------------------------------------------------------------------------------------------------------
##  Fixed Issues ##

##- Are ALL calculated Vertex normals correct? verified manually with testcubev4
##- Sphere tree paths with spaces
##- Dependencies in prb.  each line for an algorithm, one parenthesis
##- verify that your geometry have all the normals coherent. (should be as align normals has been added before meshing)
##- add some info on the help based on the manual
##- Calling external precompiled cpp, add criteria for negative radius when deleting line in sph
##- plan on msh to clu functionaly design of precompiled executable. args, paths and location,
##-     define arguments path and call exe from inside exec folder.
##- add dummy .bat to avoid error showing both unix.bat and win.bat
##- cluster visualizaton in GID pre
##- add export gidmesh as generic.msh on calculate
##- RECOMMENEDED MEDIAL - spheretree throws error if algorithm parameters are not quite good. example: small geom with high numsamples
## - Cluster visualization
##  --------------------------------------------------------------------------------------------------------------------------------------------------


## GiD events --------------------------------------------------------------------------------------------------------------------------------------------------


proc InitGIDProject { dir } {

    # Create buttons to insert algorithm values
    if { [GidUtils::IsTkDisabled] eq 0} {
        GiDMenu::Create "Sphere Cluster Creation" PRE
        #GiDMenu::InsertOption "Sphere Cluster Creation" [list "SphereTree"] 0 PRE "GidOpenConditions \"SphereTree\"" "" ""
        GiDMenu::InsertOption "Sphere Cluster Creation" [list "Define Options"] 0 PRE "GidOpenProblemData" "" ""
        GiDMenu::InsertOption "Sphere Cluster Creation" [list "Generate SPH file" ] 1 PRE [list GenerateSPHFileFromOBJFile] "" ""
        GiDMenu::InsertOption "Sphere Cluster Creation" [list "Generate CLU file" ] 2 PRE [list GenerateClusterFile] "" ""
        GiDMenu::InsertOption "Sphere Cluster Creation" [list "Generate cluster and visualize" ] 3 PRE [list OneClickGo] "" ""
        #GiDMenu::InsertOption "Sphere Cluster Creation" [list "Visualize cluster over geometry" ] 4 PRE [list ReadClusterFileintoGeometry] "" ""
        GiDMenu::InsertOption "Sphere Cluster Creation" [list "Visualize cluster over mesh" ] 4 PRE [list ReadClusterFileintoMesh] "" ""
        #GiDMenu::InsertOption "Sphere Cluster Creation" [list "Delete cluster over geometry" ] 6 PRE [list DeleteSpheresGeometry] "" ""
        GiDMenu::InsertOption "Sphere Cluster Creation" [list "Delete cluster over mesh" ] 5 PRE [list DeleteSpheresMesh] "" ""

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

    ::GidUtils::Splash [file join $::DEMClusters::ProblemTypePath splash.png] .splash 1

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
    GiD_Mesh delete element $spheres
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
    ReadClusterFileintoMesh
}


proc BeforeRunCalculation { batfilename basename dir problemtypedir gidexe args } {

    source [file join $problemtypedir OBJFile.tcl]
    set OBJOutput [GenerateOBJFile $basename $dir $problemtypedir]

    set modelname [GiD_Info Project ModelName]
    set export_msh [file join ${modelname}.gid generic.msh]
    GiD_Process Mescape Files WriteMesh $export_msh
}

proc AfterRunCalculation {basename dir problemtypedir where error errorfilename} {
    return 0
}


proc GenerateSPHFileFromOBJFile { } {
    call_SphereTree
}


namespace eval DEMClusters {
    variable ProblemName ""
    variable ProblemPath ""
    variable ProblemTypePath ""
}

proc call_SphereTree { } {

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
}

proc DEMClusters::call_TreeMedial { } {

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

    set argv "-depth $depth -branch $branch -numCover $numCover -minCover $minCover -initSpheres $initSpheres -minSpheres $minSpheres -erFact $erFact -testerLevels $testerLevels -verify -nopause -eval -expand -merge -burst -optimise balance -balExcess 0.001 -maxOptLevel 100 $genericOBJFilename"

    package require platform
    set tcl_platform [platform::generic]
    if { $tcl_platform == "linux-x86_64" } {
        # SphereTree should be installed the linux system (compiled and installed)
        set program makeTreeMedial
    } else {
        set program [file join $::DEMClusters::ProblemTypePath exec $Algorithm]
    }
    exec $program {*}$argv

    # makeTreeMedial -depth 1 -branch 100 -numCover 10000 -minCover 5 -initSpheres 1000 -minSpheres 200 -erFact 2 -testerLevels 2 -verify -nopause -eval -expand -merge -burst -optimise balance -balExcess 0.001 -maxOptLevel 100 generic.obj
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

    set argv "-depth $depth -branch $branch -numCover $numCover -minCover $minCover -testerLevels $testerLevels -verify -nopause -eval $genericOBJFilename"
    package require platform
    set tcl_platform [platform::generic]
    if { $tcl_platform == "linux-x86_64" } {
        # SphereTree should be installed the linux system (compiled and installed)
        set program makeTreeGrid
    } else {
        set program [file join $::DEMClusters::ProblemTypePath exec $Algorithm]
    }

    exec $program {*}$argv

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

    set argv "-depth $depth -branch $branch -numCover $numCover -minCover $minCover -testerLevels $testerLevels -verify -nopause -eval $genericOBJFilename"
    package require platform
    set tcl_platform [platform::generic]
    if { $tcl_platform == "linux-x86_64" } {
        # SphereTree should be installed the linux system (compiled and installed)
        set program makeTreeSpawn
    } else {
        set program [file join $::DEMClusters::ProblemTypePath exec $Algorithm]
    }

    exec $program {*}$argv

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

    set modelname [GiD_Info Project ModelName]
    set genericOBJFilename [file join ${modelname}.gid generic.obj]

    set argv "-depth $depth -branch $branch -numSamples $numSamples -minSamples $minSamples -nopause $genericOBJFilename"
    package require platform
    set tcl_platform [platform::generic]
    if { $tcl_platform == "linux-x86_64" } {
        # SphereTree should be installed the linux system (compiled and installed)
        set program makeTreeHubbard
    } else {
        set program [file join $::DEMClusters::ProblemTypePath exec $Algorithm]
    }

    exec $program {*}$argv

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
        set program makeTreeOctree
    } else {
        set program [file join $::DEMClusters::ProblemTypePath exec $Algorithm]
    }


    exec $program {*}$argv

    # makeTreeOctree ValidArgs:
    # -depth              Depth of the sphere-tree
    # -nopause            Don't pause when processing, i.e. batch mode
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

    set genericMSHFilename [file join $::DEMClusters::ProblemPath generic.msh]
    set ouputpath [file join $::DEMClusters::ProblemPath generic_cluster.clu]
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
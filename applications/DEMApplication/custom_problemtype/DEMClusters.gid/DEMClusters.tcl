## GiD events --------------------------------------------------------------------------------------------------------------------------------------------------

proc InitGIDProject { dir } {

    # Create buttons for calculate mass, center of mass and inertia later on.
    # if { [GidUtils::IsTkDisabled] eq 0} {
    #     GiDMenu::Create "Sphere Cluster Creation" PRE
    #     GiDMenu::InsertOption "Sphere Cluster Creation" [list "Inertia"] 0 PRE "GidOpenConditions \"Inertia\"" "" ""
    #     GiDMenu::UpdateMenus
    # }

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
    CalculateInertia
}


proc BeforeRunCalculation { batfilename basename dir problemtypedir gidexe args } {

    # Write file
    source [file join $problemtypedir file.tcl]
}

proc GenerateOBJFile { } {
    # Analyze the format of the OBJ and generate the file in GID
}

proc ExecuteSphereTreeAlgorithm { objfilename args} {
    #makeTreeMedial -branch NS -depth 1 -testerLevels 2 -numCover 10000 -minCover 5 -initSpheres 1000 -minSpheres 200 -erFact 2 -verify -nopause -eval -expand -merge -burst -optimise balance -balExcess 0.001 -maxOptLevel 100 file_name.obj
    # multiple args to be defined: algorithm,...

}

proc CleanSPHFile { sphfilename } {
    # Delete some lines and columns. 
    # Follow criteria to delete invalid spheres

}

proc GenerateClusterFile { sphfilename mshfilename} {
    # Look for a way to edit, compile and execute mesh_to_clu_converter.cpp with the custom filenames
    # Create cpp executable able to generate cluster file directly from inputs.
    # Or pass always the same generic sph and msh names. Generate generic cluster filename.

}


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

#-------------------------------------------------------------------------------

proc BeforeMeshGeneration {elementsize} {
    W "execute BeforeMeshGeneration"

}

proc AfterMeshGeneration {fail} {
    W "execute AfterMeshGeneration"
    CalculateInertia
}

#-------------------------------------------------------------------------------

proc BeforeRunCalculation { batfilename basename dir problemtypedir gidexe args } {

    # Write file
    source [file join $problemtypedir file.tcl]
}

#-------------------------------------------------------------------------------

proc CalculateInertia  { } { 
    W "execute CalculateInertia"
    foreach volume_id [GiD_Geometry list volume 1:end] {
        set volume_info [GiD_Info list_entities volume $volume_id]
    }
    
    set triangles [GiD_Mesh list -element_type {triangle} element]
    W $triangles

    set tetra [GiD_Mesh list -element_type {tetrahedra} element]
    W $tetra



    # 1.- calculate total volume from tetrahedras:
    # tris: number of triangles
    #     for each triangle
    #     triple producto de los vertices de cada uno.
    #     volume += Determinant(float3x3(vertices[tris[i][0]],vertices[tris[i][1]],vertices[tris[i][2]]));
    # return volume/6.0;  // since the determinant give 6 times tetra volume


}

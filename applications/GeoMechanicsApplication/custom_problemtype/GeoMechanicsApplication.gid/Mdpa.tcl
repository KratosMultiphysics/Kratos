proc WriteMdpa { basename dir problemtypedir } {

    ## Source auxiliar procedures
    source [file join $problemtypedir MdpaAuxProcs.tcl]

    ## Start MDPA file
    set filename [file join $dir ${basename}.mdpa]
    set FileVar [open $filename w]

    ## ModelPart Data
    #puts $FileVar "Begin ModelPartData"
    #puts $FileVar "  // VARIABLE_NAME value"
    #puts $FileVar "End ModelPartData"
    #puts $FileVar ""
    #puts $FileVar ""
    
    
    ## Tables
    set TableId 0
    set TableDict [dict create]
    # Solid_Displacement
    ConstraintVectorTable FileVar TableId TableDict Solid_Displacement DISPLACEMENT
    # Structural_Rotation
    ConstraintVectorTable FileVar TableId TableDict Structural_Rotation ROTATION
    # Fluid_Pressure
    PressureTable FileVar TableId TableDict Fluid_Pressure WATER_PRESSURE
    # Point_Load
    VectorTable FileVar TableId TableDict Point_Load POINT_LOAD
    # Line_Load
    VectorTable FileVar TableId TableDict Line_Load LINE_LOAD
    # Surface_Load
    VectorTable FileVar TableId TableDict Surface_Load SURFACE_LOAD
    # Normal_Load
    NormalTangentialTable FileVar TableId TableDict Normal_Load NORMAL_CONTACT_STRESS TANGENTIAL_CONTACT_STRESS
    # Normal_Fluid_Flux
    ScalarTable FileVar TableId TableDict Normal_Fluid_Flux NORMAL_FLUID_FLUX
    # Interface_Face_Load
    VectorTable FileVar TableId TableDict Interface_Face_Load FACE_LOAD
    # Interface_Normal_Fluid_Flux
    ScalarTable FileVar TableId TableDict Interface_Normal_Fluid_Flux NORMAL_FLUID_FLUX
    # Body_Acceleration
    VectorTable FileVar TableId TableDict Body_Acceleration VOLUME_ACCELERATION
    puts $FileVar ""


    ## Properties
    set PropertyId 0
    set PropertyDict [dict create]
    # Soil_two_phase part
    set Groups [GiD_Info conditions Soil_two_phase groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
        incr PropertyId
        dict set PropertyDict [lindex [lindex $Groups $i] 1] $PropertyId
        puts $FileVar "Begin Properties $PropertyId"
        puts $FileVar "End Properties"
        puts $FileVar ""
    }

    # Soil_drained part
    set Groups [GiD_Info conditions Soil_drained groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
        incr PropertyId
        dict set PropertyDict [lindex [lindex $Groups $i] 1] $PropertyId
        puts $FileVar "Begin Properties $PropertyId"
        puts $FileVar "End Properties"
        puts $FileVar ""
    }

    # Soil_undrained part
    set Groups [GiD_Info conditions Soil_undrained groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
        incr PropertyId
        dict set PropertyDict [lindex [lindex $Groups $i] 1] $PropertyId
        puts $FileVar "Begin Properties $PropertyId"
        puts $FileVar "End Properties"
        puts $FileVar ""
    }

    # Non_porous part
    set Groups [GiD_Info conditions Non_porous groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
        incr PropertyId
        dict set PropertyDict [lindex [lindex $Groups $i] 1] $PropertyId
        puts $FileVar "Begin Properties $PropertyId"
        puts $FileVar "End Properties"
        puts $FileVar ""
    }

    # Soil_Groundwater_Flow part
    set Groups [GiD_Info conditions Soil_Groundwater_Flow groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
        incr PropertyId
        dict set PropertyDict [lindex [lindex $Groups $i] 1] $PropertyId
        puts $FileVar "Begin Properties $PropertyId"
        puts $FileVar "End Properties"
        puts $FileVar ""
    }

    # Interface_Groundwater_flow part
    set Groups [GiD_Info conditions Interface_Groundwater_flow groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
        incr PropertyId
        dict set PropertyDict [lindex [lindex $Groups $i] 1] $PropertyId
        puts $FileVar "Begin Properties $PropertyId"
        puts $FileVar "End Properties"
        puts $FileVar ""
    }

    # Beam part
    set Groups [GiD_Info conditions Beam groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
            incr PropertyId
            dict set PropertyDict [lindex [lindex $Groups $i] 1] $PropertyId
            puts $FileVar "Begin Properties $PropertyId"
            puts $FileVar "End Properties"
            puts $FileVar ""
    }

    # Shell_thin_corotational part
    set Groups [GiD_Info conditions Shell_thin_corotational groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
            incr PropertyId
            dict set PropertyDict [lindex [lindex $Groups $i] 1] $PropertyId
            puts $FileVar "Begin Properties $PropertyId"
            puts $FileVar "End Properties"
            puts $FileVar ""
    }

    # Shell_thick_corotational part
    set Groups [GiD_Info conditions Shell_thick_corotational groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
            incr PropertyId
            dict set PropertyDict [lindex [lindex $Groups $i] 1] $PropertyId
            puts $FileVar "Begin Properties $PropertyId"
            puts $FileVar "End Properties"
            puts $FileVar ""
    }

    # Truss part
    set Groups [GiD_Info conditions Truss groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
            incr PropertyId
            dict set PropertyDict [lindex [lindex $Groups $i] 1] $PropertyId
            puts $FileVar "Begin Properties $PropertyId"
            puts $FileVar "End Properties"
            puts $FileVar ""
    }

    # Anchor part
    set Groups [GiD_Info conditions Anchor groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
            incr PropertyId
            dict set PropertyDict [lindex [lindex $Groups $i] 1] $PropertyId
            puts $FileVar "Begin Properties $PropertyId"
            puts $FileVar "End Properties"
            puts $FileVar ""
    }

    # Interface two_phase, drained and undrained parts
    set interface_Groups [list [GiD_Info conditions Interface_two_phase groups] [GiD_Info conditions Interface_drained groups] [GiD_Info conditions Interface_undrained groups]]
    foreach Groups $interface_Groups {
        for {set i 0} {$i < [llength $Groups]} {incr i} {
            incr PropertyId
            dict set PropertyDict [lindex [lindex $Groups $i] 1] $PropertyId
            puts $FileVar "Begin Properties $PropertyId"
            puts $FileVar "End Properties"
            puts $FileVar ""
        }
    }
    
    puts $FileVar ""

    ## Nodes
    set Nodes [GiD_Info Mesh Nodes]
    puts $FileVar "Begin Nodes"
    for {set i 0} {$i < [llength $Nodes]} {incr i 4} {
        # puts $FileVar "  [lindex $Nodes $i]  [lindex $Nodes [expr { $i+1 }]] [lindex $Nodes [expr { $i+2 }]] [lindex $Nodes [expr { $i+3 }]]"
        puts -nonewline $FileVar "  [lindex $Nodes $i]  "
        puts -nonewline $FileVar [format  "%.10f" [lindex $Nodes [expr { $i+1 }]]]
        puts -nonewline $FileVar " "
        puts -nonewline $FileVar [format  "%.10f" [lindex $Nodes [expr { $i+2 }]]]
        puts -nonewline $FileVar " "
        puts $FileVar [format  "%.10f" [lindex $Nodes [expr { $i+3 }]]]
    }
    puts $FileVar "End Nodes"
    puts $FileVar ""
    puts $FileVar ""


    ## Elements
    set FIC [GiD_AccessValue get gendata FIC_Stabilization]
    set IsQuadratic [GiD_Info Project Quadratic]
    set Dim [GiD_AccessValue get gendata Domain_Size]
    set IsMoveMesh [GiD_AccessValue get gendata Move_Mesh]
    set SolutionType [GiD_AccessValue get gendata Solution_Type]

    # Soil_two_phase part
    set Groups [GiD_Info conditions Soil_two_phase groups]
    if {$IsQuadratic eq 0} {
        if {$FIC eq false} {
            for {set i 0} {$i < [llength $Groups]} {incr i} {
                # Elements Property
                set BodyElemsProp [dict get $PropertyDict [lindex [lindex $Groups $i] 1]]

                if {$IsMoveMesh eq false} {
                    # UPwSmallStrainElement2D3N
                    WriteElements FileVar [lindex $Groups $i] triangle UPwSmallStrainElement2D3N $BodyElemsProp Triangle2D3Connectivities
                    # UPwSmallStrainElement2D4N
                    WriteElements FileVar [lindex $Groups $i] quadrilateral UPwSmallStrainElement2D4N $BodyElemsProp Quadrilateral2D4Connectivities
                    # UPwSmallStrainElement3D4N
                    WriteElements FileVar [lindex $Groups $i] tetrahedra UPwSmallStrainElement3D4N $BodyElemsProp Quadrilateral2D4Connectivities
                    # UPwSmallStrainElement3D8N
                    WriteElements FileVar [lindex $Groups $i] hexahedra UPwSmallStrainElement3D8N $BodyElemsProp Hexahedron3D8Connectivities
                } else {
                    # UPwUpdatedLagrangianElement2D3N
                    WriteElements FileVar [lindex $Groups $i] triangle UPwUpdatedLagrangianElement2D3N $BodyElemsProp Triangle2D3Connectivities
                    # UPwUpdatedLagrangianElement2D4N
                    WriteElements FileVar [lindex $Groups $i] quadrilateral UPwUpdatedLagrangianElement2D4N $BodyElemsProp Quadrilateral2D4Connectivities
                    # UPwUpdatedLagrangianElement3D4N
                    WriteElements FileVar [lindex $Groups $i] tetrahedra UPwUpdatedLagrangianElement3D4N $BodyElemsProp Quadrilateral2D4Connectivities
                    # UPwUpdatedLagrangianElement3D8N
                    WriteElements FileVar [lindex $Groups $i] hexahedra UPwUpdatedLagrangianElement3D8N $BodyElemsProp Hexahedron3D8Connectivities
                }
            }
        } else {
            for {set i 0} {$i < [llength $Groups]} {incr i} {
                # Elements Property
                set BodyElemsProp [dict get $PropertyDict [lindex [lindex $Groups $i] 1]]
                if {$IsMoveMesh eq false} {
                    # UPwSmallStrainFICElement2D3N
                    WriteElements FileVar [lindex $Groups $i] triangle UPwSmallStrainFICElement2D3N $BodyElemsProp Triangle2D3Connectivities
                    # UPwSmallStrainFICElement2D4N
                    WriteElements FileVar [lindex $Groups $i] quadrilateral UPwSmallStrainFICElement2D4N $BodyElemsProp Quadrilateral2D4Connectivities
                    # UPwSmallStrainFICElement3D4N
                    WriteElements FileVar [lindex $Groups $i] tetrahedra UPwSmallStrainFICElement3D4N $BodyElemsProp Quadrilateral2D4Connectivities
                    # UPwSmallStrainFICElement3D8N
                    WriteElements FileVar [lindex $Groups $i] hexahedra UPwSmallStrainFICElement3D8N $BodyElemsProp Hexahedron3D8Connectivities
                } else {
                    # UPwUpdatedLagrangianFICElement2D3N
                    WriteElements FileVar [lindex $Groups $i] triangle UPwUpdatedLagrangianFICElement2D3N $BodyElemsProp Triangle2D3Connectivities
                    # UPwUpdatedLagrangianFICElement2D4N
                    WriteElements FileVar [lindex $Groups $i] quadrilateral UPwUpdatedLagrangianFICElement2D4N $BodyElemsProp Quadrilateral2D4Connectivities
                    # UPwUpdatedLagrangianFICElement3D4N
                    WriteElements FileVar [lindex $Groups $i] tetrahedra UPwUpdatedLagrangianFICElement3D4N $BodyElemsProp Quadrilateral2D4Connectivities
                    # UPwUpdatedLagrangianFICElement3D8N
                    WriteElements FileVar [lindex $Groups $i] hexahedra UPwUpdatedLagrangianFICElement3D8N $BodyElemsProp Hexahedron3D8Connectivities
                }
            }
        }
    } elseif {$IsQuadratic eq 1} {
        for {set i 0} {$i < [llength $Groups]} {incr i} {
            # Elements Property
            set BodyElemsProp [dict get $PropertyDict [lindex [lindex $Groups $i] 1]]
            if {$IsMoveMesh eq false} {
                # SmallStrainUPwDiffOrderElement2D6N
                WriteElements FileVar [lindex $Groups $i] triangle SmallStrainUPwDiffOrderElement2D6N $BodyElemsProp Triangle2D6Connectivities
                # SmallStrainUPwDiffOrderElement2D8N
                WriteElements FileVar [lindex $Groups $i] quadrilateral SmallStrainUPwDiffOrderElement2D8N $BodyElemsProp Hexahedron3D8Connectivities
                # SmallStrainUPwDiffOrderElement3D10N
                WriteElements FileVar [lindex $Groups $i] tetrahedra SmallStrainUPwDiffOrderElement3D10N $BodyElemsProp Tetrahedron3D10Connectivities
                # SmallStrainUPwDiffOrderElement3D20N
                WriteElements FileVar [lindex $Groups $i] hexahedra SmallStrainUPwDiffOrderElement3D20N $BodyElemsProp Hexahedron3D20Connectivities
            } else {
                # UpdatedLagrangianUPwDiffOrderElement2D6N
                WriteElements FileVar [lindex $Groups $i] triangle UpdatedLagrangianUPwDiffOrderElement2D6N $BodyElemsProp Triangle2D6Connectivities
                # UpdatedLagrangianUPwDiffOrderElement2D8N
                WriteElements FileVar [lindex $Groups $i] quadrilateral UpdatedLagrangianUPwDiffOrderElement2D8N $BodyElemsProp Hexahedron3D8Connectivities
                # UpdatedLagrangianUPwDiffOrderElement3D10N
                WriteElements FileVar [lindex $Groups $i] tetrahedra UpdatedLagrangianUPwDiffOrderElement3D10N $BodyElemsProp Tetrahedron3D10Connectivities
                # UpdatedLagrangianUPwDiffOrderElement3D20N
                WriteElements FileVar [lindex $Groups $i] hexahedra UpdatedLagrangianUPwDiffOrderElement3D20N $BodyElemsProp Hexahedron3D20Connectivities
            }
        }
    } else {
        for {set i 0} {$i < [llength $Groups]} {incr i} {
            # Elements Property
            set BodyElemsProp [dict get $PropertyDict [lindex [lindex $Groups $i] 1]]
            if {$IsMoveMesh eq false} {
                # SmallStrainUPwDiffOrderElement2D6N
                WriteElements FileVar [lindex $Groups $i] triangle SmallStrainUPwDiffOrderElement2D6N $BodyElemsProp Triangle2D6Connectivities
                # SmallStrainUPwDiffOrderElement2D9N
                WriteElements FileVar [lindex $Groups $i] quadrilateral SmallStrainUPwDiffOrderElement2D9N $BodyElemsProp Quadrilateral2D9Connectivities
                # SmallStrainUPwDiffOrderElement3D10N
                WriteElements FileVar [lindex $Groups $i] tetrahedra SmallStrainUPwDiffOrderElement3D10N $BodyElemsProp Tetrahedron3D10Connectivities
                # SmallStrainUPwDiffOrderElement3D27N
                WriteElements FileVar [lindex $Groups $i] hexahedra SmallStrainUPwDiffOrderElement3D27N $BodyElemsProp Hexahedron3D27Connectivities
            } else {
                # UpdatedLagrangianUPwDiffOrderElement2D6N
                WriteElements FileVar [lindex $Groups $i] triangle UpdatedLagrangianUPwDiffOrderElement2D6N $BodyElemsProp Triangle2D6Connectivities
                # UpdatedLagrangianUPwDiffOrderElement2D9N
                WriteElements FileVar [lindex $Groups $i] quadrilateral UpdatedLagrangianUPwDiffOrderElement2D9N $BodyElemsProp Quadrilateral2D9Connectivities
                # UpdatedLagrangianUPwDiffOrderElement3D10N
                WriteElements FileVar [lindex $Groups $i] tetrahedra UpdatedLagrangianUPwDiffOrderElement3D10N $BodyElemsProp Tetrahedron3D10Connectivities
                # UpdatedLagrangianUPwDiffOrderElement3D27N
                WriteElements FileVar [lindex $Groups $i] hexahedra UpdatedLagrangianUPwDiffOrderElement3D27N $BodyElemsProp Hexahedron3D27Connectivities
            }
        }
    }

    # Soil_drained part
    set Groups [GiD_Info conditions Soil_drained groups]
    if {$IsQuadratic eq 0} {
        if {$FIC eq false} {
            for {set i 0} {$i < [llength $Groups]} {incr i} {
                # Elements Property
                set BodyElemsProp [dict get $PropertyDict [lindex [lindex $Groups $i] 1]]
                if {$IsMoveMesh eq false} {
                    # UPwSmallStrainElement2D3N
                    WriteElements FileVar [lindex $Groups $i] triangle UPwSmallStrainElement2D3N $BodyElemsProp Triangle2D3Connectivities
                    # UPwSmallStrainElement2D4N
                    WriteElements FileVar [lindex $Groups $i] quadrilateral UPwSmallStrainElement2D4N $BodyElemsProp Quadrilateral2D4Connectivities
                    # UPwSmallStrainElement3D4N
                    WriteElements FileVar [lindex $Groups $i] tetrahedra UPwSmallStrainElement3D4N $BodyElemsProp Quadrilateral2D4Connectivities
                    # UPwSmallStrainElement3D8N
                    WriteElements FileVar [lindex $Groups $i] hexahedra UPwSmallStrainElement3D8N $BodyElemsProp Hexahedron3D8Connectivities
                } else {
                    # UPwUpdatedLagrangianElement2D3N
                    WriteElements FileVar [lindex $Groups $i] triangle UPwUpdatedLagrangianElement2D3N $BodyElemsProp Triangle2D3Connectivities
                    # UPwUpdatedLagrangianElement2D4N
                    WriteElements FileVar [lindex $Groups $i] quadrilateral UPwUpdatedLagrangianElement2D4N $BodyElemsProp Quadrilateral2D4Connectivities
                    # UPwUpdatedLagrangianElement3D4N
                    WriteElements FileVar [lindex $Groups $i] tetrahedra UPwUpdatedLagrangianElement3D4N $BodyElemsProp Quadrilateral2D4Connectivities
                    # UPwUpdatedLagrangianElement3D8N
                    WriteElements FileVar [lindex $Groups $i] hexahedra UPwUpdatedLagrangianElement3D8N $BodyElemsProp Hexahedron3D8Connectivities
                }
            }
        } else {
            for {set i 0} {$i < [llength $Groups]} {incr i} {
                # Elements Property
                set BodyElemsProp [dict get $PropertyDict [lindex [lindex $Groups $i] 1]]
                if {$IsMoveMesh eq false} {
                    # UPwSmallStrainFICElement2D3N
                    WriteElements FileVar [lindex $Groups $i] triangle UPwSmallStrainFICElement2D3N $BodyElemsProp Triangle2D3Connectivities
                    # UPwSmallStrainFICElement2D4N
                    WriteElements FileVar [lindex $Groups $i] quadrilateral UPwSmallStrainFICElement2D4N $BodyElemsProp Quadrilateral2D4Connectivities
                    # UPwSmallStrainFICElement3D4N
                    WriteElements FileVar [lindex $Groups $i] tetrahedra UPwSmallStrainFICElement3D4N $BodyElemsProp Quadrilateral2D4Connectivities
                    # UPwSmallStrainFICElement3D8N
                    WriteElements FileVar [lindex $Groups $i] hexahedra UPwSmallStrainFICElement3D8N $BodyElemsProp Hexahedron3D8Connectivities
                } else {
                    # UPwUpdatedLagrangianFICElement2D3N
                    WriteElements FileVar [lindex $Groups $i] triangle UPwUpdatedLagrangianFICElement2D3N $BodyElemsProp Triangle2D3Connectivities
                    # UPwUpdatedLagrangianFICElement2D4N
                    WriteElements FileVar [lindex $Groups $i] quadrilateral UPwUpdatedLagrangianFICElement2D4N $BodyElemsProp Quadrilateral2D4Connectivities
                    # UPwUpdatedLagrangianFICElement3D4N
                    WriteElements FileVar [lindex $Groups $i] tetrahedra UPwUpdatedLagrangianFICElement3D4N $BodyElemsProp Quadrilateral2D4Connectivities
                    # UPwUpdatedLagrangianFICElement3D8N
                    WriteElements FileVar [lindex $Groups $i] hexahedra UPwUpdatedLagrangianFICElement3D8N $BodyElemsProp Hexahedron3D8Connectivities
                }
            }
        }
    } elseif {$IsQuadratic eq 1} {
        for {set i 0} {$i < [llength $Groups]} {incr i} {
            # Elements Property
            set BodyElemsProp [dict get $PropertyDict [lindex [lindex $Groups $i] 1]]
            if {$IsMoveMesh eq false} {
                # SmallStrainUPwDiffOrderElement2D6N
                WriteElements FileVar [lindex $Groups $i] triangle SmallStrainUPwDiffOrderElement2D6N $BodyElemsProp Triangle2D6Connectivities
                # SmallStrainUPwDiffOrderElement2D8N
                WriteElements FileVar [lindex $Groups $i] quadrilateral SmallStrainUPwDiffOrderElement2D8N $BodyElemsProp Hexahedron3D8Connectivities
                # SmallStrainUPwDiffOrderElement3D10N
                WriteElements FileVar [lindex $Groups $i] tetrahedra SmallStrainUPwDiffOrderElement3D10N $BodyElemsProp Tetrahedron3D10Connectivities
                # SmallStrainUPwDiffOrderElement3D20N
                WriteElements FileVar [lindex $Groups $i] hexahedra SmallStrainUPwDiffOrderElement3D20N $BodyElemsProp Hexahedron3D20Connectivities
            } else {
                # UpdatedLagrangianUPwDiffOrderElement2D6N
                WriteElements FileVar [lindex $Groups $i] triangle UpdatedLagrangianUPwDiffOrderElement2D6N $BodyElemsProp Triangle2D6Connectivities
                # UpdatedLagrangianUPwDiffOrderElement2D8N
                WriteElements FileVar [lindex $Groups $i] quadrilateral UpdatedLagrangianUPwDiffOrderElement2D8N $BodyElemsProp Hexahedron3D8Connectivities
                # UpdatedLagrangianUPwDiffOrderElement3D10N
                WriteElements FileVar [lindex $Groups $i] tetrahedra UpdatedLagrangianUPwDiffOrderElement3D10N $BodyElemsProp Tetrahedron3D10Connectivities
                # UpdatedLagrangianUPwDiffOrderElement3D20N
                WriteElements FileVar [lindex $Groups $i] hexahedra UpdatedLagrangianUPwDiffOrderElement3D20N $BodyElemsProp Hexahedron3D20Connectivities
            }
        }
    } else {
        for {set i 0} {$i < [llength $Groups]} {incr i} {
            # Elements Property
            set BodyElemsProp [dict get $PropertyDict [lindex [lindex $Groups $i] 1]]
            if {$IsMoveMesh eq false} {
                # SmallStrainUPwDiffOrderElement2D6N
                WriteElements FileVar [lindex $Groups $i] triangle SmallStrainUPwDiffOrderElement2D6N $BodyElemsProp Triangle2D6Connectivities
                # SmallStrainUPwDiffOrderElement2D9N
                WriteElements FileVar [lindex $Groups $i] quadrilateral SmallStrainUPwDiffOrderElement2D9N $BodyElemsProp Quadrilateral2D9Connectivities
                # SmallStrainUPwDiffOrderElement3D10N
                WriteElements FileVar [lindex $Groups $i] tetrahedra SmallStrainUPwDiffOrderElement3D10N $BodyElemsProp Tetrahedron3D10Connectivities
                # SmallStrainUPwDiffOrderElement3D27N
                WriteElements FileVar [lindex $Groups $i] hexahedra SmallStrainUPwDiffOrderElement3D27N $BodyElemsProp Hexahedron3D27Connectivities
            } else {
                # UpdatedLagrangianUPwDiffOrderElement2D6N
                WriteElements FileVar [lindex $Groups $i] triangle UpdatedLagrangianUPwDiffOrderElement2D6N $BodyElemsProp Triangle2D6Connectivities
                # UpdatedLagrangianUPwDiffOrderElement2D9N
                WriteElements FileVar [lindex $Groups $i] quadrilateral UpdatedLagrangianUPwDiffOrderElement2D9N $BodyElemsProp Quadrilateral2D9Connectivities
                # UpdatedLagrangianUPwDiffOrderElement3D10N
                WriteElements FileVar [lindex $Groups $i] tetrahedra UpdatedLagrangianUPwDiffOrderElement3D10N $BodyElemsProp Tetrahedron3D10Connectivities
                # UpdatedLagrangianUPwDiffOrderElement3D27N
                WriteElements FileVar [lindex $Groups $i] hexahedra UpdatedLagrangianUPwDiffOrderElement3D27N $BodyElemsProp Hexahedron3D27Connectivities
            }
        }
    }

    # Soil_undrained part
    set Groups [GiD_Info conditions Soil_undrained groups]
    if {$IsQuadratic eq 0} {
        if {$FIC eq false} {
            for {set i 0} {$i < [llength $Groups]} {incr i} {
                # Elements Property
                set BodyElemsProp [dict get $PropertyDict [lindex [lindex $Groups $i] 1]]
                if {$IsMoveMesh eq false} {
                    # UPwSmallStrainElement2D3N
                    WriteElements FileVar [lindex $Groups $i] triangle UPwSmallStrainElement2D3N $BodyElemsProp Triangle2D3Connectivities
                    # UPwSmallStrainElement2D4N
                    WriteElements FileVar [lindex $Groups $i] quadrilateral UPwSmallStrainElement2D4N $BodyElemsProp Quadrilateral2D4Connectivities
                    # UPwSmallStrainElement3D4N
                    WriteElements FileVar [lindex $Groups $i] tetrahedra UPwSmallStrainElement3D4N $BodyElemsProp Quadrilateral2D4Connectivities
                    # UPwSmallStrainElement3D8N
                    WriteElements FileVar [lindex $Groups $i] hexahedra UPwSmallStrainElement3D8N $BodyElemsProp Hexahedron3D8Connectivities
                } else {
                    # UPwUpdatedLagrangianElement2D3N
                    WriteElements FileVar [lindex $Groups $i] triangle UPwUpdatedLagrangianElement2D3N $BodyElemsProp Triangle2D3Connectivities
                    # UPwUpdatedLagrangianElement2D4N
                    WriteElements FileVar [lindex $Groups $i] quadrilateral UPwUpdatedLagrangianElement2D4N $BodyElemsProp Quadrilateral2D4Connectivities
                    # UPwUpdatedLagrangianElement3D4N
                    WriteElements FileVar [lindex $Groups $i] tetrahedra UPwUpdatedLagrangianElement3D4N $BodyElemsProp Quadrilateral2D4Connectivities
                    # UPwUpdatedLagrangianElement3D8N
                    WriteElements FileVar [lindex $Groups $i] hexahedra UPwUpdatedLagrangianElement3D8N $BodyElemsProp Hexahedron3D8Connectivities
                }
            }
        } else {
            for {set i 0} {$i < [llength $Groups]} {incr i} {
                # Elements Property
                set BodyElemsProp [dict get $PropertyDict [lindex [lindex $Groups $i] 1]]
                if {$IsMoveMesh eq false} {
                    # UPwSmallStrainFICElement2D3N
                    WriteElements FileVar [lindex $Groups $i] triangle UPwSmallStrainFICElement2D3N $BodyElemsProp Triangle2D3Connectivities
                    # UPwSmallStrainFICElement2D4N
                    WriteElements FileVar [lindex $Groups $i] quadrilateral UPwSmallStrainFICElement2D4N $BodyElemsProp Quadrilateral2D4Connectivities
                    # UPwSmallStrainFICElement3D4N
                    WriteElements FileVar [lindex $Groups $i] tetrahedra UPwSmallStrainFICElement3D4N $BodyElemsProp Quadrilateral2D4Connectivities
                    # UPwSmallStrainFICElement3D8N
                    WriteElements FileVar [lindex $Groups $i] hexahedra UPwSmallStrainFICElement3D8N $BodyElemsProp Hexahedron3D8Connectivities
                } else {
                    # UPwUpdatedLagrangianFICElement2D3N
                    WriteElements FileVar [lindex $Groups $i] triangle UPwUpdatedLagrangianFICElement2D3N $BodyElemsProp Triangle2D3Connectivities
                    # UPwUpdatedLagrangianFICElement2D4N
                    WriteElements FileVar [lindex $Groups $i] quadrilateral UPwUpdatedLagrangianFICElement2D4N $BodyElemsProp Quadrilateral2D4Connectivities
                    # UPwUpdatedLagrangianFICElement3D4N
                    WriteElements FileVar [lindex $Groups $i] tetrahedra UPwUpdatedLagrangianFICElement3D4N $BodyElemsProp Quadrilateral2D4Connectivities
                    # UPwUpdatedLagrangianFICElement3D8N
                    WriteElements FileVar [lindex $Groups $i] hexahedra UPwUpdatedLagrangianFICElement3D8N $BodyElemsProp Hexahedron3D8Connectivities
                }
            }
        }
    } elseif {$IsQuadratic eq 1} {
        for {set i 0} {$i < [llength $Groups]} {incr i} {
            # Elements Property
            set BodyElemsProp [dict get $PropertyDict [lindex [lindex $Groups $i] 1]]
            if {$IsMoveMesh eq false} {
                # SmallStrainUPwDiffOrderElement2D6N
                WriteElements FileVar [lindex $Groups $i] triangle SmallStrainUPwDiffOrderElement2D6N $BodyElemsProp Triangle2D6Connectivities
                # SmallStrainUPwDiffOrderElement2D8N
                WriteElements FileVar [lindex $Groups $i] quadrilateral SmallStrainUPwDiffOrderElement2D8N $BodyElemsProp Hexahedron3D8Connectivities
                # SmallStrainUPwDiffOrderElement3D10N
                WriteElements FileVar [lindex $Groups $i] tetrahedra SmallStrainUPwDiffOrderElement3D10N $BodyElemsProp Tetrahedron3D10Connectivities
                # SmallStrainUPwDiffOrderElement3D20N
                WriteElements FileVar [lindex $Groups $i] hexahedra SmallStrainUPwDiffOrderElement3D20N $BodyElemsProp Hexahedron3D20Connectivities
            } else {
                # UpdatedLagrangianUPwDiffOrderElement2D6N
                WriteElements FileVar [lindex $Groups $i] triangle UpdatedLagrangianUPwDiffOrderElement2D6N $BodyElemsProp Triangle2D6Connectivities
                # UpdatedLagrangianUPwDiffOrderElement2D8N
                WriteElements FileVar [lindex $Groups $i] quadrilateral UpdatedLagrangianUPwDiffOrderElement2D8N $BodyElemsProp Hexahedron3D8Connectivities
                # UpdatedLagrangianUPwDiffOrderElement3D10N
                WriteElements FileVar [lindex $Groups $i] tetrahedra UpdatedLagrangianUPwDiffOrderElement3D10N $BodyElemsProp Tetrahedron3D10Connectivities
                # UpdatedLagrangianUPwDiffOrderElement3D20N
                WriteElements FileVar [lindex $Groups $i] hexahedra UpdatedLagrangianUPwDiffOrderElement3D20N $BodyElemsProp Hexahedron3D20Connectivities
            }
        }
    } else {
        for {set i 0} {$i < [llength $Groups]} {incr i} {
            # Elements Property
            set BodyElemsProp [dict get $PropertyDict [lindex [lindex $Groups $i] 1]]
            if {$IsMoveMesh eq false} {
                # SmallStrainUPwDiffOrderElement2D6N
                WriteElements FileVar [lindex $Groups $i] triangle SmallStrainUPwDiffOrderElement2D6N $BodyElemsProp Triangle2D6Connectivities
                # SmallStrainUPwDiffOrderElement2D9N
                WriteElements FileVar [lindex $Groups $i] quadrilateral SmallStrainUPwDiffOrderElement2D9N $BodyElemsProp Quadrilateral2D9Connectivities
                # SmallStrainUPwDiffOrderElement3D10N
                WriteElements FileVar [lindex $Groups $i] tetrahedra SmallStrainUPwDiffOrderElement3D10N $BodyElemsProp Tetrahedron3D10Connectivities
                # SmallStrainUPwDiffOrderElement3D27N
                WriteElements FileVar [lindex $Groups $i] hexahedra SmallStrainUPwDiffOrderElement3D27N $BodyElemsProp Hexahedron3D27Connectivities
            } else {
                # UpdatedLagrangianUPwDiffOrderElement2D6N
                WriteElements FileVar [lindex $Groups $i] triangle UpdatedLagrangianUPwDiffOrderElement2D6N $BodyElemsProp Triangle2D6Connectivities
                # UpdatedLagrangianUPwDiffOrderElement2D9N
                WriteElements FileVar [lindex $Groups $i] quadrilateral UpdatedLagrangianUPwDiffOrderElement2D9N $BodyElemsProp Quadrilateral2D9Connectivities
                # UpdatedLagrangianUPwDiffOrderElement3D10N
                WriteElements FileVar [lindex $Groups $i] tetrahedra UpdatedLagrangianUPwDiffOrderElement3D10N $BodyElemsProp Tetrahedron3D10Connectivities
                # UpdatedLagrangianUPwDiffOrderElement3D27N
                WriteElements FileVar [lindex $Groups $i] hexahedra UpdatedLagrangianUPwDiffOrderElement3D27N $BodyElemsProp Hexahedron3D27Connectivities
            }
        }
    }

    # Non_porous part
    set Groups [GiD_Info conditions Non_porous groups]
    if {$IsQuadratic eq 0} {
        for {set i 0} {$i < [llength $Groups]} {incr i} {
            # Elements Property
            set BodyElemsProp [dict get $PropertyDict [lindex [lindex $Groups $i] 1]]
            if {$IsMoveMesh eq false} {
                # SmallDisplacement2D3N
                WriteElements FileVar [lindex $Groups $i] triangle SmallDisplacementElement2D3N $BodyElemsProp Triangle2D3Connectivities
                # SmallDisplacement2D4N
                WriteElements FileVar [lindex $Groups $i] quadrilateral SmallDisplacementElement2D4N $BodyElemsProp Quadrilateral2D4Connectivities
                # SmallDisplacement3D4N
                WriteElements FileVar [lindex $Groups $i] tetrahedra SmallDisplacementElement3D4N $BodyElemsProp Quadrilateral2D4Connectivities
                # SmallDisplacement3D8N
                WriteElements FileVar [lindex $Groups $i] hexahedra SmallDisplacementElement3D8N $BodyElemsProp Hexahedron3D8Connectivities
            } else {
                # UpdatedLagrangianElement2D3N
                WriteElements FileVar [lindex $Groups $i] triangle UpdatedLagrangianElement2D3N $BodyElemsProp Triangle2D3Connectivities
                # UpdatedLagrangianElement2D4N
                WriteElements FileVar [lindex $Groups $i] quadrilateral UpdatedLagrangianElement2D4N $BodyElemsProp Quadrilateral2D4Connectivities
                # UpdatedLagrangianElement3D4N
                WriteElements FileVar [lindex $Groups $i] tetrahedra UpdatedLagrangianElement3D4N $BodyElemsProp Quadrilateral2D4Connectivities
                # UpdatedLagrangianElement3D8N
                WriteElements FileVar [lindex $Groups $i] hexahedra UpdatedLagrangianElement3D8N $BodyElemsProp Hexahedron3D8Connectivities
            }
        }
    } elseif {$IsQuadratic eq 1} {
        for {set i 0} {$i < [llength $Groups]} {incr i} {
            # Elements Property
            set BodyElemsProp [dict get $PropertyDict [lindex [lindex $Groups $i] 1]]
            if {$IsMoveMesh eq false} {
                # SmallDisplacement2D6N
                WriteElements FileVar [lindex $Groups $i] triangle SmallDisplacementElement2D6N $BodyElemsProp Triangle2D6Connectivities
                # SmallDisplacement2D8N
                WriteElements FileVar [lindex $Groups $i] quadrilateral SmallDisplacementElement2D8N $BodyElemsProp Hexahedron3D8Connectivities
                # SmallDisplacement3D10N
                WriteElements FileVar [lindex $Groups $i] tetrahedra SmallDisplacementElement3D10N $BodyElemsProp Tetrahedron3D10Connectivities
                # SmallDisplacement3D20N
                WriteElements FileVar [lindex $Groups $i] hexahedra SmallDisplacementElement3D20N $BodyElemsProp Hexahedron3D20Connectivities
            } else {
                # UpdatedLagrangianElement2D6N
                WriteElements FileVar [lindex $Groups $i] triangle UpdatedLagrangianElement2D6N $BodyElemsProp Triangle2D6Connectivities
                # UpdatedLagrangianElement2D8N
                WriteElements FileVar [lindex $Groups $i] quadrilateral UpdatedLagrangianElement2D8N $BodyElemsProp Hexahedron3D8Connectivities
                # UpdatedLagrangianElement3D10N
                WriteElements FileVar [lindex $Groups $i] tetrahedra UpdatedLagrangianElement3D10N $BodyElemsProp Tetrahedron3D10Connectivities
                # UpdatedLagrangianElement3D20N
                WriteElements FileVar [lindex $Groups $i] hexahedra UpdatedLagrangianElement3D20N $BodyElemsProp Hexahedron3D20Connectivities
            }
        }
    }

    # Soil_Groundwater_Flow part
    set Groups [GiD_Info conditions Soil_Groundwater_Flow groups]
    if {$SolutionType eq "Steady_State_Groundwater_Flow"} {
        if {$IsQuadratic eq 0} {
            for {set i 0} {$i < [llength $Groups]} {incr i} {
                # Elements Property
                set BodyElemsProp [dict get $PropertyDict [lindex [lindex $Groups $i] 1]]

                # SteadyStatePwElement2D3N
                WriteElements FileVar [lindex $Groups $i] triangle SteadyStatePwElement2D3N $BodyElemsProp Triangle2D3Connectivities
                # SteadyStatePwElement2D4N
                WriteElements FileVar [lindex $Groups $i] quadrilateral SteadyStatePwElement2D4N $BodyElemsProp Quadrilateral2D4Connectivities
                # SteadyStatePwElement3D4N
                WriteElements FileVar [lindex $Groups $i] tetrahedra SteadyStatePwElement3D4N $BodyElemsProp Quadrilateral2D4Connectivities
                # SteadyStatePwElement3D8N
                WriteElements FileVar [lindex $Groups $i] hexahedra SteadyStatePwElement3D8N $BodyElemsProp Hexahedron3D8Connectivities
            }
        } else {
            for {set i 0} {$i < [llength $Groups]} {incr i} {
                # Elements Property
                set BodyElemsProp [dict get $PropertyDict [lindex [lindex $Groups $i] 1]]
                # SteadyStatePwElement2D6N
                WriteElements FileVar [lindex $Groups $i] triangle SteadyStatePwElement2D6N $BodyElemsProp Triangle2D6Connectivities
                # SteadyStatePwElement2D8N
                WriteElements FileVar [lindex $Groups $i] quadrilateral SteadyStatePwElement2D8N $BodyElemsProp Hexahedron3D8Connectivities
                # SteadyStatePwElement2D9N
                WriteElements FileVar [lindex $Groups $i] quadrilateral SteadyStatePwElement2D9N $BodyElemsProp Quadrilateral2D9Connectivities

                # SteadyStatePwElement3D10N
                WriteElements FileVar [lindex $Groups $i] tetrahedra SteadyStatePwElement3D10N $BodyElemsProp Tetrahedron3D10Connectivities
                # SteadyStatePwElement3D20N
                WriteElements FileVar [lindex $Groups $i] hexahedra SteadyStatePwElement3D20N $BodyElemsProp Hexahedron3D20Connectivities
                # SteadyStatePwElement3D27N
                WriteElements FileVar [lindex $Groups $i] hexahedra SteadyStatePwElement3D27N $BodyElemsProp Hexahedron3D27Connectivities
            }
        }
    } else {
        if {$IsQuadratic eq 0} {
            for {set i 0} {$i < [llength $Groups]} {incr i} {
                # Elements Property
                set BodyElemsProp [dict get $PropertyDict [lindex [lindex $Groups $i] 1]]

                # TransientPwElement2D3N
                WriteElements FileVar [lindex $Groups $i] triangle TransientPwElement2D3N $BodyElemsProp Triangle2D3Connectivities
                # TransientPwElement2D4N
                WriteElements FileVar [lindex $Groups $i] quadrilateral TransientPwElement2D4N $BodyElemsProp Quadrilateral2D4Connectivities
                # TransientPwElement3D4N
                WriteElements FileVar [lindex $Groups $i] tetrahedra TransientPwElement3D4N $BodyElemsProp Quadrilateral2D4Connectivities
                # TransientPwElement3D8N
                WriteElements FileVar [lindex $Groups $i] hexahedra TransientPwElement3D8N $BodyElemsProp Hexahedron3D8Connectivities
            }
        } else {
            for {set i 0} {$i < [llength $Groups]} {incr i} {
                # Elements Property
                set BodyElemsProp [dict get $PropertyDict [lindex [lindex $Groups $i] 1]]
                # TransientPwElement2D6N
                WriteElements FileVar [lindex $Groups $i] triangle TransientPwElement2D6N $BodyElemsProp Triangle2D6Connectivities
                # TransientPwElement2D8N
                WriteElements FileVar [lindex $Groups $i] quadrilateral TransientPwElement2D8N $BodyElemsProp Hexahedron3D8Connectivities
                # TransientPwElement2D9N
                WriteElements FileVar [lindex $Groups $i] quadrilateral TransientPwElement2D9N $BodyElemsProp Quadrilateral2D9Connectivities

                # TransientPwElement3D10N
                WriteElements FileVar [lindex $Groups $i] tetrahedra TransientPwElement3D10N $BodyElemsProp Tetrahedron3D10Connectivities
                # TransientPwElement3D20N
                WriteElements FileVar [lindex $Groups $i] hexahedra TransientPwElement3D20N $BodyElemsProp Hexahedron3D20Connectivities
                # TransientPwElement3D27N
                WriteElements FileVar [lindex $Groups $i] hexahedra TransientPwElement3D27N $BodyElemsProp Hexahedron3D27Connectivities
            }
        }
    }

    # Beam_Part
    set Groups [GiD_Info conditions Beam groups]
    if {$IsQuadratic eq 0} {
        for {set i 0} {$i < [llength $Groups]} {incr i} {
            # Elements Property
            set BodyElemsProp [dict get $PropertyDict [lindex [lindex $Groups $i] 1]]

            if {$Dim eq 2} {
              # CrBeamElement2D2N
              WriteElements FileVar [lindex $Groups $i] line GeoCrBeamElementLinear2D2N $BodyElemsProp Line2D2Connectivities
            } else {
              # CrBeamElement3D2N
              WriteElements FileVar [lindex $Groups $i] line GeoCrBeamElementLinear3D2N $BodyElemsProp Line2D2Connectivities
            }
        }
    } elseif {$IsQuadratic eq 1} {
        for {set i 0} {$i < [llength $Groups]} {incr i} {
            # Elements Property
            set BodyElemsProp [dict get $PropertyDict [lindex [lindex $Groups $i] 1]]

            if {$Dim eq 2} {
              # CrBeamElement2D2N twice
              WriteElementsTwoParts FileVar [lindex $Groups $i] line GeoCrBeamElementLinear2D2N $BodyElemsProp Line2D2ConnectivitiesPart1 Line2D2ConnectivitiesPart2

            } else {
              # GeoCrBeamElementLinear3D2N twice
              WriteElementsTwoParts FileVar [lindex $Groups $i] line GeoCrBeamElementLinear3D2N $BodyElemsProp Line2D2ConnectivitiesPart1 Line2D2ConnectivitiesPart2
            }
        }
    }


    # Shell_thin_corotational_Part
    set Groups [GiD_Info conditions Shell_thin_corotational groups]
    if {$IsQuadratic eq 0} {
      for {set i 0} {$i < [llength $Groups]} {incr i} {
          # Elements Property
          set BodyElemsProp [dict get $PropertyDict [lindex [lindex $Groups $i] 1]]

          # ShellThinElementCorotational3D3N
          WriteElements FileVar [lindex $Groups $i] triangle ShellThinElementCorotational3D3N $BodyElemsProp Triangle3D3Connectivities

          # ShellThinElementCorotational3D4N
          WriteElements FileVar [lindex $Groups $i] quadrilateral ShellThinElementCorotational3D4N $BodyElemsProp Quadrilateral3D4Connectivities
      }
    } elseif {$IsQuadratic eq 1} {
      for {set i 0} {$i < [llength $Groups]} {incr i} {
          # Elements Property
          set BodyElemsProp [dict get $PropertyDict [lindex [lindex $Groups $i] 1]]

          # ShellThinElementCorotational3D3N four times
          WriteElementsFourParts FileVar [lindex $Groups $i] triangle ShellThinElementCorotational3D3N $BodyElemsProp \
                                 Triangle3D3ConnectivitiesPart1\
                                 Triangle3D3ConnectivitiesPart2\
                                 Triangle3D3ConnectivitiesPart3\
                                 Triangle3D3ConnectivitiesPart4

          # ShellThinElementCorotational3D4N
          WriteElements FileVar [lindex $Groups $i] quadrilateral ShellThinElementCorotational3D4N $BodyElemsProp Quadrilateral3D4Connectivities
      }
    }

    # Shell_thick_corotational_Part
    set Groups [GiD_Info conditions Shell_thick_corotational groups]
    if {$IsQuadratic eq 0} {
      for {set i 0} {$i < [llength $Groups]} {incr i} {
          # Elements Property
          set BodyElemsProp [dict get $PropertyDict [lindex [lindex $Groups $i] 1]]

          # ShellThickElementCorotational3D3N
          WriteElements FileVar [lindex $Groups $i] triangle ShellThickElementCorotational3D3N $BodyElemsProp Triangle3D3Connectivities

          # ShellThickElementCorotational3D4N
          WriteElements FileVar [lindex $Groups $i] quadrilateral ShellThickElementCorotational3D4N $BodyElemsProp Quadrilateral3D4Connectivities
      }
    } elseif {$IsQuadratic eq 1} {
      for {set i 0} {$i < [llength $Groups]} {incr i} {
          # Elements Property
          set BodyElemsProp [dict get $PropertyDict [lindex [lindex $Groups $i] 1]]

          # ShellThickElementCorotational3D3N four times
          WriteElementsFourParts FileVar [lindex $Groups $i] triangle ShellThickElementCorotational3D3N $BodyElemsProp \
                                 Triangle3D3ConnectivitiesPart1\
                                 Triangle3D3ConnectivitiesPart2\
                                 Triangle3D3ConnectivitiesPart3\
                                 Triangle3D3ConnectivitiesPart4

          # ShellThickElementCorotational3D4N
          WriteElements FileVar [lindex $Groups $i] quadrilateral ShellThickElementCorotational3D4N $BodyElemsProp Quadrilateral3D4Connectivities
      }
    }


    # Truss_Part
    set Groups [GiD_Info conditions Truss groups]
    if {$IsQuadratic eq 0} {
        for {set i 0} {$i < [llength $Groups]} {incr i} {
            # Elements Property
            set BodyElemsProp [dict get $PropertyDict [lindex [lindex $Groups $i] 1]]

            # TrussLinearElement3D2N
            WriteElements FileVar [lindex $Groups $i] line GeoTrussLinearElement3D2N $BodyElemsProp Line2D2Connectivities
        }
    } elseif {$IsQuadratic eq 1} {
        for {set i 0} {$i < [llength $Groups]} {incr i} {
            # Elements Property
            set BodyElemsProp [dict get $PropertyDict [lindex [lindex $Groups $i] 1]]

            # TrussLinearElement3D2N twice
            WriteElementsTwoParts FileVar [lindex $Groups $i] line GeoTrussLinearElement3D2N $BodyElemsProp Line2D2ConnectivitiesPart1 Line2D2ConnectivitiesPart2
        }
    }

    # Anchor_Part
    #set AnchorId 0
    set AnchorElementDict [dict create]
    set Groups [GiD_Info conditions Anchor groups]
    if {$IsQuadratic eq 0} {
        for {set i 0} {$i < [llength $Groups]} {incr i} {
            # Elements Property
            set BodyElemsProp [dict get $PropertyDict [lindex [lindex $Groups $i] 1]]
            
            # CableElement3D2N
            #dict set AnchorElementDict $AnchorId [WriteAnchorElements FileVar [lindex $Groups $i] line GeoCableElement3D2N $BodyElemsProp Line2D2Connectivities] 
            WriteElements FileVar [lindex $Groups $i] line GeoCableElement3D2N $BodyElemsProp Line2D2Connectivities
            #incr AnchorId
        }
    } elseif {$IsQuadratic eq 1} {
        for {set i 0} {$i < [llength $Groups]} {incr i} {
            # Elements Property
            set BodyElemsProp [dict get $PropertyDict [lindex [lindex $Groups $i] 1]]

            # CableElement3D3N
            #dict set AnchorElementDict $AnchorId [WriteAnchorElements FileVar [lindex $Groups $i] line GeoCableElement3D2N $BodyElemsProp Line2D3Connectivities]
            WriteElementsTwoParts FileVar [lindex $Groups $i] line GeoCableElement3D2N $BodyElemsProp Line2D2ConnectivitiesPart1 Line2D2ConnectivitiesPart2
            #incr AnchorId
        }
    }

    # Interface two-phase, drained and undrained parts
    set interface_Groups [list [GiD_Info conditions Interface_two_phase groups] [GiD_Info conditions Interface_drained groups] [GiD_Info conditions Interface_undrained groups]]
    if {$IsQuadratic eq 0} {
        foreach Groups $interface_Groups {
            for {set i 0} {$i < [llength $Groups]} {incr i} {
                if {[lindex [lindex $Groups $i] 3] eq false} {
                    # Elements Property
                    set InterfaceElemsProp [dict get $PropertyDict [lindex [lindex $Groups $i] 1]]
                    # UPwSmallStrainInterfaceElement2D4N
                    WriteElements FileVar [lindex $Groups $i] quadrilateral UPwSmallStrainInterfaceElement2D4N $InterfaceElemsProp Quadrilateral2D4Connectivities

                    # UPwSmallStrainInterfaceElement3D6N
                    WriteElements FileVar [lindex $Groups $i] prism UPwSmallStrainInterfaceElement3D6N $InterfaceElemsProp PrismInterface3D6Connectivities
                    # UPwSmallStrainInterfaceElement3D8N
                    WriteElements FileVar [lindex $Groups $i] hexahedra UPwSmallStrainInterfaceElement3D8N $InterfaceElemsProp HexahedronInterface3D8Connectivities
                } else {
                    # Elements Property
                    set LinkInterfaceElemsProp [dict get $PropertyDict [lindex [lindex $Groups $i] 1]]
                    # UPwSmallStrainLinkInterfaceElement2D4N
                    WriteElements FileVar [lindex $Groups $i] quadrilateral UPwSmallStrainLinkInterfaceElement2D4N $LinkInterfaceElemsProp QuadrilateralInterface2D4Connectivities
                    WriteElements FileVar [lindex $Groups $i] triangle UPwSmallStrainLinkInterfaceElement2D4N $LinkInterfaceElemsProp TriangleInterface2D4Connectivities
                    # UPwSmallStrainLinkInterfaceElement3D6N
                    WriteElements FileVar [lindex $Groups $i] prism UPwSmallStrainLinkInterfaceElement3D6N $LinkInterfaceElemsProp Triangle2D6Connectivities
                    WriteElements FileVar [lindex $Groups $i] tetrahedra UPwSmallStrainLinkInterfaceElement3D6N $LinkInterfaceElemsProp TetrahedronInterface3D6Connectivities
                    # UPwSmallStrainLinkInterfaceElement3D8N
                    WriteElements FileVar [lindex $Groups $i] hexahedra UPwSmallStrainLinkInterfaceElement3D8N $LinkInterfaceElemsProp Hexahedron3D8Connectivities
                }
            }
        }
    } elseif {$IsQuadratic eq 1} {
        foreach Groups $interface_Groups {
            for {set i 0} {$i < [llength $Groups]} {incr i} {
                if {[lindex [lindex $Groups $i] 3] eq false} {
                    # Elements Property
                    set InterfaceElemsProp [dict get $PropertyDict [lindex [lindex $Groups $i] 1]]
                    # UPwSmallStrainInterfaceElement2D4N twice
                    WriteElementsTwoParts FileVar [lindex $Groups $i] quadrilateral UPwSmallStrainInterfaceElement2D4N $InterfaceElemsProp Quadrilateral2D4ConnectivitiesPart1 Quadrilateral2D4ConnectivitiesPart2

                    # UPwSmallStrainInterfaceElement3D6N twice
                    WriteElementsFourParts FileVar [lindex $Groups $i] prism UPwSmallStrainInterfaceElement3D6N $InterfaceElemsProp \
                                                      PrismInterface3D6ConnectivitiesPart1 \
                                                      PrismInterface3D6ConnectivitiesPart2 \
                                                      PrismInterface3D6ConnectivitiesPart3 \
                                                      PrismInterface3D6ConnectivitiesPart4
                    # UPwSmallStrainInterfaceElement3D8N
                    WriteElements FileVar [lindex $Groups $i] hexahedra UPwSmallStrainInterfaceElement3D8N $InterfaceElemsProp HexahedronInterface3D8Connectivities
                } else {
                    # Elements Property
                    set LinkInterfaceElemsProp [dict get $PropertyDict [lindex [lindex $Groups $i] 1]]
                    # UPwSmallStrainLinkInterfaceElement2D4N twice
                    WriteElementsTwoParts FileVar [lindex $Groups $i] quadrilateral UPwSmallStrainLinkInterfaceElement2D4N $InterfaceElemsProp Quadrilateral2D4ConnectivitiesPart1 Quadrilateral2D4ConnectivitiesPart2
                    WriteElements FileVar [lindex $Groups $i] triangle UPwSmallStrainLinkInterfaceElement2D4N $LinkInterfaceElemsProp TriangleInterface2D4Connectivities
                    # UPwSmallStrainLinkInterfaceElement3D6N
                    # UPwSmallStrainInterfaceElement3D6N twice
                    WriteElementsFourParts FileVar [lindex $Groups $i] prism UPwSmallStrainLinkInterfaceElement3D6N $InterfaceElemsProp \
                                                      PrismInterface3D6ConnectivitiesPart1 \
                                                      PrismInterface3D6ConnectivitiesPart2 \
                                                      PrismInterface3D6ConnectivitiesPart3 \
                                                      PrismInterface3D6ConnectivitiesPart4
                    WriteElements FileVar [lindex $Groups $i] tetrahedra UPwSmallStrainLinkInterfaceElement3D6N $LinkInterfaceElemsProp TetrahedronInterface3D6Connectivities
                    # UPwSmallStrainLinkInterfaceElement3D8N
                    WriteElements FileVar [lindex $Groups $i] hexahedra UPwSmallStrainLinkInterfaceElement3D8N $LinkInterfaceElemsProp Hexahedron3D8Connectivities
                }
            }
        }
    }

    #Interface_Groundwater_flow
    set Groups [GiD_Info conditions Interface_Groundwater_flow groups]
    if {$SolutionType eq "Steady_State_Groundwater_Flow"} {
        if {$IsQuadratic eq 0} {
            for {set i 0} {$i < [llength $Groups]} {incr i} {
                if {[lindex [lindex $Groups $i] 3] eq false} {
                    # Elements Property
                    set InterfaceElemsProp [dict get $PropertyDict [lindex [lindex $Groups $i] 1]]
                    # SteadyStatePwInterfaceElement2D4N
                    WriteElements FileVar [lindex $Groups $i] quadrilateral SteadyStatePwInterfaceElement2D4N $InterfaceElemsProp Quadrilateral2D4Connectivities

                    # SteadyStatePwInterfaceElement3D6N
                    WriteElements FileVar [lindex $Groups $i] prism SteadyStatePwInterfaceElement3D6N $InterfaceElemsProp PrismInterface3D6Connectivities
                    # SteadyStatePwInterfaceElement3D8N
                    WriteElements FileVar [lindex $Groups $i] hexahedra SteadyStatePwInterfaceElement3D8N $InterfaceElemsProp HexahedronInterface3D8Connectivities
                } else {
                    # Elements Property
                    set LinkInterfaceElemsProp [dict get $PropertyDict [lindex [lindex $Groups $i] 1]]
                    # SteadyStatePwLinkInterfaceElement2D4N
                    WriteElements FileVar [lindex $Groups $i] quadrilateral SteadyStatePwLinkInterfaceElement2D4N $LinkInterfaceElemsProp QuadrilateralInterface2D4Connectivities
                    WriteElements FileVar [lindex $Groups $i] triangle SteadyStatePwLinkInterfaceElement2D4N $LinkInterfaceElemsProp TriangleInterface2D4Connectivities
                    # SteadyStatePwLinkInterfaceElement3D6N
                    WriteElements FileVar [lindex $Groups $i] prism SteadyStatePwLinkInterfaceElement3D6N $LinkInterfaceElemsProp Triangle2D6Connectivities
                    WriteElements FileVar [lindex $Groups $i] tetrahedra SteadyStatePwLinkInterfaceElement3D6N $LinkInterfaceElemsProp TetrahedronInterface3D6Connectivities
                    # SteadyStatePwLinkInterfaceElement3D8N
                    WriteElements FileVar [lindex $Groups $i] hexahedra SteadyStatePwLinkInterfaceElement3D8N $LinkInterfaceElemsProp Hexahedron3D8Connectivities
                }
            }
        } elseif {$IsQuadratic eq 1} {
            for {set i 0} {$i < [llength $Groups]} {incr i} {
                if {[lindex [lindex $Groups $i] 3] eq false} {
                    # Elements Property
                    set InterfaceElemsProp [dict get $PropertyDict [lindex [lindex $Groups $i] 1]]
                    # SteadyStatePwInterfaceElement2D4N twice
                    WriteElementsTwoParts FileVar [lindex $Groups $i] quadrilateral SteadyStatePwInterfaceElement2D4N $InterfaceElemsProp Quadrilateral2D4ConnectivitiesPart1 Quadrilateral2D4ConnectivitiesPart2

                    # SteadyStatePwInterfaceElement3D6N twice
                    WriteElementsFourParts FileVar [lindex $Groups $i] prism SteadyStatePwInterfaceElement3D6N $InterfaceElemsProp \
                                                        PrismInterface3D6ConnectivitiesPart1 \
                                                        PrismInterface3D6ConnectivitiesPart2 \
                                                        PrismInterface3D6ConnectivitiesPart3 \
                                                        PrismInterface3D6ConnectivitiesPart4
                    # SteadyStatePwInterfaceElement3D8N
                    WriteElements FileVar [lindex $Groups $i] hexahedra SteadyStatePwInterfaceElement3D8N $InterfaceElemsProp HexahedronInterface3D8Connectivities
                } else {
                    # Elements Property
                    set LinkInterfaceElemsProp [dict get $PropertyDict [lindex [lindex $Groups $i] 1]]
                    # SteadyStatePwLinkInterfaceElement2D4N twice
                    WriteElementsTwoParts FileVar [lindex $Groups $i] quadrilateral SteadyStatePwLinkInterfaceElement2D4N $InterfaceElemsProp Quadrilateral2D4ConnectivitiesPart1 Quadrilateral2D4ConnectivitiesPart2
                    WriteElements FileVar [lindex $Groups $i] triangle SteadyStatePwLinkInterfaceElement2D4N $LinkInterfaceElemsProp TriangleInterface2D4Connectivities
                    # SteadyStatePwLinkInterfaceElement3D6N
                    # SteadyStatePwLinkInterfaceElement3D6N twice
                    WriteElementsFourParts FileVar [lindex $Groups $i] prism SteadyStatePwLinkInterfaceElement3D6N $InterfaceElemsProp \
                                                        PrismInterface3D6ConnectivitiesPart1 \
                                                        PrismInterface3D6ConnectivitiesPart2 \
                                                        PrismInterface3D6ConnectivitiesPart3 \
                                                        PrismInterface3D6ConnectivitiesPart4
                    WriteElements FileVar [lindex $Groups $i] tetrahedra SteadyStatePwLinkInterfaceElement3D6N $LinkInterfaceElemsProp TetrahedronInterface3D6Connectivities
                    # SteadyStatePwLinkInterfaceElement3D8N
                    WriteElements FileVar [lindex $Groups $i] hexahedra SteadyStatePwLinkInterfaceElement3D8N $LinkInterfaceElemsProp Hexahedron3D8Connectivities
                }
            }
        }
    } else {
        if {$IsQuadratic eq 0} {
            for {set i 0} {$i < [llength $Groups]} {incr i} {
                if {[lindex [lindex $Groups $i] 3] eq false} {
                    # Elements Property
                    set InterfaceElemsProp [dict get $PropertyDict [lindex [lindex $Groups $i] 1]]
                    # TransientPwInterfaceElement2D4N
                    WriteElements FileVar [lindex $Groups $i] quadrilateral TransientPwInterfaceElement2D4N $InterfaceElemsProp Quadrilateral2D4Connectivities

                    # TransientPwInterfaceElement3D6N
                    WriteElements FileVar [lindex $Groups $i] prism TransientPwInterfaceElement3D6N $InterfaceElemsProp PrismInterface3D6Connectivities
                    # TransientPwInterfaceElement3D8N
                    WriteElements FileVar [lindex $Groups $i] hexahedra TransientPwInterfaceElement3D8N $InterfaceElemsProp HexahedronInterface3D8Connectivities
                } else {
                    # Elements Property
                    set LinkInterfaceElemsProp [dict get $PropertyDict [lindex [lindex $Groups $i] 1]]
                    # TransientPwLinkInterfaceElement2D4N
                    WriteElements FileVar [lindex $Groups $i] quadrilateral TransientPwLinkInterfaceElement2D4N $LinkInterfaceElemsProp QuadrilateralInterface2D4Connectivities
                    WriteElements FileVar [lindex $Groups $i] triangle TransientPwLinkInterfaceElement2D4N $LinkInterfaceElemsProp TriangleInterface2D4Connectivities
                    # TransientPwLinkInterfaceElement3D6N
                    WriteElements FileVar [lindex $Groups $i] prism TransientPwLinkInterfaceElement3D6N $LinkInterfaceElemsProp Triangle2D6Connectivities
                    WriteElements FileVar [lindex $Groups $i] tetrahedra TransientPwLinkInterfaceElement3D6N $LinkInterfaceElemsProp TetrahedronInterface3D6Connectivities
                    # TransientPwLinkInterfaceElement3D8N
                    WriteElements FileVar [lindex $Groups $i] hexahedra TransientPwLinkInterfaceElement3D8N $LinkInterfaceElemsProp Hexahedron3D8Connectivities
                }
            }
        } elseif {$IsQuadratic eq 1} {
            for {set i 0} {$i < [llength $Groups]} {incr i} {
                if {[lindex [lindex $Groups $i] 3] eq false} {
                    # Elements Property
                    set InterfaceElemsProp [dict get $PropertyDict [lindex [lindex $Groups $i] 1]]
                    # TransientPwInterfaceElement2D4N twice
                    WriteElementsTwoParts FileVar [lindex $Groups $i] quadrilateral TransientPwInterfaceElement2D4N $InterfaceElemsProp Quadrilateral2D4ConnectivitiesPart1 Quadrilateral2D4ConnectivitiesPart2

                    # TransientPwInterfaceElement3D6N twice
                    WriteElementsFourParts FileVar [lindex $Groups $i] prism TransientPwInterfaceElement3D6N $InterfaceElemsProp \
                                                        PrismInterface3D6ConnectivitiesPart1 \
                                                        PrismInterface3D6ConnectivitiesPart2 \
                                                        PrismInterface3D6ConnectivitiesPart3 \
                                                        PrismInterface3D6ConnectivitiesPart4
                    # TransientPwInterfaceElement3D8N
                    WriteElements FileVar [lindex $Groups $i] hexahedra TransientPwInterfaceElement3D8N $InterfaceElemsProp HexahedronInterface3D8Connectivities
                } else {
                    # Elements Property
                    set LinkInterfaceElemsProp [dict get $PropertyDict [lindex [lindex $Groups $i] 1]]
                    # TransientPwLinkInterfaceElement2D4N twice
                    WriteElementsTwoParts FileVar [lindex $Groups $i] quadrilateral TransientPwLinkInterfaceElement2D4N $InterfaceElemsProp Quadrilateral2D4ConnectivitiesPart1 Quadrilateral2D4ConnectivitiesPart2
                    WriteElements FileVar [lindex $Groups $i] triangle TransientPwLinkInterfaceElement2D4N $LinkInterfaceElemsProp TriangleInterface2D4Connectivities
                    # TransientPwLinkInterfaceElement3D6N
                    # TransientPwLinkInterfaceElement3D6N twice
                    WriteElementsFourParts FileVar [lindex $Groups $i] prism TransientPwLinkInterfaceElement3D6N $InterfaceElemsProp \
                                                        PrismInterface3D6ConnectivitiesPart1 \
                                                        PrismInterface3D6ConnectivitiesPart2 \
                                                        PrismInterface3D6ConnectivitiesPart3 \
                                                        PrismInterface3D6ConnectivitiesPart4
                    WriteElements FileVar [lindex $Groups $i] tetrahedra TransientPwLinkInterfaceElement3D6N $LinkInterfaceElemsProp TetrahedronInterface3D6Connectivities
                    # TransientPwLinkInterfaceElement3D8N
                    WriteElements FileVar [lindex $Groups $i] hexahedra TransientPwLinkInterfaceElement3D8N $LinkInterfaceElemsProp Hexahedron3D8Connectivities
                }
            }
        }

    }

    # PropagationUnion (InterfaceElement)
    if {[GiD_Groups exists PropagationUnion_3d_6] eq 1} {
        # UPwSmallStrainInterfaceElement3D6N
        set PropUnionElementList [WritePropUnionElements FileVar $InterfaceElemsProp]
    }
    puts $FileVar ""

    ## Conditions
    set ConditionId 0
    set ConditionDict [dict create]
    # Point_Load
    set Groups [GiD_Info conditions Point_Load groups]
    if {$Dim eq 2} {
        # UPwForceCondition2D1N
        WriteNodalConditions FileVar ConditionId ConditionDict $Groups PointLoadCondition2D1N $BodyElemsProp
    } else {
        # UPwForceCondition3D1N
        WriteNodalConditions FileVar ConditionId ConditionDict $Groups PointLoadCondition2D1N $BodyElemsProp
    }

    # Line_Load
    set Groups [GiD_Info conditions Line_Load groups]
    if {$Dim eq 2} {
        if {$IsQuadratic eq 0} {
            # LineLoadCondition2D2N
            WriteLineConditions FileVar ConditionId ConditionDict $Groups UPwFaceLoadCondition2D2N $PropertyDict
        } else {
            # LineLoadCondition2D3N
            WriteLineConditions FileVar ConditionId ConditionDict $Groups LineLoadDiffOrderCondition2D3N $PropertyDict
        }
    } else {
        #TODO: needs to developed or adopted from SolidMechanicsApplication
        if {$IsQuadratic eq 0} {
            # LineLoadCondition3D2N
            WriteLineConditions FileVar ConditionId ConditionDict $Groups LineLoadCondition3D2N $PropertyDict
        } else {
            # LineLoadCondition3D3N
            WriteLineConditions FileVar ConditionId ConditionDict $Groups LineLoadCondition3D3N $PropertyDict
        }
    }

    # Surface_Load
    set Groups [GiD_Info conditions Surface_Load groups]
    if {$IsQuadratic eq 0} {
        for {set i 0} {$i < [llength $Groups]} {incr i} {
            set MyConditionList [list]
            # UPwFaceLoadCondition3D3N
            WriteTypeFaceConditions FileVar ConditionId MyConditionList [lindex $Groups $i] tetrahedra UPwFaceLoadCondition3D3N $PropertyDict
            WriteTypeFaceConditions FileVar ConditionId MyConditionList [lindex $Groups $i] prism UPwFaceLoadCondition3D3N $PropertyDict
            # UPwFaceLoadCondition3D4N
            WriteTypeFaceConditions FileVar ConditionId MyConditionList [lindex $Groups $i] hexahedra UPwFaceLoadCondition3D4N $PropertyDict
            dict set ConditionDict [lindex [lindex $Groups $i] 1] $MyConditionList
        }
    } elseif {$IsQuadratic eq 1} {
        for {set i 0} {$i < [llength $Groups]} {incr i} {
            set MyConditionList [list]
            # SurfaceLoadDiffOrderCondition3D6N
            WriteTypeFaceConditions FileVar ConditionId MyConditionList [lindex $Groups $i] tetrahedra SurfaceLoadDiffOrderCondition3D6N $PropertyDict
            # SurfaceLoadDiffOrderCondition3D8N
            WriteTypeFaceConditions FileVar ConditionId MyConditionList [lindex $Groups $i] hexahedra SurfaceLoadDiffOrderCondition3D8N $PropertyDict
            dict set ConditionDict [lindex [lindex $Groups $i] 1] $MyConditionList
        }
    } else {
        for {set i 0} {$i < [llength $Groups]} {incr i} {
            set MyConditionList [list]
            # SurfaceLoadDiffOrderCondition3D6N
            WriteTypeFaceConditions FileVar ConditionId MyConditionList [lindex $Groups $i] tetrahedra SurfaceLoadDiffOrderCondition3D6N $PropertyDict
            # SurfaceLoadDiffOrderCondition3D9N
            WriteTypeFaceConditions FileVar ConditionId MyConditionList [lindex $Groups $i] hexahedra SurfaceLoadDiffOrderCondition3D9N $PropertyDict
            dict set ConditionDict [lindex [lindex $Groups $i] 1] $MyConditionList
        }
    }

    # Normal_Load
    set Groups [GiD_Info conditions Normal_Load groups]
    if {$Dim eq 2} {
        if {$IsQuadratic eq 0} {
            # UPwNormalFaceLoadCondition2D2N
            WriteFaceConditions FileVar ConditionId ConditionDict $Groups UPwNormalFaceLoadCondition2D2N $PropertyDict
        } else {
            # LineNormalLoadDiffOrderCondition2D3N
            WriteFaceConditions FileVar ConditionId ConditionDict $Groups LineNormalLoadDiffOrderCondition2D3N $PropertyDict
        }
    } else {
        if {$IsQuadratic eq 0} {
            for {set i 0} {$i < [llength $Groups]} {incr i} {
                set MyConditionList [list]
                # UPwNormalFaceLoadCondition3D3N
                WriteTypeFaceConditions FileVar ConditionId MyConditionList [lindex $Groups $i] tetrahedra UPwNormalFaceLoadCondition3D3N $PropertyDict
                WriteTypeFaceConditions FileVar ConditionId MyConditionList [lindex $Groups $i] prism UPwNormalFaceLoadCondition3D3N $PropertyDict
                # UpwNormalFaceLoadCondition3D4N
                WriteTypeFaceConditions FileVar ConditionId MyConditionList [lindex $Groups $i] hexahedra UpwNormalFaceLoadCondition3D4N $PropertyDict
                dict set ConditionDict [lindex [lindex $Groups $i] 1] $MyConditionList
            }
        } elseif {$IsQuadratic eq 1} {
            for {set i 0} {$i < [llength $Groups]} {incr i} {
                set MyConditionList [list]
                # SurfaceNormalLoadDiffOrderCondition3D6N
                WriteTypeFaceConditions FileVar ConditionId MyConditionList [lindex $Groups $i] tetrahedra SurfaceNormalLoadDiffOrderCondition3D6N $PropertyDict
                # SurfaceNormalLoadDiffOrderCondition3D8N
                WriteTypeFaceConditions FileVar ConditionId MyConditionList [lindex $Groups $i] hexahedra SurfaceNormalLoadDiffOrderCondition3D8N $PropertyDict
                dict set ConditionDict [lindex [lindex $Groups $i] 1] $MyConditionList
            }
        } else {
            for {set i 0} {$i < [llength $Groups]} {incr i} {
                set MyConditionList [list]
                # SurfaceNormalLoadDiffOrderCondition3D6N
                WriteTypeFaceConditions FileVar ConditionId MyConditionList [lindex $Groups $i] tetrahedra SurfaceNormalLoadDiffOrderCondition3D6N $PropertyDict
                # SurfaceNormalLoadDiffOrderCondition3D9N
                WriteTypeFaceConditions FileVar ConditionId MyConditionList [lindex $Groups $i] hexahedra SurfaceNormalLoadDiffOrderCondition3D9N $PropertyDict
                dict set ConditionDict [lindex [lindex $Groups $i] 1] $MyConditionList
            }
        }
    }

    # Normal_Fluid_Flux
    set Groups [GiD_Info conditions Normal_Fluid_Flux groups]
    if {$Dim eq 2} {
        if {$IsQuadratic eq 0} {
            if {$FIC eq false} {
                # UPwNormalFluxCondition2D2N
                WriteFaceConditions FileVar ConditionId ConditionDict $Groups UPwNormalFluxCondition2D2N $PropertyDict
            } else {
                # UPwNormalFluxFICCondition2D2N
                WriteFaceConditions FileVar ConditionId ConditionDict $Groups UPwNormalFluxFICCondition2D2N $PropertyDict
            }
        } else {
            # LineNormalFluidFluxDiffOrderCondition2D3N
            WriteFaceConditions FileVar ConditionId ConditionDict $Groups LineNormalFluidFluxDiffOrderCondition2D3N $PropertyDict
        }
    } else {
        if {$IsQuadratic eq 0} {
            if {$FIC eq false} {
                for {set i 0} {$i < [llength $Groups]} {incr i} {
                    set MyConditionList [list]
                    # UPwNormalFluxCondition3D3N
                    WriteTypeFaceConditions FileVar ConditionId MyConditionList [lindex $Groups $i] tetrahedra UPwNormalFluxCondition3D3N $PropertyDict
                    WriteTypeFaceConditions FileVar ConditionId MyConditionList [lindex $Groups $i] prism UPwNormalFluxCondition3D3N $PropertyDict
                    # UPwNormalFluxCondition3D4N
                    WriteTypeFaceConditions FileVar ConditionId MyConditionList [lindex $Groups $i] hexahedra UPwNormalFluxCondition3D4N $PropertyDict
                    dict set ConditionDict [lindex [lindex $Groups $i] 1] $MyConditionList
                }
            } else {
                for {set i 0} {$i < [llength $Groups]} {incr i} {
                    set MyConditionList [list]
                    # UPwNormalFluxFICCondition3D3N
                    WriteTypeFaceConditions FileVar ConditionId MyConditionList [lindex $Groups $i] tetrahedra UPwNormalFluxFICCondition3D3N $PropertyDict
                    WriteTypeFaceConditions FileVar ConditionId MyConditionList [lindex $Groups $i] prism UPwNormalFluxFICCondition3D3N $PropertyDict
                    # UPwNormalFluxFICCondition3D4N
                    WriteTypeFaceConditions FileVar ConditionId MyConditionList [lindex $Groups $i] hexahedra UPwNormalFluxFICCondition3D4N $PropertyDict
                    dict set ConditionDict [lindex [lindex $Groups $i] 1] $MyConditionList
                }
            }
        } elseif {$IsQuadratic eq 1} {
            for {set i 0} {$i < [llength $Groups]} {incr i} {
                set MyConditionList [list]
                # SurfaceNormalFluidFluxDiffOrderCondition3D6N
                WriteTypeFaceConditions FileVar ConditionId MyConditionList [lindex $Groups $i] tetrahedra SurfaceNormalFluidFluxDiffOrderCondition3D6N $PropertyDict
                # SurfaceNormalFluidFluxDiffOrderCondition3D8N
                WriteTypeFaceConditions FileVar ConditionId MyConditionList [lindex $Groups $i] hexahedra SurfaceNormalFluidFluxDiffOrderCondition3D8N $PropertyDict
                dict set ConditionDict [lindex [lindex $Groups $i] 1] $MyConditionList
            }
        } else {
            for {set i 0} {$i < [llength $Groups]} {incr i} {
                set MyConditionList [list]
                # SurfaceNormalFluidFluxDiffOrderCondition3D6N
                WriteTypeFaceConditions FileVar ConditionId MyConditionList [lindex $Groups $i] tetrahedra SurfaceNormalFluidFluxDiffOrderCondition3D6N $PropertyDict
                # SurfaceNormalFluidFluxDiffOrderCondition3D9N
                WriteTypeFaceConditions FileVar ConditionId MyConditionList [lindex $Groups $i] hexahedra SurfaceNormalFluidFluxDiffOrderCondition3D9N $PropertyDict
                dict set ConditionDict [lindex [lindex $Groups $i] 1] $MyConditionList
            }
        }
    }

    # Interface_Face_Load
    set Groups [GiD_Info conditions Interface_Face_Load groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
        set MyConditionList [list]
        # UPwFaceLoadInterfaceCondition2D2N
        WriteInterfaceConditions FileVar ConditionId MyConditionList [lindex $Groups $i] linear UPwFaceLoadInterfaceCondition2D2N $InterfaceElemsProp Line2D2Connectivities
        # UPwFaceLoadInterfaceCondition3D4N
        WriteInterfaceConditions FileVar ConditionId MyConditionList [lindex $Groups $i] triangle UPwFaceLoadInterfaceCondition3D4N $InterfaceElemsProp TriangleInterface3D4Connectivities
        WriteInterfaceConditions FileVar ConditionId MyConditionList [lindex $Groups $i] quadrilateral UPwFaceLoadInterfaceCondition3D4N $InterfaceElemsProp QuadrilateralInterface3D4Connectivities
        dict set ConditionDict [lindex [lindex $Groups $i] 1] $MyConditionList
    }

    # Interface_Normal_Fluid_Flux
    set Groups [GiD_Info conditions Interface_Normal_Fluid_Flux groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
        set MyConditionList [list]
        # UPwNormalFluxInterfaceCondition2D2N
        WriteInterfaceConditions FileVar ConditionId MyConditionList [lindex $Groups $i] linear UPwNormalFluxInterfaceCondition2D2N $InterfaceElemsProp Line2D2Connectivities
        # UPwNormalFluxInterfaceCondition3D4N
        WriteInterfaceConditions FileVar ConditionId MyConditionList [lindex $Groups $i] triangle UPwNormalFluxInterfaceCondition3D4N $InterfaceElemsProp TriangleInterface3D4Connectivities
        WriteInterfaceConditions FileVar ConditionId MyConditionList [lindex $Groups $i] quadrilateral UPwNormalFluxInterfaceCondition3D4N $InterfaceElemsProp QuadrilateralInterface3D4Connectivities
        dict set ConditionDict [lindex [lindex $Groups $i] 1] $MyConditionList
    }

    # Gap_Closure
    set IsGapClosure [GiD_AccessValue get gendata Gap_Closure_Interface_Conditions]
    if {$IsGapClosure eq true} {
        set GapClosureBarsDict [dict create]

        set interface_Groups [list [GiD_Info conditions Interface_two_phase groups] [GiD_Info conditions Interface_drained groups] [GiD_Info conditions Interface_undrained groups]]
        foreach Groups $interface_Groups {
            for {set i 0} {$i < [llength $Groups]} {incr i} {
                if {[lindex [lindex $Groups $i] 135] eq true} {
                    # Elements Property
                    set InterfaceElemsProp [dict get $PropertyDict [lindex [lindex $Groups $i] 1]]
                    set ConditionList [list]
                    # InterfaceElement2D4N
                    SaveGapClosureBarsFromIE2D4N GapClosureBarsDict ConditionId ConditionList [lindex $Groups $i] $InterfaceElemsProp
                    # InterfaceElement3D6N
                    SaveGapClosureBarsFromIE3D6N GapClosureBarsDict ConditionId ConditionList [lindex $Groups $i] $InterfaceElemsProp
                    # InterfaceElement3D8N
                    SaveGapClosureBarsFromIE3D8N GapClosureBarsDict ConditionId ConditionList [lindex $Groups $i] $InterfaceElemsProp

                    dict set ConditionDict Gap_Closure_Bars_[lindex [lindex $Groups $i] 1] $ConditionList
                }
            }
        }

        # if {[dict size $GapClosureBarsDict] > 0} {
        #     puts $FileVar "Begin Conditions GapClosureCondition"
        #     dict for {Name GapClosureBar} $GapClosureBarsDict {
        #         puts $FileVar "  [dict get $GapClosureBar Id]  [dict get $GapClosureBar PropertyId]  [dict get $GapClosureBar Connectivities]"
        #     }
        #     puts $FileVar "End Conditions"
        #     puts $FileVar ""
        # }
    }


    puts $FileVar ""

    ## SubModelParts
    # Soil_two_phase part
    WriteElementSubmodelPart FileVar Soil_two_phase
    # Soil_drained part
    WriteElementSubmodelPart FileVar Soil_drained
    # Soil_undrained part
    WriteElementSubmodelPart FileVar Soil_undrained
    # Non_porous part
    WriteElementSubmodelPart FileVar Non_porous
    # Soil_Groundwater_Flow part
    WriteElementSubmodelPart FileVar Soil_Groundwater_Flow
    # Beam part
    WriteElementSubmodelPart FileVar Beam
    # Shell_thin_corotational part
    WriteElementSubmodelPart FileVar Shell_thin_corotational
    # Shell_thick_corotational part
    WriteElementSubmodelPart FileVar Shell_thick_corotational
    # Truss part
    WriteElementSubmodelPart FileVar Truss
    # Anchor part
    # WriteAnchorElementSubmodelPart FileVar Anchor $AnchorElementDict
    WriteElementSubmodelPart FileVar Anchor
    # Interface drained Part
    WriteElementSubmodelPart FileVar Interface_drained
    # Interface undrained Part
    WriteElementSubmodelPart FileVar Interface_undrained
    # Interface two_phase Part
    WriteElementSubmodelPart FileVar Interface_two_phase
    # Interface_Groundwater_flow Part
    WriteElementSubmodelPart FileVar Interface_Groundwater_flow

    # PropagationUnion (InterfaceElement)
    if {[GiD_Groups exists PropagationUnion_3d_6] eq 1} {
        WritePropUnionElementSubmodelPart FileVar $PropUnionElementList
    }
    # Solid_Displacement
    WriteConstraintSubmodelPart FileVar Solid_Displacement $TableDict
    # Structural_Rotation
    WriteConstraintSubmodelPart FileVar Structural_Rotation $TableDict
    # Fluid_Pressure
    WriteConstraintSubmodelPart FileVar Fluid_Pressure $TableDict
    # Excavation, todo add conditions to excavation
    # WriteExcavationSubmodelPartOriginal FileVar Excavation $AnchorElementDict
    WriteExcavationSubmodelPart FileVar Excavation
    # Point_Load
    WriteLoadSubmodelPart FileVar Point_Load $TableDict $ConditionDict
    # Line_Load
    WriteLoadSubmodelPart FileVar Line_Load $TableDict $ConditionDict
    # Surface_Load
    WriteLoadSubmodelPart FileVar Surface_Load $TableDict $ConditionDict
    # Normal_Load
    WriteLoadSubmodelPart FileVar Normal_Load $TableDict $ConditionDict
    # Normal_Fluid_Flux
    WriteLoadSubmodelPart FileVar Normal_Fluid_Flux $TableDict $ConditionDict
    # Interface_Face_Load
    WriteLoadSubmodelPart FileVar Interface_Face_Load $TableDict $ConditionDict
    # Interface_Normal_Fluid_Flux
    WriteLoadSubmodelPart FileVar Interface_Normal_Fluid_Flux $TableDict $ConditionDict
    # Body_Acceleration
    WriteConstraintSubmodelPart FileVar Body_Acceleration $TableDict
    # Record_DISPLACEMENT
    WriteRecordResultSubmodelPart FileVar Record_DISPLACEMENT
    # Record_VELOCITY
    WriteRecordResultSubmodelPart FileVar Record_VELOCITY
    # Record_ACCELERATION
    WriteRecordResultSubmodelPart FileVar Record_ACCELERATION
    # Record_VOLUME_ACCELERATION
    WriteRecordResultSubmodelPart FileVar Record_VOLUME_ACCELERATION
    # Record_POINT_LOAD
    WriteRecordResultSubmodelPart FileVar Record_POINT_LOAD
    # Record_LINE_LOAD
    WriteRecordResultSubmodelPart FileVar Record_LINE_LOAD
    # Record_SURFACE_LOAD
    WriteRecordResultSubmodelPart FileVar Record_SURFACE_LOAD

    # GapClosure_Bars
    if {$IsGapClosure eq true} {
        WriteGapClosureBarsSubmodelPart FileVar Interface_two_phase $ConditionDict
        WriteGapClosureBarsSubmodelPart FileVar Interface_drained $ConditionDict
        WriteGapClosureBarsSubmodelPart FileVar Interface_undrained $ConditionDict
    }

    close $FileVar

    return $TableDict
}

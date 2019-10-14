
proc WriteMdpa { basename dir problemtypedir } {
    
    ## Source auxiliar procedures
    source [file join $problemtypedir MdpaAuxProcs.tcl]
    
    ## Start MDPA file
    set filename [file join $dir ${basename}.mdpa]
    set FileVar [open $filename w]

    puts $FileVar ""
    puts $FileVar "Begin Properties 0"
    puts $FileVar "End Properties"
    puts $FileVar ""
    
    ## Properties
    set PropertyId 0
    set PropertyDict [dict create]
    
    # Solid_Displacement
    ConstraintVectorTable FileVar TableId TableDict Solid_Displacement DISPLACEMENT
    # Force
    VectorTable FileVar TableId TableDict Force FORCE
    # Face_Load
    VectorTable FileVar TableId TableDict Face_Load FACE_LOAD
    # Normal_Load
    #NormalTangentialTable FileVar TableId TableDict Normal_Load NORMAL_CONTACT_STRESS TANGENTIAL_CONTACT_STRESS

    # Surf_load
    VectorTable FileVar TableId TableDict Surf_Load SURF_LOAD

    # Pressure_Load
    NormalTangentialTable FileVar TableId TableDict Pressure_Load NORMAL_CONTACT_STRESS TANGENTIAL_CONTACT_STRESS

    # Body_Acceleration
    VectorTable FileVar TableId TableDict Body_Acceleration VOLUME_ACCELERATION

    set Groups [GiD_Info conditions Body_Part groups]


    for {set i 0} {$i < [llength $Groups]} {incr i} {
        incr PropertyId
        dict set PropertyDict [lindex [lindex $Groups $i] 1] $PropertyId
        puts $FileVar "Begin Properties $PropertyId"
        puts $FileVar "// ELASTIC PROPERTIES"
        puts $FileVar "    YOUNG_MODULUS            [lindex [lindex $Groups $i] 4]"
        puts $FileVar "    DENSITY                  [lindex [lindex $Groups $i] 5]"
        puts $FileVar "    POISSON_RATIO            [lindex [lindex $Groups $i] 6]"
		
		if {[lindex [lindex $Groups $i] 13] eq "true"} {
			puts $FileVar "    YOUNG_MODULUS_STEEL      [lindex [lindex $Groups $i] 14]"
			puts $FileVar "    DENSITY_STEEL            [lindex [lindex $Groups $i] 15]"
			puts $FileVar "    POISSON_RATIO_STEEL      [lindex [lindex $Groups $i] 16]"
			puts $FileVar "    STEEL_VOLUMETRIC_PART    [lindex [lindex $Groups $i] 17]"
			puts $FileVar "    YIELD_STRESS_C_STEEL     [lindex [lindex $Groups $i] 18]"
			puts $FileVar "    YIELD_STRESS_T_STEEL     [lindex [lindex $Groups $i] 19]"
			puts $FileVar "    FRACTURE_ENERGY_STEEL    [lindex [lindex $Groups $i] 20]"
			puts $FileVar "    HARDENING_LAW            [lindex [lindex $Groups $i] 21]"
			
			if {[lindex [lindex $Groups $i] 21] eq 3} {
				puts $FileVar "    MAXIMUM_STRESS           [lindex [lindex $Groups $i] 22]"
				puts $FileVar "    MAXIMUM_STRESS_POSITION  [lindex [lindex $Groups $i] 23]"
			}
		}
        puts $FileVar ""
        puts $FileVar "// DAMAGE PARAMETERS"
        puts $FileVar "    YIELD_SURFACE            [lindex [lindex $Groups $i] 3]"
        puts $FileVar "    YIELD_STRESS_C           [lindex [lindex $Groups $i] 7]"
        puts $FileVar "    YIELD_STRESS_T           [lindex [lindex $Groups $i] 8]"
        puts $FileVar "    FRAC_ENERGY_T            [lindex [lindex $Groups $i] 9]"
        puts $FileVar "    INTERNAL_FRICTION_ANGLE  [lindex [lindex $Groups $i] 10]"
        puts $FileVar ""
        puts $FileVar "// DYNAMIC PARAMETERS  D = alpha*M + beta*K"
        puts $FileVar "    RAYLEIGH_BETA            [lindex [lindex $Groups $i] 11]"
        puts $FileVar "    RAYLEIGH_ALPHA           [lindex [lindex $Groups $i] 12]"
        puts $FileVar "End Properties"
        puts $FileVar ""
    }
    puts $FileVar ""

    ## Nodes
    set Nodes [GiD_Info Mesh Nodes]
    puts $FileVar "Begin Nodes"
    for {set i 0} {$i < [llength $Nodes]} {incr i 4} {
        puts $FileVar "  [lindex $Nodes $i]    [format  "%.10f" [lindex $Nodes [expr { $i+1 }]]]    [format  "%.10f" [lindex $Nodes [expr { $i+2 }]]]    [format  "%.10f" [lindex $Nodes [expr { $i+3 }]]]"
    }
    puts $FileVar "End Nodes"
    puts $FileVar ""
    puts $FileVar ""
 
    ## Elements
    set Groups [GiD_Info conditions Body_Part groups]

    for {set i 0} {$i < [llength $Groups]} {incr i} {

        set ElementName ""
        if {[lindex [lindex $Groups $i] 24] eq "false"} {
            if {[lindex [lindex $Groups $i] 3] eq "ModifiedMohrCoulomb"} {
                set ElementName "SmallStrainModifiedMohrCoulombFemDemElement3D"
            } elseif {[lindex [lindex $Groups $i] 3] eq "Rankine"} {
                set ElementName "SmallStrainRankineFemDemElement3D"
            } elseif {[lindex [lindex $Groups $i] 3] eq "SimoJu"} {
                set ElementName "SmallStrainSimoJuFemDemElement3D"
            } elseif {[lindex [lindex $Groups $i] 3] eq "DruckerPrager"} {
                set ElementName "SmallStrainDruckerPragerFemDemElement3D"
            } elseif {[lindex [lindex $Groups $i] 3] eq "VonMises"} {
                set ElementName "SmallStrainVonMisesFemDemElement3D"
            } elseif {[lindex [lindex $Groups $i] 3] eq "Tresca"} {
                set ElementName "SmallStrainTrescaFemDemElement3D"
            } elseif {[lindex [lindex $Groups $i] 3] eq "MohrCoulomb"} {
                set ElementName "SmallStrainMohrCoulombFemDemElement3D"
            } else {
                set ElementName "SmallStrainModifiedMohrCoulombFemDemElement3D"
            }        
        } else {
            if {[lindex [lindex $Groups $i] 3] eq "ModifiedMohrCoulomb"} {
                set ElementName "TotalLagrangianModifiedMohrCoulombFemDemElement3D"
            } elseif {[lindex [lindex $Groups $i] 3] eq "Rankine"} {
                set ElementName "TotalLagrangianRankineFemDemElement3D"
            } elseif {[lindex [lindex $Groups $i] 3] eq "SimoJu"} {
                set ElementName "TotalLagrangianSimoJuFemDemElement3D"
            } elseif {[lindex [lindex $Groups $i] 3] eq "DruckerPrager"} {
                set ElementName "TotalLagrangianDruckerPragerFemDemElement3D"
            } elseif {[lindex [lindex $Groups $i] 3] eq "VonMises"} {
                set ElementName "TotalLagrangianVonMisesFemDemElement3D"
            } elseif {[lindex [lindex $Groups $i] 3] eq "Tresca"} {
                set ElementName "TotalLagrangianTrescaFemDemElement3D"
            } elseif {[lindex [lindex $Groups $i] 3] eq "MohrCoulomb"} {
                set ElementName "TotalLagrangianMohrCoulombFemDemElement3D"
            } else {
                set ElementName "TotalLagrangianModifiedMohrCoulombFemDemElement3D"
            }  
        }  
         # Elements Property
        set BodyElemsProp [dict get $PropertyDict [lindex [lindex $Groups $i] 1]]
		
		set MatGroups [GiD_Info conditions Body_Part groups]
		if {[lindex [lindex $MatGroups 0] 13] eq "true"} {
			WriteElements FileVar [lindex $Groups $i] tetrahedra RomFemDem3DElement $BodyElemsProp Tetrahedron3D4Connectivities
		} elseif {[GiD_AccessValue get gendata Use_Hexahedrons] eq "true"} {
            WriteElements FileVar [lindex $Groups $i] Hexahedra FemDem3DHexahedronElement $BodyElemsProp Hexahedron3D8Connectivities
		} else {
			WriteElements FileVar [lindex $Groups $i] tetrahedra $ElementName $BodyElemsProp Tetrahedron3D4Connectivities
        }
	}
    puts $FileVar ""

    ## Conditions
    set ConditionId 0
    set ConditionDict [dict create]
    set Dim [GiD_AccessValue get gendata Domain_Size]

    # Force
    set Groups [GiD_Info conditions Force groups]
    WriteNodalConditions FileVar ConditionId ConditionDict $Groups PointLoadCondition3D1N $BodyElemsProp
    
    # Face_Load
    set Groups [GiD_Info conditions Face_Load groups]
    WriteFaceConditions FileVar ConditionId ConditionDict $Groups LineLoadCondition3D2N $PropertyDict

    # Normal_Load
    set Groups [GiD_Info conditions Normal_Load groups]
    WriteFaceConditions FileVar ConditionId ConditionDict $Groups LineLoadCondition3D2N $PropertyDict

    # Surf_Load
    set Groups [GiD_Info conditions Surf_Load groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {

        set MyConditionList [list]
        WriteTypeFaceConditions FileVar ConditionId MyConditionList [lindex $Groups $i] tetrahedra SurfaceLoadCondition3D3N $PropertyDict
        dict set ConditionDict [lindex [lindex $Groups $i] 1] $MyConditionList
    }
    
    # Pressure_Load
    set Groups [GiD_Info conditions Pressure_Load groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {

        set MyConditionList [list]
        WriteTypeFaceConditions FileVar ConditionId MyConditionList [lindex $Groups $i] tetrahedra SurfaceLoadCondition3D3N $PropertyDict
        dict set ConditionDict [lindex [lindex $Groups $i] 1] $MyConditionList
    }

    puts $FileVar ""
    puts $FileVar ""
    
    ## SubModelParts
    # Body_Part
    WriteElementSubmodelPart FileVar Body_Part
    # Solid_Displacement
    WriteConstraintSubmodelPart FileVar Solid_Displacement $TableDict
    # Force
    WriteLoadSubmodelPart FileVar Force $TableDict $ConditionDict
    # Face_Load
    WriteLoadSubmodelPart FileVar Face_Load $TableDict $ConditionDict
    # Normal_Load
    #WriteLoadSubmodelPart FileVar Normal_Load $TableDict $ConditionDict
    # Body_Acceleration
    WriteConstraintSubmodelPart FileVar Body_Acceleration $TableDict

    # Surf_Load
    WriteLoadSubmodelPart FileVar Surf_Load $TableDict $ConditionDict
    # Pressure_Load
    WriteLoadSubmodelPart FileVar Pressure_Load $TableDict $ConditionDict
    
    close $FileVar
    
    return $TableDict

 
    
}
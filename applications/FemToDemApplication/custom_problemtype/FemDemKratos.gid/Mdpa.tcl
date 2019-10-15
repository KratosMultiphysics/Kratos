
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
    NormalTangentialTable FileVar TableId TableDict Normal_Load NORMAL_CONTACT_STRESS TANGENTIAL_CONTACT_STRESS
    # Body_Acceleration
    VectorTable FileVar TableId TableDict Body_Acceleration VOLUME_ACCELERATION
    set Groups [GiD_Info conditions Body_Part groups]
    #W "Grupos  [lindex [lindex $Groups 1] 3]"
    #puts $FileVar "testeo  [lindex [lindex $Groups 1] 3]"

    for {set i 0} {$i < [llength $Groups]} {incr i} {
        incr PropertyId
        dict set PropertyDict [lindex [lindex $Groups $i] 1] $PropertyId
        puts $FileVar "Begin Properties $PropertyId"
        puts $FileVar "// ELASTIC PROPERTIES"
        puts $FileVar "    YOUNG_MODULUS            [lindex [lindex $Groups $i] 4]"
        puts $FileVar "    DENSITY                  [lindex [lindex $Groups $i] 5]"
        puts $FileVar "    POISSON_RATIO            [lindex [lindex $Groups $i] 6]"
        puts $FileVar "    THICKNESS                [lindex [lindex $Groups $i] 8]"
        puts $FileVar ""
        puts $FileVar "// DAMAGE PARAMETERS"
        puts $FileVar "    YIELD_SURFACE            [lindex [lindex $Groups $i] 3]"
        puts $FileVar "    YIELD_STRESS_C           [lindex [lindex $Groups $i] 9]"
        puts $FileVar "    YIELD_STRESS_T           [lindex [lindex $Groups $i] 10]"
        puts $FileVar "    FRAC_ENERGY_T            [lindex [lindex $Groups $i] 11]"
        puts $FileVar "    INTERNAL_FRICTION_ANGLE  [lindex [lindex $Groups $i] 12]"
		#puts $FileVar "    CHARACTERISTIC_LENGTH    [lindex [lindex $Groups $i] 15]"
        puts $FileVar ""
        puts $FileVar "// DYNAMIC PARAMETERS  D = alpha*M + beta*K"
        puts $FileVar "    RAYLEIGH_BETA            [lindex [lindex $Groups $i] 13]"
        puts $FileVar "    RAYLEIGH_ALPHA           [lindex [lindex $Groups $i] 14]"
        puts $FileVar "End Properties"
        puts $FileVar ""
    }
    puts $FileVar ""
    #Warnwin "Grupos $Groups"
    #W "Grupos $Groups"

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
        if {[lindex [lindex $Groups $i] 15] eq "false"} {
            if {[lindex [lindex $Groups $i] 3] eq "ModifiedMohrCoulomb"} {
                set ElementName "SmallStrainModifiedMohrCoulombFemDemElement2D"
            } elseif {[lindex [lindex $Groups $i] 3] eq "Rankine"} {
                set ElementName "SmallStrainRankineFemDemElement2D"
            } elseif {[lindex [lindex $Groups $i] 3] eq "SimoJu"} {
                set ElementName "SmallStrainSimoJuFemDemElement2D"
            } elseif {[lindex [lindex $Groups $i] 3] eq "DruckerPrager"} {
                set ElementName "SmallStrainDruckerPragerFemDemElement2D"
            } elseif {[lindex [lindex $Groups $i] 3] eq "VonMises"} {
                set ElementName "SmallStrainVonMisesFemDemElement2D"
            } elseif {[lindex [lindex $Groups $i] 3] eq "Tresca"} {
                set ElementName "SmallStrainTrescaFemDemElement2D"
            }  elseif {[lindex [lindex $Groups $i] 3] eq "MohrCoulomb"} {
                set ElementName "SmallStrainMohrCoulombFemDemElement2D"
            } else {
                set ElementName "SmallStrainModifiedMohrCoulombFemDemElement2D"
            }        
        } else {
            if {[lindex [lindex $Groups $i] 3] eq "ModifiedMohrCoulomb"} {
                set ElementName "TotalLagrangianModifiedMohrCoulombFemDemElement2D"
            } elseif {[lindex [lindex $Groups $i] 3] eq "Rankine"} {
                set ElementName "TotalLagrangianRankineFemDemElement2D"
            } elseif {[lindex [lindex $Groups $i] 3] eq "SimoJu"} {
                set ElementName "TotalLagrangianSimoJuFemDemElement2D"
            } elseif {[lindex [lindex $Groups $i] 3] eq "DruckerPrager"} {
                set ElementName "TotalLagrangianDruckerPragerFemDemElement2D"
            } elseif {[lindex [lindex $Groups $i] 3] eq "VonMises"} {
                set ElementName "TotalLagrangianVonMisesFemDemElement2D"
            } elseif {[lindex [lindex $Groups $i] 3] eq "Tresca"} {
                set ElementName "TotalLagrangianTrescaFemDemElement2D"
            }   elseif {[lindex [lindex $Groups $i] 3] eq "MohrCoulomb"} {
                set ElementName "TotalLagrangianMohrCoulombFemDemElement2D"
            } else {
                set ElementName "TotalLagrangianModifiedMohrCoulombFemDemElement2D"
            }  
        }
         # Elements Property
        set BodyElemsProp [dict get $PropertyDict [lindex [lindex $Groups $i] 1]]
        WriteElements FileVar [lindex $Groups $i] triangle $ElementName $BodyElemsProp Triangle2D3Connectivities
    }
    puts $FileVar ""

    ## Conditions
    set ConditionId 0
    set ConditionDict [dict create]
    set Dim [GiD_AccessValue get gendata Domain_Size]
    # Force
    set Groups [GiD_Info conditions Force groups]
    WriteNodalConditions FileVar ConditionId ConditionDict $Groups PointLoadCondition2D1N $BodyElemsProp
    # Face_Load
    set Groups [GiD_Info conditions Face_Load groups]
    WriteFaceConditions FileVar ConditionId ConditionDict $Groups LineLoadCondition2D2N $PropertyDict
    # Normal_Load
    set Groups [GiD_Info conditions Normal_Load groups]
    WriteFaceConditions FileVar ConditionId ConditionDict $Groups LineLoadCondition2D2N $PropertyDict
    #puts $FileVar "End Conditions"
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
    WriteLoadSubmodelPart FileVar Normal_Load $TableDict $ConditionDict
    # Body_Acceleration
    WriteConstraintSubmodelPart FileVar Body_Acceleration $TableDict
    
    close $FileVar
    
    return $TableDict

 
    
}
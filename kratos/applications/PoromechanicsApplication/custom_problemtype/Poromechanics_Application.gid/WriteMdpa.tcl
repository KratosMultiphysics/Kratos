proc WriteMdpa { basename dir } {
    
    ## Start MDPA file
    set filename [file join $dir ${basename}.mdpa]
    set varfile [open $filename w]
    
    ## ModelPart Data
    #puts $varfile "Begin ModelPartData"
    #puts $varfile "  // VARIABLE_NAME value"
    #puts $varfile "End ModelPartData"
    #puts $varfile ""
    #puts $varfile ""
    
    ## Tables
    set TableId 0
    set AuxList [list]
    set TableList [list]
    # Solid_Displacement
    set Groups [GiD_Info conditions Solid_Displacement groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
        set SearchInList [lsearch [lindex $Groups $i] Table_Interpolation_X*]
        if {$SearchInList > -1} {
            incr TableId
            lappend AuxList $TableId
            puts $varfile "Begin Table $TableId TIME DISPLACEMENT_X"
            set Table [lindex [lindex $Groups $i] [expr { $SearchInList+1 }]]
            for {set j 2} {$j <= [lindex $Table 1]} {incr j 2} {
                puts $varfile "  [lindex $Table $j] [lindex $Table [expr { $j+1 }]]"
            }
            puts $varfile "End Table"
            puts $varfile ""
        } else {
            lappend AuxList 0
        }
        set SearchInList [lsearch [lindex $Groups $i] Table_Interpolation_Y*]
        if {$SearchInList > -1} {
            incr TableId
            lappend AuxList $TableId
            puts $varfile "Begin Table $TableId TIME DISPLACEMENT_Y"
            set Table [lindex [lindex $Groups $i] [expr { $SearchInList+1 }]]
            for {set j 2} {$j <= [lindex $Table 1]} {incr j 2} {
                puts $varfile "  [lindex $Table $j] [lindex $Table [expr { $j+1 }]]"
            }
            puts $varfile "End Table"
            puts $varfile ""
        } else {
            lappend AuxList 0
        }
        set SearchInList [lsearch [lindex $Groups $i] Table_Interpolation_Z*]
        if {$SearchInList > -1} {
            incr TableId
            lappend AuxList $TableId
            puts $varfile "Begin Table $TableId TIME DISPLACEMENT_Z"
            set Table [lindex [lindex $Groups $i] [expr { $SearchInList+1 }]]
            for {set j 2} {$j <= [lindex $Table 1]} {incr j 2} {
                puts $varfile "  [lindex $Table $j] [lindex $Table [expr { $j+1 }]]"
            }
            puts $varfile "End Table"
            puts $varfile ""
        } else {
            lappend AuxList 0
        }
        lappend TableList [lindex [lindex $Groups $i] 1] $AuxList
        set AuxList [list]
    }
    # Fluid_Pressure
    set Groups [GiD_Info conditions Fluid_Pressure groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
        set SearchInList [lsearch [lindex $Groups $i] Table_Interpolation*]
        if {$SearchInList > -1} {
            incr TableId
            lappend TableList [lindex [lindex $Groups $i] 1] $TableId
            puts $varfile "Begin Table $TableId TIME WATER_PRESSURE"
            set Table [lindex [lindex $Groups $i] [expr { $SearchInList+1 }]]
            for {set j 2} {$j <= [lindex $Table 1]} {incr j 2} {
                puts $varfile "  [lindex $Table $j] [lindex $Table [expr { $j+1 }]]"
            }
            puts $varfile "End Table"
            puts $varfile ""
        } else {
            lappend TableList [lindex [lindex $Groups $i] 1] 0
        }
    }
    # Force
    set Groups [GiD_Info conditions Force groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
        set SearchInList [lsearch [lindex $Groups $i] Table_Interpolation_X*]
        if {$SearchInList > -1} {
            incr TableId
            lappend AuxList $TableId
            puts $varfile "Begin Table $TableId TIME FORCE_X"
            set Table [lindex [lindex $Groups $i] [expr { $SearchInList+1 }]]
            for {set j 2} {$j <= [lindex $Table 1]} {incr j 2} {
                puts $varfile "  [lindex $Table $j] [lindex $Table [expr { $j+1 }]]"
            }
            puts $varfile "End Table"
            puts $varfile ""
        } else {
            lappend AuxList 0
        }
        set SearchInList [lsearch [lindex $Groups $i] Table_Interpolation_Y*]
        if {$SearchInList > -1} {
            incr TableId
            lappend AuxList $TableId
            puts $varfile "Begin Table $TableId TIME FORCE_Y"
            set Table [lindex [lindex $Groups $i] [expr { $SearchInList+1 }]]
            for {set j 2} {$j <= [lindex $Table 1]} {incr j 2} {
                puts $varfile "  [lindex $Table $j] [lindex $Table [expr { $j+1 }]]"
            }
            puts $varfile "End Table"
            puts $varfile ""
        } else {
            lappend AuxList 0
        }
        set SearchInList [lsearch [lindex $Groups $i] Table_Interpolation_Z*]
        if {$SearchInList > -1} {
            incr TableId
            lappend AuxList $TableId
            puts $varfile "Begin Table $TableId TIME FORCE_Z"
            set Table [lindex [lindex $Groups $i] [expr { $SearchInList+1 }]]
            for {set j 2} {$j <= [lindex $Table 1]} {incr j 2} {
                puts $varfile "  [lindex $Table $j] [lindex $Table [expr { $j+1 }]]"
            }
            puts $varfile "End Table"
            puts $varfile ""
        } else {
            lappend AuxList 0
        }
        lappend TableList [lindex [lindex $Groups $i] 1] $AuxList
        set AuxList [list]
    }
    # Face_Load
    set Groups [GiD_Info conditions Face_Load groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
        set SearchInList [lsearch [lindex $Groups $i] Table_Interpolation_X*]
        if {$SearchInList > -1} {
            incr TableId
            lappend AuxList $TableId
            puts $varfile "Begin Table $TableId TIME FACE_LOAD_X"
            set Table [lindex [lindex $Groups $i] [expr { $SearchInList+1 }]]
            for {set j 2} {$j <= [lindex $Table 1]} {incr j 2} {
                puts $varfile "  [lindex $Table $j] [lindex $Table [expr { $j+1 }]]"
            }
            puts $varfile "End Table"
            puts $varfile ""
        } else {
            lappend AuxList 0
        }
        set SearchInList [lsearch [lindex $Groups $i] Table_Interpolation_Y*]
        if {$SearchInList > -1} {
            incr TableId
            lappend AuxList $TableId
            puts $varfile "Begin Table $TableId TIME FACE_LOAD_Y"
            set Table [lindex [lindex $Groups $i] [expr { $SearchInList+1 }]]
            for {set j 2} {$j <= [lindex $Table 1]} {incr j 2} {
                puts $varfile "  [lindex $Table $j] [lindex $Table [expr { $j+1 }]]"
            }
            puts $varfile "End Table"
            puts $varfile ""
        } else {
            lappend AuxList 0
        }
        set SearchInList [lsearch [lindex $Groups $i] Table_Interpolation_Z*]
        if {$SearchInList > -1} {
            incr TableId
            lappend AuxList $TableId
            puts $varfile "Begin Table $TableId TIME FACE_LOAD_Z"
            set Table [lindex [lindex $Groups $i] [expr { $SearchInList+1 }]]
            for {set j 2} {$j <= [lindex $Table 1]} {incr j 2} {
                puts $varfile "  [lindex $Table $j] [lindex $Table [expr { $j+1 }]]"
            }
            puts $varfile "End Table"
            puts $varfile ""
        } else {
            lappend AuxList 0
        }
        lappend TableList [lindex [lindex $Groups $i] 1] $AuxList
        set AuxList [list]
    }
    # Normal_Load
    set Groups [GiD_Info conditions Normal_Load groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
        set SearchInList [lsearch [lindex $Groups $i] Table_Interpolation_Normal*]
        if {$SearchInList > -1} {
            incr TableId
            lappend AuxList $TableId
            puts $varfile "Begin Table $TableId TIME NORMAL_CONTACT_STRESS"
            set Table [lindex [lindex $Groups $i] [expr { $SearchInList+1 }]]
            for {set j 2} {$j <= [lindex $Table 1]} {incr j 2} {
                puts $varfile "  [lindex $Table $j] [lindex $Table [expr { $j+1 }]]"
            }
            puts $varfile "End Table"
            puts $varfile ""
        } else {
            lappend AuxList 0
        }
        set SearchInList [lsearch [lindex $Groups $i] Table_Interpolation_Tangential*]
        if {$SearchInList > -1} {
            incr TableId
            lappend AuxList $TableId
            puts $varfile "Begin Table $TableId TIME TANGENTIAL_CONTACT_STRESS"
            set Table [lindex [lindex $Groups $i] [expr { $SearchInList+1 }]]
            for {set j 2} {$j <= [lindex $Table 1]} {incr j 2} {
                puts $varfile "  [lindex $Table $j] [lindex $Table [expr { $j+1 }]]"
            }
            puts $varfile "End Table"
            puts $varfile ""
        } else {
            lappend AuxList 0
        }
        lappend TableList [lindex [lindex $Groups $i] 1] $AuxList
        set AuxList [list]
    }
    # Normal_Fluid_Flux
    set Groups [GiD_Info conditions Normal_Fluid_Flux groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
        set SearchInList [lsearch [lindex $Groups $i] Table_Interpolation*]
        if {$SearchInList > -1} {
            incr TableId
            lappend TableList [lindex [lindex $Groups $i] 1] $TableId
            puts $varfile "Begin Table $TableId TIME NORMAL_FLUID_FLUX"
            set Table [lindex [lindex $Groups $i] [expr { $SearchInList+1 }]]
            for {set j 2} {$j <= [lindex $Table 1]} {incr j 2} {
                puts $varfile "  [lindex $Table $j] [lindex $Table [expr { $j+1 }]]"
            }
            puts $varfile "End Table"
            puts $varfile ""
        } else {
            lappend TableList [lindex [lindex $Groups $i] 1] 0
        }
    }
    # Interface_Face_Load
    set Groups [GiD_Info conditions Interface_Face_Load groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
        set SearchInList [lsearch [lindex $Groups $i] Table_Interpolation_X*]
        if {$SearchInList > -1} {
            incr TableId
            lappend AuxList $TableId
            puts $varfile "Begin Table $TableId TIME FACE_LOAD_X"
            set Table [lindex [lindex $Groups $i] [expr { $SearchInList+1 }]]
            for {set j 2} {$j <= [lindex $Table 1]} {incr j 2} {
                puts $varfile "  [lindex $Table $j] [lindex $Table [expr { $j+1 }]]"
            }
            puts $varfile "End Table"
            puts $varfile ""
        } else {
            lappend AuxList 0
        }
        set SearchInList [lsearch [lindex $Groups $i] Table_Interpolation_Y*]
        if {$SearchInList > -1} {
            incr TableId
            lappend AuxList $TableId
            puts $varfile "Begin Table $TableId TIME FACE_LOAD_Y"
            set Table [lindex [lindex $Groups $i] [expr { $SearchInList+1 }]]
            for {set j 2} {$j <= [lindex $Table 1]} {incr j 2} {
                puts $varfile "  [lindex $Table $j] [lindex $Table [expr { $j+1 }]]"
            }
            puts $varfile "End Table"
            puts $varfile ""
        } else {
            lappend AuxList 0
        }
        set SearchInList [lsearch [lindex $Groups $i] Table_Interpolation_Z*]
        if {$SearchInList > -1} {
            incr TableId
            lappend AuxList $TableId
            puts $varfile "Begin Table $TableId TIME FACE_LOAD_Z"
            set Table [lindex [lindex $Groups $i] [expr { $SearchInList+1 }]]
            for {set j 2} {$j <= [lindex $Table 1]} {incr j 2} {
                puts $varfile "  [lindex $Table $j] [lindex $Table [expr { $j+1 }]]"
            }
            puts $varfile "End Table"
            puts $varfile ""
        } else {
            lappend AuxList 0
        }
        lappend TableList [lindex [lindex $Groups $i] 1] $AuxList
        set AuxList [list]
    }
    # Interface_Normal_Fluid_Flux
    set Groups [GiD_Info conditions Interface_Normal_Fluid_Flux groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
        set SearchInList [lsearch [lindex $Groups $i] Table_Interpolation*]
        if {$SearchInList > -1} {
            incr TableId
            lappend TableList [lindex [lindex $Groups $i] 1] $TableId
            puts $varfile "Begin Table $TableId TIME NORMAL_FLUID_FLUX"
            set Table [lindex [lindex $Groups $i] [expr { $SearchInList+1 }]]
            for {set j 2} {$j <= [lindex $Table 1]} {incr j 2} {
                puts $varfile "  [lindex $Table $j] [lindex $Table [expr { $j+1 }]]"
            }
            puts $varfile "End Table"
            puts $varfile ""
        } else {
            lappend TableList [lindex [lindex $Groups $i] 1] 0
        }
    }
    # Body_Acceleration
    set Groups [GiD_Info conditions Body_Acceleration groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
        set SearchInList [lsearch [lindex $Groups $i] Table_Interpolation_X*]
        if {$SearchInList > -1} {
            incr TableId
            lappend AuxList $TableId
            puts $varfile "Begin Table $TableId TIME VOLUME_ACCELERATION_X"
            set Table [lindex [lindex $Groups $i] [expr { $SearchInList+1 }]]
            for {set j 2} {$j <= [lindex $Table 1]} {incr j 2} {
                puts $varfile "  [lindex $Table $j] [lindex $Table [expr { $j+1 }]]"
            }
            puts $varfile "End Table"
            puts $varfile ""
        } else {
            lappend AuxList 0
        }
        set SearchInList [lsearch [lindex $Groups $i] Table_Interpolation_Y*]
        if {$SearchInList > -1} {
            incr TableId
            lappend AuxList $TableId
            puts $varfile "Begin Table $TableId TIME VOLUME_ACCELERATION_Y"
            set Table [lindex [lindex $Groups $i] [expr { $SearchInList+1 }]]
            for {set j 2} {$j <= [lindex $Table 1]} {incr j 2} {
                puts $varfile "  [lindex $Table $j] [lindex $Table [expr { $j+1 }]]"
            }
            puts $varfile "End Table"
            puts $varfile ""
        } else {
            lappend AuxList 0
        }
        set SearchInList [lsearch [lindex $Groups $i] Table_Interpolation_Z*]
        if {$SearchInList > -1} {
            incr TableId
            lappend AuxList $TableId
            puts $varfile "Begin Table $TableId TIME VOLUME_ACCELERATION_Z"
            set Table [lindex [lindex $Groups $i] [expr { $SearchInList+1 }]]
            for {set j 2} {$j <= [lindex $Table 1]} {incr j 2} {
                puts $varfile "  [lindex $Table $j] [lindex $Table [expr { $j+1 }]]"
            }
            puts $varfile "End Table"
            puts $varfile ""
        } else {
            lappend AuxList 0
        }
        lappend TableList [lindex [lindex $Groups $i] 1] $AuxList
        set AuxList [list]
    }
    puts $varfile ""
    
    ## Properties
    set PropertyId 0
    set PropertyList [list]
    # Body_Part
    set Groups [GiD_Info conditions Body_Part groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
        if {[lindex [lindex $Groups $i] 3]=="LinearElastic3DLaw"} {
            incr PropertyId
            lappend PropertyList [lindex [lindex $Groups $i] 1] $PropertyId
            puts $varfile "Begin Properties $PropertyId"
            puts $varfile "  CONSTITUTIVE_LAW_NAME LinearElastic3DLaw"
            puts $varfile "  YOUNG_MODULUS [lindex [lindex $Groups $i] 4]"
            puts $varfile "  POISSON_RATIO [lindex [lindex $Groups $i] 5]"
            puts $varfile "  DENSITY_SOLID [lindex [lindex $Groups $i] 6]"
            puts $varfile "  DENSITY_WATER [lindex [lindex $Groups $i] 7]"
            puts $varfile "  POROSITY [lindex [lindex $Groups $i] 8]"
            puts $varfile "  BULK_MODULUS_SOLID [lindex [lindex $Groups $i] 9]"
            puts $varfile "  BULK_MODULUS_FLUID [lindex [lindex $Groups $i] 10]"
            puts $varfile "  PERMEABILITY_XX [lindex [lindex $Groups $i] 11]"
            puts $varfile "  PERMEABILITY_YY [lindex [lindex $Groups $i] 12]"
            puts $varfile "  PERMEABILITY_ZZ [lindex [lindex $Groups $i] 13]"
            puts $varfile "  PERMEABILITY_XY [lindex [lindex $Groups $i] 14]"
            puts $varfile "  PERMEABILITY_YZ [lindex [lindex $Groups $i] 15]"
            puts $varfile "  PERMEABILITY_ZX [lindex [lindex $Groups $i] 16]"
            puts $varfile "  DYNAMIC_VISCOSITY [lindex [lindex $Groups $i] 17]"
            puts $varfile "End Properties"
            puts $varfile ""
        } elseif {[lindex [lindex $Groups $i] 3]=="LinearElasticPlaneStrain2DLaw" || [lindex [lindex $Groups $i] 3]=="LinearElasticPlaneStress2DLaw"} {
            incr PropertyId
            lappend PropertyList [lindex [lindex $Groups $i] 1] $PropertyId
            puts $varfile "Begin Properties $PropertyId"
            puts $varfile "  CONSTITUTIVE_LAW_NAME [lindex [lindex $Groups $i] 3]"
            puts $varfile "  YOUNG_MODULUS [lindex [lindex $Groups $i] 4]"
            puts $varfile "  POISSON_RATIO [lindex [lindex $Groups $i] 5]"
            puts $varfile "  DENSITY_SOLID [lindex [lindex $Groups $i] 6]"
            puts $varfile "  DENSITY_WATER [lindex [lindex $Groups $i] 7]"
            puts $varfile "  POROSITY [lindex [lindex $Groups $i] 8]"
            puts $varfile "  BULK_MODULUS_SOLID [lindex [lindex $Groups $i] 9]"
            puts $varfile "  BULK_MODULUS_FLUID [lindex [lindex $Groups $i] 10]"
            puts $varfile "  PERMEABILITY_XX [lindex [lindex $Groups $i] 11]"
            puts $varfile "  PERMEABILITY_YY [lindex [lindex $Groups $i] 12]"
            puts $varfile "  PERMEABILITY_XY [lindex [lindex $Groups $i] 14]"
            puts $varfile "  DYNAMIC_VISCOSITY [lindex [lindex $Groups $i] 17]"
            puts $varfile "  THICKNESS [lindex [lindex $Groups $i] 23]"
            puts $varfile "End Properties"
            puts $varfile ""
        } elseif {[lindex [lindex $Groups $i] 3]=="SimoJuDamage3DLaw"} {
            incr PropertyId
            lappend PropertyList [lindex [lindex $Groups $i] 1] $PropertyId
            puts $varfile "Begin Properties $PropertyId"
            if {[GiD_AccessValue get gendata Non-local_Damage]==true} {
                puts $varfile "  CONSTITUTIVE_LAW_NAME SimoJuNonlocalDamage3DLaw"
            } else {
                puts $varfile "  CONSTITUTIVE_LAW_NAME SimoJuLocalDamage3DLaw"
            }
            puts $varfile "  YOUNG_MODULUS [lindex [lindex $Groups $i] 4]"
            puts $varfile "  POISSON_RATIO [lindex [lindex $Groups $i] 5]"
            puts $varfile "  DENSITY_SOLID [lindex [lindex $Groups $i] 6]"
            puts $varfile "  DENSITY_WATER [lindex [lindex $Groups $i] 7]"
            puts $varfile "  POROSITY [lindex [lindex $Groups $i] 8]"
            puts $varfile "  BULK_MODULUS_SOLID [lindex [lindex $Groups $i] 9]"
            puts $varfile "  BULK_MODULUS_FLUID [lindex [lindex $Groups $i] 10]"
            puts $varfile "  PERMEABILITY_XX [lindex [lindex $Groups $i] 11]"
            puts $varfile "  PERMEABILITY_YY [lindex [lindex $Groups $i] 12]"
            puts $varfile "  PERMEABILITY_ZZ [lindex [lindex $Groups $i] 13]"
            puts $varfile "  PERMEABILITY_XY [lindex [lindex $Groups $i] 14]"
            puts $varfile "  PERMEABILITY_YZ [lindex [lindex $Groups $i] 15]"
            puts $varfile "  PERMEABILITY_ZX [lindex [lindex $Groups $i] 16]"
            puts $varfile "  DYNAMIC_VISCOSITY [lindex [lindex $Groups $i] 17]"
            puts $varfile "  DAMAGE_THRESHOLD [lindex [lindex $Groups $i] 18]"
            puts $varfile "  STRENGTH_RATIO [lindex [lindex $Groups $i] 19]"
            puts $varfile "  FRACTURE_ENERGY [lindex [lindex $Groups $i] 20]"
            puts $varfile "End Properties"
            puts $varfile ""
        } elseif {[lindex [lindex $Groups $i] 3]=="SimoJuDamagePlaneStrain2DLaw" || [lindex [lindex $Groups $i] 3]=="SimoJuDamagePlaneStress2DLaw"} {
            incr PropertyId
            lappend PropertyList [lindex [lindex $Groups $i] 1] $PropertyId
            puts $varfile "Begin Properties $PropertyId"
            if {[GiD_AccessValue get gendata Non-local_Damage]==true} {
                if {[lindex [lindex $Groups $i] 3]=="SimoJuDamagePlaneStrain2DLaw"} {
                    puts $varfile "  CONSTITUTIVE_LAW_NAME SimoJuNonlocalDamagePlaneStrain2DLaw"
                } else {
                    puts $varfile "  CONSTITUTIVE_LAW_NAME SimoJuNonlocalDamagePlaneStress2DLaw"
                }
            } else {
                if {[lindex [lindex $Groups $i] 3]=="SimoJuDamagePlaneStrain2DLaw"} {
                    puts $varfile "  CONSTITUTIVE_LAW_NAME SimoJuLocalDamagePlaneStrain2DLaw"
                } else {
                    puts $varfile "  CONSTITUTIVE_LAW_NAME SimoJuLocalDamagePlaneStress2DLaw"
                }
            }
            puts $varfile "  YOUNG_MODULUS [lindex [lindex $Groups $i] 4]"
            puts $varfile "  POISSON_RATIO [lindex [lindex $Groups $i] 5]"
            puts $varfile "  DENSITY_SOLID [lindex [lindex $Groups $i] 6]"
            puts $varfile "  DENSITY_WATER [lindex [lindex $Groups $i] 7]"
            puts $varfile "  POROSITY [lindex [lindex $Groups $i] 8]"
            puts $varfile "  BULK_MODULUS_SOLID [lindex [lindex $Groups $i] 9]"
            puts $varfile "  BULK_MODULUS_FLUID [lindex [lindex $Groups $i] 10]"
            puts $varfile "  PERMEABILITY_XX [lindex [lindex $Groups $i] 11]"
            puts $varfile "  PERMEABILITY_YY [lindex [lindex $Groups $i] 12]"
            puts $varfile "  PERMEABILITY_XY [lindex [lindex $Groups $i] 14]"
            puts $varfile "  DYNAMIC_VISCOSITY [lindex [lindex $Groups $i] 17]"
            puts $varfile "  DAMAGE_THRESHOLD [lindex [lindex $Groups $i] 18]"
            puts $varfile "  STRENGTH_RATIO [lindex [lindex $Groups $i] 19]"
            puts $varfile "  FRACTURE_ENERGY [lindex [lindex $Groups $i] 20]"
            puts $varfile "  THICKNESS [lindex [lindex $Groups $i] 23]"
            puts $varfile "End Properties"
            puts $varfile ""
        } elseif {[lindex [lindex $Groups $i] 3]=="ModifiedMisesDamage3DLaw"} {
            incr PropertyId
            lappend PropertyList [lindex [lindex $Groups $i] 1] $PropertyId
            puts $varfile "Begin Properties $PropertyId"
            puts $varfile "  CONSTITUTIVE_LAW_NAME ModifiedMisesNonlocalDamage3DLaw"
            puts $varfile "  YOUNG_MODULUS [lindex [lindex $Groups $i] 4]"
            puts $varfile "  POISSON_RATIO [lindex [lindex $Groups $i] 5]"
            puts $varfile "  DENSITY_SOLID [lindex [lindex $Groups $i] 6]"
            puts $varfile "  DENSITY_WATER [lindex [lindex $Groups $i] 7]"
            puts $varfile "  POROSITY [lindex [lindex $Groups $i] 8]"
            puts $varfile "  BULK_MODULUS_SOLID [lindex [lindex $Groups $i] 9]"
            puts $varfile "  BULK_MODULUS_FLUID [lindex [lindex $Groups $i] 10]"
            puts $varfile "  PERMEABILITY_XX [lindex [lindex $Groups $i] 11]"
            puts $varfile "  PERMEABILITY_YY [lindex [lindex $Groups $i] 12]"
            puts $varfile "  PERMEABILITY_ZZ [lindex [lindex $Groups $i] 13]"
            puts $varfile "  PERMEABILITY_XY [lindex [lindex $Groups $i] 14]"
            puts $varfile "  PERMEABILITY_YZ [lindex [lindex $Groups $i] 15]"
            puts $varfile "  PERMEABILITY_ZX [lindex [lindex $Groups $i] 16]"
            puts $varfile "  DYNAMIC_VISCOSITY [lindex [lindex $Groups $i] 17]"
            puts $varfile "  DAMAGE_THRESHOLD [lindex [lindex $Groups $i] 18]"
            puts $varfile "  STRENGTH_RATIO [lindex [lindex $Groups $i] 19]"
            puts $varfile "  RESIDUAL_STRENGTH [lindex [lindex $Groups $i] 21]"
            puts $varfile "  SOFTENING_SLOPE [lindex [lindex $Groups $i] 22]"
            puts $varfile "End Properties"
            puts $varfile ""
        } elseif {[lindex [lindex $Groups $i] 3]=="ModifiedMisesDamagePlaneStrain2DLaw" || [lindex [lindex $Groups $i] 3]=="ModifiedMisesDamagePlaneStress2DLaw"} {
            incr PropertyId
            lappend PropertyList [lindex [lindex $Groups $i] 1] $PropertyId
            puts $varfile "Begin Properties $PropertyId"
            if {[lindex [lindex $Groups $i] 3]=="ModifiedMisesDamagePlaneStrain2DLaw"} {
                puts $varfile "  CONSTITUTIVE_LAW_NAME ModifiedMisesNonlocalDamagePlaneStrain2DLaw"
            } else {
                puts $varfile "  CONSTITUTIVE_LAW_NAME ModifiedMisesNonlocalDamagePlaneStress2DLaw"
            }
            puts $varfile "  YOUNG_MODULUS [lindex [lindex $Groups $i] 4]"
            puts $varfile "  POISSON_RATIO [lindex [lindex $Groups $i] 5]"
            puts $varfile "  DENSITY_SOLID [lindex [lindex $Groups $i] 6]"
            puts $varfile "  DENSITY_WATER [lindex [lindex $Groups $i] 7]"
            puts $varfile "  POROSITY [lindex [lindex $Groups $i] 8]"
            puts $varfile "  BULK_MODULUS_SOLID [lindex [lindex $Groups $i] 9]"
            puts $varfile "  BULK_MODULUS_FLUID [lindex [lindex $Groups $i] 10]"
            puts $varfile "  PERMEABILITY_XX [lindex [lindex $Groups $i] 11]"
            puts $varfile "  PERMEABILITY_YY [lindex [lindex $Groups $i] 12]"
            puts $varfile "  PERMEABILITY_XY [lindex [lindex $Groups $i] 14]"
            puts $varfile "  DYNAMIC_VISCOSITY [lindex [lindex $Groups $i] 17]"
            puts $varfile "  DAMAGE_THRESHOLD [lindex [lindex $Groups $i] 18]"
            puts $varfile "  STRENGTH_RATIO [lindex [lindex $Groups $i] 19]"
            puts $varfile "  RESIDUAL_STRENGTH [lindex [lindex $Groups $i] 21]"
            puts $varfile "  SOFTENING_SLOPE [lindex [lindex $Groups $i] 22]"
            puts $varfile "  THICKNESS [lindex [lindex $Groups $i] 23]"
            puts $varfile "End Properties"
            puts $varfile ""
        }
    }
    # Interface_Part
    set Groups [GiD_Info conditions Interface_Part groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
        if {[lindex [lindex $Groups $i] 4]=="BilinearCohesive3DLaw"} {
            incr PropertyId
            lappend PropertyList [lindex [lindex $Groups $i] 1] $PropertyId
            puts $varfile "Begin Properties $PropertyId"
            puts $varfile "  CONSTITUTIVE_LAW_NAME BilinearCohesive3DLaw"
            puts $varfile "  YOUNG_MODULUS [lindex [lindex $Groups $i] 5]"
            puts $varfile "  POISSON_RATIO [lindex [lindex $Groups $i] 6]"
            puts $varfile "  DENSITY_SOLID [lindex [lindex $Groups $i] 7]"
            puts $varfile "  DENSITY_WATER [lindex [lindex $Groups $i] 8]"
            puts $varfile "  POROSITY [lindex [lindex $Groups $i] 9]"
            puts $varfile "  BULK_MODULUS_SOLID [lindex [lindex $Groups $i] 10]"
            puts $varfile "  BULK_MODULUS_FLUID [lindex [lindex $Groups $i] 11]"
            puts $varfile "  TRANSVERSAL_PERMEABILITY [lindex [lindex $Groups $i] 12]"
            puts $varfile "  DYNAMIC_VISCOSITY [lindex [lindex $Groups $i] 13]"
            puts $varfile "  DAMAGE_THRESHOLD [lindex [lindex $Groups $i] 14]"
            puts $varfile "  MINIMUM_JOINT_WIDTH [lindex [lindex $Groups $i] 15]"
            puts $varfile "  CRITICAL_DISPLACEMENT [lindex [lindex $Groups $i] 16]"
            puts $varfile "  YIELD_STRESS [lindex [lindex $Groups $i] 17]"
            puts $varfile "  FRICTION_COEFFICIENT [lindex [lindex $Groups $i] 18]"
            puts $varfile "End Properties"
            puts $varfile ""
        } elseif {[lindex [lindex $Groups $i] 4]=="BilinearCohesivePlaneStrain2DLaw" || [lindex [lindex $Groups $i] 4]=="BilinearCohesivePlaneStress2DLaw"} {
            incr PropertyId
            lappend PropertyList [lindex [lindex $Groups $i] 1] $PropertyId
            puts $varfile "Begin Properties $PropertyId"
            puts $varfile "  CONSTITUTIVE_LAW_NAME BilinearCohesive2DLaw"
            puts $varfile "  YOUNG_MODULUS [lindex [lindex $Groups $i] 5]"
            puts $varfile "  POISSON_RATIO [lindex [lindex $Groups $i] 6]"
            puts $varfile "  DENSITY_SOLID [lindex [lindex $Groups $i] 7]"
            puts $varfile "  DENSITY_WATER [lindex [lindex $Groups $i] 8]"
            puts $varfile "  POROSITY [lindex [lindex $Groups $i] 9]"
            puts $varfile "  BULK_MODULUS_SOLID [lindex [lindex $Groups $i] 10]"
            puts $varfile "  BULK_MODULUS_FLUID [lindex [lindex $Groups $i] 11]"
            puts $varfile "  TRANSVERSAL_PERMEABILITY [lindex [lindex $Groups $i] 12]"
            puts $varfile "  DYNAMIC_VISCOSITY [lindex [lindex $Groups $i] 13]"
            puts $varfile "  DAMAGE_THRESHOLD [lindex [lindex $Groups $i] 14]"
            puts $varfile "  MINIMUM_JOINT_WIDTH [lindex [lindex $Groups $i] 15]"
            puts $varfile "  CRITICAL_DISPLACEMENT [lindex [lindex $Groups $i] 16]"
            puts $varfile "  YIELD_STRESS [lindex [lindex $Groups $i] 17]"
            puts $varfile "  FRICTION_COEFFICIENT [lindex [lindex $Groups $i] 18]"
            puts $varfile "  THICKNESS [lindex [lindex $Groups $i] 19]"
            puts $varfile "End Properties"
            puts $varfile ""
        }
    }
    puts $varfile ""
    
    ## Nodes
    set Nodes [GiD_Info Mesh Nodes]
    puts $varfile "Begin Nodes"
    for {set i 0} {$i < [llength $Nodes]} {incr i 4} {
        puts $varfile "  [lindex $Nodes $i]  [lindex $Nodes [expr { $i+1 }]] [lindex $Nodes [expr { $i+2 }]] [lindex $Nodes [expr { $i+3 }]]"
    }
    puts $varfile "End Nodes"
    puts $varfile ""
    puts $varfile ""
    
    ## Elements
    set FIC [GiD_AccessValue get gendata FIC_Stabilization]
    set IsQuadratic [GiD_Info Project Quadratic]
    # Body_Part
    set Groups [GiD_Info conditions Body_Part groups]
    if {$IsQuadratic==0} {
        if {$FIC==false} {
            for {set i 0} {$i < [llength $Groups]} {incr i} {
                # Elements Material
                set ElemsMat [lindex $PropertyList [expr { [lsearch $PropertyList [lindex [lindex $Groups $i] 1]*]+1}]]
                
                # UPwSmallStrainElement2D3N
                set Entities [GiD_EntitiesGroups get [lindex [lindex $Groups $i] 1] elements -element_type triangle]
                if {[llength $Entities] > 0} {
                    puts $varfile "Begin Elements UPwSmallStrainElement2D3N"
                    for {set j 0} {$j < [llength $Entities]} {incr j} {
                        set Connectivities [Triangle2D3Connectivities [lindex $Entities $j]]
                        puts $varfile "  [lindex $Entities $j]  $ElemsMat  $Connectivities"
                    }
                    puts $varfile "End Elements"
                    puts $varfile ""
                }
                # UPwSmallStrainElement2D4N
                set Entities [GiD_EntitiesGroups get [lindex [lindex $Groups $i] 1] elements -element_type quadrilateral]
                if {[llength $Entities] > 0} {
                    puts $varfile "Begin Elements UPwSmallStrainElement2D4N"
                    for {set j 0} {$j < [llength $Entities]} {incr j} {
                        set Connectivities [Quadrilateral2D4Connectivities [lindex $Entities $j]]
                        puts $varfile "  [lindex $Entities $j]  $ElemsMat  $Connectivities"
                    }
                    puts $varfile "End Elements"
                    puts $varfile ""
                }
                # UPwSmallStrainElement3D4N
                set Entities [GiD_EntitiesGroups get [lindex [lindex $Groups $i] 1] elements -element_type tetrahedra]
                if {[llength $Entities] > 0} {
                    puts $varfile "Begin Elements UPwSmallStrainElement3D4N"
                    for {set j 0} {$j < [llength $Entities]} {incr j} {
                        set Connectivities [Quadrilateral2D4Connectivities [lindex $Entities $j]]
                        puts $varfile "  [lindex $Entities $j]  $ElemsMat  $Connectivities"
                    }
                    puts $varfile "End Elements"
                    puts $varfile ""
                }
                # UPwSmallStrainElement3D8N
                set Entities [GiD_EntitiesGroups get [lindex [lindex $Groups $i] 1] elements -element_type hexahedra]
                if {[llength $Entities] > 0} {
                    puts $varfile "Begin Elements UPwSmallStrainElement3D8N"
                    for {set j 0} {$j < [llength $Entities]} {incr j} {
                        set Connectivities [Hexahedron3D8Connectivities [lindex $Entities $j]]
                        puts $varfile "  [lindex $Entities $j]  $ElemsMat  $Connectivities"
                    }
                    puts $varfile "End Elements"
                    puts $varfile ""
                }
            }
        } else {
            for {set i 0} {$i < [llength $Groups]} {incr i} {
                # Elements Material
                set ElemsMat [lindex $PropertyList [expr { [lsearch $PropertyList [lindex [lindex $Groups $i] 1]*]+1}]]
                
                # UPwSmallStrainFICElement2D3N
                set Entities [GiD_EntitiesGroups get [lindex [lindex $Groups $i] 1] elements -element_type triangle]
                if {[llength $Entities] > 0} {
                    puts $varfile "Begin Elements UPwSmallStrainFICElement2D3N"
                    for {set j 0} {$j < [llength $Entities]} {incr j} {
                        set Connectivities [Triangle2D3Connectivities [lindex $Entities $j]]
                        puts $varfile "  [lindex $Entities $j]  $ElemsMat  $Connectivities"
                    }
                    puts $varfile "End Elements"
                    puts $varfile ""
                }
                # UPwSmallStrainFICElement2D4N
                set Entities [GiD_EntitiesGroups get [lindex [lindex $Groups $i] 1] elements -element_type quadrilateral]
                if {[llength $Entities] > 0} {
                    puts $varfile "Begin Elements UPwSmallStrainFICElement2D4N"
                    for {set j 0} {$j < [llength $Entities]} {incr j} {
                        set Connectivities [Quadrilateral2D4Connectivities [lindex $Entities $j]]
                        puts $varfile "  [lindex $Entities $j]  $ElemsMat  $Connectivities"
                    }
                    puts $varfile "End Elements"
                    puts $varfile ""
                }
                # UPwSmallStrainFICElement3D4N
                set Entities [GiD_EntitiesGroups get [lindex [lindex $Groups $i] 1] elements -element_type tetrahedra]
                if {[llength $Entities] > 0} {
                    puts $varfile "Begin Elements UPwSmallStrainFICElement3D4N"
                    for {set j 0} {$j < [llength $Entities]} {incr j} {
                        set Connectivities [Quadrilateral2D4Connectivities [lindex $Entities $j]]
                        puts $varfile "  [lindex $Entities $j]  $ElemsMat  $Connectivities"
                    }
                    puts $varfile "End Elements"
                    puts $varfile ""
                }
                # UPwSmallStrainFICElement3D8N
                set Entities [GiD_EntitiesGroups get [lindex [lindex $Groups $i] 1] elements -element_type hexahedra]
                if {[llength $Entities] > 0} {
                    puts $varfile "Begin Elements UPwSmallStrainFICElement3D8N"
                    for {set j 0} {$j < [llength $Entities]} {incr j} {
                        set Connectivities [Hexahedron3D8Connectivities [lindex $Entities $j]]
                        puts $varfile "  [lindex $Entities $j]  $ElemsMat  $Connectivities"
                    }
                    puts $varfile "End Elements"
                    puts $varfile ""
                }
            }
        }
    } elseif {$IsQuadratic==1} {
        for {set i 0} {$i < [llength $Groups]} {incr i} {
            # Elements Material
            set ElemsMat [lindex $PropertyList [expr { [lsearch $PropertyList [lindex [lindex $Groups $i] 1]*]+1}]]
            
            # SmallStrainUPwDiffOrderElement2D6N
            set Entities [GiD_EntitiesGroups get [lindex [lindex $Groups $i] 1] elements -element_type triangle]
            if {[llength $Entities] > 0} {
                puts $varfile "Begin Elements SmallStrainUPwDiffOrderElement2D6N"
                for {set j 0} {$j < [llength $Entities]} {incr j} {
                    set Connectivities [Triangle2D6Connectivities [lindex $Entities $j]]
                    puts $varfile "  [lindex $Entities $j]  $ElemsMat  $Connectivities"
                }
                puts $varfile "End Elements"
                puts $varfile ""
            }
            # SmallStrainUPwDiffOrderElement2D8N
            set Entities [GiD_EntitiesGroups get [lindex [lindex $Groups $i] 1] elements -element_type quadrilateral]
            if {[llength $Entities] > 0} {
                puts $varfile "Begin Elements SmallStrainUPwDiffOrderElement2D8N"
                
                for {set j 0} {$j < [llength $Entities]} {incr j} {
                    set Connectivities [Hexahedron3D8Connectivities [lindex $Entities $j]]
                    puts $varfile "  [lindex $Entities $j]  $ElemsMat  $Connectivities"
                }
                puts $varfile "End Elements"
                puts $varfile ""
            }
            # SmallStrainUPwDiffOrderElement3D10N
            set Entities [GiD_EntitiesGroups get [lindex [lindex $Groups $i] 1] elements -element_type tetrahedra]
            if {[llength $Entities] > 0} {
                puts $varfile "Begin Elements SmallStrainUPwDiffOrderElement3D10N"
                for {set j 0} {$j < [llength $Entities]} {incr j} {
                    set Connectivities [Tetrahedron3D10Connectivities [lindex $Entities $j]]
                    puts $varfile "  [lindex $Entities $j]  $ElemsMat  $Connectivities"
                }
                puts $varfile "End Elements"
                puts $varfile ""
            }
            # SmallStrainUPwDiffOrderElement3D20N
            set Entities [GiD_EntitiesGroups get [lindex [lindex $Groups $i] 1] elements -element_type hexahedra]
            if {[llength $Entities] > 0} {
                puts $varfile "Begin Elements SmallStrainUPwDiffOrderElement3D20N"
                for {set j 0} {$j < [llength $Entities]} {incr j} {
                    set Connectivities [Hexahedron3D20Connectivities [lindex $Entities $j]]
                    puts $varfile "  [lindex $Entities $j]  $ElemsMat  $Connectivities"
                }
                puts $varfile "End Elements"
                puts $varfile ""
            }
        }
    } else {
        for {set i 0} {$i < [llength $Groups]} {incr i} {
            # Elements Material
            set ElemsMat [lindex $PropertyList [expr { [lsearch $PropertyList [lindex [lindex $Groups $i] 1]*]+1}]]
                
            # SmallStrainUPwDiffOrderElement2D6N
            set Entities [GiD_EntitiesGroups get [lindex [lindex $Groups $i] 1] elements -element_type triangle]
            if {[llength $Entities] > 0} {
                puts $varfile "Begin Elements SmallStrainUPwDiffOrderElement2D6N"
                for {set j 0} {$j < [llength $Entities]} {incr j} {
                    set Connectivities [Triangle2D6Connectivities [lindex $Entities $j]]
                    puts $varfile "  [lindex $Entities $j]  $ElemsMat  $Connectivities"
                }
                puts $varfile "End Elements"
                puts $varfile ""
            }
            # SmallStrainUPwDiffOrderElement2D9N
            set Entities [GiD_EntitiesGroups get [lindex [lindex $Groups $i] 1] elements -element_type quadrilateral]
            if {[llength $Entities] > 0} {
                puts $varfile "Begin Elements SmallStrainUPwDiffOrderElement2D9N"
                for {set j 0} {$j < [llength $Entities]} {incr j} {
                    set Connectivities [Quadrilateral2D9Connectivities [lindex $Entities $j]]
                    puts $varfile "  [lindex $Entities $j]  $ElemsMat  $Connectivities"
                }
                puts $varfile "End Elements"
                puts $varfile ""
            }
            # SmallStrainUPwDiffOrderElement3D10N
            set Entities [GiD_EntitiesGroups get [lindex [lindex $Groups $i] 1] elements -element_type tetrahedra]
            if {[llength $Entities] > 0} {
                puts $varfile "Begin Elements SmallStrainUPwDiffOrderElement3D10N"
                for {set j 0} {$j < [llength $Entities]} {incr j} {
                    set Connectivities [Tetrahedron3D10Connectivities [lindex $Entities $j]]
                    puts $varfile "  [lindex $Entities $j]  $ElemsMat  $Connectivities"
                }
                puts $varfile "End Elements"
                puts $varfile ""
            }
            # SmallStrainUPwDiffOrderElement3D27N
            set Entities [GiD_EntitiesGroups get [lindex [lindex $Groups $i] 1] elements -element_type hexahedra]
            if {[llength $Entities] > 0} {
                puts $varfile "Begin Elements SmallStrainUPwDiffOrderElement3D27N"
                for {set j 0} {$j < [llength $Entities]} {incr j} {
                    set Connectivities [Hexahedron3D27Connectivities [lindex $Entities $j]]
                    puts $varfile "  [lindex $Entities $j]  $ElemsMat  $Connectivities"
                }
                puts $varfile "End Elements"
                puts $varfile ""
            }
        }
    }
    # Interface_Part
    set Groups [GiD_Info conditions Interface_Part groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
        # Elements Material
        set ElemsMat [lindex $PropertyList [expr { [lsearch $PropertyList [lindex [lindex $Groups $i] 1]*]+1}]]
        
        # UPwSmallStrainLinkInterfaceElement2D4N
        set Entities [GiD_EntitiesGroups get [lindex [lindex $Groups $i] 1] elements -element_type triangle]
        if {[llength $Entities] > 0} {
            puts $varfile "Begin Elements UPwSmallStrainLinkInterfaceElement2D4N"
            for {set j 0} {$j < [llength $Entities]} {incr j} {
                set Connectivities [TriangleInterface2D4Connectivities [lindex $Entities $j]]
                puts $varfile "  [lindex $Entities $j]  $ElemsMat  $Connectivities"
            }
            puts $varfile "End Elements"
            puts $varfile ""
        }
        set Entities [GiD_EntitiesGroups get [lindex [lindex $Groups $i] 1] elements -element_type quadrilateral]
        if {[llength $Entities] > 0} {
            if {[lindex [lindex $Groups $i] 3]==false} {
                # UPwSmallStrainInterfaceElement2D4N
                puts $varfile "Begin Elements UPwSmallStrainInterfaceElement2D4N"
                for {set j 0} {$j < [llength $Entities]} {incr j} {
                    set Connectivities [Quadrilateral2D4Connectivities [lindex $Entities $j]]
                    puts $varfile "  [lindex $Entities $j]  $ElemsMat  $Connectivities"
                }
                puts $varfile "End Elements"
                puts $varfile ""
            } else {
                # UPwSmallStrainLinkInterfaceElement2D4N
                puts $varfile "Begin Elements UPwSmallStrainLinkInterfaceElement2D4N"
                for {set j 0} {$j < [llength $Entities]} {incr j} {
                    set Connectivities [QuadrilateralInterface2D4Connectivities [lindex $Entities $j]]
                    puts $varfile "  [lindex $Entities $j]  $ElemsMat  $Connectivities"
                }
                puts $varfile "End Elements"
                puts $varfile ""
            }
        }
        # UPwSmallStrainLinkInterfaceElement3D6N
        set Entities [GiD_EntitiesGroups get [lindex [lindex $Groups $i] 1] elements -element_type tetrahedra]
        if {[llength $Entities] > 0} {
            puts $varfile "Begin Elements UPwSmallStrainLinkInterfaceElement3D6N"
            for {set j 0} {$j < [llength $Entities]} {incr j} {
                set Connectivities [TetrahedronInterface3D6Connectivities [lindex $Entities $j]]
                puts $varfile "  [lindex $Entities $j]  $ElemsMat  $Connectivities"
            }
            puts $varfile "End Elements"
            puts $varfile ""
        }
        set Entities [GiD_EntitiesGroups get [lindex [lindex $Groups $i] 1] elements -element_type prism]
        if {[llength $Entities] > 0} {
            if {[lindex [lindex $Groups $i] 3]==false} {
                # UPwSmallStrainInterfaceElement3D6N
                puts $varfile "Begin Elements UPwSmallStrainInterfaceElement3D6N"
                for {set j 0} {$j < [llength $Entities]} {incr j} {
                    #set Connectivities [Triangle2D6Connectivities [lindex $Entities $j]]
                    set Connectivities [PrismInterface3D6Connectivities [lindex $Entities $j]]
                    puts $varfile "  [lindex $Entities $j]  $ElemsMat  $Connectivities"
                }
                puts $varfile "End Elements"
                puts $varfile ""
            } else {
                # UPwSmallStrainLinkInterfaceElement3D6N
                puts $varfile "Begin Elements UPwSmallStrainLinkInterfaceElement3D6N"
                for {set j 0} {$j < [llength $Entities]} {incr j} {
                    set Connectivities [Triangle2D6Connectivities [lindex $Entities $j]]
                    puts $varfile "  [lindex $Entities $j]  $ElemsMat  $Connectivities"
                }
                puts $varfile "End Elements"
                puts $varfile ""
            }
        }
        set Entities [GiD_EntitiesGroups get [lindex [lindex $Groups $i] 1] elements -element_type hexahedra]
        if {[llength $Entities] > 0} {
            if {[lindex [lindex $Groups $i] 3]==false} {
                # UPwSmallStrainInterfaceElement3D8N
                puts $varfile "Begin Elements UPwSmallStrainInterfaceElement3D8N"
                for {set j 0} {$j < [llength $Entities]} {incr j} {
                    #set Connectivities [Hexahedron3D8Connectivities [lindex $Entities $j]]
                    set Connectivities [HexahedronInterface3D8Connectivities [lindex $Entities $j]]
                    puts $varfile "  [lindex $Entities $j]  $ElemsMat  $Connectivities"
                }
                puts $varfile "End Elements"
                puts $varfile ""
            } else {
                # UPwSmallStrainLinkInterfaceElement3D8N
                puts $varfile "Begin Elements UPwSmallStrainLinkInterfaceElement3D8N"
                for {set j 0} {$j < [llength $Entities]} {incr j} {
                    set Connectivities [Hexahedron3D8Connectivities [lindex $Entities $j]]
                    puts $varfile "  [lindex $Entities $j]  $ElemsMat  $Connectivities"
                }
                puts $varfile "End Elements"
                puts $varfile ""
            }
        }
        set InterfaceElemsMat $ElemsMat
    }
    puts $varfile ""
    
    ## Conditions
    set ConditionId 0
    set AuxList [list]
    set ConditionList [list]
    set Dim [GiD_AccessValue get gendata Domain_Size]
    # Force
    set Groups [GiD_Info conditions Force groups]
    if {$Dim==2} {
        # UPwForceCondition2D1N
        for {set i 0} {$i < [llength $Groups]} {incr i} {
            set Entities [GiD_EntitiesGroups get [lindex [lindex $Groups $i] 1] nodes]
            puts $varfile "Begin Conditions UPwForceCondition2D1N"
            for {set j 0} {$j < [llength $Entities]} {incr j} {
                incr ConditionId
                lappend AuxList $ConditionId
                puts $varfile "  $ConditionId  1  [lindex $Entities $j]"
            }
            puts $varfile "End Conditions"
            puts $varfile ""
            # Update ConditionList
            lappend ConditionList [lindex [lindex $Groups $i] 1] $AuxList
            set AuxList [list]
        }
    } else {
        # UPwForceCondition3D1N
        for {set i 0} {$i < [llength $Groups]} {incr i} {
            set Entities [GiD_EntitiesGroups get [lindex [lindex $Groups $i] 1] nodes]
            puts $varfile "Begin Conditions UPwForceCondition3D1N"
            for {set j 0} {$j < [llength $Entities]} {incr j} {
                incr ConditionId
                lappend AuxList $ConditionId
                puts $varfile "  $ConditionId  1  [lindex $Entities $j]"
            }
            puts $varfile "End Conditions"
            puts $varfile ""
            # Update ConditionList
            lappend ConditionList [lindex [lindex $Groups $i] 1] $AuxList
            set AuxList [list]
        }
    }
    # Face_Load
    set Groups [GiD_Info conditions Face_Load groups]
    if {$Dim==2} {
        if {$IsQuadratic==0} {
            # UPwFaceLoadCondition2D2N
            for {set i 0} {$i < [llength $Groups]} {incr i} {
                set Entities [GiD_EntitiesGroups get [lindex [lindex $Groups $i] 1] faces]
                puts $varfile "Begin Conditions UPwFaceLoadCondition2D2N"
                for {set j 0} {$j < [llength [lindex $Entities 1]]} {incr j} {
                    incr ConditionId
                    lappend AuxList $ConditionId
                    set ElementGroup [GiD_EntitiesGroups entity_groups element [lindex [lindex $Entities 0] $j]]
                    for {set k 0} {$k < [llength ElementGroup]} {incr k} {
                        set SearchInList [lsearch $PropertyList [lindex $ElementGroup $k]*]
                        if {$SearchInList > -1} {
                            set ElemsMat [lindex $PropertyList [expr { $SearchInList+1 }]]
                        }
                    }
                    set Connectivities [GiD_Mesh get element [lindex [lindex $Entities 0] $j] face [lindex [lindex $Entities 1] $j]]
                    puts $varfile "  $ConditionId  $ElemsMat  $Connectivities"
                }
                puts $varfile "End Conditions"
                puts $varfile ""
                # Update ConditionList
                lappend ConditionList [lindex [lindex $Groups $i] 1] $AuxList
                set AuxList [list]
            }
        } else {
            # LineLoadDiffOrderCondition2D3N
            for {set i 0} {$i < [llength $Groups]} {incr i} {
                set Entities [GiD_EntitiesGroups get [lindex [lindex $Groups $i] 1] faces]
                puts $varfile "Begin Conditions LineLoadDiffOrderCondition2D3N"
                for {set j 0} {$j < [llength [lindex $Entities 1]]} {incr j} {
                    incr ConditionId
                    lappend AuxList $ConditionId
                    set ElementGroup [GiD_EntitiesGroups entity_groups element [lindex [lindex $Entities 0] $j]]
                    for {set k 0} {$k < [llength ElementGroup]} {incr k} {
                        set SearchInList [lsearch $PropertyList [lindex $ElementGroup $k]*]
                        if {$SearchInList > -1} {
                            set ElemsMat [lindex $PropertyList [expr { $SearchInList+1 }]]
                        }
                    }
                    set Connectivities [GiD_Mesh get element [lindex [lindex $Entities 0] $j] face [lindex [lindex $Entities 1] $j]]
                    puts $varfile "  $ConditionId  $ElemsMat  $Connectivities"
                }
                puts $varfile "End Conditions"
                puts $varfile ""
                # Update ConditionList
                lappend ConditionList [lindex [lindex $Groups $i] 1] $AuxList
                set AuxList [list]
            }
        }
    } else {
        if {$IsQuadratic==0} {
            for {set i 0} {$i < [llength $Groups]} {incr i} {
                # UPwFaceLoadCondition3D3N
                set Entities [GiD_EntitiesGroups get [lindex [lindex $Groups $i] 1] faces -element_type tetrahedra]
                if {[llength [lindex $Entities 1]] > 0} {
                    puts $varfile "Begin Conditions UPwFaceLoadCondition3D3N"
                    for {set j 0} {$j < [llength [lindex $Entities 1]]} {incr j} {
                        incr ConditionId
                        lappend AuxList $ConditionId
                        set ElementGroup [GiD_EntitiesGroups entity_groups element [lindex [lindex $Entities 0] $j]]
                        for {set k 0} {$k < [llength ElementGroup]} {incr k} {
                            set SearchInList [lsearch $PropertyList [lindex $ElementGroup $k]*]
                            if {$SearchInList > -1} {
                                set ElemsMat [lindex $PropertyList [expr { $SearchInList+1 }]]
                            }
                        }
                        set Connectivities [GiD_Mesh get element [lindex [lindex $Entities 0] $j] face [lindex [lindex $Entities 1] $j]]
                        puts $varfile "  $ConditionId  $ElemsMat  $Connectivities"
                    }
                    puts $varfile "End Conditions"
                    puts $varfile ""
                }
                set Entities [GiD_EntitiesGroups get [lindex [lindex $Groups $i] 1] faces -element_type prism]
                if {[llength [lindex $Entities 1]] > 0} {
                    puts $varfile "Begin Conditions UPwFaceLoadCondition3D3N"
                    for {set j 0} {$j < [llength [lindex $Entities 1]]} {incr j} {
                        incr ConditionId
                        lappend AuxList $ConditionId
                        set ElementGroup [GiD_EntitiesGroups entity_groups element [lindex [lindex $Entities 0] $j]]
                        for {set k 0} {$k < [llength ElementGroup]} {incr k} {
                            set SearchInList [lsearch $PropertyList [lindex $ElementGroup $k]*]
                            if {$SearchInList > -1} {
                                set ElemsMat [lindex $PropertyList [expr { $SearchInList+1 }]]
                            }
                        }
                        set Connectivities [GiD_Mesh get element [lindex [lindex $Entities 0] $j] face [lindex [lindex $Entities 1] $j]]
                        puts $varfile "  $ConditionId  $ElemsMat  $Connectivities"
                    }
                    puts $varfile "End Conditions"
                    puts $varfile ""
                }
                # UPwFaceLoadCondition3D4N
                set Entities [GiD_EntitiesGroups get [lindex [lindex $Groups $i] 1] faces -element_type hexahedra]
                if {[llength [lindex $Entities 1]] > 0} {
                    puts $varfile "Begin Conditions UPwFaceLoadCondition3D4N"
                    for {set j 0} {$j < [llength [lindex $Entities 1]]} {incr j} {
                        incr ConditionId
                        lappend AuxList $ConditionId
                        set ElementGroup [GiD_EntitiesGroups entity_groups element [lindex [lindex $Entities 0] $j]]
                        for {set k 0} {$k < [llength ElementGroup]} {incr k} {
                            set SearchInList [lsearch $PropertyList [lindex $ElementGroup $k]*]
                            if {$SearchInList > -1} {
                                set ElemsMat [lindex $PropertyList [expr { $SearchInList+1 }]]
                            }
                        }
                        # Note: In some cases I have obtained faces with a local id = 0 !! (it should always range from 1 to the number of faces of the element)
                        set Connectivities [GiD_Mesh get element [lindex [lindex $Entities 0] $j] face [lindex [lindex $Entities 1] $j]]
                        puts $varfile "  $ConditionId  $ElemsMat  $Connectivities"
                    }
                    puts $varfile "End Conditions"
                    puts $varfile ""
                }
                # Update ConditionList
                lappend ConditionList [lindex [lindex $Groups $i] 1] $AuxList
                set AuxList [list]
            }
        } elseif {$IsQuadratic==1} {
            for {set i 0} {$i < [llength $Groups]} {incr i} {
                # SurfaceLoadDiffOrderCondition3D6N
                set Entities [GiD_EntitiesGroups get [lindex [lindex $Groups $i] 1] faces -element_type tetrahedra]
                if {[llength [lindex $Entities 1]] > 0} {
                    puts $varfile "Begin Conditions SurfaceLoadDiffOrderCondition3D6N"
                    for {set j 0} {$j < [llength [lindex $Entities 1]]} {incr j} {
                        incr ConditionId
                        lappend AuxList $ConditionId
                        set ElementGroup [GiD_EntitiesGroups entity_groups element [lindex [lindex $Entities 0] $j]]
                        for {set k 0} {$k < [llength ElementGroup]} {incr k} {
                            set SearchInList [lsearch $PropertyList [lindex $ElementGroup $k]*]
                            if {$SearchInList > -1} {
                                set ElemsMat [lindex $PropertyList [expr { $SearchInList+1 }]]
                            }
                        }
                        set Connectivities [GiD_Mesh get element [lindex [lindex $Entities 0] $j] face [lindex [lindex $Entities 1] $j]]
                        puts $varfile "  $ConditionId  $ElemsMat  $Connectivities"
                    }
                    puts $varfile "End Conditions"
                    puts $varfile ""
                }
                # SurfaceLoadDiffOrderCondition3D8N
                set Entities [GiD_EntitiesGroups get [lindex [lindex $Groups $i] 1] faces -element_type hexahedra]
                if {[llength [lindex $Entities 1]] > 0} {
                    puts $varfile "Begin Conditions SurfaceLoadDiffOrderCondition3D8N"
                    for {set j 0} {$j < [llength [lindex $Entities 1]]} {incr j} {
                        incr ConditionId
                        lappend AuxList $ConditionId
                        set ElementGroup [GiD_EntitiesGroups entity_groups element [lindex [lindex $Entities 0] $j]]
                        for {set k 0} {$k < [llength ElementGroup]} {incr k} {
                            set SearchInList [lsearch $PropertyList [lindex $ElementGroup $k]*]
                            if {$SearchInList > -1} {
                                set ElemsMat [lindex $PropertyList [expr { $SearchInList+1 }]]
                            }
                        }
                        set Connectivities [GiD_Mesh get element [lindex [lindex $Entities 0] $j] face [lindex [lindex $Entities 1] $j]]
                        puts $varfile "  $ConditionId  $ElemsMat  $Connectivities"
                    }
                    puts $varfile "End Conditions"
                    puts $varfile ""
                }
                # Update ConditionList
                lappend ConditionList [lindex [lindex $Groups $i] 1] $AuxList
                set AuxList [list]
            }
        } else {
            for {set i 0} {$i < [llength $Groups]} {incr i} {
                # SurfaceLoadDiffOrderCondition3D6N
                set Entities [GiD_EntitiesGroups get [lindex [lindex $Groups $i] 1] faces -element_type tetrahedra]
                if {[llength [lindex $Entities 1]] > 0} {
                    puts $varfile "Begin Conditions SurfaceLoadDiffOrderCondition3D6N"
                    for {set j 0} {$j < [llength [lindex $Entities 1]]} {incr j} {
                        incr ConditionId
                        lappend AuxList $ConditionId
                        set ElementGroup [GiD_EntitiesGroups entity_groups element [lindex [lindex $Entities 0] $j]]
                        for {set k 0} {$k < [llength ElementGroup]} {incr k} {
                            set SearchInList [lsearch $PropertyList [lindex $ElementGroup $k]*]
                            if {$SearchInList > -1} {
                                set ElemsMat [lindex $PropertyList [expr { $SearchInList+1 }]]
                            }
                        }
                        set Connectivities [GiD_Mesh get element [lindex [lindex $Entities 0] $j] face [lindex [lindex $Entities 1] $j]]
                        puts $varfile "  $ConditionId  $ElemsMat  $Connectivities"
                    }
                    puts $varfile "End Conditions"
                    puts $varfile ""
                }
                # SurfaceLoadDiffOrderCondition3D9N
                set Entities [GiD_EntitiesGroups get [lindex [lindex $Groups $i] 1] faces -element_type hexahedra]
                if {[llength [lindex $Entities 1]] > 0} {
                    puts $varfile "Begin Conditions SurfaceLoadDiffOrderCondition3D9N"
                    for {set j 0} {$j < [llength [lindex $Entities 1]]} {incr j} {
                        incr ConditionId
                        lappend AuxList $ConditionId
                        set ElementGroup [GiD_EntitiesGroups entity_groups element [lindex [lindex $Entities 0] $j]]
                        for {set k 0} {$k < [llength ElementGroup]} {incr k} {
                            set SearchInList [lsearch $PropertyList [lindex $ElementGroup $k]*]
                            if {$SearchInList > -1} {
                                set ElemsMat [lindex $PropertyList [expr { $SearchInList+1 }]]
                            }
                        }
                        set Connectivities [GiD_Mesh get element [lindex [lindex $Entities 0] $j] face [lindex [lindex $Entities 1] $j]]
                        puts $varfile "  $ConditionId  $ElemsMat  $Connectivities"
                    }
                    puts $varfile "End Conditions"
                    puts $varfile ""
                }
                # Update ConditionList
                lappend ConditionList [lindex [lindex $Groups $i] 1] $AuxList
                set AuxList [list]
            }
        }
    }
    # Normal_Load
    set Groups [GiD_Info conditions Normal_Load groups]
    if {$Dim==2} {
        if {$IsQuadratic==0} {
            # UPwNormalFaceLoadCondition2D2N
            for {set i 0} {$i < [llength $Groups]} {incr i} {
                set Entities [GiD_EntitiesGroups get [lindex [lindex $Groups $i] 1] faces]
                puts $varfile "Begin Conditions UPwNormalFaceLoadCondition2D2N"
                for {set j 0} {$j < [llength [lindex $Entities 1]]} {incr j} {
                    incr ConditionId
                    lappend AuxList $ConditionId
                    set ElementGroup [GiD_EntitiesGroups entity_groups element [lindex [lindex $Entities 0] $j]]
                    for {set k 0} {$k < [llength ElementGroup]} {incr k} {
                        set SearchInList [lsearch $PropertyList [lindex $ElementGroup $k]*]
                        if {$SearchInList > -1} {
                            set ElemsMat [lindex $PropertyList [expr { $SearchInList+1 }]]
                        }
                    }
                    set Connectivities [GiD_Mesh get element [lindex [lindex $Entities 0] $j] face [lindex [lindex $Entities 1] $j]]
                    puts $varfile "  $ConditionId  $ElemsMat  $Connectivities"
                }
                puts $varfile "End Conditions"
                puts $varfile ""
                # Update ConditionList
                lappend ConditionList [lindex [lindex $Groups $i] 1] $AuxList
                set AuxList [list]
            }
        } else {
            # LineNormalLoadDiffOrderCondition2D3N
            for {set i 0} {$i < [llength $Groups]} {incr i} {
                set Entities [GiD_EntitiesGroups get [lindex [lindex $Groups $i] 1] faces]
                puts $varfile "Begin Conditions LineNormalLoadDiffOrderCondition2D3N"
                for {set j 0} {$j < [llength [lindex $Entities 1]]} {incr j} {
                    incr ConditionId
                    lappend AuxList $ConditionId
                    set ElementGroup [GiD_EntitiesGroups entity_groups element [lindex [lindex $Entities 0] $j]]
                    for {set k 0} {$k < [llength ElementGroup]} {incr k} {
                        set SearchInList [lsearch $PropertyList [lindex $ElementGroup $k]*]
                        if {$SearchInList > -1} {
                            set ElemsMat [lindex $PropertyList [expr { $SearchInList+1 }]]
                        }
                    }
                    set Connectivities [GiD_Mesh get element [lindex [lindex $Entities 0] $j] face [lindex [lindex $Entities 1] $j]]
                    puts $varfile "  $ConditionId  $ElemsMat  $Connectivities"
                }
                puts $varfile "End Conditions"
                puts $varfile ""
                # Update ConditionList
                lappend ConditionList [lindex [lindex $Groups $i] 1] $AuxList
                set AuxList [list]
            }
        }
    } else {
        if {$IsQuadratic==0} {
            for {set i 0} {$i < [llength $Groups]} {incr i} {
                # UPwNormalFaceLoadCondition3D3N
                set Entities [GiD_EntitiesGroups get [lindex [lindex $Groups $i] 1] faces -element_type tetrahedra]
                if {[llength [lindex $Entities 1]] > 0} {
                    puts $varfile "Begin Conditions UPwNormalFaceLoadCondition3D3N"
                    for {set j 0} {$j < [llength [lindex $Entities 1]]} {incr j} {
                        incr ConditionId
                        lappend AuxList $ConditionId
                        set ElementGroup [GiD_EntitiesGroups entity_groups element [lindex [lindex $Entities 0] $j]]
                        for {set k 0} {$k < [llength ElementGroup]} {incr k} {
                            set SearchInList [lsearch $PropertyList [lindex $ElementGroup $k]*]
                            if {$SearchInList > -1} {
                                set ElemsMat [lindex $PropertyList [expr { $SearchInList+1 }]]
                            }
                        }
                        set Connectivities [GiD_Mesh get element [lindex [lindex $Entities 0] $j] face [lindex [lindex $Entities 1] $j]]
                        puts $varfile "  $ConditionId  $ElemsMat  $Connectivities"
                    }
                    puts $varfile "End Conditions"
                    puts $varfile ""
                }
                set Entities [GiD_EntitiesGroups get [lindex [lindex $Groups $i] 1] faces -element_type prism]
                if {[llength [lindex $Entities 1]] > 0} {
                    puts $varfile "Begin Conditions UPwNormalFaceLoadCondition3D3N"
                    for {set j 0} {$j < [llength [lindex $Entities 1]]} {incr j} {
                        incr ConditionId
                        lappend AuxList $ConditionId
                        set ElementGroup [GiD_EntitiesGroups entity_groups element [lindex [lindex $Entities 0] $j]]
                        for {set k 0} {$k < [llength ElementGroup]} {incr k} {
                            set SearchInList [lsearch $PropertyList [lindex $ElementGroup $k]*]
                            if {$SearchInList > -1} {
                                set ElemsMat [lindex $PropertyList [expr { $SearchInList+1 }]]
                            }
                        }
                        set Connectivities [GiD_Mesh get element [lindex [lindex $Entities 0] $j] face [lindex [lindex $Entities 1] $j]]
                        puts $varfile "  $ConditionId  $ElemsMat  $Connectivities"
                    }
                    puts $varfile "End Conditions"
                    puts $varfile ""
                }
                # UpwNormalFaceLoadCondition3D4N
                set Entities [GiD_EntitiesGroups get [lindex [lindex $Groups $i] 1] faces -element_type hexahedra]
                if {[llength [lindex $Entities 1]] > 0} {
                    puts $varfile "Begin Conditions UpwNormalFaceLoadCondition3D4N"
                    for {set j 0} {$j < [llength [lindex $Entities 1]]} {incr j} {
                        incr ConditionId
                        lappend AuxList $ConditionId
                        set ElementGroup [GiD_EntitiesGroups entity_groups element [lindex [lindex $Entities 0] $j]]
                        for {set k 0} {$k < [llength ElementGroup]} {incr k} {
                            set SearchInList [lsearch $PropertyList [lindex $ElementGroup $k]*]
                            if {$SearchInList > -1} {
                                set ElemsMat [lindex $PropertyList [expr { $SearchInList+1 }]]
                            }
                        }
                        set Connectivities [GiD_Mesh get element [lindex [lindex $Entities 0] $j] face [lindex [lindex $Entities 1] $j]]
                        puts $varfile "  $ConditionId  $ElemsMat  $Connectivities"
                    }
                    puts $varfile "End Conditions"
                    puts $varfile ""
                }
                # Update ConditionList
                lappend ConditionList [lindex [lindex $Groups $i] 1] $AuxList
                set AuxList [list]
            }
        } elseif {$IsQuadratic==1} {
            for {set i 0} {$i < [llength $Groups]} {incr i} {
                # SurfaceNormalLoadDiffOrderCondition3D6N
                set Entities [GiD_EntitiesGroups get [lindex [lindex $Groups $i] 1] faces -element_type tetrahedra]
                if {[llength [lindex $Entities 1]] > 0} {
                    puts $varfile "Begin Conditions SurfaceNormalLoadDiffOrderCondition3D6N"
                    for {set j 0} {$j < [llength [lindex $Entities 1]]} {incr j} {
                        incr ConditionId
                        lappend AuxList $ConditionId
                        set ElementGroup [GiD_EntitiesGroups entity_groups element [lindex [lindex $Entities 0] $j]]
                        for {set k 0} {$k < [llength ElementGroup]} {incr k} {
                            set SearchInList [lsearch $PropertyList [lindex $ElementGroup $k]*]
                            if {$SearchInList > -1} {
                                set ElemsMat [lindex $PropertyList [expr { $SearchInList+1 }]]
                            }
                        }
                        set Connectivities [GiD_Mesh get element [lindex [lindex $Entities 0] $j] face [lindex [lindex $Entities 1] $j]]
                        puts $varfile "  $ConditionId  $ElemsMat  $Connectivities"
                    }
                    puts $varfile "End Conditions"
                    puts $varfile ""
                }
                # SurfaceNormalLoadDiffOrderCondition3D8N
                set Entities [GiD_EntitiesGroups get [lindex [lindex $Groups $i] 1] faces -element_type hexahedra]
                if {[llength [lindex $Entities 1]] > 0} {
                    puts $varfile "Begin Conditions SurfaceNormalLoadDiffOrderCondition3D8N"
                    for {set j 0} {$j < [llength [lindex $Entities 1]]} {incr j} {
                        incr ConditionId
                        lappend AuxList $ConditionId
                        set ElementGroup [GiD_EntitiesGroups entity_groups element [lindex [lindex $Entities 0] $j]]
                        for {set k 0} {$k < [llength ElementGroup]} {incr k} {
                            set SearchInList [lsearch $PropertyList [lindex $ElementGroup $k]*]
                            if {$SearchInList > -1} {
                                set ElemsMat [lindex $PropertyList [expr { $SearchInList+1 }]]
                            }
                        }
                        set Connectivities [GiD_Mesh get element [lindex [lindex $Entities 0] $j] face [lindex [lindex $Entities 1] $j]]
                        puts $varfile "  $ConditionId  $ElemsMat  $Connectivities"
                    }
                    puts $varfile "End Conditions"
                    puts $varfile ""
                }
                # Update ConditionList
                lappend ConditionList [lindex [lindex $Groups $i] 1] $AuxList
                set AuxList [list]
            }
        } else {
            for {set i 0} {$i < [llength $Groups]} {incr i} {
                # SurfaceNormalLoadDiffOrderCondition3D6N
                set Entities [GiD_EntitiesGroups get [lindex [lindex $Groups $i] 1] faces -element_type tetrahedra]
                if {[llength [lindex $Entities 1]] > 0} {
                    puts $varfile "Begin Conditions SurfaceNormalLoadDiffOrderCondition3D6N"
                    for {set j 0} {$j < [llength [lindex $Entities 1]]} {incr j} {
                        incr ConditionId
                        lappend AuxList $ConditionId
                        set ElementGroup [GiD_EntitiesGroups entity_groups element [lindex [lindex $Entities 0] $j]]
                        for {set k 0} {$k < [llength ElementGroup]} {incr k} {
                            set SearchInList [lsearch $PropertyList [lindex $ElementGroup $k]*]
                            if {$SearchInList > -1} {
                                set ElemsMat [lindex $PropertyList [expr { $SearchInList+1 }]]
                            }
                        }
                        set Connectivities [GiD_Mesh get element [lindex [lindex $Entities 0] $j] face [lindex [lindex $Entities 1] $j]]
                        puts $varfile "  $ConditionId  $ElemsMat  $Connectivities"
                    }
                    puts $varfile "End Conditions"
                    puts $varfile ""
                }
                # SurfaceNormalLoadDiffOrderCondition3D9N
                set Entities [GiD_EntitiesGroups get [lindex [lindex $Groups $i] 1] faces -element_type hexahedra]
                if {[llength [lindex $Entities 1]] > 0} {
                    puts $varfile "Begin Conditions SurfaceNormalLoadDiffOrderCondition3D9N"
                    for {set j 0} {$j < [llength [lindex $Entities 1]]} {incr j} {
                        incr ConditionId
                        lappend AuxList $ConditionId
                        set ElementGroup [GiD_EntitiesGroups entity_groups element [lindex [lindex $Entities 0] $j]]
                        for {set k 0} {$k < [llength ElementGroup]} {incr k} {
                            set SearchInList [lsearch $PropertyList [lindex $ElementGroup $k]*]
                            if {$SearchInList > -1} {
                                set ElemsMat [lindex $PropertyList [expr { $SearchInList+1 }]]
                            }
                        }
                        set Connectivities [GiD_Mesh get element [lindex [lindex $Entities 0] $j] face [lindex [lindex $Entities 1] $j]]
                        puts $varfile "  $ConditionId  $ElemsMat  $Connectivities"
                    }
                    puts $varfile "End Conditions"
                    puts $varfile ""
                }
                # Update ConditionList
                lappend ConditionList [lindex [lindex $Groups $i] 1] $AuxList
                set AuxList [list]
            }
        }
    }
    # Normal_Fluid_Flux
    set Groups [GiD_Info conditions Normal_Fluid_Flux groups]
    if {$Dim==2} {
        if {$IsQuadratic==0} {
            if {$FIC==false} {
                # UPwNormalFluxCondition2D2N
                for {set i 0} {$i < [llength $Groups]} {incr i} {
                    set Entities [GiD_EntitiesGroups get [lindex [lindex $Groups $i] 1] faces]
                    puts $varfile "Begin Conditions UPwNormalFluxCondition2D2N"
                    for {set j 0} {$j < [llength [lindex $Entities 1]]} {incr j} {
                        incr ConditionId
                        lappend AuxList $ConditionId
                        set ElementGroup [GiD_EntitiesGroups entity_groups element [lindex [lindex $Entities 0] $j]]
                        for {set k 0} {$k < [llength ElementGroup]} {incr k} {
                            set SearchInList [lsearch $PropertyList [lindex $ElementGroup $k]*]
                            if {$SearchInList > -1} {
                                set ElemsMat [lindex $PropertyList [expr { $SearchInList+1 }]]
                            }
                        }
                        set Connectivities [GiD_Mesh get element [lindex [lindex $Entities 0] $j] face [lindex [lindex $Entities 1] $j]]
                        puts $varfile "  $ConditionId  $ElemsMat  $Connectivities"
                    }
                    puts $varfile "End Conditions"
                    puts $varfile ""
                    # Update ConditionList
                    lappend ConditionList [lindex [lindex $Groups $i] 1] $AuxList
                    set AuxList [list]
                }
            } else {
                # UPwNormalFluxFICCondition2D2N
                for {set i 0} {$i < [llength $Groups]} {incr i} {
                    set Entities [GiD_EntitiesGroups get [lindex [lindex $Groups $i] 1] faces]
                    puts $varfile "Begin Conditions UPwNormalFluxFICCondition2D2N"
                    for {set j 0} {$j < [llength [lindex $Entities 1]]} {incr j} {
                        incr ConditionId
                        lappend AuxList $ConditionId
                        set ElementGroup [GiD_EntitiesGroups entity_groups element [lindex [lindex $Entities 0] $j]]
                        for {set k 0} {$k < [llength ElementGroup]} {incr k} {
                            set SearchInList [lsearch $PropertyList [lindex $ElementGroup $k]*]
                            if {$SearchInList > -1} {
                                set ElemsMat [lindex $PropertyList [expr { $SearchInList+1 }]]
                            }
                        }
                        set Connectivities [GiD_Mesh get element [lindex [lindex $Entities 0] $j] face [lindex [lindex $Entities 1] $j]]
                        puts $varfile "  $ConditionId  $ElemsMat  $Connectivities"
                    }
                    puts $varfile "End Conditions"
                    puts $varfile ""
                    # Update ConditionList
                    lappend ConditionList [lindex [lindex $Groups $i] 1] $AuxList
                    set AuxList [list]
                }
            }
        } else {
            # LineNormalFluidFluxDiffOrderCondition2D3N
            for {set i 0} {$i < [llength $Groups]} {incr i} {
                set Entities [GiD_EntitiesGroups get [lindex [lindex $Groups $i] 1] faces]
                puts $varfile "Begin Conditions LineNormalFluidFluxDiffOrderCondition2D3N"
                for {set j 0} {$j < [llength [lindex $Entities 1]]} {incr j} {
                    incr ConditionId
                    lappend AuxList $ConditionId
                    set ElementGroup [GiD_EntitiesGroups entity_groups element [lindex [lindex $Entities 0] $j]]
                    for {set k 0} {$k < [llength ElementGroup]} {incr k} {
                        set SearchInList [lsearch $PropertyList [lindex $ElementGroup $k]*]
                        if {$SearchInList > -1} {
                            set ElemsMat [lindex $PropertyList [expr { $SearchInList+1 }]]
                        }
                    }
                    set Connectivities [GiD_Mesh get element [lindex [lindex $Entities 0] $j] face [lindex [lindex $Entities 1] $j]]
                    puts $varfile "  $ConditionId  $ElemsMat  $Connectivities"
                }
                puts $varfile "End Conditions"
                puts $varfile ""
                # Update ConditionList
                lappend ConditionList [lindex [lindex $Groups $i] 1] $AuxList
                set AuxList [list]
            }
        }
    } else {
        if {$IsQuadratic==0} {
            if {$FIC==false} {
                for {set i 0} {$i < [llength $Groups]} {incr i} {
                    # UPwNormalFluxCondition3D3N
                    set Entities [GiD_EntitiesGroups get [lindex [lindex $Groups $i] 1] faces -element_type tetrahedra]
                    if {[llength [lindex $Entities 1]] > 0} {
                        puts $varfile "Begin Conditions UPwNormalFluxCondition3D3N"
                        for {set j 0} {$j < [llength [lindex $Entities 1]]} {incr j} {
                            incr ConditionId
                            lappend AuxList $ConditionId
                            set ElementGroup [GiD_EntitiesGroups entity_groups element [lindex [lindex $Entities 0] $j]]
                            for {set k 0} {$k < [llength ElementGroup]} {incr k} {
                                set SearchInList [lsearch $PropertyList [lindex $ElementGroup $k]*]
                                if {$SearchInList > -1} {
                                    set ElemsMat [lindex $PropertyList [expr { $SearchInList+1 }]]
                                }
                            }
                            set Connectivities [GiD_Mesh get element [lindex [lindex $Entities 0] $j] face [lindex [lindex $Entities 1] $j]]
                            puts $varfile "  $ConditionId  $ElemsMat  $Connectivities"
                        }
                        puts $varfile "End Conditions"
                        puts $varfile ""
                    }
                    set Entities [GiD_EntitiesGroups get [lindex [lindex $Groups $i] 1] faces -element_type prism]
                    if {[llength [lindex $Entities 1]] > 0} {
                        puts $varfile "Begin Conditions UPwNormalFluxCondition3D3N"
                        for {set j 0} {$j < [llength [lindex $Entities 1]]} {incr j} {
                            incr ConditionId
                            lappend AuxList $ConditionId
                            set ElementGroup [GiD_EntitiesGroups entity_groups element [lindex [lindex $Entities 0] $j]]
                            for {set k 0} {$k < [llength ElementGroup]} {incr k} {
                                set SearchInList [lsearch $PropertyList [lindex $ElementGroup $k]*]
                                if {$SearchInList > -1} {
                                    set ElemsMat [lindex $PropertyList [expr { $SearchInList+1 }]]
                                }
                            }
                            set Connectivities [GiD_Mesh get element [lindex [lindex $Entities 0] $j] face [lindex [lindex $Entities 1] $j]]
                            puts $varfile "  $ConditionId  $ElemsMat  $Connectivities"
                        }
                        puts $varfile "End Conditions"
                        puts $varfile ""
                    }
                    # UPwNormalFluxCondition3D4N
                    set Entities [GiD_EntitiesGroups get [lindex [lindex $Groups $i] 1] faces -element_type hexahedra]
                    if {[llength [lindex $Entities 1]] > 0} {
                        puts $varfile "Begin Conditions UPwNormalFluxCondition3D4N"
                        for {set j 0} {$j < [llength [lindex $Entities 1]]} {incr j} {
                            incr ConditionId
                            lappend AuxList $ConditionId
                            set ElementGroup [GiD_EntitiesGroups entity_groups element [lindex [lindex $Entities 0] $j]]
                            for {set k 0} {$k < [llength ElementGroup]} {incr k} {
                                set SearchInList [lsearch $PropertyList [lindex $ElementGroup $k]*]
                                if {$SearchInList > -1} {
                                    set ElemsMat [lindex $PropertyList [expr { $SearchInList+1 }]]
                                }
                            }
                            set Connectivities [GiD_Mesh get element [lindex [lindex $Entities 0] $j] face [lindex [lindex $Entities 1] $j]]
                            puts $varfile "  $ConditionId  $ElemsMat  $Connectivities"
                        }
                        puts $varfile "End Conditions"
                        puts $varfile ""
                    }
                    # Update ConditionList
                    lappend ConditionList [lindex [lindex $Groups $i] 1] $AuxList
                    set AuxList [list]
                }
            } else {
                for {set i 0} {$i < [llength $Groups]} {incr i} {
                    # UPwNormalFluxFICCondition3D3N
                    set Entities [GiD_EntitiesGroups get [lindex [lindex $Groups $i] 1] faces -element_type tetrahedra]
                    if {[llength [lindex $Entities 1]] > 0} {
                        puts $varfile "Begin Conditions UPwNormalFluxFICCondition3D3N"
                        for {set j 0} {$j < [llength [lindex $Entities 1]]} {incr j} {
                            incr ConditionId
                            lappend AuxList $ConditionId
                            set ElementGroup [GiD_EntitiesGroups entity_groups element [lindex [lindex $Entities 0] $j]]
                            for {set k 0} {$k < [llength ElementGroup]} {incr k} {
                                set SearchInList [lsearch $PropertyList [lindex $ElementGroup $k]*]
                                if {$SearchInList > -1} {
                                    set ElemsMat [lindex $PropertyList [expr { $SearchInList+1 }]]
                                }
                            }
                            set Connectivities [GiD_Mesh get element [lindex [lindex $Entities 0] $j] face [lindex [lindex $Entities 1] $j]]
                            puts $varfile "  $ConditionId  $ElemsMat  $Connectivities"
                        }
                        puts $varfile "End Conditions"
                        puts $varfile ""
                    }
                    set Entities [GiD_EntitiesGroups get [lindex [lindex $Groups $i] 1] faces -element_type prism]
                    if {[llength [lindex $Entities 1]] > 0} {
                        puts $varfile "Begin Conditions UPwNormalFluxFICCondition3D3N"
                        for {set j 0} {$j < [llength [lindex $Entities 1]]} {incr j} {
                            incr ConditionId
                            lappend AuxList $ConditionId
                            set ElementGroup [GiD_EntitiesGroups entity_groups element [lindex [lindex $Entities 0] $j]]
                            for {set k 0} {$k < [llength ElementGroup]} {incr k} {
                                set SearchInList [lsearch $PropertyList [lindex $ElementGroup $k]*]
                                if {$SearchInList > -1} {
                                    set ElemsMat [lindex $PropertyList [expr { $SearchInList+1 }]]
                                }
                            }
                            set Connectivities [GiD_Mesh get element [lindex [lindex $Entities 0] $j] face [lindex [lindex $Entities 1] $j]]
                            puts $varfile "  $ConditionId  $ElemsMat  $Connectivities"
                        }
                        puts $varfile "End Conditions"
                        puts $varfile ""
                    }
                    # UPwNormalFluxFICCondition3D4N
                    set Entities [GiD_EntitiesGroups get [lindex [lindex $Groups $i] 1] faces -element_type hexahedra]
                    if {[llength [lindex $Entities 1]] > 0} {
                        puts $varfile "Begin Conditions UPwNormalFluxFICCondition3D4N"
                        for {set j 0} {$j < [llength [lindex $Entities 1]]} {incr j} {
                            incr ConditionId
                            lappend AuxList $ConditionId
                            set ElementGroup [GiD_EntitiesGroups entity_groups element [lindex [lindex $Entities 0] $j]]
                            for {set k 0} {$k < [llength ElementGroup]} {incr k} {
                                set SearchInList [lsearch $PropertyList [lindex $ElementGroup $k]*]
                                if {$SearchInList > -1} {
                                    set ElemsMat [lindex $PropertyList [expr { $SearchInList+1 }]]
                                }
                            }
                            set Connectivities [GiD_Mesh get element [lindex [lindex $Entities 0] $j] face [lindex [lindex $Entities 1] $j]]
                            puts $varfile "  $ConditionId  $ElemsMat  $Connectivities"
                        }
                        puts $varfile "End Conditions"
                        puts $varfile ""
                    }
                    # Update ConditionList
                    lappend ConditionList [lindex [lindex $Groups $i] 1] $AuxList
                    set AuxList [list]
                }
            }
        } elseif {$IsQuadratic==1} {
            for {set i 0} {$i < [llength $Groups]} {incr i} {
                # SurfaceNormalFluidFluxDiffOrderCondition3D6N
                set Entities [GiD_EntitiesGroups get [lindex [lindex $Groups $i] 1] faces -element_type tetrahedra]
                if {[llength [lindex $Entities 1]] > 0} {
                    puts $varfile "Begin Conditions SurfaceNormalFluidFluxDiffOrderCondition3D6N"
                    for {set j 0} {$j < [llength [lindex $Entities 1]]} {incr j} {
                        incr ConditionId
                        lappend AuxList $ConditionId
                        set ElementGroup [GiD_EntitiesGroups entity_groups element [lindex [lindex $Entities 0] $j]]
                        for {set k 0} {$k < [llength ElementGroup]} {incr k} {
                            set SearchInList [lsearch $PropertyList [lindex $ElementGroup $k]*]
                            if {$SearchInList > -1} {
                                set ElemsMat [lindex $PropertyList [expr { $SearchInList+1 }]]
                            }
                        }
                        set Connectivities [GiD_Mesh get element [lindex [lindex $Entities 0] $j] face [lindex [lindex $Entities 1] $j]]
                        puts $varfile "  $ConditionId  $ElemsMat  $Connectivities"
                    }
                    puts $varfile "End Conditions"
                    puts $varfile ""
                }
                # SurfaceNormalFluidFluxDiffOrderCondition3D8N
                set Entities [GiD_EntitiesGroups get [lindex [lindex $Groups $i] 1] faces -element_type hexahedra]
                if {[llength [lindex $Entities 1]] > 0} {
                    puts $varfile "Begin Conditions SurfaceNormalFluidFluxDiffOrderCondition3D8N"
                    for {set j 0} {$j < [llength [lindex $Entities 1]]} {incr j} {
                        incr ConditionId
                        lappend AuxList $ConditionId
                        set ElementGroup [GiD_EntitiesGroups entity_groups element [lindex [lindex $Entities 0] $j]]
                        for {set k 0} {$k < [llength ElementGroup]} {incr k} {
                            set SearchInList [lsearch $PropertyList [lindex $ElementGroup $k]*]
                            if {$SearchInList > -1} {
                                set ElemsMat [lindex $PropertyList [expr { $SearchInList+1 }]]
                            }
                        }
                        set Connectivities [GiD_Mesh get element [lindex [lindex $Entities 0] $j] face [lindex [lindex $Entities 1] $j]]
                        puts $varfile "  $ConditionId  $ElemsMat  $Connectivities"
                    }
                    puts $varfile "End Conditions"
                    puts $varfile ""
                }
                # Update ConditionList
                lappend ConditionList [lindex [lindex $Groups $i] 1] $AuxList
                set AuxList [list]
            }
        } else {
            for {set i 0} {$i < [llength $Groups]} {incr i} {
                # SurfaceNormalFluidFluxDiffOrderCondition3D6N
                set Entities [GiD_EntitiesGroups get [lindex [lindex $Groups $i] 1] faces -element_type tetrahedra]
                if {[llength [lindex $Entities 1]] > 0} {
                    puts $varfile "Begin Conditions SurfaceNormalFluidFluxDiffOrderCondition3D6N"
                    for {set j 0} {$j < [llength [lindex $Entities 1]]} {incr j} {
                        incr ConditionId
                        lappend AuxList $ConditionId
                        set ElementGroup [GiD_EntitiesGroups entity_groups element [lindex [lindex $Entities 0] $j]]
                        for {set k 0} {$k < [llength ElementGroup]} {incr k} {
                            set SearchInList [lsearch $PropertyList [lindex $ElementGroup $k]*]
                            if {$SearchInList > -1} {
                                set ElemsMat [lindex $PropertyList [expr { $SearchInList+1 }]]
                            }
                        }
                        set Connectivities [GiD_Mesh get element [lindex [lindex $Entities 0] $j] face [lindex [lindex $Entities 1] $j]]
                        puts $varfile "  $ConditionId  $ElemsMat  $Connectivities"
                    }
                    puts $varfile "End Conditions"
                    puts $varfile ""
                }
                # SurfaceNormalFluidFluxDiffOrderCondition3D9N
                set Entities [GiD_EntitiesGroups get [lindex [lindex $Groups $i] 1] faces -element_type hexahedra]
                if {[llength [lindex $Entities 1]] > 0} {
                    puts $varfile "Begin Conditions SurfaceNormalFluidFluxDiffOrderCondition3D9N"
                    for {set j 0} {$j < [llength [lindex $Entities 1]]} {incr j} {
                        incr ConditionId
                        lappend AuxList $ConditionId
                        set ElementGroup [GiD_EntitiesGroups entity_groups element [lindex [lindex $Entities 0] $j]]
                        for {set k 0} {$k < [llength ElementGroup]} {incr k} {
                            set SearchInList [lsearch $PropertyList [lindex $ElementGroup $k]*]
                            if {$SearchInList > -1} {
                                set ElemsMat [lindex $PropertyList [expr { $SearchInList+1 }]]
                            }
                        }
                        set Connectivities [GiD_Mesh get element [lindex [lindex $Entities 0] $j] face [lindex [lindex $Entities 1] $j]]
                        puts $varfile "  $ConditionId  $ElemsMat  $Connectivities"
                    }
                    puts $varfile "End Conditions"
                    puts $varfile ""
                }
                # Update ConditionList
                lappend ConditionList [lindex [lindex $Groups $i] 1] $AuxList
                set AuxList [list]
            }
        }
    }
    # Interface_Face_Load
    set Groups [GiD_Info conditions Interface_Face_Load groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
        # UPwFaceLoadInterfaceCondition2D2N
        set Entities [GiD_EntitiesGroups get [lindex [lindex $Groups $i] 1] elements -element_type linear]
        if {[llength $Entities] > 0} {
            puts $varfile "Begin Conditions UPwFaceLoadInterfaceCondition2D2N"
            for {set j 0} {$j < [llength $Entities]} {incr j} {
                incr ConditionId
                lappend AuxList $ConditionId
                set Connectivities [Line2D2Connectivities [lindex $Entities $j]]
                puts $varfile "  $ConditionId  $InterfaceElemsMat  $Connectivities"
            }
            puts $varfile "End Conditions"
            puts $varfile ""
        }
        # UPwFaceLoadInterfaceCondition3D4N
        set Entities [GiD_EntitiesGroups get [lindex [lindex $Groups $i] 1] elements -element_type triangle]
        if {[llength $Entities] > 0} {
            puts $varfile "Begin Conditions UPwFaceLoadInterfaceCondition3D4N"
            for {set j 0} {$j < [llength $Entities]} {incr j} {
                incr ConditionId
                lappend AuxList $ConditionId
                set Connectivities [TriangleInterface3D4Connectivities [lindex $Entities $j]]
                puts $varfile "  $ConditionId  $InterfaceElemsMat  $Connectivities"
            }
            puts $varfile "End Conditions"
            puts $varfile ""
        }
        set Entities [GiD_EntitiesGroups get [lindex [lindex $Groups $i] 1] elements -element_type quadrilateral]
        if {[llength $Entities] > 0} {
            puts $varfile "Begin Conditions UPwFaceLoadInterfaceCondition3D4N"
            for {set j 0} {$j < [llength $Entities]} {incr j} {
                incr ConditionId
                lappend AuxList $ConditionId
                set Connectivities [QuadrilateralInterface3D4Connectivities [lindex $Entities $j]]
                puts $varfile "  $ConditionId  $InterfaceElemsMat  $Connectivities"
            }
            puts $varfile "End Conditions"
            puts $varfile ""
        }
        # Update ConditionList
        lappend ConditionList [lindex [lindex $Groups $i] 1] $AuxList
        set AuxList [list]
    }
    # Interface_Normal_Fluid_Flux
    set Groups [GiD_Info conditions Interface_Normal_Fluid_Flux groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
        # UPwNormalFluxInterfaceCondition2D2N
        set Entities [GiD_EntitiesGroups get [lindex [lindex $Groups $i] 1] elements -element_type linear]
        if {[llength $Entities] > 0} {
            puts $varfile "Begin Conditions UPwNormalFluxInterfaceCondition2D2N"
            for {set j 0} {$j < [llength $Entities]} {incr j} {
                incr ConditionId
                lappend AuxList $ConditionId
                set Connectivities [Line2D2Connectivities [lindex $Entities $j]]
                puts $varfile "  $ConditionId  $InterfaceElemsMat  $Connectivities"
            }
            puts $varfile "End Conditions"
            puts $varfile ""
        }
        # UPwNormalFluxInterfaceCondition3D4N
        set Entities [GiD_EntitiesGroups get [lindex [lindex $Groups $i] 1] elements -element_type triangle]
        if {[llength $Entities] > 0} {
            puts $varfile "Begin Conditions UPwNormalFluxInterfaceCondition3D4N"
            for {set j 0} {$j < [llength $Entities]} {incr j} {
                incr ConditionId
                lappend AuxList $ConditionId
                set Connectivities [TriangleInterface3D4Connectivities [lindex $Entities $j]]
                puts $varfile "  $ConditionId  $InterfaceElemsMat  $Connectivities"
            }
            puts $varfile "End Conditions"
            puts $varfile ""
        }
        set Entities [GiD_EntitiesGroups get [lindex [lindex $Groups $i] 1] elements -element_type quadrilateral]
        if {[llength $Entities] > 0} {
            puts $varfile "Begin Conditions UPwNormalFluxInterfaceCondition3D4N"
            for {set j 0} {$j < [llength $Entities]} {incr j} {
                incr ConditionId
                lappend AuxList $ConditionId
                set Connectivities [QuadrilateralInterface3D4Connectivities [lindex $Entities $j]]
                puts $varfile "  $ConditionId  $InterfaceElemsMat  $Connectivities"
            }
            puts $varfile "End Conditions"
            puts $varfile ""
        }
        # Update ConditionList
        lappend ConditionList [lindex [lindex $Groups $i] 1] $AuxList
        set AuxList [list]
    }
    puts $varfile ""
    
    ## SubModelParts
    # Body_Part
    set Groups [GiD_Info conditions Body_Part groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
        puts $varfile "Begin SubModelPart [lindex [lindex $Groups $i] 1]"
        # Nodes
        set Entities [GiD_EntitiesGroups get [lindex [lindex $Groups $i] 1] nodes]
        puts $varfile "  Begin SubModelPartNodes"
        for {set j 0} {$j < [llength $Entities]} {incr j} {
            puts $varfile "    [lindex $Entities $j]"
        }
        puts $varfile "  End SubModelPartNodes"
        # Elements
        set Entities [GiD_EntitiesGroups get [lindex [lindex $Groups $i] 1] elements]
        puts $varfile "  Begin SubModelPartElements"
        for {set j 0} {$j < [llength $Entities]} {incr j} {
            puts $varfile "    [lindex $Entities $j]"
        }
        puts $varfile "  End SubModelPartElements"
        # Conditions
        puts $varfile "  Begin SubModelPartConditions"
        puts $varfile "  End SubModelPartConditions"
        puts $varfile "End SubModelPart"
        puts $varfile ""
    }
    # Interface_Part
    set Groups [GiD_Info conditions Interface_Part groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
        puts $varfile "Begin SubModelPart [lindex [lindex $Groups $i] 1]"
        # Nodes
        set Entities [GiD_EntitiesGroups get [lindex [lindex $Groups $i] 1] nodes]
        puts $varfile "  Begin SubModelPartNodes"
        for {set j 0} {$j < [llength $Entities]} {incr j} {
            puts $varfile "    [lindex $Entities $j]"
        }
        puts $varfile "  End SubModelPartNodes"
        # Elements
        set Entities [GiD_EntitiesGroups get [lindex [lindex $Groups $i] 1] elements]
        puts $varfile "  Begin SubModelPartElements"
        for {set j 0} {$j < [llength $Entities]} {incr j} {
            puts $varfile "    [lindex $Entities $j]"
        }
        puts $varfile "  End SubModelPartElements"
        # Conditions
        puts $varfile "  Begin SubModelPartConditions"
        puts $varfile "  End SubModelPartConditions"
        puts $varfile "End SubModelPart"
        puts $varfile ""
    }
    # Solid_Displacement
    set Groups [GiD_Info conditions Solid_Displacement groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
        puts $varfile "Begin SubModelPart [lindex [lindex $Groups $i] 1]"
        # Tables
        set SearchInList [lsearch $TableList [lindex [lindex $Groups $i] 1]*]
        set AuxList [lindex $TableList [expr { $SearchInList+1 }]]
        puts $varfile "  Begin SubModelPartTables"
        for {set j 0} {$j < [llength $AuxList]} {incr j} {
            if {[lindex $AuxList $j] > 0} {
                puts $varfile "    [lindex $AuxList $j]"
            }
        }
        puts $varfile "  End SubModelPartTables"
        # Nodes
        set Entities [GiD_EntitiesGroups get [lindex [lindex $Groups $i] 1] nodes]
        puts $varfile "  Begin SubModelPartNodes"
        for {set j 0} {$j < [llength $Entities]} {incr j} {
            puts $varfile "    [lindex $Entities $j]"
        }
        puts $varfile "  End SubModelPartNodes"
        # Elements
        puts $varfile "  Begin SubModelPartElements"
        puts $varfile "  End SubModelPartElements"
        # Conditions
        puts $varfile "  Begin SubModelPartConditions"
        puts $varfile "  End SubModelPartConditions"
        puts $varfile "End SubModelPart"
        puts $varfile ""
    }
    # Fluid_Pressure
    set Groups [GiD_Info conditions Fluid_Pressure groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
        puts $varfile "Begin SubModelPart [lindex [lindex $Groups $i] 1]"
        # Tables
        set SearchInList [lsearch $TableList [lindex [lindex $Groups $i] 1]*]
        set AuxList [lindex $TableList [expr { $SearchInList+1 }]]
        puts $varfile "  Begin SubModelPartTables"
        for {set j 0} {$j < [llength $AuxList]} {incr j} {
            if {[lindex $AuxList $j] > 0} {
                puts $varfile "    [lindex $AuxList $j]"
            }
        }
        puts $varfile "  End SubModelPartTables"
        # Nodes
        set Entities [GiD_EntitiesGroups get [lindex [lindex $Groups $i] 1] nodes]
        puts $varfile "  Begin SubModelPartNodes"
        for {set j 0} {$j < [llength $Entities]} {incr j} {
            puts $varfile "    [lindex $Entities $j]"
        }
        puts $varfile "  End SubModelPartNodes"
        # Elements
        puts $varfile "  Begin SubModelPartElements"
        puts $varfile "  End SubModelPartElements"
        # Conditions
        puts $varfile "  Begin SubModelPartConditions"
        puts $varfile "  End SubModelPartConditions"
        puts $varfile "End SubModelPart"
        puts $varfile ""
    }
    # Force
    set Groups [GiD_Info conditions Force groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
        puts $varfile "Begin SubModelPart [lindex [lindex $Groups $i] 1]"
        # Tables
        set SearchInList [lsearch $TableList [lindex [lindex $Groups $i] 1]*]
        set AuxList [lindex $TableList [expr { $SearchInList+1 }]]
        puts $varfile "  Begin SubModelPartTables"
        for {set j 0} {$j < [llength $AuxList]} {incr j} {
            if {[lindex $AuxList $j] > 0} {
                puts $varfile "    [lindex $AuxList $j]"
            }
        }
        puts $varfile "  End SubModelPartTables"
        # Nodes
        set Entities [GiD_EntitiesGroups get [lindex [lindex $Groups $i] 1] nodes]
        puts $varfile "  Begin SubModelPartNodes"
        for {set j 0} {$j < [llength $Entities]} {incr j} {
            puts $varfile "    [lindex $Entities $j]"
        }
        puts $varfile "  End SubModelPartNodes"
        # Elements
        puts $varfile "  Begin SubModelPartElements"
        puts $varfile "  End SubModelPartElements"
        #Conditions
        set SearchInList [lsearch $ConditionList [lindex [lindex $Groups $i] 1]*]
        set AuxList [lindex $ConditionList [expr { $SearchInList+1 }]]
        puts $varfile "  Begin SubModelPartConditions"
        for {set j 0} {$j < [llength $AuxList]} {incr j} {
            puts $varfile "    [lindex $AuxList $j]"
        }
        puts $varfile "  End SubModelPartConditions"
        puts $varfile "End SubModelPart"
        puts $varfile ""
    }
    # Face_Load
    set Groups [GiD_Info conditions Face_Load groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
        puts $varfile "Begin SubModelPart [lindex [lindex $Groups $i] 1]"
        # Tables
        set SearchInList [lsearch $TableList [lindex [lindex $Groups $i] 1]*]
        set AuxList [lindex $TableList [expr { $SearchInList+1 }]]
        puts $varfile "  Begin SubModelPartTables"
        for {set j 0} {$j < [llength $AuxList]} {incr j} {
            if {[lindex $AuxList $j] > 0} {
                puts $varfile "    [lindex $AuxList $j]"
            }
        }
        puts $varfile "  End SubModelPartTables"
        # Nodes
        set Entities [GiD_EntitiesGroups get [lindex [lindex $Groups $i] 1] nodes]
        puts $varfile "  Begin SubModelPartNodes"
        for {set j 0} {$j < [llength $Entities]} {incr j} {
            puts $varfile "    [lindex $Entities $j]"
        }
        puts $varfile "  End SubModelPartNodes"
        # Elements
        puts $varfile "  Begin SubModelPartElements"
        puts $varfile "  End SubModelPartElements"
        #Conditions
        set SearchInList [lsearch $ConditionList [lindex [lindex $Groups $i] 1]*]
        set AuxList [lindex $ConditionList [expr { $SearchInList+1 }]]
        puts $varfile "  Begin SubModelPartConditions"
        for {set j 0} {$j < [llength $AuxList]} {incr j} {
            puts $varfile "    [lindex $AuxList $j]"
        }
        puts $varfile "  End SubModelPartConditions"
        puts $varfile "End SubModelPart"
        puts $varfile ""
    }
    # Normal_Load
    set Groups [GiD_Info conditions Normal_Load groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
        puts $varfile "Begin SubModelPart [lindex [lindex $Groups $i] 1]"
        # Tables
        set SearchInList [lsearch $TableList [lindex [lindex $Groups $i] 1]*]
        set AuxList [lindex $TableList [expr { $SearchInList+1 }]]
        puts $varfile "  Begin SubModelPartTables"
        for {set j 0} {$j < [llength $AuxList]} {incr j} {
            if {[lindex $AuxList $j] > 0} {
                puts $varfile "    [lindex $AuxList $j]"
            }
        }
        puts $varfile "  End SubModelPartTables"
        # Nodes
        set Entities [GiD_EntitiesGroups get [lindex [lindex $Groups $i] 1] nodes]
        puts $varfile "  Begin SubModelPartNodes"
        for {set j 0} {$j < [llength $Entities]} {incr j} {
            puts $varfile "    [lindex $Entities $j]"
        }
        puts $varfile "  End SubModelPartNodes"
        # Elements
        puts $varfile "  Begin SubModelPartElements"
        puts $varfile "  End SubModelPartElements"
        #Conditions
        set SearchInList [lsearch $ConditionList [lindex [lindex $Groups $i] 1]*]
        set AuxList [lindex $ConditionList [expr { $SearchInList+1 }]]
        puts $varfile "  Begin SubModelPartConditions"
        for {set j 0} {$j < [llength $AuxList]} {incr j} {
            puts $varfile "    [lindex $AuxList $j]"
        }
        puts $varfile "  End SubModelPartConditions"
        puts $varfile "End SubModelPart"
        puts $varfile ""
    }
    # Normal_Fluid_Flux
    set Groups [GiD_Info conditions Normal_Fluid_Flux groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
        puts $varfile "Begin SubModelPart [lindex [lindex $Groups $i] 1]"
        # Tables
        set SearchInList [lsearch $TableList [lindex [lindex $Groups $i] 1]*]
        set AuxList [lindex $TableList [expr { $SearchInList+1 }]]
        puts $varfile "  Begin SubModelPartTables"
        for {set j 0} {$j < [llength $AuxList]} {incr j} {
            if {[lindex $AuxList $j] > 0} {
                puts $varfile "    [lindex $AuxList $j]"
            }
        }
        puts $varfile "  End SubModelPartTables"
        # Nodes
        set Entities [GiD_EntitiesGroups get [lindex [lindex $Groups $i] 1] nodes]
        puts $varfile "  Begin SubModelPartNodes"
        for {set j 0} {$j < [llength $Entities]} {incr j} {
            puts $varfile "    [lindex $Entities $j]"
        }
        puts $varfile "  End SubModelPartNodes"
        # Elements
        puts $varfile "  Begin SubModelPartElements"
        puts $varfile "  End SubModelPartElements"
        #Conditions
        set SearchInList [lsearch $ConditionList [lindex [lindex $Groups $i] 1]*]
        set AuxList [lindex $ConditionList [expr { $SearchInList+1 }]]
        puts $varfile "  Begin SubModelPartConditions"
        for {set j 0} {$j < [llength $AuxList]} {incr j} {
            puts $varfile "    [lindex $AuxList $j]"
        }
        puts $varfile "  End SubModelPartConditions"
        puts $varfile "End SubModelPart"
        puts $varfile ""
    }
    # Interface_Face_Load
    set Groups [GiD_Info conditions Interface_Face_Load groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
        puts $varfile "Begin SubModelPart [lindex [lindex $Groups $i] 1]"
        # Tables
        set SearchInList [lsearch $TableList [lindex [lindex $Groups $i] 1]*]
        set AuxList [lindex $TableList [expr { $SearchInList+1 }]]
        puts $varfile "  Begin SubModelPartTables"
        for {set j 0} {$j < [llength $AuxList]} {incr j} {
            if {[lindex $AuxList $j] > 0} {
                puts $varfile "    [lindex $AuxList $j]"
            }
        }
        puts $varfile "  End SubModelPartTables"
        # Nodes
        set Entities [GiD_EntitiesGroups get [lindex [lindex $Groups $i] 1] nodes]
        puts $varfile "  Begin SubModelPartNodes"
        for {set j 0} {$j < [llength $Entities]} {incr j} {
            puts $varfile "    [lindex $Entities $j]"
        }
        puts $varfile "  End SubModelPartNodes"
        # Elements
        puts $varfile "  Begin SubModelPartElements"
        puts $varfile "  End SubModelPartElements"
        #Conditions
        set SearchInList [lsearch $ConditionList [lindex [lindex $Groups $i] 1]*]
        set AuxList [lindex $ConditionList [expr { $SearchInList+1 }]]
        puts $varfile "  Begin SubModelPartConditions"
        for {set j 0} {$j < [llength $AuxList]} {incr j} {
            puts $varfile "    [lindex $AuxList $j]"
        }
        puts $varfile "  End SubModelPartConditions"
        puts $varfile "End SubModelPart"
        puts $varfile ""
    }
    # Interface_Normal_Fluid_Flux
    set Groups [GiD_Info conditions Interface_Normal_Fluid_Flux groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
        puts $varfile "Begin SubModelPart [lindex [lindex $Groups $i] 1]"
        # Tables
        set SearchInList [lsearch $TableList [lindex [lindex $Groups $i] 1]*]
        set AuxList [lindex $TableList [expr { $SearchInList+1 }]]
        puts $varfile "  Begin SubModelPartTables"
        for {set j 0} {$j < [llength $AuxList]} {incr j} {
            if {[lindex $AuxList $j] > 0} {
                puts $varfile "    [lindex $AuxList $j]"
            }
        }
        puts $varfile "  End SubModelPartTables"
        # Nodes
        set Entities [GiD_EntitiesGroups get [lindex [lindex $Groups $i] 1] nodes]
        puts $varfile "  Begin SubModelPartNodes"
        for {set j 0} {$j < [llength $Entities]} {incr j} {
            puts $varfile "    [lindex $Entities $j]"
        }
        puts $varfile "  End SubModelPartNodes"
        # Elements
        puts $varfile "  Begin SubModelPartElements"
        puts $varfile "  End SubModelPartElements"
        #Conditions
        set SearchInList [lsearch $ConditionList [lindex [lindex $Groups $i] 1]*]
        set AuxList [lindex $ConditionList [expr { $SearchInList+1 }]]
        puts $varfile "  Begin SubModelPartConditions"
        for {set j 0} {$j < [llength $AuxList]} {incr j} {
            puts $varfile "    [lindex $AuxList $j]"
        }
        puts $varfile "  End SubModelPartConditions"
        puts $varfile "End SubModelPart"
        puts $varfile ""
    }
    # Body_Acceleration
    set Groups [GiD_Info conditions Body_Acceleration groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
        puts $varfile "Begin SubModelPart [lindex [lindex $Groups $i] 1]"
        # Tables
        set SearchInList [lsearch $TableList [lindex [lindex $Groups $i] 1]*]
        set AuxList [lindex $TableList [expr { $SearchInList+1 }]]
        puts $varfile "  Begin SubModelPartTables"
        for {set j 0} {$j < [llength $AuxList]} {incr j} {
            if {[lindex $AuxList $j] > 0} {
                puts $varfile "    [lindex $AuxList $j]"
            }
        }
        puts $varfile "  End SubModelPartTables"
        # Nodes
        set Entities [GiD_EntitiesGroups get [lindex [lindex $Groups $i] 1] nodes]
        puts $varfile "  Begin SubModelPartNodes"
        for {set j 0} {$j < [llength $Entities]} {incr j} {
            puts $varfile "    [lindex $Entities $j]"
        }
        puts $varfile "  End SubModelPartNodes"
        # Elements
        puts $varfile "  Begin SubModelPartElements"
        puts $varfile "  End SubModelPartElements"
        # Conditions
        puts $varfile "  Begin SubModelPartConditions"
        puts $varfile "  End SubModelPartConditions"
        puts $varfile "End SubModelPart"
        puts $varfile ""
    }

    close $varfile
    
    return $TableList
}

#-------------------------------------------------------------------------------

proc Triangle2D3Connectivities { ElemId } {
    
    set ElementInfo [GiD_Mesh get element $ElemId]
    #ElementInfo: <layer> <elemtype> <NumNodes> <N1> <N2> ...
    for {set i 0} {$i < 3} {incr i} {
        set N($i) [lindex $ElementInfo [expr { 3+$i }]]
    }
    set Connectivities "$N(0) $N(1) $N(2)"
    
    return $Connectivities
}

#-------------------------------------------------------------------------------

proc Quadrilateral2D4Connectivities { ElemId } {
    
    #It is the same for the Tethrahedron3D4
    
    set ElementInfo [GiD_Mesh get element $ElemId]
    #ElementInfo: <layer> <elemtype> <NumNodes> <N1> <N2> ...
    for {set i 0} {$i < 4} {incr i} {
        set N($i) [lindex $ElementInfo [expr { 3+$i }]]
    }
    
    set Connectivities "$N(0) $N(1) $N(2) $N(3)"
    
    return $Connectivities
}

#-------------------------------------------------------------------------------

proc Triangle2D6Connectivities { ElemId } {
    
    #It is the same for the Prism3D6
    
    set ElementInfo [GiD_Mesh get element $ElemId]
    #ElementInfo: <layer> <elemtype> <NumNodes> <N1> <N2> ...
    for {set i 0} {$i < 6} {incr i} {
        set N($i) [lindex $ElementInfo [expr { 3+$i }]]
    }
    set Connectivities "$N(0) $N(1) $N(2) $N(3) $N(4) $N(5)"
    
    return $Connectivities
}

#-------------------------------------------------------------------------------

proc Hexahedron3D8Connectivities { ElemId } {
    
    #It is the same for Quadrilateral2D8
    
    set ElementInfo [GiD_Mesh get element $ElemId]
    #ElementInfo: <layer> <elemtype> <NumNodes> <N1> <N2> ...
    for {set i 0} {$i < 8} {incr i} {
        set N($i) [lindex $ElementInfo [expr { 3+$i }]]
    }

    set Connectivities "$N(0) $N(1) $N(2) $N(3) $N(4) $N(5) $N(6) $N(7)"
    
    return $Connectivities
}

#-------------------------------------------------------------------------------

proc Quadrilateral2D9Connectivities { ElemId } {
    
    set ElementInfo [GiD_Mesh get element $ElemId]
    #ElementInfo: <layer> <elemtype> <NumNodes> <N1> <N2> ...
    for {set i 0} {$i < 9} {incr i} {
        set N($i) [lindex $ElementInfo [expr { 3+$i }]]
    }

    set Connectivities "$N(0) $N(1) $N(2) $N(3) $N(4) $N(5) $N(6) $N(7) $N(8)"
    
    return $Connectivities
}

#-------------------------------------------------------------------------------

proc Tetrahedron3D10Connectivities { ElemId } {
        
    set ElementInfo [GiD_Mesh get element $ElemId]
    #ElementInfo: <layer> <elemtype> <NumNodes> <N1> <N2> ...
    for {set i 0} {$i < 10} {incr i} {
        set N($i) [lindex $ElementInfo [expr { 3+$i }]]
    }

    set Connectivities "$N(0) $N(1) $N(2) $N(3) $N(4) $N(5) $N(6) $N(7) $N(8) $N(9)"
    
    return $Connectivities
}

#-------------------------------------------------------------------------------

proc Hexahedron3D20Connectivities { ElemId } {
    
    set ElementInfo [GiD_Mesh get element $ElemId]
    #ElementInfo: <layer> <elemtype> <NumNodes> <N1> <N2> ...
    set Connectivities [list]
    for {set i 0} {$i < 20} {incr i} {
        set N($i) [lindex $ElementInfo [expr { 3+$i }]]
        lappend Connectivities $N($i)
    }

    return $Connectivities
}

#-------------------------------------------------------------------------------

proc Hexahedron3D27Connectivities { ElemId } {
        
    set ElementInfo [GiD_Mesh get element $ElemId]
    #ElementInfo: <layer> <elemtype> <NumNodes> <N1> <N2> ...
    set Connectivities [list]
    for {set i 0} {$i < 27} {incr i} {
        set N($i) [lindex $ElementInfo [expr { 3+$i }]]
        lappend Connectivities $N($i)
    }

    return $Connectivities
}

#-------------------------------------------------------------------------------

proc TriangleInterface2D4Connectivities { ElemId } {
    
    # Obtaining element nodes
    set ElementInfo [GiD_Mesh get element $ElemId]
    #ElementInfo: <layer> <elemtype> <NumNodes> <N1> <N2> ...
    set N1(Id) [lindex $ElementInfo 3]
    set N2(Id) [lindex $ElementInfo 4]
    set N3(Id) $N2(Id)
    set N4(Id) [lindex $ElementInfo 5]
    
    # Obtaining nodes coordinates
    set NCoord [lindex [GiD_Info Coordinates $N1(Id)] 0]
    set N1(x) [lindex $NCoord 0]
    set N1(y) [lindex $NCoord 1]
    set NCoord [lindex [GiD_Info Coordinates $N2(Id)] 0]
    set N2(x) [lindex $NCoord 0]
    set N2(y) [lindex $NCoord 1]
    #set NCoord [lindex [GiD_Info Coordinates $N3(Id)] 0]
    set N3(x) $N2(x)
    set N3(y) $N2(y)
    set NCoord [lindex [GiD_Info Coordinates $N4(Id)] 0]
    set N4(x) [lindex $NCoord 0]
    set N4(y) [lindex $NCoord 1]
    
    # Computing element lengths
    set lx [expr { sqrt( (0.5*($N2(x)+$N3(x)-$N1(x)-$N4(x)))**2 + (0.5*($N2(y)+$N3(y)-$N1(y)-$N4(y)))**2 ) }]
    set ly [expr { sqrt( (0.5*($N3(x)+$N4(x)-$N1(x)-$N2(x)))**2 + (0.5*($N3(y)+$N4(y)-$N1(y)-$N2(y)))**2 ) }]
    
    if {$ly <= $lx} {
        set Connectivities "$N1(Id) $N2(Id) $N3(Id) $N4(Id)"
    } else {
        set Connectivities "$N4(Id) $N1(Id) $N2(Id) $N3(Id)"
    }
    
    return $Connectivities
}

#-------------------------------------------------------------------------------

proc QuadrilateralInterface2D4Connectivities { ElemId } {
        
    # Obtaining element nodes
    set ElementInfo [GiD_Mesh get element $ElemId]
    #ElementInfo: <layer> <elemtype> <NumNodes> <N1> <N2> ...
    set N1(Id) [lindex $ElementInfo 3]
    set N2(Id) [lindex $ElementInfo 4]
    set N3(Id) [lindex $ElementInfo 5]
    set N4(Id) [lindex $ElementInfo 6]
    
    # Obtaining nodes coordinates
    set NCoord [lindex [GiD_Info Coordinates $N1(Id)] 0]
    set N1(x) [lindex $NCoord 0]
    set N1(y) [lindex $NCoord 1]
    set NCoord [lindex [GiD_Info Coordinates $N2(Id)] 0]
    set N2(x) [lindex $NCoord 0]
    set N2(y) [lindex $NCoord 1]
    set NCoord [lindex [GiD_Info Coordinates $N3(Id)] 0]
    set N3(x) [lindex $NCoord 0]
    set N3(y) [lindex $NCoord 1]
    set NCoord [lindex [GiD_Info Coordinates $N4(Id)] 0]
    set N4(x) [lindex $NCoord 0]
    set N4(y) [lindex $NCoord 1]
    
    # Computing element lengths
    set lx [expr { sqrt( (0.5*($N2(x)+$N3(x)-$N1(x)-$N4(x)))**2 + (0.5*($N2(y)+$N3(y)-$N1(y)-$N4(y)))**2 ) }]
    set ly [expr { sqrt( (0.5*($N3(x)+$N4(x)-$N1(x)-$N2(x)))**2 + (0.5*($N3(y)+$N4(y)-$N1(y)-$N2(y)))**2 ) }]
    
    if {$ly <= $lx} {
        set Connectivities "$N1(Id) $N2(Id) $N3(Id) $N4(Id)"
    } else {
        set Connectivities "$N4(Id) $N1(Id) $N2(Id) $N3(Id)"
    }
    
    return $Connectivities
}

#-------------------------------------------------------------------------------

proc TetrahedronInterface3D6Connectivities { ElemId } {
        
    set ElementInfo [GiD_Mesh get element $ElemId]
    #ElementInfo: <layer> <elemtype> <NumNodes> <N1> <N2> ...
    set N1(Id) [lindex $ElementInfo 3]
    set N2(Id) [lindex $ElementInfo 4]
    set N3(Id) [lindex $ElementInfo 5]
    set N4(Id) [lindex $ElementInfo 6]
    
    set Connectivities "$N1(Id) $N2(Id) $N3(Id) $N4(Id) $N2(Id) $N3(Id)"
    
    return $Connectivities
}

#-------------------------------------------------------------------------------

#proc TetrahedronInterface3D6Connectivities { ElemId } {
        
    #set ElementInfo [GiD_Mesh get element $ElemId]
    ##ElementInfo: <layer> <elemtype> <NumNodes> <N1> <N2> ...
    #set N1(Id) [lindex $ElementInfo 3]
    #set N2(Id) [lindex $ElementInfo 4]
    #set N3(Id) [lindex $ElementInfo 5]
    #set N4(Id) [lindex $ElementInfo 6]
    
    ## Obtaining nodes coordinates
    #set NCoord [lindex [GiD_Info Coordinates $N1(Id)] 0]
    #set N1(x) [lindex $NCoord 0]
    #set N1(y) [lindex $NCoord 1]
    #set N1(z) [lindex $NCoord 2]
    #set NCoord [lindex [GiD_Info Coordinates $N2(Id)] 0]
    #set N2(x) [lindex $NCoord 0]
    #set N2(y) [lindex $NCoord 1]
    #set N2(z) [lindex $NCoord 2]
    #set NCoord [lindex [GiD_Info Coordinates $N3(Id)] 0]
    #set N3(x) [lindex $NCoord 0]
    #set N3(y) [lindex $NCoord 1]
    #set N3(z) [lindex $NCoord 2]
    #set NCoord [lindex [GiD_Info Coordinates $N4(Id)] 0]
    #set N4(x) [lindex $NCoord 0]
    #set N4(y) [lindex $NCoord 1]
    #set N4(z) [lindex $NCoord 2]
    
    ## Compute vector from node 1 to node 4
    #set Vup(x) [expr { $N4(x)-$N1(x) }]
    #set Vup(y) [expr { $N4(y)-$N1(y) }]
    #set Vup(z) [expr { $N4(z)-$N1(z) }]
    
    ## Compute vector from node 1 to node 2
    #set V12(x) [expr { $N2(x)-$N1(x) }]
    #set V12(y) [expr { $N2(y)-$N1(y) }]
    #set V12(z) [expr { $N2(z)-$N1(z) }]
    ## Compute vector from node 1 to node 3
    #set V13(x) [expr { $N3(x)-$N1(x) }]
    #set V13(y) [expr { $N3(y)-$N1(y) }]
    #set V13(z) [expr { $N3(z)-$N1(z) }]
    
    ## Compute cross product between V12 and V13
    #set Vn(x) [expr { $V12(y)*$V13(z)-$V12(z)*$V13(y) }]
    #set Vn(y) [expr { $V12(z)*$V13(x)-$V12(x)*$V13(z) }]
    #set Vn(z) [expr { $V12(x)*$V13(y)-$V12(y)*$V13(x) }]
    
    ## Compute cos(angle) between Vup and Vn
    #set cosangle [expr { $Vup(x)*$Vn(x)+$Vup(y)*$Vn(y)+$Vup(z)*$Vn(z)}]
    
    ## Check Connectivities
    #if {$cosangle > -1.0e-20} {
        #set Connectivities "$N1(Id) $N2(Id) $N3(Id) $N4(Id) $N2(Id) $N3(Id)"
    #} else {
        ##W "Reordering nodes of interface element"
        #set Connectivities "$N1(Id) $N3(Id) $N2(Id) $N4(Id) $N3(Id) $N2(Id)"
    #}
    
    #return $Connectivities
#}

#-------------------------------------------------------------------------------

proc PrismInterface3D6Connectivities { ElemId } {

    # Obtaining element nodes
    set ElementInfo [GiD_Mesh get element $ElemId]
    #ElementInfo: <layer> <elemtype> <NumNodes> <N1> <N2> ...
    set N1(Id) [lindex $ElementInfo 3]
    set N2(Id) [lindex $ElementInfo 4]
    set N3(Id) [lindex $ElementInfo 5]
    set N4(Id) [lindex $ElementInfo 6]
    set N5(Id) [lindex $ElementInfo 7]
    set N6(Id) [lindex $ElementInfo 8]
    
    # Obtaining element volume
    set Volume [lindex [GiD_Info list_entities -more Elements $ElemId] 20]
    set Volume [string trimleft $Volume "Volume="]
    
    # Check Connectivities
    if {$Volume > -1.0e-10} {
        set Connectivities "$N1(Id) $N2(Id) $N3(Id) $N4(Id) $N5(Id) $N6(Id)"
    } else {
        #W "Reordering nodes of interface element"
        set Connectivities "$N1(Id) $N3(Id) $N2(Id) $N4(Id) $N6(Id) $N5(Id)"
    }
    
    return $Connectivities
}

#-------------------------------------------------------------------------------

proc HexahedronInterface3D8Connectivities { ElemId } {

    # Obtaining element nodes
    set ElementInfo [GiD_Mesh get element $ElemId]
    #ElementInfo: <layer> <elemtype> <NumNodes> <N1> <N2> ...
    set N1(Id) [lindex $ElementInfo 3]
    set N2(Id) [lindex $ElementInfo 4]
    set N3(Id) [lindex $ElementInfo 5]
    set N4(Id) [lindex $ElementInfo 6]
    set N5(Id) [lindex $ElementInfo 7]
    set N6(Id) [lindex $ElementInfo 8]
    set N7(Id) [lindex $ElementInfo 9]
    set N8(Id) [lindex $ElementInfo 10]
    
    # Obtaining element volume
    set Volume [lindex [GiD_Info list_entities -more Elements $ElemId] 22]
    set Volume [string trimleft $Volume "Volume="]
    
    # Check Connectivities
    if {$Volume > -1.0e-10} {
        set Connectivities "$N1(Id) $N2(Id) $N3(Id) $N4(Id) $N5(Id) $N6(Id) $N7(Id) $N8(Id)"
    } else {
        #W "Reordering nodes of interface element"
        set Connectivities "$N1(Id) $N4(Id) $N3(Id) $N2(Id) $N5(Id) $N8(Id) $N7(Id) $N6(Id)"
    }
    
    return $Connectivities
}

#-------------------------------------------------------------------------------

#proc HexahedronInterface3D8Connectivities { ElemId } {
    
    ## Obtaining element nodes
    #set ElementInfo [GiD_Mesh get element $ElemId]
    ##ElementInfo: <layer> <elemtype> <NumNodes> <N1> <N2> ...
    #set N1(Id) [lindex $ElementInfo 3]
    #set N2(Id) [lindex $ElementInfo 4]
    #set N3(Id) [lindex $ElementInfo 5]
    #set N4(Id) [lindex $ElementInfo 6]
    #set N5(Id) [lindex $ElementInfo 7]
    #set N6(Id) [lindex $ElementInfo 8]
    #set N7(Id) [lindex $ElementInfo 9]
    #set N8(Id) [lindex $ElementInfo 10]
    
    ## Obtaining nodes coordinates
    #set NCoord [lindex [GiD_Info Coordinates $N1(Id)] 0]
    #set N1(x) [lindex $NCoord 0]
    #set N1(y) [lindex $NCoord 1]
    #set N1(z) [lindex $NCoord 2]
    #set NCoord [lindex [GiD_Info Coordinates $N2(Id)] 0]
    #set N2(x) [lindex $NCoord 0]
    #set N2(y) [lindex $NCoord 1]
    #set N2(z) [lindex $NCoord 2]
    #set NCoord [lindex [GiD_Info Coordinates $N3(Id)] 0]
    #set N3(x) [lindex $NCoord 0]
    #set N3(y) [lindex $NCoord 1]
    #set N3(z) [lindex $NCoord 2]
    #set NCoord [lindex [GiD_Info Coordinates $N4(Id)] 0]
    #set N4(x) [lindex $NCoord 0]
    #set N4(y) [lindex $NCoord 1]
    #set N4(z) [lindex $NCoord 2]
    #set NCoord [lindex [GiD_Info Coordinates $N5(Id)] 0]
    #set N5(x) [lindex $NCoord 0]
    #set N5(y) [lindex $NCoord 1]
    #set N5(z) [lindex $NCoord 2]
    #set NCoord [lindex [GiD_Info Coordinates $N6(Id)] 0]
    #set N6(x) [lindex $NCoord 0]
    #set N6(y) [lindex $NCoord 1]
    #set N6(z) [lindex $NCoord 2]
    #set NCoord [lindex [GiD_Info Coordinates $N7(Id)] 0]
    #set N7(x) [lindex $NCoord 0]
    #set N7(y) [lindex $NCoord 1]
    #set N7(z) [lindex $NCoord 2]
    #set NCoord [lindex [GiD_Info Coordinates $N8(Id)] 0]
    #set N8(x) [lindex $NCoord 0]
    #set N8(y) [lindex $NCoord 1]
    #set N8(z) [lindex $NCoord 2]
    
    ## Computing element lengths
    #set lx [expr { sqrt( (0.25*($N2(x)+$N6(x)+$N3(x)+$N7(x)-$N1(x)-$N5(x)-$N4(x)-$N8(x)))**2 + (0.25*($N2(y)+$N6(y)+$N3(y)+$N7(y)-$N1(y)-$N5(y)-$N4(y)-$N8(y)))**2 + (0.25*($N2(z)+$N6(z)+$N3(z)+$N7(z)-$N1(z)-$N5(z)-$N4(z)-$N8(z)))**2 ) }]
    #set ly [expr { sqrt( (0.25*($N3(x)+$N4(x)+$N7(x)+$N8(x)-$N1(x)-$N2(x)-$N5(x)-$N6(x)))**2 + (0.25*($N3(y)+$N4(y)+$N7(y)+$N8(y)-$N1(y)-$N2(y)-$N5(y)-$N6(y)))**2 + (0.25*($N3(z)+$N4(z)+$N7(z)+$N8(z)-$N1(z)-$N2(z)-$N5(z)-$N6(z)))**2 ) }]
    #set lz [expr { sqrt( (0.25*($N5(x)+$N6(x)+$N7(x)+$N8(x)-$N1(x)-$N2(x)-$N3(x)-$N4(x)))**2 + (0.25*($N5(y)+$N6(y)+$N7(y)+$N8(y)-$N1(y)-$N2(y)-$N3(y)-$N4(y)))**2 + (0.25*($N5(z)+$N6(z)+$N7(z)+$N8(z)-$N1(z)-$N2(z)-$N3(z)-$N4(z)))**2 ) }]
    
    #if {$lz <= $lx} {
        #if {$lz <= $ly} {
            ## lz <= lx && lz <= ly
            #set Connectivities "$N1(Id) $N2(Id) $N3(Id) $N4(Id) $N5(Id) $N6(Id) $N7(Id) $N8(Id)"
        #} else {
            ## ly < lz <= lx
            #set Connectivities "$N1(Id) $N4(Id) $N8(Id) $N5(Id) $N2(Id) $N3(Id) $N7(Id) $N6(Id)"
        #}
    #} elseif {$ly <= $lx} {
        ## ly <= lx < lz
        #set Connectivities "$N1(Id) $N4(Id) $N8(Id) $N5(Id) $N2(Id) $N3(Id) $N7(Id) $N6(Id)"
    #} else {
        ## lx < lz && lx < ly
        #set Connectivities "$N1(Id) $N5(Id) $N6(Id) $N2(Id) $N4(Id) $N8(Id) $N7(Id) $N3(Id)"
    #}
    
    #return $Connectivities
#}

#-------------------------------------------------------------------------------

proc Line2D2Connectivities { ElemId } {
    
    set ElementInfo [GiD_Mesh get element $ElemId]
    #ElementInfo: <layer> <elemtype> <NumNodes> <N1> <N2> ...
    for {set i 0} {$i < 2} {incr i} {
        set N($i) [lindex $ElementInfo [expr { 3+$i }]]
    }
    set Connectivities "$N(0) $N(1)"
    
    return $Connectivities
}

#-------------------------------------------------------------------------------

proc TriangleInterface3D4Connectivities { ElemId } {
    
    # Obtaining element nodes
    set ElementInfo [GiD_Mesh get element $ElemId]
    #ElementInfo: <layer> <elemtype> <NumNodes> <N1> <N2> ...
    set N1(Id) [lindex $ElementInfo 3]
    set N2(Id) [lindex $ElementInfo 4]
    set N3(Id) $N2(Id)
    set N4(Id) [lindex $ElementInfo 5]
    
    # Obtaining nodes coordinates
    set NCoord [lindex [GiD_Info Coordinates $N1(Id)] 0]
    set N1(x) [lindex $NCoord 0]
    set N1(y) [lindex $NCoord 1]
    set N1(z) [lindex $NCoord 2]
    set NCoord [lindex [GiD_Info Coordinates $N2(Id)] 0]
    set N2(x) [lindex $NCoord 0]
    set N2(y) [lindex $NCoord 1]
    set N2(z) [lindex $NCoord 2]
    #set NCoord [lindex [GiD_Info Coordinates $N3(Id)] 0]
    set N3(x) $N2(x)
    set N3(y) $N2(y)
    set N3(z) $N2(z)
    set NCoord [lindex [GiD_Info Coordinates $N4(Id)] 0]
    set N4(x) [lindex $NCoord 0]
    set N4(y) [lindex $NCoord 1]
    set N4(z) [lindex $NCoord 2]
    
    # Computing element lengths
    set lx [expr { sqrt( (0.5*($N2(x)+$N3(x)-$N1(x)-$N4(x)))**2 + (0.5*($N2(y)+$N3(y)-$N1(y)-$N4(y)))**2 + (0.5*($N2(z)+$N3(z)-$N1(z)-$N4(z)))**2 ) }]
    set ly [expr { sqrt( (0.5*($N3(x)+$N4(x)-$N1(x)-$N2(x)))**2 + (0.5*($N3(y)+$N4(y)-$N1(y)-$N2(y)))**2 + (0.5*($N3(z)+$N4(z)-$N1(z)-$N2(z)))**2 ) }]
    
    if {$ly <= $lx} {
        set Connectivities "$N1(Id) $N2(Id) $N3(Id) $N4(Id)"
    } else {
        set Connectivities "$N4(Id) $N1(Id) $N2(Id) $N3(Id)"
    }
    
    return $Connectivities
}

#-------------------------------------------------------------------------------

proc QuadrilateralInterface3D4Connectivities { ElemId } {
    
    # Obtaining element nodes
    set ElementInfo [GiD_Mesh get element $ElemId]
    #ElementInfo: <layer> <elemtype> <NumNodes> <N1> <N2> ...
    set N1(Id) [lindex $ElementInfo 3]
    set N2(Id) [lindex $ElementInfo 4]
    set N3(Id) [lindex $ElementInfo 5]
    set N4(Id) [lindex $ElementInfo 6]
    
    # Obtaining nodes coordinates
    set NCoord [lindex [GiD_Info Coordinates $N1(Id)] 0]
    set N1(x) [lindex $NCoord 0]
    set N1(y) [lindex $NCoord 1]
    set N1(z) [lindex $NCoord 2]
    set NCoord [lindex [GiD_Info Coordinates $N2(Id)] 0]
    set N2(x) [lindex $NCoord 0]
    set N2(y) [lindex $NCoord 1]
    set N2(z) [lindex $NCoord 2]
    set NCoord [lindex [GiD_Info Coordinates $N3(Id)] 0]
    set N3(x) [lindex $NCoord 0]
    set N3(y) [lindex $NCoord 1]
    set N3(z) [lindex $NCoord 2]
    set NCoord [lindex [GiD_Info Coordinates $N4(Id)] 0]
    set N4(x) [lindex $NCoord 0]
    set N4(y) [lindex $NCoord 1]
    set N4(z) [lindex $NCoord 2]
    
    # Computing element lengths
    set lx [expr { sqrt( (0.5*($N2(x)+$N3(x)-$N1(x)-$N4(x)))**2 + (0.5*($N2(y)+$N3(y)-$N1(y)-$N4(y)))**2 + (0.5*($N2(z)+$N3(z)-$N1(z)-$N4(z)))**2 ) }]
    set ly [expr { sqrt( (0.5*($N3(x)+$N4(x)-$N1(x)-$N2(x)))**2 + (0.5*($N3(y)+$N4(y)-$N1(y)-$N2(y)))**2 + (0.5*($N3(z)+$N4(z)-$N1(z)-$N2(z)))**2 ) }]
    
    if {$ly <= $lx} {
        set Connectivities "$N1(Id) $N2(Id) $N3(Id) $N4(Id)"
    } else {
        set Connectivities "$N4(Id) $N1(Id) $N2(Id) $N3(Id)"
    }
    
    return $Connectivities
}

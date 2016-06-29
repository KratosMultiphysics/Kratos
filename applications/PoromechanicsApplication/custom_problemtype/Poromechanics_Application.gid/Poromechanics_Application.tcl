## GiD events --------------------------------------------------------------------------------------------------------------------------------------------------

proc InitGIDProject { dir } {
        
    GiDMenu::Create "Poromechanics Application" PRE
    GiDMenu::InsertOption "Poromechanics Application" [list "Parts"] 0 PRE "GidOpenConditions \"Parts\"" "" ""
	GiDMenu::InsertOption "Poromechanics Application" [list "Dirichlet Boundary Conditions"] 1 PRE "GidOpenConditions \"Dirichlet_Boundary_Conditions\"" "" ""
	GiDMenu::InsertOption "Poromechanics Application" [list "Loads"] 2 PRE "GidOpenConditions \"Loads\"" "" ""
    GiDMenu::InsertOption "Poromechanics Application" [list "Project Parameters"] 3 PRE "GidOpenProblemData" "" ""
	GiDMenu::UpdateMenus
}

#-------------------------------------------------------------------------------

# Pass the path and the name of the problem to the Python script
proc BeforeRunCalculation { batfilename basename dir problemtypedir gidexe args } {
    
    set TableList [::Poromechanics_Application::WriteMDPA $basename $dir]
    
    ::Poromechanics_Application::WriteProjectParameters $basename $dir $gidexe $TableList
}

#-------------------------------------------------------------------------------

#proc AfterWriteCalcFileGIDProject { file error } {
#
#    WarnWin "Stopped after writing calculation file. Continue?"
#}


## Problemtype procedures --------------------------------------------------------------------------------------------------------------------------------------

namespace eval Poromechanics_Application {
    
}

#-------------------------------------------------------------------------------

proc Poromechanics_Application::WriteMDPA { basename dir } {
    
    ## Start MDPA file
    set filename [file join $dir ${basename}.dat]
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
            puts $varfile "  THICKNESS [lindex [lindex $Groups $i] 21]"
            puts $varfile "End Properties"
            puts $varfile ""
        } elseif {[lindex [lindex $Groups $i] 3]=="RestoreSimoJu3DLaw"} {
            incr PropertyId
            lappend PropertyList [lindex [lindex $Groups $i] 1] $PropertyId
            puts $varfile "Begin Properties $PropertyId"
            puts $varfile "  CONSTITUTIVE_LAW_NAME RestoreSimoJu3DLaw"
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
        } elseif {[lindex [lindex $Groups $i] 3]=="RestoreSimoJuPlaneStrain2DLaw" || [lindex [lindex $Groups $i] 3]=="RestoreSimoJuPlaneStress2DLaw"} {
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
            puts $varfile "  DAMAGE_THRESHOLD [lindex [lindex $Groups $i] 18]"
            puts $varfile "  STRENGTH_RATIO [lindex [lindex $Groups $i] 19]"
            puts $varfile "  FRACTURE_ENERGY [lindex [lindex $Groups $i] 20]"
            puts $varfile "  THICKNESS [lindex [lindex $Groups $i] 21]"
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
            puts $varfile "  RESIDUAL_STRESS [lindex [lindex $Groups $i] 17]"
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
            puts $varfile "  RESIDUAL_STRESS [lindex [lindex $Groups $i] 17]"
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
                        set Connectivities [::Poromechanics_Application::Triangle2D3Connectivities [lindex $Entities $j]]
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
                        set Connectivities [::Poromechanics_Application::Quadrilateral2D4Connectivities [lindex $Entities $j]]
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
                        set Connectivities [::Poromechanics_Application::Quadrilateral2D4Connectivities [lindex $Entities $j]]
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
                        set Connectivities [::Poromechanics_Application::Hexahedron3D8Connectivities [lindex $Entities $j]]
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
                        set Connectivities [::Poromechanics_Application::Triangle2D3Connectivities [lindex $Entities $j]]
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
                        set Connectivities [::Poromechanics_Application::Quadrilateral2D4Connectivities [lindex $Entities $j]]
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
                        set Connectivities [::Poromechanics_Application::Quadrilateral2D4Connectivities [lindex $Entities $j]]
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
                        set Connectivities [::Poromechanics_Application::Hexahedron3D8Connectivities [lindex $Entities $j]]
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
                    set Connectivities [::Poromechanics_Application::Triangle2D6Connectivities [lindex $Entities $j]]
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
                    set Connectivities [::Poromechanics_Application::Hexahedron3D8Connectivities [lindex $Entities $j]]
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
                    set Connectivities [::Poromechanics_Application::Tetrahedron3D10Connectivities [lindex $Entities $j]]
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
                    set Connectivities [::Poromechanics_Application::Hexahedron3D20Connectivities [lindex $Entities $j]]
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
                    set Connectivities [::Poromechanics_Application::Triangle2D6Connectivities [lindex $Entities $j]]
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
                    set Connectivities [::Poromechanics_Application::Quadrilateral2D9Connectivities [lindex $Entities $j]]
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
                    set Connectivities [::Poromechanics_Application::Tetrahedron3D10Connectivities [lindex $Entities $j]]
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
                    set Connectivities [::Poromechanics_Application::Hexahedron3D27Connectivities [lindex $Entities $j]]
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
                set Connectivities [::Poromechanics_Application::TriangleInterface2D4Connectivities [lindex $Entities $j]]
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
                    set Connectivities [::Poromechanics_Application::QuadrilateralInterface2D4Connectivities [lindex $Entities $j]]
                    puts $varfile "  [lindex $Entities $j]  $ElemsMat  $Connectivities"
                }
                puts $varfile "End Elements"
                puts $varfile ""
            } else {
                # UPwSmallStrainLinkInterfaceElement2D4N
                puts $varfile "Begin Elements UPwSmallStrainLinkInterfaceElement2D4N"
                for {set j 0} {$j < [llength $Entities]} {incr j} {
                    set Connectivities [::Poromechanics_Application::QuadrilateralInterface2D4Connectivities [lindex $Entities $j]]
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
                set Connectivities [::Poromechanics_Application::TetrahedronInterface3D6Connectivities [lindex $Entities $j]]
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
                    set Connectivities [::Poromechanics_Application::Triangle2D6Connectivities [lindex $Entities $j]]
                    puts $varfile "  [lindex $Entities $j]  $ElemsMat  $Connectivities"
                }
                puts $varfile "End Elements"
                puts $varfile ""
            } else {
                # UPwSmallStrainLinkInterfaceElement3D6N
                puts $varfile "Begin Elements UPwSmallStrainLinkInterfaceElement3D6N"
                for {set j 0} {$j < [llength $Entities]} {incr j} {
                    set Connectivities [::Poromechanics_Application::Triangle2D6Connectivities [lindex $Entities $j]]
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
                    set Connectivities [::Poromechanics_Application::HexaedronInterface3D8Connectivities [lindex $Entities $j]]
                    puts $varfile "  [lindex $Entities $j]  $ElemsMat  $Connectivities"
                }
                puts $varfile "End Elements"
                puts $varfile ""
            } else {
                # UPwSmallStrainLinkInterfaceElement3D8N
                puts $varfile "Begin Elements UPwSmallStrainLinkInterfaceElement3D8N"
                for {set j 0} {$j < [llength $Entities]} {incr j} {
                    set Connectivities [::Poromechanics_Application::HexaedronInterface3D8Connectivities [lindex $Entities $j]]
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
                set Connectivities [::Poromechanics_Application::Line2D2Connectivities [lindex $Entities $j]]
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
                set Connectivities [::Poromechanics_Application::TriangleInterface3D4Connectivities [lindex $Entities $j]]
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
                set Connectivities [::Poromechanics_Application::QuadrilateralInterface3D4Connectivities [lindex $Entities $j]]
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
                set Connectivities [::Poromechanics_Application::Line2D2Connectivities [lindex $Entities $j]]
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
                set Connectivities [::Poromechanics_Application::TriangleInterface3D4Connectivities [lindex $Entities $j]]
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
                set Connectivities [::Poromechanics_Application::QuadrilateralInterface3D4Connectivities [lindex $Entities $j]]
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

proc Poromechanics_Application::Triangle2D3Connectivities { ElemId } {
    
    set ElementInfo [GiD_Mesh get element $ElemId]
    #ElementInfo: <layer> <elemtype> <NumNodes> <N1> <N2> ...
    for {set i 0} {$i < 3} {incr i} {
        set N($i) [lindex $ElementInfo [expr { 3+$i }]]
    }
    set Connectivities "$N(0) $N(1) $N(2)"
    
    return $Connectivities
}

#-------------------------------------------------------------------------------

proc Poromechanics_Application::Quadrilateral2D4Connectivities { ElemId } {
    
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

proc Poromechanics_Application::Triangle2D6Connectivities { ElemId } {
    
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

proc Poromechanics_Application::Hexahedron3D8Connectivities { ElemId } {
    
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

proc Poromechanics_Application::Quadrilateral2D9Connectivities { ElemId } {
    
    set ElementInfo [GiD_Mesh get element $ElemId]
    #ElementInfo: <layer> <elemtype> <NumNodes> <N1> <N2> ...
    for {set i 0} {$i < 9} {incr i} {
        set N($i) [lindex $ElementInfo [expr { 3+$i }]]
    }

    set Connectivities "$N(0) $N(1) $N(2) $N(3) $N(4) $N(5) $N(6) $N(7) $N(8)"
    
    return $Connectivities
}

#-------------------------------------------------------------------------------

proc Poromechanics_Application::Tetrahedron3D10Connectivities { ElemId } {
        
    set ElementInfo [GiD_Mesh get element $ElemId]
    #ElementInfo: <layer> <elemtype> <NumNodes> <N1> <N2> ...
    for {set i 0} {$i < 10} {incr i} {
        set N($i) [lindex $ElementInfo [expr { 3+$i }]]
    }

    set Connectivities "$N(0) $N(1) $N(2) $N(3) $N(4) $N(5) $N(6) $N(7) $N(8) $N(9)"
    
    return $Connectivities
}

#-------------------------------------------------------------------------------

proc Poromechanics_Application::Hexahedron3D20Connectivities { ElemId } {
    
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

proc Poromechanics_Application::Hexahedron3D27Connectivities { ElemId } {
        
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

proc Poromechanics_Application::TriangleInterface2D4Connectivities { ElemId } {
    
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

proc Poromechanics_Application::QuadrilateralInterface2D4Connectivities { ElemId } {
    
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

proc Poromechanics_Application::TetrahedronInterface3D6Connectivities { ElemId } {
        
    set ElementInfo [GiD_Mesh get element $ElemId]
    #ElementInfo: <layer> <elemtype> <NumNodes> <N1> <N2> ...
    for {set i 0} {$i < 4} {incr i} {
        set N($i) [lindex $ElementInfo [expr { 3+$i }]]
    }
    
    set Connectivities "$N(0) $N(1) $N(2) $N(3) $N(1) $N(2)"
    
    return $Connectivities
}

#-------------------------------------------------------------------------------

proc Poromechanics_Application::HexaedronInterface3D8Connectivities { ElemId } {
    
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
    set NCoord [lindex [GiD_Info Coordinates $N5(Id)] 0]
    set N5(x) [lindex $NCoord 0]
    set N5(y) [lindex $NCoord 1]
    set N5(z) [lindex $NCoord 2]
    set NCoord [lindex [GiD_Info Coordinates $N6(Id)] 0]
    set N6(x) [lindex $NCoord 0]
    set N6(y) [lindex $NCoord 1]
    set N6(z) [lindex $NCoord 2]
    set NCoord [lindex [GiD_Info Coordinates $N7(Id)] 0]
    set N7(x) [lindex $NCoord 0]
    set N7(y) [lindex $NCoord 1]
    set N7(z) [lindex $NCoord 2]
    set NCoord [lindex [GiD_Info Coordinates $N8(Id)] 0]
    set N8(x) [lindex $NCoord 0]
    set N8(y) [lindex $NCoord 1]
    set N8(z) [lindex $NCoord 2]
    
    # Computing element lengths
    set lx [expr { sqrt( (0.25*($N2(x)+$N6(x)+$N3(x)+$N7(x)-$N1(x)-$N5(x)-$N4(x)-$N8(x)))**2 + (0.25*($N2(y)+$N6(y)+$N3(y)+$N7(y)-$N1(y)-$N5(y)-$N4(y)-$N8(y)))**2 + (0.25*($N2(z)+$N6(z)+$N3(z)+$N7(z)-$N1(z)-$N5(z)-$N4(z)-$N8(z)))**2 ) }]
    set ly [expr { sqrt( (0.25*($N3(x)+$N4(x)+$N7(x)+$N8(x)-$N1(x)-$N2(x)-$N5(x)-$N6(x)))**2 + (0.25*($N3(y)+$N4(y)+$N7(y)+$N8(y)-$N1(y)-$N2(y)-$N5(y)-$N6(y)))**2 + (0.25*($N3(z)+$N4(z)+$N7(z)+$N8(z)-$N1(z)-$N2(z)-$N5(z)-$N6(z)))**2 ) }]
    set lz [expr { sqrt( (0.25*($N5(x)+$N6(x)+$N7(x)+$N8(x)-$N1(x)-$N2(x)-$N3(x)-$N4(x)))**2 + (0.25*($N5(y)+$N6(y)+$N7(y)+$N8(y)-$N1(y)-$N2(y)-$N3(y)-$N4(y)))**2 + (0.25*($N5(z)+$N6(z)+$N7(z)+$N8(z)-$N1(z)-$N2(z)-$N3(z)-$N4(z)))**2 ) }]
    
    if {$lz <= $lx} {
        if {$lz <= $ly} {
            # lz <= lx && lz <= ly
            set Connectivities "$N1(Id) $N2(Id) $N3(Id) $N4(Id) $N5(Id) $N6(Id) $N7(Id) $N8(Id)"
        } else {
            # ly < lz <= lx
            set Connectivities "$N1(Id) $N4(Id) $N8(Id) $N5(Id) $N2(Id) $N3(Id) $N7(Id) $N6(Id)"
        }
    } elseif {$ly <= $lx} {
        # ly <= lx < lz
        set Connectivities "$N1(Id) $N4(Id) $N8(Id) $N5(Id) $N2(Id) $N3(Id) $N7(Id) $N6(Id)"
    } else {
        # lx < lz && lx < ly
        set Connectivities "$N1(Id) $N5(Id) $N6(Id) $N2(Id) $N4(Id) $N8(Id) $N7(Id) $N3(Id)"
    }
    
    return $Connectivities
}

#-------------------------------------------------------------------------------

proc Poromechanics_Application::Line2D2Connectivities { ElemId } {
    
    set ElementInfo [GiD_Mesh get element $ElemId]
    #ElementInfo: <layer> <elemtype> <NumNodes> <N1> <N2> ...
    for {set i 0} {$i < 2} {incr i} {
        set N($i) [lindex $ElementInfo [expr { 3+$i }]]
    }
    set Connectivities "$N(0) $N(1)"
    
    return $Connectivities
}

#-------------------------------------------------------------------------------

proc Poromechanics_Application::TriangleInterface3D4Connectivities { ElemId } {
    
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

proc Poromechanics_Application::QuadrilateralInterface3D4Connectivities { ElemId } {
    
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

#-------------------------------------------------------------------------------

proc Poromechanics_Application::WriteProjectParameters { basename dir gidexe TableList} {
    
    ## Start ProjectParameters.json file
    set filename [file join $dir ${basename}-1.dat]
    set varfile [open $filename w]
    
    puts $varfile "\{"
    
    ## problem_data
    puts $varfile "    \"problem_data\":             \{"
    puts $varfile "        \"problem_name\":      \"$basename\","
    puts $varfile "        \"model_part_name\":   \"PorousDomain\","
    puts $varfile "        \"domain_size\":       [GiD_AccessValue get gendata Domain_Size],"
    puts $varfile "        \"time_step\":         [GiD_AccessValue get gendata Delta_Time],"
    puts $varfile "        \"start_time\":        [GiD_AccessValue get gendata Start_Time],"
    puts $varfile "        \"end_time\":          [GiD_AccessValue get gendata End_Time],"
    puts $varfile "        \"echo_level\":        [GiD_AccessValue get gendata Echo_Level],"
    puts $varfile "        \"OMP_threads\":       [GiD_AccessValue get gendata OMP_Threads],"
    puts $varfile "        \"FIC_stabilization\": [GiD_AccessValue get gendata FIC_Stabilization],"
    if {[regexp -all {\\} $gidexe] > 0} {
        regsub -all {\\} $gidexe {\\\\} gidexe
        puts $varfile "        \"gid_path\":          \"$gidexe\""
    } else {
        puts $varfile "        \"gid_path\":          \"$gidexe\""
    }
    puts $varfile "    \},"
    
    ## solver_settings
    puts $varfile "    \"solver_settings\":          \{"
    puts $varfile "        \"solver_type\":                        \"poromechanics_U_Pw_solver\","
    puts $varfile "        \"model_import_settings\":              \{"
    puts $varfile "            \"input_type\":     \"mdpa\","
    puts $varfile "            \"input_filename\": \"$basename\""
    puts $varfile "        \},"
    puts $varfile "        \"buffer_size\":                        2,"
    puts $varfile "        \"echo_level\":                         [GiD_AccessValue get gendata Echo_Level],"
    puts $varfile "        \"reform_dofs_at_each_iteration\":      false,"
    puts $varfile "        \"compute_reactions\":                  [GiD_AccessValue get gendata Write_Reactions],"
    puts $varfile "        \"move_mesh_flag\":                     true,"
    puts $varfile "        \"solution_type\":                      \"[GiD_AccessValue get gendata Solution_Type]\","
    puts $varfile "        \"scheme_type\":                        \"[GiD_AccessValue get gendata Scheme_Type]\","
    puts $varfile "        \"newmark_beta\":                       [GiD_AccessValue get gendata Newmark_Beta],"
    puts $varfile "        \"newmark_gamma\":                      [GiD_AccessValue get gendata Newmark_Gamma],"
    puts $varfile "        \"newmark_theta\":                      [GiD_AccessValue get gendata Newmark_Theta],"
    puts $varfile "        \"rayleigh_m\":                         [GiD_AccessValue get gendata Rayleigh_Mass],"
    puts $varfile "        \"rayleigh_k\":                         [GiD_AccessValue get gendata Rayleigh_Stiffness],"
    puts $varfile "        \"strategy_type\":                      \"[GiD_AccessValue get gendata Strategy_Type]\","
    puts $varfile "        \"convergence_criterion\":              \"[GiD_AccessValue get gendata Convergence_Criterion]\","
    puts $varfile "        \"displacement_relative_tolerance\":    [GiD_AccessValue get gendata Displacement_Relative_Tolerance],"
    puts $varfile "        \"displacement_absolute_tolerance\":    [GiD_AccessValue get gendata Displacement_Absolute_Tolerance],"
    puts $varfile "        \"residual_relative_tolerance\":        [GiD_AccessValue get gendata Residual_Relative_Tolerance],"
    puts $varfile "        \"residual_absolute_tolerance\":        [GiD_AccessValue get gendata Residual_Absolute_Tolerance],"
    puts $varfile "        \"max_iteration\":                      [GiD_AccessValue get gendata Max_Iterations],"
    puts $varfile "        \"desired_iterations\":                 [GiD_AccessValue get gendata Desired_Iterations],"
    puts $varfile "        \"max_radius_factor\":                  [GiD_AccessValue get gendata Max_Radius_Factor],"
    puts $varfile "        \"min_radius_factor\":                  [GiD_AccessValue get gendata Min_Radius_Factor],"
    puts $varfile "        \"builder\":                            \"[GiD_AccessValue get gendata Builder]\","
    ## linear_solver_settings
    puts $varfile "        \"linear_solver_settings\":             \{"
    if {[GiD_AccessValue get gendata Solver_Type]=="AMGCL"} {
        puts $varfile "            \"solver_type\":                    \"AMGCL\","
        puts $varfile "            \"smoother_type\":                  \"ilu0\","
        puts $varfile "            \"krylov_type\":                    \"gmres\","
        puts $varfile "            \"coarsening_type\":                \"aggregation\","
        puts $varfile "            \"max_iteration\":                  100,"
        puts $varfile "            \"provide_coordinates\":            false,"
        puts $varfile "            \"gmres_krylov_space_dimension\":   100,"
        puts $varfile "            \"verbosity\":                      [GiD_AccessValue get gendata Verbosity],"
        puts $varfile "            \"tolerance\":                      1.0e-6,"
        puts $varfile "            \"scaling\":                        [GiD_AccessValue get gendata Scaling],"
        puts $varfile "            \"block_size\":                     1,"
        puts $varfile "            \"use_block_matrices_if_possible\": true,"
        puts $varfile "            \"coarse_enough\":                  5000"
    } elseif {[GiD_AccessValue get gendata Solver_Type]=="AMGCL_NS_Solver"} {
        puts $varfile "            \"solver_type\":                    \"AMGCL_NS_Solver\","
        puts $varfile "            \"krylov_type\":                    \"gmres\","
        puts $varfile "            \"velocity_block_preconditioner\":  \{"
        puts $varfile "                \"krylov_type\":         \"bicgstab\","
        puts $varfile "                \"tolerance\":           1.0e-3,"
        puts $varfile "                \"preconditioner_type\": \"spai0\","
        puts $varfile "                \"max_iteration\":       50"
        puts $varfile "            \},"
        puts $varfile "            \"pressure_block_preconditioner\":  \{"
        puts $varfile "                \"krylov_type\":         \"cg\","
        puts $varfile "                \"tolerance\":           1.0e-2,"
        puts $varfile "                \"preconditioner_type\": \"spai0\","
        puts $varfile "                \"max_iteration\":       50"
        puts $varfile "            \},"
        puts $varfile "            \"tolerance\":                      1.0e-9,"
        puts $varfile "            \"gmres_krylov_space_dimension\":   50,"
        puts $varfile "            \"coarsening_type\":                \"aggregation\","
        puts $varfile "            \"max_iteration\":                  50,"
        puts $varfile "            \"verbosity\":                      [GiD_AccessValue get gendata Verbosity],"
        puts $varfile "            \"scaling\":                        [GiD_AccessValue get gendata Scaling],"
        puts $varfile "            \"coarse_enough\":                  5000"
    } else {
        puts $varfile "            \"solver_type\":                    \"[GiD_AccessValue get gendata Solver_Type]\","
        puts $varfile "            \"tolerance\":                      1.0e-5,"
        puts $varfile "            \"max_iteration\":                  100,"
        puts $varfile "            \"preconditioner_type\":            \"ILU0Preconditioner\","
        puts $varfile "            \"scaling\":                        [GiD_AccessValue get gendata Scaling],"
        puts $varfile "            \"verbosity\":                      [GiD_AccessValue get gendata Verbosity]"
    }
    puts $varfile "        \},"
    ## problem_domain_sub_model_part_list
    set PutStrings \[
    # Body_Part
    set Groups [GiD_Info conditions Body_Part groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
        append PutStrings \" [lindex [lindex $Groups $i] 1] \" ,
    }
    # Interface_Part
    set Groups [GiD_Info conditions Interface_Part groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
        append PutStrings \" [lindex [lindex $Groups $i] 1] \" ,
    }
    set PutStrings [string trimright $PutStrings ,]
    append PutStrings \]
    puts $varfile "        \"problem_domain_sub_model_part_list\": $PutStrings,"
    ## processes_sub_model_part_list
    set PutStrings \[
    # Solid_Displacement
    set Groups [GiD_Info conditions Solid_Displacement groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
        append PutStrings \" [lindex [lindex $Groups $i] 1] \" ,
    }
    # Fluid_Pressure
    set Groups [GiD_Info conditions Fluid_Pressure groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
        append PutStrings \" [lindex [lindex $Groups $i] 1] \" ,
    }
    # Force
    set Groups [GiD_Info conditions Force groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
        append PutStrings \" [lindex [lindex $Groups $i] 1] \" ,
    }
    # Face_Load
    set Groups [GiD_Info conditions Face_Load groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
        append PutStrings \" [lindex [lindex $Groups $i] 1] \" ,
    }
    # Normal_Load
    set Groups [GiD_Info conditions Normal_Load groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
        append PutStrings \" [lindex [lindex $Groups $i] 1] \" ,
    }
    # Normal_Fluid_Flux
    set Groups [GiD_Info conditions Normal_Fluid_Flux groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
        append PutStrings \" [lindex [lindex $Groups $i] 1] \" ,
    }
    # Interface_Face_Load
    set Groups [GiD_Info conditions Interface_Face_Load groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
        append PutStrings \" [lindex [lindex $Groups $i] 1] \" ,
    }
    # Interface_Normal_Fluid_Flux
    set Groups [GiD_Info conditions Interface_Normal_Fluid_Flux groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
        append PutStrings \" [lindex [lindex $Groups $i] 1] \" ,
    }
    # Body_Acceleration
    set Groups [GiD_Info conditions Body_Acceleration groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
        append PutStrings \" [lindex [lindex $Groups $i] 1] \" ,
    }
    set PutStrings [string trimright $PutStrings ,]
    append PutStrings \]
    puts $varfile "        \"processes_sub_model_part_list\":      $PutStrings,"
    ## loads_sub_model_part_list
    set PutStrings \[
    set iGroup 0
    # Force
    set Groups [GiD_Info conditions Force groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
        incr iGroup
        append PutStrings \" [lindex [lindex $Groups $i] 1] \" ,
    }
    # Face_Load
    set Groups [GiD_Info conditions Face_Load groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
        incr iGroup
        append PutStrings \" [lindex [lindex $Groups $i] 1] \" ,
    }
    # Normal_Load
    set Groups [GiD_Info conditions Normal_Load groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
        incr iGroup
        append PutStrings \" [lindex [lindex $Groups $i] 1] \" ,
    }
    # Normal_Fluid_Flux
    set Groups [GiD_Info conditions Normal_Fluid_Flux groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
        incr iGroup
        append PutStrings \" [lindex [lindex $Groups $i] 1] \" ,
    }
    # Interface_Face_Load
    set Groups [GiD_Info conditions Interface_Face_Load groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
        incr iGroup
        append PutStrings \" [lindex [lindex $Groups $i] 1] \" ,
    }
    # Interface_Normal_Fluid_Flux
    set Groups [GiD_Info conditions Interface_Normal_Fluid_Flux groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
        incr iGroup
        append PutStrings \" [lindex [lindex $Groups $i] 1] \" ,
    }
    # Body_Acceleration
    set Groups [GiD_Info conditions Body_Acceleration groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
        incr iGroup
        append PutStrings \" [lindex [lindex $Groups $i] 1] \" ,
    }
    if {$iGroup > 0} {
        set PutStrings [string trimright $PutStrings ,]
    }
    append PutStrings \]
    puts $varfile "        \"loads_sub_model_part_list\":          $PutStrings,"
    ## loads_variable_list
    set PutStrings \[
    # Force
    set Groups [GiD_Info conditions Force groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
        append PutStrings \" FORCE \" ,
    }
    # Face_Load
    set Groups [GiD_Info conditions Face_Load groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
        append PutStrings \" FACE_LOAD \" ,
    }
    # Normal_Load
    set Groups [GiD_Info conditions Normal_Load groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
        append PutStrings \" NORMAL_CONTACT_STRESS \" ,
    }
    # Normal_Fluid_Flux
    set Groups [GiD_Info conditions Normal_Fluid_Flux groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
        append PutStrings \" NORMAL_FLUID_FLUX \" ,
    }
    # Interface_Face_Load
    set Groups [GiD_Info conditions Interface_Face_Load groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
        append PutStrings \" FACE_LOAD \" ,
    }
    # Interface_Normal_Fluid_Flux
    set Groups [GiD_Info conditions Interface_Normal_Fluid_Flux groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
        append PutStrings \" NORMAL_FLUID_FLUX \" ,
    }
    # Body_Acceleration
    set Groups [GiD_Info conditions Body_Acceleration groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
        append PutStrings \" VOLUME_ACCELERATION \" ,
    }
    if {$iGroup > 0} {
        set PutStrings [string trimright $PutStrings ,]
    }
    append PutStrings \]
    puts $varfile "        \"loads_variable_list\":                $PutStrings"
    puts $varfile "    \},"
    
    ## constraints_process_list
    set Groups [GiD_Info conditions Solid_Displacement groups]
    set NumGroups [llength $Groups]
    set Groups [GiD_Info conditions Fluid_Pressure groups]
    incr NumGroups [llength $Groups]
    if {$NumGroups > 0} {
        set iGroup 0
        puts $varfile "    \"constraints_process_list\": \[\{"
        # Solid_Displacement
        set Groups [GiD_Info conditions Solid_Displacement groups]
        for {set i 0} {$i < [llength $Groups]} {incr i} {
            set Entities [GiD_EntitiesGroups get [lindex [lindex $Groups $i] 1] volumes]
            if {[llength $Entities] > 0} {
                incr iGroup
                puts $varfile "        \"implemented_in_file\":   \"apply_constraint_vector_table_process\","
                puts $varfile "        \"implemented_in_module\": \"KratosMultiphysics.PoromechanicsApplication\","
                puts $varfile "        \"process_name\":          \"ApplyConstraintVectorTableProcess\","
                puts $varfile "        \"Parameters\":            \{"
                puts $varfile "            \"mesh_id\":         0,"
                puts $varfile "            \"model_part_name\": \"[lindex [lindex $Groups $i] 1]\","
                puts $varfile "            \"variable_name\":   \"DISPLACEMENT\","
                if {[lindex [lindex $Groups $i] 3] == 1} {
                    set PutStrings true
                } else {
                    set PutStrings false
                }
                if {[lindex [lindex $Groups $i] 8] == 1} {
                    append PutStrings , true
                } else {
                    append PutStrings , false
                }
                if {[lindex [lindex $Groups $i] 13] == 1} {
                    append PutStrings , true
                } else {
                    append PutStrings , false
                }
                puts $varfile "            \"active\":          \[$PutStrings\],"
                if {[lindex [lindex $Groups $i] 5] == 1} {
                    set PutStrings true
                } else {
                    set PutStrings false
                }
                if {[lindex [lindex $Groups $i] 10] == 1} {
                    append PutStrings , true
                } else {
                    append PutStrings , false
                }
                if {[lindex [lindex $Groups $i] 15] == 1} {
                    append PutStrings , true
                } else {
                    append PutStrings , false
                }
                puts $varfile "            \"is_fixed\":        \[$PutStrings\],"
                puts $varfile "            \"value\":           \[[lindex [lindex $Groups $i] 4],[lindex [lindex $Groups $i] 9],[lindex [lindex $Groups $i] 14]\],"
                if {[GiD_AccessValue get gendata Strategy_Type]=="Arc-Length"} {
                    puts $varfile "            \"table\":           \[0,0,0\]"
                } else {
                    set SearchInList [lsearch $TableList [lindex [lindex $Groups $i] 1]*]
                    set AuxList [lindex $TableList [expr { $SearchInList+1 }]]
                    puts $varfile "            \"table\":           \[[lindex $AuxList 0],[lindex $AuxList 1],[lindex $AuxList 2]\]"
                }
                puts $varfile "        \}"
                if {$iGroup < $NumGroups} {
                    puts $varfile "    \},\{"
                } else {
                    puts $varfile "    \}\],"
                }
            }
        }
        for {set i 0} {$i < [llength $Groups]} {incr i} {
            set Entities [GiD_EntitiesGroups get [lindex [lindex $Groups $i] 1] surfaces]
            if {[llength $Entities] > 0} {
                incr iGroup
                puts $varfile "        \"implemented_in_file\":   \"apply_constraint_vector_table_process\","
                puts $varfile "        \"implemented_in_module\": \"KratosMultiphysics.PoromechanicsApplication\","
                puts $varfile "        \"process_name\":          \"ApplyConstraintVectorTableProcess\","
                puts $varfile "        \"Parameters\":            \{"
                puts $varfile "            \"mesh_id\":         0,"
                puts $varfile "            \"model_part_name\": \"[lindex [lindex $Groups $i] 1]\","
                puts $varfile "            \"variable_name\":   \"DISPLACEMENT\","
                if {[lindex [lindex $Groups $i] 3] == 1} {
                    set PutStrings true
                } else {
                    set PutStrings false
                }
                if {[lindex [lindex $Groups $i] 8] == 1} {
                    append PutStrings , true
                } else {
                    append PutStrings , false
                }
                if {[lindex [lindex $Groups $i] 13] == 1} {
                    append PutStrings , true
                } else {
                    append PutStrings , false
                }
                puts $varfile "            \"active\":          \[$PutStrings\],"
                if {[lindex [lindex $Groups $i] 5] == 1} {
                    set PutStrings true
                } else {
                    set PutStrings false
                }
                if {[lindex [lindex $Groups $i] 10] == 1} {
                    append PutStrings , true
                } else {
                    append PutStrings , false
                }
                if {[lindex [lindex $Groups $i] 15] == 1} {
                    append PutStrings , true
                } else {
                    append PutStrings , false
                }
                puts $varfile "            \"is_fixed\":        \[$PutStrings\],"
                puts $varfile "            \"value\":           \[[lindex [lindex $Groups $i] 4],[lindex [lindex $Groups $i] 9],[lindex [lindex $Groups $i] 14]\],"
                if {[GiD_AccessValue get gendata Strategy_Type]=="Arc-Length"} {
                    puts $varfile "            \"table\":           \[0,0,0\]"
                } else {
                    set SearchInList [lsearch $TableList [lindex [lindex $Groups $i] 1]*]
                    set AuxList [lindex $TableList [expr { $SearchInList+1 }]]
                    puts $varfile "            \"table\":           \[[lindex $AuxList 0],[lindex $AuxList 1],[lindex $AuxList 2]\]"
                }
                puts $varfile "        \}"
                if {$iGroup < $NumGroups} {
                    puts $varfile "    \},\{"
                } else {
                    puts $varfile "    \}\],"
                }
            }
        }
        for {set i 0} {$i < [llength $Groups]} {incr i} {
            set Entities [GiD_EntitiesGroups get [lindex [lindex $Groups $i] 1] lines]
            if {[llength $Entities] > 0} {
                incr iGroup
                puts $varfile "        \"implemented_in_file\":   \"apply_constraint_vector_table_process\","
                puts $varfile "        \"implemented_in_module\": \"KratosMultiphysics.PoromechanicsApplication\","
                puts $varfile "        \"process_name\":          \"ApplyConstraintVectorTableProcess\","
                puts $varfile "        \"Parameters\":            \{"
                puts $varfile "            \"mesh_id\":         0,"
                puts $varfile "            \"model_part_name\": \"[lindex [lindex $Groups $i] 1]\","
                puts $varfile "            \"variable_name\":   \"DISPLACEMENT\","
                if {[lindex [lindex $Groups $i] 3] == 1} {
                    set PutStrings true
                } else {
                    set PutStrings false
                }
                if {[lindex [lindex $Groups $i] 8] == 1} {
                    append PutStrings , true
                } else {
                    append PutStrings , false
                }
                if {[lindex [lindex $Groups $i] 13] == 1} {
                    append PutStrings , true
                } else {
                    append PutStrings , false
                }
                puts $varfile "            \"active\":          \[$PutStrings\],"
                if {[lindex [lindex $Groups $i] 5] == 1} {
                    set PutStrings true
                } else {
                    set PutStrings false
                }
                if {[lindex [lindex $Groups $i] 10] == 1} {
                    append PutStrings , true
                } else {
                    append PutStrings , false
                }
                if {[lindex [lindex $Groups $i] 15] == 1} {
                    append PutStrings , true
                } else {
                    append PutStrings , false
                }
                puts $varfile "            \"is_fixed\":        \[$PutStrings\],"
                puts $varfile "            \"value\":           \[[lindex [lindex $Groups $i] 4],[lindex [lindex $Groups $i] 9],[lindex [lindex $Groups $i] 14]\],"
                if {[GiD_AccessValue get gendata Strategy_Type]=="Arc-Length"} {
                    puts $varfile "            \"table\":           \[0,0,0\]"
                } else {
                    set SearchInList [lsearch $TableList [lindex [lindex $Groups $i] 1]*]
                    set AuxList [lindex $TableList [expr { $SearchInList+1 }]]
                    puts $varfile "            \"table\":           \[[lindex $AuxList 0],[lindex $AuxList 1],[lindex $AuxList 2]\]"
                }
                puts $varfile "        \}"
                if {$iGroup < $NumGroups} {
                    puts $varfile "    \},\{"
                } else {
                    puts $varfile "    \}\],"
                }
            }
        }
        for {set i 0} {$i < [llength $Groups]} {incr i} {
            set Entities [GiD_EntitiesGroups get [lindex [lindex $Groups $i] 1] points]
            if {[llength $Entities] > 0} {
                incr iGroup
                puts $varfile "        \"implemented_in_file\":   \"apply_constraint_vector_table_process\","
                puts $varfile "        \"implemented_in_module\": \"KratosMultiphysics.PoromechanicsApplication\","
                puts $varfile "        \"process_name\":          \"ApplyConstraintVectorTableProcess\","
                puts $varfile "        \"Parameters\":            \{"
                puts $varfile "            \"mesh_id\":         0,"
                puts $varfile "            \"model_part_name\": \"[lindex [lindex $Groups $i] 1]\","
                puts $varfile "            \"variable_name\":   \"DISPLACEMENT\","
                if {[lindex [lindex $Groups $i] 3] == 1} {
                    set PutStrings true
                } else {
                    set PutStrings false
                }
                if {[lindex [lindex $Groups $i] 8] == 1} {
                    append PutStrings , true
                } else {
                    append PutStrings , false
                }
                if {[lindex [lindex $Groups $i] 13] == 1} {
                    append PutStrings , true
                } else {
                    append PutStrings , false
                }
                puts $varfile "            \"active\":          \[$PutStrings\],"
                if {[lindex [lindex $Groups $i] 5] == 1} {
                    set PutStrings true
                } else {
                    set PutStrings false
                }
                if {[lindex [lindex $Groups $i] 10] == 1} {
                    append PutStrings , true
                } else {
                    append PutStrings , false
                }
                if {[lindex [lindex $Groups $i] 15] == 1} {
                    append PutStrings , true
                } else {
                    append PutStrings , false
                }
                puts $varfile "            \"is_fixed\":        \[$PutStrings\],"
                puts $varfile "            \"value\":           \[[lindex [lindex $Groups $i] 4],[lindex [lindex $Groups $i] 9],[lindex [lindex $Groups $i] 14]\],"
                if {[GiD_AccessValue get gendata Strategy_Type]=="Arc-Length"} {
                    puts $varfile "            \"table\":           \[0,0,0\]"
                } else {
                    set SearchInList [lsearch $TableList [lindex [lindex $Groups $i] 1]*]
                    set AuxList [lindex $TableList [expr { $SearchInList+1 }]]
                    puts $varfile "            \"table\":           \[[lindex $AuxList 0],[lindex $AuxList 1],[lindex $AuxList 2]\]"
                }
                puts $varfile "        \}"
                if {$iGroup < $NumGroups} {
                    puts $varfile "    \},\{"
                } else {
                    puts $varfile "    \}\],"
                }
            }
        }
        # Fluid_Pressure
        set Groups [GiD_Info conditions Fluid_Pressure groups]
        for {set i 0} {$i < [llength $Groups]} {incr i} {
            set Entities [GiD_EntitiesGroups get [lindex [lindex $Groups $i] 1] volumes]
            if {[llength $Entities] > 0} {
                incr iGroup
                puts $varfile "        \"implemented_in_file\":   \"apply_pore_pressure_table_process\","
                puts $varfile "        \"implemented_in_module\": \"KratosMultiphysics.PoromechanicsApplication\","
                puts $varfile "        \"process_name\":          \"ApplyPorePressureTableProcess\","
                puts $varfile "        \"Parameters\":            \{"
                puts $varfile "            \"mesh_id\":              0,"
                puts $varfile "            \"model_part_name\":      \"[lindex [lindex $Groups $i] 1]\","
                puts $varfile "            \"variable_name\":        \"WATER_PRESSURE\","
                if {[lindex [lindex $Groups $i] 8] == 1} {
                    set PutStrings true
                } else {
                    set PutStrings false
                }
                puts $varfile "            \"is_fixed\":             $PutStrings,"
                puts $varfile "            \"value\":                [lindex [lindex $Groups $i] 4],"
                if {[GiD_AccessValue get gendata Strategy_Type]=="Arc-Length"} {
                    puts $varfile "            \"table\":                0,"
                } else {
                    set SearchInList [lsearch $TableList [lindex [lindex $Groups $i] 1]*]
                    set AuxList [lindex $TableList [expr { $SearchInList+1 }]]
                    puts $varfile "            \"table\":                $AuxList,"
                }
                if {[lindex [lindex $Groups $i] 3] == "Hydrostatic"} {
                    set PutStrings true
                } else {
                    set PutStrings false
                }
                puts $varfile "            \"hydrostatic\":          $PutStrings,"
                if {[lindex [lindex $Groups $i] 5] == "Y"} {
                    set PutStrings 2
                } elseif {[lindex [lindex $Groups $i] 5] == "Z"} {
                    set PutStrings 3
                } else {
                    set PutStrings 1
                }
                puts $varfile "            \"gravity_direction\":    $PutStrings,"
                puts $varfile "            \"reference_coordinate\": [lindex [lindex $Groups $i] 6],"
                puts $varfile "            \"specific_weight\":      [lindex [lindex $Groups $i] 7]"
                puts $varfile "        \}"
                if {$iGroup < $NumGroups} {
                    puts $varfile "    \},\{"
                } else {
                    puts $varfile "    \}\],"
                }
            }
        }
        for {set i 0} {$i < [llength $Groups]} {incr i} {
            set Entities [GiD_EntitiesGroups get [lindex [lindex $Groups $i] 1] surfaces]
            if {[llength $Entities] > 0} {
                incr iGroup
                puts $varfile "        \"implemented_in_file\":   \"apply_pore_pressure_table_process\","
                puts $varfile "        \"implemented_in_module\": \"KratosMultiphysics.PoromechanicsApplication\","
                puts $varfile "        \"process_name\":          \"ApplyPorePressureTableProcess\","
                puts $varfile "        \"Parameters\":            \{"
                puts $varfile "            \"mesh_id\":              0,"
                puts $varfile "            \"model_part_name\":      \"[lindex [lindex $Groups $i] 1]\","
                puts $varfile "            \"variable_name\":        \"WATER_PRESSURE\","
                if {[lindex [lindex $Groups $i] 8] == 1} {
                    set PutStrings true
                } else {
                    set PutStrings false
                }
                puts $varfile "            \"is_fixed\":             $PutStrings,"
                puts $varfile "            \"value\":                [lindex [lindex $Groups $i] 4],"
                if {[GiD_AccessValue get gendata Strategy_Type]=="Arc-Length"} {
                    puts $varfile "            \"table\":                0,"
                } else {
                    set SearchInList [lsearch $TableList [lindex [lindex $Groups $i] 1]*]
                    set AuxList [lindex $TableList [expr { $SearchInList+1 }]]
                    puts $varfile "            \"table\":                $AuxList,"
                }
                if {[lindex [lindex $Groups $i] 3] == "Hydrostatic"} {
                    set PutStrings true
                } else {
                    set PutStrings false
                }
                puts $varfile "            \"hydrostatic\":          $PutStrings,"
                if {[lindex [lindex $Groups $i] 5] == "Y"} {
                    set PutStrings 2
                } elseif {[lindex [lindex $Groups $i] 5] == "Z"} {
                    set PutStrings 3
                } else {
                    set PutStrings 1
                }
                puts $varfile "            \"gravity_direction\":    $PutStrings,"
                puts $varfile "            \"reference_coordinate\": [lindex [lindex $Groups $i] 6],"
                puts $varfile "            \"specific_weight\":      [lindex [lindex $Groups $i] 7]"
                puts $varfile "        \}"
                if {$iGroup < $NumGroups} {
                    puts $varfile "    \},\{"
                } else {
                    puts $varfile "    \}\],"
                }
            }
        }
        for {set i 0} {$i < [llength $Groups]} {incr i} {
            set Entities [GiD_EntitiesGroups get [lindex [lindex $Groups $i] 1] lines]
            if {[llength $Entities] > 0} {
                incr iGroup
                puts $varfile "        \"implemented_in_file\":   \"apply_pore_pressure_table_process\","
                puts $varfile "        \"implemented_in_module\": \"KratosMultiphysics.PoromechanicsApplication\","
                puts $varfile "        \"process_name\":          \"ApplyPorePressureTableProcess\","
                puts $varfile "        \"Parameters\":            \{"
                puts $varfile "            \"mesh_id\":              0,"
                puts $varfile "            \"model_part_name\":      \"[lindex [lindex $Groups $i] 1]\","
                puts $varfile "            \"variable_name\":        \"WATER_PRESSURE\","
                if {[lindex [lindex $Groups $i] 8] == 1} {
                    set PutStrings true
                } else {
                    set PutStrings false
                }
                puts $varfile "            \"is_fixed\":             $PutStrings,"
                puts $varfile "            \"value\":                [lindex [lindex $Groups $i] 4],"
                if {[GiD_AccessValue get gendata Strategy_Type]=="Arc-Length"} {
                    puts $varfile "            \"table\":                0,"
                } else {
                    set SearchInList [lsearch $TableList [lindex [lindex $Groups $i] 1]*]
                    set AuxList [lindex $TableList [expr { $SearchInList+1 }]]
                    puts $varfile "            \"table\":                $AuxList,"
                }
                if {[lindex [lindex $Groups $i] 3] == "Hydrostatic"} {
                    set PutStrings true
                } else {
                    set PutStrings false
                }
                puts $varfile "            \"hydrostatic\":          $PutStrings,"
                if {[lindex [lindex $Groups $i] 5] == "Y"} {
                    set PutStrings 2
                } elseif {[lindex [lindex $Groups $i] 5] == "Z"} {
                    set PutStrings 3
                } else {
                    set PutStrings 1
                }
                puts $varfile "            \"gravity_direction\":    $PutStrings,"
                puts $varfile "            \"reference_coordinate\": [lindex [lindex $Groups $i] 6],"
                puts $varfile "            \"specific_weight\":      [lindex [lindex $Groups $i] 7]"
                puts $varfile "        \}"
                if {$iGroup < $NumGroups} {
                    puts $varfile "    \},\{"
                } else {
                    puts $varfile "    \}\],"
                }
            }
        }
        for {set i 0} {$i < [llength $Groups]} {incr i} {
            set Entities [GiD_EntitiesGroups get [lindex [lindex $Groups $i] 1] points]
            if {[llength $Entities] > 0} {
                incr iGroup
                puts $varfile "        \"implemented_in_file\":   \"apply_pore_pressure_table_process\","
                puts $varfile "        \"implemented_in_module\": \"KratosMultiphysics.PoromechanicsApplication\","
                puts $varfile "        \"process_name\":          \"ApplyPorePressureTableProcess\","
                puts $varfile "        \"Parameters\":            \{"
                puts $varfile "            \"mesh_id\":              0,"
                puts $varfile "            \"model_part_name\":      \"[lindex [lindex $Groups $i] 1]\","
                puts $varfile "            \"variable_name\":        \"WATER_PRESSURE\","
                if {[lindex [lindex $Groups $i] 8] == 1} {
                    set PutStrings true
                } else {
                    set PutStrings false
                }
                puts $varfile "            \"is_fixed\":             $PutStrings,"
                puts $varfile "            \"value\":                [lindex [lindex $Groups $i] 4],"
                if {[GiD_AccessValue get gendata Strategy_Type]=="Arc-Length"} {
                    puts $varfile "            \"table\":                0,"
                } else {
                    set SearchInList [lsearch $TableList [lindex [lindex $Groups $i] 1]*]
                    set AuxList [lindex $TableList [expr { $SearchInList+1 }]]
                    puts $varfile "            \"table\":                $AuxList,"
                }
                if {[lindex [lindex $Groups $i] 3] == "Hydrostatic"} {
                    set PutStrings true
                } else {
                    set PutStrings false
                }
                puts $varfile "            \"hydrostatic\":          $PutStrings,"
                if {[lindex [lindex $Groups $i] 5] == "Y"} {
                    set PutStrings 2
                } elseif {[lindex [lindex $Groups $i] 5] == "Z"} {
                    set PutStrings 3
                } else {
                    set PutStrings 1
                }
                puts $varfile "            \"gravity_direction\":    $PutStrings,"
                puts $varfile "            \"reference_coordinate\": [lindex [lindex $Groups $i] 6],"
                puts $varfile "            \"specific_weight\":      [lindex [lindex $Groups $i] 7]"
                puts $varfile "        \}"
                if {$iGroup < $NumGroups} {
                    puts $varfile "    \},\{"
                } else {
                    puts $varfile "    \}\],"
                }
            }
        }

    } else {
        puts $varfile "    \"constraints_process_list\": \[\],"
    }
    
    ## loads_process_list
    set Groups [GiD_Info conditions Force groups]
    set NumGroups [llength $Groups]
    set Groups [GiD_Info conditions Face_Load groups]
    incr NumGroups [llength $Groups]
    set Groups [GiD_Info conditions Normal_Load groups]
    incr NumGroups [llength $Groups]
    set Groups [GiD_Info conditions Normal_Fluid_Flux groups]
    incr NumGroups [llength $Groups]
    set Groups [GiD_Info conditions Interface_Face_Load groups]
    incr NumGroups [llength $Groups]
    set Groups [GiD_Info conditions Interface_Normal_Fluid_Flux groups]
    incr NumGroups [llength $Groups]
    set Groups [GiD_Info conditions Body_Acceleration groups]
    incr NumGroups [llength $Groups]
    if {$NumGroups > 0} {
        set iGroup 0
        puts $varfile "    \"loads_process_list\":       \[\{"
        # Force
        set Groups [GiD_Info conditions Force groups]
        for {set i 0} {$i < [llength $Groups]} {incr i} {
            incr iGroup
            puts $varfile "        \"implemented_in_file\":   \"apply_load_vector_table_process\","
            puts $varfile "        \"implemented_in_module\": \"KratosMultiphysics.PoromechanicsApplication\","
            puts $varfile "        \"process_name\":          \"ApplyLoadVectorTableProcess\","
            puts $varfile "        \"Parameters\":            \{"
            puts $varfile "            \"mesh_id\":         0,"
            puts $varfile "            \"model_part_name\": \"[lindex [lindex $Groups $i] 1]\","
            puts $varfile "            \"variable_name\":   \"FORCE\","
            if {[lindex [lindex $Groups $i] 3] == 1} {
                set PutStrings true
            } else {
                set PutStrings false
            }
            if {[lindex [lindex $Groups $i] 7] == 1} {
                append PutStrings , true
            } else {
                append PutStrings , false
            }
            if {[lindex [lindex $Groups $i] 11] == 1} {
                append PutStrings , true
            } else {
                append PutStrings , false
            }
            puts $varfile "            \"active\":          \[$PutStrings\],"
            puts $varfile "            \"value\":           \[[lindex [lindex $Groups $i] 4],[lindex [lindex $Groups $i] 8],[lindex [lindex $Groups $i] 12]\],"
            if {[GiD_AccessValue get gendata Strategy_Type]=="Arc-Length"} {
                puts $varfile "            \"table\":           \[0,0,0\]"
            } else {
                set SearchInList [lsearch $TableList [lindex [lindex $Groups $i] 1]*]
                set AuxList [lindex $TableList [expr { $SearchInList+1 }]]
                puts $varfile "            \"table\":           \[[lindex $AuxList 0],[lindex $AuxList 1],[lindex $AuxList 2]\]"
            }
            puts $varfile "        \}"
            if {$iGroup < $NumGroups} {
                puts $varfile "    \},\{"
            } else {
                puts $varfile "    \}\],"
            }
        }
        # Face_Load
        set Groups [GiD_Info conditions Face_Load groups]
        for {set i 0} {$i < [llength $Groups]} {incr i} {
            incr iGroup
            puts $varfile "        \"implemented_in_file\":   \"apply_load_vector_table_process\","
            puts $varfile "        \"implemented_in_module\": \"KratosMultiphysics.PoromechanicsApplication\","
            puts $varfile "        \"process_name\":          \"ApplyLoadVectorTableProcess\","
            puts $varfile "        \"Parameters\":            \{"
            puts $varfile "            \"mesh_id\":         0,"
            puts $varfile "            \"model_part_name\": \"[lindex [lindex $Groups $i] 1]\","
            puts $varfile "            \"variable_name\":   \"FACE_LOAD\","
            if {[lindex [lindex $Groups $i] 3] == 1} {
                set PutStrings true
            } else {
                set PutStrings false
            }
            if {[lindex [lindex $Groups $i] 7] == 1} {
                append PutStrings , true
            } else {
                append PutStrings , false
            }
            if {[lindex [lindex $Groups $i] 11] == 1} {
                append PutStrings , true
            } else {
                append PutStrings , false
            }
            puts $varfile "            \"active\":          \[$PutStrings\],"
            puts $varfile "            \"value\":           \[[lindex [lindex $Groups $i] 4],[lindex [lindex $Groups $i] 8],[lindex [lindex $Groups $i] 12]\],"
            if {[GiD_AccessValue get gendata Strategy_Type]=="Arc-Length"} {
                puts $varfile "            \"table\":           \[0,0,0\]"
            } else {
                set SearchInList [lsearch $TableList [lindex [lindex $Groups $i] 1]*]
                set AuxList [lindex $TableList [expr { $SearchInList+1 }]]
                puts $varfile "            \"table\":           \[[lindex $AuxList 0],[lindex $AuxList 1],[lindex $AuxList 2]\]"
            }
            puts $varfile "        \}"
            if {$iGroup < $NumGroups} {
                puts $varfile "    \},\{"
            } else {
                puts $varfile "    \}\],"
            }
        }
        # Normal_Load
        set Groups [GiD_Info conditions Normal_Load groups]
        for {set i 0} {$i < [llength $Groups]} {incr i} {
            incr iGroup
            puts $varfile "        \"implemented_in_file\":   \"apply_normal_load_table_process\","
            puts $varfile "        \"implemented_in_module\": \"KratosMultiphysics.PoromechanicsApplication\","
            puts $varfile "        \"process_name\":          \"ApplyNormalLoadTableProcess\","
            puts $varfile "        \"Parameters\":            \{"
            puts $varfile "            \"mesh_id\":              0,"
            puts $varfile "            \"model_part_name\":      \"[lindex [lindex $Groups $i] 1]\","
            puts $varfile "            \"variable_name\":        \"NORMAL_CONTACT_STRESS\","
            if {[lindex [lindex $Groups $i] 3] == 1} {
                set PutStrings true
            } else {
                set PutStrings false
            }
            if {[lindex [lindex $Groups $i] 11] == 1} {
                append PutStrings , true
            } else {
                append PutStrings , false
            }
            puts $varfile "            \"active\":               \[$PutStrings\],"
            puts $varfile "            \"value\":                \[[lindex [lindex $Groups $i] 5],[lindex [lindex $Groups $i] 12]\],"
            if {[GiD_AccessValue get gendata Strategy_Type]=="Arc-Length"} {
                puts $varfile "            \"table\":                \[0,0\],"
            } else {
                set SearchInList [lsearch $TableList [lindex [lindex $Groups $i] 1]*]
                set AuxList [lindex $TableList [expr { $SearchInList+1 }]]
                puts $varfile "            \"table\":                \[[lindex $AuxList 0],[lindex $AuxList 1]\],"
            }
            if {[lindex [lindex $Groups $i] 4] == "Hydrostatic"} {
                set PutStrings true
            } else {
                set PutStrings false
            }
            puts $varfile "            \"hydrostatic\":          $PutStrings,"
            if {[lindex [lindex $Groups $i] 6] == "Y"} {
                set PutStrings 2
            } elseif {[lindex [lindex $Groups $i] 6] == "Z"} {
                set PutStrings 3
            } else {
                set PutStrings 1
            }
            puts $varfile "            \"gravity_direction\":    $PutStrings,"
            puts $varfile "            \"reference_coordinate\": [lindex [lindex $Groups $i] 7],"
            puts $varfile "            \"specific_weight\":      [lindex [lindex $Groups $i] 8]"
            puts $varfile "        \}"
            if {$iGroup < $NumGroups} {
                puts $varfile "    \},\{"
            } else {
                puts $varfile "    \}\],"
            }
        }
        # Normal_Fluid_Flux
        set Groups [GiD_Info conditions Normal_Fluid_Flux groups]
        for {set i 0} {$i < [llength $Groups]} {incr i} {
            incr iGroup
            puts $varfile "        \"implemented_in_file\":   \"apply_load_scalar_table_process\","
            puts $varfile "        \"implemented_in_module\": \"KratosMultiphysics.PoromechanicsApplication\","
            puts $varfile "        \"process_name\":          \"ApplyLoadScalarTableProcess\","
            puts $varfile "        \"Parameters\":            \{"
            puts $varfile "            \"mesh_id\":         0,"
            puts $varfile "            \"model_part_name\": \"[lindex [lindex $Groups $i] 1]\","
            puts $varfile "            \"variable_name\":   \"NORMAL_FLUID_FLUX\","
            puts $varfile "            \"value\":           [lindex [lindex $Groups $i] 3],"
            if {[GiD_AccessValue get gendata Strategy_Type]=="Arc-Length"} {
                puts $varfile "            \"table\":           0"
            } else {
                set SearchInList [lsearch $TableList [lindex [lindex $Groups $i] 1]*]
                set AuxList [lindex $TableList [expr { $SearchInList+1 }]]
                puts $varfile "            \"table\":           $AuxList"
            }
            puts $varfile "        \}"
            if {$iGroup < $NumGroups} {
                puts $varfile "    \},\{"
            } else {
                puts $varfile "    \}\],"
            }
        }
        # Interface_Face_Load
        set Groups [GiD_Info conditions Interface_Face_Load groups]
        for {set i 0} {$i < [llength $Groups]} {incr i} {
            incr iGroup
            puts $varfile "        \"implemented_in_file\":   \"apply_load_vector_table_process\","
            puts $varfile "        \"implemented_in_module\": \"KratosMultiphysics.PoromechanicsApplication\","
            puts $varfile "        \"process_name\":          \"ApplyLoadVectorTableProcess\","
            puts $varfile "        \"Parameters\":            \{"
            puts $varfile "            \"mesh_id\":         0,"
            puts $varfile "            \"model_part_name\": \"[lindex [lindex $Groups $i] 1]\","
            puts $varfile "            \"variable_name\":   \"FACE_LOAD\","
            if {[lindex [lindex $Groups $i] 3] == 1} {
                set PutStrings true
            } else {
                set PutStrings false
            }
            if {[lindex [lindex $Groups $i] 7] == 1} {
                append PutStrings , true
            } else {
                append PutStrings , false
            }
            if {[lindex [lindex $Groups $i] 11] == 1} {
                append PutStrings , true
            } else {
                append PutStrings , false
            }
            puts $varfile "            \"active\":          \[$PutStrings\],"
            puts $varfile "            \"value\":           \[[lindex [lindex $Groups $i] 4],[lindex [lindex $Groups $i] 8],[lindex [lindex $Groups $i] 12]\],"
            if {[GiD_AccessValue get gendata Strategy_Type]=="Arc-Length"} {
                puts $varfile "            \"table\":           \[0,0,0\]"
            } else {
                set SearchInList [lsearch $TableList [lindex [lindex $Groups $i] 1]*]
                set AuxList [lindex $TableList [expr { $SearchInList+1 }]]
                puts $varfile "            \"table\":           \[[lindex $AuxList 0],[lindex $AuxList 1],[lindex $AuxList 2]\]"
            }
            puts $varfile "        \}"
            if {$iGroup < $NumGroups} {
                puts $varfile "    \},\{"
            } else {
                puts $varfile "    \}\],"
            }
        }
        # Interface_Normal_Fluid_Flux
        set Groups [GiD_Info conditions Interface_Normal_Fluid_Flux groups]
        for {set i 0} {$i < [llength $Groups]} {incr i} {
            incr iGroup
            puts $varfile "        \"implemented_in_file\":   \"apply_load_scalar_table_process\","
            puts $varfile "        \"implemented_in_module\": \"KratosMultiphysics.PoromechanicsApplication\","
            puts $varfile "        \"process_name\":          \"ApplyLoadScalarTableProcess\","
            puts $varfile "        \"Parameters\":            \{"
            puts $varfile "            \"mesh_id\":         0,"
            puts $varfile "            \"model_part_name\": \"[lindex [lindex $Groups $i] 1]\","
            puts $varfile "            \"variable_name\":   \"NORMAL_FLUID_FLUX\","
            puts $varfile "            \"value\":           [lindex [lindex $Groups $i] 3],"
            if {[GiD_AccessValue get gendata Strategy_Type]=="Arc-Length"} {
                puts $varfile "            \"table\":           0"
            } else {
                set SearchInList [lsearch $TableList [lindex [lindex $Groups $i] 1]*]
                set AuxList [lindex $TableList [expr { $SearchInList+1 }]]
                puts $varfile "            \"table\":           $AuxList"
            }
            puts $varfile "        \}"
            if {$iGroup < $NumGroups} {
                puts $varfile "    \},\{"
            } else {
                puts $varfile "    \}\],"
            }
        }
        # Body_Acceleration
        set Groups [GiD_Info conditions Body_Acceleration groups]
        for {set i 0} {$i < [llength $Groups]} {incr i} {
            incr iGroup
            puts $varfile "        \"implemented_in_file\":   \"apply_load_vector_table_process\","
            puts $varfile "        \"implemented_in_module\": \"KratosMultiphysics.PoromechanicsApplication\","
            puts $varfile "        \"process_name\":          \"ApplyLoadVectorTableProcess\","
            puts $varfile "        \"Parameters\":            \{"
            puts $varfile "            \"mesh_id\":         0,"
            puts $varfile "            \"model_part_name\": \"[lindex [lindex $Groups $i] 1]\","
            puts $varfile "            \"variable_name\":   \"VOLUME_ACCELERATION\","
            if {[lindex [lindex $Groups $i] 3] == 1} {
                set PutStrings true
            } else {
                set PutStrings false
            }
            if {[lindex [lindex $Groups $i] 7] == 1} {
                append PutStrings , true
            } else {
                append PutStrings , false
            }
            if {[lindex [lindex $Groups $i] 11] == 1} {
                append PutStrings , true
            } else {
                append PutStrings , false
            }
            puts $varfile "            \"active\":          \[$PutStrings\],"
            puts $varfile "            \"value\":           \[[lindex [lindex $Groups $i] 4],[lindex [lindex $Groups $i] 8],[lindex [lindex $Groups $i] 12]\],"
            if {[GiD_AccessValue get gendata Strategy_Type]=="Arc-Length"} {
                puts $varfile "            \"table\":           \[0,0,0\]"
            } else {
                set SearchInList [lsearch $TableList [lindex [lindex $Groups $i] 1]*]
                set AuxList [lindex $TableList [expr { $SearchInList+1 }]]
                puts $varfile "            \"table\":           \[[lindex $AuxList 0],[lindex $AuxList 1],[lindex $AuxList 2]\]"
            }
            puts $varfile "        \}"
            if {$iGroup < $NumGroups} {
                puts $varfile "    \},\{"
            } else {
                puts $varfile "    \}\],"
            }
        }
    } else {
        puts $varfile "    \"loads_process_list\":       \[\],"
    }

    ## output_configuration
    puts $varfile "    \"output_configuration\":     \{"
    puts $varfile "        \"result_file_configuration\": \{"
    puts $varfile "            \"gidpost_flags\":       \{"
    puts $varfile "                \"GiDPostMode\":           \"[GiD_AccessValue get gendata GiD_post_mode]\","
    puts $varfile "                \"WriteDeformedMeshFlag\": \"[GiD_AccessValue get gendata Write_deformed_mesh]\","
    puts $varfile "                \"WriteConditionsFlag\":   \"[GiD_AccessValue get gendata Write_conditions]\","
    puts $varfile "                \"MultiFileFlag\":         \"[GiD_AccessValue get gendata Multi_file_flag]\""
    puts $varfile "            \},"
    puts $varfile "            \"file_label\":          \"[GiD_AccessValue get gendata File_label]\","
    puts $varfile "            \"output_control_type\": \"[GiD_AccessValue get gendata Output_control_type]\","
    puts $varfile "            \"output_frequency\":    [GiD_AccessValue get gendata Output_frequency],"
    puts $varfile "            \"body_output\":         [GiD_AccessValue get gendata Body_output],"
    puts $varfile "            \"node_output\":         [GiD_AccessValue get gendata Node_output],"
    puts $varfile "            \"skin_output\":         [GiD_AccessValue get gendata Skin_output],"
    puts $varfile "            \"plane_output\":        \[\],"
    # nodal_results
    set PutStrings \[
    set iGroup 0
    if {[GiD_AccessValue get gendata Write_Solid_Displacement]==true} {
        incr iGroup
        append PutStrings \" DISPLACEMENT \" ,
    }
    if {[GiD_AccessValue get gendata Write_Fluid_Pressure]==true} {
        incr iGroup
        append PutStrings \" WATER_PRESSURE \" ,
    }
    if {[GiD_AccessValue get gendata Write_Reactions]==true} {
        incr iGroup
        append PutStrings \" REACTION \" , \" REACTION_WATER_PRESSURE \" ,
    }
    if {[GiD_AccessValue get gendata Write_Force]==true} {
        incr iGroup
        append PutStrings \" FORCE \" ,
    }
    if {[GiD_AccessValue get gendata Write_Face_Load]==true} {
        incr iGroup
        append PutStrings \" FACE_LOAD \" ,
    }
    if {[GiD_AccessValue get gendata Write_Normal_Load]==true} {
        incr iGroup
        append PutStrings \" NORMAL_CONTACT_STRESS \" ,
    }
    if {[GiD_AccessValue get gendata Write_Tangential_Load]==true} {
        incr iGroup
        append PutStrings \" TANGENTIAL_CONTACT_STRESS \" ,
    }
    if {[GiD_AccessValue get gendata Write_Normal_Fluid_Flux]==true} {
        incr iGroup
        append PutStrings \" NORMAL_FLUID_FLUX \" ,
    }
    if {[GiD_AccessValue get gendata Write_Body_Acceleration]==true} {
        incr iGroup
        append PutStrings \" VOLUME_ACCELERATION \" ,
    }
    if {$iGroup > 0} {
        set PutStrings [string trimright $PutStrings ,]
    }
    append PutStrings \]
    puts $varfile "            \"nodal_results\":       $PutStrings,"
    # gauss_point_results
    set PutStrings \[
    set iGroup 0
    if {[GiD_AccessValue get gendata Write_Strain]==true} {
        incr iGroup
        append PutStrings \" GREEN_LAGRANGE_STRAIN_TENSOR \" ,
    }
    if {[GiD_AccessValue get gendata Write_Effective_Stress]==true} {
        incr iGroup
        append PutStrings \" CAUCHY_STRESS_TENSOR \" ,
    }
    if {[GiD_AccessValue get gendata Write_Total_Stress]==true} {
        incr iGroup
        append PutStrings \" TOTAL_STRESS_TENSOR \" ,
    }
    if {[GiD_AccessValue get gendata Write_Von_Mises_Stress]==true} {
        incr iGroup
        append PutStrings \" VON_MISES_STRESS \" ,
    }
    if {[GiD_AccessValue get gendata Write_Fluid_Flux]==true} {
        incr iGroup
        append PutStrings \" FLUID_FLUX_VECTOR \" ,
    }
    if {[GiD_AccessValue get gendata Write_Permeability]==true} {
        incr iGroup
        append PutStrings \" PERMEABILITY_MATRIX \" ,
    }
    if {[GiD_AccessValue get gendata Write_Damage]==true} {
        incr iGroup
        append PutStrings \" DAMAGE_VARIABLE \" ,
    }
    if {[GiD_AccessValue get gendata Write_Local_Stress_Vector]==true} {
        incr iGroup
        append PutStrings \" LOCAL_STRESS_VECTOR \" ,
    }
    if {[GiD_AccessValue get gendata Write_Local_Relative_Displacement]==true} {
        incr iGroup
        append PutStrings \" LOCAL_RELATIVE_DISPLACEMENT_VECTOR \" ,
    }
    if {[GiD_AccessValue get gendata Write_Local_Fluid_Flux]==true} {
        incr iGroup
        append PutStrings \" LOCAL_FLUID_FLUX_VECTOR \" ,
    }
    if {[GiD_AccessValue get gendata Write_Local_Permeability]==true} {
        incr iGroup
        append PutStrings \" LOCAL_PERMEABILITY_MATRIX \" ,
    }
    if {$iGroup > 0} {
        set PutStrings [string trimright $PutStrings ,]
    }
    append PutStrings \]
    puts $varfile "            \"gauss_point_results\": $PutStrings"
    puts $varfile "        \},"
    puts $varfile "        \"point_data_configuration\":  \[\]"
    puts $varfile "    \},"
    
    ## restart_options
    puts $varfile "    \"restart_options\":          \{"
    puts $varfile "        \"SaveRestart\":      false,"
    puts $varfile "        \"RestartFrequency\": 0,"
    puts $varfile "        \"LoadRestart\":      false,"
    puts $varfile "        \"Restart_Step\":     0"
    puts $varfile "    \}"
    
    puts $varfile "\}"
    
    close $varfile
}
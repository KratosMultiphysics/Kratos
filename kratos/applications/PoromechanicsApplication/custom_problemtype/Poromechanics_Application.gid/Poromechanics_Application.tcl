## GiD events --------------------------------------------------------------------------------------------------------------------------------------------------

proc InitGIDProject { dir } {
    #::Poromechanics_Application::GetKratosPath
    
    GiDMenu::Create "Poromechanics Application" PRE
	GiDMenu::InsertOption "Poromechanics Application" [list "Dirichlet Boundary Conditions"] 0 PRE "GidOpenConditions \"Dirichlet_Boundary_Conditions\"" "" ""
	GiDMenu::InsertOption "Poromechanics Application" [list "Other Conditions"] 1 PRE "GidOpenConditions \"Other_Conditions\"" "" ""
	GiDMenu::InsertOption "Poromechanics Application" [list "Elements"] 2 PRE "GidOpenConditions \"Elements\"" "" ""
    GiDMenu::InsertOption "Poromechanics Application" [list "Materials"] 3 PRE "GidOpenMaterials" "" ""
    GiDMenu::InsertOption "Poromechanics Application" [list "Problem Parameters"] 4 PRE "GidOpenProblemData" "" ""
	GiDMenu::UpdateMenus
}

#-------------------------------------------------------------------------------

# Pass the path and the name of the problem to the Python script
proc BeforeRunCalculation { batfilename basename dir problemtypedir gidexe args } {
    set filename [file join $dir ${basename}-1.dat]
    set varfile [open $filename a]
    #puts $varfile "problem_name = '[file join $dir $basename]'"
    puts $varfile "problem_name = '${basename}'"
    puts $varfile "problem_path = '[file join $dir]'"
    #puts $varfile "gid_path = '${gidexe}'"
    #puts $varfile "kratos_path = '${::Poromechanics_Application::kratos_path}'"
    puts $varfile ""
    close $varfile
}

#-------------------------------------------------------------------------------

#proc AfterWriteCalcFileGIDProject { file error } {
#    WarnWin "Stopped after writing calculation file. Continue?"
#}


## Problemtype procedures --------------------------------------------------------------------------------------------------------------------------------------

namespace eval Poromechanics_Application {
    variable kratos_path ""
}

#-------------------------------------------------------------------------------

proc Poromechanics_Application::GetKratosPath { } {
    set knownpath 0
    set folder [file join $::env(HOME) Poromechanics_Application] 
    set setupfile [file join $folder "Poromechanics_Application.ini"]
    if { [file exists $setupfile] } {
        set setupdata [open $setupfile r]
        if { [gets $setupdata line] >= 0 } {
            if { [file isdirectory $line] } {
                set ::Poromechanics_Application::kratos_path $line
                set knownpath 1
        }
      }
    }

    if { $knownpath == 0 } {
        set title "Select the path to your Kratos folder"
        set ::Poromechanics_Application::kratos_path [tk_chooseDirectory -mustexist 1 -title $title ]

        set folder [file join $::env(HOME) Poromechanics_Application]
        if { ![file exists $folder] || ![file isdirectory $folder] } {
            file mkdir $folder
        }
        set filepath [file join $folder "Poromechanics_Application.ini"]
        set setupfile [open $filepath w]
        puts $setupfile $::Poromechanics_Application::kratos_path
        close $setupfile
    }
}

#-------------------------------------------------------------------------------

proc Poromechanics_Application::QuadrilateralInterface2D3Conectivities { ElemId } {
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
    
    if {$ly < $lx} {
        set Conectivities "$N1(Id) $N2(Id) $N3(Id) $N4(Id)"
    } else {
        set Conectivities "$N4(Id) $N1(Id) $N2(Id) $N3(Id)"
    }
    
    return $Conectivities
}

#-------------------------------------------------------------------------------

proc Poromechanics_Application::QuadrilateralInterface2D4Conectivities { ElemId } {
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
    
    if {$ly < $lx} {
        set Conectivities "$N1(Id) $N2(Id) $N3(Id) $N4(Id)"
    } else {
        set Conectivities "$N4(Id) $N1(Id) $N2(Id) $N3(Id)"
    }
    
    return $Conectivities
}

#-------------------------------------------------------------------------------

proc Poromechanics_Application::QuadrilateralInterface3D3Conectivities { ElemId } {
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
    
    if {$ly < $lx} {
        set Conectivities "$N1(Id) $N2(Id) $N3(Id) $N4(Id)"
    } else {
        set Conectivities "$N4(Id) $N1(Id) $N2(Id) $N3(Id)"
    }
    
    return $Conectivities
}

#-------------------------------------------------------------------------------

proc Poromechanics_Application::QuadrilateralInterface3D4Conectivities { ElemId } {
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
    
    if {$ly < $lx} {
        set Conectivities "$N1(Id) $N2(Id) $N3(Id) $N4(Id)"
    } else {
        set Conectivities "$N4(Id) $N1(Id) $N2(Id) $N3(Id)"
    }
    
    return $Conectivities
}

#-------------------------------------------------------------------------------

proc Poromechanics_Application::HexaedraInterface3D8Conectivities { ElemId } {
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
    
    if {$lz < $lx} {
        if {$lz < $ly} {
            # lz < lx && lz < ly
            set Conectivities "$N1(Id) $N2(Id) $N3(Id) $N4(Id) $N5(Id) $N6(Id) $N7(Id) $N8(Id)"
        } else {
            # ly < lz < lx
            set Conectivities "$N1(Id) $N4(Id) $N8(Id) $N5(Id) $N2(Id) $N3(Id) $N7(Id) $N6(Id)"
        }
    } elseif {$ly < $lx} {
        # ly < lx < lz
        set Conectivities "$N1(Id) $N4(Id) $N8(Id) $N5(Id) $N2(Id) $N3(Id) $N7(Id) $N6(Id)"
    } else {
        # lx < lz && lx < ly
        set Conectivities "$N1(Id) $N5(Id) $N6(Id) $N2(Id) $N4(Id) $N8(Id) $N7(Id) $N3(Id)"
    }
    
    return $Conectivities
}

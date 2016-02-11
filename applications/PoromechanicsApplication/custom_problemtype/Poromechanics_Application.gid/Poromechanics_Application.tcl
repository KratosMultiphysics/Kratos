## GiD events -----------------------------------------------------------------------------------------------

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

# Pass the path and the name of the problem to the Python script
proc BeforeRunCalculation { batfilename basename dir problemtypedir gidexe args } {
    set filename [file join $dir ${basename}-1.dat]
    set varfile [open $filename a]
    puts $varfile "problem_name = '[file join $dir $basename]'"
    puts $varfile "problem_path = '[file join $dir]'"
    #puts $varfile "gid_path = '${gidexe}'"
    #puts $varfile "kratos_path = '${::Poromechanics_Application::kratos_path}'"
    puts $varfile ""
    close $varfile
}


## Problemtype procedures -----------------------------------------------------------------------------------

namespace eval Poromechanics_Application {
    variable kratos_path ""
}

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

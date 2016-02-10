## GiD events -----------------------------------------------------------------------------------------------

proc InitGIDProject { dir } {
    #::Dam_Application::GetKratosPath
    
    GiDMenu::Create "Dam Application" PRE
	GiDMenu::InsertOption "Dam Application" [list "Dirichlet Boundary Conditions"] 0 PRE "GidOpenConditions \"Dirichlet_Boundary_Conditions\"" "" ""
	GiDMenu::InsertOption "Dam Application" [list "Load Conditions"] 1 PRE "GidOpenConditions \"Load_Conditions\"" "" ""
    GiDMenu::InsertOption "Dam Application" [list "Other Conditions"] 2 PRE "GidOpenConditions \"Other_Conditions\"" "" ""
	GiDMenu::InsertOption "Dam Application" [list "Elements"] 3 PRE "GidOpenConditions \"Elements\"" "" ""
    GiDMenu::InsertOption "Dam Application" [list "Materials"] 4 PRE "GidOpenMaterials" "" ""
    GiDMenu::InsertOption "Dam Application" [list "Problem Parameters"] 5 PRE "GidOpenProblemData" "" ""
	GiDMenu::UpdateMenus
}

# Pass the path and the name of the problem to the Python script
proc BeforeRunCalculation { batfilename basename dir problemtypedir gidexe args } {
    set filename [file join $dir ${basename}-1.dat]
    set varfile [open $filename a]
    puts $varfile "problem_name = '[file join $dir $basename]'"
    puts $varfile "problem_path = '[file join $dir]'"
    #puts $varfile "gid_path = '${gidexe}'"
    #puts $varfile "kratos_path = '${::Dam_Application::kratos_path}'"
    puts $varfile ""
    close $varfile
}


## Problemtype procedures -----------------------------------------------------------------------------------

namespace eval Dam_Application {
    variable kratos_path ""
}

proc Dam_Application::GetKratosPath { } {
    set knownpath 0
    set folder [file join $::env(HOME) Dam_Application] 
    set setupfile [file join $folder "Dam_Application.ini"]
    if { [file exists $setupfile] } {
        set setupdata [open $setupfile r]
        if { [gets $setupdata line] >= 0 } {
            if { [file isdirectory $line] } {
                set ::Dam_Application::kratos_path $line
                set knownpath 1
        }
      }
    }

    if { $knownpath == 0 } {
        set title "Select the path to your Kratos folder"
        set ::Dam_Application::kratos_path [tk_chooseDirectory -mustexist 1 -title $title ]

        set folder [file join $::env(HOME) Dam_Application]
        if { ![file exists $folder] || ![file isdirectory $folder] } {
            file mkdir $folder
        }
        set filepath [file join $folder "Dam_Application.ini"]
        set setupfile [open $filepath w]
        puts $setupfile $::Dam_Application::kratos_path
        close $setupfile
    }
}

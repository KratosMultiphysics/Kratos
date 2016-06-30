## GiD events --------------------------------------------------------------------------------------------------------------------------------------------------

proc InitGIDProject { dir } {
        
    GiDMenu::Create "Poromechanics Application" PRE
    GiDMenu::InsertOption "Poromechanics Application" [list "Parts"] 0 PRE "GidOpenConditions \"Parts\"" "" ""
	GiDMenu::InsertOption "Poromechanics Application" [list "Dirichlet Constraints"] 1 PRE "GidOpenConditions \"Dirichlet_Constraints\"" "" ""
	GiDMenu::InsertOption "Poromechanics Application" [list "Loads"] 2 PRE "GidOpenConditions \"Loads\"" "" ""
    GiDMenu::InsertOption "Poromechanics Application" [list "Project Parameters"] 3 PRE "GidOpenProblemData" "" ""
	GiDMenu::UpdateMenus
}

#-------------------------------------------------------------------------------

proc BeforeRunCalculation { batfilename basename dir problemtypedir gidexe args } {
    
    source [file join $problemtypedir WriteMdpa.tcl]
    
    set TableList [::Poromechanics_Application::WriteMDPA $basename $dir]
    
    source [file join $problemtypedir WriteProjectParameters.tcl]
    
    ::Poromechanics_Application::WriteProjectParameters $basename $dir $gidexe $TableList
}

## Problemtype procedures --------------------------------------------------------------------------------------------------------------------------------------

namespace eval Poromechanics_Application {
    
}

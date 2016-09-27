## GiD events --------------------------------------------------------------------------------------------------------------------------------------------------

proc InitGIDProject { dir } {
    
    # Initialize ProblemType Menu
    GiDMenu::Create "Poromechanics Application" PRE
    GiDMenu::InsertOption "Poromechanics Application" [list "Parts"] 0 PRE "GidOpenConditions \"Parts\"" "" ""
	GiDMenu::InsertOption "Poromechanics Application" [list "Dirichlet Constraints"] 1 PRE "GidOpenConditions \"Dirichlet_Constraints\"" "" ""
	GiDMenu::InsertOption "Poromechanics Application" [list "Loads"] 2 PRE "GidOpenConditions \"Loads\"" "" ""
    GiDMenu::InsertOption "Poromechanics Application" [list "Project Parameters"] 3 PRE "GidOpenProblemData" "" ""
	GiDMenu::UpdateMenus
    
    # Save ProblemTypePath
    set ::Poromechanics_Application::ProblemTypePath $dir
}

#-------------------------------------------------------------------------------

proc AfterReadGIDProject { filename } {
    
    # Save ProblemPath
    set projectpath $filename
    append projectpath .gid
    set ::Poromechanics_Application::ProblemPath $projectpath
    
    # Save ProblemName
    if {[regexp -all {\\} $filename] > 0} {
        # Windows
        regsub -all {\\} $filename { } filename
    } else {
        # Unix
        regsub -all {/} $filename { } filename
    }
    set filename [lreplace $filename 0 [expr { [llength $filename]-2 }]]
    set ::Poromechanics_Application::ProblemName $filename
}

#-------------------------------------------------------------------------------

proc BeforeRunCalculation { batfilename basename dir problemtypedir gidexe args } {
        
    # Write MDPA
    source [file join $problemtypedir WriteMdpa.tcl]
    set TableList [WriteMdpa $basename $dir]

    # Write ProjectParameters
    source [file join $problemtypedir WriteProjectParameters.tcl]
    WriteProjectParameters $basename $dir $TableList
        
    # Write Initial fractures data
    if {[GiD_AccessValue get gendata Fracture_Propagation]==true && [GiD_AccessValue get gendata Domain_Size]==2} {
        # Define GiDPath
        if {[regexp -all {\\} $gidexe] > 0} {
            # Windows
            regsub -all {\\} $gidexe {/} gidexe
        }
        set gidexe [string trimright $gidexe gid.exe]
        
        source [file join $problemtypedir FracturePropagation2D.tcl]
        WriteInitialFracturesData $dir $gidexe
    }
}


## Problemtype procedures --------------------------------------------------------------------------------------------------------------------------------------

namespace eval Poromechanics_Application {
    variable ProblemName ""
    variable ProblemPath ""
    variable ProblemTypePath ""
}

#-------------------------------------------------------------------------------

proc Poromechanics_Application::PropagateFractures2D { } {
    
    # Get Propagation Data
    source [file join $::Poromechanics_Application::ProblemPath PropagationData.tcl]
    set PropagationData [GetPropagationData]
    
    # Generate New Geometry and Write New FracturesData
    source [file join $::Poromechanics_Application::ProblemTypePath FracturePropagation2D.tcl]
    GenerateNewFractures $::Poromechanics_Application::ProblemPath $PropagationData
    GiD_Process Mescape Files Save
    
    # Write new MDPA
    source [file join $::Poromechanics_Application::ProblemTypePath WriteMdpa.tcl]
    set TableList [WriteMdpa $::Poromechanics_Application::ProblemName $::Poromechanics_Application::ProblemPath]

    # Write new ProjectParameters
    source [file join $::Poromechanics_Application::ProblemTypePath WriteProjectParameters.tcl]
    WriteProjectParameters $::Poromechanics_Application::ProblemName $::Poromechanics_Application::ProblemPath $TableList
    
    GiD_Process Mescape Files Save
    GiD_Process escape escape escape escape escape Quit
}
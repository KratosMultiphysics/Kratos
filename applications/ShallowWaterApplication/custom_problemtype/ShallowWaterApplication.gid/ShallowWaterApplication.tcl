## GiD events ------------------------------------------------------------------

proc InitGIDProject { dir } {
    
    # Initialize ProblemType Menu
    if { [GidUtils::IsTkDisabled] eq 0} {  
	GiDMenu::Create "Shallow Water Application" PRE
	GiDMenu::InsertOption "Shallow Water Application" [list "Parts"] 0 PRE "GidOpenConditions \"Parts\"" "" ""
	GiDMenu::InsertOption "Shallow Water Application" [list "Initial conditions"] 1 PRE "GidOpenConditions \"Initial_conditions\"" "" ""
	GiDMenu::InsertOption "Shallow Water Application" [list "Boundary conditions"] 2 PRE "GidOpenConditions \"Boundary_conditions\"" "" ""
	GiDMenu::InsertOption "Shallow Water Application" [list "Source terms"] 3 PRE "GidOpenConditions \"Source_terms\"" "" ""
	GiDMenu::InsertOption "Shallow Water Application" [list "Project Parameters"] 4 PRE "GidOpenProblemData" "" ""
	GiDMenu::UpdateMenus
    }

    # Save ProblemTypePath
    set ::ShallowWaterApplication::ProblemTypePath $dir
}

#-------------------------------------------------------------------------------

proc AfterReadGIDProject { filename } {
    
    # Save ProblemPath
    set projectpath $filename
    append projectpath .gid
    set ::ShallowWaterApplication::ProblemPath $projectpath
    
    # Save ProblemName
    # if {$::tcl_platform(platform) eq "windows"} {}
    if {[regexp -all {\\} $filename] > 0} {
	# Windows
	regsub -all {\\} $filename { } filename
    } else {
	# Unix
	regsub -all {/} $filename { } filename
    }
    set filename [lreplace $filename 0 [expr { [llength $filename]-2 }]]
    set ::ShallowWaterApplication::ProblemName $filename
}

#-------------------------------------------------------------------------------

proc BeforeRunCalculation { batfilename basename dir problemtypedir gidexe args } {
    
    # Set Parallel Configuration
    set paralleltype [GiD_AccessValue get gendata Parallel_Configuration]

    # Write MDPA
    source [file join $problemtypedir Mdpa.tcl]
    WriteMdpa $basename $dir $problemtypedir

    # Write ProjectParameters
    source [file join $problemtypedir ProjectParameters.tcl]
    WriteProjectParameters $basename $dir $problemtypedir

    # Copy python script in the problemdir
    file copy -force [file join $problemtypedir shallow_water_main.py] [file join $dir MainKratos.py]

    # Run the problem
    set run 1
    catch {
	if {$paralleltype eq "MPI"} {set run 0}
    }
    if {$run} {
	return ""
    } {
	return [list "-cancel-" [= "You have selected MPI parallelism system.\n\
		                    Input files have been written.\n\
		                    Run the case with: mpirun -np \[npartitions\] python3 MainKratos.py" ]]
    }
}


## Problemtype procedures ------------------------------------------------------

namespace eval ShallowWaterApplication {
    variable ProblemName ""
    variable ProblemPath ""
    variable ProblemTypePath ""
}

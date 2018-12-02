
## GiD events --------------------------------------------------------------------------------------------------------------------------------------------------

proc InitGIDProject { dir } {
    
# Initialize ProblemType Menu
    if { [GidUtils::IsTkDisabled] eq 0} {  
        GiDMenu::Create "FemDemKratos3D Application" PRE
        GiDMenu::InsertOption "FemDemKratos3D Application" [list "Parts"] 0 PRE "GidOpenConditions \"Parts\"" "" ""
        GiDMenu::InsertOption "FemDemKratos3D Application" [list "Dirichlet Constraints"] 1 PRE "GidOpenConditions \"Dirichlet_Constraints\"" "" ""
        GiDMenu::InsertOption "FemDemKratos3D Application" [list "Loads"] 2 PRE "GidOpenConditions \"Loads\"" "" ""
        GiDMenu::InsertOption "FemDemKratos3D Application" [list "Project Parameters"] 3 PRE "GidOpenProblemData" "" ""
        GiDMenu::InsertOption "FemDemKratos3D Application" [list "Plots"] 4 PRE "GidOpenConditions \"Plots\"" "" ""
        GiDMenu::UpdateMenus
    }

    # Save ProblemTypePath
    set ::FemDemKratos::ProblemTypePath $dir
}

#-------------------------------------------------------------------------------


proc AfterReadGIDProject { filename } {
    
    # Save ProblemPath
    set projectpath $filename
    append projectpath .gid
    set ::FemDemKratos::ProblemPath $projectpath

    # Save ProblemName
    if {[regexp -all {\\} $filename] > 0} {
        # Windows
        regsub -all {\\} $filename { } filename
    } else {
        # Unix
        regsub -all {/} $filename { } filename
    }
    set filename [lreplace $filename 0 [expr { [llength $filename]-2 }]]
    set ::FemDemKratos::ProblemName $filename
}

#-------------------------------------------------------------------------------

proc BeforeRunCalculation { batfilename basename dir problemtypedir gidexe args } {  

#---------------------------------------------------------------
    # Write MDPA
    source [file join $problemtypedir Mdpa.tcl]
    set TableDict [WriteMdpa $basename $dir $problemtypedir]

    # Write ProjectParameters
    source [file join $problemtypedir ProjectParameters.tcl]
    WriteProjectParameters $basename $dir $problemtypedir $TableDict
#---------------------------------------------------------------

    # Write materials.py
    source [file join $problemtypedir Materials.tcl]
    WriteMaterials $basename $dir $problemtypedir $TableDict
#---------------------------------------------------------------

    # If MMG remeshing is activated, we copy the .json parameters
    if {[GiD_AccessValue get gendata Activate_MMG_Remeshing] eq "true"} {
        file copy -force [file join $problemtypedir MMGParameters.json] [file join $dir MMGParameters.json]
    }

#---------------------------------------------------------------
    
    # For Coupled calculations with DEM elements
    if {[GiD_AccessValue get gendata Coupled_Calculation] eq "true"} {

        # Writes the mdpa of the discrete elements (only properties)
        source [file join $problemtypedir MdpaDEM.tcl]
        WriteMdpaDEM $basename $dir $problemtypedir

        source [file join $problemtypedir ProjectParametersDEM.tcl]
        WriteProjectParametersDEM $basename $dir $problemtypedir

        file copy -force [file join $problemtypedir KratosFemDemCoupled3DApplication.py] [file join $dir KratosFemDemCoupled3DApplication.py]
    } else {
        file copy -force [file join $problemtypedir DEM_explicit_solver_var.py] [file join $dir DEM_explicit_solver_var.py]
        file copy -force [file join $problemtypedir KratosFemDemApplication.py] [file join $dir KratosFemDemApplication.py]
    }
#---------------------------------------------------------------
    
}


## Problemtype procedures --------------------------------------------------------------------------------------------------------------------------------------

namespace eval FemDemKratos {
    variable ProblemName ""
    variable ProblemPath ""
    variable ProblemTypePath ""
}
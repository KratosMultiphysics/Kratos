## GiD events --------------------------------------------------------------------------------------------------------------------------------------------------

proc InitGIDProject { dir } {

    # Initialize ProblemType Menu
    if { [GidUtils::IsTkDisabled] eq 0} {
        GiDMenu::Create "Poromechanics Application" PRE
        GiDMenu::InsertOption "Poromechanics Application" [list "Parts"] 0 PRE "GidOpenConditions \"Parts\"" "" ""
        GiDMenu::InsertOption "Poromechanics Application" [list "Dirichlet Constraints"] 1 PRE "GidOpenConditions \"Dirichlet_Constraints\"" "" ""
        GiDMenu::InsertOption "Poromechanics Application" [list "Loads"] 2 PRE "GidOpenConditions \"Loads\"" "" ""
        GiDMenu::InsertOption "Poromechanics Application" [list "Project Parameters"] 3 PRE "GidOpenProblemData" "" ""
        GiDMenu::UpdateMenus
    }

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
    # if {$::tcl_platform(platform) eq "windows"} {}
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

    # Set Parallel Configuration
    set paralleltype [GiD_AccessValue get gendata Parallel_Configuration]

    # Write Initial fractures data
    if {([GiD_AccessValue get gendata Fracture_Propagation] eq true) && ($paralleltype ne "MPI")} {
        # Define GiDPath
        if {[regexp -all {\\} $gidexe] > 0} {
            # Windows
            regsub -all {\\} $gidexe {/} gidexe
        }
        set gidexe [string trimright $gidexe gid.exe]

        if {[GiD_AccessValue get gendata Domain_Size] eq 2} {

            source [file join $problemtypedir FracturePropagation2D.tcl]
            WriteInitialFracturesData $dir $problemtypedir $gidexe

        } else {

            source [file join $problemtypedir FracturePropagation3D.tcl]
            WriteInitialFracturesData $dir $problemtypedir $gidexe
        }
    }

    # Write MDPA
    source [file join $problemtypedir Mdpa.tcl]
    set TableDict [WriteMdpa $basename $dir $problemtypedir]

    # Write ProjectParameters
    source [file join $problemtypedir ProjectParameters.tcl]
    WriteProjectParameters $basename $dir $problemtypedir $TableDict

    # Copy python script in the problemdir
    if {[GiD_AccessValue get gendata Fracture_Propagation] eq true} {
        file copy -force [file join $problemtypedir poromechanics_fracture_main.py] [file join $dir MainKratos.py]
    } else {
        file copy -force [file join $problemtypedir KratosPoromechanics.py] [file join $dir MainKratos.py]
    }

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

    ### Measure time
    #set start_time_1 [clock clicks]
    #set end_time_1 [expr { [clock clicks]-$start_time_1 }]
    #WarnWin "Time for GenerateNewFractures: $end_time_1 clicks"
    ###
}


## Problemtype procedures --------------------------------------------------------------------------------------------------------------------------------------

namespace eval Poromechanics_Application {
    variable ProblemName ""
    variable ProblemPath ""
    variable ProblemTypePath ""
}

#-------------------------------------------------------------------------------

proc Poromechanics_Application::PropagateFractures2D { } {

    # Source Propagation Data and file
    set PropagationData [source [file join $::Poromechanics_Application::ProblemPath PropagationData.tcl]]
    source [file join $::Poromechanics_Application::ProblemTypePath FracturePropagation2D.tcl]

    # Generate New Geometry and Write New FracturesData
    GenerateNewFractures $::Poromechanics_Application::ProblemPath $::Poromechanics_Application::ProblemTypePath $PropagationData

    # Write new MDPA
    source [file join $::Poromechanics_Application::ProblemTypePath Mdpa.tcl]
    set TableDict [WriteMdpa $::Poromechanics_Application::ProblemName $::Poromechanics_Application::ProblemPath $::Poromechanics_Application::ProblemTypePath]

    # Write new ProjectParameters
    source [file join $::Poromechanics_Application::ProblemTypePath ProjectParameters.tcl]
    WriteProjectParameters $::Poromechanics_Application::ProblemName $::Poromechanics_Application::ProblemPath $::Poromechanics_Application::ProblemTypePath $TableDict

    # Quit GiD
    GiD_Process Mescape Files Save
    GiD_Process escape escape escape escape escape Quit
}

#-------------------------------------------------------------------------------

proc Poromechanics_Application::PropagateFractures3D { } {

    # Source Propagation Data and file
    set PropagationData [source [file join $::Poromechanics_Application::ProblemPath PropagationData.tcl]]
    source [file join $::Poromechanics_Application::ProblemTypePath FracturePropagation3D.tcl]

    # Generate New Geometry and Write New FracturesData
    GenerateNewFractures $::Poromechanics_Application::ProblemPath $::Poromechanics_Application::ProblemTypePath $PropagationData

    # Write new MDPA
    source [file join $::Poromechanics_Application::ProblemTypePath Mdpa.tcl]
    set TableDict [WriteMdpa $::Poromechanics_Application::ProblemName $::Poromechanics_Application::ProblemPath $::Poromechanics_Application::ProblemTypePath]

    # Write new ProjectParameters
    source [file join $::Poromechanics_Application::ProblemTypePath ProjectParameters.tcl]
    WriteProjectParameters $::Poromechanics_Application::ProblemName $::Poromechanics_Application::ProblemPath $::Poromechanics_Application::ProblemTypePath $TableDict

    # Quit GiD
    GiD_Process Mescape Files Save
    GiD_Process escape escape escape escape escape Quit
}

#--------------------------------------------------------------------------------------------------------------------------------------------------------------

proc Poromechanics_Application::CreateContactEntity { } {

    set LayerName "Layer0"

    ## Contact Surface:

    # set Line1 [list 111 1]
    # set Line2 [list 110 0]

    # # Orientation of Line1:
    # #     0: SAME1ST. The normal to the line points towards the oppsite direction of the contact surface.
    # #     1: DIFF1ST. The normal to the line points towards the contact surface.
    # # Orientation of Line2:
    # #     0: SAME1ST. The normal to the line points towards the contact surface.
    # #     1: DIFF1ST. The normal to the line points towards the oppsite direction of the contact surface.

    # GiD_Geometry create surface append contactsurface \
    #     $LayerName 2 $Line1 $Line2


    ## Contact Volume:

    # set Surf1 [list 5 1]
    # set Surf2 [list 186 0]

    # # Orientation of Surf1:
    # #     0: SAME1ST. The normal to the surface points towards the contact volume.
    # #     1: DIFF1ST. The normal to the surface points towards the oppsite direction of the contact volume.
    # # Orientation of Surf2:
    # #     0: DIFF1ST. The normal to the surface points towards the oppsite direction of the contact volume.
    # #     1: SAME1ST. The normal to the surface points towards the contact volume.

    # Transformation matrix for a contact volume between Surf1 and Surf2. R means rotation, and T translation.

    # set TransformMatrix [list Rxx Rxy Rxz Tx \
    #                           Ryx Ryy Ryz Ty \
    #                           Rzx Rzy Rzz Tz \
    #                           0.0 0.0 0.0 1.0]

    # Transformation matrix for a contact volume between two equal surfaces, one over the other (zero-thickness interface elements)

    # set TransformMatrix [list 1.0 0.0 0.0 0.0 \
    #                           0.0 1.0 0.0 0.0 \
    #                           0.0 0.0 1.0 0.0 \
    #                           0.0 0.0 0.0 1.0]

    # Transformation matrix for a contact volume generated by a translation from Surf1 to Surf2

    # set TransformMatrix [list 1.0 0.0 0.0 Tx \
    #                           0.0 1.0 0.0 Ty \
    #                           0.0 0.0 1.0 Tz \
    #                           0.0 0.0 0.0 1.0]

    # Transformation matrix for a contact volume generated by a normal translation from Surf1 to Surf2

    # # Vertex of Surf1
    # set Vertex [GiD_Geometry get point 90]
    # # Point1 of Surf1
    # set Point1 [GiD_Geometry get point 91]
    # # Point2 of Surf1
    # set Point2 [GiD_Geometry get point 93]
    # # Distance of translation
    # set Distance 0.01
    # set TransformMatrix [NormalTranslationMatrix $Vertex $Point1 $Point2 $Distance]

    # Transformation matrix for a contact volume generated by a rotation around a given axis from Surf1 to Surf2:

    # # Point at the initial position of the rotation
    # set InitPoint [GiD_Geometry get point 17]
    # # Point at the final position of the rotation
    # set FinalPoint [GiD_Geometry get point 18]
    # # Vertex of the angle of rotation
    # set Vertex [GiD_Geometry get point 19]
    # # Initial point of the rotation Axis
    # set InitAxis [GiD_Geometry get point 20]
    # # Final point of the rotation Axis
    # set FinalAxis [GiD_Geometry get point 19]
    # set TransformMatrix [RotationMatrix $InitPoint $FinalPoint $Vertex $InitAxis $FinalAxis]

    # GiD_Geometry create volume append $LayerName 2 \
    #    $Surf1 $Surf2 contactvolume $TransformMatrix
}

#-------------------------------------------------------------------------------

proc NormalTranslationMatrix {Vertex Point1 Point2 Distance} {
    # Vector in local x direction
    set Vx(0) [expr {[lindex $Point1 1]-[lindex $Vertex 1]}]
    set Vx(1) [expr {[lindex $Point1 2]-[lindex $Vertex 2]}]
    set Vx(2) [expr {[lindex $Point1 3]-[lindex $Vertex 3]}]

    # Vector in local y direction
    set Vy(0) [expr {[lindex $Point2 1]-[lindex $Vertex 1]}]
    set Vy(1) [expr {[lindex $Point2 2]-[lindex $Vertex 2]}]
    set Vy(2) [expr {[lindex $Point2 3]-[lindex $Vertex 3]}]

    # Vector in local z direction (Cross product between Vx and Vy)
    set Vz(0) [expr {$Vx(1)*$Vy(2)-$Vx(2)*$Vy(1)}]
    set Vz(1) [expr {$Vx(2)*$Vy(0)-$Vx(0)*$Vy(2)}]
    set Vz(2) [expr {$Vx(0)*$Vy(1)-$Vx(1)*$Vy(0)}]
    set InvNorm [expr {1.0/sqrt($Vz(0)*$Vz(0)+$Vz(1)*$Vz(1)+$Vz(2)*$Vz(2))}]
    set Vz(0) [expr {$Vz(0)*$InvNorm}]
    set Vz(1) [expr {$Vz(1)*$InvNorm}]
    set Vz(2) [expr {$Vz(2)*$InvNorm}]

    set Tx [expr {$Distance*Vz(0)}]
    set Ty [expr {$Distance*Vz(1)}]
    set Tz [expr {$Distance*Vz(2)}]

    return [list 1.0 0.0 0.0 $Tx \
                 0.0 1.0 0.0 $Ty \
                 0.0 0.0 1.0 $Tz \
                 0.0 0.0 0.0 1.0]

}

#-------------------------------------------------------------------------------

proc RotationMatrix {InitPoint FinalPoint Vertex InitAxis FinalAxis} {
    # Unitary vector at the initial position of the Rotation
    set Ri(0) [expr {[lindex $InitPoint 1]-[lindex $Vertex 1]}]
    set Ri(1) [expr {[lindex $InitPoint 2]-[lindex $Vertex 2]}]
    set Ri(2) [expr {[lindex $InitPoint 3]-[lindex $Vertex 3]}]
    set InvNorm [expr {1.0/sqrt($Ri(0)*$Ri(0)+$Ri(1)*$Ri(1)+$Ri(2)*$Ri(2))}]
    set Ri(0) [expr {$Ri(0)*$InvNorm}]
    set Ri(1) [expr {$Ri(1)*$InvNorm}]
    set Ri(2) [expr {$Ri(2)*$InvNorm}]
    # Unitary vector at the final position of the Rotation
    set Rf(0) [expr {[lindex $FinalPoint 1]-[lindex $Vertex 1]}]
    set Rf(1) [expr {[lindex $FinalPoint 2]-[lindex $Vertex 2]}]
    set Rf(2) [expr {[lindex $FinalPoint 3]-[lindex $Vertex 3]}]
    set InvNorm [expr {1.0/sqrt($Rf(0)*$Rf(0)+$Rf(1)*$Rf(1)+$Rf(2)*$Rf(2))}]
    set Rf(0) [expr {$Rf(0)*$InvNorm}]
    set Rf(1) [expr {$Rf(1)*$InvNorm}]
    set Rf(2) [expr {$Rf(2)*$InvNorm}]
    # Unitary rotation Axis
    set A(0) [expr {[lindex $FinalAxis 1]-[lindex $InitAxis 1]}]
    set A(1) [expr {[lindex $FinalAxis 2]-[lindex $InitAxis 2]}]
    set A(2) [expr {[lindex $FinalAxis 3]-[lindex $InitAxis 3]}]
    set InvNorm [expr {1.0/sqrt($A(0)*$A(0)+$A(1)*$A(1)+$A(2)*$A(2))}]
    set A(0) [expr {$A(0)*$InvNorm}]
    set A(1) [expr {$A(1)*$InvNorm}]
    set A(2) [expr {$A(2)*$InvNorm}]

    # Cosine of the angle of rotation
    set CosAngle [expr {$Ri(0)*$Rf(0)+$Ri(1)*$Rf(1)+$Ri(2)*$Rf(2)}]
    # Cross product between vectors Ri and Rf
    set n(0) [expr {$Ri(1)*$Rf(2)-$Ri(2)*$Rf(1)}]
    set n(1) [expr {$Ri(2)*$Rf(0)-$Ri(0)*$Rf(2)}]
    set n(2) [expr {$Ri(0)*$Rf(1)-$Ri(1)*$Rf(0)}]
    # Sine of the angle of rotation (positive angle between 0º and 90º)
    set SinAngle [expr {sqrt($n(0)*$n(0)+$n(1)*$n(1)+$n(2)*$n(2))}]

    # Transformation Matrix
    set Rxx [expr {$CosAngle+$A(0)*$A(0)*(1.0-$CosAngle)}]
    set Rxy [expr {$A(0)*$A(1)*(1.0-$CosAngle)-$A(2)*$SinAngle}]
    set Rxz [expr {$A(0)*$A(2)*(1.0-$CosAngle)+$A(1)*$SinAngle}]
    set Ryx [expr {$A(0)*$A(1)*(1.0-$CosAngle)+$A(2)*$SinAngle}]
    set Ryy [expr {$CosAngle+$A(1)*$A(1)*(1.0-$CosAngle)}]
    set Ryz [expr {$A(1)*$A(2)*(1.0-$CosAngle)-$A(0)*$SinAngle}]
    set Rzx [expr {$A(0)*$A(2)*(1.0-$CosAngle)-$A(1)*$SinAngle}]
    set Rzy [expr {$A(1)*$A(2)*(1.0-$CosAngle)+$A(0)*$SinAngle}]
    set Rzz [expr {$CosAngle+$A(2)*$A(2)*(1.0-$CosAngle)}]
    set Tx [expr {[lindex $InitAxis 1]-($Rxx*[lindex $InitAxis 1]+$Rxy*[lindex $InitAxis 2]+$Rxz*[lindex $InitAxis 3])}]
    set Ty [expr {[lindex $InitAxis 2]-($Ryx*[lindex $InitAxis 1]+$Ryy*[lindex $InitAxis 2]+$Ryz*[lindex $InitAxis 3])}]
    set Tz [expr {[lindex $InitAxis 3]-($Rzx*[lindex $InitAxis 1]+$Rzy*[lindex $InitAxis 2]+$Rzz*[lindex $InitAxis 3])}]

    return [list $Rxx $Rxy $Rxz $Tx \
                 $Ryx $Ryy $Ryz $Ty \
                 $Rzx $Rzy $Rzz $Tz \
                 0.0 0.0 0.0 1.0]
}

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
        
    # Write Initial fractures data
    if {[GiD_AccessValue get gendata Fracture_Propagation] eq true} {
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
    
    ## Contact Surface:
    
    #~ set Line1 [list 111 1]
    #~ set Line2 [list 110 0]
    
    # Orientation of Line1:
        # 0: SAME1ST. The normal to the line points towards the oppsite direction of the contact surface.
        # 1: DIFF1ST. The normal to the line points towards the contact surface.
    # Orientation of Line2:
        # 0: SAME1ST. The normal to the line points towards the contact surface.
        # 1: DIFF1ST. The normal to the line points towards the oppsite direction of the contact surface.
    
    #~ GiD_Geometry create surface append contactsurface \
        #~ Layer0 2 $Line1 $Line2
    
    
    ## Contact Volume:

    #~ set Surf1 [list 8 1]
    #~ set Surf2 [list 9 0]

    # Orientation of Surf1:
        # 0: SAME1ST. The normal to the surface points towards the contact volume.
        # 1: DIFF1ST. The normal to the surface points towards the oppsite direction of the contact volume.
    # Orientation of Surf2:
        # 0: DIFF1ST. The normal to the surface points towards the oppsite direction of the contact volume.
        # 1: SAME1ST. The normal to the surface points towards the contact volume.
    
    # Transformation matrix for a contact volume between Surf1 and Surf2. R means rotation, and T translation.
    
    #~ set TransformMatrix [list Rxx Rxy Rxz Tx \
                              #~ Ryx Ryy Ryz Ty \
                              #~ Rzx Rzy Rzz Tz \
                              #~ 0.0 0.0 0.0 1.0]
    
    # Transformation matrix for a contact volume between two equal surfaces, one over the other (zero-thickness interface elements)
    
    #~ set TransformMatrix [list 1.0 0.0 0.0 0.0 \
                              #~ 0.0 1.0 0.0 0.0 \
                              #~ 0.0 0.0 1.0 0.0 \
                              #~ 0.0 0.0 0.0 1.0]
    
    # Transformation matrix for a contact volume generated by a translation from Surf1 to Surf2
    
    #~ set TransformMatrix [list 1.0 0.0 0.0 Tx \
                              #~ 0.0 1.0 0.0 Ty \
                              #~ 0.0 0.0 1.0 Tz \
                              #~ 0.0 0.0 0.0 1.0]

    # Transformation matrix for a contact volume generated by a normal translation from Surf1 to Surf2
    
    #~ # Vertex of Surf1
    #~ set Vertex [GiD_Geometry get point 90]
    #~ # Point1 of Surf1
    #~ set Point1 [GiD_Geometry get point 91]
    #~ # Point2 of Surf1
    #~ set Point2 [GiD_Geometry get point 93]
    
    #~ # Compute unitary normal vector
    #~ set UnitNormal [ComputeUnitNormal $Vertex $Point1 $Point2]
    #~ set Distance 0.1
    #~ set Tx [expr {$Distance*[lindex $UnitNormal 0]}]
    #~ set Ty [expr {$Distance*[lindex $UnitNormal 1]}]
    #~ set Tz [expr {$Distance*[lindex $UnitNormal 2]}]
    
    #~ set TransformMatrix [list 1.0 0.0 0.0 $Tx \
                              #~ 0.0 1.0 0.0 $Ty \
                              #~ 0.0 0.0 1.0 $Tz \
                              #~ 0.0 0.0 0.0 1.0]

    # Transformation matrix for a contact volume generated by a rotation around a given axis from Surf1 to Surf2:
    
    #~ # Point at the initial position of the rotation
    #~ set InitPoint [GiD_Geometry get point 17]
    #~ set Pi(0) [lindex $InitPoint 1]
    #~ set Pi(1) [lindex $InitPoint 2]
    #~ set Pi(2) [lindex $InitPoint 3]
    #~ # Point at the final position of the rotation
    #~ set FinalPoint [GiD_Geometry get point 18]
    #~ set Pf(0) [lindex $FinalPoint 1]
    #~ set Pf(1) [lindex $FinalPoint 2]
    #~ set Pf(2) [lindex $FinalPoint 3]
    #~ # Vertex of the angle of rotation
    #~ set Vertex [GiD_Geometry get point 19]
    #~ set V(0) [lindex $Vertex 1]
    #~ set V(1) [lindex $Vertex 2]
    #~ set V(2) [lindex $Vertex 3]
    #~ # Initial point of the rotation Axis
    #~ set InitAxis [GiD_Geometry get point 20]
    #~ set Ai(0) [lindex $InitAxis 1]
    #~ set Ai(1) [lindex $InitAxis 2]
    #~ set Ai(2) [lindex $InitAxis 3]
    #~ # Final point of the rotation Axis
    #~ set FinalAxis [GiD_Geometry get point 19]
    #~ set Af(0) [lindex $FinalAxis 1]
    #~ set Af(1) [lindex $FinalAxis 2]
    #~ set Af(2) [lindex $FinalAxis 3]

    #~ # Unitary vector at the initial position of the Rotation
    #~ set Ri(0) [expr {$Pi(0)-$V(0)}]
    #~ set Ri(1) [expr {$Pi(1)-$V(1)}]
    #~ set Ri(2) [expr {$Pi(2)-$V(2)}]
    #~ set InvNorm [expr {1.0/sqrt($Ri(0)*$Ri(0)+$Ri(1)*$Ri(1)+$Ri(2)*$Ri(2))}]
    #~ set Ri(0) [expr {$Ri(0)*$InvNorm}]
    #~ set Ri(1) [expr {$Ri(1)*$InvNorm}]
    #~ set Ri(2) [expr {$Ri(2)*$InvNorm}]
    #~ # Unitary vector at the final position of the Rotation
    #~ set Rf(0) [expr {$Pf(0)-$V(0)}]
    #~ set Rf(1) [expr {$Pf(1)-$V(1)}]
    #~ set Rf(2) [expr {$Pf(2)-$V(2)}]
    #~ set InvNorm [expr {1.0/sqrt($Rf(0)*$Rf(0)+$Rf(1)*$Rf(1)+$Rf(2)*$Rf(2))}]
    #~ set Rf(0) [expr {$Rf(0)*$InvNorm}]
    #~ set Rf(1) [expr {$Rf(1)*$InvNorm}]
    #~ set Rf(2) [expr {$Rf(2)*$InvNorm}]
    #~ # Unitary rotation Axis
    #~ set A(0) [expr {$Af(0)-$Ai(0)}]
    #~ set A(1) [expr {$Af(1)-$Ai(1)}]
    #~ set A(2) [expr {$Af(2)-$Ai(2)}]
    #~ set InvNorm [expr {1.0/sqrt($A(0)*$A(0)+$A(1)*$A(1)+$A(2)*$A(2))}]
    #~ set A(0) [expr {$A(0)*$InvNorm}]
    #~ set A(1) [expr {$A(1)*$InvNorm}]
    #~ set A(2) [expr {$A(2)*$InvNorm}]
    
    #~ # Cosine of the angle of rotation
    #~ set CosAngle [expr {$Ri(0)*$Rf(0)+$Ri(1)*$Rf(1)+$Ri(2)*$Rf(2)}]
    #~ # Cross product between vectors Ri and Rf
    #~ set n(0) [expr {$Ri(1)*$Rf(2)-$Ri(2)*$Rf(1)}]
    #~ set n(1) [expr {$Ri(2)*$Rf(0)-$Ri(0)*$Rf(2)}]
    #~ set n(2) [expr {$Ri(0)*$Rf(1)-$Ri(1)*$Rf(0)}]
    #~ # Sine of the angle of rotation (positive angle between 0º and 90º)
    #~ set SinAngle [expr {sqrt($n(0)*$n(0)+$n(1)*$n(1)+$n(2)*$n(2))}]
    
    #~ # Transformation Matrix
    #~ set T(00) [expr {$CosAngle+$A(0)*$A(0)*(1.0-$CosAngle)}]
    #~ set T(01) [expr {$A(0)*$A(1)*(1.0-$CosAngle)-$A(2)*$SinAngle}]
    #~ set T(02) [expr {$A(0)*$A(2)*(1.0-$CosAngle)+$A(1)*$SinAngle}]
    #~ set T(10) [expr {$A(0)*$A(1)*(1.0-$CosAngle)+$A(2)*$SinAngle}]
    #~ set T(11) [expr {$CosAngle+$A(1)*$A(1)*(1.0-$CosAngle)}]
    #~ set T(12) [expr {$A(1)*$A(2)*(1.0-$CosAngle)-$A(0)*$SinAngle}]
    #~ set T(20) [expr {$A(0)*$A(2)*(1.0-$CosAngle)-$A(1)*$SinAngle}]
    #~ set T(21) [expr {$A(1)*$A(2)*(1.0-$CosAngle)+$A(0)*$SinAngle}]
    #~ set T(22) [expr {$CosAngle+$A(2)*$A(2)*(1.0-$CosAngle)}]
    #~ set T(03) [expr {$Ai(0)-($T(00)*$Ai(0)+$T(01)*$Ai(1)+$T(02)*$Ai(2))}]
    #~ set T(13) [expr {$Ai(1)-($T(10)*$Ai(0)+$T(11)*$Ai(1)+$T(12)*$Ai(2))}]
    #~ set T(23) [expr {$Ai(2)-($T(20)*$Ai(0)+$T(21)*$Ai(1)+$T(22)*$Ai(2))}]
    
    #~ set TransformMatrix [list $T(00) $T(01) $T(02) $T(03) \
                              #~ $T(10) $T(11) $T(12) $T(13) \
                              #~ $T(20) $T(21) $T(22) $T(23) \
                              #~ 0.0 0.0 0.0 1.0]
        
    #~ GiD_Geometry create volume append Layer0 2 \
       #~ $Surf1 $Surf2 contactvolume $TransformMatrix
}

#-------------------------------------------------------------------------------

proc ComputeUnitNormal {Vertex Point1 Point2} {
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
    
    return "$Vz(0) $Vz(1) $Vz(2)"
}
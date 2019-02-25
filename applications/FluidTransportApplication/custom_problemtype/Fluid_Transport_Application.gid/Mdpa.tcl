proc WriteMdpa { basename dir problemtypedir } {

    ## Source auxiliar procedures
    source [file join $problemtypedir MdpaAuxProcs.tcl]

    ## Start MDPA file
    set filename [file join $dir ${basename}.mdpa]
    set FileVar [open $filename w]

    ## ModelPart Data
    #puts $FileVar "Begin ModelPartData"
    #puts $FileVar "  // VARIABLE_NAME value"
    #puts $FileVar "End ModelPartData"
    #puts $FileVar ""
    #puts $FileVar ""

    ## Tables
    set TableId 0
    set TableDict [dict create]
    # Velocity
    ConstraintVectorTable FileVar TableId TableDict Velocity VELOCITY
    # Phi_Value
    set SolutionType [GiD_AccessValue get gendata Solution_Type]
    set SchemeType [GiD_AccessValue get gendata Scheme]

    if {$SolutionType eq "Steady"} {

	PressureTable FileVar TableId TableDict Phi_Value TEMPERATURE

    } else {

	if {$SchemeType eq "Implicit"} {

	    PressureTable FileVar TableId TableDict Phi_Value PHI_THETA

	} else {

	    PressureTable FileVar TableId TableDict Phi_Value TEMPERATURE

	}

    }
    # Face_Heat_Flux
    ScalarTable FileVar TableId TableDict Face_Heat_Flux FACE_HEAT_FLUX
    # Q_Source
    ScalarTable FileVar TableId TableDict Q_Source HEAT_FLUX

    puts $FileVar ""

    ## Properties
    set PropertyId 0
    set PropertyDict [dict create]
    # Body_Part
    set Groups [GiD_Info conditions Body_Part groups]

    for {set i 0} {$i < [llength $Groups]} {incr i} {
	incr PropertyId
	dict set PropertyDict [lindex [lindex $Groups $i] 1] $PropertyId
	puts $FileVar "Begin Properties $PropertyId"

	puts $FileVar "  CONDUCTIVITY [lindex [lindex $Groups $i] 3]"
	puts $FileVar "  SPECIFIC_HEAT [lindex [lindex $Groups $i] 4]"
	puts $FileVar "  DENSITY [lindex [lindex $Groups $i] 5]"
	puts $FileVar "  ABSORPTION_COEFFICIENT [lindex [lindex $Groups $i] 6]"
	puts $FileVar "End Properties"
	puts $FileVar ""
    }
    puts $FileVar ""


    ## Nodes
    set Nodes [GiD_Info Mesh Nodes]
    puts $FileVar "Begin Nodes"
    for {set i 0} {$i < [llength $Nodes]} {incr i 4} {
	# puts $FileVar "  [lindex $Nodes $i]  [lindex $Nodes [expr { $i+1 }]] [lindex $Nodes [expr { $i+2 }]] [lindex $Nodes [expr { $i+3 }]]"
	puts -nonewline $FileVar "  [lindex $Nodes $i]  "
	puts -nonewline $FileVar [format  "%.10f" [lindex $Nodes [expr { $i+1 }]]]
	puts -nonewline $FileVar " "
	puts -nonewline $FileVar [format  "%.10f" [lindex $Nodes [expr { $i+2 }]]]
	puts -nonewline $FileVar " "
	puts $FileVar [format  "%.10f" [lindex $Nodes [expr { $i+3 }]]]
    }
    puts $FileVar "End Nodes"
    puts $FileVar ""
    puts $FileVar ""

    ## Elements
    #set IsQuadratic [GiD_Info Project Quadratic]
    # Body_Part
    set Groups [GiD_Info conditions Body_Part groups]

    if {$SolutionType eq "Steady"} {

	for {set i 0} {$i < [llength $Groups]} {incr i} {
	    # Elements Property
	    set BodyElemsProp [dict get $PropertyDict [lindex [lindex $Groups $i] 1]]

	    # SteadyConvectionDiffusionFICElement2D3N
	    WriteElements FileVar [lindex $Groups $i] triangle SteadyConvectionDiffusionFICElement2D3N $BodyElemsProp Triangle2D3Connectivities
	    # SteadyConvectionDiffusionFICElement2D4N
	    WriteElements FileVar [lindex $Groups $i] quadrilateral SteadyConvectionDiffusionFICElement2D4N $BodyElemsProp Quadrilateral2D4Connectivities
	    # SteadyConvectionDiffusionFICElement3D4N
	    WriteElements FileVar [lindex $Groups $i] tetrahedra SteadyConvectionDiffusionFICElement3D4N $BodyElemsProp Quadrilateral2D4Connectivities
	    # SteadyConvectionDiffusionFICElement3D8N
	    WriteElements FileVar [lindex $Groups $i] hexahedra SteadyConvectionDiffusionFICElement3D8N $BodyElemsProp Hexahedron3D8Connectivities
	}

    } else {

	if {$SchemeType eq "Implicit"} {

	    for {set i 0} {$i < [llength $Groups]} {incr i} {
		# Elements Property
		set BodyElemsProp [dict get $PropertyDict [lindex [lindex $Groups $i] 1]]

		# TransientConvectionDiffusionFICElement2D3N
		WriteElements FileVar [lindex $Groups $i] triangle TransientConvectionDiffusionFICElement2D3N $BodyElemsProp Triangle2D3Connectivities
		# TransientConvectionDiffusionFICElement2D4N
		WriteElements FileVar [lindex $Groups $i] quadrilateral TransientConvectionDiffusionFICElement2D4N $BodyElemsProp Quadrilateral2D4Connectivities
		# TransientConvectionDiffusionFICElement3D4N
		WriteElements FileVar [lindex $Groups $i] tetrahedra TransientConvectionDiffusionFICElement3D4N $BodyElemsProp Quadrilateral2D4Connectivities
		# TransientConvectionDiffusionFICElement3D8N
		WriteElements FileVar [lindex $Groups $i] hexahedra TransientConvectionDiffusionFICElement3D8N $BodyElemsProp Hexahedron3D8Connectivities
	    }

	} else {

	    for {set i 0} {$i < [llength $Groups]} {incr i} {
		# Elements Property
		set BodyElemsProp [dict get $PropertyDict [lindex [lindex $Groups $i] 1]]

		# TransientConvectionDiffusionFICExplicitElement2D3N
		WriteElements FileVar [lindex $Groups $i] triangle TransientConvectionDiffusionFICExplicitElement2D3N $BodyElemsProp Triangle2D3Connectivities
		# TransientConvectionDiffusionFICExplicitElement2D4N
		WriteElements FileVar [lindex $Groups $i] quadrilateral TransientConvectionDiffusionFICExplicitElement2D4N $BodyElemsProp Quadrilateral2D4Connectivities
		# TransientConvectionDiffusionFICExplicitElement3D4N
		WriteElements FileVar [lindex $Groups $i] tetrahedra TransientConvectionDiffusionFICExplicitElement3D4N $BodyElemsProp Quadrilateral2D4Connectivities
		# TransientConvectionDiffusionFICExplicitElement3D8N
		WriteElements FileVar [lindex $Groups $i] hexahedra TransientConvectionDiffusionFICExplicitElement3D8N $BodyElemsProp Hexahedron3D8Connectivities
	    }

	}


    }

    puts $FileVar ""

    ## Conditions
    set ConditionId 0
    set ConditionDict [dict create]
    set Dim [GiD_AccessValue get gendata Domain_Size]
    # Face_Heat_Flux
    set Groups [GiD_Info conditions Face_Heat_Flux groups]
    if {$Dim eq 2} {
	# FluxCondition2D2N
	WriteFaceConditions FileVar ConditionId ConditionDict $Groups FluxCondition2D2N $PropertyDict
    } else {
	for {set i 0} {$i < [llength $Groups]} {incr i} {
	    set MyConditionList [list]
	    # FluxCondition3D3N
	    WriteTypeFaceConditions FileVar ConditionId MyConditionList [lindex $Groups $i] tetrahedra FluxCondition3D3N $PropertyDict
	    # FluxCondition3D4N
	    WriteTypeFaceConditions FileVar ConditionId MyConditionList [lindex $Groups $i] hexahedra FluxCondition3D4N $PropertyDict
	    dict set ConditionDict [lindex [lindex $Groups $i] 1] $MyConditionList
	}
    }

    puts $FileVar ""

    ## SubModelParts
    # Body_Part
    WriteElementSubmodelPart FileVar Body_Part
    # Velocity
    WriteConstraintSubmodelPart FileVar Velocity $TableDict
    # Phi_Value
    WriteConstraintSubmodelPart FileVar Phi_Value $TableDict
    # Face_Heat_Flux
    WriteLoadSubmodelPart FileVar Face_Heat_Flux $TableDict $ConditionDict
    # Q_Source
    WriteConstraintSubmodelPart FileVar Q_Source $TableDict

    close $FileVar

    return $TableDict
}
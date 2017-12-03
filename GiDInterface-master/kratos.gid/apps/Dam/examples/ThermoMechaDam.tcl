
proc ::Dam::examples::ThermoMechaDam {args} {
    
    DrawDamGeometry
    AssignGroupsDam
    AssignDamMeshSizes
    TreeAssignationDam

    GiD_Process 'Redraw
    GidUtils::UpdateWindow GROUPS
    GidUtils::UpdateWindow LAYER
}

proc Dam::examples::DrawDamGeometry {args} {
    
    Kratos::ResetModel
    GiD_Layers create Dam
    GiD_Layers edit to_use Dam

    # Geometry creation
    ## Points ##
    set coordinates [list 0 0 0 10 0 0 3 30 0 0 30 0 ]
    set damPoints [list ]
    foreach {x y z} $coordinates {
         lappend damPoints [GiD_Geometry create point append Dam $x $y $z]
    }

    ## Lines ##
    set damLines [list ]
    set initial [lindex $damPoints 0]
    foreach point [lrange $damPoints 1 end] {
        lappend damLines [GiD_Geometry create line append stline Dam $initial $point]
        set initial $point
    }
    lappend damLines [GiD_Geometry create line append stline Dam $initial [lindex $damPoints 0]]

     ## Surface ##
    GiD_Process Mescape Geometry Create NurbsSurface {*}$damLines escape escape

    #~ # Soil #
    GiD_Layers create Soil
    GiD_Layers edit to_use Soil
    
    #~ # Geometry creation
    #~ ## Points ##
    set soil_coordinates [list -5 0 0 -5 -5 0 15 -5 0 15 0 0 ]
    set soilPoints [list ]
    foreach {x y z} $soil_coordinates {
        lappend soilPoints [GiD_Geometry create point append Soil $x $y $z]
    }
    
    ## Lines ##
    set soilLines [list ]
    set initial [lindex $damPoints 0]
    foreach point [lrange $soilPoints 0 end] {
        lappend soilLines [GiD_Geometry create line append stline Soil $initial $point]
        set initial $point
    }
    lappend soilLines [GiD_Geometry create line append stline Soil $initial [lindex $damPoints 1]]
    
    
    lappend soilLines 1
    
    ## Surface ##
    GiD_Process Mescape Geometry Create NurbsSurface {*}$soilLines escape escape
    
    GiD_Process 'Zoom Frame
        
}

proc Dam::examples::AssignGroupsDam {args} {
    
    # Create the groups
    GiD_Groups create Dam
    GiD_Groups edit color Dam "#26d1a8ff"
    GiD_EntitiesGroups assign Dam surfaces 1
    
    GiD_Groups create Soil
    GiD_Groups edit color Soil "#e0210fff"
    GiD_EntitiesGroups assign Soil surfaces 2
    
    GiD_Groups create Displacement
    GiD_Groups edit color Displacement "#3b3b3bff"
    GiD_EntitiesGroups assign Displacement lines 7
      
    GiD_Groups create Initial
    GiD_Groups edit color Initial "#26d1a8ff"
    GiD_EntitiesGroups assign Initial surfaces {1 2}
    
    GiD_Groups create Bofang
    GiD_Groups edit color Bofang "#42eb71ff"
    GiD_EntitiesGroups assign Bofang lines {4 5}

    GiD_Groups create Uniform
    GiD_Groups edit color Uniform "#3b3b3bff"
    GiD_EntitiesGroups assign Uniform lines {3 2 9} 
    
    GiD_Groups create Thermal_Parameters_1
    GiD_Groups edit color Initial "#26d1a8ff"
    GiD_EntitiesGroups assign Thermal_Parameters_1 surfaces 1
    
    GiD_Groups create Thermal_Parameters_2
    GiD_Groups edit color Initial "#26d1a8ff"
    GiD_EntitiesGroups assign Thermal_Parameters_2 surfaces 2
    
    GiD_Groups create Hydrostatic
    GiD_Groups edit color Hydrostatic "#26d1a8fe"
    GiD_EntitiesGroups assign Hydrostatic lines {4 5}

}

proc Dam::examples::AssignDamMeshSizes {args} {
	
    set dam_mesh_size 0.25
    GiD_Process Mescape Meshing AssignSizes Surfaces $dam_mesh_size [GiD_EntitiesGroups get Dam surfaces] escape escape
    GiD_Process Mescape Meshing AssignSizes Surfaces $dam_mesh_size [GiD_EntitiesGroups get Soil surfaces] escape escape
    Kratos::BeforeMeshGeneration $dam_mesh_size
}

# Tree assign
proc Dam::examples::TreeAssignationDam {args} {

	set nd $::Model::SpatialDimension
    set root [customlib::GetBaseRoot]

    # Set Type of problem strategy set
    spdAux::SetValueOnTreeItem v "Thermo-Mechanical" DamTypeofProblem

    # Dam Part
    set damParts [spdAux::getRoute "DamParts"]
    set damNode [spdAux::AddConditionGroupOnXPath $damParts Dam]
    set props [list Element SmallDisplacementElement2D ConstitutiveLaw ThermalLinearElastic2DPlaneStrain Material "Concrete-Dam" DENSITY 2400 YOUNG_MODULUS 1.962e10 POISSON_RATIO 0.20 THERMAL_EXPANSION 1e-05]
    foreach {prop val} $props {
        set propnode [$damNode selectNodes "./value\[@n = '$prop'\]"]
        if {$propnode ne "" } {
            $propnode setAttribute v $val
        } else {
            W "Warning - Couldn't find property Dam $prop"
        }
    }
	
	#Soil Part
    set soilNode [spdAux::AddConditionGroupOnXPath $damParts Soil]
    set props_soil [list Element SmallDisplacementElement2D ConstitutiveLaw ThermalLinearElastic2DPlaneStrain Material Soil DENSITY 3000 YOUNG_MODULUS 4.9e10 POISSON_RATIO 0.25 THERMAL_EXPANSION 1e-05]
    foreach {prop val} $props_soil {
        set propnode [$soilNode selectNodes "./value\[@n = '$prop'\]"]
        if {$propnode ne "" } {
            $propnode setAttribute v $val
        } else {
            W "Warning - Couldn't find property Dam $prop"
        }
    }

	# Dirichlet Conditions
	
		# Displacements
		set damDirichletConditions [spdAux::getRoute "DamNodalConditions"]
		set displacement "$damDirichletConditions/condition\[@n='DISPLACEMENT'\]"
		set displacemnetnode [spdAux::AddConditionGroupOnXPath $displacement Displacement]
		
		# Surface Temperature 
		set initial "$damDirichletConditions/condition\[@n='INITIALTEMPERATURE'\]"
		set initialnode [spdAux::AddConditionGroupOnXPath $initial Initial]
		set props_initial [list is_fixed 0 value 7.5 ]
		foreach {prop val} $props_initial {
			 set propnode [$initialnode selectNodes "./value\[@n = '$prop'\]"]
			 if {$propnode ne "" } {
				  $propnode setAttribute v $val
			 } else {
				W "Warning - Couldn't find property Initial $prop"
			}
		}

		# Bofang Temperature
		set bofang "$damDirichletConditions/condition\[@n='BOFANGTEMPERATURE'\]"
		set bofangnode [spdAux::AddConditionGroupOnXPath $bofang Bofang]
		set props_bofang [list is_fixed 1 Gravity_Direction Y Reservoir_Bottom_Coordinate_in_Gravity_Direction 0.0 Surface_Temp 15.19 Bottom_Temp 9.35 Height_Dam 30.0 Temperature_Amplitude 6.51 Day_Ambient_Temp 201 Water_level 20.0 Outer_temp 10.0 Month 7 ]
		foreach {prop val} $props_bofang {
			 set propnode [$bofangnode selectNodes "./value\[@n = '$prop'\]"]
			 if {$propnode ne "" } {
				  $propnode setAttribute v $val
			 } else {
				W "Warning - Couldn't find property Bofang $prop"
			}
		}
		
		# Uniform Temperature
		set uniform "$damDirichletConditions/condition\[@n='INITIALTEMPERATURE'\]"
		set uniformnode [spdAux::AddConditionGroupOnXPath $uniform Uniform]
		set props_uniform [list is_fixed 1 value 10.0 ]
		foreach {prop val} $props_uniform {
			 set propnode [$uniformnode selectNodes "./value\[@n = '$prop'\]"]
			 if {$propnode ne "" } {
				  $propnode setAttribute v $val
			 } else {
				W "Warning - Couldn't find property Uniform $prop"
			}
		}
	
	
	# Thermal Load Conditions
	
		# Thermal Parameters 1
		set damThermalLoadConditions [spdAux::getRoute "DamThermalLoads"]
		set thermalparameter "$damThermalLoadConditions/condition\[@n='ThermalParameters2D'\]"
		set thermalparameternode1 [spdAux::AddConditionGroupOnXPath $thermalparameter Thermal_Parameters_1]
		
		# Thermal Parameters 2
		set thermalparameternode2 [spdAux::AddConditionGroupOnXPath $thermalparameter Thermal_Parameters_2]
		set props_thermal_2 [list ThermalDensity 3000 ]
		foreach {prop val} $props_thermal_2 {
			 set propnode [$thermalparameternode2 selectNodes "./value\[@n = '$prop'\]"]
			 if {$propnode ne "" } {
				  $propnode setAttribute v $val
			 } else {
				W "Warning - Couldn't find property Thermal_Parameters_2 $prop"
			}
		}
		
	# Load Conditions
   
		# Hydrostatic Load
        set damLoadConditions [spdAux::getRoute "DamLoads"]
		set hydro "$damLoadConditions/condition\[@n='HydroLinePressure2D'\]"
		set hydronode [spdAux::AddConditionGroupOnXPath $hydro Hydrostatic]
		set props_hydro [list Modify 0 Gravity_Direction Y Reservoir_Bottom_Coordinate_in_Gravity_Direction 0.0 Spe_weight 10000 Water_level 20.0]
		foreach {prop val} $props_hydro {
			 set propnode [$hydronode selectNodes "./value\[@n = '$prop'\]"]
			 if {$propnode ne "" } {
				  $propnode setAttribute v $val
			 } else {
				W "Warning - Couldn't find property Hydrostatic $prop"
			}
		}
		
	# Solution
	spdAux::SetValueOnTreeItem v "Days" DamTimeScale

	# Results
    set results [list REACTION No TEMPERATURE Yes POSITIVE_FACE_PRESSURE  Yes]]
    set nodal_path [spdAux::getRoute "NodalResults"]
    foreach {n v} $results {
        [$root selectNodes "$nodal_path/value\[@n = '$n'\]"] setAttribute v $v
    }
	

    spdAux::RequestRefresh


}


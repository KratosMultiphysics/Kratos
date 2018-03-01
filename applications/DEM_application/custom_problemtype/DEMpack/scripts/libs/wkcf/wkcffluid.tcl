###############################################################################
#
#    NAME: wkcffluid.tcl
#
#    PURPOSE: Useful procedures to work with the fluid application
#
#    QUANTECH ATZ-DEVELOPMENT DEPARTMENT
#
#    AUTHOR : G. Socorro
#
#    CREATED AT: 02/04/12
#
#    HISTORY:
#
#     3.3- 31/10/13-G. Socorro, change the proc GetDensityViscosityValues by the proc GetFluidMaterialProperties 
#     3.2- 14/07/13-G. Socorro, modify the proc WriteFluidIsSlipWallLawBC to write is-slip and walllaw BC
#     3.1- 17/06/13-G. Socorro, delete wmethod variable and all related procedures (*_m0,*_m1,*_m2) => now we are using only the new GiD groups
#     3.0- 12/04/13-G. Socorro, correct a bug in the proc WriteFluidInletNoSlipBC_m2 (2D case)
#     2.9- 22/03/13-G. Socorro, correct a bug in the proc WriteFluidInletNoSlipBC_m2 (using write_calc_data instead of [GiD_EntitiesGroups get $nsgroupid nodes])
#     2.8- 10/12/12-J. Garate,  PFEM PT dont need to write Density and Viscosity from WritePropertyAtNodes
#     2.7- 05/12/12-J. Garate,  PFEM Slip velocity format correction
#     2.6- 03/12/12-J. Garate,  Added Bulk Modulus on GetDensityViscosityValues return value. ::wkcf::WriteFluidProjectParameters for PFEM
#     2.5- 28/11/12-J. Garate,  Created ::wkcf::WriteFluidPFEMWallBC and ::wkcf::WriteFluidPFEMInletBC
#     2.4- 12/11/12-J. Garate,  Fixed some errors
#     2.3- 07/11/12-J. Garate,  Modification and adaptation on functions: WritePropertyAtNodes, GetDensityViscosityValues, WriteFluidBC, WriteFluidInletNoSlipBC
#                               WriteFluidFlagVariableBC, WriteFluidIsSlipBC, WriteFluidWallLawBC, WriteFluidDistanceBC, WriteOutLetPressureBC.
#                               Creation of functions using GiD_File fprintf $filechannel "%s" format
#     2.2- 05/10/12-G. Socorro, write density and viscosity variable for the LevelSet in the projectparameter.py
#     2.1- 04/10/12-G. Socorro, write variable using the format node_id 0 node_value
#     2.0- 04/10/12-G. Socorro, update the proc WritePropertyAtNodes_m1 to write the LevelSet variable at nodal level
#     1.9- 03/10/12-G. Socorro, add the proc WriteFluidDistanceBC and write free surface option in the projectparameter file
#     1.8- 28/09/12-G. Socorro, set VolumeOutput = True by default in 2D problems
#     1.7- 23/09/12-G. Socorro, update the proc WriteFluidProjectParameters to write turbulence properties in the 
#                               that wfsmethod=0
#     1.6- 22/09/12-G. Socorro, update the proc WriteFluidProjectParameters to write turbulence properties
#     1.5- 20/08/12-G. Socorro, correct a bug when write is_structure for 3D problems
#     1.4- 25/07/12-G. Socorro, add VolumeOutput to 3D problem
#     1.3- 24/07/12-G. Socorro, update the procedure WriteFluidIsSlipBC to write Y_Wall nodaldata
#     1.2- 22/07/12-G. Socorro, modify the Is-Slip BC, correct a bug when write inlet BC 
#     1.1- 08/06/12-G. Socorro, add a new variable Use_slip_conditions when use is-slip or wall-law conditions
#     1.0- 04/06/12-G. Socorro, set Laplacian form=1 for the fractional step case
#     0.9- 13/05/12-G. Socorro, modify the proc WriteFluidIsSlipBC to use the dictionary ctbclink to link conditions 
#                               with the No-Slip boundary condition
#     0.8- 10/05/12-G. Socorro, correct a bug for 3D case when write IS-SLIP boundary condition
#     0.7- 08/05/12-G. Socorro, update the condition IS-SLIP to delete the fixity value
#     0.6- 06/05/12-G. Socorro, modify the procs WriteFluidProjectParameters and WriteFluidSolvers to write using the 
#                               new and old fluid solver format, modify the proc WriteOutLetPressureBC to write using (write_calc_data)
#     0.5- 05/05/12-G. Socorro, modify the procs WriteFluidFlagVariableBC and ::wkcf::WriteFluidInletNoSlipBC_m1 to write the properties using the fast method (write_calc_data) 
#     0.4- 04/05/12-G. socorro, modify the proc WritePropertyAtNodes to write the properties using the fast method (write_calc_data) 
#     0.3- 23/04/12-G. Socorro, add the proc WriteFluidInletNoSlipBC
#     0.2- 10/04/12-G. Socorro, modify the proc WriteFluidProjectParameters
#     0.1- 02/04/12-G. Socorro, create a base source code from wkcf.tcl
#
###############################################################################

proc ::wkcf::WritePropertyAtNodes {AppId} {
    # ABSTRACT: Write some properties at the nodal level for Fluid application
    variable dprops
    variable filechannel

    # For debug
    if {!$::wkcf::pflag} {
	set inittime [clock seconds]
    }

    set cproperty "dv"
    # Free surface
    set cxpath "$AppId//c.AnalysisData//i.FreeSurface"
    set FreeSurface [::xmlutils::setXml $cxpath $cproperty]
    # wa "FreeSurface:$FreeSurface"
    
    # PFEM
    set cxpath "$AppId//c.AnalysisData//i.FluidApproach"
    set PFEM [::xmlutils::setXml $cxpath $cproperty]
    # wa "cxpath:$cxpath"
    # wa "FluidApproach:$PFEM"

    set flag [expr {($FreeSurface eq "No") && ([info exists dprops($AppId,AllKElemId)]) && ([llength $dprops($AppId,AllKElemId)]) && ($PFEM ne "PFEM-Lagrangian")}]
    # wa "flag:$flag"
    # Check for all defined kratos elements
    if {$flag} {
	set kxpath "Materials"
	# Write fluid material properties at nodal level
	set matdict [dict create]
	# Get the material properties dictionary
	set matdict [::wkcf::GetFluidMaterialProperties $AppId]
	# wa "matdict:$matdict"
	# Get the fluid type
	set NonNewtonianFluid "No"
	if {[dict exists $matdict "NonNewtonianFluid"]} {
	    set NonNewtonianFluid [dict get $matdict "NonNewtonianFluid"]
	}
	# Case: Newtoniam fluid
	set plist [list "Viscosity" "Density"]
	set kwlist [list]
	set vplist [list]
	#if {$NonNewtonianFluid eq "Yes"} {
	#    lappend plist "BinghamSmoother" "YieldStress" "PowerLawN" "PowerLawK" "GelStrength"
	#}
	# wa "plist:$plist"
	# update temporal list 
	foreach pid $plist {
	    lappend kwlist [::xmlutils::getKKWord $kxpath $pid "kkword"]
	    lappend vplist [dict get $matdict $pid]
	}
	#if {$NonNewtonianFluid eq "Yes"} {
	    # Replace viscosity value by the plastic viscosity value
	    #lset vplist 0 [dict get $matdict "PlasticViscosity"]
	#} else {
	    # Replace power law N with 1.0
	    #lset vplist 4 1.0
	    #Replace YieldStress with 0.0
	    #lset vplist 3 0.0
	    #Replace GelStrength with 0.0
	    #lset vplist 6 0.0
	#}
	# wa "kwlist:$kwlist\nvplist:$vplist"

	set cpropid "0"        

	# Write the group nodal properties
	foreach celemid $dprops($AppId,AllKElemId) {
	    # Check for all defined group identifier for this element
	    set cflag [expr {([info exists dprops($AppId,KElem,$celemid,AllGroupId)]) && ([llength $dprops($AppId,KElem,$celemid,AllGroupId)])}]
	    if {$cflag} {
		# For all defined group identifier for this element
		foreach cgroupid $dprops($AppId,KElem,$celemid,AllGroupId) {                    
		    # Write the nodal material properties
		    if {[GiD_EntitiesGroups get $cgroupid nodes -count]} {
		        foreach kw $kwlist cval $vplist {
		            GiD_File fprintf $filechannel "%s" "Begin NodalData $kw \/\/ GUI group identifier: $cgroupid"
		            foreach node_id [GiD_EntitiesGroups get $cgroupid nodes] {
		                GiD_File fprintf $filechannel "%10i %4i %20.10f" $node_id $cpropid $cval
		            }
		            GiD_File fprintf $filechannel "%s " "End NodalData"
		            GiD_File fprintf $filechannel ""
		        }
		    }
		}
	    }
	}
    }

    # Try to write the levelset properties
    set contid "PorousZones"
    set kxpath "Applications//$AppId"
  
    # For free surface
    if {$FreeSurface =="Yes"} {
	if {[info exists dprops($AppId,AllPorousZonesTypeId)] && [llength $dprops($AppId,AllPorousZonesTypeId)]} {
	    # Get the application root identifier    
	    set rootdataid $AppId
	    set cpropid "0"       
	    set cxpath "$rootdataid//c.SolutionStrategy//c.[list ${contid}]//i.UseErgunEquation"
	    set UseErgunEquation [::xmlutils::setXml $cxpath $cproperty]
	    # wa "UseErgunEquation:$UseErgunEquation"

	    # Write the group nodal properties
	    foreach czonetypeid $dprops($AppId,AllPorousZonesTypeId) {
		set cproplist [list]
		if {(($czonetypeid eq "ErgunEquationNo") && ($UseErgunEquation eq "No"))} {
		    set cproplist [list "PorosityValue" "LinearDarcyCoefficient" "NonLinearDarcyCoefficient"]
		} elseif {(($czonetypeid eq "ErgunEquationYes") && ($UseErgunEquation eq "Yes"))} {
		    set cproplist [list "PorosityValue" "DiameterValue"]
		}
		if {![llength $cproplist]} {
		    continue
		}
		# wa "czonetypeid:$czonetypeid cproplist:$cproplist"
		# Check for all defined group identifier for this zone type
		if {([info exists dprops($AppId,$contid,$czonetypeid,AllGroupId)]) && ([llength $dprops($AppId,$contid,$czonetypeid,AllGroupId)])} {
		    # For all defined group identifier for this zone type
		    foreach cgroupid $dprops($AppId,$contid,$czonetypeid,AllGroupId) {
		        # wa "cgroupid:$cgroupid"
		        if {[info exists dprops($AppId,$contid,$czonetypeid,$cgroupid,GProps)] && [llength $dprops($AppId,$contid,$czonetypeid,$cgroupid,GProps)]} {
		            foreach pid $cproplist pvalue $dprops($AppId,$contid,$czonetypeid,$cgroupid,GProps) {
		                # wa "pid:$pid pvalue:$pvalue"
		                # Write the current variable value for this group
		                if {[GiD_EntitiesGroups get $cgroupid nodes -count]} {
		                    set vkword [::xmlutils::getKKWord $kxpath "$pid" "kkword"]
		                    # wa "vkword:$vkword"
		                    GiD_File fprintf $filechannel "%s" "Begin NodalData $vkword \/\/ GUI group identifier: $cgroupid"
		                    foreach node_id [GiD_EntitiesGroups get $cgroupid nodes] {
		                        GiD_File fprintf $filechannel "%10i %4i %s" $node_id $cpropid $pvalue
		                    }
		                    GiD_File fprintf $filechannel "%s" "End NodalData"
		                    GiD_File fprintf $filechannel "%s" ""
		                }
		            }
		        }
		    }
		}
	    }
	}
    }

    # For debug
    if {!$::wkcf::pflag} {
	set endtime [clock seconds]
	set ttime [expr $endtime-$inittime]
	# WarnWinText "endtime:$endtime ttime:$ttime"
	WarnWinText "Write property at nodes: [::KUtils::Duration $ttime]"
    }
}

proc ::wkcf::GetFluidMaterialProperties {AppId {what All} {key ""}} {
    # ABSTRACT : Return the fluid material properties values
    variable dprops
    
    set cprop [list]
    
    # Set the kratos element list flag
    set flag [expr {([info exists dprops($AppId,AllKElemId)]) && ([llength $dprops($AppId,AllKElemId)])}]
    # wa "flag:$flag"
    
    if {$flag} {
	# Create the material properties dictionary
	set pdict [dict create]
	set cproplist [list "NonNewtonianFluid" "Density" "Viscosity" "BinghamSmoother" "YieldStress" "PowerLawN" "PowerLawK" "GelStrength"]
	set dfvcproplist [list No 0.0 0.0 0.0 0.0]
	# Init the dictionary
	foreach pid $cproplist dfv $dfvcproplist {
	    dict set pdict $pid $dfv
	}
	# wa "pdict:$pdict"
	# foreach PropertyId $dprops($AppId,GKProps,AllPropertyId)         
	foreach PropertyId $dprops($AppId,AllKPropertyId) {   
	    # wa "PropertyId:$PropertyId"
	    # Get the material identifier for this property 
	    set MatId $dprops($AppId,Property,$PropertyId,MatId) 
	    # Get the group identifier
	    set GroupId $dprops($AppId,Property,$PropertyId,GroupId)
	    # Get all material properties
	    set mpxpath "[::KMat::findMaterialParent $MatId]//m.[list ${MatId}]"
	    # WarnWinText "mpxpath:$mpxpath"
	    # Get the material properties
	    foreach pid $cproplist {
		if {$pid =="Density"} {
		    set xpath "c.General"
		} else {
		    set xpath "c.Fluid"
		}
		# Get the current value for this properties
		set cvalue [lindex [::KMat::getMaterialProperties "p" "$mpxpath//$xpath//p.[list $pid]"] 0 1]
		# Format some properties
		if {$pid ne "NonNewtonianFluid"} {
		    set $pid [GiD_FormatReal "%10.5e" $cvalue]
		} else {
		    set $pid $cvalue 
		}
		# Create/update the fluid material properties dictionary
		dict set pdict $pid [set $pid]
	    }
	    # Only the first property
	    break 
	}
	# Update cprop
	if {[llength [dict keys $pdict]]} {
	    switch -exact -- $what {
	    "All" {
		set cprop $pdict
	    }
	    "PropertyId" {
		if {[dict exists $pdict $key]} {
		  set cprop [dict get $pdict $key]
		} 
	    }
	    }
	    unset pdict
	}
    }
    # wa "cprop:$cprop"
    return $cprop
}


proc ::wkcf::WriteFluidBC {AppId inletvelglist noslipglist flagvariablelist kwordlist} {
    # ABSTRACT: Write the fluid boundary conditions
    variable ndime; variable dprops; variable filechannel

    # WarnWinText "inletvelglist:$inletvelglist\nnoslipglist:$noslipglist\nflagvariablelist:$flagvariablelist\nkwordlist:$kwordlist"
   
    if {([llength $inletvelglist]) || ([llength $noslipglist])} {
	# For debug
	if {!$::wkcf::pflag} {
	    set inittime [clock seconds]
	}

    # Map Inlet-NoSlip => Use no-slip values at share nodes
    set icondid "InletVelocity"; set nscondid "No-Slip"
    set cpropid "1"
    set xitem [lindex $kwordlist 0]
    set yitem [lindex $kwordlist 1]
    set zitem [lindex $kwordlist 2]
    
    if {[llength $noslipglist]} {
	# Write all no-slip condition
	
	# For each group in the no-slip condition
	foreach nsgroupid $noslipglist {
	    lassign $dprops($AppId,BC,$nscondid,$nsgroupid,GProps) nsx nsxval nsy nsyval nsz nszval

	    # X component
	    if {$nsx} {
		  if { [GiD_EntitiesGroups get $nsgroupid nodes -count] } {
		    GiD_File fprintf $filechannel "%s" "Begin NodalData $xitem \/\/ No-slip condition GUI group identifier: $nsgroupid"
		    foreach node_id [GiD_EntitiesGroups get $nsgroupid nodes] {
		        GiD_File fprintf $filechannel "%10i %8i %10.5e" $node_id $cpropid $nsxval
		    }
		    GiD_File fprintf $filechannel "%s" "End NodalData"
		    GiD_File fprintf $filechannel ""
		  }		  
	    }
	    
	    # Y component
	    if {$nsy} {
		if { [GiD_EntitiesGroups get $nsgroupid nodes -count] } {
		    GiD_File fprintf $filechannel "%s" "Begin NodalData $yitem \/\/ No-slip condition GUI group identifier: $nsgroupid"
		    foreach node_id [GiD_EntitiesGroups get $nsgroupid nodes] {
		        GiD_File fprintf $filechannel "%10i %8i %10.5e" $node_id $cpropid $nsyval
		    }
		    GiD_File fprintf $filechannel "%s" "End NodalData"
		    GiD_File fprintf $filechannel ""
		}	    
        }
	    
	    # Z component
	    if {$ndime =="3D"} {
		   if {$nsz} {
		    if { [GiD_EntitiesGroups get $nsgroupid nodes -count] } {
		    GiD_File fprintf $filechannel "%s" "Begin NodalData $zitem \/\/ No-slip condition GUI group identifier: $nsgroupid"
		    foreach node_id [GiD_EntitiesGroups get $nsgroupid nodes] {
		        GiD_File fprintf $filechannel "%10i %8i %10.5e" $node_id $cpropid $nszval
		    }
		    GiD_File fprintf $filechannel "%s" "End NodalData"
		    GiD_File fprintf $filechannel ""
		    }
		  }
	    }
	    
	    
	    variable property_number
	    
	    if {$nsx || $nsy || $nsz} {
            set property_number [expr $property_number + 1 ]
            GiD_File fprintf $filechannel "Begin Properties $property_number // GUI property identifier: $nsgroupid"
            GiD_File fprintf $filechannel "IMPOSED_PRESSURE 0"
            if {$nsx} {
                GiD_File fprintf $filechannel "IMPOSED_VELOCITY_X 1"
                GiD_File fprintf $filechannel "IMPOSED_VELOCITY_X_VALUE $nsxval"
            } else {
                GiD_File fprintf $filechannel "IMPOSED_VELOCITY_X 0"
            }
            if {$nsy} {
                GiD_File fprintf $filechannel "IMPOSED_VELOCITY_Y 1"
                GiD_File fprintf $filechannel "IMPOSED_VELOCITY_Y_VALUE $nsyval"
            } else {
                GiD_File fprintf $filechannel "IMPOSED_VELOCITY_Y 0"
            }
            if {$nsz} {
                GiD_File fprintf $filechannel "IMPOSED_VELOCITY_Z 1"
                GiD_File fprintf $filechannel "IMPOSED_VELOCITY_Z_VALUE $nszval"
            } else {
                GiD_File fprintf $filechannel "IMPOSED_VELOCITY_Z 0"
            }	    	  
            GiD_File fprintf $filechannel "IS_SLIP 0"
            GiD_File fprintf $filechannel "End Properties"
            
            GiD_File fprintf $filechannel "Begin Mesh $property_number // GUI property identifier: $nsgroupid"
            GiD_File fprintf $filechannel "Begin MeshNodes"
            foreach node_id [GiD_EntitiesGroups get $nsgroupid nodes] {
                GiD_File fprintf $filechannel "%10i" $node_id
            }                   
            GiD_File fprintf $filechannel "End MeshNodes"
            GiD_File fprintf $filechannel "End Mesh" 
            
        }   
	    
	}
    }
    
    # Use first the inlet
    if {([llength $inletvelglist]) && ([llength $noslipglist])} {
	# Check to match node identifier
	set condmatch [dict create]
	# For each group in the no-slip condition
	foreach nsgroupid $noslipglist {
	    lassign $dprops($AppId,BC,$nscondid,$nsgroupid,GProps) cx cxval cy cyval cz czval
	    # wa "nsgroupid:$nsgroupid cx:$cx cxval:$cxval cy:$cy cyval:$cyval cz:$cz czval:$czval"
	    if {[GiD_EntitiesGroups get $nsgroupid nodes -count] } {
		# For each node in the no-slip bc update the condmatch 
		foreach nsnodeid [GiD_EntitiesGroups get $nsgroupid nodes] {
		    dict set condmatch $nsnodeid [list $cx $cy $cz]
		}
	    }
	}
	
	# For all inlet velocity group identifier
	set ixcomp ""; set iycomp ""; set izcomp ""
	foreach igroupid $inletvelglist {
	    lassign $dprops($AppId,BC,$icondid,$igroupid,GProps) ix ixval iy iyval iz izval
	    # wa "igroupid:$igroupid ix:$ix iy:$iy iz:$iz"
	    # Set the inlet format dictionary
	    if {[GiD_EntitiesGroups get $igroupid nodes -count]} {
		# 3D problems
		if {$ndime =="3D"} {
		    # For each node in the inlet bc ckeck to write this node
		    foreach inodeid [GiD_EntitiesGroups get $igroupid nodes] {
		        # WarnWinText "inodeid:$inodeid"
		        # Check that this node identifier exists in the dictionary
		        if {[dict exists $condmatch $inodeid]} {
		            # Get the properties
		            lassign [dict get $condmatch $inodeid] nsx nsy nsz
		            # X component => Check x flag
		            if {($ix) && ($nsx=="0")} {
		                # Write this node identifier
		                append ixcomp "[format "%8i%8i" $inodeid $cpropid]   [GiD_FormatReal "%10.5e" $ixval]\n"
		            }
		            # Y component => Check y flag
		            if {($iy) && ($nsy=="0")} {
		                # Write this node identifier
		                append iycomp "[format "%8i%8i" $inodeid $cpropid]   [GiD_FormatReal "%10.5e" $iyval]\n"
		            }
		            # Z component => Check z flag
		            if {($iz) && ($nsz=="0")} {
		                # Write this node identifier
		                append izcomp "[format "%8i%8i" $inodeid $cpropid]   [GiD_FormatReal "%10.5e" $izval]\n"
		            }
		        } else {
		            # Write this node identifier
		            # X component => Check x flag
		            if {$ix} {
		                # Write this node identifier
		                append ixcomp "[format "%8i%8i" $inodeid $cpropid]   [GiD_FormatReal "%10.5e" $ixval]\n"
		            }
		            # Y component => Check y flag
		            if {$iy} {
		                # Write this node identifier
		                append iycomp "[format "%8i%8i" $inodeid $cpropid]   [GiD_FormatReal "%10.5e" $iyval]\n"
		            }
		            # Z component => Check z flag
		            if {$iz} {
		                # Write this node identifier
		                append izcomp "[format "%8i%8i" $inodeid $cpropid]   [GiD_FormatReal "%10.5e" $izval]\n"
		            }
		        }
		    }
		    # Write this group identifier
		    if {[string length $ixcomp]} {
		        GiD_File fprintf $filechannel "%s" "Begin NodalData $xitem \/\/ Inlet velocity condition GUI group identifier: $igroupid"
		        GiD_File fprintf $filechannel "%s" $ixcomp
		        GiD_File fprintf $filechannel "%s" "End NodalData"
		        GiD_File fprintf $filechannel ""
		    }
		    if {[string length $iycomp]} {
		        GiD_File fprintf $filechannel "%s" "Begin NodalData $yitem \/\/ Inlet velocity condition GUI group identifier: $igroupid"
		        GiD_File fprintf $filechannel "%s" $iycomp
		        GiD_File fprintf $filechannel "%s" "End NodalData"
		        GiD_File fprintf $filechannel ""
		    }
		    if {[string length $izcomp]} {
		        GiD_File fprintf $filechannel "%s" "Begin NodalData $zitem \/\/ Inlet velocity condition GUI group identifier: $igroupid"
		        GiD_File fprintf $filechannel "%s" $izcomp
		        GiD_File fprintf $filechannel "%s" "End NodalData"
		        GiD_File fprintf $filechannel ""
		    }
		    
		    if {[string length $ixcomp] || [string length $iycomp] || [string length $izcomp]} {
                set property_number [expr $property_number + 1 ]
                GiD_File fprintf $filechannel "Begin Properties $property_number // GUI property identifier: $igroupid"
                GiD_File fprintf $filechannel "IMPOSED_PRESSURE 0"
                GiD_File fprintf $filechannel "IMPOSED_VELOCITY_X 1"
                GiD_File fprintf $filechannel "IMPOSED_VELOCITY_X_VALUE $ixval"           
                GiD_File fprintf $filechannel "IMPOSED_VELOCITY_Y 1"
                GiD_File fprintf $filechannel "IMPOSED_VELOCITY_Y_VALUE $iyval"           
                GiD_File fprintf $filechannel "IMPOSED_VELOCITY_Z 1"
                GiD_File fprintf $filechannel "IMPOSED_VELOCITY_Z_VALUE $izval"
                GiD_File fprintf $filechannel "IS_SLIP 0"
                GiD_File fprintf $filechannel "End Properties"
                                    
                GiD_File fprintf $filechannel "Begin Mesh $property_number // GUI property identifier: $igroupid"
                GiD_File fprintf $filechannel "Begin MeshNodes"
                foreach node_id [GiD_EntitiesGroups get $igroupid nodes] {
                    GiD_File fprintf $filechannel "%10i" $node_id
                }                   
                GiD_File fprintf $filechannel "End MeshNodes"    
                GiD_File fprintf $filechannel "End Mesh"
		    }
		    
		    
		} elseif {$ndime =="2D"} {
		    # For each node in the inlet bc ckeck to write this node
		    foreach inodeid [GiD_EntitiesGroups get $igroupid nodes] {
		        # WarnWinText "inodeid:$inodeid"
		        # Check that this node identifier exists in the dictionary
		        if {[dict exists $condmatch $inodeid]} {
		            # Get the properties
		            lassign [dict get $condmatch $inodeid] nsx nsy nsz
		            # X component => Check x flag
		            if {($ix) && ($nsx=="0")} {
		                # Write this node identifier
		                lappend ixcomp $inodeid $cpropid $ixval
		            }
		            # Y component => Check y flag
		            if {($iy) && ($nsy=="0")} {
		                # Write this node identifier
		                lappend iycomp $inodeid $cpropid $iyval
		            }
		        } else {
		            # Write this node identifier
		            # X component => Check x flag
		            if {$ix} {
		                # Write this node identifier
		                lappend ixcomp $inodeid $cpropid $ixval
		            }
		            # Y component => Check y flag
		            if {$iy} {
		                # Write this node identifier
		                lappend iycomp $inodeid $cpropid $iyval
		            }
		        }
		    }
		    # wa "ixcomp:$ixcomp iycomp:$iycomp"
		    # Write this group identifier
		    if {[string length $ixcomp]} {
		        GiD_File fprintf $filechannel "%s" "Begin NodalData $xitem \/\/ Inlet velocity condition GUI group identifier: $igroupid"
		        foreach { inodeid cpropid ixval } $ixcomp {
		            GiD_File fprintf $filechannel "%8i %8i %10.5e" $inodeid $cpropid $ixval
		        }
		        GiD_File fprintf $filechannel "%s" "End NodalData"
		        GiD_File fprintf $filechannel ""
		    }
		    if {[string length $iycomp]} {
		        GiD_File fprintf $filechannel "%s" "Begin NodalData $yitem \/\/ Inlet velocity condition GUI group identifier: $igroupid"
		        foreach { inodeid cpropid iyval } $iycomp {
		            GiD_File fprintf $filechannel "%8i %8i %10.5e" $inodeid $cpropid $iyval
		        }
		        GiD_File fprintf $filechannel "%s" "End NodalData"
		        GiD_File fprintf $filechannel ""
		    }
		}
	    }
	    # Reset ixcomp, iycomp and zcomp
	    set ixcomp ""; set iycomp ""; set izcomp ""
	}
	
	# unset dictionary variable
	unset condmatch

    } else {
	# Write the inlet boundary condition properties

	# For each group in the inlet condition
	foreach igroupid $inletvelglist {
	    lassign $dprops($AppId,BC,$icondid,$igroupid,GProps) ix ixval iy iyval iz izval

	    # X component
	    if {$ix} {
		if { [GiD_EntitiesGroups get $igroupid nodes -count] } {
		    GiD_File fprintf $filechannel "%s" "Begin NodalData $xitem \/\/ Inlet velocity condition GUI group identifier: $igroupid"
		    foreach node_id [GiD_EntitiesGroups get $igroupid nodes] {
		        GiD_File fprintf $filechannel "%10i %8i %10.5e" $node_id $cpropid $ixval
		    }
		    GiD_File fprintf $filechannel "%s" "End NodalData"
		    GiD_File fprintf $filechannel "%s" ""
		}
	    }
	    
	    # Y component
	    if {$iy} {
		if { [GiD_EntitiesGroups get $igroupid nodes -count] } {
		    GiD_File fprintf $filechannel "%s" "Begin NodalData $yitem \/\/ Inlet velocity condition GUI group identifier: $igroupid"
		    foreach node_id [GiD_EntitiesGroups get $igroupid nodes] {
		        GiD_File fprintf $filechannel "%10i %8i %10.5e" $node_id $cpropid $iyval
		    }
		    GiD_File fprintf $filechannel "%s" "End NodalData"
		    GiD_File fprintf $filechannel "%s" ""
		}
	    }
	    
	    # Z component
	    if {$ndime =="3D"} {
		if {$iz} {
		    if { [GiD_EntitiesGroups get $igroupid nodes -count] } {
		        GiD_File fprintf $filechannel "%s" "Begin NodalData $zitem \/\/ Inlet velocity condition GUI group identifier: $igroupid"
		        foreach node_id [GiD_EntitiesGroups get $igroupid nodes] {
		            GiD_File fprintf $filechannel "%10i %8i %10.5e" $node_id $cpropid $izval
		        }
		        GiD_File fprintf $filechannel "%s" "End NodalData"
		        GiD_File fprintf $filechannel "%s" ""
		    }
		}
	    }
	    if { $ix || $iy ||$iz } {
                variable property_number
                set property_number [expr $property_number + 1 ]
                GiD_File fprintf $filechannel "Begin Properties $property_number // GUI property identifier: $igroupid"
                GiD_File fprintf $filechannel "IMPOSED_PRESSURE 0"
                GiD_File fprintf $filechannel "IMPOSED_VELOCITY_X $ix"
                GiD_File fprintf $filechannel "IMPOSED_VELOCITY_X_VALUE $ixval"           
                GiD_File fprintf $filechannel "IMPOSED_VELOCITY_Y $iy"
                GiD_File fprintf $filechannel "IMPOSED_VELOCITY_Y_VALUE $iyval"           
                GiD_File fprintf $filechannel "IMPOSED_VELOCITY_Z $iz"
                GiD_File fprintf $filechannel "IMPOSED_VELOCITY_Z_VALUE $izval"
                GiD_File fprintf $filechannel "IS_SLIP 0"
                GiD_File fprintf $filechannel "End Properties"
                                    
                GiD_File fprintf $filechannel "Begin Mesh $property_number // GUI property identifier: $igroupid"
                GiD_File fprintf $filechannel "Begin MeshNodes"
                foreach node_id [GiD_EntitiesGroups get $igroupid nodes] {
                    GiD_File fprintf $filechannel "%10i" $node_id
                }                   
                GiD_File fprintf $filechannel "End MeshNodes"    
                GiD_File fprintf $filechannel "End Mesh"
            }
	}
    }

	# For debug
	if {!$::wkcf::pflag} {
	    set endtime [clock seconds]
	    set ttime [expr $endtime-$inittime]
	    # WarnWinText "endtime:$endtime ttime:$ttime"
	    WarnWinText "Write fluid inlet-no-slip boundary conditions: [::KUtils::Duration $ttime]"
		}
		
	    }
    
    # Write Flag-variable and is_boundary nodal data conditions
    if {[llength $flagvariablelist]} {
	
	# For debug
	if {!$::wkcf::pflag} {
	    set inittime [clock seconds]
	    }
	::wkcf::WriteFluidFlagVariableBC $AppId $flagvariablelist
	
	# For debug
	if {!$::wkcf::pflag} {
	    set endtime [clock seconds]
	    set ttime [expr $endtime-$inittime]
	    # WarnWinText "endtime:$endtime ttime:$ttime"
	    WarnWinText "Write fluid flag variable boundary conditions: [::KUtils::Duration $ttime]"
	}
    }
}

proc ::wkcf::WriteFluidPFEMWallBC {AppId ccondid kwordlist} {
    variable ndime; variable dprops; variable filechannel

    # Map Inlet-NoSlip => Use no-slip values at share nodes
    set cpropid "1"
    
    # PFEM Linear velocity XYZ
    set xitem [lindex $kwordlist 0]
    set yitem [lindex $kwordlist 1]
    set zitem [lindex $kwordlist 2]
    
    # PFEM Angular velocity XYZ
    set xitem2 [lindex $kwordlist 3]
    set yitem2 [lindex $kwordlist 4]
    set zitem2 [lindex $kwordlist 5]
    
    # For each group in the inlet condition
    foreach cgroupid $dprops($AppId,BC,$ccondid,AllGroupId) {
	# wa "cgroupid:$cgroupid"
	# Get the condition properties
	lassign $dprops($AppId,BC,$ccondid,$cgroupid,GProps) ixval iyval izval ixval2 iyval2 izval2

	# X Linear Velocity
	if { [GiD_EntitiesGroups get $cgroupid nodes -count] } {
	    GiD_File fprintf $filechannel "%s" "Begin NodalData $xitem \/\/ PFEM Linear velocity condition GUI group identifier: $cgroupid"
	    foreach node_id [GiD_EntitiesGroups get $cgroupid nodes] {
		GiD_File fprintf $filechannel "%10i %8i %10.5e" $node_id $cpropid $ixval
	    }
	    GiD_File fprintf $filechannel "End NodalData"
	    GiD_File fprintf $filechannel ""
	}
	
	# Y Linear Velocity
	if { [GiD_EntitiesGroups get $cgroupid nodes -count] } {
	    GiD_File fprintf $filechannel "%s" "Begin NodalData $yitem \/\/ PFEM Linear velocity condition GUI group identifier: $cgroupid"
	    foreach node_id [GiD_EntitiesGroups get $cgroupid nodes] {
		GiD_File fprintf $filechannel "%10i %8i %10.5e" $node_id $cpropid $iyval
	    }
	    GiD_File fprintf $filechannel "End NodalData"
	    GiD_File fprintf $filechannel ""
	}
	
	# Z Linear Velocity
	if {$ndime =="3D"} {
	    if { [GiD_EntitiesGroups get $cgroupid nodes -count] } {
		GiD_File fprintf $filechannel "%s" "Begin NodalData $zitem \/\/ PFEM Linear velocity condition GUI group identifier: $cgroupid"
		foreach node_id [GiD_EntitiesGroups get $cgroupid nodes] {
		    GiD_File fprintf $filechannel "%10i %8i %10.5e" $node_id $cpropid $izval
		}
		GiD_File fprintf $filechannel "End NodalData"
		GiD_File fprintf $filechannel ""
	    }
	}
	
	# X Angular Velocity
	if { [GiD_EntitiesGroups get $cgroupid nodes -count] } {
	    GiD_File fprintf $filechannel "%s" "Begin NodalData $xitem2 \/\/ PFEM Angular velocity condition GUI group identifier: $cgroupid"
	    foreach node_id [GiD_EntitiesGroups get $cgroupid nodes] {
		GiD_File fprintf $filechannel "%10i %8i %10.5e" $node_id $cpropid $ixval2
	    }
	    GiD_File fprintf $filechannel "End NodalData"
	    GiD_File fprintf $filechannel ""
	}
	
	# Y Angular Velocity
	if { [GiD_EntitiesGroups get $cgroupid nodes -count] } {
	    GiD_File fprintf $filechannel "%s" "Begin NodalData $yitem2 \/\/ PFEM Angular velocity condition GUI group identifier: $cgroupid"
	    foreach node_id [GiD_EntitiesGroups get $cgroupid nodes] {
		GiD_File fprintf $filechannel "%10i %8i %10.5e" $node_id $cpropid $iyval2
	    }
	    GiD_File fprintf $filechannel "End NodalData"
	    GiD_File fprintf $filechannel ""
	}
	
	# Z Angular Velocity
	if {$ndime =="3D"} {
	    if { [GiD_EntitiesGroups get $cgroupid nodes -count] } {
		GiD_File fprintf $filechannel "%s" "Begin NodalData $zitem2 \/\/ PFEM Angular velocity condition GUI group identifier: $cgroupid"
		foreach node_id [GiD_EntitiesGroups get $cgroupid nodes] {
		    GiD_File fprintf $filechannel "%10i %8i %10.5e" $node_id $cpropid $izval2
		}
		GiD_File fprintf $filechannel "End NodalData"
		GiD_File fprintf $filechannel ""
	    }
	}
    }
}

proc ::wkcf::WriteFluidPFEMInletBC {AppId ccondid kwordlist} {
    variable ndime; variable dprops; variable filechannel

    set cpropid "1"
    
    # PFEM Inlet velocity XYZ
    set xitem [lindex $kwordlist 0]
    set yitem [lindex $kwordlist 1]
    set zitem [lindex $kwordlist 2]

    # For each group in the inlet condition
    foreach cgroupid $dprops($AppId,BC,$ccondid,AllGroupId) {
	# wa "cgroupid:$cgroupid"
	# Get the condition properties
	lassign $dprops($AppId,BC,$ccondid,$cgroupid,GProps) ixvar ixval iyvar iyval izvar izval
	# msg "x: $ixvar $ixval y: $iyvar $iyval z: $izvar $izval"
	# X Linear Velocity
	if {$ixvar} {
	    if { [GiD_EntitiesGroups get $cgroupid nodes -count] } {
		GiD_File fprintf $filechannel "%s" "Begin NodalData $xitem \/\/ PFEM Inlet velocity condition GUI group identifier: $cgroupid"
		foreach node_id [GiD_EntitiesGroups get $cgroupid nodes] {
		    GiD_File fprintf $filechannel "%10i %8i %10.5e" $node_id $cpropid $ixval
		}
		GiD_File fprintf $filechannel "End NodalData"
		GiD_File fprintf $filechannel ""
	    }
	}
	
	# Y Linear Velocity
	if {$iyvar} {
	    if { [GiD_EntitiesGroups get $cgroupid nodes -count] } {
		GiD_File fprintf $filechannel "%s" "Begin NodalData $yitem \/\/ PFEM Inlet velocity condition GUI group identifier: $cgroupid"
		foreach node_id [GiD_EntitiesGroups get $cgroupid nodes] {
		    GiD_File fprintf $filechannel "%10i %8i %10.5e" $node_id $cpropid $iyval
		}
		GiD_File fprintf $filechannel "End NodalData"
		GiD_File fprintf $filechannel ""
	    }
	}
	
	# Z Linear Velocity
	if {$ndime =="3D"} {
	    if {$izvar} {
		if { [GiD_EntitiesGroups get $cgroupid nodes -count] } {
		    GiD_File fprintf $filechannel "%s" "Begin NodalData $zitem \/\/ PFEM Inlet velocity condition GUI group identifier: $cgroupid"
		    foreach node_id [GiD_EntitiesGroups get $cgroupid nodes] {
		        GiD_File fprintf $filechannel "%10i %8i %10.5e" $node_id $cpropid $izval
		    }
		    GiD_File fprintf $filechannel "End NodalData"
		    GiD_File fprintf $filechannel ""
		}
	    }
	}
    }
}

proc ::wkcf::WriteFluidFlagVariableBC {AppId flagvariablelist} {
    # ABSTRACT: Write the flag variable boundary condition
    variable dprops
    variable filechannel

    # WarnWinText "flagvariablelist:$flagvariablelist"
    # For nodes with many flag variable defined flag of level two have the priority over flag of level one
    set flagvarcondid "Flag-Variable"
    set cpropid "0"
    set isbcpropid "1"

    # Write the flag condition
    set flag1 0
    set fvitem "FLAG-VARIABLE"
    set isbitem "IS_BOUNDARY"
    
    # For each group in the flag-variable condition
    # Create a dict for all nodes with flag equal to 2
    set flagvar2 [dict create]
    
    foreach cgroupid $flagvariablelist {
	lassign $dprops($AppId,BC,$flagvarcondid,$cgroupid,GProps) flagval
	# WarnWinText "cgroupid:$cgroupid flagval:$flagval"

	if {$flagval=="2"} {
	    # Write this group identifier
	    set flag1 1

	    # Flag-Variable
	    if { [GiD_EntitiesGroups get $cgroupid nodes -count] } {
		GiD_File fprintf $filechannel "%s" "Begin NodalData $fvitem \/\/ Flag-Variable condition GUI group identifier: $cgroupid"
		foreach node_id [GiD_EntitiesGroups get $cgroupid nodes] {
		    GiD_File fprintf $filechannel "%10i %8i %8i" $node_id $cpropid $flagval
		    # Update the flagvar2 dictionary
		    dict set flagvar2 $nodeid $cgroupid
		}
		GiD_File fprintf $filechannel "%s" "End NodalData"
		GiD_File fprintf $filechannel ""
	    }

	    # is_boundary
	    if { [GiD_EntitiesGroups get $cgroupid nodes -count] } {
		GiD_File fprintf $filechannel "%s" "Begin NodalData $isbitem \/\/ is_boundary associated with Flag-Variable condition GUI group identifier: $cgroupid"
		foreach node_id [GiD_EntitiesGroups get $cgroupid nodes] {
		    GiD_File fprintf $filechannel "%10i %8i %8i" $node_id $cpropid $isbcpropid
		}
		GiD_File fprintf $filechannel "%s" "End NodalData"
		GiD_File fprintf $filechannel ""
	    }
	}
    }

    # Write all group with flag-variable equal to 1
    if {$flag1} {
	set fvcomp ""; set isbcomp ""
	# For each group in the flag-variable condition
	foreach cgroupid $flagvariablelist {
	    lassign $dprops($AppId,BC,$flagvarcondid,$cgroupid,GProps) flagval
	    # WarnWinText "cgroupid:$cgroupid GProps:$GProps"
	    
	    if {$flagval=="1"} {
		if { [GiD_EntitiesGroups get $cgroupid nodes -count] } {
		    foreach nodeid [GiD_EntitiesGroups get $cgroupid nodes] {
		        if {![dict exists $flagvar2 $nodeid]} {
		            append fvcomp "[format "%10i%8i%8i" $nodeid $cpropid $flagval]\n"
		            append isbcomp "[format "%10i%8i%8i" $nodeid $cpropid $isbcpropid]\n"
		        }
		    }
		    # Write this group identifier
		    # Flag-Variable
		    if {[string length $fvcomp]} {
		        GiD_File fprintf $filechannel "%s" "Begin NodalData $fvitem \/\/ Flag-Variable condition GUI group identifier: $cgroupid"
		        GiD_File fprintf $filechannel "%s" "[string trimright ${fvcomp}]"
		        GiD_File fprintf $filechannel "%s" "End NodalData"
		        GiD_File fprintf $filechannel ""
		    }
		    # is_boundary
		    if {[string length $isbcomp]} {
		        GiD_File fprintf $filechannel "%s" "Begin NodalData $isbitem \/\/ is_boundary associated with Flag-Variable condition GUI group identifier: $cgroupid"
		        GiD_File fprintf $filechannel "%s" "[string trimright ${isbcomp}]"
		        GiD_File fprintf $filechannel "%s" "End NodalData"
		        GiD_File fprintf $filechannel ""
		    }
		    # Reset components
		    set fvcomp ""; set isbcomp ""
		}
	    }
	}
    }
   
    # unset the dict for all nodes with flag equal to 2
    unset flagvar2
}


proc ::wkcf::WriteFluidIsSlipWallLawBC {AppId ccondid kwordlist } {
    # ABSTRACT: Write is-slip/walllaw boundary conditions => Conditional data
    variable dprops;   variable ndime
    variable ctbclink; variable filechannel

    # wa "AppId:$AppId ccondid:$ccondid kwordlist:$kwordlist"
    # For debug
    if {!$::wkcf::pflag} {
	set inittime [clock seconds]
    }
    # Set the keyword values
    set isstructurekw [lindex $kwordlist 0]
    set isywallkw [lindex $kwordlist 1]

    set state 0
    # Variable to control when use slip conditions
    set dprops($AppId,UseSlipConditions) 0


    # For all defined group identifier inside this condition type
    foreach cgroupid $dprops($AppId,BC,$ccondid,AllGroupId) {
	# wa "cgroupid:$cgroupid"
	# Get the condition properties
	lassign $dprops($AppId,BC,$ccondid,$cgroupid,GProps) activateval ConstantValue
	# wa "activateval:$activateval ConstantValue:$ConstantValue"
	if {$ndime == "2D"} {
	    set GiDElemType "Linear"
	    if {[GiD_EntitiesGroups get $cgroupid elements -count -element_type $GiDElemType]} {
		set dprops($AppId,UseSlipConditions) 1
		GiD_File fprintf $filechannel "%s" "Begin ConditionalData $isstructurekw // GUI $ccondid condition group identifier: $cgroupid"
		foreach elem_id [GiD_EntitiesGroups get $cgroupid elements -element_type $GiDElemType] {
		    # set nodes [lrange [GiD_Mesh get element $elem_id] 3 end] 
		    # wa "elemid:$elem_id cfixval:$cfixval nodei:$nodei nodej:$nodej"
		    if {[dict exists $ctbclink $elem_id]} {
		    # Check that exists this element in the dictionary with the condition indentifier links
		        set condid [dict get $ctbclink $elem_id]
		        GiD_File fprintf $filechannel "%10d %10d" $condid $activateval
		    }
		}
		GiD_File fprintf $filechannel "%s" "End ConditionalData"
		GiD_File fprintf $filechannel "%s" ""
		
		# Write wall_y values
		GiD_File fprintf $filechannel "%s" " Begin NodalData $isywallkw // GUI $ccondid condition group identifier: $cgroupid"
		foreach nodeid [GiD_EntitiesGroups get $cgroupid nodes] {
		    GiD_File fprintf $filechannel "%10d %4d %10g" $nodeid $state $ConstantValue
		}
		GiD_File fprintf $filechannel "%s" "End NodalData"
		GiD_File fprintf $filechannel "" 
	
	    }
	} elseif {$ndime == "3D"} {
	    set GiDElemType "Triangle"
	    if {[GiD_EntitiesGroups get $cgroupid elements -count -element_type $GiDElemType]} {
		set dprops($AppId,UseSlipConditions) 1
		GiD_File fprintf $filechannel "%s" "Begin ConditionalData $isstructurekw // GUI $ccondid condition group identifier: $cgroupid"
		foreach elem_id [GiD_EntitiesGroups get $cgroupid elements -element_type $GiDElemType] {
		    #set nodes [lrange [GiD_Mesh get element $elem_id] 3 end] 
		    # wa "elemid:$elem_id nodei:[lindex $nodes 0] nodej:[lindex $nodes 1] nodek:[lindex $nodes 2]"
		    # Check that exists this element in the dictionary with the condition indentifier links
		    # wa [dict keys $ctbclink]
		    if {[dict exists $ctbclink $elem_id]} {
		        set condid [dict get $ctbclink $elem_id]
		        # msg "$condid $activateval"
		        GiD_File fprintf $filechannel "%10d %10d" $condid $activateval
		    }
		}
		GiD_File fprintf $filechannel "End ConditionalData"
		GiD_File fprintf $filechannel ""
		
		# Write wall_y values
		GiD_File fprintf $filechannel "%s" " Begin NodalData $isywallkw // GUI $ccondid condition group identifier: $cgroupid"
		foreach nodeid [GiD_EntitiesGroups get $cgroupid nodes] {
		    # msg "$nodeid $state $ConstantValue"
		    GiD_File fprintf $filechannel "%10d %4d %10g" $nodeid $state $ConstantValue
		}
		GiD_File fprintf $filechannel "%s" "End NodalData"
		GiD_File fprintf $filechannel ""
	    } 
	}
	variable property_number
    set property_number [expr $property_number + 1 ]
    GiD_File fprintf $filechannel "Begin Properties $property_number // GUI property identifier: $cgroupid"
    GiD_File fprintf $filechannel "IMPOSED_PRESSURE 0"
    GiD_File fprintf $filechannel "IMPOSED_VELOCITY_X 0"
    GiD_File fprintf $filechannel "IMPOSED_VELOCITY_Y 0"
    GiD_File fprintf $filechannel "IMPOSED_VELOCITY_Z 0"
    GiD_File fprintf $filechannel "IS_SLIP 1"
    GiD_File fprintf $filechannel "End Properties"
                        
    GiD_File fprintf $filechannel "Begin Mesh $property_number // GUI property identifier: $cgroupid"
    GiD_File fprintf $filechannel "Begin MeshNodes"
    foreach node_id [GiD_EntitiesGroups get $cgroupid nodes] {
        GiD_File fprintf $filechannel "%10i" $node_id
    }                   
    GiD_File fprintf $filechannel "End MeshNodes"    
    GiD_File fprintf $filechannel "End Mesh"
    }

    # For debug
    if {!$::wkcf::pflag} {
	set endtime [clock seconds]
	set ttime [expr $endtime-$inittime]
	# WarnWinText "endtime:$endtime ttime:$ttime"
	WarnWinText "Write fluid wall law/is-slip boundary conditions: [::KUtils::Duration $ttime]"
    }
}

proc ::wkcf::WritePFEMLagrangianFluidFixedWallBC {AppId ccondid kwordlist } {
    # ABSTRACT: Write fixed wall boundary conditions => Conditional data
    variable dprops;   variable ndime
    variable ctbclink; variable filechannel

    # wa "AppId:$AppId ccondid:$ccondid kwordlist:$kwordlist"
    # For debug
    if {!$::wkcf::pflag} {
	set inittime [clock seconds]
    }
    # Set the keyword values
    set isstructurekw [lindex $kwordlist 0]
    set DispX [lindex $kwordlist 1]
    set DispY [lindex $kwordlist 2]
    set DispZ [lindex $kwordlist 3]

    set state 1
    set activateval "1.0"
    set fvalue "0.0"

    # For all defined group identifier inside this condition type
    foreach cgroupid $dprops($AppId,BC,$ccondid,AllGroupId) {
	# wa "cgroupid:$cgroupid"
	# Get the condition properties
	lassign $dprops($AppId,BC,$ccondid,$cgroupid,GProps) DispXValue DispYValue DispZValue
	# wa "DispXValue:$DispXValue DispYValue:$DispYValue DispZValue:$DispZValue"
	if {$ndime == "2D"} {
	    set GiDElemType "Linear"
	    if {[GiD_EntitiesGroups get $cgroupid elements -count -element_type $GiDElemType]} {
		GiD_File fprintf $filechannel "%s" "Begin ConditionalData $isstructurekw // GUI $ccondid condition group identifier: $cgroupid"
		foreach elem_id [GiD_EntitiesGroups get $cgroupid elements -element_type $GiDElemType] {
		    if {[dict exists $ctbclink $elem_id]} {
			# Check that exists this element in the dictionary with the condition indentifier links
		        set condid [dict get $ctbclink $elem_id]
		        GiD_File fprintf $filechannel "%10d %10g %10g" $condid $activateval $fvalue
		    }
		}
		GiD_File fprintf $filechannel "%s" "End ConditionalData"
		GiD_File fprintf $filechannel "%s" ""
		
		# Write displacement values
		foreach kw [list $DispX $DispY] cvalue [list $DispXValue $DispYValue] {
		    GiD_File fprintf $filechannel "%s" " Begin NodalData $kw // GUI $ccondid condition group identifier: $cgroupid"
		    foreach nodeid [GiD_EntitiesGroups get $cgroupid nodes] {
			GiD_File fprintf $filechannel "%10d %4d %10g" $nodeid $state $cvalue
		    }
		    GiD_File fprintf $filechannel "%s" "End NodalData"
		    GiD_File fprintf $filechannel "" 
		}
	    }
	} elseif {$ndime == "3D"} {
	    set GiDElemType "Triangle"
	    if {[GiD_EntitiesGroups get $cgroupid elements -count -element_type $GiDElemType]} {
		GiD_File fprintf $filechannel "%s" "Begin ConditionalData $isstructurekw // GUI $ccondid condition group identifier: $cgroupid"
		foreach elem_id [GiD_EntitiesGroups get $cgroupid elements -element_type $GiDElemType] {
		    # Check that exists this element in the dictionary with the condition indentifier links
		    # wa [dict keys $ctbclink]
		    if {[dict exists $ctbclink $elem_id]} {
		        set condid [dict get $ctbclink $elem_id]
		        # wa "$condid $activateval"
		        GiD_File fprintf $filechannel "%10d %10g %10g" $condid $activateval $fvalue
		    }
		}
		GiD_File fprintf $filechannel "End ConditionalData"
		GiD_File fprintf $filechannel ""
	
		# Write displacement values
		foreach kw [list $DispX $DispY $DispZ] cvalue [list $DispXValue $DispYValue $DispZValue] {
		    GiD_File fprintf $filechannel "%s" " Begin NodalData $kw // GUI $ccondid condition group identifier: $cgroupid"
		    foreach nodeid [GiD_EntitiesGroups get $cgroupid nodes] {
			GiD_File fprintf $filechannel "%10d %4d %10g" $nodeid $state $cvalue
		    }
		    GiD_File fprintf $filechannel "%s" "End NodalData"
		    GiD_File fprintf $filechannel "" 
		}
	    } 
	}
    }

    # For debug
    if {!$::wkcf::pflag} {
	set endtime [clock seconds]
	set ttime [expr $endtime-$inittime]
	# WarnWinText "endtime:$endtime ttime:$ttime"
	WarnWinText "Write PFEM fixed wall boundary conditions: [::KUtils::Duration $ttime]"
    }
}

proc ::wkcf::WriteFluidDistanceBC {AppId ccondid kwordlist} {
    # ABSTRACT: Write distance boundary conditions => Nodal data
    variable dprops; variable filechannel
    
    # For debug
    if {!$::wkcf::pflag} {
	set inittime [clock seconds]
}


    set cpropid "0"       
    # For all defined group identifier inside this condition type
    foreach cgroupid $dprops($AppId,BC,$ccondid,AllGroupId) {
	# wa "cgroupid:$cgroupid"
	# Get the condition properties
	lassign $dprops($AppId,BC,$ccondid,$cgroupid,GProps) cvalue
	# wa "cvalue:$cvalue"
	if {[GiD_EntitiesGroups get $cgroupid nodes -count]} {
	    GiD_File fprintf $filechannel "%s" "Begin NodalData $kwordlist // GUI distance condition group identifier: $cgroupid"
	    foreach node_id [GiD_EntitiesGroups get $cgroupid nodes] {
		GiD_File fprintf $filechannel "%10i %5i %10.5f" $node_id $cpropid $cvalue
	    }
	    GiD_File fprintf $filechannel "%s" "End NodalData"
	    GiD_File fprintf $filechannel ""
	}
    }

    # For debug
    if {!$::wkcf::pflag} {
	set endtime [clock seconds]
	set ttime [expr $endtime-$inittime]
	# WarnWinText "endtime:$endtime ttime:$ttime"
	WarnWinText "Write fluid distance boundary conditions: [::KUtils::Duration $ttime]"
    }
}

proc ::wkcf::WriteOutLetPressureBC {AppId ccondid kwordlist} {
    # ASBTRACT: Write outlet pressure boundary condition
    variable dprops;   variable filechannel
    variable property_number

    # For debug
    if {!$::wkcf::pflag} {
	set inittime [clock seconds]
    }
    
    set kitem [lindex $kwordlist 0]
    
    # For all defined group identifier inside this condition type
    foreach cgroupid $dprops($AppId,BC,$ccondid,AllGroupId) {
	# Get the condition properties
	lassign $dprops($AppId,BC,$ccondid,$cgroupid,GProps) fixval pval
	# WarnWinText "fixval:$fixval pval:$pval"
	
	set nodes [GiD_EntitiesGroups get $cgroupid nodes]
	# Fix x
	if {$fixval =="1"} {
	    
	    GiD_File fprintf $filechannel "%s" "Begin NodalData $kitem"
	    foreach nodeid $nodes {
		GiD_File fprintf $filechannel "%10i %4i %10f" $nodeid $fixval $pval
	    }
	    GiD_File fprintf $filechannel "End NodalData"
	    GiD_File fprintf $filechannel "" 
	    
	    
	    set property_number [expr $property_number + 1 ]
        GiD_File fprintf $filechannel "Begin Properties $property_number // GUI property identifier: $cgroupid"
        GiD_File fprintf $filechannel "IMPOSED_PRESSURE 1"
        GiD_File fprintf $filechannel "PRESSURE $pval"
        GiD_File fprintf $filechannel "IMPOSED_VELOCITY_X 0"
        GiD_File fprintf $filechannel "IMPOSED_VELOCITY_Y 0"
        GiD_File fprintf $filechannel "IMPOSED_VELOCITY_Z 0"
        GiD_File fprintf $filechannel "IS_SLIP 0"
        GiD_File fprintf $filechannel "End Properties"
                            
        GiD_File fprintf $filechannel "Begin Mesh $property_number // GUI property identifier: $cgroupid"
        GiD_File fprintf $filechannel "Begin MeshNodes"
        foreach node_id [GiD_EntitiesGroups get $cgroupid nodes] {
            GiD_File fprintf $filechannel "%10i" $node_id
        }                   
        GiD_File fprintf $filechannel "End MeshNodes"    
        GiD_File fprintf $filechannel "End Mesh"
	}
    }

    # For debug
    if {!$::wkcf::pflag} {
	set endtime [clock seconds]
	set ttime [expr $endtime-$inittime]
	# WarnWinText "endtime:$endtime ttime:$ttime"
	WarnWinText "Write fluid wall law boundary conditions: [::KUtils::Duration $ttime]"
    }
}

proc ::wkcf::WriteFluidProjectParameters {AppId fileid PDir} {
    variable ndime; variable dprops
    
    # FluidApproach
    set rootid "$AppId"
    set cxpath "$rootid//c.AnalysisData//i.FluidApproach"
    set FluidApproach [::xmlutils::setXml $cxpath "dv"]
    if {$FluidApproach == "Eulerian"} {
	::wkcf::WriteEulerianFluidProjectParameters $AppId $fileid $PDir
    } elseif {$FluidApproach == "PFEM-Lagrangian"} {
	::wkcf::WritePFEMLagrangianFluidProjectParameters $AppId $fileid $PDir
    }
}

proc ::wkcf::WriteEulerianFluidProjectParameters {AppId fileid PDir} {
    variable ndime; variable dprops
    variable ActiveAppList

    set trailing_spaces  ""
    # Write fluid solver method
    # 0 => Old format
    # 1 => New format
    set wfsmethod 1
    set rootid "$AppId"
    
    # Kratos key word xpath
    set kxpath "Applications/$rootid"
    
    set domain_size 2
    if {$ndime =="2D"} {
	set domain_size 2
    } elseif {$ndime =="3D"} {
	set domain_size 3
    }
    
    # Domain size
    puts $fileid "domain_size = $domain_size"
    puts $fileid ""
    
    # Get others properties
    set cproperty "dv"
    
    # Fluid type
    if {"DEM" in $ActiveAppList} {
        set analysis_data_path [GetAnalysisDataContainer]
        set options_container [GetOptionsContainer]
    } else {
        set analysis_data_path "AnalysisData"
    }
    
    set cxpath "$rootid//c.$analysis_data_path//i.FluidType"
    set FluidType [::xmlutils::setXml $cxpath $cproperty]
    
    # WarnWinText "FluidType:$FluidType"
    if {$FluidType =="Compressible"} {
	
    } elseif {$FluidType =="Incompressible"} {
	
	# Free surface
	set cxpath "$rootid//c.AnalysisData//i.FreeSurface"
	set FreeSurface [::xmlutils::setXml $cxpath $cproperty]
	# WarnWinText "FreeSurface:$FreeSurface"
	if {$FreeSurface =="No"} {
	    # Solver type
	    set cxpath "$rootid//c.$analysis_data_path//i.SolverType"
	    set SolverType [::xmlutils::setXml $cxpath $cproperty]

	    set cxpath "$rootid//c.SolutionStrategy//c.ParallelType//i.ParallelSolutionType"
	    set parallel_type [::xmlutils::setXml $cxpath $cproperty]

	    # WarnWinText "SolverType:$SolverType"
	    # Get the kratos keyword
	    set ckword [::xmlutils::getKKWord $kxpath $SolverType]
	    # WarnWinText "ckword:$ckword"
		set trailing_spaces  "    "
		# Fluid solver configuration
		puts $fileid ""
		puts $fileid "# Fluid solver configuration"
		puts $fileid "class FluidSolverConfiguration:"

		switch -exact -- $SolverType {
		    "ElementBased" { 

		        if { ${parallel_type} == "OpenMP" } {
		          puts $fileid "${trailing_spaces}solver_type =  \"vms_fractional_step_solver\"" 
		        } else {
		          puts $fileid "${trailing_spaces}solver_type =  \"trilinos_vms_fs_fluid_solver\"" 
		        }
		    }
		    "Monolithic" { 
		        if { ${parallel_type} == "OpenMP" } {
		          puts $fileid "${trailing_spaces}solver_type =  \"vms_monolithic_solver\"" 
		        } else {
		          puts $fileid "${trailing_spaces}solver_type =  \"trilinos_vms_monolithic_solver\"" 
		        }
		    }
		}
		
		puts $fileid "    domain_size = $domain_size"
		
		# Get the turbulence properties
		set cxpath "$rootid//c.AnalysisData//i.TurbulenceModel"
		set TurbulenceModel [::xmlutils::setXml $cxpath $cproperty]
		# WarnWinText "TurbulenceModel:$TurbulenceModel"
		if {$TurbulenceModel eq "Off"} {
		    puts $fileid "${trailing_spaces}TurbulenceModel = \"None\""
		} elseif {$TurbulenceModel eq "Smagorinsky-Lilly"} {
		    puts $fileid "${trailing_spaces}TurbulenceModel = \"$TurbulenceModel\""
		    # Get the smagorinsky-lilly constant
		    set cxpath "$rootid//c.AnalysisData//i.SmagorinskyConstant"
		    set SmagorinskyConstant [::xmlutils::setXml $cxpath $cproperty]
		    # WarnWinText "SmagorinskyConstant:$SmagorinskyConstant"
		    puts $fileid "${trailing_spaces}SmagorinskyConstant = $SmagorinskyConstant"
		} elseif {$TurbulenceModel eq "Spalart-Allmaras"} {
		    puts $fileid "${trailing_spaces}TurbulenceModel = \"$TurbulenceModel\""
		# Get the value of the turbulence viscosity
		    set cxpath "$rootid//c.AnalysisData//i.TurbulentViscosity"
		    set TurbulentViscosity [::xmlutils::setXml $cxpath $cproperty]
		    # wa "TurbulentViscosity:$TurbulentViscosity"
		    puts $fileid "${trailing_spaces}TurbulentViscosity = $TurbulentViscosity"
		    
		    # Try to get the group-mesh link  
		    # SA_wall_group_ids = [1, 5, 3]
		    # Get the values
		    set basexpath "$rootid//c.AnalysisData//c.Spalart-AllmarasGroupId${ndime}"
		    set gproplist [::xmlutils::setXmlContainerIds $basexpath]
		    # wa "gproplist:$gproplist"
		    if {[llength $gproplist]} {
		        set meshidlist [list]
		        foreach cgroupid $gproplist {
		            # Get the group properties
		            set cxpath "${basexpath}//c.[list ${cgroupid}]//c.MainProperties"
		            set allgprop [::xmlutils::setXmlContainerPairs $cxpath "" "dv"]
		            # wa "allgprop:$allgprop"
		            if {[llength $allgprop]} {
		                set Activate [lindex $allgprop 0 1]
		                # wa "Activate:$Activate"
		                if {$Activate} {
		                    if {[info exists dprops($AppId,Mesh,$cgroupid,MeshIdGroup)]} {
		                        set MeshIdGroup $dprops($AppId,Mesh,$cgroupid,MeshIdGroup)
		                        # wa "MeshIdGroup:$MeshIdGroup"
		                        if {$MeshIdGroup !=""} {
		                            append meshidlist "$MeshIdGroup,"
		                        }
		                    }
		                }
		            }
		        }
		        # wa "meshidlist:$meshidlist"
		        set findcomma [string last "," $meshidlist]
		        if {$findcomma !="-1"} {
		            set meshidlist [string range $meshidlist 0 end-1]
		            append meshidlist "\]"
		            set endmeshidlist "\[${meshidlist}"
		            puts $fileid "${trailing_spaces}SA_wall_group_ids = $endmeshidlist"
		        }
		    }
		}
		
		# Monolithic,PressureSplitting,ElementBased,EdgeBased
		puts $fileid ""
		switch -exact -- $SolverType {
		    "ElementBased" {
		        # Fractional step options => ElementBased                
		        # Solution strategy
		        # Linear solvers
		        # Velocity
		        ::wkcf::WriteLinearSolvers $rootid $fileid "Velocity" $wfsmethod $trailing_spaces "velocity_linear_solver_config"
		        # Pressure
		        ::wkcf::WriteLinearSolvers $rootid $fileid "Pressure" $wfsmethod $trailing_spaces "pressure_linear_solver_config" 

		        puts $fileid "$trailing_spaces"
		        puts $fileid "${trailing_spaces}#convergence criteria settings"
		        set cxpath "$rootid//c.SolutionStrategy//c.Advanced//i.RelativeVelocityTolerance"
		        set cvalue [::xmlutils::setXml $cxpath $cproperty]
		        puts $fileid "${trailing_spaces}vel_toll = $cvalue"

		        set cxpath "$rootid//c.SolutionStrategy//c.Advanced//i.RelativePressureTolerance"
		        set cvalue [::xmlutils::setXml $cxpath $cproperty]
                        puts $fileid "${trailing_spaces}press_toll = $cvalue"		    
		    }
		    "PressureSplitting" {
		        # Pressure splitting
		        # Solution strategy
		        # Linear solvers
		        # Velocity
		        ::wkcf::WriteLinearSolvers $rootid $fileid "Velocity" $wfsmethod $trailing_spaces 
		        puts $fileid ""
		        # Pressure
		        ::wkcf::WriteLinearSolvers $rootid $fileid "Pressure" $wfsmethod $trailing_spaces
		    }
		    "Monolithic" {
		        # Monolithic
		        
		        # Solution strategy
		        # Linear solvers
		        # Velocity
		        ::wkcf::WriteLinearSolvers $rootid $fileid "Monolithic" $wfsmethod $trailing_spaces "linear_solver_config" 
		        # Write relative and absolute tolerances
		        puts $fileid "$trailing_spaces"
		        puts $fileid "${trailing_spaces}#convergence criteria settings"
		        set ctlist [list "RelativeVelocityTolerance" "AbsoluteVelocityTolerance" "RelativePressureTolerance" "AbsolutePressureTolerance"]
		        foreach cv $ctlist {
		            set cxpath "$rootid//c.SolutionStrategy//c.Advanced//i.[list ${cv}]"
		            set cvalue [::xmlutils::setXml $cxpath $cproperty]
		            set ckword [::xmlutils::getKKWord $kxpath $cv]
		            puts $fileid "${trailing_spaces}$ckword = $cvalue"
		        }
		    }
		    
		    
		}
  
	    
	    
	    
	    
	    

		
# 		# Get the turbulence properties
# 		set cxpath "$rootid//c.AnalysisData//i.TurbulenceModel"
# 		set TurbulenceModel [::xmlutils::setXml $cxpath $cproperty]
# 		# WarnWinText "TurbulenceModel:$TurbulenceModel"
# 		if {$TurbulenceModel eq "Off"} {
# 		    puts $fileid "TurbulenceModel = \"None\""
# 		} elseif {$TurbulenceModel eq "Smagorinsky-Lilly"} {
# 		    puts $fileid "TurbulenceModel = \"$TurbulenceModel\""
# 		    # Get the smagorinsky-lilly constant
# 		    set cxpath "$rootid//c.AnalysisData//i.SmagorinskyConstant"
# 		    set SmagorinskyConstant [::xmlutils::setXml $cxpath $cproperty]
# 		    # WarnWinText "SmagorinskyConstant:$SmagorinskyConstant"
# 		    puts $fileid "SmagorinskyConstant = $SmagorinskyConstant"
# 		} elseif {$TurbulenceModel eq "Spalart-Allmaras"} {
# 		    puts $fileid "TurbulenceModel = \"$TurbulenceModel\""
# 		    # Get the value of the turbulence viscosity
# 		    set cxpath "$rootid//c.AnalysisData//i.TurbulentViscosity"
# 		    set TurbulentViscosity [::xmlutils::setXml $cxpath $cproperty]
# 		    # wa "TurbulentViscosity:$TurbulentViscosity"
# 		    puts $fileid "TurbulentViscosity = $TurbulentViscosity"
# 		    
# 		    # Try to get the group-mesh link  
# 		    # SA_wall_group_ids = [1, 5, 3]
# 		    # Get the values
# 		    set basexpath "$rootid//c.AnalysisData//c.Spalart-AllmarasGroupId${ndime}"
# 		    set gproplist [::xmlutils::setXmlContainerIds $basexpath]
# 		    # wa "gproplist:$gproplist"
# 		    if {[llength $gproplist]} {
# 		        set meshidlist [list]
# 		        foreach cgroupid $gproplist {
# 		            # Get the group properties
# 		            set cxpath "${basexpath}//c.[list ${cgroupid}]//c.MainProperties"
# 		            set allgprop [::xmlutils::setXmlContainerPairs $cxpath "" "dv"]
# 		            # wa "allgprop:$allgprop"
# 		            if {[llength $allgprop]} {
# 		                set Activate [lindex $allgprop 0 1]
# 		                # wa "Activate:$Activate"
# 		                if {$Activate} {
# 		                    if {[info exists dprops($AppId,Mesh,$cgroupid,MeshIdGroup)]} {
# 		                        set MeshIdGroup $dprops($AppId,Mesh,$cgroupid,MeshIdGroup)
# 		                        # wa "MeshIdGroup:$MeshIdGroup"
# 		                        if {$MeshIdGroup !=""} {
# 		                            append meshidlist "$MeshIdGroup,"
# 		                        }
# 		                    }
# 		                }
# 		            }
# 		        }
# 		        # wa "meshidlist:$meshidlist"
# 		        set findcomma [string last "," $meshidlist]
# 		        if {$findcomma !="-1"} {
# 		            set meshidlist [string range $meshidlist 0 end-1]
# 		            append meshidlist "\]"
# 		            set endmeshidlist "\[${meshidlist}"
# 		            puts $fileid "SA_wall_group_ids = $endmeshidlist"
# 		        }
# 		    }
# 		}
	    
	    # Get the divergence clearance step
	    set cxpath "$rootid//c.ProblemParameters//i.DivergenceCleareanceStep"
	    set DivergenceCleareanceStep [::xmlutils::setXml $cxpath $cproperty]
	    # Get the kratos keyword
	    set DivergenceCleareanceStepKW [::xmlutils::getKKWord $kxpath "DivergenceCleareanceStep"]
        puts $fileid "${trailing_spaces}$DivergenceCleareanceStepKW = $DivergenceCleareanceStep"	    


	    # Use ortogonal subscales => OssSwitch
	    puts $fileid "${trailing_spaces}"
	    puts $fileid "${trailing_spaces}#other solver settings"
	    set cxpath "$rootid//c.SolutionStrategy//c.Advanced//i.OssSwitch"
	    set OssSwitch [::xmlutils::setXml $cxpath $cproperty]
	    puts $fileid "${trailing_spaces}oss_switch = $OssSwitch"
	    
	    # Calculate reactions
	    set cxpath "$rootid//c.Results//c.OnNodes//i.Reactions"
	    set Reactions [::xmlutils::setXml $cxpath $cproperty]
	    if {$Reactions =="Yes"} {
		    puts $fileid "${trailing_spaces}compute_reactions = True"
	    } else {
		    puts $fileid "${trailing_spaces}compute_reactions = False"
	    }
	    
	    
	    # Time order
	    set cxpath "$rootid//c.SolutionStrategy//c.Advanced//i.TimeOrder"
	    set TimeOrder [::xmlutils::setXml $cxpath $cproperty]
	    puts $fileid "${trailing_spaces}time_order = $TimeOrder"
	    
	    # Predictor corrector
	    set cxpath "$rootid//c.SolutionStrategy//c.Advanced//i.PredictorCorrector"
	    set PredictorCorrector [::xmlutils::setXml $cxpath $cproperty]
	    puts $fileid "${trailing_spaces}predictor_corrector = $PredictorCorrector"

	    # Use dt in stabilization => DynamicTau
	    set cxpath "$rootid//c.SolutionStrategy//c.Advanced//i.DynamicTau"
	    set DynamicTau [::xmlutils::setXml $cxpath $cproperty]
	    puts $fileid "${trailing_spaces}dynamic_tau = $DynamicTau"
	    
	    if {$SolverType in [list "ElementBased" "EdgeBased"]} {
		# Maximum velocity iterations
		set cxpath "$rootid//c.SolutionStrategy//c.Advanced//i.MaximumVelocityIterations"
		set MaximumVelocityIterations [::xmlutils::setXml $cxpath $cproperty]
		puts $fileid "${trailing_spaces}max_vel_its = $MaximumVelocityIterations"
		
		# Maximum pressure iterations
		set cxpath "$rootid//c.SolutionStrategy//c.Advanced//i.MaximumPressureIterations"
		set MaximumPressureIterations [::xmlutils::setXml $cxpath $cproperty]
		puts $fileid "${trailing_spaces}max_press_its = $MaximumPressureIterations"
		
	    } elseif {$SolverType in [list "Monolithic" "PressureSplitting"]} {
		# Maximum iterations
		set cxpath "$rootid//c.SolutionStrategy//c.Advanced//i.MaximumIterations"
		set MaximumIterations [::xmlutils::setXml $cxpath $cproperty]
		puts $fileid "${trailing_spaces}max_iteration = $MaximumIterations"
		
	    }
	    
	    # Laplacian form
	    set cxpath "$rootid//c.SolutionStrategy//c.Advanced//i.LaplacianForm"
	    set LaplacianForm [::xmlutils::setXml $cxpath $cproperty]
	    # Get the kratos keyword
	    set ckword [::xmlutils::getKKWord $kxpath $LaplacianForm]
	    if {$SolverType eq "ElementBased"} {
		# Set the default value
		puts $fileid "${trailing_spaces}laplacian_form = 1"
	    } else {
		puts $fileid "${trailing_spaces}laplacian_form = $ckword"
	    }
	    
	} else {
	    # Solver type for free surface
	    set cxpath "$rootid//c.AnalysisData//i.SolverTypeFreeSurf"
	    set SolverTypeFreeSurf [::xmlutils::setXml $cxpath $cproperty]
	    # WarnWinText "SolverTypeFreeSurf:$SolverTypeFreeSurf"
	    
	    # Check for use OpenMP
	    # Kratos key word xpath
	    set kxpath "Applications/$rootid"
	    set cxpath "$rootid//c.SolutionStrategy//c.ParallelType//i.ParallelSolutionType"
	    set ParallelSolutionType [::xmlutils::setXml $cxpath $cproperty]
	    
	    if {$ParallelSolutionType eq "OpenMP"} {
		# Write some project parameters used in the level set solver
		set pidlist [list "RedistanceFrequency" "ExtrapolationLayers" "StabdtPressureFactor" "StabdtConvectionFactor" "WallLawY" "SafetyFactor" "NumberOfInitialSteps"]
		foreach cvar $pidlist {
		    set cxpath "$rootid//c.SolutionStrategy//c.Advanced//i.[list ${cvar}]"
		    set cvarvalue [::xmlutils::setXml $cxpath $cproperty]
		    set ckword [::xmlutils::getKKWord $kxpath $cvar]
		    puts $fileid "$ckword = $cvarvalue"
		}
		
		# Write the body force properties
		set contid "LevelSetBodyForce"
		set cxpath "$rootid//c.SolutionStrategy//c.[list ${contid}]//i.GravityValue"
		set GravityValue [::xmlutils::setXml $cxpath $cproperty]
		set cxpath "$rootid//c.SolutionStrategy//c.[list ${contid}]//i.Cx"
		set Cx [::xmlutils::setXml $cxpath $cproperty]
		set cxpath "$rootid//c.SolutionStrategy//c.[list ${contid}]//i.Cy"
		set Cy [::xmlutils::setXml $cxpath $cproperty]
		set cxpath "$rootid//c.SolutionStrategy//c.[list ${contid}]//i.Cz"
		set Cz [::xmlutils::setXml $cxpath $cproperty]
		# wa "GravityValue:$GravityValue Cx:$Cx Cy:$Cy Cz:$Cz"
		set kwordlist [list "body_force_x" "body_force_y" "body_force_z"]
		set valuelist [list [expr $GravityValue*$Cx] [expr $GravityValue*$Cy] [expr $GravityValue*$Cz]]
		foreach ckword $kwordlist cvalue $valuelist {
		    puts $fileid "$ckword = $cvalue"
		}
		
		# Porous zone properties
		set contid "PorousZones"
		set cxpath "$rootid//c.SolutionStrategy//c.[list ${contid}]//i.UseErgunEquation"
		set UseErgunEquation [::xmlutils::setXml $cxpath $cproperty]
		if {$UseErgunEquation eq "Yes"} {
		    puts $fileid "UseErgun = True "
		} else {
		    puts $fileid "UseErgun = False "
		}
		
		# Material properties (Density and viscosity)
		set Density [::wkcf::GetFluidMaterialProperties $AppId "PropertyId" "Density"]
		set Viscosity [::wkcf::GetFluidMaterialProperties $AppId "PropertyId" "Viscosity"]
		# wa "Density:$Density Viscosity:$Viscosity"
		puts $fileid "Density = $Density "
		puts $fileid "Viscosity = $Viscosity "
	    }
	}
    }
    
    puts $fileid ""
    # Start time
    set cxpath "$rootid//c.ProblemParameters//i.StartTime"
    set StartTime [::xmlutils::setXml $cxpath $cproperty]
    # End time
    set cxpath "$rootid//c.ProblemParameters//i.EndTime"
    set EndTime [::xmlutils::setXml $cxpath $cproperty]
    # Delta time
    set cxpath "$rootid//c.ProblemParameters//i.DeltaTime"
    set DeltaTime [::xmlutils::setXml $cxpath $cproperty]
    
    # For use automatic delta time
    puts $fileid "#general problem settings"
    set cxpath "$rootid//c.ProblemParameters//i.UseAutomaticDeltaTime"
    set UseAutomaticDeltaTime [::xmlutils::setXml $cxpath $cproperty]
    puts $fileid "AutomaticDeltaTime = \"$UseAutomaticDeltaTime\""
    

    
    # WarnWinText "StartTime:$StartTime EndTime:$EndTime DeltaTime:$DeltaTime"
    puts $fileid "Dt = $DeltaTime"
    puts $fileid "Start_time = $StartTime"
    puts $fileid "max_time = $EndTime"
    set nsteps [expr int(double($EndTime-$StartTime)/double($DeltaTime))]
    puts $fileid "nsteps = $nsteps"
    
    puts $fileid ""

    

    
    # Check for use slip conditions
    if {[info exists dprops($AppId,UseSlipConditions)]} {
	if {$dprops($AppId,UseSlipConditions)} {
	    puts $fileid "Use_slip_conditions = True"
	}
    }
    
    # Write the group dictionary
    set arrinfo [array get dprops $AppId,Mesh,*,MeshIdGroup]
    if {[llength $arrinfo]} {
	puts $fileid ""
	puts $fileid "groups_dictionary = \{"
	foreach {name val} $arrinfo {
	    set lastchar [string last "h," $name] 
	    set firstchar [string first ",MeshIdGroup" $name] 
	    set groupid "[string range $name [expr $lastchar+2] [expr $firstchar-1]]"
	    # wa "name:$name val:$val lastchar:$lastchar firstchar:$firstchar groupid:$groupid"
	    if {$val !=""} {
		puts $fileid "        \"$groupid\" : $val," 
	    }
	}
	puts $fileid "                   \}"
    }
    

    puts $fileid "#output settings"
    # Output step 
    set cxpath "$rootid//c.Results//i.OutputDeltaTime"
    set OutputDeltaTime [::xmlutils::setXml $cxpath $cproperty]
    puts $fileid "output_time = $OutputDeltaTime"
    # WarnWinText "OutputDeltaTime:$OutputDeltaTime"
    set output_step [expr int($OutputDeltaTime/double($DeltaTime))]
    # WarnWinText "output_step:$output_step"
    puts $fileid "output_step = $output_step"
    
    # For results
    if {$ndime =="3D"} {
	#  For volumen output
	set cxpath "$rootid//c.Results//i.VolumeOutput"
	set VolumeOutput [::xmlutils::setXml $cxpath $cproperty]
	if {$VolumeOutput eq "Yes"} {
	    puts $fileid "VolumeOutput = True"
	} else {
	    puts $fileid "VolumeOutput = False"
	}
    } else {
	# Set the default value for 2D
	puts $fileid "VolumeOutput = True"
    }
    
    # On nodes results
    set cnrlist [list "Velocity" "Pressure" "Reactions" "Distance"]
    # set cnrlist [list "Velocity" "Pressure" "Reactions"]
    set nodal_results "nodal_results=\["
    foreach cnr $cnrlist {
	set cxpath "$rootid//c.Results//c.OnNodes//i.[list ${cnr}]"
	set cproperty "dv"
	set cvalue [::xmlutils::setXml $cxpath $cproperty]
	if {$cvalue =="Yes"} {
	    set cnkr [::xmlutils::getKKWord $kxpath $cnr]
	    append nodal_results "\"$cnkr\","
	}
    }
    set findcomma [string last "," $nodal_results]
    if {$findcomma !="-1"} {
	set nodal_results [string range $nodal_results 0 end-1]
	append nodal_results "\]" 
	puts $fileid "$nodal_results"
    }
    
    # Set gauss_points_results to empty
    puts $fileid "gauss_points_results=\[\]"
    
    # WarnWinText "nodal_results:$nodal_results"
    
    # GiD post mode variables
    ::wkcf::WriteGiDPostMode $AppId $fileid 
    
    puts $fileid ""
    set PName [::KUtils::GetPaths "PName"]
    puts $fileid "problem_name=\"${PName}${AppId}\"" 
    #puts $fileid "problem_path=\"$PDir\"" 
    
    # Get the kratos path 
    set cxpath "GeneralApplicationData//c.ProjectConfiguration//i.KratosPath"
    set cproperty "dv"
    set KratosPath [::xmlutils::setXml $cxpath $cproperty]
    set KratosPath [file native $KratosPath]
    
    # Write the kratos path
    puts $fileid "kratos_path=\"${KratosPath}\"" 
}

proc ::wkcf::WritePFEMLagrangianFluidProjectParameters {AppId fileid PDir} {
    variable ndime; variable dprops
    
    set rootid "$AppId"
    
    # Kratos key word xpath
    set kxpath "Applications/$rootid"

    set domain_size 2
    if {$ndime =="2D"} {
	set domain_size 2
    } elseif {$ndime =="3D"} {
	set domain_size 3
    }
   
    # Domain size
    puts $fileid "domain_size = $domain_size"
    puts $fileid ""

    # Get others properties
    set cproperty "dv"
    
    # Fluid type
    set cxpath "$rootid//c.AnalysisData//i.FluidType"
    set FluidType [::xmlutils::setXml $cxpath $cproperty]
    
    # WarnWinText "FluidType:$FluidType"
    if {$FluidType =="Compressible"} {

    } elseif {$FluidType =="Incompressible"} {

	# FSI
	puts $fileid "FSI = 0.00000e+00"        
	
	# Compute Reactions
	puts $fileid "compute_reactions = 0.00000e+00"
	
	# Dt
	set cxpath "$rootid//c.ProblemParameters//i.DeltaTime"
	set DeltaTime [::xmlutils::setXml $cxpath $cproperty]
	puts $fileid "Dt = $DeltaTime"
	
	# max_time
	set cxpath "$rootid//c.ProblemParameters//i.EndTime"
	set EndTime [::xmlutils::setXml $cxpath $cproperty]
	puts $fileid "max_time = $EndTime"
	
	# Output Step
	set cxpath "$rootid//c.Results//i.OutputDeltaTime"
	set OutStep [::xmlutils::setXml $cxpath $cproperty]
	puts $fileid "output_step = $OutStep"
	
	# Alpha Shape
	puts $fileid "alpha_shape = 1.60000e+00"
	
	# Erase Nodes
	puts $fileid "erase_nodes = 1.00000e+00"
	
	# Adaptive Refinement
	puts $fileid "adaptive_refinement = 1.00000e+00"
	
	# Delete Nodes Close to Wall
	puts $fileid "delete_nodes_close_to_wall = 1.00000e+00"

	# Material properties (Density, viscosity and Bulk Modulus)
	set Density [::wkcf::GetFluidMaterialProperties $AppId "PropertyId" "Density"]
	set Viscosity [::wkcf::GetFluidMaterialProperties $AppId "PropertyId" "Viscosity"]
	set BulkModulus [::wkcf::GetFluidMaterialProperties $AppId "PropertyId" "BulkModulus"]
	# wa "Density:$Density Viscosity:$Viscosity BulkModulus:$BulkModulus"
	puts $fileid "density = $Density "
	puts $fileid "viscosity = $Viscosity "
	puts $fileid "bulk_modulus = $BulkModulus "
	
	# Gravity
	set cxpath "$rootid//c.ProblemParameters//i.PFEMBodyForceGravity"
	set Gravity [::xmlutils::setXml $cxpath $cproperty]
	puts $fileid "with_gravity = $Gravity"
	
	# Bounding Box
	set cxpath "$rootid//c.SolutionStrategy//c.Advanced//c.Boundingbox//i.MinX"
	set B1X [::xmlutils::setXml $cxpath $cproperty]
	puts $fileid "bounding_box_corner1_x = $B1X"
	set cxpath "$rootid//c.SolutionStrategy//c.Advanced//c.Boundingbox//i.MinY"
	set B1Y [::xmlutils::setXml $cxpath $cproperty]
	puts $fileid "bounding_box_corner1_y = $B1X"
	set cxpath "$rootid//c.SolutionStrategy//c.Advanced//c.Boundingbox//i.MinZ"
	set B1Z [::xmlutils::setXml $cxpath $cproperty]
	puts $fileid "bounding_box_corner1_z = $B1X"
	set cxpath "$rootid//c.SolutionStrategy//c.Advanced//c.Boundingbox//i.MaxX"
	set B2X [::xmlutils::setXml $cxpath $cproperty]
	puts $fileid "bounding_box_corner2_x = $B2X"
	set cxpath "$rootid//c.SolutionStrategy//c.Advanced//c.Boundingbox//i.MaxY"
	set B2Y [::xmlutils::setXml $cxpath $cproperty]
	puts $fileid "bounding_box_corner2_y = $B2Y"
	set cxpath "$rootid//c.SolutionStrategy//c.Advanced//c.Boundingbox//i.MaxZ"
	set B2Z [::xmlutils::setXml $cxpath $cproperty]
	puts $fileid "bounding_box_corner2_z = $B2Z"
	
	# Lagrangian Nodes Inlet
	puts $fileid "lagrangian_nodes_inlet = 0.00000e+00"   

	# Incompressible Modified FracStep
	puts $fileid "SolverType = \"Incompressible_Modified_FracStep\""   
		
	puts $fileid ""
	
	# Comment
	puts $fileid "# Declare Python Variables"   
	
	puts $fileid ""
	set PName [::KUtils::GetPaths "PName"]
	puts $fileid "problem_name=\"${PName}${AppId}\"" 
	puts $fileid "problem_path=\"$PDir\"" 

	# Get the kratos path 
	set cxpath "GeneralApplicationData//c.ProjectConfiguration//i.KratosPath"
	set cproperty "dv"
	set KratosPath [::xmlutils::setXml $cxpath $cproperty]
	set KratosPath [file native $KratosPath]

	# Write the kratos path
	puts $fileid "kratos_path=\"${KratosPath}\"" 
    }
}

proc ::wkcf::WriteLinearSolvers {rootid fileid vartype wfsmethod trailing_spaces config_name} {
    # Write linear solver for all applications
    
    # Kratos key word xpath
    set kxpath "Applications/$rootid"
    # Set default value xml variable
    set cproperty "dv"


    puts $fileid "    # $vartype solver"
    # Define the class
    puts $fileid "    class ${config_name}:"


    set cxpath "$rootid//c.SolutionStrategy//c.SolverTypes//i.[list ${vartype}]LinearSolverType"
    set LinearSolverType [::xmlutils::setXml $cxpath $cproperty]
    if {$LinearSolverType =="Direct"} {
	# Direct solver type
	set cxpath "$rootid//c.SolutionStrategy//c.SolverTypes//i.[list ${vartype}]DirectSolverType"
	set DirectSolverType [::xmlutils::setXml $cxpath $cproperty]
	# WarnWinText "DirectSolverType:$DirectSolverType"
	set cDirectSolverType [::xmlutils::getKKWord $kxpath $DirectSolverType]
	
	    puts $fileid "${trailing_spaces}    solver_type = \"$cDirectSolverType\""
	    puts $fileid "${trailing_spaces}    scaling = False"

	
    } elseif {$LinearSolverType =="Iterative"} {
	
	# Iterative solver type 
	set cxpath "$rootid//c.SolutionStrategy//c.SolverTypes//i.[list ${vartype}]IterativeSolverType"
	set IterativeSolverType [::xmlutils::setXml $cxpath $cproperty]
	# Tolerance
	set cxpath "$rootid//c.SolutionStrategy//c.SolverTypes//i.[list ${vartype}]ISTolerance"
	set Tolerance [::xmlutils::setXml $cxpath $cproperty]
	# Maximum iteration
	set cxpath "$rootid//c.SolutionStrategy//c.SolverTypes//i.[list ${vartype}]ISMaximumIteration"
	set MaximumIteration [::xmlutils::setXml $cxpath $cproperty]
	# preconditioner type
	set cxpath "$rootid//c.SolutionStrategy//c.SolverTypes//i.[list ${vartype}]PreconditionerType"
	set PreconditionerType [::xmlutils::setXml $cxpath $cproperty]
	# WarnWinText "vartype:$vartype IterativeSolverType:$IterativeSolverType Tolerance:$Tolerance MaximumIteration:$MaximumIteration PreconditionerType:$PreconditionerType"
	
	# Solver type
	set lsolver [::xmlutils::getKKWord $kxpath $IterativeSolverType]
	# Preconditioner
	set precond [::xmlutils::getKKWord $kxpath $PreconditionerType]
	
	    puts $fileid "${trailing_spaces}    solver_type = \"$lsolver\""
	    puts $fileid "${trailing_spaces}    tolerance = $Tolerance"
	    puts $fileid "${trailing_spaces}    max_iteration = $MaximumIteration"
	    puts $fileid "${trailing_spaces}    preconditioner = \"$precond\""
	    puts $fileid "${trailing_spaces}    scaling = False"

    }
}

proc ::wkcf::WriteCutAndGraph {AppId} {
    
    # ABSTRACT: Write the cutting and point history properties
    variable dprops;  variable ActiveAppList
    variable ndime

    set filename "define_output.py"
    set PDir [::KUtils::GetPaths "PDir"]
    
    set fullname [file native [file join $PDir $filename]]
    
    # First delete the file
    if {[file exists $fullname]} {
	set res [file delete -force $fullname]
    }

    if { [catch { set fileid [open $fullname w+] }] } {
	WarnWin [= "Cannot write file %s. Permission denied" $fullname].
	return 0
    }

    # Write the output point properties
    # Kratos key word xpath
    set kxpath "Applications/$AppId"

    puts $fileid ""
    puts $fileid "def DefineOutputPoints():"
    puts $fileid "    output_points = \[\]"
    # Get the value
    set basexpath "$AppId//c.Results//c.HistoryOutputOnPoints"
    set opproplist [::xmlutils::setXmlContainerIds $basexpath]
    foreach copid $opproplist {
	# Get the output point properties
	set cxpath "${basexpath}//c.[list ${copid}]//c.HistoryGraph"
	set alloprop [::xmlutils::setXmlContainerPairs $cxpath "" "dv"]
	# wa "alloprop:$alloprop"
	if {[llength $alloprop]} {
	    set cpprop "\["
	    foreach cprop [lrange $alloprop 0 2] {
		lassign $cprop dummy val
		append cpprop "$val,"
	    }
	    set cpprop [string replace $cpprop end end "\]"]
	    
	    # Get the variable identifier
	    if {$ndime eq "2D"} {
		set varid [lindex $alloprop 3 1]
	    } elseif {$ndime eq "3D"} {
		set varid [lindex $alloprop 4 1]
	    } 
	    # Get the variable identifier from the kratos keyword
	    set varkword [::xmlutils::getKKWord $kxpath "$varid" "kkword"]
	    # wa "varkword:$varkword"
	    # Get the file identifier
	    set filename [lindex $alloprop 5 1]
	}
	# add a point with coordinates x,y,z to be written to "fileidentifier"
	puts $fileid "    output_points.append( \[ $cpprop, \"${varkword}\", \"$filename\" \] )"
    }
    # Return the output point properties
    puts $fileid "    return output_points "

    # Write the cut plane properties
    puts $fileid ""
    puts $fileid "def DefineCutPlanes():"
    puts $fileid "    cut_planes_list = \[\]"
    
    # adding a plane with normal (nx,ny,nz) passing through the point with coordinates (x,y,z) and 
    # named "plane_name"
    
    # Get the value
    set basexpath "$AppId//c.Results//c.CutOptions"
    set cutproplist [::xmlutils::setXmlContainerIds $basexpath]
    # WarnWinText "cutproplist:$cutproplist"
    foreach ccutid $cutproplist {
	# WarnWinText "ccutid:$ccutid"
	
	# Get the origin cut properties
	set cxpath "${basexpath}//c.[list ${ccutid}]//c.Origin"
	set alloprop [::xmlutils::setXmlContainerPairs $cxpath "" "dv"]
	if {[llength $alloprop]} {
	    set cpprop "\["
	    foreach cprop $alloprop {
		lassign $cprop dummy val
		append cpprop "$val,"
	    }
	    set cpprop [string replace $cpprop end end "\]"]
	}
	
	# Get the normal to plane properties
	set cxpath "${basexpath}//c.[list ${ccutid}]//c.NormalToPlane"
	set allnprop [::xmlutils::setXmlContainerPairs $cxpath "" "dv"]
	if {[llength $allnprop]} {
	    set cnprop "\["
	    foreach cprop $allnprop {
		lassign $cprop dummy val
		append cnprop "$val,"
	    }
	    set cnprop [string replace $cnprop end end "\]"]
	}
	# Append this plane
	puts $fileid "    cut_planes_list.append( \[ $cnprop, $cpprop, \"$ccutid\" \] )"
    }
    
    # Return the cut properties
    puts $fileid "    return cut_planes_list"
    
    # Write the drag forcest properties
    # Kratos key word xpath
    set kxpath "Applications/$AppId"

    puts $fileid ""
    puts $fileid "def DefineDragList():"
    puts $fileid "    drag_list = \[\]"
    # Get the value
    set basexpath "$AppId//c.Results//c.DragOptions"
    set dragproplist [::xmlutils::setXmlContainerIds $basexpath]
    # wa "dragproplist:$dragproplist"
    foreach cdpid $dragproplist {
	# Get the drag force properties
	set cxpath "${basexpath}//c.[list ${cdpid}]//c.MainProperties"
	set alldragprop [::xmlutils::setXmlContainerPairs $cxpath "" "dv"]
	# wa "alldragprop:$alldragprop"
	if {[llength $alldragprop]} {
	    if {(([info exists dprops($AppId,Mesh,$cdpid,MeshIdGroup)]) && ($dprops($AppId,Mesh,$cdpid,MeshIdGroup) !=""))} { 
		# Get the mesh identifier
		set MeshIdGroup $dprops($AppId,Mesh,$cdpid,MeshIdGroup)
		    
		# Get the file identifier
		set filename [lindex $alldragprop 0 1]
		
		# add each mesh identifier with the file when will be write the drag forces "fileidentifier"
		puts $fileid "    drag_list.append( \[ $MeshIdGroup, \"$filename\" \] )"
	    }
	}
    }

    # Return the drag forces properties
    puts $fileid "    return drag_list"

    close $fileid
}

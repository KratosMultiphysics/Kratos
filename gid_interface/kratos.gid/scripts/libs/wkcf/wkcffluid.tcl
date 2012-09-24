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
    # Write some properties at the nodal level for Fluid application
    variable wmethod

    # For debug
    if {!$::wkcf::pflag} {
	set inittime [clock seconds]
    }
    switch -exact -- $wmethod {
	"0" {
	    ::wkcf::WritePropertyAtNodes_m0 $AppId
	}
	"1" {
	    ::wkcf::WritePropertyAtNodes_m1 $AppId
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

proc ::wkcf::WritePropertyAtNodes_m1 {AppId} {
    # Write some properties at the nodal level for Fluid application
    variable dprops

    set cproplist [list "Density" "Viscosity"]
  
    # Check for all defined kratos elements
    if {([info exists dprops($AppId,AllKElemId)]) && ([llength $dprops($AppId,AllKElemId)])} {

	
	# Write viscosity and density for each node identifier
	set Density 0.0; set Viscosity 0.0
	foreach PropertyId $dprops($AppId,GKProps,AllPropertyId) {
	    # Get the material identifier for this property 
	    set MatId $dprops($AppId,Property,$PropertyId,MatId) 
	    # Get the group identifier
	    set GroupId $dprops($AppId,Property,$PropertyId,GroupId)
	    # Get all material properties
	    set mpxpath "[::KMat::findMaterialParent $MatId]//m.${MatId}"
	    # WarnWinText "mpxpath:$mpxpath"
	    # Get the material properties
	    foreach pid $cproplist {
		if {$pid =="Density"} {
		    set xpath "c.General"
		    # Get the current value for this properties
		    set cvalue [lindex [::KMat::getMaterialProperties "p" "$mpxpath//$xpath//p.$pid"] 0 1]
		    set Density [GiD_FormatReal "%10.5e" $cvalue]
		} elseif {$pid =="Viscosity"} {
		    set xpath "c.Fluid"
		    # Get the current value for this properties
		    set cvalue [lindex [::KMat::getMaterialProperties "p" "$mpxpath//$xpath//p.$pid"] 0 1]
		    set Viscosity [GiD_FormatReal "%10.5e" $cvalue]
		}
	    }
	    # Only the first property
	    break 
	}
	# WarnWinText "Density:$Density Viscosity:$Viscosity"
	
	set kxpath "Materials"
	set cpropid "0"        

	# Write the group nodal properties
	foreach celemid $dprops($AppId,AllKElemId) {
	    # Check for all defined group identifier for this element
	    if {([info exists dprops($AppId,KElem,$celemid,AllGroupId)]) && ([llength $dprops($AppId,KElem,$celemid,AllGroupId)])} {
		# For all defined group identifier for this element
		foreach cgroupid $dprops($AppId,KElem,$celemid,AllGroupId) {
		    # Group properties format for density and viscosity
		    set gprop_visco [dict create]
		    set gprop_densi [dict create]		    
		    # Write viscosity value for this group
		    set f "%10i [format "%4i" $cpropid]   $Viscosity\n"
		    set f [subst $f]
		    dict set gprop_visco $cgroupid "$f"
		    if {[write_calc_data nodes -count $gprop_visco]>0} {
			set vkword [::xmlutils::getKKWord $kxpath "Viscosity" "kkword"]
			write_calc_data puts "Begin NodalData $vkword \/\/ GUI group identifier: $cgroupid"
			write_calc_data nodes -sorted $gprop_visco
			write_calc_data puts "End NodalData"
		        write_calc_data puts ""
		    }

		    # Write density value for this group
		    set f "%10i [format "%4i" $cpropid]   $Density\n"
		    set f [subst $f]
		    dict set gprop_densi $cgroupid "$f"
		    if {[write_calc_data nodes -count $gprop_visco]>0} {
			set vkword [::xmlutils::getKKWord $kxpath "Density" "kkword"]
			write_calc_data puts "Begin NodalData $vkword \/\/ GUI group identifier: $cgroupid"
			write_calc_data nodes -sorted $gprop_densi
			write_calc_data puts "End NodalData"
		        write_calc_data puts ""
		    }
		    # Unset the group dictionary
		    unset gprop_visco
		    unset gprop_densi
		}
	    }
	}
    }
}

proc ::wkcf::WritePropertyAtNodes_m0 {AppId} {
    # Write some properties at the nodal level for Fluid application
    variable dprops

    set cproplist [list "Density" "Viscosity"]
  
    # Check for all defined kratos elements
    if {([info exists dprops($AppId,AllKElemId)]) && ([llength $dprops($AppId,AllKElemId)])} {
	
	# Create a dictionary
	set nc [dict create 0 0]
	# For all defined kratos elements        
	foreach celemid $dprops($AppId,AllKElemId) {
	    # Check for all defined group identifier for this element
	    if {([info exists dprops($AppId,KElem,$celemid,AllGroupId)]) && ([llength $dprops($AppId,KElem,$celemid,AllGroupId)])} {
		# For all defined group identifier for this element
		foreach cgroupid $dprops($AppId,KElem,$celemid,AllGroupId) {
		    # Get the GiD entity type, element type and property identifier
		    lassign $dprops($AppId,KElem,$celemid,$cgroupid,GProps) GiDEntity GiDElemType PropertyId KEKWord nDim
		    # WarnWinText "GiDEntity:$GiDEntity GiDElemType:$GiDElemType PropertyId:$PropertyId KEKWord:$KEKWord nDim:$nDim"
		    # Get all defined entities for this group identifier
		    set allelist [::KUtils::GetDefinedMeshGiDEntities $cgroupid $GiDEntity]
		    # WarnWinText "alllist:$allelist"
		    if {[llength $allelist]} {
		        foreach elemid $allelist {
			    # Get the element properties
			    foreach nodeid [lrange [GiD_Info Mesh Elements $GiDElemType $elemid] 1 end-1] {
		                dict set nc $nodeid $cgroupid 
		            }
			}
		    }
		}
	    }
	}

	# Write viscosity and density for each node identifier
	set Density 0.0; set Viscosity 0.0
	foreach PropertyId $dprops($AppId,GKProps,AllPropertyId) {
	    # Get the material identifier for this property 
	    set MatId $dprops($AppId,Property,$PropertyId,MatId) 
	    # Get the group identifier
	    set GroupId $dprops($AppId,Property,$PropertyId,GroupId)
	    # Get all material properties
	    set mpxpath "[::KMat::findMaterialParent $MatId]//m.${MatId}"
	    # WarnWinText "mpxpath:$mpxpath"
	    # Get the material properties
	    foreach pid $cproplist {
		if {$pid =="Density"} {
		    set xpath "c.General"
		    # Get the current value for this properties
		    set cvalue [lindex [::KMat::getMaterialProperties "p" "$mpxpath//$xpath//p.$pid"] 0 1]
		    set Density [GiD_FormatReal "%10.5e" $cvalue]
		} elseif {$pid =="Viscosity"} {
		    set xpath "c.Fluid"
		    # Get the current value for this properties
		    set cvalue [lindex [::KMat::getMaterialProperties "p" "$mpxpath//$xpath//p.$pid"] 0 1]
		    set Viscosity [GiD_FormatReal "%10.5e" $cvalue]
		}
	    }
	    # Only the first property
	    break 
	}
	# WarnWinText "Density:$Density Viscosity:$Viscosity"
	
	set kxpath "Materials"
	set cpropid "0"        

	# Write the group nodal properties
	foreach celemid $dprops($AppId,AllKElemId) {
	    # Check for all defined group identifier for this element
	    if {([info exists dprops($AppId,KElem,$celemid,AllGroupId)]) && ([llength $dprops($AppId,KElem,$celemid,AllGroupId)])} {
		# For all defined group identifier for this element
		foreach cgroupid $dprops($AppId,KElem,$celemid,AllGroupId) {
		    set cnodeglist [list]
		    dict for {nodeid dgroupid} $nc {
		        if {$cgroupid ==$dgroupid} {
		            lappend cnodeglist $nodeid
		        }
		    }
		    if {[llength $cnodeglist]} {
		        # Write all nodes for this group in incresing orden
		        set viscobf ""; set densibf ""
		        foreach nodeid [lsort -integer $cnodeglist] {
		            append viscobf "[format "%4i  %4i" $nodeid $cpropid]    $Viscosity\n"
		            append densibf "[format "%4i  %4i" $nodeid $cpropid]    $Density\n"
		        }
		        # Write viscosity value for this group
		        set vkword [::xmlutils::getKKWord $kxpath "Viscosity" "kkword"]
		        write_calc_data puts "Begin NodalData $vkword \/\/ GUI group identifier: $cgroupid"
		        write_calc_data puts "[string trimright $viscobf]"
		        write_calc_data puts "End NodalData"
		        write_calc_data puts ""
		        
		        # Write density value for this group 
		        set dkword [::xmlutils::getKKWord $kxpath "Density" "kkword"]
		        write_calc_data puts "Begin NodalData $dkword \/\/ GUI group identifier: $cgroupid"
		        write_calc_data puts "[string trimright $densibf]"
		        write_calc_data puts "End NodalData"
		        write_calc_data puts ""
		    }
		}
		unset cnodeglist
	    }
	}
	unset nc
    }
}

proc ::wkcf::WriteFluidBC {AppId inletvelglist noslipglist flagvariablelist kwordlist} {
    # ABSTRACT: Write the fluid boundary conditions
    variable wmethod

    # WarnWinText "inletvelglist:$inletvelglist\nnoslipglist:$noslipglist\nflagvariablelist:$flagvariablelist\nkwordlist:$kwordlist"
   
    if {([llength $inletvelglist]) || ([llength $noslipglist])} {
	# For debug
	if {!$::wkcf::pflag} {
	    set inittime [clock seconds]
	}
	switch -exact -- $wmethod {
	    "0" {
		::wkcf::WriteFluidInletNoSlipBC_m0 $AppId $inletvelglist $noslipglist $kwordlist
	    }
	    "1" {
		::wkcf::WriteFluidInletNoSlipBC_m1 $AppId $inletvelglist $noslipglist $kwordlist
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

proc ::wkcf::WriteFluidInletNoSlipBC_m1 {AppId inletvelglist noslipglist kwordlist} {
    variable ndime; variable dprops

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

	    # Group properties format 
	    set gprop_xcomp [dict create]
	    set gprop_ycomp [dict create]
	    set gprop_zcomp [dict create]
    
	    # X component
	    if {$nsx} {
		set f "%10i [format "%8i" $cpropid]   [GiD_FormatReal "%10.5e" $nsxval]\n"
		set f [subst $f]
		dict set gprop_xcomp $nsgroupid "$f"
		if {[write_calc_data nodes -count $gprop_xcomp]>0} {
		    write_calc_data puts "Begin NodalData $xitem \/\/ No-slip condition GUI group identifier: $nsgroupid"
		    write_calc_data nodes -sorted $gprop_xcomp
		    write_calc_data puts "End NodalData"
		    write_calc_data puts ""
		}
	    }
	    
	    # Y component
	    if {$nsy} {
		set f "%10i [format "%8i" $cpropid]   [GiD_FormatReal "%10.5e" $nsyval]\n"
		set f [subst $f]
		dict set gprop_ycomp $nsgroupid "$f"
		if {[write_calc_data nodes -count $gprop_ycomp]>0} {
		    write_calc_data puts "Begin NodalData $yitem \/\/ No-slip condition GUI group identifier: $nsgroupid"
		    write_calc_data nodes -sorted $gprop_ycomp
		    write_calc_data puts "End NodalData"
		    write_calc_data puts ""
		}
	    }
	    
	    # Z component
	    if {$ndime =="3D"} {
		if {$nsz} {
		    set f "%10i [format "%8i" $cpropid]   [GiD_FormatReal "%10.5e" $nszval]\n"
		    set f [subst $f]
		    dict set gprop_zcomp $nsgroupid "$f"
		    if {[write_calc_data nodes -count $gprop_zcomp]>0} {
			write_calc_data puts "Begin NodalData $zitem \/\/ No-slip condition GUI group identifier: $nsgroupid"
			write_calc_data nodes -sorted $gprop_zcomp
			write_calc_data puts "End NodalData"
			write_calc_data puts ""
		    }
		}
	    }
	    # Unset the dictionaries
	    unset gprop_xcomp 
	    unset gprop_ycomp 
	    unset gprop_zcomp
	}
    }
    
    # Use first the inlet
    if {([llength $inletvelglist]) && ([llength $noslipglist])} {
	set cformat "%10d"
	# Check to match node identifier
	set condmatch [dict create]
	# For each group in the no-slip condition
	foreach nsgroupid $noslipglist {
	    lassign $dprops($AppId,BC,$nscondid,$nsgroupid,GProps) cx cxval cy cyval cz czval
	    # WarnWinText "nsgroupid:$nsgroupid cx:$cx cxval:$cxval cy:$cy cyval:$cyval cz:$cz czval:$czval"
	    set gprop [dict create]
	    dict set gprop $nsgroupid "$cformat"
	    if {[write_calc_data nodes -count $gprop]>0} {
		# For each node in the no-slip bc update the condmatch 
		foreach nsnodeid [write_calc_data nodes -return $gprop] {
		    dict set condmatch $nsnodeid [list $cx $cy $cz]
		}
	    }
	    # Unset the dictionary
	    unset gprop 
	}
	
	# For all inlet velocity group identifier
	set ixcomp ""; set iycomp ""; set izcomp ""
	foreach igroupid $inletvelglist {
	    lassign $dprops($AppId,BC,$icondid,$igroupid,GProps) ix ixval iy iyval iz izval
	    # WarnWinText "igroupid:$igroupid ix:$ix iy:$iy iz:$iz"
	    # Set the inlet format dictionary
	    set gprop [dict create]
	    dict set gprop $igroupid "$cformat"
	    if {[write_calc_data nodes -count $gprop]>0} {
		# 3D problems
		if {$ndime =="3D"} {
		    # For each node in the inlet bc ckeck to write this node
		    foreach inodeid [write_calc_data nodes -return $gprop] {
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
			write_calc_data puts "Begin NodalData $xitem \/\/ Inlet velocity condition GUI group identifier: $igroupid"
			write_calc_data puts "[string trimright ${ixcomp}]"
			write_calc_data puts "End NodalData"
			write_calc_data puts ""
		    }
		    if {[string length $iycomp]} {
			write_calc_data puts "Begin NodalData $yitem \/\/ Inlet velocity condition GUI group identifier: $igroupid"
			write_calc_data puts "[string trimright ${iycomp}]"
			write_calc_data puts "End NodalData"
			write_calc_data puts ""
		    }
		    if {[string length $izcomp]} {
			write_calc_data puts "Begin NodalData $zitem \/\/ Inlet velocity condition GUI group identifier: $igroupid"
			write_calc_data puts "[string trimright ${izcomp}]"
			write_calc_data puts "End NodalData"
			write_calc_data puts ""
		    }
		    
		} elseif {$ndime =="2D"} {
		    # For each node in the inlet bc ckeck to write this node
		    foreach inodeid [write_calc_data nodes -return $gprop] {
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
			}
		    }
		    # Write this group identifier
		    if {[string length $ixcomp]} {
			write_calc_data puts "Begin NodalData $xitem \/\/ Inlet velocity condition GUI group identifier: $igroupid"
			write_calc_data puts "[string trimright ${ixcomp}]"
			write_calc_data puts "End NodalData"
			write_calc_data puts ""
		    }
		    if {[string length $iycomp]} {
			write_calc_data puts "Begin NodalData $yitem \/\/ Inlet velocity condition GUI group identifier: $igroupid"
			write_calc_data puts "[string trimright ${iycomp}]"
			write_calc_data puts "End NodalData"
			write_calc_data puts ""
		    }
		}
	    }
	    # Reset ixcomp, iycomp and zcomp
	    set ixcomp ""; set iycomp ""; set izcomp ""

	    # Unset the inlet format dictionary
	    unset gprop
	}
	
	# unset dictionary variable
	unset condmatch

    } else {
	# Write the inlet boundary condition properties

	# For each group in the inlet condition
	foreach igroupid $inletvelglist {
	    lassign $dprops($AppId,BC,$icondid,$igroupid,GProps) ix ixval iy iyval iz izval

	    # Group properties format 
	    set gprop_xcomp [dict create]
	    set gprop_ycomp [dict create]
	    set gprop_zcomp [dict create]   

	    # X component
	    if {$ix} {
		set f "%10i [format "%8i" $cpropid]   [GiD_FormatReal "%10.5e" $ixval]\n"
		set f [subst $f]
		dict set gprop_xcomp $igroupid "$f"
		if {[write_calc_data nodes -count $gprop_xcomp]>0} {
		    write_calc_data puts "Begin NodalData $xitem \/\/ Inlet velocity condition GUI group identifier: $igroupid"
		    write_calc_data nodes -sorted $gprop_xcomp
		    write_calc_data puts "End NodalData"
		    write_calc_data puts ""
		}
	    }
	    
	    # Y component
	    if {$iy} {
		set f "%10i [format "%8i" $cpropid]   [GiD_FormatReal "%10.5e" $iyval]\n"
		set f [subst $f]
		dict set gprop_ycomp $igroupid "$f"
		if {[write_calc_data nodes -count $gprop_ycomp]>0} {
		    write_calc_data puts "Begin NodalData $yitem \/\/ Inlet velocity condition GUI group identifier: $igroupid"
		    write_calc_data nodes -sorted $gprop_ycomp
		    write_calc_data puts "End NodalData"
		    write_calc_data puts ""
		}
	    }
	    
	    # Z component
	    if {$ndime =="3D"} {
		if {$iz} {
		    set f "%10i [format "%8i" $cpropid]   [GiD_FormatReal "%10.5e" $izval]\n"
		    set f [subst $f]
		    dict set gprop_zcomp $igroupid "$f"
		    if {[write_calc_data nodes -count $gprop_zcomp]>0} {
			write_calc_data puts "Begin NodalData $zitem \/\/ Inlet velocity condition GUI group identifier: $igroupid"
			write_calc_data nodes -sorted $gprop_zcomp
			write_calc_data puts "End NodalData"
			write_calc_data puts ""
		    }
		}
	    }
	    # Unset the dictionaries
	    unset gprop_xcomp 
	    unset gprop_ycomp 
	    unset gprop_zcomp
	}
    }
}

proc ::wkcf::WriteFluidInletNoSlipBC_m0 {AppId inletvelglist noslipglist kwordlist} {
    variable gidentitylist; variable ndime
    variable useqelem; variable dprops

    # Map Inlet-NoSlip => Use no-slip values at share nodes
    set icondid "InletVelocity"; set nscondid "No-Slip"
    set cpropid "1"
    set xitem [lindex $kwordlist 0]
    set yitem [lindex $kwordlist 1]
    set zitem [lindex $kwordlist 2]
    
    if {[llength $noslipglist]} {
	# No-slip
	foreach cgroupid $noslipglist {
	    set allnslip($cgroupid,NodeList) [list]
	    foreach GiDEntity $gidentitylist {
		# Get all defined entities for this group identifier
		switch $GiDEntity {
		    "point" {
		        set callnlist [::KUtils::GetDefinedMeshGiDEntities $cgroupid $GiDEntity "Nodes" $useqelem]
		        if {[llength $callnlist]} {
		            lappend allnslip($cgroupid,NodeList) {*}$callnlist
		        }
		    }
		    "line" - "surface" {
		        set callnlist [::KUtils::GetDefinedMeshGiDEntities $cgroupid $GiDEntity "Nodes" $useqelem]
		        if {[llength $callnlist]} {
		            lappend allnslip($cgroupid,NodeList) {*}$callnlist
		        }
		    }
		    "volume" {
		        if {$ndime =="3D"} {
		            set callnlist [::KUtils::GetDefinedMeshGiDEntities $cgroupid $GiDEntity "Nodes" $useqelem]
		            if {[llength $callnlist]} {
		                lappend allnslip($cgroupid,NodeList) {*}$callnlist
		            }
		        }
		    }
		} 
	    }
	    # WarnWinText "groupid:$cgroupid NodeList:$allnslip($cgroupid,NodeList)"
	}

	# Check to match node identifier
	# Use first the inlet
	set xcomp ""; set ycomp ""; set zcomp ""
	if {[llength $inletvelglist]} {
	    set condmatch [dict create none 0]
	    # For each group in the no-slip condition
	    foreach nsgroupid $noslipglist {
		set nsGProps $dprops($AppId,BC,$nscondid,$nsgroupid,GProps)
		# WarnWinText "nsgroupid:$nsgroupid nsGProps:$nsGProps"
		foreach nsnodeid $allnslip($nsgroupid,NodeList) {
		    set clist [list 0 0 0]
		    if {[lindex $nsGProps 0]} {
		        append xcomp "[format "%8i%8i" $nsnodeid $cpropid]   [GiD_FormatReal "%10.5e" [lindex $nsGProps 1]]\n"
		        lset clist 0 1 
		    }
		    if {[lindex $nsGProps 2]} {
		        append ycomp "[format "%8i%8i" $nsnodeid $cpropid]   [GiD_FormatReal "%10.5e" [lindex $nsGProps 3]]\n"
		        lset clist 1 1
		    }
		    if {[lindex $nsGProps 4]} {
		        append zcomp "[format "%8i%8i" $nsnodeid $cpropid]   [GiD_FormatReal "%10.5e" [lindex $nsGProps 5]]\n"
		        lset clist 2 1
		    }
		    dict set condmatch $nsnodeid $clist 
		}
		
		# Write this group identifier
		if {[string length $xcomp]} {
		    write_calc_data puts "Begin NodalData $xitem \/\/ No-slip condition GUI group identifier: $nsgroupid"
		    write_calc_data puts "[string trimright ${xcomp}]"
		    write_calc_data puts "End NodalData"
		    write_calc_data puts ""
		}
		if {[string length $ycomp]} {
		    write_calc_data puts "Begin NodalData $yitem \/\/ No-slip condition GUI group identifier: $nsgroupid"
		    write_calc_data puts "[string trimright ${ycomp}]"
		    write_calc_data puts "End NodalData"
		    write_calc_data puts ""
		}
		if {$ndime =="3D"} {
		    if {[string length $zcomp]} {
		        write_calc_data puts "Begin NodalData $zitem \/\/ No-slip condition GUI group identifier: $nsgroupid"
		        write_calc_data puts "[string trimright ${zcomp}]"
		        write_calc_data puts "End NodalData"
		        write_calc_data puts ""
		    }
		}

		# Reset xcomp, ycomp and zcomp
		set xcomp ""; set ycomp ""; set zcomp ""
	    }

	    # Get the Inlet velocity entities
	    foreach cgroupid $inletvelglist {
		set allninlet($cgroupid,NodeList) [list]
		foreach GiDEntity $gidentitylist {
		    # Get all defined entities for this group identifier
		    switch $GiDEntity {
		        "point" {
		            set callnlist [::KUtils::GetDefinedMeshGiDEntities $cgroupid $GiDEntity "Nodes" $useqelem]
		            if {[llength $callnlist]} {
		                lappend allninlet($cgroupid,NodeList) {*}$callnlist
		            }
		        }
		        "line" - "surface" {
		            set callnlist [::KUtils::GetDefinedMeshGiDEntities $cgroupid $GiDEntity "Nodes" $useqelem]
		            if {[llength $callnlist]} {
		                lappend allninlet($cgroupid,NodeList) {*}$callnlist
		            }
		        }
		        "volume" {
		            if {$ndime =="3D"} {
		                set callnlist [::KUtils::GetDefinedMeshGiDEntities $cgroupid $GiDEntity "Nodes" $useqelem]
		                if {[llength $callnlist]} {
		                    lappend allninlet($cgroupid,NodeList) {*}$callnlist
		                }
		            }
		        }
		    } 
		}
		# WarnWinText "groupid:$cgroupid NodeList:$allninlet($cgroupid,NodeList)"
	    }
	    
	    # For all inlet velocity group identifier
	    set ixcomp ""; set iycomp ""; set izcomp ""
	    foreach igroupid $inletvelglist {
		set iGProps $dprops($AppId,BC,$icondid,$igroupid,GProps)
		# WarnWinText "igroupid:$igroupid iGProps:$iGProps"
		foreach inodeid $allninlet($igroupid,NodeList) {
		    # WarnWinText "inodeid:$inodeid"
		    # Check that this node identifier exists in the dictionary
		    if {[dict exists $condmatch $inodeid]} {
		        # Get the properties
		        set nprop [dict get $condmatch $inodeid] 
		        # WarnWinText "nprop:$nprop"
		        # Check x flag
		        if {[lindex $nprop 0]=="0"} {
		            # Write this node identifier
		            if {[lindex $iGProps 0]} {
		                append ixcomp "[format "%8i%8i" $inodeid $cpropid]   [GiD_FormatReal "%10.5e" [lindex $iGProps 1]]\n"
		            }
		        }
		        # Check y flag
		        if {[lindex $nprop 1]=="0"} {
		            if {[lindex $iGProps 2]} {
		                append iycomp "[format "%8i%8i" $inodeid $cpropid]   [GiD_FormatReal "%10.5e" [lindex $iGProps 3]]\n"
		            }
		        }
		        # Check z flag
		        if {[lindex $nprop 2]=="0"} {
		            if {[lindex $iGProps 4]} {
		                append izcomp "[format "%8i%8i" $inodeid $cpropid]   [GiD_FormatReal "%10.5e" [lindex $iGProps 5]]\n"
		            }
		        }
		    } else {
		        # Write this node identifier
		        if {[lindex $iGProps 0]} {
		            append ixcomp "[format "%8i%8i" $inodeid $cpropid]   [GiD_FormatReal "%10.5e" [lindex $iGProps 1]]\n"
		        }
		        if {[lindex $iGProps 2]} {
		            append iycomp "[format "%8i%8i" $inodeid $cpropid]   [GiD_FormatReal "%10.5e" [lindex $iGProps 3]]\n"
		        }
		        if {[lindex $iGProps 4]} {
		            append izcomp "[format "%8i%8i" $inodeid $cpropid]   [GiD_FormatReal "%10.5e" [lindex $iGProps 5]]\n"
		        }
		    }
		}
		
		# Write this group identifier
		if {[string length $ixcomp]} {
		    write_calc_data puts "Begin NodalData $xitem \/\/ Inlet velocity condition GUI group identifier: $igroupid"
		    write_calc_data puts "[string trimright ${ixcomp}]"
		    write_calc_data puts "End NodalData"
		    write_calc_data puts ""
		}
		if {[string length $iycomp]} {
		    write_calc_data puts "Begin NodalData $yitem \/\/ Inlet velocity condition GUI group identifier: $igroupid"
		    write_calc_data puts "[string trimright ${iycomp}]"
		    write_calc_data puts "End NodalData"
		    write_calc_data puts ""
		}
		if {$ndime =="3D"} {
		    if {[string length $izcomp]} {
		        write_calc_data puts "Begin NodalData $zitem \/\/ Inlet velocity condition GUI group identifier: $igroupid"
		        write_calc_data puts "[string trimright ${izcomp}]"
		        write_calc_data puts "End NodalData"
		        write_calc_data puts ""
		    }
		}
		
		# Reset ixcomp, iycomp and zcomp
		set ixcomp ""; set iycomp ""; set izcomp ""
	    }
	    
	    # unset dictionary variable
	    unset condmatch
	    if {[info exists allninlet]} {
		unset allninlet
	    }
	    if {[info exists allnslip]} {
		unset allnslip
	    }
	    
	} else {

	    # Write all no-slip condition
	    # For each group in the no-slip condition
	    foreach nsgroupid $noslipglist {
		set nsGProps $dprops($AppId,BC,$nscondid,$nsgroupid,GProps)
		foreach nsnodeid $allnslip($nsgroupid,NodeList) {
		    if {[lindex $nsGProps 0]} {
		        append xcomp "[format "%8i%8i" $nsnodeid $cpropid]   [GiD_FormatReal "%10.5e" [lindex $nsGProps 1]]\n"
		    }
		    if {[lindex $nsGProps 2]} {
		        append ycomp "[format "%8i%8i" $nsnodeid $cpropid]   [GiD_FormatReal "%10.5e" [lindex $nsGProps 3]]\n"
		    }
		    if {[lindex $nsGProps 4]} {
		        append zcomp "[format "%8i%8i" $nsnodeid $cpropid]   [GiD_FormatReal "%10.5e" [lindex $nsGProps 5]]\n"
		    }
		}
		
		# Write this group identifier
		if {[string length $xcomp]} {
		    write_calc_data puts "Begin NodalData $xitem \/\/ No-slip condition GUI group identifier: $nsgroupid"
		    write_calc_data puts "[string trimright ${xcomp}]"
		    write_calc_data puts "End NodalData"
		    write_calc_data puts ""
		}
		if {[string length $ycomp]} {
		    write_calc_data puts "Begin NodalData $yitem \/\/ No-slip condition GUI group identifier: $nsgroupid"
		    write_calc_data puts "[string trimright ${ycomp}]"
		    write_calc_data puts "End NodalData"
		    write_calc_data puts ""
		}
		if {$ndime =="3D"} {
		    if {[string length $zcomp]} {
		        write_calc_data puts "Begin NodalData $zitem \/\/ No-slip condition GUI group identifier: $nsgroupid"
		        write_calc_data puts "[string trimright ${zcomp}]"
		        write_calc_data puts "End NodalData"
		        write_calc_data puts ""
		    }
		}
		
		# Reset xcomp, ycomp and zcomp
		set xcomp ""; set ycomp ""; set zcomp ""
	    }
	    if {[info exists allnslip]} {
		unset allnslip
	    }
	}
    }
}

proc ::wkcf::WriteFluidFlagVariableBC {AppId flagvariablelist} {
    # ABSTRACT: Write the flag variable boundary condition
    variable wmethod
    
    switch -exact -- $wmethod {
	"0" {
	    ::wkcf::WriteFluidFlagVariableBC_m0 $AppId $flagvariablelist
	}
	"1" {
	    ::wkcf::WriteFluidFlagVariableBC_m1 $AppId $flagvariablelist
	}
    }
}

proc ::wkcf::WriteFluidFlagVariableBC_m1 {AppId flagvariablelist} {
    variable dprops

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
	# Group dictionary properties format 
	set gprop_fv [dict create]
	set gprop_isb [dict create]

	if {$flagval=="2"} {
	    # Write this group identifier
	    set flag1 1

	    # Flag-Variable
	    set f "%10i [format "%8i%8i" $cpropid $flagval]\n"
	    set f [subst $f]
	    dict set gprop_fv $cgroupid "$f"
	    if {[write_calc_data nodes -count $gprop_fv]>0} {
		write_calc_data puts "Begin NodalData $fvitem \/\/ Flag-Variable condition GUI group identifier: $cgroupid"
		write_calc_data nodes -sorted $gprop_fv
		write_calc_data puts "End NodalData"
		write_calc_data puts ""
	
		# Update the flagvar2 dictionary
		foreach nodeid [write_calc_data nodes -return $gprop_fv] {
		    dict set flagvar2 $nodeid $cgroupid
		}
	    }

	    # is_boundary
	    set f "%10i [format "%8i%8i" $cpropid $isbcpropid]\n"
	    set f [subst $f]
	    dict set gprop_isb $cgroupid "$f"
	    if {[write_calc_data nodes -count $gprop_isb]>0} {
		write_calc_data puts "Begin NodalData $isbitem \/\/ is_boundary associated with Flag-Variable condition GUI group identifier: $cgroupid"
		write_calc_data nodes -sorted $gprop_isb
		write_calc_data puts "End NodalData"
		write_calc_data puts ""
	    }
	}
	# unset dictionaries used for format groups
	unset gprop_fv
	unset gprop_isb
    }

    # Write all group with flag-variable equal to 1
    if {$flag1} {
	set fvcomp ""; set isbcomp ""
	# For each group in the flag-variable condition
	foreach cgroupid $flagvariablelist {
	    lassign $dprops($AppId,BC,$flagvarcondid,$cgroupid,GProps) flagval
	    # WarnWinText "cgroupid:$cgroupid GProps:$GProps"
	    set gprop_fv [dict create]
	    if {$flagval=="1"} {
		set f "%10i [format "%8i%8i" $cpropid $flagval]\n"
		set f [subst $f]
		dict set gprop_fv $cgroupid "$f"
		if {[write_calc_data nodes -count $gprop_fv]>0} {
		    foreach nodeid [write_calc_data nodes -return $gprop_fv] {
			if {![dict exists $flagvar2 $nodeid]} {
			    append fvcomp "[format "%10i%8i%8i" $nodeid $cpropid $flagval]\n"
			    append isbcomp "[format "%10i%8i%8i" $nodeid $cpropid $isbcpropid]\n"
			}
		    }
		    # Write this group identifier
		    # Flag-Variable
		    if {[string length $fvcomp]} {
			write_calc_data puts "Begin NodalData $fvitem \/\/ Flag-Variable condition GUI group identifier: $cgroupid"
			write_calc_data puts "[string trimright ${fvcomp}]"
			write_calc_data puts "End NodalData"
			write_calc_data puts ""
		    }
		    # is_boundary
		    if {[string length $isbcomp]} {
			write_calc_data puts "Begin NodalData $isbitem \/\/ is_boundary associated with Flag-Variable condition GUI group identifier: $cgroupid"
			write_calc_data puts "[string trimright ${isbcomp}]"
			write_calc_data puts "End NodalData"
			write_calc_data puts ""
		    }
		    # Reset components
		    set fvcomp ""; set isbcomp ""
		}
	    }
	    # unset dictionary used for format groups
	    unset gprop_fv
	}
    }
   
    # unset the dict for all nodes with flag equal to 2
    unset flagvar2
}

proc ::wkcf::WriteFluidFlagVariableBC_m0 {AppId flagvariablelist} {
    variable gidentitylist; variable ndime
    variable useqelem; variable dprops

    #WarnWinText "flagvariablelist:$flagvariablelist"
    # For nodes with many flag variable defined flag of level two have the priority over flag of level one
    set flagvarcondid "Flag-Variable"
    set cpropid "0"
    set isbcpropid "1"

    
    # Flag-Variable
    foreach cgroupid $flagvariablelist {
	set allnflagvar($cgroupid,NodeList) [list]
	foreach GiDEntity $gidentitylist {
	    # Get all defined entities for this group identifier
	    switch $GiDEntity {
		"point" {
		    set callnlist [::KUtils::GetDefinedMeshGiDEntities $cgroupid $GiDEntity "Nodes" $useqelem]
		    if {[llength $callnlist]} {
			lappend allnflagvar($cgroupid,NodeList) {*}$callnlist
		    }
		}
		"line" - "surface" {
		    set callnlist [::KUtils::GetDefinedMeshGiDEntities $cgroupid $GiDEntity "Nodes" $useqelem]
		    if {[llength $callnlist]} {
			lappend allnflagvar($cgroupid,NodeList) {*}$callnlist
		    }
		}
		"volume" {
		    if {$ndime =="3D"} {
			set callnlist [::KUtils::GetDefinedMeshGiDEntities $cgroupid $GiDEntity "Nodes" $useqelem]
			if {[llength $callnlist]} {
			    lappend allnflagvar($cgroupid,NodeList) {*}$callnlist
			}
		    }
		}
	    } 
	}
	# WarnWinText "groupid:$cgroupid NodeList:$allnflagvar($cgroupid,NodeList)"
    }
    
    # Write the flag condition
    set flag1 0
    set fvcomp ""; set isbcomp ""
    set fvitem "FLAG-VARIABLE"
    set isbitem "IS_BOUNDARY"
    set flagvar2 [dict create]
    # For each group in the flag-variable condition
    # Create a dict for all nodes with flag equal to 2
    foreach groupid $flagvariablelist {
	set GProps $dprops($AppId,BC,$flagvarcondid,$groupid,GProps)
	# WarnWinText "groupid:$groupid GProps:$GProps"
	if {[lindex $GProps 0]=="2"} {
	    foreach nodeid $allnflagvar($groupid,NodeList) {
		append fvcomp "[format "%8i%8i%8i" $nodeid $cpropid [lindex $GProps 0]]\n"
		append isbcomp "[format "%8i%8i%8i" $nodeid $cpropid $isbcpropid]\n"
		dict set flagvar2 $nodeid $groupid
		set flag1 1
	    }
	    # Write this group identifier
	    # Flag-Variable
	    if {[string length $fvcomp]} {
		write_calc_data puts "Begin NodalData $fvitem \/\/ Flag-Variable condition GUI group identifier: $groupid"
		write_calc_data puts "[string trimright ${fvcomp}]"
		write_calc_data puts "End NodalData"
		write_calc_data puts ""
	    }
	    # is_boundary
	    if {[string length $isbcomp]} {
		write_calc_data puts "Begin NodalData $isbitem \/\/ is_boundary associated with Flag-Variable condition GUI group identifier: $groupid"
		write_calc_data puts "[string trimright ${isbcomp}]"
		write_calc_data puts "End NodalData"
		write_calc_data puts ""
	    }
	    # Reset components
	    set fvcomp "";        set isbcomp ""
	}
    }
    
    # Write all group with flag-variable equal to 1
    if {$flag1} {
	set fvcomp ""; set isbcomp ""
	# For each group in the flag-variable condition
	foreach groupid $flagvariablelist {
	    set GProps $dprops($AppId,BC,$flagvarcondid,$groupid,GProps)
	    # WarnWinText "groupid:$groupid GProps:$GProps"
	    if {[lindex $GProps 0]=="1"} {
		foreach nodeid $allnflagvar($groupid,NodeList) {
		    if {![dict exists $flagvar2 $nodeid]} {
			append fvcomp "[format "%8i%8i%8i" $nodeid $cpropid [lindex $GProps 0]]\n"
			append isbcomp "[format "%8i%8i%8i" $nodeid $cpropid $isbcpropid]\n"
		    }
		}
		# Write this group identifier
		# Flag-Variable
		if {[string length $fvcomp]} {
		    write_calc_data puts "Begin NodalData $fvitem \/\/ Flag-Variable condition GUI group identifier: $groupid"
		    write_calc_data puts "[string trimright ${fvcomp}]"
		    write_calc_data puts "End NodalData"
		    write_calc_data puts ""
		}
		# is_boundary
		if {[string length $isbcomp]} {
		    write_calc_data puts "Begin NodalData $isbitem \/\/ is_boundary associated with Flag-Variable condition GUI group identifier: $groupid"
		    write_calc_data puts "[string trimright ${isbcomp}]"
		    write_calc_data puts "End NodalData"
		    write_calc_data puts ""
		}
		# Reset components
		set fvcomp ""; set isbcomp ""
	    }
	}
    }

    # Unset temporal variables
    unset flagvar2
    if {[info exists allnflagvar]} {
	unset allnflagvar
    }
}

proc ::wkcf::WriteFluidIsSlipBC {AppId ccondid kwordlist} {
    # ASTRACT: Write is-slip boundary conditions => Conditional data
    variable dprops; variable wmethod 
    variable ndime; variable ctbclink
    
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

	# wa "activateval:$activateval"
	if {$wmethod} {
	    set gprop [dict create]
	    set f ""
	    if {$ndime == "2D"} {
		set GiDElemType "Linear"
		set f "%10d %10d %10d"
		dict set gprop $cgroupid "$f"
		if {[write_calc_data has_elements -elemtype $GiDElemType $gprop]} {
		    set dprops($AppId,UseSlipConditions) 1
		    set f "%10d [format "%4i" $activateval] %10d %10d\n"
		    set f [subst $f]
		    dict set gprop $cgroupid "$f"
		    write_calc_data puts "Begin ConditionalData $isstructurekw // GUI is-slip condition group identifier: $cgroupid"
		    # write_calc_data connectivities -sorted $gprop
		    foreach {elemid cfixval nodei nodej} [write_calc_data connectivities -return -elemtype "$GiDElemType" $gprop] {
			# wa "elemid:$elemid cfixval:$cfixval nodei:$nodei nodej:$nodej"
			# Check that exists this element in the dictionary with the condition indentifier links
			if {[dict exists $ctbclink $elemid]} {
			    set condid [dict get $ctbclink $elemid]
			    write_calc_data puts "[format "%10d %10d" $condid $activateval]"
			}
		    }
		    write_calc_data puts "End ConditionalData"
		    write_calc_data puts ""
		    unset gprop

		    # Write Y_Wall values
		    set gprop [dict create]
		    set f "%10d [format "%4d" $state] [format "%10g" $ConstantValue]\n"
		    set f [subst $f]
		    dict set gprop $cgroupid "$f"
		    write_calc_data puts " Begin NodalData $isywallkw // GUI Y-Wall condition group identifier: $cgroupid"
		    write_calc_data nodes $gprop 
		    write_calc_data puts "End NodalData"
		    write_calc_data puts ""
		
		    unset gprop
		}
	    } elseif {$ndime == "3D"} {
		set f "%10d %10d %10d %10d"
		dict set gprop $cgroupid "$f"
		set GiDElemType "Triangle"
		if {[write_calc_data has_elements -elemtype $GiDElemType $gprop]} {
		    set dprops($AppId,UseSlipConditions) 1
		    set f "%10d [format "%4i" $activateval] %10d %10d %10d\n"
		    set f [subst $f]
		    dict set gprop $cgroupid "$f"
		    write_calc_data puts "Begin ConditionalData $isstructurekw // GUI is-slip condition group identifier: $cgroupid"
		    # write_calc_data connectivities -sorted $gprop
		    foreach {elemid cfixval nodei nodej nodek} [write_calc_data connectivities -return -elemtype "$GiDElemType" $gprop] {
			# Check that exists this element in the dictionary with the condition indentifier links
			if {[dict exists $ctbclink $elemid]} {
			    set condid [dict get $ctbclink $elemid]
			    write_calc_data puts "[format "%10d %10d" $condid $cfixval]"
			}
		    }
		    write_calc_data puts "End ConditionalData"
		    write_calc_data puts ""
		    unset gprop
		    
		    # Write Y_Wall values
		    set gprop [dict create]
		    set f "%10d [format "%4d" $state] [format "%10g" $ConstantValue]\n"
		    set f [subst $f]
		    dict set gprop $cgroupid "$f"
		    write_calc_data puts " Begin NodalData $isywallkw // GUI Y-Wall condition group identifier: $cgroupid"
		    write_calc_data nodes $gprop 
		    write_calc_data puts "End NodalData"
		    write_calc_data puts ""
		
		    unset gprop
		} 
	    }
	}
    }

    # For debug
    if {!$::wkcf::pflag} {
	set endtime [clock seconds]
	set ttime [expr $endtime-$inittime]
	# WarnWinText "endtime:$endtime ttime:$ttime"
	WarnWinText "Write fluid wall is-slip boundary conditions: [::KUtils::Duration $ttime]"
    }
}

proc ::wkcf::WriteFluidWallLawBC {AppId ccondid kwordlist} {
    # ASTRACT: Write wall law boundary conditions => Nodal data
    variable dprops; variable wmethod 

    # For debug
    if {!$::wkcf::pflag} {
	set inittime [clock seconds]
    }

    # For all defined group identifier inside this condition type
    foreach cgroupid $dprops($AppId,BC,$ccondid,AllGroupId) {
	# wa "cgroupid:$cgroupid"
	# Get the condition properties
	lassign $dprops($AppId,BC,$ccondid,$cgroupid,GProps) cvalue
	# wa "cvalue:$cvalue"
	if {$wmethod} {
	    set gprop [dict create]
	    set f "%10i"
	    dict set gprop $cgroupid "$f"
	    if {[write_calc_data nodes -count $gprop]>0} {
		set f "%10i [format "%10.5f" $cvalue]\n"
		set f [subst $f]
		dict set gprop $cgroupid "$f"
		write_calc_data puts "Begin NodalData $kwordlist // GUI wall law condition group identifier: $cgroupid"
		write_calc_data nodes -sorted $gprop
		write_calc_data puts "End NodalData"
		write_calc_data puts ""
	    }
	    unset gprop
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

proc ::wkcf::WriteOutLetPressureBC {AppId ccondid kwordlist} {
    # ASBTRACT: Write outlet pressure boundary condition
    variable wmethod

    # For debug
    if {!$::wkcf::pflag} {
	set inittime [clock seconds]
    }
    switch -exact -- $wmethod {
	"0" {
	    ::wkcf::WriteOutLetPressureBC_m0 $AppId $ccondid $kwordlist
	}
	"1" {
	    ::wkcf::WriteOutLetPressureBC_m1 $AppId $ccondid $kwordlist
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

proc ::wkcf::WriteOutLetPressureBC_m1 {AppId ccondid kwordlist} {
    # Write outlet pressure boundary condition
    variable dprops

    set kitem [lindex $kwordlist 0]
    
    # For all defined group identifier inside this condition type
    foreach cgroupid $dprops($AppId,BC,$ccondid,AllGroupId) {
	# Get the condition properties
	lassign $dprops($AppId,BC,$ccondid,$cgroupid,GProps) fixval pval
	# WarnWinText "fixval:$fixval pval:$pval"
	set gprop [dict create]
	# Fix x
	if {$fixval =="1"} {
	    set f "%10i"
	    dict set gprop $cgroupid "$f"
	    if {[write_calc_data nodes -count $gprop]>0} {
		set f "%10i [format "%4i%10.5f" $fixval $pval]\n"
		set f [subst $f]
		dict set gprop $cgroupid "$f"
		# Write the pressure value
		write_calc_data puts "Begin NodalData $kitem"
		write_calc_data nodes -sorted $gprop
		write_calc_data puts "End NodalData"
		write_calc_data puts ""
	    }
	}
	unset gprop
    }
}

proc ::wkcf::WriteOutLetPressureBC_m0 {AppId ccondid kwordlist} {
    # Write outlet pressure boundary condition
    variable ndime;    variable gidentitylist
    variable useqelem; variable dprops

    set nodelist [list]
    
    # For all defined group identifier inside this condition type
    foreach cgroupid $dprops($AppId,BC,$ccondid,AllGroupId) {
	# Get the condition properties
	set GProps $dprops($AppId,BC,$ccondid,$cgroupid,GProps)
	# WarnWinText "GProps:$GProps"
	# Assign values
	lassign $GProps fixval pressureval
	# WarnWinText "fixval:$fixval pressureval:$pressureval"
	set allnlist [list]
	foreach GiDEntity $gidentitylist {
	    # Get all defined entities for this group identifier
	    switch $GiDEntity {
		"point" {
		    set callnlist [::KUtils::GetDefinedMeshGiDEntities $cgroupid $GiDEntity "Nodes" $useqelem]
		    if {[llength $callnlist]} {
			lappend allnlist $callnlist
		    }
		}
		"line" - "surface" {
		    set callnlist [::KUtils::GetDefinedMeshGiDEntities $cgroupid $GiDEntity "Nodes" $useqelem]
		    if {[llength $callnlist]} {
			lappend allnlist $callnlist
		    }
		} 
	    }
	}
	# WarnWinText "$GiDEntity alllist:$allnlist"
	foreach cprop $allnlist {
	    set cprop [lsort -integer -unique $cprop]
	    foreach nodeid $cprop {
		# Fix x
		if {$fixval =="1"} {
		    lappend nodelist "$nodeid 1 $pressureval"
		}
	    }
	}
	unset allnlist
    }
    # WarnWinText "nodelist:$nodelist"
    
    # PRESSURE
    if {[llength $nodelist]} {
	set kitem [lindex $kwordlist 0]
	write_calc_data puts "Begin NodalData $kitem"
	foreach citem $nodelist {
	    lassign $citem nodeid fix pval
	    set cf "[format "%4i%4i%10.5f" $nodeid $fix $pval]"
	    write_calc_data puts "$cf"
	}
	write_calc_data puts "End NodalData"
	write_calc_data puts ""
	unset nodelist
    }
}

proc ::wkcf::WriteFluidProjectParameters {AppId fileid PDir} {
    variable ndime; variable dprops
    
    # Write fluid solver method
    # 0 => Old format
    # 1 => New format
    set wfsmethod 0
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

	# Free surface
	set cxpath "$rootid//c.AnalysisData//i.FreeSurface"
	set FreeSurface [::xmlutils::setXml $cxpath $cproperty]
	# WarnWinText "FreeSurface:$FreeSurface"
	if {$FreeSurface =="No"} {
	    # Solver type
	    set cxpath "$rootid//c.AnalysisData//i.SolverType"
	    set SolverType [::xmlutils::setXml $cxpath $cproperty]
	    # WarnWinText "SolverType:$SolverType"
	    # Get the kratos keyword
	    set ckword [::xmlutils::getKKWord $kxpath $SolverType]
	    # WarnWinText "ckword:$ckword"
	    if {$wfsmethod=="0"} {
		puts $fileid "SolverType = \"$ckword\""
	    }

	    # Monolithic,PressureSplitting,ElementBased,EdgeBased
	    puts $fileid ""
	    switch -exact -- $SolverType {
		"ElementBased" {
		    # Fractional step options => ElementBased                
		    # Solution strategy
		    # Linear solvers
		    # Velocity
		    ::wkcf::WriteFluidSolvers $rootid $fileid "Velocity"
		    puts $fileid ""
		    # Pressure
		    ::wkcf::WriteFluidSolvers $rootid $fileid "Pressure"                 
		}
		"PressureSplitting" {
		    # Pressure splitting
		    # Solution strategy
		    # Linear solvers
		    # Velocity
		    ::wkcf::WriteFluidSolvers $rootid $fileid "Velocity"
		    puts $fileid ""
		    # Pressure
		    ::wkcf::WriteFluidSolvers $rootid $fileid "Pressure"
		}
		"Monolithic" {
		    # Monolithic

		    # Solution strategy
		    # Linear solvers
		    # Velocity
		    ::wkcf::WriteFluidSolvers $rootid $fileid "Monolithic"
		    puts $fileid ""
		}
	    }
	    
	    if {$wfsmethod} {
		# Fluid solver configuration
		puts $fileid ""
		puts $fileid "# Fluid solver configuration"
		puts $fileid "class FluidSolverConfiguration:"
		puts $fileid "    is_active = True"
		puts $fileid "    use_defaults = True"
		puts $fileid "    solving_strategy_type =  \"$ckword\""
		puts $fileid "    velocity_solver_parameters = FluidVelocityLinearSolverConfiguration()"
		puts $fileid "    pressure_solver_parameters = FluidPressureLinearSolverConfiguration()"
		
		# Get the turbulence properties
		set cxpath "$rootid//c.AnalysisData//i.TurbulenceModel"
		set TurbulenceModel [::xmlutils::setXml $cxpath $cproperty]
		# WarnWinText "TurbulenceModel:$TurbulenceModel"
		if {$TurbulenceModel eq "Off"} {
		    puts $fileid "    TurbulenceModel = \"None\""
		} elseif {$TurbulenceModel eq "Smagorinsky-Lilly"} {
		    puts $fileid "    TurbulenceModel = \"$TurbulenceModel\""
		    # Get the smagorinsky-lilly constant
		    set cxpath "$rootid//c.AnalysisData//i.SmagorinskyConstant"
		    set SmagorinskyConstant [::xmlutils::setXml $cxpath $cproperty]
		    # WarnWinText "SmagorinskyConstant:$SmagorinskyConstant"
		    puts $fileid "    SmagorinskyConstant = $SmagorinskyConstant"
		} elseif {$TurbulenceModel eq "Spalart-Allmaras"} {
		    puts $fileid "    TurbulenceModel = \"$TurbulenceModel\""
		}

	    } else {

		# Get the turbulence properties
		set cxpath "$rootid//c.AnalysisData//i.TurbulenceModel"
		set TurbulenceModel [::xmlutils::setXml $cxpath $cproperty]
		# WarnWinText "TurbulenceModel:$TurbulenceModel"
		if {$TurbulenceModel eq "Off"} {
		    puts $fileid "TurbulenceModel = \"None\""
		} elseif {$TurbulenceModel eq "Smagorinsky-Lilly"} {
		    puts $fileid "TurbulenceModel = \"$TurbulenceModel\""
		    # Get the smagorinsky-lilly constant
		    set cxpath "$rootid//c.AnalysisData//i.SmagorinskyConstant"
		    set SmagorinskyConstant [::xmlutils::setXml $cxpath $cproperty]
		    # WarnWinText "SmagorinskyConstant:$SmagorinskyConstant"
		    puts $fileid "SmagorinskyConstant = $SmagorinskyConstant"
		} elseif {$TurbulenceModel eq "Spalart-Allmaras"} {
		    puts $fileid "TurbulenceModel = \"$TurbulenceModel\""
		    # Get the value of the turbulence viscosity
		    set cxpath "$rootid//c.AnalysisData//i.TurbulentViscosity"
		    set TurbulentViscosity [::xmlutils::setXml $cxpath $cproperty]
		    # wa "TurbulentViscosity:$TurbulentViscosity"
		    puts $fileid "TurbulentViscosity = $TurbulentViscosity"
		    
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
			    set cxpath "${basexpath}//c.${cgroupid}//c.MainProperties"
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
			    puts $fileid "SA_wall_group_ids = $endmeshidlist"
			}
		    }
		}
	    }

	    # Write relative and absolute tolerances
	    puts $fileid ""
	    set ctlist [list "RelativeVelocityTolerance" "AbsoluteVelocityTolerance" "RelativePressureTolerance" "AbsolutePressureTolerance"]
	    foreach cv $ctlist {
		set cxpath "$rootid//c.SolutionStrategy//c.Advanced//i.${cv}"
		set cvalue [::xmlutils::setXml $cxpath $cproperty]
		set ckword [::xmlutils::getKKWord $kxpath $cv]
		puts $fileid "$ckword = $cvalue"
	    }
	    
	    puts $fileid ""  
	    # Time order
	    set cxpath "$rootid//c.SolutionStrategy//c.Advanced//i.TimeOrder"
	    set TimeOrder [::xmlutils::setXml $cxpath $cproperty]
	    puts $fileid "time_order = $TimeOrder"

	    # Predictor corrector
	    set cxpath "$rootid//c.SolutionStrategy//c.Advanced//i.PredictorCorrector"
	    set PredictorCorrector [::xmlutils::setXml $cxpath $cproperty]
	    puts $fileid "predictor_corrector = $PredictorCorrector"
	    
	    if {$SolverType in [list "ElementBased" "EdgeBased"]} {
		# Maximum velocity iterations
		set cxpath "$rootid//c.SolutionStrategy//c.Advanced//i.MaximumVelocityIterations"
		set MaximumVelocityIterations [::xmlutils::setXml $cxpath $cproperty]
		puts $fileid "max_vel_its = $MaximumVelocityIterations"
		
		# Maximum pressure iterations
		set cxpath "$rootid//c.SolutionStrategy//c.Advanced//i.MaximumPressureIterations"
		set MaximumPressureIterations [::xmlutils::setXml $cxpath $cproperty]
		puts $fileid "max_press_its = $MaximumPressureIterations"

	    } elseif {$SolverType in [list "Monolithic" "PressureSplitting"]} {
		# Maximum iterations
		set cxpath "$rootid//c.SolutionStrategy//c.Advanced//i.MaximumIterations"
		set MaximumIterations [::xmlutils::setXml $cxpath $cproperty]
		puts $fileid "max_iterations = $MaximumIterations"

	    }

	    # Laplacian form
	    set cxpath "$rootid//c.SolutionStrategy//c.Advanced//i.LaplacianForm"
	    set LaplacianForm [::xmlutils::setXml $cxpath $cproperty]
	    # Get the kratos keyword
	    set ckword [::xmlutils::getKKWord $kxpath $LaplacianForm]
	    if {$SolverType eq "ElementBased"} {
		# Set the default value
		puts $fileid "laplacian_form = 1"
	    } else {
		puts $fileid "laplacian_form = $ckword"
	    }
	    
	} else {
	    # Solver type for free surface
	    set cxpath "$rootid//c.AnalysisData//i.SolverTypeFreeSurf"
	    set SolverTypeFreeSurf [::xmlutils::setXml $cxpath $cproperty]
	    # WarnWinText "SolverTypeFreeSurf:$SolverTypeFreeSurf"
	}
    }

    puts $fileid ""
    # Start time
    set cxpath "$rootid//c.SolutionStrategy//i.StartTime"
    set StartTime [::xmlutils::setXml $cxpath $cproperty]
    # End time
    set cxpath "$rootid//c.SolutionStrategy//i.EndTime"
    set EndTime [::xmlutils::setXml $cxpath $cproperty]
    # Delta time
    set cxpath "$rootid//c.SolutionStrategy//i.DeltaTime"
    set DeltaTime [::xmlutils::setXml $cxpath $cproperty]
    
    # For use automatic delta time
    set cxpath "$rootid//c.SolutionStrategy//i.UseAutomaticDeltaTime"
    set UseAutomaticDeltaTime [::xmlutils::setXml $cxpath $cproperty]
    puts $fileid "AutomaticDeltaTime = \"$UseAutomaticDeltaTime\""

    # Get the divergence clearance step
    set cxpath "$rootid//c.SolutionStrategy//i.DivergenceCleareanceStep"
    set DivergenceCleareanceStep [::xmlutils::setXml $cxpath $cproperty]
    # Get the kratos keyword
    set DivergenceCleareanceStepKW [::xmlutils::getKKWord $kxpath "DivergenceCleareanceStep"]
    # wa "DivergenceCleareanceStepKW:$DivergenceCleareanceStepKW DivergenceCleareanceStep:$DivergenceCleareanceStep"
    puts $fileid "$DivergenceCleareanceStepKW = $DivergenceCleareanceStep"

    # WarnWinText "StartTime:$StartTime EndTime:$EndTime DeltaTime:$DeltaTime"
    puts $fileid "Dt = $DeltaTime"
    puts $fileid "Start_time = $StartTime"
    puts $fileid "max_time = $EndTime"
    set nsteps [expr int(double($EndTime-$StartTime)/double($DeltaTime))]
    puts $fileid "nsteps = $nsteps"

    puts $fileid ""
    # Use dt in stabilization => DynamicTau
    set cxpath "$rootid//c.SolutionStrategy//c.Advanced//i.DynamicTau"
    set DynamicTau [::xmlutils::setXml $cxpath $cproperty]
    puts $fileid "use_dt_in_stabilization = $DynamicTau"
    
    # Use ortogonal subscales => OssSwitch
    set cxpath "$rootid//c.SolutionStrategy//c.Advanced//i.OssSwitch"
    set OssSwitch [::xmlutils::setXml $cxpath $cproperty]
    puts $fileid "use_orthogonal_subscales = $OssSwitch"
    
    # Calculate reactions
    set cxpath "$rootid//c.Results//c.OnNodes//i.Reactions"
    set Reactions [::xmlutils::setXml $cxpath $cproperty]
    if {$Reactions =="Yes"} {
	puts $fileid "Calculate_reactions = True"
    } else {
	puts $fileid "Calculate_reactions = False"
    }

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

    puts $fileid ""
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
    }

    puts $fileid ""
    # On nodes results
    set cnrlist [list "Velocity" "Pressure" "Reactions" "Distance"]
    # set cnrlist [list "Velocity" "Pressure" "Reactions"]
    set nodal_results "nodal_results=\["
    foreach cnr $cnrlist {
	set cxpath "$rootid//c.Results//c.OnNodes//i.${cnr}"
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
    puts $fileid "problem_path=\"$PDir\"" 
    
    # Get the kratos path 
    set cxpath "GeneralApplicationData//c.ProjectConfiguration//i.KratosPath"
    set cproperty "dv"
    set KratosPath [::xmlutils::setXml $cxpath $cproperty]
    set KratosPath [file native $KratosPath]

    # Write the kratos path
    puts $fileid "kratos_path=\"${KratosPath}\"" 

    
}

proc ::wkcf::WriteFluidSolvers {rootid fileid vartype {wfsmethod 0}} {
    # Write fluid velocity and pressure solvers
    
    # Kratos key word xpath
    set kxpath "Applications/$rootid"
    # Set default value xml variable
    set cproperty "dv"

    puts $fileid "# $vartype solver"

    if {$wfsmethod} {
	# Define the class
	puts $fileid "class Fluid${vartype}LinearSolverConfiguration:"
    }

    set cxpath "$rootid//c.SolutionStrategy//i.${vartype}LinearSolverType"
    set LinearSolverType [::xmlutils::setXml $cxpath $cproperty]
    if {$LinearSolverType =="Direct"} {
	# Direct solver type
	set cxpath "$rootid//c.SolutionStrategy//i.${vartype}DirectSolverType"
	set DirectSolverType [::xmlutils::setXml $cxpath $cproperty]
	# WarnWinText "DirectSolverType:$DirectSolverType"
	set cDirectSolverType [::xmlutils::getKKWord $kxpath $DirectSolverType]
	
	if {$wfsmethod} {
	    puts $fileid "    solver_type = \"$cDirectSolverType\""
	} else {
	    puts $fileid "${vartype}_Linear_Solver = \"$cDirectSolverType\""
	}
	
    } elseif {$LinearSolverType =="Iterative"} {
	
	# Iterative solver type 
	set cxpath "$rootid//c.SolutionStrategy//i.${vartype}IterativeSolverType"
	set IterativeSolverType [::xmlutils::setXml $cxpath $cproperty]
	# Tolerance
	set cxpath "$rootid//c.SolutionStrategy//i.${vartype}ISTolerance"
	set Tolerance [::xmlutils::setXml $cxpath $cproperty]
	# Maximum iteration
	set cxpath "$rootid//c.SolutionStrategy//i.${vartype}ISMaximumIteration"
	set MaximumIteration [::xmlutils::setXml $cxpath $cproperty]
	# preconditioner type
	set cxpath "$rootid//c.SolutionStrategy//i.${vartype}PreconditionerType"
	set PreconditionerType [::xmlutils::setXml $cxpath $cproperty]
	# WarnWinText "vartype:$vartype IterativeSolverType:$IterativeSolverType Tolerance:$Tolerance MaximumIteration:$MaximumIteration PreconditionerType:$PreconditionerType"
	
	# Solver type
	set lsolver [::xmlutils::getKKWord $kxpath $IterativeSolverType]
	# Preconditioner
	set precond [::xmlutils::getKKWord $kxpath $PreconditionerType]
	
	if {$wfsmethod} {
	    puts $fileid "    solver_type = \"$lsolver\""
	    puts $fileid "    tolerance = $Tolerance"
	    puts $fileid "    max_iterations = $MaximumIteration"
	    puts $fileid "    preconditioner = \"$precond\""
	    puts $fileid "    use_scaling = True"
	} else {
	    puts $fileid "${vartype}_Linear_Solver = \"$lsolver\""
	    puts $fileid "${vartype}_Iterative_Tolerance = $Tolerance"
	    puts $fileid "${vartype}_Solver_Max_Iteration = $MaximumIteration"
	    puts $fileid "${vartype}_Preconditioner_type = \"$precond\""
	}
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
	set cxpath "${basexpath}//c.${copid}//c.HistoryGraph"
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
	set cxpath "${basexpath}//c.${ccutid}//c.Origin"
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
	set cxpath "${basexpath}//c.${ccutid}//c.NormalToPlane"
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
	set cxpath "${basexpath}//c.${cdpid}//c.MainProperties"
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

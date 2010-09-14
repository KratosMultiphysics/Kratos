###############################################################################
#
#    NAME: wkcfutils.tcl
#
#    PURPOSE: Some utilities procedures to write all the Kratos calculation files
#
#    QUANTECH ATZ-DEVELOPMENT DEPARTMENT
#
#    AUTHOR : G. Socorro
#
#    CREATED AT: 10/05/10
#
#    LAST MODIFICATION : correct an error when defined body force for group of element with the same property
#
#    VERSION : 0.9
#
#    HISTORY:
#
#     0.9- 07/09/10-G. Socorro, correct an error when defined body force for group of element with the same property
#     0.8- 06/09/10-G. Socorro, check for active group variable when get the element properties
#     0.7- 03/09/10-G. Socorro, correct an error with the BC for inlet and no-slip
#     0.6- 16/06/10-G. Socorro, add the elasto-plastic material model to the constitutive laws
#     0.5- 15/06/10-G. Socorro, add the damage material model to the constitutive laws
#     0.4- 11/06/10-G. Socorro, update material properties using the new constitutive equation
#     0.3- 09/06/10-G. Socorro, add FSI application
#     0.3- 31/05/10-G. Socorro, get all kratos key word from a XML file (kratos_key_words.xml)
#     0.2- 12/05/10-G. Socorro, start to add the fluid application options  
#     0.1- 10/05/10-G. Socorro, create the base source code from wkcf.tcl
#
###############################################################################

proc ::wkcf::Preprocess {} {
    # Create some global variables used when write data to the file
    # Get the problem dimension => spatial dimension
    variable ndime;    variable StructuralAnalysis
    variable gidetype; variable useqelem
    variable dprops;   variable FluidApplication
    variable FSIApplication; variable ActiveAppList

    # Check for use quadratic elements
    set useqelem [GiD_Info Project Quadratic]
    # WarnWinText "useqelem:$useqelem"

    # Get the spatial dimension
    set cxpath "GeneralApplicationData//c.Domain//i.SpatialDimension"
    set cproperty "dv"
    set ndime [::xmlutils::setXml $cxpath $cproperty]

    # Get application type
    # Structural analysis
    set cxpath "GeneralApplicationData//c.ApplicationTypes//i.StructuralAnalysis"
    set cproperty "dv"
    set StructuralAnalysis [::xmlutils::setXml $cxpath $cproperty]

    # WarnWinText "StructuralAnalysis:$StructuralAnalysis ndime:$ndime"

    # Fuild application
    set cxpath "GeneralApplicationData//c.ApplicationTypes//i.Fluid"
    set cproperty "dv"
    set FluidApplication [::xmlutils::setXml $cxpath $cproperty]

    # WarnWinText "FluidApplication:$FluidApplication"

    # FSI application
    set cxpath "GeneralApplicationData//c.ApplicationTypes//i.FluidStructureInteraction"
    set cproperty "dv"
    set FSIApplication [::xmlutils::setXml $cxpath $cproperty]

    # WarnWinText "FSIApplication:$FSIApplication"

    # Update active application list
    if {$FSIApplication =="Yes"} {
    set ActiveAppList [list "StructuralAnalysis" "Fluid"]
    # set structural analysis and fluid application to Yes
    set StructuralAnalysis "Yes"
    set FluidApplication "Yes"
    } else {
    if {$StructuralAnalysis =="Yes"} {
        set ActiveAppList [list "StructuralAnalysis"]
    } else {
        if {$FluidApplication =="Yes"} {
        set ActiveAppList [list "Fluid"]
        }
    }
    }
    
    # WarnWinText "ActiveAppList:$ActiveAppList"

    # Get the element properties
    ::wkcf::GetElementProperties
   
    # Get properties data
    ::wkcf::GetPropertiesData
   
    # Get boundary condition properties
    ::wkcf::GetBoundaryConditionProperties

    if {$StructuralAnalysis =="Yes"} {
    # Get load properties
    ::wkcf::GetLoadProperties
    }

    # Create the kratos global properties identifier
    ::wkcf::CreateKratosPropertiesIdentifier

    # ::WinUtils::PrintArray dprops
}

proc ::wkcf::CreateKratosPropertiesIdentifier {} {
    # Create the kratos global properties identifier
    variable dprops; variable StructuralAnalysis
    variable FluidApplication

    # Check for used body forces
    set usebforce "No"    
    if {$StructuralAnalysis =="Yes"} {
	set AppId "StructuralAnalysis"
	# List with all group that use body forces
	set dprops($AppId,AllBodyForceGroupId) [list]
	foreach cloadtid $dprops($AppId,AllLoadTypeId) {
	    if {$cloadtid =="BodyForce"} {
		# Check the group identifier
		if {([info exists dprops($AppId,Loads,$cloadtid,AllGroupId)]) && ([llength $dprops($AppId,Loads,$cloadtid,AllGroupId)]>0)} {
		    foreach cgroupid $dprops($AppId,Loads,$cloadtid,AllGroupId) {
			if {$cgroupid ni $dprops($AppId,AllKEGroupId)} {
			    # Error => In this version all body force must belong to the same an element group
			} else {
			    set usebforce "Yes"
			    if {$cgroupid ni $dprops($AppId,AllBodyForceGroupId)} {
				lappend dprops($AppId,AllBodyForceGroupId) $cgroupid
			    }
			}
		    }
		} else {
		    # Error: First define some group identifier for this body force
		}
		break
	    }
	}
	
	# WarnWinText "usebforce:$usebforce"

	# Create the global kratos properties list
	set dprops($AppId,GKProps,AllPropertyId) $dprops($AppId,AllKPropertyId)
	# For all defined kratos elements
	foreach celemid $dprops($AppId,AllKElemId) {
	    # WarnWinText "celemid:$celemid"
	    # For all defined group identifier for this elements
	    foreach cgroupid $dprops($AppId,KElem,$celemid,AllGroupId) {
		# Get the group properties
		lassign $dprops($AppId,KElem,$celemid,$cgroupid,GProps) GiDEntity GiDElemType PropertyId KEKWord nDim
		# WarnWinText "cgroupid:$cgroupid PropertyId:$PropertyId"
		if {![info exists dprops($AppId,GKProps,$PropertyId,AddBF)]} {
		    set dprops($AppId,GKProps,$PropertyId,AddBF) "No"
		}
		if {$usebforce =="Yes"} {
		    # Check to add body force properties
		    if {$cgroupid in $dprops($AppId,AllBodyForceGroupId)} {
			set dprops($AppId,GKProps,$PropertyId,AddBF) "Yes"
		    } 
		}
		
		# Update the global kratos property identifier 
		set GlobalPId [expr [lsearch $dprops($AppId,GKProps,AllPropertyId) $PropertyId]+1]
		set dprops($AppId,KElem,$celemid,$cgroupid,GlobalPId) $GlobalPId
	    }
	}
    }
    
    # For fluid application
    if {$FluidApplication =="Yes"} {
	set AppId "Fluid"
	# Create the global kratos properties list
	set dprops($AppId,GKProps,AllPropertyId) $dprops($AppId,AllKPropertyId)
	# For all defined kratos elements
	foreach celemid $dprops($AppId,AllKElemId) {
	    # For all defined group identifier for this elements
	    foreach cgroupid $dprops($AppId,KElem,$celemid,AllGroupId) {
		# Get the group properties
		# lassign $dprops($AppId,KElem,$celemid,$cgroupid,GProps) GiDEntity GiDElemType PropertyId KEKWord nDim
		
		# Update the global kratos property identifier 
		set GlobalPId 0
		set dprops($AppId,KElem,$celemid,$cgroupid,GlobalPId) $GlobalPId
	    }
	}
    }
}

proc ::wkcf::GetLoadProperties {} {
    # Get all load properties
    variable dprops

    # Set the application root identifier    
    set rootdataid "StructuralAnalysis"

    # Get all load properties
    set cxpath "$rootdataid//c.Loads"
    set clproplist [::xmlutils::setXmlContainerIds $cxpath]
    # WarnWinText "clproplist:$clproplist"
    # Load type list
    set dprops($rootdataid,AllLoadTypeId) [list]
    foreach cloadtid $clproplist {
    # WarnWinText "cloadtid:$cloadtid"
    # Get the group identifier defined for this load type
    set cxpath "$cxpath//c.${cloadtid}"
    set cgrouplist [::xmlutils::setXmlContainerIds $cxpath]
    # WarnWinText "current cgrouplist:$cgrouplist"
    if {[llength $cgrouplist]>0} {
        # Update load type identifier
        lappend dprops($rootdataid,AllLoadTypeId) $cloadtid
        # WarnWinText "inside cgrouplist:$cgrouplist"
        foreach cgroupid $cgrouplist {
        # Get the main properties
        set cxpath "$cxpath//c.${cgroupid}"
        set allmprop [::xmlutils::setXmlContainerPairs $cxpath "" "dv"]
        # WarnWinText "allmprop:$allmprop"
        
        # Kratos load to group link
        # Group list
        if {![info exists dprops($rootdataid,Loads,$cloadtid,AllGroupId)]} {
            set dprops($rootdataid,Loads,$cloadtid,AllGroupId) [list]
        }
        if {$cgroupid ni $dprops($rootdataid,Loads,$cloadtid,AllGroupId)} {
            lappend dprops($rootdataid,Loads,$cloadtid,AllGroupId) $cgroupid
        }
        # Group properties
        if {![info exists dprops($rootdataid,Loads,$cloadtid,$cgroupid,GProps)]} {
            set dprops($rootdataid,Loads,$cloadtid,$cgroupid,GProps) [list]
        }
        set dprops($rootdataid,Loads,$cloadtid,$cgroupid,GProps) $allmprop
        }
    }
    # Reset the path
    set cxpath "$rootdataid//c.Loads"
    }
}

proc ::wkcf::GetBoundaryConditionProperties {} {
    # Get all boundary condition properties
    variable dprops; variable ActiveAppList
    
    # For each active application
    foreach AppId $ActiveAppList {
    # Get the application root identifier    
    set rootdataid $AppId
        # Get all defined condition groups
    set cxpath "$rootdataid//c.Conditions"
    set cbcproplist [::xmlutils::setXmlContainerIds $cxpath]
    # WarnWinText "cbcproplist:$cbcproplist"
    # Boundary condition type list
    set dprops($AppId,AllBCTypeId) [list]
    foreach cbctid $cbcproplist {
        # WarnWinText "cbctid:$cbctid"
        # Get the group identifier defined for this condition
        set cxpath "${cxpath}//c.${cbctid}"
        set cbcgrouplist [::xmlutils::setXmlContainerIds $cxpath]
        # WarnWinText "cbcgrouplist:$cbcgrouplist"
        if {[llength $cbcgrouplist]>0} {
        # Update load type identifier
        lappend dprops($AppId,AllBCTypeId) $cbctid
        # WarnWinText "inside cbcgrouplist:$cbcgrouplist"
        foreach cgroupid $cbcgrouplist {
            set proplist [list]
            switch -exact -- $cbctid {
            "Displacements" - "Rotations" - "InletVelocity" - "No-Slip" {
                # Get activation properties
                foreach ca [list Ax Ay Az] cv [list Vx Vy Vz] {
                # Activation
                set acxpath "$cxpath//c.${cgroupid}//c.Activation//i.${ca}"
                # WarnWinText "Activation :cxpath:$cxpath"
                set cproperty "dv"
                set CActive [::xmlutils::setXml $acxpath $cproperty]
                # Values
                set vcxpath "$cxpath//c.${cgroupid}//c.Values//i.${cv}"
                # WarnWinText "Values :cxpath:$cxpath"
                set cproperty "dv"
                set CValue [::xmlutils::setXml $vcxpath $cproperty]
                lappend proplist $CActive $CValue
                }
            }
            "OutletPressure" {
                # Get properties
                foreach citem [list "FixPressure" "PressureValue"] {
                # set xpath
                set pcxpath "$cxpath//c.${cgroupid}//c.MainProperties//i.${citem}"
                set cproperty "dv"
                set CValue [::xmlutils::setXml $pcxpath $cproperty]
                lappend proplist $CValue
                }
            }
            }
            # WarnWinText "proplist:$proplist"
            # Kratos BC to group link
            # Group list
            if {![info exists dprops($AppId,BC,$cbctid,AllGroupId)]} {
            set dprops($AppId,BC,$cbctid,AllGroupId) [list]
            }
            if {$cgroupid ni $dprops($AppId,BC,$cbctid,AllGroupId)} {
            lappend dprops($AppId,BC,$cbctid,AllGroupId) $cgroupid
            }
            # Group properties
            if {![info exists dprops($AppId,BC,$cbctid,$cgroupid,GProps)]} {
            set dprops($AppId,BC,$cbctid,$cgroupid,GProps) [list]
            }
            set dprops($AppId,BC,$cbctid,$cgroupid,GProps) $proplist
        }
        }
        # Reset the path
        set cxpath "$rootdataid//c.Conditions"
    }
    }
}

proc ::wkcf::GetApplicationRootId {} {
    # Get the application root identifier
    variable FluidApplication;  variable StructuralAnalysis
    
    # Select the current active application
    if {$StructuralAnalysis =="Yes"} {
    set rootdataid "StructuralAnalysis"
    } else {
    if {$FluidApplication =="Yes"} {
        set rootdataid "Fluid"
    }
    }
    return $rootdataid
}

proc ::wkcf::GetPropertiesData {} {
    # Process all properties
    variable dprops; variable ActiveAppList
    variable ndime

    # For each active application
    foreach AppId $ActiveAppList {
    # Get the application root identifier    
    set rootdataid $AppId
    # Get the properties identifier 
    set cxpath "$rootdataid//c.Properties"
    set cproplist [::xmlutils::setXmlContainerIds $cxpath]
    # WarnWinText "cproplist:$cproplist"
    # All material list
    set dprops($AppId,AllMatId) [list]
    # Get the properties
    foreach propid $cproplist {
        # Material identifier
        set mxpath "$cxpath//c.${propid}//c.MainProperties//i.Material"
        set cproperty "dv"
        set MatId [::xmlutils::setXml $mxpath $cproperty]
        # WarnWinText "MatId:$MatId"
        set dprops($AppId,Property,$propid,MatId) "$MatId"
        # Get the material properties
        if {$MatId ni $dprops($AppId,AllMatId)} {
        lappend dprops($AppId,AllMatId) $MatId
        }
           
        # Thickness value
        set txpath "$cxpath//c.${propid}//c.MainProperties//i.Thickness"
        set cproperty "dv"
        set Thickness [::xmlutils::setXml $txpath $cproperty]
        # WarnWinText "Thickness:$Thickness"
        set dprops($AppId,Property,$propid,Thickness) $Thickness

        # Property type => Base element type
        set ptypexpath "$cxpath//c.${propid}//c.MainProperties//i.ElemType"
        set cproperty "dv"
        set ptype [::xmlutils::setXml $ptypexpath $cproperty]
        # WarnWinText "ptype:$ptype"
        set dprops($AppId,Property,$propid,BaseElemType) $ptype

        # Material model 
        set xpath "$cxpath//c.${propid}//c.MainProperties//i.MatModel"
        set cproperty "dv"
        set MatModel [::xmlutils::setXml $xpath $cproperty]
        # WarnWinText "MatModel:$MatModel"
        set dprops($AppId,Property,$propid,MatModel) $MatModel

        # Set fluency and behavior variables
        set dprops($AppId,Material,$MatId,UseFluency) "No"
        set dprops($AppId,Material,$MatId,Fluency) ""
        set dprops($AppId,Material,$MatId,UseBehavior) "No"
        set dprops($AppId,Material,$MatId,Behavior) ""

        # Get material properties
        switch -exact -- $MatModel {
        "Elastic-Isotropic" {
            if {($ptype=="PlaneStrain") && ($ndime =="2D")} {
            # Get the material properties
            ::wkcf::GetMaterialProperties $AppId $propid $MatId $ptype $MatModel "One"
            } elseif {($ptype=="PlaneStress") && ($ndime =="2D")} {
            set cptype "Isotropic2D"
            # Get the material properties
            ::wkcf::GetMaterialProperties $AppId $propid $MatId $cptype $MatModel "Yes"
            } elseif {(($ptype=="Solid")||($ptype=="Shell")) && ($ndime =="3D")} {
            set cptype "Isotropic3D"
            # Get the material properties
            ::wkcf::GetMaterialProperties $AppId $propid $MatId $cptype $MatModel "No"
            }
        }
        "Elastic-Orthotropic" {
        }
        "Elasto-Plastic" {
            if {($ptype=="PlaneStrain") && ($ndime =="2D")} {
            set cptype "Plasticity2D"
            # Get the material properties
            ::wkcf::GetMaterialProperties $AppId $propid $MatId $cptype $MatModel "One"
            # Get the material behavior and fluency properties
            ::wkcf::GetBehaviorFluencyProperties $AppId $MatId $MatModel $ptype $cptype
            } elseif {($ptype=="PlaneStress") && ($ndime =="2D")} {
            set cptype "Plasticity2D"
            # Get the material properties
            ::wkcf::GetMaterialProperties $AppId $propid $MatId $cptype $MatModel "Yes"
            # Get the material behavior and fluency properties
            ::wkcf::GetBehaviorFluencyProperties $AppId $MatId $MatModel $ptype $cptype
            } elseif {($ptype=="Solid") && ($ndime =="3D")} {
            set cptype "Plasticity3D"
            # Get the material properties
            ::wkcf::GetMaterialProperties $AppId $propid $MatId $cptype $MatModel "No"
            # Get the material behavior and fluency properties
            ::wkcf::GetBehaviorFluencyProperties $AppId $MatId $MatModel $ptype $cptype
            }
        }
        "Damage" {
            if {($ptype=="PlaneStrain") && ($ndime =="2D")} {
            set cptype "IsotropicDamage"
            # Get the material properties
            ::wkcf::GetMaterialProperties $AppId $propid $MatId $cptype $MatModel "One"
            # Get the material behavior and fluency properties
            ::wkcf::GetBehaviorFluencyProperties $AppId $MatId $MatModel $ptype $cptype
            } elseif {($ptype=="PlaneStress") && ($ndime =="2D")} {
            set cptype "IsotropicDamage"
            # Get the material properties
            ::wkcf::GetMaterialProperties $AppId $propid $MatId $cptype $MatModel "Yes"
            # Get the material behavior and fluency properties
            ::wkcf::GetBehaviorFluencyProperties $AppId $MatId $MatModel $ptype $cptype
            } elseif {($ptype=="Solid") && ($ndime =="3D")} {
            set cptype "IsotropicDamage3D"
            # Get the material properties
            ::wkcf::GetMaterialProperties $AppId $propid $MatId $cptype $MatModel "No"
            # Get the material behavior and fluency properties
            ::wkcf::GetBehaviorFluencyProperties $AppId $MatId $MatModel $ptype $cptype
            }
        }
        }
    }
    }
}

proc ::wkcf::GetBehaviorFluencyProperties {AppId MatId MatModel cptype ptype} {
    # Get the behavior and fluency material properties
    variable dprops
    variable ndime

    # WarnWinText "AppId:$AppId MatId:$MatId cptype:$cptype ptype:$ptype"
    # Xpath for constitutive laws
    set clxpath "CLawProperties"
    # Get all material properties
    set mpxpath "[::KMat::findMaterialParent $MatId]//m.${MatId}"
    # WarnWinText "mpxpath:$mpxpath"

    # Set fluency and behavior variables
    set dprops($AppId,Material,$MatId,UseFluency) "Yes"
    set dprops($AppId,Material,$MatId,UseBehavior) "Yes"
    # Get the softening behavior
    set mbehavior [::xmlutils::getKKWord $clxpath $ptype "mbehavior"]
    # Get the softening behavior xpath values
    set mbxpath [::xmlutils::getKKWord $clxpath $ptype "mbxpath"]
    # WarnWinText "mbehavior:$mbehavior mbxpath:$mbxpath"
    # Get the current behavior 
    set cbvalue [lindex [::KMat::getMaterialProperties "p" "$mpxpath//$mbxpath//p.$mbehavior"] 0 1]
    # Get the internal behavior properties
    set mbivalues [split [::xmlutils::getKKWord $clxpath $ptype "mbivalues"] ,]
    # Get the write behavior properties
    set mbwritev [split [::xmlutils::getKKWord $clxpath $ptype "mbwritev"] ,]
    # WarnWinText "mbwritev:$mbwritev mbivalues:$mbivalues\n$mpxpath//$mbxpath//p.$mbehavior cbvalue:$cbvalue"
    foreach mbiv $mbivalues mbwv $mbwritev {
    if {$mbiv ==$cbvalue} {
        set dprops($AppId,Material,$MatId,Behavior) "$mbwv"
        break
    }
    }
    # WarnWinText "dprops($AppId,Material,$MatId,Behavior):$dprops($AppId,Material,$MatId,Behavior)"
    if {$MatModel =="Damage"} {
    # Damage models
    # Get the energy yield function
    # Get the internal state properties
    set msivalues [split [::xmlutils::getKKWord $clxpath "MState" "msivalues"] ,]
    # Get the write behavior properties
    set mswritev [split [::xmlutils::getKKWord $clxpath "MState" "mswritev"] ,]
    # WarnWinText "mswritev:$mswritev msivalues:$msivalues"
    foreach msiv $msivalues mswv $mswritev {
        # WarnWinText "msiv:$msiv cptype:$cptype mswv:$mswv" 
        if {$msiv ==$cptype} {
        set dprops($AppId,Material,$MatId,Fluency) "EnergyYieldFunction(myState.${mswv})"
        break
        }
    }
    } elseif {$MatModel == "Elasto-Plastic"} {
    # Elasto-plastic models
    # Get the internal state properties
    set msivalues [split [::xmlutils::getKKWord $clxpath "MState" "msivalues"] ,]
    # Get the write behavior properties
    set mswritev [split [::xmlutils::getKKWord $clxpath "MState" "mswritev"] ,]
    # WarnWinText "mswritev:$mswritev msivalues:$msivalues"
    set cstate ""
    foreach msiv $msivalues mswv $mswritev {
        # WarnWinText "msiv:$msiv cptype:$cptype mswv:$mswv" 
        if {$msiv ==$cptype} {
        set cstate "myState.${mswv}"
        break
        }
    }
    # WarnWinText "cstate:$cstate"
    # Get the yield function properties
    # Get the yield criteria
    set yfid "YieldFunctions"
    set myieldcriteria [::xmlutils::getKKWord "$clxpath" $ptype "myieldcriteria"]
    # Get the yield criteria xpath values
    set mycxpath [::xmlutils::getKKWord "$clxpath" $ptype "mycxpath"]
    # Get the current yield criteria
    set cycvalue [lindex [::KMat::getMaterialProperties "p" "$mpxpath//$mycxpath//p.$myieldcriteria"] 0 1]
    # WarnWinText "myieldcriteria:$myieldcriteria mycxpath:$mycxpath cycvalue:$cycvalue"
    # Get the yield function options
    set yfivalues [split [::xmlutils::getKKWord "$clxpath//$yfid" "AvailableYieldFunction" "yfivalues"] ,]
    # Get the write yield function properties
    set yfwritev [split [::xmlutils::getKKWord "$clxpath//$yfid" "AvailableYieldFunction" "yfwritev"] ,]
    # WarnWinText "yfwritev:$yfwritev yfivalues:$yfivalues"
    set cyf ""
    foreach yfiv $yfivalues yfwv $yfwritev {
        # WarnWinText "yfiv:$yfiv cycvalue:$cycvalue yfwv:$yfwv" 
        if {$yfiv ==$cycvalue} {
        set cyf "$yfwv"
        break
        }
    }
    # WarnWinText "cyf:$cyf"

    set dprops($AppId,Material,$MatId,Fluency) "${cyf}(${cstate},myPotencialPlastic.Associated)"
    }
    # WarnWinText "dprops($AppId,Material,$MatId,Fluency):$dprops($AppId,Material,$MatId,Fluency)"
}

proc ::wkcf::GetMaterialProperties {AppId propid MatId ptype CMatModel {usethick "No"}} {
    # Get material properties
    variable dprops
    variable ndime
    
    # WarnWinText "MaterialProperties =>AppId:$AppId MatId:$MatId ptype:$ptype CMatModel:$CMatModel usethick:$usethick"
    # Xpath for constitutive laws
    set clxpath "CLawProperties"
    # Xpath for materials
    set matxpath "Materials"
  
    # Get all material properties
    set mpxpath "[::KMat::findMaterialParent $MatId]//m.${MatId}"
    # WarnWinText "mpxpath:$mpxpath"
    
    # Set the kratos model keyword base xpath
    set kmxpath "Applications//$AppId"
    
    # Material model
    set MatModel [::xmlutils::getKKWord $clxpath $ptype "matm"]
    set dprops($AppId,Material,$MatId,MatModel) "$MatModel"
    # Get the used material properties
    set mprops [split [::xmlutils::getKKWord $clxpath $ptype "mprops"] ,]
    # Get xpath values
    set mxpath [split [::xmlutils::getKKWord $clxpath $ptype "mxpath"] ,]
    # WarnWinText "MatModel:$MatModel mprops:$mprops mxpath:$mxpath"
    set matplist [list]
    foreach pid $mprops xpath $mxpath {
    # Get the kratos key word
    set kkword [::xmlutils::getKKWord $matxpath $pid "kkword"] 
    # WarnWinText "pid:$pid kkword:$kkword"
    # Get the current value for this properties
    # WarnWinText "xpath:$mpxpath//$xpath//p.$pid"
    set cvalue [lindex [::KMat::getMaterialProperties "p" "$mpxpath//$xpath//p.$pid"] 0 1]
    if {($kkword !="") && ($cvalue !="")} {
        lappend matplist [list $kkword $cvalue]
    }
    }

    # Add others properties for specific constitutive models
    if {$CMatModel == "Elasto-Plastic"} {
    # Get the yield function properties
    set yfid "YieldFunctions"
    # Get the yield criteria
    set myieldcriteria [::xmlutils::getKKWord $clxpath $ptype "myieldcriteria"]
    # Get the yield criteria xpath values
    set mycxpath [::xmlutils::getKKWord $clxpath $ptype "mycxpath"]
    # Get the current yield criteria
    set cycvalue [lindex [::KMat::getMaterialProperties "p" "$mpxpath//$mycxpath//p.$myieldcriteria"] 0 1]
    # WarnWinText "myieldcriteria:$myieldcriteria mycxpath:$mycxpath cycvalue:$cycvalue"
    # Get the yield function options
    set cyf ""
    set yfivalues [split [::xmlutils::getKKWord "$clxpath//$yfid" "AvailableYieldFunction" "yfivalues"] ,]
    foreach yfiv $yfivalues {
        if {$yfiv ==$cycvalue} {
        set cyf "$yfiv"
        break
        }
    }
    # WarnWinText "cyf:$cyf"
    # Get other properties for the specific yield function
    if {$cyf !=""} {
        # Get the yield function parameters
        set privalues [split [::xmlutils::getKKWord "$clxpath//$yfid" "$cyf" "privalues"] ,]
        # Get the write yield function parameters 
        set prxpath [split [::xmlutils::getKKWord "$clxpath//$yfid" "$cyf" "prxpath"] ,]
        # WarnWinText "$clxpath//$cyf prxpath:$prxpath privalues:$privalues"
        foreach pid $privalues xpath $prxpath {
        # Get the kratos key word
        set kkword [::xmlutils::getKKWord $matxpath $pid "kkword"] 
        # WarnWinText "pid:$pid kkword:$kkword"
        # Get the current value for this properties
        # WarnWinText "xpath:$mpxpath//$xpath//p.$pid"
        set cvalue [lindex [::KMat::getMaterialProperties "p" "$mpxpath//$xpath//p.$pid"] 0 1]
        if {($kkword !="") && ($cvalue !="")} {
            lappend matplist [list $kkword $cvalue]
        }
        }
    }
    }

    # WarnWinText "matplist:$matplist"
    set dprops($AppId,Material,$MatId,Props) $matplist
    # Get others properties
    # WarnWinText "clxpath:$clxpath ptype:$ptype"
    set cprops [::xmlutils::getKKWord $clxpath $ptype "cprops"]
    # WarnWinText "cprops:$cprops"
    # Check to use thickness value
    set dprops($AppId,Material,$MatId,CProps) [list]
    switch -exact -- $usethick {
    "Yes" {
        # Get the current value
        set Thickness $dprops($AppId,Property,$propid,Thickness)
        set kword [::xmlutils::getKKWord $kmxpath [lindex $cprops 0] "kkword"]
        # Update section properties
        set dprops($AppId,Material,$MatId,CProps) [list [list $kword $Thickness]]
        # WarnWinText "dprops($AppId,Material,$MatId,CProps):$dprops($AppId,Material,$MatId,CProps)"
    }
    "One" {
        set Thickness "1.0"
        set kword [::xmlutils::getKKWord $kmxpath [lindex $cprops 0] "kkword"]
        # Update section properties
        set dprops($AppId,Material,$MatId,CProps) [list [list $kword $Thickness]]
        # WarnWinText "dprops($AppId,Material,$MatId,CProps):$dprops($AppId,Material,$MatId,CProps)"
    }
    }
}

proc ::wkcf::GetElementProperties {} {
    variable dprops
    variable ActiveAppList

    # For each active application
    foreach AppId $ActiveAppList {
	# Get the application root identifier    
	set rootdataid $AppId
	# Get the properties element links
	# All kratos element identifier
	set dprops($AppId,AllKElemId) [list]
	# All kratos property identifier
	set dprops($AppId,AllKPropertyId) [list]
	# All Kratos element group identifier
	set dprops($AppId,AllKEGroupId) [list]
	# Get all defined element types
	set cxpath "$rootdataid//c.Elements"
	set glist [::xmlutils::setXmlContainerIds $cxpath]
	set kwxpath "Applications/$rootdataid"
	# WarnWinText "glist:$glist"
	foreach celemid $glist {
	    # Get the group identifier defined for this element 
	    set cxpath "$rootdataid//c.Elements//c.${celemid}"
	    set cgrouplist [::xmlutils::setXmlContainerIds $cxpath]
	    if {[llength $cgrouplist]>0} { 
		# WarnWinText "cgrouplist:$cgrouplist"
		foreach cgroupid $cgrouplist {
		    # Check if this group is in active state
		    # Get active propery
		    set cxpath "$rootdataid//c.Elements//c.${celemid}//c.${cgroupid}"
		    set cproperty "active"
		    set ActiveGroup [::xmlutils::setXml $cxpath $cproperty]
		    # WarnWinText "ActiveGroup:$ActiveGroup"
		    if {$ActiveGroup =="0"} {
			continue
		    }
		    # Update all kratos element identifier
		    if {$celemid ni $dprops($AppId,AllKElemId)} {
			lappend dprops($AppId,AllKElemId) $celemid
		    }
		    # Update all group identifier
		    if {$cgroupid ni $dprops($AppId,AllKEGroupId)} {
			lappend dprops($AppId,AllKEGroupId) $cgroupid
		    }
		    # Get the GiD entity type
		    set cxpath "$rootdataid//c.Elements//c.${celemid}"
		    set cproperty "GiDEntity"
		    set GiDEntity [::xmlutils::setXml $cxpath $cproperty]
		    # WarnWinText "GiDEntity:$GiDEntity"
		    # Get the GiD element type for Kratos elements type
		    set cxpath "$rootdataid//c.Elements//c.${celemid}//c.${cgroupid}//c.Properties//i.ElementType"
		    set cproperty "dv"
		    set GiDElemType [::xmlutils::setXml $cxpath $cproperty]
		    # WarnWinText "GiDElemType:$GiDElemType"
		    # Get the property identifier
		    set cxpath "$rootdataid//c.Elements//c.${celemid}//c.${cgroupid}//c.Properties//i.Property"
		    set cproperty "dv"
		    set PropertyId [::xmlutils::setXml $cxpath $cproperty]
		    # WarnWinText "PropertyId:$PropertyId"
		    # Get the Key word for the element type
		    set kelemtype [::xmlutils::getKKWord $kwxpath $celemid "kkword"]
		    # set cxpath "$rootdataid//c.Elements//c.${celemid}"
		    # set cproperty "kkword"
		    # set kelemtype [::xmlutils::setXml $cxpath $cproperty]
		    # WarnWinText "kelemtype:$kelemtype"
		    # Get ndim
		    set cxpath "$rootdataid//c.Elements//c.${celemid}"
		    set cproperty "nDim"
		    set nDim [::xmlutils::setXml $cxpath $cproperty]
		    # WarnWinText "nDim:$nDim"
		    
		    set GProps [list $GiDEntity $GiDElemType $PropertyId $kelemtype $nDim]
		    
		    # Kratos element to group link
		    # Group list
		    if {![info exists dprops($AppId,KElem,$celemid,AllGroupId)]} {
			set dprops($AppId,KElem,$celemid,AllGroupId) [list]
		    }
		    if {$cgroupid ni $dprops($AppId,KElem,$celemid,AllGroupId)} {
			lappend dprops($AppId,KElem,$celemid,AllGroupId) $cgroupid
		    }
		    # Group properties
		    if {![info exists dprops($AppId,KElem,$celemid,$cgroupid,GProps)]} {
			set dprops($AppId,KElem,$celemid,$cgroupid,GProps) [list]
		    }
		    set dprops($AppId,KElem,$celemid,$cgroupid,GProps) $GProps
		    
		    # Update AllKPropertyId list
		    if {$PropertyId ni $dprops($AppId,AllKPropertyId)} {
			lappend dprops($AppId,AllKPropertyId) $PropertyId
		    }
		    # Property to group identifier link
		    if {![info exists dprops($AppId,Property,$PropertyId,GroupId)]} {
			set dprops($AppId,Property,$PropertyId,GroupId) [list]
		    }
		    # Update the link between group and property => if the property exist add to the list
		    if {[info exists dprops($AppId,Property,$PropertyId,GroupId)]} {
			lappend dprops($AppId,Property,$PropertyId,GroupId) $cgroupid
		    } else {
			set dprops($AppId,Property,$PropertyId,GroupId) $cgroupid
		    }
		}
	    }
	}
    }
}

proc ::wkcf::UnsetLocalVariables {} {
    variable dprops

    if {[info exists dprops]} {
    unset dprops
    }
}

proc ::wkcf::GetnDimnNode {GiDElemType nDim} {
    # ABSTRACT: Get the kratos element identifier as a function of GiD element type
    variable useqelem

    set etbf ""
    set nid "3N"
    if {$nDim =="2D"} {
    if {$GiDElemType =="Triangle"} {
        if {$useqelem=="1"} {
        set nid "6N"
        } else {
        set nid "3N"
        }
    } elseif {$GiDElemType =="Quadrilateral"} {
        if {$useqelem=="1"} {
        set nid "8N"
        } else {
        set nid "4N"
        }
    }
    } elseif {$nDim =="3D"} {
    if {$GiDElemType =="Tetrahedra"} {
        if {$useqelem=="1"} {
        set nid "10N"
        } else {
        set nid "4N"
        }
    } elseif {$GiDElemType =="Hexahedra"} {
        if {$useqelem=="1"} {
        set nid "20N"
        } elseif {$useqelem=="2"} {
        set nid "27N"
        } else {
        set nid "8N"
        }
    }
    }
    set etbf "${nDim}${nid}"
    return $etbf
}

proc ::wkcf::FindBoundaries {entity} {
    # ABSTRACT: Return a list containing all boundaries entities
    # Arguments
    # entity => Intity to be processed
    #  * entity=line for models made of surfaces 
    #  * entity=surface for models made of volumes
    # Note: This procedure in the same used in the fluid_only problem type

    set boundarylist [list]
    # Generate some names
    set Entity [string toupper $entity 0 0]
    set entities [format "%ss" $entity]
    
    # Get the number of the last entity
    set instruction [format "MaxNum%ss" $Entity]
    set Max [GiD_Info Geometry $instruction]
    
    # Generate a list containing all entities and record their id and number of HigerEntities
    set EntityList [GiD_Info list_entities $entities 1:$Max]
    set candidates [regexp -all -inline {Num: ([0-9]*) HigherEntity: ([0-9]*)} $EntityList]
    
    # Find ids of entities with exactly 1 HigherEntity (this means they are in the boundary)
    for {set i 1} {$i < [llength $candidates] } {incr i 3} {
    set j [expr {$i + 1}]
    if {[lindex $candidates $j] == 1} {lappend boundarylist [lindex $candidates $i]}
    }
    return $boundarylist
}


proc ::wkcf::AlignLineNormals {direction} {
    # ABSTRACT: Makes all of boundary lines' normals point inwards or outwards
    # Arguments
    # direction => Direction option ["Inwards"|"Outwards"]
    # Note: This procedure in the same used in the fluid_only problem type
    
    switch $direction {
	Inwards {
	    set wrong_way "DIFF1ST"
	}
	Outwards {
	    set wrong_way "SAME1ST"
	} 
	default {puts "Unknown direction, line normals not aligned"}
    }
    
    set surfacelist [GiD_Geometry list surface 1:]
    
    # For each surface, we look for boundary lines oriented in the wrong direction
    set linelist [list]
    foreach surface $surfacelist {
	set surfaceinfo [GiD_Info list_entities surfaces $surface]
	set numpos [lsearch $surfaceinfo "NumLines:"]
	set numlines [lindex $surfaceinfo [expr {$numpos +1}]]
	for {set i 0} {$i < $numlines} {incr i} {
	    set orient [lindex $surfaceinfo [expr {$numpos+5+4*$i}]]
	    if {[string compare $orient $wrong_way]==0} {
		# If the normal is pointing in the wrong direction, 
		# Check if it's a contour line
		set linenum [lindex $surfaceinfo [expr {$numpos+3+4*$i}]]
		set lineinfo [GiD_Info list_entities lines $linenum]
		#set highpos [lsearch $surfinfo "HigherEntity:"]
		set higherentities [lindex $lineinfo 4]
		if {$higherentities==1} {
		    lappend linelist $linenum
		}
	    }
	}
    }
    # WarnWinText "linelist:$linelist"
    if {[llength $linelist]} {
	# If its in the contour, switch its normal
	eval GiD_Process Mescape Utilities SwapNormals Lines Select $linelist
    }

}

proc ::wkcf::AlignSurfNormals {direction} {
    # ABSTRACT: Makes all of boundary surfaces' normals point inwards or outwards
    # Arguments
    # direction => Direction option ["Inwards"|"Outwards"]
    # Note: This procedure in the same used in the fluid_only problem type
    
    switch $direction {
	Inwards {
	    set wrong_way "DIFF1ST"
	}
	Outwards {
	    set wrong_way "SAME1ST"
	} 
	default {puts "Unknown Direction, surface normals not aligned"}
    }
    
    set volumelist [GiD_Geometry list volume 1:]

    set surfacelist [list] 
    # For each volume, we look for face surfaces with oriented in the wrong direction
    foreach volume $volumelist {
	set volumeinfo [GiD_Info list_entities volumes $volume]
	set numpos [lsearch $volumeinfo "NumSurfaces:"]
	set numsurf [lindex $volumeinfo [expr {$numpos +1 }]]
	for {set i 0} {$i < $numsurf} {incr i} {
	    set orient [lindex $volumeinfo [expr {$numpos+5+4*$i}]]
	    if {[string compare $orient $wrong_way]==0} {
		# If the normal is pointing in the wrong direction,
		# Check if it's a contour surface
		set surfnum [lindex $volumeinfo [expr {$numpos+3+4*$i}]]
		set surfinfo [GiD_Info list_entities surfaces $surfnum]
		set higherentities [lindex $surfinfo 4]
		if {$higherentities==1} {
		 lappend surfacelist $surfnum
		}
	    }
	}
    }

    if {[llength $surfacelist]} {
	# If its in the contour, switch its normal
	eval GiD_Process Mescape Utilities SwapNormals Surfaces Select $surfacelist
    }
}

proc ::wkcf::WriteBatFile {} {
  
    set batfilename "Kratos.win.bat"
    set ProblemTypePath [::KUtils::GetPaths "PTDir"]
    set batfullname [file native [file join $ProblemTypePath $batfilename]]

    # First delete the file
    set res ""
    catch { set res [file delete -force $batfullname] }
    
    # Create the new file
    set f [open $batfullname w]
    # WarnWinText "batfullname:$batfullname res:$res f:$f"

    puts $f "REM @ECHO OFF"
    puts $f "REM Identification for arguments"
    puts $f "REM basename                          = %1"
    puts $f "REM Project directory                 = %2"
    puts $f "REM Problem directory                 = %3"
    puts $f " "
    puts $f "REM OutputFile: %2\\%1.info"
    puts $f "REM ErrorFile: %2\\%1.err"
    puts $f " "


    puts $f "DEL %2\\%1.info"
    puts $f "DEL %2\\%1.post.bin"
    puts $f "DEL %2\\%1.err"

    puts $f "REM Run the python script"
    puts $f "python KratosOpenMP.py > %2\\%1.info 2> %2\\%1.err"

    close $f
}

proc ::wkcf::WriteFluidSolvers {rootid fileid vartype} {
    # Write fluid velocity and pressure solvers
    
    # Kratos key word xpath
    set kxpath "Applications/$rootid"
    # Set default value xml variable
    set cproperty "dv"

    puts $fileid "# $vartype solver"
    
    set cxpath "$rootid//c.SolutionStrategy//i.${vartype}LinearSolverType"
    set LinearSolverType [::xmlutils::setXml $cxpath $cproperty]
    if {$LinearSolverType =="Direct"} {
	# Direct solver type
	set cxpath "$rootid//c.SolutionStrategy//i.${vartype}DirectSolverType"
	set DirectSolverType [::xmlutils::setXml $cxpath $cproperty]
	# WarnWinText "DirectSolverType:$DirectSolverType"
	set cDirectSolverType [::xmlutils::getKKWord $kxpath $DirectSolverType]
	puts $fileid "${vartype}_Linear_Solver = \"$cDirectSolverType\""
    
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
	# WarnWinText "IterativeSolverType:$IterativeSolverType Tolerance:$Tolerance MaximumIteration:$MaximumIteration PreconditionerType:$PreconditionerType"
    
	# Solver type
	set lsolver [::xmlutils::getKKWord $kxpath $IterativeSolverType]
	puts $fileid "${vartype}_Linear_Solver = \"$lsolver\""
	puts $fileid "${vartype}_Iterative_Tolerance = $Tolerance"
	puts $fileid "${vartype}_Solver_Max_Iteration = $MaximumIteration"
	# Preconditioner
	set precond [::xmlutils::getKKWord $kxpath $PreconditionerType]
	puts $fileid "${vartype}_Preconditioner_type = \"$precond\""
    }
 }

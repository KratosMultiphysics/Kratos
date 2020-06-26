
package require snit
package require snitwiz

namespace eval ::kwiz:: {
     variable wprops

     # Base window path
     variable bwinpath

     # Store all step snit variable step identifier
     variable stepidlist

     # Store all wiz step
     variable wizstepidlist

     # Store the wizard id
     variable wizardid

     variable editedMesh

     # Calculation percentage
     variable ccalvalue

     # Caution, editing
     variable editing
}

proc ::kwiz::Init {} {
     # ABSTRACT: Init namespace some variables

     # Init the base window path inside GiD
     variable bwinpath
     set bwinpath ".gid.wtest"

     variable wprops
     set wprops(dummy) 0

     variable editedMesh
     set editedMesh 1

     # Init step id list
     variable stepidlist [list]
     set stepidlist [list ::kwiz::step1 ::kwiz::step1bis ::kwiz::step2 ::kwiz::step3 ::kwiz::step4 ::kwiz::step5 ::kwiz::step6 ::kwiz::step7 ::kwiz::step8]

     variable wizstepidlist [list]
     set wizstepidlist [list "Init" "TestType" "ExperimentType" "GeometryMesh" "MaterialParameters" "GeneralSettings" "OutputSettings" "RunProblem" "ResultView"]

     # Calculation percentage
     variable ccalvalue 0.0



}

# Init some namespace variables
::kwiz::Init


# Initialize default values if new button is pressed
proc ::kwiz::initdefvalues { } {
     variable wprops
     #wa "Defaulteando"
     # Default values step1
     if {[info exists wprops(TestType,TestType)]} {
          if {$wprops(TestType,TestType) ==""} {
               set wprops(TestType,TestType) $::defvals::dem::wizdata::TestType
          }
     } else {
          set wprops(TestType,TestType) $::defvals::dem::wizdata::TestType
     }

     # Default values step2
     if {[info exists wprops(ExperimentType,ExperimentType)]} {
          if {$wprops(ExperimentType,ExperimentType) ==""} {
               set wprops(ExperimentType,ExperimentType) $::defvals::dem::wizdata::ExperimentType
          }
     } else {
          set wprops(ExperimentType,ExperimentType) $::defvals::dem::wizdata::ExperimentType
     }

     # Default values step3
     if {[info exists wprops(GeometryMesh,MeshType)]} {
          if {$wprops(GeometryMesh,MeshType) ==""} {
               set wprops(GeometryMesh,MeshType) $::defvals::dem::wizdata::MeshType
          }
     } else {
          set wprops(GeometryMesh,MeshType) $::defvals::dem::wizdata::MeshType
     }

      # Set default material values
     ::kwiz::initmaterialdefaultvalues

     # Default values step5
     set plist [list CalculationTime DeltaTimeSelection DeltaTime LoadingVelocity]
     set pdflist [list $::defvals::dem::wizdata::CalculationTime $::defvals::dem::wizdata::DeltaTimeSelection \
                  $::defvals::dem::wizdata::DeltaTime $::defvals::dem::wizdata::LoadingVelocity]

     foreach pid $plist dfv $pdflist {

        if {[info exists wprops(GeneralSettings,$pid)]} {
            if {$wprops(GeneralSettings,$pid) ==""} {
                set wprops(GeneralSettings,$pid) $dfv
            }
        } else {
                set wprops(GeneralSettings,$pid) $dfv
        }
     }

     # Default values step6
     set plist [list Displacement Velocity TotalForces Rhs DampForces AppliedForces ContactSigma ContactTau \
                LocalContactForce FailureCriterionState SinglePostprocessFile PrintOutputFile PrintGraphData]
     set pdflist [list $::defvals::dem::wizdata::Displacement $::defvals::dem::wizdata::Velocity \
                  $::defvals::dem::wizdata::TotalForces $::defvals::dem::wizdata::Rhs $::defvals::dem::wizdata::DampForces \
                  $::defvals::dem::wizdata::AppliedForces $::defvals::dem::wizdata::ContactSigma $::defvals::dem::wizdata::ContactTau \
                  $::defvals::dem::wizdata::LocalContactForce $::defvals::dem::wizdata::FailureCriterionState \
                  $::defvals::dem::wizdata::SinglePostprocessFile $::defvals::dem::wizdata::PrintOutputFile \
                  $::defvals::dem::wizdata::PrintGraphData]

     foreach pid $plist dfv $pdflist {

        if {[info exists wprops(OutputSettings,$pid)]} {
            if {$wprops(OutputSettings,$pid) ==""} {
                set wprops(OutputSettings,$pid) $dfv
            }
        } else {
                set wprops(OutputSettings,$pid) $dfv
        }
     }

     # Default values step7
     set plist [list SelectOMPMPI NumberOfThreads]
     set pdflist [list $::defvals::dem::wizdata::SelectOMPMPI $::defvals::dem::wizdata::NumberOfThreads]

     foreach pid $plist dfv $pdflist {

        if {[info exists wprops(RunProblem,$pid)]} {
            if {$wprops(RunProblem,$pid) ==""} {
                set wprops(RunProblem,$pid) $dfv
            }
        } else {
                set wprops(RunProblem,$pid) $dfv
        }
     }

}

proc ::kwiz::WizardAfterLoad { filename } {
     # wa "$filename"
     # wa [file tail $filename]
     # wa [winfo exists ".gid.wtest.w.layoutFrame.wiz.layout"]

     if {[file tail $filename] ne "kratos.spd"} {
          set block [::kwiz::blocking "get"]
          if {$block ne "block"} {
               if {[winfo exists ".gid.wtest.w.layoutFrame.wiz.layout"]} {
                    global KPriv
                    # wa [$KPriv(xmlWiz) asXML]
                    # Automatically go to the next step
                    set currentStepId [::kwiz::getCurrentStep]
                    set currentStepId [::kwiz::getCurrentStepInnerName $currentStepId]

                    $::kwiz::bwinpath.w gotostep $currentStepId
               }
          }
     }
}

proc ::kwiz::SetMeshLenDiam { } {
     variable wprops
     global KPriv

     set ExperimentType $wprops(ExperimentType,ExperimentType)
     set meshtype $wprops(GeometryMesh,MeshType)
     # wa "ExperimentType $ExperimentType"
     # wa "meshtype $meshtype"
     if {$ExperimentType ne "BTS"} { set ET "UCS"} else {
          set ET "BTS"
          set aux $meshtype
          set meshtype $ExperimentType
          append meshtype $aux
     }
     # wa "ET $ET"
     set xmlData $KPriv(xmlWiz)
     set path "/Wizard/MeshDB"
     set nodes [$xmlData selectNodes "$path/Item"]

     foreach node $nodes {
          if {$ET eq [$node getAttribute experimenttype ""] } {
               if {$meshtype eq [$node getAttribute id ""] } {
                   # wa [$node asXML]
                    set wprops(GeometryMesh,ProbeLength) [$node getAttribute length ""]
                    set wprops(GeometryMesh,ProbeDiameter) [$node getAttribute diameter ""]
               }
          }
     }

     # wa "$wprops(GeometryMesh,ProbeLength) $wprops(GeometryMesh,ProbeDiameter)"
}

proc ::kwiz::getCurrentStepInnerName {StepId} {
     variable wizstepidlist
     variable stepidlist

     set i [lsearch $wizstepidlist $StepId]
     set paso [lindex $stepidlist $i]
     # wa $paso
     return $paso

}

proc ::kwiz::getCurrentStep {} {
     variable wizstepidlist
     foreach stepId $wizstepidlist {
          # wa "active for $stepId ? [::kwiz::StepState $stepId "Active"]"
          if {[::kwiz::StepState $stepId "Active"]} {return $stepId}
     }
     return "Init"
}

proc ::kwiz::StepState {StepId pos {action "get"} {value 1}} {
     global KPriv
     # wa [$KPriv(xmlWiz) asXML]
     # set xmlData $KPriv(xmlWiz)
     # set xmlRoot "/Wizard/Steps/Step\[@id='$StepId'\]/Data/Item\[@idv='$pos'\]"

     set stps [$KPriv(xmlWiz) getElementsByTagName "Step"]
     foreach step $stps {
          if {[$step getAttribute id ""] eq $StepId} {
               set things [[lindex [$step childNodes] 0 ] childNodes]
               foreach node $things {
                    # wa "pre [$node getAttribute idv ""]"
                    if {[$node getAttribute idv ""] eq $pos} {
                         if {$action eq "get"} {
                              return [$node getAttribute dv ""]
                         } else {
                              $node setAttribute dv $value
                         }
                    }
               }
          }
     }

     return 0
     #set node [$xmlData selectNodes $xmlRoot]
     #wa [$node asXML]
     #if {$action eq "get"} {
     #     return [$node getAttribute dv ""]
     #} else {
     #     $node setAttribute dv $value
     #}
}

proc ::kwiz::NextEvent {stepnumber} {
     variable wizstepidlist
     set stepId [lindex $wizstepidlist [expr {$stepnumber}] ]
     ::kwiz::StepState $stepId "Active" "set" 0
     set stepId [lindex $wizstepidlist [expr {$stepnumber +1}] ]
     ::kwiz::StepState $stepId "Active" "set"
     ::kwiz::StepState $stepId "Visited" "set"

     ::kwiz::UpdateEvent $stepnumber
}
proc ::kwiz::BackEvent {stepnumber} {
     variable wizstepidlist
     set stepId [lindex $wizstepidlist [expr {$stepnumber}] ]
     ::kwiz::StepState $stepId "Active" "set" 0
     set stepId [lindex $wizstepidlist [expr {$stepnumber -1}] ]
     ::kwiz::StepState $stepId "Active" "set"

     ::kwiz::UpdateEvent $stepnumber

}
proc ::kwiz::UpdateEvent {stepnumber} {
     variable wprops
     variable wizstepidlist

     if {$stepnumber ne "0"} {
          set stepId [lindex $wizstepidlist $stepnumber]

          if {$stepId eq "GeometryMesh"} {
               ::kwiz::SetMeshLenDiam
          }

          # wa "Next from $stepId"
          ::kfiles::SaveWizardDataStep $stepId

          ::kfiles::SaveWizardToSPDStep $stepId
     }
}

proc ::kwiz::SaveActiveWizardtoSPD { filename } {
     variable wizardid

     set cxpath "DEM//c.Wizards"
     set wizards [::xmlutils::setXmlContainerIds $cxpath "Item"]

     foreach wiz $wizards {
          set cxpath "DEM//c.Wizards//i.$wiz"
          # wa "wiz $wiz"
          # wa "Id $wizardid"
          if {$wiz eq $wizardid} {
               ::xmlutils::setXml $cxpath "dv" write 1
               set wizardFile "[string range $filename 0 [expr [string length $filename] - 5]].$wizardid.wiz"
               set filen [::KUtils::GiveProjectName]
               ::xmlutils::setXml $cxpath "rpath" write "$filen.$wizardid.wiz"
          } else {
               ::xmlutils::setXml $cxpath "dv" write 0
          }
     }
     return $wizardFile
}

proc ::kwiz::ValidateDV {idv dv} {

      if {$idv eq "ExperimentType"} {
          if {$dv eq "Uniaxial compression strength"} {return "UCS"}
          if {$dv eq "Uniaxial strain compaction"} {return "Oedometric"}
     }
     return $dv
}


# Action for new and load buttons (TestType)
proc ::kwiz::wizardAction {w b} {
     variable wprops

     if { $b == 1  } {
          # Poner todos los valores por defecto
          GiD_Process Files SaveAs
          #::kwiz::initdefvalues

          if {[::KUtils::GetPaths "PName"] ne ""} {
               # Automatically go to the next step
               $::kwiz::bwinpath.w gotostep ::kwiz::step1bis
          }

     } elseif { $b == 2  } {
          # Cargar valores de archivo
          set wprops(TestType,WizardType) 2
          GiD_Process Files Read
     }

     # Set the widgets
     set wprops(TestType,WizardType) $b

     ::kwiz::CheckNextButtonStep1
     #::kwiz::checkWizardAction $w

}

proc ::kwiz::CheckProjName { } {
     global KPriv

     set root [::KUtils::GetPaths "PDir"]
     if {$root eq ""} {set root $KPriv(problemTypeDir)}

     # wa $root
     set ProjectName [GiD_Info Project "ModelName"]
     # wa $ProjectName
     if {$ProjectName eq "UNNAMED"} {
          set ProjectName ""
     }
     return $ProjectName
}


proc ::kwiz::CheckNextButtonStep1 { } {
     variable wprops
     if {[info exists wprops(TestType,TestType) ]} {
          if { $wprops(TestType,TestType) == 1 } {
               .gid.wtest.w buttonstate next normal
          } elseif { $wprops(TestType,TestType) == 2 } {
               .gid.wtest.w buttonstate next normal
          } else {
               .gid.wtest.w buttonstate next disabled
          }
     }
}


# Check variables for new and load buttons (Step1)
proc ::kwiz::checkWizardAction { w } {
     variable wprops

     if {[info exists wprops(TestType,WizardType)]} {
          if {$wprops(TestType,WizardType) ne "0"} {
               ::kwiz::HideButtonsStep1 $w
               # Set the widgets
               set img1 [::WinUtils::GetImage MaterialTest2.gif]
               set img2 [::WinUtils::GetImage CuttingTest2.gif]
               set but1 [ttk::checkbutton $w.b1 -style TMenubutton.Toolbutton -image $img1 -command "::kwiz::buttonPushedStep1 $w" -variable ::kwiz::wprops(TestType,TestType) -onvalue 1]
               set but2 [ttk::checkbutton $w.b2 -style TMenubutton.Toolbutton -image $img2 -command "::kwiz::buttonPushedStep1 $w" -variable ::kwiz::wprops(TestType,TestType) -onvalue 2]

               # Check image buttons
               ::kwiz::checkImageStep1 $w
               ::kwiz::CheckNextButtonStep1
               # Set the grid
               grid $but1 -column 1 -row 2
               grid $but2 -column 2 -row 2

          }
     }
}

# Proc for check image buttons
proc ::kwiz::checkImageStep1 { w } {
     variable wprops
     set img3 [::WinUtils::GetImage MaterialTestSel.gif]
     set img4 [::WinUtils::GetImage CuttingTestSel.gif]
     if {[info exists wprops(TestType,TestType) ]} {

          if { $wprops(TestType,TestType) == 1 } {
               set ::kwiz::wizardid "compressiontest"
               $w.b1 configure -image $img3
              .gid.wtest.w buttonstate next normal
          }
          if { $wprops(TestType,TestType) == 2 } {
               set ::kwiz::wizardid "cutting"
               $w.b2 configure -image $img4
              .gid.wtest.w buttonstate next normal
          }
     } else {
          .gid.wtest.w buttonstate next disabled
     }
}


proc ::kwiz::checkImageStep2 { w } {
    variable wprops
     set imgUCSm [::WinUtils::GetImage UCSMarked.gif]
     set imgTrim [::WinUtils::GetImage TriMarked.gif]
     set imgUSCm [::WinUtils::GetImage USCMarked.gif]
     set imgHydm [::WinUtils::GetImage TriMarked.gif]
     set imgBTSm [::WinUtils::GetImage BTSMarked.gif]

    if {[info exists wprops(ExperimentType,ExperimentType) ]} {
        if { $wprops(ExperimentType,ExperimentType) eq "UCS" } {
            $w.b1 configure -image $imgUCSm
        }
        if { $wprops(ExperimentType,ExperimentType) eq "Triaxial" } {
            $w.b2 configure -image $imgTrim
        }
        if { $wprops(ExperimentType,ExperimentType) eq "Oedometric" } {
            $w.b3 configure -image $imgUSCm
        }
        if { $wprops(ExperimentType,ExperimentType) eq "Hydrostatic" } {
            $w.b4 configure -image $imgHydm
        }
        if { $wprops(ExperimentType,ExperimentType) eq "BTS" } {
            $w.b5 configure -image $imgBTSm
        }
    }
}

proc ::kwiz::checkImageStep3 { w } {
    variable wprops

     set img1 [::WinUtils::GetImage bts_mesh_marked.gif]
     set img2 [::WinUtils::GetImage ucs_mesh_marked.gif]
     if {[info exists wprops(GeometryMesh,MeshType) ]} {
          # wa "Checando $wprops(GeometryMesh,MeshType)"
          if { $wprops(ExperimentType,ExperimentType) eq "BTS" } {
               foreach i [list 1 2 3] {
                    if { $wprops(GeometryMesh,MeshType) == $i } {
                         $w.b$i configure -image $img1
                    }
               }
          } else {
               foreach i [list 1 2 3 4 5 6] {
                    if { $wprops(GeometryMesh,MeshType) == $i } {
                         $w.b$i configure -image $img2
                    }
               }
          }
     }
}

# Change button images step 1
proc ::kwiz::buttonPushedStep1 {cpath} {
    variable wprops
    variable wizardid
    variable bwinpath

     set img1 [::WinUtils::GetImage MaterialTest2.gif]
     set img2 [::WinUtils::GetImage CuttingTest2.gif]
     set img3 [::WinUtils::GetImage MaterialTestSel.gif]
     set img4 [::WinUtils::GetImage CuttingTestSel.gif]

     if { $wprops(TestType,TestType) == 1 } {

         .gid.wtest.w buttonstate next normal
         $cpath.b1 configure -image $img3
         $cpath.b2 configure -image $img2
          set wizardid "compressiontest"
     } elseif { $wprops(TestType,TestType) == 2 } {
         .gid.wtest.w buttonstate next normal
         $cpath.b1 configure -image $img1
         $cpath.b2 configure -image $img4
         set wizardid "cutting"
     } else {
         $cpath.b1 configure -image $img1
         $cpath.b2 configure -image $img2
         .gid.wtest.w buttonstate next disabled
     }
     if { $wprops(TestType,TestType) == 0 } {
         $cpath.b1 configure -image $img1
         $cpath.b2 configure -image $img2
         .gid.wtest.w buttonstate next disabled
     }

}



proc ::kwiz::buttonPushedStep2 {cpath id} {
     variable wprops
     global dir

    variable editedMesh
    set editedMesh 1

	 # selected images
     set imgUCSs [::WinUtils::GetImage UCSSelected.gif]
     set imgTris [::WinUtils::GetImage TriSelected.gif]
     set imgUSCs [::WinUtils::GetImage USCSelected.gif]
     set imgHyds [::WinUtils::GetImage TriSelected.gif]
     set imgBTSs [::WinUtils::GetImage BTSSelected.gif]

	 # marked images
     set imgUCSm [::WinUtils::GetImage UCSMarked.gif]
     set imgTrim [::WinUtils::GetImage TriMarked.gif]
     set imgUSCm [::WinUtils::GetImage USCMarked.gif]
     set imgHydm [::WinUtils::GetImage TriMarked.gif]
     set imgBTSm [::WinUtils::GetImage BTSMarked.gif]

        # Si no hay boton apretado deselecciono todas las imagenes
     $cpath.b1 configure -image $imgUCSs
     $cpath.b2 configure -image $imgTris
     $cpath.b3 configure -image $imgUSCs
     $cpath.b4 configure -image $imgHyds
     $cpath.b5 configure -image $imgBTSs
     .gid.wtest.w buttonstate next disabled

     if { $wprops(ExperimentType,ExperimentType) ne 0 } {
     # Si hay boton apretado selecciono la imagen correspondiente

          .gid.wtest.w buttonstate next active
          switch $id {
               1 {
                    $cpath.b1 configure -image $imgUCSm
               }
               2 {
                    $cpath.b2 configure -image $imgTrim
               }
               3 {
                    $cpath.b3 configure -image $imgUSCm
               }
               4 {
                    $cpath.b4 configure -image $imgHydm
               }
               5 {
                    $cpath.b5 configure -image $imgBTSm
               }
          }
     }
}

# Change button images step 2
proc ::kwiz::buttonPushedStep2b {cpath id} {
     variable wprops
     global dir
     variable editedMesh
    set editedMesh 1

     set img1 [::WinUtils::GetImage ExpSelected.gif]
     set img2 [::WinUtils::GetImage ExpMarked.gif]


        if { $wprops(ExperimentType,ExperimentType) == 0 } {
            # Si no hay boton apretado deselecciono todas las imagenes
            .gid.wtest.w buttonstate next disabled
            foreach i [list 1 2 3 4 5] {
            $cpath.b$i configure -image $img1
            }
        } else {
            # Si hay boton apretado selecciono la imagen correspondiente
            foreach i [list 1 2 3 4 5] {

                if { $i == $id } {
                    .gid.wtest.w buttonstate next normal
                    $cpath.b$i configure -image $img2
                } else {
                    $cpath.b$i configure -image $img1
                }
            }
        }

}
# Change button images step 3
proc ::kwiz::buttonPushedStep3 {cpath} {
     variable wprops
     global dir

     variable editedMesh
     set editedMesh 1

     set img1 [::WinUtils::GetImage bts_meshgrey.gif]
     set img2 [::WinUtils::GetImage bts_mesh_marked.gif]
     set img3 [::WinUtils::GetImage ucs_mesh_marked.gif]
     set img4 [::WinUtils::GetImage ucs_meshgrey.gif]
     if { $wprops(ExperimentType,ExperimentType) eq "BTS" } {

        foreach i [list 1  3] {
            if { $wprops(GeometryMesh,MeshType) == $i } {
                .gid.wtest.w buttonstate next normal
                $cpath.b$i configure -image $img2
            } else {
                $cpath.b$i configure -image $img1
            }
        }
        if { $wprops(GeometryMesh,MeshType) ==0 } {
            .gid.wtest.w buttonstate next disabled
        }
        #$cpath.b2 configure -state disabled
     } else {

        foreach i [list 1 2 3 4 5 6] {
            if { $wprops(GeometryMesh,MeshType) == $i } {
                .gid.wtest.w buttonstate next normal
                $cpath.b$i configure -image $img3
            } else {
                $cpath.b$i configure -image $img4
            }
        }
        if { $wprops(GeometryMesh,MeshType) ==0 } {
            .gid.wtest.w buttonstate next disabled
        }
     }
}

# Set material properties default values
proc ::kwiz::initmaterialdefaultvalues {} {
    variable wprops
   set plist [list Density YoungModulus PoissonRatio ParticleFrictionAngle LCS1 LCS2 LCS3 YRC1 YRC2 YRC3 PlasticYoungModulus PlasticYieldStress \
              DamageDeformationFactor TangentialStrength NormalTensileStrength InternalFrictionAngleCoeff]
   set pdflist [list $::defvals::dem::wizmaterial::Density $::defvals::dem::wizmaterial::YoungModulus \
                $::defvals::dem::wizmaterial::PoissonRatio $::defvals::dem::wizmaterial::ParticleFrictionAngle \
                $::defvals::dem::wizmaterial::LCS1 $::defvals::dem::wizmaterial::LCS2 $::defvals::dem::wizmaterial::LCS3 \
                $::defvals::dem::wizmaterial::YRC1 $::defvals::dem::wizmaterial::YRC2 $::defvals::dem::wizmaterial::YRC3 \
                $::defvals::dem::wizmaterial::PlasticYoungModulus $::defvals::dem::wizmaterial::PlasticYieldStress \
                $::defvals::dem::wizmaterial::DamageDeformationFactor $::defvals::dem::wizmaterial::TangentialStrength \
                $::defvals::dem::wizmaterial::NormalTensileStrength $::defvals::dem::wizmaterial::InternalFrictionAngleCoeff]

    foreach pid $plist dfv $pdflist {

        if {[info exists wprops(MaterialParameters,$pid)]} {
            if {$wprops(MaterialParameters,$pid) ==""} {
                set wprops(MaterialParameters,$pid) $dfv
            }
        } else {
                set wprops(MaterialParameters,$pid) $dfv
        }
   }
}

# Combo DeltaTime GeneralSettings
proc ::kwiz::ParallelTypeSelection {w combo} {
     variable wprops
     set ParTypeSel [$combo get]

    if { $ParTypeSel == "OpenMP" } {
         set ::kwiz::wprops(RunProblem,SelectOMPMPI) "OpenMP"
    } else {
         set ::kwiz::wprops(RunProblem,SelectOMPMPI) "MPI"
    }
}

proc ::kwiz::bodystep1 { win } {
    variable wprops
    variable bwinpath

     # Set the widgets
     set img3 [::WinUtils::GetImage newprojectgrey.gif]
     set img4 [::WinUtils::GetImage loadprojectgrey.gif]
     set but3 [ttk::button $win.b3 -text "New project" -image $img3  -command "::kwiz::wizardAction $win 1"]
     set but4 [ttk::button $win.b4 -text "Load project"  -image $img4 -command "::kwiz::wizardAction $win 2"]
     .gid.wtest.w buttonstate next disabled
     # wa [::KUtils::GetPaths "PName"]
     if {[::KUtils::GetPaths "PName"] ne ""} {
          #wa "si"
          # Si tengo nombre, next
          after 1000 [.gid.wtest.w buttonstate next normal]
     } else {
          #wa "no"
          # si no tengo nombre, no next
          after 1000 [.gid.wtest.w buttonstate next disabled]
     }
     # Set the grid
     grid $but3 -column 1 -row 1
     grid $but4 -column 3 -row 1

     # Configure the grid
     grid columnconfigure $win 0 -minsize 50
     grid columnconfigure $win 1 -pad 50 -minsize 220
     grid columnconfigure $win 2 -minsize 20
     grid columnconfigure $win 3 -pad 0 -minsize 220
     grid rowconfigure $win 0 -minsize 38
}

proc ::kwiz::bodystep1bis { win } {
    variable wprops

     catch {destroy $win.b1}
     catch {destroy $win.b2}
     # Set the widgets
     set img3 [::WinUtils::GetImage MaterialTestSel.gif]
     set img4 [::WinUtils::GetImage CuttingTestSel.gif]

     set img1 [::WinUtils::GetImage MaterialTest2.gif]
     set img2 [::WinUtils::GetImage CuttingTest2.gif]
     set but1 [ttk::checkbutton $win.b1 -style TMenubutton.Toolbutton -image $img1 -command "::kwiz::buttonPushedStep1 $win" -variable ::kwiz::wprops(TestType,TestType) -onvalue 1]
     set but2 [ttk::checkbutton $win.b2 -style TMenubutton.Toolbutton -image $img2 -command "::kwiz::buttonPushedStep1 $win" -variable ::kwiz::wprops(TestType,TestType) -onvalue 2]

     # Check image buttons
     ::kwiz::checkImageStep1 $win


     set lab2 [ttk::label $win.l -text "Select test:"]
     grid $lab2 -column 1 -row 0 -sticky w

     # Set the grid
     grid $but1 -column 1 -row 2
     grid $but2 -column 3 -row 2

     # Configure the grid
     grid columnconfigure $win 0 -minsize 77
     grid columnconfigure $win 2 -minsize 50
     grid rowconfigure $win 1 -pad 40
     grid rowconfigure $win 2 -pad 40
}

proc ::kwiz::bodystep2 { win } {
     variable wprops
    if {$wprops(TestType,TestType) == 1} {
      # MATERIAL TEST

     # Set the widgets

     set imgUCSm [::WinUtils::GetImage UCSSelected.gif]
     set imgTrim [::WinUtils::GetImage TriSelected.gif]
     set imgUSCm [::WinUtils::GetImage USCSelected.gif]
     set imgHydm [::WinUtils::GetImage TriSelected.gif]
     set imgBTSm [::WinUtils::GetImage BTSSelected.gif]

     set lab [ttk::label $win.l -text "Select experiment type:"]
     set lab1 [ttk::label $win.l1 -text "Uniaxial compression strength" -wraplength 150]
     set lab2 [ttk::label $win.l2 -text "Triaxial" -wraplength 150]
     set lab3 [ttk::label $win.l3 -text "Uniaxial strain compaction" -wraplength 150]
     set lab4 [ttk::label $win.l4 -text "Hydrostatic" -wraplength 150]
     set lab5 [ttk::label $win.l5 -text "BTS" -wraplength 150]
     set but1 [ttk::checkbutton $win.b1 -style TMenubutton.Toolbutton -image $imgUCSm -command "::kwiz::buttonPushedStep2 $win 1" -variable ::kwiz::wprops(ExperimentType,ExperimentType) -onvalue "UCS"]
     set but2 [ttk::checkbutton $win.b2 -style TMenubutton.Toolbutton -image $imgTrim -command "::kwiz::buttonPushedStep2 $win 2" -variable ::kwiz::wprops(ExperimentType,ExperimentType) -onvalue "Triaxial"]
     set but3 [ttk::checkbutton $win.b3 -style TMenubutton.Toolbutton -image $imgUSCm -command "::kwiz::buttonPushedStep2 $win 3" -variable ::kwiz::wprops(ExperimentType,ExperimentType) -onvalue "Oedometric"]
     set but4 [ttk::checkbutton $win.b4 -style TMenubutton.Toolbutton -image $imgHydm -command "::kwiz::buttonPushedStep2 $win 4" -variable ::kwiz::wprops(ExperimentType,ExperimentType) -onvalue "Hydrostatic"]
     set but5 [ttk::checkbutton $win.b5 -style TMenubutton.Toolbutton -image $imgBTSm -command "::kwiz::buttonPushedStep2 $win 5" -variable ::kwiz::wprops(ExperimentType,ExperimentType) -onvalue "BTS"]


    # Set def MeshType
     if {[info exists wprops(GeometryMesh,MeshType) ]} {
          if {$wprops(GeometryMesh,MeshType) eq 0} {
               set wprops(GeometryMesh,MeshType) 1
          }
     }

     # Check image buttons
     ::kwiz::checkImageStep2 $win


     # Set the grid
     grid $lab -column 1 -row 0 -sticky w
     grid $but1 -column 1 -row 2
     grid $but2 -column 2 -row 2
     grid $but3 -column 3 -row 2
     grid $but4 -column 1 -row 4
     grid $but5 -column 2 -row 4
     grid $lab1 -column 1 -row 3
     grid $lab2 -column 2 -row 3
     grid $lab3 -column 3 -row 3
     grid $lab4 -column 1 -row 5
     grid $lab5 -column 2 -row 5

     # Configure the grid
     grid rowconfigure $win 1 -minsize 15
     grid columnconfigure $win 0 -minsize 50
     grid columnconfigure $win 1 -pad 0
     grid columnconfigure $win 2 -pad 80
     grid columnconfigure $win 3 -pad 0

    } else {
      # CUTTING TEST

          set lab2 [ttk::label $win.l1 -text "TODAVIA POR IMPLEMENTAR!"]

          pack $lab2

    }

}

proc ::kwiz::bodystep3 { win } {
     global dir
     variable wprops
     if { $wprops(ExperimentType,ExperimentType) eq "BTS" } {
        # BTS EXPERIMENT

          # Set the widgets
          set img1 [::WinUtils::GetImage bts_meshgrey.gif]
          set img2 [::WinUtils::GetImage bts_mesh_marked.gif]
          set lab0 [ttk::label $win.l0 -text "Select Geometry & Mesh:"]
          set lab1 [ttk::label $win.l1 -text "Cylinder 27000 elements\nD=150 mm, L=90 mm" -wraplength 150]
          #set lab2 [ttk::label $win.l2 -text "Cylinder 40000 elements\nD=37.8 mm, L=21.9 mm" -wraplength 150]
          set lab3 [ttk::label $win.l3 -text "Cylinder 16000 elements\nD=37.8 mm, L=21.9 mm" -wraplength 150]
          set labspace [ttk::label $win.ls -text ""]
          set but1 [ttk::checkbutton $win.b1 -style TMenubutton.Toolbutton -image $img1 -command "::kwiz::buttonPushedStep3 $win" -variable ::kwiz::wprops(GeometryMesh,MeshType) -onvalue 1]
          #set but2 [ttk::checkbutton $win.b2 -style TMenubutton.Toolbutton -image $img1 -command "::kwiz::buttonPushedStep3 $win" -variable ::kwiz::wprops(GeometryMesh,MeshType) -onvalue 2]
          set but3 [ttk::checkbutton $win.b3 -style TMenubutton.Toolbutton -image $img1 -command "::kwiz::buttonPushedStep3 $win" -variable ::kwiz::wprops(GeometryMesh,MeshType) -onvalue 3]

          ## Set def MeshType
          # if {$wprops(ExperimentType,ExperimentType) eq "BTS"} {
          #     set ::kwiz::wprops(GeometryMesh,MeshType) 1
          # }

         # Check image buttons
          ::kwiz::checkImageStep3 $win

          # Set the grid
          grid $lab0 -column 1 -row 0 -sticky w
          grid $but1 -column 1 -row 2
          #grid $but2 -column 2 -row 2
          grid $but3 -column 3 -row 2
          grid $lab1 -column 1 -row 3
          #grid $lab2 -column 2 -row 3
          grid $lab3 -column 3 -row 3

          # Grid configure
          grid rowconfigure $win 1 -minsize 15
          grid columnconfigure $win 0 -minsize 50
          grid columnconfigure $win 1 -pad 0
          grid columnconfigure $win 3 -pad 80
          grid columnconfigure $win 2 -pad 0

    } else {
        # OTHER EXPERIMENTS

     # Set the widgets
     set img1 [::WinUtils::GetImage ucs_meshgrey.gif]
     set img2 [::WinUtils::GetImage ucs_mesh_marked.gif]
     set lab0 [ttk::label $win.l0 -text "Select Geometry & Mesh:"]
     set lab1 [ttk::label $win.l1 -text "Cylinder 13000 elements\nD=150 mm, L=300 mm" -wraplength 150]
     set lab2 [ttk::label $win.l2 -text "Cylinder 36000 elements\nD=150 mm, L=300 mm" -wraplength 150]
     set lab3 [ttk::label $win.l3 -text "Cylinder 71000 elements\nD=150 mm, L=300 mm" -wraplength 150]
     set lab4 [ttk::label $win.l4 -text "Cylinder 18000 elements\nD=25.4 mm, L=48.26 mm" -wraplength 150]
     set lab5 [ttk::label $win.l5 -text "Cylinder 42000 elements\nD=25.4 mm, L=48.26 mm" -wraplength 150]
     set lab6 [ttk::label $win.l6 -text "Cylinder 65000 elements\nD=25.4 mm, L=48.26 mm" -wraplength 150]
     set labspace [ttk::label $win.ls -text ""]
     set but1 [ttk::checkbutton $win.b1 -style TMenubutton.Toolbutton -image $img1 -command "::kwiz::buttonPushedStep3 $win" -variable ::kwiz::wprops(GeometryMesh,MeshType) -onvalue 1]
     set but2 [ttk::checkbutton $win.b2 -style TMenubutton.Toolbutton -image $img1 -command "::kwiz::buttonPushedStep3 $win" -variable ::kwiz::wprops(GeometryMesh,MeshType) -onvalue 2]
     set but3 [ttk::checkbutton $win.b3 -style TMenubutton.Toolbutton -image $img1 -command "::kwiz::buttonPushedStep3 $win" -variable ::kwiz::wprops(GeometryMesh,MeshType) -onvalue 3]
     set but4 [ttk::checkbutton $win.b4 -style TMenubutton.Toolbutton -image $img1 -command "::kwiz::buttonPushedStep3 $win" -variable ::kwiz::wprops(GeometryMesh,MeshType) -onvalue 4]
     set but5 [ttk::checkbutton $win.b5 -style TMenubutton.Toolbutton -image $img1 -command "::kwiz::buttonPushedStep3 $win" -variable ::kwiz::wprops(GeometryMesh,MeshType) -onvalue 5]
     set but6 [ttk::checkbutton $win.b6 -style TMenubutton.Toolbutton -image $img1 -command "::kwiz::buttonPushedStep3 $win" -variable ::kwiz::wprops(GeometryMesh,MeshType) -onvalue 6]

     # Check image buttons
     ::kwiz::checkImageStep3 $win

     # Set the grid
     grid $lab0 -column 1 -row 0 -sticky w
     grid $but1 -column 1 -row 2
     grid $but2 -column 2 -row 2
     grid $but3 -column 3 -row 2
     grid $lab1 -column 1 -row 3
     grid $lab2 -column 2 -row 3
     grid $lab3 -column 3 -row 3
     grid $but4 -column 1 -row 4
     grid $but5 -column 2 -row 4
     grid $but6 -column 3 -row 4
     grid $lab4 -column 1 -row 5
     grid $lab5 -column 2 -row 5
     grid $lab6 -column 3 -row 5

     # Grid configure
     grid rowconfigure $win 1 -minsize 15
     grid columnconfigure $win 0 -minsize 50
     grid columnconfigure $win 1 -pad 0
     grid columnconfigure $win 2 -pad 80
     grid columnconfigure $win 3 -pad 0
    }
}

proc ::kwiz::ImportExportMaterials {{what Import}} {
     # ABSTRACT: Import/export dempack wizard material file
     variable bwinpath

     set txt1 [= "Dempack Wizard material files"]; set l1 [list "$txt1" ".kwm"]
     if {$what eq "Import"} {
          # Import
          set mfilename [Browser-ramR file read $bwinpath [= "Import material"] "" [list $l1] ".kwm" 0]
          # wa "mfilename:$mfilename"
          set cext "kwm"
          if {$mfilename != "" }  {
               set currext [string tolower [file extension $mfilename]]
               if {$currext != ".kwm" } {
                   set mfilename "${mfilename}.kwm"
               }
          } else {
               set txt [= "First select the Dempack wizard material file to be imported"]
               WarnWin "${txt}."
          }
     } elseif {$what eq "Export"} {
          # Export
          set mfilename [Browser-ramR file write $bwinpath [= "Export material"] "" [list $l1] ".kwm" 0]
          # wa "mfilename:$mfilename"
          set cext "kwm"
          if {$mfilename != "" }  {
               set currext [string tolower [file extension $mfilename]]
               if {$currext != ".kwm" } {
                   set mfilename "${mfilename}.kwm"
               }
          } else {
               set txt [= "First select the Dempack wizard material file to be exported"]
               WarnWin "${txt}."
          }
     }
     return $mfilename
}

proc ::kwiz::ExportMaterial { win } {
     # ABSTRACT: Export a material file
     global KPriv
     variable wprops

     variable wizardid
     set wizid $wizardid

     # Set the internal step identifier
     set stepId "MaterialParameters"

     # Save wizard data
     ::kfiles::SaveWizardDataStep $stepId

     # Get the file path
     set filename [::kwiz::ImportExportMaterials "Export"]

     if {$filename !=""} {
          set xmlData $KPriv(xmlWiz)
          set xmlRoot "/Wizard\[@wizardid='$wizid'\]"

          set path "/${xmlRoot}/Steps"
          set stepNode [$xmlData selectNodes "$path/Step\[@id='$stepId'\]"]
          # wa [$stepNode getAttribute id ""]
          set dataNodes [$xmlData selectNodes "$path/Step\[@id='$stepId'\]/Data"]

          foreach dataNode $dataNodes {
               # wa "stepId:$stepId"
               foreach node [$dataNode childNodes] {

                    # For nodes with no children
                    if {([$node nodeName] eq "Item") && (![$node hasChildNodes])} {
                         set idv [$node getAttribute idv ""]
                         # Assign the current variable value to the xml node
                         set dv [$node getAttribute dv ""]
                         # wwt "node:$node idv:$idv dv:$dv"

                         set path "/Kratos_Wiz_Mat"
                         set parent [$KPriv(xmlWizMat) selectNodes $path]
                         set material [$parent childNodes]
                         foreach node [$material childNodes] {
                              if {[$node getAttribute idv ""] eq $idv} {
                                   $node setAttribute dv $dv
                                   # wa "id: $idv dv: $dv"
                              }
                         }
                    }
               }
          }

          # Write the xml file
          ::xmlutils::writeFile "${filename}" $KPriv(dir) $KPriv(encrXmlWizMat) $KPriv(xmlDocWizMat) $KPriv(RDConfig) 0
     }
}


proc ::kwiz::ImportMaterial { win } {
     # ABSTRACT: Import a material file
     global KPriv
     variable wprops
     variable wizardid

     set wizid $wizardid

     # Set the step identifier
     set stepId "MaterialParameters"

     # Get the imported file
     set filename [::kwiz::ImportExportMaterials "Import"]
     if {$filename !=""} {
          # Open xml file
          set xmlArray [::xmlutils::openFile "." "$filename"]

          set KPriv(xmlWizMat) [lindex $xmlArray 0]
          set KPriv(encrXmlWizMat) [lindex $xmlArray 1]
          set KPriv(xmlDocWizMat) [lindex $xmlArray 2]

          set xmlData $KPriv(xmlWiz)
          set xmlRoot "/Wizard\[@wizardid='$wizid'\]"

          set path "/${xmlRoot}/Steps"
          set stepNode [$xmlData selectNodes "$path/Step\[@id='$stepId'\]"]
           #wa [$stepNode getAttribute id ""]
          set dataNodes [$xmlData selectNodes "$path/Step\[@id='$stepId'\]/Data"]

          foreach dataNode $dataNodes {
               # wa "stepId:$stepId"
               foreach node [$dataNode childNodes] {

                    # For nodes with no children
                    if {([$node nodeName] eq "Item") && (![$node hasChildNodes])} {
                         set idv [$node getAttribute idv ""]
                         # Assign the current variable value to the xml node

                         set path "/Kratos_Wiz_Mat"
                         set parent [$KPriv(xmlWizMat) selectNodes $path]
                         set material [$parent childNodes]
                         foreach node2 [$material childNodes] {
                              if {[$node2 getAttribute idv ""] eq $idv} {
                                   set dv [$node2 getAttribute dv ""]
                                   $node setAttribute dv $dv
                                   # wa "id: $idv dv: $dv"
                              }
                         }
                    }
               }
          }

          # Import the wizard data
          ::kfiles::ImportWizardDataStep $stepId
     }
}

proc ::kwiz::bodystep4 { win } {
    variable wprops

     # Units
     set mpatxt [= "MPa"]
     set mpatxt "\[${mpatxt}\]"
     set patxt [= "Pa"]
     set patxt "\[${patxt}\]"
     # set ulesstxt "\[-\]"
     set ulesstxt ""
     set kgmtxt [= "kg/m^3"]
     set kgmtxt "\[${kgmtxt}\]"
     set dynfrang [= "degrees"]
     set dynfrang "\[${dynfrang}\]"
     set entrywidth 8

     # Set the widgets
     set img1 [::WinUtils::GetImage MaterialPar1b.gif]

     # Left frame
     set fr1 [ttk::frame $win.fr1 -borderwidth 0]
     # Rigth frame
     set fr2 [ttk::frame $win.fr2 -borderwidth 0]

     # Import/export button
     set txt [= "Export material"]
     set but1 [ttk::button $fr1.b1 -text "${txt}..." -command [list ::kwiz::ExportMaterial $win]]
     set txt [= "Clic this button to export the current\nmaterial properties to an external file"]
     tooltip::tooltip $fr1.b1 "${txt}"

     set txt [= "Import material"]
     set but2 [ttk::button $fr1.b2 -text "${txt}..." -command [list ::kwiz::ImportMaterial $win]]
     set txt [= "Clic this button to import material\nproperties from an external file"]
     tooltip::tooltip $fr1.b2 "${txt}"

     # Labels
     set txt [= "Define material parameters"]
     set lab0 [ttk::label $fr1.l0 -text "${txt}:"]

     set labIm [ttk::label $fr2.lIm -image $img1]

     # Material frame
     set labfr1 [ttk::labelframe $fr1.lfr1 -text [= "Material"] -padding 10 -width 200 -height 200 ]

     # Density
     set txt [= "Density"]
     set lab1 [ttk::label $fr1.lfr1.l1 -text "${txt}:"]
     set ent1 [ttk::entry $fr1.lfr1.e1 -textvariable ::kwiz::wprops(MaterialParameters,Density) -width $entrywidth]
     wcb::callback $fr1.lfr1.e1 before insert wcb::checkEntryForReal
     set lab17 [ttk::label $fr1.lfr1.l5 -text "$kgmtxt"]
     set txt [= "Enter a real value for the density"]
     tooltip::tooltip $fr1.lfr1.e1 "${txt}."

     # Young modulus
     set txt [= "Young's modulus"]
     set lab2 [ttk::label $fr1.lfr1.l2 -text "${txt}:"]
     set ent2 [ttk::entry $fr1.lfr1.e2 -textvariable ::kwiz::wprops(MaterialParameters,YoungModulus) -width $entrywidth]
     wcb::callback $fr1.lfr1.e2 before insert wcb::checkEntryForReal
     set lab18 [ttk::label $fr1.lfr1.l6 -text "$patxt"]
     set txt [= "Enter a real value for the Young's modulus"]
     tooltip::tooltip $fr1.lfr1.e2 "${txt}."

     # Poisson ratio
     set txt [= "Poisson's ratio"]
     set lab3 [ttk::label $fr1.lfr1.l3 -text "${txt}:"]
     set ent3 [ttk::entry $fr1.lfr1.e3 -textvariable ::kwiz::wprops(MaterialParameters,PoissonRatio) -width $entrywidth]
     wcb::callback $fr1.lfr1.e3 before insert wcb::checkEntryForReal
     set prl [ttk::label $fr1.lfr1.ul3 -text "$ulesstxt"]
     set txt [= "Enter a real value for the Poisson's ratio"]
     tooltip::tooltip $fr1.lfr1.e3 "${txt}."

     # Dynamic friction coeff
     set txt [= "Dynamic friction angle"]
     set lab4 [ttk::label $fr1.lfr1.l4 -text "${txt}:"]
     set ent4 [ttk::entry $fr1.lfr1.e4 -textvariable ::kwiz::wprops(MaterialParameters,ParticleFrictionAngle) -width $entrywidth]
     wcb::callback $fr1.lfr1.e4 before insert wcb::checkEntryForReal
     set dfrl [ttk::label $fr1.lfr1.ul4 -text "$dynfrang"]
     set txt [= "Enter a real value for the dynamic friction angle (broken material)"]
     tooltip::tooltip $fr1.lfr1.e4 "${txt}."

     # Non-linear elasticity frame
     set labfr2 [ttk::labelframe $fr1.lfr2 -text [= "Non-linear elasticity"] -padding 10 -width 200 -height 100]

     # LCS1
     set txt [= "LCS1"]
     set lab5 [ttk::label $fr1.lfr2.l1 -text "${txt}:"]
     set ent5 [ttk::entry $fr1.lfr2.e1 -textvariable ::kwiz::wprops(MaterialParameters,LCS1) -width $entrywidth]
     wcb::callback $fr1.lfr2.e1 before insert wcb::checkEntryForReal
     set lab19 [ttk::label $fr1.lfr2.l7 -text "$mpatxt"]
     set txt [= "Enter a real value for the LCS1 parameter"]
     tooltip::tooltip $fr1.lfr2.e1 "${txt}."

     # LCS2
     set txt [= "LCS2"]
     set lab6 [ttk::label $fr1.lfr2.l2 -text "${txt}:"]
     set ent6 [ttk::entry $fr1.lfr2.e2 -textvariable ::kwiz::wprops(MaterialParameters,LCS2) -width $entrywidth]
     wcb::callback $fr1.lfr2.e2 before insert wcb::checkEntryForReal
     set lab20 [ttk::label $fr1.lfr2.l8 -text "$mpatxt"]
     set txt [= "Enter a real value for the LCS2 parameter"]
     tooltip::tooltip $fr1.lfr2.e2 "${txt}."

     # LCS3
     set txt [= "LCS3"]
     set lab7 [ttk::label $fr1.lfr2.l3 -text "${txt}:"]
     set ent7 [ttk::entry $fr1.lfr2.e3 -textvariable ::kwiz::wprops(MaterialParameters,LCS3) -width $entrywidth]
     wcb::callback $fr1.lfr2.e3 before insert wcb::checkEntryForReal
     set lab21 [ttk::label $fr1.lfr2.l9 -text "$mpatxt"]
     set txt [= "Enter a real value for the LCS3 parameter"]
     tooltip::tooltip $fr1.lfr2.e3 "${txt}."

     # YRC1
     set txt [= "YRC1"]
     set lab8 [ttk::label $fr1.lfr2.l4 -text "${txt}:"]
     set ent8 [ttk::entry $fr1.lfr2.e4 -textvariable ::kwiz::wprops(MaterialParameters,YRC1) -width $entrywidth]
     wcb::callback $fr1.lfr2.e4 before insert wcb::checkEntryForReal
     set ul8 [ttk::label $fr1.lfr2.ul4 -text "$ulesstxt"]
     set txt [= "Enter a real value for the YRC1 parameter"]
     tooltip::tooltip $fr1.lfr2.e4 "${txt}."

     # YRC2
     set txt [= "YRC2"]
     set lab9 [ttk::label $fr1.lfr2.l5 -text "${txt}:"]
     set ent9 [ttk::entry $fr1.lfr2.e5 -textvariable ::kwiz::wprops(MaterialParameters,YRC2) -width $entrywidth]
     wcb::callback $fr1.lfr2.e5 before insert wcb::checkEntryForReal
     set ul9 [ttk::label $fr1.lfr2.ul5 -text "$ulesstxt"]
     set txt [= "Enter a real value for the YRC2 parameter"]
     tooltip::tooltip $fr1.lfr2.e5 "${txt}."

     # YRC3
     set txt [= "YRC3"]
     set lab10 [ttk::label $fr1.lfr2.l6 -text "${txt}:"]
     set ent10 [ttk::entry $fr1.lfr2.e6 -textvariable ::kwiz::wprops(MaterialParameters,YRC3) -width $entrywidth]
     wcb::callback $fr1.lfr2.e6 before insert wcb::checkEntryForReal
     set ul10 [ttk::label $fr1.lfr2.ul6 -text "$ulesstxt"]
     set txt [= "Enter a real value for the YRC3 parameter"]
     tooltip::tooltip $fr1.lfr2.e6 "${txt}."

     # 1D plasticity and damage
     set labfr3 [ttk::labelframe $fr1.lfr3 -text [= "1D plasticity and damage"] -padding 10 -width 200 -height 100]

     # Unit labels
     set lab23 [ttk::label $fr1.lfr3.l5 -text "$mpatxt"]

     # Plastic young's modulus
     set txt [= "Plastic Young's modulus"]
     set lab11 [ttk::label $fr1.lfr3.l1 -text "${txt}:"]
     set ent11 [ttk::entry $fr1.lfr3.e1 -textvariable ::kwiz::wprops(MaterialParameters,PlasticYoungModulus) -width $entrywidth]
     wcb::callback $fr1.lfr3.e1 before insert wcb::checkEntryForReal
     set lab22 [ttk::label $fr1.lfr3.l4 -text "$patxt"]
     set txt [= "Enter a real value for the plastic Young's modulus"]
     tooltip::tooltip $fr1.lfr3.e1 "${txt}."

     # Plastic yield stress
     set txt [= "Plastic yield stress"]
     set lab12 [ttk::label $fr1.lfr3.l2 -text "${txt}:"]
     set ent12 [ttk::entry $fr1.lfr3.e2 -textvariable ::kwiz::wprops(MaterialParameters,PlasticYieldStress) -width $entrywidth]
     wcb::callback $fr1.lfr3.e2 before insert wcb::checkEntryForReal
     set txt [= "Enter a real value for the plastic yield stress"]
     tooltip::tooltip $fr1.lfr3.e2 "${txt}."

     # Damage deformation factor
     set txt [= "Damage deformation factor"]
     set lab13 [ttk::label $fr1.lfr3.l3 -text "${txt}:"]
     set ent13 [ttk::entry $fr1.lfr3.e3 -textvariable ::kwiz::wprops(MaterialParameters,DamageDeformationFactor) -width $entrywidth]
     wcb::callback $fr1.lfr3.e3 before insert wcb::checkEntryForReal
     set ul13 [ttk::label $fr1.lfr3.ul3 -text "$ulesstxt"]
     set txt [= "Enter a real value for the damage deformation factor"]
     tooltip::tooltip $fr1.lfr3.e3 "${txt}."

     # Failure criterion frame
     set labfr4 [ttk::labelframe $fr1.lfr4 -text [= "Failure criterion"] -padding 10 -width 200 -height 100]

     # Tangential strength
     set txt [= "Tangential strength"]
     set lab14 [ttk::label $fr1.lfr4.l1 -text "${txt}:"]
     set ent14 [ttk::entry $fr1.lfr4.e1 -textvariable ::kwiz::wprops(MaterialParameters,TangentialStrength) -width $entrywidth]
     wcb::callback $fr1.lfr4.e1 before insert wcb::checkEntryForReal
     set lab24 [ttk::label $fr1.lfr4.l4 -text "$mpatxt"]
     set txt [= "Enter a real value for the tangential strength"]
     tooltip::tooltip $fr1.lfr4.e1 "${txt}."

     # Normal tensile strength
     set txt [= "Normal tensile strength"]
     set lab15 [ttk::label $fr1.lfr4.l2 -text "${txt}:"]
     set ent15 [ttk::entry $fr1.lfr4.e2 -textvariable ::kwiz::wprops(MaterialParameters,NormalTensileStrength) -width $entrywidth]
     wcb::callback $fr1.lfr4.e2 before insert wcb::checkEntryForReal
     set lab25 [ttk::label $fr1.lfr4.l5 -text "$mpatxt"]
     set txt [= "Enter a real value for the normal tensile strength"]
     tooltip::tooltip $fr1.lfr4.e2 "${txt}."

     # Internal friction angle coeff
     set txt [= "Internal friction angle coeff"]
     set lab16 [ttk::label $fr1.lfr4.l3 -text "${txt}:"]
     set ent16 [ttk::entry $fr1.lfr4.e3 -textvariable ::kwiz::wprops(MaterialParameters,InternalFrictionAngleCoeff) -width $entrywidth]
     wcb::callback $fr1.lfr4.e3 before insert wcb::checkEntryForReal
     set ul16 [ttk::label $fr1.lfr4.ul3 -text "$ulesstxt"]
     set txt [= "Enter a real value for the internal friction angle coefficient"]
     tooltip::tooltip $fr1.lfr4.e3 "${txt}."

     # Widget geometry
     # Grid the frames
     grid $fr1 -column 1 -row 0 -sticky nw
     grid $fr2 -column 2 -row 0

     # Set the grid
     grid $lab0 -column 0 -row 0 -sticky w
     grid $labIm -column 0 -row 1 -sticky nw -rowspan 3

     # Label frames
     grid $labfr1 -column 0 -row 2 -sticky wen -ipadx 2  -columnspan 2
     grid $labfr2 -column 0 -row 3 -sticky wen -ipadx 2  -columnspan 2
     grid $labfr3 -column 0 -row 4 -sticky wen -ipadx 2 -columnspan 2
     grid $labfr4 -column 0 -row 5 -sticky wen -ipadx 2  -columnspan 2

     # Butons
     grid $but1 -column 0 -row 6 -sticky w
     grid $but2 -column 1 -row 6 -sticky e

     # Label frame 1 => Material frame
     grid $lab1 -column 1 -row 1 -sticky w
     grid $ent1 -column 2 -row 1 -sticky w
     grid $lab17 -column 3 -row 1 -sticky w
     grid $lab2 -column 1 -row 2 -sticky w
     grid $ent2 -column 2 -row 2 -sticky w
     grid $lab18 -column 3 -row 2 -sticky w
     grid $lab3 -column 1 -row 3 -sticky w
     grid $ent3 -column 2 -row 3 -sticky w
     grid $prl -column 3 -row 3 -sticky w
     grid $lab4 -column 1 -row 4 -sticky w
     grid $ent4 -column 2 -row 4 -sticky w
     grid $dfrl -column 3 -row 4 -sticky w

     # Label frame 2 => Non-linear elasticity
     grid $lab5 -column 1 -row 1 -sticky w
     grid $ent5 -column 2 -row 1 -sticky w
     grid $lab19 -column 3 -row 1 -sticky w
     grid $lab6 -column 1 -row 2 -sticky w
     grid $ent6 -column 2 -row 2 -sticky w
     grid $lab20 -column 3 -row 2 -sticky w
     grid $lab7 -column 1 -row 3 -sticky w
     grid $ent7 -column 2 -row 3 -sticky w
     grid $lab21 -column 3 -row 3 -sticky w
     grid $lab8 -column 5 -row 1 -sticky w
     grid $ent8 -column 6 -row 1 -sticky w
     grid $ul8 -column 7 -row 1 -sticky w
     grid $lab9 -column 5 -row 2 -sticky w
     grid $ent9 -column 6 -row 2 -sticky w
     grid $ul9 -column 7 -row 2 -sticky w
     grid $lab10 -column 5 -row 3 -sticky w
     grid $ent10 -column 6 -row 3 -sticky w
     grid $ul10 -column 7 -row 3 -sticky w

     # Label frame 3 => 1D plasticity and damage
     grid $lab11 -column 1 -row 1 -sticky w
     grid $ent11 -column 2 -row 1 -sticky w
     grid $lab22 -column 3 -row 1 -sticky w
     grid $lab12 -column 1 -row 2 -sticky w
     grid $ent12 -column 2 -row 2 -sticky w
     grid $lab23 -column 3 -row 2 -sticky w
     grid $lab13 -column 1 -row 3 -sticky w
     grid $ent13 -column 2 -row 3 -sticky w
     grid $ul13 -column 3 -row 3 -sticky w

     # Label frame 4 => Failure criterion
     grid $lab14 -column 1 -row 1 -sticky w
     grid $ent14 -column 2 -row 1 -sticky w
     grid $lab24 -column 3 -row 1 -sticky w
     grid $lab15 -column 1 -row 2 -sticky w
     grid $ent15 -column 2 -row 2 -sticky w
     grid $lab25 -column 3 -row 2 -sticky w
     grid $lab16 -column 1 -row 3 -sticky w
     grid $ent16 -column 2 -row 3 -sticky w
     grid $ul16 -column 3 -row 3 -sticky w

     # Grid configure
     grid columnconfigure $win 0 -minsize 10
     grid columnconfigure $win 1 -pad 10
     grid rowconfigure $fr1 5 -pad 20
     grid columnconfigure $fr1.lfr2 4 -minsize 18

}

proc ::kwiz::bodystep5 { win } {
    variable wprops

     # Units
     set mpatxt [= "MPa"]
     set mpatxt "\[${mpatxt}\]"
     set sectxt [= "sec"]
     set sectxt "\[${sectxt}\]"
     set mstxt [= "m/s"]
     set mstxt "\[${mstxt}\]"
     # set ulesstxt "\[-\]"
     set ulesstxt ""

     set entrywidth 16

     # Set the widgets
     set img1 [::WinUtils::GetImage dem-wct-settings.gif]
     set txt [= "General settings"]
     set lab0 [ttk::label $win.l0 -text "${txt}:"]

     # Simulation setting frame
     set labfr1 [ttk::labelframe $win.lf1 -text [= "Simulation settings"]]
     set labIm1 [ttk::label $win.lf1.li1 -image $img1]

     # Calculation time
     set txt [= "Calculation time"]
     set lab1 [ttk::label $win.lf1.l1 -text "${txt}:"]
     set ent1 [ttk::entry $win.lf1.e1 -textvariable ::kwiz::wprops(GeneralSettings,CalculationTime) -width $entrywidth]
     wcb::callback $win.lf1.e1 before insert wcb::checkEntryForReal
     set lab7 [ttk::label $win.lf1.l7 -text "${sectxt}"]
     set txt [= "Enter a real value for the calculation time"]
     tooltip::tooltip $win.lf1.e1 "${txt}."

     # Expected deformation
     set txt [= "Expected deformation"]
     set lab3 [ttk::label $win.lf1.l3 -text "${txt}:"]
     set result [::kwiz::ExpDeformation]
     set lab11 [ttk::label $win.lf1.17 -text $result]
     # Calculation time * Loading velocity / Specimen Lenght * 100

     # Delta time selection
     set txt [= "Delta time selection"]
     set lab2 [ttk::label $win.lf1.l2 -text "${txt}:"]
     set combo [ttk::combobox $win.lf1.c -values [list "Custom" "Predefined"] -textvariable ::kwiz::wprops(GeneralSettings,DeltaTimeSelection) \
               -state readonly -width [expr {$entrywidth-3}]]
     set txt [= "Select the delta time calculation method from the combobox"]
     tooltip::tooltip $win.lf1.c "${txt}."

     # Set default value to combo
     if {([info exists ::kwiz::wprops(GeneralSettings,DeltaTimeSelection)]) && ($::kwiz::wprops(GeneralSettings,DeltaTimeSelection) !="")} {
          $combo set $::kwiz::wprops(GeneralSettings,DeltaTimeSelection)
     } else {
          # Select the first value from the combobox
          $combo current 0
     }

     # Loading velocity
     set txt [= "Loading velocity"]
     set lab5 [ttk::label $win.lf1.l5 -text "${txt}:"]
     set ent3 [ttk::entry $win.lf1.e3 -textvariable ::kwiz::wprops(GeneralSettings,LoadingVelocity) -width $entrywidth]
     wcb::callback $win.lf1.e3 before insert wcb::checkEntryForReal
     set lab9 [ttk::label $win.lf1.l9 -text "${mstxt}"]
     set txt [= "Enter a real value for the loading velocity"]
     tooltip::tooltip $win.lf1.e3 "${txt}."

     # Pressure
     set txt [= "Pressure"]
     set lab6 [ttk::label $win.lf1.l6 -text "${txt}:"]
     set ent4 [ttk::entry $win.lf1.e4 -textvariable ::kwiz::wprops(GeneralSettings,ConfinementPressure) -width $entrywidth]
     wcb::callback $win.lf1.e4 before insert wcb::checkEntryForReal
     set lab10 [ttk::label $win.lf1.l10 -text "${mpatxt}"]
     set txt [= "Enter a real value for the pressure"]
     tooltip::tooltip $win.lf1.e4 "${txt}."

     # Set the state for delta time
     set statedts "normal"
     set DTS ""
     if {[info exists ::kwiz::wprops(GeneralSettings,DeltaTimeSelection)]} {
          set DTS $::kwiz::wprops(GeneralSettings,DeltaTimeSelection)
     }
     if {$DTS eq "Predefined" } {
         set statedts "disabled"
         # valores

     }
     # Delta time
     set txt [= "Delta time"]
     set lab4 [ttk::label $win.lf1.l4 -text "${txt}:" -state $statedts]
     set lab8 [ttk::label $win.lf1.l8 -text "${sectxt}" -state $statedts]
     set ent2 [ttk::entry $win.lf1.e2 -textvariable ::kwiz::wprops(GeneralSettings,DeltaTime) -state $statedts -width $entrywidth]
     wcb::callback $win.lf1.e2 before insert wcb::checkEntryForReal
     set txt [= "Enter a real value for delta time"]
     tooltip::tooltip $win.lf1.e2 "${txt}."

     # Set the grid
     grid $lab0 -column 1 -row 0 -sticky w

     # Label frame 1
     grid $labfr1 -column 1 -row 2 -sticky w -ipadx 15 -ipady 8

     # Widgets
     grid $labIm1 -column 0 -row 0 -rowspan 5

     # Calculation time
     grid $lab1 -column 1 -row 0 -sticky w
     grid $ent1 -column 2 -row 0 -sticky w
     grid $lab7 -column 3 -row 0 -sticky w

     # Expected deformation
     grid $lab3 -column 1 -row 1 -sticky w
     grid $lab11 -column 2 -row 1 -sticky w

     # Delta time selection
     grid $combo -column 2 -row 2 -sticky w
     grid $lab2 -column 1 -row 2 -sticky w

     # Loading velocity
     grid $lab5 -column 1 -row 4 -sticky w
     grid $ent3 -column 2 -row 4 -sticky w
     grid $lab9 -column 3 -row 4 -sticky w

     # Pressure
     grid $lab6 -column 1 -row 5 -sticky w
     grid $lab10 -column 3 -row 5 -sticky w
     grid $ent4 -column 2 -row 5 -sticky w

     # Delta time
     grid $lab4 -column 1 -row 3 -sticky w
     grid $lab8 -column 3 -row 3 -sticky w
     grid $ent2 -column 2 -row 3 -sticky w

     # Get the experiment type
     set ExperimentType ""
     if {[info exists wprops(ExperimentType,ExperimentType)]} {
          set ExperimentType $wprops(ExperimentType,ExperimentType)
     }
     set cflag [expr {($ExperimentType eq "Triaxial") || ($ExperimentType eq "Hydrostatic")}]
     if {!$cflag} {
          set wdlist [list $lab6 $lab10 $ent4]
          foreach wd $wdlist {
               if {[winfo exists $wd]} {
                    grid remove $wd
               }
          }
     }

     # Grid configure
     grid rowconfigure $win 1 -minsize 15
     grid columnconfigure $win 0 -minsize 30
     grid columnconfigure $win.lf1 0 -pad 10


     # Dependencias de Expected Deformation
     wcb::callback $win.lf1.e1 after insert ::kwiz::updateExpectedDef
     wcb::callback $ent3 after insert ::kwiz::updateExpectedDef


     # Combobox frame1 bind
     bind $combo <<ComboboxSelected>> [list ::kwiz::DeltaTimeSelection $win %W]
}

proc ::kwiz::updateExpectedDef { w idx str} {
     set labl ".gid.wtest.w.layoutFrame.wiz.layout.lf1.17"
     $labl configure -text [::kwiz::ExpDeformation]
}

proc ::kwiz::ExpDeformation { } {
	# Calculate and update the expected deformation value
	# Calculation time * Loading velocity / Specimen Lenght * 100

	variable wprops
	set len $wprops(GeometryMesh,ProbeLength)
	set CT $wprops(GeneralSettings,CalculationTime)
	set LV $wprops(GeneralSettings,LoadingVelocity)
	set result [format "%.2f" [expr {$CT * $LV / $len * 100}] ]
        append result " %"
	return $result
}

proc ::kwiz::DeltaTimeSelection {w combo} {
     # ABSTRACT: Update the combobox value for the delta time method
     variable wprops

     set ExperimentType $wprops(ExperimentType,ExperimentType)
     set meshtype $wprops(GeometryMesh,MeshType)

     set DeltaTimeSel [$combo get]
     # wa "Deltatime $DeltaTimeSel"
     set wdlist [list $w.lf1.l4 $w.lf1.e2 $w.lf1.l8]
     set cstate "normal"
     if { $DeltaTimeSel eq "Predefined" } {
          set wprops(GeneralSettings,DeltaTime) 0.5e-7
          if {$ExperimentType ne "BTS"} {
               if {$meshtype eq 1} {
                    set wprops(GeneralSettings,DeltaTime) 1.0e-7
               }
          }
          if {$ExperimentType eq "BTS"} {
               if {$meshtype eq 1} {
                    set wprops(GeneralSettings,DeltaTime) 1.0e-7
               }
          }

          set cstate "disabled"

     }
     foreach wd $wdlist {
          if {[winfo exists $wd]} {
               $wd configure -state $cstate
          }
     }
     # wa "DeltaTimeSelection:$wprops(GeneralSettings,DeltaTimeSelection)"
}

proc ::kwiz::bodystep6 { win } {
    variable wprops

     # Units
     set sectxt [= "sec"]
     set sectxt "\[${sectxt}\]"
     # set ulesstxt "\[-\]"
     set ulesstxt ""

     set entrywidth 10

     # Set the widgets
     set txt [= "Output settings"]
     set lab0 [ttk::label $win.l0 -text "${txt}:" ]
     # Nodal results frame
     set labfr1 [ttk::labelframe $win.lf1 -text [= "Nodal results"]]
     set img1 [::WinUtils::GetImage dem-wct-nodalresults.gif]
     set labIm1 [ttk::label $win.lf1.li1 -image $img1]
     set but1 [ttk::checkbutton $win.lf1.b1 -text [= "Displacement"] -offvalue "No" -onvalue "Yes" -variable ::kwiz::wprops(OutputSettings,Displacement)]
     set but2 [ttk::checkbutton $win.lf1.b2 -text [= "Velocity"] -offvalue "No" -onvalue "Yes" -variable ::kwiz::wprops(OutputSettings,Velocity)]
     set but3 [ttk::checkbutton $win.lf1.b3 -text [= "Total forces"] -offvalue "No" -onvalue "Yes" -variable ::kwiz::wprops(OutputSettings,TotalForces)]
     set but4 [ttk::checkbutton $win.lf1.b4 -text [= "Rhs"] -offvalue "No" -onvalue "Yes" -variable ::kwiz::wprops(OutputSettings,Rhs)]
     set but5 [ttk::checkbutton $win.lf1.b5 -text [= "Damp forces"] -offvalue "No" -onvalue "Yes" -variable ::kwiz::wprops(OutputSettings,DampForces)]
     set but6 [ttk::checkbutton $win.lf1.b6 -text [= "Applied forces"] -offvalue "No" -onvalue "Yes" -variable ::kwiz::wprops(OutputSettings,AppliedForces)]

     # Contact mesh results frame
     set labfr2 [ttk::labelframe $win.lf2 -text [= "Contact mesh results"]]
     set img2 [::WinUtils::GetImage dem-wct-contacts.gif]
     set labIm2 [ttk::label $win.lf2.li1 -image $img2]
     set but7 [ttk::checkbutton $win.lf2.b1 -text [= "Contact sigma"] -offvalue "No" -onvalue "Yes" -variable ::kwiz::wprops(OutputSettings,ContactSigma)]
     set but8 [ttk::checkbutton $win.lf2.b2 -text [= "Contact tau"] -offvalue "No" -onvalue "Yes" -variable ::kwiz::wprops(OutputSettings,ContactTau)]
     set but9 [ttk::checkbutton $win.lf2.b3 -text [= "Local contact force"] -offvalue "No" -onvalue "Yes" -variable ::kwiz::wprops(OutputSettings,LocalContactForce)]
     set but10 [ttk::checkbutton $win.lf2.b4 -text [= "Failure criterion state"] -offvalue "No" -onvalue "Yes" -variable ::kwiz::wprops(OutputSettings,FailureCriterionState)]

     # Postprocess definition frame
     set labfr3 [ttk::labelframe $win.lf3 -text [= "Postprocess definition"]]
     set img3 [::WinUtils::GetImage dem-wct-post.gif]
     set labIm3 [ttk::label $win.lf3.li1 -image $img3]
     set but11 [ttk::checkbutton $win.lf3.b1 -text [= "Single postprocess file?"] -offvalue "Multiples" -onvalue "Single" -variable ::kwiz::wprops(OutputSettings,SinglePostprocessFile)]

     # Print output file every
     set txt [= "Print output file every"]
     set lab1 [ttk::label $win.lf3.l1 -text "${txt}:"]
     set ent1 [ttk::entry $win.lf3.e1 -textvariable ::kwiz::wprops(OutputSettings,PrintOutputFile)]
     wcb::callback $win.lf3.e1 before insert wcb::checkEntryForReal
     set lab3 [ttk::label $win.lf3.l3 -text "${sectxt}"]
     set txt [= "Enter a real value for print output file in seconds"]
     tooltip::tooltip $win.lf3.e1 "${txt}."

     # Print graph data every
     set txt [= "Print graph data every"]
     set lab2 [ttk::label $win.lf3.l2 -text "${txt}:"]
     set ent2 [ttk::entry $win.lf3.e2 -textvariable ::kwiz::wprops(OutputSettings,PrintGraphData)]
     wcb::callback $win.lf3.e2 before insert wcb::checkEntryForReal
     set lab4 [ttk::label $win.lf3.l4 -text "${sectxt}"]
     set txt [= "Enter a real value for print graph data in seconds"]
     tooltip::tooltip $win.lf3.e2 "${txt}."


     # Set the grid
     grid $lab0 -column 0 -row 0 -sticky w -padx 30
     grid $labfr1 -row 2 -column 0 -sticky wen -padx 30 -ipadx 62
     grid $labfr2 -row 3 -column 0 -sticky wen -padx 30 -ipadx 43
     grid $labfr3 -row 4 -column 0 -sticky wen -padx 30 -ipadx 10

     # Frame 1
     grid $labIm1 -row 0 -column 0 -rowspan 5 -padx 10 -pady 10
     grid $but1 -row 1 -column 1 -sticky w
     grid $but2 -row 2 -column 1 -sticky w
     grid $but3 -row 3 -column 1 -sticky w
     grid $but4 -row 1 -column 2 -sticky w
     grid $but5 -row 2 -column 2 -sticky w
     grid $but6 -row 3 -column 2 -sticky w

     # Frame 2
     grid $labIm2 -row 0 -column 0 -rowspan 5 -padx 10 -pady 10
     grid $but7 -row 1 -column 1 -sticky w
     grid $but8 -row 2 -column 1 -sticky w
     grid $but9 -row 1 -column 2 -sticky w
     grid $but10 -row 2 -column 2 -sticky w

     # Frame 3
     grid $labIm3 -row 0 -column 0 -rowspan 5 -padx 10 -pady 10
     grid $but11 -row 1 -column 1 -sticky wn -columnspan 3
     grid $lab1 -row 2 -column 1 -sticky wn
     grid $ent1 -row 2 -column 2 -sticky wn
     grid $lab3 -row 2 -column 3 -sticky wn
     grid $lab2 -row 3 -column 1 -sticky wn
     grid $ent2 -row 3 -column 2 -sticky wn
     grid $lab4 -row 3 -column 3 -sticky wn

     # Grid configure
     grid columnconfigure $win 0 -minsize 30
     grid rowconfigure $win 1 -minsize 15
     grid rowconfigure $win 2 -pad 20
     grid rowconfigure $win 3 -pad 20
     grid rowconfigure $win 4 -pad 20
     grid columnconfigure $win.lf1 1 -pad 20
     grid columnconfigure $win.lf2 1 -pad 20
     grid columnconfigure $win.lf3 1 -pad 20
     grid rowconfigure $win.lf1 0 -minsize 10
     grid rowconfigure $win.lf2 0 -minsize 10
     grid rowconfigure $win.lf3 0 -minsize 10

}

proc ::kwiz::TopWindow { } {
    variable bwinpath
    wm attributes $bwinpath -topmost 1
}

proc ::kwiz::GetMeshPath { ExpType cmeshid } {
     global KPriv

     if {$ExpType eq "BTS"} {
          set temp "BTS"
          append temp $cmeshid
          set cmeshid $temp
     }
     set path ""
     set xmlData $KPriv(xmlWiz)
     set xmlRoot "/Wizard/MeshDB"
     set Nodes [$xmlData selectNodes $xmlRoot]
     foreach node [$Nodes childNodes] {
          if {[$node getAttribute experimenttype ""] eq $ExpType} {
               if {[$node getAttribute id ""] eq $cmeshid} {
                    set path [$node getAttribute modelpath ""]
               }
          }
     }
     # wa "meshpath $path"
     return $path
}

proc ::kwiz::blocking {what} {
     if {$what ne "get"} {
          set dir [GidUtils::GetTmp]
          set auxfile [file nativename [file join [GidUtils::GetTmp] "DemAuxFile.txt" ] ]
          set ch [open $auxfile w]
          puts $ch $what
          close $ch
     } else {
          set dir [GidUtils::GetTmp]
          set auxfile [file nativename [file join [GidUtils::GetTmp] "DemAuxFile.txt" ] ]
	 if [file exists $auxfile] {
          set ch [open $auxfile r]
          set bl [gets $ch]
          close $ch
          return $bl
	 } else {
	     return ""
	 }
     }
}

proc ::kwiz::Precalc { win } {
     ::kwiz::PrepareGeom $win
     ::kwiz::PrepareGeom $win
}

proc ::kwiz::PrepareGeom { win } {
     variable wprops
     variable editing

     ::kwiz::blocking "block"

     # ::WinUtils::PrintArray wprops
     ::kwiz::UpdateEvent 5
     ::kwiz::UpdateEvent 7

     # Vamos a buscar la malla activa
     if { $wprops(ExperimentType,ExperimentType) eq "BTS" } {
          set ExpType "BTS"
     } else {
          set ExpType "UCS"
     }
     set activemesh $wprops(GeometryMesh,MeshType)

     set modelpath [file nativename [::kwiz::GetMeshPath $ExpType $activemesh] ]

     # Get the window name
     set windowname [wm title .gid]


     # Get the model full path
     append fullpath [::KUtils::GetPaths "PFPath"] ".spd"

     # Save the model files
     ::kfiles::SaveSPD $fullpath

     # Set filetocopy model.wiz
     set filetocopy [file nativename [::kwiz::SaveActiveWizardtoSPD $fullpath] ]
     # wa "filetocopy $filetocopy"
     set filenamecorrect [lindex [split $filetocopy "."] 0]
     append filenamecorrect ".gid"
     # wa " FIlenamecorrect $filenamecorrect"

     # Temporal dir to save model.wiz
     set dir [GidUtils::GetTmp]
     # wa $dir
     if { [file isdirectory $dir] } {
          if { [file exists $filetocopy] } {
               file copy -force -- $filetocopy $dir
          }

     }

     set filedest [file join $dir [file tail $filetocopy]]
     if { ![file exists $filedest] } {
           return -1
     }

     # Congelamos GiD
      GidUtils::DisableGraphics

     # Conseguimos el path del modelo a cargar
     set cPTpath [::KUtils::GetPaths "PTDir"]
     set meshfull [file join $cPTpath $modelpath]

     #set newfil [file nativename [file join [GidUtils::GetTmp] "temp.gid"] ]
     #GiD_Process Mescape File New
     set ::KPriv(kike_no_unload_problemtype) 1
     # Cargamos el Modelo con los datos de la malla
     GiD_Process Mescape File Read $meshfull escape

     set ::KPriv(kike_no_unload_problemtype) 0
     # Eliminamos el modelo nuestro ?

     # Guardamos con el nombre que toca
     GiD_Process Mescape Files SaveAs -alsoresults:0 -- geoversion:current $filenamecorrect

    # wa "Filename $filenamecorrect"
     # Copiamos el wizard.wiz donde toca
     file copy -force -- $filedest $filenamecorrect

     # Cargamos el .wiz de disco a memoria
     ::kfiles::LoadWizardData "noload"

     # Transferimos el .wiz al ::kwiz::wprops
     ::kfiles::ImportWizardData

     # Propagar del .wiz  al arbol
     ::kfiles::SaveWizardToSPD

     # Save the model
     GiD_Process Mescape Files Save

     # Enable graphic
      GidUtils::EnableGraphics

     # Update GiD GUI
     update
     update idletasks
     if {[GiD_Info Project ViewMode] ne "MESHUSE"} {
          GiD_Process Mescape Meshing MeshView
     }

     GiD_Process 'Zoom Frame
     GiD_Process 'Render SmoothLighting
     GidUtils::UpdateWindow "LAYER"
     GidUtils::UpdateWindow "GROUPS"

     $win.lf2.b1 configure -state normal
     $win.lf2.b2 configure -state disabled

     # Reset edit mesh
     variable editedMesh
     set editedMesh 0

    # Set the window title
    # Get the window name
     wm title .gid $windowname

     ::kwiz::blocking "unlock"
}

proc ::kwiz::RunCalculate { } {

    # Save last step
    ::kwiz::UpdateEvent 7


    GiD_MustRemeshFlag set 0

    # Start the simulation
    GiD_Process Mescape Utilities Calculate

    # Start the check running process event
    ::kwiz::CheckRunningProcess

    #::kwiz::TopWindow
}

proc ::kwiz::CheckRunningProcess { } {
    # ABSTRACT:  We check if there are any project process running

    set binterval 1

    # Wait that the variable ::ccgui::IsProcessRunning is equal to 1 => pass for the event BeforeRunCalculation
    after [expr {$binterval*1000}] "::kwiz::WaitForBeforeRunCalculationEvent"

}

proc ::kwiz::WaitForBeforeRunCalculationEvent {} {
    # ABSTRACT: Wait that the variable ::KUtils::IsProcessRunning is equal to 1 => pass for the event BeforeRunCalculation

    # Check if the current project process are running
    set ProcessIsRunning [::KUtils::ProjectProcessAreRunning]

    # wa "CheckRunningProcess => ProcessIsRunning:$ProcessIsRunning"
    if {$ProcessIsRunning} {
	# Get the process running variable state
	set StartProcessControl [::KUtils::IsProcessRunningVar Get]
	# wa "StartProcessControl:$StartProcessControl"
	if {!$StartProcessControl} {
	    ::kwiz::CheckRunningProcess
	} else {
	    set interval 3
	    set firsttime 1
	    # Try to update the results
	    ::kwiz::CheckUpdateResults $interval $firsttime
	}
    }
}

proc ::kwiz::CheckUpdateResults { {interval 3} {firsttime 0}} {
    # ABSTRACT: Recursive function that is checking the results until the end of the calculation
    # Arguments
    # Interval    => Every few seconds will check the percentage calculated
    # firsttime   => In the first time the procedure is called this variables is equal to 1

     # wa "interval:$interval firsttime:$firsttime ProcessIsRunning:$ProcessIsRunning what:$what"

    # Get the state of the process running using the global variable
    set ProcessIsRunning [::KUtils::IsProcessRunningVar Get]
    # wa "CHECK UPDATE RESULTS and get ProcessIsRunning: $ProcessIsRunning"

    # Access the results file and put the percentage value in ::kwiz::ccalvalue
    # ::kwiz::UpdateLabelPercent $ProcessIsRunning $firsttime

    # Get the percentage
    # set LoadedPercentage $::kwiz::ccalvalue
    # wa "LoadedPercentage:$LoadedPercentage"

    if {$ProcessIsRunning} {

        if {$firsttime} {
          # Automatically go to the next step
          $::kwiz::bwinpath.w gotostep ::kwiz::step8
        }

	# Try to read again
	# after [expr {$interval*1000}] "::kwiz::CheckUpdateResults"

    } else {
	# Finish the simulation
	# ::kwiz::FinishSimulation $what
    }
}

proc ::kwiz::bodystep7 { win } {
     variable wprops

     # Units
     # set ulesstxt "\[-\]"
     set ulesstxt ""

     set entrywidth 16

     # Set the widgets
     # Parallel options
     set labfr1 [ttk::labelframe $win.lf1 -text [= "Parallel options"]]
     set img1 [::WinUtils::GetImage dem-wct-runsimulation.gif]
     set labIm1 [ttk::label $win.lf1.li1 -image $img1]
     set txt [= "Parallel type"]
     set lab1 [ttk::label $win.lf1.l1 -text "${txt}:"]
     set combo [ttk::combobox $win.lf1.c1 -state readonly -values [list OpenMP MPI] -width [expr {$entrywidth-3}]]
     set txt [= "Select the parallelization method from the combobox"]
     tooltip::tooltip $win.lf1.c1 "${txt}."

     # Set default value to combo
     if {[info exists wprops(RunProblem,SelectOMPMPI)]} {
          $combo set $::kwiz::wprops(RunProblem,SelectOMPMPI)
     } else {
          $combo current 0
     }

     set txt [= "Number of threads"]
     set lab2 [ttk::label $win.lf1.l2 -text "${txt}:"]
     set ent1 [ttk::entry $win.lf1.e2 -textvariable ::kwiz::wprops(RunProblem,NumberOfThreads) -width $entrywidth]
     wcb::callback $win.lf1.e2 before insert wcb::checkEntryForInt
     set ul2 [ttk::label $win.lf1.ul2 -text "${ulesstxt}"]
     set txt [= "Enter a integer value for the number of threads"]
     tooltip::tooltip $win.lf1.e2 "${txt}."

     # Run simulation
     set labfr2 [ttk::labelframe $win.lf2 -text [= "Run simulation"]]

     variable editedMesh
     if {$editedMesh} {
          set runstate "disabled"
          set preparestate "enabled"
     } else {
          set runstate "enabled"
          set preparestate "enabled"
     }

     # wa "EditedMesh $editedMesh"
     set but1 [ttk::button $win.lf2.b1 -text [= "Run simulation"] -state $runstate -command [list ::kwiz::RunCalculate] -width 20 ]
     set but2 [ttk::button $win.lf2.b2 -text [= "Prepare Data"] -state $preparestate -command [list ::kwiz::Precalc $win ] -width 20]

     # Set the grid
     grid $labfr1 -column 1 -row 1 -sticky wen -ipadx 20
     grid $labIm1 -column 0 -row 0 -rowspan 2
     grid $but1 -column 1 -row 3
     grid $but2 -column 0 -row 3


     # Frame 1
     grid $labfr2 -column 1 -row 2 -sticky wen -ipadx 20
     grid $lab1 -column 1 -row 0 -sticky w
     grid $lab2 -column 1 -row 1 -sticky w
     grid $combo -column 2 -row 0 -sticky w
     grid $ent1 -column 2 -row 1 -sticky w
     grid $ul2 -column 3 -row 1 -sticky w


     # Configure the grid
     grid columnconfigure $win 0 -minsize 30
     grid rowconfigure $win 2 -minsize 20
     grid columnconfigure $win.lf1 0 -pad 10
     grid columnconfigure $win.lf2 0 -pad 150
     grid rowconfigure $win 1 -pad 15


     # Combobox frame1 bind
     bind $combo <<ComboboxSelected>> [list ::kwiz::ParallelTypeSelection $win %W]
}

proc ::kwiz::getPoints { } {
     variable wprops
     # Get the graph file
     set dir [::KUtils::GetPaths "PFPath"]
     append dir "_Graphs"
     set name [::KUtils::GetPaths "PName"]
     set dir [file nativename [file join $dir $name]]
     # Comprobar si es BTS o no
     if { $wprops(ExperimentType,ExperimentType) eq "BTS" } {

          append dir "_bts.grf"
     } else {
          append dir "_graph.grf"
     }
     # wa $dir
     set lista [list ]
     if {[file exists $dir]} {
          set f [open $dir RDONLY]

          set i 0
          while {![eof $f]} {
               set line [string trim [gets $f] ]
               #wa $line
               if {[string length $line] > 2} {
                    if {$i eq "1"} {
                         lappend lista $line
                         set i 0
                    }
               }
               incr i
          }

          close $f
     }
     # wa $lista
     return $lista
}

proc ::kwiz::PrintGraph {win} {

     set graframe $win.graframe
     set btn $win.lf2.b1

     if {[winfo exists $btn]} {
          $btn configure -state disabled
     }

     # Check to destroy the canvas
     set cpath "$graframe.c"
     if {[winfo exists $cpath]} {
          destroy $cpath
     }

     # Get the curve points
     set pvalues [::kwiz::getPoints]

     set cpath [canvas $graframe.c -relief sunken -borderwidth 2 -width 476 ]
     set title [= "Strain Stress Graph Previsualization"]
     set x_txt [= "Axial strain"]
     set xlabel "${x_txt} \[%\]"
     set y_txt [= "Stress"]
     set ylabel "${y_txt} \[MPa\]"
     #wa $pvalues

     # Plot the graph
     ::KPlot::Plot $cpath $pvalues $title $xlabel $ylabel
     grid $cpath -in $graframe -sticky ew

     if {[winfo exists $btn]} {
          after 500 [list $btn configure -state normal]
     }
}

proc ::kwiz::bodystep8 { win } {
     variable wprops

     # Set the widgets
     set labfr1 [ttk::labelframe $win.lf1 -text [= "Result View"]]
     set labfr2 [ttk::frame $win.lf2]
     set graframe [ttk::frame $win.graframe]
     set logframe [ttk::frame $win.logframe]

     set but1 [ttk::button $win.lf2.b1 -text [= "Refresh Graph"] -command [list ::kwiz::PrintGraph $win] -width 20]
     set but2 [ttk::button $win.lf2.b2 -text [= "Terminate"] -command [list ::kwiz::Terminate] -width 20]
     set but3 [ttk::button $win.lf2.b3 -text [= "More Info"] -command [list PWViewOutput] -width 20]

     set lab1 [ttk::label $win.lf1.l1 -text [= "Select Graph"]]
     set combo [ttk::combobox $win.lf1.c1 -state readonly -values [list Graph1]]
     # Set default value to combo
     $combo set $::kwiz::wprops(ResultView,Graph1)

     catch {destroy $graframe.c}
     # Pintar el grafico
     ::kwiz::PrintGraph $win

     # JG Usar funcion ::KUtils::IsProcessRunningVar para ver el estado

     set log [tk::text $logframe.log -height 4]
     set txt [= "Calculation is running"]
     $log insert 1.0 "$txt\n"
     set txt [= "Please wait"]
     $log insert 2.0 "${txt}...\n"
     # Significa instertar este texto en la linea 2, caracter inicial 0
     $log configure -state disabled

     # Set the grid
     grid $labfr1 -column 1 -row 0 -sticky we -ipadx 20
     grid $labfr2 -column 1 -row 1 -sticky we -ipadx 20

     grid $graframe -column 1 -row 2 -sticky wes
     grid $logframe -column 1 -row 3 -sticky nsew

     # Frame 1
     grid $lab1 -column 1 -row 0 -sticky w
     grid $combo -column 2 -row 0 -sticky w

     # Frame Botones
     grid $but1 -column 1 -row 1
     grid $but2 -column 2 -row 1
     grid $but3 -column 3 -row 1

     # Frame Log
     grid $log -sticky ew

     # Configure the grid
     grid columnconfigure $win 0 -minsize 30
     grid rowconfigure $win 2 -minsize 20
     grid columnconfigure $win.lf1 0 -pad 10
     grid columnconfigure $win.lf2 0 -pad 150
     grid rowconfigure $win 1 -pad 15
      update
     ::kwiz::WhileCalculation $log
}

proc ::kwiz::WhileCalculation { log } {
     set cont 1
     #wa "While Calculation"

     set basexpath "DEM//c.DEM-SolutionStrategy//c.DEM-TimeParameters//i.DEM-ScreenInfoOutput"
     set screenOutTime [::xmlutils::setXml $basexpath "dv"]
     # wa "Screen $screenOutTime"

     if {[winfo exists $log]} {
          set infofile [::KUtils::GetPaths "PFPath"]
          set infofile [append infofile ".info"]
          if {[file exists $infofile]} {
               set f [open $infofile RDONLY]
               catch {seek $f -150 end}

               set lineuno ""
               set linedos ""
               set linetres ""
               set linecuat ""
               while {![eof $f]} {
                    set cont 1
                    update idletasks
                    #wa "Long [string length [gets $f]]"
                    set line [string trim [gets $f] ]

                    set print 0
                    #wa [string first "Real time" $line ]
                    if {[string first "Real time" $line ] ne -1} {
                         if {$lineuno ne $line} {
                              set lineuno $line
                              $log configure -state normal
                              $log delete 1.0
                              $log insert 1.0 "$line\n"
                              $log configure -state disabled
                         }

                    }
                    if {[string first "Simulation" $line ] ne -1} {
                         if {$linedos ne $line} {
                              set linedos $line
                              $log configure -state normal
                              $log delete 2.0
                              $log insert 2.0 "$line\n"
                              $log configure -state disabled
                         }

                    }
                    if {[string first "Percentage" $line] ne -1}  {
                         if {$linetres ne $line} {
                              set linetres $line
                              $log configure -state normal
                              $log delete 3.0
                              $log insert 3.0 "$line\n"
                              $log configure -state disabled
                         }

                    }
                    if {[string first "Time Step" $line] ne -1}  {
                         if {$linecuat ne $line} {
                              set linecuat $line
                              $log configure -state normal
                              $log delete 4.0
                              $log insert 4.0 "$line\n"
                              $log configure -state disabled
                         }

                    }
                    if {$line eq "ANALYSIS COMPLETED"} {
                         set cont 0
                    }
                    # condicion de salida
               }
               close $f
          }
          if {$cont} {after [lindex [split [expr {$screenOutTime *1000}] "."] 0] ::kwiz::WhileCalculation $log}

     }
     update idletasks

     update

}

proc ::kwiz::Terminate { } {
     GiD_Process Utilities CancelProcess
}

#-np- ::kwiz::CreateDEMWizard
 proc ::kwiz::CreateDEMWizard { } {
     variable bwinpath
     variable wprops
     variable stepidlist

     ::kwiz::blocking "unlock"

     # Destroy the window
     if {[winfo exists $bwinpath]} {
          destroy $bwinpath
     }
     # wa "bwinpath:$bwinpath"
     # Create the window
     InitWindow $bwinpath [= "DEM Wizard"] Pre::kwiz::CreateDEMWizardWindowGeom

     # Centrar window
     wm withdraw $bwinpath

     set x [expr {([winfo screenwidth .gid.central.s]-[winfo width .gid.central.s])/2}]
     set y [expr {([winfo screenheight .gid.central.s]-[winfo height .gid.central.s])/2}]
     if { $x < 0 } { set x 0 }
     if { $y < 0 } { set y 0 }
     WmGidGeom $bwinpath +$x+$y
     update
     wm deiconify $bwinpath

     # Tamaños !
     wm minsize $bwinpath 855 715
     #wm maxsize $bwinpath 855 715

     # First destroy all defined command (snit step data type)
     foreach cmdid [info commands ::kwiz::step*] {
          $cmdid destroy
     }

     # Create all the steps

     # New Load
     ::snitwiz::wizardstep ::kwiz::step1 -title [= "Welcome to DEM"] -body {::kwiz::bodystep1 $win}

     # Wizard type
     ::snitwiz::wizardstep ::kwiz::step1bis -title [= "Step 1: Wizard type"] -body {::kwiz::bodystep1bis $win}

     # Experiment type
     ::snitwiz::wizardstep ::kwiz::step2  -title [= "Step 2: Experiment type"] -layout basic -body {::kwiz::bodystep2 $win}

     # Geometry & Mesh
     ::snitwiz::wizardstep ::kwiz::step3  -title [= "Step 3: Geometry & Mesh"] -layout basic -body {::kwiz::bodystep3 $win}

     # Material parameters
     ::snitwiz::wizardstep ::kwiz::step4  -title [= "Step 4: Material parameters"] -layout basic -body {::kwiz::bodystep4 $win}

     # General setting
     ::snitwiz::wizardstep ::kwiz::step5  -title [= "Step 5: General setting"] -layout basic -body {::kwiz::bodystep5 $win}

     # Output settings"
     ::snitwiz::wizardstep ::kwiz::step6  -title [= "Step 6: Output settings"] -layout basic -body {::kwiz::bodystep6 $win}

     # Run simulation
     ::snitwiz::wizardstep ::kwiz::step7  -title [= "Step 7: Run simulation"] -layout basic -body {::kwiz::bodystep7 $win}

     # Result visualization
     ::snitwiz::wizardstep ::kwiz::step8  -title [= "Step 8: Result visualization"] -layout basic -body {::kwiz::bodystep8 $win}

     # Render the wizard
     ::snitwiz::wizard $bwinpath.w -steps $stepidlist

     #::snitwiz::disablenext
     $bwinpath.w buttonstate next disabled

     # Start the wizard
     $bwinpath.w start


     #::snitwiz::disablenext
      $bwinpath.w buttonstate next disabled
}

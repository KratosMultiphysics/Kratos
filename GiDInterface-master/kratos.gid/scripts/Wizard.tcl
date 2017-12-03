##################################################################################
#   This file is common for all Kratos Applications.
#   Do not change anything here unless it's strictly necessary.
##################################################################################

namespace eval Wizard {
    # Namespace variables declaration
    variable wizDoc
    variable wizardid
    variable wizwindow
    variable wprops
    variable stepidlist
}

proc Wizard::Init { } {
    package require gid_wizard
    package require wcb

    variable wizDoc
    set wizDoc ""
    
    variable wizardid
    set wizardid ""
    
    variable wizwindow
    set wizwindow ".gid.activewizard"
    
    variable wprops
    catch {unset wprops}
    set wprops(dummy) 0
    
    variable stepidlist
    set stepidlist [list ]
}

proc Wizard::LoadWizardDoc {doc_path} {
    variable wizDoc
    set wizDoc [dom parse [tDOM::xmlReadFile $doc_path] ]
    #W "Wizdoc"
    #W [$wizDoc asXML]
}

proc Wizard::CreateWindow {} {
    variable wizwindow
    variable wprops
    variable stepidlist
    variable wizardid
    # W "lista de pasos $stepidlist"
    # Destroy the window 
    if {[winfo exists $wizwindow]} {
         destroy $wizwindow
    }
    
    # Create the window
    InitWindow $wizwindow [= "Wizard"] Pre::kwiz::CreateWizardWindowGeom
   
    #set ::kwiz::wizardid "cutting"   
    
    # Centrar window
    wm withdraw $wizwindow
    
    set x [expr {([winfo screenwidth .gid.central.s]-[winfo width .gid.central.s])/2}]
    set y [expr {([winfo screenheight .gid.central.s]-[winfo height .gid.central.s])/2}]
    if { $x < 0 } { set x 0 }
    if { $y < 0 } { set y 0 }
    WmGidGeom $wizwindow +$x+$y
    update
    wm deiconify $wizwindow

    # Tamaños !
    SetWindowSize 700 500
    
    # First destroy all defined command (snit step data type)
    #foreach cmdid [info commands ::Wizard::*] {
    #     $cmdid destroy
    #}
    
    # Create all the steps
    set i 0
    set nssteplist [list ]
    foreach stepId $stepidlist {
        incr i
        catch {$stepId destroy}
        lappend nssteplist ::Wizard::${stepId}
        
        ::gid_wizard::wizardstep $stepId -title [= "Step $i: $stepId"] -layout basic -body "::${wizardid}::Wizard::$stepId \$win"
    }
    
    # Render the wizard
    # W "lista a enviar $nssteplist"
    ::gid_wizard::wizard $wizwindow.w -steps $nssteplist
         
    # Start the wizard
    $wizwindow.w start
    
    bind $wizwindow <<WizNext>> [list Wizard::NextEvent %d]
}

proc Wizard::DestroyWindow {} {
    variable wizwindow
    
    if {[winfo exists $wizwindow]} {destroy $wizwindow}
    return ""
}

proc Wizard::NextEvent {step} {
    variable wizardid
    variable stepidlist
    set stepId [lindex $stepidlist $step]
    ::${wizardid}::Wizard::Next$stepId
}


proc Wizard::ImportWizardData {} {
    variable wizDoc
    variable wprops
    variable stepidlist
    variable wizardid
    
    # ABSTRACT: Import all wizard data variables from XML 
    set xmlData $wizDoc
    #wa [$KPriv(xmlWiz) asXML]
    set wizardid [[$xmlData selectNodes "/Wizard"] getAttribute "wizardid"]
    set path "/Wizard/Steps"
    set stepNodes [$xmlData selectNodes "$path/Step"]
    set dataNodes [$xmlData selectNodes "$path/Step/Data"]
    #wa "Sn $stepNodes"
    foreach stepNode $stepNodes dataNode $dataNodes {
        set i 0
        set stepId [$stepNode getAttribute id ""]
        lappend stepidlist $stepId
        foreach node [$dataNode childNodes] {
            # For nodes with no children
            if {([$node nodeName] eq "Item") && (![$node hasChildNodes])} {
                set n [$node getAttribute n ""]
                set pn [$node getAttribute pn $n]
                set v [$node getAttribute v ""]
                set t [$node getAttribute type ""]
                #W "::kwiz::wprops($stepId,$n,value)= $v -> $i"
                set wprops($stepId,$n,value) $v
                set wprops($stepId,$n,order) $i
                set wprops($stepId,$n,type) $t
                set wprops($stepId,$n,name) $pn
                
                if {$t eq "combo"} {
                    set vs [split [$node getAttribute values ""] ","]
                    set wprops($stepId,$n,values) $vs
                }
                
                incr i                
            }
        }
    }
    #W [array names wprops]
}

proc Wizard::SetProperty { stepid propid value } {
    variable wprops
    set wprops($stepid,$propid,value) $value
}
proc Wizard::GetProperty { stepid propid } {
    variable wprops
    return $wprops($stepid,$propid)
}
proc Wizard::GetStepProperties { stepid } {
    variable wprops
    set lista [list ]
    foreach key [array names wprops] {
        if {[lindex [split $key ","] 0] eq $stepid} {if {[lindex [split $key ","] 1] ni $lista} {lappend lista [lindex [split $key ","] 1]}}
    }
    return $lista
}

proc Wizard::SetWindowSize {x y} {
    variable wizwindow
    wm minsize $wizwindow $x $y
    wm maxsize $wizwindow $x $y
    
}

Wizard::Init

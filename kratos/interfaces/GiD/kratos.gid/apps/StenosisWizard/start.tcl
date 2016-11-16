namespace eval ::StenosisWizard {
    # Variable declaration
    variable dir
}

proc ::StenosisWizard::Init { } {
    # Variable initialization
    variable dir
    
    # Init Working directory
    set dir [apps::getMyDir "StenosisWizard"]
    # We'll work on 3D space
    spdAux::SetSpatialDimmension "3D"
    # Load Fluid App
    apps::LoadAppById "Fluid"
    # Don't open the tree
    set ::spdAux::TreeVisibility 0
    
    # Enable the Wizard Module
    Kratos::LoadWizardFiles
    LoadMyFiles
}

proc ::StenosisWizard::LoadMyFiles { } {
    variable dir
    uplevel #0 [list source [file join $dir xml GetFromXML.tcl]]
    uplevel #0 [list source [file join $dir write write.tcl]]
    ::Wizard::LoadWizardDoc [file join $dir wizard Wizard_default.wiz]
    uplevel #0 [list source [file join $dir wizard Wizard_Steps.tcl]]
    Wizard::ImportWizardData
    
    
    # Init the Wizard Window
    after 600 [::StenosisWizard::StartWizardWindow]
}


proc ::StenosisWizard::StartWizardWindow { } {
    gid_groups_conds::close_all_windows
    Wizard::CreateWindow
    
}

::StenosisWizard::Init

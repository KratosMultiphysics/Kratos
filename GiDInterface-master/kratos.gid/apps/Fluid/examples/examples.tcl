namespace eval Fluid::examples {

}

proc Fluid::examples::Init { } {
    uplevel #0 [list source [file join $::Fluid::dir examples CylinderInFlow.tcl]]
    #uplevel #0 [list source [file join $::FSI::dir examples HorizontalFlexibleBar.tcl]]
    GiDMenu::InsertOption "Kratos" [list "---"] 6 PRE "" "" "" replace =
    GiDMenu::InsertOption "Kratos" [list "Cylinder in air flow" ] 7 PRE [list ::Fluid::examples::CylinderInFlow] "" "" replace =
    #GiDMenu::InsertOption "Kratos" [list "Horizontal flexible bar" ] 8 PRE [list ::FSI::examples::HorizontalFlexibleBar] "" "" replace =
    GiDMenu::UpdateMenus
}

Fluid::examples::Init
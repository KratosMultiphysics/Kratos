namespace eval Dam::examples {

}

proc Dam::examples::Init { } {
    uplevel #0 [list source [file join $::Dam::dir examples ThermoMechaDam.tcl]]
    #uplevel #0 [list source [file join $::FSI::dir examples HorizontalFlexibleBar.tcl]]
    GiDMenu::InsertOption "Kratos" [list "---"] 6 PRE "" "" "" replace =
    GiDMenu::InsertOption "Kratos" [list "Thermo-Mechanical Dam" ] 7 PRE [list ::Dam::examples::ThermoMechaDam] "" "" replace =
    #GiDMenu::InsertOption "Kratos" [list "Horizontal flexible bar" ] 8 PRE [list ::FSI::examples::HorizontalFlexibleBar] "" "" replace =
    GiDMenu::UpdateMenus
}

Dam::examples::Init

namespace eval EmbeddedFluid::examples {
    variable CylinderInFlow_Data
}

proc EmbeddedFluid::examples::Init { } {
    uplevel #0 [list source [file join $::EmbeddedFluid::dir examples CylinderInFlow.tcl]]
    GiDMenu::InsertOption "Kratos" [list "---"] 6 PRE "" "" "" replace =
    GiDMenu::InsertOption "Kratos" [list "Embedded cylinder test" ] 7 PRE [list ::EmbeddedFluid::examples::CylinderInFlow] "" "" replace =
    GiDMenu::UpdateMenus
}

EmbeddedFluid::examples::Init
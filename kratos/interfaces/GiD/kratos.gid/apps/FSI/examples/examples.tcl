namespace eval FSI::examples {

}

proc FSI::examples::Init { } {
    uplevel #0 [list source [file join $::FSI::dir examples MokChannelWithFlexibleWall.tcl]]
    uplevel #0 [list source [file join $::FSI::dir examples HorizontalFlexibleBar.tcl]]
    GiDMenu::InsertOption "Kratos" [list "---"] 6 PRE "" "" "" replace =
    GiDMenu::InsertOption "Kratos" [list "Mok - Channel with flexible wall" ] 7 PRE [list ::FSI::examples::MokChannelFlexibleWall] "" "" replace =
    GiDMenu::InsertOption "Kratos" [list "Horizontal flexible bar" ] 8 PRE [list ::FSI::examples::HorizontalFlexibleBar] "" "" replace =
    GiDMenu::UpdateMenus
}

FSI::examples::Init
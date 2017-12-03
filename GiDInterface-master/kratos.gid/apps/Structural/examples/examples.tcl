namespace eval Structural::examples {

}

proc Structural::examples::Init { } {
    uplevel #0 [list source [file join $::Structural::dir examples TrussCantilever.tcl]]
    GiDMenu::InsertOption "Kratos" [list "---"] 6 PRE "" "" "" replace =
    GiDMenu::InsertOption "Kratos" [list "Truss cantilever" ] 7 PRE [list ::Structural::examples::TrussCantilever] "" "" replace =
    GiDMenu::UpdateMenus
}

Structural::examples::Init
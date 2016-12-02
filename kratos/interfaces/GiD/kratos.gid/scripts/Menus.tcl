
proc Kratos::ToolbarAddItem {id icon code tex} {
    variable kratos_private
    if {![info exists kratos_private(MenuItems)]} {
        set kratos_private(MenuItems) [dict create]
    }
    set num [llength [dict keys $kratos_private(MenuItems)]]
    incr num
    dict set kratos_private(MenuItems) $num id $id
    dict set kratos_private(MenuItems) $num icon $icon
    dict set kratos_private(MenuItems) $num code $code
    dict set kratos_private(MenuItems) $num tex $tex
    return $num
}
proc Kratos::ToolbarDeleteItem {id} {
    variable kratos_private
    foreach num [dict keys $kratos_private(MenuItems)] {
        if {[dict get $kratos_private(MenuItems) $num id] eq $id } {
            set kratos_private(MenuItems) [dict remove $kratos_private(MenuItems) $num]
            break
        }
    }
    return $num
}

proc Kratos::ToolbarRefresh {} {
    Kratos::EndCreatePreprocessTBar
    Kratos::CreatePreprocessModelTBar
}


proc Kratos::CreatePreprocessModelTBar { {type "DEFAULT INSIDELEFT"} } {
    global KBitmapsNames KBitmapsCommands KBitmapsHelp
    variable kratos_private
    
    Kratos::ToolbarAddItem "Model" "propstree.png" [list -np- gid_groups_conds::open_conditions menu] [= "Define the model properties"]
    Kratos::ToolbarAddItem "Spacer" "" "" ""
    Kratos::ToolbarAddItem "Run" "run.png" {Utilities Calculate} [= "Run the simulation"]
    Kratos::ToolbarAddItem "Output" "output.png" [list -np- PWViewOutput] [= "View process info"]
    Kratos::ToolbarAddItem "Stop" "stop.png" {Utilities CancelProcess} [= "Cancel process"]
    Kratos::ToolbarAddItem "SpacerApp" "" "" ""
    
    set app_items_toolbar [apps::ExecuteOnCurrentApp CustomToolbarItems]
    if {$app_items_toolbar < 1} {
        Kratos::ToolbarDeleteItem "SpacerApp"
    }
    if {$app_items_toolbar ne "-1"} {
        set dir [file join $::Kratos::kratos_private(Path) images ]
        set iconslist [list ]
        set commslist [list ]
        set helpslist [list ]
        foreach item [dict keys $kratos_private(MenuItems)] {
            set icon [dict get $kratos_private(MenuItems) $item icon]
            lappend iconslist [expr {$icon ne "" ? [file join $dir $icon] : "---"}]
            lappend commslist  [dict get $kratos_private(MenuItems) $item code]
            lappend helpslist [dict get $kratos_private(MenuItems) $item tex]
        }
        
        set KBitmapsNames(0) $iconslist
        set KBitmapsCommands(0) $commslist
        set KBitmapsHelp(0) $helpslist
        
        set prefix Pre
        set name KPreprocessModelbar
        set procname ::Kratos::CreatePreprocessModelTBar
        set kratos_private(ToolBars,PreprocessModelTBar) [CreateOtherBitmaps ${name} [= "Kratos toolbar"] KBitmapsNames KBitmapsCommands KBitmapsHelp $dir $procname $type $prefix]
        
        AddNewToolbar [= "Kratos toolbar"] ${prefix}${name}WindowGeom $procname
    }
}

proc Kratos::EndCreatePreprocessTBar {} {
    variable kratos_private
   
    set name KPreprocessModelbar
    
    ReleaseToolbar ${name}
    if {[info exists kratos_private(ToolBars,PreprocessModelTBar)]} {
        destroy $kratos_private(ToolBars,PreprocessModelTBar)
    }
    if {[info exists kratos_private(MenuItems)]} {
        unset kratos_private(MenuItems)
    }
    update
}




proc Kratos::ChangeMenus { } {
    set found [GiDMenu::_FindIndex "Kratos" PRE]
    if {$found > 0} {GiDMenu::Delete "Kratos" PRE}
    GiDMenu::Create "Kratos" PRE
    variable kratos_private
    set tomode "developer mode"
    set fromode "release mode"
    if {$kratos_private(DevMode) eq "dev"} {set tomode "release mode";set fromode "developer mode"}
    GiDMenu::InsertOption "Kratos" [list "Kratos data" ] 0 PRE [list gid_groups_conds::open_conditions menu] "" "" replace =
    GiDMenu::InsertOption "Kratos" [list "---"] 1 PRE "" "" "" replace =
    GiDMenu::InsertOption "Kratos" [list "You are in $fromode" ] 2 PRE [list ] "" "" replace =
    GiDMenu::InsertOption "Kratos" [list "Switch to $tomode" ] 3 PRE [list Kratos::SwitchMode] "" "" replace =

    if {$::Kratos::kratos_private(UseWizard)} {
        GiDMenu::InsertOption "Kratos" [list "---"] 4 PRE "" "" "" replace =
        GiDMenu::InsertOption "Kratos" [list "Wizard window" ] 5 PRE [list Wizard::CreateWindow] "" "" replace =
    }
    GidChangeDataLabel "Data units" ""
    GidChangeDataLabel "Interval" ""
    GidChangeDataLabel "Conditions" ""
    GidChangeDataLabel "Materials" ""
    GidChangeDataLabel "Interval Data" ""
    GidChangeDataLabel "Problem Data" ""
    GidChangeDataLabel "Local axes" "gid_groups_conds::local_axes_menu %W"
    
    GiDMenu::InsertOption "Help" [list ---] end PREPOST {} "" "" insertafter
    GiDMenu::InsertOption "Help" [list [_ "Visit %s web" Kratos]...] end PREPOST [list VisitWeb "http://www.cimne.com/kratos"] "" "" insertafter
    
    GiDMenu::UpdateMenus
}


 

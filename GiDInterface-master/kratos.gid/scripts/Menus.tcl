
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
    Kratos::EndCreatePreprocessTBar
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
        set theme [gid_themes::GetCurrentTheme]
        set dir [file join $::Kratos::kratos_private(Path) images ]
        set iconslist [list ]
        set commslist [list ]
        set helpslist [list ]
        foreach item [dict keys $kratos_private(MenuItems)] {
            set icon [dict get $kratos_private(MenuItems) $item icon]
            if {![file exists $icon] && $icon ne ""} {
                set good_dir $dir
                if {$theme eq "GiD_black"} {
                    set good_dir [file join $dir Black]
                    if {![file exists [file join $good_dir $icon]]} {set good_dir $dir}
                }
                set icon [file join $good_dir $icon]
            }
            
            lappend iconslist [expr {$icon ne "" ? $icon : "---"}]
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
    set pos -1
    if {$kratos_private(DevMode) eq "dev"} {set tomode "release mode";set fromode "developer mode"}
    GiDMenu::InsertOption "Kratos" [list "Kratos data" ] [incr pos] PRE [list gid_groups_conds::open_conditions menu] "" "" replace =
    GiDMenu::InsertOption "Kratos" [list "Local axes window" ] [incr pos] PRE [list gid_groups_conds::local_axes_window] "" "" replace =
    GiDMenu::InsertOption "Kratos" [list "---"] [incr pos] PRE "" "" "" replace =
    GiDMenu::InsertOption "Kratos" [list "You are in $fromode" ] [incr pos] PRE [list ] "" "" replace =
    GiDMenu::InsertOption "Kratos" [list "Switch to $tomode" ] [incr pos] PRE [list Kratos::SwitchMode] "" "" replace =

    if {$::Kratos::kratos_private(UseWizard)} {
        GiDMenu::InsertOption "Kratos" [list "---"] [incr pos] PRE "" "" "" replace =
        GiDMenu::InsertOption "Kratos" [list "Wizard window" ] [incr pos] PRE [list Wizard::CreateWindow] "" "" replace =
    }
    GiDMenu::InsertOption "Kratos" [list "---"] [incr pos] PRE "" "" "" replace =
    GiDMenu::InsertOption "Kratos" [list "About Kratos" ] [incr pos] PRE [list Kratos::About] "" "" replace =
    GidChangeDataLabel "Data units" ""
    GidChangeDataLabel "Interval" ""
    GidChangeDataLabel "Conditions" ""
    GidChangeDataLabel "Materials" ""
    GidChangeDataLabel "Interval Data" ""
    GidChangeDataLabel "Problem Data" ""
    GidChangeDataLabel "Local axes" ""
    # GidChangeDataLabel "Local axes" "gid_groups_conds::local_axes_window"
    
    GiDMenu::InsertOption "Help" [list ---] end PREPOST {} "" "" insertafter
    GiDMenu::InsertOption "Help" [list [_ "Visit %s web" Kratos]] end PREPOST [list VisitWeb "http://www.cimne.com/kratos"] "" "" insertafter
    GiDMenu::InsertOption "Help" [list [_ "Do you want to develop Kratos?"]] end PREPOST [list VisitWeb "https://github.com/KratosMultiphysics"] "" "" insertafter
    
    GiDMenu::UpdateMenus
}

proc Kratos::About {} {
    Splash
}
 

proc Kratos::Splash { } {
    variable kratos_private
    set prev_splash_state [GiD_Set SplashWindow]
    GiD_Set SplashWindow 1 ;#set temporary to 1 to force show splash without take care of the GiD splash preference
    set off_x 150
    set fnt "Sans 10"
    if { $::tcl_platform(platform) == "windows" } {
        set fnt "verdana 10"
        set off_x 130
    }
    set line1 "Kratos Multiphysics Version $kratos_private(Version)"
    ::GidUtils::Splash [file join $kratos_private(Path) images splash.png] .splash 0 \
        [list $line1 $off_x 230]


    .splash.lv configure -font $fnt -background white -foreground black \
        -relief solid -borderwidth 1 -padx 12 -pady 3
    update
    
    GiD_Set SplashWindow $prev_splash_state
}

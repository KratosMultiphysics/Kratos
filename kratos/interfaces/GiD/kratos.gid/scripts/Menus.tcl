
proc Kratos::CreatePreprocessModelTBar { {type "DEFAULT INSIDELEFT"} } {
    global KBitmapsNames KBitmapsCommands KBitmapsHelp
    variable kratos_private
    set dir [file join $::Kratos::kratos_private(Path) images ]
   
    set KBitmapsNames(0) "$dir/new_props.gif --- $dir/openrunsim.gif $dir/runsimulation.png $dir/runsiminfo.gif $dir/stop.png"
            
    set KBitmapsCommands(0) [list [list -np- gid_groups_conds::open_conditions menu] "" [list -np- RunWin] {Utilities Calculate} [list -np- PWViewOutput] {Utilities CancelProcess}]

    set KBitmapsHelp(0) [list [= "Define the model properties"]  "" [= "Open the process control window"] [= "Run the simulation"] [= "View process info"] [= "Cancel process"]]         
    
    set prefix Pre
    set name KPreprocessModelbar
    set procname ::Kratos::CreatePreprocessModelTBar
    set kratos_private(ToolBars,PreprocessModelTBar) [CreateOtherBitmaps ${name} [= "Kratos toolbar"] KBitmapsNames KBitmapsCommands KBitmapsHelp $dir $procname $type $prefix]
    
    AddNewToolbar [= "Kratos toolbar"] ${prefix}${name}WindowGeom $procname
}

proc Kratos::EndCreatePreprocessTBar {} {
    variable kratos_private
   
    set name KPreprocessModelbar
    
    ReleaseToolbar ${name}
    #if {[info exists kratos_private(ToolBars,PreprocessModelTBar)]} {
        destroy $kratos_private(ToolBars,PreprocessModelTBar)
    #}
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


 

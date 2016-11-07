##################################################################################
#   This file is common for all Kratos Applications.
#   Do not change anything here unless it's strictly necessary.
##################################################################################

set ::kratos_debug 1 ;# could be 0,1,-1

##########################################################
#################### GiD Tcl events ######################
##########################################################
proc InitGIDProject { dir } {
    font configure SmallFont -size 16
    #uplevel #0 [list Kratos::InitGIDProject $dir]
    Kratos::InitGIDProject $dir
}

proc EndGIDProject {} {
    Kratos::DestroyWindows
    spdAux::EndRefreshTree
    Kratos::RegisterEnvironment
    Model::DestroyEverything
    gid_groups_conds::end_problemtype [Kratos::GiveKratosDefaultsFile]
    unset -nocomplain ::Kratos::kratos_private
}

proc ChangedLanguage { newlan } {
    Kratos::ChangeMenus
}

proc InitGIDPostProcess {} {
    gid_groups_conds::close_all_windows
    gid_groups_conds::open_post check_default
}

proc EndGIDPostProcess {} {
    gid_groups_conds::close_all_windows
    if {$::spdAux::TreeVisibility} {
        gid_groups_conds::open_conditions check_default
        gid_groups_conds::open_conditions menu
    }
}
 
# Load GiD project files (initialise XML Tdom structure)
proc LoadGIDProject { filespd } {
    if { ![file exists $filespd] } { return }
    set versionPT [gid_groups_conds::give_data_version]
    gid_groups_conds::open_spd_file $filespd
    set versionData [gid_groups_conds::give_data_version]
    if { [package vcompare $versionPT $versionData] == 1 } {
        after idle Kratos::upgrade_problemtype
    }
    #spdAux::reactiveApp
    update
    spdAux::LoadModelFiles
}

# Save GiD project files (save XML Tdom structure to spd file)
proc SaveGIDProject { filespd } {
    gid_groups_conds::save_spd_file $filespd
    Kratos::RegisterEnvironment
    FileSelector::CopyFilesIntoModel [file dirname $filespd]
}

proc BeforeTransformProblemType { file oldproblemtype newproblemtype } {

return "-cancel-"
}

proc AfterTransformProblemType { filename oldproblemtype newproblemtype } {
    set spd_file [file join $filename.gid [file tail $filename].spd]
    return [gid_groups_conds::transform_problemtype $spd_file]
}

proc AfterWriteCalcFileGIDProject { filename errorflag } {
    FileSelector::CopyFilesIntoModel [file dirname $filename]
    catch {write::Init}
    set errcode [::write::writeEvent $filename]
    if {$errcode} {return "-cancel-"}
}
proc BeforeMeshGeneration { elementsize } {
    Kratos::BeforeMeshGeneration $elementsize
}

##########################################################
#################### Kratos namespace ####################
##########################################################
namespace eval Kratos {
  variable kratos_private
}

proc Kratos::InitGIDProject { dir } {
    variable kratos_private
    unset -nocomplain kratos_private
    set kratos_private(Path) $dir ;#to know where to find the files
    set kratos_private(DevMode) "release" ; #can be dev or release
   
    array set kratos_private [ReadProblemtypeXml [file join $dir kratos.xml] Infoproblemtype {Name Version MinimumGiDVersion}]
    if { [GidUtils::VersionCmp $kratos_private(MinimumGiDVersion)] < 0 } {
        WarnWin [_ "Error: %s Interface requires GiD %s or later." $kratos_private(Name) $kratos_private(MinimumGiDVersion)]
    }
    if {[GiD_Info GiDVersion] eq "13.0-rc1"} { WarnWin "The minimum GiD version is 13.0-rc2.\n Ask the GiD Team for it."; return ""}
        
    #append to auto_path only folders that must include tcl packages (loaded on demand with package require mechanism)
    if { [lsearch -exact $::auto_path [file join $dir scripts]] == -1 } {
        lappend ::auto_path [file join $dir scripts]
    }
    # JG Sources will be in a different proc
    foreach filename {Applications.tcl Writing.tcl spdAuxiliar.tcl} {
        uplevel 1 [list source [file join $dir scripts $filename]]
    }

    # JG Sources will be in a different proc
    foreach filename {Model.tcl Entity.tcl Parameter.tcl Topology.tcl Solver.tcl ConstitutiveLaw.tcl Condition.tcl Element.tcl SolutionStrategy.tcl Process.tcl} {
        uplevel 1 [list source [file join $dir scripts Model $filename]]
    }
    # JG Sources will be in a different proc
    foreach filename {SimpleXMLViewer.tcl FileManager.tcl} {
        uplevel 1 [list source [file join $dir libs $filename]]
    }
    
    set kratos_private(UseWizard) 0
     
    Kratos::load_gid_groups_conds
    customlib::UpdateDocument
    Kratos::LoadEnvironment
    Kratos::ChangeMenus
    #set HeaderBackground [$doc selectNodes string(Infoproblemtype/Program/HeaderBackground)]
    #gid_groups_conds::SetHeaderBackground $HeaderBackground
    gid_groups_conds::SetLibDir [file join $dir exec] 
    gid_groups_conds::begin_problemtype [file join $dir kratos_default.spd] [Kratos::GiveKratosDefaultsFile]

    spdAux::processIncludes
    spdAux::parseRoutes
    update
    spdAux::LoadModelFiles
    
    after 100 [list gid_groups_conds::close_all_windows]
    after 500 [list spdAux::CreateWindow]
}

proc Kratos::DestroyWindows {} {
    spdAux::DestroyWindow
    if {$::Kratos::kratos_private(UseWizard)} {
        Wizard::DestroyWindow
    }
}

proc Kratos::LoadWizardFiles { } {
    set ::Kratos::kratos_private(UseWizard) 1
    set dir $::Kratos::kratos_private(Path)
    uplevel #0 [list source [file join $dir scripts Wizard.tcl]]
    Kratos::ChangeMenus
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
    GiDMenu::UpdateMenus
}

proc Kratos::SwitchMode {} {
    variable kratos_private
    if {$kratos_private(DevMode) eq "dev"} {
        set kratos_private(DevMode) "release"
    }  {
        set kratos_private(DevMode) "dev"
    }
    Kratos::RegisterEnvironment
    #W "Registrado $kratos_private(DevMode)"
    Kratos::ChangeMenus
    spdAux::RequestRefresh
}

proc Kratos::GetPreferencesFilePath { } {
    variable kratos_private
    set dir_name [file dirname [GiveGidDefaultsFile]]
    set file_name $kratos_private(Name)Vars.txt
    if { $::tcl_platform(platform) == "windows" } {
        return [file join $dir_name $file_name]
    } else {
        return [file join $dir_name .$file_name]
    }
}

proc Kratos::RegisterEnvironment { } {
    variable kratos_private
    set varsToSave [list DevMode]
    set preferences [dict create]
    dict set preferences DevMode $kratos_private(DevMode)
    #gid_groups_conds::set_preference DevMode $kratos_private(DevMode)
    set fp [open [Kratos::GetPreferencesFilePath] w]
    catch {set data [puts $fp [write::tcl2json $preferences]]}
    close $fp
}
proc Kratos::LoadEnvironment { } {
    variable kratos_private
    #set kratos_private(DevMode) [gid_groups_conds::get_preference DevMode releasedefault]
    set data ""
    set syspath HOME
    if {$::tcl_platform(platform) eq "windows"} {set syspath APPDATA}
    catch {
        set fp [open [Kratos::GetPreferencesFilePath] r]
        set data [read $fp]
        close $fp
    }
    foreach {k v} [write::json2dict $data] {
        set kratos_private($k) $v
    }
}

proc Kratos::load_gid_groups_conds {} {  
    package require customlib_extras ;#this require also customLib
    package require customlib_native_groups
    package require json::write
}

proc Kratos::GiveKratosDefaultsFile {} {
    variable kratos_private
    set dir_name [file dirname [GiveGidDefaultsFile]]
    set file_name $kratos_private(Name)$kratos_private(Version).ini
    if { $::tcl_platform(platform) == "windows" } {
        return [file join $dir_name $file_name]
    } else {
        return [file join $dir_name .$file_name]
    }
}

proc Kratos::upgrade_problemtype {} {
    set w [dialogwin_snit .gid._ask -title [_ "Action"] -entrytext \
            [_ "The model needs to be upgraded. Do you want to upgrade to new version?"]]
    set action [$w createwindow]
    destroy $w
    if { $action < 1 } { return }
    set project [lindex [GiD_Info Project] 0]
    GiD_Process escape escape escape escape Data Defaults TransfProblem $project
}

proc Kratos::LocalAxesMenu { menu } {    
    #if { [$menu index end] ne "none" } { return }
    $menu delete 0 end
    set local_axes [GiD_Info localaxes]
    foreach i [list Point Line Surface] name [list [_ Points] [_ Lines] [_ Surfaces]] {
        $menu add cascade -label $name -menu $menu.m$i
        destroy $menu.m$i
        set m [menu $menu.m$i -tearoff 0]
        if { [lsearch "Line Surface" $i] != -1 } {
            $m add command -label [_ "Assign Automatic"] -command [list GiD_Process \
                    escape escape escape escape Data Conditions AssignCond ${i}_Local_axes \
                    change -Automatic-]
            $m add command -label [_ "Assign Automatic alt"] -command [list GiD_Process \
                    escape escape escape escape Data Conditions AssignCond ${i}_Local_axes \
                    change -Automatic_alt-]
            $m add separator
        }
        set idx 0
        foreach j $local_axes {
            $m add command -label [_ "Assign '%s'" $j] -command [list GiD_Process \
                    escape escape escape escape Data Conditions AssignCond ${i}_Local_axes \
                    change $j]
            incr idx
        }
        if { $idx } { $m add separator }
        $m add command -label [_ "Unassign"] -command [list GiD_Process \
                escape escape escape escape Data Conditions AssignCond ${i}_Local_axes \
                Unassign]
        
        $m add separator
        $m add command -label [_ Draw] -command [list GiD_Process \
                escape escape escape escape Data Conditions DrawCond -LocalAxes- \
                ${i}_Local_axes -draw-]
    }
    set ns [list [_ "Define#C#menu"] --- [_ "Draw#C#menu"] [_ "Draw all#C#menu"] \
            --- [_ "Delete#C#menu"] [_ "Delete all#C#menu"]]
    set cs {
        "Data LocalAxes DefineLocAxes"
        {} "Data LocalAxes DrawLocAxes"
        "Data LocalAxes DrawLocAxes -All-"
        {} "Data LocalAxes DeleteLA"
        "Data LocalAxes DeleteAllLA"
    }
    $menu add separator
    foreach n $ns c $cs {
        if { $n eq "---" } {
            $menu add separator
        } else {
            $menu add command -label $n -command [concat "GiD_Process escape escape \
                    escape escape" $c]
        }
    }
}


#######################################################
#    Add items to control-1 contextual menu 
#######################################################

proc Kratos::AddDataItemsToMenu { menu whatuse type entity } {

    set DisableGraphics [.central.s disable graphics]
    set DisableWarnLine [.central.s disable warnline]
    
    .central.s disable graphics 1
    .central.s disable warnline 1

    switch $type {
        Points { set cnds point_groups }
        Lines { set cnds line_groups }
        Surfaces { set cnds surface_groups }
        Volumes { set cnds volume_groups }
        Nodes { set cnds point_groups }
        Elements { set cnds "line_groups surface_groups volume_groups" }
    }
    switch $type Points - Lines - Surfaces - Volumes { set where geometry } Nodes - \
        Elements { set where mesh }

    set doc $gid_groups_conds::doc

    set submenu ""
    set needsseparator 0
    foreach cnd $cnds {
        set ret [GiD_Info Conditions $cnd $where $entity]
        foreach i $ret {
            foreach "num face - group" $i break
            if { $num == 0 } { continue }
            if { $num eq "E" } { set num $face }
            if { $submenu == "" } {
                set j 1
                while { [winfo exists $menu.m$j] } { incr j }
                set submenu [menu $menu.m$j -tearoff 0]
                set text [_ Groups]
                $menu add cascade -label $text -menu $submenu
                set needsseparator 1
            }
            set xp [format_xpath {//group[@n=%s and ancestor::condition]} \
                    $group]
            set domNodes [$doc selectNodes $xp]

            set txt [_ "Group: %s" $group]
            if { ![llength $domNodes] } {
                ContextualEntity::AddToMenu $submenu $txt \
                    [list Kratos::_AddDataItemsToMenuItem $group]
            } else {
                set j 1
                while { [winfo exists $submenu.m$j] } { incr j }
                set subsubmenu [menu $submenu.m$j -tearoff 0]
                $submenu add cascade -label $txt -menu $subsubmenu
                ContextualEntity::AddToMenu $subsubmenu $txt \
                    [list Kratos::_AddDataItemsToMenuItem $group]
                foreach domNode $domNodes {
                    set cndNode [$domNode selectNodes ancestor::condition]
                    set blockNode [$domNode selectNodes {ancestor::blockdata[@sequence='1']}]
                    if { $blockNode ne "" } {
                        set txt [_ "Condition: %s - %s" [$blockNode @name] [$cndNode @pn]]
                    } else {
                        set txt [_ "Condition: %s" [$cndNode @pn]]
                    }
                    ContextualEntity::AddToMenu $subsubmenu $txt \
                        [list Kratos::_AddDataItemsToMenuItemD $group \
                            $domNode]
                }
            }
        }
    }
    if { $needsseparator } { $menu add separator }
    
    if { !$DisableGraphics } { .central.s disable graphics 0 }
    if { !$DisableWarnLine } { .central.s disable warnline 0 }

}

proc Kratos::_AddDataItemsToMenuItem { group } {
    gid_groups_conds::open_groups .gid window_force window
    set wg $gid_groups_conds::gid_group_var
    set tree [$wg givetreectrl]   
    $tree selection clear
    foreach item [$tree item children root] {
        if { [$tree item text $item 0] eq $group } {
            $tree selection add $item
            break
        }
    }
}

proc Kratos::_AddDataItemsToMenuItemD { group domNode } {
    set what [gid_groups_conds::open_conditions window_type]
    if { $what eq "none" } { set what window }
    set xpath [gid_groups_conds::nice_xpath $domNode]
    # W $xpath
    gid_groups_conds::open_conditions $what -select_xpath $xpath
}

proc Kratos::ResetModel { } {
    foreach layer [GiD_Info layers] {
        GiD_Process 'Layers Delete $layer Yes escape escape
    }
    foreach group [GiD_Groups list] {
        GiD_Groups delete $group
    }
}

proc Kratos::BeforeMeshGeneration {elementsize} {
    #catch {
        foreach group [GiD_Groups list] {
            GiD_Process Mescape Meshing MeshCriteria Mesh Lines {*}[GiD_EntitiesGroups get $group lines] escape escape Mescape
            GiD_Process Mescape Meshing MeshCriteria Mesh Surfaces {*}[GiD_EntitiesGroups get $group surfaces] escape escape 
        }
    #}
}

proc Kratos::PrintArray {a {pattern *}} {
    # ABSTRACT:
    # Print the content of array nicely
    
    upvar 1 $a array  
    if {![array exists array]} {
        error "\"$a\" isn't an array"
    }
    set maxl 0
    foreach name [lsort [array names array $pattern]] {
        if {[string length $name] > $maxl} {
            set maxl [string length $name]
        }
    }
    set maxl [expr {$maxl + [string length $a] + 2}]
    foreach name [lsort [array names array $pattern]] {
        set nameString [format %s(%s) $a $name]
        W "[format "%-*s = %s" $maxl $nameString $array($name)]"
    }
}

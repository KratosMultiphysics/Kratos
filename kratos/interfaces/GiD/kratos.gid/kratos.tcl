##################################################################################
#   This file is common for all Kratos Applications.
#   Do not change anything here unless it's strictly necessary.
##################################################################################

set ::kratos_debug 1 ;# could be 0,1,-1

##########################################################
#################### GiD Tcl events ######################
##########################################################
proc InitGIDProject { dir } {
    #uplevel #0 [list Kratos::InitGIDProject $dir]
    Kratos::InitGIDProject $dir
}

proc EndGIDProject {} {
    spdAux::DestroyWindow
    spdAux::EndRefreshTree
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
    gid_groups_conds::open_conditions check_default
    gid_groups_conds::open_conditions menu
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
    spdAux::reactiveApp
}

# Save GiD project files (save XML Tdom structure to spd file)
proc SaveGIDProject { filespd } {
    gid_groups_conds::save_spd_file $filespd
}

proc AfterTransformProblemType { filename oldproblemtype newproblemtype } {
    set spd_file [file join $filename.gid [file tail $filename].spd]
    return [gid_groups_conds::transform_problemtype $spd_file]
}

proc AfterWriteCalcFileGIDProject { filename errorflag } {
    set errcode [::write::writeEvent $filename]
    if {$errcode} {return "-cancel-"}
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
    array set kratos_private [ReadProblemtypeXml [file join $dir kratos.xml] Infoproblemtype {Name Version MinimumGiDVersion}]
    if { [GidUtils::VersionCmp $kratos_private(MinimumGiDVersion)] < 0 } {
        WarnWin [_ "Error: %s Interface requires GiD %s or later." $kratos_private(Name) $kratos_private(MinimumGiDVersion)]
    }    
    #append to auto_path only folders that must include tcl packages (loaded on demand with package require mechanism)
    if { [lsearch -exact $::auto_path [file join $dir scripts]] == -1 } {
        lappend ::auto_path [file join $dir scripts]
    }
    # JG Sources will be in a different proc
    foreach filename {Applications.tcl Writing.tcl spdAuxiliar.tcl } {
        uplevel 1 [list source [file join $dir scripts $filename]]
    }

    # JG Sources will be in a different proc
    foreach filename {Model.tcl Entity.tcl Parameter.tcl Topology.tcl Solver.tcl ConstitutiveLaw.tcl Condition.tcl Element.tcl SolutionStrategy.tcl Process.tcl} {
        uplevel 1 [list source [file join $dir scripts Model $filename]]
    }
    
    Kratos::ChangeMenus
     
    Kratos::load_gid_groups_conds
    #set HeaderBackground [$doc selectNodes string(Infoproblemtype/Program/HeaderBackground)]
    #gid_groups_conds::SetHeaderBackground $HeaderBackground
    gid_groups_conds::SetLibDir [file join $dir exec] 
    gid_groups_conds::begin_problemtype [file join $dir kratos_default.spd] [Kratos::GiveKratosDefaultsFile]

    spdAux::processIncludes
    spdAux::parseRoutes
    after 500 [list spdAux::CreateWindow $dir]
}

proc Kratos::ChangeMenus { } {
    GidChangeDataLabel "Data units" ""
    GidChangeDataLabel "Interval" ""
    GidChangeDataLabel "Conditions" ""
    GidChangeDataLabel "Materials" ""
    GidChangeDataLabel "Interval Data" ""
    GidChangeDataLabel "Problem Data" ""
    GidChangeDataLabel "Local axes" ""
    GidAddUserDataOptions "---" "" 3
    #GidAddUserDataOptions [_ "Groups"] [list gid_groups_conds::open_groups .gid window] 5
    GidAddUserDataOptions [_ "Data"] [list gid_groups_conds::open_conditions menu] 7
    GidAddUserDataOptions "---" "" 10
    GidAddUserDataOptionsMenu [_ "Local axes"] [list gid_groups_conds::local_axes_menu %W] 11
    GiDMenu::UpdateMenus
}

proc Kratos::load_gid_groups_conds {} {  
    package require customlib_extras ;#this require also customLib
    package require customlib_native_groups
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


package require tdom
package require snit
package require fulltktree

# if { 1 } {
    #     set f {C:\Documents and Settings\ramsan\Mis documentos\barres\newdatastyle\kk.db3}
    #     mylog::init -view_binding <Control-l> -file $f debug
    # } else {
    #     mylog::init ""
    # }

# to be inserted in tclIndex
proc boundary_conds {} {}

snit::widget boundary_conds {
    option -domdoc ""
    option -domnodes "" ; # if this is not NULL, -domdoc is not used
    option -edit_window_position down ;# up|down|full
    option -title ""
    option -title_active 1
    option -have_import_export 1
    option -modify_global_db 1
    option -can_delete_last_material 0
    option -external_update_handler ""
    option -auto_take_focus 1    
    
    variable searchframe_opener
    variable searchframe
    variable searchcombo
    variable last_searchstring
    variable searchstring
    
    variable tree
    variable vscrollbar
    variable hscrollbar
    
    variable unique_edit_id 1
    variable entryvalues
    variable entry_node_values
    variable entryvalues_units
    variable entryvalues_state
    variable entryvalues_function
    variable entryvalues_function_var
    variable entryvalues_haschanges 0
    variable entryvalues_mat_changes_stack ""
    
    variable transform_unit_value
    variable execute_proc_cancel 0    
    variable show_menubutton_about 1   
    
    delegate method * to tree
    delegate option * to tree
    delegate option -frame_bg as -bg to hull
    
    constructor { args } {                 
        gid_groups_conds::HeaderBackgroundColor
        gid_groups_conds::setHeaderBackground
        
        gid_groups_conds::setHeaderImage
        
        $win configure -frame_bg $::SystemBackground 
        
        if {[info exists gid_groups_conds::doc]} {
            set root [$gid_groups_conds::doc documentElement]
            set show_menubutton_about [[$gid_groups_conds::doc documentElement] @show_menubutton_about 1]           
        }
        
        #         set height [font metrics SystemFont -linespace]
        #         if {$height < 18} { set height 18 }   
        # Height depends on the entry combobox size
        if { [info command ::GiD_Process] eq "" } {
            set cbt [ttk::entry $win.size]
        } else {
            set cbt [ttk::entry $win.size -style Padding0.TEntry]
        }       
        set height [expr {[winfo reqheight $cbt] - 1}]        
        
        gid_groups_conds::SetHeaderName $win.f0
        set fwopen $win.f0      
        
        if { $show_menubutton_about } {
            menubutton $win.f0.b -image [cu::get_image submenu9x9] -menu $win.f0.b.m \
                -bg $::SystemBackground
            place $win.f0.b -relx 1 -x -5 -anchor ne -y 4        
            set m [menu $win.f0.b.m -tearoff 0]      
            
            $m add checkbutton -label [= "Inside graphical window"] -variable ::inside_graphical \
                -command [list gid_groups_conds::open_conditions toggle]
            
            #         $m add checkbutton -label [= "Always on top"] -variable ::always_on_top \
                #             -command [list draw_post::reopen_gui -toggle_transient 1]                
            
            $m add separator
            $m add command -label [= "About CustomLib..."] \
                -command [list gid_groups_conds::setCustomLib_window $win.f0]
        }
        
        set w_ .gid.central.boundaryconds
        
        if { ![winfo exists $w_] || [winfo class $w_] eq "Toplevel" } {
            set ::inside_graphical 0             
        } else {
            set ::inside_graphical 1
        }       
        
        #         set ::always_on_top [gid_groups_conds::get_preference draw_post::gui_transient 1]        
        
        gid_groups_conds::register_popup_help $fwopen [= "Toggle window or inline visualization"]
        bind $fwopen <1> [list draw_post::reopen_gui -toggle_frame_toplevel 1]
        
        button $win.f0.search -image [cu::get_image search-16] -bg $::SystemBackgroundHi -relief flat \
            -command [mymethod toggle_search]
        
        place $win.f0.search -relx 0 -x 0 -anchor nw -y -2
        gid_groups_conds::register_popup_help $win.f0.search [= "Search quickly in the data tree"]                                       
        
        grid columnconfigure $win.f0 0 -weight 1
        
        set fwopen $win.f0       
        
        if { !$gid_groups_conds::is_draw_post } {
            gid_groups_conds::register_popup_help $fwopen [= "Toggle window or inline visualization"]
            bind $fwopen <1> [list gid_groups_conds::open_conditions toggle]
        } else {
            gid_groups_conds::register_popup_help $fwopen [= "Open search box"]
            bind $fwopen <1> [mymethod toggle_search]
        }
        bind $fwopen <<Contextual>> [mymethod contextual_menu top "" %X %Y]
        
        
        #         set searchframe_opener [frame $win.f1 -bg $::SystemBackgroundSearchOpener \
            #                 -width 19 -height 4 -cursor hand2]
        
        #         set searchframe_opener [button $win.f1 -image [image create photo -file $gid_groups_conds::HeaderImage] \
            #                 -command [list gid_groups_conds::setheader_description_window $win.f0.f2] \
            #                 -bg $::SystemBackgroundHi -relief flat]
        
        #         bind $win.f1 <1> [mymethod toggle_search]
        #         bind $win.f1 <<Contextual>> [mymethod toggle_search]
        # 
        #         gid_groups_conds::register_popup_help $win.f1 [= "toggle search view"]
        
        set searchframe [ttk::labelframe $win.lb -text [= "Search..."]]
        #         frame $win.lb.f0 -bg $::SystemBackgroundSearchOpener -width 8 -height 4 \
            #             -cursor hand2
        #         $win.lb configure -labelwidget $win.lb.f0
        #         gid_groups_conds::register_popup_help $win.lb.f0 [= "toggle search view"]
        #         bind $win.lb.f0 <1> [mymethod toggle_search]
        #         bind $win.lb.f0 <<Contextual>> [mymethod toggle_search]
        
        button $win.lb.l1 -image [cu::get_image search-16] -relief flat \
            -command [mymethod clear_search] 
        
        install searchcombo as ttk::combobox $win.lb.cb1 -textvariable \
            [varname searchstring] -width 4
        
        gid_groups_conds::register_popup_help  $win.lb.l1 [= "Clear search"]
        gid_groups_conds::register_popup_help  $win.lb.cb1 [= "Search in tree"]
        
        set searchstring ""
        set last_searchstring ""
        grid $win.lb.l1 $win.lb.cb1 -sticky ew
        grid configure $win.lb.cb1 -padx "0 2" -pady "0 2"
        grid columnconfigure $win.lb 1 -weight 1
        
        bind $win.lb.cb1 <Return> "[mymethod search_tree full] ;break"
        set cmdT "[mymethod search_tree key] ;#"
        trace add variable [varname searchstring] write $cmdT
        bind $win.lb.cb1 <Destroy> [list trace remove variable \
                [varname searchstring] write $cmdT]
        
        frame $win.ft -bg [$win cget -frame_bg]
        install tree as treectrl $win.ft.t -highlightthickness 0 \
            -xscrollincrement 20 -showheader 0 -indent 2 \
            -font SystemFont -itemheight $height -selectmode extended -showroot no \
            -showrootbutton no -showbuttons yes -showlines yes \
            -scrollmargin 16 -xscrolldelay "500 50" -yscrolldelay "500 50" \
            -width 300 -borderwidth 0 -relief groove -doublebuffer window \
            -yscrollincrement $height
        #-backgroundimage [image create photo -file $::ProblemTypePriv(problemtypedir)/images/splash-tr.gif]
        
        bind $tree <Configure> [mymethod tree_changed_size]
        
        install vscrollbar as ttk::scrollbar $win.ft.sv -orient vertical -command [list $tree yview]
        install hscrollbar as ttk::scrollbar $win.ft.sh -orient horizontal -command [list $tree xview]
        
        #bind $win <FocusIn> [list after idle [list focus $tree]]
        
        set err [catch {
                $tree configure -openbuttonimage mac-collapse \
                    -closedbuttonimage mac-expand
            }]
        if { $err } {
            $tree configure -buttonimage [list mac-collapse open mac-expand !open ]
        }
        foreach i [list blockdata container condition groupList group \
                value disabled function functionVariable] {
            $tree state define $i
        }
        
        $tree column create
        $tree column configure 0  -text [= condition] -expand 1 
        #-itembackground {\#e0e8f0 {}} \
            
        catch { $tree configure -treecolumn 0 }
        bindtags $tree [list $tree FullTreeCtrl TreeCtrl [winfo toplevel $tree] all]
        set ::TreeCtrl::Priv(DirCnt,$tree) end
        
        $tree element create e_image image
        $tree element create e_text text -fill [list $::SystemHighlightText \
                {selected focus} $::SystemHighlightText {selected !focus disabled} \
                gray disabled darkblue {}] -lines 1       
        $tree element create e_rect rect -fill [list $::SystemHighlight {selected focus} \
                $::SystemHighlightNoFocus {selected !focus}] -showfocus yes -open we
        $tree element create f_image image
        
        set s [$tree style create s_text -orient horizontal]
        $tree style elements $s {e_rect e_image e_text f_image}
        $tree style layout $s e_image -padx {2 0} -expand ns
        $tree style layout $s e_text -padx {2 0} -squeeze x -expand ns
        $tree style layout $s e_rect -union [list e_text] \
            -iexpand nswe -ipadx 2
        $tree style layout $s f_image -padx {2 0} -expand ns
        
        $self createbindings
        
        grid $win.f0 -sticky new -row 0 -padx 2 -pady 1
        
        grid $win.lb -sticky new -padx "2 1" -pady "2 0" -row 2
        
        switch $options(-edit_window_position) {
            up {
                grid $win.ft -sticky nsew -row 4 -padx "2 1" -pady "0 2"
                grid rowconfigure $win 4 -weight 1
            }
            default {
                grid $win.ft -sticky nsew -row 3 -padx "2 1" -pady "0 2"
                grid rowconfigure $win 3 -weight 1
            }
        }
        grid $tree $win.ft.sv -sticky nsew -padx 1 -pady 1
        grid $win.ft.sh -sticky ew -row 5
        grid configure $tree -pady "2 0"
        grid columnconfigure $win.ft 0 -weight 1
        grid rowconfigure $win.ft 0 -weight 1
        grid columnconfigure $win 0 -weight 1
        
        $self configurelist $args
        
        autoscroll::autoscroll $win.ft.sv
        autoscroll::autoscroll $win.ft.sh
        
        set transform_unit_value 1
        
    }    
    onconfigure -domdoc {value} {
        set options(-domdoc) $value
        $self actualize
    }
    onconfigure -domnodes {value} {
        set options(-domnodes) $value
        
        array unset entryvalues
        array unset entry_node_values
        array unset entryvalues_units           
        
        $self actualize
    }
    onconfigure -title {value} {       
        set options(-title) $value       
        if { $value ne "" } {
            $win.f0 configure -text $value            
        }      
    }
    onconfigure -title_active {value} {
        set options(-title_active) $value
        if { $value } {
            bind $win.f0 <1> [list gid_groups_conds::open_conditions toggle]
            bind $win.f0 <<Contextual>> [list gid_groups_conds::open_conditions toggle]
            
            $win.f0 configure -cursor hand2
            
        } else {
            bind $win.f0 <1> ""
            bind $win.f0 <<Contextual>> ""
            
            $win.f0 configure -cursor ""              
            if { [winfo exists $win.f0.b] } {
                place forget $win.f0.b
            }
        }            
    }
    method add_external_frame { frame } {
        grid $frame -in $win -sticky nsew -row 5 -padx "2 1" -pady "0 2"
    }
    method info_frame_width {} {
        return [winfo width $win]
    }
    method clear_search  {} {
        if { [$searchcombo selection present] } {
            $searchcombo delete 0 end
        } else {
            $searchcombo selection range 0 end
            focus $searchcombo
        }
        $self search_tree full
    }
    method toggle_search {} {
        if { [grid info $searchframe] eq "" } {
            grid $searchframe
            focus $searchcombo
            #             grid remove $searchframe_opener
            $win.f0.search configure -relief sunken
        } else {
            grid remove $searchframe
            #             grid $searchframe_opener
            set searchstring ""
            $win.f0.search configure -relief flat
        }
    }
    method show_search {} {
        grid $searchframe
        focus $searchcombo
        #         grid remove $searchframe_opener
    }
    method hide_search {} {
        grid remove $searchframe
        #         grid $searchframe_opener
        if { $searchstring ne "" } { set searchstring "" }
    }
    method tree_changed_size {} {
        
        set height [winfo height $tree]
        if { $height < 2 } { return }
        
        set tree_prefsNode [gid_groups_conds::get_preference_node -local_global global \
                tree_preferences]
        
        set colors ""
        if { [$tree_prefsNode @background_color ""] ni [list "" auto] } {
            set c [$tree_prefsNode @background_color]
            $tree configure -background $c -backgroundimage ""
            lappend colors $c
        } elseif { [info command GiD_Set] eq "" } {
            # nothing
        } elseif { [GiD_Set BackColorType] == 0 } {
            $tree configure -background [GiDColorToTkColor [GiD_Set BackgroundColor]] \
                -backgroundimage ""
            lappend colors [GiDColorToTkColor [GiD_Set BackgroundColor]]
        } elseif { [GiD_Set BackColorSense] eq "r" } {
            $tree configure -background [GiD_Set BackColorTop] \
                -backgroundimage ""
            lappend colors [GiD_Set BackColorTop]
        } elseif { [GiD_Set BackColorSense] eq "l" } {
            $tree configure -background [GiD_Set BackColorBottom] \
                -backgroundimage ""
            lappend colors [GiD_Set BackColorBottom]
        } else {
            set img [$tree cget -backgroundimage]
            if { $img eq "" } {
                set img [image create photo]
                gid_groups_conds::add_to_internal_images_cache $img
            }
            lassign [list [GiD_Set BackColorTop] [GiD_Set BackColorBottom]] t b
            set m "#"
            set min ""
            foreach i [winfo rgb . $t] j [winfo rgb . $b] {
                append m [format %.4x [expr {round(0.5*($i+$j))}]]
            }
            switch [GiD_Set BackColorSense] {
                u { lassign [list $t $b] b t }
                ul { lassign [list $b $m] t b }
                dl { lassign [list $b $m] b t }
                ur { lassign [list $t $m] b t }
                dr { lassign [list $t $m] t b }
            }
            set h [expr {$height*2}]
            $img configure -width 2 -height $h
            $img put $t -to 0 0
            $img put $t -to 1 0
            $img put $b -to 0 [expr {$h-1}]
            $img put $b -to 1 [expr {$h-1}]
            cu::img::gradient -gradient_type double_vertical $img
            $tree configure -backgroundimage $img
            lappend colors $t $b
        }
        set minmod ""
        foreach c $colors {
            lassign [winfo rgb . $c] rc gc bc
            set mod [expr {sqrt($rc*$rc+$gc*$gc+$bc*$bc)}]
            if { $minmod eq "" || $mod < $minmod } {
                set minmod $mod
            }
        }
        
        if { [$tree_prefsNode @text_color ""] ni [list "" auto] } {
            set color [$tree_prefsNode @text_color]
            set lines_color $color
        } elseif { $minmod ne "" && $minmod < 20000 } {
            set color white
            set lines_color ""
        } else {
            set color darkblue
            set lines_color ""
        }
        $tree element configure e_text -fill [list $::SystemHighlightText \
                {selected focus} $::SystemHighlightText {selected !focus disabled} \
                gray disabled $color {}]
        if { $lines_color eq "" } {
            $tree configure -linecolor [lindex [$tree configure -linecolor] 3]
        } else {
            $tree configure -linecolor $lines_color
        }
    }
    variable search_tree_after ""
    # what can be: full and key
    method search_tree { what } {
        if { $searchstring eq $last_searchstring } { return }
        after cancel $search_tree_after
        set search_tree_after [after 200 [mymethod search_tree_do $what]]
    }
    method search_tree_do { what } {
        
        after cancel $search_tree_after
        set search_tree_after ""
        
        set str [string trim $searchstring]
        if { $str eq "" } {
            foreach item [$self selection get] {
                $self _store_tree_state_item $item
                foreach itemP [$self item ancestors $item] {
                    if { $itemP == 0 } { continue }
                    $self _store_tree_state_item $itemP
                }
            }
            $self recover_tree_state 0
            set last_searchstring $str
            return
        }
        if { $last_searchstring eq "" } {
            $self store_tree_state all 0
        }
        if { [string length $searchstring] == 1 && $last_searchstring eq "" } { return }
        
        $tree item collapse 0 -recurse
        set isexpanded 0
        foreach i [$tree item children 0] {
            set isex [$self search_expand $i [gid_groups_conds::fcompare_str \
                        $searchstring]]
            if { $isex } { set isexpanded 1 }
        }
        set last_searchstring $str
        
        set values [$searchcombo cget -values]
        
        if { [set ipos [lsearch -exact $values $str]] != -1 } {
            set values [lreplace $values $ipos $ipos]
            set replaced 1
        } else { set replaced 0 }
        
        if { $what eq "full" || $replaced } {
            set values [linsert $values 0 $str]
        } elseif { $isexpanded && [string first $str [lindex $values 0]] == -1} {
            while { [string first [lindex $values 0] $str] != -1 } {
                set values [lrange $values 1 end]
            }
            set values [linsert $values 0 $str]
        }
        set values [lrange $values 0 14]
        $searchcombo configure -values $values
    }
    method search_expand { id text { visible "" } } {
        set text_tree [gid_groups_conds::fcompare_str [$tree item text $id 0]]
        if { [string first $text $text_tree] != -1 } {
            set expand 1
            $tree item configure $id -visible 1
            set visible 1
        } else {
            set expand 0
            if { $visible != 1 } {
                $tree item configure $id -visible 0
            }
        }
        set expand_this 0
        $self expand_node $id
        foreach i [$tree item children $id] {
            set isex [$self search_expand $i $text $visible]
            if { $isex } { set expand_this 1 }
        }
        if { $expand_this } {
            if { $expand == 0 } {
                $tree item state set $id disabled
                set expand 1
            } else {
                $tree item state set $id !disabled
            }
            $tree item expand $id
            $tree item configure $id -visible 1
        } else {
            $tree item state set $id !disabled
        }
        return $expand
    }
    method expand_node { id } {
        if { [$tree item numchildren $id] > 0 } { return }       
        set domNode [$tree item element cget $id 0 e_text -data]
        gid_groups_conds::uncompress_subtree $domNode
        
        $self _add_to_tree_from_dom [$domNode childNodes] $id \
            [get_domnode_attribute $domNode state normal]
    }
    method collapse_node { id } {
        if { $id == 0 } { return }        
        set domNode [$tree item element cget $id 0 e_text -data]
        set ret [gid_groups_conds::compress_subtree $domNode]
        if { $ret } {
            foreach item [$tree item children $id] {
                $tree item delete $item
            }
        }
    }
    method set_item_state { domNode id { also_selection 1 } } {
        
        set isopen 0
        $tree item configure $id -visible 1
        set tree_stateList [split [$domNode @tree_state ""] ,]
        if { [lsearch -exact $tree_stateList open] == -1 } {
            lappend tree_stateList close
        }
        if { [lsearch -exact $tree_stateList disabled] == -1 } {
            lappend tree_stateList !disabled
        }
        if { [lsearch -exact $tree_stateList selected] == -1 } {
            lappend tree_stateList !selected
        }
        foreach state $tree_stateList {
            switch $state {
                open {
                    if { ![$tree item state get $id open] } {
                        $tree item expand $id
                    }
                    set isopen 1
                }
                close {
                    if { [$tree item state get $id open] } {
                        $tree item collapse $id
                    }
                }
                active { $tree activate $id }
                selected { if { $also_selection } { $tree selection add $id } }
                !selected {
                    if { $also_selection } {
                        $tree selection clear $id
                    }
                }
                disabled { $tree item state set $id disabled }
                !disabled { $tree item state set $id !disabled }
            }
        }
        return $isopen
    }
    method _give_icon { node } {
        if {[get_domnode_attribute $node icon ""] eq "-"} {
            set image ""
            return $image
        } elseif { [get_domnode_attribute $node icon ""] ne "" } {                      
            set err [catch {gid_groups_conds::giveimageL [get_domnode_attribute $node icon ""] 
                } image]        
            if { $err || $image == "" } {
                set image [gid_groups_conds::giveimage contents2-16]
            }
        } elseif { [get_domnode_attribute $node icon_end ""] ne "" } {
            if { [get_domnode_attribute $node icon_end ""] ne "" } {
                if {[get_domnode_attribute $node icon_end ""] eq ""} {
                    set image ""
                    return $image
                }
                set err [catch {gid_groups_conds::giveimageL [get_domnode_attribute $node icon_end ""] 
                    } image]
                if { $err || $image == "" } {
                    set image [gid_groups_conds::giveimage contents2-16]
                }
            }
        } else {
            set image [gid_groups_conds::giveimage contents2-16]
        }        
        return $image
    }
    method _add_to_tree_from_dom { childNodes tree_parent parent_state } {
        
        foreach childNode $childNodes {
            if { [$childNode nodeType] eq "COMMENT_NODE" } { continue }
            set state [get_domnode_attribute $childNode state normal]
            if { $state eq "hidden" } { continue }
            if { $parent_state eq "disabled" } { set state disabled }
            if { $state eq "disabled_on_void" && [$childNode @v ""] eq "" } {
                set state disabled
            }
            
            if { $state eq "hidden_this" } { 
                set sub_childNodes [$childNode childNode]
                foreach sub_childNode $sub_childNodes {
                    $self _add_to_tree_from_dom $sub_childNode $tree_parent $parent_state
                }
                continue
            }          
            
            switch -- [$childNode nodeName] {
                container {                                    
                    if { [$childNode selectNodes {name(..)}] eq "condition" } {
                        continue
                    }                
                    set image [$self _give_icon $childNode]                          
                    
                    set id [$self add $tree_parent [get_domnode_attribute $childNode pn] \
                            $image container \
                            $childNode]
                    set xp {.//condition/group|.//condition/groupList}
                    set boldfont 0
                    if { [llength [$childNode selectNodes $xp]] } {
                        $tree item element configure $id 0 e_text -font \
                            SystemBoldFont
                        set boldfont 1
                    } else {
                        $tree item element configure $id 0 e_text -font ""
                    } 
                    if {!$boldfont} {
                        set xp {.//blockdata}                
                        if { [llength [$childNode selectNodes $xp]] } {
                            set icount 0
                            foreach iNode [$childNode selectNodes $xp] {
                                if { [get_domnode_attribute $iNode state normal] ne "hidden" && \
                                    [$iNode @sequence 0] == "1" && [$iNode @isbold 0] == "1" } {                                                                
                                    incr icount                                 
                                }
                            }                       
                            if {$icount > 0} {
                                $tree item element configure $id 0 e_text \
                                    -font SystemBoldFont
                            } else {                            
                                $tree item element configure $id 0 e_text -font ""                            
                            }
                        }  
                    }
                    
                    set open [$self set_item_state $childNode $id]
                    set xp "container|blockdata|condition|group|value"
                    if { $open } {
                        $self _add_to_tree_from_dom [$childNode childNodes] $id \
                            $state                        
                    } elseif { [llength [$childNode selectNodes $xp]] } {
                        set is_button [$self is_button [$childNode selectNodes $xp]]
                        if { $is_button } {
                            $tree item configure $id -button 1
                        }
                    }                   
                    set title [get_domnode_attribute $childNode title ""]
                    if {$title ne ""} {
                        set font [concat [font actual SystemBoldFont] [list -underline 1]]                      
                        $tree item element configure $id 0 e_text -fill black -font [list $font]
                    }
                }
                blockdata {
                    if { [$childNode @n ""] eq "Internal data" } { continue }
                    set name [get_domnode_attribute $childNode name]
                    if { $name eq "" } {
                        set name [get_domnode_attribute $childNode pn]
                    }
                    set image [$self _give_icon $childNode]      
                    set id [$self add $tree_parent $name \
                            $image blockdata \
                            $childNode]
                    if { [$childNode @sequence 0] == 1 && [$childNode @active 1] == 1 && \
                        [$childNode @isbold 0] == 1 } {
                        $tree item element configure $id 0 e_text \
                            -font SystemBoldFont
                    }
                    set open [$self set_item_state $childNode $id]
                    if { $open } {
                        $self _add_to_tree_from_dom [$childNode childNodes] $id \
                            $state
                    } elseif { [llength [$childNode childNodes]] } {
                        set is_button [$self is_button [$childNode childNodes]]
                        if { $is_button } {
                            $tree item configure $id -button 1
                        }
                    }
                }
                condition {
                    set id [$self add $tree_parent [get_domnode_attribute $childNode pn] \
                            [gid_groups_conds::giveimage filenew-16] condition \
                            $childNode]
                    if { [llength [$childNode selectNodes group|groupList]] } {
                        $tree item element configure $id 0 e_text -font \
                            SystemBoldFont
                    }
                    set open [$self set_item_state $childNode $id]
                    if { $open } {
                        $self _add_to_tree_from_dom [$childNode childNodes] $id \
                            $state
                    } elseif { [llength [$childNode selectNodes group|groupList]] } {
                        set is_button [$self is_button [$childNode selectNodes group|groupList]]
                        if { $is_button } {
                            $tree item configure $id -button 1
                        }
                    }
                }
                groupList {
                    set id [$self add $tree_parent [= "groupList"] \
                            [gid_groups_conds::giveimage kcmdf-16] groupList \
                            $childNode]
                    set open [$self set_item_state $childNode $id]
                    if { $open } {
                        $self _add_to_tree_from_dom [$childNode childNodes] $id \
                            $state
                    } elseif { [llength [$childNode childNodes]] } {
                        set is_button [$self is_button [$childNode childNodes]]
                        if { $is_button } {
                            $tree item configure $id -button 1
                        }
                    }
                }
                group {
                    set id [$self add $tree_parent [= "group: %s" [$childNode @n]] \
                            [gid_groups_conds::giveimage kcmdf-16] group \
                            $childNode]
                    set open [$self set_item_state $childNode $id]
                    if { $open } {
                        $self _add_to_tree_from_dom [$childNode childNodes] $id \
                            $state
                    } elseif { [llength [$childNode childNodes]] } {
                        set is_button [$self is_button [$childNode childNodes]]
                        if { $is_button } {
                            $tree item configure $id -button 1
                        }
                    }
                }
                value {
                    if { [$childNode selectNodes {name(..)}] eq "condition" } { continue }                                                          
                    set v [gid_groups_conds::give_printable_value $childNode]
                    set funcNode [$childNode selectNodes function]
                    if { $funcNode ne "" } {
                        foreach fvarNode [$funcNode selectNodes functionVariable] {
                            set pn [$fvarNode @pn]
                            if { [string length $pn] > 20 } {
                                set pn [string range $pn 0 16]...
                            }
                            lappend nameList $pn
                        }
                        set v f([join $nameList ,])...
                    }
                    set name [get_domnode_attribute $childNode pn]
                    if { $name eq "" } {
                        set name [get_domnode_attribute $childNode name]
                    }
                    set txt "$name: $v"
                    
                    set id [$self add $tree_parent $txt \
                            [gid_groups_conds::giveimage editclear-16] value \
                            $childNode]
                    set open [$self set_item_state $childNode $id]
                    if { $open } {
                        $self _add_to_tree_from_dom [$childNode childNodes] $id \
                            $state
                    } elseif { $funcNode ne "" } {
                        $tree item configure $id -button 1
                    }
                }
                function {
                    set image [$self _give_icon $childNode]      
                    set id [$self add $tree_parent [= function] \
                            $image function $childNode]
                    set open [$self set_item_state $childNode $id]
                    if { $open } {
                        $self _add_to_tree_from_dom [$childNode childNodes] $id \
                            $state
                    } else {
                        $tree item configure $id -button 1
                    }
                }
                functionVariable {
                    set pn [get_domnode_attribute $childNode pn]
                    if { [string length $pn] > 20 } {
                        set pn [string range $pn 0 16]...
                    }
                    set image [$self _give_icon $childNode]      
                    set id [$self add $tree_parent $pn \
                            $image functionVariable \
                            $childNode]
                    set open [$self set_item_state $childNode $id]
                    if { $open } {
                        $self _add_to_tree_from_dom [$childNode childNodes] $id \
                            $state
                    } else {
                        $tree item configure $id -button 1
                    }
                }
            }
            if { $state eq "disabled" && [info exists id] } {
                $tree item state set $id $state
            }          
            if { [get_domnode_attribute $childNode icon ""] ne "" } {
                set err [catch { gid_groups_conds::giveimageL [get_domnode_attribute $childNode icon ""] } image]                
                if { !$err && [info exists id] } {
                    if { $image == "" } {
                        set image [gid_groups_conds::giveimage contents2-16]
                    }
                    $tree item element configure $id 0 e_image -image $image                                     
                }
                if {[get_domnode_attribute $childNode icon_end ""] ne ""} {
                    set err [catch { gid_groups_conds::giveimageL [get_domnode_attribute $childNode icon_end ""] } f_image]                    
                    if { !$err && [info exists id] } {
                        if { $f_image == "" } {
                            set f_image [gid_groups_conds::giveimage contents2-16]
                        }
                        $tree item element configure $id 0 f_image -image $f_image   
                    }
                }
            }
        }
    }
    method is_button { domNodeList } {
        set is_button 0
        foreach iNode $domNodeList {
            if { [$iNode nodeName] eq "#comment" } { 
                # nothing
                continue 
            } 
            set state [get_domnode_attribute $iNode state normal]
            if { $state eq "normal" || $state eq "disabled" } { set is_button 1; break }
        }
        return $is_button
    }
    method _actualize_attributes { childNode id } {
        
        set state [get_domnode_attribute $childNode state normal]
        if { $state eq "hidden" } { return }
        if { $state eq "disabled_on_void" && [$childNode @v ""] eq "" } {
            set state disabled
        }
        switch -- [$childNode nodeName] {
            container {
                set xp {.//condition/group|.//condition/groupList}
                if { [llength [$childNode selectNodes $xp]] } {
                    $tree item element configure $id 0 e_text -font \
                        SystemBoldFont
                } else {
                    $tree item element configure $id 0 e_text -font ""
                }
            }
            blockdata {
                if { [$childNode @sequence 0] == 1 && [$childNode @active 1] == 1 } {
                    $tree item element configure $id 0 e_text \
                        -font SystemBoldFont
                } else {
                    $tree item element configure $id 0 e_text \
                        -font ""
                }
            }
            condition {
                if { [llength [$childNode selectNodes group|groupList]] } {
                    $tree item element configure $id 0 e_text -font \
                        SystemBoldFont
                } else {
                    $tree item element configure $id 0 e_text -font ""
                }
            }
            value {
                if { [$childNode selectNodes {name(..)}] eq "condition" } { continue }
                set v [gid_groups_conds::give_printable_value $childNode]
                set funcNode [$childNode selectNodes function]
                if { $funcNode ne "" } {
                    foreach fvarNode [$funcNode selectNodes functionVariable] {
                        lappend nameList [$fvarNode @pn]
                    }
                    set v f([join $nameList ,])...
                }
                set name [get_domnode_attribute $childNode pn]
                if { $name eq "" } {
                    set name [get_domnode_attribute $childNode name]
                }
                set txt "$name: $v"
                
                $tree item element configure $id 0 e_text -text $txt
            }
        }
        
        if { $state eq "disabled" } {
            $tree item state set $id disabled
        } elseif { $state eq "normal" } {
            $tree item state set $id !disabled
        }
        if { [get_domnode_attribute $childNode icon ""] ne "" } {
            set err [catch {gid_groups_conds::giveimageL [get_domnode_attribute $childNode icon ""] 
                } image]           
            if { !$err && [info exists id] } {
                if { $image == "" } {
                    set image [gid_groups_conds::giveimage contents2-16]
                }
                $tree item element configure $id 0 e_image -image $image              
            }
            if { [get_domnode_attribute $childNode icon_end] ne "" } {
                set err [catch {gid_groups_conds::giveimageL [get_domnode_attribute $childNode icon_end ""] 
                    } f_image]                
                if { !$err && [info exists id] } {
                    if { $f_image == "" } {
                        set f_image [gid_groups_conds::giveimage contents2-16]
                    }
                    $tree item element configure $id 0 f_image -image $f_image
                }
            }
        }
    }
    method actualize { { id 0 } } {
        return [$self actualize_do -do_focus $options(-auto_take_focus) $id]
    }
    method actualize_do { args } {
        set optional {
            { -do_focus boolean 1 }
        }   
        set compulsory "id"
        parse_args $optional $compulsory $args
        
        if { $options(-domdoc) ne "" } {
            set root [$options(-domdoc) documentElement]
            set vs [$root selectNodes string(display_options/@view_conditions_search)]
            if { $vs == 1 } {
                $self show_search
            } else {
                $self hide_search
            }
            set v [$root selectNodes string(display_options/@conditions_search_values)]
            $searchcombo configure -values [split $v ,]
        } elseif { [grid info $searchframe] ne "" } {
            $self toggle_search
        }
        
        set xview [lindex [$tree xview] 0]
        set yview [lindex [$tree yview] 0]
        
        $self store_tree_state all $id
        foreach i [$tree item children $id] {
            $tree item delete $i
        }
        if { $options(-domdoc) eq "" && $options(-domnodes) eq "" } { return }
        
        if { $id == 0 } {
            if { $options(-domdoc) ne "" } {
                set domNode [$options(-domdoc) documentElement]
                set childNodes [$domNode childNodes]
                set state [$domNode @state normal]
            } else {
                set childNodes $options(-domnodes)
                set state normal
            }
        } else {           
            set domNode [$tree item element cget $id 0 e_text -data]
            $self _actualize_attributes $domNode $id
            gid_groups_conds::uncompress_subtree $domNode
            set childNodes [$domNode childNodes]
            set state [$domNode @state normal]
        }
        $self _add_to_tree_from_dom $childNodes $id $state
        
        #$self recover_tree_state
        
        if { $id != 0 } {
            if { [$tree item children $id] ne "" } {
                $tree item configure $id -button 1
            } else {
                $tree item configure $id -button 0
            }
        }
        $tree xview moveto $xview
        $tree yview moveto $yview
        if { [$tree index 1] == 1 } {
            if { [$tree index active] == 0 } { $tree activate 1 }
            if { ![$tree selection count] } { $tree selection add 1 }
        }
        if { $do_focus } {
            focus $win
        }
    }
    method actualize_domNode { args } {
        
        set nodesList [lsort -unique $args]
        if { [lindex $nodesList 0] == 0 } {
            $self actualize
            return
        }
        set itemList ""
        foreach domNode $nodesList {
            set unc [gid_groups_conds::uncompress_subtree $domNode]
            set item 0
            foreach node [lrange [$domNode selectNodes ancestor-or-self::*] 1 end] {
                set state [get_domnode_attribute $node state normal]
                if { $state eq "hidden_this" } { continue }
                set found 0
                foreach i [$tree item children $item] {                    
                    set d [$tree item element cget $i 0 e_text -data]
                    if { $d eq $node } {
                        set item $i
                        set found 1
                        break
                    }
                }
                if { !$found } {
                    set item ""
                    break
                }
            }
            if { $item ne "" } { lappend itemList $item }
            if { $unc } { gid_groups_conds::compress_subtree $domNode }
        }
        set itemList [lsort -decreasing -unique -integer $itemList]
        foreach item $itemList {
            set found 0
            foreach i [$tree item ancestors $item] {
                if { [lsearch -decreasing  -sorted -integer $itemList $i] != -1 } {
                    set found 1
                    break
                }
            }
            if { !$found } {
                $self actualize $item
            }
        }
    }
    
    delegate method _selection to tree as selection
    method selection { args } {
        #         mylog::debug "selection $args"
        return [$tree selection {*}$args]
    }
    
    method select_domNode_xpath { xpath } {
        
        if { $options(-domdoc) eq "" } { return }
        set root [$options(-domdoc) documentElement]
        set domNode [$root selectNodes $xpath]
        if { $domNode eq "" } {
            error "xpath '$xpath' returned no nodes"
        }
        return [$self select_domNode $domNode]
    }
    method select_domNode { domNode } {
        
        $tree item collapse all
        
        set ancestors [$domNode selectNodes ancestor-or-self::*]
        foreach id [$tree item children 0] {           
            set d [$tree item element cget $id 0 e_text -data]
            if { [set ipos [lsearch -exact $ancestors $d]] != -1 } {
                set id_parent $id
                $self expand_node $id
                $tree item expand $id
                set ancestors [lrange $ancestors [expr {$ipos+1}] end]
            }
        }
        set id $id_parent
        foreach node $ancestors {
            foreach id [$tree item children $id_parent] {                
                set d [$tree item element cget $id 0 e_text -data]
                if { $d eq $node } {
                    $self expand_node $id
                    $tree item expand $id
                    set id_parent $id
                    break
                }
            }
        }
        $tree selection clear all
        $tree selection add $id
        $tree activate $id
        $tree see $id
        focus $tree
        return $id
    }
    method check_if_closed_and_clear { domNode } {
        set ancestors [$domNode selectNodes ancestor-or-self::*]
        foreach id [$tree item children 0] {            
            set d [$tree item element cget $id 0 e_text -data]
            if { [set ipos [lsearch -exact $ancestors $d]] != -1 } {
                set id_parent $id
                set ancestors [lrange $ancestors [expr {$ipos+1}] end]
                break
            }
        }
        set id $id_parent
        foreach node $ancestors {
            set found 0
            foreach id [$tree item children $id_parent] {               
                set d [$tree item element cget $id 0 e_text -data]
                if { $d eq $node } {
                    set found 1
                    set id_parent $id
                    break
                }
            }
            if { !$found } { return 1 }
        }
        if { [$tree item state get $id open] } {
            return 0
        } else {
            foreach child_id [$tree item children $id] {
                $tree item delete $child_id
            }
            return 1
        }
    }
    method _store_tree_state_item { id } {
        set ss ""
        foreach state [list active open selected] {
            if { [$tree item state get $id $state] } {
                lappend ss $state
            }
        }
        if { [lsearch -exact $ss open] == -1 } { lappend ss close }   
        set domNode [$tree item element cget $id 0 e_text -data]
        if { [info command $domNode] ne "" } {
            if { ![llength $ss] } {
                if { [$domNode hasAttribute tree_state] } {
                    $domNode removeAttribute tree_state
                }
            } else {
                $domNode setAttribute tree_state [join $ss ,]
            }
        }
        foreach child_id [$tree item children $id] {
            $self _store_tree_state_item $child_id
        }
    }
    # what: all display_options
    method store_tree_state { { what all } { id 0 } } {
        if { [$tree index 1] != 1 } { return }
        if { $options(-domdoc) eq "" } { return }
        
        if { $what eq "all" } {
            if { $id != 0 } {
                $self _store_tree_state_item $id
            }
            foreach child_id [$tree item children $id] {
                $self _store_tree_state_item $child_id
            }
        }
        set root [$options(-domdoc) documentElement]
        set width [gid_groups_conds::open_conditions info_frame_width]
        
        set dispNode [$root selectNodes display_options]
        if { $dispNode eq "" } {
            set dispNode [$root appendChildTag display_options]
        }
        if { $width ne "" && $width > 1 } {
            $dispNode setAttribute frame_width $width
            $dispNode setAttribute is_frame_open 1
        } else {
            $dispNode setAttribute is_frame_open 0
        }
        if { [grid info $searchframe] eq "" } {
            $dispNode setAttribute view_conditions_search 0
        } else {
            $dispNode setAttribute view_conditions_search 1
        }
        # if one search value has a comma: problems
        $dispNode setAttribute conditions_search_values \
            [join [$searchcombo cget -values] ,]
    }
    method store_tree_state_close_frame {} {
        set root [$options(-domdoc) documentElement]
        set dispNode [$root selectNodes display_options]
        $dispNode setAttribute is_frame_open 0
    }
    method recover_tree_state { { also_selection 1 } } {
        if { [$tree index 1] != 1 } { return }
        
        foreach id [$tree range 1 last] {            
            set domNode [$tree item element cget $id 0 e_text -data]
            if { $domNode eq "" } { continue }
            $self set_item_state $domNode $id $also_selection
        }
    }
    method execbindings { what args } {
        switch $what {
            <<Contextual>> {
                foreach "x y X Y" $args break
                foreach "type id" [list "" ""] break
                foreach "type id" [$tree identify $x $y] break
                if { $type != "item" || $id == "" } { return -code break }
                focus $tree
                if { ![$tree selection includes $id] } {
                    $tree selection clear all
                    $tree selection add $id
                }
                $tree activate $id
                $self contextual_menu $x $y $X $Y
            }
            <App> {
                foreach "x y - -" [$tree item bbox active 0 e_text] break
                set X [expr {$x+[winfo rootx $tree]}]
                set Y [expr {$y+[winfo rooty $tree]}]
                $self contextual_menu $x $y $X $Y
            }
        }
    }
    method createbindings {} {
        bind $tree <KeyPress-Left> {
            #TreeCtrl::SetActiveItem %W [TreeCtrl::LeftRight %W -1]
            if { [%W item numchildren [%W index active]] } {
                %W collapse [%W index active]
            } elseif { [%W index "active parent"] != 0 } {
                %W activate "active parent"
                %W selection clear all
                %W selection add active
            }
        }
        bind $tree <KeyPress-Right> {
            #TreeCtrl::SetActiveItem %W [TreeCtrl::LeftRight %W 1]
            %W expand [%W index active]
        }
        #         bind TreeCtrl <Double-ButtonPress-1> {+
            #             set id [%W identify %x %y]
            #             if {[lindex $id 0] eq "item"} {
                #                 %W toggle [%W index active]
                #             }
            #         }
        #
        
        bind $tree <Motion> [mymethod manage_motion %x %y]
        bind $tree <Leave> [mymethod manage_motion %x %y]
        $tree notify bind $vscrollbar <Scroll-y> { %W set %l %u }
        bind $vscrollbar <ButtonPress-1> [list focus $tree]
        
        $tree notify bind $hscrollbar <Scroll-x> { %W set %l %u }
        bind $hscrollbar <ButtonPress-1> [list focus $tree]
        
        #         $tree notify install event Header
        #         $tree notify install detail Header invoke
        # 
        #         $tree notify install event Drag
        #         $tree notify install detail Drag begin
        #         $tree notify install detail Drag end
        #         $tree notify install detail Drag receive
        # 
        $tree notify install event Edit
        $tree notify install detail Edit accept
        
        #$self notify bind $self <Header-invoke> [mymethod headerinvoke %C]
        
        #         $tree notify bind DontDelete <Selection> \
            #             [mymethod cancel_cond_or_blockdata_edit]
        
        #         $tree notify bind DontDelete <Selection> [mymethod select]
        #         $tree notify bind DontDelete <Selection> {
            #             if {%c == 1} {
                #                 set selection [%T selection get]
                #                 #xmlwidget::DisplayNode %T [lindex $selection 0]
                #             }
            #         }
        
        bind $tree <KeyPress> [mymethod search_text %A]
        bind $tree <Return> "[mymethod execute_select return]"
        #bind $tree <Return> [mymethod add_group]
        bind $tree <ButtonPress-1> [mymethod popup_help_cancel]
        bind $tree <ButtonRelease-1> "TreeCtrl::Release1 %W %x %y
            [mymethod execute_select %x %y] ; break"
        
        bind $tree <Control-ButtonRelease-1> "TreeCtrl::Release1 %W %x %y
            break"
        bind $tree <Shift-ButtonRelease-1> "TreeCtrl::Release1 %W %x %y
            break"
        
        bind $tree <<Contextual>> [mymethod execbindings <<Contextual>> \
                %x %y %X %Y]
        catch {
            bind $tree <App> [mymethod execbindings <App>]
        }
        
        bind $tree <F2> "[mymethod edit_value_or_blockdata]; break"
        
        $tree notify bind $tree <Edit-accept> [mymethod edit_value_or_blockdata_ok \
                %I %t]
        
        #         $tree notify bind $tree <Drag-receive> {
            #             set self [winfo parent %W]
            #             $self end_drop %I %l %x %y
            #         }
        
        TreeCtrl::SetSensitive $tree [list [list 0 s_text e_image e_text f_image]]
        TreeCtrl::SetDragImage $tree ""
        TreeCtrl::SetEditable $tree ""
        
        #         foreach i [list edit sensitive dragimage] {
            #             set ::TreeCtrl::Priv($i,$tree) ""
            #         }
        #         lappend ::TreeCtrl::Priv(sensitive,$tree)  \
            #             [list 0 s_text e_image e_text]
        #         lappend ::TreeCtrl::Priv(dragimage,$tree) \
            #             [list 0 s_text e_text] [list 1 s_image e_image]
        #        bindtags $tree [list $tree TreeCtrlFileList TreeCtrl [winfo toplevel $tree] all]
        
        #         $tree notify bind $tree <Drag-receive> {
            #             set w [winfo parent %T]
            #             $w move_group %I
            #         }
        
        bind $tree <Delete> [mymethod delete_cond_group_or_blockdata]
        
        #bind $tree <<Cut>> [list xmlwidget::CutOrCopy $tree cut]
        #bind $tree <<Copy>> [list xmlwidget::CutOrCopy $tree copy]
        
        $tree notify bind $tree <Expand-before> [mymethod expand_node %I]
        $tree notify bind $tree <Collapse-after> [mymethod collapse_node %I]
    }
    variable manage_motion_id
    variable manage_motion_active 1
    method manage_motion { x y } {
        if { !$manage_motion_active } { return }
        
        if { [info exists manage_motion_id] } {
            after cancel $manage_motion_id
            unset -nocomplain manage_motion_id
        }
        set identify [$tree identify $x $y]
        if { [lindex $identify 0] ne "item" } { return }
        set id [lindex $identify 1]
        set manage_motion_id [after 600 [list catch [mymethod popup_help $id]]]
    }
    method popup_help_deactivate {} {
        if { [info exists manage_motion_id] } {
            after cancel $manage_motion_id
            unset -nocomplain manage_motion_id
        }
        set manage_motion_active 0
        gid_groups_conds::register_popup_help_deactivate
    }
    method popup_help_reactivate {} {
        set manage_motion_active 1
        gid_groups_conds::register_popup_help_reactivate
    }
    method popup_help_cancel {} {
        if { [info exists manage_motion_id] } {
            after cancel $manage_motion_id
            unset -nocomplain manage_motion_id
        }
    }
    method popup_help { id } {
        unset -nocomplain manage_motion_id
        
        if { [$tree item id $id] eq "" } { return }
        foreach "x1 y1 x2 y2" [$tree item bbox $id 0 e_text] break
        if { ![info exists x1] } { return }
        
        set domNode [$tree item element cget $id 0 e_text -data]
        
        set text [$self _give_node_help $domNode]              
        
        if { $text eq "" } {
            if { $x1 >= 0 && $x2 < [winfo width $tree] } { return }
            set err [catch { $tree item text $id 0 } text]
            if { $err } { return }
        }
        set x [expr {[winfo pointerx $tree]-[winfo rootx $tree]}]
        set y [expr {[winfo pointery $tree]-[winfo rooty $tree]}]
        
        if { $x < $x1 || $x > $x2 || $y < $y1 || $y > $y2 } { return }
        
        set x [expr {[winfo pointerx $tree]+15}]
        set y [expr {[winfo pointery $tree]+10}]
        
        set manage_motion_active 0
        set w [gid_groups_conds::popup_help $tree $text $x $y]
        if { ![winfo exists $w] } { return }
        bind $w <Destroy> [list set [varname manage_motion_active] 1]
    }
    method _give_node_help_image { domNode } {
        set help_image [$domNode @help_image ""]
        while { $help_image eq "" } {
            set domNode [$domNode parentNode]
            if { $domNode eq "" } { break }
            set help_image [$domNode @help_image ""]
            if { [$domNode nodeName] ne "container" } { break }
        }
        return $help_image
    }
    method _give_node_help { domNode } {
        set help [$domNode @help ""]
        while { $help eq "" } {
            set domNode [$domNode parentNode]
            if { $domNode eq "" } { break }
            set help [$domNode @help ""]
            if { [$domNode nodeName] ne "container" } { break }
        }
        return [= [subst -nocommands -novariables $help]]
    }
    method search_text { char } {
        return
        
        if { $char eq "\t" } { return }
        if { [$tree index "active bottom"] eq "" } { return }
        if { [string is wordchar -strict $char] || [string is punct -strict $char] \
            || [string is space -strict $char] } {
            if { ![info exists searchstring] || [string index $searchstring end] != $char } {
                append searchstring $char
            }
            if { [info exists searchstring_reached_end] && $searchstring_reached_end ne "" \
                && $searchstring_reached_end eq $searchstring } {
                set id "first visible"
            } elseif { [$tree compare active == "active bottom"] } {
                set id "first visible"
            } else { set id "active below" }
            set found 0
            while { $id != "" } {
                set txt [$tree item text $id 0]
                if { [string match -nocase $searchstring* $txt] } {
                    set found 1
                    break
                }
                set id [$tree index "$id below"]
            }
            if { !$found } {
                bell
                set searchstring_reached_end $searchstring
                set searchstring ""
                after 300 [list set [varname searchstring_reached_end] ""]
            } else {
                $tree activate $id
                $tree see $id
                $tree selection clear all
                $tree selection add active
                after 300 [list set [varname searchstring] ""]
            }
        }
    }
    variable execute_select_pressed ""
    variable execute_select_handler ""
    method execute_select { { x "" } { y "" } } {
        
        set id [lindex [$tree selection get] 0]
        if { $id eq "" } { return }
        if { $x eq "return" } {
            set identify "return"
        } elseif {$x ne "" } {
            set identify [$tree identify $x $y]
            if { $x ne "" && ([lindex $identify 0] ne "item" || \
                [lindex $identify 2] ne "column") } { return }
        } else { set identify "" }
        
        if { $execute_select_pressed == $id || $identify eq "return" } {
            set execute_select_pressed ""
            
            # this is double click or Return
            
            if { $identify ne "return" } {
                set ident [$tree identify $x $y]
                set id [lindex $ident 1]
                if { [lindex $ident 0] ne "item" || $id eq "" } { return }
            }
            set p $id
            while { $p != 0 && ![$tree item state get $p condition] && \
                ![$tree item state get $p blockdata] && \
                ![$tree item state get $p container]} {
                set p [$tree item parent $p]
            }
            
            set found 1
            if { [$tree item state get $id condition] } {                
                set domNodeCnd [$tree item element cget $id 0 e_text -data]
                set a0 [dict create add_group 1 edit_group 1 delete_group 1]
                set actions [dict merge $a0 [split [$domNodeCnd @actions ""] ,]]
                if { [dict get $actions add_group] } {
                    set open_window 1                    
                    if { [$domNodeCnd @open_window ""] == "0" } {
                        set open_window 0
                    }
                    if { [$domNodeCnd selectNodes {edit_command[@edit_type='exclusive' or @edit_type='exclusive_blockdata']}] ne ""} {
                        $self _edit_command_eval -item $id $domNodeCnd
                    }
                    $self apply_cond -what apply -open_window $open_window
                } else { set found 0 }
            } elseif { [$tree item state get $id blockdata] } {
                set domNode [$tree item element cget $id 0 e_text -data]
                set xp1 {value|container/value|container/container/value}
                
                if { [$domNode selectNodes {edit_command[@edit_type='exclusive' or @edit_type='exclusive_blockdata']}] ne ""} {
                    $self _edit_command_eval -item $id $domNode
                } elseif { [$domNode selectNodes $xp1] ne "" || ![$self edit_value_or_blockdata] } {
                    set del_cancel_button 0
                    if { [$domNode @del_cancel_button ""] == "1" } {
                        set del_cancel_button 1
                    }
                    set open_window 1                    
                    if { [$domNode @open_window ""] == "0" } {
                        set open_window 0
                    }                   
                    $self apply_blockdata -del_cancel_button $del_cancel_button \
                        -open_window $open_window                    
                }
            } elseif { [$tree item state get $id container] } {
                set domNode [$tree item element cget $id 0 e_text -data]
                set xp1 {value|container/value}
                if { [$domNode selectNodes {edit_command[@edit_type='exclusive']}] ne ""} {
                    if {[$domNode selectNodes {edit_command[@procmod]}] ne ""} {
                        $self _edit_command_eval -item $id -procmod 1 $domNode
                    } else {
                        $self _edit_command_eval -item $id $domNode
                    }
                } elseif { [$domNode selectNodes $xp1] ne "" } {
                    set del_cancel_button 0
                    if { [$domNode @del_cancel_button ""] == "1" } {
                        set del_cancel_button 1
                    }
                    set open_window 1                    
                    if { [$domNode @open_window ""] == "0" } {
                        set open_window 0
                    }                    
                    $self apply_blockdata -del_cancel_button $del_cancel_button \
                        -open_window $open_window                    
                } elseif { [$tree item state get $p condition] } {
                    set domNodeCnd [$tree item element cget $p 0 e_text -data]
                    set a0 [dict create add_group 1 edit_group 1 delete_group 1]
                    set actions [dict merge $a0 [split [$domNodeCnd @actions ""] ,]]
                    if { [dict get $actions add_group] } {
                        $self apply_cond -what edit
                    } else { set found 0 }
                } else { set found 0 }
            } elseif { [$tree item state get $id value] } {
                $self edit_value_or_blockdata
            } elseif { [$tree item state get $id group] } {
                set domNodeCnd [$tree item element cget $p 0 e_text -data]
                set a0 [dict create add_group 1 edit_group 1 delete_group 1]
                set actions [dict merge $a0 [split [$domNodeCnd @actions ""] ,]]
                if { [dict get $actions edit_group] } {
                    $self apply_cond -what edit
                } else { set found 0 }
            } else { set found 0 }
            if { $found } {
                if { $execute_select_handler ne "" } {
                    after cancel $execute_select_handler
                    set execute_select_handler ""
                }
            }
            return
        } else {
            set execute_select_pressed $id
            if { $execute_select_handler ne "" } {
                after cancel $execute_select_handler
                set execute_select_handler ""
            }
            #after 300 [list unset -nocomplain [varname execute_select_pressed]]
            set execute_select_handler [after 900 [list catch [mymethod \
                            execute_select_handler $id]]]
        }
        
        #         set isopen [$tree item state get $id open]
        #         if { $isopen } {
            #             $tree collapse -recurse $id
            #         } else {
            #             set sel [$tree selection get]
            #             if 0 {
                #                 $tree collapse all
                #                 foreach i [$tree item ancestors $id] {
                    #                     $tree expand $i
                    #                     if { [lsearch -exact $sel $i] != -1 } {
                        #                         $tree selection add $i
                        #                     }
                    #                 }
                #             }
            #             $tree expand $id
            #             if { [lsearch -exact $sel $id] != -1 } {
                #                 $tree selection add $id
                #             }
            #         }
    }
    method execute_select_handler { id } {
        set execute_select_pressed ""
        set execute_select_handler ""
        set isopen [$tree item state get $id open]
        
        if { $isopen } {
            #$tree collapse -recurse $id
        } else {
            #$tree expand $id
            set sel [$tree selection get]
            if { [lsearch -exact $sel $id] != -1 } {
                $tree selection add $id
            }
        }
    }
    method contextual_menu_add_cmds { m cmdsvar } {
        upvar $cmdsvar cmds
        
        set menu_types [list cascade checkbutton command radiobutton separator]
        
        foreach "img txt cmd" $cmds {
            if { $txt eq "-" } {
                $m add separator
            } elseif { [lindex $cmd 0] eq "-" } {
                set i 1
                while { [winfo exists $m.m$i] } { incr i }
                if { $img ne "" && [info command $img] eq "" } {
                    set img $gid_groups_conds::images($img-16)
                }
                $m add cascade -image $img -label $txt -compound left \
                    -menu $m.m$i
                menu $m.m$i -tearoff 0
                upvar [lindex $cmd 1] [lindex $cmd 1] 
                $self contextual_menu_add_cmds $m.m$i [lindex $cmd 1]
            } elseif { [lsearch -exact $menu_types $img] != -1 } {
                eval [list $m add $img -label $txt] $cmd
            } else {
                if { $img ne "" && [info command $img] eq "" } {
                    set img $gid_groups_conds::images($img-16)
                }
                $m add command -image $img -label $txt -compound left \
                    -command $cmd
            }
        }
    }
    method collapse_expand { collapse_expand } {
        switch $collapse_expand {
            expand {
                foreach id [$tree selection get] {
                    $tree item expand $id -recurse
                }
            }
            expandall {
                $tree item expand 0 -recurse
            }
            collapse {
                foreach id [$tree selection get] {
                    $tree item collapse $id -recurse
                }
            }
            collapseall {
                $tree item collapse 0 -recurse
            }
            collapseother {
                set id [lindex [$tree selection get] 0]
                set domNode [$tree item element cget $id 0 e_text -data]
                set item 0
                foreach node [lrange [$domNode selectNodes ancestor-or-self::*] 1 end] {
                    set state [get_domnode_attribute $node state]
                    if { $state eq "hidden_this" } { continue }
                    foreach i [$tree item children $item] {
                        set d [$tree item element cget $i 0 e_text -data]
                        if { $d ne $node } {
                            $tree item collapse $i
                        } else {
                            set item $i
                        }
                    }
                }
            }
            viewthis {
                foreach id [$tree selection get] {
                    $tree item expand $id
                }
                $self collapse_expand collapseother
            }
        }
    }
    method contextual_menu { x y X Y } {
        
        $self popup_help_deactivate
        
        if { $x ne "top" } {
            set ident [$tree identify $x $y]
            set id [lindex $ident 1]
            if { [lindex $ident 0] ne "item" || $id eq "" } {
                $self popup_help_reactivate
                return
            }
            set domNode [$tree item element cget $id 0 e_text -data]
            
            set p $id
            while { $p != 0 && ![$tree item state get $p condition] && \
                ![$tree item state get $p blockdata] && \
                ![$tree item state get $p container] } {
                set p [$tree item parent $p]
            }
        } else {
            set id ""
        }
        set cmds ""
        set draw_symbols 0
        
        if { $id eq "" } {
            if { !$gid_groups_conds::is_draw_post } {
                lappend cmds "" [= "Toggle window position"] \
                    [list gid_groups_conds::open_conditions toggle]
            }
            lappend cmds "" [= "Toggle search box"] \
                [mymethod toggle_search]
        } elseif { [$tree item state get $id condition] } {
            set open_window 1                    
            if { [$domNode @open_window ""] == "0" } {
                set open_window 0
            }
            if { $open_window } {
                lappend cmds [gid_groups_conds::giveimage bookmark_folder-16] \
                    [= "Apply to entities"] [mymethod apply_cond]
            }
            if { [llength [$domNode selectNodes group|groupList]] } {
                lappend cmds [gid_groups_conds::giveimage remove-16] \
                    [= "Unassign groups"] [mymethod delete_cond_group] \
                    - - - 
                lappend cmds [gid_groups_conds::giveimage colorize-16] \
                    [= "Draw symbols"] [mymethod draw_groups_do symbols]
                lappend cmds [gid_groups_conds::giveimage colorize-16] \
                    [= "Draw"]... [list - cmds_draw]
                if { [$domNode selectNodes symbol] ne "" } { set draw_symbols 1 }
                lappend cmds [gid_groups_conds::giveimage colorize-16] \
                    [= "List entities"] [mymethod list_entities]
            }
        } elseif { [$tree item state get $id blockdata] } {
            set msg [= "Edit"]                                  
            set xp {ancestor-or-self::*/edit_command[@edit_type='exclusive' or @edit_type='exclusive_blockdata']}
            if { [$domNode @createedittext ""] != ""} {
                set msg [= "%s" [$domNode @createedittext]]
            }   
            if {[$domNode @options 1]} {         
                if {![$domNode @edit_command_eval 1]} {
                    lappend cmds [gid_groups_conds::giveimage bookmark_folder-16] \
                        $msg [mymethod apply_blockdata]
                } elseif { [$domNode selectNodes $xp] ne "" } {
                    lappend cmds [gid_groups_conds::giveimage bookmark_folder-16] \
                        $msg [mymethod _edit_command_eval -item $id $domNode]
                } else {
                    lappend cmds [gid_groups_conds::giveimage bookmark_folder-16] \
                        $msg [mymethod apply_blockdata]
                }  
            }          
            if { [$domNode @editable_name ""] eq "1" || \
                [$domNode @editable_name ""] eq "unique" } {
                lappend cmds [gid_groups_conds::giveimage bookmark_folder-16] \
                    [= "Rename"] [mymethod edit_value_or_blockdata]
            }           
            
            if { [$domNode @sequence 0] } {
                set xp {condition/group|condition/groupList|}
                append xp {container/condition/group|} \
                    {container/condition/groupList}
                set num_groups [llength [$domNode selectNodes $xp]]
                
                set xp {/*/blockdata[@n="Internal data"]/value[@n="copy_sequence_style"]}
                set n [$domNode selectNodes $xp]
                if { $n eq "" } {
                    set copy_sequence_style copy_both
                } else {
                    set copy_sequence_style [$n @v copy_both]
                }
                if { $copy_sequence_style eq "copy_both" && !$num_groups } {
                    set copy_sequence_style copy_cond
                }
                switch $copy_sequence_style {
                    copy_both {
                        lappend cmds [gid_groups_conds::giveimage editcopy-16] \
                            [= "Copy (no assigned groups)"] [mymethod copy_block_data] \
                            [gid_groups_conds::giveimage editcopy-16] \
                            [= "Copy (with assigned groups)"] [mymethod copy_block_data \
                                -copy_cond_groups]
                    }
                    copy_cond {
                        lappend cmds [gid_groups_conds::giveimage editcopy-16] \
                            [= "Copy"] [mymethod copy_block_data]
                    }
                    copy_cond_groups {
                        lappend cmds [gid_groups_conds::giveimage editcopy-16] \
                            [= "Copy"] [mymethod copy_block_data -copy_cond_groups]
                    }
                }
                if { [$domNode @n] eq "material" } {
                    set xp {ancestor::container[@n="materials"]/blockdata[@n="material"]|}
                    append xp {ancestor::container[@n="materials"]/container/blockdata[@n="material"]}
                    set cmds_move_containers ""
                    set xp_c {ancestor::container[@n="materials"]/container}
                    foreach cNode [$domNode selectNodes $xp_c] {
                        if { [[$domNode parentNode] @pn ""] eq [$cNode @pn] } { continue }
                        lappend cmds_move_containers [$cNode @pn]
                    }
                } else {
                    set xp [format_xpath {../blockdata[@n=%s]} [$domNode @n]]
                }
                
                if { [$domNode @can_delete_last_item 0] == 1 } {
                    set options(-can_delete_last_material) 1    
                }
                
                set num [llength [$domNode selectNodes $xp]]                          
                
                if { $num > 1 } {
                    lappend cmds [gid_groups_conds::giveimage up-16] \
                        [= "Change position"] [list - cmds_move]
                    lappend cmds - - - \
                        [gid_groups_conds::giveimage remove-16] \
                        [= "Delete"] [mymethod delete_disable_blockdata]
                } elseif { $options(-can_delete_last_material) && $num == 1 } {
                    lappend cmds - - - \
                        [gid_groups_conds::giveimage remove-16] \
                        [= "Delete"] [mymethod delete_disable_blockdata]
                } elseif { [$domNode @state normal] ne "disabled" && 
                    [$domNode @sequence_type any] eq "non_void_disabled" } {
                    lappend cmds - - - \
                        [gid_groups_conds::giveimage remove-16] \
                        [= "Disable"] [mymethod delete_disable_blockdata]
                } elseif { [$domNode @active 1] && 
                    [$domNode @sequence_type any] eq "non_void_deactivated" } {
                    lappend cmds - - - \
                        [gid_groups_conds::giveimage remove-16] \
                        [= "Deactivate"] [mymethod delete_disable_blockdata]
                } elseif { [$domNode @state normal] eq "disabled" && 
                    [$domNode @sequence_type any] eq "non_void_disabled" } {
                    lappend cmds - - - \
                        [gid_groups_conds::giveimage remove-16] \
                        [= "Enable"] [mymethod enable_activate_blockdata]
                } elseif { [$domNode @active 1] == 0 && 
                    [$domNode @sequence_type any] eq "non_void_deactivated" } {
                    lappend cmds - - - \
                        [gid_groups_conds::giveimage remove-16] \
                        [= "Activate"] [mymethod enable_activate_blockdata]
                }
                if { [$domNode selectNodes .//condition] ne "" } {
                    lappend cmds [gid_groups_conds::giveimage colorize-16] \
                        [= "Draw symbols"] [mymethod draw_groups_do symbols]
                    lappend cmds - - - [gid_groups_conds::giveimage colorize-16] \
                        [= "Draw"]... [list - cmds_draw]
                    if { [$domNode selectNodes .//condition/symbol] ne "" } {
                        set draw_symbols 1
                    }
                    lappend cmds [gid_groups_conds::giveimage colorize-16] \
                        [= "List entities"] [mymethod list_entities]
                }
                
                if { $options(-have_import_export) && [$domNode @n] eq "material" } {
                    lappend cmds - - - \
                        "" [= "Import/export materials"] [mymethod import_export_materials]
                }
            }
        } elseif { [$tree item state get $id container] } {
            set xp {blockdata[@sequence='1']}
            set bd [lindex [$domNode selectNodes $xp] 0]
            if { $bd ne "" } {
                set name [$bd @pn [$bd @n]]
                set createedittext [= "Create new %s" $name]
                if { [$domNode @createedittext ""] != ""} {
                    set createedittext [$domNode @createedittext]
                }
                if { [$domNode @create_new 1] } {                   
                    lappend cmds [gid_groups_conds::giveimage editcopy-16] \
                        $createedittext [mymethod copy_block_data]
                }
            }                        
            set xp {ancestor-or-self::container[@n="materials"]}
            if { $options(-have_import_export) && [$domNode selectNodes $xp] ne "" } {
                if { $cmds ne "" } { lappend cmds - - - }
                lappend cmds "" [= "Import/export materials"] [mymethod import_export_materials]
            } 
            set msg [= "Edit"]           
            if { [$domNode @createedittext ""] != ""} {
                set msg [= "%s" [$domNode @createedittext]]
            }                 
            if { [$domNode selectNodes {edit_command[@edit_type='exclusive']}] ne "" && \
                ![$domNode @create_new 1]} {                
                lappend cmds [gid_groups_conds::giveimage bookmark_folder-16] \
                    $msg [mymethod _edit_command_eval -item $id $domNode]
            } elseif { [$tree item state get $p blockdata] && \
                [$domNode selectNodes value|container/value] ne "" && \
                ![$domNode @create_new 1]} {
                lappend cmds [gid_groups_conds::giveimage bookmark_folder-16] \
                    $msg [mymethod apply_blockdata]
            } elseif { [$tree item state get $p condition] } {
                if { $cmds ne "" } {
                    lappend cmds - - -
                }
                lappend cmds [gid_groups_conds::giveimage bookmark_folder-16] \
                    [= "Copy and apply to entities"] [mymethod apply_cond -what apply] \
                    [gid_groups_conds::giveimage bookmark_folder-16] \
                    [= "Edit"] [mymethod apply_cond -what edit]
                lappend cmds [gid_groups_conds::giveimage colorize-16] \
                    [= "Draw symbols"] [mymethod draw_groups_do symbols]
                lappend cmds - - - [gid_groups_conds::giveimage colorize-16] \
                    [= "Draw"]... [list - cmds_draw]
                set domNodeP [$tree item element cget $p 0 e_text -data]
                if { [$domNodeP selectNodes symbol] ne "" } {
                    set draw_symbols 1
                }
                lappend cmds [gid_groups_conds::giveimage colorize-16] \
                    [= "List entities"] [mymethod list_entities]
            } elseif { [$domNode selectNodes .//condition] ne "" } {
                lappend cmds [gid_groups_conds::giveimage colorize-16] \
                    [= "Draw symbols"] [mymethod draw_groups_do symbols]
                lappend cmds [gid_groups_conds::giveimage colorize-16] \
                    [= "Draw"]... [list - cmds_draw]
                if { [$domNode selectNodes .//condition/symbol] ne "" } {
                    set draw_symbols 1
                }
                lappend cmds [gid_groups_conds::giveimage colorize-16] \
                    [= "List entities"] [mymethod list_entities]
            }
        } elseif { [$tree item state get $id value] } {
            if { [$tree item state get $p condition] } {
                lappend cmds [gid_groups_conds::giveimage bookmark_folder-16] \
                    [= "Apply more"] [mymethod apply_cond -what select_more]
            }
            if { [$domNode selectNodes function] eq "" } {
                lappend cmds [gid_groups_conds::giveimage bookmark_folder-16] \
                    [= "Edit"] [mymethod edit_value_or_blockdata]
            }
            if { [$tree item state get $p condition] } {
                lappend cmds [gid_groups_conds::giveimage remove-16] \
                    [= "Unassign group"] [mymethod delete_cond_group]
                
                lappend cmds [gid_groups_conds::giveimage colorize-16] \
                    [= "Draw symbols"] [mymethod draw_groups_do symbols]
                lappend cmds - - - [gid_groups_conds::giveimage colorize-16] \
                    [= "Draw"]... [list - cmds_draw]
                set domNodeP [$tree item element cget $p 0 e_text -data]
                if { [$domNodeP selectNodes symbol] ne "" } {
                    set draw_symbols 1
                }
                lappend cmds [gid_groups_conds::giveimage colorize-16] \
                    [= "List entities"] [mymethod list_entities]
            }
        } elseif { [$tree item state get $id group] || [$tree item state get $id groupList] } {
            
            #             [gid_groups_conds::giveimage bookmark_folder-16] \
                #                 [= "Rename group"] [mymethod rename_cond_group]
            
            lappend cmds [gid_groups_conds::giveimage bookmark_folder-16] \
                [= "Copy and apply to entities"] [mymethod apply_cond -what apply] \
                [gid_groups_conds::giveimage bookmark_folder-16] \
                [= "Apply more"] [mymethod apply_cond -what select_more] \
                [gid_groups_conds::giveimage bookmark_folder-16] \
                [= "Edit"] [mymethod apply_cond -what edit] \
                [gid_groups_conds::giveimage remove-16] \
                [= "Unassign group"] [mymethod delete_cond_group] \
                [gid_groups_conds::giveimage bookmark_folder-16] \
                [= "View group in groups window"] [mymethod view_cond_group_in_group_window] \
                
            lappend cmds - - - [gid_groups_conds::giveimage colorize-16] \
                [= "Draw symbols"] [mymethod draw_groups_do symbols]
            lappend cmds [gid_groups_conds::giveimage colorize-16] \
                [= "Draw"]... [list - cmds_draw]
            lappend cmds [gid_groups_conds::giveimage up-16] \
                [= "Change position"] [list - cmds_move]
            
            set domNodeP [$tree item element cget $p 0 e_text -data]
            if { [$domNodeP selectNodes symbol] ne "" } {
                set draw_symbols 1
            }
            lappend cmds [gid_groups_conds::giveimage colorize-16] \
                [= "List entities"] [mymethod list_entities]
        }
        
        if { $cmds ne "" } { lappend cmds - - - }
        
        if { $id ne "" } {
            lappend cmds mac-expand [= "View this"] [mymethod collapse_expand viewthis]
            lappend cmds mac-expand [= "Expand"] [list - cmds_expand]
            lappend cmds - - -
        }
        lappend cmds "" [= "Tree preferences"] [mymethod tree_preferences]
        
        set cmds_draw [list [gid_groups_conds::giveimage colorize-16] \
                [= "Draw groups"] [mymethod draw_groups_do groupnames] \
                [gid_groups_conds::giveimage colorize-16] \
                [= "Draw values"] [mymethod draw_groups_do values] \
                ]
        if { $draw_symbols } {
            lappend cmds_draw [gid_groups_conds::giveimage colorize-16] \
                [= "Draw symbols"] [mymethod draw_groups_do symbols]
        }
        lappend cmds_draw - - - [gid_groups_conds::giveimage colorize-16] \
            [= "Draw all symbols"] [mymethod draw_groups_do symbols_all] \
            - - -
        
        lappend cmds_draw checkbutton [= "Auto change render"] \
            [list -image [gid_groups_conds::giveimage colorize-16] -variable \
                gid_groups_conds::do_change_opengl_rendermode -compound left]
        
        lappend cmds_draw checkbutton [= "Auto hide surfaces"] \
            [list -image [gid_groups_conds::giveimage colorize-16] -variable \
                gid_groups_conds::draw_auto_hide_surfaces -compound left]
        
        set cmds_expand [list mac-expand \
                [= "Expand"] [mymethod collapse_expand expand] \
                mac-collapse2 \
                [= "Collapse other"] [mymethod collapse_expand collapseother] \
                mac-collapse \
                [= "Collapse"] [mymethod collapse_expand collapse] \
                mac-expand2 \
                [= "Expand all"] [mymethod collapse_expand expandall] \
                mac-collapse2 \
                [= "Collapse all"] [mymethod collapse_expand collapseall] \
                ]
        
        set cmds_move [list [gid_groups_conds::giveimage up-16] \
                [= "Up"] [mymethod reorder_blockdata_or_group up] \
                [gid_groups_conds::giveimage down-16] \
                [= "Down"] [mymethod reorder_blockdata_or_group down] \
                - - - \
                [gid_groups_conds::giveimage top-16] \
                [= "First"] [mymethod reorder_blockdata_or_group first] \
                [gid_groups_conds::giveimage bottom-16] \
                [= "Last"] [mymethod reorder_blockdata_or_group last] \
                ]
        
        if { [info exists cmds_move_containers] } {
            lappend cmds_move - - -
            foreach c $cmds_move_containers {
                lappend cmds_move "" [= "Move to '%s'" $c] \
                    [mymethod reorder_blockdata_or_group moveto $c]
            }
        }
        
        destroy $win.menu
        set m [menu $win.menu -tearoff 0]
        
        $self contextual_menu_add_cmds $m cmds
        tk_popup $m $X $Y
        # WARNING: check if it works for Linux
        $self popup_help_reactivate
    }
    
    method add { parent name image state data { before "" } } {
        set treenode [$tree item create]
        if { $before eq "" } {
            $tree item lastchild $parent $treenode
        } else {
            $tree item prevsibling $before $treenode
        }
        
        $tree item state set $treenode [list $state]
        $tree item style set $treenode 0 s_text
        $tree item complex $treenode [list \
                [list e_image -image $image] [list e_text -text $name -data $data] \
                [list f_image -image ""]]       
        
        if { $parent != 0 && ![$tree item cget $parent -button] } {
            set name [$tree item text $parent 0]
            $tree item configure $parent -button yes
            #            $tree item style set $parent 0 s_text 1 s_image
            #             $tree item complex $parent \
                #                 [list [list e_text -text $name]] \
                #                 [list [list e_image -image $img]]
        }
        #set ::TreeCtrl::Priv(DirCnt,$tree) $treenode
        return $treenode
    }
    method contextual_entry { entry x y } {
        destroy $win.menu
        set m [menu $win.menu -tearoff 0]
        
        $m add command -label [= "Select all"] -command \
            [list $entry selection range 0 end]
        $m add command -label [= "Cut"] -image $gid_groups_conds::images(editcut-16) \
            -compound left -command \
            [list event generate $entry <<Cut>>]
        $m add command -label [= "Copy"] -image $gid_groups_conds::images(editcopy-16) \
            -compound left -command \
            [list event generate $entry <<Copy>>]
        $m add command -label [= "Paste"] -image $gid_groups_conds::images(editpaste-16) \
            -compound left -command \
            [list event generate $entry <<Paste>>]
        $m add separator
        $m add checkbutton -label [= "Transform unit value"] -image \
            $gid_groups_conds::images(attach-16) \
            -compound left -variable [varname transform_unit_value]
        
        tk_popup $m $x $y
    }
    #     method _setup_units_combo { n valueNode combo } {
        # 
        #         set unit_magnitude [get_domnode_attribute $valueNode unit_magnitude]
        #         if { [llength [split $unit_magnitude ,]] > 1 } {
            #             set ov $entryvalues(_selection_entity_)
            #         } else {
            #             set ov ""
            #         }
        #         foreach "value unit units unitsN" [gid_groups_conds::give_value_in_active_unit \
            #             -ov $ov $valueNode] break
        #         $combo configure -values $unitsN
        # 
        #         set unitN [gid_groups_conds::units_to_nice_units $unit]
        #         $combo state !readonly
        #         set entryvalues_units_old($n) $unitN
        #         set entryvalues($n) $value
        #         set entryvalues_units($n) $unitN
        #         $combo state readonly
        #     }
    #     method _manage_units_combo { domNode magnitude n } {
        #         if { ![info exists transform_unit_value] || !$transform_unit_value } { return }
        #         if { $entryvalues_units($n) eq $entryvalues_units_old($n) } { return }
        # 
        #         set unit_from [gid_groups_conds::nice_units_to_units $entryvalues_units_old($n)]
        #         set unit_to [gid_groups_conds::nice_units_to_units $entryvalues_units($n)]
        #         catch {
            #             set nv [gid_groups_conds::convert_unit_value $magnitude \
                #                     $entryvalues($n) $unit_from $unit_to]
            #             set entryvalues($n) $nv
            #         }
        #         set entryvalues_units_old($n) $entryvalues_units($n)
        #     }
    method _fill_data_window_function_update { funcs n button env } {
        
        if { ![winfo exists $button] } { return }
        
        if { ![info exists entryvalues_state($n)] } {
            set stateG normal
        } else {
            set stateG $entryvalues_state($n)
        }
        if { $entryvalues_function($n) ne "" } {
            set nameList ""
            set valueList ""
            foreach fvarNode [$entryvalues_function($n) selectNodes functionVariable] {
                lappend nameList [$fvarNode @pn]
                set txt [$fvarNode @pn]
                if { [$fvarNode @units ""] ne "" } {
                    append txt " ([$fvarNode @units])"
                }
                set vline [list $txt]
                foreach vNode [$fvarNode selectNodes value] {
                    set txt "[$vNode @pn]: [$vNode @v]"
                    if { [$vNode selectNodes function] ne "" } {
                        append txt "..."
                    }
                    lappend vline $txt
                }
                if { [llength $vline] > 10 } {
                    set vline [lrange $vline 0 9]
                    lappend vline "..."
                }
                lappend valueList [join $vline \n]
            }
            foreach i $env {
                lassign $i type widget
                if { $type ni "ComboBox entry entry_units checkbutton" } { continue }
                if { ![winfo exists $widget] } { continue }
                if { $type eq "entry_units" } {
                    set idx 0
                    foreach j [winfo children $widget] {
                        if { $idx == 0 } {
                            $j configure -textvariable \
                                [myvar entryvalues_function_var($n)]
                            set entryvalues_function_var($n) "f([join $nameList ,])..."
                            $j state disabled
                        }
                        set h [lrange [split [gid_groups_conds::register_popup_help $j] \n] 0 0]
                        eval lappend h $valueList
                        gid_groups_conds::register_popup_help $j [join $h \n]
                        incr idx
                    }
                    set h [lrange [split [gid_groups_conds::register_popup_help $widget] \n] 0 0]
                    lappend h {*}$valueList
                    gid_groups_conds::register_popup_help $widget [join $h \n]
                } else {
                    $widget state !disabled
                    if { $type ne "checkbutton" } {
                        $widget configure -textvariable \
                            [myvar entryvalues_function_var($n)]
                        set entryvalues_function_var($n) "f([join $nameList ,])..."
                    } else {
                        $widget configure -variable \
                            [myvar entryvalues_function_var($n)]
                        set entryvalues_function_var($n) 1
                    }
                    $widget state disabled
                    set h [lrange [split [gid_groups_conds::register_popup_help $widget] \n] 0 0]
                    eval lappend h $valueList
                    gid_groups_conds::register_popup_help $widget [join $h \n]
                }
            }
        } else {
            set nfuncs [llength $funcs]
            if {$funcs == ""} {
                set state !disabled
            } elseif { [lsearch $funcs scalar] == -1 } {
                set state disabled
            } else {
                if { $stateG ne "disabled" } {
                    set state !disabled
                } else {
                    set state disabled
                }
                incr nfuncs -1
            }
            set h [gid_groups_conds::register_popup_help $button]
            if { $nfuncs == 0 } {
                $button state disabled
                if { [llength [split $h \n]] == 1 } {
                    append h "\nThere are no functions to apply in this configuration"
                }
            } else {
                if { $stateG ne "disabled" } {
                    $button state !disabled
                }
                if { [llength [split $h \n]] == 2 } {
                    set h [lindex [split $h \n] 0]
                }
            }
            gid_groups_conds::register_popup_help $button $h
            foreach i $env {
                lassign $i type widget
                if { $type ni "ComboBox entry entry_units checkbutton" } { continue }
                if { $type eq "entry_units" } {
                    set idx 0
                    foreach j [winfo children $widget] {
                        if { $idx == 0 } {
                            $j state $state
                            $j configure -textvariable [myvar entryvalues($n)]
                        } else {
                            $j configure -textvariable [myvar entryvalues_units($n)]
                        }
                        set h [lindex [split [gid_groups_conds::register_popup_help $j] \n] 0]
                        gid_groups_conds::register_popup_help $j $h
                        incr idx
                    }
                    set h [lindex [split [gid_groups_conds::register_popup_help $widget] \n] 0]
                    gid_groups_conds::register_popup_help $widget $h
                } else {
                    $widget state $state
                    if { $type ne "checkbutton" } {
                        $widget configure -textvariable [myvar entryvalues($n)]
                    } else {
                        $widget configure -variable [myvar entryvalues($n)]
                    }
                    set h [lindex [split [gid_groups_conds::register_popup_help $widget] \n] 0]
                    gid_groups_conds::register_popup_help $widget $h
                }
            }
        }
        # this is to activate the trace
        if { [info exists entryvalues($n)] } {
            set entryvalues($n) $entryvalues($n)
        }
    }
    method has_changes {} {
        return $entryvalues_haschanges
    }
    method _fill_data_window_reset_changes { w } {
        set entryvalues_haschanges 0
        set cmd [list set [myvar entryvalues_haschanges] 1]
        foreach i [list entryvalues entryvalues_units entryvalues_function \
                entryvalues_function_var] {
            trace add variable [myvar $i] write "$cmd ;#"
            bind $w <Destroy> +[list trace remove variable [myvar $i] write "$cmd ;#"]
        }
        bind $w <Destroy> +[list set [myvar entryvalues_haschanges] 0]
        set entryvalues_mat_changes_stack ""
    }
    method _fill_data_window { w domNode title { level 0 } } {
        
        if { $level == 0 } { array unset entry_node_values }
        set containerNodes [$domNode selectNodes container]
        set valueNodes [$domNode selectNodes value|edit_command]
        if { ![llength $valueNodes] && [$domNode @isvalue 0] } {
            set valueNodes [list $domNode]
        }
        set ret ""
        if { [llength $valueNodes] } {
            set ret0 [$self _fill_data_window_values $w $domNode $valueNodes $title]
            lappend ret {*}$ret0
        }
        
        if { $level == 0 && [llength $containerNodes] } {
            if { [$domNode hasAttribute n] } { 
                set n [$domNode @n]
            } else { 
                set n ""
            }
            if { $n eq "material" || [llength $containerNodes] > 1 } {
                set ret0 [$self _fill_data_window_containersL0 $w $domNode $containerNodes]
                lappend ret {*}$ret0
                set containerNodes ""
            }
        }
        if { [llength $containerNodes] } {
            set ret0 [$self _fill_data_window_containers $w $containerNodes $level]
            lappend ret {*}$ret0
        }
        #         if { $level == 1 || $level == 0}  // Revisar para que salgan los Botones
        if { $level == 0} {
            set xp {ancestor-or-self::condition[@local_axes]}
            if { [$domNode selectNodes $xp] ne "" } {
                ttk::frame $w.fla
                ttk::label $w.fla.l -text [= "Local Axes"]: \
                    -image local_axes-16 -compound left
                cu::menubutton_button $w.fla.mbb -text [= "Open"] \
                    -menu $w.fla.mbb.m -command gid_groups_conds::local_axes_window
                menu $w.fla.mbb.m -tearoff 0
                
                gid_groups_conds::local_axes_menu -ov_variable \
                    [myvar entryvalues(_selection_entity_)] \
                    $w.fla.mbb.m
                
                set help [= "Apply or draw local axes"]
                gid_groups_conds::register_popup_help $w.fla.l \
                    $help
                gid_groups_conds::register_popup_help $w.fla.mbb \
                    $help
                
                grid $w.fla.l $w.fla.mbb -sticky w
                grid columnconfigure $w.fla 1 -weight 1
                grid $w.fla -sticky nw -padx 2 -pady "2 1"
            }
            #             set xp {edit_command|ancestor-or-self::condition/edit_command}            
            #             set edit_command [$domNode selectNodes $xp]
            #             if { $edit_command ne "" } {
                #                 set pn [get_domnode_attribute $edit_command pn]
                #                 set edit_type [get_domnode_attribute $edit_command edit_type]
                #                 ttk::button $w.fbcommand -text $pn \
                    #                     -command [mymethod _edit_command_eval_window $w $domNode $edit_type]              
                #                 grid $w.fbcommand -sticky nw -padx 2 -pady "2 1"
                #                 set help [get_domnode_attribute $edit_command help ""]
                #                 if {$help != ""} {                
                    #                     gid_groups_conds::register_popup_help $w.fbcommand \
                        #                         $help             
                    #                 }   
                #             }
            set f [lindex $ret 0]
            if { [winfo exists $f.e0.e] } {
                focus $f.e0.e
                catch { $f.e0.e selection range 0 end }
            } elseif { [winfo exists $f.e0] } {
                focus $f.e0
                catch { $f.e0 selection range 0 end }
            }
        }
        return $ret
    }
    method _fill_data_window_containersL0 { w domNode containerNodes } {
        cu::notebook $w.fc
        
        set j 0
        foreach containerNode $containerNodes {
            set xp {value|edit_command|container/value|container/edit_command}
            if { [llength [$containerNode selectNodes $xp]] == 0 } { continue }            
            set isinfo 0
            foreach icontNode [$containerNode selectNodes $xp] {
                if { [get_domnode_attribute $icontNode state] eq "normal" || \
                    [get_domnode_attribute $icontNode state] eq "disabled" || \
                    [get_domnode_attribute $icontNode state] eq "readonly" || \
                    [get_domnode_attribute $icontNode state] eq ""} { 
                    set isinfo 1; break 
                }                
            }
            if { !$isinfo } { continue }            
            set state [get_domnode_attribute $containerNode state normal]
            set f [ttk::frame $w.fc.f$j]
            $w.fc add $f -text [string tolower [get_domnode_attribute \
                        $containerNode pn]] -sticky nsew
            #set f [$w.fc insert end p$j -text [$containerNode @pn]]
            $self _fill_data_window $f $containerNode \
                [get_domnode_attribute $containerNode pn] 1
            set entry_node_values($containerNode) \
                [list [list TNotebook $w.fc $j]]
            if { $state eq "hidden" } {
                $w.fc tab $j -state $state
            }
            incr j
        }
        #$w.fc compute_size
        grid $w.fc -sticky nwes
        grid columnconfigure $w 0 -weight 1
        grid rowconfigure $w 0 -weight 1
        #$w.fc raise p0
        if { $j } {
            foreach w_i [$w.fc tabs] {
                if { [$w.fc tab $w_i -state] eq "normal" } {
                    $w.fc select $w_i
                    break
                }
            }
        }
        return $w.fc
    }
    method _fill_data_window_containers { w containerNodes level } {
        set idx 0
        foreach containerNode $containerNodes {
            #grid rowconfigure $w $idx -weight 0
            set xp {value|edit_command|container/value|container/edit_command}
            if { [llength [$containerNode selectNodes $xp]] == 0 } { continue }
            
            $self _fill_data_window $w $containerNode \
                [get_domnode_attribute $containerNode pn] \
                [expr {$level+1}]
            incr idx
        }
        #grid rowconfigure $w $idx -weight 1
    }
    method _manage_vector_entry { what args } {              
        switch $what {
            to_components {
                lassign $args n dimensions
                set list [split $entryvalues($n) ","]
                for { set i 1 } { $i <= $dimensions } { incr i } {
                    if { [set! entryvalues($n,$i)] != [lindex $list $i-1] } {
                        set entryvalues($n,$i) [lindex $list $i-1]
                    }
                }
            }
            from_components {
                lassign $args n dimensions
                set list ""
                for { set i 1 } { $i <= $dimensions } { incr i } {
                    lappend list $entryvalues($n,$i)
                }
                set v [join $list ","]
                if { $v ne $entryvalues($n) } {
                    set entryvalues($n) $v
                }
            }
            select_vector {
                lassign $args n frame 1st_coord format 2nd_coord               
                set grab ""
                if {[winfo exists $frame]} {
                    if { [grab current $frame] ne "" } {
                        set grab [grab current $frame]
                        grab release $grab
                    }                        
                }
                if { $1st_coord eq "" } {
                    set msg [= "Select origin point"]            
                } elseif {$2nd_coord eq ""} {
                    set msg [= "Select destination point"]
                }                
                set withdraw ""   
                #                 set grab ""                            
                set f [ttk::frame $frame.t]
                place $f -x 0 -y 0 -relwidth 1 -relheight 1 -anchor nw                
                ttk::label $frame.t.l -text $msg
                ttk::button $frame.t.b -image [gid_groups_conds::giveimage close17] -command \
                    [mymethod _manage_vector_entry select_vector_done $n $frame.t $grab $withdraw "" "" $format]
                grid $frame.t.l $frame.t.b -padx 2 -pady 2
                
                if {$1st_coord == ""} {
                    $self manage_buttons_set PICKCOORDINATES cmd _manage_vector_entry \
                        select_vector_done $n $frame.t $grab $withdraw $format first "" 
                } else {                                                         
                    $self manage_buttons_set PICKCOORDINATES cmd _manage_vector_entry \
                        select_vector_done $n $frame.t $grab $withdraw $format second $1st_coord                 
                }                          
            }
            select_vector_done {
                lassign $args n new_frame grab withdraw 1stcoord 2ndcoord format help 
                if { [winfo exists $new_frame] } {
                    set p_new_frame [winfo parent $new_frame]
                    set w_win [winfo parent $p_new_frame]
                }
                destroy $p_new_frame
                if { $withdraw ne "" } {
                    wm deiconify $withdraw
                }
                if { $grab ne "" } {
                    grab $grab
                }
                if { $1stcoord ne "" && $2ndcoord ne "" } {
                    set vector [m::sub_vect $2ndcoord $1stcoord]
                    
                    if {$format != ""} {
                        set ivList ""
                        foreach iv $vector {
                            set iv [format $format $iv]
                            lappend ivList $iv
                        }
                        set vector $ivList
                    }
                    set entryvalues($n) [join $vector ","]               
                    set origin $1stcoord                                         
                    
                    set gid_groups_conds::origin $origin; set gid_groups_conds::vector $vector
                    
                    gid_groups_conds::register_drawopenGLvector 
                    gid_groups_conds::drawopenGLvector 
                    gid_groups_conds::end_register_drawopenGLvector                                   
                    $w_win.b1 configure -image [gid_groups_conds::giveimage remove-16] \
                        -command [mymethod end_vector_drawing $n $w_win $format]                     
                    gid_groups_conds::register_popup_help -check_bbox $w_win.b1 \
                        [= "Finish the drawing of the vector"]
                }
            }    
            select_point {
                lassign $args n frame num_pnt format 
                set grab ""
                if { [grab current $frame] ne "" } {
                    set grab [grab current $frame]
                    grab release $grab
                }
                if { [GiD_Info Project ViewMode] eq "GEOMETRYUSE"} {
                    set msg [= "Select a point in geometry"]
                } elseif {[GiD_Info Project ViewMode] eq "MESHUSE"} {
                    set msg [= "Select a node in mesh"]
                }
                if { [info command ::draw_post::give_toplevel] ne "" } { 
                    if { ([winfo toplevel $frame] ne [draw_post::give_toplevel]) && [draw_post::give_toplevel] != "" } {
                        wm withdraw [winfo toplevel $frame]
                        set withdraw [winfo toplevel $frame]                   
                        toplevel $frame.t                   
                        wm title $frame.t $msg
                        wm transient $frame.t [draw_post::give_toplevel]
                        set f [ttk::frame $frame.t.f]
                        grid $f -sticky nsew
                    } else {
                        set withdraw ""   
                        set grab ""                    
                        set f [ttk::frame $frame.t]
                        place $f -x 0 -y 0 -relwidth 1 -relheight 1 -anchor nw
                    }
                } else {
                    set withdraw ""   
                    set grab ""                    
                    set f [ttk::frame $frame.t]
                    place $f -x 0 -y 0 -relwidth 1 -relheight 1 -anchor nw
                }
                ttk::label $frame.t.l -text $msg
                ttk::button $frame.t.b -image [gid_groups_conds::giveimage close17] -command \
                    [mymethod _manage_vector_entry select_point_done $n $frame.t $grab $withdraw "" "" $format]
                grid $frame.t.l $frame.t.b -padx 2 -pady 2
                
                if {$num_pnt != ""} {
                    if { [draw_post::is_init_postprocess] } {
                        draw_post::manage_buttons_set_graphmode \
                            PICKNODE cmd [mymethod _manage_vector_entry select_num_point_done $n $frame.t \
                                $grab $withdraw "" "" $format]  
                    } else {                                      
                        $self manage_buttons_set PICKPOINT cmd _manage_vector_entry \
                            select_num_point_done $n $frame.t $grab $withdraw $format
                    }
                } else {
                    if { [draw_post::is_init_postprocess] } {
                        draw_post::manage_buttons_set_graphmode \
                            PICKNODE cmd [mymethod _manage_vector_entry select_point_done $n $frame.t \
                                $grab $withdraw "" "" $format]  
                    } else {                                      
                        $self manage_buttons_set PICKPOINT cmd _manage_vector_entry \
                            select_point_done $n $frame.t $grab $withdraw $format
                    }
                }
            }
            select_point_done {
                lassign $args n new_frame grab withdraw num point format
                destroy $new_frame               
                if { $withdraw ne "" } {
                    wm deiconify $withdraw
                }
                if { $grab ne "" } {
                    grab $grab
                }
                if { $point ne "" } {
                    if {$format != ""} {
                        set ipntList ""
                        foreach ipnt $point {
                            set ipnt [format $format $ipnt]
                            lappend ipntList $ipnt
                        }
                        set point $ipntList
                    }
                    set entryvalues($n) [join $point ","]
                }
            }
            select_num_point_done {
                lassign $args n new_frame grab withdraw num point format
                destroy $new_frame
                if { $withdraw ne "" } {
                    wm deiconify $withdraw
                }
                if { $grab ne "" } {
                    grab $grab
                }
                if { $num ne "" } {
                    set entryvalues($n) $num
                }
            }           
        }
    } 
    method end_vector_drawing { n frame format } {
        GiD_Process 'Redraw
        $frame.b1 configure -image [gid_groups_conds::giveimage vector-16] -style Toolbutton \
            -command [mymethod _manage_vector_entry select_vector $n $frame "" $format ""] 
        gid_groups_conds::register_popup_help -check_bbox $frame.b1 \
            [= "Choose a vector on screen"]
        
    }    
    method manage_buttons_set { graphmode args } {      
        
        lassign $args cmd - type n new_frame grab withdraw format 1stor2nd 1stcoord 
        
        switch $graphmode {
            "PICKPOINT" {
                set cursor dotbox
                if { [winfo exists .gid] } {
                    set wg .gid.central.s
                    $wg configure -cursor $cursor
                }                               
                lappend args cursor $cursor 
                if { [GiD_Info Project ViewMode] eq "GEOMETRYUSE"} {             
                    set num [GidUtils::PickEntities Points single]
                } elseif {[GiD_Info Project ViewMode] eq "MESHUSE"} {
                    set num [GidUtils::PickEntities Nodes single]
                }
                set point ""
                if {$num != ""} {
                    if { $type eq "select_num_point_done" } {                        
                        set point $num   
                    } else {
                        set point {*}[GiD_Info Coordinates $num]                   
                    }                    
                }
                if { [winfo exists .gid] } {
                    set wg .gid.central.s
                    set cursor ""
                    $wg configure -cursor $cursor
                } 
                if { $type eq "select_num_point_done" } {
                    $self _manage_vector_entry \
                        select_num_point_done $n $new_frame $grab $withdraw $num $point $format 
                } else {
                    $self _manage_vector_entry \
                        select_point_done $n $new_frame $grab $withdraw $num $point $format
                } 
            } 
            "PICKCOORDINATES" {
                if { $1stor2nd eq "first" } {
                    set msg [= "Enter origin point coordinates"]
                } else {
                    set msg [= "Enter destination point coordinates"]
                }
                if { [GiD_Info Project ViewMode] eq "GEOMETRYUSE"} {             
                    set coord [GidUtils::GetCoordinates $msg NoJoin  GEOMETRYUSE]                  
                } elseif {[GiD_Info Project ViewMode] eq "MESHUSE"} {
                    set coord [GidUtils::GetCoordinates $msg NoJoin  MESHUSE]
                }    
                GidUtils::SetWarnLine [_ "Selected point coordinates: %s" $coord]            
                if { [winfo exists .gid] } {
                    set wg .gid.central.s                   
                    $wg configure -cursor ""                    
                }                                 
                if { $1stor2nd eq "first" } {
                    $self _manage_vector_entry \
                        select_vector $n $new_frame $coord $format ""                    
                } elseif {$1stor2nd eq "second" } {                                   
                    $self _manage_vector_entry \
                        select_vector_done $n $new_frame $grab $withdraw $1stcoord $coord $format                   
                }                      
            }
        }
    }
    # main function to create the data windows
    method _fill_data_window_values { w parentNode valueNodes title } {
        set idx 0
        while { [winfo exists $w.f$idx] } { incr idx }
        set f [ttk::labelframe $w.f$idx -text $title]
        grid $f -sticky nswe -padx 2 -pady 1
        grid columnconfigure $w 0 -weight 1
        #grid rowconfigure $w $idx -weight 0
        #         # left to 0 to give priority to the notebook
        #         grid rowconfigure $w 0 -weight 1
        
        set entry_node_values($parentNode) \
            [list [list TLabelframe $f]]
        set state [get_domnode_attribute $parentNode state normal]
        if { $state eq "hidden" } {
            grid remove $f
        }
        set row 0
        set idx 0
        foreach valueNode $valueNodes {
            set state [get_domnode_attribute $valueNode state normal]
            set unit_magnitude [get_domnode_attribute $valueNode unit_magnitude]
            set unit_definition [get_domnode_attribute $valueNode unit_definition]
            set unit_mesh_definition [get_domnode_attribute $valueNode unit_mesh_definition]
            
            set pn [get_domnode_attribute $valueNode pn]
            set n [$valueNode @n]
            set help [$self _give_node_help $valueNode]
            set disable_on_void [get_domnode_attribute $valueNode disable_on_void 0]
            
            set help_img [$self _give_node_help_image $valueNode] 
            if {$help_img != ""} { 
                set help_image [gid_groups_conds::giveimage $help_img]
            } else {
                set help_image ""
            }
            
            if { !$disable_on_void } {
                ttk::label $f.l$idx -text $pn:
                set env [list [list label $f.l$idx]]
            } else {
                ttk::checkbutton $f.l$idx -text $pn: -variable [myvar entryvalues_disable_on_void($n)]
                set env [list [list checkbutton $f.l$idx]]
            }
            if { $help ne "" } {               
                gid_groups_conds::register_popup_help -image_tag $help_image $f.l$idx $help 
            }           
            set dict ""
            foreach "ni vi" [split [get_domnode_attribute $valueNode dict] ,] {
                dict set dict $ni [= $vi]
            }
            set values [split [get_domnode_attribute $valueNode values] ,]
            set value [get_domnode_attribute $valueNode v]
            if { [get_domnode_attribute $valueNode units_system_definition] == 1 } {
                foreach i [gid_groups_conds::give_units_system_list] {
                    lappend values [lindex $i 0]
                    dict set dict [lindex $i 0] [= [lindex $i 1]]
                }
                lassign [gid_groups_conds::give_active_units_system] value
            }
            set isproc 0
            set proc [$valueNode @values ""]
            if { $proc != ""} {                
                if { [string index $proc 0] eq "\[" && [string index $proc end] eq "\]" } {
                    set isproc 1
                }                 
            }
            if { [$valueNode nodeName] eq "edit_command" } {
                catch { $f.l$idx configure -text "" }
                set edit_type [get_domnode_attribute $valueNode edit_type]
                ttk::button $f.e$idx -text $pn \
                    -command [mymethod _edit_command_eval_window $w $valueNode $edit_type]
                grid $f.e$idx -sticky nw -padx 2 -pady "2 1"
                if { $help ne "" } {                
                    gid_groups_conds::register_popup_help $f.e$idx $help             
                }
            } elseif { $values ne "" || $isproc } {
                if { [llength $values] == 2 && [string is boolean [lindex $values 0]] &&
                    [string is boolean [lindex $values 1]] && [$valueNode @editable 0] == 0 } {
                    if { [lindex $values 0] } {
                        foreach "onvalue offvalue" $values break
                    } else {
                        foreach "offvalue onvalue" $values break
                    }
                    ttk::checkbutton $f.e$idx -text $pn -variable \
                        [varname entryvalues($n)] \
                        -onvalue $onvalue -offvalue $offvalue \
                        -command $options(-external_update_handler)
                    destroy $f.l$idx
                    set env ""
                    lappend env [list checkbutton $f.e$idx]
                    if { $help ne "" } {
                        gid_groups_conds::register_popup_help -check_bbox $f.e$idx $help
                    }
                } else {
                    cu::combobox $f.e$idx -textvariable \
                        [varname entryvalues($n)] \
                        -state readonly -values $values \
                        -dict $dict
                    bind $f.e$idx <<ComboboxSelected>> $options(-external_update_handler)
                    if { [$valueNode @editable 0] == 1 } {
                        $f.e$idx configure -state normal
                    }
                    lappend env [list ComboBox $f.e$idx]
                    if { $help ne "" } {
                        gid_groups_conds::register_popup_help -check_bbox $f.e$idx $help
                    }
                }
            } elseif { [get_domnode_attribute $valueNode values_tree ""] ne "" } {
                cu::combobox_tree $f.e$idx -textvariable \
                    [varname entryvalues($n)] \
                    -state readonly
                bind $f.e$idx <<ComboboxSelected>> $options(-external_update_handler)
                if { [$valueNode @editable 0] == 1 } {
                    $f.e$idx configure -state normal
                }
                set last_item(-1) root
                set values_tree [split [get_domnode_attribute $valueNode values_tree] ,]
                set _img_dict_ ""
                foreach i $values_tree {
                    lassign $i level name fname icon selectable command collapse
                    set parent $last_item([expr {$level-1}])
                    if { ![dict exists $_img_dict_ $icon] } {
                        dict set _img_dict_ $icon [gid_groups_conds::giveimage $icon]
                    }
                    if { $collapse eq "" } { set collapse 0 }
                    set image [dict get $_img_dict_ $icon]
                    #                     set image [gid_groups_conds::giveimage $icon]
                    set item [$f.e$idx tree_insert -image $image -collapse $collapse \
                            -active $selectable -command $command end \
                            $name $fname $parent]
                    set last_item($level) $item
                }
                if { [lsearch -exact [$f.e$idx cget -values] $value] == -1 } {
                    set value [$f.e$idx give_first_active]
                }
                lappend env [list ComboBox $f.e$idx]
                if { $help ne "" } {
                    gid_groups_conds::register_popup_help -check_bbox $f.e$idx $help
                }
            } elseif { [get_domnode_attribute $valueNode values_check ""] ne "" } {
                set values_check [split [get_domnode_attribute $valueNode values_check] ,]
                
                cu::check_listbox $f.e$idx -values $values_check -bd 1 -relief sunken
                
                lappend env [list Fulltktree $f.e$idx]
                
                set cmd [list $f.e$idx set_selected_comma_list -varname [myvar entryvalues($n)]]
                trace add variable [myvar entryvalues($n)] write "$cmd;#"
                bind $f.e$idx <Destroy> [list trace remove variable [myvar entryvalues($n)] write "$cmd;#"]
                
                set cmd [list $f.e$idx get_selected_comma_list -varname [myvar entryvalues($n)]]
                trace add variable [myvar entryvalues($n)] read "$cmd;#"
                bind $f.e$idx <Destroy> [list trace remove variable [myvar entryvalues($n)] read "$cmd;#"]
                
                if { $help ne "" } {
                    gid_groups_conds::register_popup_help -check_bbox $f.e$idx $help
                }
            } elseif { $unit_magnitude ne "" || $unit_definition ne "" || $unit_mesh_definition == 1 } {
                if { [llength [split $unit_magnitude ,]] > 1 } {
                    set ov_variable [myvar entryvalues(_selection_entity_)]
                } else {
                    set ov_variable ""
                }
                gid_groups_conds::entry_units $f.e$idx -value_node $valueNode \
                    -value_variable [myvar entryvalues($n)] \
                    -units_variable [myvar entryvalues_units($n)] \
                    -ov_variable $ov_variable
                
                if { $unit_definition eq "" && $unit_mesh_definition != 1 } {
                    set value $entryvalues($n)
                }
                
                if { [winfo width $win] < 220 } {
                    set width 5
                } else {
                    set width 7
                }
                $f.e$idx configure -units_width $width
                
                lappend env [list entry_units $f.e$idx]
                
                if { $help ne "" } {
                    gid_groups_conds::register_popup_help -check_bbox $f.e$idx $help
                }
            } elseif { [$valueNode @fieldtype ""] eq "long text" } {
                cu::multiline_entry $f.e$idx -textvariable \
                    [myvar entryvalues($n)] -height 4
                lappend env [list entry $f.e$idx]
                if { $help ne "" } {
                    gid_groups_conds::register_popup_help -check_bbox $f.e$idx $help
                }
            } elseif { [$valueNode @fieldtype ""] eq "vector" } {
                set dimensions [$valueNode @dimensions 3]
                set format [$valueNode @format ""]
                ttk::frame $f.e$idx
                set column 0
                set entry_node_values(vector,$valueNode) ""
                for { set i 1 } { $i <= $dimensions } { incr i } {
                    if { $i > 1 } {
                        ttk::label $f.e$idx.l$i -text ","
                        grid $f.e$idx.l$i -row 0 -column $column
                        incr column
                    }
                    set entryvalues($n,$i) ""
                    ttk::entry $f.e$idx.e$i -width 6 -textvariable [myvar entryvalues($n,$i)]
                    set cmd "[mymethod _manage_vector_entry from_components $n $dimensions];#"
                    trace add variable [myvar entryvalues($n,$i)] write $cmd
                    bind $f.e$idx.e$i <Destroy> [list trace remove variable [myvar entryvalues($n,$i)] write $cmd]
                    if { $help ne "" } {
                        gid_groups_conds::register_popup_help -check_bbox $f.e$idx.e$i $help
                    }
                    if { $dimensions == "1" } {
                        grid $f.e$idx.e$i -sticky ew
                        grid configure $f.e$idx.e$i -padx "0 2" -pady "0 2"
                    } else {                                            
                        grid $f.e$idx.e$i -row 0 -column $column -sticky ew
                    }
                    grid columnconfigure $f.e$idx $column -weight 1
                    incr column
                    lappend env [list entry $f.e$idx.e$i]                    
                    if { $i > 1 } {     
                        lappend entry_node_values(vector,$valueNode) [list label \
                                $f.e$idx.l$i]
                    }
                }    
                if { [$valueNode @pick_point 0] } {
                    if { $dimensions == "1" } {
                        ttk::button $f.e$idx.b1 -image [gid_groups_conds::giveimage node16x16] -style Toolbutton \
                            -command [mymethod _manage_vector_entry select_point $n $f.e$idx 1 $format]
                    } else {
                        ttk::button $f.e$idx.b1 -image [gid_groups_conds::giveimage node16x16] -style Toolbutton \
                            -command [mymethod _manage_vector_entry select_point $n $f.e$idx "" $format]
                    }
                    if { [GiD_Info Project ViewMode] eq "GEOMETRYUSE"} {
                        set msg [= "Select point in geometry"]
                    } elseif {[GiD_Info Project ViewMode] eq "MESHUSE"} {
                        set msg [= "Select node in mesh"]
                    }
                    gid_groups_conds::register_popup_help $f.e$idx.b1 $msg
                    grid $f.e$idx.b1 -row 0 -column $column
                    incr column
                    lappend entry_node_values(vector,$valueNode) [list button \
                            $f.e$idx.b1]
                } elseif { [$valueNode @pick_vector 0] } {                    
                    ttk::button $f.e$idx.b1 -image [gid_groups_conds::giveimage vector-16] -style Toolbutton \
                        -command [mymethod _manage_vector_entry select_vector $n $f.e$idx "" $format ""]
                    grid $f.e$idx.b1 -row 0 -column $column                                       
                    incr column
                    lappend entry_node_values(vector,$valueNode) [list button \
                            $f.e$idx.b1]
                    if { $help ne "" } {
                        gid_groups_conds::register_popup_help -check_bbox $f.e$idx.b1 $help
                    }
                }
                set cmd "[mymethod _manage_vector_entry to_components $n $dimensions];#"
                trace add variable [myvar entryvalues($n)] write $cmd
                bind $f.e$idx <Destroy> [list trace remove variable [myvar entryvalues($n,$i)] write $cmd]                                                                
            } else {
                ttk::entry $f.e$idx -textvariable \
                    [myvar entryvalues($n)]
                lappend env [list entry $f.e$idx]
                if { $help ne "" } {
                    gid_groups_conds::register_popup_help -check_bbox $f.e$idx $help
                }
            }
            if { [winfo children $f.e$idx] ne "" } {
                foreach i [winfo children $f.e$idx] {
                    bind $i <<Contextual>> [mymethod _fill_data_window_values_contextual $i \
                            $valueNode %X %Y]
                }
            } else {
                bind $f.e$idx <<Contextual>> [mymethod _fill_data_window_values_contextual $f.e$idx \
                        $valueNode %X %Y]
            }
            
            if { $disable_on_void } {
                set def [get_domnode_attribute $valueNode v_default 1]
                cu::enabled_disabled_entry $f.e$idx [myvar entryvalues_disable_on_void($n)] \
                    [myvar entryvalues($n)] "" $def $value
            } else {
                #                 if {[get_domnode_attribute $valueNode values ""] ne "" && \
                    #                     [$valueNode @editable 0] == 0} {                    
                    #                     set values [split [get_domnode_attribute $valueNode values] ,]
                    #                     if { [lsearch -exact $values $entryvalues($n)] == -1 } {
                        #                         set entryvalues($n) [lindex $values 0]
                        #                     }
                    #                 } else {                
                    #                 }
                set entryvalues($n) $value
                
            }
            set cmd [mymethod _fill_data_window_dependencies \
                    -value_variables [list $valueNode]]
            append cmd ";#"
            trace add variable [myvar entryvalues($n)] write $cmd
            bind $f <Destroy> [list trace remove variable \
                    [myvar entryvalues($n)] write $cmd]
            if { $state eq "disabled" } {
                if { [winfo exists $f.l$idx] } {
                    $f.l$idx configure -state disabled
                }
                if { [winfo children $f.e$idx] eq "" } {
                    $f.e$idx configure -state disabled
                } else {
                    foreach i [winfo children $f.e$idx] {
                        $i configure -state disabled
                    }
                }
            }
            set xpo [format_xpath {ancestor::condition/value[@n=%s]} [$valueNode @n]]
            append xpo [format_xpath {|ancestor::condition/container/value[@n=%s]} [$valueNode @n]]
            set valueNodeOrig [$valueNode selectNodes $xpo]
            if { $valueNodeOrig eq "" } { set valueNodeOrig $valueNode }
            if { [$valueNode @function ""] ne "" } {
                set funcNode [$valueNode selectNodes function]
                set funcs [cu::commalist_to_list [get_domnode_attribute $valueNode function]]
                if { [lsearch $funcs matrix_func] != -1 } { 
                    set function_image matrix-x-16
                } else {
                    set function_image eq-16
                }
                
                ttk::button $f.button$idx -style Toolbutton -image \
                    [gid_groups_conds::giveimage $function_image] \
                    -command [list gid_groups_conds::function_editor \
                        -units_var [myvar entryvalues_units($n)] \
                        -ov_var [myvar entryvalues(_selection_entity_)] \
                        $win [= "Create function for %s" $pn] \
                        $valueNode [myvar entryvalues_function($n)]]
                gid_groups_conds::register_popup_help $f.button$idx \
                    [= "Create/edit a function for this field"]
                
                lappend env [list button $f.button$idx]
                set cmd [mymethod _fill_data_window_function_update \
                        $funcs $n $f.button$idx $env]
                trace add variable [myvar entryvalues_function($n)] write "$cmd;#"
                bind $f <Destroy> +[list trace remove variable \
                        [myvar entryvalues_function($n)] write "$cmd;#"]
                
                if { $funcNode ne "" } {
                    set entryvalues_function($n) [$funcNode cloneNode -deep]
                } else {
                    set entryvalues_function($n) ""
                }
            } elseif { [$valueNodeOrig selectNodes edit_command] ne "" } {
                set editNodeList [$valueNodeOrig selectNodes edit_command]
                set inum 0
                foreach editNode $editNodeList {
                    if { [get_domnode_attribute $editNode icon ""] ne "" } {
                        set image [gid_groups_conds::giveimageL [get_domnode_attribute $editNode icon ""]]
                    } else {
                        set image [gid_groups_conds::giveimage run-16]
                    }
                    set edit_type [get_domnode_attribute $editNode edit_type]
                    set proc [get_domnode_attribute $editNode proc]
                    set cmd [mymethod _edit_command_eval_window \
                            -proc $proc $win $valueNode $edit_type]
                    if { $inum == 0 } {
                        if { [llength $editNodeList] == 1 } {
                            ttk::button $f.button$idx -style Toolbutton -image \
                                $image -command $cmd
                        } else {
                            ttk::menubutton $f.button$idx -style Toolbutton -image \
                                $image -menu $f.button$idx.m
                            menu $f.button$idx.m -tearoff 0
                            $f.button$idx.m add command -label [get_domnode_attribute $editNode pn] \
                                -command $cmd
                        }
                        set help [get_domnode_attribute $editNode pn]
                        gid_groups_conds::register_popup_help $f.button$idx \
                            $help
                        lappend env [list button $f.button$idx]
                    } elseif { $proc eq "" } {
                        $f.button$idx.m add separator
                    } else {
                        $f.button$idx.m add command -label [get_domnode_attribute $editNode pn] \
                            -command $cmd
                    }
                    incr inum
                }
            } elseif { [$valueNode @local_axes ""] ne "" } {
                cu::menubutton_button $f.button$idx -style Toolbutton -image local_axes-16 \
                    -menu $f.button$idx.m -command gid_groups_conds::local_axes_window
                menu $f.button$idx.m -tearoff 0
                gid_groups_conds::local_axes_menu -ov_variable \
                    [myvar entryvalues(_selection_entity_)] \
                    $f.button$idx.m
                
                set la_state [get_domnode_attribute $valueNode local_axes]
                if { $la_state eq "normal" } { set la_state !disabled }
                $f.button$idx state $la_state
                set help [= "Apply or draw local axes"]
                gid_groups_conds::register_popup_help $f.button$idx \
                    $help
                
                #lappend env [list button $f.button$idx]
                set entry_node_values(local_axes,$valueNode) [list [list button \
                            $f.button$idx]]
            }
            set xp {name(preceding-sibling::*[1])='label'}
            if { [$valueNodeOrig selectNodes $xp] } {
                set precedingNode [$valueNodeOrig selectNodes {preceding-sibling::*[1]}]
                set row [$self _fill_data_window_values_label $f $row $idx $precedingNode]
            }
            
            set entry_node_values($valueNode) $env
            
            set cspan 3
            set col 0
            if { [winfo exists $f.l$idx] } {
                grid $f.l$idx -row $row -column 0 -sticky w
                incr cspan -1
                incr col
            }
            if { [winfo exists $f.button$idx] } {
                grid $f.button$idx -row $idx -column 2 -sticky w
                incr cspan -1
            } else {
                set isbuttonprev 0
                foreach valueNode $valueNodes {
                    if { [$valueNode @fieldtype ""] eq "vector" } {
                        set isbuttonprev 1
                    }
                }                
                if {$isbuttonprev} {
                    ttk::label $f.lb$idx -text ""
                    grid $f.lb$idx -row $idx -column 2 -sticky w 
                    incr cspan -1                    
                }
            }
            
            if {  [$valueNode @fieldtype ""] eq "long text" } {
                set col 0
                incr row 1
                incr cspan
            }          
            grid $f.e$idx -sticky new -padx 3 -pady 2 -row $row \
                -column $col -columnspan $cspan            
            
            if { $state eq "hidden" } {
                foreach i [list $f.l$idx $f.e$idx $f.button$idx $f.lb$idx] {
                    if { [winfo exists $i] } {
                        grid remove $i
                    }
                }
            }
            set entryvalues_state($n) $state
            #                 if { [winfo exists $f.un$idx] } {
                #                     grid $f.un$idx -row 0 -column 2 -sticky nw
                #                 }
            incr idx
            incr row
        }
        if { [llength $valueNodes] } {
            set xp {name(following-sibling::*[1])='label'}
            if { [$valueNodeOrig selectNodes $xp] } {
                set followingNode [$valueNodeOrig selectNodes {following-sibling::*[1]}]
                set row [$self _fill_data_window_values_label $f $row $idx $followingNode]
            }
        }
        grid columnconfigure $f 1 -weight 1
        grid rowconfigure $f $row -weight 1
        
        return $f
    }
    method _fill_data_window_values_contextual { widget domNode x y } {
        
        focus $widget
        destroy $win.vmenu
        set m [menu $win.vmenu -tearoff 0]
        
        foreach i [list Cut Copy Paste] j [list [= "Cut"] [= "Copy"] [= "Paste"]] {
            $m add command -label $j -command [list event generate $widget <<$i>>]
        }
        $m add separator
        $m add command -label [= "Select all"] -command [list $widget selection range 0 end]
        
        if { [$domNode @function ""] ne "" } {
            set n [$domNode @n]
            set pn [get_domnode_attribute $domNode pn]
            $m add separator
            $m add command -label [= "Function editor"] -command \
                [list gid_groups_conds::function_editor \
                    -units_var [myvar entryvalues_units($n)] \
                    -ov_var [myvar entryvalues(_selection_entity_)] \
                    $win [= "Create function for %s" $pn] \
                    $domNode [myvar entryvalues_function($n)]]
        } elseif { [$domNode selectNodes edit_command] ne "" } {
            set editNode [$domNode selectNodes edit_command]
            set edit_type [get_domnode_attribute $editNode edit_type]
            set proc [get_domnode_attribute $editNode proc]
            $m add separator
            $m add command -label [= "Edit"] -command \
                [mymethod _edit_command_eval_window \
                    -proc $proc $win $domNode $edit_type]
        } elseif { [$domNode @local_axes ""] ne "" } {
            $m add separator
            $m add cascade -label [= "Local axes"] -menu $m.m1
            menu $m.m1 -tearoff 0
            gid_groups_conds::local_axes_menu -ov_variable \
                [myvar entryvalues(_selection_entity_)] \
                $m.m1
            
            set la_state [get_domnode_attribute $domNode local_axes]
            if { $la_state eq "normal" } { set la_state "" }
            $m entryconfigure end -state $la_state
        }
        tk_popup $m $x $y
    }
    method _fill_data_window_values_label { f row idx n } {
        
        set img [gid_groups_conds::create_give_image [$n toXPath] [$n @image ""]]
        set text [$n @pn ""]
        
        ttk::label $f.label$idx -text $text -image $img -compound \
            [$n @compound left] -justify left
        
        # this is approximate
        $self _fill_data_window_values_label_fill $f.label$idx [$self info_frame_width]
        bind $f.label$idx <Configure> [mymethod _fill_data_window_values_label_fill %W %w]
        grid $f.label$idx -row $row -column 0 -columnspan 3 -sticky ew
        incr row
        return $row
    }
    method _fill_data_window_values_label_fill { w width } {
        if { [$w cget -image] ne "" && [$w cget -compound] in "left right" } {
            set img_width [image width [$w cget -image]]
            if { $img_width > $width-40 } {
                set width [expr {$width-40}]
            } else {
                set width [expr {$width-$img_width}]
            }
        }
        if { $width < 40 } { set width 40 }
        $w configure -wraplength $width
    }
    method _fill_data_window_in_material_c { nb domNode menu } {
        set xp {ancestor::container[@n="materials"]/blockdata[@n="material"]|}
        append xp {ancestor::container[@n="materials"]/container/blockdata[@n="material"]}
        set idx [$nb index current]
        
        foreach containerNode [array names entry_node_values] {
            set v [lindex $entry_node_values($containerNode) 0]
            if { $v eq [list TNotebook $nb $idx] } {
                break
            }
        }
        set n [$containerNode @n]
        
        foreach bdNode [$domNode selectNodes $xp] {
            if { $bdNode eq $domNode } { continue }
            if { [$bdNode selectNodes [format_xpath {container[@n=%s]} $n]] eq "" } { continue }
            
            if { ![info exists menupf] } {
                $menu add cascade -label [= "Properties from"] -menu $menu.pf
                set menupf [menu $menu.pf -tearoff 0]
            }
            $menupf add command -label [$bdNode @name] -command \
                [mymethod _fill_data_window_in_material_c_get $nb $idx $containerNode $bdNode]
        }
        if { [llength [$nb tabs]] > 1 } {
            if { [$menu index end] ne "none" } { $menu add separator }
            $menu add command -label [= "Delete"] -command \
                [mymethod _fill_data_window_in_material_c_del $nb $idx $containerNode]
        }
    }
    method _unset_entry_node_values_subtree { node } {
        unset -nocomplain entry_node_values($node)
        foreach n [$node childNodes] {
            $self _unset_entry_node_values_subtree $n
        }
    }
    method _fill_data_window_in_material_c_get { nb idx containerNode bdNode } {
        
        set n [$containerNode @n]
        set containerNode_from [$bdNode selectNodes [format_xpath {container[@n=%s]} $n]]
        set state [get_domnode_attribute $containerNode state normal]
        set f $nb.f$idx
        destroy {*}[winfo children $f]
        
        $self _unset_entry_node_values_subtree $containerNode
        
        $self _fill_data_window $f $containerNode_from \
            [get_domnode_attribute $containerNode pn] 1
        
        set entry_node_values($containerNode_from) \
            [list [list TNotebook $nb $idx]]
        
        if { $state eq "hidden" } {
            $nb tab $idx -state $state
        }
    }
    method _fill_data_window_in_material_c_del { nb idx containerNode } {
        
        lappend entryvalues_mat_changes_stack [list delete [$containerNode @n]]
        set f $nb.f$idx
        destroy $f
        unset entry_node_values($containerNode)
        
        foreach containerNode [array names entry_node_values] {
            set v [lindex $entry_node_values($containerNode) 0]
            if { [lrange $v 0 1] eq [list TNotebook $nb] } {
                set i [lindex $v 2]
                if { $i > $idx } {
                    lset entry_node_values($containerNode) 0 2 [expr {$i-1}]
                }
            }
            for { set i 0 } { $i < [llength $entry_node_values($containerNode)] } { incr i } {
                set w [lindex $entry_node_values($containerNode) $i 1]
                if { [string match $f* $w] } {
                    set entry_node_values($containerNode) [lreplace \
                            $entry_node_values($containerNode) $i $i]
                    incr i -1
                }
            }
            if { [llength $entry_node_values($containerNode)] == 0 } {
                unset entry_node_values($containerNode)
            }
        }
    }
    method _fill_data_window_in_material_c_last { nb domNode menu } {
        
        set xp {ancestor::container[@n="materials"]/blockdata[@n="material"]/container|}
        append xp {ancestor::container[@n="materials"]/container/blockdata[@n="material"]/container}
        
        lassign "" cnts dict
        foreach cNode [$domNode selectNodes $xp] {
            set state [get_domnode_attribute $cNode state normal]
            if { $state eq "hidden" } { continue }
            dict lappend cnts [$cNode @n] [$cNode selectNodes string(../@name)]
            dict set dict [$cNode @n] [get_domnode_attribute $cNode pn]
        }
        set cns ""
        foreach containerNode [$domNode selectNodes container] {
            lappend cns [$containerNode @n]
        }
        foreach i $entryvalues_mat_changes_stack {
            lassign $i op a1
            switch $op {
                add {
                    lappend cns $a1
                }
                delete {
                    set ipos [lsearch $cns $a1]
                    set cns [lreplace $cns $ipos $ipos]
                }
            }
        }
        foreach c $cns {
            set cnts [dict remove $cnts $c]
        }
        set idx 0
        dict for "n v" $cnts {
            $menu add cascade -label [= "'%s' from" [dict get $dict $n]] -menu $menu.pf$idx
            set menupf [menu $menu.pf$idx -tearoff 0]
            foreach i $v {
                $menupf add command -label $i -command \
                    [mymethod _fill_data_window_in_material_c_get_tab $nb $domNode $n $i]
            }
            incr idx
        }
    }
    method _fill_data_window_in_material_c_get_tab { nb domNode n matname } {
        
        lappend entryvalues_mat_changes_stack [list add $n $matname]
        
        set xp0 {ancestor::container[@n="materials"]//blockdata[@n="material" and @name=%s]/}
        append xp0 {container[@n=%s]}
        set xp [format_xpath $xp0 $matname $n]
        set containerNode_from [$domNode selectNodes $xp]
        
        set state [get_domnode_attribute $containerNode_from state normal]
        set idx [llength [$nb tabs]]
        set f [ttk::frame $nb.f$idx]
        $nb add $f -text [string tolower [get_domnode_attribute \
                    $containerNode_from pn]] -sticky nsew
        $self _fill_data_window $f $containerNode_from \
            [get_domnode_attribute $containerNode_from pn] 1
        set entry_node_values($containerNode_from) \
            [list [list TNotebook $nb $idx]]
        if { $state eq "hidden" } {
            $nb tab $idx -state $state
        } else {
            $nb select $idx
        }
    }
    method _edit_command_eval_window_mod { args } {
        set optional {
            { -proc proc "" }          
        }       
        set compulsory "wmain domNode edit_type id"       
        parse_args $optional $compulsory $args
        $self _edit_command_eval_window -proc $proc $wmain $domNode $edit_type
    }
    method _edit_command_eval_window { args } {
        
        set optional {
            { -proc proc "" }          
        }       
        set compulsory "wmain domNode edit_type"       
        parse_args $optional $compulsory $args
        
        gid_groups_conds::uncompress_subtree $domNode
        
        set entries [array get entryvalues]
        set entries_units [array get entryvalues_units]
        
        if { $edit_type eq "callback" } {
            set callback [mymethod _edit_command_eval_window_update $unique_edit_id]
        } else {
            set callback ""
        }
        set ret [$self _edit_command_eval -dict $entries \
                -dict_units $entries_units -proc $proc -callback $callback \
                $domNode]
        if { $edit_type eq "callback" || $ret eq "" } { return }
        
        lassign $ret dict units_dict
        
        $self apply_values_ComboBox               
        $self _edit_command_eval_window_update $unique_edit_id $dict $units_dict
    }
    method _edit_command_eval_window_update { unique_id dict units_dict } {
        
        if { $unique_id != $unique_edit_id } {
            snit_messageBox -message [= "It is not possible to update. Out of sync"]
            return
        }      
        dict for "n v" $dict {
            if { [info exists entryvalues($n)] } {
                set entryvalues($n) $v
            }
            if { [info exists entryvalues_units($n)] } {
                set unitN [gid_groups_conds::units_to_nice_units \
                        [dict get $units_dict $n]]
                set entryvalues_units($n) $unitN
            }
        }
    }
    method _edit_command_eval { args } {
        
        set optional {
            { -item item "" }
            { -dict dict "" }
            { -dict_units dict_units "" }
            { -proc proc "" }
            { -callback callback "" }
            { -procmod boolean 0 }
        }
        set compulsory "domNode"
        parse_args $optional $compulsory $args
        
        gid_groups_conds::uncompress_subtree $domNode
        
        if { $proc eq "" } {
            set xp {edit_command|ancestor-or-self::condition/edit_command}
            set edit_commandNode [$domNode selectNodes $xp]
            if { $edit_commandNode eq "" } {
                set xp {ancestor-or-self::*/edit_command[@edit_type='exclusive' or @edit_type='exclusive_blockdata']}
                set edit_commandNode [$domNode selectNodes $xp]
            }          
            if { $edit_commandNode eq "" } {
                if { $procmod } {
                    set proc [$domNode @procmod ""]                
                } else {            
                    set proc [$domNode @proc ""]
                }
            } else {
                if { $procmod } {
                    set proc [$edit_commandNode @procmod ""]                
                } else {            
                    set proc [$edit_commandNode @proc ""]
                }
            }
        }
        if { $callback ne "" } {
            lappend proc $callback
        }        
        return [gid_groups_conds::eval_proc -tree $tree -boundary_conds \
                $self -item $item -dict $dict -dict_units $dict_units \
                $proc $domNode]
    }
    method _fill_data_window_dependenciesP { args } {
        
        set optional {
            { -is_pre_dependency boolean 0 }
        }
        set compulsory "domNode"
        parse_args $optional $compulsory $args
        
        set valueNodes [$domNode selectNodes {value|container/value|container/container/value}]
        if { ![llength $valueNodes] && [$domNode @isvalue 0] } {
            set valueNodes [list $domNode]
        }
        $self _fill_data_window_dependencies -value_variables -is_pre_dependency $is_pre_dependency $valueNodes
    }
    method _fill_data_window_dependencies { args } {
        
        set optional {
            { -value_variables "" 0 }
            { -condition_variables "" 0 }
            { -is_pre_dependency boolean 0 }
        }
        set compulsory "domNodes"
        parse_args $optional $compulsory $args
        
        set dirty 1
        set iterations 0
        while { $dirty && $iterations < 20 } {
            set dirty 0
            foreach domNode $domNodes {
                set save_values ""
                if { $value_variables } {
                    foreach valueNode [array names entry_node_values] {
                        if { [string first "," $valueNode] != -1 } { continue }
                        if { [$valueNode nodeName] ne "value" } { continue }
                        if { ![$valueNode hasAttribute n] } { continue }
                        set n [$valueNode @n]
                        if { [$valueNode hasAttribute v] } {
                            dict set save_values $n v [$valueNode @v]
                        }
                        $valueNode setAttribute v $entryvalues($n)
                        if { [string trim [$valueNode @units ""] ]ne "" } {
                            dict set save_values $n units [$valueNode @units]
                            $valueNode setAttribute units $entryvalues_units($n)
                        }
                        dict set save_values $n state [$valueNode @state ""]
                        $valueNode setAttribute state $entryvalues_state($n)
                        set funcNode [$valueNode selectNodes function]
                        if { $funcNode ne "" } {
                            dict set save_values $n func $funcNode
                            $valueNode removeChild $funcNode
                        }
                        if { [info exists entryvalues_function($n)] && \
                            $entryvalues_function($n) ne "" } {
                            $valueNode appendChild $entryvalues_function($n)
                        }
                    } 
                }
                if { $condition_variables } {
                    set ov [get_domnode_attribute $domNode ov ""]                    
                    if { $ov eq "" } {
                        set ov [get_domnode_attribute $domNode ov1 ""]
                        set ovwhat ov1
                    } else {
                        set ovwhat ov
                    }       
                    set isproc 0             
                    set ovval [$domNode @$ovwhat ""] 
                    if { [string index $ovval 0] eq "\[" && [string index $ovval end] eq "\]" } {
                        set isproc 1
                    }
                    if { [llength [split $ov ,]] > 1 } {
                        set save_ov $ov
                        if { !$isproc } {
                            $domNode setAttribute $ovwhat $entryvalues(_selection_entity_)
                        }
                    }
                }
                foreach depNode [gid_groups_conds::give_dependencies_nodes \
                        $domNode] {
                    
                    set domNodeFrom $domNode
                    
                    if {![info exists entryvalues(_selection_domNodeGrp_)]} {
                        set entryvalues(_selection_domNodeGrp_) ""
                    }
                    if { [$domNode nodeName] eq "condition" && $entryvalues(_selection_domNodeGrp_) ne "" } {
                        set domNodeFrom $entryvalues(_selection_domNodeGrp_)
                    }
                    set tomodNodes [$domNodeFrom selectNodes [$depNode @node]]
                    if { [llength $tomodNodes] == 0 } {
                        if { [$depNode @check 1] == 0 } { continue }
                        error "error in _fill_data_window_dependencies. xpath='[$depNode @node]'"
                    }
                    $depNode setAttribute v [get_domnode_attribute $domNode v ""]
                    if { [$domNode hasAttribute units] } {
                        $depNode setAttribute units [get_domnode_attribute $domNode units ""]
                    }
                    foreach tomodNode $tomodNodes {
                        if { ![info exists entry_node_values($tomodNode)] } { continue }
                        set n [$tomodNode @n]
                        lassign [list 0 0 0] has_activated has_v has_units
                        lassign "" ov_default local_axes state v units ovList
                        lassign "" vector
                        if {[$tomodNode hasAttribute fieldtype]} {
                            if { [$tomodNode @fieldtype] eq "vector" } {
                                set vector 1
                            }
                        }                            
                        if {[$depNode hasAttribute att]} {
                            set iiList [list 1]
                        } else {
                            set iiList [list 1 2]
                        }
                        foreach i $iiList {
                            if { ![$depNode hasAttribute att$i] && ![$depNode hasAttribute att]} { continue }
                            if { [$depNode hasAttribute att$i] } {
                                set depattNode [$depNode @att$i]
                            } elseif { [$depNode hasAttribute att] } {
                                set depattNode [$depNode @att]
                            }
                            switch -- $depattNode {
                                "state" {
                                    if { [$depNode hasAttribute att$i] } {
                                        set state [get_domnode_attribute $depNode v$i]
                                    } else {
                                        set state [get_domnode_attribute $depNode v]
                                    }
                                    #                                 if { [info exists entryvalues_state($n)] &&
                                        #                                     [lsearch "hidden disabled" $entryvalues_state($n)] != -1 &&
                                        #                                     [lsearch "hidden disabled" $state] == -1 } {
                                        #                                     set has_activated 1
                                        #                                 }
                                    set entryvalues_state($n) $state
                                }
                                "state1" {
                                    if { [$depNode hasAttribute att$i] } {
                                        set state1 [get_domnode_attribute $depNode v$i]
                                    } else {
                                        set state1 [get_domnode_attribute $depNode v]
                                    }
                                    set entryvalues(_selection_group1_state) $state1
                                }
                                "state2" {
                                    if { [$depNode hasAttribute att$i] } {
                                        set state2 [get_domnode_attribute $depNode v$i]
                                    } else {
                                        set state2 [get_domnode_attribute $depNode v]
                                    }
                                    set entryvalues(_selection_group2_state) $state2
                                }
                                "v" {
                                    set has_v 1
                                    if { [$depNode hasAttribute att$i] } {
                                        set v [get_domnode_attribute $depNode v$i]
                                    } else {
                                        set v [get_domnode_attribute $depNode v] 
                                    }
                                }
                                "units" {
                                    set has_units 1
                                    if { [$depNode hasAttribute att$i] } {
                                        set unit [get_domnode_attribute $depNode v$i]
                                    } else {
                                        set unit [get_domnode_attribute $depNode v]
                                    }
                                }
                                "ov" - "ov1" {
                                    if { [$depNode hasAttribute att$i] } {
                                        set ovList [split [get_domnode_attribute $depNode v$i] ,]
                                    } else {
                                        set ovList [split [get_domnode_attribute $depNode v] ,]
                                    }
                                }
                                "ov_default" - "ov1_default" {
                                    if { [$depNode hasAttribute att$i] } {
                                        set ov_default [get_domnode_attribute $depNode v$i]
                                    } else {
                                        set ov_default [get_domnode_attribute $depNode v]
                                    }
                                }
                                "local_axes" {
                                    if { [$depNode hasAttribute att$i] } {
                                        set local_axes [get_domnode_attribute $depNode v$i]
                                    } else {
                                        set local_axes [get_domnode_attribute $depNode v]
                                    }
                                }
                                change_units_system {
                                    if { !$is_pre_dependency } {
                                        set xp {self::*[@unit_definition or @unit_mesh_definition='1']}
                                        if { [$tomodNode selectNodes $xp] ne "" } {
                                            if { [$depNode hasAttribute att$i] } {
                                                set units_system [get_domnode_attribute $depNode v1]
                                            } else {
                                                set units_system [get_domnode_attribute $depNode v]
                                            }
                                            set unit [gid_groups_conds::nice_units_to_units \
                                                    $entryvalues_units($n)]
                                            
                                            lassign [gid_groups_conds::give_units_actualize_value \
                                                    -for_units_system $units_system \
                                                    -check_this_unit $unit \
                                                    -perform_the_change 0 $tomodNode] unit haschanges
                                            
                                            if { $haschanges } { set has_units 1 }
                                        }
                                    }
                                }
                            }
                        }
                        if { [$tomodNode nodeName] eq "value" && $has_v } {
                            if { $v ne $entryvalues($n) || $has_activated } {
                                set entryvalues($n) $v
                                if { [$tomodNode selectNodes dependencies] ne "" } {
                                    set dirty 1
                                }
                            }
                        }
                        if { [$tomodNode nodeName] eq "value" && $has_units } {
                            set unitN [gid_groups_conds::units_to_nice_units $unit]
                            if { $unitN ne $entryvalues_units($n) } {
                                set entryvalues_units($n) $unitN
                                if { [$tomodNode selectNodes dependencies] ne "" } {
                                    set dirty 1
                                }
                            }
                        }
                        #mylog::debug "domNode=[$domNode @n] tomodNode=[$tomodNode @n] state=$state"
                        if { $state eq "" } {
                            set st [dict_getd $save_values $n state ""] 
                            if { [regexp {^\s*\[} $st] } {
                                set state [_get_attribute_value $tomodNode $st]
                            } elseif { [info exists entryvalues_state($n)] } {
                                #set state [get_domnode_attribute $tomodNode state]
                                set state $entryvalues_state($n)
                            }
                            if { $state eq "" } { set state normal }
                            #mylog::debug "state=$state"
                            if { [info exists entryvalues_state($n)] && $entryvalues_state($n) ne $state } {
                                set entryvalues_state($n) $state
                                set dirty 1
                                #mylog::debug "dirty=1"
                            }
                        }
                        set nodes ""
                        lappend nodes $state $entry_node_values($tomodNode)
                        if { $vector ne "" } {
                            lappend nodes $state $entry_node_values(vector,$tomodNode)
                        }
                        if { $local_axes ne "" } {
                            lappend nodes $local_axes $entry_node_values(local_axes,$tomodNode)
                        }
                        if { [info exists state1] } {
                            lappend nodes $state1 [set! entry_node_values(_selection_group1_,$tomodNode)]
                        }
                        if { [info exists state2] } {
                            lappend nodes $state2 [set! entry_node_values(_selection_group2_,$tomodNode)]
                        }
                        foreach "stateL nodesL" $nodes {
                            foreach i $nodesL {
                                lassign $i type w
                                switch $type {
                                    TNotebook {
                                        set idx [lindex $i 2]
                                        if { $stateL eq "" } {set state normal }
                                        $w tab $idx -state $stateL
                                    }
                                    NoteBook {
                                        #                                     foreach "idx page" [lrange $i 2 3] break
                                        #                                     if { $stateL eq "normal" || $stateL eq "disabled" } {
                                            #                                         $w itemconfigure $page -state $stateL
                                            #                                         if { [$w index $page] == -1 } {
                                                #                                             $w insert $idx $page
                                                #                                         }
                                            #                                     } elseif { $stateL eq "hidden" } {
                                            #                                         if { [$w index $page] != -1 } {
                                                #                                             $w delete $page 0
                                                #                                         }
                                            #                                     }
                                    }
                                    TLabelframe {
                                        switch $stateL {
                                            hidden { grid remove $w }
                                            default { grid $w }
                                        }
                                    }
                                    label - ComboBox - entry - entry_units - checkbutton - button {
                                        if { $stateL ne "hidden" } {
                                            if { $type eq "entry_units" } {
                                                foreach j [winfo children $w] {
                                                    switch $stateL {
                                                        normal {
                                                            catch { $j state !disabled }
                                                        }
                                                        disabled {
                                                            catch { $j state disabled }
                                                        }
                                                        default {
                                                            catch { $j state $stateL }
                                                        }
                                                    }
                                                    #$j configure -state $stateL
                                                }
                                                if { [info exists units_system] } {
                                                    $w configure -units_system $units_system
                                                }
                                            } else {
                                                switch $stateL {
                                                    normal { $w state !disabled }
                                                    disabled { $w state disabled }                                                   
                                                    default { catch { $w state $stateL } }
                                                }
                                                #$w configure -state $stateL
                                            }
                                            grid $w
                                        } else {
                                            grid remove $w
                                        }
                                        
                                        if { $type eq "ComboBox" && \
                                            [get_domnode_attribute $tomodNode values ""] ne "" && \
                                            [$tomodNode @editable 0] == 0} {
                                            set values [split [get_domnode_attribute $tomodNode values] ,]
                                            if { $values ne [$w cget -values] } {
                                                $w configure -values $values
                                                if { [lsearch -exact $values $entryvalues($n)] == -1 } {
                                                    set entryvalues($n) [lindex $values 0]
                                                }
                                            }
                                        }
                                    }
                                    TFrame {
                                        if { $ovList ne "" } {
                                            $self _apply_cond_entities_frame $w $tomodNode "" \
                                                $ovList $ov_default
                                        }
                                    }
                                }
                            }
                        }
                        if { [info exists entryvalues_function($n)] } {
                            # to activate the trace
                            set entryvalues_function($n) $entryvalues_function($n)
                        }
                    }
                    if { [$depNode hasAttribute v] } {
                        $depNode removeAttribute v
                    }
                    if { [$depNode hasAttribute units] } {
                        $depNode removeAttribute units
                    }
                }
                if { $value_variables } {
                    foreach valueNode [array names entry_node_values] {
                        if { [string first "," $valueNode] != -1 } { continue }
                        if { ![$valueNode hasAttribute n] || ![$valueNode hasAttribute v] } { continue }
                        set n [$valueNode @n]
                        
                        if { [dict exists $save_values $n v] } {
                            $valueNode setAttribute v [dict get $save_values $n v]
                        } else {
                            $valueNode removeAttribute v
                        }
                        if { [dict exists $save_values $n units] } {
                            $valueNode setAttribute units [dict get $save_values $n units]
                        } elseif { [$valueNode hasAttribute units] } {
                            $valueNode removeAttribute units
                        }
                        $valueNode setAttribute state [dict get $save_values $n state]
                        
                        if { [info exists entryvalues_function($n)] && \
                            $entryvalues_function($n) ne "" } {
                            $valueNode removeChild $entryvalues_function($n)
                        }
                        if { [dict exists $save_values $n func] } {
                            $valueNode appendChild [dict get $save_values $n func]
                        }
                    }
                }
                if { $condition_variables } {
                    if { [info exists save_ov] } {
                        set ov [get_domnode_attribute $domNode ov ""]
                        if { $ov eq "" } {
                            set ov [get_domnode_attribute $domNode ov1 ""]
                            set ovwhat ov1
                        } else {
                            set ovwhat ov
                        }
                        set isproc 0             
                        set ovval [$domNode @$ovwhat ""] 
                        if { [string index $ovval 0] eq "\[" && [string index $ovval end] eq "\]" } {
                            set isproc 1
                        }
                        if {!$isproc} {
                            $domNode setAttribute $ovwhat $save_ov
                        }
                    }
                }
            }
            incr iterations
        }
    }
    method delete_cond_group_or_blockdata {} {
        
        set ids [$tree selection get]
        if { [llength $ids] == 0 } { return }
        set id [lindex $ids 0]
        set domNode [$tree item element cget $id 0 e_text -data]
        switch [$domNode nodeName] {
            blockdata { $self delete_disable_blockdata }
            default { $self delete_cond_group }
        }
    }
    method delete_cond_group {} {
        set ids [$tree selection get]
        
        set condition_nodes ""
        set group_nodes ""
        set condition_ids ""
        foreach id $ids {
            set domNode [$tree item element cget $id 0 e_text -data]
            if { [$domNode nodeName] eq "condition" } {
                eval lappend group_nodes [$domNode selectNodes group|groupList]
            } elseif { [$domNode selectNodes ancestor-or-self::groupList] ne "" } {
                lappend group_nodes [$domNode selectNodes ancestor-or-self::groupList]
            } elseif { [$domNode selectNodes ancestor-or-self::group] ne "" } {
                lappend group_nodes [$domNode selectNodes ancestor-or-self::group]
            }
            set cndNode [$domNode selectNodes ancestor-or-self::condition]
            if { $cndNode ne "" } {
                lappend condition_nodes $cndNode
                set ancestors_ids [list $id]
                eval lappend ancestors_ids [$tree item ancestors $id]
                foreach i $ancestors_ids {
                    if { [$tree item state get $i condition] } {
                        lappend condition_ids $i
                        break
                    }
                }
                set a0 [dict create add_group 1 edit_group 1 delete_group 1]
                set actions [dict merge $a0 [split [$cndNode @actions ""] ,]]
                if { ![dict get $actions delete_group] } {
                    set pn [get_domnode_attribute $cndNode pn]
                    snit_messageBox -message [= "It is not possible to delete groups from %s" $pn]
                    return
                }
            }
        }
        set group_nodes [lsort -unique $group_nodes]
        if { [llength $group_nodes] == 0 } { return }
        
        $self cancel_cond_or_blockdata_edit
        
        set condition_nodes [lsort -unique $condition_nodes]
        set condition_ids [lsort -unique $condition_ids]
        set cndName ""
        
        if { [llength $group_nodes] == 1 } {
            set domNode [lindex $group_nodes 0]
            if { [$domNode nodeName] eq "groupList" } {
                set ns ""
                foreach i [$domNode selectNodes group] {
                    lappend ns [$i @n]
                }
                set n [join $ns ,]
            } else {
                set n [$domNode @n]
            }
            set cndName [get_domnode_attribute [$domNode parentNode] pn]
            set txt [= "Are you sure to unassign group '%s' from condition '%s'?" $n $cndName]
        } else {
            set groups ""
            foreach domNode $group_nodes {
                if { [$domNode nodeName] eq "groupList" } {
                    set ns ""
                    foreach i [$domNode selectNodes group] {
                        lappend ns [$i @n]
                    }
                    set n ([join $ns ,])
                } else {
                    set n [$domNode @n]
                }
                lappend groups $n
            }
            set groupsP [join $groups ,]
            if { [string length $groupsP] > 100 } {
                set groupsP [string range $groupsP 0 96]
                regexp {^(.*),[^,]$} $groupsP {} groupsP
                append groupsP ...
            }
            if { [llength $condition_nodes] == 1 } {
                set cndName [get_domnode_attribute [lindex $condition_nodes 0] pn]
                set txt [= "Are you sure to unassign groups '%s' from condition '%s'?" \
                        $groupsP $cndName]
            } else {
                set txt [= "Are you sure to unassign groups '%s' from their condition?" \
                        $groupsP]
            }
        }
        set w [dialogwin_snit $win._ask -title [= "Unassign or delete"] -entrytext \
                $txt]
        set f [$w giveframe]
        
        ttk::radiobutton $f.cb1 -text [= "Keep the group/s"] -variable \
            [$w give_uservar deletegroups 0] -value 0
        ttk::radiobutton $f.cb2 -text [= "Fully delete the group/s (geometry/mesh entities not modified)"] -variable \
            [$w give_uservar deletegroups] -value 1
        
        grid $f.cb1 -sticky w -pady 2 -padx "5 60"
        grid $f.cb2 -sticky w -pady 2 -padx "5 60"
        
        bind $w <Return> [list $w invokeok]
        set action [$w createwindow]
        set deletegroups [$w give_uservar_value deletegroups]
        destroy $w
        if { $action < 1 } { return }
        
        if { $cndName ne "" } {
            set label [_ "Unassign groups from '%s'" $cndName]
        } else {
            se label [_ "Unassign groups"]
        }
        gid_groups_conds::start_group_undo_batch $label
        
        set groups ""
        foreach domNode $group_nodes {
            if { [$domNode nodeName] eq "groupList" } {
                foreach i [$domNode selectNodes group] {
                    lappend groups [$i @n]
                }
            } else {
                lappend groups [$domNode @n]
            }
            gid_groups_conds::deleteN $domNode
        }
        foreach id $condition_ids {
            $self actualize $id
            
            set pc [$tree item parent $id]
            while { $pc ne "" } {
                if { [$tree item state get $pc container] } {
                    set domNode [$tree item element cget $pc 0 e_text -data]
                    set xp {condition/group|condition/groupList}
                    if { [llength [$domNode selectNodes $xp]] } {
                        $tree item element configure $pc 0 e_text -font \
                            SystemBoldFont
                    } else {
                        $tree item element configure $pc 0 e_text -font ""
                    }
                }
                set pc [$tree item parent $pc]
            }
        }
        if { !$deletegroups } {
            gid_groups_conds::check_groups_auto_from_bc $groups
        } else {
            foreach group $groups {
                gid_groups_conds::delete_group -noactualizegroups $group
            }
            gid_groups_conds::actualize_groups_window
        }
        
        set actualize_nodesList ""
        foreach domNodeCnd $condition_nodes {
            set ret [gid_groups_conds::check_condition_dependencies1 $domNodeCnd]
            lappend actualize_nodesList {*}$ret
        }
        $self actualize_domNode {*}$actualize_nodesList
        gid_groups_conds::check_draw_nodes_state
        gid_groups_conds::end_group_undo_batch
    }
    method view_cond_group_in_group_window {} {
        
        set ids [$tree selection get]
        
        foreach id $ids {
            set domNode [$tree item element cget $id 0 e_text -data]
            if { [$domNode nodeName] ne "group" } { continue }
            if { $gid_groups_conds::is_draw_post } {
                draw_post::change_pane_state_name -focus groups tree
                set gr [ggc::manage_geometry_groups id [$domNode @n]]
                draw_post::create_groups_cmdT set_active_layer_num $gr
            } else {
                gid_groups_conds::open_groups .gid window_force
                gid_groups_conds::groups_cmd_if_exists select_group [$domNode @n]
            }
            break
        }
    }
    method rename_cond_group {} {
        set ids [$tree selection get]
        
        foreach id $ids {
            set domNode [$tree item element cget $id 0 e_text -data]
            if { [$domNode nodeName] ne "group" } { continue }
            gid_groups_conds::open_groups .gid window_force
            gid_groups_conds::groups_cmd_if_exists rename_group_start_name [$domNode @n]
            break
        }
    }
    proc recursively_change_state { w state } {
        catch { $w configure -state $state }
        catch { $w configure -disabledforeground grey60 }
        
        foreach i [winfo children $w] {
            recursively_change_state $i $state
        }
    }
    method _apply_cond_entities_frame { frame domNodeCnd domNodeGrp ovList ov_default } {
        
        foreach i [winfo children $frame] { destroy $i }
        set grp [lindex $domNodeGrp 0]
        
        switch [gid_groups_conds::is_geometry_or_mesh] {
            geometry {
                set ents [list [= Points] point point [= Lines] line line \
                        [= Surfaces] surface coonsurf [= Volumes] volume volume \
                        [= Volumes] surface_as_volume volume]
            }
            mesh {
                set ents [list [= Nodes] point node [= "L Elements"] line element \
                        [= "S Elements"] surface element \
                        [= "V Elements"] volume element]
            }
        }
        set dict ""
        foreach "ni vi" [split [get_domnode_attribute $domNodeCnd dict] ,] {
            dict set dict $ni [= $vi]
        }
        set idx 0
        foreach "name n img" $ents {
            if { [lsearch -exact $ovList $n] == -1 } { continue }
            if { [dict exists $dict $n] } { set name [dict get $dict $n] }
            ttk::radiobutton $frame.r$idx \
                -text [string tolower $name] \
                -image $::GidPriv(image_${img}24) -compound top \
                -style ToolbuttonT \
                -variable [varname entryvalues(_selection_entity_)] \
                -value $n
            gid_groups_conds::register_popup_help $frame.r$idx \
                [= "Apply the condition to %s" $name]
            grid $frame.r$idx -sticky w -row 0 -column $idx
            if { $grp ne "" && [$grp @ov ""] ne "" } {
                if { $n ne [$grp @ov] } {
                    $frame.r$idx state disabled
                }
            }
            incr idx
        }
        $self _apply_cond_selection_entity $domNodeGrp $ovList $ov_default
    }
    method _apply_cond_selection_entity { domNodeGrp ovList ov_default } {
        set grp [lindex $domNodeGrp 0]
        if { $grp ne "" && [$grp @ov ""] ne "" && [lsearch -exact $ovList [$grp @ov]] != -1 } {
            set entryvalues(_selection_entity_) [$grp @ov]
        } elseif { [lsearch -exact $ovList $ov_default] != -1 } {
            set entryvalues(_selection_entity_) $ov_default
        } else {
            set entryvalues(_selection_entity_) [lindex $ovList 0]
        }
    }
    # what can be apply or edit or select_more
    method apply_cond { args } {
        
        set optional {
            { -what apply|select_more|edit apply }
            { -show_buttons boolean 1 }
            { -open_window boolean 1 }
        }
        set compulsory ""
        parse_args $optional $compulsory $args
        
        if {!$open_window} { return }
        
        set ids [$tree selection get]
        if { [llength $ids] != 1 } {
            tk_messageBox -message [= "It is necessary to select one condition"]
            return
        }
        set id [lindex $ids 0]
        set cndid $id
        while { ![$tree item state get $cndid condition] } {
            set cndid [$tree item parent $cndid]
            if { $cndid == 0 } { return }
        }
        set domNodeCnd [$tree item element cget $cndid 0 e_text -data]
        set domNode [$tree item element cget $id 0 e_text -data]
        if { [$domNode selectNodes ancestor-or-self::groupList] ne "" } {
            set domNode [$domNode selectNodes ancestor-or-self::groupList]
            while { ![$tree item state get $id groupList] } {
                set id [$tree item parent $id]
                if { $id == 0 } { return }
            }
        }
        set domNodeGrp [$domNode selectNodes ancestor-or-self::groupList/group]
        if { $domNodeGrp eq "" } {
            set domNodeGrp [$domNode selectNodes ancestor-or-self::group]
        }
        if { $what eq "apply" && $domNodeGrp ne "" } {
            if { [$domNode selectNodes ancestor-or-self::groupList] ne "" } {
                set domNode [$domNode selectNodes ancestor-or-self::groupList]
            } else {
                set domNode $domNodeGrp
            }
        }
        if { $what eq "select_more"} {
            if { $domNodeGrp eq "" } { return }
            set domNode $domNodeGrp
        }
        if { $what eq "edit" && $domNodeGrp eq "" } { return }
        
        #         if { $what ne "apply" } {
            #             set entryvalues(_selection_domNodeGrp_) $domNodeGrp
            #         } else {
            #             set entryvalues(_selection_domNodeGrp_) ""
            #         }
        
        if { $what ne "apply" } {
            set entryvalues(_selection_domNodeGrp_) $domNodeGrp
        } 
        
        $self cancel_cond_edit
        #         $self cancel_cond_or_blockdata_edit
        
        if { [winfo exists $win.data] } { return }
        
        set pn [get_domnode_attribute $domNodeCnd pn]
        
        foreach i [list "" 1 2] {
            set entryvalues(_selection_created_group${i}_) ""
        }
        
        set a0 [dict create add_group 1 edit_group 1 delete_group 1]
        set actions [dict merge $a0 [split [$domNodeCnd @actions ""] ,]]
        if { $what eq "apply" && ![dict get $actions add_group] } {
            snit_messageBox -message [= "It is not possible to apply new groups to %s" $pn]
            return
        }
        if { $what eq "edit" && ![dict get $actions edit_group] } {
            snit_messageBox -message [= "It is not possible to edit group data for %s" $pn]
            return
        }
        
        set ov [get_domnode_attribute $domNodeCnd ov ""]
        set ov_default [get_domnode_attribute $domNodeCnd ov_default ""]
        if { $ov eq "" } {
            set ov [get_domnode_attribute $domNodeCnd ov1 ""]
            set ov_default [get_domnode_attribute $domNodeCnd ov1_default ""]
        }
        set ov0 $ov
        set dict ""
        foreach "ni vi" [split [get_domnode_attribute $domNodeCnd dict] ,] {
            dict set dict $ni [= $vi]
        }
        set ovList [split $ov ,]
        
        incr unique_edit_id
        
        ttk::frame $win.data
        bind $win.data <Destroy> [mymethod apply_cond_cancel]
        
        set pn [get_domnode_attribute $domNodeCnd pn]
        switch $what {
            apply {
                set title [= "Apply %s" $pn]
            }
            edit {
                set title [= "Edit %s" $pn]
            }
            select_more {
                set title [= "Select more %s" $pn]
            }
            default {
                set title [= "%s %s" $what $pn]
            }
        }
        
        if { [llength $ovList] > 1 } {
            $self _apply_cond_selection_entity $domNodeGrp $ovList $ov_default
        } else {
            set entryvalues(_selection_entity_) $ov
        }
        
        gid_groups_conds::start_group_undo_batch [_ "Apply condition '%s'" $pn]
        
        $self _fill_data_window $win.data $domNode $title
        
        if { [llength $ovList] > 1 } {
            ttk::frame $win.data.ents
            grid $win.data.ents -sticky nw
            $self _apply_cond_entities_frame $win.data.ents $domNodeCnd $domNodeGrp \
                $ovList $ov_default
            if { $domNodeGrp eq "" } {
                set entry_node_values($domNodeCnd) \
                    [list [list TFrame $win.data.ents]]
            }
        }
        
        #frame $win.data.fs -width 80 -height 2 -bd 4 -relief groove
        ttk::separator $win.data.fs -orient horizontal
        
        set groups [gid_groups_conds::give_groups]
        lappend groups ""
        
        set f [ttk::frame $win.data.fg]
        
        if { $what eq "select_more" } {
            recursively_change_state $win.data disabled
            recursively_change_state $f normal
        }
        
        set grp_help [= "Select an existing group or auto create a new one in order to apply the data to it"]
        set existing_ov ""
        foreach i [list "" 1 2] {
            set ov$i [get_domnode_attribute $domNodeCnd ov$i ""]
            if { [set ov$i] eq "" } { continue }
            lappend existing_ov $i
            set new_grp_name "[$domNodeCnd @pn] ${i}Auto"
            set ovp [get_domnode_attribute $domNodeCnd ov${i}p [= "Group"]]
            ttk::label $f.l$i -textvariable \
                [varname entryvalues(_selection_label_${i}_)]
            set txt [lindex [split $ovp ,] 0]
            set entryvalues(_selection_label_${i}_) "$txt:"
            gid_groups_conds::register_popup_help $f.l$i $grp_help
            
            set combo [cu::combobox_tree $f.cb$i -state readonly -textvariable \
                    [myvar entryvalues(_selection_group${i}_)]]
            gid_groups_conds::register_popup_help $f.cb$i $grp_help
            
            foreach group $groups {
                if { $group eq "" } {
                    set name [_ "(New group)"]
                } else {
                    set name $group
                }
                $combo tree_insert end $name $group 0
            }
            
            ttk::button $f.li$i -image $::GidPriv(image_select) -text [= "Select"] -compound left \
                -command [mymethod apply_cond_select_entities $i $domNodeCnd $f.li$i \
                    $combo $dict $what $new_grp_name]
            $f.li$i configure -width [expr {[string length [= "Select"]]-1}]
            gid_groups_conds::register_popup_help $f.li$i \
                [= "Create a new group and select entities into it"]
            
            grid $f.l$i $f.cb$i $f.li$i -sticky nw -padx 1
            grid configure $f.cb$i -sticky new
            grid columnconfigure $f 1 -weight 1
            
            set entry_node_values(_selection_group${i}_,$domNodeCnd) [list \
                    "label $f.l$i" "ComboBox $combo" "button $f.li$i"]
            
            set entryvalues(_selection_group${i}_state) "normal"
            
            if { $what eq "edit" || $what eq "select_more" } {
                if { $i eq "" } {
                    set grp $domNodeGrp
                } else {
                    set grp [lindex $domNodeGrp [expr {$i-1}]]
                }
                set entryvalues(_selection_group_is_new) 0
                set entryvalues(_selection_group${i}_) [$grp @n]
            } else {
                set entryvalues(_selection_group_is_new) 1
                set entryvalues(_selection_group${i}_) ""
            }
            switch $what {
                "edit" {
                    $f.l$i configure -state disabled
                    $f.cb$i configure -state disabled
                    $f.li$i configure -state disabled
                }
                "select_more" {
                    $f.cb$i configure -state disabled 
                }
                "apply" {
                    set vname [myvar entryvalues(_selection_group${i}_)]
                    set cmd [myproc deactivate_button_on_variable_void $vname $f.li$i]
                    trace add variable $vname write "$cmd;#"
                    set remove_cmd [list trace remove variable $vname write "$cmd;#"]
                    bind $f.cb$i <Destroy> $remove_cmd
                    $f.li$i configure -command "$remove_cmd ; [$f.li$i cget -command]"
                }
            }
            bind $combo <Return> "[list $win.data.b.b1 invoke] ; break"
        }
        if { [llength $ovList] > 1 } {
            set cmd [mymethod _apply_cond_change_sel_entity $domNodeCnd]
            set v [varname entryvalues(_selection_entity_)]
            trace add variable $v write "$cmd;#"
            bind  $win.data.ents.r0 <Destroy> [list trace remove variable \
                    $v write "$cmd;#"]
            uplevel #0 $cmd
        }
        
        $self _fill_data_window_dependenciesP -is_pre_dependency 1 $domNode
        $self _fill_data_window_reset_changes $win.data
        
        ttk::frame $win.data.b
        
        switch $what {
            select_more {
                set txt1 ""
                set txt2 [= Close]
            }
            default {
                set txt1 [= OK]
                set txt2 [= Cancel]
            }
        }
        
        set len 0
        foreach i [list $txt1 $txt2] {
            if { [string length $i] > $len } {
                set len [string length $i]
            }
        }
        if { $txt1 ne "" } {
            ttk::button $win.data.b.b1 -text $txt1 -image [gid_groups_conds::giveimage ok-16]\
                -compound left -width $len -command \
                [mymethod apply_cond_ok $existing_ov $cndid $id $what]
            gid_groups_conds::register_popup_help $win.data.b.b1 $grp_help
            
            set cmd [mymethod _apply_cond_change_group $existing_ov $win.data.b.b1]
            foreach i $existing_ov {
                set v [varname entryvalues(_selection_group${i}_)]
                trace add variable $v write "$cmd;#"
                bind $win.data.b.b1 <Destroy> [list trace remove variable \
                        $v write "$cmd;#"]
            }
            $self _apply_cond_change_group $existing_ov $win.data.b.b1
        }        
        ttk::button $win.data.b.b2 -text $txt2 -image [gid_groups_conds::giveimage remove-16] \
            -compound left -width $len -command [list destroy $win.data]
        
        if { [winfo exists $win.data.b.b1] } {
            grid $win.data.b.b1 $win.data.b.b2 -padx 2 -pady 3
        } else {
            grid $win.data.b.b2 -padx 2 -pady 3
        }
        
        if { $show_buttons } {
            grid $win.data.fs -sticky nwe -padx "2 40" -pady "10 3"
            grid $f -sticky new
            grid $win.data.b
        }
        grid $win.data -sticky ewns -padx 2 -pady 2
        grid columnconfigure $win.data 0 -weight 1
        
        switch $options(-edit_window_position) {
            up {
                grid $win.data - -sticky ewns -padx 2 -pady 2 -row 1
            }
            down {
                grid $win.data - -sticky ewns -padx 2 -pady 2
            }
            full {
                place $win.data -x 0 -y 0 -relwidth 1 -relheight 1 -anchor nw
            }
        }
        if { [winfo exists $win.data.f.e0] } {
            tk::TabToWindow $win.data.f.e0
        }
        $self correct_pane_height $win.data
        update
        if { [llength [$tree selection get]] } {
            $tree see [lindex [$tree selection get] 0]
        }
    }
    proc deactivate_button_on_variable_void { variablename button } {
        upvar #0 $variablename v
        if { $v eq "" } {
            $button state !disabled
        } else {
            $button state disabled
        }
    }
    method _apply_cond_set_group { ovi group } {
        set entryvalues(_selection_group${ovi}_) $group
    }
    method _apply_cond_change_sel_entity { domNodeCnd } {
        
        set ov [get_domnode_attribute $domNodeCnd ov ""]
        if { $ov eq "" } { set ov [get_domnode_attribute $domNodeCnd ov1 ""] }
        set ipos [lsearch [split $ov ,] $entryvalues(_selection_entity_)]
        
        foreach i [list "" 1 2] {
            set ov$i [get_domnode_attribute $domNodeCnd ov$i ""]
            if { [set ov$i] eq "" } { continue }
            set ovp [get_domnode_attribute $domNodeCnd ov${i}p [= "Group"]]
            set txt [lindex [split $ovp ,] $ipos]
            if { $txt ne "" } {
                set entryvalues(_selection_label_${i}_) "$txt:"
            }
        }
        $self _fill_data_window_dependencies -condition_variables \
            [list $domNodeCnd]
    }
    method _apply_cond_change_group { existing_ov button } {
        
        set isfull 1
        foreach i $existing_ov {
            if { $entryvalues(_selection_group${i}_state) in "disabled hidden" } { continue }
            if { [string trim $entryvalues(_selection_group${i}_)] eq "" } {
                set isfull 0
                break
            }
        }
        if { !$isfull } {
            $button state disabled
        } else {
            $button state !disabled
        }
    }
    method apply_cond_select_entities { ipos domNodeCnd button combo dict what group_name } {
        
        array set dictA $dict
        #         if { [string match *,* $ov] } {
            #             destroy $win.menu
            #             set m [menu $win.menu -tearoff 0]
            #             set ovList [split $ov ,]
            #             switch [lindex [GiD_Info Project] 4] {
                #                 GEOMETRYUSE {
                    #                     set ents [list [= Points] point point [= Lines] line line \
                        #                             [= Surfaces] surface coonsurf [= Volumes] volume volume \
                        #                             [= Volumes] surface_as_volume volume]
                    #                 }
                #                 MESHUSE {
                    #                     set ents [list [= Nodes] point node [= "L Elements"] line element \
                        #                             [= "S Elements"] surface element \
                        #                             [= "V Elements"] volume element]
                    #                 }
                #             }
            #             foreach "name n img" $ents {
                #                 if { [lsearch -exact $ovList $n] == -1 } { continue }
                #                 if { [info exists dictA($n)] } { set name $dictA($n) }
                #                 $m add command -label $name -image $::GidPriv(image_${img}16) \
                    #                     -compound left -command \
                    #                     [mymethod apply_cond_select_entities $ipos $button $combo \
                    #                         $n $dict $what $group_name]
                #             }
            #             set x [expr {[winfo rootx $button]}]
            #             set y [expr {[winfo rooty $button]+[winfo height $button]}]
            #             tk_popup $m $x $y
            #             return
            #         }
        
        set ov0 [get_domnode_attribute $domNodeCnd ov ""]
        if { $ov0 eq "" } { set ov0 [get_domnode_attribute $domNodeCnd ov1 ""] }
        set ovList [split $ov0 ,]
        set ov [get_domnode_attribute $domNodeCnd ov$ipos ""]
        if { [llength $ovList] > 1 } {
            set ov0 $entryvalues(_selection_entity_)
            set is [lsearch $ovList $ov0]
            set ov [lindex [split $ov ,] $is]
        }
        if { $what ne "select_more" } {
            lassign [gid_groups_conds::add_new_group $group_name] group
            set entryvalues(_selection_group${ipos}_) $group
            set entryvalues(_selection_created_group${ipos}_) $group
            $combo configure -values [gid_groups_conds::give_groups]
        } else {
            set group $entryvalues(_selection_group${ipos}_)
        }
        $combo state disabled
        if { [info exists dictA($ov)] } {
            set ovName $dictA($ov)
        } else { set ovName $ov }
        gid_groups_conds::apply_entities_to_group $group $ov $ovName $button $self \
            [mymethod end_apply_cond_select_entities $combo $button $what]
        
    }
    method end_apply_cond_select_entities { combo button what } {
        tk::TabToWindow $combo
        if { $what ne "select_more" } {
            $combo state "!disabled !readonly"
            $combo selection range 0 end
            $combo icursor end
        }
        $button state !disabled
        #         if { $what eq "select_more" } {
            #             destroy $win.data
            #         }
    } 
    
    method cancel_cond_edit {} {
        if { ![winfo exists $win.data] } { return }
        set group ""
        if { [set! entryvalues(_selection_group_is_new)] == 1 } {
            foreach i [array names entryvalues _selection_group] {
                set group [string trim $entryvalues($i)]
                if { $group ne "" } { break }
            }
        }
        if { $group ne "" || $entryvalues_haschanges } {
            set ret [snit_messageBox -type okcancel -default ok \
                    -message [= "Are you sure to dismiss the data that is being edited?"] \
                    -parent $win.data]
            if { $ret eq "cancel" } { return }
        }
        destroy $win.data
    }
    
    method cancel_cond_or_blockdata_edit {} {
        if { ![winfo exists $win.data] } { return }
        set group ""
        if { [set! entryvalues(_selection_group_is_new)] == 1 } {
            foreach i [array names entryvalues _selection_group*] {
                set group [string trim $entryvalues($i)]
                if { $group ne "" } { break }
            }
        }
        if { $group ne "" || $entryvalues_haschanges } {
            set ret [snit_messageBox -type okcancel -default ok \
                    -message [= "Are you sure to dismiss the data that is being edited?"] \
                    -parent $win.data]
            if { $ret eq "cancel" } { return }
        }
        destroy $win.data
    }
    method apply_cond_cancel {} {
        
        foreach i [array names entryvalues _selection_created_group*] {
            if { $entryvalues($i) ne "" } {
                gid_groups_conds::delete_group $entryvalues($i)
            }
        }
        foreach "n v" [array get entryvalues_function] {
            if { $v ne "" } { $v delete }
        }
        foreach i [list entryvalues entry_node_values entryvalues_units \
                entryvalues_state entryvalues_function] {
            array unset $i
        }
        set entryvalues_mat_changes_stack ""
        incr unique_edit_id
        gid_groups_conds::end_group_undo_batch
    }
    # what can be apply or edit
    method apply_cond_ok { existing_ov cndid id what } {
        
        set err [catch { $self apply_cond_ok_do $existing_ov $cndid $id $what } errstring]
        if { $err } {
            if { $errstring ne "" } {
                snit_messageBox -message $errstring -parent $win
            }
        }
    }
    method apply_cond_ok_do { args } {
        set optional {
            { -actualize_tree boolean 1 }
        }
        set compulsory "existing_ov cndid id what"
        parse_args $optional $compulsory $args        
        $tree item expand $id
        
        set ret [$self check_edited_values]
        if { $ret == -1 } { error "" }
        set domNodeCnd [$tree item element cget $cndid 0 e_text -data]
        set domNode [$tree item element cget $id 0 e_text -data]
        set groups [gid_groups_conds::give_groups]
        
        if { $what eq "apply" } {
            set groupList ""
            foreach i $existing_ov {
                set group [string trim $entryvalues(_selection_group${i}_)]
                if { $entryvalues(_selection_group${i}_state) in "disabled hidden" } {
                    if { $entryvalues(_selection_created_group${i}_) ne "" } {
                        gid_groups_conds::delete_group -noactualizegroups \
                            $entryvalues(_selection_created_group${i}_)
                    }
                    set group ""
                } elseif { $entryvalues(_selection_created_group${i}_) ne "" } {
                    if { $group eq "" } {
                        error [= "Cannot rename selected group to void"]
                    }
                    if { $group ne $entryvalues(_selection_created_group${i}_) } {
                        if { [lsearch -exact $groups $group] != -1 } {
                            error [= "Cannot rename group '%s' to '%s'. Already exists" \
                                    $entryvalues(_selection_created_group${i}_) $group]
                        }
                        gid_groups_conds::rename_group \
                            $entryvalues(_selection_created_group${i}_) \
                            $group
                        set ids [$tree selection get]
                        if { [llength $ids] != 1 } {
                            error [= "It is necessary to select one condition"]                            
                        }
                        set id [lindex $ids 0]
                        set cndid $id
                        while { ![$tree item state get $cndid condition] } {
                            set cndid [$tree item parent $cndid]
                            if { $cndid == 0 } { error "" }
                        } 
                        set domNodeCnd [$tree item element cget $cndid 0 e_text -data]
                        set domNode [$tree item element cget $id 0 e_text -data]
                    }
                } else {
                    if { $group eq "" } {
                        error [= "Select a group to apply the condition to it"] 
                    }
                    if { [lsearch -exact $groups $group] == -1 } {
                        error [= "Group '%s' does not exist" $group]
                    }
                    set gr [$domNodeCnd selectNodes [format_xpath {group[@n=%s]} $group]]
                    if { $gr ne "" } {
                        error [= "It is not possible to apply the same group several times to the same condition"]
                    }
                }
                lappend groupList $group
            }            
        }
        set isopen [$tree item state get $cndid open]
        set actualize_nodesList ""
        
        if { $what eq "apply" } {
            set ov0 [get_domnode_attribute $domNodeCnd ov ""]
            if { $ov0 eq "" } { set ov0 [get_domnode_attribute $domNodeCnd ov1 ""] }
            set ovList [split $ov0 ,]
            
            if { [llength $existing_ov] > 1 } {
                set gn [gid_groups_conds::addN $domNodeCnd groupList ""]
                set id_group [$self add $cndid [= "groupList"] \
                        [gid_groups_conds::giveimage kcmdf-16] groupList $gn]
                set idx 0
                foreach i $existing_ov {
                    set group [lindex $groupList $idx]
                    if { $group eq "" } { continue }
                    if { [llength $ovList] > 1 } {
                        set ov0 $entryvalues(_selection_entity_)
                        set is [lsearch $ovList $ov0]
                        set ov [get_domnode_attribute $domNodeCnd ov$i ""]
                        set ov [lindex [split $ov ,] $is]
                        set att [list n $group ov $ov]
                    } else {
                        set att [list n $group]
                    }
                    set g [gid_groups_conds::addN $gn group $att]
                    $self add $id_group [= "group: %s" $group] \
                        [gid_groups_conds::giveimage kcmdf-16] group $g
                    incr idx
                }
            } else {
                set ov [get_domnode_attribute $domNodeCnd ov ""]
                set ovList [split $ov ,]
                if { [llength $ovList] > 1 } {
                    set ov $entryvalues(_selection_entity_)
                    set att [list n $group ov $ov]
                } else {
                    set att [list n $group]
                }
                set gn [gid_groups_conds::addN $domNodeCnd group $att]
                set id_group [$self add $cndid [= "group: %s" $group] \
                        [gid_groups_conds::giveimage kcmdf-16] group $gn]
                $self item collapse $id_group
            }
            set ret [$self copy_cond_values $domNodeCnd $gn $id_group normal]
            lappend actualize_nodesList {*}$ret
            foreach group $groupList {
                if { $group eq "" } { continue }
                set newtype [gid_groups_conds::give_group_name_auto_from_bc \
                        $domNodeCnd]
                if { $newtype ne "" } {
                    gid_groups_conds::change_group_att $group type $newtype
                }
            }
        } elseif { $what eq "edit" } {
            lappend actualize_nodesList {*}[$self apply_edited_values $id]
        } else {
            error "error in apply_cond_ok"
        }
        if { !$isopen } {
            $tree item collapse $id
        }
        
        set p $cndid
        $tree item element configure $cndid 0 e_text -font \
            SystemBoldFont
        
        foreach i [array names entryvalues _selection_created_group*_] {
            set entryvalues($i) ""
        }
        
        set pc [$tree item parent $p]
        while { $pc ne "" } {
            if { [$tree item state get $pc container] } {
                $tree item element configure $pc 0 e_text -font \
                    SystemBoldFont
            }
            set pc [$tree item parent $pc]
        }
        while { $p != 0 && ![$tree item state get $p blockdata] } {
            set p [$tree item parent $p]
        }
        if { $p != 0 } {
            set domNodeBD [$tree item element cget $p 0 e_text -data]
            if { [$domNodeBD @state normal] eq "disabled" &&
                [$domNodeBD @sequence_type any] eq "non_void_disabled" } {
                gid_groups_conds::setAttributesN $domNodeBD [list state normal]
                lappend actualize_nodesList $domNodeBD
            } elseif { [$domNodeBD @active 1] == 0 } {
                gid_groups_conds::setAttributesN $domNodeBD [list active 1]
                if { [$domNodeBD @sequence_type any] eq "non_void_deactivated" } {
                    lappend actualize_nodesList $domNodeBD
                }
            }
        }
        
        destroy $win.data
        
        set ret [gid_groups_conds::check_condition_dependencies1 $domNodeCnd]
        lappend actualize_nodesList {*}$ret
        
        if { $actualize_tree } {
            $self actualize_domNode {*}$actualize_nodesList
        }
        gid_groups_conds::check_draw_nodes_state
        
        set update_proc [$domNode @update_proc ""]
        if {$update_proc != ""} {
            if { [string index $update_proc 0] eq "\[" && [string index $update_proc end] eq "\]" } {
                set update_proc [string range $update_proc 1 end-1]
            } 
            set ret [gid_groups_conds::eval_proc $update_proc $domNode]
            if { $ret == -1 } { return }            
        }
    }
    method apply_blockdata { args } {  
        variable execute_proc_cancel      
        set optional {
            { -del_cancel_button boolean 0 }
            { -open_window boolean 1 }
            { -create_new boolean 0}
        }
        set compulsory ""
        parse_args $optional $compulsory $args
        
        if {!$open_window} { return }
        if { $create_new } {  set execute_proc_cancel 1 }
        
        set ids [$tree selection get]
        if { [llength $ids] != 1 } {
            tk_messageBox -message [= "It is necessary to select a blockdata item"]
            return
        }
        set id [lindex $ids 0]
        
        set domNode [$tree item element cget $id 0 e_text -data]
        gid_groups_conds::uncompress_subtree $domNode
        
        $self cancel_cond_or_blockdata_edit
        if { [winfo exists $win.data] } { return }
        
        incr unique_edit_id
        
        set pn [get_domnode_attribute $domNode pn]
        gid_groups_conds::start_group_undo_batch [_ "Edit data '%s'" $pn]
        
        destroy $win.data
        ttk::frame $win.data
        bind $win.data <Destroy> [mymethod apply_blockdata_cancel $id]
        
        set ret [$self _fill_data_window $win.data $domNode $pn]
        set notebook [lindex $ret 0]
        $self _fill_data_window_dependenciesP -is_pre_dependency 1 $domNode
        $self _fill_data_window_reset_changes $win.data
        
        set len 0
        foreach i [list [= OK] [= Cancel] [= More]...] {
            if { [string length $i] > $len } {
                set len [string length $i]
            }
        }
        ttk::frame $win.data.b
        ttk::button $win.data.b.b1 -text [= OK] -image [gid_groups_conds::giveimage ok-16]\
            -compound left -width $len -command \
            [mymethod apply_blockdata_ok $id]
        if { !$del_cancel_button } {
            ttk::button $win.data.b.b2 -text [= Cancel] -image [gid_groups_conds::giveimage remove-16]\
                -compound left -width $len -command [list destroy $win.data]
        }
        if { [$domNode @n] eq "material" } {
            set cmd [list $notebook configure -menu_callback \
                    [mymethod _fill_data_window_in_material_c $notebook $domNode] \
                    -last_menu_callback \
                    [mymethod _fill_data_window_in_material_c_last $notebook $domNode]]
            append cmd ";" [list $win.data.b.b3 state disabled]
            if {[$domNode @morebutton 1]} {
                ttk::button $win.data.b.b3 -text [= More]... -image [gid_groups_conds::giveimage up-16]\
                    -compound left -width $len -command $cmd
            }
        }
        
        if { $del_cancel_button } {
            if { [winfo exists $win.data.b.b3] } {
                grid $win.data.b.b1 $win.data.b.b3 -padx 2 -pady 3
            } else {                    
                grid $win.data.b.b1 -padx 2 -pady 3
            }
        } else {
            if { [winfo exists $win.data.b.b3] } {
                grid $win.data.b.b1 $win.data.b.b2 $win.data.b.b3 -padx 2 -pady 3
            } else {                    
                grid $win.data.b.b1 $win.data.b.b2 -padx 2 -pady 3
            }                        
        }
        grid $win.data - -sticky new
        grid $win.data.b
        grid columnconfigure $win.data 0 -weight 1
        
        switch $options(-edit_window_position) {
            up {
                grid $win.data -sticky ewns -padx 2 -pady 2 -row 1
            }
            down {
                grid $win.data -sticky ewns -padx 2 -pady 2
            }
            full {
                place $win.data -x 0 -y 0 -relwidth 1 -relheight 1 -anchor nw
            }
        }
        if { [winfo exists $win.data.e0] } {
            tk::TabToWindow $win.data.f.e0
        }
        $self correct_pane_height $win.data
    }
    method correct_pane_height { data_window } {
        
        if { $gid_groups_conds::is_draw_post } {
            update idletasks
            set d [dict create [winfo parent $win] [expr {[winfo reqheight $data_window]+35}]]
            draw_post::correct_panes_heights -min_heights_dict $d
        }
    }
    method apply_blockdata_cancel { id } {
        variable execute_proc_cancel
        
        if { ![winfo exists $tree] } { return }
        set p $id
        if {$p in [$tree item range 0 end]} {
            while { $p != 0 && ![$tree item state get $p condition] && \
                ![$tree item state get $p blockdata] && \
                ![$tree item state get $p container] } {
                set p [$tree item parent $p]
            }
        } else {
            return
        }
        if { $p != 0 } {
            set domNodeBD [$tree item element cget $p 0 e_text -data]
            set proc_cancel [$domNodeBD @proc_cancel ""]
            if {$proc_cancel != "" } { 
                if { $execute_proc_cancel == "1" } {
                    if { [string index $proc_cancel 0] eq "\[" && [string index $proc_cancel end] eq "\]" } {
                        set proc_cancel [string range $proc_cancel 1 end-1]
                        set ret [gid_groups_conds::eval_proc $proc_cancel $domNodeBD]               
                    }        
                } else {
                    gid_groups_conds::actualize_conditions_window
                }
            }
        }
        
        foreach "n v" [array get entryvalues_function] {
            if { $v ne "" } { $v delete }
        }
        foreach i [list entryvalues entry_node_values entryvalues_units \
                entryvalues_state entryvalues_function] {
            array unset $i
        }
        set entryvalues_mat_changes_stack ""
        incr unique_edit_id
        gid_groups_conds::end_group_undo_batch
    }
    method apply_blockdata_ok { id } {
        variable execute_proc_cancel
        
        set ret [$self check_edited_values]
        if { $ret == -1 } { return }
        
        set p $id        
        while { $p != 0 && ![$tree item state get $p condition] && \
            ![$tree item state get $p blockdata] && \
            ![$tree item state get $p container] } {
            set p [$tree item parent $p]
        }        
        set domNodeBD [$tree item element cget $p 0 e_text -data]
        
        set isopen [$tree item state get $p open]
        
        set actualize_nodesList [$self apply_edited_values $p]
        
        if { !$isopen } {
            $tree item collapse $p
        }
        
        set check_values [$domNodeBD @check_values ""]
        if {$check_values != ""} {
            if { [string index $check_values 0] eq "\[" && [string index $check_values end] eq "\]" } {
                set check_values [string range $check_values 1 end-1]
                set ret [gid_groups_conds::eval_proc $check_values $domNodeBD]
                if { $ret == -1 } { return }
            }        
        }
        
        if { [$domNodeBD @actualize_parents 0] == "1" } {
            lappend actualize_nodesList [[$domNodeBD parentNode] parentNode]
        }       
        if { [$domNodeBD @state normal] eq "disabled" &&
            [$domNodeBD @sequence_type any] eq "non_void_disabled" } {
            gid_groups_conds::setAttributesN $domNodeBD [list state normal]
            lappend actualize_nodesList $domNodeBD
        } elseif { [$domNodeBD @active 1] == 0 } {
            gid_groups_conds::setAttributesN $domNodeBD [list active 1]
            #[$domNodeBD @sequence_type any] eq "non_void_deactivated"
            lappend actualize_nodesList $domNodeBD
        }
        $self actualize_domNode {*}$actualize_nodesList
        
        if { [$domNodeBD @actualize_tree 0] == "1" } {
            $self actualize
        }
        set execute_proc_cancel 0
        destroy $win.data  
        
        set update_proc [$domNodeBD @update_proc ""]
        if {$update_proc != ""} {
            if { [string index $update_proc 0] eq "\[" && [string index $update_proc end] eq "\]" } {
                set update_proc [string range $update_proc 1 end-1]
            } 
            set ret [gid_groups_conds::eval_proc $update_proc $domNodeBD]
            if { $ret == -1 } { return }            
        }     
    }
    #     method move_block_data { where } {
        # 
        #         set ids [$tree selection get]
        #         if { [llength $ids] == 0 } {
            #             tk_messageBox -message [= "It is necessary to select a blockdata item"]
            #             return
            #         }
        #         if { $where eq "bottom" || $where eq "down" } {
            #             foreach "tmpids ids" [list $ids ""] break
            #             foreach i $tmpids {
                #                 set ids [linsert $ids 0 $i]
                #             }
            #         }
        #         set domNode0 [$tree item element cget [lindex $ids 0] 0 e_text -data]
        #         set leveldomNodes [$domNode0 selectNodes ../*]
        # 
        #         switch $where {
            #             up {
                #                 set before_after before
                #                 set ipos [lsearch -exact $leveldomNodes $domNode0]
                #                 if { $ipos == 0 } { return }
                #                 set domNodeRef [lindex $leveldomNodes [expr {$ipos-1}]]
                #             }
            #             top {
                #                 set before_after before
                #                 set ipos [lsearch -exact $leveldomNodes $domNode0]
                #                 if { $ipos == 0 } { return }
                #                 set domNodeRef [lindex $leveldomNodes 0]
                #             }
            #             down {
                #                 set before_after after
                #                 set ipos [lsearch -exact $leveldomNodes $domNode0]
                #                 if { $ipos == [llength $leveldomNodes]-1 } { return }
                #                 set domNodeRef [lindex $leveldomNodes [expr {$ipos+1}]]
                #             }
            #             bottom {
                #                 set before_after after
                #                 set ipos [lsearch -exact $leveldomNodes $domNode0]
                #                 if { $ipos == [llength $leveldomNodes]-1 } { return }
                #                 set domNodeRef [lindex $leveldomNodes end]
                #             }
            #         }
        #         foreach id $ids {
            #             set domNode [$tree item element cget $id 0 e_text -data]
            #             gid_groups_conds::moveNode [nice_xpath $domNodeRef] \
                #                 [nice_xpath $domNode] $before_after
            #         }
        #         $self actualize
        #     }
    method copy_block_data { args } {
        
        set optional {
            { -copy_cond_groups "" 0 }
            { -domNode node "" }
            { -newname name "" }
        }
        set compulsory ""
        parse_args $optional $compulsory $args
        
        if { $domNode eq "" } {
            set ids [$tree selection get]
            if { [llength $ids] != 1} {
                snit_messageBox -message \
                    [= "It is necessary to select one blockdata item"]
                return
            }
            set id [lindex $ids 0]
            set domNode [$tree item element cget $id 0 e_text -data]
        } else {
            set id ""
        }
        
        gid_groups_conds::uncompress_subtree $domNode
        
        if { [$domNode nodeName] eq "container" } {            
            set xp {blockdata[@sequence='1']}
            set bdList [$domNode selectNodes $xp] 
            if { $bdList eq "" } { return }            
            
            if { $id ne "" } {
                $self expand_node $id
                $tree item expand $id
                foreach id_child [$tree item children $id] {
                    foreach bNode $bdList {
                        if { [$tree item element cget $id_child 0 e_text -data] eq $bNode } {
                            set id $id_child
                            set domNode $bNode
                            break
                        }
                    }
                }
            }
        } 
        if { ![$domNode @sequence 0] } {
            snit_messageBox -message [= "Blockdata item cannot be copied"]
            return    
        }
        
        set actualize_nodesList ""
        if { [$domNode @state normal] eq "disabled" &&
            [$domNode @sequence_type any] eq "non_void_disabled" } {
            gid_groups_conds::setAttributesN $domNode [list state normal]
            lappend actualize_nodesList $domNode
        } elseif { [$domNode @active 1] == 0 } {
            # [$domNode @sequence_type any] eq "non_void_deactivated"
            gid_groups_conds::setAttributesN $domNode [list active 1]
            lappend actualize_nodesList $domNode
        }
        set newNode [gid_groups_conds::copyNodeN $domNode [$domNode parentNode]]
        
        if { !$copy_cond_groups } {
            set xp {condition/group|condition/groupList|}
            append xp {container/condition/group|} \
                {container/condition/groupList}
            foreach d [$newNode selectNodes $xp] {
                gid_groups_conds::deleteN $d
            }
        }
        gid_groups_conds::setAttributesN $newNode [list state normal tree_state close]
        if { $newname ne "" } {
            gid_groups_conds::setAttributesN $newNode [list name $newname]
            gid_groups_conds::actualize_conditions_window
        } elseif { [$domNode @editable_name ""] eq "unique" } {
            set values ""
            foreach node [$domNode selectNodes "../blockdata"] {
                if { $node eq $newNode } { continue }
                lappend values [get_domnode_attribute $node name]
            }
            set basename [get_domnode_attribute $newNode name]
            set newname $basename
            if { ![regexp {(.*[^\d])(\d+)$} $newname {} basename idx] } {
                set idx ""
            }
            while { [lsearch -exact $values $newname] != -1 } {
                if { $idx eq "" } { set idx 2 } else { incr idx }
                set newname $basename$idx
            }
            if { $newname ne $basename } {
                gid_groups_conds::setAttributesN $newNode [list name $newname]
            }
        }
        
        if { $id ne "" } {
            set name [get_domnode_attribute $newNode name]
            if { $name eq "" } {
                set name [get_domnode_attribute $newNode pn]
            }
            set tree_parent [$tree item parent $id]
            set image [$self _give_icon $newNode]
            set newid [$self add $tree_parent $name \
                    $image blockdata \
                    $newNode]
            if { [$newNode @sequence 0] == 1 } {
                $tree item element configure $newid 0 e_text \
                    -font SystemBoldFont
            }         
            if { [get_domnode_attribute $newNode icon ""] ne "" } {
                set image [gid_groups_conds::giveimageL \
                        [get_domnode_attribute $newNode icon ""]]
                $tree item element configure $newid 0 e_image -image $image
                if {[get_domnode_attribute $newNode icon_end ""] ne ""} {
                    set f_image [gid_groups_conds::giveimageL \
                            [get_domnode_attribute $newNode icon_end ""]]
                    $tree item element configure $newid 0 f_image -image $f_image
                }
            }
            
            $tree item configure $newid -button yes
            $tree item collapse $newid
        } else {
            set newid ""
        }       
        
        set ret [gid_groups_conds::check_condition_dependencies1 "" $newNode]
        lappend actualize_nodesList {*}$ret
        $self actualize_domNode {*}$actualize_nodesList
        
        if { $id ne "" } {
            $tree selection clear all
            $tree selection add $newid
            $tree activate $newid
            $tree see $newid
            focus $tree
            update
            event generate $tree <F2>
        }
        
        if { [set proc [$domNode @update_proc ""]] ne "" } {
            gid_groups_conds::eval_proc -tree $tree -boundary_conds \
                $self -item $newid $proc $newNode
        }
        #         $self actualize
    }
    method delete_disable_blockdata {} {
        
        set ids [$tree selection get]
        if { [llength $ids] == 0 } {
            tk_messageBox -message [= "It is necessary to select a blockdata item"]
            return
        }
        
        $self cancel_cond_or_blockdata_edit
        
        set update_data ""
        # delete from end to beginning so as to let the first disabled
        set nids ""
        foreach id $ids {
            set domNode [$tree item element cget $id 0 e_text -data]
            if { [$domNode nodeName] ne "blockdata" || ![$domNode @sequence 0] } {
                continue
            }
            gid_groups_conds::uncompress_subtree $domNode
            
            set before_update_proc [$domNode @before_update_proc ""]
            if { $before_update_proc ne "" } {
                set ret [catch { gid_groups_conds::eval_proc -tree $tree -boundary_conds \
                            $self -item $id $before_update_proc $domNode } errstring]
                if { $ret } {
                    if { $errstring ne "" } {
                        tk_messageBox -message $errstring
                    }
                    return
                }
            }
            set nids [linsert $nids 0 $id]
        }
        if { [llength $nids] == 0 } {
            snit_messageBox -message [= "Some blockdata items cannot be deleted"] \
                -parent $win
        }
        set confirm 1
        foreach id $nids {
            set domNode [$tree item element cget $id 0 e_text -data]
            if { $update_data eq "" && [$domNode @update_proc ""] ne "" } {
                lappend update_data [$domNode @update_proc] \
                    [$domNode parentNode] [$tree item parent $id]
            }
            if { [$domNode @n] eq "material" } {
                set xp {ancestor::container[@n="materials"]/blockdata[@n="material"]|}
                append xp {ancestor::container[@n="materials"]/container/blockdata[@n="material"]}
            } else {
                set xp [format_xpath {../blockdata[@n=%s]} [$domNode @n]]
            }
            set nds [$domNode selectNodes $xp]
            if { [llength $nds] > 1 || ([llength $nds]==1 && $options(-can_delete_last_material)) } {
                set domNode [$tree item element cget $id 0 e_text -data]
                set name [get_domnode_attribute $domNode name]
                if { $name eq "" } {
                    set name [get_domnode_attribute $domNode pn]
                }
                if { $confirm } {
                    set msg [= "Are you sure to delete '%s'?" $name]
                    if {[$domNode @mod_msg ""] != ""} {
                        set msg [get_domnode_attribute $domNode mod_msg]
                    }                        
                    set w [dialogwin_snit $win._ask -title [= "Delete"] -entrytext $msg]                    
                    set f [$w giveframe]
                    
                    if {[$domNode @mod_msg ""] == ""} {
                        ttk::checkbutton $f.cb1 -text [= "Repeat my answer"] -variable \
                            [$w give_uservar repeat 0]                                            
                        grid $f.cb1 -sticky w -pady 2 -padx "5 60"
                    } else {
                        $w set_uservar_value repeat 0
                    }
                    
                    bind [winfo toplevel $f] <Return> [list $w invokeok]
                    set action [$w createwindow]
                    if { [$w give_uservar_value repeat] } { set confirm 0 }
                    destroy $w
                    if { $action < 1 } { return }
                }                
                gid_groups_conds::deleteN $domNode               
                $tree item delete $id             
                $self actualize           
            } elseif { [$domNode @state normal] ne "disabled" && 
                [$domNode @sequence_type any] eq "non_void_disabled" } {
                if { [$domNode selectNodes {.//group}] ne "" } {
                    if { $confirm } {
                        tk_messageBox -message [= "Block cannot be disabled. It contains assigned groups"]
                    }
                    return
                }
                gid_groups_conds::setAttributesN $domNode [list state disabled]
                $self actualize $id
            } elseif { [$domNode @active 1] == 1 } {
                # [$domNode @sequence_type any] eq "non_void_deactivated"
                if { [$domNode selectNodes {.//group}] ne "" } {
                    if { $confirm } {
                        snit_messageBox -message [= "Block cannot be deactivated. It contains assigned groups"]
                    }
                    return
                }
                gid_groups_conds::setAttributesN $domNode [list active 0]
                $self actualize $id
            }
        }
        if { $update_data ne "" } {
            foreach "proc domNode item" $update_data break
            gid_groups_conds::eval_proc -tree $tree -boundary_conds \
                $self -item $item $proc $domNode
        }
    }
    method enable_activate_blockdata {} {
        set ids [$tree selection get]
        if { [llength $ids] == 0 } {
            snit_messageBox -parent $win \
                -message [= "It is necessary to select a blockdata item"]
            return
        }
        foreach id $ids {
            set domNode [$tree item element cget $id 0 e_text -data]
            gid_groups_conds::uncompress_subtree $domNode
            
            if { [$domNode @state normal] eq "disabled" && 
                [$domNode @sequence_type any] eq "non_void_disabled" } {
                gid_groups_conds::setAttributesN $domNode [list state normal]
            } elseif { [$domNode @active 1] == 0 } {
                #[$domNode @sequence_type any] eq "non_void_deactivated"
                gid_groups_conds::setAttributesN $domNode [list active 1]
            }
            $self actualize $id
        }
    }
    method reorder_blockdata_or_group { where args } {
        set needs_actualizeList ""
        set ids [$tree selection get]
        if { $where eq "down" || $where eq "last" } {
            set ids [lreverse $ids]
        }
        set update_data ""
        set idx 0
        foreach id $ids {
            set domNode [$tree item element cget $id 0 e_text -data]
            gid_groups_conds::uncompress_subtree $domNode
            set isgood 0
            if { [$domNode nodeName] eq "blockdata" && [$domNode @sequence 0] } {
                set isgood 1
            }
            if { [$domNode nodeName] eq "group" && [[$domNode parentNode] nodeName] eq "groupList" } {
                set domNode [$domNode parentNode]
            }
            if { [$domNode nodeName] eq "group" && [[$domNode parentNode] nodeName] eq "condition" } {
                set isgood 1
            }
            if { [$domNode nodeName] eq "groupList" && [[$domNode parentNode] nodeName] eq "condition" } {
                set isgood 1
            }
            if { !$isgood } { continue }
            
            set before_update_proc [$domNode @before_update_proc ""]
            if { $before_update_proc ne "" } {
                set ret [catch { gid_groups_conds::eval_proc -tree $tree -boundary_conds \
                            $self -item $id $before_update_proc $domNode } errstring]
                if { $ret } {
                    if { $errstring ne "" } {
                        tk_messageBox -message $errstring
                    }
                    return
                }
            }
            
            set name [$domNode nodeName]
            switch $where {
                up {
                    set xp "preceding-sibling::$name\[1\]"
                    set domNodeRef [$domNode selectNodes $xp]
                    set before_after before
                }
                first {
                    set xp "preceding-sibling::$name\[position()=last()\]"
                    set domNodeRef [$domNode selectNodes $xp]
                    set before_after before
                }
                down {
                    set xp "following-sibling::$name\[1\]"
                    set domNodeRef [$domNode selectNodes $xp]
                    set before_after after
                }
                last {
                    set xp "following-sibling::$name\[position()=last()\]"
                    set domNodeRef [$domNode selectNodes $xp]
                    set before_after after
                }
                moveto {
                    set c [lindex $args 0]
                    set xp [format_xpath {ancestor::container[@n="materials"]/container[@pn=%s]} $c]
                    set domNodeRef [$domNode selectNodes $xp]
                    set before_after child
                }
            }
            if { $domNodeRef ne "" } {
                lappend needs_actualizeList [$domNode parentNode]
                if { $before_after eq "child" } {
                    lappend needs_actualizeList $domNodeRef
                }
                gid_groups_conds::moveNodeN $domNodeRef $domNode $before_after
                
                if { [$domNode @update_proc ""] ne "" } {
                    lappend update_data [$domNode @update_proc] \
                        $domNode $id
                }
            }
            incr idx
        }
        if { $idx == 0 } {
            snit_messageBox -parent $win \
                -message [= "It is necessary to select a blockdata or group item"]
            return
        }
        $self actualize_domNode {*}$needs_actualizeList
        if { $update_data ne "" } {
            foreach "proc domNode item" $update_data {
                gid_groups_conds::eval_proc -tree $tree -boundary_conds \
                    $self -item $item $proc $domNode
            }
        }
    }
    # these funcs are only valid for simple values with or without units (no functions)
    method set_edited_values { d } {
        foreach n [array names entryvalues] {
            if { [string match "_*" $n] } { continue }
            if { [dict exists $d $n] } {
                set entryvalues($n) [dict get $d $n v]
                if { [info exists entryvalues_units($n)] } {
                    set entryvalues_units($n) [gid_groups_conds::units_to_nice_units \
                            [dict get $d $n units]]
                }
            }
        }
    }
    method give_edited_values {} {
        set d ""
        foreach n [array names entryvalues] {
            if { [string match "_*" $n] } { continue }
            dict set d $n v $entryvalues($n)
            if { [info exists entryvalues_units($n)] } {
                set newunits [gid_groups_conds::nice_units_to_units $entryvalues_units($n)]
                dict set d $n units $newunits
            }
        }
        return $d
    }
    method check_edited_value { args } {
        
        set optional {
            { -entry_units widget "" }
            { -check_value boolean 1 }
            { -value number "" }
            { -display_message boolean 1 }
            { -string_is string_is ""}
        }
        set compulsory "parent domNode magnitude unit"
        parse_args $optional $compulsory $args
        
        if { ![gid_groups_conds::is_unit_correct $magnitude $unit] } {
            if { $display_message } {
                set pn [$domNode @pn]
                snit_messageBox -message [= "Unit for property '%s' is not correct" $pn] \
                    -parent $parent
            }
            if { $entry_units ne "" } { $entry_units traverse_in units }
            return -1
        }
        if { $check_value && $string_is != "" } {
            set check_value 0
            if {$string_is eq "integer_or_void"} {
                if {$value ne ""} { set string_is "integer" } 
            } 
            if {$string_is eq "double_or_void"} {
                if {$value ne ""} { set string_is "double" } 
            } 
            if { $string_is eq "double" || $string_is eq "integer" } {
                if { ![string is $string_is -strict $value] } {
                    if { [info command ::draw_post::resolve_parametric] ne "" } {
                        set err [catch { draw_post::resolve_parametric -parent $win -ask_user 1 \
                                    -- $value } ret]
                        if { !$err } {
                            set value $ret
                        }
                    }
                }
                if { ![string is $string_is -strict $value] } {
                    if { $display_message } {
                        set pn [$domNode @pn]
                        snit_messageBox -message [= "Value for property '%s' is not correct" $pn] \
                            -parent $parent
                    }
                    return -1
                }              
            }
            if { $string_is eq "%" } {
                set isok 1
                if { ![string is double -strict $value] } {
                    set isok 0
                }
                if { !$isok || $value > 100 } {
                    if { $display_message } {
                        set pn [$domNode @pn]
                        snit_messageBox -message [= "Percentage for property '%s' is not correct" $pn] \
                            -parent $parent
                    }
                    return -1         
                }           
            }
            if { $string_is eq "list_of_double" } {               
                foreach val $value {
                    if { ![string is double -strict $val] } {
                        if { $display_message } {
                            set pn [$domNode @pn]
                            snit_messageBox -message [= "Value for property '%s' is not correct" $pn] \
                                -parent $parent
                        }
                        return -1
                    }
                }                
            }
        }    
        set err 0    
        if { $check_value && ![string is double -strict $value] } {
            if { [info command ::draw_post::resolve_parametric] ne "" } {
                set err [catch { draw_post::resolve_parametric -parent $win -ask_user 1 \
                            -unit_magnitude $magnitude -units $unit -- $value } ret]
            } else {
                set err 1
            }
        }
        if { $err } {
            if { $display_message } {
                set pn [$domNode @pn]
                snit_messageBox -message [= "Value for property '%s' is not correct" $pn] \
                    -parent $parent
            }
            if { $entry_units ne "" } { $entry_units traverse_in entry }
            return -1
        }
    }
    method check_edited_values { args } {
        
        set optional {
            { -display_message boolean 1 }
        }
        set compulsory ""
        parse_args $optional $compulsory $args
        
        if { [info exists entryvalues(_selection_entity_)] } {
            set ov $entryvalues(_selection_entity_)
        } else {
            set ov ""
        }
        foreach n [array names entryvalues_units] {
            if { $entryvalues_state($n) eq "disabled" } { continue }
            if { [info exists entryvalues_function($n)] && $entryvalues_function($n) ne "" } {
                continue
            }
            set found 0
            foreach domNode [array names entry_node_values] {
                if { [string first "," $domNode] != -1 } { continue }
                if { [$domNode @n ""] eq $n } { 
                    set found 1
                    break
                }
            }            
            if {!$found} { continue }            
            if { [$domNode @n ""] ne $n } { continue }
            set unit_definition [get_domnode_attribute $domNode unit_definition]
            set unit_mesh_definition [get_domnode_attribute $domNode unit_mesh_definition]
            set disable_on_void [get_domnode_attribute $domNode disable_on_void 0]
            
            if { $unit_definition ne "" } {
                set magnitude $unit_definition
                set check_value 0
            } elseif { $unit_mesh_definition == 1 } {
                set magnitude L
                set check_value 0
            } else {
                set magnitude [gid_groups_conds::give_magnitude [$domNode @unit_magnitude] $ov]
                set check_value 1
            }
            set newunits [gid_groups_conds::nice_units_to_units $entryvalues_units($n)]
            set entry_units [lindex $entry_node_values($domNode) 1 1]
            
            if { $disable_on_void && [string trim $entryvalues($n)] eq "" } {
                set ret 0
            } else {
                set ret [$self check_edited_value -entry_units $entry_units -check_value $check_value \
                        -value $entryvalues($n) -display_message $display_message -string_is [$domNode @string_is ""] \
                        $entry_units $domNode $magnitude $newunits ]
            }
            if { $ret == -1 } { return -1 }
        }
        foreach domNode [array names entry_node_values] {
            if { [string first "," $domNode] != -1 } { continue }
            set n [$domNode @n ""]
            if { $n in [array names entryvalues_units] } { continue }
            if {[$domNode @string_is ""] != ""} {
                set check_value 1
                lassign [list "" "" ""] entry_units newunits magnitude
                set ret [$self check_edited_value -entry_units $entry_units -check_value $check_value \
                        -value $entryvalues($n) -display_message $display_message -string_is [$domNode @string_is ""] \
                        $entry_units $domNode $magnitude $newunits]
                if { $ret == -1 } { return -1 }
            }
        }
        return 0
    }
    method apply_edited_values { id } {
        set domNode [$tree item element cget $id 0 e_text -data]
        
        set actualize_nodesList ""
        
        foreach i $entryvalues_mat_changes_stack {
            set args [lassign $i op]
            switch $op {
                add {
                    lassign $args n matname
                    set xp0 {ancestor::container[@n="materials"]//}
                    append xp0 {blockdata[@n="material" and @name=%s]/}
                    append xp0 {container[@n=%s]}
                    set xp [format_xpath $xp0 $matname $n]
                    set containerNode_from [$domNode selectNodes $xp]
                    $domNode appendChild [$containerNode_from cloneNode -deep]
                    lappend actualize_nodesList $domNode
                }
                delete {
                    set n [lindex $args 0]
                    set xp [format_xpath {.//container[@n=%s]} $n]
                    set node [$domNode selectNodes $xp]
                    lappend actualize_nodesList [$node parentNode]
                    $node delete
                }
            }
        }
        set nal [$self apply_edited_values_do $domNode $id]
        lappend actualize_nodesList {*}$nal
        return [lsort -unique $actualize_nodesList]
    }
    method apply_edited_values_do { domNode id } {
        
        set actualize_nodesList ""
        
        foreach valueNode [$domNode childNodes] {
            set d ""
            if { $id ne "" } {
                if {$id in [$tree item range 0 end]} {
                    foreach id_son [$tree item children $id] {
                        set d [$tree item element cget $id_son 0 e_text -data]
                        if { $valueNode eq $d } { break }
                    }
                } else {                   
                    continue
                }
            }
            if { $valueNode ne $d } { set id_son "" }
            switch -- [$valueNode nodeName] {
                value {
                    set pn [get_domnode_attribute $valueNode pn]
                    set n [$valueNode @n]
                    set unit_magnitude [get_domnode_attribute $valueNode unit_magnitude]
                    set units_system_definition [get_domnode_attribute $valueNode units_system_definition]
                    set unit_definition [get_domnode_attribute $valueNode unit_definition]
                    set unit_mesh_definition [get_domnode_attribute $valueNode unit_mesh_definition]
                    if { ![info exists entryvalues_state($n)] } {
                        set state normal
                    } else {
                        set state $entryvalues_state($n)
                    }
                    
                    set dict ""
                    foreach "ni vi" [split [get_domnode_attribute $valueNode dict] ,] {
                        dict set dict $ni [= $vi]
                    }
                    if { $units_system_definition == 1 } {
                        foreach i [gid_groups_conds::give_units_system_list] {
                            dict set dict [lindex $i 0] [= [lindex $i 1]]
                        }
                    }
                    if { $id_son ne "" } {
                        set txt_old [$tree item text $id_son 0]
                    } else {
                        set txt_old [gid_groups_conds::give_printable_value $valueNode]
                    } 
                    set newvalue ""                  
                    if { [info exists entryvalues($n)] } { 
                        set newvalue $entryvalues($n)         
                    }          
                    set newvalueT [dict_getd $dict $newvalue $newvalue]
                    
                    if { $unit_definition ne "" || $unit_mesh_definition == 1 } {
                        set unitN $entryvalues_units($n)
                        set newunits [gid_groups_conds::nice_units_to_units $unitN]
                        set txt "$pn: $unitN"
                        
                        if { $unit_definition ne "" } {
                            set magnitude $unit_definition
                        } else {
                            set magnitude L
                        }
                        gid_groups_conds::add_to_units_preferences $magnitude $newunits
                    } elseif { $unit_magnitude ne "" } {
                        set unitN $entryvalues_units($n)
                        set newunits [gid_groups_conds::nice_units_to_units $unitN]
                        set txt "$pn: $newvalueT $unitN"
                        gid_groups_conds::add_to_units_preferences $unit_magnitude $newunits
                    } else {
                        set txt "$pn: $newvalueT"
                    }
                    if { $txt_old ne $txt || [$valueNode @function ""] ne "" } {
                        if { $unit_definition ne "" } {
                            gid_groups_conds::set_active_unit $unit_definition $newunits
                        } elseif { $unit_mesh_definition ne "" } {
                            gid_groups_conds::set_mesh_unit $newunits
                        } elseif { $units_system_definition == 1 } {
                            gid_groups_conds::set_active_units_system $newvalue
                        } else {
                            set attList [list v $newvalue]
                            if { $unit_magnitude ne "" } {
                                lappend attList units $newunits
                                set resolve_parametric 1
                            } elseif { [$valueNode @string_is ""] in "integer_or_void double_or_void integer double" } {
                                set resolve_parametric 1
                            } else {
                                set resolve_parametric 0
                            }
                            gid_groups_conds::setAttributesN -resolve_parametric \
                                $resolve_parametric $valueNode $attList
                        }
                        if { [$valueNode @function ""] ne "" } {
                            set funcNode [$valueNode selectNodes function]
                            if { $funcNode ne "" } {
                                gid_groups_conds::deleteN $funcNode
                            }
                            if { $entryvalues_function($n) ne "" } {
                                gid_groups_conds::copyNode_recursive \
                                    $entryvalues_function($n) $valueNode
                                set entryvalues_function($n) ""
                            }
                            lappend actualize_nodesList $valueNode
                        }
                        if { $state ne "hidden" } {
                            if { $id_son ne "" } {
                                $tree item text $id_son 0 $txt
                            }
                            if { [$valueNode @actualize_tree 0] } {
                                lappend actualize_nodesList 0
                            }
                        }
                        set hc [gid_groups_conds::check_node_dependencies $valueNode]
                        lappend actualize_nodesList $valueNode {*}$hc
                        
                        #                         foreach node [$valueNode selectNodes .//value] {
                            #                             set hc [gid_groups_conds::check_node_dependencies $node]
                            #                             if { $hc && $hc > $needs_actualize && $state ne "hidden" } {
                                #                                 set needs_actualize $hc
                                #                             }
                            #                         }
                    }
                }
                container {
                    set xp {value|container/value}
                    if { [llength [$valueNode selectNodes $xp]] == 0 } { continue }
                    set state [get_domnode_attribute $valueNode state normal]
                    if { $state eq "hidden" } { continue }
                    
                    set nal [$self apply_edited_values_do $valueNode $id_son]
                    lappend actualize_nodesList {*}$nal
                }
            }
        }
        return [lsort -unique $actualize_nodesList]
    }
    method copy_cond_values { domNode_from domNode_to id_to parent_state { level 0 } } {
        
        set actualize_nodesList ""
        set valueNodes [$domNode_from selectNodes value|container|edit_command]
        if { ![llength $valueNodes] && [$domNode_from @isvalue 0] } {
            lappend valueNodes $domNode_from
        }
        foreach valueNode $valueNodes {
            set state [get_domnode_attribute $valueNode state normal]
            if { $parent_state eq "hidden" } { set state hidden }
            switch [$valueNode nodeName] {
                edit_command {
                    set attList ""
                    foreach i [$valueNode selectNodes @*] {
                        lappend attList {*}$i
                    }
                    set v [gid_groups_conds::addN $domNode_to edit_command $attList]
                }
                value - condition {
                    set n [$valueNode @n]
                    set state $entryvalues_state($n)
                    if { $parent_state eq "hidden" } { set state hidden }
                    set pn [get_domnode_attribute $valueNode pn]
                    set unit_magnitude [get_domnode_attribute $valueNode unit_magnitude]
                    set units_system_definition [get_domnode_attribute $valueNode units_system_definition]
                    set unit_definition [get_domnode_attribute $valueNode unit_definition]
                    set unit_mesh_definition [get_domnode_attribute $valueNode unit_mesh_definition]
                    
                    set unitN ""
                    set dict ""
                    foreach "ni vi" [split [get_domnode_attribute $valueNode dict] ,] {
                        dict set dict $ni [= $vi]
                    }
                    if { $units_system_definition == 1 } {
                        foreach i [gid_groups_conds::give_units_system_list] {
                            dict set dict [lindex $i 0] [= [lindex $i 1]]
                        }
                    }
                    set newvalue ""                  
                    if { [info exists entryvalues($n)] } { 
                        set newvalue $entryvalues($n)         
                    }   
                    if { $unit_magnitude ne "" || $unit_definition ne "" || $unit_mesh_definition == 1 } {
                        set unitN $entryvalues_units($n)
                        set newunits [gid_groups_conds::nice_units_to_units $unitN]
                        
                        if { $unit_definition ne "" } {
                            set magnitude $unit_definition
                        } elseif { $unit_mesh_definition == 1 } {
                            set magnitude L
                        } else {
                            set magnitude $unit_magnitude
                        }
                        gid_groups_conds::add_to_units_preferences $magnitude $newunits
                    }
                    
                    set newvalueT [dict_getd $dict $newvalue $newvalue]
                    
                    set attList ""
                    foreach i [$valueNode selectNodes @*] {
                        lappend attList {*}$i
                    }
                    foreach i [list v units] {
                        dict unset attList $i
                    }
                    set resolve_parametric 0
                    
                    if { $unit_definition ne "" } {
                        gid_groups_conds::set_active_unit $unit_definition $newunits
                    } elseif { $unit_mesh_definition ne "" } {
                        gid_groups_conds::set_mesh_unit $newunits
                    } elseif { $units_system_definition == 1 } {
                        gid_groups_conds::set_active_units_system $newvalue
                    } else {
                        dict set attList v $newvalue
                        if { [llength [split $unit_magnitude ,]] > 1 } {
                            set ov $entryvalues(_selection_entity_)
                            dict set attList unit_magnitude [dict get [split $unit_magnitude ,] $ov]
                        }
                        if { $unit_magnitude ne "" } {
                            dict set attList units $newunits
                            set resolve_parametric 1
                        } elseif { [$valueNode @string_is ""] in "integer_or_void double_or_void integer double" } {
                            set resolve_parametric 1
                        }
                    }
                    foreach i [list ov isvalue] {
                        dict unset attList $i
                    }
                    set v [gid_groups_conds::addN -resolve_parametric $resolve_parametric \
                            $domNode_to value $attList]
                    
                    set actualize_func 0
                    if { [$valueNode @function ""] ne "" } {
                        if { $entryvalues_function($n) ne "" } {
                            set funcNode $entryvalues_function($n)
                            gid_groups_conds::copyNode_recursive \
                                $funcNode $v
                            
                            set nameList ""
                            foreach fvarNode [$funcNode selectNodes functionVariable] {
                                lappend nameList [$fvarNode @pn]
                            }
                            set newvalueT "f([join $nameList ,])..."
                            set actualize_func 1
                        }
                    }
                    set txt "$pn: $newvalueT"
                    if { $unitN ne "" } { append txt " $unitN" }
                    
                    if { $state ne "hidden" } {
                        set id [$self add $id_to $txt \
                                [gid_groups_conds::giveimage editclear-16] \
                                value $v]
                        $self item collapse $id
                        if { [$valueNode @actualize_tree 0] } {
                            lappend actualize_nodesList 0
                        } elseif { $actualize_func } {
                            lappend actualize_nodesList $valueNode
                        }
                    }
                }
                container {
                    set attList ""
                    foreach i [$valueNode selectNodes @*] { eval lappend attList $i }
                    set c [gid_groups_conds::addN $domNode_to container $attList]
                    
                    if { $state ne "hidden" } {
                        set image [$self _give_icon $valueNode]
                        set id_c [$self add $id_to [get_domnode_attribute \
                                    $valueNode pn] \
                                $image container $c]
                        $self item collapse $id_c
                    } else {
                        set id_c ""
                    }
                    set nal [$self copy_cond_values $valueNode $c $id_c $state \
                            [expr {$level+1}]]
                    
                    lappend actualize_nodesList {*}$nal
                }
            }
        }
        if { $level == 0 } {
            foreach valueNode [$domNode_to selectNodes ".//value"] {
                set hc [gid_groups_conds::check_node_dependencies $valueNode]
                lappend actualize_nodesList {*}$hc
            }
        }
        return [lsort -unique $actualize_nodesList]
    }
    proc nice_xpath { domNode } {
        set doc [$domNode ownerDocument]
        set root [$doc documentElement]
        
        set nodeList ""
        while { $domNode ne $root } {
            set parentNode [$domNode parentNode]
            if { [$domNode @n ""] ne "" } {
                set xp [$domNode nodeName]\[@n='[$domNode @n]'\]
            } else {
                set xp [$domNode nodeName]
            }
            if { [llength [$parentNode selectNodes $xp]] > 1 } {
                set ipos [lsearch -exact [$parentNode selectNodes $xp] $domNode]
                append xp \[[expr {$ipos+1}]\]
            }
            set nodeList [linsert $nodeList 0 $xp]
            set domNode $parentNode
        }
        return [join $nodeList /]
    }
    method edit_value_or_blockdata {} {
        set ids [$tree selection get]
        if { [llength $ids] != 1 } {
            return 0
        }
        set id [lindex $ids 0]
        
        set domNode [$tree item element cget $id 0 e_text -data]
        
        $self cancel_cond_or_blockdata_edit
        
        if { [$domNode nodeName] eq "blockdata" } {
            if { [$domNode @editable_name ""] eq "1" || \
                [$domNode @editable_name ""] eq "unique" } {
                set values ""
                if { [$domNode @n] eq "material" } {
                    set xp {ancestor::container[@n="materials"]/blockdata[@n="material"]|}
                    append xp {ancestor::container[@n="materials"]/container/blockdata[@n="material"]}
                } else {
                    set xp [format_xpath {../blockdata[@n=%s]} [$domNode @n]]
                }
                foreach node [$domNode selectNodes $xp] {
                    lappend values [get_domnode_attribute $node name]
                }
                ::TreeCtrl::ComboOpen $tree $id 0 e_text 1 \
                    $values [get_domnode_attribute $domNode name]               
                return 1
            } else {
                set edit_values 0
                foreach ichild [$domNode childNodes] {
                    if { [$ichild nodeType] ne "condition"} {
                        set edit_values 1; break
                    }
                }
                if { $edit_values } { return 1 } else { return 0 }                
            } 
        } elseif { [$domNode nodeName] eq "value" } {
            if { [get_domnode_attribute $domNode state] eq "disabled" } { return 0 }
            if { [$domNode selectNodes function] ne "" } { return }
            set function [get_domnode_attribute $domNode function ""]
            if { $function ne "" && [lsearch [split $function ,] "scalar"] == -1 } {
                return
            }
            set editable [$domNode @editable 0]
            set value [get_domnode_attribute $domNode v]
            array unset dict
            foreach "ni vi" [split [get_domnode_attribute $domNode dict] ,] {
                set dict($ni) [= $vi]
            }
            if { [info exists dict($value)] } {
                set value $dict($value)
            }
            set values [split [get_domnode_attribute $domNode values] ,]
            set values_tree [split [get_domnode_attribute $domNode values_tree] ,]
            set values_check [split [get_domnode_attribute $domNode values_check] ,]
            
            if { [get_domnode_attribute $domNode units_system_definition] == 1 } {
                foreach i [gid_groups_conds::give_units_system_list] {
                    lappend values [lindex $i 0]
                    set dict([lindex $i 0]) [= [lindex $i 1]]
                }
                lassign [gid_groups_conds::give_active_units_system] - value
            } 
            
            if { [llength $values] } {
                for { set i 0 } { $i < [llength $values] } { incr i } {
                    if { [info exists dict([lindex $values $i])] } {
                        lset values $i $dict([lindex $values $i])
                    }
                }
            } elseif { [llength $values_tree] } {
                set idx 0
                foreach i $values_tree {
                    lset values_tree $idx 3 [gid_groups_conds::giveimage [lindex $i 3]]
                    incr idx
                }
            } else {
                set editable 1
                lappend values $value
            }
            set unit_magnitude [get_domnode_attribute $domNode unit_magnitude]
            set unit_definition [get_domnode_attribute $domNode unit_definition]
            set unit_mesh_definition [get_domnode_attribute $domNode unit_mesh_definition]
            if { $unit_magnitude ne "" || $unit_definition ne "" || $unit_mesh_definition == 1 } {
                destroy $tree.entry_units
                gid_groups_conds::entry_units $tree.entry_units -value_node $domNode
                TreeCtrl::AnyWidgetOpen $tree $id 0 e_text $tree.entry_units
            } elseif { [llength $values_tree] } {
                ::TreeCtrl::ComboTreeOpen $tree $id 0 e_text $editable \
                    $values_tree $value
            } elseif { [llength $values_check] } {
                destroy $tree.entry_tree_check
                
                if { [info command cu::menubutton_check_listbox1] eq "" } {
                    snit::widgetadaptor cu::menubutton_check_listbox1 {
                        delegate method * to hull
                        delegate option * to hull
                        constructor args {
                            installhull [cu::menubutton_check_listbox $self]
                            $self configurelist $args
                        }
                        method get_text {} { return [$self get_selected_comma_list] }
                    }
                }
                cu::menubutton_check_listbox1 $tree.entry_tree_check -values $values_check \
                    -label_prefix [= "Number of %s" [$domNode @pn]]
                $tree.entry_tree_check set_selected_comma_list $value
                TreeCtrl::AnyWidgetOpen $tree $id 0 e_text $tree.entry_tree_check
            } else {
                ::TreeCtrl::ComboOpen $tree $id 0 e_text $editable \
                    $values $value
            }
            return 1
        }
    }
    method edit_value_or_blockdata_ok { id newvalue } {      
        if { [lsearch [$tree range 1 last] $id] == "-1" } {
            set id [lindex [$tree selection get] 0]           
            if { $id eq "" } { return }
        }
        set domNode [$tree item element cget $id 0 e_text -data]
        if { [$domNode nodeName] eq "container" } { return }
        if { [$domNode nodeName] eq "value" } {
            set text ""
            if {$id != 0} {               
                set text [$tree item element cget $id 0 e_text -text]
            }
            regexp {^([^:]+):} $text {} before
            
            if { ![info exists before] } { set before "" }
            
            set unit_magnitude [get_domnode_attribute $domNode unit_magnitude]
            set unit_definition [get_domnode_attribute $domNode unit_definition]
            set units_system_definition [get_domnode_attribute $domNode units_system_definition]
            set unit_mesh_definition [get_domnode_attribute $domNode unit_mesh_definition]
            set disable_on_void [get_domnode_attribute $domNode disable_on_void 0]
            set string_is [get_domnode_attribute $domNode string_is]
            
            if { $unit_definition ne "" || $unit_mesh_definition == 1 } {
                set unit [gid_groups_conds::nice_units_to_units $newvalue]
                if { $unit_definition ne "" } {
                    set magnitude $unit_definition
                } else {
                    set magnitude L  
                }
                set ret [$self check_edited_value -check_value 0 -string_is $string_is $tree $domNode $magnitude $unit]
                if { $ret == -1 } { return }
                set txt "$before: $newvalue"
                
                gid_groups_conds::add_to_units_preferences $magnitude $unit
            } elseif { $unit_magnitude ne "" } {
                lassign $newvalue newvalue unitN
                set unit [gid_groups_conds::nice_units_to_units $unitN]
                
                if { $disable_on_void && [string trim $newvalue] eq "" } {
                    set ret 0
                } else {
                    set ret [$self check_edited_value -check_value 1 -value $newvalue -string_is $string_is \
                            $tree $domNode $unit_magnitude $unit]
                }
                if { $ret == -1 } { return }
                set txt "$before: $newvalue $unitN"
                gid_groups_conds::add_to_units_preferences $unit_magnitude $unit
            } elseif { $string_is ne ""} {
                set unit ""
                set ret [$self check_edited_value -check_value 1 \
                        -value $newvalue -string_is $string_is \
                        $tree $domNode $unit_magnitude $unit]
                if { $ret == -1 } { return }
                set txt "$before: $newvalue"                
            } else {
                set txt "$before: $newvalue"
            }
            $tree item element configure $id 0 e_text -text $txt
            
            set dict ""
            foreach "ni vi" [split [get_domnode_attribute $domNode dict] ,] {
                dict set dict $ni [= $vi]
            }
            
            if { $units_system_definition == 1 } {
                foreach i [gid_groups_conds::give_units_system_list] {
                    dict set dict [lindex $i 0] [= [lindex $i 1]]
                }
            }
            set dict_inv [dict_inverse $dict]
            if { [dict exists $dict_inv $newvalue] } {
                set newvalue [dict get $dict_inv $newvalue]
            }
            
            set actualize_nodesList ""
            
            if { $unit_definition ne "" } {
                gid_groups_conds::set_active_unit $unit_definition $unit
            } elseif { $unit_mesh_definition ne "" } {
                gid_groups_conds::set_mesh_unit $unit
            } elseif { $units_system_definition == 1 } {
                gid_groups_conds::set_active_units_system $newvalue
            } else {
                set attList [list v $newvalue]
                if { $unit_magnitude ne "" } {
                    lappend attList units $unit
                    set resolve_parametric 1
                } elseif { [$domNode @string_is ""] in "integer_or_void double_or_void integer double" } {
                    set resolve_parametric 1
                } else {
                    set resolve_parametric 0
                }
                gid_groups_conds::setAttributesN -resolve_parametric $resolve_parametric \
                    $domNode $attList
                lappend actualize_nodesList $domNode
            }
            if { [$domNode @actualize_tree 0] } {
                lappend actualize_nodesList 0
            }
            set hc [gid_groups_conds::check_node_dependencies $domNode]
            lappend actualize_nodesList {*}$hc
            
            set p $id
            while { $p != 0 && ![$tree item state get $p blockdata] } {
                set p [$tree item parent $p]
            }
            if { $p != 0 } {
                set domNodeBD [$tree item element cget $p 0 e_text -data]
                
                if { [$domNodeBD @state normal] eq "disabled" &&
                    [$domNodeBD @sequence_type any] eq "non_void_disabled" } {
                    gid_groups_conds::setAttributesN $domNodeBD \
                        [list state normal]
                    lappend actualize_nodesList $domNodeBD
                } elseif { [$domNodeBD @active 1] == 0 } {
                    #[$domNodeBD @sequence_type any] eq "non_void_deactivated"
                    gid_groups_conds::setAttributesN $domNodeBD \
                        [list active 1]
                    lappend actualize_nodesList $domNodeBD
                }
            }
            if { [$domNode @state normal] eq "disabled_on_void" && \
                [$domNode @v ""] eq "" } {
                gid_groups_conds::setAttributesN $domNode [list state disabled]
            }
            $self actualize_domNode {*}$actualize_nodesList
            
            if { [$domNode selectNodes {ancestor::condition}] ne "" } {
                gid_groups_conds::check_draw_nodes_state
            }
            if { [set proc [$domNode @update_proc ""]] ne "" } {
                gid_groups_conds::eval_proc -tree $tree -boundary_conds \
                    $self -item $id $proc $domNode
            }
        } elseif { [$domNode nodeName] eq "blockdata" } {
            if { $newvalue eq [get_domnode_attribute $domNode name] } { return }
            set values ""
            if { [$domNode @n] eq "material" } {
                set xp {ancestor::container[@n="materials"]/blockdata[@n="material"]|}
                append xp {ancestor::container[@n="materials"]/container/blockdata[@n="material"]}
            } else {
                set xp [format_xpath {../blockdata[@n=%s]} [$domNode @n]]
            }
            foreach node [$domNode selectNodes $xp] {
                lappend values [get_domnode_attribute $node name]
            }
            if { [lsearch -exact $values $newvalue] != -1 } {
                snit_messageBox -type ok  -parent $win \
                    -message [= "Identifier '%s' already exists" $newvalue]
                return
            }
            set before_update_proc [$domNode @before_update_proc ""]
            if { $before_update_proc ne "" } {
                set ret [catch { gid_groups_conds::eval_proc -tree $tree -boundary_conds \
                            $self -item $id $before_update_proc $domNode } errstring]
                if { $ret } {
                    if { $errstring ne "" } {
                        tk_messageBox -message $errstring
                    }
                    return
                }
            } 
            $tree item element configure $id 0 e_text -text "$newvalue"
            gid_groups_conds::setAttributesN $domNode [list name $newvalue]
            
            if { [set proc [$domNode @update_proc ""]] ne "" } {
                gid_groups_conds::eval_proc -tree $tree -boundary_conds \
                    $self -item $id $proc $domNode
            }
        } else { error "error in edit_value_or_blockdata_ok" }
    }
    proc valuesList { groupNode } {
        set list ""
        foreach valueNode [$groupNode selectNodes .//value] {
            lappend list [list [$valueNode @n] [$valueNode @v] \
                    [$valueNode @units ""]]
        }
        if { ![llength $list] } {
            set cndNode [$groupNode parentNode]
            lappend list [list [$cndNode @n] [get_domnode_attribute \
                        $cndNode @pn] ""]
        }
        return $list
    }
    proc valuesNamesList { conditionNode } {
        set list ""
        foreach valueNode [$conditionNode selectNodes value|container/value] {
            lappend list [$valueNode @pn]
        }
        if { ![llength $list] } {
            set list [list [get_domnode_attribute $conditionNode pn]]
        }
        return $list
    }
    method list_entities {} {
        set domNodes ""
        foreach id [$tree selection get] {
            set domNode [$tree item element cget $id 0 e_text -data]
            gid_groups_conds::uncompress_subtree $domNode
            lappend domNodes $domNode
        }
        gid_groups_conds::list_entities_cnd $win $domNodes
    }
    # what: groupnames, values, symbols, symbols_all
    method draw_groups_do { args } {
        
        set optional {
            { -domNodes domNodes "" }
        }
        set compulsory "what"
        parse_args $optional $compulsory $args
        
        set b $win.finish.b
        if { [winfo exists $b] } {
            $b invoke
        } else {
            $self end_draw_groups_do "" "" finish
        }
        if { $domNodes eq "" } {
            foreach id [$tree selection get] {
                lappend domNodes [$tree item element cget $id 0 e_text -data]
            }
        }
        set groups [gid_groups_conds::draw_nodes -domNodes $domNodes $what]
        
        set t $win
        destroy $t.finish
        set f [frame $t.finish -bd 1 -relief groove]
        set scr "[list destroy $f]; raise $t"
        set b [ttk::button $f.b -text [= Finish] -width 10 -command \
                [mymethod end_draw_groups_do $scr $f.b finish]]
        if { $what eq "symbols" || $what eq "symbols_all" } {
            set b2 [ttk::button $f.b2 -text [= "Keep drawing"] -command \
                    [mymethod end_draw_groups_do $scr $f.b maintain]]
            gid_groups_conds::register_popup_help $b2 \
                [= "Continue drawing the symbols after leaving the drawing function"]
        }
        
        if 0 {
            set text [text $f.t -bd [$f cget -bd] \
                    -tabs "10 [expr {[winfo width $t]-20}] right"]
            
            $text insert end "[= {Drawing groups}]:\n\n"
            $text insert end "\t[= Group]\t[= Values]\n" bold
            $text tag configure bold -font "[font configure [$text cget -font]] -weight bold"
            foreach i $groups {
                $text insert end "\t[lindex $i 0]\t[lindex $i 2]\n"
            }
            $text configure -state disabled
            
            if { ![info exists b2] } {
                grid $b - -padx 2 -pady 10 -sticky nw
            } else {
                grid $b $b2 -padx 2 -pady 10 -sticky nw
            }
            grid $text - -padx 2 -pady 10 -sticky nwes
        } else {
            set columns [list \
                    [list 15 [= "Group"] left text 1] \
                    ]
            if { $what eq "groupnames" } {
                lset columns 0 0 30
            } else {
                set valuesNames ""
                foreach i $groups {
                    set v ""
                    foreach j [lindex $i 2] { lappend v [lindex $j 0] }
                    if { $valuesNames eq "" } {
                        set valuesNames $v
                    } elseif { $valuesNames ne $v } {
                        set valuesNames ""
                        break
                    }
                }
                if { $valuesNames ne "" } {
                    foreach i $valuesNames {
                        lappend columns [list 10 $i left text 1]
                    }
                } else {
                    lappend columns [list 10 [= Values] left text 1]
                    lset columns 1 0 60
                }
            }
            package require fulltktree
            fulltktree $f.tree -width 500 -height 200 \
                -columns $columns -expand 0 -showbuttons 0 \
                -selectmode extended -showlines 0 -indent 0
            foreach i $groups {
                set list [list [lindex $i 0]]
                set l ""
                foreach j [lindex $i 2] {
                    lappend l [lindex $j 1]
                }
                if { $what eq "groupnames" } {
                    #nothing
                } elseif { $valuesNames ne "" } {
                    eval lappend list $l
                } else {
                    lappend list $l
                }
                $f.tree insert end $list
            }
            if { ![info exists b2] } {
                grid $b - -padx 2 -pady 10 -sticky nw
            } else {
                grid $b $b2 -padx 2 -pady 10 -sticky nw
            }
            grid $f.tree - -padx 2 -pady 10 -sticky nwes
        }
        grid columnconfigure $f 1 -weight 1
        grid rowconfigure $f 1 -weight 1
        
        $self _keep_frame_on_top $t $f
        
        #grab $b
        bind $b <Escape> "[list $b invoke]; break"
        bind $b <Return> [list $b invoke]
        focus $b
        if { !$gid_groups_conds::is_draw_post } {
            delayedop changefunc [list $b invoke]
        } else {
            draw_post::manage_buttons_set_graphmode \
                ENDCMD end_cmd [list $b invoke]   
        }
    }
    method end_draw_groups_do { scr b what } {
        if { $b ne "" } {
            if { !$gid_groups_conds::is_draw_post } {
                set found [delayedop -cancel changefunc [list $b invoke]]
                if { $found } { GiD_Process escape escape }
            } else {
                draw_post::end_draw_graphmode   
            }
        }
        eval $scr
        if { $what eq "finish" } {
            gid_groups_conds::end_drawall
        } else {
            gid_groups_conds::end_drawall_maintain
        }
    }
    method _keep_frame_on_top { t f } {
        if { ![winfo exists $f] } {
            bind $t <Configure> ""
            return
        }
        place $f -in $t -x 0 -y 0 -anchor nw -width \
            [winfo width $t] -height [winfo height $t]
        bind $t <Configure> [mymethod _keep_frame_on_top $t $f]
    }
    proc _tree_preferences_give_defaults {} {
        return [dict create background_color auto text_color auto]
    }
    method tree_preferences { args } {
        
        set optional {
            { -to_defaults "" 0 }
        }
        set compulsory ""
        parse_args $optional $compulsory $args
        
        set w $win.tree_prefs
        destroy $w
        dialogwin_snit $w -title [= "Tree preferences"] -callback \
            [mymethod tree_preferences_do] -grab 0 -transient 1 \
            -morebuttons [list [= "Apply"] [= "Defaults"]]
        set f [$w giveframe]
        
        set tree_prefsNode [gid_groups_conds::get_preference_node -local_global global \
                tree_preferences]
        
        set values [_tree_preferences_give_defaults]
        if { !$to_defaults } {
            dict for "n v" $values {
                if { [$tree_prefsNode hasAttribute $n] } {
                    dict set values $n [$tree_prefsNode @$n]
                }
            }
        }
        set xml "<formulae usableviewmode='1' tabstyle='notebook' updatebutton='0'/>"
        set docf [dom parse $xml]
        set root [$docf documentElement]
        dom createNodeCmd text t
        foreach node [list title description container param] {
            dom createNodeCmd element $node
        }
        
        $root appendFromScript {
            title { t [= "Tree preferences"] }
            description {}
            container n "graphical" pns [= "Graphical"] pn [= "graphical"] { \
                param n "background_color" pn [= "Background color"] \
                field_type color \
                value [dict get $values background_color] \
                default_value white \
                disabled_value [list auto white] \
                help [= "Background color of the tree. If deactivated, it is chosen automatically"] \
                {}
                param n "text_color" pn [= "Text color"] \
                field_type color \
                value [dict get $values text_color] \
                default_value darkblue \
                disabled_value [list auto darkblue] \
                help [= "Color of the text of the tree. If deactivated, it is chosen automatically"] \
                {}
            }
        }
        set f [$w giveframe]
        set key $w
        set f_in [formulae::create_windowD -state disabled \
                $f $key $root]
        grid $f_in -sticky nsew
        grid rowconfigure $f 0 -weight 1
        grid columnconfigure $f 0 -weight 1
        
        bind $w <Return> [list $w invokeok]
        set action [$w createwindow]
        $docf delete
    }
    method tree_preferences_do { w } {
        
        if { [$w giveaction] < 1 } {
            destroy $w
            return
        }
        if { [$w giveaction] == 3 } {
            $self tree_preferences -to_defaults
            return
        }
        set key $w
        set d [formulae::give_params_dict $key]
        
        set tree_prefsNode [gid_groups_conds::get_preference_node -local_global global \
                tree_preferences]
        set values [_tree_preferences_give_defaults]
        foreach n [dict keys $values] {
            set v [dict get $d $n]
            $tree_prefsNode setAttribute $n $v
        }
        $self tree_changed_size
        if { [$w giveaction] == 1 } { destroy $w }
    }
    method import_export_materials {} {
        
        set id [lindex [$tree selection get] 0]
        set domNode [$tree item element cget $id 0 e_text -data]
        set cNode [$domNode selectNodes {ancestor-or-self::container[@n="materials"]}]
        set xp [$cNode toXPath]
        gid_groups_conds::import_export_materials -boundary_conds_sel $self \
            $win $xp
    }
    method apply_values_ComboBox {} {
        
        foreach containerNode [array names entry_node_values] {
            if { [string first "," $containerNode] != -1 } { continue }
            set proc [$containerNode @values ""]
            if { $proc == ""} { continue }
            if { [string index $proc 0] eq "\[" && [string index $proc end] eq "\]"} {                      
                set values [split [get_domnode_attribute $containerNode values ""] ,]  
                lassign [lindex $entry_node_values($containerNode) 1] type widget 
                if { $type eq "ComboBox" } {
                    $widget configure -values $values
                }                                              
            }                             
        }  
    }
}

















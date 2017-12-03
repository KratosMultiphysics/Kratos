# Usage:
#     source {E:\PROYECTOS\FileWizard\NativeFileManager\NativeFileOpener.tcl}
#     FileSelector::InitWindow W
#     FileSelector::InitWindow puts
#     W $::FileSelector::files_to_model
#     puts W $::FileSelector::files_to_model
#     FileSelector::CopyFilesIntoModel {E:\TEMP\aaaa}

namespace eval FileSelector {
    variable selected_file
    variable save_to_model
    variable result_proc_name
    variable result_proc_args
    variable w
    variable files_to_model
    variable files_list
}

proc FileSelector::Start {} {
    variable selected_file
    set selected_file ""
    
    variable result_proc_name
    set result_proc_name ""
    variable result_proc_args
    set result_proc_args ""
    
    variable w
    set w .gid.fileSelector
    
    variable files_to_model
    set files_to_model [list ]
    
    variable files_list
    set files_list [list ]
}
FileSelector::Start

# PUBLIC FUNCTIONS
proc FileSelector::InitWindow {updateProcname args} {
    set ::FileSelector::selected_file ""
    set ::FileSelector::save_to_model 0
    variable result_proc_name
    variable result_proc_args
    set result_proc_name $updateProcname
    set result_proc_args $args
    FileSelector::_OpenFileSelector
}
proc FileSelector::FinishWindow {result} {
    variable result_proc_name
    variable result_proc_args
    
    if {$result} {
        variable save_to_model
        variable files_to_model
        variable selected_file
        variable files_list
        
        if {$save_to_model} {
            lappend files_to_model $selected_file
            set selected_file [file join . [file tail $selected_file] ]
        }
        lappend files_list $selected_file
        $result_proc_name $selected_file $result_proc_args
    } else {
        $result_proc_name "" $result_proc_args
    }
    catch {variable w; destroy $w}
}

proc FileSelector::CopyFilesIntoModel { dir } {
    variable files_to_model
    variable files_list
    foreach f $files_to_model {
        set files_list [lsearch -all -inline -not -exact $files_list $f]
        file copy -force $f $dir
        lappend files_list [file join $dir $f]
    }
    set files_to_model [list ]
}

proc FileSelector::GetAllFiles { } {
   variable files_list
   return $files_list
}
proc FileSelector::AddFile { fileid } {
   variable files_list
   if {$fileid ne "" && $fileid ni $files_list} {lappend files_list $fileid}
}


# PRIVATE FUNCTIONS
proc FileSelector::_OpenFileSelector2 { } {
    variable w
    catch {destroy $w}
    toplevel $w
    wm title $w "Select a file"
    wm minsize $w 400 200
    wm resizable $w 0 0
    # Top frame
    set fr1 [ttk::frame $w.fr1 -borderwidth 10]
    
    # Label
    set lab1 [ttk::label $fr1.lab1 -text {Filename: } -justify left -anchor w ]
    grid $lab1 -row 0 -column 0 -sticky ew 
    
    # Entry
    grid [ttk::entry $fr1.ent1 -width 40 -textvariable ::FileSelector::selected_file]  -column 1 -row 0 -sticky we; # -state readonly
    
    # Button browse
    grid [ttk::button $fr1.browse -text "Browse" -command "set ::FileSelector::selected_file \[tk_getOpenFile\]" ]  -column 2 -row 0 -sticky we
    
    # Checkbutton
    grid [ttk::checkbutton $fr1.check -text "Save file into model?" -variable ::FileSelector::save_to_model] -column 0 -row 1 -columnspan 3 -sticky we
    grid $fr1 -column 0 -row 0 -sticky nw
    
    # Bottom frame
    set fr2 [ttk::frame $w.fr2 -borderwidth 10]
    
    # Button browse
    grid [ttk::button $fr2.cancel -text "Cancel" -command [list FileSelector::FinishWindow 0 ]]  -column 0 -row 0 -sticky w
    grid [ttk::button $fr2.ok -text "Ok" -command [list FileSelector::FinishWindow 1] ] -column 1 -row 0 -sticky e
    
    grid $fr2 -column 0 -row 1 -sticky swe
}

proc FileSelector::_OpenFileSelector { } {
    variable w
    ::InitWindow $w [_ "Select a file"] PreFileSelectorWindowGeom FileSelector
    if { ![winfo exists $w] } return ;# windows disabled || usemorewindows == 0

     # Top frame
    set fr1 [ttk::frame $w.fr1 -borderwidth 10]
    
    # Label
    set lab1 [ttk::label $fr1.lab1 -text {Filename: } -justify left -anchor w ]
    grid $lab1 -row 0 -column 0 -sticky ew 
    
    # Entry
    grid [ttk::entry $fr1.ent1 -width 40 -textvariable ::FileSelector::selected_file]  -column 1 -row 0 -sticky we; # -state readonly
    
    # Button browse
    grid [ttk::button $fr1.browse -text "Browse" -command "set ::FileSelector::selected_file \[tk_getOpenFile\]" ]  -column 2 -row 0 -sticky we
    
    # Checkbutton
    grid [ttk::checkbutton $fr1.check -text "Save file into model?" -variable ::FileSelector::save_to_model] -column 0 -row 1 -columnspan 3 -sticky we
    
    grid $fr1 -column 0 -row 0 -sticky nw
   
    ttk::frame $w.but -style BottomFrame.TFrame   
    ttk::button $w.but.accept -text [_ "Apply"] -command "[list FileSelector::FinishWindow 1 ]"  -underline 0 -style BottomFrame.TButton   
    ttk::button $w.but.close -text [_ "Close"] -command "[list FileSelector::FinishWindow 0 ]" -underline 0 -style BottomFrame.TButton   


    grid columnconfigure $w.fr1 1 -weight 1

    grid $w.but.accept -row 1 -column 1 -padx 5 -pady 6
    grid $w.but.close -row 1 -column 3 -padx 5 -pady 6
    grid $w.but -row 4 -column 0  -sticky ews -columnspan 7
    if { $::tcl_version >= 8.5 } { grid anchor $w.but center }

    grid $w.but -row 3 -sticky ews -columnspan 7

    grid columnconfigure $w 0 -weight 1
    grid rowconfigure $w 3 -weight 1
    #
    ## Resize behavior management
    #wm minsize $w 180 200

    focus $w.but.accept
    bind $w <Alt-c> "$w.but.close invoke"
    bind $w <Escape> "$w.but.close invoke"
}
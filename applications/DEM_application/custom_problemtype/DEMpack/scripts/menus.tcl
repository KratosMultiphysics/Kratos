##############################################################################
#
#    NAME: menus.tcl
#
#    PURPOSE: TCL script to work with menus in preprocess and postprocess in the 
#             Kratos problem type
#
#    QUANTECH ATZ-DEVELOPMENT DEPARTMENT
#
#    AUTHOR : G. Socorro
#
#    CREATED AT: 29/03/06
#
#    HISTORY:
#
#     1.8-09/07/13-G. Socorro, add the Toolbar argument when call the proc OpenGiDGroupTab 
#     1.7-17/06/13-G. Socorro, modify the procs AddMenuToPreprocessMenu and CreatePreprocessModelTBar to use only the new GiD group options
#     1.6-12/02/13-G. Socorro, modify the command to access to the group window for version 11.0.x and 11.1.x
#     1.5-10/10/12-G. Socorro, change arc.gif by curves.gif
#     1.4-08/10/12-J. Garate, adapted to New GiD Groups
#     1.3-03/10/12-G. Socorro, update some menu option to use the new curve module
#     1.2-01/10/12-J. Garate, Enable/disable Curves Module
#     1.1-20/09/12-J. Garate, Add Curve's menu button
#     1.0-26/04/12-G. Socorro, change GiD_Groups by Cond_Groups
#     0.9-29/03/12-J. Garate, icon path change to adapt to GiD Themes
#     0.8-29/03/12-G. Socorro, update the call to the model and material window using "::KMProps::StartBaseWindow"
#     0.7-22/03/12-J. Gárate, Cambio a funciones públicas de los grupos de GiD
#     0.6-12/03/12-J. Gárate, Adaptacion a los nuevos grupos de GiD
#     0.4-01/03/12-J. Gárate, Añadidos los botones de recarga del PT, en modo Debugger
#     0.3-27/02/12-J. Gárate, Cambio el boton de Toolbar y del menu de "Grupos" a los grupos de GiD
#     0.2-08/06/10-G. Socorro, add more options to the Kratos preprocess toolbar
#     0.1-01/02/10-G. Socorro, create a base source code
#
###############################################################################

# kmtb => Kratos Menus and ToolBars

namespace eval ::kmtb:: {
}

proc ::kmtb::ChangePreprocessMenu {dir} {   
    # Disabled all base options
    GidChangeDataLabel "Data units" ""
    GidChangeDataLabel "Interval" ""
    GidChangeDataLabel "Conditions" ""
    GidChangeDataLabel "Materials" ""
    GidChangeDataLabel "Intervals" ""
    GidChangeDataLabel "Problem Data" ""
    GidChangeDataLabel "Local axes" ""
    
    GiDMenu::Create [= "Kratos#C#menu"] PRE -1 =
    GiDMenu::InsertOption [= "Kratos#C#menu"] [list [= "Model properties#C#menu"]] 0 PRE KMProps::StartBaseWindow "" "" replace =
    GiDMenu::InsertOption [= "Kratos#C#menu"] [list [= "Material database#C#menu"]] 1 PRE {KMProps::StartBaseWindow Materials} "" "" replace =
    GiDMenu::InsertOption [= "Kratos#C#menu"] "---" 2 PRE "" "" "" replace =
    GiDMenu::InsertOption [= "Kratos#C#menu"] [list [= "Project settings#C#menu"]] 3 PRE kps::InitSettingWindow "" "" replace =    
    
    GiDMenu::InsertOption "Help" [list ---] end PREPOST {} "" "" insertafter
    GiDMenu::InsertOption "Help" [list [_ "Visit %s web" $::kipt::ProgramName]...] end PREPOST [list VisitWeb $::kipt::Web] "" "" insertafter
    GiDMenu::InsertOption "Help" [list [concat [_ "About"] " " $::kipt::ProgramName]...] end PREPOST kipt::About "" "" insertafter
        
    # Update the menu properties
    ::GiDMenu::UpdateMenus
}

proc ::kmtb::CreatePreprocessModelTBar {{type "DEFAULT INSIDELEFT"}} {
    
    global KBitmapsNames KBitmapsCommands KBitmapsHelp 
    global KPriv
    set dir $::KPriv(dir)
    catch {unset KBitmapsNames KBitmapsCommands KBitmapsHelp}
    if {$KPriv(what_dempack_package) eq "C-DEMPack"} {     
        if {[kipt::CurvesModule]} {

            set KBitmapsNames(0) "$KPriv(imagesdir)/groups.gif $KPriv(imagesdir)/new_props.gif $KPriv(imagesdir)/maticon.gif $KPriv(imagesdir)/nonlinear.gif \
                                  --- $KPriv(imagesdir)/openrunsim.gif $KPriv(imagesdir)/runsimulation.gif $KPriv(imagesdir)/runsiminfo.gif \
                                  $KPriv(imagesdir)/stop.gif $KPriv(imagesdir)/wizard.gif"
            
            set KBitmapsCommands(0) [list \
                                    [list -np- ::WinUtils::OpenGiDGroupTab Toolbar] \
                                    [list -np- ::KMProps::StartBaseWindow] \
                                    [list -np- ::KMProps::StartBaseWindow Materials] \
                                    [list -np- ::KMProps::StartBaseWindow Curve] \
                                    "" \
                                    [list -np- RunWin] \
                                    {Utilities Calculate} \
                                    [list -np- PWViewOutput] \
                                    {Utilities CancelProcess} \
                                    [list -np- ::kmtb::OpenWiz]]

            set KBitmapsHelp(0) [list [= "Define the group properties using the group editor"] \
                                      [= "Define the model properties"] \
                                      [= "Define the material properties"] \
                                      [= "Define the curve properties"] \
                                      "" \
                                      [= "Open the process control window"] \
                                      [= "Run the simulation"] \
                                      [= "View process info"] \
                                      [= "Cancel process"] \
                                      [= "Launch Wizard"]] 
        } else {

            set KBitmapsNames(0) "$KPriv(imagesdir)/groups.gif $KPriv(imagesdir)/new_props.gif $KPriv(imagesdir)/maticon.gif \
                                 --- $KPriv(imagesdir)/openrunsim.gif $KPriv(imagesdir)/runsimulation.png $KPriv(imagesdir)/runsiminfo.gif \
                                  $KPriv(imagesdir)/stop.png $KPriv(imagesdir)/wizard.gif"
			
            set KBitmapsCommands(0) [list \
                                    [list -np- ::WinUtils::OpenGiDGroupTab Toolbar] \
                                    [list -np- ::KMProps::StartBaseWindow] \
                                    [list -np- ::KMProps::StartBaseWindow Materials] \
                                    "" \
                                    [list -np- RunWin] \
                                    {Utilities Calculate} \
                                    [list -np- PWViewOutput] \
                                    {Utilities CancelProcess}\
                                    [list -np- ::kmtb::OpenWiz]]
            
            set KBitmapsHelp(0) [list [= "Define the group properties using the group editor"] \
                                      [= "Define the model properties"] \
                                      [= "Define the material properties"] \
                                      "" \
                                      [= "Open the process control window"] \
                                      [= "Run the simulation"] \
                                      [= "View process info"] \
                                      [= "Cancel process"] \
                                      [= "Launch Wizard"]] 
        }
    } else {
        if {[kipt::CurvesModule]} {
            set KBitmapsNames(0) "$KPriv(imagesdir)/groups.gif $KPriv(imagesdir)/new_props.gif $KPriv(imagesdir)/maticon.gif $KPriv(imagesdir)/nonlinear.gif \
                                  --- $KPriv(imagesdir)/openrunsim.gif $KPriv(imagesdir)/runsimulation.png $KPriv(imagesdir)/runsiminfo.gif \
                                  $KPriv(imagesdir)/stop.png $KPriv(imagesdir)"
            
            set KBitmapsCommands(0) [list \
                                    [list -np- ::WinUtils::OpenGiDGroupTab Toolbar] \
                                    [list -np- ::KMProps::StartBaseWindow] \
                                    [list -np- ::KMProps::StartBaseWindow Materials] \
                                    [list -np- ::KMProps::StartBaseWindow Curve] \
                                    "" \
                                    [list -np- RunWin] \
                                    {Utilities Calculate} \
                                    [list -np- PWViewOutput] \
                                    {Utilities CancelProcess} ]

            set KBitmapsHelp(0) [list [= "Define the group properties using the group editor"] \
                                      [= "Define the model properties"] \
                                      [= "Define the material properties"] \
                                      [= "Define the curve properties"] \
                                      "" \
                                      [= "Open the process control window"] \
                                      [= "Run the simulation"] \
                                      [= "View process info"] \
                                      [= "Cancel process"] ]
        } else {

            set KBitmapsNames(0) "$KPriv(imagesdir)/groups.gif $KPriv(imagesdir)/new_props.gif $KPriv(imagesdir)/maticon.gif \
                                  --- $KPriv(imagesdir)/openrunsim.gif $KPriv(imagesdir)/runsimulation.png $KPriv(imagesdir)/runsiminfo.gif \
                                  $KPriv(imagesdir)/stop.png"
                
            set KBitmapsCommands(0) [list \
                                    [list -np- ::WinUtils::OpenGiDGroupTab Toolbar] \
                                    [list -np- ::KMProps::StartBaseWindow] \
                                    [list -np- ::KMProps::StartBaseWindow Materials] \
                                    "" \
                                    [list -np- RunWin] \
                                    {Utilities Calculate} \
                                    [list -np- PWViewOutput] \
                                    {Utilities CancelProcess}]
            
            set KBitmapsHelp(0) [list [= "Define the group properties using the group editor"] \
                                      [= "Define the model properties"] \
                                      [= "Define the material properties"] \
                                      "" \
                                      [= "Open the process control window"] \
                                      [= "Run the simulation"] \
                                      [= "View process info"] \
                                      [= "Cancel process"] ]
        }
    }
    # prefix values:
    # Pre        Only active in the preprocessor
    # Post       Only active in the postprocessor
    # PrePost    Active Always
    
    set prefix Pre
    set name KPreprocessModelbar
    set procname ::kmtb::CreatePreprocessModelTBar
    set KPriv(ToolBars,PreprocessModelTBar) [CreateOtherBitmaps $name [= "Model definition toolbar"] KBitmapsNames KBitmapsCommands KBitmapsHelp $dir $procname $type $prefix]
    
    AddNewToolbar [= "Model definition toolbar"] ${prefix}${name}WindowGeom $procname 
}

proc ::kmtb::EndCreatePreprocessTBar {} {
    global KPriv
   
    set name KPreprocessModel
    set winname ::kmtb::CreatePreprocessModelTBar
    catch { 
        ReleaseToolbar ${name}bar
        rename $winname ""
        destroy $KPriv(ToolBars,PreprocessModelTBar)
        unset KPriv(ToolBars,PreprocessModelTBar)
    }
}

proc ::kmtb::UpdateVersion {} {
    return [GiD_Process MEscape data defaults TransfProblem]
}

proc ::kmtb::OpenWiz {} {
    ::kwiz::CreateDEMWizard
}

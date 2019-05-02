###############################################################################
#
#        NAME: kmaterials.tcl
#
#        PURPOSE: Main window to manage material properties
#
#        QUANTECH ATZ-DEVELOPMENT DEPARTMENT
#
#        AUTHOR : G. Socorro
#
#        CREATED AT: 25/02/2010
#
#        HISTORY:
#
#   1.9- 17/07/13-G. Socorro, modify the proc FillTreeMat to take into account the current state of the property (normal,hidden or disabled)
#   1.8- 15/07/13-G. Socorro, create a local variable createframeafteredit to disable/enable the bottom frame creation after rename a material
#   1.7- 15/07/13-G. Socorro, update the proc insertXml and insertXmlCopy
#        1.6- 11/07/13-G. Socorro, correct the bug in the proc splitNode from the last Adrià modification
#        1.5- 20/09/12-J. Garate, Minor bug fixing, 
#        1.4- 24/07/12-J. Garate, Minor bug fixing, function comments 
#        1.3- 20/07/12-J. Garate, Materials Tab Frame Combos with multiple values 
#        1.2- 19/07/12-J. Garate, Materials Tab Frame and service functions 
#        1.1- 12/06/12-J. Garate, Materials Tab Frame and service functions 
#                                        (::KMat::CreateTabTree ::KMat::Dict_create ::KMat::GetParentNB 
#                                        ::KMat::IsNum ::KMat::GetNBPath ::KMat::addTabtoNb ::KMat::Combos)
#        1.0- 14/05/12-J. Garate, corrected bug, material combobox
#        0.9- 30/03/12-G. Socorro, add the variables TreeMatsPath and NbMatsPath
#        0.8- 22/06/11-G. Socorro, delete snit, tdom and xmlstruct from the package require
#        0.7- 07/06/11 GS, add composite and plastic material to the structural analysis application
#        0.6- 27/09/10 LC, Correct bugs inserting new materials and double click in containers
#        0.5- 08/06/10 KS, New materials database structure + template for add material.
#        0.4- 11/05/10 KS, varius bugs fixed.
#        0.3- 26/04/10 KS, Material Groups.
#        0.2- 20/04/10 KS, Read, Add, Remove, Rename OK.
#        0.1- 25/02/10 LCA, create a base source code from the kegroups.tcl script
#
###############################################################################

package require treectrl
package require tooltip
package provide KMat 1.0 


# Create a base namespace KMat
namespace eval ::KMat:: {
    
    # Path of the base window 
    variable WinPath ".gid.kmprops"
    # Tree material properties path
    variable TreeMatsPath
    # Notebook material properties path
    variable NbMatsPath
    
    variable WinLayout 
    variable SystemHighlight
    variable SystemHighlightText
    
    variable lastSelected {}
    variable abdlist
    variable xml ""
    
    # Se inicializan las clases dinámicamente leyendo del xml
    variable visibilityVars {}
    
    
}

set ::KMat::UseFullTkTree 1

# ERROR: Using fulltktree gives an error when hidding the tree again!

if { $::KMat::UseFullTkTree} {
    package require fulltktree
}

proc ::KMat::Init {} {
    
    variable WinLayout;        variable SystemHighlight
    variable SystemHighlightText; variable abdlist
    
    # Get default colors
    set w [listbox .listbox]
    set SystemHighlight [$w cget -selectbackground]
    set SystemHighlightText [$w cget -selectforeground]
    destroy $w
    if { $::tcl_platform(platform) == "unix" } {
        # I hate that gray selection color
        set SystemHighlight #316ac5
        set SystemHighlightText White
    }
    
    set WinLayout "OUTSIDE"
    
    set abdlist [list A11 A12 A13 A21 A22 A23 A31 A32 A33 \
            B11 B12 B13 B21 B22 B23 B31 B32 B33 \
            D11 D12 D13 D21 D22 D23 D31 D32 D33 ]
    
    global KPriv
    if {[info exists KPriv(xmlMat)]} {
        set ::KMat::xml $KPriv(xmlMat)
    }
    if {[info exists ::MaterialTrees)]} {
        ::struct::tree ::MaterialTrees
    }
    
}

# Init KMat namespace
::KMat::Init

proc ::KMat::CreateTreeAndToolbar { w } {
    
    # Create the treectrl properties
    set mdf [ttk::frame $w.middle ]
    set T [::KMat::CreateTreeProperties $w]
    
    grid $w.middle -sticky wens
    
    # Create the frame where set the properties
    set f [ttk::frame $w.fBottom -borderwidth 0 ]
    
    # Grid for toolbar
    grid $f -row 2 -column 0 -sticky wes
    
    # Create the toolbar frame
    set tbf [ttk::frame $w.tbar -borderwidth 0]
    
    ttk::button $tbf.newmaterial -image [::WinUtils::GetImage add.gif]  -command [list ::KMat::CreateNewMaterial $T] -style Toolbutton
    tooltip::tooltip $tbf.newmaterial [= "Create new material"]
    
    ttk::button $tbf.deleteMaterialsId -image [::WinUtils::GetImage delete_icon.gif]  -command [list ::KMat::DeleteMaterial $T] -style Toolbutton
    tooltip::tooltip $tbf.deleteMaterialsId [= "Delete material"]
    
    ttk::button $tbf.refreshTree -image [::WinUtils::GetImage refresh.gif]  -command [list ::KMat::refreshTree] -style Toolbutton
    tooltip::tooltip $tbf.refreshTree [= "Refresh Tree"]
    
    #        ttk::button $tbf.testbutton -image [::WinUtils::GetImage uc.gif]  -command [list ::KMat::testbutton] -style Toolbutton
    #        tooltip::tooltip $tbf.testbutton [= "tests"]
    
    
    #set tf [ttk::frame $w.close] 
    grid $tbf -sticky ews
    grid anchor $tbf w
    #grid [ttk::button $tf.bClose -text [= "Close"] -command [list destroy $w]]  -sticky ew -padx 5 -pady 3
    
    # Grid for toolbar 
    #grid $w.tbar -row 2 -column 0 -sticky wes
    grid $tbf.newmaterial -sticky we -row 0 -column 0
    grid $tbf.deleteMaterialsId -sticky we -row 0 -column 1
    grid $tbf.refreshTree -sticky we -row 0 -column 2
    #        grid $tbf.testbutton -sticky we -row 0 -column 3
    
    #grid rowconfigure $w 2 -weight 1
    #grid columnconfigure $w 0 -weight 1
    
    focus $T
    return $T
}

proc ::KMat::CreateTreeProperties {w} {
    
    variable SystemHighlight
    variable SystemHighlightText
    
    # Scrollbars
    set vsb $w.vsb1
    set hsb $w.hsb1
    
    
    # Create the treectrl and set the scrollbars
    if { $::KMat::UseFullTkTree} {
        if { [ catch {set T [fulltktree $w.t]}]} {
            set T [fulltktree $w.t]
        }
    } else {
        set T [treectrl $w.t -xscrollcommand [list $hsb set] -yscrollcommand [list $vsb set]]
        ttk::scrollbar $vsb -orient vertical -command [list $T yview]
        ttk::scrollbar $hsb -orient horizontal -command [list $T xview]
    }
    
    
    # Set the height
    set height [font metrics [$T cget -font] -linespace]
    if {$height < 19} {
        set height 19
    }
    
    # Configure the treectrl
    $T configure -indent 15 -itemheight $height -selectmode browse \
        -showroot 0 -showrootbutton 0 -showbuttons 1 -showlines 1 \
        -highlightthickness 0 -borderwidth 0 -height 300 \
        -xscrollincrement 20 -yscrollincrement 20 
    
    $T column create -text [= "Materials"] -tags C0 -weight 0
    
    # Configure the column weight and arrow
    $T column configure C0 -weight 1 -arrow up
    
    # Configure the column that have the tree
    $T configure -treecolumn C0
    
    # Intento de desactivar el drag
    #$T column dragconfigure -enable no
    
    # Create elements
    $T element create elemImgAny image
    $T element create elemTxtRead text -fill [list $SystemHighlightText {selected focus}] -lines 1
    $T element create elemRectSel rect -fill [list $SystemHighlight {selected focus} gray {selected !focus}] -showfocus yes
    $T element create eWindow window
    
    # Create styles using the elements
    set S [$T style create styAnyRead]
    $T style elements $S {elemImgAny elemRectSel elemTxtRead}
    $T style layout $S elemImgAny -expand ns
    $T style layout $S elemTxtRead -padx 4 -expand ns -squeeze x
    $T style layout $S elemRectSel -union [list elemTxtRead] -iexpand ns -ipadx 2
    
    set S [$T style create styFrame -orient horizontal]
    $T style elements $S {eWindow}
    $T style layout $S eWindow -squeeze x -padx {0 1} -pady {2 2}
    
    # Items
    set item root
    $T item configure $item -button yes
    $T item style set $item C0 styAnyRead
    $T item element configure $item C0 elemTxtRead -text [= "Root"]
    
    # List of lists: {column style element ...} specifying text elements the user can edit
    TreeCtrl::SetEditable $T {
        
    }
    
    # List of lists: {column style element ...} specifying elements
    # the user can click on or select with the selection rectangle
    TreeCtrl::SetSensitive $T {
        {C0 styAnyRead elemTxtRead}
    }
    
    # List of lists: {column style element ...} specifying elements
    # added to the drag image when dragging selected items
    TreeCtrl::SetDragImage $T {
        
    }
    
    # Some notify install
    #$T notify install <Drag-receive>
    $T notify install <Edit-accept>
    
    # Notify bind
    # TODO
    
    #                $T notify bind DragTag <Drag-receive> { ::KMProps::ReceiveDragGroups  %T %l %I }
    $T notify bind EditTag <Edit-accept> { ::KMat::SetMatToRename %T %I %t }
    
    bind $T <Button-1> [list ::KMat::ClickTree %x %y $T]
    bind $T <Double-Button-1> [list ::KMat::DoubleClickTree %x %y $T]
    #         bind $T <Return> [list SetLayersTo TOUSE $T]
    #         bind $T <Key-Delete> [list SetLayersToDelete $T]
    #         bind $T <Alt_L> [list InvertSelectionTableList $T]
    #         bind $T <Alt_R> [list InvertSelectionTableList $T]
    #         bind $T <Meta_L> [list InvertSelectionTableList $T]
    #         bind $T <Meta_R> [list InvertSelectionTableList $T]
    bind $T <F2> [list ::KMat::BeginEditMaterial $T]
    
    bind $T <Button-3> "[list ::KMat::MenuContextualGroup %W %x %y] ; break"
    # MouseWheel
    if {[string equal "x11" [tk windowingsystem]]} {
        # Support for mousewheels on Linux/Unix commonly comes through mapping
        # the wheel to the extended buttons.  If you have a mousewheel, find
        # Linux configuration info at:
        #        http://www.inria.fr/koala/colas/mouse-wheel-scroll/
        bind $T <4> [list $T yview scroll -3 units ]
        bind $T <5> [list $T yview scroll 3 units ]
    } elseif {[string equal [tk windowingsystem] "aqua"]} {
        bind $T <MouseWheel> [subst -nocommands { $T yview scroll {eval [expr - (%D)]} units } ]
    } else {
        bind $T <MouseWheel> [subst -nocommands { $T yview scroll [expr - (%D / 120) * 3] units } ]    
    }
    # Grid the tree
    if { $::KMat::UseFullTkTree} {
        grid $T -sticky nsew
    } else {
        grid $T $vsb -sticky nsew
        grid configure $vsb -sticky ns
        grid $hsb -sticky ew
        grid remove $vsb $hsb
        bind $T <Configure> [list ConfigureListScrollbars $T $hsb $vsb]
    }
    grid columnconfigure $w 0 -weight 1
    grid rowconfigure $w 0 -weight 1
    
    bindtags $T [list $T TreeCtrlFileList TreeCtrl [winfo toplevel $T] all]
    
    return $T
}

#---------------------------------------------------------------------------------------------- 
# Lee el xml y carga el árbol de propiedades de forma recursiva (sin maximo de niveles)
#----------------------------------------------------------------------------------------------

proc ::KMat::FillTreeMatRecursive { } {
    variable abdlist; variable TreeMatsPath
    global KPriv
    
    set KPriv(materialsId) {}
    
    set T $TreeMatsPath
    set nodes [$::KMat::xml selectNodes "/Kratos_KMat_DB/Materials/MaterialGroup\[@id\]"]
    
    foreach node $nodes {
        ::KMat::FillNextLevel 0 $node $T $node
    }
}

proc ::KMat::FillNextLevel { depth acumPath T node } {
    
    global KPriv
    
    set ParentPath "root"
    set parent ""
    
    if {$depth != 0} {
        set parent [lrange $acumPath 0 0]
        set ParentPath [$parent getAttribute id ""]
        foreach elem [lrange $acumPath 1 end] {
            set split [::KMat::splitNode $elem]
            append ParentPath "//$split"
        }
        lappend acumPath $node
        
        set item [::KMat::InsertNewItem [$node getAttribute pid ""] [$node getAttribute id ""] $T $parent $ParentPath]
        set nodes [$node childNodes]
        incr depth 1
        if {$item != -1} {
            foreach node2 $nodes {
                ::KMat::FillNextLevel $depth $acumPath $T $node2
            }
        }
    } else {
        set item [::KMat::InsertNewItem [$node getAttribute pid ""] [$node getAttribute id ""] $T "" "root" [$node hasChildNodes] [::KMat::stateNode $node] [$node getAttribute open "0"]]
        set nodes2 [$node childNodes]
        foreach node2 $nodes2 {
            ::KMat::FillNextLevel 1 $acumPath $T $node2
        }
    }
}

#---------------------------------------------------------------------------------------------- 
# Lee el xml y carga el árbol de propiedades de forma iterativa como máximo hasta 7 niveles -> Recursiva (Proximamente)
#----------------------------------------------------------------------------------------------
proc ::KMat::FillTreeMat { } {
    
    variable abdlist; variable TreeMatsPath
    global KPriv
    
    set KPriv(materialsId) {}
    
    set T $TreeMatsPath
    set nodes [$::KMat::xml selectNodes "/Kratos_KMat_DB/Materials/MaterialGroup\[@id\]"]
    
    foreach node $nodes {                          
        # Nos guardamos todos los Id
        set pid        [$node getAttribute pid ""]
        set id_l1 [$node getAttribute id ""]
        set open [$node getAttribute open "0"]
        set state [::KMat::stateNode $node]
        # wa "Level 1: pid:$pid id_l1:$id_l1 open:$open state:$state\n"
        set item [::KMat::InsertNewItem $pid $id_l1 $T "" "root" [$node hasChildNodes] $state $open]
        # wa "item:$item"
        foreach node2 [$node childNodes] {
            set pid [$node2 getAttribute pid ""]
            set id [$node2 getAttribute id ""]
            set open [$node2 getAttribute open "0"]
            set state [::KMat::stateNode $node2]
            set node2split [::KMat::splitNode $node2]
            set cpath "$id_l1//"
            # wa "Level 2: pid:$pid id:$id open:$open state:$state node2split:$node2split\n"
            set item2 [::KMat::InsertNewItem $pid $node2split $T $cpath "$item" [$node2 hasChildNodes] $state $open]
            # wa "item2:$item2"
            if {$item2 != -1} {                                        
                lappend KPriv(materialsId) $pid
                
                # Seleccionamos los hijos (3º nivel)
                foreach node3 [$node2 childNodes] {
                    set pid [$node3 getAttribute pid ""]
                    set id [$node3 getAttribute id ""]
                    set open [$node3 getAttribute open "0"]
                    set state [::KMat::stateNode $node3]
                    set node3split [::KMat::splitNode $node3]
                    set cpath "$id_l1//$node2split//"
                    # wa "Level 3: pid:$pid id:$id open:$open state:$state node3split:$node3split cpath:$cpath"
                    set item3 -1
                    if {$pid ni $abdlist } {
                        set item3 [::KMat::InsertNewItem $pid $node3split $T $cpath "$item2" [$node3 hasChildNodes] $state $open]
                    }
                    if {$item3 != -1} {
                        foreach node4 [$node3 childNodes] {
                            set pid [$node4 getAttribute pid ""]
                            set id [$node4 getAttribute id ""]
                            set open [$node4 getAttribute open "0"]
                            set state [::KMat::stateNode $node4]
                            set node4split [::KMat::splitNode $node4]
                            set cpath "$id_l1//$node2split//$node3split//"
                            # wa "Level 4: pid:$pid id:$id open:$open state:$state node4split:$node4split cpath:$cpath"
                            set item4 [::KMat::InsertNewItem $pid $node4split $T "$cpath" "$item3" [$node4 hasChildNodes] $state $open]
                            if {$item4 != -1} {
                                foreach node5 [$node4 childNodes] {
                                    set pid [$node5 getAttribute pid ""]
                                    set id [$node5 getAttribute id ""]
                                    set open [$node5 getAttribute open "0"]
                                    set state [::KMat::stateNode $node5]
                                    set node5split [::KMat::splitNode $node5]
                                    set cpath "$id_l1//$node2split//$node3split//$node4split//"
                                    # wa "Level 5: pid:$pid id:$id open:$open state:$state node5split:$node5split cpath:$cpath"
                                    set item5 [::KMat::InsertNewItem $pid $node5split $T "$cpath" "$item4" [$node5 hasChildNodes] $state $open]
                                    if {$item5 != -1} {
                                        foreach node6 [$node5 childNodes] {
                                            set pid [$node6 getAttribute pid ""]
                                            set id [$node6 getAttribute id ""]
                                            set open [$node6 getAttribute open "0"]
                                            set state [::KMat::stateNode $node6]
                                            set node6split [::KMat::splitNode $node6]
                                            set cpath "$id_l1//$node2split//$node3split//$node4split//$node5split//"
                                            # wa "Level 6: pid:$pid id:$id open:$open state:$state node6split:$node6split cpath:$cpath"                
                                            set item6 [::KMat::InsertNewItem $pid $node6split $T "$cpath" "$item5" [$node6 hasChildNodes] $state $open]
                                            if {$item6 != -1} {
                                                foreach node7 [$node6 childNodes] {                        
                                                    set pid [$node7 getAttribute pid ""]
                                                    set id [$node7 getAttribute id ""]
                                                    set open [$node7 getAttribute open "0"]
                                                    set state [::KMat::stateNode $node7]
                                                    set node7split [::KMat::splitNode $node7]
                                                    set cpath "$id_l1//$node2split//$node3split//$node4split//$node5split//$node6split//"
                                                    # wa "Level 7: pid:$pid id:$id open:$open state:$state node7split:$node7split cpath:$cpath"
                                                    set item7 [::KMat::InsertNewItem $pid $node7split $T "$cpath" "$item6" [$node7 hasChildNodes] $state $open]
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}


proc ::KMat::DoubleClickTree {x y T {item ""}} {
    # Llegamos aqui cuando el usuario hace doble click en el arbol
    variable lastSelected
    
    #Si no llega directamente el item, miramos cual ha sido pulsado
    if { $item == "" } {  
        set info [$T identify $x $y]
        if { [lindex $info 0] == "item" && [llength $info] >= 4 } {
            set item [lindex $info 1]
        } else {
            return ""
        }
    }
    
    # Conseguimos el path, el id y el valor
    set fullname [DecodeName [$T item tag names $item]]
    set idFull [string map { "." "" "//" ""} $fullname]        
    set id [::KMat::setXml $fullname id ]
    set value [::KMat::setXml $fullname value "" ]
    
    set f "$T.f$idFull"
    if {[winfo exists $f]} {
        destroy $f
    }
    
    if { $value != "" && [llength $lastSelected] > 0 && [lindex $lastSelected 0] > $item} {
        # Si se trata de un Item (Property) entramos aqui
        ::KMat::cmbSelectChange [lindex $lastSelected 0] $T 1 "anterior"
    }
    set f "$T.f$idFull"
    if {[winfo exists $f]} {
        destroy $f
    }
    
    set f [::KMat::buildFrame $T $item]
    
    #Solo es necesario seguir si se trata de un item editable
    if { $f == "" } {
        return ""
    }
    
    $T item style set $item C0 styFrame
    $T item element configure $item C0 eWindow -window $f                
    ::KMat::insertIcon $item $T $fullname
    
    #Nos guardamos este item y el valor seleccionado
    set f "$T.f$idFull.cmb"
    set xpath "[::KMat::setXPath $fullname]"
    set selCombo [::xmlutils::getComboValue $::KMat::xml $xpath $f]
    set selComboText [::xmlutils::getComboValue $::KMat::xml $xpath $f "text"]
    set lastSelected [list $item $selCombo $selComboText]
    
    set path [DecodeName [$T item tag names $item]]   
    set splitted [::KEGroups::split2 $path //]        
    
    if { [llength $splitted] == 3 } {
        # Si tenemos que construir el arbol de Tabs, entraremos por aqui
        set clase "Tab"
        ::KMat::buildTabFrame $T $item $f $clase
        return ""
    }
    
    
}

proc ::KMat::buildfatherFrame { T item {class "Tab"} } {
    
    variable abdlist
    
    set fullname [DecodeName [$T item tag names $item]]
    set id [::KMat::setXml $fullname id "" ]
    set pid [::KMat::setXml $fullname pid "" ]
    set f [::KMat::iniFrameBottom]
    # Miramos los descendientes directos y si son container ponemos un tab por cada uno q tenga items
    set children [$T item children $item]
    
    set listTabs {}
    set listItems {}
    set acceptItems {}
    
    set nb ${f}.nb
    
    grid [ttk::notebook $nb ] -row 0 -column 0 -columnspan 2 -padx 0 -sticky nswe -in $f
    
    # Declaramos un tab para el material
    set fTab ${nb}.f$id
    $nb add [ttk::labelframe $fTab -text "[= Properties]" -padding {10 0 10 10}] \
        -text "[string range $pid 0 20]"
    
    set count 0
    foreach itemChild $children {
        set fullname [DecodeName [$T item tag names $itemChild]]
        set nodeName [::xmlutils::getXmlNodeName $fullname "mat"]
        
        set xpath "[::KMat::setXPath $fullname]"
        set comboList [::xmlutils::getValues $::KMat::xml $xpath]
        
        lappend listItems $itemChild
        
        set id [::KMat::setXml $fullname id "" ]
        set pid [::KMat::setXml $fullname pid "" ]
        set value [::KMat::setXml $fullname value "" ]
        set help [::KMat::setXml $fullname help "" ]
        
        grid [ttk::label $fTab.label$id -text [= $pid] ] \
            -row $count -column 0 -pady 5 -sticky nw -in $fTab
        grid [ttk::combobox $fTab.cmb$id -state normal -values $comboList -textvariable "::KMat::cmb$id" ] \
            -row $count -column 1 -padx 3 -pady 5 -sticky nw -in $fTab
        #bind $fTab.cmb$id <<ComboboxSelected>> "::KMat::cmbSelectChange $item $T"
        
        set ::KMat::cmb$id $value
        incr count
    }
    
    set clase [lappend fTab $listItems]
    #$f.bBottomOk and $f.bBottomCancel names need to popup click ok cancel before change item
    grid [ttk::button $f.bBottomOk -text "Ok"  -command "[list ::KMat::acceptTabFrame $T $acceptItems $clase $item]" ] \
        -row $count -column 0 -sticky sw  -pady 3 -padx 20  -in $f
    tooltip::tooltip $f.bBottomOk [= "Confirm values"]
    
    grid [ttk::button $f.bBottomCancel -text "Cancel"  -command "::KMat::cancelBottom" ] \
        -row $count -column 0 -sticky sw  -pady 3 -padx 100 -in $f
    tooltip::tooltip $f.bBottomCancel [= "Cancel assignation"]
}

proc ::KMat::Tree_Return { } {
    ::KMProps::DestroyBottomFrame     
}

proc ::KMat::buildTabFrame { T item f {class "Tab"} } {
    # Crea el arbol de Tabs
    
    global KPriv
    
    # Miramos los descendientes directos y si son container ponemos un tab por cada uno q tenga items
    set children [$T item children $item]
    set listTabs {}
    set acceptItems {}
    if { [info exist ::Tabdict] } {
        unset ::Tabdict
    }
    if { [info exist ::Itemsdict] } {
        unset ::Itemsdict
    }
    set ::Tabdict [dict create]
    set ::Itemsdict [dict create]
    set listTabs [::KMat::Dict_create $T $item]
    
    # listTabs contiene int a int b string c
    # Donde a es el subnivel que tiene
    # Donde b es el identificador de su padre en el subnivel superior
    # Donde c es su nombre
    
    # Tabdict contiene el id del tab, su path del xml
    # Itemsdict contiene el id del Item y su path en el arbol
    
    #Reseteamos la variable que nos indica si estamos en una propiedad
    set ::KMat::propertyName ""
    if { [llength $listTabs] >= 1 } {
        # Si tenemos Tabs que añadir
        set row 0
        
        set rootname [DecodeName [$T item tag names $item]]
        
        set xpath [::xmlutils::setXPath $rootname "mats"]
        set node [$KPriv(xmlMat) selectNodes $xpath]
        
        # Vamos a crear el arbol de tabs, con los Items y los combobox
        set row [::KMat::CreateTabTree $node $listTabs $T $f]
        # Aqui ya se han creado los Notebooks y los combos
        
        
        # Si pulsan Esc también forzamos la salida del Tab
        bind $T <Return> [list ::KMat::Tree_Return]
        
        incr row 1
        
        #Finalmente colocamos los Botones de OK y de Cancel debajo de todo
        #$f.bBottomOk and $f.bBottomCancel names need to popup click ok cancel before change item
        grid [ttk::button $f.bBottomOk -text "Ok"  -command "[list ::KMat::acceptTabFrame $T $acceptItems $class $item]" ] \
            -row $row -column 0 -sticky sw  -pady 3 -padx 20  -in $f
        tooltip::tooltip $f.bBottomOk [= "Confirm values"]
        
        grid [ttk::button $f.bBottomCancel -text "Cancel"  -command "::KMat::DestroyBottomFrame" ] \
            -row $row -column 0 -sticky sw  -pady 3 -padx 100  -in $f
        tooltip::tooltip $f.bBottomCancel [= "Cancel assignation"]
    }
}

proc ::KMat::CreateTabTree {node listTabs T f} {
    # Crea los Notebooks con tabs y coloca los Combos dentro de cada tab
    
    set stopid [$node getAttribute id ""]
    set createdLevels {}
    # stopid es el id del nodo que hemos seleccionado
    
    foreach { level parent id } $listTabs {
        # Nos recorremos la lista de tabs y los vamos colocando en los Notebooks correspondientes
        
        if { $level ni $createdLevels} {
            # Si el notebook aun no ha sido creado lo creamos
            lappend createdLevels $level
            if {$level == 0} {
                # Tratamos el nb0, que es especial ya que se acopla sobre f
                set nb ${f}.nb0
                if {[ winfo exists $nb] } {
                    destroy $nb
                }
                set nbPath [ttk::notebook $nb ]
                grid $nbPath -row $level -sticky nsew
                ttk::notebook::enableTraversal $nb
                
            } else {
                # Tratamos los demas nb, que se acoplan sobre el tab correspondiente del nb anterior
                set nb [::KMat::GetNBPath $f $level $parent $listTabs $stopid $id]
                
                if {[ winfo exists $nb] } {
                    destroy $nb
                }
                grid [ttk::notebook $nb ] -row $level -sticky nsew -in [::KMat::GetParentNB $nb 1]
                ttk::notebook::enableTraversal $nb
            }
        }
        # Finalmente añadimos el tab
        ::KMat::addTabtoNb $f $level $parent $id $listTabs $stopid
    }
    return 3
}

proc ::KMat::Dict_create {T item {i 0} {pi 0} } {
    # Llena los diccionarios y devuelve la Lista de Tabs, con la siguiente estructura
    # listTabs contiene int a int b string c
    # Donde a es el subnivel que tiene
    # Donde b es el identificador de su padre en el subnivel superior
    # Donde c es su nombre
    
    # Tabdict contiene el id del tab, su path del xml
    # Itemsdict contiene el id del Item y su path en el arbol
    set children [$T item children $item]
    set listTabs ""
    set mi 0
    foreach itemChild $children {
        set fullname [DecodeName [$T item tag names $itemChild]]
        set nodeName [::xmlutils::getXmlNodeName $fullname "mats"]
        #Miramos si cada hijo es container
        if { $nodeName == "Container" } {
            
            if { [$T item numchildren $itemChild] > 0 } {
                set Tabid [::xmlutils::setXml $fullname id "read" "" "mat"]
                
                lappend listTabs $i $pi $Tabid
                set aux [::KMat::Dict_create $T $itemChild [expr $i+1] $mi]
                foreach {aux1 aux2 aux3} $aux {
                    lappend listTabs $aux1 $aux2 $aux3
                }
                dict set ::Tabdict $Tabid $fullname
            }
        } else {
            set Tabid [::xmlutils::setXml $fullname id "read" "" "mat"]
            dict set ::Itemsdict $Tabid $fullname
        }
        incr mi 1
    }
    return $listTabs
}

proc ::KMat::GetParentNB { nb {cut 2} } {
    # Devuelve la direccion del Notebook del padre de $nb
    set id [string range $nb end end]
    set anb $nb
    if {$id > 0} {
        set anb [::KMat::CheckNB $anb]
        for { set i 1} {$i <= $cut} {incr i 1} {
            while { [string range $anb end end] != "." } {
                set anb [string range $anb 0 end-1]
            }
            set anb [string range $anb 0 end-1]
        }
    }
    return $anb
}

proc ::KMat::CheckNB { nb } {
    # Comprueba si $nb es un path valido para un Notebook
    
    set id [string range $nb end end]
    set anb $nb
    if {$id > 0} {
        set i 0
        while { ![::KMat::IsNum [string range $anb end end]] } {
            incr i 1
            set anb [string range $anb 0 end-1]
        }
        set row [string range $anb end end]
        set i1 $i
        set anb [string range $anb 0 end-1]
        while { ![::KMat::IsNum [string range $anb end end]] } {
            incr i 1
            set anb [string range $anb 0 end-1]
        }
        
        if { $row == [string range $anb end end] } {
            set nb [string range $nb 0 end-[expr $i1 +4]]
        }
    }
    
    if {![::KMat::IsNum [string range $nb end end]]} {
        set i 0
        set anb $nb
        while { [string range $anb end end] != "." } {
            incr i 1
            set anb [string range $anb 0 end-1]
        }
        set i1 $i
        set anb [string range $anb 0 end-1]
        if { ![::KMat::IsNum [string range $anb end end]]} {
            while { [string range $anb end end] != "." } {
                set anb [string range $anb 0 end-1]
            }
            set anb [string range $anb 0 end-1]
            set apended [string range $nb end-$i1 end]
            append anb $apended
            set nb $anb
        }
    }
    return $nb
}

proc ::KMat::IsNum { n } {
    # Devuelve 1 si $n es un numero, 0 si no
    
    if { $n >= 0 } {
        if { $n <= 9 } {
            return 1
        }
    }
    return 0
}

proc ::KMat::GetNBPath { f level parent listTabs stopid tid} {
    # Dado el id de un tab, devuelve su path dentro de los Notebook
    
    set fullname [dict get $::Tabdict $tid]
    set splitted [::KMProps::split2 $fullname //]
    set i 0
    set j 0
    set nb "${f}.nb0"
    
    foreach check $splitted {
        if { [string range $check 2 end] == $stopid} {
            break
        }
        incr j 1
    }
    incr j 1
    while {$i < $level } {
        set pid [string range [lrange $splitted [expr $j+$i] [expr $j+$i]] 2 end]
        set aux [append nb ".f$pid"]
        incr i 1
        set nb [append aux ".nb$i"]
    }
    return $nb
}

proc ::KMat::addTabtoNb { f level parent tid listTabs stopid} {
    # Busca en que nb tiene que ir el tid (tab a colocar) y lo inserta
    global KPriv
    
    set nb [::KMat::GetNBPath $f $level $parent $listTabs $stopid $tid]
    # Busca en que Notebook tenemos que insertar el tab
    
    if {$tid in [dict keys $::Tabdict]} {
        # Nos aseguramos que el tab está en el diccionario antes creado
        
        set fullname [dict get $::Tabdict $tid]
        set pid [::xmlutils::setXml $fullname pid "read" "" "mat"]
        set help [::xmlutils::setXml $fullname help "read" "" "mat"]
        set kids 0
        # Obtenemos los datos para añadir los tab
        # kids nos indicará si hay que añadir combos o no
        
        set fTab ${nb}.f
        # Preparamos el nombre del Tab, tiene que empezar por minuscula, con lo que ponemos una f delante
        
        set xpath [::xmlutils::setXPath $fullname "mats"]
        set nodeaux [$KPriv(xmlMat) selectNodes $xpath]
        set i 0
        foreach chnode $nodeaux {
            if {$pid != [$chnode getAttribute id ""] } {
                incr i 1
            } else {
                break
            }
        }
        set aux [lindex [$nodeaux childNodes ] $i]
        
        if { [$aux hasChildNodes] } {
            # Si dentro de este tab tenemos que meter otro Notebook, añadimos el ttk:frame
            set fTab [append fTab $tid]
            if {[winfo exists $fTab]} {
                destroy $fTab
            }
            set framePath [ttk::frame $fTab -padding {10 10 10 10} ]
            
        } else {
            # Si dentro de este tab tenemos que meter los combobox, añadimos el ttk:labelframe
            set fTab [append fTab $tid]
            if {[winfo exists $fTab]} {
                destroy $fTab
            }
            set framePath [ttk::labelframe $fTab -padding {10 10 10 10} -text [= "Properties"] -labelanchor nw]
            set kids 1
        }
        
        $nb add $framePath -text "[string range [= $pid] 0 20]" -sticky nsew
        # Añadimos el tab
        if { $kids } {
            # Añadimos los Combos
            ::KMat::Combos $chnode $level $framePath
        }
    }
}

proc ::KMat::Combos {node row framePath} {
    # Cuando un tab es padre (sus hijos son Items, no containers" creamos los combos
    
    global KPriv
    incr row 1
    
    set listContainer [$node childNodes]
    # Los combos a añadir son los hijos del nodo
    
    for {set i 0} { $i < [llength $listContainer] } {incr i} {
        # Para cada hijo, obtenemos sus datos
        set id [[lindex $listContainer $i] getAttribute id ""]
        set pid [[lindex $listContainer $i] getAttribute pid ""]
        set value [[lindex $listContainer $i] getAttribute value ""]
        set help [[lindex $listContainer $i] getAttribute help ""]
        
        if { [winfo exists $framePath.lbl$id ] } {
            destroy $framePath.lbl$id
        }
        
        # Para cada item añadimos label y combo
        grid [ttk::label $framePath.lbl$id -text [= $pid]: ] \
            -row [expr $i+2] -column 0 -pady 2 -sticky nw -in $framePath
        
        # Añadir al diccionario el path del combo
        if {$id in [dict keys $::Itemsdict]} {
            set fullname [dict get $::Itemsdict $id]
            set lista {}
            lappend lista $fullname
            lappend lista $framePath.cmb$id
            dict set ::Itemsdict $id $lista
        }
        
        # Obtenemos la lista de valores para el combo si existe
        
        set icomboList ""
        
        set comboList [::xmlutils::getXMLValues $fullname "" "" "" "" "mat"]
        
        set CBState [::xmlutils::setXml $fullname CBState "read" "" "mat"]
        
        if { $CBState == "normal" } {
            set values $comboList
            set comboList {}
        } else {
            set values {}
        }
        
        if { [llength $comboList] > 0 } {
            
            if { [winfo exists $framePath.cmb$id ] } {
                destroy $framePath.cmb$id
            }
            grid [ttk::combobox $framePath.cmb$id -values $comboList -width [::KMProps::getCmbWidth $comboList] -textvariable "::KMat::cmb$id"] \
                -row [expr $i+2] -column 1 -padx 3 -pady 2 -sticky nw -in $framePath
            tooltip::tooltip $framePath.cmb$id [= $help]
            
            ::xmlutils::setComboDv $framePath.cmb$id $fullname $value
            
        } else { 
            if { [winfo exists $framePath.cmb$id ] } {
                destroy $framePath.cmb$id
            }
            grid [ttk::combobox $framePath.cmb$id -values $values -textvariable "::KMat::cmb$id" -width [::KMProps::getCmbWidth $comboList]] \
                -row [expr $i+2] -column 1 -padx 3 -pady 2 -sticky nw -in $framePath
            tooltip::tooltip $framePath.cmb$id [= $help]
        }
        set ::KMat::cmb$id $value
        
    }
    # Se han añadido los label y los combos
    update idletasks
    return [expr $row +1]
}

proc ::KMat::Combobox_Escape { item T } {
    ::KMat::cmbCancel $item $T
}

proc ::KMat::Combobox_Return { item T } {
    ::KMat::cmbSelectChange $item $T 1 
}

proc ::KMat::Combobox_FocusOut { item T } {
    ::KMat::cmbSelectChange $item $T 0
}

proc ::KMat::buildFrame { T item } {
    
    global KPriv
    set fullname [DecodeName [$T item tag names $item]]
    set idFull [string map { "." "" "//" ""} $fullname]                
    #Comprobamos que sea un item
    
    if { [::KMProps::itemType $fullname] == "p" } {
        set id [::KMat::setXml $fullname id "" ]
        set pid [::KMat::setXml $fullname pid "" ]
        set unit [::KMat::setXml $fullname unit "" ]
        set icon [::KMat::setXml $fullname icon "" ]
        set state [::KMat::setXml $fullname state "" ]
        set help [::KMat::setXml $fullname help "" ]
        
        set xpath "[::KMat::setXPath $fullname]"
        set comboList [::xmlutils::getValues $::KMat::xml $xpath]
        set value [::xmlutils::getValueText $::KMat::xml $xpath "value"]
        
        #
        #---------------------------#---------------------------#
        # Configurar frame en función del XML
        #---------------------------#---------------------------#
        #
        set bg "#F8F8F8"                                        
        set f [frame "$T.f$idFull" -borderwidth 0 -background $bg]
        if { [llength $comboList] > 0 } {
            grid [ttk::combobox $f.cmb -values $comboList -state readonly -textvariable "::KMat::cmb$idFull"] \
                -row 0 -column 0 -padx 3 -sticky ne -in $f
            
            ::xmlutils::setComboValue $::KMat::xml $xpath $f.cmb $value
            #set selected [::xmlutils::getSelected $value $comboList]
            #$f.cmb current $selected
            bind $f.cmb <<ComboboxSelected>> [list KMat::cmbSelectChange $item $T 0]
        } else {
            grid [ttk::combobox $f.cmb -state normal -textvariable ::KMat::cmb$idFull] \
                -row 0 -column 0 -padx 3 -sticky nw -in $f
            set ::KMat::cmb$idFull $value
            bind $f.cmb <FocusOut> [list ::KMat::Combobox_FocusOut $item $T]
            bind $f.cmb <Escape> [list ::KMat::Combobox_Escape $item $T]
        }
        # Si pulsan intro o Esc también forzamos la salida del combo (por probar)
        bind $f.cmb <Return> [list ::KMat::Combobox_Return $item $T]
        
        return $f                  
    } elseif { [::KMProps::itemType $fullname] == "c"} {
        # Falta por implementar ::KMProps::itemType $fullname == m, por si hacen click en un material, o en un grupo de materiales
        if { [::KMat::ContainerType $fullname] == "padre" } {
            ::KMat::CreateBottomFrame $T $item
            return ""
        } elseif { [::KMat::ContainerType $fullname] == "abuelo" } {
            ::KMat::CreateGrandpaFrame $T $item
            return ""
        }
    }
}

proc ::KMat::CreateBottomFrame { T item } {
    # ABSTRACT: Create the botton frame
    
    # If exists the bottom frame destroy it
    set f [::KMat::DestroyBottomFrame]
    
    # Create the frame where set the properties
    ttk::frame $f -borderwidth 0
    # Grid for toolbar
    grid $f -row 2 -column 0 -sticky wes
    ::KMat::buildfatherFrame $T $item
    
    return $f
    
}

proc ::KMat::CreateGrandpaFrame { T item } {
    
    set f [::KMat::DestroyBottomFrame]
    ttk::frame $f -borderwidth 0
    grid $f -row 2 -column 0 -sticky wes
    ::KMat::buildTabFrame $T $item $f
    return $f
}

proc ::KMat::DestroyBottomFrame { } {
    # ABSTRACT: Destroy the botton frame
    variable NbMatsPath
    
    set f ${NbMatsPath}.fBottom    
    if {[winfo exists $f]} {
        foreach w [winfo children $f] {
            destroy $w
        }
        destroy $f
    }
    return $f
}

proc ::KMat::ContainerType { fullname } {
    
    global KPriv
    set path [::xmlutils::setXPath $fullname "mat"]
    set node [$KPriv(xmlMat) selectNodes $path]
    set childs [$node childNodes]
    set contcontrol 0
    foreach child $childs {
        if { [$child hasChildNodes] } {
            set contcontrol 1
            break
        }
    }
    if { $contcontrol == 0 } {
        return "padre"        
    } else {
        return "abuelo"
    }
}

proc ::KMat::buildTabFrameABD { T item {class "Tab"} } {
    
    variable abdlist
    
    set fullname [DecodeName [$T item tag names $item]]
    set id [::KMat::setXml $fullname id "" ]
    set pid [::KMat::setXml $fullname pid "" ]
    
    set f [::KMat::iniFrameBottom]
    
    # Miramos los descendientes directos y si son container ponemos un tab por cada uno q tenga items
    set children [$T item children $item]
    
    set listTabs {}
    set listItems {}
    set acceptItems {}
    
    set nb ${f}.nb
    #        grid [ttk::notebook $nb ] -row 0 -sticky ewn
    grid [ttk::notebook $nb ] -row 0 -column 0 -columnspan 2 -padx 0 -sticky nw -in $f
    
    # declaramos un tab para el material
    set fTab ${nb}.f$id
    $nb add [ttk::labelframe $fTab -text "[= Properties]" -padding {10 0 10 10}] \
        -text "[string range $pid 0 20]"
    
    set count 1
    # foreach itemChild $children {
        #                 set fullname [DecodeName [$T item tag names $itemChild]]
        #                 set nodeName [::xmlutils::getXmlNodeName $fullname "mat"]
        #                 set comboList [::xmlutils::getXMLValues $fullname]
        
        #                 Lappend listItems $itemChild
        
        #                 set id [::KMat::setXml $fullname id ""]
        #                 set pid [::KMat::setXml $fullname pid "" ]
        #                 set value [::KMat::setXml $fullname value "" ]
        #                 set help [::KMat::setXml $fullname help "" ]
        
        #                 if { $id ni $abdlist } {
            #                         grid [ttk::label $fTab.label$id -text "$pid" ] \
            #                                 -row $count -column 0 -pady 5 -sticky nw -in $fTab
            #                         grid [ttk::combobox $fTab.cmb$id -state normal -textvariable "::KMat::cmb$id" ] \
            #                                 -row $count -column 1 -padx 3 -pady 5 -sticky nw -in $fTab
            #                         set ::KMat::cmb$id $value
            #                         incr count
            #                 }
        # }
    
    #$f.bBottomOk and $f.bBottomCancel names need to popup click ok cancel before change item
    grid [ttk::button $f.bBottomOk -text "Ok"  -command "::KMat::acceptTabFrameABD $T $listItems $class $item" ] \
        -row $count -column 0 -sticky sw  -pady 3 -padx 20  -in $f
    tooltip::tooltip $f.bBottomOk [= "Confirm values"]
    
    grid [ttk::button $f.bBottomCancel -text "Cancel"  -command "::KMat::cancelBottom" ] \
        -row $count -column 0 -sticky sw  -pady 3 -padx 100  -in $f
    tooltip::tooltip $f.bBottomCancel [= "Cancel assignation"]
    
}

proc ::KMat::ClickTree { x y T } {
    
    variable lastSelected
    set info [$T identify $x $y]
    
    if { [lindex $info 0] == "item" && [llength $info] >= 4 } {
        set item [lindex $info 1]
        set col [lindex $info 3]
        set fullname [DecodeName [$T item tag names $item]]
        
        #test getMaterialProperties
        #set matlist [::KMat::getMaterialProperties "p" $fullname]
        
        set id [::KMat::setXml $fullname id "" ]
        
        #Eliminamos el anterior combo, si aun está visible
        if {[llength $lastSelected] > 0 && [lindex $lastSelected 0] != $item} { 
            
            ::KMat::cmbSelectChange [lindex $lastSelected 0] $T 1 "anterior"
        }
        
        #Si pinchan en un item con hijos lo abrimos
        #                                $T item toggle $item
        
    } elseif { [lindex $info 0] == "header" && [lindex $info 1] == "0" } {
        if { [$T column cget C0 -arrow] == "up" } {
            $T column configure C0 -arrow down
            $T item sort 0 -dictionary -increasing
        } else {
            $T column configure C0 -arrow up
            $T item sort 0 -dictionary -decreasing
        }
        return ""
    } else {
        return ""
    }
    
    if { $col != 0 } {
        if { ![$T selection includes $item] } {
            $T selection clear
            $T selection add $item
        }
    }
    if { $col == 0 } {
        #SetLayersTo TOUSE $T
    } elseif { $col == 1 } {                
        set parent [winfo parent $T]
    }
    
    if { $col != 0 } {
        return -code break
    }
    return ""
}

proc ::KMat::CreateNewMaterial { {T ""} {name ""} } {
    
    # insertamos en el grupo adecuado. Si no, no insertamos
    set item [$T selection get 0]        
    if {$item == "" } {
        set txt [= "No materials group selected"]
        WarnWin "${txt}."
        return ""
    }
    
    # wa "item:$item"
    
    # Insertamos siempre el nuevo material en el padre del item seleccionado
    set aux 99
    set aux2 [$T item parent $item]
    while {$aux != 1 && $aux2 != 0} {
        set padre [$T item parent $item]
        set item $padre                
        set path [DecodeName [$T item tag names $item]]   
        set splitted [::KEGroups::split2 $path //]
        set aux [llength $splitted]
    }
    
    set path [DecodeName [$T item tag names $item]]   
    
    if { $name == "" } {
        set name [::KMat::GetAutomaticMatName "" $path]
    } else {
        if { ![::KEGroups::isValidGroupName $name] } {
            WarnWin [= "Bad material name, start or end by '//' is not allowed"]
            return ""
        }
    }
    
    ::KMat::insertXml "$path" $name 1 Generic
    
    ::KMat::refreshTree $T
    
    return $name   
}

proc ::KMat::DeleteMaterial { {T ""} {name ""} } {
    
    global KPriv
    
    set items [$T selection get]
    
    if { $items == "" } {
        WarnWin [_ "No material selected."]
        return ""
    }
    
    #Obtenemos el path del nodo
    set path [DecodeName [$T item tag names $items]]   
    
    # si no estamos en el segundo nivel (materiales) no borramos nada
    set splitted [::KEGroups::split2 $path //]                
    
    #                set splitted [::KEGroups::split2 $path //]
    #                set aux [llength $splitted]
    
    if { [llength $splitted] == 1 } {
        WarnWin [_ "No material selected. Can not delete a material group."]
        return ""
    }
    if { [llength $splitted] == 3 } {
        set padre [$T item parent $items]
        set items $padre                
        set path [DecodeName [$T item tag names $items]]   
        set splitted [::KEGroups::split2 $path //]                
    }
    
    if { [llength $splitted] != 2 } {
        WarnWin [_ "Error deleting Material."]
        return ""                
    }
    
    if {[llength $items] > 0 } {
        
        # Para avisar de qué items se van a borrar
        set tuttoItem {}
        foreach it $items {
            lappend tuttoItem "[lindex [$T item text $it] 0]"
        }
        set aviso "Are you sure you want to delete $tuttoItem ?"
        set confirmado [::WinUtils::confirmBox "." "$aviso"]
        if { $confirmado == "ok" } {
            
            #Buscamos todos los descendientes del item a eliminar
            set completList {}
            foreach item $items {                                                
                lappend completList $item
                foreach i [$T item descendants $item] {
                    set idx [lsearch -exact $items $i]
                    if { $idx == -1 } {
                        lappend completList $i                
                    }
                }
            }
            for {set i 0} { $i < [llength $completList] } {incr i} {
                set MaterialId [$T item text [lindex $completList $i] 0]
                set KPriv(materialsId) [::KEGroups::listReplace $KPriv(materialsId) $MaterialId]
            }
            
            for {set i 0} { $i < [llength $items] } {incr i} {
                #Elimina el grupo del xml
                catch { ::xmlutils::unsetXml [DecodeName [$T item tag names [lindex $items $i]]] "mat" }                           
                
                #Resetea todas las ocurrencias del material en las propiedades del .spd
                ::KMProps::checkMaterials [$T item text $item 0]
                
                #Elimina el grupo del árbol
                ::KMProps::deleteItem $T $item
                
            }
        } else {
            #WarnWin [= "No material selected"]
        }
    }
}


proc ::KMat::DeleteTree { {T ""} } {
    variable TreeMatsPath
    
    if { $T == "" } {
        set T $TreeMatsPath
    }
    if { [winfo exists $::KMProps::winpath] } {
        foreach item [$T item range 0 end] {          
            # Elimina el item del árbol
            catch {$T item delete $item}                                                
        }
    }                
}


proc ::KMat::MenuContextualGroup { T x y } {
    
    set w $T.menucontextualgroup
    if { [winfo exists $w] } {
        destroy $w
    }
    menu $w
    
    $w add command -label [= "New Material"] -command [list ::KMat::CreateNewMaterial $T]
    $w add command -label [= "Delete Material"] -command [list ::KMat::DeleteMaterial $T]
    
    set item [$T selection get 0]        
    if {$item != "" } {                
        set path [DecodeName [$T item tag names $item]]   
        set splitted [::KEGroups::split2 $path //]                
        if { [llength $splitted] == 2 } {
            $w add command -label [= "Rename material"] -command [list ::KMat::BeginEditMaterial $T]                                
            $w add command -label [= "Copy material"] -command [list ::KMat::CopyMaterial $T]                                
        } else {
            $w add command -label [= "Rename material"] -state disabled                                
            $w add command -label [= "Copy material"] -state disabled                                
        }
    }
    
    
    #                $w add command -label [= "Collapse All"] -command [list $T collapse -recurse "$item"] -state normal
    #                $w add command -label [= "Expand All"] -command [list $T expand -recurse "$item"] -state normal
    
    set x [expr [winfo rootx $T]+$x+2]
    set y [expr [winfo rooty $T]+$y]
    GiD_PopupMenu $w $x $y
}

proc ::KMat::insertXml { path id state type } {
    global KPriv
    
    # wa "path:$path id:$id state:$state type:$type"
    if { $path == "root" } { 
        set xpath "/Kratos_KMat_DB/Materials"
    } else {
        set xpath "[::KMat::setMatXPath $path]"
    }
    # wa "xpath:$xpath"
    set templatePath "/Kratos_KMat_DB/Templates/Template"
    if { $path == "Metal" } {
        set maticon "grey.gif"
    } elseif { $path == "Fluid" } {
        set maticon "blue.gif"                
    } elseif { $path == "Plastic" } {
        set maticon "red.gif"                
    } elseif { $path == "Composite" } {
        set maticon "green.gif"                
    } elseif { $path == "DEMMaterial" } {
        set maticon "grey.gif"                
    }
    
    set CurrentTemplateId "NewMaterial"
    if {$path == "DEMMaterial"} {
        set CurrentTemplateId "NewDEMMaterial"
    }
    
    set attributesArray [list id=\'$id\' pid=\'$id\' icon=\'$maticon\' help=\'$id\' open=\'1\']
    # wa "attributesArray:$attributesArray"
    ::xmlutils::copyTemplate $::KMat::xml $xpath $templatePath $CurrentTemplateId "Material" $attributesArray 
    
    set xmlArray [::xmlutils::replaceTemplate $::KMat::xml $xpath]
    
    set KPriv(xmlDocMat) [lindex $xmlArray 0]
    set KPriv(xmlMat) [lindex $xmlArray 1]
    
    set ::KMat::xml $KPriv(xmlMat)
}

proc ::KMat::insertXmlCopy { path id state type sourcematname} {
    global KPriv
    
    if { $path == "root" } { 
        set xpath "/Kratos_KMat_DB/Materials"
    } else {
        set xpath "[::KMat::setMatXPath $path]"
    }
    
    set templatePath "$xpath/Material"
    if { $path == "Metal" } {
        set maticon "grey.gif"
    } elseif { $path == "Fluid" } {
        set maticon "blue.gif"                
    } elseif { $path == "Plastic" } {
        set maticon "red.gif"                
    } elseif { $path == "Composite" } {
        set maticon "green.gif"                
    } elseif { $path == "DEMMaterial" } {
        set maticon "grey.gif"                
    }
    
    set attributesArray [list id=\'$id\' pid=\'$id\' icon=\'$maticon\' help=\'$id\' open=\'1\']
    
    # wa "sourcematname:$sourcematname"
    ::xmlutils::copyTemplate  $::KMat::xml $xpath $templatePath "$sourcematname" "Material" $attributesArray 
    
    set xmlArray [::xmlutils::replaceTemplate $::KMat::xml $xpath]
    
    set KPriv(xmlDocMat) [lindex $xmlArray 0]
    set KPriv(xmlMat) [lindex $xmlArray 1]
    
    set ::KMat::xml $KPriv(xmlMat)
}


# Prepara la query para utilizar las funciones de domNOde
proc ::KMat::setMatXPath { path } {
    set splitted [::KEGroups::split2 $path //]        
    # wa "splitted:$splitted path:$path"        
    if { [llength $splitted] >= 1 } {                
        set xpath "/Kratos_KMat_DB/Materials/MaterialGroup\[@id='[lindex $splitted 0]'\]"
    } 
    if { [llength $splitted] >= 2 } {                
        set xpath "$xpath/Group\[@id='[lindex $splitted 1]'\]"
    } 
    if { [llength $splitted] >= 3 } {                
        set xpath "$xpath/Group\[@id='[lindex $splitted 2]'\]"
    }
    if { [llength $splitted] >= 4 } {                
        set xpath "$xpath/Group\[@id='[lindex $splitted 3]'\]"
    }
    if { [llength $splitted] >= 5 } {                
        set xpath "$xpath/Group\[@id='[lindex $splitted 4]'\]"
    }
    return $xpath
}



proc ::KMat::InsertNewMaterial { MatName T {state 1} {type "Generic"} {parent ""} {parentitem root} {childs true} } {
    
    if { $parent != "" } {
        set fullname $parent$MatName
    } else {
        set fullname $MatName
    }
    
    set item [$T item create -button yes -tags [EncodeName $fullname]]
    $T item lastchild $parentitem $item
    $T item style set $item C0 styAnyRead
    $T item element configure $item C0 elemTxtRead -text "$MatName"
    
    set item2 [::KMat::InsertNewItem [= "Density"] "p.Density" $T "$fullname" "$item" 0]        
    set item3 [::KMat::InsertNewItem [= "Young_Modulus"] "p.Young_Modulus" $T "$fullname" "$item" 0]        
    set item4 [::KMat::InsertNewItem [= "Poisson_Ratio"] "p.Poisson_Ratio" $T "$fullname" "$item" 0]        
    set item5 [::KMat::InsertNewItem [= "Viscosity"] "p.Viscosity" $T "$fullname" "$item" 0]        
    set item6 [::KMat::InsertNewItem [= "Description"] "p.Description" $T "$fullname" "$item" 0]        
    
    set item7 [::KMat::InsertNewItem [= "Description2"] "c.General//p.Description2" $T "$fullname" "$item" 0]
    set item8 [::KMat::InsertNewItem [= "Description"] "p.Description" $T "$fullname" "$item" 0]
    set item8 [::KMat::InsertNewItem [= "Density"] "p.Density" $T "$fullname" "$item" 0]
    
    ::KMat::refreshTree $T
    return $item
}


proc ::KMat::BeginEditMaterial { T } {
    set I [$T selection get 0]
    set C 0
    set E elemTxtRead
    
    set path [DecodeName [$T item tag names $I]]   
    set splitted [::KEGroups::split2 $path //]                
    if { [llength $splitted] == 2 } {
        ::TreeCtrl::FileListEdit $T $I $C $E
    }
    
}

proc ::KMat::CopyMaterial { {T ""} {name ""} } {
    
    # insertamos en el grup adecuado. Si no, no insertamos
    set item [$T selection get 0]        
    if {$item == "" } {
        set txt [= "No materials group selected"]
        WarnWin "${txt}."
        return ""
    }
    
    set fullname [DecodeName [$T item tag names $item]]
    set sourcematname [::KMat::setXml $fullname id "" ]
    
    #Insertamos siempre el nuevo material en el padre del item seleccionado
    set aux 99
    set aux2 [$T item parent $item]
    while {$aux != 1 && $aux2 != 0} {
        set padre [$T item parent $item]
        set item $padre                
        set path [DecodeName [$T item tag names $item]]   
        set splitted [::KEGroups::split2 $path //]
        set aux [llength $splitted]
    }
    
    set path [DecodeName [$T item tag names $item]]   
    
    if { $name == "" } {
        set name [::KMat::GetAutomaticMatName "" $sourcematname]
    } else {
        if { ![::KEGroups::isValidGroupName $name] } {
            msg [= "Bad material name, start or end by '//' is not allowed"]
            return ""
        }
    }
    
    ::KMat::insertXmlCopy "$path" $name 1 Generic $sourcematname
    
    
    ::KMat::refreshTree $T
    
    return $name   
    
}


proc ::KMat::SetMatToRename { T item newtext } {
    
    global KPriv
    if { $newtext == ""} {  return  }
    if { $item == 0 } {
        msg [= "Root folder can't be edited"]
        return
    }
    set oldId [$T item text $item 0]
    #Controlamos q el nombre no esté ya en el árbol (a no ser q no se haya cambiado)
    if { $oldId != $newtext && $newtext in $KPriv(materialsId) } {
        msg "The group name '$s' already exist.\nChoose another, please." $newtext
        return
    }
    #Validamos q el nombre no tenga carácteres que vulneran la seguridad y quitamos espacios
    set newtext [::KUtils::parseTreeStr $newtext]
    if { $newtext == -1 } {
        msg "You can't use some reservate chars like:\n  :   /   $   .   \\  %  "
        return
    }
    
    #Guardamos el nivel del item a renombrar
    set fullname [DecodeName [$T item tag names $item]]
    set splitted [::KEGroups::split2 $fullname //]
    set whereRename [llength $splitted]
    
    #Recorremos toda su familia                cambiando cada path
    set items [$T item descendants $item]
    
    foreach i $items {                
        set fullname [DecodeName [$T item tag names $i]]
        set splitted [::KEGroups::split2 $fullname //]
        
        #Sustituimos la posición apropiada del path por el nuevo nombre
        set splitted [::KEGroups::listReplace $splitted [lindex $splitted [expr $whereRename - 1]] m.[list $newtext]]
        
        #Reconstruimos el nuevo path
        set fullNewName ""
        foreach iSplit $splitted {                
            set fullNewName "$fullNewName$iSplit//"
        }
        set fullNewName [string range $fullNewName 0 end-2]
        
        #Cambiamos la etiqueta (el path) en el arbol
        $T item tag remove $i [list names [$T item tag names $i]]
        $T item tag add $i [EncodeName $fullNewName]
    }
    
    #Cambiamos el item a renombrar 
    set fullname [DecodeName [$T item tag names $item]]
    
    ::KMat::editTag $T $item $fullname $newtext
    
    #Renombramos también todas las ocurrencias del material en el .spd de propiedades
    ::KMProps::checkMaterials $oldId $newtext
    
    return ""
}


#
# Renombra un item, modificando su path y reconstruyendo el combo
#
proc ::KMat::editTag { T item fullname newtext } {                
    global KPriv                
    
    set parts [::KEGroups::split2 $fullname //]
    lset parts end m.[list $newtext]
    set newPath [join $parts //]
    
    # Renombra en la lista de ID's
    set idItem [lindex [$T item text $item] 0]
    set KPriv(materialsId) [::KEGroups::listReplace $KPriv(materialsId) $idItem m.[list $newtext]]
    
    # Cambiar nombre en el árbol
    $T item tag remove $item [list names [$T item tag names $item]]
    $T item tag add $item [EncodeName $newPath]
    $T item element configure $item C0 elemTxtRead -text $newtext
    
    # Cambiar nombre en el XML
    ::KMat::setXml $fullname pid $newtext 
    ::KMat::setXml $fullname id  $newtext
    
    set createframeafteredit 0
    if {$createframeafteredit} {
        set childs [$T item children $item]
        foreach child $childs {                
            set fullname [DecodeName [$T item tag names $child]]
            set idFull [string map { "." "" "//" ""} $fullname]
            destroy "$T.f$idFull"                
            set f [::KMat::buildFrame $T $child]                                        
        }
    }
    
    return $newPath
}

proc ::KMat::GetAutomaticMatName { {auto ""} { startname "" } } {                
    global KPriv
    
    set name ""
    
    set i 0
    foreach grup $KPriv(materialsId) {
        incr $i
    }
    if { [llength $KPriv(materialsId)] > 0 } {
        for {set i 1} {$i<10000} {incr i} {
            # set name ${auto}Material${i}
            set name ${auto}$startname${i}
            if { [lsearch -exact $KPriv(materialsId) $name] == -1 } { break }
        }
    } else {
        if { $auto == "" } {
            set name "Material1"
        } else {
            set name "${auto}Material1"
        }
    }
    return $name
}

#
# Separa cada node en "inicialNombreTag.idNodo"
#
proc ::KMat::splitNode { node } {
    
    set id [$node getAttribute id ""]                
    if { [$node tagName] == "Container"} {
        return "c.[list $id]"
    } elseif { [$node tagName] == "Item"} {
        return "i.[list $id]"
    } elseif { [$node tagName] == "Property"} {
        return "p.[list $id]"
    } elseif { [$node tagName] == "Material"} {
        return "m.[list $id]"
    } else {
        return "NoTree"
    }
}

proc ::KMat::InsertNewItem { propName id T {parent ""} {parentitem root} {childs true} {state "normal"} {open 0}} {
    set propName [= $propName]
    # wa "propName:$propName id:$id T:$T parent:$parent parentitem:$parentitem childs:$childs state:$state open:$open"
    if { $state == "hidden" } {
        #No inserta este item pero se sigue el proceso
        return $parentitem
    } elseif { $state == -1 } {
        #No inserta este item ni su descendencia
        return -1
    } 
    if {[string first "NoTree" $id] != -1 } {
        #Solo insertamos determinados nodos (container, item..)
        return "-1"
    }
    if { $parent != "" } {
        # set fullname "$parent//$id"
        set fullname "$parent$id"
    } else {
        set fullname $id
    }
    # wa "fullname:$fullname"
    set xpath "[::KMat::setXPath $fullname]"
    set tooltip [::xmlutils::getValueText $::KMat::xml $xpath "help"]
    
    if { $childs } {                
        set item [$T item create -button yes -tags [EncodeName $fullname] -open $open]
        if { $::KMat::UseFullTkTree} {
            $T popup_enter_help $item [= $tooltip]
        }
        $T item lastchild $parentitem $item
        $T item style set $item C0 styAnyRead
        $T item element configure $item C0 elemTxtRead -text "$propName"                                                
    } else {
        set item [$T item create -button no -tags [EncodeName $fullname] -open $open]
        if { $::KMat::UseFullTkTree} {
            # Looks like this is the one that creates the annoying tooltip
            $T popup_enter_help $item [= $tooltip]
        }
        $T item lastchild $parentitem $item
        $T item style set $item C0 styFrame
        
        set xpath "[::KMat::setXPath $fullname]"
        # wa "xpath:$xpath"
        set value [::xmlutils::getValueText $::KMat::xml $xpath "value"]
        
        $T item style set $item C0 styAnyRead
        if { $parentitem == "root" } {
            $T item element configure $item C0 elemTxtRead -text "$propName $value"
        } else {
            $T item element configure $item C0 elemTxtRead -text "$propName: $value"
        }
    }
    
    # Consultamos el icono en el xml, y si existe en nuestro directorio se lo añadimos                
    set icon [::KMat::setXml $fullname icon "" ]
    set imagen [::WinUtils::GetImage $icon]
    if { $imagen != -1 } {
        $T item image $item C0 $imagen
    } else {
    }        
    
    # Miramos si el item tiene que estar a disabled (viene del proc ::KMProps::stateNode)
    if {$state == "disabled"} {
        $T item enabled $item 0
    }        
    return $item
}




proc ::KMat::cmbSelectChange { item T {remove 1} {selectVal current}  } {   
    
    set fullname [DecodeName [$T item tag names $item]]
    set idFull [string map { "." "" "//" ""} $fullname]
    
    set id [::KMat::setXml $fullname id "" ]
    set pid [::KMat::setXml $fullname pid "" ]
    
    set xpath "[::KMat::setXPath $fullname]"
    set value [::xmlutils::getValueText $::KMat::xml $xpath "value"]
    
    #Antes solo se hacía esto:
    #set selCombo [set ::KMat::cmb$idFull]                        
    
    #Gestión de ivalues
    set xpath "[::KMat::setXPath $fullname]"
    set comboState [::xmlutils::getComboState $::KMat::xml $xpath]
    if { $comboState == "normal" } {
        
        set selCombo [set ::KMat::cmb$idFull]
        set selComboText $selCombo
        if {$selectVal != "current"} {
            set ::KMat::lastSelected {}
        }
    } else {
        
        if {$selectVal == "current"} {
            
            set f "$T.f${idFull}.cmb"
            set selCombo [::xmlutils::getComboValue $::KMat::xml $xpath $f]
            set selComboText [::xmlutils::getComboValue $::KMat::xml $xpath $f "text"]
            set ::KMat::lastSelected [list $item $selCombo $selComboText]
            
        } else {
            
            set selCombo [lindex $::KFun::lastSelected 1]
            set selComboText [lindex $::KFun::lastSelected 2]
            set ::KMat::lastSelected {}
        }
    }
    
    if { $remove } {
        
        set f "$T.f$idFull"
        if {[winfo exists $f]} {
            destroy $f
        }
        $T item style set $item C0 styAnyRead
        $T item element configure $item C0 elemTxtRead -text "$pid: $selComboText"
    }
    
    #Guarda el nuevo valor en el xml
    #  "selcombo = $selCombo"
    #msg "fullname = $fullname"
    ::KMat::setXml $fullname value $selCombo
    
    ::KMat::insertIcon $item $T $fullname
    
    #Volvemos a cargar el árbol de propiedades
    if { !$remove } {
        #::KMat::refreshTree $T
    }
    ::KMat::refreshTree $T
}

proc ::KMat::cmbCancel { item T  } {
    
    set fullname [DecodeName [$T item tag names $item]]
    set idFull [string map { "." "" "//" ""} $fullname]
    
    set pid [::KMat::setXml $fullname pid ""]
    set id [::KMat::setXml $fullname id ""]
    set value [::KMat::setXml $fullname value ""]
    
    set xpath "[::KMat::setXPath $fullname]"
    set valueText [::xmlutils::getValueText $::KMat::xml $xpath "value"]
    
    set ::KMat::cmb$idFull "$value"
    set f "$T.f$idFull"
    if {[winfo exists $f]} {
        destroy $f
    }
    $T item style set $item C0 styAnyRead
    $T item element configure $item C0 elemTxtRead -text "$pid: $valueText"
    
    ::KMat::insertIcon $item $T $fullname
    set ::KMat::lastSelected ""
    ::KMat::refreshTree
}

proc ::KMat::insertIcon { item T fullname {type "mat"}} {
    
    # Consultamos el icono en el xml, y si existe en nuestro directorio se lo añadimos                
    if {$type == "mat"} {
        set icon [::KMat::setXml $fullname icon ""]
    } else {
        set icon [::KProps::setXml $fullname icon]
    }
    
    set imagen [::WinUtils::GetImage $icon]
    if { $imagen != -1 } {
        $T item image $item C0 $imagen
    }                        
}


proc ::KMat::refreshTree { {T ""} } {
    variable TreeMatsPath; variable lastSelected
    
    if {$T == ""} {
        set T $TreeMatsPath
    }
    
    #Primero hay que asegurarse de que exista la ventana
    if { [winfo exists $::KMProps::WinPath] } {
        
        foreach item [$T item range 0 end] {          
            
            set fullname [DecodeName [$T item tag names $item]]                                                
            if {$fullname != "" } {                
                catch {
                    ::KMat::setXml $fullname open [$T item isopen $item]
                }
            }
        }
    }
    ::KMat::DeleteTree                
    ::KMat::FillTreeMat
    
    set lastSelected {}
}

proc ::KMat::findMaterialParent { matid } {
    
    global KPriv
    
    set nodes [$::KMat::xml selectNodes "/Kratos_KMat_DB/Materials/MaterialGroup\[@id\]"]
    
    foreach node $nodes {                          
        set nodes2 [$node childNodes]
        foreach node2 $nodes2 {  
            set aux [$node2 getAttribute id ""]
            if { $aux == $matid} {
                set parent [$node getAttribute id ""]
                return $parent
            }
        }                
    }
}

proc ::KMat::getMaterials {{application ""}} {
    variable TreeMatsPath
    global KPriv
    
    set KPriv(materialsList) {}
    set T $TreeMatsPath
    
    set nodes ""
    if { $::KMat::xml !="" } {
        set nodes [$::KMat::xml selectNodes "/Kratos_KMat_DB/Materials/MaterialGroup\[@id\]"]
    }
    
    if {$application == ""} {
        foreach node $nodes {                
            set nodes2 [$node childNodes]
            foreach node2 $nodes2 {  
                set aux [$node2 getAttribute id ""]
                lappend KPriv(materialsList) $aux                                
            }                
        }
    } else {
        if {$application == "StructuralAnalysis"} {
            set smatlist [list "Metal" "Composite" "Plastic"]
            foreach node $nodes {                
                if {[$node getAttribute id ""] in $smatlist} {
                    set nodes2 [$node childNodes]
                    foreach node2 $nodes2 {  
                        set aux [$node2 getAttribute id ""]
                        lappend KPriv(materialsList) $aux                                
                    }                
                }
            }
        } elseif {($application == "Fluid")||($application == "PFEM")} {
            foreach node $nodes {                
                if { [$node getAttribute id ""] == "Fluid"} {
                    set nodes2 [$node childNodes]
                    foreach node2 $nodes2 {  
                        set aux [$node2 getAttribute id ""]
                        lappend KPriv(materialsList) $aux                                
                    }                
                }
            }
        } elseif {$application == "FluidStructureInteraction"} {
            foreach node $nodes {                
                if { [$node getAttribute id ""] == "Metal" || [$node getAttribute id ""] == "Fluid"} {
                    set nodes2 [$node childNodes]
                    foreach node2 $nodes2 {  
                        set aux [$node2 getAttribute id ""]
                        lappend KPriv(materialsList) $aux                                
                    }                
                }
            }
        } elseif {$application == "ConvectionDiffusion"} {
            foreach node $nodes {                
                if { [$node getAttribute id ""] == "Metal" || [$node getAttribute id ""] == "Composite"} {
                    set nodes2 [$node childNodes]
                    foreach node2 $nodes2 {  
                        set aux [$node2 getAttribute id ""]
                        lappend KPriv(materialsList) $aux                                
                    }                
                }
            }
        } elseif {$application == "DEM"} {
            foreach node $nodes {                
                if { [$node getAttribute id ""] == "DEMMaterial"} {
                    set nodes2 [$node childNodes]
                    foreach node2 $nodes2 {  
                        set aux [$node2 getAttribute id ""]
                        lappend KPriv(materialsList) $aux                                
                    }                
                }
            }
        } elseif {$application == "DSOLID"} {
            foreach node $nodes {
                if { [$node getAttribute id ""] == "SolidMaterial" } {
                    set nodes2 [$node childNodes]
                    foreach node2 $nodes2 {
                        set aux [$node2 getAttribute id ""]
                        lappend KPriv(materialsList) $aux
                    }
                }
            }
        } else {
            foreach node $nodes {                
                set nodes2 [$node childNodes]
                foreach node2 $nodes2 {  
                    set aux [$node2 getAttribute id ""]
                    lappend KPriv(materialsList) $aux                                
                }                
            }                                        
        }
    }
    
    return $KPriv(materialsList)
}

proc ::KMat::getMaterialProperty { matid prop } {
    
    global KPriv
    set KPriv(materialsList) {}
    
    set xpath "/Kratos_KMat_DB/Materials/MaterialGroup/Material\[@id=\"$matid\"]/Property\[@id=\"$prop\"]"
    set node [$::KMat::xml selectNodes $xpath]
    
    set value [$node @value]
    return $value
}

proc ::KMat::getMaterialProperties { pathtype path } {
    
    # pathtype = c if you pass a container path
    # pathtype = p if you pass a property path
    global KPriv
    if { $pathtype == "c" } {
        set xpath "[::KMat::setXPath $path]/Property\[@id\]"
    } elseif { $pathtype == "p" } {
        set xpath "[::KMat::setXPath $path]"
    }
    
    set matPropsList {}
    
    set nodes [$::KMat::xml selectNodes "$xpath"]
    
    foreach node $nodes {                          
        set aux "[$node getAttribute id ""] [$node @value]"
        lappend matPropsList $aux                                
    }
    
    return $matPropsList
    
}


proc ::KMat::testbutton { } {
    
    set matlist [::KMat::getMaterialProperties]
}


#
# Editar o extraer propiedades del xml en memoria
#
proc ::KMat::setXml { path property {value ""} } {
    
    global KPriv
    
    set xpath "[::KMat::setXPath $path]"
    
    if { $value == "" } {
        set value [$::KMat::xml set "$xpath/@$property" ]                                
        
        
        #Cuando hay espacios el xml devuelve una lista y si la imprimes tal cual aparecen corchetes
        if { [llength $value] == 1 } {
            set value [lindex $value 0]
        }
        return $value
    } else {
        $::KMat::xml set "$xpath/@$property" "$value"
        return "1"
    }
}

#
# Crea el frame inferior 
#
proc ::KMat::iniFrameBottom { } {
    
    #Destruye el frame inferior si existía
    set f [::KMat::cancelBottom]
    
    # Create the frame where set the properties
    ttk::frame $f -borderwidth 0
    # Grid for toolbar
    grid $f -row 2 -column 0 -sticky wes
    
    return $f
}

#
# Destruye el frame inferior 
#
proc ::KMat::cancelBottom { } {
    variable NbMatsPath
    set f ${NbMatsPath}.fBottom
    if { [winfo exists $f]} {
        foreach w [winfo children $f] {
            destroy $w
        }
    }
    destroy $f
    
    return $f
}




proc ::KMat::acceptTabFrameABD { T listItems class {itemSel ""}} {
    WarnWin [_ "Accept Tab Frame ABD, save values in XML"]
}


proc ::KMat::acceptTabFrame { T listItems class {itemSel ""} } { 
    
    # Boton OK debajo de los Notebookç
    if { $class == "Tab" } {
        foreach {item lista} $::Itemsdict {
            
            set fullname [lindex $lista 0]
            set combopath [lindex $lista 1]
            
            # Leer el valor del combo
            set valor [$combopath get]
            set fullname [lindex $fullname 0]
            # Escribir en el combo
            ::KMat::setXml $fullname value $valor
        }
        
    } else {
        set fTab [lindex $class 0]
        set listItems [lindex $class 1]
        
        foreach item $listItems {
            set fullname [DecodeName [$T item tag names $item]]
            set id [::KMat::setXml $fullname id "" ]
            set combopath $fTab.cmb$id
            set valor [$combopath get]
            ::KMat::setXml $fullname value $valor
        }
    }
    #Volvemos a cargar el árbol para q el path de los items sea correcto
    ::KMat::cancelBottom
    ::KMat::refreshTree $T
}

proc ::KMat::initVisibilityClass { } {
    
    variable visibilityVars
    set visibilityVars {}
    global KPriv
    
    set classes ""
    if { $::KMat::xml !=""}  {
        set classes [$::KMat::xml set "/ Kratos_KMat_DB/Materials/ClassConfiguration/Class" ]
    }
    set classes [split $classes ","]
    
    foreach class $classes {
        lappend visibilityVars $class
        set ::KMat::$class ""
    }
}


#
# Valida varias cosas para cada nodo
#
proc ::KMat::stateNode { node } {
    
    #Validamos para cada nodo si tiene que estar visible 
    #(en función de los valores elegidos en algunos combos)
    # wa "nodeName:[$node nodeName]"
    if { [$node nodeName] == "Property" } {
        
        # Leemos la class del nodo para ver si requiere de acciones especiales (ocultar nodos)
        set value [$node getAttribute value ""]
        set class [$node getAttribute class ""]
        
        # wa "value:$value class:$class"
        # Equivalente a Switch $class
        foreach var $::KMat::visibilityVars {        
            if {$var == $class} {                                
                #Caso especial para el solver de fluidos
                if { $var == "fluidSolvTyp" } {   
                    #El solver de fluidos está duplicado dependiendo de una variable prebia
                    set freeYesOrNo [$node getAttribute freeSurf ""]
                    #msg "FREEEEEEEEEEEEEE:  $freeYesOrNo == \"\" || $freeYesOrNo == $::KMProps::freeSurf"
                    if { $freeYesOrNo == "" || $freeYesOrNo == $::KMProps::freeSurf } {
                        set ::KMProps::fluidSolvTyp "$value"
                    }
                } else {
                    # Caso general
                    set ::KMat::$var $value
                }
            }
        }
    } else {
        
        # Caso especial para application
        set class [$node getAttribute class ""]
        if { $class == "application" } {                        
            set apliState [$node getAttribute state ""]
            #msg "apliState$apliState"
            if {$apliState != "hiddenAll" } {
                set ::KMProps::application [$node getAttribute id ""]
            }
        }
    }
    
    set state [$node getAttribute state "normal"]
    # wa "\n[$node getAttribute id ""]                                state:$state"
    
    #Si el estado es hiddenAll se oculta el nodo y toda su descendencia
    if {$state == "hiddenAll"} {
        return "-1"
    }
    
    foreach var $::KMat::visibilityVars {
        set globalVar [set ::KMat::$var]                
        set nodeValuesVar [split [$node getAttribute $var ""] ","]                
        #Si el nodo tiene alguna restriccion de clase (p.ej. del tipo strucType=Shell)
        # y no coincide con el valor seleccionado, ocultamos el nodo
        if { $nodeValuesVar != "" && !($globalVar in $nodeValuesVar) } {                        
            if {$var == "strucType" } {
                if { $globalVar != "Generic"} { 
                    #msg "$nodeValuesVar \"\" !=  && $globalVar in $nodeValuesVar"
                    return -1
                }
            } else {
                #msg "$nodeValuesVar \"\" !=  && $globalVar in $nodeValuesVar"
                return -1
            }
        }
    }
    # Devolvemos el estado del nodo ("normal" por defecto)
    return $state        
}

#
#
#
# Prepara la query para utilizar las funciones de domNOde
proc ::KMat::setXPath { path } {
    
    set splitted [::KMProps::split2 $path //]
    
    set i 0
    foreach itemId $splitted {
        if { $i == 0 } {                                
            #set xpath "/Kratos_KMat_DB/Materials/Material\[@id='$itemId'\]"                  
            set xpath "/Kratos_KMat_DB/Materials/MaterialGroup\[@id='$itemId'\]"                 
            
        } else {
            if { [string index $itemId 0] == "p" } {
                set xpath "$xpath/Property\[@id='[string range $itemId 2 end]'\]"
            }
            if { [string index $itemId 0] == "m" } {
                set xpath "$xpath/Material\[@id='[string range $itemId 2 end]'\]"
            }
            if { [string index $itemId 0] == "c" } {
                set xpath "$xpath/Container\[@id='[string range $itemId 2 end]'\]"
            }
        }                
        incr i
    }
    return $xpath
}

proc ::KMat::FromNodetoItem { node } {
    
    set path [::xmlutils::getPathFromNode $node 1 mats]
    #msg "path $path"
    set item [::KMat::FromPathtoItem $path] 
    
    #msg "item = $item"
    return $item
}

proc ::KMat::FromPathtoItem { path } {
    variable TreeMatsPath
    
    set Item 1
    set splitted [::KMProps::split2 $path / ]
    #        msg "splitted $splitted"
    
    set materialgroup [lindex $splitted 3]
    #        msg "materialgroup $materialgroup"
    set splitted2 [::KMProps::split2 $path ' ]
    #        msg "splitted2 $splitted2"
    
    set compid [lindex $splitted2 1]
    #        msg "compid $compid"
    set T $TreeMatsPath
    
    set treenodes [$T item range 1 end]
    foreach item $treenodes {
        #                msg "item = $item"
        set fullname [DecodeName [$T item tag names $item]]
        #msg "fullname = $fullname"
        set id [::xmlutils::setXml $fullname id "read" "" mat]
        #                msg "id $id"
        if {$id == $compid } {
            set Item $item
            break
        }
    }
    
    set compid [lindex $splitted2 3]
    #        msg "compid $compid"
    set T $TreeMatsPath
    
    set treenodes [$T item range $Item end]
    foreach item $treenodes {
        #                msg "item = $item"
        set fullname [DecodeName [$T item tag names $item]]
        #msg "fullname = $fullname"
        set id [::xmlutils::setXml $fullname id "read" "" mat]
        #                msg "id $id"
        if {$id == $compid } {
            set Item $item
            break
        }
    }
    set splitted2 [lrange $splitted2 5 end]
    #        msg "splitted2 $splitted2"
    
    foreach {cont rubish} $splitted2 {
        set compid [lindex $splitted2 0]
        #                msg "compid $compid"
        set T $TreeMatsPath
        
        set treenodes [$T item range $Item end]
        foreach item $treenodes {
            #                        msg "item = $item"
            set fullname [DecodeName [$T item tag names $item]]
            #msg "fullname = $fullname"
            set id [::xmlutils::setXml $fullname id "read" "" mat]
            #                        msg "id $id"
            if {$id == $compid } {
                set Item $item
                break
            }
        }
        set splitted2 [lrange $splitted2 2 end]
    }
    return $Item
}

proc ::KMat::NewcmbSelectChange { item combopos {remove 1} {selectVal current} } {   
    variable TreeMatsPath
    set T $TreeMatsPath
    set fullname [DecodeName [$T item tag names $item]]
    set idFull [string map { "." "" "//" ""} $fullname]
    
    set id [::KMat::setXml $fullname id "" ]
    set pid [::KMat::setXml $fullname pid "" ]
    
    set xpath "[::KMat::setXPath $fullname]"
    set value [::xmlutils::getValueText $::KMat::xml $xpath "value"]
    
    #Antes solo se hacía esto:
    #set selCombo [set ::KMat::cmb$idFull]                        
    
    #Gestión de ivalues
    set xpath "[::KMat::setXPath $fullname]"
    set comboState [::xmlutils::getComboState $::KMat::xml $xpath]
    if { $comboState == "normal" } {
        
        set selCombo [set ::KMat::cmb$idFull]
        set selComboText $selCombo
        if {$selectVal != "current"} {
            set ::KMat::lastSelected {}
        }
    } else {
        
        if {$selectVal == "current"} {
            
            set f "$T.f${idFull}.cmb"
            set selCombo [::xmlutils::getComboValue $::KMat::xml $xpath $f]
            set selComboText [::xmlutils::getComboValue $::KMat::xml $xpath $f "text"]
            set ::KMat::lastSelected [list $item $selCombo $selComboText]
            
        } else {
            
            set selCombo [lindex $::KFun::lastSelected 1]
            set selComboText [lindex $::KFun::lastSelected 2]
            set ::KMat::lastSelected {}
        }
    }
    
    if { $remove } {
        
        set f "$T.f$idFull"
        if {[winfo exists $f]} {
            destroy $f
        }
        $T item style set $item C0 styAnyRead
        $T item element configure $item C0 elemTxtRead -text "$pid: $selComboText"
    }
    
    #Guarda el nuevo valor en el xml
    #  "selcombo = $selCombo"
    #msg "fullname = $fullname"
    ::KMat::setXml $fullname value $selCombo
    
    ::KMat::insertIcon $item $T $fullname
    
    #Volvemos a cargar el árbol de propiedades
    if { !$remove } {
        #::KMat::refreshTree $T
    }
    ::KMat::refreshTree $T
}

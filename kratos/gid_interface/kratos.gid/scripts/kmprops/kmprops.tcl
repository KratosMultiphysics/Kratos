###############################################################################
#
#  NAME: kmprop.tcl
#
#  PURPOSE: Main window to manage model properties
#
#  QUANTECH ATZ-DEVELOPMENT DEPARTMENT
#
#  AUTHOR : L. Calvo
#
#  CREATED AT: 25/02/2010
#
#  LAST MODIFICATION : 
#
#  VERSION : 0.1
#
#  HISTORY:
#   0.4- 05/05/10 LCA: Se han habilitado combos con idioma añadiendo al spd la propiedad "ivalues"
#   0.3- 26/04/10 LCA: Estructura correcta del árbol, con filtros por dimensión(2D,3D) 
#                y por Types (Structural, Solution y Analysis). Funciona link entre propiedades y materiales.
#   0.2- 29/03/10 LCA: Hasta este momento tenemos funcionando el tab Properties
#                (el de Materials está en proceso). Se pueden editar items y asignar grupos a condiciones.
#   0.1- 25/02/10 LCA, create a base source code from the kegroups.tcl script
#
###############################################################################

package require treectrl
package require tooltip
package provide KMprops 1.0 
package require snit
package require tdom
package require xmlstruct 1.0

# Create a base namespace KEGroups
namespace eval ::KMProps:: {
	
	# Group properties array
	variable Props
	# Path of the base window 
	variable WinPath ".gid.kmprops"
	# Window layout ["OUTSIDE"|"LEFT"|"RIGHT"] 
	variable WinLayout 
	variable SystemHighlight
	variable SystemHighlightText
	variable TreePropsPath
	variable ngroups
	variable lastSelected {}
	variable selectedEntity ""
	
	variable application "StrucutalAnalysis"
	
	#Se inicializan las clases dinámicamente leyendo del xml
	variable visibilityVars {}
	#variable visibilityVars {nDim strucType soluType analysType linearSolvTyp solverType fluidType freeSurf}
}

proc ::KMProps::Pruebas { } {
	
	#::KMValid::ValidateModel
	#msg "f[::xmlutils::getKKWord StructuralAnalysis TotalLagrangian2]f"
	return 0
}

proc ::KMProps::Init {} {
	
	variable WinLayout
	variable SystemHighlight
	variable SystemHighlightText
	variable ngroups
	variable Props
	
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

	set WinnLayout "OUTSIDE"
	set ngroups 10000
}

proc ::KMProps::InitBaseWindow {{whattab "Model"} {what "OUTSIDE"} } {
	
	#PRUEBAS - Fución para hacer pruebas. Si devuelve True se no carga la ventana
	if {[::KMProps::Pruebas]} { return 0 }
	
	variable WinLayout

	# Init KEGroups namaspace variables
	::KMProps::Init
	::KMat::Init
	
	set w "$::KMProps::WinPath"

	if {$what == "OUTSIDE"} {
		
		# Open the window outside
		::KMProps::OpenWindowOutside $w
		
		if {![winfo exists $w]} return 
		
		#TAB PROPIEDADES y MATERIALES para separar los dos árboles
		set nb "$w.nb"
		pack [ttk::notebook $nb] -expand 1 -fill both
		
		#
		# PROPERTIES TREE 
		#
		set fProp ${nb}.fProp
		$nb add [ttk::frame $fProp -padding {1 1 1 1}] -text "[= Model]"
		
		set ::KMProps::TreePropsPath [::KMProps::CreateTreeAndToolbar $fProp]
		
		#
		# MATERIALS TREE
		#
		set fMat ${nb}.fMat
		$nb add [ttk::frame $fMat -padding {1 1 1 1}] -text "[= Materials]"
		
		set ::KMat::TreeMatPath [::KMat::CreateTreeAndToolbar $fMat]
		
		::KMat::initVisibilityClass
		::KMat::FillTreeMat
		
		::KMProps::initVisibilityClass
		::KMProps::FillTreeProps
		
		# Select active tab
		switch -exact -- $whattab {
		    "Model" {
		        $nb select "$nb.fProp"
		    }
		    "Materials" {
		        $nb select "$nb.fMat"
		    }
		}
	
		# Binding
		bind $w <Alt-c> "destroy $w"
		#bind $w <Escape> "destroy $w"
		
	} else {
		
	}
}

proc ::KMProps::CreateTreeAndToolbar { w } {
	
	# Create the treectrl properties 
	set mdf [ttk::frame $w.middle]
	set T [::KMProps::CreateTreeProperties $w]

	grid $w.middle -sticky wens

	#grid rowconfigure $w 2 -weight 1
	#grid columnconfigure $w 0 -weight 1
	
	## For lower buttons
	#set tf [ttk::frame $w.close] 
	#grid $tf -sticky ews
	#grid anchor $tf center
	#grid [ttk::button $tf.bClose -text [= "Close"] -command [list destroy $w]]  -sticky ew -padx 5 -pady 3

	focus $T

	return $T
}

#
# Crea el frame inferior 
#
proc ::KMProps::iniFrameBottom { } {
	
	#Destruye el frame inferior si existía
	set f [::KMProps::cancelBottom]
	
	# Create the frame where set the properties
	ttk::frame $f -borderwidth 0
	
	# Grid for toolbar
	grid $f -row 2 -column 0 -sticky wes
	
	return $f

}

#
# Destruye el frame inferior 
#
proc ::KMProps::cancelBottom { } {
	
	set f ${::KMProps::WinPath}.nb.fProp.fBottom
	if { [winfo exists $f]} {
		foreach w [winfo children $f] {
		        destroy $w
		}
	}
	destroy $f
	
	return $f
}

proc ::KMProps::initVisibilityClass { } {
	
	variable visibilityVars
	set visibilityVars {}
	global KPriv

	set classes [$KPriv(xml) set "/Kratos_Data/ClassConfiguration/Class" ]
	
	set classes [split $classes ","]
	
	foreach class $classes {
		lappend visibilityVars $class
		set ::KMProps::$class ""
	}
}

#---------------------------------------------------------------------------------------------- 
# Lee el xml y carga el árbol de propiedades de forma iterativa como máximo hasta 7 niveles
#----------------------------------------------------------------------------------------------
proc ::KMProps::FillTreeProps { } {
	
	variable dimension

	global KPriv
	
	#Obtenemos los grupos si aun no han sido cargados (aun no han cargado su ventana)
	if { [llength $KPriv(groupsId)] == 0 } {
		::KEGroups::getXmlGroupsId
	}
	
	set T $::KMProps::TreePropsPath
	
	#Seleccionamos todos los nodos del primer nivel
	set nodes [$KPriv(xml) selectNodes "/Kratos_Data/RootData\[@id\]"]
	
	foreach node $nodes {
		# Insertamos cada RootData de 1er nivel en el árbol
		set item [::KMProps::InsertNewProp $node [$node getAttribute id ""] $T "" "root" [$node hasChildNodes] [::KMProps::stateNode $node] [$node getAttribute open "0"]]
		if {$item != -1} {
		        
		        set nodes2 [$node childNodes]
		        foreach node2 $nodes2 {
		        set item2 [::KMProps::InsertNewProp $node2 [::KMProps::splitNode $node2] $T "[$node getAttribute id 0]//" "$item" [$node2 hasChildNodes] [::KMProps::stateNode $node2] [$node2 getAttribute open "0"]]
		        if {$item2 != -1} {
		                
		                #Seleccionamos los hijos (3º nivel)
		                set nodes3 [$node2 childNodes]
		                foreach node3 $nodes3 {
		                
		                set item3 [::KMProps::InsertNewProp $node3 [::KMProps::splitNode $node3] $T "[$node getAttribute id 0]//[::KMProps::splitNode $node2]//" "$item2" [$node3 hasChildNodes] [::KMProps::stateNode $node3] [$node3 getAttribute open "0"]]
		                if {$item3 != -1} {
		                        #Seleccionamos los hijos (4º nivel)
		                        set nodes4 [$node3 childNodes]                 
		                        foreach node4 $nodes4 {
		                        
		                        set item4 [::KMProps::InsertNewProp $node4 [::KMProps::splitNode $node4] $T "[$node getAttribute id 0]//[::KMProps::splitNode $node2]//[::KMProps::splitNode $node3]//" "$item3" [$node4 hasChildNodes] [::KMProps::stateNode $node4] [$node4 getAttribute open "0"]]
		                        if {$item4 != -1} {
		                                #Seleccionamos los hijos (5º nivel)
		                                set nodes5 [$node4 childNodes]
		                                foreach node5 $nodes5 {
		                                
		                                set item5 [::KMProps::InsertNewProp $node5 [::KMProps::splitNode $node5] $T "[$node getAttribute id 0]//[::KMProps::splitNode $node2]//[::KMProps::splitNode $node3]//[::KMProps::splitNode $node4]//" "$item4" [$node5 hasChildNodes] [::KMProps::stateNode $node5] [$node5 getAttribute open "0"]]
		                                if {$item5 != -1} {
		                                        #Seleccionamos los hijos (6º nivel)
		                                        set nodes6 [$node5 childNodes]   
		                                        foreach node6 $nodes6 {
		                                        
		                                        set item6 [::KMProps::InsertNewProp $node6 [::KMProps::splitNode $node6] $T "[$node getAttribute id 0]//[::KMProps::splitNode $node2]//[::KMProps::splitNode $node3]//[::KMProps::splitNode $node4]//[::KMProps::splitNode $node5]//" "$item5" [$node6 hasChildNodes] [::KMProps::stateNode $node6] [$node6 getAttribute open "0"]]
		                                        if {$item6 != -1} {
		                                                #Seleccionamos los hijos (6º nivel)
		                                                set nodes7 [$node6 childNodes]
		                                                
		                                                foreach node7 $nodes7 {
		                                                set item7 [::KMProps::InsertNewProp $node7 [::KMProps::splitNode $node7] $T "[$node getAttribute id 0]//[::KMProps::splitNode $node2]//[::KMProps::splitNode $node3]//[::KMProps::splitNode $node4]//[::KMProps::splitNode $node5]//[::KMProps::splitNode $node6]//" "$item6" [$node7 hasChildNodes] [::KMProps::stateNode $node7] [$node7 getAttribute open "0"]]
		                                                }
		                                        }}
		                                }}
		                        }}
		                }}
		        }}
		}}
	
	return ""
}

#
# Valida varias cosas para cada nodo
#  DIMENSION (2D / 3D)
#  STATE (normal, hidden, disabled)
#
proc ::KMProps::stateNode { node } {
	
	#Validamos para cada nodo si tiene que estar visible 
	#(en función de los valores elegidos en algunos combos)
	if { [$node nodeName] == "Item" } {
		
		#Salvedad para no mostrar la propiedad Thickness en algunos casos
		set id [$node getAttribute id ""]
		
		if { $id == "ElemType" } {
		        
		        set ::KMProps::ElemTypeThickness [$node getAttribute dv ""]
		        
		} elseif { $id == "Thickness" } {
		        
		        if { ![::KMProps::showThickness]} {
		                
		                return "hidden"
		        }
		}
	}
	
	set state [$node getAttribute state "normal"]
	#msg "\n[$node getAttribute id ""]        state:$state"
	
	#Si el estado es hiddenAll se oculta el nodo y toda su descendencia
	if {$state == "hiddenAll"} {
		return "-1"
	}
	if { [$node nodeName] == "Item" } {
		
		#Salvedad para no mostrar la propiedad Thickness en algunos casos
		set id [$node getAttribute id ""]
		if { $id == "PressureValue" } {
		        set id $id
		}
	}
		
	foreach var $::KMProps::visibilityVars {
		
		set globalVar [set ::KMProps::$var]
		
		set nodeValuesVar [split [$node getAttribute $var ""] ","]

		
		#Si el nodo tiene alguna restriccion de clase (p.ej. del tipo strucType=Shell)
		# y no coincide con el valor seleccionado, ocultamos el nodo        
		if { $nodeValuesVar != "" } {
		        
		        if { !($globalVar in $nodeValuesVar) } {
		        
		        #msg "$state:  --------- > g: $globalVar in nodeVals: $nodeValuesVar"
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
	}
	
	#Leemos la class del nodo para ver si requiere de acciones especiales (ocultar nodos)
	#Si no se llega aquí es porque el nodo no era visible, y en ese caso no se tiene que hacer
	set value [$node getAttribute dv ""]
	
	set class [$node getAttribute class ""]
	#Equivalente a Switch $class
	foreach var $::KMProps::visibilityVars {
		
		if {$var == $class} {
		
		#Caso especial para el solver de fluidos
		if { $var == "fluidSolvTyp" } {
		        
		        #El solver de fluidos está duplicado dependiendo de una variable prebia
		        set freeYesOrNo [$node getAttribute freeSurf ""]
		        
		        if { $freeYesOrNo == "" || $freeYesOrNo == $::KMProps::freeSurf } {
		        set ::KMProps::fluidSolvTyp "$value"
		        }
		} else {
		        #Caso general
		        set ::KMProps::$var $value
		}
		}
	}
		
	
	#Devolvemos el estado del nodo ("normal" por defecto)
	return $state
}

#
# Valida varias cosas para cada nodo
#  DIMENSION (2D / 3D)
#  STATE (normal, hidden, disabled)
#
proc ::KMProps::stateNode2 { node } {
	
	#Validamos para cada nodo si tiene que estar visible 
	#(en función de los valores elegidos en algunos combos)
	if { [$node nodeName] == "Item" } {
		
		#Salvedad para no mostrar la propiedad Thickness en algunos casos
		set id [$node getAttribute id ""]
		
		if { $id == "ElemType" } {
		        
		        set ::KMProps::ElemTypeThickness [$node getAttribute dv ""]
		        
		} elseif { $id == "Thickness" } {
		        
		        if { ![::KMProps::showThickness]} {
		                
		                return "hidden"
		        }
		}
		
		#Leemos la class del nodo para ver si requiere de acciones especiales (ocultar nodos)
		
		set value [$node getAttribute dv ""]
		
		set class [$node getAttribute class ""]
		#Equivalente a Switch $class
		foreach var $::KMProps::visibilityVars {
		        
		        if {$var == $class} {
		        
		        #Caso especial para el solver de fluidos
		        if { $var == "fluidSolvTyp" } {
		                
		                #El solver de fluidos está duplicado dependiendo de una variable prebia
		                set freeYesOrNo [$node getAttribute freeSurf ""]
		                
		                if { $freeYesOrNo == "" || $freeYesOrNo == $::KMProps::freeSurf } {
		                set ::KMProps::fluidSolvTyp "$value"
		                }
		        } else {
		                #Caso general
		                set ::KMProps::$var $value
		        }
		        }
		}
	} else {
		
		#Caso especial para application
		#set class [$node getAttribute class ""]
		#if { $class == "application" } {
		
		#set apliState [$node getAttribute state ""]
		
		#if {$apliState != "hiddenAll" } {
		#set ::KMProps::application [$node getAttribute id ""]
		#}
		#}
	}
	
	set state [$node getAttribute state "normal"]
	
	#Si el estado es hiddenAll se oculta el nodo y toda su descendencia
	if {$state == "hiddenAll"} {
		return "-1"
	}
	if { [$node nodeName] == "Item" } {
		
		#Salvedad para no mostrar la propiedad Thickness en algunos casos
		set id [$node getAttribute id ""]
		if { $id == "PressureValue" } {
		        set id $id
		}
	}
		
	foreach var $::KMProps::visibilityVars {
		
		set globalVar [set ::KMProps::$var]
		
		set nodeValuesVar [split [$node getAttribute $var ""] ","]

		
		#Si el nodo tiene alguna restriccion de clase (p.ej. del tipo strucType=Shell)
		# y no coincide con el valor seleccionado, ocultamos el nodo        
		if { $nodeValuesVar != "" } {
		        
		        if { !($globalVar in $nodeValuesVar) } {
		        
		        if {$var == "strucType" } {
		                if { $globalVar != "Generic"} { 
		                
		                return -1
		                }
		        } else {
		                
		                return -1
		        }
		        }
		}
	}
	
	#Devolvemos el estado del nodo ("normal" por defecto)
	return $state
}

#
# Separa cada node en "inicialNombreTag.idNodo"
#
proc ::KMProps::splitNode { node } {
	
	set id [$node getAttribute id ""]
	
	if { [$node tagName] == "Container"} {
		return "c.$id"
	} elseif { [$node tagName] == "Item"} {
		return "i.$id"
	} elseif { [$node tagName] == "Property"} {
		return "p.[$node getAttribute id ""]"
	} elseif { [$node tagName] == "Material"} {
		return "m.[$node getAttribute id ""]"
	} else {
		return "NoTree"
	}
}

#---------------------------------------------------------------------------------------------- 
# Lee RECURSIVAMENTE el xml y carga el árbol de propiedades
#----------------------------------------------------------------------------------------------
#proc ::KMProps::FillTree { T node path item} {
#                                
#                                set path "$path//[::KMProps::splitNode $node]"
#                                
#                                
#                                if { [$node hasChildNodes] } {
#                                                
#                                        set item [::KMProps::FillTree $T [lindex [$node childNodes] 0] $path $item]
#                                                
#                                } else {
#                                        set item [::KMProps::InsertProp [$node getAttribute pid ""] [::KMProps::splitNode $node] $T "" "$path" "$item" [$node hasChildNodes] ]
#                                        
#                                        set parentNode [$node parentNode]
#                                        $node delete
#                                        
#                                        set item [::KMProps::FillTree $T $parentNode $path $item]
#                                        return $item
#                                }
#                
#                return ""
#}

proc ::KMProps::refreshTree { {T ""} {onlySave 0} } {
	
	if {$T == ""} {
		set T $::KMProps::TreePropsPath
	}
	#Primero hay que asegurarse de que exista el árbol
	if { [winfo exists $::KMProps::WinPath] } {
		
		foreach item [$T item range 0 end] {
		        
		        set fullname [DecodeName [$T item tag names $item]]
		        
		        if {$fullname != "" && [::xmlutils::getXmlNodeName $fullname] != "Item"} {
		        
		        catch {
		                ::xmlutils::setXml $fullname open "write" [$T item isopen $item]
		        }
		       
		        }
		}
	}
	
	if { $onlySave == 0 } {
		
		set selectedItem [$T selection get]
		
		#Vuelve a crear la ventana y el tree  
		#::KMProps::InitBaseWindow
		
		#Vacía el árbol y lo vuelve a llenar
		::KMat::DeleteTree $T
		::KMProps::FillTreeProps
		
		#Si está disponible añadimos la selección anterior
		catch {
		        $T selection add $selectedItem
		        $T see [$T selection get]
		}
	}
	set ::KMProps::lastSelected {}
}

proc ::KMProps::refreshTreeOld { T } {
	
	array unset stateitems
	array unset selectionitems
	#try to mantain the previous state: open, selected
	foreach item [$T item range 0 end] {
		set fullnameencoded [$T item tag names $item]
		
		set stateitems($fullnameencoded) [$T item isopen $item]
		set selectionitems($fullnameencoded) [$T selection includes $item]
	}
	$T item delete all ;#coult try not to delete and rebuild all...
	
	::KMProps::InitBaseWindow
	
	foreach item [$T item range 0 end] {
		set fullnameencoded [$T item tag names $item]
		
		if { [info exists stateitems($fullnameencoded)] &&
		         !$stateitems($fullnameencoded) } {
		        $T item collapse $item
		}
		if { [info exists selectionitems($fullnameencoded)] &&
		         $selectionitems($fullnameencoded) } {
		        $T selection add $item
		}
	}
	set id [$T item id "tag [EncodeName [GiD_Info Project LayerToUse]]"]
	if { $id != "" } {
		$T activate $id
	}
	if { [$T selection count] == 1 } {
		$T see [$T selection get]
	}
}

proc ::KMProps::CreateTreeProperties {w} {
	
	variable SystemHighlight
	variable SystemHighlightText

	# Scrollbars
	set vsb $w.vsb1
	set hsb $w.hsb1

	# Create the treectrl and set the scrollbars
	set T [treectrl $w.t -xscrollcommand [list $hsb set] -yscrollcommand [list $vsb set]]
	ttk::scrollbar $vsb -orient vertical   -command [list $T yview]
	ttk::scrollbar $hsb -orient horizontal -command [list $T xview]
	
	# Set the height
	set height [font metrics [$T cget -font] -linespace]
	if {$height < 18} {
		set height 18
	}
	
	# Configure the treectrl
	$T configure -indent 15 -itemheight $height -selectmode browse \
		-showroot 0 -showrootbutton 0 -showbuttons 1 -showlines 1 \
		-highlightthickness 0 -borderwidth 0 -height 400 \
		-xscrollincrement 20 -yscrollincrement 20
	
	# Create the column identifier list
	#set collistid [list [= "Properties"] [= "C1"]]
	#set i 0
	#foreach cid $collistid {
	# $T column create -text $cid -tags C$i -weight 0
	# incr i
	#}
	$T column create -text [= " "] -tags C0 -weight 0
	
	# Configure the column weight and arrow
	$T column configure C0 -weight 1 -arrow up

	# Configure the column that have the tree
	$T configure -treecolumn C0
	
	# Intento de desactivar el drag
	$T column dragconfigure -enable no

	# Create elements
	$T element create elemImgAny image
	$T element create elemTxtRead text -fill [list $SystemHighlightText {selected focus}] -lines 1
	$T element create elemRectSel rect -fill [list $SystemHighlight {selected focus} gray {selected !focus}] -showfocus yes
	$T element create eWindow window
	
	
	# Create styles using the elements
	set S [$T style create styAnyRead]
	$T style elements $S {elemImgAny elemRectSel elemTxtRead }
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

	# List of lists: {column style element ...} specifying text elements
	# the user can edit
	TreeCtrl::SetEditable $T {

	}
	
	# List of lists: {column style element ...} specifying elements
	# the user can click on or select with the selection rectangle
	TreeCtrl::SetSensitive $T {
		{C0 styAnyRead elemTxtRead}
	}

	# List of lists: {column style element ...} specifying elements
	# added to the drag image when dragging selected items
	#TreeCtrl::SetDragImage $T {
	
	#}
	#{C0 styAnyRead elemTxtRead}
	
	# Some notify install
	$T notify install <Drag-receive>
	#$T notify install <Edit-accept>

	# Notify bind
	# TODO

	$T notify bind DragTag <Drag-receive> { ::KMProps::ReceiveDragGroups  %T %l %I }
	#$T notify bind EditTag <Edit-accept> { ::KMProps::SetPropsToRename %T %I %t }

	bind $T <Button-1> [list ::KMProps::ClickTree %x %y $T]
	bind $T <Double-Button-1> [list ::KMProps::DoubleClickTree %x %y $T]
	bind $T <Return> [list ::KMProps::IntroEvent $T]
	#   bind $T <Key-Delete> [list SetLayersToDelete $T]
	#   bind $T <Alt_L> [list InvertSelectionTableList $T]
	#   bind $T <Alt_R> [list InvertSelectionTableList $T]
	#   bind $T <Meta_L> [list InvertSelectionTableList $T]
	#   bind $T <Meta_R> [list InvertSelectionTableList $T]
	#                bind $T <F2> [list ::KMProps::BeginEditGroups $T]

	bind $T <Button-3> "[list ::KMProps::MenuContextual %W %x %y] ; break"

	# Grid the tree
	grid $T $vsb -sticky nsew
	grid configure $vsb -sticky ns
	grid $hsb -sticky ew
	grid remove $vsb $hsb
	bind $T <Configure> [list ConfigureListScrollbars $T $hsb $vsb]
	
	grid columnconfigure $w 0 -weight 1
	grid rowconfigure $w 0 -weight 1

	bindtags $T [list $T TreeCtrlFileList TreeCtrl [winfo toplevel $T] all]

	return $T

}

proc ::TreeCtrl::FileListMotion1 { a b c } {
	
}

proc ::KMProps::CreateTreeMaterials { w } {
	
	global KPriv
	variable TreeMatPath
	
	set ::KMProps::TreeMatPath "$w.treeMat"
}

#
###################################################################################################
#--------------------------------------------------------------------------------------------------
# Funciones para manejar el XML de propiedades (las genéricas están en /lib/xml/xmlutils.tcl
#--------------------------------------------------------------------------------------------------
###################################################################################################
#
# Prepara la query para utilizar las funciones de domNOde
proc ::KMProps::setXPath { path  } {
	
	set splitted [::KMProps::split2 $path //]
	
	set i 0
	
	foreach itemId $splitted {
		
		if { $i == 0 } {                
		        #El primer elemento será siempre del nivel 'RootData'
		        set xpath "/Kratos_Data/RootData\[@id='$itemId'\]"
		} else {
		        
		        if { [string index $itemId 0] == "c" } {
		        set xpath "$xpath/Container\[@id='[string range $itemId 2 end]'\]"
		        } else {
		        set xpath "$xpath/Item\[@id='[string range $itemId 2 end]'\]"
		        }
		}
		incr i
	}
	
	return $xpath
}

proc ::KMProps::getTemplateStructure { id {type "containers"}} {
	
	global KPriv
	
	set T $::KMProps::TreePropsPath
	
	set listTemplate {}
	
	if { $type == "containers" } {
		
		#Seleccionamos todos los containers del template y sus respectivos items
		set nodes [$KPriv(xml) selectNodes "/Kratos_Data/Templates/Template\[@id='$id'\]/Container"]
		
		foreach node $nodes {
		        set listChilds [$node getAttribute id ""]
		        foreach item [$node childNodes] {
		        
		        set listChilds [lappend listChilds [$item getAttribute id ""]]
		        }
		        set listTemplate [lappend listTemplate $listChilds]
		}
	} else { # "Items"
		set listTemplate [$KPriv(xml) selectNodes "/Kratos_Data/Templates/Template\[@id='$id'\]/Items"]
	}
	
	return $listTemplate
}

proc ::KMProps::copyTemplate { idTemplate fullname groupId clase } {
	
	global KPriv
	
	#set idTemplate "Forces"
	set id [::KMProps::getPropTemplate "$idTemplate" id]
	set pid [::KMProps::getPropTemplate "$idTemplate" pid]
	set help [::KMProps::getPropTemplate "$idTemplate" help]
	set icon [::KMProps::getPropTemplate "$idTemplate" icon]
	
	set xpath "[::xmlutils::setXPath $fullname]"
	
	set template [$KPriv(xml) set "/Kratos_Data/Templates/Template\[@id='$idTemplate'\]"]
	
	#No se puede insertar en el xml un fragmento con mas de un nodo, por eso utilizamos 
	#las etiquetas auxiliares "<gouptemplate>$template</gouptemplate>"
	set template "<groupTemplate>$template</groupTemplate>"
	
	if { $clase == "OnlyGetText" } {
		return $template
	}
	
	$KPriv(xml) lappend "$xpath/Container id=\"$groupId\" pid=\"$groupId\" class=\"$clase\" icon=\"$icon\" help=\"$help\" open=\"1\"" $template
	
	#Eliminamos ahora las etiquetas auxiliares "<gouptemplate></gouptemplate>"
	::KMProps::replaceTemplate
	
	return $template
}

proc ::KMProps::replaceTemplate { } {
	
	global KPriv
	
	set xmlText [$KPriv(xml) asXML]
	
	set xmlText [string map {"<groupTemplate>" "" "</groupTemplate>" ""} $xmlText]
	
	set KPriv(xmlDoc) [dom parse $xmlText]
	
	set KPriv(xml) [$KPriv(xmlDoc) documentElement]
}

proc ::KMProps::getPropTemplate {id property {templatePath ""}} {
	
	global KPriv
	
	set xpath "/Kratos_Data/Templates/Template\[@id='$id'\]"
	
	set splitted [::KMProps::split2 $templatePath //]
	if {[llength $splitted] >= 1} {

		set idContainer [lindex $splitted 0]
		if { $idContainer != "OnlyItems" } {
		        set xpath "$xpath/Container\[@id='$idContainer'\]"
		}
	}
	if { [llength $splitted] >= 2 } {
		set idItem [lindex $splitted 1]
		set xpath "$xpath/Item\[@id='$idItem'\]"
	}

	set value [$KPriv(xml) set "$xpath/@$property" ]
	
	return [lindex $value 0]
	
}

proc ::KMProps::itemType {fullname} {
	
	set splitted [::KMProps::split2 $fullname //]
	
	if { [llength $splitted] < 2 } {
		return "r"
	} else {
		set letra [string index [lindex $splitted end] 0]
		#Devuelve identificador de item (en esta caso "c" o "i" de container o item)
		return $letra
	}
}


# - END XML FUNCTIONS

###################################################################################################
#
# Dado el path de un item (fullname) te devuelve la letra inicial del último elemento 
#

proc ::KMProps::DoubleClickTree { x y T {item ""}} {
	
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
	#set col [lindex $info 3]
	
	set fullname [DecodeName [$T item tag names $item]]
	set idFull [string map { "." "" "//" ""} $fullname]
	
	if { [::xmlutils::setXml $fullname state] == "disabled" } {
		return ""
	}
	
	set id [::xmlutils::setXml $fullname id]
	set dv [::xmlutils::setXml $fullname dv]
	
	#Eliminamos el anterior combo, si aun está visible
	if {[llength $lastSelected] > 0 && [lindex $lastSelected 0] != $item} { 
		
		#set state [::xmlutils::getComboBoxState $fullname]
		::KMProps::cmbSelectChange [lindex $lastSelected 0] $T 1 "anterior"
	}
	
	#Destruimos el frame inferior
	::KMProps::cancelBottom
	
	if { [::xmlutils::getXmlNodeName $fullname]  == "Item" } {
		
		set f [::KMProps::buildFrame $T $item]  
		$T item style set $item C0 styFrame
		$T item element configure $item C0 eWindow -window $f
		
		#Nos guardamos este item como y el valor seleccionado (es el último visible)
		set f "$T.f$idFull.cmb"
		set selCombo [::xmlutils::getComboDv $f $fullname]
		set selComboText [::xmlutils::getComboDv $f $fullname "text"]
		set lastSelected [list $item $selCombo $selComboText]
		
	} else {
		
		set clase [::xmlutils::setXml $fullname class]
		set idTemplate [::xmlutils::setXml $fullname idTemplate]
		
		
		if {$clase == "Groups"} {
		        #Comprobamos si en esta asignación de Grupos necesitarán propiedades
		        set template [::KMProps::copyTemplate ${idTemplate} $fullname "Nothing" "OnlyGetText"]
		        
		        #Validamos que exista algún combo de "properties" en el template
		        if {[string match "*GCV=\"Properties*" $template]} {
		        
		        #Miramos si hay alguna propiedad dada de alta (si la hay,el combo no puede estar vacío)
		        set props [::KMProps::getProps $fullname]
		        if { [llength $props] == 0 } {
		                WarnWin [= "You must define a valid Property for this element type."]
		                
		                #Abrir la edición de propiedades
		                $T selection clear
		                set parentItem [$T item parent [$T item parent $item]]
		                foreach i [$T item children $parentItem] {
		                if { [$T item text $i 0] == "Properties" } {
		                        #$T selection add $i
		                        ::KMProps::DoubleClickTree 0 0 $T $i
		                }
		                }
		                
		                return  ""
		        }
		        }
		        ::KMProps::buildGroupsFrame $T $idTemplate $item $fullname
		        
		} elseif {$clase == "Properties" } { ;# || $clase == "Property"
		        
		        ::KMProps::buildPropertyFrame $T $idTemplate $item $fullname
		        
		} elseif {$clase == "Tab" || $clase == "Property" || $clase == "Group"} {
		        
		        ::KMProps::buildTabFrame $T $item $clase
		}
	}
}

proc ::KMProps::ClickTree { x y T } {
	
	variable lastSelected
	
	set info [$T identify $x $y]
	
	if { [lindex $info 0] == "item" && [llength $info] >= 4 } {
		
		set item [lindex $info 1]
		set col [lindex $info 3]
		
		set fullname [DecodeName [$T item tag names $item]]
		set id [::xmlutils::setXml $fullname id]
		
		#Eliminamos el anterior combo, si aun está visible
		if {[llength $lastSelected] > 0 && [lindex $lastSelected 0] != $item} { 
		        
		        ::KMProps::cmbSelectChange [lindex $lastSelected 0] $T 1 "anterior"
		}
		
		#Obtenemos la clase de nodo para actuar en consecuencia
		set clase [::xmlutils::setXml $fullname class]
		
		#Miramos si se ha pulsado container o item
		set nodeName [::xmlutils::getXmlNodeName $fullname]
		if { $nodeName == "Container" } {
		        #Abrimos el item si tiene hijos y no es de alta de propiedades o grupos
		        if { $clase == "Groups" || $clase == "Properties" } {
		        
		        #En este caso desplegamos el frame de alta de propiedades o grupos
		        ::KMProps::DoubleClickTree $x $y $T
		        } else {
		        #$T item toggle $item
		        }
		} elseif {$nodeName == "Item"} {
		        
		}
		
		
	} elseif { [lindex $info 0] == "header" && [lindex $info 1] == "0" } {
		
		if { [$T column cget C0 -arrow] == "up" } {
		        $T column configure C0 -arrow down
		        #$T item sort 0 -dictionary -increasing
		} else {
		        $T column configure C0 -arrow up
		        #$T item sort 0 -dictionary -decreasing
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
		WarnWin "column 1 ???"
	}
	
	if { $col != 0 } {
		return -code break
	}
	return ""
	
}

proc ::KMProps::IntroEvent { T } {
	
	::KMProps::DoubleClickTree 0 0 $T
}

#
# Construye un frame con un tab por cada container que cuelgue de "item"
# y dentro de cada tab, etiquetas y combos para cada item (si los hay)
#
proc ::KMProps::buildGroupsFrame { T idTemplate item fullname} {
	
	global KPriv
	
	set f [::KMProps::iniFrameBottom]
	
	#El combo de propiedades es importante resetearlo
	# (posteriormente se utiliza para validar si hay alguna seleccionada)
	if {[info exists ::KMProps::cmbProperty]} {
		unset ::KMProps::cmbProperty
	}
	#Parametro para utilizar o no los combos "Activation"
	set activation 0
	set listT [::KMProps::getTemplateStructure $idTemplate]
	
	if {[llength $listT] >= 1 } {
		set nb ${f}.nb
		grid [ttk::notebook $nb ] -row 0 -column 0 -columnspan 2 -padx 0 -sticky nw -in $f
		
		
		#Lista de listas con formato {idContainer idItem1 idItem2...}
		foreach listContainer $listT {
		        
		        #Si tiene como mínimo el container y un item ponemos tab y dentro label-combos
		        if {[llength $listContainer] >= 2} {
		        
		        set idContainer [lindex $listContainer 0]
		        
		        #Si solo hay un container le damos el nombre del item pulsado, no del template
		        if {[llength $listT] == 1} {
		                set pid "[$T item text $item 0]"
		        } else {
		                set pid [= [::KMProps::getPropTemplate $idTemplate pid $idContainer]]
		                #Lo limitamos por si es demasiado largo
		                set pid "[string range $pid 0 20]"
		        }
		        
		        # Para cada container declaramos un tab
		        set fTab ${nb}.f$idContainer
		        $nb add [ttk::labelframe $fTab -text "[= Properties]" -padding {10 0 10 10}] \
		                -text "$pid"
		        
		        for {set i 1} { $i < [llength $listContainer] } {incr i} {
		                
		                set id [lindex $listContainer $i]
		                
		                #Si no coincide la dimensión 2D/3D no lo ponemos
		                set nDim [::KMProps::getPropTemplate $idTemplate nDim "$idContainer//$id"]
		                if { $nDim == "" || $nDim == $::KMProps::nDim } {
		                
		                #Comprobamos el estado para ver si está activo o no
		                set state [::KMProps::getPropTemplate $idTemplate state "$idContainer//$id"]
		                
		                if {$state != "hidden"} {
		                        
		                        set pid [::KMProps::getPropTemplate $idTemplate pid "$idContainer//$id"]
		                        set tooltip [::KMProps::getPropTemplate $idTemplate tooltip "$idContainer//$id"]
		                        set dv [::KMProps::getPropTemplate $idTemplate dv "$idContainer//$id"]
		                        set function [::KMProps::getPropTemplate $idTemplate function "$idContainer//$id"]
		                        
		                        #Obtenemos la lista de valores para el combo si existe
		                        set comboList [::xmlutils::getXMLValues "$idContainer//$id" $idTemplate "" $fullname]
		                        set icomboList [::xmlutils::getXMLValues "$idContainer//$id" $idTemplate "ivalues" $fullname]
		                        
		                        set CBState [::KMProps::getPropTemplate $idTemplate CBState "$idContainer//$id"]
								if { $CBState == "normal" } {
									set values $comboList
									set comboList {}
								} else {
									set values {}
								}
		                        
		                        #Para cada item añadimos label y combo
		                        grid [ttk::label $fTab.lbl$id -text "$pid:" ] \
		                        -row $i -column 0 -pady 2 -sticky nw -in $fTab
		                        
		                        if { [llength $comboList] > 0 } {
		                        
		                                if {$state != "disable"} {
		                                        set state "readonly"
		                                }
		                                if {[string length $id] == 2 && ([string index $id end] == "x" || [string index $id end] == "y" || [string index $id end] == "z") } { 
		                                        set width 15
		                                } else {
		                                        set width 20
		                                }
		                                grid [ttk::combobox $fTab.cmb$id -values $comboList -state $state -width $width -textvariable "::KMProps::cmb$id"] \
		                                        -row $i -column 1 -padx 3 -pady 2 -sticky ne -in $fTab
		                                tooltip::tooltip $fTab.cmb$id [= "%s" $tooltip]
		                                
		                                if {$id == "Ax" || $id == "Ay" || $id == "Az"} {
		                                        bind $fTab.cmb$id <<ComboboxSelected>> "::KMProps::cmbDisable $fullname $f.nb $id"
		                                        set activation 1
		                                }
		                                ::xmlutils::setComboDv $fTab.cmb$id $fullname $dv $idTemplate 
		                                
		                        } else {
		                        
		                        if {$id == "Vx" || $id == "Vy" || $id == "Vz"} {
		                                if { $activation } {
		                                        set activeId "A[string range $id 1 1]"
		                                        set value [set "::KMProps::cmb$activeId"]
		                                        if { $value == 0 } {
		                                                set state "disabled"
		                                        }
		                                }
		                        }
		                        
		                        grid [ttk::combobox $fTab.cmb$id -width 20 -state $state -values $values -textvariable "::KMProps::cmb$id"] \
		                                -row $i -column 1 -padx 3 -pady 2 -sticky nw -in $fTab
		                        tooltip::tooltip $fTab.cmb$id [= "%s" $tooltip]
		                        
		                        #set dv [::KMProps::getPropTemplate $idTemplate dvText "$idContainer//$id"]
		                        set ::KMProps::cmb$id $dv
		                        
		                        }
		                        if {$function != "" } {
		                                grid [ttk::button $fTab.funct$id -text "functions" -command "KFun::InitBaseWindow $fTab $id" -style TMenubutton.Toolbutton] \
		                                -row $i -column 2 -sticky nw  -pady 0 -padx 3 -in $fTab
		                                tooltip::tooltip $fTab.funct$id [= "Function manager"]
		                                set img [::WinUtils::GetImage "functions.gif"]
		                                if { $img != -1 } { $fTab.funct$id configure -image $img }
		                                
		                                grid [ttk::button $fTab.deleteFunct$id -text "delete" -command "::KMProps::unassignFunction $fTab $id $fullname" -style TMenubutton.Toolbutton] \
		                                -row $i -column 3 -sticky nw  -pady 0 -padx 2 -in $fTab
		                                tooltip::tooltip $fTab.deleteFunct$id [= "Unassign function"]
		                                set img [::WinUtils::GetImage "delete_icon.gif"]
		                                if { $img != -1 } { $fTab.deleteFunct$id configure -image $img }
		                                
		                        }
		                }
		                }                                                
		        }
		        }                                                
		}
		
		#
		# PUNTOS, LINEAS, SUPERFICIES Y VOLUMENES
		#
		set whatuse [GiD_Info Project ViewMode]
		set geomlist [GetGeometryImageFiles]
		set meshlist {node.gif element.gif}
		set col 0
		
		#Miramos si hay restricción de entidades
		set entityList [split [::xmlutils::setXml $fullname GiDEntity] "," ]
		
		#set entityList [split [::KMProps::getPropTemplate $idTemplate GiDEntity] "," ]
		
		switch $whatuse {
		        
		        GEOMETRYUSE {
		        foreach i $geomlist {
		                
		                set command [file rootname $i]
		                
		                if {$command in $entityList || $entityList == "" } {
		                if {$col == 0} {
		                        #Por defecto dejamos marcada la primera
		                        set i "[string range $i 0 [expr [string length $i] - 5]]_sel.gif"
		                        set ::KMProps::selectedEntity $command
		                }
		                
		                set fb "${f}.b$command"
		                grid [ttk::button $fb -text "$i" -command "::KMProps::changeImage $command $i $f" -style TMenubutton.Toolbutton] \
		                        -row 1 -column 0 -sticky nw  -pady 3 -padx [expr (50 * $col) + 15] -in $f
		                tooltip::tooltip $fb [= "Entity %s" $command]
		                set img [::WinUtils::GetImage $i]
		                if { $img != -1 } {
		                        $fb configure -image $img
		                }
		                incr col
		                }
		        }
		        }
		        MESHUSE {
		        foreach i $meshlist {
		                set command [file rootname $i]
		                if {$command in $entityList || $entityList == ""  } {
		                if {$col == 0} {
		                        #Por defecto dejamos marcada la primera
		                        set i "[string range $i 0 [expr [string length $i] - 5]]_sel.gif"
		                        set ::KMProps::selectedEntity $command
		                }
		                
		                set fb "${f}.b$command"
		                grid [ttk::button $fb -text "$i" -command "::KMProps::changeImage $command $i $f" -style TMenubutton.Toolbutton] \
		                        -row 1 -column 0 -sticky nw  -pady 3 -padx [expr (50 * $col) + 15] -in $f
		                tooltip::tooltip $fb [= "Entity $command"]
		                set img [::WinUtils::GetImage $i]
		                if { $img != -1 } {
		                        $fb configure -image $img
		                }
		                incr col
		                }
		        }
		        }
		}
		
		grid [ttk::label $f.lGroups -text "Group:" ] \
		        -row 2 -column 0 -pady 3 -padx 5 -sticky nw -in $f
		
		#COMBO DE GRUPOS
		set fGroups $f.cGroups
		set filterGroups [::KMProps::getGroups $entityList]
		grid [ttk::combobox $fGroups -state readonly -values "$filterGroups" -textvariable "::KMProps::selGroup"  -postcommand "::KMProps::changeGroups [list $entityList] $fGroups" -width 20] \
		        -row 2 -column 0 -pady 3 -padx 55 -sticky nw -in $f
		
		set ::KMProps::selGroup ""
		if { [llength $filterGroups] > 0 } {
		        $f.cGroups current 0
		}
		bind $f.cGroups <<ComboboxSelected>> "::KMProps::cmbChangeCheckGroups $f"
		
		# BOTON A LA DERECHA DE LOS GRUPOS (CREAR GRUPO AUTOMATICAMENTE)
		grid [ttk::button $f.iGroups -text [= "newGroup"] -command "::KMProps::autoNewGroup $id" ] \
		        -row 2 -column 0 -sticky nw  -pady 3 -padx 230 -in $f
		tooltip::tooltip $f.iGroups [= "Create automatic new group"]
		$f.iGroups configure -image [::WinUtils::GetImage "newAutoGroup.gif" ]
		
		grid [ttk::button $f.bPropOk -text [= "Ok"]  -command "::KMProps::acceptGroups $T $idTemplate $fullname $item {$listT} {$entityList} $fGroups" ] \
		        -row 3 -column 0 -sticky nw  -pady 3 -padx 20  -in $f
		tooltip::tooltip $f.bPropOk [= "Assign condition to the selected group"]
		
		grid [ttk::button $f.bPropCancel -text [= "Cancel"]  -command "::KMProps::cancelBottom" ] \
		        -row 3 -column 0 -sticky nw  -pady 3 -padx 100  -in $f
		tooltip::tooltip $f.bPropCancel [= "Cancel assignation"]

	}
	
	bind $T <KeyPress> "if { %k == 27   } { ::KMProps::cancelBottom }"
}

proc ::KMProps::changeGroups { entityList f {fullname ""} } {
	
	set valores [::KMProps::getGroups $entityList $fullname]
	$f configure -values $valores
	
	if { !($::KMProps::selGroup in $valores) } {
		
		if {[string range $::KMProps::selGroup 0 8] == "AutoGroup"} {
		        WarnWin [= "The new group '%s' has not any usefull entity assigned." $::KMProps::selGroup]
		}
		set ::KMProps::selGroup ""
	} 
}

proc ::KMProps::changeCmbValues { f path {idTemplate ""} {elemTypeDv ""} {noTemplateFullname ""} } { 
	
	set values [::xmlutils::getXMLValues $path $idTemplate "" "$noTemplateFullname" "$elemTypeDv"]
	
	$f configure -values $values
	
	if { [llength $values] > 0 && [$f current] == -1 } {
		$f current 0        
	}
}

proc ::KMProps::cmbChangeCheckGroups { f } {
	
	global KPriv
	
	if { [winfo exists $f.cGroups] } {
		
		if { [llength $KPriv(groupsId)] > 0 } {
		        
		        $f.cGroups configure -values $KPriv(groupsId)
		        if { $::KMProps::selGroup ni $KPriv(groupsId) } {
		        set ::KMProps::selGroup [lindex $KPriv(groupsId) 0]
		        }
		} else {
		        $f.cGroups configure -values {}
		        set ::KMProps::selGroup ""
		}
	}
}

proc ::KMProps::unassignFunction { f id fullname } {
	
	set fcmb "$f.cmb$id"
	
	if {[winfo exists $fcmb]} {
		$fcmb configure -state normal
	}
	
	set ::KMProps::cmb$id "0.0"
	
	#set fFun "${::KMProps::WinPath}.functions"
	#if {[winfo exists $fFun]} {        
	#        destroy $fFun
	#}
	
}

proc ::KMProps::assignFunction { f id idFunction } {
	
	set fcmb "$f.cmb$id"
	
	if {[winfo exists $fcmb]} {
		$fcmb configure -state disabled
	}
	set ::KMProps::cmb$id "$idFunction"
}

proc ::KMProps::getGroups { entityList {fullname ""}} {
	
	global KPriv
	
	# Switch state
	if { $entityList == "" } {
		set PState [GiD_Info Project ViewMode]
		if {($PState == "GEOMETRYUSE")} {
		        set entityList {point line surface volume}
		} else {
		        set entityList {nodes elements}
		}
	}
	
	set grupos {}
	
	foreach groupId $KPriv(groupsId) {
		foreach entity $entityList {
		        if { [::KEGroups::getGroupGiDEntities $groupId $entity "hasEntities"] } {
		        if { !( $groupId in $grupos) } {
		                lappend grupos $groupId
		        }
		        }
		}
	}
	
	#foreach groupId $KPriv(groupsId) {
	#set gEntities [::KEGroups::getAssignedGiDEntities $groupId]
	#set geomEntities [lindex $gEntities 0]
	#foreach "points lines surfaces volumes" $geomEntities {
	#        
	#        # Si el grupo tiene alguna de esas entidades asignadas lo añadirá al combo
	#        if { "point" in $entityList && [llength $points] } {
	#        lappend grupos $groupId
	#        } elseif {"line" in $entityList && [llength $lines] } {
	#        lappend grupos $groupId
	#        } elseif {"surface" in $entityList && [llength $surfaces] } {
	#        lappend grupos $groupId
	#        } elseif {"volume" in $entityList && [llength $volumes] } {
	#        lappend grupos $groupId
	#        }
	#}
	#}
	
	if {$fullname != ""} {
		
		#Eliminamos de la lista los grupos ya asignados a esta propiedad
		set assignedGroups [::xmlutils::setXmlContainerIds $fullname]
		foreach g $assignedGroups {
		        if { $g != $::KMProps::selGroup } {
		        set grupos [::KEGroups::listReplace $grupos $g]
		        }
		}
	}
	return $grupos
}


proc ::KMProps::cmbElemTypeChange { f fullname {idTemplate ""} } {
	
	global KPriv
	
	if { [info exists ::KMProps::cmbElemType] } {
		
		set dv [::xmlutils::getComboDv $f $fullname "id" $idTemplate]
		
		set xpath "Kratos_KWords/ElementCLaws/Item\[@id='ElementTypes'\]"
		::xmlutils::getAttribute $KPriv(xmlDocKKW) $xpath dv $dv
		
		#Actualizamos los valores de Material Mode en función del nuevo dv
		set fMat [string map {"ElemType" "MatModel"} $f]
		
		if { $idTemplate == "" } {
		        
		        set fullnameMat [string map {"ElemType" "MatModel"} $fullname]
		        ::KMProps::changeCmbValues $fMat "$fullnameMat" $idTemplate "NoSearchFullname"
		} else {
		        
		        ::KMProps::changeCmbValues $fMat "MainProperties//MatModel" $idTemplate
		}
		
		set fCmbThick [string map {"ElemType" "Thickness"} $f]
		set fLblThick [string map {"cmbThickness" "lblThickness"} $fCmbThick]
		if { [::KMProps::showThickness] } {

		        if { [winfo exists $fCmbThick] } { 
		                grid $fCmbThick 
		        }
		        if { [winfo exists $fLblThick] } { 
		                grid $fLblThick 
		        }
		} else {
		        if { [winfo exists $fCmbThick] } { 
		                grid remove $fCmbThick 
		        }
		        if { [winfo exists $fLblThick] } { 
		                grid remove $fLblThick 
		        }
		
		}
		
	}
}

#
#Consulta el xml KKWord para ver si se tiene que mostrar el elemento del Tickness
#
proc ::KMProps::showThickness { } {
	
	#return 1
	
	global KPriv
	
	#Leemos el xml si aun no se ha leído
	#::xmlutils::initKKWord
	
	if { [info exists ::KMProps::ElemTypeThickness] } {
		set dv $::KMProps::ElemTypeThickness
		unset ::KMProps::ElemTypeThickness
	} else {
		set xpath "Kratos_KWords/ElementCLaws/Item\[@id='ElementTypes'\]"
		set dv [::xmlutils::getAttribute $KPriv(xmlDocKKW) $xpath dv]
	}
	
	set xpath "Kratos_KWords/ElementCLaws/Item\[@id='Tickness$::KMProps::nDim'\]"
	set ListElementType [split [::xmlutils::getAttribute $KPriv(xmlDocKKW) $xpath elementType] ","]
	
	if { $dv in $ListElementType } {
		return 1
	} else {
		return 0
	}
}


#
# Construye un frame consultando el template de propiedades correspondiente,
#  con un tab por cada container y dentro de cada tab, 
#  etiquetas y combos para cada item (si los hay)
#
proc ::KMProps::buildPropertyFrame { T idTemplate item fullname } { 
	
	global KPriv
	
	set f [::KMProps::iniFrameBottom]
	
	set listT [::KMProps::getTemplateStructure $idTemplate]
	
	if {[llength $listT] >= 1 } {
		set nb ${f}.nb
		grid [ttk::notebook $nb ] -row 0 -column 0 -columnspan 2 -padx 0 -sticky nw -in $f
		
		#Lista de listas con formato {idContainer idItem1 idItem2...}
		foreach listContainer $listT {
		        
		        #Si tiene como mínimo el container(1er elemento) y un item ponemos tab y dentro label-combos
		        if {[llength $listContainer] >= 2} {
		        
		        set idContainer [lindex $listContainer 0]
		        
		        set pid [::KMProps::getPropTemplate $idTemplate pid $idContainer]
		        
		        # Para cada container declaramos un tab
		        set fTab ${nb}.f$idContainer
		        $nb add [ttk::labelframe $fTab -text "[= Properties]" -padding {10 0 10 10}] \
		                -text "[string range $pid 0 20]"
		        
		        #En el caso del primer tab forzamos el item de "Nombre de propiedad" y 2campos mas
		        if { $idContainer == [lindex $listT 0 0] } {
		                
		                #Nombre de propiedad
		                grid [ttk::label $fTab.lblName -text "[= Property Name:]" ] \
		                -row 0 -column 0 -pady 5 -sticky nw -in $fTab
		                grid [ttk::combobox $fTab.cmbPropertyName -state normal -textvariable "::KMProps::propertyName" ] \
		                -row 0 -column 1 -padx 3 -pady 5 -sticky nw -in $fTab
		                tooltip::tooltip $fTab.cmbPropertyName [= Choose a new property name. ]
		                
		                set ::KMProps::propertyName "[::KEGroups::GetAutomaticPropertyName $fullname]"
		                focus $fTab.cmbPropertyName
		        }
		        
		        for {set i 1} { $i < [llength $listContainer] } {incr i} {
		                
		                set id [lindex $listContainer $i]
		                
		                #Los nodos ocultos no se deben mostrar
		                set state [::KMProps::getPropTemplate $idTemplate state "$idContainer//$id"]
		                if {$state != "hidden" } {
		                        
		                #Si no coincide la dimensión 2D/3D no lo ponemos
		                set nDim [::KMProps::getPropTemplate $idTemplate nDim "$idContainer//$id"]
		                if { $nDim == "" || $nDim == $::KMProps::nDim } {
		                        
		                        set pid [::KMProps::getPropTemplate $idTemplate pid "$idContainer//$id"]
		                        set tooltip [::KMProps::getPropTemplate $idTemplate tooltip "$idContainer//$id"]
		                        set dv [::KMProps::getPropTemplate $idTemplate dv "$idContainer//$id"]
		                        
		                        set CBState [::KMProps::getPropTemplate $idTemplate CBState "$idContainer//$id"]
								if { $CBState == "normal" } {
									set values $comboList
									set comboList {}
								} else {
									set values {}
								}
		                        
		                        # Es importante que ElemType sea el primer combo para actualizar el "dv" en el xml
		                        if { $id == "ElemType" } {
		                                set comboList [::xmlutils::getXMLValues "$idContainer//$id" $idTemplate "" "" $dv]
		                        } else {
		                                set comboList [::xmlutils::getXMLValues "$idContainer//$id" $idTemplate "" "$fullname"]
		                        }
		                        
		                        #Para cada item añadimos label y combo
		                        grid [ttk::label $fTab.lbl$id -text "$pid:" ] \
		                                -row $i -column 0 -pady 5 -sticky nw -in $fTab
		                        
		                        if { [llength $comboList] > 0 } {
		                                grid [ttk::combobox $fTab.cmb$id -values $comboList -state readonly \
		                                          -textvariable "::KMProps::cmb$id" \
		                                          -postcommand [list ::KMProps::changeCmbValues "$fTab.cmb$id" "$idContainer//$id" "$idTemplate" "" "$fullname"] ] \
		                                -row $i -column 1 -padx 5 -pady 2 -sticky ne -in $fTab
		                                
		                                ::xmlutils::setComboDv $fTab.cmb$id $fullname $dv $idTemplate
		                                
		                                #En este caso se tendrá qué recargar si existe el combo de "Material Model"
		                                if { $id == "ElemType" } {
		                                        bind $fTab.cmb$id <<ComboboxSelected>> "::KMProps::cmbElemTypeChange $fTab.cmb$id $idContainer//$id $idTemplate"
		                                }
		                                
		                        } else {
		                                grid [ttk::combobox $fTab.cmb$id -state normal -values $values -textvariable "::KMProps::cmb$id"] \
		                                -row $i -column 1 -padx 5 -pady 2 -sticky nw -in $fTab
		                                set ::KMProps::cmb$id $dv
		                        }
		                        tooltip::tooltip $fTab.cmb$id $tooltip
		                        
		                        #En el caso de ElemType actualizamos la variable de filtrado del Thickness
		                        if { $id == "ElemType" } {
		                                set ::KMProps::ElemTypeThickness $dv
		                        
		                        #En el caso del Thickness, no siempre se mostrará
		                        } elseif { $id == "Thickness" && ![::KMProps::showThickness]} {
		                                
		                                grid remove $fTab.lbl$id
		                                grid remove $fTab.cmb$id
		                        }
		                }
		                }
		        }
		        }                                                
		}
	}

	grid [ttk::button $f.bPropOk -text [= "Ok"]  -command "::KMProps::acceptProperty $T $idTemplate $fullname $item {$listT}" ] \
		-row 3 -column 0 -sticky nw  -pady 3 -padx 20  -in $f
	tooltip::tooltip $f.bPropOk [= "Assign condition to the selected group"]
	
	grid [ttk::button $f.bPropCancel -text [= "Cancel"]  -command "::KMProps::cancelBottom" ] \
		-row 3 -column 0 -sticky nw  -pady 3 -padx 100  -in $f
	tooltip::tooltip $f.bPropCancel [= "Cancel assignation"]
	
	bind $T <KeyPress> "if { %k == 27   } { ::KMProps::cancelBottom }"

}

proc ::KEGroups::GetAutomaticPropertyName { fullname } {
	
	set name ""
	
	set propIds [::xmlutils::setXmlContainerIds "[::KMProps::getApplication $fullname]//c.Properties"]
	set i 0
	
	if { [llength propIds] > 0 } {
		
		for {set i 1} {$i<10000} {incr i} {
		        
		        set name "Property${i}"
		        if { [lsearch -exact $propIds $name] == -1 } { break }
		}
	} else {
		set name "Property1"
	}
	return $name
}

#
# Click en botón OK de Property
#
proc ::KMProps::acceptProperty { T idTemplate fullname item listT} {
	
	set property $::KMProps::propertyName
	
	#Validamos que la propiedad no tenga carácteres extraños
	set property [::KUtils::parseTreeStr $property]
	if { $property == -1 } {
		
		WarnWin [= "You can't use some reservate chars like:\n  :   /   $   .   \\  %  "]
		set ::KMProps::propertyName ""
		return ""
	}
	
	#Comprobamos q el nombre no sea vacío
	if { $property == "" } {
		WarnWin [= "The property name can not be empty"]
		return ""
	}
	
	#Comprobamos q la propiedad aun no exista
	set id [::xmlutils::setXml "$fullname//c.$property" id]

	if { $id != "" } {
		
		WarnWin [= "This name property it is already assigned."]
		return
	}

	::KMProps::copyTemplate ${idTemplate} $fullname "$property" "Property"
	
	#
	# Ahora debemos actualizar todos los valores en el xml
	#
	if {[llength $listT] >= 1 } {
		
		#Lista de listas con formato {idContainer idItem1 idItem2...}
		foreach listContainer $listT {
		        
		        #Si tiene como mínimo el container y un item entramos
		        if {[llength $listContainer] >= 2} {
		        
		        set idContainer [lindex $listContainer 0]
		        
		        #Recorremos los items
		        for {set i 1} { $i < [llength $listContainer] } {incr i} {
		                
		                set id [lindex $listContainer $i]
		                
		                #if { $id != "Thickness" || ($id == "Thickness" && [::KMProps::showThickness])} {
		                
		                        set fullNombre "$fullname//c.$property//c.$idContainer//i.$id"
		                        
		                        #Los nodos ocultos no existían en el formulario
		                                set state [::xmlutils::setXml $fullNombre state]
		                                if {$state != "hidden" } {
		                                
		                                #Validamos la dimensión de cada elemento
		                                set nDim [::xmlutils::setXml $fullNombre nDim]
		                                if { $nDim == "" || $nDim == $::KMProps::nDim } {
		                                
		                                        set value [set ::KMProps::cmb$id]
		                                        
		                                        ::xmlutils::setXml $fullNombre dv "write" $value
		                                }
		                        }
		                #}
		        }
		        }
		}
	}

	#Destruimos el frame inferior
	::KMProps::cancelBottom
	
	#Recargar tree
	::KMProps::refreshTree $T
}

proc ::KMProps::changeImage {entity img path} {
	
	variable selectedEntity
	
	if {$selectedEntity != "" } {
		
		set f "$path.b$selectedEntity"
		
		$f configure -image [::WinUtils::GetImage "${selectedEntity}.gif"]
	}
	
	"$path.b$entity" configure -image [::WinUtils::GetImage "${entity}_sel.gif"]

	set selectedEntity $entity
}

proc ::KMProps::autoNewGroup { id } {
	
	global KPriv
	
	set GroupId [::KEGroups::GetAutomaticGroupName "Auto"]
	
	#lappend KPriv(groupsId) $GroupId
	
	set color [::KEGroups::randomColor]
	::KEGroups::insertXml "root" $GroupId $color 1 Generic
	
	set f ${::KMProps::WinPath}.nb.fProp.fBottom
	
	#$f.cGroups configure -values $KPriv(groupsId)
	#msg "selectedEntity:$::KMProps::selectedEntity id:$id groupsId:$KPriv(groupsId)"
	set ::KMProps::selGroup $GroupId 
	
	#Cerramos la ventana de grupos
	destroy $::KEGroups::WinPath
	
	::KEGroups::SelectionAssign $::KMProps::selectedEntity $GroupId $::KMProps::WinPath
	
	#Referescamos la ventana de grupos
	::KEGroups::InitBaseWindow
	
	#Cerramos la ventana de grupos
	destroy $::KEGroups::WinPath
	
	#Ponemos el foco en la ventana de propiedades
	focus $::KMProps::WinPath
	
	#Lo ejecuta demasiado pronto
	#set postcommand [$f.cGroups cget -postcommand]
	#set args [lrange $postcommand 1 end]
	#::KMProps::changeGroups $args
}

proc ::KMProps::getApplication { fullname } {
	
	set application [lindex [split $fullname "//"] 0]
	
	return $application
}

proc ::KMProps::getProps { fullname } {
	
	global KPriv
	
	set application [::KMProps::getApplication $fullname]
	
	#Miramos si hay alguna propiedad dada de alta (si la hay,obligatoriamente estará seleccionada)
	set props [::xmlutils::setXmlContainerIds "${application}//c.Properties"]
	
	#Miramos si es necesario filtrar por Constitutive Laws
	
	#Extraemos del fullname, el id del elemento pulsado
	set lSplit [split $fullname "//"]
	set lRange [lindex $lSplit 4] 
	set id [string range $lRange 2 end]
	
	#msg "$fullname \n $lSplit \n $lRange \n $id \n\n"
	
	set cLawsList {}
	set xPath "Kratos_KWords/ElementValidCLaws/Item\[@id='$id'\]"
	set active [::xmlutils::getAttribute $KPriv(xmlDocKKW) $xPath "active"]
	if { $active == 1 } {
		set cLawsList [split [::xmlutils::getAttribute $KPriv(xmlDocKKW) $xPath claws] ","]   
		
		if { $cLawsList != "" } {
		        
		        set FilterProps {} 
		        foreach idProp $props {
		        
		        set dv [::xmlutils::setXml "${application}//c.Properties//c.${idProp}//c.MainProperties//i.MatModel" dv]
		        
		        if {$dv in $cLawsList } {
		                
		                lappend FilterProps $idProp
		        }
		        }
		        return $FilterProps
		}
	}
	
	return $props
}

proc ::KMProps::acceptGroups { T idTemplate fullname item listT entityList fGroups} {
	
	set grupo $::KMProps::selGroup
	
	if { $grupo == "" } {
		WarnWin [= "You have to choose one group\n (you can create a new one pushing the button on the right)"]
	} else {
		
		#Primero comprobamos q el grupo aun no exista
		set id [::xmlutils::setXml "$fullname//c.$grupo" id]
		
		if { $id != "" } {
		        
		        WarnWin [= "This group it is already assigned to this condition."]
		} else {
		        
		        #Validamos que haya alguna propiedad seleccionada
		        if {[info exists ::KMProps::cmbProperty]} {
		        if {$::KMProps::cmbProperty == "" } {
		                WarnWin [= "You must define a Property before!"]
		                return ""
		        }
		        }
		        
		        #Comprobamos que el grupo no sea un AutoGroup sin entidades
		        ::KMProps::changeGroups $entityList $fGroups
		        if { $::KMProps::selGroup == "" } {
		                return ""
		        }
		        
		        set template [::KMProps::copyTemplate ${idTemplate} $fullname "$grupo" "Group"]
		        #Validamos que haya algún combo de "properties" en el template
		        if {[string match "*GCV=\"Properties*" $template]} {
		        
		                #Miramos si hay alguna propiedad dada de alta (si la hay,obligatoriamente estará seleccionada)
		                set props [::xmlutils::setXmlContainerIds "[::KMProps::getApplication $fullname]//c.Properties"]
		                if { [llength $props] < 1 } {
		                        WarnWin [= "You must define a Property before."]
		                        #::KMProps::deleteProps $T $itemSel $newPropertyName
		                        return  ""
		                }
		        }
		        
		        
		        #
		        # Ahora debemos actualizar todos los valores en el xml
		        #
		        set fBottom ${::KMProps::WinPath}.nb.fProp.fBottom
		        
		        if {[llength $listT] >= 1 } {
		        
		        #Lista de listas con formato {idContainer idItem1 idItem2...}
		        foreach listContainer $listT {
		                
		                #Si tiene como mínimo el container y un item entramos
		                if {[llength $listContainer] >= 2} {
		                
		                set idContainer [lindex $listContainer 0]
		                
		                #set id [::KMProps::getPropTemplate $idTemplate id $idContainer]
		                
		                #Recorremos los items
		                for {set i 1} { $i < [llength $listContainer] } {incr i} {
		                        
		                        set id [lindex $listContainer $i]
		                        
		                        #Cuando tengamos nodos ocultos o incompatibles la variable no existirá
		                        if {[info exists ::KMProps::cmb$id]} {
		                        
		                        set fullNombre "$fullname//c.$grupo//c.$idContainer//i.$id"
		                        
		                        set f "${fBottom}.nb.f${idContainer}.cmb${id}"
		                        
		                        if { [winfo exists $f] } {
		                                
		                                if { [$f cget -state] == "readonly" } {
		                                
		                                        set value [::xmlutils::getComboDv $f $fullNombre]
		                                } else {
		                                
		                                        set value [set ::KMProps::cmb$id]                                                                          
		                                }
		                                
		                                if {$id == "Vx" || $id == "Vy" || $id == "Vz"} {
		                                
		                                        set activeId "A[string range $id 1 1]"
		                                        set fullActive "$fullname//c.$grupo//c.Activation//i.$activeId"
		                                        set active [::xmlutils::setXml $fullActive dv "read"]
		                                        
		                                        if { $active == 0 } {
		                                                ::xmlutils::setXml $fullNombre state "write" "disabled"
		                                        } else {
		                                                ::xmlutils::setXml $fullNombre state "write" "normal"
		                                        }
		                                }
		                                #Comprobamos si el combo tiene una función asignada
		                                set function [::xmlutils::setXml $fullNombre function]
		                                if { $function != "" && [$f cget -state] == "disabled"} {
		                                        
		                                        ::xmlutils::setXml $fullNombre function "write" 1
		                                        ::xmlutils::setXml $fullNombre state "write" "disabled"
		                                } 
		                                
		                                #Guarda el nuevo valor en el xml
		                                ::xmlutils::setXml $fullNombre dv "write" $value
		                        }
		                        }
		                }
		                }
		        }
		        }

		        #Destruimos el frame inferior
		        ::KMProps::cancelBottom
		        
		        ::KMProps::refreshTree $T
		        
		        $T selection add $item
		        $T item expand $item
		}
	}
}

#
# Construye un frame con un tab por cada container que cuelgue de "item"
# y dentro de cada tab, etiquetas y combos para cada item (si los hay)
# 
proc ::KMProps::buildTabFrame { T item {class "Tab"} } {

	set f [::KMProps::iniFrameBottom]
	
	# Miramos los descendientes directos y si son container ponemos un tab por cada uno q tenga items
	set children [$T item children $item]

	set listTabs {}
	set listItems {}
	set acceptItems {}
	foreach itemChild $children {
		
		set fullname [DecodeName [$T item tag names $itemChild]]
		
		set nodeName [::xmlutils::getXmlNodeName $fullname]
		
		#Miramos si cada hijo es container o item
		if { $nodeName == "Container" } {
		        #Si no tiene items no agregamos el tab
		        if { [$T item numchildren $itemChild] > 0 } {
		        lappend listTabs $fullname
		        }
		        #msg "numchild:[$T item numchildren $itemChild]"
		} elseif { $nodeName == "Item" } {
		        lappend listItems $itemChild
		}
	}
	
	#Reseteamos la variable que nos indica si estamos en una propiedad
	set ::KMProps::propertyName ""
	
	# Si no tiene containers pero tiene items, utilizamos como tab el elemento padre seleccionado
	if { [llength $listTabs] == 0 && [llength $listItems] >= 1 } {
		
		set fullname [DecodeName [$T item tag names $item]]
		lappend listTabs $fullname
	} else {
		set listItems {}
	}
	
	if { [llength $listTabs] >= 1 } {
		
		set nb ${f}.nb
		grid [ttk::notebook $nb ] -row 0 -sticky ewn
		
		set i 0
		foreach fullname $listTabs {
		        
		        set id [::xmlutils::setXml $fullname id]
		        set pid [::xmlutils::setXml $fullname pid]                                
		        set help [::xmlutils::setXml $fullname help]
		        
		        #Para cada container declaramos un tab
		        set fTab ${nb}.f$id
		        
		        $nb add [ttk::labelframe $fTab -text "[= Properties]" -padding {10 10 10 10}] -text "[string range [= $pid] 0 20]"
		        
		        if {[llength $listItems] > 0 } {
		        #En este caso en realidad los nietos son los propios hijos del item pulsado
		        set nietos [$T item children $item]
		        } else {
		        #Miramos los hijos de cada container (nietos del pulsado originalmente)
		        set itemChild [lindex $children $i]
		        set nietos [$T item children $itemChild]
		        }
		        
		        ######################
		        # NOMBRE PROPIEDAD - Añadimos "Nombre de propiedad"
		        if {$class == "Property" && $i == 0 } {
		        
		        #Nombre de propiedad
		        grid [ttk::label $fTab.lblName -text "[= Property Name:]" ] \
		                -row 0 -column 0 -pady 5 -sticky nw -in $fTab
		        
		        grid [ttk::combobox $fTab.cmbPropertyName -state normal -textvariable "::KMProps::propertyName"] \
		                -row 0 -column 1 -padx 3 -pady 5 -sticky nw -in $fTab
		        
		        set ::KMProps::propertyName "[$T item text $item 0]"
		        
		        } elseif { $class == "Group" && $i == 0 } {
		        
		        ######################
		        #COMBO DE GRUPOS
		        set fullname [DecodeName [$T item tag names [$T item parent $item]]]
		        
		        set ::KMProps::selGroup "[$T item text $item 0]"
		        
		        set entityList [split [::xmlutils::setXml $fullname GiDEntity] ","]
		        
		        set filterGroups [::KMProps::getGroups $entityList $fullname]
		        
		        grid [ttk::label $fTab.lblName -text "[= Group:]" ] \
		                -row 0 -column 0 -pady 5 -sticky nw -in $fTab
		        
		        set fGroups $fTab.cGroups
		        grid [ttk::combobox $fGroups -state readonly -values "$filterGroups" -textvariable "::KMProps::selGroup"  -postcommand "::KMProps::changeGroups [list $entityList] $fGroups $fullname" -width 20] \
		                -row 0 -column 1 -pady 5 -sticky nw -in $fTab
		        
		        
		        
		        bind $fTab.cGroups <<ComboboxSelected>> "::KMProps::cmbChangeCheckGroups $fTab"
		        ######################
		        }
		        
		        set row 1
		        foreach nieto $nietos {
		        
		        set fullname [DecodeName [$T item tag names $nieto]]
		        set dv [::xmlutils::setXml $fullname dv]
		        #Comprobamos q sea un item
		        if { [::xmlutils::getXmlNodeName $fullname]  == "Item" } {
		                
		                set state [::xmlutils::setXml $fullname state]
		                
		                if {$state != "hidden"} {
		                
		                        set id [::xmlutils::setXml $fullname id]
		                        set pid [::xmlutils::setXml $fullname pid]
		                        set tooltip [::xmlutils::setXml $fullname help]
		                        set function [::xmlutils::setXml $fullname function]
								
		                        #Nos guardamos cada item para actualizar sus valores al final
		                        lappend acceptItems $nieto $fTab
		                        
		                        ## Es importante que ElemType sea el primer combo para actualizar el "dv" en el xml
		                        #if { $id == "ElemType" } {
		                        #        set comboList [::xmlutils::getXMLValues "$idContainer//$id" $idTemplate "" "" $dv]
		                        #} else {
		                        #        set comboList [::xmlutils::getXMLValues "$idContainer//$id" $idTemplate]
		                        #}
		                        
		                        set comboList [::xmlutils::getXMLValues "$fullname"]
		                        
		                        set CBState [::xmlutils::setXml $fullname CBState]
								if { $CBState == "normal" } {
									set values $comboList
									set comboList {}
								} else {
									set values {}
								}
								
		                        #Para cada item añadimos label y combo
		                        grid [ttk::label $fTab.lbl$id -text "${pid}:" ] \
		                                -row $row -column 0 -padx 3 -pady 5 -sticky nw -in $fTab
		                        
		                        if { [llength $comboList] > 0 } {
		                                
		                                if {$state != "disabled"} {
		                                set state "readonly"
		                                }
		                                grid [ttk::combobox $fTab.cmb$id -values $comboList -state $state \
		                                          -textvariable "::KMProps::cmb$id" ]\
		                                -row $row -column 1 -padx 3 -pady 5 -sticky ne -in $fTab
		                                
		                                ::xmlutils::setComboDv $fTab.cmb$id $fullname $dv
		                                #set selected [::xmlutils::getSelected $dv $comboList]
		                                #$fTab.cmb$id current $selected
		                                
		                                if {$id == "Ax" || $id == "Ay" || $id == "Az"} {
		                                ::KMProps::cmbDisable $fullname $f.nb
		                                bind $fTab.cmb$id <<ComboboxSelected>> "::KMProps::cmbDisable $fullname $f.nb"
		                                } elseif { $id == "ElemType" } {
		                                #En este caso se tendrá qué recargar si existe el combo de "Material Model"
		                                bind $fTab.cmb$id <<ComboboxSelected>> "::KMProps::cmbElemTypeChange $fTab.cmb$id $fullname"
		                                }
		                                
		                        } else {
		                                if {$state != "disabled"} {
		                                set state "normal"
		                                }
		                                
		                                grid [ttk::combobox $fTab.cmb$id -state $state -values $values -textvariable "::KMProps::cmb$id"] \
		                                -row $row -column 1 -padx 3 -pady 5 -sticky nw -in $fTab
		                                set ::KMProps::cmb$id $dv
		                        }
		                        
		                        if {$function != "" } {
		                                grid [ttk::button $fTab.funct$id -text "functions" -command "KFun::InitBaseWindow $fTab $id" -style TMenubutton.Toolbutton] \
		                                -row $row -column 2 -sticky nw  -pady 0 -padx 3 -in $fTab
		                                tooltip::tooltip $fTab.funct$id [= "Function manager"]
		                                set img [::WinUtils::GetImage "functions.gif"]
		                                if { $img != -1 } { $fTab.funct$id configure -image $img }
		                                
		                                grid [ttk::button $fTab.deleteFunct$id -text "delete" -command "::KMProps::unassignFunction $fTab $id $fullname" -style TMenubutton.Toolbutton] \
		                                -row $row -column 3 -sticky nw  -pady 0 -padx 2 -in $fTab
		                                tooltip::tooltip $fTab.deleteFunct$id [= "Unassign function"]
		                                set img [::WinUtils::GetImage "delete_icon.gif"]
		                                if { $img != -1 } { $fTab.deleteFunct$id configure -image $img }
		                                
		                        }
		                        
		                        #En el caso de ElemType actualizamos la variable de filtrado del Thickness
		                        if { $id == "ElemType" } {
		                                set ::KMProps::ElemTypeThickness $dv
		                        
		                        #En el caso del Thickness, no siempre se mostrará
		                        } elseif { $id == "Thickness" && ![::KMProps::showThickness]} {
		                                
		                                grid remove $fTab.lbl$id
		                                grid remove $fTab.cmb$id
		                        }
		                        
		                        incr row
		                }
		        }
		        }
		        incr i
		}
		# Si pulsan Esc también forzamos la salida del Tab
		bind $T <KeyPress> "if { %k == 27   } { ::KMProps::cancelBottom }"
		#bind $T <KeyPress> "if { %k == 13   } {  [list ::KMProps::acceptTabFrame $T $acceptItems $class $item] }"
		#bind $::KMProps::WinPath <KeyPress> "if { %k == 13   } { [list ::KMProps::acceptTabFrame $T $acceptItems $class $item] }"
		
		if { [llength $acceptItems] } {
		        grid [ttk::button $f.bPropOk -text "Ok"  -command "[list ::KMProps::acceptTabFrame $T $acceptItems $class $item]" ] \
		        -row 1 -column 0 -sticky sw  -pady 3 -padx 20  -in $f
		        tooltip::tooltip $f.bPropOk [= "Confirm values"]
		        
		        grid [ttk::button $f.bPropCancel -text "Cancel"  -command "::KMProps::cancelBottom" ] \
		        -row 1 -column 0 -sticky sw  -pady 3 -padx 100  -in $f
		        tooltip::tooltip $f.bPropCancel [= "Cancel assignation"]
		}
	}
}

#
# Deshabilita el widget de values si el flag activation está a 0
#
proc ::KMProps::cmbDisable { fullname f {id "" }} {
	
	if {$id == "" } {
		set iniIdEmpty ""
		set id [::xmlutils::setXml $fullname id]
	} else {
		set iniIdEmpty "$id"
	}
	
	set selCombo [set ::KMProps::cmb$id]
	
	set whatDisable [string map {A V} $id]
	
	set fullname [string map [list "Activation//i.$id" "Values//i.$whatDisable"] $fullname]
	
	#msg "state?[$f.fValues.cmb$whatDisable  "
	if {$selCombo == 0 } {
		
		if {[winfo exists $f.fValues.cmb$whatDisable]} {
		        $f.fValues.cmb$whatDisable configure -state disabled
		}
		if { $iniIdEmpty == "" } {
		        #set ::KMProps::cmb$whatDisable "Disabled"
		        ::xmlutils::setXml $fullname state "write" disabled
		        
		}
	} else {
		if {[winfo exists $f.fValues.cmb$whatDisable]} {
		        $f.fValues.cmb$whatDisable configure -state normal
		}
		if { $iniIdEmpty == "" } {
		        set dvChild [::xmlutils::setXml $fullname dv]
		        set ::KMProps::cmb$whatDisable $dvChild                
		        ::xmlutils::setXml $fullname state "write" normal
		}
	}
}

proc ::KMProps::getCmbWidth { comboList } {

	set width 15
	#Validamos el tamaño de los string del combo para ponerle uno o otro tamaño
	foreach c $comboList {
		if { [string length $c] > 12 } {
		        set width 25
		}
		if { [string length $c] > 20 } {
		        set width 32
		}
	}
	return $width
}                
#
# Contruye un combo dinámico en el item pulsado del árbol y posteriormente se elimina
#
proc ::KMProps::buildFrame { T item { type "props" } } {
	
	set fullname [DecodeName [$T item tag names $item]]
	set idFull [string map { "." "" "//" ""} $fullname]
	
	#Comprobamos que sea un item
	if { [::xmlutils::getXmlNodeName $fullname]  == "Item" } {

		set id [::xmlutils::setXml $fullname id]
		set pid [::xmlutils::setXml $fullname pid]
		set dv [::xmlutils::setXml $fullname dv]
		set tooltip [::xmlutils::setXml $fullname help]
		set comboList [::xmlutils::getXMLValues $fullname]
		
		set CBState [::xmlutils::setXml $fullname CBState]
		if { $CBState == "normal" } {
		        set values $comboList
		        set comboList {}
		} else {
		        set values {}
		}
		
		#---------------------------#---------------------------#
		# Configurar frame en función del XML
		#---------------------------#---------------------------#
		set bg "#F8F8F8"
		set f "$T.f$idFull"
		
		if { [winfo exists $f] } {
		        return
		}
		set f [frame "$T.f$idFull" -borderwidth 0 -background $bg]
		
		if { [llength $comboList] > 0 } {
		        
		        grid [ttk::combobox $f.cmb -values $comboList -state readonly -width [::KMProps::getCmbWidth $comboList] -textvariable "::KMProps::cmb$idFull"] \
		        -row 0 -column 0 -padx 3 -sticky ne -in $f
		        
		        ::xmlutils::setComboDv $f.cmb $fullname $dv
		        #set selected [::xmlutils::getSelected $dv $comboList]
		        #$f.cmb current $selected
		        
		        bind $f.cmb <<ComboboxSelected>> "::KMProps::cmbSelectChange $item $T 1 current"
		} else {
		        
		        grid [ttk::combobox $f.cmb -state normal -values $values -textvariable "::KMProps::cmb$idFull"] \
		        -row 0 -column 0 -padx 3 -sticky nw -in $f
		        
		        set ::KMProps::cmb$idFull $dv
		        
		        bind $f.cmb <Leave> "::KMProps::cmbSelectChange $item $T 0 current"
		        bind $f.cmb <FocusOut> "::KMProps::cmbSelectChange $item $T 1 current"
		        bind $f.cmb <Escape> "::KMProps::cmbCancel $item $T"
		}
		# Si pulsan intro o Esc también forzamos la salida del combo
		bind $f.cmb <KeyPress> "if { %k == 13  } { ::KMProps::cmbSelectChange $item $T 1 current}"
		bind $T <KeyPress> "if { %k == 27   } { ::KMProps::cmbCancel $item $T }"
		
		tooltip::tooltip $f.cmb "$tooltip"
		
		return $f
		
	} else {
		
		return ""
	}
}

proc ::KMProps::cmbCancel { item T  } {
	
	set fullname [DecodeName [$T item tag names $item]]
	set idFull [string map { "." "" "//" ""} $fullname]
	
	set pid [::xmlutils::setXml $fullname pid]
	set id [::xmlutils::setXml $fullname id]
	set dv [::xmlutils::setXml $fullname dv]
	set dvText [::xmlutils::setXml $fullname dvText]
	
	set ::KMProps::cmb$idFull "$dv"
	
	set f "$T.f$idFull"
	if {[winfo exists $f]} {
		destroy $f
	}
	#Miramos si tiene algun estilo especial
	$T item style set $item C0 styAnyRead
	if { [::xmlutils::setXml $fullname style] == "*" } {
		$T item element configure $item C0 elemTxtRead -text "$pid* : $dvText"
	} else {
		$T item element configure $item C0 elemTxtRead -text "$pid: $dvText"
	}
	
	set ::KMProps::lastSelected {}
}
#$id $T 1 current
proc ::KMProps::cmbSelectChange { item T {remove 1} {selectVal "current"} } {
	
	#msg "cmbselectChange $remove"
	set refresh 0
	set fullname [DecodeName [$T item tag names $item]]
	set idFull [string map { "." "" "//" ""} $fullname]
	
	set id [::xmlutils::setXml $fullname id]
		
	set comboState [::xmlutils::getComboBoxState $fullname]
	if { $comboState == "normal" } {
		
		#Si han pulsado el combo en modo editar no tiene que hacer nada (antes desaparecía)
		set f "${T}.f${idFull}.cmb"
		if { [winfo exists $f] && [$f state] == "pressed" } {
			return
		}
		
		set selCombo [set ::KMProps::cmb$idFull]
		set selComboText $selCombo
		if {$selectVal != "current"} {
		        set ::KMProps::lastSelected {}
		}
	} else {
		
		if {$selectVal == "current"} {
		        
		        set f "$T.f$idFull.cmb"
		        set selCombo [::xmlutils::getComboDv $f $fullname]
		        set selComboText [::xmlutils::getComboDv $f $fullname "text"]
		        set ::KMProps::lastSelected [list $item $selCombo $selComboText]
		        
		} else {  
		        
		        set selCombo [lindex $::KMProps::lastSelected 1]
		        set selComboText [lindex $::KMProps::lastSelected 2]
		        set ::KMProps::lastSelected {}
		}
	}
	
	if {$id == "Ax" || $id == "Ay" || $id == "Az"} { 
		
		set whatDisable [string map {A V} $id]
		set abuelo [$T item parent [$T item parent $item]]
		foreach itemChild [$T item descendants $abuelo] {
		        
		        set fullnameChild [DecodeName [$T item tag names $itemChild]]
		        set idChild [::xmlutils::setXml $fullnameChild id]
		        if { $idChild == $whatDisable } {
		        if {$selCombo == 0 } {
		                $T item enabled $itemChild 0
		                $T item element configure $itemChild C0 elemTxtRead -fill { gray }
		                #$T item element configure $item C0 elemTxtRead -text "[$T item text $itemChild 0]: Disabled" -fill { gray }
		                ::xmlutils::setXml $fullnameChild state "write" disabled
		        } else {
		                $T item enabled $itemChild 1
		                set dvChild [::xmlutils::setXml $fullnameChild dv]
		                $T item element configure $itemChild C0 elemTxtRead -fill { black }
		                ::xmlutils::setXml $fullnameChild state "write" normal
		        }
		        }
		}
	} elseif { $id == "ElemType" } {
		
		global KPriv        
		set matFullname [string map {"ElemType" "MatModel"} $fullname]
		
		::xmlutils::setXml $matFullname dv
		
		##Si se cambia el ElementType hay que cambiar el dv de Material Model
		set xpath "Kratos_KWords/ElementCLaws/Item\[@id='ElementTypes'\]"
		::xmlutils::getAttribute $KPriv(xmlDocKKW) $xpath dv "$selCombo"
		set values [::xmlutils::getXMLValues "$matFullname" "" "" "" "NoElementFilter"]
		
		set dvMat [::xmlutils::setXml $matFullname dv]
		
		if { $dvMat ni $values } {
		        ::xmlutils::setXml $matFullname dv "write" [lindex $values 0]
		        #msg "$$matFullname \n $selCombo --> $values dvMat:$dvMat dv:[::xmlutils::setXml $matFullname dv] values0:[lindex $values 0]"
		        
		}
		#Hacemos un refresh para refrescar los cambios, tanto en matModel com en Thickness
		set refresh 1
		
		#
		#
		#
		#set values [::xmlutils::getXMLValues "$matFullname" "" "" "" "NoElementFilter"]
		#
		#if { $dvMat ni $values } {
		#        ::xmlutils::setXml $matFullname dv [lindex $values 0]
		#        ::KMProps::refreshTree $T
		#}
	}
	
	set pid [::xmlutils::setXml $fullname pid]
	
	#msg "id:$id rem:$remove sel:$selCombo"
	if { $remove } {
		
		#Eliminamos el frame
		set f "$T.f$idFull"
		if {[winfo exists $f]} {
		        destroy $f
		}
		
		#Miramos si tiene algun estilo especial
		$T item style set $item C0 styAnyRead
		if { [::xmlutils::setXml $fullname style] == "*" } {
		        $T item element configure $item C0 elemTxtRead -text "$pid* : $selComboText"
		} else {
		        $T item element configure $item C0 elemTxtRead -text "$pid: $selComboText"
		}
	}
	
	#Vuelve a poner la imagen por si no recargamos el árbol
	set icon [::xmlutils::setXml $fullname icon]
	set imagen [::WinUtils::GetImage $icon]
	if { $imagen != -1 } {
		$T item image $item C0 $imagen
	}
	
	#Miramos si ha cambiado para actualizar en consecuencia
	if { [::xmlutils::setXml $fullname dv] != $selCombo } {
		
		#Guarda el nuevo valor en el xml
		::xmlutils::setXml $fullname dv "write" $selCombo
		
		#Si el combo tiene un class especial habrá q reconstruir el árbol
		set clase [::xmlutils::setXml $fullname class]
		if { $clase != "" } {
		        
		        #Comprueba si el cambio va a afectar a la estructura del arbol
		        ::KMProps::specialComboAction $T $clase $selCombo $item $id
		        
		        #Volvemos a cargar el árbol de propiedades
		        set refresh 1
		}
	}
	
	if {$refresh} {
		::KMProps::refreshTree $T
	}
}

#
# Si el atributo class de un nodo requiere algún cambio en el xml
# se prepara todo para q el refresh tree cambie el treeCtl
#
proc ::KMProps::specialComboAction { T clase selCombo item id } {
	
	global KPriv 
	if { $clase == "application" } {
		
		if {$selCombo == "Yes"} {
		        
		        #Este combo solo permite una aplicación activa al mismo tiempo
		        #Así que desactivamos el resto
		        set padre [$T item parent $item]
		        set fullParent [DecodeName [$T item tag names $padre]]
		        set parentId [::xmlutils::setXml $fullParent id]
		        set applications [::xmlutils::setXmlContainerIds $fullParent "Item"]
		        foreach aplicId $applications {
		        if { $aplicId != $id } {
		                ::xmlutils::setXml "${fullParent}//i.${aplicId}" dv "write" "No"
		                ::xmlutils::setXml "$aplicId" state "write" "hiddenAll"
		        }
		        }
		        
		        #Activamos la que corresponda
		        if { $id == "FluidStructureInteraction" } {
		        
		        #Activamos Fluidos y Estructuras
		        ::xmlutils::setXml "StructuralAnalysis" state "write" "normal"
		        ::xmlutils::setXml "StructuralAnalysis" open "write" 1
		        ::xmlutils::setXml "Fluid" state "write" "normal"
		        ::xmlutils::setXml "Fluid" open "write" 1
		        
		        #Forzamos el tipo de fluido a "Incompressible"
		        set fluidType [::xmlutils::setXml "Fluid//c.AnalysisData//i.FluidType" dv]
		        if { $fluidType == "Compressible"} {
		                ::xmlutils::setXml "Fluid//c.AnalysisData//i.FluidType" dv "write" "Incompressible"
		        }
		        } else {
		        #Activamos la aplicación seleccionada
		        ::xmlutils::setXml "$id" state "write" "normal"
		        #La desplegamos
		        ::xmlutils::setXml "$id" open "write" 1
		        }
		        
		} else {
		        
		        if { $id == "FluidStructureInteraction" } {
		        #if{ [::xmlutils::setXml "${fullParent}//i.${aplicId}" dv] } {}
		        ::xmlutils::setXml "StructuralAnalysis" state "write" "hiddenAll"
		        ::xmlutils::setXml "Fluid" state "write" "hiddenAll"
		        } else {
		        ::xmlutils::setXml "$id" state "write" "hiddenAll"
		        }
		}
	}

	foreach var $::KMProps::visibilityVars {
		
		if {$var == $clase} {
		        
		        #Caso general
		        if { [set ::KMProps::$var] != "$selCombo" } {
		        set ::KMProps::$var $selCombo
		        #msg "var$var"
		        }
		        
		        #Casos especiales
		        #if { $var == "fluidType" } {
		        #        if {$selCombo == "Compressible" } {
		        #                WarnWin [= "Compressible Fluids are not still available."]           
		        #                set fullname [DecodeName [$T item tag names $item]]
		        #                ::xmlutils::setXml "$fullname" dv "write" "Incompressible"
		        #        } 
		        
		        #} 
		}
	}
	return ""
}

proc ::KMProps::acceptTabFrame { T listItems class {itemSel ""}} {
	
	#
	# Actualiza el resto de valores del xml con el contenido de cada combo
	#
	foreach {item f} $listItems {
		
		set fullname [DecodeName [$T item tag names $item]]
		set idFull [string map { "." "" "//" ""} $fullname]
		set id [::xmlutils::setXml $fullname id]
		set pid [::xmlutils::setXml $fullname pid]
		set estado [::xmlutils::setXml $fullname state]
		
		
		#msg "id: $id estado:[$f.cmb$id cget -state]"
		
		set cmbState [$f.cmb$id cget -state]
		
		if { $cmbState == "readonly" } {
		        
		        set selCombo [::xmlutils::getComboDv "$f.cmb$id" $fullname]
		        set selComboText [::xmlutils::getComboDv "$f.cmb$id" $fullname "text"]
		} else {
		        
		        set selCombo [set ::KMProps::cmb$id]
		        set selComboText $selCombo 
		        
		}
		#$cmbState == "readonly"
		
		
		set function [::xmlutils::setXml $fullname function]
		if { $function != "" } {
		        if { $cmbState == "disabled" } {
		                ::xmlutils::setXml $fullname state "write" "disabled"
		                ::xmlutils::setXml $fullname function "write" 1
		        } else {
		                ::xmlutils::setXml $fullname state "write" "normal"
		        }
		}
		
		#Guarda el nuevo valor en el xml
		::xmlutils::setXml $fullname dv "write" $selCombo
		$T item style set $item C0 styAnyRead
		
		#Al final ya se hace un refresh
		#if {$cmbState == "disabled"} {
		#  $T item element configure $item C0 elemTxtRead -fill { gray }
		#} else {
		#  $T item element configure $item C0 elemTxtRead -text "$pid: $selComboText" -fill { black }
		#}
		
		#Guarda el nuevo valor en el xml (si está desabilitado deja el que había)
		#if {$cmbState == "disabled"} {
		#  ::xmlutils::setXml $fullname dv "write" $selCombo
		#  ::xmlutils::setXml $fullname state "write" "disabled"
		#} else {
		#  ::xmlutils::setXml $fullname state "write" "normal"
		#}
		
		
		#Si el combo era especial se tendrá q reconstruir el árbol
		set clase [::xmlutils::setXml $fullname class]
		
		#Comprueba si el cambio va a afectar a la estructura del arbol
		::KMProps::specialComboAction $T $clase $selCombo $item $id
	}
	
	#
	# Caso especial para las propiedades
	#
	if {$class == "Property" && $itemSel != "" } {
		
		set newPropertyName $::KMProps::propertyName
		
		#Validamos que la propiedad no tenga carácteres extraños
		set newPropertyName [::KUtils::parseTreeStr $newPropertyName]
		if { $newPropertyName == -1 } {
		        
		        WarnWin [= "You can't use some reservate chars like:\n  :   /   $   .   \\  %  "]
		        set ::KMProps::propertyName ""
		        return
		}
		
		#Comprobamos q el nombre no sea vacío
		if { $newPropertyName == "" } {
		        WarnWin [= "The property name can not be empty"]
		        return
		}
		
		#Comprobamos q si la propiedad ha cambiado no exista el nuevo nombre
		set fullname [DecodeName [$T item tag names $itemSel]]
		set oldId [::xmlutils::setXml "$fullname" id]
		set repeatId [::xmlutils::setXml "[::KMProps::getApplication $fullname]//c.Properties//c.$newPropertyName" id]
		
		if { $repeatId != "" && $oldId != $newPropertyName } {
		        WarnWin [= "This property name already exist."]
		        return
		}
		
		#Si todo es correcto, renombramos la propiedad en todo el árbol y luego lo recargamos
		#(deleteProps hace un rename si llega el 3er argumento "newPropertyName")
		::KMProps::deleteProps $T $itemSel $newPropertyName
		
		::xmlutils::setXml $fullname pid "write" $newPropertyName
		::xmlutils::setXml $fullname id "write" $newPropertyName
		
	} elseif { $class == "Group" } {

		#
		# Caso especial para los Grupos
		#
		set fullname [DecodeName [$T item tag names $itemSel]]
		
		::xmlutils::setXml $fullname pid "write" $::KMProps::selGroup
		::xmlutils::setXml $fullname id "write" $::KMProps::selGroup
	}
	
	#Volvemos a cargar el árbol para q el path de los items sea correcto
	::KMProps::refreshTree $T
	
	::KMProps::cancelBottom
}



proc ::KMProps::editTag { T item fullname newtext { newPath "" }  } {
	
	global KPriv
	
	if { $newPath == "" } {
		set parts [::KMProps::split2 $fullname //]
		lset parts end $newtext
		set newPath [join $parts //]
	}
	
	if { $fullname != $newPath } {
		
		#msg "$fullname != $newPath"
		
		#Cambiar nombre en el árbol
		$T item tag remove $item [list names [$T item tag names $item]]
		
		$T item tag add $item [EncodeName $newPath]
		$T item element configure $item C0 elemTxtRead -text $newtext
		
		#Cambiar nombre en el XML
		#::xmlutils::setXml $fullname id "write" $newtext
		
		return $newPath
		
	} else {
		
		return $fullname
	}
	
}

#separator can be a multicharacter, like //
proc ::KMProps::split2 { x separator } {
	set splitchar \uFFFF ;#forbidden unicode character, x must never contain it
	return [split [string map "$separator $splitchar" $x] $splitchar]
}

proc ::KMProps::BeginEditGroups { T } {
	set I [$T item id active]
	set C 0
	set E elemTxtRead
	::TreeCtrl::FileListEdit $T $I $C $E
}

proc ::KMProps::ReceiveDragGroups { T dragged_list dst } {
	
	#set dstname [DecodeName [$T item tag names $dst]]
	#if { [IsItemFolderOfLayers $T $dst] } {
	##folder of layers
	#set cmd [GetOldAndNewLayernamesRecursive $T $dragged_list $dstname ""]
	#GiD_Process Layers {*}$cmd escape
	#} else {
	#WarnWin "$dstname is not a folder"
	#}
}

proc ::KMProps::WriteGeomToVar { w what geomname {InitComm ""}} {
	variable WinLayout
	global GidPriv

	set trans 1
	update idletasks
	if { $WinLayout == "RIGHT" } {
		# maintain all dimensions except the new width
		set prevgeom [lindex $GidPriv($geomname) 1]
		foreach {width height x y} [split $prevgeom x+] break
		set width [winfo width $w]
	} else {
		# "OUTSIDE"
		set width [winfo width $w]
		set height [winfo height $w]
		set x [winfo x $w]
		set y [winfo y $w]
	}
	set geom ${width}x${height}+${x}+${y}
	set GidPriv($geomname) "$what $geom $trans $InitComm"
	#  msg $GidPriv($geomname)
}

proc ::KMProps::OpenWindowOutside { w } {
	
	if { [winfo exists $w] } {
		::KMProps::CloseWindowOutside $w
	}

	# Init the window
	set title [= "Project properties"]
	InitWindow $w $title KMPropsWindowGeom ::KMProps::InitBaseWindow
	
	wm protocol $w WM_DELETE_WINDOW "[list ::KMProps::refreshTree "" 1];destroy $w"
	
	return $w
}

proc ::KMProps::CloseWindowOutside { w } {
	
	::KMProps::refreshTree "" 1
	
	destroy $w
}

proc ::KMProps::InsertNewProp { node id T {parent ""} {parentitem root} {childs true} {state "normal"} {open 0}} {
	
	set propName [$node getAttribute pid ""]
	
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
	
	if { $id == "" } {
		WarnWin [$node getAttribute id ""]
		#set id [$node getAttribute gpid ""]
	}
	
	if { $parent != "" } {
		set fullname $parent$id
	} else {
		set fullname $id
	}
	
	
	if { $childs } {
		
		set item [$T item create -button yes -tags [EncodeName $fullname] -open $open]
		$T item lastchild $parentitem $item
		$T item style set $item C0 styAnyRead
		$T item element configure $item C0 elemTxtRead -text "$propName"
		
	} else {
		
		set item [$T item create -button no -tags [EncodeName $fullname]]
		$T item lastchild $parentitem $item
		
		if { [string range $id 0 1] == "i." } {
		        
		        #set f [::KMProps::buildFrame $T $item]
		        
		        $T item style set $item C0 styFrame                                                
		        #$T item element configure $item C0 eWindow -window $f
		        
		        set dv [::xmlutils::setXml $fullname dvText]
		        
		        $T item style set $item C0 styAnyRead
		        if {$state != "disabled"} {
		        
		        #Miramos si tiene algun estilo especial
		        if { [::xmlutils::setXml $fullname style] == "*" } {
		                $T item element configure $item C0 elemTxtRead -text "$propName* : $dv"
		        } else {
		                $T item element configure $item C0 elemTxtRead -text "$propName: $dv"
		        }
		        
		        } else {
		        $T item element configure $item C0 elemTxtRead -text "$propName: $dv" -fill {gray }
		        }
		        
		        
		        
		} elseif { [string range $id 0 1] == "c." } {                                                                                                                   
		        
		        $T item style set $item C0 styAnyRead
		        $T item element configure $item C0 elemTxtRead -text "$propName"
		}
	}
	
	# Consultamos el icono en el xml, y si existe en nuestro directorio se lo añadimos                
	set icon [::xmlutils::setXml $fullname icon]
	set imagen [::WinUtils::GetImage $icon]
	if { $imagen != -1 } {
		$T item image $item C0 $imagen
	} else {
		# Si no encuentra la imagen no mostramos ninguna (descomentando mostramos defaultIcon)
		#$T item image $item C0 [::WinUtils::GetImage "defaultIcon.gif"]
	}
	
	#Miramos si el item tiene que estar a disabled (viene del proc ::KMProps::stateNode)
	if {$state == "disabled"} {
		$T item enabled $item 0
	}
	
	return $item
}

proc ::KMProps::changeTypeEvent {T fullname item {changeTo ""} } {
	
	if { $changeTo == "" } {
		set tipo [ set ::KMProps::type${item}]
	} else {
		set tipo "$changeTo"
	}
	
	::xmlutils::setXml $fullname type "write" $tipo
	
	#Cambiamos también el tipo de todos los item seleccionados
	set items [$T selection get]
	
	for {set i 0} { $i < [llength $items] } {incr i} {
		
		set item [lindex $items $i]
		set ::KMProps::type$item $tipo
		
		set fullname [DecodeName [$T item tag names $item]]
		
		::xmlutils::setXml $fullname type "write" $tipo
	}
}

proc ::KMProps::randomColor { } {

	set r [expr {int (255 * rand())}]
	set g [expr {int (255 * rand())}]
	set b [expr {int (255 * rand())}]                
	
	return [format "\#%02x%02x%02x" $r $g $b]
}

proc ::KMProps::MenuContextual { T x y } {
	
	set w $T.menucontextual
	if { [winfo exists $w] } {
		destroy $w
	}

	menu $w

	#$w add command -label [= "New Group#C#layer"] -command [list ::KMProps::CreateNewGroupId $T]
	
	set item [$T selection get]
	#msg "item$item"
	if { $item != "" } {
		set fullname [DecodeName [$T item tag names $item]]
		set class [::xmlutils::setXml $fullname class]
		
		if {$class == "Groups" || $class == "Properties"} {
		        $w add command -label [= "New"] -command [list ::KMProps::DoubleClickTree 0 0 $T $item] -state normal
		        $w add separator
		}

		if {$class == "Group"} {
		        
		        $w add command -label [= "Delete"] -command [list ::KMProps::deleteGroupCondition $T $item] -state normal
		        
		} elseif {$class == "Property"} {
		        
		        $w add command -label [= "Delete"] -command [list ::KMProps::deleteProps $T $item] -state normal
		        
		} else {                
		        $w add command -label [= "Delete"] -state disabled
		}
		
		
		if {$class == "Group" || $class == "Property"} {
		        $w add command -label [= "Edit"] -command [list ::KMProps::DoubleClickTree 0 0 $T $item] -state normal
		} else {
		        $w add command -label [= "Edit"] -state disabled
		}
		
		$w add separator
	} else {
		set item "root"                
	}                
	
	
	
	$w add command -label [= "Collapse All"] -command [list $T collapse -recurse "$item"] -state normal
	$w add command -label [= "Expand All"] -command [list $T expand -recurse "$item"] -state normal
	
	set x [expr [winfo rootx $T]+$x+2]
	set y [expr [winfo rooty $T]+$y]
	GiD_PopupMenu $w $x $y
}

#
# Borrar el grupo de la condición preguntando si lo queremos borrar completamente
#
proc ::KMProps::deleteGroupCondition { T item } {

	set fullname [DecodeName [$T item tag names $item]]
	set GroupId [$T item text $item 0]
	
	set aviso "Removing group $GroupId. Please, choose the properly option:\n\n\n Yes: Desassign this group (recomended)\n\n No: Complete group removing (with all its descendants).\n\n Cancel: Keep the group assigned."
	set answer [::WinUtils::confirmBox "." "$aviso" "yesnocancel"]
	if { $answer == "yes" } {
		
		#Desasigna de la gemoetría el item seleccionado
		#::KEGroups::UnAssignCondition $GroupId
		
		#Elimina el grupo del xml
		::xmlutils::unsetXml [DecodeName [$T item tag names $item]]
		
		#Elimina el grupo del árbol
		::KMProps::deleteItem $T $item
		
	} elseif {$answer == "no" } {
		
		#Consultamos el grupo a eliminar
		set GroupId [$T item text $item 0]
		
		set visibleGroups "[winfo exists $::KEGroups::WinPath]"
		#Iniciamos la ventana de grupos
		::KEGroups::InitBaseWindow 
		
		#Seleccionamos el item correspondiente de grupos
		set TG $::KEGroups::TreePath
		
		$TG selection clear
		
		set grupos [$TG item descendants "root"]
		foreach grup $grupos {
		        set gId [$TG item text $grup 0]
		        if { $gId == $GroupId } {
		        $TG selection add $grup
		        }
		}
		
		#Borramos recursivamente el grupo seleccionado y sus hijos
		::KEGroups::DeleteGroupsId $TG
		
		#Esto ya se hace desde Grupos y para todo el árbol de propiedades
		#Elimina el grupo del xml
		#::xmlutils::unsetXml [DecodeName [$T item tag names $item]]
		
		#Elimina el grupo del árbol
		#$T item delete $item
		
		#Si la ventana no estaba visible la volvemos a ocultar (destruir)
		if { !$visibleGroups } {
		        
		        destroy $::KEGroups::WinPath                 
		} else {
		        focus $T
		}
		
	} else {
	}
}

#
# Elimina todas las asignaciones de grupos que existan en el árbol de propiedades
#
proc ::KMProps::deleteGroups { GroupId {editName ""}} {
	
	#Primero hay que asegurarse de que exista el árbol
	if { [winfo exists $::KMProps::WinPath] } {
		
		#Si está el frameBottom activo nos lo cargamos
		::KMProps::cancelBottom
		
		set T $::KMProps::TreePropsPath
		
		# Recorremos tooodo el árbol
		set grupos [$T item descendants "root"]
		foreach item $grupos {
		        
		        #Comprobamos q el item aun exista
		        if { [$T item id $item] != "" } {
		        set itemId [$T item text $item 0]
		        if { $itemId == $GroupId } {
		                
		                set fullname [DecodeName [$T item tag names $item]]                
		                if {$editName == ""} {
		                #Elimina el grupo del xml
		                ::xmlutils::unsetXml $fullname
		                
		                #Elimina el grupo del árbol
		                ::KMProps::deleteItem $T $item
		                } else {
		                #Renombra el nodo en el xml
		                ::xmlutils::setXml $fullname pid "write" $editName
		                ::xmlutils::setXml $fullname id "write" $editName
		                }
		                
		        }
		        }
		}
		#Cierra la ventana de propiedades y la vuelve a cargar para q se renombren
		# los items automáticamente
		::KMProps::refreshTree $T
		
	} else {
		
		global KPriv
		
		#Si no está abierta la ventana, deberemos recorrer el xml a mano
		set nodes [$KPriv(xml) getElementsByTagName "Container"]
		
		foreach node $nodes {
		        
		        catch {
		        #Si encuentra el grupo intenta borrar el nodo y su descendencia
		        if { [$node getAttribute id ""] == "$GroupId" } {
		                if { [$node getAttribute class ""] == "Group" } {
		                
		                if {$editName == ""} {
		                        #Borra el grupo y su descencencia
		                        $node delete
		                } else {
		                        # Rename the node
		                        $node setAttribute pid $editName
		                        $node setAttribute id $editName

		                }
		                }
		        }
		        }
		}
	}
}

proc ::KMProps::deleteItem { T item } {
	
	#Comprobamos q el item exista
	if { [$T item id $item] != "" } {
		
		#msg "T:$T item$item"
		
		#Si el padre se queda sin hijos le quita la imagen
		set padre [$T item parent $item]
		
		#Elimina el item del árbol
		$T item delete $item 
		
		#Si se ha borrado bien, miramos si el padre sigue teniendo hijos
		if { [$T item numchildren $padre] == 0 } {
		        #Si no quedan hijos eliminamos el botón de desplegar
		        $T item configure $padre -button 0
		}
	}
}

#
# Elimina todas las asignaciones de grupos que existan en el árbol de propiedades
#
proc ::KMProps::deleteProps { T itemProp {editName ""}} {

	set propsToDelete {}
	
	# Recorremos tooodo el árbol
	set items [$T item descendants "root"]
	foreach item $items {
		
		#Comprobamos q el item aun exista
		if { [$T item id $item] != "" } {
		        
		        set fullname [DecodeName [$T item tag names $item]]
		        set itemDv [::xmlutils::setXml $fullname dv]
		        set propId [$T item text $itemProp 0]
		        
		        if { $itemDv == $propId } {
		        
		        
		        if {[::xmlutils::setXml $fullname GCV] == "Properties"} {
		                if {$editName == ""} {
		                
		                lappend propsToDelete $item
		                } else {
		                #Renombra el nodo en el xml
		                ::xmlutils::setXml $fullname dv "write" $editName
		                #::xmlutils::setXml $fullname id "write" $editName
		                }
		        }
		        }
		}
	}
	if { $editName == "" } {
		if { [llength $propsToDelete] > 0 } {
		        set aviso [= "This property it is assigned to %s groups.\nDo you want to delete the property and all this groups?" [llength $propsToDelete]]
		        set answer [::WinUtils::confirmBox "." "$aviso"]
		        if { $answer == "ok" } {
		        
		        #Borramos todos los grupos
		        foreach item $propsToDelete {        
		                
		                #No hay que borrar la propiedad, sino el grupo que la incluye
		                set itemGroup [$T item parent $item]
		                set fullname [DecodeName [$T item tag names $itemGroup]]
		                if {[::xmlutils::setXml $fullname class] == "Group"} {
		                
		                #Elimina el grupo del xml
		                ::xmlutils::unsetXml $fullname
		                
		                #Elimina el grupo del árbol (ya no hace falta, se hace refresh)
		                }
		        }
		        
		        #Borramos la propiedad
		        set fullProp [DecodeName [$T item tag names $itemProp]]
		        ::xmlutils::unsetXml $fullProp
		        }
		} else {
		        set aviso "Do you want to remove the property '[$T item text $itemProp 0]'?" 
		        set answer [::WinUtils::confirmBox "." "$aviso"]
		        if { $answer == "ok" } {
		        #Borramos la propiedad
		        set fullProp [DecodeName [$T item tag names $itemProp]]
		        ::xmlutils::unsetXml $fullProp
		        }
		}
	}
	
	#Cierra la ventana de propiedades y la vuelve a cargar para q se renombren
	# los items automáticamente
	::KMProps::refreshTree $T
}

#
# Elimina todas las asignaciones de grupos que existan en el árbol de propiedades
#
proc ::KMProps::chekMaterials { materialId {editName ""} } {
	
	global KPriv
	
	#Recorremos el xml
	set nodes [$KPriv(xml) getElementsByTagName "Item"]
	
	foreach node $nodes {
		
		#Si encuentra el grupo intenta borrar el nodo y su descendencia
		if { [$node getAttribute GCV ""] == "Materials" } {
		        if { [$node getAttribute dv ""] == $materialId } {
		        
		        $node setAttribute dv $editName
		        }
		}
	}
	
	::KMProps::refreshTree
}

#
#
#
proc ::KMProps::checkFunctions { functionId {newName ""} {action "execute"}} {

	global KPriv
	set numChanges 0
	
	#Si no está abierta la ventana, deberemos recorrer el xml a mano
	set nodes [$KPriv(xml) getElementsByTagName "Item"]
	
	foreach node $nodes {
		
		#Si encuentra el grupo intenta borrar el nodo y su descendencia
		if { [$node getAttribute function ""] == 1 } {
		        if { [$node getAttribute dv ""] == "$functionId" } {
		                
		                if { $action == "execute" } {
		                        if {$newName == ""} {
		                                #Resetea la función
		                                $node setAttribute state ""
		                                $node setAttribute function 0
		                                $node setAttribute dv "0.0"
		                        } else {
		                                # Rename the function
		                                $node setAttribute dv $newName
		                        }
		                }
		                set numChanges [expr $numChanges + 1]
		        }
		}
	}
	if { $numChanges > 0 && $action == "execute"} {
		if { [winfo exists $::KMProps::WinPath] } {
		        ::KMProps::cancelBottom
		        ::KMProps::refreshTree
		}
	}
	return $numChanges
}
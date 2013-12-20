###############################################################################
#
#	NAME: projectSettings.tcl
#
#	PURPOSE: KRATOS PROJECT SETTINGS
#
#	QUANTECH ATZ-DEVELOPMENT DEPARTMENT
#
#	AUTHORS : Luis C.A.  => LCA
#                 G. Socorro => GSM
#
#	CREATED AT: 08/07/2010
#
#	HISTORY:
#
#       0.3- 19/12/13- GSM, modify the proc updateKratosIniFile to include the widget state attribute
#       0.2- 15/07/13- GSM, update to use translation, disable kratos_debug.ini file
#	0.1- 08/07/10- LCA, create the base source code
#
###############################################################################

package provide kps 1.0 


# Create a base namespace kps (KRATOS PROJECT SETTINGS)
namespace eval ::kps:: {
	
	# Window path
	variable winpath
}

proc ::kps::Init {} {
	
	variable winpath

	set winpath ".gid.settingWin"
	
}

proc ::kps::InitSettingWindow {{w .gid.settingWin}} {
	
	# Init some properties
	::kps::Init

	# create report window
	set result [ ::kps::CreateSettingsWindow $w]
	
	return $result
	
}

proc ::kps::CreateSettingsWindow {w} {
	
	global KPriv
	
	# Init the window
	set title [= "Kratos global settings"]
	InitWindow $w $title ::kps::CreateSettingsWindowGeom ::kps::CreateSettingsWindow

	set f $w.psettings
	ttk::frame $f -padding {3 3 3 3}
	
	grid $f -sticky ewn
		
	# Lower buttons
	set def_back [$w cget -background]
	ttk::frame $w.frmButtons
	
	# Close button
	ttk::button $w.frmButtons.btnclose -text [= "Save and Close"] \
		-command "::kps::WindowbClose $w" \
		-underline 0 -width 14
	
	tooltip::tooltip $w.frmButtons.btnclose [= "Close the Kratos global settinh window"]
	
	# Geometry manager
	grid $w.frmButtons \
		-sticky ews \
		-columnspan 7
	
	grid $w.frmButtons.btnclose \
		-in $w.frmButtons \
		-row 0 -columnspan 7 \
		-sticky n \
		-padx 4m -pady 2m

	
	# For w
	wm protocol $w WM_DELETE_WINDOW "::kps::WindowbClose $w"

	grid rowconfigure $w 0 -weight 1
	grid columnconf $w 0 -weight 1
	
	# For button
	grid rowconfigure $w.frmButtons 0 -weight 1
	grid columnconf $w.frmButtons {0 1 2 3 4 5 6 7} -weight 1
	
	# Set the focus to the close button   
	focus $w.frmButtons.btnclose
	
	# Binding
	bind $w <Alt-c> "tkButtonInvoke $w.frmButtons.btnclose"
	bind $w <Escape> "tkButtonInvoke $w.frmButtons.btnclose"
	
	::kps::updateKratosIniFile $f
}

#
# Crea el archivo Kratos.ini, si aun no existe e inicializa
#  las variables globales KPriv(xmlIni)...
#
proc ::kps::initializeConfigFile { } {
	
    global KPriv
    
    set file [::kfiles::giveConfigFile]
    # wa "file:$file"
    
    # Open the spd default => the block GlobalProjectSettings
    
    set xmlFileDefault "$KPriv(dir)/kratos_default.spd"
    set xmlArray [::xmlutils::openFile "." "$xmlFileDefault"]
    set xmlDef [lindex $xmlArray 0]
    set xmlDoc [lindex $xmlArray 2]

    set nodeDefSettings [$xmlDef selectNodes "/Kratos_Data/GlobalProjectSettings"]
	
    # Create the file in case that the file does not exist 
    if { ![file exists $file] } {
	set fileid [open $file w]
	
	puts $fileid "<KratosConfig>"
	puts $fileid "[$nodeDefSettings asXML]"
	puts $fileid "</KratosConfig>"
	
	close $fileid
    }
	
    if { $KPriv(xmlIni) == "" } {
	
	set xmlArray [::xmlutils::openFile "." "$file"]
		
	set KPriv(xmlIni) [lindex $xmlArray 0]
	set KPriv(xmlDocIni) [lindex $xmlArray 2]
    }
    
    return $nodeDefSettings
}

#
# Devuelve el valor que tiene la variable "var" en el .ini (no en el default)
#
proc ::kps::getConfigValue { var } {

	global KPriv
	
	set xpath "/KratosConfig/GlobalProjectSettings/Item\[@id='$var'\]"
	
	set dv [::xmlutils::getAttribute $KPriv(xmlIni) $xpath dv]
	
	return $dv
}

#Elimina los nodos obsoletos
proc ::kps::deleteOldNodes { nodeDefSettings globalIniNode } {
	
	set childrenDef [$nodeDefSettings childNodes]
	
	set defIds {}
	
	foreach node $childrenDef {
		lappend defIds [$node getAttribute id ""]
	}
	
	set childrenIni [$globalIniNode childNodes]
	foreach node $childrenIni {
		
		set id [$node getAttribute id ""]
		if { !($id in $defIds) } {
			$node delet
		}
	}
}

#
# Si falta algun par�metro se a�ade ( actualiza el .ini compar�ndolo con el default)
# Si llega una w (est� la ventana abierta) construimos etiquetas y widgets)
#
proc ::kps::updateKratosIniFile { {w ""} } {
    
    global KPriv
    
    # Initialize the configuration file
    set nodeDefSettings [::kps::initializeConfigFile]
    
    set xmlIni $KPriv(xmlIni)
    
    set path "/KratosConfig/GlobalProjectSettings"
    set globalIniNode [$xmlIni selectNodes $path]
    
    ::kps::deleteOldNodes $nodeDefSettings $globalIniNode
    
    set i 0
    foreach node [$nodeDefSettings childNodes] {
	
	set id [$node getAttribute id ""]
	set xpath "${path}/Item\[@id='$id'\]"
	
	set exists [$xmlIni selectNodes $xpath]
	if {$exists == ""} {
	    
	    $globalIniNode appendChild $node
	}
	
	set dv [::xmlutils::getAttribute $xmlIni $xpath dv]
	
	set pid [= [::xmlutils::getAttribute $xmlIni $xpath pid]]
	
	# Get the widget type
	set widgettype [$node getAttribute widget ""]
	# Get the widget state
	set widgetstate [$node getAttribute state ""]
	
	if { $w != "" } {
	    
	    set fcmb ${w}.cmb$id
	    
	    if {$widgettype == "check" } {
		
		grid [ttk::checkbutton $fcmb -variable kps::cmb$id -text $pid] \
		    -row $i -column 0 -padx 3 -pady 10 -sticky nw -in $w
		
		set kps::cmb$id $dv
		
		# Update the state
		if {($widgetstate eq "hidden")||($widgetstate eq "hiddenAll")} {
			grid remove $fcmb	
		} elseif {($widgetstate eq "disabled")} {
			$fcmb configure -state disable
		}
	    } else {
		
		grid [ttk::label ${w}.lbl$id -text "$pid" -padding {10 10 10 10}] \
		    -row $i -column 0 -sticky nw -pady 2 -in $w
		
		# Para sacar los valores del combo o confirmar que es un text si no los hay
		set comboList [split [::xmlutils::getAttribute $xmlIni $xpath values] ","]
		
		if { [llength $comboList] } {
		    set state "readonly"
		} else {
		    set state "normal"
		}
		grid [ttk::combobox $fcmb -values $comboList -state $state -width 20 -textvariable "::kps::cmb$id"] \
		    -row $i -column 1 -padx 3 -pady 10 -sticky nw -in $w
		
		if {$state == "readonly" } {
		    ::xmlutils::setComboValue $xmlIni $xpath $fcmb $dv
		} else {
		    set kps::cmb$id $dv
		}
		
		# Update the state
		if {($widgetstate eq "hidden")||($widgetstate eq "hiddenAll")} {
			grid remove ${w}.lbl$id
			grid remove $fcmb
		} elseif {($widgetstate eq "disabled")} {
			$fcmb configure -state disable
			${w}.lbl$id configure -state disable
		}
	    }				 
	}
	
	incr i
	
    }
    
    # Lo guardamos por si se ha a�adido alg�n nodo
    set iniDir [file dirname [::kfiles::giveConfigFile]]
    ::xmlutils::writeFile "kratos.ini" $iniDir "utf-8" $xmlIni 0 0
}


proc ::kps::WindowbClose {{w .gid.settingWin}} {
	
	global KPriv 

	set xml $KPriv(xmlIni)
	
	#Recorremos cada nodo para actualizar los cambios en el xml
	set path "/KratosConfig/GlobalProjectSettings"
	set globalIniNode [$xml selectNodes $path]
	
	foreach node [$globalIniNode childNodes] {
		
		set id [$node getAttribute id ""]
		set xpath "${path}/Item\[@id='$id'\]"
		
		set selCombo [set ::kps::cmb$id]
		set comboState [::xmlutils::getComboState $xml $xpath]
		
		if { $comboState == "normal" || [$node getAttribute widget ""] == "check" } {
			
			set selCombo [set ::kps::cmb$id]

		} else {
			
			set f "$w.psettings.cmb${id}"
			set selCombo [::xmlutils::getComboValue $xml $xpath $f]
		}
		
		::xmlutils::getAttribute $xml $xpath dv $selCombo
	}
		
	# Al cerrar tenemos que actualizar los cambios en el .ini
	set dir [file dirname [::kfiles::giveConfigFile]]
	::xmlutils::writeFile "kratos.ini" $dir "utf-8" $KPriv(xmlIni) 0 0
	
	
	# Destroy the window widget	
	if {[winfo exists $w]} { 
	    destroy $w
	}
}

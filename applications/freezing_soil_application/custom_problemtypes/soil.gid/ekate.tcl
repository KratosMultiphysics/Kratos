# kratos0.1.tcl  -*- TCL -*-
# Kratos Team - 2003
# http://www.cimne.upc.es/kratos/           http://www.kratos.cimne.upc.es
#---------------------------------------------------------------------------
# This file is written in TCL lenguage 
# For more information about TCL look at: http://www.sunlabs.com/research/tcl/
#
# At least two procs must be in this file:
#
#    InitGIDProject dir - Will be called whenever a project is begun to be used.
#           where dir is the project's directory
#    EndGIDProject - Will be called whenever a project ends to be used.
#
# For more information about GID internals, check the program scripts.
#---------------------------------------------------------------------------


proc MyBitmaps { dir { type "DEFAULT INSIDELEFT"} } {
    global MyBitmapsNames MyBitmapsCommands MyBitmapsHelp kratosPriv



    set MyBitmapsNames(0) "images/Thermal.gif --- \
            images/Displacement.gif \
            images/Pressure.gif \
            images/Source.gif \
            images/BoundCond.gif \
            images/Materials.gif \
            images/Data.gif --- \
            images/Mesh.gif \
            images/Compute.gif --- \
            images/k.gif"
    set MyBitmapsCommands(0) [list \
            [ list -np- HelpWindow "CUSTOM_HELP_REFERENCE" \
		   [file join $dir html Contents.html] \
		   [file join $dir html index.html]] \
            [ list ""] \
	    [list -np- GidOpenConditions Fixed_Displacement] \
	    [list -np- GidOpenConditions Face_Pressure] \
	    [list -np- GidOpenMaterials Solids] \
	    [list -np- GidOpenProblemData] \
	    [ list ""] \
	    "Meshing generate" \
	    "Utilities Calculate" \
	    [ list ""] \
	    [list -np- WebPageKratos $kratosPriv(Web)]
    ]
    set MyBitmapsHelp(0) { "About Thermal . . ." "" \
	"Fixed Displacements" \
        "Face Pressure" \
	"Heat Sources" \
        "Velocity" \
        "Materials" \
	"Problem Data" "" \
        "Generate Mesh" \
        "Calculate" ""\
        "Kratos WebPage"}
    
    # prefix values:
    #          Pre        Only active in the preprocessor
    #          Post       Only active in the postprocessor
    #          PrePost    Active Always

    set prefix Pre

    set kratosPriv(toolbarwin) [CreateOtherBitmaps MyBar "My toolbar" \
	    MyBitmapsNames MyBitmapsCommands \
	    MyBitmapsHelp $dir "MyBitmaps [list $dir]" $type $prefix]
    AddNewToolbar "kratos bar" ${prefix}MyBarWindowGeom \
	    "MyBitmaps [list $dir]"
}

proc EndMyBitmaps {} {
    global kratosPriv
    
    ReleaseToolbar "kratos bar"
    rename MyBitmaps ""
    
    catch { destroy $kratosPriv(toolbarwin) }
}

proc InitGIDProject { dir } {
    global MenuNames MenuEntries MenuCommands MenuAcceler
    global MenuNamesP MenuEntriesP MenuCommandsP MenuAccelerP
    global kratosPriv
    global GidUtils

    set kratosPriv(ProgramName) "ekate" 
    set kratosPriv(VersionNumber) 0.1 
    set kratosPriv(Web) "http://www.cimne.upc.es/kratos/"

    #set GidVersion [.gid.central.s info GiDVersion]
    #set GidVersion [string trim $GidVersion]
    #set GidVersion [string index $GidVersion 0]

    #if { $GidVersion < 7 } {
	#WarnWin [_ "%s v%s is not compatible with GiD version lower than 7" \
	#	$kratosPriv(ProgramName) $kratosPriv(VersionNumber)]
	return
    }
    
    #    WarnWinText $dir
    
    Splash $dir 2501
    
    set ipos [lsearch $MenuNames Help]
    if { $ipos != -1 } {
	set MenuEntries($ipos) [linsert $MenuEntries($ipos) 0 \
                [_ "%s v%s Help" $kratosPriv(ProgramName) $kratosPriv(VersionNumber)] --- ]
	set MenuCommands($ipos) [linsert $MenuCommands($ipos) 0 \
		[list -np- HelpWindow "CUSTOM_HELP_REFERENCE" \
	        [file join $dir html Contents.html] \
		[file join $dir html index.html]] ""]
	set MenuAcceler($ipos) [linsert $MenuAcceler($ipos) 0 ""]
	
	lappend MenuEntries($ipos) \
	        [_ "About %s v%s" $kratosPriv(ProgramName) $kratosPriv(VersionNumber)] --- \
	        [_ "%s WebPage" $kratosPriv(ProgramName)]
	lappend MenuCommands($ipos) [list -np- Splash $dir] "" [list -np- WebPageKratos $kratosPriv(Web)]
	lappend MenuAcceler($ipos) ""
	
    	UpdateMenus
    }
	
    set ipos [lsearch $MenuNamesP Help]
    if { $ipos != -1 } {
	set MenuEntriesP($ipos) [linsert $MenuEntriesP($ipos) 0 \
                [_ "%s v%s Help" $kratosPriv(ProgramName) $kratosPriv(VersionNumber)] --- ]
	set MenuCommandsP($ipos) [linsert $MenuCommandsP($ipos) 0 \
		[list -np- HelpWindow "CUSTOM_HELP_REFERENCE" \
	        [file join $dir html Contents.html] \
		[file join $dir html index.html]] ""]
	set MenuAccelerP($ipos) [linsert $MenuAccelerP($ipos) 0 ""]
	
	lappend MenuEntriesP($ipos) \
	        [_ "About %s v%s" $kratosPriv(ProgramName) $kratosPriv(VersionNumber)]
	lappend MenuCommandsP($ipos) [list -np- Splash $dir]
	lappend MenuAccelerP($ipos) ""
	
    	UpdateMenus
    }
    
    MyBitmaps $dir

    GidChangeDataLabel "Data units" ""
    GidChangeDataLabel "Interval" ""
    GidChangeDataLabel "Conditions" ""
    GidChangeDataLabel "Local Axes" ""
    GidChangeDataLabel "Materials" ""
    
    GidAddUserDataOptions "---" "" 1    
    GidAddUserDataOptions "Fixed Displacement" "GidOpenConditions Fixed_Displacement" 2
    GidAddUserDataOptions "Forces" "GidOpenConditions Forces" 3
    GidAddUserDataOptions "---" "" 4
	GidAddUserDataOptions "Fixed Pressures" "GidOpenConditions Fixed_Pressures" 5
	GidAddUserDataOptions "Surface Pressure" "GidOpenConditions Face_Pressure" 6
    GidAddUserDataOptions "---" "" 7
	GidAddUserDataOptions "Initial Conditions" "GidOpenConditions Initial_Conditions" 8
    GidAddUserDataOptions "---" "" 9
    GidAddUserDataOptions "Model Boundaries" "GidOpenConditions Model_Boundaries" 10
    GidAddUserDataOptions "Contact" "GidOpenConditions Contact" 11
    GidAddUserDataOptions "Tying" "GidOpenConditions Tying" 12
	GidAddUserDataOptions "Hydraulic Jacks" "GidOpenConditions Hydraulic_Jacks" 13
    GidAddUserDataOptions "Activation" "GidOpenConditions Activation" 14
    GidAddUserDataOptions "Groups" "GidOpenConditions Groups" 15
    GidAddUserDataOptions "---" "" 16
    GidAddUserDataOptions "Solids" "GidOpenMaterials Solids" 17
	GidAddUserDataOptions "---" "" 18
    GidAddUserDataOptions "Select Boundaries" GetBoundary 19
	GidAddUserDataOptions "Cut Model" CutModel 20
    
}

proc EndGIDProject {} {
    EndMyBitmaps
}

proc WebPageKratos { dir } {
    global kratosPriv

    VisitWeb $dir
}


proc HelpOnkratos { dir } {
    global kratosPriv

    WarnWin [_ "To obtain help for %s v%s, check the lates news in %s" \
			$kratosPriv(ProgramName) $kratosPriv(VersionNumber) $kratosPriv(Web)] 
}

#  Modified Routine for Kratos 2003
#  February 3rd, 2003
#  Program Logo and Version

proc Splash { dir { TimeOut 0 } } {
    global GIDDEFAULT
   
    if { [.gid.central.s disable windows] } { return }

    if { [ winfo exist .splash]} {
	destroy .splash
	update
    }

    toplevel .splash
        
    set im [image create photo -file [ file join $dir images/splash.gif]]
    set x [expr [winfo screenwidth .splash]/2-[image width $im]/2]
    set y [expr [winfo screenheight .splash]/2-[image height $im]/2]

    wm geom .splash +$x+$y
    wm transient .splash .gid
    wm overrideredirect .splash 1
    pack [label .splash.l -image $im -relief ridge -bd 2]
    
    bind .splash <1> "destroy .splash"
    bind .splash <KeyPress> "destroy .splash"
    raise .splash .gid
    grab .splash
    focus .splash
    update

    if { $TimeOut > 0 } {
        after $TimeOut "if { [ winfo exist .splash] } { destroy .splash  }"
    }
}

proc CutModel {} {
	#domain file path
	set domain_file [GiD_Info Project ModelName]
	set domain_file [split $domain_file "/"]
	set domain_file [lrange $domain_file 0 [expr [llength $domain_file] -2]]
	lappend domain_file "domain.sat"
	set domain_file [join $domain_file "/"]
	#WarnWin [_ $domain_file]
	#get list of all existing volumes
	set i [GiD_Info Geometry NumVolumes]
	set stop 1
	set volume_list {}
	while {$i >= $stop} {
		#determine id of new volume
		set last_volume [GiD_Info Geometry NumVolumes]
		#get layer of current volume
		set str [GiD_Info list_entities volumes $i]
		set end [string first "NumSurfaces" $str]
		set begin [string first "LAYER" $str]
		set layername [string range $str [expr $begin+7] [expr $end-2] ]
		#change to active layer
		GiD_Process Mescape View Layers ToUse $layername
		GiD_Process Mescape
		#read domain volume
		GiD_Process Mescape Files ACISRead $domain_file
		#intersect volume and domain volume
		GiD_Process Mescape Geometry Create IntSolid3D Intersect $i [expr $last_volume+1]
		#WarnWin [_ $volume_id ]
		incr i -1
	}
	GiD_Process Mescape
	#clean up geometry
	GiD_Process Mescape Geometry Delete surface InvertSelection Mescape
	GiD_Process Mescape Geometry Delete line InvertSelection Mescape
	GiD_Process Mescape Geometry Delete point InvertSelection Mescape
}

proc GetBoundary {} {
	set i 1
	set stop [GiD_Info Geometry NumSurfaces]
	#bounding box
	set bbcoords [GiD_Info layers -bbox]
	#WarnWin [_ $bbcoords ]
	set bbcoords [lindex $bbcoords 0]
	set bbcoords [split $bbcoords ]
	set bottom [expr [expr [lindex $bbcoords 2] < [lindex $bbcoords 5]] ? [lindex $bbcoords 2] : [lindex $bbcoords 5]]
	while {$i <= $stop} {
		set str [GiD_Info list_entities surfaces $i]
		set end [string first "conditions" $str]
		set begin [string first "HigherEntity" $str]
		set numVol [string range $str [expr $begin+14] [expr $end-2] ]
		#if surface is boundary
		#WarnWin [_ $str ]
		if { [string compare $numVol "1"] == 0 } {
			#GiD_Process Mescape Data Conditions AssignCond IsBoundary Change 1 $i Mescape
			#check for position of surface
			set aend [string first "Normal" $str]
			set abegin [string first "Center" $str]
			set center [string range $str [expr $abegin+8] [expr $aend-2] ]
			set center [split $center]
			#for gid_8.1.1b
			#set bend [string first "END" $str]
			#for gid_8.1.7
			set bbegin $aend
			#set normal [string range $str [expr $bbegin+8] [expr $bend-2] ]
			set normal [string range $str [expr $bbegin+8] [string length $str] ]
			set normal [split $normal]
			set znormal [expr abs([lindex $normal 2]) - 0.00001]
			if { [expr $bottom == [lindex $center 2] ] } then {
				#WarnWin [_ "Surface is bottom surface" ]
				GiD_Process Mescape Data Conditions AssignCond IsBottom Change 1 $i Mescape
			} elseif { $znormal < 0.0 } then {
					#WarnWin [_ "Surface is side surface"]
					GiD_Process Mescape Data Conditions AssignCond IsSide Change 1 $i Mescape
					#GiD_Process Mescape Data Conditions AssignCond Surface_Water_Pressure Change 1 0.0 $i Mescape
			} else {
				#WarnWin [_ "neither Bottom nor side: must be top"]
				GiD_Process Mescape Data Conditions AssignCond IsTop Change 1 $i Mescape
			}
		}
		incr i
	}
}

proc GetLayerNodes { layername } {
        set surfaces [GiD_Info layers -entities surfaces $layername ]
        #WarnWin [_ $surfaces]
	#get all lines that are attached to the surfaces belonging to the specified layer
        set lines {}
	foreach surf $surfaces {
	        set lineinfo [GiD_Info list_entities surfaces $surf]
		WarnWin [_ $lineinfo]
                set lineinfo [split $lineinfo "\n"]
                foreach line $lineinfo {
			set line [split $line]
			if {[lindex $line 0] == "Line:"} {
				lappend lines [lindex $line 1]
			}
		}
                
	}
        set lines [lsort -unique $lines]
	#WarnWin [_ $lines]
	#get all points belonging to the specified lines
	set points {}
	foreach line $lines {
		set lineinfo [GiD_Info list_entities lines $line]
		set lineinfo [split $lineinfo "\n"]
		foreach item $lineinfo {
			set item [split $item]
			if {[lindex $item 0] == "Points:"} {
				lappend points [lindex $item 1]
				lappend points [lindex $item 2]
			}
		}
	}
	set points [lsort -unique $points]
	#WarnWin [_ $points]
	set nodes [GiD_Info list_entities PreNodes 1]
	#WarnWin [_ $nodes]
	return $points
}

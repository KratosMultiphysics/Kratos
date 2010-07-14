###############################################################################
#
#    NAME: utils.tcl
#
#    PURPOSE: Kratos problem type utilities procedures
#
#    QUANTECH ATZ-DEVELOPMENT DEPARTMENT
#
#    AUTHOR : G. Socorro
#
#    CREATED AT: 01/11/09
#
#    LAST MODIFICATION : add the procedure GetDefinedMeshGiDEntities to get the GiD mesh entities that belong to a group identifier
#
#    VERSION : 0.2
#
#    HISTORY:
#
#     0.2- 22/04/10-G. Socorro, add the procedure GetDefinedMeshGiDEntities to get the GiD mesh entities that belong to a group identifier
#     0.1- 01/11/09-G. Socorro, create a base source code
#
###############################################################################

#  Create a base namespace KUtils

namespace eval ::KUtils:: {
 
}

proc ::KUtils::GetDefinedMeshGiDEntities {groupid {etype point} {what Elements} {isquadratic 0}} {

    # WarnWinText "groupid:$groupid etype:$etype what:$what isquadratic:$isquadratic"
    set gidetype [list Linear Triangle Quadrilateral Tetrahedra Hexahedra]

    ::GidUtils::DisableGraphics
    set frozen_layers [::KEGroups::UnFreezeLayers]
    set infoProj [GiD_Info Project ViewMode]
    # WarnWinText "frozen_layers:$frozen_layers"
    set meshelist [list]
    if {$infoProj == "GEOMETRYUSE"} {
	GiD_Process MEscape Meshing MeshView
    }
    set condid "${etype}_groups"
    # WarnWinText "condid:$condid"
    if {$what == "Elements"} {
	# Get all defined elements for this group
	foreach cprop [GiD_Info conditions $condid mesh] {
	    lassign $cprop cid eid - cgroupid
	    # WarnWinText "cid:$cid eid:$eid cgroupid:$cgroupid"
	    if {($cgroupid == $groupid) && ($cid =="E")} {
		lappend meshelist $eid
	    }
	}
	::KEGroups::FreezeLayers $frozen_layers
	::GidUtils::EnableGraphics
	
	return $meshelist
    } elseif {$what == "Nodes"} {
	set nlist [list]
	# Get all defined nodes for this group
	if {$etype !="point"} { 
	    foreach cprop [GiD_Info conditions $condid mesh] {
		lassign $cprop cid eid - cgroupid
		# WarnWinText "cid:$cid eid:$eid cgroupid:$cgroupid"
		if {($cgroupid == $groupid) && ($cid !="E") && ($cid !="N")} {
		    # Get the element properties
		    foreach GiDElemType $gidetype {
			# WarnWinText "GiDElemType:$GiDElemType"
			switch $GiDElemType {
			    "Triangle" {
				foreach nodeid [::KUtils::GetTriangleEdgeNodes $isquadratic $cid $eid] {
				    if {$nodeid ni $nlist} {
					lappend nlist $nodeid
				    }  
				}
			    }
			    "Quadrilateral" {
				foreach nodeid [::KUtils::GetQuadrilateralEdgeNodes $isquadratic $cid $eid] {
				    if {$nodeid ni $nlist} {
					lappend nlist $nodeid
				    }  
				}
			    }
			    "Tetrahedra" {
				# WarnWinText "cid:$cid eid:$eid cgroupid:$cgroupid nodeidlist: [::KUtils::GetTetrahedraFaceNodes $isquadratic $cid $eid]"
				foreach nodeid [::KUtils::GetTetrahedraFaceNodes $isquadratic $cid $eid] {
				    if {$nodeid ni $nlist} {
					lappend nlist $nodeid
				    }  
				}
			    }
			    "Hexahedra" {
				foreach nodeid [::KUtils::GetHexahedraFaceNodes $isquadratic $cid $eid] {
				    if {$nodeid ni $nlist} {
					lappend nlist $nodeid
				    }  
				}
			    }
			}
		    }
		}
	    }
	} else {
	    foreach cprop [GiD_Info conditions $condid mesh] {
		lassign $cprop cid eid - cgroupid
		if {($cgroupid == $groupid) && ($cid =="N")} {
		    lappend nlist $eid
		}
	    }
	}

	::KEGroups::FreezeLayers $frozen_layers
	::GidUtils::EnableGraphics
	
	return $nlist
    }
}

proc ::KUtils::GetTriangleEdgeNodes {isquadratic elemid edgeid} {
    # WarnWinText "elemid:$elemid edgeid:$edgeid"
    set nprop [lrange [GiD_Info Mesh Elements Triangle $elemid] 1 end-1]
    set nodelist [list]
    if {[llength $nprop]>0} {
	switch $isquadratic {
	    "0" {
		switch $edgeid {
		    "1" {
			foreach lnode [list 0 1] {
			    lappend nodelist [lindex $nprop $lnode]
			}
		    }
		    "2" {
			foreach lnode [list 1 2] {
			    lappend nodelist [lindex $nprop $lnode]
			}
		    }
		    "3" {
			foreach lnode [list 2 0] {
			    lappend nodelist [lindex $nprop $lnode]
			}
		    }
		}
	    }
	    "1" {
		switch $edgeid {
		    "1" {
			foreach lnode [list 0 3 1] {
			    lappend nodelist [lindex $nprop $lnode]
			}
		    }
		    "2" {
			foreach lnode [list 1 4 2] {
			    lappend nodelist [lindex $nprop $lnode]
			}
		    }
		    "3" {
			foreach lnode [list 2 5 0] {
			    lappend nodelist [lindex $nprop $lnode]
			}
		    }
		}
	    }
	}
    }
    return $nodelist
}

proc ::KUtils::GetQuadrilateralEdgeNodes {isquadratic elemid edgeid} {

    set nprop [lrange [GiD_Info Mesh Elements Quadrilateral $elemid] 1 end-1]
    set nodelist [list]
    if {[llength $nprop]>0} {
	switch $isquadratic {
	    "0" {
		switch $edgeid {
		    "1" {
			foreach lnode [list 0 1] {
			    lappend nodelist [lindex $nprop $lnode]
			}
		    }
		    "2" {
			foreach lnode [list 1 2] {
			    lappend nodelist [lindex $nprop $lnode]
			}
		    }
		    "3" {
			foreach lnode [list 2 3] {
			    lappend nodelist [lindex $nprop $lnode]
			}
		    }
		    "4" {
			foreach lnode [list 3 0] {
			    lappend nodelist [lindex $nprop $lnode]
			}
		    }
		}
	    }
	    "1" - "2" {
		switch $edgeid {
		    "1" {
			foreach lnode [list 0 4 1] {
			    lappend nodelist [lindex $nprop $lnode]
			}
		    }
		    "2" {
			foreach lnode [list 1 5 2] {
			    lappend nodelist [lindex $nprop $lnode]
			}
		    }
		    "3" {
			foreach lnode [list 2 6 3] {
			    lappend nodelist [lindex $nprop $lnode]
			}
		    }
		    "4" {
			foreach lnode [list 3 7 0] {
			    lappend nodelist [lindex $nprop $lnode]
			}
		    }
		}
	    }
	}
    }
    return $nodelist
}

proc ::KUtils::GetTetrahedraFaceNodes {isquadratic elemid edgeid} {
    # WarnWinText "elemid:$elemid edgeid:$edgeid"
    set nprop [lrange [GiD_Info Mesh Elements Tetrahedra $elemid] 1 end-1]
    set nodelist [list]
    if {[llength $nprop]>0} {
	switch $isquadratic {
	    "0" {
		switch $edgeid {
		    "1" {
			foreach lnode [list 0 1 2] {
			    lappend nodelist [lindex $nprop $lnode]
			}
		    }
		    "2" {
			foreach lnode [list 1 3 2] {
			    lappend nodelist [lindex $nprop $lnode]
			}
		    }
		    "3" {
			foreach lnode [list 2 3 0] {
			    lappend nodelist [lindex $nprop $lnode]
			}
		    }
		    "4" {
			foreach lnode [list 3 1 0] {
			    lappend nodelist [lindex $nprop $lnode]
			}
		    }
		}
	    }
	    "1" {
		switch $edgeid {
		  "1" {
		      foreach lnode [list 0 1 2 4 5 6] {
			  lappend nodelist [lindex $nprop $lnode]
		      }
		    }
		    "2" {
			foreach lnode [list 1 3 2 8 9 5] {
			    lappend nodelist [lindex $nprop $lnode]
			}
		    }
		    "3" {
			foreach lnode [list 2 3 0 9 7 6] {
			    lappend nodelist [lindex $nprop $lnode]
			}
		    }
		    "4" {
			foreach lnode [list 3 1 0 8 4 7] {
			    lappend nodelist [lindex $nprop $lnode]
			}
		    }
		}
	    }
	}
    }
    return $nodelist
}

proc ::KUtils::GetHexahedraFaceNodes {isquadratic elemid edgeid} {

    set nprop [lrange [GiD_Info Mesh Elements Hexahedra $elemid] 1 end-1]
    set nodelist [list]
    if {[llength $nprop]>0} {
	switch $isquadratic {
	    "0" {
		switch $edgeid {
		    "1" {
			foreach lnode [list 0 1 2 3] {
			    lappend nodelist [lindex $nprop $lnode]
			}
		    }
		    "2" {
			foreach lnode [list 0 3 7 4] {
			    lappend nodelist [lindex $nprop $lnode]
			}
		    }
		    "3" {
			foreach lnode [list 0 4 5 1] {
			    lappend nodelist [lindex $nprop $lnode]
			}
		    }
		    "4" {
			foreach lnode [list 1 5 6 2] {
			    lappend nodelist [lindex $nprop $lnode]
			}
		    }
		    "5" {
			foreach lnode [list 2 6 7 3] {
			    lappend nodelist [lindex $nprop $lnode]
			}
		    }
		    "6" {
			foreach lnode [list 4 7 6 5] {
			    lappend nodelist [lindex $nprop $lnode]
			}
		    }
		}
	    }
	    "1" {
		switch $edgeid {
		    "1" {
			foreach lnode [list 0 1 2 3 8 9 10 11] {
			    lappend nodelist [lindex $nprop $lnode]
			}
		    }
		    "2" {
			foreach lnode [list 0 3 7 4 11 15 19 12] {
			    lappend nodelist [lindex $nprop $lnode]
			}
		    }
		    "3" {
			foreach lnode [list 0 4 5 1 12 16 13 8] {
			    lappend nodelist [lindex $nprop $lnode]
			}
		    }
		    "4" {
			foreach lnode [list 1 5 6 2 13 17 14 9] {
			    lappend nodelist [lindex $nprop $lnode]
			}
		    }
		    "5" {
			foreach lnode [list 2 6 7 3 14 18 15 10] {
			    lappend nodelist [lindex $nprop $lnode]
			}
		    }
		    "6" {
			foreach lnode [list 4 7 6 5 19 18 17 16] {
			    lappend nodelist [lindex $nprop $lnode]
			}
		    }
		}
	    }
	    "2" {
		switch $edgeid {
		    "1" {
			foreach lnode [list 0 1 2 3 8 9 10 11 20] {
			    lappend nodelist [lindex $nprop $lnode]
			}
		    }
		    "2" {
			foreach lnode [list 0 3 7 4 11 15 19 11 24] {
			    lappend nodelist [lindex $nprop $lnode]
			}
		    }
		    "3" {
			foreach lnode [list 0 4 5 1 12 16 15 8 21] {
			    lappend nodelist [lindex $nprop $lnode]
			}
		    }
		    "4" {
			foreach lnode [list 1 5 6 2 13 17 14 9 22] {
			    lappend nodelist [lindex $nprop $lnode]
			}
		    }
		    "5" {
			foreach lnode [list 2 6 7 3 14 18 15 10 23] {
			    lappend nodelist [lindex $nprop $lnode]
			}
		    }
		    "6" {
			foreach lnode [list 4 7 6 5 19 18 17 16 25] {
			    lappend nodelist [lindex $nprop $lnode]
			}
		    }
		}
	    }
	}
    }
    return $nodelist
}

proc ::KUtils::AutoMkIndex {} {

    set bpath "D:/GiD/GiD 9.2.7b/problemtypes/kratos.gid"

    # For scripts directory
    set allsdirname [list]
    set scriptspath [file join $dir/scripts/]
    lappend allsdirname "$scriptspath"
    cd $scriptspath
    set dirlist [glob *]
    # WarnWinText "dirlist:$dirlist\n\n"
    foreach cdir $dirlist {
	if {[file isdirectory $cdir]} {
	    lappend allsdirname "$cdir/*.tcl"
	}
    }
    if {[llength $allsdirname]>0} {
	auto_mkindex "$allsdirname" 
    }
    
}  

proc ::KUtils::ProjectInfo {what} {
    # ABSTRACT:
    # This procedure return some project information
    # ProblemType|ModelName|AreChanges|LayerToUse|ViewMode|Quadratic|RenderMode|ExistPost|Debug
    # |TmpDirectory|MustReMesh|LastElementSize|BackgroundFilename

    set result ""
    catch { set result [GiD_Info Project $what]}
    return $result
    
}

proc ::KUtils::ShowAllNamespace {} {
    set allnames [namespace children ::]
    foreach name $allnames {
	WarnWinText "$name\n"
    }
}

proc ::KUtils::CheckProjectSave {} {
    # ABSTRACT:
    # This procedure check if project have been saved
    
    # Error
    set val "0" 
    
    # Obtain project information
    set ProjectName ""
    catch { set ProjectName [::sadUtils::GiveProjectName] }
    if {$ProjectName==""} {
	WarnWin [= "The project must be saved"].
	set val "0"
    } else {
	set val "1"
    }
    return $val
}


proc ::KUtils::DeleteNamespaces {} {
    
    set nplist [list ""]
    foreach np $nplist {
	catch { namespace delete $np }
    }
}



proc ::KUtils::ChangeMatrixEntries {fpath matrixid MatrixValues} {
    
    set count 0; set ccol 0
    foreach elem $MatrixValues {
	incr count 1; incr ccol 1
	set crow 1
	if {($count>3)&&($count<=6)} {
	    set crow 2
	} elseif {$count>6} { 
	    set crow 3
	} 
	if {($count==4)||($count==7)} {
	    set ccol 1
	}
	# WarnWinText "count:$count ccol:$ccol crow:$crow"
	set endpath $fpath.e${matrixid}${crow}${ccol}
	set empty [$endpath delete 0 end]
	set empty [$endpath insert 0 $elem]
    }
}

proc ::KUtils::GeometryMode {} {
    
    set mode [lindex [GiD_Info Project] 4]
    ;# WarnWin "mode:$mode"
    if {$mode == "MESHUSE"} {
	GiD_Process Mescape Geometry
    }
}


proc ::KUtils::GiveTempDir {} {
    # Get temporal directory
    global env
    
    if { [info exists env(TEMP)] } {
	return $env(TEMP)
    }
    if { [info exists env(TMP)] } {
	return $env(TMP)
    }
    if { [info exists env(WINDIR) && [file isdir [file join $env(WINDIR) temp]]] } {
	return [file join $env(WINDIR) temp]
    }
    if { [file isdir /tmp] } {
	return /tmp
    }
    return ""
}

proc ::KUtils::GiveProjectName {} {
    # Give the project name
    
    set ProjectName [GiD_Info Project "ModelName"]
    if { [file extension $ProjectName] == ".gid" } {
	set ProjectName [file root $ProjectName]
    }
    
    if { $ProjectName == "UNNAMED" } {
	return ""
    } else {  
	set basename [file tail $ProjectName]  
	return $basename
    }
}

proc ::KUtils::GiveProjectDir {} {
    # Give the project directory
    
    set ProjectName [GiD_Info Project "ModelName"]
    if { [file extension $ProjectName] == ".gid" } {
	set ProjectName [file root $ProjectName]
    }
    
    if { $ProjectName == "UNNAMED" } {
	return ""
    }
    
    set directory $ProjectName.gid
    if { [file pathtype $directory] == "relative" } {
	set directory [file join [pwd] $directory]
    }
    return $directory
}

proc ::KUtils::GetPaths {what} {
    # ABSTRACT:
    # Utilities procedure to get some path and name 
    # What can be:
    # PDir   : Return the project directory
    # PTDir  : Return the problem type directory
    # TDir   : Return a temporal directory for windows
    # PName  : Return the project name
    # PFPath : Return the project directory including the project name in the path

    set rpath ""
    switch -exact -- $what {
	"PDir" {
	    set ProblemDirectory [::KUtils::GiveProjectDir]
	    set rpath $ProblemDirectory
	}
	"PTDir" {
	    set ptdir [GiD_Info problemtypepath] 
	    set rpath $ptdir
	}
	"TDir" {
	    set tempdir [::KUtils::GiveTempDir]
	    set rpath $tempdir
	}
	"PName" {
	    set ProjectName [::KUtils::GiveProjectName]
	    set rpath $ProjectName
	}
	"PFPath" {
	    set ProjectName [::KUtils::GiveProjectName]
	    set ProblemDirectory [::KUtils::GiveProjectDir]
	    set basename [file tail $ProjectName]
	    set fullname [file native [file join $ProblemDirectory $basename]]
	    set rpath $fullname
	}
    }
    return $rpath
}

proc ::KUtils::ExecCmd {command} {
    GiD_Process MEscape $command  
}

proc ::KUtils::GetDate {} {
    set t [clock seconds]
    return [clock format $t -format "%a %d %b %Y"]
}

proc ::KUtils::GetHour {} {
    set t [clock seconds]
    return [clock format $t -format "%H:%M:%S"]
}

proc ::KUtils::Duration { int_time } {
     set timeList [list]
     foreach div {86400 3600 60 1} mod {0 24 60 60} name {day hr min sec} {
         set n [expr {$int_time / $div}]
         if {$mod > 0} {set n [expr {$n % $mod}]}
         if {$n > 1} {
             lappend timeList "$n ${name}s"
         } elseif {$n == 1} {
             lappend timeList "$n $name"
         }
     }
     return [join $timeList]
 }


proc ::KUtils::parseTreeStr { texto } {
	
	set Ltexto [split $texto ""]
	
	if { ":" in $Ltexto || "/" in $Ltexto || "$" in $Ltexto  || "." in $Ltexto || "\\" in $Ltexto  || "%" in $Ltexto } {
		return -1
	}
	#Para reducir el número de espacios
	set texto [string trim $texto]
	
	#Sustituye espacios por underscore '_'
	set texto [string map { " " "_" } $texto]
	#[::KUtils::underscore $texto]
	
	return $texto
}

proc ::KUtils::underscore { texto {put "yes" } } {
	
	if { $put == "yes" } {
		[string map { " " "_" } $texto]
	} else {
		[string map { "_" " " } $texto]
	}
}

#e.g. from {a {b {c d}} e} returns {a b c d e}
proc ::KUtils::FlatEmbeddedLists { a } {
    set res {}
    foreach i $a {
    set i [string trim $i]
    if { [llength $i] > 1 } {
        eval lappend res [FlatEmbeddedLists $i]
    } else {
        lappend res $i
    }
    }
    return $res
}



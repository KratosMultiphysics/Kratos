###############################################################################
#
#    NAME: utils.tcl
#
#    PURPOSE: Kratos problem type utilities procedures
#
#    QUANTECH ATZ-DEVELOPMENT DEPARTMENT
#
#    AUTHOR : G. Socorro => GSM
#
#    CREATED AT: 01/11/09
#
#    HISTORY:
#
#     1.4- 06/03/13- GSM, add the proc IsProcessRunningVar, ProjectProcessAreRunning and ReadProjectProcessRunningState
#     1.3- 17/06/13- GSM, delete the proc DeleteAllGroupIdentifier
#     1.2- 22/10/12- J.Garate, Support for new GiD Groups
#     1.1- 26/09/12- J.Garate, update Parsing function to allow spacing when renaming
#     1.0- 20/07/12-GSM, update the proc ElemenDosListas and rename it to TwoListRepeatedItems
#     0.9- 13/07/12- J. Garate, ElemenDosListas Compara 2 listas a ver si hay coincidencias
#     0.8- 07/05/12- J. Garate, update renaming groups restrictions
#     0.7- 03/05/12-GS, add the procedure DeleteAllGroupIdentifier
#     0.6- 22/06/11-GS, create the proc CreateTBEFiles and DeleteTBEFiles to work with the TBE files
#     0.5- 09/12/10-GS, modify the procedure GetDefinedMeshGiDEntities to get all nodes that defined a surface when using shell elements
#     0.4- 14/09/10-GS, modify ReadResultsFromFiles procedure to load lst GiD files ${appid}.post.lst
#     0.3- 08/09/10-GS, add the procedure ReadResultsFromFiles to check and read kratos result files in the GiD postprocess
#     0.2- 22/04/10-GS, add the procedure GetDefinedMeshGiDEntities to get the GiD mesh entities that belong to a group identifier
#     0.1- 01/11/09-GS, create a base source code
#
###############################################################################

#  Create a base namespace KUtils

namespace eval ::KUtils:: {
 variable IsProcessRunning 0
}


proc ::KUtils::CreateTBEFiles {} {
    
    # Get the problem type subdirectory
    set dir [::KUtils::GetPaths "PTDir"]
    
    # For scripts directory
    set scriptspath [file join $dir/scripts/]
    # WarnWinText "scriptspath:$scriptspath"
    cd $scriptspath
    set dirlist [glob *]
    # WarnWinText "dirlist:$dirlist\n\n"
    foreach cdir $dirlist {
        # WarnWinText "cdir:$cdir"
        if {[file isdirectory $cdir]} {
            cd $cdir
            set tcllist ""
            catch { set tcllist [glob *]}
            # WarnWinText "tcllist:$tcllist\n"
            if {[llength $tcllist]} {
                foreach tcl_cl1 $tcllist {
                    # WarnWinText "tcl_cl1:$tcl_cl1"
                    if {[file isdirectory $tcl_cl1]} {
                        cd $tcl_cl1
                        set tcllist_cl2 ""
                        catch { set tcllist_cl2 [glob *]}
                        # WarnWinText "tcllist_cl2:$tcllist_cl2"
                        if {[llength $tcllist_cl2]} {
                            foreach tcl_cl2 $tcllist_cl2 {
                                if {[file isdirectory $tcl_cl2]} {
                                    cd $tcl_cl2
                                    set tcllist_cl3 ""
                                    catch { set tcllist_cl3 [glob *] }
                                    # WarnWinText "tcllist_cl3:$tcllist_cl3\n"
                                    if {[llength $tcllist_cl3]} {
                                        foreach tcl_cl3 $tcllist_cl3 {
                                            if {[file extension $tcl_cl3]==".tcl"} {
                                                loadtbefile -create $tcl_cl3
                                            }
                                        }
                                    }
                                    cd ..
                                } else {
                                    if {[file extension $tcl_cl2]==".tcl"} {
                                        loadtbefile -create $tcl_cl2
                                    }
                                }
                            }
                        }
                        cd ..
                    } else {
                        if {[file extension $tcl_cl1]==".tcl"} {
                            loadtbefile -create $tcl_cl1
                        }            
                    }
                    }
            }
            cd ..
        } else {
            # WarnWin "cdir:$cdir"
            if {[file extension $cdir]==".tcl"} {
                loadtbefile -create $cdir
            }
        }
    }
}



proc ::KUtils::DeleteTBEFiles {} {
    
    # Get the problem type subdirectory
    set dir [::KUtils::GetPaths "PTDir"]
    
    # For scripts directory
    set scriptspath [file join $dir/scripts/]
    cd $scriptspath
    set dirlist [glob *]
    # WarnWinText "dirlist:$dirlist\n\n"
    foreach cdir $dirlist {
        if {[file isdirectory $cdir]} {
            cd $cdir
            set tbelist ""
            catch { set tbelist [glob *] }
            # WarnWinText "cdir:$cdir => tbelist:$tbelist\n"
            if {[llength $tbelist]} {
                foreach tbe_level1 $tbelist {
                    # WarnWinText "tbe_level1:$tbe_level1"
                    if {[file isdirectory $tbe_level1]} {
                        cd $tbe_level1
                        set tbe_level2list ""
                        catch { set tbe_level2list [glob *]}
                        # WarnWinText "tbe_level2list:$tbe_level2list"
                        if {[llength $tbe_level2list]} {
                            foreach tbe_level2 $tbe_level2list {
                                # WarnWinText "tbe_level2:$tbe_level2"
                                if {[file extension $tbe_level2]==".tbe"} {
                                    catch { file delete $tbe_level2}
                                }
                            }
                        }
                        cd ..
                    } else {
                        if {[file extension $tbe_level1]==".tbe"} {
                            catch { file delete $tbe_level1}
                        }
                    }
                }
            }
            cd ..
        } else {
            if {[file extension $cdir]==".tbe"} {
                # WarnWinText "cdir:$cdir =>[file extension $cdir]\n"
                catch { file delete $cdir}
            }
        }
    }
}

proc ::KUtils::ReadResultsFromFiles {appid rtype pmode {what CheckRFiles}} {
    # ABSTRACT
    # Check and read the Kratos result files to the GiD postprocess
    # Arguments
    # appid  => Application identifier
    # rtype  => Result type ["Single"|"Multiples"]
    # pmode  => Postprocess mode ["Ascii"|"Binary"]

    # ok
    set ok 1 
    # WarnWinText "appid:$appid rtype:$rtype pmode:$pmode what:$what"
    set fpath [::KUtils::GetPaths "PFPath"]
    # WarnWinText "fpath:$fpath"
    if {$appid =="Fluid"} {
        # For fluid application load the lst file
        set ext ".post.lst"
        # End file name
        set efn "${fpath}${ext}"
        if {$what =="CheckRFiles"} {
            if {[file exists $efn] == 1} {
                return 1
            } else {
                return 0
            }
        } elseif {$what =="ReadRFiles"} {
            if {[file exists $efn] == 1} {
                cd [::KUtils::GetPaths "PDir"]
                GiD_Process MEscape files read [file tail $efn]
            }
        }

    } else {
        if {$rtype =="Single"} {

            # Single file
            if {$pmode =="Ascii"} {
                set ext "${appid}_0.post.res"
            } elseif {$pmode =="Binary"} {
                set ext "${appid}.post.bin"
            }
            
            # End file name
            set efn "${fpath}${ext}"
            if {$what =="CheckRFiles"} {
                if {[file exists $efn] == 1} {
                    return 1
                } else {
                    return 0
                }
            } elseif {$what =="ReadRFiles"} {
                if {[file exists $efn] == 1} {
                    cd [::KUtils::GetPaths "PDir"]
                    GiD_Process MEscape files read [file tail $efn]
                }
            }
            
        } elseif {$rtype =="Multiples"} {
            # Get the project name
            set pname [::KUtils::GetPaths "PName"]
            if {$pmode =="Ascii"} {
                set ext "${appid}.post.res"
                set fp "${pname}${appid}_*.post.res"
            } elseif {$pmode =="Binary"} {
                set ext "${appid}.post.bin"
                set fp "${pname}${appid}_*.post.bin"
            }
            
            if {$what =="CheckRFiles"} {
		
		set efn "${fpath}${ext}"
                if {[file exists $efn] == 1} {
		    return 1
		} else {
		    cd [::KUtils::GetPaths "PDir"]
		    set flist [glob $fp]
		    # WarnWinText "flist:$flist\n\n"
		    if {[llength $flist]} {
			# Check the first file
			set checkfname [lindex $flist 0]
			if {[file exists $checkfname] == 1} {
			    return 1
			} else {
			    return 0
			}
		    } else {
			return 0
		    }
		}
            } elseif {$what =="ReadRFiles"} {
                # Get the file list
                set efn "${fpath}${ext}"
                if {[file exists $efn] == 1} {
                    cd [::KUtils::GetPaths "PDir"]
		    GiD_Process MEscape Files ReadMultiple $ext
 	        } else {
		    cd [::KUtils::GetPaths "PDir"]
		    set flist [glob $fp]
		    if {[llength $flist]} {
			set endflist [list]
			foreach fileid $flist {
			    if {[file exists $fileid] == 1} {
				lappend endflist $fileid
			    }
			}
			if {[llength $endflist]} {
			    # Read all file in the GiD postprocess
			    # Multiple files
			    GiD_Process MEscape Files ReadMultiple $endflist
			}

		    }   
		} 
	    }
	}
    }
    return $ok
}

proc ::KUtils::GetDefinedMeshGiDEntities {groupid {etype point} {what Elements} {isquadratic 0}} {
    WarnWinText "KUtils::GetDefinedMeshGiDEntities MUST NOT BE CALLED, it use old condiciones instead of current groups!!"
        
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
                } elseif {($cgroupid == $groupid) && ($cid =="E")} {
                    # Condition applied over element by we need only the node identifier of this elements
                     # Get the element properties
                    foreach GiDElemType $gidetype {
                        # WarnWinText "GiDElemType:$GiDElemType"
                        switch $GiDElemType {
                            "Triangle" {
                                foreach nodeid [lrange [GiD_Info Mesh Elements Triangle $eid] 1 end-1] {
                                    if {$nodeid ni $nlist} {
                                        lappend nlist $nodeid
                                    }  
                                }
                            }
                            "Quadrilateral" {
                                foreach nodeid [lrange [GiD_Info Mesh Elements Quadrilateral $eid] 1 end-1] {
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


proc ::KUtils::parseTreeStr { texto { allowspacing 0} } {
        
        set Ltexto [split $texto ""]
    
    if { $allowspacing == 0 } {
        if {"\ " in $Ltexto} { return -1 }
        }
        if { ":" in $Ltexto || "$" in $Ltexto || "@" in $Ltexto  || "." in $Ltexto || "\\" in $Ltexto  || "%" in $Ltexto } {
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

proc ::KUtils::TwoListRepeatedItems { lista1 lista2 {ccase "nodes"} } {

    # wa "lista1:$lista1 lista2:$lista2"
    if { 1 } {
        
        set retlist ""
        if {$ccase == "nodes"} {
            foreach item $lista1 {
                if {$item in $lista2} {
                    lappend retlist $item
                }
            }
        } elseif {$ccase == "groups"} {
            foreach item $lista1 {
                set dif [lsearch $item $lista2]
                if {$dif != "-1"} {
                    lappend retlist $item
                }
            }
        }
        
        return $retlist
        
    } else {
        
        
        set retlist ""
        if {$elem == "nodes"} {
            foreach item $lista1 {
                set dif [lsearch $item $lista2]
                # msg "dif $dif"
                if {$dif != "-1"} {
                    lappend retlist $item
                }
            }
        } elseif {$elem == "groups"} {
            foreach item $lista1 {
                set dif [lsearch $item $lista2]
                if {$dif != "-1"} {
                    lappend retlist $item
                }
            }
        }
        return $retlist
    }
}

proc ::KUtils::IsProcessRunningVar {what {cvalue 0}} {
    # ABSTRACT: Get/set the process control variable
    # wa "whta:$what cvalue:$cvalue"
    set rval 0 
    switch -exact -- $what {
	"Get" {
	    if {[info exists ::KUtils::IsProcessRunning]} {
		set rval $::KUtils::IsProcessRunning
	    }
	}
	"Set" {
	    set ::KUtils::IsProcessRunning $cvalue
	}
    }
    # wa "::KUtils:IsProcessRunning:$::KUtils::IsProcessRunning"
    return $rval
}

proc ::KUtils::ProjectProcessAreRunning { } {
    # ABSTRACT:  We check if there are any process running for this project
    
    # Get the current project name
    set currentPName [::KUtils::GetPaths "PName"]
    #wwt "currentPName:$currentPName   exists RunProcInfo:[info exists ::RunProcInfo]"

    # No process running
    set IsRunning 0
    if {[info exists ::RunProcInfo] } {
	foreach processArray $::RunProcInfo {
	    # wwt "processArray:$processArray"
	    set proces [lindex $processArray 0]
	    # wwt "proces: $proces == currentPName:$currentPName"
	    if { $proces == $currentPName} {
		set IsRunning 1
	    }
	}
    } 
    
    return $IsRunning
}

proc ::KUtils::ReadProjectProcessRunningState { } {
    # ABSTRACT: Read the process running state from the file projectname.info
 
    set pdir [::KUtils::GetPaths "PDir"] 
    # wwt "pdir:$pdir"
    if {$pdir eq ""} return ""
 
    # Get the current project name
    set currentPName [::KUtils::GetPaths "PName"]
    
    # startkratos startmainloop endkratos
    set cstate [list 0 0 0]
    set data [list] 

    # Get the file path
    set fpath [file nativename [file join ${pdir} currentPName.info]]
    # wwt "fpath:$fpath"
    if {[file exists ${fpath}]} {
	set cval [catch { set resfile [open "$fpath" "r"] }]
	# All is ok
	if {$cval==0} {
	    set data [read $resfile]
	    # wwt data:$data
	    close $resfile
	} else {
	    return "$cstate"
	}
	# wwt data:$data
    }
    set findstartkratos [string last "Importing" $data]
    if {$findstartkratos !="-1"} {
	lset cstate 0 1
    }
    set findstartmainloop [string last "Real time calculation" $data]
    if {$findstartmainloop !="-1"} {
	lset cstate 1 1
    }
    
    # Ok
    set findendkratosok [string last "ANALYSIS COMPLETED" $data]
    # Error
    set findendkratoserror [string last "KRATOS TERMINATED WITH ERROR" $data]
    if {($findendkratosok !="-1")|| ($findendkratoserror !="-1")} {
	lset cstate 2 1
    } 
	
    wa "findstartkratos:$findstartkratos\nfindstartmainloop:$findstartmainloop\nfindendkratos:$findendkratos"   
  
    return $cstate
}

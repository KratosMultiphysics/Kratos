proc InitGIDProject {dir} {
	# Custom Menu
}

proc BeforeMeshGeneration {elementsize} {

	global surf_elemtype_check
	global vol_elemtype_check
	set surf_elemtype_check 0
	set vol_elemtype_check 0

	# Look for Elements with custom ElemTypes

	# Reset Automatic Conditions from previous executions
	GiD_Process Mescape Meshing MeshCriteria DefaultMesh Lines 1:end
	GiD_Process Mescape Meshing MeshCriteria DefaultMesh Surfaces 1:end escape

	if { $surf_elemtype_check > 0 } {
		GiD_Process Mescape Meshing ElemType Default Surfaces 1:end escape
		set surf_elemtype_check 0
	}
	if { $vol_elemtype_check > 0 } {
		GiD_Process Mescape Meshing ElemType Default Volumes 1:end escape
		set vol_elemtype_check 0
	}

	# End Reset Block
	
	# Volume Model Parts
	set volumelist [GiD_Geometry list volume 1:]
	
	# Assign Volume Parts

	# Generate lists and assign conditions
	# End Volume Block

	# Surface Model Parts
	set surfacelist [GiD_Geometry list surface 1:]

	if {[GiD_Info Geometry NumVolumes]>0} {
		set boundary [findboundary surface]
		# 3D Boundary Section
	}
	
	# Inherited Surface Parts
	
	# Edit Surface Parts
	
	# Generate lists and assign conditions
	# End Surface Block

	# Line Model Parts
	set linelist [GiD_Geometry list line 1:]

	if {[GiD_Info Geometry NumVolumes]==0} {
		set boundary [findboundary line]
		# 2D Boundary Section
	}

	# Inherited Line Parts

	# Edit Line Parts

	# Generate lists and assign conditions
	# End Line Block

	# Point Model Parts
	set pointlist [GiD_Geometry list point 1:]
	
	# Inherited Point Parts
	
	# Edit Point Parts

	# Generate lists and assign conditions
	# End Point Block

	# Select Elements from Options

	# Assign Non-Default Mesh Criteria to Entities
	# End Meshing Block
}

proc AfterMeshGeneration {fail} {
	# After Mesh Generation
}

proc domainsize { } {
	# Returns 2 for 2D problems and 3 for 3d problems
	set bbox [GiD_Info Layers -bbox]
	set bbox [lindex $bbox 0]
	set dx [expr {[lindex $bbox 0]-[lindex $bbox 3]}]
	set dy [expr {[lindex $bbox 1]-[lindex $bbox 4]}]
	set dz [expr {[lindex $bbox 2]-[lindex $bbox 5]}]
	set depth [expr {$dx*$dy*$dz}]
	if {$depth == 0} {set ds 2} else {set ds 3}
	return $ds
}

# Makes all of boundary lines' normals point inwards or outwards
proc alignlinenormals {Direction} {
	switch $Direction {
	Inwards {set wrong_way "DIFF1ST"}
	Outwards {set wrong_way "SAME1ST"} 
	default {puts "Unknown Direction, line normals not aligned"}
	}

	set surfacelist [GiD_Geometry list surface 1:]

	# For each surface, we look for boundary lines oriented in the wrong direction
	foreach surface $surfacelist {
		set surfaceinfo [GiD_Info list_entities surfaces $surface]
		set numpos [lsearch $surfaceinfo "NumLines:"]
		set numlines [lindex $surfaceinfo [expr {$numpos +1}]]
		for {set i 0} {$i < $numlines} {incr i} {
			set orient [lindex $surfaceinfo [expr {$numpos+5+4*$i}]]
			if {[string compare $orient $wrong_way]==0} {
				# If the normal is pointing in the wrong direction, 
				# Check if it's a contour line
				set linenum [lindex $surfaceinfo [expr {$numpos+3+4*$i}]]
				set lineinfo [GiD_Info list_entities lines $linenum]
				#set highpos [lsearch $surfinfo "HigherEntity:"]
				set higherentities [lindex $lineinfo 4]
				if {$higherentities==1} {
					# If its in the contour, switch its normal
					GiD_Process Mescape Utilities SwapNormals Lines Select $linenum
				}
			}
		}
	}
}

# Makes all of boundary surfaces' normals point inwards or outwards
proc alignsurfnormals {Direction} {
	switch $Direction {
	Inwards {set wrong_way "DIFF1ST"}
	Outwards {set wrong_way "SAME1ST"} 
	default {puts "Unknown Direction, surface normals not aligned"}
	}

	set volumelist [GiD_Geometry list volume 1:]

	# For each volume, we look for face surfaces with oriented in the wrong direction
	foreach volume $volumelist {
		set volumeinfo [GiD_Info list_entities volumes $volume]
		set numpos [lsearch $volumeinfo "NumSurfaces:"]
		set numsurf [lindex $volumeinfo [expr {$numpos +1 }]]
		for {set i 0} {$i < $numsurf} {incr i} {
			set orient [lindex $volumeinfo [expr {$numpos+5+4*$i}]]
			if {[string compare $orient $wrong_way]==0} {
				# If the normal is pointing in the wrong direction,
				# Check if it's a contour surface
				set surfnum [lindex $volumeinfo [expr {$numpos+3+4*$i}]]
				set surfinfo [GiD_Info list_entities surfaces $surfnum]
				#set highpos [lsearch $surfinfo "HigherEntity:"]
				set higherentities [lindex $surfinfo 4]
				if {$higherentities==1} {
					# If its in the contour, switch its normal
					GiD_Process Mescape Utilities SwapNormals Surfaces Select $surfnum
				}
			}
		}
	}
}

# Assigns a condition (or a model part) from a line to its ends
proc cond_linetopoint {Condition} {
	set linename "line_$Condition"
	set pointname "point_$Condition"
	# Identify points ending lines with "Condition" set and apply "Condition" to them
	set condlines [GiD_Info conditions $linename geometry]
	foreach line $condlines {
		set id [lindex $line 1]
		set new_auto [expr {[lindex $line 3] + 1}]
		set val [lrange $line 4 end]
		set l_info [GiD_Info list_entities lines $id]
		set ppos [lsearch $l_info "Points:"]
		set points {}
		lappend points [lindex $l_info [expr {$ppos+1}]]
		lappend points [lindex $l_info [expr {$ppos+2}]]
		foreach point $points {
			# Look for an existing condition of this type assigned to this entity. If there is one, store its "automatic" rating, if there isn't, store a high enough value.
			set current_info [lindex [GiD_Info conditions $pointname geometry $point] 0]
			if {[string compare [lindex $current_info 0] E]== 0} {
				set current_auto [lindex $current_info 3]
			} else {
				set current_auto 10
			}
			if {$new_auto <= $current_auto} {
				set valuestring "$new_auto $val"
				GiD_AssignData Condition $pointname points $valuestring $point
				# Assigns the conditions with "Auto" in the "automatic" question
				# so it will be recognised as a condition assigned via tcl
			}
		}
	}
}

# Assigns a condition (or a model part) from a surface to its boundary lines
proc cond_surfacetoline {Condition} {
	set surfname "surface_$Condition"	
	set linename "line_$Condition"
	# Identify lines limiting surfaces with "Condition" set and apply "Condition" to them
	set condsurf [GiD_Info conditions $surfname geometry]
	foreach surf $condsurf {
		set linelist {}
		set id [lindex $surf 1]
		set new_auto [expr {[lindex $surf 3] + 1}]
		set val [lrange $surf 4 end]
		set s_info [GiD_Info list_entities surfaces $id]
		set d [lsearch $s_info "NumLines:"]
		set numlines [lindex $s_info [expr {$d+1}]]
		for {set i 0} {$i < $numlines} {incr i} {
			set line [lindex $s_info [expr {$d+3+4*$i}]]
			lappend linelist $line
		}
		foreach line $linelist {
			# Look for an existing condition of this type assigned to this entity. If there is one, store its "automatic" rating, if there isn't, store a high enough value.
			set current_info [lindex [GiD_Info conditions $linename geometry $line] 0]
			if {[string compare [lindex $current_info 0] E]== 0} {
				set current_auto [lindex $current_info 3]
			} else {
				set current_auto 100
			}
			if {$new_auto <= $current_auto} {
				set valuestring "$new_auto $val"
				GiD_AssignData Condition $linename lines $valuestring $line
				# Assigns the conditions with "Auto" in the "automatic" question
				# so it will be recognised as a condition assigned via tcl
			}
		}
	}
}

# Assigns a condition (or a model part) from a volume to its faces
proc cond_volumetosurface {Condition} {
	set volname "volume_$Condition"
	set surfname "surface_$Condition"
	# Identify surfaces limiting volumes with "Condition" set and apply "Condition" to them
	set condvol [GiD_Info conditions $volname geometry]
	foreach vol $condvol {
		set surflist {}
		set id [lindex $vol 1]
		set new_auto [expr {[lindex $vol 3] + 1}]
		set val [lrange $vol 4 end]
		set v_info [GiD_Info list_entities volumes $id]
		set d [lsearch $v_info "NumSurfaces:"]
		set numsurfs [lindex $v_info [expr {$d+1}]]
		for {set i 0} {$i < $numsurfs} {incr i} {
			set surf [lindex $v_info [expr {$d+3+4*$i}]]
			lappend surflist $surf
		}
		foreach surf $surflist {
			# Look for an existing condition of this type assigned to this entity. If there is one, store its "automatic" rating, if there isn't, store a high enough value.
			set current_info [lindex [GiD_Info conditions $surfname geometry $surf] 0]
			if {[string compare [lindex $current_info 0] E]== 0} {
				set current_auto [lindex $current_info 3]
			} else {
				set current_auto 100
			}
			if {$new_auto <= $current_auto} {
				set valuestring "$new_auto $val"
				GiD_AssignData Condition $surfname surfaces $valuestring $surf
				# Assigns the conditions with "Auto" in the "automatic" question
				# so it will be recognised as a condition assigned via tcl
			}
		}
	}
}

# Assigns all materials from a surface to its boundary lines (unless the line has a material given by the user)
proc assign_materials {} {
	set materiallist [GiD_Info materials]
	set entitylist [list "volume" "surface" "line" "point"]
	for {set i 0} {$i < 3} {incr i} {
		set big_ent [lindex $entitylist $i]
		set small_ent [lindex $entitylist [expr {$i+1}]]
		set ent_list [GiD_Geometry list $big_ent 1:]
		foreach entity $ent_list {
			set entinfo [GiD_Info list_entities "${big_ent}s" $entity]
			set mat_pos [lsearch $entinfo "material:"]
			set mat_num [lindex $entinfo [expr {$mat_pos +1}]]
			if {$mat_num > 0} {
				set contour_list {}
				switch $big_ent {
				"volume" {
					set d [lsearch $entinfo "NumSurfaces:"]
					set numsurfs [lindex $entinfo [expr {$d+1}]]
					for {set j 0} {$j < $numsurfs} {incr j} {
						set surf [lindex $entinfo [expr {$d+3+4*$j}]]
						lappend contour_list $surf
					}
				}
				"surface" {
					set d [lsearch $entinfo "NumLines:"]
					set numlines [lindex $entinfo [expr {$d+1}]]
					for {set j 0} {$j < $numlines} {incr j} {
						set line [lindex $entinfo [expr {$d+3+4*$j}]]
						lappend contour_list $line
					}
				}
				"line" {
					set ppos [lsearch $entinfo "Points:"]
					lappend contour_list [lindex $entinfo [expr {$ppos+1}]]
					lappend contour_list [lindex $entinfo [expr {$ppos+2}]]
				}}
				set mat_name [lindex $materiallist [expr {$mat_num-1}]]
				GiD_AssignData Material $mat_name ${small_ent}s $contour_list
			}
		}
	}
}

# Assigns a condition defined via a model part. args contains the ids of the entities to edit
proc condfrompart {Part Condition mode entity args} {
	set entlist [GiD_Info conditions ${entity}_${Part} geometry $args]
	switch $mode {
	"always" {set modestring ""}
	"only2D" {set modestring "2D"}
	"only3D" {set modestring "3D"}
	}
	foreach item $entlist {
		set id [lindex $item 1]
		set auto [lindex $item 3]
		set new_auto [expr {$auto +1}]
		# Verify if this entity already has a condition of the same type assigned
		# If there is one, check its priority compared to the one we want to assign
		set current_info [lindex [GiD_Info conditions ${entity}_${Condition} geometry $id] 0]
		if {[string compare [lindex $current_info 0] E]== 0} {
			set current_auto [lindex $current_info 3]
		} else {
			set current_auto 100
		}
		if {$new_auto <= $current_auto} {
			set exp [format {BEGIN%s%s (.*?) END%s%s} $modestring $Condition $modestring $Condition]
			# matches everything (non-greedy) between 'BEGIN...COND_NAME ' and ' END...COND_NAME'
			set val [regexp -linestop -all -inline $exp $item]
			set val [lindex $val 1]
			set val "$new_auto $val"
			GiD_AssignData Condition ${entity}_${Condition} ${entity}s $val $id
		}
	}
}

# Mesh all entities with "Element" assigned (use to mesh with boundary elements)
proc meshelement {Element Entity} {
	set elemlist [createlist $Entity $Element]
	foreach id $elemlist {
	GiD_Process Mescape Meshing MeshCriteria Mesh ${Entity}s $id
	}
	GiD_Process Mescape
}

# Mesh all entities with "Elemtype" elements
proc meshtype {Element Entity Elemtype} {
	global surf_elemtype_check
	global vol_elemtype_check
	if {[string compare $Elemtype "None"]!= 0} {
		set elemlist [createlist $Entity $Element]
		if {[llength $elemlist] > 0} {
			switch $Entity {
				surface {set surf_elemtype_check 1}
				volume {set vol_elemtype_check 1}
			}
		}
		GiD_Process Mescape Meshing ElemType $Elemtype "$elemlist"
		GiD_Process Mescape
	}
}

# Check if elements with a non default ElemType have been used
proc check_elemtype {Element Entity Elemtype} {
	if {[string compare $Elemtype "None"]!= 0} {
		switch $Entity {
			surface {
				global surf_elemtype_check
				set elemlist [createlist surface $Element]
				if {[llength $elemlist] > 0} {set surf_elemtype_check 1}
			} volume {
				global vol_elemtype_check
				set elemlist [createlist volume $Element]
				if {[llength $elemlist] > 0} {set vol_elemtype_check 1}
			}
		}
	}
}

# Make a list containing all boundary entities
# use entity=line for models made of surfaces and entity=surface for models made of volumes
proc findboundary {entity} {
	set boundarylist {}
	# Generate some names
	set Entity [string toupper $entity 0 0]
	set entities [format "%ss" $entity]
	
	# Get the number of the last entity
	set instruction [format "MaxNum%ss" $Entity]
	set Max [GiD_Info Geometry $instruction]
	
	# Generate a list containing all entities and record their id and number of HigerEntities
	set EntityList [GiD_Info list_entities $entities 1:$Max]
	set candidates [regexp -all -inline {Num: ([0-9]*) HigherEntity: ([0-9]*)} $EntityList]

	# Find ids of entities with exactly 1 HigherEntity (this means they are in the boundary)
	for {set i 1} {$i < [llength $candidates] } {incr i 3} {
		set j [expr {$i + 1}]
		if {[lindex $candidates $j] == 1} {lappend boundarylist [lindex $candidates $i]}
	}
	return $boundarylist
}

# Given an entity (point, line, surface, volume) and Model Part or a Condition, returns a list of all entities with that condition assigned
proc createlist {entity Part} {
	set partlist {}
	foreach item [GiD_Info conditions ${entity}_${Part} geometry] {
		lappend partlist [lindex $item 1]
	}
	return $partlist
}

# Unassigns automatically assigned GiD Conditions (Model Parts, Conditions, Elements) from given entity types
proc cleanautomatic {Condition  args} {
	foreach entity $args {
		set autolist {}
		set infolist [GiD_Info conditions ${entity}_${Condition} geometry]
		foreach item $infolist {
			set id [regexp -inline {^E ([0-9]*) - ([0-9]*)} $item]
			# If it's "automatic" value is >0 (its always 0 for user-assigned data), store its id
			if {[lindex $id 2] > 0} {lappend autolist [lindex $id 1]}
		}
		GiD_UnAssignData Condition ${entity}_${Condition} ${entity}s $autolist
	}
}

# Assign a given single-value condition to all entities that don't have it assigned. The value given is taken from General Data
proc assigndefault {Condition Property args} {
	set defvalue [GiD_AccessValue get gendata $Property]
	foreach entity $args {
		set applylist {}
		set entitylist [GiD_Geometry list $entity 1:]
		set condlist [createlist $entity $Condition]
		foreach item $entitylist {
			if {[lsearch $condlist $item] == -1} then {
				lappend applylist $item
			}
		}
		set newval "1 $defvalue"
		GiD_AssignData Condition ${entity}_${Condition} ${entity}s $newval $applylist
	}
}

## Mesh Point elements
#proc create_point_elems { Condition } {
#	set materiallist [GiD_Info materials]
#	foreach node [GiD_Info conditions point_${Condition} mesh] {
#		set nodenum [lindex $node 1]
#		set pointnum [lindex $node 4]
#		set pointinfo [GiD_Info list_entities Points $pointnum]
#		regexp {material: ([0-9]+)} $pointinfo {} matnum
#		if { $matnum > 0 } {
#			set matname [lindex $materiallist [expr {$matnum-1}]]
#			GiD_Mesh create element append point 1 $nodenum $matname
#		} else {
#			GiD_Mesh create element append point 1 $nodenum
#		}
#	}
#}

# Mesh Point elements
proc create_point_elems { Condition } {
	set materiallist [GiD_Info materials]
	foreach node [GiD_Info conditions point_${Condition} mesh] {
		set nodenum [lindex $node 1]
		set pointnum [lindex $node 4]
		set pointinfo [GiD_Info list_entities Points $pointnum]
		regexp {material: ([0-9]+)} $pointinfo {} matnum
		if { $matnum > 0 } {
			set matname [lindex $materiallist [expr {$matnum-1}]]
			GiD_Mesh create element append point 1 $nodenum $matname
		} else {
			GiD_Mesh create element append point 1 $nodenum
		}
		# Assign an element condition over the point
		set elemnum [GiD_Info Mesh MaxNumElements]
		# Note that the element we just created will be the one with the highest number
		GiD_AssignData Condition element_${Condition} body_elements "" $elemnum
	}
}

# Assign an element from an option
proc assign_element_choice { Option entity args } {
	set user_elements {}
	foreach possible_element $args {
		lappend user_elements [createlist $entity $possible_element]
	}
	foreach entinfo [GiD_Info conditions ${entity}_${Option} geometry] {
		set entnum [lindex $entinfo 1]
		if {[lsearch $user_elements $entnum] < 0} {
			set new_auto [expr {[lindex $entinfo 3]+1}]
			set choice [lindex $entinfo 4]
			if {[string compare $choice "Use_Default"] == 0} {set choice [GiD_AccessValue get gendata ${Option}]}
			# Get the default values for the condition and drop the 'Automatic' one
			set defvals [GiD_Info Condition ${entity}_${choice}]
			set valstring ""
			for {set i 5} {$i < [llength $defvals] } {incr i 2} {
				set val [lindex $defvals $i]
				set valstring "$valstring $val"
			}
			set valstring "$new_auto $valstring"
			GiD_AssignData Condition ${entity}_${choice} ${entity}s $valstring $entnum
			# It is always assigned using default values (but elements should't have values, so it shouldn't be a problem)
		}
	}
}

# Generate a report listing several "interesting" facts about the problem
proc cond_report {} {
# Check that there is a model
set have_geometry [GiD_Info Geometry NumPoints]
if {$have_geometry>0} then {
	# Conditions
	set cond_books [GiD_Info conditions BOOKS]
	# list: {bookname1 1 bookname2 2 bookname3 3 ...}
	set book_names {}
	set used_books {}
	set book_strings {}
	# Generate 3 lists containing: book names; 0 if the book is used in the project, 1 otherwise; a string to display in the message box (only if the book was used)
	for {set i 0} {$i < [llength $cond_books]} {incr i 2} {
		set bookname [lindex $cond_books $i]
		set book_names [concat $book_names [list "$bookname"]]
		set used_books [concat $used_books 0]
		set book_strings [concat $book_strings [list "Used $bookname:\n"]]
	}
	# For each entity: Generate a list of conditions and add all used ones to our string
	set namelist { {point ovpnt {point_(.*)} } {line ovline {line_(.*)} } {surface ovsurf {surface_(.*)} } {volume ovvol {volume_(.*)} } }
	foreach lst $namelist {
		set entity [lindex $lst 0]
		set ov [lindex $lst 1]
		set regstring [lindex $lst 2]
		set condlist [GiD_Info conditions $ov]
		# Check each condition: if it was used in the model, add it to its book's string
		foreach condition $condlist {
			set cond_info [GiD_Info conditions $condition geometry]
			set one [lindex $cond_info 0]
			set found ""
			regexp {^E [0-9]* - [0-9]*} $one found
			# Try to find condition information: if some is found, "found" variable will be ocverwritten (with something beginning with 'E')
			if {[string compare $found ""] != 0} then {
				# condition was assigned
				set shortname ""
				regexp $regstring $condition everything shortname
				# shortname is the name of the condition without the point_, line_, ... prefix
				# find the book containing the condition and check if a condition with the same shortname has been used before (ex: for line_CONDITION, check if point_CONDITION was used)
				# The following line avoids breaking this code if the condition does not follow the problem type generator's standard naming convention of entity_Name
				if {[string compare $shortname ""] == 0} {set shortname $condition}
				set condbook [GiD_Info conditions $condition BOOK]
				if {[string compare $condbook "Default"]!=0} {
					set book_index [lsearch $book_names $condbook]
					set bookstring [lindex $book_strings $book_index]
					set cond_used unused
					regexp $shortname $bookstring cond_used
					if { [string compare $cond_used unused] == 0} {
						# The condition is not in the string, add it
						lset book_strings $book_index "${bookstring}\t$shortname\n"
						lset used_books $book_index 1
					}
				}
			}
		}
	}
	# Build a result string, to display in the message window
	set result_string ""
	set have_conditions 0
	for {set i 0} {$i < [llength $used_books]} {incr i} {
		if {[lindex $used_books $i] == 1} then {
			set str [lindex $book_strings $i]
			set result_string "${result_string}\n$str"
			set have_conditions 1
		}
	}
	if {$have_conditions == 0} {set result_string "No conditions have been assigned yet\n"}
	# Materials
	# Find out if we are in a 2D or a 3D problem
	if {[GiD_Info Geometry NumVolumes]>0} then {
		set maxent volume
		set otherents "Surfaces, lines and points"
		set check_list [list volume surface line]
	} else {
		set maxent surface
		set otherents "Lines and points"
		set check_list [list surface line]
	}
	# Check if lower entities will recive a material automatically
	if {[GiD_AccessValue get gendata Transfer_materials_to_lower_entities] == 1} {
		set ent_list  [GiD_Geometry list $maxent 1:]
		# Check that all highest level entities have an assigned material
		set have_mat 1
		foreach ent $ent_list {
			set ent_info [GiD_Info list_entities ${maxent}s $ent]
			regexp {material: ([0-9]+)} $ent_info everything matnumber
			if {$matnumber == 0} then {set have_mat 0}
		}
		if {$have_mat == 0} then {
			set result_string "${result_string}\nSome ${maxent}s don't have any material assigned\n"
		} else {
			set result_string "${result_string}\nAll ${maxent}s have a material assigned\n"
		}
		set result_string "${result_string}$otherents will inherit a material from their parent ${maxent}s\n"
	} else {
		foreach ent_name $check_list {
			set ent_list  [GiD_Geometry list $ent_name 1:]
			set have_mat 1
			foreach ent $ent_list {
				set ent_info [GiD_Info list_entities ${maxent}s $ent]
				regexp {material: ([0-9]+)} $ent_info everything matnumber
				if {$matnumber == 0} then {set have_mat 0}
			}
			if {$have_mat == 0} then {
				set result_string "${result_string}\nSome ${ent_name}s don't have any material assigned\n"
			} else {
				set result_string "${result_string}\nAll ${ent_name}s have a material assigned\n"
			}
		}
	}
	# Mesh
	set meshinfo [GiD_Info Mesh]
	if {$meshinfo == 0} then {set result_string "${result_string}\nA mesh must be generated"}
} else { set result_string "You should draw your model first" }
# Display the results in a window
	set w .gid.win_example

	InitWindow $w "Model Status" Modelstatus ""	"" 1

	frame $w.top
    label $w.top.results -text $result_string -justify left

    frame $w.bottom
    button $w.bottom.ok -text "OK" -height 1 -width 5 -command "destroy $w"

    pack $w.top.results -side left

    pack $w.bottom.ok -side left -anchor center

    pack $w.top -padx 15 -pady 5
    pack $w.bottom -side top -padx 6 -ipady 2
}


#to be used by TKWIDGET to easily add a select file button
#e.g.
#QUESTION: your_question
#VALUE: your_filename
#TKWIDGET: GidUtils::TkwidgetGetFilenameButton

proc TkwidgetFilePath { event args } {
    global tkwidgedprivfilepaths tkwidgedprivfilepathbuttons
    switch $event {
	INIT {
	    set PARENT [lindex $args 0]
	    upvar [lindex $args 1] ROW
	    set GDN [lindex $args 2]
	    set STRUCT [lindex $args 3]
	    set QUESTION [lindex $args 4]
	    #initialize variable to current field value
	    set tkwidgedprivfilepaths($QUESTION) [DWLocalGetValue $GDN $STRUCT $QUESTION]
	    #set entry $PARENT.e$ROW
	    set entry ""
	    foreach item [grid slaves $PARENT -row [expr $ROW-1]] {
		if { [winfo class $item] == "Entry"  || [winfo class $item] == "TEntry" } {
		    #assumed that it is the only entry of this row
		    set entry $item
		    break
		}
	    }
	    if { $entry != "" } {
		set tkwidgedprivfilepathbuttons($QUESTION) [button $PARENT.bfolder$QUESTION \
		        -image [GetImage "folder.gif"] \
		        -command [list GidUtils::_GetFilenameCmd tkwidgedprivfilepaths($QUESTION) $entry 0]]
		grid $tkwidgedprivfilepathbuttons($QUESTION) -row [expr $ROW-1] -column 2 -sticky w
		grid configure $entry -sticky ew
	    }
	}
	SYNC {
	    #set GDN [lindex $args 0]
	    #set STRUCT [lindex $args 1]
	    #set QUESTION [lindex $args 2]
	    #DWLocalSetValue $GDN $STRUCT $QUESTION $tkwidgedprivfilepaths($QUESTION)
	}
	DEPEND {
	    #set GDN [lindex $args 0]
	    #set STRUCT [lindex $args 1]
	    set QUESTION [lindex $args 2]
	    set ACTION [lindex $args 3]
	    #set value [lindex $args 4]
	    if { [info exists tkwidgedprivfilepathbuttons($QUESTION)] && \
		    [winfo exists $tkwidgedprivfilepathbuttons($QUESTION)] } {
		if { $ACTION == "HIDE" } {
		    grid remove $tkwidgedprivfilepathbuttons($QUESTION)
		} else { 
		    #RESTORE
		    grid $tkwidgedprivfilepathbuttons($QUESTION)
		}
	    } else {

	    }
	}
	CLOSE {
	    array unset tkwidgedprivfilepaths
	    array unset tkwidgedprivfilepathbuttons
	}
	default {
	    return [list ERROR [= "Unexpected tkwidget event"]]
	}
    }
    return ""
}

proc getmatnum {matname} {
	return [ expr {[lsearch [GiD_Info materials] $matname] + 1} ]
}

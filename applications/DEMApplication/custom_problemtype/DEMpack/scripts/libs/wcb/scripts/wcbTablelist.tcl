#==============================================================================
# Contains Wcb procedures for tablelist widgets.
#
# Copyright (c) 1999-2010  Csaba Nemethi (E-mail: csaba.nemethi@t-online.de)
#==============================================================================

#
# Private procedure
# =================
#

#------------------------------------------------------------------------------
# wcb::tablelistWidgetCmd
#
# Processes the Tcl command corresponding to a tablelist widget w with
# registered callbacks.  In this procedure, the execution of the commands
# activate, selection set, selection clear, activatecell, cellselection set,
# and cellselection clear is preceded by calls to the corresponding
# before-callbacks and followed by calls to the corresponding after-callbacks,
# in the global scope.
#------------------------------------------------------------------------------
proc wcb::tablelistWidgetCmd {w argList} {
    set orig [list ::_$w]

    set argCount [llength $argList]
    if {$argCount == 0} {
	# Let Tk report the error
	return [uplevel 2 $orig $argList]
    }

    set option [lindex $argList 0]
    set opLen [string length $option]
    set opArgs [lrange $argList 1 end]

    if {[string compare $option "activate"] == 0} {
	if {$argCount == 2} {
	    return [wcb::processCmd $w activate activate $opArgs]
	} else {
	    # Let Tk report the error
	    return [uplevel 2 $orig $argList]
	}

    } elseif {[string first $option "selection"] == 0 && $opLen >= 3} {
	set selOption [lindex $opArgs 0]

	if {[string first $selOption "set"] == 0} {
	    if {$argCount == 3 || $argCount == 4} {
		set selOpArgs [lrange $opArgs 1 end]
		return [wcb::processCmd $w selset "selection set" \
			$selOpArgs]
	    } else {
		# Let Tk report the error
		return [uplevel 2 $orig $argList]
	    }
	} elseif {[string first $selOption "clear"] == 0} {
	    if {$argCount == 3 || $argCount == 4} {
		set selOpArgs [lrange $opArgs 1 end]
		return [wcb::processCmd $w selclear "selection clear" \
			$selOpArgs]
	    } else {
		# Let Tk report the error
		return [uplevel 2 $orig $argList]
	    }
	} else {
	    return [uplevel 2 $orig $argList]
	}

    } elseif {[string first $option "activatecell"] == 0 && $opLen >= 9} {
	if {$argCount == 2} {
	    return [wcb::processCmd $w activatecell activatecell $opArgs]
	} else {
	    # Let Tk report the error
	    return [uplevel 2 $orig $argList]
	}

    } elseif {[string first $option "cellselection"] == 0 && $opLen >= 3} {
	set selOption [lindex $opArgs 0]

	if {[string first $selOption "set"] == 0} {
	    if {$argCount == 3 || $argCount == 4} {
		set selOpArgs [lrange $opArgs 1 end]
		return [wcb::processCmd $w cellselset "cellselection set" \
			$selOpArgs]
	    } else {
		# Let Tk report the error
		return [uplevel 2 $orig $argList]
	    }
	} elseif {[string first $selOption "clear"] == 0} {
	    if {$argCount == 3 || $argCount == 4} {
		set selOpArgs [lrange $opArgs 1 end]
		return [wcb::processCmd $w cellselclear "cellselection clear" \
			$selOpArgs]
	    } else {
		# Let Tk report the error
		return [uplevel 2 $orig $argList]
	    }
	} else {
	    return [uplevel 2 $orig $argList]
	}

    } else {
	return [uplevel 2 $orig $argList]
    }
}

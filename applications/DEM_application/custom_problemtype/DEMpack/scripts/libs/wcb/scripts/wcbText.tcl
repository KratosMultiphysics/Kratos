#==============================================================================
# Contains Wcb procedures for text and ctext widgets.
#
# REMARK: Everything stated below for text widgets is valid for ctext widgets,
#	  too.
#
# Copyright (c) 1999-2010  Csaba Nemethi (E-mail: csaba.nemethi@t-online.de)
#==============================================================================

#
# Namespace initialization
# ========================
#

namespace eval wcb {
    #
    # Some regexp patterns:
    #
    if {$tk_version >= 8.1} {
	variable alphaOrNlPat	{^[[:alpha:]\n]*$}
	variable digitOrNlPat	{^[[:digit:]\n]*$}
	variable alnumOrNlPat	{^[[:alnum:]\n]*$}
    } else {
	# Ugly because of the \n:
	variable alphaOrNlPat	"^\[A-Za-z\n]*$"
	variable digitOrNlPat	"^\[0-9\n]*$"
	variable alnumOrNlPat	"^\[A-Za-z0-9\n]*$"
    }
}

#
# Simple before-insert callback routines for text widgets
# =======================================================
#

#------------------------------------------------------------------------------
# wcb::checkStrsForRegExp
#
# Checks whether the strings to be inserted into the text widget w, contained
# in the list args of the form "string ?tagList string tagList ...?", are
# matched by the regular expression exp; if not, it cancels the insert
# operation.
#------------------------------------------------------------------------------
proc wcb::checkStrsForRegExp {exp w idx args} {
    foreach {str tagList} $args {
	if {![regexp -- $exp $str]} {
	    cancel
	    return ""
	}
    }
}

#------------------------------------------------------------------------------
# wcb::checkStrsForAlpha
#
# Checks whether the strings to be inserted into the text widget w, contained
# in the list args of the form "string ?tagList string tagList ...?", are
# alphabetic; if not, it cancels the insert operation.
#------------------------------------------------------------------------------
proc wcb::checkStrsForAlpha {w idx args} {
    variable alphaOrNlPat
    checkStrsForRegExp $alphaOrNlPat $w $idx $args
}

#------------------------------------------------------------------------------
# wcb::checkStrsForNum
#
# Checks whether the strings to be inserted into the text widget w, contained
# in the list args of the form "string ?tagList string tagList ...?", are
# numeric; if not, it cancels the insert operation.
#------------------------------------------------------------------------------
proc wcb::checkStrsForNum {w idx args} {
    variable digitOrNlPat
    checkStrsForRegExp $digitOrNlPat $w $idx $args
}

#------------------------------------------------------------------------------
# wcb::checkStrsForAlnum
#
# Checks whether the strings to be inserted into the text widget w, contained
# in the list args of the form "string ?tagList string tagList ...?", are
# alphanumeric; if not, it cancels the insert operation.
#------------------------------------------------------------------------------
proc wcb::checkStrsForAlnum {w idx args} {
    variable alnumOrNlPat
    checkStrsForRegExp $alnumOrNlPat $w $idx $args
}

#------------------------------------------------------------------------------
# wcb::convStrsToUpper
#
# Replaces the strings to be inserted into the text widget w, contained in the
# list args of the form "string ?tagList string tagList ...?", with their
# uppercase equivalents.
#------------------------------------------------------------------------------
proc wcb::convStrsToUpper {w idx args} {
    set n 1
    foreach {str tagList} $args {
	replace $n $n [string toupper $str]
	incr n 2
    }
}

#------------------------------------------------------------------------------
# wcb::convStrsToLower
#
# Replaces the strings to be inserted into the text widget w, contained in the
# list args of the form "string ?tagList string tagList ...?", with their
# lowercase equivalents.
#------------------------------------------------------------------------------
proc wcb::convStrsToLower {w idx args} {
    set n 1
    foreach {str tagList} $args {
	replace $n $n [string tolower $str]
	incr n 2
    }
}

#
# Private procedure
# =================
#

#------------------------------------------------------------------------------
# wcb::textWidgetCmd
#
# Processes the Tcl command corresponding to a text widget w with registered
# callbacks.  In this procedure, the execution of the commands insert, delete,
# and mark set insert is preceded by calls to the corresponding before-
# callbacks and followed by calls to the corresponding after-callbacks, in the
# global scope.
#------------------------------------------------------------------------------
proc wcb::textWidgetCmd {w argList} {
    set orig [list ::_$w]

    set argCount [llength $argList]
    if {$argCount == 0} {
	# Let Tk report the error
	return [uplevel 2 $orig $argList]
    }

    set option [lindex $argList 0]
    set opLen [string length $option]
    set opArgs [lrange $argList 1 end]

    if {[string first $option "insert"] == 0 && $opLen >= 3} {
	if {$argCount >= 2} {
	    return [wcb::processCmd $w insert insert $opArgs]
	} else {
	    # Let Tk report the error
	    return [uplevel 2 $orig $argList]
	}

    } elseif {[string first $option "delete"] == 0 && $opLen >= 3} {
	if {$argCount == 2 || $argCount == 3} {
	    return [wcb::processCmd $w delete delete $opArgs]
	} else {
	    # Let Tk report the error
	    return [uplevel 2 $orig $argList]
	}

    } elseif {[string first $option "replace"] == 0} {
	if {$argCount >= 4} {
	    return [wcb::processCmd $w replace replace $opArgs]
	} else {
	    # Let Tk report the error
	    return [uplevel 2 $orig $argList]
	}

    } elseif {[string first $option "mark"] == 0} {
	set markOption [lindex $opArgs 0]

	if {[string first $markOption "set"] == 0} {
	    if {$argCount == 4 &&
		[string compare [lindex $opArgs 1] "insert"] == 0} {
		set markOpArgs [lrange $opArgs 2 end]
		return [wcb::processCmd $w motion "mark set insert" \
			$markOpArgs]
	    } else {
		return [uplevel 2 $orig $argList]
	    }
	} else {
	    return [uplevel 2 $orig $argList]
	}

    } elseif {[string first $option "tag"] == 0} {
	set tagOption [lindex $opArgs 0]

	if {[string first $tagOption "add"] == 0} {
	    if {$argCount >= 4 &&
		[string compare [lindex $opArgs 1] "sel"] == 0} {
		set selOpArgs [lrange $opArgs 2 end]
		return [wcb::processCmd $w selset "tag add sel" \
			$selOpArgs]
	    } else {
		return [uplevel 2 $orig $argList]
	    }
	} elseif {[string first $tagOption "remove"] == 0} {
	    if {$argCount >= 4 &&
		[string compare [lindex $opArgs 1] "sel"] == 0} {
		set selOpArgs [lrange $opArgs 2 end]
		return [wcb::processCmd $w selclear "tag remove sel" \
			$selOpArgs]
	    } else {
		return [uplevel 2 $orig $argList]
	    }
	} else {
	    return [uplevel 2 $orig $argList]
	}

    } else {
	return [uplevel 2 $orig $argList]
    }
}

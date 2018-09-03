#==============================================================================
# Contains Wcb procedures for Tk or tile entry, BWidget Entry, Tk or tile
# spinbox, and tile combobox widgets.
#
# REMARK: Everything stated below for entry widgets is valid for tile entry and
#         BWidget Entry widgets, too.  Similarly, everything stated below for
#         spinbox widgets is valid for tile spinbox widgets, too.
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
	variable alphaPat	{^[[:alpha:]]*$}
	variable digitPat	{^[[:digit:]]*$}
	variable alnumPat	{^[[:alnum:]]*$}
    } else {
	variable alphaPat	{^[A-Za-z]*$}
	variable digitPat	{^[0-9]*$}
	variable alnumPat	{^[A-Za-z0-9]*$}
    }
}

#
# Utility procedures for entry, spinbox, and tile combobox widgets
# ================================================================
#

#------------------------------------------------------------------------------
# wcb::changeEntryText
#
# Replaces the text of the entry, spinbox, or tile combobox widget w with the
# string str, by using the delete and insert operations.  If one of these
# subcommands is canceled by some before-callback then the procedure keeps the
# original entry, spinbox, or tile combobox string and returns 0, otherwise it
# returns 1.
#------------------------------------------------------------------------------
proc wcb::changeEntryText {w str} {
    set oldStr [$w get]
    set oldPos [$w index insert]

    $w delete 0 end
    if {[canceled $w delete]} {
	return 0
    }

    $w insert 0 $str
    if {[canceled $w insert]} {
	$w insert 0 $oldStr
	set result 0
    } else {
	set result 1
    }
    $w icursor $oldPos
    return $result
}

#------------------------------------------------------------------------------
# wcb::postInsertEntryLen
#
# Returns the length of the text that would be contained in the entry, spinbox,
# or tile combobox widget w after inserting the string str.
#------------------------------------------------------------------------------
proc wcb::postInsertEntryLen {w str} {
    return [expr {[$w index end] + [string length $str]}]
}

#------------------------------------------------------------------------------
# wcb::postInsertEntryText
#
# Returns the text that would be contained in the entry, spinbox, or tile
# combobox widget w after inserting the string str before the character
# indicated by the index idx.
#------------------------------------------------------------------------------
proc wcb::postInsertEntryText {w idx str} {
    set oldText [$w get]
    set idx [$w index $idx]

    append newText [string range $oldText 0 [expr {$idx - 1}]] \
		   $str \
		   [string range $oldText $idx end]
    return $newText
}

#------------------------------------------------------------------------------
# wcb::postDeleteEntryText
#
# Returns the text that would be contained in the entry, spinbox, or tile
# combobox widget w after deleting the range of characters starting with the
# index given by from and stopping just before the one given by the first
# element of args (if any).
#------------------------------------------------------------------------------
proc wcb::postDeleteEntryText {w from args} {
    set first [$w index $from]

    if {[llength $args] == 0} {
	set last $first
    } else {
	set to [lindex $args 0]
	set last [expr {[$w index $to] - 1}]
    }

    return [string replace [$w get] $first $last]
}

#
# Simple before-insert callback routines for
# entry, spinbox, and tile combobox widgets
# ==========================================
#

#------------------------------------------------------------------------------
# wcb::checkStrForRegExp
#
# Checks whether the string str to be inserted into the entry, spinbox, or tile
# combobox widget w is matched by the regular expression exp; if not, it
# cancels the insert operation.
#------------------------------------------------------------------------------
proc wcb::checkStrForRegExp {exp w idx str} {
    if {![regexp -- $exp $str]} {
	cancel
    }
}

#------------------------------------------------------------------------------
# wcb::checkStrForAlpha
#
# Checks whether the string str to be inserted into the entry, spinbox, or tile
# combobox widget w is alphabetic; if not, it cancels the insert operation.
#------------------------------------------------------------------------------
proc wcb::checkStrForAlpha {w idx str} {
    variable alphaPat
    checkStrForRegExp $alphaPat $w $idx $str
}

#------------------------------------------------------------------------------
# wcb::checkStrForNum
#
# Checks whether the string str to be inserted into the entry, spinbox, or tile
# combobox widget w is numeric; if not, it cancels the insert operation.
#------------------------------------------------------------------------------
proc wcb::checkStrForNum {w idx str} {
    variable digitPat
    checkStrForRegExp $digitPat $w $idx $str
}

#------------------------------------------------------------------------------
# wcb::checkStrForAlnum
#
# Checks whether the string str to be inserted into the entry, spinbox, or tile
# combobox widget w is alphanumeric; if not, it cancels the insert operation.
#------------------------------------------------------------------------------
proc wcb::checkStrForAlnum {w idx str} {
    variable alnumPat
    checkStrForRegExp $alnumPat $w $idx $str
}

#------------------------------------------------------------------------------
# wcb::convStrToUpper
#
# Replaces the string str to be inserted into the entry, spinbox, or tile
# combobox widget w with its uppercase equivalent.
#------------------------------------------------------------------------------
proc wcb::convStrToUpper {w idx str} {
    replace 1 1 [string toupper $str]
    return ""
}

#------------------------------------------------------------------------------
# wcb::convStrToLower
#
# Replaces the string str to be inserted into the entry, spinbox, or tile
# combobox widget w with its lowercase equivalent.
#------------------------------------------------------------------------------
proc wcb::convStrToLower {w idx str} {
    replace 1 1 [string tolower $str]
    return ""
}

#
# Further before-insert callback routines for
# entry, spinbox, and tile combobox widgets
# ===========================================
#

#------------------------------------------------------------------------------
# wcb::checkEntryForInt
#
# Checks whether the text contained in the entry, spinbox, or tile combobox
# widget w after inserting the string str before the character indicated by the
# index idx would represent (the starting part of) an integer number; if not,
# it cancels the insert operation.
#------------------------------------------------------------------------------
proc wcb::checkEntryForInt {w idx str} {
    set newText [postInsertEntryText $w $idx $str]
    if {![regexp {^[+-]?[0-9]*$} $newText]} {
	cancel
    }
}

#------------------------------------------------------------------------------
# wcb::checkEntryForUInt
#
# Checks whether the text contained in the entry, spinbox, or tile combobox
# widget w after inserting the string str before the character indicated by the
# index idx would represent (the starting part of) an unsigned integer no
# greater than max; if not, it cancels the insert operation.  The value * for
# max means: no upper bound.
#------------------------------------------------------------------------------
proc wcb::checkEntryForUInt {max w idx str} {
    set newText [postInsertEntryText $w $idx $str]
    if {![regexp {^[0-9]*$} $newText]} {
	cancel
    } elseif {[string compare $max *] != 0} {
	scan $newText "%d" val
	if {$val > $max} {
	    cancel
	}
    }
}

#------------------------------------------------------------------------------
# wcb::checkEntryForReal
#
# Checks whether the text contained in the entry, spinbox, or tile combobox
# widget w after inserting the string str before the character indicated by the
# index idx would represent (the starting part of) a real number; if not, it
# cancels the insert operation.
#------------------------------------------------------------------------------
proc wcb::checkEntryForReal {w idx str} {
    set newText [postInsertEntryText $w $idx $str]
    if {![regexp {^[+-]?[0-9]*\.?[0-9]*([0-9]\.?[eE][+-]?[0-9]*)?$} $newText]} {
	cancel
    }
}

#------------------------------------------------------------------------------
# wcb::checkEntryForFixed
#
# Checks whether the text contained in the entry, spinbox, or tile combobox
# widget w after inserting the string str before the character indicated by the
# index idx would represent (the starting part of) a real number with at most
# cnt digits after the decimal point; if not, it cancels the insert operation.
# The value * for cnt means: unlimited number of digits after the decimal
# point.
#------------------------------------------------------------------------------
proc wcb::checkEntryForFixed {cnt w idx str} {
    set pattern {^[+-]?[0-9]*\.?}
    if {[string compare $cnt "*"] == 0} {
	append pattern {[0-9]*$}
    } else {
	for {set n 0} {$n < $cnt} {incr n} {
	    append pattern {[0-9]?}
	}
	append pattern $
    }

    set newText [postInsertEntryText $w $idx $str]
    if {![regexp $pattern $newText]} {
	cancel
    }
}

#------------------------------------------------------------------------------
# wcb::checkEntryLen
#
# Checks whether the length of the text contained in the entry, spinbox, or
# tile combobox widget w after inserting the string str would be greater than
# len; if yes, it cancels the insert operation.
#------------------------------------------------------------------------------
proc wcb::checkEntryLen {len w idx str} {
    if {[postInsertEntryLen $w $str] > $len} {
	cancel
    }
}

#
# Private procedure
# =================
#

#------------------------------------------------------------------------------
# wcb::entryWidgetCmd
#
# Processes the Tcl command corresponding to an entry, spinbox, or tile
# combobox widget w with registered callbacks.  In this procedure, the
# execution of the commands insert, delete, and icursor is preceded by calls to
# the corresponding before-callbacks and followed by calls to the corresponding
# after-callbacks, in the global scope.
#------------------------------------------------------------------------------
proc wcb::entryWidgetCmd {w argList} {
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
	if {$argCount == 3} {
	    return [wcb::processCmd $w insert insert $opArgs]
	} else {
	    # Let Tk report the error
	    return [uplevel 2 $orig $argList]
	}

    } elseif {[string first $option "delete"] == 0} {
	if {$argCount == 2 || $argCount == 3} {
	    return [wcb::processCmd $w delete delete $opArgs]
	} else {
	    # Let Tk report the error
	    return [uplevel 2 $orig $argList]
	}

    } elseif {[string first $option "icursor"] == 0 && $opLen >= 2} {
	if {$argCount == 2} {
	    return [wcb::processCmd $w motion icursor $opArgs]
	} else {
	    # Let Tk report the error
	    return [uplevel 2 $orig $argList]
	}

    } else {
	return [uplevel 2 $orig $argList]
    }
}

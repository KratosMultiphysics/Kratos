#==============================================================================
# Contains common Wcb procedures.
#
# Copyright (c) 1999-2010  Csaba Nemethi (E-mail: csaba.nemethi@t-online.de)
#==============================================================================

#
# Namespace initialization
# ========================
#

namespace eval wcb {
    #
    # Bind some cleanup operations to the <Destroy> event
    # for all widgets having the binding tag WcbCleanup
    #
    bind WcbCleanup <Destroy> {
	wcb::cleanup %W
    }
}

#
# Basic procedures
# ================
#

#------------------------------------------------------------------------------
# wcb::callback
#
# Retrieves, sets, or removes the callbacks for the widget w, the argument
# when, and the command corresponding to option.  when can be "before" or
# "after", and option can take one of the following values:
#
#   - "insert", "delete", or "motion",	for a Tk or tile entry, BWidget Entry,
#					Tk or tile spinbox, tile combobox,
#				 	text, or ctext widget;
#   - "replace",			for a text or ctext widget;
#   - "activate",			for a listbox or tablelist widget;
#   - "selset" or "selclear",		for a listbox, tablelist, text, or
#					ctext widget;
#   - "activatecell", "cellselset", or
#     "cellselclear",			for a tablelist widget.
#
# If no arguments after the option parameter are specified, then the procedure
# just returns the current before- or after-callback list, respectively, for
# the given widget operation.
#
# Otherwise, if at least one of the arguments following the option parameter is
# a nonempty string, then:
#
#   - if called for the first time for this widget with at least one nonempty
#     argument following the option parameter, then it replaces the Tcl command
#     w with a new procedure in which the execution of the widget operations
#     associated with the above values of option is preceded by invocations of
#     the corresponding before-callbacks and followed by calls to the
#     corresponding after-callbacks, in the global scope;
#   - it sets the callback list to the one built from these arguments and
#     returns the new list.
#
# Otherwise (i.e. if all arguments following the option parameter are empty),
# then the procedure unregisters all the corresponding callbacks for the given
# widget and returns an empty string.
#
# When a callback is invoked, the name of the original Tcl command for the
# widget w as well as the command arguments are automatically appended to it as
# parameters.
#------------------------------------------------------------------------------
proc wcb::callback {w when option args} {
    if {![winfo exists $w]} {
	return -code error "bad window path name \"$w\""
    }

    if {[string first $when "before"] == 0} {
	set when before
    } elseif {[string first $when "after"] == 0} {
	set when after
    } else {
	return -code error "bad argument \"$when\": must be before or after"
    }

    if {[catch {fullCallbackOpt $w $option} result] != 0} {
	return -code error $result
    }
    set option $result

    variable data
    if {[llength $args] == 0} {
	if {[info exists data($w-$when-$option)]} {
	    return $data($w-$when-$option)
	} else {
	    return {}
	}
    } elseif {[areAllEmptyStrings $args]} {
	catch {unset data($w-$when-$option)}
	return ""
    } else {
	switch [winfo class $w] {
	    Entry -
	    TEntry -
	    Spinbox -
	    TSpinbox -
	    TCombobox {
		set widgetCmd entryWidgetCmd
	    }
	    Listbox {
		set widgetCmd listboxWidgetCmd
	    }
	    Tablelist {
		set widgetCmd tablelistWidgetCmd
	    }
	    Text -
	    Ctext {
		set widgetCmd textWidgetCmd
	    }
	}
	redefWidgetCmd $w $widgetCmd
	return [set data($w-$when-$option) $args]
    }
}

#------------------------------------------------------------------------------
# wcb::cbappend
#
# Appends the arguments represented by args to the current before- or after-
# callback list, respectively, for the given widget operation.
#------------------------------------------------------------------------------
proc wcb::cbappend {w when option args} {
    if {[catch {callback $w $when $option} result] != 0} {
	return -code error $result
    }

    eval lappend result $args
    return [eval [list callback $w $when $option] $result]
}

#------------------------------------------------------------------------------
# wcb::cbprepend
#
# Prepends the arguments represented by args to the current before- or after-
# callback list, respectively, for the given widget operation.
#------------------------------------------------------------------------------
proc wcb::cbprepend {w when option args} {
    if {[catch {callback $w $when $option} result] != 0} {
	return -code error $result
    }

    set result [eval [list linsert $result 0] $args]
    return [eval [list callback $w $when $option] $result]
}

#------------------------------------------------------------------------------
# wcb::cancel
#
# If invoked from a before-callback for a widget command, this procedure
# cancels the execution of that command and of the remaining callbacks, and
# calls script in the global scope.
#------------------------------------------------------------------------------
proc wcb::cancel {{script bell}} {
    variable data
    set data(canceled-[info level 1]) 1

    if {[string compare $script ""] != 0} {
	uplevel #0 $script
    }
}

#------------------------------------------------------------------------------
# wcb::canceled
#
# Returns 1 if the most recent invocation of the widget operation correspondig
# to w and option has been aborted by some before-callback, and 0 otherwise.
#------------------------------------------------------------------------------
proc wcb::canceled {w option} {
    if {![winfo exists $w]} {
	return -code error "bad window path name \"$w\""
    }

    if {[catch {fullCallbackOpt $w $option} result] != 0} {
	return -code error $result
    }
    set option $result

    variable data
    if {[info exists data($w-canceled-$option)]} {
	return $data($w-canceled-$option)
    } else {
	return 0
    }
}

#------------------------------------------------------------------------------
# wcb::extend
#
# If invoked from a before-callback for a widget command, this procedure
# appends the values given in args to the argument list of that command.  The
# new argument list will be passed to the remaining callbacks for that command,
# too.
#------------------------------------------------------------------------------
proc wcb::extend args {
    variable data
    upvar 0 data(args-[info level 1]) var
    eval lappend var $args
}

#------------------------------------------------------------------------------
# wcb::replace
#
# If invoked from a before-callback for a widget command, this procedure
# replaces the arguments having the indices first through last of that command
# with the values given in args.  The new argument list will be passed to the
# remaining callbacks for that command, too.  The arguments are numbered from 0.
#------------------------------------------------------------------------------
proc wcb::replace {first last args} {
    variable data
    upvar 0 data(args-[info level 1]) var
    set var [eval [list lreplace $var $first $last] $args]
}

#------------------------------------------------------------------------------
# wvb::pathname
#
# Returns the path name of the widget corresponding to the Tcl command origCmd
# (which is supposed to be of the form "::_pathName").
#------------------------------------------------------------------------------
proc wcb::pathname origCmd {
    return [string range $origCmd 3 end]
}

#
# Private procedures
# ==================
#

#------------------------------------------------------------------------------
# wcb::cleanup
#
# Unregisters all callbacks defined for w and deletes the Tcl command w.
#------------------------------------------------------------------------------
proc wcb::cleanup w {
    variable data
    foreach when {before after canceled} {
	foreach option {insert delete motion activate selset selclear} {
	    catch {unset data($w-$when-$option)}
	}
    }

    catch {rename ::$w ""}
    catch {rename ::_$w ""}		;# necessary for tablelist widgets
}

#------------------------------------------------------------------------------
# wcb::fullCallbackOpt
#
# Returns the full callback option corresponding to the possibly abbreviated
# option opt.
#------------------------------------------------------------------------------
proc wcb::fullCallbackOpt {w opt} {
    set opLen [string length $opt]

    switch [winfo class $w] {
	Entry -
	TEntry -
	Spinbox -
	TSpinbox -
	TCombobox {
	    if {[string first $opt "insert"] == 0} {
		set opt insert
	    } elseif {[string first $opt "delete"] == 0} {
		set opt delete
	    } elseif {[string first $opt "motion"] == 0} {
		set opt motion
	    } else {
		return -code error \
		       "bad option \"$opt\": must be insert, delete, or motion"
	    }
	}

	Listbox {
	    if {[string first $opt "activate"] == 0} {
		set opt activate
	    } elseif {[string first $opt "selset"] == 0 && $opLen >= 4} {
		set opt selset
	    } elseif {[string first $opt "selclear"] == 0 && $opLen >= 4} {
		set opt selclear
	    } else {
		return -code error \
		       "bad option \"$opt\": must be activate, selset, or\
			selclear"
	    }
	}

	Tablelist {
	    if {[string compare $opt "activate"] == 0} {
		set opt activate
	    } elseif {[string first $opt "selset"] == 0 && $opLen >= 4} {
		set opt selset
	    } elseif {[string first $opt "selclear"] == 0 && $opLen >= 4} {
		set opt selclear
	    } elseif {[string first $opt "activatecell"] == 0 && $opLen >= 9} {
		set opt activatecell
	    } elseif {[string first $opt "cellselset"] == 0 && $opLen >= 8} {
		set opt cellselset
	    } elseif {[string first $opt "cellselclear"] == 0 && $opLen >= 8} {
		set opt cellselclear
	    } else {
		return -code error \
		       "bad option \"$opt\": must be activate, selset,\
			selclear, activatecell, cellselset, or cellselclear"
	    }
	}

	Text -
	Ctext {
	    if {[string first $opt "insert"] == 0} {
		set opt insert
	    } elseif {[string first $opt "delete"] == 0} {
		set opt delete
	    } elseif {[string first $opt "replace"] == 0} {
		set opt replace
	    } elseif {[string first $opt "motion"] == 0} {
		set opt motion
	    } elseif {[string first $opt "selset"] == 0 && $opLen >= 4} {
		set opt selset
	    } elseif {[string first $opt "selclear"] == 0 && $opLen >= 4} {
		set opt selclear
	    } else {
		return -code error \
		       "bad option \"$opt\": must be insert, delete, replace,\
		        motion, selset, or selclear"
	    }
	}

	default {
	    return -code error \
		   "window \"$w\" is not a Tk or tile entry,\
		    BWidget Entry, Tk or tile spinbox, tile combobox,\
		    listbox, tablelist, text, or ctext widget"
	}
    }

    return $opt
}

#------------------------------------------------------------------------------
# wcb::areAllEmptyStrings
#
# Returns 1 if all elements of the list lst are empty strings and 0 otherwise.
#------------------------------------------------------------------------------
proc wcb::areAllEmptyStrings lst {
    foreach elem $lst {
	if {[string compare $elem ""] != 0} {
	    return 0
	}
    }

    return 1
}

#------------------------------------------------------------------------------
# wcb::redefWidgetCmd
#
# Renames the Tcl command w to _w, builds a new widget procedure w that invokes
# cmd, and appends WcbCleanup to the list of binding tags of the widget w.
#------------------------------------------------------------------------------
proc wcb::redefWidgetCmd {w cmd} {
    if {[catch {rename ::$w ::_$w}] != 0} {
	return ""
    }

    #
    # If the command within the catch below returns an error, we
    # will substitute all occurrences of ::_$w with $w.  To this
    # end we need a version of $w in which the characters |, *, +,
    # ?, (, ., ^, $, \, [, {, }, ,, :, =, and ! are escaped, and
    # another version in which the characters & and \ are escaped
    # (to suppress the special processing of &, \0, \1, ..., \9).
    #
    regsub -all {\||\*|\+|\?|\(|\.|\^|\$|\\|\[|\{|\}|\,|\:|\=\!} $w {\\\0} w1
    regsub -all {&|\\} $w {\\\0} w2

    proc ::$w args [format {
	if {[catch {wcb::%s %s $args} result] == 0} {
	    return $result
	} else {
	    regsub -all -- %s $result %s result
	    return -code error $result
	}
    } $cmd [list $w] [list ::_$w1] [list $w2]]

    bindtags $w [linsert [bindtags $w] end WcbCleanup]
}

#------------------------------------------------------------------------------
# wcb::processCmd
#
# Invokes the before-callbacks registered for the widget w and the command
# corresponding to wcbOp, then executes the script "::_w cmdOp argList", and
# finally invokes the after-callbacks.
#------------------------------------------------------------------------------
proc wcb::processCmd {w wcbOp cmdOp argList} {
    variable data
    set data($w-canceled-$wcbOp) 0
    set orig [list ::_$w]

    #
    # Invoke the before-callbacks
    #
    if {[info exists data($w-before-$wcbOp)]} {
	foreach cb $data($w-before-$wcbOp) {
	    if {[string compare $cb ""] != 0} {
		#
		# Set the two array elements that might be changed
		# by cancel, extend, or replace, invoked (directly
		# or indirectly) from within the callback
		#
		set cb [eval list $cb]
		set cbScript [concat $cb $orig $argList]
		set data(canceled-$cbScript) 0
		set data(args-$cbScript) $argList

		#
		# Invoke the callback and get the new
		# values of the two array elements
		#
		uplevel #0 $cb $orig $argList
		set data($w-canceled-$wcbOp) $data(canceled-$cbScript)
		set argList $data(args-$cbScript)

		#
		# Remove the two array elements
		#
		unset data(canceled-$cbScript)
		unset data(args-$cbScript)

		if {$data($w-canceled-$wcbOp)} {
		    return ""
		}
	    }
	}
    }

    #
    # Execute the widget command
    #
    eval $orig $cmdOp $argList

    #
    # Invoke the after-callbacks
    #
    if {[info exists data($w-after-$wcbOp)]} {
	foreach cb $data($w-after-$wcbOp) {
	    if {[string compare $cb ""] != 0} {
		uplevel #0 $cb $orig $argList
	    }
	}
    }
}

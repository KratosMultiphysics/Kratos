#==============================================================================
# Wcb package index file.
#
# Copyright (c) 1999-2011  Csaba Nemethi (E-mail: csaba.nemethi@t-online.de)
#==============================================================================

#
# Regular package:
#
package ifneeded wcb 3.4 [list source [file join $dir wcb.tcl]]

#
# Alias:
#
package ifneeded Wcb 3.4 { package require -exact wcb 3.4 }

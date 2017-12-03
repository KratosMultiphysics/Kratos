# AniGif Package written in pure Tcl/Tk
#
# anigif.tcl v1.3 2002-09-09 (c) 2001-2002 Ryan Casey
#
# AniGif is distributed under the same license as Tcl/Tk.  As of
# AniGif 1.3, this license is applicable to all previous versions.
#
# ###############################  USAGE  #################################
#
#  ::anigif::anigif FILENAME WIDGET INDEX
#    FILENAME: appropriate path and file to use for the animated gif
#    WIDGET:   a label widget to place the animated gif into
#    DELAY:    how long to wait before loading next image (Default: 100ms)
#    INDEX:    what image to begin on (first image is 0) (Default: 0)
#
#  ::anigif::stop WIDGET
#  ::anigif::restart WIDGET INDEX
#    INDEX:    defaults to next index in loop
#  ::anigif::destroy WIDGET
#
#  NOTES:
#    There is currently a -zoom and -subsample hack to keep transparency.
#    Anigif does not handle interlaced gifs properly.  The image will look
#      distorted.
#    A delay of 0 renders as fast as possible, per the GIF specification.
#      This is currently set to 40 ms to approximate the IE default.
#    If you experience a problem with a compressed gif, try uncompressing
#      it. Search the web for gifsicle.    
#
# ############################## HISTORY #################################
#
#  1.3: Fixed error in disposal flag handling.
#       Added handling for non-valid comment/graphic blocks.
#       Searches for actual loop control block.  If it extists, loops.
#       Added more comments.
#  1.2: Now handles single playthrough gifs or gifs with partial images
#       Fixed bug in delay time (unsigned int was being treated as signed)
#  1.1: Reads default timing instead of 100 ms or user-defined.
#       You can no longer set the delay manually.
#  1.0: Moved all anigif variables to the anigif namespace
#  0.9: Initial release
# 

namespace eval anigif {

  proc anigif2 {w list delay {idx 0}} {
    if { ![winfo exists $w] } {
#Cleanup
#???
      destroy $w
      return
    } else {
      if { $idx >= [llength $list]  } {
        set idx 0
        if { [set ::anigif::${w}(repeat)] == 0} {
          # Non-repeating GIF
          ::anigif::stop $w
          return
        }
      } 
      set dispflag [lindex [set ::anigif::${w}(disposal)] $idx]
      #W "dispfalg $dispflag"
      switch "$dispflag" {
        "000" {
              # Do nothing
            }
        "001" {
              # Do not dispose
            }
        "100" {
              # Restore to background
              [set ::anigif::${w}(curimage)] blank
            }
        "101" {
              # Restore to previous - not supported
              # As recommended, since this is not supported, it is set to blank
              [set ::anigif::${w}(curimage)] blank
            }
        default { #W "no match: $dispflag" }
      }
      [set ::anigif::${w}(curimage)] copy [lindex $list $idx] -subsample 2 2
      if { [lindex $delay $idx] == 0 } {
        ::anigif::stop $w
        return
      }
      update
      set ::anigif::${w}(asdf) "::anigif::anigif2 $w [list $list]"
      set ::anigif::${w}(loop) [after [lindex $delay $idx] "[set ::anigif::${w}(asdf)] [list $delay] [expr {$idx + 1}]"]
      set ::anigif::${w}(idx) [incr idx]
    }
  }

  proc anigif {fnam w {idx 0}} {
    set n 0
    set images {}
    set delay {}
    set fnam [string map  {\\ /} $fnam]
    set fin [open $fnam r]
    fconfigure $fin -translation binary
    set data [read $fin [file size $fnam]]
    close $fin

    # Find Loop Record
    set start [string first "\x21\xFF\x0B" $data]
    if {$start < 0} {
      set repeat 0
    } else {
      set repeat 1
    }

    # Find Control Records
    set start [string first "\x21\xF9\x04" $data]
#catch "image create photo xpic$n$w \
#      -file ${fnam} \
#      -format \{gif89 -index $n\}" err
#      W $err
#      
    while {![catch "image create photo xpic$n$w \
      -file ${fnam} \
      -format \{gif89 -index $n\}"]} {
        set stop [string first "\x00" $data [expr {$start + 1}]]
        #W "index $n"
        if {$stop < $start} {
          break
        }
        set record [string range $data $start $stop]
        binary scan $record @4c1 thisdelay
        if {[info exists thisdelay]} {

          # Change to unsigned integer
          set thisdelay [expr {$thisdelay & 0xFF}];

          binary scan $record @2b3b3b1b1 -> disposalval userinput transflag

          lappend images pic$n$w
          image create photo pic$n$w
          pic$n$w copy xpic$n$w -zoom 2 2
          image delete xpic$n$w
          lappend disposal $disposalval

          # Convert hundreths to thousandths for after
          set thisdelay [expr {$thisdelay * 10}]

          # If 0, set to fastest (25 ms min to seem to match browser default)
          if {$thisdelay == 0} {set thisdelay 40}

          lappend delay $thisdelay
          unset thisdelay

          incr n
        }

        if {($start >= 0) && ($stop >= 0)} {
          set start [string first "\x21\xF9\x04" $data [expr {$stop + 1}]]
        } else {
          break
        }
    }
    set ::anigif::${w}(repeat) $repeat
    set ::anigif::${w}(delay) $delay
    #set disposal "100"
    set ::anigif::${w}(disposal) $disposal
    set ::anigif::${w}(curimage) [image create photo]
    [set ::anigif::${w}(curimage)] blank
    [set ::anigif::${w}(curimage)] copy pic0${w} -subsample 2 2
    $w configure -image [set ::anigif::${w}(curimage)]

    anigif2 $w $images $delay $idx
  }

  proc stop {w} {
    catch {
      after cancel [set ::anigif::${w}(loop)]
    }
  }

  proc restart {w {idx -1}} {
    if {$idx == -1} {
      if { [lindex ::anigif::${w}(delay) $idx] < 0 } {
        set idx 0
      } else {
        set idx [set ::anigif::${w}(idx)]
      }
    }
    catch {
      ::anigif::stop $w
      eval "[set ::anigif::${w}(asdf)] [list [set ::anigif::${w}(delay)]] $idx"
    }
  }

  proc destroy {w} {
    catch {
      ::anigif::stop $w
      set wlength [string length $w]
      foreach imagename [image names] {
        if {[string equal [string range $imagename [string first "." $imagename] end] $w]} {
          image delete $imagename
        }
      }
     unset ::anigif::${w}
    }
  }

}

package provide anigif 1.3

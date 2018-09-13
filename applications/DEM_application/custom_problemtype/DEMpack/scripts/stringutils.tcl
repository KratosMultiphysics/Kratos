###############################################################################
#
#    NAME: kstringutils.tcl
#
#    PURPOSE: Utilities procedures to work with string data
#
#    QUANTECH - DEVELOPMENT DEPARTMENT

#    AUTHOR : G. Socorro
#
#    CREATED AT: 01/11/09
#
#    LAST MODIFICATION : create the base source code
#
#    VERSION : 0.1
#
#    HISTORY:
#
#     0.1-01/11/09-G. Socorro, create the base source code
#
###############################################################################

#******************************************************************************
#                  BEGIN KStrUtils NAMESPACE
#******************************************************************************

namespace eval KStrUtils {
    
}

# List of procedure that is inside the KStrUtils namespace
# 1. 
# 2.

proc KStrUtils::PutUnderScoreSpace { str } {
    regsub -all { } $str {_} str
    return $str
}

proc KStrUtils::PutUnderScoreSpaceToList { list } {
    
    set newlist ""
    foreach elem $list {
	regsub -all { } $elem {_} clem
	;# WarnWinText "elem:$elem => clen:$clem\n"
	lappend newlist $clem 
    }
    
    return $newlist
}

proc KStrUtils::DeleteUnderScoreSpace { str } {
    regsub -all {_} $str { } str
    return $str
} 

proc KStrUtils::DeleteUnderScoreSpaceToList { list } {
    
    set newlist ""
    foreach elem $list {
	regsub -all {_} $elem { } clem
	lappend newlist $clem 
    }
    
    return $newlist
}

proc KStrUtils::WriteAsciiCode {} {
    
    set max 1000
    for {set ii 0} {$ii<$max} {incr ii} {
	set aa [format %c $ii]
	;# puts "ii:$ii val:$aa"
	;# WarnWinText "ii:$ii val:$aa"
    }
}

proc KStrUtils::StringCap {str {idx -1}} {
    # Capitalize a string, or one char in it
    # Arguments:
    #   str input string
    #   idx  idx to capitalize
    # Results:  Returns string with specified capitalization
    
    if {$i>-1} {
	if {[string length $str]>$i} {
	    return $str
	} else {
	}
    } else {
	return [string toupper [string index $str 0]][string tolower \
							  [string range $str 1 end]]
    }
}


proc KStrUtils::StringReverse {s} {
    # Reverses input string
    # Arguments:  s : input string to reverse
    # Returns:  string with chars reversed
    
    if {[set ii [string len $s]]} {
      while {$ii} {append r [string index $s [incr ii -1]]}
      return $r
    }
}
  

proc KStrUtils::Obfuscate { s } {
    # If I describe it, it ruins it...
    # Arguments:  s -> input string
    # Returns: output

    if {[set len [string len $s]]} {
      set ii -1
      while {[incr ii]<$len} {
          set c [string index $s $ii]
          if {[regexp "\[\]\\\[ \{\}\t\n\"\]" $c]} {
            append r $c
          } else {
            scan $c %c c
            append r \\[format %0.3o $c]
          }
      }
      return $r
    }
}

proc KStrUtils::Untabify {str {tablen 8}} {
    # Removes tabs from a string, replacing with appropriate number of spaces
    # Arguments:
    #   str     input string
    #   tablen  tab length, defaults to 8
    # Returns: string sans tabs

    set out {}
    while {[set ii [string first "\t" $str]] != -1} {
      set jj [expr {$tablen-($i%$tablen)}]
      append out [string range $str 0 [incr ii -1]][format %*s $jj { }]
      set str [string range $str [incr ii 2] end]
    }
    return $out$str
}

proc KStrUtils::Tabify {str {tablen 8}} {
    # Converts excess spaces to tab chars
    # Arguments:
    #   str     input string
    #   tablen  tab length, defaults to 8
    # Returns:  string with tabs replacing excess space where appropriate

    ## We must first untabify so that \t is not interpreted to be one char
    set str [untabify $str]
    set out {}
    while {[set i [string first { } $str]] != -1} {
      ## Align i to the upper tablen boundary
      set i [expr {$i+$tablen-($i%$tablen)-1}]
      set s [string range $str 0 $i]
      if {[string match {* } $s]} {
          append out [string trimright $s { }]\t
      } else {
          append out $s
      }
      set str [string range $str [incr i] end]
    }
    return $out$str
}

proc KStrUtils::WrapLines "txt {len 75} {P \n\n} {P2 \254}" {
    #   Wraps text to a specific max line length
    # Arguments:
    #   txt           input text
    #   len           desired max line length+1, defaults to 75
    #   P       paragraph boundary chars, defaults to \n\n
    #   P2            substitute for $P while processing, defaults to \254
    #           this char must not be in the input text
    # Returns:
    #   text with lines no longer than $len, except where a single word
    #   is longer than $len chars.  does not preserve paragraph boundaries.

    regsub -all $P $txt $P2 txt
    regsub -all "\n" $txt { } txt
    incr len -1
    set out {}
    while {[string len $txt]>$len} {
      set i [string last { } [string range $txt 0 $len]]
      if {$i == -1 && [set i [string first { } $txt]] == -1} break
      append out [string trim [string range $txt 0 [incr i -1]]]\n
      set txt [string range $txt [incr i 2] end]
    }
    regsub -all $P2 $out$txt $P txt
    return $txt
}

proc KStrUtils::RomanToDecimal {x} {
    # Converts a roman numeral to decimal
    # Arguments:  x ->number in roman numeral format
    # Returns:  decimal number

    set result ""
    foreach elem {
	{ 1000      m }   { 900 cm }    
	{ 500 d }   { 400 id }    
	{ 100 c }   { 90  ic }    
	{ 50  l }    
	{ 10  x }   { 9   ix }    
	{ 5   v }   { 4   iv }    
	{ 1   i }
    } {
	set digit [lindex $elem 0]
	set roman [lindex $elem 1]
	while {$x >= $digit} {
	    append result $roman
	    incr x -$digit
	}
    }
    return $result
}

proc KStrUtils::BinaryToHexadecimal {bin} {
    #  Converts binary to hex number
    # Arguments: bin -> number in binary format
    # Returns: hexadecimal number

    ## No sanity checking is done
    array set t {
      0000 0 0001 1 0010 2 0011 3 0100 4
      0101 5 0110 6 0111 7 1000 8 1001 9
      1010 a 1011 b 1100 c 1101 d 1110 e 1111 f
    }
    set diff [expr {4-[string length $bin]%4}]
    if {$diff != 4} {
        set bin [format %0${diff}d$bin 0]
    }
    regsub -all .... $bin {$t(&)} hex
    return [subst $hex]
}

proc KStrUtils::HexadecimalToBinary {hex} {
    # Converts hex number to bin
    # Arguments: hex -> number in hex format
    # Returns: binary number (in chars, not binary format)

    array set t {
	0 0000 1 0001 2 0010 3 0011 4 0100
	5 0101 6 0110 7 0111 8 1000 9 1001
	a 1010 b 1011 c 1100 d 1101 e 1110 f 1111
	A 1010 B 1011 C 1100 D 1101 E 1110 F 1111
    }
    regsub {^0[xX]} $hex {} hex
    regsub -all . $hex {$t(&)} bin
    return [subst $bin]
}

#******************************************************************************
#                  END KStrUtils NAMESPACE
#******************************************************************************


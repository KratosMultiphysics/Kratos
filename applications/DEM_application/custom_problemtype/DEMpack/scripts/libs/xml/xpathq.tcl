###############################################################################
#
#    NAME: xpath.tcl
#
#    PURPOSE: Utilities procedures to work with XPATH query (see http://www.w3.org/TR/xpath)
#
#    QUANTECH ATZ-DEVELOPMENT DEPARTMENT
#
#    AUTHOR : G. Socorro
#
#    CREATED AT: 01/11/09
#
#    LAST MODIFICATION : create a base source code
#
#    VERSION : 0.1
#
#    HISTORY:
#
#     0.2- 23/12/09-G. Socorro, add format and str from 
#                               http://objectmix.com/tcl/34751-1-xpath-construction-tcl-2-anyone-know-better-string-handling-method.html by Ramon Ribo
#                               at 07-27-2007, 02:40 AM 
#     0.1- 01/11/09-G. Socorro, create a base source code
#
###############################################################################
# Some documentation
# 2006-05-22 From linux-magazine article issue 20
# Table 3: XPATH examples
#
# Query                Description
#  /option              The option element directly below the root node
#  //option             All elements in the document called option
#  //option[3]          The third option element
#  /table/*             All elements below table, where table must be located directly below the root node
#  //table[1]           The first table element in a document
#  //table[last()]      The last table element
#  //@colspan           All colspan attributes in a document
#  //td[@colspan]       All td elements with the attribute colspan
#  //table[@width]      All table elements that have a width attribute
#  //table[@width=690]  All table elements with a width attribute that has avalue of 690
#  //*[count(tr)=2]     All elements with two tr child nodes
#  //tr/td|th           All td and th elements contained within a tr element
#  //table//img         All img elements contained within a table-element
#  //table[1]//img[2]   Second img element in the first table
#  //status[.='hit']/.. All elements that contain a <status>hit</status> element

# From wiki: http://wiki.tcl.tk/1957
# Two notes:
# 1). The brackets '[ ]' have syntactical meaning both in tcl and in XPath expressions. Don't forget to protect the brackets in your XPath expressions.
# 2). The XPath expression //a is not the best example one could choose. The '//' (which is the abbreviation for /descendant-or-self::node()/) is one of the 
#     most expensive XPath location steps for almost all known XPath engines. It means, that the XPath engine has to scan the whole tree beneath the node. 
#     Avoiding // - of course, if possible - could amazingly speed up your XPath queries or your XSLT stylesheets. Rolf Ade.

# Create a base namespace xpathq (Xpath query)

package provide xpathq 1.0

namespace eval ::xpathq:: {
 
}

proc ::xpathq::format {string args } {
    set cmd [list format $string]
    foreach i $args {
	lappend cmd [xpath_str $i]
    }
    return [eval $cmd]
}

proc ::xpathq::String { str } {

    foreach "strList type pos" [list "" "" 0] break
    while 1 {
	switch $type {
	    "" {
		set ret [regexp -start $pos -indices {['"]} $str idxs]
                if { !$ret } {
                   lappend strList "\"[string range $str $pos end]\""
                   break
                 }
                set idx [lindex $idxs 0]
                switch -- [string index $str $idx] {
                   ' { set type apostrophe }
                   \" { set type quote }
                 }
              }
           apostrophe {
                  set ret [regexp -start $pos -indices {["]} $str idxs]
                  if { !$ret } {
		      lappend strList "\"[string range $str $pos end]\""
		      break
		  }
		set idx [lindex $idxs 0]
		lappend strList "\"[string range $str $pos [expr {$idx-1}]]\""
		set type quote
		set pos $idx
	    }
	    quote {
		set ret [regexp -start $pos -indices {[']} $str idxs]
		if { !$ret } {
		    lappend strList "'[string range $str $pos end]'"
		    break
		}
		set idx [lindex $idxs 0]
		lappend strList "'[string range $str $pos [expr {$idx-1}]]'"
		set type apostrophe
		set pos $idx
	    }
	}
    }

    if { [llength $strList] > 1 } {
	return "concat([join $strList ,])"
    } else {
	return [lindex $strList 0]
    }
}











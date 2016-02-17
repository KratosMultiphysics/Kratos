package require TclOO
package require tdom
 
# This is a metaclass, a class that defines the behavior of other classes
catch {singleton destroy}
oo::class create singleton {
    superclass oo::class
    variable object
    variable root
    unexport create ;# Doesn't make sense to have named singletons
    method new args {
        if {![info exists object]} {
           set object [next {*}$args]
        }
        return $object
    }
}
 
singleton create xmlGet {
    variable root
    
    method init {pat} {
        set root [my openxml $pat]
    }
    
    method getElements {tag val} {
        
        set nodes [$root getElementsByTagName "ElementItem"]
        set lista [list ]
        foreach n $nodes {
            if {[my cumple $n $tag $val]} {lappend lista $n}
        }
        
        return $lista
        
    }
    method getElementsPretty {tag val pn} {
        set nodes [my getElements $tag $val]
    
        set lista ""
	  foreach n $nodes {
	    append lista [[$n getElementsByTagName $pn] text] ","
	  }
	  return [string range $lista 0 end-1]
        
    }
    
    method openxml { pat} {
        set xmlfile $pat
        set f [open $xmlfile]
        set doc [dom parse [read $f]]
        close $f
        return [$doc documentElement]
    }
    
    method cumple {node tag val} {
        set item [lindex [$node getElementsByTagName $tag] 0]
        if {[$item text] == $val} {return 1}
        return 0
    }
}
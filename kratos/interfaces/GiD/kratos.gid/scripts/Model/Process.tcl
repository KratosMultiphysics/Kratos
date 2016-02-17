namespace eval Model {
# Clase Process
catch {Process destroy}
oo::class create Process {
    superclass Entity
    
    constructor {n} {
        next $n
    }
}
}

proc Model::GetProcess { id } {
    variable Processes
    
    foreach s $Processes {
        if {[$s getName] eq $id} {return $s}
    }
    
}
proc Model::ParseProcesses { doc } {
    variable Processes
    
    set ProcNodeList [$doc getElementsByTagName Process]
    foreach ProcNode $ProcNodeList {
        lappend Processes [ParseProcessNode $ProcNode]
    }
}

proc Model::ParseProcessNode { node } {
    set name [$node getAttribute n]
    
    set sl [::Model::Process new $name]
    $sl setPublicName [$node getAttribute pn]
    
    foreach attr [$node attributes] {
        $sl setAttribute $attr [$node getAttribute $attr]
    }
    
    foreach in [[$node getElementsByTagName inputs] getElementsByTagName parameter] {
        set sl [ParseInputParamNode $sl $in]
    }
    return $sl
}

proc Model::getAllProcs { } {
    variable Processes
    return $Processes
	
}


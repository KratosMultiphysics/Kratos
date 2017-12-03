##################################################################################
#   This file is common for all Kratos Applications.
#   Do not change anything here unless it's strictly necessary.
##################################################################################

namespace eval Model {
catch {Material destroy}
oo::class create Material {
    superclass Entity
    
    variable MaterialType
    
    constructor {n} {
        next $n
        variable MaterialType
        
        set MaterialType ""
    }
    method getMaterialType { } {variable MaterialType; return $MaterialType}
    method setMaterialType { mt } {variable MaterialType; set MaterialType $mt}
}
}


proc Model::ParseMaterials { doc } {
    variable Materials
    
    set MatNodeList [$doc getElementsByTagName Material]
    foreach MatNode $MatNodeList {
        lappend Materials [ParseMatNode $MatNode]
    }
}

proc Model::ParseMatNode { node } {
    set name [$node getAttribute n]
    
    set mat [::Model::Material new $name]
    $mat setPublicName [$node getAttribute n]
    $mat setMaterialType [$node getAttribute MaterialType]
    $mat setHelp [$node getAttribute help]
    
    foreach att [$node attributes] {
        $mat setAttribute $att [split [$node getAttribute $att] ","]
    }
    foreach in [[$node getElementsByTagName inputs] getElementsByTagName parameter]  {
        set mat [ParseInputParamNode $mat $in]
    }
    
    return $mat
}

proc Model::GetMaterials {args} { 
    variable Materials
    #W "Get elements $args"
    set cumplen [list ]
    foreach mat $Materials {
        if {[$mat cumple {*}$args]} { lappend cumplen $mat}
    }
    #W "Elementos buenos $cumplen"
    return $cumplen
}

proc Model::getMaterial {mid} { 
    variable Materials

    foreach mat $Materials {
        if {[$mat getName] eq $mid} { return $mat}
    }
    return ""
}


proc Model::ForgetMaterials { } {
    variable Materials
    set Materials [list ]
}

proc Model::ForgetMaterial { mid } {
    variable Materials
    set Materials2 [list ]
    foreach material $Materials {
        if {[$material getName] ne $mid} {
            lappend Materials2 $material
        }
    }
    set Materials $Materials2
}
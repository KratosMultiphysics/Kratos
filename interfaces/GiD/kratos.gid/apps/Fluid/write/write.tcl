namespace eval Fluid::write {
    # Namespace variables declaration
}

proc Fluid::write::Init { } {
    # Namespace variables inicialization
 
}

# Events
proc Fluid::write::writeModelPartEvent { } {
    write::initWriteData "FLParts" "FLMaterials"
    write::writeModelPartData
    writeProperties
    write::writeMaterials
    write::writeNodalCoordinates
    write::writeElementConnectivities
    writeConditions
    writeMeshes
}

proc Fluid::write::writeParametersEvent { } {
    write::WriteString "Project Parameters"
    write::WriteString "domain_size = [write::getValue nDim]"
    
    write::WriteString "End Project Parameters"
}

proc Fluid::write::writeCustomFilesEvent { } {

}

# MDPA Blocks
proc Fluid::write::writeProperties { } {
    # Begin Properties
    write::WriteString "Begin Properties 0"
    write::WriteString "End Properties"
    write::WriteString ""
}

proc Fluid::write::writeConditions { } {
    writeInitialConditions
    writeBoundaryConditions
}

proc Fluid::write::writeBoundaryConditions { } {
    set uniqueNames [list "FLInlet"  "FLOutletPressure"]
    foreach uniqueName $uniqueNames {
        write::writeNodalData $uniqueName
    }
    
    # IS Slip
    
}

proc Fluid::write::writeInitialConditions { } {
    #writeGravity "SMGravity"
    
    set nodalCondList {"FLIniVel" "FLIniPres"}
    foreach cond $nodalCondList {
        write::writeNodalData $cond
    }
    
}

proc Fluid::write::writeMeshes { } {
    
    
}



Fluid::write::Init

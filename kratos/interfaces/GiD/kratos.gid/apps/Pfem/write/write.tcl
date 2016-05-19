namespace eval Pfem::write {
}

proc Pfem::write::Init { } {
    # Namespace variables inicialization
}

# Project Parameters
proc Pfem::write::writeParametersEvent { } {
    write::WriteString "Project parameters file"

}

# Model Part Blocks
proc Pfem::write::writeModelPartEvent { } {
    write::initWriteData "SMParts" "SMMaterials"
    write::WriteString "Begin ModelPartData"
    write::WriteString "End ModelPartData"
    write::WriteString "Begin Properties 0"
    write::WriteString "End Properties"
    
    write::WriteString "Begin Nodes"
    set dim [write::getValue "PFEM_dimensions"]
    #if { $dim eq "2D" || $dim eq "2Da"}   
    customlib::WriteCoordinates "%5d %14.5e %14.5e %14.5e%.0s\n"           
    write::WriteString "End Nodes"
    
    write::WriteString "Begin Elements TwoStepUpdatedLagrangianVPFluidElement2D"
    set elements_conditions [list "Fluids2D"]
    #customlib::InitMaterials $elements_conditions
    #set element_formats [list {"%10d" "element" "id"} {"%10d" "element" "connectivities"} {"%10d" "material" "MID"}]
    set element_formats [list {"%10d" "element" "id"} {"%10d" "element" "connectivities"}]
    customlib::WriteConnectivities $elements_conditions $element_formats 
    write::WriteString "End Elements"
}

# Custom files (Copy python scripts, write materials file...)
proc Pfem::write::writeCustomFilesEvent { } {

}


Pfem::write::Init

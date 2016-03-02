namespace eval ::Fluid {
    # Variable declaration
    variable dir
    variable prefix
}

proc ::Fluid::Init { } {
    # Variable initialization
    variable dir
    variable prefix
    
    set dir [apps::getMyDir "Fluid"]
    set prefix FL
}

proc ::Fluid::LoadMyFiles { } {
    variable dir
    
    uplevel #0 [list source [file join $dir xml GetFromXML.tcl]]
    uplevel #0 [list source [file join $dir write write.tcl]]
    uplevel #0 [list source [file join $dir write writeProjectParameters.tcl]]
}

proc ::Fluid::DrawAutoGeom { } {
    # Draw
    GiD_Process Mescape Geometry Create Object Prism 4 1 1 0 0.0,0.0,1.0 1 3 escape
    # Groups
    GiD_Process 'Groups Create Group0 escape 'Groups Edit Rename Group0 Body escape escape escape escape Utilities EntitiesGroups Assign Body Volumes 1 escape Mescape 'Groups Create Group0 escape 'Groups Edit Rename Group0 Inlet escape escape escape escape Utilities EntitiesGroups Assign Inlet Surfaces 6 escape Mescape 'Groups Create Group0 escape 'Groups Edit Rename Group0 Outlet escape escape escape escape Utilities EntitiesGroups Assign Outlet Surfaces 1 escape Mescape 'Groups Create Group0 escape 'Groups Edit Rename Group0 Slip escape escape escape escape Utilities EntitiesGroups Assign Slip Surfaces 2 2 3 5 2 escape Mescape 'Groups Create Group0 escape 'Groups Edit Rename Group0 NoSlip escape escape escape escape Utilities EntitiesGroups Assign NoSlip Surfaces 4 escape Mescape
    
    GiD_Process "-tcl- gid_groups_conds::addF {container\[@n='Fluid'\]/condition\[@n='Parts'\]} group {n Body}" "-tcl- gid_groups_conds::addF {container\[@n='Fluid'\]/condition\[@n='Parts'\]/group\[@n='Body'\]} value {n Element pn Element dict {FractionalStep3D,Fractional Step 3D} actualize_tree 1 values {\[GetElements\]} state hidden v FractionalStep3D}" "-tcl- gid_groups_conds::addF {container\[@n='Fluid'\]/condition\[@n='Parts'\]/group\[@n='Body'\]} value {n ConstitutiveLaw pn {Fluid type} actualize_tree 1 values {\[GetConstitutiveLaws\]} dict {Newtonian,Newtonian,HerschelBulkley,Herschel Bulkley} state normal v Newtonian}" "-tcl- gid_groups_conds::addF {container\[@n='Fluid'\]/condition\[@n='Parts'\]/group\[@n='Body'\]} value {n Material pn Material editable 0 help {Choose a material from the database} values_tree {\[give_materials_list\]} value Steel actualize_tree 1 state normal v Air}" "-tcl- gid_groups_conds::addF -resolve_parametric 1 {container\[@n='Fluid'\]/condition\[@n='Parts'\]/group\[@n='Body'\]} value {n DENSITY pn Density state {\[PartParamState\]} unit_magnitude Density help {} v 1.225 units kg/m^3}" "-tcl- gid_groups_conds::addF {container\[@n='Fluid'\]/condition\[@n='Parts'\]/group\[@n='Body'\]} value {n VISCOSITY pn Viscosity state {\[PartParamState\]} unit_magnitude {} help {Expressed in Pa*s} v 1.81e-5}" "-tcl- gid_groups_conds::addF {container\[@n='Fluid'\]/condition\[@n='Parts'\]/group\[@n='Body'\]} value {n YIELD_STRESS pn {Yield stress (t)} state {\[PartParamState\]} unit_magnitude {} help {} v 0.0}" "-tcl- gid_groups_conds::addF {container\[@n='Fluid'\]/condition\[@n='Parts'\]/group\[@n='Body'\]} value {n CONSISTENCY_INDEX pn {Consistency index (k)} state {\[PartParamState\]} unit_magnitude {} help {} v 0.0}" "-tcl- gid_groups_conds::addF {container\[@n='Fluid'\]/condition\[@n='Parts'\]/group\[@n='Body'\]} value {n FLOW_INDEX pn {Flow index (?)} state {\[PartParamState\]} unit_magnitude {} help {} v 0.0}" 
    GiD_Process "-tcl- gid_groups_conds::addF {container\[@n='Fluid'\]/container\[@n='InitialConditions'\]/condition\[@n='INITIAL_VELOCITY'\]} group {n Inlet}" "-tcl- gid_groups_conds::addF -resolve_parametric 1 {container\[@n='Fluid'\]/container\[@n='InitialConditions'\]/condition\[@n='INITIAL_VELOCITY'\]/group\[@n='Inlet'\]} value {n ValX wn {INITIAL_VELOCITY _X} pn {X Value} help {Initial Velocity} unit_magnitude Velocity state {} v 0.0 units m/s}" "-tcl- gid_groups_conds::addF -resolve_parametric 1 {container\[@n='Fluid'\]/container\[@n='InitialConditions'\]/condition\[@n='INITIAL_VELOCITY'\]/group\[@n='Inlet'\]} value {n ValY wn {INITIAL_VELOCITY _Y} pn {Y Value} help {Initial Velocity} unit_magnitude Velocity state {} v 0.0 units m/s}" "-tcl- gid_groups_conds::addF -resolve_parametric 1 {container\[@n='Fluid'\]/container\[@n='InitialConditions'\]/condition\[@n='INITIAL_VELOCITY'\]/group\[@n='Inlet'\]} value {n ValZ wn {INITIAL_VELOCITY _Z} pn {Z Value} help {Initial Velocity} unit_magnitude Velocity state {} v 0.0 units m/s}" 
    GiD_Process "-tcl- gid_groups_conds::addF {container\[@n='Fluid'\]/container\[@n='InitialConditions'\]/condition\[@n='INITIAL_PRESSURE'\]} group {n Inlet}" "-tcl- gid_groups_conds::addF -resolve_parametric 1 {container\[@n='Fluid'\]/container\[@n='InitialConditions'\]/condition\[@n='INITIAL_PRESSURE'\]/group\[@n='Inlet'\]} value {n Value pn Value unit_magnitude P help {Initial Pressure} state {} v 1 units Pa}" "-tcl- gid_groups_conds::addF {container\[@n='Fluid'\]/container\[@n='BoundaryConditions'\]/condition\[@n='Inlet3D'\]} group {n Inlet}" "-tcl- gid_groups_conds::addF {container\[@n='Fluid'\]/container\[@n='BoundaryConditions'\]/condition\[@n='Inlet3D'\]/group\[@n='Inlet'\]} value {n FixX pn {X Imposed} values 1,0 help {} state {} v 1}" "-tcl- gid_groups_conds::addF {container\[@n='Fluid'\]/container\[@n='BoundaryConditions'\]/condition\[@n='Inlet3D'\]/group\[@n='Inlet'\]} value {n FixY pn {Y Imposed} values 1,0 help {} state {} v 1}" "-tcl- gid_groups_conds::addF {container\[@n='Fluid'\]/container\[@n='BoundaryConditions'\]/condition\[@n='Inlet3D'\]/group\[@n='Inlet'\]} value {n FixZ pn {Z Imposed} values 1,0 help {} state {} v 1}" "-tcl- gid_groups_conds::addF -resolve_parametric 1 {container\[@n='Fluid'\]/container\[@n='BoundaryConditions'\]/condition\[@n='Inlet3D'\]/group\[@n='Inlet'\]} value {n ValX wn {Inlet3D _X} pn {Velocity X} help {} unit_magnitude Velocity state {} v 0.0 units m/s}" "-tcl- gid_groups_conds::addF -resolve_parametric 1 {container\[@n='Fluid'\]/container\[@n='BoundaryConditions'\]/condition\[@n='Inlet3D'\]/group\[@n='Inlet'\]} value {n ValY wn {Inlet3D _Y} pn {Velocity Y} help {} unit_magnitude Velocity state {} v 0.0 units m/s}" "-tcl- gid_groups_conds::addF -resolve_parametric 1 {container\[@n='Fluid'\]/container\[@n='BoundaryConditions'\]/condition\[@n='Inlet3D'\]/group\[@n='Inlet'\]} value {n ValZ wn {Inlet3D _Z} pn {Velocity Z} help {} unit_magnitude Velocity state {\[CheckDimension 3D\]} v 0.0 units m/s}" "-tcl- gid_groups_conds::addF {container\[@n='Fluid'\]/container\[@n='BoundaryConditions'\]/condition\[@n='Outlet3D'\]} group {n Outlet}" "-tcl- gid_groups_conds::addF -resolve_parametric 1 {container\[@n='Fluid'\]/container\[@n='BoundaryConditions'\]/condition\[@n='Outlet3D'\]/group\[@n='Outlet'\]} value {n Pressure pn Pressure unit_magnitude P help {} state {} v 0 units Pa}" "-tcl- gid_groups_conds::addF {container\[@n='Fluid'\]/container\[@n='BoundaryConditions'\]/condition\[@n='Slip3D'\]} group {n Slip}" "-tcl- gid_groups_conds::addF {container\[@n='Fluid'\]/container\[@n='BoundaryConditions'\]/condition\[@n='NoSlip3D'\]} group {n NoSlip}" 

    GiD_Process Mescape Meshing Generate 1 MeshingParametersFrom=Preferences Mescape Meshing MeshView 
    
    gid_groups_conds::actualize_conditions_window
}

::Fluid::Init

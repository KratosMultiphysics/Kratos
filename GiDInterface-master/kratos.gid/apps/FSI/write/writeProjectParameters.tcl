# Project Parameters
proc FSI::write::getParametersDict { } {
   set projectParametersDict [dict create]

   # FSI section
   set FSIParametersDict [dict create]
   # Problem data
   set problemDataDict [dict create]
   set paralleltype [write::getValue ParallelType]
   dict set problemDataDict parallel_type $paralleltype
   dict set FSIParametersDict problem_data $problemDataDict
   # Solver settings
   set solverSettingsDict [dict create]
   set currentStrategyId [write::getValue FSISolStrat]
   set strategy_write_name [[::Model::GetSolutionStrategy $currentStrategyId] getAttribute "ImplementedInPythonFile"]
   dict set solverSettingsDict coupling_scheme $currentStrategyId
   dict set solverSettingsDict solver_type $strategy_write_name
   set solverSettingsDict [dict merge $solverSettingsDict [write::getSolutionStrategyParametersDict] ]
   set solverSettingsDict [dict merge $solverSettingsDict [write::getSolversParametersDict FSI] ]

   dict set solverSettingsDict mesh_solver [write::getValue FSIALEParams MeshSolver]

   set solidInterfacesList [write::GetMeshFromCondition STLoads StructureInterface2D]
   lappend solidInterfacesList {*}[write::GetMeshFromCondition STLoads StructureInterface3D]
   dict set solverSettingsDict structure_interfaces_list $solidInterfacesList

   set fluid_interface_UniqueName FluidNoSlipInterface$::Model::SpatialDimension
   set fluidInterfacesList [write::GetMeshFromCondition FLBC $fluid_interface_UniqueName]
   dict set solverSettingsDict fluid_interfaces_list $fluidInterfacesList

   dict set FSIParametersDict solver_settings $solverSettingsDict
   dict set FSIParametersDict mapper_settings [GetMappingSettingsList]

   # Structural section
   UpdateUniqueNames Structural
   apps::setActiveAppSoft Structural
   Structural::write::ApplyConfiguration

   set StructuralParametersDict [Structural::write::getParametersEvent]
   set current [dict get $StructuralParametersDict solver_settings model_import_settings input_filename]
   dict set StructuralParametersDict solver_settings model_import_settings input_filename "${current}_Structural"

   # Fluid section
   UpdateUniqueNames Fluid
   apps::setActiveAppSoft Fluid
   write::SetConfigurationAttribute parts_un FLParts
   
   set FluidParametersDict [Fluid::write::getParametersDict]
   set current [dict get $FluidParametersDict solver_settings model_import_settings input_filename]
   dict set FluidParametersDict solver_settings model_import_settings input_filename "${current}_Fluid"
   set nodalresults [dict get $FluidParametersDict output_configuration result_file_configuration nodal_results]
   lappend nodalresults "MESH_DISPLACEMENT"
   dict set FluidParametersDict output_configuration result_file_configuration nodal_results $nodalresults
   UpdateUniqueNames FSI
   apps::setActiveAppSoft FSI

   dict set projectParametersDict structure_solver_settings $StructuralParametersDict
   dict set projectParametersDict fluid_solver_settings $FluidParametersDict
   dict set projectParametersDict coupling_solver_settings $FSIParametersDict

   return $projectParametersDict
}

proc FSI::write::writeParametersEvent { } {
   set projectParametersDict [getParametersDict]
   write::SetParallelismConfiguration
   write::WriteJSON $projectParametersDict
}

proc FSI::write::UpdateUniqueNames { appid } {
    set unList [list "Results"]
    foreach un $unList {
         set current_un [apps::getAppUniqueName $appid $un]
         spdAux::setRoute $un [spdAux::getRoute $current_un]
    }
}

proc FSI::write::GetMappingSettingsList { } {
    set mappingsList [list ]

    set fluid_interface_name FluidNoSlipInterface$::Model::SpatialDimension
    set structural_interface_name StructureInterface$::Model::SpatialDimension
    set structuralInterface [lindex [write::GetMeshFromCondition STLoads $structural_interface_name] 0]
    foreach fluid_interface [[customlib::GetBaseRoot] selectNodes "[spdAux::getRoute FLBC]/condition\[@n = '$fluid_interface_name'\]/group" ] {
        set map [dict create]
        set mapper_face [write::getValueByNode [$fluid_interface selectNodes ".//value\[@n='mapper_face']"] ]
        dict set map mapper_face $mapper_face
        dict set map fluid_interface_submodelpart_name [write::getMeshId $fluid_interface_name [get_domnode_attribute $fluid_interface n]]
        dict set map structure_interface_submodelpart_name $structuralInterface
        lappend mappingsList $map
    }

    return $mappingsList
}

# {
#     "mapper_face" : "Unique" (otherwise "Positive" or "Negative")
#     "fluid_interface_submodelpart_name" : "FluidNoSlipInterface2D_FluidInterface",
#     "structure_interface_submodelpart_name" : "StructureInterface2D_StructureInterface"
# }
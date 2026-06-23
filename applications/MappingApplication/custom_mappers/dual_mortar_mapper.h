//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#pragma once

// System includes

// External includes

// Project includes
#include "mappers/mapper.h"
#include "mappers/mapper_flags.h"
#include "factories/linear_solver_factory.h"
#include "processes/simple_mortar_mapper_wrapper_process.h"

namespace Kratos
{
///@name Kratos Classes
///@{

/**
 * @class DualMortarMapper
 * @ingroup MappingApplication
 * @brief Mapper based on the (dual) mortar formulation, working between non-matching surface meshes
 * @details This mapper transfers nodal values between two non-matching interfaces using the mortar
 * method. Thanks to the use of dual Lagrange multiplier (dual shape function) spaces the mapping
 * operator is (block) diagonal, so no global system needs to be assembled and solved (an explicit
 * lumped resolution is enough), unless a linear solver is provided through "linear_solver_settings".
 *
 * It only works between surfaces (2D lines, 3D triangles/quadrilaterals); volumetric geometries are
 * not supported.
 *
 * The actual mortar operator is the proven, heavily templated core implementation
 * (@ref SimpleMortarMapperProcess), which is reused here through
 * @ref SimpleMortarMapperProcessWrapper::Create. The legacy core process
 * (Python: KratosMultiphysics.SimpleMortarMapperProcess) remains fully self-contained and
 * independent of the MappingApplication; this Mapper is a standalone port that shares the same
 * operator, so both produce identical results.
 * @author Vicente Mataix Ferrandiz
 */
template<class TSparseSpace, class TDenseSpace>
class DualMortarMapper
    : public Mapper<TSparseSpace, TDenseSpace>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of DualMortarMapper
    KRATOS_CLASS_POINTER_DEFINITION(DualMortarMapper);

    typedef Mapper<TSparseSpace, TDenseSpace> BaseType;
    typedef typename BaseType::MapperUniquePointerType MapperUniquePointerType;
    typedef typename BaseType::TMappingMatrixType TMappingMatrixType;

    /// The linear solver type used by the underlying mortar process
    typedef SimpleMortarMapperProcessWrapper::LinearSolverType LinearSolverType;
    typedef LinearSolverFactory<TSparseSpace, TDenseSpace> LinearSolverFactoryType;

    ///@}
    ///@name Life Cycle
    ///@{

    // Default constructor, needed for registration
    DualMortarMapper(
        ModelPart& rModelPartOrigin,
        ModelPart& rModelPartDestination
        ) : mrModelPartOrigin(rModelPartOrigin),
            mrModelPartDestination(rModelPartDestination)
    {}

    DualMortarMapper(
        ModelPart& rModelPartOrigin,
        ModelPart& rModelPartDestination,
        Parameters JsonParameters
        ) : mrModelPartOrigin(rModelPartOrigin),
            mrModelPartDestination(rModelPartDestination),
            mMapperSettings(JsonParameters)
    {
        KRATOS_TRY

        // We validate against our own (superset) defaults, keeping track of the user-defined entities
        const bool origin_are_conditions_is_defined = mMapperSettings.Has("origin_are_conditions");
        const bool destination_are_conditions_is_defined = mMapperSettings.Has("destination_are_conditions");
        mMapperSettings.ValidateAndAssignDefaults(GetMapperDefaultSettings());

        // Automatic detection of the entity types (reusing the core helper)
        SimpleMortarMapperProcessWrapper::DetectEntities(mrModelPartOrigin, mrModelPartDestination, mMapperSettings, origin_are_conditions_is_defined, destination_are_conditions_is_defined);

        // The (optional) linear solver for the implicit resolution
        if (mMapperSettings["linear_solver_settings"].Has("solver_type")) {
            mpLinearSolver = LinearSolverFactoryType().Create(mMapperSettings["linear_solver_settings"]);
        }

        KRATOS_CATCH("")
    }

    /// Destructor.
    ~DualMortarMapper() override = default;

    ///@}
    ///@name Operations
    ///@{

    void UpdateInterface(
        Kratos::Flags MappingOptions,
        double SearchRadius) override
    {
        // The mortar operator (search included) is rebuilt on every Map call, so there is nothing to do here
    }

    void Map(
        const Variable<double>& rOriginVariable,
        const Variable<double>& rDestinationVariable,
        Kratos::Flags MappingOptions) override
    {
        MapInternal(rOriginVariable, rDestinationVariable, MappingOptions, false);
    }

    void Map(
        const Variable< array_1d<double, 3> >& rOriginVariable,
        const Variable< array_1d<double, 3> >& rDestinationVariable,
        Kratos::Flags MappingOptions) override
    {
        MapInternal(rOriginVariable, rDestinationVariable, MappingOptions, false);
    }

    void InverseMap(
        const Variable<double>& rOriginVariable,
        const Variable<double>& rDestinationVariable,
        Kratos::Flags MappingOptions) override
    {
        MapInternal(rOriginVariable, rDestinationVariable, MappingOptions, true);
    }

    void InverseMap(
        const Variable< array_1d<double, 3> >& rOriginVariable,
        const Variable< array_1d<double, 3> >& rDestinationVariable,
        Kratos::Flags MappingOptions) override
    {
        MapInternal(rOriginVariable, rDestinationVariable, MappingOptions, true);
    }

    MapperUniquePointerType Clone(
        ModelPart& rModelPartOrigin,
        ModelPart& rModelPartDestination,
        Parameters JsonParameters) const override
    {
        return Kratos::make_unique<DualMortarMapper<TSparseSpace, TDenseSpace>>(
            rModelPartOrigin, rModelPartDestination, JsonParameters);
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "DualMortarMapper";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "DualMortarMapper";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        BaseType::PrintData(rOStream);
    }

    ///@}
private:
    ///@name Member Variables
    ///@{

    ModelPart& mrModelPartOrigin;
    ModelPart& mrModelPartDestination;
    Parameters mMapperSettings;
    typename LinearSolverType::Pointer mpLinearSolver = nullptr;

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief Performs the actual mapping by building and running the core mortar operator
     * @param rOriginVariable The variable on the (mapper) origin
     * @param rDestinationVariable The variable on the (mapper) destination
     * @param MappingOptions The mapping flags
     * @param Inverse If true the mapping is performed from destination to origin (InverseMap)
     */
    template<class TDataType>
    void MapInternal(
        const Variable<TDataType>& rOriginVariable,
        const Variable<TDataType>& rDestinationVariable,
        Kratos::Flags MappingOptions,
        const bool Inverse)
    {
        KRATOS_TRY

        // Resolve the data flow direction (the underlying process always maps origin -> destination)
        ModelPart& r_process_origin = Inverse ? mrModelPartDestination : mrModelPartOrigin;
        ModelPart& r_process_destination = Inverse ? mrModelPartOrigin : mrModelPartDestination;
        const std::string process_origin_variable = Inverse ? rDestinationVariable.Name() : rOriginVariable.Name();
        const std::string process_destination_variable = Inverse ? rOriginVariable.Name() : rDestinationVariable.Name();
        const bool process_origin_are_conditions = Inverse ? mMapperSettings["destination_are_conditions"].GetBool() : mMapperSettings["origin_are_conditions"].GetBool();
        const bool process_destination_are_conditions = Inverse ? mMapperSettings["origin_are_conditions"].GetBool() : mMapperSettings["destination_are_conditions"].GetBool();

        // The historical flags refer to the data flow (source/target), so they map the same way in both directions
        const bool process_origin_historical = !MappingOptions.Is(MapperFlags::FROM_NON_HISTORICAL);
        const bool process_destination_historical = !MappingOptions.Is(MapperFlags::TO_NON_HISTORICAL);

        // Build the configuration for the core mortar process
        Parameters process_parameters = mMapperSettings.Clone();
        process_parameters.RemoveValue("linear_solver_settings");
        process_parameters["origin_variable"].SetString(process_origin_variable);
        process_parameters["destination_variable"].SetString(process_destination_variable);
        process_parameters["origin_are_conditions"].SetBool(process_origin_are_conditions);
        process_parameters["destination_are_conditions"].SetBool(process_destination_are_conditions);
        process_parameters["origin_variable_historical"].SetBool(process_origin_historical);
        process_parameters["destination_variable_historical"].SetBool(process_destination_historical);

        // Build and run the core mortar operator (reuses the proven, templated implementation)
        SimpleMortarMapperProcessWrapper::Create(r_process_origin, r_process_destination, process_parameters, mpLinearSolver)->Execute();

        KRATOS_CATCH("")
    }

    /**
     * @brief Returns the default settings of the mapper (the mortar process defaults plus the linear solver settings)
     */
    Parameters GetMapperDefaultSettings() const
    {
        Parameters default_settings = SimpleMortarMapperProcessWrapper::StaticGetDefaultParameters();
        default_settings.AddValue("linear_solver_settings", Parameters(R"({})"));
        return default_settings;
    }

    ///@}

}; // Class DualMortarMapper

///@}

}  // namespace Kratos.

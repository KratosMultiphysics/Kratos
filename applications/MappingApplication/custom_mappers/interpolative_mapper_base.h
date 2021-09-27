//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher, Jordi Cotela
//
// See Master-Thesis P.Bucher
// "Development and Implementation of a Parallel
//  Framework for Non-Matching Grid Mapping"

#if !defined(KRATOS_INTERPOLATIVE_MAPPER_BASE_H_INCLUDED )
#define  KRATOS_INTERPOLATIVE_MAPPER_BASE_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/kratos_filesystem.h"
#include "input_output/vtk_output.h"
#include "utilities/variable_utils.h"

#include "mappers/mapper.h"
#include "mapping_application_variables.h"
#include "custom_searching/interface_communicator.h"
#include "custom_utilities/interface_vector_container.h"
#include "mappers/mapper_flags.h"
#include "custom_utilities/mapper_local_system.h"
#include "custom_utilities/mapping_matrix_utilities.h"
#include "custom_utilities/mapper_utilities.h"

namespace Kratos
{
///@name Kratos Classes
///@{

template<class TSparseSpace, class TDenseSpace, class TMapperBackend>
class KRATOS_API(MAPPING_APPLICATION) InterpolativeMapperBase : public Mapper<TSparseSpace, TDenseSpace>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of InterpolativeMapperBase
    KRATOS_CLASS_POINTER_DEFINITION(InterpolativeMapperBase);

    typedef Mapper<TSparseSpace, TDenseSpace> BaseType;

    typedef typename TMapperBackend::InterfaceCommunicatorType InterfaceCommunicatorType;
    typedef Kratos::unique_ptr<InterfaceCommunicator> InterfaceCommunicatorPointerType;
    typedef typename InterfaceCommunicator::MapperInterfaceInfoUniquePointerType MapperInterfaceInfoUniquePointerType;

    typedef Kratos::unique_ptr<MapperLocalSystem> MapperLocalSystemPointer;
    typedef std::vector<MapperLocalSystemPointer> MapperLocalSystemPointerVector;

    typedef InterfaceVectorContainer<TSparseSpace, TDenseSpace> InterfaceVectorContainerType;
    typedef Kratos::unique_ptr<InterfaceVectorContainerType> InterfaceVectorContainerPointerType;

    typedef std::size_t IndexType;

    typedef typename BaseType::MapperUniquePointerType MapperUniquePointerType;
    typedef typename BaseType::TMappingMatrixType TMappingMatrixType;
    typedef Kratos::unique_ptr<TMappingMatrixType> TMappingMatrixUniquePointerType;

    typedef Variable<double> ComponentVariableType;

    ///@}
    ///@name Life Cycle
    ///@{

    // Default constructor, needed for registration
    InterpolativeMapperBase(ModelPart& rModelPartOrigin,
                         ModelPart& rModelPartDestination)
                         : mrModelPartOrigin(rModelPartOrigin),
                           mrModelPartDestination(rModelPartDestination) {}

    InterpolativeMapperBase(ModelPart& rModelPartOrigin,
                         ModelPart& rModelPartDestination,
                         Parameters JsonParameters)
                         : mrModelPartOrigin(rModelPartOrigin),
                           mrModelPartDestination(rModelPartDestination),
                           mMapperSettings(JsonParameters)
    {
        mpInterfaceVectorContainerOrigin = Kratos::make_unique<InterfaceVectorContainerType>(rModelPartOrigin);
        mpInterfaceVectorContainerDestination = Kratos::make_unique<InterfaceVectorContainerType>(rModelPartDestination);
    }

    /// Destructor.
    ~InterpolativeMapperBase() override = default;

    ///@}
    ///@name Operations
    ///@{

    void UpdateInterface(
        Kratos::Flags MappingOptions,
        double SearchRadius) override
    {
        KRATOS_TRY;

        KRATOS_WARNING_IF("Mapper", mMapperSettings["use_initial_configuration"].GetBool()) << "Updating the interface while using the initial configuration for mapping!" << std::endl;

        // Set the Flags according to the type of remeshing
        if (MappingOptions.Is(MapperFlags::REMESHED)) {
            InitializeInterface(MappingOptions);
        }
        else {
            BuildMappingMatrix(MappingOptions);
        }

        if (mpInverseMapper) {
            mpInverseMapper->UpdateInterface(MappingOptions,
                                             SearchRadius);
        }

        KRATOS_CATCH("");
    }

    void Map(
        const Variable<double>& rOriginVariable,
        const Variable<double>& rDestinationVariable,
        Kratos::Flags MappingOptions) override
    {
        KRATOS_TRY;

        if (MappingOptions.Is(MapperFlags::USE_TRANSPOSE)) {
            GetInverseMapper().MapInternalTranspose(rDestinationVariable, rOriginVariable, MappingOptions);
        } else {
            MapInternal(rOriginVariable, rDestinationVariable, MappingOptions);
        }

        KRATOS_CATCH("");
    }

    void Map(
        const Variable< array_1d<double, 3> >& rOriginVariable,
        const Variable< array_1d<double, 3> >& rDestinationVariable,
        Kratos::Flags MappingOptions) override
    {
        KRATOS_TRY;

        if (MappingOptions.Is(MapperFlags::USE_TRANSPOSE)) {
            GetInverseMapper().MapInternalTranspose(rDestinationVariable, rOriginVariable, MappingOptions);
        } else {
            MapInternal(rOriginVariable, rDestinationVariable, MappingOptions);
        }

        KRATOS_CATCH("");
    }

    void InverseMap(
        const Variable<double>& rOriginVariable,
        const Variable<double>& rDestinationVariable,
        Kratos::Flags MappingOptions) override
    {
        KRATOS_TRY;

        if (MappingOptions.Is(MapperFlags::USE_TRANSPOSE)) {
            MapInternalTranspose(rOriginVariable, rDestinationVariable, MappingOptions);
        } else {
            GetInverseMapper().Map(rDestinationVariable, rOriginVariable, MappingOptions);
        }

        KRATOS_CATCH("");
    }

    void InverseMap(
        const Variable< array_1d<double, 3> >& rOriginVariable,
        const Variable< array_1d<double, 3> >& rDestinationVariable,
        Kratos::Flags MappingOptions) override
    {
        KRATOS_TRY;

        if (MappingOptions.Is(MapperFlags::USE_TRANSPOSE)) {
            MapInternalTranspose(rOriginVariable, rDestinationVariable, MappingOptions);
        } else {
            GetInverseMapper().Map(rDestinationVariable, rOriginVariable, MappingOptions);
        }

        KRATOS_CATCH("");
    }

    ///@}
    ///@name Access
    ///@{

    TMappingMatrixType& GetMappingMatrix() override
    {
        return *mpMappingMatrix;
    }

    ///@}
    ///@name Inquiry
    ///@{

    int AreMeshesConforming() const override
    {
        KRATOS_WARNING_ONCE("Mapper") << "Developer-warning: \"AreMeshesConforming\" is deprecated and will be removed in the future" << std::endl;
        return mpIntefaceCommunicator->AreMeshesConforming();
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "InterpolativeMapperBase";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "InterpolativeMapperBase";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        BaseType::PrintData(rOStream);
    }

protected:

   /**
    * @brief Initializing the Mapper
    * This has to be called in the constructor of the
    * derived classes, since it involves calls to
    * pure virtual functions
    */
    void Initialize()
    {
        KRATOS_TRY;

        InitializeInterfaceCommunicator();
        InitializeInterface();

        KRATOS_CATCH("");
    }

    void ValidateInput()
    {
        MapperUtilities::CheckInterfaceModelParts(0);

        Parameters mapper_default_settings(GetMapperDefaultSettings());
        mMapperSettings.ValidateAndAssignDefaults(mapper_default_settings);

        if (mMapperSettings["search_radius"].GetDouble() < 0.0) {
            const double search_radius = MapperUtilities::ComputeSearchRadius(
                                            mrModelPartOrigin,
                                            mrModelPartDestination,
                                            mMapperSettings["echo_level"].GetInt());
            mMapperSettings["search_radius"].SetDouble(search_radius);
        }
    }

private:
    ///@name Member Variables
    ///@{

    ModelPart& mrModelPartOrigin;
    ModelPart& mrModelPartDestination;

    Parameters mMapperSettings;

    MapperUniquePointerType mpInverseMapper = nullptr;

    TMappingMatrixUniquePointerType mpMappingMatrix;

    MapperLocalSystemPointerVector mMapperLocalSystems;

    InterfaceCommunicatorPointerType mpIntefaceCommunicator;
    InterfaceVectorContainerPointerType mpInterfaceVectorContainerOrigin;
    InterfaceVectorContainerPointerType mpInterfaceVectorContainerDestination;

    ///@}
    ///@name Private Operations
    ///@{

    void InitializeInterfaceCommunicator()
    {
        KRATOS_TRY;

        mpIntefaceCommunicator = Kratos::make_unique<InterfaceCommunicatorType>(
            mrModelPartOrigin,
            mMapperLocalSystems,
            mMapperSettings);

        KRATOS_CATCH("");
    }

    void InitializeInterface(Kratos::Flags MappingOptions = Kratos::Flags())
    {
        KRATOS_TRY;

        CreateMapperLocalSystems(mrModelPartDestination.GetCommunicator(),
                                mMapperLocalSystems);

        BuildMappingMatrix(MappingOptions);

        KRATOS_CATCH("");
    }

    void BuildMappingMatrix(Kratos::Flags MappingOptions = Kratos::Flags())
    {
        KRATOS_TRY;

        const bool use_initial_configuration = mMapperSettings["use_initial_configuration"].GetBool();

        if (use_initial_configuration) {
            MapperUtilities::SaveCurrentConfiguration(mrModelPartOrigin);
            MapperUtilities::SaveCurrentConfiguration(mrModelPartDestination);

            VariableUtils().UpdateCurrentToInitialConfiguration(mrModelPartOrigin.Nodes());
            VariableUtils().UpdateCurrentToInitialConfiguration(mrModelPartDestination.Nodes());
        }

        AssignInterfaceEquationIds(); // Has to be done ever time in case of overlapping interfaces!

        KRATOS_ERROR_IF_NOT(mpIntefaceCommunicator) << "mpIntefaceCommunicator is a nullptr!" << std::endl;

        const MapperInterfaceInfoUniquePointerType p_ref_interface_info = GetMapperInterfaceInfo();

        mpIntefaceCommunicator->ExchangeInterfaceData(mrModelPartDestination.GetCommunicator(),
                                                    MappingOptions,
                                                    p_ref_interface_info);

        const int echo_level = mMapperSettings["echo_level"].GetInt();

        MappingMatrixUtilities<TSparseSpace, TDenseSpace>::BuildMappingMatrix(
            mpMappingMatrix,
            mpInterfaceVectorContainerOrigin->pGetVector(),
            mpInterfaceVectorContainerDestination->pGetVector(),
            mpInterfaceVectorContainerOrigin->GetModelPart(),
            mpInterfaceVectorContainerDestination->GetModelPart(),
            mMapperLocalSystems,
            echo_level);

        if (echo_level > 0) {
            PrintPairingInfo(echo_level);
        }

        if (use_initial_configuration) {
            MapperUtilities::RestoreCurrentConfiguration(mrModelPartOrigin);
            MapperUtilities::RestoreCurrentConfiguration(mrModelPartDestination);
        }

        KRATOS_CATCH("");
    }

    void AssignInterfaceEquationIds()
    {
        MapperUtilities::AssignInterfaceEquationIds(mrModelPartOrigin.GetCommunicator());
        MapperUtilities::AssignInterfaceEquationIds(mrModelPartDestination.GetCommunicator());
    }

    void MapInternal(const Variable<double>& rOriginVariable,
                     const Variable<double>& rDestinationVariable,
                     Kratos::Flags MappingOptions)
    {
        KRATOS_TRY;

        mpInterfaceVectorContainerOrigin->UpdateSystemVectorFromModelPart(rOriginVariable, MappingOptions);

        TSparseSpace::Mult(
            *mpMappingMatrix,
            mpInterfaceVectorContainerOrigin->GetVector(),
            mpInterfaceVectorContainerDestination->GetVector()); // rQd = rMdo * rQo

        mpInterfaceVectorContainerDestination->UpdateModelPartFromSystemVector(rDestinationVariable, MappingOptions);

        KRATOS_CATCH("");
    }

    void MapInternalTranspose(const Variable<double>& rOriginVariable,
                              const Variable<double>& rDestinationVariable,
                              Kratos::Flags MappingOptions)
    {
        KRATOS_TRY;

        mpInterfaceVectorContainerDestination->UpdateSystemVectorFromModelPart(rDestinationVariable, MappingOptions);

        TSparseSpace::TransposeMult(
            *mpMappingMatrix,
            mpInterfaceVectorContainerDestination->GetVector(),
            mpInterfaceVectorContainerOrigin->GetVector()); // rQo = rMdo^T * rQd

        mpInterfaceVectorContainerOrigin->UpdateModelPartFromSystemVector(rOriginVariable, MappingOptions);

        KRATOS_CATCH("");
    }

    void MapInternal(const Variable<array_1d<double, 3>>& rOriginVariable,
                     const Variable<array_1d<double, 3>>& rDestinationVariable,
                     Kratos::Flags MappingOptions)
    {
        KRATOS_TRY;

        for (const auto var_ext : {"_X", "_Y", "_Z"}) {
            const auto& var_origin = KratosComponents<ComponentVariableType>::Get(rOriginVariable.Name() + var_ext);
            const auto& var_destination = KratosComponents<ComponentVariableType>::Get(rDestinationVariable.Name() + var_ext);

            MapInternal(var_origin, var_destination, MappingOptions);
        }

        KRATOS_CATCH("");
    }

    void MapInternalTranspose(const Variable<array_1d<double, 3>>& rOriginVariable,
                              const Variable<array_1d<double, 3>>& rDestinationVariable,
                              Kratos::Flags MappingOptions)
    {
        KRATOS_TRY;

        for (const auto var_ext : {"_X", "_Y", "_Z"}) {
            const auto& var_origin = KratosComponents<ComponentVariableType>::Get(rOriginVariable.Name() + var_ext);
            const auto& var_destination = KratosComponents<ComponentVariableType>::Get(rDestinationVariable.Name() + var_ext);

            MapInternalTranspose(var_origin, var_destination, MappingOptions);
        }

        KRATOS_CATCH("");
    }

    void PrintPairingInfo(const int EchoLevel)
    {
        std::stringstream warning_msg;

        if (EchoLevel > 1) {
            // Initialize the values for printing later
            VariableUtils().SetNonHistoricalVariable(PAIRING_STATUS, 1, mrModelPartDestination.Nodes());
        }

        for (const auto& rp_local_sys : mMapperLocalSystems) {
            const auto pairing_status = rp_local_sys->GetPairingStatus();

            if (pairing_status != MapperLocalSystem::PairingStatus::InterfaceInfoFound) {
                warning_msg << rp_local_sys->PairingInfo(EchoLevel);

                if (pairing_status == MapperLocalSystem::PairingStatus::Approximation)
                    warning_msg << " is using an approximation";
                else if (pairing_status == MapperLocalSystem::PairingStatus::NoInterfaceInfo)
                    warning_msg << " has not found a neighbor";

                KRATOS_WARNING_ALL_RANKS("Mapper") << warning_msg.str() << std::endl; // TODO use data-comm of the destination-MP

                // reset the stringstream
                warning_msg.str( std::string() );
                warning_msg.clear();
            }
        }

        if (mMapperSettings["print_pairing_status_to_file"].GetBool()) {
            // print a debug ModelPart to check the pairing

            const std::string pairing_status_file_path = mMapperSettings["pairing_status_file_path"].GetString();

            filesystem::create_directories(pairing_status_file_path);

            const std::string file_name = FilesystemExtensions::JoinPaths({
                pairing_status_file_path,
                std::string(Info() + "_PairingStatus_O_" + mrModelPartOrigin.FullName() + "_D_" + mrModelPartDestination.FullName())
            });

            KRATOS_INFO("Mapper") << "Printing file with PAIRING_STATUS: " << file_name << ".vtk" << std::endl;

            Parameters vtk_params( R"({
                "file_format"                        : "binary",
                "save_output_files_in_folder"        : false,
                "nodal_data_value_variables"         : ["PAIRING_STATUS"]
            })");

            VtkOutput(mrModelPartDestination, vtk_params).PrintOutput(file_name);
        }
    }

    // functions for customizing the behavior of this Mapper
    virtual void CreateMapperLocalSystems(
        const Communicator& rModelPartCommunicator,
        std::vector<Kratos::unique_ptr<MapperLocalSystem>>& rLocalSystems) = 0;

    virtual MapperInterfaceInfoUniquePointerType GetMapperInterfaceInfo() const = 0;

    virtual Parameters GetMapperDefaultSettings() const = 0;

    ///@}
    ///@name Private  Access
    ///@{

    InterpolativeMapperBase& GetInverseMapper()
    {
        KRATOS_TRY;

        if (!mpInverseMapper) {
            InitializeInverseMapper();
        }
        return *(static_cast<InterpolativeMapperBase*>(mpInverseMapper.get()));

        KRATOS_CATCH("");
    }

    void InitializeInverseMapper()
    {
        KRATOS_TRY;

        mpInverseMapper = this->Clone(mrModelPartDestination,
                                      mrModelPartOrigin,
                                      mMapperSettings);

        KRATOS_CATCH("");
    }

}; // Class InterpolativeMapperBase

}  // namespace Kratos.

#endif // KRATOS_INTERPOLATIVE_MAPPER_BASE_H_INCLUDED  defined

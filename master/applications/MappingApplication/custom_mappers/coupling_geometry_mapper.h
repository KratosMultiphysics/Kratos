//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher
//                   Tobias Teschemachen
//

#pragma once

// System includes

// External includes

// Project includes
#include "mappers/mapper.h"
#include "custom_utilities/interface_vector_container.h"
#include "custom_utilities/mapper_local_system.h"

#include "custom_utilities/mapping_intersection_utilities.h"
#include "custom_modelers/mapping_geometries_modeler.h"
#include "modeler/modeler_factory.h"

#include "linear_solvers/linear_solver.h"

namespace Kratos
{
///@name Kratos Classes
///@{

class CouplingGeometryLocalSystem : public MapperLocalSystem
{
public:

    explicit CouplingGeometryLocalSystem(GeometryPointerType pGeom,
                                         const bool IsProjection,
                                         const bool IsDualMortar,
                                         const bool IsDestinationIsSlave
                                         )
        : mpGeom(pGeom),
          mIsProjection(IsProjection),
          mIsDualMortar(IsDualMortar),
          mIsDestinationIsSlave(IsDestinationIsSlave)
        {}

    void CalculateAll(MatrixType& rLocalMappingMatrix,
                      EquationIdVectorType& rOriginIds,
                      EquationIdVectorType& rDestinationIds,
                      MapperLocalSystem::PairingStatus& rPairingStatus) const override;

    CoordinatesArrayType& Coordinates() const override
    {
        KRATOS_DEBUG_ERROR_IF_NOT(mpGeom) << "Members are not initialized!" << std::endl;
        // return mpGeom->Center(); // check why not compiling...
        KRATOS_ERROR << "not implemented, needs checking" << std::endl;
    }

    MapperLocalSystemUniquePointer Create(GeometryPointerType pGeometry) const override
    {
        return Kratos::make_unique<CouplingGeometryLocalSystem>(pGeometry, mIsProjection, mIsDualMortar, mIsDestinationIsSlave);
    }

    /// Turn back information as a string.
    void PairingInfo(std::ostream& rOStream, const int EchoLevel) const override {KRATOS_ERROR << "Not implemented!"<<std::endl;}

private:
    GeometryPointerType mpGeom;
    bool mIsProjection; // Set to true is we are projecting the master onto the slave.
                        // Set to false if we are projecting the slave onto the slave.
    bool mIsDualMortar = false;
    bool mIsDestinationIsSlave = true;

};

// CouplingGeometryMapper
//
// The mapper always forward maps from the master to the slave.
// Normally:
//      master  =   interface origin
//      slave   =   interface destination
//
// However, this can be reversed by setting 'destination_is_slave' = false.
// This yields:
//      master  =   interface destination
//      slave   =   interface origin

template<class TSparseSpace, class TDenseSpace>
class CouplingGeometryMapper : public Mapper<TSparseSpace, TDenseSpace>
{
public:

    ///@name Type Definitions
    ///@{

    ///@}
    ///@name Pointer Definitions
    /// Pointer definition of CouplingGeometryMapper
    KRATOS_CLASS_POINTER_DEFINITION(CouplingGeometryMapper);

    typedef Mapper<TSparseSpace, TDenseSpace> BaseType;

    typedef Kratos::unique_ptr<MapperLocalSystem> MapperLocalSystemPointer;
    typedef std::vector<MapperLocalSystemPointer> MapperLocalSystemPointerVector;

    typedef InterfaceVectorContainer<TSparseSpace, TDenseSpace> InterfaceVectorContainerType;
    typedef Kratos::unique_ptr<InterfaceVectorContainerType> InterfaceVectorContainerPointerType;

    typedef std::size_t IndexType;

    typedef typename BaseType::MapperUniquePointerType MapperUniquePointerType;
    typedef typename BaseType::TMappingMatrixType MappingMatrixType;
    typedef Kratos::unique_ptr<MappingMatrixType> MappingMatrixUniquePointerType;

    typedef LinearSolver<TSparseSpace, TDenseSpace> LinearSolverType;
    typedef Kratos::shared_ptr<LinearSolverType> LinearSolverSharedPointerType;

    typedef typename TSparseSpace::VectorType TSystemVectorType;
    typedef Kratos::unique_ptr<TSystemVectorType> TSystemVectorUniquePointerType;

    ///@}
    ///@name Life Cycle
    ///@{

    // Default constructor, needed for registration
    CouplingGeometryMapper(ModelPart& rModelPartOrigin,
                         ModelPart& rModelPartDestination)
                        : mrModelPartOrigin(rModelPartOrigin),
                          mrModelPartDestination(rModelPartDestination){}


    CouplingGeometryMapper(ModelPart& rModelPartOrigin,
                         ModelPart& rModelPartDestination,
                         Parameters JsonParameters);

    /// Destructor.
    ~CouplingGeometryMapper() override = default;

    ///@}
    ///@name Operations
    ///@{

    void UpdateInterface(
        Kratos::Flags MappingOptions,
        double SearchRadius) override
    {
        mpModeler->PrepareGeometryModel();

        AssignInterfaceEquationIds();

        KRATOS_ERROR << "Not implemented!" << std::endl;
    }

    void Map(
        const Variable<double>& rOriginVariable,
        const Variable<double>& rDestinationVariable,
        Kratos::Flags MappingOptions) override
    {
        if (MappingOptions.Is(MapperFlags::USE_TRANSPOSE)) {
            MappingOptions.Reset(MapperFlags::USE_TRANSPOSE);
            MappingOptions.Set(MapperFlags::INTERNAL_USE_TRANSPOSE, true);
            GetInverseMapper()->Map(rDestinationVariable, rOriginVariable, MappingOptions);
        }
        else if (MappingOptions.Is(MapperFlags::INTERNAL_USE_TRANSPOSE)) {
            MapInternalTranspose(rOriginVariable, rDestinationVariable, MappingOptions);
        }
        else {
            MapInternal(rOriginVariable, rDestinationVariable, MappingOptions);
        }
    }

    void Map(
        const Variable< array_1d<double, 3> >& rOriginVariable,
        const Variable< array_1d<double, 3> >& rDestinationVariable,
        Kratos::Flags MappingOptions) override
    {
        if (MappingOptions.Is(MapperFlags::USE_TRANSPOSE)) {
            MappingOptions.Reset(MapperFlags::USE_TRANSPOSE);
            MappingOptions.Set(MapperFlags::INTERNAL_USE_TRANSPOSE, true);
            GetInverseMapper()->Map(rDestinationVariable, rOriginVariable, MappingOptions);
        }
        else if (MappingOptions.Is(MapperFlags::INTERNAL_USE_TRANSPOSE)) {
            MapInternalTranspose(rOriginVariable, rDestinationVariable, MappingOptions);
        }
        else {
            MapInternal(rOriginVariable, rDestinationVariable, MappingOptions);
        }
    }

    void InverseMap(
        const Variable<double>& rOriginVariable,
        const Variable<double>& rDestinationVariable,
        Kratos::Flags MappingOptions) override
    {
        if (MappingOptions.Is(MapperFlags::USE_TRANSPOSE)) {
            MapInternalTranspose(rOriginVariable, rDestinationVariable, MappingOptions);
        }
        else {
            GetInverseMapper()->Map(rDestinationVariable, rOriginVariable, MappingOptions);
        }
    }

    void InverseMap(
        const Variable< array_1d<double, 3> >& rOriginVariable,
        const Variable< array_1d<double, 3> >& rDestinationVariable,
        Kratos::Flags MappingOptions) override
    {
        if (MappingOptions.Is(MapperFlags::USE_TRANSPOSE)) {
            MapInternalTranspose(rOriginVariable, rDestinationVariable, MappingOptions);
        }
        else {
            GetInverseMapper()->Map(rDestinationVariable, rOriginVariable, MappingOptions);
        }
    }

    ///@}
    ///@name Access
    ///@{

    MappingMatrixType& GetMappingMatrix() override
    {
        if (mMapperSettings["precompute_mapping_matrix"].GetBool() || mMapperSettings["dual_mortar"].GetBool()) return *(mpMappingMatrix.get());
        else KRATOS_ERROR << "'precompute_mapping_matrix' or 'dual_mortar' must be 'true' in your parameters to retrieve the computed mapping matrix!" << std::endl;
    }

    MapperUniquePointerType Clone(ModelPart& rModelPartOrigin,
                                  ModelPart& rModelPartDestination,
                                  Parameters JsonParameters) const override
    {
        return Kratos::make_unique<CouplingGeometryMapper<TSparseSpace, TDenseSpace>>(
            rModelPartOrigin,
            rModelPartDestination,
            JsonParameters);
    }

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "CouplingGeometryMapper";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "CouplingGeometryMapper";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        BaseType::PrintData(rOStream);
    }

    // Get values
    // Always returns the true origin/destination regardless of 'destination_is_slave'
    ModelPart& GetInterfaceModelPartOrigin() override
    {

        return mpCouplingMP->GetSubModelPart("interface_origin");
    }

    ModelPart& GetInterfaceModelPartDestination() override
    {
        return mpCouplingMP->GetSubModelPart("interface_destination");
    }

private:

    ///@name Private Operations
    ///@{
    typename Modeler::Pointer mpModeler = nullptr;

    ModelPart& mrModelPartOrigin;
    ModelPart& mrModelPartDestination;
    ModelPart* mpCouplingMP = nullptr;
    ModelPart* mpCouplingInterfaceMaster = nullptr;
    ModelPart* mpCouplingInterfaceSlave = nullptr;

    Parameters mMapperSettings;

    MapperUniquePointerType mpInverseMapper = nullptr;

    MappingMatrixUniquePointerType mpMappingMatrix;
    MappingMatrixUniquePointerType mpMappingMatrixProjector;
    MappingMatrixUniquePointerType mpMappingMatrixSlave;

    TSystemVectorUniquePointerType mpTempVector;

    MapperLocalSystemPointerVector mMapperLocalSystemsProjector;
    MapperLocalSystemPointerVector mMapperLocalSystemsSlave;

    InterfaceVectorContainerPointerType mpInterfaceVectorContainerMaster;
    InterfaceVectorContainerPointerType mpInterfaceVectorContainerSlave;

    LinearSolverSharedPointerType mpLinearSolver = nullptr;


    void InitializeInterface(Kratos::Flags MappingOptions = Kratos::Flags());

    void AssignInterfaceEquationIds()
    {
        MapperUtilities::AssignInterfaceEquationIds(mpCouplingInterfaceSlave->GetCommunicator());
        MapperUtilities::AssignInterfaceEquationIds(mpCouplingInterfaceMaster->GetCommunicator());
    }

    void MapInternal(const Variable<double>& rOriginVariable,
                     const Variable<double>& rDestinationVariable,
                     Kratos::Flags MappingOptions);

    void MapInternalTranspose(const Variable<double>& rOriginVariable,
                              const Variable<double>& rDestinationVariable,
                              Kratos::Flags MappingOptions);

    void MapInternal(const Variable<array_1d<double, 3>>& rOriginVariable,
                     const Variable<array_1d<double, 3>>& rDestinationVariable,
                     Kratos::Flags MappingOptions);

    void MapInternalTranspose(const Variable<array_1d<double, 3>>& rOriginVariable,
                              const Variable<array_1d<double, 3>>& rDestinationVariable,
                              Kratos::Flags MappingOptions);

    void EnforceConsistencyWithScaling(
        const MappingMatrixType& rInterfaceMatrixSlave,
        MappingMatrixType& rInterfaceMatrixProjected,
        const double scalingLimit = 1.1);

    void CreateLinearSolver();

    void CalculateMappingMatrixWithSolver(MappingMatrixType& rConsistentInterfaceMatrix, MappingMatrixType& rProjectedInterfaceMatrix);

    Parameters GetMapperDefaultSettings() const
    {
        return Parameters(R"({
            "echo_level"                    : 0,
            "dual_mortar"                   : false,
            "precompute_mapping_matrix"     : false,
            "modeler_name"                  : "UNSPECIFIED",
            "modeler_parameters"            : {},
            "consistency_scaling"           : true,
            "row_sum_tolerance"             : 1e-12,
            "destination_is_slave"          : true,
            "linear_solver_settings"        : {}
        })");
    }

    ///@}
    ///@name Private  Access
    ///@{

    MapperUniquePointerType& GetInverseMapper()
    {
        if (!mpInverseMapper) {
            InitializeInverseMapper();
        }
        return mpInverseMapper;
    }

    void InitializeInverseMapper()
    {
        KRATOS_ERROR << "Inverse Mapping is not supported yet!" << std::endl;
        mpInverseMapper = this->Clone(mrModelPartDestination,
                                      mrModelPartOrigin,
                                      mMapperSettings);
    }

    ///@}

}; // Class CouplingGeometryMapper

///@} addtogroup block
}  // namespace Kratos.
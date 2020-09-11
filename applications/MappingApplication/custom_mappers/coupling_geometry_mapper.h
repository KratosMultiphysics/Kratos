//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher
//                   Tobias Teschemachen
//

#if !defined(KRATOS_COUPLING_GEOMETRY_MAPPER_H_INCLUDED)
#define  KRATOS_COUPLING_GEOMETRY_MAPPER_H_INCLUDED

// System includes

// External includes

// Project includes
#include "mapper.h"
#include "custom_utilities/interface_vector_container.h"
#include "custom_utilities/mapper_flags.h"
#include "custom_utilities/mapper_local_system.h"

#include "custom_utilities/mapping_intersection_utilities.h"
#include "custom_modelers/mapping_geometries_modeler.h"
#include "modeler/modeler_factory.h"

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
                                         const bool IsConsistent = false
                                         )
        : mpGeom(pGeom),
          mIsProjection(IsProjection),
          mIsDualMortar(IsDualMortar),
          mIsConsistentMortar(IsConsistent)
    {
        KRATOS_ERROR_IF(mIsConsistentMortar && mIsProjection)
            << "CouplingGeometryLocalSystem incorrectly created with mIsConsistentMortar == true and mIsProjection == true."
            << "\nThese two options are incompatible.\n";
    }

    void CalculateAll(MatrixType& rLocalMappingMatrix,
                      EquationIdVectorType& rOriginIds,
                      EquationIdVectorType& rDestinationIds,
                      MapperLocalSystem::PairingStatus& rPairingStatus) const override;

    CoordinatesArrayType& Coordinates() const override
    {
        KRATOS_DEBUG_ERROR_IF_NOT(mpGeom) << "Members are not intitialized!" << std::endl;
        // return mpGeom->Center(); // check why not compiling...
        KRATOS_ERROR << "not implemented, needs checking" << std::endl;
    }

    /// Turn back information as a string.
    std::string PairingInfo(const int EchoLevel) const override;

private:
    GeometryPointerType mpGeom;
    bool mIsProjection; // Set to true is we are projecting the master onto the slave.
                        // Set to false if we are projecting the slave onto the slave.
    bool mIsConsistentMortar = false; // by default we use conservative mapping
    bool mIsDualMortar = false;

};

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
    typedef typename BaseType::TMappingMatrixType TMappingMatrixType;
    typedef Kratos::unique_ptr<TMappingMatrixType> TMappingMatrixUniquePointerType;

    typedef Matrix DenseMappingMatrixType;
    typedef Kratos::unique_ptr<DenseMappingMatrixType> DenseMappingMatrixUniquePointerType;

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
                         Parameters JsonParameters)
                        : mrModelPartOrigin(rModelPartOrigin),
                          mrModelPartDestination(rModelPartDestination),
                          mMapperSettings(JsonParameters)
    {
        mpModeler = (ModelerFactory::Create(
            mMapperSettings["modeler_name"].GetString(),
            rModelPartOrigin.GetModel(),
            mMapperSettings["modeler_parameters"]));

        // adds destination model part
        mpModeler->GenerateNodes(rModelPartDestination);

        mpModeler->SetupGeometryModel();
        mpModeler->PrepareGeometryModel();

        // here use whatever ModelPart(s) was created by the Modeler
        mpCouplingMP = &(rModelPartOrigin.GetModel().GetModelPart("coupling"));
        mpCouplingInterfaceOrigin = mpCouplingMP->pGetSubModelPart("interface_origin");
        mpCouplingInterfaceDestination = mpCouplingMP->pGetSubModelPart("interface_destination");

        mpInterfaceVectorContainerOrigin = Kratos::make_unique<InterfaceVectorContainerType>(*mpCouplingInterfaceOrigin);
        mpInterfaceVectorContainerDestination = Kratos::make_unique<InterfaceVectorContainerType>(*mpCouplingInterfaceDestination);

        mpCouplingMP->GetMesh().SetValue(IS_DUAL_MORTAR, mMapperSettings["dual_mortar"].GetBool());

        mpCouplingMP->GetMesh().SetValue(IS_CONSISTENT_MORTAR, false);
        if (mMapperSettings.Has("consistent_mortar"))
        {
            mpCouplingMP->GetMesh().SetValue(IS_CONSISTENT_MORTAR, mMapperSettings["consistent_mortar"].GetBool());
        }

        this->InitializeInterface();
    }

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

private:

    ///@name Private Operations
    ///@{
    typename Modeler::Pointer mpModeler = nullptr;

    ModelPart& mrModelPartOrigin;
    ModelPart& mrModelPartDestination;
    ModelPart* mpCouplingMP = nullptr;
    ModelPart* mpCouplingInterfaceOrigin = nullptr;
    ModelPart* mpCouplingInterfaceDestination = nullptr;

    Parameters mMapperSettings;

    MapperUniquePointerType mpInverseMapper = nullptr;

    DenseMappingMatrixUniquePointerType mpMappingMatrix;
    DenseMappingMatrixUniquePointerType mpMappingMatrix_consistent_force;

    MapperLocalSystemPointerVector mMapperLocalSystems;

    InterfaceVectorContainerPointerType mpInterfaceVectorContainerOrigin;
    InterfaceVectorContainerPointerType mpInterfaceVectorContainerDestination;


    void InitializeInterface(Kratos::Flags MappingOptions = Kratos::Flags());

    void AssignInterfaceEquationIds()
    {
        MapperUtilities::AssignInterfaceEquationIds(mpCouplingInterfaceDestination->GetCommunicator());
        MapperUtilities::AssignInterfaceEquationIds(mpCouplingInterfaceOrigin->GetCommunicator());
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

    void CreateMapperLocalSystems(
        const Communicator& rModelPartCommunicator,
        std::vector<Kratos::unique_ptr<MapperLocalSystem>>& rLocalSystems)
    {
        MapperUtilities::CreateMapperLocalSystemsFromGeometries<CouplingGeometryLocalSystem>(
            rModelPartCommunicator,
            rLocalSystems);
    }

    void EnforceConsistencyWithScaling(const Matrix& rInterfaceMatrixSlave,
        Matrix& rInterfaceMatrixProjected,
        const double scalingLimit = 1.1);

    void CheckMappingMatrixConsistency(Matrix& rMappingMatrix)
    {
        for (size_t row = 0; row < rMappingMatrix.size1(); ++row) {
            double row_sum = 0.0;
            for (size_t col = 0; col < rMappingMatrix.size2(); ++col) row_sum += rMappingMatrix(row, col);
            if (std::abs(row_sum - 1.0) > 1e-12) {
                KRATOS_WATCH(rMappingMatrix)
                KRATOS_WATCH(row_sum)
                KRATOS_ERROR << "mapping matrix is not consistent\n";
            }
        }
    }

    DenseMappingMatrixUniquePointerType CreateInterfaceMappingMatrix(
        Matrix& rConsistentComponent, Matrix& rProjectedComponent)
    {
        // get total interface mapping matrix
        Matrix inv_interface_matrix_slave(rConsistentComponent.size1(), rConsistentComponent.size2(), 0.0);
        if (mMapperSettings["dual_mortar"].GetBool())
        {
            for (size_t i = 0; i < rConsistentComponent.size1(); ++i) {
                if (rConsistentComponent(i, i) > std::numeric_limits<double>::epsilon()) {
                    inv_interface_matrix_slave(i, i) = 1.0 / rConsistentComponent(i, i);
                }
            }
        }
        else
        {
            double aux_det_slave = 0;
            MathUtils<double>::InvertMatrix(rConsistentComponent, inv_interface_matrix_slave, aux_det_slave);
        }

        return Kratos::make_unique<DenseMappingMatrixType>(prod(inv_interface_matrix_slave, rProjectedComponent));
    }

    Parameters GetMapperDefaultSettings() const
    {
        // @tobiasteschemachen
        return Parameters( R"({
            "echo_level" : 0
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

#endif // KRATOS_COUPLING_GEOMETRY_MAPPER_H_INCLUDED defined
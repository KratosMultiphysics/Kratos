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


namespace Kratos
{
///@name Kratos Classes
///@{

class CouplingGeometryLocalSystem : public MapperLocalSystem
{
public:

    explicit CouplingGeometryLocalSystem(GeometryPointerType pGeom) : mpGeom(pGeom) {}

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

};

template<class TSparseSpace, class TDenseSpace>
class CouplingGeomteryMapper : public Mapper<TSparseSpace, TDenseSpace>
{
public:

    ///@name Type Definitions
    ///@{

    ///@}
    ///@name Pointer Definitions
    /// Pointer definition of CouplingGeomteryMapper
    KRATOS_CLASS_POINTER_DEFINITION(CouplingGeomteryMapper);

    typedef Mapper<TSparseSpace, TDenseSpace> BaseType;

    typedef Kratos::unique_ptr<MapperLocalSystem> MapperLocalSystemPointer;
    typedef std::vector<MapperLocalSystemPointer> MapperLocalSystemPointerVector;

    typedef InterfaceVectorContainer<TSparseSpace, TDenseSpace> InterfaceVectorContainerType;
    typedef Kratos::unique_ptr<InterfaceVectorContainerType> InterfaceVectorContainerPointerType;

    typedef std::size_t IndexType;

    typedef typename BaseType::MapperUniquePointerType MapperUniquePointerType;
    typedef typename BaseType::TMappingMatrixType TMappingMatrixType;
    typedef Kratos::unique_ptr<TMappingMatrixType> TMappingMatrixUniquePointerType;

    typedef VariableComponent< VectorComponentAdaptor<array_1d<double, 3> > > ComponentVariableType;

    ///@}
    ///@name Life Cycle
    ///@{

    // Default constructor, needed for registration
    CouplingGeomteryMapper(ModelPart& rModelPartOrigin,
                         ModelPart& rModelPartDestination)
                         : mrModelPartOrigin(rModelPartOrigin),
                           mrModelPartDestination(rModelPartDestination) {}

    CouplingGeomteryMapper(ModelPart& rModelPartOrigin,
                         ModelPart& rModelPartDestination,
                         Parameters JsonParameters)
                         : mrModelPartOrigin(rModelPartOrigin),
                           mrModelPartDestination(rModelPartDestination),
                           mMapperSettings(JsonParameters)
    {
        mpInterfaceVectorContainerOrigin = Kratos::make_unique<InterfaceVectorContainerType>(rModelPartOrigin);
        mpInterfaceVectorContainerDestination = Kratos::make_unique<InterfaceVectorContainerType>(rModelPartDestination);

        this->InitializeInterface();
    }

    /// Destructor.
    ~CouplingGeomteryMapper() override = default;

    ///@}
    ///@name Operations
    ///@{

    void UpdateInterface(
        Kratos::Flags MappingOptions,
        double SearchRadius) override
    {
        KRATOS_ERROR << "Not implemented!" << std::endl;
    }

    void Map(
        const Variable<double>& rOriginVariable,
        const Variable<double>& rDestinationVariable,
        Kratos::Flags MappingOptions) override
    {
        MapInternal(rOriginVariable, rDestinationVariable, MappingOptions);
    }

    void Map(
        const Variable< array_1d<double, 3> >& rOriginVariable,
        const Variable< array_1d<double, 3> >& rDestinationVariable,
        Kratos::Flags MappingOptions) override
    {
        MapInternal(rOriginVariable, rDestinationVariable, MappingOptions);
    }

    void InverseMap(
        const Variable<double>& rOriginVariable,
        const Variable<double>& rDestinationVariable,
        Kratos::Flags MappingOptions) override
    {
        KRATOS_ERROR << "Not implemented!" << std::endl;
    }

    void InverseMap(
        const Variable< array_1d<double, 3> >& rOriginVariable,
        const Variable< array_1d<double, 3> >& rDestinationVariable,
        Kratos::Flags MappingOptions) override
    {
        KRATOS_ERROR << "Not implemented!" << std::endl;
    }

    ///@}
    ///@name Access
    ///@{

    TMappingMatrixType* pGetMappingMatrix() override
    {
        KRATOS_ERROR << "Not implemented!" << std::endl;
    }

    MapperUniquePointerType Clone(ModelPart& rModelPartOrigin,
                                  ModelPart& rModelPartDestination,
                                  Parameters JsonParameters) const override
    {
        return Kratos::make_unique<CouplingGeomteryMapper<TSparseSpace, TDenseSpace>>(
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
        return "CouplingGeomteryMapper";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "CouplingGeomteryMapper";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        BaseType::PrintData(rOStream);
    }

private:

    ///@name Private Operations
    ///@{

    ModelPart& mrModelPartOrigin;
    ModelPart& mrModelPartDestination;

    Parameters mMapperSettings;

    MapperUniquePointerType mpInverseMapper = nullptr;

    TMappingMatrixUniquePointerType mpMappingMatrix;

    MapperLocalSystemPointerVector mMapperLocalSystems;

    InterfaceVectorContainerPointerType mpInterfaceVectorContainerOrigin;
    InterfaceVectorContainerPointerType mpInterfaceVectorContainerDestination;


    void InitializeInterface(Kratos::Flags MappingOptions = Kratos::Flags());

    void BuildMappingMatrix(Kratos::Flags MappingOptions = Kratos::Flags());

    void AssignInterfaceEquationIds()
    {
        MapperUtilities::AssignInterfaceEquationIds(mrModelPartOrigin.GetCommunicator());
        MapperUtilities::AssignInterfaceEquationIds(mrModelPartDestination.GetCommunicator());
    }

    void MapInternal(const Variable<double>& rOriginVariable,
                     const Variable<double>& rDestinationVariable,
                     Kratos::Flags MappingOptions);

    void MapInternal(const Variable<array_1d<double, 3>>& rOriginVariable,
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

    Parameters GetMapperDefaultSettings() const
    {
        // @tobiasteschemachen
        return Parameters( R"({
            "echo_level" : 0
        })");
    }

    ///@}

}; // Class CouplingGeomteryMapper

///@} addtogroup block
}  // namespace Kratos.

#endif // KRATOS_COUPLING_GEOMETRY_MAPPER_H_INCLUDED defined
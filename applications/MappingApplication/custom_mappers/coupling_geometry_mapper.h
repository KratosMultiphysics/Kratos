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

#if !defined(KRATOS_COUPLING_GEOMETRY_MAPPER_H_INCLUDED)
#define  KRATOS_COUPLING_GEOMETRY_MAPPER_H_INCLUDED

// System includes

// External includes

// Project includes
#include "mapper.h"


namespace Kratos
{
///@name Kratos Classes
///@{

class CouplingGeomteryLocalSystem : public MapperLocalSystem
{
public:

    explicit CouplingGeomteryLocalSystem(GeometryPointerType pGeom) : mpGeom(pGeom) {}

    void CalculateAll(MatrixType& rLocalMappingMatrix,
                      EquationIdVectorType& rOriginIds,
                      EquationIdVectorType& rDestinationIds,
                      MapperLocalSystem::PairingStatus& rPairingStatus) const override;

    CoordinatesArrayType& Coordinates() const override
    {
        KRATOS_DEBUG_ERROR_IF_NOT(mpGeom) << "Members are not intitialized!" << std::endl;
        // return mpGeom->Center();
    }

    /// Turn back information as a string.
    std::string PairingInfo(const int EchoLevel) const override;

private:
    GeometryPointerType mpGeom;

};

/// Nearest Neighbor Mapper
/** This class implements the Nearest Neighbor Mapping technique.
* Each node on the destination side gets assigned is's closest neighbor on the other side of the interface.
* In the mapping phase every node gets assigned the value of it's neighbor
* For information abt the available echo_levels and the JSON default-parameters
* look into the class description of the MapperCommunicator
*/
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
    typedef typename BaseType::MapperUniquePointerType MapperUniquePointerType;
    // typedef typename BaseType::MapperInterfaceInfoUniquePointerType MapperInterfaceInfoUniquePointerType;

    typedef typename BaseType::TMappingMatrixType TMappingMatrixType;

    ///@}
    ///@name Life Cycle
    ///@{

    // Default constructor, needed for registration
    CouplingGeomteryMapper(ModelPart& rModelPartOrigin,
                          ModelPart& rModelPartDestination){}

    CouplingGeomteryMapper(ModelPart& rModelPartOrigin,
                          ModelPart& rModelPartDestination,
                          Parameters JsonParameters)
    {
        // this->ValidateInput();
        // this->Initialize();
    }

    /// Destructor.
    ~CouplingGeomteryMapper() override = default;

    ///@}
    ///@name Operations
    ///@{

    virtual void UpdateInterface(
        Kratos::Flags MappingOptions,
        double SearchRadius) override {}

    /**
    * @brief Mapping from Origin to Destination, Scalar Variable
    * Data is exchanged on the Interface, from the Origin-Modelpart
    * to the Destination-ModelPart (the modelparts were specified in the
    * construction Phase of the Mapper)
    * @param rOriginVariable Variable on the Origin-ModelPart
    * @param rDestinationVariable Variable on the Destination-ModelPart
    * @param MappingOptions flags used to specify options for the mapping
    * @see InverseMap
    */
    virtual void Map(
        const Variable<double>& rOriginVariable,
        const Variable<double>& rDestinationVariable,
        Kratos::Flags MappingOptions) override {}

    /**
    * @brief Mapping from Origin to Destination, Vector Variable
    * Same as Map, but maps an array3-variable
    * @see Map
    */
    virtual void Map(
        const Variable< array_1d<double, 3> >& rOriginVariable,
        const Variable< array_1d<double, 3> >& rDestinationVariable,
        Kratos::Flags MappingOptions) override {}

    /**
    * @brief Mapping from Destination to Origin, Scalar Variable
    * Data is exchanged on the Interface, from the Destination-Modelpart
    * to the Origin-ModelPart (the modelparts were specified in the
    * construction Phase of the Mapper)
    * It does the opposite of Map
    * @param rOriginVariable Variable on the Origin-ModelPart
    * @param rDestinationVariable Variable on the Destination-ModelPart
    * @param MappingOptions flags used to specify options for the mapping
    * @see Map
    */
    virtual void InverseMap(
        const Variable<double>& rOriginVariable,
        const Variable<double>& rDestinationVariable,
        Kratos::Flags MappingOptions) override {}

    /**
    * @brief Mapping from Destination to Origin, Vector Variable
    * Same as InveseMap, but maps an array3-variable
    * @see InverseMap
    */
    virtual void InverseMap(
        const Variable< array_1d<double, 3> >& rOriginVariable,
        const Variable< array_1d<double, 3> >& rDestinationVariable,
        Kratos::Flags MappingOptions)override
        {}

    ///@}
    ///@name Access
    ///@{

    /**
     * @brief This method returns the mapping-matrix
     * @return The mapping-matrix
     */
    virtual TMappingMatrixType* pGetMappingMatrix() override {return nullptr;}

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

    void CreateMapperLocalSystems(
        const Communicator& rModelPartCommunicator,
        std::vector<Kratos::unique_ptr<MapperLocalSystem>>& rLocalSystems)
    {
        // MapperUtilities::CreateMapperLocalSystemsFromNodes<CouplingGeomteryLocalSystem>(
        //     rModelPartCommunicator,
        //     rLocalSystems);
    }

    // MapperInterfaceInfoUniquePointerType GetMapperInterfaceInfo() const override
    // {
    //     return Kratos::make_unique<CouplingGeomteryInterfaceInfo>();
    // }

    Parameters GetMapperDefaultSettings() const
    {
        return Parameters( R"({
            "search_radius"            : -1.0,
            "search_iterations"        : 3,
            "echo_level"               : 0
        })");
    }

    ///@}

}; // Class CouplingGeomteryMapper

///@} addtogroup block
}  // namespace Kratos.

#endif // KRATOS_COUPLING_GEOMETRY_MAPPER_H_INCLUDED  defined
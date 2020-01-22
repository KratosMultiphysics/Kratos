//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Tobias Teschemacher
//                   Pooyan Dadvand
//                   Philipp Bucher
//

#if !defined(KRATOS_COUPLING_GEOMETRY_H_INCLUDED )
#define  KRATOS_COUPLING_GEOMETRY_H_INCLUDED

// System includes

// External includes

// Project includes
#include "geometry.h"


namespace Kratos
{
///@name Kratos Classes
///@{

/**
 * @class CouplingGeometry
 * @ingroup KratosCore
 * @brief The CouplingGeometry can connect different geometries, those
 can be of different entities but have to coincide with the dimension.
 */
template<class TPointType> class CouplingGeometry
    : public Geometry<TPointType>
{
public:
    ///@}
    ///@name Type Definitions
    ///@{

    /// Geometry as base class.
    typedef Geometry<TPointType> BaseType;
    typedef Geometry<TPointType> GeometryType;

    typedef typename GeometryType::Pointer GeometryPointer;
    typedef std::vector<GeometryPointer> GeometryPointerVector;

    /// Pointer definition of CouplingGeometry
    KRATOS_CLASS_POINTER_DEFINITION( CouplingGeometry );

    typedef TPointType PointType;

    typedef typename BaseType::IndexType IndexType;
    typedef typename BaseType::SizeType SizeType;

    typedef typename BaseType::PointsArrayType PointsArrayType;

    ///@}
    ///@name Public Static Members
    ///@{

    static constexpr IndexType Master = 0;
    static constexpr IndexType Slave = 1;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor for coupling one master to one slave geometry.
    CouplingGeometry(
        GeometryPointer pMasterGeometry,
        GeometryPointer pSlaveGeometry)
        : BaseType(PointsArrayType(), &(pMasterGeometry->GetGeometryData()))
    {
        KRATOS_DEBUG_ERROR_IF(pMasterGeometry->Dimension() != pSlaveGeometry->Dimension())
            << "Geometries of different dimensional size!" << std::endl;

        mpGeometries.resize(2);

        mpGeometries[0] = pMasterGeometry;
        mpGeometries[1] = pSlaveGeometry;
    }

    /// Constructor for coupling multiple points.
    CouplingGeometry(
        GeometryPointerVector GeometryPointerVector)
        : BaseType(PointsArrayType(), &(GeometryPointerVector[0]->GetGeometryData()))
        , mpGeometries(GeometryPointerVector)
    {
    }

    explicit CouplingGeometry()
        : BaseType()
    {
    }

    /**
     * Copy constructor.
     * Construct this geometry as a copy of given geometry.
     *
     * @note This copy constructor does not copy the points and new
     * geometry shares points with given source geometry. It is
     * obvious that any change to this new geometry's point affect
     * source geometry's points too.
     */
    CouplingGeometry( CouplingGeometry const& rOther )
        : BaseType( rOther )
        , mpGeometries(rOther.mpGeometries)
    {
    }

    /**
     * Copy constructor from a geometry with other point type.
     * Construct this geometry as a copy of given geometry which
     * has different type of points. The given goemetry's
     * TOtherPointType* must be implicity convertible to this
     * geometry PointType.
     *
     * @note This copy constructor does not copy the points and new
     * geometry shares points with given source geometry. It is
     * obvious that any change to this new geometry's point affect
     * source geometry's points too.
     */
    template<class TOtherPointType> explicit CouplingGeometry(
        CouplingGeometry<TOtherPointType> const& rOther )
        : BaseType( rOther )
        , mpGeometries( rOther.mpGeometries )
    {
    }

    /// Destructor
    ~CouplingGeometry() override = default;

    ///@}
    ///@name Operators
    ///@{

    /**
     * Assignment operator.
     *
     * @note This operator don't copy the points and this
     * geometry shares points with given source geometry. It's
     * obvious that any change to this geometry's point affect
     * source geometry's points too.
     *
     * @see Clone
     * @see ClonePoints
     */
    CouplingGeometry& operator=( const CouplingGeometry& rOther )
    {
        BaseType::operator=( rOther );
        mpGeometries = rOther.mpGeometries;
        return *this;
    }

    /**
     * @brief Assignment operator for geometries with different point type.
     *
     * @note This operator don't copy the points and this
     * geometry shares points with given source geometry. It's
     * obvious that any change to this geometry's point affect
     * source geometry's points too.
     *
     * @see Clone
     * @see ClonePoints
     */
    template<class TOtherPointType>
    CouplingGeometry& operator=(
        CouplingGeometry<TOtherPointType> const & rOther )
    {
        BaseType::operator=( rOther );
        mpGeometries = rOther.mpGeometries;
        return *this;
    }

    ///@}
    ///@name Operations
    ///@{

    typename BaseType::Pointer Create(
        PointsArrayType const& ThisPoints ) const override
    {
        return Kratos::make_shared<CouplingGeometry>();
    }

    ///@}
    ///@name Input and output
    ///@{

    /**
    * @brief This function returns the pointer of the geometry part
    *        which is corresponding to the index.
    *        Checks if index is available only in debug.
    * @param Index: 0 -> master, all bigger than 0 -> slaves.
    * @return pointer of geometry, corresponding to the index.
    */
    GeometryPointer pGetGeometryPart(IndexType Index) override
    {
        KRATOS_DEBUG_ERROR_IF(mpGeometries.size() <= Index) << "Index "
            << Index << " out of range. CouplingGeometry #" << this->Id()
            << " has " << mpGeometries.size() << " geometries." << std::endl;

        return mpGeometries[Index];
    }

    /**
    * @brief This function returns the const pointer of the geometry part
    *        which is corresponding to the index.
    *        Checks if index is available only in debug.
    * @param Index: 0 -> master, all bigger than 0 -> slaves.
    * @return const pointer of geometry, corresponding to the index.
    */
    const GeometryPointer pGetGeometryPart(IndexType Index) const override
    {
        KRATOS_DEBUG_ERROR_IF(mpGeometries.size() <= Index) << "Index \""
            << Index << "\" out of range. CouplingGeometry #" << this->Id() 
            << " has " << mpGeometries.size() << " geometries." << std::endl;

        return mpGeometries[Index];
    }

    /**
    * @brief ONLY for coupling_geometry. Not necessary in base class.
    * @details Allows to exchange certain geometries.
    * @param Index of the geometry part. 0->Master; 1->Slave
     */
    void SetGeometryPart(IndexType Index, GeometryPointer pGeometry)
    {
        KRATOS_DEBUG_ERROR_IF(mpGeometries.size() <= Index) << "Index out of range: "
            << Index << " composite contains only of: "
            << mpGeometries.size() << " geometries." << std::endl;

        if (0 == Index){
            if (mpGeometries.size() > 1) {
                KRATOS_ERROR_IF(pGeometry->Dimension() != mpGeometries[1]->Dimension())
                    << "Dimension of new master geometry does not coincide with the other geometries. "
                    << "Dimension of new geometry: " << pGeometry->Dimension()
                    << ", dimension of coupling geometry: " << mpGeometries[1]->Dimension() << std::endl;
            }

            this->SetGeometryData(&(pGeometry->GetGeometryData()));
        }

        KRATOS_ERROR_IF(pGeometry->Dimension() != mpGeometries[0]->Dimension())
            << "Dimension of new entity does not coincide with this coupling geometry. "
            << "Dimension of new geometry: " << pGeometry->Dimension()
            << ", dimension of coupling geometry: " << this->Dimension() << std::endl;

        mpGeometries[Index] = pGeometry;
    }

    /**
    * @brief ONLY for coupling_geometry. Not necessary in base class.
    * @details Allows to enhance the coupling geometry, with another geometry.
    * @param Index of the geometry part. 0->Master; 1->Slave
     */
    IndexType AddGeometryPart(GeometryPointer pGeometry)
    {
        KRATOS_DEBUG_ERROR_IF(mpGeometries[0]->Dimension() != pGeometry->Dimension())
            << "Geometries of different dimensional size!" << std::endl;

        KRATOS_ERROR_IF(pGeometry->Dimension() != mpGeometries[0]->Dimension())
            << "Dimension of new entity does not coincide with this coupling geometry. "
            << "Dimension of new geometry: " << pGeometry->Dimension()
            << ", dimension of coupling geometry: " << this->Dimension() << std::endl;

        IndexType new_index = mpGeometries.size();

        mpGeometries.push_back(pGeometry);

        return new_index;
    }

    /**
     * @brief The number of geometry part
     * @return The number of geometry parts that this geometry contains.
     */
    SizeType NumberOfGeometryParts() const override
    {
        return mpGeometries.size();
    }

    ///@}
    ///@name Geometry Information
    ///@{

    /**
     * @brief Returns the domain size of the master geometry.
     * @return The domaon size of the master geometry
     */
    double DomainSize() const override
    {
        KRATOS_DEBUG_ERROR_IF(mpGeometries.size() == 0) << "No master assigned. Geometry vector of size 0." << std::endl;

        return mpGeometries[0]->DomainSize();
    }

    /**
     * @brief Returns the center of the master geometry.
     * @return The center of the master geometry
     */
    Point Center() const override
    {
        KRATOS_DEBUG_ERROR_IF(mpGeometries.size() == 0) << "No master assigned. Geometry vector of size 0." << std::endl;

        return mpGeometries[0]->Center();
    }

    ///@}
    ///@name Information
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "Coupling geometry that holds a master and a set of slave geometries.";
    }

    /// Print information about this object.
    void PrintInfo( std::ostream& rOStream ) const override
    {
        rOStream << "Coupling geometry that holds a master and a set of slave geometries.";
    }

    /// Print object's data.
    void PrintData( std::ostream& rOStream ) const override
    {
        BaseType::PrintData( rOStream );
        std::cout << std::endl;
        rOStream << "    CouplingGeometry with " << mpGeometries.size() << " geometries.";
    }

    ///@}
private:
    ///@name Private Member Variables
    ///@{

    GeometryPointerVector mpGeometries;

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save( Serializer& rSerializer ) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, BaseType );
        rSerializer.save("Geometries", mpGeometries);
    }

    void load( Serializer& rSerializer ) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BaseType );
        rSerializer.load("Geometries", mpGeometries);
    }

    ///@}
    ///@name Private Friends
    ///@{

    template<class TOtherPointType> friend class CouplingGeometry;

    ///@}
}; // Class Geometry

///@}
///@name Input and output
///@{
/**
 * input stream functions
 */
template<class TPointType> inline std::istream& operator >> (
    std::istream& rIStream,
    CouplingGeometry<TPointType>& rThis );
/**
 * output stream functions
 */
template<class TPointType> inline std::ostream& operator << (
    std::ostream& rOStream,
    const CouplingGeometry<TPointType>& rThis )
{
    rThis.PrintInfo( rOStream );
    rOStream << std::endl;
    rThis.PrintData( rOStream );
    return rOStream;
}

///@}
}// namespace Kratos.

#endif // KRATOS_COUPLING_GEOMETRY_H_INCLUDED  defined

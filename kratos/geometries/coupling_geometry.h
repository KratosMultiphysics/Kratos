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
    GeometryPointer pGetGeometryPart(const IndexType Index) override
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
    const GeometryPointer pGetGeometryPart(const IndexType Index) const override
    {
        KRATOS_DEBUG_ERROR_IF(mpGeometries.size() <= Index) << "Index \""
            << Index << "\" out of range. CouplingGeometry #" << this->Id()
            << " has " << mpGeometries.size() << " geometries." << std::endl;

        return mpGeometries[Index];
    }

    /**
     * @brief Allows to exchange certain geometries.
     * @param Index of the geometry part. 0->Master; 1->Slave
     * @param pGeometry The new geometry to add
     */
    void SetGeometryPart(
        const IndexType Index,
        GeometryPointer pGeometry
        ) override
    {
        KRATOS_DEBUG_ERROR_IF(mpGeometries.size() <= Index) << "Index out of range: "
            << Index << " composite contains only of: "
            << mpGeometries.size() << " geometries." << std::endl;

        if (Index == 0){
            this->SetGeometryData(&(pGeometry->GetGeometryData()));
        }

        mpGeometries[Index] = pGeometry;
    }

    /**
     * @brief Allows to enhance the coupling geometry, with another geometry.
     * @param pGeometry The new geometry to add
     */
    IndexType AddGeometryPart(GeometryPointer pGeometry) override
    {
        const IndexType new_index = mpGeometries.size();

        mpGeometries.push_back(pGeometry);

        return new_index;
    }

    /**
     * @brief Removes a geometry part
     * @param pGeometry The new geometry to remove
     */
    void RemoveGeometryPart(GeometryPointer pGeometry) override
    {
        const IndexType geometry_id = pGeometry->Id();
        IndexType to_remove_id = 0;
        for (const auto& p_geom : mpGeometries) {
            if (p_geom->Id() == geometry_id) {
                break;
            }
            ++to_remove_id;
        }

        RemoveGeometryPart(to_remove_id);
    }

    /**
     * @brief Removes a geometry part
     * @param Index of the geometry part.
     */
    void RemoveGeometryPart(const IndexType Index) override
    {
        KRATOS_ERROR_IF(Index == 0) << "Master geometry should not be removed from the CouplingGeometry" << std::endl;

        const SizeType number_of_geometries = NumberOfGeometryParts();
        for (IndexType i = Index; i < number_of_geometries - 1; ++i) {
            mpGeometries[i] = mpGeometries[i + 1];
        }
        mpGeometries[number_of_geometries - 1] = nullptr;
        mpGeometries.erase(mpGeometries.begin() + number_of_geometries - 1);
    }

    /**
     * @brief Use to check if certain Indexed object is within the geometry parts of this geometry.
     * @param Index of the geometry part. This index can be used differently within the derived classes.
     * @return true if has geometry part
     */
    bool HasGeometryPart(const IndexType Index) const override
    {
        if (Index < NumberOfGeometryParts()) {
            return true;
        } else {
            return false;
        }
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

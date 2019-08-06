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
//  contributors:    Pooyan Dadvand
//                   Philipp Bucher
//

#if !defined(KRATOS_COUPLING_MASTER_SLAVE_H_INCLUDED )
#define  KRATOS_COUPLING_MASTER_SLAVE_H_INCLUDED

// System includes
#include <vector>

// External includes

// Project includes
#include "includes/define.h"
#include "geometry.h"


namespace Kratos
{
///@name Kratos Classes
///@{

/**
 * @class CouplingMasterSlave
 * @ingroup KratosCore
 * @brief The CouplingMasterSlave can connect .
 */
template<class TPointType> class CouplingMasterSlave
    : public Geometry<TPointType>
{
public:
    ///@}
    ///@name Type Definitions
    ///@{

    /**
     * Geometry as base class.
     */
    typedef Geometry<TPointType> BaseType;
    typedef Geometry<TPointType> GeometryType;

    /**
     * Pointer definition of CouplingMasterSlave
     */
    KRATOS_CLASS_POINTER_DEFINITION( CouplingMasterSlave );

    typedef TPointType PointType;

    typedef GeometryData::IntegrationMethod IntegrationMethod;

    typedef typename BaseType::GeometriesArrayType GeometriesArrayType;
    typedef typename BaseType::IndexType IndexType;
    typedef typename BaseType::SizeType SizeType;
    typedef typename BaseType::PointsArrayType PointsArrayType;
    typedef typename BaseType::CoordinatesArrayType CoordinatesArrayType;

    typedef typename BaseType::IntegrationPointType IntegrationPointType;
    typedef typename BaseType::IntegrationPointsArrayType IntegrationPointsArrayType;
    typedef typename BaseType::IntegrationPointsContainerType IntegrationPointsContainerType;

    typedef typename BaseType::ShapeFunctionsValuesContainerType ShapeFunctionsValuesContainerType;
    typedef typename BaseType::ShapeFunctionsLocalGradientsContainerType ShapeFunctionsLocalGradientsContainerType;
    typedef typename BaseType::ShapeFunctionsGradientsType ShapeFunctionsGradientsType;
    typedef typename BaseType::ShapeFunctionsSecondDerivativesType ShapeFunctionsSecondDerivativesType;
    typedef typename BaseType::ShapeFunctionsThirdDerivativesType ShapeFunctionsThirdDerivativesType;

    ///@}
    ///@name Life Cycle
    ///@{

    CouplingMasterSlave(
        typename GeometryType::Pointer pMasterGeometry,
        typename GeometryType::Pointer pSlaveGeometry,
        bool IsGeometryDataInitialized = false)
        : BaseType(PointsArrayType(), &(pMasterGeometry->GetGeometryData())),
        mIsGeometryDataInitialized(IsGeometryDataInitialized)
    {
        KRATOS_DEBUG_ERROR_IF(pMasterGeometry->Dimension() != pSlaveGeometry->Dimension())
            << "Geometries of different dimensional size!" << std::endl;

        mpGeometries.resize(2);

        mpGeometries[0] = pMasterGeometry;
        mpGeometries[1] = pSlaveGeometry;
    }

    explicit CouplingMasterSlave(const PointsArrayType& ThisPoints)
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
    CouplingMasterSlave( CouplingMasterSlave const& rOther )
        : BaseType( rOther )
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
    template<class TOtherPointType> explicit CouplingMasterSlave( CouplingMasterSlave<TOtherPointType> const& rOther )
        : BaseType( rOther )
    {
    }

    /**
     * Destructor. Does nothing!!!
     */
    ~CouplingMasterSlave() override {}

    GeometryData::KratosGeometryFamily GetGeometryFamily() const override
    {
        return GeometryData::Kratos_generic_family;
    }

    GeometryData::KratosGeometryType GetGeometryType() const override
    {
        return GeometryData::Kratos_generic_type;
    }

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
    CouplingMasterSlave& operator=( const CouplingMasterSlave& rOther )
    {
        BaseType::operator=( rOther );
        return *this;
    }

    /**
     * Assignment operator for geometries with different point type.
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
    CouplingMasterSlave& operator=( CouplingMasterSlave<TOtherPointType> const & rOther )
    {
        BaseType::operator=( rOther );
        return *this;
    }

    ///@}
    ///@name Operations
    ///@{

    typename BaseType::Pointer Create( PointsArrayType const& ThisPoints ) const override
    {
        return typename BaseType::Pointer( new CouplingMasterSlave( ThisPoints) );
    }

    ///@}

    ///@name Shape Function
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    typename GeometryType& GetGeometryPart(IndexType Index) const override
    {
        KRATOS_DEBUG_ERROR_IF(mpGeometries.size() < Index) << "Index out of range: "
            << Index << " composite contains only of: "
            << mpGeometries.size() << " geometries." << std::endl;

        return *mpGeometries[Index];
    }

    ///@}

    /** Returns the domain size of the master geometry.
    */
    double DomainSize() const override
    {
        return mpGeometries[0]->DomainSize();
    }

    /** Returns the center of the master geometry.
    */
    Point Center() const override
    {
        return mpGeometries[0]->Center();
    }

    ///@name Input and output
    ///@{

    /**
     * Turn back information as a string.
     *
     * @return String contains information about this geometry.
     * @see PrintData()
     * @see PrintInfo()
     */
    std::string Info() const override
    {
        return "Composite geometry that holds a master and a set of slave geometries.";
    }

    /**
     * Print information about this object.
     * @param rOStream Stream to print into it.
     * @see PrintData()
     * @see Info()
     */
    void PrintInfo( std::ostream& rOStream ) const override
    {
        rOStream << "Composite geometry that holds a master and a set of slave geometries.";
    }

    /**
     * Print geometry's data into given stream.
     * Prints it's points
     * by the order they stored in the geometry and then center
     * point of geometry.
     *
     * @param rOStream Stream to print into it.
     * @see PrintInfo()
     * @see Info()
     */
    void PrintData( std::ostream& rOStream ) const override
    {
        BaseType::PrintData( rOStream );
        std::cout << std::endl;
        rOStream << "    CouplingGeometry with " << mpGeometries << " geometries.";
    }

    ///@}
    ///@name Friends
    ///@{

    ///@}

protected:
    /**
     * There are no protected members in class CouplingMasterSlave
     */

private:
    ///@name Static Member Variables
    ///@{

    std::vector<typename GeometryType::Pointer> mpGeometries;

    bool mIsGeometryDataInitialized;

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save( Serializer& rSerializer ) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, BaseType );
        rSerializer.save("IsGeometryDataInitialized", mIsGeometryDataInitialized);
    }

    void load( Serializer& rSerializer ) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BaseType );
        rSerializer.load("IsGeometryDataInitialized", mIsGeometryDataInitialized);
    }

    CouplingMasterSlave() : BaseType() {}

    ///@}
    ///@name Member Variables
    ///@{


    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Private Friends
    ///@{

    template<class TOtherPointType> friend class CouplingMasterSlave;

    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}
}; // Class Geometry

///@}
///@name Type Definitions
///@{


///@}
///@name Input and output
///@{
/**
 * input stream functions
 */
template<class TPointType> inline std::istream& operator >> (
    std::istream& rIStream,
    CouplingMasterSlave<TPointType>& rThis );
/**
 * output stream functions
 */
template<class TPointType> inline std::ostream& operator << (
    std::ostream& rOStream,
    const CouplingMasterSlave<TPointType>& rThis )
{
    rThis.PrintInfo( rOStream );
    rOStream << std::endl;
    rThis.PrintData( rOStream );
    return rOStream;
}

///@}
}// namespace Kratos.

#endif // KRATOS_COUPLING_MASTER_SLAVE_H_INCLUDED  defined

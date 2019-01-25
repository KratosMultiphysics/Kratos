//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//                   Janosch Stascheit
//                   Felix Nagel
//  contributors:    Hoang Giang Bui
//                   Josep Maria Carbonell
//

#if !defined(KRATOS_LINE_3D_N_H_INCLUDED )
#define  KRATOS_LINE_3D_N_H_INCLUDED

// System includes

// External includes

// Project includes
#include "geometries/geometry.h"
#include "integration/line_gauss_legendre_integration_points.h"


namespace Kratos
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

template<class TPointType>

class Line3DN : public Geometry<TPointType>
{
public:
    ///@}
    ///@name Type Definitions
    ///@{

    /// Geometry as base class.
    typedef Geometry<TPointType> BaseType;

    /// Pointer definition of Line3D3
    KRATOS_CLASS_POINTER_DEFINITION( Line3DN );

    /** Integration methods implemented in geometry.
    */
    typedef GeometryData::IntegrationMethod IntegrationMethod;

    /** A Vector of counted pointers to Geometries. Used for
    returning edges of the geometry.
     */
    typedef typename BaseType::GeometriesArrayType GeometriesArrayType;

    /** Redefinition of template parameter TPointType.
     */
    typedef TPointType PointType;

    /** Type used for indexing in geometry class.std::size_t used for indexing
    point or integration point access methods and also all other
    methods which need point or integration point index.
    */
    typedef typename BaseType::IndexType IndexType;


    /** This typed used to return size or dimension in
    geometry. Dimension, WorkingDimension, PointsNumber and
    ... return this type as their results.
    */
    typedef typename BaseType::SizeType SizeType;

    /** Array of counted pointers to point. This type used to hold
    geometry's points.
    */
    typedef  typename BaseType::PointsArrayType PointsArrayType;

    /** This type used for representing an integration point in
    geometry. This integration point is a point with an
    additional weight component.
    */
    typedef typename BaseType::IntegrationPointType IntegrationPointType;

    /** A Vector of IntegrationPointType which used to hold
    integration points related to an integration
    method. IntegrationPoints functions used this type to return
    their results.
    */
    typedef typename BaseType::IntegrationPointsArrayType IntegrationPointsArrayType;

    /** A Vector of IntegrationPointsArrayType which used to hold
    integration points related to different integration method
    implemented in geometry.
    */
    typedef typename BaseType::IntegrationPointsContainerType IntegrationPointsContainerType;

    /** A third order tensor used as shape functions' values
    continer.
    */
    typedef typename BaseType::ShapeFunctionsValuesContainerType ShapeFunctionsValuesContainerType;

    /** A fourth order tensor used as shape functions' local
    gradients container in geometry.
    */
    typedef typename BaseType::ShapeFunctionsLocalGradientsContainerType ShapeFunctionsLocalGradientsContainerType;

    /** A third order tensor to hold jacobian matrices evaluated at
    integration points. Jacobian and InverseOfJacobian functions
    return this type as their result.
    */
    typedef typename BaseType::JacobiansType JacobiansType;

    /** A third order tensor to hold shape functions' local
    gradients. ShapefunctionsLocalGradients function return this
    type as its result.
    */
    typedef typename BaseType::ShapeFunctionsGradientsType ShapeFunctionsGradientsType;

    /** Type of the normal vector used for normal to edges in geomety.
     */
    typedef typename BaseType::NormalType NormalType;

    /**
     * Type of coordinates array
     */
    typedef typename BaseType::CoordinatesArrayType CoordinatesArrayType;

    ///@}
    ///@name Life Cycle
    ///@{

//     Line3D3( const PointType& FirstPoint, const PointType& SecondPoint, const PointType& ThirdPoint )
//         : BaseType( PointsArrayType(), &msGeometryData )
//     {
//         BaseType::Points().push_back( typename PointType::Pointer( new PointType( FirstPoint ) ) );
//         BaseType::Points().push_back( typename PointType::Pointer( new PointType( SecondPoint ) ) );
//         BaseType::Points().push_back( typename PointType::Pointer( new PointType( ThirdPoint ) ) );
//     }


    Line3DN( const PointsArrayType& ThisPoints )
        : BaseType( ThisPoints, &msGeometryData )
    {
        KRATOS_INFO("line3DN") << "created with " << BaseType::PointsNumber() << " nodes" << std::endl;
    }

    /** Copy constructor.
     * Construct this geometry as a copy of given geometry.
     * @note This copy constructor don't copy the points and new
     * geometry shares points with given source geometry. It's
     * obvious that any change to this new geometry's point affect
     * source geometry's points too.
     */
    Line3DN( Line3DN const& rOther )
        : BaseType( rOther )
    {
    }

    /** Copy constructor from a geometry with other point type.
     * Construct this geometry as a copy of given geometry which
     * has different type of points. The given goemetry's
     * TOtherPointType* must be implicity convertible to this
     * geometry PointType.
     * @note This copy constructor don't copy the points and new
     * geometry shares points with given source geometry. It's
     * obvious that any change to this new geometry's point affect
     * source geometry's points too.
     */
    template<class TOtherPointType> Line3DN( Line3DN<TOtherPointType> const& rOther )
        : BaseType( rOther )
    {
    }

    /// Destructor. Do nothing!!!
    ~Line3DN() override {}

    GeometryData::KratosGeometryFamily GetGeometryFamily() const override
    {
        return GeometryData::Kratos_Linear;
    }

    GeometryData::KratosGeometryType GetGeometryType() const override
    {
        return GeometryData::Kratos_Line3DN;
    }

    ///@}
    ///@name Operators
    ///@{

    /** Assignment operator.

    @note This operator don't copy the points and this
    geometry shares points with given source geometry. It's
    obvious that any change to this geometry's point affect
    source geometry's points too.

    @see Clone
    @see ClonePoints
    */
    Line3DN& operator=( const Line3DN& rOther )
    {
        BaseType::operator=( rOther );
        return *this;
    }

    /** Assignment operator for geometries with different point type.
     * @note This operator don't copy the points and this
     * geometry shares points with given source geometry. It's
     * obvious that any change to this geometry's point affect
     * source geometry's points too.
     * @see Clone
     * @see ClonePoints
     */
    template<class TOtherPointType>
    Line3DN& operator=( Line3DN<TOtherPointType> const & rOther )
    {
        BaseType::operator=( rOther );
        return *this;
    }

    ///@}
    ///@name Operations
    ///@{

    typename BaseType::Pointer Create( PointsArrayType const& ThisPoints ) const override
    {
      // line 3 conectivities order 1-3-2
      return typename BaseType::Pointer( new Line3DN( ThisPoints ) );
    }


    std::string Info() const override
    {
        return "1 dimensional line with n nodes in 3D space";
    }

    /** Print information about this object.

    @param rOStream Stream to print into it.
    @see PrintData()
    @see Info()
    */
    void PrintInfo( std::ostream& rOStream ) const override
    {
        rOStream << "1 dimensional line with n nodes in 3D space";
    }


protected:


private:
    ///@name Static Member Variables
    ///@{

    static const GeometryData msGeometryData;

    ///@}
    ///@name Member Variables
    ///@{


    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save( Serializer& rSerializer ) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, BaseType );
    }

    void load( Serializer& rSerializer ) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BaseType );
    }

    Line3DN(): BaseType( PointsArrayType(), &msGeometryData ) {}


    template<class TOtherPointType> friend class Line3DN;

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


/// input stream function
template<class TPointType>
inline std::istream& operator >> ( std::istream& rIStream,
                                   Line3DN<TPointType>& rThis );

/// output stream function
template<class TPointType>
inline std::ostream& operator << ( std::ostream& rOStream,
                                   const Line3DN<TPointType>& rThis )
{
    rThis.PrintInfo( rOStream );
    rOStream << std::endl;
    rThis.PrintData( rOStream );

    return rOStream;
}

}  // namespace Kratos.

#endif // KRATOS_LINE_3D_3_H_INCLUDED  defined

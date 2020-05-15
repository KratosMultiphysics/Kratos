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
//

#if !defined(KRATOS_QUADRATURE_POINT_PARTITIONED_GEOMETRY_H_INCLUDED )
#define  KRATOS_QUADRATURE_POINT_PARTITIONED_GEOMETRY_H_INCLUDED

// System includes

// External includes

// Project includes
#include "geometries/geometry.h"
#include "geometries/geometry_dimension.h"
#include "geometries/quadrature_point_geometry.h"


namespace Kratos
{
/**
 * @class QuadraturePointGeometry
 * @ingroup KratosCore
 * @brief A sinlge quadrature point, that can be used for geometries without
 *        a predefined integration scheme, i.e. they can handle material point elements,
 *        isogeometric analysis elements or standard finite elements which are defined
 *        at a single quadrature point.
 *        Shape functions and integration types are precomputed and are set from
 *        from outside.
 *        The parent pointer can provide the adress to the owner of this quadrature point.
 */
template<class TPointType,
    int TWorkingSpaceDimension,
    int TLocalSpaceDimension = TWorkingSpaceDimension,
    int TDimension = TLocalSpaceDimension>
class QuadraturePointPartitionedGeometry
    : public QuadraturePointGeometry<TPointType, TWorkingSpaceDimension, TLocalSpaceDimension, TDimension>
{
public:

    /// Pointer definition of QuadraturePointGeometry
    KRATOS_CLASS_POINTER_DEFINITION(QuadraturePointPartitionedGeometry);

    typedef QuadraturePointGeometry<TPointType, TWorkingSpaceDimension, TLocalSpaceDimension, TDimension> BaseType;
    typedef Geometry<TPointType> GeometryType;

    typedef typename GeometryType::IndexType IndexType;
    typedef typename GeometryType::SizeType SizeType;

    typedef typename GeometryType::PointsArrayType PointsArrayType;
    typedef typename GeometryType::CoordinatesArrayType CoordinatesArrayType;

    typedef typename GeometryType::IntegrationPointType IntegrationPointType;
    typedef typename GeometryType::IntegrationPointsArrayType IntegrationPointsArrayType;

    typedef typename GeometryData::ShapeFunctionsGradientsType ShapeFunctionsGradientsType;

    typedef GeometryShapeFunctionContainer<GeometryData::IntegrationMethod> GeometryShapeFunctionContainerType;

    typedef typename GeometryType::IntegrationPointsContainerType IntegrationPointsContainerType;
    typedef typename GeometryType::ShapeFunctionsValuesContainerType ShapeFunctionsValuesContainerType;
    typedef typename GeometryType::ShapeFunctionsLocalGradientsContainerType ShapeFunctionsLocalGradientsContainerType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor with points and geometry shape function container
    /*
    QuadraturePointPartitionedGeometry(
        const PointsArrayType& ThisPoints,
        GeometryShapeFunctionContainerType& ThisGeometryShapeFunctionContainer)
        : BaseType(ThisPoints, &mGeometryData)
        , mGeometryData(
            &msGeometryDimension,
            ThisGeometryShapeFunctionContainer)
    {
    }
    */
    QuadraturePointPartitionedGeometry(
        const PointsArrayType& ThisPoints,
        GeometryShapeFunctionContainerType& ThisGeometryShapeFunctionContainer)
        : QuadraturePointGeometry(ThisPoints, ThisGeometryShapeFunctionContainer)
    {
    }

    /// Constructor with points, geometry shape function container, parent
    /*
    QuadraturePointPartitionedGeometry(
        const PointsArrayType& ThisPoints,
        GeometryShapeFunctionContainerType& ThisGeometryShapeFunctionContainer,
        GeometryType* pGeometryParent)
        : BaseType(ThisPoints, &mGeometryData)
        , mGeometryData(
            &msGeometryDimension,
            ThisGeometryShapeFunctionContainer)
        , mpGeometryParent(pGeometryParent)
    {
    }
    */
    QuadraturePointPartitionedGeometry(
    const PointsArrayType& ThisPoints,
    GeometryShapeFunctionContainerType& ThisGeometryShapeFunctionContainer,
    GeometryType* pGeometryParent)
    : QuadraturePointGeometry(ThisPoints, ThisGeometryShapeFunctionContainer, pGeometryParent)
    {
    }

    /// Destructor.
    ~QuadraturePointPartitionedGeometry() override = default;

    /// Copy constructor.
    QuadraturePointPartitionedGeometry(
        QuadraturePointPartitionedGeometry const& rOther )
        : BaseType( rOther )
        , mGeometryData(rOther.mGeometryData)
        , mpGeometryParent(rOther.mpGeometryParent)
    {
    }

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    QuadraturePointPartitionedGeometry& operator=(
        const QuadraturePointPartitionedGeometry& rOther )
    {
        BaseType::operator=( rOther );

        mGeometryData = rOther.mGeometryData;
        mpGeometryParent = rOther.mpGeometryParent;

        return *this;
    }

    ///@}
    ///@name Operations
    ///@{
    /*
    typename BaseType::Pointer Create( PointsArrayType const& ThisPoints ) const override
    {
        KRATOS_ERROR << "QuadraturePointPartitionedGeometry cannot be created with 'PointsArrayType const& ThisPoints'. "
            << "This constructor is not allowed as it would remove the evaluated shape functions as the ShapeFunctionContainer is not being copied."
            << std::endl;
    }
    */
    ///@}
    ///@name Input and output
    ///@{
    /*
    /// Turn back information as a string.
    std::string Info() const override
    {
        return "Quadrature point templated by local space dimension and working space dimension.";
    }

    /// Print information about this object.
    void PrintInfo( std::ostream& rOStream ) const override
    {
        rOStream << "Quadrature point templated by local space dimension and working space dimension.";
    }

    /// Print object's data.
    void PrintData( std::ostream& rOStream ) const override
    {
    }
    */
    ///@}

protected:

    ///@name Constructor
    ///@{

    /// Standard Constructor
    QuadraturePointPartitionedGeometry()
        : BaseType(
            PointsArrayType(),
            &mGeometryData)
        , mGeometryData(
            &msGeometryDimension,
            GeometryData::GI_GAUSS_1,
            {}, {}, {})
    {
    }

    ///@}

private:
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

    ///@}
}; // Class Geometry

///@name Input and output
///@{

/// input stream function
template<class TPointType,
    int TWorkingSpaceDimension,
    int TLocalSpaceDimension,
    int TDimension>
inline std::istream& operator >> (
    std::istream& rIStream,
    QuadraturePointPartitionedGeometry<TPointType, TWorkingSpaceDimension, TLocalSpaceDimension, TDimension>& rThis );

/// output stream function
template<class TPointType,
    int TWorkingSpaceDimension,
    int TLocalSpaceDimension,
    int TDimension>
inline std::ostream& operator << (
    std::ostream& rOStream,
    const QuadraturePointPartitionedGeometry<TPointType, TWorkingSpaceDimension, TLocalSpaceDimension, TDimension>& rThis )
{
    rThis.PrintInfo( rOStream );
    rOStream << std::endl;
    rThis.PrintData( rOStream );

    return rOStream;
}

///@}

}  // namespace Kratos.

#endif // KRATOS_QUADRATURE_POINT_PARTITIONED_GEOMETRY_H_INCLUDED  defined

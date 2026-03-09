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

    using BaseType::ShapeFunctionsValues;
    using BaseType::ShapeFunctionsLocalGradients;
    using BaseType::Jacobian;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor with points and geometry shape function container
    QuadraturePointPartitionedGeometry(
        const PointsArrayType& ThisPoints,
        GeometryShapeFunctionContainerType& ThisGeometryShapeFunctionContainer)
        : BaseType(ThisPoints, ThisGeometryShapeFunctionContainer)
    {
    }

    /// Constructor with points, geometry shape function container, parent
    QuadraturePointPartitionedGeometry(
    const PointsArrayType& ThisPoints,
    GeometryShapeFunctionContainerType& ThisGeometryShapeFunctionContainer,
    GeometryType* pGeometryParent)
    : BaseType(ThisPoints, ThisGeometryShapeFunctionContainer, pGeometryParent)
    {
    }

    /// Destructor.
    ~QuadraturePointPartitionedGeometry() override = default;

    /// Copy constructor.
    QuadraturePointPartitionedGeometry(
        QuadraturePointPartitionedGeometry const& rOther )
        : BaseType( rOther )
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

        return *this;
    }

    Matrix& Jacobian(Matrix& rResult, IndexType IntegrationPointIndex, GeometryData::IntegrationMethod ThisMethod) const override
    {
        const SizeType working_space_dimension = this->WorkingSpaceDimension();
        const SizeType local_space_dimension = this->LocalSpaceDimension();
        if (rResult.size1() != working_space_dimension || rResult.size2() != local_space_dimension)
            rResult.resize(working_space_dimension, local_space_dimension, false);

        const Matrix& r_shape_functions_gradient_in_integration_point = ShapeFunctionsLocalGradients(ThisMethod)[IntegrationPointIndex];
        const Matrix& r_N = ShapeFunctionsValues();

        rResult.clear();
        const SizeType points_number = this->PointsNumber();
        IndexType active_point_index = 0;
        for (IndexType i = 0; i < points_number; ++i) {
            if (r_N(IntegrationPointIndex,i) >= 0.0) { // only use active nodes
                const array_1d<double, 3>& r_coordinates = (*this)[i].Coordinates();
                for (IndexType k = 0; k < working_space_dimension; ++k) {
                    const double value = r_coordinates[k];
                    for (IndexType m = 0; m < local_space_dimension; ++m)
                        rResult(k, m) += value * r_shape_functions_gradient_in_integration_point(active_point_index, m);
                }
                active_point_index += 1; // increment active node counter
            }
        }
        return rResult;
    }



protected:

    ///@name Constructor
    ///@{

    /// Standard Constructor
    QuadraturePointPartitionedGeometry()
        : BaseType()
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

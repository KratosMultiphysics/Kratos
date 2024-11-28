//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Manuel Messmer
//

#if !defined(KRATOS_QUADRATURE_POINT_SURFACE_IN_VOLUME_GEOMETRY_H_INCLUDED )
#define  KRATOS_QUADRATURE_POINT_SURFACE_IN_VOLUME_GEOMETRY_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/variables.h"
#include "geometries/quadrature_point_geometry.h"

namespace Kratos
{
/**
 * @class QuadraturePointSurfaceInVolumeGeometry
 * @ingroup KratosCore
 * @brief A sinlge quadrature point, that can be used for geometries without
 *        a predefined integration scheme, i.e. they can handle material point elements,
 *        isogeometric analysis elements or standard finite elements which are defined
 *        at a single quadrature point.
 *        This point defines a surface segment described in a underlying volume.
 *        Shape functions and integration types have to be precomputed and are set from outside.
 *        The parent pointer can provide the adress to the owner of this quadrature point.
 */
template<class TPointType>
class QuadraturePointSurfaceInVolumeGeometry
    : public QuadraturePointGeometry<TPointType, 3, 3, 2>
{
public:

    /// Pointer definition of QuadraturePointGeometry
    KRATOS_CLASS_POINTER_DEFINITION(QuadraturePointSurfaceInVolumeGeometry);

    typedef QuadraturePointGeometry<TPointType, 3, 3, 2> BaseType;
    typedef Geometry<TPointType> GeometryType;

    typedef typename GeometryType::IndexType IndexType;
    typedef typename GeometryType::SizeType SizeType;

    typedef typename GeometryType::PointsArrayType PointsArrayType;
    typedef typename GeometryType::CoordinatesArrayType CoordinatesArrayType;

    typedef GeometryShapeFunctionContainer<GeometryData::IntegrationMethod> GeometryShapeFunctionContainerType;
    typedef BoundedMatrix<double,3,2> TangentMatrixType;

    /// using base class functions
    using BaseType::Jacobian;
    using BaseType::Calculate;
    using BaseType::DeterminantOfJacobian;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor with points and geometry shape function container
    QuadraturePointSurfaceInVolumeGeometry(
        const PointsArrayType& ThisPoints,
        GeometryShapeFunctionContainerType& ThisGeometryShapeFunctionContainer,
        const TangentMatrixType& LocalTangents )
        : BaseType(ThisPoints, ThisGeometryShapeFunctionContainer)
    {
        KRATOS_ERROR_IF( mLocalTangents.size1() != LocalTangents.size1() || mLocalTangents.size2() != LocalTangents.size2() )
            << "QuadraturePointSurfaceInVolumeGeometry :: Dimensions of LocalTangents do not match (3,2). "
            << "Given Dimensions are: (" << LocalTangents.size1() << "," << LocalTangents.size2() <<"). " << std::endl;
        mLocalTangents = LocalTangents;
    }

    // SBM
    QuadraturePointSurfaceInVolumeGeometry(
        const PointsArrayType& ThisPoints,
        GeometryShapeFunctionContainerType& ThisGeometryShapeFunctionContainer,
        const TangentMatrixType& LocalTangents,
        GeometryType* pGeometryParent,
        Vector& Normal)
        : BaseType(ThisPoints, ThisGeometryShapeFunctionContainer, pGeometryParent)
    {
        KRATOS_ERROR_IF( mLocalTangents.size1() != LocalTangents.size1() || mLocalTangents.size2() != LocalTangents.size2() )
            << "QuadraturePointSurfaceInVolumeGeometry :: Dimensions of LocalTangents do not match (3,2). "
            << "Given Dimensions are: (" << LocalTangents.size1() << "," << LocalTangents.size2() <<"). " << std::endl;
        mLocalTangents = LocalTangents;
        mNormal = Normal;
    }

    /// Constructor with points, geometry shape function container and parent
    QuadraturePointSurfaceInVolumeGeometry(
        const PointsArrayType& ThisPoints,
        GeometryShapeFunctionContainerType& ThisGeometryShapeFunctionContainer,
        const TangentMatrixType& LocalTangents,
        GeometryType* pGeometryParent)
        : BaseType(ThisPoints, ThisGeometryShapeFunctionContainer, pGeometryParent)
    {
        KRATOS_ERROR_IF( mLocalTangents.size1() != LocalTangents.size1() || mLocalTangents.size2() != LocalTangents.size2() )
            << "QuadraturePointSurfaceInVolumeGeometry :: Dimensions of LocalTangents do not match (3,2). "
            << "Given Dimensions are: (" << LocalTangents.size1() << "," << LocalTangents.size2() <<"). " << std::endl;
        mLocalTangents = LocalTangents;
    }

    /// Destructor.
    ~QuadraturePointSurfaceInVolumeGeometry() override = default;

    /// Copy constructor.
    QuadraturePointSurfaceInVolumeGeometry(QuadraturePointSurfaceInVolumeGeometry const& rOther )
        : BaseType( rOther )
        , mLocalTangents(rOther.mLocalTangents)
    {
    }

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    QuadraturePointSurfaceInVolumeGeometry& operator=( const QuadraturePointSurfaceInVolumeGeometry& rOther )
    {
        BaseType::operator=( rOther );

        mLocalTangents = rOther.mLocalTangents;

        return *this;
    }

    ///@}
    ///@name Dynamic access to internals
    ///@{

    /// Assign with Matrix
    void Assign(
        const Variable<Matrix >& rVariable,
        const Matrix& rInput) override
    {
        KRATOS_ERROR_IF( rInput.size1() != mLocalTangents.size1() || rInput.size2() != mLocalTangents.size2() )
            << "QuadraturePointSurfaceInVolume :: Input Matrix does not have the same size as LocalTangentMatrix. Size must be (3,2)." << std::endl;

        if (rVariable == LOCAL_TANGENT_MATRIX) {
            mLocalTangents = rInput;
        }
    }

    /// Calculate with Matrix
    void Calculate(
        const Variable<Matrix >& rVariable,
        Matrix& rOutput) const override
    {
        if (rVariable == LOCAL_TANGENT_MATRIX) {
            if( rOutput.size1() != mLocalTangents.size1() || rOutput.size2() != mLocalTangents.size2() ){
                rOutput.resize(mLocalTangents.size1(),mLocalTangents.size2());
            }
            rOutput = mLocalTangents;
        }
    } 
    
    
    /// Calculate with Matrix
    void Calculate(
        const Variable<Vector>& rVariable,
        Vector& rOutput) const override
    {
        if (rVariable == NORMAL) {
            if( rOutput.size() != mNormal.size()){
                rOutput.resize(mNormal.size());
            }

            if (mNormal.size() == 0) KRATOS_ERROR << "[QUADRATURE_POINT_IN_SURFACE_GEOMETRY]:: Normal not defined" << std::endl;
            rOutput = mNormal;
        }
    }
    
    void Calculate(
        const Variable<array_1d<double, 3>>& rVariable,
        array_1d<double, 3>& rOutput) const override
    {
        if (rVariable == NORMAL) {
            if( rOutput.size() != mNormal.size()){
                rOutput.resize(mNormal.size());
            }

            if (mNormal.size() == 0) KRATOS_ERROR << "[QUADRATURE_POINT_IN_SURFACE_GEOMETRY]:: Normal not defined" << std::endl;
            rOutput = mNormal;
        }
    }

    ///@}
    ///@name Normal
    ///@{

    /**
    * @brief Computes the normal of the surface
    *        lying on the underlying volume.
    * @param rResult Normal to the surface lying in the volume in global coordinates.
    * @param IntegrationPointIndex Index should be always 0.
    **/
    CoordinatesArrayType Normal(
        IndexType IntegrationPointIndex,
        GeometryData::IntegrationMethod ThisMethod) const override
    {
        KRATOS_DEBUG_ERROR_IF(IntegrationPointIndex != 0)
            << "Trying to access Normal of QuadraturePointCurveOnSurface "
            << "with an integration point index != 0." << std::endl;

        Matrix J;
        this->Jacobian(J, IntegrationPointIndex, ThisMethod);

        // Map tangents into global space
        TangentMatrixType global_tangents = prod(J, mLocalTangents);

        Vector g1 = column(global_tangents,0);
        Vector g2 = column(global_tangents,1);

        // Compute normal
        CoordinatesArrayType normal = MathUtils<double>::CrossProduct( g1, g2 );

        return normal;
    }

    ///@}
    ///@name Jacobian
    ///@{

    /**
    * @brief Returns the respective surface area of this
    *        quadrature point in global space.
    * @param IntegrationPointIndex Index should be always 0.
    **/
    double DeterminantOfJacobian(
        IndexType IntegrationPointIndex,
        GeometryData::IntegrationMethod ThisMethod) const override
    {
        KRATOS_DEBUG_ERROR_IF(IntegrationPointIndex != 0)
            << "Trying to access DeterminantOfJacobian of QuadraturePointSurfaceInVolume "
            << "with an integration point index != 0." << std::endl;

        Matrix J;
        this->Jacobian(J, IntegrationPointIndex, ThisMethod);
        TangentMatrixType global_tangents = prod(J, mLocalTangents);

        return MathUtils<double>::GeneralizedDet( global_tangents );
    }

    /**
     * @brief Returns the determinant of the jacobian at this
     *        quadrature point.
     * @param rResult Vector of results.
     */
    Vector& DeterminantOfJacobian(
        Vector& rResult,
        GeometryData::IntegrationMethod ThisMethod) const override
    {
        if (rResult.size() != 1)
            rResult.resize(1, false);

        rResult[0] = this->DeterminantOfJacobian(0, ThisMethod);

        return rResult;
    }

    ///@}
    ///@name Kratos Geometry Families
    ///@{

    /**
     * @brief Gets the geometry family.
     * @details This function returns the family type of the geometry. The geometry family categorizes the geometry into a broader classification, aiding in its identification and processing.
     * @return GeometryData::KratosGeometryFamily The geometry family.
     */
    GeometryData::KratosGeometryFamily GetGeometryFamily() const override
    {
        return GeometryData::KratosGeometryFamily::Kratos_Quadrature_Geometry;
    }

    /**
     * @brief Gets the geometry type.
     * @details This function returns the specific type of the geometry. The geometry type provides a more detailed classification of the geometry.
     * @return GeometryData::KratosGeometryType The specific geometry type.
     */
    GeometryData::KratosGeometryType GetGeometryType() const override
    {
        return GeometryData::KratosGeometryType::Kratos_Quadrature_Point_Surface_In_Volume_Geometry;
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "Quadrature point geometry for a surface in a nurbs volume with Id: "
            << std::to_string(this->Id()) << ", containing: " << std::to_string(this->size()) << " points." << std::endl;

        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo( std::ostream& rOStream ) const override
    {
        rOStream << "Quadrature point geometry for a surface in a nurbs volume with Id: "
            << this->Id() << ", containing: " << this->size() << " points." << std::endl;
    }

    /// Print object's data.
    void PrintData( std::ostream& rOStream ) const override
    {
    }
    ///@}

private:
    ///@name Member Variables
    ///@{

    TangentMatrixType mLocalTangents;

    Vector mNormal;

    ///@}
    ///@name Serialization
    ///@{

    /// Default constructor for serializer
    QuadraturePointSurfaceInVolumeGeometry()
        : BaseType()
    {
    }

    friend class Serializer;

    void save( Serializer& rSerializer ) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseType);
        rSerializer.save("LocalTangents", mLocalTangents);
    }

    void load( Serializer& rSerializer ) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BaseType );
        rSerializer.load("LocalTangents", mLocalTangents);
    }

    ///@}
}; // Class Geometry

}  // namespace Kratos.

#endif // KRATOS_QUADRATURE_POINT_SURFACE_IN_VOLUME_GEOMETRY_H_INCLUDED  defined

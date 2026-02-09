//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ricky Aristio
//

#pragma once

// System includes

// External includes

// Project includes
#include "includes/variables.h"
#include "geometries/geometry.h"
#include "geometries/geometry_dimension.h"
#include "utilities/integration_utilities.h"

namespace Kratos
{
/**
 * @class KnotSpanGeometry
 * @ingroup KratosCore
 * @brief Multiple quadrature points obtained by non-zero knot span, that can be used for geometries without
 *        a predefined integration scheme, i.e. they can handle isogeometric analysis elements.
 *        Shape functions and integration types are precomputed and are set from
 *        from outside.
 *        The parent pointer can provide the adress to the owner of this quadrature point.
 */
template<class TPointType,
    int TWorkingSpaceDimension,
    int TLocalSpaceDimension = TWorkingSpaceDimension,
    int TDimension = TLocalSpaceDimension>
class KnotSpanGeometry
    : public QuadraturePointGeometry<TPointType, 3, 2, 2>
{
public:

    /// Pointer definition of KnotSpanGeometry
    KRATOS_CLASS_POINTER_DEFINITION( KnotSpanGeometry );

    typedef QuadraturePointGeometry<TPointType, 3, 2, 2> BaseType;
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

    /// using base class functions
    using BaseType::Jacobian;
    using BaseType::DeterminantOfJacobian;
    using BaseType::ShapeFunctionsValues;
    using BaseType::ShapeFunctionsLocalGradients;
    using BaseType::InverseOfJacobian;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor with points and all shape function containers separately
    // KnotSpanGeometry(
    //     const PointsArrayType& ThisPoints,
    //     const IntegrationPointsContainerType& rIntegrationPoints,
    //     const ShapeFunctionsValuesContainerType& rShapeFunctionValues,
    //     const ShapeFunctionsLocalGradientsContainerType& rShapeFunctionsDerivativesVector,
    //     const std::size_t PointsInU,
    //     const std::size_t PointsInV,
    //     const double KnotSpanIntervalUBegin,
    //     const double KnotSpanIntervalUEnd,
    //     const double KnotSpanIntervalVBegin,
    //     const double KnotSpanIntervalVEnd)
    //     : BaseType(ThisPoints, rIntegrationPoints, rShapeFunctionValues, rShapeFunctionsDerivativesVector)
    //     , mPointsInU(PointsInU)
    //     , mPointsInV(PointsInV)
    //     , mKnotSpanIntervalUBegin(KnotSpanIntervalUBegin)
    //     , mKnotSpanIntervalUEnd(KnotSpanIntervalUEnd)
    //     , mKnotSpanIntervalVBegin(KnotSpanIntervalVBegin)
    //     , mKnotSpanIntervalVEnd(KnotSpanIntervalVEnd)
    // {
    // }

    /// Constructor with points and all shape function containers separately including the parent
    // KnotSpanGeometry(
    //     const PointsArrayType& ThisPoints,
    //     const IntegrationPointsContainerType& rIntegrationPoints,
    //     const ShapeFunctionsValuesContainerType& rShapeFunctionValues,
    //     const ShapeFunctionsLocalGradientsContainerType& rShapeFunctionsDerivativesVector,
    //     GeometryType* pGeometryParent,
    //     const std::size_t PointsInU,
    //     const std::size_t PointsInV,
    //     const double KnotSpanIntervalUBegin,
    //     const double KnotSpanIntervalUEnd,
    //     const double KnotSpanIntervalVBegin,
    //     const double KnotSpanIntervalVEnd)
    //     : BaseType(ThisPoints, rIntegrationPoints, rShapeFunctionValues, rShapeFunctionsDerivativesVector, pGeometryParent)
    //     , mPointsInU(PointsInU)
    //     , mPointsInV(PointsInV)
    //     , mKnotSpanIntervalUBegin(KnotSpanIntervalUBegin)
    //     , mKnotSpanIntervalUEnd(KnotSpanIntervalUEnd)
    //     , mKnotSpanIntervalVBegin(KnotSpanIntervalVBegin)
    //     , mKnotSpanIntervalVEnd(KnotSpanIntervalVEnd)
    // {
    // }

    /// Constructor with points and geometry shape function container
    KnotSpanGeometry(
        const PointsArrayType& ThisPoints,
        const GeometryShapeFunctionContainerType& ThisGeometryShapeFunctionContainer,
        const std::size_t PointsInU,
        const std::size_t PointsInV,
        const double KnotSpanIntervalUBegin,
        const double KnotSpanIntervalUEnd,
        const double KnotSpanIntervalVBegin,
        const double KnotSpanIntervalVEnd)
        : BaseType(ThisPoints, ThisGeometryShapeFunctionContainer)
        , mPointsInU(PointsInU)
        , mPointsInV(PointsInV)
        , mKnotSpanIntervalUBegin(KnotSpanIntervalUBegin)
        , mKnotSpanIntervalUEnd(KnotSpanIntervalUEnd)
        , mKnotSpanIntervalVBegin(KnotSpanIntervalVBegin)
        , mKnotSpanIntervalVEnd(KnotSpanIntervalVEnd)
    {
    }

    /// Constructor with points, geometry shape function container, parent
    KnotSpanGeometry(
        const PointsArrayType& ThisPoints,
        const GeometryShapeFunctionContainerType& ThisGeometryShapeFunctionContainer,
        GeometryType* pGeometryParent,
        const std::size_t PointsInU,
        const std::size_t PointsInV,
        const double KnotSpanIntervalUBegin,
        const double KnotSpanIntervalUEnd,
        const double KnotSpanIntervalVBegin,
        const double KnotSpanIntervalVEnd)
        : BaseType(ThisPoints, ThisGeometryShapeFunctionContainer, pGeometryParent)
        , mPointsInU(PointsInU)
        , mPointsInV(PointsInV)
        , mKnotSpanIntervalUBegin(KnotSpanIntervalUBegin)
        , mKnotSpanIntervalUEnd(KnotSpanIntervalUEnd)
        , mKnotSpanIntervalVBegin(KnotSpanIntervalVBegin)
        , mKnotSpanIntervalVEnd(KnotSpanIntervalVEnd)
    {
    }

    /// Constructor with points, N, Vector<DN_De, ...>
    // KnotSpanGeometry(
    //     const PointsArrayType& ThisPoints,
    //     const IntegrationPointType& ThisIntegrationPoint,
    //     const Matrix& ThisShapeFunctionsValues,
    //     const DenseVector<Matrix>& ThisShapeFunctionsDerivatives)
    //     : BaseType(ThisPoints, &mGeometryData)
    //     , mGeometryData(
    //         &msGeometryDimension,
    //         GeometryShapeFunctionContainerType(
    //             GeometryData::IntegrationMethod::GI_GAUSS_1,
    //             ThisIntegrationPoint,
    //             ThisShapeFunctionsValues,
    //             ThisShapeFunctionsDerivatives))
    // {
    // }

    /// Constructor with points, N, Vector<DN_De, ...>, parent
    // KnotSpanGeometry(
    //     const PointsArrayType& ThisPoints,
    //     const IntegrationPointType& ThisIntegrationPoint,
    //     const Matrix& ThisShapeFunctionsValues,
    //     const DenseVector<Matrix>& ThisShapeFunctionsDerivatives,
    //     GeometryType* pGeometryParent)
    //     : BaseType(ThisPoints, &mGeometryData)
    //     , mGeometryData(
    //         &msGeometryDimension,
    //         GeometryShapeFunctionContainerType(
    //             GeometryData::IntegrationMethod::GI_GAUSS_1,
    //             ThisIntegrationPoint,
    //             ThisShapeFunctionsValues,
    //             ThisShapeFunctionsDerivatives))
    //     , mpGeometryParent(pGeometryParent)
    // {
    // }

    // /// Constructor.
    // KnotSpanGeometry(
    //     const PointsArrayType& ThisPoints) = delete;

    // /// Constructor with Geometry Id
    // KnotSpanGeometry(
    //     const IndexType GeometryId,
    //     const PointsArrayType& ThisPoints
    // ) : BaseType( GeometryId, ThisPoints, &mGeometryData )
    //     , mGeometryData(
    //         &msGeometryDimension,
    //         GeometryData::IntegrationMethod::GI_GAUSS_1,
    //         {}, {}, {})
    // {
    // }

    // /// Constructor with Geometry Name
    // KnotSpanGeometry(
    //     const std::string& GeometryName,
    //     const PointsArrayType& ThisPoints
    // ) : BaseType( GeometryName, ThisPoints, &mGeometryData )
    //     , mGeometryData(
    //         &msGeometryDimension,
    //         GeometryData::IntegrationMethod::GI_GAUSS_1,
    //         {}, {}, {})
    // {
    // }

    /// Destructor.
    ~KnotSpanGeometry() override = default;

    /// Copy constructor.
    KnotSpanGeometry(
        KnotSpanGeometry const& rOther )
        : BaseType( rOther )
        , mPointsInU(rOther.PointsInU)
        , mPointsInV(rOther.PointsInV)
        , mKnotSpanIntervalUBegin(rOther.KnotSpanIntervalUBegin)
        , mKnotSpanIntervalUEnd(rOther.KnotSpanIntervalUEnd)
        , mKnotSpanIntervalVBegin(rOther.KnotSpanIntervalVBegin)
        , mKnotSpanIntervalVEnd(rOther.KnotSpanIntervalVEnd)
    {
    }

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    KnotSpanGeometry& operator=(
        const KnotSpanGeometry& rOther )
    {
        BaseType::operator=( rOther );

        mKnotSpanIntervalUBegin = rOther.KnotSpanIntervalUBegin;
        mKnotSpanIntervalUEnd = rOther.KnotSpanIntervalUEnd;
        mKnotSpanIntervalVBegin = rOther.KnotSpanIntervalVBegin;
        mKnotSpanIntervalVEnd = rOther.KnotSpanIntervalVEnd;

        return *this;
    }

    ///@}
    ///@name Operations
    ///@{

    // /**
    //  * @brief Creates a new geometry pointer
    //  * @param ThisPoints the nodes of the new geometry
    //  * @return Pointer to the new geometry
    //  */
    // typename BaseType::Pointer Create( PointsArrayType const& ThisPoints ) const override
    // {
    //     KRATOS_ERROR << "KnotSpanGeometry cannot be created with 'PointsArrayType const& PointsArrayType'. "
    //         << "This constructor is not allowed as it would remove the evaluated shape functions as the ShapeFunctionContainer is not being copied."
    //         << std::endl;
    // }

    // /**
    //  * @brief Creates a new geometry pointer
    //  * @param NewGeometryId the ID of the new geometry
    //  * @param rThisPoints the nodes of the new geometry
    //  * @return Pointer to the new geometry
    //  */
    // typename BaseType::Pointer Create(
    //     const IndexType NewGeometryId,
    //     PointsArrayType const& rThisPoints
    //     ) const override
    // {
    //     return typename BaseType::Pointer( new KnotSpanGeometry( NewGeometryId, rThisPoints ) );
    // }

    // /**
    //  * @brief Creates a new geometry pointer
    //  * @param NewGeometryId the ID of the new geometry
    //  * @param rGeometry reference to an existing geometry
    //  * @return Pointer to the new geometry
    //  */
    // typename BaseType::Pointer Create(
    //     const IndexType NewGeometryId,
    //     const BaseType& rGeometry
    // ) const override
    // {
    //     auto p_geometry = typename BaseType::Pointer( new KnotSpanGeometry( NewGeometryId, rGeometry.Points() ) );
    //     p_geometry->SetData(rGeometry.GetData());
    //     return p_geometry;
    // }

    ///@}
    ///@name Dynamic access to internals
    ///@{

    /// Calculate with array_1d<double, 3>
    // void Calculate(
    //     const Variable<array_1d<double, 3>>& rVariable,
    //     array_1d<double, 3>& rOutput) const override
    // {
    //     if (rVariable == CHARACTERISTIC_GEOMETRY_LENGTH)
    //     {
    //         rOutput = this->IntegrationPoints()[0];
    //         mpGeometryParent->Calculate(rVariable, rOutput);
    //     }
    // }

    // ///@}
    // ///@name  Geometry Shape Function Container
    // ///@{

    // /* @brief SetGeometryShapeFunctionContainer updates the GeometryShapeFunctionContainer within
    //  *        the GeometryData. This function works only for geometries with a non-const GeometryData.
    //  */
    // void SetGeometryShapeFunctionContainer(
    //     const GeometryShapeFunctionContainer<GeometryData::IntegrationMethod>& rGeometryShapeFunctionContainer) override
    // {
    //     mGeometryData.SetGeometryShapeFunctionContainer(rGeometryShapeFunctionContainer);
    // }

    // ///@}
    // ///@name Parent
    // ///@{

    // GeometryType& GetGeometryParent(IndexType Index) const override
    // {
    //     return *mpGeometryParent;
    // }

    // void SetGeometryParent(GeometryType* pGeometryParent) override
    // {
    //     mpGeometryParent = pGeometryParent;
    // }

    ///@}
    ///@name Dynamic access to internals
    ///@{

    /// Calculate with Vector
    void Calculate(
        const Variable<Vector>& rVariable,
        Vector& rOutput) const override
    {
        // if (rVariable == DETERMINANTS_OF_JACOBIAN_PARENT) {
        //     DeterminantOfJacobianParent(rOutput);
        // }
    }

    ///@}
    ///@name Geometrical Operations
    ///@{

    /// Returns the domain size of this quadrature point.
    double DomainSize() const override
    {
        return IntegrationUtilities::ComputeDomainSize(*this);
    }

    /**
    * @brief Calculates global location of this integration point.
    * \f[
    * c_i = \sum_j^N(x_j)*x_i
    * \f]
    * j is the index of the node and i the global direction (x,y,z).
    *
    * @return Point which is the location of this quadrature point.
    */
    Point Center() const override
    {
        const std::size_t node_number = this->PointsNumber();

        Point point(0.0, 0.0, 0.0);
        const Matrix& r_N = this->ShapeFunctionsValues();

        for (IndexType point_number = 0; point_number < this->IntegrationPointsNumber(); ++point_number) {
            for (IndexType i = 0; i < node_number; ++i) {
                point += (*this)[i] * r_N(point_number, i);
            }
        }
        return point;
    }

    ///@}
    ///@name Mathematical Informations
    ///@{

    // /// Returns the polynomial degree of the parent geometry
    // SizeType PolynomialDegree(IndexType LocalDirectionIndex) const override
    // {
    //     KRATOS_DEBUG_ERROR_IF_NOT(mpGeometryParent)
    //         << "Trying to call PolynomialDegree(LocalDirectionIndex) from quadrature point. "
    //         << "Pointer to parent is not assigned." << std::endl;

    //     return mpGeometryParent->PolynomialDegree(LocalDirectionIndex);
    // }

    ///@}
    ///@name Coordinates
    ///@{

    /**
    * This method provides the global coordinates corresponding to the
    * local coordinates. The mapping is done on the parent.
    * Error if Parent is not assigned is only thrown in debug mode.
    */
    // CoordinatesArrayType& GlobalCoordinates(
    //     CoordinatesArrayType& rResult,
    //     CoordinatesArrayType const& LocalCoordinates
    // ) const override
    // {
    //     KRATOS_DEBUG_ERROR_IF(mpGeometryParent == nullptr)
    //         << "Trying to call GlobalCoordinates(LocalCoordinates) from quadrature point. "
    //         << "Pointer to parent is not assigned." << std::endl;

    //     return mpGeometryParent->GlobalCoordinates(rResult, LocalCoordinates);
    // }

    /* @brief returns the respective segment domain size of this
     *        quadrature point, computed on the parent of this geometry.
     *        Required for reduced quadrature point geometries (Not all
     *        nodes are part of this geometry - used for mapping).
     * @param rResult vector of results of this quadrature point.
     */
    // CoordinatesArrayType& GlobalCoordinates(
    //     CoordinatesArrayType& rResult,
    //     CoordinatesArrayType const& LocalCoordinates,
    //     Matrix& DeltaPosition
    // ) const override
    // {
    //     KRATOS_DEBUG_ERROR_IF(mpGeometryParent == nullptr)
    //         << "Trying to call GlobalCoordinates(LocalCoordinates, DeltaPosition) from quadrature point. "
    //         << "Pointer to parent is not assigned." << std::endl;

    //     return mpGeometryParent->GlobalCoordinates(rResult, LocalCoordinates, DeltaPosition);
    // }

    ///@}
    ///@name Jacobian
    ///@{

    /**
    * @brief Jacobian in given point. Computed on parent
    *        geometry.
    * Error if Parent is not assigned is only thrown in debug mode.
    * @param rCoordinates point which jacobians has to
    *        be calculated in it.
    * @return jacobian matrix \f$ J \f$ in given point.
    */
    // Matrix& Jacobian(
    //     Matrix& rResult,
    //     const CoordinatesArrayType& rCoordinates
    // ) const override
    // {
    //     KRATOS_DEBUG_ERROR_IF(mpGeometryParent == nullptr)
    //         << "Trying to call Jacobian(LocalCoordinates) from quadrature point. "
    //         << "Pointer to parent is not assigned." << std::endl;

    //     return mpGeometryParent->Jacobian(rResult, rCoordinates);
    // }

    /** Determinant of jacobian in given point. Computed on parent
    *        geometry.
    * Error if Parent is not assigned is only thrown in debug mode.
    * @param rPoint point which determinant of jacobians has to
    *        be calculated in it.
    * @return Determinamt of jacobian matrix \f$ |J| \f$ in given
    *         point.
    */
    // double DeterminantOfJacobian(
    //     const CoordinatesArrayType& rPoint
    // ) const override
    // {
    //     KRATOS_DEBUG_ERROR_IF(mpGeometryParent == nullptr)
    //         << "Trying to call DeterminantOfJacobian(rPoint) from quadrature point. "
    //         << "Pointer to parent is not assigned." << std::endl;

    //     return mpGeometryParent->DeterminantOfJacobian(rPoint);
    // }


    /* @brief returns the respective segment length of this
     *        quadrature point. Length of vector always 1.
     * @param rResult vector of results of this quadrature point.
     */
    // Vector& DeterminantOfJacobianParent(
    //     Vector& rResult) const
    // {
    //     if (rResult.size() != 1)
    //         rResult.resize(1, false);

    //     rResult[0] = this->GetGeometryParent(0).DeterminantOfJacobian(this->IntegrationPoints()[0]);

    //     return rResult;
    // }

    // Matrix& InverseOfJacobian(
    //     Matrix& rResult,
    //     const CoordinatesArrayType& rCoordinates
    // ) const override
    // {
    //     KRATOS_DEBUG_ERROR_IF(mpGeometryParent == nullptr)
    //         << "Trying to call InverseOfJacobian(rPoint) from quadrature point. "
    //         << "Pointer to parent is not assigned." << std::endl;

    //     return mpGeometryParent->InverseOfJacobian(rResult, rCoordinates);
    // }

    ///@}
    ///@name Shape Function
    ///@{

    // Vector& ShapeFunctionsValues(
    //     Vector &rResult,
    //     const CoordinatesArrayType& rCoordinates
    // ) const override
    // {
    //     KRATOS_DEBUG_ERROR_IF(mpGeometryParent == nullptr)
    //         << "Trying to call ShapeFunctionsValues(rCoordinates) from quadrature point. "
    //         << "Pointer to parent is not assigned." << std::endl;

    //     return mpGeometryParent->ShapeFunctionsValues(rResult, rCoordinates);
    // }

    // virtual Matrix& ShapeFunctionsLocalGradients(
    //     Matrix& rResult,
    //     const CoordinatesArrayType& rPoint
    // ) const override
    // {
    //     KRATOS_DEBUG_ERROR_IF(mpGeometryParent == nullptr)
    //         << "Trying to call ShapeFunctionsLocalGradients(rPoint) from quadrature point. "
    //         << "Pointer to parent is not assigned." << std::endl;

    //     return mpGeometryParent->ShapeFunctionsLocalGradients(rResult, rPoint);
    // }

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
        return GeometryData::KratosGeometryType::Kratos_Knot_Span_Geometry;
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "Knot span templated by local space dimension and working space dimension.";
    }

    /// Print information about this object.
    void PrintInfo( std::ostream& rOStream ) const override
    {
        rOStream << "Knot span templated by local space dimension and working space dimension.";
    }

    /// Print object's data.
    void PrintData( std::ostream& rOStream ) const override
    {
    }
    ///@}

protected:

    ///@name Constructor
    ///@{

    // /// Standard Constructor
    // KnotSpanGeometry()
    //     : BaseType(
    //         PointsArrayType(),
    //         &mGeometryData)
    //     , mGeometryData(
    //         &msGeometryDimension,
    //         GeometryData::IntegrationMethod::GI_GAUSS_1,
    //         {}, {}, {})
    // {
    // }

    ///@}

private:
    ///@name Static Member Variables
    ///@{

    std::size_t mPointsInU;
    std::size_t mPointsInV; 
    double mKnotSpanIntervalUBegin;
    double mKnotSpanIntervalUEnd;
    double mKnotSpanIntervalVBegin;
    double mKnotSpanIntervalVEnd;

    // static const GeometryDimension msGeometryDimension;

    // ///@}
    // ///@name Member Variables
    // ///@{

    // GeometryData mGeometryData;

    // // quatrature point can be related to a parent geometry. To keep the connection,
    // // this geometry is related to the integration point.
    // GeometryType* mpGeometryParent = nullptr;

    ///@}
    ///@name Serialization
    ///@{

    KnotSpanGeometry()
        : BaseType()
    {
    }

    friend class Serializer;

    void save( Serializer& rSerializer ) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, BaseType );
        rSerializer.save("mPointsInU", mPointsInU);
        rSerializer.save("mPointsInV", mPointsInV);
        rSerializer.save("KnotSpanIntervalUBegin", mKnotSpanIntervalUBegin);
        rSerializer.save("KnotSpanIntervalUEnd", mKnotSpanIntervalUEnd);
        rSerializer.save("KnotSpanIntervalVBegin", mKnotSpanIntervalVBegin);
        rSerializer.save("KnotSpanIntervalVEnd", mKnotSpanIntervalVEnd);

        // rSerializer.save("IntegrationPoints", mGeometryData.IntegrationPoints());
        // rSerializer.save("ShapeFunctionsValues", mGeometryData.ShapeFunctionsValues());
        // rSerializer.save("ShapeFunctionsLocalGradients", mGeometryData.ShapeFunctionsLocalGradients());
    }

    void load( Serializer& rSerializer ) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BaseType );
        rSerializer.load("mPointsInU", mPointsInU);
        rSerializer.load("mPointsInV", mPointsInV);
        rSerializer.load("KnotSpanIntervalUBegin", mKnotSpanIntervalUBegin);
        rSerializer.load("KnotSpanIntervalUEnd", mKnotSpanIntervalUEnd);
        rSerializer.load("KnotSpanIntervalVBegin", mKnotSpanIntervalVBegin);
        rSerializer.load("KnotSpanIntervalVEnd", mKnotSpanIntervalVEnd);

        // IntegrationPointsContainerType integration_points;
        // ShapeFunctionsValuesContainerType shape_functions_values;
        // ShapeFunctionsLocalGradientsContainerType shape_functions_local_gradients;

        // rSerializer.load("IntegrationPoints", integration_points[static_cast<int>(GeometryData::IntegrationMethod::GI_GAUSS_1)]);
        // rSerializer.load("ShapeFunctionsValues", shape_functions_values[static_cast<int>(GeometryData::IntegrationMethod::GI_GAUSS_1)]);
        // rSerializer.load("ShapeFunctionsLocalGradients", shape_functions_local_gradients[static_cast<int>(GeometryData::IntegrationMethod::GI_GAUSS_1)]);

        // mGeometryData.SetGeometryShapeFunctionContainer(GeometryShapeFunctionContainer<GeometryData::IntegrationMethod>(
        //     GeometryData::IntegrationMethod::GI_GAUSS_1,
        //     integration_points,
        //     shape_functions_values,
        //     shape_functions_local_gradients));

    }

    ///@}
}; // Class Geometry

///@name Input and output
///@{

/// input stream function
// template<class TPointType,
//     int TWorkingSpaceDimension,
//     int TLocalSpaceDimension,
//     int TDimension>
// inline std::istream& operator >> (
//     std::istream& rIStream,
//     KnotSpanGeometry<TPointType, TWorkingSpaceDimension, TLocalSpaceDimension, TDimension>& rThis );

/// output stream function
// template<class TPointType,
//     int TWorkingSpaceDimension,
//     int TLocalSpaceDimension,
//     int TDimension>
// inline std::ostream& operator << (
//     std::ostream& rOStream,
//     const KnotSpanGeometry<TPointType, TWorkingSpaceDimension, TLocalSpaceDimension, TDimension>& rThis )
// {
//     rThis.PrintInfo( rOStream );
//     rOStream << std::endl;
//     rThis.PrintData( rOStream );

//     return rOStream;
// }

///@}
///@name Type Dimension Definition
///@{

// template<class TPointType,
//     int TWorkingSpaceDimension,
//     int TLocalSpaceDimension,
//     int TDimension>
// const GeometryDimension KnotSpanGeometry<
//     TPointType,
//     TWorkingSpaceDimension,
//     TLocalSpaceDimension,
//     TDimension>::msGeometryDimension(
//         TWorkingSpaceDimension,
//         TLocalSpaceDimension);

///@}

}  // namespace Kratos.

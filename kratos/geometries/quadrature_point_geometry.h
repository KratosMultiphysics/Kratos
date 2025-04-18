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
 * @class QuadraturePointGeometry
 * @ingroup KratosCore
 * @brief A single quadrature point, that can be used for geometries without
 *        a predefined integration scheme, i.e. they can handle material point elements,
 *        isogeometric analysis elements or standard finite elements which are defined
 *        at a single quadrature point.
 *        Shape functions and integration types are precomputed and are set from
 *        from outside.
 *        The parent pointer can provide the address to the owner of this quadrature point.
 */
template<class TPointType,
    int TWorkingSpaceDimension,
    int TLocalSpaceDimension = TWorkingSpaceDimension,
    int TDimension = TLocalSpaceDimension>
class QuadraturePointGeometry
    : public Geometry<TPointType>
{
public:

    /// Pointer definition of QuadraturePointGeometry
    KRATOS_CLASS_POINTER_DEFINITION( QuadraturePointGeometry );

    typedef Geometry<TPointType> BaseType;
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
    QuadraturePointGeometry(
        const PointsArrayType& ThisPoints,
        const IntegrationPointsContainerType& rIntegrationPoints,
        const ShapeFunctionsValuesContainerType& rShapeFunctionValues,
        const ShapeFunctionsLocalGradientsContainerType& rShapeFunctionsDerivativesVector)
        : BaseType(ThisPoints, &mGeometryData)
        , mGeometryData(
            &msGeometryDimension,
            GeometryData::IntegrationMethod::GI_GAUSS_1,
            rIntegrationPoints,
            rShapeFunctionValues,
            rShapeFunctionsDerivativesVector)
    {
    }

    /// Constructor with points and all shape function containers separately including the parent
    QuadraturePointGeometry(
        const PointsArrayType& ThisPoints,
        const IntegrationPointsContainerType& rIntegrationPoints,
        const ShapeFunctionsValuesContainerType& rShapeFunctionValues,
        const ShapeFunctionsLocalGradientsContainerType& rShapeFunctionsDerivativesVector,
        GeometryType* pGeometryParent)
        : BaseType(ThisPoints, &mGeometryData)
        , mGeometryData(
            &msGeometryDimension,
            GeometryData::IntegrationMethod::GI_GAUSS_1,
            rIntegrationPoints,
            rShapeFunctionValues,
            rShapeFunctionsDerivativesVector)
        , mpGeometryParent(pGeometryParent)
    {
    }

    /// Constructor with points and geometry shape function container
    QuadraturePointGeometry(
        const PointsArrayType& ThisPoints,
        const GeometryShapeFunctionContainerType& ThisGeometryShapeFunctionContainer)
        : BaseType(ThisPoints, &mGeometryData)
        , mGeometryData(
            &msGeometryDimension,
            ThisGeometryShapeFunctionContainer)
    {
    }

    /// Constructor with points, geometry shape function container, parent
    QuadraturePointGeometry(
        const PointsArrayType& ThisPoints,
        const GeometryShapeFunctionContainerType& ThisGeometryShapeFunctionContainer,
        GeometryType* pGeometryParent)
        : BaseType(ThisPoints, &mGeometryData)
        , mGeometryData(
            &msGeometryDimension,
            ThisGeometryShapeFunctionContainer)
        , mpGeometryParent(pGeometryParent)
    {
    }

    /// Constructor with points, N, Vector<DN_De, ...>
    QuadraturePointGeometry(
        const PointsArrayType& ThisPoints,
        const IntegrationPointType& ThisIntegrationPoint,
        const Matrix& ThisShapeFunctionsValues,
        const DenseVector<Matrix>& ThisShapeFunctionsDerivatives)
        : BaseType(ThisPoints, &mGeometryData)
        , mGeometryData(
            &msGeometryDimension,
            GeometryShapeFunctionContainerType(
                GeometryData::IntegrationMethod::GI_GAUSS_1,
                ThisIntegrationPoint,
                ThisShapeFunctionsValues,
                ThisShapeFunctionsDerivatives))
    {
    }

    /// Constructor with points, N, Vector<DN_De, ...>, parent
    QuadraturePointGeometry(
        const PointsArrayType& ThisPoints,
        const IntegrationPointType& ThisIntegrationPoint,
        const Matrix& ThisShapeFunctionsValues,
        const DenseVector<Matrix>& ThisShapeFunctionsDerivatives,
        GeometryType* pGeometryParent)
        : BaseType(ThisPoints, &mGeometryData)
        , mGeometryData(
            &msGeometryDimension,
            GeometryShapeFunctionContainerType(
                GeometryData::IntegrationMethod::GI_GAUSS_1,
                ThisIntegrationPoint,
                ThisShapeFunctionsValues,
                ThisShapeFunctionsDerivatives))
        , mpGeometryParent(pGeometryParent)
    {
    }

    /// Constructor.
    QuadraturePointGeometry(
        const PointsArrayType& ThisPoints) = delete;

    /// Constructor with Geometry Id
    QuadraturePointGeometry(
        const IndexType GeometryId,
        const PointsArrayType& ThisPoints
    ) : BaseType( GeometryId, ThisPoints, &mGeometryData )
        , mGeometryData(
            &msGeometryDimension,
            GeometryData::IntegrationMethod::GI_GAUSS_1,
            {}, {}, {})
    {
    }

    /// Constructor with Geometry Name
    QuadraturePointGeometry(
        const std::string& GeometryName,
        const PointsArrayType& ThisPoints
    ) : BaseType( GeometryName, ThisPoints, &mGeometryData )
        , mGeometryData(
            &msGeometryDimension,
            GeometryData::IntegrationMethod::GI_GAUSS_1,
            {}, {}, {})
    {
    }

    /// Destructor.
    ~QuadraturePointGeometry() override = default;

    /// Copy constructor.
    QuadraturePointGeometry(
        QuadraturePointGeometry const& rOther )
        : BaseType( rOther )
        , mGeometryData(rOther.mGeometryData)
        , mpGeometryParent(rOther.mpGeometryParent)
    {
    }

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    QuadraturePointGeometry& operator=(
        const QuadraturePointGeometry& rOther )
    {
        BaseType::operator=( rOther );

        mGeometryData = rOther.mGeometryData;
        mpGeometryParent = rOther.mpGeometryParent;

        return *this;
    }

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Creates a new geometry pointer
     * @param ThisPoints the nodes of the new geometry
     * @return Pointer to the new geometry
     */
    typename BaseType::Pointer Create( PointsArrayType const& ThisPoints ) const override
    {
        KRATOS_ERROR << "QuadraturePointGeometry cannot be created with 'PointsArrayType const& PointsArrayType'. "
            << "This constructor is not allowed as it would remove the evaluated shape functions as the ShapeFunctionContainer is not being copied."
            << std::endl;
    }

    /**
     * @brief Creates a new geometry pointer
     * @param NewGeometryId the ID of the new geometry
     * @param rThisPoints the nodes of the new geometry
     * @return Pointer to the new geometry
     */
    typename BaseType::Pointer Create(
        const IndexType NewGeometryId,
        PointsArrayType const& rThisPoints
        ) const override
    {
        return typename BaseType::Pointer( new QuadraturePointGeometry( NewGeometryId, rThisPoints ) );
    }

    /**
     * @brief Creates a new geometry pointer
     * @param NewGeometryId the ID of the new geometry
     * @param rGeometry reference to an existing geometry
     * @return Pointer to the new geometry
     */
    typename BaseType::Pointer Create(
        const IndexType NewGeometryId,
        const BaseType& rGeometry
    ) const override
    {
        auto p_geometry = typename BaseType::Pointer( new QuadraturePointGeometry( NewGeometryId, rGeometry.Points() ) );
        p_geometry->SetData(rGeometry.GetData());
        return p_geometry;
    }

    ///@}
    ///@name Dynamic access to internals
    ///@{

    /// Calculate with array_1d<double, 3>
    void Calculate(
        const Variable<array_1d<double, 3>>& rVariable,
        array_1d<double, 3>& rOutput) const override
    {
        if (rVariable == CHARACTERISTIC_GEOMETRY_LENGTH)
        {
            rOutput = this->IntegrationPoints()[0];
            mpGeometryParent->Calculate(rVariable, rOutput);
        }
    }

    ///@}
    ///@name  Geometry Shape Function Container
    ///@{

    /* @brief SetGeometryShapeFunctionContainer updates the GeometryShapeFunctionContainer within
     *        the GeometryData. This function works only for geometries with a non-const GeometryData.
     */
    void SetGeometryShapeFunctionContainer(
        const GeometryShapeFunctionContainer<GeometryData::IntegrationMethod>& rGeometryShapeFunctionContainer) override
    {
        mGeometryData.SetGeometryShapeFunctionContainer(rGeometryShapeFunctionContainer);
    }

    ///@}
    ///@name Parent
    ///@{

    GeometryType& GetGeometryParent(IndexType Index) const override
    {
        return *mpGeometryParent;
    }

    void SetGeometryParent(GeometryType* pGeometryParent) override
    {
        mpGeometryParent = pGeometryParent;
    }

    ///@}
    ///@name Dynamic access to internals
    ///@{

    /// Calculate with Vector
    void Calculate(
        const Variable<Vector>& rVariable,
        Vector& rOutput) const override
    {
        if (rVariable == DETERMINANTS_OF_JACOBIAN_PARENT) {
            DeterminantOfJacobianParent(rOutput);
        }
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

    /// Returns the polynomial degree of the parent geometry
    SizeType PolynomialDegree(IndexType LocalDirectionIndex) const override
    {
        KRATOS_DEBUG_ERROR_IF_NOT(mpGeometryParent)
            << "Trying to call PolynomialDegree(LocalDirectionIndex) from quadrature point. "
            << "Pointer to parent is not assigned." << std::endl;

        return mpGeometryParent->PolynomialDegree(LocalDirectionIndex);
    }

    ///@}
    ///@name Coordinates
    ///@{

    /**
    * This method provides the global coordinates corresponding to the
    * local coordinates. The mapping is done on the parent.
    * Error if Parent is not assigned is only thrown in debug mode.
    */
    CoordinatesArrayType& GlobalCoordinates(
        CoordinatesArrayType& rResult,
        CoordinatesArrayType const& LocalCoordinates
    ) const override
    {
        KRATOS_DEBUG_ERROR_IF(mpGeometryParent == nullptr)
            << "Trying to call GlobalCoordinates(LocalCoordinates) from quadrature point. "
            << "Pointer to parent is not assigned." << std::endl;

        return mpGeometryParent->GlobalCoordinates(rResult, LocalCoordinates);
    }

    /* @brief returns the respective segment domain size of this
     *        quadrature point, computed on the parent of this geometry.
     *        Required for reduced quadrature point geometries (Not all
     *        nodes are part of this geometry - used for mapping).
     * @param rResult vector of results of this quadrature point.
     */
    CoordinatesArrayType& GlobalCoordinates(
        CoordinatesArrayType& rResult,
        CoordinatesArrayType const& LocalCoordinates,
        Matrix& DeltaPosition
    ) const override
    {
        KRATOS_DEBUG_ERROR_IF(mpGeometryParent == nullptr)
            << "Trying to call GlobalCoordinates(LocalCoordinates, DeltaPosition) from quadrature point. "
            << "Pointer to parent is not assigned." << std::endl;

        return mpGeometryParent->GlobalCoordinates(rResult, LocalCoordinates, DeltaPosition);
    }

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
    Matrix& Jacobian(
        Matrix& rResult,
        const CoordinatesArrayType& rCoordinates
    ) const override
    {
        KRATOS_DEBUG_ERROR_IF(mpGeometryParent == nullptr)
            << "Trying to call Jacobian(LocalCoordinates) from quadrature point. "
            << "Pointer to parent is not assigned." << std::endl;

        return mpGeometryParent->Jacobian(rResult, rCoordinates);
    }

    /** Determinant of jacobian in given point. Computed on parent
    *        geometry.
    * Error if Parent is not assigned is only thrown in debug mode.
    * @param rPoint point which determinant of jacobians has to
    *        be calculated in it.
    * @return Determinamt of jacobian matrix \f$ |J| \f$ in given
    *         point.
    */
    double DeterminantOfJacobian(
        const CoordinatesArrayType& rPoint
    ) const override
    {
        KRATOS_DEBUG_ERROR_IF(mpGeometryParent == nullptr)
            << "Trying to call DeterminantOfJacobian(rPoint) from quadrature point. "
            << "Pointer to parent is not assigned." << std::endl;

        return mpGeometryParent->DeterminantOfJacobian(rPoint);
    }


    /* @brief returns the respective segment length of this
     *        quadrature point. Length of vector always 1.
     * @param rResult vector of results of this quadrature point.
     */
    Vector& DeterminantOfJacobianParent(
        Vector& rResult) const
    {
        if (rResult.size() != 1)
            rResult.resize(1, false);

        rResult[0] = this->GetGeometryParent(0).DeterminantOfJacobian(this->IntegrationPoints()[0]);

        return rResult;
    }

    Matrix& InverseOfJacobian(
        Matrix& rResult,
        const CoordinatesArrayType& rCoordinates
    ) const override
    {
        KRATOS_DEBUG_ERROR_IF(mpGeometryParent == nullptr)
            << "Trying to call InverseOfJacobian(rPoint) from quadrature point. "
            << "Pointer to parent is not assigned." << std::endl;

        return mpGeometryParent->InverseOfJacobian(rResult, rCoordinates);
    }

    ///@}
    ///@name Shape Function
    ///@{

    Vector& ShapeFunctionsValues(
        Vector &rResult,
        const CoordinatesArrayType& rCoordinates
    ) const override
    {
        KRATOS_DEBUG_ERROR_IF(mpGeometryParent == nullptr)
            << "Trying to call ShapeFunctionsValues(rCoordinates) from quadrature point. "
            << "Pointer to parent is not assigned." << std::endl;

        return mpGeometryParent->ShapeFunctionsValues(rResult, rCoordinates);
    }

    virtual Matrix& ShapeFunctionsLocalGradients(
        Matrix& rResult,
        const CoordinatesArrayType& rPoint
    ) const override
    {
        KRATOS_DEBUG_ERROR_IF(mpGeometryParent == nullptr)
            << "Trying to call ShapeFunctionsLocalGradients(rPoint) from quadrature point. "
            << "Pointer to parent is not assigned." << std::endl;

        return mpGeometryParent->ShapeFunctionsLocalGradients(rResult, rPoint);
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
        return GeometryData::KratosGeometryType::Kratos_Quadrature_Point_Geometry;
    }

    ///@}
    ///@name Input and output
    ///@{

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
    ///@}

protected:

    ///@name Constructor
    ///@{

    /// Standard Constructor
    QuadraturePointGeometry()
        : BaseType(
            PointsArrayType(),
            &mGeometryData)
        , mGeometryData(
            &msGeometryDimension,
            GeometryData::IntegrationMethod::GI_GAUSS_1,
            {}, {}, {})
    {
    }

    ///@}

private:
    ///@name Static Member Variables
    ///@{

    static const GeometryDimension msGeometryDimension;

    ///@}
    ///@name Member Variables
    ///@{

    GeometryData mGeometryData;

    // quatrature point can be related to a parent geometry. To keep the connection,
    // this geometry is related to the integration point.
    GeometryType* mpGeometryParent = nullptr;

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save( Serializer& rSerializer ) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, BaseType );

        rSerializer.save("IntegrationPoints", mGeometryData.IntegrationPoints());
        rSerializer.save("ShapeFunctionsValues", mGeometryData.ShapeFunctionsValues());
        rSerializer.save("ShapeFunctionsLocalGradients", mGeometryData.ShapeFunctionsLocalGradients());
    }

    void load( Serializer& rSerializer ) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BaseType );

        IntegrationPointsContainerType integration_points;
        ShapeFunctionsValuesContainerType shape_functions_values;
        ShapeFunctionsLocalGradientsContainerType shape_functions_local_gradients;

        rSerializer.load("IntegrationPoints", integration_points[static_cast<int>(GeometryData::IntegrationMethod::GI_GAUSS_1)]);
        rSerializer.load("ShapeFunctionsValues", shape_functions_values[static_cast<int>(GeometryData::IntegrationMethod::GI_GAUSS_1)]);
        rSerializer.load("ShapeFunctionsLocalGradients", shape_functions_local_gradients[static_cast<int>(GeometryData::IntegrationMethod::GI_GAUSS_1)]);

        mGeometryData.SetGeometryShapeFunctionContainer(GeometryShapeFunctionContainer<GeometryData::IntegrationMethod>(
            GeometryData::IntegrationMethod::GI_GAUSS_1,
            integration_points,
            shape_functions_values,
            shape_functions_local_gradients));


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
    QuadraturePointGeometry<TPointType, TWorkingSpaceDimension, TLocalSpaceDimension, TDimension>& rThis );

/// output stream function
template<class TPointType,
    int TWorkingSpaceDimension,
    int TLocalSpaceDimension,
    int TDimension>
inline std::ostream& operator << (
    std::ostream& rOStream,
    const QuadraturePointGeometry<TPointType, TWorkingSpaceDimension, TLocalSpaceDimension, TDimension>& rThis )
{
    rThis.PrintInfo( rOStream );
    rOStream << std::endl;
    rThis.PrintData( rOStream );

    return rOStream;
}

///@}
///@name Type Dimension Definition
///@{

template<class TPointType,
    int TWorkingSpaceDimension,
    int TLocalSpaceDimension,
    int TDimension>
const GeometryDimension QuadraturePointGeometry<
    TPointType,
    TWorkingSpaceDimension,
    TLocalSpaceDimension,
    TDimension>::msGeometryDimension(
        TWorkingSpaceDimension,
        TLocalSpaceDimension);

///@}

}  // namespace Kratos.

//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//

#pragma once

// System includes

// External includes

// Project includes
#include "geometries/geometry.h"
#include "geometries/nurbs_curve_geometry.h"
#include "geometries/nurbs_shape_function_utilities/nurbs_interval.h"


namespace Kratos
{
///@name Kratos Classes
///@{

/**
 * @class BrepCurve
 * @ingroup KratosCore
 * @brief The BrepCurve acts as topology for curves. Those
 *        can be enclosed by a certain set of points.
 */
template<class TContainerPointType, class TContainerPointEmbeddedType = TContainerPointType>
class BrepCurve
    : public Geometry<typename TContainerPointType::value_type>
{
public:
    ///@}
    ///@name Type Definitions
    ///@{

    /// Pointer definition of BrepCurve
    KRATOS_CLASS_POINTER_DEFINITION( BrepCurve );

    typedef typename TContainerPointType::value_type PointType;

    typedef Geometry<typename TContainerPointType::value_type> BaseType;
    typedef Geometry<typename TContainerPointType::value_type> GeometryType;
    typedef typename GeometryType::Pointer GeometryPointer;

    typedef GeometryData::IntegrationMethod IntegrationMethod;

    typedef NurbsCurveGeometry<3, TContainerPointType> NurbsCurveType;

    typedef typename BaseType::GeometriesArrayType GeometriesArrayType;

    typedef typename BaseType::IndexType IndexType;
    typedef typename BaseType::SizeType SizeType;

    typedef typename BaseType::PointsArrayType PointsArrayType;
    typedef typename BaseType::CoordinatesArrayType CoordinatesArrayType;
    typedef typename BaseType::IntegrationPointsArrayType IntegrationPointsArrayType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor for untrimmed curve
    BrepCurve(
        typename NurbsCurveType::Pointer pCurve)
        : BaseType(PointsArrayType(), &msGeometryData)
        , mpNurbsCurve(pCurve)
    {
        mIsTrimmed = false;
    }

    /// Constructor default
    explicit BrepCurve(const PointsArrayType& ThisPoints)
        : BaseType(ThisPoints, &msGeometryData)
    {
    }

    /// Copy constructor.
    BrepCurve(BrepCurve const& rOther)
        : BaseType(rOther)
        , mpNurbsCurve(rOther.mpNurbsCurve)
        , mIsTrimmed(rOther.mIsTrimmed)
    {
    }

    /// Copy constructor from a geometry with different point type.
    template<class TOtherContainerPointType, class TOtherContainerPointEmbeddedType>
    explicit BrepCurve(
        BrepCurve<TOtherContainerPointType, TOtherContainerPointEmbeddedType> const& rOther )
        : BaseType(rOther)
        , mpNurbsCurve(rOther.mpNurbsCurve)
        , mIsTrimmed(rOther.mIsTrimmed)
    {
    }

    /// Destructor
    ~BrepCurve() override = default;

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    BrepCurve& operator=( const BrepCurve& rOther )
    {
        BaseType::operator=( rOther );
        mpNurbsCurve = rOther.mpNurbsCurve;
        mIsTrimmed = rOther.mIsTrimmed;
        return *this;
    }

    /// Assignment operator for geometries with different point type.
    template<class TOtherContainerPointType, class TOtherContainerPointEmbeddedType>
    BrepCurve& operator=( BrepCurve<TOtherContainerPointType, TOtherContainerPointEmbeddedType> const & rOther )
    {
        BaseType::operator=( rOther );
        mpNurbsCurve = rOther.mpNurbsCurve;
        mIsTrimmed = rOther.mIsTrimmed;
        return *this;
    }

    ///@}
    ///@name Operations
    ///@{

    typename BaseType::Pointer Create(PointsArrayType const& ThisPoints) const override
    {
        return typename BaseType::Pointer(new BrepCurve(ThisPoints));
    }

    ///@}
    ///@name Access to Geometry Parts
    ///@{

    /**
    * @brief This function returns the pointer of the geometry
    *        which is corresponding to the index.
    *        Curve of the geometry is accessible with GeometryType::BACKGROUND_GEOMETRY_INDEX.
    * @param Index: trim_index or GeometryType::BACKGROUND_GEOMETRY_INDEX.
    * @return pointer of geometry, corresponding to the index.
    */
    GeometryPointer pGetGeometryPart(const IndexType Index) override
    {
        const auto& const_this = *this;
        return std::const_pointer_cast<GeometryType>(
            const_this.pGetGeometryPart(Index));
    }

    /**
    * @brief This function returns the pointer of the geometry
    *        which is corresponding to the trim index.
    *        Curve of the geometry is accessible with GeometryType::BACKGROUND_GEOMETRY_INDEX.
    * @param Index: trim_index or GeometryType::BACKGROUND_GEOMETRY_INDEX.
    * @return pointer of geometry, corresponding to the index.
    */
    const GeometryPointer pGetGeometryPart(const IndexType Index) const override
    {
        if (Index == GeometryType::BACKGROUND_GEOMETRY_INDEX)
            return mpNurbsCurve;

        KRATOS_ERROR << "Index " << Index << " not existing in BrepCurve: "
            << this->Id() << std::endl;
    }

    /**
    * @brief This function is used to check if this BrepCurve
    *        has certain trim or curve object.
    * @param Index of the geometry part.
    * @return true if has trim or curve
    */
    bool HasGeometryPart(const IndexType Index) const override
    {
        if (Index == GeometryType::BACKGROUND_GEOMETRY_INDEX)
            return true;

        return false;
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
            mpNurbsCurve->Calculate(rVariable, rOutput);
        }
    }

    ///@}
    ///@name Mathematical Informations
    ///@{

    /// Return polynomial degree of the nurbs surface
    SizeType PolynomialDegree(IndexType LocalDirectionIndex) const override
    {
        return mpNurbsCurve->PolynomialDegree(LocalDirectionIndex);
    }

    ///@}
    ///@name Information
    ///@{

    /*
    * @brief checks if the BrepCurve has any boundary trim information.
    * @return true if has no trimming.
    */
    bool IsTrimmed() const
    {
        return mIsTrimmed;
    }

    /// Returns number of points of NurbsSurface.
    SizeType PointsNumberInDirection(IndexType DirectionIndex) const override
    {
        return mpNurbsCurve->PointsNumberInDirection(DirectionIndex);
    }

    ///@}
    ///@name Geometrical Operations
    ///@{

    /// Provides the center of the underlying surface
    Point Center() const override
    {
        return mpNurbsCurve->Center();
    }

    /**
    * @brief Calls projection of its nurbs curve.
    *        Projects a certain point on the geometry, or finds
    *        the closest point, depending on the provided
    *        initial guess. The external point does not necessary
    *        lay on the geometry.
    *        It shall deal as the interface to the mathematical
    *        projection function e.g. the Newton-Raphson.
    *        Thus, the breaking criteria does not necessarily mean
    *        that it found a point on the surface, if it is really
    *        the closest or not. It shows only if the breaking
    *        criteria, defined by the tolerance is reached.
    *
    *        This function requires an initial guess, provided by
    *        rProjectedPointLocalCoordinates.
    *
    *        This function can be a very costly operation.
    *
    * @param rPointGlobalCoordinates the point to which the
    *        projection has to be found.
    * @param rProjectedPointLocalCoordinates the location of the
    *        projection in local coordinates.
    *        The variable is as initial guess!
    * @param Tolerance accepted of orthogonal error to projection.
    * @return It is chosen to take an int as output parameter to
    *         keep more possibilities within the interface.
    *         0 -> failed
    *         1 -> converged
    */
    int ProjectionPointGlobalToLocalSpace(
        const CoordinatesArrayType& rPointGlobalCoordinates,
        CoordinatesArrayType& rProjectedPointLocalCoordinates,
        const double Tolerance = std::numeric_limits<double>::epsilon()
        ) const override
    {
        return mpNurbsCurve->ProjectionPointGlobalToLocalSpace(
            rPointGlobalCoordinates, rProjectedPointLocalCoordinates, Tolerance);
    }

    /*
    * @brief This method maps from dimension space to working space.
    * @param rResult array_1d<double, 3> with the coordinates in working space
    * @param LocalCoordinates The local coordinates in dimension space
    * @return array_1d<double, 3> with the coordinates in working space
    * @see PointLocalCoordinates
    */
    CoordinatesArrayType& GlobalCoordinates(
        CoordinatesArrayType& rResult,
        const CoordinatesArrayType& rLocalCoordinates
    ) const override
     {
        mpNurbsCurve->GlobalCoordinates(rResult, rLocalCoordinates);

        return rResult;
    }

    ///@}
    ///@name Integration Info
    ///@{

    /// Provides the default integration dependent on the polynomial degree.
    IntegrationInfo GetDefaultIntegrationInfo() const override
    {
        return mpNurbsCurve->GetDefaultIntegrationInfo();
    }

    ///@}
    ///@name Integration Points
    ///@{

    /* Creates integration points on the nurbs surface of this geometry.
     * @param return integration points.
     */
    void CreateIntegrationPoints(
        IntegrationPointsArrayType& rIntegrationPoints,
        IntegrationInfo& rIntegrationInfo) const override
    {
        std::vector<double> spans;
        mpNurbsCurve->SpansLocalSpace(spans);

        IntegrationPointUtilities::CreateIntegrationPoints1D(
            rIntegrationPoints, spans, rIntegrationInfo);
    }

    ///@}
    ///@name Quadrature Point Geometries
    ///@{

    /* @brief calls function of underlying nurbs curve.
     *
     * @param rResultGeometries list of quadrature point geometries.
     * @param NumberOfShapeFunctionDerivatives the number of evaluated
     *        derivatives of shape functions at the quadrature point geometries.
     * @param rIntegrationPoints the corresponding integration points.
     *
     * @see quadrature_point_geometry.h
     */
    void CreateQuadraturePointGeometries(
        GeometriesArrayType& rResultGeometries,
        IndexType NumberOfShapeFunctionDerivatives,
        const IntegrationPointsArrayType& rIntegrationPoints,
        IntegrationInfo& rIntegrationInfo) override
    {
        mpNurbsCurve->CreateQuadraturePointGeometries(
            rResultGeometries, NumberOfShapeFunctionDerivatives, rIntegrationPoints, rIntegrationInfo);

        for (IndexType i = 0; i < rResultGeometries.size(); ++i) {
            rResultGeometries(i)->SetGeometryParent(this);
        }
    }

    ///@}
    ///@name Shape Function
    ///@{

    Vector& ShapeFunctionsValues(
        Vector &rResult,
        const CoordinatesArrayType& rCoordinates) const override
    {
        mpNurbsCurve->ShapeFunctionsValues(rResult, rCoordinates);

        return rResult;
    }

    Matrix& ShapeFunctionsLocalGradients(
        Matrix& rResult,
        const CoordinatesArrayType& rCoordinates) const override
    {
        mpNurbsCurve->ShapeFunctionsLocalGradients(rResult, rCoordinates);

        return rResult;
    }

    ///@}
    ///@name Geometry Classification
    ///@{

    /**
     * @brief Gets the geometry family.
     * @details This function returns the family type of the geometry. The geometry family categorizes the geometry into a broader classification, aiding in its identification and processing.
     * @return GeometryData::KratosGeometryFamily The geometry family.
     */
    GeometryData::KratosGeometryFamily GetGeometryFamily() const override
    {
        return GeometryData::KratosGeometryFamily::Kratos_Brep;
    }

    /**
     * @brief Gets the geometry type.
     * @details This function returns the specific type of the geometry. The geometry type provides a more detailed classification of the geometry.
     * @return GeometryData::KratosGeometryType The specific geometry type.
     */
    GeometryData::KratosGeometryType GetGeometryType() const override
    {
        return GeometryData::KratosGeometryType::Kratos_Brep_Curve;
    }

    ///@}
    ///@name Information
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "BrepCurve";
    }

    /// Print information about this object.
    void PrintInfo( std::ostream& rOStream ) const override
    {
        rOStream << "BrepCurve";
    }

    /// Print object's data.
    void PrintData( std::ostream& rOStream ) const override
    {
        BaseType::PrintData( rOStream );
        std::cout << std::endl;
        rOStream << "    BrepCurve " << std::endl;
    }

    ///@}

private:
    ///@name Static Member Variables
    ///@{

    static const GeometryData msGeometryData;

    static const GeometryDimension msGeometryDimension;

    ///@}
    ///@name Member Variables
    ///@{

    typename NurbsCurveType::Pointer mpNurbsCurve;

    /** IsTrimmed is used to optimize processes as
    *   e.g. creation of integration domain.
    */
    bool mIsTrimmed;

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save( Serializer& rSerializer ) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, BaseType );
        rSerializer.save("NurbsCurve", mpNurbsCurve);
        rSerializer.save("IsTrimmed", mIsTrimmed);
    }

    void load( Serializer& rSerializer ) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BaseType );
        rSerializer.load("NurbsCurve", mpNurbsCurve);
        rSerializer.load("IsTrimmed", mIsTrimmed);
    }

    BrepCurve()
        : BaseType( PointsArrayType(), &msGeometryData )
    {}

    ///@}

}; // Class BrepCurve

///@name Input and output
///@{

/// input stream functions
template<class TContainerPointType, class TContainerPointEmbeddedType = TContainerPointType> inline std::istream& operator >> (
    std::istream& rIStream,
    BrepCurve<TContainerPointType, TContainerPointEmbeddedType>& rThis );

/// output stream functions
template<class TContainerPointType, class TContainerPointEmbeddedType = TContainerPointType> inline std::ostream& operator << (
    std::ostream& rOStream,
    const BrepCurve<TContainerPointType, TContainerPointEmbeddedType>& rThis )
{
    rThis.PrintInfo( rOStream );
    rOStream << std::endl;
    rThis.PrintData( rOStream );
    return rOStream;
}

///@}
///@name Static Type Declarations
///@{

template<class TContainerPointType, class TContainerPointEmbeddedType> const
GeometryData BrepCurve<TContainerPointType, TContainerPointEmbeddedType>::msGeometryData(
    &msGeometryDimension,
    GeometryData::IntegrationMethod::GI_GAUSS_1,
    {}, {}, {});

template<class TContainerPointType, class TContainerPointEmbeddedType>
const GeometryDimension BrepCurve<TContainerPointType, TContainerPointEmbeddedType>::msGeometryDimension(3, 1);

///@}
}// namespace Kratos.

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
//                   Andreas Apostolatos
//                   Pooyan Dadvand
//                   Philipp Bucher
//

#if !defined(KRATOS_BREP_CURVE_ON_SURFACE_3D_H_INCLUDED )
#define  KRATOS_BREP_CURVE_ON_SURFACE_3D_H_INCLUDED

// System includes

// External includes

// Project includes
#include "geometries/geometry.h"
#include "geometries/nurbs_shape_function_utilities/nurbs_interval.h"

#include "geometries/nurbs_curve_on_surface_geometry.h"

#include "utilities/nurbs_utilities/projection_nurbs_geometry_utilities.h"

namespace Kratos
{
///@name Kratos Classes
///@{

/**
 * @class BrepCurveOnSurface
 * @ingroup KratosCore
 * @brief The BrepCurveOnSurface acts as topology for curves on surfaces.
 */
template<class TContainerPointType, class TContainerPointEmbeddedType = TContainerPointType>
class BrepCurveOnSurface
    : public Geometry<typename TContainerPointType::value_type>
{
public:
    ///@}
    ///@name Type Definitions
    ///@{

    /** Pointer definition of BrepCurveOnSurface */
    KRATOS_CLASS_POINTER_DEFINITION( BrepCurveOnSurface );

    typedef typename TContainerPointType::value_type PointType;

    typedef Geometry<typename TContainerPointType::value_type> BaseType;
    typedef Geometry<typename TContainerPointType::value_type> GeometryType;

    typedef GeometryData::IntegrationMethod IntegrationMethod;

    typedef NurbsSurfaceGeometry<3, TContainerPointType> NurbsSurfaceType;
    typedef NurbsCurveGeometry<2, TContainerPointEmbeddedType> NurbsCurveType;

    typedef NurbsCurveOnSurfaceGeometry<3, TContainerPointEmbeddedType, TContainerPointType> NurbsCurveOnSurfaceType;

    typedef typename NurbsCurveOnSurfaceType::Pointer NurbsCurveOnSurfacePointerType;

    typedef typename BaseType::GeometriesArrayType GeometriesArrayType;

    typedef typename BaseType::IndexType IndexType;
    typedef typename BaseType::SizeType SizeType;

    typedef typename BaseType::PointsArrayType PointsArrayType;
    typedef typename BaseType::CoordinatesArrayType CoordinatesArrayType;
    typedef typename BaseType::IntegrationPointsArrayType IntegrationPointsArrayType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// constructor for untrimmed surface
    BrepCurveOnSurface( 
        typename NurbsSurfaceType::Pointer pSurface,
        typename NurbsCurveType::Pointer pCurve,
        bool SameCurveDirection = true)
        : BaseType(PointsArrayType(), &msGeometryData)
        , mpCurveOnSurface(
            Kratos::make_shared<NurbsCurveOnSurfaceType>(
                pSurface, pCurve))
        , mCurveNurbsInterval(pCurve->DomainInterval())
        , mSameCurveDirection(SameCurveDirection)
    {
    }

    /// constructor for trimmed surface
    BrepCurveOnSurface(
        typename NurbsSurfaceType::Pointer pSurface,
        typename NurbsCurveType::Pointer pCurve,
        NurbsInterval CurveNurbsInterval,
        bool SameCurveDirection = true)
        : BaseType(PointsArrayType(), &msGeometryData)
        , mpCurveOnSurface(
            Kratos::make_shared<NurbsCurveOnSurfaceType>(
                pSurface, pCurve))
        , mCurveNurbsInterval(CurveNurbsInterval)
        , mSameCurveDirection(SameCurveDirection)
    {
    }

    /// constructor for untrimmed surface with curve on surface
    BrepCurveOnSurface(
        NurbsCurveOnSurfacePointerType pNurbsCurveOnSurface,
        bool SameCurveDirection = true)
        : BaseType(PointsArrayType(), &msGeometryData)
        , mpCurveOnSurface(pNurbsCurveOnSurface)
        , mCurveNurbsInterval(pNurbsCurveOnSurface->DomainInterval())
        , mSameCurveDirection(SameCurveDirection)
    {
    }

    /// constructor for trimmed surface with curve on surface
    BrepCurveOnSurface(
        NurbsCurveOnSurfacePointerType pNurbsCurveOnSurface,
        NurbsInterval CurveNurbsInterval,
        bool SameCurveDirection = true)
        : BaseType(PointsArrayType(), &msGeometryData)
        , mpCurveOnSurface(pNurbsCurveOnSurface)
        , mCurveNurbsInterval(CurveNurbsInterval)
        , mSameCurveDirection(SameCurveDirection)
    {
    }

    BrepCurveOnSurface()
        : BaseType(PointsArrayType(), &msGeometryData)
    {}

    explicit BrepCurveOnSurface(const PointsArrayType& ThisPoints)
        : BaseType(ThisPoints, &msGeometryData)
    {
    }

    /// Copy constructor.
    BrepCurveOnSurface( BrepCurveOnSurface const& rOther )
        : BaseType( rOther )
        , mpCurveOnSurface(rOther.mpCurveOnSurface)
        , mCurveNurbsInterval(rOther.mCurveNurbsInterval)
        , mSameCurveDirection(rOther.mSameCurveDirection)
    {
    }

    /// Copy constructor from a geometry with different point type.
    template<class TOtherContainerPointType, class TOtherContainerPointEmbeddedType>
    explicit BrepCurveOnSurface(
        BrepCurveOnSurface<TOtherContainerPointType, TOtherContainerPointEmbeddedType> const& rOther )
        : BaseType( rOther )
        , mpCurveOnSurface(rOther.mpCurveOnSurface)
        , mCurveNurbsInterval(rOther.mCurveNurbsInterval)
        , mSameCurveDirection(rOther.mSameCurveDirection)
    {
    }

    /// Destructor
    ~BrepCurveOnSurface() override = default;

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator
    BrepCurveOnSurface& operator=( const BrepCurveOnSurface& rOther )
    {
        BaseType::operator=( rOther );
        mpCurveOnSurface = rOther.mpCurveOnSurface;
        mCurveNurbsInterval = rOther.mCurveNurbsInterval;
        mSameCurveDirection = rOther.mSameCurveDirection;
        return *this;
    }

    /// Assignment operator with different point type
    template<class TOtherContainerPointType, class TOtherContainerPointEmbeddedType>
    BrepCurveOnSurface& operator=( BrepCurveOnSurface<TOtherContainerPointType, TOtherContainerPointEmbeddedType> const & rOther )
    {
        BaseType::operator=( rOther );
        mpCurveOnSurface = rOther.mpCurveOnSurface;
        mCurveNurbsInterval = rOther.mCurveNurbsInterval;
        mSameCurveDirection = rOther.mSameCurveDirection;
        return *this;
    }

    ///@}
    ///@name Operations
    ///@{

    typename BaseType::Pointer Create( PointsArrayType const& ThisPoints ) const override
    {
        return typename BaseType::Pointer( new BrepCurveOnSurface( ThisPoints ) );
    }

    ///@}
    ///@name Mathematical Informations
    ///@{

    /// Return polynomial degree of the nurbs curve on surface
    SizeType PolynomialDegree(IndexType LocalDirectionIndex) const override
    {
        return mpCurveOnSurface->PolynomialDegree(LocalDirectionIndex);
    }

    ///@}
    ///@name Set/ Get functions
    ///@{

    /*
    * @brief Indicates if the NURBS-curve is pointing in the same direction
    *        as the B-Rep curve.
    * @return true -> brep curve and nurbs curve point in same direction.
    *        false -> brep curve and nurbs curve point in opposite directions.
    */
    bool HasSameCurveDirection()
    {
        return mSameCurveDirection;
    }

    /// Returns the const NurbsCurveOnSurface::Pointer of this brep.
    const NurbsInterval DomainInterval() const
    {
        return mCurveNurbsInterval;
    }

    /// Returns the NurbsCurveOnSurface::Pointer of this brep.
    NurbsCurveOnSurfacePointerType pGetCurveOnSurface()
    {
        return mpCurveOnSurface;
    }

    /// Returns the const NurbsCurveOnSurface::Pointer of this brep.
    const NurbsCurveOnSurfacePointerType pGetCurveOnSurface() const
    {
        return mpCurveOnSurface;
    }

    /// Returns number of points of NurbsCurveOnSurface.
    SizeType PointsNumberInDirection(IndexType DirectionIndex) const override
    {
        return mpCurveOnSurface->PointsNumberInDirection(DirectionIndex);
    }

    ///@}
    ///@name Curve Properties
    ///@{

    /* @brief Provides intersections of the nurbs curve with the knots of the surface,
     *         using the interval of this curve.
     * @param vector of span intervals.
     * @param index of chosen direction, for curves always 0.
     */
    void Spans(std::vector<double>& rSpans, IndexType DirectionIndex = 0) const override
    {
        if (rSpans.size() > 0) {
            double external_start = rSpans[0];
            double external_end = rSpans[rSpans.size() - 1];
            mCurveNurbsInterval.IsInside(external_start);
            mCurveNurbsInterval.IsInside(external_end);

            NurbsInterval trimmed_external_interval(external_start, external_end);

            std::vector<double> new_spans;
            for (IndexType i = 0; i < rSpans.size(); i++) {
                double temp = rSpans[i];
                if (trimmed_external_interval.IsInside(temp)) {
                    new_spans.push_back(temp);
                }
            }
            rSpans = new_spans;
        }
        else {
            rSpans.resize(2);
            rSpans[0] = mCurveNurbsInterval.GetT0();
            rSpans[1] = mCurveNurbsInterval.GetT1();
        }

        mpCurveOnSurface->Spans(rSpans);
    }

    ///@}
    ///@name Geometrical Operations
    ///@{

    /// Provides the center of the underlying curve on surface
    Point Center() const override
    {
        return mpCurveOnSurface->Center();
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
        return mpCurveOnSurface->GlobalCoordinates(rResult, rLocalCoordinates);
    }

    /* @brief This method maps from dimension space to working space and computes the
     *        number of derivatives at the dimension parameter.
     * From Piegl and Tiller, The NURBS Book, Algorithm A3.2/ A4.2
     * @param LocalCoordinates The local coordinates in dimension space
     * @param Derivative Number of computed derivatives
     * @return std::vector<array_1d<double, 3>> with the coordinates in working space
     * @see PointLocalCoordinates
     */
    void GlobalSpaceDerivatives(
        std::vector<CoordinatesArrayType>& rGlobalSpaceDerivatives,
        const CoordinatesArrayType& rLocalCoordinates,
        const SizeType DerivativeOrder) const override
    {
        return mpCurveOnSurface->GlobalSpaceDerivatives(
            rGlobalSpaceDerivatives, rLocalCoordinates, DerivativeOrder);
    }

    ///@}
    ///@name Projection
    ///@{

    /* Makes projection of rPointGlobalCoordinates to
     * the closest point rProjectedPointGlobalCoordinates on the curve,
     * with local coordinates rProjectedPointLocalCoordinates.
     *
     * Condiders limits of this BrepCurveOnSurface as borders.
     *
     * @param Tolerance is the breaking criteria.
     * @return 1 -> projection succeeded
     *         0 -> projection failed
     */
    int ProjectionPoint(
        const CoordinatesArrayType& rPointGlobalCoordinates,
        CoordinatesArrayType& rProjectedPointGlobalCoordinates,
        CoordinatesArrayType& rProjectedPointLocalCoordinates,
        const double Tolerance = std::numeric_limits<double>::epsilon()
        ) const override
    {
        const bool success = ProjectionNurbsGeometryUtilities::NewtonRaphsonCurve(
            rProjectedPointLocalCoordinates,
            rPointGlobalCoordinates,
            rProjectedPointGlobalCoordinates,
            *this,
            20, Tolerance);

        return (success)
            ? 1
            : 0;
    }

    ///@}
    ///@name IsInside
    ///@{

    /// returns if rPointLocalCoordinates[0] is inside -> 1 or ouside -> 0
    int IsInsideLocalSpace(
        const CoordinatesArrayType& rPointLocalCoordinates,
        const double Tolerance = std::numeric_limits<double>::epsilon()
        ) const override
    {
        const double min_parameter = mCurveNurbsInterval.MinParameter();
        if (rPointLocalCoordinates[0] < min_parameter) {
            return 0;
        }

        const double max_parameter = mCurveNurbsInterval.MaxParameter();
        if (rPointLocalCoordinates[0] > max_parameter) {
            return 0;
        }

        return 1;
    }

    /* Returns if rPointLocalCoordinates[0] is inside -> 1 or ouside -> 0
     * and sets it to the closest border */
    int SetInsideLocalSpace(
        CoordinatesArrayType& rPointLocalCoordinates,
        const double Tolerance = std::numeric_limits<double>::epsilon()
        ) const override
    {
        return mCurveNurbsInterval.IsInside(rPointLocalCoordinates[0]);
    }


    ///@}
    ///@name Integration Points
    ///@{

    /* Creates integration points on the nurbs surface of this geometry
     * with the domain limits of this brep curve on surface.
     * @param return integration points.
     */
    void CreateIntegrationPoints(
        IntegrationPointsArrayType& rIntegrationPoints,
        IntegrationInfo& rIntegrationInfo) const override
    {
        std::vector<double> spans;
        if (!rIntegrationInfo.HasSpansInDirection(0)) {
            std::vector<double> spans;
            Spans(spans);
            rIntegrationInfo.SetSpans(spans, 0);
        }

        mpCurveOnSurface->CreateIntegrationPoints(
            rIntegrationPoints, rIntegrationInfo);
    }

    ///@}
    ///@name Quadrature Point Geometries
    ///@{

    /* @brief creates a list of quadrature point geometries
     *        from a list of integration points on the
     *        curve on surface of this geometry.
     *
     * @param rResultGeometries list of quadrature point geometries.
     * @param rIntegrationPoints list of provided integration points.
     * @param NumberOfShapeFunctionDerivatives the number of evaluated
     *        derivatives of shape functions at the quadrature point geometries.
     *
     * @see quadrature_point_geometry.h
     */
    void CreateQuadraturePointGeometries(
        GeometriesArrayType& rResultGeometries,
        IndexType NumberOfShapeFunctionDerivatives,
        const IntegrationPointsArrayType& rIntegrationPoints) override
    {
        mpCurveOnSurface->CreateQuadraturePointGeometries(
            rResultGeometries, NumberOfShapeFunctionDerivatives, rIntegrationPoints);

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
        return mpCurveOnSurface->ShapeFunctionsValues(rResult, rCoordinates);
    }

    Matrix& ShapeFunctionsLocalGradients(
        Matrix& rResult,
        const CoordinatesArrayType& rCoordinates) const override
    {
        return mpCurveOnSurface->ShapeFunctionsLocalGradients(rResult, rCoordinates);
    }

    ///@}
    ///@name Input and output
    ///@{


    ///@}
    ///@name Information
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "Brep face curve";
    }

    /// Print information about this object.
    void PrintInfo( std::ostream& rOStream ) const override
    {
        rOStream << "Brep face curve";
    }

    /// Print object's data.
    void PrintData( std::ostream& rOStream ) const override
    {
        BaseType::PrintData( rOStream );
        std::cout << std::endl;
        rOStream << "    Brep face curve : " << std::endl;
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

    NurbsCurveOnSurfacePointerType mpCurveOnSurface;

    NurbsInterval mCurveNurbsInterval;

    /** true-> brep curve and nurbs curve point in same direction.
    *  false-> brep curve and nurbs curve point in opposite directions. */
    bool mSameCurveDirection;

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save( Serializer& rSerializer ) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, BaseType );
        rSerializer.save("CurveOnSurface", mpCurveOnSurface);
        rSerializer.save("NurbsInterval", mCurveNurbsInterval);
        rSerializer.save("SameCurveDirection", mSameCurveDirection);
    }

    void load( Serializer& rSerializer ) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BaseType );
        rSerializer.load("CurveOnSurface", mpCurveOnSurface);
        rSerializer.load("NurbsInterval", mCurveNurbsInterval);
        rSerializer.load("SameCurveDirection", mSameCurveDirection);
    }

    ///@}
}; // Class BrepCurveOnSurface

///@name Input and output
///@{

/// input stream functions
template<class TContainerPointType, class TContainerPointEmbeddedType = TContainerPointType> inline std::istream& operator >> (
    std::istream& rIStream,
    BrepCurveOnSurface<TContainerPointType, TContainerPointEmbeddedType>& rThis );

/// output stream functions
template<class TContainerPointType, class TContainerPointEmbeddedType = TContainerPointType> inline std::ostream& operator << (
    std::ostream& rOStream,
    const BrepCurveOnSurface<TContainerPointType, TContainerPointEmbeddedType>& rThis )
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
GeometryData BrepCurveOnSurface<TContainerPointType, TContainerPointEmbeddedType>::msGeometryData(
    &msGeometryDimension,
    GeometryData::GI_GAUSS_1,
    {}, {}, {});

template<class TContainerPointType, class TContainerPointEmbeddedType>
const GeometryDimension BrepCurveOnSurface<TContainerPointType, TContainerPointEmbeddedType>::msGeometryDimension(
    1, 3, 2);

///@}
}// namespace Kratos.

#endif // KRATOS_BREP_CURVE_ON_SURFACE_3D_H_INCLUDED  defined

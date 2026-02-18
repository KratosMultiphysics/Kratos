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
//                   Nicolò Antonelli
//                   Andrea Gorgi
//

#if !defined(KRATOS_BREP_FACE_3D_H_INCLUDED )
#define  KRATOS_BREP_FACE_3D_H_INCLUDED

// System includes

// External includes

// Project includes
#include "geometries/geometry.h"
#include "geometries/brep_curve_on_surface.h"
#include "geometries/nurbs_shape_function_utilities/nurbs_interval.h"

// trimming integration
#include "utilities/geometry_utilities/brep_trimming_utilities.h"
#include "utilities/geometry_utilities/brep_sbm_utilities.h"

namespace Kratos
{
///@name Kratos Classes
///@{

/**
 * @class BrepSurface
 * @ingroup KratosCore
 * @brief The BrepSurface acts as topology for faces. Those
 *        can be enclosed by a certain set of brep face curves.
 * @tparam TShiftedBoundary Boolean flag indicating whether is 
 *        defined with shifted boundary conditions.
 */
template<class TContainerPointType, bool TShiftedBoundary, class TContainerPointEmbeddedType = TContainerPointType>
class BrepSurface
    : public Geometry<typename TContainerPointType::value_type>
{
public:
    ///@}
    ///@name Type Definitions
    ///@{

    /** Pointer definition of BrepSurface */
    KRATOS_CLASS_POINTER_DEFINITION( BrepSurface );

    typedef typename TContainerPointType::value_type PointType;

    typedef Geometry<typename TContainerPointType::value_type> BaseType;
    typedef Geometry<typename TContainerPointType::value_type> GeometryType;
    typedef typename GeometryType::Pointer GeometryPointer;
    using GeometrySurrogateArrayType = DenseVector<GeometryPointer>;

    typedef GeometryData::IntegrationMethod IntegrationMethod;

    typedef NurbsSurfaceGeometry<3, TContainerPointType> NurbsSurfaceType;
    typedef BrepCurveOnSurface<TContainerPointType, TShiftedBoundary, TContainerPointEmbeddedType> BrepCurveOnSurfaceType;

    typedef BrepTrimmingUtilities<TShiftedBoundary> BrepTrimmingUtilitiesType;
    typedef BrepSbmUtilities<Node> BrepSbmUtilitiesType;

    typedef DenseVector<typename BrepCurveOnSurfaceType::Pointer> BrepCurveOnSurfaceArrayType;
    typedef DenseVector<typename BrepCurveOnSurfaceType::Pointer> BrepCurveOnSurfaceLoopType;
    typedef DenseVector<DenseVector<typename BrepCurveOnSurfaceType::Pointer>> BrepCurveOnSurfaceLoopArrayType;

    typedef typename BaseType::GeometriesArrayType GeometriesArrayType;

    typedef typename BaseType::IndexType IndexType;
    typedef typename BaseType::SizeType SizeType;

    typedef typename BaseType::PointsArrayType PointsArrayType;
    typedef typename BaseType::CoordinatesArrayType CoordinatesArrayType;
    typedef typename BaseType::IntegrationPointsArrayType IntegrationPointsArrayType;

    typedef GeometryShapeFunctionContainer<GeometryData::IntegrationMethod>::IntegrationPointsContainerType IntegrationPointsContainerType;
    typedef GeometryShapeFunctionContainer<GeometryData::IntegrationMethod>::ShapeFunctionsValuesContainerType ShapeFunctionsValuesContainerType;
    typedef GeometryShapeFunctionContainer<GeometryData::IntegrationMethod>::ShapeFunctionsLocalGradientsContainerType ShapeFunctionsLocalGradientsContainerType;
    typedef GeometryShapeFunctionContainer<GeometryData::IntegrationMethod>::ShapeFunctionsDerivativesContainerType ShapeFunctionsDerivativesContainerType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor for untrimmed patch
    BrepSurface(
        typename NurbsSurfaceType::Pointer pSurface)
        : BaseType(PointsArrayType(), &msGeometryData)
        , mpNurbsSurface(pSurface)
    {
        mIsTrimmed = false;
    }

    /// Constructor for trimmed patch
    BrepSurface(
        typename NurbsSurfaceType::Pointer pSurface,
        BrepCurveOnSurfaceLoopArrayType& BrepOuterLoopArray,
        BrepCurveOnSurfaceLoopArrayType& BrepInnerLoopArray)
        : BaseType(PointsArrayType(), &msGeometryData)
        , mpNurbsSurface(pSurface)
        , mOuterLoopArray(BrepOuterLoopArray)
        , mInnerLoopArray(BrepInnerLoopArray)
    {
        mIsTrimmed = !(mOuterLoopArray.size() == 0 && mInnerLoopArray.size() == 0);
    }

    /// Constructor for trimmed patch including IsTrimmed
    BrepSurface(
        typename NurbsSurfaceType::Pointer pSurface,
        BrepCurveOnSurfaceLoopArrayType& BrepOuterLoopArray,
        BrepCurveOnSurfaceLoopArrayType& BrepInnerLoopArray,
        bool IsTrimmed)
        : BaseType(PointsArrayType(), &msGeometryData)
        , mpNurbsSurface(pSurface)
        , mOuterLoopArray(BrepOuterLoopArray)
        , mInnerLoopArray(BrepInnerLoopArray)
        , mIsTrimmed(IsTrimmed)
    {
    }

    explicit BrepSurface(const PointsArrayType& ThisPoints)
        : BaseType(ThisPoints, &msGeometryData)
    {
    }

    /// Copy constructor.
    BrepSurface( BrepSurface const& rOther )
        : BaseType( rOther )
        , mpNurbsSurface(rOther.mpNurbsSurface)
        , mOuterLoopArray(rOther.mOuterLoopArray)
        , mInnerLoopArray(rOther.mInnerLoopArray)
        , mEmbeddedEdgesArray(rOther.mEmbeddedEdgesArray)
        , mIsTrimmed(rOther.mIsTrimmed)
    {
    }

    /// Copy constructor from a geometry with different point type.
    template<class TOtherContainerPointType, class TOtherContainerPointEmbeddedType>
    explicit BrepSurface(
        BrepSurface<TOtherContainerPointType, TShiftedBoundary, TOtherContainerPointEmbeddedType> const& rOther )
        : BaseType( rOther )
        , mpNurbsSurface(rOther.mpNurbsSurface)
        , mOuterLoopArray(rOther.mOuterLoopArray)
        , mInnerLoopArray(rOther.mInnerLoopArray)
        , mEmbeddedEdgesArray(rOther.mEmbeddedEdgesArray)
        , mIsTrimmed(rOther.mIsTrimmed)
    {
    }

    /// Destructor
    ~BrepSurface() override = default;

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    BrepSurface& operator=( const BrepSurface& rOther )
    {
        BaseType::operator=( rOther );
        mpNurbsSurface = rOther.mpNurbsSurface;
        mOuterLoopArray = rOther.mOuterLoopArray;
        mInnerLoopArray = rOther.mInnerLoopArray;
        mEmbeddedEdgesArray = rOther.mEmbeddedEdgesArray;
        mIsTrimmed = rOther.mIsTrimmed;
        mpSurrogateInnerLoopGeometries = rOther.mpSurrogateInnerLoopGeometries;
        mpSurrogateOuterLoopGeometries = rOther.mpSurrogateOuterLoopGeometries;
        return *this;
    }

    /// Assignment operator for geometries with different point type.
    template<class TOtherContainerPointType, class TOtherContainerPointEmbeddedType>
    BrepSurface& operator=( BrepSurface<TOtherContainerPointType, TShiftedBoundary, TOtherContainerPointEmbeddedType> const & rOther )
    {
        BaseType::operator=( rOther );
        mpNurbsSurface = rOther.mpNurbsSurface;
        mOuterLoopArray = rOther.mOuterLoopArray;
        mInnerLoopArray = rOther.mInnerLoopArray;
        mEmbeddedEdgesArray = rOther.mEmbeddedEdgesArray;
        mIsTrimmed = rOther.mIsTrimmed;
        mpSurrogateInnerLoopGeometries = rOther.mpSurrogateInnerLoopGeometries;
        mpSurrogateOuterLoopGeometries = rOther.mpSurrogateOuterLoopGeometries;
        return *this;
    }

    ///@}
    ///@name Operations
    ///@{

    typename BaseType::Pointer Create( PointsArrayType const& ThisPoints ) const override
    {
        return typename BaseType::Pointer( new BrepSurface( ThisPoints ) );
    }

    ///@}
    ///@name Access to Geometry Parts
    ///@{

    /**
    * @brief This function returns the pointer of the geometry
    *        which is corresponding to the trim index.
    *        Surface of the geometry is accessible with SURFACE_INDEX.
    * @param Index: trim_index or SURFACE_INDEX.
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
    *        Surface of the geometry is accessible with GeometryType::BACKGROUND_GEOMETRY_INDEX.
    * @param Index: trim_index or GeometryType::BACKGROUND_GEOMETRY_INDEX.
    * @return pointer of geometry, corresponding to the index.
    */
    const GeometryPointer pGetGeometryPart(const IndexType Index) const override
    {
        if (Index == GeometryType::BACKGROUND_GEOMETRY_INDEX)
            return mpNurbsSurface;

        for (IndexType i = 0; i < mOuterLoopArray.size(); ++i)
        {
            for (IndexType j = 0; j < mOuterLoopArray[i].size(); ++j)
            {
                if (mOuterLoopArray[i][j]->Id() == Index)
                    return mOuterLoopArray[i][j];
            }
        }

        for (IndexType i = 0; i < mInnerLoopArray.size(); ++i)
        {
            for (IndexType j = 0; j < mInnerLoopArray[i].size(); ++j)
            {
                if (mInnerLoopArray[i][j]->Id() == Index)
                    return mInnerLoopArray[i][j];
            }
        }

        for (IndexType i = 0; i < mEmbeddedEdgesArray.size(); ++i)
        {
            if (mEmbeddedEdgesArray[i]->Id() == Index)
                return mEmbeddedEdgesArray[i];
        }

        KRATOS_ERROR << "Index " << Index << " not existing in BrepSurface: "
            << this->Id() << std::endl;
    }

    /**
    * @brief This function is used to check if this BrepSurface
    *        has certain trim or surface object.
    * @param Index of the geometry part.
    * @return true if has trim or surface
    */
    bool HasGeometryPart(const IndexType Index) const override
    {
        if (Index == GeometryType::BACKGROUND_GEOMETRY_INDEX)
            return true;

        for (IndexType i = 0; i < mOuterLoopArray.size(); ++i)
        {
            for (IndexType j = 0; j < mOuterLoopArray[i].size(); ++j)
            {
                if (mOuterLoopArray[i][j]->Id() == Index)
                    return true;
            }
        }

        for (IndexType i = 0; i < mInnerLoopArray.size(); ++i)
        {
            for (IndexType j = 0; j < mInnerLoopArray[i].size(); ++j)
            {
                if (mInnerLoopArray[i][j]->Id() == Index)
                    return true;
            }
        }

        for (IndexType i = 0; i < mEmbeddedEdgesArray.size(); ++i)
        {
            if (mEmbeddedEdgesArray[i]->Id() == Index)
                return true;
        }

        return false;
    }

    /// @brief Used to add the embedded edges to the brep surface.
    void AddEmbeddedEdges(BrepCurveOnSurfaceArrayType EmbeddedEdges)
    {
        mEmbeddedEdgesArray = EmbeddedEdges;
    }

    /// Access the nested loop of outer loops.
    const BrepCurveOnSurfaceLoopArrayType& GetOuterLoops() const {
        return mOuterLoopArray;
    }

    /// Access the nested loop of inner loops.
    const BrepCurveOnSurfaceLoopArrayType& GetInnerLoops() const {
        return mInnerLoopArray;
    }

    /// Access the array of embedded edges.
    const BrepCurveOnSurfaceArrayType& GetEmbeddedEdges() const {
        return mEmbeddedEdgesArray;
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
            mpNurbsSurface->Calculate(rVariable, rOutput);
        }
    }

    ///@}
    ///@name Mathematical Informations
    ///@{

    /// Return polynomial degree of the nurbs surface
    SizeType PolynomialDegree(IndexType LocalDirectionIndex) const override
    {
        return mpNurbsSurface->PolynomialDegree(LocalDirectionIndex);
    }

    ///@}
    ///@name Information
    ///@{

    /*
    * @brief checks if the BrepSurface has any boundary trim information.
    * @return true if has no trimming.
    */
    bool IsTrimmed() const
    {
        return mIsTrimmed;
    }

    /// Returns number of points of NurbsSurface.
    SizeType PointsNumberInDirection(IndexType DirectionIndex) const override
    {
        return mpNurbsSurface->PointsNumberInDirection(DirectionIndex);
    }

    ///@}
    ///@name Geometrical Operations
    ///@{

    /// Provides the center of the underlying surface
    Point Center() const override
    {
        return mpNurbsSurface->Center();
    }

    /**
    * @brief Calls projection of its nurbs surface.
    *        Projects a certain point on the geometry, or finds
    *        the closest point, depending on the provided
    *        initial guess. The external point does not necessary
    *        lay on the geometry.
    *        It shall deal as the interface to the mathematical
    *        projection function e.g. the Newton-Raphson.
    *        Thus, the breaking criteria does not necessarily mean
    *        that it found a point on the surface, if it is really
    *        the closest if or not. It shows only if the breaking
    *        criteria, defined by the tolerance is reached.
    *
    *        This function requires an initial guess, provided by
    *        rProjectedPointLocalCoordinates.
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
        return mpNurbsSurface->ProjectionPointGlobalToLocalSpace(
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
        mpNurbsSurface->GlobalCoordinates(rResult, rLocalCoordinates);

        return rResult;
    }

    ///@}
    ///@name Integration Info
    ///@{

    /// Provides the default integration dependent on the polynomial degree.
    IntegrationInfo GetDefaultIntegrationInfo() const override
    {
        return mpNurbsSurface->GetDefaultIntegrationInfo();
    }

    ///@}
    ///@name Integration Points
    ///@{

    /* Creates integration points on the nurbs surface of this geometry.
     * Accounting for whether the surface is trimmed or untrimmed, and whether shifted 
     * boundary conditions are used
     * 
     * - **Untrimmed Surface**: -> Non-cutting case
     *   - If `TShiftedBoundary` is true, the method prepares for the shifted boundary method (SBM) 
     *   - Otherwise, it directly uses `CreateIntegrationPoints` from the underlying NURBS surface.
     * - **Trimmed Surface**: -> Cutting case
     *   - It calls `BrepTrimmingUtilities::CreateBrepSurfaceTrimmingIntegrationPoints` 
     *     to generate integration points that conform to the trimming curve.
     * 
     * @param return integration points.
     */
    void CreateIntegrationPoints(
        IntegrationPointsArrayType& rIntegrationPoints,
        IntegrationInfo& rIntegrationInfo) const override
    {
        if (!mIsTrimmed) {
            // sbm case
            if constexpr (TShiftedBoundary) {

                std::vector<double> spans_u;
                std::vector<double> spans_v;
                mpNurbsSurface->SpansLocalSpace(spans_u, 0);
                mpNurbsSurface->SpansLocalSpace(spans_v, 1);
                
                // Call  "BrepSBMUtilities::CreateBrepSurfaceSBMIntegrationPoints"
                BrepSbmUtilitiesType::CreateBrepSurfaceSbmIntegrationPoints(
                    spans_u, 
                    spans_v,
                    *mpSurrogateOuterLoopGeometries,
                    *mpSurrogateInnerLoopGeometries,
                    rIntegrationPoints,
                    rIntegrationInfo);
            }
            // body-fitted case
            else {
                mpNurbsSurface->CreateIntegrationPoints(
                    rIntegrationPoints, rIntegrationInfo);
            }
        }
        // trimmed case
        else
        {
            std::vector<double> spans_u;
            std::vector<double> spans_v;
            mpNurbsSurface->SpansLocalSpace(spans_u, 0);
            mpNurbsSurface->SpansLocalSpace(spans_v, 1);

            BrepTrimmingUtilitiesType::CreateBrepSurfaceTrimmingIntegrationPoints(
                rIntegrationPoints,
                mOuterLoopArray, mInnerLoopArray,
                spans_u, spans_v,
                rIntegrationInfo);
        }
    }

    ///@}
    ///@name Quadrature Point Geometries
    ///@{

    /* @brief calls function of underlying nurbs surface and updates
     *        the parent to itself.
     *
     * @param rResultGeometries list of quadrature point geometries.
     * @param NumberOfShapeFunctionDerivatives the number of evaluated
     *        derivatives of shape functions at the quadrature point geometries.
     * @param rIntegrationPoints list of provided integration points.
     *
     * @see quadrature_point_geometry.h
     */
    void CreateQuadraturePointGeometries(
        GeometriesArrayType& rResultGeometries,
        IndexType NumberOfShapeFunctionDerivatives,
        const IntegrationPointsArrayType& rIntegrationPoints,
        IntegrationInfo& rIntegrationInfo) override
    {
        mpNurbsSurface->CreateQuadraturePointGeometries(
            rResultGeometries, NumberOfShapeFunctionDerivatives, rIntegrationPoints, rIntegrationInfo);

        for (IndexType i = 0; i < rResultGeometries.size(); ++i) {
            rResultGeometries(i)->SetGeometryParent(this);
        }
    }


    /**
     * @brief Creates a list of quadrature point geometries
     *        from a list of integration points.
     * @param rResultGeometries List of quadrature point geometries.
     * @param NumberOfShapeFunctionDerivatives the number of evaluated
     *        derivatives of shape functions at the quadrature point geometries.
     * @param rIntegrationInfo.
     * @see quadrature_point_geometry.h
     */
    void CreateQuadraturePointGeometries(
        GeometriesArrayType& rResultGeometries,
        IndexType NumberOfShapeFunctionDerivatives,
        IntegrationInfo& rIntegrationInfo) override
    {
        // Option 1: One QuadraturePointGeometry is created containing all integration points. This should be used
        // when all rIntegrationPoints are located inside the same element.
        if(  IntegrationInfo::QuadratureMethod::GAUSS == rIntegrationInfo.GetQuadratureMethod(0) )
        {
            const SizeType points_in_u = mpNurbsSurface->PolynomialDegreeU() + 1;
            const SizeType points_in_v = mpNurbsSurface->PolynomialDegreeV() + 1;

            auto knot_span_intervals_u = mpNurbsSurface->KnotSpanIntervalsU();
            auto knot_span_intervals_v = mpNurbsSurface->KnotSpanIntervalsV();

            const SizeType number_of_elements =
                knot_span_intervals_u.size() * knot_span_intervals_v.size();

            if (rResultGeometries.size() != number_of_elements) {
                rResultGeometries.resize(number_of_elements);
            }

            IntegrationPointsArrayType IntegrationPoints;
            const SizeType number_of_points_per_span =
                points_in_u * points_in_v;

            if (IntegrationPoints.size() != number_of_points_per_span) {
                IntegrationPoints.resize(number_of_points_per_span);
            }

            for (IndexType i = 0; i < knot_span_intervals_u.size(); ++i) {
                for (IndexType j = 0; j < knot_span_intervals_v.size(); ++j) {

                    // KRATOS_WATCH(i)
                    // KRATOS_WATCH(j)

                    GeometriesArrayType result_geometries_per_span;

                    // ////////////////////////////Testing///////////////////////////////////////////
                    // // shape function container.
                    // NurbsSurfaceShapeFunction shape_function_container(
                    //     mpNurbsSurface->PolynomialDegreeU(), mpNurbsSurface->PolynomialDegreeV(), NumberOfShapeFunctionDerivatives);

                    // auto default_method = this->GetDefaultIntegrationMethod();
                    // SizeType num_nonzero_cps = shape_function_container.NumberOfNonzeroControlPoints();
                    // const SizeType num_points = IntegrationPoints.size();
                    // KRATOS_ERROR_IF(num_points < 1) << "List of integration points is empty.\n";

                    // // Initialize containers.
                    // IntegrationPointsContainerType integration_points;
                    // ShapeFunctionsValuesContainerType shape_function_values;
                    // ShapeFunctionsDerivativesContainerType shape_function_derivatives;

                    // integration_points[0] = IntegrationPoints;
                    // shape_function_values[0].resize(IntegrationPoints.size(), num_nonzero_cps);
                    // shape_function_derivatives[0].resize(NumberOfShapeFunctionDerivatives - 1);

                    // for (IndexType i = 0; i < NumberOfShapeFunctionDerivatives - 1; ++i) {
                    //     shape_function_derivatives[0][i].resize(num_points);
                    // }
                    // for (IndexType i = 0; i < NumberOfShapeFunctionDerivatives - 1; ++i) {
                    //     for( IndexType i_point = 0; i_point < num_points; ++i_point){
                    //         shape_function_derivatives[0][i][i_point].resize(num_nonzero_cps, i + 2);
                    //     }
                    // }

                    // for (IndexType i_point = 0; i_point < num_points; ++i_point)
                    // {
                    //     if (mpNurbsSurface->IsRational()) {
                    //         shape_function_container.ComputeNurbsShapeFunctionValues(
                    //             mpNurbsSurface->KnotsU(), mpNurbsSurface->KnotsV(), mpNurbsSurface->Weights(), IntegrationPoints[i_point][0], IntegrationPoints[i_point][1]);
                    //     }
                    //     else {
                    //         shape_function_container.ComputeBSplineShapeFunctionValues(
                    //             mpNurbsSurface->KnotsU(), mpNurbsSurface->KnotsV(), IntegrationPoints[i_point][0], IntegrationPoints[i_point][1]);
                    //     }

                    //     /// Get Shape Functions.
                    //     for (IndexType j = 0; j < num_nonzero_cps; ++j) {
                    //         shape_function_values[0](i_point, j) = shape_function_container(j, 0);
                    //     }

                    //     if (NumberOfShapeFunctionDerivatives > 0) {
                    //         IndexType shape_derivative_index = 1;
                    //         for (IndexType n = 0; n < NumberOfShapeFunctionDerivatives - 1; n++) {
                    //             for (IndexType k = 0; k < n + 2; k++) {
                    //                 for (IndexType j = 0; j < num_nonzero_cps; j++) {
                    //                 shape_function_derivatives[0][n][i_point](j, k) = shape_function_container(j, shape_derivative_index + k);
                    //                 }
                    //             }
                    //             shape_derivative_index += n + 2;
                    //         }
                    //     }
                    // }

                    // /// Get List of Control Points
                    // PointsArrayType nonzero_control_points(num_nonzero_cps);
                    // auto cp_indices = shape_function_container.ControlPointIndices(
                    //     mpNurbsSurface->NumberOfControlPointsU(), mpNurbsSurface->NumberOfControlPointsV());
                    // for (IndexType j = 0; j < num_nonzero_cps; j++) {
                    //     nonzero_control_points(j) = mpNurbsSurface->pGetPoint(cp_indices[j]);
                    // }

                    // KRATOS_WATCH(nonzero_control_points)

                    // // Idea:
                    // // KnotSpanGeometry knot_span_geo(points_in_u, points_in_v,
                    // //     knot_span_intervals_u[i], knot_span_intervals_v[j], mpNurbsSurface);

                    // // knot_span_geo.CreateIntegrationPoints()
                    // // knot_span_geo.CreateQuadraturePoints()

                    // ////////////////////////////Testing///////////////////////////////////////////


                    if (!mIsTrimmed) {
                        // sbm case
                        if constexpr (TShiftedBoundary) {
                            // TODO: Next PR -> Call  "BrepSBMUtilities::CreateBrepSurfaceSBMIntegrationPoints"
                            mpNurbsSurface->CreateIntegrationPoints(
                                IntegrationPoints, points_in_u, points_in_v,
                                knot_span_intervals_u[i], knot_span_intervals_v[j]);
                        }
                        // body-fitted case
                        else {
                            mpNurbsSurface->CreateIntegrationPoints(
                                IntegrationPoints, points_in_u, points_in_v,
                                knot_span_intervals_u[i], knot_span_intervals_v[j]);
                        }
                    }
                    // trimmed case
                    else
                    {
                        BrepTrimmingUtilitiesType::CreateBrepSurfaceTrimmingIntegrationPointsPerKnot(
                            IntegrationPoints,
                            mOuterLoopArray, mInnerLoopArray,
                            knot_span_intervals_u[i], knot_span_intervals_v[j],
                            rIntegrationInfo);
                    }

                    // KRATOS_WATCH(IntegrationPoints.size())

                    mpNurbsSurface->CreateQuadraturePointGeometries(
                        result_geometries_per_span,
                        NumberOfShapeFunctionDerivatives,
                        IntegrationPoints,
                        rIntegrationInfo);

                    IndexType knot_span_index = i* knot_span_intervals_v.size()+ j;

                    rResultGeometries(knot_span_index) = result_geometries_per_span(0);
                }
            }

            for (IndexType i = 0; i < rResultGeometries.size(); ++i) {
                rResultGeometries(i)->SetGeometryParent(this);
            }
        }
        // Option 2: A list of QuadraturePointGeometry is created, one for each integration points.
        else if ( IntegrationInfo::QuadratureMethod::CUSTOM == rIntegrationInfo.GetQuadratureMethod(0) )
        {
            IntegrationPointsArrayType IntegrationPoints;
            CreateIntegrationPoints(IntegrationPoints, rIntegrationInfo);

            this->CreateQuadraturePointGeometries(
                rResultGeometries,
                NumberOfShapeFunctionDerivatives,
                IntegrationPoints,
                rIntegrationInfo);
        }
        else {
            KRATOS_ERROR << "Integration method not available.\n";
        }
    }


    ///@}
    ///@name Shape Function
    ///@{

    Vector& ShapeFunctionsValues(
        Vector &rResult,
        const CoordinatesArrayType& rCoordinates) const override
    {
        mpNurbsSurface->ShapeFunctionsValues(rResult, rCoordinates);

        return rResult;
    }

    Matrix& ShapeFunctionsLocalGradients(
        Matrix& rResult,
        const CoordinatesArrayType& rCoordinates) const override
    {
        mpNurbsSurface->ShapeFunctionsLocalGradients(rResult, rCoordinates);

        return rResult;
    }

    ///@}
    ///@name Geometry Family
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
        return GeometryData::KratosGeometryType::Kratos_Brep_Surface;
    }

    /**
     * @brief Set the Surrogate Outer Loop Geometries object
     * @param pSurrogateOuterLoopArray 
     */
    void SetSurrogateOuterLoopGeometries(Kratos::shared_ptr<GeometrySurrogateArrayType> pSurrogateOuterLoopArray)
    {
        mpSurrogateOuterLoopGeometries = pSurrogateOuterLoopArray;
    }
    
    /**
     * @brief Set the Surrogate Inner Loop Geometries object
     * @param pSurrogateInnerLoopArray 
     */
    void SetSurrogateInnerLoopGeometries(Kratos::shared_ptr<GeometrySurrogateArrayType> pSurrogateInnerLoopArray)
    {
        mpSurrogateInnerLoopGeometries = pSurrogateInnerLoopArray;
    }

    /**
     * @brief Get the Surrogate Inner Loop Geometries object
     * @return GeometrySurrogateArrayType 
     */
    GeometrySurrogateArrayType& GetSurrogateInnerLoopGeometries()
    {
        return *mpSurrogateInnerLoopGeometries;
    }

    /**
     * @brief Get the Surrogate Outer Loop Geometries object
     * @return GeometrySurrogateArrayType 
     */
    GeometrySurrogateArrayType& GetSurrogateOuterLoopGeometries()
    {
        return *mpSurrogateOuterLoopGeometries;
    }

    ///@}
    ///@name Information
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "Brep surface";
    }

    /// Print information about this object.
    void PrintInfo( std::ostream& rOStream ) const override
    {
        rOStream << "Brep surface";
    }

    /// Print object's data.
    void PrintData( std::ostream& rOStream ) const override
    {
        BaseType::PrintData( rOStream );
        std::cout << std::endl;
        rOStream << "    Brep surface " << std::endl;
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

    typename NurbsSurfaceType::Pointer mpNurbsSurface;

    BrepCurveOnSurfaceLoopArrayType mOuterLoopArray;
    BrepCurveOnSurfaceLoopArrayType mInnerLoopArray;

    BrepCurveOnSurfaceArrayType mEmbeddedEdgesArray;

    // For SBM
    Kratos::shared_ptr<GeometrySurrogateArrayType> mpSurrogateOuterLoopGeometries = nullptr;
    Kratos::shared_ptr<GeometrySurrogateArrayType> mpSurrogateInnerLoopGeometries = nullptr;
    
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
        rSerializer.save("NurbsSurface", mpNurbsSurface);
        rSerializer.save("OuterLoopArray", mOuterLoopArray);
        rSerializer.save("InnerLoopArray", mInnerLoopArray);
        rSerializer.save("EmbeddedEdgesArray", mEmbeddedEdgesArray);
        rSerializer.save("IsTrimmed", mIsTrimmed);
        rSerializer.save("SurrogateInnerLoopGeometries", mpSurrogateInnerLoopGeometries);
        rSerializer.save("SurrogateOuterLoopGeometries", mpSurrogateOuterLoopGeometries);
    }

    void load( Serializer& rSerializer ) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BaseType );
        rSerializer.load("NurbsSurface", mpNurbsSurface);
        rSerializer.load("OuterLoopArray", mOuterLoopArray);
        rSerializer.load("InnerLoopArray", mInnerLoopArray);
        rSerializer.load("EmbeddedEdgesArray", mEmbeddedEdgesArray);
        rSerializer.load("IsTrimmed", mIsTrimmed);
        rSerializer.save("SurrogateInnerLoopGeometries", mpSurrogateInnerLoopGeometries);
        rSerializer.save("SurrogateOuterLoopGeometries", mpSurrogateOuterLoopGeometries);
    }

    BrepSurface()
        : BaseType( PointsArrayType(), &msGeometryData )
    {}

    ///@}

}; // Class BrepSurface

///@name Input and output
///@{

/// input stream functions
template<class TContainerPointType, bool TShiftedBoundary, class TContainerPointEmbeddedType = TContainerPointType> inline std::istream& operator >> (
    std::istream& rIStream,
    BrepSurface<TContainerPointType, TShiftedBoundary, TContainerPointEmbeddedType>& rThis );

/// output stream functions
template<class TContainerPointType, bool TShiftedBoundary, class TContainerPointEmbeddedType = TContainerPointType> inline std::ostream& operator << (
    std::ostream& rOStream,
    const BrepSurface<TContainerPointType, TShiftedBoundary, TContainerPointEmbeddedType>& rThis )
{
    rThis.PrintInfo( rOStream );
    rOStream << std::endl;
    rThis.PrintData( rOStream );
    return rOStream;
}

///@}
///@name Static Type Declarations
///@{

template<class TContainerPointType, bool TShiftedBoundary, class TContainerPointEmbeddedType> const
GeometryData BrepSurface<TContainerPointType, TShiftedBoundary, TContainerPointEmbeddedType>::msGeometryData(
    &msGeometryDimension,
    GeometryData::IntegrationMethod::GI_GAUSS_1,
    {}, {}, {});

template<class TContainerPointType, bool TShiftedBoundary, class TContainerPointEmbeddedType>
const GeometryDimension BrepSurface<TContainerPointType, TShiftedBoundary, TContainerPointEmbeddedType>::msGeometryDimension(3, 2);

///@}
}// namespace Kratos.

#endif // KRATOS_BREP_FACE_3D_H_INCLUDED  defined

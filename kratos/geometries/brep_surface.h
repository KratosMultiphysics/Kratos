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

#if !defined(KRATOS_BREP_FACE_3D_H_INCLUDED )
#define  KRATOS_BREP_FACE_3D_H_INCLUDED

// System includes

// External includes
#include "geometries/clipper.hpp"

// Project includes
#include "geometries/geometry.h"
#include "geometries/brep_curve_on_surface.h"
#include "geometries/nurbs_shape_function_utilities/nurbs_interval.h"

// includes for trimming integration
#include "utilities/tessellation_utilities/curve_tessellation.h"

namespace Kratos
{
///@name Kratos Classes
///@{

using namespace ClipperLib;

/**
 * @class BrepSurface
 * @ingroup KratosCore
 * @brief The BrepSurface acts as topology for faces. Those
 *        can be enclosed by a certain set of brep face curves.
 */
template<class TContainerPointType, class TContainerPointEmbeddedType = TContainerPointType>
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

    typedef GeometryData::IntegrationMethod IntegrationMethod;

    typedef NurbsSurfaceGeometry<3, TContainerPointType> NurbsSurfaceType;
    typedef BrepCurveOnSurface<TContainerPointType, TContainerPointEmbeddedType> BrepCurveOnSurfaceType;

    typedef DenseVector<typename BrepCurveOnSurfaceType::Pointer> BrepCurveOnSurfaceLoopType;
    typedef DenseVector<DenseVector<typename BrepCurveOnSurfaceType::Pointer>> BrepCurveOnSurfaceLoopArrayType;

    typedef typename BaseType::GeometriesArrayType GeometriesArrayType;

    typedef typename BaseType::IndexType IndexType;
    typedef typename BaseType::SizeType SizeType;

    typedef typename BaseType::PointsArrayType PointsArrayType;
    typedef typename BaseType::CoordinatesArrayType CoordinatesArrayType;
    typedef typename BaseType::IntegrationPointsArrayType IntegrationPointsArrayType;

    typedef CurveTessellation<PointerVector<typename TContainerPointType::value_type>> CurveTesselationType;
    typedef std::vector<std::pair<double, CoordinatesArrayType>> TessellationType;

    static constexpr IndexType SURFACE_INDEX = -1;

    typedef signed long long cInt;

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
        , mIsTrimmed(rOther.mIsTrimmed)
    {
    }

    /// Copy constructor from a geometry with different point type.
    template<class TOtherContainerPointType, class TOtherContainerPointEmbeddedType>
    explicit BrepSurface(
        BrepSurface<TOtherContainerPointType, TOtherContainerPointEmbeddedType> const& rOther )
        : BaseType( rOther )
        , mpNurbsSurface(rOther.mpNurbsSurface)
        , mOuterLoopArray(rOther.mOuterLoopArray)
        , mInnerLoopArray(rOther.mInnerLoopArray)
        , mIsTrimmed(rOther.mIsTrimmed)
    {
    }

    /// Destructor
    ~BrepSurface() override = default;

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator
    BrepSurface& operator=( const BrepSurface& rOther )
    {
        BaseType::operator=( rOther );
        mpNurbsSurface = rOther.mpNurbsSurface;
        mOuterLoopArray = rOther.mOuterLoopArray;
        mInnerLoopArray = rOther.mInnerLoopArray;
        mIsTrimmed = rOther.mIsTrimmed;
        return *this;
    }

    /// Assignment operator with different point type
    template<class TOtherContainerPointType, class TOtherContainerPointEmbeddedType>
    BrepSurface& operator=( BrepSurface<TOtherContainerPointType, TOtherContainerPointEmbeddedType> const & rOther )
    {
        BaseType::operator=( rOther );
        mpNurbsSurface = rOther.mpNurbsSurface;
        mOuterLoopArray = rOther.mOuterLoopArray;
        mInnerLoopArray = rOther.mInnerLoopArray;
        mIsTrimmed = rOther.mIsTrimmed;
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
    ///@name Point Access
    ///@{

    PointType& operator[](const SizeType& i) override
    {
        return (*mpNurbsSurface)[i];
    }

    PointType const& operator[](const SizeType& i) const override
    {
        return (*mpNurbsSurface)[i];
    }

    typename PointType::Pointer& operator()(const SizeType& i) override
    {
        return (*mpNurbsSurface)(i);
    }

    const typename PointType::Pointer& operator()(const SizeType& i) const override
    {
        return (*mpNurbsSurface)(i);
    }

    SizeType size() const override
    {
        return mpNurbsSurface->size();
    }

    ///@}
    ///@name Access to Geometry Parts
    ///@{

    /**
    * @brief This function returns the pointer of the geometry
    *        which is corresponding to the trim index.
    *        Surface of the geometry is accessable with SURFACE_INDEX.
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
    *        Surface of the geometry is accessable with SURFACE_INDEX.
    * @param Index: trim_index or SURFACE_INDEX.
    * @return pointer of geometry, corresponding to the index.
    */
    const GeometryPointer pGetGeometryPart(const IndexType Index) const override
    {
        if (Index == SURFACE_INDEX)
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
        if (Index == SURFACE_INDEX)
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

        return false;
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

    /* @brief Provides Spans of the underlying nurbs surface.
     * @param vector of span intervals.
     * @param index of chosen direction, possible direction 0 and 1.
     */
    void Spans(std::vector<double>& rSpans, IndexType DirectionIndex) const override
    {
        mpNurbsSurface->Spans(rSpans, DirectionIndex);
    }

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
    * @param rProjectedPointGlobalCoordinates the location of the
    *        projection in global coordinates.
    * @param rProjectedPointLocalCoordinates the location of the
    *        projection in local coordinates.
    *        The variable is as initial guess!
    * @param Tolerance accepted of orthogonal error to projection.
    * @return It is chosen to take an int as output parameter to
    *         keep more possibilities within the interface.
    *         0 -> failed
    *         1 -> converged
    */
    int ProjectionPoint(
        const CoordinatesArrayType& rPointGlobalCoordinates,
        CoordinatesArrayType& rProjectedPointGlobalCoordinates,
        CoordinatesArrayType& rProjectedPointLocalCoordinates,
        const double Tolerance = std::numeric_limits<double>::epsilon()
        ) const override
    {
        return mpNurbsSurface->ProjectionPoint(
            rPointGlobalCoordinates, rProjectedPointGlobalCoordinates, rProjectedPointLocalCoordinates, Tolerance);
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
    ///@name Integration Points
    ///@{

    /* Creates integration points on the nurbs surface of this geometry.
     * @param return integration points.
     */
    void CreateIntegrationPoints(
        IntegrationPointsArrayType& rIntegrationPoints,
        IntegrationInfo& rIntegrationInfo) const override
    {
        Paths outer_loops(mOuterLoopArray.size()), inner_loops(mInnerLoopArray.size()), solution;
        const double factor = 1e-10; 

        for (IndexType i = 0; i < mOuterLoopArray.size(); ++i) {
            for (IndexType j = 0; j < mOuterLoopArray[i].size(); ++j) {               
                CurveTesselationType curve_tesselation;
                auto geometry_outer = *(mOuterLoopArray[i][j].get());
                curve_tesselation.Tessellate(
                    geometry_outer, 0.01, 1, true);
                auto tesselation = curve_tesselation.GetTessellation();
                for (IndexType u = 0; u < tesselation.size(); ++u) {
                    outer_loops[i] << ToIntPoint(std::get<1>(tesselation[u])[0], std::get<1>(tesselation[u])[1], factor);
                }
            }
        }

        for (IndexType i = 0; i < mInnerLoopArray.size(); ++i) {
            for (IndexType j = 0; j < mInnerLoopArray[i].size(); ++j) {
                CurveTesselationType curve_tesselation;
                auto geometry_inner = *(mInnerLoopArray[i][j].get());
                curve_tesselation.Tessellate(
                    geometry_inner, 0.01, 1);
                auto tesselation = curve_tesselation.GetTessellation();
                for (IndexType u = 0; u < tesselation.size(); ++u) {
                    inner_loops[i] << ToIntPoint(std::get<1>(tesselation[u])[0], std::get<1>(tesselation[u])[1], factor);
                }
            }
        }

        std::vector<double> spans_u;
        std::vector<double> spans_v;
        mpNurbsSurface->Spans(spans_u, 0);
        mpNurbsSurface->Spans(spans_v, 1);
        
        //perform intersection
        Clipper c;

        auto knot_span_intervals_u = mpNurbsSurface->KnotSpanIntervalsU();
        auto knot_span_intervals_v = mpNurbsSurface->KnotSpanIntervalsV();

        //c.Execute(ctIntersection, solution, pftNonZero, pftNonZero);
        for (IndexType i = 0; i < knot_span_intervals_u.size(); ++i) {
            for (IndexType j = 0; j < knot_span_intervals_v.size(); ++j) {

                c.AddPaths(outer_loops, ptSubject, true);
                c.AddPaths(inner_loops, ptSubject, true);

                Paths span(1);
                span[0] <<
                    ToIntPoint(knot_span_intervals_u[i].GetT0(), knot_span_intervals_v[j].GetT0(), factor) << ToIntPoint(knot_span_intervals_u[i].GetT1(), knot_span_intervals_v[j].GetT0(), factor) <<
                    ToIntPoint(knot_span_intervals_u[i].GetT1(), knot_span_intervals_v[j].GetT1(), factor) << ToIntPoint(knot_span_intervals_u[i].GetT0(), knot_span_intervals_v[j].GetT1(), factor);

                c.AddPaths(span, ptClip, true);
                c.Execute(ctIntersection, solution, pftNonZero, pftNonZero);

                const SizeType points_in_u = mpNurbsSurface->PolynomialDegreeU() + 1;
                const SizeType points_in_v = mpNurbsSurface->PolynomialDegreeV() + 1;

                if (solution.size() == 0) {

                    continue;
                }
                else if (std::abs((std::abs(Area(solution[0])) - std::abs(Area(span[0])))*(factor*factor)) < 1e-6) {

                    const SizeType number_of_integration_points =
                        points_in_u * points_in_v;

                    SizeType initial_integration_size = rIntegrationPoints.size();

                    if (rIntegrationPoints.size() != initial_integration_size+number_of_integration_points) {
                        rIntegrationPoints.resize(initial_integration_size+number_of_integration_points);
                    }

                    typename IntegrationPointsArrayType::iterator integration_point_iterator = rIntegrationPoints.begin();
                    advance(integration_point_iterator, initial_integration_size);
                    

                    IntegrationPointUtilities::IntegrationPoints2D(
                        integration_point_iterator,
                        points_in_u, points_in_v,
                        mpNurbsSurface->KnotSpanIntervalsU()[i].GetT0(), mpNurbsSurface->KnotSpanIntervalsU()[i].GetT1(),
                        mpNurbsSurface->KnotSpanIntervalsV()[j].GetT0(), mpNurbsSurface->KnotSpanIntervalsV()[j].GetT1());
                }
                else {

                    std::vector<Matrix> triangles;
                    Triangulate_OPT(solution[0], triangles, factor);
                    
                    int number_of_triangles = triangles.size();
                    int number_of_points = std::max(points_in_u, points_in_v)+1;
                    int number_of_gauss_points;
                    if (number_of_points == 1) {number_of_gauss_points = 1;}
                    else if (number_of_points == 2) {number_of_gauss_points = 3;}
                    else if (number_of_points == 3) {number_of_gauss_points = 4;}
                    else if (number_of_points == 4) {number_of_gauss_points = 6;}
                    else if (number_of_points == 5) {number_of_gauss_points = 7;}
                    else if (number_of_points == 6) {number_of_gauss_points = 12;}
                    else if (number_of_points == 7) {number_of_gauss_points = 13;}

                    const SizeType number_of_integration_points = number_of_triangles * number_of_gauss_points;

                    SizeType initial_integration_size = rIntegrationPoints.size();

                    if (rIntegrationPoints.size() != initial_integration_size+number_of_integration_points) {
                        rIntegrationPoints.resize(initial_integration_size+number_of_integration_points);
                    }

                    typename IntegrationPointsArrayType::iterator integration_point_iterator = rIntegrationPoints.begin();
                    advance(integration_point_iterator, initial_integration_size);
                    
                    for (IndexType i = 0; i < triangles.size(); ++i)
                    {
                        IntegrationPointUtilities::IntegrationPointsTriangle(
                        integration_point_iterator,
                        number_of_points,
                        triangles[i]);
                    }              
                }
                c.Clear();
            }
        }
        SizeType final_integration_size = rIntegrationPoints.size();
    }

    ///@}
    ///@name Quadrature Point Geometries
    ///@{

    /* @brief calls function of undelying nurbs surface and updates
     *        the parent to itself.
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
        mpNurbsSurface->CreateQuadraturePointGeometries(
            rResultGeometries, NumberOfShapeFunctionDerivatives, rIntegrationPoints);

        for (IndexType i = 0; i < rResultGeometries.size(); ++i) {
            rResultGeometries(i)->SetGeometryParent(this);
        }
    }
    struct DPState {
        bool visible;
        double weight;
        long bestvertex;
    };
    struct Diagonal {
        long index1;
        long index2;
    };

    //Triangulation 
    void Triangulate_OPT(const Path &polygon, std::vector<Matrix>& triangles, const double factor) const
    {
        array_1d<double, 2> p1, p2, p3, p4;
        int bestvertex;
        double weight, minweight, d1, d2;
        Diagonal diagonal, newdiagonal;
        std::list<Diagonal> diagonals;
        bool ret = true;

        int n = polygon.size();
        std::vector< IntPoint > const& points = polygon;
        matrix<DPState> dpstates(n, n);

        //init states and visibility
        for (unsigned int i = 0; i<(n - 1); i++) {
            p1 = FromIntPoint(points[i], factor);
            for (unsigned int j = i + 1; j<n; j++) {
                dpstates(j, i).visible = true;
                dpstates(j, i).weight = 0;
                dpstates(j, i).bestvertex = -1;
                if (j != (i + 1)) {
                    p2 = FromIntPoint(points[j], factor);

                    //visibility check
                    if (i == 0) p3 = FromIntPoint(points[n - 1], factor);
                        else p3 = FromIntPoint(points[i - 1], factor);
                    if (i == (n - 1)) p4 = FromIntPoint(points[0], factor);
                        else p4 = FromIntPoint(points[i + 1], factor);
                    if (!InCone(p3, p1, p4, p2)) {
                        dpstates(j, i).visible = false;
                        continue;
                    }

                    if (j == 0) p3 = FromIntPoint(points[n - 1], factor);
                        else p3 = FromIntPoint(points[j - 1], factor);
                    if (j == (n - 1)) p4 = FromIntPoint(points[0], factor);
                        else p4 = FromIntPoint(points[j + 1], factor);
                    if (!InCone(p3, p2, p4, p1)) {
                        dpstates(j, i).visible = false;
                        continue;
                    }

                    for (unsigned int k = 0; k<n; k++) {
                        p3 = FromIntPoint(points[k], factor);
                        if (k == (n - 1)) p4 = FromIntPoint(points[0], factor);
                            else p4 = FromIntPoint(points[k + 1], factor);
                        if (Intersects(p1, p2, p3, p4)) {
                            dpstates(j, i).visible = false;
                            break;
                        }
                    }
                }
            }
        }

        dpstates(n - 1, 0).visible = true;
        dpstates(n - 1, 0).weight = 0;
        dpstates(n - 1, 0).bestvertex = -1;

        for (unsigned int gap = 2; gap<n; gap++) {
            for (unsigned int i = 0; i<(n - gap); i++) {
                int j = i + gap;
                if (!dpstates(j, i).visible) continue;
                    bestvertex = -1;
                for (unsigned int k = (i + 1); k<j; k++) {
                    if (!dpstates(k, i).visible) continue;
                    if (!dpstates(j, k).visible) continue;

                    if (k <= (i + 1)) d1 = 0;
                        else d1 = Distance(FromIntPoint(points[i], factor), FromIntPoint(points[k], factor));
                    if (j <= (k + 1)) d2 = 0;
                        else d2 = Distance(FromIntPoint(points[k], factor), FromIntPoint(points[j], factor));

                    weight = dpstates(k, i).weight + dpstates(j, k).weight + d1 + d2;

                    if ((bestvertex == -1) || (weight<minweight)) {
                        bestvertex = k;
                        minweight = weight;
                    }
                }
                if (bestvertex == -1) {
                    KRATOS_THROW_ERROR(std::runtime_error, "Triangulate: No points in polygon.", std::endl);
                }

                dpstates(j, i).bestvertex = bestvertex;
                dpstates(j, i).weight = minweight;
            }
        }

        newdiagonal.index1 = 0;
        newdiagonal.index2 = n - 1;
        diagonals.push_back(newdiagonal);

        while (!diagonals.empty()) {
            diagonal = *(diagonals.begin());
            diagonals.pop_front();

            bestvertex = dpstates(diagonal.index2, diagonal.index1).bestvertex;
            if (bestvertex == -1) {
                ret = false;
                break;
            }
            Matrix triangle(3, 2);
            triangle(0, 0) = FromIntPoint(points[diagonal.index1], factor)[0];
            triangle(0, 1) = FromIntPoint(points[diagonal.index1], factor)[1];
            triangle(1, 0) = FromIntPoint(points[bestvertex], factor)[0];
            triangle(1, 1) = FromIntPoint(points[bestvertex], factor)[1];
            triangle(2, 0) = FromIntPoint(points[diagonal.index2], factor)[0];
            triangle(2, 1) = FromIntPoint(points[diagonal.index2], factor)[1];

            if (GetAreaOfTriangle(triangle)>1e-5)
                triangles.push_back(triangle);
            else
            {
                std::cout << "triangle with zero area" << GetAreaOfTriangle(triangle) << std::endl;
            }
            if (bestvertex > (diagonal.index1 + 1)) {
                newdiagonal.index1 = diagonal.index1;
                newdiagonal.index2 = bestvertex;
                diagonals.push_back(newdiagonal);
            }
            if (diagonal.index2 > (bestvertex + 1)) {
                newdiagonal.index1 = bestvertex;
                newdiagonal.index2 = diagonal.index2;
                diagonals.push_back(newdiagonal);
            }
        }
    }

    bool InCone(array_1d<double, 2> &p1, array_1d<double, 2> &p2,
        array_1d<double, 2> &p3, array_1d<double, 2> &p) const
    {
        if (IsConvex(p1, p2, p3)) {
            if (!IsConvex(p1, p2, p)) return false;
            if (!IsConvex(p2, p3, p)) return false;
            return true;
        }
        else {
            if (IsConvex(p1, p2, p)) return true;
            if (IsConvex(p2, p3, p)) return true;
        return false;
        }
    }

    bool IsConvex(
        const array_1d<double, 2>& p1, const array_1d<double, 2>& p2, 
        const array_1d<double, 2>& p3) const
    {
        double tmp;
        tmp = (p3[1] - p1[1])*(p2[0] - p1[0]) - (p3[0] - p1[0])*(p2[1] - p1[1]);
        if (tmp>0) return true;
        else return false;
    }

    double Distance(array_1d<double, 2> point_1, array_1d<double, 2> point_2) const
    {
        return sqrt(point_1[0] * point_2[0] + point_1[1] * point_2[1]);
    }

    double GetAreaOfTriangle(const Matrix& triangle) const
    {
        double area = (triangle(0, 0)*(triangle(1, 1) - triangle(2, 1))
        + triangle(1, 0)*(triangle(2, 1) - triangle(0, 1))
        + triangle(2, 0)*(triangle(0, 1) - triangle(1, 1))) / 2;

        return area;
    }

    //checks if two lines intersect
    bool Intersects(array_1d<double, 2> &p11, array_1d<double, 2> &p12,
        array_1d<double, 2> &p21, array_1d<double, 2> &p22) const
    {
        if ((p11[0] == p21[0]) && (p11[1] == p21[1])) return false;
        if ((p11[0] == p22[0]) && (p11[1] == p22[1])) return false;
        if ((p12[0] == p21[0]) && (p12[1] == p21[1])) return false;
        if ((p12[0] == p22[0]) && (p12[1] == p22[1])) return false;

        array_1d<double, 2> v1ort, v2ort, v;
        double dot11, dot12, dot21, dot22;

        v1ort[0] = p12[1] - p11[1];
        v1ort[1] = p11[0] - p12[0];

        v2ort[0] = p22[1] - p21[1];
        v2ort[1] = p21[0] - p22[0];

        v[0] = p21[0] - p11[0];
        v[1] = p21[1] - p11[1];
        dot21 = v[0] * v1ort[0] + v[1] * v1ort[1];
        v[0] = p22[0] - p11[0];
        v[1] = p22[1] - p11[1];
        dot22 = v[0] * v1ort[0] + v[1] * v1ort[1];

        v[0] = p11[0] - p21[0];
        v[1] = p11[1] - p21[1];
        dot11 = v[0] * v2ort[0] + v[1] * v2ort[1];
        v[0] = p12[0] - p21[0];
        v[1] = p12[1] - p21[1];
        dot12 = v[0] * v2ort[0] + v[1] * v2ort[1];

        if (dot11*dot12>0) return false;
        if (dot21*dot22>0) return false;

        return true;
    }

    inline ClipperLib::IntPoint
    ToIntPoint(
        const double x,
        const double y,
        const double factor) const
    {
        ClipperLib::IntPoint intPoint;

        intPoint.X = static_cast<cInt>(x / factor);
        intPoint.Y = static_cast<cInt>(y / factor);

        return intPoint;
    }

    inline array_1d<double, 2>
    FromIntPoint(
        const ClipperLib::IntPoint& intPoint,
        const double factor) const
    {
        array_1d<double, 2> point;

        point[0] = intPoint.X * factor;
        point[1] = intPoint.Y * factor;

        return point;
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
        rSerializer.save("IsTrimmed", mIsTrimmed);
    }

    void load( Serializer& rSerializer ) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BaseType );
        rSerializer.load("NurbsSurface", mpNurbsSurface);
        rSerializer.load("OuterLoopArray", mOuterLoopArray);
        rSerializer.load("InnerLoopArray", mInnerLoopArray);
        rSerializer.load("IsTrimmed", mIsTrimmed);
    }

    BrepSurface()
        : BaseType( PointsArrayType(), &msGeometryData )
    {}

    ///@}

}; // Class BrepSurface

///@name Input and output
///@{

/// input stream functions
template<class TContainerPointType, class TContainerPointEmbeddedType = TContainerPointType> inline std::istream& operator >> (
    std::istream& rIStream,
    BrepSurface<TContainerPointType, TContainerPointEmbeddedType>& rThis );

/// output stream functions
template<class TContainerPointType, class TContainerPointEmbeddedType = TContainerPointType> inline std::ostream& operator << (
    std::ostream& rOStream,
    const BrepSurface<TContainerPointType, TContainerPointEmbeddedType>& rThis )
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
GeometryData BrepSurface<TContainerPointType, TContainerPointEmbeddedType>::msGeometryData(
    &msGeometryDimension,
    GeometryData::GI_GAUSS_1,
    {}, {}, {});

template<class TContainerPointType, class TContainerPointEmbeddedType>
const GeometryDimension BrepSurface<TContainerPointType, TContainerPointEmbeddedType>::msGeometryDimension(
    2, 3, 2);

///@}
}// namespace Kratos.

#endif // KRATOS_BREP_FACE_3D_H_INCLUDED  defined

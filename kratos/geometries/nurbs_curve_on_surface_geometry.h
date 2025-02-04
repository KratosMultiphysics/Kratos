//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Thomas Oberbichler
//                   Tobias Teschemacher
//                   Andreas Apostolatos
//
//  Ported from the ANurbs library (https://github.com/oberbichler/ANurbs)
//

#if !defined(KRATOS_NURBS_CURVE_ON_SURFACE_H_INCLUDED )
#define  KRATOS_NURBS_CURVE_ON_SURFACE_H_INCLUDED

// Project includes
#include "geometries/geometry.h"

#include "geometries/nurbs_curve_geometry.h"
#include "geometries/nurbs_surface_geometry.h"
#include "geometries/nurbs_shape_function_utilities/nurbs_curve_shape_functions.h"
#include "geometries/nurbs_shape_function_utilities/nurbs_interval.h"

#include "utilities/curve_axis_intersection.h"

#include "utilities/nurbs_utilities/projection_nurbs_geometry_utilities.h"

namespace Kratos {

template <int TWorkingSpaceDimension, class TCurveContainerPointType, class TSurfaceContainerPointType>
class NurbsCurveOnSurfaceGeometry : public Geometry<typename TSurfaceContainerPointType::value_type>
{
public:
    ///@name Type Definitions
    ///@{

    typedef typename TSurfaceContainerPointType::value_type NodeType;
    typedef typename TCurveContainerPointType::value_type CurveNodeType;

    typedef Geometry<NodeType> GeometryType;
    typedef typename GeometryType::Pointer GeometryPointer;

    typedef Geometry<NodeType> BaseType;

    typedef typename BaseType::IndexType IndexType;
    typedef typename BaseType::SizeType SizeType;

    typedef NurbsSurfaceGeometry<3, TSurfaceContainerPointType> NurbsSurfaceType;
    typedef NurbsCurveGeometry<2, TCurveContainerPointType> NurbsCurveType;
    typedef typename NurbsCurveType::Pointer NurbsCurvePointerType;

    typedef typename BaseType::CoordinatesArrayType CoordinatesArrayType;
    typedef typename BaseType::PointsArrayType PointsArrayType;
    typedef typename BaseType::GeometriesArrayType GeometriesArrayType;
    typedef typename BaseType::IntegrationPointsArrayType IntegrationPointsArrayType;

    // using base class functionalities.
    using BaseType::CreateQuadraturePointGeometries;
    using BaseType::pGetPoint;
    using BaseType::GetPoint;

    /// Counted pointer of NurbsCurveOnSurfaceGeometry
    KRATOS_CLASS_POINTER_DEFINITION(NurbsCurveOnSurfaceGeometry);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    NurbsCurveOnSurfaceGeometry(
        typename NurbsSurfaceType::Pointer pSurface,
        typename NurbsCurveType::Pointer pCurve)
        : BaseType(PointsArrayType(), &msGeometryData)
        , mpNurbsSurface(pSurface)
        , mpNurbsCurve(pCurve)
    {
    }

    /// Default constructor
    NurbsCurveOnSurfaceGeometry()
        : BaseType(PointsArrayType(), &msGeometryData)
    {};

    /// Copy constructor
    NurbsCurveOnSurfaceGeometry(NurbsCurveOnSurfaceGeometry const& rOther)
        : BaseType(rOther)
        , mpNurbsSurface(rOther.mpNurbsSurface)
        , mpNurbsCurve(rOther.mpNurbsCurve)
    {
    }

    /// Copy constructor, with different point type.
    template<class TOtherCurveContainerPointType, class TOtherSurfaceContainerPointType> NurbsCurveOnSurfaceGeometry(
        NurbsCurveOnSurfaceGeometry<TWorkingSpaceDimension, TOtherCurveContainerPointType, TOtherSurfaceContainerPointType> const& rOther)
        : BaseType(rOther, &msGeometryData)
        , mpNurbsSurface(rOther.mpNurbsSurface)
        , mpNurbsCurve(rOther.mpNurbsCurve)
    {
    }

    /// Destructor
    ~NurbsCurveOnSurfaceGeometry() override = default;

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator
    NurbsCurveOnSurfaceGeometry& operator=(const NurbsCurveOnSurfaceGeometry& rOther)
    {
        BaseType::operator=(rOther);
        mpNurbsSurface = rOther.mpNurbsSurface;
        mpNurbsCurve = rOther.mpNurbsCurve;
        return *this;
    }

    /// @brief Assignment operator for geometries with different point type.
    template<class TOtherCurveContainerPointType, class TOtherSurfaceContainerPointType>
    NurbsCurveOnSurfaceGeometry& operator=(
        NurbsCurveOnSurfaceGeometry<TWorkingSpaceDimension, TOtherCurveContainerPointType, TOtherSurfaceContainerPointType> const & rOther)
    {
        BaseType::operator=(rOther);
        mpNurbsSurface = rOther.mpNurbsSurface;
        mpNurbsCurve = rOther.mpNurbsCurve;
        return *this;
    }

    ///@}
    ///@name Operations
    ///@{

    typename BaseType::Pointer Create(
        TSurfaceContainerPointType const& ThisPoints) const override
    {
        KRATOS_ERROR << "NurbsCurveOnSurfaceGeometry cannot be created with 'PointsArrayType const& ThisPoints'. "
            << "'Create' is not allowed as it would not contain the required pointers to the surface and the curve."
            << std::endl;
    }

    ///@}
    ///@name Access to Geometry Parts
    ///@{

    /**
    * @brief This function returns the pointer of the geometry
    *        which is corresponding to the trim index.
    *        Possible indices are:
    *        SURFACE_INDEX or EMBEDDED_CURVE_INDEX.
    * @param Index: SURFACE_INDEX or EMBEDDED_CURVE_INDEX.
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
    *        Possible indices are:
    *        SURFACE_INDEX or EMBEDDED_CURVE_INDEX.
    * @param Index: SURFACE_INDEX or EMBEDDED_CURVE_INDEX.
    * @return pointer of geometry, corresponding to the index.
    */
    const GeometryPointer pGetGeometryPart(const IndexType Index) const override
    {
        if (Index == GeometryType::BACKGROUND_GEOMETRY_INDEX)
            return mpNurbsSurface;

        KRATOS_ERROR << "Index " << Index << " not existing in NurbsCurveOnSurface: "
            << this->Id() << std::endl;
    }

    /**
    * @brief This function is used to check if the index is
    *        GeometryType::BACKGROUND_GEOMETRY_INDEX.
    * @param Index of the geometry part.
    * @return true if GeometryType::BACKGROUND_GEOMETRY_INDEX.
    */
    bool HasGeometryPart(const IndexType Index) const override
    {
        return Index == GeometryType::BACKGROUND_GEOMETRY_INDEX;
    }

    ///@}
    ///@name Set / Calculate access
    ///@{

    void Calculate(
        const Variable<array_1d<double, 3>>& rVariable,
        array_1d<double, 3>& rOutput
        ) const override
    {
        if (rVariable == PARAMETER_2D_COORDINATES) {
            array_1d<double, 3>& local_parameter = rOutput;
            mpNurbsCurve->GlobalCoordinates(rOutput, local_parameter);
        }
    }

    ///@}
    ///@name Mathematical Information
    ///@{

    /// Return polynomial degree of the curve
    SizeType PolynomialDegree(IndexType LocalDirectionIndex) const override
    {
        return mpNurbsSurface->PolynomialDegree(0) + mpNurbsSurface->PolynomialDegree(1);
    }

    ///@}
    ///@name Set/ Get functions
    ///@{
    /// Returns the NurbsCurve::Pointer of this CurveOnSurface.
    NurbsCurvePointerType pGetCurve()
    {
        return mpNurbsCurve;
    }

    /// Returns the const NurbsCurveOnSurface::Pointer of this brep.
    const NurbsCurvePointerType pGetCurve() const
    {
        return mpNurbsCurve;
    }

    ///@}
    ///@name Curve Properties
    ///@{

    /// Returns number of points of NurbsCurve.
    SizeType PointsNumberInDirection(IndexType LocalDirectionIndex) const override
    {
        return mpNurbsCurve->PointsNumberInDirection(LocalDirectionIndex);
    }

    /* @brief Provides intersections of the nurbs curve with the knots of the surface,
     *         using the interval of this curve.
     * @param vector of span intervals.
     * @param index of chosen direction, for curves always 0.
     */
    void SpansLocalSpace(
        std::vector<double>& rSpans, 
        IndexType DirectionIndex = 0) const override
    {
        auto interval = mpNurbsCurve->DomainInterval();
        this->SpansLocalSpace(rSpans, interval.GetT0(), interval.GetT1());
    }

    /* @brief Provides intersections of the nurbs curve with the knots of the surface,
     *         using the interval of this curve in the SBM case.
     * @param vector of span intervals.
     * @param index of chosen direction, for curves always 0.
     */
    void SpansLocalSpaceSBM(
        std::vector<double>& rSpans, 
        IndexType DirectionIndex = 0) const
    {
        auto interval = mpNurbsCurve->DomainInterval();
        this->SpansLocalSpaceSBM(rSpans, interval.GetT0(), interval.GetT1());
    }

    /* @brief  Provides intersections of the nurbs curve with the knots of the surface.
     * @return vector of interval limitations.
     */
    void SpansLocalSpace(
        std::vector<double>& rSpans,
        double Start, 
        double End) const
    {
        std::vector<double> surface_spans_u;
        std::vector<double> surface_spans_v;
        mpNurbsSurface->SpansLocalSpace(surface_spans_u, 0);
        mpNurbsSurface->SpansLocalSpace(surface_spans_v, 1);

        // trimming case: compute also inner knot spans intersections
        CurveAxisIntersection<CurveNodeType>::ComputeAxisIntersection(
            rSpans,
            *(mpNurbsCurve.get()), Start, End,
            surface_spans_u, surface_spans_v,
            1e-6);
    }

    /* @brief  Provides intersections of the nurbs curve with the knots of the surface in the SBM case.
     * @return vector of interval limitations.
     */
    void SpansLocalSpaceSBM(
        std::vector<double>& rSpans,
        double Start, 
        double End) const
    {
        std::vector<double> surface_spans_u;
        std::vector<double> surface_spans_v;
        mpNurbsSurface->SpansLocalSpace(surface_spans_u, 0);
        mpNurbsSurface->SpansLocalSpace(surface_spans_v, 1);

        // SBM case: compute only knot spans intersections
        ComputeAxisIntersectionSBM(
            rSpans,
            Start, End,
            surface_spans_u, surface_spans_v);
    }
    
     /**
     * @brief Computes the knot span intersections for the shifted boundary method (SBM) case.
     * This function calculates the intersections between a NURBS curve segment and the knot spans
     * of a surface in the local parameter space, specifically for SBM-based geometries. The computed
     * intersections are returned in the `rSpans` vector.
     * 
     * * @details 
     * The function operates as follows:
     * - Computes the physical coordinates of the start (`Start`) and end (`End`) points of the curve segment.
     * - Determines the segment's orientation (horizontal or vertical) based on the difference in coordinates.
     * - For each direction (U or V):
     *   - Calculates the relevant knot interval.
     *   - Iterates over the surface spans in the corresponding direction, identifying intersections within the interval.
     *   - Maps the intersection values back to the curve's parameter space and stores them in `rSpans`.
     * 
     * @param rSpans 
     * @param Start 
     * @param End 
     * @param surface_spans_u 
     * @param surface_spans_v 
     * @param Tolerance 
     * @return int 
     */
   void ComputeAxisIntersectionSBM(
            std::vector<double>& rSpans,
            double Start,
            double End,
            const std::vector<double>& surface_spans_u,
            const std::vector<double>& surface_spans_v) const 
    {
        CoordinatesArrayType physical_coord_1 = ZeroVector(3);
        CoordinatesArrayType local_coord_1 = ZeroVector(3);
        local_coord_1[0] = Start;
        mpNurbsCurve->GlobalCoordinates(physical_coord_1, local_coord_1);

        CoordinatesArrayType physical_coord_2 = ZeroVector(3);
        CoordinatesArrayType local_coord_2 = ZeroVector(3);
        local_coord_2[0] = End;
        mpNurbsCurve->GlobalCoordinates(physical_coord_2, local_coord_2);

        // 2D need to re-order
        bool reverse_order = false;
        if ((physical_coord_2[0] - physical_coord_1[0] < 0 ) || 
            (physical_coord_2[1] - physical_coord_1[1] < 0 )) {
            reverse_order = true;
        }
        std::vector<double> tempSpans;
        const double tolerance_orientation = 1e-10;
        const double tolerance_intersection = 1e-15;
        
        // Scale factor between surface coordinate and curve one
        double physical_length_segment = norm_2(physical_coord_1-physical_coord_2);
        double parameter_length_segment = norm_2(local_coord_1-local_coord_2);
        double scale_factor = parameter_length_segment/physical_length_segment;

        std::vector<double> knot_interval(2);

        // Compute constant coordinate -> understand segment orientation (vertical/horizontal)
        if (std::abs(physical_coord_1[0]-physical_coord_2[0]) > tolerance_orientation) {
            // horizontal case
            knot_interval[0] = physical_coord_1[0]; 
            knot_interval[1] = physical_coord_2[0];
            std::sort(knot_interval.begin(), knot_interval.end());

            knot_interval[0] -= tolerance_intersection; 
            knot_interval[1] += tolerance_intersection;
            // Compare with surface_spans_u
            for (IndexType i = 0; i < surface_spans_u.size(); i++) {
                double curr_knot_value = surface_spans_u[i];
                if (curr_knot_value < knot_interval[0]) {continue;}
                if (std::abs(curr_knot_value - knot_interval[0]) < tolerance_intersection*10) knot_interval[0] = curr_knot_value;
                if (curr_knot_value > knot_interval[1]) {break;}
                if (std::abs(curr_knot_value - knot_interval[1]) < tolerance_intersection*10) knot_interval[1] = curr_knot_value;
                double knot_value_in_curve_parameter = Start + (curr_knot_value-knot_interval[0]) * scale_factor;

                tempSpans.push_back(knot_value_in_curve_parameter);
            }
            
        } else if (std::abs(physical_coord_1[1]-physical_coord_2[1]) > tolerance_orientation) {
            // vertical case
            knot_interval[0] = physical_coord_1[1]; 
            knot_interval[1] = physical_coord_2[1];
            std::sort(knot_interval.begin(), knot_interval.end());

            knot_interval[0] -= tolerance_intersection; 
            knot_interval[1] += tolerance_intersection;
            // Compare with surface_spans_v
            for (IndexType i = 0; i < surface_spans_v.size(); i++) {
                double curr_knot_value = surface_spans_v[i];
                if (curr_knot_value < knot_interval[0]) {continue;}
                if (std::abs(curr_knot_value - knot_interval[0]) < tolerance_intersection*10) knot_interval[0] = curr_knot_value;
                if (curr_knot_value > knot_interval[1]) {break;}
                if (std::abs(curr_knot_value - knot_interval[1]) < tolerance_intersection*10) knot_interval[1] = curr_knot_value;
                double knot_value_in_curve_parameter = Start + (curr_knot_value-knot_interval[0]) * scale_factor;
                
                tempSpans.push_back(knot_value_in_curve_parameter);
            }
        } else {
            KRATOS_ERROR << "ComputeAxisIntersectionSBM :: brep not parallel to any surface knot vector";
        }

        rSpans.resize(tempSpans.size());
        if (tempSpans.size() > 2){
            if (reverse_order) {
                for (IndexType i = 0; i < tempSpans.size(); i++) {
                    rSpans[i] = End - tempSpans[tempSpans.size()-1-i];
                }
            } else {
                rSpans = tempSpans;
            }
        } else rSpans = tempSpans;
        
        KRATOS_ERROR_IF(rSpans.size()<2) << "ComputeAxisIntersectionSBM :: Wrong number of intersection found (<2)" << std::endl;
    }

    /* @brief Provides the nurbs boundaries of the NURBS/B-Spline curve.
     * @return domain interval.
     */
    NurbsInterval DomainInterval() const
    {
        return mpNurbsCurve->DomainInterval();
    }

    ///@}
    ///@name Projection Point
    ///@{

    /* @brief Makes projection of rPointGlobalCoordinates to
     *       the closest point on the curve, with
     *       local coordinates rProjectedPointLocalCoordinates.
     *
     * @param Tolerance is the breaking criteria.
     * @return 1 -> projection succeeded
     *         0 -> projection failed
     */
    int ProjectionPointGlobalToLocalSpace(
        const CoordinatesArrayType& rPointGlobalCoordinates,
        CoordinatesArrayType& rProjectedPointLocalCoordinates,
        const double Tolerance = std::numeric_limits<double>::epsilon()
    ) const override
    {
        CoordinatesArrayType point_global_coordinates;

        return ProjectionNurbsGeometryUtilities::NewtonRaphsonCurve(
            rProjectedPointLocalCoordinates,
            rPointGlobalCoordinates,
            point_global_coordinates,
            *this,
            20, Tolerance);
    }

    ///@}
    ///@name IsInside
    ///@{

    /* @brief checks and returns if local coordinate rPointLocalCoordinates[0]
     *        is inside the local/parameter space.
     * @return inside -> 1
     *         outside -> 0
     */
    int IsInsideLocalSpace(
        const CoordinatesArrayType& rPointLocalCoordinates,
        const double Tolerance = std::numeric_limits<double>::epsilon()
    ) const override
    {
        return mpNurbsCurve->IsInsideLocalSpace(rPointLocalCoordinates, Tolerance);
    }

    ///@}
    ///@name ClosestPoint
    ///@{

    /* @brief Makes a check if the provided parameter rPointLocalCoordinates[0]
     *        is inside the curve, or on the boundary or if it lays outside.
     *        If it is outside, it is set to the boundary which is closer to it.
     * @return if rPointLocalCoordinates[0] was before the projection:
     *         on boundary -> 2 - meaning that it is equal to start or end point.
     *         inside -> 1
     *         outside -> 0
     */
    int ClosestPointLocalToLocalSpace(
        const CoordinatesArrayType& rPointLocalCoordinates,
        CoordinatesArrayType& rClosestPointLocalCoordinates,
        const double Tolerance = std::numeric_limits<double>::epsilon()
    ) const override
    {
        return mpNurbsCurve->ClosestPointLocalToLocalSpace(rPointLocalCoordinates, rClosestPointLocalCoordinates, Tolerance);
    }


    ///@}
    ///@name Geometrical Information
    ///@{

    /// Provides the center of the underlying surface
    Point Center() const override
    {
        return mpNurbsSurface->Center();
    }

    /// Computes the length of a nurbs curve
    double Length() const override
    {
        IntegrationPointsArrayType integration_points;
        IntegrationInfo integration_info = GetDefaultIntegrationInfo();
        CreateIntegrationPoints(integration_points, integration_info);

        double length = 0.0;
        for (IndexType i = 0; i < integration_points.size(); ++i) {
            const double determinant_jacobian = DeterminantOfJacobian(integration_points[i]);
            length += integration_points[i].Weight() * determinant_jacobian;
        }
        return length;
    }

    ///@}
    ///@name Jacobian
    ///@{

    double DeterminantOfJacobian(
        const CoordinatesArrayType& rPoint) const override
    {
        std::vector<CoordinatesArrayType> global_space_derivatives(2);
        this->GlobalSpaceDerivatives(
            global_space_derivatives, rPoint, 1);
        return norm_2(global_space_derivatives[1]);
    }

    ///@}
    ///@name Integration Info
    ///@{

    /// Provides the default integration dependent on the polynomial degree of the underlying surface.
    IntegrationInfo GetDefaultIntegrationInfo() const override
    {
        IndexType p = mpNurbsSurface->PolynomialDegreeU() + mpNurbsSurface->PolynomialDegreeV() + 1;
        return IntegrationInfo(1, p, IntegrationInfo::QuadratureMethod::GAUSS);
    }

    ///@}
    ///@name Integration Points
    ///@{

    /* Creates integration points according to the knot intersections
     * of the underlying nurbs surface.
     * @param result integration points.
     */
    void CreateIntegrationPoints(
        IntegrationPointsArrayType& rIntegrationPoints,
        IntegrationInfo& rIntegrationInfo) const override
    {
        std::vector<double> spans;
        SpansLocalSpace(spans);

        IntegrationPointUtilities::CreateIntegrationPoints1D(
            rIntegrationPoints, spans, rIntegrationInfo);
    }

    ///@}
    ///@name Quadrature Point Geometries
    ///@{

    /**
     * @brief This method creates a list of quadrature point geometries
     *        from a list of integration points.
     *
     * @param rResultGeometries list of quadrature point geometries.
     * @param rIntegrationPoints list of integration points.
     * @param NumberOfShapeFunctionDerivatives the number provided
     *        derivatives of shape functions in the system.
     *
     * @see quadrature_point_geometry.h
     */
    void CreateQuadraturePointGeometries(
        GeometriesArrayType& rResultGeometries,
        IndexType NumberOfShapeFunctionDerivatives,
        const IntegrationPointsArrayType& rIntegrationPoints,
        IntegrationInfo& rIntegrationInfo) override
    {
        // shape function container.
        NurbsSurfaceShapeFunction shape_function_container(
            mpNurbsSurface->PolynomialDegreeU(), mpNurbsSurface->PolynomialDegreeV(),
            NumberOfShapeFunctionDerivatives);

        // Resize containers.
        if (rResultGeometries.size() != rIntegrationPoints.size())
            rResultGeometries.resize(rIntegrationPoints.size());

        auto default_method = this->GetDefaultIntegrationMethod();
        SizeType num_nonzero_cps = shape_function_container.NumberOfNonzeroControlPoints();

        Matrix N(1, num_nonzero_cps);
        DenseVector<Matrix> shape_function_derivatives(NumberOfShapeFunctionDerivatives - 1);
        for (IndexType i = 0; i < NumberOfShapeFunctionDerivatives - 1; i++) {
            shape_function_derivatives[i].resize(num_nonzero_cps, i + 2);
        }

        for (IndexType i = 0; i < rIntegrationPoints.size(); ++i)
        {
            std::vector<CoordinatesArrayType> global_space_derivatives(2);
            mpNurbsCurve->GlobalSpaceDerivatives(
                global_space_derivatives,
                rIntegrationPoints[i],
                1);

            if (mpNurbsSurface->IsRational()) {
                shape_function_container.ComputeNurbsShapeFunctionValues(
                    mpNurbsSurface->KnotsU(), mpNurbsSurface->KnotsV(), mpNurbsSurface->Weights(),
                    global_space_derivatives[0][0], global_space_derivatives[0][1]);
            }
            else {
                shape_function_container.ComputeBSplineShapeFunctionValues(
                    mpNurbsSurface->KnotsU(), mpNurbsSurface->KnotsV(),
                    global_space_derivatives[0][0], global_space_derivatives[0][1]);
            }

            /// Get List of Control Points
            PointsArrayType nonzero_control_points(num_nonzero_cps);
            auto cp_indices = shape_function_container.ControlPointIndices(
                mpNurbsSurface->NumberOfControlPointsU(), mpNurbsSurface->NumberOfControlPointsV());
            for (IndexType j = 0; j < num_nonzero_cps; j++) {
                nonzero_control_points(j) = mpNurbsSurface->pGetPoint(cp_indices[j]);
            }
            /// Get Shape Functions N
            for (IndexType j = 0; j < num_nonzero_cps; j++) {
                N(0, j) = shape_function_container(j, 0);
            }

            /// Get Shape Function Derivatives DN_De, ...
            if (NumberOfShapeFunctionDerivatives > 0) {
                IndexType shape_derivative_index = 1;
                for (IndexType n = 0; n < NumberOfShapeFunctionDerivatives - 1; n++) {
                    for (IndexType k = 0; k < n + 2; k++) {
                        for (IndexType j = 0; j < num_nonzero_cps; j++) {
                            shape_function_derivatives[n](j, k) = shape_function_container(j, shape_derivative_index + k);
                        }
                    }
                    shape_derivative_index += n + 2;
                }
            }

            GeometryShapeFunctionContainer<GeometryData::IntegrationMethod> data_container(
                default_method, rIntegrationPoints[i],
                N, shape_function_derivatives);

            rResultGeometries(i) = CreateQuadraturePointsUtility<NodeType>::CreateQuadraturePointCurveOnSurface(
                data_container, nonzero_control_points,
                global_space_derivatives[1][0], global_space_derivatives[1][1], this);
        }
    }


    /**
     * @brief This method creates a list of quadrature point geometries
     *        from a list of integration points in the SBM case.
     *
     * * @details
     * 1. Checks whether the BRep is internal (aligned with specific knot spans) or external.
     *    - Internal boundaries are associated with a fixed span (knot span edge) in the U and V directions.
     *    - External boundaries are associated with body-fitted cases.
     * 2. Iterates over each integration point:
     *    - Computes global space derivatives for the NURBS curve at the integration point.
     *    - Gathers nonzero control points and associated shape function values.
     * 3. Creates a geometry shape function container and assigns it to the result geometries array.
     * 
     * @param rResultGeometries list of quadrature point geometries.
     * @param rIntegrationPoints list of integration points.
     * @param NumberOfShapeFunctionDerivatives the number provided
     *        derivatives of shape functions in the system.
     *
     * @see quadrature_point_geometry.h
     */
    void CreateQuadraturePointGeometriesSBM(
        GeometriesArrayType& rResultGeometries,
        IndexType NumberOfShapeFunctionDerivatives,
        const IntegrationPointsArrayType& rIntegrationPoints,
        IntegrationInfo& rIntegrationInfo) 
    {
        // shape function container.
        NurbsSurfaceShapeFunction shape_function_container(
            mpNurbsSurface->PolynomialDegreeU(), mpNurbsSurface->PolynomialDegreeV(),
            NumberOfShapeFunctionDerivatives);

        // Resize containers.
        if (rResultGeometries.size() != rIntegrationPoints.size())
            rResultGeometries.resize(rIntegrationPoints.size());

        auto default_method = this->GetDefaultIntegrationMethod();
        SizeType num_nonzero_cps = shape_function_container.NumberOfNonzeroControlPoints();

        Matrix N(1, num_nonzero_cps);
        DenseVector<Matrix> shape_function_derivatives(NumberOfShapeFunctionDerivatives - 1);
        for (IndexType i = 0; i < NumberOfShapeFunctionDerivatives - 1; i++) {
            shape_function_derivatives[i].resize(num_nonzero_cps, i + 2);
        }

        // the brep surface in which the brep curve is defined has been marked as sbm in the nurbs_modeler_sbm
        bool is_brep_internal = false;
        IndexType InternalBrepSpanU;
        IndexType InternalBrepSpanV;

        // SBM -> check to which axis the brep is alligned
        std::vector<CoordinatesArrayType> first_integration_point(2); // first integration point of the brep in the parameter space
        std::vector<CoordinatesArrayType> last_integration_point(2); // last integration point of the brep in the parameter space
        mpNurbsCurve->GlobalSpaceDerivatives(first_integration_point,rIntegrationPoints[0],1); 
        mpNurbsCurve->GlobalSpaceDerivatives(last_integration_point,rIntegrationPoints[rIntegrationPoints.size()-1],1); 

        // check if the brep is internal (external breps do not need to be computed in a different knot span) 
        is_brep_internal = true; 
        const double tolerance = 1e-14;
        // At this point all the true are boundary GPs
        if (std::abs(first_integration_point[0][0]-mpNurbsSurface->KnotsU()[0]) < tolerance)
            is_brep_internal = false;
        if (std::abs(first_integration_point[0][0]-mpNurbsSurface->KnotsU()[mpNurbsSurface->KnotsU().size()-1]) < tolerance) 
            is_brep_internal = false;
        if (std::abs(first_integration_point[0][1]-mpNurbsSurface->KnotsV()[0]) < tolerance) 
            is_brep_internal = false;
        if (std::abs(first_integration_point[0][1]-mpNurbsSurface->KnotsV()[mpNurbsSurface->KnotsV().size()-1]) < tolerance) 
            is_brep_internal = false;

        // If is_brep_internal -> the current GP is a surrogate GP otherwise it is an external GP
        if (is_brep_internal) {
            // brep = knot spans edge -> same SpanU and SpanV for all its integration points
            InternalBrepSpanU = NurbsUtilities::GetLowerSpan(mpNurbsSurface->PolynomialDegreeU(), mpNurbsSurface->KnotsU(), first_integration_point[0][0]);
            InternalBrepSpanV = NurbsUtilities::GetLowerSpan(mpNurbsSurface->PolynomialDegreeV(), mpNurbsSurface->KnotsV(), first_integration_point[0][1]);

            // Understand if the knot-boundary is along U or V
            bool gp_is_along_V;
            if (std::abs(first_integration_point[0][0] - last_integration_point[0][0]) < tolerance) {
                gp_is_along_V = true;
            }
            else {gp_is_along_V = false;}

            if (gp_is_along_V && first_integration_point[0][1] > last_integration_point[0][1]) {
                // Is decreasing along y -> Add 1 at InternalBrepSpanU
                InternalBrepSpanU++;
            }
            if (!gp_is_along_V && first_integration_point[0][0] < last_integration_point[0][0]) {
                // Is increasing along x -> Add 1 at InternalBrepSpanV
                InternalBrepSpanV++;
            }
        }

        // Loop over the integration points of the brep
        for (IndexType i = 0; i < rIntegrationPoints.size(); ++i) 
        {
            std::vector<CoordinatesArrayType> global_space_derivatives(2);
            mpNurbsCurve->GlobalSpaceDerivatives(
                global_space_derivatives,
                rIntegrationPoints[i],
                1);

            if (mpNurbsSurface->IsRational()) { 
                shape_function_container.ComputeNurbsShapeFunctionValues(
                    mpNurbsSurface->KnotsU(), mpNurbsSurface->KnotsV(), mpNurbsSurface->Weights(),
                    global_space_derivatives[0][0], global_space_derivatives[0][1]);
            }
            else {
                if (is_brep_internal) {
                    // SBM case on internal brep 
                    shape_function_container.ComputeBSplineShapeFunctionValuesAtSpan(
                        mpNurbsSurface->KnotsU(),
                        mpNurbsSurface->KnotsV(),
                        InternalBrepSpanU,
                        InternalBrepSpanV,
                        global_space_derivatives[0][0],
                        global_space_derivatives[0][1]);
                }
                else {
                    // SBM case on external brep
                    shape_function_container.ComputeBSplineShapeFunctionValues(
                        mpNurbsSurface->KnotsU(), mpNurbsSurface->KnotsV(),
                        global_space_derivatives[0][0], global_space_derivatives[0][1]);
                } 
            }

            /// Get List of Control Points
            PointsArrayType nonzero_control_points(num_nonzero_cps);
            auto cp_indices = shape_function_container.ControlPointIndices(
                mpNurbsSurface->NumberOfControlPointsU(), mpNurbsSurface->NumberOfControlPointsV());
            for (IndexType j = 0; j < num_nonzero_cps; j++) {
                nonzero_control_points(j) = mpNurbsSurface->pGetPoint(cp_indices[j]);
            }
            /// Get Shape Functions N
            for (IndexType j = 0; j < num_nonzero_cps; j++) {
                N(0, j) = shape_function_container(j, 0);
            }

            /// Get Shape Function Derivatives DN_De, ...
            if (NumberOfShapeFunctionDerivatives > 0) {
                IndexType shape_derivative_index = 1;
                for (IndexType n = 0; n < NumberOfShapeFunctionDerivatives - 1; n++) {
                    for (IndexType k = 0; k < n + 2; k++) {
                        for (IndexType j = 0; j < num_nonzero_cps; j++) {
                            shape_function_derivatives[n](j, k) = shape_function_container(j, shape_derivative_index + k);
                        }
                    }
                    shape_derivative_index += n + 2;
                }
            }

            GeometryShapeFunctionContainer<GeometryData::IntegrationMethod> data_container(
                default_method, rIntegrationPoints[i],
                N, shape_function_derivatives);

            rResultGeometries(i) = CreateQuadraturePointsUtility<NodeType>::CreateQuadraturePointCurveOnSurface(
                data_container, nonzero_control_points,
                global_space_derivatives[1][0], global_space_derivatives[1][1], this);
        }
    }

    ///@}
    ///@name Operation within Global Space
    ///@{

    /*
    * @brief This method maps from dimension space to working space.
    * From Piegl and Tiller, The NURBS Book, Algorithm A3.1/ A4.1
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
        // Compute the coordinates of the embedded curve in the parametric space of the surface
        CoordinatesArrayType result_local = mpNurbsCurve->GlobalCoordinates(rResult, rLocalCoordinates);

        // Compute and return the coordinates of the surface in the geometric space
        return mpNurbsSurface->GlobalCoordinates(rResult, result_local);
    }

    /**
    * @brief This method maps from dimension space to working space and computes the
    *        number of derivatives at the dimension parameter.
    * From ANurbs library (https://github.com/oberbichler/ANurbs)
    * @param LocalCoordinates The local coordinates in dimension space
    * @param Derivative Number of computed derivatives
    * @return std::vector<array_1d<double, 3>> with the coordinates in working space
    * @see PointLocalCoordinates
    */
    void GlobalSpaceDerivatives(
        std::vector<CoordinatesArrayType>& rGlobalSpaceDerivatives,
        const CoordinatesArrayType& rCoordinates,
        const SizeType DerivativeOrder) const override
    {
        // Check size of output
        if (rGlobalSpaceDerivatives.size() != DerivativeOrder + 1) {
            rGlobalSpaceDerivatives.resize(DerivativeOrder + 1);
        }

        // Compute the gradients of the embedded curve in the parametric space of the surface
        std::vector<array_1d<double, 3>> curve_derivatives;
        mpNurbsCurve->GlobalSpaceDerivatives(curve_derivatives, rCoordinates, DerivativeOrder);

        // Compute the gradients of the surface in the geometric space
        array_1d<double, 3> surface_coordinates =  ZeroVector(3);
        surface_coordinates[0] = curve_derivatives[0][0];
        surface_coordinates[1] = curve_derivatives[0][1];
        std::vector<array_1d<double, 3>> surface_derivatives;
        mpNurbsSurface->GlobalSpaceDerivatives(surface_derivatives, surface_coordinates, DerivativeOrder);

        std::function<array_1d<double, 3>(int, int, int)> c;
        c = [&](int DerivativeOrder, int i, int j) -> array_1d<double, 3> {
            if (DerivativeOrder > 0) {
                array_1d<double, 3> result = ZeroVector(3);

                for (int a = 1; a <= DerivativeOrder; a++) {
                    result += (
                        c(DerivativeOrder - a, i + 1, j) * curve_derivatives[a][0] +
                        c(DerivativeOrder - a, i, j + 1) * curve_derivatives[a][1]
                        ) * NurbsUtilities::GetBinomCoefficient(DerivativeOrder - 1, a - 1);
                }

                return result;
            }
            else {
                const int index = NurbsSurfaceShapeFunction::IndexOfShapeFunctionRow(i, j);
                return surface_derivatives[index];
            }
        };
        for (SizeType i = 0; i <= DerivativeOrder; i++) {
            rGlobalSpaceDerivatives[i] = c(i, 0, 0);
        }
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
        return GeometryData::KratosGeometryFamily::Kratos_Nurbs;
    }

    /**
     * @brief Gets the geometry type.
     * @details This function returns the specific type of the geometry. The geometry type provides a more detailed classification of the geometry.
     * @return GeometryData::KratosGeometryType The specific geometry type.
     */
    GeometryData::KratosGeometryType GetGeometryType() const override
    {
        return GeometryData::KratosGeometryType::Kratos_Nurbs_Curve_On_Surface;
    }

    ///@}
    ///@name Information
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "2 dimensional nurbs curve on 3D surface.";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "2 dimensional nurbs curve on 3D surface.";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }

    ///@}
private:
    ///@name Private Static Member Variables
    ///@{

    static const GeometryData msGeometryData;

    static const GeometryDimension msGeometryDimension;

    ///@}
    ///@name Private Member Variables
    ///@{

    typename NurbsSurfaceType::Pointer mpNurbsSurface;
    typename NurbsCurveType::Pointer mpNurbsCurve;

    ///@}
    ///@name Private Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseType);
        rSerializer.save("pNurbsSurface", mpNurbsSurface);
        rSerializer.save("pNurbsCurve", mpNurbsCurve);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseType);
        rSerializer.load("pNurbsSurface", mpNurbsSurface);
        rSerializer.load("pNurbsCurve", mpNurbsCurve);
    }

    ///@}

}; // class NurbsCurveOnSurfaceGeometry

template<int TWorkingSpaceDimension, class TCurveContainerPointType, class TSurfaceContainerPointType>
const GeometryData NurbsCurveOnSurfaceGeometry<TWorkingSpaceDimension, TCurveContainerPointType, TSurfaceContainerPointType>::msGeometryData(
    &msGeometryDimension,
    GeometryData::IntegrationMethod::GI_GAUSS_1,
    {}, {}, {});

template<int TWorkingSpaceDimension, class TCurveContainerPointType, class TSurfaceContainerPointType>
const GeometryDimension NurbsCurveOnSurfaceGeometry<TWorkingSpaceDimension, TCurveContainerPointType, TSurfaceContainerPointType>::msGeometryDimension(TWorkingSpaceDimension, 1);

} // namespace Kratos

#endif // KRATOS_NURBS_CURVE_ON_SURFACE_H_INCLUDED defined

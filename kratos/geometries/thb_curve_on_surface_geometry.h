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

#if !defined(KRATOS_THB_NURBS_CURVE_ON_SURFACE_H_INCLUDED )
#define  KRATOS_THB_NURBS_CURVE_ON_SURFACE_H_INCLUDED

// External includes
#include <gismo.h>
using namespace gismo;

// Project includes
#include "geometries/geometry.h"

#include "geometries/nurbs_curve_geometry.h"
#include "geometries/thb_surface_geometry.h"
#include "geometries/nurbs_shape_function_utilities/nurbs_curve_shape_functions.h"
#include "geometries/nurbs_shape_function_utilities/nurbs_interval.h"

#include "utilities/curve_axis_intersection.h"

#include "utilities/nurbs_utilities/projection_nurbs_geometry_utilities.h"

namespace Kratos {

template <int TWorkingSpaceDimension, class TCurveContainerPointType, class TSurfaceContainerPointType>
class THBCurveOnSurfaceGeometry : public Geometry<typename TSurfaceContainerPointType::value_type>
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

    typedef THBSurfaceGeometry<3, TSurfaceContainerPointType> THBSurfaceType;
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

    /// Counted pointer of THBCurveOnSurfaceGeometry
    KRATOS_CLASS_POINTER_DEFINITION(THBCurveOnSurfaceGeometry);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    THBCurveOnSurfaceGeometry(
        typename THBSurfaceType::Pointer pSurface,
        typename NurbsCurveType::Pointer pCurve)
        : BaseType(PointsArrayType(), &msGeometryData)
        , mpNurbsSurface(pSurface)
        , mpNurbsCurve(pCurve)
    {
    }

    /// Default constructor
    THBCurveOnSurfaceGeometry()
        : BaseType(PointsArrayType(), &msGeometryData)
    {};

    /// Copy constructor
    THBCurveOnSurfaceGeometry(THBCurveOnSurfaceGeometry const& rOther)
        : BaseType(rOther)
        , mpNurbsSurface(rOther.mpNurbsSurface)
        , mpNurbsCurve(rOther.mpNurbsCurve)
    {
    }

    /// Copy constructor, with different point type.
    template<class TOtherCurveContainerPointType, class TOtherSurfaceContainerPointType> THBCurveOnSurfaceGeometry(
        THBCurveOnSurfaceGeometry<TWorkingSpaceDimension, TOtherCurveContainerPointType, TOtherSurfaceContainerPointType> const& rOther)
        : BaseType(rOther, &msGeometryData)
        , mpNurbsSurface(rOther.mpNurbsSurface)
        , mpNurbsCurve(rOther.mpNurbsCurve)
    {
    }

    /// Destructor
    ~THBCurveOnSurfaceGeometry() override = default;

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator
    THBCurveOnSurfaceGeometry& operator=(const THBCurveOnSurfaceGeometry& rOther)
    {
        BaseType::operator=(rOther);
        mpNurbsSurface = rOther.mpNurbsSurface;
        mpNurbsCurve = rOther.mpNurbsCurve;
        return *this;
    }

    /// @brief Assignment operator for geometries with different point type.
    template<class TOtherCurveContainerPointType, class TOtherSurfaceContainerPointType>
    THBCurveOnSurfaceGeometry& operator=(
        THBCurveOnSurfaceGeometry<TWorkingSpaceDimension, TOtherCurveContainerPointType, TOtherSurfaceContainerPointType> const & rOther)
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
        KRATOS_ERROR << "THBCurveOnSurfaceGeometry cannot be created with 'PointsArrayType const& ThisPoints'. "
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
    ///@name Mathematical Informations
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
    void SpansLocalSpace(std::vector<double>& rSpans, IndexType DirectionIndex = 0) const override
    {
        auto interval = mpNurbsCurve->DomainInterval();
        this->SpansLocalSpace(rSpans, interval.GetT0(), interval.GetT1());
    }

    /* @brief  Provides intersections of the nurbs curve with the knots of the surface.
     * @return vector of interval limitations.
     */
    void SpansLocalSpace(std::vector<double>& rSpans,
        double Start, double End) const
    {
        std::vector<double> surface_spans_u;
        std::vector<double> surface_spans_v;
        mpNurbsSurface->SpansLocalSpace(surface_spans_u, 0);
        mpNurbsSurface->SpansLocalSpace(surface_spans_v, 1);

        CurveAxisIntersection<CurveNodeType>::ComputeAxisIntersection(
            rSpans,
            *(mpNurbsCurve.get()), Start, End,
            surface_spans_u, surface_spans_v,
            1e-6);
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

    /* @brief Makes a check if the provided paramater rPointLocalCoordinates[0]
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
    ///@name Geometrical Informations
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

        gismo::gsMatrix<> ThbResult;
        gismo::gsMatrix<> ThbResultDeriv;
        gismo::gsMatrix<> ThbResultDeriv2;
        gismo::gsMatrix<int> ThbIndices;
        gismo::gsMatrix<int> ThbIndicesDeriv;
        gismo::gsMatrix<int> ThbIndicesDeriv2;
        gismo::gsMatrix<> rLocalCoordinateThb(2,1);

        for (IndexType i = 0; i < rIntegrationPoints.size(); ++i)
        {
            std::vector<CoordinatesArrayType> global_space_derivatives(2);
            mpNurbsCurve->GlobalSpaceDerivatives(
                global_space_derivatives,
                rIntegrationPoints[i],
                1);

            //G+Smo
            gismo::gsTHBSpline<2> thb_geom = mpNurbsSurface->ToGismo();
            gismo::gsTHBSplineBasis<2> thb = thb_geom.basis();

            rLocalCoordinateThb(0,0) = std::floor(global_space_derivatives[0][0] * 1e12) / 1e12;
            rLocalCoordinateThb(1,0) = std::floor(global_space_derivatives[0][1] * 1e12) / 1e12;
            
            thb.active_into(rLocalCoordinateThb,ThbIndices);
            thb.eval_into(rLocalCoordinateThb,ThbResult);

            thb.active_into(rLocalCoordinateThb,ThbIndicesDeriv);
            thb.deriv_into(rLocalCoordinateThb,ThbResultDeriv);

            thb.active_into(rLocalCoordinateThb,ThbIndicesDeriv2);
            thb.deriv2_into(rLocalCoordinateThb,ThbResultDeriv2);

            SizeType num_nonzero_cps = ThbIndices.size();

            Matrix N(1, num_nonzero_cps);
            DenseVector<Matrix> shape_function_derivatives(NumberOfShapeFunctionDerivatives - 1);
            for (IndexType i = 0; i < NumberOfShapeFunctionDerivatives - 1; i++) {
                shape_function_derivatives[i].resize(num_nonzero_cps, i + 2);
            }

            /// Get List of Control Points
            PointsArrayType nonzero_control_points(num_nonzero_cps);
            for (IndexType j = 0; j < num_nonzero_cps; j++) {
                nonzero_control_points(j) =  mpNurbsSurface->pGetPoint(ThbIndices(j,0));
            }
            /// Get Shape Functions N
            for (IndexType j = 0; j < num_nonzero_cps; j++) {
                N(0, j) = ThbResult(j, 0);
            }
            
            /// Get Shape Function Derivatives DN_De, ...
            if (NumberOfShapeFunctionDerivatives > 0) {
                for (IndexType n = 0; n < NumberOfShapeFunctionDerivatives - 1; n++) {
                    if(n == 0){
                    for (IndexType k = 0; k < n + 2; k++) {
                        for (IndexType j = 0; j < num_nonzero_cps; j++) {
                            shape_function_derivatives[n](j, k) = ThbResultDeriv(j*2 + k , 0);
                        }
                    }
                    }
                    else if(n == 1){
                        for (IndexType j = 0; j < num_nonzero_cps; j++) {
                            shape_function_derivatives[n](j, 0) = ThbResultDeriv2(j*3, 0);
                        }
                        for (IndexType j = 0; j < num_nonzero_cps; j++) {
                            shape_function_derivatives[n](j, 1) = ThbResultDeriv2(j*3 + 2 , 0);
                        }
                        for (IndexType j = 0; j < num_nonzero_cps; j++) {
                            shape_function_derivatives[n](j, 2) = ThbResultDeriv2(j*3 + 1 , 0);
                        }
                    }
                    else
                    {
                        KRATOS_THROW_ERROR( std::invalid_argument, "third derivatives and more are not yet defined!", "" )
                    }
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

    GeometryData::KratosGeometryFamily GetGeometryFamily() const override
    {
        return GeometryData::KratosGeometryFamily::Kratos_Nurbs;
    }

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

    typename THBSurfaceType::Pointer mpNurbsSurface;
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

}; // class THBCurveOnSurfaceGeometry

template<int TWorkingSpaceDimension, class TCurveContainerPointType, class TSurfaceContainerPointType>
const GeometryData THBCurveOnSurfaceGeometry<TWorkingSpaceDimension, TCurveContainerPointType, TSurfaceContainerPointType>::msGeometryData(
    &msGeometryDimension,
    GeometryData::IntegrationMethod::GI_GAUSS_1,
    {}, {}, {});

template<int TWorkingSpaceDimension, class TCurveContainerPointType, class TSurfaceContainerPointType>
const GeometryDimension THBCurveOnSurfaceGeometry<TWorkingSpaceDimension, TCurveContainerPointType, TSurfaceContainerPointType>::msGeometryDimension(TWorkingSpaceDimension, 1);

} // namespace Kratos

#endif // KRATOS_NURBS_CURVE_ON_SURFACE_H_INCLUDED defined

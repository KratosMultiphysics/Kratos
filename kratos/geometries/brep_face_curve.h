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
//  contributors:    Pooyan Dadvand
//                   Philipp Bucher
//

#if !defined(KRATOS_BREP_FACE_CURVE_3D_H_INCLUDED )
#define  KRATOS_BREP_FACE_CURVE_3D_H_INCLUDED

// System includes

// External includes

// Project includes
#include "geometries/geometry.h"
#include "geometries/geometry_dimension.h"


namespace Kratos
{
///@name Kratos Classes
///@{

/**
 * @class BrepFaceCurve
 * @ingroup KratosCore
 * @brief The BrepFaceCurve acts as topology for faces having a pointer
 *        to a surface geometry and a list of BrepEdge s. Those edges
 *        define the trimming of the respective surface.
 */
template<class TPointType> class BrepFaceCurve
    : public Geometry<TPointType>
{
public:
    ///@}
    ///@name Type Definitions
    ///@{

    /** Pointer definition of BrepFaceCurve */
    KRATOS_CLASS_POINTER_DEFINITION( BrepFaceCurve );

    typedef TPointType PointType;

    typedef Geometry<TPointType> BaseType;
    typedef Geometry<TPointType> GeometryType;

    typedef GeometryData::IntegrationMethod IntegrationMethod;

    typedef NurbsSurfaceGeometry<3, TPointType> NurbsSurfaceType;
    typedef NurbsCurveGeometry<2, TPointType> NurbsCurveType;

    typedef typename BaseType::GeometriesArrayType GeometriesArrayType;
    typedef typename BaseType::IndexType IndexType;
    typedef typename BaseType::SizeType SizeType;
    typedef typename BaseType::PointsArrayType PointsArrayType;
    typedef typename BaseType::CoordinatesArrayType CoordinatesArrayType;

    typedef typename BaseType::IntegrationPointType IntegrationPointType;
    typedef typename BaseType::IntegrationPointsArrayType IntegrationPointsArrayType;
    typedef typename BaseType::IntegrationPointsContainerType IntegrationPointsContainerType;

    typedef typename BaseType::ShapeFunctionsValuesContainerType ShapeFunctionsValuesContainerType;
    typedef typename BaseType::ShapeFunctionsLocalGradientsContainerType ShapeFunctionsLocalGradientsContainerType;
    typedef typename BaseType::ShapeFunctionsGradientsType ShapeFunctionsGradientsType;
    typedef typename BaseType::ShapeFunctionsSecondDerivativesType ShapeFunctionsSecondDerivativesType;
    typedef typename BaseType::ShapeFunctionsThirdDerivativesType ShapeFunctionsThirdDerivativesType;

    typedef typename BaseType::JacobiansType JacobiansType;
    typedef typename BaseType::NormalType NormalType;

    // New Shape function typedefs
    //typedef GeometryData::ShapeFunctionsDerivativesType ShapeFunctionsDerivativesType;
    //typedef GeometryData::ShapeFunctionsDerivativesIntegrationPointsType ShapeFunctionsDerivativesIntegrationPointsType;
    //typedef GeometryData::ShapeFunctionsDerivativesVectorType ShapeFunctionsDerivativesVectorType;


    ///@}
    ///@name Life Cycle
    ///@{

    BrepFaceCurve( 
        typename NurbsSurfaceType::Pointer pSurface,
        typename NurbsCurveType::Pointer pCurve)
        : BaseType(PointsArrayType(), &mGeometryData)
        , mpNurbsSurface(pSurface)
        , mpNurbsCurve(pCurve)
        , mGeometryData(
            &msGeometryDimension,
            GeometryData::GI_GAUSS_1,
            {}, {}, {})
    {
    }

    explicit BrepFaceCurve(const PointsArrayType& ThisPoints)
        : BaseType(ThisPoints, &mGeometryData)
        , mGeometryData(
            &msGeometryDimension,
            GeometryData::GI_GAUSS_1,
            {}, {}, {})
    {
    }

    /**
     * Copy constructor.
     * Construct this geometry as a copy of given geometry.
     *
     * @note This copy constructor does not copy the points and new
     * geometry shares points with given source geometry. It is
     * obvious that any change to this new geometry's point affect
     * source geometry's points too.
     */
    BrepFaceCurve( BrepFaceCurve const& rOther )
        : BaseType( rOther )
    {
    }

    /**
     * Copy constructor from a geometry with other point type.
     * Construct this geometry as a copy of given geometry which
     * has different type of points. The given goemetry's
     * TOtherPointType* must be implicity convertible to this
     * geometry PointType.
     *
     * @note This copy constructor does not copy the points and new
     * geometry shares points with given source geometry. It is
     * obvious that any change to this new geometry's point affect
     * source geometry's points too.
     */
    template<class TOtherPointType> explicit BrepFaceCurve( BrepFaceCurve<TOtherPointType> const& rOther )
        : BaseType( rOther )
    {
    }

    /**
     * Destructor. Does nothing!!!
     */
    ~BrepFaceCurve() override {}

    GeometryData::KratosGeometryFamily GetGeometryFamily() const override
    {
        return GeometryData::Kratos_generic_family;
    }

    GeometryData::KratosGeometryType GetGeometryType() const override
    {
        return GeometryData::Kratos_generic_type;
    }

    ///@}
    ///@name Operators
    ///@{

    /**
     * Assignment operator.
     *
     * @note This operator don't copy the points and this
     * geometry shares points with given source geometry. It's
     * obvious that any change to this geometry's point affect
     * source geometry's points too.
     *
     * @see Clone
     * @see ClonePoints
     */
    BrepFaceCurve& operator=( const BrepFaceCurve& rOther )
    {
        BaseType::operator=( rOther );
        return *this;
    }

    /**
     * Assignment operator for geometries with different point type.
     *
     * @note This operator don't copy the points and this
     * geometry shares points with given source geometry. It's
     * obvious that any change to this geometry's point affect
     * source geometry's points too.
     *
     * @see Clone
     * @see ClonePoints
     */
    template<class TOtherPointType>
    BrepFaceCurve& operator=( BrepFaceCurve<TOtherPointType> const & rOther )
    {
        BaseType::operator=( rOther );
        return *this;
    }

    ///@}
    ///@name Operations
    ///@{

    typename BaseType::Pointer Create( PointsArrayType const& ThisPoints ) const override
    {
        return typename BaseType::Pointer( new BrepFaceCurve( ThisPoints ) );
    }

    ///@}

    /**
    * @brief Calculates the boundingbox of the geometry.
    * @details Corresponds with the highest and lowest point in space
    * @param rLowPoint  Lower point of the boundingbox.
    * @param rHighPoint Higher point of the boundingbox.
    */
    //bool IsInsideBoundingBox(
    //    const TPointType& rPoint,
    //    double BoundaryTolerance = std::numeric_limits<double>::epsilon) const
    //{
    //    const SizeType working_space_dimension = WorkingSpaceDimension();

    //    TPointType LowPoint;
    //    TPointType HighPoint;

    //    BoundingBox(LowPoint, HighPoint);

    //    for (IndexType i = 0; i < working_space_dimension; ++i)
    //    {
    //        if ((LowPoint[i] - BoundaryTolerance > rPoint[i])
    //            || (HighPoint[i] + BoundaryTolerance < rPoint[i])
    //            return false;
    //    }
    //    return true;
    //}

    /**
    * @brief Calculates the boundingbox of the geometry.
    * @details Corresponds with the highest and lowest point in space
    * @param rLowPoint  Lower point of the boundingbox.
    * @param rHighPoint Higher point of the boundingbox.
    */
    //override void BoundingBox(
    //    TPointType& rLowPoint,
    //    TPointType& rHighPoint
    //) const
    //{
    //    const SizeType working_space_dimension = WorkingSpaceDimension();

    //    for (auto point in mLinearizedPolygon) { //The first node is already assigned, so we can start from 1
    //        const auto& r_point = this->GetPoint(point);
    //        for (IndexType i = 0; i < working_space_dimension; ++i) {
    //            rHighPoint[i] = (rHighPoint[i] < r_point[i]) ? r_point[i] : rHighPoint[i];
    //            rLowPoint[i] = (rLowPoint[i]  > r_point[i]) ? r_point[i] : rLowPoint[i];
    //        }
    //    }
    //}

    //bool GetClosestPoint()
    //{

    //}

    //std::vector<TPointType> GetPointsInRadius(
    //    const TPointType& rPoint,
    //    const double Radius,
    //    const int MaxIterations,
    //    const double Accuracy
    //    )
    //{
    //    std::vector<double> parameters;
    //    if (IsInsideBoundingBox(rPoint, Radius))
    //    {
    //        SearchRadius = Radius * 2;
    //        for (auto point in mLinearizedPolygon)
    //        {
    //            if (norm_2(rPoint - point) < SearchRadius)
    //            {
    //                double new_parameter = 0.0;
    //                if (NewtonRaphson<3>(
    //                    new_parameter,
    //                    point.t,
    //                    mp_nurbs_curve,
    //                    MaxIterations,
    //                    Radius,
    //                    Accuracy))
    //                {
    //                    for (int t = 0; t < parameters.size(); ++t)
    //                    {
    //                        if (parameters[t] - new_parameter < Accuracy)
    //                        {
    //                            parameters.push_back(new_parameter);
    //                            break;
    //                        }
    //                    }
    //                }
    //            }
    //        }
    //    }
    //    std::vector<TPointType> new_points;
    //    for (int t = 0; t < parameters.size(); ++t)
    //    {
    //        new_points.push_back(mp_curve->PointAt(parameters[t]));
    //    }
    //    return new_points;
    //}

    //GeometryType::Pointer GetIntegrationPointSurface(const double CurveParameter)
    //{

    //}

    /**
    * @brief Returns the local coordinates of a given arbitrary point
    * @param rResult The vector containing the local coordinates of the point
    * @param rPoint The point in global coordinates
    * @return The vector containing the local coordinates of the point
    */
    virtual CoordinatesArrayType& PointLocalCoordinates(
        CoordinatesArrayType& rResult,
        const CoordinatesArrayType& rPoint
    ) const
    {
        // ONLY PART OF THE UTILITY

        // This function would be better with a true false return

        // Get initial guess from any point on surface - see the point list.

        //external utility
        iga_surface_utilities::NewtonRaphson(
            const std::shared_ptr<NodeSurfaceGeometry3D>& pSurface,
            const array_1d<double, 3>& rPoint,
            const double InitialU,
            const double InitialV,
            double& rU,
            double& rV,
            const double Accuracy,
            const double Tolerance,
            const double MaxIterations
        );
    }

    /** This method maps from dimension space to working space.
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
        pCurve->GlobalCoordinates(rResult, rLocalCoordinates);

        pSurface->GlobalCoordinates(rResult, rLocalCoordinates);

        return rResult;
    }

    ///@name Shape Function
    ///@{

    //override ShapeFunctionsDerivativesType& ShapeFunctionsDerivative(
    //    IndexType DerivativeOrder,
    //    Vector &rResult,
    //    const CoordinatesArrayType& rCoordinates) const
    //{
    //    return pSurface->ShapeFunctionsDerivative(DerivativeOrder, rCoordinates);
    //}

    //override ShapeFunctionsDerivativesVectorType& ShapeFunctionsDerivatives(
    //    IndexType DerivativeOrder,
    //    Vector &rResult,
    //    const CoordinatesArrayType& rCoordinates) const
    //{
    //    return pSurface->ShapeFunctionsDerivatives(DerivativeOrder, rCoordinates);
    //}

    ///@}
    ///@name Input and output
    ///@{

    //std::vector<geometry::Pointer> GetIntegrationPointGeomtries()
    //{
    //    if (!mIsGeometryDataInitialized)
    //        ComputeGeometryData();

    //    SizeType number_of_integration_points = mpGeometryData->IntegrationPointsNumber();

    //    IntegrationPointsArrayType integration_points = mpGeometryData->IntegrationPoints();
    //    ShapeFunctionsDerivativesVectorType shape_functions_derivatives_all = mpGeometryData->ShapeFunctionsDerivativesAll();

    //    SizeType number_of_derivatives = shape_functions_derivatives_all.size();

    //    std::vector<geometry> IntegrationPointList(number_of_integration_points);
    //    for (i = 0.0; i < number_of_integration_points; ++i)
    //    {
    //        ShapeFunctionsDerivativesVectorType shape_functions_derivatives_all_one_integration_point(number_of_derivatives);

    //        for (j = 0.0; j < number_of_integration_points; ++j)
    //        {
    //            shape_functions_derivatives_all_one_integration_point[j] = shape_functions_derivatives_all[j][i];
    //        }

    //        IntegrationPointSurface3d(integration_points[i], shape_functions_derivatives_all_one_integration_point);
    //    }
    //}

    //std::vector<geometry::Pointer> GetLineGeometries()
    //{
    //    if (!mIsPolygonInitialized)
    //        ComputeGeometryData();

    //    SizeType number_of_integration_points = mpGeometryData->IntegrationPointsNumber();

    //    IntegrationPointsArrayType integration_points = mpGeometryData->IntegrationPoints();
    //    ShapeFunctionsDerivativesVectorType shape_functions_derivatives_all = mpGeometryData->ShapeFunctionsDerivativesAll();

    //    SizeType number_of_derivatives = shape_functions_derivatives_all.size();

    //    std::vector<geometry> IntegrationPointList(number_of_integration_points);
    //    for (i = 0.0; i < number_of_integration_points; ++i)
    //    {
    //        ShapeFunctionsDerivativesVectorType shape_functions_derivatives_all_one_integration_point(number_of_derivatives);

    //        for (j = 0.0; j < number_of_integration_points; ++j)
    //        {
    //            shape_functions_derivatives_all_one_integration_point[j] = shape_functions_derivatives_all[j][i];
    //        }

    //        IntegrationPointSurface3d(integration_points[i], shape_functions_derivatives_all_one_integration_point);
    //    }
    //}


    //LinearizedPolygon GetLinearizedPolygon()
    //{
    //    if (!mIsLinearizedPolygonInitialized)
    //    {
    //        this->InitializeLinearizedPolygon();
    //    }
    //}


    /**
     * Turn back information as a string.
     *
     * @return String contains information about this geometry.
     * @see PrintData()
     * @see PrintInfo()
     */
    std::string Info() const override
    {
        return "Brep face curve";
    }

    /**
     * Print information about this object.
     * @param rOStream Stream to print into it.
     * @see PrintData()
     * @see Info()
     */
    void PrintInfo( std::ostream& rOStream ) const override
    {
        rOStream << "Brep face curve";
    }

    /**
     * Print geometry's data into given stream.
     * Prints it's points
     * by the order they stored in the geometry and then center
     * point of geometry.
     *
     * @param rOStream Stream to print into it.
     * @see PrintInfo()
     * @see Info()
     */
    void PrintData( std::ostream& rOStream ) const override
    {
        BaseType::PrintData( rOStream );
        std::cout << std::endl;
        rOStream << "    Brep face curve : " << std::endl;
    }

    ///@}

protected:

private:
    ///@name Static Member Variables
    ///@{

    static const GeometryDimension msGeometryDimension;

    ///@}
    ///@name Member Variables
    ///@{

    GeometryData mGeometryData;

    typename NurbsSurfaceType::Pointer mpNurbsSurface;
    typename NurbsCurveType::Pointer mpNurbsCurve;

    //bool m_is_geometry_data_initialized;
    //GeometryData mGeometryData;

    //bool m_is_linearized_polygon_initialized;
    //LinearizedPolygon mLinearizedPolygon;

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save( Serializer& rSerializer ) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, BaseType );
        rSerializer.save("NurbsSurface", mpNurbsSurface);
        rSerializer.save("NurbsCurve", mpNurbsCurve);
        rSerializer.save("GeometryData", mGeometryData);
    }

    void load( Serializer& rSerializer ) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BaseType );
        rSerializer.load("NurbsSurface", mpNurbsSurface);
        rSerializer.load("NurbsCurve", mpNurbsCurve);
        rSerializer.load("GeometryData", mGeometryData);
    }

    BrepFaceCurve()
        : BaseType( PointsArrayType(), &mGeometryData )
        , mGeometryData(
            &msGeometryDimension,
            GeometryData::GI_GAUSS_1,
            {}, {}, {})
    {}

    ///@}
    ///@name Private Operations
    ///@{

    //void ComputeGeometryData(IntegrationMethod)
    //{
    //    if (mIsBoundaryPolygonInitialized)
    //        CreateBoundaryPolygon();

    //    //Make Clipper or equal

    //    //Create Integration Points

    //    int degree_u = pSurface->DegreeU();
    //    int degree_v = pSurface->DegreeV();

    //    int integration_degree = std::max(degree_u, degree_v) + 1;

    //    for (int i = 0; i < rClipper.NbSpansU(); ++i)
    //    {
    //        for (int j = 0; j < rClipper.NbSpansU(); ++j)
    //        {
    //            if (rClipper.SpanTrimType(i, j) == ANurbs::Empty)
    //            {
    //                continue;
    //            }
    //            else if (rClipper.SpanTrimType(i, j) == ANurbs::Full)
    //            {
    //                auto integration_points = ANurbs::IntegrationPoints<double>::Points2(
    //                    degree_u + 1,
    //                    degree_v + 1,
    //                    rClipper.SpanU(i),
    //                    rClipper.SpanV(j));

    //                for (int i = 0; i < integration_points.size(); ++i)
    //                {
    //                    mpNurbsSurface->Compute(
    //                        integration_point_polygon.IntegrationPoint(i).u,
    //                        integration_point_polygon.IntegrationPoint(i).v);
    //                }
    //            }
    //            else if (rClipper.SpanTrimType(i, j) == ANurbs::Trimmed)
    //            {
    //                auto polygons = rClipper.SpanPolygons(i, j);

    //                for (int p = 0; p < polygons.size(); ++p)
    //                {
    //                    auto integration_point_polygon = ANurbs::PolygonIntegrationPoints<Kratos::array_1d<double, 2>>();

    //                    integration_point_polygon.Compute(degree, polygons[p]);

    //                    for (int i = 0; i < integration_point_polygon.NbIntegrationPoints(); ++i)
    //                    {
    //                        mpNurbsSurface->Compute(
    //                            integration_point_polygon.IntegrationPoint(i).u,
    //                            integration_point_polygon.IntegrationPoint(i).v);

    //                        int number_of_non_zero_cps = shape.NonzeroPoleIndices().size();
    //                    }
    //                }
    //            }
    //        }
    //    }
    //    mGeometryData.SetData(...);
    //}


    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Private Friends
    ///@{

    template<class TOtherPointType> friend class BrepFaceCurve;

    ///@}
    ///@name Un accessible methods
    ///@{



    ///@}
}; // Class Geometry

///@}
///@name Type Definitions
///@{


///@}
///@name Input and output
///@{
/**
 * input stream functions
 */
template<class TPointType> inline std::istream& operator >> (
    std::istream& rIStream,
    BrepFaceCurve<TPointType>& rThis );
/**
 * output stream functions
 */
template<class TPointType> inline std::ostream& operator << (
    std::ostream& rOStream,
    const BrepFaceCurve<TPointType>& rThis )
{
    rThis.PrintInfo( rOStream );
    rOStream << std::endl;
    rThis.PrintData( rOStream );
    return rOStream;
}

template<class TPointType>
const GeometryDimension BrepFaceCurve<TPointType>::msGeometryDimension(
    1,
    3,
    2);

///@}
}// namespace Kratos.

#endif // KRATOS_BREP_FACE_CURVE_3D_H_INCLUDED  defined

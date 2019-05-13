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

#if !defined(KRATOS_BREP_CURVE_3D_H_INCLUDED )
#define  KRATOS_BREP_CURVE_3D_H_INCLUDED

// System includes

// External includes

// Project includes


namespace Kratos
{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/**
 * @class BrepCurve
 * @ingroup KratosCore
 * @brief The BrepCurve acts as topology for faces having a pointer
 *        to a surface geometry and a list of BrepEdge s. Those edges
 *        define the trimming of the respective surface.
 */
template<class TPointType> class BrepCurve
    : public Geometry<TPointType>
{
public:
    ///@}
    ///@name Type Definitions
    ///@{

    /**
     * Geometry as base class.
     */
    typedef Geometry<TPointType> BaseType;
    typedef Geometry<TPointType> GeometryType;

    /**
     * Pointer definition of BrepCurve
     */
    KRATOS_CLASS_POINTER_DEFINITION( BrepCurve );

    typedef TPointType PointType;

    typedef GeometryData::IntegrationMethod IntegrationMethod;

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
    typedef GeometryData::ShapeFunctionsDerivativesType ShapeFunctionsDerivativesType;
    typedef GeometryData::ShapeFunctionsDerivativesIntegrationPointsType ShapeFunctionsDerivativesIntegrationPointsType;
    typedef GeometryData::ShapeFunctionsDerivativesVectorType ShapeFunctionsDerivativesVectorType;

    // Own typedefs
    /**
    * Type of one brep loop.
    */
    typedef typename std::vector<BrepEdge::Pointer> BrepLoopType;


    ///@}
    ///@name Life Cycle
    ///@{

    BrepCurve( 
        typename NurbsSurface::Pointer pSurface,
        typename BrepLoopType InnerBrepLoopVector,
        typename BrepLoopType OuterBrepLoopVector)
        : BaseType(pSurface->mPointerVector(), nullptr ),
        mpNurbsSurface(pSurface),
        mInnerLoops(InnerBrepLoopVector),
        mOuterLoops(OuterBrepLoopVector)
    {
        mIsGeometryDataInitialized = false;
        mIsBoundaryPolygonInitialized = false;
        mIsSurfacePointCloudInitialized = false;
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
    BrepCurve( BrepCurve const& rOther )
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
    template<class TOtherPointType> explicit BrepCurve( BrepCurve<TOtherPointType> const& rOther )
        : BaseType( rOther )
    {
    }

    /**
     * Destructor. Does nothing!!!
     */
    ~BrepCurve() override {}

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
    BrepCurve& operator=( const BrepCurve& rOther )
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
    BrepCurve& operator=( BrepCurve<TOtherPointType> const & rOther )
    {
        BaseType::operator=( rOther );
        return *this;
    }

    ///@}
    ///@name Operations
    ///@{

    typename BaseType::Pointer Create( PointsArrayType const& ThisPoints ) const override
    {
        return typename BaseType::Pointer( new BrepCurve( ThisPoints ) );
    }

    ///@}

    /**
    * @brief Calculates the boundingbox of the geometry.
    * @details Corresponds with the highest and lowest point in space
    * @param rLowPoint  Lower point of the boundingbox.
    * @param rHighPoint Higher point of the boundingbox.
    */
    override void BoundingBox(
        TPointType& rLowPoint,
        TPointType& rHighPoint
    ) const
    {
        const SizeType dim = WorkingSpaceDimension();

        for (each point in mBoundaryPolygon and mSurfacePointCloud) { //The first node is already assigned, so we can start from 1
            const auto& r_point = this->GetPoint(point);
            for (IndexType i = 0; i < dim; ++i) {
                rHighPoint[i] = (rHighPoint[i] < r_point[i]) ? r_point[i] : rHighPoint[i];
                rLowPoint[i] = (rLowPoint[i]  > r_point[i]) ? r_point[i] : rLowPoint[i];
            }
        }
    }

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

    /**
    * Returns whether given arbitrary point is inside the Geometry and the respective
    * local point for the given global point
    * @param rPoint The point to be checked if is inside o note in global coordinates
    * @param rResult The local coordinates of the point
    * @param Tolerance The  tolerance that will be considered to check if the point is inside or not
    * @return True if the point is inside, false otherwise
    */
    virtual bool IsInside(
        const CoordinatesArrayType& rPoint,
        CoordinatesArrayType& rResult,
        const double Tolerance = std::numeric_limits<double>::epsilon()
    )
    {
        // boundarypolygon should be used. Boost can do that, there are more libraries...

        KRATOS_ERROR << "Calling base class IsInside method instead of derived class one. Please check the definition of derived class. " << *this << std::endl;
        return false;
    }

    ///@name Shape Function
    ///@{

    override ShapeFunctionsDerivativesType& ShapeFunctionsDerivative(
        IndexType DerivativeOrder,
        Vector &rResult,
        const CoordinatesArrayType& rCoordinates) const
    {
        return pSurface->ShapeFunctionsDerivative(DerivativeOrder, rCoordinates);
    }

    override ShapeFunctionsDerivativesVectorType& ShapeFunctionsDerivatives(
        IndexType DerivativeOrder,
        Vector &rResult,
        const CoordinatesArrayType& rCoordinates) const
    {
        return pSurface->ShapeFunctionsDerivatives(DerivativeOrder, rCoordinates);
    }

    ///@}
    ///@name Input and output
    ///@{

    std::vector<geometry::Pointer> GetIntegrationPointGeomtries()
    {
        if (!mIsGeometryDataInitialized)
            ComputeGeometryData();

        SizeType number_of_integration_points = mpGeometryData->IntegrationPointsNumber();

        IntegrationPointsArrayType integration_points = mpGeometryData->IntegrationPoints();
        ShapeFunctionsDerivativesVectorType shape_functions_derivatives_all = mpGeometryData->ShapeFunctionsDerivativesAll();

        SizeType number_of_derivatives = shape_functions_derivatives_all.size();

        std::vector<geometry> IntegrationPointList(number_of_integration_points);
        for (i = 0.0; i < number_of_integration_points; ++i)
        {
            ShapeFunctionsDerivativesVectorType shape_functions_derivatives_all_one_integration_point(number_of_derivatives);

            for (j = 0.0; j < number_of_integration_points; ++j)
            {
                shape_functions_derivatives_all_one_integration_point[j] = shape_functions_derivatives_all[j][i];
            }

            IntegrationPointSurface3d(integration_points[i], shape_functions_derivatives_all_one_integration_point);
        }
    }

    /**
     * Turn back information as a string.
     *
     * @return String contains information about this geometry.
     * @see PrintData()
     * @see PrintInfo()
     */
    std::string Info() const override
    {
        return "2 dimensional triangle with three nodes in 3D space";
    }

    /**
     * Print information about this object.
     * @param rOStream Stream to print into it.
     * @see PrintData()
     * @see Info()
     */
    void PrintInfo( std::ostream& rOStream ) const override
    {
        rOStream << "2 dimensional triangle with three nodes in 3D space";
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
    /**
     * :TODO: needs to be reviewed because it is not properly implemented yet
     * (comment by janosch)
     */
    void PrintData( std::ostream& rOStream ) const override
    {
        BaseType::PrintData( rOStream );
        std::cout << std::endl;
        Matrix jacobian;
        Jacobian( jacobian, PointType() );
        rOStream << "    Jacobian in the origin\t : " << jacobian;
    }

    ///@}
    ///@name Friends
    ///@{

    ///@}

protected:
    /**
     * There are no protected members in class BrepCurve
     */

    BPolygonGetBoundaryPolygon()
    {
        if (!mIsBoundaryPolygonInitialized) {
            this->InitializeBoundaryPolygon();
        }

            )
}

private:
    ///@name Static Member Variables
    ///@{

    static const GeometryData msGeometryData;

    NurbsSurfacePointer mpNurbsSurface;

    BrepLoopVectorType mInnerLoops;
    BrepLoopVectorType mOuterLoops;

    bool mIsGeometryDataInitialized;
    GeometryData mGeometryData;

    bool mIsBoundaryPolygonInitialized;
    BoundaryPolygon mBoundaryPolygon;

    bool mIsSurfacePointCloudInitialized;
    SurfaceCload mSurfacePointCloud;

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save( Serializer& rSerializer ) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, BaseType );
    }

    void load( Serializer& rSerializer ) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BaseType );
    }

    BrepCurve(): BaseType( PointsArrayType(), &mGeometryData ) {}

    ///@}
    ///@name Member Variables
    ///@{


    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    void ComputeGeometryData(IntegrationMethod)
    {
        if (mIsBoundaryPolygonInitialized)
            CreateBoundaryPolygon();

        //Make Clipper or equal

        //Create Integration Points

        int degree_u = pSurface->DegreeU();
        int degree_v = pSurface->DegreeV();

        int integration_degree = std::max(degree_u, degree_v) + 1;

        for (int i = 0; i < rClipper.NbSpansU(); ++i)
        {
            for (int j = 0; j < rClipper.NbSpansU(); ++j)
            {
                if (rClipper.SpanTrimType(i, j) == ANurbs::Empty)
                {
                    continue;
                }
                else if (rClipper.SpanTrimType(i, j) == ANurbs::Full)
                {
                    auto integration_points = ANurbs::IntegrationPoints<double>::Points2(
                        degree_u + 1,
                        degree_v + 1,
                        rClipper.SpanU(i),
                        rClipper.SpanV(j));

                    for (int i = 0; i < integration_points.size(); ++i)
                    {
                        mpNurbsSurface->Compute(
                            integration_point_polygon.IntegrationPoint(i).u,
                            integration_point_polygon.IntegrationPoint(i).v);
                    }
                }
                else if (rClipper.SpanTrimType(i, j) == ANurbs::Trimmed)
                {
                    auto polygons = rClipper.SpanPolygons(i, j);

                    for (int p = 0; p < polygons.size(); ++p)
                    {
                        auto integration_point_polygon = ANurbs::PolygonIntegrationPoints<Kratos::array_1d<double, 2>>();

                        integration_point_polygon.Compute(degree, polygons[p]);

                        for (int i = 0; i < integration_point_polygon.NbIntegrationPoints(); ++i)
                        {
                            mpNurbsSurface->Compute(
                                integration_point_polygon.IntegrationPoint(i).u,
                                integration_point_polygon.IntegrationPoint(i).v);

                            int number_of_non_zero_cps = shape.NonzeroPoleIndices().size();
                        }
                    }
                }
            }
        }
        mGeometryData.SetData(...);
    }


    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Private Friends
    ///@{

    template<class TOtherPointType> friend class BrepCurve;

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
    BrepCurve<TPointType>& rThis );
/**
 * output stream functions
 */
template<class TPointType> inline std::ostream& operator << (
    std::ostream& rOStream,
    const BrepCurve<TPointType>& rThis )
{
    rThis.PrintInfo( rOStream );
    rOStream << std::endl;
    rThis.PrintData( rOStream );
    return rOStream;
}

template<class TPointType>
const GeometryData BrepCurve<TPointType>::msGeometryData(1,
    3,
    1,
    GeometryData::GI_GAUSS_1,
    {},
    {},
    {});

///@}
}// namespace Kratos.

#endif // KRATOS_BREP_CURVE_3D_H_INCLUDED  defined

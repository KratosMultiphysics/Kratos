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

#if !defined(KRATOS_BREP_SURFACE_ON_VOLUME_H_INCLUDED )
#define  KRATOS_BREP_SURFACE_ON_VOLUME_H_INCLUDED

// System includes

// External includes

// Project includes
#include "geometries/geometry.h"
#include "geometries/nurbs_shape_function_utilities/nurbs_interval.h"

#include "geometries/nurbs_surface_on_volume_geometry.h"

#include "utilities/nurbs_utilities/projection_nurbs_geometry_utilities.h"

namespace Kratos
{
///@name Kratos Classes
///@{

/**
 * @class BrepSurfaceOnVolume
 * @ingroup KratosCore
 * @brief The BrepSurfaceOnVolume acts as topology for curves on surfaces.
 */
template<class TVolumeContainerPointType, class TContainerPointEmbeddedType = TVolumeContainerPointType>
class BrepSurfaceOnVolume
    : public Geometry<typename TVolumeContainerPointType::value_type>
{
public:
    ///@}
    ///@name Type Definitions
    ///@{

    /** Pointer definition of BrepSurfaceOnVolume */
    KRATOS_CLASS_POINTER_DEFINITION( BrepSurfaceOnVolume );

    typedef typename TVolumeContainerPointType::value_type PointType;

    typedef Geometry<typename TVolumeContainerPointType::value_type> BaseType;
    typedef Geometry<typename TVolumeContainerPointType::value_type> GeometryType;
    typedef typename GeometryType::Pointer GeometryPointer;

    typedef GeometryData::IntegrationMethod IntegrationMethod;

    // template<int TWorkingSpaceDimension, class TVolumeContainerPointType>

    typedef NurbsVolumeGeometry<TVolumeContainerPointType> NurbsVolumeType;

    typedef NurbsSurfaceGeometry<3, TVolumeContainerPointType> NurbsSurfaceType;

    typedef NurbsSurfaceOnVolumeGeometry<3,  TVolumeContainerPointType> NurbsSurfaceOnVolumeType;

    typedef typename NurbsSurfaceOnVolumeType::Pointer NurbsSurfaceOnVolumePointerType;

    typedef typename BaseType::GeometriesArrayType GeometriesArrayType;

    typedef typename BaseType::IndexType IndexType;
    typedef typename BaseType::SizeType SizeType;

    typedef typename BaseType::PointsArrayType PointsArrayType;
    typedef typename BaseType::CoordinatesArrayType CoordinatesArrayType;
    typedef typename BaseType::IntegrationPointsArrayType IntegrationPointsArrayType;

    static constexpr IndexType SURFACE_ON_VOLUME_INDEX = std::numeric_limits<IndexType>::max() - 2;

    ///@}
    ///@name Life Cycle
    ///@{

    /// constructor for untrimmed surface
    BrepSurfaceOnVolume(
        typename NurbsVolumeType::Pointer pVolume,
        typename NurbsSurfaceType::Pointer pSurface)
        : BaseType(PointsArrayType(), &msGeometryData)
        , mpSurfaceOnVolume(
            Kratos::make_shared<NurbsSurfaceOnVolumeType>(
                pVolume, pSurface))
        , mSurfaceNurbsIntervalU(pSurface->DomainIntervalU())
        , mSurfaceNurbsIntervalV(pSurface->DomainIntervalV())
    {
    }

    BrepSurfaceOnVolume()
        : BaseType(PointsArrayType(), &msGeometryData)
    {}

    explicit BrepSurfaceOnVolume(const PointsArrayType& ThisPoints)
        : BaseType(ThisPoints, &msGeometryData)
    {
    }

    /// Copy constructor.
    BrepSurfaceOnVolume( BrepSurfaceOnVolume const& rOther )
        : BaseType( rOther )
        , mpSurfaceOnVolume(rOther.mpSurfaceOnVolume)
        , mSurfaceNurbsIntervalU(rOther.mSurfaceNurbsIntervalU)
        , mSurfaceNurbsIntervalV(rOther.mSurfaceNurbsIntervalV)
    {
    }

    /// Copy constructor from a geometry with different point type.
    template<class TOtherContainerPointType, class TOtherContainerPointEmbeddedType>
    explicit BrepSurfaceOnVolume(
        BrepSurfaceOnVolume<TOtherContainerPointType, TOtherContainerPointEmbeddedType> const& rOther )
        : BaseType( rOther )
        , mpSurfaceOnVolume(rOther.mpSurfaceOnVolume)
        , mSurfaceNurbsIntervalU(rOther.mSurfaceNurbsIntervalU)
        , mSurfaceNurbsIntervalV(rOther.mSurfaceNurbsIntervalV)
    {
    }

    /// Destructor
    ~BrepSurfaceOnVolume() override = default;

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    BrepSurfaceOnVolume& operator=( const BrepSurfaceOnVolume& rOther )
    {
        BaseType::operator=( rOther );
        mpSurfaceOnVolume = rOther.mpSurfaceOnVolume;
        mSurfaceNurbsIntervalU = rOther.mSurfaceNurbsIntervalU;
        mSurfaceNurbsIntervalV = rOther.mSurfaceNurbsIntervalV;
        return *this;
    }

    /// Assignment operator for geometries with different point type.
    template<class TOtherContainerPointType, class TOtherContainerPointEmbeddedType>
    BrepSurfaceOnVolume& operator=( BrepSurfaceOnVolume<TOtherContainerPointType, TOtherContainerPointEmbeddedType> const & rOther )
    {
        BaseType::operator=( rOther );
        mpSurfaceOnVolume = rOther.mpSurfaceOnVolume;
        mSurfaceNurbsIntervalU = rOther.mSurfaceNurbsIntervalU;
        mSurfaceNurbsIntervalV = rOther.mSurfaceNurbsIntervalV;
        return *this;
    }

    ///@}
    ///@name Operations
    ///@{

    typename BaseType::Pointer Create( PointsArrayType const& ThisPoints ) const override
    {
        return typename BaseType::Pointer( new BrepSurfaceOnVolume( ThisPoints ) );
    }

    ///@}
    ///@name Access to Geometry Parts
    ///@{

    /**
    * @return pointer of geometry, corresponding to the index.
    */
    GeometryPointer pGetGeometryPart(const IndexType Index) override
    {
        const auto& const_this = *this;
        return std::const_pointer_cast<GeometryType>(
            const_this.pGetGeometryPart(Index));
    }


    const GeometryPointer pGetGeometryPart(const IndexType Index) const override
    {
        if (Index == GeometryType::BACKGROUND_GEOMETRY_INDEX)
            return mpSurfaceOnVolume->pGetGeometryPart(GeometryType::BACKGROUND_GEOMETRY_INDEX);

        if (Index == SURFACE_ON_VOLUME_INDEX)
            return mpSurfaceOnVolume;

        KRATOS_ERROR << "Index " << Index << " not existing in BrepSurfaceOnVolume: "
            << this->Id() << std::endl;
    }

    /**
    * @brief This function is used to check if the index is either
    *        GeometryType::BACKGROUND_GEOMETRY_INDEX or SURFACE_ON_VOLUME_INDEX.
    * @param Index of the geometry part.
    * @return true if GeometryType::BACKGROUND_GEOMETRY_INDEX or SURFACE_ON_VOLUME_INDEX.
    */
    bool HasGeometryPart(const IndexType Index) const override
    {
        return (Index == GeometryType::BACKGROUND_GEOMETRY_INDEX || Index == SURFACE_ON_VOLUME_INDEX);
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
            mpSurfaceOnVolume->Calculate(rVariable, rOutput);
        }
    }

    ///@}
    ///@name Mathematical Informations
    ///@{

    SizeType PolynomialDegree(IndexType LocalDirectionIndex) const override
    {
        return mpSurfaceOnVolume->PolynomialDegree(LocalDirectionIndex);
    }

    ///@}
    ///@name Set/ Get functions
    ///@{


    NurbsInterval DomainIntervalU() const
    {
        return mSurfaceNurbsIntervalU;
    }
    NurbsInterval DomainIntervalV() const
    {
        return mSurfaceNurbsIntervalV;
    }


    /// Returns the NurbsSurfaceOnVolume::Pointer of this brep.
    NurbsSurfaceOnVolumePointerType pGetSurfaceOnVolume()
    {
        return mpSurfaceOnVolume;
    }

    /// Returns the const NurbsSurfaceOnVolume::Pointer of this brep.
    const NurbsSurfaceOnVolumePointerType pGetSurfaceOnVolume() const
    {
        return mpSurfaceOnVolume;
    }

    /// Returns number of points of NurbsSurfaceOnVolume.
    SizeType PointsNumberInDirection(IndexType DirectionIndex) const override
    {
        return mpSurfaceOnVolume->PointsNumberInDirection(DirectionIndex);
    }


    /* @brief Provides intersections of the nurbs curve with the knots of the surface,
     *         using the interval of this curve.
     * @param vector of span intervals.
     * @param index of chosen direction, for curves always 0.
     */

    void SpansLocalSpace(std::vector<double>& rSpansU, std::vector<double>& rSpansV) const
    {
        
        mpSurfaceOnVolume->SpansLocalSpace(rSpansU, rSpansV,
            mSurfaceNurbsIntervalU.GetT0(), mSurfaceNurbsIntervalU.GetT1(),
            mSurfaceNurbsIntervalV.GetT0(), mSurfaceNurbsIntervalV.GetT1());
    }

    /* @brief checks and returns if local coordinate rPointLocalCoordinates[0]
     *        is inside the local/parameter space.
     * @return on boundary -> 2 - meaning that it is equal to start or end point.
     *         inside -> 1
     *         outside -> 0
     */

    // int IsInsideLocalSpace(
    //     const CoordinatesArrayType& rPointLocalCoordinates,
    //     const double Tolerance = std::numeric_limits<double>::epsilon()
    // ) const override
    // {
    //     const double min_parameter = mSurfaceNurbsInterval.MinParameter();
    //     if (rPointLocalCoordinates[0] < min_parameter) {
    //         return 0;
    //     } else if (std::abs(rPointLocalCoordinates[0] - min_parameter) < Tolerance) {
    //         return 2;
    //     }

    //     const double max_parameter = mSurfaceNurbsInterval.MaxParameter();
    //     if (rPointLocalCoordinates[0] > max_parameter) {
    //         return 0;
    //     } else if (std::abs(rPointLocalCoordinates[0] - max_parameter) < Tolerance) {
    //         return 2;
    //     }

    //     return 1;
    // }

    ///@}
    ///@name ClosestPoint
    ///@{

    /* @brief Makes a check if the provided paramater rPointLocalCoordinates[0]
     *        is inside the curve, or on the boundary or if it lays outside.
     *        If it is outside, it is set to the boundary which is closer to it.
     * @return if rPointLocalCoordinates[0] was before the projection:
     *         inside -> 1
     *         outside -> 0
     *         on boundary -> 2 - meaning that it is equal to start or end point.
     */
    // virtual int ClosestPointLocalToLocalSpace(
    //     const CoordinatesArrayType& rPointLocalCoordinates,
    //     CoordinatesArrayType& rClosestPointLocalCoordinates,
    //     const double Tolerance = std::numeric_limits<double>::epsilon()
    // ) const override
    // {
    //     const double min_parameter = mSurfaceNurbsInterval.MinParameter();
    //     if (rPointLocalCoordinates[0] < min_parameter) {
    //         rClosestPointLocalCoordinates[0] = min_parameter;
    //         return 0;
    //     } else if (std::abs(rPointLocalCoordinates[0] - min_parameter) < Tolerance) {
    //         rClosestPointLocalCoordinates[0] = rPointLocalCoordinates[0];
    //         return 2;
    //     }

    //     const double max_parameter = mSurfaceNurbsInterval.MaxParameter();
    //     if (rPointLocalCoordinates[0] > max_parameter) {
    //         rClosestPointLocalCoordinates[0] = min_parameter;
    //         return 0;
    //     } else if (std::abs(rPointLocalCoordinates[0] - max_parameter) < Tolerance) {
    //         rClosestPointLocalCoordinates[0] = rPointLocalCoordinates[0];
    //         return 2;
    //     }

    //     rClosestPointLocalCoordinates[0] = rPointLocalCoordinates[0];
    //     return 1;
    // }

    ///@}
    ///@name Projection
    ///@{

     /* @brief Makes projection of rPointGlobalCoordinates to
      *       the closest point on the curve, with
      *       local coordinates rProjectedPointLocalCoordinates.
      *
      * @param Tolerance is the breaking criteria.
      * @return 1 -> projection succeeded
      *         0 -> projection failed
      */
    // int ProjectionPointGlobalToLocalSpace(
    //     const CoordinatesArrayType& rPointGlobalCoordinates,
    //     CoordinatesArrayType& rProjectedPointLocalCoordinates,
    //     const double Tolerance = std::numeric_limits<double>::epsilon()
    // ) const override
    // {
    //     CoordinatesArrayType point_global_coordinates;

    //     return ProjectionNurbsGeometryUtilities::NewtonRaphsonSurface(
    //         rProjectedPointLocalCoordinates,
    //         rPointGlobalCoordinates,
    //         point_global_coordinates,
    //         *this,
    //         20, Tolerance);
    // }

    ///@}
    ///@name Geometrical Operations
    ///@{

    /// Provides the center of the underlying curve on surface
    Point Center() const override
    {
        return mpSurfaceOnVolume->Center();
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
        return mpSurfaceOnVolume->GlobalCoordinates(rResult, rLocalCoordinates);
    }

    /**
    * @brief This method maps from local space to global/working space and computes the
    *        number of derivatives at the underlying nurbs curve on surface
    *        at the parameter rLocalCoordinates[0].
    *
    * @param LocalCoordinates The local coordinates in paramater space
    * @param Derivative Number of computed derivatives
    *        0 -> Location = PointLocalCoordinates
    *        1 -> Tangent
    *        2 -> Curvature
    *        ...
    * @return std::vector<array_1d<double, 3>> with the global space derivatives
    */
    void GlobalSpaceDerivatives(
        std::vector<CoordinatesArrayType>& rGlobalSpaceDerivatives,
        const CoordinatesArrayType& rLocalCoordinates,
        const SizeType DerivativeOrder) const override
    {
        return mpSurfaceOnVolume->GlobalSpaceDerivatives(rGlobalSpaceDerivatives, rLocalCoordinates, DerivativeOrder);
    }

    /**
    * Returns whether given arbitrary point is inside the Geometry and the respective
    * local point for the given global point
    * @param rPoint The point to be checked if is inside o note in global coordinates
    * @param rResult The local coordinates of the point
    * @param Tolerance The  tolerance that will be considered to check if the point is inside or not
    * @return True if the point is inside, false otherwise
    */
    bool IsInside(
        const CoordinatesArrayType& rPoint,
        CoordinatesArrayType& rResult,
        const double Tolerance = std::numeric_limits<double>::epsilon()
    ) const override
    {
        KRATOS_ERROR << "IsInside is not yet implemented within the BrepSurfaceOnVolume";
    }

    ///@}
    ///@name Geometrical Informations
    ///@{

    /// Computes the length of a nurbs curve
    double Length() const override
    {
        IntegrationPointsArrayType integration_points;
        IntegrationInfo integration_info = GetDefaultIntegrationInfo();
        CreateIntegrationPoints(integration_points, integration_info);

        double length = 0.0;
        for (IndexType i = 0; i < integration_points.size(); ++i) {
            const double determinant_jacobian = mpSurfaceOnVolume->DeterminantOfJacobian(integration_points[i]);
            length += integration_points[i].Weight() * determinant_jacobian;
        }
        return length;
    }

    ///@}
    ///@name Integration Info
    ///@{

    /// Provides the default integration dependent on the polynomial degree.
    IntegrationInfo GetDefaultIntegrationInfo() const override
    {
        return mpSurfaceOnVolume->GetDefaultIntegrationInfo();
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
        std::vector<double> spansU; std::vector<double> spansV;
        SpansLocalSpace(spansU, spansV);

        IntegrationPointUtilities::CreateIntegrationPoints2D(
            rIntegrationPoints, spansU, spansV, rIntegrationInfo);

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
        const IntegrationPointsArrayType& rIntegrationPoints,
        IntegrationInfo& rIntegrationInfo) override
    {
        mpSurfaceOnVolume->CreateQuadraturePointGeometries(
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
        return mpSurfaceOnVolume->ShapeFunctionsValues(rResult, rCoordinates);
    }

    Matrix& ShapeFunctionsLocalGradients(
        Matrix& rResult,
        const CoordinatesArrayType& rCoordinates) const override
    {
        return mpSurfaceOnVolume->ShapeFunctionsLocalGradients(rResult, rCoordinates);
    }

    ///@}
    ///@name Geometry Classification
    ///@{

    GeometryData::KratosGeometryFamily GetGeometryFamily() const override
    {
        return GeometryData::KratosGeometryFamily::Kratos_Brep;
    }

    GeometryData::KratosGeometryType GetGeometryType() const override
    {
        return GeometryData::KratosGeometryType::Kratos_Brep_Surface_On_Volume;
    }

    ///@}
    ///@name Information
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "Brep face surface";
    }

    /// Print information about this object.
    void PrintInfo( std::ostream& rOStream ) const override
    {
        rOStream << "Brep face surface";
    }

    /// Print object's data.
    void PrintData( std::ostream& rOStream ) const override
    {
        BaseType::PrintData( rOStream );
        std::cout << std::endl;
        rOStream << "    Brep face surface : " << std::endl;
    }

    void SetIsExitingDirectionSBM(int isExitingDirection) 
    {
        /* isExitingDirection :
                            -1 -> entering surface/external body fitted surface
                             0 -> exiting surface in x direction
                             1 -> exiting surface in y direction
                             2 -> exiting surface in z direction
                            */
        KRATOS_ERROR_IF(isExitingDirection > 2 || isExitingDirection < -1) 
                     << "WRONG EXITING DIRECTION: not external face and not parallel to any direction";
        mpSurfaceOnVolume->SetIsExitingDirectionSBM(isExitingDirection);
    }

    void SetNormalSBM(Vector &Normal) 
    {
        mpSurfaceOnVolume->SetNormalSBM(Normal);
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

    NurbsSurfaceOnVolumePointerType mpSurfaceOnVolume;

    NurbsInterval mSurfaceNurbsIntervalU;

    NurbsInterval mSurfaceNurbsIntervalV;


    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save( Serializer& rSerializer ) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, BaseType );
        rSerializer.save("SurfaceOnVolume", mpSurfaceOnVolume);
        rSerializer.save("NurbsIntervalU", mSurfaceNurbsIntervalU);
        rSerializer.save("NurbsIntervalV", mSurfaceNurbsIntervalV);
    }

    void load( Serializer& rSerializer ) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BaseType );
        rSerializer.load("SurfaceOnVolume", mpSurfaceOnVolume);
        rSerializer.load("NurbsIntervalU", mSurfaceNurbsIntervalU);
        rSerializer.load("NurbsIntervalV", mSurfaceNurbsIntervalV);
    }

    ///@}
}; // Class BrepSurfaceOnVolume

///@name Input and output
///@{

/// input stream functions
template<class TContainerPointType, class TContainerPointEmbeddedType = TContainerPointType> inline std::istream& operator >> (
    std::istream& rIStream,
    BrepSurfaceOnVolume<TContainerPointType, TContainerPointEmbeddedType>& rThis );

/// output stream functions
template<class TContainerPointType, class TContainerPointEmbeddedType = TContainerPointType> inline std::ostream& operator << (
    std::ostream& rOStream,
    const BrepSurfaceOnVolume<TContainerPointType, TContainerPointEmbeddedType>& rThis )
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
GeometryData BrepSurfaceOnVolume<TContainerPointType, TContainerPointEmbeddedType>::msGeometryData(
    &msGeometryDimension,
    GeometryData::IntegrationMethod::GI_GAUSS_1,
    {}, {}, {});

template<class TContainerPointType, class TContainerPointEmbeddedType>
const GeometryDimension BrepSurfaceOnVolume<TContainerPointType, TContainerPointEmbeddedType>::msGeometryDimension(3, 1);

///@}
}// namespace Kratos.

#endif // KRATOS_BREP_CURVE_ON_SURFACE_3D_H_INCLUDED  defined

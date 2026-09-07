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

#pragma once

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
 * @tparam TShiftedBoundary Boolean flag indicating whether is defined with shifted boundary conditions.
 */
template<class TContainerPointType, bool TShiftedBoundary, class TContainerPointEmbeddedType = TContainerPointType >
class BrepSurfaceOnVolume
    : public Geometry<typename TContainerPointType::value_type>
{
public:
    ///@}
    ///@name Type Definitions
    ///@{

    /** Pointer definition of BrepSurfaceOnVolume */
    KRATOS_CLASS_POINTER_DEFINITION( BrepSurfaceOnVolume );

    typedef typename TContainerPointType::value_type PointType;

    typedef Geometry<typename TContainerPointType::value_type> BaseType;
    typedef Geometry<typename TContainerPointType::value_type> GeometryType;
    typedef typename GeometryType::Pointer GeometryPointer;

    typedef GeometryData::IntegrationMethod IntegrationMethod;

    typedef NurbsVolumeGeometry<TContainerPointType> NurbsVolumeType;
    typedef NurbsSurfaceGeometry<3, TContainerPointEmbeddedType> NurbsSurfaceType;

    typedef NurbsSurfaceOnVolumeGeometry<3, TContainerPointType> NurbsSurfaceOnVolumeType;

    typedef typename NurbsSurfaceOnVolumeType::Pointer NurbsSurfaceOnVolumePointerType;

    typedef typename BaseType::GeometriesArrayType GeometriesArrayType;

    typedef typename BaseType::IndexType IndexType;
    typedef typename BaseType::SizeType SizeType;

    typedef typename BaseType::PointsArrayType PointsArrayType;
    typedef typename BaseType::CoordinatesArrayType CoordinatesArrayType;
    typedef typename BaseType::IntegrationPointsArrayType IntegrationPointsArrayType;

    static constexpr IndexType SURFACE_ON_VOLUME_INDEX = std::numeric_limits<IndexType>::max() - 2;

    using GeometryType::SpansLocalSpace;

    ///@}
    ///@name Life Cycle
    ///@{

    /// constructor for untrimmed volume
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
        BrepSurfaceOnVolume<TOtherContainerPointType, TShiftedBoundary, TOtherContainerPointEmbeddedType> const& rOther )
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
    BrepSurfaceOnVolume& operator=( BrepSurfaceOnVolume<TOtherContainerPointType, TShiftedBoundary, TOtherContainerPointEmbeddedType> const & rOther )
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
    ///@name Mathematical Informations
    ///@{

    /// Return polynomial degree of the nurbs curve on volume
    SizeType PolynomialDegree(IndexType LocalDirectionIndex) const override
    {
        return mpSurfaceOnVolume->PolynomialDegree(LocalDirectionIndex);
    }

    ///@}
    ///@name Set/ Get functions
    ///@{

    /* @brief Provides the natural boundaries of the NURBS/B-Spline curve.
     * @return domain interval.
     */
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

    ///@}
    ///@name Surface Properties
    ///@{

    /* @brief Provides intersections of the nurbs curve with the knots of the volume,
     *         using the interval of this curve.
     * @param vector of span intervals.
     * @param index of chosen direction, for curves always 0.
     */
    void SpansLocalSpace(std::vector<double>& rSpansU, std::vector<double>& rSpansV) const
    {
        /* When the `TShiftedBoundary` template parameter is true, the spans are computed using 
        *  the shifted boundary method by invoking `SpansLocalSpaceSBM`.
        *  Otherwise, the spans are computed using the standard method by invoking `SpansLocalSpace`
        */
        if constexpr (TShiftedBoundary) {
            
            mpSurfaceOnVolume->SpansLocalSpace(rSpansU, rSpansV,
                mSurfaceNurbsIntervalU.GetT0(), mSurfaceNurbsIntervalU.GetT1(),
                mSurfaceNurbsIntervalV.GetT0(), mSurfaceNurbsIntervalV.GetT1());

        } else {
            mpSurfaceOnVolume->SpansLocalSpace(rSpansU, rSpansV,
                mSurfaceNurbsIntervalU.GetT0(), mSurfaceNurbsIntervalU.GetT1(),
                mSurfaceNurbsIntervalV.GetT0(), mSurfaceNurbsIntervalV.GetT1());
        }
        
    }

    ///@}
    ///@name Geometrical Operations
    ///@{

    /// Provides the center of the underlying curve on volume
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
    *        number of derivatives at the underlying nurbs curve on volume
    *        at the parameter rLocalCoordinates[0].
    *
    * @param LocalCoordinates The local coordinates in parameter space
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

    /* Creates integration points on the nurbs volume of this geometry.
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
     *        curve on volume of this geometry.
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
        // If `TShiftedBoundary` is true, the method `CreateQuadraturePointGeometriesSbm` is called.
        // Otherwise, the method `CreateQuadraturePointGeometries` is used for standard processing.
        if constexpr (TShiftedBoundary) {
            mpSurfaceOnVolume->CreateQuadraturePointGeometriesSbm(
                rResultGeometries, NumberOfShapeFunctionDerivatives, rIntegrationPoints, rIntegrationInfo);
        } else {
            mpSurfaceOnVolume->CreateQuadraturePointGeometries(
                rResultGeometries, NumberOfShapeFunctionDerivatives, rIntegrationPoints, rIntegrationInfo);
        }
        
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
     * @details This function returns the type of the geometry. The geometry type provides a more specific classification within the geometry family, allowing for finer categorization and processing.
     * @return GeometryData::KratosGeometryType The geometry type.
     */
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

    void SetIsExitingDirectionSbm(int isExitingDirection) 
    {
        /* isExitingDirection :
                            -1 -> entering surface/external body fitted surface
                             0 -> exiting surface in x direction
                             1 -> exiting surface in y direction
                             2 -> exiting surface in z direction
                            */
        KRATOS_ERROR_IF(isExitingDirection > 2 || isExitingDirection < -1) 
                     << "WRONG EXITING DIRECTION: not external face and not parallel to any direction";
        mpSurfaceOnVolume->SetIsExitingDirectionSbm(isExitingDirection);  
    }

    void SetNormalSbm(Vector &Normal) 
    {
        mpSurfaceOnVolume->SetNormalSbm(Normal);
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
template<class TContainerPointType, bool TShiftedBoundary, class TContainerPointEmbeddedType = TContainerPointType> inline std::istream& operator >> (
    std::istream& rIStream,
    BrepSurfaceOnVolume<TContainerPointType, TShiftedBoundary, TContainerPointEmbeddedType>& rThis );

/// output stream functions
template<class TContainerPointType, bool TShiftedBoundary, class TContainerPointEmbeddedType = TContainerPointType> inline std::ostream& operator << (
    std::ostream& rOStream,
    const BrepSurfaceOnVolume<TContainerPointType, TShiftedBoundary, TContainerPointEmbeddedType>& rThis )
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
GeometryData BrepSurfaceOnVolume<TContainerPointType, TShiftedBoundary, TContainerPointEmbeddedType>::msGeometryData(
    &msGeometryDimension,
    GeometryData::IntegrationMethod::GI_GAUSS_1,
    {}, {}, {});

template<class TContainerPointType, bool TShiftedBoundary, class TContainerPointEmbeddedType>
const GeometryDimension BrepSurfaceOnVolume<TContainerPointType, TShiftedBoundary, TContainerPointEmbeddedType>::msGeometryDimension(3, 1);

///@}
}// namespace Kratos.
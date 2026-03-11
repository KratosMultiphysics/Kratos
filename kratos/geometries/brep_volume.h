//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Nicolò Antonelli
//                   Andrea Gorgi
//

#pragma once

// System includes

// External includes

// Project includes
#include "geometries/geometry.h"
#include "geometries/brep_surface_on_volume.h"
#include "geometries/nurbs_shape_function_utilities/nurbs_interval.h"

// integration
#include "utilities/geometry_utilities/brep_sbm_utilities.h"

namespace Kratos
{
///@name Kratos Classes
///@{

/**
 * @class BrepVolume
 * @ingroup KratosCore
 * @brief The BrepVolume acts as topology for faces. Those
 *        can be enclosed by a certain set of brep face curves.
 * @tparam TShiftedBoundary Boolean flag indicating whether is 
 *        defined with shifted boundary conditions.
 */
template<class TContainerPointType, bool TShiftedBoundary, class TContainerPointEmbeddedType = TContainerPointType>
class BrepVolume
    : public Geometry<typename TContainerPointType::value_type>
{
public:
    ///@}
    ///@name Type Definitions
    ///@{

    /** Pointer definition of BrepVolume */
    KRATOS_CLASS_POINTER_DEFINITION( BrepVolume );

    typedef typename TContainerPointType::value_type PointType;

    typedef Geometry<typename TContainerPointType::value_type> BaseType;
    typedef Geometry<typename TContainerPointType::value_type> GeometryType;
    typedef typename GeometryType::Pointer GeometryPointer;
    using GeometrySurrogateArrayType = DenseVector<GeometryPointer>;

    typedef GeometryData::IntegrationMethod IntegrationMethod;

    typedef NurbsVolumeGeometry<TContainerPointType> NurbsVolumeType;
    typedef BrepSurfaceOnVolume<TContainerPointType, TShiftedBoundary, TContainerPointEmbeddedType> BrepSurfaceOnVolumeType;

    typedef BrepSbmUtilities<Node> BrepSbmUtilitiesType;

    typedef DenseVector<typename BrepSurfaceOnVolumeType::Pointer> BrepSurfaceOnVolumeArrayType;
    typedef DenseVector<typename BrepSurfaceOnVolumeType::Pointer> BrepSurfaceOnVolumeLoopType;
    typedef DenseVector<DenseVector<typename BrepSurfaceOnVolumeType::Pointer>> BrepSurfaceOnVolumeLoopArrayType;

    typedef typename BaseType::GeometriesArrayType GeometriesArrayType;

    typedef typename BaseType::IndexType IndexType;
    typedef typename BaseType::SizeType SizeType;

    typedef typename BaseType::PointsArrayType PointsArrayType;
    typedef typename BaseType::CoordinatesArrayType CoordinatesArrayType;
    typedef typename BaseType::IntegrationPointsArrayType IntegrationPointsArrayType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor for untrimmed patch
    BrepVolume(
        typename NurbsVolumeType::Pointer pVolume)
        : BaseType(PointsArrayType(), &msGeometryData)
        , mpNurbsVolume(pVolume)
    {
        mIsTrimmed = false;
    }

    /// Constructor for trimmed patch
    BrepVolume(
        typename NurbsVolumeType::Pointer pVolume,
        BrepSurfaceOnVolumeLoopArrayType& BrepOuterLoopArray,
        BrepSurfaceOnVolumeLoopArrayType& BrepInnerLoopArray)
        : BaseType(PointsArrayType(), &msGeometryData)
        , mpNurbsVolume(pVolume)
        , mOuterLoopArray(BrepOuterLoopArray)
        , mInnerLoopArray(BrepInnerLoopArray)
    {
        mIsTrimmed = !(mOuterLoopArray.size() == 0 && mInnerLoopArray.size() == 0);
    }

    /// Constructor for trimmed patch including IsTrimmed
    BrepVolume(
        typename NurbsVolumeType::Pointer pVolume,
        BrepSurfaceOnVolumeLoopArrayType& BrepOuterLoopArray,
        BrepSurfaceOnVolumeLoopArrayType& BrepInnerLoopArray,
        bool IsTrimmed)
        : BaseType(PointsArrayType(), &msGeometryData)
        , mpNurbsVolume(pVolume)
        , mOuterLoopArray(BrepOuterLoopArray)
        , mInnerLoopArray(BrepInnerLoopArray)
        , mIsTrimmed(IsTrimmed)
    {
    }

    explicit BrepVolume(const PointsArrayType& ThisPoints)
        : BaseType(ThisPoints, &msGeometryData)
    {
    }

    /// Copy constructor.
    BrepVolume( BrepVolume const& rOther )
        : BaseType( rOther )
        , mpNurbsVolume(rOther.mpNurbsVolume)
        , mOuterLoopArray(rOther.mOuterLoopArray)
        , mInnerLoopArray(rOther.mInnerLoopArray)
        , mEmbeddedEdgesArray(rOther.mEmbeddedEdgesArray)
        , mIsTrimmed(rOther.mIsTrimmed)
    {
    }

    /// Destructor
    ~BrepVolume() override = default;

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    BrepVolume& operator=( const BrepVolume& rOther )
    {
        BaseType::operator=( rOther );
        mpNurbsVolume = rOther.mpNurbsVolume;
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
    BrepVolume& operator=( BrepVolume<TOtherContainerPointType, TShiftedBoundary, TOtherContainerPointEmbeddedType> const & rOther )
    {
        BaseType::operator=( rOther );
        mpNurbsVolume = rOther.mpNurbsVolume;
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
        return typename BaseType::Pointer( new BrepVolume( ThisPoints ) );
    }


    ///@}
    ///@name Mathematical Informations
    ///@{

    /// Return polynomial degree of the nurbs surface
    SizeType PolynomialDegree(IndexType LocalDirectionIndex) const override
    {
        return mpNurbsVolume->PolynomialDegree(LocalDirectionIndex);
    }

    ///@}
    ///@name Geometrical Operations
    ///@{

    /// Provides the center of the underlying surface
    Point Center() const override
    {
        return mpNurbsVolume->Center();
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
        return mpNurbsVolume->ProjectionPointGlobalToLocalSpace(
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
        mpNurbsVolume->GlobalCoordinates(rResult, rLocalCoordinates);

        return rResult;
    }

    ///@}
    ///@name Integration Info
    ///@{

    /// Provides the default integration dependent on the polynomial degree.
    IntegrationInfo GetDefaultIntegrationInfo() const override
    {
        return mpNurbsVolume->GetDefaultIntegrationInfo();
    }

    ///@}
    ///@name Integration Points
    ///@{

    /* Creates integration points on the nurbs volume of this geometry.
     * Accounting for whether the surface is trimmed or untrimmed, and whether shifted 
     * boundary conditions are used
     * 
     * - **Untrimmed Volume**: -> Non-cutting case
     *   - If `TShiftedBoundary` is true, the method prepares for the shifted boundary method (SBM) 
     *   - Otherwise, it directly uses `CreateIntegrationPoints` from the underlying NURBS surface.
     * - **Trimmed Volume**: -> Cutting case
     *   - It calls `BrepTrimmingUtilities::CreateBrepVolumeTrimmingIntegrationPoints` 
     *     to generate integration points that conform to the trimming surface.
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
                std::vector<double> spans_w;
                mpNurbsVolume->SpansLocalSpace(spans_u, 0);
                mpNurbsVolume->SpansLocalSpace(spans_v, 1);
                mpNurbsVolume->SpansLocalSpace(spans_w, 2);
                
                // Call  "BrepSBMUtilities::CreateBrepVolumeSbmIntegrationPoints"
                BrepSbmUtilitiesType::CreateBrepVolumeSbmIntegrationPoints(
                    spans_u, 
                    spans_v,
                    spans_w,
                    *mpSurrogateOuterLoopGeometries,
                    *mpSurrogateInnerLoopGeometries,
                    rIntegrationPoints,
                    rIntegrationInfo); 
            }
            // body-fitted case
            else {
                mpNurbsVolume->CreateIntegrationPoints(
                    rIntegrationPoints, rIntegrationInfo);
            }
        }
        // trimmed case
        else
        {
            KRATOS_ERROR << "Trimmed BrepVolume is not yet implemented" << std::endl;
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
        mpNurbsVolume->CreateQuadraturePointGeometries(
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
        mpNurbsVolume->ShapeFunctionsValues(rResult, rCoordinates);

        return rResult;
    }

    Matrix& ShapeFunctionsLocalGradients(
        Matrix& rResult,
        const CoordinatesArrayType& rCoordinates) const override
    {
        mpNurbsVolume->ShapeFunctionsLocalGradients(rResult, rCoordinates);

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
        return GeometryData::KratosGeometryType::Kratos_Brep_Volume;
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
        return "Brep volume";
    }

    /// Print information about this object.
    void PrintInfo( std::ostream& rOStream ) const override
    {
        rOStream << "Brep volume";
    }

    /// Print object's data.
    void PrintData( std::ostream& rOStream ) const override
    {
        BaseType::PrintData( rOStream );
        std::cout << std::endl;
        rOStream << "    Brep volume " << std::endl;
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

    typename NurbsVolumeType::Pointer mpNurbsVolume;

    BrepSurfaceOnVolumeLoopArrayType mOuterLoopArray;
    BrepSurfaceOnVolumeLoopArrayType mInnerLoopArray;

    BrepSurfaceOnVolumeArrayType mEmbeddedEdgesArray;

    // For Sbm
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
        rSerializer.save("NurbsSurface", mpNurbsVolume);
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
        rSerializer.load("NurbsSurface", mpNurbsVolume);
        rSerializer.load("OuterLoopArray", mOuterLoopArray);
        rSerializer.load("InnerLoopArray", mInnerLoopArray);
        rSerializer.load("EmbeddedEdgesArray", mEmbeddedEdgesArray);
        rSerializer.load("IsTrimmed", mIsTrimmed);
        rSerializer.load("SurrogateInnerLoopGeometries", mpSurrogateInnerLoopGeometries);
        rSerializer.load("SurrogateOuterLoopGeometries", mpSurrogateOuterLoopGeometries);
    }

    BrepVolume()
        : BaseType( PointsArrayType(), &msGeometryData )
    {}

    ///@}

}; // Class BrepVolume

///@name Input and output
///@{

/// input stream functions
template<class TContainerPointType, bool TShiftedBoundary, class TContainerPointEmbeddedType = TContainerPointType> inline std::istream& operator >> (
    std::istream& rIStream,
    BrepVolume<TContainerPointType, TShiftedBoundary, TContainerPointEmbeddedType>& rThis );

/// output stream functions
template<class TContainerPointType, bool TShiftedBoundary, class TContainerPointEmbeddedType = TContainerPointType> inline std::ostream& operator << (
    std::ostream& rOStream,
    const BrepVolume<TContainerPointType, TShiftedBoundary, TContainerPointEmbeddedType>& rThis )
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
GeometryData BrepVolume<TContainerPointType, TShiftedBoundary, TContainerPointEmbeddedType>::msGeometryData(
    &msGeometryDimension,
    GeometryData::IntegrationMethod::GI_GAUSS_1,
    {}, {}, {});

template<class TContainerPointType, bool TShiftedBoundary, class TContainerPointEmbeddedType>
const GeometryDimension BrepVolume<TContainerPointType, TShiftedBoundary, TContainerPointEmbeddedType>::msGeometryDimension(3, 2);

///@}
}// namespace Kratos.
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

#if !defined(KRATOS_BREP_VOLUME_H_INCLUDED )
#define  KRATOS_BREP_VOLUME_H_INCLUDED

// System includes

// External includes

// Project includes
#include "geometries/geometry.h"
#include "geometries/brep_surface_on_volume.h"
#include "geometries/nurbs_shape_function_utilities/nurbs_interval.h"

// trimming integration
// #include "utilities/geometry_utilities/brep_trimming_utilities.h"

// SBM integration
#include "utilities/geometry_utilities/brep_sbm_utilities.h"

namespace Kratos
{
///@name Kratos Classes
///@{

/**
 * @class BrepVolume
 * @ingroup KratosCore
 * @brief The BrepVolume acts as topology for volumes. Those
 *        can be enclosed by a certain set of brep volume faces.
 */
template<class TContainerPointType, class TContainerPointEmbeddedType = TContainerPointType>
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

    typedef GeometryData::IntegrationMethod IntegrationMethod;

    typedef NurbsVolumeGeometry<TContainerPointType> NurbsVolumeType;
    typedef BrepSurfaceOnVolume<TContainerPointType, TContainerPointEmbeddedType> BrepSurfaceOnVolumeType;

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


    // Constructor for SBM
    /// Constructor for trimmed patch
    BrepVolume(
        typename NurbsVolumeType::Pointer pVolume, 
        ModelPart& rSurrogateModelPart_inner, ModelPart& rSurrogateModelPart_outer)
        : BaseType(PointsArrayType(), &msGeometryData)
        , mpNurbsVolume(pVolume) 
        , mpSurrogateModelPart_inner(&rSurrogateModelPart_inner)
        , mpSurrogateModelPart_outer(&rSurrogateModelPart_outer)
    {
    }

    explicit BrepVolume(const PointsArrayType& ThisPoints)
        : BaseType(ThisPoints, &msGeometryData)
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
        return *this;
    }

    /// Assignment operator for geometries with different point type.
    template<class TOtherContainerPointType, class TOtherContainerPointEmbeddedType>
    BrepVolume& operator=( BrepVolume<TOtherContainerPointType, TOtherContainerPointEmbeddedType> const & rOther )
    {
        BaseType::operator=( rOther );
        mpNurbsVolume = rOther.mpNurbsVolume;
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
    *        Surface of the geometry is accessable with GeometryType::BACKGROUND_GEOMETRY_INDEX.
    * @param Index: trim_index or GeometryType::BACKGROUND_GEOMETRY_INDEX.
    * @return pointer of geometry, corresponding to the index.
    */
    const GeometryPointer pGetGeometryPart(const IndexType Index) const override
    {
        if (Index == GeometryType::BACKGROUND_GEOMETRY_INDEX)
            return mpNurbsVolume;

        KRATOS_ERROR << "Index " << Index << " not existing in BrepVolume: "
            << this->Id() << std::endl;
    }

    /**
    * @brief This function is used to check if this BrepVolume
    *        has certain trim or surface object.
    * @param Index of the geometry part.
    * @return true if has trim or surface
    */
    bool HasGeometryPart(const IndexType Index) const override
    {
        if (Index == GeometryType::BACKGROUND_GEOMETRY_INDEX)
            return true;

        return false;
    }


    /// Return polynomial degree of the nurbs surface
    SizeType PolynomialDegree(IndexType LocalDirectionIndex) const override
    {
        return mpNurbsVolume->PolynomialDegree(LocalDirectionIndex);
    }

    ///@}
    ///@name Information
    ///@{


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

    /* Creates integration points on the nurbs surface of this geometry.
     * @param return integration points.
     */
    void CreateIntegrationPoints(
        IntegrationPointsArrayType& rIntegrationPoints,
        IntegrationInfo& rIntegrationInfo) const override
    {   
        std::vector<double> spans_u;
        std::vector<double> spans_v;
        std::vector<double> spans_w;
        mpNurbsVolume->SpansLocalSpace(spans_u, 0);
        mpNurbsVolume->SpansLocalSpace(spans_v, 1);
        mpNurbsVolume->SpansLocalSpace(spans_w, 2);

        BrepSBMUtilities::CreateBrepVolumeSBMIntegrationPoints(
            rIntegrationPoints,
            spans_u, spans_v, spans_w,
            *mpSurrogateModelPart_inner, 
            *mpSurrogateModelPart_outer,
            rIntegrationInfo);
        
    }

    ///@}
    ///@name Quadrature Point Geometries
    ///@{

    /* @brief calls function of undelying nurbs surface and updates
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
        KRATOS_WATCH("Start CreateQuadraturePointGeometries")
        mpNurbsVolume->CreateQuadraturePointGeometries(
            rResultGeometries, NumberOfShapeFunctionDerivatives, rIntegrationPoints, rIntegrationInfo);

        for (IndexType i = 0; i < rResultGeometries.size(); ++i) {
            rResultGeometries(i)->SetGeometryParent(this);
        }
        KRATOS_WATCH("End CreateQuadraturePointGeometries")
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

    GeometryData::KratosGeometryFamily GetGeometryFamily() const override
    {
        return GeometryData::KratosGeometryFamily::Kratos_Brep;
    }

    GeometryData::KratosGeometryType GetGeometryType() const override
    {
        return GeometryData::KratosGeometryType::Kratos_Brep_Volume;
    }

    ///@}
    ///@name Information
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "Brep Volume";
    }

    /// Print information about this object.
    void PrintInfo( std::ostream& rOStream ) const override
    {
        rOStream << "Brep Volume";
    }

    /// Print object's data.
    void PrintData( std::ostream& rOStream ) const override
    {
        BaseType::PrintData( rOStream );
        std::cout << std::endl;
        rOStream << "    Brep Volume " << std::endl;
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

    ModelPart* mpSurrogateModelPart_inner = nullptr;
    ModelPart* mpSurrogateModelPart_outer = nullptr;


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
        rSerializer.save("NurbsVolume", mpNurbsVolume);
        rSerializer.save("OuterLoopArray", mOuterLoopArray);
        rSerializer.save("InnerLoopArray", mInnerLoopArray);
        rSerializer.save("EmbeddedEdgesArray", mEmbeddedEdgesArray);
        rSerializer.save("IsTrimmed", mIsTrimmed);
    }

    void load( Serializer& rSerializer ) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BaseType );
        rSerializer.load("NurbsVolume", mpNurbsVolume);
        rSerializer.load("OuterLoopArray", mOuterLoopArray);
        rSerializer.load("InnerLoopArray", mInnerLoopArray);
        rSerializer.load("EmbeddedEdgesArray", mEmbeddedEdgesArray);
        rSerializer.load("IsTrimmed", mIsTrimmed);
    }

    BrepVolume()
        : BaseType( PointsArrayType(), &msGeometryData )
    {}

    ///@}

}; // Class BrepVolume

///@name Input and output
///@{

/// input stream functions
template<class TContainerPointType, class TContainerPointEmbeddedType = TContainerPointType> inline std::istream& operator >> (
    std::istream& rIStream,
    BrepVolume<TContainerPointType, TContainerPointEmbeddedType>& rThis );

/// output stream functions
template<class TContainerPointType, class TContainerPointEmbeddedType = TContainerPointType> inline std::ostream& operator << (
    std::ostream& rOStream,
    const BrepVolume<TContainerPointType, TContainerPointEmbeddedType>& rThis )
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
GeometryData BrepVolume<TContainerPointType, TContainerPointEmbeddedType>::msGeometryData(
    &msGeometryDimension,
    GeometryData::IntegrationMethod::GI_GAUSS_1,
    {}, {}, {});

template<class TContainerPointType, class TContainerPointEmbeddedType>
const GeometryDimension BrepVolume<TContainerPointType, TContainerPointEmbeddedType>::msGeometryDimension(3, 2);

///@}
}// namespace Kratos.

#endif // KRATOS_BREP_FACE_3D_H_INCLUDED  defined

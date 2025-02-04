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

// Project includes
#include "geometries/geometry.h"
#include "geometries/brep_curve_on_surface.h"
#include "geometries/nurbs_shape_function_utilities/nurbs_interval.h"

// trimming integration
#include "utilities/geometry_utilities/brep_trimming_utilities.h"

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
    typedef DenseVector<GeometryPointer> GeometrySurrogateArrayType;

    typedef GeometryData::IntegrationMethod IntegrationMethod;

    typedef NurbsSurfaceGeometry<3, TContainerPointType> NurbsSurfaceType;
    typedef BrepCurveOnSurface<TContainerPointType, TShiftedBoundary, TContainerPointEmbeddedType> BrepCurveOnSurfaceType;

    typedef DenseVector<typename BrepCurveOnSurfaceType::Pointer> BrepCurveOnSurfaceArrayType;
    typedef DenseVector<typename BrepCurveOnSurfaceType::Pointer> BrepCurveOnSurfaceLoopType;
    typedef DenseVector<DenseVector<typename BrepCurveOnSurfaceType::Pointer>> BrepCurveOnSurfaceLoopArrayType;

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

    // Constructor for SBM
    BrepSurface(
        typename NurbsSurfaceType::Pointer pSurface, 
        BrepCurveOnSurfaceLoopArrayType& BrepOuterLoopArray,
        BrepCurveOnSurfaceLoopArrayType& BrepInnerLoopArray,
        GeometrySurrogateArrayType& rSurrogateInnerLoopGeometries, 
        GeometrySurrogateArrayType& rSurrogateOuterLoopGeometries)
        : BaseType(PointsArrayType(), &msGeometryData)
        , mpNurbsSurface(pSurface) 
        , mOuterLoopArray(BrepOuterLoopArray)
        , mInnerLoopArray(BrepInnerLoopArray)
        , mpSurrogateInnerLoopGeometries(rSurrogateInnerLoopGeometries)
        , mpSurrogateOuterLoopGeometries(rSurrogateOuterLoopGeometries)
    {
        mIsTrimmed = false;
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
                // TODO: Next PR -> Call  "BrepSBMUtilities::CreateBrepSurfaceSBMIntegrationPoints"
                mpNurbsSurface->CreateIntegrationPoints(
                    rIntegrationPoints, rIntegrationInfo);
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

            BrepTrimmingUtilities::CreateBrepSurfaceTrimmingIntegrationPoints(
                rIntegrationPoints,
                mOuterLoopArray, mInnerLoopArray,
                spans_u, spans_v,
                rIntegrationInfo);
        }
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
        mpNurbsSurface->CreateQuadraturePointGeometries(
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


    GeometrySurrogateArrayType mpSurrogateInnerLoopGeometries;
    GeometrySurrogateArrayType mpSurrogateOuterLoopGeometries;

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
    }

    void load( Serializer& rSerializer ) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BaseType );
        rSerializer.load("NurbsSurface", mpNurbsSurface);
        rSerializer.load("OuterLoopArray", mOuterLoopArray);
        rSerializer.load("InnerLoopArray", mInnerLoopArray);
        rSerializer.load("EmbeddedEdgesArray", mEmbeddedEdgesArray);
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

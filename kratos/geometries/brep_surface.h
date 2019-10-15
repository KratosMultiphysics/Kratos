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


namespace Kratos
{
///@name Kratos Classes
///@{

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

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor for untrimmed patch
    BrepSurface( 
        typename NurbsSurfaceType::Pointer pSurface)
        : BaseType(PointsArrayType(), &msGeometryData)
        , mpNurbsSurface(pSurface)
    {
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
    {
    }

    /// Destructor
    ~BrepSurface() override = default;

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
    BrepSurface& operator=( const BrepSurface& rOther )
    {
        BaseType::operator=( rOther );
        mpNurbsSurface = rOther.mpNurbsSurface;
        mOuterLoopArray = rOther.mOuterLoopArray;
        mInnerLoopArray = rOther.mInnerLoopArray;
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
    template<class TOtherContainerPointType, class TOtherContainerPointEmbeddedType>
    BrepSurface& operator=( BrepSurface<TOtherContainerPointType, TOtherContainerPointEmbeddedType> const & rOther )
    {
        BaseType::operator=( rOther );
        mpNurbsSurface = rOther.mpNurbsSurface;
        mOuterLoopArray = rOther.mOuterLoopArray;
        mInnerLoopArray = rOther.mInnerLoopArray;
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
    ///@name Information
    ///@{

    /*
    * @brief checks if the Brep Face has any boundary trim information
    * @return true if face has no boundary trimming curves
    */
    bool IsTrimmed() const
    {
        if (mOuterLoopArray.size() == 0 && mInnerLoopArray.size() == 0)
            return false;
    }

    ///@}
    ///@name Geometrical Operations
    ///@{

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
    ///@name Input and output
    ///@{


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
    }

    void load( Serializer& rSerializer ) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BaseType );
        rSerializer.load("NurbsSurface", mpNurbsSurface);
        rSerializer.load("OuterLoopArray", mOuterLoopArray);
        rSerializer.load("InnerLoopArray", mInnerLoopArray);
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

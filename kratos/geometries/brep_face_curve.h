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

#if !defined(KRATOS_BREP_FACE_CURVE_3D_H_INCLUDED )
#define  KRATOS_BREP_FACE_CURVE_3D_H_INCLUDED

// System includes

// External includes

// Project includes
#include "geometries/geometry.h"
#include "geometries/nurbs_shape_function_utilities/nurbs_interval.h"


namespace Kratos
{
///@name Kratos Classes
///@{

/**
 * @class BrepFaceCurve
 * @ingroup KratosCore
 * @brief The BrepFaceCurve acts as topology for faces. Those
 *        can be enclosed by a certain set of boundary curves.
 */
template<class TPointType, class TPointEmbeddedType = TPointType>
class BrepFaceCurve
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
    typedef NurbsCurveGeometry<2, TPointEmbeddedType> NurbsCurveType;

    typedef typename BaseType::GeometriesArrayType GeometriesArrayType;

    typedef typename BaseType::IndexType IndexType;
    typedef typename BaseType::SizeType SizeType;

    typedef typename BaseType::PointsArrayType PointsArrayType;
    typedef typename BaseType::CoordinatesArrayType CoordinatesArrayType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// constructor for untrimmed surface
    BrepFaceCurve( 
        typename NurbsSurfaceType::Pointer pSurface)
        : BaseType(PointsArrayType(), &msGeometryData)
        , mpNurbsSurface(pSurface)
    {
    }

    BrepFaceCurve(
        typename NurbsSurfaceType::Pointer pSurface,
        typename NurbsCurveType::Pointer pCurve,
        Interval CurveInterval)
        : BaseType(PointsArrayType(), &msGeometryData)
        , mpNurbsSurface(pSurface)
        , mpNurbsCurve(pCurve)
        , mCurveInterval(CurveInterval)
    {
    }

    explicit BrepFaceCurve(const PointsArrayType& ThisPoints)
        : BaseType(ThisPoints, &msGeometryData)
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
    template<class TOtherPointType, class TOtherPointEmbeddedType> explicit BrepFaceCurve(
        BrepFaceCurve<TOtherPointType, TOtherPointEmbeddedType> const& rOther )
        : BaseType( rOther )
    {
    }

    /// Destructor
    ~BrepFaceCurve() override = default;

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
    template<class TOtherPointType, class TOtherPointEmbeddedType>
    BrepFaceCurve& operator=( BrepFaceCurve<TOtherPointType, TOtherPointEmbeddedType> const & rOther )
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
        CoordinatesArrayType surface_coordinates(3, 0.0);
        mpNurbsCurve->GlobalCoordinates(surface_coordinates, rLocalCoordinates);
        mpNurbsSurface->GlobalCoordinates(rResult, surface_coordinates);

        return rResult;
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
    ) const
    {
        rResult = rPoint;
        return mCurveInterval.IsInside(rResult[0]);
    }

    ///@name Shape Function
    ///@{

    Vector& ShapeFunctionsValues(
        Vector &rResult,
        const CoordinatesArrayType& rCoordinates) const override
    {
        CoordinatesArrayType surface_coordinates(3, 0.0);
        mpNurbsCurve->GlobalCoordinates(surface_coordinates, rCoordinates);
        mpNurbsSurface->ShapeFunctionsValues(rResult, surface_coordinates);

        return rResult;
    }

    Matrix& ShapeFunctionsLocalGradients(
        Matrix& rResult,
        const CoordinatesArrayType& rCoordinates) const override
    {
        CoordinatesArrayType surface_coordinates(3, 0.0);
        mpNurbsCurve->GlobalCoordinates(surface_coordinates, rCoordinates);
        mpNurbsSurface->ShapeFunctionsLocalGradients(rResult, surface_coordinates);

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
        return "Brep face curve";
    }

    /// Print information about this object.
    void PrintInfo( std::ostream& rOStream ) const override
    {
        rOStream << "Brep face curve";
    }

    /// Print object's data.
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

    static const GeometryData msGeometryData;

    ///@}
    ///@name Member Variables
    ///@{

    typename NurbsSurfaceType::Pointer mpNurbsSurface;
    typename NurbsCurveType::Pointer mpNurbsCurve;

    Interval mCurveInterval;

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save( Serializer& rSerializer ) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, BaseType );
        rSerializer.save("NurbsSurface", mpNurbsSurface);
        rSerializer.save("NurbsCurve", mpNurbsCurve);
    }

    void load( Serializer& rSerializer ) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BaseType );
        rSerializer.load("NurbsSurface", mpNurbsSurface);
        rSerializer.load("NurbsCurve", mpNurbsCurve);
    }

    BrepFaceCurve()
        : BaseType( PointsArrayType(), &msGeometryData )
    {}

    ///@}
    ///@name Private Friends
    ///@{

    template<class TOtherPointType, class TOtherPointEmbeddedType> friend class BrepFaceCurve;

    ///@}
}; // Class BrepFaceCurve

///@name Input and output
///@{

/// input stream functions
template<class TPointType, class TPointEmbeddedType> inline std::istream& operator >> (
    std::istream& rIStream,
    BrepFaceCurve<TPointType, TPointEmbeddedType>& rThis );

/// output stream functions
template<class TPointType, class TPointEmbeddedType> inline std::ostream& operator << (
    std::ostream& rOStream,
    const BrepFaceCurve<TPointType, TPointEmbeddedType>& rThis )
{
    rThis.PrintInfo( rOStream );
    rOStream << std::endl;
    rThis.PrintData( rOStream );
    return rOStream;
}

///@}
///@name Static Type Declarations
///@{

template<class TPointType, class TPointEmbeddedType> const
GeometryData BrepFaceCurve<TPointType, TPointEmbeddedType>::msGeometryData(
    1, 3, 2,
    GeometryData::GI_GAUSS_1,
    {}, {}, {});

///@}
}// namespace Kratos.

#endif // KRATOS_BREP_FACE_CURVE_3D_H_INCLUDED  defined

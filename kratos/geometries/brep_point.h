//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//

#if !defined(KRATOS_BREP_POINT_H_INCLUDED )
#define  KRATOS_BREP_POINT_H_INCLUDED

// System includes

// External includes

// Project includes
#include "geometries/geometry.h"

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
template<class TContainerPointType, int TLocalSpaceDimension, int TWorkingSpaceDimension>
class BrepPoint
    : public Geometry<typename TContainerPointType::value_type>
{
public:
    ///@}
    ///@name Type Definitions
    ///@{

    /// Pointer definition of BrepPoint
    KRATOS_CLASS_POINTER_DEFINITION( BrepPoint );

    typedef typename TContainerPointType::value_type PointType;

    typedef Geometry<typename TContainerPointType::value_type> BaseType;
    typedef Geometry<typename TContainerPointType::value_type> GeometryType;
    typedef typename GeometryType::Pointer GeometryPointer;

    typedef GeometryData::IntegrationMethod IntegrationMethod;

    typedef typename BaseType::IndexType IndexType;
    typedef typename BaseType::SizeType SizeType;

    typedef typename BaseType::PointsArrayType PointsArrayType;
    typedef typename BaseType::CoordinatesArrayType CoordinatesArrayType;
    typedef typename BaseType::GeometriesArrayType GeometriesArrayType;
    typedef typename BaseType::IntegrationPointsArrayType IntegrationPointsArrayType;

    static constexpr IndexType BACKGROUND_GEOMETRY_INDEX = -1;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor for untrimmed patch
    BrepPoint(
        const CoordinatesArrayType LocalCoordinates,
        typename GeometryType::Pointer pBackgroundGeometry)
        : BaseType(PointsArrayType(), &msGeometryData)
        , mLocalCoordinates(LocalCoordinates)
        , mpBackgroundGeometry(pBackgroundGeometry)
    {
        KRATOS_ERROR_IF(pBackgroundGeometry->LocalSpaceDimension() != this->LocalSpaceDimension())
            << "Local space dimension of background geometry (" << pBackgroundGeometry->LocalSpaceDimension()
            << ") does not coincide with LocalSpaceDimension of this geometry (" << this->LocalSpaceDimension() << ")."
            << std::endl;
    }

    explicit BrepPoint(const PointsArrayType& ThisPoints)
        : BaseType(ThisPoints, &msGeometryData)
    {
    }

    /// Copy constructor.
    BrepPoint(BrepPoint const& rOther )
        : BaseType( rOther )
        , mLocalCoordinates(rOther.mLocalCoordinates)
        , mpBackgroundGeometry(rOther.mpBackgroundGeometry)
    {
    }

    /// Destructor
    ~BrepPoint() override = default;

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    BrepPoint& operator=( const BrepPoint& rOther )
    {
        BaseType::operator=( rOther );
        mLocalCoordinates = rOther.mLocalCoordinates;
        mpBackgroundGeometry = rOther.mpBackgroundGeometry;
        return *this;
    }

    ///@}
    ///@name Operations
    ///@{

    typename BaseType::Pointer Create( PointsArrayType const& ThisPoints ) const override
    {
        return typename BaseType::Pointer( new BrepPoint( ThisPoints ) );
    }

    ///@}
    ///@name Access to Geometry Parts
    ///@{

    /**
    * @brief This function returns the pointer of the background geometry
    *        of this point. Accesseable via BACKGROUND_GEOMETRY_INDEX.
    * @param Index: only possible with BACKGROUND_GEOMETRY_INDEX.
    * @return pointer of geometry, corresponding to the index.
    */
    GeometryPointer pGetGeometryPart(IndexType Index) override
    {
        const auto& const_this = *this;
        return std::const_pointer_cast<GeometryType>(
            const_this.pGetGeometryPart(Index));
    }

    /**
    * @brief This function returns the pointer of the background geometry
    *        of this point. Accesseable via BACKGROUND_GEOMETRY_INDEX.
    * @param Index: only possible with BACKGROUND_GEOMETRY_INDEX.
    * @return pointer of geometry, corresponding to the index.
    */
    const GeometryPointer pGetGeometryPart(IndexType Index) const override
    {
        if (Index == BACKGROUND_GEOMETRY_INDEX)
            return mpBackgroundGeometry;

        KRATOS_ERROR << "Index " << Index << " not existing as geometry part in BrepPoint #"
            << this->Id() << std::endl;
    }

    /* @brief only true if BACKGROUND_GEOMETRY_INDEX is used.
     *        Brep Points have no other geometries.
     */
    bool HasGeometryPart(IndexType Index) const override
    {
        if (Index == BACKGROUND_GEOMETRY_INDEX)
            return true;

        return false;
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
        mpBackgroundGeometry->GlobalCoordinates(rResult, mLocalCoordinates);

        return rResult;
    }

    ///@}
    ///@name Integration Points
    ///@{

    /* Creates a integration point on the background geometry.
     * @param return integration points list with one entry.
     */
    void CreateIntegrationPoints(
        IntegrationPointsArrayType& rIntegrationPoints) const override
    {
        rIntegrationPoints.resize(1);
        rIntegrationPoints[0][0] = mLocalCoordinates[0];
        rIntegrationPoints[0][1] = mLocalCoordinates[1];
        rIntegrationPoints[0][2] = mLocalCoordinates[2];
        rIntegrationPoints[0].Weight() = 1.0;
    }

    ///@}
    ///@name Quadrature Point Geometries
    ///@{

    /* @brief creates a list of quadrature point geometry
     *        at the specific location of this BrepPoint on the background geometry.
     *
     * @param rResultGeometries list of ONE quadrature point geometries.
     * @param rIntegrationPoints list of ONE integration point.
     * @param NumberOfShapeFunctionDerivatives the number of evaluated
     *        derivatives of shape functions at the quadrature point geometries.
     *
     * @see quadrature_point_geometry.h
     */
    void CreateQuadraturePointGeometries(
        GeometriesArrayType& rResultGeometries,
        IndexType NumberOfShapeFunctionDerivatives) override
    {
        IntegrationPointsArrayType integration_point(1);
        this->CreateIntegrationPoints(integration_point);

        mpBackgroundGeometry->CreateQuadraturePointGeometries(
            rResultGeometries, NumberOfShapeFunctionDerivatives, integration_point);
    }

    ///@}
    ///@name Shape Function
    ///@{

    Vector& ShapeFunctionsValues(
        Vector &rResult,
        const CoordinatesArrayType& rCoordinates) const override
    {
        mpBackgroundGeometry->ShapeFunctionsValues(rResult, mLocalCoordinates);

        return rResult;
    }

    Matrix& ShapeFunctionsLocalGradients(
        Matrix& rResult,
        const CoordinatesArrayType& rCoordinates) const override
    {
        mpBackgroundGeometry->ShapeFunctionsLocalGradients(rResult, mLocalCoordinates);

        return rResult;
    }

    ///@}
    ///@name Information
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "BrepPoint";
    }

    /// Print information about this object.
    void PrintInfo( std::ostream& rOStream ) const override
    {
        rOStream << "BrepPoint";
    }

    /// Print object's data.
    void PrintData( std::ostream& rOStream ) const override
    {
        BaseType::PrintData( rOStream );
        std::cout << std::endl;
        rOStream << "    BrepPoint" << std::endl;
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

    CoordinatesArrayType mLocalCoordinates;

    typename GeometryType::Pointer mpBackgroundGeometry;

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save( Serializer& rSerializer ) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, BaseType );
        rSerializer.save("LocalCoordinates", mLocalCoordinates);
        rSerializer.save("pBackgroundGeometry", mpBackgroundGeometry);
    }

    void load( Serializer& rSerializer ) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BaseType );
        rSerializer.load("LocalCoordinates", mLocalCoordinates);
        rSerializer.load("pBackgroundGeometry", mpBackgroundGeometry);
    }

    BrepPoint()
        : BaseType( PointsArrayType(), &msGeometryData )
    {}

    ///@}

}; // Class BrepPoint

///@name Input and output
///@{

/// input stream functions
template<class TContainerPointType, int TLocalSpaceDimension, int TWorkingSpaceDimension> inline std::istream& operator >> (
    std::istream& rIStream,
    BrepPoint<TContainerPointType, TLocalSpaceDimension, TWorkingSpaceDimension>& rThis );

/// output stream functions
template<class TContainerPointType, int TLocalSpaceDimension, int TWorkingSpaceDimension> inline std::ostream& operator << (
    std::ostream& rOStream,
    const BrepPoint<TContainerPointType, TLocalSpaceDimension, TWorkingSpaceDimension>& rThis )
{
    rThis.PrintInfo( rOStream );
    rOStream << std::endl;
    rThis.PrintData( rOStream );
    return rOStream;
}

///@}
///@name Static Type Declarations
///@{

template<class TContainerPointType, int TLocalSpaceDimension, int TWorkingSpaceDimension> const
GeometryData BrepPoint<TContainerPointType, TLocalSpaceDimension, TWorkingSpaceDimension>::msGeometryData(
    &msGeometryDimension,
    GeometryData::GI_GAUSS_1,
    {}, {}, {});

template<class TContainerPointType, int TLocalSpaceDimension, int TWorkingSpaceDimension>
const GeometryDimension BrepPoint<TContainerPointType, TLocalSpaceDimension, TWorkingSpaceDimension>::msGeometryDimension(
    TLocalSpaceDimension, TWorkingSpaceDimension, 0);

///@}
}// namespace Kratos.

#endif // KRATOS_BREP_FACE_3D_H_INCLUDED  defined

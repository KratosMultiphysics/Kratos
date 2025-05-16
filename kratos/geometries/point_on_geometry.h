//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//

#if !defined(KRATOS_POINT_ON_GEOMETRY_H_INCLUDED )
#define  KRATOS_POINT_ON_GEOMETRY_H_INCLUDED

// System includes

// External includes

// Project includes
#include "geometries/geometry.h"

#include "quadrature_point_geometry.h"

namespace Kratos
{
///@name Kratos Classes
///@{

/* @class PointOnGeometry
 * @ingroup KratosCore
 * @brief The PointOnGeometry acts as topology for points on various types
 *        of geometries as faces or curves.
 */
template<class TContainerPointType, int TWorkingSpaceDimension, int TLocalSpaceDimensionOfBackground>
class PointOnGeometry
    : public Geometry<typename TContainerPointType::value_type>
{
public:
    ///@}
    ///@name Type Definitions
    ///@{

    /// Pointer definition of PointOnGeometry
    KRATOS_CLASS_POINTER_DEFINITION( PointOnGeometry );

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

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor with the local coordinates and background geometry
    PointOnGeometry(
        const CoordinatesArrayType LocalCoordinates,
        typename GeometryType::Pointer pBackgroundGeometry)
        : BaseType(PointsArrayType(), &msGeometryData)
        , mLocalCoordinates(LocalCoordinates)
        , mpBackgroundGeometry(pBackgroundGeometry)
    {
        KRATOS_ERROR_IF(pBackgroundGeometry->WorkingSpaceDimension() != this->WorkingSpaceDimension())
            << "Working space dimension of background geometry (" << pBackgroundGeometry->WorkingSpaceDimension()
            << ") does not coincide with WorkingSpaceDimension of this geometry (" << this->WorkingSpaceDimension() << ")."
            << std::endl;
        KRATOS_ERROR_IF(pBackgroundGeometry->LocalSpaceDimension() != this->LocalSpaceDimension())
            << "Local space dimension of background geometry (" << pBackgroundGeometry->LocalSpaceDimension()
            << ") does not coincide with LocalSpaceDimension of this geometry (" << this->LocalSpaceDimension() << ")."
            << std::endl;
    }

    /// Copy constructor.
    PointOnGeometry(PointOnGeometry const& rOther )
        : BaseType( rOther )
        , mLocalCoordinates(rOther.mLocalCoordinates)
        , mpBackgroundGeometry(rOther.mpBackgroundGeometry)
    {
    }

    /// Destructor
    ~PointOnGeometry() override = default;

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    PointOnGeometry& operator=( const PointOnGeometry& rOther )
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
        return typename BaseType::Pointer( new PointOnGeometry( ThisPoints ) );
    }

    ///@}
    ///@name Access to Geometry Parts
    ///@{

    /**
    * @brief This function returns the pointer of the background geometry
    *        of this point. Accesseable via GeometryType::BACKGROUND_GEOMETRY_INDEX.
    * @param Index: only possible with GeometryType::BACKGROUND_GEOMETRY_INDEX.
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
    *        of this point. Accesseable via GeometryType::BACKGROUND_GEOMETRY_INDEX.
    * @param Index: only possible with GeometryType::BACKGROUND_GEOMETRY_INDEX.
    * @return pointer of geometry, corresponding to the index.
    */
    const GeometryPointer pGetGeometryPart(IndexType Index) const override
    {
        if (Index == GeometryType::BACKGROUND_GEOMETRY_INDEX)
            return mpBackgroundGeometry;

        KRATOS_ERROR << "Index " << Index << " not existing as geometry part in PointOnGeometry #"
            << this->Id() << std::endl;
    }

    /* @brief only true if BACKGROUND_GEOMETRY_INDEX is used.
     *        Brep Points have no other geometries.
     */
    bool HasGeometryPart(IndexType Index) const override
    {
        return (Index == GeometryType::BACKGROUND_GEOMETRY_INDEX);
    }

    ///@}
    ///@name Mathematical Informations
    ///@{

    /// Return polynomial degree of the background geometry
    SizeType PolynomialDegree(IndexType LocalDirectionIndex) const override
    {
        return mpBackgroundGeometry->PolynomialDegree(LocalDirectionIndex);
    }

    ///@}
    ///@name Geometrical Operations
    ///@{

    /* @brief Calculates the global location of this point on the background geometry.
     *
     * @return Point which is the location of this quadrature point.
     */
    Point Center() const override
    {
        CoordinatesArrayType global_location;
        mpBackgroundGeometry->GlobalCoordinates(global_location, mLocalCoordinates);
        return Point(global_location);
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
        KRATOS_ERROR << "Trying to access the GlobalCoordinates via PointOnGeometry. "
            << "However, this is not possible. Try to access the global coordinates of"
            << " any local coordinates via accessing the background geometry through: "
            << "pGetGeometryPart(GeometryType::BACKGROUND_GEOMETRY_INDEX)." << std::endl;

        return rResult;
    }

    ///@}
    ///@name Integration Info
    ///@{

    /// Provides the default integration dependent on the polynomial degree.
    IntegrationInfo GetDefaultIntegrationInfo() const override
    {
        return IntegrationInfo(0, 0, IntegrationInfo::QuadratureMethod::Default);
    }

    ///@}
    ///@name Integration Points
    ///@{

    /* Creates a integration point on the background geometry.
     * @param return integration points list with one entry.
     */
    void CreateIntegrationPoints(
        IntegrationPointsArrayType& rIntegrationPoints,
        IntegrationInfo& rIntegrationInfo) const override
    {
        if (rIntegrationPoints.size() != 1) { rIntegrationPoints.resize(1); }
        rIntegrationPoints[0][0] = mLocalCoordinates[0];
        rIntegrationPoints[0][1] = mLocalCoordinates[1];
        rIntegrationPoints[0][2] = mLocalCoordinates[2];
        rIntegrationPoints[0].Weight() = 1.0;
    }

    ///@}
    ///@name Quadrature Point Geometries
    ///@{

    /* @brief creates a list of quadrature point geometry
     *        at the specific location of this PointOnGeometry on the background geometry.
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
        IndexType NumberOfShapeFunctionDerivatives,
        IntegrationInfo& rIntegrationInfo) override
    {
        IntegrationPointsArrayType integration_point(1);
        this->CreateIntegrationPoints(integration_point, rIntegrationInfo);

        GeometriesArrayType rQuadraturePointGeometries(1);
        mpBackgroundGeometry->CreateQuadraturePointGeometries(
            rQuadraturePointGeometries, NumberOfShapeFunctionDerivatives, integration_point, rIntegrationInfo);

        if (rResultGeometries.size() != 1) { rResultGeometries.resize(1); }
        // assignment operator for quadrature point geometry with Dimension being 0.
        rResultGeometries(0) = Kratos::make_shared<
            QuadraturePointGeometry<PointType, TWorkingSpaceDimension, TLocalSpaceDimensionOfBackground, 0>>(
                std::move(rQuadraturePointGeometries(0)->Points()),
                rQuadraturePointGeometries(0)->GetGeometryData().GetGeometryShapeFunctionContainer(),
                this);
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
    ///@name Geometry Family
    ///@{

    /**
     * @brief Gets the geometry family.
     * @details This function returns the family type of the geometry. The geometry family categorizes the geometry into a broader classification, aiding in its identification and processing.
     * @return GeometryData::KratosGeometryFamily The geometry family.
     */
    GeometryData::KratosGeometryFamily GetGeometryFamily() const override
    {
        return GeometryData::KratosGeometryFamily::Kratos_Point;
    }

    ///@}
    ///@name Information
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "PointOnGeometry";
    }

    /// Print information about this object.
    void PrintInfo( std::ostream& rOStream ) const override
    {
        rOStream << "PointOnGeometry";
    }

    /// Print object's data.
    void PrintData( std::ostream& rOStream ) const override
    {
        BaseType::PrintData( rOStream );
        std::cout << std::endl;
        rOStream << "    PointOnGeometry" << std::endl;
    }

    ///@}

private:
    static const GeometryData msGeometryData;

    static const GeometryDimension msGeometryDimension;

    ///@name Member Variables
    ///@{

    CoordinatesArrayType mLocalCoordinates;

    typename GeometryType::Pointer mpBackgroundGeometry;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor for serializer.
    explicit PointOnGeometry(const PointsArrayType& ThisPoints)
        : BaseType(ThisPoints, &msGeometryData)
    {
    }

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

    PointOnGeometry()
        : BaseType( PointsArrayType(), &msGeometryData )
    {}

    ///@}

}; // Class PointOnGeometry

///@name Input and output
///@{

/// input stream functions
template<class TContainerPointType, int TWorkingSpaceDimension, int TLocalSpaceDimensionOfBackground> inline std::istream& operator >> (
    std::istream& rIStream,
    PointOnGeometry<TContainerPointType, TWorkingSpaceDimension, TLocalSpaceDimensionOfBackground>& rThis );

/// output stream functions
template<class TContainerPointType, int TWorkingSpaceDimension, int TLocalSpaceDimensionOfBackground> inline std::ostream& operator << (
    std::ostream& rOStream,
    const PointOnGeometry<TContainerPointType, TWorkingSpaceDimension, TLocalSpaceDimensionOfBackground>& rThis )
{
    rThis.PrintInfo( rOStream );
    rOStream << std::endl;
    rThis.PrintData( rOStream );
    return rOStream;
}

///@}
///@name Static Type Declarations
///@{

template<class TContainerPointType, int TWorkingSpaceDimension, int TLocalSpaceDimensionOfBackground> const
GeometryData PointOnGeometry<TContainerPointType, TWorkingSpaceDimension, TLocalSpaceDimensionOfBackground>::msGeometryData(
    &msGeometryDimension,
    GeometryData::IntegrationMethod::GI_GAUSS_1,
    {}, {}, {});

template<class TContainerPointType, int TWorkingSpaceDimension, int TLocalSpaceDimensionOfBackground>
const GeometryDimension PointOnGeometry<TContainerPointType, TWorkingSpaceDimension, TLocalSpaceDimensionOfBackground>::msGeometryDimension(TWorkingSpaceDimension, TLocalSpaceDimensionOfBackground);

///@}
}// namespace Kratos.

#endif // KRATOS_POINT_ON_GEOMETRY_H_INCLUDED  defined

//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Andrea Gorgi
//                   
//

#if !defined(KRATOS_NURBS_COUPLING_GEOMETRY_2d_H_INCLUDED )
#define  KRATOS_NURBS_COUPLING_GEOMETRY_2d_H_INCLUDED

// System includes

// External includes

// Project includes
#include "geometries/geometry.h"
#include "geometries/nurbs_shape_function_utilities/nurbs_interval.h"

#include "geometries/nurbs_curve_on_surface_geometry.h"
#include "geometries/brep_curve_on_surface.h"

#include "utilities/nurbs_utilities/projection_nurbs_geometry_utilities.h"

#include "includes/io.h"

namespace Kratos
{
///@name Kratos Classes
///@{

/**
 * @class NurbsCouplingGeometry2D
 * @ingroup KratosCore
 * @brief The NurbsCouplingGeometry2D acts as topology for nurbs coupled curves on surfaces.
 */
template<class TPointType, class TSurfaceContainerPointType> class NurbsCouplingGeometry2D
    : public Geometry<TPointType>
{
public:
    //@name Type Definitions
    ///@{

    /// Geometry as base class.
    typedef Geometry<TPointType> BaseType;
    typedef Geometry<TPointType> GeometryType;

    typedef typename GeometryType::Pointer GeometryPointer;
    typedef std::vector<GeometryPointer> GeometryPointerVector;

    /// Pointer definition of CouplingGeometry
    KRATOS_CLASS_POINTER_DEFINITION( NurbsCouplingGeometry2D );

    typedef TPointType PointType;

    typedef typename BaseType::IndexType IndexType;
    typedef typename BaseType::SizeType SizeType;

    typedef typename BaseType::IntegrationPointsArrayType IntegrationPointsArrayType;
    typedef typename BaseType::PointsArrayType PointsArrayType;
    typedef typename BaseType::CoordinatesArrayType CoordinatesArrayType;
    typedef std::vector<CoordinatesArrayType> CoordinatesArrayVectorType;
    typedef PointerVector<GeometryType> GeometriesArrayType;

    typedef NurbsSurfaceGeometry<3, TSurfaceContainerPointType> NurbsSurfaceType;
    typedef typename TSurfaceContainerPointType::value_type NodeType;

    typedef typename NurbsSurfaceType::Pointer NurbsSurfaceTypePointer;

    // MODIFIED --------------------------------------
    typedef PointerVector<Node> ContainerNodeType;
    typedef PointerVector<Point> ContainerEmbeddedNodeType;

    typedef BrepCurveOnSurface<ContainerNodeType, ContainerEmbeddedNodeType> BrepCurveOnSurfaceType;

    typedef DenseVector<typename BrepCurveOnSurfaceType::Pointer> BrepCurveOnSurfaceArrayType;
    //--------------------------------------------------------------------------------------

    static constexpr IndexType MasterIndex = 0;
    static constexpr IndexType SlaveIndex = 1;

    static constexpr IndexType CURVE_ON_SURFACE_INDEX = std::numeric_limits<IndexType>::max() - 2;
    static constexpr IndexType EMBEDDED_CURVE_INDEX = std::numeric_limits<IndexType>::max() - 3;


    ///@}
    ///@name Life Cycle
    ///@{

    /// constructor for untrimmed surface
    /// constructor for untrimmed surface
NurbsCouplingGeometry2D(
    GeometriesArrayType& slaveGeometryList,
    GeometriesArrayType& masterGeometryList)
{
    // Assicurati di ottenere il puntatore alle superfici corrette.
    mpNurbsSurfaceMaster = std::dynamic_pointer_cast<NurbsSurfaceType>(
        masterGeometryList[0].pGetGeometryPart(GeometryType::BACKGROUND_GEOMETRY_INDEX));
    mpNurbsSurfaceSlave = std::dynamic_pointer_cast<NurbsSurfaceType>(
        slaveGeometryList[0].pGetGeometryPart(GeometryType::BACKGROUND_GEOMETRY_INDEX));

    // Ridimensiona le liste
    mSlaveGeometryList.resize(slaveGeometryList.size());
    mMasterGeometryList.resize(masterGeometryList.size());

    // Itera su ogni elemento per fare il cast da Geometry a BrepCurveOnSurface
    for (std::size_t i = 0; i < slaveGeometryList.size(); ++i) {
        auto slave_part = std::dynamic_pointer_cast<BrepCurveOnSurfaceType>(slaveGeometryList(i));
        if (!slave_part) {
            KRATOS_ERROR << "Element " << i << " in slaveGeometryList is not of type BrepCurveOnSurface." << std::endl;
        }
        mSlaveGeometryList[i] = slave_part;
    }

    for (std::size_t i = 0; i < masterGeometryList.size(); ++i) {
        auto master_part = std::dynamic_pointer_cast<BrepCurveOnSurfaceType>(masterGeometryList(i));
        if (!master_part) {
            KRATOS_ERROR << "Element " << i << " in masterGeometryList is not of type BrepCurveOnSurface." << std::endl;
        }
        mMasterGeometryList[i] = master_part;
    }
}


    explicit NurbsCouplingGeometry2D()
        : BaseType()
    {
    }

    /// Copy constructor.
    NurbsCouplingGeometry2D( NurbsCouplingGeometry2D const& rOther)
        : mSlaveGeometryList(rOther.mSlaveGeometryList)
        , mMasterGeometryList(rOther.mMasterGeometryList)
    {
    }

    /// Destructor
    ~NurbsCouplingGeometry2D() override = default;

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    NurbsCouplingGeometry2D& operator=( const NurbsCouplingGeometry2D& rOther )
    {
        mSlaveGeometryList = rOther.mSlaveGeometryList;
        mMasterGeometryList = rOther.mMasterGeometryList;
        return *this;
    }

    typename BaseType::Pointer Create(
        PointsArrayType const& ThisPoints ) const override
    {
        return Kratos::make_shared<NurbsCouplingGeometry2D>();
    }


    ///@}
    ///@name Operations
    ///@{


    ///@}
    ///@name Access to Geometry Parts
    ///@{

    // /**
    // * @brief This function returns the pointer of the geometry
    // *        which is corresponding to the index.
    // *        Possible indices are:
    // *        SURFACE_INDEX, EMBEDDED_CURVE_INDEX or CURVE_ON_SURFACE_INDEX.
    // * @param Index: SURFACE_INDEX, EMBEDDED_CURVE_INDEX or CURVE_ON_SURFACE_INDEX.
    // * @return pointer of geometry, corresponding to the index.
    // */
    // GeometryPointer pGetGeometryPart(const IndexType Index) override
    // {
    //     const auto& const_this = *this;
    //     return std::const_pointer_cast<GeometryType>(
    //         const_this.pGetGeometryPart(Index));
    // }

    // /**
    // * @brief This function returns the pointer of the geometry
    // *        which is corresponding to the index.
    // *        Possible indices are:
    // *        GeometryType::BACKGROUND_GEOMETRY_INDEX, EMBEDDED_CURVE_INDEX or CURVE_ON_SURFACE_INDEX.
    // * @param Index: SURFACE_INDEX, EMBEDDED_CURVE_INDEX or CURVE_ON_SURFACE_INDEX.
    // * @return pointer of geometry, corresponding to the index.
    // */
    // const GeometryPointer pGetGeometryPart(const IndexType Index) const override
    // {
    //     if (Index == GeometryType::BACKGROUND_GEOMETRY_INDEX)
    //         return mpCurveOnSurface->pGetGeometryPart(GeometryType::BACKGROUND_GEOMETRY_INDEX);

    //     if (Index == CURVE_ON_SURFACE_INDEX)
    //         return mpCurveOnSurface;

    //     KRATOS_ERROR << "Index " << Index << " not existing in BrepCurveOnSurface: "
    //         << this->Id() << std::endl;
    // }

    // /**
    // * @brief This function is used to check if the index is either
    // *        GeometryType::BACKGROUND_GEOMETRY_INDEX or CURVE_ON_SURFACE_INDEX.
    // * @param Index of the geometry part.
    // * @return true if GeometryType::BACKGROUND_GEOMETRY_INDEX or CURVE_ON_SURFACE_INDEX.
    // */
    // bool HasGeometryPart(const IndexType Index) const override
    // {
    //     return (Index == GeometryType::BACKGROUND_GEOMETRY_INDEX || Index == CURVE_ON_SURFACE_INDEX);
    // }

    ///@}
    ///@name Set / Calculate access
    ///@{


    ///@}
    ///@name Mathematical Informations
    ///@{


    ///@}
    ///@name Set/ Get functions
    ///@{


    /// Returns the Slave Geometry list of this coupling geometry.
    BrepCurveOnSurfaceArrayType pGetSlaveGeometry()
    {
        return mSlaveGeometryList;
    }

    /// Returns the Master Geometry list of this coupling geometry.
    BrepCurveOnSurfaceArrayType pGetMasterGeometry()
    {
        return mMasterGeometryList;
    }


    ///@}
    ///@name Curve Properties
    ///@{

    /* @brief Provides intersections of the nurbs curve with the knots of the surface,
     *         using the interval of this curve.
     * @param vector of span intervals.
     * @param index of chosen direction, for curves always 0.
     */
    void SpansLocalSpaceForParentIntegration(NurbsSurfaceTypePointer rpNurbsSurfaceParent,
                                             NurbsSurfaceTypePointer rpNurbsSurfacePaired, 
                                             std::vector<std::vector<double>>& integration_edges_on_parameter_parent_list,
                                             std::vector<std::vector<double>>& spans_parent_list,
                                             std::vector<std::vector<double>>& spans_paired_list,
                                             BrepCurveOnSurfaceArrayType& rParentGeometryList,
                                             BrepCurveOnSurfaceArrayType& rPairedGeometryList);

    ///@}

    ///@}
    ///@name Integration Info
    ///@{

    /// Provides the default integration dependent on the polynomial degree.
    IntegrationInfo GetDefaultIntegrationInfo() const override
    {
        // IndexType numberOfIntegrationPoint_slave = (mSlaveGeometryList[0].GetDefaultIntegrationInfo()).mNumberOfIntegrationPointsPerSpanVector[0];
        // IndexType numberOfIntegrationPoint_master = (mMasterGeometryList[0].GetDefaultIntegrationInfo()).mNumberOfIntegrationPointsPerSpanVector[0];
        // return IntegrationInfo(1, std::max(numberOfIntegrationPoint_master, numberOfIntegrationPoint_slave), IntegrationInfo::QuadratureMethod::GAUSS);
        const SizeType local_space_dimension = this->LocalSpaceDimension();

        std::vector<SizeType> number_of_points_per_span_per_direction(local_space_dimension);
        std::vector<IntegrationInfo::QuadratureMethod> quadrature_method_per_direction(local_space_dimension);
        for (IndexType i = 0; i < local_space_dimension; ++i) {
            SizeType max_p = 0;

            max_p = std::max(mSlaveGeometryList[0]->PolynomialDegree(i) + 1, max_p);
            max_p = std::max(mMasterGeometryList[0]->PolynomialDegree(i) + 1, max_p);

            number_of_points_per_span_per_direction[i] = max_p;
            quadrature_method_per_direction[i] = IntegrationInfo::QuadratureMethod::GAUSS;
        }

        return IntegrationInfo(number_of_points_per_span_per_direction, quadrature_method_per_direction);
    }

    ///@}
    ///@name Integration Points
    ///@{

    /* Creates integration points on the nurbs surface of this geometry.
     * @param return integration points.
     */
    void CreateIntegrationPoints(
        NurbsSurfaceTypePointer rpNurbsSurfaceParent,
        NurbsSurfaceTypePointer rpNurbsSurfacePaired,
        std::vector<IntegrationPointsArrayType>& rIntegrationPoints,
        IntegrationInfo& rIntegrationInfo,
        std::vector<std::vector<double>>& spans_parent_list,
        std::vector<std::vector<double>>& spans_paired_list,
        BrepCurveOnSurfaceArrayType& rParentGeometryList,
        BrepCurveOnSurfaceArrayType& rPairedGeometryList
        );
    
    ///@}
    ///@name Integration Points
    ///@{

    ///@}
    ///@name Create Quadrature Points Geometries
    ///@{

    /* Creates integration points on the nurbs surface of this geometry.
     * @param return integration points.
     */

    void CreateQuadraturePointGeometries(
        GeometriesArrayType& rResultGeometries,
        IndexType NumberOfShapeFunctionDerivatives,
        IntegrationInfo& rIntegrationInfo) override;


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
    void CreateQuadraturePointGeometriesOnParent(
        GeometriesArrayType& rResultGeometries,
        IndexType NumberOfShapeFunctionDerivatives,
        const std::vector<IntegrationPointsArrayType>& rIntegrationPoints,
        IntegrationInfo& rIntegrationInfo,
        NurbsSurfaceTypePointer rpNurbsSurfaceParent,
        NurbsSurfaceTypePointer rpNurbsSurfacePaired,
        BrepCurveOnSurfaceArrayType& rParentGeometryList,
        BrepCurveOnSurfaceArrayType& rPairedGeometryList,
        const IndexType& integrationDomain,
        SizeType& quadraturePointId);


    void CreateQuadraturePointGeometriesSlave(
        GeometriesArrayType& rResultGeometries,
        IndexType NumberOfShapeFunctionDerivatives,
        const IntegrationPointsArrayType& rIntegrationPoints,
        IntegrationInfo& rIntegrationInfo);

    ///@}




    ///@}
    
    GeometryData::KratosGeometryType GetGeometryType() const override
    {
        return GeometryData::KratosGeometryType::Kratos_Nurbs_Coupling_Geometry_2d;
    }


    ///@}
    ///@name Information
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "NURBS COUPLING GEOMETRIES 2D";
    }

    /// Print information about this object.
    void PrintInfo( std::ostream& rOStream ) const override
    {
        rOStream << "NURBS COUPLING GEOMETRIES 2D";
    }

    /// Print object's data.
    void PrintData( std::ostream& rOStream ) const override
    {
        BaseType::PrintData( rOStream );
        std::cout << std::endl;
        rOStream << "    NURBS COUPLING GEOMETRIES 2D : " << std::endl;
    }

    ///@}

private:
    ///@name Static Member Variables
    ///@{

    static const GeometryData msGeometryData;

    static const GeometryDimension msGeometryDimension;

    BrepCurveOnSurfaceArrayType mMasterGeometryList;
    BrepCurveOnSurfaceArrayType mSlaveGeometryList;

    NurbsSurfaceTypePointer mpNurbsSurfaceMaster;
    NurbsSurfaceTypePointer mpNurbsSurfaceSlave;
    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save( Serializer& rSerializer ) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, BaseType );
        rSerializer.save("SlaveGeometryList", mSlaveGeometryList);
        rSerializer.save("MasterGeometryList", mMasterGeometryList);
    }

    void load( Serializer& rSerializer ) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BaseType );
        rSerializer.load("SlaveGeometryList", mSlaveGeometryList);
        rSerializer.load("MasterGeometryList", mMasterGeometryList);
    }

    ///@}
}; // Class NurbsCouplingGeometry2D

template<class TPointType, class TSurfaceContainerPointType> inline std::istream& operator >> (
    std::istream& rIStream,
    NurbsCouplingGeometry2D<TPointType, TSurfaceContainerPointType>& rThis );
/**
 * output stream functions
 */
template<class TPointType, class TSurfaceContainerPointType> inline std::ostream& operator << (
    std::ostream& rOStream,
    const NurbsCouplingGeometry2D<TPointType, TSurfaceContainerPointType>& rThis )
{
    rThis.PrintInfo( rOStream );
    rOStream << std::endl;
    rThis.PrintData( rOStream );
    return rOStream;
}

///@}
}// namespace Kratos.

#endif // KRATOS_BREP_CURVE_ON_SURFACE_3D_H_INCLUDED  defined

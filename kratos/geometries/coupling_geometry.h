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
//                   Pooyan Dadvand
//                   Philipp Bucher
//

#if !defined(KRATOS_COUPLING_GEOMETRY_H_INCLUDED )
#define  KRATOS_COUPLING_GEOMETRY_H_INCLUDED

// System includes

// External includes

// Project includes
#include "geometry.h"

#include "integration/integration_point_utilities.h"
#include "utilities/tessellation_utilities/curve_tessellation.h"

namespace Kratos
{
///@name Kratos Classes
///@{

/**
 * @class CouplingGeometry
 * @ingroup KratosCore
 * @brief The CouplingGeometry connects two or more geometries
 *        of different types and entities.
 */
template<class TPointType> class CouplingGeometry
    : public Geometry<TPointType>
{
public:
    ///@name Type Definitions
    ///@{

    /// Geometry as base class.
    typedef Geometry<TPointType> BaseType;
    typedef Geometry<TPointType> GeometryType;

    typedef typename GeometryType::Pointer GeometryPointer;
    typedef std::vector<GeometryPointer> GeometryPointerVector;

    /// Pointer definition of CouplingGeometry
    KRATOS_CLASS_POINTER_DEFINITION( CouplingGeometry );

    typedef TPointType PointType;

    typedef typename BaseType::IndexType IndexType;
    typedef typename BaseType::SizeType SizeType;

    typedef typename BaseType::IntegrationPointsArrayType IntegrationPointsArrayType;
    typedef typename BaseType::PointsArrayType PointsArrayType;
    typedef typename BaseType::CoordinatesArrayType CoordinatesArrayType;
    typedef std::vector<CoordinatesArrayType> CoordinatesArrayVectorType;
    typedef PointerVector<GeometryType> GeometriesArrayType;

    ///@}
    ///@name Public Static Members
    ///@{

    static constexpr IndexType Master = 0;
    static constexpr IndexType Slave = 1;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor for coupling one master to one slave geometry.
    CouplingGeometry(
        GeometryPointer pMasterGeometry,
        GeometryPointer pSlaveGeometry)
        : BaseType(PointsArrayType(), &(pMasterGeometry->GetGeometryData()))
    {
        mpGeometries.resize(2);

        mpGeometries[0] = pMasterGeometry;
        mpGeometries[1] = pSlaveGeometry;
    }

    /// Constructor for coupling multiple points.
    CouplingGeometry(
        GeometryPointerVector GeometryPointerVector)
        : BaseType(PointsArrayType(), &(GeometryPointerVector[0]->GetGeometryData()))
        , mpGeometries(GeometryPointerVector)
    {
    }

    explicit CouplingGeometry()
        : BaseType()
    {
    }

    /// Copy constructor.
    CouplingGeometry( CouplingGeometry const& rOther )
        : BaseType( rOther )
        , mpGeometries(rOther.mpGeometries)
    {
    }

    /// Copy constructor with other point type.
    template<class TOtherPointType> explicit CouplingGeometry(
        CouplingGeometry<TOtherPointType> const& rOther )
        : BaseType( rOther )
        , mpGeometries( rOther.mpGeometries )
    {
    }

    /// Destructor
    ~CouplingGeometry() override = default;

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    CouplingGeometry& operator=( const CouplingGeometry& rOther )
    {
        BaseType::operator=( rOther );
        mpGeometries = rOther.mpGeometries;
        return *this;
    }

    /// @brief Assignment operator with different point type.
    template<class TOtherPointType>
    CouplingGeometry& operator=(
        CouplingGeometry<TOtherPointType> const & rOther )
    {
        BaseType::operator=( rOther );
        mpGeometries = rOther.mpGeometries;
        return *this;
    }

    ///@}
    ///@name Operations
    ///@{

    typename BaseType::Pointer Create(
        PointsArrayType const& ThisPoints ) const override
    {
        return Kratos::make_shared<CouplingGeometry>();
    }

    ///@}
    ///@name Input and output
    ///@{

    /**
    * @brief This function returns the pointer of the geometry part
    *        which is corresponding to the index.
    *        Checks if index is available only in debug.
    * @param Index: 0 -> master, all bigger than 0 -> slaves.
    * @return pointer of geometry, corresponding to the index.
    */
    GeometryPointer pGetGeometryPart(const IndexType Index) override
    {
        KRATOS_DEBUG_ERROR_IF(mpGeometries.size() <= Index) << "Index "
            << Index << " out of range. CouplingGeometry #" << this->Id()
            << " has " << mpGeometries.size() << " geometries." << std::endl;

        return mpGeometries[Index];
    }

    /**
    * @brief This function returns the const pointer of the geometry part
    *        which is corresponding to the index.
    *        Checks if index is available only in debug.
    * @param Index: 0 -> master, all bigger than 0 -> slaves.
    * @return const pointer of geometry, corresponding to the index.
    */
    const GeometryPointer pGetGeometryPart(const IndexType Index) const override
    {
        KRATOS_DEBUG_ERROR_IF(mpGeometries.size() <= Index) << "Index \""
            << Index << "\" out of range. CouplingGeometry #" << this->Id()
            << " has " << mpGeometries.size() << " geometries." << std::endl;

        return mpGeometries[Index];
    }

    /**
     * @brief Allows to exchange certain geometries.
     * @param Index of the geometry part. 0->Master; 1->Slave
     * @param pGeometry The new geometry to add
     */
    void SetGeometryPart(
        const IndexType Index,
        GeometryPointer pGeometry
        ) override
    {
        KRATOS_DEBUG_ERROR_IF(mpGeometries.size() <= Index) << "Index out of range: "
            << Index << " composite contains only of: "
            << mpGeometries.size() << " geometries." << std::endl;

        if (Index == 0){
            this->SetGeometryData(&(pGeometry->GetGeometryData()));
        }

        mpGeometries[Index] = pGeometry;
    }

    /**
     * @brief Allows to enhance the coupling geometry, with another geometry.
     * @param pGeometry The new geometry to add
     */
    IndexType AddGeometryPart(GeometryPointer pGeometry) override
    {
        const IndexType new_index = mpGeometries.size();

        mpGeometries.push_back(pGeometry);

        return new_index;
    }

    /**
     * @brief Removes a geometry part
     * @param pGeometry The new geometry to remove
     */
    void RemoveGeometryPart(GeometryPointer pGeometry) override
    {
        const IndexType geometry_id = pGeometry->Id();
        IndexType to_remove_id = 0;
        for (const auto& p_geom : mpGeometries) {
            if (p_geom->Id() == geometry_id) {
                break;
            }
            ++to_remove_id;
        }

        RemoveGeometryPart(to_remove_id);
    }

    /**
     * @brief Removes a geometry part
     * @param Index of the geometry part.
     */
    void RemoveGeometryPart(const IndexType Index) override
    {
        KRATOS_ERROR_IF(Index == 0) << "Master geometry should not be removed from the CouplingGeometry" << std::endl;

        const SizeType number_of_geometries = NumberOfGeometryParts();
        for (IndexType i = Index; i < number_of_geometries - 1; ++i) {
            mpGeometries[i] = mpGeometries[i + 1];
        }
        mpGeometries[number_of_geometries - 1] = nullptr;
        mpGeometries.erase(mpGeometries.begin() + number_of_geometries - 1);
    }

    /**
     * @brief Use to check if certain Indexed object is within the geometry parts of this geometry.
     * @param Index of the geometry part. This index can be used differently within the derived classes.
     * @return true if has geometry part
     */
    bool HasGeometryPart(const IndexType Index) const override
    {
        if (Index < NumberOfGeometryParts()) {
            return true;
        } else {
            return false;
        }
    }

    /**
     * @brief The number of geometry part
     * @return The number of geometry parts that this geometry contains.
     */
    SizeType NumberOfGeometryParts() const override
    {
        return mpGeometries.size();
    }

    ///@}
    ///@name Geometry Information
    ///@{

    /**
     * @brief Returns the domain size of the master geometry.
     * @return The domaon size of the master geometry
     */
    double DomainSize() const override
    {
        KRATOS_DEBUG_ERROR_IF(mpGeometries.size() == 0) << "No master assigned. Geometry vector of size 0." << std::endl;

        return mpGeometries[0]->DomainSize();
    }

    /**
     * @brief Returns the center of the master geometry.
     * @return The center of the master geometry
     */
    Point Center() const override
    {
        KRATOS_DEBUG_ERROR_IF(mpGeometries.size() == 0) << "No master assigned. Geometry vector of size 0." << std::endl;

        return mpGeometries[0]->Center();
    }

    ///@}
    ///@name Spans
    ///@{

    /* @brief Provides the combined spans of all geometry parts of this geometry
     *        in local parameter coordinates of the geometry
     *        according to its direction from LocalDirectionIndex.
     *
     * @param resulting vector of span intervals.
     * @param LocalDirectionIndex of chosen direction, for curves always 0.
     */
    virtual void SpansLocalSpace(
        std::vector<double>& rSpans,
        IndexType LocalDirectionIndex = 0) const override
    {
        const double model_tolerance = 1e-2;

        if (this->LocalSpaceDimension() == 1) {
            std::vector<double> master_span_intersections_in_master_local_space;
            std::vector<double> slave_span_intersections_in_master_local_space;

            mpGeometries[0]->SpansLocalSpace(master_span_intersections_in_master_local_space);

            // Create tessellation for estimation of initial guesses.
            CurveTessellation<PointerVector<TPointType>> curve_tessellation_master;
            curve_tessellation_master.Tessellate(
                *(mpGeometries[0].get()),
                master_span_intersections_in_master_local_space,
                1e-2, mpGeometries[0]->PolynomialDegree(0), false);

            CoordinatesArrayType local_coords_span_intersection_on_slave = ZeroVector(3);
            CoordinatesArrayType global_coords_span_intersection_on_slave = ZeroVector(3);
            CoordinatesArrayType local_coords_master = ZeroVector(3);
            CoordinatesArrayType global_coords_master;
            for (IndexType i = 1; i < mpGeometries.size(); ++i) {
                std::vector<double> intersection_slave_spans;
                mpGeometries[i]->SpansLocalSpace(intersection_slave_spans);

                for (IndexType j = 0; j < intersection_slave_spans.size(); ++j) {
                    local_coords_span_intersection_on_slave[0] = intersection_slave_spans[j];
                    // Get global coordinates of span intersection on slave.
                    mpGeometries[i]->GlobalCoordinates(
                        global_coords_span_intersection_on_slave, local_coords_span_intersection_on_slave);
                    // Get initial guess for projection on master curve.
                    curve_tessellation_master.GetClosestPoint(
                        global_coords_span_intersection_on_slave, global_coords_master, local_coords_master);
                    // Projection on master curve.
                    int success = mpGeometries[0]->ProjectionPointGlobalToLocalSpace(
                        global_coords_span_intersection_on_slave, local_coords_master);

                    #ifdef KRATOS_DEBUG
                        mpGeometries[0]->GlobalCoordinates(global_coords_master, local_coords_master);
                    #endif
                    KRATOS_DEBUG_ERROR_IF((success > 1 && (norm_2(global_coords_span_intersection_on_slave - local_coords_master) > model_tolerance))
                        || (success == 0 && (norm_2(global_coords_span_intersection_on_slave - local_coords_master) < model_tolerance)))
                        << "Projection of intersection spans failed. Global Coordinates on slave: "
                        << global_coords_span_intersection_on_slave << ", and global coordinates on master: "
                        << global_coords_master << ". Difference: " << norm_2(global_coords_span_intersection_on_slave - global_coords_master)
                        << " larger than model tolerance: " << model_tolerance << std::endl;

                    // If success == 0, it is considered that the projection is on one of the boundaries.
                    slave_span_intersections_in_master_local_space.push_back(local_coords_master[0]);
                }
            }

            MergeSpans(rSpans, master_span_intersections_in_master_local_space, slave_span_intersections_in_master_local_space);
        }
    }

    ///@}
    ///@name Inquiry
    ///@{

    /// Provides the default integration per geometry.
    IntegrationInfo GetDefaultIntegrationInfo() const override
    {
        const SizeType local_space_dimension = this->LocalSpaceDimension();

        std::vector<SizeType> number_of_points_per_span_per_direction(local_space_dimension);
        std::vector<IntegrationInfo::QuadratureMethod> quadrature_method_per_direction(local_space_dimension);
        for (IndexType i = 0; i < local_space_dimension; ++i) {
            SizeType max_p = 0;
            for (IndexType j = 0; j < NumberOfGeometryParts(); ++j) {
                max_p = std::max(mpGeometries[j]->PolynomialDegree(i) + 1, max_p);
            }
            number_of_points_per_span_per_direction[i] = max_p;
            quadrature_method_per_direction[i] = IntegrationInfo::QuadratureMethod::GAUSS;
        }

        return IntegrationInfo(number_of_points_per_span_per_direction, quadrature_method_per_direction);
    }

    ///@}
    ///@name Integration Points
    ///@{

    /* Creates integration points on the master considering all slave intersections.
     * @return integration points.
     */
    void CreateIntegrationPoints(
        IntegrationPointsArrayType& rIntegrationPoints,
        IntegrationInfo& rIntegrationInfo) const override
    {
        if (this->LocalSpaceDimension() == 1) {
            std::vector<double> spans;
            this->SpansLocalSpace(spans, 0);

            IntegrationPointUtilities::CreateIntegrationPoints1D(
                rIntegrationPoints, spans, rIntegrationInfo);
        }
    }

    /* @brief This method creates a list of quadrature point geometries
     *        from a list of integration points.
     *        Initially, all integration points are transformed to
     *        quadrature point geometries on the master. Then all points
     *        are mapped to all slave geometries. Finally, coupling geometries,
     *        containing quadrature point geometries per integration point on
     *        master and all slaves are created.
     *
     * @param rResultGeometries list of quadrature point geometries.
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
        const double model_tolerance = 1e-2;

        const SizeType num_integration_points = rIntegrationPoints.size();

        if (rResultGeometries.size() != num_integration_points) {
            rResultGeometries.resize(num_integration_points);
        }

        // Create quadrature points on master
        GeometriesArrayType master_quadrature_points(num_integration_points);
        mpGeometries[0]->CreateQuadraturePointGeometries(
            master_quadrature_points,
            NumberOfShapeFunctionDerivatives,
            rIntegrationPoints,
            rIntegrationInfo);

        // Compute vector of location
        CoordinatesArrayVectorType integration_points_global_coords_vector(num_integration_points);
        for (SizeType i = 0; i < num_integration_points; ++i) {
            integration_points_global_coords_vector[i] = master_quadrature_points[i].Center();
        }

        // First slave
        IntegrationPointsArrayType integration_points_slave = rIntegrationPoints;
        CoordinatesArrayType local_slave_coords = ZeroVector(3);
        CoordinatesArrayType global_slave_coords = ZeroVector(3);

        if (!rIntegrationInfo.Is(IntegrationInfo::DO_NOT_CREATE_TESSELLATION_ON_SLAVE)) {
            if (this->LocalSpaceDimension() == 1) {
                CurveTessellation<PointerVector<TPointType>> curve_tesselation;
                curve_tesselation.Tessellate(*(mpGeometries[1].get()), 1e-2, mpGeometries[1]->PolynomialDegree(0));

                for (SizeType j = 0; j < num_integration_points; ++j) {
                    curve_tesselation.GetClosestPoint(
                        integration_points_global_coords_vector[j],
                        global_slave_coords,
                        local_slave_coords);

                    mpGeometries[1]->ProjectionPointGlobalToLocalSpace(
                        integration_points_global_coords_vector[j],
                        local_slave_coords);

                    integration_points_slave[j][0] = local_slave_coords[0];
                    integration_points_slave[j][1] = local_slave_coords[1];
                    integration_points_slave[j][2] = local_slave_coords[2];
                }
            }
            else {
                KRATOS_ERROR << "Tessellation for " << this->LocalSpaceDimension()
                    << "-dimensional objects is not implemented." << std::endl;
            }
        }
        else {
            for (SizeType j = 0; j < num_integration_points; ++j) {
                mpGeometries[1]->ProjectionPointGlobalToLocalSpace(
                    integration_points_global_coords_vector[j],
                    local_slave_coords);

                integration_points_slave[j][0] = local_slave_coords[0];
                integration_points_slave[j][1] = local_slave_coords[1];
                integration_points_slave[j][2] = local_slave_coords[2];
            }
        }

        GeometriesArrayType slave_quadrature_points(num_integration_points);
        mpGeometries[1]->CreateQuadraturePointGeometries(
            slave_quadrature_points,
            NumberOfShapeFunctionDerivatives,
            integration_points_slave,
            rIntegrationInfo);

        for (SizeType i = 0; i < num_integration_points; ++i) {
            KRATOS_DEBUG_ERROR_IF(norm_2(master_quadrature_points(i)->Center() - slave_quadrature_points(i)->Center()) > model_tolerance)
                << "Difference between master and slave coordinates above model tolerance of " << model_tolerance
                << ". Location of master: " << master_quadrature_points(i)->Center() << ", location of slave: "
                << slave_quadrature_points(i)->Center() << ". Distance: "
                << norm_2(master_quadrature_points(i)->Center() - slave_quadrature_points(i)->Center()) << std::endl;

            rResultGeometries(i) = Kratos::make_shared<CouplingGeometry<PointType>>(
                master_quadrature_points(i), slave_quadrature_points(i));
        }

        KRATOS_ERROR_IF(mpGeometries.size() > 2)
            << "CreateQuadraturePointGeometries not implemented for coupling of more than 2 geometries. "
            << mpGeometries.size() << " are given." << std::endl;
    }

    void CreateQuadraturePointGeometries(
        GeometriesArrayType& rResultGeometries,
        IndexType NumberOfShapeFunctionDerivatives,
        IntegrationInfo& rIntegrationInfo) override
    {
        const double model_tolerance = 1e-3;

        if (this->LocalSpaceDimension() != 0) {
            BaseType::CreateQuadraturePointGeometries(rResultGeometries, NumberOfShapeFunctionDerivatives, rIntegrationInfo);
        }
        else {
            rResultGeometries.resize(1);

            GeometriesArrayType master_quadrature_points(1);
            mpGeometries[0]->CreateQuadraturePointGeometries(
                master_quadrature_points,
                NumberOfShapeFunctionDerivatives,
                rIntegrationInfo);

            GeometriesArrayType slave_quadrature_points(1);
            mpGeometries[1]->CreateQuadraturePointGeometries(
                slave_quadrature_points,
                NumberOfShapeFunctionDerivatives,
                rIntegrationInfo);

            KRATOS_DEBUG_ERROR_IF(norm_2(master_quadrature_points(0)->Center() - slave_quadrature_points(0)->Center()) > model_tolerance)
                << "Difference between master and slave coordinates above model tolerance of " << model_tolerance
                << ". Location of master: " << master_quadrature_points(0)->Center() << ", location of slave: "
                << slave_quadrature_points(0)->Center() << ". Distance: "
                << norm_2(master_quadrature_points(0)->Center() - slave_quadrature_points(0)->Center()) << std::endl;

            rResultGeometries(0) = Kratos::make_shared<CouplingGeometry<PointType>>(
                master_quadrature_points(0), slave_quadrature_points(0));

            if (mpGeometries.size() > 2) {
                for (IndexType i = 2; i < mpGeometries.size(); ++i) {
                    GeometriesArrayType more_slave_quadrature_points(1);
                    mpGeometries[i]->CreateQuadraturePointGeometries(
                        more_slave_quadrature_points,
                        NumberOfShapeFunctionDerivatives,
                        rIntegrationInfo);

                    KRATOS_DEBUG_ERROR_IF(norm_2(master_quadrature_points(0)->Center() - more_slave_quadrature_points(0)->Center()) > model_tolerance)
                        << "Difference between master and slave coordinates above model tolerance of " << model_tolerance
                        << ". Location of master: " << master_quadrature_points(0)->Center() << ", location of slave: "
                        << more_slave_quadrature_points(0)->Center() << ". Distance: "
                        << norm_2(master_quadrature_points(0)->Center() - more_slave_quadrature_points(0)->Center()) << std::endl;

                    rResultGeometries(0)->AddGeometryPart(more_slave_quadrature_points(0));
                }
            }
        }
    }

    ///@}
    ///@name Span Utilities
    ///@{

    static void MergeSpans(
        std::vector<double>& rResultSpans,
        const std::vector<double>& rSpans1,
        const std::vector<double>& rSpans2,
        double Tolerance = 1e-6)
    {
        NurbsInterval interval_1(rSpans1[0], rSpans1[rSpans1.size() - 1]);
        NurbsInterval interval_2(rSpans2[0], rSpans2[rSpans2.size() - 1]);

        for (IndexType i = 0; i < rSpans1.size(); ++i) {
            double temp = rSpans1[i];
            interval_2.IsInside(temp);
            rResultSpans.push_back(temp);
        }
        for (IndexType i = 0; i < rSpans2.size(); ++i) {
            double temp = rSpans2[i];
            interval_1.IsInside(temp);
            rResultSpans.push_back(temp);
        }

        SortUnique(rResultSpans, Tolerance);
    }

    static void SortUnique(
        std::vector<double>& rResultSpans,
        const double Tolerance)
    {
        std::sort(std::begin(rResultSpans), std::end(rResultSpans));

        auto last = std::unique(std::begin(rResultSpans), std::end(rResultSpans),
            [Tolerance](double a, double b) { return b - a < Tolerance; });

        auto nb_unique = std::distance(std::begin(rResultSpans), last);

        rResultSpans.resize(nb_unique);
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
        return GeometryData::KratosGeometryFamily::Kratos_Composite;
    }

    /**
     * @brief Gets the geometry type.
     * @details This function returns the specific type of the geometry. The geometry type provides a more detailed classification of the geometry.
     * @return GeometryData::KratosGeometryType The specific geometry type.
     */
    GeometryData::KratosGeometryType GetGeometryType() const override
    {
        return GeometryData::KratosGeometryType::Kratos_Coupling_Geometry;
    }

    ///@}
    ///@name Information
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "Coupling geometry that holds a master and a set of slave geometries.";
    }

    /// Print information about this object.
    void PrintInfo( std::ostream& rOStream ) const override
    {
        rOStream << "Coupling geometry that holds a master and a set of slave geometries.";
    }

    /// Print object's data.
    void PrintData( std::ostream& rOStream ) const override
    {
        BaseType::PrintData( rOStream );
        std::cout << std::endl;
        rOStream << "    CouplingGeometry with " << mpGeometries.size() << " geometries.";
    }

    ///@}
private:
    ///@name Private Member Variables
    ///@{

    GeometryPointerVector mpGeometries;

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save( Serializer& rSerializer ) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, BaseType );
        rSerializer.save("Geometries", mpGeometries);
    }

    void load( Serializer& rSerializer ) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BaseType );
        rSerializer.load("Geometries", mpGeometries);
    }

    ///@}
    ///@name Private Friends
    ///@{

    template<class TOtherPointType> friend class CouplingGeometry;

    ///@}
}; // Class Geometry

///@}
///@name Input and output
///@{
/**
 * input stream functions
 */
template<class TPointType> inline std::istream& operator >> (
    std::istream& rIStream,
    CouplingGeometry<TPointType>& rThis );
/**
 * output stream functions
 */
template<class TPointType> inline std::ostream& operator << (
    std::ostream& rOStream,
    const CouplingGeometry<TPointType>& rThis )
{
    rThis.PrintInfo( rOStream );
    rOStream << std::endl;
    rThis.PrintData( rOStream );
    return rOStream;
}

///@}
}// namespace Kratos.

#endif // KRATOS_COUPLING_GEOMETRY_H_INCLUDED  defined

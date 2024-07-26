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

    static constexpr IndexType MasterIndex = 0;
    static constexpr IndexType SlaveIndex = 1;


    ///@}
    ///@name Life Cycle
    ///@{

    /// constructor for untrimmed surface
    NurbsCouplingGeometry2D(
        GeometriesArrayType& slaveGeometryList,
        GeometriesArrayType& masterGeometryList)
        : mSlaveGeometryList(slaveGeometryList)
        , mMasterGeometryList(masterGeometryList)
    {

        mpNurbsSurfaceMaster = std::dynamic_pointer_cast<NurbsSurfaceType>(mMasterGeometryList[0].pGetGeometryPart(GeometryType::BACKGROUND_GEOMETRY_INDEX));
        mpNurbsSurfaceSlave = std::dynamic_pointer_cast<NurbsSurfaceType>(mSlaveGeometryList[0].pGetGeometryPart(GeometryType::BACKGROUND_GEOMETRY_INDEX));
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
    GeometriesArrayType pGetSlaveGeometry()
    {
        return mSlaveGeometryList;
    }

    /// Returns the Master Geometry list of this coupling geometry.
    GeometriesArrayType pGetMasterGeometry()
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
    static void SpansLocalSpaceForParentIntegration(std::vector<std::vector<double>>& integration_edges_on_parameter_parent_list,
                                             std::vector<std::vector<double>>& spans_parent_list,
                                             std::vector<std::vector<double>>& spans_paired_list,
                                             GeometriesArrayType& rParentGeometryList,
                                             GeometriesArrayType& rPairedGeometryList)
    {
        for (int i_brep_parent = 0; i_brep_parent < spans_parent_list.size(); i_brep_parent++) {
            integration_edges_on_parameter_parent_list[i_brep_parent] = spans_parent_list[i_brep_parent];
        }

        for (int i_brep_paired = 0; i_brep_paired < spans_paired_list.size(); i_brep_paired++) {
            for (int i = 0; i < spans_paired_list[i_brep_paired].size(); i++) {
                
                CoordinatesArrayType locCoord(3); locCoord[0] = spans_paired_list[i_brep_paired][i];
                
                double incumb_distance = 1e16;
                double projection_parameter_in_best_parent_brep;
                int best_parent_brep_index = -1; //initialize at impossible value
                for (int i_brep_parent = 0; i_brep_parent < spans_parent_list.size(); i_brep_parent++) {
                    Vector interval; 
                    rParentGeometryList[i_brep_parent].DomainInterval(interval);


                    double current_distance;
                    double current_projection_parameter_in_parent;
                    bool hasFoundProjection = GetProjection(locCoord, rPairedGeometryList[i_brep_paired], 
                                                            rParentGeometryList[i_brep_parent], current_projection_parameter_in_parent, current_distance);

                    if (hasFoundProjection && current_distance < incumb_distance &&
                        interval[0] <= current_projection_parameter_in_parent && current_projection_parameter_in_parent <= interval[1]) {
                        incumb_distance = current_distance;
                        best_parent_brep_index = i_brep_parent;
                        projection_parameter_in_best_parent_brep = current_projection_parameter_in_parent;
                    }
                } 
                
                if (best_parent_brep_index > -1) integration_edges_on_parameter_parent_list[best_parent_brep_index].push_back(projection_parameter_in_best_parent_brep);
            }
        }

        // maybe we can do better
        for (int i_brep_parent = 0; i_brep_parent < spans_parent_list.size(); i_brep_parent++) {
            std::sort(integration_edges_on_parameter_parent_list[i_brep_parent].begin(), integration_edges_on_parameter_parent_list[i_brep_parent].end());
            auto last = std::unique(integration_edges_on_parameter_parent_list[i_brep_parent].begin(), integration_edges_on_parameter_parent_list[i_brep_parent].end());
            integration_edges_on_parameter_parent_list[i_brep_parent].erase(last, integration_edges_on_parameter_parent_list[i_brep_parent].end());
        }
            
    }

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

            max_p = std::max(mSlaveGeometryList[0].PolynomialDegree(i) + 1, max_p);
            max_p = std::max(mMasterGeometryList[0].PolynomialDegree(i) + 1, max_p);

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
    static void CreateIntegrationPoints(
        std::vector<IntegrationPointsArrayType>& rIntegrationPoints,
        IntegrationInfo& rIntegrationInfo,
        std::vector<std::vector<double>>& spans_parent_list,
        std::vector<std::vector<double>>& spans_paired_list,
        GeometriesArrayType& rParentGeometryList,
        GeometriesArrayType& rPairedGeometryList
        )
    {
        std::vector<std::vector<double>> integration_edges_on_parameter_parent_list(spans_parent_list.size());
        SpansLocalSpaceForParentIntegration(integration_edges_on_parameter_parent_list, spans_parent_list, spans_paired_list, rParentGeometryList, rPairedGeometryList);
        
        for (int i_brep_parent = 0; i_brep_parent < spans_parent_list.size(); i_brep_parent++){
            IntegrationPointUtilities::CreateIntegrationPoints1D(
                rIntegrationPoints[i_brep_parent], integration_edges_on_parameter_parent_list[i_brep_parent], rIntegrationInfo);
        }

    }
    
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
        IntegrationInfo& rIntegrationInfo) override
    {
        std::vector<std::vector<double>> spans_slave_list(mSlaveGeometryList.size());
        for (int i = 0; i < mSlaveGeometryList.size(); i++) mSlaveGeometryList[i].SpansLocalSpace(spans_slave_list[i]);

        std::vector<std::vector<double>> spans_master_list(mMasterGeometryList.size());
        for (int i = 0; i < mMasterGeometryList.size(); i++) mMasterGeometryList[i].SpansLocalSpace(spans_master_list[i]);


        // for master
        std::vector<IntegrationPointsArrayType> IntegrationPointsMaster(mMasterGeometryList.size());
        CreateIntegrationPoints(IntegrationPointsMaster, rIntegrationInfo, spans_master_list, spans_slave_list, 
                                mMasterGeometryList, mSlaveGeometryList);


         // // for slave
        std::vector<IntegrationPointsArrayType> IntegrationPointsSlave(mSlaveGeometryList.size());
        CreateIntegrationPoints(IntegrationPointsSlave, rIntegrationInfo, spans_slave_list, spans_master_list, 
                                mSlaveGeometryList, mMasterGeometryList);


        // // Resize containers.
        int numberIntegrationPointsMaster = 0; 
        int numberIntegrationPointsSlave = 0; 
        for (int i = 0; i < mMasterGeometryList.size(); i++) numberIntegrationPointsMaster += IntegrationPointsMaster[i].size();  
        for (int i = 0; i < mSlaveGeometryList.size(); i++) numberIntegrationPointsSlave += IntegrationPointsSlave[i].size(); 

        const int numberIntegrationPoints = numberIntegrationPointsMaster+numberIntegrationPointsSlave;
        if (rResultGeometries.size() != numberIntegrationPoints)
            rResultGeometries.resize(numberIntegrationPoints);


        SizeType quadraturePointId = 0;
        this->CreateQuadraturePointGeometriesOnParent(
            rResultGeometries,
            NumberOfShapeFunctionDerivatives,
            IntegrationPointsMaster,
            rIntegrationInfo,
            mpNurbsSurfaceMaster,
            mpNurbsSurfaceSlave,
            mMasterGeometryList,
            mSlaveGeometryList,
            MasterIndex,
            quadraturePointId);

       
        GeometriesArrayType rResultGeometries2;
        this->CreateQuadraturePointGeometriesOnParent(
            rResultGeometries,
            NumberOfShapeFunctionDerivatives,
            IntegrationPointsSlave,
            rIntegrationInfo,
            mpNurbsSurfaceSlave,
            mpNurbsSurfaceMaster,
            mSlaveGeometryList,
            mMasterGeometryList,
            SlaveIndex,
            quadraturePointId);

        rResultGeometries.resize(quadraturePointId);
        

    }


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
        GeometriesArrayType& rParentGeometryList,
        GeometriesArrayType& rPairedGeometryList,
        const IndexType& integrationDomain,
        SizeType& quadraturePointId) 
    {
        // FOR EACH OF THE INTEGRATION POINTS FIND THE PROJECTION ON THE SLAVE, AND SAVE ALL THE VALUES OF THE BASIS FUNCTIONS 
        
        // shape function container.
        NurbsSurfaceShapeFunction shape_function_container_master(
            rpNurbsSurfaceParent->PolynomialDegreeU(), rpNurbsSurfaceParent->PolynomialDegreeV(),
            NumberOfShapeFunctionDerivatives);

        NurbsSurfaceShapeFunction shape_function_container_slave(
            rpNurbsSurfacePaired->PolynomialDegreeU(), rpNurbsSurfacePaired->PolynomialDegreeV(),
            NumberOfShapeFunctionDerivatives);

        
        auto default_method = this->GetDefaultIntegrationMethod();
        SizeType num_nonzero_cps_master = shape_function_container_master.NumberOfNonzeroControlPoints();

        Matrix N_m(1, num_nonzero_cps_master);
        DenseVector<Matrix> shape_function_derivatives_m(NumberOfShapeFunctionDerivatives - 1);
        for (IndexType i = 0; i < NumberOfShapeFunctionDerivatives - 1; i++) {
            shape_function_derivatives_m[i].resize(num_nonzero_cps_master, i + 2);
        }
        //--------------
        SizeType num_nonzero_cps_slave = shape_function_container_slave.NumberOfNonzeroControlPoints();

        Matrix N_s(1, num_nonzero_cps_slave);
        DenseVector<Matrix> shape_function_derivatives_s(NumberOfShapeFunctionDerivatives - 1);
        for (IndexType i = 0; i < NumberOfShapeFunctionDerivatives - 1; i++) {
            shape_function_derivatives_s[i].resize(num_nonzero_cps_slave, i + 2);
        }

        for (IndexType i_brep_m = 0; i_brep_m < rIntegrationPoints.size(); i_brep_m++) {
            
            std::vector<CoordinatesArrayType> global_space_first(2);
            std::vector<CoordinatesArrayType> global_space_second(2);
            rParentGeometryList[i_brep_m].GlobalSpaceDerivatives(global_space_first ,rIntegrationPoints[i_brep_m][0], 1);  // i = 0
            rParentGeometryList[i_brep_m].GlobalSpaceDerivatives(global_space_second,rIntegrationPoints[i_brep_m][1], 1); // i = 1

            //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
            for (IndexType i = 0; i < rIntegrationPoints[i_brep_m].size(); ++i)
            {  
                // MASTER
                //ÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑ
                std::vector<CoordinatesArrayType> global_space_derivatives_master(2);
                rParentGeometryList[i_brep_m].GlobalSpaceDerivatives(
                        global_space_derivatives_master,
                        rIntegrationPoints[i_brep_m][i],
                        1);

                //***+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
                std::ofstream outputFile("txt_files/Gauss_Point_Contact_coordinates.txt", std::ios::app);
                if (!outputFile.is_open())
                {
                    std::cerr << "Failed to open the file for writing." << std::endl;
                    return;
                }
                outputFile << std::setprecision(14); // Set precision to 10^-14
                outputFile << global_space_derivatives_master[0][0] << "  " << global_space_derivatives_master[0][1]   <<"\n";
                outputFile.close();
                //++++++++++++++++++++++++++++++++++++++++
                if (rpNurbsSurfaceParent->IsRational()) {
                    shape_function_container_master.ComputeNurbsShapeFunctionValues(
                        rpNurbsSurfaceParent->KnotsU(), rpNurbsSurfaceParent->KnotsV(), rpNurbsSurfaceParent->Weights(),
                        global_space_derivatives_master[0][0], global_space_derivatives_master[0][1]);
                }
                else {
                    // IN ORDER TO CHECK IF YOU ARE USING TRIM OR SBM APPROACH
                    std::ifstream file("txt_files/input_data.txt");
                    std::string line;
                    int SBM_technique;
                    std::getline(file, line);
                    std::getline(file, line); // Read the second line
                    SBM_technique = std::stoi(line);
                    file.close();

                    // MODIFIED
                    bool is_surrogate_boundary = true; 
                    const double threshold = 1e-14 ;
                    // At this point all the true are boundary GPs
                    if (std::abs(global_space_derivatives_master[0][0]-rpNurbsSurfaceParent->KnotsU()[0]) < threshold) {is_surrogate_boundary = false;} // External boundary
                    if (std::abs(global_space_derivatives_master[0][0]-rpNurbsSurfaceParent->KnotsU()[rpNurbsSurfaceParent->KnotsU().size()-1]) < threshold) {is_surrogate_boundary = false;} // External boundary
                    if (std::abs(global_space_derivatives_master[0][1]-rpNurbsSurfaceParent->KnotsV()[0]) < threshold) {is_surrogate_boundary = false;} // External boundary
                    if (std::abs(global_space_derivatives_master[0][1]-rpNurbsSurfaceParent->KnotsV()[rpNurbsSurfaceParent->KnotsV().size()-1]) < threshold) {is_surrogate_boundary = false;} // External boundary
                    
                    if (SBM_technique==1) {is_surrogate_boundary=false;}

                    // If it is 1 -> the current GP is a surrogate GP otherwise it is an external GP
                    if (is_surrogate_boundary) {
                        // KRATOS_WATCH("SURROGATE BOUNDARY GP")
                        IndexType SpanU = NurbsUtilities::GetLowerSpan(rpNurbsSurfaceParent->PolynomialDegreeU(), rpNurbsSurfaceParent->KnotsU(), global_space_derivatives_master[0][0]);
                        IndexType SpanV = NurbsUtilities::GetLowerSpan(rpNurbsSurfaceParent->PolynomialDegreeV(), rpNurbsSurfaceParent->KnotsV(), global_space_derivatives_master[0][1]);

                        // Understand if the knot-boundary is along U or V
                        bool gp_is_along_U;
                        if (std::abs(global_space_first[0][0] - global_space_second[0][0]) < threshold) {
                            gp_is_along_U = true;
                        }
                        else {gp_is_along_U = false; }

                        if (gp_is_along_U && global_space_first[0][1] > global_space_second[0][1]) {
                            // Is decreasing along y -> Add 1 at SpanU
                            SpanU++;
                        }
                        if (!gp_is_along_U && global_space_first[0][0] < global_space_second[0][0]) {
                            // Is increasing along x -> Add 1 at SpanV
                            SpanV++;
                        }

                        shape_function_container_master.ComputeBSplineShapeFunctionValuesAtSpan(
                            rpNurbsSurfaceParent->KnotsU(),
                            rpNurbsSurfaceParent->KnotsV(),
                            SpanU,
                            SpanV,
                            global_space_derivatives_master[0][0],
                            global_space_derivatives_master[0][1]);

                    }
                    else {
                        // 'EXTERNAL GP'
                        shape_function_container_master.ComputeBSplineShapeFunctionValues(
                            rpNurbsSurfaceParent->KnotsU(), rpNurbsSurfaceParent->KnotsV(),
                            global_space_derivatives_master[0][0], global_space_derivatives_master[0][1]);
                    }
                }

                /// Get List of Control Points
                PointsArrayType nonzero_control_points_master(num_nonzero_cps_master);
                auto cp_indices_master = shape_function_container_master.ControlPointIndices(
                    rpNurbsSurfaceParent->NumberOfControlPointsU(), rpNurbsSurfaceParent->NumberOfControlPointsV());
                for (IndexType j = 0; j < num_nonzero_cps_master; j++) {
                    nonzero_control_points_master(j) = rpNurbsSurfaceParent->pGetPoint(cp_indices_master[j]);
                }
                /// Get Shape Functions N
                for (IndexType j = 0; j < num_nonzero_cps_master; j++) {
                    N_m(0, j) = shape_function_container_master(j, 0);
                }

                /// Get Shape Function Derivatives DN_De, ...
                if (NumberOfShapeFunctionDerivatives > 0) {
                    IndexType shape_derivative_index = 1;
                    for (IndexType n = 0; n < NumberOfShapeFunctionDerivatives - 1; n++) {
                        for (IndexType k = 0; k < n + 2; k++) {
                            for (IndexType j = 0; j < num_nonzero_cps_master; j++) {
                                shape_function_derivatives_m[n](j, k) = shape_function_container_master(j, shape_derivative_index + k);
                            }
                        }
                        shape_derivative_index += n + 2;
                    }
                }

                GeometryShapeFunctionContainer<GeometryData::IntegrationMethod> data_container_master(
                    default_method, rIntegrationPoints[i_brep_m][i],
                    N_m, shape_function_derivatives_m);



                // SLAVE
                //ÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑ    

                // project point 

                
                CoordinatesArrayType master_quadrature_point(3);
                master_quadrature_point[0] = global_space_derivatives_master[0][0]; master_quadrature_point[1] = global_space_derivatives_master[0][1];
                master_quadrature_point[2] = 0.0;
                double best_distance = 1e16;
                CoordinatesArrayType best_projected_point_on_slave_local;
                int best_bred_id_slave = -1;
                bool isConverged;
                for (IndexType i_brep_s = 0; i_brep_s < rIntegrationPoints.size(); i_brep_s++) {
                    
                    CoordinatesArrayType local_coord_projected_on_slave; //first trial
                    isConverged = rPairedGeometryList[i_brep_s].ProjectionPointGlobalToLocalSpace(master_quadrature_point, local_coord_projected_on_slave, 1e-12);

                    if (isConverged) {
                        CoordinatesArrayType point_projected_on_slave(3);
                        
                        std::vector<CoordinatesArrayType> global_space_derivatives_slave(2);
                        rPairedGeometryList[i_brep_s].GlobalCoordinates(point_projected_on_slave, local_coord_projected_on_slave);

                        double current_distance = norm_2(master_quadrature_point-point_projected_on_slave);

                        if (current_distance < best_distance) {
                            best_distance = current_distance;
                            best_projected_point_on_slave_local = local_coord_projected_on_slave;
                            best_bred_id_slave = i_brep_s;
                        }
                    } else {
                        Vector interval; 
                        rPairedGeometryList[i_brep_s].DomainInterval(interval);
                        local_coord_projected_on_slave[0] = interval[0];
                        CoordinatesArrayType point_projected_on_slave(3);        
                        rPairedGeometryList[i_brep_s].GlobalCoordinates(point_projected_on_slave, local_coord_projected_on_slave);

                        double current_distance = norm_2(master_quadrature_point-point_projected_on_slave);

                        if (current_distance < best_distance) {
                            best_distance = current_distance;
                            best_projected_point_on_slave_local = local_coord_projected_on_slave;
                            best_bred_id_slave = i_brep_s;
                        }

                        local_coord_projected_on_slave[0] = interval[1];
                        rPairedGeometryList[i_brep_s].GlobalCoordinates(point_projected_on_slave, local_coord_projected_on_slave);

                        current_distance = norm_2(master_quadrature_point-point_projected_on_slave);


                        if (current_distance < best_distance) {
                            best_distance = current_distance;
                            best_projected_point_on_slave_local = local_coord_projected_on_slave;
                            best_bred_id_slave = i_brep_s;
                        }
                    }

                }
                if (!isConverged) continue;
                // WARNING ! REMOVE

                IntegrationPoint<1> integrationPointSlave(best_projected_point_on_slave_local[0]);

                std::vector<CoordinatesArrayType> global_space_derivatives_slave(2);
                rPairedGeometryList[best_bred_id_slave].GlobalSpaceDerivatives(
                        global_space_derivatives_slave,
                        integrationPointSlave,
                        1);

                //***+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
                std::ofstream outputFile2("txt_files/Gauss_Point_Contact_coordinates.txt", std::ios::app);
                if (!outputFile2.is_open())
                {
                    std::cerr << "Failed to open the file for writing." << std::endl;
                    return;
                }
                outputFile2 << std::setprecision(14); // Set precision to 10^-14
                outputFile2 << global_space_derivatives_slave[0][0] << "  " << global_space_derivatives_slave[0][1]   <<"\n";
                outputFile2.close();
                //++++++++++++++++++++++++++++++++++++++++
                if (rpNurbsSurfacePaired->IsRational()) {
                    shape_function_container_slave.ComputeNurbsShapeFunctionValues(
                        rpNurbsSurfacePaired->KnotsU(), rpNurbsSurfacePaired->KnotsV(), rpNurbsSurfacePaired->Weights(),
                        global_space_derivatives_slave[0][0], global_space_derivatives_slave[0][1]);
                }
                else {
                    // IN ORDER TO CHECK IF YOU ARE USING TRIM OR SBM APPROACH
                    std::ifstream file("txt_files/input_data.txt");
                    std::string line;
                    int SBM_technique;
                    std::getline(file, line);
                    std::getline(file, line); // Read the second line
                    SBM_technique = std::stoi(line);
                    file.close();

                    // MODIFIED
                    bool is_surrogate_boundary = true; 
                    const double threshold = 1e-14 ;
                    // At this point all the true are boundary GPs
                    if (std::abs(global_space_derivatives_slave[0][0]-rpNurbsSurfacePaired->KnotsU()[0]) < threshold) {is_surrogate_boundary = false;} // External boundary
                    if (std::abs(global_space_derivatives_slave[0][0]-rpNurbsSurfacePaired->KnotsU()[rpNurbsSurfacePaired->KnotsU().size()-1]) < threshold) {is_surrogate_boundary = false;} // External boundary
                    if (std::abs(global_space_derivatives_slave[0][1]-rpNurbsSurfacePaired->KnotsV()[0]) < threshold) {is_surrogate_boundary = false;} // External boundary
                    if (std::abs(global_space_derivatives_slave[0][1]-rpNurbsSurfacePaired->KnotsV()[rpNurbsSurfacePaired->KnotsV().size()-1]) < threshold) {is_surrogate_boundary = false;} // External boundary
                    
                    if (SBM_technique==1) {is_surrogate_boundary=false;}

                    // If it is 1 -> the current GP is a surrogate GP otherwise it is an external GP
                    if (is_surrogate_boundary) {
                        // KRATOS_WATCH("SURROGATE BOUNDARY GP")
                        IndexType SpanU = NurbsUtilities::GetLowerSpan(rpNurbsSurfacePaired->PolynomialDegreeU(), rpNurbsSurfacePaired->KnotsU(), global_space_derivatives_slave[0][0]);
                        IndexType SpanV = NurbsUtilities::GetLowerSpan(rpNurbsSurfacePaired->PolynomialDegreeV(), rpNurbsSurfacePaired->KnotsV(), global_space_derivatives_slave[0][1]);

                        // Understand if the knot-boundary is along U or V
                        bool gp_is_along_U;
                        if (std::abs(global_space_first[0][0] - global_space_second[0][0]) < threshold) {
                            gp_is_along_U = true;
                        }
                        else {gp_is_along_U = false; }

                        if (gp_is_along_U && global_space_first[0][1] > global_space_second[0][1]) {
                            // Is decreasing along y -> Add 1 at SpanU
                            SpanU++;
                        }
                        if (!gp_is_along_U && global_space_first[0][0] < global_space_second[0][0]) {
                            // Is increasing along x -> Add 1 at SpanV
                            SpanV++;
                        }

                        shape_function_container_slave.ComputeBSplineShapeFunctionValuesAtSpan(
                            rpNurbsSurfacePaired->KnotsU(),
                            rpNurbsSurfacePaired->KnotsV(),
                            SpanU,
                            SpanV,
                            global_space_derivatives_slave[0][0],
                            global_space_derivatives_slave[0][1]);

                    }
                    else {
                        // 'EXTERNAL GP'
                        shape_function_container_slave.ComputeBSplineShapeFunctionValues(
                            rpNurbsSurfacePaired->KnotsU(), rpNurbsSurfacePaired->KnotsV(),
                            global_space_derivatives_slave[0][0], global_space_derivatives_slave[0][1]);
                    }
                }

                /// Get List of Control Points
                PointsArrayType nonzero_control_points_slave(num_nonzero_cps_slave);
                auto cp_indices_slave = shape_function_container_slave.ControlPointIndices(
                    rpNurbsSurfacePaired->NumberOfControlPointsU(), rpNurbsSurfacePaired->NumberOfControlPointsV());
                for (IndexType j = 0; j < num_nonzero_cps_slave; j++) {
                    nonzero_control_points_slave(j) = rpNurbsSurfacePaired->pGetPoint(cp_indices_slave[j]);
                }
                /// Get Shape Functions N
                for (IndexType j = 0; j < num_nonzero_cps_slave; j++) {
                    N_s(0, j) = shape_function_container_slave(j, 0);
                }

                /// Get Shape Function Derivatives DN_De, ...
                if (NumberOfShapeFunctionDerivatives > 0) {
                    IndexType shape_derivative_index = 1;
                    for (IndexType n = 0; n < NumberOfShapeFunctionDerivatives - 1; n++) {
                        for (IndexType k = 0; k < n + 2; k++) {
                            for (IndexType j = 0; j < num_nonzero_cps_slave; j++) {
                                shape_function_derivatives_s[n](j, k) = shape_function_container_slave(j, shape_derivative_index + k);
                            }
                        }
                        shape_derivative_index += n + 2;
                    }
                }

                GeometryShapeFunctionContainer<GeometryData::IntegrationMethod> data_container_slave(
                    default_method, integrationPointSlave,
                    N_s, shape_function_derivatives_s);

                rResultGeometries(quadraturePointId) = CreateQuadraturePointsUtility<NodeType>::CreateQuadraturePointCouplingGeometry2D(
                    data_container_master, data_container_slave, 
                    nonzero_control_points_master, nonzero_control_points_slave,
                    global_space_derivatives_master[1][0], global_space_derivatives_master[1][1], 
                    global_space_derivatives_slave[1][0], global_space_derivatives_slave[1][1],
                    &rParentGeometryList[i_brep_m], &rPairedGeometryList[best_bred_id_slave], this);

                // set value of the background integration domain
                rResultGeometries(quadraturePointId)->SetValue(ACTIVATION_LEVEL, integrationDomain);
                quadraturePointId++;
            }
        }
    }


    void CreateQuadraturePointGeometriesSlave(
        GeometriesArrayType& rResultGeometries,
        IndexType NumberOfShapeFunctionDerivatives,
        const IntegrationPointsArrayType& rIntegrationPoints,
        IntegrationInfo& rIntegrationInfo) 
    {
        // mpCurveOnSurface->CreateQuadraturePointGeometries(
        //     rResultGeometries, NumberOfShapeFunctionDerivatives, rIntegrationPoints, rIntegrationInfo);

        // for (IndexType i = 0; i < rResultGeometries.size(); ++i) {
        //     rResultGeometries(i)->SetGeometryParent(this);
        // }
    }

    ///@}

    ///@name Projection functionalaties
    ///@{
    static bool GetProjection(CoordinatesArrayType& slavePointLocalCoord, GeometryType &slave_geometry, GeometryType &master_geometry, double& localProjection, double& distance) {
        
        CoordinatesArrayType slavePoint;
        slave_geometry.GlobalCoordinates(slavePoint, slavePointLocalCoord);

        // find slave normal in that point
        CoordinatesArrayType normal = slave_geometry.UnitNormal(slavePointLocalCoord);

        // KRATOS_WATCH(slavePoint)
        // KRATOS_WATCH(normal)

        // ANOTHER WAY -> MINIMUM DISTANCE OF ALL POINTS IN MASTER FROM POINT IN SLAVE, AND NOT THE OPPOSITE
        // CoordinatesArrayType result; //first trial
        // bool isConverged = master_geometry.ProjectionPointGlobalToLocalSpace(slavePoint, result);
        // localProjection = result[0];

        // return isConverged;
        
        // find closest point in the master brep-------------------------------------------------
        // get the corresponding point on the line

        const double toll = 1e-10;
        double res = toll + 1;
        int it = 0;
        const int itMax = 2e2;

        // starting point NEWTON RAPHSON
        CoordinatesArrayType globCoordMaster;
        CoordinatesArrayType t(3); t[0] = 0; //first trial

        int larger_direction = 1;
        int smaller_direction = 0;
        if (std::abs(normal[0]) > std::abs(normal[1])) { larger_direction = 0; smaller_direction = 1;}

        master_geometry.GlobalCoordinates(globCoordMaster, t);
        const double s_line = (globCoordMaster[larger_direction] - slavePoint[larger_direction])/normal[larger_direction];
        const double smaller_direction_line_current = s_line*normal[smaller_direction] + slavePoint[smaller_direction];
        res = smaller_direction_line_current - globCoordMaster[smaller_direction];

        while (std::abs(res) > toll && it < itMax) {
            //----------------
            std::vector<CoordinatesArrayType> rGlobalSpaceDerivatives;

            master_geometry.GlobalSpaceDerivatives(rGlobalSpaceDerivatives, t, 1);

            // KRATOS_WATCH(rGlobalSpaceDerivatives)


            const double f_der = normal[smaller_direction]/normal[larger_direction] * rGlobalSpaceDerivatives[1][larger_direction] - rGlobalSpaceDerivatives[1][smaller_direction];

            t[0] -= res/f_der;

            master_geometry.GlobalCoordinates(globCoordMaster, t);

            const double s_line = (globCoordMaster[larger_direction] - slavePoint[larger_direction])/normal[larger_direction];

            const double smaller_direction_line_current = s_line*normal[smaller_direction] + slavePoint[smaller_direction];

            res = smaller_direction_line_current - globCoordMaster[smaller_direction];
        }

        localProjection = t[0];

        if (it >= itMax) {KRATOS_WARNING("NR for contact projection has not converged, res:") << res;
        return 0;
        }

        // KRATOS_WATCH(t)
        // KRATOS_WATCH(globCoordMaster)

        distance = norm_2(slavePoint-globCoordMaster);
        return 1;
    }    



    void ComputeProjectionOnSlave() {

        KRATOS_WATCH("ssssssssssssssi")
        // CoordinatesArrayType master_quadrature_point(3);
        // master_quadrature_point[0] = global_space_derivatives_master[0][0]; master_quadrature_point[1] = global_space_derivatives_master[0][1];
        // master_quadrature_point[2] = 0.0;
        // double best_distance = 1e16;
        // CoordinatesArrayType best_projected_point_on_slave_local;
        // int best_bred_id_slave = -1;
        // bool isConverged;
        // for (IndexType i_brep_s = 0; i_brep_s < rIntegrationPoints.size(); i_brep_s++) {
            
        //     CoordinatesArrayType local_coord_projected_on_slave; //first trial
        //     isConverged = rPairedGeometryList[i_brep_s].ProjectionPointGlobalToLocalSpace(master_quadrature_point, local_coord_projected_on_slave, 1e-12);

        //     if (isConverged) {
        //         CoordinatesArrayType point_projected_on_slave(3);
                
        //         std::vector<CoordinatesArrayType> global_space_derivatives_slave(2);
        //         rPairedGeometryList[i_brep_s].GlobalCoordinates(point_projected_on_slave, local_coord_projected_on_slave);

        //         double current_distance = norm_2(master_quadrature_point-point_projected_on_slave);

        //         if (current_distance < best_distance) {
        //             best_distance = current_distance;
        //             best_projected_point_on_slave_local = local_coord_projected_on_slave;
        //             best_bred_id_slave = i_brep_s;
        //         }
        //     } else {
        //         Vector interval; 
        //         rPairedGeometryList[i_brep_s].DomainInterval(interval);
        //         local_coord_projected_on_slave[0] = interval[0];
        //         CoordinatesArrayType point_projected_on_slave(3);        
        //         rPairedGeometryList[i_brep_s].GlobalCoordinates(point_projected_on_slave, local_coord_projected_on_slave);

        //         double current_distance = norm_2(master_quadrature_point-point_projected_on_slave);

        //         if (current_distance < best_distance) {
        //             best_distance = current_distance;
        //             best_projected_point_on_slave_local = local_coord_projected_on_slave;
        //             best_bred_id_slave = i_brep_s;
        //         }

        //         local_coord_projected_on_slave[0] = interval[1];
        //         rPairedGeometryList[i_brep_s].GlobalCoordinates(point_projected_on_slave, local_coord_projected_on_slave);

        //         current_distance = norm_2(master_quadrature_point-point_projected_on_slave);


        //         if (current_distance < best_distance) {
        //             best_distance = current_distance;
        //             best_projected_point_on_slave_local = local_coord_projected_on_slave;
        //             best_bred_id_slave = i_brep_s;
        //         }
        //     }

        // }
        // if (!isConverged) continue;
        // // WARNING ! REMOVE

        // IntegrationPoint<1> integrationPointSlave(best_projected_point_on_slave_local[0]);

        // std::vector<CoordinatesArrayType> global_space_derivatives_slave(2);
        // rPairedGeometryList[best_bred_id_slave].GlobalSpaceDerivatives(
        //         global_space_derivatives_slave,
        //         integrationPointSlave,
        //         1);
    }

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

    GeometriesArrayType mMasterGeometryList;
    GeometriesArrayType mSlaveGeometryList;

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

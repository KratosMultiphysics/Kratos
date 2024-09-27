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

// System includes

// External includes

// Project includes
#include "geometries/nurbs_coupling_geometry_2d.h"
#include "includes/node.h"
#include "utilities/nurbs_utilities/projection_nurbs_contact_utilities.h"

namespace Kratos
{
    ///@name Curve Properties
    ///@{

    /* @brief Provides intersections of the nurbs curve with the knots of the surface,
     *         using the interval of this curve.
     * @param vector of span intervals.
     * @param index of chosen direction, for curves always 0.
     */
    template<class TPointType, class TSurfaceContainerPointType> 
    void NurbsCouplingGeometry2D<TPointType, TSurfaceContainerPointType>::SpansLocalSpaceForParentIntegration(
                                             NurbsSurfaceTypePointer rpNurbsSurfaceParent,
                                             NurbsSurfaceTypePointer rpNurbsSurfacePaired,
                                             std::vector<std::vector<double>>& integration_edges_on_parameter_parent_list,
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
                    bool hasFoundProjection = ProjectionNurbsContactUtilities<TPointType, TSurfaceContainerPointType>::GetProjection(rpNurbsSurfaceParent, rpNurbsSurfacePaired, locCoord, rPairedGeometryList[i_brep_paired], 
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

        // for (int i_brep_parent = 0; i_brep_parent < spans_parent_list.size(); i_brep_parent++) 
        // {
        //     Vector interval; 
        //     rParentGeometryList[i_brep_parent].DomainInterval(interval);
        //     int n_gauss_points = 500;
        //     for (int i = 0; i < n_gauss_points; i++) 
        //     {
        //         double new_gauss_point_loc_pos = interval[0] + (interval[1] - interval[0])/(n_gauss_points-1) * i;
        //         integration_edges_on_parameter_parent_list[i_brep_parent].push_back(new_gauss_point_loc_pos);
        //     }   

        //     std::sort(integration_edges_on_parameter_parent_list[i_brep_parent].begin(), integration_edges_on_parameter_parent_list[i_brep_parent].end());
        //     auto last = std::unique(integration_edges_on_parameter_parent_list[i_brep_parent].begin(), integration_edges_on_parameter_parent_list[i_brep_parent].end());
        //     integration_edges_on_parameter_parent_list[i_brep_parent].erase(last, integration_edges_on_parameter_parent_list[i_brep_parent].end());
        // }
            
    }

    ///@}

    ///@}

    ///@}
    ///@name Integration Points
    ///@{

    /* Creates integration points on the nurbs surface of this geometry.
     * @param return integration points.
     */
    template<class TPointType, class TSurfaceContainerPointType> 
    void NurbsCouplingGeometry2D<TPointType, TSurfaceContainerPointType>::CreateIntegrationPoints(
        NurbsSurfaceTypePointer rpNurbsSurfaceParent,
        NurbsSurfaceTypePointer rpNurbsSurfacePaired,
        std::vector<IntegrationPointsArrayType>& rIntegrationPoints,
        IntegrationInfo& rIntegrationInfo,
        std::vector<std::vector<double>>& spans_parent_list,
        std::vector<std::vector<double>>& spans_paired_list,
        GeometriesArrayType& rParentGeometryList,
        GeometriesArrayType& rPairedGeometryList
        )
    {
        std::vector<std::vector<double>> integration_edges_on_parameter_parent_list(spans_parent_list.size());
        SpansLocalSpaceForParentIntegration(rpNurbsSurfaceParent, rpNurbsSurfacePaired, integration_edges_on_parameter_parent_list, 
                                            spans_parent_list, spans_paired_list, rParentGeometryList, rPairedGeometryList);
        
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
    template<class TPointType, class TSurfaceContainerPointType> 
    void NurbsCouplingGeometry2D<TPointType, TSurfaceContainerPointType>::CreateQuadraturePointGeometries(
        GeometriesArrayType& rResultGeometries,
        IndexType NumberOfShapeFunctionDerivatives,
        IntegrationInfo& rIntegrationInfo)
    {

        std::string name_output_file1 = "txt_files/projection_contact_slave.txt";
        std::string name_output_file2 = "txt_files/projection_contact_master.txt";
        std::ofstream outputFile1(name_output_file1);  outputFile1.close();
        std::ofstream outputFile2(name_output_file2);  outputFile2.close();
        //***************************************************** */

        std::vector<std::vector<double>> spans_slave_list(mSlaveGeometryList.size());
        for (int i = 0; i < mSlaveGeometryList.size(); i++) mSlaveGeometryList[i].SpansLocalSpace(spans_slave_list[i]);

        std::vector<std::vector<double>> spans_master_list(mMasterGeometryList.size());
        for (int i = 0; i < mMasterGeometryList.size(); i++) mMasterGeometryList[i].SpansLocalSpace(spans_master_list[i]);


        // for master
        std::vector<IntegrationPointsArrayType> IntegrationPointsMaster(mMasterGeometryList.size());
        CreateIntegrationPoints(mpNurbsSurfaceMaster, mpNurbsSurfaceSlave, IntegrationPointsMaster, rIntegrationInfo, spans_master_list, spans_slave_list, 
                                mMasterGeometryList, mSlaveGeometryList);


         // // for slave
        std::vector<IntegrationPointsArrayType> IntegrationPointsSlave(mSlaveGeometryList.size());
        CreateIntegrationPoints(mpNurbsSurfaceSlave, mpNurbsSurfaceMaster, IntegrationPointsSlave, rIntegrationInfo, spans_slave_list, spans_master_list, 
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
    template<class TPointType, class TSurfaceContainerPointType> 
    void NurbsCouplingGeometry2D<TPointType, TSurfaceContainerPointType>::CreateQuadraturePointGeometriesOnParent(
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

                // nonzero_control_points_master[0].GetSolutionStepValue(DISPLACEMENT_X);
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
                bool isConvergedAtLeastOnce = false;
                CoordinatesArrayType best_projected_on_slave;
                for (IndexType i_brep_s = 0; i_brep_s < rPairedGeometryList.size(); i_brep_s++) {
                    CoordinatesArrayType local_coord_projected_on_slave; //first trial
                    local_coord_projected_on_slave[0] = 0;
                    CoordinatesArrayType rProjectedPointGlobalCoordinates;
                    double current_distance;
                    // bool isConverged = rPairedGeometryList[i_brep_s].ProjectionPointGlobalToLocalSpace(master_quadrature_point, local_coord_projected_on_slave, 1e-12);
                    
                    bool isConverged = ProjectionNurbsContactUtilities<TPointType, TSurfaceContainerPointType>::NewtonRaphsonCurveOnDeformed(rpNurbsSurfaceParent, rpNurbsSurfacePaired,
                                                    local_coord_projected_on_slave,
                                                    rIntegrationPoints[i_brep_m][i],
                                                    rProjectedPointGlobalCoordinates,
                                                    rParentGeometryList[i_brep_m], rPairedGeometryList[i_brep_s], 
                                                    current_distance,
                                                    10,
                                                    20, 1e-12);

                    // bool isConverged = ProjectionNurbsContactUtilities<TPointType, TSurfaceContainerPointType>:: NewtonRaphsonCurve(
                    //                                             local_coord_projected_on_slave,
                    //                                             master_quadrature_point,
                    //                                             rProjectedPointGlobalCoordinates,
                    //                                             rPairedGeometryList[i_brep_s],
                    //                                             20,
                    //                                             1e-12);
                    if (isConverged) {
                        isConvergedAtLeastOnce = true;
                        // CoordinatesArrayType point_projected_on_slave(3);
                        
                        // std::vector<CoordinatesArrayType> global_space_derivatives_slave(2);
                        // rPairedGeometryList[i_brep_s].GlobalCoordinates(point_projected_on_slave, local_coord_projected_on_slave);

                        // //
                        // best_projected_on_slave = point_projected_on_slave;
                        // //
                        
                        // double current_distance = norm_2(master_quadrature_point-point_projected_on_slave);

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

                        // NEW
                        CoordinatesArrayType displacement_on_projection; CoordinatesArrayType displacement_on_master_quadrature_point;
                        ProjectionNurbsContactUtilities<TPointType, TSurfaceContainerPointType>::GetDisplacement(rpNurbsSurfacePaired, point_projected_on_slave, displacement_on_projection);
                        ProjectionNurbsContactUtilities<TPointType, TSurfaceContainerPointType>::GetDisplacement(rpNurbsSurfaceParent, master_quadrature_point, displacement_on_master_quadrature_point);

                        CoordinatesArrayType point_projected_on_slave_deformed = point_projected_on_slave + displacement_on_projection;
                        CoordinatesArrayType master_quadrature_point_deformed = master_quadrature_point + displacement_on_master_quadrature_point;
                        
                        current_distance = norm_2(master_quadrature_point_deformed-point_projected_on_slave_deformed);

                        // double current_distance = norm_2(master_quadrature_point-point_projected_on_slave);

                        if (current_distance < best_distance) {
                            best_distance = current_distance;
                            best_projected_point_on_slave_local = local_coord_projected_on_slave;
                            best_bred_id_slave = i_brep_s;
                        }

                        local_coord_projected_on_slave[0] = interval[1];
                        rPairedGeometryList[i_brep_s].GlobalCoordinates(point_projected_on_slave, local_coord_projected_on_slave);

                         // NEW
                        ProjectionNurbsContactUtilities<TPointType, TSurfaceContainerPointType>::GetDisplacement(rpNurbsSurfacePaired, point_projected_on_slave, displacement_on_projection);
                        point_projected_on_slave_deformed = point_projected_on_slave + displacement_on_projection;

                        current_distance = norm_2(master_quadrature_point_deformed-point_projected_on_slave_deformed);

                        // current_distance = norm_2(master_quadrature_point-point_projected_on_slave);


                        if (current_distance < best_distance) {
                            best_distance = current_distance;
                            best_projected_point_on_slave_local = local_coord_projected_on_slave;
                            best_bred_id_slave = i_brep_s;
                        }
                    }

                }
                // if (master_quadrature_point[1] < 0.5) {
                //     KRATOS_WATCH(master_quadrature_point)
                //     KRATOS_WATCH(best_projected_on_slave)
                // }
                if (!isConvergedAtLeastOnce) continue;


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


    template<class TPointType, class TSurfaceContainerPointType> 
    void NurbsCouplingGeometry2D<TPointType, TSurfaceContainerPointType>::CreateQuadraturePointGeometriesSlave(
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
   
    template class NurbsCouplingGeometry2D<Kratos::Node, Kratos::PointerVector<Kratos::Node, Kratos::intrusive_ptr<Kratos::Node>, std::vector<Kratos::intrusive_ptr<Kratos::Node>, std::allocator<Kratos::intrusive_ptr<Kratos::Node>>>>>;




///@}
}// namespace Kratos.


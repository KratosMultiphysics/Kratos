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
                                             BrepCurveOnSurfaceArrayType& rParentGeometryList,
                                             BrepCurveOnSurfaceArrayType& rPairedGeometryList)
    {
        for (int i_brep_parent = 0; i_brep_parent < spans_parent_list.size(); i_brep_parent++) {
            integration_edges_on_parameter_parent_list[i_brep_parent] = spans_parent_list[i_brep_parent];
        }

        for (int i_brep_paired = 0; i_brep_paired < spans_paired_list.size(); i_brep_paired++) {
            for (int i = 0; i < spans_paired_list[i_brep_paired].size(); i++) {
                
                CoordinatesArrayType locCoord(3); locCoord[0] = spans_paired_list[i_brep_paired][i]+1e-6;
                
                double incumb_distance = 1e16;
                double projection_parameter_in_best_parent_brep;
                int best_parent_brep_index = -1; //initialize at impossible value
                for (int i_brep_parent = 0; i_brep_parent < spans_parent_list.size(); i_brep_parent++) {
                    Vector interval; 
                    rParentGeometryList[i_brep_parent]->DomainInterval(interval);


                    double current_distance;
                    double current_projection_parameter_in_parent;
                    bool hasFoundProjection = ProjectionNurbsContactUtilities<TPointType, TSurfaceContainerPointType>::
                                                    GetProjection(rpNurbsSurfacePaired, rpNurbsSurfaceParent, locCoord, *rPairedGeometryList[i_brep_paired], 
                                                                  *rParentGeometryList[i_brep_parent], current_projection_parameter_in_parent, current_distance);

                    if (hasFoundProjection && current_distance < incumb_distance &&
                        interval[0] <= current_projection_parameter_in_parent && current_projection_parameter_in_parent <= interval[1]) {
                        incumb_distance = current_distance;
                        best_parent_brep_index = i_brep_parent;
                        projection_parameter_in_best_parent_brep = current_projection_parameter_in_parent;
                    }
                } 
                
                if (best_parent_brep_index > -1) 
                {
                    // check if valuye is already present
                    bool is_present = false;
                    for (IndexType i_vertex = 0; i_vertex < integration_edges_on_parameter_parent_list[best_parent_brep_index].size(); i_vertex++)
                    {
                        if (std::abs(projection_parameter_in_best_parent_brep - integration_edges_on_parameter_parent_list[best_parent_brep_index][i_vertex]) < 1e-4) 
                        {
                            is_present = true;
                            break;
                        }
                    }
                    if (is_present) continue;
                    integration_edges_on_parameter_parent_list[best_parent_brep_index].push_back(projection_parameter_in_best_parent_brep);
                }
            }
        }

        // maybe we can do better
        for (int i_brep_parent = 0; i_brep_parent < spans_parent_list.size(); i_brep_parent++) {
            std::sort(integration_edges_on_parameter_parent_list[i_brep_parent].begin(), integration_edges_on_parameter_parent_list[i_brep_parent].end());
            // auto last = std::unique(integration_edges_on_parameter_parent_list[i_brep_parent].begin(), integration_edges_on_parameter_parent_list[i_brep_parent].end());
            // integration_edges_on_parameter_parent_list[i_brep_parent].erase(last, integration_edges_on_parameter_parent_list[i_brep_parent].end());
        }

        
        // for (int i_brep_parent = 0; i_brep_parent < spans_parent_list.size(); i_brep_parent++) 
        // {
        //     Vector interval; 
        //     rParentGeometryList[i_brep_parent]->DomainInterval(interval);

        //     int n_gauss_points = 5;
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

    template<class TPointType, class TSurfaceContainerPointType> 
    void NurbsCouplingGeometry2D<TPointType, TSurfaceContainerPointType>::FilterSpansForProjection(NurbsSurfaceTypePointer rpNurbsSurfaceParent,
        NurbsSurfaceTypePointer rpNurbsSurfacePaired,
        std::vector<std::vector<double>>& integration_edges_on_parameter_parent_list,
        std::vector<std::vector<std::vector<double>>>& integration_edges_on_parameter_parent_list_filtered,
        BrepCurveOnSurfaceArrayType& rParentGeometryList,
        BrepCurveOnSurfaceArrayType& rPairedGeometryList)
    {
        SizeType count_active_spans = 0;
        bool last_is_converged = false;
        integration_edges_on_parameter_parent_list_filtered.resize(rParentGeometryList.size());
        for (IndexType i_brep_parent = 0; i_brep_parent < rParentGeometryList.size(); i_brep_parent++)
        {
            for (IndexType i = 0; i < integration_edges_on_parameter_parent_list[i_brep_parent].size(); i++)
            {
                CoordinatesArrayType vertex_edge_parent_geometry = ZeroVector(3);
                vertex_edge_parent_geometry[0] = integration_edges_on_parameter_parent_list[i_brep_parent][i];

                CoordinatesArrayType projection_on_slave;
                int best_brep_slave_index;
                bool is_converged = ProjectionNurbsContactUtilities<TPointType, TSurfaceContainerPointType>::GetProjectionOnPairedGeometry(rpNurbsSurfaceParent,
                                                rpNurbsSurfacePaired,
                                                *rParentGeometryList[i_brep_parent],
                                                rPairedGeometryList,
                                                vertex_edge_parent_geometry,
                                                projection_on_slave,
                                                best_brep_slave_index,
                                                25,
                                                50,
                                                1e-9);  


                if ((is_converged && last_is_converged)) //i>0 ||
                {
                    int n_additional_subdivision = 1;
                    double V0 = integration_edges_on_parameter_parent_list[i_brep_parent][i-1];
                    double V1 = integration_edges_on_parameter_parent_list[i_brep_parent][i];
                    double original_span_length = V1-V0;

                    for (IndexType i_div = 0; i_div < n_additional_subdivision; i_div++) {
                        std::vector<double> current_span(2);
                        current_span[0] = V0 + original_span_length/n_additional_subdivision*(i_div); 
                        current_span[1] = V0 + original_span_length/n_additional_subdivision*(i_div+1); 
                        integration_edges_on_parameter_parent_list_filtered[i_brep_parent].push_back(current_span);
                    }
                    // std::vector<double> current_span(2);
                    // current_span[0] = integration_edges_on_parameter_parent_list[i_brep_parent][i-1]; 
                    // current_span[1] = integration_edges_on_parameter_parent_list[i_brep_parent][i]; 
                    // integration_edges_on_parameter_parent_list_filtered[i_brep_parent].push_back(current_span);

                }
                last_is_converged = is_converged;
            }
            
        }
        // KRATOS_WATCH(integration_edges_on_parameter_parent_list)
        // KRATOS_WATCH(integration_edges_on_parameter_parent_list_filtered)

        // exit(0);
        
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
        std::vector<std::vector<IntegrationPointsArrayType>>& rIntegrationPoints,
        IntegrationInfo& rIntegrationInfo,
        std::vector<std::vector<double>>& spans_parent_list,
        std::vector<std::vector<double>>& spans_paired_list,
        BrepCurveOnSurfaceArrayType& rParentGeometryList,
        BrepCurveOnSurfaceArrayType& rPairedGeometryList,
        std::vector<std::vector<int>>& projected_paired_geometry_id
        )
    {   

        std::vector<std::vector<double>> integration_edges_on_parameter_parent_list(spans_parent_list.size());
        SpansLocalSpaceForParentIntegration(rpNurbsSurfaceParent, rpNurbsSurfacePaired, integration_edges_on_parameter_parent_list, 
                                            spans_parent_list, spans_paired_list, rParentGeometryList, rPairedGeometryList);

        std::vector<std::vector<std::vector<double>>> integration_edges_on_parameter_parent_list_filtered;
        FilterSpansForProjection(rpNurbsSurfaceParent, rpNurbsSurfacePaired, integration_edges_on_parameter_parent_list,
                                 integration_edges_on_parameter_parent_list_filtered, rParentGeometryList, rPairedGeometryList);

        projected_paired_geometry_id.resize(spans_parent_list.size());
        
        for (int i_brep_parent = 0; i_brep_parent < spans_parent_list.size(); i_brep_parent++){
            // IntegrationPointUtilities::CreateIntegrationPoints1D(
            //     rIntegrationPoints[i_brep_parent], integration_edges_on_parameter_parent_list[i_brep_parent], rIntegrationInfo);
            rIntegrationPoints[i_brep_parent].resize(2);
            CreateIntegrationPoints1DGauss(
                rIntegrationPoints[i_brep_parent][0], integration_edges_on_parameter_parent_list_filtered[i_brep_parent], rIntegrationInfo);

            // project the integration points and check if everything went correctly
            int n_integration_points = rIntegrationPoints[i_brep_parent][0].size();

            // initialize slave integration points array
            rIntegrationPoints[i_brep_parent][1].resize(n_integration_points);
            projected_paired_geometry_id[i_brep_parent].resize(n_integration_points);

            for (IndexType i = 0; i < rIntegrationPoints[i_brep_parent][0].size(); ++i)
            { 
                // ###########################################################################################################3
                // ###########################################################################################################3

                CoordinatesArrayType projection_on_slave;
                int best_brep_slave_index;
                bool is_converged = ProjectionNurbsContactUtilities<TPointType, TSurfaceContainerPointType>::GetProjectionOnPairedGeometry(rpNurbsSurfaceParent,
                                              rpNurbsSurfacePaired,
                                              *rParentGeometryList[i_brep_parent],
                                              rPairedGeometryList,
                                              rIntegrationPoints[i_brep_parent][0][i],
                                              projection_on_slave,
                                              best_brep_slave_index
                                              ); 
                
                projected_paired_geometry_id[i_brep_parent][i] = best_brep_slave_index;

                IntegrationPoint<1> integrationPointSlave(projection_on_slave[0]);
                rIntegrationPoints[i_brep_parent][1][i] = integrationPointSlave;
                // ###########################################################################################################3
            }
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

        std::ofstream outputFile3("txt_files/Gauss_Point_Contact_coordinates.txt"); outputFile3.close();
        //***************************************************** */

        std::vector<std::vector<double>> spans_slave_list(mSlaveGeometryList.size());
        for (int i = 0; i < mSlaveGeometryList.size(); i++) mSlaveGeometryList[i]->SpansLocalSpace(spans_slave_list[i]);

        std::vector<std::vector<double>> spans_master_list(mMasterGeometryList.size());
        for (int i = 0; i < mMasterGeometryList.size(); i++) mMasterGeometryList[i]->SpansLocalSpace(spans_master_list[i]);


        // for master
        std::vector<std::vector<IntegrationPointsArrayType>> IntegrationPoints(mMasterGeometryList.size());
        std::vector<std::vector<int>> projected_paired_geometry_id;
        CreateIntegrationPoints(mpNurbsSurfaceMaster, mpNurbsSurfaceSlave, IntegrationPoints, rIntegrationInfo, spans_master_list, spans_slave_list, 
                                mMasterGeometryList, mSlaveGeometryList, projected_paired_geometry_id);

        // // Resize containers.
        int numberIntegrationPointsMaster = 0; 
        for (int i = 0; i < mMasterGeometryList.size(); i++) numberIntegrationPointsMaster += IntegrationPoints[i][0].size();  

        const int numberIntegrationPoints = numberIntegrationPointsMaster;
        if (rResultGeometries.size() != numberIntegrationPoints)
            rResultGeometries.resize(numberIntegrationPoints);


        SizeType quadraturePointId = 0;

        this->CreateQuadraturePointGeometriesOnParent(
            rResultGeometries,
            NumberOfShapeFunctionDerivatives,
            IntegrationPoints,
            projected_paired_geometry_id,
            rIntegrationInfo,
            mpNurbsSurfaceMaster,
            mpNurbsSurfaceSlave,
            mMasterGeometryList,
            mSlaveGeometryList,
            MasterIndex,
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
        const std::vector<std::vector<IntegrationPointsArrayType>>& rIntegrationPoints,
        const std::vector<std::vector<int>>& projected_paired_geometry_id,
        IntegrationInfo& rIntegrationInfo,
        NurbsSurfaceTypePointer rpNurbsSurfaceParent,
        NurbsSurfaceTypePointer rpNurbsSurfacePaired,
        BrepCurveOnSurfaceArrayType& rParentGeometryList,
        BrepCurveOnSurfaceArrayType& rPairedGeometryList,
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

        bool is_sbm = rpNurbsSurfaceParent->GetValue(IS_SBM) ||  rpNurbsSurfacePaired->GetValue(IS_SBM);
        bool is_brep_internal = false;
        IndexType InternalBrepSpanU;
        IndexType InternalBrepSpanV;

        for (IndexType i_brep_m = 0; i_brep_m < rIntegrationPoints.size(); i_brep_m++) {

            // SBM CHECK FOR INTEGRATION
            if (is_sbm)
            {
                // IF IS_SBM -> check to which axis the brep is alligned
                std::vector<CoordinatesArrayType> first_integration_point(2); // first integration point of the brep in the parameter space
                std::vector<CoordinatesArrayType> last_integration_point(2); // last integration point of the brep in the parameter space
                rParentGeometryList[i_brep_m]->LocalSpaceDerivatives(first_integration_point,rIntegrationPoints[i_brep_m][0][0],1); 
                rParentGeometryList[i_brep_m]->LocalSpaceDerivatives(last_integration_point,rIntegrationPoints[i_brep_m][0][rIntegrationPoints.size()-1],1); 

                // check if the brep is internal (external breps do not need to be computed in a different knot span) 
                is_brep_internal = true; 
                const double tolerance = 1e-14;
                // At this point all the true are boundary GPs
                if (std::abs(first_integration_point[0][0]-rpNurbsSurfaceParent->KnotsU()[0]) < tolerance)
                    is_brep_internal = false;
                if (std::abs(first_integration_point[0][0]-rpNurbsSurfaceParent->KnotsU()[rpNurbsSurfaceParent->KnotsU().size()-1]) < tolerance) 
                    is_brep_internal = false;
                if (std::abs(first_integration_point[0][1]-rpNurbsSurfaceParent->KnotsV()[0]) < tolerance) 
                    is_brep_internal = false;
                if (std::abs(first_integration_point[0][1]-rpNurbsSurfaceParent->KnotsV()[rpNurbsSurfaceParent->KnotsV().size()-1]) < tolerance) 
                    is_brep_internal = false;

                // If is_brep_internal -> the current GP is a surrogate GP otherwise it is an external GP
                if (is_brep_internal) {
                    // brep = knot spans edge -> same SpanU and SpanV for all its integration points
                    InternalBrepSpanU = NurbsUtilities::GetLowerSpan(rpNurbsSurfaceParent->PolynomialDegreeU(), rpNurbsSurfaceParent->KnotsU(), first_integration_point[0][0]);
                    InternalBrepSpanV = NurbsUtilities::GetLowerSpan(rpNurbsSurfaceParent->PolynomialDegreeV(), rpNurbsSurfaceParent->KnotsV(), first_integration_point[0][1]);

                    // Understand if the knot-boundary is along U or V
                    bool gp_is_along_U;
                    if (std::abs(first_integration_point[0][0] - last_integration_point[0][0]) < tolerance) {
                        gp_is_along_U = true;
                    }
                    else {gp_is_along_U = false; }

                    if (gp_is_along_U && first_integration_point[0][1] > last_integration_point[0][1]) {
                        // Is decreasing along y -> Add 1 at InternalBrepSpanU
                        InternalBrepSpanU++;
                    }
                    if (!gp_is_along_U && first_integration_point[0][0] < last_integration_point[0][0]) {
                        // Is increasing along x -> Add 1 at InternalBrepSpanV
                        InternalBrepSpanV++;
                    }
                }
            }
            //--------------------------------------------------
            //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
            for (IndexType i = 0; i < rIntegrationPoints[i_brep_m][0].size(); ++i)
            {  
                // MASTER
                //ÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑ
                std::vector<CoordinatesArrayType> parameter_space_derivatives_master(2);
                std::vector<CoordinatesArrayType> physical_space_derivatives_master(2);

                rParentGeometryList[i_brep_m]->LocalSpaceDerivatives(
                        parameter_space_derivatives_master,
                        rIntegrationPoints[i_brep_m][0][i],
                        1);

                rParentGeometryList[i_brep_m]->GlobalSpaceDerivatives(
                        physical_space_derivatives_master,
                        rIntegrationPoints[i_brep_m][0][i],
                        1);


                
                // std::ofstream outputFile("txt_files/Gauss_Point_Contact_coordinates.txt", std::ios::app);
                // if (!outputFile.is_open())
                // {
                //     std::cerr << "Failed to open the file for writing." << std::endl;
                //     return;
                // }
                // outputFile << std::setprecision(14); // Set precision to 10^-14
                // outputFile << parameter_space_derivatives_master[0][0] << "  " << parameter_space_derivatives_master[0][1]   <<"\n";
                // outputFile.close();

                
                //++++++++++++++++++++++++++++++++++++++++
                if (rpNurbsSurfaceParent->IsRational()) {
                    shape_function_container_master.ComputeNurbsShapeFunctionValues(
                        rpNurbsSurfaceParent->KnotsU(), rpNurbsSurfaceParent->KnotsV(), rpNurbsSurfaceParent->Weights(),
                        parameter_space_derivatives_master[0][0], parameter_space_derivatives_master[0][1]);
                }
                else {
                    if (is_sbm && is_brep_internal) {
                    // SBM case on internal brep 
                    shape_function_container_master.ComputeBSplineShapeFunctionValuesAtSpan(
                        rpNurbsSurfaceParent->KnotsU(),
                        rpNurbsSurfaceParent->KnotsV(),
                        InternalBrepSpanU,
                        InternalBrepSpanV,
                        parameter_space_derivatives_master[0][0],
                        parameter_space_derivatives_master[0][1]);
                    }
                    else {
                        // Trimming case or external brep
                        shape_function_container_master.ComputeBSplineShapeFunctionValues(
                            rpNurbsSurfaceParent->KnotsU(), rpNurbsSurfaceParent->KnotsV(),
                            parameter_space_derivatives_master[0][0], parameter_space_derivatives_master[0][1]);
                    } 
                }

                //***+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
                


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
                    default_method, rIntegrationPoints[i_brep_m][0][i],
                    N_m, shape_function_derivatives_m);
                    


                int projected_brep_id = projected_paired_geometry_id[i_brep_m][i];

                IntegrationPoint<1> integrationPointSlave = rIntegrationPoints[i_brep_m][1][i];

                std::vector<CoordinatesArrayType> parameter_space_derivatives_slave(2);
                rPairedGeometryList[projected_brep_id]->LocalSpaceDerivatives(
                        parameter_space_derivatives_slave,
                        integrationPointSlave,
                        1);



                //***+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
                CoordinatesArrayType physical_integration_point_slave(2);
                rPairedGeometryList[projected_brep_id]->GlobalCoordinates(
                        physical_integration_point_slave,
                        integrationPointSlave);

                Vector displacement_on_slave_integration_point;

                ProjectionNurbsContactUtilities<TPointType, TSurfaceContainerPointType>::
                                            GetDisplacement(rpNurbsSurfacePaired, parameter_space_derivatives_slave[0], displacement_on_slave_integration_point);

                CoordinatesArrayType physical_integration_point_slave_deformed = physical_integration_point_slave + displacement_on_slave_integration_point;

                Vector displacement_on_projection; 
                ProjectionNurbsContactUtilities<TPointType, TSurfaceContainerPointType>::GetDisplacement(rpNurbsSurfaceParent, parameter_space_derivatives_master[0], displacement_on_projection);

                std::ofstream outputFile2("txt_files/Gauss_Point_Contact_coordinates.txt", std::ios::app);
                if (!outputFile2.is_open())
                {
                    std::cerr << "Failed to open the file for writing." << std::endl;
                    return;
                }

                rParentGeometryList[i_brep_m]->LocalSpaceDerivatives(
                        parameter_space_derivatives_master,
                        rIntegrationPoints[i_brep_m][0][i],
                        1);

                rParentGeometryList[i_brep_m]->GlobalSpaceDerivatives(
                        physical_space_derivatives_master,
                        rIntegrationPoints[i_brep_m][0][i],
                        1);
                
                CoordinatesArrayType master_quadrature_point_parameter(3);
                master_quadrature_point_parameter[0] = parameter_space_derivatives_master[0][0]; master_quadrature_point_parameter[1] = parameter_space_derivatives_master[0][1];
                master_quadrature_point_parameter[2] = 0.0;

                CoordinatesArrayType master_quadrature_point_physical(3);
                master_quadrature_point_physical[0] = physical_space_derivatives_master[0][0]; master_quadrature_point_physical[1] = physical_space_derivatives_master[0][1];
                master_quadrature_point_physical[2] = 0.0;

                Vector displacement_on_master_quadrature_point;
                ProjectionNurbsContactUtilities<TPointType, TSurfaceContainerPointType>::GetDisplacement(rpNurbsSurfaceParent, master_quadrature_point_parameter, displacement_on_master_quadrature_point);
                CoordinatesArrayType master_quadrature_point_deformed = master_quadrature_point_physical + displacement_on_master_quadrature_point;         
                outputFile2 << std::setprecision(14); // Set precision to 10^-14
                outputFile2 << master_quadrature_point_deformed[0] << "  " << master_quadrature_point_deformed[1]  <<" " 
                            << physical_integration_point_slave_deformed[0]  << "  " << physical_integration_point_slave_deformed[1]   <<"\n";
                outputFile2.close();

                //++++++++++++++++++++++++++++++++++++++++
                //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                if (rpNurbsSurfacePaired->IsRational()) {
                    shape_function_container_slave.ComputeNurbsShapeFunctionValues(
                        rpNurbsSurfacePaired->KnotsU(), rpNurbsSurfacePaired->KnotsV(), rpNurbsSurfacePaired->Weights(),
                        parameter_space_derivatives_slave[0][0], parameter_space_derivatives_slave[0][1]);
                }
                else { // TO DO: SBM case
                    
                    shape_function_container_slave.ComputeBSplineShapeFunctionValues(
                        rpNurbsSurfacePaired->KnotsU(), rpNurbsSurfacePaired->KnotsV(),
                        parameter_space_derivatives_slave[0][0], parameter_space_derivatives_slave[0][1]);
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
                    parameter_space_derivatives_master[1][0], parameter_space_derivatives_master[1][1], 
                    parameter_space_derivatives_slave[1][0], parameter_space_derivatives_slave[1][1],
                    &(*rParentGeometryList[i_brep_m]), &(*rPairedGeometryList[projected_brep_id]));

                // set value of the background integration domain
                quadraturePointId++;
            }
        }
    }

    ///@}

    ///@name Projection functionalaties
    ///@{
   
    template class NurbsCouplingGeometry2D<Kratos::Node, Kratos::PointerVector<Kratos::Node, Kratos::intrusive_ptr<Kratos::Node>, std::vector<Kratos::intrusive_ptr<Kratos::Node>, std::allocator<Kratos::intrusive_ptr<Kratos::Node>>>>>;




///@}
}// namespace Kratos.

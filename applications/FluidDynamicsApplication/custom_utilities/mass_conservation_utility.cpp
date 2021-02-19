// --- mass_conservation_utility.cpp ---- Fri, May 22, 2020 11:26:17 AM ----
//  Altair Manufacturing Solver
//
//  Author: Simon Wenczowski, ddiez --- Maintained by: ddiez
//  Copyright: Altair Engineering, Inc. 2015 - 2020
// ************************************************************

// System includes

// External includes

// Project includes
#include "containers/model.h"
#include "modified_shape_functions/tetrahedra_3d_4_modified_shape_functions.h"
#include "modified_shape_functions/triangle_2d_3_modified_shape_functions.h"

// Application includes
#include "mass_conservation_utility.h"
// #include "custom_utilities/compute_filled_volume_utility.h"


namespace Kratos
{

/* Public functions *******************************************************/

/// constructor
MassConservationUtility::MassConservationUtility(
    Model& rModel,
    Parameters Settings)
    : mrModelPart(rModel.GetModelPart(Settings["model_part_name"].GetString())) {

    this->ValidateInputAndInitialize(Settings);
}



/// constructor (direct input of settings)
MassConservationUtility::MassConservationUtility(
    ModelPart& rModelPart,
    Parameters Settings)
    : mrModelPart(rModelPart) {

    this->ValidateInputAndInitialize(Settings);
}

void MassConservationUtility::ValidateInputAndInitialize(
    Parameters Settings) {
    const Parameters default_parameters = GetDefaultParameters();
    Settings.ValidateAndAssignDefaults(default_parameters);
    mEchoLevel = Settings["echo_level"].GetInt();
//     mAverageEdge = ComputeFilledVolumeUtility(mrModelPart).ComputeAverageEdge();
 }

/// Initialization function to find the initial volumes and print first lines in the log-file
/// NOTE THAT WE DO THIS IN THE CalculateInitialVolume TO ENSURE THAT THE LEVEL SET IS ALREADY SET.
/// IT MIGHT HAPPEN THAT THE LEVEL SET IS NOT INITIALIZED AT THIS POINT
// std::string MassConservationUtility::Initialize(){
//     double neg_vol = 0.0;
//     double inter_area = 0.0;

//     // ComputeVolumesAndInterface(neg_vol, inter_area );

//     // mInitialNegativeVolume = neg_vol;
//     std::string output_line =   "------ Initial values ----------------- \n";
//     output_line +=              "  negative volume (water) = " + std::to_string(this->mInitialNegativeVolume) + "\n";
//     output_line +=              "------ Time step values --------------- \n";
//     output_line +=              " negative volume simulation (water) \t\terror_relative \t\tvolumen de agua teorico \t\tvolumen de agua teorico dentro del tanque \t\tincremento de volumen teorico \t\ttime  \n";
//     return output_line;
// }
// Esto calcularlo antes de la simulación del paso de tiempo que toca
double MassConservationUtility::CalculateWaterVolume()
{
    double neg_vol=0.0;
    double inter_area=0.0;
    ComputeVolumesAndInterface(neg_vol, inter_area);
    return neg_vol;
}

std::string MassConservationUtility::CalculateInitialVolume()
{
    mInitialNegativeVolume = CalculateWaterVolume();

    std::string output_line =   "------ Initial values ----------------- \n";
    output_line +=              "  negative volume (water) = " + std::to_string(mInitialNegativeVolume) + "\n";
    output_line +=              "------ Time step values --------------- \n";
    output_line +=              " negative volume simulation (water) \t\terror_relative \t\tvolumen de agua teorico \t\tvolumen de agua teorico dentro del tanque \t\tincremento de volumen teorico \t\ttime  \n";
    return output_line;
}

std::string MassConservationUtility::ComputeBalancedVolume(){

    double neg_vol = 0.0;
    double inter_area = 0.0;

    ComputeVolumesAndInterface(neg_vol, inter_area);
    double net_inflow_inlet = ComputeFlowOverBoundary(INLET);
    double net_inflow_outlet = ComputeFlowOverBoundary(OUTLET);
    mInterfaceArea = inter_area;
    const double current_time = mrModelPart.GetProcessInfo()[TIME];
    const double current_dt = mrModelPart.GetProcessInfo()[DELTA_TIME];

    mQNet0 = net_inflow_inlet + net_inflow_outlet;
    mDeltaTheoreticalNegativeVolume += current_dt * mQNet0;

    mTheoreticalNegativeVolume = mInitialNegativeVolume + mDeltaTheoreticalNegativeVolume;

    mVolumeError = neg_vol - mTheoreticalNegativeVolume;

    KRATOS_WATCH(mTheoreticalNegativeVolume)
    KRATOS_WATCH(mVolumeError)

    KRATOS_INFO_IF("MassConservationUtility", mEchoLevel > 0) << "Theoretical Negative Volume: " << mTheoreticalNegativeVolume << std::endl;
    KRATOS_INFO_IF("MassConservationUtility", mEchoLevel > 0) << "Volume error: " << mVolumeError << std::endl;

    // assembly of the log message
    std::string output_line_timestep =   std::to_string(neg_vol) + "\t\t";
    output_line_timestep +=            std::to_string(mVolumeError) + "\t\t";
    output_line_timestep +=              std::to_string( mTheoreticalNegativeVolume) + "\t\t";
    output_line_timestep +=              std::to_string(mInitialNegativeVolume) + "\t\t";
    output_line_timestep +=             std::to_string(mDeltaTheoreticalNegativeVolume) + "\t\t";
    output_line_timestep +=             std::to_string(current_time) + "\n";


    return output_line_timestep;
}



// double MassConservationUtility::ComputeTimeStepForConvectionSign(double& rOrthogonalFlow){
//     const double tol = 1e-14;
//     double time_step_for_convection = 0.0;
//     double corrected_time_step= 0.0;
//     const double current_dt = mrModelPart.GetProcessInfo()[TIME];
//     if ( mVolumeError < 0.0 ){
//         // case: water volume was lost by mistake
//         mFluidVolumeConservation = FluidVolumeConservation::VOLUME_LOST;
//         double water_outflow_over_boundary = rOrthogonalFlow;
//         KRATOS_WATCH(mVolumeError)
//         KRATOS_WATCH(water_outflow_over_boundary)
//         // checking if flow is sufficient (avoid division by 0)
//         if ( std::abs(water_outflow_over_boundary) > tol ){
//             time_step_for_convection = -mVolumeError / water_outflow_over_boundary;
//              KRATOS_WATCH(time_step_for_convection )
//             corrected_time_step= current_dt + std::abs(time_step_for_convection);


//         }
//     }
//     else if ( mVolumeError > 0.0) {
//         // case: water volume was gained by mistake
//         mFluidVolumeConservation = FluidVolumeConservation::VOLUME_GAINED;
//         double water_inflow_over_boundary = rOrthogonalFlow;
//         KRATOS_WATCH(mVolumeError)
//         KRATOS_WATCH(water_inflow_over_boundary)
//         // checking if flow is sufficient (avoid division by 0)
//         if ( std::abs(water_inflow_over_boundary) > tol ){
//             time_step_for_convection =- mVolumeError / water_inflow_over_boundary;
//             KRATOS_WATCH(time_step_for_convection )
//             corrected_time_step= current_dt - std::abs(time_step_for_convection);

//         }
//     }
//     else {
//         // case: Exactly the correct volume of water is present
//         mFluidVolumeConservation = FluidVolumeConservation::EXACT_VOLUME;
//     }
//     KRATOS_WARNING_IF("MassConservationUtility", time_step_for_convection < 0.0) << "A time step smaller than 0.0 was computed." << std::endl;

//     return time_step_for_convection;
// }










double MassConservationUtility::ComputeTimeStepForConvection(double& rOrthogonalFlow)
{
    const double limit_factor = 0.8;
    const double current_dt = mrModelPart.GetProcessInfo()[DELTA_TIME];

    KRATOS_WATCH(mVolumeError)
    KRATOS_WATCH(rOrthogonalFlow)

    const double tol = 1e-14;
    double time_step_for_convection = 0.0;
    if ( mVolumeError < 0.0 ){

        KRATOS_WATCH("We're negative volume")
        // case: water volume was lost by mistake
        mFluidVolumeConservation = FluidVolumeConservation::VOLUME_LOST;
        // checking if flow is sufficient (avoid division by 0)
        if (std::abs(rOrthogonalFlow) > tol) {
            // If negative orthogonal flow (from air to water) negative delta time to go "backwards" and loose mass
            if (rOrthogonalFlow < 0.0) {
                time_step_for_convection = -std::abs(mVolumeError / rOrthogonalFlow);
            // If positive orthogonal flow (from water to air) positive delta time to go "frontwards" and earn mass
            } else {
                time_step_for_convection = std::abs(mVolumeError / rOrthogonalFlow);
            }
        }
    }
    else if ( mVolumeError > 0.0) {
        KRATOS_WATCH("We're positive volume")
        // case: water volume was gained by mistake
        mFluidVolumeConservation = FluidVolumeConservation::VOLUME_GAINED;
        // checking if flow is sufficient (avoid division by 0)
        if ( std::abs(rOrthogonalFlow) > tol ){
            // If negative orthogonal flow (from air to water) positive delta time to go "backwards" and loose mass
            if (rOrthogonalFlow < 0.0) {
                time_step_for_convection = std::abs(mVolumeError / rOrthogonalFlow);
            // If positive orthogonal flow (from water to air) negative delta time to go "frontwards" and earn mass
            } else {
                time_step_for_convection = -std::abs(mVolumeError / rOrthogonalFlow);
            }
        }
    }
    else {
        // case: Exactly the correct volume of water is present
        mFluidVolumeConservation = FluidVolumeConservation::EXACT_VOLUME;
    }

   

    const double limit_dt = limit_factor * current_dt;
    if (std::abs(time_step_for_convection) > limit_dt){
        time_step_for_convection = (time_step_for_convection < 0.0) ? -limit_dt : limit_dt;
    }
    
    KRATOS_WATCH(time_step_for_convection)

    return time_step_for_convection;
}







void MassConservationUtility::RestoreDistanceValues( const Variable<double>& rAuxDistVar ){

    const int number_nodes = mrModelPart.NumberOfNodes();

    #pragma omp parallel for
    for (int i_node = 0; i_node < number_nodes; ++i_node){
        auto it_node = mrModelPart.NodesBegin() + i_node;
        it_node->GetSolutionStepValue(DISTANCE) = it_node->GetValue(rAuxDistVar);
    }
}



void MassConservationUtility::ReCheckTheMassConservation(){
    double neg_vol = 0.0;
    double inter_area = 0.0;
    ComputeVolumesAndInterface(neg_vol, inter_area );

    mVolumeError = neg_vol - mTheoreticalNegativeVolume;

    mInterfaceArea = inter_area;
    KRATOS_INFO_IF("MassConservationUtility", mEchoLevel > 0) << "WATERVOLUME: " << neg_vol<< std::endl;
    KRATOS_INFO_IF("MassConservationUtility", mEchoLevel > 0) << "Volume error after correction: " << mVolumeError << std::endl;
}

void MassConservationUtility::ComputeVolumesAndInterface(double& rNegativeVolume, double& rInterfaceArea ){
    // initalisation (necessary because no reduction for type reference)
    // double pos_vol = 0.0;
    double neg_vol = 0.0;
    double int_area = 0.0;

    #pragma omp parallel for reduction(+:neg_vol, int_area)
    for (int i_elem = 0; i_elem < static_cast<int>(mrModelPart.NumberOfElements()); ++i_elem){
        // iteration over all elements
        const auto it_elem = mrModelPart.ElementsBegin() + i_elem;

        Matrix shape_functions;
        GeometryType::ShapeFunctionsGradientsType shape_derivatives;

        auto& rGeom = it_elem->GetGeometry();
        unsigned int pt_count_pos = 0;
        unsigned int pt_count_neg = 0;

        // instead of using data.isCut()
        for (unsigned int pt = 0; pt < rGeom.Points().size(); pt++){
            if ( rGeom[pt].FastGetSolutionStepValue(DISTANCE) > 0.0 ){
                pt_count_pos++;
            } else {
                pt_count_neg++;
            }
        }

        // if ( pt_count_pos == rGeom.PointsNumber() ){
        //     // all nodes are positive (pointer is necessary to maintain polymorphism of DomainSize())
        //     pos_vol += it_elem->pGetGeometry()->DomainSize();
        // }
        if ( pt_count_neg == rGeom.PointsNumber() ){
            // all nodes are negative (pointer is necessary to maintain polymorphism of DomainSize())
            neg_vol += it_elem->pGetGeometry()->DomainSize();
        }
        else if ( 0 < pt_count_neg && 0 < pt_count_pos ){
            // element is cut by the surface (splitting)
            Kratos::unique_ptr<ModifiedShapeFunctions> p_modified_sh_func = nullptr;
            Vector w_gauss_pos_side(3, 0.0);
            Vector w_gauss_neg_side(3, 0.0);
            Vector w_gauss_interface(3, 0.0);

            Vector Distance( rGeom.PointsNumber(), 0.0 );
            for (unsigned int i = 0; i < rGeom.PointsNumber(); i++){
                // Control mechanism to avoid 0.0 ( is necessary because "distance_modification" possibly not yet executed )
                double& r_dist = rGeom[i].FastGetSolutionStepValue(DISTANCE);
                if (std::abs(r_dist) < 1.0e-12) {
                    const double aux_dist = 1.0e-6* rGeom[i].GetValue(NODAL_H);
                    if (r_dist > 0.0) {
                        #pragma omp critical
                        r_dist = aux_dist;
                    } else {
                        #pragma omp critical
                        r_dist = -aux_dist;
                    }
                }
                Distance[i] = rGeom[i].FastGetSolutionStepValue(DISTANCE);
            }

            if ( rGeom.PointsNumber() == 3 ){ p_modified_sh_func = Kratos::make_unique<Triangle2D3ModifiedShapeFunctions>(it_elem->pGetGeometry(), Distance); }
            else if ( rGeom.PointsNumber() == 4 ){ p_modified_sh_func = Kratos::make_unique<Tetrahedra3D4ModifiedShapeFunctions>(it_elem->pGetGeometry(), Distance); }
            else { KRATOS_ERROR << "The process can not be applied on this kind of element" << std::endl; }

            // Call the positive side modified shape functions calculator (Gauss weights woulb be enough)
            // Object p_modified_sh_func has full knowledge of slit geometry
            // p_modified_sh_func->ComputePositiveSideShapeFunctionsAndGradientsValues(
            //         shape_functions,                    // N
            //         shape_derivatives,                  // DN
            //         w_gauss_pos_side,                   // includes the weights of the GAUSS points (!!!)
            //         GeometryData::GI_GAUSS_1);          // first order Gauss integration (1 point per triangle)

            // for ( unsigned int i = 0; i < w_gauss_pos_side.size(); i++){
            //     pos_vol += w_gauss_pos_side[i];
            // }

            // Call the negative side modified shape functions calculator
            // Object p_modified_sh_func has full knowledge of slit geometry
            p_modified_sh_func->ComputeNegativeSideShapeFunctionsAndGradientsValues(
                    shape_functions,                    // N
                    shape_derivatives,                  // DN
                    w_gauss_neg_side,                   // includes the weights of the GAUSS points (!!!)
                    GeometryData::GI_GAUSS_1);          // first order Gauss integration

            for ( unsigned int i = 0; i < w_gauss_neg_side.size(); i++){
                neg_vol += w_gauss_neg_side[i];
            }

            // Concerning their area, the positive and negative side of the interface are equal
            p_modified_sh_func->ComputeInterfacePositiveSideShapeFunctionsAndGradientsValues(
                    shape_functions,                    // N
                    shape_derivatives,                  // DN
                    w_gauss_interface,                  // includes the weights of the GAUSS points (!!!)
                    GeometryData::GI_GAUSS_1);          // first order Gauss integration

            for ( unsigned int i = 0; i < w_gauss_interface.size(); i++){
                int_area += std::abs( w_gauss_interface[i] );
            }
        }
    }
    // assigning the values to the arguments of type reference
    rNegativeVolume = neg_vol;
    rInterfaceArea = int_area;
}


double MassConservationUtility::OrthogonalFlowIntoAir()
{
    KRATOS_WATCH("HOLA")
    double outflow = 0.0;
    #pragma omp parallel for reduction(+: outflow)
    for (int i_elem = 0; i_elem < static_cast<int>(mrModelPart.NumberOfElements()); ++i_elem){
        const auto it_elem = mrModelPart.ElementsBegin() + i_elem;

        Matrix shape_functions;
        GeometryType::ShapeFunctionsGradientsType shape_derivatives;

        auto r_geom = it_elem->GetGeometry();
        unsigned int pt_count_pos = 0;
        unsigned int pt_count_neg = 0;
        const bool geom_is_cut = IsGeometryCut(r_geom, pt_count_neg, pt_count_pos);
        // KRATOS_WATCH(pt_count_neg)
        // KRATOS_WATCH(pt_count_pos)
        if (geom_is_cut){
            // KRATOS_WATCH("IS INTERSECTED")
            // element is cut by the surface (splitting)
            Vector w_gauss_interface(3, 0.0);

            Vector distance( r_geom.PointsNumber(), 0.0 );
            for (unsigned int i = 0; i < r_geom.PointsNumber(); i++){
                distance[i] = r_geom[i].GetSolutionStepValue(DISTANCE);
            }

            const auto p_modified_sh_func = GetModifiedShapeFunctions(it_elem->pGetGeometry(), distance);

            // Concerning their area, the positive and negative side of the interface are equal
            p_modified_sh_func->ComputeInterfaceNegativeSideShapeFunctionsAndGradientsValues(
                    shape_functions,                    // N
                    shape_derivatives,                  // DN
                    w_gauss_interface,                  // includes the weights of the GAUSS points (!!!)
                    GeometryData::GI_GAUSS_2);          // second order Gauss integration

            // negative side outwards area normal vector values for the Gauss pts. of given quadrature
            std::vector<Vector> normal_vectors;
            p_modified_sh_func->ComputeNegativeSideInterfaceAreaNormals(
                normal_vectors,
                GeometryData::GI_GAUSS_2
            );

            // iteration over all 3 or 6 integration points
            for ( unsigned int i_gauss = 0; i_gauss < w_gauss_interface.size(); i_gauss++)
            {
                const double& r_weight = w_gauss_interface[i_gauss];
                const auto& r_N = row(shape_functions, i_gauss);

                auto& r_normal = normal_vectors[i_gauss];
                if (norm_2( r_normal ) > 1e-15) {
                    r_normal /= norm_2( r_normal );
                }

                array_1d<double,3> interpolated_velocity = ZeroVector(3);
                for (unsigned int n_node = 0; n_node < r_geom.PointsNumber(); n_node++){
                    noalias( interpolated_velocity ) += r_N[n_node] * (r_geom[n_node].FastGetSolutionStepValue(VELOCITY)
                        - r_geom[n_node].FastGetSolutionStepValue(MESH_VELOCITY));
                }

                // KRATOS_WATCH(interpolated_velocity)
                // KRATOS_WATCH(r_normal)
                const double r_orthogonal_flow = inner_prod( r_normal, interpolated_velocity );
                // KRATOS_WATCH(r_orthogonal_flow)
                outflow += r_weight * r_orthogonal_flow;
            }
        }
    }
    const double global_outflow = mrModelPart.GetCommunicator().GetDataCommunicator().SumAll(outflow);
    // KRATOS_WATCH(global_outflow)
    return global_outflow;
}


double MassConservationUtility::ComputeFlowOverBoundary( const Kratos::Flags BoundaryFlag ){

    // Convention: "mass" is considered as "water", meaning the volumes with a negative distance is considered
    double inflow_over_boundary = 0.0;
    const double epsilon = 1.0e-8;

    #pragma omp parallel for reduction(+: inflow_over_boundary)
    for (int i_cond = 0; i_cond < static_cast<int>(mrModelPart.NumberOfConditions()); ++i_cond){

        // iteration over all conditions (pointer to condition)
        const auto p_condition = mrModelPart.ConditionsBegin() + i_cond;

        if ( p_condition->Is( BoundaryFlag ) ){

            auto& r_geometry = p_condition->GetGeometry();
            Vector distance( r_geometry.PointsNumber(), 0.0 );

            unsigned int neg_count = 0;
            unsigned int pos_count = 0;
            const bool is_geometry_cut = IsGeometryCut(r_geometry, neg_count, pos_count);

            // leave the current iteration of the condition is completely on positive side
            if ( pos_count == r_geometry.PointsNumber() ){ continue; }
            array_1d<double, 3> normal;
            this->CalculateNormal3D( normal, r_geometry);
            if( norm_2( normal ) < epsilon ){ continue; }
            else { normal /= norm_2( normal ); }

            // --- the condition is completely on the negative side (3D)
            if ( neg_count == r_geometry.PointsNumber() ){

                const GeometryType::IntegrationPointsArrayType& IntegrationPoints = r_geometry.IntegrationPoints(GeometryData::GI_GAUSS_2);
                const unsigned int num_gauss = IntegrationPoints.size();
                Vector gauss_pts_det_jabobian = ZeroVector(num_gauss);
                r_geometry.DeterminantOfJacobian(gauss_pts_det_jabobian, GeometryData::GI_GAUSS_2);
                const Matrix n_container = r_geometry.ShapeFunctionsValues( GeometryData::GI_GAUSS_2 );

                for (unsigned int i_gauss = 0; i_gauss < num_gauss; i_gauss++){
                    const auto& N = row(n_container, i_gauss);
                    double const wGauss = gauss_pts_det_jabobian[i_gauss] * IntegrationPoints[i_gauss].Weight();
                    array_1d<double,3> interpolated_velocity = ZeroVector(3);
                    for (unsigned int n_node = 0; n_node < r_geometry.PointsNumber(); n_node++){
                        noalias( interpolated_velocity ) += N[n_node] * (r_geometry[n_node].FastGetSolutionStepValue(VELOCITY)
                            - r_geometry[n_node].FastGetSolutionStepValue(MESH_VELOCITY));
                    }
                    inflow_over_boundary -= wGauss * inner_prod( normal, interpolated_velocity );
                }

            // --- the condition is cut
            } else if ( is_geometry_cut ){

                Matrix r_shape_functions;
                GeometryType::ShapeFunctionsGradientsType r_shape_derivatives;
                Vector w_gauss_neg_side;

                for (unsigned int i = 0; i < r_geometry.PointsNumber(); i++) {
                    distance[i] = r_geometry[i].GetSolutionStepValue( DISTANCE );
                }

                // generating an auxiliary Triangle2D3 geometry the "Triangle2D3ModifiedShapeFunctions" can work with
                const auto aux_2D_triangle = GenerateAuxTriangle( r_geometry );
                // passing the auxiliary triangle
                const auto p_modified_sh_func = Kratos::make_unique<Triangle2D3ModifiedShapeFunctions>( aux_2D_triangle, distance);

                p_modified_sh_func->ComputeNegativeSideShapeFunctionsAndGradientsValues(
                    r_shape_functions,                  // N
                    r_shape_derivatives,                // DN
                    w_gauss_neg_side,                   // includes the weights of the GAUSS points (!!!)
                    GeometryData::GI_GAUSS_2);          // second order Gauss integration

                auto i_points = aux_2D_triangle->IntegrationPoints( GeometryData::GI_GAUSS_2 );

                // interating velocity over the negative area of the condition
                for ( unsigned int i_gauss = 0; i_gauss < w_gauss_neg_side.size(); i_gauss++){
                    const array_1d<double,3>& N = row(r_shape_functions, i_gauss);
                    array_1d<double,3> interpolated_velocity = ZeroVector(3);
                    for (unsigned int n_node = 0; n_node < r_geometry.PointsNumber(); n_node++){
                        noalias( interpolated_velocity ) += N[n_node] * r_geometry[n_node].FastGetSolutionStepValue(VELOCITY);
                    }
                    // abs() is necessary because the auxiliary Triangle2D3 geometry could possibly be inverted
                    // the normal still comes from the oiginal triangle
                    inflow_over_boundary -= std::abs( w_gauss_neg_side[i_gauss] ) * inner_prod( normal, interpolated_velocity );
                }
            }
        }
    }
    const double global_inflow_over_boundary = mrModelPart.GetCommunicator().GetDataCommunicator().SumAll(inflow_over_boundary);
    return global_inflow_over_boundary;
}

void MassConservationUtility::CalculateNormal3D(array_1d<double,3>& An, const Geometry<Node<3> >& pGeometry){

    array_1d<double,3> v1,v2;
    v1[0] = pGeometry[1].X() - pGeometry[0].X();
    v1[1] = pGeometry[1].Y() - pGeometry[0].Y();
    v1[2] = pGeometry[1].Z() - pGeometry[0].Z();

    v2[0] = pGeometry[2].X() - pGeometry[0].X();
    v2[1] = pGeometry[2].Y() - pGeometry[0].Y();
    v2[2] = pGeometry[2].Z() - pGeometry[0].Z();

    MathUtils<double>::CrossProduct(An,v1,v2);
    An *= 0.5;
}



/// Function to convert Triangle3D3N into Triangle2D3N which can be handled by the splitting utilitity
Triangle2D3<Node<3>>::Pointer MassConservationUtility::GenerateAuxTriangle( const Geometry<Node<3> >& rGeom ){

    // Generating auxiliary "Triangle2D3" because the original geometry is "Triangle3D3"

    // vectors that will form the rows of the rotation matrix
    array_1d<double,3> vec_u = rGeom[1].Coordinates() - rGeom[0].Coordinates();
    vec_u /= norm_2( vec_u );

    array_1d<double,3> vec_w;
    MathUtils<double>::CrossProduct(vec_w, vec_u, ( rGeom[2].Coordinates() - rGeom[0].Coordinates() ) );
    vec_w /= norm_2( vec_w );

    array_1d<double,3> vec_v;
    MathUtils<double>::CrossProduct(vec_v, vec_u, vec_w );

    // assembly of the rotation matrix
    Matrix rot_mat = ZeroMatrix(3,3);
    for (unsigned int i = 0; i < 3; i++){
        rot_mat(0,i) = vec_u[i];
        rot_mat(1,i) = vec_v[i];
        rot_mat(2,i) = vec_w[i];
    }

    // rotating the original geometry into a position parallel to the x-y plane by applying the rotation matrix
    // the z-coordinates are then equal and will be neglected afterwards (check is performed before)
    array_1d<double,3> coord1_transformed = prod( rot_mat, rGeom[0].Coordinates() );
    array_1d<double,3> coord2_transformed = prod( rot_mat, rGeom[1].Coordinates() );
    array_1d<double,3> coord3_transformed = prod( rot_mat, rGeom[2].Coordinates() );
    KRATOS_DEBUG_ERROR_IF_NOT( std::abs(coord1_transformed[2] - coord2_transformed[2])<1.0e-7 &&
                            std::abs(coord1_transformed[2] - coord3_transformed[2])<1.0e-7 );

    // creating auxiliary nodes based on the transformed position
    Node<3UL>::Pointer node1 = Kratos::make_intrusive<Kratos::Node<3UL>>( mrModelPart.Nodes().size() + 2, coord1_transformed[0], coord1_transformed[1] );
    Node<3UL>::Pointer node2 = Kratos::make_intrusive<Kratos::Node<3UL>>( mrModelPart.Nodes().size() + 3, coord2_transformed[0], coord2_transformed[1] );
    Node<3UL>::Pointer node3 = Kratos::make_intrusive<Kratos::Node<3UL>>( mrModelPart.Nodes().size() + 4, coord3_transformed[0], coord3_transformed[1] );

    // finally creating the desired Triangle2D3 based on the nodes
    Triangle2D3<Node<3>>::Pointer aux_2D_triangle = Kratos::make_shared< Triangle2D3<Node<3> > >( node1, node2, node3 );
    return aux_2D_triangle;
}

const Parameters MassConservationUtility::GetDefaultParameters()
{
    const Parameters default_parameters = Parameters(R"(
    {
        "model_part_name"   : "define_model_part_name",
        "correct_backwards" : true,
        "echo_level"        : 0
    })" );
    return default_parameters;
}

bool MassConservationUtility::IsGeometryCut(
    const GeometryType &rGeometry,
    unsigned int &PtCountNeg,
    unsigned int &PtCountPos)
{
    PtCountNeg = 0;
    PtCountPos = 0;
    for (unsigned int pt = 0; pt < rGeometry.Points().size(); pt++){
        if ( rGeometry[pt].GetSolutionStepValue(DISTANCE) > 0.0 ){
            PtCountPos++;
        } else {
            PtCountNeg++;
        }
    }
    return (PtCountNeg > 0 && PtCountPos > 0);
}

ModifiedShapeFunctions::Pointer MassConservationUtility::GetModifiedShapeFunctions(
    GeometryType::Pointer pGeometry,
    const Vector& rNodalDistances) {

    const GeometryData::KratosGeometryType geometry_type = pGeometry->GetGeometryType();
    switch (geometry_type){
        case GeometryData::KratosGeometryType::Kratos_Triangle2D3:
            return Kratos::make_unique<Triangle2D3ModifiedShapeFunctions>(pGeometry, rNodalDistances);
        case GeometryData::KratosGeometryType::Kratos_Tetrahedra3D4:
            return Kratos::make_unique<Tetrahedra3D4ModifiedShapeFunctions>(pGeometry, rNodalDistances);
        default:
            KRATOS_ERROR << "Asking for a non-implemented modified shape functions geometry.";
    }
}

void MassConservationUtility::RevertVelocityDirection()
{
    const int number_nodes = mrModelPart.NumberOfNodes();
    #pragma omp parallel for
    for (int i = 0; i < number_nodes; ++i) {
        auto inode = mrModelPart.NodesBegin() + i;
        inode->GetSolutionStepValue(VELOCITY) *= -1.0;
        inode->GetSolutionStepValue(VELOCITY, 1) *= -1.0;
    }

}

double MassConservationUtility::GetVolumeError()
{
    return mVolumeError;
}

double MassConservationUtility::GetTheoreticalVolume()
{
    return mTheoreticalNegativeVolume;
}


// double MassConservationUtility::GetDistanceDelta()
// {
//      const int number_nodes = mrModelPart.NumberOfNodes();
//      for (int i = 0; i < number_nodes; ++i) {
//          auto inode = mrModelPart.NodesBegin() + i;
//          inode->GetSolutionStepValue(DISTANCE,1);
//          ;
//          inode->GetSolutionStepValue(AUX_DISTANCE)= inode->GetSolutionStepValue(DISTANCE)-inode->GetSolutionStepValue(DISTANCE,1)

//          mDeltaDistance= inode->GetSolutionStepValue(AUX_DISTANCE)
//         }
//     return mDeltaDistance
//     }

};  // namespace Kratos.
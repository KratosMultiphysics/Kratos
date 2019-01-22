//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Simon Wenczowski
//
//

// System includes

// External includes

// Project includes
#include "utilities/openmp_utils.h"
#include "utilities/body_distance_calculation_utils.h"
#include "utilities/parallel_levelset_distance_calculator.h"
#include "utilities/variable_utils.h"
#include "includes/deprecated_variables.h"

// Application includes
#include "two_fluids_inlet_process.h"
#include "fluid_dynamics_application_variables.h"


namespace Kratos
{

/* Public functions *******************************************************/

TwoFluidsInletProcess::TwoFluidsInletProcess(
    ModelPart& rModelPart,
    Parameters& rParameters,
    Process::Pointer dist_proc )
    : Process(), mrInletModelPart(rModelPart) {

    KRATOS_CHECK_VARIABLE_KEY( DISTANCE )
    KRATOS_CHECK_VARIABLE_KEY( AUX_DISTANCE )
    KRATOS_CHECK_VARIABLE_KEY( NODAL_AREA )

    // checking all parameters that define the interface
    Parameters default_parameters( R"(
    {
        "interface_normal"          : [0.0,1.0,0.0],
        "point_on_interface"        : [0.0,0.25,0.0],
        "inlet_transition_radius"   : 0.05
    }  )" );

    rParameters.ValidateAndAssignDefaults(default_parameters);

    // finding the complete model of the inlet model part
    ModelPart& rRootModelPart = mrInletModelPart.GetRootModelPart();

    // setting the parameters to the private data members of the class
    mInterfaceNormal = rParameters["interface_normal"].GetVector();
    mInterfacePoint = rParameters["point_on_interface"].GetVector();
    mInletRadius = rParameters["inlet_transition_radius"].GetDouble();

    // normalization of itnerface normal vector
    if ( norm_2( mInterfaceNormal ) > 1.0e-7 ){
        mInterfaceNormal = mInterfacePoint / norm_2( mInterfaceNormal );
    } else {
        KRATOS_ERROR << "Error thrown in TwoFluidsInletProcess: 'interface_normal' in 'interface_settings' must not have a norm of 0.0." << std::endl;
    }
    // checking to avoid 0 transition radius (necessary for division)
    if ( mInletRadius < 1.0e-7 ){
        KRATOS_ERROR << "Error thrown in TwoFluidsInletProcess: 'inlet_transition_radius' in 'interface_settings' must not have a value smaller or equal 0.0." << std::endl;
    }

    // setting flags for the inlet on nodes and conditions
    #pragma omp parallel for
    for (int i_node = 0; i_node < static_cast<int>( mrInletModelPart.NumberOfNodes() ); ++i_node){
        // iteration over all nodes
        auto it_nodes = mrInletModelPart.NodesBegin() + i_node;
        it_nodes->Set( INLET, true );
    }
    #pragma omp parallel for
    for (int i_cond = 0; i_cond < static_cast<int>( mrInletModelPart.NumberOfConditions() ); ++i_cond){
        // iteration over all conditions
        auto it_cond = mrInletModelPart.ConditionsBegin() + i_cond;
        it_cond->Set( INLET, true );
    }

    // Comment: The historical DISTANCE variable is used to compute a distance field that is then stored in a non-historical variable AUX_DISTANCE
    //          The functions for the distance computation do not work on non-historical variables - this is reason for the following procedure (*)

    // (*) temporally storing the distance field as an older version of itself (it can be assured that nothing is over-written at the start)
    #pragma omp parallel for
    for (int i_node = 0; i_node < static_cast<int>( rRootModelPart.NumberOfNodes() ); ++i_node){
        // iteration over all nodes
        auto it_node = rRootModelPart.NodesBegin() + i_node;
        it_node->GetSolutionStepValue(DISTANCE, 2) = it_node->GetSolutionStepValue(DISTANCE, 0);
    }

    // Preparing the distance computation
    // VariableUtils var_utils;
    // var_utils.SetNonHistoricalVariable( IS_VISITED, 0.0, rRootModelPart.Nodes() );
    // var_utils.SetNonHistoricalVariable( IS_VISITED, 1.0, mrInletModelPart.Nodes() );
    // var_utils.SetScalarVar( DISTANCE, 0.0, rRootModelPart.Nodes() );
    // var_utils.SetScalarVar( DISTANCE, 0.0, mrInletModelPart.Nodes() );

    // var_utils.SetScalarVar( NODAL_AREA, 0.0, rRootModelPart.Nodes() );
    // var_utils.SetScalarVar( NODAL_AREA, 1.0, mrInletModelPart.Nodes() );

    // for (int i_node = 0; i_node < static_cast<int>( rRootModelPart.NumberOfNodes() ); ++i_node){
    //     // iteration over all nodes to ensure the existence of the non-hist. variable (OMP parallel leads to problems)
    //     auto it_node = rRootModelPart.NodesBegin() + i_node;
    //     it_node->GetValue(DISTANCE) = 0.0;
    // }

    // const unsigned int dim = rRootModelPart.GetProcessInfo()[DOMAIN_SIZE];
    // // MPI Version used here
    // if ( dim == 2 ){
    //     ParallelDistanceCalculator<2> parallel_distance_util;
    //     parallel_distance_util.CalculateDistancesLagrangianSurface(rRootModelPart, DISTANCE, NODAL_AREA, 100, mInletRadius);
    // } else if ( dim == 3 ){
    //     ParallelDistanceCalculator<3> parallel_distance_util;
    //     parallel_distance_util.CalculateDistancesLagrangianSurface(rRootModelPart, DISTANCE, NODAL_AREA, 100, mInletRadius);
    // } else {
    //     KRATOS_ERROR << "Error thrown in TwoFluidsInletProcess: Dimension not valid." << std::endl;
    // }

    // Preparation for variational distance calculation process
    // setting distance of inlet nodes to 0.0
    // setting rest of distances to 1.0
    #pragma omp parallel for
    for (int i_node = 0; i_node < static_cast<int>( rRootModelPart.NumberOfNodes() ); ++i_node){
        // iteration over all nodes
        auto it_node = rRootModelPart.NodesBegin() + i_node;
        it_node->GetSolutionStepValue(DISTANCE, 0) = 1.0;
    }

    #pragma omp parallel for
    for (int i_node = 0; i_node < static_cast<int>( mrInletModelPart.NumberOfNodes() ); ++i_node){
        // iteration over all nodes
        auto it_node = mrInletModelPart.NodesBegin() + i_node;
        it_node->GetSolutionStepValue(DISTANCE, 0) = -1.0e-7;
    }

    // Variational distance calculation process is executed to calculate distance from inlet
    KRATOS_WATCH( "DITANCE from the inlet is being calcualted ")
    dist_proc->Execute();
    KRATOS_WATCH( "DITANCE from the inlet is being calcualted end ")

    // scaling the distance values such that 1.0 is reached at the inlet
    const double scaling_factor = 1.0 / mInletRadius;
    #pragma omp parallel for
    for (int i_node = 0; i_node < static_cast<int>( rRootModelPart.NumberOfNodes() ); ++i_node){
        // iteration over all nodes
        auto it_node = rRootModelPart.NodesBegin() + i_node;

        if ( (mInletRadius - it_node->GetSolutionStepValue(DISTANCE, 0)) >= 0 ){
            // inside the transition radius (from 1 to 0)
            it_node->GetSolutionStepValue(DISTANCE, 0) = scaling_factor * (mInletRadius - it_node->GetSolutionStepValue(DISTANCE, 0));
        } else {
            // outside the transition radius
            it_node->GetSolutionStepValue(DISTANCE, 0) = 0.0;
        }
    }

    // saving the value of DISTANCE to the non-historical variable AUX_DISTANCE
    VariableUtils var_utils;
    var_utils.SaveScalarVar( DISTANCE, AUX_DISTANCE, rRootModelPart.Nodes() );

    // (*) restoring the original distance field from its stored version
    #pragma omp parallel for
    for (int i_node = 0; i_node < static_cast<int>( rRootModelPart.NumberOfNodes() ); ++i_node){
        // iteration over all nodes
        auto it_node = rRootModelPart.NodesBegin() + i_node;
        it_node->GetSolutionStepValue(DISTANCE, 0) = it_node->GetSolutionStepValue(DISTANCE, 2);
    }

    // subdividing the inlet into two sub_model_part
    mrInletModelPart.CreateSubModelPart("fluid_1_inlet");
    mrInletModelPart.CreateSubModelPart("fluid_2_inlet");
    ModelPart& rWaterInlet = mrInletModelPart.GetSubModelPart("fluid_1_inlet");
    ModelPart& rAirInlet = mrInletModelPart.GetSubModelPart("fluid_2_inlet");

    // classifying nodes (no OMP parallel possible)
    std::vector<IndexType> index_node_water;
    std::vector<IndexType> index_node_air;
    for (int i_node = 0; i_node < static_cast<int>( mrInletModelPart.NumberOfNodes() ); ++i_node){
        auto it_node = mrInletModelPart.NodesBegin() + i_node;
        const double inlet_dist = ComputeNodalDistanceInInletDistanceField(it_node);
        if (inlet_dist <= 0.0){
            index_node_water.push_back( it_node->GetId() );
        } else {
            index_node_air.push_back( it_node->GetId() );
        }
    }
    rWaterInlet.AddNodes( index_node_water );
    rAirInlet.AddNodes( index_node_air );

    // classifying conditions (no OMP parallel possible)
    std::vector<IndexType> index_cond_water;
    std::vector<IndexType> index_cond_air;
    for (int i_cond = 0; i_cond < static_cast<int>( mrInletModelPart.NumberOfConditions() ); ++i_cond){
        auto it_cond = mrInletModelPart.ConditionsBegin() + i_cond;
        unsigned int pos_counter = 0;
        unsigned int neg_counter = 0;
        for (int i_node = 0; i_node < static_cast<int>(it_cond->GetGeometry().PointsNumber()); i_node++){
            const Node<3>& rNode = (it_cond->GetGeometry())[i_node];
            const double inlet_dist = ComputeNodalDistanceInInletDistanceField( rNode );
            if ( inlet_dist > 0 ){ pos_counter++; }
            if ( inlet_dist <= 0 ){ neg_counter++; }
        }
        // the conditions cut by the interface are neither assigned to both the positive and the negative side
        if( pos_counter > 0 ){
            index_cond_air.push_back( it_cond->GetId() );
        }
        if( neg_counter > 0 ){
            index_cond_water.push_back( it_cond->GetId() );
        }
    }
    rWaterInlet.AddConditions( index_cond_water );
    rAirInlet.AddConditions( index_cond_air );

    rRootModelPart.GetCommunicator().Barrier();
}



void TwoFluidsInletProcess::SmoothDistanceField(){

    ModelPart& rRootModelPart = mrInletModelPart.GetRootModelPart();

    // #pragma omp parallel for
    for (int i_node = 0; i_node < static_cast<int>( rRootModelPart.NumberOfNodes() ); ++i_node){
        // iteration over all nodes
        auto it_node = rRootModelPart.NodesBegin() + i_node;

        // check if node is inside "inlet_transition_radius"
        if ( std::abs(it_node->GetValue(AUX_DISTANCE)) > 1.0e-5 ){
            // finding distance value of the node in the inlet field
            const double inlet_dist = ComputeNodalDistanceInInletDistanceField(it_node);
            // finding distance value for the node in the domain distance field
            const double domain_dist = it_node->FastGetSolutionStepValue( DISTANCE, 0 );
            // introducing a smooth transition based in the distance from the inlet stored in "AUX_DISTANCE"
            const double weighting_factor_inlet_field = it_node->GetValue(AUX_DISTANCE);
            const double weighting_factor_domain_field = 1.0 - weighting_factor_inlet_field;

            const double smoothed_dist = weighting_factor_inlet_field * inlet_dist + weighting_factor_domain_field * domain_dist;
            it_node->FastGetSolutionStepValue( DISTANCE, 0 ) = smoothed_dist;
        }

    }

}


/* Private functions ****************************************************/

double TwoFluidsInletProcess::ComputeNodalDistanceInInletDistanceField( const ModelPart::NodesContainerType::iterator node ){
    const array_1d<double,3> distance = node->Coordinates() - mInterfacePoint;
    const array_1d<double,3>& normal = mInterfaceNormal;
    const double inlet_distance =   distance[0]*normal[0] +
                                    distance[1]*normal[1] +
                                    distance[2]*normal[2];
    return inlet_distance;
}

double TwoFluidsInletProcess::ComputeNodalDistanceInInletDistanceField( const Node<3>& node ){
    const array_1d<double,3> distance = node.Coordinates() - mInterfacePoint;
    const array_1d<double,3>& normal = mInterfaceNormal;
    const double inlet_distance =   distance[0]*normal[0] +
                                    distance[1]*normal[1] +
                                    distance[2]*normal[2];
    return inlet_distance;
}


};  // namespace Kratos.
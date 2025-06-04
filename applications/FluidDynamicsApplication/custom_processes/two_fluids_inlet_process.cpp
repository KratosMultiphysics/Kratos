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
#include "utilities/variable_utils.h"

// Application includes
#include "two_fluids_inlet_process.h"
#include "fluid_dynamics_application_variables.h"


namespace Kratos
{

/* Public functions *******************************************************/

TwoFluidsInletProcess::TwoFluidsInletProcess(
    ModelPart& rModelPart,
    Parameters& rParameters,
    Process::Pointer pDistanceProcess )
    : Process(), mrInletModelPart(rModelPart) {

    // checking all parameters that define the interface
    Parameters default_parameters( R"(
    {
        "interface_normal"          : [0.0,1.0,0.0],
        "point_on_interface"        : [0.0,0.25,0.0],
        "inlet_transition_radius"   : 0.05
    }  )" );

    rParameters.ValidateAndAssignDefaults(default_parameters);

    // finding the complete model of the inlet model part
    ModelPart& r_root_model_part = mrInletModelPart.GetRootModelPart();

    // setting the parameters to the private data members of the class
    mInterfaceNormal = rParameters["interface_normal"].GetVector();
    mInterfacePoint = rParameters["point_on_interface"].GetVector();
    mInletRadius = rParameters["inlet_transition_radius"].GetDouble();

    // normalization of interface normal vector
    if ( norm_2( mInterfaceNormal ) > 1.0e-7 ){
        mInterfaceNormal /= norm_2( mInterfaceNormal );
    } else {
        KRATOS_ERROR << "Error thrown in TwoFluidsInletProcess: 'interface_normal' in 'interface_settings' must not have a norm of 0.0." << std::endl;
    }
    // checking to avoid 0 transition radius (necessary for division)
    if ( mInletRadius < 1.0e-7 ){
        KRATOS_ERROR << "Error thrown in TwoFluidsInletProcess: 'inlet_transition_radius' in 'interface_settings' must not have a value smaller or equal 0.0." << std::endl;
    }

    // setting flags for the inlet on nodes and conditions
    VariableUtils().SetFlag( INLET, true, mrInletModelPart.Nodes() );
    VariableUtils().SetFlag( INLET, true, mrInletModelPart.Conditions() );

    // Comment: The historical DISTANCE variable is used to compute a distance field that is then stored in a non-historical variable AUX_DISTANCE
    //          The functions for the distance computation do not work on non-historical variables - this is reason for the following procedure (*)

    const int buffer_size = r_root_model_part.GetBufferSize();
    KRATOS_ERROR_IF( buffer_size < 2 ) << "TwoFluidsInletProcess: There is no space for an intermediate storage" << std::endl;

    // (*) temporally storing the distance field as an older version of itself (it can be assured that nothing is over-written at the start)
    #pragma omp parallel for
    for (int i_node = 0; i_node < static_cast<int>( r_root_model_part.NumberOfNodes() ); ++i_node){
        // iteration over all nodes
        auto it_node = r_root_model_part.NodesBegin() + i_node;
        it_node->GetSolutionStepValue(DISTANCE, (buffer_size-1) ) = it_node->GetSolutionStepValue(DISTANCE, 0);
    }

    // Preparation for variational distance calculation process
    // setting distance of inlet nodes to 0.0 (bit below)
    // setting rest of distances to 1.0
    VariableUtils().SetVariable( DISTANCE, 1.0, r_root_model_part.Nodes() );
    VariableUtils().SetVariable( DISTANCE, -1.0e-7, mrInletModelPart.Nodes() );

    // Variational distance calculation process is executed to calculate distance from inlet
    pDistanceProcess->Execute();

    // scaling the distance values such that 1.0 is reached at the inlet
    const double scaling_factor = 1.0 / mInletRadius;
    #pragma omp parallel for
    for (int i_node = 0; i_node < static_cast<int>( r_root_model_part.NumberOfNodes() ); ++i_node){
        // iteration over all nodes
        auto it_node = r_root_model_part.NodesBegin() + i_node;
        auto& r_dist = it_node->GetSolutionStepValue(DISTANCE, 0);

        if ( (mInletRadius - r_dist) >= 0 ){
            // inside the transition radius (from 1 to 0)
            r_dist = scaling_factor * (mInletRadius - r_dist);
        } else {
            // outside the transition radius
            r_dist = 0.0;
        }
    }

    // Saving the value of DISTANCE to the non-historical variable AUX_DISTANCE
    VariableUtils().SaveVariable<Variable<double>>( DISTANCE, AUX_DISTANCE, r_root_model_part.Nodes() );

    // (*) restoring the original distance field from its stored version
    #pragma omp parallel for
    for (int i_node = 0; i_node < static_cast<int>( r_root_model_part.NumberOfNodes() ); ++i_node){
        // iteration over all nodes
        auto it_node = r_root_model_part.NodesBegin() + i_node;
        it_node->GetSolutionStepValue(DISTANCE, 0) = it_node->GetSolutionStepValue(DISTANCE, (buffer_size-1) );
    }

    // subdividing the inlet into two sub_model_part
    mrInletModelPart.CreateSubModelPart("fluid_1_inlet");
    mrInletModelPart.CreateSubModelPart("fluid_2_inlet");
    ModelPart& r_fluid_1_inlet = mrInletModelPart.GetSubModelPart("fluid_1_inlet");
    ModelPart& r_fluid_2_inlet = mrInletModelPart.GetSubModelPart("fluid_2_inlet");

    // classifying nodes (no OMP parallel possible)
    std::vector<IndexType> index_node_fluid1;
    std::vector<IndexType> index_node_fluid2;
    for (int i_node = 0; i_node < static_cast<int>( mrInletModelPart.NumberOfNodes() ); ++i_node){
        auto it_node = mrInletModelPart.NodesBegin() + i_node;
        const double inlet_dist = inner_prod( ( it_node->Coordinates() - mInterfacePoint ), mInterfaceNormal );
        if (inlet_dist <= 0.0){
            index_node_fluid1.push_back( it_node->GetId() );
        } else {
            index_node_fluid2.push_back( it_node->GetId() );
        }
    }
    r_fluid_1_inlet.AddNodes( index_node_fluid1 );
    r_fluid_2_inlet.AddNodes( index_node_fluid2 );

    // classifying conditions (no OMP parallel possible)
    std::vector<IndexType> index_cond_fluid1;
    std::vector<IndexType> index_cond_fluid2;
    for (int i_cond = 0; i_cond < static_cast<int>( mrInletModelPart.NumberOfConditions() ); ++i_cond){
        auto it_cond = mrInletModelPart.ConditionsBegin() + i_cond;
        unsigned int pos_counter = 0;
        unsigned int neg_counter = 0;
        for (int i_node = 0; i_node < static_cast<int>(it_cond->GetGeometry().PointsNumber()); i_node++){
            const Node& rNode = (it_cond->GetGeometry())[i_node];
            // const double inlet_dist = ComputeNodalDistanceInInletDistanceField( rNode );
            const double inlet_dist = inner_prod( ( rNode.Coordinates() - mInterfacePoint ), mInterfaceNormal );

            if ( inlet_dist > 0 ){ pos_counter++; }
            if ( inlet_dist <= 0 ){ neg_counter++; }
        }
        // the conditions cut by the interface are neither assigned to both the positive and the negative side
        if( pos_counter > 0 ){
            index_cond_fluid2.push_back( it_cond->GetId() );
        }
        if( neg_counter > 0 ){
            index_cond_fluid1.push_back( it_cond->GetId() );
        }
    }
    r_fluid_1_inlet.AddConditions( index_cond_fluid1 );
    r_fluid_2_inlet.AddConditions( index_cond_fluid2 );

    r_root_model_part.GetCommunicator().GetDataCommunicator().Barrier();
}



void TwoFluidsInletProcess::SmoothDistanceField(){

    ModelPart& r_root_model_part = mrInletModelPart.GetRootModelPart();

    // #pragma omp parallel for
    for (int i_node = 0; i_node < static_cast<int>( r_root_model_part.NumberOfNodes() ); ++i_node){

        // iteration over all nodes
        auto it_node = r_root_model_part.NodesBegin() + i_node;

        // check if node is inside "inlet_transition_radius"
        if ( std::abs(it_node->GetValue(AUX_DISTANCE)) > 1.0e-5 ){

            // finding distance value of the node in the inlet field
            const double& r_inlet_dist = inner_prod( ( it_node->Coordinates() - mInterfacePoint ), mInterfaceNormal );

            // finding distance value for the node in the domain distance field
            double& r_domain_dist = it_node->FastGetSolutionStepValue( DISTANCE, 0 );
            // introducing a smooth transition based in the distance from the inlet stored in "AUX_DISTANCE"
            const double& r_weighting_factor_inlet_field = it_node->GetValue(AUX_DISTANCE);

            // setting the smoothed distance field
            r_domain_dist = r_weighting_factor_inlet_field * r_inlet_dist + ( 1.0 - r_weighting_factor_inlet_field ) * r_domain_dist;
        }

    }

}

};  // namespace Kratos.

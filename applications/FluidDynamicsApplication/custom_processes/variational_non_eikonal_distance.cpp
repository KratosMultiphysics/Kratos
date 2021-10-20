//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    KratosAppGenerator
//
//

// System includes
// See the header file

// External includes
// See the header file

// Project includes
// See the header file

// Application includes
#include "variational_non_eikonal_distance.h"
// See the header file

namespace Kratos
{

/* Public functions *******************************************************/
VariationalNonEikonalDistance::VariationalNonEikonalDistance(
    ModelPart& rModelPart,
    TLinearSolver::Pointer plinear_solver)
    : Process(),
      mrModelPart(rModelPart)
{
    // Member variables initialization
    // Nothing!

    // Generate an auxilary model part and populate it by elements of type MySimpleElement
    CreateAuxModelPart();

    auto p_builder_solver = Kratos::make_shared<ResidualBasedBlockBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver> >(plinear_solver);
    //auto p_builder_solver = Kratos::make_unique<ResidualBasedBlockBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver> >(plinear_solver);

    InitializeSolutionStrategy(plinear_solver, p_builder_solver);

    mpGradientCalculator = Kratos::make_unique<ComputeGradientProcessType>(
        mrModelPart,
        DISTANCE_AUX2,
        DISTANCE_GRADIENT,
        NODAL_AREA,
        false);
}

void VariationalNonEikonalDistance::CreateAuxModelPart()
{
    ModelPart::NodesContainerType& r_nodes = mrModelPart.Nodes();
    ModelPart::ElementsContainerType& r_elems = mrModelPart.Elements();

    Model& current_model = mrModelPart.GetModel();
    if(current_model.HasModelPart( mAuxModelPartName ))
        current_model.DeleteModelPart( mAuxModelPartName );

    // Adding DISTANCE_AUX2 and DISTANCE_GRADIENT_X,Y,Z to the solution variables is not needed if it is already a solution variable of the problem
    mrModelPart.AddNodalSolutionStepVariable(DISTANCE_AUX2);
    mrModelPart.AddNodalSolutionStepVariable(DISTANCE_GRADIENT);

    // Ensure that the nodes have distance as a DOF
    VariableUtils().AddDof<Variable<double> >(DISTANCE_AUX2, mrModelPart);

    // Ensure that the nodes have DISTANCE_AUX2 and DISTANCE_GRADIENT_X,Y,Z as a DOF !!!! NOT NEEDED HERE IF IT IS DONE IN THE PYTHON SCRIPT !!!!
    //#pragma omp parallel for
       for (int k = 0; k < static_cast<int>(mrModelPart.NumberOfNodes()); ++k) {
            auto it_node = mrModelPart.NodesBegin() + k;
            it_node->AddDof(DISTANCE_AUX2);
            /* it_node->AddDof(DISTANCE_GRADIENT_X);
            it_node->AddDof(DISTANCE_GRADIENT_Y);
            it_node->AddDof(DISTANCE_GRADIENT_Z); */
        }

    // Generate AuxModelPart
    ModelPart& r_distance_model_part = current_model.CreateModelPart( mAuxModelPartName );

    Element::Pointer p_distance_element = Kratos::make_intrusive<VariationalNonEikonalDistanceElement>();

    ConnectivityPreserveModeler modeler;
    modeler.GenerateModelPart(mrModelPart, r_distance_model_part, *p_distance_element);

    const double delta_time = mrModelPart.pGetProcessInfo()->GetValue(DELTA_TIME);
    r_distance_model_part.pGetProcessInfo()->SetValue(DELTA_TIME, delta_time);
}

void VariationalNonEikonalDistance::Execute()
{
    KRATOS_TRY;

    const unsigned int NumNodes = mrModelPart.NumberOfNodes();
    const unsigned int NumElements = mrModelPart.NumberOfElements();

    const unsigned int num_nodes = 4;

    Model& current_model = mrModelPart.GetModel();
    auto& r_redistancing_model_part = current_model.GetModelPart( mAuxModelPartName );

    double distance_min = 1.0e10;
    unsigned int node_nearest = 0;
    //#pragma omp parallel for
    for (unsigned int k = 0; k < NumNodes; ++k) {
        auto it_node = mrModelPart.NodesBegin() + k;

        const double distance = it_node->FastGetSolutionStepValue(DISTANCE);

        if (std::abs(distance) < 1.0e-12){
            it_node->Fix(DISTANCE_AUX2);
        } else {
            it_node->Free(DISTANCE_AUX2);
        }

        /* if (abs(distance) < distance_min){
            #pragma omp critical
            distance_min = abs(distance);
            #pragma omp critical
            node_nearest = k;
        } */
    }

    double h_min = 1.0e12;
    #pragma omp parallel for
    for(unsigned int i = 0; i < mrModelPart.Elements().size(); i++){
        auto it_elem = mrModelPart.ElementsBegin() + i;
        auto& geom = it_elem->GetGeometry();

        const double he = ElementSizeCalculator<3,4>::AverageElementSize(geom);
        if (he < h_min){
            #pragma omp critical
            h_min = he;
        }

        unsigned int n_pos = 0;
        for (int geom_node=0; geom_node<num_nodes; ++geom_node){
            auto& distance_i = geom[geom_node].FastGetSolutionStepValue(DISTANCE);
            if (distance_i > 0){
                n_pos++;
            }
        }
        /* if (n_pos < num_nodes && n_pos > 0){
            for (int geom_node=0; geom_node<num_nodes; ++geom_node){
                #pragma omp critical
                geom[geom_node].Fix(DISTANCE_AUX2);
            }
        } */
    }

    #pragma omp parallel for
    for (unsigned int k = 0; k < NumNodes; ++k) {
        auto it_node = mrModelPart.NodesBegin() + k;
        it_node->SetValue(NODAL_AREA, 0.0);

        const double distance = it_node->FastGetSolutionStepValue(DISTANCE);

        auto it_node_redistancing = r_redistancing_model_part.NodesBegin() + k;

        if (!it_node->IsFixed(DISTANCE_AUX2)){
            if (distance > 0.0){
                it_node->FastGetSolutionStepValue(DISTANCE_AUX2) = 1.0*h_min;
            } else{
                it_node->FastGetSolutionStepValue(DISTANCE_AUX2) = -1.0*h_min;
            }
        }

        it_node->Set(BOUNDARY,false);
    }

    //KRATOS_WATCH(node_nearest)
    //KRATOS_WATCH(distance_min)

    //auto it_node = mrModelPart.NodesBegin() + node_nearest;
    //it_node->Fix(DISTANCE_AUX2);

    //#pragma omp parallel for
    for (auto it_cond = mrModelPart.ConditionsBegin(); it_cond != mrModelPart.ConditionsEnd(); ++it_cond){
        Geometry< Node<3> >& geom = it_cond->GetGeometry();
        for(unsigned int i=0; i<geom.size(); i++){
            geom[i].Set(BOUNDARY,true);
        }
    }

    r_redistancing_model_part.pGetProcessInfo()->SetValue(FRACTIONAL_STEP,1);

    mpGradientCalculator->Execute(); // To provide the initial condition for DISTANCE_GRADIENT

    KRATOS_INFO("VariationalNonEikonalDistance") << "Reconstruction of levelset, about to solve the LSE" << std::endl;
    mp_solving_strategy->Solve();
    KRATOS_INFO("VariationalNonEikonalDistance") << "Reconstruction of levelset, LSE is solved" << std::endl;

    r_redistancing_model_part.pGetProcessInfo()->SetValue(FRACTIONAL_STEP,2);

    /* #pragma omp parallel for
    for (unsigned int k = 0; k < NumNodes; ++k) {
        auto it_node = mrModelPart.NodesBegin() + k;
        auto it_node_redistancing = r_redistancing_model_part.NodesBegin() + k;

        if (it_node->IsFixed(DISTANCE)){//(it_node->GetValue(IS_STRUCTURE) == 1.0){//
            it_node_redistancing->Fix(DISTANCE_AUX2);
            it_node->Fix(DISTANCE_AUX2);
            it_node->FastGetSolutionStepValue(FLAG_VARIABLE) = 543;
        } else {
            it_node_redistancing->Free(DISTANCE_AUX2);
            it_node->Free(DISTANCE_AUX2);
            it_node->FastGetSolutionStepValue(FLAG_VARIABLE) = -345;
        }
    } */

    mpGradientCalculator->Execute(); // To provide the initial condition for DISTANCE_GRADIENT

    unsigned int iteration = 0;
    double max_grad_norm_deviation = 1.0e2;
    double norm_grad_norm_deviation = 0.0;
    //for (unsigned iter = 0; iter < 50; ++iter){
    while (max_grad_norm_deviation > 2.0e-1 && iteration < 5){
        KRATOS_INFO("VariationalNonEikonalDistance") << "Redistancing, about to solve the LSE" << std::endl;
        mp_solving_strategy->Solve();
        KRATOS_INFO("VariationalNonEikonalDistance") << "Redistancing, LSE is solved" << std::endl;

        mpGradientCalculator->Execute();

        max_grad_norm_deviation = 0.0;
        norm_grad_norm_deviation = 0.0;
        #pragma omp parallel for
        for (unsigned int k = 0; k < NumNodes; ++k) {
            auto it_node = mrModelPart.NodesBegin() + k;
            const double grad_norm_dev = std::abs(
                    norm_2( it_node->GetValue(DISTANCE_GRADIENT) ) - 1.0);
            if ( grad_norm_dev > max_grad_norm_deviation ){
                #pragma omp critical
                max_grad_norm_deviation = grad_norm_dev;
            }
            norm_grad_norm_deviation += grad_norm_dev;
        }
        iteration++;
        KRATOS_INFO("Maximum deviation in the norm of distance gradient") <<
            norm_grad_norm_deviation/static_cast<double>(NumNodes) << std::endl;
    }
    if (max_grad_norm_deviation > 2.0e-1){
        KRATOS_INFO("VariationalNonEikonalDistance") << "Convergence is not achieved." << std::endl;
    }

    /* #pragma omp parallel for
    for (unsigned int k = 0; k < NumNodes; ++k) {
        auto it_node = mrModelPart.NodesBegin() + k;
        auto it_node_redistancing = r_redistancing_model_part.NodesBegin() + k;

        if (it_node->IsFixed(DISTANCE)){//if (it_node->GetValue(IS_STRUCTURE) == 1.0){//
            it_node_redistancing->Fix(DISTANCE_AUX2);
            it_node->Fix(DISTANCE_AUX2);
            it_node->FastGetSolutionStepValue(FLAG_VARIABLE) = 543;
        } else {
            it_node_redistancing->Free(DISTANCE_AUX2);
            it_node->Free(DISTANCE_AUX2);
        }
    } */

    KRATOS_CATCH("")
}

void VariationalNonEikonalDistance::ExecuteInitialize()
{
    KRATOS_TRY;
    KRATOS_CATCH("");
}

void VariationalNonEikonalDistance::ExecuteBeforeSolutionLoop() {
    this->ExecuteInitializeSolutionStep();
    this->ExecuteFinalizeSolutionStep();
}

void VariationalNonEikonalDistance::ExecuteInitializeSolutionStep() {
//Nothing
}

void VariationalNonEikonalDistance::ExecuteFinalizeSolutionStep() {
//Nothing
}

/* Protected functions ****************************************************/

}; // namespace Kratos

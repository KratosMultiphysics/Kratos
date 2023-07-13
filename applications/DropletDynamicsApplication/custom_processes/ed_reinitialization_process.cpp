//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Mohammad R. Hashemi
//
//

// System includes
// See the header file

// External includes
// See the header file

// Project includes
// See the header file

// Application includes
#include "ed_reinitialization_process.h"
// See the header file

namespace Kratos
{

/* Public functions *******************************************************/
EDReinitializationProcess::EDReinitializationProcess(
    ModelPart& rModelPart,
    TLinearSolver::Pointer plinear_solver)
    : Process(),
      mrModelPart(rModelPart)
{
    mAuxModelPartIsCreated = false;

    mpBlockBuilderSolver = Kratos::make_shared<ResidualBasedBlockBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver> >(plinear_solver);

    mpGradientCalculator = Kratos::make_unique<ComputeGradientProcessType>(
        mrModelPart,
        DISTANCE_AUX2,
        DISTANCE_GRADIENT,
        NODAL_AREA,
        false);
}

void EDReinitializationProcess::CreateAuxModelPart()
{
    ModelPart::NodesContainerType& r_nodes = mrModelPart.Nodes();
    ModelPart::ElementsContainerType& r_elems = mrModelPart.Elements();

    Model& current_model = mrModelPart.GetModel();
    if(current_model.HasModelPart( mAuxModelPartName ))
        current_model.DeleteModelPart( mAuxModelPartName );

    ModelPart* p_aux_model_part = &(current_model.CreateModelPart(mAuxModelPartName));

    // Generate
    p_aux_model_part->Nodes().clear();
    p_aux_model_part->Conditions().clear();
    p_aux_model_part->Elements().clear();

    p_aux_model_part->SetProcessInfo(mrModelPart.pGetProcessInfo());
    p_aux_model_part->SetBufferSize(mrModelPart.GetBufferSize());
    p_aux_model_part->SetProperties(mrModelPart.pProperties());
    p_aux_model_part->Tables() = mrModelPart.Tables();

    // Assigning the nodes to the new model part
    p_aux_model_part->Nodes() = mrModelPart.Nodes();

    // Adding DISTANCE_AUX2 to the solution variables is not needed if it is already a solution variable of the problem
    mrModelPart.AddNodalSolutionStepVariable(DISTANCE_AUX2);
    mrModelPart.AddNodalSolutionStepVariable(DISTANCE_GRADIENT);

    // Ensure that the nodes have distance as a DOF
    VariableUtils().AddDof<Variable<double> >(DISTANCE_AUX2, mrModelPart);

    //#pragma omp parallel for
    for (int k = 0; k < static_cast<int>(mrModelPart.NumberOfNodes()); ++k) {
        auto it_node = mrModelPart.NodesBegin() + k;
        it_node->AddDof(DISTANCE_AUX2);
    }

    p_aux_model_part->Elements().reserve(mrModelPart.NumberOfElements());
    for (auto it_elem = mrModelPart.ElementsBegin(); it_elem != mrModelPart.ElementsEnd(); ++it_elem){
        Element::Pointer p_element = Kratos::make_intrusive< EDReinitializationElement >(
            it_elem->Id(),
            it_elem->pGetGeometry(),
            it_elem->pGetProperties());

        // Assign EXACTLY THE SAME GEOMETRY, so that memory is saved!!
        p_element->pGetGeometry() = it_elem->pGetGeometry();

        p_aux_model_part->Elements().push_back(p_element);
    }

    Communicator::Pointer pComm = mrModelPart.GetCommunicator().Create();
    p_aux_model_part->SetCommunicator(pComm);

    const double delta_time = mrModelPart.pGetProcessInfo()->GetValue(DELTA_TIME);
    p_aux_model_part->pGetProcessInfo()->SetValue(DELTA_TIME, delta_time);
    
    mAuxModelPartIsCreated = true;

}

void EDReinitializationProcess::Execute()
{
    KRATOS_TRY;
    KRATOS_INFO("EllipticDistanceReinitialization") << "Start." << std::endl;
    if(mAuxModelPartIsCreated == false){
        CreateAuxModelPart();
        InitializeSolutionStrategy(mpBlockBuilderSolver);
    }
    KRATOS_INFO("EllipticDistanceReinitialization") << "CreateAuxModelPart." << std::endl;
    const unsigned int NumNodes = mrModelPart.NumberOfNodes();
    const unsigned int NumElements = mrModelPart.NumberOfElements();

    const unsigned int num_nodes = 4;

    Model& current_model = mrModelPart.GetModel();
    auto& r_redistancing_model_part = current_model.GetModelPart( mAuxModelPartName );
    KRATOS_INFO("EllipticDistanceReinitialization") << "CreateAuxModelPart." << std::endl;
    // Finding Min element size
    // ToDo 
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
    }
    KRATOS_INFO("EllipticDistanceReinitialization") << "Minimum element size is " << h_min << std::endl;
    // Make DISTANCE_AUX2 a fixed variable for the approxiately interfacial nodes
    // ToDo 
    double distance_min = 1.0e10;
    unsigned int node_nearest = 0;
    //#pragma omp parallel for
    for (unsigned int k = 0; k < NumNodes; ++k) {
        auto it_node = mrModelPart.NodesBegin() + k;

        const double distance = it_node->FastGetSolutionStepValue(DISTANCE);

        it_node->FastGetSolutionStepValue(DISTANCE_AUX2) = distance;

        if ( std::abs(distance) < 1.0e-12 || (it_node->IsFixed(DISTANCE) && !it_node->Is(INLET)) ){
            it_node->Fix(DISTANCE_AUX2);
        } else {
            it_node->Free(DISTANCE_AUX2);
        }
    }
    KRATOS_INFO("EllipticDistanceReinitialization") << "Make DistanceAux2 fixied on surface." << std::endl;
    #pragma omp parallel for
    for (unsigned int k = 0; k < NumNodes; ++k) {

        auto it_node = mrModelPart.NodesBegin() + k;

        it_node->SetValue(NODAL_AREA, 0.0);

        it_node->Set(BOUNDARY,false);
    }
 
    //#pragma omp parallel for
    for (auto it_cond = mrModelPart.ConditionsBegin(); it_cond != mrModelPart.ConditionsEnd(); ++it_cond){
        Geometry< Node >& geom = it_cond->GetGeometry();
        for(unsigned int i=0; i<geom.size(); i++){
            geom[i].Set(BOUNDARY,true);
        }
    }
   KRATOS_INFO("EllipticDistanceReinitialization") << "Set all conditions boundary." << std::endl;
    //************************************************************************************
    //************************************************************************************
    // Doing redistancing for a few iterations: better DISTANCE_GRADIENT for BC in reinitialization step
    /* r_redistancing_model_part.pGetProcessInfo()->SetValue(FRACTIONAL_STEP,2);

    mpGradientCalculator->Execute(); // To provide the initial condition for DISTANCE_GRADIENT

    unsigned int iteration = 0;
    double max_grad_norm_deviation = 1.0e2;
    double norm_grad_norm_deviation = 0.0;

    const unsigned int time_step = mrModelPart.pGetProcessInfo()->GetValue(STEP);

    while (iteration < 5){

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
        KRATOS_INFO("Deviation in the norm of distance gradient") <<
            norm_grad_norm_deviation/static_cast<double>(NumNodes) << std::endl;
    } */

    //if (time_step % 10 == 0)
    {

    //************************************************************************************
    //************************************************************************************
    // Reinitialization
    r_redistancing_model_part.pGetProcessInfo()->SetValue(FRACTIONAL_STEP,0);
    KRATOS_INFO("EllipticDistanceReinitialization") << "FRACTIONAL_STEP is set to 0." << std::endl;
    mpGradientCalculator->Execute(); // To provide the initial condition for DISTANCE_GRADIENT
    KRATOS_INFO("EllipticDistanceReinitialization") << "Gradient Calculation is done." << std::endl;
    #pragma omp parallel for
    for (unsigned int k = 0; k < NumNodes; ++k) {

        auto it_node = mrModelPart.NodesBegin() + k;

        const double distance = it_node->FastGetSolutionStepValue(DISTANCE);

        if (!it_node->IsFixed(DISTANCE_AUX2)){
            if (distance > 0.0){
                it_node->FastGetSolutionStepValue(DISTANCE_AUX2) = 1.0*h_min;
            } else{
                it_node->FastGetSolutionStepValue(DISTANCE_AUX2) = -1.0*h_min;
            }
        }
    }

    KRATOS_INFO("EllipticDistanceReinitialization") << "Reconstruction of levelset, about to solve the LSE" << std::endl;   
    KRATOS_INFO("EllipticDistanceReinitialization") << mp_solving_strategy->Info() << std::endl; 
    mp_solving_strategy->Solve();
    KRATOS_INFO("EllipticDistanceReinitialization") << "Reconstruction of levelset, LSE is solved" << std::endl;

    //************************************************************************************
    //************************************************************************************
    // Reinitialization: step 2 to better estimate norm(grad_phi) on boundaries (seems not essential)
    /* r_redistancing_model_part.pGetProcessInfo()->SetValue(FRACTIONAL_STEP,1);

    mpGradientCalculator->Execute(); // To provide the better approximation of the BC (new gradient)

    KRATOS_INFO("VariationalNonEikonalDistance") << "Reconstruction of levelset, about to solve the LSE" << std::endl;
    mp_solving_strategy->Solve();
    KRATOS_INFO("VariationalNonEikonalDistance") << "Reconstruction of levelset, LSE is solved" << std::endl; */

    //************************************************************************************
    //************************************************************************************
    // Redistancing
    r_redistancing_model_part.pGetProcessInfo()->SetValue(FRACTIONAL_STEP,2);

    unsigned int iteration = 0;
    double max_grad_norm_deviation = 1.0e2;
    double norm_grad_norm_deviation = 0.0;

    std::ofstream gradient_error_file;
    gradient_error_file.open ("gradient_error_file.txt");

    while (max_grad_norm_deviation > 1.0e-3 && iteration < 15){
        mpGradientCalculator->Execute();
        // Calculating max_grad_norm_dev & norm of grad_norm_dev
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
        norm_grad_norm_deviation /= static_cast<double>(NumNodes);
        KRATOS_INFO("Deviation in the norm of distance gradient") <<
            norm_grad_norm_deviation << std::endl;

        gradient_error_file << iteration << "\t" << norm_grad_norm_deviation << std::endl;

        KRATOS_INFO("EllipticDistanceReinitialization") << "Redistancing, about to solve the LSE" << std::endl;
        mp_solving_strategy->Solve();
        KRATOS_INFO("EllipticDistanceReinitialization") << "Redistancing, LSE is solved" << std::endl;

        iteration++;
    }

    auto NumFreeNodes = NumNodes;
    for (unsigned int k = 0; k < NumNodes; ++k) {
        auto it_node = mrModelPart.NodesBegin() + k;
        const double dist = it_node->FastGetSolutionStepValue(DISTANCE_AUX2);
        if (dist > 3.0*h_min || dist < -3.0*h_min){
            it_node->Fix(DISTANCE_AUX2);
            NumFreeNodes -= 1;
        }
    }

    iteration = 0;
    max_grad_norm_deviation = 1.0e2;
    norm_grad_norm_deviation = 0.0;

    while (max_grad_norm_deviation > 1.0e-2 && iteration < 1){

        mpGradientCalculator->Execute();

        max_grad_norm_deviation = 0.0;
        norm_grad_norm_deviation = 0.0;
        #pragma omp parallel for
        for (unsigned int k = 0; k < NumNodes; ++k) {
            auto it_node = mrModelPart.NodesBegin() + k;
            if (!it_node->IsFixed(DISTANCE_AUX2)){
                const double grad_norm_dev = std::abs(
                    norm_2( it_node->GetValue(DISTANCE_GRADIENT) ) - 1.0);
                if ( grad_norm_dev > max_grad_norm_deviation ){
                    #pragma omp critical
                    max_grad_norm_deviation = grad_norm_dev;
                }
                norm_grad_norm_deviation += grad_norm_dev;
            }
        }
        norm_grad_norm_deviation /= static_cast<double>(NumNodes);
        KRATOS_INFO("Deviation in the norm of distance gradient") <<
            norm_grad_norm_deviation << std::endl;

        gradient_error_file << iteration << "\t" << norm_grad_norm_deviation << std::endl;

        KRATOS_INFO("EllipticDistanceReinitialization") << "Redistancing, about to solve the LSE" << std::endl;
        mp_solving_strategy->Solve();
        KRATOS_INFO("EllipticDistanceReinitialization") << "Redistancing, LSE is solved" << std::endl;

        iteration++;
    }

    gradient_error_file.close();

    if (max_grad_norm_deviation > 1.0e-2){
        KRATOS_INFO("EllipticDistanceReinitialization") << "Convergence is not achieved." << std::endl;
    }

    }

    KRATOS_CATCH("")
}

void EDReinitializationProcess::ExecuteInitialize()
{
    KRATOS_TRY;
    KRATOS_CATCH("");
}

void EDReinitializationProcess::ExecuteBeforeSolutionLoop() {
    this->ExecuteInitializeSolutionStep();
    this->ExecuteFinalizeSolutionStep();
}

void EDReinitializationProcess::ExecuteInitializeSolutionStep() {
//Nothing
}

void EDReinitializationProcess::ExecuteFinalizeSolutionStep() {
//Nothing
}

/* Protected functions ****************************************************/

}; // namespace Kratos

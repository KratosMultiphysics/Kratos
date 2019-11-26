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
#include "surface_smoothing_process.h"
// See the header file

namespace Kratos
{

/* Public functions *******************************************************/
SurfaceSmoothingProcess::SurfaceSmoothingProcess(
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

    InitializeSolutionStrategy(plinear_solver, p_builder_solver);
}

void SurfaceSmoothingProcess::CreateAuxModelPart()
{
    ModelPart::NodesContainerType& r_nodes = mrModelPart.Nodes();
    ModelPart::ElementsContainerType& r_elems = mrModelPart.Elements();

    // Adding DISTANCE to the solution variables is not needed if it is already a solution variable of the problem
    mrModelPart.AddNodalSolutionStepVariable(DISTANCE_AUX);
    
    Model& current_model = mrModelPart.GetModel();
    if(current_model.HasModelPart( mAuxModelPartName ))
        current_model.DeleteModelPart( mAuxModelPartName );

    // Ensure that the nodes have DISTANCE_AUX as a DOF !!!! NOT NEEDED HERE IF IT IS DONE IN THE PYTHON SCRIPT !!!!
    #pragma omp parallel for
       for (int k = 0; k < static_cast<int>(mrModelPart.NumberOfNodes()); ++k) {
            auto it_node = mrModelPart.NodesBegin() + k;
            it_node->AddDof(DISTANCE_AUX);
        }

    // Generate AuxModelPart
    ModelPart& r_smoothing_model_part = current_model.CreateModelPart( mAuxModelPartName );

    Element::Pointer p_smoothing_element = Kratos::make_intrusive<SurfaceSmoothingElement>();

    ConnectivityPreserveModeler modeler;
    modeler.GenerateModelPart(mrModelPart, r_smoothing_model_part, *p_smoothing_element);
}

void SurfaceSmoothingProcess::Execute()
{
    #pragma omp parallel for
       for (int k = 0; k < static_cast<int>(mrModelPart.NumberOfNodes()); ++k) {
            auto it_node = mrModelPart.NodesBegin() + k;
            it_node->FastGetSolutionStepValue(DISTANCE_AUX) = it_node->FastGetSolutionStepValue(DISTANCE);
        }

    KRATOS_INFO("SurfaceSmoothingProcess") << "About to solve the LSE" << std::endl;
    mp_solving_strategy->Solve();
    KRATOS_INFO("SurfaceSmoothingProcess") << "LSE is solved" << std::endl;
}

void SurfaceSmoothingProcess::ExecuteInitialize()
{
    KRATOS_TRY;
    KRATOS_CATCH("");
}

void SurfaceSmoothingProcess::ExecuteBeforeSolutionLoop() {
    this->ExecuteInitializeSolutionStep();
    this->ExecuteFinalizeSolutionStep();
}

void SurfaceSmoothingProcess::ExecuteInitializeSolutionStep() {
//Nothing
}

void SurfaceSmoothingProcess::ExecuteFinalizeSolutionStep() {
//Nothing
}

/* Protected functions ****************************************************/

//void SurfaceSmoothingProcess::CheckAndStoreVariablesList(const std::vector<std::string>& rVariableStringArray){
//Nothing
//}

}; // namespace Kratos
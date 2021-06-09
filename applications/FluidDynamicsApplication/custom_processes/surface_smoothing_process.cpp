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

    /* const auto it_element_begin = mrModelPart.ElementsBegin();
    for(unsigned int i = 0; i < rModelPart.Elements().size(); i++){
        auto it_elem = it_element_begin + i;
        const auto& faces = it_elem->pGetGeometry()->Faces();
        const auto& neighbour_elems = it_elem->GetValue(NEIGHBOUR_ELEMENTS);

        for (unsigned int i_face = 0; i_face < faces.size(); i_face++) {
        //if (neighbour_elems[ i_face ].Id() == this->Id() ){
            const auto& r_face = faces[i_face];
            unsigned int contact_node = 0;

            const unsigned int num_face_nodes = 4 - 1;
            for (unsigned int j=0; j < num_face_nodes; ++j){
                if ( r_face[j].GetValue(IS_STRUCTURE) == 1.0 ){
                    contact_node++;
                }
            }

            if (contact_node == num_face_nodes){
                for (unsigned int k=0; k < it_elem->GetGeometry().size(); ++k){
                    KRATOS_INFO("Smoothing Process") << "nodes IS_STRUCTURE "
                        << it_elem->GetGeometry()[k].GetValue(IS_STRUCTURE) << std::endl;
                }
                for (unsigned int j=0; j < num_face_nodes; ++j){
                    KRATOS_INFO("Smoothing Process") << "face nodes Z "
                        << r_face[j].Z() << std::endl;
                }
                KRATOS_INFO("Smoothing Process") << "i_face " << i_face << std::endl;
                for (unsigned int ii = 0; ii < faces.size(); ii++){
                    KRATOS_INFO("Smoothing Process") << it_elem->Id() << " " <<
                        neighbour_elems[ ii ].Id() << std::endl;
                }
            }
        }
    } */

    // Generate an auxilary model part and populate it by elements of type MySimpleElement
    CreateAuxModelPart();

    //auto p_builder_solver = Kratos::make_shared<ResidualBasedBlockBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver> >(plinear_solver);
    //auto p_builder_solver = Kratos::make_unique<ResidualBasedBlockBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver> >(plinear_solver);

    InitializeSolutionStrategy(plinear_solver);//, p_builder_solver);
}

void SurfaceSmoothingProcess::CreateAuxModelPart()
{
    ModelPart::NodesContainerType& r_nodes = mrModelPart.Nodes();
    ModelPart::ElementsContainerType& r_elems = mrModelPart.Elements();
    
    Model& current_model = mrModelPart.GetModel();
    if(current_model.HasModelPart( mAuxModelPartName ))
        current_model.DeleteModelPart( mAuxModelPartName );

    // Adding DISTANCE to the solution variables is not needed if it is already a solution variable of the problem
    mrModelPart.AddNodalSolutionStepVariable(DISTANCE_AUX);

    // Ensure that the nodes have distance as a DOF
    VariableUtils().AddDof<Variable<double> >(DISTANCE_AUX, mrModelPart);

    // Ensure that the nodes have DISTANCE_AUX as a DOF !!!! NOT NEEDED HERE IF IT IS DONE IN THE PYTHON SCRIPT !!!!
    /* #pragma omp parallel for
    for (int k = 0; k < static_cast<int>(mrModelPart.NumberOfNodes()); ++k) {
        auto it_node = mrModelPart.NodesBegin() + k;
        it_node->AddDof(DISTANCE_AUX);
    } */

    // Generate AuxModelPart
    ModelPart& r_smoothing_model_part = current_model.CreateModelPart( mAuxModelPartName );

    Element::Pointer p_smoothing_element = Kratos::make_intrusive<SurfaceSmoothingElement>();

    ConnectivityPreserveModeler modeler;
    modeler.GenerateModelPart(mrModelPart, r_smoothing_model_part, *p_smoothing_element);

    const double delta_time = mrModelPart.pGetProcessInfo()->GetValue(DELTA_TIME);
    r_smoothing_model_part.pGetProcessInfo()->SetValue(DELTA_TIME, delta_time);
        
}

void SurfaceSmoothingProcess::Execute()
{
    KRATOS_TRY;
    
    const unsigned int NumNodes = mrModelPart.NumberOfNodes();
    const unsigned int NumElements = mrModelPart.NumberOfElements();

    std::vector<double> DistDiff(NumNodes, 0.0);
    std::vector<double> DistDiffAvg(NumNodes, 0.0);
    std::vector<double> NumNeighbors(NumNodes, 0.0);

    /* Model& current_model = mrModelPart.GetModel();
    auto& r_smoothing_model_part = current_model.GetModelPart( mAuxModelPartName ); */

    #pragma omp parallel for
    for (unsigned int k = 0; k < NumNodes; ++k) {
        auto it_node = mrModelPart.NodesBegin() + k;
        it_node->SetValue(NODAL_AREA, 0.0);

        //it_node->Free(DISTANCE_AUX);
        const double distance = it_node->FastGetSolutionStepValue(DISTANCE);
        it_node->FastGetSolutionStepValue(DISTANCE_AUX) = distance;

        /* auto it_node_smoothing = r_smoothing_model_part.NodesBegin() + k;
        if (it_node->IsFixed(DISTANCE)){
            it_node_smoothing->Fix(DISTANCE_AUX);
        } else {
            it_node_smoothing->Free(DISTANCE_AUX);
        } */

        if ( it_node->IsFixed(DISTANCE) ) { //GetValue(IS_STRUCTURE) == 1.0 ){
            it_node->Fix(DISTANCE_AUX);
        } else {
            it_node->Free(DISTANCE_AUX);
        }
    }

    //Model& current_model = mrModelPart.GetModel();
    //ModelPart& r_smoothing_model_part = current_model.GetModelPart( mAuxModelPartName );
    //r_smoothing_model_part.pGetProcessInfo()->SetValue(FRACTIONAL_STEP,1);

    KRATOS_INFO("SurfaceSmoothingProcess") << "About to solve the LSE" << std::endl;
    mp_solving_strategy->Solve();
    KRATOS_INFO("SurfaceSmoothingProcess") << "LSE is solved" << std::endl;

    //#pragma omp parallel for
    //for (unsigned int k = 0; k < NumNodes; ++k) {
    //    auto it_node = mrModelPart.NodesBegin() + k;
    //    const double distance = it_node->FastGetSolutionStepValue(DISTANCE_AUX);
    //    it_node->FastGetSolutionStepValue(DISTANCE) = distance;
    //}

    //r_smoothing_model_part.pGetProcessInfo()->SetValue(FRACTIONAL_STEP,2);

    //KRATOS_INFO("SurfaceSmoothingProcess") << "About to solve the LSE" << std::endl;
    //mp_solving_strategy->Solve();
    //KRATOS_INFO("SurfaceSmoothingProcess") << "LSE is solved" << std::endl;

////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////

    #pragma omp parallel for
    for (unsigned int k = 0; k < NumNodes; ++k) {
        auto it_node = mrModelPart.NodesBegin() + k;
        DistDiff[it_node->Id()-1] = it_node->FastGetSolutionStepValue(DISTANCE_AUX) - it_node->FastGetSolutionStepValue(DISTANCE);
    }

    const int num_dim  = 3;
    // const int num_nodes  = num_dim + 1;

    // for (unsigned int k = 0; k < NumElements; ++k) {
    //     auto it_elem = mrModelPart.ElementsBegin() + k;
    //     auto geom = it_elem->pGetGeometry();

    //     for (unsigned int i=0; i < num_nodes; i++){
    //         for (unsigned int j=i+1; j < num_nodes; j++){
    //             const int iId = (*geom)[i].Id() - 1;
    //             const int jId = (*geom)[j].Id() - 1;
    //             NumNeighbors[iId] += 1.0;
    //             NumNeighbors[jId] += 1.0;

    //             DistDiffAvg[iId] += DistDiff[jId];
    //             DistDiffAvg[jId] += DistDiff[iId];
    //         }
    //     }
    // }

    Vector iPosition(num_dim);
    Vector jPosition(num_dim);

    KRATOS_INFO("SurfaceSmoothingProcess") << "Correction For Loop" << std::endl;

    #pragma omp parallel for firstprivate(iPosition, jPosition)
    for (unsigned int k = 0; k < NumNodes; ++k) {
        auto it_node = mrModelPart.NodesBegin() + k;
        const unsigned int iId = it_node->Id() - 1;
        {
            iPosition[0] = it_node->X();
            iPosition[1] = it_node->Y();
            iPosition[2] = it_node->Z();
        }
        auto& n_nodes = it_node->GetValue(NEIGHBOUR_NODES);
        //KRATOS_INFO("SurfaceSmoothingProcess, iId") << iId << std::endl;
        //KRATOS_INFO("SurfaceSmoothingProcess, size") << n_nodes.size() << std::endl;
        for (unsigned int j = 0; j < n_nodes.size(); ++j) {
            if ((n_nodes[j].GetValue(IS_STRUCTURE) == 1.0 && it_node->GetValue(IS_STRUCTURE) == 1.0)
                    || (n_nodes[j].GetValue(IS_STRUCTURE) == 0.0 && it_node->GetValue(IS_STRUCTURE) == 0.0)){
                const unsigned int jId = n_nodes[j].Id() - 1;
                //KRATOS_INFO("SurfaceSmoothingProcess, jId") << jId << std::endl;
                {
                    jPosition[0] = n_nodes[j].X();
                    jPosition[1] = n_nodes[j].Y();
                    jPosition[2] = n_nodes[j].Z();
                }
                const double dist = norm_2(iPosition - jPosition);
                NumNeighbors[iId] += 1.0/dist;
                DistDiffAvg[iId] += DistDiff[jId]/dist;
            }
        }
    }

    #pragma omp parallel for
    for (unsigned int k = 0; k < NumNodes; ++k) {
        auto it_node = mrModelPart.NodesBegin() + k;
        //KRATOS_INFO("SurfaceSmoothingProcess, Nodes Id") << it_node->Id() << std::endl;
        if (NumNeighbors[it_node->Id()-1] != 0.0 && !it_node->IsFixed(DISTANCE)/* && it_node->GetValue(IS_STRUCTURE) == 0.0 */){
             it_node->FastGetSolutionStepValue(DISTANCE_AUX) =
                it_node->FastGetSolutionStepValue(DISTANCE_AUX) - 1.0/NumNeighbors[it_node->Id()-1]*DistDiffAvg[it_node->Id()-1];
        }

    //     it_node->Free(DISTANCE_AUX);
    }

////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////

    // // Auxiliar containers
    // Vector N, values;
    // Matrix InvJ0, J0;
    // double detJ0 = 0.0;

    // // First element iterator
    // const auto it_element_begin = mrModelPart.ElementsBegin();

    // // Current domain size
    // const std::size_t dimension = mrModelPart.GetProcessInfo()[DOMAIN_SIZE];

    // // Initial resize
    // const auto& r_first_element_geometry = it_element_begin->GetGeometry();
    // const std::size_t number_of_nodes_first_element = r_first_element_geometry.PointsNumber();

    // // Iterate over the elements
    // #pragma omp parallel for firstprivate(N, J0, InvJ0, detJ0)
    // for(int i_elem=0; i_elem<static_cast<int>(mrModelPart.Elements().size()); ++i_elem) {
    //     auto it_elem = it_element_begin + i_elem;
    //     auto& r_geometry = it_elem->GetGeometry();

    //     // Current geometry information
    //     const std::size_t number_of_nodes = r_geometry.PointsNumber();

    //     // Resize if needed
    //     if (N.size() != number_of_nodes)
    //         N.resize(number_of_nodes);

    //     // The integration points
    //     const auto& r_integration_method = r_geometry.GetDefaultIntegrationMethod();
    //     const auto& r_integration_points = r_geometry.IntegrationPoints(r_integration_method);
    //     const std::size_t number_of_integration_points = r_integration_points.size();

    //     // The containers of the shape functions and the local gradients
    //     const Matrix& rNcontainer = r_geometry.ShapeFunctionsValues(r_integration_method);

    //     for ( IndexType point_number = 0; point_number < number_of_integration_points; ++point_number ) {
    //         // Getting the shape functions
    //         noalias(N) = row(rNcontainer, point_number);

    //         // Getting the jacobians and local gradients
    //         GeometryUtils::JacobianOnInitialConfiguration(r_geometry, r_integration_points[point_number], J0);
    //         MathUtils<double>::GeneralizedInvertMatrix(J0, InvJ0, detJ0);

    //         const double gauss_point_volume = r_integration_points[point_number].Weight() * detJ0;

    //         double dist_diff = 0.0;
    //         for(std::size_t i_node=0; i_node<number_of_nodes; ++i_node) {
    //             dist_diff += N[i_node] * (r_geometry[i_node].FastGetSolutionStepValue(DISTANCE_AUX)
    //                 - r_geometry[i_node].FastGetSolutionStepValue(DISTANCE));
    //         }

    //         for(std::size_t i_node=0; i_node<number_of_nodes; ++i_node) {
    //             for(std::size_t k=0; k<dimension; ++k) {
    //                 #pragma omp atomic
    //                 DistDiff[ r_geometry[i_node].Id() - 1 ] += N[i_node] * gauss_point_volume * dist_diff;
    //             }

    //             double& vol = r_geometry[i_node].GetValue(NODAL_AREA);

    //             #pragma omp atomic
    //             vol += N[i_node] * gauss_point_volume;
    //         }
    //     }
    // }

    // const auto it_node_begin = mrModelPart.NodesBegin();

    // #pragma omp parallel for
    // for(int i = 0; i < static_cast<int>(mrModelPart.Nodes().size()); ++i) {
    //     auto it_node = it_node_begin + i;
    //     it_node->FastGetSolutionStepValue(DISTANCE_AUX) -=
    //         ( DistDiff[ it_node->Id() - 1 ] / it_node->GetValue(NODAL_AREA) );
    // }

////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////

    KRATOS_CATCH("");
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
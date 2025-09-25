//
//   Project Name:        Kratos
//   Last Modified by:    $Author: gcasas $
//   Date:                $Date: 2016-03-08 08:56:42 $
//
//
#include "utilities/math_utils.h"
#include "derivative_recovery.h"
#include "linear_solvers/amgcl_solver.h"
#include "swimming_DEM_application.h"
#include "swimming_dem_application_variables.h"
#include "custom_utilities/fields/ethier_flow_field.h"
#include "omp.h"

namespace Kratos
{
template <std::size_t TDim>
void DerivativeRecovery<TDim>::RecoverGradientOfAScalar(const VariableData& origin_variable, const VariableData& destination_variable)
{
    for (int i = 0; i < (int)mModelPart.Nodes().size(); ++i){
        NodeIteratorType i_particle = mModelPart.NodesBegin() + i;
        Node::Pointer p_node = *(i_particle.base());
        array_1d<double, 3>& gradient = p_node->FastGetSolutionStepValue(TORQUE);
        gradient[0] = 0.0;
        gradient[1] = 0.0;
        gradient[2] = 99.0;
    }
}
//***************************************************************************************************************
//***************************************************************************************************************
template <std::size_t TDim>
void DerivativeRecovery<TDim>::AddTimeDerivative(ModelPart& r_model_part, Variable<array_1d<double, 3> >& material_derivative_container)
{
    const double delta_time_inv = 1.0 / r_model_part.GetProcessInfo()[DELTA_TIME];
    for (NodeIteratorType inode = r_model_part.NodesBegin(); inode != r_model_part.NodesEnd(); ++inode){
        array_1d <double, 3>& material_derivative = inode->FastGetSolutionStepValue(material_derivative_container);
        const array_1d <double, 3> eulerian_rate_of_change = inode->FastGetSolutionStepValue(ACCELERATION);
        noalias(material_derivative) += eulerian_rate_of_change;
    }
}
//***************************************************************************************************************
//***************************************************************************************************************
template <std::size_t TDim>
void DerivativeRecovery<TDim>::AddTimeDerivativeComponent(ModelPart& r_model_part, Variable<array_1d<double, 3> >& material_derivative_container, const int i_component)
{
    const double delta_time_inv = 1.0 / r_model_part.GetProcessInfo()[DELTA_TIME];
    for (NodeIteratorType inode = r_model_part.NodesBegin(); inode != r_model_part.NodesEnd(); ++inode){
        const double eulerian_rate_of_change = delta_time_inv * (inode->FastGetSolutionStepValue(VELOCITY)[i_component] - inode->FastGetSolutionStepValue(VELOCITY, 1)[i_component]);
        array_1d <double, 3>& material_derivative = inode->FastGetSolutionStepValue(material_derivative_container);
        material_derivative[i_component] += eulerian_rate_of_change;
    }
}
//***************************************************************************************************************
//***************************************************************************************************************
template <std::size_t TDim>
void DerivativeRecovery<TDim>::CalculateLocalMassMatrix(const unsigned& number_of_dofs, const ModelPart::ElementsContainerType::iterator& it_elem, Matrix& local_mass_matrix)
{
    GeometryData::IntegrationMethod integration_method = it_elem->GetIntegrationMethod();
    Geometry<Node>& r_geometry = it_elem->GetGeometry();
    unsigned local_size = r_geometry.size() * number_of_dofs;
    local_mass_matrix.resize(local_size, local_size);

    unsigned int NumNodes = r_geometry.size();
    const std::vector<IntegrationPoint<3>> r_integrations_points = r_geometry.IntegrationPoints( integration_method );
    unsigned int r_number_integration_points = r_geometry.IntegrationPointsNumber(integration_method);
    Matrix NContainer = r_geometry.ShapeFunctionsValues(integration_method);

    Vector detJ_vector(r_number_integration_points);
    r_geometry.DeterminantOfJacobian(detJ_vector, integration_method);

    DenseVector<Matrix> shape_derivatives;
    r_geometry.ShapeFunctionsIntegrationPointsGradients(shape_derivatives, detJ_vector, integration_method);

    for (unsigned int g = 0; g < r_number_integration_points; g++)
    {
        double Weight = r_integrations_points[g].Weight() * detJ_vector[g];
        for (unsigned int a = 0; a < NumNodes; ++a)
        {
            for (unsigned int b = 0; b < NumNodes; ++b)
            {
                for(unsigned int i = 0; i < TDim; i++)
                {
                    unsigned local_pos_a = a * TDim + i;
                    unsigned local_pos_b = b * TDim + i;
                    // for(unsigned int j = 0; j < number_of_dofs; j++)
                    // {
                    // unsigned local_pos_b = b * number_of_dofs + j;
                    local_mass_matrix(local_pos_a, local_pos_b) += Weight * NContainer(g, a) * NContainer(g, b);
                    // }
                }
            }
        }
    }
}
//***************************************************************************************************************
//***************************************************************************************************************
template <std::size_t TDim>
void DerivativeRecovery<TDim>::AssembleMassMatrix(SparseSpaceType::MatrixType& global_matrix, const Matrix& local_lhs, std::map<unsigned, unsigned>& element_id_map)
{
    const SizeType local_size = local_lhs.size1();

    for (unsigned i_local = 0; i_local < local_size; i_local++)
    {
        unsigned i_global = element_id_map[i_local];
        double* values_vector = global_matrix.value_data().begin();

        std::size_t* index1_vector = global_matrix.index1_data().begin();
        std::size_t* index2_vector = global_matrix.index2_data().begin();
        size_t left_limit = index1_vector[i_global];
        //find the first entry
        size_t last_pos;
        unsigned int position = left_limit;
        while(element_id_map[0] != index2_vector[position])position++;

        last_pos = position;
        size_t last_found = element_id_map[0];
        double& r_a = values_vector[last_pos];
        const double& v_a = local_lhs(i_local,0);
        AtomicAdd(r_a,  v_a);
        size_t posit = 0;
        unsigned int pos;
        for (unsigned int j=1; j<element_id_map.size(); j++) {
            unsigned int id_to_find = element_id_map[j];
            if(id_to_find > last_found) {
              pos = last_pos+1;
              while(id_to_find != index2_vector[pos]) pos++;
              posit = pos;
            } else if(id_to_find < last_found) {
              pos = last_pos-1;
              while(id_to_find != index2_vector[pos]) pos--;
              posit = pos;
            } else {
                posit = last_pos;
            }
            double& r = values_vector[posit];
            const double& v = local_lhs(i_local,j);
            AtomicAdd(r,  v);
            last_found = id_to_find;
            last_pos = posit;
        }
    }
}
//***************************************************************************************************************
//***************************************************************************************************************
template <std::size_t TDim>
void DerivativeRecovery<TDim>::ComputeNonZeroMassMatrixIndex(ModelPart& r_model_part, const std::vector<std::map<unsigned, unsigned>>& elements_id_map, std::vector<std::unordered_set<unsigned>>& inds)
{
    unsigned system_size = inds.size();
    unsigned number_of_elements = elements_id_map.size();
    for (unsigned int e = 0; e < number_of_elements; e++)
    {
        ModelPart::ElementsContainerType::iterator it_elem = r_model_part.ElementsBegin() + e;
        Geometry<Node>& r_geometry = it_elem->GetGeometry();
        unsigned NumNodes = r_geometry.size();

        // Get over all pair of nodes of each element
        std::map<unsigned, unsigned> element_inds_map = elements_id_map[e];
        unsigned number_of_dofs = TDim;


        // Node in row
        for(unsigned n = 0; n < NumNodes; n++)
        {
            // Node in col
            for(unsigned m = 0; m < NumNodes; m++)
            {
                for(unsigned d1 = 0; d1 < number_of_dofs; d1++)
                {
                    unsigned local_dof_1_pos  = n * number_of_dofs + d1;
                    for(unsigned d2 = 0; d2 < number_of_dofs; d2++)
                    {
                        unsigned local_dof_2_pos  = m * number_of_dofs + d2;
                        unsigned global_dof_1_pos = element_inds_map[local_dof_1_pos];
                        unsigned global_dof_2_pos = element_inds_map[local_dof_2_pos];
                        inds[global_dof_1_pos].insert(global_dof_2_pos);
                    }
                }
            }
        }
    }
}
//***************************************************************************************************************
//***************************************************************************************************************
template <std::size_t TDim>
void DerivativeRecovery<TDim>::ConstructMassMatrixStructure(ModelPart& r_model_part, const unsigned& system_size, const std::vector<std::map<unsigned, unsigned>>& elements_id_map, SparseSpaceType::MatrixType& global_mass_matrix)
{
    // See `mapping_matrix_utilities.cpp`

    // Indices for which mass matrix is non-zero
    std::vector<std::unordered_set<unsigned>> indices(system_size);
    ComputeNonZeroMassMatrixIndex(r_model_part, elements_id_map, indices);

    // computing the number of non-zero entries
    SizeType num_non_zero_entries = 0;
    for (const auto& r_row_indices : indices) { // looping the indices per row
        num_non_zero_entries += r_row_indices.size(); // adding the number of col-indices
    }
    global_mass_matrix = SparseSpaceType::MatrixType(indices.size(), indices.size(), num_non_zero_entries);

    double* mass_matrix_values = global_mass_matrix.value_data().begin();
    IndexType* mass_matrix_row_indices = global_mass_matrix.index1_data().begin();
    IndexType* mass_matrix_col_indices = global_mass_matrix.index2_data().begin();

    mass_matrix_row_indices[0] = 0;

    // #pragma omp parallel for
    for (IndexType i=0; i < system_size; ++i) {
        mass_matrix_row_indices[i+1] = mass_matrix_row_indices[i] + indices[i].size();  // mass_matrix_row_indices[i + 1] ???
    }

    #pragma omp parallel for
    for (IndexType i = 0; i < system_size; ++i) {

        const IndexType row_begin = mass_matrix_row_indices[i];
        const IndexType row_end = mass_matrix_row_indices[i + 1];
        IndexType k = row_begin;

        std::unordered_set<unsigned>::iterator row_inds_it;
        // for (const auto index : indices[i]) {
        for(row_inds_it = indices[i].begin(); row_inds_it != indices[i].end(); row_inds_it++){
            // mass_matrix_col_indices[k] = index;
            mass_matrix_col_indices[k] = *row_inds_it;
            mass_matrix_values[k] = 0.0;
            ++k;
        }
        std::sort(&mass_matrix_col_indices[row_begin], &mass_matrix_col_indices[row_end]);
    }

    global_mass_matrix.set_filled(indices.size() + 1, num_non_zero_entries);
}
//***************************************************************************************************************
//***************************************************************************************************************
template <std::size_t TDim>
void DerivativeRecovery<TDim>::CalculateFieldL2Projection(ModelPart& r_model_part,
                                                                 VelocityField::Pointer flow_field,
                                                                 Variable<array_1d<double, 3> >& vector_container,
                                                                 Variable<array_1d<double, 3> >& vector_rate_container,
                                                                 Variable<array_1d<double, 3> >& material_derivative_container)
{
    KRATOS_INFO("SwimmingDEM") << "Constructing the material derivative for a given flow field..." << std::endl;

    // System size
    const unsigned number_of_elements = r_model_part.NumberOfElements();
    const unsigned number_of_nodes = r_model_part.NumberOfNodes();
    const unsigned number_of_dofs = TDim;  // vx, vy, vz (vi,j for i, j = 0, ..., TDim - 1)
    const SizeType system_size = number_of_nodes * number_of_dofs;

    // Map the ids of the nodes to their position in the model part
    std::map <std::size_t, unsigned int> id_to_position;
    unsigned int node_pos = 0;
    for (NodeIteratorType inode = r_model_part.NodesBegin(); inode != r_model_part.NodesEnd(); ++inode){
        noalias(inode->FastGetSolutionStepValue(material_derivative_container)) = ZeroVector(3);
        id_to_position[inode->Id()] = node_pos;
        ++node_pos;
    }

    // Compute the mass matrix of the system (only once)
    if(!mMassMatrixAlreadyComputed)
    {
        // Map of local dof position to global dof position
        std::vector<std::map<unsigned, unsigned>> elements_id_map(number_of_elements);
        for (unsigned int e = 0; e < number_of_elements; e++)
        {
            ModelPart::ElementsContainerType::iterator it_elem = r_model_part.ElementsBegin() + e;
            Geometry<Node>& r_geometry = it_elem->GetGeometry();

            std::map<unsigned, unsigned> element_map;
            for(unsigned n = 0; n < r_geometry.size(); n++)
            {
                const SizeType node_pos = id_to_position[r_geometry[n].Id()];
                for (unsigned i = 0; i < number_of_dofs; i++)
                {
                    unsigned local_dof_position = n * number_of_dofs + i;  // Each node has a local position in the element
                    unsigned global_dof_position = node_pos * number_of_dofs + i;  // Each DOF in the model part has a unique position (we assume the nodes in r_geometry are always at the same order)
                    element_map[local_dof_position] = global_dof_position;  // This maps the local DOF position of the element to the global position of this DOF
                }
            }
            elements_id_map[e] = element_map;
        }

        // Construct the structure of the sparse mass matrix
        ConstructMassMatrixStructure(r_model_part, system_size, elements_id_map, m_global_mass_matrix);

        // #pragma omp parallel
        // {
        Matrix local_lhs;
        ModelPart::ElementsContainerType::iterator el_begin = r_model_part.ElementsBegin();
        // #pragma omp for schedule(guided)
        for (unsigned e = 0; e < number_of_elements; e++){
            ModelPart::ElementsContainerType::iterator it_elem = el_begin + e;
            // Node id to global index map
            CalculateLocalMassMatrix(number_of_dofs, it_elem, local_lhs);
            AssembleMassMatrix(m_global_mass_matrix, local_lhs, elements_id_map[e]);
            local_lhs = ZeroMatrix(local_lhs.size1(), local_lhs.size2());
        }
        // }

        mMassMatrixAlreadyComputed = true;
    }

    // Compute the projection of grad(u_j) for each j (RHS)
    // unsigned num_threads = omp_get_num_threads();

    unsigned num_threads = omp_get_max_threads();
    std::vector<SparseSpaceType::VectorType> thread_local_rhs(num_threads, ZeroVector(system_size));
    SparseSpaceType::VectorType ProjectionRHS = ZeroVector(system_size);
    SparseSpaceType::VectorType ProjectionCoefficients = ZeroVector(system_size);
    for(unsigned j = 0; j < TDim; j++)
    {
        // Reset thread local rhs
        for(unsigned t = 0; t < num_threads; t++)
            thread_local_rhs[t] = ZeroVector(system_size);
        ProjectionRHS = ZeroVector(system_size);
        ProjectionCoefficients = ZeroVector(system_size);

        // for (unsigned int e = 0; e < number_of_elements; e++)
        #pragma omp parallel
        {
            // Parallel variables
            auto& rLocalMesh = r_model_part.GetCommunicator().LocalMesh();
            ModelPart::ElementIterator elements_begin = rLocalMesh.ElementsBegin();

            unsigned tid = omp_get_thread_num();
            auto& local_rhs = thread_local_rhs[tid];

            #pragma omp for
            for(unsigned e = 0; e<static_cast<int>(rLocalMesh.NumberOfElements()); ++e)
            {
                const ModelPart::ElementIterator it_elem = elements_begin + e;
                // ModelPart::ElementsContainerType::iterator it_elem = r_model_part.ElementsBegin() + e;

                Geometry<Node>& r_geometry = it_elem->GetGeometry();
                unsigned int NumNodes = r_geometry.size();
                GeometryData::IntegrationMethod integration_method = it_elem->GetIntegrationMethod();
                const std::vector<IntegrationPoint<3>> r_integrations_points = r_geometry.IntegrationPoints( integration_method );
                unsigned int r_number_integration_points = r_geometry.IntegrationPointsNumber(integration_method);
                Vector detJ_vector(r_number_integration_points);
                r_geometry.DeterminantOfJacobian(detJ_vector, integration_method);
                Matrix NContainer = r_geometry.ShapeFunctionsValues(integration_method);
                DenseVector<Matrix> shape_derivatives;
                r_geometry.ShapeFunctionsIntegrationPointsGradients(shape_derivatives, detJ_vector, integration_method);

                for (unsigned int g = 0; g < r_number_integration_points; g++)
                {
                    double Weight = r_integrations_points[g].Weight() * detJ_vector[g];

                    // Calculate the gauss global coordinates
                    array_1d<double, 3> gauss_coords = ZeroVector(3);
                    for(unsigned n = 0; n < NumNodes; n++)
                    {
                        Vector node_coords = r_geometry[n].Coordinates();
                        double shape_function_at_gauss_point = NContainer(g, n);
                        for(unsigned d = 0; d < TDim; d++)
                        {
                            gauss_coords[d] += node_coords[d] * shape_function_at_gauss_point;
                        }
                    }

                    // Calculate the exact gradient at x = xg
                    array_1d<array_1d<double, 3>, 3> grad_exact;
                    flow_field->CalculateGradient(0.0, gauss_coords, grad_exact, omp_get_thread_num());
                    for (unsigned int a = 0; a < NumNodes; ++a){
                        double shape_function_at_gauss_point = NContainer(g, a);
                        for (unsigned d = 0; d < TDim; d++)
                        {
                            unsigned int a_global = TDim * id_to_position[r_geometry[a].Id()] + d;
                            double projection_rhs_value = Weight * grad_exact[j][d] * shape_function_at_gauss_point;
                            local_rhs[a_global] += projection_rhs_value;
                        }
                    }
                }
            }
        }
        // Reduction outside the parallel region
        for (unsigned tid = 0; tid < num_threads; tid++)
        {
            for (int i = 0; i < system_size; ++i) {
                ProjectionRHS[i] += thread_local_rhs[tid][i];
            }
        }

        // Invert mass matrix
        AMGCLSolver<SparseSpaceType, DenseSpaceType> LinearSolver;
        LinearSolver.Solve(m_global_mass_matrix, ProjectionCoefficients, ProjectionRHS);

        // Add values of material derivatives
        #pragma omp parallel for
        for (unsigned int i = 0; i < number_of_nodes; i++)
        {
            ModelPart::NodesContainerType::iterator inode = r_model_part.NodesBegin() + i;
            array_1d <double, 3>& material_derivative_node = inode->FastGetSolutionStepValue(material_derivative_container);
            array_1d <double, 3>& elemental_values = inode->FastGetSolutionStepValue(vector_container);
            for (unsigned int d = 0; d < TDim; d++)
            {
                unsigned int index = TDim * id_to_position[inode->Id()] + d;
                material_derivative_node[j] += elemental_values[d] * ProjectionCoefficients[index];

                // Store the gradient of u_j, i.e. \partial_d u_j, d = 0, 1, 2
                if (mStoreFullGradient)
                {
                    if (j == 0)
                    {
                        array_1d<double, 3>& gradient = inode->FastGetSolutionStepValue(VELOCITY_X_GRADIENT);
                        gradient[d] = ProjectionCoefficients[index];
                    } else if (j == 1)
                    {
                        array_1d<double, 3>& gradient = inode->FastGetSolutionStepValue(VELOCITY_Y_GRADIENT);
                        gradient[d] = ProjectionCoefficients[index];
                    } else if (j == 2)
                    {
                        array_1d<double, 3>& gradient = inode->FastGetSolutionStepValue(VELOCITY_Z_GRADIENT);
                        gradient[d] = ProjectionCoefficients[index];
                    }
                }
            }
        }

        // Reset variables
        // ProjectionRHS = ZeroVector(system_size);
        // ProjectionCoefficients = ZeroVector(system_size);
    }
    KRATOS_INFO("SwimmingDEM") << "Finished constructing the material derivative for a given flow field." << std::endl;

}
//***************************************************************************************************************
//***************************************************************************************************************
template <std::size_t TDim>
void DerivativeRecovery<TDim>::CalculateVectorMaterialDerivativeExactL2Parallel(ModelPart& r_model_part,
                                                                 Variable<array_1d<double, 3> >& vector_container,
                                                                 Variable<array_1d<double, 3> >& vector_rate_container,
                                                                 Variable<array_1d<double, 3> >& material_derivative_container)
{
    KRATOS_INFO("SwimmingDEM") << "Constructing the material derivative by computing the L2 projection..." << std::endl;

    // System size
    const unsigned number_of_elements = r_model_part.NumberOfElements();
    const unsigned number_of_nodes = r_model_part.NumberOfNodes();
    const unsigned number_of_dofs = TDim;  // vx, vy, vz (vi,j for i, j = 0, ..., TDim - 1)
    const SizeType system_size = number_of_nodes * number_of_dofs;

    // Map the ids of the nodes to their position in the model part
    std::map <std::size_t, unsigned int> id_to_position;
    unsigned int node_pos = 0;
    for (NodeIteratorType inode = r_model_part.NodesBegin(); inode != r_model_part.NodesEnd(); ++inode){
        noalias(inode->FastGetSolutionStepValue(material_derivative_container)) = ZeroVector(3);
        id_to_position[inode->Id()] = node_pos;
        ++node_pos;
    }

    // Compute the mass matrix of the system (only once)
    if(!mMassMatrixAlreadyComputed)
    {
        // Map of local dof position to global dof position
        std::vector<std::map<unsigned, unsigned>> elements_id_map(number_of_elements);
        for (unsigned int e = 0; e < number_of_elements; e++)
        {
            ModelPart::ElementsContainerType::iterator it_elem = r_model_part.ElementsBegin() + e;
            Geometry<Node>& r_geometry = it_elem->GetGeometry();

            std::map<unsigned, unsigned> element_map;
            for(unsigned n = 0; n < r_geometry.size(); n++)
            {
                const SizeType node_pos = id_to_position[r_geometry[n].Id()];
                for (unsigned i = 0; i < number_of_dofs; i++)
                {
                    unsigned local_dof_position = n * number_of_dofs + i;  // Each node has a local position in the element
                    unsigned global_dof_position = node_pos * number_of_dofs + i;  // Each DOF in the model part has a unique position (we assume the nodes in r_geometry are always at the same order)
                    element_map[local_dof_position] = global_dof_position;  // This maps the local DOF position of the element to the global position of this DOF
                }
            }
            elements_id_map[e] = element_map;
        }

        // Construct the structure of the sparse mass matrix
        ConstructMassMatrixStructure(r_model_part, system_size, elements_id_map, m_global_mass_matrix);

        // #pragma omp parallel
        // {
        Matrix local_lhs;
        ModelPart::ElementsContainerType::iterator el_begin = r_model_part.ElementsBegin();
        // #pragma omp for schedule(guided)
        for (unsigned e = 0; e < number_of_elements; e++){
            ModelPart::ElementsContainerType::iterator it_elem = el_begin + e;
            // Node id to global index map
            CalculateLocalMassMatrix(number_of_dofs, it_elem, local_lhs);
            AssembleMassMatrix(m_global_mass_matrix, local_lhs, elements_id_map[e]);
            local_lhs = ZeroMatrix(local_lhs.size1(), local_lhs.size2());
        }
        // }

        mMassMatrixAlreadyComputed = true;
    }

    // Compute the projection of grad(u_j) for each j (RHS)
    // unsigned num_threads = omp_get_num_threads();

    unsigned num_threads = omp_get_max_threads();
    std::vector<SparseSpaceType::VectorType> thread_local_rhs(num_threads, ZeroVector(system_size));
    SparseSpaceType::VectorType ProjectionRHS = ZeroVector(system_size);
    SparseSpaceType::VectorType ProjectionCoefficients = ZeroVector(system_size);
    for(unsigned j = 0; j < TDim; j++)
    {
        // Reset thread local rhs
        for(unsigned t = 0; t < num_threads; t++)
            thread_local_rhs[t] = ZeroVector(system_size);
        ProjectionRHS = ZeroVector(system_size);
        ProjectionCoefficients = ZeroVector(system_size);

        // for (unsigned int e = 0; e < number_of_elements; e++)
        #pragma omp parallel
        {
            // Parallel variables
            auto& rLocalMesh = r_model_part.GetCommunicator().LocalMesh();
            ModelPart::ElementIterator elements_begin = rLocalMesh.ElementsBegin();

            unsigned tid = omp_get_thread_num();
            auto& local_rhs = thread_local_rhs[tid];

            #pragma omp for
            for(unsigned e = 0; e<static_cast<int>(rLocalMesh.NumberOfElements()); ++e)
            {
                const ModelPart::ElementIterator it_elem = elements_begin + e;
                // ModelPart::ElementsContainerType::iterator it_elem = r_model_part.ElementsBegin() + e;

                Geometry<Node>& r_geometry = it_elem->GetGeometry();
                unsigned int NumNodes = r_geometry.size();
                GeometryData::IntegrationMethod integration_method = it_elem->GetIntegrationMethod();
                const std::vector<IntegrationPoint<3>> r_integrations_points = r_geometry.IntegrationPoints( integration_method );
                unsigned int r_number_integration_points = r_geometry.IntegrationPointsNumber(integration_method);
                Vector detJ_vector(r_number_integration_points);
                r_geometry.DeterminantOfJacobian(detJ_vector, integration_method);
                Matrix NContainer = r_geometry.ShapeFunctionsValues(integration_method);
                DenseVector<Matrix> shape_derivatives;
                r_geometry.ShapeFunctionsIntegrationPointsGradients(shape_derivatives, detJ_vector, integration_method);

                for (unsigned int g = 0; g < r_number_integration_points; g++)
                {
                    double Weight = r_integrations_points[g].Weight() * detJ_vector[g];

                    // Use nodal values
                    for (unsigned int a = 0; a < NumNodes; ++a){
                        for (unsigned int b = 0; b < NumNodes; ++b){
                            for (unsigned int d = 0; d < TDim; d++)
                            {
                                // Mass matrix (LHS)
                                unsigned int a_global = TDim * id_to_position[r_geometry[a].Id()] + d;
                                double nodal_value_j = r_geometry[b].FastGetSolutionStepValue(vector_container)[j];
                                double projection_rhs_value = Weight * nodal_value_j * NContainer(g, a) * shape_derivatives[g](b, d);
                                local_rhs[a_global] += projection_rhs_value;
                            }
                        }
                    }
                }
            }
        }
        // Reduction outside the parallel region
        for (unsigned tid = 0; tid < num_threads; tid++)
        {
            for (int i = 0; i < system_size; ++i) {
                ProjectionRHS[i] += thread_local_rhs[tid][i];
            }
        }

        // Invert mass matrix
        AMGCLSolver<SparseSpaceType, DenseSpaceType> LinearSolver;
        LinearSolver.Solve(m_global_mass_matrix, ProjectionCoefficients, ProjectionRHS);

        // Add values of material derivatives
        #pragma omp parallel for
        for (unsigned int i = 0; i < number_of_nodes; i++)
        {
            ModelPart::NodesContainerType::iterator inode = r_model_part.NodesBegin() + i;
            array_1d <double, 3>& material_derivative_node = inode->FastGetSolutionStepValue(material_derivative_container);
            array_1d <double, 3>& elemental_values = inode->FastGetSolutionStepValue(vector_container);
            for (unsigned int d = 0; d < TDim; d++)
            {
                unsigned int index = TDim * id_to_position[inode->Id()] + d;
                material_derivative_node[j] += elemental_values[d] * ProjectionCoefficients[index];

                // Store the gradient of u_j, i.e. \partial_d u_j, d = 0, 1, 2
                if (mStoreFullGradient)
                {
                    if (j == 0)
                    {
                        array_1d<double, 3>& gradient = inode->FastGetSolutionStepValue(VELOCITY_X_GRADIENT);
                        gradient[d] = ProjectionCoefficients[index];
                    } else if (j == 1)
                    {
                        array_1d<double, 3>& gradient = inode->FastGetSolutionStepValue(VELOCITY_Y_GRADIENT);
                        gradient[d] = ProjectionCoefficients[index];
                    } else if (j == 2)
                    {
                        array_1d<double, 3>& gradient = inode->FastGetSolutionStepValue(VELOCITY_Z_GRADIENT);
                        gradient[d] = ProjectionCoefficients[index];
                    }
                }
            }
        }

        // Reset variables
        // ProjectionRHS = ZeroVector(system_size);
        // ProjectionCoefficients = ZeroVector(system_size);
    }
    KRATOS_INFO("SwimmingDEM") << "Finished constructing the material derivative by computing the L2 projection." << std::endl;
}
//***************************************************************************************************************
//***************************************************************************************************************
template <std::size_t TDim>
void DerivativeRecovery<TDim>::CalculateVectorMaterialDerivativeExactL2(ModelPart& r_model_part,
                                                                 Variable<array_1d<double, 3> >& vector_container,
                                                                 Variable<array_1d<double, 3> >& vector_rate_container,
                                                                 Variable<array_1d<double, 3> >& material_derivative_container)
{
    KRATOS_INFO("SwimmingDEM") << "Constructing the material derivative by computing the L2 projection (old)..." << std::endl;
    std::map <std::size_t, unsigned int> id_to_position;
    unsigned int entry = 0;
    for (NodeIteratorType inode = r_model_part.NodesBegin(); inode != r_model_part.NodesEnd(); ++inode){
        noalias(inode->FastGetSolutionStepValue(material_derivative_container)) = ZeroVector(3);
        id_to_position[inode->Id()] = entry;
        ++entry;
    }
    unsigned int number_of_nodes = r_model_part.NumberOfNodes();
    unsigned int number_of_elements = r_model_part.NumberOfElements();
    if (number_of_nodes != entry)
    {
        KRATOS_ERROR << "Found discrepancy on the number of nodes in the L2 projection." << std::endl;
    }
    std::cout << "Preparing to compute L2" << std::endl;

    // For each u_j compute the gradient
    Matrix massMatrix = ZeroMatrix(number_of_nodes * TDim, number_of_nodes * TDim);
    Vector rhs = ZeroVector(number_of_nodes * TDim);
    
    // bool mass_matrix_computed = false;
    for (unsigned int j = 0; j < TDim; j++) {
        Vector L2Projection = ZeroVector(number_of_nodes * TDim);  // grad_proj_component_j

        std::cout << "Building mass matrix for j = " << j << std::endl;
        for (unsigned int e = 0; e < number_of_elements; e++)
        {
            // std::cout << "  Computing element e = " << e << std::endl;
            ModelPart::ElementsContainerType::iterator it_elem = r_model_part.ElementsBegin() + e;
            Geometry<Node>& r_geometry = it_elem->GetGeometry();
            unsigned int NumNodes = r_geometry.size();
            GeometryData::IntegrationMethod integration_method = it_elem->GetIntegrationMethod();
            const std::vector<IntegrationPoint<3>> r_integrations_points = r_geometry.IntegrationPoints( integration_method );
            unsigned int r_number_integration_points = r_geometry.IntegrationPointsNumber(integration_method);
            Vector detJ_vector(r_number_integration_points);
            r_geometry.DeterminantOfJacobian(detJ_vector, integration_method);
            Matrix NContainer = r_geometry.ShapeFunctionsValues(integration_method);
            DenseVector<Matrix> shape_derivatives;
            r_geometry.ShapeFunctionsIntegrationPointsGradients(shape_derivatives, detJ_vector, integration_method);

            Matrix localLHS = ZeroMatrix(NumNodes * TDim, NumNodes * TDim);
            for (unsigned int g = 0; g < r_number_integration_points; g++)
            {
                double Weight = r_integrations_points[g].Weight() * detJ_vector[g];
                for (unsigned int a = 0; a < NumNodes; ++a){
                    for (unsigned int b = 0; b < NumNodes; ++b){
                        for (unsigned int d = 0; d < TDim; d++)
                        {
                            // Mass matrix (LHS)
                            unsigned int a_global = TDim * id_to_position[r_geometry[a].Id()] + d;
                            unsigned int b_global = TDim * id_to_position[r_geometry[b].Id()] + d;
                            double nodal_value_j = r_geometry[b].FastGetSolutionStepValue(vector_container)[j];
                            massMatrix(a_global, b_global) += Weight * NContainer(g, a) * NContainer(g, b);
                            rhs[a_global] += Weight * nodal_value_j * NContainer(g, a) * shape_derivatives[g](b, d);

                            // localLHS(a * TDim + d, b * TDim + d) += Weight * NContainer(g, a) * NContainer(g, b);
                        }
                    }
                }
            }
            // Print mass matrix
            // if (e == 89)
            // {
            //     for (unsigned a = 0; a < NumNodes; a++)
            //     {
            //         for(unsigned b = 0; b < NumNodes; b++)
            //         {
            //             for(unsigned d = 0; d < TDim; d++)
            //             {
            //                 unsigned local_ind_a = a * TDim + d;
            //                 unsigned local_ind_b = b * TDim + d;
            //                 for (unsigned int g = 0; g < r_number_integration_points; g++)
            //                 {
            //                     double Weight = r_integrations_points[g].Weight() * detJ_vector[g];
            //                     localLHS(local_ind_a, local_ind_b) += Weight * NContainer(g, a) * NContainer(g, b);
            //                 }
            //                 unsigned int a_global = TDim * id_to_position[r_geometry[a].Id()] + d;
            //                 unsigned int b_global = TDim * id_to_position[r_geometry[b].Id()] + d;
            //                 std::cout << "localMassMatrix[" << local_ind_a << ", " << local_ind_b << "] = " << localLHS(local_ind_a, local_ind_b) << " will be send to (a_global, b_glboal) = (" << a_global << ", " << b_global << ")" << std::endl;
            //             }
            //         }
            //     }
            //     exit(0);
            // }
        }
        // for(unsigned k = 0; k < number_of_nodes * TDim; k++)
        // {
        //     // std::cout << "rhs[" << k << "] = " << rhs[k] << std::endl;
        //     std::cout << "massMatrix[" << k << ", "<< k << "] = " << massMatrix(k, k) << std::endl;
        // }
        // exit(0);

        // Solve system
        std::cout << "Solving for j = " << j << std::endl;
        MathUtils<double>::Solve(massMatrix, L2Projection, rhs);
        std::cout << "Solved for j = " << j << std::endl;

        // Add values of material derivatives
        for (unsigned int i = 0; i < number_of_nodes; i++) {
            ModelPart::NodesContainerType::iterator inode = r_model_part.NodesBegin() + i;
            array_1d <double, 3>& material_derivative_node = inode->FastGetSolutionStepValue(material_derivative_container);
            array_1d <double, 3>& elemental_values = inode->FastGetSolutionStepValue(vector_container);
            for (unsigned int d = 0; d < TDim; d++)
            {
                unsigned int index = TDim * id_to_position[inode->Id()] + d;
                material_derivative_node[j] += elemental_values[d] * L2Projection(index);

                // Store the gradient of u_j, i.e. \partial_d u_j, d = 0, 1, 2
                if (mStoreFullGradient)
                {
                    if (j == 0)
                    {
                        array_1d<double, 3>& gradient = inode->FastGetSolutionStepValue(VELOCITY_X_GRADIENT);
                        gradient[d] = L2Projection(index);
                    } else if (j == 1)
                    {
                        array_1d<double, 3>& gradient = inode->FastGetSolutionStepValue(VELOCITY_Y_GRADIENT);
                        gradient[d] = L2Projection(index);
                    } else if (j == 2)
                    {
                        array_1d<double, 3>& gradient = inode->FastGetSolutionStepValue(VELOCITY_Z_GRADIENT);
                        gradient[d] = L2Projection(index);
                    }
                }
            }
        }

        // Reset variables
        // if (!mass_matrix_computed)
        // {
        massMatrix = ZeroMatrix(TDim * number_of_nodes, TDim * number_of_nodes);
        // mass_matrix_computed = true;
        // }
        rhs = ZeroVector(TDim * number_of_nodes);
    }

    AddTimeDerivative(r_model_part, material_derivative_container);
    KRATOS_INFO("SwimmingDEM") << "Finished constructing the material derivative by computing the L2 projection..." << std::endl;

}
//***************************************************************************************************************
//***************************************************************************************************************
template <std::size_t TDim>
void DerivativeRecovery<TDim>::CalculateVectorMaterialDerivative(ModelPart& r_model_part,
                                                                 Variable<array_1d<double, 3> >& vector_container,
                                                                 Variable<array_1d<double, 3> >& vector_rate_container,
                                                                 Variable<array_1d<double, 3> >& material_derivative_container)
{
    KRATOS_INFO("SwimmingDEM") << "Constructing the material derivative by derivating nodal averages..." << std::endl;
    std::map <std::size_t, unsigned int> id_to_position;
    unsigned int entry = 0;
    for (NodeIteratorType inode = r_model_part.NodesBegin(); inode != r_model_part.NodesEnd(); ++inode){
        noalias(inode->FastGetSolutionStepValue(material_derivative_container)) = ZeroVector(3);
        id_to_position[inode->Id()] = entry;
        ++entry;
    }
    std::vector<array_1d <double, 3> > convective_contributions_to_the_derivative;
    convective_contributions_to_the_derivative.resize(entry);
    array_1d <double, 3> grad = ZeroVector(3);
    // array_1d <double, TDim + 1 > elemental_values;
    // array_1d <double, TDim + 1 > N; // shape functions vector
    // BoundedMatrix<double, TDim + 1, TDim> DN_DX;
    for (unsigned int j = 0; j < TDim; ++j){ // for each component of the original vector value
        // for each element, constructing the gradient contribution (to its nodes) of the component v_j and storing it in material_derivative_container
        for (ModelPart::ElementIterator ielem = r_model_part.ElementsBegin(); ielem != r_model_part.ElementsEnd(); ++ielem){
            // computing the shape function derivatives
            Geometry<Node >& geom = ielem->GetGeometry();
            // FROM HERE, I TRY TO DO AN IMPLEMENTATION FOR ANY ELEMENT
            GeometryData::IntegrationMethod integration_method = ielem->GetIntegrationMethod();
            auto number_integration_points = geom.IntegrationPointsNumber(integration_method);
            const std::vector<IntegrationPoint<3>>& IntegrationPoints = geom.IntegrationPoints(integration_method);
            double NumNodes = geom.size();
            Vector gauss_weights = ZeroVector(number_integration_points);
            Matrix shape_functions = ZeroMatrix(number_integration_points,NumNodes);
            DenseVector<Matrix> shape_derivatives;
            Matrix DN_DX = ZeroMatrix(NumNodes, TDim);
            Vector N = ZeroVector(NumNodes);
            Vector elemental_values = ZeroVector(NumNodes);

            Vector DetJ;
            geom.ShapeFunctionsIntegrationPointsGradients(shape_derivatives,DetJ,integration_method);
            shape_functions = geom.ShapeFunctionsValues(integration_method);

            for (unsigned int g = 0; g < number_integration_points; g++){
                gauss_weights[g] = DetJ[g] * IntegrationPoints[g].Weight();
                for (unsigned int i = 0; i < NumNodes; ++i){
                    for (unsigned int l = 0; l < NumNodes; ++l){
                        for (unsigned int d = 0; d < TDim; ++d)
                            DN_DX(l,d) += gauss_weights[g] * shape_functions(g,i) * shape_derivatives[g](l,d);
                    }
                }
            }

            //GeometryUtils::CalculateGeometryData(geom, DN_DX, N, Volume);
            for (unsigned int i = 0; i < NumNodes; ++i){
                elemental_values[i] = geom[i].FastGetSolutionStepValue(vector_container)[j];
            }
            array_1d <double, 3> grad_aux = prod(trans(DN_DX), elemental_values); // its dimension may be 2
            for (unsigned int i = 0; i < TDim; ++i){
                grad[i] = grad_aux[i];
            }
            // double nodal_area = Volume / static_cast<double>(TDim + 1);
            // grad *= nodal_area;
            for (unsigned int i = 0; i < NumNodes; ++i){
                geom[i].FastGetSolutionStepValue(material_derivative_container) += grad / NumNodes; // we use material_derivative_container to store the gradient of one component at a time
            }
        }

        // normalizing the constributions to the gradient and getting the j-component of the material derivative
        for (NodeIteratorType inode = r_model_part.NodesBegin(); inode != r_model_part.NodesEnd(); ++inode){
            array_1d <double, 3>& stored_gradient_of_component_j = inode->FastGetSolutionStepValue(material_derivative_container);
            stored_gradient_of_component_j /= inode->FastGetSolutionStepValue(NODAL_AREA);
            if (mStoreFullGradient){
                if (j == 0){
                    array_1d <double, 3>& gradient = inode->FastGetSolutionStepValue(VELOCITY_X_GRADIENT);
                    noalias(gradient) = stored_gradient_of_component_j;
                }
                else if (j == 1){
                    array_1d <double, 3>& gradient = inode->FastGetSolutionStepValue(VELOCITY_Y_GRADIENT);
                    noalias(gradient) = stored_gradient_of_component_j;
                }
                else {
                    array_1d <double, 3>& gradient = inode->FastGetSolutionStepValue(VELOCITY_Z_GRADIENT);
                    noalias(gradient) = stored_gradient_of_component_j;
                }
            }
            const array_1d <double, 3>& velocity = inode->FastGetSolutionStepValue(VELOCITY);
            convective_contributions_to_the_derivative[id_to_position[inode->Id()]][j] = DEM_INNER_PRODUCT_3(velocity, stored_gradient_of_component_j);
            stored_gradient_of_component_j = ZeroVector(3);
        }
    }
    // Adding convective part
    for (NodeIteratorType inode = r_model_part.NodesBegin(); inode != r_model_part.NodesEnd(); ++inode){
        const array_1d <double, 3>& stored_convective_contribution = convective_contributions_to_the_derivative[id_to_position[inode->Id()]];
        array_1d <double, 3>& material_derivative = inode->FastGetSolutionStepValue(material_derivative_container);
        material_derivative = stored_convective_contribution;
    }
    // Adding Eulerian time derivative contribution
    AddTimeDerivative(r_model_part, material_derivative_container);
    KRATOS_INFO("SwimmingDEM") << "Finished constructing the material derivative by derivating nodal averages..." << std::endl;
}
//***************************************************************************************************************
//***************************************************************************************************************
template <std::size_t TDim>
void DerivativeRecovery<TDim>::RecoverLagrangianAcceleration(ModelPart& r_model_part)
{
    for (NodeIteratorType inode = r_model_part.NodesBegin(); inode != r_model_part.NodesEnd(); ++inode){
        const array_1d <double, 3>& lagrangian_acceleration = inode->FastGetSolutionStepValue(ACCELERATION);
        array_1d <double, 3>& material_acceleration = inode->FastGetSolutionStepValue(MATERIAL_ACCELERATION);
        noalias(material_acceleration) = lagrangian_acceleration;
    }
}
//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
// This function modifies the material derivative using the pre-computed value of the gradient
template <std::size_t TDim>
void DerivativeRecovery<TDim>::CalculateVectorMaterialDerivativeFromGradient(ModelPart& r_model_part,
                                                                             Variable<array_1d<double, 3> >& vector_gradient_container_x,
                                                                             Variable<array_1d<double, 3> >& vector_gradient_container_y,
                                                                             Variable<array_1d<double, 3> >& vector_gradient_container_z,
                                                                             Variable<array_1d<double, 3> >& vector_rate_container,
                                                                             Variable<array_1d<double, 3>  >& material_derivative_container)
{
    #pragma omp parallel for
    for (int i = 0; i < (int)r_model_part.Nodes().size(); ++i){
        NodeIteratorType i_particle = r_model_part.NodesBegin() + i;
        Node::Pointer p_node = *(i_particle.base());
        const array_1d <double, 3>& velocity = p_node->FastGetSolutionStepValue(VELOCITY);
        array_1d <double, 3>& material_derivative = p_node->FastGetSolutionStepValue(material_derivative_container);
        material_derivative[0] = DEM_INNER_PRODUCT_3(velocity, p_node->FastGetSolutionStepValue(vector_gradient_container_x));
        material_derivative[1] = DEM_INNER_PRODUCT_3(velocity, p_node->FastGetSolutionStepValue(vector_gradient_container_y));
        material_derivative[2] = DEM_INNER_PRODUCT_3(velocity, p_node->FastGetSolutionStepValue(vector_gradient_container_z));
    }
    AddTimeDerivative(r_model_part, material_derivative_container);
}
//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
// This function modifies one component of the material derivative using the pre-computed value of the gradient of the first component of the velocity
template <std::size_t TDim>
void DerivativeRecovery<TDim>::CalculateVectorMaterialDerivativeComponent(ModelPart& r_model_part, Variable<array_1d<double, 3> >& vector_component_gradient_container, Variable<array_1d<double, 3> >& vector_rate_container, Variable<array_1d<double, 3>  >& material_derivative_container)
{
    int current_component = r_model_part.GetProcessInfo()[CURRENT_COMPONENT];
    KRATOS_ERROR_IF(current_component != 0 && current_component != 1 && current_component != 2) << "The value of CURRENT_COMPONENT passed to the ComputeComponentGradientSimplex element is not 0, 1 or 2, but " << current_component << std::endl;
    for (NodeIteratorType inode = r_model_part.NodesBegin(); inode != r_model_part.NodesEnd(); ++inode){
        const array_1d <double, 3>& gradient_of_component = inode->FastGetSolutionStepValue(vector_component_gradient_container);
        const array_1d <double, 3>& velocity = inode->FastGetSolutionStepValue(VELOCITY);
        array_1d <double, 3>& material_derivative = inode->FastGetSolutionStepValue(material_derivative_container);
        material_derivative[current_component] = DEM_INNER_PRODUCT_3(velocity, gradient_of_component);
    }
    AddTimeDerivativeComponent(r_model_part, material_derivative_container, current_component);
}
//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
template <std::size_t TDim>
void DerivativeRecovery<TDim>::RecoverSuperconvergentMatDeriv(ModelPart& r_model_part, Variable<array_1d<double, 3> >& vector_container, Variable<array_1d<double, 3> >& vector_rate_container, Variable<array_1d<double, 3> >& material_derivative_container, unsigned int& ord)
{
    mCalculatingTheGradient = true;
    if (mFirstGradientRecovery){
        KRATOS_INFO("SwimmingDEM") << "Constructing first-step neighbour clouds for material derivative..." << std::endl;
        SetNeighboursAndWeights(r_model_part, ord);
        mFirstGradientRecovery = false;
        KRATOS_INFO("SwimmingDEM") << "Finished constructing neighbour clouds for material derivative." << std::endl;
    }
    if (mSomeCloudsDontWork){ // a default value is necessary in the cases where recovery is not possible
        CalculateVectorMaterialDerivative(r_model_part, vector_container, vector_rate_container, material_derivative_container);
        std::cout << "Error in RecoverSuperConvergentMatDeriv: mSomeCloudsDontWork = " << mSomeCloudsDontWork << std::endl;
        exit(1);
    }
    // Solving least squares problem (Zhang, 2006)
    unsigned int n_relevant_terms = 3;
    std::vector<array_1d <double, 3> > polynomial_coefficients; // vector to store, for each node, the corresponding values of the polynomial coefficients relevant for the calculation of the material derivative
    polynomial_coefficients.resize(n_relevant_terms);
    for (NodeIteratorType inode = r_model_part.NodesBegin(); inode != r_model_part.NodesEnd(); ++inode){
        GlobalPointersVector<Node >& neigh_nodes = inode->GetValue(NEIGHBOUR_NODES);
        unsigned int n_neigh = neigh_nodes.size();
        if (!n_neigh){ // then we keep the defualt value
            continue;
        }
        for (unsigned int i = 0; i < n_relevant_terms; ++i){ // resetting polynomial_coefficients to 0
            polynomial_coefficients[i] = ZeroVector(3);
        }
        const Vector& nodal_weights = inode->FastGetSolutionStepValue(NODAL_WEIGHTS);
        for (unsigned int k = 0; k < TDim; ++k){
            for (unsigned int i_neigh = 0; i_neigh < n_neigh; ++i_neigh){
                const array_1d<double, 3>& neigh_nodal_value = neigh_nodes[i_neigh].FastGetSolutionStepValue(vector_container);
                for (unsigned int d = 0; d < n_relevant_terms; ++d){
                    polynomial_coefficients[d][k] += nodal_weights[n_relevant_terms * i_neigh + d] * neigh_nodal_value[k];
                }
            }
        }
        array_1d <double, 3>& recovered_mat_deriv = inode->FastGetSolutionStepValue(material_derivative_container);
        const array_1d <double, 3>& velocity = inode->FastGetSolutionStepValue(vector_container);
        recovered_mat_deriv[0] = velocity[0] * polynomial_coefficients[0][0] + velocity[1] * polynomial_coefficients[1][0] + velocity[2] * polynomial_coefficients[2][0];
        recovered_mat_deriv[1] = velocity[0] * polynomial_coefficients[0][1] + velocity[1] * polynomial_coefficients[1][1] + velocity[2] * polynomial_coefficients[2][1];
        recovered_mat_deriv[2] = velocity[0] * polynomial_coefficients[0][2] + velocity[1] * polynomial_coefficients[1][2] + velocity[2] * polynomial_coefficients[2][2];
    }
    AddTimeDerivative(r_model_part, material_derivative_container);
    mCalculatingTheGradient = false;
}
//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
// This function adds the contribution to the vorticity corresponding to the pre-computed gradient of one velocity component
template <std::size_t TDim>
void DerivativeRecovery<TDim>::CalculateVorticityFromGradient(ModelPart& r_model_part, Variable<array_1d<double, 3> >& vector_gradient_container_x, Variable<array_1d<double, 3> >& vector_gradient_container_y, Variable<array_1d<double, 3> >& vector_gradient_container_z, Variable<array_1d<double, 3>  >& vorticity_container)
{
    #pragma omp parallel for
    for (int i = 0; i < (int)r_model_part.Nodes().size(); ++i){
        NodeIteratorType i_particle = r_model_part.NodesBegin() + i;
        Node::Pointer p_node = *(i_particle.base());
        array_1d<double, 3>& vorticity = p_node->FastGetSolutionStepValue(vorticity_container);
        const array_1d<double, 3>& gradient_x = p_node->FastGetSolutionStepValue(vector_gradient_container_x);
        const array_1d<double, 3>& gradient_y = p_node->FastGetSolutionStepValue(vector_gradient_container_y);
        const array_1d<double, 3>& gradient_z = p_node->FastGetSolutionStepValue(vector_gradient_container_z);
        vorticity[0] = gradient_z[1] - gradient_y[2];
        vorticity[1] = gradient_x[2] - gradient_z[0];
        vorticity[2] = gradient_y[0] - gradient_x[1];
    }
}
//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
// This function adds the contribution to the vorticity corresponding to the pre-computed gradient of one velocity component
template <std::size_t TDim>
void DerivativeRecovery<TDim>::CalculateVorticityContributionOfTheGradientOfAComponent(ModelPart& r_model_part, Variable<array_1d<double, 3> >& vector_component_gradient_container, Variable<array_1d<double, 3>  >& vorticity_container)
{
    int current_component = r_model_part.GetProcessInfo()[CURRENT_COMPONENT];
    KRATOS_ERROR_IF(current_component != 0 && current_component != 1 && current_component != 2) << "The value of CURRENT_COMPONENT passed to the ComputeComponentGradientSimplex element is not 0, 1 or 2, but " << current_component << std::endl;
    for (NodeIteratorType inode = r_model_part.NodesBegin(); inode != r_model_part.NodesEnd(); ++inode){
        const array_1d <double, 3>& gradient_of_component = inode->FastGetSolutionStepValue(vector_component_gradient_container);
        array_1d <double, 3>& vorticity = inode->FastGetSolutionStepValue(vorticity_container);
        if (current_component == 0){
            vorticity[1] += gradient_of_component[2];
            vorticity[2] -= gradient_of_component[1];
        }
        else if (current_component == 1){
            vorticity[0] -= gradient_of_component[2];
            vorticity[2] += gradient_of_component[0];
        }
        else {
            vorticity[0] += gradient_of_component[1];
            vorticity[1] -= gradient_of_component[0];
        }
    }
}
//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
template <std::size_t TDim>
template <class TScalarVariable>
void DerivativeRecovery<TDim>::CalculateGradient(ModelPart& r_model_part, TScalarVariable& scalar_container, Variable<array_1d<double, 3> >& gradient_container)
{
    for (NodeIteratorType inode = r_model_part.NodesBegin(); inode != r_model_part.NodesEnd(); ++inode){
        noalias(inode->FastGetSolutionStepValue(gradient_container)) = ZeroVector(3);
    }
    array_1d <double, 3> grad = ZeroVector(3); // its dimension is always 3
    //array_1d <double, TDim + 1 > elemental_values;
    //array_1d <double, TDim + 1 > N; // shape functions vector
    //BoundedMatrix<double, TDim + 1, TDim> DN_DX;
    for (ModelPart::ElementIterator ielem = r_model_part.ElementsBegin(); ielem != r_model_part.ElementsEnd(); ++ielem){
        // computing the shape function derivatives
        Geometry<Node >& geom = ielem->GetGeometry();
        // FROM HERE, I TRY TO DO AN IMPLEMENTATION FOR ANY ELEMENT
        GeometryData::IntegrationMethod integration_method = ielem->GetIntegrationMethod();
        auto number_integration_points = geom.IntegrationPointsNumber(integration_method);
        const std::vector<IntegrationPoint<3>>& IntegrationPoints = geom.IntegrationPoints(integration_method);

        double NumNodes = geom.size();
        Vector gauss_weights = ZeroVector(number_integration_points);
        Matrix shape_functions = ZeroMatrix(number_integration_points,NumNodes);
        DenseVector<Matrix> shape_derivatives;
        Matrix DN_DX = ZeroMatrix(NumNodes, TDim);
        Vector N = ZeroVector(NumNodes);
        Vector elemental_values = ZeroVector(NumNodes);
        // computing DN_DX values for the strain rate
        //ielem->CalculateGeometryData(gauss_weights, shape_functions, shape_derivatives);
        Vector DetJ;
        geom.ShapeFunctionsIntegrationPointsGradients(shape_derivatives,DetJ,integration_method);
        shape_functions = geom.ShapeFunctionsValues(integration_method);

        //double Volume = geom.Volume();
        for (unsigned int g = 0; g < number_integration_points; g++){
            gauss_weights[g] = DetJ[g] * IntegrationPoints[g].Weight();
            for (unsigned int i = 0; i < NumNodes; ++i){
                for (unsigned int j = 0; j < NumNodes; ++j){
                    for (unsigned int d = 0; d < TDim; ++d)
                        DN_DX(j,d) += gauss_weights[g] * shape_functions(g,i) * shape_derivatives[g](j,d);
                }
            }
        }

        //GeometryUtils::CalculateGeometryData(geom, DN_DX, N, Volume);
        // getting the gradients;
        for (unsigned int i = 0; i < NumNodes; ++i){
            elemental_values[i] = geom[i].FastGetSolutionStepValue(scalar_container);
        }
        array_1d <double, TDim> grad_aux = prod(trans(DN_DX), elemental_values); // its dimension may be 2
        for (unsigned int d = 0; d < TDim; ++d){
            grad[d] = grad_aux[d];
        }
        // AVERAGING
        for (unsigned int i = 0; i < NumNodes; ++i){
            geom[i].FastGetSolutionStepValue(gradient_container) += grad/NumNodes;
        }
    }
    for (NodeIteratorType inode = r_model_part.NodesBegin(); inode != r_model_part.NodesEnd(); ++inode){
        inode->FastGetSolutionStepValue(gradient_container) /= inode->FastGetSolutionStepValue(NODAL_AREA);
    }
}
//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
template <std::size_t TDim>
void DerivativeRecovery<TDim>::SmoothVectorField(ModelPart& r_model_part, Variable<array_1d<double, 3> >& vector_field, Variable<array_1d<double, 3> >& auxiliary_veriable)
{
    for (NodeIteratorType inode = r_model_part.NodesBegin(); inode != r_model_part.NodesEnd(); ++inode){
        noalias(inode->FastGetSolutionStepValue(auxiliary_veriable)) = ZeroVector(3);
    }
    array_1d <double, TDim + 1 > N; // shape functions vector
    BoundedMatrix<double, TDim + 1, TDim> DN_DX;
    for (ModelPart::ElementIterator ielem = r_model_part.ElementsBegin(); ielem != r_model_part.ElementsEnd(); ++ielem){
        // computing the shape function derivatives
        Geometry<Node >& geom = ielem->GetGeometry();
        double Volume;
        GeometryUtils::CalculateGeometryData(geom, DN_DX, N, Volume);
        array_1d <double, 3> average = ZeroVector(3); // its dimension is always 3
        for (unsigned int i = 0; i < TDim; ++i){
            noalias(average) += geom[i].FastGetSolutionStepValue(vector_field);
        }
        double nodal_area = Volume / static_cast<double>(TDim + 1);
        average *= nodal_area;
        for (unsigned int i = 0; i < TDim + 1; ++i){
            geom[i].FastGetSolutionStepValue(auxiliary_veriable) += average;
        }
    }
    for (NodeIteratorType inode = r_model_part.NodesBegin(); inode != r_model_part.NodesEnd(); ++inode){
        noalias(inode->FastGetSolutionStepValue(vector_field)) = inode->FastGetSolutionStepValue(auxiliary_veriable) / (3 * inode->FastGetSolutionStepValue(NODAL_AREA));
    }
}
//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
template <std::size_t TDim>
template <class TScalarVariable>
void DerivativeRecovery<TDim>::RecoverSuperconvergentGradient(ModelPart& r_model_part, TScalarVariable& scalar_container, Variable<array_1d<double, 3> >& gradient_container, unsigned int& ord)
{
    mCalculatingTheGradient = true;
    if (mFirstGradientRecovery){
        KRATOS_INFO("SwimmingDEM") << "Constructing first-step neighbour clouds for gradient recovery..." << std::endl;
        SetEdgeNodesAndWeights(r_model_part);
        SetNeighboursAndWeights(r_model_part, ord);
        mFirstGradientRecovery = false;
        KRATOS_INFO("SwimmingDEM") << "Finished constructing neighbour clouds for gradient recovery." << std::endl;
    }
    if (mSomeCloudsDontWork){ // a default value is necessary in the cases where recovery is not possible
        CalculateGradient(r_model_part, scalar_container, gradient_container);
    }

    // Solving least squares problem (Zhang, 2006)
    for (NodeIteratorType inode = r_model_part.NodesBegin(); inode != r_model_part.NodesEnd(); ++inode){
        GlobalPointersVector<Node >& neigh_nodes = inode->GetValue(NEIGHBOUR_NODES);
        unsigned int n_neigh = neigh_nodes.size();
        if (!n_neigh){ // we keep the defualt value
            continue;
        }

        array_1d <double, 3>& recovered_gradient = inode->FastGetSolutionStepValue(gradient_container);
        // recovered_gradient = ZeroVector(3);
        Vector gradient = ZeroVector(3);
        const Vector& nodal_weights = inode->GetValue(NODAL_WEIGHTS);
        std::cout << "Node with Id = " << inode->Id() << std::endl;
        for (unsigned int i_neigh = 0; i_neigh < n_neigh; ++i_neigh){
            // If the node is an edge node compute the gradient differently
            if(inode->GetValue(IS_EDGE_NODE))
            {
                continue;
            } else {
                const double& neigh_nodal_value = neigh_nodes[i_neigh].FastGetSolutionStepValue(scalar_container);
                for (unsigned int d = 0; d < TDim; ++d){
                    gradient[d] += nodal_weights[TDim * i_neigh + d] * neigh_nodal_value;
                    std::cout << "  Neighbour with ID = " << neigh_nodes[i_neigh].Id() << " has w(d = " << d << ") = " << nodal_weights[TDim * i_neigh + d] << std::endl;
                }
            }
        }
        exit(0);
        // std::cout << "gradient = " << gradient << std::endl;
        // exit(0);

        for (unsigned int d = 0; d < TDim; ++d)
        {
            recovered_gradient[d] = gradient[d];
        }
    }
    mCalculatingTheGradient = false;
}
//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
template <std::size_t TDim>
void DerivativeRecovery<TDim>::RecoverSuperconvergentLaplacian(ModelPart& r_model_part, Variable<array_1d<double, 3> >& vector_container, Variable<array_1d<double, 3> >& laplacian_container, unsigned int& ord)
{
    std::cout << "Entering RecoverSuperconvergentLaplacian" << std::endl;

    mCalculatingTheLaplacian = true;
    if (mFirstLaplacianRecovery){
        KRATOS_INFO("SwimmingDEM") << "Constructing first-step neighbour clouds for laplacian recovery..." << std::endl;
        SetNeighboursAndWeights(r_model_part, ord);
        mFirstLaplacianRecovery = false;
        KRATOS_INFO("SwimmingDEM") << "Finished constructing neighbour clouds for laplacian recovery." << std::endl;
    }
    if (mSomeCloudsDontWork){ // a default value is necessary in the cases where recovery is not possible
        CalculateVectorLaplacian(r_model_part, vector_container, laplacian_container);
    }
    // Solving a least squares problem for every node (Zhang, 2006)
    unsigned int n_relevant_terms = 6;
    std::vector<array_1d <double, 3> > polynomial_coefficients; // vector to store, for each node, the corresponding values of the polynomial coefficients relevant for the calculation of the laplacian
    polynomial_coefficients.resize(n_relevant_terms);
    for (NodeIteratorType inode = r_model_part.NodesBegin(); inode != r_model_part.NodesEnd(); ++inode){
        GlobalPointersVector<Node >& neigh_nodes = inode->GetValue(NEIGHBOUR_NODES);
        unsigned int n_neigh = neigh_nodes.size();
        if (!n_neigh){ // then we keep the defualt value
            continue;
        }
        for (unsigned int i = 0; i < n_relevant_terms; ++i){ // resetting polynomial_coefficients to 0
            polynomial_coefficients[i] = ZeroVector(3);
        }
        const Vector& nodal_weights = inode->FastGetSolutionStepValue(NODAL_WEIGHTS);
        array_1d <double, 3>& recovered_laplacian = inode->FastGetSolutionStepValue(laplacian_container);
        noalias(recovered_laplacian) = ZeroVector(3);
        for (unsigned int k = 0; k < TDim; ++k){
           for (unsigned int i_neigh = 0; i_neigh < n_neigh; ++i_neigh){
                const array_1d<double, 3>& neigh_nodal_value = neigh_nodes[i_neigh].FastGetSolutionStepValue(vector_container);
                for (unsigned int d = 0; d < n_relevant_terms; ++d){
                    polynomial_coefficients[d][k] += nodal_weights[n_relevant_terms * i_neigh + d] * neigh_nodal_value[k];
                }
            }
        }
        recovered_laplacian[0] = 2 * (polynomial_coefficients[3][0] + polynomial_coefficients[4][0] + polynomial_coefficients[5][0]);
        recovered_laplacian[1] = 2 * (polynomial_coefficients[3][1] + polynomial_coefficients[4][1] + polynomial_coefficients[5][1]);
        recovered_laplacian[2] = 2 * (polynomial_coefficients[3][2] + polynomial_coefficients[4][2] + polynomial_coefficients[5][2]);
//        recovered_laplacian[0] = polynomial_coefficients[0][0] + polynomial_coefficients[1][0] + 2 * polynomial_coefficients[3][0];
//        recovered_laplacian[1] = polynomial_coefficients[0][0] + polynomial_coefficients[2][0] + 2 * polynomial_coefficients[4][0];
//        recovered_laplacian[2] = polynomial_coefficients[1][0] + polynomial_coefficients[2][0] + 2 * polynomial_coefficients[5][0];
    }
    mCalculatingTheLaplacian = false;

    std::cout << "Exiting RecoverSuperconvergentLaplacian" << std::endl;
}
//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
template <std::size_t TDim>
void DerivativeRecovery<TDim>::RecoverSuperconvergentVelocityLaplacianFromGradient(ModelPart& r_model_part, Variable<array_1d<double, 3> >& vector_container, Variable<array_1d<double, 3> >& laplacian_container, unsigned int& ord)
{
    std::cout << "Entering RecoverSuperconvergentVelocityLaplacianFromGradient" << std::endl;

    mCalculatingTheGradient = true;
    if (mFirstLaplacianRecovery){
        KRATOS_INFO("SwimmingDEM") << "Finished constructing neighbour clouds for laplacian recovery." << std::endl;
        SetNeighboursAndWeights(r_model_part, ord);
        mFirstLaplacianRecovery = false;
        KRATOS_INFO("SwimmingDEM") << "Finished constructing neighbour clouds for laplacian recovery." << std::endl;
    }
    if (mSomeCloudsDontWork){ // a default value is necessary in the cases where recovery is not possible
        CalculateVectorLaplacian(r_model_part, vector_container, laplacian_container);
    }
    // Solving least squares problem (Zhang, 2006)
    unsigned int n_relevant_terms = 3;
    std::vector<array_1d <double, 3> > polynomial_coefficients; // vector to store, for each node, the corresponding values of the polynomial coefficients relevant for the calculation of the laplacian
    polynomial_coefficients.resize(n_relevant_terms);
    array_1d< array_1d<double, 3>, 3> gradient;
    for (NodeIteratorType inode = r_model_part.NodesBegin(); inode != r_model_part.NodesEnd(); ++inode){
        GlobalPointersVector<Node >& neigh_nodes = inode->GetValue(NEIGHBOUR_NODES);
        unsigned int n_neigh = neigh_nodes.size();
        if (!n_neigh){ // then we keep the defualt value
            continue;
        }
        for (unsigned int i = 0; i < n_relevant_terms; ++i){ // resetting polynomial_coefficients to 0
            polynomial_coefficients[i] = ZeroVector(3);
        }
        const Vector& nodal_weights = inode->FastGetSolutionStepValue(NODAL_WEIGHTS);
        array_1d <double, 3>& recovered_laplacian = inode->FastGetSolutionStepValue(laplacian_container);
        recovered_laplacian =  ZeroVector(3);
        for (unsigned int i_neigh = 0; i_neigh < n_neigh; ++i_neigh){
            noalias(gradient[0]) = neigh_nodes[i_neigh].FastGetSolutionStepValue(VELOCITY_X_GRADIENT);
            noalias(gradient[1]) = neigh_nodes[i_neigh].FastGetSolutionStepValue(VELOCITY_Y_GRADIENT);
            noalias(gradient[2]) = neigh_nodes[i_neigh].FastGetSolutionStepValue(VELOCITY_Z_GRADIENT);
            recovered_laplacian[0] += nodal_weights[n_relevant_terms * i_neigh] * gradient[0][0];
            recovered_laplacian[0] += nodal_weights[n_relevant_terms * i_neigh + 1] * gradient[0][1];
            recovered_laplacian[0] += nodal_weights[n_relevant_terms * i_neigh + 2] * gradient[0][2];
            recovered_laplacian[1] += nodal_weights[n_relevant_terms * i_neigh] * gradient[1][0];
            recovered_laplacian[1] += nodal_weights[n_relevant_terms * i_neigh + 1] * gradient[1][1];
            recovered_laplacian[1] += nodal_weights[n_relevant_terms * i_neigh + 2] * gradient[1][2];
            recovered_laplacian[2] += nodal_weights[n_relevant_terms * i_neigh] * gradient[2][0];
            recovered_laplacian[2] += nodal_weights[n_relevant_terms * i_neigh + 1] * gradient[2][1];
            recovered_laplacian[2] += nodal_weights[n_relevant_terms * i_neigh + 2] * gradient[2][2];
        }
    }
    mCalculatingTheGradient = false;
    std::cout << "Exiting RecoverSuperconvergentVelocityLaplacianFromGradient" << std::endl;
}
//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
template <std::size_t TDim>
void DerivativeRecovery<TDim>::RecoverSuperconvergentMatDerivAndLaplacian(ModelPart& r_model_part, Variable<array_1d<double, 3> >& vector_container, Variable<array_1d<double, 3> >& vector_rate_container, Variable<array_1d<double, 3> >& material_derivative_container, Variable<array_1d<double, 3> >& laplacian_container, unsigned int& ord)
{
    std::cout << "Entering RecoverSuperconvergentMatDerivAndLaplacian" << std::endl;

    mCalculatingGradientAndLaplacian = true;
    if (mFirstLaplacianRecovery){
        KRATOS_INFO("SwimmingDEM") << "Constructing first-step neighbour clouds for material derivative and laplacian recovery..." << std::endl;
        SetNeighboursAndWeights(r_model_part, ord);
        mFirstLaplacianRecovery = false;
        KRATOS_INFO("SwimmingDEM") << "Finished constructing neighbour clouds for material derivative and laplacian recovery." << std::endl;
    }
    if (mSomeCloudsDontWork){ // a default value is necessary in the cases where recovery is not possible
        CalculateVectorLaplacian(r_model_part, vector_container, laplacian_container);
        CalculateVectorMaterialDerivative(r_model_part, vector_container, vector_rate_container, material_derivative_container);
    }
    // Solving least squares problem (Zhang, 2006)
    unsigned int n_relevant_terms = 9;
    std::vector<array_1d <double, 3> > polynomial_coefficients; // vector to store, for each node, the corresponding values of the polynomial coefficients relevant for the calculation of the laplacian and the material derivative
    polynomial_coefficients.resize(n_relevant_terms);
    for (NodeIteratorType inode = r_model_part.NodesBegin(); inode != r_model_part.NodesEnd(); ++inode){
        GlobalPointersVector<Node >& neigh_nodes = inode->GetValue(NEIGHBOUR_NODES);
        unsigned int n_neigh = neigh_nodes.size();
        if (!n_neigh){ // then we keep the defualt value
            continue;
        }
        for (unsigned int i = 0; i < n_relevant_terms; ++i){ // resetting polynomial_coefficients to 0
            polynomial_coefficients[i] = ZeroVector(3);
        }
        const Vector& nodal_weights = inode->FastGetSolutionStepValue(NODAL_WEIGHTS);
        for (unsigned int k = 0; k < TDim; ++k){
            for (unsigned int i_neigh = 0; i_neigh < n_neigh; ++i_neigh){
                const array_1d<double, 3>& neigh_nodal_value = neigh_nodes[i_neigh].FastGetSolutionStepValue(vector_container);
                for (unsigned int d = 0; d < n_relevant_terms; ++d){
                    polynomial_coefficients[d][k] += nodal_weights[n_relevant_terms * i_neigh + d] * neigh_nodal_value[k];
                }
            }
        }
        array_1d <double, 3>& recovered_laplacian = inode->FastGetSolutionStepValue(laplacian_container);
        array_1d <double, 3>& recovered_mat_deriv = inode->FastGetSolutionStepValue(material_derivative_container);
        const array_1d <double, 3>& velocity = inode->FastGetSolutionStepValue(vector_container);
        recovered_mat_deriv[0] = velocity[0] * polynomial_coefficients[0][0] + velocity[1] * polynomial_coefficients[1][0] + velocity[2] * polynomial_coefficients[2][0];
        recovered_mat_deriv[1] = velocity[0] * polynomial_coefficients[0][1] + velocity[1] * polynomial_coefficients[1][1] + velocity[2] * polynomial_coefficients[2][1];
        recovered_mat_deriv[2] = velocity[0] * polynomial_coefficients[0][2] + velocity[1] * polynomial_coefficients[1][2] + velocity[2] * polynomial_coefficients[2][2];
        recovered_laplacian[0] = 2 * (polynomial_coefficients[6][0] + polynomial_coefficients[7][0] + polynomial_coefficients[8][0]);
        recovered_laplacian[1] = 2 * (polynomial_coefficients[6][1] + polynomial_coefficients[7][1] + polynomial_coefficients[8][1]);
        recovered_laplacian[2] = 2 * (polynomial_coefficients[6][2] + polynomial_coefficients[7][2] + polynomial_coefficients[8][2]);
    }
    AddTimeDerivative(r_model_part, material_derivative_container);
    mCalculatingGradientAndLaplacian = false;
    std::cout << "Exiting RecoverSuperconvergentMatDerivAndLaplacian" << std::endl;
}
//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
template <std::size_t TDim>
void DerivativeRecovery<TDim>::CalculateVectorLaplacian(ModelPart& r_model_part, Variable<array_1d<double, 3> >& vector_container, Variable<array_1d<double, 3> >& laplacian_container)
{
    KRATOS_INFO("SwimmingDEM") << "Constructing the Laplacian by derivating nodal averages..." << std::endl;
    std::map <std::size_t, unsigned int> id_to_position;
    unsigned int entry = 0;
    for (NodeIteratorType inode = r_model_part.NodesBegin(); inode != r_model_part.NodesEnd(); ++inode){
        noalias(inode->FastGetSolutionStepValue(laplacian_container)) = ZeroVector(3);
        id_to_position[inode->Id()] = entry;
        ++entry;
    }
    std::vector<array_1d <double, 3> > laplacians;
    laplacians.resize(entry);
    std::fill(laplacians.begin(), laplacians.end(), ZeroVector(3));
    array_1d <double, 3> grad = ZeroVector(3);
    array_1d <double, TDim + 1 > elemental_values;
    array_1d <double, TDim + 1 > N; // shape functions vector
    BoundedMatrix<double, TDim + 1, TDim> DN_DX;
    BoundedMatrix<double, TDim + 1, TDim> elemental_vectors; // They carry the nodal gradients of the corresponding component v_j
    const double nodal_area_share = 1.0 / static_cast<double>(TDim + 1);
    for (unsigned int j = 0; j < TDim; ++j){ // for each component of the original vector value
        // for each element, constructing the gradient contribution (to its nodes) of the component v_j and storing it in laplacian_container
        for (ModelPart::ElementIterator ielem = r_model_part.ElementsBegin(); ielem != r_model_part.ElementsEnd(); ++ielem){
            // computing the shape function derivatives
            Geometry<Node >& geom = ielem->GetGeometry();
            double Volume;
            GeometryUtils::CalculateGeometryData(geom, DN_DX, N, Volume);
            for (unsigned int i = 0; i < TDim + 1; ++i){
                elemental_values[i] = geom[i].FastGetSolutionStepValue(vector_container)[j];
            }
            array_1d <double, TDim> grad_aux = prod(trans(DN_DX), elemental_values); // its dimension may be 2
            for (unsigned int i = 0; i < TDim; ++i){
                grad[i] = grad_aux[i];
            }
            double nodal_area = Volume * nodal_area_share;
            grad *= nodal_area;
            for (unsigned int i = 0; i < TDim + 1; ++i){
                geom[i].FastGetSolutionStepValue(laplacian_container) += grad; // we use laplacian_container to store the gradient of one component at a time
            }
        }
        // normalizing the contributions to the gradient
        for (NodeIteratorType inode = r_model_part.NodesBegin(); inode != r_model_part.NodesEnd(); ++inode){
            inode->FastGetSolutionStepValue(laplacian_container) /= inode->FastGetSolutionStepValue(NODAL_AREA);
        }
        // for each element, constructing the divergence contribution (to its nodes) of the gradient of component v_j
        for (ModelPart::ElementIterator ielem = r_model_part.ElementsBegin(); ielem != r_model_part.ElementsEnd(); ++ielem){
            Geometry<Node >& geom = ielem->GetGeometry();
            double Volume;
            GeometryUtils::CalculateGeometryData(geom, DN_DX, N, Volume);
            for (unsigned int i = 0; i < TDim + 1; ++i){
                for (unsigned int k = 0; k < TDim; ++k){
                    elemental_vectors(i, k) = geom[i].FastGetSolutionStepValue(laplacian_container)[k]; // it is actually the gradient of component v_j
                }
            }
            BoundedMatrix<double, TDim, TDim> grad_aux = prod(trans(DN_DX), elemental_vectors); // its dimension may be 2
            double divergence_of_vi = 0.0;
            for (unsigned int k = 0; k < TDim; ++k){ // the divergence is the trace of the gradient
                divergence_of_vi += grad_aux(k, k);
            }
            double nodal_area = Volume / static_cast<double>(TDim + 1);
            divergence_of_vi *= nodal_area;
            for (unsigned int i = 0; i < TDim + 1; ++i){
                laplacians[id_to_position[geom[i].Id()]][j] += divergence_of_vi; // adding the contribution of the elemental divergence to each of its nodes
            }
        }
        // clearing the values stored in laplacian_container for the next component
        for (NodeIteratorType inode = r_model_part.NodesBegin(); inode != r_model_part.NodesEnd(); ++inode){
            array_1d <double, 3>& current_gradient_of_vi = inode->FastGetSolutionStepValue(laplacian_container);
            noalias(current_gradient_of_vi) = ZeroVector(3);
        }
    } // for each component of the vector value
    for (NodeIteratorType inode = r_model_part.NodesBegin(); inode != r_model_part.NodesEnd(); ++inode){
        array_1d <double, 3>& stored_laplacian = laplacians[id_to_position[inode->Id()]];
        array_1d <double, 3>& laplacian = inode->FastGetSolutionStepValue(laplacian_container);
        noalias(laplacian) = stored_laplacian / inode->FastGetSolutionStepValue(NODAL_AREA);
    }
    KRATOS_INFO("SwimmingDEM") << "Finished constructing the Laplacian by derivating nodal averages..." << std::endl;
}
//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
template <std::size_t TDim>
void DerivativeRecovery<TDim>::CalculateVelocityLaplacianRate(ModelPart& r_model_part)
{
    double delta_t_inv = 1.0 / r_model_part.GetProcessInfo()[DELTA_TIME];
    DenseVector<unsigned int> nodes_partition;
    OpenMPUtils::CreatePartition(ParallelUtilities::GetNumThreads(), r_model_part.Nodes().size(), nodes_partition);

    #pragma omp parallel for
    for (int k = 0; k < ParallelUtilities::GetNumThreads(); ++k){
        NodesArrayType& pNodes = r_model_part.GetCommunicator().LocalMesh().Nodes();
        NodeIteratorType node_begin = pNodes.ptr_begin() + nodes_partition[k];
        NodeIteratorType node_end   = pNodes.ptr_begin() + nodes_partition[k + 1];
        for (ModelPart::NodesContainerType::iterator inode = node_begin; inode != node_end; ++inode){
            array_1d <double, 3>& laplacian_rate = inode->FastGetSolutionStepValue(VELOCITY_LAPLACIAN_RATE);
            array_1d <double, 3>& laplacian      = inode->FastGetSolutionStepValue(VELOCITY_LAPLACIAN);
            array_1d <double, 3>& old_laplacian  = inode->FastGetSolutionStepValue(VELOCITY_LAPLACIAN, 1);
            noalias(laplacian_rate) = delta_t_inv * (laplacian - old_laplacian);
        }
    }
}
//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
//                  P R I V A T E   M E T H O D S
//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
template <std::size_t TDim>
bool DerivativeRecovery<TDim>::IsEdgeNode(Geometry<Node>::GeometriesArrayType& edge, Node::Pointer& p_node)
{
    return false;
}
//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
template <std::size_t TDim>
void DerivativeRecovery<TDim>::SetEdgeNodesAndWeights(ModelPart& r_model_part)
{
    // If quadratic elements are used in the superconvergent recovery,
    // then we have to be careful about the edge nodes. For an edge
    // node, its recovered gradient is the weighted mean of the recovery
    // of its vertex nodes.

    std::vector<int> ids;
    std::vector<int>::iterator it;
    for(ElementIteratorType ielem = r_model_part.ElementsBegin(); ielem != r_model_part.ElementsEnd(); ++ielem)
    {
        Geometry<Node >& geom = ielem->GetGeometry();
        Geometry<Node>::GeometriesArrayType edges = geom.GenerateEdges();

        for(unsigned int n = 0; n < edges.size(); n++)
        {
            Geometry<Node>& edge = edges[n];
            if (edge.size() == 2)
                continue;

            Node::NodeType& edge_node = edge[edge.size() - 1];  // Last node is the edge node
            edge_node.SetValue(IS_EDGE_NODE, true);
        }
    }
}
//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
template <std::size_t TDim>
void DerivativeRecovery<TDim>::SetNeighboursAndWeights(ModelPart& r_model_part, const unsigned int& ord)
{
    // Finding elements concurrent to each node. The nodes of these elements will form the initial cloud of points
    FindNodalNeighboursProcess neighbour_finder(r_model_part);
    neighbour_finder.Execute();
    const unsigned int n_max_iterations = 100;
    unsigned int i = 0;

    for (NodeIteratorType inode = r_model_part.NodesBegin(); inode != r_model_part.NodesEnd(); ++inode){
        // The weights only make sense if the node is a vertex node
        if (inode->GetValue(IS_EDGE_NODE)) continue;

        bool the_cloud_of_neighbours_is_successful = SetInitialNeighboursAndWeights(r_model_part, *(inode.base()), ord);
        GlobalPointersVector<Node >& neigh_nodes = inode->GetValue(NEIGHBOUR_NODES);
        unsigned int iteration = 0;
        while (!the_cloud_of_neighbours_is_successful && iteration < n_max_iterations){
            the_cloud_of_neighbours_is_successful = SetNeighboursAndWeights(r_model_part, *(inode.base()), ord);
            ++iteration;
        }
        if (iteration >= n_max_iterations){ // giving up on this method, settling for the default method
            mSomeCloudsDontWork = true;
            neigh_nodes.clear();
            inode->FastGetSolutionStepValue(NODAL_WEIGHTS).clear();
            KRATOS_WARNING("SwimmingDEM") << "Warning!, for the node with id " << inode->Id() << " it has not been possible to form an adequate cloud of neighbours" << std::endl;
            KRATOS_WARNING("SwimmingDEM") << "for the gradient recovery. A lower accuracy method has been employed for this node." << std::endl;
        }
        ++i;
    }
}
//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
template <std::size_t TDim>
void DerivativeRecovery<TDim>::SetNeighboursAndWeightsForTheLaplacian(ModelPart& r_model_part, const unsigned int& ord)
{
    // Finding elements concurrent to each node. The nodes of these elements will form the initial cloud of points
    FindNodalNeighboursProcess neighbour_finder(r_model_part);
    neighbour_finder.Execute();
    const unsigned int n_max_iterations = 100;
    unsigned int i = 0;
    for (NodeIteratorType inode = r_model_part.NodesBegin(); inode != r_model_part.NodesEnd(); ++inode){
        bool the_cloud_of_neighbours_is_successful = SetInitialNeighboursAndWeights(r_model_part, *(inode.base()), ord);
        GlobalPointersVector<Node >& neigh_nodes = inode->GetValue(NEIGHBOUR_NODES);
        unsigned int iteration = 0;
        while (!the_cloud_of_neighbours_is_successful && iteration < n_max_iterations){
            the_cloud_of_neighbours_is_successful = SetNeighboursAndWeights(r_model_part, *(inode.base()), ord);
            ++iteration;
        }
        if (iteration >= n_max_iterations){ // giving up on this method, settling for the default method
            mSomeCloudsDontWork = true;
            neigh_nodes.clear();
            inode->FastGetSolutionStepValue(NODAL_WEIGHTS).clear();
            KRATOS_WARNING("SwimmingDEM") << "Warning!, for the node with id " << inode->Id() << " it has not been possible to form an adequate cloud of neighbours" << std::endl;
            KRATOS_WARNING("SwimmingDEM") << "for the gradient recovery. A lower accuracy method has been employed for this node." << std::endl;
        }
        ++i;
    }
}
//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
template <std::size_t TDim>
void DerivativeRecovery<TDim>::OrderByDistance(Node::Pointer &p_node, GlobalPointersVector<Node >& neigh_nodes)
{
    const unsigned int n_nodes = neigh_nodes.size();
    std::vector<double> distances_squared;
    distances_squared.resize(n_nodes);
    const array_1d <double, 3>& origin = p_node->Coordinates();
    for (unsigned int i = 0; i < n_nodes; ++i){
        const array_1d <double, 3> rel_coordinates = (neigh_nodes[i] - origin);
        distances_squared[i] = DEM_INNER_PRODUCT_3(rel_coordinates, rel_coordinates);
    }
    std::vector <std::pair<unsigned int, double> > ordering;
    ordering.resize(n_nodes);
    for (unsigned int i = 0; i < n_nodes; ++i){
        ordering[i] = std::make_pair(i, distances_squared[i]);
    }
    std::sort(ordering.begin(), ordering.end(), IsCloser());
    GlobalPointersVector<Node > ordered_neighbours;
    for (unsigned int i = 0; i < n_nodes; ++i){
        Node::WeakPointer& p_neigh = neigh_nodes(ordering[i].first);
        ordered_neighbours.push_back(p_neigh);
    }
    ordered_neighbours.swap(neigh_nodes);
}
//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
template <std::size_t TDim>
bool DerivativeRecovery<TDim>::SetInitialNeighboursAndWeights(ModelPart& r_model_part, Node::Pointer &p_node, const unsigned int& ord)
{
    GlobalPointersVector<Element>& neigh_elems = p_node->GetValue(NEIGHBOUR_ELEMENTS);
    GlobalPointersVector<Node >& neigh_nodes = p_node->GetValue(NEIGHBOUR_NODES);
    std::map<std::size_t, std::size_t> ids; // map to keep track of all different ids corresponding to already added neighbours to avoid repetition
    ids[p_node->Id()] = p_node->Id();
    for (unsigned int i_el = 0; i_el < neigh_elems.size(); ++i_el){
        Geometry<Node >& geom = neigh_elems[i_el].GetGeometry();
        unsigned int jj = 0; // index of the node in geom corresponding to neighbour neigh_elems[i_el]
        if (geom[jj].Id() == p_node->Id()){ // skipping itself
            ++jj;
        }
        for (unsigned int j = 0; j < TDim; ++j){
            Node::Pointer p_neigh = geom(jj);
            if (ids.find(p_neigh->Id()) == ids.end()){
                neigh_nodes.push_back(p_neigh);
                ids[p_neigh->Id()] = p_neigh->Id();
            }
        }
    }
    OrderByDistance(p_node, neigh_nodes);
    unsigned int min_neighbours = Factorial(TDim + ord) / (Factorial(TDim) * Factorial(ord));
    if (neigh_nodes.size() < min_neighbours){ // Not worthwhile checking, since there are 10 independent coefficients to be determined
        return false;
    }
    else {
        return SetWeightsAndRunLeastSquaresTest(r_model_part, p_node, ord);
    }
}
//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
template <std::size_t TDim>
bool DerivativeRecovery<TDim>::SetNeighboursAndWeights(ModelPart& r_model_part, Node::Pointer& p_node, const unsigned int& ord)
{
    GlobalPointersVector<Node >& neigh_nodes = p_node->GetValue(NEIGHBOUR_NODES);
    const unsigned int node_increase_per_neighbour = 1;
    const unsigned int node_increase_overall = 1;
    std::map<std::size_t, std::size_t> ids;
    ids[p_node->Id()] = p_node->Id();
    for (unsigned int i = 0; i < (unsigned int)neigh_nodes.size(); ++i){
        auto p_neigh = neigh_nodes(i);
        ids[p_neigh->Id()] = p_neigh->Id();
    }
    const unsigned int n_neigh = neigh_nodes.size();
    for (unsigned int i = 0; i < n_neigh; ++i){
        auto p_neigh = neigh_nodes(i);
        GlobalPointersVector<Node >& neigh_neigh_nodes = p_neigh->GetValue(NEIGHBOUR_NODES);
        unsigned int n_new_nodes = 0;
        for (unsigned int j = 0; j < (unsigned int)neigh_neigh_nodes.size(); ++j){
            auto p_neigh_neigh = neigh_neigh_nodes(j);
            if (ids.find(p_neigh_neigh->Id()) == ids.end()){
                neigh_nodes.push_back(p_neigh_neigh);
                ids[p_neigh_neigh->Id()] = p_neigh_neigh->Id();
                ++n_new_nodes;
            }
            if (n_new_nodes >= node_increase_per_neighbour){
                break;
            }
        }
    }
    OrderByDistance(p_node, neigh_nodes);
    const unsigned int new_size = std::min(n_neigh + node_increase_overall, (unsigned int)neigh_nodes.size());
    neigh_nodes.resize(new_size); // keeping only nearest nodes

    unsigned int min_neighbours = Factorial(TDim + ord) / (Factorial(TDim) * Factorial(ord));
    if (neigh_nodes.size() < min_neighbours){ // it is not worthwhile checking, since there are 10 independent coefficients to be determined
        return false;
    }
    else {
        return SetWeightsAndRunLeastSquaresTest(r_model_part, p_node, ord);
    }
}
//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
template <std::size_t TDim>
double DerivativeRecovery<TDim>::SecondDegreeTestPolynomial(const array_1d <double, 3>& coordinates)
{
    const double x = coordinates[0];
    const double y = coordinates[1];
    double z = 0.0;
    if (TDim == 3) z = coordinates[2];
    return 1.0 + x + y + z + x * y + x * z + y * z + x * x + y * y + z * z;
}
//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
template <std::size_t TDim>
double DerivativeRecovery<TDim>::SecondDegreeGenericPolynomial(DenseMatrix<double> C, const array_1d <double, 3>& coordinates)
{
    const double x = coordinates[0];
    const double y = coordinates[1];
    if (TDim == 3){
        const double z = coordinates[2];
        return C(0,0) + C(1,0) * x + C(2,0) * y + C(3,0) * z + C(4,0) * x * y + C(5,0) * x * z + C(6,0) * y * z + C(7,0) * x * x + C(8,0) * y * y + C(9,0) * z * z;
    }
    else{
        return C(0,0) + C(1,0) * x + C(2,0) * y + C(3,0) * x * y + C(4,0) * x * x + C(5,0) * y * y;
    }

}
//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
template <std::size_t TDim>
inline int DerivativeRecovery<TDim>:: Factorial(const unsigned int n){
    if (n == 0){
        return 1;
    }
    unsigned int k = n;
    for (unsigned int i = n - 1; i > 0; --i){
        k *= i;
    }
    return k;
}
//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
template <std::size_t TDim>
bool DerivativeRecovery<TDim>::SetWeightsAndRunLeastSquaresTest(ModelPart& r_model_part, Node::Pointer& p_node, const unsigned int& ord)
{
    unsigned int n_poly_terms = Factorial(TDim + ord) / (Factorial(ord) * Factorial(TDim)); // 2 is the polynomial order
    //KRATOS_ERROR_IF(TDim == 2) << "Gradient recovery not implemented yet in 2D!)" << std::endl;
    GlobalPointersVector<Node >& neigh_nodes = p_node->GetValue(NEIGHBOUR_NODES);
    unsigned int n_nodal_neighs = (unsigned int)neigh_nodes.size();
    const double h_inv = 1.0 / CalculateTheMaximumDistanceToNeighbours(p_node); // we use it as a scaling parameter to improve stability
    const array_1d <double, 3> origin = p_node->Coordinates();
    DenseMatrix<double> TestNodalValues(n_nodal_neighs, 1);
    DenseMatrix<double> A(n_nodal_neighs, n_poly_terms);

    // Write the A matrix from Zhang (2005), each row being (1 xi_i eta_i ... xi_i^k eta_i^l ... eta_i^(k+1))
    // where the i-th row correspoinds to the ith neighbor (in relative coords)
    // unsigned int min_neighbours = Factorial(TDim + ord) / (Factorial(TDim) * Factorial(ord));
    // std::cout << "Setting matrix A for node " << p_node->Id() << std::endl;
    for (unsigned int i = 0; i < n_nodal_neighs; ++i){
        A(i, 0) = 1.0;
        if constexpr (TDim == 3){
            Node& neigh = neigh_nodes[i];
            const array_1d <double, 3> rel_coordinates = (neigh.Coordinates() - origin) * h_inv;
            TestNodalValues(i, 0) = SecondDegreeTestPolynomial(rel_coordinates);

            for (unsigned int k = 1; k < n_poly_terms; ++k)
            {
                // Terms of order 1
                if (ord >= 1)
                {
                    if (k < 4)
                    {
                        A(i, k) = rel_coordinates[k - 1]; // x, y and z
                    }
                }

                // Terms of order 2
                if (ord >= 2)
                {
                    if (k == 4)
                    {
                        A(i, k) = rel_coordinates[0] * rel_coordinates[1]; // x y
                    }
                    else if (k == 5)
                    {
                        A(i, k) = rel_coordinates[0] * rel_coordinates[2]; // x z
                    }
                    else if (k == 6)
                    {
                        A(i, k) = rel_coordinates[1] * rel_coordinates[2]; // y z
                    }
                    else if (k == 7)
                    {
                        A(i, k) = rel_coordinates[0] * rel_coordinates[0]; // x^2
                    }
                    else if (k == 8)
                    {
                        A(i, k) = rel_coordinates[1] * rel_coordinates[1]; // y^2
                    }
                    else if (k == 9)
                    {
                        A(i, k) = rel_coordinates[2] * rel_coordinates[2]; // z^2
                    }
                }

                // Terms of order 3
                if (ord == 3)
                {
                    if (k == 10)
                    {
                        A(i, k) == rel_coordinates[0] * rel_coordinates[1] * rel_coordinates[2];  // x y z
                    }
                    else if (k == 11)
                    {
                        A(i, k) = rel_coordinates[0] * rel_coordinates[0] * rel_coordinates[1];  // x^2 y
                    }
                    else if (k == 12)
                    {
                        A(i, k) = rel_coordinates[0] * rel_coordinates[0] * rel_coordinates[2];  // x^2 z
                    }
                    else if (k == 13)
                    {
                        A(i, k) = rel_coordinates[0] * rel_coordinates[1] * rel_coordinates[1];  // x y^2
                    }
                    else if (k == 14)
                    {
                        A(i, k) = rel_coordinates[0] * rel_coordinates[2] * rel_coordinates[2];  // x z^2
                    }
                    else if (k == 15)
                    {
                        A(i, k) = rel_coordinates[1] * rel_coordinates[1] * rel_coordinates[2];  // y^2 z
                    }
                    else if (k == 16)
                    {
                        A(i, k) = rel_coordinates[1] * rel_coordinates[2] * rel_coordinates[2];  // y z^2
                    }
                    else if (k == 17)
                    {
                        A(i, k) = rel_coordinates[0] * rel_coordinates[0] * rel_coordinates[0];  // x^3
                    }
                    else if (k == 18)
                    {
                        A(i, k) = rel_coordinates[1] * rel_coordinates[1] * rel_coordinates[1]; // y^3
                    }
                    else if (k == 19)
                    {
                        A(i, k) = rel_coordinates[2] * rel_coordinates[2] * rel_coordinates[2];  // z^3
                    }
                }

                if (ord > 3)
                {
                    KRATOS_ERROR << "Superconvergent gradient recovery of order " << ord << "not implemented yet in 3D!" << std::endl;
                }
            }
        }
        if constexpr (TDim == 2){
            Node& neigh = neigh_nodes[i];
            const array_1d <double, 3> rel_coordinates = (neigh.Coordinates() - origin) * h_inv;
            TestNodalValues(i, 0) = SecondDegreeTestPolynomial(rel_coordinates);
            for (unsigned int d = 1; d < 6; ++d){
                if (d < 3){
                    A(i, d) = rel_coordinates[d - 1];
                }
                else if (d == 3){
                    A(i, d) = rel_coordinates[0] * rel_coordinates[1];
                }
                else {
                    A(i, d) = rel_coordinates[d - 4] * rel_coordinates[d - 4];
                }
            }
        }
        //else {
        //    KRATOS_ERROR << "Gradient recovery not implemented yet in 2D!)" << std::endl;
        //}
    }
    // std::cout << "Matrix A with size (" << A.size1() << ", " << A.size2() << ") succesfully computed for node " << p_node->Id() << std::endl;

    DenseMatrix<double>AtransA(n_poly_terms, n_poly_terms);
    // std::cout << "Computing A^T A for node " << p_node->Id() << std::endl;
    noalias(AtransA) = prod(trans(A), A);
    // std::cout << "Computed A^T A with size (" << AtransA.size1() << ", " << AtransA.size2() << ") for node " << p_node->Id() << std::endl;
    if (std::abs(MathUtils<double>::Det(AtransA)) < 0.01){
        // std::cout << "Determinant of AtransA not ok for node " << p_node->Id() << std::endl;
        return false;
    }
    // std::cout << "Determinant of AtransA ok for node " << p_node->Id() << std::endl;
        
    // unsigned int n_relevant_terms = 0;
    // if (mCalculatingTheGradient){
    //     n_relevant_terms = TDim;
    // }
    // else if (mCalculatingTheLaplacian){
    //     n_relevant_terms = n_poly_terms - (TDim + 1);
    // }
    // else {
    //     n_relevant_terms = n_poly_terms - 1;
    // }

    // std::vector<unsigned int> relevant_terms;
    // relevant_terms.resize(n_relevant_terms);
    // double normalization =  1.0;
    // if (mCalculatingTheGradient){
    //     normalization = h_inv;
    //     relevant_terms[0] = 1;
    //     relevant_terms[1] = 2;
    //     if (TDim == 3) relevant_terms[2] = 3;
    // }
    // else if (mCalculatingTheLaplacian){
    //     normalization = h_inv * h_inv;
    //     relevant_terms[0] = 4;
    //     relevant_terms[1] = 5;
    //     relevant_terms[2] = 6;
    //     relevant_terms[3] = 7;
    //     relevant_terms[4] = 8;
    //     relevant_terms[5] = 9;
    // }
    // else {
    //     for (unsigned int i = 1; i < n_poly_terms; ++i){
    //         relevant_terms[i - 1] = i;
    //     }
    // }

    double normalization = h_inv;  // Review this!
    unsigned int n_relevant_terms = TDim;
    std::vector<unsigned int> relevant_terms(n_relevant_terms);
    for (unsigned int i = 0; i < n_relevant_terms; ++i){
        relevant_terms[i] = i + 1;
    }

    p_node->SetLock();
    p_node->SetValue(NODAL_WEIGHTS, ZeroVector(n_relevant_terms * n_nodal_neighs));
    Vector& nodal_weights = p_node->GetValue(NODAL_WEIGHTS);
    p_node->UnSetLock();

    const DenseMatrix<double> AtransAinv = mMyCustomFunctions.Inverse(AtransA);
    DenseMatrix<double>AtransAinvAtrans(n_poly_terms, n_nodal_neighs);
    noalias(AtransAinvAtrans) = prod(AtransAinv, trans(A));

    // std::cout << "AtransAinvAtrans computed for node " << p_node->Id() << std::endl;

    for (unsigned int i = 0; i < n_nodal_neighs; ++i){
        for (unsigned int d = 0; d < n_relevant_terms; ++d){
            if (mCalculatingGradientAndLaplacian){
                if (d > 2){
                    normalization = h_inv * h_inv;
                }
                else {
                    normalization = h_inv;
                }
            }
            nodal_weights(n_relevant_terms * i + d) = AtransAinvAtrans(relevant_terms[d], i) * normalization;
        }
    }
    // std::cout << "Nodal weights have been set for node " << p_node->Id() << std::endl;

    DenseMatrix<double> C = ZeroMatrix(n_nodal_neighs, 1);
    C = prod(AtransAinvAtrans, TestNodalValues);
    double abs_difference = 0.0;
    for (unsigned int i = 0; i < n_nodal_neighs; ++i){
        const array_1d <double, 3>& rel_coordinates = (neigh_nodes[i].Coordinates() - origin) * h_inv;
        abs_difference += std::abs(SecondDegreeGenericPolynomial(C, rel_coordinates) - SecondDegreeTestPolynomial(rel_coordinates));
    }
//        abs_difference = abs(C(7, 0) - 1.0) + abs(C(8, 0) - 1.0) + abs(C(8, 0) - 1.0);
    const double tolerance = 0.001; // recommended by E Ortega
    return (abs_difference > tolerance? false : true);
}
//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
template <std::size_t TDim>
unsigned int DerivativeRecovery<TDim>::GetNumberOfUniqueNeighbours(const int my_id, const GlobalPointersVector<Element>& my_neighbour_elements)
{
    std::vector<int> ids;
    ids.push_back(my_id);
    for (unsigned int i_el = 0; i_el < my_neighbour_elements.size(); ++i_el){
        const Geometry<Node >& geom = my_neighbour_elements[i_el].GetGeometry();
        for (unsigned int jj = 0; jj < TDim + 1; ++jj){
            int id = (int)geom[jj].Id();
            std::vector<int>::iterator it;
            it = find(ids.begin(), ids.end(), id);
            if (it >= ids.end()){
                ids.push_back(id);
            }
        }
    }
    return (int)ids.size();
}
//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
template <std::size_t TDim>
double DerivativeRecovery<TDim>::CalculateTheMaximumDistanceToNeighbours(Node::Pointer& p_node)
{
    double max_distance_yet = 0.0;
    const array_1d <double, 3>& coors = p_node->Coordinates();
    GlobalPointersVector<Node >& neigh_nodes = p_node->GetValue(NEIGHBOUR_NODES);
    for (unsigned int i = 0; i < (unsigned int)neigh_nodes.size(); ++i){
        array_1d <double, 3> delta = neigh_nodes[i].Coordinates() - coors;
        double distance_2 = DEM_INNER_PRODUCT_3(delta, delta);
        max_distance_yet = max_distance_yet > distance_2 ? max_distance_yet : distance_2;
    }
    return std::sqrt(max_distance_yet);
}
//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
template <std::size_t TDim>
double DerivativeRecovery<TDim>::CalculateTheMaximumEdgeLength(ModelPart& r_model_part)
{
    double max_distance_yet = 0.0;
    for (ModelPart::ElementIterator ielem = r_model_part.ElementsBegin(); ielem != r_model_part.ElementsEnd(); ++ielem){
        Geometry<Node >& geom = ielem->GetGeometry();
        unsigned int n_nodes = geom.size();
        for (unsigned int k = 1; k < n_nodes - 1; ++k){
            for (unsigned int i = k; i < n_nodes; ++i){
                array_1d <double, 3> delta_i = geom[k - 1] - geom[i];
                double distance_2 = DEM_INNER_PRODUCT_3(delta_i, delta_i);
                max_distance_yet = max_distance_yet > distance_2 ? max_distance_yet : distance_2;
            }
        }
    }
    return std::sqrt(max_distance_yet);
}
//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
template <std::size_t TDim>
double DerivativeRecovery<TDim>::CalculateTheMinumumEdgeLength(ModelPart& r_model_part)
{
    double min_distance_yet = 0.0;
    bool first_node = true;
    for (ModelPart::ElementIterator ielem = r_model_part.ElementsBegin(); ielem != r_model_part.ElementsEnd(); ++ielem){
        Geometry<Node >& geom = ielem->GetGeometry();
        if (first_node){ // assign the distance (squared) between any two nodes to min_distance_yet
            array_1d <double, 3> delta = geom[0] - geom[1];
            double distance_2 = DEM_INNER_PRODUCT_3(delta, delta);
            min_distance_yet = distance_2;
        }
        unsigned int n_nodes = geom.size();
        // unsigned int n_nodes = static_cast<unsigned int>(TDim + 1);
        for (unsigned int k = 1; k < n_nodes - 1; ++k){
            for (unsigned int i = k; i < n_nodes; ++i){
                array_1d <double, 3> delta_i = geom[k - 1] - geom[i];
                double distance_2 = DEM_INNER_PRODUCT_3(delta_i, delta_i);
                min_distance_yet = min_distance_yet < distance_2 ? min_distance_yet : distance_2;
            }
        }
    }
    return std::sqrt(min_distance_yet);
}
//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
template <std::size_t TDim>
void DerivativeRecovery<TDim>::ComputeGradientForVertexNode(Node::Pointer& inode, Vector& gradient, Variable<double>& scalar_container)
{
    GlobalPointersVector<Node >& neigh_nodes = inode->GetValue(NEIGHBOUR_NODES);
    unsigned int n_neigh = neigh_nodes.size();
    if (!n_neigh){ // we keep the defualt value
        return;
    }

    const Vector& nodal_weights = inode->GetValue(NODAL_WEIGHTS);
    for (unsigned int i_neigh = 0; i_neigh < n_neigh; ++i_neigh){
        const double& neigh_nodal_value = neigh_nodes[i_neigh].FastGetSolutionStepValue(scalar_container);
        for (unsigned int d = 0; d < TDim; ++d){
            gradient[d] += nodal_weights[TDim * i_neigh + d] * neigh_nodal_value;
        }
    }
}
//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
template <std::size_t TDim>
void DerivativeRecovery<TDim>::ComputeGradientForEdgeNode(Node::Pointer& inode, Vector& gradient, Variable<double>& scalar_container)
{
    // GlobalPointersVector<Node >& neigh_nodes = inode->GetValue(CONTIGUOUS_NODES);
    // unsigned int n_neigh = neigh_nodes.size();
    // if (!n_neigh){ // we keep the defualt value
    //     return;
    // }
    // Vector dr = neigh_nodes[1].Coordinates() - neigh_nodes[0].Coordinates();


    // // Compute the coefficients matrix and the P vector
    // DenseMatrix<double> P(TDim, n_poly_terms);
    // DenseMatrix<double> coefficientsMatrix(n_poly_terms, n_neigh);
    // ComputeDerivativeMonomialsVector(z_i, *(inode.base()), ord, P);

    // bool is_matrix_succesfully_computed = true;
    // ComputeCoefficientsMatrix(*(inode.base()), coefficientsMatrix, ord, is_matrix_succesfully_computed);
    // if(!is_matrix_succesfully_computed)
    // {
    //     KRATOS_ERROR << "Unable to compute the recovery for node with ID = " << std::endl;
    // }

    // array_1d<double, 3> z_i = inode->Coordinates();
    // bool is_matrix_succesfully_computed = true;
    // // std::cout << "The node " << inode->Coordinates() << " is not edge node! It has " << n_neigh << " neighbours!" << std::endl;

    // // Evaluate the gradient centered at z_i at z = z_i
    // for (unsigned int i_neigh = 0; i_neigh < n_neigh; ++i_neigh)
    // {
    //     const double& neigh_nodal_value = neigh_nodes[i_neigh].FastGetSolutionStepValue(scalar_container);
    //     for (unsigned int d = 0; d < TDim; ++d)
    //     {
    //         double weight = 0.0;
    //         for (unsigned int n = 0; n < n_poly_terms; n++)
    //         {
    //             weight += P(d, n) * coefficientsMatrix(n, i_neigh);
    //             // std::cout << "    - P(" << d << ", " << n << ") = " << P(d, n) << ", C(" << n << ", " << i_neigh << ") = " << coefficientsMatrix(i_neigh, n) << std::endl;
    //         }
    //         gradient[d] += weight * neigh_nodal_value;
    //     }
    // }
}
//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
template <std::size_t TDim>
template <class TScalarVariable>
void DerivativeRecovery<TDim>::RecoverSuperconvergentGradientAlt(ModelPart& r_model_part, TScalarVariable& scalar_container, Variable<array_1d<double, 3> >& gradient_container, unsigned int& ord)
{
    // Construct the clouds
    mCalculatingTheGradient = true;
    if (mFirstGradientRecovery){
        KRATOS_INFO("SwimmingDEM") << "Constructing first-step neighbour clouds for gradient recovery..." << std::endl;
        ClassifyEdgeNodes(r_model_part);
        SetNeighboursAndWeights(r_model_part, ord);
        mFirstGradientRecovery = false;
        KRATOS_INFO("SwimmingDEM") << "Finished constructing neighbour clouds for gradient recovery." << std::endl;
    }
    KRATOS_ERROR_IF(mSomeCloudsDontWork) << "Unable to construct neighbour cloudes for super-convergent gradient recovery" << std::endl;

    // Number of terms of the polynomial
    // unsigned int n_poly_terms = Factorial(TDim + ord) / (Factorial(ord) * Factorial(TDim)); // 2 is the polynomial order

    // Compute the coefficients for each node (Zhang, 2005)
    for (NodeIteratorType inode = r_model_part.NodesBegin(); inode != r_model_part.NodesEnd(); ++inode)
    {

        Vector gradient = ZeroVector(3);
        if (inode->GetValue(IS_EDGE_NODE))
        {
            // GlobalPointersVector<Node> &edge_nodes = inode->GetValue(EDGE_NODES);
        }
        else
        {
            // GlobalPointersVector<Node> &neigh_nodes = inode->GetValue(NEIGHBOUR_NODES);
            // unsigned int n_neigh = neigh_nodes.size();

            // Compute the coefficients matrix and the P vector
            // DenseMatrix<double> P(TDim, n_poly_terms);
            // DenseMatrix<double> coefficientsMatrix(n_poly_terms, n_neigh);

            // array_1d<double, 3> z_i = inode->Coordinates();
            // ComputeDerivativeMonomialsVector(z_i, *(inode.base()), ord, P);
            // ComputeDerivativeMonomialsVector(*(inode.base()), ord, P);
            // bool is_matrix_succesfully_computed = true;
            // ComputeCoefficientsMatrix(*(inode.base()), coefficientsMatrix, ord, is_matrix_succesfully_computed);
            // // std::cout << "The node " << inode->Coordinates() << " is not edge node! It has " << n_neigh << " neighbours!" << std::endl;

            // // Evaluate the gradient centered at z_i at z = z_i
            // for (unsigned int i_neigh = 0; i_neigh < n_neigh; ++i_neigh)
            // {
            //     const double& neigh_nodal_value = neigh_nodes[i_neigh].FastGetSolutionStepValue(scalar_container);
            //     for (unsigned int d = 0; d < TDim; ++d)
            //     {
            //         double weight = 0.0;
            //         for (unsigned int n = 0; n < n_poly_terms; n++)
            //         {
            //             weight += P(d, n) * coefficientsMatrix(n, i_neigh);
            //             // std::cout << "    - P(" << d << ", " << n << ") = " << P(d, n) << ", C(" << n << ", " << i_neigh << ") = " << coefficientsMatrix(i_neigh, n) << std::endl;
            //         }
            //         gradient[d] += weight * neigh_nodal_value;
            //     }
            // }
            ComputeGradientForVertexNode(*(inode.base()), gradient, scalar_container);
        }

        // Copy gradient values to the nodal gradient container variable
        array_1d <double, 3>& recovered_gradient = inode->FastGetSolutionStepValue(gradient_container);
        for(unsigned int i = 0; i < TDim; i++) recovered_gradient[i] = gradient[i];
    }
}
//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
template <std::size_t TDim>
void DerivativeRecovery<TDim>::ClassifyEdgeNodes(ModelPart& r_model_part)
{
    // If quadratic elements are used in the superconvergent recovery,
    // then we have to be careful about the edge nodes. For an edge
    // node, its recovered gradient is the weighted mean of the recovery
    // of its vertex nodes.

    std::vector<int> ids;
    std::vector<int>::iterator it;
    for(ElementIteratorType ielem = r_model_part.ElementsBegin(); ielem != r_model_part.ElementsEnd(); ++ielem)
    {
        Geometry<Node >& geom = ielem->GetGeometry();
        Geometry<Node>::GeometriesArrayType edges = geom.GenerateEdges();

        for(unsigned int n = 0; n < edges.size(); n++)
        {
            Geometry<Node>& edge = edges[n];
            if (edge.size() == 2)
                continue;

            Node::NodeType& edge_node = edge[edge.size() - 1];  // Last node is the edge node
            edge_node.SetValue(IS_EDGE_NODE, true);
        }
    }
}
//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
template <std::size_t TDim>
void DerivativeRecovery<TDim>::ComputeCoefficientsMatrix(Node::Pointer& p_node, DenseMatrix<double>& CoeffsMatrix, const unsigned int& ord, bool& is_matrix_successfully_computed)
{
    unsigned int n_poly_terms = Factorial(TDim + ord) / (Factorial(ord) * Factorial(TDim)); // 2 is the polynomial order

    GlobalPointersVector<Node >& neigh_nodes = p_node->GetValue(NEIGHBOUR_NODES);
    unsigned int n_nodal_neighs = (unsigned int)neigh_nodes.size();
    const double h_inv = 1.0 / CalculateTheMaximumDistanceToNeighbours(p_node); // we use it as a scaling parameter to improve stability
    // const double h_inv = 1.0;
    const array_1d <double, 3> origin = p_node->Coordinates();
    DenseMatrix<double> A(n_nodal_neighs, n_poly_terms);

    // Write the A matrix from Zhang (2005), each row being (1 xi_i eta_i ... xi_i^k eta_i^l ... eta_i^(k+1))
    // where the i-th row correspoinds to the ith neighbor (in relative coords)
    // unsigned int min_neighbours = Factorial(TDim + ord) / (Factorial(TDim) * Factorial(ord));
    for (unsigned int i = 0; i < n_nodal_neighs; ++i){
        A(i, 0) = 1.0;
        if constexpr (TDim == 3){
            Node& neigh = neigh_nodes[i];
            const array_1d <double, 3> rel_coordinates = (neigh.Coordinates() - origin) * h_inv;
            for (unsigned int k = 1; k < n_poly_terms; ++k)
            {
                // Terms of order 1
                if (ord >= 1)
                {
                    if (k < 4)
                    {
                        A(i, k) = rel_coordinates[k - 1]; // x, y and z
                    }
                }

                // Terms of order 2
                if (ord >= 2)
                {
                    if (k == 4)
                    {
                        A(i, k) = rel_coordinates[0] * rel_coordinates[1]; // x y
                    }
                    else if (k == 5)
                    {
                        A(i, k) = rel_coordinates[0] * rel_coordinates[2]; // x z
                    }
                    else if (k == 6)
                    {
                        A(i, k) = rel_coordinates[1] * rel_coordinates[2]; // y z
                    }
                    else if (k == 7)
                    {
                        A(i, k) = rel_coordinates[0] * rel_coordinates[0]; // x^2
                    }
                    else if (k == 8)
                    {
                        A(i, k) = rel_coordinates[1] * rel_coordinates[1]; // y^2
                    }
                    else if (k == 9)
                    {
                        A(i, k) = rel_coordinates[2] * rel_coordinates[2]; // z^2
                    }
                }

                // Terms of order 3
                if (ord == 3)
                {
                    if (k == 10)
                    {
                        A(i, k) == rel_coordinates[0] * rel_coordinates[1] * rel_coordinates[2];  // x y z
                    }
                    else if (k == 11)
                    {
                        A(i, k) = rel_coordinates[0] * rel_coordinates[0] * rel_coordinates[1];  // x^2 y
                    }
                    else if (k == 12)
                    {
                        A(i, k) = rel_coordinates[0] * rel_coordinates[0] * rel_coordinates[2];  // x^2 z
                    }
                    else if (k == 13)
                    {
                        A(i, k) = rel_coordinates[0] * rel_coordinates[1] * rel_coordinates[1];  // x y^2
                    }
                    else if (k == 14)
                    {
                        A(i, k) = rel_coordinates[0] * rel_coordinates[2] * rel_coordinates[2];  // x z^2
                    }
                    else if (k == 15)
                    {
                        A(i, k) = rel_coordinates[1] * rel_coordinates[1] * rel_coordinates[2];  // y^2 z
                    }
                    else if (k == 16)
                    {
                        A(i, k) = rel_coordinates[1] * rel_coordinates[2] * rel_coordinates[2];  // y z^2
                    }
                    else if (k == 17)
                    {
                        A(i, k) = rel_coordinates[0] * rel_coordinates[0] * rel_coordinates[0];  // x^3
                    }
                    else if (k == 18)
                    {
                        A(i, k) = rel_coordinates[1] * rel_coordinates[1] * rel_coordinates[1]; // y^3
                    }
                    else if (k == 19)
                    {
                        A(i, k) = rel_coordinates[2] * rel_coordinates[2] * rel_coordinates[2];  // z^3
                    }
                }

                if (ord > 3)
                {
                    KRATOS_ERROR << "Superconvergent gradient recovery of order " << ord << "not implemented yet in 3D!" << std::endl;
                }
            }
        }
    }

    DenseMatrix<double>AtransA(n_poly_terms, n_poly_terms);
    noalias(AtransA) = prod(trans(A), A);
    if (std::abs(MathUtils<double>::Det(AtransA)) < 0.01){
        is_matrix_successfully_computed = false;
    }
    const DenseMatrix<double> AtransAinv = mMyCustomFunctions.Inverse(AtransA);

    // Assign the values to the coeffs matrix
    DenseMatrix<double> AtransAinvAtrans(n_poly_terms, n_nodal_neighs);
    noalias(AtransAinvAtrans) = prod(AtransAinv, trans(A));

    for (unsigned int i = 0; i < AtransAinvAtrans.size1(); i++)
    {
        for (unsigned int j = 0; j < AtransAinvAtrans.size2(); j++)
        {
            CoeffsMatrix(i, j) = AtransAinvAtrans(i, j) * h_inv;
        }
    }

    // Print matricies
    if(p_node->Id() == -1)
    {
        std::cout << "(new) A =" << std::endl;
        for (unsigned int i = 0; i < A.size1(); i++)
        {
            for (unsigned int j = 0; j < A.size2(); j++)
            {
                std::cout << "  A(" << i << ", " << j << ") = " << A(i, j) << std::endl;
            }
        }
        std::cout << "(new) AtransAinvAtrans =" << std::endl;
        for (unsigned int i = 0; i < AtransAinvAtrans.size1(); i++)
        {
            for (unsigned int j = 0; j < AtransAinvAtrans.size2(); j++)
            {
                std::cout << "  AtrAinvAtr(" << i << ", " << j << ") = " << AtransAinvAtrans(i, j) << std::endl;
            }
        }
        std::cout << "for node " << p_node->Coordinates() << std::endl;
    }

}
//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
template <std::size_t TDim>
void DerivativeRecovery<TDim>::ComputeDerivativeMonomialsVector(Node::Pointer& p_node, unsigned int& ord, DenseMatrix<double>& result)
{
    // If the gradient is computed at the node itself, then
    // the computation can be simplified

    KRATOS_ERROR_IF(ord > 3) << "Order 3 method has not been implemented for the ZG recovery!" << std::endl;

    for (unsigned int i = 0; i < result.size1(); i++)
    {
        for (unsigned int j = 0; j < result.size2(); j++)
        {
            result(i, j) = 0.0;
        }
    }

    // Constant terms
    // result(0, 0) = 0.0;  // d( 1 ) / dx
    // result(1, 0) = 0.0;  // d( 1 ) / dy
    // result(2, 0) = 0.0;  // d( 1 ) / dz

    // Terms of order 1
    // dP / dx
    result(0, 1) = 1.0; // d( x ) / dx
    // result(0, 2) = 0.0;  // d( y ) / dx
    // result(0, 3) = 0.0;  // d( z ) / dx

    // dP / dy
    // result(1, 1) = 0.0;  // d( x ) / dy
    result(1, 2) = 1.0; // d( y ) / dy
    // result(1, 3) = 0.0;  // d( z ) / dy

    // dP / dz
    // result(2, 1) = 0.0;  // d( x ) / dz
    // result(2, 2) = 0.0;  // d( y ) / dz
    result(2, 3) = 1.0; // d( z ) / dz
}
//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
template <std::size_t TDim>
void DerivativeRecovery<TDim>::ComputeDerivativeMonomialsVector(array_1d<double, 3>& pos, Node::Pointer& p_node, unsigned int& ord, DenseMatrix<double>& result)
{
    // Be aware that the derivatives here must be consistent with the definition provided
    // in the computation of the coefficients matrix

    KRATOS_ERROR_IF(ord > 3) << "Order 3 method has not been implemented for the ZG recovery!" << std::endl;

    // Constant terms and linear terms are constant (they do not depend on the relative coordinates)
    ComputeDerivativeMonomialsVector(p_node, ord, result);

    // const double h_inv = 1.0 / CalculateTheMaximumDistanceToNeighbours(p_node); // we use it as a scaling parameter to improve stability
    const double h_inv = 1.0;
    const array_1d <double, 3> origin = p_node->Coordinates();
    const array_1d<double, 3> rel_coordinates = h_inv * (pos - origin);

    // Terms of order 2
    if(ord > 1)
    {
        // dP / dx
        result(0, 4) = rel_coordinates[1];      // d ( x*y ) / dx 
        result(0, 5) = rel_coordinates[2];      // d ( x*z ) / dx
        result(0, 6) = 0.0;                     // d ( y*z ) / dx
        result(0, 7) = 2. * rel_coordinates[0]; // d ( x^2 ) / dx
        result(0, 8) = 0.0;                     // d ( y^2 ) / dx
        result(0, 9) = 0.0;                     // d ( z^2 ) / dx

        // dP / dy
        result(1, 4) = rel_coordinates[0];      // d ( x*y ) / dy 
        result(1, 5) = 0.0;                     // d ( x*z ) / dy
        result(1, 6) = rel_coordinates[2];      // d ( y*z ) / dy
        result(1, 7) = 0.0;                     // d ( x^2 ) / dy
        result(1, 8) = 2. * rel_coordinates[1]; // d ( y^2 ) / dy
        result(1, 9) = 0.0;                     // d ( z^2 ) / dy

        // dP / dz
        result(2, 4) = 0.0;                     // d ( x*y ) / dz 
        result(2, 5) = rel_coordinates[0];      // d ( x*z ) / dz
        result(2, 6) = rel_coordinates[1];      // d ( y*z ) / dz
        result(2, 7) = 0.0;                     // d ( x^2 ) / dz
        result(2, 8) = 0.0;                     // d ( y^2 ) / dz
        result(2, 9) = 2. * rel_coordinates[2]; // d ( z^2 ) / dz
    }

    // Terms of order 3
    if(ord > 2)
    {
        // dP / dx
        result(0, 10) = rel_coordinates[1] * rel_coordinates[2];       // d ( x*y*z ) / dx
        result(0, 11) = 2. * rel_coordinates[0] * rel_coordinates[1];  // d ( x^2*y ) / dx
        result(0, 12) = 2. * rel_coordinates[0] * rel_coordinates[2];  // d ( x^2*z ) / dx
        result(0, 13) = rel_coordinates[1] * rel_coordinates[1];       // d ( x*y^2 ) / dx
        result(0, 14) = rel_coordinates[2] * rel_coordinates[2];       // d ( x*z^2 ) / dx
        result(0, 15) = 0.0;                                           // d ( y^2*z ) / dx
        result(0, 16) = 0.0;                                           // d ( y*z^2 ) / dx
        result(0, 17) = 3. * rel_coordinates[0] * rel_coordinates[0];  // d ( x^3 ) / dx
        result(0, 18) = 0.0;                                           // d ( y^3 ) / dx
        result(0, 19) = 0.0;                                           // d ( z^3 ) / dx

        // dP / dy
        result(1, 10) = rel_coordinates[0] * rel_coordinates[2];       // d ( x*y*z ) / dy
        result(1, 11) = rel_coordinates[0] * rel_coordinates[0];       // d ( x^2*y ) / dy
        result(1, 12) = 0.0;                                           // d ( x^2*z ) / dy
        result(1, 13) = 2. * rel_coordinates[0] * rel_coordinates[1];  // d ( x*y^2 ) / dy
        result(1, 14) = 0.0;                                           // d ( x*z^2 ) / dy
        result(1, 15) = 2. * rel_coordinates[1] * rel_coordinates[2];  // d ( y^2*z ) / dy
        result(1, 16) = rel_coordinates[2] * rel_coordinates[2];       // d ( y*z^2 ) / dy
        result(1, 17) = 0.0;                                           // d ( x^3 ) / dy
        result(1, 18) = 3. * rel_coordinates[1] * rel_coordinates[1];  // d ( y^3 ) / dy
        result(1, 19) = 0.0;                                           // d ( z^3 ) / dy

        // dP / dz
        result(2, 10) = rel_coordinates[0] * rel_coordinates[1];       // d ( x*y*z ) / dz
        result(2, 11) = 0.0;                                           // d ( x^2*y ) / dz
        result(2, 12) = rel_coordinates[0] * rel_coordinates[0];       // d ( x^2*z ) / dz
        result(2, 13) = 0.0;                                           // d ( x*y^2 ) / dz
        result(2, 14) = 2. * rel_coordinates[0] * rel_coordinates[2];  // d ( x*z^2 ) / dz
        result(2, 15) = rel_coordinates[1] * rel_coordinates[1];       // d ( y^2*z ) / dz
        result(2, 16) = 2. * rel_coordinates[1] * rel_coordinates[2];  // d ( y*z^2 ) / dz
        result(2, 17) = 0.0;                                           // d ( x^3 ) / dz
        result(2, 18) = 0.0;                                           // d ( y^3 ) / dz
        result(2, 19) = 3. * rel_coordinates[2] * rel_coordinates[2];  // d ( z^3 ) / dz
    }                       

    // Normalization factor
    // for(unsigned d = 0; d < TDim; d++)
    // {
    //     for (unsigned int i = 0; i < result.size2(); i++)
    //     {
    //         // result(d, i) *= h_inv;
    //         result(d, i) *= h_inv;
    //     }
    // }
}
//**************************************************************************************************************************************************
//**************************************************************************************************************************************************

// Explicit instantiations
template class DerivativeRecovery<2>;
template class DerivativeRecovery<3>;
template void KRATOS_API(SWIMMING_DEM_APPLICATION) DerivativeRecovery<2>::RecoverSuperconvergentGradient< Variable<double> >(ModelPart&,  Variable<double>&, Variable<array_1d<double, 3> >&, unsigned int&);
template void KRATOS_API(SWIMMING_DEM_APPLICATION) DerivativeRecovery<3>::RecoverSuperconvergentGradient< Variable<double> >(ModelPart&,  Variable<double>&, Variable<array_1d<double, 3> >&, unsigned int&);
template void KRATOS_API(SWIMMING_DEM_APPLICATION) DerivativeRecovery<3>::RecoverSuperconvergentGradientAlt< Variable<double> >(ModelPart&,  Variable<double>&, Variable<array_1d<double, 3> >&, unsigned int&);
template void KRATOS_API(SWIMMING_DEM_APPLICATION) DerivativeRecovery<2>::CalculateGradient< Variable<double> >(ModelPart&,  Variable<double>&, Variable<array_1d<double, 3> >&);
template void KRATOS_API(SWIMMING_DEM_APPLICATION) DerivativeRecovery<3>::CalculateGradient< Variable<double> >(ModelPart&,  Variable<double>&, Variable<array_1d<double, 3> >&);
}  // namespace Kratos.

//
//   Project Name:        Kratos
//   Last Modified by:    $Author: gcasas $
//   Date:                $Date: 2016-03-08 08:56:42 $
//
//

#include "derivative_recovery.h"

namespace Kratos
{
template <std::size_t TDim>
void DerivativeRecovery<TDim>::RecoverGradientOfAScalar(const VariableData& origin_variable, const VariableData& destination_variable)
{
    for (int i = 0; i < (int)mModelPart.Nodes().size(); i++){
        NodeIteratorType i_particle = mModelPart.NodesBegin() + i;
        Node<3>::Pointer p_node = *(i_particle.base());
        array_1d<double, 3>& gradient = p_node->FastGetSolutionStepValue(TORQUE);
        gradient[0] = 0.0;
        gradient[1] = 0.0;
        gradient[2] = 99.0;
    }
}
//***************************************************************************************************************
//***************************************************************************************************************
template <std::size_t TDim>
void DerivativeRecovery<TDim>::CalculateVectorMaterialDerivative(ModelPart& r_model_part, Variable<array_1d<double, 3> >& vector_container, Variable<array_1d<double, 3> >& vector_rate_container, Variable<array_1d<double, 3> >& material_derivative_container)
{
    std::cout << "Constructing the material derivative by derivating nodal averages...\n";
    std::map <std::size_t, unsigned int> id_to_position;
    unsigned int entry = 0;
    double delta_time_inv = 1.0 / r_model_part.GetProcessInfo()[DELTA_TIME];

    for (NodeIteratorType inode = r_model_part.NodesBegin(); inode != r_model_part.NodesEnd(); inode++){
        noalias(inode->FastGetSolutionStepValue(material_derivative_container)) = ZeroVector(3);
        id_to_position[inode->Id()] = entry;
        entry++;
    }

    std::vector<array_1d <double, 3> > convective_contributions_to_the_derivative;
    convective_contributions_to_the_derivative.resize(entry);

    array_1d <double, 3> grad = ZeroVector(3);
    array_1d <double, TDim + 1 > elemental_values;
    array_1d <double, TDim + 1 > N; // shape functions vector
    boost::numeric::ublas::bounded_matrix<double, TDim + 1, TDim> DN_DX;

    for (unsigned int j = 0; j < TDim; ++j){ // for each component of the original vector value

        // for each element, constrno_balls_small_cellular_Post_Filesucting the gradient contribution (to its nodes) of the component v_j and storing it in material_derivative_container

        for (ModelPart::ElementIterator ielem = r_model_part.ElementsBegin(); ielem != r_model_part.ElementsEnd(); ielem++){
            // computing the shape function derivatives
            Geometry<Node<3> >& geom = ielem->GetGeometry();
            double Volume;
            GeometryUtils::CalculateGeometryData(geom, DN_DX, N, Volume);

            for (unsigned int i = 0; i < TDim + 1; ++i){
                elemental_values[i] = geom[i].FastGetSolutionStepValue(vector_container)[j];
            }

            array_1d <double, 3> grad_aux = prod(trans(DN_DX), elemental_values); // its dimension may be 2

            for (unsigned int i = 0; i < TDim; ++i){
                grad[i] = grad_aux[i];
            }

            double nodal_area = Volume / static_cast<double>(TDim + 1);
            grad *= nodal_area;

            for (unsigned int i = 0; i < TDim + 1; ++i){
                geom[i].FastGetSolutionStepValue(material_derivative_container) += grad; // we use material_derivative_container to store the gradient of one component at a time
            }
        }

        // normalizing the constributions to the gradient and getting the j-component of the material derivative

        for (NodeIteratorType inode = r_model_part.NodesBegin(); inode != r_model_part.NodesEnd(); inode++){
            array_1d <double, 3>& stored_gradient_of_component_j = inode->FastGetSolutionStepValue(material_derivative_container);
            stored_gradient_of_component_j /= inode->FastGetSolutionStepValue(NODAL_AREA);
            const array_1d <double, 3>& velocity = inode->FastGetSolutionStepValue(VELOCITY);
            convective_contributions_to_the_derivative[id_to_position[inode->Id()]][j] = DEM_INNER_PRODUCT_3(velocity, stored_gradient_of_component_j);
            stored_gradient_of_component_j = ZeroVector(3);
        }

    }

    for (NodeIteratorType inode = r_model_part.NodesBegin(); inode != r_model_part.NodesEnd(); inode++){
        const array_1d <double, 3>& stored_convective_contribution = convective_contributions_to_the_derivative[id_to_position[inode->Id()]];
        const array_1d <double, 3>& eulerian_rate_of_change = delta_time_inv * (inode->FastGetSolutionStepValue(VELOCITY) - inode->FastGetSolutionStepValue(VELOCITY, 1));
        array_1d <double, 3>& material_derivative = inode->FastGetSolutionStepValue(material_derivative_container);
        material_derivative = eulerian_rate_of_change + stored_convective_contribution;
    }

    std::cout << "Finished constructing the material derivative by derivating nodal averages...\n";
}
//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
template <std::size_t TDim>
void DerivativeRecovery<TDim>::CalculateGradient(ModelPart& r_model_part, Variable<double>& scalar_container, Variable<array_1d<double, 3> >& gradient_container)
{
    for (NodeIteratorType inode = r_model_part.NodesBegin(); inode != r_model_part.NodesEnd(); inode++){
        noalias(inode->FastGetSolutionStepValue(gradient_container)) = ZeroVector(3);
    }

    array_1d <double, 3> grad = ZeroVector(3); // its dimension is always 3
    array_1d <double, TDim + 1 > elemental_values;
    array_1d <double, TDim + 1 > N; // shape functions vector
    boost::numeric::ublas::bounded_matrix<double, TDim + 1, TDim> DN_DX;

    for (ModelPart::ElementIterator ielem = r_model_part.ElementsBegin(); ielem != r_model_part.ElementsEnd(); ielem++){
        // computing the shape function derivatives
        Geometry<Node<3> >& geom = ielem->GetGeometry();
        double Volume;

        GeometryUtils::CalculateGeometryData(geom, DN_DX, N, Volume);

        // getting the gradients;

        for (unsigned int i = 0; i < TDim + 1; ++i){
            elemental_values[i] = geom[i].FastGetSolutionStepValue(scalar_container);
        }


        array_1d <double, TDim> grad_aux = prod(trans(DN_DX), elemental_values); // its dimension may be 2

        for (unsigned int i = 0; i < TDim; ++i){
            grad[i] = grad_aux[i];
        }

        double nodal_area = Volume / static_cast<double>(TDim + 1);
        grad *= nodal_area;

        for (unsigned int i = 0; i < TDim + 1; ++i){
            geom[i].FastGetSolutionStepValue(gradient_container) += grad;
        }
    }

    for (NodeIteratorType inode = r_model_part.NodesBegin(); inode != r_model_part.NodesEnd(); inode++){
        inode->FastGetSolutionStepValue(gradient_container) /= inode->FastGetSolutionStepValue(NODAL_AREA);
    }
}
//**************************************************************************************************************************************************
//**************************************************************************************************************************************************

template <std::size_t TDim>
void DerivativeRecovery<TDim>::CalculateVectorLaplacian(ModelPart& r_model_part, Variable<array_1d<double, 3> >& vector_container, Variable<array_1d<double, 3> >& laplacian_container)
{
    std::cout << "Constructing the Laplacian by derivating nodal averages...\n";
    std::map <std::size_t, unsigned int> id_to_position;
    unsigned int entry = 0;

    for (NodeIteratorType inode = r_model_part.NodesBegin(); inode != r_model_part.NodesEnd(); inode++){
        noalias(inode->FastGetSolutionStepValue(laplacian_container)) = ZeroVector(3);
        id_to_position[inode->Id()] = entry;
        entry++;
    }

    std::vector<array_1d <double, 3> > laplacians;
    laplacians.resize(entry);

    array_1d <double, 3> grad = ZeroVector(3);
    array_1d <double, TDim + 1 > elemental_values;
    array_1d <double, TDim + 1 > N; // shape functions vector
    boost::numeric::ublas::bounded_matrix<double, TDim + 1, TDim> DN_DX;
    boost::numeric::ublas::bounded_matrix<double, TDim + 1, TDim> elemental_vectors; // They carry the nodal gradients of the corresponding component v_j


    for (unsigned int j = 0; j < TDim; ++j){ // for each component of the original vector value

        // for each element, constructing the gradient contribution (to its nodes) of the component v_j and storing it in laplacian_container

        for (ModelPart::ElementIterator ielem = r_model_part.ElementsBegin(); ielem != r_model_part.ElementsEnd(); ielem++){
            // computing the shape function derivatives
            Geometry<Node<3> >& geom = ielem->GetGeometry();
            double Volume;
            GeometryUtils::CalculateGeometryData(geom, DN_DX, N, Volume);

            for (unsigned int i = 0; i < TDim + 1; ++i){
                elemental_values[i] = geom[i].FastGetSolutionStepValue(vector_container)[j];
            }

            array_1d <double, TDim> grad_aux = prod(trans(DN_DX), elemental_values); // its dimension may be 2

            for (unsigned int i = 0; i < TDim; ++i){
                grad[i] = grad_aux[i];
            }

            double nodal_area = Volume / static_cast<double>(TDim + 1);
            grad *= nodal_area;

            for (unsigned int i = 0; i < TDim + 1; ++i){
                geom[i].FastGetSolutionStepValue(laplacian_container) += grad; // we use laplacian_container to store the gradient of one component at a time
            }
        }

        // normalizing the constributions to the gradient

        for (NodeIteratorType inode = r_model_part.NodesBegin(); inode != r_model_part.NodesEnd(); inode++){
            inode->FastGetSolutionStepValue(laplacian_container) /= inode->FastGetSolutionStepValue(NODAL_AREA);
        }

        // for each element, constructing the divergence contribution (to its nodes) of the gradient of component v_j

        for (ModelPart::ElementIterator ielem = r_model_part.ElementsBegin(); ielem != r_model_part.ElementsEnd(); ielem++){
            Geometry<Node<3> >& geom = ielem->GetGeometry();
            double Volume;
            GeometryUtils::CalculateGeometryData(geom, DN_DX, N, Volume);

            for (unsigned int i = 0; i < TDim + 1; ++i){
                for (unsigned int k = 0; k < TDim; ++k){
                    elemental_vectors(i, k) = geom[i].FastGetSolutionStepValue(laplacian_container)[k]; // it is actually the gradient of component v_j
                }
            }

            boost::numeric::ublas::bounded_matrix<double, TDim, TDim> grad_aux = prod(trans(DN_DX), elemental_vectors); // its dimension may be 2
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

        for (NodeIteratorType inode = r_model_part.NodesBegin(); inode != r_model_part.NodesEnd(); inode++){
            array_1d <double, 3>& current_gradient_of_vi = inode->FastGetSolutionStepValue(laplacian_container);
            current_gradient_of_vi = ZeroVector(3);
        }

    } // for each component of the vector value

    for (NodeIteratorType inode = r_model_part.NodesBegin(); inode != r_model_part.NodesEnd(); inode++){
        array_1d <double, 3>& stored_laplacian = laplacians[id_to_position[inode->Id()]];
        array_1d <double, 3>& laplacian = inode->FastGetSolutionStepValue(laplacian_container);
        laplacian = stored_laplacian / inode->FastGetSolutionStepValue(NODAL_AREA);
    }

    std::cout << "Finished constructing the Laplacian by derivating nodal averages...\n";
}

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************

template <std::size_t TDim>
void DerivativeRecovery<TDim>::RecoverSuperconvergentLaplacian(ModelPart& r_model_part, Variable<array_1d<double, 3> >& vector_container, Variable<array_1d<double, 3> >& laplacian_container)
{
    mCalculatingTheLaplacian = true;

    if (mFirstLaplacianRecovery){
        std::cout << "Constructing first-step neighbour clouds for laplacian recovery...\n";
        SetNeighboursAndWeights(r_model_part);
        mFirstLaplacianRecovery = false;
        std::cout << "Finished constructing neighbour clouds for laplacian recovery.\n";
    }

    if (mSomeCloudsDontWork){ // a default value is necessary in the cases where recovery is not possible
        CalculateVectorLaplacian(r_model_part, vector_container, laplacian_container);
    }

    // Solving least squares problem (Zhang, 2006)
    unsigned int n_relevant_terms = 6;

    std::vector<array_1d <double, 3> > polynomial_coefficients; // vector to store, for each node, the corresponding values of the polynomial coefficients relevant for the calculation of the laplacian
    polynomial_coefficients.resize(n_relevant_terms);

    for (NodeIteratorType inode = r_model_part.NodesBegin(); inode != r_model_part.NodesEnd(); inode++){
        WeakPointerVector<Node<3> >& neigh_nodes = inode->GetValue(NEIGHBOUR_NODES);
        unsigned int n_neigh = neigh_nodes.size();

        if (!n_neigh){ // then we keep the defualt value
            continue;
        }

        for (unsigned int i = 0; i < n_relevant_terms; i++){ // resetting polynomial_coefficients to 0
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
        recovered_laplacian[0] = 2 * (polynomial_coefficients[3][0] + polynomial_coefficients[4][0] + polynomial_coefficients[5][0]);
        recovered_laplacian[1] = 2 * (polynomial_coefficients[3][1] + polynomial_coefficients[4][1] + polynomial_coefficients[5][1]);
        recovered_laplacian[2] = 2 * (polynomial_coefficients[3][2] + polynomial_coefficients[4][2] + polynomial_coefficients[5][2]);
//        recovered_laplacian[0] = polynomial_coefficients[0][0] + polynomial_coefficients[1][0] + 2 * polynomial_coefficients[3][0];
//        recovered_laplacian[1] = polynomial_coefficients[0][0] + polynomial_coefficients[2][0] + 2 * polynomial_coefficients[4][0];
//        recovered_laplacian[2] = polynomial_coefficients[1][0] + polynomial_coefficients[2][0] + 2 * polynomial_coefficients[5][0];
    }

    mCalculatingTheLaplacian = false;
}


// Explicit instantiations
template class DerivativeRecovery<2>;
template class DerivativeRecovery<3>;

}  // namespace Kratos.

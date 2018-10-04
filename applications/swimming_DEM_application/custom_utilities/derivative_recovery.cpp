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
    for (int i = 0; i < (int)mModelPart.Nodes().size(); ++i){
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
void DerivativeRecovery<TDim>::AddTimeDerivative(ModelPart& r_model_part, Variable<array_1d<double, 3> >& material_derivative_container)
{
    const double delta_time_inv = 1.0 / r_model_part.GetProcessInfo()[DELTA_TIME];

    for (NodeIteratorType inode = r_model_part.NodesBegin(); inode != r_model_part.NodesEnd(); ++inode){
        array_1d <double, 3>& material_derivative = inode->FastGetSolutionStepValue(material_derivative_container);
        const array_1d <double, 3> eulerian_rate_of_change = delta_time_inv * (inode->FastGetSolutionStepValue(VELOCITY) - inode->FastGetSolutionStepValue(VELOCITY, 1));
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
void DerivativeRecovery<TDim>::CalculateVectorMaterialDerivative(ModelPart& r_model_part,
                                                                 Variable<array_1d<double, 3> >& vector_container,
                                                                 Variable<array_1d<double, 3> >& vector_rate_container,
                                                                 Variable<array_1d<double, 3> >& material_derivative_container)
{
    KRATOS_INFO("DEM-FLUID") << "Constructing the material derivative by derivating nodal averages..." << std::endl;
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
    array_1d <double, TDim + 1 > elemental_values;
    array_1d <double, TDim + 1 > N; // shape functions vector
    BoundedMatrix<double, TDim + 1, TDim> DN_DX;

    for (unsigned int j = 0; j < TDim; ++j){ // for each component of the original vector value

        // for each element, constructing the gradient contribution (to its nodes) of the component v_j and storing it in material_derivative_container

        for (ModelPart::ElementIterator ielem = r_model_part.ElementsBegin(); ielem != r_model_part.ElementsEnd(); ++ielem){
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

    KRATOS_INFO("DEM-FLUID") << "Finished constructing the material derivative by derivating nodal averages..." << std::endl;
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
        Node<3>::Pointer p_node = *(i_particle.base());
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

    if (current_component != 0 && current_component != 1 && current_component != 2){
        KRATOS_THROW_ERROR(std::invalid_argument,"The value of CURRENT_COMPONENT passed to the ComputeComponentGradientSimplex element is not 0, 1 or 2, but ", current_component);
    }

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
void DerivativeRecovery<TDim>::RecoverSuperconvergentMatDeriv(ModelPart& r_model_part, Variable<array_1d<double, 3> >& vector_container, Variable<array_1d<double, 3> >& vector_rate_container, Variable<array_1d<double, 3> >& material_derivative_container)
{
    mCalculatingTheGradient = true;

    if (mFirstGradientRecovery){
        KRATOS_INFO("DEM-FLUID") << "Constructing first-step neighbour clouds for material derivative..." << std::endl;
        SetNeighboursAndWeights(r_model_part);
        mFirstGradientRecovery = false;
        KRATOS_INFO("DEM-FLUID") << "Finished constructing neighbour clouds for material derivative." << std::endl;
    }

    if (mSomeCloudsDontWork){ // a default value is necessary in the cases where recovery is not possible
        CalculateVectorMaterialDerivative(r_model_part, vector_container, vector_rate_container, material_derivative_container);
    }

    // Solving least squares problem (Zhang, 2006)
    unsigned int n_relevant_terms = 3;

    std::vector<array_1d <double, 3> > polynomial_coefficients; // vector to store, for each node, the corresponding values of the polynomial coefficients relevant for the calculation of the material derivative
    polynomial_coefficients.resize(n_relevant_terms);

    for (NodeIteratorType inode = r_model_part.NodesBegin(); inode != r_model_part.NodesEnd(); ++inode){
        WeakPointerVector<Node<3> >& neigh_nodes = inode->GetValue(NEIGHBOUR_NODES);
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
        Node<3>::Pointer p_node = *(i_particle.base());
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

    if (current_component != 0 && current_component != 1 && current_component != 2){
        KRATOS_THROW_ERROR(std::invalid_argument,"The value of CURRENT_COMPONENT passed to the ComputeComponentGradientSimplex element is not 0, 1 or 2, but ", current_component);
    }

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
    array_1d <double, TDim + 1 > elemental_values;
    array_1d <double, TDim + 1 > N; // shape functions vector
    BoundedMatrix<double, TDim + 1, TDim> DN_DX;

    for (ModelPart::ElementIterator ielem = r_model_part.ElementsBegin(); ielem != r_model_part.ElementsEnd(); ++ielem){

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

        Geometry<Node<3> >& geom = ielem->GetGeometry();
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
void DerivativeRecovery<TDim>::RecoverSuperconvergentGradient(ModelPart& r_model_part, TScalarVariable& scalar_container, Variable<array_1d<double, 3> >& gradient_container)
{
    mCalculatingTheGradient = true;

    if (mFirstGradientRecovery){
        KRATOS_INFO("DEM-FLUID") << "Constructing first-step neighbour clouds for gradient recovery..." << std::endl;
        SetNeighboursAndWeights(r_model_part);
        mFirstGradientRecovery = false;
        KRATOS_INFO("DEM-FLUID") << "Finished constructing neighbour clouds for gradient recovery." << std::endl;
    }

    if (mSomeCloudsDontWork){ // a default value is necessary in the cases where recovery is not possible
        CalculateGradient(r_model_part, scalar_container, gradient_container);
    }

    // Solving least squares problem (Zhang, 2006)

    for (NodeIteratorType inode = r_model_part.NodesBegin(); inode != r_model_part.NodesEnd(); ++inode){
        WeakPointerVector<Node<3> >& neigh_nodes = inode->GetValue(NEIGHBOUR_NODES);
        unsigned int n_neigh = neigh_nodes.size();

        if (!n_neigh){ // we keep the defualt value
            continue;
        }

        array_1d <double, 3>& recovered_gradient = inode->FastGetSolutionStepValue(gradient_container);
        recovered_gradient = ZeroVector(3);
        const Vector& nodal_weights = inode->FastGetSolutionStepValue(NODAL_WEIGHTS);

        for (unsigned int i_neigh = 0; i_neigh < n_neigh; ++i_neigh){
            const double& neigh_nodal_value = neigh_nodes[i_neigh].FastGetSolutionStepValue(scalar_container);

            for (unsigned int d = 0; d < TDim; ++d){
                recovered_gradient[d] += nodal_weights[3 * i_neigh + d] * neigh_nodal_value;
            }
        }
    }

    mCalculatingTheGradient = false;
}
//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
//template <std::size_t TDim>
//void DerivativeRecovery<TDim>::RecoverSuperconvergentGradient(ModelPart& r_model_part, const VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > >& scalar_container, Variable<array_1d<double, 3> >& gradient_container)
//{
//    mCalculatingTheGradient = true;

//    if (mFirstGradientRecovery){
//        std::cout << "Constructing first-step neighbour clouds for gradient recovery...\n";
//        SetNeighboursAndWeights(r_model_part);
//        mFirstGradientRecovery = false;
//        std::cout << "Finished constructing neighbour clouds for gradient recovery.\n";
//    }

//    if (mSomeCloudsDontWork){ // a default value is necessary in the cases where recovery is not possible
//        CalculateGradient(r_model_part, scalar_container, gradient_container);
//    }

//    // Solving least squares problem (Zhang, 2006)

//    for (NodeIteratorType inode = r_model_part.NodesBegin(); inode != r_model_part.NodesEnd(); ++inode){
//        WeakPointerVector<Node<3> >& neigh_nodes = inode->GetValue(NEIGHBOUR_NODES);
//        unsigned int n_neigh = neigh_nodes.size();

//        if (!n_neigh){ // we keep the defualt value
//            continue;
//        }

//        array_1d <double, 3>& recovered_gradient = inode->FastGetSolutionStepValue(gradient_container);
//        recovered_gradient = ZeroVector(3);
//        const Vector& nodal_weights = inode->FastGetSolutionStepValue(NODAL_WEIGHTS);

//        for (unsigned int i_neigh = 0; i_neigh < n_neigh; ++i_neigh){
//            const double& neigh_nodal_value = neigh_nodes[i_neigh].FastGetSolutionStepValue(scalar_container);

//            for (unsigned int d = 0; d < TDim; ++d){
//                recovered_gradient[d] += nodal_weights[3 * i_neigh + d] * neigh_nodal_value;
//            }
//        }
//    }

//    mCalculatingTheGradient = false;
//}
//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
template <std::size_t TDim>
void DerivativeRecovery<TDim>::RecoverSuperconvergentLaplacian(ModelPart& r_model_part, Variable<array_1d<double, 3> >& vector_container, Variable<array_1d<double, 3> >& laplacian_container)
{
    mCalculatingTheLaplacian = true;

    if (mFirstLaplacianRecovery){
        KRATOS_INFO("DEM-FLUID") << "Constructing first-step neighbour clouds for laplacian recovery..." << std::endl;
        SetNeighboursAndWeights(r_model_part);
        mFirstLaplacianRecovery = false;
        KRATOS_INFO("DEM-FLUID") << "Finished constructing neighbour clouds for laplacian recovery." << std::endl;
    }

    if (mSomeCloudsDontWork){ // a default value is necessary in the cases where recovery is not possible
        CalculateVectorLaplacian(r_model_part, vector_container, laplacian_container);
    }

    // Solving a least squares problem for every node (Zhang, 2006)
    unsigned int n_relevant_terms = 6;

    std::vector<array_1d <double, 3> > polynomial_coefficients; // vector to store, for each node, the corresponding values of the polynomial coefficients relevant for the calculation of the laplacian
    polynomial_coefficients.resize(n_relevant_terms);

    for (NodeIteratorType inode = r_model_part.NodesBegin(); inode != r_model_part.NodesEnd(); ++inode){
        WeakPointerVector<Node<3> >& neigh_nodes = inode->GetValue(NEIGHBOUR_NODES);
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
}
//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
template <std::size_t TDim>
void DerivativeRecovery<TDim>::RecoverSuperconvergentVelocityLaplacianFromGradient(ModelPart& r_model_part, Variable<array_1d<double, 3> >& vector_container, Variable<array_1d<double, 3> >& laplacian_container)
{
    mCalculatingTheGradient = true;

    if (mFirstLaplacianRecovery){
        KRATOS_INFO("DEM-FLUID") << "Finished constructing neighbour clouds for laplacian recovery." << std::endl;
        SetNeighboursAndWeights(r_model_part);
        mFirstLaplacianRecovery = false;
        KRATOS_INFO("DEM-FLUID") << "Finished constructing neighbour clouds for laplacian recovery." << std::endl;
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
        WeakPointerVector<Node<3> >& neigh_nodes = inode->GetValue(NEIGHBOUR_NODES);
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
}
//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
template <std::size_t TDim>
void DerivativeRecovery<TDim>::RecoverSuperconvergentMatDerivAndLaplacian(ModelPart& r_model_part, Variable<array_1d<double, 3> >& vector_container, Variable<array_1d<double, 3> >& vector_rate_container, Variable<array_1d<double, 3> >& material_derivative_container, Variable<array_1d<double, 3> >& laplacian_container)
{
    mCalculatingGradientAndLaplacian = true;

    if (mFirstLaplacianRecovery){
        KRATOS_INFO("DEM-FLUID") << "Constructing first-step neighbour clouds for material derivative and laplacian recovery..." << std::endl;
        SetNeighboursAndWeights(r_model_part);
        mFirstLaplacianRecovery = false;
        KRATOS_INFO("DEM-FLUID") << "Finished constructing neighbour clouds for material derivative and laplacian recovery." << std::endl;
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
        WeakPointerVector<Node<3> >& neigh_nodes = inode->GetValue(NEIGHBOUR_NODES);
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
}
//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
template <std::size_t TDim>
void DerivativeRecovery<TDim>::CalculateVectorLaplacian(ModelPart& r_model_part, Variable<array_1d<double, 3> >& vector_container, Variable<array_1d<double, 3> >& laplacian_container)
{
    KRATOS_INFO("DEM-FLUID") << "Constructing the Laplacian by derivating nodal averages..." << std::endl;
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
            Geometry<Node<3> >& geom = ielem->GetGeometry();
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

    KRATOS_INFO("DEM-FLUID") << "Finished constructing the Laplacian by derivating nodal averages..." << std::endl;
}
//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
template <std::size_t TDim>
void DerivativeRecovery<TDim>::CalculateVelocityLaplacianRate(ModelPart& r_model_part)
{
    double delta_t_inv = 1.0 / r_model_part.GetProcessInfo()[DELTA_TIME];
    DenseVector<unsigned int> nodes_partition;
    OpenMPUtils::CreatePartition(OpenMPUtils::GetNumThreads(), r_model_part.Nodes().size(), nodes_partition);

    #pragma omp parallel for
    for (int k = 0; k < OpenMPUtils::GetNumThreads(); ++k){
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
void DerivativeRecovery<TDim>::SetNeighboursAndWeights(ModelPart& r_model_part)
{
    // Finding elements concurrent to each node. The nodes of these elements will form the initial cloud of points
    FindNodalNeighboursProcess neighbour_finder = FindNodalNeighboursProcess(r_model_part);
    neighbour_finder.Execute();
    const unsigned int n_max_iterations = 100;

    unsigned int i = 0;
    for (NodeIteratorType inode = r_model_part.NodesBegin(); inode != r_model_part.NodesEnd(); ++inode){
        bool the_cloud_of_neighbours_is_successful = SetInitialNeighboursAndWeights(r_model_part, *(inode.base()));
        WeakPointerVector<Node<3> >& neigh_nodes = inode->GetValue(NEIGHBOUR_NODES);

        unsigned int iteration = 0;
        while (!the_cloud_of_neighbours_is_successful && iteration < n_max_iterations){
            the_cloud_of_neighbours_is_successful = SetNeighboursAndWeights(r_model_part, *(inode.base()));
            ++iteration;
        }

        if (iteration >= n_max_iterations){ // giving up on this method, settling for the default method
            mSomeCloudsDontWork = true;
            neigh_nodes.clear();
            inode->FastGetSolutionStepValue(NODAL_WEIGHTS).clear();
            KRATOS_WARNING("DEM-FLUID") << "Warning!, for the node with id " << inode->Id() << " it has not been possible to form an adequate cloud of neighbours" << std::endl;
            KRATOS_WARNING("DEM-FLUID") << "for the gradient recovery. A lower accuracy method has been employed for this node." << std::endl;
        }

        ++i;
    }
}
//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
template <std::size_t TDim>
void DerivativeRecovery<TDim>::SetNeighboursAndWeightsForTheLaplacian(ModelPart& r_model_part)
{
    // Finding elements concurrent to each node. The nodes of these elements will form the initial cloud of points
    FindNodalNeighboursProcess neighbour_finder = FindNodalNeighboursProcess(r_model_part);
    neighbour_finder.Execute();
    const unsigned int n_max_iterations = 100;

    unsigned int i = 0;
    for (NodeIteratorType inode = r_model_part.NodesBegin(); inode != r_model_part.NodesEnd(); ++inode){
        bool the_cloud_of_neighbours_is_successful = SetInitialNeighboursAndWeights(r_model_part, *(inode.base()));
        WeakPointerVector<Node<3> >& neigh_nodes = inode->GetValue(NEIGHBOUR_NODES);

        unsigned int iteration = 0;
        while (!the_cloud_of_neighbours_is_successful && iteration < n_max_iterations){
            the_cloud_of_neighbours_is_successful = SetNeighboursAndWeights(r_model_part, *(inode.base()));
            ++iteration;
        }

        if (iteration >= n_max_iterations){ // giving up on this method, settling for the default method
            mSomeCloudsDontWork = true;
            neigh_nodes.clear();
            inode->FastGetSolutionStepValue(NODAL_WEIGHTS).clear();
            KRATOS_WARNING("DEM-FLUID") << "Warning!, for the node with id " << inode->Id() << " it has not been possible to form an adequate cloud of neighbours" << std::endl;
            KRATOS_WARNING("DEM-FLUID") << "for the gradient recovery. A lower accuracy method has been employed for this node." << std::endl;
        }

        ++i;
    }
}
//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
template <std::size_t TDim>
void DerivativeRecovery<TDim>::OrderByDistance(Node<3>::Pointer &p_node, WeakPointerVector<Node<3> >& neigh_nodes)
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
    WeakPointerVector<Node<3> > ordered_neighbours;

    for (unsigned int i = 0; i < n_nodes; ++i){
        Node<3>::WeakPointer& p_neigh = neigh_nodes(ordering[i].first);
        ordered_neighbours.push_back(p_neigh);
    }

    ordered_neighbours.swap(neigh_nodes);
}
//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
template <std::size_t TDim>
bool DerivativeRecovery<TDim>::SetInitialNeighboursAndWeights(ModelPart& r_model_part, Node<3>::Pointer &p_node)
{
    WeakPointerVector<Element>& neigh_elems = p_node->GetValue(NEIGHBOUR_ELEMENTS);
    WeakPointerVector<Node<3> >& neigh_nodes = p_node->GetValue(NEIGHBOUR_NODES);
    std::map<std::size_t, std::size_t> ids; // map to keep track of all different ids corresponding to already added neighbours to avoid repetition
    ids[p_node->Id()] = p_node->Id();

    unsigned int i = 0;

    for (unsigned int i_el = 0; i_el < neigh_elems.size(); ++i_el){
        Geometry<Node<3> >& geom = neigh_elems[i_el].GetGeometry();

        unsigned int jj = 0; // index of the node in geom corresponding to neighbour neigh_elems[i_el]
        if (geom[jj].Id() == p_node->Id()){ // skipping itself
            ++jj;
        }

        for (unsigned int j = 0; j < TDim; ++j){
            Node<3>::Pointer p_neigh = geom(jj);

            if (ids.find(p_neigh->Id()) == ids.end()){
                neigh_nodes.push_back(p_neigh);
                ids[p_neigh->Id()] = p_neigh->Id();
            }
        }

        i += TDim;
    }

    OrderByDistance(p_node, neigh_nodes);

    if (neigh_nodes.size() < 10){ // Not worthwhile checking, since there are 10 independent coefficients to be determined
        return false;
    }

    else {
        return SetWeightsAndRunLeastSquaresTest(r_model_part, p_node);
    }
}
//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
template <std::size_t TDim>
bool DerivativeRecovery<TDim>::SetNeighboursAndWeights(ModelPart& r_model_part, Node<3>::Pointer& p_node)
{
    WeakPointerVector<Node<3> >& neigh_nodes = p_node->GetValue(NEIGHBOUR_NODES);
    const unsigned int node_increase_per_neighbour = 1;
    const unsigned int node_increase_overall = 1;
    std::map<std::size_t, std::size_t> ids;
    ids[p_node->Id()] = p_node->Id();

    for (unsigned int i = 0; i < (unsigned int)neigh_nodes.size(); ++i){
        Node<3>::Pointer p_neigh = neigh_nodes(i).lock();
        ids[p_neigh->Id()] = p_neigh->Id();
    }

    const unsigned int n_neigh = neigh_nodes.size();

    for (unsigned int i = 0; i < n_neigh; ++i){
        Node<3>::Pointer p_neigh = neigh_nodes(i).lock();
        WeakPointerVector<Node<3> >& neigh_neigh_nodes = p_neigh->GetValue(NEIGHBOUR_NODES);
        unsigned int n_new_nodes = 0;
        for (unsigned int j = 0; j < (unsigned int)neigh_neigh_nodes.size(); ++j){
            Node<3>::Pointer p_neigh_neigh = neigh_neigh_nodes(j).lock();
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

    if (neigh_nodes.size() < 10){ // it is not worthwhile checking, since there are 10 independent coefficients to be determined
        return false;
    }

    else {
        return SetWeightsAndRunLeastSquaresTest(r_model_part, p_node);
    }
}
//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
template <std::size_t TDim>
double DerivativeRecovery<TDim>::SecondDegreeTestPolynomial(const array_1d <double, 3>& coordinates)
{
    const double x = coordinates[0];
    const double y = coordinates[1];
    const double z = coordinates[2];
    return 1.0 + x + y + z + x * y + x * z + y * z + x * x + y * y + z * z;
}
//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
template <std::size_t TDim>
double DerivativeRecovery<TDim>::SecondDegreeGenericPolynomial(DenseMatrix<double> C, const array_1d <double, 3>& coordinates)
{
    const double x = coordinates[0];
    const double y = coordinates[1];
    const double z = coordinates[2];
    return C(0,0) + C(1,0) * x + C(2,0) * y + C(3,0) * z + C(4,0) * x * y + C(5,0) * x * z + C(6,0) * y * z + C(7,0) * x * x + C(8,0) * y * y + C(9,0) * z * z;
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
bool DerivativeRecovery<TDim>::SetWeightsAndRunLeastSquaresTest(ModelPart& r_model_part, Node<3>::Pointer& p_node)
{
    unsigned int n_poly_terms = Factorial(TDim + 2) / (2 * Factorial(TDim)); // 2 is the polynomial order

    if (TDim == 2){
        KRATOS_THROW_ERROR(std::runtime_error, "Gradient recovery not implemented yet in 2D!)","");
    }

    WeakPointerVector<Node<3> >& neigh_nodes = p_node->GetValue(NEIGHBOUR_NODES);
    unsigned int n_nodal_neighs = (unsigned int)neigh_nodes.size();
    const double h_inv = 1.0 / CalculateTheMaximumDistanceToNeighbours(p_node); // we use it as a scaling parameter to improve stability
    const array_1d <double, 3> origin = p_node->Coordinates();
    DenseMatrix<double> TestNodalValues(n_nodal_neighs, 1);
    DenseMatrix<double> A(n_nodal_neighs, n_poly_terms);

    for (unsigned int i = 0; i < n_nodal_neighs; ++i){
        A(i, 0) = 1.0;

        if (TDim == 3){
            Node<3>& neigh = neigh_nodes[i];
            const array_1d <double, 3> rel_coordinates = (neigh.Coordinates() - origin) * h_inv;
            TestNodalValues(i, 0) = SecondDegreeTestPolynomial(rel_coordinates);

            for (unsigned int d = 1; d < 10; ++d){
                if (d < 4){
                    A(i, d) = rel_coordinates[d - 1];
                }
                else if (d == 4){
                    A(i, d) = rel_coordinates[0] * rel_coordinates[1];
                }
                else if (d == 5){
                    A(i, d) = rel_coordinates[0] * rel_coordinates[2];
                }
                else if (d == 6){
                    A(i, d) = rel_coordinates[1] * rel_coordinates[2];
                }
                else {
                    A(i, d) = rel_coordinates[d - 7] * rel_coordinates[d - 7];
                }
            }
        }

        else {
            KRATOS_THROW_ERROR(std::runtime_error,"Gradient recovery not implemented yet in 2D!)","");
        }
    }

    DenseMatrix<double>AtransA(n_poly_terms, n_poly_terms);
    noalias(AtransA) = prod(trans(A), A);

    if (std::abs(mMyCustomFunctions.template determinant< DenseMatrix<double> >(AtransA)) < 0.01){
        return false;
    }

    else {
        unsigned int n_relevant_terms = 0;

        if (mCalculatingTheGradient){
            n_relevant_terms = TDim;
        }

        else if (mCalculatingTheLaplacian){
            n_relevant_terms = n_poly_terms - (TDim + 1);
        }

        else {
            n_relevant_terms = n_poly_terms - 1;
        }

        std::vector<unsigned int> relevant_terms;
        relevant_terms.resize(n_relevant_terms);
        double normalization =  1.0;

        if (mCalculatingTheGradient){
            normalization = h_inv;
            relevant_terms[0] = 1;
            relevant_terms[1] = 2;
            relevant_terms[2] = 3;
        }

        else if (mCalculatingTheLaplacian){
            normalization = h_inv * h_inv;
            relevant_terms[0] = 4;
            relevant_terms[1] = 5;
            relevant_terms[2] = 6;
            relevant_terms[3] = 7;
            relevant_terms[4] = 8;
            relevant_terms[5] = 9;
        }

        else {
            for (unsigned int i = 1; i < n_poly_terms; ++i){
                relevant_terms[i - 1] = i;
            }
        }

        Vector& nodal_weights = p_node->FastGetSolutionStepValue(NODAL_WEIGHTS);
        nodal_weights.resize(n_relevant_terms * n_nodal_neighs);
        const DenseMatrix<double> AtransAinv = mMyCustomFunctions.Inverse(AtransA);
//        for (unsigned i = 0; i < n_poly_terms; i++){
//            for (unsigned j = 0; j < n_poly_terms; j++){
//                if (abs(AtransAinv(i,j))>1e6){
//                    return false;
//                }
//            }
//        }
        DenseMatrix<double>AtransAinvAtrans(n_poly_terms, n_nodal_neighs);
        noalias(AtransAinvAtrans) = prod(AtransAinv, trans(A));

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

        DenseMatrix<double> C(n_nodal_neighs, 1);
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
}
//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
template <std::size_t TDim>
unsigned int DerivativeRecovery<TDim>::GetNumberOfUniqueNeighbours(const int my_id, const WeakPointerVector<Element>& my_neighbour_elements)
{
    std::vector<int> ids;
    ids.push_back(my_id);

    for (unsigned int i_el = 0; i_el < my_neighbour_elements.size(); ++i_el){
        const Geometry<Node<3> >& geom = my_neighbour_elements[i_el].GetGeometry();
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
double DerivativeRecovery<TDim>::CalculateTheMaximumDistanceToNeighbours(Node<3>::Pointer& p_node)
{
    double max_distance_yet = 0.0;
    const array_1d <double, 3>& coors = p_node->Coordinates();
    WeakPointerVector<Node<3> >& neigh_nodes = p_node->GetValue(NEIGHBOUR_NODES);

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
        Geometry<Node<3> >& geom = ielem->GetGeometry();
        unsigned int n_nodes = static_cast<unsigned int>(TDim + 1);

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
        Geometry<Node<3> >& geom = ielem->GetGeometry();

        if (first_node){ // assign the distance (squared) between any two nodes to min_distance_yet
            array_1d <double, 3> delta = geom[0] - geom[1];
            double distance_2 = DEM_INNER_PRODUCT_3(delta, delta);
            min_distance_yet = distance_2;
        }

        unsigned int n_nodes = static_cast<unsigned int>(TDim + 1);

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


// Explicit instantiations
template class DerivativeRecovery<2>;
template class DerivativeRecovery<3>;

template void KRATOS_API(SWIMMING_DEM_APPLICATION) DerivativeRecovery<2>::RecoverSuperconvergentGradient< Variable<double> >(ModelPart&,  Variable<double>&, Variable<array_1d<double, 3> >&);
template void KRATOS_API(SWIMMING_DEM_APPLICATION) DerivativeRecovery<3>::RecoverSuperconvergentGradient< Variable<double> >(ModelPart&,  Variable<double>&, Variable<array_1d<double, 3> >&);
template void KRATOS_API(SWIMMING_DEM_APPLICATION) DerivativeRecovery<2>::RecoverSuperconvergentGradient< VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > >& >(ModelPart&,  VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > >&, Variable<array_1d<double, 3> >&);
template void KRATOS_API(SWIMMING_DEM_APPLICATION) DerivativeRecovery<3>::RecoverSuperconvergentGradient< VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > >& >(ModelPart&,  VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > >&, Variable<array_1d<double, 3> >&);


template void KRATOS_API(SWIMMING_DEM_APPLICATION) DerivativeRecovery<2>::CalculateGradient< Variable<double> >(ModelPart&,  Variable<double>&, Variable<array_1d<double, 3> >&);
template void KRATOS_API(SWIMMING_DEM_APPLICATION) DerivativeRecovery<3>::CalculateGradient< Variable<double> >(ModelPart&,  Variable<double>&, Variable<array_1d<double, 3> >&);
template void KRATOS_API(SWIMMING_DEM_APPLICATION) DerivativeRecovery<2>::CalculateGradient< VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > >& >(ModelPart&,  VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > >&, Variable<array_1d<double, 3> >&);
template void KRATOS_API(SWIMMING_DEM_APPLICATION) DerivativeRecovery<3>::CalculateGradient< VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > >& >(ModelPart&,  VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > >&, Variable<array_1d<double, 3> >&);

}  // namespace Kratos.

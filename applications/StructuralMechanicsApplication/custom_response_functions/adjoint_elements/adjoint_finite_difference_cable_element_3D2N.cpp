// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Martin Fusseder, https://github.com/MFusseder
//


#include "adjoint_finite_difference_cable_element_3D2N.h"
#include "structural_mechanics_application_variables.h"
#include "custom_response_functions/response_utilities/stress_response_definitions.h"
#include "custom_elements/truss_elements/cable_element_3D2N.hpp"
#include "custom_response_functions/response_utilities/finite_difference_utility.h"


namespace Kratos
{

template <class TPrimalElement>
void AdjointFiniteDifferenceCableElement<TPrimalElement>::CalculateSensitivityMatrix(const Variable<double>& rDesignVariable, Matrix& rOutput,
                                            const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    // Get perturbation size
    const double delta = this->GetPerturbationSize(rDesignVariable, rCurrentProcessInfo);
    const SizeType number_of_nodes = this->mpPrimalElement->GetGeometry().PointsNumber();
    const SizeType dimension = rCurrentProcessInfo.GetValue(DOMAIN_SIZE);
    const SizeType num_dofs_per_node = this->mHasRotationDofs ?  2 * dimension : dimension;
    const SizeType local_size = number_of_nodes * num_dofs_per_node;

    if( rDesignVariable == TEMPERATURE )
    {
        rOutput.resize( number_of_nodes, local_size, false);

        Vector RHS;
        Vector derived_RHS;
        
        this->pGetPrimalElement()->CalculateRightHandSide(RHS, rCurrentProcessInfo);

        for(IndexType i_node = 0; i_node < this->mpPrimalElement->GetGeometry().PointsNumber(); ++i_node)
        {
            // Get pseudo-load contribution from utility
            FiniteDifferenceUtility::CalculateRightHandSideDerivative(*(this->pGetPrimalElement()), RHS, rDesignVariable,
                                                                      this->mpPrimalElement->GetGeometry()[i_node], delta, 
                                                                      derived_RHS, rCurrentProcessInfo);

            KRATOS_ERROR_IF_NOT(derived_RHS.size() == local_size) << "Size of the pseudo-load does not fit! [ derived_RHS.size() = " << derived_RHS.size() << ", local_size = " << local_size << " ]." << std::endl;

            for(IndexType i = 0; i < derived_RHS.size(); ++i)
                rOutput(i_node, i) = derived_RHS[i];
        }
    }
    // else if( rDesignVariable == TRUSS_PRESTRESS_PK2 )
    // {
    //     Vector RHS;
    //     this->pGetPrimalElement()->CalculateRightHandSide(RHS, rCurrentProcessInfo);

    //     // Get pseudo-load from utility
    //     //FiniteDifferenceUtility::CalculateRightHandSideDerivative(*pGetPrimalElement(), RHS, rDesignVariable, delta, rOutput, rCurrentProcessInfo);
    //     rOutput.resize(1,RHS.size(), false);

    //     std::stringstream filename;
    //     filename << "sensitivity_element_" << this->mpPrimalElement->Id() << ".dat";

    //     std::ifstream inFile(filename.str());
    //     if (!inFile.is_open()) {
    //         KRATOS_ERROR << "Could not open file " << filename.str() << " for reading" << std::endl;
    //     }

    //     std::vector<double> values;
    //     double val;
    //     while (inFile >> val) {
    //         values.push_back(val);
    //     }
    //     inFile.close();

    //     // Create 1-row matrix
        
    //     for (std::size_t j = 0; j < values.size(); ++j) {
    //         rOutput(0, j) = values[j];
    //     }

    
    // }
    else
    {
        Vector RHS;
        this->pGetPrimalElement()->CalculateRightHandSide(RHS, rCurrentProcessInfo);
        //KRATOS_WATCH(RHS)
        //KRATOS_WATCH(this->pGetPrimalElement())

        // Get pseudo-load from utility
        FiniteDifferenceUtility::CalculateRightHandSideDerivative(*(this->pGetPrimalElement()), RHS, rDesignVariable, delta, rOutput, rCurrentProcessInfo);
    }

    if (rOutput.size1() == 0 || rOutput.size2() == 0)
    {
        //std::cout << "Passing 0 of size 1 to sensitivity matrix for variable: " << rDesignVariable.Name() << std::endl;
        rOutput = ZeroMatrix(1, local_size);
    }

    KRATOS_CATCH("")
}

template <class TPrimalElement>
void AdjointFiniteDifferenceCableElement<TPrimalElement>::CalculateSensitivityMatrix(const Variable<array_1d<double,3>>& rDesignVariable, Matrix& rOutput,
                                            const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    
    const SizeType number_of_nodes = this->mpPrimalElement->GetGeometry().PointsNumber();
    const SizeType dimension = rCurrentProcessInfo.GetValue(DOMAIN_SIZE);
    const SizeType num_dofs_per_node = (this->mHasRotationDofs) ?  2 * dimension : dimension;
    const SizeType local_size = number_of_nodes * num_dofs_per_node;
    
    if( rDesignVariable == SHAPE_SENSITIVITY )
    {
        const double delta = this->GetPerturbationSize(rDesignVariable, rCurrentProcessInfo);
        const std::vector<const FiniteDifferenceUtility::array_1d_component_type*> coord_directions = {&SHAPE_SENSITIVITY_X, &SHAPE_SENSITIVITY_Y, &SHAPE_SENSITIVITY_Z};
        Vector derived_RHS;

        if ( (rOutput.size1() != dimension * number_of_nodes) || (rOutput.size2() != local_size ) )
            rOutput.resize(dimension * number_of_nodes, local_size, false);

        IndexType index = 0;

        Vector RHS;
        this->pGetPrimalElement()->CalculateRightHandSide(RHS, rCurrentProcessInfo);
        for(auto& node_i : this->mpPrimalElement->GetGeometry())
        {
            for(IndexType coord_dir_i = 0; coord_dir_i < dimension; ++coord_dir_i)
            {
                // Get pseudo-load contribution from utility
                FiniteDifferenceUtility::CalculateRightHandSideDerivative(*(this->pGetPrimalElement()), RHS, *coord_directions[coord_dir_i],
                                                                            node_i, delta, derived_RHS, rCurrentProcessInfo);

                KRATOS_ERROR_IF_NOT(derived_RHS.size() == local_size) << "Size of the pseudo-load does not fit!" << std::endl;

                for(IndexType i = 0; i < derived_RHS.size(); ++i)
                    rOutput( (coord_dir_i + index*dimension), i) = derived_RHS[i];
            }
            index++;
        }
    }
    else if( rDesignVariable == PRE_STRESS )
    {
        Vector delta_vector;
        if (this->mpPrimalElement->pGetProperties()->Has(PRESTRESS_VECTOR)){

            delta_vector = this->GetPerturbationSize(PRESTRESS_VECTOR, rCurrentProcessInfo);
        }
        else{
            delta_vector = ZeroVector(3);
        }
        //KRATOS_WATCH(delta_vector)
        //const std::vector<const FiniteDifferenceUtility::array_1d_component_type*> components = {&PRE_STRESS_SENSITIVITY_XX, &PRE_STRESS_SENSITIVITY_YY, &PRE_STRESS_SENSITIVITY_XY};
        Matrix derived_RHS;
        const SizeType dimension = rCurrentProcessInfo.GetValue(DOMAIN_SIZE);
        
        KRATOS_ERROR_IF_NOT(dimension > 1) << "CalculateSensitivityMatrix for Vector variables is only available for 2 and 3D!" << std::endl;
        
        SizeType vector_size;
        if (dimension == 2){
            vector_size = 1;
        }
        else if (dimension == 3){
            vector_size = 3;
        }

        if ( (rOutput.size1() != vector_size) || (rOutput.size2() != local_size ) )
            rOutput.resize(vector_size , local_size, false);

        IndexType index = 0;

        Vector RHS;
        this->pGetPrimalElement()->CalculateRightHandSide(RHS, rCurrentProcessInfo);
        
        // Get pseudo-load contribution from utility
        FiniteDifferenceUtility::CalculateRightHandSideDerivative(*(this->pGetPrimalElement()), RHS, PRESTRESS_VECTOR,
                                                                    delta_vector, derived_RHS, rCurrentProcessInfo);

        KRATOS_ERROR_IF_NOT(derived_RHS.size2() == local_size) << "Size of the pseudo-load does not fit!" << std::endl;

        for(IndexType i = 0; i < derived_RHS.size1(); ++i)
            for(IndexType j = 0; j < derived_RHS.size2(); ++j)
                if (i==2){
                    rOutput(i,j) = 0.0;
                }
                else {
                    rOutput( i , j) = derived_RHS(i,j);
                }

    }
    else{
        //std::cout << "Passing 0 of size 3 to sensitivity matrix for variable: " << rDesignVariable.Name() << std::endl;
        rOutput = ZeroMatrix(3, local_size);
    }
    KRATOS_CATCH("")
}

// template <class TPrimalElement>
// void AdjointFiniteDifferenceCableElement<TPrimalElement>::CalculateOnIntegrationPoints(const Variable<array_1d<double, 3 > >& rVariable,
// 					      std::vector< array_1d<double, 3 > >& rOutput,
// 					      const ProcessInfo& rCurrentProcessInfo)
// {
//     KRATOS_TRY

//     if (rVariable == ADJOINT_STRAIN) {
//         std::vector<Vector> strain_vector;
//         this->CalculateAdjointFieldOnIntegrationPoints(STRAIN, strain_vector, rCurrentProcessInfo);
//         if (rOutput.size() != strain_vector.size()) {
//             rOutput.resize(strain_vector.size());
//         }

//         KRATOS_ERROR_IF(strain_vector[0].size() != 3) << "Dimension of strain vector not as expected!" << std::endl;

//         for(IndexType i = 0; i < strain_vector.size(); ++i) {
//             for (IndexType j = 0; j < 3 ; ++j) {
//                 rOutput[i][j] = strain_vector[i][j];
//             }
//         }
//     } else {
//         this->CalculateAdjointFieldOnIntegrationPoints(rVariable, rOutput, rCurrentProcessInfo);
//     }

//     KRATOS_CATCH("")
// }


// template <class TPrimalElement>
// void AdjointFiniteDifferenceCableElement<TPrimalElement>::CalculateStressDisplacementDerivative(const Variable<Vector>& rStressVariable,
//                                     Matrix& rOutput, const ProcessInfo& rCurrentProcessInfo)
// {
//     KRATOS_TRY;

//     TracedStressType traced_stress_type = static_cast<TracedStressType>(this->GetValue(TRACED_STRESS_TYPE));

//     if(traced_stress_type  == TracedStressType::FX)
//     {
//         // ensure that adjoint load is determined without influence of pre-stress
//         // pre-stress does not cancel out when computing this derivative with unit-displacements!
//         Properties::Pointer p_global_properties = this->mpPrimalElement->pGetProperties();

//         Properties::Pointer p_local_property(Kratos::make_shared<Properties>(*p_global_properties));
//         this->mpPrimalElement->SetProperties(p_local_property);

//         p_local_property->SetValue(TRUSS_PRESTRESS_PK2, 0.0);

//         AdjointFiniteDifferencingBaseElement<TPrimalElement>::CalculateStressDisplacementDerivative(rStressVariable,
//                                            rOutput, rCurrentProcessInfo);

//         this->mpPrimalElement->SetProperties(p_global_properties);
//     }
//     else
//         AdjointFiniteDifferencingBaseElement<TPrimalElement>::CalculateStressDisplacementDerivative(rStressVariable,
//                                    rOutput, rCurrentProcessInfo);

//     KRATOS_CATCH("")
// }

template <class TPrimalElement>
void AdjointFiniteDifferenceCableElement<TPrimalElement>::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
                                       const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;
    // Matrix primal_lhs;
    // this->mpPrimalElement->CalculateLeftHandSide(primal_lhs,
    //                                            rCurrentProcessInfo);

    //const double delta = this->GetPerturbationSize(rDesignVariable, rCurrentProcessInfo);
    const double delta = 1e-8;

    const SizeType number_of_nodes = this->mpPrimalElement->GetGeometry().PointsNumber();
    const SizeType dimension = rCurrentProcessInfo.GetValue(DOMAIN_SIZE);
    const SizeType num_dofs_per_node = (this->mHasRotationDofs) ?  2 * dimension : dimension;
    const SizeType local_size = number_of_nodes * num_dofs_per_node;

    Vector RHS;
    this->mpPrimalElement->CalculateRightHandSide(RHS, rCurrentProcessInfo);

    //const std::vector<const FiniteDifferenceUtility::array_1d_component_type*> coord_directions = {&SHAPE_SENSITIVITY_X, &SHAPE_SENSITIVITY_Y, &SHAPE_SENSITIVITY_Z};
    // Build vector of variables containing the DOF-variables of the primal problem
    std::vector<const Variable<double>*> primal_solution_variable_list {&DISPLACEMENT_X, &DISPLACEMENT_Y, &DISPLACEMENT_Z};

    // Store primal results and initialize deformation
    Vector initial_state_variables;
    initial_state_variables.resize(local_size, false);

    for (IndexType i = 0; i < number_of_nodes; ++i)
    {
        const IndexType index = i * num_dofs_per_node;
        for(IndexType j = 0; j < primal_solution_variable_list.size(); ++j)
        {
            initial_state_variables[index + j] = this->mpPrimalElement->GetGeometry()[i].FastGetSolutionStepValue(*primal_solution_variable_list[j]);
            //this->mpPrimalElement->GetGeometry()[i].FastGetSolutionStepValue(*primal_solution_variable_list[j]) = 0.0;
        }
    }

    // Compute gradient using unit deformation states
    rLeftHandSideMatrix.resize(local_size, local_size, false);
    rLeftHandSideMatrix.clear();

    //std::vector<Matrix> partial_stress_derivatives;

    for (IndexType i = 0; i < number_of_nodes; ++i)
    {
        const IndexType index = i * num_dofs_per_node;

        for(IndexType j = 0; j < primal_solution_variable_list.size(); ++j)
        {
            this->mpPrimalElement->GetGeometry()[i].FastGetSolutionStepValue(*primal_solution_variable_list[j]) = initial_state_variables[index + j] + delta;
            //this->mpPrimalElement->CalculateOnIntegrationPoints(PK2_STRESS_TENSOR, partial_stress_derivatives, rCurrentProcessInfo);
            Vector rhs_perturbed;
            this->mpPrimalElement->CalculateRightHandSide(rhs_perturbed, rCurrentProcessInfo);

            Vector sensitivity = (rhs_perturbed - RHS)/ delta;

            for(IndexType k = 0; k < RHS.size(); ++k)
            {
                rLeftHandSideMatrix((index + j), k)= -1.0* sensitivity[k]; // -1 becasue in RHS K is on right side with a minus, so R2-R1 is actually K1-K2
            }

            this->mpPrimalElement->GetGeometry()[i].FastGetSolutionStepValue(*primal_solution_variable_list[j]) = initial_state_variables[index + j] ;
        }
    }

    //compare if the primal_lhs is different from partial_R/partial_u
    // Matrix diff;
    // noalias(diff) = primal_lhs - rLeftHandSideMatrix;
    // std::cout << "difference between lhs and fd partial R/partial_u " << diff << std::endl;

    // // Recall primal solution
    // for (IndexType i = 0; i < number_of_nodes; ++i)
    // {
    //     const IndexType index = i * num_dofs_per_node;
    //     for(IndexType j = 0; j < primal_solution_variable_list.size(); ++j)
    //         this->mpPrimalElement->GetGeometry()[i].FastGetSolutionStepValue(*primal_solution_variable_list[j]) = initial_state_variables[index + j];
    // }

    // Vector derived_RHS;

    // // if ( (rOutput.size1() != dimension * number_of_nodes) || (rOutput.size2() != local_size ) )
    // //     rOutput.resize(dimension * number_of_nodes, local_size, false);

    // // IndexType index = 0;

    // Vector RHS;
    // this->mpPrimalElement->CalculateRightHandSide(RHS, rCurrentProcessInfo);
    // // for(auto& node_i : mpPrimalElement->GetGeometry())
    // // {
    // //     for(IndexType coord_dir_i = 0; coord_dir_i < dimension; ++coord_dir_i)
    // //     {
    // //         // Get pseudo-load contribution from utility
    // //         FiniteDifferenceUtility::CalculateRightHandSideDerivative(*(this->mpPrimalElement), RHS, *coord_directions[coord_dir_i],
    // //                                                                     node_i, delta, derived_RHS, rCurrentProcessInfo);

    // //         KRATOS_ERROR_IF_NOT(derived_RHS.size() == local_size) << "Size of the pseudo-load does not fit!" << std::endl;

    // //         for(IndexType i = 0; i < derived_RHS.size(); ++i)
    // //             rOutput( (coord_dir_i + index*dimension), i) = derived_RHS[i];
    // //     }
    // //     index++;
    // // }
   
    // rLeftHandSideMatrix = ZeroMatrix(0, local_size);

    KRATOS_CATCH("")
}

template <class TPrimalElement>
void AdjointFiniteDifferenceCableElement<TPrimalElement>::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseType);
}

template <class TPrimalElement>
void AdjointFiniteDifferenceCableElement<TPrimalElement>::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseType);
}

template class AdjointFiniteDifferenceCableElement<CableElement3D2N>;

} // namespace Kratos.



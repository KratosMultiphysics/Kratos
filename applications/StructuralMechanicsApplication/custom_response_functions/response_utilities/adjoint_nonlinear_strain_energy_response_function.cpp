// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Mahmoud Sesa, https://github.com/mahmoudsesa
//

// System includes

// External includes
#include "utilities/variable_utils.h"

// Project includes
#include "adjoint_nonlinear_strain_energy_response_function.h"

namespace Kratos
{
    AdjointNonlinearStrainEnergyResponseFunction::AdjointNonlinearStrainEnergyResponseFunction(ModelPart& rModelPart, Parameters ResponseSettings)
    : AdjointStructuralResponseFunction(rModelPart, ResponseSettings)
    {
        for(auto cond_it = rModelPart.ConditionsBegin(); cond_it != rModelPart.ConditionsEnd(); ++cond_it)
        {
            mConditions[cond_it->Id()] = cond_it;
            SizeType number_of_nodes = cond_it->GetGeometry().size();
            SizeType dimension = cond_it->GetGeometry().WorkingSpaceDimension();
            SizeType vec_size = number_of_nodes * dimension;
            mConditionsRHS[cond_it->Id()] = ZeroVector(vec_size);
            mResponseGradient_0[cond_it->Id()] = ZeroVector(vec_size);
        }
    }

    AdjointNonlinearStrainEnergyResponseFunction::~AdjointNonlinearStrainEnergyResponseFunction()
    {
    }

    double AdjointNonlinearStrainEnergyResponseFunction::CalculateValue(ModelPart& rModelPart)
    {
        KRATOS_TRY;

        this->CheckForBodyForces(rModelPart);

        ProcessInfo &r_current_process_info = rModelPart.GetProcessInfo();
        double response_increment_value = 0.0;

        // Check if there adjoint elements, because response function calculation is done for primal analysis
        KRATOS_ERROR_IF( r_current_process_info.Has(IS_ADJOINT) && r_current_process_info[IS_ADJOINT] )
        << "Calculate value for strain energy response is only available when using primal elements" << std::endl;

        // sum all elemental strain energy increment values calculated by trapezoidal rule: E = 0.5 * (f_ext_i - f_ext_i-1) * (u_i - u_i-1)


        //const int num_elem = static_cast<int> (rModelPart.NumberOfElements());

        // TODO Mahmoud: Calculation using the exact value for the external force at the last time step, not just an approximation
        #pragma omp parallel
        {
        Matrix LHS;
        Vector RHS;
        Vector disp;
        Vector external_force;
        Vector external_force_previous_step;

        Vector disp_previous_step;
        Vector disp_increment;
        Vector average_load;

        #pragma omp for reduction(+:response_increment_value)
        for (int i = 0; i< static_cast<int> (rModelPart.NumberOfElements()); i++)
        {
            auto elem_it = rModelPart.ElementsBegin() + i;
            elem_it->GetValuesVector(disp,0);
            elem_it->GetValuesVector(disp_previous_step, 1);
            elem_it->CalculateLocalSystem(LHS, RHS, r_current_process_info);

            disp_increment = disp - disp_previous_step;
            external_force = -1.0 * RHS;
            external_force_previous_step = external_force - prod(LHS , disp_increment);
            average_load = 0.5 * (external_force + external_force_previous_step);

            response_increment_value += inner_prod(average_load , disp_increment);
        }
        }
        return response_increment_value;
        KRATOS_CATCH("");
    }

    void AdjointNonlinearStrainEnergyResponseFunction::CalculateGradient(const Element& rAdjointElement,
                                    const Matrix& rResidualGradient,
                                    Vector& rResponseGradient,
                                    const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;

        if (rResponseGradient.size() != rResidualGradient.size1())
            rResponseGradient.resize(rResidualGradient.size1(), false);

        rResponseGradient.clear();

        KRATOS_CATCH("");
    }

    // \frac{1}{2} (f_{ext}_i - f_{ext}_i-1) + frac{1}{2} (u^T_i - u^T_i-1) \cdot ( \frac{\partial f_{ext}_i}{\partial u} + frac{\partial f_{ext}_i-1}{\partial u})
    void AdjointNonlinearStrainEnergyResponseFunction::CalculateGradient(const Condition& rAdjointCondition,
                                    const Matrix& rResidualGradient,
                                    Vector& rResponseGradient,
                                    const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;

        ProcessInfo process_info = rProcessInfo;

        const SizeType number_of_nodes = rAdjointCondition.GetGeometry().size();
        const SizeType dimension = rAdjointCondition.GetGeometry().WorkingSpaceDimension();
        const SizeType vec_size = number_of_nodes * dimension;

        Vector displacement = ZeroVector(vec_size);
        Vector displacement_previous_step = ZeroVector(vec_size);
        double delta = rAdjointCondition.GetValue(PERTURBATION_SIZE);
        Vector RHS;

        // getting the displacements for the current and the previous time steps
        int i = 0;
        for (auto& node_i : rAdjointCondition.GetGeometry())
        {
            project(displacement, range(i, dimension)) = rAdjointCondition.GetGeometry()[node_i].FastGetSolutionStepValue(DISPLACEMENT);
            project(displacement_previous_step , range(i , dimension))= rAdjointCondition.GetGeometry()[node_i].FastGetSolutionStepValue(DISPLACEMENT,1);
            i += 3;
        }

        mConditions[rAdjointCondition.Id()]->CalculateRightHandSide(RHS, process_info);

        //calculation of partial f_{ext}_i}{\partial u}
        Matrix partial_derivative_matrix = ZeroMatrix(vec_size, vec_size);
        int i_2 = 0;
        for (auto& node_i : rAdjointCondition.GetGeometry())
        {
            Vector perturbed_RHS = Vector(0);

            // Pertubation, gradient analysis and recovery of x
            node_i.X() += delta;
            node_i.FastGetSolutionStepValue(DISPLACEMENT_X) += delta;
            mConditions[rAdjointCondition.Id()]->CalculateRightHandSide(perturbed_RHS, process_info);
            row(partial_derivative_matrix, i_2) = (perturbed_RHS - RHS) / delta;
            node_i.X() -= delta;
            node_i.FastGetSolutionStepValue(DISPLACEMENT_X) -= delta;

            // Reset the pertubed vector
            perturbed_RHS = Vector(0);

            // Pertubation, gradient analysis and recovery of y
            node_i.Y() += delta;
            node_i.FastGetSolutionStepValue(DISPLACEMENT_Y) += delta;
            mConditions[rAdjointCondition.Id()]->CalculateRightHandSide(perturbed_RHS, process_info);
            row(partial_derivative_matrix, i_2 + 1) = (perturbed_RHS - RHS) / delta;
            node_i.Y() -= delta;
            node_i.FastGetSolutionStepValue(DISPLACEMENT_Y) -= delta;

            // Reset the pertubed vector
            perturbed_RHS = Vector(0);

            // Pertubation, gradient analysis and recovery of z
            node_i.Z() += delta;
            node_i.FastGetSolutionStepValue(DISPLACEMENT_Z) += delta;
            mConditions[rAdjointCondition.Id()]->CalculateRightHandSide(perturbed_RHS, process_info);
            row(partial_derivative_matrix, i_2 + 2) = (perturbed_RHS - RHS) / delta;
            node_i.Z() -= delta;
            node_i.FastGetSolutionStepValue(DISPLACEMENT_Z) -= delta;

            i_2 += 3;
        }

        if(mExternalForceDisplacementDerivative.size1() == 0)
            mExternalForceDisplacementDerivative = ZeroMatrix(vec_size, vec_size);

        // Summing up the partial derivative terms
        rResponseGradient = 0.50 * (RHS + mConditionsRHS[rAdjointCondition.Id()] +
                            prod(displacement - displacement_previous_step, mExternalForceDisplacementDerivative + partial_derivative_matrix));

        // storing the RHS for each condition
        mConditionsRHS[rAdjointCondition.Id()] = RHS;

        // storing the partial derivative partial f_{ext}_i}{\partial u}
        mExternalForceDisplacementDerivative = partial_derivative_matrix;

        // This stores rResponseGradient values for each condition
        mResponseGradient_1[rAdjointCondition.Id()] = rResponseGradient;

        KRATOS_CATCH("");
    }

    // TODO Mahmoud: This is not the correct place to implement the scaling factor, another function should be added
    // to the base class to be used for calculating the scaling factor
    // This calculate a scaling factor for the current response gradient w.r.t the response gradient of the last step
    // Later this should be replaced with a more general approach that is suitable for follower loads
    // a scaling factor would give accurate results for non-follower loads, however it should be used carefully for follower loads
    // frac{\partial E_i}{\partial u_i} / frac{\partial E_i-1}{\partial u_i-1}
    void AdjointNonlinearStrainEnergyResponseFunction::CalculateFirstDerivativesGradient(const Condition& rAdjointCondition,
                                                   const Matrix& rResidualGradient,
                                                   Vector& rResponseGradient,
                                                   const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY

        Vector response_gradient_0;
        Vector response_gradient_1;
        response_gradient_0 = mResponseGradient_0[rAdjointCondition.Id()];
        response_gradient_1 = mResponseGradient_1[rAdjointCondition.Id()];
        rResponseGradient = ZeroVector(1);
        const double tolerance = 1e-5;

        // Here a scalar value is calculated (lambda_1 / lambda_0)
        for (IndexType i = 0; i < response_gradient_1.size(); i++)
        {
            if( std::abs(response_gradient_0[i]) > tolerance )
            {
                KRATOS_WATCH(response_gradient_0[i])
                rResponseGradient[0] = response_gradient_1[i] / response_gradient_0[i];
                break;
            }
        }

        mResponseGradient_0[rAdjointCondition.Id()] = response_gradient_1;
        KRATOS_CATCH("");
    }

    void AdjointNonlinearStrainEnergyResponseFunction::CalculatePartialSensitivity(Element& rAdjointElement,
                                                const Variable<array_1d<double, 3>>& rVariable,
                                                const Matrix& rSensitivityMatrix,
                                                Vector& rSensitivityGradient,
                                                const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY

        if (rSensitivityGradient.size() != 0)
            rSensitivityGradient.resize(0, false);

        KRATOS_CATCH("");
    }


    void AdjointNonlinearStrainEnergyResponseFunction::CalculatePartialSensitivity(Element& rAdjointElement,
                                                const Variable<double>& rVariable,
                                                const Matrix& rSensitivityMatrix,
                                                Vector& rSensitivityGradient,
                                                const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;

        if (rSensitivityGradient.size() != 0)
            rSensitivityGradient.resize(0, false);

        KRATOS_CATCH("");
    }

    void AdjointNonlinearStrainEnergyResponseFunction::CalculatePartialSensitivity(Condition& rAdjointCondition,
                                                const Variable<double>& rVariable,
                                                const Matrix& rSensitivityMatrix,
                                                Vector& rSensitivityGradient,
                                                const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;

        if (rSensitivityGradient.size() != 0)
            rSensitivityGradient.resize(0, false);



        KRATOS_CATCH("");
    }


    // Calculates the derivate of the response function with respect to the design parameters
    // Computation of \frac{1}{2} (u^T_i - u^T_i-1)  \cdot ( \frac{\partial f_{ext}_i}{\partial x} + frac{\partial f_{ext}_i-1}{\partial x})
    void AdjointNonlinearStrainEnergyResponseFunction::CalculatePartialSensitivity(Condition& rAdjointCondition,
                                                const Variable<array_1d<double, 3>>& rVariable,
                                                const Matrix& rSensitivityMatrix,
                                                Vector& rSensitivityGradient,
                                                const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;
        // TODO Mahmoud: remove this assignment
        ProcessInfo process_info = rProcessInfo;

        const SizeType number_of_nodes = rAdjointCondition.GetGeometry().size();
        const SizeType dimension = rAdjointCondition.GetGeometry().WorkingSpaceDimension();
        const SizeType mat_size = number_of_nodes * dimension;

        // Matrix for partial f_{ext}_i}{\partial x}
        Matrix partial_derivative_matrix = ZeroMatrix(mat_size, mat_size);

        if(mExternalForceDesignVariableDerivative.size1() == 0)
            mExternalForceDesignVariableDerivative = ZeroMatrix(mat_size, mat_size);

        if (rSensitivityGradient.size() != mat_size)
            rSensitivityGradient.resize(mat_size, false);

        double Delta = rAdjointCondition.GetValue(PERTURBATION_SIZE);

        Vector displacement = ZeroVector(mat_size);
        Vector displacement_previous_step = ZeroVector(mat_size);
        Vector RHS;
        Matrix LHS;

        // Accessing nodal displacements at each node for the previous the current and the previous time steps
        int i_1 = 0;
        for (auto& node_i : rAdjointCondition.GetGeometry())
        {
            project(displacement, range(i_1 , i_1 + dimension)) = rAdjointCondition.GetGeometry()[node_i].FastGetSolutionStepValue(DISPLACEMENT);
            project(displacement_previous_step , range(i_1 , i_1 + dimension))= rAdjointCondition.GetGeometry()[node_i].FastGetSolutionStepValue(DISPLACEMENT,1);
            i_1 += 3;
        }

        rAdjointCondition.CalculateRightHandSide(RHS, process_info);

        // TODO, use rVariable instead of the nodal degrees of freedom
        int i_2 = 0;
        for (auto& node_i : rAdjointCondition.GetGeometry())
        {
            Vector perturbed_RHS = Vector(0);

            // Pertubation, gradient analysis and recovery of x
            node_i.X0() += Delta;
            node_i.X() += Delta;
            rAdjointCondition.CalculateRightHandSide(perturbed_RHS, process_info);
            row(partial_derivative_matrix, i_2) = (perturbed_RHS - RHS) / Delta;
            node_i.X0() -= Delta;
            node_i.X() -= Delta;

            // Reset pertubed vector
            perturbed_RHS = Vector(0);

            // Pertubation, gradient analysis and recovery of y
            node_i.Y0() += Delta;
            node_i.Y() += Delta;
            rAdjointCondition.CalculateRightHandSide(perturbed_RHS, process_info);
            row(partial_derivative_matrix, i_2 + 1) = (perturbed_RHS - RHS) / Delta;
            node_i.Y0() -= Delta;
            node_i.Y() -= Delta;

            // Reset pertubed vector
            perturbed_RHS = Vector(0);

            // Pertubation, gradient analysis and recovery of z
            node_i.Z0() += Delta;
            node_i.Z() += Delta;
            rAdjointCondition.CalculateRightHandSide(perturbed_RHS, process_info);
            row(partial_derivative_matrix, i_2 + 2) = (perturbed_RHS - RHS) / Delta;
            node_i.Z0() -= Delta;
            node_i.Z() -= Delta;

            i_2 += 3;
        }

        // Summing up the terms for the sensitivity calculation
        rSensitivityGradient = 0.50 * prod(displacement - displacement_previous_step , mExternalForceDesignVariableDerivative + partial_derivative_matrix);

        // passing the value of the partial derivative matrix to a member variable
        mExternalForceDesignVariableDerivative = partial_derivative_matrix;

        KRATOS_CATCH("");
    }

    void AdjointNonlinearStrainEnergyResponseFunction::CheckForBodyForces(ModelPart& rModelPart)
    {
        const double numerical_limit = std::numeric_limits<double>::epsilon();

        for (auto& elem_i : rModelPart.Elements())
        {
        Vector acc = ZeroVector(3);
        if (elem_i.GetProperties().Has( VOLUME_ACCELERATION ))
            acc+= elem_i.GetProperties()[VOLUME_ACCELERATION];

        if( elem_i.GetGeometry()[0].SolutionStepsDataHas(VOLUME_ACCELERATION))
            acc += elem_i.GetGeometry()[0].FastGetSolutionStepValue(VOLUME_ACCELERATION);

        KRATOS_ERROR_IF( norm_2(acc)>numerical_limit )
            << "Calculation of nonlinear strain energy response is untested for self weight!" << std::endl;
        }
    }

} // namespace Kratos
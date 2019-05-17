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
        mpModelPart = &rModelPart;
        for(auto cond_it = rModelPart.ConditionsBegin(); cond_it != rModelPart.ConditionsEnd(); ++cond_it)
        {
            // ToDO Mahmoud: the size of condition RHS is going to change for some condtion
            mConditions[cond_it->Id()] = cond_it;
            SizeType number_of_nodes = cond_it->GetGeometry().size();
            SizeType dimension = cond_it->GetGeometry().WorkingSpaceDimension();
            SizeType vec_size = number_of_nodes * dimension;
            mConditionsRHS[cond_it->Id()] = ZeroVector(vec_size);
            mResponseGradient_0[cond_it->Id()] = ZeroVector(vec_size);
            mResponseGradient_1[cond_it->Id()] = ZeroVector(vec_size);
        }
    }

    AdjointNonlinearStrainEnergyResponseFunction::~AdjointNonlinearStrainEnergyResponseFunction()
    {
    }

    void AdjointNonlinearStrainEnergyResponseFunction::FinalizeSolutionStep()
    {
        auto condition_pointer = mpModelPart->Conditions().begin();
        Matrix residual_gradient;
        Vector adjoint_values;
        ProcessInfo r_process_info = mpModelPart->GetProcessInfo();
        double load_factor_ratio = this->CalculateAdjointScalingFactor(*condition_pointer, residual_gradient, r_process_info);
        mpModelPart->GetProcessInfo().SetValue(ADJOINT_CORRECTION_FACTOR, load_factor_ratio);
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
        // TODO Mahmoud: Calculation using the exact value for the external force at the last time step, not just an approximation
        #pragma omp parallel
        {
        Matrix LHS;
        Vector RHS;
        Vector disp;
        Vector disp_0;
        Vector disp_increment;
        Vector average_load;

        #pragma omp for reduction(+:response_increment_value)
        for (int i = 0; i< static_cast<int> (rModelPart.NumberOfConditions()); i++)
        {
            auto cond_it = rModelPart.ConditionsBegin() + i;
            cond_it->GetValuesVector(disp,0);
            cond_it->GetValuesVector(disp_0,1);
            disp_increment = disp - disp_0;

            cond_it->CalculateRightHandSide(RHS, r_current_process_info);
            average_load = 0.5 * (RHS + mConditionsRHS[cond_it->Id()]);
            response_increment_value += inner_prod(average_load , disp_increment);
            mConditionsRHS[cond_it->Id()] = RHS;
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
    // TODO Mahmoud: The derivative for shape sensitivity different from element properties sensitivities,
    // equation for element sensitivities is: \frac{1}{2} (f_{ext}_i - f_{ext}_i-1)
    // TODO Mahmoud: function calculate gradient should have information about response variable as in CalculatePartialSensitivity()
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
        double delta = rProcessInfo[PERTURBATION_SIZE];
        Vector RHS;

        // This ensures that stored RHS equals zero at the initial time step
        if(rProcessInfo(STEP) == 1)
            mConditionsRHS[rAdjointCondition.Id()] = ZeroVector(vec_size);

        // getting the displacements for the current and the previous time steps
        int i = 0;
        for (auto& node_i : rAdjointCondition.GetGeometry())
        {
            for(IndexType j = 0; j < dimension; j++)
            {
                // TODO Mahmoud: This way returns zero displacement in case of surface loads, while for point loads it gives a value
                //displacement[i+j] = rAdjointCondition.GetGeometry()[node_i].FastGetSolutionStepValue(DISPLACEMENT)[j];
                //displacement_previous_step[i+j] = rAdjointCondition.GetGeometry()[node_i].FastGetSolutionStepValue(DISPLACEMENT,1)[j];
                displacement[i+j] = node_i.FastGetSolutionStepValue(DISPLACEMENT)[j];
                displacement_previous_step[i+j] = node_i.FastGetSolutionStepValue(DISPLACEMENT,1)[j];
            }
            i += 3;
        }

        mConditions[rAdjointCondition.Id()]->CalculateRightHandSide(RHS, process_info);

        // //calculation of partial f_{ext}_i}{\partial u}
        // //TODO Mahmoud: pertubing the displacement has no effect on the RHS
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

        if(mExternalForceDisplacementDerivative[rAdjointCondition.Id()].size1() == 0)
            mExternalForceDisplacementDerivative[rAdjointCondition.Id()] = ZeroMatrix(vec_size, vec_size);

        // Summing up the partial derivative terms
        rResponseGradient = 0.50 * (RHS + mConditionsRHS[rAdjointCondition.Id()]
                          + prod(displacement - displacement_previous_step,  mExternalForceDisplacementDerivative[rAdjointCondition.Id()] + partial_derivative_matrix ) );

        //rResponseGradient = 0.50 * (RHS + mConditionsRHS[rAdjointCondition.Id()]);


        // Storing the RHS for each condition
        mConditionsRHS[rAdjointCondition.Id()] = RHS;

        // storing the partial derivative partial f_{ext}_i}{\partial u}
       // mExternalForceDisplacementDerivative[rAdjointCondition.Id()] = partial_derivative_matrix;

        // This stores rResponseGradient values for each condition
        mResponseGradient_1[rAdjointCondition.Id()] = rResponseGradient;

        KRATOS_CATCH("");
    }

    // This calculate a scaling factor for the current response gradient w.r.t the response gradient of the last step
    // Later this should be replaced with a more general approach that is suitable for follower loads
    // a scaling factor would give accurate results for non-follower loads, however it should be used carefully for follower loads
    // frac{\partial E_i}{\partial u_i} / frac{\partial E_i-1}{\partial u_i-1}
    double AdjointNonlinearStrainEnergyResponseFunction::CalculateAdjointScalingFactor(const Condition& rAdjointCondition,
                                                   const Matrix& rResidualGradient,
                                                   const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY

        Vector response_gradient_0;
        Vector response_gradient_1;
        double l2_norm_1;
        double l2_norm_0;
        double scaling_factor = 0;
        response_gradient_0 = mResponseGradient_0[rAdjointCondition.Id()];
        response_gradient_1 = mResponseGradient_1[rAdjointCondition.Id()];
        const double tolerance = 1e-5;

        l2_norm_1 = norm_2(response_gradient_1);
        l2_norm_0 = norm_2(response_gradient_0);
        scaling_factor = l2_norm_1 / l2_norm_0;

        mResponseGradient_0[rAdjointCondition.Id()] = response_gradient_1;
        return scaling_factor;

        KRATOS_CATCH("");
    }

    void AdjointNonlinearStrainEnergyResponseFunction::CalculatePartialSensitivity(Element& rAdjointElement,
                                                const Variable<array_1d<double, 3>>& rVariable,
                                                const Matrix& rSensitivityMatrix,
                                                Vector& rSensitivityGradient,
                                                const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY

        rSensitivityGradient = ZeroVector(rSensitivityMatrix.size1());

        KRATOS_CATCH("");
    }


    void AdjointNonlinearStrainEnergyResponseFunction::CalculatePartialSensitivity(Element& rAdjointElement,
                                                const Variable<double>& rVariable,
                                                const Matrix& rSensitivityMatrix,
                                                Vector& rSensitivityGradient,
                                                const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;

        rSensitivityGradient = ZeroVector(rSensitivityMatrix.size1());

        KRATOS_CATCH("");
    }

    void AdjointNonlinearStrainEnergyResponseFunction::CalculatePartialSensitivity(Condition& rAdjointCondition,
                                                const Variable<double>& rVariable,
                                                const Matrix& rSensitivityMatrix,
                                                Vector& rSensitivityGradient,
                                                const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;

        rSensitivityGradient = ZeroVector(rSensitivityMatrix.size1());

        KRATOS_CATCH("");
    }


    // ToDo Mahmoud: this function only work for shape sensitivities, it should modified to work in general for
    // other types of sensitivity
    // Calculates the derivate of the response function with respect to the design parameters
    // Computation of  ( \frac{\partial f_{ext}_i}{\partial x} + frac{\partial f_{ext}_i-1}{\partial x})  \cdot \frac{1}{2} (u^T_i - u^T_i-1)
    void AdjointNonlinearStrainEnergyResponseFunction::CalculatePartialSensitivity(Condition& rAdjointCondition,
                                                const Variable<array_1d<double, 3>>& rVariable,
                                                const Matrix& rSensitivityMatrix,
                                                Vector& rSensitivityGradient,
                                                const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;

        KRATOS_ERROR_IF( rVariable != SHAPE_SENSITIVITY )
            << "Calculation of nonlinear strain energy sensitivity is available only for shape variables!" << std::endl;

        // TODO Mahmoud: This copy assignment is expensive
        ProcessInfo process_info = rProcessInfo;

        const SizeType number_of_nodes = rAdjointCondition.GetGeometry().size();
        const SizeType dimension = rAdjointCondition.GetGeometry().WorkingSpaceDimension();
        const SizeType mat_size = number_of_nodes * dimension;

        // Matrix for partial f_{ext}_i}{\partial x}
        Matrix partial_derivative_matrix = ZeroMatrix(mat_size, mat_size);

        if(mExternalForceDesignVariableDerivative[rAdjointCondition.Id()].size1() == 0)
            mExternalForceDesignVariableDerivative[rAdjointCondition.Id()] = ZeroMatrix(mat_size, mat_size);

        if (rSensitivityGradient.size() != mat_size)
            rSensitivityGradient.resize(mat_size, false);

        double Delta = rProcessInfo[PERTURBATION_SIZE];
        Vector displacement = ZeroVector(mat_size);
        Vector displacement_previous_step = ZeroVector(mat_size);
        Vector RHS;
        Matrix LHS;

        // Accessing nodal displacements at each node for the previous the current and the previous time steps
        int i_1 = 0;
        for (auto& node_i : rAdjointCondition.GetGeometry())
        {
            for(IndexType j = 0; j < dimension; j++)
            {
                // TODO Mahmoud: check why accessing nodal displacements by this way does not work for surface loads
                displacement[i_1+j] = node_i.FastGetSolutionStepValue(DISPLACEMENT)[j];
                displacement_previous_step[i_1+j] = node_i.FastGetSolutionStepValue(DISPLACEMENT,1)[j];
            }
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
        rSensitivityGradient = 0.50 * prod(mExternalForceDesignVariableDerivative[rAdjointCondition.Id()] + partial_derivative_matrix, displacement - displacement_previous_step );
        // passing the value of the partial derivative matrix to a member variable
        mExternalForceDesignVariableDerivative[rAdjointCondition.Id()] = partial_derivative_matrix;

        KRATOS_CATCH("");
    }

    void AdjointNonlinearStrainEnergyResponseFunction::CheckForBodyForces(ModelPart& rModelPart)
    {
        KRATOS_TRY;

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

        KRATOS_CATCH("");
    }

} // namespace Kratos
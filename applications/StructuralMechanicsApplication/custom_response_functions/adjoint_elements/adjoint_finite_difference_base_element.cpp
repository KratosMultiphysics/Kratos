// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Armin Geiser, https://github.com/armingeiser
//

// System includes

// External includes

// Project includes
#include "adjoint_finite_difference_base_element.h"
#include "custom_response_functions/response_utilities/stress_response_definitions.h"
#include "custom_response_functions/response_utilities/finite_difference_utility.h"
#include "includes/checks.h"
#include "custom_elements/shell_thin_element_3D3N.hpp"
#include "custom_elements/cr_beam_element_linear_3D2N.hpp"
#include "custom_elements/truss_element_3D2N.hpp"
#include "custom_elements/truss_element_linear_3D2N.hpp"
#include "custom_elements/small_displacement.h"
#include "custom_elements/spring_damper_element_3D2N.hpp"


namespace Kratos
{

template <class TPrimalElement>
void AdjointFiniteDifferencingBaseElement<TPrimalElement>::EquationIdVector(EquationIdVectorType& rResult,
    const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY
    const GeometryType& geom = this->GetGeometry();

    const SizeType number_of_nodes = geom.PointsNumber();
    const SizeType dimension = geom.WorkingSpaceDimension();
    const SizeType num_dofs_per_node = (mHasRotationDofs) ?  2 * dimension : dimension;
    const SizeType num_dofs = number_of_nodes * num_dofs_per_node;

    if(rResult.size() != num_dofs)
        rResult.resize(num_dofs, false);

    for(IndexType i = 0; i < geom.size(); ++i)
    {
        const IndexType index = i * num_dofs_per_node;
        const NodeType& iNode = geom[i];

        rResult[index]     = iNode.GetDof(ADJOINT_DISPLACEMENT_X).EquationId();
        rResult[index + 1] = iNode.GetDof(ADJOINT_DISPLACEMENT_Y).EquationId();
        rResult[index + 2] = iNode.GetDof(ADJOINT_DISPLACEMENT_Z).EquationId();

        if(mHasRotationDofs)
        {
            rResult[index + 3] = iNode.GetDof(ADJOINT_ROTATION_X).EquationId();
            rResult[index + 4] = iNode.GetDof(ADJOINT_ROTATION_Y).EquationId();
            rResult[index + 5] = iNode.GetDof(ADJOINT_ROTATION_Z).EquationId();
        }
    }
    KRATOS_CATCH("")
}

template <class TPrimalElement>
void AdjointFiniteDifferencingBaseElement<TPrimalElement>::GetDofList(DofsVectorType& rElementalDofList,
    const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    const GeometryType & geom = this->GetGeometry();

    const SizeType number_of_nodes = geom.PointsNumber();
    const SizeType dimension =  geom.WorkingSpaceDimension();
    const SizeType num_dofs_per_node = (mHasRotationDofs) ?  2 * dimension : dimension;
    const SizeType num_dofs = number_of_nodes * num_dofs_per_node;

    if (rElementalDofList.size() != num_dofs)
        rElementalDofList.resize(num_dofs);

    for (IndexType i = 0; i < number_of_nodes; ++i)
    {
        const IndexType index = i * num_dofs_per_node ;
        rElementalDofList[index    ] = this->GetGeometry()[i].pGetDof(ADJOINT_DISPLACEMENT_X);
        rElementalDofList[index + 1] = this->GetGeometry()[i].pGetDof(ADJOINT_DISPLACEMENT_Y);
        rElementalDofList[index + 2] = this->GetGeometry()[i].pGetDof(ADJOINT_DISPLACEMENT_Z);

        if(mHasRotationDofs)
        {
            rElementalDofList[index + 3] = this->GetGeometry()[i].pGetDof(ADJOINT_ROTATION_X);
            rElementalDofList[index + 4] = this->GetGeometry()[i].pGetDof(ADJOINT_ROTATION_Y);
            rElementalDofList[index + 5] = this->GetGeometry()[i].pGetDof(ADJOINT_ROTATION_Z);
        }
    }
    KRATOS_CATCH("")
}

template <class TPrimalElement>
void AdjointFiniteDifferencingBaseElement<TPrimalElement>::GetValuesVector(Vector& rValues, int Step) const
{
    KRATOS_TRY

    const GeometryType & geom = this->GetGeometry();

    const SizeType number_of_nodes = geom.PointsNumber();
    const SizeType dimension =  geom.WorkingSpaceDimension();
    const SizeType num_dofs_per_node = (mHasRotationDofs) ?  2 * dimension : dimension;
    const SizeType num_dofs = number_of_nodes * num_dofs_per_node;

    if(rValues.size() != num_dofs)
        rValues.resize(num_dofs, false);

    for (IndexType i = 0; i < number_of_nodes; ++i)
    {
        const NodeType & iNode = geom[i];
        const array_1d<double,3>& disp = iNode.FastGetSolutionStepValue(ADJOINT_DISPLACEMENT, Step);

        const IndexType index = i * num_dofs_per_node;
        rValues[index]     = disp[0];
        rValues[index + 1] = disp[1];
        rValues[index + 2] = disp[2];

        if(mHasRotationDofs)
        {
            const array_1d<double,3>& rot = iNode.FastGetSolutionStepValue(ADJOINT_ROTATION, Step);
            rValues[index + 3] = rot[0];
            rValues[index + 4] = rot[1];
            rValues[index + 5] = rot[2];
        }
    }
    KRATOS_CATCH("")
}

template <class TPrimalElement>
void AdjointFiniteDifferencingBaseElement<TPrimalElement>::Calculate(const Variable<Matrix >& rVariable, Matrix& rOutput,
                                                const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    if(rVariable == STRESS_DISP_DERIV_ON_GP)
    {
        this->CalculateStressDisplacementDerivative(STRESS_ON_GP, rOutput, rCurrentProcessInfo);
    }
    else if(rVariable == STRESS_DISP_DERIV_ON_NODE)
    {
        this->CalculateStressDisplacementDerivative(STRESS_ON_NODE, rOutput, rCurrentProcessInfo);
    }
    else if(rVariable == STRESS_DESIGN_DERIVATIVE_ON_GP)
    {
        const std::string& design_variable_name = this->GetValue( DESIGN_VARIABLE_NAME );

        if (KratosComponents<Variable<double>>::Has(design_variable_name))
        {
            const Variable<double>& r_variable =
                KratosComponents<Variable<double>>::Get(design_variable_name);
            this->CalculateStressDesignVariableDerivative(r_variable, STRESS_ON_GP, rOutput, rCurrentProcessInfo);
        }
        else if (KratosComponents<Variable<array_1d<double, 3>>>::Has(design_variable_name))
        {
            const Variable<array_1d<double, 3>>& r_variable =
                KratosComponents<Variable<array_1d<double, 3>>>::Get(design_variable_name);
            this->CalculateStressDesignVariableDerivative(r_variable, STRESS_ON_GP, rOutput, rCurrentProcessInfo);
        }
    }
    else if(rVariable == STRESS_DESIGN_DERIVATIVE_ON_NODE)
    {
        std::string& design_variable_name = this->GetValue( DESIGN_VARIABLE_NAME );

        if (KratosComponents<Variable<double>>::Has(design_variable_name))
        {
            const Variable<double>& r_variable =
                KratosComponents<Variable<double>>::Get(design_variable_name);
            this->CalculateStressDesignVariableDerivative(r_variable, STRESS_ON_NODE, rOutput, rCurrentProcessInfo);
        }
        else if (KratosComponents<Variable<array_1d<double, 3>>>::Has(design_variable_name))
        {
            const Variable<array_1d<double, 3>>& r_variable =
                KratosComponents<Variable<array_1d<double, 3>>>::Get(design_variable_name);
            this->CalculateStressDesignVariableDerivative(r_variable, STRESS_ON_NODE, rOutput, rCurrentProcessInfo);
        }
    }
    else if(rVariable == LOCAL_ELEMENT_ORIENTATION)
    {
        this->pGetPrimalElement()->Calculate(rVariable, rOutput, rCurrentProcessInfo);
    }
    else
    {
        KRATOS_WARNING("AdjointFiniteDifferencingBaseElement") << "Calculate function called for unknown variable: " << rVariable << std::endl;
        rOutput.clear();
    }

    KRATOS_CATCH("")
}

template <class TPrimalElement>
void AdjointFiniteDifferencingBaseElement<TPrimalElement>::CalculateOnIntegrationPoints(const Variable<double>& rVariable,
                    std::vector<double>& rValues,
                    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    if(this->Has(rVariable))
    {
        // Get result value for output
        const double& output_value = this->GetValue(rVariable);

        // Resize Output
        const SizeType  gauss_points_number = this->GetGeometry()
            .IntegrationPointsNumber(this->GetIntegrationMethod());
        if (rValues.size() != gauss_points_number)
            rValues.resize(gauss_points_number);

        // Write scalar result value on all Gauss-Points
        for(IndexType i = 0; i < gauss_points_number; ++i)
            rValues[i] = output_value;
    }
    else
        KRATOS_ERROR << "Unsupported output variable." << std::endl;

    KRATOS_CATCH("")
}

template <class TPrimalElement>
void AdjointFiniteDifferencingBaseElement<TPrimalElement>::CalculateOnIntegrationPoints(
        const Variable<array_1d<double, 3 > >& rVariable, std::vector< array_1d<double, 3 > >& rOutput, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    if(this->Has(rVariable)) {
        // Get result value for output
        const auto& output_value = this->GetValue(rVariable);

        // Resize Output
        const SizeType gauss_points_number = this->GetGeometry()
            .IntegrationPointsNumber(this->GetIntegrationMethod());
        if (rOutput.size() != gauss_points_number) {
            rOutput.resize(gauss_points_number);
        }

        // Write scalar result value on all Gauss-Points
        for(IndexType i = 0; i < gauss_points_number; ++i) {
            rOutput[i] = output_value;
        }

    } else {
        KRATOS_ERROR << "Unsupported output variable." << std::endl;
    }

    KRATOS_CATCH("")
}

template <class TPrimalElement>
int AdjointFiniteDifferencingBaseElement<TPrimalElement>::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    int return_value = Element::Check(rCurrentProcessInfo);

    KRATOS_ERROR_IF_NOT(mpPrimalElement) << "Primal element pointer is nullptr!" << std::endl;

    const GeometryType& r_geom = this->GetGeometry();

    // TODO generic way of doing these checks without checking the dofs..

    // Check dofs
    for (IndexType i = 0; i < r_geom.size(); ++i) {
        const auto& r_node = r_geom[i];

        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DISPLACEMENT, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(ADJOINT_DISPLACEMENT, r_node);

        KRATOS_CHECK_DOF_IN_NODE(ADJOINT_DISPLACEMENT_X, r_node);
        KRATOS_CHECK_DOF_IN_NODE(ADJOINT_DISPLACEMENT_Y, r_node);
        KRATOS_CHECK_DOF_IN_NODE(ADJOINT_DISPLACEMENT_Z, r_node);

        if(mHasRotationDofs) {
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(ROTATION, r_node);
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(ADJOINT_ROTATION, r_node);
            KRATOS_CHECK_DOF_IN_NODE(ADJOINT_ROTATION_X, r_node);
            KRATOS_CHECK_DOF_IN_NODE(ADJOINT_ROTATION_Y, r_node);
            KRATOS_CHECK_DOF_IN_NODE(ADJOINT_ROTATION_Z, r_node);
        }
    }

    return return_value;

    KRATOS_CATCH("")
}

// Sensitivity functions

template <class TPrimalElement>
void AdjointFiniteDifferencingBaseElement<TPrimalElement>::CalculateSensitivityMatrix(const Variable<double>& rDesignVariable, Matrix& rOutput,
                                            const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    // Get perturbation size
    const double delta = this->GetPerturbationSize(rDesignVariable, rCurrentProcessInfo);

    Vector RHS;
    this->pGetPrimalElement()->CalculateRightHandSide(RHS, rCurrentProcessInfo);

    // Get pseudo-load from utility
    FiniteDifferenceUtility::CalculateRightHandSideDerivative(*pGetPrimalElement(), RHS, rDesignVariable, delta, rOutput, rCurrentProcessInfo);

    if (rOutput.size1() == 0 || rOutput.size2() == 0)
    {
        const SizeType number_of_nodes = mpPrimalElement->GetGeometry().PointsNumber();
        const SizeType dimension = rCurrentProcessInfo.GetValue(DOMAIN_SIZE);
        const SizeType num_dofs_per_node = (mHasRotationDofs) ?  2 * dimension : dimension;
        const SizeType local_size = number_of_nodes * num_dofs_per_node;
        rOutput = ZeroMatrix(0, local_size);
    }

    KRATOS_CATCH("")
}

template <class TPrimalElement>
void AdjointFiniteDifferencingBaseElement<TPrimalElement>::CalculateSensitivityMatrix(const Variable<array_1d<double,3>>& rDesignVariable, Matrix& rOutput,
                                            const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    const double delta = this->GetPerturbationSize(rDesignVariable, rCurrentProcessInfo);

    const SizeType number_of_nodes = mpPrimalElement->GetGeometry().PointsNumber();
    const SizeType dimension = rCurrentProcessInfo.GetValue(DOMAIN_SIZE);
    const SizeType num_dofs_per_node = (mHasRotationDofs) ?  2 * dimension : dimension;
    const SizeType local_size = number_of_nodes * num_dofs_per_node;

    if( rDesignVariable == SHAPE_SENSITIVITY )
    {
        const std::vector<const FiniteDifferenceUtility::array_1d_component_type*> coord_directions = {&SHAPE_SENSITIVITY_X, &SHAPE_SENSITIVITY_Y, &SHAPE_SENSITIVITY_Z};
        Vector derived_RHS;

        if ( (rOutput.size1() != dimension * number_of_nodes) || (rOutput.size2() != local_size ) )
            rOutput.resize(dimension * number_of_nodes, local_size, false);

        IndexType index = 0;

        Vector RHS;
        pGetPrimalElement()->CalculateRightHandSide(RHS, rCurrentProcessInfo);
        for(auto& node_i : mpPrimalElement->GetGeometry())
        {
            for(IndexType coord_dir_i = 0; coord_dir_i < dimension; ++coord_dir_i)
            {
                // Get pseudo-load contribution from utility
                FiniteDifferenceUtility::CalculateRightHandSideDerivative(*pGetPrimalElement(), RHS, *coord_directions[coord_dir_i],
                                                                            node_i, delta, derived_RHS, rCurrentProcessInfo);

                KRATOS_ERROR_IF_NOT(derived_RHS.size() == local_size) << "Size of the pseudo-load does not fit!" << std::endl;

                for(IndexType i = 0; i < derived_RHS.size(); ++i)
                    rOutput( (coord_dir_i + index*dimension), i) = derived_RHS[i];
            }
            index++;
        }
    }
    else
        rOutput = ZeroMatrix(0, local_size);

    KRATOS_CATCH("")
}

template <class TPrimalElement>
void AdjointFiniteDifferencingBaseElement<TPrimalElement>::CalculateStressDisplacementDerivative(const Variable<Vector>& rStressVariable,
                                            Matrix& rOutput, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    const SizeType num_nodes = mpPrimalElement->GetGeometry().PointsNumber();
    const SizeType dimension = mpPrimalElement->GetGeometry().WorkingSpaceDimension();
    const SizeType num_dofs_per_node = (mHasRotationDofs) ?  2 * dimension : dimension;
    const SizeType num_dofs = num_nodes * num_dofs_per_node;
    Vector initial_state_variables;
    Vector stress_derivatives_vector;

    // TODO first calculation only to get the size of the stress vector
    TracedStressType traced_stress_type = static_cast<TracedStressType>(this->GetValue(TRACED_STRESS_TYPE));
    if (rStressVariable == STRESS_ON_GP)
        StressCalculation::CalculateStressOnGP(*pGetPrimalElement(), traced_stress_type, stress_derivatives_vector, rCurrentProcessInfo);
    else
        StressCalculation::CalculateStressOnNode(*pGetPrimalElement(), traced_stress_type, stress_derivatives_vector, rCurrentProcessInfo);
    rOutput.resize(num_dofs, stress_derivatives_vector.size(), false);
    rOutput.clear();
    initial_state_variables.resize(num_dofs, false);

    // Build vector of variables containing the DOF-variables of the primal problem
    std::vector<const Variable<double>*> primal_solution_variable_list;
    primal_solution_variable_list.reserve(num_dofs_per_node);
    primal_solution_variable_list.push_back(&DISPLACEMENT_X);
    primal_solution_variable_list.push_back(&DISPLACEMENT_Y);
    primal_solution_variable_list.push_back(&DISPLACEMENT_Z);

    if(mHasRotationDofs)
    {
        primal_solution_variable_list.push_back(&ROTATION_X);
        primal_solution_variable_list.push_back(&ROTATION_Y);
        primal_solution_variable_list.push_back(&ROTATION_Z);
    }

    // TODO Find a better way of doing this check
    KRATOS_ERROR_IF(rCurrentProcessInfo.Has(NL_ITERATION_NUMBER))
        << "This stress displacement derivative computation is only usable for linear cases!" << std::endl;

    for (IndexType i = 0; i < num_nodes; ++i)
    {
        const IndexType index = i * num_dofs_per_node;
        for(IndexType j = 0; j < primal_solution_variable_list.size(); ++j)
        {
            initial_state_variables[index + j] = mpPrimalElement->GetGeometry()[i].FastGetSolutionStepValue(*primal_solution_variable_list[j]);
            mpPrimalElement->GetGeometry()[i].FastGetSolutionStepValue(*primal_solution_variable_list[j]) = 0.0;
        }
    }
    for (IndexType i = 0; i < num_nodes; ++i)
    {
        const IndexType index = i * num_dofs_per_node;
        for(IndexType j = 0; j < primal_solution_variable_list.size(); ++j)
        {
            mpPrimalElement->GetGeometry()[i].FastGetSolutionStepValue(*primal_solution_variable_list[j]) = 1.0;

            TracedStressType traced_stress_type = static_cast<TracedStressType>(this->GetValue(TRACED_STRESS_TYPE));
            if (rStressVariable == STRESS_ON_GP)
                StressCalculation::CalculateStressOnGP(*pGetPrimalElement(), traced_stress_type, stress_derivatives_vector, rCurrentProcessInfo);
            else
                StressCalculation::CalculateStressOnNode(*pGetPrimalElement(), traced_stress_type, stress_derivatives_vector, rCurrentProcessInfo);

            for(IndexType k = 0; k < stress_derivatives_vector.size(); ++k)
                rOutput(index+j, k) = stress_derivatives_vector[k];

            stress_derivatives_vector.clear();

            mpPrimalElement->GetGeometry()[i].FastGetSolutionStepValue(*primal_solution_variable_list[j]) = 0.0;
        }
    }
    for (IndexType i = 0; i < num_nodes; ++i)
    {
        const IndexType index = i * num_dofs_per_node;
        for(IndexType j = 0; j < primal_solution_variable_list.size(); ++j)
            mpPrimalElement->GetGeometry()[i].FastGetSolutionStepValue(*primal_solution_variable_list[j]) = initial_state_variables[index + j];
    }

    KRATOS_CATCH("")
}

template <class TPrimalElement>
void AdjointFiniteDifferencingBaseElement<TPrimalElement>::CalculateStressDesignVariableDerivative(const Variable<double>& rDesignVariable,
                                                const Variable<Vector>& rStressVariable, Matrix& rOutput,
                                                const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    // Define working variables
    Vector stress_vector_undist;
    Vector stress_vector_dist;

    // Compute stress on GP before perturbation
    TracedStressType traced_stress_type = static_cast<TracedStressType>(this->GetValue(TRACED_STRESS_TYPE));
    if (rStressVariable == STRESS_ON_GP)
        StressCalculation::CalculateStressOnGP(*pGetPrimalElement(), traced_stress_type, stress_vector_undist, rCurrentProcessInfo);
    else
        StressCalculation::CalculateStressOnNode(*pGetPrimalElement(), traced_stress_type, stress_vector_undist, rCurrentProcessInfo);

    const SizeType stress_vector_size = stress_vector_undist.size();

    // Get perturbation size
    const double delta = this->GetPerturbationSize(rDesignVariable, rCurrentProcessInfo);

    if( mpPrimalElement->GetProperties().Has(rDesignVariable) )
    {
        rOutput.resize(1, stress_vector_size, false);

        // Save property pointer
        Properties::Pointer p_global_properties = mpPrimalElement->pGetProperties();

        // Create new property and assign it to the element
        Properties::Pointer p_local_property(Kratos::make_shared<Properties>(Properties(*p_global_properties)));
        mpPrimalElement->SetProperties(p_local_property);

        // perturb the design variable
        const double current_property_value = mpPrimalElement->GetProperties()[rDesignVariable];
        p_local_property->SetValue(rDesignVariable, (current_property_value + delta));

        // Compute stress on GP after perturbation
        TracedStressType traced_stress_type = static_cast<TracedStressType>(this->GetValue(TRACED_STRESS_TYPE));
        if (rStressVariable == STRESS_ON_GP)
            StressCalculation::CalculateStressOnGP(*pGetPrimalElement(), traced_stress_type, stress_vector_dist, rCurrentProcessInfo);
        else
            StressCalculation::CalculateStressOnNode(*pGetPrimalElement(), traced_stress_type, stress_vector_dist, rCurrentProcessInfo);

        // Compute derivative of stress w.r.t. design variable with finite differences
        for(size_t j = 0; j < stress_vector_size; ++j)
            rOutput(0, j) = (stress_vector_dist[j]-stress_vector_undist[j])/delta;

        // Give element original properties back
        mpPrimalElement->SetProperties(p_global_properties);
    }
    else
    {
        rOutput = ZeroMatrix(0, stress_vector_size);
    }

    KRATOS_CATCH("")
}

template <class TPrimalElement>
void AdjointFiniteDifferencingBaseElement<TPrimalElement>::CalculateStressDesignVariableDerivative(const Variable<array_1d<double,3>>& rDesignVariable,
                                            const Variable<Vector>& rStressVariable,
                                            Matrix& rOutput, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    // define working variables
    Vector stress_vector_undist;
    Vector stress_vector_dist;

    // Compute stress on GP before perturbation
    TracedStressType traced_stress_type = static_cast<TracedStressType>(this->GetValue(TRACED_STRESS_TYPE));
    if (rStressVariable == STRESS_ON_GP)
        StressCalculation::CalculateStressOnGP(*pGetPrimalElement(), traced_stress_type, stress_vector_undist, rCurrentProcessInfo);
    else
        StressCalculation::CalculateStressOnNode(*pGetPrimalElement(), traced_stress_type, stress_vector_undist, rCurrentProcessInfo);

    const SizeType stress_vector_size = stress_vector_undist.size();

    // Get perturbation size
    const double delta = this->GetPerturbationSize(rDesignVariable, rCurrentProcessInfo);

    if(rDesignVariable == SHAPE_SENSITIVITY)
    {
        const SizeType number_of_nodes = mpPrimalElement->GetGeometry().PointsNumber();
        const SizeType dimension = rCurrentProcessInfo.GetValue(DOMAIN_SIZE);


        rOutput.resize(dimension * number_of_nodes, stress_vector_size, false);

        IndexType index = 0;
        //TODO: look that this works also for parallel computing
        for(auto& node_i : mpPrimalElement->GetGeometry())
        {
            for(IndexType coord_dir_i = 0; coord_dir_i < dimension; ++coord_dir_i)
            {
                // perturb the design variable
                node_i.GetInitialPosition()[coord_dir_i] += delta;
                node_i[coord_dir_i] += delta;

                // Compute stress on GP after perturbation
                TracedStressType traced_stress_type = static_cast<TracedStressType>(this->GetValue(TRACED_STRESS_TYPE));
                if (rStressVariable == STRESS_ON_GP)
                    StressCalculation::CalculateStressOnGP(*pGetPrimalElement(), traced_stress_type, stress_vector_dist, rCurrentProcessInfo);
                else
                    StressCalculation::CalculateStressOnNode(*pGetPrimalElement(), traced_stress_type, stress_vector_dist, rCurrentProcessInfo);

                // Compute derivative of stress w.r.t. design variable with finite differences
                for(IndexType i = 0; i < stress_vector_size; ++i)
                    rOutput( (coord_dir_i + index*dimension), i) = (stress_vector_dist[i]-stress_vector_undist[i])/delta;

                // Reset pertubed vector
                stress_vector_dist = Vector(0);

                // unperturb the design variable
                node_i.GetInitialPosition()[coord_dir_i] -= delta;
                node_i[coord_dir_i] -= delta;

                }
            index++;
        }// end loop over element nodes
    }
    else
    {
        rOutput = ZeroMatrix(0, stress_vector_size);
    }

    KRATOS_CATCH("")
}

// private
template <class TPrimalElement>
double AdjointFiniteDifferencingBaseElement<TPrimalElement>::GetPerturbationSize(const Variable<double>& rDesignVariable, const ProcessInfo& rCurrentProcessInfo) const
{
    double delta = rCurrentProcessInfo[PERTURBATION_SIZE];
    if (rCurrentProcessInfo[ADAPT_PERTURBATION_SIZE]) {
            delta *= this->GetPerturbationSizeModificationFactor(rDesignVariable);
    }

    KRATOS_DEBUG_ERROR_IF_NOT(delta > 0) << "The perturbation size is not > 0!";
    return delta;
}

template <class TPrimalElement>
double AdjointFiniteDifferencingBaseElement<TPrimalElement>::GetPerturbationSize(const Variable<array_1d<double,3>>& rDesignVariable, const ProcessInfo& rCurrentProcessInfo) const
{
    double delta = rCurrentProcessInfo[PERTURBATION_SIZE];
    if (rCurrentProcessInfo[ADAPT_PERTURBATION_SIZE]) {
            delta *= this->GetPerturbationSizeModificationFactor(rDesignVariable);
    }

    KRATOS_DEBUG_ERROR_IF_NOT(delta > 0) << "The perturbation size is not > 0!";
    return delta;
}

template <class TPrimalElement>
double AdjointFiniteDifferencingBaseElement<TPrimalElement>::GetPerturbationSizeModificationFactor(const Variable<double>& rDesignVariable) const
{
    KRATOS_TRY;

    if ( mpPrimalElement->GetProperties().Has(rDesignVariable) )
    {
        const double variable_value = mpPrimalElement->GetProperties()[rDesignVariable];
        return variable_value;
    }
    else
        return 1.0;

    KRATOS_CATCH("")
}

template <class TPrimalElement>
double AdjointFiniteDifferencingBaseElement<TPrimalElement>::GetPerturbationSizeModificationFactor(const Variable<array_1d<double,3>>& rDesignVariable) const
{
    KRATOS_TRY;

    // For shape derivatives the size of the element (length, area, ...) is used as default perturbation size modification factor.
    // Later on this value is multiplied with a user defined factor. This product is then used as final perturbation size for computing
    // derivatives with finite differences.
    if(rDesignVariable == SHAPE_SENSITIVITY)
    {
        const double domain_size = mpPrimalElement->GetGeometry().DomainSize();
        KRATOS_DEBUG_ERROR_IF(domain_size <= 0.0)
            << "Pertubation size for shape derivatives of element" << this->Id() << "<= 0.0" << std::endl;
        return domain_size;
    }
    else
        return 1.0;

    KRATOS_CATCH("")
}

template <class TPrimalElement>
void AdjointFiniteDifferencingBaseElement<TPrimalElement>::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer,  Element );
    rSerializer.save("mpPrimalElement", mpPrimalElement);
    rSerializer.save("mHasRotationDofs", mHasRotationDofs);
}

template <class TPrimalElement>
void AdjointFiniteDifferencingBaseElement<TPrimalElement>::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Element );
    rSerializer.load("mpPrimalElement", mpPrimalElement);
    rSerializer.load("mHasRotationDofs", mHasRotationDofs);

}

template class AdjointFiniteDifferencingBaseElement<ShellThinElement3D3N>;
template class AdjointFiniteDifferencingBaseElement<CrBeamElementLinear3D2N>;
template class AdjointFiniteDifferencingBaseElement<TrussElement3D2N>;
template class AdjointFiniteDifferencingBaseElement<TrussElementLinear3D2N>;
template class AdjointFiniteDifferencingBaseElement<SmallDisplacement>;
template class AdjointFiniteDifferencingBaseElement<SpringDamperElement3D2N>;

} // namespace Kratos


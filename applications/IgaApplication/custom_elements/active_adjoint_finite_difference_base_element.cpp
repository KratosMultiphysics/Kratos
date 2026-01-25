// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Armin Geiser, https://github.com/armingeiser
//                   Leonhard Rieder (for active adjoint element implementation)

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "active_adjoint_finite_difference_base_element.h"

#include "custom_response_functions/response_utilities/stress_response_definitions.h"
#include "custom_response_functions/response_utilities/finite_difference_utility.h"
#include "includes/checks.h"

#include "custom_elements/active_adjoint_finite_difference_base_element.h"
#include "custom_elements/active_shell_3p_element.h"

#include "structural_mechanics_application_variables.h"


namespace Kratos
{

namespace AdjointFiniteDifferenceBaseElementHelperUtils
{

template <class TData>
void CalculateOnIntegrationPoints(
    Element& rPrimalElement,
    const Element& rAdjointElement,
    const Variable<TData>& rVariable,
    std::vector<TData>& rValues,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    if (rAdjointElement.Has(rVariable)) {
        // Get result value for output
        const auto& output_value = rAdjointElement.GetValue(rVariable);

        // Resize Output
        const std::size_t gauss_points_number = rAdjointElement.GetGeometry().IntegrationPointsNumber(rAdjointElement.GetIntegrationMethod());
        if (rValues.size() != gauss_points_number) {
            rValues.resize(gauss_points_number);
        }

        // Write scalar result value on all Gauss-Points
        for (IndexType i = 0; i < gauss_points_number; ++i) {
            rValues[i] = output_value;
        }
    }
    else {
        rPrimalElement.CalculateOnIntegrationPoints(rVariable, rValues, rCurrentProcessInfo);
    }

    KRATOS_CATCH("");
}
} // namespace AdjointFiniteDifferenceBaseElementHelperUtils

template <class TPrimalElement>
void ActiveAdjointFiniteDifferencingBaseElement<TPrimalElement>::EquationIdVector(EquationIdVectorType& rResult,
    const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY
    const SizeType number_of_control_points = GetGeometry().size();

    if (rResult.size() != 3 * number_of_control_points + 6)
        rResult.resize(3 * number_of_control_points + 6, false);

    const IndexType pos = this->GetGeometry()[0].GetDofPosition(ADJOINT_DISPLACEMENT_X);

    for (IndexType i = 0; i < number_of_control_points; ++i) {
        const IndexType index = i * 3;
        rResult[index]     = GetGeometry()[i].GetDof(ADJOINT_DISPLACEMENT_X, pos).EquationId();
        rResult[index + 1] = GetGeometry()[i].GetDof(ADJOINT_DISPLACEMENT_Y, pos + 1).EquationId();
        rResult[index + 2] = GetGeometry()[i].GetDof(ADJOINT_DISPLACEMENT_Z, pos + 2).EquationId();
    }

    // Global actuation-dofs from the parent-geometry actuation node
    const NodeType& r_global_node = GetGeometry().GetGeometryParent(0).GetValue(ACTIVE_SHELL_NODE_GP)[0];
    IndexType offset = 3 * number_of_control_points;
    rResult[offset + 0] = r_global_node.GetDof(ADJOINT_ACTIVE_SHELL_ALPHA).EquationId();
    rResult[offset + 1] = r_global_node.GetDof(ADJOINT_ACTIVE_SHELL_BETA).EquationId();
    rResult[offset + 2] = r_global_node.GetDof(ADJOINT_ACTIVE_SHELL_GAMMA).EquationId();
    rResult[offset + 3] = r_global_node.GetDof(ADJOINT_ACTIVE_SHELL_KAPPA_1).EquationId();
    rResult[offset + 4] = r_global_node.GetDof(ADJOINT_ACTIVE_SHELL_KAPPA_2).EquationId();
    rResult[offset + 5] = r_global_node.GetDof(ADJOINT_ACTIVE_SHELL_KAPPA_12).EquationId();

    KRATOS_CATCH("")
}

template <class TPrimalElement>
void ActiveAdjointFiniteDifferencingBaseElement<TPrimalElement>::GetDofList(DofsVectorType& rElementalDofList,
    const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY
    const SizeType number_of_control_points = GetGeometry().size();

    rElementalDofList.resize(0);
    rElementalDofList.reserve(3 * number_of_control_points + 6);

    // add adjoint displacements, but no rotations (-> Kirchhoff-Love shell theory)
    for (IndexType i = 0; i < number_of_control_points; ++i) {
        rElementalDofList.push_back(GetGeometry()[i].pGetDof(ADJOINT_DISPLACEMENT_X));
        rElementalDofList.push_back(GetGeometry()[i].pGetDof(ADJOINT_DISPLACEMENT_Y));
        rElementalDofList.push_back(GetGeometry()[i].pGetDof(ADJOINT_DISPLACEMENT_Z));
    }

    // add global actuation-dofs from parent-geometry actuation node:
    const NodeType& r_global_node = GetGeometry().GetGeometryParent(0).GetValue(ACTIVE_SHELL_NODE_GP)[0];
    rElementalDofList.push_back(r_global_node.pGetDof(ADJOINT_ACTIVE_SHELL_ALPHA));
    rElementalDofList.push_back(r_global_node.pGetDof(ADJOINT_ACTIVE_SHELL_BETA));
    rElementalDofList.push_back(r_global_node.pGetDof(ADJOINT_ACTIVE_SHELL_GAMMA));
    rElementalDofList.push_back(r_global_node.pGetDof(ADJOINT_ACTIVE_SHELL_KAPPA_1));
    rElementalDofList.push_back(r_global_node.pGetDof(ADJOINT_ACTIVE_SHELL_KAPPA_2));
    rElementalDofList.push_back(r_global_node.pGetDof(ADJOINT_ACTIVE_SHELL_KAPPA_12));

    KRATOS_CATCH("")
}

template <class TPrimalElement>
void ActiveAdjointFiniteDifferencingBaseElement<TPrimalElement>::GetValuesVector(Vector& rValues, int Step) const
{
    KRATOS_TRY

    const SizeType number_of_control_points = GetGeometry().size();

    rValues.resize(3 * number_of_control_points + 6);

    // set values to adjoint displacements, but no rotations (-> Kirchhoff-Love shell theory)
    IndexType local_index = 0;
    for (IndexType i = 0; i < number_of_control_points; ++i) {
        rValues[local_index++] = GetGeometry()[i].FastGetSolutionStepValue(ADJOINT_DISPLACEMENT_X);
        rValues[local_index++] = GetGeometry()[i].FastGetSolutionStepValue(ADJOINT_DISPLACEMENT_Y);
        rValues[local_index++] = GetGeometry()[i].FastGetSolutionStepValue(ADJOINT_DISPLACEMENT_Z);
    }

    //set values for adjoint actuation dofs
    const NodeType& r_global_node = GetGeometry().GetGeometryParent(0).GetValue(ACTIVE_SHELL_NODE_GP)[0];
    rValues[local_index++] = r_global_node.FastGetSolutionStepValue(ADJOINT_ACTIVE_SHELL_ALPHA);
    rValues[local_index++] = r_global_node.FastGetSolutionStepValue(ADJOINT_ACTIVE_SHELL_BETA);
    rValues[local_index++] = r_global_node.FastGetSolutionStepValue(ADJOINT_ACTIVE_SHELL_GAMMA);
    rValues[local_index++] = r_global_node.FastGetSolutionStepValue(ADJOINT_ACTIVE_SHELL_KAPPA_1);
    rValues[local_index++] = r_global_node.FastGetSolutionStepValue(ADJOINT_ACTIVE_SHELL_KAPPA_2);
    rValues[local_index++] = r_global_node.FastGetSolutionStepValue(ADJOINT_ACTIVE_SHELL_KAPPA_12);

    KRATOS_CATCH("")
}

template <class TPrimalElement>
void ActiveAdjointFiniteDifferencingBaseElement<TPrimalElement>::Calculate(const Variable<Matrix >& rVariable, Matrix& rOutput,
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
        KRATOS_WARNING("ActiveAdjointFiniteDifferencingBaseElement") << "Calculate function called for unknown variable: " << rVariable << std::endl;
        rOutput.clear();
    }

    KRATOS_CATCH("")
}

template <class TPrimalElement>
void ActiveAdjointFiniteDifferencingBaseElement<TPrimalElement>::CalculateOnIntegrationPoints(
    const Variable<bool>& rVariable,
    std::vector<bool>& rOutput,
    const ProcessInfo& rCurrentProcessInfo)
{
    AdjointFiniteDifferenceBaseElementHelperUtils::CalculateOnIntegrationPoints(*mpPrimalElement, *this, rVariable, rOutput, rCurrentProcessInfo);
}

template <class TPrimalElement>
void ActiveAdjointFiniteDifferencingBaseElement<TPrimalElement>::CalculateOnIntegrationPoints(
    const Variable<double>& rVariable,
    std::vector<double>& rOutput,
    const ProcessInfo& rCurrentProcessInfo)
{
    AdjointFiniteDifferenceBaseElementHelperUtils::CalculateOnIntegrationPoints(*mpPrimalElement, *this, rVariable, rOutput, rCurrentProcessInfo);
}

template <class TPrimalElement>
void ActiveAdjointFiniteDifferencingBaseElement<TPrimalElement>::CalculateOnIntegrationPoints(
    const Variable<array_1d<double, 3>>& rVariable,
    std::vector<array_1d<double, 3>>& rOutput,
    const ProcessInfo& rCurrentProcessInfo)
{
    AdjointFiniteDifferenceBaseElementHelperUtils::CalculateOnIntegrationPoints(*mpPrimalElement, *this, rVariable, rOutput, rCurrentProcessInfo);
}

template <class TPrimalElement>
void ActiveAdjointFiniteDifferencingBaseElement<TPrimalElement>::CalculateOnIntegrationPoints(
    const Variable<array_1d<double, 4>>& rVariable,
    std::vector<array_1d<double, 4>>& rOutput,
    const ProcessInfo& rCurrentProcessInfo)
{
    AdjointFiniteDifferenceBaseElementHelperUtils::CalculateOnIntegrationPoints(*mpPrimalElement, *this, rVariable, rOutput, rCurrentProcessInfo);
}

template <class TPrimalElement>
void ActiveAdjointFiniteDifferencingBaseElement<TPrimalElement>::CalculateOnIntegrationPoints(
    const Variable<array_1d<double, 6>>& rVariable,
    std::vector<array_1d<double, 6>>& rOutput,
    const ProcessInfo& rCurrentProcessInfo)
{
    AdjointFiniteDifferenceBaseElementHelperUtils::CalculateOnIntegrationPoints(*mpPrimalElement, *this, rVariable, rOutput, rCurrentProcessInfo);
}

template <class TPrimalElement>
void ActiveAdjointFiniteDifferencingBaseElement<TPrimalElement>::CalculateOnIntegrationPoints(
    const Variable<array_1d<double, 9>>& rVariable,
    std::vector<array_1d<double, 9>>& rOutput,
    const ProcessInfo& rCurrentProcessInfo)
{
    AdjointFiniteDifferenceBaseElementHelperUtils::CalculateOnIntegrationPoints(*mpPrimalElement, *this, rVariable, rOutput, rCurrentProcessInfo);
}

template <class TPrimalElement>
void ActiveAdjointFiniteDifferencingBaseElement<TPrimalElement>::CalculateOnIntegrationPoints(
    const Variable<Vector>& rVariable,
    std::vector<Vector>& rOutput,
    const ProcessInfo& rCurrentProcessInfo)
{
    AdjointFiniteDifferenceBaseElementHelperUtils::CalculateOnIntegrationPoints(*mpPrimalElement, *this, rVariable, rOutput, rCurrentProcessInfo);
}

template <class TPrimalElement>
void ActiveAdjointFiniteDifferencingBaseElement<TPrimalElement>::CalculateOnIntegrationPoints(
    const Variable<Matrix>& rVariable,
    std::vector<Matrix>& rOutput,
    const ProcessInfo& rCurrentProcessInfo)
{
    AdjointFiniteDifferenceBaseElementHelperUtils::CalculateOnIntegrationPoints(*mpPrimalElement, *this, rVariable, rOutput, rCurrentProcessInfo);
}

template <class TPrimalElement>
int ActiveAdjointFiniteDifferencingBaseElement<TPrimalElement>::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    int return_value = Element::Check(rCurrentProcessInfo);

    //Temporarily commented out Check function to facilitate debugging. Re-enable once debugging is complete.
    // KRATOS_ERROR_IF_NOT(mpPrimalElement) << "Primal element pointer is nullptr!" << std::endl;

    // const GeometryType& r_geom = this->GetGeometry();

    // // TODO generic way of doing these checks without checking the dofs..

    // // Check dofs
    // for (IndexType i = 0; i < r_geom.size(); ++i) {
    //     const auto& r_node = r_geom[i];

    //     KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DISPLACEMENT, r_node);
    //     KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(ADJOINT_DISPLACEMENT, r_node);

    //     KRATOS_CHECK_DOF_IN_NODE(ADJOINT_DISPLACEMENT_X, r_node);
    //     KRATOS_CHECK_DOF_IN_NODE(ADJOINT_DISPLACEMENT_Y, r_node);
    //     KRATOS_CHECK_DOF_IN_NODE(ADJOINT_DISPLACEMENT_Z, r_node);

    // }

    return return_value;

    KRATOS_CATCH("")
}

// Sensitivity functions

template <class TPrimalElement>
void ActiveAdjointFiniteDifferencingBaseElement<TPrimalElement>::CalculateSensitivityMatrix(const Variable<double>& rDesignVariable, Matrix& rOutput,
                                            const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    // Get perturbation size
    const double delta = this->GetPerturbationSize(rDesignVariable, rCurrentProcessInfo);
    const SizeType number_of_nodes = mpPrimalElement->GetGeometry().PointsNumber();
    const SizeType dimension = rCurrentProcessInfo.GetValue(DOMAIN_SIZE);
    const SizeType num_dofs_per_node = mHasRotationDofs ?  2 * dimension : dimension;
    const SizeType local_size = number_of_nodes * num_dofs_per_node;

    if( rDesignVariable == TEMPERATURE )
    {
        rOutput.resize( number_of_nodes, local_size, false);

        Vector RHS;
        Vector derived_RHS;
        
        pGetPrimalElement()->CalculateRightHandSide(RHS, rCurrentProcessInfo);

        for(IndexType i_node = 0; i_node < mpPrimalElement->GetGeometry().PointsNumber(); ++i_node)
        {
            // Get pseudo-load contribution from utility
            FiniteDifferenceUtility::CalculateRightHandSideDerivative(*pGetPrimalElement(), RHS, rDesignVariable,
                                                                      mpPrimalElement->GetGeometry()[i_node], delta, 
                                                                      derived_RHS, rCurrentProcessInfo);

            KRATOS_ERROR_IF_NOT(derived_RHS.size() == local_size) << "Size of the pseudo-load does not fit! [ derived_RHS.size() = " << derived_RHS.size() << ", local_size = " << local_size << " ]." << std::endl;

            for(IndexType i = 0; i < derived_RHS.size(); ++i)
                rOutput(i_node, i) = derived_RHS[i];
        }
    }

    // FD derivative of primal RHS wrt global actuation variable on ACTIVE_SHELL_NODE_GP[0] (perturb +delta, recompute RHS, difference/ delta, restore).
    else if (
        rDesignVariable == ACTIVE_SHELL_ALPHA || 
        rDesignVariable == ACTIVE_SHELL_BETA || 
        rDesignVariable == ACTIVE_SHELL_GAMMA ||
        rDesignVariable == ACTIVE_SHELL_KAPPA_1 ||
        rDesignVariable == ACTIVE_SHELL_KAPPA_2 ||
        rDesignVariable == ACTIVE_SHELL_KAPPA_12)
    {
        Vector ref_RHS, perturbed_RHS;
        this->pGetPrimalElement()->CalculateRightHandSide(ref_RHS, rCurrentProcessInfo);

        NodeType& r_global_node = GetGeometry().GetGeometryParent(0).GetValue(ACTIVE_SHELL_NODE_GP)[0];

        r_global_node.FastGetSolutionStepValue(rDesignVariable) += delta;
        this->pGetPrimalElement()->CalculateRightHandSide(perturbed_RHS, rCurrentProcessInfo);
        r_global_node.FastGetSolutionStepValue(rDesignVariable) -= delta;

        noalias(perturbed_RHS) = (perturbed_RHS - ref_RHS) / delta;

        rOutput.resize(1, local_size, false);
        std::copy(perturbed_RHS.begin(), perturbed_RHS.end(), rOutput.data().begin());
    }
    else
    {
        Vector RHS;
        this->pGetPrimalElement()->CalculateRightHandSide(RHS, rCurrentProcessInfo);

        // Get pseudo-load from utility
        FiniteDifferenceUtility::CalculateRightHandSideDerivative(*pGetPrimalElement(), RHS, rDesignVariable, delta, rOutput, rCurrentProcessInfo);
    }

    if (rOutput.size1() == 0 || rOutput.size2() == 0)
    {
        rOutput = ZeroMatrix(0, local_size);
    }

    KRATOS_CATCH("")
}

template <class TPrimalElement>
void ActiveAdjointFiniteDifferencingBaseElement<TPrimalElement>::CalculateSensitivityMatrix(const Variable<array_1d<double,3>>& rDesignVariable, Matrix& rOutput,
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
void ActiveAdjointFiniteDifferencingBaseElement<TPrimalElement>::CalculateStressDisplacementDerivative(const Variable<Vector>& rStressVariable,
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
void ActiveAdjointFiniteDifferencingBaseElement<TPrimalElement>::CalculateStressDesignVariableDerivative(const Variable<double>& rDesignVariable,
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
void ActiveAdjointFiniteDifferencingBaseElement<TPrimalElement>::CalculateStressDesignVariableDerivative(const Variable<array_1d<double,3>>& rDesignVariable,
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
double ActiveAdjointFiniteDifferencingBaseElement<TPrimalElement>::GetPerturbationSize(const Variable<double>& rDesignVariable, const ProcessInfo& rCurrentProcessInfo) const
{
    double delta = rCurrentProcessInfo[PERTURBATION_SIZE];
    if (rCurrentProcessInfo[ADAPT_PERTURBATION_SIZE]) {
            delta *= this->GetPerturbationSizeModificationFactor(rDesignVariable);
    }

    KRATOS_DEBUG_ERROR_IF_NOT(delta > 0) << "The perturbation size is not > 0!";
    return delta;
}

template <class TPrimalElement>
double ActiveAdjointFiniteDifferencingBaseElement<TPrimalElement>::GetPerturbationSize(const Variable<array_1d<double,3>>& rDesignVariable, const ProcessInfo& rCurrentProcessInfo) const
{
    double delta = rCurrentProcessInfo[PERTURBATION_SIZE];
    if (rCurrentProcessInfo[ADAPT_PERTURBATION_SIZE]) {
            delta *= this->GetPerturbationSizeModificationFactor(rDesignVariable);
    }

    KRATOS_DEBUG_ERROR_IF_NOT(delta > 0) << "The perturbation size is not > 0!";
    return delta;
}

template <class TPrimalElement>
double ActiveAdjointFiniteDifferencingBaseElement<TPrimalElement>::GetPerturbationSizeModificationFactor(const Variable<double>& rDesignVariable) const
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
double ActiveAdjointFiniteDifferencingBaseElement<TPrimalElement>::GetPerturbationSizeModificationFactor(const Variable<array_1d<double,3>>& rDesignVariable) const
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
void ActiveAdjointFiniteDifferencingBaseElement<TPrimalElement>::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer,  Element );
    rSerializer.save("mpPrimalElement", mpPrimalElement);
    rSerializer.save("mHasRotationDofs", mHasRotationDofs);
}

template <class TPrimalElement>
void ActiveAdjointFiniteDifferencingBaseElement<TPrimalElement>::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Element );
    rSerializer.load("mpPrimalElement", mpPrimalElement);
    rSerializer.load("mHasRotationDofs", mHasRotationDofs);

}

// template instantiations - other templates where removed (see adjoint-base element)
template class ActiveAdjointFiniteDifferencingBaseElement<ActiveShell3pElement>;


} // namespace Kratos


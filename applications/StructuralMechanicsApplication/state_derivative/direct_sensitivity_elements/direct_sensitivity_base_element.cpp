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
#include "direct_sensitivity_base_element.h"
#include "structural_mechanics_application_variables.h"
#include "custom_response_functions/response_utilities/element_finite_difference_utility.h"
#include "includes/checks.h"



namespace Kratos
{

DirectSensitivityBaseElement::DirectSensitivityBaseElement(Element::Pointer pPrimalElement, DirectSensitivityVariable& rDirectSensitivityVariable)
                    : Element(pPrimalElement->Id(), pPrimalElement->pGetGeometry(), pPrimalElement->pGetProperties())
                    , mpPrimalElement(pPrimalElement)
                    , mrDirectSensitivityVariable(rDirectSensitivityVariable)
{
}

DirectSensitivityBaseElement::DirectSensitivityBaseElement(Element::Pointer pPrimalElement, bool HasRotationDofs, DirectSensitivityVariable& rDirectSensitivityVariable)
                    : Element(pPrimalElement->Id(), pPrimalElement->pGetGeometry(), pPrimalElement->pGetProperties())
                    , mpPrimalElement(pPrimalElement)
                    , mHasRotationDofs(HasRotationDofs)
                    , mrDirectSensitivityVariable(rDirectSensitivityVariable)
{
}

DirectSensitivityBaseElement::~DirectSensitivityBaseElement() {}

void DirectSensitivityBaseElement::EquationIdVector(EquationIdVectorType& rResult,
    ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    GeometryType& geom = this->GetGeometry();

    const SizeType number_of_nodes = geom.PointsNumber();
    const SizeType dimension = geom.WorkingSpaceDimension();
    const SizeType num_dofs_per_node = (mHasRotationDofs) ?  2 * dimension : dimension;
    const SizeType num_dofs = number_of_nodes * num_dofs_per_node;

    if(rResult.size() != num_dofs)
        rResult.resize(num_dofs, false);

    for(IndexType i = 0; i < geom.size(); ++i)
    {
        const IndexType index = i * num_dofs_per_node;
        NodeType& iNode = geom[i];

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

void DirectSensitivityBaseElement::GetDofList(DofsVectorType& rElementalDofList,
    ProcessInfo& rCurrentProcessInfo)
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
        rElementalDofList[index    ] = GetGeometry()[i].pGetDof(ADJOINT_DISPLACEMENT_X);
        rElementalDofList[index + 1] = GetGeometry()[i].pGetDof(ADJOINT_DISPLACEMENT_Y);
        rElementalDofList[index + 2] = GetGeometry()[i].pGetDof(ADJOINT_DISPLACEMENT_Z);

        if(mHasRotationDofs)
        {
            rElementalDofList[index + 3] = GetGeometry()[i].pGetDof(ADJOINT_ROTATION_X);
            rElementalDofList[index + 4] = GetGeometry()[i].pGetDof(ADJOINT_ROTATION_Y);
            rElementalDofList[index + 5] = GetGeometry()[i].pGetDof(ADJOINT_ROTATION_Z);
        }
    }
    KRATOS_CATCH("")
}

void DirectSensitivityBaseElement::GetValuesVector(Vector& rValues, int Step)
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

void DirectSensitivityBaseElement::Calculate(const Variable<Matrix >& rVariable, Matrix& rOutput,
                                                const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    std::vector<TracedStressType> TracedStressTypeVector;
    this->GetListOfTracedStressTypes(TracedStressTypeVector);
    
    if(rVariable == STRESS_DISP_DERIV_ON_GP)
        this->AssembleStressDisplacementDerivativeMatrix(STRESS_ON_GP, TracedStressTypeVector, rOutput, rCurrentProcessInfo);
    else if (rVariable == STRESS_DISP_DERIV_ON_NODE)    
        this->AssembleStressDisplacementDerivativeMatrix(STRESS_ON_NODE, TracedStressTypeVector, rOutput, rCurrentProcessInfo);     
    else if(rVariable == STRESS_DESIGN_DERIVATIVE_ON_GP)
        this->CalculateStressDesignVariableDerivative(STRESS_ON_GP, TracedStressTypeVector, rOutput, rCurrentProcessInfo);
    
    else if(rVariable == STRESS_DESIGN_DERIVATIVE_ON_NODE)
        this->CalculateStressDesignVariableDerivative(STRESS_ON_NODE, TracedStressTypeVector, rOutput, rCurrentProcessInfo);
    else
    {
        KRATOS_WARNING("DirectSensitivityBaseElement") << "Calculate function called for unknown variable: " << rVariable << std::endl;
        rOutput.clear();
    }

    KRATOS_CATCH("")
}


// Sensitivity functions

void DirectSensitivityBaseElement::AssembleStressDisplacementDerivativeMatrix(const Variable<Vector>& rStressVariable,  
                                                std::vector<TracedStressType>& rTracedStressTypeVector,
                                                Matrix& rOutput,
                                                const ProcessInfo& rCurrentProcessInfo)
{
    const SizeType num_stress = rTracedStressTypeVector.size();
    const SizeType num_gauss_points = GetGeometry().IntegrationPointsNumber(this->GetIntegrationMethod()); 
    const SizeType num_nodes = mpPrimalElement->GetGeometry().PointsNumber();
    const SizeType dimension = mpPrimalElement->GetGeometry().WorkingSpaceDimension();
    const SizeType num_dofs_per_node = (mHasRotationDofs) ?  2 * dimension : dimension;
    const SizeType num_dofs = num_nodes * num_dofs_per_node;
    const SizeType num_derivatives = num_dofs * num_stress;
    SizeType num_points = 0;

    if(rStressVariable == STRESS_ON_GP)
        num_points = num_gauss_points;
    else
        num_points = num_nodes;

    if (rOutput.size1() != num_derivatives || rOutput.size2() != num_points)
        rOutput.resize(num_derivatives, num_points);
        
    Matrix DerivativeMatrix;
    
    for (IndexType i = 0; i < num_stress; i++)
    {   
        IndexType index = i * num_dofs;

        TracedStressType traced_stress_type = rTracedStressTypeVector[i];
        
        this->CalculateStressDisplacementDerivative(rStressVariable, traced_stress_type, DerivativeMatrix, rCurrentProcessInfo);
            
        for (IndexType j = 0; j < num_dofs; j++)
            for (IndexType k = 0; k < num_points; j++)
                rOutput(index + j, k) = DerivativeMatrix(j, k); 
    }
}



void DirectSensitivityBaseElement::CalculateStressDisplacementDerivative(const Variable<Vector>& rStressVariable, 
                                            TracedStressType& rTracedStressType,
                                            Matrix& rOutput, 
                                            const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    const SizeType num_nodes = mpPrimalElement->GetGeometry().PointsNumber();
    const SizeType dimension = mpPrimalElement->GetGeometry().WorkingSpaceDimension();
    const SizeType num_dofs_per_node = (mHasRotationDofs) ?  2 * dimension : dimension;
    const SizeType num_dofs = num_nodes * num_dofs_per_node;
    Vector initial_state_variables;
    Vector stress_derivatives_vector;

    // TODO first calculation only to get the size of the stress vector
    if (rStressVariable == STRESS_ON_GP)
        StressCalculation::CalculateStressOnGP(*pGetPrimalElement(), rTracedStressType, stress_derivatives_vector, rCurrentProcessInfo);
    else
        StressCalculation::CalculateStressOnNode(*pGetPrimalElement(), rTracedStressType, stress_derivatives_vector, rCurrentProcessInfo);
    rOutput.resize(num_dofs, stress_derivatives_vector.size() );
    rOutput.clear();
    initial_state_variables.resize(num_dofs);

    // Build vector of variables containing the DOF-variables of the primal problem
    std::vector<VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>>> primal_solution_variable_list;
    primal_solution_variable_list.reserve(num_dofs_per_node);
    primal_solution_variable_list.push_back(DISPLACEMENT_X);
    primal_solution_variable_list.push_back(DISPLACEMENT_Y);
    primal_solution_variable_list.push_back(DISPLACEMENT_Z);

    if(mHasRotationDofs)
    {
        primal_solution_variable_list.push_back(ROTATION_X);
        primal_solution_variable_list.push_back(ROTATION_Y);
        primal_solution_variable_list.push_back(ROTATION_Z);
    }

    // TODO Find a better way of doing this check
    KRATOS_ERROR_IF(rCurrentProcessInfo.Has(NL_ITERATION_NUMBER))
        << "This stress displacement derivative computation is only usable for linear cases!" << std::endl;

    for (IndexType i = 0; i < num_nodes; ++i)
    {
        const IndexType index = i * num_dofs_per_node;
        for(IndexType j = 0; j < primal_solution_variable_list.size(); ++j)
        {
            initial_state_variables[index + j] = mpPrimalElement->GetGeometry()[i].FastGetSolutionStepValue(primal_solution_variable_list[j]);
            mpPrimalElement->GetGeometry()[i].FastGetSolutionStepValue(primal_solution_variable_list[j]) = 0.0;
        }
    }
    for (IndexType i = 0; i < num_nodes; ++i)
    {
        const IndexType index = i * num_dofs_per_node;
        for(IndexType j = 0; j < primal_solution_variable_list.size(); ++j)
        {
            mpPrimalElement->GetGeometry()[i].FastGetSolutionStepValue(primal_solution_variable_list[j]) = 1.0;

            if (rStressVariable == STRESS_ON_GP)
                StressCalculation::CalculateStressOnGP(*pGetPrimalElement(), rTracedStressType, stress_derivatives_vector, rCurrentProcessInfo);
            else
                StressCalculation::CalculateStressOnNode(*pGetPrimalElement(), rTracedStressType, stress_derivatives_vector, rCurrentProcessInfo);

            for(IndexType k = 0; k < stress_derivatives_vector.size(); ++k)
                rOutput(index+j, k) = stress_derivatives_vector[k];

            stress_derivatives_vector.clear();

            mpPrimalElement->GetGeometry()[i].FastGetSolutionStepValue(primal_solution_variable_list[j]) = 0.0;
        }
    }
    for (IndexType i = 0; i < num_nodes; ++i)
    {
        const IndexType index = i * num_dofs_per_node;
        for(IndexType j = 0; j < primal_solution_variable_list.size(); ++j)
            mpPrimalElement->GetGeometry()[i].FastGetSolutionStepValue(primal_solution_variable_list[j]) = initial_state_variables[index + j];
    }

    KRATOS_CATCH("")
}


void DirectSensitivityBaseElement::CalculateStressDesignVariableDerivative( const Variable<Vector>& rStressVariable,
                                                std::vector<TracedStressType>& rTracedStressTypeVector, 
                                                Matrix& rOutput,
                                                const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    const SizeType num_stress = rTracedStressTypeVector.size(); 
    const SizeType num_gauss_points = GetGeometry().IntegrationPointsNumber(this->GetIntegrationMethod()); 
    const SizeType num_nodes = mpPrimalElement->GetGeometry().PointsNumber();  

    if (rStressVariable == STRESS_ON_GP)
        rOutput.resize(num_stress, num_gauss_points);
    else
        rOutput.resize(num_stress, num_nodes); 
    
    // Define working variables
    Vector stress_vector_undist;
    Vector stress_vector_dist;

    rOutput.clear();

    // Sizing of rOutput
    if( mpPrimalElement->Id() == mrDirectSensitivityVariable.GetTracedElementId() )
    { 
        for(IndexType i = 0; i < num_stress; ++i)
        {  
            TracedStressType traced_stress_type = rTracedStressTypeVector[i];

            // Compute stress before perturbation
            if (rStressVariable == STRESS_ON_GP)
                StressCalculation::CalculateStressOnGP(*pGetPrimalElement(), traced_stress_type, stress_vector_undist, rCurrentProcessInfo);
            else
                StressCalculation::CalculateStressOnNode(*pGetPrimalElement(), traced_stress_type, stress_vector_undist, rCurrentProcessInfo);

            // Get size of stress_vector
            const SizeType stress_vector_size = stress_vector_undist.size();            
            
            // Get perturbation size
            const double delta = mrDirectSensitivityVariable.GetPerturbationSize();
             
            // Perturb the properties of the element
            mrDirectSensitivityVariable.PerturbDesignVariable(*mpPrimalElement);
        
            // Compute stress after perturbation        
            if (rStressVariable == STRESS_ON_GP)
                StressCalculation::CalculateStressOnGP(*pGetPrimalElement(), traced_stress_type, stress_vector_dist, rCurrentProcessInfo);
            else
                StressCalculation::CalculateStressOnNode(*pGetPrimalElement(), traced_stress_type, stress_vector_dist, rCurrentProcessInfo);

            // Compute derivative of stress w.r.t. design variable with finite differences        
            for(IndexType j = 0; j < stress_vector_size; ++j)
                rOutput(i, j) = (stress_vector_dist[j] - stress_vector_undist[j]) / delta;

            // Give element original properties back
            mrDirectSensitivityVariable.UnperturbDesignVariable(*mpPrimalElement);
        }
    }
    else
        rOutput.resize(0, 0, false);

    KRATOS_CATCH("")
}


void DirectSensitivityBaseElement::CalculateOnIntegrationPoints(const Variable<double>& rVariable,
                    std::vector<double>& rValues,
                    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    if(this->Has(rVariable))
    {
        // Get result value for output
        const double& output_value = this->GetValue(rVariable);

        // Resize Output
        const SizeType  write_points_number = GetGeometry()
            .IntegrationPointsNumber(this->GetIntegrationMethod());
        if (rValues.size() != write_points_number)
            rValues.resize(write_points_number);

        // Write scalar result value on all Gauss-Points
        for(IndexType i = 0; i < write_points_number; ++i)
            rValues[i] = output_value;
    }
    else
        KRATOS_ERROR << "Unsupported output variable." << std::endl;

    KRATOS_CATCH("")
}


int DirectSensitivityBaseElement::Check(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    int return_value = Element::Check(rCurrentProcessInfo);

    KRATOS_ERROR_IF_NOT(mpPrimalElement) << "Primal element pointer is nullptr!" << std::endl;

    GeometryType& r_geom = GetGeometry();

    // verify that the variables are correctly initialized
    KRATOS_CHECK_VARIABLE_KEY(DISPLACEMENT);
    KRATOS_CHECK_VARIABLE_KEY(VELOCITY);
    KRATOS_CHECK_VARIABLE_KEY(ACCELERATION);
    KRATOS_CHECK_VARIABLE_KEY(DENSITY);
    KRATOS_CHECK_VARIABLE_KEY(CONSTITUTIVE_LAW);
    KRATOS_CHECK_VARIABLE_KEY(ADJOINT_DISPLACEMENT);

    if(mHasRotationDofs)
    {
        KRATOS_CHECK_VARIABLE_KEY(ROTATION);
        KRATOS_CHECK_VARIABLE_KEY(ADJOINT_ROTATION);
    }

    // TODO generic way of doing these checks without checking the dofs..

    // Check dofs
    for (IndexType i = 0; i < r_geom.size(); ++i)
    {
        auto& r_node = r_geom[i];

        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DISPLACEMENT, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(ADJOINT_DISPLACEMENT, r_node);

        KRATOS_CHECK_DOF_IN_NODE(ADJOINT_DISPLACEMENT_X, r_node);
        KRATOS_CHECK_DOF_IN_NODE(ADJOINT_DISPLACEMENT_Y, r_node);
        KRATOS_CHECK_DOF_IN_NODE(ADJOINT_DISPLACEMENT_Z, r_node);

        if(mHasRotationDofs)
        {
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

// Sensitivity Matrix

void DirectSensitivityBaseElement::CalculateSensitivityMatrix(const Variable<double>& rDesignVariable, Matrix& rOutput,
                                            const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    // Get perturbation size
    const double delta = this->GetPerturbationSize(rDesignVariable);
    ProcessInfo process_info = rCurrentProcessInfo;

    Vector RHS;
    pGetPrimalElement()->CalculateRightHandSide(RHS, process_info);

    // Get pseudo-load from utility
    ElementFiniteDifferenceUtility::CalculateRightHandSideDerivative(*pGetPrimalElement(), RHS, rDesignVariable, delta, rOutput, process_info);

    KRATOS_CATCH("")
}

void DirectSensitivityBaseElement::CalculateSensitivityMatrix(const Variable<array_1d<double,3>>& rDesignVariable, Matrix& rOutput,
                                            const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    const double delta = this->GetPerturbationSize(rDesignVariable);
    ProcessInfo process_info = rCurrentProcessInfo;

    if( rDesignVariable == SHAPE )
    {
        const SizeType number_of_nodes = mpPrimalElement->GetGeometry().PointsNumber();
        const SizeType dimension = rCurrentProcessInfo.GetValue(DOMAIN_SIZE);
        const SizeType num_dofs_per_node = (mHasRotationDofs) ?  2 * dimension : dimension;
        const SizeType local_size = number_of_nodes * num_dofs_per_node;
        const std::vector<ElementFiniteDifferenceUtility::array_1d_component_type> coord_directions = {SHAPE_X, SHAPE_Y, SHAPE_Z};
        Vector derived_RHS;

        if ( (rOutput.size1() != dimension * number_of_nodes) || (rOutput.size2() != local_size ) )
            rOutput.resize(dimension * number_of_nodes, local_size);

        IndexType index = 0;
        
        Vector RHS;
        pGetPrimalElement()->CalculateRightHandSide(RHS, process_info);
        for(auto& node_i : mpPrimalElement->GetGeometry())
        {
            for(IndexType coord_dir_i = 0; coord_dir_i < dimension; ++coord_dir_i)
            {
                // Get pseudo-load contribution from utility
                ElementFiniteDifferenceUtility::CalculateRightHandSideDerivative(*pGetPrimalElement(), RHS, coord_directions[coord_dir_i],
                                                                            node_i, delta, derived_RHS, process_info);

                KRATOS_ERROR_IF_NOT(derived_RHS.size() == local_size) << "Size of the pseudo-load does not fit!" << std::endl;

                for(IndexType i = 0; i < derived_RHS.size(); ++i)
                    rOutput( (coord_dir_i + index*dimension), i) = derived_RHS[i];
            }
            index++;
        }
    }

    KRATOS_CATCH("")
}

void DirectSensitivityBaseElement::GetListOfTracedStressTypes(std::vector<TracedStressType>& rTracedStressTypeVector)
    {
        KRATOS_TRY;
        
        if (mHasRotationDofs)
        {
            std::vector<std::string> traced_stress = {"FX", "FY", "FZ", "MX", "MY", "MZ"};
            
            SizeType num_stress = traced_stress.size();

            if (rTracedStressTypeVector.size() != num_stress)
                rTracedStressTypeVector.resize(num_stress);
        
            for( IndexType stress_id = 0; stress_id < num_stress; stress_id++ )
                rTracedStressTypeVector[stress_id] = StressResponseDefinitions::ConvertStringToTracedStressType(traced_stress[stress_id]);
        }

        else
        {
            std::vector<std::string> traced_stress = {"FX", "FY", "FZ"};

            SizeType num_stress = traced_stress.size();     
        
            if (rTracedStressTypeVector.size() != num_stress)
                rTracedStressTypeVector.resize(num_stress);
        
            for( IndexType stress_id = 0; stress_id < num_stress; stress_id++ )
                rTracedStressTypeVector[stress_id] = StressResponseDefinitions::ConvertStringToTracedStressType(traced_stress[stress_id]);
        }

        KRATOS_CATCH("");
    }  


// private
double DirectSensitivityBaseElement::GetPerturbationSize(const Variable<double>& rDesignVariable)
{
    const double correction_factor = this->GetPerturbationSizeModificationFactor(rDesignVariable);
    const double delta = this->GetValue(PERTURBATION_SIZE) * correction_factor;
    KRATOS_DEBUG_ERROR_IF_NOT(delta > 0) << "The perturbation size is not > 0!";
    return delta;
}

double DirectSensitivityBaseElement::GetPerturbationSize(const Variable<array_1d<double,3>>& rDesignVariable)
{
    const double correction_factor = this->GetPerturbationSizeModificationFactor(rDesignVariable);
    const double delta = this->GetValue(PERTURBATION_SIZE) * correction_factor;
    KRATOS_DEBUG_ERROR_IF_NOT(delta > 0) << "The perturbation size is not > 0!";
    return delta;
}

double DirectSensitivityBaseElement::GetPerturbationSizeModificationFactor(const Variable<double>& rDesignVariable)
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

double DirectSensitivityBaseElement::GetPerturbationSizeModificationFactor(const Variable<array_1d<double,3>>& rDesignVariable)
{
    KRATOS_TRY;

    // For shape derivatives the size of the element (length, area, ...) is used as default perturbation size modification factor.
    // Later on this value is multiplied with a user defined factor. This product is then used as final perturbation size for computing
    // derivatives with finite differences.
    if(rDesignVariable == SHAPE)
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

void DirectSensitivityBaseElement::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer,  Element );
    rSerializer.save("mpPrimalElement", mpPrimalElement);
}

void DirectSensitivityBaseElement::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Element );
    rSerializer.load("mpPrimalElement", mpPrimalElement);

}

} // namespace Kratos


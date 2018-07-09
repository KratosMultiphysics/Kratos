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
#include "structural_mechanics_application_variables.h"
#include "custom_response_functions/response_utilities/response_data.h"
#include "includes/checks.h"


namespace Kratos
{

AdjointFiniteDifferencingBaseElement::AdjointFiniteDifferencingBaseElement(Element::Pointer pPrimalElement)
                    : Element(pPrimalElement->Id(), pPrimalElement->pGetGeometry(), pPrimalElement->pGetProperties()),
                        mpPrimalElement(pPrimalElement)
{
}

AdjointFiniteDifferencingBaseElement::~AdjointFiniteDifferencingBaseElement() {}

Element::Pointer AdjointFiniteDifferencingBaseElement::Create(Element::Pointer pPrimalElement) const
{
    return Kratos::make_shared<AdjointFiniteDifferencingBaseElement>(pPrimalElement);
}

void AdjointFiniteDifferencingBaseElement::EquationIdVector(EquationIdVectorType& rResult,
    ProcessInfo& rCurrentProcessInfo)
{
    const SizeType num_dofs = 6 * GetGeometry().PointsNumber(); // 6 dofs per node
    if(rResult.size() != num_dofs)
        rResult.resize(num_dofs, false);

    GeometryType & geom = this->GetGeometry();

    for(SizeType i = 0; i < geom.size(); i++)
    {
        int index = i * 6;
        NodeType & iNode = geom[i];

        rResult[index]     = iNode.GetDof(ADJOINT_DISPLACEMENT_X).EquationId();
        rResult[index + 1] = iNode.GetDof(ADJOINT_DISPLACEMENT_Y).EquationId();
        rResult[index + 2] = iNode.GetDof(ADJOINT_DISPLACEMENT_Z).EquationId();

        rResult[index + 3] = iNode.GetDof(ADJOINT_ROTATION_X).EquationId();
        rResult[index + 4] = iNode.GetDof(ADJOINT_ROTATION_Y).EquationId();
        rResult[index + 5] = iNode.GetDof(ADJOINT_ROTATION_Z).EquationId();
    }

}

void AdjointFiniteDifferencingBaseElement::GetDofList(DofsVectorType& rElementalDofList,
    ProcessInfo& rCurrentProcessInfo) {

    const SizeType num_dofs = 6 * GetGeometry().PointsNumber(); // 6 dofs per node
    rElementalDofList.resize(0);
    rElementalDofList.reserve(num_dofs);

    GeometryType & geom = this->GetGeometry();

    for (SizeType i = 0; i < geom.size(); i++)
    {
        NodeType & iNode = geom[i];

        rElementalDofList.push_back(iNode.pGetDof(ADJOINT_DISPLACEMENT_X));
        rElementalDofList.push_back(iNode.pGetDof(ADJOINT_DISPLACEMENT_Y));
        rElementalDofList.push_back(iNode.pGetDof(ADJOINT_DISPLACEMENT_Z));

        rElementalDofList.push_back(iNode.pGetDof(ADJOINT_ROTATION_X));
        rElementalDofList.push_back(iNode.pGetDof(ADJOINT_ROTATION_Y));
        rElementalDofList.push_back(iNode.pGetDof(ADJOINT_ROTATION_Z));
    }
}

void AdjointFiniteDifferencingBaseElement::GetValuesVector(Vector& rValues, int Step) {

    const SizeType num_dofs = 6 * GetGeometry().PointsNumber(); // 6 dofs per node
    if(rValues.size() != num_dofs)
        rValues.resize(num_dofs, false);

    const GeometryType & geom = this->GetGeometry();
    const int dimension = geom.WorkingSpaceDimension();

    for (SizeType i = 0; i < geom.size(); i++)
    {
        const NodeType & iNode = geom[i];
        const array_1d<double,3>& disp = iNode.FastGetSolutionStepValue(ADJOINT_DISPLACEMENT, Step);
        const array_1d<double,3>& rot = iNode.FastGetSolutionStepValue(ADJOINT_ROTATION, Step);

        const SizeType index = i * dimension * 2;
        rValues[index]     = disp[0];
        rValues[index + 1] = disp[1];
        rValues[index + 2] = disp[2];

        rValues[index + 3] = rot[0];
        rValues[index + 4] = rot[1];
        rValues[index + 5] = rot[2];
    }
}

void AdjointFiniteDifferencingBaseElement::Calculate(const Variable<Matrix >& rVariable, Matrix& rOutput,
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
        const std::string design_varible_name = this->GetValue( DESIGN_VARIABLE_NAME );

        if (KratosComponents<Variable<double>>::Has(design_varible_name) == true)
        {
            const Variable<double>& r_variable =
                KratosComponents<Variable<double>>::Get(design_varible_name);
            this->CalculateStressDesignVariableDerivative(r_variable, STRESS_ON_GP, rOutput, rCurrentProcessInfo);
        }
        else if (KratosComponents<Variable<array_1d<double, 3>>>::Has(design_varible_name) == true)
        {
            const Variable<array_1d<double, 3>>& r_variable =
                KratosComponents<Variable<array_1d<double, 3>>>::Get(design_varible_name);
            this->CalculateStressDesignVariableDerivative(r_variable, STRESS_ON_GP, rOutput, rCurrentProcessInfo);
        }
    }
    else if(rVariable == STRESS_DESIGN_DERIVATIVE_ON_NODE)
    {
        std::string design_varible_name = this->GetValue( DESIGN_VARIABLE_NAME );

        if (KratosComponents<Variable<double>>::Has(design_varible_name) == true)
        {
            const Variable<double>& r_variable =
                KratosComponents<Variable<double>>::Get(design_varible_name);
            this->CalculateStressDesignVariableDerivative(r_variable, STRESS_ON_NODE, rOutput, rCurrentProcessInfo);
        }
        else if (KratosComponents<Variable<array_1d<double, 3>>>::Has(design_varible_name) == true)
        {
            const Variable<array_1d<double, 3>>& r_variable =
                KratosComponents<Variable<array_1d<double, 3>>>::Get(design_varible_name);
            this->CalculateStressDesignVariableDerivative(r_variable, STRESS_ON_NODE, rOutput, rCurrentProcessInfo);
        }
    }
        else
    {
        rOutput.clear();
    }

    KRATOS_CATCH("")
}

void AdjointFiniteDifferencingBaseElement::GetValueOnIntegrationPoints(const Variable<double>& rVariable,
                    std::vector<double>& rValues,
                    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    if(this->Has(rVariable))
    {
        // Get result value for output
        double output_value = this->GetValue(rVariable);

        // Resize Output
        const unsigned int&  write_points_number = GetGeometry()
            .IntegrationPointsNumber(this->GetIntegrationMethod());
        if (rValues.size() != write_points_number)
        {
            rValues.resize(write_points_number);
        }

        // Write scalar result value on all Gauss-Points
        for(unsigned int i = 0; i < write_points_number; ++i)
        {
            rValues[i] = output_value;
        }
    }
    else
        KRATOS_ERROR << "Unsupported output variable." << std::endl;

    KRATOS_CATCH("")
}

int AdjointFiniteDifferencingBaseElement::Check(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    GeometryType& r_geom = GetGeometry();

    // verify that the variables are correctly initialized
    KRATOS_CHECK_VARIABLE_KEY(DISPLACEMENT);
    KRATOS_CHECK_VARIABLE_KEY(ROTATION);
    KRATOS_CHECK_VARIABLE_KEY(VELOCITY);
    KRATOS_CHECK_VARIABLE_KEY(ACCELERATION);
    KRATOS_CHECK_VARIABLE_KEY(DENSITY);
    KRATOS_CHECK_VARIABLE_KEY(CONSTITUTIVE_LAW);
    KRATOS_CHECK_VARIABLE_KEY(ADJOINT_DISPLACEMENT);
    KRATOS_CHECK_VARIABLE_KEY(ADJOINT_ROTATION);

    // TODO generic way of doing these checks without checking the dofs..

    // Check dofs
    for (unsigned int i = 0; i < r_geom.size(); i++)
    {
        auto& r_node = r_geom[i];

        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DISPLACEMENT, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(ROTATION, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(ADJOINT_DISPLACEMENT, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(ADJOINT_ROTATION, r_node);

        KRATOS_CHECK_DOF_IN_NODE(ADJOINT_DISPLACEMENT_X, r_node);
        KRATOS_CHECK_DOF_IN_NODE(ADJOINT_DISPLACEMENT_Y, r_node);
        KRATOS_CHECK_DOF_IN_NODE(ADJOINT_DISPLACEMENT_Z, r_node);
        KRATOS_CHECK_DOF_IN_NODE(ADJOINT_ROTATION_X, r_node);
        KRATOS_CHECK_DOF_IN_NODE(ADJOINT_ROTATION_Y, r_node);
        KRATOS_CHECK_DOF_IN_NODE(ADJOINT_ROTATION_Z, r_node);
    }

    return 0;

    KRATOS_CATCH("")
}

// Sensitivity functions

void AdjointFiniteDifferencingBaseElement::CalculateSensitivityMatrix(const Variable<double>& rDesignVariable, Matrix& rOutput,
                                            const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    // define working variables
    Vector RHS_undist;
    Vector RHS_dist;
    ProcessInfo copy_process_info = rCurrentProcessInfo;

    // Compute RHS before disturbing
    mpPrimalElement->CalculateRightHandSide(RHS_undist, copy_process_info);
    rOutput.resize(1,RHS_undist.size());

    // Get disturbance measure
    double delta = this->GetValue(DISTURBANCE_MEASURE);
    double correction_factor = this->GetDisturbanceMeasureCorrectionFactor(rDesignVariable);
    delta *= correction_factor;

    if ( mpPrimalElement->GetProperties().Has(rDesignVariable) )
    {
        // Save properties and its pointer
        Properties& r_global_property = mpPrimalElement->GetProperties();
        Properties::Pointer p_global_properties = mpPrimalElement->pGetProperties();

        // Create new property and assign it to the element
        Properties::Pointer p_local_property(new Properties(r_global_property));
        mpPrimalElement->SetProperties(p_local_property);

        // Disturb the design variable
        const double current_property_value = mpPrimalElement->GetProperties()[rDesignVariable];
        p_local_property->SetValue(rDesignVariable, (current_property_value + delta));

        this->AfterPerturbation(rDesignVariable, rCurrentProcessInfo);

        // Compute RHS after disturbance
        mpPrimalElement->CalculateRightHandSide(RHS_dist, copy_process_info);

        // Compute derivative of RHS w.r.t. design variable with finite differences
        noalias(RHS_dist) -= RHS_undist;
        RHS_dist /= delta;
        for(unsigned int i = 0; i < RHS_dist.size(); i++)
            rOutput(0, i) = RHS_dist[i];

        // Give element original properties back
        mpPrimalElement->SetProperties(p_global_properties);

        this->AfterPerturbation(rDesignVariable, rCurrentProcessInfo);

        //call one last time to make sure everything is as it was before TODO improve this..
        mpPrimalElement->CalculateRightHandSide(RHS_dist, copy_process_info);
    }
    else
        rOutput.clear();

    KRATOS_CATCH("")

}

void AdjointFiniteDifferencingBaseElement::CalculateSensitivityMatrix(const Variable<array_1d<double,3>>& rDesignVariable, Matrix& rOutput,
                                            const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    // define working variables
    Vector RHS_undist;
    Vector RHS_dist;
    ProcessInfo copy_process_info = rCurrentProcessInfo;

    // Get disturbance measure
    double delta= this->GetValue(DISTURBANCE_MEASURE);
    double correction_factor = this->GetDisturbanceMeasureCorrectionFactor(rDesignVariable);
    delta *= correction_factor;

    if(rDesignVariable == SHAPE)
    {
        const int number_of_nodes = mpPrimalElement->GetGeometry().PointsNumber();
        const unsigned int dimension = rCurrentProcessInfo.GetValue(DOMAIN_SIZE);
        const int local_size = number_of_nodes * dimension * 2;

        rOutput.resize(dimension * number_of_nodes, local_size);

        // compute RHS before disturbing
        mpPrimalElement->CalculateRightHandSide(RHS_undist, copy_process_info);

        int index = 0;
        //TODO: look that this works also for parallel computing
        for(auto& node_i : mpPrimalElement->GetGeometry())
        {
            for(std::size_t coord_dir_i = 0; coord_dir_i < dimension; coord_dir_i++)
            {
                // disturb the design variable
                node_i.GetInitialPosition()[coord_dir_i] += delta;

                this->AfterPerturbation(rDesignVariable, rCurrentProcessInfo);

                // compute RHS after disturbance
                mpPrimalElement->CalculateRightHandSide(RHS_dist, copy_process_info);

                //compute derivative of RHS w.r.t. design variable with finite differences
                noalias(RHS_dist) -= RHS_undist;
                RHS_dist /= delta;
                for(unsigned int i = 0; i < RHS_dist.size(); i++)
                    rOutput( (coord_dir_i + index*dimension), i) = RHS_dist[i];

                // Reset perturbed vector
                noalias(RHS_dist) = ZeroVector(RHS_dist.size());

                // undisturb the design variable
                node_i.GetInitialPosition()[coord_dir_i] -= delta;

                this->AfterPerturbation(rDesignVariable, rCurrentProcessInfo);
            }
            index++;

            //call one last time to make sure everything is as it was before TODO improve this..
            mpPrimalElement->CalculateRightHandSide(RHS_dist, copy_process_info);

        }// end loop over element nodes

    }
    else
        KRATOS_ERROR << "Unsupported design variable!" << std::endl;

    KRATOS_CATCH("")
}

void AdjointFiniteDifferencingBaseElement::CalculateStressDisplacementDerivative(const Variable<Vector>& rStressVariable,
                                            Matrix& rOutput, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    const int num_nodes = mpPrimalElement->GetGeometry().PointsNumber();
    const int dimension = mpPrimalElement->GetGeometry().WorkingSpaceDimension();
    const int num_dofs_per_node = dimension * 2; // TODO generalize
    const int num_dofs = num_nodes * num_dofs_per_node;
    Vector initial_state_variables;
    Vector stress_derivatives_vector;

    // TODO first calculation only to get the size of the stress vector
    this->Calculate(rStressVariable, stress_derivatives_vector, rCurrentProcessInfo);
    rOutput.resize(num_dofs, stress_derivatives_vector.size() );
    rOutput.clear();
    initial_state_variables.resize(num_dofs);

    // Build vector of variables containing the DOF-variables of the primal problem
    std::vector<VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>>> primal_solution_variable_list;
    primal_solution_variable_list.reserve(num_dofs_per_node);
    primal_solution_variable_list.push_back(DISPLACEMENT_X);
    primal_solution_variable_list.push_back(DISPLACEMENT_Y);
    primal_solution_variable_list.push_back(DISPLACEMENT_Z);
    primal_solution_variable_list.push_back(ROTATION_X);
    primal_solution_variable_list.push_back(ROTATION_Y);
    primal_solution_variable_list.push_back(ROTATION_Z);

    // TODO Armin
    KRATOS_ERROR_IF(rCurrentProcessInfo.Has(NL_ITERATION_NUMBER))
        << "Stress displacement derivative computation is currently only available for linear cases !" << std::endl;

    for (int i = 0; i < num_nodes; i++)
    {
        int index = i * num_dofs_per_node;
        for(unsigned int j = 0; j < primal_solution_variable_list.size(); j++)
        {
            initial_state_variables[index + j] = mpPrimalElement->GetGeometry()[i].FastGetSolutionStepValue(primal_solution_variable_list[j]);
            mpPrimalElement->GetGeometry()[i].FastGetSolutionStepValue(primal_solution_variable_list[j]) = 0.0;
        }
    }
    for (int i = 0; i < num_nodes; i++)
    {
        int index = i * num_dofs_per_node;
        for(unsigned int j = 0; j < primal_solution_variable_list.size(); j++)
        {
            mpPrimalElement->GetGeometry()[i].FastGetSolutionStepValue(primal_solution_variable_list[j]) = 1.0;

            this->Calculate(rStressVariable, stress_derivatives_vector, rCurrentProcessInfo);

            for(unsigned int k = 0; k < stress_derivatives_vector.size(); k++)
                rOutput(index+j, k) = stress_derivatives_vector[k];

            stress_derivatives_vector.clear();

            mpPrimalElement->GetGeometry()[i].FastGetSolutionStepValue(primal_solution_variable_list[j]) = 0.0;
        }
    }
    for (int i = 0; i < num_nodes; i++)
    {
        int index = i * num_dofs_per_node;
        for(unsigned int j = 0; j < primal_solution_variable_list.size(); j++)
            mpPrimalElement->GetGeometry()[i].FastGetSolutionStepValue(primal_solution_variable_list[j]) = initial_state_variables[index + j];
    }

    KRATOS_CATCH("")
}

void AdjointFiniteDifferencingBaseElement::CalculateStressDesignVariableDerivative(const Variable<double>& rDesignVariable,
                                                const Variable<Vector>& rStressVariable, Matrix& rOutput,
                                                const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    // Define working variables
    Vector stress_vector_undist;
    Vector stress_vector_dist;

    // Compute stress on GP before disturbance
    this->Calculate(rStressVariable, stress_vector_undist, rCurrentProcessInfo);

    // Get disturbance measure
    double delta= this->GetValue(DISTURBANCE_MEASURE);
    double correction_factor = this->GetDisturbanceMeasureCorrectionFactor(rDesignVariable);
    delta *= correction_factor;

    const SizeType stress_vector_size = stress_vector_undist.size();
    rOutput.resize(1, stress_vector_size);

    if( mpPrimalElement->GetProperties().Has(rDesignVariable) )
    {
        // Save properties and its pointer
        Properties& r_global_property = mpPrimalElement->GetProperties();
        Properties::Pointer p_global_properties = mpPrimalElement->pGetProperties();

        // Create new property and assign it to the element
        Properties::Pointer p_local_property(new Properties(r_global_property));
        mpPrimalElement->SetProperties(p_local_property);

        // Disturb the design variable
        const double current_property_value = mpPrimalElement->GetProperties()[rDesignVariable];
        p_local_property->SetValue(rDesignVariable, (current_property_value + delta));

        this->AfterPerturbation(rDesignVariable, rCurrentProcessInfo);

        // Compute stress on GP after disturbance
        this->Calculate(rStressVariable, stress_vector_dist, rCurrentProcessInfo);

        // Compute derivative of stress w.r.t. design variable with finite differences
        noalias(stress_vector_dist)  -= stress_vector_undist;
        stress_vector_dist  /= delta;

        for(size_t j = 0; j < stress_vector_size; j++)
            rOutput(0, j) = stress_vector_dist[j];

        // Give element original properties back
        mpPrimalElement->SetProperties(p_global_properties);

        this->AfterPerturbation(rDesignVariable, rCurrentProcessInfo);
    }
    else
        rOutput.clear();

    KRATOS_CATCH("")
}

void AdjointFiniteDifferencingBaseElement::CalculateStressDesignVariableDerivative(const Variable<array_1d<double,3>>& rDesignVariable,
                                            const Variable<Vector>& rStressVariable,
                                            Matrix& rOutput, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    // define working variables
    Vector stress_vector_undist;
    Vector stress_vector_dist;

    // Get disturbance measure
    double delta= this->GetValue(DISTURBANCE_MEASURE);
    double correction_factor = this->GetDisturbanceMeasureCorrectionFactor(rDesignVariable);
    delta *= correction_factor;

    if(rDesignVariable == SHAPE)
    {
        const int number_of_nodes = mpPrimalElement->GetGeometry().PointsNumber();
        const unsigned int dimension = rCurrentProcessInfo.GetValue(DOMAIN_SIZE);

        // Compute stress on GP before disturbance
        this->Calculate(rStressVariable, stress_vector_undist, rCurrentProcessInfo);

        const SizeType stress_vector_size = stress_vector_undist.size();
        rOutput.resize(dimension * number_of_nodes, stress_vector_size);

        int index = 0;
        //TODO: look that this works also for parallel computing
        for(auto& node_i : mpPrimalElement->GetGeometry())
        {
            for(std::size_t coord_dir_i = 0; coord_dir_i < dimension; coord_dir_i++)
            {
                // disturb the design variable
                node_i.GetInitialPosition()[coord_dir_i] += delta;

                this->AfterPerturbation(rDesignVariable, rCurrentProcessInfo);

                // Compute stress on GP after disturbance
                this->Calculate(rStressVariable, stress_vector_dist, rCurrentProcessInfo);

                // Compute derivative of stress w.r.t. design variable with finite differences
                noalias(stress_vector_dist)  -= stress_vector_undist;
                stress_vector_dist  /= delta;

                for(size_t i = 0; i < stress_vector_size; i++)
                    rOutput( (coord_dir_i + index*dimension), i) = stress_vector_dist[i];

                // Reset pertubed vector
                stress_vector_dist = Vector(0);

                // undisturb the design variable
                node_i.GetInitialPosition()[coord_dir_i] -= delta;

                this->AfterPerturbation(rDesignVariable, rCurrentProcessInfo);
            }
            index++;
        }// end loop over element nodes
    }
    else
        KRATOS_ERROR << "Unsupported design variable!" << std::endl;

    KRATOS_CATCH("")
}

// private

double AdjointFiniteDifferencingBaseElement::GetDisturbanceMeasureCorrectionFactor(const Variable<double>& rDesignVariable)
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

double AdjointFiniteDifferencingBaseElement::GetDisturbanceMeasureCorrectionFactor(const Variable<array_1d<double,3>>& rDesignVariable)
{
    KRATOS_TRY;

    if(rDesignVariable == SHAPE)
    {
        KRATOS_ERROR << "GetDisturbanceMeasureCorrectionFactor NOT_IMPLEMENTED" << std::endl;
    }
    else
        return 1.0;

    KRATOS_CATCH("")
}

void AdjointFiniteDifferencingBaseElement::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer,  Element );
}

void AdjointFiniteDifferencingBaseElement::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Element );

}

} // namespace Kratos


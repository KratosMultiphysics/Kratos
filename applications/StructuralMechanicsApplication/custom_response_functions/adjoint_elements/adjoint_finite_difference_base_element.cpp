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

AdjointFiniteDifferencingBaseElement::AdjointFiniteDifferencingBaseElement(IndexType NewId,
                        GeometryType::Pointer pGeometry)
                        : Element(NewId, pGeometry) {}

AdjointFiniteDifferencingBaseElement::AdjointFiniteDifferencingBaseElement(IndexType NewId,
                        GeometryType::Pointer pGeometry,
                        PropertiesType::Pointer pProperties, Element::Pointer pPrimalElement)
                    : Element(NewId, pGeometry, pProperties) {

    mpPrimalElement = pPrimalElement;
    mpPrimalShellElement = dynamic_pointer_cast<ShellThinElement3D3N>(pPrimalElement);
    }

AdjointFiniteDifferencingBaseElement::~AdjointFiniteDifferencingBaseElement() {}

Element::Pointer AdjointFiniteDifferencingBaseElement::Create(IndexType NewId,
                NodesArrayType const& rThisNodes, PropertiesType::Pointer pProperties,
                Element::Pointer pPrimalElement) const
{
    const GeometryType& rGeom = this->GetGeometry();
    return Kratos::make_shared<AdjointFiniteDifferencingBaseElement>(
        NewId, rGeom.Create(rThisNodes), pProperties, pPrimalElement);
}

void AdjointFiniteDifferencingBaseElement::EquationIdVector(EquationIdVectorType& rResult,
    ProcessInfo& rCurrentProcessInfo)
{
    const SizeType num_dofs = mpPrimalShellElement->GetNumberOfDofs();
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

    const SizeType num_dofs = mpPrimalShellElement->GetNumberOfDofs();
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

    const SizeType num_dofs = mpPrimalShellElement->GetNumberOfDofs();
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

void AdjointFiniteDifferencingBaseElement::Calculate(const Variable<Vector >& rVariable,
                           Vector& rOutput,
                           const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    const SizeType num_gps = mpPrimalShellElement->GetNumberOfGPs();

    if(rVariable == STRESS_ON_GP)
    {
        TracedStressType traced_stress_type = static_cast<TracedStressType>(this->GetValue(TRACED_STRESS_TYPE));

        int direction_1 = 0;
        int direction_2 = 0;
        std::vector<Matrix> stress_vector;
        bool stress_is_moment = true;

        switch (traced_stress_type)
        {
            case TracedStressType::MXX:
            {
                direction_1 = 0;
                direction_2 = 0;
                break;
            }
            case TracedStressType::MXY:
            {
                direction_1 = 0;
                direction_2 = 1;
                break;
            }
            case TracedStressType::MXZ:
            {
                direction_1 = 0;
                direction_2 = 2;
                break;
            }
            case TracedStressType::MYX:
            {
                direction_1 = 1;
                direction_2 = 0;
                break;
            }
            case TracedStressType::MYY :
            {
                direction_1 = 1;
                direction_2 = 1;
                break;
            }
            case TracedStressType::MYZ:
            {
                direction_1 = 1;
                direction_2 = 2;
                break;
            }
            case TracedStressType::MZX:
            {
                direction_1 = 2;
                direction_2 = 0;
                break;
            }
            case TracedStressType::MZY:
            {
                direction_1 = 2;
                direction_2 = 1;
                break;
            }
            case TracedStressType::MZZ :
            {
                direction_1 = 2;
                direction_2 = 2;
                break;
            }
            case TracedStressType::FXX :
            {
                direction_1 = 0;
                direction_2 = 0;
                stress_is_moment = false;
                break;
            }
            case TracedStressType::FXY:
            {
                direction_1 = 0;
                direction_2 = 1;
                stress_is_moment = false;
                break;
            }
            case TracedStressType::FXZ:
            {
                direction_1 = 0;
                direction_2 = 2;
                stress_is_moment = false;
                break;
            }
            case TracedStressType::FYX:
            {
                direction_1 = 1;
                direction_2 = 0;
                stress_is_moment = false;
                break;
            }
            case TracedStressType::FYY:
            {
                direction_1 = 1;
                direction_2 = 1;
                stress_is_moment = false;
                break;
            }
            case TracedStressType::FYZ:
            {
                direction_1 = 1;
                direction_2 = 2;
                stress_is_moment = false;
                break;
            }
            case TracedStressType::FZX:
            {
                direction_1 = 2;
                direction_2 = 0;
                stress_is_moment = false;
                break;
            }
            case TracedStressType::FZY:
            {
                direction_1 = 2;
                direction_2 = 1;
                stress_is_moment = false;
                break;
            }
            case TracedStressType::FZZ:
            {
                direction_1 = 2;
                direction_2 = 2;
                stress_is_moment = false;
                break;
            }
            default:
                KRATOS_ERROR << "Invalid stress type! Stress type not supported for this element!" << std::endl;
        }

        if(stress_is_moment)
            mpPrimalElement->GetValueOnIntegrationPoints(SHELL_MOMENT_GLOBAL, stress_vector, rCurrentProcessInfo);
        else
            mpPrimalElement->GetValueOnIntegrationPoints(SHELL_FORCE_GLOBAL, stress_vector, rCurrentProcessInfo);

        rOutput.resize(num_gps);
        for(size_t i = 0; i < num_gps; i++)
        {
            rOutput(i) = stress_vector[i](direction_1, direction_2);
        }

    }
    else
    {
        rOutput.resize(num_gps);
        rOutput.clear();
    }

    KRATOS_CATCH("")
}

void AdjointFiniteDifferencingBaseElement::Calculate(const Variable<Matrix >& rVariable, Matrix& rOutput,
                                                const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    if(rVariable == STRESS_DISP_DERIV_ON_GP)
    {
       this->CalculateStressDisplacementDerivative(STRESS_ON_GP, rOutput, rCurrentProcessInfo);
    }
    else if(rVariable == STRESS_DESIGN_DERIVATIVE_ON_GP)
    {
        std::string design_variable_name = this->GetValue( DESIGN_VARIABLE_NAME );

        if (KratosComponents<Variable<double>>::Has(design_variable_name) == true)
        {
            const Variable<double>& r_variable =
                KratosComponents<Variable<double>>::Get(design_variable_name);
            this->CalculateStressDesignVariableDerivative(r_variable, STRESS_ON_GP, rOutput, rCurrentProcessInfo);
        }
        else if (KratosComponents<Variable<array_1d<double, 3>>>::Has(design_variable_name) == true)
        {
            const Variable<array_1d<double, 3>>& r_variable =
                KratosComponents<Variable<array_1d<double, 3>>>::Get(design_variable_name);
            this->CalculateStressDesignVariableDerivative(r_variable, STRESS_ON_GP, rOutput, rCurrentProcessInfo);
        }
    }
    else
    {
        rOutput.clear();
    }

    KRATOS_CATCH("")
}

void AdjointFiniteDifferencingBaseElement::CalculateOnIntegrationPoints(const Variable<double>& rVariable,
                          std::vector<double>& rOutput,
                          const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    if(this->Has(rVariable))
    {
        // Get result value for output
        double output_value = this->GetValue(rVariable);
        const SizeType num_gps = mpPrimalShellElement->GetNumberOfGPs(); //TODO?

        // Resize Output
        if(rOutput.size() != num_gps)
            rOutput.resize(num_gps);

        // Write scalar result value on all Gauss-Points
        for(unsigned int i = 0; i < num_gps; i++)
            rOutput[i] = output_value;

        //OPT_INTERPOLATE_RESULTS_TO_STANDARD_GAUSS_POINTS(rOutput);
    }
    else
        KRATOS_ERROR << "Unsupported output variable." << std::endl;



    KRATOS_CATCH("")

}

void AdjointFiniteDifferencingBaseElement::GetValueOnIntegrationPoints(const Variable<double>& rVariable,
                    std::vector<double>& rValues,
                    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;
    this->CalculateOnIntegrationPoints(rVariable, rValues, rCurrentProcessInfo);
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
    KRATOS_CHECK_VARIABLE_KEY(SHELL_CROSS_SECTION);
    KRATOS_CHECK_VARIABLE_KEY(THICKNESS);
    KRATOS_CHECK_VARIABLE_KEY(CONSTITUTIVE_LAW);
    KRATOS_CHECK_VARIABLE_KEY(ADJOINT_DISPLACEMENT);
    KRATOS_CHECK_VARIABLE_KEY(ADJOINT_ROTATION);

    // check properties
    KRATOS_ERROR_IF(this->pGetProperties() == nullptr) << "Properties not provided for element " << this->Id() << std::endl;

    const PropertiesType & props = this->GetProperties();

    if(props.Has(SHELL_CROSS_SECTION)) // if the user specified a cross section ...
    {
        const ShellCrossSection::Pointer & section = props[SHELL_CROSS_SECTION];
        KRATOS_ERROR_IF(section == nullptr) << "SHELL_CROSS_SECTION not provided for element " << this->Id() << std::endl;

        section->Check(props, r_geom, rCurrentProcessInfo);
    }
    else // ... allow the automatic creation of a homogeneous section from a material and a thickness
    {
        KRATOS_ERROR_IF_NOT(props.Has(CONSTITUTIVE_LAW)) << "CONSTITUTIVE_LAW not provided for element " << this->Id() << std::endl;
        const ConstitutiveLaw::Pointer& claw = props[CONSTITUTIVE_LAW];
        KRATOS_ERROR_IF(claw == nullptr) << "CONSTITUTIVE_LAW not provided for element " << this->Id() << std::endl;

        KRATOS_ERROR_IF_NOT(props.Has(THICKNESS)) <<  "THICKNESS not provided for element " <<  this->Id() << std::endl;
        KRATOS_ERROR_IF(props[THICKNESS] <= 0.0) << "wrong THICKNESS value provided for element " << this->Id() << std::endl;

        ShellCrossSection::Pointer dummySection = ShellCrossSection::Pointer(new ShellCrossSection());
        dummySection->BeginStack();
        dummySection->AddPly(props[THICKNESS], 0.0, 5, this->pGetProperties());
        dummySection->EndStack();
        dummySection->SetSectionBehavior(ShellCrossSection::Thin);
        dummySection->Check(props, r_geom, rCurrentProcessInfo);
    }

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
        // TODO ActualizeElementWithDisturbance

        // Save properties and its pointer
        Properties& r_global_property = mpPrimalElement->GetProperties();
        Properties::Pointer p_global_properties = mpPrimalElement->pGetProperties();

        // Create new property and assign it to the element
        Properties::Pointer p_local_property(new Properties(r_global_property));
        mpPrimalElement->SetProperties(p_local_property);

        // Disturb the design variable
        const double current_property_value = mpPrimalElement->GetProperties()[rDesignVariable];
        p_local_property->SetValue(rDesignVariable, (current_property_value + delta));

        // TODO ActualizeElementWithDisturbance
        mpPrimalShellElement->ResetSections(); // TODO move to separate function at derived element
        mpPrimalElement->Initialize(); // TODO move to separate function at derived element

        // Compute RHS after disturbance
        mpPrimalElement->CalculateRightHandSide(RHS_dist, copy_process_info);

        // Compute derivative of RHS w.r.t. design variable with finite differences
        noalias(RHS_dist) -= RHS_undist;
        RHS_dist /= delta;
        for(unsigned int i = 0; i < RHS_dist.size(); i++)
            rOutput(0, i) = RHS_dist[i];

        // Give element original properties back
        mpPrimalElement->SetProperties(p_global_properties);
        mpPrimalShellElement->ResetSections();
        mpPrimalElement->Initialize();
        mpPrimalElement->CalculateRightHandSide(RHS_dist, copy_process_info);

        // TODO FinalizeDisturbance
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
                    // TODO InitializeDisturbance

                    // disturb the design variable
                    node_i.GetInitialPosition()[coord_dir_i] += delta;

                    // TODO ActualizeElementWithDisturbance

                    // compute RHS after disturbance
                    mpPrimalElement->CalculateRightHandSide(RHS_dist, copy_process_info);

                    //compute derivative of RHS w.r.t. design variable with finite differences
                    noalias(RHS_dist) -= RHS_undist;
                    RHS_dist /= delta;
                    for(unsigned int i = 0; i < RHS_dist.size(); i++)
                        rOutput( (coord_dir_i + index*dimension), i) = RHS_dist[i];

                    // Reset pertubed vector
                    RHS_dist = Vector(0);

                    // undisturb the design variable
                    node_i.GetInitialPosition()[coord_dir_i] -= delta;

                    // TODO FinalizeDisturbance

                }
                index++;

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
    const SizeType num_dofs = mpPrimalShellElement->GetNumberOfDofs();
    const SizeType num_gps = mpPrimalShellElement->GetNumberOfGPs();
    ProcessInfo copy_process_info = rCurrentProcessInfo;
    Vector initial_state_variables;
    Vector stress_derivatives_vector;

    rOutput.resize(num_dofs, num_gps);
    rOutput.clear();
    initial_state_variables.resize(num_dofs);

    // Built vector of variables containing the DOF-variables of the primal problem
    std::vector<VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>>> primal_solution_variable_list;
    primal_solution_variable_list.push_back(DISPLACEMENT_X);
    primal_solution_variable_list.push_back(DISPLACEMENT_Y);
    primal_solution_variable_list.push_back(DISPLACEMENT_Z);
    primal_solution_variable_list.push_back(ROTATION_X);
    primal_solution_variable_list.push_back(ROTATION_Y);
    primal_solution_variable_list.push_back(ROTATION_Z);

    // TODO Armin
    KRATOS_ERROR_IF(rCurrentProcessInfo.Has(NL_ITERATION_NUMBER))
        << "Stress displacement derivative computation is currently only for linear cases availible!" << std::endl;

    for (int i = 0; i < num_nodes; i++)
    {
        int index = i * dimension * 2;
        for(unsigned int j = 0; j < primal_solution_variable_list.size(); j++)
        {
            initial_state_variables[index + j] = mpPrimalElement->GetGeometry()[i].FastGetSolutionStepValue(primal_solution_variable_list[j]);
            mpPrimalElement->GetGeometry()[i].FastGetSolutionStepValue(primal_solution_variable_list[j]) = 0.0;
        }
    }
    for (int i = 0; i < num_nodes; i++)
    {
        int index = i * dimension * 2;
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
        int index = i * dimension * 2;
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

    const SizeType num_gps = mpPrimalShellElement->GetNumberOfGPs();
    rOutput.resize(1, num_gps);

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

        // TODO
        mpPrimalShellElement->ResetSections();
        mpPrimalElement->Initialize();

        // Compute stress on GP after disturbance
        this->Calculate(rStressVariable, stress_vector_dist, rCurrentProcessInfo);

        // Compute derivative of stress w.r.t. design variable with finite differences
        noalias(stress_vector_dist)  -= stress_vector_undist;
        stress_vector_dist  /= delta;

        for(size_t j = 0; j < num_gps; j++)
            rOutput(0, j) = stress_vector_dist[j];

        // Give element original properties back
        mpPrimalElement->SetProperties(p_global_properties);

        // TODO
        mpPrimalShellElement->ResetSections();
        mpPrimalElement->Initialize();
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

        const SizeType num_gps = mpPrimalShellElement->GetNumberOfGPs();
        rOutput.resize(dimension * number_of_nodes, num_gps);

        // Compute stress on GP before disturbance
        this->Calculate(rStressVariable, stress_vector_undist, rCurrentProcessInfo);

        int index = 0;
        //TODO: look that this works also for parallel computing
        for(auto& node_i : mpPrimalElement->GetGeometry())
        {
            for(std::size_t coord_dir_i = 0; coord_dir_i < dimension; coord_dir_i++)
            {
                // disturb the design variable
                node_i.GetInitialPosition()[coord_dir_i] += delta;

                // Compute stress on GP after disturbance
                this->Calculate(rStressVariable, stress_vector_dist, rCurrentProcessInfo);

                // Compute derivative of stress w.r.t. design variable with finite differences
                noalias(stress_vector_dist)  -= stress_vector_undist;
                stress_vector_dist  /= delta;

                for(size_t i = 0; i < num_gps; i++)
                    rOutput( (coord_dir_i + index*dimension), i) = stress_vector_dist[i];

                // Reset pertubed vector
                stress_vector_dist = Vector(0);

                // undisturb the design variable
                node_i.GetInitialPosition()[coord_dir_i] -= delta;
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
        double dx, dy, dz, L = 0.0;

        const GeometryType& geometry = mpPrimalElement->GetGeometry();

        dx = geometry[1].X0() - geometry[0].X0();
        dy = geometry[1].Y0() - geometry[0].Y0();
        dz = geometry[1].Z0() - geometry[0].Z0();
        L += sqrt(dx*dx + dy*dy + dz*dz);
        dx = geometry[2].X0() - geometry[1].X0();
        dy = geometry[2].Y0() - geometry[1].Y0();
        dz = geometry[2].Z0() - geometry[1].Z0();
        L += sqrt(dx*dx + dy*dy + dz*dz);
        dx = geometry[2].X0() - geometry[0].X0();
        dy = geometry[2].Y0() - geometry[0].Y0();
        dz = geometry[2].Z0() - geometry[0].Z0();
        L += sqrt(dx*dx + dy*dy + dz*dz);
        L /= 3.0;

        return L;
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


// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Martin Fusseder, https://github.com/MFusseder
//

#include "shell_thin_adjoint_element_3D3N.hpp"
#include "custom_utilities/shellt3_corotational_coordinate_transformation.hpp"
#include "structural_mechanics_application_variables.h"
#include "custom_response_functions/response_utilities/response_data.h"
#include "includes/checks.h"

//#include <string>
//#include <iomanip>

//----------------------------------------
// preprocessors for the integration
// method used by this element.

//#define OPT_1_POINT_INTEGRATION

//----------------------------------------
// preprocessors to handle the output
// in case of 3 integration points

//#define OPT_USES_INTERIOR_GAUSS_POINTS

#ifdef OPT_1_POINT_INTEGRATION
#define OPT_INTERPOLATE_RESULTS_TO_STANDARD_GAUSS_POINTS(X)
#else
#ifdef OPT_USES_INTERIOR_GAUSS_POINTS
#define OPT_INTERPOLATE_RESULTS_TO_STANDARD_GAUSS_POINTS(X)
#else
#define OPT_INTERPOLATE_RESULTS_TO_STANDARD_GAUSS_POINTS(X) ShellUtilities::InterpToStandardGaussPoints(X)
#endif // OPT_USES_INTERIOR_GAUSS_POINTS
#endif // OPT_1_POINT_INTEGRATION

namespace Kratos
{

// =====================================================================================
//
// Class ShellThinAdjointElement3D3N
//
// =====================================================================================

ShellThinAdjointElement3D3N::ShellThinAdjointElement3D3N(IndexType NewId,
        GeometryType::Pointer pGeometry,
        bool NLGeom)
    : ShellThinElement3D3N(NewId, pGeometry, NLGeom)
{

}

ShellThinAdjointElement3D3N::ShellThinAdjointElement3D3N(IndexType NewId,
        GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties,
        bool NLGeom)
    : ShellThinElement3D3N(NewId, pGeometry, pProperties, NLGeom)
{

}

ShellThinAdjointElement3D3N::ShellThinAdjointElement3D3N(IndexType NewId,
        GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties,
        CoordinateTransformationBasePointerType pCoordinateTransformation)
    : ShellThinElement3D3N(NewId, pGeometry, pProperties, pCoordinateTransformation)
{

}

ShellThinAdjointElement3D3N::~ShellThinAdjointElement3D3N()
{
}

Element::Pointer ShellThinAdjointElement3D3N::Create(IndexType NewId,
	NodesArrayType const& ThisNodes,
	PropertiesType::Pointer pProperties) const
{
	GeometryType::Pointer newGeom(GetGeometry().Create(ThisNodes));
	return Kratos::make_shared< ShellThinAdjointElement3D3N >(NewId, newGeom,
		pProperties, ShellThinElement3D3N::mpCoordinateTransformation->Create(newGeom));
}

void ShellThinAdjointElement3D3N::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo)
{
    const SizeType num_dofs = GetNumberOfDofs();
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

void ShellThinAdjointElement3D3N::GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& CurrentProcessInfo)
{
    const SizeType num_dofs = GetNumberOfDofs();
    ElementalDofList.resize(0);
    ElementalDofList.reserve(num_dofs);

    GeometryType & geom = this->GetGeometry();

    for (SizeType i = 0; i < geom.size(); i++)
    {
        NodeType & iNode = geom[i];

        ElementalDofList.push_back(iNode.pGetDof(ADJOINT_DISPLACEMENT_X));
        ElementalDofList.push_back(iNode.pGetDof(ADJOINT_DISPLACEMENT_Y));
        ElementalDofList.push_back(iNode.pGetDof(ADJOINT_DISPLACEMENT_Z));

        ElementalDofList.push_back(iNode.pGetDof(ADJOINT_ROTATION_X));
        ElementalDofList.push_back(iNode.pGetDof(ADJOINT_ROTATION_Y));
        ElementalDofList.push_back(iNode.pGetDof(ADJOINT_ROTATION_Z));
    }
}

int ShellThinAdjointElement3D3N::Check(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    GeometryType& r_geom = GetGeometry();

    this->CheckVariables();
    KRATOS_CHECK_VARIABLE_KEY(ADJOINT_DISPLACEMENT);
    KRATOS_CHECK_VARIABLE_KEY(ADJOINT_ROTATION);

    this->CheckProperties(rCurrentProcessInfo);

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

void ShellThinAdjointElement3D3N::GetValuesVector(Vector& values, int Step)
{
    const SizeType num_dofs = GetNumberOfDofs();
    if(values.size() != num_dofs)
        values.resize(num_dofs, false);

    const GeometryType & geom = GetGeometry();
    const int dimension = geom.WorkingSpaceDimension();

    for (SizeType i = 0; i < geom.size(); i++)
    {
        const NodeType & iNode = geom[i];
        const array_1d<double,3>& disp = iNode.FastGetSolutionStepValue(ADJOINT_DISPLACEMENT, Step);
        const array_1d<double,3>& rot = iNode.FastGetSolutionStepValue(ADJOINT_ROTATION, Step);

        const SizeType index = i * dimension * 2;
        values[index]     = disp[0];
        values[index + 1] = disp[1];
        values[index + 2] = disp[2];

        values[index + 3] = rot[0];
        values[index + 4] = rot[1];
        values[index + 5] = rot[2];
    }
}

void ShellThinAdjointElement3D3N::CalculateSensitivityMatrix(const Variable<double>& rDesignVariable, Matrix& rOutput,
                                            const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    // define working variables
    Vector RHS_undist;
    Vector RHS_dist;
    ProcessInfo copy_process_info = rCurrentProcessInfo;

    // Compute RHS before disturbing
    this->CalculateRightHandSide(RHS_undist, copy_process_info);
    rOutput.resize(1,RHS_undist.size());

    // Get disturbance measure
    double delta = this->GetValue(DISTURBANCE_MEASURE);
    double correction_factor = this->GetDisturbanceMeasureCorrectionFactor(rDesignVariable);
    delta *= correction_factor;

    if ( this->GetProperties().Has(rDesignVariable) )
    {
        // Save properties and its pointer
        Properties& r_global_property = this->GetProperties();
        Properties::Pointer p_global_properties = this->pGetProperties();

        // Create new property and assign it to the element
        Properties::Pointer p_local_property(new Properties(r_global_property));
        this->SetProperties(p_local_property);

        // Disturb the design variable
        const double current_property_value = this->GetProperties()[rDesignVariable];
        p_local_property->SetValue(rDesignVariable, (current_property_value + delta));

        // Compute RHS after disturbance
        this->CalculateRightHandSide(RHS_dist, copy_process_info);

        // Compute derivative of RHS w.r.t. design variable with finite differences
        noalias(RHS_dist) -= RHS_undist;
        RHS_dist /= delta;
        for(unsigned int i = 0; i < RHS_dist.size(); i++)
            rOutput(0, i) = RHS_dist[i];

        // Give element original properties back
        this->SetProperties(p_global_properties);
        this->CalculateRightHandSide(RHS_dist, copy_process_info);
    }
    else
        rOutput.clear();

    KRATOS_CATCH("")

}

void ShellThinAdjointElement3D3N::CalculateSensitivityMatrix(const Variable<array_1d<double,3>>& rDesignVariable, Matrix& rOutput,
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
            const int number_of_nodes = GetGeometry().PointsNumber();
            const unsigned int dimension = rCurrentProcessInfo.GetValue(DOMAIN_SIZE);
            const int local_size = number_of_nodes * dimension * 2;

            rOutput.resize(dimension * number_of_nodes, local_size);

            // compute RHS before disturbing
            this->CalculateRightHandSide(RHS_undist, copy_process_info);

            int index = 0;
            //TODO: look that this works also for parallel computing
            for(auto& node_i : this->GetGeometry())
            {
                for(std::size_t coord_dir_i = 0; coord_dir_i < dimension; coord_dir_i++)
                {
                    // disturb the design variable
                    node_i.GetInitialPosition()[coord_dir_i] += delta;

                    // compute RHS after disturbance
                    this->CalculateRightHandSide(RHS_dist, copy_process_info);

                    //compute derivative of RHS w.r.t. design variable with finite differences
                    noalias(RHS_dist) -= RHS_undist;
                    RHS_dist /= delta;
                    for(unsigned int i = 0; i < RHS_dist.size(); i++)
                        rOutput( (coord_dir_i + index*dimension), i) = RHS_dist[i];

                    // Reset pertubed vector
                    RHS_dist = Vector(0);

                    // undisturb the design variable
                    node_i.GetInitialPosition()[coord_dir_i] -= delta;
                }
                index++;

                this->CalculateRightHandSide(RHS_dist, copy_process_info);

            }// end loop over element nodes
        }
        else
            KRATOS_ERROR << "Unsupported design variable!" << std::endl;

        KRATOS_CATCH("")

}

void ShellThinAdjointElement3D3N::Calculate(const Variable<Vector >& rVariable,
                           Vector& rOutput,
                           const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    const SizeType num_gps = GetNumberOfGPs();

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
            ShellThinElement3D3N::GetValueOnIntegrationPoints(SHELL_MOMENT_GLOBAL, stress_vector, rCurrentProcessInfo);
        else
            ShellThinElement3D3N::GetValueOnIntegrationPoints(SHELL_FORCE_GLOBAL, stress_vector, rCurrentProcessInfo);

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

void ShellThinAdjointElement3D3N::Calculate(const Variable<Matrix >& rVariable, Matrix& rOutput,
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

void ShellThinAdjointElement3D3N::CalculateLeftHandSide( MatrixType& rLeftHandSideMatrix,
        ProcessInfo& rCurrentProcessInfo)
{
    Vector dummy;
    ShellThinElement3D3N::CalculateLocalSystem(rLeftHandSideMatrix, dummy, rCurrentProcessInfo);
}
// =====================================================================================
//
// Class ShellThinAdjointElement3D3N - Results on Gauss Points
//
// =====================================================================================

void ShellThinAdjointElement3D3N::CalculateOnIntegrationPoints(const Variable<double>& rVariable,
                          std::vector<double>& rOutput,
                          const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    if(this->Has(rVariable))
    {
        // Get result value for output
        double output_value = this->GetValue(rVariable);
        const SizeType num_gps = GetNumberOfGPs();

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

void ShellThinAdjointElement3D3N::GetValueOnIntegrationPoints(const Variable<double>& rVariable,
                    std::vector<double>& rValues,
                    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;
    this->CalculateOnIntegrationPoints(rVariable, rValues, rCurrentProcessInfo);
    KRATOS_CATCH("")
}

double ShellThinAdjointElement3D3N::GetDisturbanceMeasureCorrectionFactor(const Variable<double>& rDesignVariable)
{
    KRATOS_TRY;

    if ( this->GetProperties().Has(rDesignVariable) )
    {
        const double variable_value = this->GetProperties()[rDesignVariable];
        return variable_value;
    }
    else
        return 1.0;

    KRATOS_CATCH("")
}

double ShellThinAdjointElement3D3N::GetDisturbanceMeasureCorrectionFactor(const Variable<array_1d<double,3>>& rDesignVariable)
{
    KRATOS_TRY;

    if(rDesignVariable == SHAPE)
    {
        double dx, dy, dz, L = 0.0;

        dx = this->GetGeometry()[1].X0() - this->GetGeometry()[0].X0();
        dy = this->GetGeometry()[1].Y0() - this->GetGeometry()[0].Y0();
        dz = this->GetGeometry()[1].Z0() - this->GetGeometry()[0].Z0();
        L += sqrt(dx*dx + dy*dy + dz*dz);
        dx = this->GetGeometry()[2].X0() - this->GetGeometry()[1].X0();
        dy = this->GetGeometry()[2].Y0() - this->GetGeometry()[1].Y0();
        dz = this->GetGeometry()[2].Z0() - this->GetGeometry()[1].Z0();
        L += sqrt(dx*dx + dy*dy + dz*dz);
        dx = this->GetGeometry()[2].X0() - this->GetGeometry()[0].X0();
        dy = this->GetGeometry()[2].Y0() - this->GetGeometry()[0].Y0();
        dz = this->GetGeometry()[2].Z0() - this->GetGeometry()[0].Z0();
        L += sqrt(dx*dx + dy*dy + dz*dz);
        L /= 3.0;

        return L;
    }
    else
        return 1.0;

    KRATOS_CATCH("")
}

void ShellThinAdjointElement3D3N::CalculateStressDisplacementDerivative(const Variable<Vector>& rStressVariable,
                                            Matrix& rOutput, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    const int num_nodes = this->GetGeometry().PointsNumber();
    const int dimension = this->GetGeometry().WorkingSpaceDimension();
    const SizeType num_dofs = GetNumberOfDofs();
    const SizeType num_gps = GetNumberOfGPs();
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

    KRATOS_ERROR_IF(rCurrentProcessInfo.Has(NL_ITERATION_NUMBER))
        << "Stress displacement derivative computation is currently only for linear cases available!" << std::endl;

    for (int i = 0; i < num_nodes; i++)
    {
        int index = i * dimension * 2;
        for(unsigned int j = 0; j < primal_solution_variable_list.size(); j++)
        {
            initial_state_variables[index + j] = this->GetGeometry()[i].FastGetSolutionStepValue(primal_solution_variable_list[j]);
            this->GetGeometry()[i].FastGetSolutionStepValue(primal_solution_variable_list[j]) = 0.0;
        }
    }
    for (int i = 0; i < num_nodes; i++)
    {
        int index = i * dimension * 2;
        for(unsigned int j = 0; j < primal_solution_variable_list.size(); j++)
        {
            this->GetGeometry()[i].FastGetSolutionStepValue(primal_solution_variable_list[j]) = 1.0;

            this->Calculate(rStressVariable, stress_derivatives_vector, rCurrentProcessInfo);

            for(unsigned int k = 0; k < stress_derivatives_vector.size(); k++)
                rOutput(index+j, k) = stress_derivatives_vector[k];

            stress_derivatives_vector.clear();

            this->GetGeometry()[i].FastGetSolutionStepValue(primal_solution_variable_list[j]) = 0.0;
        }
    }
    for (int i = 0; i < num_nodes; i++)
    {
        int index = i * dimension * 2;
        for(unsigned int j = 0; j < primal_solution_variable_list.size(); j++)
            this->GetGeometry()[i].FastGetSolutionStepValue(primal_solution_variable_list[j]) = initial_state_variables[index + j];
    }

    KRATOS_CATCH("")
}

void ShellThinAdjointElement3D3N::CalculateStressDesignVariableDerivative(const Variable<double>& rDesignVariable,
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

        const SizeType num_gps = GetNumberOfGPs();
        rOutput.resize(1, num_gps);

        if( this->GetProperties().Has(rDesignVariable) )
        {
            // Save properties and its pointer
            Properties& r_global_property = this->GetProperties();
            Properties::Pointer p_global_properties = this->pGetProperties();

            // Create new property and assign it to the element
            Properties::Pointer p_local_property(new Properties(r_global_property));
            this->SetProperties(p_local_property);

            // Disturb the design variable
            const double current_property_value = this->GetProperties()[rDesignVariable];
            p_local_property->SetValue(rDesignVariable, (current_property_value + delta));

            // Compute stress on GP after disturbance
            this->Calculate(rStressVariable, stress_vector_dist, rCurrentProcessInfo);

            // Compute derivative of stress w.r.t. design variable with finite differences
            noalias(stress_vector_dist)  -= stress_vector_undist;
            stress_vector_dist  /= delta;

            for(size_t j = 0; j < num_gps; j++)
                rOutput(0, j) = stress_vector_dist[j];

            // Give element original properties back
            this->SetProperties(p_global_properties);
        }
        else
            rOutput.clear();

    KRATOS_CATCH("")
}

void ShellThinAdjointElement3D3N::CalculateStressDesignVariableDerivative(const Variable<array_1d<double,3>>& rDesignVariable,
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
        const int number_of_nodes = GetGeometry().PointsNumber();
        const unsigned int dimension = rCurrentProcessInfo.GetValue(DOMAIN_SIZE);

        const SizeType num_gps = GetNumberOfGPs();
        rOutput.resize(dimension * number_of_nodes, num_gps);

        // Compute stress on GP before disturbance
        this->Calculate(rStressVariable, stress_vector_undist, rCurrentProcessInfo);

        int index = 0;
        //TODO: look that this works also for parallel computing
        for(auto& node_i : this->GetGeometry())
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

// =====================================================================================
//
// Class ShellThinAdjointElement3D3N - Serialization
//
// =====================================================================================

void ShellThinAdjointElement3D3N::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer,  ShellThinElement3D3N );
}

void ShellThinAdjointElement3D3N::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, ShellThinElement3D3N );

}

}



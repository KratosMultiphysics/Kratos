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
    }

    AdjointNonlinearStrainEnergyResponseFunction::~AdjointNonlinearStrainEnergyResponseFunction()
    {
    }

// Calculates the response value increment for one time step during the primal analysis
// It is not tested for problems involving self weight
// TODO changing the variable names to follow the guidline
void AdjointNonlinearStrainEnergyResponseFunction::CalculateResponseIncrement(ModelPart& rModelPart)
{
    KRATOS_TRY;

    ProcessInfo &r_current_process_info = rModelPart.GetProcessInfo();
    double response_increment_value = 0.0;

    // Check if there is adjoint elements, because response function calculation is done for primal analysis
    KRATOS_ERROR_IF( r_current_process_info.Has(IS_ADJOINT) && r_current_process_info[IS_ADJOINT] )
    << "Calculate value for strain energy response is only available when using primal elements" << std::endl;

    // sum all elemental strain energy increment values calculated by trapezoidal rule: E = 0.5 * (f_ext_i - f_ext_i-1) * (u_i - u_i-1)
    Matrix LHS;
    Vector RHS;
    Vector disp;
    Vector external_force;
    Vector external_force_previous_step;
    Vector disp_previous_step;
    Vector disp_increment;
    Vector average_load;

    for (auto& elem_i : rModelPart.Elements())
    {
        elem_i.GetValuesVector(disp,0);
        elem_i.GetValuesVector(disp_previous_step, 1);
        elem_i.CalculateLocalSystem(LHS, RHS, r_current_process_info);

        disp_increment = disp - disp_previous_step;
        external_force = -1.0 * RHS;
        external_force_previous_step = external_force - prod(LHS , disp_increment);
        average_load = 0.5 * (external_force + external_force_previous_step);

        response_increment_value += inner_prod(average_load , disp_increment);
    }

    m_response_value += response_increment_value; 

    KRATOS_CATCH("");
}


// TODO, CalculateValue should calculate the response function at each increment,
// and the total value should be calculated in the python script
double AdjointNonlinearStrainEnergyResponseFunction::CalculateValue(ModelPart& rModelPart)
{
    KRATOS_TRY;

    this->CheckForBodyForces(rModelPart);

    return m_response_value;

    KRATOS_CATCH("");
}

void AdjointNonlinearStrainEnergyResponseFunction::CalculateSensitivityGradient(Element& rAdjointElem,
                                              const Variable<array_1d<double,3>>& rVariable,
                                              const Matrix& rDerivativesMatrix,
                                              Vector& rResponseGradient,
                                              ProcessInfo& rProcessInfo)
{
    KRATOS_TRY

    if (rResponseGradient.size() != 0)
        rResponseGradient.resize(0, false);

    KRATOS_CATCH("");
}


void AdjointNonlinearStrainEnergyResponseFunction::CalculateSensitivityGradient(Element& rAdjointElem,
                                              const Variable<double>& rVariable,
                                              const Matrix& rDerivativesMatrix,
                                              Vector& rResponseGradient,
                                              ProcessInfo& rProcessInfo)
{
    KRATOS_TRY;

    if (rResponseGradient.size() != 0)
        rResponseGradient.resize(0, false);

    KRATOS_CATCH("");
}


// Calculates the derivate of the response function with respect to the design parameters
// Computation of \frac{1}{2} (u^T_i - u^T_i-1)  \cdot ( \frac{\partial f_{ext}_i}{\partial x} + frac{\partial f_{ext}_i-1}{\partial x})
void AdjointNonlinearStrainEnergyResponseFunction::CalculateSensitivityGradient(Condition& rAdjointCondition,
                                                                const Variable<double>& rVariable,
                                                                const Matrix& rDerivativesMatrix,
                                                                Vector& rResponseGradient,
                                                                ProcessInfo& rProcessInfo)
{
    KRATOS_TRY;   
    
    KRATOS_CATCH("");
}

void AdjointNonlinearStrainEnergyResponseFunction::CalculateSensitivityGradient(Condition& rAdjointCondition,
                                              const Variable<array_1d<double,3>>& rVariable,
                                              const Matrix& rDerivativesMatrix,
                                              Vector& rResponseGradient,
                                              ProcessInfo& rProcessInfo)
{
    KRATOS_TRY;

    const SizeType number_of_nodes = rAdjointCondition.GetGeometry().size();
    const SizeType dimension = rAdjointCondition.GetGeometry().WorkingSpaceDimension();
    const SizeType mat_size = number_of_nodes * dimension;

    // Matrix for the load partial derivative
    Matrix partial_derivative_matrix = ZeroMatrix(mat_size, mat_size);

    if (rResponseGradient.size() != mat_size)
        rResponseGradient.resize(mat_size, false);

    // step size is temporarily hard coded, later it should be read from json file
    double step_size = 0.01;

    Vector displacement = ZeroVector(mat_size);
    Vector displacement_previous_step = ZeroVector(mat_size);
	Vector RHS;
    Matrix LHS;
    //std::vector<Vector> partial_derivative;

    // This is temporarily hard coded
     double Delta;
     Delta = step_size;

    // accessing nodal displacements for the previous the current and the previous time steps
    // TODO, think about loads where the element is defined on elements, not nodes
    int i_1 = 0;
    // there is a warning here
    for (auto& node_i : rAdjointCondition.GetGeometry())
    {
        project(displacement, range(i_1 , i_1 + dimension)) = rAdjointCondition.GetGeometry()[node_i].FastGetSolutionStepValue(DISPLACEMENT);
        project(displacement_previous_step , range(i_1 , i_1 + dimension))= rAdjointCondition.GetGeometry()[node_i].FastGetSolutionStepValue(DISPLACEMENT,1);   
    }
        std::cout << "Displacement:" << displacement << "::" << displacement_previous_step << std::endl;

    // this does not return the displacement vector
    // rAdjointCondition.GetValuesVector(displacement_previous_step,0);

    // Testing how to access point load for the previous time steps
    // Vector force_1;
    // Vector force_2;

	// force_1 = rAdjointCondition.GetGeometry()[0].FastGetSolutionStepValue(POINT_LOAD);
    // force_2 = rAdjointCondition.GetGeometry()[0].FastGetSolutionStepValue(POINT_LOAD,1);
    // int size = rAdjointCondition.GetGeometry().size();


	//Semi-analytic computation of partial derivative for the force vector w.r.t. node coordinates at the current time step
	rAdjointCondition.CalculateRightHandSide(RHS, rProcessInfo);
    int i_2 = 0;
    // This returns a vector of partial derivative vectors
    // TODO, use rVariable instead of the nodal degrees of freedom
	for (auto& node_i : rAdjointCondition.GetGeometry())
	{
	Vector perturbed_RHS = Vector(0);

	// Pertubation, gradient analysis and recovery of x
	node_i.X0() += Delta;
	rAdjointCondition.CalculateRightHandSide(perturbed_RHS, rProcessInfo);
    column(partial_derivative_matrix, i_2) = (perturbed_RHS - RHS) / Delta;
	node_i.X0() -= Delta;

	// Reset pertubed vector
	perturbed_RHS = Vector(0);

	// Pertubation, gradient analysis and recovery of y
	node_i.Y0() += Delta;
	rAdjointCondition.CalculateRightHandSide(perturbed_RHS, rProcessInfo);
    column(partial_derivative_matrix, i_2 + 1) = (perturbed_RHS - RHS) / Delta;
	node_i.Y0() -= Delta;

	// Reset pertubed vector
	perturbed_RHS = Vector(0);

	// Pertubation, gradient analysis and recovery of z
	node_i.Z0() += Delta;
	rAdjointCondition.CalculateRightHandSide(perturbed_RHS, rProcessInfo);
    column(partial_derivative_matrix, i_2 + 2) = (perturbed_RHS - RHS) / Delta;
	node_i.Z0() -= Delta;

    i_2 += 3;
	}

    // TODO, calculating 0.5 * delta(u) * (partial(f_i)/partial(x)  + partial(f_i-1)/partial(x))
    rResponseGradient = 0.50 * prod(displacement - displacement_previous_step , m_partial_derivative_0 + partial_derivative_matrix);
    std::cout << "response gradient" << rResponseGradient << std::endl;


    // passing the value of the partial derivative matrix to a member variable
    m_partial_derivative_0 = partial_derivative_matrix;

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

// TODO, Discussion of the implementation of the response gradient
// Test the effieciency of calculatevalue() implementation,
// Implementation of the sensitivities, this should be independent from the response function implementation

} // namespace Kratos
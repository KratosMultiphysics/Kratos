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
        // initializing mConditionsRHS 
        for(auto& cond_i : rModelPart.Conditions())
        {
            SizeType number_of_nodes = cond_i.GetGeometry().size();
            SizeType dimension = cond_i.GetGeometry().WorkingSpaceDimension();
            SizeType vec_size = number_of_nodes * dimension;
            mConditionsRHS.push_back(ZeroVector(vec_size));

            // This saves a vector of pointers to the model condition as a member variable because I can not use const condition in 
            // CalculateGradient()
            mConditions.push_back(rModelPart.pGetCondition(cond_i.Id()));
        }
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



// TODO Mahmoud: change the base function so it doesn't const anymore
// calculates the derivative of the response function with respect to nodal displacements
// \frac{1}{2} (f_{ext}_i - f_{ext}_i-1) + frac{1}{2} (u^T_i - u^T_i-1) \cdot ( \frac{\partial f_{ext}_i}{\partial u} + frac{\partial f_{ext}_i-1}{\partial u})
void AdjointNonlinearStrainEnergyResponseFunction::CalculateGradient(const Condition& rAdjointCondition,
                                   const Matrix& rResidualGradient,
                                   Vector& rResponseGradient,
                                   const ProcessInfo& rProcessInfo)
{
    KRATOS_TRY;

    ProcessInfo process_info = rProcessInfo;

    // cloning the condition data, it does not work
    //auto cond_i = rAdjointCondition.Clone(rAdjointCondition.Id(), rAdjointCondition.GetGeometry());
    //auto cond_i = rAdjointCondition.Create(rAdjointCondition.Id(), rAdjointCondition.GetGeometry(), rAdjointCondition.GetProperties() );
    // Deep copy elemental data does not work on const variable
    //cond_i->Data() = rAdjointCondition.Data();

    const SizeType number_of_nodes = rAdjointCondition.GetGeometry().size();
    const SizeType dimension = rAdjointCondition.GetGeometry().WorkingSpaceDimension();
    const SizeType vec_size = number_of_nodes * dimension;

    Vector displacement = ZeroVector(vec_size);
    Vector displacement_previous_step = ZeroVector(vec_size);

    // getting the displacements for the current and the previous time steps
    int i = 0;
    for (auto& node_i : rAdjointCondition.GetGeometry())
    {
        project(displacement, range(i, dimension)) = rAdjointCondition.GetGeometry()[node_i].FastGetSolutionStepValue(DISPLACEMENT);
        project(displacement_previous_step , range(i , dimension))= rAdjointCondition.GetGeometry()[node_i].FastGetSolutionStepValue(DISPLACEMENT,1);   
        i += 3;
    }
    // perturbation value , hard coded value
    Vector delta = 0.001 * displacement;
    std::cout << "displacement" << displacement << std::endl;
    
    // getting RHS for previous time step
    Vector RHS_previous_step = mConditionsRHS[rAdjointCondition.Id() - 1];

    // Getting RHS for the current time step
    Vector RHS;

    // TODO Mahmoud: change this later to use the condition directory
    mConditions[rAdjointCondition.Id() - 1]->CalculateRightHandSide(RHS, process_info);
    

    //calculation of the perturbed RHS 
    Matrix partial_derivative_matrix = ZeroMatrix(vec_size, vec_size);
    int i_2 = 0;
	for (auto& node_i : rAdjointCondition.GetGeometry())
	{
	Vector perturbed_RHS = Vector(0);
	// Pertubation, gradient analysis and recovery of x
	node_i.X() += delta[i_2]; 
    if(delta[i_2] != 0)
    {
        mConditions[rAdjointCondition.Id() - 1]->CalculateRightHandSide(perturbed_RHS, process_info);
        row(partial_derivative_matrix, i_2) = (perturbed_RHS - RHS) / delta[i_2];
	    node_i.X() -= delta[i_2];
    }
	
	// Reset pertubed vector
	perturbed_RHS = Vector(0);
	// Pertubation, gradient analysis and recovery of y
    if(delta[i_2 + 1] != 0) 
    {
        node_i.Y() += delta[i_2 + 1];
	    mConditions[rAdjointCondition.Id() - 1]->CalculateRightHandSide(perturbed_RHS, process_info);
        std::cout << "RHS" << RHS << "perturbed RHS" << perturbed_RHS << std::endl;
        row(partial_derivative_matrix, i_2 + 1) = (perturbed_RHS - RHS) / 0.001;
	    node_i.Y() -= delta[i_2 + 1];
    }

	// Reset pertubed vector
	perturbed_RHS = Vector(0);
	// Pertubation, gradient analysis and recovery of z
    if(delta[i_2 + 2] != 0) 
    {
        node_i.Z() += delta[i_2 + 2];
        mConditions[rAdjointCondition.Id() - 1]->CalculateRightHandSide(perturbed_RHS, process_info);
        row(partial_derivative_matrix, i_2 + 2) = (perturbed_RHS - RHS) / delta[i_2 + 2];
        node_i.Z() -= delta[i_2 + 2];
    }
    i_2 += 3;
	}

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
    ProcessInfo process_info = rProcessInfo;

    const SizeType number_of_nodes = rAdjointCondition.GetGeometry().size();
    const SizeType dimension = rAdjointCondition.GetGeometry().WorkingSpaceDimension();
    const SizeType mat_size = number_of_nodes * dimension;

    // Matrix for the load partial derivative
    Matrix partial_derivative_matrix = ZeroMatrix(mat_size, mat_size);

    if(m_partial_derivative_0.size1() == 0)
        m_partial_derivative_0 = ZeroMatrix(mat_size, mat_size);

    if (rSensitivityGradient.size() != mat_size)
        rSensitivityGradient.resize(mat_size, false);

    // step size is temporarily hard coded, later it should be read from json file
    double Delta = rAdjointCondition.GetValue(PERTURBATION_SIZE);

    Vector displacement = ZeroVector(mat_size);
    Vector displacement_previous_step = ZeroVector(mat_size);
	Vector RHS;
    Matrix LHS;
    //std::vector<Vector> partial_derivative

    // accessing nodal displacements for the previous the current and the previous time steps
    // TODO, think about loads where the condition is defined on elements, not nodes
    int i_1 = 0;
    // there is a warning here
    for (auto& node_i : rAdjointCondition.GetGeometry())
    {
        project(displacement, range(i_1 , i_1 + dimension)) = rAdjointCondition.GetGeometry()[node_i].FastGetSolutionStepValue(DISPLACEMENT);
        project(displacement_previous_step , range(i_1 , i_1 + dimension))= rAdjointCondition.GetGeometry()[node_i].FastGetSolutionStepValue(DISPLACEMENT,1); 
        i_1 += 3; 
    }
        //std::cout << "Displacement:" << displacement << "::" << displacement_previous_step << std::endl;

	//Semi-analytic computation of partial derivative for the force vector w.r.t. node coordinates at the current time step
	rAdjointCondition.CalculateRightHandSide(RHS, process_info);
    // This returns a vector of partial derivative vectors
    // TODO, use rVariable instead of the nodal degrees of freedom
    int i_2 = 0;
	for (auto& node_i : rAdjointCondition.GetGeometry())
	{
	Vector perturbed_RHS = Vector(0);

	// Pertubation, gradient analysis and recovery of x
	node_i.X() += Delta;
	rAdjointCondition.CalculateRightHandSide(perturbed_RHS, process_info);
    row(partial_derivative_matrix, i_2) = (perturbed_RHS - RHS) / Delta;
	node_i.X() -= Delta;

	// Reset pertubed vector
	perturbed_RHS = Vector(0);

	// Pertubation, gradient analysis and recovery of y
	node_i.Y() += Delta;
	rAdjointCondition.CalculateRightHandSide(perturbed_RHS, process_info);
    row(partial_derivative_matrix, i_2 + 1) = (perturbed_RHS - RHS) / Delta;
	node_i.Y() -= Delta;

	// Reset pertubed vector
	perturbed_RHS = Vector(0);

	// Pertubation, gradient analysis and recovery of z
	node_i.Z() += Delta;
	rAdjointCondition.CalculateRightHandSide(perturbed_RHS, process_info);
    row(partial_derivative_matrix, i_2 + 2) = (perturbed_RHS - RHS) / Delta;
	node_i.Z() -= Delta;

    i_2 += 3;
	}

    // TODO, calculating 0.5 * delta(u) * (partial(f_i)/partial(x)  + partial(f_i-1)/partial(x))
    rSensitivityGradient = 0.50 * prod(displacement - displacement_previous_step , m_partial_derivative_0 + partial_derivative_matrix);
   // std::cout << "response gradient" << rSensitivityGradient << std::endl;

    // passing the value of the partial derivative matrix to a member variable
    //TODO, this way is wrong theoratically, the RHS should be stored, not its derivative
    m_partial_derivative_0 = partial_derivative_matrix;
    
    // storing the RHS for each condition
    mConditionsRHS[rAdjointCondition.Id() - 1] = RHS;

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
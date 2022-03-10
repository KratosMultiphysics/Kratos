// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:     BSD License
//	             license: structural_mechanics_application/license.txt
//
//  Main authors:    Baumgaertner Daniel, https://github.com/dbaumgaertner
//                   Geiser Armin, https://github.com/armingeiser
//

#ifndef STRAIN_ENERGY_RESPONSE_FUNCTION_UTILITY_H
#define STRAIN_ENERGY_RESPONSE_FUNCTION_UTILITY_H

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------
#include <iostream>
#include <string>

// ------------------------------------------------------------------------------
// External includes
// ------------------------------------------------------------------------------

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "includes/define.h"
#include "includes/kratos_parameters.h"
#include "includes/element.h"
#include "utilities/integration_utilities.h"
#include "utilities/geometry_utilities.h"
#include "includes/model_part.h"
#include "utilities/variable_utils.h"
#include "processes/find_nodal_neighbours_process.h"
#include "finite_difference_utility.h"

// ==============================================================================

namespace Kratos
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/// Short class definition.
/** Detail class definition.

 */

class StrainEnergyResponseFunctionUtility
{
public:
	///@name Type Definitions
	///@{

	typedef array_1d<double, 3> array_3d;

	/// Pointer definition of StrainEnergyResponseFunctionUtility
	KRATOS_CLASS_POINTER_DEFINITION(StrainEnergyResponseFunctionUtility);

	///@}
	///@name Life Cycle
	///@{

	/// Default constructor.
	StrainEnergyResponseFunctionUtility(ModelPart& model_part, Parameters responseSettings)
	: mrModelPart(model_part)
	{
		// Set gradient mode
		std::string gradient_mode = responseSettings["gradient_settings"]["gradient_mode"].GetString();

		if (gradient_mode.compare("semi_analytic") == 0)
		{
			double delta = responseSettings["gradient_settings"]["step_size"].GetDouble();
			mDelta = delta;
		}
		else
			KRATOS_ERROR << "Specified gradient_mode '" << gradient_mode << "' not recognized. The only option is: semi_analytic" << std::endl;

	}

	/// Destructor.
	virtual ~StrainEnergyResponseFunctionUtility()
	{
	}

	///@}
	///@name Operators
	///@{

	///@}
	///@name Operations
	///@{

	// ==============================================================================
	void Initialize()
	{

		for (auto& elem_i : mrModelPart.Elements())
		{
			Properties& elem_i_prop = elem_i.GetProperties();
			Properties::Pointer model_part_new_prop = mrModelPart.CreateNewProperties(mrModelPart.NumberOfProperties()+1);
			*model_part_new_prop = elem_i_prop;
			elem_i.SetProperties(model_part_new_prop);
		}

	}

	// --------------------------------------------------------------------------
	double CalculateValue()
	{
		KRATOS_TRY;

		const ProcessInfo &CurrentProcessInfo = mrModelPart.GetProcessInfo();
		double strain_energy = 0.0;

		// Sum all elemental strain energy values calculated as: W_e = u_e^T K_e u_e
		for (auto& elem_i : mrModelPart.Elements())
		{
			const bool element_is_active = elem_i.IsDefined(ACTIVE) ? elem_i.Is(ACTIVE) : true;
			if(element_is_active)
			{
				Matrix LHS;
				Vector RHS;
				Vector u;

				// Get state solution relevant for energy calculation
				const auto& rConstElemRef = elem_i;
				rConstElemRef.GetValuesVector(u,0);

				elem_i.CalculateLocalSystem(LHS,RHS,CurrentProcessInfo);

				// Compute strain energy
				strain_energy += 0.5 * inner_prod(u,prod(LHS,u));
			}
		}

		return strain_energy;

		KRATOS_CATCH("");
	}

	// --------------------------------------------------------------------------
	void CalculateGradient()
	{
		KRATOS_TRY;

		// Formula computed in general notation:
		// \frac{dF}{dx} = \frac{\partial F}{\partial x}
		//                 - \lambda^T \frac{\partial R_S}{\partial x}

		// Adjoint variable:
		// \lambda       = \frac{1}{2} u^T

		// Correspondingly in specific notation for given response function
		// \frac{dF}{dx} = \frac{1}{2} u^T \cdot \frac{\partial f_{ext}}{\partial x}
		//                 + \frac{1}{2} u^T \cdot ( \frac{\partial f_{ext}}{\partial x} - \frac{\partial K}{\partial x} )

		// First gradients are initialized
		VariableUtils().SetHistoricalVariableToZero(SHAPE_SENSITIVITY, mrModelPart.Nodes());
		VariableUtils().SetHistoricalVariableToZero(THICKNESS_SENSITIVITY, mrModelPart.Nodes());
		VariableUtils().SetHistoricalVariableToZero(YOUNG_MODULUS_SENSITIVITY, mrModelPart.Nodes());


		// Gradient calculation
		// 1st step: Calculate partial derivative of response function w.r.t. node coordinates
		// 2nd step: Calculate adjoint field
		// 3rd step: Calculate partial derivative of state equation w.r.t. node coordinates and multiply with adjoint field

		// Semi analytic sensitivities
		CalculateResponseDerivativePartByFiniteDifferencing();
		CalculateAdjointField();
		CalculateStateDerivativePartByFiniteDifferencing();

		KRATOS_CATCH("");
	}

	// ==============================================================================

	///@}
	///@name Access
	///@{

	///@}
	///@name Inquiry
	///@{

	///@}
	///@name Input and output
	///@{

	/// Turn back information as a string.
	std::string Info() const
	{
		return "StrainEnergyResponseFunctionUtility";
	}

	/// Print information about this object.
	virtual void PrintInfo(std::ostream &rOStream) const
	{
		rOStream << "StrainEnergyResponseFunctionUtility";
	}

	/// Print object's data.
	virtual void PrintData(std::ostream &rOStream) const
	{
	}

	///@}
	///@name Friends
	///@{

	///@}

protected:
	///@name Protected static Member Variables
	///@{

	///@}
	///@name Protected member Variables
	///@{

	///@}
	///@name Protected Operators
	///@{

	///@}
	///@name Protected Operations
	///@{

	// ==============================================================================
	void CalculateAdjointField()
	{
		KRATOS_TRY;

		// Adjoint field may be directly obtained from state solution

		KRATOS_CATCH("");
	}

	// --------------------------------------------------------------------------
	void CalculateStateDerivativePartByFiniteDifferencing()
	{
		KRATOS_TRY;

		// Working variables
		const ProcessInfo &CurrentProcessInfo = mrModelPart.GetProcessInfo();

		// Computation of: \frac{1}{2} u^T \cdot ( - \frac{\partial K}{\partial x} )
		for (auto& elem_i : mrModelPart.Elements())
		{
			const bool element_is_active = elem_i.IsDefined(ACTIVE) ? elem_i.Is(ACTIVE) : true;
			if(element_is_active)
			{
				Vector u;
				Vector lambda;
				Vector RHS;

				// Get state solution
				const auto& rConstElemRef = elem_i;
				rConstElemRef.GetValuesVector(u,0);

				// Get adjoint variables (Corresponds to 1/2*u)
				lambda = 0.5*u;

				// Semi-analytic computation of partial derivative of state equation w.r.t. node coordinates
				elem_i.CalculateRightHandSide(RHS, CurrentProcessInfo);
				for (auto& node_i : elem_i.GetGeometry())
				{
					array_3d gradient_contribution(3, 0.0);
					Vector derived_RHS = Vector(0);

					// x-direction
					FiniteDifferenceUtility::CalculateRightHandSideDerivative(elem_i, RHS, SHAPE_SENSITIVITY_X, node_i, mDelta, derived_RHS, CurrentProcessInfo);
					gradient_contribution[0] = inner_prod(lambda, derived_RHS);

					// y-direction
					FiniteDifferenceUtility::CalculateRightHandSideDerivative(elem_i, RHS, SHAPE_SENSITIVITY_Y, node_i, mDelta, derived_RHS, CurrentProcessInfo);
					gradient_contribution[1] = inner_prod(lambda, derived_RHS);

					// z-direction
					FiniteDifferenceUtility::CalculateRightHandSideDerivative(elem_i, RHS, SHAPE_SENSITIVITY_Z, node_i, mDelta, derived_RHS, CurrentProcessInfo);
					gradient_contribution[2] = inner_prod(lambda, derived_RHS);

					// Assemble sensitivity to node
					noalias(node_i.FastGetSolutionStepValue(SHAPE_SENSITIVITY)) += gradient_contribution;
				}
				//now compute YOUNG_MODULUS_SENSITIVITY
				Properties& elem_i_prop = elem_i.GetProperties();
				if(elem_i_prop.Has(YOUNG_MODULUS))
				{
					auto original_young_modul = elem_i_prop.GetValue(YOUNG_MODULUS);
					Vector ANAL_RHS_SENS;
					// elem_i_prop.SetValue(YOUNG_MODULUS,1);
					elem_i.CalculateRightHandSide(ANAL_RHS_SENS, CurrentProcessInfo);
					double anal_sens = inner_prod(lambda, ANAL_RHS_SENS);
					const auto& r_geom = elem_i.GetGeometry();	
					const auto& integration_method = r_geom.GetDefaultIntegrationMethod();
					const auto& integration_points = r_geom.IntegrationPoints(integration_method);
					const unsigned int NumGauss = integration_points.size();
					Vector GaussPtsJDet = ZeroVector(NumGauss);
					r_geom.DeterminantOfJacobian(GaussPtsJDet, integration_method);
					const auto& Ncontainer = r_geom.ShapeFunctionsValues(integration_method); 
					for (auto& node_i : r_geom){
						auto node_weight = node_i.GetValue(NUMBER_OF_NEIGHBOUR_ELEMENTS);
						node_i.FastGetSolutionStepValue(YOUNG_MODULUS_SENSITIVITY) += anal_sens/node_weight ;
					}					
					// for(std::size_t i_point = 0; i_point<integration_points.size(); ++i_point)
					// {
					// 	const double IntToReferenceWeight = integration_points[i_point].Weight() * GaussPtsJDet[i_point];
					// 	const auto& rN = row(Ncontainer,i_point);
					// 	int node_index = 0;
					// 	for (auto& node_i : r_geom){
					// 		node_i.FastGetSolutionStepValue(YOUNG_MODULUS_SENSITIVITY) += anal_sens * IntToReferenceWeight * rN[node_index];
					// 		node_index++;
					// 	}						
					// }		
				}
				if(elem_i_prop.Has(THICKNESS))
				{
					auto original_thickness = elem_i_prop.GetValue(THICKNESS);
					Vector ANAL_RHS_SENS;
					elem_i.CalculateRightHandSide(ANAL_RHS_SENS, CurrentProcessInfo);
					double anal_sens = inner_prod(lambda, ANAL_RHS_SENS);
					const auto& r_geom = elem_i.GetGeometry();	
					const auto& integration_method = r_geom.GetDefaultIntegrationMethod();
					const auto& integration_points = r_geom.IntegrationPoints(integration_method);
					const unsigned int NumGauss = integration_points.size();
					Vector GaussPtsJDet = ZeroVector(NumGauss);
					r_geom.DeterminantOfJacobian(GaussPtsJDet, integration_method);
					const auto& Ncontainer = r_geom.ShapeFunctionsValues(integration_method); 
					for (auto& node_i : r_geom){
						auto node_weight = node_i.GetValue(NUMBER_OF_NEIGHBOUR_ELEMENTS);
						node_i.FastGetSolutionStepValue(THICKNESS_SENSITIVITY) += anal_sens/node_weight ;
					}					
					// for(std::size_t i_point = 0; i_point<integration_points.size(); ++i_point)
					// {
					// 	const double IntToReferenceWeight = integration_points[i_point].Weight() * GaussPtsJDet[i_point];
					// 	const auto& rN = row(Ncontainer,i_point);
					// 	int node_index = 0;
					// 	for (auto& node_i : r_geom){
					// 		node_i.FastGetSolutionStepValue(THICKNESS_SENSITIVITY) += anal_sens * IntToReferenceWeight * rN[node_index];
					// 		node_index++;
					// 	}						
					// }		
				}				
				
				
				
				
				
				// auto thickness


				// double pert = 1e-6*young_modul;

				// //forward_pert
				// elem_i_prop.SetValue(YOUNG_MODULUS,young_modul+pert);
				// Vector RHS_F;
				// elem_i.CalculateRightHandSide(RHS_F, CurrentProcessInfo);

				// //backward_pert
				// elem_i_prop.SetValue(YOUNG_MODULUS,young_modul-pert);
				// Vector RHS_B;
				// elem_i.CalculateRightHandSide(RHS_B, CurrentProcessInfo);

				// Vector FD_RHS_SENS=(RHS_F-RHS_B)/(2*pert);
				// double fd_sens = inner_prod(lambda, FD_RHS_SENS);



				// elem_i_prop.SetValue(YOUNG_MODULUS,young_modul);

				// std::cout<<"fd_sens : "<<fd_sens<<", anal_sens : "<<anal_sens<<std::endl;

				// Properties& elem_i_prop = elem_i.GetProperties();
				// std::cout<<"elem_i_prop[YOUNG_MODULUS] : "<<elem_i_prop.GetValue(YOUNG_MODULUS)<<std::endl;
				// elem_i_prop.SetValue(YOUNG_MODULUS,1.0);
			}
		}

		// Computation of \frac{1}{2} u^T \cdot ( \frac{\partial f_{ext}}{\partial x} )
		for (auto& cond_i : mrModelPart.Conditions())
		{
			//detect if the condition is active or not. If the user did not make any choice the element
			//is active by default
			const bool condition_is_active = cond_i.IsDefined(ACTIVE) ? cond_i.Is(ACTIVE) : true;
			if (condition_is_active)
			{
				Vector u;
				Vector lambda;
				Vector RHS;

				// Get adjoint variables (Corresponds to 1/2*u)
				const auto& rConstCondRef = cond_i;
				rConstCondRef.GetValuesVector(u,0);
				lambda = 0.5*u;

				// Semi-analytic computation of partial derivative of force vector w.r.t. node coordinates
				cond_i.CalculateRightHandSide(RHS, CurrentProcessInfo);
				for (auto& node_i : cond_i.GetGeometry())
				{
					array_3d gradient_contribution(3, 0.0);
					Vector perturbed_RHS = Vector(0);

					// Pertubation, gradient analysis and recovery of x
					node_i.X0() += mDelta;
					cond_i.CalculateRightHandSide(perturbed_RHS, CurrentProcessInfo);
					gradient_contribution[0] = inner_prod(lambda, (perturbed_RHS - RHS) / mDelta);
					node_i.X0() -= mDelta;

					// Reset pertubed vector
					perturbed_RHS = Vector(0);

					// Pertubation, gradient analysis and recovery of y
					node_i.Y0() += mDelta;
					cond_i.CalculateRightHandSide(perturbed_RHS, CurrentProcessInfo);
					gradient_contribution[1] = inner_prod(lambda, (perturbed_RHS - RHS) / mDelta);
					node_i.Y0() -= mDelta;

					// Reset pertubed vector
					perturbed_RHS = Vector(0);

					// Pertubation, gradient analysis and recovery of z
					node_i.Z0() += mDelta;
					cond_i.CalculateRightHandSide(perturbed_RHS, CurrentProcessInfo);
					gradient_contribution[2] = inner_prod(lambda, (perturbed_RHS - RHS) / mDelta);
					node_i.Z0() -= mDelta;

					// Assemble shape gradient to node
					noalias(node_i.FastGetSolutionStepValue(SHAPE_SENSITIVITY)) += gradient_contribution;
				}
			}
		}

		KRATOS_CATCH("");
	}

	// --------------------------------------------------------------------------
	void CalculateResponseDerivativePartByFiniteDifferencing()
	{
		KRATOS_TRY;

		// Working variables
		const ProcessInfo &CurrentProcessInfo = mrModelPart.GetProcessInfo();

		// Computation of \frac{1}{2} u^T \cdot ( \frac{\partial f_{ext}}{\partial x} )
		for (auto& cond_i : mrModelPart.Conditions())
		{
			//detect if the condition is active or not. If the user did not make any choice the element
			//is active by default
			const bool condition_is_active = cond_i.IsDefined(ACTIVE) ? cond_i.Is(ACTIVE) : true;
			if (condition_is_active)
			{
				Vector u;
				Vector RHS;

				// Get state solution
				const auto& rConstCondRef = cond_i;
				rConstCondRef.GetValuesVector(u,0);

				// Perform finite differencing of RHS vector
				cond_i.CalculateRightHandSide(RHS, CurrentProcessInfo);
				for (auto& node_i : cond_i.GetGeometry())
				{
					array_3d gradient_contribution(3, 0.0);
					Vector perturbed_RHS = Vector(0);

					// Pertubation, gradient analysis and recovery of x
					node_i.X0() += mDelta;
					cond_i.CalculateRightHandSide(perturbed_RHS, CurrentProcessInfo);
					gradient_contribution[0] = inner_prod(0.5*u, (perturbed_RHS - RHS) / mDelta);
					node_i.X0() -= mDelta;

					// Reset pertubed vector
					perturbed_RHS = Vector(0);

					// Pertubation, gradient analysis and recovery of y
					node_i.Y0() += mDelta;
					cond_i.CalculateRightHandSide(perturbed_RHS, CurrentProcessInfo);
					gradient_contribution[1] = inner_prod(0.5*u, (perturbed_RHS - RHS) / mDelta);
					node_i.Y0() -= mDelta;

					// Reset pertubed vector
					perturbed_RHS = Vector(0);

					// Pertubation, gradient analysis and recovery of z
					node_i.Z0() += mDelta;
					cond_i.CalculateRightHandSide(perturbed_RHS, CurrentProcessInfo);
					gradient_contribution[2] = inner_prod(0.5*u, (perturbed_RHS - RHS) / mDelta);
					node_i.Z0() -= mDelta;

					// Assemble shape gradient to node
					noalias(node_i.FastGetSolutionStepValue(SHAPE_SENSITIVITY)) += gradient_contribution;
				}
			}
		}
		KRATOS_CATCH("");
	}

	// ==============================================================================

	///@}
	///@name Protected  Access
	///@{

	///@}
	///@name Protected Inquiry
	///@{

	///@}
	///@name Protected LifeCycle
	///@{

	///@}

private:
	///@name Static Member Variables
	///@{

	///@}
	///@name Member Variables
	///@{

	ModelPart &mrModelPart;
	double mDelta;

	///@}
///@name Private Operators
	///@{

	///@}
	///@name Private Operations
	///@{

	///@}
	///@name Private  Access
	///@{

	///@}
	///@name Private Inquiry
	///@{

	///@}
	///@name Un accessible methods
	///@{

	/// Assignment operator.
	//      StrainEnergyResponseFunctionUtility& operator=(StrainEnergyResponseFunctionUtility const& rOther);

	/// Copy constructor.
	//      StrainEnergyResponseFunctionUtility(StrainEnergyResponseFunctionUtility const& rOther);

	///@}

}; // Class StrainEnergyResponseFunctionUtility

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

} // namespace Kratos.

#endif // STRAIN_ENERGY_RESPONSE_FUNCTION_UTILITY_H
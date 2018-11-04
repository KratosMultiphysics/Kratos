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
#include "includes/model_part.h"
#include "utilities/variable_utils.h"
#include "processes/find_nodal_neighbours_process.h"
#include "element_finite_difference_utility.h"

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
		std::string gradient_mode = responseSettings["gradient_mode"].GetString();

		if (gradient_mode.compare("semi_analytic") == 0)
		{
			double delta = responseSettings["step_size"].GetDouble();
			mDelta = delta;
		}
		else
			KRATOS_ERROR << "Specified gradient_mode '" << gradient_mode << "' not recognized. The only option is: semi_analytic" << std::endl;

		mConsiderDiscretization =  responseSettings["consider_discretization"].GetBool();
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
	}

	// --------------------------------------------------------------------------
	double CalculateValue()
	{
		KRATOS_TRY;

		ProcessInfo &CurrentProcessInfo = mrModelPart.GetProcessInfo();
		double strain_energy = 0.0;

		// Sum all elemental strain energy values calculated as: W_e = u_e^T K_e u_e
		for (auto& elem_i : mrModelPart.Elements())
		{
			Matrix LHS;
			Vector RHS;
			Vector u;

			// Get state solution relevant for energy calculation
			elem_i.GetValuesVector(u,0);

			elem_i.CalculateLocalSystem(LHS,RHS,CurrentProcessInfo);

			// Compute strain energy
			strain_energy += 0.5 * inner_prod(u,prod(LHS,u));
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
		VariableUtils().SetToZero_VectorVar(SHAPE_SENSITIVITY, mrModelPart.Nodes());


		// Gradient calculation
		// 1st step: Calculate partial derivative of response function w.r.t. node coordinates
		// 2nd step: Calculate adjoint field
		// 3rd step: Calculate partial derivative of state equation w.r.t. node coordinates and multiply with adjoint field

		// Semi analytic sensitivities
		CalculateResponseDerivativePartByFiniteDifferencing();
		CalculateAdjointField();
		CalculateStateDerivativePartByFiniteDifferencing();

		if (mConsiderDiscretization)
			this->ConsiderDiscretization();

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
		ProcessInfo &CurrentProcessInfo = mrModelPart.GetProcessInfo();

		// Computation of: \frac{1}{2} u^T \cdot ( - \frac{\partial K}{\partial x} )
		for (auto& elem_i : mrModelPart.Elements())
		{
			Vector u;
			Vector lambda;
			Vector RHS;

			// Get state solution
			elem_i.GetValuesVector(u,0);

			// Get adjoint variables (Corresponds to 1/2*u)
			lambda = 0.5*u;

			// Semi-analytic computation of partial derivative of state equation w.r.t. node coordinates
			elem_i.CalculateRightHandSide(RHS, CurrentProcessInfo);
			for (auto& node_i : elem_i.GetGeometry())
			{
				array_3d gradient_contribution(3, 0.0);
				Vector derived_RHS = Vector(0);

				// x-direction
				ElementFiniteDifferenceUtility::CalculateRightHandSideDerivative(elem_i, SHAPE_X, node_i, mDelta, derived_RHS, CurrentProcessInfo);
				gradient_contribution[0] = inner_prod(lambda, derived_RHS);

                // y-direction
				ElementFiniteDifferenceUtility::CalculateRightHandSideDerivative(elem_i, SHAPE_Y, node_i, mDelta, derived_RHS, CurrentProcessInfo);
				gradient_contribution[1] = inner_prod(lambda, derived_RHS);

                // z-direction
				ElementFiniteDifferenceUtility::CalculateRightHandSideDerivative(elem_i, SHAPE_Z, node_i, mDelta, derived_RHS, CurrentProcessInfo);
				gradient_contribution[2] = inner_prod(lambda, derived_RHS);

				// Assemble sensitivity to node
				noalias(node_i.FastGetSolutionStepValue(SHAPE_SENSITIVITY)) += gradient_contribution;
			}
		}

		// Computation of \frac{1}{2} u^T \cdot ( \frac{\partial f_{ext}}{\partial x} )
		for (auto& cond_i : mrModelPart.Conditions())
		{
			//detect if the condition is active or not. If the user did not make any choice the element
			//is active by default
			bool condition_is_active = true;
			if (cond_i.IsDefined(ACTIVE))
				condition_is_active = cond_i.Is(ACTIVE);

			if (condition_is_active)
			{
				Vector u;
				Vector lambda;
				Vector RHS;

				// Get adjoint variables (Corresponds to 1/2*u)
				cond_i.GetValuesVector(u,0);
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
		ProcessInfo &CurrentProcessInfo = mrModelPart.GetProcessInfo();

		// Computation of \frac{1}{2} u^T \cdot ( \frac{\partial f_{ext}}{\partial x} )
		for (auto& cond_i : mrModelPart.Conditions())
		{
			//detect if the condition is active or not. If the user did not make any choice the element
			//is active by default
			bool condition_is_active = true;
			if (cond_i.IsDefined(ACTIVE))
				condition_is_active = cond_i.Is(ACTIVE);

			if (condition_is_active)
			{
				Vector u;
				Vector RHS;

				// Get state solution
				cond_i.GetValuesVector(u,0);

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

	// --------------------------------------------------------------------------
  	virtual void ConsiderDiscretization()
	{


		// Start process to identify element neighbors for every node
		FindNodalNeighboursProcess neigbhorFinder = FindNodalNeighboursProcess(mrModelPart, 10, 10);
		neigbhorFinder.Execute();

		std::cout<< "> Considering discretization size!" << std::endl;
		for(auto& node_i : mrModelPart.Nodes())
		{
			WeakPointerVector<Element >& ng_elem = node_i.GetValue(NEIGHBOUR_ELEMENTS);

			double scaling_factor = 0.0;
			for(std::size_t i = 0; i < ng_elem.size(); i++)
			{
				Kratos::Element& ng_elem_i = ng_elem[i];
				Element::GeometryType& element_geometry = ng_elem_i.GetGeometry();

				scaling_factor += element_geometry.DomainSize();
			}

			node_i.FastGetSolutionStepValue(SHAPE_SENSITIVITY) /= scaling_factor;
		}
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
	bool mConsiderDiscretization;

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

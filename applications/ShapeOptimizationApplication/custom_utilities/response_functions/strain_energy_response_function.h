// ==============================================================================
//  KratosShapeOptimizationApplication
//
//  License:         BSD License
//                   license: ShapeOptimizationApplication/license.txt
//
//  Main authors:    Baumgaertner Daniel, https://github.com/dbaumgaertner
//                   Geiser Armin, https://github.com/armingeiser
//
// ==============================================================================

#ifndef STRAIN_ENERGY_RESPONSE_FUNCTION_H
#define STRAIN_ENERGY_RESPONSE_FUNCTION_H

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------
#include <iostream>
#include <string>
#include <algorithm>

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "../../kratos/includes/define.h"
#include "../../kratos/processes/process.h"
#include "../../kratos/includes/node.h"
#include "../../kratos/includes/element.h"
#include "../../kratos/includes/model_part.h"
#include "../../kratos/includes/kratos_flags.h"
#include "processes/find_nodal_neighbours_process.h"
#include "response_function.h"

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

class StrainEnergyResponseFunction : ResponseFunction
{
public:
	///@name Type Definitions
	///@{

	typedef array_1d<double, 3> array_3d;

	/// Pointer definition of StrainEnergyResponseFunction
	KRATOS_CLASS_POINTER_DEFINITION(StrainEnergyResponseFunction);

	///@}
	///@name Life Cycle
	///@{

	/// Default constructor.
	StrainEnergyResponseFunction(ModelPart& model_part, Parameters responseSettings)
	: mr_model_part(model_part)
	{
		// Set gradient mode
		std::string gradientMode = responseSettings["gradient_mode"].GetString();

		// Mode 1: semi-analytic sensitivities
		if (gradientMode.compare("semi_analytic") == 0)
		{
			mGradientMode = 1;
			double delta = responseSettings["step_size"].GetDouble();
			mDelta = delta;
		}
		// Throw error message in case of wrong specification
		else
			KRATOS_ERROR << "Specified gradient_mode not recognized. Options are: semi_analytic. Specified gradient_mode: " << gradientMode << std::endl;

		mConsiderDiscretization =  responseSettings["consider_discretization"].GetBool();

		// Initialize member variables to NULL
		m_strain_energy = 0.0;
	}

	/// Destructor.
	virtual ~StrainEnergyResponseFunction()
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
		//not needed because only semi-analytical sensitivity analysis is implemented yet
	}

	// --------------------------------------------------------------------------
	void CalculateValue()
	{
		KRATOS_TRY;

		ProcessInfo &CurrentProcessInfo = mr_model_part.GetProcessInfo();
		m_strain_energy = 0.0;

		// Sum all elemental strain energy values calculated as: W_e = u_e^T K_e u_e
		for (ModelPart::ElementIterator elem_i = mr_model_part.ElementsBegin(); elem_i != mr_model_part.ElementsEnd(); ++elem_i)
		{
			Matrix LHS;
			Vector RHS;
			Vector u;

			// Get state solution relevant for energy calculation
			elem_i->GetValuesVector(u,0);

			elem_i->CalculateLocalSystem(LHS,RHS,CurrentProcessInfo);

			// Compute strain energy
			m_strain_energy += 0.5 * inner_prod(u,prod(LHS,u));
		}

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
		array_3d zeros_array(3, 0.0);
		for (ModelPart::NodeIterator node_i = mr_model_part.NodesBegin(); node_i != mr_model_part.NodesEnd(); ++node_i)
			noalias(node_i->FastGetSolutionStepValue(STRAIN_ENERGY_SHAPE_GRADIENT)) = zeros_array;

		// Gradient calculation separated in analytic and semi-analytic approaches
		// In each approach, three following steps are performed:
		// 1st step: Calculate partial derivative of response function w.r.t. node coordinates
		// 2nd step: Calculate adjoint field
		// 3rd step: Calculate partial derivative of state equation w.r.t. node coordinates and multiply with adjoint field

		switch (mGradientMode)
		{
		// Semi analytic sensitivities
		case 1:
		{
			CalculateResponseDerivativePartByFiniteDifferencing();
			CalculateAdjointField();
			CalculateStateDerivativePartByFiniteDifferencing();
			break;
		}
		}

		if (mConsiderDiscretization)
			this->ConsiderDiscretization();

		KRATOS_CATCH("");
	}

	// --------------------------------------------------------------------------
	double GetValue()
	{
		KRATOS_TRY;

		return m_strain_energy;

		KRATOS_CATCH("");
	}

	// --------------------------------------------------------------------------
	pybind11::dict GetGradient()
	{
		KRATOS_TRY;

		// Dictionary to store all sensitivities along with Ids of corresponding nodes
		pybind11::dict dFdX;

		// Fill dictionary with gradient information
		for (ModelPart::NodeIterator node_i = mr_model_part.NodesBegin(); node_i != mr_model_part.NodesEnd(); ++node_i)
			dFdX[pybind11::cast(node_i->Id())] = node_i->FastGetSolutionStepValue(STRAIN_ENERGY_SHAPE_GRADIENT);

		return dFdX;

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
	virtual std::string Info() const
	{
		return "StrainEnergyResponseFunction";
	}

	/// Print information about this object.
	virtual void PrintInfo(std::ostream &rOStream) const
	{
		rOStream << "StrainEnergyResponseFunction";
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
		ProcessInfo &CurrentProcessInfo = mr_model_part.GetProcessInfo();

		// Computation of: \frac{1}{2} u^T \cdot ( - \frac{\partial K}{\partial x} )
		for (ModelPart::ElementIterator elem_i = mr_model_part.ElementsBegin(); elem_i != mr_model_part.ElementsEnd(); ++elem_i)
		{
			Vector u;
			Vector lambda;
			Vector RHS;

			// Get state solution
			elem_i->GetValuesVector(u,0);

			// Get adjoint variables (Corresponds to 1/2*u)
			lambda = 0.5*u;

			// Semi-analytic computation of partial derivative of state equation w.r.t. node coordinates
			elem_i->CalculateRightHandSide(RHS, CurrentProcessInfo);
			for (ModelPart::NodeIterator node_i = elem_i->GetGeometry().begin(); node_i != elem_i->GetGeometry().end(); ++node_i)
			{
				array_3d gradient_contribution(3, 0.0);
				Vector perturbed_RHS = Vector(0);

				// Pertubation, gradient analysis and recovery of x
				node_i->X0() += mDelta;
				elem_i->CalculateRightHandSide(perturbed_RHS, CurrentProcessInfo);
				gradient_contribution[0] = inner_prod(lambda, (perturbed_RHS - RHS) / mDelta);
				node_i->X0() -= mDelta;

				// Reset pertubed vector
				perturbed_RHS = Vector(0);

				// Pertubation, gradient analysis and recovery of y
				node_i->Y0() += mDelta;
				elem_i->CalculateRightHandSide(perturbed_RHS, CurrentProcessInfo);
				gradient_contribution[1] = inner_prod(lambda, (perturbed_RHS - RHS) / mDelta);
				node_i->Y0() -= mDelta;

				// Reset pertubed vector
				perturbed_RHS = Vector(0);

				// Pertubation, gradient analysis and recovery of z
				node_i->Z0() += mDelta;
				elem_i->CalculateRightHandSide(perturbed_RHS, CurrentProcessInfo);
				gradient_contribution[2] = inner_prod(lambda, (perturbed_RHS - RHS) / mDelta);
				node_i->Z0() -= mDelta;

				// Assemble sensitivity to node
				noalias(node_i->FastGetSolutionStepValue(STRAIN_ENERGY_SHAPE_GRADIENT)) += gradient_contribution;
			}
		}

		// Computation of \frac{1}{2} u^T \cdot ( \frac{\partial f_{ext}}{\partial x} )
		const int nconditions = static_cast<int>(mr_model_part.Conditions().size());
		ModelPart::ConditionsContainerType::iterator cond_begin = mr_model_part.ConditionsBegin();
		for (int k = 0; k < nconditions; k++)
		{
			ModelPart::ConditionsContainerType::iterator cond_i = cond_begin + k;

			//detect if the condition is active or not. If the user did not make any choice the element
			//is active by default
			bool condition_is_active = true;
			if ((cond_i)->IsDefined(ACTIVE))
				condition_is_active = (cond_i)->Is(ACTIVE);

			if (condition_is_active)
			{
				Vector u;
				Vector lambda;
				Vector RHS;

				// Get adjoint variables (Corresponds to 1/2*u)
				cond_i->GetValuesVector(u,0);
				lambda = 0.5*u;

				// Semi-analytic computation of partial derivative of force vector w.r.t. node coordinates
				cond_i->CalculateRightHandSide(RHS, CurrentProcessInfo);
				for (ModelPart::NodeIterator node_i = cond_i->GetGeometry().begin(); node_i != cond_i->GetGeometry().end(); ++node_i)
				{
					array_3d gradient_contribution(3, 0.0);
					Vector perturbed_RHS = Vector(0);

					// Pertubation, gradient analysis and recovery of x
					node_i->X0() += mDelta;
					cond_i->CalculateRightHandSide(perturbed_RHS, CurrentProcessInfo);
					gradient_contribution[0] = inner_prod(lambda, (perturbed_RHS - RHS) / mDelta);
					node_i->X0() -= mDelta;

					// Reset pertubed vector
					perturbed_RHS = Vector(0);

					// Pertubation, gradient analysis and recovery of y
					node_i->Y0() += mDelta;
					cond_i->CalculateRightHandSide(perturbed_RHS, CurrentProcessInfo);
					gradient_contribution[1] = inner_prod(lambda, (perturbed_RHS - RHS) / mDelta);
					node_i->Y0() -= mDelta;

					// Reset pertubed vector
					perturbed_RHS = Vector(0);

					// Pertubation, gradient analysis and recovery of z
					node_i->Z0() += mDelta;
					cond_i->CalculateRightHandSide(perturbed_RHS, CurrentProcessInfo);
					gradient_contribution[2] = inner_prod(lambda, (perturbed_RHS - RHS) / mDelta);
					node_i->Z0() -= mDelta;

					// Assemble shape gradient to node
					noalias(node_i->FastGetSolutionStepValue(STRAIN_ENERGY_SHAPE_GRADIENT)) += gradient_contribution;
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
		ProcessInfo &CurrentProcessInfo = mr_model_part.GetProcessInfo();

		// Computation of \frac{1}{2} u^T \cdot ( \frac{\partial f_{ext}}{\partial x} )
		const int nconditions = static_cast<int>(mr_model_part.Conditions().size());
		ModelPart::ConditionsContainerType::iterator cond_begin = mr_model_part.ConditionsBegin();
		for (int k = 0; k < nconditions; k++)
		{
			ModelPart::ConditionsContainerType::iterator cond_i = cond_begin + k;

			//detect if the condition is active or not. If the user did not make any choice the element
			//is active by default
			bool condition_is_active = true;
			if ((cond_i)->IsDefined(ACTIVE))
				condition_is_active = (cond_i)->Is(ACTIVE);

			if (condition_is_active)
			{
				Vector u;
				Vector RHS;

				// Get state solution
				cond_i->GetValuesVector(u,0);

				// Perform finite differencing of RHS vector
				cond_i->CalculateRightHandSide(RHS, CurrentProcessInfo);
				for (ModelPart::NodeIterator node_i = cond_i->GetGeometry().begin(); node_i != cond_i->GetGeometry().end(); ++node_i)
				{
					array_3d gradient_contribution(3, 0.0);
					Vector perturbed_RHS = Vector(0);

					// Pertubation, gradient analysis and recovery of x
					node_i->X0() += mDelta;
					cond_i->CalculateRightHandSide(perturbed_RHS, CurrentProcessInfo);
					gradient_contribution[0] = inner_prod(0.5*u, (perturbed_RHS - RHS) / mDelta);
					node_i->X0() -= mDelta;

					// Reset pertubed vector
					perturbed_RHS = Vector(0);

					// Pertubation, gradient analysis and recovery of y
					node_i->Y0() += mDelta;
					cond_i->CalculateRightHandSide(perturbed_RHS, CurrentProcessInfo);
					gradient_contribution[1] = inner_prod(0.5*u, (perturbed_RHS - RHS) / mDelta);
					node_i->Y0() -= mDelta;

					// Reset pertubed vector
					perturbed_RHS = Vector(0);

					// Pertubation, gradient analysis and recovery of z
					node_i->Z0() += mDelta;
					cond_i->CalculateRightHandSide(perturbed_RHS, CurrentProcessInfo);
					gradient_contribution[2] = inner_prod(0.5*u, (perturbed_RHS - RHS) / mDelta);
					node_i->Z0() -= mDelta;

					// Assemble shape gradient to node
					noalias(node_i->FastGetSolutionStepValue(STRAIN_ENERGY_SHAPE_GRADIENT)) += gradient_contribution;
				}
			}
		}
		KRATOS_CATCH("");
	}

	// --------------------------------------------------------------------------
  	virtual void ConsiderDiscretization(){


		// Start process to identify element neighbors for every node
		FindNodalNeighboursProcess neigbhorFinder = FindNodalNeighboursProcess(mr_model_part, 10, 10);
		neigbhorFinder.Execute();

		std::cout<< "> Considering discretization size!" << std::endl;
		for(ModelPart::NodeIterator node_i=mr_model_part.NodesBegin(); node_i!=mr_model_part.NodesEnd(); node_i++)
		{
			WeakPointerVector<Element >& ng_elem = node_i->GetValue(NEIGHBOUR_ELEMENTS);

			double scaling_factor = 0.0;
			for(unsigned int i = 0; i < ng_elem.size(); i++)
			{
				Kratos::Element& ng_elem_i = ng_elem[i];
				Element::GeometryType& element_geometry = ng_elem_i.GetGeometry();

				if( isElementOfTypeShell(element_geometry) )
					scaling_factor += element_geometry.Area();
				else
					scaling_factor += element_geometry.Volume();
			}

			node_i->FastGetSolutionStepValue(STRAIN_ENERGY_SHAPE_GRADIENT) /= scaling_factor;
		}
	}

	// --------------------------------------------------------------------------
	bool isElementOfTypeShell( Element::GeometryType& given_element_geometry )
	{
		if(given_element_geometry.WorkingSpaceDimension() != given_element_geometry.LocalSpaceDimension())
			return true;
		else
		    return false;
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

	ModelPart &mr_model_part;
	unsigned int mGradientMode;
	double m_strain_energy;
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
	//      StrainEnergyResponseFunction& operator=(StrainEnergyResponseFunction const& rOther);

	/// Copy constructor.
	//      StrainEnergyResponseFunction(StrainEnergyResponseFunction const& rOther);

	///@}

}; // Class StrainEnergyResponseFunction

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

} // namespace Kratos.

#endif // STRAIN_ENERGY_RESPONSE_FUNCTION_H

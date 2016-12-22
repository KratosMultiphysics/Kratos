// ==============================================================================
/*
 KratosShapeOptimizationApplication
 A library based on:
 Kratos
 A General Purpose Software for Multi-Physics Finite Element Analysis
 (Released on march 05, 2007).

 Copyright (c) 2016: Daniel Baumgaertner
                     daniel.baumgaertner@tum.de
                     Chair of Structural Analysis
                     Technische Universitaet Muenchen
                     Arcisstrasse 21 80333 Munich, Germany

 Permission is hereby granted, free  of charge, to any person obtaining
 a  copy  of this  software  and  associated  documentation files  (the
 "Software"), to  deal in  the Software without  restriction, including
 without limitation  the rights to  use, copy, modify,  merge, publish,
 distribute,  sublicense and/or  sell copies  of the  Software,  and to
 permit persons to whom the Software  is furnished to do so, subject to
 the following condition:

 Distribution of this code for  any  commercial purpose  is permissible
 ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

 The  above  copyright  notice  and  this permission  notice  shall  be
 included in all copies or substantial portions of the Software.

 THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
 EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
 MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
 CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
 TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
 SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */
//==============================================================================
//
//   Project Name:        KratosShape                            $
//   Created by:          $Author:    daniel.baumgaertner@tum.de $
//                        $Author:           armin.geiser@tum.de $
//   Date:                $Date:                   December 2016 $
//   Revision:            $Revision:                         0.0 $
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
// External includes
// ------------------------------------------------------------------------------
#include <boost/python.hpp>
#include <boost/numeric/ublas/io.hpp>

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "../../kratos/includes/define.h"
#include "../../kratos/processes/process.h"
#include "../../kratos/includes/node.h"
#include "../../kratos/includes/element.h"
#include "../../kratos/includes/model_part.h"
#include "../../kratos/includes/kratos_flags.h"
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
	StrainEnergyResponseFunction(ModelPart &model_part, boost::python::dict response_settings)
	: mr_model_part(model_part)
	{
		// Set gradient mode
		boost::python::extract<const char *> grad_mode(response_settings["gradient_mode"]);

		// Mode 1: analytic sensitivities
		if (std::strcmp(grad_mode, "analytic") == 0)
			m_gradient_mode = 1;

		// Mode 2: semi-analytic sensitivities
		else if (std::strcmp(grad_mode, "semi_analytic") == 0)
		{
			m_gradient_mode = 2;
			boost::python::extract<double> delta(response_settings["step_size"]);
			m_delta = delta;
		}

		// Throw error message in case of wrong specification
		else
			KRATOS_THROW_ERROR(std::invalid_argument, "Specified gradient_mode not recognized. Options are: analytic , semi_analytic. Specified gradient_mode: ", grad_mode);

		// Initialize member variables to NULL
		m_initial_value = 0.0;
		m_initial_value_defined = false;
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
	void initialize()
	{
		// In case of analytic sensitivity analysis, check if specified elements in model_part provide necessary sensitivity information.
		// To check, we compare if the class type of all the given elements in the model_part is among the elements
		// that provide the required sensitivity information (reference elements)
		// The reference class type is: "SmallDisplacementAnalyticSensitivityElement"

		if(m_gradient_mode==1)
		{
			const char element_name[] = "SmallDisplacementAnalyticSensitivityElement3D4N";
			Element const &reference_element = KratosComponents<Element>::Get(element_name);

			bool sensitivity_analysis_implemented = true;
			for (ModelPart::ElementsContainerType::iterator elem_i = mr_model_part.ElementsBegin(); elem_i != mr_model_part.ElementsEnd(); ++elem_i)
			{
				Element const &given_element = mr_model_part.Elements()[elem_i->Id()];

				if (typeid(given_element) != typeid(reference_element))
				{
					sensitivity_analysis_implemented = false;
					break;
				}
			}

			if (!sensitivity_analysis_implemented)
				KRATOS_THROW_ERROR(std::logic_error, "Analytic sensitivity analysis for given element type not implemented. Please choose for complete model part elements that support an analytic sensitivity analysis", "");
		}
	}

	// --------------------------------------------------------------------------
	void calculate_value()
	{
		KRATOS_TRY;

		// Working variables / vectors
		ProcessInfo &CurrentProcessInfo = mr_model_part.GetProcessInfo();
		m_strain_energy = 0.0;
		Vector u;
		Vector RHS;

		// Computation of strain_energy = 1/2*u*f
		const int nconditions = static_cast<int>(mr_model_part.Conditions().size());
		ModelPart::ConditionsContainerType::iterator cond_begin = mr_model_part.ConditionsBegin();
		for (int k = 0; k < nconditions; k++)
		{
			ModelPart::ConditionsContainerType::iterator cond_i = cond_begin + k;

			// Detect if the condition is active or not. If the user did not make any choice the element
			// Is active by default
			bool condition_is_active = true;
			if ((cond_i)->IsDefined(ACTIVE))
				condition_is_active = (cond_i)->Is(ACTIVE);

			if (condition_is_active)
			{
				// Get state solution relevant for energy calculation
				cond_i->GetValuesVector(u,0);

				// Calculate RHS
				cond_i->CalculateRightHandSide(RHS, CurrentProcessInfo);

				// Compute strain energy
				m_strain_energy += 0.5 * inner_prod(RHS, u);
			}
		}

		// Set initial value if not done yet
		if(!m_initial_value_defined)
		{
			m_initial_value = m_strain_energy;
			m_initial_value_defined = true;
		}


		KRATOS_CATCH("");
	}

	// --------------------------------------------------------------------------
	void calculate_gradient()
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

		switch (m_gradient_mode)
		{
		// analytic sensitivities
		case 1:
		{
			std::cout << "WARNING: in StrainEnergyResponseFunction::calculate_gradient()!!!! No variation of external force considerd yet in analytical sensitivity analysis" << std::endl;

			// calculate_response_derivative_part_analytically();
			calculate_adjoint_field();
			calculate_state_derivative_part_analytically();
			break;
		}
		// Semi analytic sensitivities
		case 2:
		{
			calculate_response_derivative_part_by_finite_differencing();
			calculate_adjoint_field();
			calculate_state_derivative_part_by_finite_differencing();
			break;
		}
		}

		KRATOS_CATCH("");
	}
	// --------------------------------------------------------------------------
	double get_initial_value()
	{
		KRATOS_TRY;

		if(!m_initial_value_defined)
			KRATOS_THROW_ERROR(std::logi:error, "Initial value not yet defined! First compute it by calling \"calculate_value()\"", m_initial_value_defined);

		return m_initial_value;

		KRATOS_CATCH("");
	}

	// --------------------------------------------------------------------------
	double get_value()
	{
		KRATOS_TRY;

		return m_strain_energy;

		KRATOS_CATCH("");
	}

	// --------------------------------------------------------------------------
	boost::python::dict get_gradient()
	{
		KRATOS_TRY;

		// Dictionary to store all sensitivities along with Ids of corresponding nodes
		boost::python::dict dFdX;

		// Fill dictionary with gradient information
		for (ModelPart::NodeIterator node_i = mr_model_part.NodesBegin(); node_i != mr_model_part.NodesEnd(); ++node_i)
			dFdX[node_i->Id()] = node_i->FastGetSolutionStepValue(STRAIN_ENERGY_SHAPE_GRADIENT);

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
	void calculate_adjoint_field()
	{
		KRATOS_TRY;

		// Adjoint field may be directly obtained from state solution

		KRATOS_CATCH("");
	}

	// --------------------------------------------------------------------------
	void calculate_state_derivative_part_analytically()
	{
		KRATOS_TRY;

		// Working variables
		ProcessInfo &CurrentProcessInfo = mr_model_part.GetProcessInfo();

		// Computation of: \frac{1}{2} u^T \cdot ( - \frac{\partial K}{\partial x} )
		for (ModelPart::ElementIterator elem_i = mr_model_part.ElementsBegin(); elem_i != mr_model_part.ElementsEnd(); ++elem_i)
		{
			Vector u;
			Vector lambda;

			// Get state solution
			elem_i->GetValuesVector(u,0);

			// Get adjoint variables (Corresponds to 1/2*u)
			lambda = 0.5*u;

			// Analytic computation of partial derivative of state equation w.r.t. node coordinates
			for (ModelPart::NodeIterator node_i = elem_i->GetGeometry().begin(); node_i != elem_i->GetGeometry().end(); ++node_i)
			{
				array_3d gradient_contribution(3, 0.0);
				int node_index = node_i - elem_i->GetGeometry().begin();

				// Specify node for which DKDXU (including DKDXU_X,DKDXU_Y,DKDXU_Z) shall be computed
				elem_i->SetValue(ACTIVE_NODE_INDEX, node_index);
				elem_i->Calculate(DKDXU, u, CurrentProcessInfo);

				gradient_contribution[0] = -inner_prod(lambda, elem_i->GetValue(DKDXU_X));
				gradient_contribution[1] = -inner_prod(lambda, elem_i->GetValue(DKDXU_Y));
				gradient_contribution[2] = -inner_prod(lambda, elem_i->GetValue(DKDXU_Z));

				// Assemble gradient to node
				noalias(node_i->FastGetSolutionStepValue(STRAIN_ENERGY_SHAPE_GRADIENT)) += gradient_contribution;
			}
		}

		KRATOS_CATCH("");
	}

	// --------------------------------------------------------------------------
	void calculate_state_derivative_part_by_finite_differencing()
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
				node_i->X0() += m_delta;
				elem_i->CalculateRightHandSide(perturbed_RHS, CurrentProcessInfo);
				gradient_contribution[0] = inner_prod(lambda, (perturbed_RHS - RHS) / m_delta);
				node_i->X0() -= m_delta;

				// Reset pertubed vector
				perturbed_RHS = Vector(0);

				// Pertubation, gradient analysis and recovery of y
				node_i->Y0() += m_delta;
				elem_i->CalculateRightHandSide(perturbed_RHS, CurrentProcessInfo);
				gradient_contribution[1] = inner_prod(lambda, (perturbed_RHS - RHS) / m_delta);
				node_i->Y0() -= m_delta;

				// Reset pertubed vector
				perturbed_RHS = Vector(0);

				// Pertubation, gradient analysis and recovery of z
				node_i->Z0() += m_delta;
				elem_i->CalculateRightHandSide(perturbed_RHS, CurrentProcessInfo);
				gradient_contribution[2] = inner_prod(lambda, (perturbed_RHS - RHS) / m_delta);
				node_i->Z0() -= m_delta;

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
					node_i->X0() += m_delta;
					cond_i->CalculateRightHandSide(perturbed_RHS, CurrentProcessInfo);
					gradient_contribution[0] = inner_prod(lambda, (perturbed_RHS - RHS) / m_delta);
					node_i->X0() -= m_delta;

					// Reset pertubed vector
					perturbed_RHS = Vector(0);

					// Pertubation, gradient analysis and recovery of y
					node_i->Y0() += m_delta;
					cond_i->CalculateRightHandSide(perturbed_RHS, CurrentProcessInfo);
					gradient_contribution[1] = inner_prod(lambda, (perturbed_RHS - RHS) / m_delta);
					node_i->Y0() -= m_delta;

					// Reset pertubed vector
					perturbed_RHS = Vector(0);

					// Pertubation, gradient analysis and recovery of z
					node_i->Z0() += m_delta;
					cond_i->CalculateRightHandSide(perturbed_RHS, CurrentProcessInfo);
					gradient_contribution[2] = inner_prod(lambda, (perturbed_RHS - RHS) / m_delta);
					node_i->Z0() -= m_delta;

					// Assemble shape gradient to node
					noalias(node_i->FastGetSolutionStepValue(STRAIN_ENERGY_SHAPE_GRADIENT)) += gradient_contribution;
				}
			}
		}

		KRATOS_CATCH("");
	}

	// --------------------------------------------------------------------------
	void calculate_response_derivative_part_by_finite_differencing()
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
					node_i->X0() += m_delta;
					cond_i->CalculateRightHandSide(perturbed_RHS, CurrentProcessInfo);
					gradient_contribution[0] = inner_prod(0.5*u, (perturbed_RHS - RHS) / m_delta);
					node_i->X0() -= m_delta;

					// Reset pertubed vector
					perturbed_RHS = Vector(0);

					// Pertubation, gradient analysis and recovery of y
					node_i->Y0() += m_delta;
					cond_i->CalculateRightHandSide(perturbed_RHS, CurrentProcessInfo);
					gradient_contribution[1] = inner_prod(0.5*u, (perturbed_RHS - RHS) / m_delta);
					node_i->Y0() -= m_delta;

					// Reset pertubed vector
					perturbed_RHS = Vector(0);

					// Pertubation, gradient analysis and recovery of z
					node_i->Z0() += m_delta;
					cond_i->CalculateRightHandSide(perturbed_RHS, CurrentProcessInfo);
					gradient_contribution[2] = inner_prod(0.5*u, (perturbed_RHS - RHS) / m_delta);
					node_i->Z0() -= m_delta;

					// Assemble shape gradient to node
					noalias(node_i->FastGetSolutionStepValue(STRAIN_ENERGY_SHAPE_GRADIENT)) += gradient_contribution;
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

	ModelPart &mr_model_part;
	unsigned int m_gradient_mode;
	double m_strain_energy;
	double m_delta;
	double m_initial_value;
	bool m_initial_value_defined;

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

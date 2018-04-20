// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Fusseder Martin
//

#ifndef EIGENFREQUENCY_RESPONSE_FUNCTION_LIN_SCAL_UTILITY_H
#define EIGENFREQUENCY_RESPONSE_FUNCTION_LIN_SCAL_UTILITY_H

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

//template<class TDenseSpace>

class EigenfrequencyResponseFunctionLinScalUtility
{
public:
	///@name Type Definitions
	///@{

	/// Pointer definition of EigenfrequencyResponseFunctionLinScalUtility
	KRATOS_CLASS_POINTER_DEFINITION(EigenfrequencyResponseFunctionLinScalUtility);

	///@}
	///@name Life Cycle
	///@{

	/// Default constructor.
	EigenfrequencyResponseFunctionLinScalUtility(ModelPart& model_part, Parameters responseSettings)
	: mrModelPart(model_part)
	{
		// Set gradient mode
		const std::string gradient_mode = responseSettings["gradient_mode"].GetString();

		// Mode 1: semi-analytic sensitivities
		if (gradient_mode.compare("semi_analytic") == 0)
		{
			mGradientMode = 1;
			mDelta = responseSettings["step_size"].GetDouble();
		}
		else
			KRATOS_ERROR << "Specified gradient_mode '" << gradient_mode << "' not recognized. The only option is: semi_analytic" << std::endl;


        // Get array of numbers of the eigenfrequencies which have to be traced by this response function
		mTracedEigenfrequencies.resize(responseSettings["traced_eigenfrequency"].size(),false);

		for(std::size_t i = 0; i < mTracedEigenfrequencies.size(); i++)
			mTracedEigenfrequencies[i] = responseSettings["traced_eigenfrequency"][i].GetInt();

		// Ask for weighting factors and check their number
		KRATOS_ERROR_IF(responseSettings["weighting_factors"].size() != mTracedEigenfrequencies.size()) << "The number of chosen eigenvalues does not fit to the number of weighting factors!" << std::endl;
		mWeightingFactors.resize(responseSettings["weighting_factors"].size(),false);

		// Read weighting factors and sum them up
		double test_sum_weight_facs = 0.0;
		for(std::size_t i = 0; i < mWeightingFactors.size(); i++)
		{
			mWeightingFactors[i] = responseSettings["weighting_factors"][i].GetDouble();
			test_sum_weight_facs += mWeightingFactors[i];
		}

		// Check the weighting factors and modify them for the case that their sum is unequal to one
		if(std::abs(test_sum_weight_facs - 1.0) > 1e-12)
		{
			for(std::size_t i = 0; i < mTracedEigenfrequencies.size(); i++)
				mWeightingFactors[i] /= test_sum_weight_facs;

			std::cout << "> The sum of the chosen weighting factors is unequal to one. A scaling process was executed for them!" << std::endl;
		}
	}

	/// Destructor.
	virtual ~EigenfrequencyResponseFunctionLinScalUtility()
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
	double CalculateValue()
	{
		KRATOS_TRY;

		double m_resp_function_value = 0.0;

		// Compute response function by weighting with linear scalarization
		for(std::size_t i = 0; i < mTracedEigenfrequencies.size(); i++)
			m_resp_function_value += mWeightingFactors[i] * get_single_eigenvalue(mTracedEigenfrequencies[i]);

		return m_resp_function_value;

		KRATOS_CATCH("");
	}

	// --------------------------------------------------------------------------
	double get_single_eigenvalue(int id_eigenvalue)
	{
		KRATOS_TRY;

		const int num_of_computed_eigenvalues = (mrModelPart.GetProcessInfo()[EIGENVALUE_VECTOR]).size();

		KRATOS_ERROR_IF(num_of_computed_eigenvalues < id_eigenvalue) << "The chosen eigenvalue was not solved by the eigenvalue analysis!" << std::endl;

		return (mrModelPart.GetProcessInfo()[EIGENVALUE_VECTOR])[id_eigenvalue-1];

		KRATOS_CATCH("");
	}

	// --------------------------------------------------------------------------
	Vector get_eigenvector_of_element(ModelPart::ElementType& traced_element, int id_eigenvalue, int size_of_eigenvector)
	{
		KRATOS_TRY;

		Vector eigenvector_of_element;
		eigenvector_of_element.resize(size_of_eigenvector,false);

		// Get eigenvector of element
		int k = 0;
		const int NumNodeDofs = size_of_eigenvector/traced_element.GetGeometry().size();
		for (auto& node_i : traced_element.GetGeometry())
		{
			Matrix& rNodeEigenvectors = node_i.GetValue(EIGENVECTOR_MATRIX);

			for (int i = 0; i < NumNodeDofs; i++)
				eigenvector_of_element(i+NumNodeDofs*k) = rNodeEigenvectors((id_eigenvalue-1),i);

			k++;
		}

		return eigenvector_of_element;

		KRATOS_CATCH("");
	}

	// --------------------------------------------------------------------------
	void CalculateGradient()
	{
		KRATOS_TRY;

		// First gradients are initialized
		VariableUtils().SetToZero_VectorVar(SHAPE_SENSITIVITY, mrModelPart.Nodes());

		// Gradient calculation is done by a semi-analytic approaches
		// The gradient is computed in one step

		switch (mGradientMode)
		{
		// Semi-analytic sensitivities
		case 1:
		{
			calculate_complete_response_derivative_by_finite_differencing();
			break;
		}

		}// End switch mGradientMode

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
		return "EigenfrequencyResponseFunctionLinScalUtility";
	}

	/// Print information about this object.
	virtual void PrintInfo(std::ostream &rOStream) const
	{
		rOStream << "EigenfrequencyResponseFunctionLinScalUtility";
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
	void calculate_complete_response_derivative_by_finite_differencing()
	{
		KRATOS_TRY;

		// Working variables
		ProcessInfo &CurrentProcessInfo = mrModelPart.GetProcessInfo();

		for (auto& elem_i : mrModelPart.Elements())
		{
			Matrix mass_matrix_org;
			Matrix LHS_org;
			Vector eigenvector_of_element = Vector(0);
			Vector dummy;
			Matrix aux_matrix = Matrix(0,0);
			Vector aux_vector = Vector(0);
			elem_i.CalculateMassMatrix(mass_matrix_org, CurrentProcessInfo);
			elem_i.CalculateLocalSystem(LHS_org, dummy ,CurrentProcessInfo);
			double traced_eigenvalue = 0.0;

			// Get size of element eigenvector and initialize eigenvector
			int num_dofs_element = mass_matrix_org.size1();
			eigenvector_of_element.resize(num_dofs_element,false);

			// Semi-analytic computation of partial derivative of state equation w.r.t. node coordinates
			for (auto& node_i : elem_i.GetGeometry())
			{
				Vector gradient_contribution(3, 0.0);
				Matrix perturbed_LHS = Matrix(0,0);
				Matrix perturbed_mass_matrix = Matrix(0,0);

				// Derivative of response w.r.t. x-coord ------------------------
				node_i.X0() += mDelta;

				elem_i.CalculateMassMatrix(perturbed_mass_matrix, CurrentProcessInfo);
				perturbed_mass_matrix = (perturbed_mass_matrix - mass_matrix_org) / mDelta;

				elem_i.CalculateLocalSystem(perturbed_LHS, dummy ,CurrentProcessInfo);
				perturbed_LHS = (perturbed_LHS - LHS_org) / mDelta;

				// Loop over eigenvalues
				for(std::size_t i = 0; i < mTracedEigenfrequencies.size(); i++)
				{
					traced_eigenvalue = get_single_eigenvalue(mTracedEigenfrequencies[i]);
					eigenvector_of_element = get_eigenvector_of_element(elem_i, mTracedEigenfrequencies[i], num_dofs_element);

					aux_matrix = perturbed_LHS;
					aux_matrix -= (perturbed_mass_matrix * traced_eigenvalue);
					aux_vector = prod(aux_matrix , eigenvector_of_element);
					gradient_contribution[0] += inner_prod(eigenvector_of_element , aux_vector) * mWeightingFactors[i];

					eigenvector_of_element = Vector(0);
					aux_matrix = Matrix(0,0);
					aux_vector = Vector(0);
				}

				node_i.X0() -= mDelta;

				// End derivative of response w.r.t. x-coord --------------------

				// Reset pertubed RHS and mass matrix
				perturbed_LHS = Matrix(0,0);
				perturbed_mass_matrix = Matrix(0,0);


				// Derivative of response w.r.t. y-coord ------------------------
				node_i.Y0() += mDelta;

				elem_i.CalculateMassMatrix(perturbed_mass_matrix, CurrentProcessInfo);
				perturbed_mass_matrix = (perturbed_mass_matrix - mass_matrix_org) / mDelta;

				elem_i.CalculateLocalSystem(perturbed_LHS, dummy ,CurrentProcessInfo);
				perturbed_LHS = (perturbed_LHS - LHS_org) / mDelta;

				// Loop over eigenvalues
				for(std::size_t i = 0; i < mTracedEigenfrequencies.size(); i++)
				{
					traced_eigenvalue = get_single_eigenvalue(mTracedEigenfrequencies[i]);
					eigenvector_of_element = get_eigenvector_of_element(elem_i, mTracedEigenfrequencies[i], num_dofs_element);

					aux_matrix = perturbed_LHS;
					aux_matrix -= (perturbed_mass_matrix * traced_eigenvalue);
					aux_vector = prod(aux_matrix , eigenvector_of_element);
					gradient_contribution[1] += inner_prod(eigenvector_of_element , aux_vector) * mWeightingFactors[i];

					eigenvector_of_element = Vector(0);
				    aux_matrix = Matrix(0,0);
					aux_vector = Vector(0);;
				}

				node_i.Y0() -= mDelta;
				// End derivative of response w.r.t. y-coord --------------------

				// Reset pertubed RHS and mass matrix
				perturbed_LHS = Matrix(0,0);
				perturbed_mass_matrix= Matrix(0,0);

				// Derivative of response w.r.t. z-coord ------------------------
				node_i.Z0() += mDelta;

				elem_i.CalculateMassMatrix(perturbed_mass_matrix, CurrentProcessInfo);
				perturbed_mass_matrix = (perturbed_mass_matrix - mass_matrix_org) / mDelta;

				elem_i.CalculateLocalSystem(perturbed_LHS, dummy ,CurrentProcessInfo);
				perturbed_LHS = (perturbed_LHS - LHS_org) / mDelta;

				// Loop over eigenvalues
				for(std::size_t i = 0; i < mTracedEigenfrequencies.size(); i++)
				{
					traced_eigenvalue = get_single_eigenvalue(mTracedEigenfrequencies[i]);
					eigenvector_of_element = get_eigenvector_of_element(elem_i, mTracedEigenfrequencies[i], num_dofs_element);

					aux_matrix = perturbed_LHS;
					aux_matrix -= (perturbed_mass_matrix * traced_eigenvalue);
					aux_vector = prod(aux_matrix , eigenvector_of_element);
					gradient_contribution[2] += inner_prod(eigenvector_of_element , aux_vector) * mWeightingFactors[i];

					eigenvector_of_element = Vector(0);
					aux_matrix = Matrix(0,0);
					aux_vector = Vector(0);
				}

				node_i.Z0() -= mDelta;
				// End derivative of response w.r.t. z-coord --------------------

				// Assemble sensitivity to node
				noalias(node_i.FastGetSolutionStepValue(SHAPE_SENSITIVITY)) += gradient_contribution;

			}// End loop over nodes of element

		}// End loop over elements
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
	unsigned int mGradientMode;
	double mDelta;
	std::vector<int> mTracedEigenfrequencies;
	std::vector<double> mWeightingFactors;

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
	//      EigenfrequencyResponseFunctionLinScalUtility& operator=(EigenfrequencyResponseFunctionLinScalUtility const& rOther);

	/// Copy constructor.
	//      EigenfrequencyResponseFunctionLinScalUtility(EigenfrequencyResponseFunctionLinScalUtility const& rOther);

	///@}

}; // Class EigenfrequencyResponseFunctionLinScalUtility

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

} // namespace Kratos.

#endif // EIGENFREQUENCY_RESPONSE_FUNCTION_LIN_SCAL_UTILITY_H

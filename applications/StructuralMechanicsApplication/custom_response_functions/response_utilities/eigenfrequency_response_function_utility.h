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

#ifndef EIGENFREQUENCY_RESPONSE_FUNCTION_UTILITY_H
#define EIGENFREQUENCY_RESPONSE_FUNCTION_UTILITY_H

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

class EigenfrequencyResponseFunctionUtility
{
public:
	///@name Type Definitions
	///@{

	/// Pointer definition of EigenfrequencyResponseFunctionUtility
	KRATOS_CLASS_POINTER_DEFINITION(EigenfrequencyResponseFunctionUtility);

	///@}
	///@name Life Cycle
	///@{

	/// Default constructor.
	EigenfrequencyResponseFunctionUtility(ModelPart& model_part, Parameters responseSettings)
	: mrModelPart(model_part)
	{
		std::string gradient_mode = responseSettings["gradient_mode"].GetString();
		if (gradient_mode.compare("semi_analytic") == 0)
		{
			mDelta = responseSettings["step_size"].GetDouble();
		}
		else
			KRATOS_ERROR << "Specified gradient_mode '" << gradient_mode << "' not recognized. The only option is: semi_analytic" << std::endl;

		// Get number of eigenfrequency for which the structure has to be optimized
		mTracedEigenValue = responseSettings["traced_eigenfrequency"].GetInt();
	}

	/// Destructor.
	virtual ~EigenfrequencyResponseFunctionUtility()
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

		double eigenvalue = 0.0;

		int num_of_computed_eigenvalues = (mrModelPart.GetProcessInfo()[EIGENVALUE_VECTOR]).size();

		KRATOS_ERROR_IF(num_of_computed_eigenvalues < mTracedEigenValue) << "The chosen eigenvalue was not solved by the eigenvalue analysis!" << std::endl;

		eigenvalue = (mrModelPart.GetProcessInfo()[EIGENVALUE_VECTOR])[mTracedEigenValue - 1];

		return eigenvalue;

		KRATOS_CATCH("");
	}

	// --------------------------------------------------------------------------
	void CalculateGradient()
	{
		KRATOS_TRY;

		// Formula computed in general notation:
		// \frac{dF}{dx} = eigenvector^T\cdot (frac{\partial RHS}{\partial x} -
		//				   eigenvalue frac{\partial mass_matrix}{\partial x})\cdot eigenvector

		// First gradients are initialized
		VariableUtils().SetToZero_VectorVar(SHAPE_SENSITIVITY, mrModelPart.Nodes());

		// Gradient calculation is done by a semi-analytic approach
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
		return "EigenfrequencyResponseFunctionUtility";
	}

	/// Print information about this object.
	virtual void PrintInfo(std::ostream &rOStream) const
	{
		rOStream << "EigenfrequencyResponseFunctionUtility";
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

		const double eigenvalue = (mrModelPart.GetProcessInfo()[EIGENVALUE_VECTOR])[mTracedEigenValue - 1];

		// Computation of: \frac{dF}{dx} = eigenvector^T\cdot (frac{\partial RHS}{\partial x} -
		//				                   eigenvalue frac{\partial mass_matrix}{\partial x})\cdot eigenvector
		for (auto& elem_i : mrModelPart.Elements())
		{
			Matrix mass_matrix_org;
			Matrix LHS_org;
			Vector eigenvector_of_element = Vector(0);
			Vector dummy;
			elem_i.CalculateMassMatrix(mass_matrix_org, CurrentProcessInfo);
			elem_i.CalculateLocalSystem(LHS_org, dummy ,CurrentProcessInfo);

			// Get size of element eigenvector and initialize eigenvector
			int num_dofs_element = mass_matrix_org.size1();
			eigenvector_of_element.resize(num_dofs_element,false);

			// Get eigenvector of element
			int k = 0;
			const int NumNodeDofs = num_dofs_element/elem_i.GetGeometry().size();
			for (auto& node_i : elem_i.GetGeometry())
			{
				Matrix& rNodeEigenvectors = node_i.GetValue(EIGENVECTOR_MATRIX);

				for (int i = 0; i < NumNodeDofs; i++)
                    eigenvector_of_element(i+NumNodeDofs*k) = rNodeEigenvectors((mTracedEigenValue-1),i);

				k++;
			}

			// Semi-analytic computation of partial derivative of state equation w.r.t. node coordinates
			for (auto& node_i : elem_i.GetGeometry())
			{

				Vector gradient_contribution(3, 0.0);
				Matrix perturbed_LHS = Matrix(0,0);
				Matrix perturbed_mass_matrix = Matrix(0,0);
				Vector aux = Vector(0);

				// Derivative of response w.r.t. x-coord ------------------------
				node_i.X0() += mDelta;

				elem_i.CalculateMassMatrix(perturbed_mass_matrix, CurrentProcessInfo);
				perturbed_mass_matrix = (perturbed_mass_matrix - mass_matrix_org) / mDelta;
				perturbed_mass_matrix *= eigenvalue;

				elem_i.CalculateLocalSystem(perturbed_LHS, dummy ,CurrentProcessInfo);
				perturbed_LHS = (perturbed_LHS - LHS_org) / mDelta;

				perturbed_LHS -= perturbed_mass_matrix;

				aux = prod(perturbed_LHS , eigenvector_of_element);
				gradient_contribution[0] = inner_prod(eigenvector_of_element , aux);

				node_i.X0() -= mDelta;
				// End derivative of response w.r.t. x-coord --------------------

				// Reset pertubed RHS and mass matrix
				perturbed_LHS = Matrix(0,0);
				perturbed_mass_matrix = Matrix(0,0);
				aux = Vector(0);

				// Derivative of response w.r.t. y-coord ------------------------
				node_i.Y0() += mDelta;

				elem_i.CalculateMassMatrix(perturbed_mass_matrix, CurrentProcessInfo);

				perturbed_mass_matrix = (perturbed_mass_matrix - mass_matrix_org) / mDelta;
				perturbed_mass_matrix *= eigenvalue;

				elem_i.CalculateLocalSystem(perturbed_LHS, dummy ,CurrentProcessInfo);

				perturbed_LHS = (perturbed_LHS - LHS_org) / mDelta;

				perturbed_LHS -= perturbed_mass_matrix;

				aux = prod(perturbed_LHS , eigenvector_of_element);
				gradient_contribution[1] = inner_prod(eigenvector_of_element , aux);
				node_i.Y0() -= mDelta;
				// End derivative of response w.r.t. y-coord --------------------

				// Reset pertubed RHS and mass matrix
				perturbed_LHS = Matrix(0,0);
				perturbed_mass_matrix= Matrix(0,0);
				aux = Vector(0);

				// Derivative of response w.r.t. z-coord ------------------------
				node_i.Z0() += mDelta;

				elem_i.CalculateMassMatrix(perturbed_mass_matrix, CurrentProcessInfo);
				perturbed_mass_matrix = (perturbed_mass_matrix - mass_matrix_org) / mDelta;
				perturbed_mass_matrix *= eigenvalue;

				elem_i.CalculateLocalSystem(perturbed_LHS, dummy ,CurrentProcessInfo);
				perturbed_LHS = (perturbed_LHS - LHS_org) / mDelta;

				perturbed_LHS -= perturbed_mass_matrix;

				aux = prod(perturbed_LHS , eigenvector_of_element);
				gradient_contribution[2] =  inner_prod(eigenvector_of_element , aux);
				node_i.Z0() -= mDelta;
				// End derivative of response w.r.t. z-coord --------------------

				// Assemble sensitivity to node
				noalias(node_i.FastGetSolutionStepValue(SHAPE_SENSITIVITY)) += gradient_contribution;
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
	unsigned int mGradientMode;
	unsigned int mWeightingMethod;
	double mDelta;
	int mTracedEigenValue;


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
	//      EigenfrequencyResponseFunctionUtility& operator=(SEigenfrequencyResponseFunctionUtility const& rOther);

	/// Copy constructor.
	//      EigenfrequencyResponseFunctionUtility(EigenfrequencyResponseFunctionUtility const& rOther);

	///@}

}; // Class EigenfrequencyResponseFunctionUtility

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

} // namespace Kratos.

#endif // EIGENFRQUENCY_RESPONSE_FUNCTION_UTILITY_H

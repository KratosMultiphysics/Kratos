// ==============================================================================
//  KratosShapeOptimizationApplication
//
//  License:         BSD License
//                   license: ShapeOptimizationApplication/license.txt
//
//  Main authors:    Fusseder Martin
//                   martin.fusseder@tum.de
//
// ==============================================================================

#ifndef EIGENFREQUENCY_RESPONSE_FUNCTION_H
#define EIGENFREQUENCY_RESPONSE_FUNCTION_H

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
#include "includes/define.h"
#include "processes/process.h"
#include "includes/node.h"
#include "includes/element.h"
#include "includes/model_part.h"
#include "includes/kratos_flags.h"
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

//template<class TDenseSpace>

class EigenfrequencyResponseFunction : public ResponseFunction
{
public:
	///@name Type Definitions
	///@{

	// TODO solve this via template or how to get this from Eigensolverstrategy
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;

	typedef array_1d<double, 3> array_3d;
    typedef LocalSpaceType::VectorType DenseVectorType;
	typedef LocalSpaceType::MatrixType DenseMatrixType;
	typedef Variable<DenseVectorType> VariableDenseVectorType;
	typedef Variable<DenseMatrixType> VariableDenseMatrixType;

	/// Pointer definition of EigenfrequencyResponseFunction
	KRATOS_CLASS_POINTER_DEFINITION(EigenfrequencyResponseFunction);

	///@}
	///@name Life Cycle
	///@{

	/// Default constructor.
	EigenfrequencyResponseFunction(ModelPart& model_part, Parameters& responseSettings)
	: mrModelPart(model_part)
	{
		std::string gradientMode = responseSettings["gradient_mode"].GetString();
		if (gradientMode.compare("semi_analytic") == 0)
		{
			mDelta = responseSettings["step_size"].GetDouble();
		}
		else
			KRATOS_ERROR << "Specified gradient_mode not recognized. The only option is: semi_analytic. Specified gradient_mode: " << gradientMode << std::endl;

		// Get number of eigenfrequency for which the structure has to be optimized
		mTracedEigenValue = responseSettings["traced_eigenfrequency"].GetInt();
	}

	/// Destructor.
	virtual ~EigenfrequencyResponseFunction()
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

		const VariableDenseVectorType& rEIGENVALUE_VECTOR =
            KratosComponents<VariableDenseVectorType>::Get("EIGENVALUE_VECTOR");

		int num_of_computed_eigenvalues = (mrModelPart.GetProcessInfo()[rEIGENVALUE_VECTOR]).size();

		KRATOS_ERROR_IF(num_of_computed_eigenvalues < mTracedEigenValue) << "The chosen eigenvalue was not solved by the eigenvalue analysis!" << std::endl;

		eigenvalue = (mrModelPart.GetProcessInfo()[rEIGENVALUE_VECTOR])[mTracedEigenValue - 1];

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
		array_3d zeros_array(3, 0.0);
		for (auto& node_i : mrModelPart.Nodes())
			noalias(node_i.FastGetSolutionStepValue(SHAPE_SENSITIVITY) )= zeros_array;

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
	virtual std::string Info() const
	{
		return "EigenfrequencyResponseFunction";
	}

	/// Print information about this object.
	virtual void PrintInfo(std::ostream &rOStream) const
	{
		rOStream << "EigenfrequencyResponseFunction";
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

		const VariableDenseVectorType& rEIGENVALUE_VECTOR =
            KratosComponents<VariableDenseVectorType>::Get("EIGENVALUE_VECTOR");
		const double eigenvalue = (mrModelPart.GetProcessInfo()[rEIGENVALUE_VECTOR])[mTracedEigenValue - 1];

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

			const VariableDenseMatrixType& rEIGENVECTOR_MATRIX =
           	      KratosComponents<VariableDenseMatrixType>::Get("EIGENVECTOR_MATRIX");

			// Get eigenvector of element
			int k = 0;
			const int NumNodeDofs = num_dofs_element/elem_i.GetGeometry().size();
			for (auto& node_i : elem_i.GetGeometry())
			{
				Matrix& rNodeEigenvectors = node_i.GetValue(rEIGENVECTOR_MATRIX);

				for (int i = 0; i < NumNodeDofs; i++)
                    eigenvector_of_element(i+NumNodeDofs*k) = rNodeEigenvectors((mTracedEigenValue-1),i);

				k++;
			}

			// Semi-analytic computation of partial derivative of state equation w.r.t. node coordinates
			for (auto& node_i : elem_i.GetGeometry())
			{

				array_3d gradient_contribution(3, 0.0);
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

	virtual void ConsiderDiscretization()
	{
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
	//      EigenfrequencyResponseFunction& operator=(SEigenfrequencyResponseFunction const& rOther);

	/// Copy constructor.
	//      EigenfrequencyResponseFunction(EigenfrequencyResponseFunction const& rOther);

	///@}

}; // Class EigenfrequencyResponseFunction

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

} // namespace Kratos.

#endif // EIGENFRQUENCY_RESPONSE_FUNCTION_H

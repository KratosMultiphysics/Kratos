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

//template<class TDenseSpace>

class EigenfrequencyResponseFunction : ResponseFunction
{
public:
	///@name Type Definitions
	///@{

	// TODO solve this via template or how to get this from Eigensolverstrategy 
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;

	typedef array_1d<double, 3> array_3d;
    typedef typename LocalSpaceType::VectorType DenseVectorType;
	typedef typename LocalSpaceType::MatrixType DenseMatrixType;
	typedef Variable<DenseVectorType> VariableDenseVectorType;
	typedef Variable<DenseMatrixType> VariableDenseMatrixType;
	

	/// Pointer definition of EigenfrequencyResponseFunction
	KRATOS_CLASS_POINTER_DEFINITION(EigenfrequencyResponseFunction);

	///@}
	///@name Life Cycle
	///@{

	/// Default constructor.
	EigenfrequencyResponseFunction(ModelPart& model_part, Parameters& responseSettings)
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
		else
			KRATOS_THROW_ERROR(std::invalid_argument, "Specified gradient_mode not recognized. The only option is: semi_analytic. Specified gradient_mode: ", gradientMode);

		// Get number of eigenfrequency for which the structure has to be optimized
		m_traced_eigenvalue = responseSettings["traced_eigenfrequency"].GetInt();


		// Initialize member variables to NULL
		m_initial_value = 0.0;
		m_initial_value_defined = false;
		m_eigenvalue = 0.0;

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
	void initialize()
	{
		//not needed because only semi-analytical sensitivity analysis is implemented yet
	}

	// --------------------------------------------------------------------------
	void calculate_value()
	{
		KRATOS_TRY;

		m_eigenvalue = 0.0;

		const VariableDenseVectorType& rEIGENVALUE_VECTOR =
            KratosComponents<VariableDenseVectorType>::Get("EIGENVALUE_VECTOR");

		int num_of_computed_eigenvalues = (mr_model_part.GetProcessInfo()[rEIGENVALUE_VECTOR]).size();

		if(num_of_computed_eigenvalues < m_traced_eigenvalue)
			KRATOS_THROW_ERROR(std::runtime_error, "The chosen eigenvalue was not solved by the eigenvalue analysis!", "");

		m_eigenvalue = 	(mr_model_part.GetProcessInfo()[rEIGENVALUE_VECTOR])[m_traced_eigenvalue - 1]; 

		// Change sign of response: only maximization makes sense in case of eigenfrequency optimization		
		m_eigenvalue *= (-1.0);

		std::cout << ("I use eigenvalue: ") << m_eigenvalue << std::endl;

		// Set initial value if not done yet
		if(!m_initial_value_defined)
		{
			m_initial_value = m_eigenvalue;
			m_initial_value_defined = true;
		}


		KRATOS_CATCH("");
	}

	// --------------------------------------------------------------------------
	void calculate_gradient()
	{
		KRATOS_TRY;

		// Formula computed in general notation:
		// \frac{dF}{dx} = eigenvector^T\cdot (frac{\partial RHS}{\partial x} - 
		//				   eigenvalue frac{\partial mass_matrix}{\partial x})\cdot eigenvector

		// First gradients are initialized
		array_3d zeros_array(3, 0.0);
		for (ModelPart::NodeIterator node_i = mr_model_part.NodesBegin(); node_i != mr_model_part.NodesEnd(); ++node_i)
			noalias(node_i->FastGetSolutionStepValue(EIGENFREQUENCY_SHAPE_GRADIENT) )= zeros_array;

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

		return m_eigenvalue;

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
			dFdX[node_i->Id()] = node_i->FastGetSolutionStepValue(EIGENFREQUENCY_SHAPE_GRADIENT);

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
		ProcessInfo &CurrentProcessInfo = mr_model_part.GetProcessInfo();

		// Computation of: \frac{dF}{dx} = eigenvector^T\cdot (frac{\partial RHS}{\partial x} - 
		//				                   eigenvalue frac{\partial mass_matrix}{\partial x})\cdot eigenvector
		for (ModelPart::ElementIterator elem_i = mr_model_part.ElementsBegin(); elem_i != mr_model_part.ElementsEnd(); ++elem_i)
		//for (auto& i_elem : mr_model_part.Elements())
		{

	

			Matrix mass_matrix_org;
			Matrix LHS_org;
			Vector eigenvector_of_element = Vector(0);
			Vector dummy;
			elem_i->CalculateMassMatrix(mass_matrix_org, CurrentProcessInfo);
			elem_i->CalculateLocalSystem(LHS_org, dummy ,CurrentProcessInfo);
		
			// Get size of element eigenvector and initialize eigenvector
			int num_dofs_element = mass_matrix_org.size1();
			eigenvector_of_element.resize(num_dofs_element,false);

			const VariableDenseMatrixType& rEIGENVECTOR_MATRIX =
           	      KratosComponents<VariableDenseMatrixType>::Get("EIGENVECTOR_MATRIX");
		
			int k = 0;
			// Get eigenvector of element
			for (ModelPart::NodeIterator node_i = elem_i->GetGeometry().begin(); node_i != elem_i->GetGeometry().end(); ++node_i)
			{
				Matrix& rNodeEigenvectors = node_i->GetValue(rEIGENVECTOR_MATRIX);

				ModelPart::NodeType::DofsContainerType& NodeDofs = node_i->GetDofs();
    
				int NumNodeDofs = NodeDofs.size();

				for (int i = 0; i < NumNodeDofs; i++)
                {
                    eigenvector_of_element(i+NumNodeDofs*k) = rNodeEigenvectors((m_traced_eigenvalue-1),i);
                }
				k++;
			}
			
			// Semi-analytic computation of partial derivative of state equation w.r.t. node coordinates
			for (ModelPart::NodeIterator node_i = elem_i->GetGeometry().begin(); node_i != elem_i->GetGeometry().end(); ++node_i)
			{

				array_3d gradient_contribution(3, 0.0);
				Matrix perturbed_LHS = Matrix(0,0);
				Matrix perturbed_mass_matrix = Matrix(0,0);
				Vector aux = Vector(0);

				// Derivative of response w.r.t. x-coord ------------------------
				node_i->X0() += mDelta;

				elem_i->CalculateMassMatrix(perturbed_mass_matrix, CurrentProcessInfo);
				perturbed_mass_matrix = (perturbed_mass_matrix - mass_matrix_org) / mDelta;
				perturbed_mass_matrix *= m_eigenvalue;
				
				elem_i->CalculateLocalSystem(perturbed_LHS, dummy ,CurrentProcessInfo);
				perturbed_LHS = (perturbed_LHS - LHS_org) / mDelta;

				perturbed_LHS -= perturbed_mass_matrix;

				aux = prod(perturbed_LHS , eigenvector_of_element);
				gradient_contribution[0] = inner_prod(eigenvector_of_element , aux);

				node_i->X0() -= mDelta;
				// End derivative of response w.r.t. x-coord --------------------

				// Reset pertubed RHS and mass matrix
				perturbed_LHS = Matrix(0,0);
				perturbed_mass_matrix = Matrix(0,0);
				aux = Vector(0);

				// Derivative of response w.r.t. y-coord ------------------------
				node_i->Y0() += mDelta;
				
				elem_i->CalculateMassMatrix(perturbed_mass_matrix, CurrentProcessInfo);

				perturbed_mass_matrix = (perturbed_mass_matrix - mass_matrix_org) / mDelta;
				perturbed_mass_matrix *= m_eigenvalue;

				elem_i->CalculateLocalSystem(perturbed_LHS, dummy ,CurrentProcessInfo);

				perturbed_LHS = (perturbed_LHS - LHS_org) / mDelta;

				perturbed_LHS -= perturbed_mass_matrix;
				
				aux = prod(perturbed_LHS , eigenvector_of_element);
				gradient_contribution[1] = inner_prod(eigenvector_of_element , aux);
				node_i->Y0() -= mDelta;
				// End derivative of response w.r.t. y-coord --------------------

				// Reset pertubed RHS and mass matrix
				perturbed_LHS = Matrix(0,0);
				perturbed_mass_matrix= Matrix(0,0);
				aux = Vector(0);

				// Derivative of response w.r.t. z-coord ------------------------
				node_i->Z0() += mDelta;

				elem_i->CalculateMassMatrix(perturbed_mass_matrix, CurrentProcessInfo);
				perturbed_mass_matrix = (perturbed_mass_matrix - mass_matrix_org) / mDelta;
				perturbed_mass_matrix *= m_eigenvalue;

				elem_i->CalculateLocalSystem(perturbed_LHS, dummy ,CurrentProcessInfo);
				perturbed_LHS = (perturbed_LHS - LHS_org) / mDelta;

				perturbed_LHS -= perturbed_mass_matrix;
				
				aux = prod(perturbed_LHS , eigenvector_of_element);
				gradient_contribution[2] =  inner_prod(eigenvector_of_element , aux);
				node_i->Z0() -= mDelta;
				// End derivative of response w.r.t. z-coord --------------------

				// Change sign of gradient: only maximization makes sense in case of eigenfrequency optimization	
				gradient_contribution[0] *= (-1.0);
				gradient_contribution[1] *= (-1.0);
				gradient_contribution[2] *= (-1.0);

				// Assemble sensitivity to node
				noalias(node_i->FastGetSolutionStepValue(EIGENFREQUENCY_SHAPE_GRADIENT)) += gradient_contribution;
				
			}
		}
		KRATOS_CATCH("");
	}

	virtual void consider_discretization(){
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
	unsigned int mWeightingMethod;
	double m_eigenvalue; 
	double mDelta;
	double m_initial_value;
	bool m_initial_value_defined;
	int m_traced_eigenvalue;


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

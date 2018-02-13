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

#ifndef EIGENFREQUENCY_RESPONSE_FUNCTION_LIN_SCAL_H
#define EIGENFREQUENCY_RESPONSE_FUNCTION_LIN_SCAL_H

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

class EigenfrequencyResponseFunctionLinScal : ResponseFunction
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


	/// Pointer definition of EigenfrequencyResponseFunctionLinScal
	KRATOS_CLASS_POINTER_DEFINITION(EigenfrequencyResponseFunctionLinScal);

	///@}
	///@name Life Cycle
	///@{

	/// Default constructor.
	EigenfrequencyResponseFunctionLinScal(ModelPart& model_part, Parameters& responseSettings)
	: mr_model_part(model_part)
	{
		// Set gradient mode
		std::string gradientMode = responseSettings["gradient_mode"].GetString();

		// Mode 1: semi-analytic sensitivities
		if (gradientMode.compare("semi_analytic") == 0)
		{
			mGradientMode = 1;
			mDelta = responseSettings["step_size"].GetDouble();
		}
		else
			KRATOS_THROW_ERROR(std::invalid_argument, "Specified gradient_mode not recognized. The only option is: semi_analytic. Specified gradient_mode: ", gradientMode);


        // Get array of numbers of the eigenfrequencies which have to be traced by this response function
		m_num_eigenvalues =  responseSettings["traced_eigenfrequency"].size();
		m_vector_ev.resize(m_num_eigenvalues,false);

		for(int i = 0; i < m_num_eigenvalues; i++)
			m_vector_ev[i] = responseSettings["traced_eigenfrequency"][i].GetInt();

		// Ask for weighting factors and check their number
		m_num_weight_fac =  responseSettings["weighting_factors"].size();
		if(m_num_weight_fac != m_num_eigenvalues)
			KRATOS_THROW_ERROR(std::logic_error, "The number of chosen eigenvalues does not fit to the number of weighting factors!", "");
		m_vector_weight_fac.resize(m_num_weight_fac,false);

		// Read weighting factors and sum them up
		double test_sum_weight_facs = 0.0;
		for(int i = 0; i < m_num_weight_fac; i++)
		{
			m_vector_weight_fac[i] = responseSettings["weighting_factors"][i].GetDouble();
			test_sum_weight_facs += m_vector_weight_fac[i];
		}

		// Check the weighting factors and modify them for the case that their sum is unequal to one
		if(test_sum_weight_facs != 1.0)
		{
			for(int i = 0; i < m_num_eigenvalues; i++)
				m_vector_weight_fac[i] /= test_sum_weight_facs;

			std::cout << "> The sum of the chosen weighting factors is unequal to one. A scaling process was executed for them!" << std::endl;
		}

		// Initialize member variables to NULL
		m_initial_value = 0.0;
		m_initial_value_defined = false;
		m_resp_function_value = 0.0;
	}

	/// Destructor.
	virtual ~EigenfrequencyResponseFunctionLinScal()
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

		m_resp_function_value = 0.0;

		// Compute response function by weighting with linear scalarization
		for(int i = 0; i < m_num_eigenvalues; i++)
			m_resp_function_value += m_vector_weight_fac[i] * get_single_eigenvalue(m_vector_ev[i]);

		// Set initial value if not done yet
		if(!m_initial_value_defined)
		{
			m_initial_value = m_resp_function_value;
			m_initial_value_defined = true;
		}

		KRATOS_CATCH("");
	}

	// --------------------------------------------------------------------------
	double get_single_eigenvalue(int id_eigenvalue)
	{
		KRATOS_TRY;

		const VariableDenseVectorType& rEIGENVALUE_VECTOR =
            KratosComponents<VariableDenseVectorType>::Get("EIGENVALUE_VECTOR");

		int num_of_computed_eigenvalues = (mr_model_part.GetProcessInfo()[rEIGENVALUE_VECTOR]).size();

		if(num_of_computed_eigenvalues < id_eigenvalue)
			KRATOS_THROW_ERROR(std::runtime_error, "The chosen eigenvalue was not solved by the eigenvalue analysis!", "");

		return (mr_model_part.GetProcessInfo()[rEIGENVALUE_VECTOR])[id_eigenvalue-1];

		KRATOS_CATCH("");
	}

	// --------------------------------------------------------------------------
	Vector get_eigenvector_of_element(ModelPart::ElementType& traced_element, int id_eigenvalue, int size_of_eigenvector)
	{
		KRATOS_TRY;

		Vector eigenvector_of_element;
		eigenvector_of_element.resize(size_of_eigenvector,false);

		const VariableDenseMatrixType& rEIGENVECTOR_MATRIX =
           	  KratosComponents<VariableDenseMatrixType>::Get("EIGENVECTOR_MATRIX");

		// Get eigenvector of element
		int k = 0;
		const int NumNodeDofs = size_of_eigenvector/traced_element.GetGeometry().size();
		for (auto& node_i : traced_element.GetGeometry())
		{
			Matrix& rNodeEigenvectors = node_i.GetValue(rEIGENVECTOR_MATRIX);

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
		array_3d zeros_array(3, 0.0);
		for (auto& node_i : mr_model_part.Nodes())
			noalias(node_i.FastGetSolutionStepValue(EIGENFREQUENCY_SHAPE_GRADIENT) ) = zeros_array;

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
	double GetInitialValue()
	{
		KRATOS_TRY;

		if(!m_initial_value_defined)
			KRATOS_THROW_ERROR(std::logi:error, "Initial value not yet defined! First compute it by calling \"CalculateValue()\"", m_initial_value_defined);

		return m_initial_value;

		KRATOS_CATCH("");
	}

	// --------------------------------------------------------------------------
	double GetValue()
	{
		KRATOS_TRY;

		return m_resp_function_value;

		KRATOS_CATCH("");
	}

	// --------------------------------------------------------------------------
	boost::python::dict GetGradient()
	{
		KRATOS_TRY;

		// Dictionary to store all sensitivities along with Ids of corresponding nodes
		boost::python::dict dFdX;

		// Fill dictionary with gradient information
		for (auto& node_i : mr_model_part.Nodes())
			dFdX[node_i.Id()] = node_i.FastGetSolutionStepValue(EIGENFREQUENCY_SHAPE_GRADIENT);

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
		return "EigenfrequencyResponseFunctionLinScal";
	}

	/// Print information about this object.
	virtual void PrintInfo(std::ostream &rOStream) const
	{
		rOStream << "EigenfrequencyResponseFunctionLinScal";
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

		for (auto& elem_i : mr_model_part.Elements())
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
				array_3d gradient_contribution(3, 0.0);
				Matrix perturbed_LHS = Matrix(0,0);
				Matrix perturbed_mass_matrix = Matrix(0,0);

				// Derivative of response w.r.t. x-coord ------------------------
				node_i.X0() += mDelta;

				elem_i.CalculateMassMatrix(perturbed_mass_matrix, CurrentProcessInfo);
				perturbed_mass_matrix = (perturbed_mass_matrix - mass_matrix_org) / mDelta;

				elem_i.CalculateLocalSystem(perturbed_LHS, dummy ,CurrentProcessInfo);
				perturbed_LHS = (perturbed_LHS - LHS_org) / mDelta;

				// Loop over eigenvalues
				for(int i = 0; i < m_num_eigenvalues; i++)
				{
					traced_eigenvalue = get_single_eigenvalue(m_vector_ev[i]);
					eigenvector_of_element = get_eigenvector_of_element(elem_i, m_vector_ev[i], num_dofs_element);

					aux_matrix = perturbed_LHS;
					aux_matrix -= (perturbed_mass_matrix * traced_eigenvalue);
					aux_vector = prod(aux_matrix , eigenvector_of_element);
					gradient_contribution[0] += inner_prod(eigenvector_of_element , aux_vector) * m_vector_weight_fac[i];

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
				for(int i = 0; i < m_num_eigenvalues; i++)
				{
					traced_eigenvalue = get_single_eigenvalue(m_vector_ev[i]);
					eigenvector_of_element = get_eigenvector_of_element(elem_i, m_vector_ev[i], num_dofs_element);

					aux_matrix = perturbed_LHS;
					aux_matrix -= (perturbed_mass_matrix * traced_eigenvalue);
					aux_vector = prod(aux_matrix , eigenvector_of_element);
					gradient_contribution[1] += inner_prod(eigenvector_of_element , aux_vector) * m_vector_weight_fac[i];

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
				for(int i = 0; i < m_num_eigenvalues; i++)
				{
					traced_eigenvalue = get_single_eigenvalue(m_vector_ev[i]);
					eigenvector_of_element = get_eigenvector_of_element(elem_i, m_vector_ev[i], num_dofs_element);

					aux_matrix = perturbed_LHS;
					aux_matrix -= (perturbed_mass_matrix * traced_eigenvalue);
					aux_vector = prod(aux_matrix , eigenvector_of_element);
					gradient_contribution[2] += inner_prod(eigenvector_of_element , aux_vector) * m_vector_weight_fac[i];

					eigenvector_of_element = Vector(0);
					aux_matrix = Matrix(0,0);
					aux_vector = Vector(0);
				}

				node_i.Z0() -= mDelta;
				// End derivative of response w.r.t. z-coord --------------------

				// Assemble sensitivity to node
				noalias(node_i.FastGetSolutionStepValue(EIGENFREQUENCY_SHAPE_GRADIENT)) += gradient_contribution;

			}// End loop over nodes of element

		}// End loop over elements
		KRATOS_CATCH("");
	}

	virtual void ConsiderDiscretization(){
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
	double m_resp_function_value;
	double mDelta;
	double m_initial_value;
	bool m_initial_value_defined;
	int m_num_eigenvalues;
	std::vector<int> m_vector_ev;
	int m_num_weight_fac;
	std::vector<double> m_vector_weight_fac;

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
	//      EigenfrequencyResponseFunctionLinScal& operator=(EigenfrequencyResponseFunctionLinScal const& rOther);

	/// Copy constructor.
	//      EigenfrequencyResponseFunctionLinScal(EigenfrequencyResponseFunctionLinScal const& rOther);

	///@}

}; // Class EigenfrequencyResponseFunctionLinScal

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

} // namespace Kratos.

#endif // EIGENFRQUENCY_RESPONSE_FUNCTION_H

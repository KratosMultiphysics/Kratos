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

#ifndef EIGENFREQUENCY_RESPONSE_FUNCTION_KS_H
#define EIGENFREQUENCY_RESPONSE_FUNCTION_KS_H

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

class EigenfrequencyResponseFunctionKS : ResponseFunction
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
	

	/// Pointer definition of EigenfrequencyResponseFunctionKS
	KRATOS_CLASS_POINTER_DEFINITION(EigenfrequencyResponseFunctionKS);

	///@}
	///@name Life Cycle
	///@{

	/// Default constructor.
	EigenfrequencyResponseFunctionKS(ModelPart& model_part, Parameters& responseSettings)
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


        // Get array of numbers of eigenfrequencies which have to be traced by this response function
		m_num_eigenvalues =  responseSettings["traced_eigenfrequency"].size();
		m_vector_ev.resize(m_num_eigenvalues,false);
		
		for(int i = 0; i < m_num_eigenvalues; i++)
			m_vector_ev[i] = responseSettings["traced_eigenfrequency"][i].GetInt();
		
		// Set KS-Parameter
		m_KS_parameter= responseSettings["KS_parameter"].GetDouble();
	
		// Initialize member variables to NULL
		m_initial_value = 0.0;
		m_initial_value_defined = false;
		m_resp_function_value = 0.0;
		m_scal_fac_computed = false;


  		// Check validity of KS-Parameter
  		if (m_KS_parameter <= 0)
  		{
			KRATOS_THROW_ERROR(std::invalid_argument, "KS-parameter has to be positive. Your choice is: ", m_KS_parameter);
  		}
  		else if (m_KS_parameter > 100.0 || m_KS_parameter < 1.0)
  		{
			KRATOS_THROW_ERROR(std::invalid_argument, "KS-parameter is out of the recommended range of [1, 100]. Your choice is: ", m_KS_parameter);
  		}


	}

	/// Destructor.
	virtual ~EigenfrequencyResponseFunctionKS()
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

		m_resp_function_value = 0.0;

		// Compute response function by weighting with Kreisselmeier-Steinhauser
		for(int i = 0; i < m_num_eigenvalues; i++)
			m_resp_function_value += compute_single_ks_value(i);

		m_resp_function_value  = log( m_resp_function_value );	
		m_resp_function_value /= m_KS_parameter;

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

		double current_eigenvalue = 0.0;

		const VariableDenseVectorType& rEIGENVALUE_VECTOR =
            KratosComponents<VariableDenseVectorType>::Get("EIGENVALUE_VECTOR");

		int num_of_computed_eigenvalues = (mr_model_part.GetProcessInfo()[rEIGENVALUE_VECTOR]).size();

		if(num_of_computed_eigenvalues < id_eigenvalue)
			KRATOS_THROW_ERROR(std::runtime_error, "The chosen eigenvalue was not solved by the eigenvalue analysis!", "");

		current_eigenvalue = (mr_model_part.GetProcessInfo()[rEIGENVALUE_VECTOR])[id_eigenvalue-1]; 

		return current_eigenvalue;
		
		KRATOS_CATCH("");

	}

	// --------------------------------------------------------------------------
	double compute_single_ks_value(int id_eigenvalue_in_vector)
	{

		KRATOS_TRY;

		double lower_bound_mapping = 0.1;
    	double upper_bound_mapping = 0.5;
    	double shift= 80.0;

		// Perform tranformation of eigenvalues into [0.1 , 0.5] (only performed once)
    	if (!m_scal_fac_computed)
		{

        	m_scaling_factors_eigenvalues.resize(m_num_eigenvalues);
        	double smallest_ev = get_single_eigenvalue(m_vector_ev[0]);
        	double largest_ev =  get_single_eigenvalue(m_vector_ev[m_num_eigenvalues - 1]);
        	double lambda_scaled;
			double current_eigenvalue = 0.0;

        	if (m_num_eigenvalues == 1)
            	m_scaling_factors_eigenvalues[0] = upper_bound_mapping / smallest_ev;
        	else
        	{
          		// Loop over traced eigenvalues
          		for (int i = 0; i < m_num_eigenvalues; ++i)
          		{
					current_eigenvalue = get_single_eigenvalue(m_vector_ev[i]);
					lambda_scaled = lower_bound_mapping  + (current_eigenvalue - smallest_ev) *
                	(upper_bound_mapping - lower_bound_mapping ) / (largest_ev - smallest_ev);
					m_scaling_factors_eigenvalues[i] = lambda_scaled / current_eigenvalue;
          		}
			}
        
        	// Change boolean
        	m_scal_fac_computed = true;
   		}

		double traced_eigenvalue = get_single_eigenvalue(m_vector_ev[id_eigenvalue_in_vector]);

    	// Calculate KS-term for one eigenvalue
   		double KS_value = (std::exp( -1.0 * m_KS_parameter *
       					  std::sqrt(traced_eigenvalue * m_scaling_factors_eigenvalues[id_eigenvalue_in_vector]) + shift));
		//double KS_value = std::exp( -1.0 * m_KS_parameter * current_eigenvalue * m_scaling_factors_eigenvalues[id_eigenvalue_in_vector]  + shift  );

    	return KS_value;

		KRATOS_CATCH("");
	}

	// --------------------------------------------------------------------------
	Vector get_eigenvector_of_element(ModelPart::ElementIterator traced_element, int id_eigenvalue, int size_of_eigenvector)
	{

		KRATOS_TRY;

		Vector eigenvector_of_element;
		eigenvector_of_element.resize(size_of_eigenvector,false);

		const VariableDenseMatrixType& rEIGENVECTOR_MATRIX =
           	  KratosComponents<VariableDenseMatrixType>::Get("EIGENVECTOR_MATRIX");
		
		int k = 0;
		for (ModelPart::NodeIterator node_i = traced_element->GetGeometry().begin(); node_i != traced_element->GetGeometry().end(); ++node_i)
		{
			Matrix& rNodeEigenvectors = node_i->GetValue(rEIGENVECTOR_MATRIX);

			ModelPart::NodeType::DofsContainerType& NodeDofs = node_i->GetDofs();

			int NumNodeDofs = NodeDofs.size();

			for (int i = 0; i < NumNodeDofs; i++)
            {
                eigenvector_of_element(i+NumNodeDofs*k) = rNodeEigenvectors((id_eigenvalue-1),i);
            }
			k++;
		}

		return eigenvector_of_element;

		KRATOS_CATCH("");

	}

	// --------------------------------------------------------------------------
	void calculate_gradient()
	{
		KRATOS_TRY;

		// First gradients are initialized
		array_3d zeros_array(3, 0.0);
		for (ModelPart::NodeIterator node_i = mr_model_part.NodesBegin(); node_i != mr_model_part.NodesEnd(); ++node_i)
			noalias(node_i->FastGetSolutionStepValue(EIGENFREQUENCY_SHAPE_GRADIENT) ) = zeros_array;

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

		return m_resp_function_value;

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
		return "EigenfrequencyResponseFunctionKS";
	}

	/// Print information about this object.
	virtual void PrintInfo(std::ostream &rOStream) const
	{
		rOStream << "EigenfrequencyResponseFunctionKS";
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

		for (ModelPart::ElementIterator elem_i = mr_model_part.ElementsBegin(); elem_i != mr_model_part.ElementsEnd(); ++elem_i)
		{
			Matrix mass_matrix_org;
			Matrix LHS_org;
			Vector eigenvector_of_element = Vector(0);
			Vector dummy;
			Matrix aux_matrix = Matrix(0,0);
			Vector aux_vector = Vector(0);
			elem_i->CalculateMassMatrix(mass_matrix_org, CurrentProcessInfo);
			elem_i->CalculateLocalSystem(LHS_org, dummy ,CurrentProcessInfo);
			double traced_eigenvalue = 0.0;
		
			// Get size of element eigenvector and initialize eigenvector
			int num_dofs_element = mass_matrix_org.size1();
			eigenvector_of_element.resize(num_dofs_element,false);

			Vector ks_values = Vector(0);
			ks_values.resize(m_num_eigenvalues,false);
			double ks_values_sum = 0.0;

			for(int i = 0; i < m_num_eigenvalues; i++)
			{
				ks_values[i] = compute_single_ks_value(i);
				ks_values_sum += ks_values[i]; 
			}

			// Semi-analytic computation of partial derivative of state equation w.r.t. node coordinates
			for (ModelPart::NodeIterator node_i = elem_i->GetGeometry().begin(); node_i != elem_i->GetGeometry().end(); ++node_i)
			{
				array_3d gradient_contribution(3, 0.0);
				Matrix perturbed_LHS = Matrix(0,0);
				Matrix perturbed_mass_matrix = Matrix(0,0);
				
				// Derivative of response w.r.t. x-coord ------------------------
				node_i->X0() += mDelta;

				elem_i->CalculateMassMatrix(perturbed_mass_matrix, CurrentProcessInfo);
				perturbed_mass_matrix = (perturbed_mass_matrix - mass_matrix_org) / mDelta;
				
				elem_i->CalculateLocalSystem(perturbed_LHS, dummy ,CurrentProcessInfo);
				perturbed_LHS = (perturbed_LHS - LHS_org) / mDelta;

				// Loop over eigenvalues
				for(int i = 0; i < m_num_eigenvalues; i++)
				{
					traced_eigenvalue = get_single_eigenvalue(m_vector_ev[i]);
					eigenvector_of_element = get_eigenvector_of_element(elem_i, m_vector_ev[i], num_dofs_element);

					aux_matrix = perturbed_LHS;
					aux_matrix -= (perturbed_mass_matrix * traced_eigenvalue);
					aux_vector = prod(aux_matrix , eigenvector_of_element);
					
					gradient_contribution[0] += (inner_prod(eigenvector_of_element , aux_vector) *
												ks_values[i] * m_scaling_factors_eigenvalues[i] / 
												(2.0 * std::sqrt(m_scaling_factors_eigenvalues[i] * traced_eigenvalue)));							

					eigenvector_of_element = Vector(0);
					aux_matrix = Matrix(0,0);
					aux_vector = Vector(0);

				}
								
				gradient_contribution[0] /= ks_values_sum;
				gradient_contribution[0] *= -1.0;
			
				node_i->X0() -= mDelta;

				// End derivative of response w.r.t. x-coord --------------------

				// Reset pertubed RHS and mass matrix
				perturbed_LHS = Matrix(0,0);
				perturbed_mass_matrix = Matrix(0,0);
		

				// Derivative of response w.r.t. y-coord ------------------------
				node_i->Y0() += mDelta;

				elem_i->CalculateMassMatrix(perturbed_mass_matrix, CurrentProcessInfo);
				perturbed_mass_matrix = (perturbed_mass_matrix - mass_matrix_org) / mDelta;

				elem_i->CalculateLocalSystem(perturbed_LHS, dummy ,CurrentProcessInfo);
				perturbed_LHS = (perturbed_LHS - LHS_org) / mDelta;

				// Loop over eigenvalues
				for(int i = 0; i < m_num_eigenvalues; i++)
				{
					traced_eigenvalue = get_single_eigenvalue(m_vector_ev[i]);
					eigenvector_of_element = get_eigenvector_of_element(elem_i, m_vector_ev[i], num_dofs_element);

					aux_matrix = perturbed_LHS;
					aux_matrix -= (perturbed_mass_matrix * traced_eigenvalue);
					aux_vector = prod(aux_matrix , eigenvector_of_element);
					gradient_contribution[1] += (inner_prod(eigenvector_of_element , aux_vector) *
												ks_values[i] * m_scaling_factors_eigenvalues[i] / 
												(2.0 * std::sqrt(m_scaling_factors_eigenvalues[i] * traced_eigenvalue)));
					
					eigenvector_of_element = Vector(0);
				    aux_matrix = Matrix(0,0);
					aux_vector = Vector(0);
				}

				gradient_contribution[1] /= ks_values_sum;
				gradient_contribution[1] *= -1.0;
		
				node_i->Y0() -= mDelta;
				// End derivative of response w.r.t. y-coord --------------------

				// Reset pertubed RHS and mass matrix
				perturbed_LHS = Matrix(0,0);
				perturbed_mass_matrix= Matrix(0,0);

				// Derivative of response w.r.t. z-coord ------------------------
				node_i->Z0() += mDelta;

				elem_i->CalculateMassMatrix(perturbed_mass_matrix, CurrentProcessInfo);
				perturbed_mass_matrix = (perturbed_mass_matrix - mass_matrix_org) / mDelta;


				elem_i->CalculateLocalSystem(perturbed_LHS, dummy ,CurrentProcessInfo);
				perturbed_LHS = (perturbed_LHS - LHS_org) / mDelta;

				// Loop over eigenvalues
				for(int i = 0; i < m_num_eigenvalues; i++)
				{
					traced_eigenvalue = get_single_eigenvalue(m_vector_ev[i]);
					eigenvector_of_element = get_eigenvector_of_element(elem_i, m_vector_ev[i], num_dofs_element);

					aux_matrix = perturbed_LHS;
					aux_matrix -= (perturbed_mass_matrix * traced_eigenvalue);
					aux_vector = prod(aux_matrix , eigenvector_of_element);
					gradient_contribution[2] += (inner_prod(eigenvector_of_element , aux_vector) *
												ks_values[i] * m_scaling_factors_eigenvalues[i] / 
												(2.0 * std::sqrt(m_scaling_factors_eigenvalues[i] * traced_eigenvalue)));

					eigenvector_of_element = Vector(0);
					aux_matrix = Matrix(0,0);
					aux_vector = Vector(0);
				}

				gradient_contribution[2] /= ks_values_sum;
				gradient_contribution[2] *= -1.0;
	
				node_i->Z0() -= mDelta;
				// End derivative of response w.r.t. z-coord --------------------

				// Assemble sensitivity to node
				noalias(node_i->FastGetSolutionStepValue(EIGENFREQUENCY_SHAPE_GRADIENT)) += gradient_contribution;

			}// End loop over nodes of element
	
		}// End loop over elements
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
	double m_resp_function_value;
	double mDelta;
	double m_initial_value;
	bool m_initial_value_defined;
	int m_num_eigenvalues;
	std::vector<int> m_vector_ev;
	double m_KS_parameter;
	Vector m_scaling_factors_eigenvalues;
	bool m_scal_fac_computed;

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
	//      EigenfrequencyResponseFunctionKS& operator=(EigenfrequencyResponseFunctionKS const& rOther);

	/// Copy constructor.
	//      EigenfrequencyResponseFunctionKS(EigenfrequencyResponseFunctionKS const& rOther);

	///@}

}; // Class EigenfrequencyResponseFunctionKS

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

} // namespace Kratos.

#endif // EIGENFRQUENCY_RESPONSE_FUNCTION_KS_H

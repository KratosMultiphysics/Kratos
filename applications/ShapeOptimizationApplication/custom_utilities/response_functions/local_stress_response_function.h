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

#ifndef LOCAL_STRESS_RESPONSE_FUNCTION_H
#define LOCAL_STRESS_RESPONSE_FUNCTION_H

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
#include "structural_response_function.h"


#include "includes/kratos_parameters.h"
#include "includes/ublas_interface.h"
#include "utilities/openmp_utils.h"

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

class LocalStressResponseFunction : StructuralResponseFunction
{
public:
	///@name Type Definitions
	///@{

	typedef StructuralResponseFunction BaseType;
	typedef array_1d<double, 3> array_3d;

	

	/// Pointer definition of LocalStressResponseFunction
	KRATOS_CLASS_POINTER_DEFINITION(LocalStressResponseFunction);

	///@}
	///@name Life Cycle
	///@{

	/// Default constructor.
	LocalStressResponseFunction(ModelPart& model_part, Parameters& responseSettings)
	: StructuralResponseFunction(model_part, responseSettings)
	{
		ModelPart& r_model_part = this->GetModelPart();

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
	

		//get traced element
		m_id_of_traced_element = responseSettings["traced_element"].GetInt();
		m_traced_pElement = r_model_part.pGetElement(m_id_of_traced_element);

		//give stress location to traced element
		m_id_of_location = responseSettings["stress_location"].GetInt();
		m_traced_pElement->SetValue(LOCATION_OF_TRACED_STRESS, m_id_of_location);

		//tell traced element the stress type 
		m_traced_stress_type = responseSettings["stress_type"].GetString();
		m_traced_pElement->SetValue(TRACED_STRESS_TYPE, m_traced_stress_type);

		m_stress_treatment = responseSettings["stress_treatment"].GetString();
		m_traced_pElement->SetValue(STRESS_TREATMENT, m_stress_treatment);

		// Initialize member variables to NULL
		m_initial_value = 0.0;
		m_initial_value_defined = false;
		m_stress_value = 0.0;	

	}

	/// Destructor.
	virtual ~LocalStressResponseFunction()
	{
	}

	///@}
	///@name Operators
	///@{

	///@}
	///@name Operations
	///@{

	void Initialize() override
	{
		KRATOS_TRY;

		BaseType::Initialize();

		KRATOS_CATCH("");
	}

	// ==============================================================================
	void CalculateValue()
	{
		KRATOS_TRY;

		ModelPart& r_model_part = this->GetModelPart();

		// Working variables
		ProcessInfo &CurrentProcessInfo = r_model_part.GetProcessInfo();

		m_traced_pElement->Calculate(STRESS_VALUE, m_stress_value, CurrentProcessInfo);

		//just for testing the computation of the adjoint load--> erase this later!!!!!!!!!!!!!!!!!
		std::cout << ("Response Function value= ") << m_stress_value << std::endl;

		Vector adjoint_load;
		Vector zero_adjoint_load;
		m_traced_pElement->Calculate(ADJOINT_LOAD, adjoint_load, CurrentProcessInfo);

		m_traced_pElement->Calculate(ZERO_ADJOINT_LOAD, zero_adjoint_load, CurrentProcessInfo);

		int size_load = adjoint_load.size();

		for(int i = 0; i < size_load; i++)
		{
			//std::cout << ("adjoint_load = ") << adjoint_load[i] << std::endl;	
			std::cout << ("zero_adjoint_load = ") << zero_adjoint_load[i] << std::endl;	
		}

		this->UpdateSensitivities();
		//-------------------------------------------------------------------------!!!!!!!!!!!!!!!!!!

		// Set initial value if not done yet
		if(!m_initial_value_defined)
		{
			m_initial_value = m_stress_value;
			m_initial_value_defined = true;
		}


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

		return m_stress_value;

		KRATOS_CATCH("");
	}

	// --------------------------------------------------------------------------
	/*boost::python::dict get_gradient()
	{
		KRATOS_TRY;

		// Dictionary to store all sensitivities along with Ids of corresponding nodes
		boost::python::dict dFdX;

		ModelPart& r_model_part = this->GetModelPart();

		// Fill dictionary with gradient information
		for (ModelPart::NodeIterator node_i = r_model_part.NodesBegin(); node_i != r_model_part.NodesEnd(); ++node_i)
			dFdX[node_i->Id()] = node_i->FastGetSolutionStepValue(LOCAL_STRESS_GRADIENT);

		return dFdX;

		KRATOS_CATCH("");
	}*/

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
		return "LocalStressResponseFunction";
	}

	/// Print information about this object.
	virtual void PrintInfo(std::ostream &rOStream) const
	{
		rOStream << "LocalStressResponseFunction";
	}

	/// Print object's data.
	virtual void PrintData(std::ostream &rOStream) const
	{
	}

	///@}
	///@name Friends
	///@{

	///@}

	// ==============================================================================
	void CalculateGradient(const Element& rElem, const Matrix& rLHS, 
								Vector& rOutput, ProcessInfo& rProcessInfo) override
	{
		if(rElem.Id() == m_id_of_traced_element)
		{
			//computes adjoint load in global direction. Or is it in local direction required???????????????
			m_traced_pElement->Calculate(ADJOINT_LOAD, rOutput, rProcessInfo);
		}
		else
		{
			// There is only a contribution of the traced element to the adjoint load case
			int num_DOFs_element = rLHS.size1();
			rOutput.resize(num_DOFs_element);
			for(int i = 0; i < num_DOFs_element; i++){ rOutput[i] = 0.0; }
		}

	}

	void UpdateSensitivities() override
	{
		KRATOS_TRY;

		BaseType::UpdateSensitivities();

		KRATOS_CATCH("");
	}


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
	void CalculateSensitivityGradient(Element& rElem,
                                	  const Variable<array_1d<double,3>>& rVariable,
                                      const Matrix& rDerivativesMatrix,
                                      Vector& rRHSContribution,
                                      ProcessInfo& rProcessInfo) //--> auf base class anpassen
    {
      	KRATOS_TRY

      	if (rRHSContribution.size() != rDerivativesMatrix.size1())
          	rRHSContribution.resize(rDerivativesMatrix.size1(), false);

		Vector zero_adjoint_vector;	  
		zero_adjoint_vector  = ZeroVector(rDerivativesMatrix.size1());

		if(rElem.Id() == m_id_of_traced_element)
		{
			rElem.Calculate(ZERO_ADJOINT_LOAD, zero_adjoint_vector, rProcessInfo);
			noalias(rRHSContribution) = prod(rDerivativesMatrix, zero_adjoint_vector);
		}
		else
		{
			 noalias(rRHSContribution) = zero_adjoint_vector;
		}	  

      KRATOS_CATCH("")
	}

	// ==============================================================================
	void CalculateSensitivityGradient(Element& rElem,
                                          const Variable<double>& rVariable,
                                          const Matrix& rDerivativesMatrix,
                                          Vector& rRHSContribution,
                                          ProcessInfo& rProcessInfo) //--> auf base class anpassen
    {
      	KRATOS_TRY

      	if (rRHSContribution.size() != rDerivativesMatrix.size1())
          	rRHSContribution.resize(rDerivativesMatrix.size1(), false);

		Vector zero_adjoint_vector;	  
		zero_adjoint_vector  = ZeroVector(rDerivativesMatrix.size1());

		if(rElem.Id() == m_id_of_traced_element)
		{
			rElem.Calculate(ZERO_ADJOINT_LOAD, zero_adjoint_vector, rProcessInfo);
			noalias(rRHSContribution) = prod(rDerivativesMatrix, zero_adjoint_vector);;
		}
		else
		{
			 noalias(rRHSContribution) = zero_adjoint_vector;
		}	  

       KRATOS_CATCH("")
	}

	// ==============================================================================
	void CalculateSensitivityGradient(Condition& rCond,
                                	  const Variable<array_1d<double,3>>& rVariable,
                                      const Matrix& rDerivativesMatrix,
                                      Vector& rRHSContribution,
                                      ProcessInfo& rProcessInfo)  //--> auf base class anpassen
    {
		KRATOS_TRY;

        // ----> insert code

		KRATOS_CATCH("");
	}

	// ==============================================================================
	void CalculateSensitivityGradient(Condition& rCond,
                                          const Variable<double>& rVariable,
                                          const Matrix& rDerivativesMatrix,
                                          Vector& rRHSContribution,
                                          ProcessInfo& rProcessInfo) //--> auf base class anpassen
    {
		KRATOS_TRY;

         // ----> insert code

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
	unsigned int mGradientMode;
	double m_stress_value; 
	double mDelta;
	double m_initial_value;
	bool m_initial_value_defined;
	unsigned int m_id_of_traced_element;
	int m_id_of_location;
	Element::Pointer m_traced_pElement;
	std::string m_traced_stress_type;
	std::string m_stress_treatment;

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
	//      LocalStressResponseFunction& operator=(SLocalStressResponseFunction const& rOther);

	/// Copy constructor.
	//      LocalStressResponseFunction(LocalStressResponseFunction const& rOther);

	///@}

}; // Class LocalStressResponseFunction

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

} // namespace Kratos.

#endif // LOCAL_STRESS_RESPONSE_FUNCTION_H

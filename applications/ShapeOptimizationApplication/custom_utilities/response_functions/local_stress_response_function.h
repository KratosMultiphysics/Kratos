// ==============================================================================
//  KratosShapeOptimizationApplication
//
//  License:         BSD License
//                   license: ShapeOptimizationApplication/license.txt
//
//  Main authors:    Fusseder Martin   
//                   martin.fusseder@tum.de
//	TODO: Check that this response function  works in a correct way for all conditions
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

class LocalStressResponseFunction : public StructuralResponseFunction
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

	LocalStressResponseFunction(ModelPart& rModelPart, Parameters& rParameters)
	: StructuralResponseFunction(rModelPart, rParameters)
	{
		ModelPart& r_model_part = this->GetModelPart();

		// Set gradient mode
		std::string gradientMode = rParameters["gradient_mode"].GetString();

		// Mode 1: semi-analytic sensitivities
		if (gradientMode.compare("semi_analytic") == 0)
		{
			mGradientMode = 1;
			double delta = rParameters["step_size"].GetDouble();
			mDelta = delta;
		}
		else
			KRATOS_THROW_ERROR(std::invalid_argument, "Specified gradient_mode not recognized. The only option is: semi_analytic. Specified gradient_mode: ", gradientMode);
	

		//get traced element
		m_id_of_traced_element = rParameters["traced_element"].GetInt();
		m_traced_pElement = r_model_part.pGetElement(m_id_of_traced_element);

		//tell traced element the stress type 
		m_traced_stress_type = rParameters["stress_type"].GetString();
		m_traced_pElement->SetValue(TRACED_STRESS_TYPE, m_traced_stress_type);

		// get info how and where to treat the stress
		m_stress_treatment = rParameters["stress_treatment"].GetString();
		if(m_stress_treatment != "mean" && m_stress_treatment != "GP" && m_stress_treatment != "node")
			KRATOS_ERROR << "Chosen option for stress treatmeant is not availible! Chose 'GP','node' or 'mean'!" << std::endl; 
		if(m_stress_treatment == "GP" || m_stress_treatment == "node")
			m_id_of_location = rParameters["stress_location"].GetInt();
		if(m_id_of_location < 1)
			KRATOS_THROW_ERROR(std::invalid_argument, "Chose a 'stress_location' > 0. Specified 'stress_location': ", m_id_of_location);
			
		// Initialize member variables to NULL
		//m_initial_value = 0.0;
		//m_initial_value_defined = false;
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
	double CalculateValue(ModelPart& rModelPart) override
	{
		KRATOS_TRY;

		// Working variables
		ProcessInfo &CurrentProcessInfo = rModelPart.GetProcessInfo();

		Vector element_stress;
		if(m_stress_treatment == "mean" || m_stress_treatment == "GP")
			m_traced_pElement->Calculate(STRESS_ON_GP, element_stress, CurrentProcessInfo);
		else 
			m_traced_pElement->Calculate(STRESS_ON_NODE, element_stress, CurrentProcessInfo);

		int stress_vec_size = element_stress.size();

		if(m_stress_treatment == "mean")
		{
			for(int i = 0; i < stress_vec_size; i++)
				m_stress_value += element_stress[i];

			m_stress_value /= stress_vec_size;	
		}
		else if(m_stress_treatment == "GP")
		{
			if(stress_vec_size >= m_id_of_location)
				m_stress_value = element_stress[m_id_of_location - 1];
			else
				KRATOS_ERROR << "Chosen Gauss-Point is not availible. Chose 'stress_location' between 1 and " << 
								stress_vec_size  << "!"<< std::endl; 
		}
		else if(m_stress_treatment == "node")
		{
			const int num_ele_nodes = m_traced_pElement->GetGeometry().PointsNumber();
			if(num_ele_nodes >= m_id_of_location)
				m_stress_value = element_stress[m_id_of_location - 1];
			else
				KRATOS_ERROR << "Chosen Node is not availible. The element has only " << 
								num_ele_nodes  << " nodes."<< std::endl; 

		}

		// Set initial value if not done yet
		/*if(!m_initial_value_defined)
		{
			m_initial_value = m_stress_value;
			m_initial_value_defined = true;
		}*/

		return m_stress_value;

		KRATOS_CATCH("");
	}
	// --------------------------------------------------------------------------
	/*double GetInitialValue() old function, don't needed for sensitivity analysis
	{
		KRATOS_TRY;

		if(!m_initial_value_defined)
			KRATOS_THROW_ERROR(std::logi:error, "Initial value not yet defined! First compute it by calling \"CalculateValue()\"", m_initial_value_defined);

		return m_initial_value;

		KRATOS_CATCH("");
	}

	// --------------------------------------------------------------------------
	double GetValue() old function, don't needed for sensitivity analysis
	{
		KRATOS_TRY;

		return m_stress_value;

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
	/*virtual std::string Info() const            old function, don't needed for sensitivity analysis
	{
		return "LocalStressResponseFunction";
	}

	/// Print information about this object.
	virtual void PrintInfo(std::ostream &rOStream) const    old function, don't needed for sensitivity analysis
	{
		rOStream << "LocalStressResponseFunction";
	}

	/// Print object's data.
	virtual void PrintData(std::ostream &rOStream) const    old function, don't needed for sensitivity analysis
	{
	}*/

	///@}
	///@name Friends
	///@{

	///@}

	// ==============================================================================
	double GetDisturbanceMeasure() const override
	{ 
		return mDelta; 
	}

	// ==============================================================================
	void CalculateGradient(const Element& rAdjointElem, const Matrix& rAdjointMatrix,
                                   Vector& rResponseGradient,
                                   ProcessInfo& rProcessInfo) override
	{
		rResponseGradient.resize(rAdjointMatrix.size1());
		rResponseGradient.clear();
		
		if(rAdjointElem.Id() == m_id_of_traced_element)
		{
			Matrix stress_displ_deriv;
			if(m_stress_treatment == "mean" || m_stress_treatment == "GP")
				m_traced_pElement->Calculate(STRESS_DISP_DERIV_ON_GP, stress_displ_deriv, rProcessInfo);
			else
				m_traced_pElement->Calculate(STRESS_DISP_DERIV_ON_NODE, stress_displ_deriv, rProcessInfo);	

			int num_of_dofs = stress_displ_deriv.size1();
			int num_of_deriv = stress_displ_deriv.size2();
			double stress_displ_deriv_value = 0.0;

			if(rResponseGradient.size() != stress_displ_deriv.size1())
				KRATOS_ERROR << "Size of stress displacement derivative does not fit!" << std::endl; 

			for (int dof_it = 0 ; dof_it < num_of_dofs; dof_it++)
			{
				if(m_stress_treatment == "mean")
				{
					for(int GP_it = 0; GP_it < num_of_deriv; GP_it++)
						stress_displ_deriv_value += stress_displ_deriv(dof_it, GP_it);

					stress_displ_deriv_value /= num_of_deriv;	
				}
				else if(m_stress_treatment == "GP")
				{
					if(num_of_deriv >= m_id_of_location)
						stress_displ_deriv_value = stress_displ_deriv(dof_it, (m_id_of_location-1));
					else
						KRATOS_ERROR << "Chosen Gauss-Point is not availible. Chose 'stress_location' between 1 and " << 
									num_of_deriv  << "!"<< std::endl; 
				}
				else if(m_stress_treatment == "node")
				{
					if(num_of_deriv >= m_id_of_location)
						stress_displ_deriv_value = stress_displ_deriv(dof_it, (m_id_of_location-1));
					else
						KRATOS_ERROR << "Chosen node is not availible. The element has only " << 
									num_of_deriv  << " nodes."<< std::endl; 

				}
				rResponseGradient[dof_it] = (-1) * stress_displ_deriv_value;
				stress_displ_deriv_value = 0.0;	
			}
		}
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
	void CalculateSensitivityGradient(Element& rAdjointElem, 
                                      const Variable<double>& rVariable,
                                      const Matrix& rDerivativesMatrix,
                                      Vector& rResponseGradient,
                                      ProcessInfo& rProcessInfo) override
    {
      	KRATOS_TRY


		if(rAdjointElem.Id() == m_id_of_traced_element)
		{
			BaseType::CalculateSensitivityGradient(rAdjointElem, rVariable, rDerivativesMatrix, rResponseGradient, rProcessInfo);
		}
		else
		{
			if (rResponseGradient.size() != rDerivativesMatrix.size1())
          			rResponseGradient.resize(rDerivativesMatrix.size1(), false);
			rResponseGradient.clear();
		}
		
      /*	if (rResponseGradient.size() != rDerivativesMatrix.size1())
          	rResponseGradient.resize(rDerivativesMatrix.size1(), false);
		
		Vector zero_adjoint_vector;	  
		zero_adjoint_vector  = ZeroVector(rDerivativesMatrix.size1());

		if(rAdjointElem.Id() == m_id_of_traced_element)
		{
			rAdjointElem.Calculate(ZERO_ADJOINT_LOAD, zero_adjoint_vector, rProcessInfo);
			noalias(rResponseGradient) = prod(rDerivativesMatrix, zero_adjoint_vector);
		}
		else
		{
			 noalias(rResponseGradient) = zero_adjoint_vector;
		}*/

       KRATOS_CATCH("")
	}

	// ==============================================================================
	void CalculateSensitivityGradient(Condition& rAdjointCondition, 
                                     const Variable<double>& rVariable,
                                     const Matrix& rDerivativesMatrix,
                                     Vector& rResponseGradient,
                                     ProcessInfo& rProcessInfo) override
	{										  
		KRATOS_TRY;

		//TODO: Rework this. Maybe not always zero vetor.
		if (rResponseGradient.size() != rDerivativesMatrix.size1())
          		rResponseGradient.resize(rDerivativesMatrix.size1(), false);
		rResponseGradient.clear();	

		KRATOS_CATCH("");
	}

	// ==============================================================================
	void CalculateSensitivityGradient(Element& rAdjointElem, 
                                      const Variable<array_1d<double,3>>& rVariable,
                                      const Matrix& rDerivativesMatrix,
                                      Vector& rResponseGradient,
                                      ProcessInfo& rProcessInfo) override
    {
      	KRATOS_TRY

		//Vector zero_adjoint_vector;	  
		//zero_adjoint_vector  = ZeroVector(rDerivativesMatrix.size1());

		// Do it only for traced element in order to prevent e.g. double disturbance if DV are the nodal coordinates
		if(rAdjointElem.Id() == m_id_of_traced_element) 
		{
			BaseType::CalculateSensitivityGradient(rAdjointElem, rVariable, rDerivativesMatrix, rResponseGradient, rProcessInfo);
			//rAdjointElem.Calculate(ZERO_ADJOINT_LOAD, zero_adjoint_vector, rProcessInfo);
			//noalias(rResponseGradient) = prod(rDerivativesMatrix, zero_adjoint_vector);
		}
		else
		{
			if (rResponseGradient.size() != rDerivativesMatrix.size1())
          		rResponseGradient.resize(rDerivativesMatrix.size1(), false);
			rResponseGradient.clear();	  
		}	

      KRATOS_CATCH("")
	}

	// ==============================================================================
	void CalculateSensitivityGradient(Condition& rAdjointCondition, 
                                      const Variable<array_1d<double,3>>& rVariable,
                                      const Matrix& rDerivativesMatrix,
                                      Vector& rResponseGradient,
                                      ProcessInfo& rProcessInfo)
    {
		KRATOS_TRY;
		BaseType::CalculateSensitivityGradient(rAdjointCondition, rVariable, rDerivativesMatrix, rResponseGradient, rProcessInfo);
        //TODO: Rework this. A zero vector is not the general case. It is valid for e.g. point loads
		/*
			if (rResponseGradient.size() != rDerivativesMatrix.size1())
          		rResponseGradient.resize(rDerivativesMatrix.size1(), false);
			rResponseGradient.Clear();*/

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
	//double m_initial_value;
	//bool m_initial_value_defined;
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

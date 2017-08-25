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

#ifndef NODAL_DISPLACEMENT_RESPONSE_FUNCTION_H
#define NODAL_DISPLACEMENT_RESPONSE_FUNCTION_H

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

class NodalDisplacementResponseFunction : public StructuralResponseFunction
{
public:
	///@name Type Definitions
	///@{

	typedef StructuralResponseFunction BaseType;
	//typedef array_1d<double, 3> array_3d;
	typedef Element::DofsVectorType DofsVectorType;
	typedef Node<3>::Pointer PointTypePointer; //try to ensure that this response function also works for 2D
	typedef VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>> VariableComponentType;
	

	

	/// Pointer definition of NodalDisplacementResponseFunction
	KRATOS_CLASS_POINTER_DEFINITION(NodalDisplacementResponseFunction);

	///@}
	///@name Life Cycle
	///@{

	/// Default constructor.
	NodalDisplacementResponseFunction(ModelPart& model_part, Parameters& responseSettings)
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
	

		// Get id of node where a displacement should be traced
		m_id_of_traced_node = responseSettings["traced_node"].GetInt();

		// Get the corresponding dof to the displacement which should be traced
		// By this response function e.g. DISPLACEMENT_X, ROTATION_X,...
		m_traced_dof_label = responseSettings["traced_dof"].GetString();


		// Get pointer to traced node
		m_traced_pNode = r_model_part.pGetNode(m_id_of_traced_node); 

		// Check if variable for traced dof is valid
		if( !(KratosComponents< VariableComponent< VectorComponentAdaptor<array_1d<double, 3> > > >::Has(m_traced_dof_label)) ) 
		{
			KRATOS_THROW_ERROR(std::invalid_argument, "Specified traced DOF is not availible. Specified DOF: ", m_traced_dof_label);
		}	
		else
		{
			const VariableComponentType& rTRACED_DOF =
            	KratosComponents<VariableComponentType>::Get(m_traced_dof_label);	
			if (m_traced_pNode->SolutionStepsDataHas(rTRACED_DOF) == false)
				KRATOS_THROW_ERROR(std::invalid_argument, "Specified DOF is not availible at traced node.", "");		
		}		

		// Check if variable for traced adjoint dof is valid
		if( !(KratosComponents< VariableComponent< VectorComponentAdaptor<array_1d<double, 3> > > >::Has(std::string("ADJOINT_") + m_traced_dof_label)) ) 
		{
			KRATOS_THROW_ERROR(std::invalid_argument, "Specified traced adjoint DOF is not availible.", "");
		}
		else
		{
			const VariableComponentType& rTRACED_ADJOINT_DOF =
            	KratosComponents<VariableComponentType>::Get(std::string("ADJOINT_") + m_traced_dof_label);		
			if (m_traced_pNode->SolutionStepsDataHas(rTRACED_ADJOINT_DOF) == false)	
				KRATOS_THROW_ERROR(std::invalid_argument, "Specified adjoint DOF is not availible at traced node.", "");		
		}
		
		m_displacement_value = 0.0;	
		this->GetNeighboringElementPointer();

	}

	/// Destructor.
	virtual ~NodalDisplacementResponseFunction()
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



		// Initialize a condition which is for this response function the adjoint load
		/*std::string label_adjoint_load_condition;

		if(m_traced_displacement_type == "DISPLACEMENT")
			label_adjoint_load_condition = POINT_LOAD;
		else if(m_traced_displacement_type == "ROTATION")
			label_adjoint_load_condition = POINT_MOMENT;
		else
			KRATOS_ERROR << "Plaese choose between DISPLACEMENT and ROTATION. Your choice is " 
			<< m_traced_displacement_type << "." << std::endl;

		if (KratosComponents<Variable<double>>::Has(label_adjoint_load_condition) == true)
        {
            const Variable<double>& r_adjoint_load_variable =
                KratosComponents<Variable<double>>::Get(label_adjoint_load_condition);
		}	
		m_traced_pNode.FastGetSolutionStepValue(r_adjoint_load_variable) = r_adjoint_load_variable.Zero();*/

		
	
		KRATOS_CATCH("");	
	}

	// ==============================================================================
	void GetNeighboringElementPointer() 
	{
		KRATOS_TRY;	

		ModelPart& r_model_part = this->GetModelPart();
		bool neighboring_element_found = false;

		for (ModelPart::ElementIterator it = r_model_part.ElementsBegin(); it != r_model_part.ElementsEnd(); ++it)
        {
			const unsigned int NumberOfNodes = it->GetGeometry().size();
            for(unsigned int i = 0; i < NumberOfNodes; ++i)
			{
				int current_node_id = it->GetGeometry()[i].Id();
				if(current_node_id == m_id_of_traced_node)
				{
					m_neighboring_pElement = r_model_part.pGetElement(it->Id());
					neighboring_element_found = true;
					break;
				}
			}
			if(neighboring_element_found) { break; }
		}	
		if(!neighboring_element_found)
			KRATOS_ERROR << "No neighboring element is availible for the traced node." << std::endl;

		KRATOS_CATCH("");
		
	}

	// ==============================================================================
	double CalculateValue(ModelPart& rModelPart) override  //TODO: test this function
	{
		KRATOS_TRY;

		const VariableComponentType& rTRACED_DOF =
            KratosComponents<VariableComponentType>::Get(m_traced_dof_label);	

		m_displacement_value = m_traced_pNode->FastGetSolutionStepValue(rTRACED_DOF);
		
		return m_displacement_value;

		//std::cout << ("current displacement value = ") << m_displacement_value << std::endl;

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

		return m_displacement_value;

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
		return "NodalDisplacementResponseFunction";
	}

	/// Print information about this object.
	virtual void PrintInfo(std::ostream &rOStream) const    old function, don't needed for sensitivity analysis
	{
		rOStream << "NodalDisplacementResponseFunction";
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
		KRATOS_TRY;
		       
		if (rResponseGradient.size() != rAdjointMatrix.size1())
            rResponseGradient.resize(rAdjointMatrix.size1(), false);

        rResponseGradient.clear();

		if(rAdjointElem.Id() == m_neighboring_pElement->Id() )
		{
			DofsVectorType dofs_of_lement; 
			m_neighboring_pElement->GetDofList(dofs_of_lement,rProcessInfo);
			const VariableComponentType& rTRACED_ADJOINT_DOF =
            	KratosComponents<VariableComponentType>::Get(std::string("ADJOINT_") + m_traced_dof_label);	
			int index = 0;
			for(unsigned int i = 0; i < dofs_of_lement.size(); ++i)
			{
				if(m_traced_pNode->pGetDof(rTRACED_ADJOINT_DOF) == dofs_of_lement[i])
				{
					rResponseGradient[i] = -1; //TODO: Check for correct sign!
				}
				index++;
			}
		}

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
	void CalculateSensitivityGradient(Element& rAdjointElem, 
                                      const Variable<double>& rVariable,
                                      const Matrix& rDerivativesMatrix,
                                      Vector& rResponseGradient,
                                      ProcessInfo& rProcessInfo) override
    {
      	KRATOS_TRY

		if (rResponseGradient.size() != rDerivativesMatrix.size1())
        	rResponseGradient.resize(rDerivativesMatrix.size1(), false);

        rResponseGradient.clear();


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

		if (rResponseGradient.size() != rDerivativesMatrix.size1())
        	rResponseGradient.resize(rDerivativesMatrix.size1(), false);

        rResponseGradient.clear();

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

        if (rResponseGradient.size() != rDerivativesMatrix.size1())
        	rResponseGradient.resize(rDerivativesMatrix.size1(), false);

        rResponseGradient.clear();

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
	double m_displacement_value; 
	double mDelta;
	int m_id_of_traced_node;
	std::string m_traced_dof_label;
	PointTypePointer  m_traced_pNode;
	Element::Pointer m_neighboring_pElement;
	
	
	//std::string m_traced_displacement_type;
	//std::string m_traced_direction;

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
	//      NodalDisplacementResponseFunction& operator=(SNodalDisplacementResponseFunction const& rOther);

	/// Copy constructor.
	//      NodalDisplacementResponseFunction(NodalDisplacementResponseFunction const& rOther);

	///@}

}; // Class NodalDisplacementResponseFunction

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

} // namespace Kratos.

#endif // NODAL_DISPLACEMENT_RESPONSE_FUNCTION_H

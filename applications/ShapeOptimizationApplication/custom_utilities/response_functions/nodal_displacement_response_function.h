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
	//typedef VariableComponent< VectorComponentAdaptor<array_1d<double, 3> > > component_type;	
	

	

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
	

		//get id of node where a diplacement should be traced
		m_id_of_traced_node = responseSettings["traced_node"].GetInt();

		//get the corresponding dof to the displacement which should be traced
		//by this response function e.g. DISPLACEMENT_X or ROTATION_X
		m_traced_dof_label = responseSettings["traced_dof"].GetString();
		m_traced_dof_label = "ADJOINT_" + m_traced_dof_label;

		// get variable for traced DOF
		/*if( KratosComponents< VariableComponent< VectorComponentAdaptor<array_1d<double, 3> > > >::Has(m_traced_dof_label) ) 
        	m_var_component = KratosComponents< component_type >::Get(m_traced_dof_label);
		else
			KRATOS_THROW_ERROR(std::invalid_argument, "Specified DOF is not availible. Specified DOF: ", m_traced_dof_label);	*/

		m_traced_pNode = r_model_part.pGetNode(m_id_of_traced_node); 
	
		m_displacement_value = 0.0;	
		m_adjoint_load_is_set = false;

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
	double CalculateValue(ModelPart& rModelPart) override
	{
		KRATOS_TRY;

		// Working variables
		//ProcessInfo &CurrentProcessInfo = rModelPart.GetProcessInfo();

		//TODO: test this
		m_displacement_value = m_traced_pNode->FastGetSolutionStepValue(
								KratosComponents< VariableComponent< VectorComponentAdaptor<array_1d<double, 3> > > >::Get(m_traced_dof_label));
		

		return m_displacement_value;

		std::cout << ("current displacement value = ") << m_displacement_value << std::endl;

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
	void CalculateGradient(const Element& rAdjointElem, const Matrix& rAdjointMatrix,
                                   Vector& rResponseGradient,
                                   ProcessInfo& rProcessInfo) override
	{
		KRATOS_TRY;
		       
		if (rResponseGradient.size() != rAdjointMatrix.size1())
            rResponseGradient.resize(rAdjointMatrix.size1(), false);

        rResponseGradient.clear();

		//TODO: find a more elegant solution for this
		if(!m_adjoint_load_is_set)
		{
			const unsigned int NumberOfNodes = rAdjointElem.GetGeometry().size();
            for(unsigned int i = 0; i < NumberOfNodes; ++i)
			{
				int current_node_id = rAdjointElem.GetGeometry()[i].Id();
				if(current_node_id == m_id_of_traced_node)
				{
					DofsVectorType dofs_of_lement; 
					ModelPart& r_model_part = this->GetModelPart();
					Element::Pointer traced_pElement = r_model_part.pGetElement(rAdjointElem.Id());
					traced_pElement->GetDofList(dofs_of_lement,rProcessInfo);
					int index = 0;
					for(unsigned int i = 0; i < dofs_of_lement.size(); ++i)
					{
						if(m_traced_pNode->pGetDof(KratosComponents< VariableComponent< VectorComponentAdaptor<array_1d<double, 3> > > >::Get(m_traced_dof_label)) == dofs_of_lement[i])
						{
							rResponseGradient[i] = -1; //TODO: Check for correct sign!
							m_adjoint_load_is_set = true;
						}
						index++;
					}// loop dofs of element
				}// enter if one of the element nodes is the traced node
			}// loop nodes of element  	
		}// enter only if the adjoint load was not given by a neighboring element

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
	bool m_adjoint_load_is_set;
	PointTypePointer  m_traced_pNode;
	
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

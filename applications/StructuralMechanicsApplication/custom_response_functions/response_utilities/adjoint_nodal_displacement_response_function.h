// ==============================================================================
//  KratosStructuralMechanicsApplication
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Fusseder Martin
//                   martin.fusseder@tum.de
//
// ==============================================================================

#ifndef ADJOINT_NODAL_DISPLACEMENT_RESPONSE_FUNCTION_H
#define ADJOINT_NODAL_DISPLACEMENT_RESPONSE_FUNCTION_H

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
#include "adjoint_structural_response_function.h"


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

class AdjointNodalDisplacementResponseFunction : public AdjointStructuralResponseFunction
{
public:
	///@name Type Definitions
	///@{

	typedef AdjointStructuralResponseFunction BaseType;
	//typedef array_1d<double, 3> array_3d;
	typedef Element::DofsVectorType DofsVectorType;
	typedef Node<3>::Pointer PointTypePointer; //try to ensure that this response function also works for 2D
	typedef VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>> VariableComponentType;




	/// Pointer definition of AdjointNodalDisplacementResponseFunction
	KRATOS_CLASS_POINTER_DEFINITION(AdjointNodalDisplacementResponseFunction);

	///@}
	///@name Life Cycle
	///@{

	/// Default constructor.
	AdjointNodalDisplacementResponseFunction(ModelPart& model_part, Parameters& responseSettings)
	: AdjointStructuralResponseFunction(model_part, responseSettings)
	{
		ModelPart& r_model_part = this->GetModelPart();

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
	virtual ~AdjointNodalDisplacementResponseFunction()
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

		m_displacement_value = rModelPart.Nodes()[(m_traced_pNode->Id())].FastGetSolutionStepValue(rTRACED_DOF, 0);

		return m_displacement_value;

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
		return "AdjointNodalDisplacementResponseFunction";
	}

	/// Print information about this object.
	virtual void PrintInfo(std::ostream &rOStream) const    old function, don't needed for sensitivity analysis
	{
		rOStream << "AdjointNodalDisplacementResponseFunction";
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

		if( rAdjointElem.Id() == m_neighboring_pElement->Id() )
		{
			DofsVectorType dofs_of_lement;
			m_neighboring_pElement->GetDofList(dofs_of_lement,rProcessInfo);

			const VariableComponentType& rTRACED_ADJOINT_DOF =
            	KratosComponents<VariableComponentType>::Get(std::string("ADJOINT_") + m_traced_dof_label);

			for(unsigned int i = 0; i < dofs_of_lement.size(); ++i)
			{
				if(m_traced_pNode->pGetDof(rTRACED_ADJOINT_DOF) == dofs_of_lement[i])
				{
					rResponseGradient[i] = -1;
				}
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
	double m_displacement_value;
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
	//      AdjointNodalDisplacementResponseFunction& operator=(SAdjointNodalDisplacementResponseFunction const& rOther);

	/// Copy constructor.
	//      AdjointNodalDisplacementResponseFunction(AdjointNodalDisplacementResponseFunction const& rOther);

	///@}

}; // Class AdjointNodalDisplacementResponseFunction

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

} // namespace Kratos.

#endif // ADJOINT_NODAL_DISPLACEMENT_RESPONSE_FUNCTION_H

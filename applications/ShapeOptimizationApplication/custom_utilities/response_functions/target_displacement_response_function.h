// ==============================================================================
//  KratosShapeOptimizationApplication
//
//  License:         BSD License
//                   license: ShapeOptimizationApplication/license.txt
//
//  Main authors:    Kayis AliNuri   
//                   alinuri.kayis@tum.de
//
// ==============================================================================

#ifndef TARGET_DISPLACEMENT_RESPONSE_FUNCTION_H
#define TARGET_DISPLACEMENT_RESPONSE_FUNCTION_H

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
#include "../../kratos/includes/condition.h"
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

class TargetDisplacementResponseFunction : public StructuralResponseFunction
{
public:
	///@name Type Definitions
	///@{

	typedef StructuralResponseFunction BaseType;
	typedef std::size_t IndexType;
	//typedef Condition ConditionType;

	/// Pointer definition of TargetDisplacementResponseFunction
	KRATOS_CLASS_POINTER_DEFINITION(TargetDisplacementResponseFunction);

	///@}
	///@name Life Cycle
	///@{

	/// Default constructor.
	TargetDisplacementResponseFunction(ModelPart& model_part, Parameters& responseSettings)
	: StructuralResponseFunction(model_part, responseSettings)
	{
		// Get name of condition to check while assigning stock sensitivity values to RHS 
		m_condition_name = 
			responseSettings["condition_name"].GetString();
		// Get SubModelPart name for the nodes where new conditions are created
		mStockSensitivityModelPartName = 
			responseSettings["stock_sentivity_model_part_name"].GetString();
		
		Parameters nodal_stock_sensitivity_variables = 
			responseSettings["nodal_stock_sensitivity_variables"];
        mNodalStockSensitivityVariables.resize(nodal_stock_sensitivity_variables.size());
        for (unsigned int i = 0; i < nodal_stock_sensitivity_variables.size(); ++i)
			mNodalStockSensitivityVariables[i] = nodal_stock_sensitivity_variables[i].GetString();
			
		

	}

	/// Destructor.
	virtual ~TargetDisplacementResponseFunction()
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
		ModelPart& r_model_part = this->GetModelPart();

		if (r_model_part.HasSubModelPart(mStockSensitivityModelPartName) == false)
		KRATOS_ERROR << "No sub model part \"" << mStockSensitivityModelPartName
					 << "\"" << std::endl;

		ModelPart& r_sub_model_part = r_model_part.GetSubModelPart(mStockSensitivityModelPartName);

		//getting the array of the conditions
		const int nconditions = static_cast<int>(r_model_part.Conditions().size());

		// creating new conditions and adding them into submodelpart
		ModelPart::NodeIterator nodes_begin = r_sub_model_part.Nodes().begin();
		ModelPart::NodeIterator nodes_end = r_sub_model_part.Nodes().end();
		unsigned int index = 1;
		for (auto it = nodes_begin; it != nodes_end; ++it)
		{
			std::vector<IndexType> nodeIds;
			nodeIds.push_back(it->Id());
			std::vector< ModelPart::IndexType > ConditionsIds;
			ConditionsIds.push_back(nconditions + index);
			unsigned int elem_id = nconditions + index;
			unsigned int prop_id = 0;
			r_model_part.CreateNewCondition(m_condition_name,elem_id,nodeIds,r_model_part.pGetProperties(prop_id));
			r_sub_model_part.AddConditions(ConditionsIds);
			index++;
		}

		/*
	
		BaseType::Initialize();
		ModelPart& r_model_part = this->GetModelPart();
		//getting the array of the conditions
		const int nconditions = static_cast<int>(r_model_part.Conditions().size());
	#pragma omp parallel
    {
		ModelPart::NodeIterator nodes_begin;
		ModelPart::NodeIterator nodes_end;
		OpenMPUtils::PartitionedIterators(
                r_model_part.GetSubModelPart(mStockSensitivityModelPartName).Nodes(),
                nodes_begin, nodes_end);
		unsigned int index = 1;
		for (auto it = nodes_begin; it != nodes_end; ++it)
		{
			std::vector<IndexType> nodeIds;
			nodeIds.push_back(it->Id());
			std::vector< ModelPart::IndexType > ConditionsIds;
			ConditionsIds.push_back(nconditions + index);
			r_model_part.CreateNewCondition(m_condition_name,ConditionsIds,nodeIds,mp.pProperties()[0]);
			r_model_part.GetSubModelPart(mStockSensitivityModelPartName).AddConditions(ConditionIds);
			index++;
		}
		
            ModelPart::ConditionIterator conditions_begin;
            ModelPart::ConditionIterator conditions_end;
            OpenMPUtils::PartitionedIterators(
                r_model_part.GetSubModelPart(mSensitivityModelPartName).Conditions(),
                conditions_begin, conditions_end);
            for (auto it = conditions_begin; it != conditions_end; ++it)
                it->SetValue(UPDATE_SENSITIVITIES, true);   
                
    }  

		*/
	
	
		KRATOS_CATCH("");	
	}

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
		return "TargetDisplacementResponseFunction";
	}

	/// Print information about this object.
	virtual void PrintInfo(std::ostream &rOStream) const    old function, don't needed for sensitivity analysis
	{
		rOStream << "TargetDisplacementResponseFunction";
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
	void CalculateGradient(const Condition& rAdjointCondition, const Matrix& rAdjointMatrix,
                                   Vector& rResponseGradient,
                                   ProcessInfo& rProcessInfo) override
	{
		KRATOS_TRY;

		if (rResponseGradient.size() != rAdjointMatrix.size1())
			rResponseGradient.resize(rAdjointMatrix.size1(), false);

		rResponseGradient.clear();

		for (auto label : mNodalStockSensitivityVariables)
        {
				     
            if (KratosComponents<Variable<array_1d<double, 3>>>::Has(label) == true)
            {
                const Variable<array_1d<double, 3>>& r_variable =
                    KratosComponents<Variable<array_1d<double, 3>>>::Get(label);
				// Checking can be done by variable in Condition::Calculate()	
				if (rAdjointCondition.Info().compare(m_condition_name) == 0)
				{
					const unsigned int dim = rAdjointCondition.GetGeometry().WorkingSpaceDimension();
					const unsigned int NumberOfNodes = rAdjointCondition.GetGeometry().size();
					for(unsigned int i = 0; i < NumberOfNodes; ++i)
					{
						for (unsigned int ii = 0; ii < dim; ++ii)
						{
							rResponseGradient[ii] =
								-1 * rAdjointCondition.GetGeometry()[i].FastGetSolutionStepValue(r_variable, 0)[ii];
						}
					}
				}

            }
            else
                KRATOS_ERROR << "Unsupported variable: " << label << "." << std::endl;
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
	std::string m_condition_name;
	std::string mStockSensitivityModelPartName;
	std::vector<std::string> mNodalStockSensitivityVariables;
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
	//      TargetDisplacementResponseFunction& operator=(STargetDisplacementResponseFunction const& rOther);

	/// Copy constructor.
	//      TargetDisplacementResponseFunction(TargetDisplacementResponseFunction const& rOther);

	///@}

}; // Class TargetDisplacementResponseFunction

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

} // namespace Kratos.

#endif // TARGET_DISPLACEMENT_RESPONSE_FUNCTION_H



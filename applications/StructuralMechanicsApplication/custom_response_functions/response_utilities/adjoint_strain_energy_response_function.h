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

#ifndef ADJOINT_STRAIN_ENERGY_RESPONSE_FUNCTION_H
#define ADJOINT_STRAIN_ENERGY_RESPONSE_FUNCTION_H

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

class AdjointStrainEnergyResponseFunction : public AdjointStructuralResponseFunction
{
public:
	///@name Type Definitions
	///@{

	typedef AdjointStructuralResponseFunction BaseType;
	typedef array_1d<double, 3> array_3d;



	/// Pointer definition of AdjointStrainEnergyResponseFunction
	KRATOS_CLASS_POINTER_DEFINITION(AdjointStrainEnergyResponseFunction);

	///@}
	///@name Life Cycle
	///@{

	/// Default constructor.
	AdjointStrainEnergyResponseFunction(ModelPart& model_part, Parameters& responseSettings)
	: AdjointStructuralResponseFunction(model_part, responseSettings)
	{

		// Initialize member variables to NULL
		m_current_response_value = 0.0;

	}

	/// Destructor.
	virtual ~AdjointStrainEnergyResponseFunction()
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

		// It is necessary to initialize the elements since no adjoint problem is solved for this response type.
		// For this response type the elements are only created!
		ModelPart& r_model_part = this->GetModelPart();
	#pragma omp parallel
        {
            ModelPart::ElementIterator elements_begin;
            ModelPart::ElementIterator elements_end;
            OpenMPUtils::PartitionedIterators(r_model_part.Elements(), elements_begin, elements_end);
            for (auto it = elements_begin; it != elements_end; ++it)
                it->Initialize();
        }
		// TODO: Check if initialization is also necessary for conditions!

		KRATOS_CATCH("");
	}

	// ==============================================================================
	double CalculateValue(ModelPart& rModelPart) override
	{
		KRATOS_TRY;

		ModelPart& r_model_part = rModelPart; 
		ProcessInfo &CurrentProcessInfo = r_model_part.GetProcessInfo();
		m_current_response_value = 0.0;

		// Check if there are at the time of calling adjoint or primal elements
		// TODO: ist there a smarter solution to get this information??
		std::string element_name = r_model_part.Elements()[r_model_part.ElementsBegin()->Id()].Info();
		if( element_name.find("Adjoint", 0) != std::string::npos)
			KRATOS_ERROR << "Calculate value for strain energy response is not availibe when using adjoint elements" << std::endl;
			
		// Sum all elemental strain energy values calculated as: W_e = u_e^T K_e u_e
		for (auto& elem_i : r_model_part.Elements())
		{
			Matrix LHS;
			Vector RHS;
			Vector disp;

			// Get state solution relevant for energy calculation
			elem_i.GetValuesVector(disp,0);

			elem_i.CalculateLocalSystem(LHS,RHS,CurrentProcessInfo);

			// Compute strain energy
			m_current_response_value += 0.5 * inner_prod(disp,prod(LHS,disp));
 		}

		return m_current_response_value;

		KRATOS_CATCH("");
	}
	// --------------------------------------------------------------------------
	/*double GetInitialValue()
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

		return m_current_response_value;

		KRATOS_CATCH("");
	}*/

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
		return "AdjointStrainEnergyResponseFunction";
	}

	/// Print information about this object.
	virtual void PrintInfo(std::ostream &rOStream) const
	{
		rOStream << "AdjointStrainEnergyResponseFunction";
	}

	/// Print object's data.
	virtual void PrintData(std::ostream &rOStream) const
	{
	}

	///@}
	///@name Friends
	///@{

	///@}



	// =============================================================================
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

		// There will be a mistake, if body forces are considered. Because the elements are responsible for the body forces!

     	 KRATOS_CATCH("")
	}

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

		// There will be a mistake, if body forces are considered. Because the elements are responsible for the body forces!

        KRATOS_CATCH("")
	}

	// ==============================================================================
	void CalculateSensitivityGradient(Condition& rAdjointCondition,
                                              const Variable<array_1d<double,3>>& rVariable,
                                              const Matrix& rDerivativesMatrix,
                                              Vector& rResponseGradient,
                                              ProcessInfo& rProcessInfo) override
    {
		KRATOS_TRY;

		Vector adjoint_variables;

		rAdjointCondition.GetValuesVector(adjoint_variables); // = 0.5*u

		if (adjoint_variables.size() != rDerivativesMatrix.size2())
			KRATOS_ERROR << "Size of adjoint vector does not fit to the size of the pseudo load!" << std::endl;

		if (rResponseGradient.size() != rDerivativesMatrix.size2())
			rResponseGradient.resize(adjoint_variables.size(), false);

		noalias(rResponseGradient) = prod(rDerivativesMatrix, adjoint_variables);

		KRATOS_CATCH("");
	}

	// ==============================================================================
	void CalculateSensitivityGradient(Condition& rAdjointCondition,
                                              const Variable<double>& rVariable,
                                              const Matrix& rDerivativesMatrix,
                                              Vector& rResponseGradient,
                                              ProcessInfo& rProcessInfo) override
    {
		KRATOS_TRY;

		Vector adjoint_variables;

		rAdjointCondition.GetValuesVector(adjoint_variables); // = 0.5*u

		if (adjoint_variables.size() != rDerivativesMatrix.size2())
			KRATOS_ERROR << "Size of adjoint vector does not fit to the size of the pseudo load!" << std::endl;

		if (rResponseGradient.size() != rDerivativesMatrix.size2())
			rResponseGradient.resize(adjoint_variables.size(), false);

		noalias(rResponseGradient) = prod(rDerivativesMatrix, adjoint_variables);

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

	double m_current_response_value;

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
	//      AdjointStrainEnergyResponseFunction& operator=(AdjointStrainEnergyResponseFunction const& rOther);

	/// Copy constructor.
	//      AdjointStrainEnergyResponseFunction(AdjointStrainEnergyResponseFunction const& rOther);

	///@}

}; // Class AdjointStrainEnergyResponseFunction

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

} // namespace Kratos.

#endif // ADJOINT_STRAIN_ENERGY_RESPONSE_FUNCTION_H

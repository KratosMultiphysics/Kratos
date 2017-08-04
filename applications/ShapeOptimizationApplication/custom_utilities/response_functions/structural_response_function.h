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

#ifndef STRUCTURAL_RESPONSE_FUNCTION_H
#define STRUCTURAL_RESPONSE_FUNCTION_H

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


#include "includes/kratos_parameters.h"
#include "includes/ublas_interface.h"
#include "utilities/openmp_utils.h"

#include "shape_optimization_application.h"

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

class StructuralResponseFunction
{
public:
	///@name Type Definitions
	///@{


	
	typedef array_1d<double, 3> array_3d;

	

	/// Pointer definition of StructuralResponseFunction
	KRATOS_CLASS_POINTER_DEFINITION(StructuralResponseFunction);

	///@}
	///@name Life Cycle
	///@{

	/// Default constructor.
	StructuralResponseFunction(ModelPart& model_part, Parameters& responseSettings)
	: mr_model_part(model_part)
	{
		KRATOS_TRY;
	
		//mSensitivityModelPartName = responseSettings["sensitivity_model_part_name"].GetString();

        Parameters nodal_sensitivity_variables = responseSettings["nodal_sensitivity_variables"];
        mNodalSensitivityVariables.resize(nodal_sensitivity_variables.size());
        for (unsigned int i = 0; i < nodal_sensitivity_variables.size(); ++i)
			mNodalSensitivityVariables[i] = nodal_sensitivity_variables[i].GetString();

		Parameters element_sensitivity_variables = responseSettings["element_sensitivity_variables"];
        mElementSensitivityVariables.resize(element_sensitivity_variables.size());
        for (unsigned int i = 0; i < element_sensitivity_variables.size(); ++i)
			mElementSensitivityVariables[i] = element_sensitivity_variables[i].GetString();	

		Parameters condition_sensitivity_variables = responseSettings["condition_sensitivity_variables"];
        mConditionSensitivityVariables.resize(condition_sensitivity_variables.size());
        for (unsigned int i = 0; i < condition_sensitivity_variables.size(); ++i)
			mConditionSensitivityVariables[i] = condition_sensitivity_variables[i].GetString();	

		KRATOS_CATCH("");
	}

	/// Destructor.
	virtual ~StructuralResponseFunction()
	{
	}

	///@}
	///@name Operators
	///@{

	///@}
	///@name Operations
	///@{

	ModelPart& GetModelPart()
    {
      return mr_model_part;
    }

    const ModelPart& GetModelPart() const
    {
      return mr_model_part;
	}

	// ==============================================================================
	virtual void Initialize()
	{
		KRATOS_TRY;

		std::cout << ("I was in base class to initialize") << std::endl;
	//move this to base class---------------------------------------------------------------
        // Set sensitivity variables to zero.
        for (auto label : mNodalSensitivityVariables)
        {
            if (KratosComponents<Variable<double>>::Has(label) == true)
            {
                const Variable<double>& r_variable =
                    KratosComponents<Variable<double>>::Get(label);

#pragma omp parallel
                {
                    ModelPart::NodeIterator nodes_begin;
                    ModelPart::NodeIterator nodes_end;
                    OpenMPUtils::PartitionedIterators(mr_model_part.Nodes(),
                                                      nodes_begin, nodes_end);
                    for (auto it = nodes_begin; it != nodes_end; ++it)
                        it->FastGetSolutionStepValue(r_variable) = r_variable.Zero();
                }
            }
            else if (KratosComponents<Variable<array_1d<double, 3>>>::Has(label) == true)
            {
                const Variable<array_1d<double, 3>>& r_variable =
                    KratosComponents<Variable<array_1d<double, 3>>>::Get(label);

#pragma omp parallel
                {
                    ModelPart::NodeIterator nodes_begin;
                    ModelPart::NodeIterator nodes_end;
                    OpenMPUtils::PartitionedIterators(mr_model_part.Nodes(),
                                                      nodes_begin, nodes_end);
                    for (auto it = nodes_begin; it != nodes_end; ++it)
                        it->FastGetSolutionStepValue(r_variable) = r_variable.Zero();
                }
            }
            else
                KRATOS_ERROR << "Unsupported variable: " << label << "." << std::endl;

	}// end loop over mNodalSensitivityVariables



		KRATOS_CATCH("");
	}

	// --------------------------------------------------------------------------
	virtual void CalculateValue()  = 0;
	
	// --------------------------------------------------------------------------
	virtual double GetInitialValue() = 0;

	// --------------------------------------------------------------------------
	virtual double GetValue() = 0;

	// --------------------------------------------------------------------------
	//virtual boost::python::dict get_gradient() = 0;

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
		return "StructuralResponseFunction";
	}

	/// Print information about this object.
	virtual void PrintInfo(std::ostream &rOStream) const
	{
		rOStream << "StructuralResponseFunction";
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
	virtual void CalculateGradient(const Element& rElem, const Matrix& rLHS, Vector& rOutput, 
									ProcessInfo& rProcessInfo)
	{
		KRATOS_TRY;

        if (rOutput.size() != rLHS.size1())
            rOutput.resize(rLHS.size1(), false);

        for (unsigned int k = 0; k < rOutput.size(); ++k)
            rOutput[k] = 0.0;

		KRATOS_CATCH("");

	}

	// ==============================================================================
	virtual void UpdateSensitivities()
	{
        KRATOS_TRY;

		// Nodal sensitivity variables e.g. node coordinates---------------------------------
        for (auto label : mNodalSensitivityVariables)
        {
            if (KratosComponents<Variable<double>>::Has(label) == true)
            {
                const Variable<double>& r_variable =
                    KratosComponents<Variable<double>>::Get(label);
                this->UpdateNodalSensitivities(r_variable);
            }
            else if (KratosComponents<Variable<array_1d<double,3>>>::Has(label) == true)
            {
                const Variable<array_1d<double,3>>& r_variable =
                    KratosComponents<Variable<array_1d<double,3>>>::Get(label);
                this->UpdateNodalSensitivities(r_variable);
            }
            else
                KRATOS_ERROR << "Unsupported node variable: " << label << "." << std::endl;
        }

		// Elemental sensitivity variables e.g. 2nd moment of inertia---------------------------------
		for (auto label : mElementSensitivityVariables)
        {
            if (KratosComponents<Variable<double>>::Has(label) == true)
            {
                const Variable<double>& r_variable =
                    KratosComponents<Variable<double>>::Get(label);
                this->UpdateElementSensitivities(r_variable);
            }
            else if (KratosComponents<Variable<array_1d<double,3>>>::Has(label) == true)
            {
                const Variable<array_1d<double,3>>& r_variable =
                    KratosComponents<Variable<array_1d<double,3>>>::Get(label);
                this->UpdateElementSensitivities(r_variable);
            }
            else
                KRATOS_ERROR << "Unsupported element variable: " << label << "." << std::endl;
        }

		std::cout << ("I looped all design variables and all elements") << std::endl;

		// Add condition sensitivity variables e.g. point load---------------------------------------
		for (auto label : mConditionSensitivityVariables)
        {
            if (KratosComponents<Variable<double>>::Has(label) == true)
            {
                const Variable<double>& r_variable =
                    KratosComponents<Variable<double>>::Get(label);
                this->UpdateConditionSensitivities(r_variable);
            }
            else if (KratosComponents<Variable<array_1d<double,3>>>::Has(label) == true)
            {
                const Variable<array_1d<double,3>>& r_variable =
                    KratosComponents<Variable<array_1d<double,3>>>::Get(label);
                this->UpdateConditionSensitivities(r_variable);
            }
            else
                KRATOS_ERROR << "Unsupported condition variable: " << label << "." << std::endl;
        }
      
        KRATOS_CATCH("");
	}


	protected:
	///@name Protected static Member Variables
	///@{

	///@}
	///@name Protected member Variables
	///@{

	std::vector<std::string> mNodalSensitivityVariables;
	std::vector<std::string> mElementSensitivityVariables;
	std::vector<std::string> mConditionSensitivityVariables;	

	///@}
	///@name Protected Operators
	///@{

	///@}
	///@name Protected Operations
	///@{


	// ==============================================================================
	template <typename TDataType>
    void UpdateNodalSensitivities(Variable<TDataType> const& rSensitivityVariable) 
	{
		KRATOS_TRY;

		ModelPart& r_model_part = this->GetModelPart();
        ProcessInfo& r_process_info = r_model_part.GetProcessInfo();
        Vector sensitivity_vector;
    	Vector response_gradient;
        Vector adjoint_vector;
        Matrix sensitivity_matrix;

		// Loop elements
		for (ModelPart::ElementIterator it = r_model_part.ElementsBegin(); it != r_model_part.ElementsEnd(); ++it)
        {

			Element::GeometryType& r_geom = it->GetGeometry();	

            // Compute the pseudo load
            it->CalculateSensitivityMatrix(
                rSensitivityVariable, sensitivity_matrix, r_process_info);

            // This part of the sensitivity is computed from the objective
            // with primal variables treated as constant.
            this->CalculateSensitivityGradient(
                *it, rSensitivityVariable, sensitivity_matrix,
                response_gradient, r_process_info);

			// Get the adjoint displacement field
			this->GetAdjointVariables(*it, adjoint_vector);	
            //it->GetValuesVector(adjoint_vector);

            if (sensitivity_vector.size() != sensitivity_matrix.size1())
                sensitivity_vector.resize(sensitivity_matrix.size1(), false);

			// Compute the whole sensitivity
            noalias(sensitivity_vector) =
                    			(prod(sensitivity_matrix, adjoint_vector) +
                                    response_gradient);

            this->AssembleNodalSensitivityContribution(
                rSensitivityVariable, sensitivity_vector, r_geom);  //----> check for correct output
        }

		// Reset 
		/*sensitivity_vector =  Vector(0);
		response_gradient = Vector(0);
        adjoint_vector = Vector(0);
        sensitivity_matrix = Matrix (0,0);

		// Loop conditions
		const int nconditions = static_cast<int>(mr_model_part.Conditions().size());
		ModelPart::ConditionsContainerType::iterator cond_begin = mr_model_part.ConditionsBegin();
		for (int k = 0; k < nconditions; k++)
		{
			ModelPart::ConditionsContainerType::iterator cond_i = cond_begin + k;

			Condition::GeometryType& r_geom = cond_i->GetGeometry();	

            // Compute the pseudo load
            cond_i->CalculateSensitivityMatrix(
                rSensitivityVariable, sensitivity_matrix, r_process_info);

            // This part of the sensitivity is computed from the objective
            // with primal variables treated as constant.
            this->CalculateSensitivityGradient(
                *cond_i, rSensitivityVariable, sensitivity_matrix,
                response_gradient, r_process_info);

			// Get the adjoint displacement field
			this->GetAdjointVariables(*cond_i, adjoint_vector);	
            //it->GetValuesVector(adjoint_vector);

            if (sensitivity_vector.size() != sensitivity_matrix.size1())
                sensitivity_vector.resize(sensitivity_matrix.size1(), false);

			// Compute the whole sensitivity
            noalias(sensitivity_vector) =
                    			(prod(sensitivity_matrix, adjoint_vector) +
                                    response_gradient);

            this->AssembleNodalSensitivityContribution(
                rSensitivityVariable, sensitivity_vector, r_geom);	//----> check for correct output
        }*/
    
        r_model_part.GetCommunicator().AssembleCurrentData(rSensitivityVariable);

		KRATOS_CATCH("");
	}

	// ==============================================================================
	template <typename TDataType>
	void UpdateElementSensitivities(Variable<TDataType> const& rSensitivityVariable)  
	{
		KRATOS_TRY;

		ModelPart& r_model_part = this->GetModelPart();
        ProcessInfo& r_process_info = r_model_part.GetProcessInfo();
        Vector sensitivity_vector;
    	Vector response_gradient;
        Vector adjoint_vector;
        Matrix sensitivity_matrix;

		for (ModelPart::ElementIterator it = r_model_part.ElementsBegin(); it != r_model_part.ElementsEnd(); ++it)
        {

            // Compute the pseudo load
            it->CalculateSensitivityMatrix(
            	rSensitivityVariable, sensitivity_matrix, r_process_info);

            // This part of the sensitivity is computed from the objective
            // with primal variables treated as constant.
            this->CalculateSensitivityGradient(
                    *it, rSensitivityVariable, sensitivity_matrix,
                    response_gradient, r_process_info);
	
			// Get the adjoint displacement field
			this->GetAdjointVariables(*it, adjoint_vector);			
            //it->GetValuesVector(adjoint_vector);

            if (sensitivity_vector.size() != sensitivity_matrix.size1())
               sensitivity_vector.resize(sensitivity_matrix.size1(), false);


			// Compute the whole sensitivity
            noalias(sensitivity_vector) =
                				(prod(sensitivity_matrix, adjoint_vector) +
                                 response_gradient);

			std::cout << ("element sensitivty = ") << sensitivity_vector[0] << std::endl;
               // this->AssembleNodalSensitivityContribution(
                  //  rSensitivityVariable, sensitivity_vector[k], r_geom);		//----> check for correct output
        }
       

        r_model_part.GetCommunicator().AssembleCurrentData(rSensitivityVariable);

		KRATOS_CATCH("");
	}

	// ==============================================================================
	template <typename TDataType>
    void UpdateConditionSensitivities(Variable<TDataType> const& rSensitivityVariable) 
	{
		KRATOS_TRY;

		ModelPart& r_model_part = this->GetModelPart();
		ProcessInfo& r_process_info = r_model_part.GetProcessInfo();
		Vector sensitivity_vector;
    	Vector response_gradient;
        Vector adjoint_vector;
        Matrix sensitivity_matrix;

		const int nconditions = static_cast<int>(mr_model_part.Conditions().size());
		ModelPart::ConditionsContainerType::iterator cond_begin = mr_model_part.ConditionsBegin();
		for (int k = 0; k < nconditions; k++)
		{
			ModelPart::ConditionsContainerType::iterator cond_i = cond_begin + k;

			Condition::GeometryType& r_geom = cond_i->GetGeometry();	

            // Compute the pseudo load
            cond_i->CalculateSensitivityMatrix(
                rSensitivityVariable, sensitivity_matrix, r_process_info);

            // This part of the sensitivity is computed from the objective
            // with primal variables treated as constant.
            this->CalculateSensitivityGradient(
                *cond_i, rSensitivityVariable, sensitivity_matrix,
                response_gradient, r_process_info);

			// Get the adjoint displacement field
			this->GetAdjointVariables(*cond_i, adjoint_vector);	
            //it->GetValuesVector(adjoint_vector);

            if (sensitivity_vector.size() != sensitivity_matrix.size1())
                sensitivity_vector.resize(sensitivity_matrix.size1(), false);

			// Compute the whole sensitivity
            noalias(sensitivity_vector) =
                    			(prod(sensitivity_matrix, adjoint_vector) +
                                    response_gradient);

            this->AssembleNodalSensitivityContribution(
                rSensitivityVariable, sensitivity_vector, r_geom); //----> check for correct output
        }
    
        r_model_part.GetCommunicator().AssembleCurrentData(rSensitivityVariable);


		KRATOS_CATCH("");
	}

	// ==============================================================================
	virtual void CalculateSensitivityGradient(Element& rElem,
                                	  const Variable<array_1d<double,3>>& rVariable,
                                      const Matrix& rDerivativesMatrix,
                                      Vector& rRHSContribution,
                                      ProcessInfo& rProcessInfo)
    {
		KRATOS_TRY;

        KRATOS_ERROR << "This should be implemented in the derived class." << std::endl;

		KRATOS_CATCH("");
	}

	// ==============================================================================
	virtual void CalculateSensitivityGradient(Element& rElem,
                                          const Variable<double>& rVariable,
                                          const Matrix& rDerivativesMatrix,
                                          Vector& rRHSContribution,
                                          ProcessInfo& rProcessInfo)
    {
		KRATOS_TRY;

        KRATOS_ERROR << "This should be implemented in the derived class." << std::endl;

		KRATOS_CATCH("");
	}

	// ==============================================================================
	virtual void CalculateSensitivityGradient(Condition& rCond,
                                	  const Variable<array_1d<double,3>>& rVariable,
                                      const Matrix& rDerivativesMatrix,
                                      Vector& rRHSContribution,
                                      ProcessInfo& rProcessInfo)
    {
		KRATOS_TRY;

        KRATOS_ERROR << "This should be implemented in the derived class." << std::endl;

		KRATOS_CATCH("");
	}

	// ==============================================================================
	virtual void CalculateSensitivityGradient(Condition& rCond,
                                          const Variable<double>& rVariable,
                                          const Matrix& rDerivativesMatrix,
                                          Vector& rRHSContribution,
                                          ProcessInfo& rProcessInfo)
    {
		KRATOS_TRY;

        KRATOS_ERROR << "This should be implemented in the derived class." << std::endl;

		KRATOS_CATCH("");
	}

	// ==============================================================================
	virtual void GetAdjointVariables(Element& rElem, Vector& rValues)
    {
		KRATOS_TRY;

		rElem.GetValuesVector(rValues);

		KRATOS_CATCH("");
	}

	// ==============================================================================
	virtual void GetAdjointVariables(Condition& rCond, Vector& rValues)
    {
		KRATOS_TRY;

		rCond.GetValuesVector(rValues);

		KRATOS_CATCH("");
	}

	// ==============================================================================
    void AssembleNodalSensitivityContribution(Variable<double> const& rSensitivityVariable,
                                              Vector const& rSensitivityVector,
                                              Element::GeometryType& rGeom) 
    {
        unsigned int index = 0;
        for (unsigned int i_node = 0; i_node < rGeom.PointsNumber(); ++i_node)
        {
                double& r_sensitivity =
                    rGeom[i_node].FastGetSolutionStepValue(rSensitivityVariable);
                rGeom[i_node].SetLock();
                r_sensitivity += rSensitivityVector[index++];
                rGeom[i_node].UnSetLock();
        }
    }

	// ==============================================================================
    void AssembleNodalSensitivityContribution(Variable<array_1d<double, 3>> const& rSensitivityVariable,
                                              Vector const& rSensitivityVector,
                                              Element::GeometryType& rGeom) 
    {
        unsigned int index = 0;
        for (unsigned int i_node = 0; i_node < rGeom.PointsNumber(); ++i_node)
        {
            array_1d<double, 3>& r_sensitivity =
                rGeom[i_node].FastGetSolutionStepValue(rSensitivityVariable);
            rGeom[i_node].SetLock();
            for (unsigned int d = 0; d < rGeom.WorkingSpaceDimension(); ++d)
                r_sensitivity[d] += rSensitivityVector[index++];
            rGeom[i_node].UnSetLock();
        }
	}

	// ==============================================================================
    void AssembleElementSensitivityContribution(Variable<double> const& rSensitivityVariable,
                                              Vector const& rSensitivityVector,
                                              Element::GeometryType& rGeom) 
    {
		//-->add code
	}

	// ==============================================================================
    void AssembleElementSensitivityContribution(Variable<array_1d<double, 3>> const& rSensitivityVariable,
                                              Vector const& rSensitivityVector,
                                              Element::GeometryType& rGeom)
    {
		//-->add code
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
    //std::string mSensitivityModelPartName;

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
	//      StructuralResponseFunction& operator=(SStructuralResponseFunction const& rOther);

	/// Copy constructor.
	//      StructuralResponseFunction(StructuralResponseFunction const& rOther);

	///@}

}; // Class StructuralResponseFunction

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

} // namespace Kratos.

#endif // STRUCTURAL_RESPONSE_FUNCTION_H

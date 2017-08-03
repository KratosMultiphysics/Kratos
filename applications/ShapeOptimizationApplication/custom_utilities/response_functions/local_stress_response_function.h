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
#include "response_function.h"


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

class LocalStressResponseFunction : ResponseFunction
{
public:
	///@name Type Definitions
	///@{



	typedef array_1d<double, 3> array_3d;

	

	/// Pointer definition of LocalStressResponseFunction
	KRATOS_CLASS_POINTER_DEFINITION(LocalStressResponseFunction);

	///@}
	///@name Life Cycle
	///@{

	/// Default constructor.
	LocalStressResponseFunction(ModelPart& model_part, Parameters& responseSettings)
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
	

		//get traced element
		m_id_of_traced_element = responseSettings["traced_element"].GetInt();
		m_traced_pElement = mr_model_part.pGetElement(m_id_of_traced_element);

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

		//move this to base class---------------------------------------------------------------
		//mSensitivityModelPartName = responseSettings["sensitivity_model_part_name"].GetString();

        Parameters nodal_sensitivity_variables = responseSettings["nodal_sensitivity_variables"];
        mNodalSensitivityVariables.resize(nodal_sensitivity_variables.size());
        for (unsigned int i = 0; i < nodal_sensitivity_variables.size(); ++i)
			mNodalSensitivityVariables[i] = nodal_sensitivity_variables[i].GetString();

		Parameters element_sensitivity_variables = responseSettings["element_sensitivity_variables"];
        mElementSensitivityVariables.resize(element_sensitivity_variables.size());
        for (unsigned int i = 0; i < element_sensitivity_variables.size(); ++i)
			mElementSensitivityVariables[i] = element_sensitivity_variables[i].GetString();	

		// maybe add also sensitivity variables	for conditions
		//----------------------------------------------------------------------------------------	

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

	// ==============================================================================
	void initialize()
	{
		KRATOS_TRY;

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
	void calculate_value() //---------------------------------------------------------------> rename it into CalculateValue()
	{
		KRATOS_TRY;

		// Working variables
		ProcessInfo &CurrentProcessInfo = mr_model_part.GetProcessInfo();

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
	void calculate_gradient()
	{
		KRATOS_TRY;

		// Formula computed in general notation:
		// \frac{dF}{dx} = eigenvector^T\cdot (frac{\partial RHS}{\partial x} - 
		//				   eigenvalue frac{\partial mass_matrix}{\partial x})\cdot eigenvector

		// First gradients are initialized
		array_3d zeros_array(3, 0.0);
		for (ModelPart::NodeIterator node_i = mr_model_part.NodesBegin(); node_i != mr_model_part.NodesEnd(); ++node_i)
			noalias(node_i->FastGetSolutionStepValue(LOCAL_STRESS_GRADIENT) )= zeros_array;

		// Gradient calculation is done by a semi-analytic approaches
		// The gradient is computed in one step

		switch (mGradientMode)
		{
		// Semi-analytic sensitivities
		case 1:
		{
			
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

		return m_stress_value;

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
			dFdX[node_i->Id()] = node_i->FastGetSolutionStepValue(LOCAL_STRESS_GRADIENT);

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
	void CalculateGradient(const Element& rElem, const Matrix& rLHS, Vector& rOutput, ProcessInfo& rProcessInfo)
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

	// ==============================================================================
	virtual void UpdateSensitivities() //---------------------------------------------move this to base class
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

		// Add condition sensitivity variables e.g. line load---------------------------------------

		//move this to base class
      
        KRATOS_CATCH("");
	}
	// ==============================================================================
	template <typename TDataType>
    void UpdateNodalSensitivities(Variable<TDataType> const& rSensitivityVariable) //---->move this to base class
	{
		KRATOS_TRY;

        ProcessInfo& r_process_info = mr_model_part.GetProcessInfo();
     
        Vector sensitivity_vector;
    	Vector response_gradient;
        Vector adjoint_vector;
        Matrix sensitivity_matrix;

		for (ModelPart::ElementIterator it = mr_model_part.ElementsBegin(); it != mr_model_part.ElementsEnd(); ++it)
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
            it->GetValuesVector(adjoint_vector);

            if (sensitivity_vector.size() != sensitivity_matrix.size1())
                sensitivity_vector.resize(sensitivity_matrix.size1(), false);

			// Compute the whole sensitivity
            noalias(sensitivity_vector) =
                    			(prod(sensitivity_matrix, adjoint_vector) +
                                    response_gradient);

            this->AssembleNodalSensitivityContribution(
                rSensitivityVariable, sensitivity_vector, r_geom);
        }

		//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		// also loop over conditions since e.g. load entries in RHS are changing when node coordinates are disturbed
		// and add sensitivity_vector from element and condition loop
		//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    
        mr_model_part.GetCommunicator().AssembleCurrentData(rSensitivityVariable);



		KRATOS_CATCH("");
	}

	// ==============================================================================
	template <typename TDataType>
	void UpdateElementSensitivities(Variable<TDataType> const& rSensitivityVariable)  //---->move this to base class
	{
		KRATOS_TRY;
		
        ProcessInfo& r_process_info = mr_model_part.GetProcessInfo();
     
        Vector sensitivity_vector;
    	Vector response_gradient;
        Vector adjoint_vector;
        Matrix sensitivity_matrix;

		for (ModelPart::ElementIterator it = mr_model_part.ElementsBegin(); it != mr_model_part.ElementsEnd(); ++it)
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
            it->GetValuesVector(adjoint_vector);

            if (sensitivity_vector.size() != sensitivity_matrix.size1())
               sensitivity_vector.resize(sensitivity_matrix.size1(), false);


			// Compute the whole sensitivity
            noalias(sensitivity_vector) =
                				(prod(sensitivity_matrix, adjoint_vector) +
                                 response_gradient);

               // this->AssembleNodalSensitivityContribution(
                  //  rSensitivityVariable, sensitivity_vector[k], r_geom);
        }
       

        mr_model_part.GetCommunicator().AssembleCurrentData(rSensitivityVariable);

		KRATOS_CATCH("");
	}

	// ==============================================================================
	void CalculateSensitivityGradient(Element& rElem,
                                	  const Variable<array_1d<double,3>>& rVariable,
                                      const Matrix& rDerivativesMatrix,
                                      Vector& rRHSContribution,
                                      ProcessInfo& rProcessInfo)
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
                                          ProcessInfo& rProcessInfo)
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
    void AssembleNodalSensitivityContribution(Variable<double> const& rSensitivityVariable,
                                              Vector const& rSensitivityVector,
                                              Element::GeometryType& rGeom) //----------------------------> move this to base class
    {
        unsigned int index = 0;
        for (unsigned int i_node = 0; i_node < rGeom.PointsNumber(); ++i_node)
        {
           // if (rGeom[i_node].GetValue(UPDATE_SENSITIVITIES) == true)
            //{
                double& r_sensitivity =
                    rGeom[i_node].FastGetSolutionStepValue(rSensitivityVariable);
                rGeom[i_node].SetLock();
                r_sensitivity += rSensitivityVector[index++];
                rGeom[i_node].UnSetLock();
            //}
            //else
               // ++index;
        }
    }

	// ==============================================================================
    void AssembleNodalSensitivityContribution(Variable<array_1d<double, 3>> const& rSensitivityVariable,
                                              Vector const& rSensitivityVector,
                                              Element::GeometryType& rGeom) //----------------------------> move this to base class
    {
        unsigned int index = 0;
        for (unsigned int i_node = 0; i_node < rGeom.PointsNumber(); ++i_node)
        {
           // if (rGeom[i_node].GetValue(UPDATE_SENSITIVITIES) == true)
            //{
                array_1d<double, 3>& r_sensitivity =
                    rGeom[i_node].FastGetSolutionStepValue(rSensitivityVariable);
                rGeom[i_node].SetLock();
                for (unsigned int d = 0; d < rGeom.WorkingSpaceDimension(); ++d)
                    r_sensitivity[d] += rSensitivityVector[index++];
                rGeom[i_node].UnSetLock();
            //}
            //else
              //  index += rGeom.WorkingSpaceDimension();
        }
	}

	// ==============================================================================
    void AssembleElementSensitivityContribution(Variable<double> const& rSensitivityVariable,
                                              Vector const& rSensitivityVector,
                                              Element::GeometryType& rGeom) //----------------------------> move this to base class
    {
	}

	// ==============================================================================
    void AssembleElementSensitivityContribution(Variable<array_1d<double, 3>> const& rSensitivityVariable,
                                              Vector const& rSensitivityVector,
                                              Element::GeometryType& rGeom) //----------------------------> move this to base class
    {
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
	double m_stress_value; 
	double mDelta;
	double m_initial_value;
	bool m_initial_value_defined;
	unsigned int m_id_of_traced_element;
	int m_id_of_location;
	Element::Pointer m_traced_pElement;
	std::string m_traced_stress_type;
	std::string m_stress_treatment;

	//move this to base class-----------------------------
    std::string mSensitivityModelPartName;
    std::vector<std::string> mNodalSensitivityVariables;
	std::vector<std::string> mElementSensitivityVariables;
	//-----------------------------------------------------



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

#endif // EIGENFRQUENCY_RESPONSE_FUNCTION_H

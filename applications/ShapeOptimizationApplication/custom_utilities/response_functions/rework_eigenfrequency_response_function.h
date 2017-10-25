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

#ifndef REWORK_EIGENFREQUENCY_RESPONSE_FUNCTION_H
#define REWORK_EIGENFREQUENCY_RESPONSE_FUNCTION_H

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

class ReworkEigenfrequencyResponseFunction : public StructuralResponseFunction
{
public:
	///@name Type Definitions
	///@{

	typedef StructuralResponseFunction BaseType;
	typedef array_1d<double, 3> array_3d;

	// TODO solve this via template or how to get this from Eigensolverstrategy 
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;

	typedef typename LocalSpaceType::VectorType DenseVectorType;
	typedef typename LocalSpaceType::MatrixType DenseMatrixType;
	typedef Variable<DenseVectorType> VariableDenseVectorType;
	typedef Variable<DenseMatrixType> VariableDenseMatrixType;

	

	/// Pointer definition of ReworkEigenfrequencyResponseFunction
	KRATOS_CLASS_POINTER_DEFINITION(ReworkEigenfrequencyResponseFunction);

	///@}
	///@name Life Cycle
	///@{

	/// Default constructor.
	ReworkEigenfrequencyResponseFunction(ModelPart& model_part, Parameters& responseSettings)
	: StructuralResponseFunction(model_part, responseSettings)
	{
	
		// Initialize member variables to NULL
		m_initial_value = 0.0;
		m_initial_value_defined = false;
		m_current_response_value = 0.0;

		// Get number of eigenfrequency for which the structure has to be optimized
		m_traced_eigenvalue = responseSettings["traced_eigenfrequency"].GetInt();

	}

	/// Destructor.
	virtual ~ReworkEigenfrequencyResponseFunction()
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

		ModelPart& r_model_part = this->GetModelPart();

		const VariableDenseVectorType& rEIGENVALUE_VECTOR =
            KratosComponents<VariableDenseVectorType>::Get("EIGENVALUE_VECTOR");

		int num_of_computed_eigenvalues = (r_model_part.GetProcessInfo()[rEIGENVALUE_VECTOR]).size();

		if(num_of_computed_eigenvalues < m_traced_eigenvalue)
			KRATOS_THROW_ERROR(std::runtime_error, "The chosen eigenvalue was not solved by the eigenvalue analysis!", "");

		m_current_response_value = 	(r_model_part.GetProcessInfo()[rEIGENVALUE_VECTOR])[m_traced_eigenvalue - 1]; 

		// Change sign of response: only maximization makes sense in case of eigenfrequency optimization		
		m_current_response_value *= (-1.0); // TODO: do this also in SA?

		// Set initial value if not done yet
		if(!m_initial_value_defined)
		{
			m_initial_value = m_current_response_value;
			m_initial_value_defined = true;
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
		return "ReworkEigenfrequencyResponseFunction";
	}

	/// Print information about this object.
	virtual void PrintInfo(std::ostream &rOStream) const
	{
		rOStream << "ReworkEigenfrequencyResponseFunction";
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
	/*void UpdateSensitivities() override
	{
		KRATOS_TRY;

		BaseType::UpdateSensitivities();

		KRATOS_CATCH("");
	}*/

	void CalculateGradient(const Element& rAdjointElem,
                                   const Matrix& rAdjointMatrix,
                                   Vector& rResponseGradient,
                                   ProcessInfo& rProcessInfo) override
    {
        KRATOS_TRY;

		KRATOS_ERROR << "I am wrong here. There is no adjoint problem to solve for eigenfrequency as response." << std::endl;

        KRATOS_CATCH("");
    }

    void CalculateGradient(const Condition& rAdjointCondition,
                                   const Matrix& rAdjointMatrix,
                                   Vector& rResponseGradient,
                                   ProcessInfo& rProcessInfo) override
    {
        KRATOS_TRY;

		KRATOS_ERROR << "I am wrong here. There is no adjoint problem to solve for eigenfrequency as response." << std::endl;

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

	template <typename TDataType>
    void UpdateNodalSensitivities(Variable<TDataType> const& rSensitivityVariable) override
	{
		KRATOS_TRY;

	    ModelPart& r_model_part = this->GetModelPart();
        ProcessInfo& r_process_info = r_model_part.GetProcessInfo();
        //double delta_time = -r_process_info[DELTA_TIME];
        const int num_threads = OpenMPUtils::GetNumThreads();
        std::vector<Vector> sensitivity_vector(num_threads);
        std::vector<Vector> response_gradient(num_threads);
        std::vector<Vector> adjoint_vector(num_threads);
        std::vector<Matrix> sensitivity_matrix(num_threads);

        Communicator& r_comm = r_model_part.GetCommunicator();
        if (r_comm.TotalProcesses() > 1)
        {
            // here we make sure we only add the old sensitivity once
            // when we assemble.
#pragma omp parallel
            {
                ModelPart::NodeIterator nodes_begin;
                ModelPart::NodeIterator nodes_end;
                OpenMPUtils::PartitionedIterators(r_model_part.Nodes(), nodes_begin, nodes_end);
                for (auto it = nodes_begin; it != nodes_end; ++it)
                    if (it->FastGetSolutionStepValue(PARTITION_INDEX) != r_comm.MyPID())
                        it->FastGetSolutionStepValue(rSensitivityVariable) =
                            rSensitivityVariable.Zero();//------------------------------->What is here done??
            }
        }
        // Assemble element contributions.
#pragma omp parallel
        {
            ModelPart::ElementIterator elements_begin;
            ModelPart::ElementIterator elements_end;
            OpenMPUtils::PartitionedIterators(r_model_part.Elements(),
                                              elements_begin, elements_end);
            int k = OpenMPUtils::ThisThread();

            for (auto it = elements_begin; it != elements_end; ++it)
            {
                //std::cout << ("I compute now sensitivities of element #") << it->Id() << std::endl;
                Element::GeometryType& r_geom = it->GetGeometry();
                bool update_sensitivities = false;
                for (unsigned int i_node = 0; i_node < r_geom.PointsNumber(); ++i_node)
                    if (r_geom[i_node].GetValue(UPDATE_SENSITIVITIES) == true)
                    {
                        update_sensitivities = true;
                        break;
                    }

                if (update_sensitivities == false) // true for most elements
                    continue;

                // --> calculate here sensitivities

                this->AssembleNodalSensitivityContribution(
                    rSensitivityVariable, sensitivity_vector[k], r_geom);  //----> check for correct output
            }
        }

//         Assemble condition contributions.
 #pragma omp parallel
        {
            ModelPart::ConditionIterator conditions_begin;
            ModelPart::ConditionIterator conditions_end;
            OpenMPUtils::PartitionedIterators(r_model_part.Conditions(),
                                               conditions_begin, conditions_end);
            int k = OpenMPUtils::ThisThread();

            for (auto it = conditions_begin; it != conditions_end; ++it)
            {
                Condition::GeometryType& r_geom = it->GetGeometry();
                bool update_sensitivities = false;
                for (unsigned int i_node = 0; i_node < r_geom.PointsNumber(); ++i_node)
                    if (r_geom[i_node].GetValue(UPDATE_SENSITIVITIES) == true)
                    {
                        update_sensitivities = true;
                        break;
                    }

                if (update_sensitivities == false)
                    continue;

			    // --> calculate here sensitivities

                this->AssembleNodalSensitivityContribution(
                    rSensitivityVariable, sensitivity_vector[k], r_geom);	//----> check for correct output
            }        
        }
    
        r_model_part.GetCommunicator().AssembleCurrentData(rSensitivityVariable);

		KRATOS_CATCH("");
	}

	// ==============================================================================
	template <typename TDataType>
	void UpdateElementSensitivities(Variable<TDataType> const& rSensitivityVariable,
	Variable<TDataType> const& rOutputVariable)  override
	{
		KRATOS_TRY;

        ModelPart& r_model_part = this->GetModelPart();
        ProcessInfo& r_process_info = r_model_part.GetProcessInfo();
        //double delta_time = -r_process_info[DELTA_TIME];
        const int num_threads = OpenMPUtils::GetNumThreads();
        std::vector<Vector> sensitivity_vector(num_threads);
        std::vector<Vector> response_gradient(num_threads);
        std::vector<Vector> adjoint_vector(num_threads);
        std::vector<Matrix> sensitivity_matrix(num_threads);

  
        //std::cout << ("I compute now element sensitivities") << std::endl;
	#pragma omp parallel
        {
            ModelPart::ElementIterator elements_begin;
            ModelPart::ElementIterator elements_end;
            OpenMPUtils::PartitionedIterators(r_model_part.Elements(),
                                              elements_begin, elements_end);
            int k = OpenMPUtils::ThisThread();

            //std::cout << ("Inertia sensitivities:")  << std::endl;

            for (auto it = elements_begin; it != elements_end; ++it)
            {
                if (it->GetValue(UPDATE_SENSITIVITIES) == true)
                {
         
					// --> calculate here sensitivities
			    
                    this->AssembleElementSensitivityContribution(
                  	        rOutputVariable, sensitivity_vector[k], *it);		//----> check for correct output
              
                }
            }
        }
       
        r_model_part.GetCommunicator().AssembleCurrentData(rSensitivityVariable);

		KRATOS_CATCH("");
	}

	// ==============================================================================
	template <typename TDataType>
    void UpdateConditionSensitivities(Variable<TDataType> const& rSensitivityVariable, 
	Variable<TDataType> const& rOutputVariable) override
	{
		KRATOS_TRY;

	    ModelPart& r_model_part = this->GetModelPart();
        ProcessInfo& r_process_info = r_model_part.GetProcessInfo();
        //double delta_time = -r_process_info[DELTA_TIME];
        const int num_threads = OpenMPUtils::GetNumThreads();
        std::vector<Vector> sensitivity_vector(num_threads);
        std::vector<Vector> response_gradient(num_threads);
        std::vector<Vector> adjoint_vector(num_threads);
        std::vector<Matrix> sensitivity_matrix(num_threads);

	//         Assemble condition contributions.
#pragma omp parallel
        {
            ModelPart::ConditionIterator conditions_begin;
            ModelPart::ConditionIterator conditions_end;
            OpenMPUtils::PartitionedIterators(r_model_part.Conditions(),
                                               conditions_begin, conditions_end);
            int k = OpenMPUtils::ThisThread();

            for (auto it = conditions_begin; it != conditions_end; ++it)
            {
            
                if (it->GetValue(UPDATE_SENSITIVITIES) == true)
                {   
					// --> calculate here sensitivities

		            Condition::GeometryType& r_geom = it->GetGeometry();
                    this->AssembleConditionSensitivityContribution(
               		         rOutputVariable, sensitivity_vector[k], r_geom); //----> check for correct output

                }
            }
        }
    
        r_model_part.GetCommunicator().AssembleCurrentData(rSensitivityVariable);

		KRATOS_CATCH("");
	}	

	// ==============================================================================
	void CalculateSensitivityGradient(Element& rAdjointElem,
                                              const Variable<array_1d<double,3>>& rVariable,
                                              const Matrix& rDerivativesMatrix,
                                              Vector& rResponseGradient,
                                              ProcessInfo& rProcessInfo) override
    {
      	KRATOS_TRY;

      	KRATOS_ERROR << "I am wrong here. There is no partial derivative to compute." << std::endl;  

     	 KRATOS_CATCH("");
	}

	// ==============================================================================
	void CalculateSensitivityGradient(Element& rAdjointElem,
                                              const Variable<double>& rVariable,
                                              const Matrix& rDerivativesMatrix,
                                              Vector& rResponseGradient,
                                              ProcessInfo& rProcessInfo) override
    {
      	KRATOS_TRY;

		KRATOS_ERROR << "I am wrong here. There is no partial derivative to compute." << std::endl;  

        KRATOS_CATCH("");
	}

	// ==============================================================================
	void CalculateSensitivityGradient(Condition& rAdjointCondition,
                                              const Variable<array_1d<double,3>>& rVariable,
                                              const Matrix& rDerivativesMatrix,
                                              Vector& rResponseGradient,
                                              ProcessInfo& rProcessInfo) override
    {
		KRATOS_TRY;

		KRATOS_ERROR << "I am wrong here. There is no partial derivative to compute." << std::endl;  

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

		KRATOS_ERROR << "I am wrong here. There is no partial derivative to compute." << std::endl;  

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
	double m_initial_value;
	bool m_initial_value_defined;
	int m_traced_eigenvalue;

	///@}
///@name Private Operators
	///@{

	///@}
	///@name Private Operations
	///@{
// --------------------------------------------------------------------------
	double GetEigenvalue(int id_eigenvalue)
	{

		KRATOS_TRY;

		ModelPart& r_model_part = this->GetModelPart();

		double current_eigenvalue = 0.0;

		const VariableDenseVectorType& rEIGENVALUE_VECTOR =
            KratosComponents<VariableDenseVectorType>::Get("EIGENVALUE_VECTOR");

		int num_of_computed_eigenvalues = (r_model_part.GetProcessInfo()[rEIGENVALUE_VECTOR]).size();

		if(num_of_computed_eigenvalues < id_eigenvalue)
			KRATOS_THROW_ERROR(std::runtime_error, "The chosen eigenvalue was not solved by the eigenvalue analysis!", "");

		current_eigenvalue = (r_model_part.GetProcessInfo()[rEIGENVALUE_VECTOR])[id_eigenvalue-1]; 

		return current_eigenvalue;
		
		KRATOS_CATCH("");

	}

	// --------------------------------------------------------------------------
	Vector GetEigenvectorOfElement(ModelPart::ElementIterator traced_element, int id_eigenvalue, int size_of_eigenvector)
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
	//      ReworkEigenfrequencyResponseFunction& operator=(ReworkEigenfrequencyResponseFunction const& rOther);

	/// Copy constructor.
	//      ReworkEigenfrequencyResponseFunction(ReworkEigenfrequencyResponseFunction const& rOther);

	///@}

}; // Class ReworkEigenfrequencyResponseFunction

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

} // namespace Kratos.

#endif // REWORK_EIGENFREQUENCY_RESPONSE_FUNCTION_H

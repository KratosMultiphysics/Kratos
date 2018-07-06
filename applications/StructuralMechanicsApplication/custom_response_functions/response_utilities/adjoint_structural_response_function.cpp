// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Martin Fusseder, https://github.com/MFusseder 
//   



// System includes

// External includes

// Project includes
#include "utilities/openmp_utils.h"
#include "adjoint_structural_response_function.h"
#include "utilities/variable_utils.h"


// Application includes

namespace Kratos
{

    /// Constructor.
    AdjointStructuralResponseFunction::AdjointStructuralResponseFunction(ModelPart& rModelPart, Parameters ResponseSettings)
      : mrModelPart(rModelPart)
    {
        KRATOS_TRY;

        mSensitivityModelPartName = ResponseSettings["sensitivity_model_part_name"].GetString();

        this->ReadDesignVariables(mNodalSensitivityScalarVariables, mNodalSensitivityVectorVariables, ResponseSettings["nodal_sensitivity_variables"]);
        this->ReadDesignVariables(mElementSensitivityScalarVariables, mElementSensitivityVectorVariables, ResponseSettings["element_sensitivity_variables"]);
        this->ReadDesignVariables(mConditionSensitivityScalarVariables, mConditionSensitivityVectorVariables, ResponseSettings["condition_sensitivity_variables"]);
        
        // Set gradient mode
        const std::string gradient_mode = ResponseSettings["gradient_mode"].GetString();

        // Mode 1: semi-analytic sensitivities
        if (gradient_mode == "semi_analytic")
        {
            mGradientMode = 1;
            double delta = ResponseSettings["step_size"].GetDouble();
            mDelta = delta;
        }
        else
            KRATOS_ERROR << "Specified gradient_mode not recognized. The only option is: semi_analytic. Specified gradient_mode: " <<  gradient_mode << std::endl;


        KRATOS_CATCH("");
    }

    /// Destructor.
    AdjointStructuralResponseFunction::~AdjointStructuralResponseFunction()
    {
    }

    ModelPart& AdjointStructuralResponseFunction::GetModelPart()
    {
      return mrModelPart;
    }

    ModelPart& AdjointStructuralResponseFunction::GetModelPart() const 
    {
      return mrModelPart;
    }

    void AdjointStructuralResponseFunction::Initialize()
    {
        KRATOS_TRY;

        this->Clear();
        this->Check();

        ModelPart& r_model_part = this->GetModelPart();
        ModelPart& r_sensitivity_model_part = r_model_part.GetSubModelPart(mSensitivityModelPartName);

        // Initialize flags.
        VariableUtils().SetNonHistoricalVariable(UPDATE_SENSITIVITIES, true, r_sensitivity_model_part.Nodes());
        VariableUtils().SetNonHistoricalVariable(UPDATE_SENSITIVITIES, true, r_sensitivity_model_part.Elements());
        VariableUtils().SetNonHistoricalVariable(UPDATE_SENSITIVITIES, true, r_sensitivity_model_part.Conditions());

    if(mGradientMode == 1)
    {
        VariableUtils().SetNonHistoricalVariable(DISTURBANCE_MEASURE, mDelta, r_model_part.Elements());
        VariableUtils().SetNonHistoricalVariable(DISTURBANCE_MEASURE, mDelta, r_model_part.Conditions());
    }

        KRATOS_CATCH("");
    }

    void AdjointStructuralResponseFunction::InitializeSolutionStep(){}

    void AdjointStructuralResponseFunction::FinalizeSolutionStep()
    {
        KRATOS_TRY;

        this->UpdateSensitivities();

        KRATOS_CATCH("");
    }

    void AdjointStructuralResponseFunction::Check()
    {
        KRATOS_TRY;

        ModelPart& r_model_part = this->GetModelPart();

        KRATOS_ERROR_IF_NOT(r_model_part.HasSubModelPart(mSensitivityModelPartName))
            << "No sub model part \"" << mSensitivityModelPartName << "\"" << std::endl;

        KRATOS_CATCH("");
    }

    void AdjointStructuralResponseFunction::Clear()
    {
        KRATOS_TRY;

        ModelPart& r_model_part = this->GetModelPart();

        // Reset flags.
        VariableUtils().SetNonHistoricalVariable(UPDATE_SENSITIVITIES, false, r_model_part.Nodes());
        VariableUtils().SetNonHistoricalVariable(UPDATE_SENSITIVITIES, false, r_model_part.Elements());
        VariableUtils().SetNonHistoricalVariable(UPDATE_SENSITIVITIES, false, r_model_part.Conditions());

        // Set nodal sensitivity result variables to zero.
        for (const auto& variable_pair : mNodalSensitivityScalarVariables)
            VariableUtils().SetToZero_ScalarVar(variable_pair[1], r_model_part.Nodes()); 
        for (const auto& variable_pair : mNodalSensitivityVectorVariables)
            VariableUtils().SetToZero_VectorVar(variable_pair[1], r_model_part.Nodes());
        // Set elemental sensitivity result variables to zero.
        for (const auto& variable_pair : mElementSensitivityScalarVariables)
        {
            #pragma omp parallel for
            for (int i = 0; i< static_cast<int> (r_model_part.Elements().size()); ++i) 
            {
                ElementsContainerType::iterator it = r_model_part.ElementsBegin() + i;
                it->SetValue(variable_pair[1], variable_pair[1].Zero());
            }
        }
        for (const auto& variable_pair : mElementSensitivityVectorVariables)
        {
             #pragma omp parallel for
            for (int i = 0; i< static_cast<int> (r_model_part.Elements().size()); ++i) 
            {
                ElementsContainerType::iterator it = r_model_part.ElementsBegin() + i;
                it->SetValue(variable_pair[1], variable_pair[1].Zero());
            }   
        }
        // Set conditional sensitivity result variables to zero.
        for (const auto& variable_pair : mConditionSensitivityScalarVariables)
        {
            #pragma omp parallel for
            for (int i = 0; i< static_cast<int> (r_model_part.Conditions().size()); ++i) 
            {
                ConditionsContainerType::iterator it = r_model_part.ConditionsBegin() + i;
                const unsigned int number_of_nodes = it->GetGeometry().size();
                for(unsigned int j = 0; j < number_of_nodes; ++j)
                    it->GetGeometry()[j].FastGetSolutionStepValue(variable_pair[1]) = variable_pair[1].Zero();
            }
        }
        for (const auto& variable_pair : mConditionSensitivityVectorVariables)
        {
            #pragma omp parallel for
            for (int i = 0; i< static_cast<int> (r_model_part.Conditions().size()); ++i) 
            {
                ConditionsContainerType::iterator it = r_model_part.ConditionsBegin() + i;
                const unsigned int number_of_nodes = it->GetGeometry().size();
                for(unsigned int j = 0; j < number_of_nodes; ++j)
                    it->GetGeometry()[j].FastGetSolutionStepValue(variable_pair[1]) = variable_pair[1].Zero();
            }  
        }

        KRATOS_CATCH("");
    }

    /// Calculate the local gradient w.r.t. primal solution.
    /**
     * @param[in]     rAdjointElem      the adjoint element.
     * @param[in]     rAdjointMatrix    the transposed gradient of the
     *                                  element's residual w.r.t. primal.
     * @param[out]    rResponseGradient the gradient of the response function.
     * @param[in]     rProcessInfo      the current process info.
     */
    void AdjointStructuralResponseFunction::CalculateGradient(const Element& rAdjointElem,
                                   const Matrix& rAdjointMatrix,
                                   Vector& rResponseGradient,
                                   ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;

        if (rResponseGradient.size() != rAdjointMatrix.size1())
            rResponseGradient.resize(rAdjointMatrix.size1(), false);

        rResponseGradient.clear();

        KRATOS_CATCH("");
    }

    void AdjointStructuralResponseFunction::CalculateGradient(const Condition& rAdjointCondition,
                                   const Matrix& rAdjointMatrix,
                                   Vector& rResponseGradient,
                                   ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;

        if (rResponseGradient.size() != rAdjointMatrix.size1())
            rResponseGradient.resize(rAdjointMatrix.size1(), false);

        rResponseGradient.clear();

        KRATOS_CATCH("");
    }

    /// Calculate the local gradient w.r.t. first derivatives of primal solution.
    /**
     * @param[in]     rAdjointElem      the adjoint element.
     * @param[in]     rAdjointMatrix    the transposed gradient of the
     *                                  element's residual w.r.t. first derivatives.
     * @param[out]    rResponseGradient the gradient of the response function.
     * @param[in]     rProcessInfo      the current process info.
     */
    void AdjointStructuralResponseFunction::CalculateFirstDerivativesGradient(const Element& rAdjointElem,
                                                   const Matrix& rAdjointMatrix,
                                                   Vector& rResponseGradient,
                                                   ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;

        if (rResponseGradient.size() != rAdjointMatrix.size1())
            rResponseGradient.resize(rAdjointMatrix.size1(), false);

        rResponseGradient.clear();

        KRATOS_CATCH("");
    }

    void AdjointStructuralResponseFunction::CalculateFirstDerivativesGradient(const Condition& rAdjointCondition,
                                                   const Matrix& rAdjointMatrix,
                                                   Vector& rResponseGradient,
                                                   ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;
        if (rResponseGradient.size() != rAdjointMatrix.size1())
            rResponseGradient.resize(rAdjointMatrix.size1(), false);
        rResponseGradient.clear();
        KRATOS_CATCH("");
    }

    /// Calculate the local gradient w.r.t. second derivatives of primal solution.
    /**
     * @param[in]     rAdjointElem      the adjoint element.
     * @param[in]     rAdjointMatrix    the transposed gradient of the
     *                                  element's residual w.r.t. second derivatives.
     * @param[out]    rResponseGradient the gradient of the response function.
     * @param[in]     rProcessInfo      the current process info.
     */
    void AdjointStructuralResponseFunction::CalculateSecondDerivativesGradient(const Element& rAdjointElem,
                                                    const Matrix& rAdjointMatrix,
                                                    Vector& rResponseGradient,
                                                    ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;

        if (rResponseGradient.size() != rAdjointMatrix.size1())
            rResponseGradient.resize(rAdjointMatrix.size1(), false);

        rResponseGradient.clear();

        KRATOS_CATCH("");
    }

    void AdjointStructuralResponseFunction::CalculateSecondDerivativesGradient(const Condition& rAdjointCondition,
                                                    const Matrix& rAdjointMatrix,
                                                    Vector& rResponseGradient,
                                                    ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;
        if (rResponseGradient.size() != rAdjointMatrix.size1())
            rResponseGradient.resize(rAdjointMatrix.size1(), false);
        rResponseGradient.clear();
        KRATOS_CATCH("");
    }

    void AdjointStructuralResponseFunction::UpdateSensitivities()
    {
        KRATOS_TRY;

        for (const auto& variable_pair : mNodalSensitivityScalarVariables)
            this->UpdateNodalSensitivities(variable_pair[0], variable_pair[1]);
        for (const auto& variable_pair : mNodalSensitivityVectorVariables)
            this->UpdateNodalSensitivities(variable_pair[0], variable_pair[1]);
        for (const auto& variable_pair : mElementSensitivityScalarVariables)
            this->UpdateElementSensitivities(variable_pair[0], variable_pair[1]);
        for (const auto& variable_pair : mElementSensitivityVectorVariables)
            this->UpdateElementSensitivities(variable_pair[0], variable_pair[1]);
        for (const auto& variable_pair : mConditionSensitivityScalarVariables)
            this->UpdateConditionSensitivities(variable_pair[0], variable_pair[1]);
        for (const auto& variable_pair : mConditionSensitivityVectorVariables)
            this->UpdateConditionSensitivities(variable_pair[0], variable_pair[1]);

        KRATOS_CATCH("");
    }


    /// Calculate the scalar valued response function
    double AdjointStructuralResponseFunction::CalculateValue(ModelPart& rModelPart)
    {
        return 0.0;
    }

    template <typename TDataType>
    void AdjointStructuralResponseFunction::UpdateNodalSensitivities(Variable<TDataType> const& rSensitivityVariable, Variable<TDataType> const& rOutputVariable)
    {
        KRATOS_TRY;

        ModelPart& r_model_part = this->GetModelPart();
        ProcessInfo& r_process_info = r_model_part.GetProcessInfo();
        const int num_threads = 1; 
        std::vector<Vector> sensitivity_vector(num_threads);
        std::vector<Vector> response_gradient(num_threads);
        std::vector<Vector> adjoint_vector(num_threads);
        std::vector<Matrix> sensitivity_matrix(num_threads);

        int k = 0;

        for (ModelPart::ElementIterator it = r_model_part.ElementsBegin(); it != r_model_part.ElementsEnd(); ++it)
        {
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

            // Compute the pseudo load
            it->CalculateSensitivityMatrix(
                rSensitivityVariable, sensitivity_matrix[k], r_process_info);

            if(sensitivity_matrix[k].size1() > 0)
            {  
                // This part of the sensitivity is computed from the objective
                // with primal variables treated as constant.
                this->CalculateSensitivityGradient(
                    *it, rSensitivityVariable, sensitivity_matrix[k],
                    response_gradient[k], r_process_info);

                // Get the adjoint displacement field
                it->GetValuesVector(adjoint_vector[k]);

                if (sensitivity_vector[k].size() != sensitivity_matrix[k].size1())
                    sensitivity_vector[k].resize(sensitivity_matrix[k].size1(), false);

                // Compute the whole sensitivity
                noalias(sensitivity_vector[k]) = (prod(sensitivity_matrix[k], adjoint_vector[k]) +
                                response_gradient[k]);

                this->AssembleNodalSensitivityContribution(
                    rOutputVariable, sensitivity_vector[k], r_geom); 
            }
        }
    
        // Assemble condition contributions.
        for (ModelPart::ConditionIterator it = r_model_part.ConditionsBegin(); it != r_model_part.ConditionsEnd(); ++it)
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

            // This is multiplied with the adjoint to compute sensitivity
            // contributions from the condition.
            it->CalculateSensitivityMatrix(
                rSensitivityVariable, sensitivity_matrix[k], r_process_info);

            if(sensitivity_matrix[k].size1() > 0)
            {  
                // This part of the sensitivity is computed from the objective
                // with primal variables treated as constant.
                this->CalculateSensitivityGradient(
                    *it, rSensitivityVariable, sensitivity_matrix[k],
                    response_gradient[k], r_process_info);

                // Get the adjoint displacement field
                it->GetValuesVector(adjoint_vector[k]);

                if (sensitivity_vector[k].size() != sensitivity_matrix[k].size1())
                    sensitivity_vector[k].resize(sensitivity_matrix[k].size1(), false);

                // Compute the whole sensitivity
                noalias(sensitivity_vector[k]) = (prod(sensitivity_matrix[k], adjoint_vector[k]) + response_gradient[k]);

                this->AssembleNodalSensitivityContribution(rOutputVariable, sensitivity_vector[k], r_geom);	
            }
        }

        r_model_part.GetCommunicator().AssembleCurrentData(rSensitivityVariable);

        KRATOS_CATCH("");
    }

    // ==============================================================================
    template <typename TDataType>
    void AdjointStructuralResponseFunction::UpdateElementSensitivities(Variable<TDataType> const& rSensitivityVariable, Variable<TDataType> const& rOutputVariable)
    {
        KRATOS_TRY;

        ModelPart& r_model_part = this->GetModelPart();
        ProcessInfo& r_process_info = r_model_part.GetProcessInfo();
        const int num_threads = OpenMPUtils::GetNumThreads();
        std::vector<Vector> sensitivity_vector(num_threads);
        std::vector<Vector> response_gradient(num_threads);
        std::vector<Vector> adjoint_vector(num_threads);
        std::vector<Matrix> sensitivity_matrix(num_threads);

        #pragma omp parallel for
        for (int i = 0; i< static_cast<int> (r_model_part.Elements().size()); ++i) 
        {
            const unsigned int k = OpenMPUtils::ThisThread();
            ElementsContainerType::iterator it = r_model_part.ElementsBegin() + i;

            if (it->GetValue(UPDATE_SENSITIVITIES) == true)
            {
                // Compute the pseudo load
                it->CalculateSensitivityMatrix(
                    rSensitivityVariable, sensitivity_matrix[k], r_process_info);

                if(sensitivity_matrix[k].size1() > 0)
                {    
                    // This part of the sensitivity is computed from the objective
                    // with primal variables treated as constant.
                    this->CalculateSensitivityGradient(
                        *it, rSensitivityVariable, sensitivity_matrix[k],
                            response_gradient[k], r_process_info);
                    
                  
                    // Get the adjoint displacement field
                    it->GetValuesVector(adjoint_vector[k]);

                    if (sensitivity_vector[k].size() != sensitivity_matrix[k].size1())
                        sensitivity_vector[k].resize(sensitivity_matrix[k].size1(), false);

                    // Compute the whole sensitivity
                    noalias(sensitivity_vector[k]) = (prod(sensitivity_matrix[k], adjoint_vector[k]) +
                                    response_gradient[k]);

                    this->AssembleElementSensitivityContribution(
                                rOutputVariable, sensitivity_vector[k], *it);		
                }
            }
            
        }

        r_model_part.GetCommunicator().AssembleCurrentData(rSensitivityVariable);

        KRATOS_CATCH("");
    }

    // ==============================================================================
    template <typename TDataType>
    void AdjointStructuralResponseFunction::UpdateConditionSensitivities(Variable<TDataType> const& rSensitivityVariable, Variable<TDataType> const& rOutputVariable)
    {
        KRATOS_TRY;

        ModelPart& r_model_part = this->GetModelPart();
        ProcessInfo& r_process_info = r_model_part.GetProcessInfo();
        const int num_threads = OpenMPUtils::GetNumThreads();
        std::vector<Vector> sensitivity_vector(num_threads);
        std::vector<Vector> response_gradient(num_threads);
        std::vector<Vector> adjoint_vector(num_threads);
        std::vector<Matrix> sensitivity_matrix(num_threads);

        //  Assemble condition contributions.
        #pragma omp parallel for
        for (int i = 0; i< static_cast<int> (r_model_part.Conditions().size()); ++i) 
        {
            const unsigned int k = OpenMPUtils::ThisThread();
            ConditionsContainerType::iterator it = r_model_part.ConditionsBegin() + i;

            if (it->GetValue(UPDATE_SENSITIVITIES) == true)
            {
                // Compute the pseudo load
                it->CalculateSensitivityMatrix(
                    rSensitivityVariable, sensitivity_matrix[k], r_process_info);

                if(sensitivity_matrix[k].size1() > 0)
                {      
                    // This part of the sensitivity is computed from the objective
                    // with primal variables treated as constant.
                    this->CalculateSensitivityGradient(
                        *it, rSensitivityVariable, sensitivity_matrix[k],
                        response_gradient[k], r_process_info);

                    // Get the adjoint displacement field
                    it->GetValuesVector(adjoint_vector[k]);

                    if (sensitivity_vector[k].size() != sensitivity_matrix[k].size1())
                        sensitivity_vector[k].resize(sensitivity_matrix[k].size1(), false);

                    // Compute the whole sensitivity
                    noalias(sensitivity_vector[k]) = (prod(sensitivity_matrix[k], adjoint_vector[k]) +
                                    response_gradient[k]);

                    Condition::GeometryType& r_geom = it->GetGeometry();
                    this->AssembleConditionSensitivityContribution(
                                rOutputVariable, sensitivity_vector[k], r_geom); 
                }
            }
        }

        r_model_part.GetCommunicator().AssembleCurrentData(rSensitivityVariable);

        KRATOS_CATCH("");
    }



    /// Calculate the local gradient of response function w.r.t. the sensitivity variable.
    /**
     * @param[in]     rAdjointElem       the adjoint element.
     * @param[in]     rVariable          the sensitivity variable.
     * @param[in]     rDerivativesMatrix the transposed gradient of the element's
     *                                   residual w.r.t. the sensitivity variable.
     * @param[out]    rResponseGradient  the gradient of the response function.
     * @param[in,out] rProcessInfo       the current process info.
     */
    void AdjointStructuralResponseFunction::CalculateSensitivityGradient(Element& rAdjointElem, 
                                              const Variable<double>& rVariable,
                                              const Matrix& rDerivativesMatrix,
                                              Vector& rResponseGradient,
                                              ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;

        KRATOS_ERROR << "This should be implemented in the derived class." << std::endl;

        KRATOS_CATCH("");
    }

    void AdjointStructuralResponseFunction::CalculateSensitivityGradient(Condition& rAdjointCondition,
                                              const Variable<double>& rVariable,
                                              const Matrix& rDerivativesMatrix,
                                              Vector& rResponseGradient,
                                              ProcessInfo& rProcessInfo)
    {
         KRATOS_TRY;

         KRATOS_ERROR << "This should be implemented in the derived class." << std::endl;

         KRATOS_CATCH("");
    }

    /// Calculate the local gradient of response function w.r.t. the sensitivity variable.
    /**
     * @param[in]     rAdjointElem       the adjoint element.
     * @param[in]     rVariable          the sensitivity variable.
     * @param[in]     rDerivativesMatrix the transposed gradient of the element's
     *                                   residual w.r.t. the sensitivity variable.
     * @param[out]    rResponseGradient  the gradient of the response function.
     * @param[in,out] rProcessInfo       the current process info.
     */
    void AdjointStructuralResponseFunction::CalculateSensitivityGradient(Element& rAdjointElem, 
                                              const Variable<array_1d<double,3>>& rVariable,
                                              const Matrix& rDerivativesMatrix,
                                              Vector& rResponseGradient,
                                              ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;

        KRATOS_ERROR << "This should be implemented in the derived class." << std::endl;

        KRATOS_CATCH("");
    }

    void AdjointStructuralResponseFunction::CalculateSensitivityGradient(Condition& rAdjointCondition, 
                                              const Variable<array_1d<double,3>>& rVariable,
                                              const Matrix& rDerivativesMatrix,
                                              Vector& rResponseGradient,
                                              ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;

        KRATOS_ERROR << "This should be implemented in the derived class." << std::endl;

        KRATOS_CATCH("");
    }

    void AdjointStructuralResponseFunction::AssembleNodalSensitivityContribution(Variable<double> const& rSensitivityVariable,
                                              Vector const& rSensitivityVector,
                                              Element::GeometryType& rGeom)
    {
        unsigned int index = 0;
        for (unsigned int i_node = 0; i_node < rGeom.PointsNumber(); ++i_node)
        {
            if (rGeom[i_node].GetValue(UPDATE_SENSITIVITIES) == true)
            {
                double& r_sensitivity =
                    rGeom[i_node].FastGetSolutionStepValue(rSensitivityVariable);
                rGeom[i_node].SetLock();
                r_sensitivity += rSensitivityVector[index++];
                rGeom[i_node].UnSetLock();
            }
            else
                ++index;
        }
    }

    void AdjointStructuralResponseFunction::AssembleNodalSensitivityContribution(Variable<array_1d<double, 3>> const& rSensitivityVariable,
                                              Vector const& rSensitivityVector,
                                              Element::GeometryType& rGeom)
    {
        unsigned int index = 0;
        for (unsigned int i_node = 0; i_node < rGeom.PointsNumber(); ++i_node)
        {
            if (rGeom[i_node].GetValue(UPDATE_SENSITIVITIES) == true)
            {
                array_1d<double, 3>& r_sensitivity =
                    rGeom[i_node].FastGetSolutionStepValue(rSensitivityVariable);
                rGeom[i_node].SetLock();
                for (unsigned int d = 0; d < rGeom.WorkingSpaceDimension(); ++d)
                    r_sensitivity[d] += rSensitivityVector[index++];
                rGeom[i_node].UnSetLock();
            }
            else
                index += rGeom.WorkingSpaceDimension();
        }
    }

    void AdjointStructuralResponseFunction::AssembleElementSensitivityContribution(Variable<double> const& rSensitivityVariable,
                                                Vector const& rSensitivityVector,
                                                Element& rElem)
    {
        rElem.SetValue(rSensitivityVariable , rSensitivityVector[0]);
        // attention: one has to ensure that element is able to print the variable type later on his Gauss-Points
    }

    void AdjointStructuralResponseFunction::AssembleElementSensitivityContribution(Variable<array_1d<double, 3>> const& rSensitivityVariable,
                                                Vector const& rSensitivityVector,
                                                Element& rElem)
    {
        rElem.SetValue(rSensitivityVariable , rSensitivityVector);
        // attention: one has to ensure that element is able to print the variable type later on his Gauss-Points
    }

    void AdjointStructuralResponseFunction::AssembleConditionSensitivityContribution(Variable<double> const& rSensitivityVariable,
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

    void AdjointStructuralResponseFunction::AssembleConditionSensitivityContribution(Variable<array_1d<double, 3>> const& rSensitivityVariable,
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

    void AdjointStructuralResponseFunction::ReadDesignVariables(std::vector<std::vector<Variable<double>>>& rScalarDesignVariables, 
        std::vector<std::vector<Variable<array_1d<double,3>>>>& rVectorDesignVariables, Parameters DesignVariableSettings)
    {
        for (unsigned int i = 0; i < DesignVariableSettings.size(); ++i)
        {
            const std::string variable_label = DesignVariableSettings[i].GetString();
            const std::string output_variable_label = variable_label + "_SENSITIVITY";
            std::vector<Variable<double>> helper_scalar_variables;
            std::vector<Variable<array_1d<double,3>>> helper_vector_variables;

            if (KratosComponents<Variable<double>>::Has(variable_label))
            {
                const Variable<double>& r_variable =
                    KratosComponents<Variable<double>>::Get(variable_label);

                helper_scalar_variables.push_back(r_variable); 

                if (KratosComponents<Variable<double>>::Has(output_variable_label))
                {
                    const Variable<double>& r_output_variable =
                        KratosComponents<Variable<double>>::Get(output_variable_label); 
                    helper_scalar_variables.push_back(r_output_variable); 
                }
                else
                    KRATOS_ERROR << "Unsupported output variable: " << output_variable_label << "." << std::endl;  

                rScalarDesignVariables.push_back(helper_scalar_variables);           
            }
            else if (KratosComponents<Variable<array_1d<double,3>>>::Has(variable_label))
            {
                const Variable<array_1d<double, 3>>& r_variable =
                    KratosComponents<Variable<array_1d<double, 3>>>::Get(variable_label);

                helper_vector_variables.push_back(r_variable); 

                if (KratosComponents<Variable<array_1d<double,3>>>::Has(output_variable_label))
                {
                    const Variable<array_1d<double, 3>>& r_output_variable =
                        KratosComponents<Variable<array_1d<double, 3>>>::Get(output_variable_label); 
                    helper_vector_variables.push_back(r_output_variable); 
                }  
                else
                    KRATOS_ERROR << "Unsupported output variable: " << output_variable_label << "." << std::endl; 

                rVectorDesignVariables.push_back(helper_vector_variables);      
            }
            else
                KRATOS_ERROR << "Unsupported variable: " << variable_label << "." << std::endl;  
        }
    }

};



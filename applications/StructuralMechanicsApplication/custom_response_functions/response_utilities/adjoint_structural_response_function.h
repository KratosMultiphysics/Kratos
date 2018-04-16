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

#ifndef ADJOINT_STRUCTURAL_RESPONSE_FUNCTION_H
#define ADJOINT_STRUCTURAL_RESPONSE_FUNCTION_H

// System includes
#include <vector>
#include <string>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/condition.h"
#include "includes/process_info.h"
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"
#include "includes/ublas_interface.h"
#include "utilities/openmp_utils.h"


// Application includes

namespace Kratos
{
///@addtogroup StructuralMechanicsApplication
///@{

///@name Kratos Classes
///@{

/// A base class for response functions.
class AdjointStructuralResponseFunction
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(AdjointStructuralResponseFunction);

    typedef VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>> VariableComponentType;
    typedef Variable<array_1d<double, 3>> VariableWithComponentsType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    AdjointStructuralResponseFunction(ModelPart& rModelPart, Parameters& rParameters)
      : mrModelPart(rModelPart)
    {
        KRATOS_TRY;

        mSensitivityModelPartName =
           rParameters["sensitivity_model_part_name"].GetString();

        Parameters nodal_sensitivity_variables = rParameters["nodal_sensitivity_variables"];
        mNodalSensitivityVariables.resize(nodal_sensitivity_variables.size());
        for (unsigned int i = 0; i < nodal_sensitivity_variables.size(); ++i)
            mNodalSensitivityVariables[i] = nodal_sensitivity_variables[i].GetString();

        Parameters element_sensitivity_variables = rParameters["element_sensitivity_variables"];
        mElementSensitivityVariables.resize(element_sensitivity_variables.size());
        for (unsigned int i = 0; i < element_sensitivity_variables.size(); ++i)
            mElementSensitivityVariables[i] = element_sensitivity_variables[i].GetString();

        Parameters condition_sensitivity_variables = rParameters["condition_sensitivity_variables"];
        mConditionSensitivityVariables.resize(condition_sensitivity_variables.size());
        for (unsigned int i = 0; i < condition_sensitivity_variables.size(); ++i)
            mConditionSensitivityVariables[i] = condition_sensitivity_variables[i].GetString();

        // Set gradient mode
        std::string gradient_mode = rParameters["gradient_mode"].GetString();

        // Mode 1: semi-analytic sensitivities
        if (gradient_mode.compare("semi_analytic") == 0)
        {
            mGradientMode = 1;
            double delta = rParameters["step_size"].GetDouble();
            mDelta = delta;
        }
        else
            KRATOS_ERROR << "Specified gradient_mode not recognized. The only option is: semi_analytic. Specified gradient_mode: " <<  gradient_mode << std::endl;


        KRATOS_CATCH("");
    }

    /// Destructor.
    virtual ~AdjointStructuralResponseFunction()
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
      return mrModelPart;
    }

    const ModelPart& GetModelPart() const
    {
      return mrModelPart;
    }

    virtual void Initialize()
    {
        KRATOS_TRY;

        this->Clear();
        this->Check();

        ModelPart& r_model_part = this->GetModelPart();

        // Initialize flags.
#pragma omp parallel
        {
            ModelPart::NodeIterator nodes_begin;
            ModelPart::NodeIterator nodes_end;
            OpenMPUtils::PartitionedIterators(
                r_model_part.GetSubModelPart(mSensitivityModelPartName).Nodes(),
                nodes_begin, nodes_end);
            for (auto it = nodes_begin; it != nodes_end; ++it)
                it->SetValue(UPDATE_SENSITIVITIES, true);
        }
#pragma omp parallel
        {
            ModelPart::ElementIterator elements_begin;
            ModelPart::ElementIterator elements_end;
            OpenMPUtils::PartitionedIterators(
                r_model_part.GetSubModelPart(mSensitivityModelPartName).Elements(),
                elements_begin, elements_end);
            for (auto it = elements_begin; it != elements_end; ++it)
                it->SetValue(UPDATE_SENSITIVITIES, true);

        }
#pragma omp parallel
        {
            ModelPart::ConditionIterator conditions_begin;
            ModelPart::ConditionIterator conditions_end;
            OpenMPUtils::PartitionedIterators(
                r_model_part.GetSubModelPart(mSensitivityModelPartName).Conditions(),
                conditions_begin, conditions_end);
            for (auto it = conditions_begin; it != conditions_end; ++it)
                it->SetValue(UPDATE_SENSITIVITIES, true);

        }

    if(mGradientMode == 1)
    {
#pragma omp parallel
        {
            ModelPart::ElementIterator elements_begin;
            ModelPart::ElementIterator elements_end;
            OpenMPUtils::PartitionedIterators(r_model_part.Elements(), elements_begin, elements_end);
            for (auto it = elements_begin; it != elements_end; ++it)
                it->SetValue(DISTURBANCE_MEASURE, mDelta);
        }
#pragma omp parallel
        {
            ModelPart::ConditionIterator conditions_begin;
            ModelPart::ConditionIterator conditions_end;
            OpenMPUtils::PartitionedIterators(r_model_part.Conditions(), conditions_begin, conditions_end);
            for (auto it = conditions_begin; it != conditions_end; ++it)
                it->SetValue(DISTURBANCE_MEASURE, mDelta);
        }
    }

        KRATOS_CATCH("");
    }

    virtual void InitializeSolutionStep()
    {
    }

    virtual void FinalizeSolutionStep()
    {
        KRATOS_TRY;

        this->UpdateSensitivities();

        KRATOS_CATCH("");
    }

    virtual void Check()
    {
        KRATOS_TRY;

        ModelPart& r_model_part = this->GetModelPart();

        KRATOS_ERROR_IF_NOT(r_model_part.HasSubModelPart(mSensitivityModelPartName))
            << "No sub model part \"" << mSensitivityModelPartName << "\"" << std::endl;

        KRATOS_CATCH("");
    }

    virtual void Clear()
    {
        KRATOS_TRY;

        ModelPart& r_model_part = this->GetModelPart();

        // Reset flags.
#pragma omp parallel
        {
            ModelPart::NodeIterator nodes_begin;
            ModelPart::NodeIterator nodes_end;
            OpenMPUtils::PartitionedIterators(r_model_part.Nodes(), nodes_begin, nodes_end);
            for (auto it = nodes_begin; it != nodes_end; ++it)
                it->SetValue(UPDATE_SENSITIVITIES, false);
        }
#pragma omp parallel
        {
            ModelPart::ElementIterator elements_begin;
            ModelPart::ElementIterator elements_end;
            OpenMPUtils::PartitionedIterators(r_model_part.Elements(), elements_begin, elements_end);
            for (auto it = elements_begin; it != elements_end; ++it)
                it->SetValue(UPDATE_SENSITIVITIES, false);
        }
#pragma omp parallel
        {
            ModelPart::ConditionIterator conditions_begin;
            ModelPart::ConditionIterator conditions_end;
            OpenMPUtils::PartitionedIterators(r_model_part.Conditions(), conditions_begin, conditions_end);
            for (auto it = conditions_begin; it != conditions_end; ++it)
                it->SetValue(UPDATE_SENSITIVITIES, false);
        }


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
                    OpenMPUtils::PartitionedIterators(r_model_part.Nodes(),
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
                    OpenMPUtils::PartitionedIterators(r_model_part.Nodes(),
                                                      nodes_begin, nodes_end);
                    for (auto it = nodes_begin; it != nodes_end; ++it)
                        it->FastGetSolutionStepValue(r_variable) = r_variable.Zero();
                }
            }
            else
                KRATOS_ERROR << "Unsupported variable: " << label << "." << std::endl;
        }

        // Set sensitivity variables to zero.
        for (auto label : mElementSensitivityVariables)
        {
            std::string output_variable_label = label + "_SENSITIVITY";
            if (KratosComponents<Variable<double>>::Has(output_variable_label) == true)
            {
                const Variable<double>& r_variable =
                    KratosComponents<Variable<double>>::Get(output_variable_label);

#pragma omp parallel
                {
                    ModelPart::ElementIterator elements_begin;
                    ModelPart::ElementIterator elements_end;
                    OpenMPUtils::PartitionedIterators(r_model_part.Elements(),
                                                      elements_begin, elements_end);
                    for (auto it = elements_begin; it != elements_end; ++it)
                        it->SetValue(r_variable, r_variable.Zero());

                }
            }
            else if (KratosComponents<Variable<array_1d<double, 3>>>::Has(output_variable_label) == true)
            {
                const Variable<array_1d<double, 3>>& r_variable =
                    KratosComponents<Variable<array_1d<double, 3>>>::Get(output_variable_label);

#pragma omp parallel
                {
                    ModelPart::ElementIterator elements_begin;
                    ModelPart::ElementIterator elements_end;
                    OpenMPUtils::PartitionedIterators(r_model_part.Elements(),
                                                      elements_begin, elements_end);
                    for (auto it = elements_begin; it != elements_end; ++it)
                        it->SetValue(r_variable, r_variable.Zero());
                }
            }
            else
                KRATOS_ERROR << "Unsupported variable: " << output_variable_label << "." << std::endl;
        }

        // Set sensitivity variables to zero.
        for (auto label : mConditionSensitivityVariables)
        {
            std::string output_variable_label = label + "_SENSITIVITY";
            if (KratosComponents<Variable<double>>::Has(output_variable_label) == true)
            {
                const Variable<double>& r_variable =
                    KratosComponents<Variable<double>>::Get(output_variable_label);

#pragma omp parallel
                {
                    ModelPart::ConditionIterator conditions_begin;
                    ModelPart::ConditionIterator conditions_end;
                    OpenMPUtils::PartitionedIterators(r_model_part.Conditions(),
                                                      conditions_begin, conditions_end);
                    for (auto it = conditions_begin; it != conditions_end; ++it)
                    {
                        //loop nodes
                        const unsigned int NumberOfNodes = it->GetGeometry().size();
                        for(unsigned int i = 0; i < NumberOfNodes; ++i)
                            it->GetGeometry()[i].FastGetSolutionStepValue(r_variable) = r_variable.Zero();
                    }

                }
            }
            else if (KratosComponents<Variable<array_1d<double, 3>>>::Has(output_variable_label) == true)
            {
                const Variable<array_1d<double, 3>>& r_variable =
                    KratosComponents<Variable<array_1d<double, 3>>>::Get(output_variable_label);

#pragma omp parallel
                {
                    ModelPart::ConditionIterator conditions_begin;
                    ModelPart::ConditionIterator conditions_end;
                    OpenMPUtils::PartitionedIterators(r_model_part.Conditions(),
                                                      conditions_begin, conditions_end);
                    for (auto it = conditions_begin; it != conditions_end; ++it)
                    {
                        //loop nodes
                        const unsigned int NumberOfNodes = it->GetGeometry().size();
                        for(unsigned int i = 0; i < NumberOfNodes; ++i)
                            it->GetGeometry()[i].FastGetSolutionStepValue(r_variable) = r_variable.Zero();
                    }
                }
            }
            else
                KRATOS_ERROR << "Unsupported variable: " << label << "." << std::endl;
        }
        //TODO: Check if it is a good idea to write the conditional sensitivities to the nodes?!

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
    virtual void CalculateGradient(const Element& rAdjointElem,
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

    virtual void CalculateGradient(const Condition& rAdjointCondition,
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
    virtual void CalculateFirstDerivativesGradient(const Element& rAdjointElem,
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

    // virtual void CalculateFirstDerivativesGradient(const Condition& rAdjointCondition,
    //                                                const Matrix& rAdjointMatrix,
    //                                                Vector& rResponseGradient,
    //                                                ProcessInfo& rProcessInfo)
    // {
    //     KRATOS_TRY;

    //     if (rResponseGradient.size() != rAdjointMatrix.size1())
    //         rResponseGradient.resize(rAdjointMatrix.size1(), false);

    //     rResponseGradient.clear();

    //     KRATOS_CATCH("");
    // }

    /// Calculate the local gradient w.r.t. second derivatives of primal solution.
    /**
     * @param[in]     rAdjointElem      the adjoint element.
     * @param[in]     rAdjointMatrix    the transposed gradient of the
     *                                  element's residual w.r.t. second derivatives.
     * @param[out]    rResponseGradient the gradient of the response function.
     * @param[in]     rProcessInfo      the current process info.
     */
    virtual void CalculateSecondDerivativesGradient(const Element& rAdjointElem,
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

    // virtual void CalculateSecondDerivativesGradient(const Condition& rAdjointCondition,
    //                                                 const Matrix& rAdjointMatrix,
    //                                                 Vector& rResponseGradient,
    //                                                 ProcessInfo& rProcessInfo)
    // {
    //     KRATOS_TRY;

    //     if (rResponseGradient.size() != rAdjointMatrix.size1())
    //         rResponseGradient.resize(rAdjointMatrix.size1(), false);

    //     rResponseGradient.clear();

    //     KRATOS_CATCH("");
    // }

    virtual void UpdateSensitivities()
    {
        KRATOS_TRY;

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
                KRATOS_ERROR << "Unsupported variable: " << label << "." << std::endl;
        }

        // Elemental sensitivity variables e.g. 2nd moment of inertia---------------------------------
        for (auto label : mElementSensitivityVariables)
        {
            //create label for output variable
            std::string output_label = label + "_SENSITIVITY";

            if (KratosComponents<Variable<double>>::Has(label) == true)
            {
                const Variable<double>& r_variable =
                    KratosComponents<Variable<double>>::Get(label);

                //check for output variable to save later the computed sensitivity
                if( KratosComponents<Variable<double>>::Has(output_label) == true )
                {
                    const Variable<double>& r_output_variable =
                        KratosComponents<Variable<double>>::Get(output_label);
                    this->UpdateElementSensitivities(r_variable, r_output_variable);
                }
                else
                    KRATOS_ERROR << "Unsupported element variable for output: " << output_label << "." << std::endl;


            }
            else if (KratosComponents<Variable<array_1d<double,3>>>::Has(label) == true)
            {
                const Variable<array_1d<double,3>>& r_variable =
                    KratosComponents<Variable<array_1d<double,3>>>::Get(label);

                //check for output variable to save later the computed sensitivity
                if (KratosComponents<Variable<array_1d<double,3>>>::Has(output_label) == true)
                {
                    const Variable<array_1d<double,3>>& r_output_variable =
                        KratosComponents<Variable<array_1d<double,3>>>::Get(output_label);
                    this->UpdateElementSensitivities(r_variable, r_output_variable);
                }
                else
                    KRATOS_ERROR << "Unsupported element variable for output: " << output_label << "." << std::endl;


            }
            else
                KRATOS_ERROR << "Unsupported element variable: " << label << "." << std::endl;
        }

        for (auto label : mConditionSensitivityVariables)
        {
            //create label for output variable
            std::string output_label = label + "_SENSITIVITY";

            if (KratosComponents<Variable<double>>::Has(label) == true)
            {
                const Variable<double>& r_variable =
                    KratosComponents<Variable<double>>::Get(label);

                //check for output variable to save later the computed sensitivity
                if( KratosComponents<Variable<double>>::Has(output_label) == true )
                {
                    const Variable<double>& r_output_variable =
                        KratosComponents<Variable<double>>::Get(output_label);
                    this->UpdateConditionSensitivities(r_variable, r_output_variable);
                }
                else
                    KRATOS_ERROR << "Unsupported condition variable for output: " << output_label << "." << std::endl;
            }
            else if (KratosComponents<Variable<array_1d<double,3>>>::Has(label) == true)
            {
                const Variable<array_1d<double,3>>& r_variable =
                    KratosComponents<Variable<array_1d<double,3>>>::Get(label);

                //check for output variable to save later the computed sensitivity
                if (KratosComponents<Variable<array_1d<double,3>>>::Has(output_label) == true)
                {
                    const Variable<array_1d<double,3>>& r_output_variable =
                        KratosComponents<Variable<array_1d<double,3>>>::Get(output_label);
                    this->UpdateConditionSensitivities(r_variable, r_output_variable);
                }
                else
                    KRATOS_ERROR << "Unsupported condition variable for output: " << output_label << "." << std::endl;

            }
            else
                KRATOS_ERROR << "Unsupported condition variable: " << label << "." << std::endl;
        }

        KRATOS_CATCH("");
    }


    /// Calculate the scalar valued response function
    virtual double CalculateValue(ModelPart& rModelPart)
    {
        return 0.0;
    }

    ///@}

protected:
    ///@name Protected member Variables
    ///@{

    ModelPart& mrModelPart;

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    /// Calculate and add the sensitivity from the current adjoint and primal solutions.
    /*
     * This function updates the accumulated sensitivities for a single time step.
     * The accumulated sensitivity is defined as:
     *
     * \f[
     * d_{\mathbf{s}}\bar{J} = \Sigma_{n=1}^N
     *   (\partial_{\mathbf{s}}J^n + \lambda^{nT}\partial_{\mathbf{s}}\mathbf{r}^n)
     *    \Delta t
     * \f]
     *
     * with \f$\mathbf{r}^n\f$ the residual of the governing partial differential
     * equation for the current step. This function should be called once for
     * each design variable per time step.
     */
   /* template <typename TDataType>
    void UpdateNodalSensitivities(Variable<TDataType> const& rSensitivityVariable)
    {
        KRATOS_TRY

        ModelPart& r_model_part = this->GetModelPart();
        ProcessInfo& r_process_info = r_model_part.GetProcessInfo();
        double delta_time = -r_process_info[DELTA_TIME];
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
                            rSensitivityVariable.Zero();
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

                // This is multiplied with the adjoint to compute sensitivity
                // contributions from the element.
                it->CalculateSensitivityMatrix(
                    rSensitivityVariable, sensitivity_matrix[k], r_process_info);

                // This is the sensitivity contribution coming from the response
                // function with primal variables treated as constant.
                this->CalculateSensitivityGradient(
                    *it, rSensitivityVariable, sensitivity_matrix[k],
                    response_gradient[k], r_process_info);

                // Get the element adjoint vector.
                it->GetValuesVector(adjoint_vector[k]);

                if (sensitivity_vector[k].size() != sensitivity_matrix[k].size1())
                    sensitivity_vector[k].resize(sensitivity_matrix[k].size1(), false);

                // Calculated the total sensitivity contribution for the element.
                noalias(sensitivity_vector[k]) =
                    delta_time * (prod(sensitivity_matrix[k], adjoint_vector[k]) +
                                  response_gradient[k]);

                this->AssembleNodalSensitivityContribution(
                    rSensitivityVariable, sensitivity_vector[k], r_geom);
            }
        }*/
//         // Assemble condition contributions.
// #pragma omp parallel
//         {
//             ModelPart::ConditionIterator conditions_begin;
//             ModelPart::ConditionIterator conditions_end;
//             OpenMPUtils::PartitionedIterators(r_model_part.Conditions(),
//                                               conditions_begin, conditions_end);
//             int k = OpenMPUtils::ThisThread();

//             for (auto it = conditions_begin; it != conditions_end; ++it)
//             {
//                 Condition::GeometryType& r_geom = it->GetGeometry();
//                 bool update_sensitivities = false;
//                 for (unsigned int i_node = 0; i_node < r_geom.PointsNumber(); ++i_node)
//                     if (r_geom[i_node].GetValue(UPDATE_SENSITIVITIES) == true)
//                     {
//                         update_sensitivities = true;
//                         break;
//                     }

//                 if (update_sensitivities == false)
//                     continue;

//                 // This is multiplied with the adjoint to compute sensitivity
//                 // contributions from the condition.
//                 it->CalculateSensitivityMatrix(
//                     rSensitivityVariable, sensitivity_matrix[k], r_process_info);

//                 // This is the sensitivity contribution coming from the response
//                 // function with primal variables treated as constant.
//                 this->CalculateSensitivityGradient(
//                     *it, rSensitivityVariable, sensitivity_matrix[k],
//                     response_gradient[k], r_process_info);

//                 // Get the element adjoint vector.
//                 it->GetValuesVector(adjoint_vector[k]);

//                 if (sensitivity_vector[k].size() != sensitivity_matrix[k].size1())
//                     sensitivity_vector[k].resize(sensitivity_matrix[k].size1(), false);

//                 // Calculated the total sensitivity contribution for the condition.
//                 noalias(sensitivity_vector[k]) =
//                     delta_time * (prod(sensitivity_matrix[k], adjoint_vector[k]) +
//                                   response_gradient[k]);

//                 this->AssembleNodalSensitivityContribution(
//                     rSensitivityVariable, sensitivity_vector[k], r_geom);
//             }
//         }

 //       r_model_part.GetCommunicator().AssembleCurrentData(rSensitivityVariable);

   //     KRATOS_CATCH("")
    //}
    template <typename TDataType>
    void UpdateNodalSensitivities(Variable<TDataType> const& rSensitivityVariable)
    {
        KRATOS_TRY;

        ModelPart& r_model_part = this->GetModelPart();
        ProcessInfo& r_process_info = r_model_part.GetProcessInfo();
        //double delta_time = -r_process_info[DELTA_TIME];
        const int num_threads = 1; //OpenMPUtils::GetNumThreads();
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
/*#pragma omp parallel
        {
            ModelPart::ElementIterator elements_begin;
            ModelPart::ElementIterator elements_end;
            OpenMPUtils::PartitionedIterators(r_model_part.Elements(),
                                              elements_begin, elements_end);
            int k = OpenMPUtils::ThisThread();*/

            int k = 0;

            //for (auto it = elements_begin; it != elements_end; ++it)
            for (ModelPart::ElementIterator it = r_model_part.ElementsBegin(); it != r_model_part.ElementsEnd(); ++it)
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
                    noalias(sensitivity_vector[k]) =
                                    /*delta_time * */(prod(sensitivity_matrix[k], adjoint_vector[k]) +
                                    response_gradient[k]);

                    this->AssembleNodalSensitivityContribution(
                        rSensitivityVariable, sensitivity_vector[k], r_geom);  //----> check for correct output
                }
            }
        //}

//         Assemble condition contributions.
/*#pragma omp parallel
        {
            ModelPart::ConditionIterator conditions_begin;
            ModelPart::ConditionIterator conditions_end;
            OpenMPUtils::PartitionedIterators(r_model_part.Conditions(),
                                               conditions_begin, conditions_end);
            int k = OpenMPUtils::ThisThread();*/

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
                    noalias(sensitivity_vector[k]) =
                                    /*delta_time * */(prod(sensitivity_matrix[k], adjoint_vector[k]) +
                                        response_gradient[k]);


                    /*for(unsigned i = 0; i < sensitivity_vector.size(); i++)
                    {
                        std::cout << ("SA result conditions = ") << i << sensitivity_vector[i] << std::endl;
                    }*/

                    this->AssembleNodalSensitivityContribution(
                        rSensitivityVariable, sensitivity_vector[k], r_geom);	//----> check for correct output
                }
            }
        //}

        r_model_part.GetCommunicator().AssembleCurrentData(rSensitivityVariable);

        KRATOS_CATCH("");
    }

    // ==============================================================================
    template <typename TDataType>
    void UpdateElementSensitivities(Variable<TDataType> const& rSensitivityVariable, Variable<TDataType> const& rOutputVariable)
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
                        noalias(sensitivity_vector[k]) =
                                    /*delta_time * */(prod(sensitivity_matrix[k], adjoint_vector[k]) +
                                    response_gradient[k]);

                        //std::cout << ("element sensitivty = ") << sensitivity_vector[0] << std::endl;
                        this->AssembleElementSensitivityContribution(
                                  rOutputVariable, sensitivity_vector[k], *it);		//----> check for correct output

                        //std::cout <<  sensitivity_vector[k][0] << std::endl;
                    }
                }
            }
        }

        r_model_part.GetCommunicator().AssembleCurrentData(rSensitivityVariable);

        KRATOS_CATCH("");
    }

    // ==============================================================================
    template <typename TDataType>
    void UpdateConditionSensitivities(Variable<TDataType> const& rSensitivityVariable, Variable<TDataType> const& rOutputVariable)
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
                        noalias(sensitivity_vector[k]) =
                                /*delta_time * */(prod(sensitivity_matrix[k], adjoint_vector[k]) +
                                    response_gradient[k]);

                        /*for(unsigned i = 0; i < sensitivity_vector.size(); i++)
                        {
                            std::cout << ("SA result = ") << sensitivity_vector[i] << std::endl;
                        }*/

                        Condition::GeometryType& r_geom = it->GetGeometry();
                        this->AssembleConditionSensitivityContribution(
                                   rOutputVariable, sensitivity_vector[k], r_geom); //----> check for correct output
                    }
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
    virtual void CalculateSensitivityGradient(Element& rAdjointElem, //make const Element?
                                              const Variable<double>& rVariable,
                                              const Matrix& rDerivativesMatrix,
                                              Vector& rResponseGradient,
                                              ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;

        KRATOS_ERROR << "This should be implemented in the derived class." << std::endl;

        KRATOS_CATCH("");
    }

    virtual void CalculateSensitivityGradient(Condition& rAdjointCondition, //make const Condition?
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
    virtual void CalculateSensitivityGradient(Element& rAdjointElem, //make const Element?
                                              const Variable<array_1d<double,3>>& rVariable,
                                              const Matrix& rDerivativesMatrix,
                                              Vector& rResponseGradient,
                                              ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;

        KRATOS_ERROR << "This should be implemented in the derived class." << std::endl;

        KRATOS_CATCH("");
    }

    virtual void CalculateSensitivityGradient(Condition& rAdjointCondition, //make const Condition?
                                              const Variable<array_1d<double,3>>& rVariable,
                                              const Matrix& rDerivativesMatrix,
                                              Vector& rResponseGradient,
                                              ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;

        KRATOS_ERROR << "This should be implemented in the derived class." << std::endl;

        KRATOS_CATCH("");
    }

    void AssembleNodalSensitivityContribution(Variable<double> const& rSensitivityVariable,
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

    void AssembleNodalSensitivityContribution(Variable<array_1d<double, 3>> const& rSensitivityVariable,
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

    void AssembleElementSensitivityContribution(Variable<double> const& rSensitivityVariable,
                                                Vector const& rSensitivityVector,
                                                Element& rElem)
    {
        rElem.SetValue(rSensitivityVariable , rSensitivityVector[0]);
        // attention: one has to ensure that element is able to print the variable type later on his Gauss-Points
    }

    void AssembleElementSensitivityContribution(Variable<array_1d<double, 3>> const& rSensitivityVariable,
                                                Vector const& rSensitivityVector,
                                                Element& rElem)
    {
        rElem.SetValue(rSensitivityVariable , rSensitivityVector);
        // attention: one has to ensure that element is able to print the variable type later on his Gauss-Points
    }

    void AssembleConditionSensitivityContribution(Variable<double> const& rSensitivityVariable,
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

    void AssembleConditionSensitivityContribution(Variable<array_1d<double, 3>> const& rSensitivityVariable,
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

    ///@}

private:
    ///@name Member Variables
    ///@{

    std::string mSensitivityModelPartName;
    std::vector<std::string> mNodalSensitivityVariables;
    std::vector<std::string> mElementSensitivityVariables;
    std::vector<std::string> mConditionSensitivityVariables;
    unsigned int mGradientMode;
    double mDelta;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    ///@}
};

///@} // Kratos Classes

///@} //Structural Mechanics Application group

} /* namespace Kratos.*/

#endif /* KRATOS_STRUCTURAL_RESPONSE_FUNCTION defined */

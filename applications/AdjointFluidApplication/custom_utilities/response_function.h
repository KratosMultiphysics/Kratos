//  KratosAdjointFluidApplication
//
//  License:		 BSD License
//					 license: AdjointFluidApplication/license.txt
//
//  Main authors:    Michael Andre, https://github.com/msandre
//

#if !defined(KRATOS_RESPONSE_FUNCTION)
#define KRATOS_RESPONSE_FUNCTION

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
///@addtogroup AdjointFluidApplication
///@{

///@name Kratos Classes
///@{

/// A base class for response functions.
class ResponseFunction
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(ResponseFunction);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    ResponseFunction(ModelPart& rModelPart, Parameters& rParameters)
      : mrModelPart(rModelPart)
    {
        KRATOS_TRY;

        mSensitivityModelPartName =
            rParameters["sensitivity_model_part_name"].GetString();
        Parameters nodal_sensitivity_variables = rParameters["nodal_sensitivity_variables"];
        mNodalSensitivityVariables.resize(nodal_sensitivity_variables.size());
        for (unsigned int i = 0; i < nodal_sensitivity_variables.size(); ++i)
            mNodalSensitivityVariables[i] = nodal_sensitivity_variables[i].GetString();
        
        mResponseName = "";
        mOutputFilename = "";

        KRATOS_CATCH("");
    }

    /// Destructor.
    virtual ~ResponseFunction()
    {
        if (mIsOutputFileOpened)
            mOutputFileStream.close();
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

    void SetOutputFilename(const std::string OutputFilename)
    {
        this->mOutputFilename = OutputFilename;
    }

    void SetResponseName(const std::string ResponseName)
    {
        this->mResponseName = ResponseName;
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

        if (mOutputFilename.compare("") != 0) 
        {
            mOutputFileStream.open(mOutputFilename + ".data");
            mOutputFileStream<<"#time       "<<mResponseName<<std::endl;
            mIsOutputFileOpened = true;
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

        double response_value = this->CalculateValue( this->GetModelPart() );
        
        if (mIsOutputFileOpened) 
        {
            ProcessInfo& rProcessInfo = this->GetModelPart().GetProcessInfo();
            mOutputFileStream.precision(5);
            mOutputFileStream<<std::scientific<<rProcessInfo[TIME]<<" ";
            mOutputFileStream.precision(15);
            mOutputFileStream<<std::scientific<<response_value<<std::endl;        
        }

        KRATOS_CATCH("");
    }

    virtual void Check()
    {
        KRATOS_TRY;

        ModelPart& r_model_part = this->GetModelPart();

        if (r_model_part.HasSubModelPart(mSensitivityModelPartName) == false)
            KRATOS_ERROR << "No sub model part \"" << mSensitivityModelPartName
                         << "\"" << std::endl;

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

    // virtual void CalculateGradient(const Condition& rAdjointCondition,
    //                                const Matrix& rAdjointMatrix,
    //                                Vector& rResponseGradient,
    //                                ProcessInfo& rProcessInfo)
    // {
    //     KRATOS_TRY;

    //     if (rResponseGradient.size() != rAdjointMatrix.size1())
    //         rResponseGradient.resize(rAdjointMatrix.size1(), false);

    //     rResponseGradient.clear();

    //     KRATOS_CATCH("");
    // }

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
    template <typename TDataType>
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
        }
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

        r_model_part.GetCommunicator().AssembleCurrentData(rSensitivityVariable);

        KRATOS_CATCH("")
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
    virtual void CalculateSensitivityGradient(const Element& rAdjointElem,
                                              const Variable<double>& rVariable,
                                              const Matrix& rDerivativesMatrix,
                                              Vector& rResponseGradient,
                                              ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;

        KRATOS_ERROR << "This should be implemented in the derived class." << std::endl;

        KRATOS_CATCH("");
    }

    // virtual void CalculateSensitivityGradient(const Condition& rAdjointCondition,
    //                                           const Variable<double>& rVariable,
    //                                           const Matrix& rDerivativesMatrix,
    //                                           Vector& rResponseGradient,
    //                                         mIsOutputFileOpened  ProcessInfo& rProcessInfo)
    // {
    //     KRATOS_TRY;

    //     KRATOS_ERROR << "This should be implemented in the derived class." << std::endl;

    //     KRATOS_CATCH("");
    // }

    /// Calculate the local gradient of response function w.r.t. the sensitivity variable.
    /**
     * @param[in]     rAdjointElem       the adjoint element.
     * @param[in]     rVariable          the sensitivity variable.
     * @param[in]     rDerivativesMatrix the transposed gradient of the element's
     *                                   residual w.r.t. the sensitivity variable.
     * @param[out]    rResponseGradient  the gradient of the response function.
     * @param[in,out] rProcessInfo       the current process info.
     */
    virtual void CalculateSensitivityGradient(const Element& rAdjointElem,
                                              const Variable<array_1d<double,3>>& rVariable,
                                              const Matrix& rDerivativesMatrix,
                                              Vector& rResponseGradient,
                                              ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;

        KRATOS_ERROR << "This should be implemented in the derived class." << std::endl;

        KRATOS_CATCH("");
    }

    // virtual void CalculateSensitivityGradient(const Condition& rAdjointCondition,
    //                                           const Variable<array_1d<double,3>>& rVariable,
    //                                           const Matrix& rDerivativesMatrix,
    //                                           Vector& rResponseGradient,
    //                                           ProcessInfo& rProcessInfo)
    // {
    //     KRATOS_TRY;

    //     KRATOS_ERROR << "This should be implemented in the derived class." << std::endl;

    //     KRATOS_CATCH("");
    // }

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

    ///@}

private:
    ///@name Member Variables
    ///@{

    std::string mSensitivityModelPartName;
    std::vector<std::string> mNodalSensitivityVariables;
    std::string mOutputFilename;
    std::ofstream mOutputFileStream; 
    std::string mResponseName;
    bool mIsOutputFileOpened;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    ///@}
};

///@} // Kratos Classes

///@} // Adjoint Fluid Application group

} /* namespace Kratos.*/

#endif /* KRATOS_RESPONSE_FUNCTION defined */

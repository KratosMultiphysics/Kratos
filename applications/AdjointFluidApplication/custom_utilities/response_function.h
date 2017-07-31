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
        mBoundaryModelPartName = rParameters["boundary_model_part_name"].GetString();
        Parameters nodal_sensitivity_variables = rParameters["nodal_sensitivity_variables"];
        mNodalSensitivityVariables.resize(nodal_sensitivity_variables.size());
        for (unsigned int i = 0; i < nodal_sensitivity_variables.size(); ++i)
            mNodalSensitivityVariables[i] = nodal_sensitivity_variables[i].GetString();
    }

    /// Destructor.
    virtual ~ResponseFunction()
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
        ModelPart& r_model_part = this->GetModelPart();

        if (r_model_part.HasSubModelPart(mBoundaryModelPartName) == false)
            KRATOS_ERROR << "No sub model part \"" << mBoundaryModelPartName
                         << "\"" << std::endl;
        
        for (auto &node : r_model_part.Nodes())
            node.Set(BOUNDARY, false);

        for (auto &node : r_model_part.GetSubModelPart(mBoundaryModelPartName).Nodes())
            node.Set(BOUNDARY, true);

        for (auto label : mNodalSensitivityVariables)
        {
            if (KratosComponents<Variable<double>>::Has(label) == true)
            {
                const Variable<double>& rVariable =
                    KratosComponents<Variable<double>>::Get(label);
                for (auto& node : r_model_part.Nodes())
                    node.FastGetSolutionStepValue(rVariable) = rVariable.Zero();
            }
            else if (KratosComponents<Variable<array_1d<double,3>>>::Has(label) == true)
            {
                const Variable<array_1d<double,3>>& rVariable =
                    KratosComponents<Variable<array_1d<double,3>>>::Get(label);
                for (auto& node : r_model_part.Nodes())
                    node.FastGetSolutionStepValue(rVariable) = rVariable.Zero();
            }
            else
                KRATOS_ERROR << "Unsupported variable: " << label << "." << std::endl;
        }
    }

    virtual void InitializeSolutionStep()
    {
    }

    virtual void FinalizeSolutionStep()
    {
      this->UpdateSensitivities();
    }

    virtual void Check()
    {
    }

    /// Calculate the local gradient w.r.t. velocity
    /**
     * @param[in]     rElem            the local adjoint element.
     * @param[in]     rAdjointMatrix   the transposed gradient of the local
     *                                 element's residual w.r.t. velocity.
     * @param[out]    rRHSContribution the gradient of the response function.
     * @param[in,out] rProcessInfo     the current process info.
     */
    virtual void CalculateFirstDerivativesGradient(const Element& rElem,
                                                   const Matrix& rAdjointMatrix,
                                                   Vector& rRHSContribution,
                                                   ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY

        if (rRHSContribution.size() != rAdjointMatrix.size1())
            rRHSContribution.resize(rAdjointMatrix.size1(), false);

        for (unsigned int k = 0; k < rRHSContribution.size(); ++k)
            rRHSContribution[k] = 0.0;

        KRATOS_CATCH("")
    }

    /// Calculate the local gradient w.r.t. acceleration
    /**
     * @param[in]     rElem            the local adjoint element.
     * @param[in]     rAdjointMatrix   the transposed gradient of the local
     *                                 element's residual w.r.t. acceleration.
     * @param[out]    rRHSContribution the gradient of the response function.
     * @param[in,out] rProcessInfo     the current process info.
     */
    virtual void CalculateSecondDerivativesGradient(const Element& rElem,
                                                    const Matrix& rAdjointMatrix,
                                                    Vector& rRHSContribution,
                                                    ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY

        if (rRHSContribution.size() != rAdjointMatrix.size1())
            rRHSContribution.resize(rAdjointMatrix.size1(), false);

        for (unsigned int k = 0; k < rRHSContribution.size(); ++k)
            rRHSContribution[k] = 0.0;

        KRATOS_CATCH("")
    }


    virtual void UpdateSensitivities()
    {
      KRATOS_TRY;

        for (auto label : mNodalSensitivityVariables)
        {
            if (KratosComponents<Variable<double>>::Has(label) == true)
            {
                const Variable<double>& rVariable =
                    KratosComponents<Variable<double>>::Get(label);
                this->UpdateNodalSensitivities(rVariable);
            }
            else if (KratosComponents<Variable<array_1d<double,3>>>::Has(label) == true)
            {
                const Variable<array_1d<double,3>>& rVariable =
                    KratosComponents<Variable<array_1d<double,3>>>::Get(label);
                this->UpdateNodalSensitivities(rVariable);
            }
            else
                KRATOS_ERROR << "Unsupported variable: " << label << "." << std::endl;
        }
      
      KRATOS_CATCH("");
    }


    /// Calculate the scalar valued response function
    virtual double CalculateValue(ModelPart& rModelPart)
    {
        KRATOS_TRY

        return 0.0;

        KRATOS_CATCH("")
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
            for (auto it = r_model_part.NodesBegin(); it != r_model_part.NodesEnd(); ++it)
                if (it->FastGetSolutionStepValue(PARTITION_INDEX) != r_comm.MyPID())
                    it->FastGetSolutionStepValue(rSensitivityVariable) =
                        rSensitivityVariable.Zero();
        }

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
                bool is_boundary = false;
                for (unsigned int i_node = 0; i_node < r_geom.PointsNumber(); ++i_node)
                    if (r_geom[i_node].Is(BOUNDARY) == true)
                    {
                        is_boundary = true;
                        break;
                    }

                if (is_boundary == false) // true for most elements
                    continue;

                it->CalculateSensitivityMatrix(
                    rSensitivityVariable, sensitivity_matrix[k], r_process_info);

                this->CalculateSensitivityContribution(
                    *it, rSensitivityVariable, sensitivity_matrix[k],
                    response_gradient[k], r_process_info);

                it->GetValuesVector(adjoint_vector[k]);

                if (sensitivity_vector[k].size() != sensitivity_matrix[k].size1())
                    sensitivity_vector[k].resize(sensitivity_matrix[k].size1(), false);

                noalias(sensitivity_vector[k]) =
                    delta_time * (prod(sensitivity_matrix[k], adjoint_vector[k]) +
                                  response_gradient[k]);

                this->AssembleNodalSensitivityContribution(
                    rSensitivityVariable, sensitivity_vector[k], r_geom);
            }
        }

        r_model_part.GetCommunicator().AssembleCurrentData(rSensitivityVariable);

        KRATOS_CATCH("")
    }

    /// Calculate the local gradient w.r.t. design variable.
    /**
     * @param[in]     rElem              the local adjoint element.
     * @param[in]     rVariable          the sensitivity variable.
     * @param[in]     rDerivativesMatrix the transposed gradient of the local
     *                                   element's residual.
     * @param[out]    rRHSContribution   the gradient of the response function.
     * @param[in,out] rProcessInfo       the current process info.
     */
    virtual void CalculateSensitivityContribution(const Element& rElem,
                                                  const Variable<double>& rVariable,
                                                  const Matrix& rDerivativesMatrix,
                                                  Vector& rRHSContribution,
                                                  ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY

        if (rRHSContribution.size() != rDerivativesMatrix.size1())
            rRHSContribution.resize(rDerivativesMatrix.size1(), false);

        for (unsigned int k = 0; k < rRHSContribution.size(); ++k)
            rRHSContribution[k] = 0.0;

        KRATOS_CATCH("")
    }

    /// Calculate the local gradient w.r.t. design variable.
    /**
     * @param[in]     rElem              the local adjoint element.
     * @param[in]     rVariable          the sensitivity variable.
     * @param[in]     rDerivativesMatrix the transposed gradient of the local
     *                                   element's residual.
     * @param[out]    rRHSContribution   the gradient of the response function.
     * @param[in,out] rProcessInfo       the current process info.
     */
    virtual void CalculateSensitivityContribution(const Element& rElem,
                                                  const Variable<array_1d<double,3>>& rVariable,
                                                  const Matrix& rDerivativesMatrix,
                                                  Vector& rRHSContribution,
                                                  ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY

        if (rRHSContribution.size() != rDerivativesMatrix.size1())
            rRHSContribution.resize(rDerivativesMatrix.size1(), false);

        for (unsigned int k = 0; k < rRHSContribution.size(); ++k)
            rRHSContribution[k] = 0.0;

        KRATOS_CATCH("")
    }

    void AssembleNodalSensitivityContribution(Variable<double> const& rSensitivityVariable,
                                              Vector const& rSensitivityVector,
                                              Element::GeometryType& rGeom)
    {
        unsigned int index = 0;
        for (unsigned int i_node = 0; i_node < rGeom.PointsNumber(); ++i_node)
        {
            if (rGeom[i_node].Is(BOUNDARY) == true)
            {
                double& rSensitivity =
                    rGeom[i_node].FastGetSolutionStepValue(rSensitivityVariable);
                rGeom[i_node].SetLock();
                rSensitivity += rSensitivityVector[index++];
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
            if (rGeom[i_node].Is(BOUNDARY) == true)
            {
                array_1d<double, 3>& rSensitivity =
                    rGeom[i_node].FastGetSolutionStepValue(rSensitivityVariable);
                rGeom[i_node].SetLock();
                for (unsigned int d = 0; d < rGeom.WorkingSpaceDimension(); ++d)
                    rSensitivity[d] += rSensitivityVector[index++];
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

    std::string mBoundaryModelPartName;
    std::vector<std::string> mNodalSensitivityVariables;

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

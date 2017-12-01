//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Michael Andre, https://github.com/msandre
//

#if !defined(KRATOS_RESPONSE_FUNCTION_SENSITIVITY_BUILDER_UTILITY_H_INCLUDED)
#define KRATOS_RESPONSE_FUNCTION_SENSITIVITY_BUILDER_UTILITY_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"

namespace Kratos
{
///@name Kratos Classes
///@{

/// Implementation for assembling sensitivities from local contributions.
/**
 * The utility is intended to be inherited by a concrete response function,
 * which customizes CalculatePartialSensitivity and calls the provided
 * build functions to assemble its sensitivity variables.
 * 
 */
class ResponseFunctionSensitivityBuilderUtility
{
protected:
    ///@name Type Definitions
    ///@{

    typedef Geometry<Node<3>> GeometryType;

    typedef ModelPart::NodesContainerType NodesContainerType;

    typedef ModelPart::ElementsContainerType ElementsContainerType;

    typedef ModelPart::ConditionsContainerType ConditionsContainerType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Destructor.
    ~ResponseFunctionSensitivityBuilderUtility() = default;

    ///@}
    ///@name Protected Operations
    ///@{

    /// Calculate the local gradient w.r.t. the sensitivity variable.
    /**
     * Defines the local contribution from an element for a sensitivity variable.
     * This is overridden by the derived response function and called in the
     * assembly loop.
     * 
     * @param[in]     rVariable             sensitivity variable.
     * @param[in]     rElement              local adjoint element.
     * @param[in]     rSensitivityMatrix    transposed gradient of the residual
     *                                      w.r.t. the sensitivity variable.
     * @param[out]    rPartialSensitivity   gradient of the response function
     *                                      w.r.t. the sensitivity variable.
     * @param[in]     rProcessInfo          current process info.
     */
    virtual void CalculatePartialSensitivity(Variable<double> const& rVariable,
                                             Element const& rElement,
                                             Matrix const& rSensitivityMatrix,
                                             Vector& rPartialSensitivity,
                                             ProcessInfo const& rProcessInfo) const
    {
        KRATOS_TRY;

        KRATOS_ERROR << "Calling the base class method." << std::endl;

        KRATOS_CATCH("");
    }

    virtual void CalculatePartialSensitivity(Variable<double> const& rVariable,
                                             Condition const& rCondition,
                                             Matrix const& rSensitivityMatrix,
                                             Vector& rPartialSensitivity,
                                             ProcessInfo const& rProcessInfo) const
    {
        KRATOS_TRY;

        KRATOS_ERROR << "Calling the base class method." << std::endl;

        KRATOS_CATCH("");
    }

    virtual void CalculatePartialSensitivity(Variable<array_1d<double, 3>> const& rVariable,
                                             Element const& rElement,
                                             Matrix const& rSensitivityMatrix,
                                             Vector& rPartialSensitivity,
                                             ProcessInfo const& rProcessInfo) const
    {
        KRATOS_TRY;

        KRATOS_ERROR << "Calling the base class method." << std::endl;

        KRATOS_CATCH("");
    }

    virtual void CalculatePartialSensitivity(Variable<array_1d<double, 3>> const& rVariable,
                                             Condition const& rCondition,
                                             Matrix const& rSensitivityMatrix,
                                             Vector& rPartialSensitivity,
                                             ProcessInfo const& rProcessInfo) const
    {
        KRATOS_TRY;

        KRATOS_ERROR << "Calling the base class method." << std::endl;

        KRATOS_CATCH("");
    }

    void BuildNodalSolutionStepSensitivities(std::string const& rVariable,
                                             ModelPart& rModelPart,
                                             double ScalingFactor = 1.0)
    {
        KRATOS_TRY;

        if (KratosComponents<Variable<double>>::Has(rVariable) == true)
        {
            const Variable<double>& r_variable =
                KratosComponents<Variable<double>>::Get(rVariable);
            BuildNodalSolutionStepSensitivities(r_variable, rModelPart, ScalingFactor);
        }
        else if (KratosComponents<Variable<array_1d<double, 3>>>::Has(rVariable) == true)
        {
            const Variable<array_1d<double, 3>>& r_variable =
                KratosComponents<Variable<array_1d<double, 3>>>::Get(rVariable);
            BuildNodalSolutionStepSensitivities(r_variable, rModelPart, ScalingFactor);
        }
        else
            KRATOS_ERROR << "Unsupported variable: " << rVariable << "." << std::endl;

        KRATOS_CATCH("");
    }

    template <class TDataType>
    void BuildNodalSolutionStepSensitivities(Variable<TDataType> const& rVariable,
                                             ModelPart& rModelPart,
                                             double ScalingFactor = 1.0)
    {
        KRATOS_TRY;

        Communicator& r_comm = rModelPart.GetCommunicator();
        if (r_comm.TotalProcesses() > 1)
        {
// Make sure we only add the old sensitivity once when we assemble.
#pragma omp parallel
            {
                ModelPart::NodeIterator nodes_begin;
                ModelPart::NodeIterator nodes_end;
                OpenMPUtils::PartitionedIterators(rModelPart.Nodes(), nodes_begin, nodes_end);
                for (auto it = nodes_begin; it != nodes_end; ++it)
                    if (it->FastGetSolutionStepValue(PARTITION_INDEX) != r_comm.MyPID())
                        it->FastGetSolutionStepValue(rVariable) = rVariable.Zero();
            }
        }

        BuildNodalSolutionStepElementContributions(
            rVariable, ScalingFactor, rModelPart.Elements(), rModelPart.GetProcessInfo());

        BuildNodalSolutionStepConditionContributions(rVariable, ScalingFactor,
                                                     rModelPart.Conditions(),
                                                     rModelPart.GetProcessInfo());

        r_comm.AssembleCurrentData(rVariable);

        KRATOS_CATCH("");
    }

    ///@}

private:
    ///@name Private Operations
    ///@{

    template <class TDataType>
    void BuildNodalSolutionStepElementContributions(Variable<TDataType> const& rVariable,
                                                    double ScalingFactor,
                                                    ElementsContainerType& rElements,
                                                    ProcessInfo const& rProcessInfo)
    {
        KRATOS_TRY;

        const int num_threads = OpenMPUtils::GetNumThreads();
        std::vector<Vector> local_sensitivity(num_threads);
        std::vector<Vector> partial_sensitivity(num_threads);
        std::vector<Vector> adjoint_vector(num_threads);
        std::vector<Matrix> sensitivity_matrix(num_threads);

#pragma omp parallel
        {
            ModelPart::ElementIterator elements_begin;
            ModelPart::ElementIterator elements_end;
            OpenMPUtils::PartitionedIterators(rElements, elements_begin, elements_end);
            int k = OpenMPUtils::ThisThread();

            for (auto it = elements_begin; it != elements_end; ++it)
            {
                Element::GeometryType& r_geom = it->GetGeometry();

                if (HasActiveNodes(r_geom) == false)
                    continue;

                it->CalculateSensitivityMatrix(rVariable, sensitivity_matrix[k], rProcessInfo);

                CalculatePartialSensitivity(rVariable, *it, sensitivity_matrix[k],
                                            partial_sensitivity[k], rProcessInfo);

                it->GetValuesVector(adjoint_vector[k]);

                if (local_sensitivity[k].size() != sensitivity_matrix[k].size1())
                    local_sensitivity[k].resize(sensitivity_matrix[k].size1(), false);

                noalias(local_sensitivity[k]) =
                    ScalingFactor * (prod(sensitivity_matrix[k], adjoint_vector[k]) +
                                     partial_sensitivity[k]);

                AssembleNodalSolutionStepSensitivityContribution(
                    rVariable, local_sensitivity[k], r_geom);
            }
        }

        KRATOS_CATCH("");
    }

    template <class TDataType>
    void BuildNodalSolutionStepConditionContributions(Variable<TDataType> const& rVariable,
                                                      double ScalingFactor,
                                                      ConditionsContainerType& rConditions,
                                                      ProcessInfo const& rProcessInfo)
    {
        KRATOS_TRY;

        const int num_threads = OpenMPUtils::GetNumThreads();
        std::vector<Vector> local_sensitivity(num_threads);
        std::vector<Vector> partial_sensitivity(num_threads);
        std::vector<Vector> adjoint_vector(num_threads);
        std::vector<Matrix> sensitivity_matrix(num_threads);

#pragma omp parallel
        {
            ModelPart::ConditionIterator conditions_begin;
            ModelPart::ConditionIterator conditions_end;
            OpenMPUtils::PartitionedIterators(rConditions, conditions_begin, conditions_end);
            int k = OpenMPUtils::ThisThread();

            for (auto it = conditions_begin; it != conditions_end; ++it)
            {
                Condition::GeometryType& r_geom = it->GetGeometry();

                if (HasActiveNodes(r_geom) == false)
                    continue;

                it->CalculateSensitivityMatrix(rVariable, sensitivity_matrix[k], rProcessInfo);

                CalculatePartialSensitivity(rVariable, *it, sensitivity_matrix[k],
                                            partial_sensitivity[k], rProcessInfo);

                it->GetValuesVector(adjoint_vector[k]);

                if (local_sensitivity[k].size() != sensitivity_matrix[k].size1())
                    local_sensitivity[k].resize(sensitivity_matrix[k].size1(), false);

                noalias(local_sensitivity[k]) =
                    ScalingFactor * (prod(sensitivity_matrix[k], adjoint_vector[k]) +
                                     partial_sensitivity[k]);

                AssembleNodalSolutionStepSensitivityContribution(
                    rVariable, local_sensitivity[k], r_geom);
            }
        }

        KRATOS_CATCH("");
    }

    bool HasActiveNodes(GeometryType const& rGeom)
    {
        bool result = false;
        for (unsigned i_node = 0; i_node < rGeom.PointsNumber(); ++i_node)
            if (rGeom[i_node].GetValue(UPDATE_SENSITIVITIES))
            {
                result = true;
                break;
            }
        return result;
    }

    void AssembleNodalSolutionStepSensitivityContribution(Variable<double> const& rVariable,
                                                          Vector const& rSensitivityVector,
                                                          GeometryType& rGeom)
    {
        unsigned int index = 0;
        for (unsigned int i_node = 0; i_node < rGeom.PointsNumber(); ++i_node)
        {
            if (rGeom[i_node].GetValue(UPDATE_SENSITIVITIES))
            {
                double& r_sensitivity = rGeom[i_node].FastGetSolutionStepValue(rVariable);
                rGeom[i_node].SetLock();
                r_sensitivity += rSensitivityVector[index++];
                rGeom[i_node].UnSetLock();
            }
            else
                ++index;
        }
    }

    void AssembleNodalSolutionStepSensitivityContribution(Variable<array_1d<double, 3>> const& rVariable,
                                                          Vector const& rSensitivityVector,
                                                          GeometryType& rGeom)
    {
        unsigned int index = 0;
        for (unsigned int i_node = 0; i_node < rGeom.PointsNumber(); ++i_node)
        {
            if (rGeom[i_node].GetValue(UPDATE_SENSITIVITIES))
            {
                array_1d<double, 3>& r_sensitivity =
                    rGeom[i_node].FastGetSolutionStepValue(rVariable);
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

};

///@} // Kratos Classes
} /* namespace Kratos.*/

#endif /* KRATOS_RESPONSE_FUNCTION_SENSITIVITY_BUILDER_UTILITY_H_INCLUDED defined */

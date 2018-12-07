//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    
//

#if !defined(KRATOS_ELEMENT_SENSITIVITY_BUILDER_H_INCLUDED)
#define KRATOS_ELEMENT_SENSITIVITY_BUILDER_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "utilities/openmp_utils.h"
#include "response_functions/adjoint_response_function.h"

namespace Kratos
{
///@name Kratos Classes
///@{

class ElementSensitivityBuilder
{
public:
    ///@name Type Definitions
    ///@{

    typedef Geometry<Node<3>> GeometryType;

    ///@}
    ///@name Life Cycle
    ///@{

    ElementSensitivityBuilder(ModelPart& rModelPart, AdjointResponseFunction::Pointer pResponseFunction)
        : mrModelPart(rModelPart), mpResponseFunction(pResponseFunction)
    {
    }

    ///@}
    ///@name Operations
    ///@{

    void BuildElementSensitivities(std::string const& rVariable, double ScalingFactor = 1.0)
    {
        KRATOS_TRY;

        if (KratosComponents<Variable<double>>::Has(rVariable) == true)
        {
            const Variable<double>& r_variable =
                KratosComponents<Variable<double>>::Get(rVariable);
            BuildElementSensitivities(r_variable, ScalingFactor);
        }
        else if (KratosComponents<Variable<array_1d<double, 3>>>::Has(rVariable) == true)
        {
            const Variable<array_1d<double, 3>>& r_variable =
                KratosComponents<Variable<array_1d<double, 3>>>::Get(rVariable);
            BuildElementSensitivities(r_variable, ScalingFactor);
        }
        else
            KRATOS_ERROR << "Unsupported variable: " << rVariable << "." << std::endl;

        KRATOS_CATCH("");
    }

    ///@}

private:

    ///@name Private Operations
    ///@{

    template <class TDataType>
    void BuildElementSensitivities(Variable<TDataType> const& rVariable, double ScalingFactor)
    {
        KRATOS_TRY;

        auto& r_elements = mrModelPart.Elements();
        const auto& r_process_info = mrModelPart.GetProcessInfo();

        const int num_threads = OpenMPUtils::GetNumThreads();
        std::vector<Vector> local_sensitivity(num_threads);
        std::vector<Vector> partial_sensitivity(num_threads);
        std::vector<Vector> adjoint_vector(num_threads);
        std::vector<Matrix> sensitivity_matrix(num_threads);

#pragma omp parallel
        {
            ModelPart::ElementIterator elements_begin;
            ModelPart::ElementIterator elements_end;
            OpenMPUtils::PartitionedIterators(r_elements, elements_begin, elements_end);
            int k = OpenMPUtils::ThisThread();

            for (auto it = elements_begin; it != elements_end; ++it)
            {
                it->CalculateSensitivityMatrix(rVariable, sensitivity_matrix[k], r_process_info);

                mpResponseFunction->CalculatePartialSensitivity(
                    *it, rVariable, sensitivity_matrix[k],
                    partial_sensitivity[k], r_process_info);

                it->GetValuesVector(adjoint_vector[k]);

                if (local_sensitivity[k].size() != sensitivity_matrix[k].size1())
                    local_sensitivity[k].resize(sensitivity_matrix[k].size1(), false);

                noalias(local_sensitivity[k]) =
                    ScalingFactor * (prod(sensitivity_matrix[k], adjoint_vector[k]) +
                                     partial_sensitivity[k]);

                AssembleElementSensitivityContribution(
                    rVariable, local_sensitivity[k], *it);
            }
        }

        KRATOS_CATCH("");
    }

    ///@}

private:
    ///@name Member Variables
    ///@{

    ModelPart& mrModelPart;
    AdjointResponseFunction::Pointer mpResponseFunction;

    ///@}
    ///@name Private Operations
    ///@{

    void AssembleElementSensitivityContribution(Variable<double> const& rVariable,
                                                Vector const& rSensitivityVector,
                                                Element& rElement)
    {
        KRATOS_DEBUG_ERROR_IF(rSensitivityVector.size() != 1) << "rSensitivityVector.size() = " << rSensitivityVector.size() << std::endl;
        rElement.GetValue(rVariable) += rSensitivityVector[0];
    }

    void AssembleElementSensitivityContribution(Variable<array_1d<double, 3>> const& rVariable,
                                                Vector const& rSensitivityVector,
                                                Element& rElement)
    {
        array_1d<double, 3>& r_sensitivity = rElement.GetValue(rVariable);
        const auto ws_dim = rElement.GetGeometry().WorkingSpaceDimension();
        KRATOS_DEBUG_ERROR_IF(rSensitivityVector.size() != ws_dim) << "rSensitivityVector.size() = " << rSensitivityVector.size() << std::endl;
        for (unsigned d = 0; d < ws_dim; ++d)
            r_sensitivity[d] += rSensitivityVector[d];
    }

    ///@}

};

///@} // Kratos Classes
} /* namespace Kratos.*/

#endif /* KRATOS_ELEMENT_SENSITIVITY_BUILDER_H_INCLUDED defined */

// KRATOS ___ ___  _  ___   __   ___ ___ ___ ___
//       / __/ _ \| \| \ \ / /__|   \_ _| __| __|
//      | (_| (_) | .` |\ V /___| |) | || _|| _|
//       \___\___/|_|\_| \_/    |___/___|_| |_|  APPLICATION
//
//  License:       BSD License
//                 Kratos default license: kratos/license.txt
//
//  Main authors:  Jordi Cotela
//

#ifndef KRATOS_LOCAL_TEMPERATURE_AVERAGE_RESPONSE_FUNCTION_H_INCLUDED
#define KRATOS_LOCAL_TEMPERATURE_AVERAGE_RESPONSE_FUNCTION_H_INCLUDED

#include "includes/kratos_flags.h"
#include "includes/model_part.h"
#include "utilities/variable_utils.h"
#include "response_functions/adjoint_response_function.h"

namespace Kratos {
///@addtogroup ConvectionDiffusionApplication
///@{

///@name Kratos Classes
///@{

class LocalTemperatureAverageResponseFunction: public AdjointResponseFunction
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(LocalTemperatureAverageResponseFunction);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    LocalTemperatureAverageResponseFunction(Parameters Settings, ModelPart& rModelPart)
    {
        KRATOS_TRY;

        mTargetModelPartName = Settings["model_part_name"].GetString();

        auto& r_target_model_part = GetTargetModelPart(rModelPart, mTargetModelPartName);
        auto& r_nodes = r_target_model_part.Nodes();
        mNumNodes = r_nodes.size();

        VariableUtils variable_utils;
        variable_utils.SetFlag(STRUCTURE,true,r_nodes);

        // Note: this should not be parallel, the operation is not threadsafe if the variable is uninitialized
        for (auto& r_node : r_nodes)
        {
            r_node.SetValue(NUMBER_OF_NEIGHBOUR_ELEMENTS,0);
        }

        mNumNodes = rModelPart.GetCommunicator().GetDataCommunicator().SumAll(mNumNodes);

        auto& r_elements = rModelPart.Elements();
        const int num_elements = r_elements.size();

        #pragma omp parallel for
        for (int i = 0; i < num_elements; i++)
        {
            auto i_elem = r_elements.begin() + i;
            auto& r_geom = i_elem->GetGeometry();
            for (unsigned int i = 0; i < r_geom.PointsNumber(); i++)
            {
                auto& r_node = r_geom[i];
                if (r_node.Is(STRUCTURE))
                {
                    r_node.SetLock();
                    r_node.GetValue(NUMBER_OF_NEIGHBOUR_ELEMENTS) += 1;
                    r_node.UnSetLock();
                }
            }
        }

        rModelPart.GetCommunicator().AssembleNonHistoricalData(NUMBER_OF_NEIGHBOUR_ELEMENTS);

        KRATOS_CATCH("");
    }

    /// Destructor.
    ~LocalTemperatureAverageResponseFunction() override
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

        KRATOS_CATCH("");
    }

    void CalculateGradient(const Element& rAdjointElement,
                           const Matrix& rResidualGradient,
                           Vector& rResponseGradient,
                           const ProcessInfo& rProcessInfo) override
    {
        ComputePointTemperatureSensitivityContribution(rResidualGradient, rAdjointElement.GetGeometry().Points(),rResponseGradient);
    }

    void CalculateGradient(const Condition& rAdjointCondition,
                           const Matrix& rResidualGradient,
                           Vector& rResponseGradient,
                           const ProcessInfo& rProcessInfo) override
    {
        noalias(rResponseGradient) = ZeroVector(rResidualGradient.size1());
    }

    void CalculateFirstDerivativesGradient(const Element& rAdjointElement,
                                           const Matrix& rResidualGradient,
                                           Vector& rResponseGradient,
                                           const ProcessInfo& rProcessInfo) override
    {
        ComputePointTemperatureSensitivityContribution(rResidualGradient, rAdjointElement.GetGeometry().Points(),rResponseGradient);
    }

    void CalculateSecondDerivativesGradient(const Element& rAdjointElement,
                                            const Matrix& rResidualGradient,
                                            Vector& rResponseGradient,
                                            const ProcessInfo& rProcessInfo) override
    {
        ComputePointTemperatureSensitivityContribution(rResidualGradient, rAdjointElement.GetGeometry().Points(),rResponseGradient);
    }

    void CalculatePartialSensitivity(Element& rAdjointElement,
                                     const Variable<array_1d<double, 3>>& rVariable,
                                     const Matrix& rSensitivityMatrix,
                                     Vector& rSensitivityGradient,
                                     const ProcessInfo& rProcessInfo) override
    {
        if (rSensitivityGradient.size() != rSensitivityMatrix.size1())
            rSensitivityGradient.resize(rSensitivityMatrix.size1(), false);
        noalias(rSensitivityGradient) = ZeroVector(rSensitivityMatrix.size1());
    }

    void CalculatePartialSensitivity(Condition& rAdjointElement,
                                     const Variable<array_1d<double, 3>>& rVariable,
                                     const Matrix& rSensitivityMatrix,
                                     Vector& rSensitivityGradient,
                                     const ProcessInfo& rProcessInfo) override
    {
        if (rSensitivityGradient.size() != rSensitivityMatrix.size1())
            rSensitivityGradient.resize(rSensitivityMatrix.size1(), false);
        noalias(rSensitivityGradient) = ZeroVector(rSensitivityMatrix.size1());
    }

    double CalculateValue(ModelPart& rModelPart) override
    {
        KRATOS_TRY;
        const ModelPart& r_target_model_part =
            GetTargetModelPart(rModelPart, mTargetModelPartName);

        const double domain_aggregated_temperature =
            VariableUtils().SumHistoricalVariable<double>(TEMPERATURE, r_target_model_part);

        const Communicator& r_communicator = r_target_model_part.GetCommunicator();
        const int number_of_nodes = r_communicator.LocalMesh().NumberOfNodes();
        const int total_nodes = r_communicator.GetDataCommunicator().SumAll(number_of_nodes);

        return domain_aggregated_temperature / static_cast<double>(total_nodes);

        KRATOS_CATCH("");
    }

    ///@}

protected:
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    ///@}

private:
    ///@name Member Variables
    ///@{

	int mNumNodes = 0;
	std::string mTargetModelPartName;

    ///@}
    ///@name Private Operators
    ///@{

    void ComputePointTemperatureSensitivityContribution(
        const Matrix& rDerivativesOfResidual,
        const Element::NodesArrayType& rNodes,
        Vector& rLocalSensitivityContribution) const
    {
        if (rLocalSensitivityContribution.size() != rDerivativesOfResidual.size1())
            rLocalSensitivityContribution.resize(rDerivativesOfResidual.size1(), false);

        noalias(rLocalSensitivityContribution) = ZeroVector(rLocalSensitivityContribution.size());

        const unsigned int num_nodes = rNodes.size();

        for (unsigned int i = 0; i < num_nodes; i++)
        {
            if (rNodes[i].Is(STRUCTURE))
            {
                double factor = 1.0 / (rNodes[i].GetValue(NUMBER_OF_NEIGHBOUR_ELEMENTS)*mNumNodes);
                rLocalSensitivityContribution[i] = factor;
            }
        }
    }

    ModelPart& GetTargetModelPart(ModelPart& rModelPart, const std::string& rTargetModelPartName)
    {
        KRATOS_TRY;
        if (rModelPart.Name() == rTargetModelPartName)
		{
            return rModelPart;
		}
		else if (rModelPart.HasSubModelPart(rTargetModelPartName))
		{
			return rModelPart.GetSubModelPart(rTargetModelPartName);
		}
		else
		{
			KRATOS_ERROR << "Unknown ModelPart " << rTargetModelPartName << "." << std::endl;
		}
        KRATOS_CATCH("")
        return rModelPart;
    }

    ///@}
    ///@name Private Operations
    ///@{

    ///@}
};

///@} // Kratos Classes

///@} // ConvectionDiffusionApplication group

}

#endif // KRATOS_LOCAL_TEMPERATURE_AVERAGE_RESPONSE_FUNCTION_H_INCLUDED

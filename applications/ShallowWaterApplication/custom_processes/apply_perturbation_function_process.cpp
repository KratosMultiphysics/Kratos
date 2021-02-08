//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Maso Sotomayor
//


// System includes


// External includes


// Project includes
#include "includes/checks.h"
#include "apply_perturbation_function_process.h"


namespace Kratos
{

template<class TVarType>
ApplyPerturbationFunctionProcess<TVarType>::ApplyPerturbationFunctionProcess(
    ModelPart& rThisModelPart,
    NodePointerType pNode,
    TVarType& rThisVariable,
    Parameters& rThisParameters
) : Process(), mrModelPart(rThisModelPart), mrVariable(rThisVariable)
{
    ValidateParameters(rThisParameters);
    mSourcePoints.push_back(pNode);
}


template<class TVarType>
ApplyPerturbationFunctionProcess<TVarType>::ApplyPerturbationFunctionProcess(
    ModelPart& rThisModelPart,
    NodesArrayType& rSourcePoints,
    TVarType& rThisVariable,
    Parameters& rThisParameters
) : Process(), mrModelPart(rThisModelPart), mrVariable(rThisVariable)
{
    ValidateParameters(rThisParameters);
    mSourcePoints = rSourcePoints;
}


template<class TVarType>
int ApplyPerturbationFunctionProcess<TVarType>::Check()
{
    if (mrModelPart.Nodes().size() != 0) {
        const auto& r_node = *mrModelPart.NodesBegin();
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(mrVariable, r_node);
    }
    KRATOS_CHECK(mInfluenceDistance >= std::numeric_limits<double>::epsilon());
    return 0;
}


template<class TVarType>
void ApplyPerturbationFunctionProcess<TVarType>::Execute()
{
    ExecuteBeforeSolutionLoop();
}


template<class TVarType>
void ApplyPerturbationFunctionProcess<TVarType>::ExecuteBeforeSolutionLoop()
{
    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(mrModelPart.Nodes().size()); i++)
    {
        auto i_node = mrModelPart.NodesBegin() + i;
        double distance = ComputeDistance(*i_node.base());
        double& r_value = i_node->FastGetSolutionStepValue(mrVariable);
        r_value = ComputeInitialValue(distance);
    }
}


template<class TVarType>
double ApplyPerturbationFunctionProcess<TVarType>::ComputeDistance(NodePointerType pNode)
{
    array_1d<double, 3>& coord = pNode->Coordinates();
    double sqr_distance = std::pow(mInfluenceDistance, 2) + 1.0;
    for (IndexType i = 0; i < mSourcePoints.size(); i++)
    {
        auto search_node = mSourcePoints.begin() + i;
        array_1d<double, 3>& source_coord = search_node->Coordinates();
        sqr_distance = std::min(sqr_distance, PointPointSquareDistance(coord, source_coord));
    }
    return std::sqrt(sqr_distance);
}


template<class TVarType>
double ApplyPerturbationFunctionProcess<TVarType>::PointPointSquareDistance(array_1d<double, 3>& rCoordA, array_1d<double, 3>& rCoordB)
{
    return std::pow(rCoordA[0] - rCoordB[0], 2) + std::pow(rCoordA[1] - rCoordB[1], 2) + std::pow(rCoordA[2] - rCoordB[2], 2);
}


template<class TVarType>
double ApplyPerturbationFunctionProcess<TVarType>::ComputeInitialValue(double& rDistance)
{
    double result = mDefaultValue;
    if (rDistance < mInfluenceDistance)
        result += 0.5 * mPerturbation * (1 + std::cos(mHalfWaveNumber * rDistance));
    return result;
}


template<class TVarType>
void ApplyPerturbationFunctionProcess<TVarType>::ValidateParameters(Parameters& rParameters)
{
    // default parameters
    Parameters default_parameters = Parameters(R"(
    {
        "default_value"              : 0.0,
        "distance_of_influence"      : 1.0,
        "maximum_perturbation_value" : 1.0
    })");
    rParameters.ValidateAndAssignDefaults(default_parameters);

    // Initialization of member variables
    mDefaultValue = rParameters["default_value"].GetDouble();
    mInfluenceDistance = rParameters["distance_of_influence"].GetDouble();
    mPerturbation = rParameters["maximum_perturbation_value"].GetDouble();
    mHalfWaveNumber = std::acos(-1) / mInfluenceDistance;
}

template class ApplyPerturbationFunctionProcess<Variable<double>>;

}  // namespace Kratos.

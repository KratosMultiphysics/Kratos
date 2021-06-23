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
#include "move_mesh_utility.h"
#include "utilities/parallel_utilities.h"


namespace Kratos
{

MoveMeshUtility::MoveMeshUtility(
    ModelPart& rLagrangianModelPart,
    ModelPart& rEulerianModelPart,
    Parameters ThisParameters)
    : mrLagrangianModelPart(rLagrangianModelPart)
    , mrEulerianModelPart(rEulerianModelPart)
    , mLagrangianSearchStructure(mrLagrangianModelPart)
    , mEulerianSearchStructure(mrEulerianModelPart)
{
    ThisParameters.ValidateAndAssignDefaults(GetDefaultParameters());
    mMaxResults = ThisParameters["maximum_results"].GetDouble();

    FillVariablesList(mScalarVariablesToLagrangian, ThisParameters["map_variables_to_lagrangian"]);
    FillVariablesList(mVectorVariablesToLagrangian, ThisParameters["map_variables_to_lagrangian"]);
    FillVariablesList(mScalarVariablesToEulerian, ThisParameters["map_variables_to_eulerian"]);
    FillVariablesList(mVectorVariablesToEulerian, ThisParameters["map_variables_to_eulerian"]);
}

const Parameters MoveMeshUtility::GetDefaultParameters() const
{
    const Parameters default_parameters = Parameters(R"({
        "map_variables_to_lagrangian" : ["TOPOGRAPHY","MANNING"],
        "map_variables_to_eulerian"   : ["HEIGHT","VELOCITY"],
        "maximum_results"             : 10000
    })" );
    return default_parameters;
}

int MoveMeshUtility::Check()
{
    KRATOS_ERROR_IF(mrEulerianModelPart.NumberOfNodes() == 0) << "Move mesh utility: The Eulerian model part is empty\n" << mrEulerianModelPart << std::endl;
    return 0;
}

void MoveMeshUtility::Initialize()
{
    mEulerianSearchStructure.UpdateSearchDatabase();
}

void MoveMeshUtility::MoveMesh()
{
    std::size_t num_nodes = mrEulerianModelPart.ElementsBegin()->GetGeometry().size();
    double dt = mrLagrangianModelPart.GetProcessInfo()[DELTA_TIME];

    struct TLS {
        Vector N;
        ResultContainerType results;
    };
    TLS tls;
    tls.N.resize(num_nodes);
    tls.results.resize(mMaxResults);

    block_for_each(mrLagrangianModelPart.Nodes(), tls, [&](NodeType& rNode, TLS& rTLS){
        Element::Pointer p_element;
        auto results_begin = rTLS.results.begin();
        bool is_found = MoveNode(rNode, dt, rTLS.N, p_element, results_begin);
        if (is_found)
            MapToLagrangian(rNode, rTLS.N, p_element);
    });
}

void MoveMeshUtility::MapResults()
{
    mLagrangianSearchStructure.UpdateSearchDatabase();

    std::size_t num_nodes;
    if (mrLagrangianModelPart.NumberOfNodes() != 0)
        num_nodes = mrLagrangianModelPart.ElementsBegin()->GetGeometry().size();

    struct TLS {
        Vector N;
        ResultContainerType results;
    };
    TLS tls;
    tls.N.resize(num_nodes);
    tls.results.resize(mMaxResults);

    block_for_each(mrEulerianModelPart.Nodes(), tls, [&](NodeType& rNode, TLS& rTLS){
        array_1d<double,3>& r_pos = rNode.Coordinates();
        Element::Pointer p_element;
        auto results_begin = rTLS.results.begin();
        bool is_found = mLagrangianSearchStructure.FindPointOnMesh(r_pos, rTLS.N, p_element, results_begin, mMaxResults);
        MapToEulerian(rNode, rTLS.N, p_element, is_found);
    });
}

bool MoveMeshUtility::MoveNode(
    NodeType& rNode,
    double Dt,
    Vector& rN,
    Element::Pointer& pElement,
    ResultIteratorType& rResultBegin)
{
    array_1d<double,3>& r_pos = rNode.Coordinates();
    const array_1d<double,3>& vel = rNode.FastGetSolutionStepValue(VELOCITY);
    const array_1d<double,3>& acc = rNode.FastGetSolutionStepValue(ACCELERATION);
    r_pos += vel * Dt + 0.5 * acc * Dt * Dt;
    rNode.FastGetSolutionStepValue(DISPLACEMENT) = r_pos;
    bool is_found = mEulerianSearchStructure.FindPointOnMesh(r_pos, rN, pElement, rResultBegin, mMaxResults);
    return is_found;
}

void MoveMeshUtility::MapToLagrangian(
    NodeType& rNode,
    const Vector& rN,
    const Element::Pointer pElement)
{
    GeometryType r_geom = pElement->GetGeometry();
    for (std::size_t v = 0; v != mScalarVariablesToLagrangian.size(); ++v)
    {
        const Variable<double>& r_var = *(mScalarVariablesToLagrangian[v]);
        InterpolateVariable(rNode, rN, r_geom, r_var);
    }
    for (std::size_t v = 0; v != mVectorVariablesToLagrangian.size(); ++v)
    {
        const Variable<array_1d<double,3>>& r_var = *(mVectorVariablesToLagrangian[v]);
        InterpolateVariable(rNode, rN, r_geom, r_var);
    }
}

void MoveMeshUtility::MapToEulerian(
    NodeType& rNode,
    const Vector& rN,
    const Element::Pointer pElement,
    const bool IsFound)
{
    if (IsFound) {
        GeometryType r_geom = pElement->GetGeometry();
        for (std::size_t v = 0; v != mScalarVariablesToEulerian.size(); ++v)
        {
            const Variable<double>& r_var = *(mScalarVariablesToEulerian[v]);
            InterpolateVariable(rNode, rN, r_geom, r_var);
        }
        for (std::size_t v = 0; v != mVectorVariablesToEulerian.size(); ++v)
        {
            const Variable<array_1d<double,3>>& r_var = *(mVectorVariablesToEulerian[v]);
            InterpolateVariable(rNode, rN, r_geom, r_var);
        }
    } else {
        for (std::size_t v = 0; v != mScalarVariablesToEulerian.size(); ++v)
        {
            const Variable<double>& r_var = *(mScalarVariablesToEulerian[v]);
            rNode.FastGetSolutionStepValue(r_var) = 0.0;
        }
        for (std::size_t v = 0; v != mVectorVariablesToEulerian.size(); ++v)
        {
            const Variable<array_1d<double,3>>& r_var = *(mVectorVariablesToEulerian[v]);
            rNode.FastGetSolutionStepValue(r_var) = ZeroVector(3);
        }
    }
}

template<class TDataType>
void MoveMeshUtility::FillVariablesList(
    std::vector<const Variable<TDataType>*>& rVariablesList,
    const Parameters VariablesNames)
{
    const auto variables_names = VariablesNames.GetStringArray();
    for (const std::string& variable_name : variables_names)
    {
        if (KratosComponents<Variable<TDataType>>::Has(variable_name))
        {
            const auto& r_variable = KratosComponents<Variable<TDataType>>::Get(variable_name);
            rVariablesList.push_back(&r_variable);
        }
    }
}

template<class TDataType>
void MoveMeshUtility::InterpolateVariable(
    NodeType& rNode,
    const Vector& rN,
    const GeometryType& rGeometry,
    const Variable<TDataType>& rVariable)
{
    TDataType& r_value = rNode.FastGetSolutionStepValue(rVariable);
    r_value = rN[0] * rGeometry[0].FastGetSolutionStepValue(rVariable);
    for (std::size_t i = 1; i != rGeometry.size(); ++i)
    {
        r_value += rN[i] * rGeometry[i].FastGetSolutionStepValue(rVariable);
    }
}

template void MoveMeshUtility::FillVariablesList(std::vector<const Variable<double>*>&, Parameters);
template void MoveMeshUtility::FillVariablesList(std::vector<const Variable<array_1d<double,3>>*>&, Parameters);

template void MoveMeshUtility::InterpolateVariable(NodeType&, const Vector&, const GeometryType&, const Variable<double>&);
template void MoveMeshUtility::InterpolateVariable(NodeType&, const Vector&, const GeometryType&, const Variable<array_1d<double,3>>&);

}  // namespace Kratos.

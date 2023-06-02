//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

#pragma once

// System incldues
#include <variant>

// Project includes
#include "containers/array_1d.h"
#include "containers/container_expression/container_expression.h"
#include "containers/container_expression/container_data_io.h"
#include "containers/variable.h"
#include "expression_io.h"
#include "includes/define.h"
#include "includes/model_part.h"

namespace Kratos {


class KRATOS_API(KRATOS_CORE) VariableExpressionIO: public ExpressionIO
{
public:
    ///@name Type definitions
    ///@{

    using VariableType = std::variant<
                                const Variable<int>*,
                                const Variable<double>*,
                                const Variable<array_1d<double, 3>>*,
                                const Variable<array_1d<double, 4>>*,
                                const Variable<array_1d<double, 6>>*,
                                const Variable<array_1d<double, 9>>*,
                                const Variable<Vector>*,
                                const Variable<Matrix>*>;

    KRATOS_CLASS_POINTER_DEFINITION(VariableExpressionIO);

    ///}
    ///@name Public enums
    ///@{

    enum ContainerType
    {
        NodalHistorical,
        NodalNonHistorical,
        ConditionNonHistorical,
        ElementNonHistorical
    };

    enum MeshType
    {
        Local,
        Interface,
        Ghost
    };

    ///@}
    ///@name Life cycle
    ///@{

    template<class TDataType>
    VariableExpressionIO(
        ModelPart& rModelPart,
        const Variable<TDataType>& rVariable,
        const ContainerType& rContainerType,
        const MeshType& rMeshType = MeshType::Local);

    template <class TDataType, class TMeshType = Kratos::MeshType::Local>
    VariableExpressionIO(
        ContainerExpression<ModelPart::NodesContainerType, TMeshType>& rContainer,
        const Variable<TDataType>& rVariable,
        const bool IsHistorical = false);

    template <class TContainerType, class TDataType, class TMeshType = Kratos::MeshType::Local>
    VariableExpressionIO(
        ContainerExpression<TContainerType, TMeshType>& rContainer,
        const Variable<TDataType>& rVariable);

    ~VariableExpressionIO() override = default;

    ///@}
    ///@name Public operations
    ///@{

    Expression::Pointer Read() override;

    void Write(const Expression& rExpression) override;

    ///@}

private:
    ///@name Private member variables
    ///@{

    ModelPart& mrModelPart;

    const VariableType mpVariable;

    const ContainerType mContainerType;

    const MeshType mMeshType;

    ///@}
    ///@name Private operations
    ///@{

    ModelPart::MeshType& GetMesh();

    ///@}
}; // class ExpressionIO


} // namespace Kratos

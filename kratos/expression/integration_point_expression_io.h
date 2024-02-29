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
#include "expression/container_expression.h"
#include "expression/traits.h"
#include "containers/variable.h"
#include "expression_io.h"
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/global_variables.h"

namespace Kratos {


class KRATOS_API(KRATOS_CORE) IntegrationPointExpressionIO
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


    ///}
    ///@name Public classes
    ///@{

    class KRATOS_API(KRATOS_CORE) Input : public ExpressionInput
    {
    public:
        ///@name Type definitions
        ///@{

        KRATOS_CLASS_POINTER_DEFINITION(Input);

        ///@}
        ///@name Life cycle
        ///@{

        Input(
            ModelPart& rModelPart,
            const VariableType& rVariable,
            Globals::DataLocation CurrentLocation,
            MeshType CurrentMeshType = MeshType::Local);

        ~Input() override = default;

        ///@}
        ///@name Public operations
        ///@{

        Expression::Pointer Execute() const override;

        ///@}

    private:
        ///@name Private member variables
        ///@{

        ModelPart* mpModelPart;

        VariableType mpVariable;

        Globals::DataLocation mDataLocation;

        MeshType mMeshType;

        ///@}

    };

    class KRATOS_API(KRATOS_CORE) Output : public ExpressionOutput
    {
    public:
        ///@name Type definitions
        ///@{

        KRATOS_CLASS_POINTER_DEFINITION(Output);

        ///@}
        ///@name Life cycle
        ///@{

        Output(
            ModelPart& rModelPart,
            const VariableType& rVariable,
            Globals::DataLocation CurrentLocation,
            MeshType CurrentMeshType = MeshType::Local);

        ~Output() override = default;

        ///@}
        ///@name Public operations
        ///@{

        void Execute(const Expression& rExpression) override;

        ///@}

    private:
        ///@name Private member variables
        ///@{

        ModelPart* mpModelPart;

        VariableType mpVariable;

        Globals::DataLocation mDataLocation;

        MeshType mMeshType;

        ///@}

    };

    ///@}
    ///@name Public static operations
    ///@{

    template<class TContainerType, MeshType TMeshType>
    KRATOS_API(KRATOS_CORE) static void Read(
        ContainerExpression<TContainerType, TMeshType>& rContainerExpression,
        const VariableType& rVariable);

    template<class TContainerType, MeshType TMeshType>
    KRATOS_API(KRATOS_CORE) static void Write(
        const ContainerExpression<TContainerType, TMeshType>& rContainerExpression,
        const VariableType& rVariable);

    ///@}

}; // class ExpressionIO


} // namespace Kratos

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

// Project includes
#include "expression/container_expression.h"
#include "expression/traits.h"
#include "expression_io.h"
#include "includes/define.h"
#include "includes/global_variables.h"
#include "includes/model_part.h"

namespace Kratos {


class KRATOS_API(KRATOS_CORE) NodalPositionExpressionIO
{
public:
    ///@name Type definitions
    ///@{

    using Configuration = Kratos::Globals::Configuration;

    ///@}
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
            const ModelPart& rModelPart,
            const Configuration& rConfiguration,
            const MeshType& rMeshType = MeshType::Local);

        ~Input() override = default;

        ///@}
        ///@name Public operations
        ///@{

        Expression::Pointer Execute() const override;

        ///@}

    private:
        ///@name Private member variables
        ///@{

        const ModelPart& mrModelPart;

        const Configuration mConfiguration;

        const MeshType mMeshType;

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
            const Configuration& rConfiguration,
            const MeshType& rMeshType = MeshType::Local);

        ~Output() override = default;

        ///@}
        ///@name Public operations
        ///@{

        void Execute(const Expression& rExpression) override;

        ///@}

    private:
        ///@name Private member variables
        ///@{

        ModelPart& mrModelPart;

        const Configuration mConfiguration;

        const MeshType mMeshType;

        ///@}

    };

    ///@}
    ///@name Public static operations
    ///@{

    template<MeshType TMeshType = MeshType::Local>
    KRATOS_API(KRATOS_CORE) static void Read(
        ContainerExpression<ModelPart::NodesContainerType, TMeshType>& rContainerExpression,
        const Configuration& rConfiguration);

    template<MeshType TMeshType = MeshType::Local>
    KRATOS_API(KRATOS_CORE) static void Write(
        const ContainerExpression<ModelPart::NodesContainerType, TMeshType>& rContainerExpression,
        const Configuration& rConfiguration);

    ///@}

}; // class ExpressionIO


} // namespace Kratos

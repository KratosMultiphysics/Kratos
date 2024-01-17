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
#include "includes/model_part.h"

namespace Kratos {


class KRATOS_API(KRATOS_CORE) EntityDomainSizeExpressionIO
{
public:
    ///@name Type definitions
    ///@{

    using Configuration = Kratos::Globals::Configuration;

    ///@}
    ///@name Public classes
    ///@{

    class KRATOS_API(KRATOS_CORE) EntityDomainSizeExpressionInput : public ExpressionInput
    {
    public:
        ///@name Type definitions
        ///@{

        KRATOS_CLASS_POINTER_DEFINITION(EntityDomainSizeExpressionInput);

        ///@}
        ///@name Life cycle
        ///@{

        EntityDomainSizeExpressionInput(
            const ModelPart& rModelPart,
            const ContainerType& rContainerType);

        ~EntityDomainSizeExpressionInput() override = default;

        ///@}
        ///@name Public operations
        ///@{

        Expression::Pointer Execute() const override;

        ///@}

    private:
        ///@name Private member variables
        ///@{

        const ModelPart& mrModelPart;

        const ContainerType mContainerType;

        ///@}

    };

    ///@}
    ///@name Public static operations
    ///@{

    template<class TContainerType, MeshType TMeshType>
    KRATOS_API(KRATOS_CORE) static void Read(ContainerExpression<TContainerType, TMeshType>& rContainerExpression);

    ///@}

}; // class ExpressionIO


} // namespace Kratos

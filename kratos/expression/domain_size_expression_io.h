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
#include "includes/global_variables.h"

namespace Kratos {

/**
 * @brief Utility class for domain size IO
 * @details This utility is used to compute domain sizes of each entity
 *          (such as of conditions/elements). The domain size for conditions
 *          and elements differs depending on the underlying geometry.
 *          1. If it is a line geometry, then domain size represents the line length.
 *          2. If it is a surface geometry, then domain size represents the surface area.
 *          3. If is is a volume geometry, then domain size represents the volume.
 * @see Geometry::DomainSize
 */
class KRATOS_API(KRATOS_CORE) DomainSizeExpressionIO
{
public:
    ///@name Public classes
    ///@{

    /**
     * @brief Construct a new domain expression input
     * @details This constructs an input object to read in domain sizes of each entity.
     *
     */
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
            Globals::DataLocation CurrentLocation);

        ///@}
        ///@name Public operations
        ///@{

        Expression::Pointer Execute() const override;

        ///@}

    private:
        ///@name Private member variables
        ///@{

        ModelPart const * mpModelPart;

        Globals::DataLocation mDataLocation;

        ///@}

    };

    ///@}
    ///@name Public static operations
    ///@{

    /// @tparam TContainerType @ref ModelPart::ElementsContainerType or @ref ModelPart::ConditionsContainerType
    template<class TContainerType, MeshType TMeshType>
    KRATOS_API(KRATOS_CORE) static void Read(ContainerExpression<TContainerType, TMeshType>& rContainerExpression);

    ///@}

}; // class ExpressionIO


} // namespace Kratos

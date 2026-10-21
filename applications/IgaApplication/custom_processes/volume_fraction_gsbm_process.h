//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Mario Caballero
//

#pragma once

// System includes

// External includes

// Project includes
#include "containers/model.h"
#include "includes/define.h"
#include "processes/process.h"

namespace Kratos
{

///@name Kratos Classes
///@{

class KRATOS_API(IGA_APPLICATION) VolumeFractionGsbmProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(VolumeFractionGsbmProcess);

    using SizeType = std::size_t;
    using IndexType = std::size_t;
    using GeometryType = Geometry<Node>;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    VolumeFractionGsbmProcess(
        Model& rModel,
        Parameters ThisParameters)
        : Process()
        , mpModel(&rModel)
        , mParameters(ThisParameters)
    {
        mParameters.ValidateAndAssignDefaults(GetDefaultParameters());
    }

    /// Destructor.
    ~VolumeFractionGsbmProcess() override = default;

    // Get the default parameters
    const Parameters GetDefaultParameters() const override;

    ///@}
    ///@name Operations
    ///@{

    /// Assigns VOLUME_FRACTION to the elements and facet sides. Each cut facet
    /// gets NEIGHBOUR_NODES (a node at its closest skin point) and IDENTIFIER,
    /// as in SBM. The node holds the skin's NORMAL, LOCAL_TANGENT, CURVATURE
    /// and CURVATURE_DERIVATIVES there. NURBS skins are tessellated into
    /// <skin>.tessellation, and projected onto exactly with SBM's
    /// ProjectToNurbsSkin. Skin loops must be counterclockwise, holes
    /// included, as in SBM.
    /// Assumes uv = xy
    void ExecuteInitialize() override;

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "VolumeFractionGsbmProcess";
    }

    ///@}

private:
    ///@name Member Variables
    ///@{

    Model* mpModel = nullptr;
    Parameters mParameters;

    ///@}
};

///@}

} // namespace Kratos

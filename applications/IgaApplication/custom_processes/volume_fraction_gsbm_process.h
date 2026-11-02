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
    /// With weight_approach "binary", VOLUME_FRACTION is 1 on the knot spans
    /// fully inside and 0 on the others: the SBM as a GSBM.
    /// Entities that add nothing are set not ACTIVE: the elements of empty knot
    /// spans, the facets without a jump in VOLUME_FRACTION and, with
    /// ghost_penalty_model_part_name, the ghost penalty conditions with
    /// \beta_2 = 0. So skeleton_model_part_name must hold the conditions of the
    /// shifted boundary terms, as elements must be scaled by VOLUME_FRACTION.
    /// Assumes uv = xy
    void ExecuteInitialize() override;

    /// With update_each_step, for moving skins: removes the projection nodes,
    /// the cut facets' NEIGHBOUR_NODES and IDENTIFIER, and the tessellations,
    /// then runs ExecuteInitialize for the skins where they are now
    void ExecuteInitializeSolutionStep() override;

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

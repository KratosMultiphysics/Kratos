//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Nicolò Antonelli
//

#pragma once

// Project includes
#include "modeler/surrogate_boundary_modeler.h"

namespace Kratos {

/**
 * @class GapSurrogateBoundaryModeler
 * @brief Adds a quadrature representation of the SBM gap.
 * @details The base SurrogateBoundaryModeler owns classification, active-mesh
 *          generation, projection and surrogate-face topology. This derived
 *          modeler uses those results to create, in 2D:
 *          - at most one quadrature-point gap element per surrogate face;
 *          - boundary quadrature on the reconstructed physical skin;
 *          - one coupling condition on every non-trivial interface between
 *            adjacent gap patches whose basis extensions come from different
 *            parent elements.
 *
 *          Gap entities contain the nodes of their active parent element. No
 *          fictitious nodes or DOFs are introduced: the parent FEM basis is
 *          evaluated outside the parent cell at every gap quadrature point.
 *          The physics remains application-owned through registered element
 *          and condition names supplied in the settings.
 *
 *          A surrogate face can collapse to zero gap measure at an exact
 *          grid/skin tangency or in the fitted limit. Its physical-boundary
 *          condition is retained (with its reconstructed quadrature, possibly
 *          of zero line measure), while no zero-volume gap element is created.
 *          The face still participates in topology so adjacent positive-area
 *          patches are coupled to the correct active parent element.
 *
 * @note The first implementation supports a 2D Quad4/Line2 workflow.
 *       It reconstructs each true-boundary piece as the straight segment
 *       joining the closest-point projections of the surrogate endpoints.
 */
class KRATOS_API(KRATOS_CORE) GapSurrogateBoundaryModeler
    : public SurrogateBoundaryModeler
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(GapSurrogateBoundaryModeler);

    GapSurrogateBoundaryModeler() = default;

    GapSurrogateBoundaryModeler(
        Model& rModel,
        Parameters ModelerParameters = Parameters());

    ~GapSurrogateBoundaryModeler() override = default;

    Modeler::Pointer Create(
        Model& rModel,
        const Parameters ModelParameters) const override
    {
        return Kratos::make_shared<GapSurrogateBoundaryModeler>(
            rModel, ModelParameters);
    }

    void SetupModelPart() override;

    /**
     * @brief Rebuilds the active, surrogate and gap views after skin motion.
     * @details The persistent Cartesian background is owned by the base
     *          modeler. Gap quadrature elements and their interface conditions
     *          are transient topology and are recreated at every update.
     */
    void PrepareGeometryModel() override;

    const Parameters GetDefaultParameters() const override;

    std::string Info() const override
    {
        return "GapSurrogateBoundaryModeler";
    }

private:
    Model* mpGapModel = nullptr;
    Parameters mGapParameters;

    void ClearGapModelPart2D();
    void SetupGapModelPart2D();
};

} // namespace Kratos

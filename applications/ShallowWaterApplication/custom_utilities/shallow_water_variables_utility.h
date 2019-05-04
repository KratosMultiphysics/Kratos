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

#ifndef KRATOS_SHALLOW_WATER_VARIABLES_UTILITY_H_INCLUDED
#define KRATOS_SHALLOW_WATER_VARIABLES_UTILITY_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/model_part.h"
#include "shallow_water_application.h"

namespace Kratos
{

typedef std::size_t IndexType;

typedef std::vector<IndexType> IndexVectorType;

typedef std::map<int, Properties::Pointer> PropertiesMapType;

class ShallowWaterVariablesUtility
{
public:

    KRATOS_CLASS_POINTER_DEFINITION(ShallowWaterVariablesUtility);

    ShallowWaterVariablesUtility(ModelPart& rModelPart, const double& rDryHeight = 1e-3)
        : mrModelPart(rModelPart)
    {
        KRATOS_TRY

        std::cout << "Initializing shallow water variables utility" << std::endl;
        mWaterHeightConvert = mrModelPart.GetProcessInfo()[WATER_HEIGHT_UNIT_CONVERTER];
        mDryHeight = rDryHeight;
        mZeroValue = 1e-8;

        KRATOS_CATCH("")
    }

    ~ShallowWaterVariablesUtility()
    {}

    /**
     * This method computes the free surface elevation as HEIGHT + BATHYMETRY
     */
    void ComputeFreeSurfaceElevation()
    {
        KRATOS_TRY

        const int nnodes = static_cast<int>(mrModelPart.Nodes().size());
        ModelPart::NodesContainerType::iterator nodes_begin = mrModelPart.NodesBegin();
        #pragma omp parallel for
        for(int i = 0; i < nnodes; i++)
        {
            auto node = nodes_begin + i;
            node->FastGetSolutionStepValue(FREE_SURFACE_ELEVATION) = node->FastGetSolutionStepValue(HEIGHT) + (node->FastGetSolutionStepValue(BATHYMETRY) / mWaterHeightConvert);
        }

        KRATOS_CATCH("")
    }

    /**
     * This method computes the water height as FREE_SURFACE_ELEVATION - BATHYMETRY
     */
    void ComputeHeightFromFreeSurface()
    {
        KRATOS_TRY

        const int nnodes = static_cast<int>(mrModelPart.Nodes().size());
        ModelPart::NodesContainerType::iterator nodes_begin = mrModelPart.NodesBegin();

        #pragma omp parallel for
        for(int i = 0; i < nnodes; i++)
        {
            auto node = nodes_begin + i;
            node->FastGetSolutionStepValue(HEIGHT) = node->FastGetSolutionStepValue(FREE_SURFACE_ELEVATION) - (node->FastGetSolutionStepValue(BATHYMETRY) / mWaterHeightConvert);
        }

        KRATOS_CATCH("")
    }

    /**
     * This method computes the velocity as the MOMENTUM / HEIGHT
     */
    void ComputeVelocity()
    {
        KRATOS_TRY

        const int nnodes = static_cast<int>(mrModelPart.Nodes().size());
        ModelPart::NodesContainerType::iterator nodes_begin = mrModelPart.NodesBegin();
        #pragma omp parallel for
        for(int i = 0; i < nnodes; i++)
        {
            auto node = nodes_begin + i;
            node->GetSolutionStepValue(VELOCITY) = node->FastGetSolutionStepValue(MOMENTUM) / (node->FastGetSolutionStepValue(HEIGHT) * mWaterHeightConvert);
        }

        KRATOS_CATCH("")
    }

    /**
     * This method computes the momentum as the VELOCITY * HEIGHT
     */
    void ComputeMomentum()
    {
        KRATOS_TRY

        const int nnodes = static_cast<int>(mrModelPart.Nodes().size());
        ModelPart::NodesContainerType::iterator nodes_begin = mrModelPart.NodesBegin();

        #pragma omp parallel for
        for(int i = 0; i < nnodes; i++)
        {
            auto node = nodes_begin + i;
            node->GetSolutionStepValue(MOMENTUM) = node->FastGetSolutionStepValue(VELOCITY) * (node->FastGetSolutionStepValue(HEIGHT) * mWaterHeightConvert);
        }

        KRATOS_CATCH("")
    }

    /**
     * Set the dry nodes to zero
     */
    void CheckDryConservedVariables()
    {
        KRATOS_TRY

        const int nnodes = static_cast<int>(mrModelPart.Nodes().size());
        ModelPart::NodesContainerType::iterator nodes_begin = mrModelPart.NodesBegin();

        #pragma omp parallel for
        for(int i = 0; i < nnodes; i++)
        {
            auto node = nodes_begin + i;

            if (node->FastGetSolutionStepValue(HEIGHT) < mDryHeight &&
                node->FastGetSolutionStepValue(RAIN)   < mDryHeight )
            {
                node->FastGetSolutionStepValue(HEIGHT)     = mZeroValue;
                node->FastGetSolutionStepValue(MOMENTUM_X) = 0;
                node->FastGetSolutionStepValue(MOMENTUM_Y) = 0;
            }
        }
        KRATOS_CATCH("")
    }

    /**
     * Set the dry nodes to zero
     */
    void CheckDryPrimitiveVariables()
    {
        KRATOS_TRY

        const int nnodes = static_cast<int>(mrModelPart.Nodes().size());
        ModelPart::NodesContainerType::iterator nodes_begin = mrModelPart.NodesBegin();

        #pragma omp parallel for
        for(int i = 0; i < nnodes; i++)
        {
            auto node = nodes_begin + i;

            if (node->FastGetSolutionStepValue(HEIGHT) < mDryHeight &&
                node->FastGetSolutionStepValue(RAIN)   < mDryHeight )
            {
                node->FastGetSolutionStepValue(HEIGHT)     = mZeroValue;
                node->FastGetSolutionStepValue(VELOCITY_X) = 0;
                node->FastGetSolutionStepValue(VELOCITY_Y) = 0;
            }
        }
        KRATOS_CATCH("")
    }

    /**
     * This method looks for the dry elements and sets the proper flag
     */
    void SetDryWetState()
    {
        KRATOS_TRY

        // Getting the elements from the model
        const int nelem = static_cast<int>(mrModelPart.Elements().size());
        ModelPart::ElementsContainerType::iterator elem_begin = mrModelPart.ElementsBegin();

        // And now, if an element has all nodes dry, it is not active
        #pragma omp parallel for
        for(int k = 0; k < nelem; k++)
        {
            auto elem = elem_begin + k;
            bool wet_node = false;                     // The nodal flag
            for(Node<3>& node : elem->GetGeometry())
            {
                double dry_height = 0.2 * std::sqrt(node.FastGetSolutionStepValue(NODAL_AREA)) * norm_2(node.FastGetSolutionStepValue(TOPOGRAPHY_GRADIENT));
                dry_height = std::max(dry_height, mDryHeight);
                if (node.FastGetSolutionStepValue(HEIGHT) >= mDryHeight ||
                    node.FastGetSolutionStepValue(RAIN)   >= mDryHeight )
                    wet_node = true;  // It means there is almost a wet node
            }

            if (wet_node)
                elem->Set(FLUID, true);
            else
                elem->Set(FLUID, false);
        }

        KRATOS_CATCH("")
    }

    void DeactivateDryElements()
    {
        KRATOS_TRY

        // Getting the elements from the model
        const int nelem = static_cast<int>(mrModelPart.Elements().size());
        ModelPart::ElementsContainerType::iterator elem_begin = mrModelPart.ElementsBegin();

        #pragma omp parallel for
        for (int k = 0; k < nelem; k++)
        {
            auto elem = elem_begin + k;

            if (elem->Is(FLUID))
                elem->Set(ACTIVE, true);
            else
                elem->Set(ACTIVE, false);
        }

        KRATOS_CATCH("")
    }

protected:

private:

    ModelPart& mrModelPart;
    double mWaterHeightConvert;
    double mDryHeight;
    double mZeroValue;

}; // class ShallowWaterVariablesUtility

} // namespace Kratos

# endif // KRATOS_SHALLOW_WATER_VARIABLES_UTILITY_H_INCLUDED

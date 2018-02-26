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

#if !defined(KRATOS_SHALLOW_WATER_VARIABLES_UTILITY_H_INCLUDED)
#define  KRATOS_SHALLOW_WATER_VARIABLES_UTILITY_H_INCLUDED

// System includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "utilities/openmp_utils.h"
#include "processes/find_nodal_neighbours_process.h"

// External includes 

// Project includes
#include "shallow_water_application.h"

namespace Kratos
{

    class ShallowWaterVariablesUtility
    {
    public:

        KRATOS_CLASS_POINTER_DEFINITION(ShallowWaterVariablesUtility);

        ShallowWaterVariablesUtility(ModelPart& model_part) :
            mrModelPart(model_part)  
        {
            KRATOS_TRY
            
            std::cout << "Initializing shallow water variables utility" << std::endl; 
            mWaterHeightConvert = mrModelPart.GetProcessInfo()[WATER_HEIGHT_UNIT_CONVERTER];
            mThreshold = 1e-3;
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

            ModelPart::NodesContainerType& r_nodes = mrModelPart.Nodes();
            #pragma omp parallel for
            for(unsigned int i = 0; i < static_cast<unsigned int>(r_nodes.size()); i++)
            {
                ModelPart::NodesContainerType::iterator inode = r_nodes.begin() + i;
                inode->FastGetSolutionStepValue(FREE_SURFACE_ELEVATION) = inode->FastGetSolutionStepValue(HEIGHT) + (inode->FastGetSolutionStepValue(BATHYMETRY) / mWaterHeightConvert);
            }

            KRATOS_CATCH("")
        }

        /**
         * This method computes the water height as FREE_SURFACE_ELEVATION - BATHYMETRY
         */
        void ComputeHeightFromInitialFreeSurface()
        {
            KRATOS_TRY

            ModelPart::NodesContainerType& r_nodes = mrModelPart.Nodes();
            #pragma omp parallel for
            for(unsigned int i = 0; i < static_cast<unsigned int>(r_nodes.size()); i++)
            {
                ModelPart::NodesContainerType::iterator inode = r_nodes.begin() + i;
                inode->FastGetSolutionStepValue(HEIGHT) = inode->FastGetSolutionStepValue(FREE_SURFACE_ELEVATION) - (inode->FastGetSolutionStepValue(BATHYMETRY) / mWaterHeightConvert);
            }

            KRATOS_CATCH("")
        }

        /** 
         * This method computes the velocity as the MOMENTUM / HEIGHT
         */
        void ComputeVelocity()
        {
            KRATOS_TRY

            ModelPart::NodesContainerType& r_nodes = mrModelPart.Nodes();
            #pragma omp parallel for
            for(unsigned int i = 0; i < static_cast<unsigned int>(r_nodes.size()); i++)
            {
                ModelPart::NodesContainerType::iterator inode = r_nodes.begin() + i;
                inode->GetSolutionStepValue(VELOCITY) = inode->FastGetSolutionStepValue(MOMENTUM) / (inode->FastGetSolutionStepValue(HEIGHT) * mWaterHeightConvert);
            }

            KRATOS_CATCH("")
        }

        /** 
         * This method computes the momentum as the VELOCITY * HEIGHT
         */
        void ComputeMomentum()
        {
            KRATOS_TRY

            ModelPart::NodesContainerType& r_nodes = mrModelPart.Nodes();
            #pragma omp parallel for
            for(unsigned int i = 0; i < static_cast<unsigned int>(r_nodes.size()); i++)
            {
                ModelPart::NodesContainerType::iterator inode = r_nodes.begin() + i;
                inode->GetSolutionStepValue(MOMENTUM) = inode->FastGetSolutionStepValue(VELOCITY) * (inode->FastGetSolutionStepValue(HEIGHT) * mWaterHeightConvert);
            }

            KRATOS_CATCH("")
        }

        void CheckDryConservedVariables()
        {
            KRATOS_TRY

            ModelPart::NodesContainerType& r_nodes = mrModelPart.Nodes();
            #pragma omp parallel for
            for(unsigned int i = 0; i < static_cast<unsigned int>(r_nodes.size()); i++)
            {
                ModelPart::NodesContainerType::iterator inode = r_nodes.begin() + i;
                if (inode->FastGetSolutionStepValue(HEIGHT) < mThreshold &&
                    inode->FastGetSolutionStepValue(RAIN)   < mThreshold )
                {
                    inode->FastGetSolutionStepValue(HEIGHT)     = mZeroValue;
                    inode->FastGetSolutionStepValue(MOMENTUM_X) = 0;
                    inode->FastGetSolutionStepValue(MOMENTUM_Y) = 0;
                }
            }
            KRATOS_CATCH("")
        }

        void CheckDryPrimitiveVariables()
        {
            KRATOS_TRY

            ModelPart::NodesContainerType& r_nodes = mrModelPart.Nodes();
            #pragma omp parallel for
            for(unsigned int i = 0; i < static_cast<unsigned int>(r_nodes.size()); i++)
            {
                ModelPart::NodesContainerType::iterator inode = r_nodes.begin() + i;
                if (inode->FastGetSolutionStepValue(HEIGHT) < mThreshold &&
                    inode->FastGetSolutionStepValue(RAIN)   < mThreshold )
                {
                    inode->FastGetSolutionStepValue(HEIGHT)     = mZeroValue;
                    inode->FastGetSolutionStepValue(VELOCITY_X) = 0;
                    inode->FastGetSolutionStepValue(VELOCITY_Y) = 0;
                }
            }
            KRATOS_CATCH("")
        }

        void SetDryWetState()
        {
            KRATOS_TRY
            
            //~ int expected_neigh_elems = 6;
            //~ int expected_neigh_nodes = 6;
            //~ FindNodalNeighboursProcess find_nodes_process(mrModelPart, expected_neigh_elems, expected_neigh_nodes);
            //~ find_nodes_process.Execute();

            //~ ModelPart::NodesContainerType& r_nodes = mrModelPart.Nodes();
            //~ // We loop all the nodes to check if they are dry
            //~ #pragma omp parallel for
            //~ for(unsigned int i = 0; i < static_cast<unsigned int>(r_nodes.size()); i++)
            //~ {
                //~ ModelPart::NodesContainerType::iterator inode = r_nodes.begin() + i;
                //~ // If current node is dry, is candidate to be inactive
                //~ if (inode->FastGetSolutionStepValue(HEIGHT) < mThreshold && 
                    //~ inode->FastGetSolutionStepValue(RAIN)   < mThreshold )
                //~ {
                    //~ WeakPointerVector< Node<3> >& rneigh = inode->GetValue(NEIGHBOUR_NODES);
                    //~ // We loop all the neighbour nodes to check if they are dry
                    //~ // If a neighbour node is wet, current node is candidate to be wet, so it is active
                    //~ bool neigh_wet = false;
                    //~ for( WeakPointerVector<Node<3> >::iterator jnode = rneigh.begin(); jnode!=rneigh.end(); jnode++)
                    //~ {
                        //~ if (jnode->FastGetSolutionStepValue(HEIGHT) >= mThreshold ||
                            //~ jnode->FastGetSolutionStepValue(RAIN)   >= mThreshold )
                            //~ neigh_wet = true;
                    //~ }
                    //~ if (neigh_wet)
                        //~ inode->Set(ACTIVE, true);
                    //~ else
                        //~ inode->Set(ACTIVE, false);
                //~ }
                //~ // If current element is wet, set active
                //~ else
                    //~ inode->Set(ACTIVE, true);
            //~ }
            
            // Way B: elements
            
            // Getting the elements from the model
            const unsigned int nelements = static_cast<int>(mrModelPart.Elements().size());
            int nnodes;
            bool wet_node;
            
            // And now, if an element has all nodes dry, it is not active
            #pragma omp parallel for
            for(unsigned int k = 0; k < nelements; k++)
            {
                ModelPart::ElementsContainerType::iterator it = mrModelPart.ElementsBegin() + k;
                nnodes = it->GetGeometry().size();
                wet_node = false;
                for(int l = 0; l < nnodes; l++)
                {
                    if (it->GetGeometry()[l].FastGetSolutionStepValue(HEIGHT) >= mThreshold ||
                        it->GetGeometry()[l].FastGetSolutionStepValue(RAIN)   >= mThreshold )
                        wet_node = true;  // It means there is almost a wet node
                }
                if (wet_node)
                    it->Set(ACTIVE, true);
                else
                    it->Set(ACTIVE, false);
            }
            
            KRATOS_CATCH("")
        }

    protected:

    private:

        ModelPart& mrModelPart;
        double mWaterHeightConvert;
        double mThreshold;
        double mZeroValue;

    }; // class ShallowWaterVariablesUtility

} // namespace Kratos

# endif // KRATOS_SHALLOW_WATER_VARIABLES_UTILITY_H_INCLUDED

/*
==============================================================================
KratosShallowWaterApplication 
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu 
rrossi@cimne.upc.edu
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain


Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE. 

==============================================================================
*/
 
//   
//   Project Name:        Kratos       
//   Last Modified by:    Miguel Masó Sotomayor
//   Date:                September 15th 2017
//   Revision:            1.0
//
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

/*
==============================================================================
KratosIncompressibleFluidApplication
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
//   Last Modified by:    $Author: rrossi $
//   Date:                $Date: 2008-10-13 06:58:23 $
//   Revision:            $Revision: 1.4 $
//
//


#if !defined(KRATOS_MARK_FOR_REFINEMENT )
#define  KRATOS_MARK_FOR_REFINEMENT



// System includes
#include <string>
#include <iostream>
#include <algorithm>

// External includes


// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/node.h"
#include "utilities/geometry_utilities.h"
#include "utilities/openmp_utils.h"
#include "incompressible_fluid_application.h"

namespace Kratos
{
    ///@addtogroup IncompressibleFluidApplication
    ///@{

    ///@name Kratos Classes
    ///@{

    /// A class for dynamic mesh refinement.
    /**
     Provides several utility functions to refine the mesh during the solution
     procedure. The main function is MarkForRefinement, which chooses which elements
     to refine based on the value of an elemental variable and a user-defined tolerance
     */
    class RefinementUtilities
    {
    public:
        ///@name Life Cycle
        ///@{

        /// Default constructor
        RefinementUtilities()
        {}

        /// Destructor
        ~RefinementUtilities(){}

        ///@}
        ///@name Operations
        ///@{

        /// Identify the elements that need to be refined.
        /**
         This class mark elements for refinement based on the elemental value of the variable
         given as input. If the value exceeds the given tolerance, the element will be marked for
         refinement. If an element is refined several times, this procedure may choose to refine
         its neighbours too, to avoid creating elements with poor shapes.
         @param rVariable The variable used to check which elements must be refined. It is assumed
         that another function will write, to allow for different refinement criteria
         @param ThisModelPart The model part containig the elements we want to refine
         @param DomainSize Number of spatial dimensions of the problem (2 or 3)
         @param admissible_ratio The refinement tolerance. All elements with a greater value of
         rVariable will be refined
         @param admissible_area Minimum area (or volume) for the new elements. Elements won't be
         refined if this results in elements smaller than this size.
         @param max_levels Maximum number of successive refinements on the same element
         */
        void MarkForRefinement(Variable<double>& rVariable,
                               ModelPart& ThisModelPart,
                               unsigned int DomainSize,
                               double admissible_ratio,
                               double admissible_area,
                               int max_levels)
        {
            KRATOS_TRY;

            //reset the nodal values
            for(ModelPart::NodesContainerType::iterator it=ThisModelPart.NodesBegin(); it!=ThisModelPart.NodesEnd(); it++)
            {
                it->GetValue(REFINEMENT_LEVEL) = 0;
            }

            //mark elements for splitting depending on the desired error ratio
            unsigned int number_of_splitted_elements = 0;
            for(ModelPart::ElementsContainerType::iterator it=ThisModelPart.ElementsBegin(); it!=ThisModelPart.ElementsEnd(); it++)
            {
                double ratio = it->GetValue(rVariable);

                if(ratio > admissible_ratio && it->GetValue(REFINEMENT_LEVEL) <= max_levels )
                {
                    //mark for splitting
                    it->GetValue(SPLIT_ELEMENT) = true;
                    int& current_level = it->GetValue(REFINEMENT_LEVEL);
                    current_level += 1;

                    number_of_splitted_elements++;

                    //mark all of the nodes with the refinement level
                    Geometry< Node<3> >& geom = it->GetGeometry();

                    for(unsigned int i=0; i<geom.size(); i++)
                    {
                        if(geom[i].GetValue(REFINEMENT_LEVEL) < current_level)
                        {
                            geom[i].GetValue(REFINEMENT_LEVEL) = current_level;
                        }
                    }
                }
            }

            bool is_ok = false;
            int nit = 0;
            while(is_ok==false && nit<2*max_levels)
            {
                is_ok = true;

                ModelPart::ElementsContainerType aux_elem_list;

                //fill a list with all of the elements that have to be refined to obtain a correct gradient
                for(ModelPart::ElementsContainerType::iterator it=ThisModelPart.ElementsBegin(); it!=ThisModelPart.ElementsEnd(); it++)
                {
                    //mark all of the nodes with the refinement level
                    Geometry< Node<3> >& geom = it->GetGeometry();

                    //determine for each element the maximum level of refinement (basically if the neighbours have been refined)
                    int level = geom[0].GetValue(REFINEMENT_LEVEL);
//                    int min_level = level;
                    int max_level = level;
                    for(unsigned int i=1; i<geom.size(); i++)
                    {
                        level = geom[i].GetValue(REFINEMENT_LEVEL);
//                        if(level < min_level) min_level = level;
                        if(level > max_level) max_level = level;
                    }

                    const int& current_level = it->GetValue(REFINEMENT_LEVEL);
                    //if there is a difference of level of refinement greater than 1, then refine to have a smoother grading of elements
                    if(max_level > current_level + 1 && current_level < max_levels)
                    {
                        //more iterations of the overall algorithm are needed
                        is_ok = false;

                        aux_elem_list.push_back( *(it.base()) );
                    }
                }

                //now signal such elements for spltting and color their nodes correctly
                for(ModelPart::ElementsContainerType::iterator it=aux_elem_list.begin(); it!=aux_elem_list.end(); it++)
                {
                    int& current_level = it->GetValue(REFINEMENT_LEVEL);
                    //mark for splitting
                    it->GetValue(SPLIT_ELEMENT) = true;
                    current_level += 1; //here we increase the level

                    number_of_splitted_elements++;

                    //mark all of the nodes with the refinement level
                    Geometry< Node<3> >& geom = it->GetGeometry();

                    for(unsigned int i=0; i<geom.size(); i++)
                    {
                        if(geom[i].GetValue(REFINEMENT_LEVEL) < current_level)
                        {
                            geom[i].GetValue(REFINEMENT_LEVEL) = current_level;
                        }
                    }
                }

                //increase iterations of the grading algorithm
                nit += 1;
            }

            double MinSize,ElemSize;

            if (DomainSize == 2)
            {
                MinSize = admissible_area * 4.0;

                for( ModelPart::ElementIterator itElem = ThisModelPart.ElementsBegin(); itElem != ThisModelPart.ElementsEnd(); ++itElem)
                {
                    ElemSize = GeometryUtils::CalculateVolume2D(itElem->GetGeometry());
                    if (ElemSize < MinSize)
                    {
                        itElem->SetValue(SPLIT_ELEMENT,false);
                        number_of_splitted_elements -= 1;
                    }
                }
            }
            else // (DomainSize == 3)
            {
                MinSize = admissible_area * 8.0;

                for( ModelPart::ElementIterator itElem = ThisModelPart.ElementsBegin(); itElem != ThisModelPart.ElementsEnd(); ++itElem)
                {
                    ElemSize = GeometryUtils::CalculateVolume3D(itElem->GetGeometry());
                    if (ElemSize < MinSize)
                    {
                        itElem->SetValue(SPLIT_ELEMENT,false);
                        number_of_splitted_elements -= 1;
                    }
                }
            }

            KRATOS_WATCH("***********************************************************");
            std::cout << "total number of refined elements = " << number_of_splitted_elements << std::endl;

            KRATOS_CATCH("")
        }

        /// Calculate the subscale velocity on the model's mesh and use it to mark elements for refinement.
        /**
         This function is intended to work with VMS elements (and derived classes). It will call the element's
         Calculate method to estimate the error based on the relative magnitude of the subscale.
         The result will be stored in the variable ERROR_RATIO for each element.

         This refinement criteria is based on G. Hauke, M. Doweidar, M. Miana, The multiscale approach to error
         estimation and adaptivity, Computer Methods in Applied Mechanics and Engineering,
         Volume 195, Issues 13-16, 2006
         @param rModelPart The model part containing the elements
         @see VMS for the calculation of the error ratio at the elemental level
         */
        void SubscaleErrorEstimate(ModelPart& rModelPart)
        {
            // Partitioning
            const int NumThreads = OpenMPUtils::GetNumThreads();
            OpenMPUtils::PartitionVector ElementPartition;
            OpenMPUtils::DivideInPartitions(rModelPart.Elements().size(),NumThreads,ElementPartition);

            // Initialize error values
            for( ModelPart::ElementIterator itElem = rModelPart.ElementsBegin(); itElem != rModelPart.ElementsEnd(); ++itElem)
                itElem->SetValue(ERROR_RATIO,0.0);

            //Compute average kinetic energy
            double Atot = 0.0;
            double avg_vel = 0.0;
            array_1d<double,3> vgauss;

            const double NodeFactor = 1.0 / static_cast<double>(rModelPart.ElementsBegin()->GetGeometry().size());
            const int ElementsEnd = static_cast<unsigned int>(rModelPart.Elements().size());
            
            #pragma omp parallel for reduction(+:Atot,avg_vel)
            for(int k=0; k< ElementsEnd; k++)
            {
                ModelPart::ElementsContainerType::iterator itElem = rModelPart.ElementsBegin()+k;

                Geometry<Node<3> >& geom = itElem->GetGeometry();
                double Area = geom.Area();

                noalias(vgauss) = geom[0].FastGetSolutionStepValue(VELOCITY);
                for(unsigned int i=1; i<geom.size(); i++)
                    noalias(vgauss) += geom[i].FastGetSolutionStepValue(VELOCITY);
                double norm_v = norm_2(vgauss) * NodeFactor;

                avg_vel += Area*norm_v;
                Atot += Area;
            }
            avg_vel/=Atot;

            if(Atot < 1e-10)
                KRATOS_ERROR(std::logic_error,"area can not be zero!!","")

            #pragma omp parallel
            {
                int k = OpenMPUtils::ThisThread();

                // Initialize the iterator boundaries for this thread
                ModelPart::ElementsContainerType::iterator ElemBegin = rModelPart.ElementsBegin() + ElementPartition[k];
                ModelPart::ElementsContainerType::iterator ElemEnd = rModelPart.ElementsBegin() + ElementPartition[k+1];

                ProcessInfo& rProcessInfo = rModelPart.GetProcessInfo();
                double Error;

                // Ask each element to calculate its Error ratio
                for( ModelPart::ElementsContainerType::iterator itElem = ElemBegin; itElem != ElemEnd; ++itElem)
                {
                    itElem->Calculate(ERROR_RATIO,Error,rProcessInfo);

                    Error/=avg_vel;
                    itElem->SetValue(ERROR_RATIO,Error);
                }
            }
        }

        ///@} //Operators

	private:

	};

        ///@} //Kratos Classes

        ///@}

}  // namespace Kratos.

#endif // KRATOS_MARK_FOR_REFINEMENT  defined



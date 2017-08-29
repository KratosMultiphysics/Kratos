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
//   Date:                August 6th 2017
//   Revision:            1.0
//
//

#if !defined(KRATOS_DRY_BED_UTILITY_H_INCLUDED)
#define  KRATOS_DRY_BED_UTILITY_H_INCLUDED

// System includes
#include <string>
#include <iostream> 
#include <algorithm>

// External includes 

// Project includes
#include "shallow_water_application.h"

namespace Kratos
{

    class DryBedUtility
    {
    public:
        
        KRATOS_CLASS_POINTER_DEFINITION(DryBedUtility);
        
        DryBedUtility(ModelPart& model_part)
            : mr_model_part(model_part)  
        {
            KRATOS_TRY
            std::cout << "Initializing dry/wet state utility" << std::endl; 
            KRATOS_CATCH("")
        }
        
        ~DryBedUtility()
        {}
        
        void CheckConservedVariables(ModelPart::NodesContainerType& rNodes)
        {
            KRATOS_TRY
            ModelPart::NodesContainerType::iterator inodebegin = rNodes.begin();
            vector<unsigned int> node_partition;
            #ifdef _OPENMP
                int number_of_threads = omp_get_max_threads();
            #else
                int number_of_threads = 1;
            #endif
            OpenMPUtils::CreatePartition(number_of_threads, rNodes.size(), node_partition);
            
            #pragma omp parallel for
            for(int kkk=0; kkk<number_of_threads; kkk++)
            {
                for(unsigned int ii=node_partition[kkk]; ii<node_partition[kkk+1]; ii++)
                {
                    ModelPart::NodesContainerType::iterator inode = inodebegin+ii;
                    if (inode->FastGetSolutionStepValue(HEIGHT) < 1e-3)
                    {
                        inode->GetSolutionStepValue(HEIGHT)     = 1e-3;
                        inode->GetSolutionStepValue(MOMENTUM_X) = 0;
                        inode->GetSolutionStepValue(MOMENTUM_Y) = 0;
                    }
                    
                }
            }
            
            KRATOS_CATCH("")
        }
        
        void CheckPrimitiveVariables(ModelPart::NodesContainerType& rNodes)
        {
            KRATOS_TRY
            ModelPart::NodesContainerType::iterator inodebegin = rNodes.begin();
            vector<unsigned int> node_partition;
            #ifdef _OPENMP
                int number_of_threads = omp_get_max_threads();
            #else
                int number_of_threads = 1;
            #endif
            OpenMPUtils::CreatePartition(number_of_threads, rNodes.size(), node_partition);
            
            #pragma omp parallel for
            for(int kkk=0; kkk<number_of_threads; kkk++)
            {
                for(unsigned int ii=node_partition[kkk]; ii<node_partition[kkk+1]; ii++)
                {
                    ModelPart::NodesContainerType::iterator inode = inodebegin+ii;
                    if (inode->FastGetSolutionStepValue(HEIGHT) < 1e-6)
                    {
                        inode->GetSolutionStepValue(HEIGHT)     = 1e-3;
                        inode->GetSolutionStepValue(VELOCITY_X) = 0;
                        inode->GetSolutionStepValue(VELOCITY_Y) = 0;
                    }
                    
                }
            }
            
            KRATOS_CATCH("")
        }
        
    protected:
    
    private:
    
        ModelPart& mr_model_part;
        
    }; // class DryBedUtility

}  // namespace Kratos

# endif // KRATOS_DRY_BED_UTILITY_H_INCLUDED

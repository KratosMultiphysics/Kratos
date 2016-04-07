#ifndef PRE_UTILITES_H
#define PRE_UTILITES_H

// Project includes
#include "utilities/timer.h"
#include "includes/variables.h"
#include "DEM_application.h"

/* System includes */
#include <limits>
#include <iostream>
#include <iomanip>

/* External includes */
#ifdef _OPENMP
#include <omp.h>
#endif

/* Project includes */
#include "includes/define.h"
#include "utilities/openmp_utils.h"

namespace Kratos
{

class PreUtilities
    {
     public:

     typedef ModelPart::ElementsContainerType                         ElementsArrayType;
     typedef ModelPart::NodesContainerType::ContainerType             NodesContainerType;
     typedef WeakPointerVector<Element>                               ParticleWeakVectorType;
     typedef WeakPointerVector<Element>::iterator                     ParticleWeakIteratorType;

     KRATOS_CLASS_POINTER_DEFINITION(PreUtilities);

      /// Default constructor

      PreUtilities(ModelPart& rModelPart)
      {
          //mInitialCenterOfMassAndMass = CalculateCenterOfMass(rModelPart);
          //mInitialMass                = CalculateTotalMass(rModelPart);
      }

      /// Destructor

      virtual ~PreUtilities() {}

      void MeasureTopHeight(ModelPart& rModelPart, double& subtotal, double& weight )
      {
        /*
        ElementsArrayType& pElements        = rModelPart.Elements();

        for (ElementsArrayType::iterator it= pElements.begin(); it!=pElements.end(); ++it)
        {
                 
            if( it->GetGeometry()[0].FastGetSolutionStepValue(GROUP_ID) == 1 )
            {
              
              ParticleWeakVectorType& mrNeighbours = it->GetValue(NEIGHBOUR_ELEMENTS);            
              
              for(ParticleWeakIteratorType ineighbour = mrNeighbours.begin();  
              ineighbour != mrNeighbours.end(); ineighbour++)
              {
                
                if( ineighbour->GetGeometry()[0].FastGetSolutionStepValue(GROUP_ID) != 1 )
                {
                
                    subtotal += it->GetGeometry()[0].Coordinates()[1]*it->GetGeometry()[0].FastGetSolutionStepValue(RADIUS);
                    weight += it->GetGeometry()[0].FastGetSolutionStepValue(RADIUS);
                    
                    break;
                }
                  
              }
              
            }
            
         }
        */
         
      }
          
     void MeasureBotHeight(ModelPart& rModelPart, double& subtotal, double& weight )
      {
          /*
            ElementsArrayType& pElements        = rModelPart.Elements();

            for (ElementsArrayType::iterator it= pElements.begin(); it!=pElements.end(); ++it)
            {

                if( it->GetGeometry()[0].FastGetSolutionStepValue(GROUP_ID) == 2 )
                {

                  ParticleWeakVectorType& mrNeighbours = it->GetValue(NEIGHBOUR_ELEMENTS);

                  for(ParticleWeakIteratorType ineighbour = mrNeighbours.begin();
                  ineighbour != mrNeighbours.end(); ineighbour++)
                  {

                    if( ineighbour->GetGeometry()[0].FastGetSolutionStepValue(GROUP_ID) != 2 )
                    {

                        subtotal += it->GetGeometry()[0].Coordinates()[1]*it->GetGeometry()[0].FastGetSolutionStepValue(RADIUS);
                        weight += it->GetGeometry()[0].FastGetSolutionStepValue(RADIUS);

                        break;
                    }

                  }

                }

             }
         */
      }

        array_1d<double, 3> GetInitialCenterOfMass()
        {
            return mInitialCenterOfMassAndMass;
        }

        /// Turn back information as a stemplate<class T, std::size_t dim> tring.

        virtual std::string Info() const
        {
            return "";
        }

        /// Print information about this object.

        virtual void PrintInfo(std::ostream& rOStream) const
        {
        }

        /// Print object's data.

        virtual void PrintData(std::ostream& rOStream) const
        {
        }

        vector<unsigned int>&    GetElementPartition(){return (mElementPartition);};

    protected:

        vector<unsigned int>                mElementPartition;

    private:

        array_1d<double, 3> mInitialCenterOfMassAndMass;
        double mInitialMass;

        /// Assignment operator
        PreUtilities & operator=(PreUtilities const& rOther);

    }; // Class PreUtilities

/// output stream function
// 	template<std::size_t TDim>
// 	inline std::ostream& operator << (std::ostream& rOStream)
// 	{
// 		rThis.PrintInfo(rOStream);
// 		rOStream << std::endl;
// 		rThis.PrintData(rOStream);
//
// 		return rOStream;
// 	}

} // namespace Kratos.

#endif // PRE_UTILITES_H

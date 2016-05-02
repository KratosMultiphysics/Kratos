#ifndef PRE_UTILITES_H
#define PRE_UTILITES_H

// Project includes
#include "utilities/timer.h"
#include "includes/variables.h"
//#include "DEM_application.h"

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
#include "cluster_information.h"

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
      
      void SetClusterInformationInProperties(std::string const& name,
                                            boost::python::list& list_of_coordinates, 
                                            boost::python::list& list_of_radii, 
                                            double size, 
                                            double volume, 
                                            boost::python::list& inertias, 
                                            Properties::Pointer& p_properties) {
          
          ClusterInformation cl_info;
          
          cl_info.mName = name;
          
          array_1d<double,3> coords(3,0.0);
          for(int i=0; i<boost::python::len(list_of_coordinates); i++){
            boost::python::list list(list_of_coordinates[i]);
            coords[0] =  boost::python::extract<double>(list[0]);
            coords[1] =  boost::python::extract<double>(list[1]);
            coords[2] =  boost::python::extract<double>(list[2]);
            cl_info.mListOfCoordinates.push_back(coords);
          }
          for(int i=0; i<boost::python::len(list_of_radii); i++){
            cl_info.mListOfRadii.push_back(boost::python::extract<double>(list_of_radii[i]));
          }
          //TODO: check the sizes (should be the same)
          cl_info.mSize = size;
          cl_info.mVolume = volume;
          cl_info.mInertias[0] = boost::python::extract<double>(inertias[0]);
          cl_info.mInertias[1] = boost::python::extract<double>(inertias[1]);
          cl_info.mInertias[2] = boost::python::extract<double>(inertias[2]);
          
          p_properties->SetValue(CLUSTER_INFORMATION, cl_info);
          
      }

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

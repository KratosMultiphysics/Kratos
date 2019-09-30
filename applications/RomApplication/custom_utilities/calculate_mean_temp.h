//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//  Kratos default license: kratos/license.txt
//
//  Main authors:    YOUR_NAME_HERE
//

#if !defined(KRATOS_CALCUALTE_MEAN_TEMP_UTILITY_INCLUDED )
#define  KRATOS_CALCUALTE_MEAN_TEMP_UTILITY_INCLUDED

// System includes
#include <string>
#include <iostream> 
#include <algorithm>

// Project includes 
#include "includes/define.h"
#include "includes/variables.h" 
#include "includes/model_part.h"
#include "input_output/logger.h"

namespace Kratos
{
    // This class is to be modified by the user to customize the interpolation process
    class CalculateMeanTemperature 
    {
    public:

        KRATOS_CLASS_POINTER_DEFINITION(CalculateMeanTemperature);

        CalculateMeanTemperature(ModelPart& rModelPart)
            : mrModelPart(rModelPart) //mrModelPart is saved as private variable (declared at the end of the file)  
        {
            KRATOS_TRY
            KRATOS_INFO("CalculateMeanTemperature") << "Hello, I am the constructor of the Utility" << std::endl; 
            KRATOS_CATCH("")
        }

        ~CalculateMeanTemperature()
        {}

        void Calculate()
        {
            KRATOS_TRY

            double area;              // We create the needed variables
            double sum_areas=0.0;
            double sum_temperatures=0.0;
            double one_third=1.0/3.0;

            // Getting data for the given geometry
            for(auto it_elem = mrModelPart.ElementsBegin(); it_elem!=mrModelPart.ElementsEnd(); ++it_elem) // Loop the elements
            {
                Geometry<Node<3> >& geom = it_elem->GetGeometry(); 
                area = CalculateArea(geom);               // We call CalculateArea (private function)  
                sum_areas += area;
                for (unsigned int k = 0; k < 3; ++k) {
                    sum_temperatures += geom[k].FastGetSolutionStepValue(TEMPERATURE) * one_third * area;
                }
            }
            const double mean_temperature = sum_temperatures / sum_areas;
            KRATOS_INFO("CalculateMeanTemperature") << "Finished, the mean temperature is " << mean_temperature << std::endl;   //we print the result  

            KRATOS_CATCH("")
        } 
    
    protected:


    private:
        double CalculateArea(Element::GeometryType& geom)
        {
            return 0.5 * ((geom[1].X() - geom[0].X())*(geom[2].Y() - geom[0].Y())- (geom[1].Y() - geom[0].Y())*(geom[2].X() - geom[0].X()));
        }

        ModelPart& mrModelPart;

    };

}  // namespace Kratos.

#endif // KRATOS_CALCUALTE_MEAN_TEMP_UTILITY_INCLUDED  defined
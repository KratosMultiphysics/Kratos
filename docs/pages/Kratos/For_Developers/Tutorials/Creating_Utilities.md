---
title: Creating an Utility
keywords: 
tags: [Tutorial-Creating-an-Utility.md]
sidebar: kratos_for_developers
summary: 
---

## Overview

Utilities are useful tools to perform arbitrary actions. Furthermore, it is possible to create the solver to your specific problems and avoid using the builder and solver, for example to create a new edge based builder. For educational purposes, we will create here a simple utility that calculates the mean temperature inside the domain once the solution has been obtained. 

## Modifying add_custom_utilites_to_python.cpp

To begin with we have to tell Kratos that we have an Utility. Doing so is straightforward since we only need to edit the file **custom_python/add_custom_utilities_to_python.cpp** , including the source of our utility and some lines that serve as an interfase between the **C++** functions and python. Once you've done this Kratos knows you have a new utility and the functions you've declared in this file can be called directly from python. Have in mind that in the constructor you must write the input arguments, while in all the other functions you just type the name of them. 

```cpp
//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics 
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    YOUR_NAME_HERE
//

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/define_python.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_utilities/calculate_mean_temp.h"

namespace Kratos
{

namespace Python
{
using namespace pybind11;

  void  AddCustomUtilitiesToPython(pybind11::module& m)
  {
    class_<CalculateMeanTemperature, typename CalculateMeanTemperature::Pointer>(m, "CalculateMeanTemperature")
    .def(init<ModelPart&>()) // The input parameters is a model part 
    .def("Execute",&CalculateMeanTemperature::Calculate) // When we call "Execute" in python, Calculate is called in C++. Notice we don't write the input parameters here 
    ;

  }

}  // namespace Python.

} // Namespace Kratos

```

## Creating custom_utilities/calculate_mean_temp.h

```cpp
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
                Geometry<Node >& geom = it_elem->GetGeometry(); 
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
```

And that's it, this is all you need, to use it you must first call the constructor and then execute it. 
